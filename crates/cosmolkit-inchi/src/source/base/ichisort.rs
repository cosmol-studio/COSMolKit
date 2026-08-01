use crate::source::base::ichimap2::DifferentiateRanks2;
use crate::source::base::util::{inchi_calloc, inchi_free};
use crate::source_types::{
    AT_NUMB, AT_RANK, BOND_DOUBLE, CANON_GLOBALS, CANON_STAT, CT_OUT_OF_RAM, NEIGH_LIST, S_CHAR,
    SourceConstPointer, SourceHeap, SourceHeapError, SourceMutPointer, T_GROUP_INFO, sp_ATOM,
    tagAtInvariantIndexes_AT_INV_BREAK1, tagAtInvariantIndexes_AT_INV_LENGTH,
};

#[allow(non_snake_case)]
pub(crate) enum CompNeighborsATNumberContext<'a> {
    Globals {
        heap: &'a SourceHeap,
        globals: &'a CANON_GLOBALS,
    },
    Slices {
        neighbors: &'a [AT_NUMB],
        ranks: &'a [AT_RANK],
    },
}

#[allow(non_snake_case)]
pub(crate) fn CompNeighborsAT_NUMBER(
    a1: AT_NUMB,
    a2: AT_NUMB,
    context: CompNeighborsATNumberContext<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:453 CompNeighborsAT_NUMBER
    // INCHI✔️✔️: complete active CT_NEIGH_INCREASE source frame follows verbatim.
    /*
    int CompNeighborsAT_NUMBER( const void* a1, const void* a2, void *p )
    {
        CANON_GLOBALS *pCG = (CANON_GLOBALS *) p;
    #ifdef CT_NEIGH_INCREASE
        return (int) pCG->m_pn_RankForSort[pCG->m_pNeighborsForSort[( int )*(const AT_NUMB*) a1]] -
               (int) pCG->m_pn_RankForSort[pCG->m_pNeighborsForSort[( int )*(const AT_NUMB*) a2]];
    #else
        return (int) ( (CANON_GLOBALS *) pCG )->m_pn_RankForSort[pNeighborsForSort[( int )*(const AT_NUMB*) a2]] -
               (int) ( (CANON_GLOBALS *) pCG )->m_pn_RankForSort[pNeighborsForSort[( int )*(const AT_NUMB*) a1]];
    #endif
    }
    */
    // END INCHI C FUNCTION: CompNeighborsAT_NUMBER

    let (neighbors, ranks) = match context {
        CompNeighborsATNumberContext::Globals { heap, globals } => (
            heap.slice(globals.m_pNeighborsForSort)?,
            heap.slice(globals.m_pn_RankForSort)?,
        ),
        CompNeighborsATNumberContext::Slices { neighbors, ranks } => (neighbors, ranks),
    };
    let first_neighbor = *neighbors
        .get(usize::from(a1))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second_neighbor = *neighbors
        .get(usize::from(a2))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let first_rank = *ranks
        .get(usize::from(first_neighbor))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second_rank = *ranks
        .get(usize::from(second_neighbor))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    Ok(i32::from(first_rank) - i32::from(second_rank))
}

#[allow(non_snake_case)]
pub(crate) fn comp_AT_RANK(a1: AT_RANK, a2: AT_RANK) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:467 comp_AT_RANK
    // INCHI✔️✔️: int comp_AT_RANK( const void* a1, const void* a2, void *p )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     return ( int )*(const AT_RANK*) a1 - ( int )*(const AT_RANK*) a2;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: comp_AT_RANK

    i32::from(a1) - i32::from(a2)
}

#[allow(non_snake_case)]
pub(crate) fn CompAtomInvariants2Only(
    heap: &SourceHeap,
    a1: AT_RANK,
    a2: AT_RANK,
    pCG: &CANON_GLOBALS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:502 CompAtomInvariants2Only
    // INCHI✔️✔️: int CompAtomInvariants2Only( const void* a1, const void* a2, void *p )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     CANON_GLOBALS *pCG = (CANON_GLOBALS *) p;
    // INCHI✔️✔️:     const ATOM_INVARIANT2 *pAI1 = pCG->m_pAtomInvariant2ForSort + ( int )*(const AT_RANK*) a1;
    // INCHI✔️✔️:     const ATOM_INVARIANT2 *pAI2 = pCG->m_pAtomInvariant2ForSort + ( int )*(const AT_RANK*) a2;
    // INCHI✔️✔️:     int i;
    // INCHI✔️✔️:     for (i = 0; i < AT_INV_BREAK1; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (pAI1->val[i] == pAI2->val[i])
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         return  (int) pAI1->val[i] - (int) pAI2->val[i];
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (pAI1->iso_sort_key != pAI2->iso_sort_key)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return ( pAI1->iso_sort_key > pAI2->iso_sort_key ) ? 1 : -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     for (; i < AT_INV_LENGTH; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (pAI1->val[i] != pAI2->val[i])
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         return  (int) pAI1->val[i] - (int) pAI2->val[i];
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (pAI1->iso_aux_key != pAI2->iso_aux_key)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return ( pAI1->iso_aux_key > pAI2->iso_aux_key ) ? 1 : -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: CompAtomInvariants2Only

    let invariants = heap.slice(pCG.m_pAtomInvariant2ForSort)?;
    let first = invariants
        .get(usize::from(a1))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second = invariants
        .get(usize::from(a2))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut i = 0_usize;
    while i < tagAtInvariantIndexes_AT_INV_BREAK1 as usize {
        if first.val[i] != second.val[i] {
            return Ok(i32::from(first.val[i]) - i32::from(second.val[i]));
        }
        i += 1;
    }
    if first.iso_sort_key != second.iso_sort_key {
        return Ok(if first.iso_sort_key > second.iso_sort_key {
            1
        } else {
            -1
        });
    }
    while i < tagAtInvariantIndexes_AT_INV_LENGTH as usize {
        if first.val[i] != second.val[i] {
            i += 1;
            continue;
        }
        return Ok(i32::from(first.val[i]) - i32::from(second.val[i]));
    }
    if first.iso_aux_key != second.iso_aux_key {
        return Ok(if first.iso_aux_key > second.iso_aux_key {
            1
        } else {
            -1
        });
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CompAtomInvariants2(
    heap: &SourceHeap,
    a1: AT_RANK,
    a2: AT_RANK,
    pCG: &CANON_GLOBALS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:536 CompAtomInvariants2
    // INCHI✔️✔️: int CompAtomInvariants2( const void* a1, const void* a2, void *p )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     /*  Warning: the following line may be compiler implementation dependent */
    // INCHI✔️✔️:     int ret = CompAtomInvariants2Only( a1, a2, p );
    // INCHI✔️✔️:     if (!ret)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ret = ( int )*(const AT_RANK*) a1 - ( int )*(const AT_RANK*) a2;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return ret;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: CompAtomInvariants2

    let ret = CompAtomInvariants2Only(heap, a1, a2, pCG)?;
    Ok(if ret == 0 {
        i32::from(a1) - i32::from(a2)
    } else {
        ret
    })
}

#[allow(non_snake_case)]
pub(crate) fn CompRank(
    heap: &SourceHeap,
    a1: AT_RANK,
    a2: AT_RANK,
    pCG: &CANON_GLOBALS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:475 CompRank
    // INCHI✔️✔️: int CompRank( const void* a1, const void* a2, void *p )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     CANON_GLOBALS *pCG = (CANON_GLOBALS *) p;
    // INCHI✔️✔️:     int ret = (int) pCG->m_pn_RankForSort[( int )*(const AT_RANK*) a1] -
    // INCHI✔️✔️:               (int) pCG->m_pn_RankForSort[( int )*(const AT_RANK*) a2];
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return ret;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: CompRank

    let ranks = heap.slice(pCG.m_pn_RankForSort)?;
    let first = ranks
        .get(usize::from(a1))
        .copied()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second = ranks
        .get(usize::from(a2))
        .copied()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    Ok(i32::from(first) - i32::from(second))
}

#[allow(non_snake_case)]
pub(crate) fn CompRanksOrd(
    heap: &SourceHeap,
    a1: AT_RANK,
    a2: AT_RANK,
    pCG: &CANON_GLOBALS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:486 CompRanksOrd
    // INCHI✔️✔️: int CompRanksOrd( const void* a1, const void* a2, void *p )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int ret;
    // INCHI✔️✔️:     CANON_GLOBALS *pCG = (CANON_GLOBALS *) p;
    // INCHI✔️✔️:     ret = (int) pCG->m_pn_RankForSort[( int )*(const AT_RANK*) a1] -
    // INCHI✔️✔️:           (int) pCG->m_pn_RankForSort[( int )*(const AT_RANK*) a2];
    // INCHI✔️✔️:     if (!ret)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ret = ( int )*(const AT_RANK*) a1 - ( int )*(const AT_RANK*) a2;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return ret;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: CompRanksOrd

    let ranks = heap.slice(pCG.m_pn_RankForSort)?;
    let first = ranks
        .get(usize::from(a1))
        .copied()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second = ranks
        .get(usize::from(a2))
        .copied()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let ret = i32::from(first) - i32::from(second);
    Ok(if ret == 0 {
        i32::from(a1) - i32::from(a2)
    } else {
        ret
    })
}

#[allow(non_snake_case)]
pub(crate) fn insertions_sort_AT_NUMBERS(
    heap: &mut SourceHeap,
    base: SourceMutPointer<AT_NUMB>,
    num: i32,
    compare: &mut dyn FnMut(&SourceHeap, AT_NUMB, AT_NUMB) -> Result<i32, SourceHeapError>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:331 insertions_sort_AT_NUMBERS
    // INCHI✔️❌: int insertions_sort_AT_NUMBERS( void *pCG,
    // INCHI✔️❌:                                 AT_NUMB *base,
    // INCHI✔️❌:                                 int num,
    // INCHI✔️❌:                                 int( *compare )( const void *e1, const void *e2, void * ) )
    // INCHI✔️❌: {
    // INCHI✔️❌:     AT_NUMB *i, *j, *pk, tmp;
    // INCHI✔️❌:     int  k, num_trans = 0;
    // INCHI✔️❌:     for (k = 1, pk = base; k < num; k++, pk++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (j = ( i = pk ) + 1, tmp = *j; j > base && ( *compare )( i, &tmp, pCG ) > 0; j = i, i--)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *j = *i;
    // INCHI✔️❌:             num_trans++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         *j = tmp;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return num_trans;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: insertions_sort_AT_NUMBERS

    let mut num_trans = 0_i32;
    let mut k = 1_i32;
    while k < num {
        let tmp = source_at_number_get(heap, base, k)?;
        let mut j = k;
        let mut i = k.wrapping_sub(1);
        while j > 0 {
            let left = source_at_number_get(heap, base, i)?;
            if compare(heap, left, tmp)? <= 0 {
                break;
            }
            source_at_number_set(heap, base, j, left)?;
            num_trans = num_trans.wrapping_add(1);
            j = i;
            i = i.wrapping_sub(1);
        }
        source_at_number_set(heap, base, j, tmp)?;
        k = k.wrapping_add(1);
    }
    Ok(num_trans)
}

#[allow(non_snake_case)]
pub(crate) fn CompChemElemLex(a1: &[u8], a2: &[u8]) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:551 CompChemElemLex
    // INCHI✔️✔️: int CompChemElemLex( const void *a1, const void *a2 )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     return memcmp( a1, a2, 2 );
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: CompChemElemLex

    let first = a1.get(..2).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second = a2.get(..2).ok_or(SourceHeapError::PointerOutOfBounds)?;
    for index in 0..2 {
        if first[index] != second[index] {
            return Ok(i32::from(first[index]) - i32::from(second[index]));
        }
    }
    Ok(0)
}

fn source_at_number_get(
    heap: &SourceHeap,
    pointer: SourceMutPointer<AT_NUMB>,
    index: i32,
) -> Result<AT_NUMB, SourceHeapError> {
    heap.slice(pointer.as_const())?
        .get(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .copied()
        .ok_or(SourceHeapError::PointerOutOfBounds)
}

fn source_at_number_set(
    heap: &mut SourceHeap,
    pointer: SourceMutPointer<AT_NUMB>,
    index: i32,
    value: AT_NUMB,
) -> Result<(), SourceHeapError> {
    heap.slice_mut(pointer)?
        .get_mut(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)
        .map(|slot| *slot = value)
}

#[allow(non_snake_case)]
pub(crate) fn CompareNeighListLex(
    heap: &SourceHeap,
    pp1: NEIGH_LIST,
    pp2: NEIGH_LIST,
    nRank: SourceMutPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:560 CompareNeighListLex
    // INCHI✔️❌: int CompareNeighListLex( NEIGH_LIST pp1, NEIGH_LIST pp2, const AT_RANK *nRank )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int len1 = (int) *pp1++;
    // INCHI✔️❌:     int len2 = (int) *pp2++;
    // INCHI✔️❌:     int len = inchi_min( len1, len2 );
    // INCHI✔️❌:     int diff = 0;
    // INCHI✔️❌:     int ret;
    // INCHI✔️❌:     /* djb-rwth: fixing oss-fuzz issue #25642 */
    // INCHI✔️❌:     while ((len > 0) && !diff)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         len--;
    // INCHI✔️❌:         diff = (int)nRank[*pp1++] - (int)nRank[*pp2++];
    // INCHI✔️❌:     };
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = diff ? diff : (len1 - len2);
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompareNeighListLex

    let len1 = i32::from(source_at_number_get(heap, pp1, 0)?);
    let len2 = i32::from(source_at_number_get(heap, pp2, 0)?);
    let mut len = len1.min(len2);
    let mut diff = 0_i32;
    let mut index = 1_i32;
    while len > 0 && diff == 0 {
        len = len.wrapping_sub(1);
        let atom1 = source_at_number_get(heap, pp1, index)?;
        let atom2 = source_at_number_get(heap, pp2, index)?;
        diff = i32::from(source_at_number_get(heap, nRank, i32::from(atom1))?)
            - i32::from(source_at_number_get(heap, nRank, i32::from(atom2))?);
        index = index.wrapping_add(1);
    }

    Ok(if diff != 0 { diff } else { len1 - len2 })
}

#[allow(non_snake_case)]
pub(crate) fn CompareNeighListLexUpToMaxRank(
    heap: &SourceHeap,
    pp1: NEIGH_LIST,
    pp2: NEIGH_LIST,
    nRank: SourceMutPointer<AT_RANK>,
    nMaxAtNeighRank: AT_RANK,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:582 CompareNeighListLexUpToMaxRank
    // INCHI✔️❌: int CompareNeighListLexUpToMaxRank( NEIGH_LIST pp1, NEIGH_LIST pp2, const AT_RANK *nRank, AT_RANK nMaxAtNeighRank )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int len1 = (int) *pp1++;
    // INCHI✔️❌:     int len2 = (int) *pp2++;
    // INCHI✔️❌:     int diff = 0;
    // INCHI✔️❌:     int len;
    // INCHI✔️❌:     while (0 < len1 && nRank[pp1[len1 - 1]] > nMaxAtNeighRank)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         len1--;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     while (0 < len2 && nRank[pp2[len2 - 1]] > nMaxAtNeighRank)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         len2--;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     len = inchi_min( len1, len2 );
    // INCHI✔️❌:     while (len-- > 0 && !( diff = (int) nRank[*pp1++] - (int) nRank[*pp2++] ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return diff ? diff : ( len1 - len2 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompareNeighListLexUpToMaxRank

    let mut len1 = i32::from(source_at_number_get(heap, pp1, 0)?);
    let mut len2 = i32::from(source_at_number_get(heap, pp2, 0)?);
    while len1 > 0 {
        let atom = source_at_number_get(heap, pp1, len1)?;
        if source_at_number_get(heap, nRank, i32::from(atom))? <= nMaxAtNeighRank {
            break;
        }
        len1 = len1.wrapping_sub(1);
    }
    while len2 > 0 {
        let atom = source_at_number_get(heap, pp2, len2)?;
        if source_at_number_get(heap, nRank, i32::from(atom))? <= nMaxAtNeighRank {
            break;
        }
        len2 = len2.wrapping_sub(1);
    }

    let mut len = len1.min(len2);
    let mut diff = 0_i32;
    let mut index = 1_i32;
    while len > 0 {
        len = len.wrapping_sub(1);
        let atom1 = source_at_number_get(heap, pp1, index)?;
        let atom2 = source_at_number_get(heap, pp2, index)?;
        diff = i32::from(source_at_number_get(heap, nRank, i32::from(atom1))?)
            - i32::from(source_at_number_get(heap, nRank, i32::from(atom2))?);
        index = index.wrapping_add(1);
        if diff != 0 {
            break;
        }
    }
    Ok(if diff != 0 { diff } else { len1 - len2 })
}

#[allow(non_snake_case)]
pub(crate) fn compare_NeighLists(
    heap: &SourceHeap,
    op1: NEIGH_LIST,
    op2: NEIGH_LIST,
    pCG: &CANON_GLOBALS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:607 compare_NeighLists
    // INCHI✔️❌: int compare_NeighLists( const NEIGH_LIST *op1, const NEIGH_LIST *op2, void *p )
    // INCHI✔️❌: {
    // INCHI✔️❌:     CANON_GLOBALS *pCG = (CANON_GLOBALS *) p;
    // INCHI✔️❌:     return CompareNeighListLex( *op1, *op2, pCG->m_pn_RankForSort );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: compare_NeighLists

    CompareNeighListLex(heap, op1, op2, pCG.m_pn_RankForSort.as_mut())
}

#[allow(non_snake_case)]
pub(crate) fn CompNeighListRanks(
    heap: &SourceHeap,
    a1: AT_RANK,
    a2: AT_RANK,
    pCG: &CANON_GLOBALS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:615 CompNeighListRanks
    // INCHI✔️❌: int CompNeighListRanks( const void* a1, const void* a2, void *p )
    // INCHI✔️❌: {
    // INCHI✔️❌:     CANON_GLOBALS *pCG = (CANON_GLOBALS *) p;
    // INCHI✔️❌:     int ret;
    // INCHI✔️❌:     ret = (int) pCG->m_pn_RankForSort[*( (const AT_RANK*) a1 )] -
    // INCHI✔️❌:           (int) pCG->m_pn_RankForSort[*( (const AT_RANK*) a2 )];
    // INCHI✔️❌:     if (!ret)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = compare_NeighLists( pCG->m_pNeighList_RankForSort + *( (const AT_RANK*) a1 ),
    // INCHI✔️❌:                                      pCG->m_pNeighList_RankForSort + *( (const AT_RANK*) a2 ), p );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompNeighListRanks

    let ranks = heap.slice(pCG.m_pn_RankForSort)?;
    let first_index = usize::from(a1);
    let second_index = usize::from(a2);
    let first_rank = *ranks
        .get(first_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second_rank = *ranks
        .get(second_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let ret = i32::from(first_rank) - i32::from(second_rank);
    if ret != 0 {
        return Ok(ret);
    }

    let neighbor_lists = heap.slice(pCG.m_pNeighList_RankForSort)?;
    let first = *neighbor_lists
        .get(first_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second = *neighbor_lists
        .get(second_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    compare_NeighLists(heap, first, second, pCG)
}

#[allow(non_snake_case)]
pub(crate) fn CompNeighLists(
    heap: &SourceHeap,
    a1: AT_RANK,
    a2: AT_RANK,
    pCG: &CANON_GLOBALS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:632 CompNeighLists
    // INCHI✔️❌: int CompNeighLists( const void* a1, const void* a2, void *p )
    // INCHI✔️❌: {
    // INCHI✔️❌:     CANON_GLOBALS *pCG = (CANON_GLOBALS *) p;
    // INCHI✔️❌:     int ret;
    // INCHI✔️❌:     ret = compare_NeighLists( pCG->m_pNeighList_RankForSort + *( (const AT_RANK*) a1 ),
    // INCHI✔️❌:                               pCG->m_pNeighList_RankForSort + *( (const AT_RANK*) a2 ), p );
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompNeighLists
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: CompNeighLists
    // INCHI✔️❌: typedef unsigned short AT_RANK;
    // INCHI✔️❌: typedef AT_RANK  *NEIGH_LIST;
    // INCHI✔️❌:     const NEIGH_LIST      *m_pNeighList_RankForSort;
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: CompNeighLists

    let neighbor_lists = heap.slice(pCG.m_pNeighList_RankForSort)?;
    let first = *neighbor_lists
        .get(usize::from(a1))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second = *neighbor_lists
        .get(usize::from(a2))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    compare_NeighLists(heap, first, second, pCG)
}

#[allow(non_snake_case)]
pub(crate) fn CompNeighListsUpToMaxRank(
    heap: &SourceHeap,
    a1: AT_RANK,
    a2: AT_RANK,
    pCG: &CANON_GLOBALS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:644 CompNeighListsUpToMaxRank
    // INCHI✔️❌: int CompNeighListsUpToMaxRank( const void* a1, const void* a2, void *p )
    // INCHI✔️❌: {
    // INCHI✔️❌:     CANON_GLOBALS *pCG = (CANON_GLOBALS *) p;
    // INCHI✔️❌:     int ret;
    // INCHI✔️❌:     ret = CompareNeighListLexUpToMaxRank( pCG->m_pNeighList_RankForSort[*( (const AT_RANK*) a1 )],
    // INCHI✔️❌:                                           pCG->m_pNeighList_RankForSort[*( (const AT_RANK*) a2 )],
    // INCHI✔️❌:                                           pCG->m_pn_RankForSort, pCG->m_nMaxAtNeighRankForSort );
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompNeighListsUpToMaxRank
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: CompNeighListsUpToMaxRank
    // INCHI✔️❌: typedef unsigned short AT_RANK;
    // INCHI✔️❌: typedef AT_RANK  *NEIGH_LIST;
    // INCHI✔️❌:     const NEIGH_LIST      *m_pNeighList_RankForSort;
    // INCHI✔️❌:     const AT_RANK         *m_pn_RankForSort;
    // INCHI✔️❌:     AT_RANK m_nMaxAtNeighRankForSort;
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: CompNeighListsUpToMaxRank

    let neighbor_lists = heap.slice(pCG.m_pNeighList_RankForSort)?;
    let first = *neighbor_lists
        .get(usize::from(a1))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second = *neighbor_lists
        .get(usize::from(a2))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    CompareNeighListLexUpToMaxRank(
        heap,
        first,
        second,
        pCG.m_pn_RankForSort.as_mut(),
        pCG.m_nMaxAtNeighRankForSort,
    )
}

#[allow(non_snake_case)]
pub(crate) fn CompNeighListRanksOrd(
    heap: &SourceHeap,
    a1: AT_RANK,
    a2: AT_RANK,
    pCG: &CANON_GLOBALS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:657 CompNeighListRanksOrd
    // INCHI✔️❌: int CompNeighListRanksOrd( const void* a1, const void* a2, void *p )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = CompNeighListRanks( a1, a2, p );
    // INCHI✔️❌:     if (!ret)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = ( int )*( (const AT_RANK*) a1 ) - ( int )*( (const AT_RANK*) a2 ); /*  keep original order if identical */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompNeighListRanksOrd

    let ret = CompNeighListRanks(heap, a1, a2, pCG)?;
    Ok(if ret == 0 {
        i32::from(a1) - i32::from(a2)
    } else {
        ret
    })
}

#[allow(non_snake_case)]
pub(crate) fn CompRanksInvOrd(a1: AT_RANK, a2: AT_RANK) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:670 CompRanksInvOrd
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    int CompRanksInvOrd( const void* a1, const void* a2, void *p )
    {
        return ( int )*(const AT_RANK*) a2 - ( int )*(const AT_RANK*) a1;
    }
    */
    // END INCHI C FUNCTION: CompRanksInvOrd

    i32::from(a2) - i32::from(a1)
}

#[allow(non_snake_case)]
pub(crate) fn CompNeighborsRanksCountEql(
    a1: AT_RANK,
    a2: AT_RANK,
    ranks: &[AT_RANK],
    pCG: &mut CANON_GLOBALS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:677 CompNeighborsRanksCountEql
    // INCHI✔️✔️: complete source frame follows verbatim; CT_NEIGH_INCREASE is active.
    /*
    int CompNeighborsRanksCountEql( const void* a1, const void* a2, void *p )
    {
        CANON_GLOBALS *pCG = (CANON_GLOBALS *) p;
    #ifdef CT_NEIGH_INCREASE
        int ret = (int) pCG->m_pn_RankForSort[( int )*(const AT_RANK*) a1] -
                  (int) pCG->m_pn_RankForSort[( int )*(const AT_RANK*) a2];
    #else
        int ret = (int) pCG->m_pn_RankForSort[( int )*(const AT_RANK*) a2] -
                  (int) pCG->m_pn_RankForSort[( int )*(const AT_RANK*) a1];
    #endif
        pCG->m_nNumCompNeighborsRanksCountEql += !ret;

        return ret;
    }
    */
    // END INCHI C FUNCTION: CompNeighborsRanksCountEql

    let first = *ranks
        .get(usize::from(a1))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second = *ranks
        .get(usize::from(a2))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let ret = i32::from(first) - i32::from(second);
    pCG.m_nNumCompNeighborsRanksCountEql = pCG
        .m_nNumCompNeighborsRanksCountEql
        .wrapping_add(i32::from(ret == 0));
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn insertions_sort_NeighList_AT_NUMBERS(
    heap: &mut SourceHeap,
    base: NEIGH_LIST,
    nRank: SourceMutPointer<AT_RANK>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:355 insertions_sort_NeighList_AT_NUMBERS
    // INCHI✔️❌: void insertions_sort_NeighList_AT_NUMBERS( NEIGH_LIST base, AT_RANK *nRank )
    // INCHI✔️❌: {
    // INCHI✔️❌:     AT_NUMB *i, *j, *pk, tmp;
    // INCHI✔️❌:     AT_RANK rj; /* optimization */
    // INCHI✔️❌:     int k, num = (int) *base++;
    // INCHI✔️❌:     for (k = 1, pk = base; k < num; k++, pk++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (j = ( i = pk ) + 1, rj = nRank[(int) *j]; j > base && nRank[(int) *i] > rj; j = i, i--)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             tmp = *i;
    // INCHI✔️❌:             *i = *j;
    // INCHI✔️❌:             *j = tmp;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: insertions_sort_NeighList_AT_NUMBERS

    let num = i32::from(source_at_number_get(heap, base, 0)?);
    let mut k = 1_i32;
    while k < num {
        let mut j = k.wrapping_add(1);
        let mut i = k;
        let current = source_at_number_get(heap, base, j)?;
        let rj = source_at_number_get(heap, nRank, i32::from(current))?;
        while j > 1 {
            let previous = source_at_number_get(heap, base, i)?;
            if source_at_number_get(heap, nRank, i32::from(previous))? <= rj {
                break;
            }
            let tmp = previous;
            let current = source_at_number_get(heap, base, j)?;
            source_at_number_set(heap, base, i, current)?;
            source_at_number_set(heap, base, j, tmp)?;
            j = i;
            i = i.wrapping_sub(1);
        }
        k = k.wrapping_add(1);
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn insertions_sort_NeighList_AT_NUMBERS3(
    heap: &mut SourceHeap,
    base: NEIGH_LIST,
    nRank: SourceMutPointer<AT_RANK>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:396 insertions_sort_NeighList_AT_NUMBERS3
    // INCHI✔️❌: int insertions_sort_NeighList_AT_NUMBERS3( NEIGH_LIST base, AT_RANK *nRank )
    // INCHI✔️❌: {
    // INCHI✔️❌:     AT_NUMB *i, *j, *pk, tmp;
    // INCHI✔️❌:     AT_RANK rj;
    // INCHI✔️❌:     int k, n, num = (int) *base++;
    // INCHI✔️❌:     for (k = 1, pk = base, n = 0; k < num; k++, pk++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (j = ( i = pk ) + 1, rj = nRank[(int) ( tmp = *j )];
    // INCHI✔️❌:              j > base && nRank[(int) *i] > rj;
    // INCHI✔️❌:              j = i, i--)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *j = *i;
    // INCHI✔️❌:             n++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         *j = tmp;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return n;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: insertions_sort_NeighList_AT_NUMBERS3

    let num = i32::from(source_at_number_get(heap, base, 0)?);
    let mut n = 0_i32;
    let mut k = 1_i32;
    while k < num {
        let tmp = source_at_number_get(heap, base, k.wrapping_add(1))?;
        let rj = source_at_number_get(heap, nRank, i32::from(tmp))?;
        let mut j = k.wrapping_add(1);
        let mut i = k;
        while j > 1 {
            let previous = source_at_number_get(heap, base, i)?;
            if source_at_number_get(heap, nRank, i32::from(previous))? <= rj {
                break;
            }
            source_at_number_set(heap, base, j, previous)?;
            n = n.wrapping_add(1);
            j = i;
            i = i.wrapping_sub(1);
        }
        source_at_number_set(heap, base, j, tmp)?;
        k = k.wrapping_add(1);
    }
    Ok(n)
}

#[allow(non_snake_case)]
pub(crate) fn insertions_sort_NeighListBySymmAndCanonRank(
    heap: &mut SourceHeap,
    base: NEIGH_LIST,
    nSymmRank: SourceMutPointer<AT_RANK>,
    nCanonRank: SourceMutPointer<AT_RANK>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:421 insertions_sort_NeighListBySymmAndCanonRank
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    void insertions_sort_NeighListBySymmAndCanonRank( NEIGH_LIST base,
                                                      const AT_RANK *nSymmRank,
                                                      const AT_RANK *nCanonRank )
    {
        AT_NUMB *i, *j, *pk, tmp;
        int  diff;
        int k, num = (int) *base++;
        for (k = 1, pk = base; k < num; k++, pk++)
        {
            for (j = ( i = pk ) + 1;
                 j > base &&    /*  always j > i */
                 ( 0 > ( diff = (int) nSymmRank[(int) *i] - (int) nSymmRank[(int) *j] ) ||
                 (!diff && nCanonRank[(int) *i] < nCanonRank[(int) *j]) ); /* djb-rwth: addressing LLVM warning */
                 j = i, i--)
            {
                tmp = *i;
                *i = *j;
                *j = tmp;
            }
        }

    }
    */
    // END INCHI C FUNCTION: insertions_sort_NeighListBySymmAndCanonRank

    let num = i32::from(source_at_number_get(heap, base, 0)?);
    let mut k = 1_i32;
    while k < num {
        let mut i = k;
        let mut j = k.wrapping_add(1);
        while j > 1 {
            let left = source_at_number_get(heap, base, i)?;
            let right = source_at_number_get(heap, base, j)?;
            let diff = i32::from(source_at_number_get(heap, nSymmRank, i32::from(left))?)
                - i32::from(source_at_number_get(heap, nSymmRank, i32::from(right))?);
            let should_swap = diff < 0
                || (diff == 0
                    && source_at_number_get(heap, nCanonRank, i32::from(left))?
                        < source_at_number_get(heap, nCanonRank, i32::from(right))?);
            if !should_swap {
                break;
            }
            source_at_number_set(heap, base, i, right)?;
            source_at_number_set(heap, base, j, left)?;
            j = i;
            i = i.wrapping_sub(1);
        }
        k = k.wrapping_add(1);
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn insertions_sort_AT_RANK(
    base: &mut [AT_RANK],
    num: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:375 insertions_sort_AT_RANK
    // INCHI✔️✔️: int insertions_sort_AT_RANK( AT_RANK *base, int num )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     AT_RANK *i, *j, *pk, tmp;
    // INCHI✔️✔️:     int  k, num_trans = 0;
    // INCHI✔️✔️:     for (k = 1, pk = base; k < num; k++, pk++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         for (j = ( i = pk ) + 1, tmp = *j; j > base && *i > tmp; j = i, i--)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             *j = *i;
    // INCHI✔️✔️:             num_trans++;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         *j = tmp;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return num_trans;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: insertions_sort_AT_RANK

    let count = usize::try_from(num.max(0)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if base.len() < count {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let mut num_trans = 0_i32;
    for k in 1..count {
        let temporary = base[k];
        let mut j = k;
        while j > 0 && base[j - 1] > temporary {
            base[j] = base[j - 1];
            num_trans = num_trans.wrapping_add(1);
            j -= 1;
        }
        base[j] = temporary;
    }
    Ok(num_trans)
}

pub(crate) fn inchi_qsort(
    base: &mut [u8],
    num: usize,
    width: usize,
    compare: &mut dyn FnMut(&[u8], &[u8]) -> Result<i32, SourceHeapError>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:66 inchi_qsort
    // INCHI✔️✔️: void inchi_qsort( void *pParam,
    // INCHI✔️✔️:                   void *base,
    // INCHI✔️✔️:                   size_t num,
    // INCHI✔️✔️:                   size_t width,
    // INCHI✔️✔️:                   int( *comp )( const void *, const void *, void * ) )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     char *lo, *hi;              /* ends of sub-array currently sorting */
    // INCHI✔️✔️:     char *mid;                  /* points to middle of subarray */
    // INCHI✔️✔️:     char *loguy, *higuy;        /* traveling pointers for partition step */
    // INCHI✔️✔️:     size_t size;                /* size of the sub-array */
    // INCHI✔️✔️:     char *lostk[STKSIZ], *histk[STKSIZ];
    // INCHI✔️✔️:     int stkptr;                 /* stack for saving sub-array to be processed */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (num < 2)
    // INCHI✔️✔️:         return;                 /* nothing to do */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     stkptr = 0;                 /* initialize stack */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     lo = (char *) base;
    // INCHI✔️✔️:     hi = (char *) base + width * ( num - 1 );        /* initialize limits */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* This entry point is for pseudo-recursion calling: setting
    // INCHI✔️✔️:     lo and hi and jumping to here is like recursion, but stkptr is
    // INCHI✔️✔️:     preserved, locals aren't, so we preserve stuff on the stack */
    // INCHI✔️✔️: recurse:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     size = ( hi - lo ) / width + 1;        /* number of el's to sort */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* First we pick a partitioning element.  The efficiency of the
    // INCHI✔️✔️:     algorithm demands that we find one that is approximately the median
    // INCHI✔️✔️:     of the values, but also that we select one fast.  We choose the
    // INCHI✔️✔️:     median of the first, middle, and last elements, to avoid bad
    // INCHI✔️✔️:     performance in the face of already sorted data, or data that is made
    // INCHI✔️✔️:     up of multiple sorted runs appended together.  Testing shows that a
    // INCHI✔️✔️:     median-of-three algorithm provides better performance than simply
    // INCHI✔️✔️:     picking the middle element for the latter case. */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     mid = lo + ( size / 2 ) * width;      /* find middle element */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Sort the first, middle, last elements into order */
    // INCHI✔️✔️:     if (comp( lo, mid, pParam ) > 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inchi_swap( lo, mid, width );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (comp( lo, hi, pParam ) > 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inchi_swap( lo, hi, width );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (comp( mid, hi, pParam ) > 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inchi_swap( mid, hi, width );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* We now wish to partition the array into three pieces, one consisting
    // INCHI✔️✔️:     of elements <= partition element, one of elements equal to the
    // INCHI✔️✔️:     partition element, and one of elements > than it.  This is done
    // INCHI✔️✔️:     below; comments indicate conditions established at every step. */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     loguy = lo;
    // INCHI✔️✔️:     higuy = hi;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Note that higuy decreases and loguy increases on every iteration,
    // INCHI✔️✔️:     so loop must terminate. */
    // INCHI✔️✔️:     for (;;)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* lo <= loguy < hi, lo < higuy <= hi,
    // INCHI✔️✔️:         A[i] <= A[mid] for lo <= i <= loguy,
    // INCHI✔️✔️:         A[i] > A[mid] for higuy <= i < hi,
    // INCHI✔️✔️:         A[hi] >= A[mid] */
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* The doubled loop is to avoid calling comp(mid,mid), since some
    // INCHI✔️✔️:         existing comparison funcs don't work when passed the same
    // INCHI✔️✔️:         value for both pointers. */
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (mid > loguy)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             do
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 loguy += width;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             while (loguy < mid && comp( loguy, mid, pParam ) <= 0);
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (mid <= loguy)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             do
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 loguy += width;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             while (loguy <= hi && comp( loguy, mid, pParam ) <= 0);
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         do
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             higuy -= width;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         while (higuy > mid && comp( higuy, mid, pParam ) > 0);
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (higuy < loguy)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         inchi_swap( loguy, higuy, width );
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (mid == higuy)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             mid = loguy;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     higuy += width;
    // INCHI✔️✔️:     if (mid < higuy)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         do
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             higuy -= width;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         while (higuy > mid && comp( higuy, mid, pParam ) == 0);
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (mid >= higuy)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         do
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             higuy -= width;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         while (higuy > lo && comp( higuy, mid, pParam ) == 0);
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (higuy - lo >= hi - loguy)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (lo < higuy)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             lostk[stkptr] = lo;
    // INCHI✔️✔️:             histk[stkptr] = higuy;
    // INCHI✔️✔️:             ++stkptr;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (loguy < hi)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             lo = loguy;
    // INCHI✔️✔️:             goto recurse;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (loguy < hi)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             lostk[stkptr] = loguy;
    // INCHI✔️✔️:             histk[stkptr] = hi;
    // INCHI✔️✔️:             ++stkptr;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (lo < higuy)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             hi = higuy;
    // INCHI✔️✔️:             goto recurse;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     --stkptr;
    // INCHI✔️✔️:     if (stkptr >= 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         lo = lostk[stkptr];
    // INCHI✔️✔️:         hi = histk[stkptr];
    // INCHI✔️✔️:         goto recurse;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: inchi_qsort
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: inchi_qsort
    // INCHI✔️✔️: #define STKSIZ (8*sizeof(void*) - 2)
    // INCHI✔️✔️: GCC/Linux LP64: STKSIZ == 62
    // END INCHI ACTIVE MACRO CONFIGURATION: inchi_qsort

    if num < 2 {
        return Ok(());
    }
    if width == 0 {
        return Err(SourceHeapError::UnsupportedSourceBehavior);
    }
    let required = num
        .checked_mul(width)
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    if required > base.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let width = isize::try_from(width).map_err(|_| SourceHeapError::PointerOffsetOverflow)?;
    let mut lo = 0_isize;
    let mut hi = isize::try_from(required - width as usize)
        .map_err(|_| SourceHeapError::PointerOffsetOverflow)?;
    let mut low_stack = [0_isize; 62];
    let mut high_stack = [0_isize; 62];
    let mut stack_pointer = 0_usize;

    let compare_at = |bytes: &[u8],
                      first: isize,
                      second: isize,
                      compare: &mut dyn FnMut(&[u8], &[u8]) -> Result<i32, SourceHeapError>|
     -> Result<i32, SourceHeapError> {
        let first = usize::try_from(first).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let second = usize::try_from(second).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let width = usize::try_from(width).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let first_slice = bytes
            .get(first..first + width)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let second_slice = bytes
            .get(second..second + width)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        compare(first_slice, second_slice)
    };

    loop {
        let size = (hi - lo) / width + 1;
        let mut mid = lo + (size / 2) * width;
        if compare_at(base, lo, mid, compare)? > 0 {
            inchi_swap(base, lo as usize, mid as usize, width as usize)?;
        }
        if compare_at(base, lo, hi, compare)? > 0 {
            inchi_swap(base, lo as usize, hi as usize, width as usize)?;
        }
        if compare_at(base, mid, hi, compare)? > 0 {
            inchi_swap(base, mid as usize, hi as usize, width as usize)?;
        }

        let mut low_guy = lo;
        let mut high_guy = hi;
        loop {
            if mid > low_guy {
                loop {
                    low_guy += width;
                    if !(low_guy < mid && compare_at(base, low_guy, mid, compare)? <= 0) {
                        break;
                    }
                }
            }
            if mid <= low_guy {
                loop {
                    low_guy += width;
                    if !(low_guy <= hi && compare_at(base, low_guy, mid, compare)? <= 0) {
                        break;
                    }
                }
            }
            loop {
                high_guy -= width;
                if !(high_guy > mid && compare_at(base, high_guy, mid, compare)? > 0) {
                    break;
                }
            }
            if high_guy < low_guy {
                break;
            }
            inchi_swap(base, low_guy as usize, high_guy as usize, width as usize)?;
            if mid == high_guy {
                mid = low_guy;
            }
        }

        high_guy += width;
        if mid < high_guy {
            loop {
                high_guy -= width;
                if !(high_guy > mid && compare_at(base, high_guy, mid, compare)? == 0) {
                    break;
                }
            }
        }
        if mid >= high_guy {
            loop {
                high_guy -= width;
                if !(high_guy > lo && compare_at(base, high_guy, mid, compare)? == 0) {
                    break;
                }
            }
        }

        let mut recurse = false;
        if high_guy - lo >= hi - low_guy {
            if lo < high_guy {
                low_stack[stack_pointer] = lo;
                high_stack[stack_pointer] = high_guy;
                stack_pointer += 1;
            }
            if low_guy < hi {
                lo = low_guy;
                recurse = true;
            }
        } else {
            if low_guy < hi {
                low_stack[stack_pointer] = low_guy;
                high_stack[stack_pointer] = hi;
                stack_pointer += 1;
            }
            if lo < high_guy {
                hi = high_guy;
                recurse = true;
            }
        }
        if recurse {
            continue;
        }
        if stack_pointer == 0 {
            return Ok(());
        }
        stack_pointer -= 1;
        lo = low_stack[stack_pointer];
        hi = high_stack[stack_pointer];
    }
}

pub(crate) fn inchi_swap(
    bytes: &mut [u8],
    first: usize,
    second: usize,
    width: usize,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:286 inchi_swap
    // INCHI✔️✔️: void inchi_swap( char *a, char *b, size_t width )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     char tmp;
    // INCHI✔️✔️:     if (a != b)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         while (width--)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             tmp = *a;
    // INCHI✔️✔️:             *a++ = *b;
    // INCHI✔️✔️:             *b++ = tmp;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: inchi_swap

    if first == second {
        return Ok(());
    }
    first
        .checked_add(width)
        .filter(|end| *end <= bytes.len())
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    second
        .checked_add(width)
        .filter(|end| *end <= bytes.len())
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for offset in 0..width {
        let temporary = bytes[first + offset];
        bytes[first + offset] = bytes[second + offset];
        bytes[second + offset] = temporary;
    }
    Ok(())
}

pub(crate) fn insertions_sort(
    base: &mut [u8],
    number: usize,
    width: usize,
    compare: &mut dyn FnMut(&[u8], &[u8]) -> Result<i32, SourceHeapError>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:304 insertions_sort
    // INCHI✔️✔️: int insertions_sort( void *pCG,
    // INCHI✔️✔️:                      void *base,
    // INCHI✔️✔️:                      size_t num, size_t width,
    // INCHI✔️✔️:                      int( *compare )( const void *, const void *, void * ) ) /* djb-rwth: types of variables are sufficient */
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     char *i, *j, *pk = (char*) base;
    // INCHI✔️✔️:     int  num_trans = 0;
    // INCHI✔️✔️:     size_t k;
    // INCHI✔️✔️:     for (k = 1; k < num; k++, pk += width)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /*for( i = pk, j = pk + width; j > (char*)base && (*compare)(i,j) > 0; j=i, i -= width )*/
    // INCHI✔️✔️:         for (i = j = pk + width;
    // INCHI✔️✔️:              j > ( char* )base && ( i -= width, ( *compare )( i, j, pCG ) ) > 0;
    // INCHI✔️✔️:              j = i)        /* changed to keep BoundsChecker happy 2007-09-24 DT */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             inchi_swap( i, j, width );
    // INCHI✔️✔️:             num_trans++;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return num_trans;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: insertions_sort

    let required = number
        .checked_mul(width)
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    if required > base.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let mut transactions = 0_i32;
    for element in 1..number {
        let mut right = element * width;
        while right > 0 {
            let left = right - width;
            let ordering = {
                let first = &base[left..left + width];
                let second = &base[right..right + width];
                compare(first, second)?
            };
            if ordering <= 0 {
                break;
            }
            inchi_swap(base, left, right, width)?;
            transactions = transactions.wrapping_add(1);
            right = left;
        }
    }
    Ok(transactions)
}

#[allow(non_snake_case)]
pub(crate) fn CreateNeighListFromLinearCT(
    heap: &mut SourceHeap,
    linear_ct: SourceConstPointer<AT_NUMB>,
    length_ct: i32,
    number_of_atoms: i32,
) -> Result<SourceMutPointer<NEIGH_LIST>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:701 CreateNeighListFromLinearCT
    // INCHI✔️✔️: NEIGH_LIST *CreateNeighListFromLinearCT( AT_NUMB *LinearCT, int nLenCT, int num_atoms )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     /* Atom numbers in LinearCT are canonical numbers
    // INCHI✔️✔️:      * order: parent[i] > neigh[i][0] < neigh[i][1]...<neigh[i][n] < parent[i+1] > neigh[i+1][0] < ...
    // INCHI✔️✔️:      *        parent[i] < parent[i+1]
    // INCHI✔️✔️:      */
    // INCHI✔️✔️:     int i, j;
    // INCHI✔️✔️:     S_CHAR     *valence = NULL;
    // INCHI✔️✔️:     NEIGH_LIST *pp = NULL;
    // INCHI✔️✔️:     AT_NUMB    *pAtList = NULL;
    // INCHI✔️✔️:     AT_RANK     n_vertex, n_neigh;
    // INCHI✔️✔️:     int err = 1, num_bonds;
    // INCHI✔️✔️:     int length, start;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if ((int) LinearCT[0] > num_atoms)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         goto exit_function;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     valence = (S_CHAR*)inchi_calloc((long long)num_atoms + 1, sizeof(valence[0]));
    // INCHI✔️✔️:     if (!valence) /* djb-rwth: cast operator added */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         goto exit_function;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i = 1, num_bonds = 0, n_vertex = LinearCT[0]; i < nLenCT; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (( n_neigh = LinearCT[i] ) < n_vertex)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             valence[n_neigh] ++;
    // INCHI✔️✔️:             valence[n_vertex] ++;
    // INCHI✔️✔️:             num_bonds += 2;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if ((int) ( n_vertex = n_neigh ) > num_atoms)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 goto exit_function;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if ((int) n_vertex != num_atoms)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         goto exit_function;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     length = num_bonds + num_atoms + 1;
    // INCHI✔️✔️:     pp = (NEIGH_LIST*)inchi_calloc(((long long)num_atoms + 1), sizeof(NEIGH_LIST));
    // INCHI✔️✔️:     pAtList = (AT_NUMB*)inchi_malloc(length * sizeof(AT_NUMB));
    // INCHI✔️✔️:     if (pp && pAtList) /* djb-rwth: cast operator added; addressing LLVM warning */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /*  Create empty connection table */
    // INCHI✔️✔️:         for (i = 1, length = 0; i <= num_atoms; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             start = length;
    // INCHI✔️✔️:             length += ( valence[i] + 1 );
    // INCHI✔️✔️:             pp[i - 1] = pAtList + start;
    // INCHI✔️✔️:             pp[i - 1][0] = 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         /*  Fill out the CT */
    // INCHI✔️✔️:         for (i = 1, n_vertex = LinearCT[0] - 1; i < nLenCT; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (( n_neigh = LinearCT[i] - 1 ) < n_vertex)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 /*  Vertex - neighbor connection */
    // INCHI✔️✔️:                 j = (int) ( ++pp[(int) n_vertex][0] );
    // INCHI✔️✔️:                 pp[(int) n_vertex][j] = n_neigh;
    // INCHI✔️✔️:                 /*  neighbor - vertex connection */
    // INCHI✔️✔️:                 j = (int) ( ++pp[(int) n_neigh][0] );
    // INCHI✔️✔️:                 pp[(int) n_neigh][j] = n_vertex;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if ((int) ( n_vertex = n_neigh ) >= num_atoms)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     goto exit_function;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         err = 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️: exit_function:
    // INCHI✔️✔️:     if (valence)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inchi_free( valence );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (err) /* djb-rwth: ignoring LLVM warning */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (pAtList)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             inchi_free( pAtList );
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (pp)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             inchi_free( pp );
    // INCHI✔️✔️:             pp = NULL;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return pp; /* djb-rwth: ignoring LLVM warning: since a pointer is returned, memory should be freed in a function which calls *CreateNeighListFromLinearCT */
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: CreateNeighListFromLinearCT

    let available_ct_length = heap.slice(linear_ct)?.len();
    let first = *heap
        .slice(linear_ct)?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if i32::from(first) > number_of_atoms {
        return Ok(SourceMutPointer::null());
    }
    let atom_count =
        usize::try_from(number_of_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let valence = match inchi_calloc::<S_CHAR>(heap, atom_count as u64 + 1, 1) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return Ok(SourceMutPointer::null()),
        Err(error) => return Err(error),
    };
    let ct_length = if length_ct <= 0 {
        0
    } else {
        usize::try_from(length_ct).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
    };
    if available_ct_length < ct_length {
        inchi_free(heap, valence)?;
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let mut vertex = first;
    let mut number_of_bonds = 0_i32;
    let mut valid = true;
    for index in 1..ct_length {
        let neighbor = heap.slice(linear_ct)?[index];
        if neighbor < vertex {
            let values = heap.slice_mut(valence)?;
            let neighbor_index = usize::from(neighbor);
            let vertex_index = usize::from(vertex);
            if neighbor_index >= values.len() || vertex_index >= values.len() {
                valid = false;
                break;
            }
            values[neighbor_index] = values[neighbor_index]
                .checked_add(1)
                .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
            values[vertex_index] = values[vertex_index]
                .checked_add(1)
                .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
            number_of_bonds = number_of_bonds.wrapping_add(2);
        } else {
            vertex = neighbor;
            if i32::from(vertex) > number_of_atoms {
                valid = false;
                break;
            }
        }
    }
    if !valid || i32::from(vertex) != number_of_atoms {
        inchi_free(heap, valence)?;
        return Ok(SourceMutPointer::null());
    }
    let list_length = number_of_bonds
        .wrapping_add(number_of_atoms)
        .wrapping_add(1);
    let pointer_list = match inchi_calloc::<NEIGH_LIST>(heap, atom_count as u64 + 1, 1) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
        Err(error) => {
            inchi_free(heap, valence)?;
            return Err(error);
        }
    };
    let atom_list = match usize::try_from(list_length) {
        Ok(length) => match heap.allocate(vec![0_u16; length]) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => {
                inchi_free(heap, valence)?;
                inchi_free(heap, pointer_list)?;
                return Err(error);
            }
        },
        Err(_) => {
            inchi_free(heap, valence)?;
            inchi_free(heap, pointer_list)?;
            return Ok(SourceMutPointer::null());
        }
    };

    if pointer_list.is_null() || atom_list.is_null() {
        inchi_free(heap, valence)?;
        inchi_free(heap, atom_list)?;
        inchi_free(heap, pointer_list)?;
        return Ok(SourceMutPointer::null());
    }

    let mut position = 0_usize;
    for atom in 1..=atom_count {
        let start = position;
        let atom_valence = heap.slice(valence.as_const())?[atom];
        position = position.wrapping_add((i32::from(atom_valence) + 1) as usize);
        let pointer = atom_list
            .offset(i64::try_from(start).map_err(|_| SourceHeapError::PointerOffsetOverflow)?)?;
        heap.slice_mut(pointer_list)?[atom - 1] = pointer;
        heap.slice_mut(pointer)?[0] = 0;
    }
    vertex = first.wrapping_sub(1);
    for index in 1..ct_length {
        let raw_neighbor = heap.slice(linear_ct)?[index];
        let neighbor = raw_neighbor.wrapping_sub(1);
        if neighbor < vertex {
            for (owner, adjacent) in [(vertex, neighbor), (neighbor, vertex)] {
                let owner_pointer = heap.slice(pointer_list.as_const())?[usize::from(owner)];
                let owner_list = heap.slice_mut(owner_pointer)?;
                owner_list[0] = owner_list[0].wrapping_add(1);
                let entry = usize::from(owner_list[0]);
                owner_list[entry] = adjacent;
            }
        } else {
            vertex = neighbor;
            if i32::from(vertex) >= number_of_atoms {
                inchi_free(heap, valence)?;
                inchi_free(heap, atom_list)?;
                inchi_free(heap, pointer_list)?;
                return Ok(SourceMutPointer::null());
            }
        }
    }
    inchi_free(heap, valence)?;
    Ok(pointer_list)
}

#[allow(non_snake_case)]
pub(crate) fn CreateNeighList(
    heap: &mut SourceHeap,
    num_atoms: i32,
    num_at_tg: i32,
    at: SourceConstPointer<sp_ATOM>,
    bDoubleBondSquare: i32,
    t_group_info: Option<&T_GROUP_INFO>,
) -> Result<SourceMutPointer<NEIGH_LIST>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:810 CreateNeighList
    // INCHI✔️✔️: NEIGH_LIST *CreateNeighList( int num_atoms,
    // INCHI✔️✔️:                              int num_at_tg,
    // INCHI✔️✔️:                              sp_ATOM* at,
    // INCHI✔️✔️:                              int bDoubleBondSquare,
    // INCHI✔️✔️:                              T_GROUP_INFO *t_group_info )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     /*  +1 to add NULL termination */
    // INCHI✔️✔️:     NEIGH_LIST *pp = (NEIGH_LIST *) inchi_calloc( ( (long long)num_at_tg + 1 ), sizeof( NEIGH_LIST ) ); /* djb-rwth: cast operator added */
    // INCHI✔️✔️:     T_GROUP   *t_group = NULL;
    // INCHI✔️✔️:     AT_NUMB   *nEndpointAtomNumber = NULL;
    // INCHI✔️✔️:     int        num_t_groups = 0;
    // INCHI✔️✔️:     int        nFirstEndpointAtNoPos;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     AT_NUMB   *pAtList = NULL;
    // INCHI✔️✔️:     int        length, start, val, i, j;
    // INCHI✔️✔️:     if (pp)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (num_at_tg > num_atoms)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             t_group = t_group_info->t_group;
    // INCHI✔️✔️:             num_t_groups = t_group_info->num_t_groups;
    // INCHI✔️✔️:             nEndpointAtomNumber = t_group_info->nEndpointAtomNumber;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (!bDoubleBondSquare)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             for (i = 0, length = 0; i < num_atoms; i++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 length += (int) at[i].valence + ( num_t_groups && at[i].endpoint );
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             length += num_atoms;
    // INCHI✔️✔️:             for (i = 0; i < num_t_groups; i++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 length += (int) t_group[i].nNumEndpoints;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             length += num_t_groups;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             for (i = 0, length = 0; i < num_atoms; i++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 val = (int) at[i].valence;
    // INCHI✔️✔️:                 for (j = 0; j < val; j++)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     length += 1 + ( bDoubleBondSquare && BOND_DOUBLE == at[i].bond_type[j] );
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 length += ( num_t_groups && at[i].endpoint );
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             length += num_atoms;
    // INCHI✔️✔️:             for (i = 0; i < num_t_groups; i++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 length += (int) t_group[i].nNumEndpoints;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             length += num_t_groups;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         length++; /*  +1 to save number of neighbors */
    // INCHI✔️✔️:         pAtList = (AT_NUMB*)inchi_malloc(length * sizeof(AT_NUMB));
    // INCHI✔️✔️:         if (pAtList) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (!bDoubleBondSquare)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 for (i = 0, length = 0; i < num_atoms; i++)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     val = at[i].valence;
    // INCHI✔️✔️:                     start = length++;
    // INCHI✔️✔️:                     for (j = 0; j < val; j++)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         pAtList[length++] = at[i].neighbor[j];
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     /*  add endpoint */
    // INCHI✔️✔️:                     if (num_t_groups && at[i].endpoint)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         pAtList[length++] = num_atoms + (int) at[i].endpoint - 1;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     pAtList[start] = length - start - 1;  /*  number of neighbors before the list of neighbors */
    // INCHI✔️✔️:                     pp[i] = pAtList + start;              /*  pointer to the <num.neigh.><list of neigh> */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 for (i = 0, length = 0; i < num_atoms; i++)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     val = at[i].valence;
    // INCHI✔️✔️:                     start = length++;
    // INCHI✔️✔️:                     for (j = 0; j < val; j++)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         pAtList[length++] = at[i].neighbor[j]; /* djb-rwth: buffer overrun avoided implicitly */
    // INCHI✔️✔️:                         if (bDoubleBondSquare && BOND_DOUBLE == at[i].bond_type[j])
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             pAtList[length++] = at[i].neighbor[j]; /*  a list of neighbor orig. numbers */
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     /*  Add endpoint */
    // INCHI✔️✔️:                     if (num_t_groups && at[i].endpoint)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         pAtList[length++] = num_atoms + (int) at[i].endpoint - 1;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     pAtList[start] = length - start - 1;  /*  number of neighbors before the list of neighbors */
    // INCHI✔️✔️:                     pp[i] = pAtList + start;              /*  pointer to the <num.neigh.><list of neigh> */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:
    // INCHI✔️✔️:             /*  Add t-groups */
    // INCHI✔️✔️:             for (i = 0; i < num_t_groups; i++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 val = (int) t_group[i].nNumEndpoints;
    // INCHI✔️✔️:                 start = length++;
    // INCHI✔️✔️:                 nFirstEndpointAtNoPos = (int) t_group[i].nFirstEndpointAtNoPos;
    // INCHI✔️✔️:                 for (j = 0; j < val; j++)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     pAtList[length++] = nEndpointAtomNumber[nFirstEndpointAtNoPos + j];
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 pAtList[start] = length - start - 1;  /*  number of neighbors before the list of neighbors */
    // INCHI✔️✔️:                 pp[num_atoms + i] = pAtList + start;    /*  pointer to the <num.neigh.><list of neigh> */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             inchi_free(pAtList); /* djb-rwth: fixing coverity ID #499598 */
    // INCHI✔️✔️:             inchi_free( pp );
    // INCHI✔️✔️:             return NULL;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     } /* djb-rwth: ignoring LLVM warning */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* djb-rwth: fixing coverity ID #499598 -- pp uses pAtList values */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return pp; /* djb-rwth: ignoring LLVM warning: since a pointer is returned, memory should be freed in a function which calls *CreateNeighList */
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: CreateNeighList

    let pointer_count = i64::from(num_at_tg)
        .checked_add(1)
        .and_then(|value| u64::try_from(value).ok())
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let pointer_list = match inchi_calloc::<NEIGH_LIST>(heap, pointer_count, 8) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return Ok(SourceMutPointer::null()),
        Err(error) => return Err(error),
    };
    let atom_count = match usize::try_from(num_atoms) {
        Ok(count) => count,
        Err(_) => {
            inchi_free(heap, pointer_list)?;
            return Err(SourceHeapError::SourceIntegerOverflow);
        }
    };

    let (num_t_groups, group_pointer, endpoint_pointer) = if num_at_tg > num_atoms {
        let Some(info) = t_group_info else {
            inchi_free(heap, pointer_list)?;
            return Err(SourceHeapError::NullPointer);
        };
        (info.num_t_groups, info.t_group, info.nEndpointAtomNumber)
    } else {
        (0, SourceMutPointer::null(), SourceMutPointer::null())
    };
    let group_count =
        usize::try_from(num_t_groups.max(0)).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;

    let mut list_length = i64::from(num_atoms) + i64::from(num_t_groups);
    for index in 0..atom_count {
        let atom = heap
            .slice(at)?
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let valence = i32::from(atom.valence);
        if valence < 0 {
            inchi_free(heap, pointer_list)?;
            return Err(SourceHeapError::UnsupportedSourceBehavior);
        }
        let valence =
            usize::try_from(valence).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if valence > atom.neighbor.len() {
            inchi_free(heap, pointer_list)?;
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        if bDoubleBondSquare == 0 {
            list_length += valence as i64;
        } else {
            for bond in 0..valence {
                list_length += 1 + i64::from(atom.bond_type[bond] == BOND_DOUBLE as u8);
            }
        }
        list_length += i64::from(num_t_groups != 0 && atom.endpoint != 0);
    }
    for index in 0..group_count {
        let group = heap
            .slice(group_pointer.as_const())?
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        list_length += i64::from(group.nNumEndpoints);
    }
    list_length = list_length
        .checked_add(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let list_length = match usize::try_from(list_length) {
        Ok(length) => length,
        Err(_) => {
            inchi_free(heap, pointer_list)?;
            return Err(SourceHeapError::SourceIntegerOverflow);
        }
    };
    let atom_list = match heap.allocate(vec![0_u16; list_length]) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => {
            inchi_free(heap, pointer_list)?;
            return Ok(SourceMutPointer::null());
        }
        Err(error) => {
            inchi_free(heap, pointer_list)?;
            return Err(error);
        }
    };

    let mut position = 0_usize;
    for index in 0..atom_count {
        let (valence, neighbors, bond_types, endpoint) = {
            let atom = heap
                .slice(at)?
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            (
                usize::try_from(atom.valence)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                atom.neighbor,
                atom.bond_type,
                atom.endpoint,
            )
        };
        let start = position;
        position += 1;
        for bond in 0..valence {
            heap.slice_mut(atom_list)?[position] = neighbors[bond];
            position += 1;
            if bDoubleBondSquare != 0 && bond_types[bond] == BOND_DOUBLE as u8 {
                heap.slice_mut(atom_list)?[position] = neighbors[bond];
                position += 1;
            }
        }
        if num_t_groups != 0 && endpoint != 0 {
            heap.slice_mut(atom_list)?[position] =
                num_atoms.wrapping_add(i32::from(endpoint)).wrapping_sub(1) as AT_NUMB;
            position += 1;
        }
        heap.slice_mut(atom_list)?[start] = position.wrapping_sub(start).wrapping_sub(1) as AT_NUMB;
        heap.slice_mut(pointer_list)?[index] = atom_list.offset(start as i64)?;
    }

    for index in 0..group_count {
        let (valence, first_endpoint) = {
            let group = heap
                .slice(group_pointer.as_const())?
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            (
                usize::from(group.nNumEndpoints),
                usize::from(group.nFirstEndpointAtNoPos),
            )
        };
        let start = position;
        position += 1;
        for endpoint in 0..valence {
            let value = *heap
                .slice(endpoint_pointer.as_const())?
                .get(first_endpoint + endpoint)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            heap.slice_mut(atom_list)?[position] = value;
            position += 1;
        }
        heap.slice_mut(atom_list)?[start] = position.wrapping_sub(start).wrapping_sub(1) as AT_NUMB;
        let pointer_index = atom_count + index;
        let pointers = heap.slice_mut(pointer_list)?;
        if pointer_index >= pointers.len() {
            inchi_free(heap, atom_list)?;
            inchi_free(heap, pointer_list)?;
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        pointers[pointer_index] = atom_list.offset(start as i64)?;
    }

    Ok(pointer_list)
}

#[allow(non_snake_case)]
pub(crate) fn FreeNeighList(
    heap: &mut SourceHeap,
    pointer_list: SourceMutPointer<NEIGH_LIST>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:941 FreeNeighList
    // INCHI✔️✔️: void FreeNeighList( NEIGH_LIST *pp )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     if (pp)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (pp[0])
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             inchi_free( pp[0] );
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         inchi_free( pp );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: FreeNeighList

    if !pointer_list.is_null() {
        let atom_list = heap.slice(pointer_list.as_const())?[0];
        if !atom_list.is_null() {
            inchi_free(heap, atom_list)?;
        }
        inchi_free(heap, pointer_list)?;
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn BreakAllTies(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    num_atoms: i32,
    num_max: i32,
    pRankStack: SourceMutPointer<SourceMutPointer<AT_RANK>>,
    NeighList: SourceMutPointer<NEIGH_LIST>,
    nTempRank: SourceMutPointer<AT_RANK>,
    pCS: &mut CANON_STAT,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:956 BreakAllTies
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int BreakAllTies( CANON_GLOBALS *pCG,
                      int num_atoms,
                      int num_max,
                      AT_RANK **pRankStack,
                      NEIGH_LIST *NeighList,
                      AT_RANK *nTempRank,
                      CANON_STAT *pCS )
    {
        int i, nRet = -1, nNumRanks = 1 /* value does not matter*/;

        AT_RANK *nPrevRank = *pRankStack++;
        AT_RANK *nPrevAtomNumber = *pRankStack++;

        AT_RANK *nNewRank = NULL;
        AT_RANK *nNewAtomNumber = NULL;

        if (!pRankStack[0])
        {
            pRankStack[0] = (AT_RANK *) inchi_malloc( num_max * sizeof( *nNewRank ) );
        }
        if (!pRankStack[1])
        {
            pRankStack[1] = (AT_RANK *) inchi_malloc( num_max * sizeof( *nNewAtomNumber ) );
        }
        if (!pRankStack[0] || !pRankStack[1])
        {
            return CT_OUT_OF_RAM;  /*   <BRKPT> */
        }
        nNewRank = pRankStack[0];
        nNewAtomNumber = pRankStack[1];

        if (nNewRank && nNewAtomNumber)
        {
            memcpy(nNewAtomNumber, nPrevAtomNumber, num_atoms * sizeof(nNewAtomNumber[0]));
            memcpy(nNewRank, nPrevRank, num_atoms * sizeof(nNewRank[0]));
            for (i = 1, nRet = 0; i < num_atoms; i++)
            {
                /*  12-12-2001: replaced Prev... with New... */
                if (nNewRank[(int) nNewAtomNumber[i - 1]] == nNewRank[(int) nNewAtomNumber[i]])
                {
                    nNewRank[nNewAtomNumber[i - 1]] = (AT_RANK) i;
                    nNumRanks = DifferentiateRanks2( pCG, num_atoms, NeighList,
                                                     nNumRanks, nNewRank, nTempRank,
                                                     nNewAtomNumber,
                                                     &pCS->lNumNeighListIter, 1 );
                    pCS->lNumBreakTies++;
                    nRet++;
                }
            }
        }

        return nRet;
    }
    */
    // END INCHI C FUNCTION: BreakAllTies

    let count = usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let capacity = usize::try_from(num_max).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let (previous_rank, previous_atom_number) = {
        let stack = heap.slice(pRankStack.as_const())?;
        (
            *stack.first().ok_or(SourceHeapError::PointerOutOfBounds)?,
            *stack.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
    };

    for slot in 2..=3_usize {
        let current = *heap
            .slice(pRankStack.as_const())?
            .get(slot)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if current.is_null() {
            let allocated = match heap.allocate(vec![0_u16; capacity]) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                Err(error) => return Err(error),
            };
            heap.slice_mut(pRankStack)?[slot] = allocated;
        }
    }

    let (new_rank, new_atom_number) = {
        let stack = heap.slice(pRankStack.as_const())?;
        (
            *stack.get(2).ok_or(SourceHeapError::PointerOutOfBounds)?,
            *stack.get(3).ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
    };
    if new_rank.is_null() || new_atom_number.is_null() {
        return Ok(CT_OUT_OF_RAM);
    }

    let copy_prefix = |heap: &mut SourceHeap,
                       destination: SourceMutPointer<AT_RANK>,
                       source: SourceMutPointer<AT_RANK>|
     -> Result<(), SourceHeapError> {
        if destination == source {
            if heap.slice(source.as_const())?.len() < count {
                return Err(SourceHeapError::PointerOutOfBounds);
            }
            return Ok(());
        }
        heap.with_slice_mut_and_heap(destination, |destination, heap| {
            let source = heap
                .slice(source.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let destination = destination
                .get_mut(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            destination.copy_from_slice(source);
            Ok(())
        })
    };
    copy_prefix(heap, new_atom_number, previous_atom_number)?;
    copy_prefix(heap, new_rank, previous_rank)?;

    let mut result = 0_i32;
    let mut number_of_ranks = 1_i32;
    let mut index = 1_i32;
    while index < num_atoms {
        let (previous_atom, current_atom) = {
            let atom_numbers = heap.slice(new_atom_number.as_const())?;
            let previous_index = usize::try_from(index.wrapping_sub(1))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let current_index =
                usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            (
                *atom_numbers
                    .get(previous_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                *atom_numbers
                    .get(current_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
        };
        let tied = {
            let ranks = heap.slice(new_rank.as_const())?;
            *ranks
                .get(usize::from(previous_atom))
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                == *ranks
                    .get(usize::from(current_atom))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
        };
        if tied {
            heap.slice_mut(new_rank)?[usize::from(previous_atom)] = index as AT_RANK;
            number_of_ranks = DifferentiateRanks2(
                heap,
                pCG,
                num_atoms,
                NeighList,
                number_of_ranks,
                new_rank,
                nTempRank,
                new_atom_number,
                &mut pCS.lNumNeighListIter,
                1,
            )?;
            pCS.lNumBreakTies = pCS.lNumBreakTies.wrapping_add(1);
            result = result.wrapping_add(1);
        }
        index = index.wrapping_add(1);
    }
    Ok(result)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn iisort(
    heap: &mut SourceHeap,
    list: SourceMutPointer<i32>,
    number: i32,
) -> Result<SourceMutPointer<i32>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:1014 iisort
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int * iisort( int *list, int num )
{
    int i;
    for (i = 1; i < num; i++)
    {
        int tmp = list[i];
        int j = i - 1;
        while (j >= 0 && list[j] > tmp)
        {
            list[j + 1] = list[j];
            j--;
        }
        list[j + 1] = tmp;
    }

    return list;
}
    */
    // END INCHI C FUNCTION: iisort
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: iisort
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; signed int comparison and stable insertion order are preserved.
    // INCHI✔️❌: Rust performs the same O(num^2) in-place insertion sort without a temporary allocation or clone.
    // INCHI✔️❌: SourceHeap's BTreeMap allocation lookup is materially more expensive than the source's direct pointer access.
    // END INCHI ACTIVE MACRO CONFIGURATION: iisort

    if number <= 1 {
        return Ok(list);
    }
    let number =
        usize::try_from(number).map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    let values = heap.slice_mut(list)?;
    if number > values.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    for index in 1..number {
        let temporary = values[index];
        let mut previous = index;
        while previous > 0 && values[previous - 1] > temporary {
            values[previous] = values[previous - 1];
            previous -= 1;
        }
        values[previous] = temporary;
    }
    Ok(list)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{ATOM_INVARIANT2, T_GROUP};

    #[test]
    fn source_port__ichisort__breakallties__line_956() {
        fn empty_neighbors(heap: &mut SourceHeap, count: usize) -> SourceMutPointer<NEIGH_LIST> {
            let atoms = heap
                .allocate_model_storage(vec![sp_ATOM::default(); count])
                .unwrap();
            CreateNeighList(heap, count as i32, count as i32, atoms.as_const(), 0, None).unwrap()
        }

        let mut heap = SourceHeap::default();
        let previous_rank = heap.allocate_model_storage(vec![1_u16, 2, 3]).unwrap();
        let previous_order = heap.allocate_model_storage(vec![0_u16, 1, 2]).unwrap();
        let stack = heap
            .allocate_model_storage(vec![
                previous_rank,
                previous_order,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ])
            .unwrap();
        let neighbors = empty_neighbors(&mut heap, 3);
        let temporary = heap.allocate_model_storage(vec![0_u16; 3]).unwrap();
        let mut stats = CANON_STAT {
            lNumBreakTies: 7,
            lNumNeighListIter: 11,
            ..CANON_STAT::default()
        };
        assert_eq!(
            BreakAllTies(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                3,
                3,
                stack,
                neighbors,
                temporary,
                &mut stats,
            ),
            Ok(0)
        );
        let stack_values = heap.slice(stack.as_const()).unwrap();
        assert_eq!(heap.slice(stack_values[2].as_const()).unwrap(), &[1, 2, 3]);
        assert_eq!(heap.slice(stack_values[3].as_const()).unwrap(), &[0, 1, 2]);
        assert_eq!((stats.lNumBreakTies, stats.lNumNeighListIter), (7, 11));

        let mut tied_heap = SourceHeap::default();
        let previous_rank = tied_heap.allocate_model_storage(vec![2_u16, 2]).unwrap();
        let previous_order = tied_heap.allocate_model_storage(vec![0_u16, 1]).unwrap();
        let new_rank = tied_heap
            .allocate_model_storage(vec![99_u16, 98, 97])
            .unwrap();
        let new_order = tied_heap
            .allocate_model_storage(vec![96_u16, 95, 94])
            .unwrap();
        let stack = tied_heap
            .allocate_model_storage(vec![previous_rank, previous_order, new_rank, new_order])
            .unwrap();
        let neighbors = empty_neighbors(&mut tied_heap, 2);
        let temporary = tied_heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
        let mut stats = CANON_STAT::default();
        assert_eq!(
            BreakAllTies(
                &mut tied_heap,
                &mut CANON_GLOBALS::default(),
                2,
                3,
                stack,
                neighbors,
                temporary,
                &mut stats,
            ),
            Ok(1)
        );
        assert_eq!(tied_heap.slice(new_rank.as_const()).unwrap(), &[1, 2, 97]);
        assert_eq!(tied_heap.slice(new_order.as_const()).unwrap(), &[0, 1, 94]);
        assert_eq!(stats.lNumBreakTies, 1);
        assert!(stats.lNumNeighListIter > 0);

        let mut first_failure = SourceHeap::default();
        let first_rank = first_failure.allocate_model_storage(vec![1_u16]).unwrap();
        let first_order = first_failure.allocate_model_storage(vec![0_u16]).unwrap();
        let first_stack = first_failure
            .allocate_model_storage(vec![
                first_rank,
                first_order,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ])
            .unwrap();
        first_failure.fail_after_allocations(0);
        assert_eq!(
            BreakAllTies(
                &mut first_failure,
                &mut CANON_GLOBALS::default(),
                1,
                1,
                first_stack,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut CANON_STAT::default(),
            ),
            Ok(CT_OUT_OF_RAM)
        );
        assert!(first_failure.slice(first_stack.as_const()).unwrap()[2].is_null());
        assert!(!first_failure.slice(first_stack.as_const()).unwrap()[3].is_null());

        let mut second_failure = SourceHeap::default();
        let second_rank = second_failure.allocate_model_storage(vec![1_u16]).unwrap();
        let second_order = second_failure.allocate_model_storage(vec![0_u16]).unwrap();
        let second_stack = second_failure
            .allocate_model_storage(vec![
                second_rank,
                second_order,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ])
            .unwrap();
        second_failure.fail_after_allocations(1);
        assert_eq!(
            BreakAllTies(
                &mut second_failure,
                &mut CANON_GLOBALS::default(),
                1,
                1,
                second_stack,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut CANON_STAT::default(),
            ),
            Ok(CT_OUT_OF_RAM)
        );
        assert!(!second_failure.slice(second_stack.as_const()).unwrap()[2].is_null());
        assert!(second_failure.slice(second_stack.as_const()).unwrap()[3].is_null());
    }

    #[test]
    fn source_port__ichisort__comp_at_rank__line_467() {
        assert_eq!(comp_AT_RANK(0, 0), 0);
        assert_eq!(comp_AT_RANK(0, u16::MAX), -65_535);
        assert_eq!(comp_AT_RANK(u16::MAX, 0), 65_535);
        assert_eq!(comp_AT_RANK(1234, 1235), -1);
        assert_eq!(comp_AT_RANK(1235, 1234), 1);
    }

    #[test]
    fn source_port__ichisort__compranksinvord__line_670() {
        assert_eq!(CompRanksInvOrd(0, 0), 0);
        assert_eq!(CompRanksInvOrd(0, u16::MAX), 65_535);
        assert_eq!(CompRanksInvOrd(u16::MAX, 0), -65_535);
        assert_eq!(CompRanksInvOrd(1234, 1235), 1);
        assert_eq!(CompRanksInvOrd(1235, 1234), -1);
    }

    #[test]
    fn source_port__ichisort__compneighborsrankscounteql__line_677() {
        let mut heap = SourceHeap::default();
        let ranks = heap
            .allocate_model_storage(vec![u16::MAX, 0, 77, 77])
            .unwrap();
        let mut globals = CANON_GLOBALS {
            m_pn_RankForSort: ranks.as_const(),
            m_nNumCompNeighborsRanksCountEql: i32::MAX,
            ..CANON_GLOBALS::default()
        };
        assert_eq!(
            CompNeighborsRanksCountEql(0, 1, heap.slice(ranks.as_const()).unwrap(), &mut globals),
            Ok(65_535)
        );
        assert_eq!(globals.m_nNumCompNeighborsRanksCountEql, i32::MAX);
        assert_eq!(
            CompNeighborsRanksCountEql(1, 0, heap.slice(ranks.as_const()).unwrap(), &mut globals),
            Ok(-65_535)
        );
        assert_eq!(globals.m_nNumCompNeighborsRanksCountEql, i32::MAX);
        assert_eq!(
            CompNeighborsRanksCountEql(2, 3, heap.slice(ranks.as_const()).unwrap(), &mut globals),
            Ok(0)
        );
        assert_eq!(globals.m_nNumCompNeighborsRanksCountEql, i32::MIN);
        assert_eq!(
            CompNeighborsRanksCountEql(3, 3, heap.slice(ranks.as_const()).unwrap(), &mut globals),
            Ok(0)
        );
        assert_eq!(globals.m_nNumCompNeighborsRanksCountEql, i32::MIN + 1);

        assert_eq!(
            CompNeighborsRanksCountEql(4, 0, heap.slice(ranks.as_const()).unwrap(), &mut globals),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(globals.m_nNumCompNeighborsRanksCountEql, i32::MIN + 1);
        assert_eq!(
            CompNeighborsRanksCountEql(0, 0, &[], &mut globals),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(globals.m_nNumCompNeighborsRanksCountEql, i32::MIN + 1);
    }

    #[test]
    fn source_port__ichisort__compneighborsat_number__line_453() {
        fn compare(
            heap: &SourceHeap,
            left: AT_NUMB,
            right: AT_NUMB,
            globals: &CANON_GLOBALS,
        ) -> Result<i32, SourceHeapError> {
            CompNeighborsAT_NUMBER(
                left,
                right,
                CompNeighborsATNumberContext::Globals { heap, globals },
            )
        }

        let mut heap = SourceHeap::default();
        let neighbors = heap.allocate_model_storage(vec![2_u16, 0, 3, 1]).unwrap();
        let ranks = heap
            .allocate_model_storage(vec![50_u16, 10, 30, 30])
            .unwrap();
        let globals = CANON_GLOBALS {
            m_pNeighborsForSort: neighbors.as_const(),
            m_pn_RankForSort: ranks.as_const(),
            ..CANON_GLOBALS::default()
        };

        assert_eq!(compare(&heap, 0, 1, &globals), Ok(-20));
        assert_eq!(compare(&heap, 1, 0, &globals), Ok(20));
        assert_eq!(compare(&heap, 0, 2, &globals), Ok(0));
        assert_eq!(compare(&heap, 3, 3, &globals), Ok(0));
        assert_eq!(
            compare(&heap, 4, 0, &globals),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let invalid_neighbor = CANON_GLOBALS {
            m_pNeighborsForSort: heap
                .allocate_model_storage(vec![99_u16])
                .unwrap()
                .as_const(),
            m_pn_RankForSort: ranks.as_const(),
            ..CANON_GLOBALS::default()
        };
        assert_eq!(
            compare(&heap, 0, 0, &invalid_neighbor),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            compare(&heap, 0, 0, &CANON_GLOBALS::default()),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(
            CompNeighborsAT_NUMBER(
                0,
                1,
                CompNeighborsATNumberContext::Slices {
                    neighbors: &[2, 0, 3, 1],
                    ranks: &[50, 10, 30, 30],
                },
            ),
            Ok(-20)
        );
    }

    #[test]
    fn source_port__ichisort__compatominvariants2only__line_502() {
        fn invariant(values: [AT_NUMB; 10], isotope: i64, auxiliary: i8) -> ATOM_INVARIANT2 {
            ATOM_INVARIANT2 {
                val: values,
                iso_sort_key: isotope,
                iso_aux_key: auxiliary,
            }
        }

        let mut heap = SourceHeap::default();
        let pointer = heap
            .allocate_model_storage(vec![
                invariant([1, 2, 3, 4, 5, 6, 7, 10, 20, 30], 5, -1),
                invariant([1, 2, 3, 4, 5, 6, 9, 40, 50, 60], 5, 1),
                invariant([1, 2, 3, 4, 5, 6, 7, 40, 20, 60], 5, 1),
                invariant([1, 2, 3, 4, 5, 6, 7, 40, 50, 60], 4, 1),
                invariant([1, 2, 3, 4, 5, 6, 7, 40, 50, 60], 5, 1),
                invariant([1, 2, 3, 4, 5, 6, 7, 41, 51, 61], 5, -1),
            ])
            .unwrap();
        let globals = CANON_GLOBALS {
            m_pAtomInvariant2ForSort: pointer.as_const(),
            ..CANON_GLOBALS::default()
        };

        assert_eq!(CompAtomInvariants2Only(&heap, 0, 1, &globals), Ok(-2));
        assert_eq!(CompAtomInvariants2Only(&heap, 0, 3, &globals), Ok(1));
        assert_eq!(CompAtomInvariants2Only(&heap, 3, 0, &globals), Ok(-1));
        // The source's second loop returns zero at its first equal slot.
        assert_eq!(CompAtomInvariants2Only(&heap, 0, 2, &globals), Ok(0));
        assert_eq!(CompAtomInvariants2Only(&heap, 4, 5, &globals), Ok(1));
        assert_eq!(CompAtomInvariants2Only(&heap, 5, 4, &globals), Ok(-1));
        assert_eq!(CompAtomInvariants2Only(&heap, 4, 4, &globals), Ok(0));
        assert_eq!(
            CompAtomInvariants2Only(&heap, 4, 6, &globals),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            CompAtomInvariants2Only(&heap, 0, 0, &CANON_GLOBALS::default()),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichisort__compatominvariants2__line_536() {
        let mut heap = SourceHeap::default();
        let mut different = ATOM_INVARIANT2::default();
        different.val[0] = 2;
        let pointer = heap
            .allocate_model_storage(vec![
                ATOM_INVARIANT2::default(),
                different,
                ATOM_INVARIANT2::default(),
            ])
            .unwrap();
        let globals = CANON_GLOBALS {
            m_pAtomInvariant2ForSort: pointer.as_const(),
            ..CANON_GLOBALS::default()
        };

        assert_eq!(CompAtomInvariants2(&heap, 0, 1, &globals), Ok(-2));
        assert_eq!(CompAtomInvariants2(&heap, 2, 0, &globals), Ok(2));
        assert_eq!(CompAtomInvariants2(&heap, 0, 2, &globals), Ok(-2));
        assert_eq!(CompAtomInvariants2(&heap, 2, 2, &globals), Ok(0));
        assert_eq!(
            CompAtomInvariants2(&heap, 3, 0, &globals),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichisort__comprank__line_475() {
        let mut heap = SourceHeap::default();
        let ranks = heap.allocate_model_storage(vec![5_u16, 2, 5, 9]).unwrap();
        let globals = CANON_GLOBALS {
            m_pn_RankForSort: ranks.as_const(),
            ..CANON_GLOBALS::default()
        };

        assert_eq!(CompRank(&heap, 0, 1, &globals), Ok(3));
        assert_eq!(CompRank(&heap, 1, 3, &globals), Ok(-7));
        assert_eq!(CompRank(&heap, 0, 2, &globals), Ok(0));
        assert_eq!(CompRank(&heap, 2, 0, &globals), Ok(0));
        assert_eq!(CompRank(&heap, 2, 2, &globals), Ok(0));
        assert_eq!(
            CompRank(&heap, 4, 0, &globals),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            CompRank(&heap, 0, 0, &CANON_GLOBALS::default()),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichisort__compranksord__line_486() {
        let mut heap = SourceHeap::default();
        let ranks = heap.allocate_model_storage(vec![5_u16, 2, 5, 9]).unwrap();
        let globals = CANON_GLOBALS {
            m_pn_RankForSort: ranks.as_const(),
            ..CANON_GLOBALS::default()
        };
        assert_eq!(CompRanksOrd(&heap, 0, 1, &globals), Ok(3));
        assert_eq!(CompRanksOrd(&heap, 1, 3, &globals), Ok(-7));
        assert_eq!(CompRanksOrd(&heap, 0, 2, &globals), Ok(-2));
        assert_eq!(CompRanksOrd(&heap, 2, 0, &globals), Ok(2));
        assert_eq!(CompRanksOrd(&heap, 2, 2, &globals), Ok(0));
        assert_eq!(
            CompRanksOrd(&heap, 4, 0, &globals),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            CompRanksOrd(&heap, 0, 0, &CANON_GLOBALS::default()),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichisort__inchi_qsort__line_66() {
        let mut records = vec![3_u8, 30, 1, 10, 2, 20, 2, 21, 1, 11];
        let mut calls = Vec::new();
        inchi_qsort(&mut records, 5, 2, &mut |left, right| {
            calls.push((left[0], right[0]));
            Ok(i32::from(left[0]) - i32::from(right[0]))
        })
        .unwrap();
        assert_eq!(records, vec![1, 11, 1, 10, 2, 20, 2, 21, 3, 30]);
        assert_eq!(calls.first(), Some(&(3, 2)));

        let mut equal = vec![7_u8; 64];
        inchi_qsort(&mut equal, 32, 2, &mut |left, right| {
            Ok(i32::from(left[0]) - i32::from(right[0]))
        })
        .unwrap();
        assert_eq!(equal, vec![7_u8; 64]);

        let mut no_op = Vec::new();
        assert_eq!(inchi_qsort(&mut no_op, 0, 0, &mut |_, _| Ok(0)), Ok(()));
        assert_eq!(
            inchi_qsort(&mut [1_u8, 2], 2, 0, &mut |_, _| Ok(0)),
            Err(SourceHeapError::UnsupportedSourceBehavior)
        );
        assert_eq!(
            inchi_qsort(&mut [1_u8], 2, 1, &mut |_, _| Ok(0)),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut partial = vec![3_u8, 1, 2];
        assert_eq!(
            inchi_qsort(&mut partial, 3, 1, &mut |_, _| {
                Err(SourceHeapError::UnsupportedSourceBehavior)
            }),
            Err(SourceHeapError::UnsupportedSourceBehavior)
        );
        assert_eq!(partial, vec![3, 1, 2]);
    }

    #[test]
    fn source_port__ichisort__compchemelemlex__line_551() {
        assert_eq!(CompChemElemLex(b"C ", b"C "), Ok(0));
        assert_eq!(CompChemElemLex(b"Br", b"Cl"), Ok(-1));
        assert_eq!(CompChemElemLex(b"Cl", b"Br"), Ok(1));
        assert_eq!(CompChemElemLex(b"Ca", b"Cl"), Ok(-11));
        assert_eq!(CompChemElemLex(&[0xff, 0], &[0x01, 0]), Ok(254));
        assert_eq!(
            CompChemElemLex(b"C", b"C "),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            CompChemElemLex(b"C ", b"C"),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichisort__compare_neighlists__line_607() {
        let mut heap = SourceHeap::default();
        let ranks = heap.allocate_model_storage(vec![10_u16, 20, 30]).unwrap();
        let first = heap.allocate_model_storage(vec![2_u16, 0, 1]).unwrap();
        let second = heap.allocate_model_storage(vec![2_u16, 0, 2]).unwrap();
        let equal = heap.allocate_model_storage(vec![2_u16, 0, 1]).unwrap();
        let short = heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
        let globals = CANON_GLOBALS {
            m_pn_RankForSort: ranks.as_const(),
            ..CANON_GLOBALS::default()
        };

        assert_eq!(compare_NeighLists(&heap, first, second, &globals), Ok(-10));
        assert_eq!(compare_NeighLists(&heap, second, first, &globals), Ok(10));
        assert_eq!(compare_NeighLists(&heap, first, equal, &globals), Ok(0));
        assert_eq!(compare_NeighLists(&heap, short, first, &globals), Ok(-1));
        assert_eq!(compare_NeighLists(&heap, first, short, &globals), Ok(1));
        assert_eq!(
            compare_NeighLists(&heap, SourceMutPointer::null(), first, &globals),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(
            compare_NeighLists(&heap, first, second, &CANON_GLOBALS::default()),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichisort__compneighlistranks__line_615() {
        let mut heap = SourceHeap::default();
        let first = heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
        let second = heap.allocate_model_storage(vec![1_u16, 1]).unwrap();
        let equal = heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
        let neighbor_lists = heap
            .allocate_model_storage(vec![first, second, equal, second])
            .unwrap();
        let ranks = heap.allocate_model_storage(vec![7_u16, 3, 7, 7]).unwrap();
        let globals = CANON_GLOBALS {
            m_pNeighList_RankForSort: neighbor_lists.as_const(),
            m_pn_RankForSort: ranks.as_const(),
            ..CANON_GLOBALS::default()
        };

        assert_eq!(CompNeighListRanks(&heap, 0, 1, &globals), Ok(4));
        assert_eq!(CompNeighListRanks(&heap, 1, 0, &globals), Ok(-4));
        assert_eq!(CompNeighListRanks(&heap, 0, 2, &globals), Ok(0));
        assert_eq!(CompNeighListRanks(&heap, 2, 0, &globals), Ok(0));
        assert_eq!(CompNeighListRanks(&heap, 0, 3, &globals), Ok(4));
        assert_eq!(CompNeighListRanks(&heap, 3, 0, &globals), Ok(-4));
        assert_eq!(CompNeighListRanks(&heap, 0, 0, &globals), Ok(0));
        assert_eq!(
            CompNeighListRanks(&heap, 4, 0, &globals),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let missing_lists = CANON_GLOBALS {
            m_pn_RankForSort: ranks.as_const(),
            ..CANON_GLOBALS::default()
        };
        assert_eq!(
            CompNeighListRanks(&heap, 0, 2, &missing_lists),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(
            CompNeighListRanks(&heap, 0, 0, &CANON_GLOBALS::default()),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichisort__compneighlists__line_632() {
        let mut heap = SourceHeap::default();
        let ranks = heap.allocate_model_storage(vec![10_u16, 20, 30]).unwrap();
        let first = heap.allocate_model_storage(vec![2_u16, 0, 1]).unwrap();
        let second = heap.allocate_model_storage(vec![2_u16, 0, 2]).unwrap();
        let equal = heap.allocate_model_storage(vec![2_u16, 0, 1]).unwrap();
        let short = heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
        let neighbor_lists = heap
            .allocate_model_storage(vec![first, second, equal, short])
            .unwrap();
        let globals = CANON_GLOBALS {
            m_pNeighList_RankForSort: neighbor_lists.as_const(),
            m_pn_RankForSort: ranks.as_const(),
            ..CANON_GLOBALS::default()
        };

        assert_eq!(CompNeighLists(&heap, 0, 1, &globals), Ok(-10));
        assert_eq!(CompNeighLists(&heap, 1, 0, &globals), Ok(10));
        assert_eq!(CompNeighLists(&heap, 0, 2, &globals), Ok(0));
        assert_eq!(CompNeighLists(&heap, 3, 0, &globals), Ok(-1));
        assert_eq!(CompNeighLists(&heap, 0, 3, &globals), Ok(1));
        assert_eq!(
            CompNeighLists(&heap, 4, 0, &globals),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            CompNeighLists(&heap, 0, 0, &CANON_GLOBALS::default()),
            Err(SourceHeapError::NullPointer)
        );

        let missing_ranks = CANON_GLOBALS {
            m_pNeighList_RankForSort: neighbor_lists.as_const(),
            ..CANON_GLOBALS::default()
        };
        assert_eq!(
            CompNeighLists(&heap, 0, 1, &missing_ranks),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn official_c_oracle__compneighlists__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/official_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--comp-neigh-lists-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "CompNeighLists");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_array = |field: &str| {
                official["input"][field]
                    .as_array()
                    .unwrap_or_else(|| panic!("{case_id}: {field} must be an array"))
                    .iter()
                    .map(|value| {
                        u16::try_from(
                            value
                                .as_u64()
                                .unwrap_or_else(|| panic!("{case_id}: {field} must contain u16")),
                        )
                        .unwrap_or_else(|_| panic!("{case_id}: {field} value exceeds u16"))
                    })
                    .collect::<Vec<_>>()
            };
            let first_values = parse_array("first_list");
            let second_values = parse_array("second_list");
            let ranks_values = parse_array("ranks");
            let first_len = usize::from(first_values[0]);
            let second_len = usize::from(second_values[0]);
            let first_index = u16::try_from(
                official["input"]["first_index"]
                    .as_u64()
                    .expect("first_index must be u16"),
            )
            .expect("first_index exceeds u16");
            let second_index = u16::try_from(
                official["input"]["second_index"]
                    .as_u64()
                    .expect("second_index must be u16"),
            )
            .expect("second_index exceeds u16");

            let mut heap = SourceHeap::default();
            let ranks = heap.allocate_model_storage(ranks_values).unwrap();
            let first = heap
                .allocate_model_storage(first_values[..=first_len].to_vec())
                .unwrap();
            let second = heap
                .allocate_model_storage(second_values[..=second_len].to_vec())
                .unwrap();
            let neighbor_lists = heap.allocate_model_storage(vec![first, second]).unwrap();
            let globals = CANON_GLOBALS {
                m_pNeighList_RankForSort: neighbor_lists.as_const(),
                m_pn_RankForSort: ranks.as_const(),
                ..CANON_GLOBALS::default()
            };
            let first_before = heap.slice(first.as_const()).unwrap().to_vec();
            let second_before = heap.slice(second.as_const()).unwrap().to_vec();
            let ranks_before = heap.slice(ranks.as_const()).unwrap().to_vec();
            let lists_before = heap.slice(neighbor_lists.as_const()).unwrap().to_vec();
            let globals_before = globals.clone();
            let rust = CompNeighLists(&heap, first_index, second_index, &globals).unwrap();
            let expected = i32::try_from(
                official["output"]["result"]
                    .as_i64()
                    .expect("result must be an integer"),
            )
            .expect("official result must fit i32");
            assert_eq!(rust, expected, "{case_id}");
            assert_eq!(heap.slice(first.as_const()).unwrap(), first_before);
            assert_eq!(heap.slice(second.as_const()).unwrap(), second_before);
            assert_eq!(heap.slice(ranks.as_const()).unwrap(), ranks_before);
            assert_eq!(heap.slice(neighbor_lists.as_const()).unwrap(), lists_before);
            assert_eq!(globals, globals_before);
            for field in [
                "first_list_unchanged",
                "second_list_unchanged",
                "ranks_unchanged",
                "list_pointers_unchanged",
                "indices_unchanged",
                "globals_unchanged",
            ] {
                assert_eq!(official["output"][field], true, "{case_id}: {field}");
            }
            record_count += 1;
        }
        assert_eq!(record_count, 7);
    }

    #[test]
    fn source_port__ichisort__compneighlistsuptomaxrank__line_644() {
        let mut heap = SourceHeap::default();
        let ranks = heap
            .allocate_model_storage(vec![1_u16, 2, 3, 9, 10, u16::MAX])
            .unwrap();
        let first = heap.allocate_model_storage(vec![3_u16, 0, 1, 3]).unwrap();
        let second = heap.allocate_model_storage(vec![3_u16, 0, 2, 4]).unwrap();
        let short = heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
        let long = heap.allocate_model_storage(vec![2_u16, 0, 1]).unwrap();
        let boundary_high = heap.allocate_model_storage(vec![1_u16, 5]).unwrap();
        let boundary_low = heap.allocate_model_storage(vec![1_u16, 4]).unwrap();
        let neighbor_lists = heap
            .allocate_model_storage(vec![
                first,
                second,
                short,
                long,
                boundary_high,
                boundary_low,
            ])
            .unwrap();
        let mut globals = CANON_GLOBALS {
            m_pNeighList_RankForSort: neighbor_lists.as_const(),
            m_pn_RankForSort: ranks.as_const(),
            m_nMaxAtNeighRankForSort: 10,
            ..CANON_GLOBALS::default()
        };

        assert_eq!(CompNeighListsUpToMaxRank(&heap, 0, 1, &globals), Ok(-1));
        assert_eq!(CompNeighListsUpToMaxRank(&heap, 1, 0, &globals), Ok(1));

        globals.m_nMaxAtNeighRankForSort = 1;
        assert_eq!(CompNeighListsUpToMaxRank(&heap, 0, 1, &globals), Ok(0));

        globals.m_nMaxAtNeighRankForSort = 2;
        assert_eq!(CompNeighListsUpToMaxRank(&heap, 2, 3, &globals), Ok(-1));
        assert_eq!(CompNeighListsUpToMaxRank(&heap, 3, 2, &globals), Ok(1));

        globals.m_nMaxAtNeighRankForSort = u16::MAX;
        assert_eq!(
            CompNeighListsUpToMaxRank(&heap, 4, 5, &globals),
            Ok(i32::from(u16::MAX) - 10)
        );
        assert_eq!(
            CompNeighListsUpToMaxRank(&heap, 6, 0, &globals),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            CompNeighListsUpToMaxRank(&heap, 0, 0, &CANON_GLOBALS::default()),
            Err(SourceHeapError::NullPointer)
        );

        let missing_ranks = CANON_GLOBALS {
            m_pNeighList_RankForSort: neighbor_lists.as_const(),
            m_nMaxAtNeighRankForSort: u16::MAX,
            ..CANON_GLOBALS::default()
        };
        assert_eq!(
            CompNeighListsUpToMaxRank(&heap, 0, 1, &missing_ranks),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn official_c_oracle__compneighlistsuptomaxrank__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/official_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--comp-neigh-lists-up-to-max-rank-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "CompNeighListsUpToMaxRank");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_array = |field: &str| {
                official["input"][field]
                    .as_array()
                    .unwrap_or_else(|| panic!("{case_id}: {field} must be an array"))
                    .iter()
                    .map(|value| {
                        u16::try_from(
                            value
                                .as_u64()
                                .unwrap_or_else(|| panic!("{case_id}: {field} must contain u16")),
                        )
                        .unwrap_or_else(|_| panic!("{case_id}: {field} value exceeds u16"))
                    })
                    .collect::<Vec<_>>()
            };
            let first_values = parse_array("first_list");
            let second_values = parse_array("second_list");
            let ranks_values = parse_array("ranks");
            let first_len = usize::from(first_values[0]);
            let second_len = usize::from(second_values[0]);
            let first_index = u16::try_from(
                official["input"]["first_index"]
                    .as_u64()
                    .expect("first_index must be u16"),
            )
            .expect("first_index exceeds u16");
            let second_index = u16::try_from(
                official["input"]["second_index"]
                    .as_u64()
                    .expect("second_index must be u16"),
            )
            .expect("second_index exceeds u16");
            let max_rank = u16::try_from(
                official["input"]["max_rank"]
                    .as_u64()
                    .expect("max_rank must be u16"),
            )
            .expect("max_rank exceeds u16");

            let mut heap = SourceHeap::default();
            let ranks = heap.allocate_model_storage(ranks_values).unwrap();
            let first = heap
                .allocate_model_storage(first_values[..=first_len].to_vec())
                .unwrap();
            let second = heap
                .allocate_model_storage(second_values[..=second_len].to_vec())
                .unwrap();
            let neighbor_lists = heap.allocate_model_storage(vec![first, second]).unwrap();
            let globals = CANON_GLOBALS {
                m_pNeighList_RankForSort: neighbor_lists.as_const(),
                m_pn_RankForSort: ranks.as_const(),
                m_nMaxAtNeighRankForSort: max_rank,
                ..CANON_GLOBALS::default()
            };
            let first_before = heap.slice(first.as_const()).unwrap().to_vec();
            let second_before = heap.slice(second.as_const()).unwrap().to_vec();
            let ranks_before = heap.slice(ranks.as_const()).unwrap().to_vec();
            let lists_before = heap.slice(neighbor_lists.as_const()).unwrap().to_vec();
            let globals_before = globals.clone();
            let rust =
                CompNeighListsUpToMaxRank(&heap, first_index, second_index, &globals).unwrap();
            let expected = i32::try_from(
                official["output"]["result"]
                    .as_i64()
                    .expect("result must be an integer"),
            )
            .expect("official result must fit i32");
            assert_eq!(rust, expected, "{case_id}");
            assert_eq!(heap.slice(first.as_const()).unwrap(), first_before);
            assert_eq!(heap.slice(second.as_const()).unwrap(), second_before);
            assert_eq!(heap.slice(ranks.as_const()).unwrap(), ranks_before);
            assert_eq!(heap.slice(neighbor_lists.as_const()).unwrap(), lists_before);
            assert_eq!(globals, globals_before);
            for field in [
                "first_list_unchanged",
                "second_list_unchanged",
                "ranks_unchanged",
                "list_pointers_unchanged",
                "indices_unchanged",
                "globals_unchanged",
            ] {
                assert_eq!(official["output"][field], true, "{case_id}: {field}");
            }
            record_count += 1;
        }
        assert_eq!(record_count, 8);
    }

    #[test]
    fn source_port__ichisort__compneighlistranksord__line_657() {
        let mut heap = SourceHeap::default();
        let first = heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
        let second = heap.allocate_model_storage(vec![1_u16, 1]).unwrap();
        let equal = heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
        let neighbor_lists = heap
            .allocate_model_storage(vec![first, second, equal, second])
            .unwrap();
        let ranks = heap.allocate_model_storage(vec![7_u16, 3, 7, 7]).unwrap();
        let globals = CANON_GLOBALS {
            m_pNeighList_RankForSort: neighbor_lists.as_const(),
            m_pn_RankForSort: ranks.as_const(),
            ..CANON_GLOBALS::default()
        };

        assert_eq!(CompNeighListRanksOrd(&heap, 0, 1, &globals), Ok(4));
        assert_eq!(CompNeighListRanksOrd(&heap, 1, 0, &globals), Ok(-4));
        assert_eq!(CompNeighListRanksOrd(&heap, 0, 2, &globals), Ok(-2));
        assert_eq!(CompNeighListRanksOrd(&heap, 2, 0, &globals), Ok(2));
        assert_eq!(CompNeighListRanksOrd(&heap, 0, 3, &globals), Ok(4));
        assert_eq!(CompNeighListRanksOrd(&heap, 3, 0, &globals), Ok(-4));
        assert_eq!(CompNeighListRanksOrd(&heap, 0, 0, &globals), Ok(0));
        assert_eq!(
            CompNeighListRanksOrd(&heap, 4, 0, &globals),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let rank_short_circuit = CANON_GLOBALS {
            m_pn_RankForSort: ranks.as_const(),
            ..CANON_GLOBALS::default()
        };
        assert_eq!(
            CompNeighListRanksOrd(&heap, 0, 1, &rank_short_circuit),
            Ok(4)
        );
        assert_eq!(
            CompNeighListRanksOrd(&heap, 0, 2, &rank_short_circuit),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichisort__insertions_sort_at_numbers__line_331() {
        let mut heap = SourceHeap::default();
        let allocation = heap
            .allocate_model_storage(vec![99_u16, 21, 11, 22, 12, 88])
            .unwrap();
        let base = allocation.offset(1).unwrap();
        let mut calls = Vec::new();
        assert_eq!(
            insertions_sort_AT_NUMBERS(&mut heap, base, 4, &mut |_, left, right| {
                calls.push((left, right));
                Ok(i32::from(left / 10) - i32::from(right / 10))
            }),
            Ok(3)
        );
        assert_eq!(
            heap.slice(allocation.as_const()).unwrap(),
            &[99, 11, 12, 21, 22, 88]
        );
        assert_eq!(
            calls,
            vec![(21, 11), (21, 22), (22, 12), (21, 12), (11, 12)]
        );

        let null = SourceMutPointer::<AT_NUMB>::null();
        let mut call_count = 0;
        assert_eq!(
            insertions_sort_AT_NUMBERS(&mut heap, null, 0, &mut |_, _, _| {
                call_count += 1;
                Ok(0)
            }),
            Ok(0)
        );
        assert_eq!(
            insertions_sort_AT_NUMBERS(&mut heap, null, -1, &mut |_, _, _| {
                call_count += 1;
                Ok(0)
            }),
            Ok(0)
        );
        assert_eq!(call_count, 0);

        let partial = heap.allocate_model_storage(vec![3_u16, 2]).unwrap();
        assert_eq!(
            insertions_sort_AT_NUMBERS(&mut heap, partial, 3, &mut |_, left, right| {
                Ok(i32::from(left) - i32::from(right))
            }),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(partial.as_const()).unwrap(), &[2, 3]);

        let descending = heap.allocate_model_storage(vec![1_u16, 2, 3, 4]).unwrap();
        assert_eq!(
            insertions_sort_AT_NUMBERS(&mut heap, descending, 4, &mut |_, left, right| {
                Ok(i32::from(right) - i32::from(left))
            }),
            Ok(6)
        );
        assert_eq!(heap.slice(descending.as_const()).unwrap(), &[4, 3, 2, 1]);
    }

    #[test]
    fn source_port__ichisort__insertions_sort_neighlist_at_numbers__line_355() {
        let mut heap = SourceHeap::default();
        let ranks = heap
            .allocate_model_storage(vec![50_u16, 10, 40, 20, 30, 20])
            .unwrap();
        let allocation = heap
            .allocate_model_storage(vec![99_u16, 5, 0, 1, 2, 3, 4, 88])
            .unwrap();
        let list = allocation.offset(1).unwrap();
        assert_eq!(
            insertions_sort_NeighList_AT_NUMBERS(&mut heap, list, ranks),
            Ok(())
        );
        assert_eq!(
            heap.slice(allocation.as_const()).unwrap(),
            &[99, 5, 1, 3, 4, 2, 0, 88]
        );

        let stable = heap
            .allocate_model_storage(vec![4_u16, 5, 3, 1, 0])
            .unwrap();
        assert_eq!(
            insertions_sort_NeighList_AT_NUMBERS(&mut heap, stable, ranks),
            Ok(())
        );
        assert_eq!(heap.slice(stable.as_const()).unwrap(), &[4, 1, 5, 3, 0]);

        let empty = heap.allocate_model_storage(vec![0_u16]).unwrap();
        let singleton = heap.allocate_model_storage(vec![1_u16, 4]).unwrap();
        assert_eq!(
            insertions_sort_NeighList_AT_NUMBERS(&mut heap, empty, SourceMutPointer::null()),
            Ok(())
        );
        assert_eq!(
            insertions_sort_NeighList_AT_NUMBERS(&mut heap, singleton, SourceMutPointer::null()),
            Ok(())
        );

        let malformed = heap.allocate_model_storage(vec![3_u16, 0, 1]).unwrap();
        assert_eq!(
            insertions_sort_NeighList_AT_NUMBERS(&mut heap, malformed, ranks),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(malformed.as_const()).unwrap(), &[3, 1, 0]);

        let invalid_atom = heap.allocate_model_storage(vec![2_u16, 0, 99]).unwrap();
        assert_eq!(
            insertions_sort_NeighList_AT_NUMBERS(&mut heap, invalid_atom, ranks),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(invalid_atom.as_const()).unwrap(), &[2, 0, 99]);
    }

    #[test]
    fn source_port__ichisort__compareneighlistlex__line_560() {
        let mut heap = SourceHeap::default();
        let ranks = heap
            .allocate_model_storage(vec![5_u16, 2, 9, 2, u16::MAX, 0])
            .unwrap();
        let first = heap.allocate_model_storage(vec![3_u16, 1, 0, 4]).unwrap();
        let second = heap.allocate_model_storage(vec![3_u16, 3, 0, 4]).unwrap();
        assert_eq!(CompareNeighListLex(&heap, first, second, ranks), Ok(0));

        heap.slice_mut(second).unwrap()[2] = 2;
        assert_eq!(CompareNeighListLex(&heap, first, second, ranks), Ok(-4));
        assert_eq!(CompareNeighListLex(&heap, second, first, ranks), Ok(4));

        let prefix = heap.allocate_model_storage(vec![2_u16, 1, 0]).unwrap();
        let longer = heap.allocate_model_storage(vec![3_u16, 3, 0, 1]).unwrap();
        assert_eq!(CompareNeighListLex(&heap, prefix, longer, ranks), Ok(-1));
        assert_eq!(CompareNeighListLex(&heap, longer, prefix, ranks), Ok(1));

        let empty1 = heap.allocate_model_storage(vec![0_u16]).unwrap();
        let empty2 = heap.allocate_model_storage(vec![0_u16, 99]).unwrap();
        assert_eq!(CompareNeighListLex(&heap, empty1, empty2, ranks), Ok(0));

        let max_rank = heap
            .allocate_model_storage(vec![2_u16, 4, u16::MAX])
            .unwrap();
        let min_rank = heap
            .allocate_model_storage(vec![2_u16, 5, u16::MAX])
            .unwrap();
        assert_eq!(
            CompareNeighListLex(&heap, max_rank, min_rank, ranks),
            Ok(i32::from(u16::MAX))
        );

        let short_circuit1 = heap.allocate_model_storage(vec![2_u16, 2]).unwrap();
        let short_circuit2 = heap.allocate_model_storage(vec![2_u16, 0]).unwrap();
        assert_eq!(
            CompareNeighListLex(&heap, short_circuit1, short_circuit2, ranks),
            Ok(4)
        );

        let malformed1 = heap.allocate_model_storage(vec![2_u16, 1]).unwrap();
        let malformed2 = heap.allocate_model_storage(vec![2_u16, 3]).unwrap();
        assert_eq!(
            CompareNeighListLex(&heap, malformed1, malformed2, ranks),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichisort__compareneighlistlexuptomaxrank__line_582() {
        let mut heap = SourceHeap::default();
        let ranks = heap
            .allocate_model_storage(vec![1_u16, 3, 5, 7, u16::MAX])
            .unwrap();
        let first = heap
            .allocate_model_storage(vec![4_u16, 0, 1, 2, 4])
            .unwrap();
        let second = heap
            .allocate_model_storage(vec![4_u16, 0, 1, 3, 4])
            .unwrap();

        assert_eq!(
            CompareNeighListLexUpToMaxRank(&heap, first, second, ranks, 5),
            Ok(1)
        );
        assert_eq!(
            CompareNeighListLexUpToMaxRank(&heap, first, second, ranks, 7),
            Ok(-2)
        );
        assert_eq!(
            CompareNeighListLexUpToMaxRank(&heap, first, second, ranks, u16::MAX),
            Ok(-2)
        );

        let all_trimmed1 = heap.allocate_model_storage(vec![2_u16, 0, 1]).unwrap();
        let all_trimmed2 = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
        assert_eq!(
            CompareNeighListLexUpToMaxRank(&heap, all_trimmed1, all_trimmed2, ranks, 0),
            Ok(0)
        );

        let interior_high = heap.allocate_model_storage(vec![2_u16, 4, 0]).unwrap();
        let low = heap.allocate_model_storage(vec![2_u16, 0, 0]).unwrap();
        assert_eq!(
            CompareNeighListLexUpToMaxRank(&heap, interior_high, low, ranks, 5),
            Ok(i32::from(u16::MAX) - 1)
        );

        let short_circuit1 = heap.allocate_model_storage(vec![2_u16, 3, 0]).unwrap();
        let short_circuit2 = heap.allocate_model_storage(vec![2_u16, 2, 0]).unwrap();
        assert_eq!(
            CompareNeighListLexUpToMaxRank(&heap, short_circuit1, short_circuit2, ranks, u16::MAX),
            Ok(2)
        );

        let malformed = heap.allocate_model_storage(vec![2_u16, 0]).unwrap();
        assert_eq!(
            CompareNeighListLexUpToMaxRank(&heap, malformed, second, ranks, 5),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichisort__insertions_sort_neighlist_at_numbers3__line_396() {
        let mut heap = SourceHeap::default();
        let ranks = heap
            .allocate_model_storage(vec![3_u16, 1, 2, 1, 3])
            .unwrap();
        let list = heap
            .allocate_model_storage(vec![5_u16, 0, 1, 2, 3, 4, 99])
            .unwrap();
        assert_eq!(
            insertions_sort_NeighList_AT_NUMBERS3(&mut heap, list, ranks),
            Ok(4)
        );
        assert_eq!(
            heap.slice(list.as_const()).unwrap(),
            &[5, 1, 3, 2, 0, 4, 99]
        );

        let descending_ranks = heap
            .allocate_model_storage(vec![5_u16, 4, 3, 2, 1])
            .unwrap();
        let descending = heap
            .allocate_model_storage(vec![5_u16, 0, 1, 2, 3, 4])
            .unwrap();
        assert_eq!(
            insertions_sort_NeighList_AT_NUMBERS3(&mut heap, descending, descending_ranks),
            Ok(10)
        );
        assert_eq!(
            heap.slice(descending.as_const()).unwrap(),
            &[5, 4, 3, 2, 1, 0]
        );

        let empty = heap.allocate_model_storage(vec![0_u16]).unwrap();
        let one = heap.allocate_model_storage(vec![1_u16, u16::MAX]).unwrap();
        assert_eq!(
            insertions_sort_NeighList_AT_NUMBERS3(&mut heap, empty, SourceMutPointer::null()),
            Ok(0)
        );
        assert_eq!(
            insertions_sort_NeighList_AT_NUMBERS3(&mut heap, one, SourceMutPointer::null()),
            Ok(0)
        );

        let partial_ranks = heap.allocate_model_storage(vec![2_u16, 1]).unwrap();
        let partial = heap.allocate_model_storage(vec![3_u16, 0, 1]).unwrap();
        assert_eq!(
            insertions_sort_NeighList_AT_NUMBERS3(&mut heap, partial, partial_ranks),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(partial.as_const()).unwrap(), &[3, 1, 0]);
    }

    #[test]
    fn source_port__ichisort__insertions_sort_neighlistbysymmandcanonrank__line_421() {
        let mut heap = SourceHeap::default();
        let symmetry = heap
            .allocate_model_storage(vec![2_u16, 1, 2, 3, 2])
            .unwrap();
        let canonical = heap
            .allocate_model_storage(vec![1_u16, 9, 5, 0, 5])
            .unwrap();
        let list = heap
            .allocate_model_storage(vec![5_u16, 0, 1, 2, 3, 4, 99])
            .unwrap();
        assert_eq!(
            insertions_sort_NeighListBySymmAndCanonRank(&mut heap, list, symmetry, canonical),
            Ok(())
        );
        assert_eq!(
            heap.slice(list.as_const()).unwrap(),
            &[5, 3, 2, 4, 0, 1, 99]
        );

        let extrema_symmetry = heap
            .allocate_model_storage(vec![0_u16, u16::MAX, u16::MAX])
            .unwrap();
        let extrema_canonical = heap
            .allocate_model_storage(vec![u16::MAX, 0, u16::MAX])
            .unwrap();
        let extrema = heap.allocate_model_storage(vec![3_u16, 0, 1, 2]).unwrap();
        assert_eq!(
            insertions_sort_NeighListBySymmAndCanonRank(
                &mut heap,
                extrema,
                extrema_symmetry,
                extrema_canonical,
            ),
            Ok(())
        );
        assert_eq!(heap.slice(extrema.as_const()).unwrap(), &[3, 2, 1, 0]);

        let empty = heap.allocate_model_storage(vec![0_u16]).unwrap();
        let one = heap.allocate_model_storage(vec![1_u16, u16::MAX]).unwrap();
        assert_eq!(
            insertions_sort_NeighListBySymmAndCanonRank(
                &mut heap,
                empty,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(())
        );
        assert_eq!(
            insertions_sort_NeighListBySymmAndCanonRank(
                &mut heap,
                one,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(())
        );

        let partial_symmetry = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
        let partial_canonical = heap.allocate_model_storage(vec![0_u16, 0]).unwrap();
        let partial = heap.allocate_model_storage(vec![4_u16, 0, 1]).unwrap();
        assert_eq!(
            insertions_sort_NeighListBySymmAndCanonRank(
                &mut heap,
                partial,
                partial_symmetry,
                partial_canonical,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(partial.as_const()).unwrap(), &[4, 1, 0]);

        let invalid_atom = heap.allocate_model_storage(vec![2_u16, 0, 9]).unwrap();
        assert_eq!(
            insertions_sort_NeighListBySymmAndCanonRank(
                &mut heap,
                invalid_atom,
                symmetry,
                canonical,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichisort__insertions_sort_at_rank__line_375() {
        let mut negative = [3, 2, 1];
        assert_eq!(insertions_sort_AT_RANK(&mut negative, -1), Ok(0));
        assert_eq!(negative, [3, 2, 1]);

        let mut ranks: [AT_RANK; 8] = [5, 1, 3, 3, 0, u16::MAX, u16::MIN, 2];
        let original = ranks;
        let expected_transitions = original
            .iter()
            .enumerate()
            .map(|(left, value)| {
                original[left + 1..]
                    .iter()
                    .filter(|right| value > right)
                    .count()
            })
            .sum::<usize>() as i32;
        assert_eq!(
            insertions_sort_AT_RANK(&mut ranks, original.len() as i32),
            Ok(expected_transitions)
        );
        assert_eq!(ranks, [u16::MIN, 0, 1, 2, 3, 3, 5, u16::MAX]);

        let mut prefix = [4, 2, 3, 1, 99, 98];
        assert_eq!(insertions_sort_AT_RANK(&mut prefix, 4), Ok(5));
        assert_eq!(prefix, [1, 2, 3, 4, 99, 98]);
        assert_eq!(
            insertions_sort_AT_RANK(&mut prefix, 7),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichisort__iisort__line_1014() {
        fn permutations(values: &mut [i32], start: usize, output: &mut Vec<Vec<i32>>) {
            if start == values.len() {
                output.push(values.to_vec());
                return;
            }
            for index in start..values.len() {
                values.swap(start, index);
                permutations(values, start + 1, output);
                values.swap(start, index);
            }
        }

        let mut heap = SourceHeap::default();
        let values = heap
            .allocate_model_storage(vec![3_i32, i32::MAX, -1, 3, i32::MIN, 91])
            .unwrap();
        assert_eq!(iisort(&mut heap, values, 5), Ok(values));
        assert_eq!(
            heap.slice(values.as_const()).unwrap(),
            &[i32::MIN, -1, 3, 3, i32::MAX, 91]
        );

        let descending = heap
            .allocate_model_storage(vec![5_i32, 4, 3, 2, 1])
            .unwrap();
        assert_eq!(iisort(&mut heap, descending, 5), Ok(descending));
        assert_eq!(heap.slice(descending.as_const()).unwrap(), &[1, 2, 3, 4, 5]);

        let mut source_values = [i32::MIN, 3, -1, 3, i32::MAX];
        let mut cases = Vec::new();
        permutations(&mut source_values, 0, &mut cases);
        for case in cases {
            let pointer = heap.allocate_model_storage(case).unwrap();
            assert_eq!(iisort(&mut heap, pointer, 5), Ok(pointer));
            assert_eq!(
                heap.slice(pointer.as_const()).unwrap(),
                &[i32::MIN, -1, 3, 3, i32::MAX]
            );
        }

        let interior_allocation = heap
            .allocate_model_storage(vec![91_i32, 5, 4, 3, 2, 1, 92])
            .unwrap();
        let interior = interior_allocation.offset(1).unwrap();
        assert_eq!(iisort(&mut heap, interior, 5), Ok(interior));
        assert_eq!(
            heap.slice(interior_allocation.as_const()).unwrap(),
            &[91, 1, 2, 3, 4, 5, 92]
        );

        let null = SourceMutPointer::null();
        assert_eq!(iisort(&mut heap, null, i32::MIN), Ok(null));
        assert_eq!(iisort(&mut heap, null, 0), Ok(null));
        assert_eq!(iisort(&mut heap, null, 1), Ok(null));
        assert_eq!(
            iisort(&mut heap, null, 2),
            Err(SourceHeapError::NullPointer)
        );
        let before_error = heap.slice(descending.as_const()).unwrap().to_vec();
        assert_eq!(
            iisort(&mut heap, descending, 6),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(descending.as_const()).unwrap(), before_error);
    }

    #[test]
    fn source_port__ichisort__inchi_swap__line_286() {
        let mut bytes = *b"abcdefgh";
        assert_eq!(inchi_swap(&mut bytes, 0, 4, 4), Ok(()));
        assert_eq!(&bytes, b"efghabcd");

        let mut bytes = *b"abcdef";
        assert_eq!(inchi_swap(&mut bytes, 1, 3, 3), Ok(()));
        assert_eq!(&bytes, b"adefcb");

        let mut bytes = *b"abcd";
        assert_eq!(inchi_swap(&mut bytes, 2, 2, usize::MAX), Ok(()));
        assert_eq!(&bytes, b"abcd");
        assert_eq!(inchi_swap(&mut bytes, 0, 4, 0), Ok(()));
        assert_eq!(&bytes, b"abcd");
        assert_eq!(
            inchi_swap(&mut bytes, 0, 3, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(&bytes, b"abcd");
    }

    #[test]
    fn source_port__ichisort__insertions_sort__line_304() {
        let mut records = [2_u8, b'a', 1, b'b', 2, b'c', 1, b'd'];
        let mut comparisons = 0;
        assert_eq!(
            insertions_sort(&mut records, 4, 2, &mut |left, right| {
                comparisons += 1;
                Ok(i32::from(left[0]) - i32::from(right[0]))
            }),
            Ok(3)
        );
        assert_eq!(&records, &[1, b'b', 1, b'd', 2, b'a', 2, b'c']);
        assert_eq!(comparisons, 5);

        let mut descending = [1_u8, 2, 3, 4];
        assert_eq!(
            insertions_sort(&mut descending, 4, 1, &mut |left, right| {
                Ok(i32::from(right[0]) - i32::from(left[0]))
            }),
            Ok(6)
        );
        assert_eq!(&descending, &[4, 3, 2, 1]);

        let mut empty = [];
        assert_eq!(insertions_sort(&mut empty, 0, 4, &mut |_, _| Ok(0)), Ok(0));
        let mut one = [7_u8, 8];
        assert_eq!(insertions_sort(&mut one, 1, 2, &mut |_, _| Ok(1)), Ok(0));
        assert_eq!(
            insertions_sort(&mut one, 2, 2, &mut |_, _| Ok(0)),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(&one, &[7, 8]);
        assert_eq!(insertions_sort(&mut one, 3, 0, &mut |_, _| Ok(0)), Ok(0));
    }

    #[test]
    fn source_port__ichisort__createneighlistfromlinearct__line_701() {
        let mut heap = SourceHeap::default();
        let ct = heap
            .allocate_model_storage(vec![2_u16, 1, 3, 1, 2])
            .unwrap();
        let neighbors = CreateNeighListFromLinearCT(&mut heap, ct.as_const(), 5, 3).unwrap();
        assert!(!neighbors.is_null());
        let pointers = heap.slice(neighbors.as_const()).unwrap();
        assert_eq!(pointers.len(), 4);
        assert!(pointers[3].is_null());
        assert_eq!(pointers[0].difference(pointers[0]).unwrap(), 0);
        assert_eq!(pointers[1].difference(pointers[0]).unwrap(), 3);
        assert_eq!(pointers[2].difference(pointers[0]).unwrap(), 6);
        assert_eq!(
            &heap.slice(pointers[0].as_const()).unwrap()[..3],
            &[2, 1, 2]
        );
        assert_eq!(
            &heap.slice(pointers[1].as_const()).unwrap()[..3],
            &[2, 0, 2]
        );
        assert_eq!(
            &heap.slice(pointers[2].as_const()).unwrap()[..3],
            &[2, 0, 1]
        );

        for (values, length, atoms) in [
            (vec![4_u16], 1, 3),
            (vec![2_u16, 1], 2, 3),
            (vec![2_u16, 1, 4], 3, 3),
        ] {
            let ct = heap.allocate_model_storage(values).unwrap();
            assert!(
                CreateNeighListFromLinearCT(&mut heap, ct.as_const(), length, atoms)
                    .unwrap()
                    .is_null()
            );
        }

        for (successful, expected_calls) in [(0, 1), (1, 3), (2, 3)] {
            let mut failing_heap = SourceHeap::default();
            let ct = failing_heap
                .allocate_model_storage(vec![2_u16, 1, 3, 1, 2])
                .unwrap();
            failing_heap.fail_after_allocations(successful);
            assert!(
                CreateNeighListFromLinearCT(&mut failing_heap, ct.as_const(), 5, 3)
                    .unwrap()
                    .is_null()
            );
            assert_eq!(failing_heap.source_allocation_calls(), expected_calls);
        }
    }

    #[test]
    fn source_port__ichisort__createneighlist__line_810() {
        fn atom(neighbors: &[AT_NUMB], bond_types: &[u8], endpoint: AT_NUMB) -> sp_ATOM {
            let mut atom = sp_ATOM {
                valence: neighbors.len() as i8,
                endpoint,
                ..sp_ATOM::default()
            };
            atom.neighbor[..neighbors.len()].copy_from_slice(neighbors);
            atom.bond_type[..bond_types.len()].copy_from_slice(bond_types);
            atom
        }

        let atoms = vec![
            atom(&[1, 2], &[1, BOND_DOUBLE as u8], 0),
            atom(&[0], &[1], 0),
            atom(&[0], &[BOND_DOUBLE as u8], 0),
        ];
        let mut heap = SourceHeap::default();
        let atoms_pointer = heap.allocate_model_storage(atoms.clone()).unwrap();
        let neighbors =
            CreateNeighList(&mut heap, 3, 3, atoms_pointer.as_const(), 0, None).unwrap();
        let pointers = heap.slice(neighbors.as_const()).unwrap().to_vec();
        assert_eq!(pointers.len(), 4);
        assert!(pointers[3].is_null());
        assert_eq!(pointers[1].difference(pointers[0]).unwrap(), 3);
        assert_eq!(pointers[2].difference(pointers[0]).unwrap(), 5);
        assert_eq!(
            &heap.slice(pointers[0].as_const()).unwrap()[..7],
            &[2, 1, 2, 1, 0, 1, 0]
        );
        assert_eq!(heap.live_allocation_count(), 3);
        FreeNeighList(&mut heap, neighbors).unwrap();
        assert_eq!(heap.live_allocation_count(), 1);

        let squared = CreateNeighList(&mut heap, 3, 3, atoms_pointer.as_const(), 1, None).unwrap();
        let squared_pointers = heap.slice(squared.as_const()).unwrap().to_vec();
        assert_eq!(squared_pointers[1].difference(squared_pointers[0]), Ok(4));
        assert_eq!(squared_pointers[2].difference(squared_pointers[0]), Ok(6));
        assert!(squared_pointers[3].is_null());
        assert_eq!(
            &heap.slice(squared_pointers[0].as_const()).unwrap()[..9],
            &[3, 1, 2, 2, 1, 0, 2, 0, 0]
        );
        FreeNeighList(&mut heap, squared).unwrap();

        let mut taut_heap = SourceHeap::default();
        let taut_atoms = taut_heap
            .allocate_model_storage(vec![atom(&[], &[], 1), atom(&[], &[], 2)])
            .unwrap();
        let groups = taut_heap
            .allocate_model_storage(vec![
                T_GROUP {
                    nNumEndpoints: 2,
                    nFirstEndpointAtNoPos: 0,
                    ..T_GROUP::default()
                },
                T_GROUP {
                    nNumEndpoints: 1,
                    nFirstEndpointAtNoPos: 2,
                    ..T_GROUP::default()
                },
            ])
            .unwrap();
        let endpoint_numbers = taut_heap.allocate_model_storage(vec![0_u16, 1, 1]).unwrap();
        let info = T_GROUP_INFO {
            t_group: groups,
            nEndpointAtomNumber: endpoint_numbers,
            num_t_groups: 2,
            ..T_GROUP_INFO::default()
        };
        let taut_neighbors =
            CreateNeighList(&mut taut_heap, 2, 4, taut_atoms.as_const(), 0, Some(&info)).unwrap();
        let taut_pointers = taut_heap.slice(taut_neighbors.as_const()).unwrap().to_vec();
        assert_eq!(taut_pointers.len(), 5);
        assert!(taut_pointers[4].is_null());
        assert_eq!(taut_pointers[1].difference(taut_pointers[0]), Ok(2));
        assert_eq!(taut_pointers[2].difference(taut_pointers[0]), Ok(4));
        assert_eq!(taut_pointers[3].difference(taut_pointers[0]), Ok(7));
        assert_eq!(
            &taut_heap.slice(taut_pointers[0].as_const()).unwrap()[..9],
            &[1, 2, 1, 3, 2, 0, 1, 1, 1]
        );
        FreeNeighList(&mut taut_heap, taut_neighbors).unwrap();
        assert_eq!(taut_heap.live_allocation_count(), 3);

        let mut boundary_heap = SourceHeap::default();
        let boundary_atom = atom(&(0_u16..20).collect::<Vec<_>>(), &[1_u8; 20], 0);
        let boundary_atoms = boundary_heap
            .allocate_model_storage(vec![boundary_atom])
            .unwrap();
        let boundary =
            CreateNeighList(&mut boundary_heap, 1, 1, boundary_atoms.as_const(), 0, None).unwrap();
        let boundary_pointer = boundary_heap.slice(boundary.as_const()).unwrap()[0];
        assert_eq!(
            boundary_heap.slice(boundary_pointer.as_const()).unwrap()[0],
            20
        );
        assert_eq!(
            &boundary_heap.slice(boundary_pointer.as_const()).unwrap()[1..21],
            &(0_u16..20).collect::<Vec<_>>()
        );
        FreeNeighList(&mut boundary_heap, boundary).unwrap();

        let mut zero_heap = SourceHeap::default();
        zero_heap.trace_source_allocations();
        let zero =
            CreateNeighList(&mut zero_heap, 0, 0, SourceConstPointer::null(), 0, None).unwrap();
        assert!(!zero.is_null());
        assert!(zero_heap.slice(zero.as_const()).unwrap()[0].is_null());
        assert_eq!(zero_heap.source_allocation_calls(), 2);

        for successful_allocations in [0, 1] {
            let mut failure_heap = SourceHeap::default();
            let failure_atoms = failure_heap.allocate_model_storage(atoms.clone()).unwrap();
            let baseline = failure_heap.live_allocation_count();
            failure_heap.fail_after_allocations(successful_allocations);
            assert!(
                CreateNeighList(&mut failure_heap, 3, 3, failure_atoms.as_const(), 0, None,)
                    .unwrap()
                    .is_null()
            );
            assert_eq!(
                failure_heap.source_allocation_calls(),
                successful_allocations + 1
            );
            assert_eq!(failure_heap.live_allocation_count(), baseline);
        }
    }

    #[test]
    fn source_port__ichisort__freeneighlist__line_941() {
        let mut heap = SourceHeap::default();
        assert_eq!(FreeNeighList(&mut heap, SourceMutPointer::null()), Ok(()));

        let pointer_list = heap
            .allocate_model_storage(vec![SourceMutPointer::<AT_NUMB>::null()])
            .unwrap();
        assert_eq!(FreeNeighList(&mut heap, pointer_list), Ok(()));
        assert_eq!(
            heap.slice(pointer_list.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let atom_list = heap.allocate_model_storage(vec![1_u16, 2, 3]).unwrap();
        let pointer_list = heap
            .allocate_model_storage(vec![atom_list, atom_list.offset(2).unwrap()])
            .unwrap();
        assert_eq!(FreeNeighList(&mut heap, pointer_list), Ok(()));
        assert_eq!(
            heap.slice(atom_list.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(pointer_list.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
    }
}
