use std::mem::size_of;
use std::sync::atomic::{AtomicU16, Ordering};

use crate::source::base::util::inchi_free;
use crate::source_types::{
    AT_NUMB, CANON_GLOBALS, NodeSet, SourceHeap, SourceHeapError, SourceMutPointer, Vertex, bitWord,
};

static RANK_MARK_BIT: AtomicU16 = AtomicU16::new(0);
static RANK_MASK_BIT: AtomicU16 = AtomicU16::new(u16::MAX);
const SOURCE_SIZEOF_POINTER: u64 = 8;

pub(crate) fn rank_mark_bit() -> AT_NUMB {
    RANK_MARK_BIT.load(Ordering::Relaxed)
}

pub(crate) fn rank_mask_bit() -> AT_NUMB {
    RANK_MASK_BIT.load(Ordering::Relaxed)
}

#[allow(non_snake_case)]
pub(crate) fn NodeSetCreate(
    heap: &mut SourceHeap,
    pCG: &CANON_GLOBALS,
    pSet: &mut NodeSet,
    n: i32,
    L: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:718 NodeSetCreate
    // INCHI✔️✔️: int NodeSetCreate( struct tagCANON_GLOBALS *pCG,
    // INCHI✔️✔️:                    NodeSet *pSet,
    // INCHI✔️✔️:                    int n,
    // INCHI✔️✔️:                    int L )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, len;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     len = ( n + pCG->m_num_bit - 1 ) / pCG->m_num_bit;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     pSet->bitword = (bitWord**) inchi_calloc( L, sizeof( pSet->bitword[0] ) );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (!pSet->bitword)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     pSet->bitword[0] = (bitWord*) inchi_calloc( (long long)len * (long long)L, sizeof( pSet->bitword[0][0] ) ); /* djb-rwth: cast operators added */
    // INCHI✔️✔️:     if (!pSet->bitword[0])
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* Cleanup */
    // INCHI✔️✔️:         inchi_free( pSet->bitword );
    // INCHI✔️✔️:         pSet->bitword = NULL;
    // INCHI✔️✔️:         return 0; /* failed */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     for (i = 1; i < L; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         pSet->bitword[i] = pSet->bitword[i - 1] + len;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     pSet->len_set = len;
    // INCHI✔️✔️:     pSet->num_set = L;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 1;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: NodeSetCreate

    let len = n
        .wrapping_add(pCG.m_num_bit)
        .wrapping_sub(1)
        .wrapping_div(pCG.m_num_bit);

    pSet.bitword = match crate::source::base::util::inchi_calloc::<SourceMutPointer<bitWord>>(
        heap,
        L as u64,
        SOURCE_SIZEOF_POINTER,
    ) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => {
            pSet.bitword = SourceMutPointer::null();
            return Ok(0);
        }
        Err(error) => return Err(error),
    };

    let storage_count = i64::from(len)
        .checked_mul(i64::from(L))
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    let storage = match u64::try_from(storage_count)
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)
        .and_then(|count| {
            crate::source::base::util::inchi_calloc::<bitWord>(
                heap,
                count,
                size_of::<bitWord>() as u64,
            )
        }) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => {
            inchi_free(heap, pSet.bitword)?;
            pSet.bitword = SourceMutPointer::null();
            return Ok(0);
        }
        Err(error) => return Err(error),
    };

    {
        let rows = heap.slice_mut(pSet.bitword)?;
        let row_count =
            usize::try_from(L).map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        rows.get_mut(0)
            .ok_or(SourceHeapError::PointerOutOfBounds)
            .map(|slot| *slot = storage)?;
        let mut i = 1_i32;
        while i < L {
            let previous = rows
                .get(usize::try_from(i - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let current = previous.offset(i64::from(len))?;
            rows.get_mut(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)
                .map(|slot| *slot = current)?;
            i = i.wrapping_add(1);
        }
        let _ = rows
            .get(..row_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
    }

    pSet.len_set = len;
    pSet.num_set = L;
    Ok(1)
}

#[allow(non_snake_case)]
pub(crate) fn NodeSetFree(
    heap: &mut SourceHeap,
    _pCG: &CANON_GLOBALS,
    pSet: &mut NodeSet,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:754 NodeSetFree
    // INCHI✔️✔️: void NodeSetFree( struct tagCANON_GLOBALS *pCG, NodeSet *pSet )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     if (pSet && pSet->bitword)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (pSet->bitword[0])
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             inchi_free( pSet->bitword[0] );
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         inchi_free( pSet->bitword );
    // INCHI✔️✔️:         pSet->bitword = NULL;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: NodeSetFree

    if pSet.bitword.is_null() {
        return Ok(());
    }
    let first = pointer_array_get(heap, pSet.bitword, 0)?;
    if !first.is_null() {
        inchi_free(heap, first)?;
    }
    inchi_free(heap, pSet.bitword)?;
    pSet.bitword = SourceMutPointer::null();
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn SetBitFree(
    heap: &mut SourceHeap,
    globals: &mut CANON_GLOBALS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3259 SetBitFree
    // INCHI✔️❌: int SetBitFree( CANON_GLOBALS *pCG )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (pCG->m_bBit)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( pCG->m_bBit );
    // INCHI✔️❌:         pCG->m_bBit = NULL;
    // INCHI✔️❌:         INCHI_HEAPCHK
    // INCHI✔️❌:         return 1; /* success */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0; /* already destroyed */
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: SetBitFree

    if !globals.m_bBit.is_null() {
        inchi_free(heap, globals.m_bBit)?;
        globals.m_bBit = SourceMutPointer::null();
        return Ok(1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn SetBitCreate(
    heap: &mut SourceHeap,
    globals: &mut CANON_GLOBALS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3205 SetBitCreate
    // INCHI✔️✔️: int SetBitCreate( CANON_GLOBALS *pCG )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     bitWord b1 = 1, b2;
    // INCHI✔️✔️:     AT_NUMB n1, n2;
    // INCHI✔️✔️: #ifdef INCHI_CANON_USE_HASH
    // INCHI✔️❌:     CtHash  h1, h2;
    // INCHI✔️❌: #endif
    // INCHI✔️✔️:     int    i;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (pCG->m_bBitInitialized)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         INCHI_HEAPCHK
    // INCHI✔️✔️:         return 0; /* already created */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* djb-rwth: removing redundant code */
    // INCHI✔️✔️:     pCG->m_num_bit = 1;
    // INCHI✔️✔️:     for (b1 = 1, pCG->m_num_bit = 1; b1 < ( b2 = (bitWord) ( ( b1 << 1 )& BIT_WORD_MASK ) ); b1 = b2, pCG->m_num_bit++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     pCG->m_bBit = (bitWord*) inchi_calloc( pCG->m_num_bit, sizeof( bitWord ) );
    // INCHI✔️✔️:     if (!pCG->m_bBit)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         INCHI_HEAPCHK
    // INCHI✔️✔️:             return -1; /* failed */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     for (i = 0, b1 = 1; i < pCG->m_num_bit; i++, b1 <<= 1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         pCG->m_bBit[i] = b1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (n1 = 1; n1 < ( n2 = (AT_RANK) ( ( n1 << 1 )& AT_RANK_MASK ) ); n1 = n2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     rank_mark_bit = n1;
    // INCHI✔️✔️:     rank_mask_bit = ~n1;
    // INCHI✔️✔️:
    // INCHI✔️❌: #ifdef INCHI_CANON_USE_HASH
    // INCHI✔️❌:     for (h1 = 1; h1 < ( h2 = ( h1 << 1 ) ); h1 = h2)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     hash_mark_bit = h1;
    // INCHI✔️❌: #endif
    // INCHI✔️✔️:     pCG->m_bBitInitialized = 1;
    // INCHI✔️✔️:     INCHI_HEAPCHK
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 1;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: SetBitCreate

    if globals.m_bBitInitialized != 0 {
        return Ok(0);
    }

    globals.m_num_bit = 1;
    let mut b1: bitWord = 1;
    loop {
        let b2 = b1.wrapping_shl(1) & u16::MAX;
        if b1 >= b2 {
            break;
        }
        globals.m_num_bit += 1;
        b1 = b2;
    }

    let created = match crate::source::base::util::inchi_calloc::<bitWord>(
        heap,
        globals.m_num_bit as u64,
        size_of::<bitWord>() as u64,
    ) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => {
            globals.m_bBit = SourceMutPointer::null();
            return Ok(-1);
        }
        Err(error) => return Err(error),
    };
    globals.m_bBit = created;

    {
        let bits = heap.slice_mut(globals.m_bBit)?;
        let mut mask: bitWord = 1;
        for index in 0..usize::try_from(globals.m_num_bit)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?
        {
            bits[index] = mask;
            mask = mask.wrapping_shl(1);
        }
    }

    let mut n1: AT_NUMB = 1;
    loop {
        let n2 = n1.wrapping_shl(1) & u16::MAX;
        if n1 >= n2 {
            break;
        }
        n1 = n2;
    }
    RANK_MARK_BIT.store(n1, Ordering::Relaxed);
    RANK_MASK_BIT.store(!n1, Ordering::Relaxed);

    globals.m_bBitInitialized = 1;
    Ok(1)
}

fn pointer_array_get<T: 'static>(
    heap: &SourceHeap,
    pointer: SourceMutPointer<SourceMutPointer<T>>,
    index: i32,
) -> Result<SourceMutPointer<T>, SourceHeapError> {
    let offset = usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    heap.slice(pointer.as_const())?
        .get(offset)
        .copied()
        .ok_or(SourceHeapError::PointerOutOfBounds)
}

fn source_get<T: Copy + 'static>(
    heap: &SourceHeap,
    pointer: SourceMutPointer<T>,
    index: i32,
) -> Result<T, SourceHeapError> {
    let offset = usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    heap.slice(pointer.as_const())?
        .get(offset)
        .copied()
        .ok_or(SourceHeapError::PointerOutOfBounds)
}

#[allow(non_snake_case)]
pub(crate) fn NodeSetFromRadEndpoints(
    heap: &mut SourceHeap,
    pCG: &CANON_GLOBALS,
    cur_nodes: &NodeSet,
    k: i32,
    RadEndpoints: SourceMutPointer<Vertex>,
    num_v: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1138 NodeSetFromRadEndpoints
    // INCHI✔️❌: void NodeSetFromRadEndpoints( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                               NodeSet *cur_nodes,
    // INCHI✔️❌:                               int k, /*Node *v*/
    // INCHI✔️❌:                               Vertex RadEndpoints[],
    // INCHI✔️❌:                               int num_v )
    // INCHI✔️❌: {
    // INCHI✔️❌:     bitWord *Bits = cur_nodes->bitword[k];
    // INCHI✔️❌:     int      len = cur_nodes->len_set * sizeof( bitWord );
    // INCHI✔️❌:     int      i, j;
    // INCHI✔️❌:
    // INCHI✔️❌:     memset( Bits, 0, len ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 1; i < num_v; i += 2)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         j = (int) RadEndpoints[i];
    // INCHI✔️❌:         Bits[j / pCG->m_num_bit] |= pCG->m_bBit[j % pCG->m_num_bit];
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: NodeSetFromRadEndpoints
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: NodeSetFromRadEndpoints
    // INCHI✔️❌: typedef U_SHORT bitWord
    // END INCHI ACTIVE MACRO CONFIGURATION: NodeSetFromRadEndpoints

    let bits = pointer_array_get(heap, cur_nodes.bitword, k)?;
    let len_set =
        usize::try_from(cur_nodes.len_set).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    {
        let target = heap.slice_mut(bits)?;
        let clear_len = target
            .get_mut(..len_set)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        clear_len.fill(0);
    }

    let mut i = 1_i32;
    while i < num_v {
        let j = source_get(heap, RadEndpoints, i)?;
        let word_index = j
            .checked_div(pCG.m_num_bit)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let bit_index = j
            .checked_rem(pCG.m_num_bit)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mask = source_get(heap, pCG.m_bBit, bit_index)?;
        let word_offset =
            usize::try_from(word_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let target = heap.slice_mut(bits)?;
        let word = target
            .get_mut(word_offset)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        *word |= mask;
        i = i.wrapping_add(2);
    }

    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn RemoveFromNodeSet(
    heap: &mut SourceHeap,
    pCG: &CANON_GLOBALS,
    cur_nodes: &NodeSet,
    k: i32,
    v: SourceMutPointer<Vertex>,
    num_v: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1159 RemoveFromNodeSet
    // INCHI✔️❌: void RemoveFromNodeSet( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                         NodeSet *cur_nodes, int k, Vertex v[], int num_v )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (cur_nodes->bitword)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         bitWord *Bits = cur_nodes->bitword[k];
    // INCHI✔️❌:         /*int      len  = cur_nodes->len_set*sizeof(bitWord);*/
    // INCHI✔️❌:         int      i, j;
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = 0; i < num_v; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             j = (int) v[i];
    // INCHI✔️❌:             Bits[j / pCG->m_num_bit] &= ~pCG->m_bBit[j % pCG->m_num_bit];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: RemoveFromNodeSet
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: RemoveFromNodeSet
    // INCHI✔️❌: typedef U_SHORT bitWord
    // END INCHI ACTIVE MACRO CONFIGURATION: RemoveFromNodeSet

    if cur_nodes.bitword.is_null() {
        return Ok(());
    }
    let bits = pointer_array_get(heap, cur_nodes.bitword, k)?;
    let mut i = 0_i32;
    while i < num_v {
        let j = source_get(heap, v, i)?;
        let word_index = j
            .checked_div(pCG.m_num_bit)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let bit_index = j
            .checked_rem(pCG.m_num_bit)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mask = source_get(heap, pCG.m_bBit, bit_index)?;
        let word_offset =
            usize::try_from(word_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let target = heap.slice_mut(bits)?;
        let word = target
            .get_mut(word_offset)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        *word &= !mask;
        i = i.wrapping_add(1);
    }

    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn DoNodeSetsIntersect(
    heap: &SourceHeap,
    cur_nodes: &NodeSet,
    k1: i32,
    k2: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1178 DoNodeSetsIntersect
    // INCHI✔️❌: int DoNodeSetsIntersect( NodeSet *cur_nodes, int k1, int k2 )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (cur_nodes->bitword)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         bitWord *Bits1 = cur_nodes->bitword[k1];
    // INCHI✔️❌:         bitWord *Bits2 = cur_nodes->bitword[k2];
    // INCHI✔️❌:         int      len = cur_nodes->len_set;
    // INCHI✔️❌:         int      i;
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = 0; i < len; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (Bits1[i] & Bits2[i])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: DoNodeSetsIntersect

    if cur_nodes.bitword.is_null() {
        return Ok(0);
    }
    let bits1 = pointer_array_get(heap, cur_nodes.bitword, k1)?;
    let bits2 = pointer_array_get(heap, cur_nodes.bitword, k2)?;
    let len =
        usize::try_from(cur_nodes.len_set).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let slice1 = heap.slice(bits1.as_const())?;
    let slice2 = heap.slice(bits2.as_const())?;
    for i in 0..len {
        if slice1
            .get(i)
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            & slice2
                .get(i)
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
            != 0
        {
            return Ok(1);
        }
    }

    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn IsNodeSetEmpty(
    heap: &SourceHeap,
    cur_nodes: &NodeSet,
    k: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1201 IsNodeSetEmpty
    // INCHI✔️❌: int IsNodeSetEmpty( NodeSet *cur_nodes, int k )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (cur_nodes->bitword)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         bitWord *Bits = cur_nodes->bitword[k];
    // INCHI✔️❌:         int      len = cur_nodes->len_set;
    // INCHI✔️❌:         int      i;
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = 0; i < len; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (Bits[i])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 1;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: IsNodeSetEmpty

    if cur_nodes.bitword.is_null() {
        return Ok(1);
    }
    let bits = pointer_array_get(heap, cur_nodes.bitword, k)?;
    let len =
        usize::try_from(cur_nodes.len_set).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let slice = heap.slice(bits.as_const())?;
    for i in 0..len {
        if slice
            .get(i)
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            != 0
        {
            return Ok(0);
        }
    }

    Ok(1)
}

#[allow(non_snake_case)]
pub(crate) fn AddNodeSet2ToNodeSet1(
    heap: &mut SourceHeap,
    cur_nodes: &NodeSet,
    k1: i32,
    k2: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1223 AddNodeSet2ToNodeSet1
    // INCHI✔️❌: void AddNodeSet2ToNodeSet1( NodeSet *cur_nodes, int k1, int k2 )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (cur_nodes->bitword)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         bitWord *Bits1 = cur_nodes->bitword[k1];
    // INCHI✔️❌:         bitWord *Bits2 = cur_nodes->bitword[k2];
    // INCHI✔️❌:         int      len = cur_nodes->len_set;
    // INCHI✔️❌:         int      i;
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = 0; i < len; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             Bits1[i] |= Bits2[i];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: AddNodeSet2ToNodeSet1

    if cur_nodes.bitword.is_null() {
        return Ok(());
    }
    let bits1 = pointer_array_get(heap, cur_nodes.bitword, k1)?;
    let bits2 = pointer_array_get(heap, cur_nodes.bitword, k2)?;
    let len = cur_nodes.len_set;
    let mut i = 0_i32;
    while i < len {
        let source = source_get(heap, bits2, i)?;
        let offset = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let target = heap.slice_mut(bits1)?;
        let word = target
            .get_mut(offset)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        *word |= source;
        i = i.wrapping_add(1);
    }

    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn AddNodesToRadEndpoints(
    heap: &mut SourceHeap,
    pCG: &CANON_GLOBALS,
    cur_nodes: &NodeSet,
    k: i32,
    RadEndpoints: SourceMutPointer<Vertex>,
    vRad: Vertex,
    nStart: i32,
    nLen: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1241 AddNodesToRadEndpoints
    // INCHI✔️❌: int AddNodesToRadEndpoints( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                             NodeSet *cur_nodes,
    // INCHI✔️❌:                             int k,
    // INCHI✔️❌:                             Vertex RadEndpoints[],
    // INCHI✔️❌:                             Vertex vRad,
    // INCHI✔️❌:                             int nStart, int nLen )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int n = nStart;
    // INCHI✔️❌:     if (cur_nodes->bitword)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         bitWord *Bits = cur_nodes->bitword[k];
    // INCHI✔️❌:         int      len = cur_nodes->len_set;
    // INCHI✔️❌:         int      i, j;
    // INCHI✔️❌:         Vertex   v;
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = 0, v = 0; i < len; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (Bits[i])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (j = 0; j < pCG->m_num_bit; j++, v++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (Bits[i] & pCG->m_bBit[j])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (n >= nLen)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             return -1; /* overflow */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         RadEndpoints[n++] = vRad;
    // INCHI✔️❌:                         RadEndpoints[n++] = v;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 v += pCG->m_num_bit;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return n;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: AddNodesToRadEndpoints
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: AddNodesToRadEndpoints
    // INCHI✔️❌: typedef U_SHORT bitWord
    // END INCHI ACTIVE MACRO CONFIGURATION: AddNodesToRadEndpoints

    let mut n = nStart;
    if cur_nodes.bitword.is_null() {
        return Ok(n);
    }

    let bits = pointer_array_get(heap, cur_nodes.bitword, k)?;
    let len = cur_nodes.len_set;
    let mut i = 0_i32;
    let mut v: Vertex = 0;
    while i < len {
        let word = source_get(heap, bits, i)?;
        if word != 0 {
            let mut j = 0_i32;
            while j < pCG.m_num_bit {
                let mask = source_get(heap, pCG.m_bBit, j)?;
                if word & mask != 0 {
                    if n >= nLen {
                        return Ok(-1);
                    }
                    let first_offset =
                        usize::try_from(n).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    heap.slice_mut(RadEndpoints)?
                        .get_mut(first_offset)
                        .ok_or(SourceHeapError::PointerOutOfBounds)
                        .map(|slot| *slot = vRad)?;
                    n = n.wrapping_add(1);
                    let second_offset =
                        usize::try_from(n).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    heap.slice_mut(RadEndpoints)?
                        .get_mut(second_offset)
                        .ok_or(SourceHeapError::PointerOutOfBounds)
                        .map(|slot| *slot = v)?;
                    n = n.wrapping_add(1);
                }
                j = j.wrapping_add(1);
                v = v.wrapping_add(1);
            }
        } else {
            v = v.wrapping_add(pCG.m_num_bit);
        }
        i = i.wrapping_add(1);
    }

    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_port__ichican2__nodesetcreate__line_718() {
        let globals = CANON_GLOBALS {
            m_num_bit: 16,
            ..CANON_GLOBALS::default()
        };
        let mut heap = SourceHeap::default();
        let mut set = NodeSet::default();

        assert_eq!(NodeSetCreate(&mut heap, &globals, &mut set, 33, 3), Ok(1));
        assert_eq!(set.len_set, 3);
        assert_eq!(set.num_set, 3);
        let rows = heap.slice(set.bitword.as_const()).unwrap().to_vec();
        assert_eq!(rows.len(), 3);
        assert_eq!(heap.slice(rows[0].as_const()).unwrap(), &[0; 9]);
        assert_eq!(heap.slice(rows[1].as_const()).unwrap(), &[0; 6]);
        assert_eq!(heap.slice(rows[2].as_const()).unwrap(), &[0; 3]);
        heap.slice_mut(rows[1]).unwrap()[0] = 0xabcd;
        heap.slice_mut(rows[2]).unwrap()[2] = 0x1234;
        assert_eq!(
            heap.slice(rows[0].as_const()).unwrap(),
            &[0, 0, 0, 0xabcd, 0, 0, 0, 0, 0x1234]
        );

        let mut first_failure_heap = SourceHeap::default();
        first_failure_heap.fail_after_allocations(0);
        let mut first_failure = NodeSet {
            len_set: 7,
            num_set: 8,
            ..NodeSet::default()
        };
        assert_eq!(
            NodeSetCreate(&mut first_failure_heap, &globals, &mut first_failure, 33, 3),
            Ok(0)
        );
        assert_eq!(first_failure_heap.source_allocation_calls(), 1);
        assert!(first_failure.bitword.is_null());
        assert_eq!(first_failure.len_set, 7);
        assert_eq!(first_failure.num_set, 8);

        let mut second_failure_heap = SourceHeap::default();
        second_failure_heap.fail_after_allocations(1);
        let mut second_failure = NodeSet {
            len_set: 9,
            num_set: 10,
            ..NodeSet::default()
        };
        assert_eq!(
            NodeSetCreate(
                &mut second_failure_heap,
                &globals,
                &mut second_failure,
                33,
                3
            ),
            Ok(0)
        );
        assert_eq!(second_failure_heap.source_allocation_calls(), 2);
        assert!(second_failure.bitword.is_null());
        assert_eq!(second_failure.len_set, 9);
        assert_eq!(second_failure.num_set, 10);
    }

    #[test]
    fn source_port__ichican2__nodesetfree__line_754() {
        let globals = CANON_GLOBALS {
            m_num_bit: 16,
            ..CANON_GLOBALS::default()
        };
        let mut heap = SourceHeap::default();
        let mut set = NodeSet::default();
        assert_eq!(NodeSetCreate(&mut heap, &globals, &mut set, 20, 2), Ok(1));
        let outer = set.bitword;
        let rows = heap.slice(set.bitword.as_const()).unwrap().to_vec();
        assert_eq!(rows.len(), 2);

        assert_eq!(NodeSetFree(&mut heap, &globals, &mut set), Ok(()));
        assert!(set.bitword.is_null());
        assert_eq!(set.len_set, 2);
        assert_eq!(set.num_set, 2);
        assert_eq!(
            heap.slice(outer.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(rows[0].as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );

        let null_row_outer = heap
            .allocate(vec![SourceMutPointer::<bitWord>::null()])
            .unwrap();
        let mut null_row_set = NodeSet {
            bitword: null_row_outer,
            len_set: 5,
            num_set: 1,
        };
        assert_eq!(NodeSetFree(&mut heap, &globals, &mut null_row_set), Ok(()));
        assert!(null_row_set.bitword.is_null());
        assert_eq!(
            heap.slice(null_row_outer.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );

        let mut empty = NodeSet {
            len_set: 3,
            num_set: 4,
            ..NodeSet::default()
        };
        assert_eq!(NodeSetFree(&mut heap, &globals, &mut empty), Ok(()));
        assert!(empty.bitword.is_null());
        assert_eq!(empty.len_set, 3);
        assert_eq!(empty.num_set, 4);
    }

    #[test]
    fn source_port__ichican2__setbitcreate__line_3205() {
        let mut heap = SourceHeap::default();
        let mut globals = CANON_GLOBALS::default();

        assert_eq!(SetBitCreate(&mut heap, &mut globals), Ok(1));
        assert_eq!(globals.m_bBitInitialized, 1);
        assert_eq!(globals.m_num_bit, 16);
        assert_eq!(
            heap.slice(globals.m_bBit.as_const()).unwrap(),
            &[
                0x0001, 0x0002, 0x0004, 0x0008, 0x0010, 0x0020, 0x0040, 0x0080, 0x0100, 0x0200,
                0x0400, 0x0800, 0x1000, 0x2000, 0x4000, 0x8000,
            ]
        );
        assert_eq!(rank_mark_bit(), 0x8000);
        assert_eq!(rank_mask_bit(), 0x7fff);

        let created = globals.m_bBit;
        heap.trace_source_allocations();
        assert_eq!(SetBitCreate(&mut heap, &mut globals), Ok(0));
        assert_eq!(heap.source_allocation_calls(), 0);
        assert_eq!(globals.m_bBit, created);
        assert_eq!(globals.m_num_bit, 16);
        assert_eq!(globals.m_bBitInitialized, 1);

        let mut failing_heap = SourceHeap::default();
        failing_heap.fail_after_allocations(0);
        let mut failing_globals = CANON_GLOBALS {
            m_bBit: created,
            ..CANON_GLOBALS::default()
        };
        assert_eq!(
            SetBitCreate(&mut failing_heap, &mut failing_globals),
            Ok(-1)
        );
        assert_eq!(failing_heap.source_allocation_calls(), 1);
        assert!(failing_globals.m_bBit.is_null());
        assert_eq!(failing_globals.m_num_bit, 16);
        assert_eq!(failing_globals.m_bBitInitialized, 0);
    }

    #[test]
    fn source_port__ichican2__setbitfree__line_3259() {
        let mut heap = SourceHeap::default();
        let bits = heap.allocate(vec![1_u16, 2, 4, 8]).unwrap();
        let mut globals = CANON_GLOBALS {
            m_bBit: bits,
            m_bBitInitialized: 7,
            m_num_bit: 4,
            ..CANON_GLOBALS::default()
        };
        assert_eq!(SetBitFree(&mut heap, &mut globals), Ok(1));
        assert!(globals.m_bBit.is_null());
        assert_eq!(globals.m_bBitInitialized, 7);
        assert_eq!(globals.m_num_bit, 4);
        assert_eq!(
            heap.slice(bits.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(SetBitFree(&mut heap, &mut globals), Ok(0));
        assert_eq!(globals.m_bBitInitialized, 7);
        assert_eq!(globals.m_num_bit, 4);
    }

    #[test]
    fn source_port__ichican2__nodesetfromradendpoints__line_1138() {
        let mut heap = SourceHeap::default();
        let bit_masks = heap
            .allocate((0..16).map(|bit| 1_u16 << bit).collect())
            .unwrap();
        let globals = CANON_GLOBALS {
            m_bBit: bit_masks,
            m_num_bit: 16,
            ..CANON_GLOBALS::default()
        };
        let set0 = heap.allocate(vec![0xaaaa_u16, 0xbbbb]).unwrap();
        let set1 = heap.allocate(vec![0xffff_u16, 0xffff]).unwrap();
        let bitword = heap.allocate(vec![set0, set1]).unwrap();
        let cur_nodes = NodeSet {
            bitword,
            num_set: 2,
            len_set: 2,
        };
        let endpoints = heap.allocate(vec![99, 0, 88, 17, 77, 31]).unwrap();

        assert_eq!(
            NodeSetFromRadEndpoints(&mut heap, &globals, &cur_nodes, 1, endpoints, 6),
            Ok(())
        );
        assert_eq!(heap.slice(set0.as_const()).unwrap(), &[0xaaaa, 0xbbbb]);
        assert_eq!(heap.slice(set1.as_const()).unwrap(), &[0x0001, 0x8002]);

        let clear_only = heap.allocate(vec![12, 13]).unwrap();
        assert_eq!(
            NodeSetFromRadEndpoints(&mut heap, &globals, &cur_nodes, 1, clear_only, 1),
            Ok(())
        );
        assert_eq!(heap.slice(set1.as_const()).unwrap(), &[0x0000, 0x0000]);
    }

    #[test]
    fn source_port__ichican2__removefromnodeset__line_1159() {
        let mut heap = SourceHeap::default();
        let bit_masks = heap
            .allocate((0..16).map(|bit| 1_u16 << bit).collect())
            .unwrap();
        let globals = CANON_GLOBALS {
            m_bBit: bit_masks,
            m_num_bit: 16,
            ..CANON_GLOBALS::default()
        };
        let set0 = heap.allocate(vec![0xffff_u16, 0xffff]).unwrap();
        let set1 = heap.allocate(vec![0xffff_u16, 0xffff]).unwrap();
        let bitword = heap.allocate(vec![set0, set1]).unwrap();
        let cur_nodes = NodeSet {
            bitword,
            num_set: 2,
            len_set: 2,
        };
        let atoms = heap.allocate(vec![0, 17, 31]).unwrap();

        assert_eq!(
            RemoveFromNodeSet(&mut heap, &globals, &cur_nodes, 1, atoms, 3),
            Ok(())
        );
        assert_eq!(heap.slice(set1.as_const()).unwrap(), &[0xfffe, 0x7ffd]);
        assert_eq!(heap.slice(set0.as_const()).unwrap(), &[0xffff, 0xffff]);

        let empty = NodeSet::default();
        assert_eq!(
            RemoveFromNodeSet(&mut heap, &globals, &empty, 0, atoms, 3),
            Ok(())
        );
    }

    #[test]
    fn source_port__ichican2__donodesetsintersect__line_1178() {
        let mut heap = SourceHeap::default();
        let set0 = heap.allocate(vec![0x0001_u16, 0x0002]).unwrap();
        let set1 = heap.allocate(vec![0x0004_u16, 0x0002]).unwrap();
        let set2 = heap.allocate(vec![0x0008_u16, 0x0010]).unwrap();
        let bitword = heap.allocate(vec![set0, set1, set2]).unwrap();
        let cur_nodes = NodeSet {
            bitword,
            num_set: 3,
            len_set: 2,
        };

        assert_eq!(DoNodeSetsIntersect(&heap, &cur_nodes, 0, 1), Ok(1));
        assert_eq!(DoNodeSetsIntersect(&heap, &cur_nodes, 0, 2), Ok(0));
        assert_eq!(DoNodeSetsIntersect(&heap, &NodeSet::default(), 0, 1), Ok(0));
    }

    #[test]
    fn source_port__ichican2__isnodesetempty__line_1201() {
        let mut heap = SourceHeap::default();
        let set0 = heap.allocate(vec![0x0000_u16, 0x0000]).unwrap();
        let set1 = heap.allocate(vec![0x0001_u16, 0x0000]).unwrap();
        let bitword = heap.allocate(vec![set0, set1]).unwrap();
        let cur_nodes = NodeSet {
            bitword,
            num_set: 2,
            len_set: 2,
        };

        assert_eq!(IsNodeSetEmpty(&heap, &cur_nodes, 0), Ok(1));
        assert_eq!(IsNodeSetEmpty(&heap, &cur_nodes, 1), Ok(0));
        assert_eq!(IsNodeSetEmpty(&heap, &NodeSet::default(), 0), Ok(1));
    }

    #[test]
    fn source_port__ichican2__addnodeset2tonodeset1__line_1223() {
        let mut heap = SourceHeap::default();
        let set0 = heap.allocate(vec![0x0001_u16, 0x00f0]).unwrap();
        let set1 = heap.allocate(vec![0x0004_u16, 0x0f00]).unwrap();
        let set2 = heap.allocate(vec![0x8000_u16, 0x0008]).unwrap();
        let bitword = heap.allocate(vec![set0, set1, set2]).unwrap();
        let cur_nodes = NodeSet {
            bitword,
            num_set: 3,
            len_set: 2,
        };

        assert_eq!(AddNodeSet2ToNodeSet1(&mut heap, &cur_nodes, 0, 1), Ok(()));
        assert_eq!(heap.slice(set0.as_const()).unwrap(), &[0x0005, 0x0ff0]);
        assert_eq!(heap.slice(set1.as_const()).unwrap(), &[0x0004, 0x0f00]);
        assert_eq!(heap.slice(set2.as_const()).unwrap(), &[0x8000, 0x0008]);

        assert_eq!(
            AddNodeSet2ToNodeSet1(&mut heap, &NodeSet::default(), 0, 1),
            Ok(())
        );
    }

    #[test]
    fn source_port__ichican2__addnodestoradendpoints__line_1241() {
        let mut heap = SourceHeap::default();
        let bit_masks = heap
            .allocate((0..16).map(|bit| 1_u16 << bit).collect())
            .unwrap();
        let globals = CANON_GLOBALS {
            m_bBit: bit_masks,
            m_num_bit: 16,
            ..CANON_GLOBALS::default()
        };
        let set0 = heap.allocate(vec![0x0000_u16, 0x8003]).unwrap();
        let set1 = heap.allocate(vec![0x0004_u16, 0x0000]).unwrap();
        let bitword = heap.allocate(vec![set0, set1]).unwrap();
        let cur_nodes = NodeSet {
            bitword,
            num_set: 2,
            len_set: 2,
        };
        let endpoints = heap.allocate(vec![-7; 10]).unwrap();

        assert_eq!(
            AddNodesToRadEndpoints(&mut heap, &globals, &cur_nodes, 0, endpoints, 42, 2, 10),
            Ok(8)
        );
        assert_eq!(
            heap.slice(endpoints.as_const()).unwrap(),
            &[-7, -7, 42, 16, 42, 17, 42, 31, -7, -7]
        );

        assert_eq!(
            AddNodesToRadEndpoints(&mut heap, &globals, &cur_nodes, 1, endpoints, 77, 8, 10),
            Ok(10)
        );
        assert_eq!(
            heap.slice(endpoints.as_const()).unwrap(),
            &[-7, -7, 42, 16, 42, 17, 42, 31, 77, 2]
        );

        assert_eq!(
            AddNodesToRadEndpoints(&mut heap, &globals, &cur_nodes, 1, endpoints, 99, 10, 10),
            Ok(-1)
        );
        assert_eq!(
            AddNodesToRadEndpoints(
                &mut heap,
                &globals,
                &NodeSet::default(),
                0,
                endpoints,
                99,
                12,
                10
            ),
            Ok(12)
        );
    }
}
