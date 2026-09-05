use std::mem::{MaybeUninit, size_of};
use std::sync::atomic::{AtomicU16, Ordering};

use crate::source::base::ichisort::{
    AtomInvariant2SortWorkspace, CompChemElemLex, CreateNeighList, FreeNeighList, inchi_qsort,
};
use crate::source::base::{
    ichicano::FixCanonEquivalenceInfo,
    ichiisot::make_iso_sort_key,
    ichimap2::{DifferentiateRanks2, DifferentiateRanks4},
    util::{inchi_calloc, inchi_free},
};
use crate::source_types::local_ichican2::{
    BASE_H_NUMBER, Cell, ConTable, EMPTY_CT, EMPTY_H_NUMBER, INCHI_CANON_INFINITY, MAX_LAYERS,
    MAX_SET_SIZE, NORMALLY_ALLOWED_MAX_SET_SIZE, Node, Transposition, UnorderedPartition, kLeast,
};
use crate::source_types::{
    AT_FLAG_ISO_H_POINT, AT_ISO_SORT_KEY, AT_NUMB, AT_RANK, ATOM_INVARIANT2, ATOM_SIZES, BCN,
    CANON_COUNTS, CANON_DATA, CANON_FLAG_ISO_FIXED_H_DIFF, CANON_FLAG_ISO_ONLY_NON_TAUT_DIFF,
    CANON_FLAG_ISO_TAUT_DIFF, CANON_FLAG_NO_H_RECANON, CANON_FLAG_NO_TAUT_H_DIFF, CANON_GLOBALS,
    CT_CANON_ERR, CT_ERR_MAX, CT_ERR_MIN, CT_ISOCOUNT_ERR, CT_OUT_OF_RAM, CT_TIMEOUT_ERR,
    CT_UNKNOWN_ERR, CT_USER_QUIT_ERR, FLAG_NORM_CONSIDER_TAUT, INCHI_CLOCK, MAXVAL, NEIGH_LIST,
    NUM_CHEM_ELEMENTS, NUM_H, NodeSet, Partition, SourceConstPointer, SourceHeap, SourceHeapError,
    SourceMutPointer, StableSourceSlice, T_GROUP, T_GROUP_INFO, T_NUM_ISOTOPIC, T_NUM_NO_ISOTOPIC,
    TAUT_NON, TAUT_NUM, TAUT_YES, TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE, TG_FLAG_FOUND_ISOTOPIC_H_DONE,
    USER_ACTION_QUIT, Vertex, bitWord, clock_t, inchiTime, sp_ATOM,
    tagAtInvariantIndexes_AT_INV_HILL_ORDER, tagAtInvariantIndexes_AT_INV_NUM_CONNECTIONS,
    tagAtInvariantIndexes_AT_INV_NUM_H, tagAtInvariantIndexes_AT_INV_NUM_H_FIX,
    tagAtInvariantIndexes_AT_INV_NUM_TG_ENDPOINTS, tagAtInvariantIndexes_AT_INV_TAUT_ISO,
    tagAtInvariantIndexes_AT_INV_TG_NUMBERS,
};

static RANK_MARK_BIT: AtomicU16 = AtomicU16::new(0);
static RANK_MASK_BIT: AtomicU16 = AtomicU16::new(u16::MAX);
const SOURCE_SIZEOF_POINTER: u64 = 8;
const MAX_LAYERS_LEN: usize = MAX_LAYERS as usize;
type KLeastLayers = [kLeast; MAX_LAYERS_LEN];

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
    // INCHI✔️✔️:     len = ( n + pCG->m_num_bit - 1 ) / pCG->m_num_bit;
    let len = n
        .wrapping_add(pCG.m_num_bit)
        .wrapping_sub(1)
        .wrapping_div(pCG.m_num_bit);

    // INCHI✔️✔️:     pSet->bitword = (bitWord**) inchi_calloc( L, sizeof( pSet->bitword[0] ) );
    pSet.bitword = match crate::source::base::util::inchi_calloc::<SourceMutPointer<bitWord>>(
        heap,
        L as u64,
        SOURCE_SIZEOF_POINTER,
    ) {
        Ok(pointer) => pointer,
        // INCHI✔️✔️:     if (!pSet->bitword)
        // INCHI✔️✔️:     {
        // INCHI✔️✔️:         return 0;
        // INCHI✔️✔️:     }
        Err(SourceHeapError::AllocationFailed) => {
            pSet.bitword = SourceMutPointer::null();
            return Ok(0);
        }
        Err(error) => return Err(error),
    };

    let storage_count = i64::from(len)
        .checked_mul(i64::from(L))
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    // INCHI✔️✔️:     pSet->bitword[0] = (bitWord*) inchi_calloc(
    // INCHI✔️✔️:         (long long)len * (long long)L, sizeof( pSet->bitword[0][0] ) );
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
        // INCHI✔️✔️:     if (!pSet->bitword[0])
        // INCHI✔️✔️:     {
        // INCHI✔️✔️:         /* Cleanup */
        // INCHI✔️✔️:         inchi_free( pSet->bitword );
        // INCHI✔️✔️:         pSet->bitword = NULL;
        // INCHI✔️✔️:         return 0; /* failed */
        // INCHI✔️✔️:     }
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
        let rows = rows
            .get_mut(..row_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if rows.is_empty() {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let rows_pointer = rows.as_mut_ptr();

        // INCHI✔️✔️:     pSet->bitword[0] = (bitWord*) allocated storage;
        // SAFETY: the complete `L`-row prefix was validated above and remains
        // exclusively borrowed for the source loop below.
        unsafe { *rows_pointer = storage };
        let mut previous = storage;
        let row_length =
            u64::try_from(len).map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        // INCHI✔️✔️:     for (i = 1; i < L; i++)
        // INCHI✔️✔️:     {
        let mut i = 1_usize;
        while i < row_count {
            // INCHI✔️✔️:         pSet->bitword[i] = pSet->bitword[i - 1] + len;
            // SAFETY: storage_count == len * L was checked before allocating
            // the contiguous storage, and i remains smaller than L.
            previous = unsafe { previous.add_unchecked(row_length) };
            // SAFETY: `i < row_count`, and `rows_pointer` addresses the
            // validated, exclusively borrowed prefix.
            unsafe { *rows_pointer.add(i) = previous };
            i += 1;
        }
        // INCHI✔️✔️:     }
    }

    // INCHI✔️✔️:     pSet->len_set = len;
    pSet.len_set = len;
    // INCHI✔️✔️:     pSet->num_set = L;
    pSet.num_set = L;
    // INCHI✔️✔️:     return 1;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: NodeSetCreate
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

#[inline(always)]
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

fn source_clone<T: Clone + 'static>(
    heap: &SourceHeap,
    pointer: SourceMutPointer<T>,
    index: i32,
) -> Result<T, SourceHeapError> {
    let offset = usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    heap.slice(pointer.as_const())?
        .get(offset)
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)
}

#[inline(always)]
fn source_set<T: Copy + 'static>(
    heap: &mut SourceHeap,
    pointer: SourceMutPointer<T>,
    index: i32,
    value: T,
) -> Result<(), SourceHeapError> {
    let offset = usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    heap.slice_mut(pointer)?
        .get_mut(offset)
        .ok_or(SourceHeapError::PointerOutOfBounds)
        .map(|slot| *slot = value)
}

#[allow(non_snake_case)]
pub(crate) fn NodeSetFromVertices(
    heap: &mut SourceHeap,
    pCG: &CANON_GLOBALS,
    cur_nodes: &NodeSet,
    l: i32,
    v: SourceMutPointer<Node>,
    num_v: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1052 NodeSetFromVertices
    // INCHI✔️❌: void NodeSetFromVertices( CANON_GLOBALS *pCG, NodeSet *cur_nodes, int l, Node *v, int num_v )
    // INCHI✔️❌: {
    // INCHI✔️❌:     bitWord *Bits = cur_nodes->bitword[l - 1];
    // INCHI✔️❌:     int      len = cur_nodes->len_set * sizeof( bitWord );
    // INCHI✔️❌:     int      i, j;
    // INCHI✔️❌:
    // INCHI✔️❌:     memset( Bits, 0, len ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < num_v; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         j = (int) v[i] - 1;
    // INCHI✔️❌:         Bits[j / pCG->m_num_bit] |= pCG->m_bBit[j % pCG->m_num_bit];
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: NodeSetFromVertices
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: NodeSetFromVertices
    // INCHI✔️❌: typedef AT_NUMB Node
    // INCHI✔️❌: typedef U_SHORT bitWord
    // END INCHI ACTIVE MACRO CONFIGURATION: NodeSetFromVertices

    let bits = pointer_array_get(heap, cur_nodes.bitword, l.wrapping_sub(1))?;
    let len_set =
        usize::try_from(cur_nodes.len_set).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    heap.slice_mut(bits)?
        .get_mut(..len_set)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .fill(0);

    // SAFETY: CanonGraph owns the vertex array, global mask table, and NodeSet
    // row as independent fixed-size allocations for the complete source call.
    let vertices = unsafe { heap.stable_slice(v.as_const())? };
    let masks = unsafe { heap.stable_slice(pCG.m_bBit.as_const())? };
    let mut bit_words = unsafe { heap.stable_slice_mut(bits)? };
    // INCHI✔️✔️:     for (i = 0; i < num_v; i++)
    let mut i = 0_i32;
    while i < num_v {
        // INCHI✔️✔️:         j = (int)v[i] - 1;
        let j = i32::from(
            *vertices.get(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)?,
        )
        .wrapping_sub(1);
        let word_index = j
            .checked_div(pCG.m_num_bit)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let bit_index = j
            .checked_rem(pCG.m_num_bit)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        // INCHI✔️✔️:         Bits[j / pCG->m_num_bit] |= pCG->m_bBit[j % pCG->m_num_bit];
        let mask = *masks
            .get(usize::try_from(bit_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)?;
        let word = bit_words.get_mut(
            usize::try_from(word_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?,
        )?;
        *word |= mask;
        i = i.wrapping_add(1);
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn AllNodesAreInSet(
    heap: &SourceHeap,
    cur_nodes: &NodeSet,
    lcur_nodes: i32,
    set: &NodeSet,
    lset: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1071 AllNodesAreInSet
    // INCHI✔️❌: int AllNodesAreInSet( NodeSet *cur_nodes, int lcur_nodes, NodeSet *set, int lset )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i;
    // INCHI✔️❌:     int n = cur_nodes->len_set;
    // INCHI✔️❌:
    // INCHI✔️❌:     bitWord *BitsNode = cur_nodes->bitword[lcur_nodes - 1];
    // INCHI✔️❌:     bitWord *BitsSet = set->bitword[lset - 1];
    // INCHI✔️❌:
    // INCHI✔️❌:     /* find any BitsNode[i] bit not in BitsSet[i] */
    // INCHI✔️❌:     for (i = 0; i < n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (BitsNode[i] & ~BitsSet[i])
    // INCHI✔️❌:         {
    // INCHI✔️❌:
    // INCHI✔️❌:             INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:     return 1;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: AllNodesAreInSet
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: AllNodesAreInSet
    // INCHI✔️❌: typedef U_SHORT bitWord
    // END INCHI ACTIVE MACRO CONFIGURATION: AllNodesAreInSet

    let bits_node = pointer_array_get(heap, cur_nodes.bitword, lcur_nodes.wrapping_sub(1))?;
    let bits_set = pointer_array_get(heap, set.bitword, lset.wrapping_sub(1))?;
    // SAFETY: both NodeSet rows remain live and immutable for this source loop.
    let node_words = unsafe { heap.stable_slice(bits_node.as_const())? };
    let set_words = unsafe { heap.stable_slice(bits_set.as_const())? };
    // INCHI✔️✔️:     for (i = 0; i < n; i++)
    let mut i = 0_i32;
    while i < cur_nodes.len_set {
        let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let node_word = *node_words.get(index)?;
        let set_word = *set_words.get(index)?;
        // INCHI✔️✔️:         if (BitsNode[i] & ~BitsSet[i])
        if node_word & !set_word != 0 {
            // INCHI✔️✔️:             return 0;
            return Ok(0);
        }
        i = i.wrapping_add(1);
    }
    Ok(1)
}

#[allow(non_snake_case)]
pub(crate) fn PartitionGetMcrAndFixSet(
    heap: &mut SourceHeap,
    pCG: &CANON_GLOBALS,
    p: &Partition,
    Mcr: &NodeSet,
    Fix: &NodeSet,
    n: i32,
    l: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1098 PartitionGetMcrAndFixSet
    // INCHI✔️❌: void PartitionGetMcrAndFixSet( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                                Partition *p,
    // INCHI✔️❌:                                NodeSet *Mcr,
    // INCHI✔️❌:                                NodeSet *Fix,
    // INCHI✔️❌:                                int n,
    // INCHI✔️❌:                                int l )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j1, j2;
    // INCHI✔️❌:     AT_RANK r, r1;
    // INCHI✔️❌:     bitWord *McrBits = Mcr->bitword[l - 1];
    // INCHI✔️❌:     bitWord *FixBits = Fix->bitword[l - 1];
    // INCHI✔️❌:     int     len = Mcr->len_set * sizeof( bitWord );
    // INCHI✔️❌:
    // INCHI✔️❌:     memset( McrBits, 0, len ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( FixBits, 0, len ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     for (i = 0, r = 1; i < n; i++, r++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (r == ( r1 = ( rank_mask_bit&p->Rank[j1 = (int) p->AtNumber[i]] ) ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             FixBits[j1 / pCG->m_num_bit] |= pCG->m_bBit[j1 % pCG->m_num_bit];
    // INCHI✔️❌:             McrBits[j1 / pCG->m_num_bit] |= pCG->m_bBit[j1 % pCG->m_num_bit];
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (r = r1; i + 1 < n && r == ( rank_mask_bit&p->Rank[j2 = (int) p->AtNumber[i + 1]] ); i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (j1 > j2)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     j1 = j2;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             McrBits[j1 / pCG->m_num_bit] |= pCG->m_bBit[j1 % pCG->m_num_bit];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: PartitionGetMcrAndFixSet
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: PartitionGetMcrAndFixSet
    // INCHI✔️❌: typedef U_SHORT AT_RANK
    // INCHI✔️❌: typedef U_SHORT bitWord
    // END INCHI ACTIVE MACRO CONFIGURATION: PartitionGetMcrAndFixSet

    let mcr_bits = pointer_array_get(heap, Mcr.bitword, l.wrapping_sub(1))?;
    let fix_bits = pointer_array_get(heap, Fix.bitword, l.wrapping_sub(1))?;
    let len = usize::try_from(Mcr.len_set).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    heap.slice_mut(mcr_bits)?
        .get_mut(..len)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .fill(0);
    heap.slice_mut(fix_bits)?
        .get_mut(..len)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .fill(0);

    // SAFETY: CanonGraph allocates the atom order, ranks, mask table, Mcr row,
    // and Fix row independently and keeps all five fixed-size buffers live for
    // this call. No helper is invoked while the two writable views are live.
    let atoms = unsafe { heap.stable_slice(p.AtNumber.as_const())? };
    let ranks = unsafe { heap.stable_slice(p.Rank.as_const())? };
    let masks = unsafe { heap.stable_slice(pCG.m_bBit.as_const())? };
    let mut mcr_words = unsafe { heap.stable_slice_mut(mcr_bits)? };
    let mut fix_words = unsafe { heap.stable_slice_mut(fix_bits)? };
    let rank_mask = rank_mask_bit();
    // INCHI✔️✔️:     for (i = 0, r = 1; i < n; i++, r++)
    let mut i = 0_i32;
    let mut r: AT_RANK = 1;
    while i < n {
        let atom_index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        // INCHI✔️✔️:         if (r == ( r1 = ( rank_mask_bit&p->Rank[j1 = (int) p->AtNumber[i]] ) ))
        let mut j1 = i32::from(*atoms.get(atom_index)?);
        let r1 = rank_mask
            & *ranks.get(usize::try_from(j1).map_err(|_| SourceHeapError::PointerOutOfBounds)?)?;
        if r == r1 {
            // INCHI✔️✔️:             FixBits[j1 / pCG->m_num_bit] |= pCG->m_bBit[j1 % pCG->m_num_bit];
            let word_index = j1
                .checked_div(pCG.m_num_bit)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let bit_index = j1
                .checked_rem(pCG.m_num_bit)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mask = *masks.get(
                usize::try_from(bit_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?,
            )?;
            let word_offset =
                usize::try_from(word_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            *fix_words.get_mut(word_offset)? |= mask;
            // INCHI✔️✔️:             McrBits[j1 / pCG->m_num_bit] |= pCG->m_bBit[j1 % pCG->m_num_bit];
            *mcr_words.get_mut(word_offset)? |= mask;
        } else {
            // INCHI✔️✔️:             for (r = r1; i + 1 < n && r == ( rank_mask_bit&p->Rank[j2 = (int) p->AtNumber[i + 1]] ); i++)
            r = r1;
            while i.wrapping_add(1) < n {
                let next_index = usize::try_from(i.wrapping_add(1))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let j2 = i32::from(*atoms.get(next_index)?);
                let next_rank = rank_mask
                    & *ranks.get(
                        usize::try_from(j2).map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )?;
                if r != next_rank {
                    break;
                }
                // INCHI✔️✔️:                 if (j1 > j2)
                if j1 > j2 {
                    // INCHI✔️✔️:                     j1 = j2;
                    j1 = j2;
                }
                i = i.wrapping_add(1);
            }
            // INCHI✔️✔️:             McrBits[j1 / pCG->m_num_bit] |= pCG->m_bBit[j1 % pCG->m_num_bit];
            let word_index = j1
                .checked_div(pCG.m_num_bit)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let bit_index = j1
                .checked_rem(pCG.m_num_bit)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mask = *masks.get(
                usize::try_from(bit_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?,
            )?;
            *mcr_words.get_mut(
                usize::try_from(word_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?,
            )? |= mask;
        }
        i = i.wrapping_add(1);
        r = r.wrapping_add(1);
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn PartitionGetTransposition(
    heap: &mut SourceHeap,
    pFrom: &Partition,
    pTo: &Partition,
    n: i32,
    gamma: &Transposition,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1285 PartitionGetTransposition
    // INCHI✔️❌: void PartitionGetTransposition( Partition *pFrom,
    // INCHI✔️❌:                                 Partition *pTo,
    // INCHI✔️❌:                                 int n,
    // INCHI✔️❌:                                 Transposition *gamma )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         gamma->nAtNumb[(int) pFrom->AtNumber[i]] = pTo->AtNumber[i];
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: PartitionGetTransposition

    if n <= 0 {
        return Ok(());
    }
    // SAFETY: the two partition orders and transposition are independent
    // fixed-size CanonGraph allocations and remain live throughout the loop.
    let from_atoms = unsafe { heap.stable_slice(pFrom.AtNumber.as_const())? };
    let to_atoms = unsafe { heap.stable_slice(pTo.AtNumber.as_const())? };
    let mut transposition = unsafe { heap.stable_slice_mut(gamma.nAtNumb)? };
    // INCHI✔️✔️:     for (i = 0; i < n; i++)
    let mut i = 0_i32;
    while i < n {
        let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let target = usize::from(*from_atoms.get(index)?);
        // INCHI✔️✔️:         gamma->nAtNumb[(int) pFrom->AtNumber[i]] = pTo->AtNumber[i];
        *transposition.get_mut(target)? = *to_atoms.get(index)?;
        i = i.wrapping_add(1);
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn nGetMcr2(
    heap: &mut SourceHeap,
    nEqArray: SourceMutPointer<AT_RANK>,
    n: AT_RANK,
) -> Result<AT_RANK, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1305 nGetMcr2
    // INCHI✔️❌: AT_RANK nGetMcr2( AT_RANK *nEqArray, AT_RANK n )
    // INCHI✔️❌: {
    // INCHI✔️❌:     AT_RANK n1, n2, mcr; /*  recursive version is much shorter. */
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:     n1 = nEqArray[(int) n];
    // INCHI✔️❌:
    // INCHI✔️❌:     if (n == n1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return n;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  1st pass: find mcr */
    // INCHI✔️❌:     while (n1 != ( n2 = nEqArray[(int) n1] ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         n1 = n2;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  2nd pass: copy mcr to each element of the set starting from nEqArray[n] */
    // INCHI✔️❌:     mcr = n1;
    // INCHI✔️❌:     n1 = n;
    // INCHI✔️❌:     while ( /*n1*/ mcr != ( n2 = nEqArray[(int) n1] ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nEqArray[(int) n1] = mcr;
    // INCHI✔️❌:         n1 = n2;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( mcr );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: nGetMcr2

    // SAFETY: the equivalence array is one fixed-size allocation. This
    // function is its sole accessor while path compression is in progress.
    let mut equivalences = unsafe { heap.stable_slice_mut(nEqArray)? };
    // INCHI✔️✔️:     n1 = nEqArray[(int) n];
    let mut n1 = *equivalences.get(usize::from(n))?;
    // INCHI✔️✔️:     if (n == n1)
    if n == n1 {
        return Ok(n);
    }

    // INCHI✔️✔️:     while (n1 != ( n2 = nEqArray[(int) n1] ))
    loop {
        let n2 = *equivalences.get(usize::from(n1))?;
        if n1 == n2 {
            break;
        }
        n1 = n2;
    }

    // INCHI✔️✔️:     mcr = n1;
    let mcr = n1;
    // INCHI✔️✔️:     n1 = n;
    n1 = n;
    // INCHI✔️✔️:     while (mcr != ( n2 = nEqArray[(int) n1] ))
    loop {
        let index = usize::from(n1);
        let n2 = *equivalences.get(index)?;
        if mcr == n2 {
            break;
        }
        // INCHI✔️✔️:         nEqArray[(int) n1] = mcr;
        *equivalences.get_mut(index)? = mcr;
        n1 = n2;
    }
    Ok(mcr)
}

#[allow(non_snake_case)]
pub(crate) fn nJoin2Mcrs2(
    heap: &mut SourceHeap,
    nEqArray: SourceMutPointer<AT_RANK>,
    mut n1: AT_RANK,
    mut n2: AT_RANK,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1342 nJoin2Mcrs2
    // INCHI✔️❌: int nJoin2Mcrs2( AT_RANK *nEqArray, AT_RANK n1, AT_RANK n2 )
    // INCHI✔️❌: {
    // INCHI✔️❌:     n1 = nGetMcr2( nEqArray, n1 );
    // INCHI✔️❌:     n2 = nGetMcr2( nEqArray, n2 );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (n1 < n2)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nEqArray[n2] = n1;
    // INCHI✔️❌:
    // INCHI✔️❌:         INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:         return 1; /*  a change has been made */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (n2 < n1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nEqArray[n1] = n2;
    // INCHI✔️❌:
    // INCHI✔️❌:         INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:         return 1; /*  a change has been made */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0; /*  no changes */
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: nJoin2Mcrs2

    n1 = nGetMcr2(heap, nEqArray, n1)?;
    n2 = nGetMcr2(heap, nEqArray, n2)?;
    if n1 < n2 {
        source_set(heap, nEqArray, i32::from(n2), n1)?;
        return Ok(1);
    }
    if n2 < n1 {
        source_set(heap, nEqArray, i32::from(n1), n2)?;
        return Ok(1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn GetUnorderedPartitionMcrNode(
    heap: &mut SourceHeap,
    p1: &UnorderedPartition,
    v: Node,
) -> Result<Node, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1372 GetUnorderedPartitionMcrNode
    // INCHI✔️❌: Node GetUnorderedPartitionMcrNode( UnorderedPartition *p1, Node v )
    // INCHI✔️❌: {
    // INCHI✔️❌:     Node ret = (Node) ( 1 + nGetMcr2( p1->equ2, (AT_RANK) ( v - 1 ) ) );
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: GetUnorderedPartitionMcrNode

    let zero_based = v.wrapping_sub(1) as AT_RANK;
    let ret = nGetMcr2(heap, p1.equ2, zero_based)?.wrapping_add(1) as Node;
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn UnorderedPartitionJoin(
    heap: &mut SourceHeap,
    p1: &UnorderedPartition,
    p2: &UnorderedPartition,
    n: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1385 UnorderedPartitionJoin
    // INCHI✔️❌: int UnorderedPartitionJoin( UnorderedPartition *p1,
    // INCHI✔️❌:                             UnorderedPartition *p2,
    // INCHI✔️❌:                             int n )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j;
    // INCHI✔️❌:     int nNumChanges = 0;
    // INCHI✔️❌:     for (i = 0; i < n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (( j = (int) p1->equ2[i] ) == i || p2->equ2[(int) i] == p2->equ2[(int) j])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         nNumChanges += nJoin2Mcrs2( p2->equ2, (AT_NUMB) i, (AT_NUMB) j );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:     return nNumChanges;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: UnorderedPartitionJoin

    let mut i = 0_i32;
    let mut nNumChanges = 0_i32;
    while i < n {
        let j = i32::from(source_get(heap, p1.equ2, i)?);
        if j == i || source_get(heap, p2.equ2, i)? == source_get(heap, p2.equ2, j)? {
            i = i.wrapping_add(1);
            continue;
        }
        nNumChanges =
            nNumChanges.wrapping_add(nJoin2Mcrs2(heap, p2.equ2, i as AT_NUMB, j as AT_NUMB)?);
        i = i.wrapping_add(1);
    }
    Ok(nNumChanges)
}

#[allow(non_snake_case)]
pub(crate) fn PartitionSatisfiesLemma_2_25(
    heap: &SourceHeap,
    p: &Partition,
    n: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1407 PartitionSatisfiesLemma_2_25
    // INCHI✔️❌: int PartitionSatisfiesLemma_2_25( Partition *p, int n )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int nPartitionSize = 0;
    // INCHI✔️❌:     int nNumNonTrivialCells = 0;
    // INCHI✔️❌:     AT_RANK r;
    // INCHI✔️❌:     int i, num;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = num = 0, r = 1; i < n; i++, r++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (( rank_mask_bit & p->Rank[(int) p->AtNumber[i]] ) == r)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nPartitionSize++;
    // INCHI✔️❌:             if (num)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* num+1 = cell size > 1 */
    // INCHI✔️❌:                 nNumNonTrivialCells++;
    // INCHI✔️❌:                 num = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* check Lemma_2_25 conditions */
    // INCHI✔️❌:     if (n <= nPartitionSize + 4 ||
    // INCHI✔️❌:          n == nPartitionSize + nNumNonTrivialCells ||
    // INCHI✔️❌:          n == nPartitionSize + nNumNonTrivialCells + 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: PartitionSatisfiesLemma_2_25

    if n <= 0 {
        return Ok(1);
    }
    // SAFETY: both partition arrays are fixed-size, independent allocations
    // owned by CanonGraph and remain immutable for this source loop.
    let atoms = unsafe { heap.stable_slice(p.AtNumber.as_const())? };
    let ranks = unsafe { heap.stable_slice(p.Rank.as_const())? };
    let mut nPartitionSize = 0_i32;
    let mut nNumNonTrivialCells = 0_i32;
    // INCHI✔️✔️:     for (i = num = 0, r = 1; i < n; i++, r++)
    let mut i = 0_i32;
    let mut num = 0_i32;
    let mut r: AT_RANK = 1;
    let rank_mask = rank_mask_bit();
    while i < n {
        let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let atom = usize::from(*atoms.get(index)?);
        // INCHI✔️✔️:         if (( rank_mask_bit & p->Rank[(int) p->AtNumber[i]] ) == r)
        if rank_mask & *ranks.get(atom)? == r {
            // INCHI✔️✔️:             nPartitionSize++;
            nPartitionSize = nPartitionSize.wrapping_add(1);
            // INCHI✔️✔️:             if (num)
            if num != 0 {
                // INCHI✔️✔️:                 nNumNonTrivialCells++;
                nNumNonTrivialCells = nNumNonTrivialCells.wrapping_add(1);
                // INCHI✔️✔️:                 num = 0;
                num = 0;
            }
        } else {
            // INCHI✔️✔️:             num++;
            num = num.wrapping_add(1);
        }
        i = i.wrapping_add(1);
        r = r.wrapping_add(1);
    }

    // INCHI✔️✔️:     if (n <= nPartitionSize + 4 ||
    // INCHI✔️✔️:          n == nPartitionSize + nNumNonTrivialCells ||
    // INCHI✔️✔️:          n == nPartitionSize + nNumNonTrivialCells + 1)
    if n <= nPartitionSize.wrapping_add(4)
        || n == nPartitionSize.wrapping_add(nNumNonTrivialCells)
        || n == nPartitionSize
            .wrapping_add(nNumNonTrivialCells)
            .wrapping_add(1)
    {
        // INCHI✔️✔️:         return 1;
        return Ok(1);
    }
    // INCHI✔️✔️:     return 0;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn PartitionCopy(
    heap: &mut SourceHeap,
    To: &Partition,
    From: &Partition,
    n: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1445 PartitionCopy
    // INCHI✔️✔️: void PartitionCopy( Partition *To, Partition *From, int n )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i;
    // INCHI✔️✔️:     memcpy(To->AtNumber, From->AtNumber, n * sizeof(To->AtNumber[0]));
    // INCHI✔️✔️:     memcpy(To->Rank, From->Rank, n * sizeof(To->AtNumber[0]));
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i = 0; i < n; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         To->Rank[i] &= rank_mask_bit;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     INCHI_HEAPCHK
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: PartitionCopy

    if n <= 0 {
        return Ok(());
    }
    let count = usize::try_from(n).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let source_atoms_are_bounded =
        heap.has_proven_index_bound(From.AtNumber.as_const(), count, count);
    // SAFETY: CanonGraph partitions own distinct arrays that remain live and
    // fixed-size throughout this source call, matching memcpy's C precondition.
    let source_atoms = unsafe { heap.stable_slice(From.AtNumber.as_const())? };
    // INCHI✔️✔️:     memcpy(To->AtNumber, From->AtNumber, n * sizeof(To->AtNumber[0]));
    heap.slice_mut(To.AtNumber)?
        .get_mut(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .copy_from_slice(source_atoms.prefix(count)?);
    drop(source_atoms);
    if source_atoms_are_bounded {
        // The source memcpy copies the complete proved prefix without changing
        // any value, so the destination has exactly the same index bound.
        heap.record_index_bound(To.AtNumber, count, count)?;
    }

    // SAFETY: Rank has the same distinct-allocation and lifetime guarantees.
    let source_ranks = unsafe { heap.stable_slice(From.Rank.as_const())? };
    // INCHI✔️✔️:     memcpy(To->Rank, From->Rank, n * sizeof(To->AtNumber[0]));
    heap.slice_mut(To.Rank)?
        .get_mut(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .copy_from_slice(source_ranks.prefix(count)?);
    drop(source_ranks);

    // INCHI✔️✔️:     for (i = 0; i < n; i++)
    let rank_mask = rank_mask_bit();
    let ranks = heap
        .slice_mut(To.Rank)?
        .get_mut(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut i = 0_usize;
    while i < count {
        // INCHI✔️✔️:         To->Rank[i] &= rank_mask_bit;
        ranks[i] &= rank_mask;
        i += 1;
    }
    Ok(())
}

fn allocate_partition_storage(
    heap: &mut SourceHeap,
    count: usize,
) -> Result<SourceMutPointer<AT_RANK>, SourceHeapError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| SourceHeapError::AllocationFailed)?;
    values.resize(count, 0);
    heap.allocate(values)
}

#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn PartitionColorVertex(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    G: SourceMutPointer<NEIGH_LIST>,
    p: &mut [Partition],
    v: Node,
    n: i32,
    n_tg: i32,
    n_max: i32,
    bDigraph: i32,
    nNumPrevRanks: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1464 PartitionColorVertex
    // INCHI✔️✔️: int PartitionColorVertex( CANON_GLOBALS *pCG,
    // INCHI✔️✔️:                           Graph *G,
    // INCHI✔️✔️:                           Partition *p,
    // INCHI✔️✔️:                           Node v,
    // INCHI✔️✔️:                           int n,
    // INCHI✔️✔️:                           int n_tg,
    // INCHI✔️✔️:                           int n_max,
    // INCHI✔️✔️:                           int bDigraph,
    // INCHI✔️✔️:                           int nNumPrevRanks )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int     nNumNewRanks, i, j;
    // INCHI✔️✔️:     long    lNumNeighListIter = 0;
    // INCHI✔️✔️:     AT_RANK rv, r;
    // INCHI✔️✔️:     AT_NUMB s, sv;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i = 1; i <= 2; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (!p[i].AtNumber)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             p[i].AtNumber = (AT_NUMB *) inchi_malloc( n_max * sizeof( p[0].AtNumber[0] ) );
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (!p[i].Rank)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             p[i].Rank = (AT_RANK *) inchi_malloc( n_max * sizeof( p[0].Rank[0] ) );
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (!p[i].AtNumber || !p[i].Rank)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:
    // INCHI✔️✔️:             INCHI_HEAPCHK
    // INCHI✔️✔️:
    // INCHI✔️✔️:             return CT_OUT_OF_RAM;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     PartitionCopy( p + 1, p, n_tg );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     sv = v - 1;          /* atom number we are looking for */
    // INCHI✔️✔️:     if (sv >= (AT_NUMB) n_tg)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:
    // INCHI✔️✔️:         INCHI_HEAPCHK
    // INCHI✔️✔️:
    // INCHI✔️✔️:         return CT_CANON_ERR; /* !!! severe program error: sv not found !!! */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     rv = p[1].Rank[(int) sv];  /* rank of this atom */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* second, locate sv among all vertices that have same rank as v */
    // INCHI✔️✔️:     s = n_max + 1; /* always greater than sv; this initialization is needed only to keep the compiler happy */
    // INCHI✔️✔️:     for (j = (int) rv - 1; 0 <= j && rv == ( p[1].Rank[(int) ( s = p[1].AtNumber[j] )] ) && s != sv; j--) /* djb-rwth: removing redundant code */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (s != sv)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         INCHI_HEAPCHK
    // INCHI✔️✔️:         return CT_CANON_ERR; /* !!! severe program error: sv not found !!! */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* shift preceding atom numbers to the right to fill the gap after removing sv */
    // INCHI✔️✔️:     r = rv - 1; /* initialization only to keep compiler happy */
    // INCHI✔️✔️:     for (i = j--; 0 <= j && rv == ( r = p[1].Rank[(int) ( s = p[1].AtNumber[j] )] ); i = j, j--)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         p[1].AtNumber[i] = s;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     r = ( i > 0 ) ? ( r + 1 ) : 1;  /* new reduced rank = (next lower rank)+1 or 1 */
    // INCHI✔️✔️:     /* insert sv and adjust its rank */
    // INCHI✔️✔️:     p[1].AtNumber[i] = sv;
    // INCHI✔️✔️:     p[1].Rank[(int) sv] = r;
    // INCHI✔️✔️:
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* make equitable partition */
    // INCHI✔️✔️:     if (bDigraph)
    // INCHI✔️✔️:     {
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         nNumNewRanks = DifferentiateRanks2( pCG, n_tg, G,
    // INCHI✔️❌:                                          nNumPrevRanks+1, p[1].Rank, p[2].Rank,
    // INCHI✔️❌:                                          p[1].AtNumber, &lNumNeighListIter, 1 );
    // INCHI✔️❌:         */
    // INCHI✔️✔️:         nNumNewRanks = DifferentiateRanks4( pCG, n_tg, G,
    // INCHI✔️✔️:                                          nNumPrevRanks + 1, p[1].Rank, p[2].Rank /* temp array */,
    // INCHI✔️✔️:                                          p[1].AtNumber, (AT_RANK) n, &lNumNeighListIter );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️❌:          /*
    // INCHI✔️❌:          nNumNewRanks = DifferentiateRanks2( pCG, n_tg, G,
    // INCHI✔️❌:                                           nNumPrevRanks+1, p[1].Rank, p[2].Rank,
    // INCHI✔️❌:                                           p[1].AtNumber, &lNumNeighListIter, 1 );
    // INCHI✔️❌:          */
    // INCHI✔️✔️:         nNumNewRanks = DifferentiateRanks3( pCG, n_tg, G,
    // INCHI✔️✔️:                                          nNumPrevRanks + 1, p[1].Rank, p[2].Rank /* temp array */,
    // INCHI✔️✔️:                                          p[1].AtNumber, &lNumNeighListIter );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     INCHI_HEAPCHK
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return nNumNewRanks;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: PartitionColorVertex

    if p.len() < 3 {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let allocation_count = usize::try_from(n_max).ok();
    for index in 1..=2 {
        if p[index].AtNumber.is_null() {
            let allocation = allocation_count
                .ok_or(SourceHeapError::AllocationFailed)
                .and_then(|count| allocate_partition_storage(heap, count));
            match allocation {
                Ok(pointer) => p[index].AtNumber = pointer,
                Err(SourceHeapError::AllocationFailed) => {}
                Err(error) => return Err(error),
            }
        }
        if p[index].Rank.is_null() {
            let allocation = allocation_count
                .ok_or(SourceHeapError::AllocationFailed)
                .and_then(|count| allocate_partition_storage(heap, count));
            match allocation {
                Ok(pointer) => p[index].Rank = pointer,
                Err(SourceHeapError::AllocationFailed) => {}
                Err(error) => return Err(error),
            }
        }
        if p[index].AtNumber.is_null() || p[index].Rank.is_null() {
            return Ok(CT_OUT_OF_RAM);
        }
    }

    let (source_storage, remaining) = p.split_at_mut(1);
    let source = source_storage
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let (destination_storage, temporary_storage) = remaining.split_at_mut(1);
    let destination = destination_storage
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let temporary = temporary_storage
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    PartitionCopy(heap, destination, source, n_tg)?;

    let sv = v.wrapping_sub(1);
    if sv >= n_tg as AT_NUMB {
        return Ok(CT_CANON_ERR);
    }

    // SAFETY: PartitionCreate gives AtNumber and Rank distinct allocations of
    // `n_max` elements, and this function neither frees nor resizes them. When
    // PartitionCopy propagated the source permutation bound, every write below
    // only moves an existing member inside that same prefix.
    let atom_count = usize::try_from(n_tg).ok();
    let atoms_are_bounded = atom_count.is_some_and(|count| {
        heap.has_proven_index_bound(destination.AtNumber.as_const(), count, count)
    });
    let mut atoms = if atoms_are_bounded {
        let count = atom_count.expect("a proved non-negative prefix has a valid count");
        unsafe { heap.stable_index_bounded_slice_mut(destination.AtNumber, count, count)? }
    } else {
        unsafe { heap.stable_slice_mut(destination.AtNumber)? }
    };
    // SAFETY: see above; Rank is the other independently owned partition array.
    let mut ranks = unsafe { heap.stable_slice_mut(destination.Rank)? };

    // INCHI✔️✔️:     rv = p[1].Rank[sv];
    let rv = *ranks.get(usize::from(sv))?;
    // INCHI✔️✔️:     for (s = n_max + 1, j = rv - 1;
    let mut s = n_max.wrapping_add(1) as AT_NUMB;
    let mut j = i32::from(rv).wrapping_sub(1);
    // INCHI✔️✔️:          0 <= j && rv == p[1].Rank[(int)(s = p[1].AtNumber[j])] && s != sv;
    // INCHI✔️✔️:          j--)
    while j >= 0 {
        s = *atoms.get(usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?)?;
        if rv != *ranks.get(usize::from(s))? || s == sv {
            break;
        }
        j = j.wrapping_sub(1);
    }
    // INCHI✔️✔️:     if (s != sv)
    if s != sv {
        // INCHI✔️✔️:         return CT_CANON_ERR;
        return Ok(CT_CANON_ERR);
    }

    // INCHI✔️✔️:     r = rv - 1;
    let mut r = rv.wrapping_sub(1);
    // INCHI✔️✔️:     for (i = j--; 0 <= j && rv == (r = p[1].Rank[(int)(s = p[1].AtNumber[j])]); i = j, j--)
    let mut i = j;
    j = j.wrapping_sub(1);
    while j >= 0 {
        s = *atoms.get(usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?)?;
        r = *ranks.get(usize::from(s))?;
        if rv != r {
            break;
        }
        // INCHI✔️✔️:         p[1].AtNumber[i] = s;
        *atoms.get_mut(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)? = s;
        i = j;
        j = j.wrapping_sub(1);
    }
    // INCHI✔️✔️:     r = (i > 0) ? (r + 1) : 1;
    r = if i > 0 { r.wrapping_add(1) } else { 1 };
    // INCHI✔️✔️:     p[1].AtNumber[i] = sv;
    *atoms.get_mut(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)? = sv;
    // INCHI✔️✔️:     p[1].Rank[(int)sv] = r;
    *ranks.get_mut(usize::from(sv))? = r;
    drop(atoms);
    drop(ranks);

    let mut lNumNeighListIter = 0_i64;
    if bDigraph != 0 {
        crate::source::base::ichimap2::DifferentiateRanks4(
            heap,
            pCG,
            n_tg,
            G,
            nNumPrevRanks.wrapping_add(1),
            destination.Rank,
            temporary.Rank,
            destination.AtNumber,
            n as AT_RANK,
            &mut lNumNeighListIter,
        )
    } else {
        crate::source::base::ichimap2::DifferentiateRanks3(
            heap,
            pCG,
            n_tg,
            G,
            nNumPrevRanks.wrapping_add(1),
            destination.Rank,
            temporary.Rank,
            destination.AtNumber,
            &mut lNumNeighListIter,
        )
    }
}

#[allow(non_snake_case)]
pub(crate) fn CellGetMinNode(
    heap: &SourceHeap,
    p: &Partition,
    W: &Cell,
    v: Node,
    pCD: Option<&CANON_DATA>,
) -> Result<Node, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1586 CellGetMinNode
    // INCHI✔️❌: Node CellGetMinNode( Partition *p, Cell *W, Node v, CANON_DATA *pCD )
    // INCHI✔️❌: {
    // INCHI✔️❌:     AT_NUMB i;
    // INCHI✔️❌:     AT_NUMB uCurAtNumb, uMinAtNumb = INCHI_CANON_INFINITY;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* in case of emty cell:  (W->first=INCHI_CANON_INFINITY) > (W->next=0); returns INCHI_CANON_INFINITY */
    // INCHI✔️❌:     if (W->first > W->next)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return INCHI_CANON_INFINITY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_AUX_RANKING == 1 )
    // INCHI✔️❌:     if (pCD && pCD->nAuxRank)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         AT_RANK uMinAuxRank, uCurAuxRank;
    // INCHI✔️❌:         int     nCurAtNumb;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_AUX_RANKING_ALL == 1 )
    // INCHI✔️❌:         AT_RANK uInpAuxRank;
    // INCHI✔️❌:         int     nInpAtNumb, nMinAtNumb;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = W->first; i < W->next; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             uCurAtNumb = p->AtNumber[(int) i];
    // INCHI✔️❌:             if (!( p->Rank[(int) uCurAtNumb] & rank_mark_bit ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 break; /* found the first unmarked yet node */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i == W->next)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return INCHI_CANON_INFINITY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_AUX_RANKING_ALL == 1 )
    // INCHI✔️❌:         /*==== vertex ordering definition ===
    // INCHI✔️❌:          * vertex v1 < v2 <=> (AuxRank(v1)==AuxRank(v2) && AtNumb(v1) < AtNumb(v2)) || (AuxRank(v1) < AuxRank(v2))
    // INCHI✔️❌:          * vertex v1 > v2 <=> (AuxRank(v1)==AuxRank(v2) && AtNumb(v1) > AtNumb(v2)) || (AuxRank(v1) > AuxRank(v2))
    // INCHI✔️❌:          * vertex v1 = v2 <=> (AuxRank(v1)==AuxRank(v2) && AtNumb(v1) == AtNumb(v2))
    // INCHI✔️❌:          */
    // INCHI✔️❌:
    // INCHI✔️❌:         /* set initial vMin so that vMin > any vertex */
    // INCHI✔️❌:         uMinAuxRank = INCHI_CANON_INFINITY;
    // INCHI✔️❌:         nMinAtNumb = INCHI_CANON_INFINITY;
    // INCHI✔️❌:         /* set vInp */
    // INCHI✔️❌:         if (v)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nInpAtNumb = (int) v - 1; /* previous vertex */
    // INCHI✔️❌:             uInpAuxRank = pCD->nAuxRank[nInpAtNumb];
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nInpAtNumb = -1; /* less than any vertex */
    // INCHI✔️❌:             uInpAuxRank = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* find vMin = min { vCur : (vCur > vInp) && (vCur in W) } */
    // INCHI✔️❌:         for (; i < W->next; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nCurAtNumb = (int) p->AtNumber[(int) i];
    // INCHI✔️❌:             if (!( p->Rank[nCurAtNumb] & rank_mark_bit ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* vertex nCurAtNumb is not marked, find whether it fits the conditions */
    // INCHI✔️❌:                 uCurAuxRank = pCD->nAuxRank[nCurAtNumb];
    // INCHI✔️❌:                 if ((uCurAuxRank == uInpAuxRank && nCurAtNumb > nInpAtNumb) || uCurAuxRank > uInpAuxRank) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* here vCur > vInp */
    // INCHI✔️❌:                     if (uCurAuxRank == uMinAuxRank && nCurAtNumb < nMinAtNumb)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* vCur < vMin (1) */
    // INCHI✔️❌:                         nMinAtNumb = nCurAtNumb;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else if (uCurAuxRank < uMinAuxRank)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* vCur < vMin (2) */
    // INCHI✔️❌:                         uMinAuxRank = uCurAuxRank;
    // INCHI✔️❌:                         nMinAtNumb = nCurAtNumb;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         uMinAtNumb = ( nMinAtNumb == INCHI_CANON_INFINITY ) ? INCHI_CANON_INFINITY : (AT_NUMB) nMinAtNumb;
    // INCHI❌❌: #else
    // INCHI❌❌:
    // INCHI❌❌:         if (v)
    // INCHI❌❌:         {
    // INCHI❌❌:             nCurAtNumb = (int) v - 1;
    // INCHI❌❌:             /* any valid found node must have nAuxRank == uMinAuxRank */
    // INCHI❌❌:             uMinAuxRank = pCD->nAuxRank[nCurAtNumb];
    // INCHI❌❌:         }
    // INCHI❌❌:         else
    // INCHI❌❌:         {
    // INCHI❌❌:             /* any valid found node must have minimal uMinAuxRank from pCD->nAuxRank[] */
    // INCHI❌❌:             uMinAuxRank = INCHI_CANON_INFINITY; /* undefined */
    // INCHI❌❌:         }
    // INCHI❌❌:
    // INCHI❌❌:         for (; i < W->next; i++)
    // INCHI❌❌:         {
    // INCHI❌❌:             uCurAtNumb = p->AtNumber[(int) i];
    // INCHI❌❌:             nCurAtNumb = (int) uCurAtNumb;
    // INCHI❌❌:             if (uCurAtNumb >= v && !( p->Rank[nCurAtNumb] & rank_mark_bit ))
    // INCHI❌❌:             {
    // INCHI❌❌:                 uCurAuxRank = pCD->nAuxRank[nCurAtNumb];
    // INCHI❌❌:                 if (v)
    // INCHI❌❌:                 {
    // INCHI❌❌:                     /* get next node */
    // INCHI❌❌:                     /* find node with smallest uCurAtNumb among nodes with aux. ranks equal to uMinAuxRank */
    // INCHI❌❌:                     if (uCurAuxRank == uMinAuxRank && uCurAtNumb < uMinAtNumb)
    // INCHI❌❌:                     {
    // INCHI❌❌:                         uMinAtNumb = uCurAtNumb;
    // INCHI❌❌:                     }
    // INCHI❌❌:                 }
    // INCHI❌❌:                 else
    // INCHI❌❌:                 {
    // INCHI❌❌:                     /* get first node */
    // INCHI❌❌:                     /* find node with smallest smallest uCurAtNumb among nodes with smallest aux. ranks */
    // INCHI❌❌:                     if (uMinAuxRank > uCurAuxRank)
    // INCHI❌❌:                     {
    // INCHI❌❌:                         uMinAuxRank = uCurAuxRank;
    // INCHI❌❌:                         uMinAtNumb = uCurAtNumb;
    // INCHI❌❌:                     }
    // INCHI❌❌:                     else
    // INCHI❌❌:                     {
    // INCHI❌❌:                         if (uMinAuxRank == uCurAuxRank && uCurAtNumb < uMinAtNumb)
    // INCHI❌❌:                         {
    // INCHI❌❌:                             uMinAtNumb = uCurAtNumb;
    // INCHI❌❌:                         }
    // INCHI❌❌:                     }
    // INCHI❌❌:                 }
    // INCHI❌❌:             }
    // INCHI❌❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     else
    // INCHI✔️❌: #endif /* } USE_AUX_RANKING */
    // INCHI✔️❌:
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = W->first; i < W->next; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             uCurAtNumb = p->AtNumber[(int) i];
    // INCHI✔️❌:             if (uCurAtNumb >= v && !( p->Rank[(int) uCurAtNumb] & rank_mark_bit ) && uCurAtNumb < uMinAtNumb)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 uMinAtNumb = uCurAtNumb;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (uMinAtNumb != INCHI_CANON_INFINITY)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         uMinAtNumb++;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:     return uMinAtNumb;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CellGetMinNode

    let infinity = INCHI_CANON_INFINITY as AT_NUMB;
    if W.first > W.next {
        return Ok(infinity);
    }

    let mut minimum_atom = infinity;
    let mut i = W.first as AT_NUMB;
    if i32::from(i) >= W.next {
        return Ok(infinity);
    }
    // SAFETY: the partition order and ranks are independent fixed-size
    // allocations kept live and immutable throughout this cell scan.
    let atoms = unsafe { heap.stable_slice(p.AtNumber.as_const())? };
    let ranks = unsafe { heap.stable_slice(p.Rank.as_const())? };
    if let Some(canon_data) = pCD.filter(|canon_data| !canon_data.nAuxRank.is_null()) {
        // SAFETY: nAuxRank is another immutable CanonGraph allocation with the
        // same vertex-domain lifetime as the partition arrays.
        let aux_ranks = unsafe { heap.stable_slice(canon_data.nAuxRank.as_const())? };
        // INCHI✔️✔️:         for (i = W->first; i < W->next; i++)
        while i32::from(i) < W.next {
            let atom = *atoms.get(usize::from(i))?;
            // INCHI✔️✔️:             if (!( p->Rank[(int) uCurAtNumb] & rank_mark_bit ))
            if *ranks.get(usize::from(atom))? & rank_mark_bit() == 0 {
                break;
            }
            i = i.wrapping_add(1);
        }
        // INCHI✔️✔️:         if (i == W->next)
        if i32::from(i) == W.next {
            // INCHI✔️✔️:             return INCHI_CANON_INFINITY;
            return Ok(infinity);
        }

        let mut minimum_aux_rank = infinity as AT_RANK;
        let mut minimum_atom_number = i32::from(infinity);
        let (input_atom_number, input_aux_rank) = if v != 0 {
            let atom_number = i32::from(v).wrapping_sub(1);
            (
                atom_number,
                *aux_ranks.get(
                    usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )?,
            )
        } else {
            (-1, 0)
        };
        // INCHI✔️✔️:         for (; i < W->next; i++)
        while i32::from(i) < W.next {
            let current_atom_number = i32::from(*atoms.get(usize::from(i))?);
            // INCHI✔️✔️:             if (!( p->Rank[nCurAtNumb] & rank_mark_bit ))
            if *ranks.get(
                usize::try_from(current_atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
            )? & rank_mark_bit()
                == 0
            {
                let current_aux_rank = *aux_ranks.get(
                    usize::try_from(current_atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )?;
                if (current_aux_rank == input_aux_rank && current_atom_number > input_atom_number)
                    || current_aux_rank > input_aux_rank
                {
                    if current_aux_rank == minimum_aux_rank
                        && current_atom_number < minimum_atom_number
                    {
                        minimum_atom_number = current_atom_number;
                    } else if current_aux_rank < minimum_aux_rank {
                        minimum_aux_rank = current_aux_rank;
                        minimum_atom_number = current_atom_number;
                    }
                }
            }
            i = i.wrapping_add(1);
        }
        minimum_atom = if minimum_atom_number == i32::from(infinity) {
            infinity
        } else {
            minimum_atom_number as AT_NUMB
        };
    } else {
        // INCHI✔️✔️:         for (i = W->first; i < W->next; i++)
        while i32::from(i) < W.next {
            let current_atom = *atoms.get(usize::from(i))?;
            // INCHI✔️✔️:             if (uCurAtNumb >= v && !( p->Rank[(int) uCurAtNumb] & rank_mark_bit ) && uCurAtNumb < uMinAtNumb)
            if current_atom >= v
                && *ranks.get(usize::from(current_atom))? & rank_mark_bit() == 0
                && current_atom < minimum_atom
            {
                minimum_atom = current_atom;
            }
            i = i.wrapping_add(1);
        }
    }
    // INCHI✔️✔️:     if (uMinAtNumb != INCHI_CANON_INFINITY)
    if minimum_atom != infinity {
        // INCHI✔️✔️:         uMinAtNumb++;
        minimum_atom = minimum_atom.wrapping_add(1);
    }
    Ok(minimum_atom)
}

#[allow(non_snake_case)]
pub(crate) fn CellGetNumberOfNodes(
    heap: &SourceHeap,
    p: &Partition,
    W: &Cell,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1748 CellGetNumberOfNodes
    // INCHI✔️❌: int CellGetNumberOfNodes( Partition *p, Cell *W )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int first = W->first;
    // INCHI✔️❌:     int next = W->next;
    // INCHI✔️❌:     int i, num;
    // INCHI✔️❌:     for (i = first, num = 0; i < next; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (!( rank_mark_bit & p->Rank[(int) p->AtNumber[i]] ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:     return num;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CellGetNumberOfNodes

    if W.first >= W.next {
        return Ok(0);
    }
    // SAFETY: both partition arrays remain live, immutable, and independently
    // allocated throughout the cell scan.
    let atoms = unsafe { heap.stable_slice(p.AtNumber.as_const())? };
    let ranks = unsafe { heap.stable_slice(p.Rank.as_const())? };
    // INCHI✔️✔️:     for (i = first, num = 0; i < next; i++)
    let mut i = W.first;
    let mut number = 0_i32;
    while i < W.next {
        let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let atom = usize::from(*atoms.get(index)?);
        // INCHI✔️✔️:         if (!( rank_mark_bit & p->Rank[(int) p->AtNumber[i]] ))
        if rank_mark_bit() & *ranks.get(atom)? == 0 {
            // INCHI✔️✔️:             num++;
            number = number.wrapping_add(1);
        }
        i = i.wrapping_add(1);
    }
    Ok(number)
}

#[allow(non_snake_case)]
pub(crate) fn CellIntersectWithSet(
    heap: &mut SourceHeap,
    pCG: &CANON_GLOBALS,
    p: &Partition,
    W: &Cell,
    Mcr: &NodeSet,
    l: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1768 CellIntersectWithSet
    // INCHI✔️❌: int CellIntersectWithSet( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                           Partition *p, Cell *W, NodeSet *Mcr, int l )
    // INCHI✔️❌: {
    // INCHI✔️❌:     bitWord *McrBits = Mcr->bitword[l - 1];
    // INCHI✔️❌:     int first = W->first;
    // INCHI✔️❌:     int next = W->next;
    // INCHI✔️❌:     int i, j, k;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (first >= next)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* for testing only */
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (i = first, k = 0; i < next; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         j = (int) p->AtNumber[i];
    // INCHI✔️❌:         if (!( McrBits[j / pCG->m_num_bit] & pCG->m_bBit[j % pCG->m_num_bit] ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* BC: reading uninit memory ???-not examined yet */
    // INCHI✔️❌:             k += !( p->Rank[j] & rank_mark_bit ); /* for testing only */
    // INCHI✔️❌:             p->Rank[j] |= rank_mark_bit;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:     return k;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CellIntersectWithSet

    let mcr_bits = pointer_array_get(heap, Mcr.bitword, l.wrapping_sub(1))?;
    let first = W.first;
    let next = W.next;
    if first >= next {
        return Ok(0);
    }

    // SAFETY: CanonGraph owns the atom order, rank array, and mask table as
    // independent fixed-size allocations. NodeSet owns its row in another
    // fixed-size allocation. All four remain live without resize or free for
    // this call, and only the rank allocation is written by the source loop.
    let atoms = unsafe { heap.stable_slice(p.AtNumber.as_const())? };
    let mcr_words = unsafe { heap.stable_slice(mcr_bits.as_const())? };
    let masks = unsafe { heap.stable_slice(pCG.m_bBit.as_const())? };
    let mut ranks = unsafe { heap.stable_slice_mut(p.Rank)? };
    // INCHI✔️✔️:     for (i = first, k = 0; i < next; i++)
    let mut i = first;
    let mut k = 0_i32;
    while i < next {
        let atom_index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        // INCHI✔️✔️:         j = (int) p->AtNumber[i];
        let j = i32::from(*atoms.get(atom_index)?);
        let word_index = j
            .checked_div(pCG.m_num_bit)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let bit_index = j
            .checked_rem(pCG.m_num_bit)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        // INCHI✔️✔️:         if (!( McrBits[j / pCG->m_num_bit] &
        // INCHI✔️✔️:                pCG->m_bBit[j % pCG->m_num_bit] ))
        let word = *mcr_words
            .get(usize::try_from(word_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)?;
        let mask = *masks
            .get(usize::try_from(bit_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)?;
        if word & mask == 0 {
            let rank_index = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            // INCHI✔️✔️:             k += !( p->Rank[j] & rank_mark_bit );
            let rank = *ranks.get(rank_index)?;
            k = k.wrapping_add(i32::from(rank & rank_mark_bit() == 0));
            // INCHI✔️✔️:             p->Rank[j] |= rank_mark_bit;
            *ranks.get_mut(rank_index)? = rank | rank_mark_bit();
        }
        i = i.wrapping_add(1);
    }
    Ok(k)
}

#[allow(non_snake_case)]
pub(crate) fn CtPartClear(
    heap: &mut SourceHeap,
    Ct: &mut ConTable,
    k: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1798 CtPartClear
    // INCHI✔️✔️: void CtPartClear( ConTable *Ct, int k )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int start;
    // INCHI✔️✔️:     int len;
    // INCHI✔️✔️:     /* connection table */
    // INCHI✔️✔️:     start = k > 1 ? Ct->nextCtblPos[k - 1] : 0;
    // INCHI✔️✔️:     len = Ct->lenCt - start;
    // INCHI✔️✔️:     if (len > 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         memset( Ct->Ctbl + start, 0, ( (long long)Ct->lenCt - (long long)start ) * sizeof( Ct->Ctbl[0] ) ); /* djb-rwth: cast operators added; memset_s C11/Annex K variant? */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     Ct->lenCt = start;
    // INCHI✔️✔️:     Ct->lenPos = k;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     INCHI_HEAPCHK
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: CtPartClear

    let start = if k > 1 {
        i32::from(source_get(heap, Ct.nextCtblPos, k.wrapping_sub(1))?)
    } else {
        0
    };
    let len = Ct.lenCt.wrapping_sub(start);
    if len > 0 {
        let start_index =
            usize::try_from(start).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let end_index =
            usize::try_from(Ct.lenCt).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        heap.slice_mut(Ct.Ctbl)?
            .get_mut(start_index..end_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .fill(0);
    }
    Ct.lenCt = start;
    Ct.lenPos = k;
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn insertions_sort_NeighList_AT_NUMBERS2(
    heap: &mut SourceHeap,
    base: NEIGH_LIST,
    nRank: SourceMutPointer<AT_RANK>,
    max_rj: AT_RANK,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1819 insertions_sort_NeighList_AT_NUMBERS2
    // INCHI✔️✔️: void insertions_sort_NeighList_AT_NUMBERS2( NEIGH_LIST base,
    // INCHI✔️✔️:                                             AT_RANK *nRank,
    // INCHI✔️✔️:                                             AT_RANK max_rj )
    // INCHI✔️✔️: {

    // SAFETY: every graph row is a fixed-size allocation distinct from the
    // partition rank array and both remain live throughout canonicalization.
    let mut neighbors = unsafe { heap.stable_slice_mut(base)? };
    // INCHI✔️✔️:     int k, num = (int) *base++;
    let num = i32::from(*neighbors.get(0)?);
    if num <= 1 {
        return Ok(());
    }
    let ranks = unsafe { heap.stable_slice(nRank.as_const())? };
    // INCHI✔️✔️:     for (k = 1, pk = base; k < num; k++, pk++)
    let mut k = 1_i32;
    while k < num {
        // INCHI✔️✔️:         i = pk;
        let mut i = k;
        // INCHI✔️✔️:         j = i + 1;
        let mut j = i.wrapping_add(1);
        let j_index = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let j_atom = *neighbors.get(j_index)?;
        // INCHI✔️✔️:         rj = ( rank_mask_bit & nRank[(int) *j] );
        let rj = rank_mask_bit() & *ranks.get(usize::from(j_atom))?;
        // INCHI✔️✔️:         if (rj < max_rj)
        if rj < max_rj {
            // INCHI✔️✔️:             while (j > base && rj < ( rank_mask_bit & nRank[(int) *i] ))
            while j > 1 {
                let i_index =
                    usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let i_atom = *neighbors.get(i_index)?;
                if rj >= rank_mask_bit() & *ranks.get(usize::from(i_atom))? {
                    break;
                }
                // INCHI✔️✔️:                 tmp = *i;
                // INCHI✔️✔️:                 *i = *j;
                *neighbors.get_mut(i_index)? = j_atom;
                // INCHI✔️✔️:                 *j = tmp;
                *neighbors.get_mut(
                    usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )? = i_atom;
                // INCHI✔️✔️:                 j = i--;
                j = i;
                i = i.wrapping_sub(1);
            }
        }
        k = k.wrapping_add(1);
    }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: insertions_sort_NeighList_AT_NUMBERS2
    Ok(())
}

/// Sorts one source-created graph row after proving the complete row and rank
/// index domain once, and returns the same writable row for `CtPartFill`.
fn try_sort_neighbor_row_source_layout(
    heap: &mut SourceHeap,
    base: NEIGH_LIST,
    ranks: &[AT_RANK],
    max_rj: AT_RANK,
    rank_mask: AT_RANK,
) -> Result<Option<StableSourceSlice<AT_NUMB>>, SourceHeapError> {
    let mut neighbors = unsafe { heap.stable_slice_mut(base)? };
    let Some(&num) = neighbors.get(0).ok() else {
        return Ok(None);
    };
    let row_len = usize::from(num)
        .checked_add(1)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let Ok(row) = neighbors.prefix_mut(row_len) else {
        return Ok(None);
    };
    if row[1..]
        .iter()
        .any(|&atom| usize::from(atom) >= ranks.len())
    {
        return Ok(None);
    }

    // INCHI✔️✔️:     for (k = 1, pk = base; k < num; k++, pk++)
    let mut k = 1_usize;
    while k < usize::from(num) {
        // INCHI✔️✔️:         i = pk;
        let mut i = k;
        // INCHI✔️✔️:         j = i + 1;
        let mut j = i + 1;
        // SAFETY: `row_len` proves every row index used by the source loop,
        // and the domain scan above proves every rank index before mutation.
        let j_atom = unsafe { *row.get_unchecked(j) };
        // INCHI✔️✔️:         rj = ( rank_mask_bit & nRank[(int) *j] );
        let rj = rank_mask & unsafe { *ranks.get_unchecked(usize::from(j_atom)) };
        // INCHI✔️✔️:         if (rj < max_rj)
        if rj < max_rj {
            // INCHI✔️✔️:             while (j > base && rj < ( rank_mask_bit & nRank[(int) *i] ))
            while j > 1 {
                let i_atom = unsafe { *row.get_unchecked(i) };
                if rj >= rank_mask & unsafe { *ranks.get_unchecked(usize::from(i_atom)) } {
                    break;
                }
                // INCHI✔️✔️:                 tmp = *i;
                // INCHI✔️✔️:                 *i = *j;
                unsafe { *row.get_unchecked_mut(i) = j_atom };
                // INCHI✔️✔️:                 *j = tmp;
                unsafe { *row.get_unchecked_mut(j) = i_atom };
                // INCHI✔️✔️:                 j = i--;
                j = i;
                i -= 1;
            }
        }
        k += 1;
    }
    Ok(Some(neighbors))
}

#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn CtPartFill(
    heap: &mut SourceHeap,
    G: SourceMutPointer<NEIGH_LIST>,
    pCD: &CANON_DATA,
    p: &Partition,
    Ct: Option<&mut ConTable>,
    mut k: i32,
    n: i32,
    n_tg: i32,
    n_max: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1849 CtPartFill
    // INCHI✔️❌: void CtPartFill( Graph *G,
    // INCHI✔️❌:                  CANON_DATA *pCD,
    // INCHI✔️❌:                  Partition *p,
    // INCHI✔️❌:                  ConTable *Ct,
    // INCHI✔️❌:                  int k,
    // INCHI✔️❌:                  int n,
    // INCHI✔️❌:                  int n_tg,
    // INCHI✔️❌:                  int n_max
    // INCHI✔️❌:                 )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*  k = (new index in Ct->nextAtRank[] and Ct->nextCtblPos[]) + 1 */
    // INCHI✔️❌:
    // INCHI✔️❌:     int     startCtbl;
    // INCHI✔️❌:     int     startAtOrd;
    // INCHI✔️❌:     AT_RANK r, rj, nn, j, rj_prev; /* djb-rwth: ignoring LLVM warning as the variable is used */
    // INCHI✔️❌:     int     i, m, an_sao;
    // INCHI✔️❌:
    // INCHI❌❌: #ifdef INCHI_CANON_USE_HASH
    // INCHI❌❌:     CtHash  hash = 0;
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:     k--;
    // INCHI✔️❌:     if (Ct && k) /* djb-rwth: fixing oss-fuzz issue #69612 */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         startCtbl = Ct->nextCtblPos[k - 1];
    // INCHI✔️❌:         startAtOrd = Ct->nextAtRank[k - 1] - 1;  /* here  p->Rank[p->AtNumber[r-1]] = r */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         startCtbl = 0;
    // INCHI✔️❌:         startAtOrd = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /******* Well-defined (by fixed ranks) part of the connection table ************/
    // INCHI✔️❌:     /* djb-rwth: fixing oss-fuzz issue #69612 */
    // INCHI✔️❌:     if (startAtOrd < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* djb-rwth: fixing oss-fuzz issue #391043585 */
    // INCHI✔️❌:     if (startAtOrd < n_max && Ct)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         an_sao = (int)p->AtNumber[startAtOrd];
    // INCHI✔️❌:         r = (rank_mask_bit & p->Rank[an_sao]);
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = startAtOrd; i < n_tg && r == (rank_mask_bit & p->Rank[m = (int)p->AtNumber[i]]); i++, r++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             Ct->Ctbl[startCtbl++] = r;
    // INCHI✔️❌:             insertions_sort_NeighList_AT_NUMBERS2(G[m], p->Rank, r);
    // INCHI✔️❌:             nn = G[m][0];   /* number of neighbors */
    // INCHI✔️❌:             rj_prev = 0;    /* debug only */ /* djb-rwth: ignoring LLVM warning as the variable is used */
    // INCHI✔️❌:
    // INCHI❌❌: #ifdef INCHI_CANON_USE_HASH
    // INCHI❌❌:             hash = add2crc32(hash, (AT_NUMB)(r + n));
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:             for (j = 1; j <= nn && (rj = (rank_mask_bit & p->Rank[(int)G[m][j]])) < r; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 Ct->Ctbl[startCtbl++] = rj;
    // INCHI✔️❌:
    // INCHI❌❌: #ifdef INCHI_CANON_USE_HASH
    // INCHI❌❌:                 hash = add2crc32(hash, rj);
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI❌❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI❌❌:                 /* debug only */
    // INCHI❌❌:                 if (rj < rj_prev)
    // INCHI❌❌:                 {
    // INCHI❌❌:                     int stop = 1;   /* <BRKPT> */
    // INCHI❌❌:                 }
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:                 rj_prev = rj; /* djb-rwth: ignoring LLVM warning as the variable is used */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:     /****************** Well-defined part of base hydrogen atoms *******************/
    // INCHI✔️❌:     if (Ct)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (pCD->NumH && Ct->NumH)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nn = inchi_min(n, i);
    // INCHI✔️❌:             for (j = startAtOrd; j < nn; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* atoms */
    // INCHI✔️❌:                 Ct->NumH[j] = pCD->NumH[p->AtNumber[j]];
    // INCHI✔️❌:             }
    // INCHI✔️❌:             for (; j < i; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* t-groups */
    // INCHI✔️❌:                 int data_pos = n + T_NUM_NO_ISOTOPIC * ((int)p->AtNumber[j] - n);
    // INCHI✔️❌:                 for (m = 0; m < T_NUM_NO_ISOTOPIC; m++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     Ct->NumH[nn++] = pCD->NumH[data_pos++];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             Ct->lenNumH = nn;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             Ct->lenNumH = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:             /****************** Well-defined part of fixed hydrogen atoms *******************/
    // INCHI✔️❌:             if (pCD->NumHfixed && Ct->NumHfixed)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nn = inchi_min(n, i);
    // INCHI✔️❌:                 for (j = startAtOrd; j < nn; j++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     Ct->NumHfixed[j] = pCD->NumHfixed[p->AtNumber[j]];
    // INCHI✔️❌:
    // INCHI✔️❌:                     INCHI_HEAPCHK
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* Ct->lenNumHfixed = nn; */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;/* Ct->lenNumHfixed = 0; */
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:         INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:             /****************** Well-defined part of isotopic keys ***************************/
    // INCHI✔️❌:             if (pCD->iso_sort_key && Ct->iso_sort_key)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (j = startAtOrd; j < i; j++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     Ct->iso_sort_key[j] = pCD->iso_sort_key[p->AtNumber[j]];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 Ct->len_iso_sort_key = i;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 Ct->len_iso_sort_key = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:         INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:             /****************** Well-defined part of isotopic iso_exchg_atnos ***************************/
    // INCHI✔️❌:             if (pCD->iso_exchg_atnos && Ct->iso_exchg_atnos)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (j = startAtOrd; j < i; j++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     Ct->iso_exchg_atnos[j] = pCD->iso_exchg_atnos[p->AtNumber[j]];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 Ct->len_iso_exchg_atnos = i;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 Ct->len_iso_exchg_atnos = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:         INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI❌❌:             /******** Well-defined part of isotopic keys for fixed hydrogen atoms ************/
    // INCHI❌❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI❌❌:             if (pCD->iso_sort_key_Hfixed && Ct->iso_sort_key_Hfixed)
    // INCHI❌❌:             {
    // INCHI❌❌:                 nn = inchi_min(n, i);
    // INCHI❌❌:                 for (j = startAtOrd; j < nn; j++)
    // INCHI❌❌:                 {
    // INCHI❌❌:                     Ct->iso_sort_key_Hfixed[j] = pCD->iso_sort_key_Hfixed[p->AtNumber[j]];
    // INCHI❌❌:                 }
    // INCHI❌❌:                 Ct->len_iso_sort_key_Hfixed = nn;
    // INCHI❌❌:             }
    // INCHI❌❌:             else
    // INCHI❌❌:             {
    // INCHI❌❌:                 Ct->len_iso_sort_key_Hfixed = 0;
    // INCHI❌❌:             }
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:             Ct->lenCt = startCtbl; /* not aways increases */
    // INCHI✔️❌:         Ct->nextCtblPos[k] = startCtbl;
    // INCHI✔️❌:         Ct->nextAtRank[k] = r;
    // INCHI✔️❌:         Ct->lenPos = k + 1;
    // INCHI✔️❌:
    // INCHI✔️❌:         /* The rest of the CTable */
    // INCHI✔️❌:
    // INCHI❌❌: #ifdef INCHI_CANON_USE_HASH
    // INCHI❌❌:         while (i < n)
    // INCHI❌❌:         {
    // INCHI❌❌:             r = (rank_mask_bit & p->Rank[m = (int)p->AtNumber[i]]);
    // INCHI❌❌:             hash = add2crc32(hash, (AT_NUMB)(r + n));
    // INCHI❌❌:             r++;
    // INCHI❌❌:             insertions_sort_NeighList_AT_NUMBERS2(G[m], p->Rank, r);
    // INCHI❌❌:             nn = G[m][0];
    // INCHI❌❌:             rj_prev = 0; /* debug only */
    // INCHI❌❌:             for (j = 1; j <= nn && (rj = (rank_mask_bit & p->Rank[(int)G[m][j]])) < r; j++)
    // INCHI❌❌:             {
    // INCHI❌❌:                 hash = add2crc32(hash, rj);
    // INCHI❌❌:             }
    // INCHI❌❌:             i++;
    // INCHI❌❌:         }
    // INCHI❌❌:         Ct->hash[k] = hash;
    // INCHI❌❌: #endif
    // INCHI✔️❌:     }
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CtPartFill
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CtPartFill
    // INCHI✔️❌: #define bRELEASE_VERSION  1    /* 1=> release version; comment out to disable */
    // INCHI✔️❌: #define USE_ISO_SORT_KEY_HFIXED  0
    // INCHI✔️❌: /* #define INCHI_CANON_USE_HASH */
    // INCHI✔️❌: #define INCHI_HEAPCHK
    // END INCHI ACTIVE MACRO CONFIGURATION: CtPartFill

    k = k.wrapping_sub(1);
    let Some(ct) = Ct else {
        return Ok(());
    };
    let (mut start_ctbl, start_at_ord) = if k != 0 {
        (
            i32::from(source_get(heap, ct.nextCtblPos, k.wrapping_sub(1))?),
            i32::from(source_get(heap, ct.nextAtRank, k.wrapping_sub(1))?).wrapping_sub(1),
        )
    } else {
        (0, 0)
    };
    if start_at_ord < 0 || start_at_ord >= n_max {
        return Ok(());
    }

    // SAFETY: CanonGraph constructs the connection table, graph pointer table,
    // graph rows, atom order, and ranks as distinct fixed-size allocations.
    // The nested sort only mutates one graph row and never frees or resizes it.
    let graph_storage = unsafe { heap.stable_slice(G.as_const())? };
    let graph = graph_storage.prefix(graph_storage.len())?;
    let atom_number_storage = unsafe { heap.stable_slice(p.AtNumber.as_const())? };
    let atom_numbers = atom_number_storage.prefix(atom_number_storage.len())?;
    let rank_storage = unsafe { heap.stable_slice(p.Rank.as_const())? };
    let ranks = rank_storage.prefix(rank_storage.len())?;
    let mut connection_table_storage = unsafe { heap.stable_slice_mut(ct.Ctbl)? };
    let connection_table_len = connection_table_storage.len();
    let connection_table = connection_table_storage.prefix_mut(connection_table_len)?;
    let rank_mask = rank_mask_bit();

    // INCHI✔️✔️:         an_sao = (int)p->AtNumber[startAtOrd];
    let first_atom_index = start_at_ord as usize;
    if first_atom_index >= atom_numbers.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    // SAFETY: `start_at_ord >= 0` above and the explicit length check prove
    // the exact source access before it occurs.
    let first_atom = i32::from(unsafe { *atom_numbers.get_unchecked(first_atom_index) });
    // INCHI✔️✔️:         r = (rank_mask_bit & p->Rank[an_sao]);
    let first_rank_index = first_atom as usize;
    if first_rank_index >= ranks.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    // SAFETY: atom numbers are unsigned and the rank length check proves the
    // same indexed read as `p->Rank[an_sao]`.
    let mut r = rank_mask & unsafe { *ranks.get_unchecked(first_rank_index) };
    // INCHI✔️✔️:         for (i = startAtOrd; i < n_tg &&
    // INCHI✔️✔️:              r == (rank_mask_bit & p->Rank[m = (int)p->AtNumber[i]]);
    // INCHI✔️✔️:              i++, r++)
    let mut i = start_at_ord;
    while i < n_tg {
        let index = i as usize;
        if i < 0 || index >= atom_numbers.len() {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        // SAFETY: the immediately preceding check proves this source access.
        let m = unsafe { *atom_numbers.get_unchecked(index) };
        let m_index = usize::from(m);
        if m_index >= ranks.len() {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        // SAFETY: `m_index` was checked against the stable rank allocation.
        if r != rank_mask & unsafe { *ranks.get_unchecked(m_index) } {
            break;
        }
        // INCHI✔️✔️:             Ct->Ctbl[startCtbl++] = r;
        let output_index = start_ctbl as usize;
        if start_ctbl < 0 || output_index >= connection_table.len() {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        // SAFETY: the output index is checked at the source write point.
        unsafe { *connection_table.get_unchecked_mut(output_index) = r };
        start_ctbl = start_ctbl.wrapping_add(1);
        if m_index >= graph.len() {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        // SAFETY: graph and rank use the same checked unsigned atom index.
        let neighbors = unsafe { *graph.get_unchecked(m_index) };
        // INCHI✔️✔️:             insertions_sort_NeighList_AT_NUMBERS2(G[m], p->Rank, r);
        let mut neighbor_storage =
            match try_sort_neighbor_row_source_layout(heap, neighbors, ranks, r, rank_mask)? {
                Some(storage) => storage,
                None => {
                    insertions_sort_NeighList_AT_NUMBERS2(heap, neighbors, p.Rank, r)?;
                    unsafe { heap.stable_slice_mut(neighbors)? }
                }
            };
        let neighbor_len = neighbor_storage.len();
        let neighbors = neighbor_storage.prefix_mut(neighbor_len)?;
        // INCHI✔️✔️:             nn = G[m][0];
        if neighbors.is_empty() {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        // SAFETY: the graph row was checked nonempty at the source read point.
        let nn = unsafe { *neighbors.get_unchecked(0) };
        // INCHI✔️✔️:             for (j = 1; j <= nn &&
        // INCHI✔️✔️:                  (rj = (rank_mask_bit & p->Rank[(int)G[m][j]])) < r; j++)
        let mut j: AT_RANK = 1;
        while j <= nn {
            let neighbor_index = usize::from(j);
            if neighbor_index >= neighbors.len() {
                return Err(SourceHeapError::PointerOutOfBounds);
            }
            // SAFETY: the graph-row index is checked immediately above.
            let neighbor = unsafe { *neighbors.get_unchecked(neighbor_index) };
            let neighbor_rank_index = usize::from(neighbor);
            if neighbor_rank_index >= ranks.len() {
                return Err(SourceHeapError::PointerOutOfBounds);
            }
            // SAFETY: the neighbor atom index is checked against ranks.
            let rj = rank_mask & unsafe { *ranks.get_unchecked(neighbor_rank_index) };
            if rj >= r {
                break;
            }
            // INCHI✔️✔️:                 Ct->Ctbl[startCtbl++] = rj;
            let output_index = start_ctbl as usize;
            if start_ctbl < 0 || output_index >= connection_table.len() {
                return Err(SourceHeapError::PointerOutOfBounds);
            }
            // SAFETY: the output index is checked at the source write point.
            unsafe { *connection_table.get_unchecked_mut(output_index) = rj };
            start_ctbl = start_ctbl.wrapping_add(1);
            j = j.wrapping_add(1);
        }
        i = i.wrapping_add(1);
        r = r.wrapping_add(1);
    }

    // INCHI✔️✔️:         if (pCD->NumH && Ct->NumH)
    if !pCD.NumH.is_null() && !ct.NumH.is_null() {
        // SAFETY: CTableCreate gives each output field its own fixed-size
        // allocation; CANON_DATA owns distinct immutable source fields.
        let source_h_storage = unsafe { heap.stable_slice(pCD.NumH.as_const())? };
        let source_h = source_h_storage.prefix(source_h_storage.len())?;
        let mut target_h_storage = unsafe { heap.stable_slice_mut(ct.NumH)? };
        let target_h_len = target_h_storage.len();
        let target_h = target_h_storage.prefix_mut(target_h_len)?;
        // INCHI✔️✔️:             nn = inchi_min(n, i);
        let mut nn = n.min(i) as AT_RANK;
        // INCHI✔️✔️:             for (j = startAtOrd; j < nn; j++)
        let mut j = start_at_ord as AT_RANK;
        while j < nn {
            let atom = atom_numbers
                .get(usize::from(j))
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            // INCHI✔️✔️:                 Ct->NumH[j] = pCD->NumH[p->AtNumber[j]];
            *target_h
                .get_mut(usize::from(j))
                .ok_or(SourceHeapError::PointerOutOfBounds)? = source_h
                .get(usize::from(atom))
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            j = j.wrapping_add(1);
        }
        // INCHI✔️✔️:             for (; j < i; j++)
        while i32::from(j) < i {
            let atom = i32::from(
                atom_numbers
                    .get(usize::from(j))
                    .copied()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            // INCHI✔️✔️:                 int data_pos = n + T_NUM_NO_ISOTOPIC * ((int)p->AtNumber[j] - n);
            let mut data_pos = n.wrapping_add(
                (crate::source_types::T_NUM_NO_ISOTOPIC as i32).wrapping_mul(atom.wrapping_sub(n)),
            );
            // INCHI✔️✔️:                 for (m = 0; m < T_NUM_NO_ISOTOPIC; m++)
            let mut m = 0_i32;
            while m < crate::source_types::T_NUM_NO_ISOTOPIC as i32 {
                let source_index =
                    usize::try_from(data_pos).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                // INCHI✔️✔️:                     Ct->NumH[nn++] = pCD->NumH[data_pos++];
                *target_h
                    .get_mut(usize::from(nn))
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = source_h
                    .get(source_index)
                    .copied()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                nn = nn.wrapping_add(1);
                data_pos = data_pos.wrapping_add(1);
                m = m.wrapping_add(1);
            }
            j = j.wrapping_add(1);
        }
        // INCHI✔️✔️:             Ct->lenNumH = nn;
        ct.lenNumH = i32::from(nn);
    } else {
        // INCHI✔️✔️:             Ct->lenNumH = 0;
        ct.lenNumH = 0;
    }

    // INCHI✔️✔️:             if (pCD->NumHfixed && Ct->NumHfixed)
    if !pCD.NumHfixed.is_null() && !ct.NumHfixed.is_null() {
        let source_fixed_h_storage = unsafe { heap.stable_slice(pCD.NumHfixed.as_const())? };
        let source_fixed_h = source_fixed_h_storage.prefix(source_fixed_h_storage.len())?;
        let mut target_fixed_h_storage = unsafe { heap.stable_slice_mut(ct.NumHfixed)? };
        let target_fixed_h_len = target_fixed_h_storage.len();
        let target_fixed_h = target_fixed_h_storage.prefix_mut(target_fixed_h_len)?;
        // INCHI✔️✔️:                 nn = inchi_min(n, i);
        let nn = n.min(i) as AT_RANK;
        // INCHI✔️✔️:                 for (j = startAtOrd; j < nn; j++)
        let mut j = start_at_ord as AT_RANK;
        while j < nn {
            let atom = atom_numbers
                .get(usize::from(j))
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            // INCHI✔️✔️:                     Ct->NumHfixed[j] = pCD->NumHfixed[p->AtNumber[j]];
            *target_fixed_h
                .get_mut(usize::from(j))
                .ok_or(SourceHeapError::PointerOutOfBounds)? = source_fixed_h
                .get(usize::from(atom))
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            j = j.wrapping_add(1);
        }
    }

    // INCHI✔️✔️:             if (pCD->iso_sort_key && Ct->iso_sort_key)
    if !pCD.iso_sort_key.is_null() && !ct.iso_sort_key.is_null() {
        let source_iso_storage = unsafe { heap.stable_slice(pCD.iso_sort_key.as_const())? };
        let source_iso = source_iso_storage.prefix(source_iso_storage.len())?;
        let mut target_iso_storage = unsafe { heap.stable_slice_mut(ct.iso_sort_key)? };
        let target_iso_len = target_iso_storage.len();
        let target_iso = target_iso_storage.prefix_mut(target_iso_len)?;
        // INCHI✔️✔️:                 for (j = startAtOrd; j < i; j++)
        let mut j = start_at_ord as AT_RANK;
        while i32::from(j) < i {
            let atom = atom_numbers
                .get(usize::from(j))
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            // INCHI✔️✔️:                     Ct->iso_sort_key[j] = pCD->iso_sort_key[p->AtNumber[j]];
            *target_iso
                .get_mut(usize::from(j))
                .ok_or(SourceHeapError::PointerOutOfBounds)? = source_iso
                .get(usize::from(atom))
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            j = j.wrapping_add(1);
        }
        // INCHI✔️✔️:                 Ct->len_iso_sort_key = i;
        ct.len_iso_sort_key = i;
    } else {
        // INCHI✔️✔️:                 Ct->len_iso_sort_key = 0;
        ct.len_iso_sort_key = 0;
    }

    // INCHI✔️✔️:             if (pCD->iso_exchg_atnos && Ct->iso_exchg_atnos)
    if !pCD.iso_exchg_atnos.is_null() && !ct.iso_exchg_atnos.is_null() {
        let source_exchange_storage = unsafe { heap.stable_slice(pCD.iso_exchg_atnos.as_const())? };
        let source_exchange = source_exchange_storage.prefix(source_exchange_storage.len())?;
        let mut target_exchange_storage = unsafe { heap.stable_slice_mut(ct.iso_exchg_atnos)? };
        let target_exchange_len = target_exchange_storage.len();
        let target_exchange = target_exchange_storage.prefix_mut(target_exchange_len)?;
        // INCHI✔️✔️:                 for (j = startAtOrd; j < i; j++)
        let mut j = start_at_ord as AT_RANK;
        while i32::from(j) < i {
            let atom = atom_numbers
                .get(usize::from(j))
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            // INCHI✔️✔️:                     Ct->iso_exchg_atnos[j] = pCD->iso_exchg_atnos[p->AtNumber[j]];
            *target_exchange
                .get_mut(usize::from(j))
                .ok_or(SourceHeapError::PointerOutOfBounds)? = source_exchange
                .get(usize::from(atom))
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            j = j.wrapping_add(1);
        }
        // INCHI✔️✔️:                 Ct->len_iso_exchg_atnos = i;
        ct.len_iso_exchg_atnos = i;
    } else {
        // INCHI✔️✔️:                 Ct->len_iso_exchg_atnos = 0;
        ct.len_iso_exchg_atnos = 0;
    }

    ct.lenCt = start_ctbl;
    source_set(heap, ct.nextCtblPos, k, start_ctbl as AT_NUMB)?;
    source_set(heap, ct.nextAtRank, k, r)?;
    ct.lenPos = k.wrapping_add(1);
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn CtPartINCHI_CANON_INFINITY(
    heap: &mut SourceHeap,
    Ct: &ConTable,
    cmp: SourceMutPointer<crate::source_types::S_CHAR>,
    mut k: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2065 CtPartINCHI_CANON_INFINITY
    // INCHI✔️❌: void CtPartINCHI_CANON_INFINITY( ConTable *Ct, S_CHAR *cmp, int k )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int     startCtbl;
    // INCHI✔️❌:     /*int     startAtOrd;*/
    // INCHI✔️❌:     k--;
    // INCHI✔️❌:     if (k)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         startCtbl = Ct->nextCtblPos[k - 1];
    // INCHI✔️❌:         /*startAtOrd = Ct->nextAtRank[k-1]-1;*/  /* here  p->Rank[p->AtNumber[r-1]] = r */
    // INCHI✔️❌:         if (cmp)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             memset( cmp, 0, k * sizeof( cmp[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         startCtbl = 0;
    // INCHI✔️❌:         /*startAtOrd = 0;*/
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (!startCtbl || Ct->Ctbl[startCtbl - 1] != EMPTY_CT)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Ct->Ctbl[startCtbl] = EMPTY_CT;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CtPartINCHI_CANON_INFINITY

    k = k.wrapping_sub(1);
    let start_ctbl = if k != 0 {
        let start = i32::from(source_get(heap, Ct.nextCtblPos, k.wrapping_sub(1))?);
        if !cmp.is_null() {
            let count = usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            heap.slice_mut(cmp)?
                .get_mut(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .fill(0);
        }
        start
    } else {
        0
    };
    if start_ctbl == 0
        || source_get(heap, Ct.Ctbl, start_ctbl.wrapping_sub(1))? != EMPTY_CT as AT_RANK
    {
        source_set(heap, Ct.Ctbl, start_ctbl, EMPTY_CT as AT_RANK)?;
    }
    Ok(())
}

#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn CtPartCompare(
    heap: &mut SourceHeap,
    Ct1: &ConTable,
    Ct2: &ConTable,
    cmp: SourceMutPointer<crate::source_types::S_CHAR>,
    mut kLeastForLayer: Option<&mut KLeastLayers>,
    mut k: i32,
    bOnlyCommon: i32,
    bSplitTautCompare: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2123 CtPartCompare
    // INCHI✔️❌: int CtPartCompare( ConTable *Ct1,
    // INCHI✔️❌:                    ConTable *Ct2,
    // INCHI✔️❌:                    S_CHAR *cmp,
    // INCHI✔️❌:                    kLeast *kLeastForLayer,
    // INCHI✔️❌:                    int k,
    // INCHI✔️❌:                    int bOnlyCommon,
    // INCHI✔️❌:                    int bSplitTautCompare )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int     startCt1, endCt1, startCt2, endCt2; /*endCt,*/
    // INCHI✔️❌:     int     startAt1, endAt1, startAt2, endAt2; /*endCt,*/
    // INCHI✔️❌:     int     midCt /* end of atoms only Ct */, midNumH = 0 /* end of atoms only NumH */, maxVert;
    // INCHI✔️❌:     int     diff, i, k1, k2, lenNumH, /*mid_iso_sort_key,*/ midAt; /* djb-rwth: ignoring LLVM warning: variables used */
    // INCHI✔️❌:     int     nLayer = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     k--;
    // INCHI✔️❌:     i = -1;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* set kLeastForLayer[nLayer].k = (k+1) or -(k+1)
    // INCHI✔️❌:            kLeastForLayer[nLayer].i = iDiff
    // INCHI✔️❌:         if all the conditions are met:
    // INCHI✔️❌:         1) kLeastForLayer[nLayer].k = 0
    // INCHI✔️❌:         2) diff==0 for all layers < nLayer
    // INCHI✔️❌:
    // INCHI✔️❌:         sign:
    // INCHI✔️❌:         if the final diff < 0 then kLeastForLayer[nLayer].k = -(k+1) else
    // INCHI✔️❌:         if the final diff > 0 then kLeastForLayer[nLayer].k = +(k+1)
    // INCHI✔️❌:
    // INCHI✔️❌:         k+1 instead of k takes into account k--; statememt above)
    // INCHI✔️❌:
    // INCHI✔️❌:         meaning:
    // INCHI✔️❌:         ========
    // INCHI✔️❌:         abs(kLeastForLayer[nLayer].k) is the greatest level k at which
    // INCHI✔️❌:         difference at layer nLayer are zeroes of hidden by differences in smaller nLayer.
    // INCHI✔️❌:
    // INCHI✔️❌:         "Hidden by difference in smaller level" means that nLayer of comparison
    // INCHI✔️❌:         has not been reached because the difference was discovered at a previous layer.
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         Lambda vs zf_zeta comparison
    // INCHI✔️❌:         =============================================
    // INCHI✔️❌:         accept only diff == 0
    // INCHI✔️❌:
    // INCHI✔️❌:         Lambda vs pzb_rho and pzb_rho_fix comparison
    // INCHI✔️❌:         =============================================
    // INCHI✔️❌:         Maintain kLeastForLayer[] and kLeastForLayerFix[]
    // INCHI✔️❌:
    // INCHI✔️❌:         The algorithm provides that pzb_rho(m-1) < pzb_rho(m) <= pzb_rho_fix
    // INCHI✔️❌:
    // INCHI✔️❌:         Definition: pzb_rho(m-1) < pzb_rho(m) means that
    // INCHI✔️❌:         -----------------------------------------------
    // INCHI✔️❌:         pzb_rho(m-1)[nLayerCurr] == pzb_rho(m)[nLayerCurr] for nLayerCurr = 0..nLayerDiff-1
    // INCHI✔️❌:         pzb_rho(m-1)[nLayerDiff] <  pzb_rho(m)[nLayerDiff]
    // INCHI✔️❌:
    // INCHI✔️❌:         Definition: pzb_rho(m-1)[nLayerDiff] <  pzb_rho(m)[nLayerDiff] means that
    // INCHI✔️❌:         -------------------------------------------------------------------------
    // INCHI✔️❌:         pzb_rho(m-1)[nLayerDiff][i]     == pzb_rho(m)[nLayerDiff][i] for i=0..iDdiff-1
    // INCHI✔️❌:         pzb_rho(m-1)[nLayerDiff][iDdiff] < pzb_rho(m)[nLayerDiff][iDdiff]
    // INCHI✔️❌:
    // INCHI✔️❌:         This defines nLayerDiff(pzb1, pzb2) where pszb1 = pzb_rho(a), pzb2=pzb_rho(b) (a<b) or pzb_rho_fix
    // INCHI✔️❌:                and   iDdiff    (pzb1, pzb2).
    // INCHI✔️❌:         In case pzb_rho(m)[nLayerCurr] == pzb_rho_fix[nLayerCurr] for all non-NULL nLayerCurr in pzb_rho_fix,
    // INCHI✔️❌:            nLayerDiff(pzb_rho(m), pzb_rho_fix) = the first layer in pzb_rho(m) not present in pzb_rho_fix
    // INCHI✔️❌:            iDdiff    (pzb_rho(m), pzb_rho_fix) = -1
    // INCHI✔️❌:         Case when such a layer does not exist means program error
    // INCHI✔️❌:
    // INCHI✔️❌:         Suppose L_rho = nLayerDiff(Lambda, pzb_rho(m))
    // INCHI✔️❌:                 L_fix = nLayerDiff(Lambda, pzb_rho_fix)
    // INCHI✔️❌:                 I_rho = iDdiff    (Lambda, pzb_rho(m))
    // INCHI✔️❌:                 I_fix = iDdiff    (Lambda, pzb_rho_fix)
    // INCHI✔️❌:                 kLeastForLayer determined from Lambda vs pzb_rho(m) comparison
    // INCHI✔️❌:         Then:
    // INCHI✔️❌:
    // INCHI✔️❌:         1. Comparison Lambda vs pzb_rho_fix before reaching discrete partition
    // INCHI✔️❌:         ----------------------------------------------------------------------
    // INCHI✔️❌:         a)    0 < abs(kLeastForLayerFix[L_fix].k) <= k-1 (* in this case I_fix >= 0 *)  &&
    // INCHI✔️❌:               ((L_fix < L_rho) || (L_fix == L_rho && I_fix < I_rho))
    // INCHI✔️❌:               =>
    // INCHI✔️❌:               qzb_rho_fix = kLeastForLayerFix[L_fix].k if prevoiusly qzb_rho_fix == 0
    // INCHI✔️❌:
    // INCHI✔️❌:         b)    otherwise do not change qzb_rho_fix, except the following:
    // INCHI✔️❌:
    // INCHI✔️❌:         c)    Special case L_rho == L_fix && I_rho == I_fix. Let L=L_rho, I = I_rho.
    // INCHI✔️❌:
    // INCHI✔️❌:               Compare 3 valirs: Lambda[L][I], pzb_rho(m)[L][I], pzb_rho_fix[L][I]
    // INCHI✔️❌:               The algorithm provides pzb_rho(m)[L][I] < pzb_rho_fix[L][I]
    // INCHI✔️❌:               (pzb_rho(m)[L][I]==pzb_rho_fix[L][I] <=> pzb_rho(m)[L][I]==pzb_rho_fix[L][I]
    // INCHI✔️❌:                is impossible by construction)
    // INCHI✔️❌:               There are 3 possibilities:
    // INCHI✔️❌:               c1) Lambda[L][I]     < pzb_rho(m)[L][I]  < pzb_rho_fix[L][I] <=>
    // INCHI✔️❌:                   kLeastForLayer[L].k  < 0 && kLeastForLayerFix[L].k < 0
    // INCHI✔️❌:                   => qzb_rho := kLeastForLayer[L].k, reject too small Lambda
    // INCHI✔️❌:               c2) pzb_rho(m)[L][I] < Lambda[L][I]      < pzb_rho_fix[L][I]
    // INCHI✔️❌:                   kLeastForLayer[L].k  > 0 && kLeastForLayerFix[L].k < 0
    // INCHI✔️❌:                   => qzb_rho := kLeastForLayer[L].k, accept Lambda, rho:=nu
    // INCHI✔️❌:               c3) pzb_rho(m)[L][I] < pzb_rho_fix[L][I] < Lambda[L][I]
    // INCHI✔️❌:                   kLeastForLayer[L].k  > 0 && kLeastForLayerFix[L].k > 0
    // INCHI✔️❌:                   => qzb_rho_fix := kLeastForLayerFix[L].k, reject too big Lambda
    // INCHI✔️❌:
    // INCHI✔️❌:               Case
    // INCHI✔️❌:                   kLeastForLayer[L].k  < 0 && kLeastForLayerFix[L].k > 0 is impossible
    // INCHI✔️❌:                   because it means
    // INCHI✔️❌:                   pzb_rho_fix < Lambda < pzb_rho(m) <=> pzb_rho_fix < pzb_rho(m)
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:             Case (c3) occurs in case of (a)
    // INCHI✔️❌:             Case (c1)
    // INCHI✔️❌:
    // INCHI✔️❌:         2. Comparison Lambda vs pzb_rho before reaching discrete partition
    // INCHI✔️❌:         ----------------------------------------------------------------------
    // INCHI✔️❌:         a) (L_rho < L_fix) || (L_rho == L_fix && I_rho < I_fix)  =>
    // INCHI✔️❌:
    // INCHI✔️❌:            Lambda differs from pzb_rho(m) in the part of pzb_rho(m) that will never change
    // INCHI✔️❌:            qzb_rho = kLeastForLayer[L_rho].k; reject Labmda or accept pzb_rho(m+1):=Labmda
    // INCHI✔️❌:
    // INCHI✔️❌:         b) (L_rho == L_fix && I_rho > I_fix) && kLeastForLayer[L_rho].k < 0
    // INCHI✔️❌:            Lambda < pzb_rho(m), therefore
    // INCHI✔️❌:            qzb_rho = kLeastForLayer[L_rho].k; reject Labmda
    // INCHI✔️❌:
    // INCHI✔️❌:         c) (L_rho > L_fix) =>
    // INCHI✔️❌:            qzb_rho := 0 because more significant difference may be discovered
    // INCHI✔️❌:            in layer < L_rho later. The final comparison may be needed at the
    // INCHI✔️❌:            level of discrete partition.
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (cmp)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = 0; i <= k && !cmp[i]; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i < k)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cmp[k] = cmp[i];
    // INCHI✔️❌:             return (int) cmp[i];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     k1 = Ct1->lenPos - 1; /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:     k2 = Ct2->lenPos - 1; /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:
    // INCHI❌❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI❌❌:     if (k > k1 || k > k2)
    // INCHI❌❌:     {
    // INCHI❌❌:         int stop = 1;
    // INCHI❌❌:     }
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     diff = 0; /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (k)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         startCt1 = Ct1->nextCtblPos[k - 1];
    // INCHI✔️❌:         startCt2 = Ct2->nextCtblPos[k - 1];
    // INCHI✔️❌:         startAt1 = Ct1->nextAtRank[k - 1] - 1;
    // INCHI✔️❌:         startAt2 = Ct2->nextAtRank[k - 1] - 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         startCt1 = startCt2 = 0;
    // INCHI✔️❌:         startAt1 = startAt2 = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     endCt1 = Ct1->nextCtblPos[k];
    // INCHI✔️❌:     endCt2 = Ct2->nextCtblPos[k];
    // INCHI✔️❌:     endAt1 = (int) Ct1->nextAtRank[k] - 1;
    // INCHI✔️❌:     endAt2 = (int) Ct2->nextAtRank[k] - 1;
    // INCHI✔️❌:
    // INCHI✔️❌:     maxVert = inchi_min( Ct1->maxVert, Ct2->maxVert );
    // INCHI✔️❌:
    // INCHI❌❌: #ifdef INCHI_CANON_USE_HASH
    // INCHI❌❌:     if (!diff)
    // INCHI❌❌:     {
    // INCHI❌❌:         if (Ct1->hash[k] > Ct2->hash[k])
    // INCHI❌❌:             diff = 1;
    // INCHI❌❌:         else
    // INCHI❌❌:             if (Ct1->hash[k] < Ct2->hash[k])
    // INCHI❌❌:                 diff = -1;
    // INCHI❌❌:     }
    // INCHI❌❌:     if (diff)
    // INCHI❌❌:     {
    // INCHI❌❌:         goto done;
    // INCHI❌❌:     }
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     /************************** lengths **************************************************/
    // INCHI✔️❌:     if ((diff = -( startCt1 - startCt2 ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* comparing two INCHI_CANON_INFINITY terminations */
    // INCHI✔️❌:         if (bOnlyCommon &&
    // INCHI✔️❌:              startCt1 >= Ct1->nLenCTAtOnly && startCt2 >= Ct2->nLenCTAtOnly &&
    // INCHI✔️❌:              Ct1->Ctbl[startCt1] == EMPTY_CT && Ct2->Ctbl[startCt2] == EMPTY_CT)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (bOnlyCommon)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             startCt1 = inchi_min( startCt1, startCt2 ); /* djb-rwth: removing redundant code */
    // INCHI✔️❌:             startAt1 = startAt2 = inchi_min( startAt1, startAt2 );
    // INCHI✔️❌:             if (Ct1->lenCt == Ct2->lenCt)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 endCt1 = endCt2 = inchi_max( endCt1, endCt2 );
    // INCHI✔️❌:                 endAt1 = endAt2 = inchi_max( endAt1, endAt2 );
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI❌❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI❌❌:             else
    // INCHI❌❌:             {
    // INCHI❌❌:                 int stop = 1;
    // INCHI❌❌:             }
    // INCHI❌❌: #endif
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         /* comparing (taut tail) vs INCHI_CANON_INFINITY termination -- ??? */
    // INCHI✔️❌:             if (startCt1 > startCt2 &&
    // INCHI✔️❌:                  Ct1->maxVert > Ct2->maxVert &&
    // INCHI✔️❌:                  startAt2 == Ct2->maxVert)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto done;
    // INCHI✔️❌:             }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     lenNumH = Ct1->lenNumH;
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:
    // INCHI✔️❌:     if ((diff = -( endCt1 - endCt2 ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* negative sign reproduces results for NSC=28393 */
    // INCHI✔️❌:         if (bOnlyCommon)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             endCt1 = inchi_min( endCt1, endCt2 ); /* djb-rwth: removing redundant code */
    // INCHI✔️❌:             endAt1 = endAt2 = inchi_min( endAt1, endAt2 );
    // INCHI✔️❌:             lenNumH = inchi_min( Ct1->lenNumH, Ct2->lenNumH );
    // INCHI✔️❌:             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* take care of case when comparing tautomeric vs non-tautomeric:
    // INCHI✔️❌:             since (taut)->maxVert > (non-taut)->maxVert, --???
    // INCHI✔️❌:             (taut)->maxlenCt  > (non-taut)->maxlenCt     --!!!
    // INCHI✔️❌:             compare up to min out of the two, ignoring INCHI_CANON_INFINITY in the last position */
    // INCHI✔️❌:             if (endCt1 > endCt2 && Ct1->maxlenCt > Ct2->maxlenCt)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (endAt2 == Ct2->maxVert + 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* remove INCHI_CANON_INFINITY termination of the shorter CT */
    // INCHI✔️❌:                     /* should never happen */
    // INCHI✔️❌:                     endAt2--;
    // INCHI✔️❌:                     lenNumH = endAt1 = endAt2; /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                     endCt2--;
    // INCHI✔️❌:                     endCt1 = endCt2;
    // INCHI✔️❌:                     diff = 0; /* djb-rwth: ignoring LLVM warning: value used? */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else if (endAt2 == Ct2->maxVert)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* remove INCHI_CANON_INFINITY termination of CT */
    // INCHI✔️❌:                     lenNumH = endAt1 = endAt2; /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                     endCt1 = endCt2;
    // INCHI✔️❌:                     diff = 0; /* djb-rwth: ignoring LLVM warning: value used? */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     goto done;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto done;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (bSplitTautCompare)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         midCt = inchi_min( Ct1->nLenCTAtOnly, Ct2->nLenCTAtOnly );
    // INCHI✔️❌:         if (midCt > endCt1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             midCt = endCt1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         midAt = inchi_min( maxVert, endAt1 );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         midCt = endCt1;
    // INCHI✔️❌:         midAt = endAt1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*endCt   = min(endCt1, endCt2);*/
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:     /************ layer 0: connection table without tautomeric groups ********/
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = startCt1; i < midCt && Ct1->Ctbl[i] == Ct2->Ctbl[i]; i++)
    // INCHI✔️❌:     /*for ( i = startCt1; i < endCt && !(diff = (int)Ct1->Ctbl[i] - (int)Ct2->Ctbl[i]); i ++ )*/
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i < midCt)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         diff = (int) Ct1->Ctbl[i] - (int) Ct2->Ctbl[i];
    // INCHI✔️❌:         goto done;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:     /******** layer 1 NumH: H atoms without tautomeric H *********************/
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:     nLayer++;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*============= check limits for consistency  ==========*/
    // INCHI✔️❌:     if ((diff = -( startAt1 - startAt2 ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto done;   /* should not happen */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if ((diff = -( endAt1 - endAt2 ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto done;   /* should not happen */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*============= comparison =============================*/
    // INCHI✔️❌:     if (Ct1->NumH && Ct2->NumH)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (endAt1 < maxVert)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             midNumH = lenNumH = endAt1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (bSplitTautCompare)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             midNumH = maxVert;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             midNumH = lenNumH;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* lenNumH = (endAt2 >= maxVert)? lenNumH : endAt2; */
    // INCHI✔️❌:         /* endAt1 = (endAt2 == n)? lenNumH : endAt2; */
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = startAt1; i < midNumH && Ct1->NumH[i] == Ct2->NumH[i]; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i < midNumH)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             diff = (int) Ct1->NumH[i] - (int) Ct2->NumH[i];
    // INCHI✔️❌:             goto done;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:     /************** layer 2: tautomeric part of CT and tautomeric H **********/
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:     nLayer++;
    // INCHI✔️❌:     for (i = midCt; i < endCt1 && Ct1->Ctbl[i] == Ct2->Ctbl[i]; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ; /* compare tautomeric groups part of CT */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (i < endCt1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         diff = (int) Ct1->Ctbl[i] - (int) Ct2->Ctbl[i];
    // INCHI✔️❌:         goto done;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (Ct1->NumH && Ct2->NumH)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = midNumH; i < lenNumH && Ct1->NumH[i] == Ct2->NumH[i]; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ; /* compare tautomeric H */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i < lenNumH)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             diff = (int) Ct1->NumH[i] - (int) Ct2->NumH[i];
    // INCHI✔️❌:             i += endCt1 - midCt;
    // INCHI✔️❌:             goto done;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:     /************** layer 3: Fixed H atoms ***********************************/
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:     nLayer++;
    // INCHI✔️❌:     if (Ct1->NumHfixed && Ct2->NumHfixed)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = startAt1; i < midAt && Ct1->NumHfixed[i] == Ct2->NumHfixed[i]; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i < midAt)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             diff = (int) Ct1->NumHfixed[i] - (int) Ct2->NumHfixed[i];
    // INCHI✔️❌:             goto done;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:     /************** layer 4: isotopic atoms H, incl. tautomeric **************/
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:     nLayer++;
    // INCHI✔️❌:     if (Ct1->iso_sort_key && Ct2->iso_sort_key)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = startAt1; i < endAt1 && Ct1->iso_sort_key[i] == Ct2->iso_sort_key[i]; i++)
    // INCHI✔️❌:             ;
    // INCHI✔️❌:         if (i < endAt1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             diff = Ct1->iso_sort_key[i] > Ct2->iso_sort_key[i] ? 1 : -1;
    // INCHI✔️❌:             goto done;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (Ct1->iso_exchg_atnos && Ct2->len_iso_exchg_atnos)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = startAt1; i < endAt1 && Ct1->iso_exchg_atnos[i] == Ct2->iso_exchg_atnos[i]; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i < endAt1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             diff = Ct1->iso_exchg_atnos[i] > Ct2->iso_exchg_atnos[i] ? 1 : -1;
    // INCHI✔️❌:             goto done;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI❌❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI❌❌:
    // INCHI❌❌:     /*************************************************************************/
    // INCHI❌❌:     /************** layer 6: Fixed isotopic H atoms **************************/
    // INCHI❌❌:     /*************************************************************************/
    // INCHI❌❌:     nLayer++;
    // INCHI❌❌:     if (Ct1->iso_sort_key_Hfixed && Ct2->iso_sort_key_Hfixed)
    // INCHI❌❌:     {
    // INCHI❌❌:         for (i = startAt1; i < midAt && Ct1->iso_sort_key_Hfixed[i] == Ct2->iso_sort_key_Hfixed[i]; i++)
    // INCHI❌❌:         {
    // INCHI❌❌:             ;
    // INCHI❌❌:         }
    // INCHI❌❌:         if (i < midAt)
    // INCHI❌❌:         {
    // INCHI❌❌:             diff = Ct1->iso_sort_key_Hfixed[i] > Ct2->iso_sort_key_Hfixed[i] ? 1 : -1;
    // INCHI❌❌:             goto done;
    // INCHI❌❌:         }
    // INCHI❌❌:     }
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: done:
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef INCHI_CANON_MIN
    // INCHI✔️❌:     diff = -diff;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     if (diff)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         diff = ( diff > 0 ) ? ( nLayer + 1 ) : -( nLayer + 1 ); /* return the discovered difference layer number >= 1 */
    // INCHI✔️❌:         if (kLeastForLayer)
    // INCHI✔️❌:         {
    // INCHI✔️❌:
    // INCHI❌❌: #if ( bRELEASE_VERSION != 1 )
    // INCHI❌❌:             if (abs( kLeastForLayer[nLayer].k ) > k + 1)
    // INCHI❌❌:             { /* for debug only */
    // INCHI❌❌:                 int stop = 1; /* <BRKPT> */
    // INCHI❌❌:             }
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:             if (!kLeastForLayer[nLayer].k)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 kLeastForLayer[nLayer].k = ( diff > 0 ) ? ( k + 1 ) : -( k + 1 );
    // INCHI✔️❌:                 kLeastForLayer[nLayer].i = i;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (nLayer /* && !bOnlyCommon */)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 diff = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI❌❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI❌❌:     else
    // INCHI❌❌:     {
    // INCHI❌❌:         int stop = 1;  /* for debug only */
    // INCHI❌❌:     }
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     if (cmp)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         cmp[k] = ( diff > 0 ) ? 1 : ( diff < 0 ) ? -1 : 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return diff;
    // INCHI✔️❌: }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // END INCHI C FUNCTION: CtPartCompare
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CtPartCompare
    // INCHI✔️❌: #define INCHI_CANON_MIN
    // INCHI✔️❌: #define bRELEASE_VERSION  1
    // INCHI✔️❌: #define USE_ISO_SORT_KEY_HFIXED  0
    // INCHI✔️❌: /* #define INCHI_CANON_USE_HASH */
    // END INCHI ACTIVE MACRO CONFIGURATION: CtPartCompare

    // INCHI✔️✔️:     k--;
    k = k.wrapping_sub(1);
    // INCHI✔️✔️:     i = -1;
    let mut i = -1_i32;
    let mut layer = 0_i32;

    // SAFETY: CanonGraph allocates cmp independently from every CTable field
    // and keeps it fixed-size and live throughout this comparison.
    let mut comparisons = if cmp.is_null() {
        None
    } else {
        Some(unsafe { heap.stable_slice_mut(cmp)? })
    };
    // INCHI✔️✔️:     if (cmp)
    if let Some(values) = comparisons.as_mut() {
        // INCHI✔️✔️:         for (i = 0; i <= k && !cmp[i]; i++)
        i = 0;
        while i <= k {
            let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            if *values.get(index)? != 0 {
                break;
            }
            i = i.wrapping_add(1);
        }
        // INCHI✔️✔️:         if (i < k)
        if i < k {
            let source_index =
                usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let target_index =
                usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            // INCHI✔️✔️:             cmp[k] = cmp[i];
            let value = *values.get(source_index)?;
            *values.get_mut(target_index)? = value;
            // INCHI✔️✔️:             return (int) cmp[i];
            return Ok(i32::from(value));
        }
    }
    // INCHI✔️✔️:     k1 = Ct1->lenPos - 1;
    let _k1 = Ct1.lenPos.wrapping_sub(1);
    // INCHI✔️✔️:     k2 = Ct2->lenPos - 1;
    let _k2 = Ct2.lenPos.wrapping_sub(1);

    // SAFETY: CTableCreate owns each table field in a separate fixed-size
    // allocation. CtPartCompare only reads those fields, and CanonGraph keeps
    // them live without resize or free for the complete call.
    let ct_positions1 = unsafe { heap.stable_slice(Ct1.nextCtblPos.as_const())? };
    let ct_positions2 = unsafe { heap.stable_slice(Ct2.nextCtblPos.as_const())? };
    let atom_positions1 = unsafe { heap.stable_slice(Ct1.nextAtRank.as_const())? };
    let atom_positions2 = unsafe { heap.stable_slice(Ct2.nextAtRank.as_const())? };
    let (mut sc1, mut sa1, sc2, mut sa2) = if k != 0 {
        let previous =
            usize::try_from(k.wrapping_sub(1)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        // INCHI✔️✔️:         startCt1 = Ct1->nextCtblPos[k - 1];
        let start_ct1 = i32::from(*ct_positions1.get(previous)?);
        // INCHI✔️✔️:         startCt2 = Ct2->nextCtblPos[k - 1];
        let start_ct2 = i32::from(*ct_positions2.get(previous)?);
        // INCHI✔️✔️:         startAt1 = Ct1->nextAtRank[k - 1] - 1;
        let start_at1 = i32::from(*atom_positions1.get(previous)?).wrapping_sub(1);
        // INCHI✔️✔️:         startAt2 = Ct2->nextAtRank[k - 1] - 1;
        let start_at2 = i32::from(*atom_positions2.get(previous)?).wrapping_sub(1);
        (start_ct1, start_at1, start_ct2, start_at2)
    } else {
        // INCHI✔️✔️:         startCt1 = startCt2 = 0;
        // INCHI✔️✔️:         startAt1 = startAt2 = 0;
        (0, 0, 0, 0)
    };
    let current = usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    // INCHI✔️✔️:     endCt1 = Ct1->nextCtblPos[k];
    let mut ec1 = i32::from(*ct_positions1.get(current)?);
    // INCHI✔️✔️:     endCt2 = Ct2->nextCtblPos[k];
    let mut ec2 = i32::from(*ct_positions2.get(current)?);
    // INCHI✔️✔️:     endAt1 = (int) Ct1->nextAtRank[k] - 1;
    let mut ea1 = i32::from(*atom_positions1.get(current)?).wrapping_sub(1);
    // INCHI✔️✔️:     endAt2 = (int) Ct2->nextAtRank[k] - 1;
    let mut ea2 = i32::from(*atom_positions2.get(current)?).wrapping_sub(1);
    // INCHI✔️✔️:     maxVert = inchi_min( Ct1->maxVert, Ct2->maxVert );
    let max_vert = Ct1.maxVert.min(Ct2.maxVert);
    let mut len_h = Ct1.lenNumH;
    let mut mid_h = 0_i32;
    let mut table_views = None;

    let (mut diff, done_i, done_layer) = 'compare: {
        // INCHI✔️✔️:     if ((diff = -( startCt1 - startCt2 )))
        let mut diff = sc2.wrapping_sub(sc1);
        if diff != 0 {
            // INCHI✔️✔️:         if (bOnlyCommon && startCt1 >= Ct1->nLenCTAtOnly &&
            // INCHI✔️✔️:             startCt2 >= Ct2->nLenCTAtOnly && Ct1->Ctbl[startCt1] == EMPTY_CT &&
            // INCHI✔️✔️:             Ct2->Ctbl[startCt2] == EMPTY_CT)
            if bOnlyCommon != 0 && sc1 >= Ct1.nLenCTAtOnly && sc2 >= Ct2.nLenCTAtOnly {
                let tables1 = unsafe { heap.stable_slice(Ct1.Ctbl.as_const())? };
                let tables2 = unsafe { heap.stable_slice(Ct2.Ctbl.as_const())? };
                let start1 =
                    usize::try_from(sc1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let start2 =
                    usize::try_from(sc2).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if *tables1.get(start1)? == EMPTY_CT as AT_RANK
                    && *tables2.get(start2)? == EMPTY_CT as AT_RANK
                {
                    // INCHI✔️✔️:             return 0;
                    return Ok(0);
                }
                table_views = Some((tables1, tables2));
            }
            // INCHI✔️✔️:         if (bOnlyCommon)
            if bOnlyCommon != 0 {
                sc1 = sc1.min(sc2);
                sa1 = sa1.min(sa2);
                sa2 = sa1;
                if Ct1.lenCt == Ct2.lenCt {
                    ec1 = ec1.max(ec2);
                    ec2 = ec1;
                    ea1 = ea1.max(ea2);
                    ea2 = ea1;
                }
            } else if sc1 > sc2 && Ct1.maxVert > Ct2.maxVert && sa2 == Ct2.maxVert {
                return Ok(0);
            } else {
                break 'compare (diff, i, layer);
            }
        }
        // INCHI✔️✔️:     if ((diff = -( endCt1 - endCt2 )))
        diff = ec2.wrapping_sub(ec1);
        if diff != 0 {
            if bOnlyCommon != 0 {
                ec1 = ec1.min(ec2);
                ea1 = ea1.min(ea2);
                ea2 = ea1;
                len_h = Ct1.lenNumH.min(Ct2.lenNumH);
            } else if ec1 > ec2 && Ct1.maxlenCt > Ct2.maxlenCt {
                if ea2 == Ct2.maxVert.wrapping_add(1) {
                    ea2 = ea2.wrapping_sub(1);
                    ea1 = ea2;
                    len_h = ea2;
                    ec2 = ec2.wrapping_sub(1);
                    ec1 = ec2;
                    diff = 0;
                } else if ea2 == Ct2.maxVert {
                    ea1 = ea2;
                    len_h = ea2;
                    ec1 = ec2;
                    diff = 0;
                } else {
                    break 'compare (diff, i, layer);
                }
            } else {
                break 'compare (diff, i, layer);
            }
        }
        if table_views.is_none() {
            table_views = Some(
                (unsafe { heap.stable_slice(Ct1.Ctbl.as_const())? }, unsafe {
                    heap.stable_slice(Ct2.Ctbl.as_const())?
                }),
            );
        }
        let (tables1, tables2) = table_views
            .as_ref()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let (mid_ct, mid_at) = if bSplitTautCompare != 0 {
            (
                Ct1.nLenCTAtOnly.min(Ct2.nLenCTAtOnly).min(ec1),
                max_vert.min(ea1),
            )
        } else {
            (ec1, ea1)
        };
        // INCHI✔️✔️:     for (i = startCt1; i < midCt && Ct1->Ctbl[i] == Ct2->Ctbl[i]; i++)
        i = sc1;
        while i < mid_ct {
            let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            if *tables1.get(index)? != *tables2.get(index)? {
                break;
            }
            i = i.wrapping_add(1);
        }
        // INCHI✔️✔️:     if (i < midCt)
        if i < mid_ct {
            let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            // INCHI✔️✔️:         diff = (int) Ct1->Ctbl[i] - (int) Ct2->Ctbl[i];
            diff = i32::from(*tables1.get(index)?).wrapping_sub(i32::from(*tables2.get(index)?));
            break 'compare (diff, i, layer);
        }
        layer = layer.wrapping_add(1);
        diff = sa2.wrapping_sub(sa1);
        if diff != 0 {
            break 'compare (diff, i, layer);
        }
        diff = ea2.wrapping_sub(ea1);
        if diff != 0 {
            break 'compare (diff, i, layer);
        }
        let mut hydrogen_views = None;
        // INCHI✔️✔️:     if (Ct1->NumH && Ct2->NumH)
        if !Ct1.NumH.is_null() && !Ct2.NumH.is_null() {
            if ea1 < max_vert {
                mid_h = ea1;
                len_h = ea1;
            } else if bSplitTautCompare != 0 {
                mid_h = max_vert;
            } else {
                mid_h = len_h;
            }
            let hydrogens1 = unsafe { heap.stable_slice(Ct1.NumH.as_const())? };
            let hydrogens2 = unsafe { heap.stable_slice(Ct2.NumH.as_const())? };
            // INCHI✔️✔️:         for (i = startAt1; i < midNumH && Ct1->NumH[i] == Ct2->NumH[i]; i++)
            i = sa1;
            while i < mid_h {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if *hydrogens1.get(index)? != *hydrogens2.get(index)? {
                    break;
                }
                i = i.wrapping_add(1);
            }
            // INCHI✔️✔️:         if (i < midNumH)
            if i < mid_h {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                // INCHI✔️✔️:             diff = (int) Ct1->NumH[i] - (int) Ct2->NumH[i];
                diff = i32::from(*hydrogens1.get(index)?)
                    .wrapping_sub(i32::from(*hydrogens2.get(index)?));
                break 'compare (diff, i, layer);
            }
            hydrogen_views = Some((hydrogens1, hydrogens2));
        }
        layer = layer.wrapping_add(1);
        // INCHI✔️✔️:     for (i = midCt; i < endCt1 && Ct1->Ctbl[i] == Ct2->Ctbl[i]; i++)
        i = mid_ct;
        while i < ec1 {
            let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            if *tables1.get(index)? != *tables2.get(index)? {
                break;
            }
            i = i.wrapping_add(1);
        }
        // INCHI✔️✔️:     if (i < endCt1)
        if i < ec1 {
            let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            // INCHI✔️✔️:         diff = (int) Ct1->Ctbl[i] - (int) Ct2->Ctbl[i];
            diff = i32::from(*tables1.get(index)?).wrapping_sub(i32::from(*tables2.get(index)?));
            break 'compare (diff, i, layer);
        }
        if let Some((hydrogens1, hydrogens2)) = &hydrogen_views {
            // INCHI✔️✔️:         for (i = midNumH; i < lenNumH && Ct1->NumH[i] == Ct2->NumH[i]; i++)
            i = mid_h;
            while i < len_h {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if *hydrogens1.get(index)? != *hydrogens2.get(index)? {
                    break;
                }
                i = i.wrapping_add(1);
            }
            // INCHI✔️✔️:         if (i < lenNumH)
            if i < len_h {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                // INCHI✔️✔️:             diff = (int) Ct1->NumH[i] - (int) Ct2->NumH[i];
                diff = i32::from(*hydrogens1.get(index)?)
                    .wrapping_sub(i32::from(*hydrogens2.get(index)?));
                // INCHI✔️✔️:             i += endCt1 - midCt;
                i = i.wrapping_add(ec1.wrapping_sub(mid_ct));
                break 'compare (diff, i, layer);
            }
        }
        layer = layer.wrapping_add(1);
        // INCHI✔️✔️:     if (Ct1->NumHfixed && Ct2->NumHfixed)
        if !Ct1.NumHfixed.is_null() && !Ct2.NumHfixed.is_null() {
            let fixed_hydrogens1 = unsafe { heap.stable_slice(Ct1.NumHfixed.as_const())? };
            let fixed_hydrogens2 = unsafe { heap.stable_slice(Ct2.NumHfixed.as_const())? };
            // INCHI✔️✔️:         for (i = startAt1; i < midAt && Ct1->NumHfixed[i] == Ct2->NumHfixed[i]; i++)
            i = sa1;
            while i < mid_at {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if *fixed_hydrogens1.get(index)? != *fixed_hydrogens2.get(index)? {
                    break;
                }
                i = i.wrapping_add(1);
            }
            // INCHI✔️✔️:         if (i < midAt)
            if i < mid_at {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                // INCHI✔️✔️:             diff = (int) Ct1->NumHfixed[i] - (int) Ct2->NumHfixed[i];
                diff = i32::from(*fixed_hydrogens1.get(index)?)
                    .wrapping_sub(i32::from(*fixed_hydrogens2.get(index)?));
                break 'compare (diff, i, layer);
            }
        }
        layer = layer.wrapping_add(1);
        // INCHI✔️✔️:     if (Ct1->iso_sort_key && Ct2->iso_sort_key)
        if !Ct1.iso_sort_key.is_null() && !Ct2.iso_sort_key.is_null() {
            let isotope_keys1 = unsafe { heap.stable_slice(Ct1.iso_sort_key.as_const())? };
            let isotope_keys2 = unsafe { heap.stable_slice(Ct2.iso_sort_key.as_const())? };
            // INCHI✔️✔️:         for (i = startAt1; i < endAt1 && Ct1->iso_sort_key[i] == Ct2->iso_sort_key[i]; i++)
            i = sa1;
            while i < ea1 {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if *isotope_keys1.get(index)? != *isotope_keys2.get(index)? {
                    break;
                }
                i = i.wrapping_add(1);
            }
            // INCHI✔️✔️:         if (i < endAt1)
            if i < ea1 {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                // INCHI✔️✔️:             diff = Ct1->iso_sort_key[i] > Ct2->iso_sort_key[i] ? 1 : -1;
                diff = if *isotope_keys1.get(index)? > *isotope_keys2.get(index)? {
                    1
                } else {
                    -1
                };
                break 'compare (diff, i, layer);
            }
        }
        // INCHI✔️✔️:     if (Ct1->iso_exchg_atnos && Ct2->len_iso_exchg_atnos)
        if !Ct1.iso_exchg_atnos.is_null() && Ct2.len_iso_exchg_atnos != 0 {
            let exchange_atoms1 = unsafe { heap.stable_slice(Ct1.iso_exchg_atnos.as_const())? };
            let exchange_atoms2 = unsafe { heap.stable_slice(Ct2.iso_exchg_atnos.as_const())? };
            // INCHI✔️✔️:         for (i = startAt1; i < endAt1 && Ct1->iso_exchg_atnos[i] == Ct2->iso_exchg_atnos[i]; i++)
            i = sa1;
            while i < ea1 {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if *exchange_atoms1.get(index)? != *exchange_atoms2.get(index)? {
                    break;
                }
                i = i.wrapping_add(1);
            }
            // INCHI✔️✔️:         if (i < endAt1)
            if i < ea1 {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                // INCHI✔️✔️:             diff = Ct1->iso_exchg_atnos[i] > Ct2->iso_exchg_atnos[i] ? 1 : -1;
                diff = if *exchange_atoms1.get(index)? > *exchange_atoms2.get(index)? {
                    1
                } else {
                    -1
                };
                break 'compare (diff, i, layer);
            }
        }
        (0, i, layer)
    };
    i = done_i;
    layer = done_layer;

    diff = diff.wrapping_neg();
    if diff != 0 {
        diff = if diff > 0 {
            layer.wrapping_add(1)
        } else {
            layer.wrapping_add(1).wrapping_neg()
        };
        if let Some(values) = kLeastForLayer.as_mut() {
            let index = usize::try_from(layer).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let value = values
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if value.k == 0 {
                *values
                    .get_mut(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = kLeast {
                    k: if diff > 0 {
                        k.wrapping_add(1)
                    } else {
                        k.wrapping_add(1).wrapping_neg()
                    },
                    i,
                };
            }
            if layer != 0 {
                diff = 0;
            }
        }
    }
    // INCHI✔️✔️:     if (cmp)
    if let Some(values) = comparisons.as_mut() {
        let index = usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        // INCHI✔️✔️:         cmp[k] = ( diff > 0 ) ? 1 : ( diff < 0 ) ? -1 : 0;
        *values.get_mut(index)? = if diff > 0 {
            1
        } else if diff < 0 {
            -1
        } else {
            0
        };
    }
    Ok(diff)
}

#[allow(non_snake_case)]
pub(crate) fn CtFullCompare(
    heap: &mut SourceHeap,
    Ct1: &ConTable,
    Ct2: &ConTable,
    bOnlyCommon: i32,
    bSplitTautCompare: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2623 CtFullCompare
    // INCHI✔️❌: int CtFullCompare( ConTable *Ct1,
    // INCHI✔️❌:                    ConTable *Ct2,
    // INCHI✔️❌:                    int bOnlyCommon,
    // INCHI✔️❌:                    int bSplitTautCompare )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int     startCt1, endCt1, endCt2; /*endCt,*/
    // INCHI✔️❌:     int     startAt1, endAt1, endAt2; /*endCt,*/
    // INCHI✔️❌:     int     midCt   /* end of atoms only in Ctbl */,
    // INCHI✔️❌:         midNumH = 0 /* end of atoms only NumH */,
    // INCHI✔️❌:         midAt   /* end of atoms only */; /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:     int     diff = 0, i, k1, k2, lenNumH1, lenNumH2, maxVert /* min num atoms */;
    // INCHI✔️❌:     int     len_iso_sort_key1, len_iso_sort_key2 /*, mid_iso_sort_key*/;
    // INCHI✔️❌:     int     nLayer = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     k1 = Ct1->lenPos - 1;
    // INCHI✔️❌:     k2 = Ct2->lenPos - 1;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:
    // INCHI✔️❌:     startCt1 = 0; /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     startAt1 = 0; /* djb-rwth: removing redundant code */
    // INCHI✔️❌:
    // INCHI✔️❌:     endCt1 = Ct1->nextCtblPos[k1];
    // INCHI✔️❌:     endCt2 = Ct2->nextCtblPos[k2];
    // INCHI✔️❌:     endAt1 = (int) Ct1->nextAtRank[k1] - 1;
    // INCHI✔️❌:     endAt2 = (int) Ct2->nextAtRank[k2] - 1;
    // INCHI✔️❌:
    // INCHI✔️❌:     maxVert = inchi_min( Ct1->maxVert, Ct2->maxVert );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (bOnlyCommon)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         endCt1 = inchi_min( endCt1, endCt2 );
    // INCHI✔️❌:         endCt1 = endCt2 = inchi_min( endCt1, Ct1->lenCt );
    // INCHI✔️❌:         endAt1 = inchi_min( endAt1, endAt2 ); /* djb-rwth: removing redundant code */
    // INCHI✔️❌:
    // INCHI✔️❌:         if (Ct1->Ctbl[endCt1] == EMPTY_CT || Ct2->Ctbl[endCt1] == EMPTY_CT) /* djb-rwth: redundant conditions removed */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             endCt1 = endCt2 = endCt1 - 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         lenNumH1 =
    // INCHI✔️❌:         lenNumH2 = inchi_min( Ct1->lenNumH, Ct2->lenNumH );
    // INCHI✔️❌:
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         len_iso_sort_key1 =
    // INCHI✔️❌:         len_iso_sort_key2 = inchi_min( Ct1->len_iso_sort_key, Ct1->len_iso_sort_key );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Ct1->Ctbl[endCt1 - 1] == EMPTY_CT)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             endCt1--;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (Ct2->Ctbl[endCt2 - 1] == EMPTY_CT)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             endCt2--;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         lenNumH1 = Ct1->lenNumH;
    // INCHI✔️❌:         lenNumH2 = Ct2->lenNumH;
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:
    // INCHI✔️❌:         len_iso_sort_key1 = Ct1->len_iso_sort_key;
    // INCHI✔️❌:         len_iso_sort_key2 = Ct2->len_iso_sort_key;
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if ((diff = -( endCt1 - endCt2 ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* negative sign reproduces results for NSC=28393 */
    // INCHI✔️❌:         goto done;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (bSplitTautCompare)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         midCt = inchi_min( Ct1->nLenCTAtOnly, Ct2->nLenCTAtOnly );
    // INCHI✔️❌:         if (midCt > endCt1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             midCt = endCt1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         midAt = inchi_min( maxVert, endAt1 ); /* djb-rwth: ignoring LLVM warning: variable used? */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         midCt = endCt1;
    // INCHI✔️❌:         midAt = endAt1; /* djb-rwth: ignoring LLVM warning: variable used? */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:     /************ layer 0: connection table without tautomeric groups ********/
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:     for (i = startCt1; i < midCt && Ct1->Ctbl[i] == Ct2->Ctbl[i]; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (i < midCt)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         diff = (int) Ct1->Ctbl[i] - (int) Ct2->Ctbl[i];
    // INCHI✔️❌:         goto done;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:     /************* layer 1: H atoms without tautomeric H *********************/
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:     nLayer++;
    // INCHI✔️❌:     if (Ct1->NumH && Ct2->NumH)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((diff = -( lenNumH1 - lenNumH2 ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* negative sign reproduces results for NSC=28393 */
    // INCHI✔️❌:             goto done;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (endAt1 < maxVert)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             midNumH = lenNumH1 = endAt1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (bSplitTautCompare)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             midNumH = maxVert;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             midNumH = lenNumH1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = startAt1; i < midNumH && Ct1->NumH[i] == Ct2->NumH[i]; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i < midNumH)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             diff = (int) Ct1->NumH[i] - (int) Ct2->NumH[i];
    // INCHI✔️❌:             goto done;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:     /************** layer 2: tautomeric part of CT and tautomeric H **********/
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:     nLayer++;
    // INCHI✔️❌:     for (i = midCt; i < endCt1 && Ct1->Ctbl[i] == Ct2->Ctbl[i]; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ; /* compare tautomeric groups part of CT */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (i < endCt1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         diff = (int) Ct1->Ctbl[i] - (int) Ct2->Ctbl[i];
    // INCHI✔️❌:         goto done;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (Ct1->NumH && Ct2->NumH)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = midNumH; i < lenNumH1 && Ct1->NumH[i] == Ct2->NumH[i]; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ; /* compare tautomeric H */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i < lenNumH1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             diff = (int) Ct1->NumH[i] - (int) Ct2->NumH[i];
    // INCHI✔️❌:             goto done;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:     /************** layer 3: Fixed H atoms ***********************************/
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:     nLayer++;
    // INCHI✔️❌:     if (Ct1->NumHfixed && Ct2->NumHfixed)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = startAt1; i < endAt1 && Ct1->NumHfixed[i] == Ct2->NumHfixed[i]; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i < endAt1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             diff = (int) Ct1->NumHfixed[i] - (int) Ct2->NumHfixed[i];
    // INCHI✔️❌:             goto done;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:     /************** layer 4: isotopic atoms, H and isotopic taut H ***********/
    // INCHI✔️❌:     /*************************************************************************/
    // INCHI✔️❌:     nLayer++;
    // INCHI✔️❌:     if (Ct1->iso_sort_key && Ct2->iso_sort_key)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((diff = -( len_iso_sort_key1 - len_iso_sort_key2 ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* negative sign reproduces results for NSC=28393 */
    // INCHI✔️❌:             goto done;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         for (i = startAt1; i < endAt1 && Ct1->iso_sort_key[i] == Ct2->iso_sort_key[i]; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i < endAt1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             diff = Ct1->iso_sort_key[i] > Ct2->iso_sort_key[i] ? 1 : -1;
    // INCHI✔️❌:             goto done;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI❌❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI❌❌:
    // INCHI❌❌:
    // INCHI❌❌:     /*************************************************************************/
    // INCHI❌❌:     /************** layer 6: Fixed isotopic H atoms **************************/
    // INCHI❌❌:     /*************************************************************************/
    // INCHI❌❌:     nLayer++;
    // INCHI❌❌:     if (Ct1->iso_sort_key_Hfixed && Ct2->iso_sort_key_Hfixed)
    // INCHI❌❌:     {
    // INCHI❌❌:         for (i = startAt1; i < midAt && Ct1->iso_sort_key_Hfixed[i] == Ct2->iso_sort_key_Hfixed[i]; i++)
    // INCHI❌❌:             ;
    // INCHI❌❌:         if (i < midAt)
    // INCHI❌❌:         {
    // INCHI❌❌:             diff = Ct1->iso_sort_key_Hfixed[i] > Ct2->iso_sort_key_Hfixed[i] ? 1 : -1;
    // INCHI❌❌:             goto done;
    // INCHI❌❌:         }
    // INCHI❌❌:     }
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: done:
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef INCHI_CANON_MIN
    // INCHI✔️❌:     diff = -diff;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     if (diff)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         diff = ( diff > 0 ) ? ( nLayer + 1 ) : -( nLayer + 1 ); /* return the discovered difference layer number >= 1 */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return diff;
    // INCHI✔️❌: }
    // INCHI✔️❌:
    // END INCHI C FUNCTION: CtFullCompare
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CtFullCompare
    // INCHI✔️❌: #define INCHI_CANON_MIN
    // INCHI✔️❌: #define USE_ISO_SORT_KEY_HFIXED  0
    // END INCHI ACTIVE MACRO CONFIGURATION: CtFullCompare

    // INCHI✔️✔️:     k1 = Ct1->lenPos - 1;
    let k1 = Ct1.lenPos.wrapping_sub(1);
    // INCHI✔️✔️:     k2 = Ct2->lenPos - 1;
    let k2 = Ct2.lenPos.wrapping_sub(1);

    // SAFETY: CTableCreate owns each table field in a fixed-size allocation.
    // CtFullCompare neither writes, frees, nor resizes any of them. Views are
    // created in the same order in which the C function first dereferences
    // the corresponding fields.
    let ct_positions1 = unsafe { heap.stable_slice(Ct1.nextCtblPos.as_const())? };
    // INCHI✔️✔️:     endCt1 = Ct1->nextCtblPos[k1];
    let mut end_ct1 = i32::from(
        *ct_positions1
            .get(usize::try_from(k1).map_err(|_| SourceHeapError::PointerOutOfBounds)?)?,
    );
    let ct_positions2 = unsafe { heap.stable_slice(Ct2.nextCtblPos.as_const())? };
    // INCHI✔️✔️:     endCt2 = Ct2->nextCtblPos[k2];
    let mut end_ct2 = i32::from(
        *ct_positions2
            .get(usize::try_from(k2).map_err(|_| SourceHeapError::PointerOutOfBounds)?)?,
    );
    let atom_positions1 = unsafe { heap.stable_slice(Ct1.nextAtRank.as_const())? };
    // INCHI✔️✔️:     endAt1 = (int) Ct1->nextAtRank[k1] - 1;
    let mut end_at1 = i32::from(
        *atom_positions1
            .get(usize::try_from(k1).map_err(|_| SourceHeapError::PointerOutOfBounds)?)?,
    )
    .wrapping_sub(1);
    let atom_positions2 = unsafe { heap.stable_slice(Ct2.nextAtRank.as_const())? };
    // INCHI✔️✔️:     endAt2 = (int) Ct2->nextAtRank[k2] - 1;
    let end_at2 = i32::from(
        *atom_positions2
            .get(usize::try_from(k2).map_err(|_| SourceHeapError::PointerOutOfBounds)?)?,
    )
    .wrapping_sub(1);
    // INCHI✔️✔️:     maxVert = inchi_min( Ct1->maxVert, Ct2->maxVert );
    let max_vert = Ct1.maxVert.min(Ct2.maxVert);

    let table1 = unsafe { heap.stable_slice(Ct1.Ctbl.as_const())? };
    let table2 = unsafe { heap.stable_slice(Ct2.Ctbl.as_const())? };
    let (mut len_h1, len_h2, len_iso1, len_iso2);
    // INCHI✔️✔️:     if (bOnlyCommon)
    if bOnlyCommon != 0 {
        end_ct1 = end_ct1.min(end_ct2);
        end_ct1 = end_ct1.min(Ct1.lenCt);
        end_ct2 = end_ct1;
        end_at1 = end_at1.min(end_at2);
        let end_index =
            usize::try_from(end_ct1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        // INCHI✔️✔️:         if (Ct1->Ctbl[endCt1] == EMPTY_CT || Ct2->Ctbl[endCt1] == EMPTY_CT)
        if *table1.get(end_index)? == EMPTY_CT as AT_RANK
            || *table2.get(end_index)? == EMPTY_CT as AT_RANK
        {
            end_ct1 = end_ct1.wrapping_sub(1);
            end_ct2 = end_ct1;
        }
        len_h1 = Ct1.lenNumH.min(Ct2.lenNumH);
        len_h2 = len_h1;
        len_iso1 = Ct1.len_iso_sort_key.min(Ct1.len_iso_sort_key);
        len_iso2 = len_iso1;
    } else {
        let end_index1 = usize::try_from(end_ct1.wrapping_sub(1))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        // INCHI✔️✔️:         if (Ct1->Ctbl[endCt1 - 1] == EMPTY_CT)
        if *table1.get(end_index1)? == EMPTY_CT as AT_RANK {
            end_ct1 = end_ct1.wrapping_sub(1);
        }
        let end_index2 = usize::try_from(end_ct2.wrapping_sub(1))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        // INCHI✔️✔️:         if (Ct2->Ctbl[endCt2 - 1] == EMPTY_CT)
        if *table2.get(end_index2)? == EMPTY_CT as AT_RANK {
            end_ct2 = end_ct2.wrapping_sub(1);
        }
        len_h1 = Ct1.lenNumH;
        len_h2 = Ct2.lenNumH;
        len_iso1 = Ct1.len_iso_sort_key;
        len_iso2 = Ct2.len_iso_sort_key;
    }
    let mut i = 0_i32;
    let mut layer = 0_i32;
    let mut mid_h = 0_i32;
    let mut diff = 'compare: {
        // INCHI✔️✔️:     if ((diff = -( endCt1 - endCt2 )))
        let mut value = end_ct2.wrapping_sub(end_ct1);
        if value != 0 {
            break 'compare value;
        }
        let (mid_ct, _mid_at) = if bSplitTautCompare != 0 {
            (
                Ct1.nLenCTAtOnly.min(Ct2.nLenCTAtOnly).min(end_ct1),
                max_vert.min(end_at1),
            )
        } else {
            (end_ct1, end_at1)
        };
        // INCHI✔️✔️:     for (i = startCt1; i < midCt && Ct1->Ctbl[i] == Ct2->Ctbl[i]; i++)
        i = 0;
        while i < mid_ct {
            let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            if *table1.get(index)? != *table2.get(index)? {
                break;
            }
            i = i.wrapping_add(1);
        }
        // INCHI✔️✔️:     if (i < midCt)
        if i < mid_ct {
            let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            // INCHI✔️✔️:         diff = (int) Ct1->Ctbl[i] - (int) Ct2->Ctbl[i];
            value = i32::from(*table1.get(index)?).wrapping_sub(i32::from(*table2.get(index)?));
            break 'compare value;
        }
        layer = layer.wrapping_add(1);
        let mut hydrogen_views = None;
        // INCHI✔️✔️:     if (Ct1->NumH && Ct2->NumH)
        if !Ct1.NumH.is_null() && !Ct2.NumH.is_null() {
            // INCHI✔️✔️:         if ((diff = -( lenNumH1 - lenNumH2 )))
            value = len_h2.wrapping_sub(len_h1);
            if value != 0 {
                break 'compare value;
            }
            if end_at1 < max_vert {
                mid_h = end_at1;
                len_h1 = end_at1;
            } else if bSplitTautCompare != 0 {
                mid_h = max_vert;
            } else {
                mid_h = len_h1;
            }
            let hydrogens1 = unsafe { heap.stable_slice(Ct1.NumH.as_const())? };
            let hydrogens2 = unsafe { heap.stable_slice(Ct2.NumH.as_const())? };
            // INCHI✔️✔️:         for (i = startAt1; i < midNumH && Ct1->NumH[i] == Ct2->NumH[i]; i++)
            i = 0;
            while i < mid_h {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if *hydrogens1.get(index)? != *hydrogens2.get(index)? {
                    break;
                }
                i = i.wrapping_add(1);
            }
            if i < mid_h {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                value = i32::from(*hydrogens1.get(index)?)
                    .wrapping_sub(i32::from(*hydrogens2.get(index)?));
                break 'compare value;
            }
            hydrogen_views = Some((hydrogens1, hydrogens2));
        }
        layer = layer.wrapping_add(1);
        // INCHI✔️✔️:     for (i = midCt; i < endCt1 && Ct1->Ctbl[i] == Ct2->Ctbl[i]; i++)
        i = mid_ct;
        while i < end_ct1 {
            let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            if *table1.get(index)? != *table2.get(index)? {
                break;
            }
            i = i.wrapping_add(1);
        }
        if i < end_ct1 {
            let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            value = i32::from(*table1.get(index)?).wrapping_sub(i32::from(*table2.get(index)?));
            break 'compare value;
        }
        if let Some((hydrogens1, hydrogens2)) = &hydrogen_views {
            // INCHI✔️✔️:         for (i = midNumH; i < lenNumH1 && Ct1->NumH[i] == Ct2->NumH[i]; i++)
            i = mid_h;
            while i < len_h1 {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if *hydrogens1.get(index)? != *hydrogens2.get(index)? {
                    break;
                }
                i = i.wrapping_add(1);
            }
            if i < len_h1 {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                value = i32::from(*hydrogens1.get(index)?)
                    .wrapping_sub(i32::from(*hydrogens2.get(index)?));
                break 'compare value;
            }
        }
        layer = layer.wrapping_add(1);
        // INCHI✔️✔️:     if (Ct1->NumHfixed && Ct2->NumHfixed)
        if !Ct1.NumHfixed.is_null() && !Ct2.NumHfixed.is_null() {
            let fixed_hydrogens1 = unsafe { heap.stable_slice(Ct1.NumHfixed.as_const())? };
            let fixed_hydrogens2 = unsafe { heap.stable_slice(Ct2.NumHfixed.as_const())? };
            // INCHI✔️✔️:         for (i = startAt1; i < endAt1 && Ct1->NumHfixed[i] == Ct2->NumHfixed[i]; i++)
            i = 0;
            while i < end_at1 {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if *fixed_hydrogens1.get(index)? != *fixed_hydrogens2.get(index)? {
                    break;
                }
                i = i.wrapping_add(1);
            }
            if i < end_at1 {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                value = i32::from(*fixed_hydrogens1.get(index)?)
                    .wrapping_sub(i32::from(*fixed_hydrogens2.get(index)?));
                break 'compare value;
            }
        }
        layer = layer.wrapping_add(1);
        // INCHI✔️✔️:     if (Ct1->iso_sort_key && Ct2->iso_sort_key)
        if !Ct1.iso_sort_key.is_null() && !Ct2.iso_sort_key.is_null() {
            // INCHI✔️✔️:         if ((diff = -( len_iso_sort_key1 - len_iso_sort_key2 )))
            value = len_iso2.wrapping_sub(len_iso1);
            if value != 0 {
                break 'compare value;
            }
            let isotope_keys1 = unsafe { heap.stable_slice(Ct1.iso_sort_key.as_const())? };
            let isotope_keys2 = unsafe { heap.stable_slice(Ct2.iso_sort_key.as_const())? };
            // INCHI✔️✔️:         for (i = startAt1; i < endAt1 && Ct1->iso_sort_key[i] == Ct2->iso_sort_key[i]; i++)
            i = 0;
            while i < end_at1 {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if *isotope_keys1.get(index)? != *isotope_keys2.get(index)? {
                    break;
                }
                i = i.wrapping_add(1);
            }
            if i < end_at1 {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                value = if *isotope_keys1.get(index)? > *isotope_keys2.get(index)? {
                    1
                } else {
                    -1
                };
                break 'compare value;
            }
        }
        0
    };
    diff = diff.wrapping_neg();
    if diff != 0 {
        diff = if diff > 0 {
            layer.wrapping_add(1)
        } else {
            layer.wrapping_add(1).wrapping_neg()
        };
    }
    Ok(diff)
}

#[allow(non_snake_case)]
pub(crate) fn CtFullCompareLayers(kLeastForLayer: &KLeastLayers) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2863 CtFullCompareLayers
    // INCHI✔️✔️: int CtFullCompareLayers( kLeast *kLeastForLayer )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int iLayer;
    // INCHI✔️✔️:     /* check for the rejection condition: Lambda > zb_rho_fix */
    // INCHI✔️✔️:     for (iLayer = 0; iLayer < MAX_LAYERS; iLayer++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (kLeastForLayer[iLayer].k)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return ( kLeastForLayer[iLayer].k > 0 ) ? ( iLayer + 1 ) : -( iLayer + 1 );
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: CtFullCompareLayers

    for (layer, value) in kLeastForLayer.iter().enumerate() {
        if value.k != 0 {
            let layer = i32::try_from(layer)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?
                .wrapping_add(1);
            return Ok(if value.k > 0 {
                layer
            } else {
                layer.wrapping_neg()
            });
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CtCompareLayersGetFirstDiff(
    kLeast_rho: Option<&KLeastLayers>,
    nOneAdditionalLayer: i32,
    L_rho: &mut i32,
    I_rho: &mut i32,
    k_rho: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2880 CtCompareLayersGetFirstDiff
    // INCHI✔️✔️: int CtCompareLayersGetFirstDiff( kLeast *kLeast_rho,
    // INCHI✔️✔️:                                  int nOneAdditionalLayer,
    // INCHI✔️✔️:                                  int *L_rho,
    // INCHI✔️✔️:                                  int *I_rho,
    // INCHI✔️✔️:                                  int *k_rho )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int iLayer;
    // INCHI✔️✔️:     if (kLeast_rho)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         for (iLayer = 0; iLayer < MAX_LAYERS; iLayer++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (kLeast_rho[iLayer].k)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 *L_rho = iLayer;
    // INCHI✔️✔️:                 *I_rho = kLeast_rho[iLayer].i;
    // INCHI✔️✔️:                 *k_rho = kLeast_rho[iLayer].k;
    // INCHI✔️✔️:                 break;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (iLayer == MAX_LAYERS)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (nOneAdditionalLayer)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 *L_rho = nOneAdditionalLayer;  /* ??? subtract 1 ??? */
    // INCHI✔️✔️:                 *I_rho = -1;
    // INCHI✔️✔️:                 *k_rho = 0;
    // INCHI✔️✔️:                 return 0; /* difference may be in the first additional layer */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 *L_rho = INCHI_CANON_INFINITY;
    // INCHI✔️✔️:                 *I_rho = -1;
    // INCHI✔️✔️:                 *k_rho = 0;
    // INCHI✔️✔️:                 return 0; /* no difference found */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return 1; /* difference in a real layer */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1; /* no input, should not happen */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: CtCompareLayersGetFirstDiff

    let Some(values) = kLeast_rho else {
        return Ok(-1);
    };
    for (layer, value) in values.iter().enumerate() {
        if value.k != 0 {
            *L_rho = i32::try_from(layer).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            *I_rho = value.i;
            *k_rho = value.k;
            return Ok(1);
        }
    }
    *L_rho = if nOneAdditionalLayer != 0 {
        nOneAdditionalLayer
    } else {
        INCHI_CANON_INFINITY as i32
    };
    *I_rho = -1;
    *k_rho = 0;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CtPartCompareLayers(
    kLeast_rho: Option<&KLeastLayers>,
    L_rho_fix_prev: i32,
    nOneAdditionalLayer: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2929 CtPartCompareLayers
    // INCHI✔️✔️: int CtPartCompareLayers( kLeast *kLeast_rho,
    // INCHI✔️✔️:                          int L_rho_fix_prev,
    // INCHI✔️✔️:                          int nOneAdditionalLayer )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int L_rho, I_rho, k_rho;
    // INCHI✔️✔️:     if (0 < CtCompareLayersGetFirstDiff( kLeast_rho, nOneAdditionalLayer, &L_rho, &I_rho, &k_rho ) &&
    // INCHI✔️✔️:          /* differences has been found in a real layer or all real layers are identical */
    // INCHI✔️✔️:          L_rho <= L_rho_fix_prev)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* in this layer pzb_rho == pzb_rho_fix or in the previous real layer */
    // INCHI✔️✔️:         return k_rho > 0 ? ( L_rho + 1 ) : -( L_rho + 1 );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: CtPartCompareLayers

    let (mut layer, mut item, mut k) = (0_i32, 0_i32, 0_i32);
    if CtCompareLayersGetFirstDiff(
        kLeast_rho,
        nOneAdditionalLayer,
        &mut layer,
        &mut item,
        &mut k,
    )? > 0
        && layer <= L_rho_fix_prev
    {
        let result = layer.wrapping_add(1);
        return Ok(if k > 0 { result } else { result.wrapping_neg() });
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn UpdateCompareLayers(
    kLeastForLayer: Option<&mut KLeastLayers>,
    hzz: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2947 UpdateCompareLayers
    // INCHI✔️✔️: void UpdateCompareLayers( kLeast kLeastForLayer[], int hzz )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i;
    // INCHI✔️✔️:     if (kLeastForLayer)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         for (i = 0; i < MAX_LAYERS; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (abs( kLeastForLayer[i].k ) >= hzz)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 kLeastForLayer[i].k = 0;
    // INCHI✔️✔️:                 kLeastForLayer[i].i = 0;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: UpdateCompareLayers

    let Some(values) = kLeastForLayer else {
        return Ok(());
    };
    for value in values.iter_mut() {
        if value.k.abs() >= hzz {
            value.k = 0;
            value.i = 0;
        }
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn TranspositionGetMcrAndFixSetAndUnorderedPartition(
    heap: &mut SourceHeap,
    pCG: &CANON_GLOBALS,
    gamma: &mut Transposition,
    McrSet: &NodeSet,
    FixSet: &NodeSet,
    n: i32,
    l: i32,
    p: &mut UnorderedPartition,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3125 TranspositionGetMcrAndFixSetAndUnorderedPartition
    // INCHI✔️❌: void TranspositionGetMcrAndFixSetAndUnorderedPartition( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                                                         Transposition *gamma,
    // INCHI✔️❌:                                                         NodeSet *McrSet,
    // INCHI✔️❌:                                                         NodeSet *FixSet,
    // INCHI✔️❌:                                                         int n,
    // INCHI✔️❌:                                                         int l,
    // INCHI✔️❌:                                                         UnorderedPartition *p )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, k, mcr;
    // INCHI✔️❌:     AT_RANK next;
    // INCHI✔️❌:     bitWord *McrBits = McrSet->bitword[l - 1];
    // INCHI✔️❌:     bitWord *FixBits = FixSet->bitword[l - 1];
    // INCHI✔️❌:     int     len = McrSet->len_set * sizeof( bitWord );
    // INCHI✔️❌:
    // INCHI✔️❌:     memset( McrBits, 0, len ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( FixBits, 0, len ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     for (i = 0; i < n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         p->equ2[i] = INCHI_CANON_INFINITY; /* for debug only */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         j = (int) ( next = gamma->nAtNumb[i] );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (j == i)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             FixBits[i / pCG->m_num_bit] |= pCG->m_bBit[i % pCG->m_num_bit];
    // INCHI✔️❌:             McrBits[i / pCG->m_num_bit] |= pCG->m_bBit[i % pCG->m_num_bit];
    // INCHI✔️❌:             /* p->next[i] = INCHI_CANON_INFINITY; */ /* no link to same orbit points */
    // INCHI✔️❌:             p->equ2[i] = next;  /* fixed point */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (!( rank_mark_bit & next ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             gamma->nAtNumb[i] |= rank_mark_bit;
    // INCHI✔️❌:             mcr = inchi_min( j, i );
    // INCHI✔️❌:             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:             /* mark all nodes in the cycle to ignore later; find mcr */
    // INCHI✔️❌:             while (!( rank_mark_bit & ( next = gamma->nAtNumb[j] ) ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 gamma->nAtNumb[j] |= rank_mark_bit;
    // INCHI✔️❌:                 j = (int) next;
    // INCHI✔️❌:                 if (mcr > j)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     mcr = j;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* djb-rwth: removing redundant code */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             McrBits[mcr / pCG->m_num_bit] |= pCG->m_bBit[mcr % pCG->m_num_bit]; /* save mcr */
    // INCHI✔️❌:             /* fill out the unordered partition, the mcr first, other in the cycle after that */
    // INCHI✔️❌:             p->equ2[mcr] = mcr;
    // INCHI✔️❌:             for (k = mcr; mcr != ( j = (int) ( rank_mask_bit & gamma->nAtNumb[k] ) ); k = j)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 p->equ2[j] = mcr;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI❌❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI❌❌:     /* for debug only */
    // INCHI❌❌:     for (i = 0; i < n; i++)
    // INCHI❌❌:     {
    // INCHI❌❌:         if (p->equ2[i] >= n)
    // INCHI❌❌:         {
    // INCHI❌❌:             int stop = 1;
    // INCHI❌❌:         }
    // INCHI❌❌:     }
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     /* remove the marks */
    // INCHI✔️❌:     for (i = 0; i < n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         gamma->nAtNumb[i] &= rank_mask_bit;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: TranspositionGetMcrAndFixSetAndUnorderedPartition
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: TranspositionGetMcrAndFixSetAndUnorderedPartition
    // INCHI✔️❌: #define bRELEASE_VERSION  1
    // END INCHI ACTIVE MACRO CONFIGURATION: TranspositionGetMcrAndFixSetAndUnorderedPartition

    let mcr_bits = pointer_array_get(heap, McrSet.bitword, l.wrapping_sub(1))?;
    let fix_bits = pointer_array_get(heap, FixSet.bitword, l.wrapping_sub(1))?;
    let set_len =
        usize::try_from(McrSet.len_set).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    // SAFETY: CanonGraph creates the McrSet and FixSet rows as distinct,
    // fixed-size allocations and keeps both live for this complete call.
    let mut mcr_words = unsafe { heap.stable_slice_mut(mcr_bits)? };
    // INCHI✔️✔️:     memset( McrBits, 0, len );
    mcr_words.prefix_mut(set_len)?.fill(0);
    let mut fix_words = unsafe { heap.stable_slice_mut(fix_bits)? };
    // INCHI✔️✔️:     memset( FixBits, 0, len );
    fix_words.prefix_mut(set_len)?.fill(0);

    // The C loops do not dereference gamma, equ2, or the mask table when n is
    // non-positive; the two memset operations above still execute.
    if n <= 0 {
        return Ok(());
    }

    // SAFETY: TranspositionCreate, UnorderedPartitionCreate, and CanonGraph
    // own these three arrays independently from each other and from both set
    // rows. No allocation is resized or freed while these views are live.
    let mut permutation = unsafe { heap.stable_slice_mut(gamma.nAtNumb)? };
    let mut equivalences = unsafe { heap.stable_slice_mut(p.equ2)? };
    let masks = unsafe { heap.stable_slice(pCG.m_bBit.as_const())? };

    // INCHI✔️✔️:     for (i = 0; i < n; i++)
    let mut i = 0_i32;
    while i < n {
        let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        // INCHI✔️✔️:         p->equ2[i] = INCHI_CANON_INFINITY;
        *equivalences.get_mut(index)? = INCHI_CANON_INFINITY as AT_NUMB;
        i = i.wrapping_add(1);
    }

    // INCHI✔️✔️:     for (i = 0; i < n; i++)
    i = 0;
    while i < n {
        let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        // INCHI✔️✔️:         j = (int) ( next = gamma->nAtNumb[i] );
        let next = *permutation.get(index)?;
        let mut j = i32::from(next);
        // INCHI✔️✔️:         if (j == i)
        if j == i {
            let word = i.wrapping_div(pCG.m_num_bit);
            let bit_index = i.wrapping_rem(pCG.m_num_bit);
            let word_index =
                usize::try_from(word).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let mask_index =
                usize::try_from(bit_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let bit = *masks.get(mask_index)?;
            // INCHI✔️✔️:             FixBits[i / pCG->m_num_bit] |= pCG->m_bBit[i % pCG->m_num_bit];
            *fix_words.get_mut(word_index)? |= bit;
            // INCHI✔️✔️:             McrBits[i / pCG->m_num_bit] |= pCG->m_bBit[i % pCG->m_num_bit];
            *mcr_words.get_mut(word_index)? |= bit;
            // INCHI✔️✔️:             p->equ2[i] = next;
            *equivalences.get_mut(index)? = next;
        // INCHI✔️✔️:         else if (!( rank_mark_bit & next ))
        } else if rank_mark_bit() & next == 0 {
            // INCHI✔️✔️:             gamma->nAtNumb[i] |= rank_mark_bit;
            *permutation.get_mut(index)? = next | rank_mark_bit();
            // INCHI✔️✔️:             mcr = inchi_min( j, i );
            let mut mcr = j.min(i);
            // INCHI✔️✔️:             while (!( rank_mark_bit & ( next = gamma->nAtNumb[j] ) ))
            loop {
                let cycle_index =
                    usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let cycle_next = *permutation.get(cycle_index)?;
                if rank_mark_bit() & cycle_next != 0 {
                    break;
                }
                // INCHI✔️✔️:                 gamma->nAtNumb[j] |= rank_mark_bit;
                *permutation.get_mut(cycle_index)? = cycle_next | rank_mark_bit();
                // INCHI✔️✔️:                 j = (int) next;
                j = i32::from(cycle_next);
                // INCHI✔️✔️:                 if (mcr > j)
                if mcr > j {
                    // INCHI✔️✔️:                     mcr = j;
                    mcr = j;
                }
            }
            let word = mcr.wrapping_div(pCG.m_num_bit);
            let bit_index = mcr.wrapping_rem(pCG.m_num_bit);
            let word_index =
                usize::try_from(word).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let mask_index =
                usize::try_from(bit_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let bit = *masks.get(mask_index)?;
            // INCHI✔️✔️:             McrBits[mcr / pCG->m_num_bit] |= pCG->m_bBit[mcr % pCG->m_num_bit];
            *mcr_words.get_mut(word_index)? |= bit;
            let mcr_index =
                usize::try_from(mcr).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            // INCHI✔️✔️:             p->equ2[mcr] = mcr;
            *equivalences.get_mut(mcr_index)? = mcr as AT_NUMB;
            // INCHI✔️✔️:             for (k = mcr; mcr != (j = (int)(rank_mask_bit & gamma->nAtNumb[k])); k = j)
            let mut k = mcr;
            loop {
                let cycle_index =
                    usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                j = i32::from(rank_mask_bit() & *permutation.get(cycle_index)?);
                if mcr == j {
                    break;
                }
                // INCHI✔️✔️:                 p->equ2[j] = mcr;
                *equivalences.get_mut(
                    usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )? = mcr as AT_NUMB;
                k = j;
            }
        }
        i = i.wrapping_add(1);
    }

    // INCHI✔️✔️:     for (i = 0; i < n; i++)
    i = 0;
    while i < n {
        let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        // INCHI✔️✔️:         gamma->nAtNumb[i] &= rank_mask_bit;
        *permutation.get_mut(index)? &= rank_mask_bit();
        i = i.wrapping_add(1);
    }
    Ok(())
}

#[allow(non_snake_case)]
#[allow(non_snake_case)]
pub(crate) fn GetOneAdditionalLayer(
    pCD: Option<&CANON_DATA>,
    pzb_rho_fix: Option<&ConTable>,
) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3566 GetOneAdditionalLayer
    // INCHI✔️✔️: int GetOneAdditionalLayer( CANON_DATA *pCD, ConTable *pzb_rho_fix )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int nLastLayer = -1, nNumLast = 0, nLayer = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (!pCD || !pzb_rho_fix)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     nLayer++; /* 1 */
    // INCHI✔️✔️:     if (pCD->NumH && !pzb_rho_fix->NumH)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nLastLayer = nLayer;
    // INCHI✔️✔️:         nNumLast++;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     nLayer++; /* 2 */
    // INCHI✔️✔️:     if (pCD->nLenCTAtOnly < pCD->nLenLinearCT && pzb_rho_fix->nLenCTAtOnly == pzb_rho_fix->lenCt)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nLastLayer = nLayer;
    // INCHI✔️✔️:         nNumLast++;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     nLayer++; /* 3 */
    // INCHI✔️✔️:     if (pCD->NumHfixed && !pzb_rho_fix->NumHfixed)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nLastLayer = nLayer;
    // INCHI✔️✔️:         nNumLast++;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     nLayer++; /* 4 */
    // INCHI✔️✔️:     if (pCD->iso_sort_key && !pzb_rho_fix->iso_sort_key)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nLastLayer = nLayer;
    // INCHI✔️✔️:         nNumLast++;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*
    // INCHI✔️✔️:     nLayer ++; // 5
    // INCHI✔️✔️:     if ( pCD->nLenCTAtOnly < pCD->nLenLinearCT && pCD->iso_sort_key &&
    // INCHI✔️✔️:         (pzb_rho_fix->nLenCTAtOnly == pzb_rho_fix->lenCt || !pzb_rho_fix->iso_sort_key ) )
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:         nLastLayer = nLayer;
    // INCHI✔️✔️:         nNumLast ++;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:
    // INCHI❌❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI❌❌:
    // INCHI❌❌:     nLayer++; /* 6 */
    // INCHI❌❌:     if (pCD->iso_sort_key_Hfixed && !pzb_rho_fix->iso_sort_key_Hfixed)
    // INCHI❌❌:     {
    // INCHI❌❌:         nLastLayer = nLayer;
    // INCHI❌❌:         nNumLast++;
    // INCHI❌❌:     }
    // INCHI❌❌: #endif
    // INCHI✔️✔️:     if (1 == nNumLast)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return nLastLayer;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: GetOneAdditionalLayer
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetOneAdditionalLayer
    // INCHI✔️✔️: #define USE_ISO_SORT_KEY_HFIXED  0
    // END INCHI ACTIVE MACRO CONFIGURATION: GetOneAdditionalLayer

    let (Some(data), Some(fixed)) = (pCD, pzb_rho_fix) else {
        return 0;
    };
    let mut last_layer = -1_i32;
    let mut number_last = 0_i32;
    let mut layer = 0_i32;

    layer = layer.wrapping_add(1);
    if !data.NumH.is_null() && fixed.NumH.is_null() {
        last_layer = layer;
        number_last = number_last.wrapping_add(1);
    }
    layer = layer.wrapping_add(1);
    if data.nLenCTAtOnly < data.nLenLinearCT && fixed.nLenCTAtOnly == fixed.lenCt {
        last_layer = layer;
        number_last = number_last.wrapping_add(1);
    }
    layer = layer.wrapping_add(1);
    if !data.NumHfixed.is_null() && fixed.NumHfixed.is_null() {
        last_layer = layer;
        number_last = number_last.wrapping_add(1);
    }
    layer = layer.wrapping_add(1);
    if !data.iso_sort_key.is_null() && fixed.iso_sort_key.is_null() {
        last_layer = layer;
        number_last = number_last.wrapping_add(1);
    }
    if number_last == 1 { last_layer } else { 0 }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum CanonGraphLabel {
    InitialDescent,
    L2,
    L6,
    L7,
    L8,
    L9,
    L10,
    L11,
    L12,
    L13,
    L14,
    L15,
    L16,
    L17,
    ExitFunction,
}

fn partition_at(pi: &[Partition], k: i32) -> Result<&Partition, SourceHeapError> {
    pi.get(usize::try_from(k.wrapping_sub(1)).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn CanonGraph(
    heap: &mut SourceHeap,
    ic: &mut INCHI_CLOCK,
    pCG: &mut CANON_GLOBALS,
    n: i32,
    n_tg: i32,
    n_max: i32,
    bDigraph: i32,
    G: SourceMutPointer<NEIGH_LIST>,
    pi: &mut [Partition],
    nSymmRank: SourceMutPointer<AT_RANK>,
    nCanonRank: SourceMutPointer<AT_RANK>,
    nAtomNumberCanon: SourceMutPointer<AT_NUMB>,
    pCD: &mut CANON_DATA,
    pCC: &mut CANON_COUNTS,
    pzb_rho_inp: SourceMutPointer<ConTable>,
    mut pp_zb_rho_out: Option<&mut SourceMutPointer<ConTable>>,
    LargeMolecules: i32,
    user_action: Option<fn() -> u32>,
    console_quit: Option<fn() -> i32>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3662 CanonGraph
    // INCHI✔️❌: int CanonGraph( INCHI_CLOCK *ic,
    // INCHI✔️❌:                 CANON_GLOBALS *pCG,
    // INCHI✔️❌:                 int n,
    // INCHI✔️❌:                 int n_tg,
    // INCHI✔️❌:                 int n_max,
    // INCHI✔️❌:                 int bDigraph,
    // INCHI✔️❌:                 Graph *G,
    // INCHI✔️❌:                 Partition pi[],
    // INCHI✔️❌:                 AT_RANK *nSymmRank,
    // INCHI✔️❌:                 AT_RANK *nCanonRank,
    // INCHI✔️❌:                 AT_NUMB *nAtomNumberCanon,
    // INCHI✔️❌:                 CANON_DATA *pCD,
    // INCHI✔️❌:                 CANON_COUNTS *pCC,
    // INCHI✔️❌:                 ConTable **pp_zb_rho_inp,
    // INCHI✔️❌:                 ConTable **pp_zb_rho_out,
    // INCHI✔️❌:                 int LargeMolecules )
    // INCHI✔️❌: {
    // INCHI✔️❌:
    // INCHI✔️❌:     /* bDigraph != 0
    // INCHI✔️❌:        means consider edges from atoms to t-groups
    // INCHI✔️❌:        as directed, that is, do not include
    // INCHI✔️❌:        t-group ranks in comparing neighbors
    // INCHI✔️❌:        when refining partition
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Always set
    // INCHI✔️❌:        lab = true
    // INCHI✔️❌:        dig = true
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* in the comments:
    // INCHI✔️❌:         m = |zeta|
    // INCHI✔️❌:         r = |rho|
    // INCHI✔️❌:
    // INCHI✔️❌:         m < n or r < n in case pi[k] in P (i.e. satisfies Lemma 2.25)
    // INCHI✔️❌:
    // INCHI✔️❌:         Just after passing point B:
    // INCHI✔️❌:         ===========================
    // INCHI✔️❌:         K = k-1
    // INCHI✔️❌:         wi = v[i], i = 1..K
    // INCHI✔️❌:         Gamma(0) = Gamma = Aut(G)pi
    // INCHI✔️❌:         Gamma(i) = Gamma(w1,w2,...,wi) pointwise stabilizer for i=1..K
    // INCHI✔️❌:         zeta is a terminal node =>
    // INCHI✔️❌:           the coarsest equitable partition that fixes w1,...,wK is discrete =>
    // INCHI✔️❌:             Gamma(K)=1
    // INCHI✔️❌:         At point A only:
    // INCHI✔️❌:             index = |Gamma(k-1)|/|Gamma(k)|
    // INCHI✔️❌:         At points A and B:
    // INCHI✔️❌:             size  = |Gamma(k-1)|
    // INCHI✔️❌:             theta = theta(Gamma(k-1));
    // INCHI✔️❌:             Gamma(k-1) = <Y>, where Y is the set of all automprhisms output up
    // INCHI✔️❌:               to the present stage (in Step 10 = L10 )
    // INCHI✔️❌:             |Y| <= n - |theta|
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     AT_RANK    *pCt = pCD->LinearCT;
    // INCHI✔️❌:     /*int         nMaxLenCt        =  pCD->nMaxLenLinearCT;*/
    // INCHI✔️❌:     int        *nLenCt = &pCD->nLenLinearCT;
    // INCHI✔️❌:     CANON_DATA *pCD1 = pCD;
    // INCHI✔️❌:
    // INCHI✔️❌:     int i, k, k2, index, l, ok, ret = 0, res;
    // INCHI✔️❌:     int t_Lemma;   /* hh: if pi[k-1] satisfies Lemma 2-25 then
    // INCHI✔️❌:                           t_Lemma = min{i| i=1..k && pi[i-1] satisfies Lemma 2-25}*/
    // INCHI✔️❌:                    /* otherwise t_Lemma = k --> here this is always the case */
    // INCHI✔️❌:     int t_eq_zeta; /* ht: min{i|i=1..m && all terminal modes descended from or equal
    // INCHI✔️❌:                       to zeta(i) have been shown to be equivalent}. */
    // INCHI✔️❌:     int h_zeta;    /* h: the longest common ancestor of zeta and nu is nu(h_zeta) */
    // INCHI✔️❌:     int h_rho;     /* hb: the longest common ancestor of rho and nu is nu(h_rho) */
    // INCHI✔️❌:     int hz_rho;    /* hzb: max{i|i=1..min(k,r) && Lambda(G,pi,nu(i)) == Lambda(G,pi,rho(i))} */
    // INCHI✔️❌:     int hz_zeta;   /* hzf: max{i|i=1..min(k,m) && Lambda(G,pi,nu(i)) == Lambda(G,pi,zeta(i))} */
    // INCHI✔️❌:     int qzb_rho;   /* Ct(Lambda[k]) - Ct(rho[k]) */
    // INCHI✔️❌:     double size;   /* |Aut(G)| */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     int  nNumLayers = ( NULL != pCD->NumH ) + ( NULL != pCD->NumHfixed ) +
    // INCHI✔️❌:                       /* (bDigraph && pCD->nLenLinearCT > pCD->nLenCTAtOnly)*/ /* ??? tautomeric */
    // INCHI✔️❌:                       ( NULL != pCD->iso_sort_key )
    // INCHI❌❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI❌❌:                     + ( NULL != pCD->iso_sort_key_Hfixed )
    // INCHI❌❌: #endif
    // INCHI✔️❌:                     ;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     int  dig = ( bDigraph || nNumLayers );
    // INCHI✔️❌:     int  bSplitTautCompare = ( bDigraph || nNumLayers ); /* compare taut. H and tgroups connections after H */
    // INCHI✔️❌:                         /* digraph: 1=>do not use Lemma 2.25, 0 => use */
    // INCHI✔️❌:     int  lab = 1;       /* label: 1=>find canonical numbering;
    // INCHI✔️❌:                             0=>do not find canonical numbering, do not use rho */
    // INCHI✔️❌:     int  r;             /* |rho| */
    // INCHI✔️❌:     int  bZetaEqRho = lab;
    // INCHI✔️❌:     int  bZetaIsomorph;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*const int L = MAX_SET_SIZE;*/
    // INCHI✔️❌:     int L_curr_max_set_size = MAX_SET_SIZE;
    // INCHI✔️❌:
    // INCHI✔️❌:     UnorderedPartition theta, theta_from_gamma;
    // INCHI✔️❌:
    // INCHI✔️❌:     Cell *W;  /* W[i] is the first non-trivial cell of pi[i+1] */
    // INCHI✔️❌:     Node *v;  /* v[i] is in W[i] to create T(G,pi,nu[i+1]) */
    // INCHI✔️❌:     Node tvc, tvh;
    // INCHI✔️❌:     S_CHAR *e, *qzb = NULL;   /* qzb = NULL always */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* current node CT */
    // INCHI✔️❌:     ConTable Lambda;
    // INCHI✔️❌:     /* first leaf CT */
    // INCHI✔️❌:     ConTable zf_zeta; /* Ct for zeta,  the first discovered terminal node */
    // INCHI✔️❌:     /* best leaf/node CT: find the greatest pzb_rho possibly subject to pzb_rho[k] <= pzb_rho_fix[k] condition */
    // INCHI✔️❌:     ConTable *pzb_rho = NULL;  /* Ct for rho,   the best discovered terminal node */
    // INCHI✔️❌:     /* fixed input CT: for all k pzb_rho[k] <= pzb_rho_fix[k]; at the end pzb_rho == pzb_rho_fix */
    // INCHI✔️❌:     ConTable *pzb_rho_fix = ( pp_zb_rho_inp && *pp_zb_rho_inp ) ? *pp_zb_rho_inp : NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     NodeSet Omega; /* MAX_SET_SIZE */
    // INCHI✔️❌:     NodeSet Phi;   /* MAX_SET_SIZE */
    // INCHI✔️❌:     NodeSet cur_nodes;   /* 1 each */
    // INCHI✔️❌:
    // INCHI✔️❌:     Transposition gamma;
    // INCHI✔️❌:
    // INCHI✔️❌:     Partition zeta;      /* the first discovered terminal node */
    // INCHI✔️❌:     Partition rho;       /* the best discovered terminal node */
    // INCHI✔️❌:
    // INCHI✔️❌:     int nNumFoundGenerators = 0;
    // INCHI✔️❌:     int qzb_rho_fix = 0;
    // INCHI✔️❌:     int hzb_rho_fix = 0;
    // INCHI✔️❌:     int bRhoIsDiscrete = 1;
    // INCHI✔️❌:     kLeast kLeast_rho[MAX_LAYERS];
    // INCHI✔️❌:     kLeast kLeast_rho_fix[MAX_LAYERS];
    // INCHI✔️❌:     int nOneAdditionalLayer;
    // INCHI✔️❌:     int pzb_rho_fix_reached = 0;
    // INCHI✔️❌:     int L_rho_fix_prev = 0, I_rho_fix_prev = -1, k_rho_fix_prev = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!LargeMolecules)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         L_curr_max_set_size = NORMALLY_ALLOWED_MAX_SET_SIZE;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Note: Layered comparison should be consistent, especially in layer numbers.
    // INCHI✔️❌:              Layered comparison is implemented in:
    // INCHI✔️❌:                      CtFullCompare()
    // INCHI✔️❌:                      CtPartCompare()
    // INCHI✔️❌:                      GetOneAdditionalLayer()
    // INCHI✔️❌:
    // INCHI✔️❌:              The partial comparison results in kLeast[] are used in
    // INCHI✔️❌:                      CtFullCompareLayers()
    // INCHI✔️❌:                      CtPartCompareLayers()
    // INCHI✔️❌:                      CtCompareLayersGetFirstDiff()
    // INCHI✔️❌:                      UpdateCompareLayers()
    // INCHI✔️❌:      */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     nOneAdditionalLayer = GetOneAdditionalLayer( pCD1, pzb_rho_fix );
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* next 2 lines for debug only */
    // INCHI✔️❌:     /* num_g++;  */
    // INCHI✔️❌:     /* WriteGraph( G, n_tg, num_g, "V:\\IChI_v10\\Gordon-Graphs\\hard\\k06g08v312-alt.dre", "a+" ); */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* memory allocation */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (0 > SetBitCreate( pCG ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (pzb_rho_fix && pzb_rho_fix->nLenCTAtOnly != pCD->nLenCTAtOnly)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* consistency check */
    // INCHI✔️❌:         return -2;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     ok = 1;
    // INCHI✔️❌:
    // INCHI✔️❌:     ok &= UnorderedPartitionCreate( &theta, n_tg );
    // INCHI✔️❌:     ok &= UnorderedPartitionCreate( &theta_from_gamma, n_tg );
    // INCHI✔️❌:
    // INCHI✔️❌:     /* djb-rwth: fixing dereferencing NULL pointer and buffer overflows */
    // INCHI✔️❌:     W = (Cell*) inchi_calloc( n_tg, sizeof( W[0] ) );
    // INCHI✔️❌:     v = (Node*) inchi_calloc( n_tg, sizeof( v[0] ) );
    // INCHI✔️❌:     e = (S_CHAR*) inchi_calloc( n_tg, sizeof( e[0] ) );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!W || !v || !e)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ok &= 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: /*
    // INCHI✔️❌:     ok &= ( NULL != ( W = (Cell*) inchi_calloc( n_tg, sizeof( W[0] ) ) ) );
    // INCHI✔️❌:     ok &= ( NULL != ( v = (Node*) inchi_calloc( n_tg, sizeof( v[0] ) ) ) );
    // INCHI✔️❌:     ok &= ( NULL != ( e = (S_CHAR*) inchi_calloc( n_tg, sizeof( e[0] ) ) ) );
    // INCHI✔️❌: */
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     ok &= (NULL != (v   = (Node*)inchi_calloc( n_tg, sizeof(W[0]))));
    // INCHI✔️❌:     ok &= (NULL != (e   = (S_CHAR*)inchi_calloc( n_tg, sizeof(W[0]))));
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*    ok &= (NULL != (qzb = (S_CHAR*)inchi_calloc( n_tg, sizeof(W[0])))); */
    // INCHI✔️❌:     ok &= CTableCreate( &Lambda, n, pCD );
    // INCHI✔️❌:     ok &= CTableCreate( &zf_zeta, n, pCD );
    // INCHI✔️❌:     ok &= ( ( pzb_rho = (ConTable *) inchi_calloc( 1, sizeof( *pzb_rho ) ) ) &&
    // INCHI✔️❌:             CTableCreate( pzb_rho, n, pCD ) );
    // INCHI✔️❌:
    // INCHI✔️❌:     ok &= NodeSetCreate( pCG, &Omega, n_tg, L_curr_max_set_size );
    // INCHI✔️❌:     ok &= NodeSetCreate( pCG, &Phi, n_tg, L_curr_max_set_size );
    // INCHI✔️❌:     ok &= NodeSetCreate( pCG, &cur_nodes, n_tg, 1 );
    // INCHI✔️❌:
    // INCHI✔️❌:     ok &= PartitionCreate( &zeta, n_tg );
    // INCHI✔️❌:     ok &= PartitionCreate( &rho, n_tg );
    // INCHI✔️❌:     ok &= TranspositionCreate( &gamma, n_tg );
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌: /*L1:*/
    // INCHI✔️❌:     k = 1;
    // INCHI✔️❌:     size = 1.0;
    // INCHI✔️❌:     hz_rho = index = l = 0; /* djb-rwth: removing redundant code */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!ok)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function; /* initialization failed */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     UnorderedPartitionMakeDiscrete( &theta, n_tg );
    // INCHI✔️❌:     t_Lemma = 2;
    // INCHI✔️❌:
    // INCHI✔️❌:     pCC->lNumBreakTies = 0;
    // INCHI✔️❌:     pCC->lNumDecreasedCT = 0;
    // INCHI✔️❌:     pCC->lNumRejectedCT = 0;
    // INCHI✔️❌:     pCC->lNumEqualCT = 1;
    // INCHI✔️❌:     pCC->lNumTotCT = 0;
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:
    // INCHI✔️❌:     hzb_rho_fix = 1;
    // INCHI✔️❌:
    // INCHI✔️❌:     memset( kLeast_rho, 0, sizeof( kLeast_rho ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( kLeast_rho_fix, 0, sizeof( kLeast_rho_fix ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (PartitionIsDiscrete( &pi[k - 1], n_tg ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* added the following 3 lines to the original to create Ct */
    // INCHI✔️❌:         PartitionCopy( &rho, &pi[k - 1], n_tg );
    // INCHI✔️❌:         CtPartFill( G, pCD, &pi[k - 1], pzb_rho, 1, n, n_tg, n_max );
    // INCHI✔️❌:         CtPartINCHI_CANON_INFINITY( pzb_rho, qzb, 2 );
    // INCHI✔️❌:         pCC->lNumTotCT++;
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         /* goto L18; */
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!dig && PartitionSatisfiesLemma_2_25( &pi[0], n ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         t_Lemma = 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     PartitionGetFirstCell( &pi[k-1], &W[k-1], k, n );
    // INCHI✔️❌:     v[k-1] = CellGetMinNode( &pi[k-1], &W[k-1], 0, pCD1 );
    // INCHI✔️❌:     CtPartClear( &Lambda, 1 );
    // INCHI✔️❌:     e[k-1] = 0;
    // INCHI✔️❌:     */
    // INCHI✔️❌:     CtPartClear( &Lambda, 1 );
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌: /* L2: reach the first leaf and save it in zeta and rho */
    // INCHI✔️❌:     while (k)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* the two next lines intentionally switched */
    // INCHI✔️❌:         /* Create equitable partition in pi[k]  */
    // INCHI✔️❌:         PartitionGetFirstCell( &pi[k - 1], W, k, n );
    // INCHI✔️❌:         v[k - 1] = CellGetMinNode( &pi[k - 1], &W[k - 1], 0, pCD1 );
    // INCHI✔️❌:         e[k - 1] = 0;
    // INCHI✔️❌:         if (dig || !PartitionSatisfiesLemma_2_25( &pi[k - 1], n ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             t_Lemma = k + 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* e[k-1] = 0; */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             Node vv = v[k - 1];
    // INCHI✔️❌:             if (0 > ( ret = PartitionColorVertex( pCG, G, &pi[k - 1], vv /*v[k-1]*/, n, n_tg, n_max, bDigraph, 0 ) ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto exit_error;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         pCC->lNumBreakTies++;
    // INCHI✔️❌:         k++;
    // INCHI✔️❌:
    // INCHI✔️❌:         CtPartFill( G, pCD, &pi[k - 1], &Lambda, k - 1, n, n_tg, n_max);
    // INCHI✔️❌:
    // INCHI✔️❌:         /* return -1; *//* debug only */
    // INCHI✔️❌:         /* if(h_zeta==0)goto L5; L5: */
    // INCHI✔️❌:         /* the first terminal node has not been reached yet */
    // INCHI✔️❌:         /* search for the predefined numbering */
    // INCHI✔️❌:         if (pzb_rho_fix && QZFIX_OK( qzb_rho_fix ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             qzb_rho_fix = CtPartCompare( &Lambda, pzb_rho_fix, qzb, kLeast_rho_fix, k - 1, 1, bSplitTautCompare );
    // INCHI✔️❌:             if (QZFIX_OK( qzb_rho_fix ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 hzb_rho_fix = k;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (lab && QZFIX_OK( qzb_rho_fix ))  /* DCh */
    // INCHI✔️❌:             CtPartCopy( pzb_rho, &Lambda, k - 1 );
    // INCHI✔️❌:         CtPartCopy( &zf_zeta, &Lambda, k - 1 );
    // INCHI✔️❌:         /*goto L4; L4:*/
    // INCHI✔️❌:         if (PartitionIsDiscrete( &pi[k - 1], n ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             break;  /* goto L7; */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* goto L2; */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     pCC->lNumTotCT++;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* L7; L7: */
    // INCHI✔️❌:     /* if ( h_zeta == 0 ) goto L18; L18:*/
    // INCHI✔️❌:     h_zeta = t_eq_zeta = hz_zeta = k;
    // INCHI✔️❌:     CtPartINCHI_CANON_INFINITY( &zf_zeta, NULL, k );
    // INCHI✔️❌:     /******************** <<<===== B **************************/
    // INCHI✔️❌:     PartitionCopy( &zeta, &pi[k - 1], n_tg );
    // INCHI✔️❌:     if (lab)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (pzb_rho_fix)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (0 == qzb_rho_fix)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 qzb_rho_fix = CtFullCompare( &Lambda, pzb_rho_fix, 1, bSplitTautCompare );
    // INCHI✔️❌:                 if (qzb_rho_fix > 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     hzb_rho_fix = 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (hzb_rho_fix > 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 PartitionCopy( &rho, &pi[hzb_rho_fix - 1], n_tg );
    // INCHI✔️❌:                 /*CtPartINCHI_CANON_INFINITY( pzb_rho, qzb, k );*/
    // INCHI✔️❌:             }
    // INCHI✔️❌:             hz_rho = h_rho = hzb_rho_fix;
    // INCHI✔️❌:             bRhoIsDiscrete = ( hzb_rho_fix == k );
    // INCHI✔️❌:             if (bRhoIsDiscrete)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 CtPartINCHI_CANON_INFINITY( pzb_rho, qzb, k );
    // INCHI✔️❌:                 pzb_rho_fix_reached = !qzb_rho_fix;
    // INCHI✔️❌:                 CtCompareLayersGetFirstDiff( kLeast_rho_fix, nOneAdditionalLayer,
    // INCHI✔️❌:                                  &L_rho_fix_prev, &I_rho_fix_prev, &k_rho_fix_prev );
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI❌❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI❌❌:             else
    // INCHI❌❌:             {
    // INCHI❌❌:                 int stop = 1;
    // INCHI❌❌:             }
    // INCHI❌❌: #endif
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             PartitionCopy( &rho, &pi[k - 1], n_tg );
    // INCHI✔️❌:             hz_rho = h_rho = k;
    // INCHI✔️❌:             CtPartINCHI_CANON_INFINITY( pzb_rho, qzb, k );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         qzb_rho = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     r = k;
    // INCHI✔️❌:     v[k - 1] = INCHI_CANON_INFINITY;     /* DCh */
    // INCHI✔️❌:     CellMakeEmpty( W, k ); /* DCh */
    // INCHI✔️❌:     k--;
    // INCHI✔️❌:     goto L13;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: L2:
    // INCHI✔️❌:     /* the two next lines intentionally switched */
    // INCHI✔️❌:     /* Create equitable partition in pi[k]  */
    // INCHI✔️❌:     if (0 > ( ret = PartitionColorVertex( pCG, G, &pi[k - 1], v[k - 1],
    // INCHI✔️❌:         n, n_tg, n_max, bDigraph, 0 ) ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_error;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     pCC->lNumBreakTies++;
    // INCHI✔️❌:     k++;
    // INCHI✔️❌:     CtPartFill( G, pCD, &pi[k - 1], &Lambda, k - 1, n, n_tg, n_max );
    // INCHI✔️❌:     e[k - 1] = 0;         /* moved  */
    // INCHI✔️❌:     v[k - 1] = INCHI_CANON_INFINITY;  /* added by DCh. */
    // INCHI✔️❌:     CellMakeEmpty( W, k ); /* DCh */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (hz_zeta == k - 1 &&
    // INCHI✔️❌:          0 == CtPartCompare( &Lambda, &zf_zeta, NULL, NULL, k - 1, 0, bSplitTautCompare ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         hz_zeta = k; /* max{k|Lambda(G,pi,nu(k))==Lambda(G,pi,zeta) }  */
    // INCHI✔️❌:     } /* added */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* -- old code ---
    // INCHI✔️❌:     if ( pzb_rho_fix && QZFIX_OK(qzb_rho_fix) ) {
    // INCHI✔️❌:         qzb_rho_fix = CtPartCompare( &Lambda, pzb_rho_fix, qzb, kLeast_rho_fix, k-1, 1, bSplitTautCompare );
    // INCHI✔️❌:         if ( QZFIX_OK(qzb_rho_fix) ) {
    // INCHI✔️❌:             hzb_rho_fix = k;
    // INCHI✔️❌:         } else {
    // INCHI✔️❌:             pCC->lNumRejectedCT ++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* --- new code ---*/
    // INCHI✔️❌:     if (pzb_rho_fix && !qzb_rho_fix)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         qzb_rho_fix = CtPartCompare( &Lambda, pzb_rho_fix, qzb, kLeast_rho_fix, k - 1, 1, bSplitTautCompare );
    // INCHI✔️❌:         if (!qzb_rho_fix && bRhoIsDiscrete)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             qzb_rho_fix = CtPartCompareLayers( kLeast_rho_fix, L_rho_fix_prev, nOneAdditionalLayer );
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( FIX_ChCh_CONSTIT_CANON_BUG == 1 )
    // INCHI✔️❌:             if (qzb_rho_fix)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 int L_rho_fix_diff = abs( qzb_rho_fix ) - 1;
    // INCHI✔️❌:                 if (L_rho_fix_diff < L_rho_fix_prev ||
    // INCHI✔️❌:                      (L_rho_fix_diff == (L_rho_fix_prev && kLeast_rho_fix[L_rho_fix_diff].i < I_rho_fix_prev))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     qzb_rho_fix = L_rho_fix_diff + 1; /* positive difference will be rejected */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI❌❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI❌❌:             if (qzb_rho_fix)
    // INCHI❌❌:             {
    // INCHI❌❌:                 int stop = 1; /* debug only */
    // INCHI❌❌:             }
    // INCHI❌❌: #endif
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (!QZFIX_OK( qzb_rho_fix ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pCC->lNumRejectedCT++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (pzb_rho_fix && QZFIX_OK( qzb_rho_fix ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         hzb_rho_fix = k;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* if (!lab) goto L3; */
    // INCHI✔️❌:     if (lab && QZFIX_OK( qzb_rho_fix ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* once the difference has been found it is meaningful as long as k increments */
    // INCHI✔️❌:         /* cur_qzb2 = CtPartCompare( &Lambda, pzb_rho, qzb, k-1 ); */ /* rho compare */
    // INCHI✔️❌:         if (hz_rho == k - 1 && !qzb_rho && bRhoIsDiscrete)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int qzb_rho_temp = 0; /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:             qzb_rho = CtPartCompare( &Lambda, pzb_rho, qzb, kLeast_rho, k - 1, 0, bSplitTautCompare );
    // INCHI✔️❌:             /* old code */
    // INCHI✔️❌:             if (!qzb_rho && pzb_rho_fix_reached &&
    // INCHI✔️❌:                   nOneAdditionalLayer && 0 > kLeast_rho[nOneAdditionalLayer].k)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 qzb_rho_temp = -( nOneAdditionalLayer + 1 ); /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:                 /* qzb_rho = -(nOneAdditionalLayer+1); *//* early rejection */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* new code */
    // INCHI✔️❌:             if (!qzb_rho && bRhoIsDiscrete)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 qzb_rho = CtPartCompareLayers( kLeast_rho, L_rho_fix_prev, 0 );
    // INCHI✔️❌: #if ( FIX_ChCh_CONSTIT_CANON_BUG == 1 )
    // INCHI✔️❌:                 if (qzb_rho)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     int L_rho_diff = abs( qzb_rho ) - 1;
    // INCHI✔️❌:                     if (L_rho_diff < L_rho_fix_prev ||
    // INCHI✔️❌:                          (L_rho_diff == (L_rho_fix_prev && kLeast_rho[L_rho_diff].i < I_rho_fix_prev))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         qzb_rho = -( L_rho_diff + 1 ); /* negative difference will be rejected */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* compare old results to new */
    // INCHI❌❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI❌❌:             if (qzb_rho_temp && qzb_rho_temp != qzb_rho)
    // INCHI❌❌:             {
    // INCHI❌❌:                 int stop = 1; /* <BRKPT> */
    // INCHI❌❌:             }
    // INCHI❌❌: #endif
    // INCHI✔️❌:             if (!qzb_rho)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 hz_rho = k;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (qzb_rho < 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     pCC->lNumRejectedCT++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (qzb_rho > 0 || (!qzb_rho && !bRhoIsDiscrete)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* found better rho */
    // INCHI✔️❌:             if (!nNumLayers)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 CtPartCopy( pzb_rho, &Lambda, k - 1 );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: /*L3:*/
    // INCHI✔️❌:     /*if ( hz_rho == k || (lab && qzb_rho >= 0 ) )*/
    // INCHI✔️❌:     /*if ( hz_zeta == k || hz_rho == k || (lab && qzb_rho >= 0 ) ) goto L4; else goto L6;*/
    // INCHI✔️❌:     if (hz_zeta == k || hz_rho == k ||
    // INCHI✔️❌:         ( lab && qzb_rho >= 0 && QZFIX_OK( qzb_rho_fix ) ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*L4: check for possible isomorphism or found a better rho */
    // INCHI✔️❌:         if (PartitionIsDiscrete( &pi[k - 1], n ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pCC->lNumTotCT++;
    // INCHI✔️❌:             goto L7;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         PartitionGetFirstCell( &pi[k - 1], W, k, n );
    // INCHI✔️❌:         v[k - 1] = CellGetMinNode( &pi[k - 1], &W[k - 1], 0, pCD1 );
    // INCHI✔️❌:         if (!dig && PartitionSatisfiesLemma_2_25( &pi[k - 1], n ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ; /* found additional isomprphism */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             t_Lemma = k + 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         e[k - 1] = 0;  /* created new cell W[k-1] */
    // INCHI✔️❌:         goto L2;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: L6:
    // INCHI✔️❌:     /* a better rho or no good node was found at this level; return to smaller k */
    // INCHI✔️❌:     k2 = k;
    // INCHI✔️❌:     k = inchi_min( t_Lemma - 1, inchi_max( t_eq_zeta - 1, hz_rho ) );
    // INCHI✔️❌:     if (k2 == t_Lemma)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto L13;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* store isomorphism found from Lemma 2.25. should be dig=0 !!! */
    // INCHI✔️❌:     if (dig)
    // INCHI✔️❌:     {
    // INCHI❌❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI❌❌:         int stop = 1;
    // INCHI❌❌: #endif
    // INCHI✔️❌:         goto L13;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     l = inchi_min( l + 1, L_curr_max_set_size );
    // INCHI✔️❌:     PartitionGetMcrAndFixSet( pCG, &pi[t_Lemma - 1], &Omega, &Phi, n_tg, l );
    // INCHI✔️❌:     goto L12;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: L7:
    // INCHI✔️❌:     /* from L4: pi[k-1] is discrete */
    // INCHI✔️❌:     if (h_zeta == 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*goto L18;*/  /* error. the first T(nu) leaf was found */
    // INCHI✔️❌:         ret = CT_CANON_ERR;
    // INCHI✔️❌:         goto exit_error;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (k != hz_zeta)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto L8;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*  here zeta^gamma == nu */
    // INCHI✔️❌:     /*  if ( G^gamma == G ) goto L10; */
    // INCHI✔️❌:     if (0 == ( res = CtFullCompare( &Lambda, &zf_zeta, 0, bSplitTautCompare ) )) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         PartitionGetTransposition( &zeta, &pi[k - 1], n_tg, &gamma );
    // INCHI✔️❌:         bZetaIsomorph = 1; /* for testing only */
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         goto L10;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI❌❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI❌❌:     {
    // INCHI❌❌:         int stop = 1;
    // INCHI❌❌:     }
    // INCHI❌❌: #endif
    // INCHI✔️❌:     /* !!! we should never come here !!! */
    // INCHI✔️❌:     if (!nNumLayers)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = -2;
    // INCHI✔️❌:         goto exit_error;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: L8: /* here nu is discrete: check rho for being a bettere leaf or isomorphism */
    // INCHI✔️❌:     /*if ( !lab || qzb_rho < 0 || !QZFIX_OK(qzb_rho_fix) )*/
    // INCHI✔️❌:     if (!lab || (qzb_rho < 0 && ( !pzb_rho_fix || qzb_rho_fix > 0 ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto L6;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (pzb_rho_fix && kLeast_rho_fix && 0 == qzb_rho_fix) /* djb-rwth: addressing of array kLeast_rho_fix will always evaluate to true? */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* check for the rejection condition: Lambda > zb_rho_fix */
    // INCHI✔️❌:         if (kLeast_rho_fix) /* djb-rwth: addressing of array kLeast_rho_fix will always evaluate to true? */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int qzb_rho_fix_alt;
    // INCHI✔️❌:             qzb_rho_fix = CtFullCompareLayers( kLeast_rho_fix );
    // INCHI✔️❌:             /* for debug only */
    // INCHI✔️❌:             qzb_rho_fix_alt = CtFullCompare( &Lambda, pzb_rho_fix, 1, bSplitTautCompare );
    // INCHI✔️❌:             if (qzb_rho_fix != qzb_rho_fix_alt)
    // INCHI✔️❌:             {
    // INCHI❌❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI❌❌:                 int stop = 1;
    // INCHI❌❌: #endif
    // INCHI✔️❌:                 qzb_rho_fix = qzb_rho_fix_alt;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* end debug */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             qzb_rho_fix = CtFullCompare( &Lambda, pzb_rho_fix, 1, bSplitTautCompare );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (!pzb_rho_fix_reached)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pzb_rho_fix_reached = !qzb_rho_fix;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (0 < qzb_rho_fix)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* Lambda > pzb_rho_fix, ignore this node */
    // INCHI✔️❌:             /* hzb_rho_fix = min( hzb_rho_fix, hz_rho ); */ /* ??? */
    // INCHI✔️❌:             qzb_rho_fix = 0;
    // INCHI✔️❌:             goto L6;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         qzb_rho_fix = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (qzb_rho < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto L6;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (qzb_rho > 0 || !bRhoIsDiscrete)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto L9; /* note: p67 says k > PartitionSize( &rho, n ) */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (k < r)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto L9; /* cannot understand it... */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* !!! we should never come here if G(nu) != G(rho): CtPartCompare must be enough !!! */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* if ( G(nu) > G(rho) ) goto L9; */
    // INCHI✔️❌:     if (kLeast_rho) /* djb-rwth: addressing of array kLeast_rho will always evaluate to true? */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int cur_qzb_alt;
    // INCHI✔️❌:         qzb_rho = CtFullCompareLayers( kLeast_rho );
    // INCHI✔️❌:         /* for debug only */
    // INCHI✔️❌:         cur_qzb_alt = CtFullCompare( &Lambda, pzb_rho, 0, bSplitTautCompare );
    // INCHI✔️❌:         if (qzb_rho != cur_qzb_alt)
    // INCHI✔️❌:         {
    // INCHI❌❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI❌❌:             int stop = 1;
    // INCHI❌❌: #endif
    // INCHI✔️❌:             qzb_rho = cur_qzb_alt;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* end debug */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         qzb_rho = CtFullCompare( &Lambda, pzb_rho, 0, bSplitTautCompare );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* qzb_rho difference can be due to layers 1..MAX_LAYERS-1 only */
    // INCHI✔️❌:     if (0 < qzb_rho)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* CtFullCompare( &Lambda, pzb_rho, 0, bSplitTautCompare ); */
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         goto L9;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* if ( G(nu) < G(rho) ) goto L6; */
    // INCHI✔️❌:     if (0 > qzb_rho)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         qzb_rho = 0;
    // INCHI✔️❌:         goto L6;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* nu^gamma == rho */
    // INCHI✔️❌:     if (r != k)
    // INCHI✔️❌:     {  /* if() is for debug only */
    // INCHI✔️❌:         r = k;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     PartitionGetTransposition( &pi[k - 1], &rho, n_tg, &gamma );
    // INCHI✔️❌:     bZetaIsomorph = 0; /* DCh */
    // INCHI✔️❌:     pCC->lNumEqualCT++;
    // INCHI✔️❌:     goto L10;
    // INCHI✔️❌:
    // INCHI✔️❌: L9:
    // INCHI✔️❌:     /* rho := nu; */
    // INCHI✔️❌:     PartitionCopy( &rho, &pi[k - 1], n_tg );
    // INCHI✔️❌:     if (nNumLayers)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         CtFullCopy( pzb_rho, &Lambda );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     bZetaEqRho = 0;
    // INCHI✔️❌:     qzb_rho = 0;
    // INCHI✔️❌:     CtCompareLayersGetFirstDiff( kLeast_rho_fix, nOneAdditionalLayer,
    // INCHI✔️❌:                                  &L_rho_fix_prev, &I_rho_fix_prev,
    // INCHI✔️❌:                                  &k_rho_fix_prev );
    // INCHI✔️❌:     memset( kLeast_rho, 0, sizeof( kLeast_rho ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     h_rho = hz_rho = k;
    // INCHI✔️❌:     CtPartINCHI_CANON_INFINITY( pzb_rho, qzb, k );
    // INCHI✔️❌:     pCC->lNumDecreasedCT++;
    // INCHI✔️❌:     pCC->lNumEqualCT = 1;
    // INCHI✔️❌:     bRhoIsDiscrete = 1;
    // INCHI✔️❌:     goto L6;
    // INCHI✔️❌:
    // INCHI✔️❌: L10: /* discrete pi[k-1] && G^gamma == G */
    // INCHI✔️❌:     pCC->lNumEqualCT += bZetaEqRho || !( bZetaIsomorph || qzb_rho );
    // INCHI✔️❌:     l = inchi_min( l + 1, L_curr_max_set_size );
    // INCHI✔️❌:     /* Omega[l] := mcr(gamma);
    // INCHI✔️❌:        Phi[l]   := fix(gamma);
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     TranspositionGetMcrAndFixSetAndUnorderedPartition( pCG, &gamma, &Omega, &Phi,
    // INCHI✔️❌:                                                        n_tg, l, &theta_from_gamma );
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( theta(gamma) <= theta ) goto L11;
    // INCHI✔️❌:     theta := theta v theta(gamma);
    // INCHI✔️❌:     UnorderedPartitionJoin() returns 0 if theta_from_gamma is finer than theta,
    // INCHI✔️❌:     which means no changes in theta: theta_from_gamma ^ theta == theta.
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if (!UnorderedPartitionJoin( &theta_from_gamma, &theta, n_tg ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto L11; /* no new isomorphism found */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  Output gamma (it is the Aut(G) generator) -- omitted -- */
    // INCHI✔️❌:     nNumFoundGenerators++;
    // INCHI✔️❌:     /* if ( tvc in mcr(theta) ) goto L11; */
    // INCHI✔️❌:     if (tvc == GetUnorderedPartitionMcrNode( &theta, tvc ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto L11;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     k = h_zeta;
    // INCHI✔️❌:     goto L13;
    // INCHI✔️❌:
    // INCHI✔️❌: L11:
    // INCHI✔️❌:     k = lab ? h_rho : h_zeta; /***Changed*** originally was k = h_rho; */
    // INCHI✔️❌:
    // INCHI✔️❌: L12:
    // INCHI✔️❌:     /* if ( e[k-1] == 1 ) */
    // INCHI❌❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI❌❌:     if (e[k - 1] == 1 && v[k - 1] == INCHI_CANON_INFINITY)
    // INCHI❌❌:     {
    // INCHI❌❌:         int stop = 1;          /* <BRKPT> testing only */
    // INCHI❌❌:     }
    // INCHI❌❌: #endif
    // INCHI✔️❌:     if (e[k - 1] == 1 && v[k - 1] != INCHI_CANON_INFINITY)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* INCHI_CANON_INFINITY for testing only */
    // INCHI✔️❌:         CellIntersectWithSet( pCG, &pi[k - 1], &W[k - 1], &Omega, l );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: L13:
    // INCHI✔️❌:
    // INCHI✔️❌:     if ((UserAction && USER_ACTION_QUIT == ( *UserAction )( )) ||
    // INCHI✔️❌:          (ConsoleQuit && ( *ConsoleQuit )( ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = CT_USER_QUIT_ERR;
    // INCHI✔️❌:         goto exit_error;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (bInchiTimeIsOver( ic, pCD->ulTimeOutTime ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = CT_TIMEOUT_ERR;
    // INCHI✔️❌:         goto exit_error;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (k == 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function; /* stop */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (lab && k < h_rho)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /***Added***/
    // INCHI✔️❌:         h_rho = k;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (k > h_zeta)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (v[k - 1] == INCHI_CANON_INFINITY)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*** Added by DCh for testing only ****/
    // INCHI✔️❌:             k--;
    // INCHI✔️❌:             goto L13;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         goto L17;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (k == h_zeta)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto L14;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     h_zeta = k;
    // INCHI✔️❌:     tvc = tvh = CellGetMinNode( &pi[k - 1], &W[k - 1], 0, pCD1 );
    // INCHI✔️❌:
    // INCHI✔️❌: L14:
    // INCHI✔️❌:     /* if v[k] and tvh are in the same cell of theta then index ++ */
    // INCHI✔️❌:     if (GetUnorderedPartitionMcrNode( &theta, v[k - 1] ) ==
    // INCHI✔️❌:          GetUnorderedPartitionMcrNode( &theta, tvh ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         index++;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     v[k - 1] = CellGetMinNode( &pi[k - 1], &W[k - 1], v[k - 1], pCD1 );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (v[k - 1] == INCHI_CANON_INFINITY)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto L16;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (v[k - 1] != GetUnorderedPartitionMcrNode( &theta, v[k - 1] ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto L14;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: L15:
    // INCHI✔️❌:     t_Lemma = inchi_min( t_Lemma, k + 1 );
    // INCHI✔️❌:     hz_zeta = inchi_min( hz_zeta, k );
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( lab && hz_rho >= k ) {
    // INCHI✔️❌:         hz_rho = k;
    // INCHI✔️❌:         qzb_rho = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if (lab)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (hz_rho >= k /*-1*/)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             qzb_rho = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (hz_rho > k)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             hz_rho = k;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         UpdateCompareLayers( kLeast_rho, hz_rho );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (pzb_rho_fix)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (hzb_rho_fix >= k /*-1*/)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             qzb_rho_fix = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (hzb_rho_fix > k)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             hzb_rho_fix = k;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         UpdateCompareLayers( kLeast_rho_fix, hzb_rho_fix );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     goto L2;
    // INCHI✔️❌:
    // INCHI✔️❌: L16:
    // INCHI✔️❌:     if (t_eq_zeta == k + 1 && index == CellGetNumberOfNodes( &pi[k - 1], &W[k - 1] ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         t_eq_zeta = k;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     size *= (double) index;
    // INCHI✔️❌:     /******************** <<<===== A **************************/
    // INCHI✔️❌:     /* passed K times after passing point A. At these passes
    // INCHI✔️❌:        k = K, K-1, ..., 1 in this order
    // INCHI✔️❌:     */
    // INCHI✔️❌:     index = 0;
    // INCHI✔️❌:     k--;
    // INCHI✔️❌:     goto L13;
    // INCHI✔️❌:
    // INCHI✔️❌: L17:
    // INCHI✔️❌:     /* if ( e[k-1] == 0 ) */
    // INCHI❌❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI❌❌:     if (e[k - 1] == 0 && v[k - 1] == INCHI_CANON_INFINITY)
    // INCHI❌❌:     {
    // INCHI❌❌:         /* testing only */
    // INCHI❌❌:         int stop = 1;  /* <BRKPT> */
    // INCHI❌❌:     }
    // INCHI❌❌: #endif
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( e[k] == 0 set W[k] = Intersection(W[k], Omega[i]) for each i = 1..l,
    // INCHI✔️❌:          such that {v[1]..v[k-1]} in Phi[i]
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if (e[k - 1] == 0 && v[k - 1] != INCHI_CANON_INFINITY) /* Added v[k-1]!=... DCh */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         NodeSetFromVertices( pCG, &cur_nodes, 1, v, k - 1 );
    // INCHI✔️❌:         for (i = 1; i <= l; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (AllNodesAreInSet( &cur_nodes, 1, &Phi, i ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 CellIntersectWithSet( pCG, &pi[k - 1], &W[k - 1], &Omega, i );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     e[k - 1] = 1;
    // INCHI✔️❌:     v[k - 1] = CellGetMinNode( &pi[k - 1], &W[k - 1], v[k - 1], pCD1 );
    // INCHI✔️❌:     if (v[k - 1] != INCHI_CANON_INFINITY)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto L15;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     k--;
    // INCHI✔️❌:     goto L13;
    // INCHI✔️❌: /* L18: see above */
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     /* CtPartFill( G, pCD, &rho, pzb_rho, 1, n, n_tg ); */
    // INCHI✔️❌:     if (!bRhoIsDiscrete)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = CT_CANON_ERR;
    // INCHI✔️❌:         goto exit_error;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (pzb_rho_fix)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         qzb_rho_fix = CtFullCompare( pzb_rho_fix, pzb_rho, 1, bSplitTautCompare );
    // INCHI✔️❌:         if (qzb_rho_fix)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = CT_CANON_ERR;
    // INCHI✔️❌:             goto exit_error;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* SymmRank */
    // INCHI✔️❌:     memset( nSymmRank, 0, n_tg * sizeof( nSymmRank[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     for (i = 0; i < n_tg; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         k = rho.AtNumber[i];
    // INCHI✔️❌:         k2 = (int) GetUnorderedPartitionMcrNode( &theta, (AT_NUMB) ( k + 1 ) ) - 1;
    // INCHI✔️❌:         if (!nSymmRank[k2] || nSymmRank[k2] > rho.Rank[k])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nSymmRank[k2] = rho.Rank[k];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (i = 0; i < n_tg; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         k = rho.AtNumber[i];
    // INCHI✔️❌:         k2 = (int) GetUnorderedPartitionMcrNode( &theta, (AT_NUMB) ( k + 1 ) ) - 1;
    // INCHI✔️❌:         nSymmRank[k] = nSymmRank[k2];
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* CanonRank, nAtomNumberCanon */
    // INCHI✔️❌:     memcpy(nCanonRank, rho.Rank, n_tg * sizeof(nCanonRank[0]));
    // INCHI✔️❌:     memcpy(nAtomNumberCanon, rho.AtNumber, n_tg * sizeof(nAtomNumberCanon[0]));
    // INCHI✔️❌:     /* LinearCT */
    // INCHI✔️❌:     if (pzb_rho)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *nLenCt = pzb_rho->lenCt - 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (pCt && pzb_rho && (*nLenCt > 0)) /* djb-rwth: GHI #164 fix */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         memcpy(pCt, pzb_rho->Ctbl, *nLenCt * sizeof(pCt[0]));
    // INCHI✔️❌:     }
    // INCHI✔️❌:     pCC->lNumTotCT = pCC->lNumDecreasedCT + pCC->lNumRejectedCT + pCC->lNumEqualCT;
    // INCHI✔️❌:     pCC->dGroupSize = size;
    // INCHI✔️❌:     pCC->lNumGenerators = nNumFoundGenerators;
    // INCHI✔️❌:     pCC->lNumStoredIsomorphisms = l;
    // INCHI✔️❌:     /* Note: check nNumFoundGenerators */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (pp_zb_rho_out && !*pp_zb_rho_out)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *pp_zb_rho_out = pzb_rho;
    // INCHI✔️❌:         pzb_rho = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: exit_error:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:     UnorderedPartitionFree( &theta );
    // INCHI✔️❌:     UnorderedPartitionFree( &theta_from_gamma );
    // INCHI✔️❌:     if (W)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( W );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (v)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( v );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (e)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( e );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (qzb)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( qzb );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     CTableFree( &Lambda );
    // INCHI✔️❌:     CTableFree( &zf_zeta );
    // INCHI✔️❌:     if (pzb_rho)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         CTableFree( pzb_rho );
    // INCHI✔️❌:         inchi_free( pzb_rho );
    // INCHI✔️❌:         pzb_rho = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* CTableFree( &zf_zeta2 ); */
    // INCHI✔️❌:     NodeSetFree( pCG, &Omega );
    // INCHI✔️❌:     NodeSetFree( pCG, &Phi );
    // INCHI✔️❌:     /* NodeSetFree( &mcr_theta, n, 1 ); */
    // INCHI✔️❌:     NodeSetFree( pCG, &cur_nodes );
    // INCHI✔️❌:     PartitionFree( &zeta );
    // INCHI✔️❌:     /* PartitionFree( &zeta2 ); */
    // INCHI✔️❌:     PartitionFree( &rho );
    // INCHI✔️❌:     TranspositionFree( &gamma );
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CanonGraph
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CanonGraph
    // INCHI✔️❌: #define COMPILE_ANSI_ONLY
    // INCHI✔️❌: #define TARGET_API_LIB
    // INCHI✔️❌: GCC/Linux branch
    // INCHI✔️❌: #define bRELEASE_VERSION 1
    // INCHI✔️❌: #define USE_ISO_SORT_KEY_HFIXED 0
    // INCHI✔️❌: #define FIX_ChCh_CONSTIT_CANON_BUG 1
    // INCHI❌❌: INCHI_CANON_USE_HASH is undefined
    // END INCHI ACTIVE MACRO CONFIGURATION: CanonGraph
    let pzb_rho_fix_storage = if pzb_rho_inp.is_null() {
        None
    } else {
        // SAFETY: CanonGraph only reads the caller-owned input ConTable, and
        // the allocation remains live for the complete source call.
        Some(unsafe { heap.stable_slice(pzb_rho_inp.as_const())? })
    };
    let pzb_rho_fix = if let Some(storage) = pzb_rho_fix_storage.as_ref() {
        Some(storage.get(0)?)
    } else {
        None
    };
    if SetBitCreate(heap, pCG)? < 0 {
        return Ok(-1);
    }
    if pzb_rho_fix.is_some_and(|fixed| fixed.nLenCTAtOnly != pCD.nLenCTAtOnly) {
        return Ok(-2);
    }
    let (mut theta, mut theta_from_gamma) =
        (UnorderedPartition::default(), UnorderedPartition::default());
    let (mut Lambda, mut zf_zeta) = (ConTable::default(), ConTable::default());
    let (mut Omega, mut Phi, mut cur_nodes) =
        (NodeSet::default(), NodeSet::default(), NodeSet::default());
    let (mut zeta, mut rho) = (Partition::default(), Partition::default());
    let mut gamma = Transposition::default();
    // INCHI✔️✔️:     kLeast kLeast_rho[MAX_LAYERS];
    let mut kLeast_rho_storage = MaybeUninit::<KLeastLayers>::uninit();
    // INCHI✔️✔️:     kLeast kLeast_rho_fix[MAX_LAYERS];
    let mut kLeast_rho_fix_storage = MaybeUninit::<KLeastLayers>::uninit();
    let mut ok = 1_i32;
    ok &= UnorderedPartitionCreate(heap, &mut theta, n_tg)?;
    ok &= UnorderedPartitionCreate(heap, &mut theta_from_gamma, n_tg)?;
    let W_pointer = crate::source::base::util::inchi_calloc::<Cell>(
        heap,
        n_tg as u64,
        size_of::<Cell>() as u64,
    )
    .unwrap_or_else(|_| SourceMutPointer::null());
    let v = crate::source::base::util::inchi_calloc::<Node>(
        heap,
        n_tg as u64,
        size_of::<Node>() as u64,
    )
    .unwrap_or_else(|_| SourceMutPointer::null());
    let e = crate::source::base::util::inchi_calloc::<crate::source_types::S_CHAR>(
        heap,
        n_tg as u64,
        size_of::<crate::source_types::S_CHAR>() as u64,
    )
    .unwrap_or_else(|_| SourceMutPointer::null());
    if W_pointer.is_null() || v.is_null() || e.is_null() {
        ok = 0;
    }
    let mut W = if W_pointer.is_null() {
        None
    } else {
        // SAFETY: this is the unique writable view of the source `Cell *W`
        // allocation. CanonGraph keeps it live until the source cleanup block.
        Some(unsafe { heap.stable_slice_mut(W_pointer)? })
    };
    ok &= CTableCreate(heap, &mut Lambda, n, pCD)?;
    ok &= CTableCreate(heap, &mut zf_zeta, n, pCD)?;
    let mut pzb_rho_pointer =
        crate::source::base::util::inchi_calloc::<ConTable>(heap, 1, size_of::<ConTable>() as u64)
            .unwrap_or_else(|_| SourceMutPointer::null());
    let mut pzb_rho_storage = if pzb_rho_pointer.is_null() {
        ok = 0;
        None
    } else {
        // SAFETY: this allocation is distinct from every earlier CanonGraph
        // allocation and is retained until the source cleanup or output
        // ownership-transfer branch.
        let mut storage = unsafe { heap.stable_slice_mut(pzb_rho_pointer)? };
        ok &= CTableCreate(heap, storage.get_mut(0)?, n, pCD)?;
        Some(storage)
    };
    let pzb_rho_is_usable = pzb_rho_storage.as_ref().is_some_and(|storage| {
        storage
            .get(0)
            .is_ok_and(|table| con_table_is_usable(table, pCD))
    });
    if !con_table_is_usable(&Lambda, pCD)
        || !con_table_is_usable(&zf_zeta, pCD)
        || (!pzb_rho_pointer.is_null() && !pzb_rho_is_usable)
    {
        ok = 0;
    }
    let max_set_size = if LargeMolecules == 0 {
        NORMALLY_ALLOWED_MAX_SET_SIZE as i32
    } else {
        MAX_SET_SIZE as i32
    };
    ok &= NodeSetCreate(heap, pCG, &mut Omega, n_tg, max_set_size)?;
    ok &= NodeSetCreate(heap, pCG, &mut Phi, n_tg, max_set_size)?;
    ok &= NodeSetCreate(heap, pCG, &mut cur_nodes, n_tg, 1)?;
    ok &= PartitionCreate(heap, &mut zeta, n_tg)?;
    ok &= PartitionCreate(heap, &mut rho, n_tg)?;
    ok &= TranspositionCreate(heap, &mut gamma, n_tg)?;

    let execution = (|| -> Result<i32, SourceHeapError> {
        let nNumLayers = i32::from(!pCD.NumH.is_null())
            + i32::from(!pCD.NumHfixed.is_null())
            + i32::from(!pCD.iso_sort_key.is_null());
        let dig = i32::from(bDigraph != 0 || nNumLayers != 0);
        let bSplitTautCompare = dig;
        let mut bZetaEqRho = 1_i32;
        let mut bZetaIsomorph = 0_i32;
        let nOneAdditionalLayer = GetOneAdditionalLayer(Some(pCD), pzb_rho_fix);
        let mut ret = 0_i32;
        let (mut k, mut index, mut l) = (1_i32, 0_i32, 0_i32);
        let mut t_Lemma;
        let (mut t_eq_zeta, mut h_zeta, mut h_rho, mut hz_rho, mut hz_zeta) = (0, 0, 0, 0, 0);
        let mut qzb_rho = 0_i32;
        let mut size = 1.0_f64;
        let mut r = 0_i32;
        let (mut tvc, mut tvh) = (0 as Node, 0 as Node);
        let mut nNumFoundGenerators = 0_i32;
        let (mut qzb_rho_fix, mut hzb_rho_fix, mut bRhoIsDiscrete) = (0_i32, 1_i32, 1_i32);
        let mut pzb_rho_fix_reached = 0_i32;
        let (mut L_rho_fix_prev, mut I_rho_fix_prev, mut k_rho_fix_prev) = (0_i32, -1_i32, 0_i32);
        if ok == 0 {
            return Ok(CT_OUT_OF_RAM);
        }
        let pzb_rho = pzb_rho_storage
            .as_mut()
            .ok_or(SourceHeapError::NullPointer)?
            .get_mut(0)?;
        let W = W.as_mut().ok_or(SourceHeapError::NullPointer)?;
        let W_len = W.len();
        UnorderedPartitionMakeDiscrete(heap, &theta, n_tg)?;
        t_Lemma = 2;
        pCC.lNumBreakTies = 0;
        pCC.lNumDecreasedCT = 0;
        pCC.lNumRejectedCT = 0;
        pCC.lNumEqualCT = 1;
        pCC.lNumTotCT = 0;
        // INCHI✔️✔️:     memset( kLeast_rho, 0, sizeof( kLeast_rho ) );
        let kLeast_rho = kLeast_rho_storage.write(std::array::from_fn(|_| kLeast::default()));
        // INCHI✔️✔️:     memset( kLeast_rho_fix, 0, sizeof( kLeast_rho_fix ) );
        let kLeast_rho_fix =
            kLeast_rho_fix_storage.write(std::array::from_fn(|_| kLeast::default()));
        let mut label = CanonGraphLabel::InitialDescent;
        'machine: loop {
            match label {
                CanonGraphLabel::InitialDescent => {
                    if PartitionIsDiscrete(heap, partition_at(pi, k)?, n_tg)? != 0 {
                        PartitionCopy(heap, &rho, partition_at(pi, k)?, n_tg)?;
                        CtPartFill(
                            heap,
                            G,
                            pCD,
                            partition_at(pi, k)?,
                            Some(&mut *pzb_rho),
                            1,
                            n,
                            n_tg,
                            n_max,
                        )?;
                        CtPartINCHI_CANON_INFINITY(
                            heap,
                            &mut *pzb_rho,
                            SourceMutPointer::null(),
                            2,
                        )?;
                        pCC.lNumTotCT += 1;
                        label = CanonGraphLabel::ExitFunction;
                        continue;
                    }
                    if dig == 0 && PartitionSatisfiesLemma_2_25(heap, partition_at(pi, 1)?, n)? != 0
                    {
                        t_Lemma = 1;
                    }
                    CtPartClear(heap, &mut Lambda, 1)?;
                    while k != 0 {
                        PartitionGetFirstCell(
                            heap,
                            partition_at(pi, k)?,
                            W.prefix_mut(W_len)?,
                            k,
                            n,
                        )?;
                        let wi = usize::try_from(k - 1)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let minimum =
                            CellGetMinNode(heap, partition_at(pi, k)?, W.get(wi)?, 0, Some(pCD))?;
                        source_set(heap, v, k - 1, minimum)?;
                        source_set(heap, e, k - 1, 0)?;
                        if dig != 0
                            || PartitionSatisfiesLemma_2_25(heap, partition_at(pi, k)?, n)? == 0
                        {
                            t_Lemma = k + 1;
                        }
                        let vv = source_get(heap, v, k - 1)?;
                        let result = PartitionColorVertex(
                            heap,
                            pCG,
                            G,
                            pi.get_mut(wi..)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?,
                            vv,
                            n,
                            n_tg,
                            n_max,
                            bDigraph,
                            0,
                        )?;
                        if result < 0 {
                            ret = result;
                            break 'machine;
                        }
                        pCC.lNumBreakTies += 1;
                        k += 1;
                        CtPartFill(
                            heap,
                            G,
                            pCD,
                            partition_at(pi, k)?,
                            Some(&mut Lambda),
                            k - 1,
                            n,
                            n_tg,
                            n_max,
                        )?;
                        if let Some(fixed) = pzb_rho_fix.filter(|_| qzb_rho_fix <= 0) {
                            qzb_rho_fix = CtPartCompare(
                                heap,
                                &Lambda,
                                fixed,
                                SourceMutPointer::null(),
                                Some(kLeast_rho_fix),
                                k - 1,
                                1,
                                bSplitTautCompare,
                            )?;
                            if qzb_rho_fix <= 0 {
                                hzb_rho_fix = k;
                            }
                        }
                        if qzb_rho_fix <= 0 {
                            CtPartCopy(heap, &mut *pzb_rho, &Lambda, k - 1)?;
                        }
                        CtPartCopy(heap, &mut zf_zeta, &Lambda, k - 1)?;
                        if PartitionIsDiscrete(heap, partition_at(pi, k)?, n)? != 0 {
                            break;
                        }
                    }
                    pCC.lNumTotCT += 1;
                    h_zeta = k;
                    t_eq_zeta = k;
                    hz_zeta = k;
                    CtPartINCHI_CANON_INFINITY(heap, &mut zf_zeta, SourceMutPointer::null(), k)?;
                    PartitionCopy(heap, &zeta, partition_at(pi, k)?, n_tg)?;
                    if let Some(fixed) = pzb_rho_fix {
                        if qzb_rho_fix == 0 {
                            qzb_rho_fix =
                                CtFullCompare(heap, &Lambda, fixed, 1, bSplitTautCompare)?;
                            if qzb_rho_fix > 0 {
                                hzb_rho_fix = 1;
                            }
                        }
                        if hzb_rho_fix > 1 {
                            PartitionCopy(heap, &rho, partition_at(pi, hzb_rho_fix)?, n_tg)?;
                        }
                        hz_rho = hzb_rho_fix;
                        h_rho = hzb_rho_fix;
                        bRhoIsDiscrete = i32::from(hzb_rho_fix == k);
                        if bRhoIsDiscrete != 0 {
                            CtPartINCHI_CANON_INFINITY(
                                heap,
                                &mut *pzb_rho,
                                SourceMutPointer::null(),
                                k,
                            )?;
                            pzb_rho_fix_reached = i32::from(qzb_rho_fix == 0);
                            CtCompareLayersGetFirstDiff(
                                Some(kLeast_rho_fix),
                                nOneAdditionalLayer,
                                &mut L_rho_fix_prev,
                                &mut I_rho_fix_prev,
                                &mut k_rho_fix_prev,
                            )?;
                        }
                    } else {
                        PartitionCopy(heap, &rho, partition_at(pi, k)?, n_tg)?;
                        hz_rho = k;
                        h_rho = k;
                        CtPartINCHI_CANON_INFINITY(
                            heap,
                            &mut *pzb_rho,
                            SourceMutPointer::null(),
                            k,
                        )?;
                    }
                    qzb_rho = 0;
                    r = k;
                    source_set(heap, v, k - 1, INCHI_CANON_INFINITY as Node)?;
                    CellMakeEmpty(W.prefix_mut(W_len)?, k)?;
                    k -= 1;
                    label = CanonGraphLabel::L13;
                }
                CanonGraphLabel::L2 => {
                    let vv = source_get(heap, v, k - 1)?;
                    let wi =
                        usize::try_from(k - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let result = PartitionColorVertex(
                        heap,
                        pCG,
                        G,
                        pi.get_mut(wi..)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                        vv,
                        n,
                        n_tg,
                        n_max,
                        bDigraph,
                        0,
                    )?;
                    if result < 0 {
                        ret = result;
                        break 'machine;
                    }
                    pCC.lNumBreakTies += 1;
                    k += 1;
                    CtPartFill(
                        heap,
                        G,
                        pCD,
                        partition_at(pi, k)?,
                        Some(&mut Lambda),
                        k - 1,
                        n,
                        n_tg,
                        n_max,
                    )?;
                    source_set(heap, e, k - 1, 0)?;
                    source_set(heap, v, k - 1, INCHI_CANON_INFINITY as Node)?;
                    CellMakeEmpty(W.prefix_mut(W_len)?, k)?;
                    if hz_zeta == k - 1
                        && CtPartCompare(
                            heap,
                            &Lambda,
                            &zf_zeta,
                            SourceMutPointer::null(),
                            None,
                            k - 1,
                            0,
                            bSplitTautCompare,
                        )? == 0
                    {
                        hz_zeta = k;
                    }
                    if let Some(fixed) = pzb_rho_fix.filter(|_| qzb_rho_fix == 0) {
                        qzb_rho_fix = CtPartCompare(
                            heap,
                            &Lambda,
                            fixed,
                            SourceMutPointer::null(),
                            Some(kLeast_rho_fix),
                            k - 1,
                            1,
                            bSplitTautCompare,
                        )?;
                        if qzb_rho_fix == 0 && bRhoIsDiscrete != 0 {
                            qzb_rho_fix = CtPartCompareLayers(
                                Some(kLeast_rho_fix),
                                L_rho_fix_prev,
                                nOneAdditionalLayer,
                            )?;
                            if qzb_rho_fix != 0 {
                                let layer = qzb_rho_fix.abs() - 1;
                                let value = kLeast_rho_fix
                                    .get(
                                        usize::try_from(layer)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                if layer < L_rho_fix_prev
                                    || layer
                                        == i32::from(
                                            L_rho_fix_prev != 0 && value.i < I_rho_fix_prev,
                                        )
                                {
                                    qzb_rho_fix = layer + 1;
                                }
                            }
                        }
                        if qzb_rho_fix > 0 {
                            pCC.lNumRejectedCT += 1;
                        }
                    }
                    if pzb_rho_fix.is_some() && qzb_rho_fix <= 0 {
                        hzb_rho_fix = k;
                    }
                    if qzb_rho_fix <= 0 {
                        if hz_rho == k - 1 && qzb_rho == 0 && bRhoIsDiscrete != 0 {
                            qzb_rho = CtPartCompare(
                                heap,
                                &Lambda,
                                &*pzb_rho,
                                SourceMutPointer::null(),
                                Some(kLeast_rho),
                                k - 1,
                                0,
                                bSplitTautCompare,
                            )?;
                            if qzb_rho == 0 && bRhoIsDiscrete != 0 {
                                qzb_rho = CtPartCompareLayers(Some(kLeast_rho), L_rho_fix_prev, 0)?;
                                if qzb_rho != 0 {
                                    let layer = qzb_rho.abs() - 1;
                                    let value = kLeast_rho
                                        .get(
                                            usize::try_from(layer)
                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                        )
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    if layer < L_rho_fix_prev
                                        || layer
                                            == i32::from(
                                                L_rho_fix_prev != 0 && value.i < I_rho_fix_prev,
                                            )
                                    {
                                        qzb_rho = -(layer + 1);
                                    }
                                }
                            }
                            if qzb_rho == 0 {
                                hz_rho = k;
                            } else if qzb_rho < 0 {
                                pCC.lNumRejectedCT += 1;
                            }
                        }
                        if (qzb_rho > 0 || (qzb_rho == 0 && bRhoIsDiscrete == 0)) && nNumLayers == 0
                        {
                            CtPartCopy(heap, &mut *pzb_rho, &Lambda, k - 1)?;
                        }
                    }
                    if hz_zeta == k || hz_rho == k || (qzb_rho >= 0 && qzb_rho_fix <= 0) {
                        if PartitionIsDiscrete(heap, partition_at(pi, k)?, n)? != 0 {
                            pCC.lNumTotCT += 1;
                            label = CanonGraphLabel::L7;
                        } else {
                            PartitionGetFirstCell(
                                heap,
                                partition_at(pi, k)?,
                                W.prefix_mut(W_len)?,
                                k,
                                n,
                            )?;
                            let wi = usize::try_from(k - 1)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let minimum = CellGetMinNode(
                                heap,
                                partition_at(pi, k)?,
                                W.get(wi)?,
                                0,
                                Some(pCD),
                            )?;
                            source_set(heap, v, k - 1, minimum)?;
                            if dig != 0
                                || PartitionSatisfiesLemma_2_25(heap, partition_at(pi, k)?, n)? == 0
                            {
                                t_Lemma = k + 1;
                            }
                            source_set(heap, e, k - 1, 0)?;
                            label = CanonGraphLabel::L2;
                        }
                    } else {
                        label = CanonGraphLabel::L6;
                    }
                }
                CanonGraphLabel::L6 => {
                    let k2 = k;
                    k = (t_Lemma - 1).min((t_eq_zeta - 1).max(hz_rho));
                    if k2 == t_Lemma || dig != 0 {
                        label = CanonGraphLabel::L13;
                    } else {
                        l = (l + 1).min(max_set_size);
                        PartitionGetMcrAndFixSet(
                            heap,
                            pCG,
                            partition_at(pi, t_Lemma)?,
                            &Omega,
                            &Phi,
                            n_tg,
                            l,
                        )?;
                        label = CanonGraphLabel::L12;
                    }
                }
                CanonGraphLabel::L7 => {
                    if h_zeta == 0 {
                        ret = CT_CANON_ERR;
                        break 'machine;
                    }
                    if k != hz_zeta {
                        label = CanonGraphLabel::L8;
                        continue;
                    }
                    if CtFullCompare(heap, &Lambda, &zf_zeta, 0, bSplitTautCompare)? == 0 {
                        PartitionGetTransposition(heap, &zeta, partition_at(pi, k)?, n_tg, &gamma)?;
                        bZetaIsomorph = 1;
                        label = CanonGraphLabel::L10;
                    } else if nNumLayers == 0 {
                        ret = -2;
                        break 'machine;
                    } else {
                        label = CanonGraphLabel::L8;
                    }
                }
                CanonGraphLabel::L8 => {
                    if qzb_rho < 0 && (pzb_rho_fix.is_none() || qzb_rho_fix > 0) {
                        label = CanonGraphLabel::L6;
                        continue;
                    }
                    if let Some(fixed) = pzb_rho_fix.filter(|_| qzb_rho_fix == 0) {
                        qzb_rho_fix = CtFullCompareLayers(kLeast_rho_fix)?;
                        let alt = CtFullCompare(heap, &Lambda, fixed, 1, bSplitTautCompare)?;
                        if qzb_rho_fix != alt {
                            qzb_rho_fix = alt;
                        }
                        if pzb_rho_fix_reached == 0 {
                            pzb_rho_fix_reached = i32::from(qzb_rho_fix == 0);
                        }
                        if qzb_rho_fix > 0 {
                            qzb_rho_fix = 0;
                            label = CanonGraphLabel::L6;
                            continue;
                        }
                        qzb_rho_fix = 0;
                    }
                    if qzb_rho < 0 {
                        label = CanonGraphLabel::L6;
                    } else if qzb_rho > 0 || bRhoIsDiscrete == 0 || k < r {
                        label = CanonGraphLabel::L9;
                    } else {
                        qzb_rho = CtFullCompareLayers(kLeast_rho)?;
                        let alt = CtFullCompare(heap, &Lambda, &*pzb_rho, 0, bSplitTautCompare)?;
                        if qzb_rho != alt {
                            qzb_rho = alt;
                        }
                        if qzb_rho > 0 {
                            label = CanonGraphLabel::L9;
                        } else if qzb_rho < 0 {
                            qzb_rho = 0;
                            label = CanonGraphLabel::L6;
                        } else {
                            r = k;
                            PartitionGetTransposition(
                                heap,
                                partition_at(pi, k)?,
                                &rho,
                                n_tg,
                                &gamma,
                            )?;
                            bZetaIsomorph = 0;
                            pCC.lNumEqualCT += 1;
                            label = CanonGraphLabel::L10;
                        }
                    }
                }
                CanonGraphLabel::L9 => {
                    PartitionCopy(heap, &rho, partition_at(pi, k)?, n_tg)?;
                    if nNumLayers != 0 {
                        CtFullCopy(heap, &mut *pzb_rho, &Lambda)?;
                    }
                    bZetaEqRho = 0;
                    qzb_rho = 0;
                    CtCompareLayersGetFirstDiff(
                        Some(kLeast_rho_fix),
                        nOneAdditionalLayer,
                        &mut L_rho_fix_prev,
                        &mut I_rho_fix_prev,
                        &mut k_rho_fix_prev,
                    )?;
                    kLeast_rho.fill(kLeast::default());
                    h_rho = k;
                    hz_rho = k;
                    CtPartINCHI_CANON_INFINITY(heap, &mut *pzb_rho, SourceMutPointer::null(), k)?;
                    pCC.lNumDecreasedCT += 1;
                    pCC.lNumEqualCT = 1;
                    bRhoIsDiscrete = 1;
                    label = CanonGraphLabel::L6;
                }
                CanonGraphLabel::L10 => {
                    pCC.lNumEqualCT +=
                        i64::from(bZetaEqRho != 0 || !(bZetaIsomorph != 0 || qzb_rho != 0));
                    l = (l + 1).min(max_set_size);
                    TranspositionGetMcrAndFixSetAndUnorderedPartition(
                        heap,
                        pCG,
                        &mut gamma,
                        &Omega,
                        &Phi,
                        n_tg,
                        l,
                        &mut theta_from_gamma,
                    )?;
                    if UnorderedPartitionJoin(heap, &theta_from_gamma, &theta, n_tg)? == 0 {
                        label = CanonGraphLabel::L11;
                    } else {
                        nNumFoundGenerators += 1;
                        if tvc == GetUnorderedPartitionMcrNode(heap, &theta, tvc)? {
                            label = CanonGraphLabel::L11;
                        } else {
                            k = h_zeta;
                            label = CanonGraphLabel::L13;
                        }
                    }
                }
                CanonGraphLabel::L11 => {
                    k = h_rho;
                    label = CanonGraphLabel::L12;
                }
                CanonGraphLabel::L12 => {
                    if source_get(heap, e, k - 1)? == 1
                        && source_get(heap, v, k - 1)? != INCHI_CANON_INFINITY as Node
                    {
                        let wi = usize::try_from(k - 1)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        CellIntersectWithSet(
                            heap,
                            pCG,
                            partition_at(pi, k)?,
                            W.get(wi)?,
                            &Omega,
                            l,
                        )?;
                    }
                    label = CanonGraphLabel::L13;
                }
                CanonGraphLabel::L13 => {
                    if user_action.is_some_and(|callback| callback() == USER_ACTION_QUIT)
                        || console_quit.is_some_and(|callback| callback() != 0)
                    {
                        ret = CT_USER_QUIT_ERR;
                        break 'machine;
                    }
                    let timeout: Option<inchiTime> = if pCD.ulTimeOutTime.is_null() {
                        None
                    } else {
                        Some(source_clone(heap, pCD.ulTimeOutTime, 0)?)
                    };
                    if crate::source::base::ichicano::bInchiTimeIsOver(
                        ic,
                        timeout.as_ref(),
                        clock_result,
                    ) != 0
                    {
                        ret = CT_TIMEOUT_ERR;
                        break 'machine;
                    }
                    if k == 0 {
                        label = CanonGraphLabel::ExitFunction;
                        continue;
                    }
                    if k < h_rho {
                        h_rho = k;
                    }
                    if k > h_zeta {
                        if source_get(heap, v, k - 1)? == INCHI_CANON_INFINITY as Node {
                            k -= 1;
                        } else {
                            label = CanonGraphLabel::L17;
                        }
                        continue;
                    }
                    if k == h_zeta {
                        label = CanonGraphLabel::L14;
                    } else {
                        h_zeta = k;
                        let wi = usize::try_from(k - 1)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        tvh = CellGetMinNode(heap, partition_at(pi, k)?, W.get(wi)?, 0, Some(pCD))?;
                        tvc = tvh;
                        label = CanonGraphLabel::L14;
                    }
                }
                CanonGraphLabel::L14 => {
                    let current = source_get(heap, v, k - 1)?;
                    if GetUnorderedPartitionMcrNode(heap, &theta, current)?
                        == GetUnorderedPartitionMcrNode(heap, &theta, tvh)?
                    {
                        index += 1;
                    }
                    let wi =
                        usize::try_from(k - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let minimum =
                        CellGetMinNode(heap, partition_at(pi, k)?, W.get(wi)?, current, Some(pCD))?;
                    source_set(heap, v, k - 1, minimum)?;
                    if minimum == INCHI_CANON_INFINITY as Node {
                        label = CanonGraphLabel::L16;
                    } else if minimum != GetUnorderedPartitionMcrNode(heap, &theta, minimum)? {
                        label = CanonGraphLabel::L14;
                    } else {
                        label = CanonGraphLabel::L15;
                    }
                }
                CanonGraphLabel::L15 => {
                    t_Lemma = t_Lemma.min(k + 1);
                    hz_zeta = hz_zeta.min(k);
                    if hz_rho >= k {
                        qzb_rho = 0;
                    }
                    if hz_rho > k {
                        hz_rho = k;
                    }
                    UpdateCompareLayers(Some(kLeast_rho), hz_rho)?;
                    if pzb_rho_fix.is_some() {
                        if hzb_rho_fix >= k {
                            qzb_rho_fix = 0;
                        }
                        if hzb_rho_fix > k {
                            hzb_rho_fix = k;
                        }
                        UpdateCompareLayers(Some(kLeast_rho_fix), hzb_rho_fix)?;
                    }
                    label = CanonGraphLabel::L2;
                }
                CanonGraphLabel::L16 => {
                    let wi =
                        usize::try_from(k - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    if t_eq_zeta == k + 1
                        && index == CellGetNumberOfNodes(heap, partition_at(pi, k)?, W.get(wi)?)?
                    {
                        t_eq_zeta = k;
                    }
                    size *= f64::from(index);
                    index = 0;
                    k -= 1;
                    label = CanonGraphLabel::L13;
                }
                CanonGraphLabel::L17 => {
                    if source_get(heap, e, k - 1)? == 0
                        && source_get(heap, v, k - 1)? != INCHI_CANON_INFINITY as Node
                    {
                        NodeSetFromVertices(heap, pCG, &cur_nodes, 1, v, k - 1)?;
                        let wi = usize::try_from(k - 1)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let mut set_index = 1_i32;
                        while set_index <= l {
                            if AllNodesAreInSet(heap, &cur_nodes, 1, &Phi, set_index)? != 0 {
                                CellIntersectWithSet(
                                    heap,
                                    pCG,
                                    partition_at(pi, k)?,
                                    W.get(wi)?,
                                    &Omega,
                                    set_index,
                                )?;
                            }
                            set_index += 1;
                        }
                    }
                    source_set(heap, e, k - 1, 1)?;
                    let wi =
                        usize::try_from(k - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let current = source_get(heap, v, k - 1)?;
                    let minimum =
                        CellGetMinNode(heap, partition_at(pi, k)?, W.get(wi)?, current, Some(pCD))?;
                    source_set(heap, v, k - 1, minimum)?;
                    if minimum != INCHI_CANON_INFINITY as Node {
                        label = CanonGraphLabel::L15;
                    } else {
                        k -= 1;
                        label = CanonGraphLabel::L13;
                    }
                }
                CanonGraphLabel::ExitFunction => {
                    if bRhoIsDiscrete == 0 {
                        ret = CT_CANON_ERR;
                        break 'machine;
                    }
                    if let Some(fixed) = pzb_rho_fix {
                        if CtFullCompare(heap, fixed, &*pzb_rho, 1, bSplitTautCompare)? != 0 {
                            ret = CT_CANON_ERR;
                            break 'machine;
                        }
                    }
                    let count =
                        usize::try_from(n_tg).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    heap.slice_mut(nSymmRank)?
                        .get_mut(..count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .fill(0);
                    let mut i = 0_i32;
                    while i < n_tg {
                        let atom = i32::from(source_get(heap, rho.AtNumber, i)?);
                        let representative = i32::from(GetUnorderedPartitionMcrNode(
                            heap,
                            &theta,
                            (atom + 1) as AT_NUMB,
                        )?) - 1;
                        let old = source_get(heap, nSymmRank, representative)?;
                        let rank = source_get(heap, rho.Rank, atom)?;
                        if old == 0 || old > rank {
                            source_set(heap, nSymmRank, representative, rank)?;
                        }
                        i += 1;
                    }
                    i = 0;
                    while i < n_tg {
                        let atom = i32::from(source_get(heap, rho.AtNumber, i)?);
                        let representative = i32::from(GetUnorderedPartitionMcrNode(
                            heap,
                            &theta,
                            (atom + 1) as AT_NUMB,
                        )?) - 1;
                        source_set(
                            heap,
                            nSymmRank,
                            atom,
                            source_get(heap, nSymmRank, representative)?,
                        )?;
                        i += 1;
                    }
                    i = 0;
                    while i < n_tg {
                        source_set(heap, nCanonRank, i, source_get(heap, rho.Rank, i)?)?;
                        source_set(
                            heap,
                            nAtomNumberCanon,
                            i,
                            source_get(heap, rho.AtNumber, i)?,
                        )?;
                        i += 1;
                    }
                    pCD.nLenLinearCT = pzb_rho.lenCt - 1;
                    if !pCD.LinearCT.is_null() && pCD.nLenLinearCT > 0 {
                        i = 0;
                        while i < pCD.nLenLinearCT {
                            source_set(heap, pCD.LinearCT, i, source_get(heap, pzb_rho.Ctbl, i)?)?;
                            i += 1;
                        }
                    }
                    pCC.lNumTotCT = pCC.lNumDecreasedCT + pCC.lNumRejectedCT + pCC.lNumEqualCT;
                    pCC.dGroupSize = size;
                    pCC.lNumGenerators = i64::from(nNumFoundGenerators);
                    pCC.lNumStoredIsomorphisms = i64::from(l);
                    if let Some(output) = pp_zb_rho_out.as_deref_mut() {
                        if output.is_null() {
                            *output = pzb_rho_pointer;
                            pzb_rho_pointer = SourceMutPointer::null();
                        }
                    }
                    break 'machine;
                }
            }
        }
        Ok(ret)
    })();
    let mut cleanup_error = None;
    macro_rules! cleanup {
        ($expression:expr) => {
            if let Err(error) = $expression {
                if cleanup_error.is_none() {
                    cleanup_error = Some(error);
                }
            }
        };
    }
    cleanup!(UnorderedPartitionFree(heap, &mut theta));
    cleanup!(UnorderedPartitionFree(heap, &mut theta_from_gamma));
    if !W_pointer.is_null() {
        cleanup!(inchi_free(heap, W_pointer));
    }
    if !v.is_null() {
        cleanup!(inchi_free(heap, v));
    }
    if !e.is_null() {
        cleanup!(inchi_free(heap, e));
    }
    cleanup!(CTableFree(heap, Some(&mut Lambda)));
    cleanup!(CTableFree(heap, Some(&mut zf_zeta)));
    if !pzb_rho_pointer.is_null() {
        if let Some(storage) = pzb_rho_storage.as_mut() {
            // SAFETY: the source allocation contains exactly one ConTable and
            // ownership was not transferred when pzb_rho_pointer is non-null.
            let pzb_rho = unsafe { storage.get_unchecked_mut(0) };
            cleanup!(CTableFree(heap, Some(pzb_rho)));
        }
        cleanup!(inchi_free(heap, pzb_rho_pointer));
    }
    cleanup!(NodeSetFree(heap, pCG, &mut Omega));
    cleanup!(NodeSetFree(heap, pCG, &mut Phi));
    cleanup!(NodeSetFree(heap, pCG, &mut cur_nodes));
    cleanup!(PartitionFree(heap, Some(&mut zeta)));
    cleanup!(PartitionFree(heap, Some(&mut rho)));
    cleanup!(TranspositionFree(heap, Some(&mut gamma)));
    match execution {
        Err(error) => Err(error),
        Ok(_) if cleanup_error.is_some() => Err(cleanup_error.expect("checked")),
        Ok(result) => Ok(result),
    }
}

#[allow(non_snake_case)]
#[allow(non_snake_case)]
pub(crate) fn CTableCreate(
    heap: &mut SourceHeap,
    Ct: &mut ConTable,
    n: i32,
    pCD: &CANON_DATA,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:769 CTableCreate
    // INCHI✔️❌: int CTableCreate( ConTable *Ct, int n, CANON_DATA *pCD )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int maxlenCt = pCD->nMaxLenLinearCT + 1; /* add one element for CtPartINCHI_CANON_INFINITY() */
    // INCHI✔️❌:     int maxlenNumH = pCD->NumH ? ( pCD->maxlenNumH + 1 ) : 0;
    // INCHI✔️❌:     int maxlenNumHfixed = pCD->NumHfixed ? ( pCD->maxlenNumHfixed + 1 ) : 0;
    // INCHI✔️❌:     int maxlenIso = pCD->maxlen_iso_sort_key ? ( pCD->maxlen_iso_sort_key + 1 ) : 0;
    // INCHI✔️❌:     int maxlenIsoExchg = pCD->iso_exchg_atnos ? ( pCD->maxlen_iso_exchg_atnos + 1 ) : 0;
    // INCHI✔️❌:
    // INCHI❌❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI❌❌:     int maxlenIsoHfixed = pCD->maxlen_iso_sort_key_Hfixed ? ( pCD->maxlen_iso_sort_key_Hfixed + 1 ) : 0;
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     memset( Ct, 0, sizeof( Ct[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:     Ct->maxVert = n;
    // INCHI✔️❌:
    // INCHI✔️❌:     n++;
    // INCHI✔️❌:
    // INCHI✔️❌:     Ct->Ctbl = (AT_RANK*) inchi_calloc( maxlenCt, sizeof( Ct->Ctbl[0] ) );
    // INCHI✔️❌:     Ct->nextCtblPos = (AT_NUMB*) inchi_calloc( n, sizeof( Ct->nextCtblPos[0] ) );
    // INCHI✔️❌:     Ct->nextAtRank = (AT_RANK*) inchi_calloc( n, sizeof( Ct->nextAtRank[0] ) );
    // INCHI✔️❌:     if (maxlenNumH)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Ct->NumH = (NUM_H *) inchi_calloc( maxlenNumH, sizeof( Ct->NumH[0] ) );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (maxlenNumHfixed)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Ct->NumHfixed = (NUM_H *) inchi_calloc( maxlenNumHfixed, sizeof( Ct->NumH[0] ) );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (maxlenIso)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Ct->iso_sort_key = (AT_ISO_SORT_KEY *) inchi_calloc( maxlenIso, sizeof( Ct->iso_sort_key[0] ) );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (maxlenIsoExchg)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Ct->iso_exchg_atnos = (S_CHAR *) inchi_calloc( maxlenIsoExchg, sizeof( Ct->iso_exchg_atnos[0] ) );
    // INCHI✔️❌:     }
    // INCHI❌❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI❌❌:     if (maxlenIsoHfixed)
    // INCHI❌❌:     {
    // INCHI❌❌:         Ct->iso_sort_key_Hfixed = (AT_ISO_SORT_KEY *) inchi_calloc( maxlenIsoHfixed, sizeof( Ct->iso_sort_key_Hfixed[0] ) );
    // INCHI❌❌:     }
    // INCHI❌❌: #endif
    // INCHI❌❌: #ifdef INCHI_CANON_USE_HASH
    // INCHI❌❌:     Ct->hash = (CtHash*) inchi_calloc( n, sizeof( Ct->hash[0] ) );
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     Ct->lenCt = 0;
    // INCHI✔️❌:     Ct->nLenCTAtOnly = pCD->nLenCTAtOnly;
    // INCHI✔️❌:     Ct->maxlenCt = maxlenCt;
    // INCHI✔️❌:     Ct->lenNumH = 0;
    // INCHI✔️❌:     Ct->maxlenNumH = maxlenNumH;
    // INCHI✔️❌:     Ct->len_iso_sort_key = 0;
    // INCHI✔️❌:     Ct->maxlen_iso_sort_key = maxlenIso;
    // INCHI✔️❌:     Ct->len_iso_exchg_atnos = 0;
    // INCHI✔️❌:     Ct->maxlen_iso_exchg_atnos = maxlenIso;
    // INCHI✔️❌:
    // INCHI❌❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI❌❌:     Ct->len_iso_sort_key_Hfixed = 0;
    // INCHI❌❌:     Ct->maxlen_iso_sort_key_Hfixed = maxlenIsoHfixed;
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     Ct->maxPos = n;
    // INCHI✔️❌:     Ct->lenPos = 0;
    // INCHI✔️❌:     /* djb-rwth: fixing a NULL pointer dereferences */
    // INCHI✔️❌:     if (Ct->nextAtRank)
    // INCHI✔️❌:         Ct->nextAtRank[0] = 0;
    // INCHI✔️❌:     if (Ct->nextCtblPos)
    // INCHI✔️❌:         Ct->nextCtblPos[0] = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (Ct->Ctbl && Ct->nextCtblPos &&
    // INCHI✔️❌:         ( !maxlenNumH || Ct->NumH ) &&
    // INCHI✔️❌:         ( !maxlenNumHfixed || Ct->NumHfixed ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CTableCreate
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CTableCreate
    // INCHI✔️❌: #define USE_ISO_SORT_KEY_HFIXED  0
    // INCHI✔️❌: /* #define INCHI_CANON_USE_HASH */
    // END INCHI ACTIVE MACRO CONFIGURATION: CTableCreate

    let max_len_ct = pCD.nMaxLenLinearCT.wrapping_add(1);
    let max_len_h = if !pCD.NumH.is_null() {
        pCD.maxlenNumH.wrapping_add(1)
    } else {
        0
    };
    let max_len_fixed_h = if !pCD.NumHfixed.is_null() {
        pCD.maxlenNumHfixed.wrapping_add(1)
    } else {
        0
    };
    let max_len_iso = if pCD.maxlen_iso_sort_key != 0 {
        pCD.maxlen_iso_sort_key.wrapping_add(1)
    } else {
        0
    };
    let max_len_iso_exchange = if !pCD.iso_exchg_atnos.is_null() {
        pCD.maxlen_iso_exchg_atnos.wrapping_add(1)
    } else {
        0
    };

    *Ct = ConTable::default();
    Ct.maxVert = n;
    let allocation_n = n.wrapping_add(1);

    Ct.Ctbl = crate::source::base::util::inchi_calloc::<AT_RANK>(
        heap,
        max_len_ct as u64,
        size_of::<AT_RANK>() as u64,
    )
    .unwrap_or_else(|_| SourceMutPointer::null());
    Ct.nextCtblPos = crate::source::base::util::inchi_calloc::<AT_NUMB>(
        heap,
        allocation_n as u64,
        size_of::<AT_NUMB>() as u64,
    )
    .unwrap_or_else(|_| SourceMutPointer::null());
    Ct.nextAtRank = crate::source::base::util::inchi_calloc::<AT_RANK>(
        heap,
        allocation_n as u64,
        size_of::<AT_RANK>() as u64,
    )
    .unwrap_or_else(|_| SourceMutPointer::null());
    if max_len_h != 0 {
        Ct.NumH = crate::source::base::util::inchi_calloc(
            heap,
            max_len_h as u64,
            size_of::<NUM_H>() as u64,
        )
        .unwrap_or_else(|_| SourceMutPointer::null());
    }
    if max_len_fixed_h != 0 {
        Ct.NumHfixed = crate::source::base::util::inchi_calloc(
            heap,
            max_len_fixed_h as u64,
            size_of::<NUM_H>() as u64,
        )
        .unwrap_or_else(|_| SourceMutPointer::null());
    }
    if max_len_iso != 0 {
        Ct.iso_sort_key = crate::source::base::util::inchi_calloc(
            heap,
            max_len_iso as u64,
            size_of::<AT_ISO_SORT_KEY>() as u64,
        )
        .unwrap_or_else(|_| SourceMutPointer::null());
    }
    if max_len_iso_exchange != 0 {
        Ct.iso_exchg_atnos = crate::source::base::util::inchi_calloc(
            heap,
            max_len_iso_exchange as u64,
            size_of::<i8>() as u64,
        )
        .unwrap_or_else(|_| SourceMutPointer::null());
    }

    Ct.lenCt = 0;
    Ct.nLenCTAtOnly = pCD.nLenCTAtOnly;
    Ct.maxlenCt = max_len_ct;
    Ct.lenNumH = 0;
    Ct.maxlenNumH = max_len_h;
    Ct.len_iso_sort_key = 0;
    Ct.maxlen_iso_sort_key = max_len_iso;
    Ct.len_iso_exchg_atnos = 0;
    Ct.maxlen_iso_exchg_atnos = max_len_iso;
    Ct.maxPos = allocation_n;
    Ct.lenPos = 0;
    if !Ct.nextAtRank.is_null() {
        source_set(heap, Ct.nextAtRank, 0, 0)?;
    }
    if !Ct.nextCtblPos.is_null() {
        source_set(heap, Ct.nextCtblPos, 0, 0)?;
    }

    Ok(
        if !Ct.Ctbl.is_null()
            && !Ct.nextCtblPos.is_null()
            && (max_len_h == 0 || !Ct.NumH.is_null())
            && (max_len_fixed_h == 0 || !Ct.NumHfixed.is_null())
        {
            1
        } else {
            0
        },
    )
}

#[allow(non_snake_case)]
#[allow(non_snake_case)]
pub(crate) fn CTableFree(
    heap: &mut SourceHeap,
    Ct: Option<&mut ConTable>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:851 CTableFree
    // INCHI✔️❌: void CTableFree( ConTable *Ct )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (Ct)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Ct->Ctbl)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( Ct->Ctbl );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (Ct->nextCtblPos)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( Ct->nextCtblPos );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (Ct->nextAtRank)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( Ct->nextAtRank );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (Ct->NumH)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( Ct->NumH );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (Ct->NumHfixed)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( Ct->NumHfixed );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (Ct->iso_sort_key)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( Ct->iso_sort_key );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (Ct->iso_exchg_atnos)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( Ct->iso_exchg_atnos );
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI❌❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI❌❌:         if (Ct->iso_sort_key_Hfixed)
    // INCHI❌❌:         {
    // INCHI❌❌:             inchi_free( Ct->iso_sort_key_Hfixed );
    // INCHI❌❌:         }
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI❌❌: #ifdef INCHI_CANON_USE_HASH
    // INCHI❌❌:         if (Ct->hash)
    // INCHI❌❌:         {
    // INCHI❌❌:             inchi_free( Ct->hash );
    // INCHI❌❌:         }
    // INCHI❌❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         memset( Ct, 0, sizeof( Ct[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     }
    // END INCHI C FUNCTION: CTableFree
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CTableFree
    // INCHI✔️❌: #define USE_ISO_SORT_KEY_HFIXED  0
    // INCHI✔️❌: /* #define INCHI_CANON_USE_HASH */
    // END INCHI ACTIVE MACRO CONFIGURATION: CTableFree

    let Some(Ct) = Ct else {
        return Ok(());
    };
    if !Ct.Ctbl.is_null() {
        inchi_free(heap, Ct.Ctbl)?;
    }
    if !Ct.nextCtblPos.is_null() {
        inchi_free(heap, Ct.nextCtblPos)?;
    }
    if !Ct.nextAtRank.is_null() {
        inchi_free(heap, Ct.nextAtRank)?;
    }
    if !Ct.NumH.is_null() {
        inchi_free(heap, Ct.NumH)?;
    }
    if !Ct.NumHfixed.is_null() {
        inchi_free(heap, Ct.NumHfixed)?;
    }
    if !Ct.iso_sort_key.is_null() {
        inchi_free(heap, Ct.iso_sort_key)?;
    }
    if !Ct.iso_exchg_atnos.is_null() {
        inchi_free(heap, Ct.iso_exchg_atnos)?;
    }
    *Ct = ConTable::default();
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn TranspositionCreate(
    heap: &mut SourceHeap,
    p: &mut Transposition,
    n: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:694 TranspositionCreate
    // INCHI✔️✔️: int TranspositionCreate( Transposition *p, int n )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     p->nAtNumb = (AT_NUMB*) inchi_calloc( n, sizeof( p->nAtNumb[0] ) );
    // INCHI✔️✔️:     if (p->nAtNumb)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: TranspositionCreate

    p.nAtNumb = match crate::source::base::util::inchi_calloc::<AT_NUMB>(
        heap,
        n as u64,
        size_of::<AT_NUMB>() as u64,
    ) {
        Ok(pointer) => pointer,
        Err(_) => SourceMutPointer::null(),
    };
    Ok(if p.nAtNumb.is_null() { 0 } else { 1 })
}

#[allow(non_snake_case)]
pub(crate) fn TranspositionFree(
    heap: &mut SourceHeap,
    p: Option<&mut Transposition>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:707 TranspositionFree
    // INCHI✔️✔️: void TranspositionFree( Transposition *p )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     if (p && p->nAtNumb)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inchi_free( p->nAtNumb );
    // INCHI✔️✔️:         p->nAtNumb = NULL;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: TranspositionFree

    let Some(p) = p else {
        return Ok(());
    };
    if !p.nAtNumb.is_null() {
        inchi_free(heap, p.nAtNumb)?;
        p.nAtNumb = SourceMutPointer::null();
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn UnorderedPartitionCreate(
    heap: &mut SourceHeap,
    p: &mut UnorderedPartition,
    n: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:904 UnorderedPartitionCreate
    // INCHI✔️✔️: int UnorderedPartitionCreate( UnorderedPartition *p, int n )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     p->equ2 = (AT_NUMB*) inchi_calloc( n, sizeof( p->equ2[0] ) );
    // INCHI✔️✔️:     /* p->next = (AT_NUMB*)inchi_calloc( n, sizeof(p->next[0])); */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (p->equ2 /*&& p->next*/)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: UnorderedPartitionCreate

    p.equ2 = crate::source::base::util::inchi_calloc::<AT_NUMB>(
        heap,
        n as u64,
        size_of::<AT_NUMB>() as u64,
    )
    .unwrap_or_else(|_| SourceMutPointer::null());
    Ok(if p.equ2.is_null() { 0 } else { 1 })
}

#[allow(non_snake_case)]
pub(crate) fn UnorderedPartitionFree(
    heap: &mut SourceHeap,
    p: &mut UnorderedPartition,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:919 UnorderedPartitionFree
    // INCHI✔️✔️: void UnorderedPartitionFree( UnorderedPartition *p )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     if (p->equ2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inchi_free( p->equ2 );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* if (p->next) inchi_free( p->next); */
    // INCHI✔️✔️:     p->equ2 = NULL;
    // INCHI✔️✔️:     /* p->next = NULL; */
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: UnorderedPartitionFree

    if !p.equ2.is_null() {
        inchi_free(heap, p.equ2)?;
    }
    p.equ2 = SourceMutPointer::null();
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn UnorderedPartitionMakeDiscrete(
    heap: &mut SourceHeap,
    p: &UnorderedPartition,
    n: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:933 UnorderedPartitionMakeDiscrete
    // INCHI✔️✔️: void UnorderedPartitionMakeDiscrete( UnorderedPartition *p, int n )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i;
    // INCHI✔️✔️:     for (i = 0; i < n; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         p->equ2[i] = (AT_NUMB) i;
    // INCHI✔️✔️:         /* p->next[i] = INCHI_CANON_INFINITY; */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     INCHI_HEAPCHK
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: UnorderedPartitionMakeDiscrete

    let mut i = 0_i32;
    while i < n {
        source_set(heap, p.equ2, i, i as AT_NUMB)?;
        i = i.wrapping_add(1);
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn PartitionCreate(
    heap: &mut SourceHeap,
    p: &mut Partition,
    n: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:946 PartitionCreate
    // INCHI✔️✔️: int PartitionCreate( Partition *p, int n )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     p->AtNumber = (AT_NUMB*) inchi_calloc( n, sizeof( p->AtNumber[0] ) );
    // INCHI✔️✔️:     p->Rank = (AT_RANK*) inchi_calloc( n, sizeof( p->Rank[0] ) );
    // INCHI✔️✔️:     if (p->AtNumber && p->Rank)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: PartitionCreate

    p.AtNumber = crate::source::base::util::inchi_calloc::<AT_NUMB>(
        heap,
        n as u64,
        size_of::<AT_NUMB>() as u64,
    )
    .unwrap_or_else(|_| SourceMutPointer::null());
    p.Rank = crate::source::base::util::inchi_calloc::<AT_RANK>(
        heap,
        n as u64,
        size_of::<AT_RANK>() as u64,
    )
    .unwrap_or_else(|_| SourceMutPointer::null());
    Ok(if !p.AtNumber.is_null() && !p.Rank.is_null() {
        1
    } else {
        0
    })
}

#[allow(non_snake_case)]
pub(crate) fn PartitionIsDiscrete(
    heap: &SourceHeap,
    p: &Partition,
    n: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:978 PartitionIsDiscrete
    // INCHI✔️✔️: int PartitionIsDiscrete( Partition *p, int n )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i;
    // INCHI✔️✔️:     AT_RANK r;
    // INCHI✔️✔️:     for (i = 0, r = 1; i < n; i++, r++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (r != ( rank_mask_bit & p->Rank[p->AtNumber[i]] ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             INCHI_HEAPCHK
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     INCHI_HEAPCHK
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 1;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: PartitionIsDiscrete

    if n <= 0 {
        return Ok(1);
    }
    // SAFETY: both partition arrays are independent fixed-size allocations
    // owned by CanonGraph and remain immutable for this source loop.
    let atoms = unsafe { heap.stable_slice(p.AtNumber.as_const())? };
    let ranks = unsafe { heap.stable_slice(p.Rank.as_const())? };
    let count = usize::try_from(n).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atoms = atoms.prefix(count)?;
    let ranks = ranks.prefix(ranks.len())?;
    let rank_mask = rank_mask_bit();
    // INCHI✔️✔️:     for (i = 0, r = 1; i < n; i++, r++)
    let mut i = 0_usize;
    let mut rank = 1 as AT_RANK;
    while i < count {
        let atom = usize::from(atoms[i]);
        // INCHI✔️✔️:         if (r != ( rank_mask_bit & p->Rank[p->AtNumber[i]] ))
        let atom_rank = ranks
            .get(atom)
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if rank != rank_mask & atom_rank {
            // INCHI✔️✔️:             return 0;
            return Ok(0);
        }
        i += 1;
        rank = rank.wrapping_add(1);
    }
    // INCHI✔️✔️:     return 1;
    Ok(1)
}

#[allow(non_snake_case)]
pub(crate) fn PartitionGetFirstCell(
    heap: &SourceHeap,
    p: &Partition,
    baseW: &mut [Cell],
    k: i32,
    n: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:998 PartitionGetFirstCell
    // INCHI✔️✔️: int PartitionGetFirstCell( Partition *p, Cell *baseW, int k, int n )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i;
    // INCHI✔️✔️:     AT_RANK r;
    // INCHI✔️✔️:     Cell *W = baseW + k - 1;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     i = ( k > 1 ) ? baseW[k - 2].first + 1 : 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (i < n)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* bypass single vertex cells */
    // INCHI✔️✔️:         for (r = (AT_RANK) ( i + 1 ); i < n && r == ( rank_mask_bit & p->Rank[(int) p->AtNumber[i]] ); i++, r++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             ;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (i < n)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         W->first = i;
    // INCHI✔️✔️:         for (r = ( rank_mask_bit & p->Rank[(int) p->AtNumber[i]] ), i++;
    // INCHI✔️✔️:               i < n && r == ( rank_mask_bit & p->Rank[(int) p->AtNumber[i]] );
    // INCHI✔️✔️:               i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             ;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         W->next = i;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         INCHI_HEAPCHK
    // INCHI✔️✔️:
    // INCHI✔️✔️:         return ( W->next - W->first );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     W->first = INCHI_CANON_INFINITY;
    // INCHI✔️✔️:     W->next = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     INCHI_HEAPCHK
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: PartitionGetFirstCell

    let target =
        usize::try_from(k.wrapping_sub(1)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let _ = baseW
        .get(target)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut i = if k > 1 {
        let previous =
            usize::try_from(k.wrapping_sub(2)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        baseW
            .get(previous)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .first
            .wrapping_add(1)
    } else {
        0
    };
    if i < n {
        let mut rank = i.wrapping_add(1) as AT_RANK;
        while i < n {
            let atom = source_get(heap, p.AtNumber, i)?;
            if rank != rank_mask_bit() & source_get(heap, p.Rank, i32::from(atom))? {
                break;
            }
            i = i.wrapping_add(1);
            rank = rank.wrapping_add(1);
        }
    }
    if i < n {
        let first = i;
        baseW
            .get_mut(target)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .first = first;
        let atom = source_get(heap, p.AtNumber, i)?;
        let rank = rank_mask_bit() & source_get(heap, p.Rank, i32::from(atom))?;
        i = i.wrapping_add(1);
        while i < n {
            let atom = source_get(heap, p.AtNumber, i)?;
            if rank != rank_mask_bit() & source_get(heap, p.Rank, i32::from(atom))? {
                break;
            }
            i = i.wrapping_add(1);
        }
        let cell = baseW
            .get_mut(target)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        cell.next = i;
        return Ok(cell.next.wrapping_sub(first));
    }
    let cell = baseW
        .get_mut(target)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    cell.first = INCHI_CANON_INFINITY as i32;
    cell.next = 0;
    Ok(0)
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

#[allow(non_snake_case)]
pub(crate) fn CtPartCopy(
    heap: &mut SourceHeap,
    Ct1: &mut ConTable,
    Ct2: &ConTable,
    mut k: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2965 CtPartCopy
    // INCHI✔️❌: void CtPartCopy( ConTable *Ct1 /* to */,
    // INCHI✔️❌:                  ConTable *Ct2 /* from */,
    // INCHI✔️❌:                  int k )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int     startCt1, startCt2, endCt2;
    // INCHI✔️❌:     int     len2, len2H, len2Hfixed = 0, len2iso_sort_key, len2iso_exchg_atnos, i;
    // INCHI✔️❌:     int     startAt1, startAt2, endAt2; /*endCt,*/
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:     int     len2iso_sort_key_Hfixed;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     k--;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (k)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         startCt1 = Ct1->nextCtblPos[k - 1];
    // INCHI✔️❌:         startCt2 = Ct2->nextCtblPos[k - 1];
    // INCHI✔️❌:         startAt1 = Ct1->nextAtRank[k - 1] - 1;
    // INCHI✔️❌:         startAt2 = Ct2->nextAtRank[k - 1] - 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         startCt1 = startCt2 = 0;
    // INCHI✔️❌:         startAt1 = startAt2 = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     endCt2 = Ct2->nextCtblPos[k];
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     endAt2 = (int) Ct2->nextAtRank[k] - 1;
    // INCHI✔️❌:
    // INCHI✔️❌:     len2 = endCt2 - startCt2;
    // INCHI✔️❌:     /* len    = min(len1, len2); */
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI✔️❌:     if (startCt1 != startCt2 || startAt1 != startAt2)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int stop = 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     /* copy connection table: Ctbl */
    // INCHI✔️❌:     for (i = 0; i < len2; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Ct1->Ctbl[startCt1 + i] = Ct2->Ctbl[startCt2 + i];
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* copy number of H: NumH */
    // INCHI✔️❌:     len2H = 0;
    // INCHI✔️❌:     if (Ct1->NumH && Ct2->NumH)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         len2H = endAt2 - startAt2;
    // INCHI✔️❌:         if (endAt2 > Ct2->maxVert)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             len2H = Ct2->lenNumH - startAt2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         for (i = 0; i < len2H; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             Ct1->NumH[startAt1 + i] = Ct2->NumH[startAt2 + i];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* copy number of fixed H */
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     if (Ct1->NumHfixed && Ct2->NumHfixed)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         len2Hfixed = endAt2 - startAt2;
    // INCHI✔️❌:         for (i = 0; i < len2Hfixed; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             Ct1->NumHfixed[startAt1 + i] = Ct2->NumHfixed[startAt2 + i];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* copy isotopic keys */
    // INCHI✔️❌:     len2iso_sort_key = 0;
    // INCHI✔️❌:     if (Ct1->iso_sort_key && Ct2->iso_sort_key)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         len2iso_sort_key = endAt2 - startAt2;
    // INCHI✔️❌:         for (i = 0; i < len2iso_sort_key; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             Ct1->iso_sort_key[startAt1 + i] = Ct2->iso_sort_key[startAt2 + i];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     len2iso_exchg_atnos = 0;
    // INCHI✔️❌:     if (Ct1->iso_exchg_atnos && Ct2->iso_exchg_atnos)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         len2iso_exchg_atnos = endAt2 - startAt2;
    // INCHI✔️❌:         for (i = 0; i < len2iso_exchg_atnos; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             Ct1->iso_exchg_atnos[startAt1 + i] = Ct2->iso_exchg_atnos[startAt2 + i];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:     len2iso_sort_key_Hfixed = 0;
    // INCHI✔️❌:     if (Ct1->iso_sort_key_Hfixed && Ct2->iso_sort_key_Hfixed)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         len2iso_sort_key_Hfixed = endAt2 - startAt2;
    // INCHI✔️❌:         for (i = 0; i < len2iso_sort_key; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             Ct1->iso_sort_key_Hfixed[startAt1 + i] = Ct2->iso_sort_key_Hfixed[startAt2 + i];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     Ct1->lenCt = startCt1 + len2;
    // INCHI✔️❌:     Ct1->nextCtblPos[k] = startCt1 + len2;
    // INCHI✔️❌:     Ct1->nextAtRank[k] = Ct2->nextAtRank[k];
    // INCHI✔️❌:     if (len2H)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Ct1->lenNumH = startAt1 + len2H;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( len2Hfixed )
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Ct1->lenNumHfixed   = startAt1 + len2Hfixed;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (len2iso_sort_key)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Ct1->len_iso_sort_key = startAt1 + len2iso_sort_key;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (len2iso_exchg_atnos)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Ct1->len_iso_exchg_atnos = startAt1 + len2iso_exchg_atnos;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:     if (len2iso_sort_key_Hfixed)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Ct1->len_iso_sort_key_Hfixed = startAt1 + len2iso_sort_key_Hfixed;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef INCHI_CANON_USE_HASH
    // INCHI✔️❌:     Ct1->hash[k] = Ct2->hash[k];
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     Ct1->lenPos = k + 1;
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CtPartCopy
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CtPartCopy
    // INCHI✔️❌: #define bRELEASE_VERSION  1    /* 1=> release version; comment out to disable */
    // INCHI✔️❌: #define USE_ISO_SORT_KEY_HFIXED  0  /* 0=> normal mode: merge isotopic taut H to isotopic atom sorting key in
    // INCHI✔️❌:                                            taut H-fixed canonicalization;
    // INCHI✔️❌:                                        1=> add one more "string" iso_sort_Hfixed to the canonicalization */
    // INCHI✔️❌: /*
    // INCHI✔️❌: #define INCHI_CANON_USE_HASH
    // INCHI✔️❌: */ /* djb-rwth: constant has not been defined? */
    // INCHI✔️❌: #define INCHI_HEAPCHK          /* default: no explicit heap checking during the execution */
    // END INCHI ACTIVE MACRO CONFIGURATION: CtPartCopy

    k = k.wrapping_sub(1);
    let (start_ct1, start_ct2, start_at1, start_at2) = if k != 0 {
        (
            i32::from(source_get(heap, Ct1.nextCtblPos, k.wrapping_sub(1))?),
            i32::from(source_get(heap, Ct2.nextCtblPos, k.wrapping_sub(1))?),
            i32::from(source_get(heap, Ct1.nextAtRank, k.wrapping_sub(1))?).wrapping_sub(1),
            i32::from(source_get(heap, Ct2.nextAtRank, k.wrapping_sub(1))?).wrapping_sub(1),
        )
    } else {
        (0, 0, 0, 0)
    };

    let end_ct2 = i32::from(source_get(heap, Ct2.nextCtblPos, k)?);
    let end_at2 = i32::from(source_get(heap, Ct2.nextAtRank, k)?).wrapping_sub(1);
    let len2 = end_ct2.wrapping_sub(start_ct2);

    // SAFETY: CtPartCopy is called with two ConTables created independently by
    // CTableCreate. Corresponding source and destination fields are distinct,
    // fixed-size allocations that remain live for the complete field loop.
    let source_ctbl = unsafe { heap.stable_slice(Ct2.Ctbl.as_const())? };
    let mut target_ctbl = unsafe { heap.stable_slice_mut(Ct1.Ctbl)? };
    // INCHI✔️✔️:     for (i = 0; i < len2; i++)
    let mut i = 0_i32;
    while i < len2 {
        let source_index = usize::try_from(start_ct2.wrapping_add(i))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let target_index = usize::try_from(start_ct1.wrapping_add(i))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        // INCHI✔️✔️:         Ct1->Ctbl[startCt1 + i] = Ct2->Ctbl[startCt2 + i];
        *target_ctbl.get_mut(target_index)? = *source_ctbl.get(source_index)?;
        i = i.wrapping_add(1);
    }
    drop(source_ctbl);
    drop(target_ctbl);

    // INCHI✔️✔️:     len2H = 0;
    let mut len2_h = 0_i32;
    // INCHI✔️✔️:     if (Ct1->NumH && Ct2->NumH)
    if !Ct1.NumH.is_null() && !Ct2.NumH.is_null() {
        // INCHI✔️✔️:         len2H = endAt2 - startAt2;
        len2_h = end_at2.wrapping_sub(start_at2);
        // INCHI✔️✔️:         if (endAt2 > Ct2->maxVert)
        if end_at2 > Ct2.maxVert {
            // INCHI✔️✔️:             len2H = Ct2->lenNumH - startAt2;
            len2_h = Ct2.lenNumH.wrapping_sub(start_at2);
        }
        let source_h = unsafe { heap.stable_slice(Ct2.NumH.as_const())? };
        let mut target_h = unsafe { heap.stable_slice_mut(Ct1.NumH)? };
        // INCHI✔️✔️:         for (i = 0; i < len2H; i++)
        i = 0;
        while i < len2_h {
            let source_index = usize::try_from(start_at2.wrapping_add(i))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let target_index = usize::try_from(start_at1.wrapping_add(i))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            // INCHI✔️✔️:             Ct1->NumH[startAt1 + i] = Ct2->NumH[startAt2 + i];
            *target_h.get_mut(target_index)? = *source_h.get(source_index)?;
            i = i.wrapping_add(1);
        }
    }

    // INCHI✔️✔️:     if (Ct1->NumHfixed && Ct2->NumHfixed)
    if !Ct1.NumHfixed.is_null() && !Ct2.NumHfixed.is_null() {
        // INCHI✔️✔️:         len2Hfixed = endAt2 - startAt2;
        let len2_h_fixed = end_at2.wrapping_sub(start_at2);
        let source_fixed_h = unsafe { heap.stable_slice(Ct2.NumHfixed.as_const())? };
        let mut target_fixed_h = unsafe { heap.stable_slice_mut(Ct1.NumHfixed)? };
        // INCHI✔️✔️:         for (i = 0; i < len2Hfixed; i++)
        i = 0;
        while i < len2_h_fixed {
            let source_index = usize::try_from(start_at2.wrapping_add(i))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let target_index = usize::try_from(start_at1.wrapping_add(i))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            // INCHI✔️✔️:             Ct1->NumHfixed[startAt1 + i] = Ct2->NumHfixed[startAt2 + i];
            *target_fixed_h.get_mut(target_index)? = *source_fixed_h.get(source_index)?;
            i = i.wrapping_add(1);
        }
    }

    // INCHI✔️✔️:     len2iso_sort_key = 0;
    let mut len2_iso_sort_key = 0_i32;
    // INCHI✔️✔️:     if (Ct1->iso_sort_key && Ct2->iso_sort_key)
    if !Ct1.iso_sort_key.is_null() && !Ct2.iso_sort_key.is_null() {
        // INCHI✔️✔️:         len2iso_sort_key = endAt2 - startAt2;
        len2_iso_sort_key = end_at2.wrapping_sub(start_at2);
        let source_iso = unsafe { heap.stable_slice(Ct2.iso_sort_key.as_const())? };
        let mut target_iso = unsafe { heap.stable_slice_mut(Ct1.iso_sort_key)? };
        // INCHI✔️✔️:         for (i = 0; i < len2iso_sort_key; i++)
        i = 0;
        while i < len2_iso_sort_key {
            let source_index = usize::try_from(start_at2.wrapping_add(i))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let target_index = usize::try_from(start_at1.wrapping_add(i))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            // INCHI✔️✔️:             Ct1->iso_sort_key[startAt1 + i] = Ct2->iso_sort_key[startAt2 + i];
            *target_iso.get_mut(target_index)? = *source_iso.get(source_index)?;
            i = i.wrapping_add(1);
        }
    }

    // INCHI✔️✔️:     len2iso_exchg_atnos = 0;
    let mut len2_iso_exchg_atnos = 0_i32;
    // INCHI✔️✔️:     if (Ct1->iso_exchg_atnos && Ct2->iso_exchg_atnos)
    if !Ct1.iso_exchg_atnos.is_null() && !Ct2.iso_exchg_atnos.is_null() {
        // INCHI✔️✔️:         len2iso_exchg_atnos = endAt2 - startAt2;
        len2_iso_exchg_atnos = end_at2.wrapping_sub(start_at2);
        let source_exchange = unsafe { heap.stable_slice(Ct2.iso_exchg_atnos.as_const())? };
        let mut target_exchange = unsafe { heap.stable_slice_mut(Ct1.iso_exchg_atnos)? };
        // INCHI✔️✔️:         for (i = 0; i < len2iso_exchg_atnos; i++)
        i = 0;
        while i < len2_iso_exchg_atnos {
            let source_index = usize::try_from(start_at2.wrapping_add(i))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let target_index = usize::try_from(start_at1.wrapping_add(i))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            // INCHI✔️✔️:             Ct1->iso_exchg_atnos[startAt1 + i] = Ct2->iso_exchg_atnos[startAt2 + i];
            *target_exchange.get_mut(target_index)? = *source_exchange.get(source_index)?;
            i = i.wrapping_add(1);
        }
    }

    // INCHI✔️✔️:     Ct1->lenCt = startCt1 + len2;
    Ct1.lenCt = start_ct1.wrapping_add(len2);
    // INCHI✔️✔️:     Ct1->nextCtblPos[k] = startCt1 + len2;
    source_set(heap, Ct1.nextCtblPos, k, Ct1.lenCt as AT_NUMB)?;
    let next_at_rank = source_get(heap, Ct2.nextAtRank, k)?;
    // INCHI✔️✔️:     Ct1->nextAtRank[k] = Ct2->nextAtRank[k];
    source_set(heap, Ct1.nextAtRank, k, next_at_rank)?;
    // INCHI✔️✔️:     if (len2H)
    if len2_h != 0 {
        // INCHI✔️✔️:         Ct1->lenNumH = startAt1 + len2H;
        Ct1.lenNumH = start_at1.wrapping_add(len2_h);
    }
    // INCHI✔️✔️:     if (len2iso_sort_key)
    if len2_iso_sort_key != 0 {
        // INCHI✔️✔️:         Ct1->len_iso_sort_key = startAt1 + len2iso_sort_key;
        Ct1.len_iso_sort_key = start_at1.wrapping_add(len2_iso_sort_key);
    }
    // INCHI✔️✔️:     if (len2iso_exchg_atnos)
    if len2_iso_exchg_atnos != 0 {
        // INCHI✔️✔️:         Ct1->len_iso_exchg_atnos = startAt1 + len2iso_exchg_atnos;
        Ct1.len_iso_exchg_atnos = start_at1.wrapping_add(len2_iso_exchg_atnos);
    }
    // INCHI✔️✔️:     Ct1->lenPos = k + 1;
    Ct1.lenPos = k.wrapping_add(1);
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn CtFullCopy(
    heap: &mut SourceHeap,
    Ct1: &mut ConTable,
    Ct2: &ConTable,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3113 CtFullCopy
    // INCHI✔️❌: void CtFullCopy( ConTable *Ct1, ConTable *Ct2 )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /* Ct1 does not have INCHI_CANON_INFINITY termination */
    // INCHI✔️❌:     int k;
    // INCHI✔️❌:     for (k = 0; k < Ct2->lenPos; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         CtPartCopy( Ct1 /* to */, Ct2 /* from */, k + 1 );
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CtFullCopy

    let mut k = 0_i32;
    while k < Ct2.lenPos {
        CtPartCopy(heap, Ct1, Ct2, k.wrapping_add(1))?;
        k = k.wrapping_add(1);
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn CellMakeEmpty(baseW: &mut [Cell], mut k: i32) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1040 CellMakeEmpty
    // INCHI✔️✔️: void CellMakeEmpty( Cell *baseW, int k )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     k--;
    // INCHI✔️✔️:     baseW[k].first = INCHI_CANON_INFINITY;
    // INCHI✔️✔️:     baseW[k].next = 0;
    // INCHI✔️✔️:     baseW[k].prev = -1;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     INCHI_HEAPCHK
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: CellMakeEmpty

    k = k.wrapping_sub(1);
    let cell = baseW
        .get_mut(usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    cell.first = INCHI_CANON_INFINITY as i32;
    cell.next = 0;
    cell.prev = -1;
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn PartitionFree(
    heap: &mut SourceHeap,
    p: Option<&mut Partition>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:959 PartitionFree
    // INCHI✔️❌: void PartitionFree( Partition *p )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (p)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (p->AtNumber)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( p->AtNumber );
    // INCHI✔️❌:             p->AtNumber = NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (p->Rank)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( p->Rank );
    // INCHI✔️❌:             p->Rank = NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: PartitionFree

    let Some(p) = p else {
        return Ok(());
    };
    if !p.AtNumber.is_null() {
        inchi_free(heap, p.AtNumber)?;
        p.AtNumber = SourceMutPointer::null();
    }
    if !p.Rank.is_null() {
        inchi_free(heap, p.Rank)?;
        p.Rank = SourceMutPointer::null();
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn SetInitialRanks2(
    heap: &mut SourceHeap,
    num_atoms: i32,
    pAtomInvariant2: SourceMutPointer<ATOM_INVARIANT2>,
    nNewRank: SourceMutPointer<AT_RANK>,
    nAtomNumber: SourceMutPointer<AT_RANK>,
    pCG: &mut CANON_GLOBALS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4674 SetInitialRanks2
    // INCHI✔️✔️: int SetInitialRanks2( int num_atoms,
    // INCHI✔️✔️:                       ATOM_INVARIANT2* pAtomInvariant2,
    // INCHI✔️✔️:                       AT_RANK *nNewRank,
    // INCHI✔️✔️:                       AT_RANK *nAtomNumber,
    // INCHI✔️✔️:                       CANON_GLOBALS *pCG )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, nNumDiffRanks;
    // INCHI✔️✔️:     AT_RANK nCurrentRank;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nAtomNumber[i] = (AT_RANK) i;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* global for qsort */
    // INCHI✔️✔️:     pCG->m_pAtomInvariant2ForSort = pAtomInvariant2;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     inchi_qsort( pCG, nAtomNumber, num_atoms, sizeof( nAtomNumber[0] ), CompAtomInvariants2 );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* nNewRank[i]: non-decreading order; do not increment nCurrentRank */
    // INCHI✔️✔️:     /*           if consecutive sorted atom invariants are identical */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* djb-rwth: fixing oss-fuzz issue #69315 */
    // INCHI✔️✔️:     nNumDiffRanks = 1;
    // INCHI✔️✔️:     if (num_atoms > 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nCurrentRank = (AT_RANK)num_atoms;
    // INCHI✔️✔️:         nNewRank[nAtomNumber[num_atoms - 1]] = nCurrentRank;
    // INCHI✔️✔️:         for (i = num_atoms - 1; i > 0; i--)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* Note: CompAtomInvariants2Only() in following line implicitly reads pAtomInvariant2 pointed by pAtomInvariant2ForSort */
    // INCHI✔️✔️:             if (CompAtomInvariants2Only(&nAtomNumber[i - 1], &nAtomNumber[i], pCG))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 nNumDiffRanks++;
    // INCHI✔️✔️:                 nCurrentRank = (AT_RANK)i;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             nNewRank[nAtomNumber[i - 1]] = nCurrentRank;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return nNumDiffRanks;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: SetInitialRanks2

    let count = usize::try_from(num_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if count > 0 {
        let atom_numbers = heap.slice_mut(nAtomNumber)?;
        let atom_numbers = atom_numbers
            .get_mut(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for (index, atom_number) in atom_numbers.iter_mut().enumerate() {
            *atom_number = index as AT_RANK;
        }
    }

    pCG.m_pAtomInvariant2ForSort = pAtomInvariant2.as_const();

    let invariant_storage = if count > 1 {
        // SAFETY: ATOM_INVARIANT2 storage cannot alias either AT_RANK
        // allocation by SourceHeap's allocation type invariant. This function
        // neither resizes nor frees it while sorting atom numbers and writing
        // ranks, so the source array remains stable for both comparator loops.
        Some(unsafe { heap.stable_slice(pAtomInvariant2.as_const())? })
    } else {
        None
    };
    let invariant_workspace = invariant_storage
        .as_ref()
        .map(|invariants| AtomInvariant2SortWorkspace::new(invariants.prefix(count)?, count))
        .transpose()?;

    if count > 1 {
        let invariant_workspace = invariant_workspace
            .as_ref()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        heap.with_slice_mut_and_heap(nAtomNumber, |atom_numbers, _heap| {
            let atom_numbers = atom_numbers
                .get_mut(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let byte_len = count
                .checked_mul(size_of::<AT_RANK>())
                .ok_or(SourceHeapError::AllocationSizeOverflow)?;
            // AT_RANK has no padding; this exposes the same contiguous records passed to C qsort.
            let bytes = unsafe {
                std::slice::from_raw_parts_mut(atom_numbers.as_mut_ptr().cast::<u8>(), byte_len)
            };
            inchi_qsort(bytes, count, size_of::<AT_RANK>(), &mut |first, second| {
                let first = AT_RANK::from_ne_bytes(
                    first
                        .try_into()
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                );
                let second = AT_RANK::from_ne_bytes(
                    second
                        .try_into()
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                );
                // SAFETY: the source initialized every record to 0..count and
                // inchi_qsort only permutes those records.
                Ok(unsafe { invariant_workspace.compare_in_bounds(first, second) })
            })
        })?;
    }
    if count > 0 {
        // The source initializes 0..num_atoms and qsort only permutes it.
        heap.record_index_bound(nAtomNumber, count, count)?;
    }

    let mut nNumDiffRanks = 1_i32;
    if count > 0 {
        heap.with_slice_mut_and_heap(nNewRank, |new_ranks, heap| {
            let atom_numbers = heap.slice(nAtomNumber.as_const())?;
            let atom_numbers = atom_numbers
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mut nCurrentRank = num_atoms as AT_RANK;
            let last_atom = usize::from(atom_numbers[count - 1]);
            *new_ranks
                .get_mut(last_atom)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = nCurrentRank;
            let mut i = count - 1;
            while i > 0 {
                // SAFETY: the source initialized 0..count and qsort only
                // permuted that complete array above.
                if unsafe {
                    invariant_workspace
                        .as_ref()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .compare_only_in_bounds(atom_numbers[i - 1], atom_numbers[i])
                } != 0
                {
                    nNumDiffRanks = nNumDiffRanks.wrapping_add(1);
                    nCurrentRank = i as AT_RANK;
                }
                let atom = usize::from(atom_numbers[i - 1]);
                *new_ranks
                    .get_mut(atom)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = nCurrentRank;
                i -= 1;
            }
            Ok(())
        })?;
    }

    Ok(nNumDiffRanks)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn FillOutAtomInvariant2(
    heap: &mut SourceHeap,
    at: SourceConstPointer<sp_ATOM>,
    num_atoms: i32,
    num_at_tg: i32,
    pAtomInvariant: SourceMutPointer<ATOM_INVARIANT2>,
    bIgnoreIsotopic: i32,
    bHydrogensInRanks: i32,
    bHydrogensFixedInRanks: i32,
    bDigraph: i32,
    bTautGroupsOnly: i32,
    t_group_info: Option<&T_GROUP_INFO>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4719 FillOutAtomInvariant2
    // INCHI✔️✔️: void FillOutAtomInvariant2( sp_ATOM* at,
    // INCHI✔️✔️:                             int num_atoms,
    // INCHI✔️✔️:                             int num_at_tg,
    // INCHI✔️✔️:                             ATOM_INVARIANT2* pAtomInvariant,
    // INCHI✔️✔️:                             int bIgnoreIsotopic,
    // INCHI✔️✔️:                             int bHydrogensInRanks,
    // INCHI✔️✔️:                             int bHydrogensFixedInRanks,
    // INCHI✔️✔️:                             int bDigraph,
    // INCHI✔️✔️:                             int bTautGroupsOnly,
    // INCHI✔️✔️:                             T_GROUP_INFO *t_group_info )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, k, j, i_t_group;
    // INCHI✔️✔️:     /* tautomers */
    // INCHI✔️✔️:     T_GROUP          *t_group = NULL;
    // INCHI✔️✔️:     int               num_t_groups = 0;
    // INCHI✔️✔️:     int               num_tautomer_iso = 0;
    // INCHI✔️✔️: #define ELEM_NAME_LEN  2
    // INCHI✔️✔️:     char ChemElements[ELEM_NAME_LEN*NUM_CHEM_ELEMENTS + ELEM_NAME_LEN];
    // INCHI✔️✔️:     char CurElement[ELEM_NAME_LEN + ELEM_NAME_LEN], *pCurElem;
    // INCHI✔️✔️:     int  nNumChemElements = 0;
    // INCHI✔️✔️:     int  nNumHydrogenAtoms = 0;
    // INCHI✔️✔️:     int  nNumCarbonAtoms = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     ChemElements[0] = '\0'; /* djb-rwth: initialisation prevents empty string comparison */
    // INCHI✔️✔️:     memset( ChemElements, 0, sizeof( ChemElements ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:     memset( CurElement, 0, sizeof( CurElement ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:     nNumChemElements = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (num_at_tg > num_atoms && t_group_info)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         t_group = t_group_info->t_group;
    // INCHI✔️✔️:         num_t_groups = t_group_info->num_t_groups;
    // INCHI✔️✔️:         num_tautomer_iso = t_group_info->bIgnoreIsotopic ? 0 : T_NUM_ISOTOPIC;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (!bTautGroupsOnly)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:
    // INCHI✔️✔️:         for (i = 0; i < num_atoms; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (!strcmp( at[i].elname, "C" ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 nNumCarbonAtoms++;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else if (!strcmp( at[i].elname, "H" ) ||
    // INCHI✔️✔️:                  !strcmp( at[i].elname, "D" ) ||
    // INCHI✔️✔️:                  !strcmp( at[i].elname, "T" ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 nNumHydrogenAtoms++;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 CurElement[0] = at[i].elname[0];
    // INCHI✔️✔️:                 CurElement[1] = at[i].elname[1] ? at[i].elname[1] : ' ';
    // INCHI✔️✔️:                 if (!( pCurElem = strstr( ChemElements, CurElement ) )) /* djb-rwth: ignoring LLVM warning: variable used for function return value */
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     strcat(ChemElements, CurElement);
    // INCHI✔️✔️:                     nNumChemElements++;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (nNumChemElements > 1)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             qsort( ChemElements, nNumChemElements, ELEM_NAME_LEN, CompChemElemLex );
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (nNumCarbonAtoms)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (nNumChemElements)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 memmove(ChemElements + ELEM_NAME_LEN, ChemElements, (long long)ELEM_NAME_LEN * (long long)nNumChemElements); /* djb-rwth: cast operators added */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             ChemElements[0] = 'C';
    // INCHI✔️✔️:             ChemElements[1] = ' ';
    // INCHI✔️✔️:             nNumChemElements++;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (nNumHydrogenAtoms)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             ChemElements[ELEM_NAME_LEN*nNumChemElements] = 'H';
    // INCHI✔️✔️:             ChemElements[ELEM_NAME_LEN*nNumChemElements + 1] = ' ';
    // INCHI✔️✔️:             nNumChemElements++;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /* general */
    // INCHI✔️✔️:         for (i = 0; i < num_atoms; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             memset( &pAtomInvariant[i], 0, sizeof( pAtomInvariant[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:             CurElement[0] = at[i].elname[0];
    // INCHI✔️✔️:             CurElement[1] = at[i].elname[1] ? at[i].elname[1] : ' ';
    // INCHI✔️✔️:             pCurElem = strstr( ChemElements, CurElement );
    // INCHI✔️✔️:             if (pCurElem)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 j = (int) ( pCurElem - ChemElements ) / ELEM_NAME_LEN + 1;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 j = nNumChemElements; /* must be D or T */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             /* at[i].hill_type = (U_CHAR) j; */
    // INCHI✔️✔️:             pAtomInvariant[i].val[AT_INV_HILL_ORDER] = j;
    // INCHI✔️✔️:
    // INCHI✔️✔️:             pAtomInvariant[i].val[AT_INV_NUM_CONNECTIONS] = at[i].valence;
    // INCHI✔️✔️:             if (bHydrogensInRanks)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 pAtomInvariant[i].val[AT_INV_NUM_H] = ( ( t_group && at[i].endpoint > 0 ) ? 0 : at[i].num_H );
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (bHydrogensFixedInRanks)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 pAtomInvariant[i].val[AT_INV_NUM_H_FIX] = ( ( t_group && at[i].endpoint > 0 ) ? at[i].num_H : 0 );
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (!bDigraph &&  t_group && ( i_t_group = (int) at[i].endpoint - 1 ) >= 0 && i_t_group < num_t_groups)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 pAtomInvariant[i].val[AT_INV_NUM_TG_ENDPOINTS] = t_group[i_t_group].nNumEndpoints;
    // INCHI✔️✔️:                 for (j = 0; j < T_NUM_NO_ISOTOPIC; j++)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     pAtomInvariant[i].val[AT_INV_TG_NUMBERS + j] = t_group[i_t_group].num[j];
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 for (j = 0; j < num_tautomer_iso; j++)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     pAtomInvariant[i].val[AT_INV_TAUT_ISO + j] = t_group[i_t_group].num[j + T_NUM_NO_ISOTOPIC];
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             pAtomInvariant[i].iso_sort_key = bIgnoreIsotopic ? 0 : at[i].iso_sort_key;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* fill tautomeric groups only */
    // INCHI✔️✔️:         memset( pAtomInvariant, 0, num_at_tg * sizeof( pAtomInvariant[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /**************************************/
    // INCHI✔️✔️:     /*          tautomeric groups         */
    // INCHI✔️✔️:     /**************************************/
    // INCHI✔️✔️:     for (i = num_atoms; i < num_at_tg; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:
    // INCHI✔️✔️:         k = i - num_atoms;
    // INCHI✔️✔️:         memset( &pAtomInvariant[i], 0, sizeof( pAtomInvariant[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:         if (!t_group)
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         /* make sure ranks of t-groups are larger than that of any atom */
    // INCHI✔️✔️:          /* greater than for any real atom */
    // INCHI✔️✔️:         pAtomInvariant[i].val[AT_INV_HILL_ORDER] = bTautGroupsOnly ? num_at_tg : nNumChemElements + 1;
    // INCHI✔️✔️:         /* greater than for any real atom */
    // INCHI✔️✔️:         pAtomInvariant[i].val[AT_INV_NUM_CONNECTIONS] = MAXVAL + 1;
    // INCHI✔️✔️:         if (k < num_t_groups)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             pAtomInvariant[i].val[AT_INV_NUM_TG_ENDPOINTS] = t_group[k].nNumEndpoints;
    // INCHI✔️✔️:             for (j = 0; j < T_NUM_NO_ISOTOPIC; j++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 pAtomInvariant[i].val[AT_INV_TAUT_ISO + j] = t_group[k].num[j];
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             for (j = 0; j < num_tautomer_iso; j++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 pAtomInvariant[i].val[AT_INV_TAUT_ISO + j] = t_group[k].num[j + T_NUM_NO_ISOTOPIC];
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: FillOutAtomInvariant2
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FillOutAtomInvariant2
    // INCHI✔️✔️: #define ELEM_NAME_LEN 2
    // INCHI✔️✔️: #define NUM_CHEM_ELEMENTS 127
    // INCHI✔️✔️: #define MAXVAL 20
    // INCHI✔️✔️: #define T_NUM_NO_ISOTOPIC 2
    // INCHI✔️✔️: #define T_NUM_ISOTOPIC 3
    // END INCHI ACTIVE MACRO CONFIGURATION: FillOutAtomInvariant2

    const ELEM_NAME_LEN: usize = 2;
    let atom_count = usize::try_from(num_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let total_count =
        usize::try_from(num_at_tg).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if total_count < atom_count {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    heap.with_slice_mut_and_heap(pAtomInvariant, |invariants, heap| {
        let invariants = invariants
            .get_mut(..total_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let atoms = if atom_count == 0 {
            &[][..]
        } else {
            heap.slice(at)?
                .get(..atom_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        };
        let (t_group_pointer, num_t_groups, num_tautomer_iso) = if num_at_tg > num_atoms
            && let Some(info) = t_group_info
        {
            (
                Some(info.t_group),
                info.num_t_groups,
                if info.bIgnoreIsotopic != 0 {
                    0
                } else {
                    T_NUM_ISOTOPIC as i32
                },
            )
        } else {
            (None, 0, 0)
        };
        let t_groups: Option<&[T_GROUP]> = match t_group_pointer {
            Some(pointer) if !pointer.is_null() => Some(heap.slice(pointer.as_const())?),
            _ => None,
        };

        let mut chem_elements = [0_u8; NUM_CHEM_ELEMENTS as usize * ELEM_NAME_LEN + ELEM_NAME_LEN];
        let mut n_num_chem_elements = 0_usize;
        let mut n_num_hydrogen_atoms = 0_i32;
        let mut n_num_carbon_atoms = 0_i32;
        let find_element = |elements: &[u8], pair: [u8; 2]| -> Option<usize> {
            let end = elements
                .iter()
                .position(|byte| *byte == 0)
                .unwrap_or(elements.len());
            elements[..end].windows(2).position(|window| window == pair)
        };

        if bTautGroupsOnly == 0 {
            for atom in atoms {
                let first = atom.elname[0] as u8;
                let second = atom.elname[1] as u8;
                if first == b'C' && second == 0 {
                    n_num_carbon_atoms = n_num_carbon_atoms.wrapping_add(1);
                } else if matches!(first, b'H' | b'D' | b'T') && second == 0 {
                    n_num_hydrogen_atoms = n_num_hydrogen_atoms.wrapping_add(1);
                } else {
                    let pair = [first, if second != 0 { second } else { b' ' }];
                    if find_element(&chem_elements, pair).is_none() {
                        let start = n_num_chem_elements
                            .checked_mul(ELEM_NAME_LEN)
                            .ok_or(SourceHeapError::AllocationSizeOverflow)?;
                        chem_elements
                            .get_mut(start..start + ELEM_NAME_LEN)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .copy_from_slice(&pair);
                        n_num_chem_elements += 1;
                    }
                }
            }
            if n_num_chem_elements > 1 {
                let byte_count = n_num_chem_elements * ELEM_NAME_LEN;
                inchi_qsort(
                    &mut chem_elements[..byte_count],
                    n_num_chem_elements,
                    ELEM_NAME_LEN,
                    &mut CompChemElemLex,
                )?;
            }
            if n_num_carbon_atoms != 0 {
                if n_num_chem_elements != 0 {
                    let byte_count = n_num_chem_elements * ELEM_NAME_LEN;
                    chem_elements.copy_within(0..byte_count, ELEM_NAME_LEN);
                }
                chem_elements[..2].copy_from_slice(b"C ");
                n_num_chem_elements += 1;
            }
            if n_num_hydrogen_atoms != 0 {
                let start = ELEM_NAME_LEN * n_num_chem_elements;
                chem_elements[start..start + 2].copy_from_slice(b"H ");
                n_num_chem_elements += 1;
            }

            for (index, atom) in atoms.iter().enumerate() {
                invariants[index] = ATOM_INVARIANT2::default();
                let pair = [
                    atom.elname[0] as u8,
                    if atom.elname[1] != 0 {
                        atom.elname[1] as u8
                    } else {
                        b' '
                    },
                ];
                let hill_order = find_element(&chem_elements, pair)
                    .map(|offset| offset / ELEM_NAME_LEN + 1)
                    .unwrap_or(n_num_chem_elements);
                invariants[index].val[tagAtInvariantIndexes_AT_INV_HILL_ORDER as usize] =
                    hill_order as AT_NUMB;
                invariants[index].val[tagAtInvariantIndexes_AT_INV_NUM_CONNECTIONS as usize] =
                    atom.valence as AT_NUMB;
                if bHydrogensInRanks != 0 {
                    invariants[index].val[tagAtInvariantIndexes_AT_INV_NUM_H as usize] =
                        if t_groups.is_some() && atom.endpoint > 0 {
                            0
                        } else {
                            atom.num_H as AT_NUMB
                        };
                }
                if bHydrogensFixedInRanks != 0 {
                    invariants[index].val[tagAtInvariantIndexes_AT_INV_NUM_H_FIX as usize] =
                        if t_groups.is_some() && atom.endpoint > 0 {
                            atom.num_H as AT_NUMB
                        } else {
                            0
                        };
                }
                let group_index = i32::from(atom.endpoint).wrapping_sub(1);
                if bDigraph == 0
                    && let Some(groups) = t_groups
                    && group_index >= 0
                    && group_index < num_t_groups
                {
                    let group = groups
                        .get(group_index as usize)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    invariants[index].val[tagAtInvariantIndexes_AT_INV_NUM_TG_ENDPOINTS as usize] =
                        group.nNumEndpoints;
                    for j in 0..T_NUM_NO_ISOTOPIC as usize {
                        invariants[index].val
                            [tagAtInvariantIndexes_AT_INV_TG_NUMBERS as usize + j] = group.num[j];
                    }
                    for j in 0..num_tautomer_iso as usize {
                        invariants[index].val[tagAtInvariantIndexes_AT_INV_TAUT_ISO as usize + j] =
                            group.num[j + T_NUM_NO_ISOTOPIC as usize];
                    }
                }
                invariants[index].iso_sort_key = if bIgnoreIsotopic != 0 {
                    0
                } else {
                    atom.iso_sort_key
                };
            }
        } else {
            invariants.fill(ATOM_INVARIANT2::default());
        }

        for index in atom_count..total_count {
            let group_index = index - atom_count;
            invariants[index] = ATOM_INVARIANT2::default();
            let Some(groups) = t_groups else {
                continue;
            };
            invariants[index].val[tagAtInvariantIndexes_AT_INV_HILL_ORDER as usize] =
                if bTautGroupsOnly != 0 {
                    num_at_tg as AT_NUMB
                } else {
                    (n_num_chem_elements + 1) as AT_NUMB
                };
            invariants[index].val[tagAtInvariantIndexes_AT_INV_NUM_CONNECTIONS as usize] =
                (MAXVAL + 1) as AT_NUMB;
            if (group_index as i32) < num_t_groups {
                let group = groups
                    .get(group_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                invariants[index].val[tagAtInvariantIndexes_AT_INV_NUM_TG_ENDPOINTS as usize] =
                    group.nNumEndpoints;
                for j in 0..T_NUM_NO_ISOTOPIC as usize {
                    invariants[index].val[tagAtInvariantIndexes_AT_INV_TAUT_ISO as usize + j] =
                        group.num[j];
                }
                for j in 0..num_tautomer_iso as usize {
                    invariants[index].val[tagAtInvariantIndexes_AT_INV_TAUT_ISO as usize + j] =
                        group.num[j + T_NUM_NO_ISOTOPIC as usize];
                }
            }
        }
        Ok(())
    })
}

#[allow(non_snake_case)]
pub(crate) fn CleanNumH(
    heap: &mut SourceHeap,
    NumH: SourceMutPointer<NUM_H>,
    len: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4881 CleanNumH
    // INCHI✔️✔️: void CleanNumH( NUM_H *NumH, int len )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i;
    // INCHI✔️✔️:     if (NumH)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         for (i = 0; i < len; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (NumH[i] == EMPTY_H_NUMBER)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 NumH[i] = 0;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 NumH[i] -= BASE_H_NUMBER;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: CleanNumH
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CleanNumH
    // INCHI✔️✔️: #define EMPTY_H_NUMBER (INCHI_CANON_INFINITY-1) /* 32766 */
    // INCHI✔️✔️: #define BASE_H_NUMBER ((INCHI_CANON_INFINITY-1)/2) /* 16383 */
    // INCHI✔️✔️: typedef signed short NUM_H;
    // END INCHI ACTIVE MACRO CONFIGURATION: CleanNumH

    if NumH.is_null() || len <= 0 {
        return Ok(());
    }
    let count = usize::try_from(len).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let values = heap.slice_mut(NumH)?;
    let values = values
        .get_mut(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for value in values {
        *value = if *value == EMPTY_H_NUMBER as NUM_H {
            0
        } else {
            value.wrapping_sub(BASE_H_NUMBER as NUM_H)
        };
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn CleanCt(
    heap: &mut SourceHeap,
    Ct: SourceMutPointer<AT_RANK>,
    len: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4902 CleanCt
    // INCHI✔️✔️: int CleanCt( AT_RANK *Ct, int len )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     if (Ct && Ct[len] == EMPTY_CT)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         Ct[len] = 0;
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: CleanCt
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CleanCt
    // INCHI✔️✔️: #define EMPTY_CT 0
    // END INCHI ACTIVE MACRO CONFIGURATION: CleanCt

    if Ct.is_null() {
        return Ok(0);
    }
    let index = usize::try_from(len).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let value = heap
        .slice_mut(Ct)?
        .get_mut(index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if *value == EMPTY_CT as AT_RANK {
        *value = 0;
        return Ok(1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CleanIsoSortKeys(
    heap: &mut SourceHeap,
    isk: SourceMutPointer<AT_ISO_SORT_KEY>,
    len: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4915 CleanIsoSortKeys
    // INCHI✔️✔️: void CleanIsoSortKeys( AT_ISO_SORT_KEY * isk, int len )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i;
    // INCHI✔️✔️:     if (isk)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         for (i = 0; i < len; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (isk[i] == EMPTY_ISO_SORT_KEY)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 isk[i] = 0;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: CleanIsoSortKeys
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CleanIsoSortKeys
    // INCHI✔️✔️: #define EMPTY_ISO_SORT_KEY LONG_MAX
    // INCHI✔️✔️: GCC/Linux LP64: LONG_MAX == 9223372036854775807
    // END INCHI ACTIVE MACRO CONFIGURATION: CleanIsoSortKeys

    if isk.is_null() || len <= 0 {
        return Ok(());
    }
    let count = usize::try_from(len).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let values = heap.slice_mut(isk)?;
    let values = values
        .get_mut(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for value in values {
        if *value == AT_ISO_SORT_KEY::MAX {
            *value = 0;
        }
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn DeAllocBCN(
    heap: &mut SourceHeap,
    pBCN: Option<&mut BCN>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4960 DeAllocBCN
    // INCHI✔️❌: void DeAllocBCN( BCN *pBCN )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int    i, k;
    // INCHI✔️❌:     FTCN  *ftcn;
    // INCHI✔️❌:     if (!pBCN)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (pBCN->pRankStack)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = 0; i < pBCN->nMaxLenRankStack; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             FREE_ARRAY( pBCN->pRankStack[i] )
    // INCHI✔️❌:         }
    // INCHI✔️❌:         FREE_ARRAY( pBCN->pRankStack )
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (k = 0; k < TAUT_NUM; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ftcn = pBCN->ftcn + k;
    // INCHI✔️❌:         FreeNeighList( ftcn->NeighList );
    // INCHI✔️❌:
    // INCHI✔️❌:         FREE_ARRAY( ftcn->LinearCt )
    // INCHI✔️❌:
    // INCHI✔️❌:         PartitionFree( &ftcn->PartitionCt );
    // INCHI✔️❌:
    // INCHI✔️❌:         FREE_ARRAY( ftcn->nSymmRankCt )
    // INCHI✔️❌:         FREE_ARRAY( ftcn->nNumHOrig )
    // INCHI✔️❌:         FREE_ARRAY( ftcn->nNumH )
    // INCHI✔️❌:         FREE_ARRAY( ftcn->nNumHOrigFixH )
    // INCHI✔️❌:         FREE_ARRAY( ftcn->nNumHFixH )
    // INCHI✔️❌:
    // INCHI✔️❌:         PartitionFree( &ftcn->PartitionCtIso );
    // INCHI✔️❌:
    // INCHI✔️❌:         FREE_ARRAY( ftcn->nSymmRankCtIso )
    // INCHI✔️❌:         FREE_ARRAY( ftcn->iso_sort_keys )
    // INCHI✔️❌:         FREE_ARRAY( ftcn->iso_sort_keysOrig )
    // INCHI✔️❌:         FREE_ARRAY( ftcn->iso_exchg_atnos )
    // INCHI✔️❌:         FREE_ARRAY( ftcn->iso_exchg_atnosOrig )
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: DeAllocBCN
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: DeAllocBCN
    // INCHI✔️❌: #define FREE_ARRAY( X) if (X) inchi_free( X);
    // INCHI✔️❌: #define TAUT_NUM  2
    // END INCHI ACTIVE MACRO CONFIGURATION: DeAllocBCN

    let Some(pBCN) = pBCN else {
        return Ok(());
    };
    if !pBCN.pRankStack.is_null() {
        let mut i = 0_i32;
        while i < pBCN.nMaxLenRankStack {
            let rank = pointer_array_get(heap, pBCN.pRankStack, i)?;
            if !rank.is_null() {
                inchi_free(heap, rank)?;
            }
            i = i.wrapping_add(1);
        }
        inchi_free(heap, pBCN.pRankStack)?;
    }
    let mut k = 0_usize;
    while k < TAUT_NUM as usize {
        let ftcn = &mut pBCN.ftcn[k];
        FreeNeighList(heap, ftcn.NeighList)?;

        if !ftcn.LinearCt.is_null() {
            inchi_free(heap, ftcn.LinearCt)?;
        }

        PartitionFree(heap, Some(&mut ftcn.PartitionCt))?;

        if !ftcn.nSymmRankCt.is_null() {
            inchi_free(heap, ftcn.nSymmRankCt)?;
        }
        if !ftcn.nNumHOrig.is_null() {
            inchi_free(heap, ftcn.nNumHOrig)?;
        }
        if !ftcn.nNumH.is_null() {
            inchi_free(heap, ftcn.nNumH)?;
        }
        if !ftcn.nNumHOrigFixH.is_null() {
            inchi_free(heap, ftcn.nNumHOrigFixH)?;
        }
        if !ftcn.nNumHFixH.is_null() {
            inchi_free(heap, ftcn.nNumHFixH)?;
        }

        PartitionFree(heap, Some(&mut ftcn.PartitionCtIso))?;

        if !ftcn.nSymmRankCtIso.is_null() {
            inchi_free(heap, ftcn.nSymmRankCtIso)?;
        }
        if !ftcn.iso_sort_keys.is_null() {
            inchi_free(heap, ftcn.iso_sort_keys)?;
        }
        if !ftcn.iso_sort_keysOrig.is_null() {
            inchi_free(heap, ftcn.iso_sort_keysOrig)?;
        }
        if !ftcn.iso_exchg_atnos.is_null() {
            inchi_free(heap, ftcn.iso_exchg_atnos)?;
        }
        if !ftcn.iso_exchg_atnosOrig.is_null() {
            inchi_free(heap, ftcn.iso_exchg_atnosOrig)?;
        }
        k = k.wrapping_add(1);
    }
    Ok(())
}

fn source_prefix_equal<T: PartialEq + 'static>(
    heap: &SourceHeap,
    left: SourceMutPointer<T>,
    right: SourceMutPointer<T>,
    len: i32,
) -> Result<bool, SourceHeapError> {
    let len = usize::try_from(len).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let left = heap
        .slice(left.as_const())?
        .get(..len)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let right = heap
        .slice(right.as_const())?
        .get(..len)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    Ok(left == right)
}

fn con_table_is_usable(table: &ConTable, data: &CANON_DATA) -> bool {
    !table.Ctbl.is_null()
        && !table.nextCtblPos.is_null()
        && !table.nextAtRank.is_null()
        && (data.NumH.is_null() || !table.NumH.is_null())
        && (data.NumHfixed.is_null() || !table.NumHfixed.is_null())
        && (data.maxlen_iso_sort_key == 0 || !table.iso_sort_key.is_null())
        && (data.iso_exchg_atnos.is_null() || !table.iso_exchg_atnos.is_null())
}

fn free_con_table(
    heap: &mut SourceHeap,
    pointer: SourceMutPointer<ConTable>,
) -> Result<(), SourceHeapError> {
    if pointer.is_null() {
        return Ok(());
    }
    let mut table = source_clone(heap, pointer, 0)?;
    CTableFree(heap, Some(&mut table))?;
    inchi_free(heap, pointer)
}

fn sync_partition_stack(
    heap: &mut SourceHeap,
    stack: SourceMutPointer<SourceMutPointer<AT_RANK>>,
    partitions: &[Partition],
) -> Result<(), SourceHeapError> {
    for (index, partition) in partitions.iter().enumerate() {
        source_set(heap, stack, (index * 2) as i32, partition.Rank)?;
        source_set(heap, stack, (index * 2 + 1) as i32, partition.AtNumber)?;
    }
    Ok(())
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn GetBaseCanonRanking(
    heap: &mut SourceHeap,
    ic: &mut INCHI_CLOCK,
    num_atoms: i32,
    mut num_at_tg: i32,
    at: [SourceMutPointer<sp_ATOM>; TAUT_NUM as usize],
    mut t_group_info: Option<&mut T_GROUP_INFO>,
    sizes: &[ATOM_SIZES; TAUT_NUM as usize],
    pBCN: &mut BCN,
    ulTimeOutTime: SourceMutPointer<inchiTime>,
    pCG: &mut CANON_GLOBALS,
    bFixIsoFixedH: i32,
    LargeMolecules: i32,
    user_action: Option<fn() -> u32>,
    console_quit: Option<fn() -> i32>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:5105 GetBaseCanonRanking
    // INCHI✔️❌: int GetBaseCanonRanking( INCHI_CLOCK *ic,
    // INCHI✔️❌:                          int num_atoms,
    // INCHI✔️❌:                          int num_at_tg,
    // INCHI✔️❌:                          sp_ATOM* at[],
    // INCHI✔️❌:                          T_GROUP_INFO *t_group_info,
    // INCHI✔️❌:                          ATOM_SIZES s[],
    // INCHI✔️❌:                          BCN *pBCN,
    // INCHI✔️❌:                          struct tagInchiTime *ulTimeOutTime,
    // INCHI✔️❌:                          CANON_GLOBALS *pCG,
    // INCHI✔️❌:                          int bFixIsoFixedH,
    // INCHI✔️❌:                          int LargeMolecules )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0;
    // INCHI✔️❌:     int iBase;                   /* base structure index, always valid; = TAUT_YES except special fully non-taut mode */
    // INCHI✔️❌:     int iOther;                  /* other than basic structure index, usually non-taut; may be = iBase */
    // INCHI✔️❌:     int bReqNonTaut;             /* 1 => requested non-tautomeric results */
    // INCHI✔️❌:     int bReqTaut;                /* 1 => requested tautomeric results and the base structure is tautomeric */
    // INCHI✔️❌:     int bChanged;
    // INCHI✔️❌:
    // INCHI✔️❌:     sp_ATOM *at_base = NULL;
    // INCHI✔️❌:     sp_ATOM *at_other = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     int bTautIgnoreIsotopic = 0;
    // INCHI✔️❌:     /*int bIgnoreIsotopic     = 0;*/
    // INCHI✔️❌:     int nNumCurrRanks = 0;
    // INCHI✔️❌:     int nMaxLenRankStack = 0;
    // INCHI✔️❌:     int num_max = num_at_tg;
    // INCHI✔️❌:     long lCount;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* local allocations */
    // INCHI✔️❌:     ATOM_INVARIANT2 *pAtomInvariant = NULL;
    // INCHI✔️❌:     NEIGH_LIST     *NeighList[TAUT_NUM];
    // INCHI✔️❌:     ConTable *Ct_Temp = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* initial partition for canonicalization */
    // INCHI✔️❌:     AT_RANK *nRank = NULL;
    // INCHI✔️❌:     AT_NUMB *nAtomNumber = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* canonicalization output */
    // INCHI✔️❌:
    // INCHI✔️❌:     ConTable *Ct_NoH = NULL;
    // INCHI✔️❌:     AT_RANK *nCanonRankNoH = NULL;
    // INCHI✔️❌:     AT_NUMB *nAtomNumberCanonNoH = NULL;
    // INCHI✔️❌:     AT_RANK *nSymmRankNoH = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     ConTable *Ct_NoTautH = NULL;
    // INCHI✔️❌:     AT_RANK *nSymmRankNoTautH = NULL;
    // INCHI✔️❌:     AT_RANK *nCanonRankNoTautH = NULL;
    // INCHI✔️❌:     AT_NUMB *nAtomNumberCanonNoTautH = NULL;
    // INCHI✔️❌:     NUM_H   *numHNoTautH = NULL;
    // INCHI✔️❌:     int      lenNumHNoTautH;
    // INCHI✔️❌:     int      maxlenNumHNoTautH;
    // INCHI✔️❌:
    // INCHI✔️❌:     ConTable *Ct_Base = NULL;
    // INCHI✔️❌:     AT_RANK *nSymmRankBase = NULL;
    // INCHI✔️❌:     AT_RANK *nCanonRankBase = NULL;
    // INCHI✔️❌:     AT_NUMB *nAtomNumberCanonBase = NULL;
    // INCHI✔️❌:     NUM_H   *numH = NULL;
    // INCHI✔️❌:     int      lenNumH = num_atoms;
    // INCHI✔️❌:     int      maxlenNumH = 0;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_AUX_RANKING == 1 )
    // INCHI✔️❌:     AT_RANK *nRankAux = NULL;
    // INCHI✔️❌:     AT_NUMB *nAtomNumberAux = NULL;
    // INCHI✔️❌:     ATOM_INVARIANT2 *pAtomInvariantAux = NULL;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     ConTable *Ct_FixH = NULL;
    // INCHI✔️❌:     AT_RANK *nSymmRankFixH = NULL;
    // INCHI✔️❌:     AT_RANK *nCanonRankFixH = NULL;
    // INCHI✔️❌:     AT_NUMB *nAtomNumberCanonFixH = NULL;
    // INCHI✔️❌:     NUM_H   *NumHfixed = NULL;
    // INCHI✔️❌:     int      maxlenNumHfixed;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* isotopic canonicalization */
    // INCHI✔️❌:
    // INCHI✔️❌:     ConTable *Ct_NoTautHIso = NULL;
    // INCHI✔️❌:     AT_RANK *nSymmRankNoTautHIso = NULL;
    // INCHI✔️❌:     AT_RANK *nCanonRankNoTautHIso = NULL;
    // INCHI✔️❌:     AT_NUMB *nAtomNumberCanonNoTautHIso = NULL;
    // INCHI✔️❌:     AT_ISO_SORT_KEY *iso_sort_key_NoTautH = NULL;
    // INCHI✔️❌:     int              maxlen_iso_sort_key_NoTautH = 0;
    // INCHI✔️❌:     int              len_iso_sort_key_NoTautH = 0;
    // INCHI✔️❌:     int num_iso_NoTautH = 0, num_iso_NoAuxBase;
    // INCHI✔️❌:
    // INCHI✔️❌:     ConTable *Ct_BaseIso = NULL;
    // INCHI✔️❌:     AT_RANK *nSymmRankBaseIso = NULL;
    // INCHI✔️❌:     AT_RANK *nCanonRankBaseIso = NULL;
    // INCHI✔️❌:     AT_NUMB *nAtomNumberCanonBaseIso = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     AT_ISO_SORT_KEY *iso_sort_keyBase = NULL;
    // INCHI✔️❌:     int              maxlen_iso_sort_keyBase = 0;
    // INCHI✔️❌:     int              len_iso_sort_keyBase = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     int              bUseIsoAuxBase[TAUT_NUM];
    // INCHI✔️❌:     S_CHAR          *iso_exchg_atnos = NULL;
    // INCHI✔️❌:     int              len_iso_exchg_atnos = 0;
    // INCHI✔️❌:     int              maxlen_iso_exchg_atnos = 0;
    // INCHI✔️❌:     int num_iso_Base = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     AT_ISO_SORT_KEY  iso_sort_key;
    // INCHI✔️❌:
    // INCHI✔️❌:     ConTable *Ct_FixHIso = NULL;
    // INCHI✔️❌:     AT_RANK *nSymmRankFixHIso = NULL;
    // INCHI✔️❌:     AT_RANK *nCanonRankFixHIso = NULL;
    // INCHI✔️❌:     AT_NUMB *nAtomNumberCanonFixHIso = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:     AT_ISO_SORT_KEY  iso_sort_key2;
    // INCHI✔️❌:     AT_ISO_SORT_KEY *iso_sort_key_Hfixed = NULL;
    // INCHI✔️❌:     int              maxlen_iso_sort_key_Hfixed;
    // INCHI✔️❌:     int              len_iso_sort_key_Hfixed;
    // INCHI✔️❌:     int num_iso_Hfixed;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     AT_RANK *nTempRank = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     CANON_DATA    pCD[3]; /* = &CanonData; */
    // INCHI✔️❌:     CANON_COUNTS  CanonCounts;
    // INCHI✔️❌:     CANON_COUNTS *pCC = &CanonCounts;
    // INCHI✔️❌:
    // INCHI✔️❌:     int i, j, k, m;
    // INCHI✔️❌:     int nCanonFlags[2];
    // INCHI✔️❌:
    // INCHI✔️❌:     /* */
    // INCHI✔️❌:     int iflag;
    // INCHI✔️❌:
    // INCHI✔️❌:     memset( pCD, 0, sizeof( pCD ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( pCC, 0, sizeof( pCC[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( bUseIsoAuxBase, 0, sizeof( bUseIsoAuxBase ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( nCanonFlags, 0, sizeof( nCanonFlags ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     NeighList[TAUT_NON] = NULL;
    // INCHI✔️❌:     NeighList[TAUT_YES] = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* select base structure, find whether it is tautomeric or not */
    // INCHI✔️❌:     /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     if (at[TAUT_YES] && s[TAUT_YES].nLenCT &&
    // INCHI✔️❌:          t_group_info && ( s[TAUT_YES].nLenLinearCTTautomer > 0 && /* ordinary tautomerism */
    // INCHI✔️❌:              ((t_group_info->t_group && t_group_info->num_t_groups > 0) ||
    // INCHI✔️❌:              /* protons have been moved */
    // INCHI✔️❌:              ( t_group_info->tni.bNormalizationFlags & FLAG_NORM_CONSIDER_TAUT ) ||
    // INCHI✔️❌:              /* tautomerism due to possible isotopic proton exchange */
    // INCHI✔️❌:              (t_group_info->nNumIsotopicEndpoints > 1 &&
    // INCHI✔️❌:              ( t_group_info->bTautFlagsDone & ( TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE ) ))) ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* tautomeric: (1) has tautomeric atoms OR
    // INCHI✔️❌:                        (2) H-atoms have been rearranged due to proton addition/removal OR
    // INCHI✔️❌:                        (3) Found isotopic H-atoms on tautomeric or hetero atoms
    // INCHI✔️❌:          */
    // INCHI✔️❌:         iBase = TAUT_YES;
    // INCHI✔️❌:         bReqTaut = 1;
    // INCHI✔️❌:         bUseIsoAuxBase[iBase] = ( s[iBase].nLenIsotopicEndpoints > 1 ) &&
    // INCHI✔️❌:             ( t_group_info->bTautFlagsDone & ( TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE ) );
    // INCHI✔️❌:         if (at[TAUT_NON] && s[TAUT_NON].nLenCT)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             iOther = TAUT_NON; /* tautomeric and non-tautomeric */
    // INCHI✔️❌:             bReqNonTaut = 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             iOther = iBase; /* tautomeric only */
    // INCHI✔️❌:             bReqNonTaut = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     else if (at[TAUT_NON] && s[TAUT_NON].nLenCT)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* force pure non-tautomeric processing; happens for testing only */
    // INCHI✔️❌:         iBase = TAUT_NON;
    // INCHI✔️❌:         bReqTaut = 0;
    // INCHI✔️❌:         iOther = iBase;
    // INCHI✔️❌:         bReqNonTaut = 1;
    // INCHI✔️❌:         num_at_tg = num_atoms;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     else if (at[TAUT_YES] && s[TAUT_YES].nLenCT)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* although the user requested tautomeric processing, tautomerism has not been found */
    // INCHI✔️❌:         /* however, the results should be saved in the TAUT_YES elements of the arrays */
    // INCHI✔️❌:         iBase = TAUT_YES;
    // INCHI✔️❌:         bReqTaut = 0;
    // INCHI✔️❌:         bUseIsoAuxBase[iBase] = ( s[iBase].nLenIsotopicEndpoints > 1 );
    // INCHI✔️❌:         iOther = iBase;
    // INCHI✔️❌:         bReqNonTaut = 1;
    // INCHI✔️❌:         num_at_tg = num_atoms;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = CT_UNKNOWN_ERR;
    // INCHI✔️❌:         goto exit_error;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (bReqTaut)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* save "process isotopic" mark; temporarily set it to NO */
    // INCHI✔️❌:         bTautIgnoreIsotopic = t_group_info->bIgnoreIsotopic;
    // INCHI✔️❌:         t_group_info->bIgnoreIsotopic = 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:     num_iso_Hfixed =
    // INCHI✔️❌:         len_iso_sort_key_Hfixed =
    // INCHI✔️❌:         maxlen_iso_sort_key_Hfixed = 0;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     /* prepare initial data */
    // INCHI✔️❌:     at_base = at[iBase];
    // INCHI✔️❌:     at_other = at[iOther];
    // INCHI✔️❌:     pAtomInvariant = (ATOM_INVARIANT2 *) inchi_calloc( num_max, sizeof( pAtomInvariant[0] ) );
    // INCHI✔️❌:     nSymmRankNoH = (AT_RANK *) inchi_calloc( num_max, sizeof( nSymmRankNoH[0] ) );
    // INCHI✔️❌:     nCanonRankNoH = (AT_RANK *) inchi_calloc( num_max, sizeof( nCanonRankNoH[0] ) );
    // INCHI✔️❌:     nAtomNumberCanonNoH = (AT_NUMB *) inchi_calloc( num_max, sizeof( nAtomNumberCanonNoH[0] ) );
    // INCHI✔️❌:     nRank = (AT_RANK *) inchi_calloc( num_max, sizeof( nRank[0] ) );
    // INCHI✔️❌:     nAtomNumber = (AT_NUMB *) inchi_calloc( num_max, sizeof( nAtomNumber[0] ) );
    // INCHI✔️❌:     nTempRank = (AT_RANK *) inchi_calloc( num_max, sizeof( nTempRank[0] ) );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!pAtomInvariant ||
    // INCHI✔️❌:          !nSymmRankNoH || !nCanonRankNoH || !nAtomNumberCanonNoH ||
    // INCHI✔️❌:          !nRank || !nAtomNumber || !nTempRank)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_error_alloc;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_AUX_RANKING == 1 )
    // INCHI✔️❌:     nRankAux = (AT_RANK *) inchi_calloc( num_max, sizeof( nRankAux[0] ) );
    // INCHI✔️❌:     nAtomNumberAux = (AT_NUMB *) inchi_calloc( num_max, sizeof( nAtomNumberAux[0] ) );
    // INCHI✔️❌:     pAtomInvariantAux = (ATOM_INVARIANT2 *) inchi_malloc( num_max * sizeof( pAtomInvariantAux[0] ) );
    // INCHI✔️❌:     if (!nRankAux || !nAtomNumberAux || !pAtomInvariantAux)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_error_alloc;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     if (bReqTaut)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (!( NeighList[TAUT_YES] =
    // INCHI✔️❌:                CreateNeighList( num_atoms, num_at_tg, at_base, 0, t_group_info ) ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_error_alloc;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* needed for the hydrogenless structure */
    // INCHI✔️❌:         if (!( NeighList[TAUT_NON] =
    // INCHI✔️❌:                CreateNeighList( num_atoms, num_atoms, at_base, 0, NULL ) ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_error_alloc;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (!( NeighList[TAUT_NON] =
    // INCHI✔️❌:                CreateNeighList( num_atoms, num_atoms, at_base, 0, NULL ) ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_error_alloc;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         NeighList[TAUT_YES] = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:         INCHI_HEAPCHK
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* avoid memory leaks in case of error */
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     pBCN->ftcn[TAUT_NON].NeighList          = NeighList[TAUT_NON];
    // INCHI✔️❌:     pBCN->ftcn[TAUT_YES].NeighList          = NeighList[TAUT_YES];
    // INCHI✔️❌:     */
    // INCHI✔️❌:     pBCN->nMaxLenRankStack = 0;
    // INCHI✔️❌:     pBCN->num_max = num_max;        /* allocated nRank[] arrays lengths in pRankStack */
    // INCHI✔️❌:     pBCN->num_at_tg = num_at_tg;  /* all of the following arrays have this length */
    // INCHI✔️❌:     pBCN->num_atoms = num_atoms;
    // INCHI✔️❌:     pBCN->ulTimeOutTime = ulTimeOutTime;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* initial partitioning of a hydrogenless skeleton: fill out the inveriant */
    // INCHI✔️❌:     FillOutAtomInvariant2( at_base,
    // INCHI✔️❌:                            num_atoms,
    // INCHI✔️❌:                            num_atoms,
    // INCHI✔️❌:                            pAtomInvariant,
    // INCHI✔️❌:                            1 /*bIgnoreIsotopic*/,
    // INCHI✔️❌:                            0 /*bHydrogensInRanks*/,
    // INCHI✔️❌:                            0 /*bHydrogensFixedInRanks*/,
    // INCHI✔️❌:                            0 /*bTaut=bDigraph*/,
    // INCHI✔️❌:                            0 /* bTautGroupsOnly */,
    // INCHI✔️❌:                            NULL /*t_group_info*/ );
    // INCHI✔️❌:
    // INCHI✔️❌:     /* initial partitioning of a hydrogenless skeleton: create equitable partition (assign initial ranks) */
    // INCHI✔️❌:     nNumCurrRanks = SetInitialRanks2( num_atoms, pAtomInvariant, nRank, nAtomNumber, pCG );
    // INCHI✔️❌:
    // INCHI✔️❌:     lCount = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* make equitable partition in pBCN->pRankStack[0,1] */
    // INCHI✔️❌:     nNumCurrRanks = DifferentiateRanks2( pCG,
    // INCHI✔️❌:                                          num_atoms,
    // INCHI✔️❌:                                          NeighList[TAUT_NON],
    // INCHI✔️❌:                                          nNumCurrRanks, nRank,
    // INCHI✔️❌:                                          nTempRank,
    // INCHI✔️❌:                                          nAtomNumber,
    // INCHI✔️❌:                                          &lCount,
    // INCHI✔️❌:                                          0 /* 0 means use qsort */ );
    // INCHI✔️❌:
    // INCHI✔️❌:     /* allocate partition stack */
    // INCHI✔️❌:     nMaxLenRankStack = 2 * ( num_at_tg - nNumCurrRanks ) + 8;  /* was 2*(...) + 6 */
    // INCHI✔️❌:     pBCN->pRankStack = (AT_RANK **) inchi_calloc( nMaxLenRankStack, sizeof( pBCN->pRankStack[0] ) );
    // INCHI✔️❌:     if (!pBCN->pRankStack)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         pBCN->nMaxLenRankStack = 0; /* avoid memory leaks in case of error */
    // INCHI✔️❌:         goto exit_error_alloc;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     pBCN->nMaxLenRankStack = nMaxLenRankStack; /* avoid memory leaks in case of error */
    // INCHI✔️❌:     /* init partition stack */
    // INCHI✔️❌:     pBCN->pRankStack[0] = nRank;
    // INCHI✔️❌:     pBCN->pRankStack[1] = nAtomNumber;
    // INCHI✔️❌:
    // INCHI✔️❌:     /********************************************************************************************/
    // INCHI✔️❌:     /* get NoH/no taut groups  canonical numbering, connection table, and equivalence partition */
    // INCHI✔️❌:     /********************************************************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:     /* pointers */
    // INCHI✔️❌:     pCD[iOther].LinearCT = NULL;
    // INCHI✔️❌:     pCD[iOther].NumH = NULL;
    // INCHI✔️❌:     pCD[iOther].NumHfixed = NULL;
    // INCHI✔️❌:     pCD[iOther].iso_sort_key = NULL;
    // INCHI✔️❌:     pCD[iOther].iso_exchg_atnos = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:     pCD[iOther].iso_sort_key_Hfixed = NULL;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     /* variables - unchanged */
    // INCHI✔️❌:     pCD[iOther].ulTimeOutTime = pBCN->ulTimeOutTime;
    // INCHI✔️❌:     pCD[iOther].nMaxLenLinearCT = s[iOther].nLenCTAtOnly + 1;
    // INCHI✔️❌:     /* return values & input/output */
    // INCHI✔️❌:     pCD[iOther].nLenLinearCT = s[iOther].nLenCTAtOnly;
    // INCHI✔️❌:     pCD[iOther].nLenCTAtOnly = s[iOther].nLenCTAtOnly;
    // INCHI✔️❌:     pCD[iOther].lenNumH = 0;
    // INCHI✔️❌:     pCD[iOther].lenNumHfixed = 0;
    // INCHI✔️❌:     pCD[iOther].len_iso_sort_key = 0;
    // INCHI✔️❌:     pCD[iOther].maxlen_iso_sort_key = 0;
    // INCHI✔️❌:     pCD[iOther].len_iso_exchg_atnos = 0;
    // INCHI✔️❌:     pCD[iOther].maxlen_iso_exchg_atnos = 0;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:     pCD[iOther].len_iso_sort_key_Hfixed = 0;
    // INCHI✔️❌:     pCD[iOther].maxlen_iso_sort_key_Hfixed = 0;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = CanonGraph01( ic, pCG, num_atoms, num_atoms, num_max, 0,
    // INCHI✔️❌:                         NeighList[TAUT_NON], (Partition *) pBCN->pRankStack,
    // INCHI✔️❌:                         nSymmRankNoH, nCanonRankNoH, nAtomNumberCanonNoH,
    // INCHI✔️❌:                         pCD + iOther, pCC, NULL, &Ct_NoH, LargeMolecules );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ret < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_error;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* update initial partitioning */
    // INCHI✔️❌:     nNumCurrRanks = FixCanonEquivalenceInfo( pCG, num_atoms, nSymmRankNoH, nRank, nTempRank, nAtomNumber, &bChanged ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* repartition if necessary */
    // INCHI✔️❌:     if (bChanged & 3)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Ct_NoH)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             CTableFree( Ct_NoH );
    // INCHI✔️❌:             inchi_free( Ct_NoH );
    // INCHI✔️❌:             Ct_NoH = NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         pCD[iOther].nCanonFlags |= CANON_FLAG_NO_H_RECANON;
    // INCHI✔️❌:
    // INCHI✔️❌:         ret = CanonGraph02( ic, pCG, num_atoms, num_atoms, num_max, 0,
    // INCHI✔️❌:                             NeighList[TAUT_NON], (Partition *) pBCN->pRankStack,
    // INCHI✔️❌:                             nSymmRankNoH, nCanonRankNoH, nAtomNumberCanonNoH,
    // INCHI✔️❌:                             pCD + iOther, pCC, NULL, &Ct_NoH, LargeMolecules );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (ret < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_error;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /********************************************************************************/
    // INCHI✔️❌:     /* get NoTautH canonical numbering, connection table, and equivalence partition */
    // INCHI✔️❌:     /********************************************************************************/
    // INCHI✔️❌:     maxlenNumHNoTautH = num_atoms + 1;
    // INCHI✔️❌:     nSymmRankNoTautH = (AT_RANK *) inchi_calloc( num_max, sizeof( nSymmRankNoTautH[0] ) );
    // INCHI✔️❌:     nCanonRankNoTautH = (AT_RANK *) inchi_calloc( num_max, sizeof( nCanonRankNoTautH[0] ) );
    // INCHI✔️❌:     nAtomNumberCanonNoTautH = (AT_NUMB *) inchi_calloc( num_max, sizeof( nAtomNumberCanonNoTautH[0] ) );
    // INCHI✔️❌:     numHNoTautH = (NUM_H *) inchi_calloc( maxlenNumHNoTautH, sizeof( numHNoTautH[0] ) );
    // INCHI✔️❌:     if (!numHNoTautH || !nSymmRankNoTautH || !nCanonRankNoTautH || !nAtomNumberCanonNoTautH)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_error_alloc;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* find number of H atoms attached to not-a-tautomeric-endpoint atoms */
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         numHNoTautH[i] = ( !at_base[i].endpoint && at_base[i].num_H ) ? at_base[i].num_H + BASE_H_NUMBER : EMPTY_H_NUMBER;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* pointers */
    // INCHI✔️❌:     pCD[iOther].LinearCT = NULL;
    // INCHI✔️❌:     pCD[iOther].NumH = numHNoTautH;
    // INCHI✔️❌:     pCD[iOther].NumHfixed = NULL;
    // INCHI✔️❌:     pCD[iOther].iso_sort_key = NULL;
    // INCHI✔️❌:     pCD[iOther].iso_exchg_atnos = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:     pCD[iOther].iso_sort_key_Hfixed = NULL;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     /* variables - unchanged */
    // INCHI✔️❌:     pCD[iOther].ulTimeOutTime = pBCN->ulTimeOutTime;
    // INCHI✔️❌:     pCD[iOther].nMaxLenLinearCT = s[iOther].nLenCTAtOnly + 1;
    // INCHI✔️❌:     pCD[iOther].maxlenNumH = maxlenNumHNoTautH;
    // INCHI✔️❌:     /* return values & input/output */
    // INCHI✔️❌:     pCD[iOther].nLenLinearCT = s[iOther].nLenCTAtOnly;
    // INCHI✔️❌:     pCD[iOther].nLenCTAtOnly = s[iOther].nLenCTAtOnly;
    // INCHI✔️❌:     pCD[iOther].lenNumH = lenNumHNoTautH = num_atoms;
    // INCHI✔️❌:     pCD[iOther].lenNumHfixed = 0;
    // INCHI✔️❌:     pCD[iOther].len_iso_sort_key = 0;
    // INCHI✔️❌:     pCD[iOther].maxlen_iso_sort_key = 0;
    // INCHI✔️❌:     pCD[iOther].len_iso_exchg_atnos = 0;
    // INCHI✔️❌:     pCD[iOther].maxlen_iso_exchg_atnos = 0;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:     pCD[iOther].len_iso_sort_key_Hfixed = 0;
    // INCHI✔️❌:     pCD[iOther].maxlen_iso_sort_key_Hfixed = 0;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     pCD[iOther].nAuxRank = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* check whether we need NoTautH cononicalization */
    // INCHI✔️❌:     memset( nTempRank, 0, num_max * sizeof( nTempRank[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (nTempRank[nSymmRankNoH[i] - 1] < i)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nTempRank[nSymmRankNoH[i] - 1] = i; /* greatest class representative */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (numHNoTautH[i] != numHNoTautH[nTempRank[nSymmRankNoH[i] - 1]])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pCD[iOther].nCanonFlags |= CANON_FLAG_NO_TAUT_H_DIFF;
    // INCHI✔️❌:             break; /* atoms so far found to be equivalent have different number of H; the canonicalization is needed */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* i = 0; *//* debug: force to call the canonicalization */
    // INCHI✔️❌:     if (i < num_atoms)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* needs canonicalization */
    // INCHI✔️❌:         /* get aux canonical ranking of the structure with attached H */
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_AUX_RANKING == 1 )
    // INCHI✔️❌:         /* refine no-H partition according to not-a-taut-H distribution */
    // INCHI✔️❌:         memset( pAtomInvariantAux, 0, num_max * sizeof( pAtomInvariantAux[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pAtomInvariantAux[i].val[0] = nSymmRankNoH[i];
    // INCHI✔️❌:             pAtomInvariantAux[i].val[1] = numHNoTautH[i]; /* additional differentiation: not-a-taut-H distribution */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* initial partitioning */
    // INCHI✔️❌:         nNumCurrRanks = SetInitialRanks2( num_atoms, pAtomInvariantAux, nRankAux, nAtomNumberAux, pCG );
    // INCHI✔️❌:
    // INCHI✔️❌:         /* make equitable partition */
    // INCHI✔️❌:         nNumCurrRanks = DifferentiateRanks2( pCG, num_atoms, NeighList[TAUT_NON],
    // INCHI✔️❌:                                             nNumCurrRanks, nRankAux,
    // INCHI✔️❌:                                             nTempRank, nAtomNumberAux, &lCount, 0 /* 0 means use qsort */ ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:
    // INCHI✔️❌:         /* to accelerate do not call CanonGraph() to find really equivalent atoms */
    // INCHI✔️❌:         pCD[iOther].nAuxRank = nRankAux;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         ret = CanonGraph03( ic, pCG, num_atoms, num_atoms, num_max, 1 /* digraph?? was 0 */,
    // INCHI✔️❌:                             NeighList[TAUT_NON], (Partition *) pBCN->pRankStack,
    // INCHI✔️❌:                             nSymmRankNoTautH, nCanonRankNoTautH, nAtomNumberCanonNoTautH,
    // INCHI✔️❌:                             pCD + iOther, pCC, &Ct_NoH, &Ct_NoTautH, LargeMolecules );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (ret < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_error;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* in case of non-tautomeric structure the final results are in:
    // INCHI✔️❌:
    // INCHI✔️❌:                    nSymmRankNoTautH
    // INCHI✔️❌:                    nCanonRankNoTautH
    // INCHI✔️❌:                    nAtomNumberCanonNoTautH
    // INCHI✔️❌:                    Ct_NoTautH
    // INCHI✔️❌:                    numHNoTautH (original H positions)
    // INCHI✔️❌:         */
    // INCHI✔️❌:     } /* if ( i < num_atoms )  */
    // INCHI✔️❌:
    // INCHI✔️❌:     else
    // INCHI✔️❌:
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         /* copy the results of the previous (no H) canonicalization */
    // INCHI✔️❌:         /* in this case numHNoTautH[] is not needed for the next canonicalization(s) */
    // INCHI✔️❌:         if (( Ct_Temp = (ConTable *) inchi_calloc( 1, sizeof( *Ct_Temp ) ) ) &&
    // INCHI✔️❌:              CTableCreate( Ct_Temp, num_atoms, pCD + iOther ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             CtFullCopy( Ct_Temp, Ct_NoH );
    // INCHI✔️❌:             /* since Ct_NoH does not have Ct_NoH->NumH we have to fill out Ct_Temp->NumH separately */
    // INCHI✔️❌:             for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 Ct_Temp->NumH[nCanonRankNoH[i] - 1] = numHNoTautH[i];
    // INCHI✔️❌:                 /*Ct_Temp->NumH[i] = numHNoTautH[nAtomNumberCanonNoH[i]]; -- alternative */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             Ct_Temp->lenNumH = num_atoms;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_error_alloc;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         Ct_NoTautH = Ct_Temp;
    // INCHI✔️❌:         Ct_Temp = NULL;
    // INCHI✔️❌:         /* djb-rwth: functions replaced with their safe C11 variants */
    // INCHI✔️❌:         memcpy(nSymmRankNoTautH, nSymmRankNoH, num_atoms * sizeof(nSymmRankNoTautH[0]));
    // INCHI✔️❌:         memcpy(nCanonRankNoTautH, nCanonRankNoH, num_atoms * sizeof(nCanonRankNoTautH[0]));
    // INCHI✔️❌:         memcpy(nAtomNumberCanonNoTautH, nAtomNumberCanonNoH, num_atoms * sizeof(nAtomNumberCanonNoTautH[0]));
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* in case of non-tautomeric component this is the final result */
    // INCHI✔️❌:     /* i = CtFullCompare( Ct_NoTautH, Ct_Temp, num_atoms, 0, 0 );*/
    // INCHI✔️❌:
    // INCHI✔️❌:     /*******************************************************************************************/
    // INCHI✔️❌:     /* If only Isotopic atoms and isotopic H, tautomerism has not been found:                  */
    // INCHI✔️❌:     /* get isotopic canonical numbering, connection table, and equivalence partition           */
    // INCHI✔️❌:     /*******************************************************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:     if (s[iOther].num_isotopic_atoms && !s[iOther].bIgnoreIsotopic && !bReqTaut && bReqNonTaut)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         maxlen_iso_sort_key_NoTautH = num_atoms + 1;
    // INCHI✔️❌:         nSymmRankNoTautHIso = (AT_RANK *) inchi_calloc( num_max, sizeof( nSymmRankNoTautHIso[0] ) );
    // INCHI✔️❌:         nCanonRankNoTautHIso = (AT_RANK *) inchi_calloc( num_max, sizeof( nCanonRankNoTautHIso[0] ) );
    // INCHI✔️❌:         nAtomNumberCanonNoTautHIso = (AT_NUMB *) inchi_calloc( num_max, sizeof( nAtomNumberCanonNoTautHIso[0] ) );
    // INCHI✔️❌:         iso_sort_key_NoTautH = (AT_ISO_SORT_KEY *) inchi_calloc( maxlen_iso_sort_key_NoTautH, sizeof( iso_sort_key_NoTautH[0] ) );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!nSymmRankNoTautHIso || !nCanonRankNoTautHIso || !nAtomNumberCanonNoTautHIso || !iso_sort_key_NoTautH)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_error_alloc;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* fill out isotopic non-tautomeric keys */
    // INCHI✔️❌:         num_iso_NoTautH = 0;
    // INCHI✔️❌:         for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (at_base[i].endpoint)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* should not happen */
    // INCHI✔️❌:                 iso_sort_key = make_iso_sort_key( at_base[i].iso_atw_diff, 0, 0, 0 );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 iso_sort_key = make_iso_sort_key( at_base[i].iso_atw_diff, at_base[i].num_iso_H[0], at_base[i].num_iso_H[1], at_base[i].num_iso_H[2] );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (iso_sort_key)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 iso_sort_key_NoTautH[i] = iso_sort_key;
    // INCHI✔️❌:                 num_iso_NoTautH++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 iso_sort_key_NoTautH[i] = EMPTY_ISO_SORT_KEY;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* pointers */
    // INCHI✔️❌:         pCD[iOther].LinearCT = NULL; /* LinearCT; */
    // INCHI✔️❌:         pCD[iOther].NumH = numHNoTautH;
    // INCHI✔️❌:         pCD[iOther].NumHfixed = NULL;
    // INCHI✔️❌:         pCD[iOther].iso_sort_key = iso_sort_key_NoTautH;
    // INCHI✔️❌:         pCD[iOther].iso_exchg_atnos = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:         pCD[iOther].iso_sort_key_Hfixed = NULL;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         /* variables - unchanged */
    // INCHI✔️❌:         pCD[iOther].ulTimeOutTime = pBCN->ulTimeOutTime;
    // INCHI✔️❌:         pCD[iOther].nMaxLenLinearCT = s[iOther].nLenCTAtOnly + 1;
    // INCHI✔️❌:         pCD[iOther].maxlenNumH = maxlenNumHNoTautH;
    // INCHI✔️❌:         /* return values & input/output */
    // INCHI✔️❌:         pCD[iOther].nLenLinearCT = s[iOther].nLenCTAtOnly;
    // INCHI✔️❌:         pCD[iOther].nLenCTAtOnly = s[iOther].nLenCTAtOnly;
    // INCHI✔️❌:         pCD[iOther].lenNumH = lenNumHNoTautH /*= num_atoms*/;
    // INCHI✔️❌:         pCD[iOther].lenNumHfixed = 0;
    // INCHI✔️❌:         pCD[iOther].len_iso_sort_key = len_iso_sort_key_NoTautH = num_atoms;
    // INCHI✔️❌:         pCD[iOther].maxlen_iso_sort_key = maxlen_iso_sort_key_NoTautH;
    // INCHI✔️❌:         pCD[iOther].len_iso_exchg_atnos = 0;
    // INCHI✔️❌:         pCD[iOther].maxlen_iso_exchg_atnos = 0;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:         pCD[iOther].len_iso_sort_key_Hfixed = 0;
    // INCHI✔️❌:         pCD[iOther].maxlen_iso_sort_key_Hfixed = 0;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         pCD[iOther].nAuxRank = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:         if (num_iso_NoTautH)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* check whether we need NoTautH cononicalization */
    // INCHI✔️❌:             memset( nTempRank, 0, num_max * sizeof( nTempRank[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:             for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (nTempRank[nSymmRankNoTautH[i] - 1] < i)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nTempRank[nSymmRankNoTautH[i] - 1] = i; /* greatest class representative */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (iso_sort_key_NoTautH[i] != iso_sort_key_NoTautH[nTempRank[nSymmRankNoTautH[i] - 1]])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     pCD[iOther].nCanonFlags |= CANON_FLAG_ISO_ONLY_NON_TAUT_DIFF;
    // INCHI✔️❌:                     break; /* atoms so far found to be equivalent differ in isotopes; the canonicalization is needed */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             i = num_atoms;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* i = 0; *//* debug: force to call the canonicalization */
    // INCHI✔️❌:         if (i < num_atoms)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* we need canonicalization */
    // INCHI✔️❌:             /* get aux canonical ranking of the structure with isotopic non-tautomeric H */
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_AUX_RANKING == 1 )
    // INCHI✔️❌:             /* refine no-taut-H partition according to non-taut H isotopic distribution */
    // INCHI✔️❌:             memset( pAtomInvariantAux, 0, num_max * sizeof( pAtomInvariantAux[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:             for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pAtomInvariantAux[i].val[0] = nSymmRankNoTautH[i];
    // INCHI✔️❌:                 pAtomInvariantAux[i].iso_sort_key = iso_sort_key_NoTautH[i]; /* additional differentiation */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* initial ranks for non-taut H isotopic distribution */
    // INCHI✔️❌:             nNumCurrRanks = SetInitialRanks2( num_atoms, pAtomInvariantAux, nRankAux, nAtomNumberAux, pCG );
    // INCHI✔️❌:             /* make equitable */
    // INCHI✔️❌:             nNumCurrRanks = DifferentiateRanks2( pCG, num_atoms, NeighList[TAUT_NON],
    // INCHI✔️❌:                                                 nNumCurrRanks, nRankAux,
    // INCHI✔️❌:                                                 nTempRank, nAtomNumberAux, &lCount, 0 /* 0 means use qsort */ ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:             /* to accelerate do not call CanonGraph() to find really equivalent atoms */
    // INCHI✔️❌:             pCD[iOther].nAuxRank = nRankAux;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:             ret = CanonGraph04( ic, pCG,
    // INCHI✔️❌:                                 num_atoms, num_atoms, num_max, 1 /* digraph?? was 0 */,
    // INCHI✔️❌:                                 NeighList[TAUT_NON], (Partition *) pBCN->pRankStack,
    // INCHI✔️❌:                                 nSymmRankNoTautHIso, nCanonRankNoTautHIso, nAtomNumberCanonNoTautHIso,
    // INCHI✔️❌:                                 pCD + iOther, pCC, &Ct_NoTautH, &Ct_NoTautHIso, LargeMolecules );
    // INCHI✔️❌:             if (ret < 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto exit_error;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* in case of non-tautomeric structure the final results are in:
    // INCHI✔️❌:
    // INCHI✔️❌:                        nSymmRankNoTautHIso
    // INCHI✔️❌:                        nCanonRankNoTautHIso
    // INCHI✔️❌:                        nAtomNumberCanonNoTautHIso
    // INCHI✔️❌:                        Ct_NoTautHIso
    // INCHI✔️❌:                        iso_sort_key_NoTautH (original isotopic atom positions)
    // INCHI✔️❌:             */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* copy the results of the previous (no taut H) canonicalization */
    // INCHI✔️❌:             /* in this case numHNoTautH[] is not needed for the next canonicalization(s) */
    // INCHI✔️❌:             if (( Ct_Temp = (ConTable *) inchi_calloc( 1, sizeof( *Ct_Temp ) ) ) &&
    // INCHI✔️❌:                  CTableCreate( Ct_Temp, num_atoms, pCD + iOther ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 CtFullCopy( Ct_Temp, Ct_NoTautH );
    // INCHI✔️❌:                 /* since Ct_NoTautH does not have Ct_NoTautH->iso_sort_key we have to fill out Ct_Temp->iso_sort_key separately */
    // INCHI✔️❌:                 for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     Ct_Temp->iso_sort_key[nCanonRankNoTautH[i] - 1] = iso_sort_key_NoTautH[i];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 Ct_Temp->len_iso_sort_key = num_atoms;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto exit_error_alloc;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             Ct_NoTautHIso = Ct_Temp;
    // INCHI✔️❌:             Ct_Temp = NULL;
    // INCHI✔️❌:             memcpy(nSymmRankNoTautHIso, nSymmRankNoTautH, num_atoms * sizeof(nSymmRankNoTautHIso[0]));
    // INCHI✔️❌:             memcpy(nCanonRankNoTautHIso, nCanonRankNoTautH, num_atoms * sizeof(nCanonRankNoTautHIso[0]));
    // INCHI✔️❌:             memcpy(nAtomNumberCanonNoTautHIso, nAtomNumberCanonNoTautH, num_atoms * sizeof(nAtomNumberCanonNoTautHIso[0]));
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* in case of non-tautomeric component this is the final result */
    // INCHI✔️❌:         /* i = CtFullCompare( Ct_NoTautHIso, Ct_Temp, num_atoms, 0, 0 );*/
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (bReqTaut)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*****************************************************************************/
    // INCHI✔️❌:         /* Tautomeric Structure Canonicalizaton:                                     */
    // INCHI✔️❌:         /* get base canonical numbering, connection table, and equivalence partition */
    // INCHI✔️❌:         /*****************************************************************************/
    // INCHI✔️❌:         /* find H atoms attached to non-tautomeric-endpoints and to tautomeric endpoints */
    // INCHI✔️❌:         maxlenNumH = num_atoms + T_NUM_NO_ISOTOPIC*( num_at_tg - num_atoms ) + 1; /* including negative charges */
    // INCHI✔️❌:         nSymmRankBase = (AT_RANK *) inchi_calloc( num_max, sizeof( nSymmRankBase[0] ) );
    // INCHI✔️❌:         nCanonRankBase = (AT_RANK *) inchi_calloc( num_max, sizeof( nCanonRankBase[0] ) );
    // INCHI✔️❌:         nAtomNumberCanonBase = (AT_NUMB *) inchi_calloc( num_max, sizeof( nAtomNumberCanonBase[0] ) );
    // INCHI✔️❌:         numH = (NUM_H *) inchi_calloc( maxlenNumH, sizeof( numH[0] ) );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!numH || !nSymmRankBase || !nCanonRankBase || !nAtomNumberCanonBase)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_error_alloc;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* non-tautomeric H counts */
    // INCHI✔️❌:         for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             numH[i] = ( !at_base[i].endpoint && at_base[i].num_H ) ? at_base[i].num_H + BASE_H_NUMBER : EMPTY_H_NUMBER;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* tautomeric H and negative charge counts */
    // INCHI✔️❌:         for (i = k = num_atoms; i < num_at_tg; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             m = i - num_atoms;
    // INCHI✔️❌:             for (j = 0; j < T_NUM_NO_ISOTOPIC; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* non-zeroes for j=1 are negative charge counts; T_NUM_NO_ISOTOPIC=2 entry per t-group */
    // INCHI✔️❌:                 numH[k++] = t_group_info->t_group[m].num[j] ? t_group_info->t_group[m].num[j] + BASE_H_NUMBER : EMPTY_H_NUMBER;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* pointers */
    // INCHI✔️❌:         pCD[iBase].LinearCT = NULL;
    // INCHI✔️❌:         pCD[iBase].NumH = numH; /* num_atoms non-tautomeric H; num_tg pairs of H and (-) in t-groups */
    // INCHI✔️❌:         pCD[iBase].NumHfixed = NULL;
    // INCHI✔️❌:         pCD[iBase].iso_sort_key = NULL;
    // INCHI✔️❌:         pCD[iBase].iso_exchg_atnos = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:         pCD[iBase].iso_sort_key_Hfixed = NULL;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         /* variables - unchanged */
    // INCHI✔️❌:         pCD[iBase].ulTimeOutTime = pBCN->ulTimeOutTime;
    // INCHI✔️❌:         pCD[iBase].nMaxLenLinearCT = s[iBase].nLenCT + 1;
    // INCHI✔️❌:         pCD[iBase].maxlenNumH = maxlenNumH;
    // INCHI✔️❌:         /* return values & input/output */
    // INCHI✔️❌:         pCD[iBase].nLenLinearCT = s[iBase].nLenCT;
    // INCHI✔️❌:         pCD[iBase].nLenCTAtOnly = s[iBase].nLenCTAtOnly;
    // INCHI✔️❌:         pCD[iBase].lenNumH = lenNumH = k;
    // INCHI✔️❌:         pCD[iBase].lenNumHfixed = 0;
    // INCHI✔️❌:         pCD[iBase].len_iso_sort_key = 0;
    // INCHI✔️❌:         pCD[iBase].maxlen_iso_sort_key = 0;
    // INCHI✔️❌:         pCD[iBase].len_iso_exchg_atnos = 0;
    // INCHI✔️❌:         pCD[iBase].maxlen_iso_exchg_atnos = 0;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:         pCD[iBase].len_iso_sort_key_Hfixed = 0;
    // INCHI✔️❌:         pCD[iBase].maxlen_iso_sort_key_Hfixed = 0;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         pCD[iBase].nAuxRank = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:         /* make sure the initial partition is equitable (at this point t-groups do not have ranks yet) */
    // INCHI✔️❌:         FillOutAtomInvariant2( at_base,
    // INCHI✔️❌:                                num_atoms,
    // INCHI✔️❌:                                num_at_tg,
    // INCHI✔️❌:                                pAtomInvariant,
    // INCHI✔️❌:                                1 /*bIgnoreIsotopic*/,
    // INCHI✔️❌:                                0 /*bHydrogensInRanks*/,
    // INCHI✔️❌:                                0 /*bHydrogensFixedInRanks*/,
    // INCHI✔️❌:                                1 /*bTaut=bDigraph*/,
    // INCHI✔️❌:                                1 /* bTautGroupsOnly */,
    // INCHI✔️❌:                                t_group_info );
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pAtomInvariant[i].val[0] = pBCN->pRankStack[0][i];
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* initial ranks for t-group(s) only */
    // INCHI✔️❌:         nNumCurrRanks = SetInitialRanks2( num_at_tg, pAtomInvariant, nRank, nAtomNumber, pCG );
    // INCHI✔️❌:
    // INCHI✔️❌:         /* make equitable, call digraph procedure;
    // INCHI✔️❌:            pBCN->pRankStack[0] is nRank, pBCN->pRankStack[1] is nAtomNumber
    // INCHI✔️❌:            This should only split ranks of tautomeric groups */
    // INCHI✔️❌:         nNumCurrRanks = DifferentiateRanks4( pCG, num_at_tg, NeighList[TAUT_YES],
    // INCHI✔️❌:                                          nNumCurrRanks, pBCN->pRankStack[0], nTempRank /* temp array */,
    // INCHI✔️❌:                                          pBCN->pRankStack[1], (AT_RANK) num_atoms, &lCount ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_AUX_RANKING == 1 )
    // INCHI✔️❌:         /* refine no-H partition according to non-taut H distribution */
    // INCHI✔️❌:         memset( pAtomInvariantAux, 0, num_max * sizeof( pAtomInvariantAux[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pAtomInvariantAux[i].val[0] = nSymmRankNoTautH[i];
    // INCHI✔️❌:             pAtomInvariantAux[i].val[1] = numH[i]; /* additional differentiation */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         * djb-rwth: original badly written loop
    // INCHI✔️❌:         for (j = i; i < num_at_tg; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pAtomInvariantAux[i].val[0] = nRank[i];
    // INCHI✔️❌:         }
    // INCHI✔️❌:         */
    // INCHI✔️❌:         for (j = i; j < num_at_tg; j++) /* djb-rwth: corrected loop */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pAtomInvariantAux[j].val[0] = nRank[j];
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* initial ranks for t-group(s) */
    // INCHI✔️❌:         nNumCurrRanks = SetInitialRanks2( num_at_tg, pAtomInvariantAux, nRankAux, nAtomNumberAux, pCG );
    // INCHI✔️❌:         /* make equitable, call digraph procedure */
    // INCHI✔️❌:         nNumCurrRanks = DifferentiateRanks4( pCG, num_at_tg, NeighList[TAUT_YES],
    // INCHI✔️❌:                                          nNumCurrRanks, nRankAux, nTempRank /* temp array */,
    // INCHI✔️❌:                                          nAtomNumberAux, (AT_RANK) num_atoms, &lCount ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:         /* to accelerate do not call CanonGraph() to find really equivalent atoms */
    // INCHI✔️❌:         pCD[iBase].nAuxRank = nRankAux;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         ret = CanonGraph05( ic, pCG, num_atoms, num_at_tg, num_max, 1 /* digraph*/,
    // INCHI✔️❌:                             NeighList[TAUT_YES], (Partition *) pBCN->pRankStack,
    // INCHI✔️❌:                             nSymmRankBase, nCanonRankBase, nAtomNumberCanonBase,
    // INCHI✔️❌:                             pCD + iBase, pCC, &Ct_NoTautH, &Ct_Base, LargeMolecules );
    // INCHI✔️❌:         if (ret < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_error;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* tautomeric isotopic structure */
    // INCHI✔️❌:         /**************************************************************************************/
    // INCHI✔️❌:         /* Isotopic atoms and isotopic H atoms and isotopic tautomeric groups                 */
    // INCHI✔️❌:         /* get isotopic canonical numbering, connection table, and equivalence partition      */
    // INCHI✔️❌:         /**************************************************************************************/
    // INCHI✔️❌:         if ((s[iBase].num_isotopic_atoms && !s[iBase].bIgnoreIsotopic) ||
    // INCHI✔️❌:              (s[iBase].bHasIsotopicTautGroups && !bTautIgnoreIsotopic) ||
    // INCHI✔️❌:              (bUseIsoAuxBase[iBase] && !bTautIgnoreIsotopic)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:
    // INCHI✔️❌:             t_group_info->bIgnoreIsotopic = bTautIgnoreIsotopic;
    // INCHI✔️❌:
    // INCHI✔️❌:             nSymmRankBaseIso = (AT_RANK *) inchi_calloc( num_max, sizeof( nSymmRankBaseIso[0] ) );
    // INCHI✔️❌:             nCanonRankBaseIso = (AT_RANK *) inchi_calloc( num_max, sizeof( nCanonRankBaseIso[0] ) );
    // INCHI✔️❌:             nAtomNumberCanonBaseIso = (AT_NUMB *) inchi_calloc( num_max, sizeof( nAtomNumberCanonBaseIso[0] ) );
    // INCHI✔️❌:             if (bUseIsoAuxBase[iBase])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 maxlen_iso_exchg_atnos = num_max + 1;
    // INCHI✔️❌:                 iso_exchg_atnos = (S_CHAR  *) inchi_calloc( maxlen_iso_exchg_atnos, sizeof( iso_exchg_atnos[0] ) );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             maxlen_iso_sort_keyBase = num_max + 1; /* num_at_tg+1;*/
    // INCHI✔️❌:             iso_sort_keyBase = (AT_ISO_SORT_KEY *) inchi_calloc( maxlen_iso_sort_keyBase, sizeof( iso_sort_keyBase[0] ) );
    // INCHI✔️❌:             if (!nSymmRankBaseIso || !nCanonRankBaseIso || !nAtomNumberCanonBaseIso ||
    // INCHI✔️❌:                  !iso_sort_keyBase ||
    // INCHI✔️❌:                  (maxlen_iso_exchg_atnos && !iso_exchg_atnos)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto exit_error_alloc;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* atoms */
    // INCHI✔️❌:             num_iso_NoTautH = 0;
    // INCHI✔️❌:             num_iso_NoAuxBase = 0;
    // INCHI✔️❌:             if (iso_exchg_atnos)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 len_iso_exchg_atnos = num_at_tg;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at_base[i].endpoint || (iso_exchg_atnos && ( at_base[i].cFlags & AT_FLAG_ISO_H_POINT ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* tautomeric or may have exchangeable isotopic H */
    // INCHI✔️❌:                     iso_sort_key = make_iso_sort_key( at_base[i].iso_atw_diff, 0, 0, 0 );
    // INCHI✔️❌:                     if (iso_exchg_atnos)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         num_iso_NoAuxBase += !at_base[i].endpoint; /* these non-taut atom may exchange isotopic H as tautomeric atoms do */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* non-mobile H */
    // INCHI✔️❌:                     iso_sort_key = make_iso_sort_key( at_base[i].iso_atw_diff, at_base[i].num_iso_H[0], at_base[i].num_iso_H[1], at_base[i].num_iso_H[2] );
    // INCHI✔️❌:                     if (iso_exchg_atnos)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         iso_exchg_atnos[i] = 1; /* atom cannot have exchangable isotopic H atom(s) */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (iso_sort_key)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     num_iso_NoTautH++;
    // INCHI✔️❌:                     iso_sort_keyBase[i] = iso_sort_key;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     iso_sort_keyBase[i] = EMPTY_ISO_SORT_KEY;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* check marking and count of non-taut atoms that may exchange isotopic H -- debug only */
    // INCHI✔️❌:             if (iso_exchg_atnos)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (num_iso_NoAuxBase != t_group_info->nIsotopicEndpointAtomNumber[0])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ret = CT_ISOCOUNT_ERR;
    // INCHI✔️❌:                     goto exit_error;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 for (i = 1; i <= num_iso_NoAuxBase; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     j = t_group_info->nIsotopicEndpointAtomNumber[i];
    // INCHI✔️❌:                     if (at_base[j].endpoint || !( at_base[j].cFlags & AT_FLAG_ISO_H_POINT ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ret = CT_ISOCOUNT_ERR;
    // INCHI✔️❌:                         goto exit_error;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* t-groups */
    // INCHI✔️❌:             num_iso_Base = 0;
    // INCHI✔️❌:             if (iso_exchg_atnos)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (i = num_atoms; i < num_at_tg; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     iso_sort_keyBase[i] = EMPTY_ISO_SORT_KEY; /* new mode: do not provide info about isotopic tautomeric H */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (i = num_atoms; i < num_at_tg; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* should not happen anymore */
    // INCHI✔️❌:                     m = i - num_atoms;
    // INCHI✔️❌:                     if ((iso_sort_key = t_group_info->t_group[m].iWeight)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* old approach: each t-group has its own isotopic "weight" */
    // INCHI✔️❌:                         num_iso_Base++;
    // INCHI✔️❌:                         iso_sort_keyBase[i] = iso_sort_key;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         iso_sort_keyBase[i] = EMPTY_ISO_SORT_KEY;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (!num_iso_NoAuxBase && iso_exchg_atnos)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* all atoms that may exchange isotopic H are either tautomeric or not present */
    // INCHI✔️❌:                 inchi_free( iso_exchg_atnos );
    // INCHI✔️❌:                 iso_exchg_atnos = NULL;
    // INCHI✔️❌:                 len_iso_exchg_atnos = 0;
    // INCHI✔️❌:                 maxlen_iso_exchg_atnos = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (!num_iso_NoTautH && !num_iso_Base && iso_sort_keyBase)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* no isotopic atoms present */
    // INCHI✔️❌:                 inchi_free( iso_sort_keyBase );
    // INCHI✔️❌:                 iso_sort_keyBase = NULL;
    // INCHI✔️❌:                 maxlen_iso_sort_keyBase = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 len_iso_sort_keyBase = num_at_tg;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (!iso_exchg_atnos && !iso_sort_keyBase)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* no isotopic part at all or only tautomeric groups */
    // INCHI✔️❌:                 inchi_free( nSymmRankBaseIso );        nSymmRankBaseIso = NULL;
    // INCHI✔️❌:                 inchi_free( nCanonRankBaseIso );       nCanonRankBaseIso = NULL;
    // INCHI✔️❌:                 inchi_free( nAtomNumberCanonBaseIso ); nAtomNumberCanonBaseIso = NULL;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* proceed with tautomeric isotopic canonicalization */
    // INCHI✔️❌:                 /* pointers */
    // INCHI✔️❌:                 pCD[iBase].LinearCT = NULL;
    // INCHI✔️❌:                 pCD[iBase].NumH = numH; /* num_atoms non-tautomeric H; num_tg pairs of H and (-) in t-groups */
    // INCHI✔️❌:                 pCD[iBase].NumHfixed = NULL;
    // INCHI✔️❌:                 pCD[iBase].iso_sort_key = iso_sort_keyBase;
    // INCHI✔️❌:                 pCD[iBase].iso_exchg_atnos = iso_exchg_atnos;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:                 pCD[iBase].iso_sort_key_Hfixed = NULL;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* variables - unchanged */
    // INCHI✔️❌:                 pCD[iBase].ulTimeOutTime = pBCN->ulTimeOutTime;
    // INCHI✔️❌:                 pCD[iBase].nMaxLenLinearCT = s[iBase].nLenCT + 1;
    // INCHI✔️❌:                 pCD[iBase].maxlenNumH = maxlenNumH;
    // INCHI✔️❌:                 /* return values & input/output */
    // INCHI✔️❌:                 pCD[iBase].nLenLinearCT = s[iBase].nLenCT;
    // INCHI✔️❌:                 pCD[iBase].nLenCTAtOnly = s[iBase].nLenCTAtOnly;
    // INCHI✔️❌:                 pCD[iBase].lenNumH = lenNumH /* = k */;
    // INCHI✔️❌:                 pCD[iBase].lenNumHfixed = 0;
    // INCHI✔️❌:                 pCD[iBase].len_iso_sort_key = len_iso_sort_keyBase;
    // INCHI✔️❌:                 pCD[iBase].maxlen_iso_sort_key = maxlen_iso_sort_keyBase;
    // INCHI✔️❌:                 pCD[iBase].len_iso_exchg_atnos = len_iso_exchg_atnos;
    // INCHI✔️❌:                 pCD[iBase].maxlen_iso_exchg_atnos = maxlen_iso_exchg_atnos;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:                 pCD[iBase].len_iso_sort_key_Hfixed = 0;
    // INCHI✔️❌:                 pCD[iBase].maxlen_iso_sort_key_Hfixed = 0;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:                 pCD[iBase].nAuxRank = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (num_iso_NoTautH || num_iso_Base || num_iso_NoAuxBase)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* check whether we need actual canonicalization */
    // INCHI✔️❌:                     memset( nTempRank, 0, num_max * sizeof( nTempRank[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:                     for (i = 0; i < num_at_tg; i++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (nTempRank[nSymmRankBase[i] - 1] < i)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             nTempRank[nSymmRankBase[i] - 1] = i; /* greatest class representative */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     for (i = 0; i < num_at_tg; i++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (( iso_sort_keyBase ? ( iso_sort_keyBase[i] != iso_sort_keyBase[nTempRank[nSymmRankBase[i] - 1]] ) : 0 ) ||
    // INCHI✔️❌:                             ( iso_exchg_atnos ? ( iso_exchg_atnos[i] != iso_exchg_atnos[nTempRank[nSymmRankBase[i] - 1]] ) : 0 ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             pCD[iBase].nCanonFlags |= CANON_FLAG_ISO_TAUT_DIFF;
    // INCHI✔️❌:                             break; /* atoms so far found to be equivalent have different number of H; the canonicalization is needed */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     i = num_at_tg; /* should not happen */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* i = 0; *//* debug: force to call the canonicalization */
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (i < num_at_tg)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* we need canonicalization */
    // INCHI✔️❌:                     /* get aux canonical ranking of the structure with isotopic non-tautomeric H */
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_AUX_RANKING == 1 )
    // INCHI✔️❌:                 /* refine no-taut-H partition according to non-taut H + t-groups isotopic distribution */
    // INCHI✔️❌:                     memset( pAtomInvariantAux, 0, num_max * sizeof( pAtomInvariantAux[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:                     for (i = 0; i < num_at_tg; i++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         pAtomInvariantAux[i].val[0] = nSymmRankBase[i];
    // INCHI✔️❌:                         pAtomInvariantAux[i].iso_sort_key = iso_sort_keyBase ? iso_sort_keyBase[i] : 0; /* additional differentiation */
    // INCHI✔️❌:                         pAtomInvariantAux[i].iso_aux_key = iso_exchg_atnos ? iso_exchg_atnos[i] : 0;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /* initial ranks for non-taut H isotopic distribution */
    // INCHI✔️❌:                     nNumCurrRanks = SetInitialRanks2( num_at_tg, pAtomInvariantAux, nRankAux, nAtomNumberAux, pCG );
    // INCHI✔️❌:                     /* make equitable, not a digraph procedure */
    // INCHI✔️❌:                     nNumCurrRanks = DifferentiateRanks2( pCG, num_at_tg, NeighList[TAUT_YES],
    // INCHI✔️❌:                                                         nNumCurrRanks, nRankAux,
    // INCHI✔️❌:                                                         nTempRank, nAtomNumberAux, &lCount, 0 /* 0 means first use qsort */ ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                     /* to accelerate do not call CanonGraph() to find really equivalent atoms */
    // INCHI✔️❌:                     pCD[iBase].nAuxRank = nRankAux;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:                     ret = CanonGraph06( ic, pCG, num_atoms, num_at_tg, num_max, 1 /* digraph */,
    // INCHI✔️❌:                                         NeighList[TAUT_YES], (Partition *) pBCN->pRankStack,
    // INCHI✔️❌:                                         nSymmRankBaseIso, nCanonRankBaseIso, nAtomNumberCanonBaseIso,
    // INCHI✔️❌:                                         pCD + iBase, pCC, &Ct_Base, &Ct_BaseIso, LargeMolecules );
    // INCHI✔️❌:                     if (ret < 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         goto exit_error;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /* in case of a tautomeric structure the final results are in:
    // INCHI✔️❌:
    // INCHI✔️❌:                                nSymmRankBaseIso
    // INCHI✔️❌:                                nCanonRankBaseIso
    // INCHI✔️❌:                                nAtomNumberCanonBaseIso
    // INCHI✔️❌:                                Ct_BaseIso
    // INCHI✔️❌:                                iso_sort_keyBase (original isotopic atom & t-group positions)
    // INCHI✔️❌:                                Ct_BaseIso->iso_exchg_atnos: 0=>can exchange isotopic H, including tautomeric atoms
    // INCHI✔️❌:                                iso_exchg_atnos            : same, in order of t_group_info->nIsotopicEndpointAtomNumber[]
    // INCHI✔️❌:                     */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* copy the results of the previous (no taut H) canonicalization */
    // INCHI✔️❌:                     /* in this case numHNoTautH[] is not needed for the next canonicalization(s) */
    // INCHI✔️❌:                     if (( Ct_Temp = (ConTable *) inchi_calloc( 1, sizeof( *Ct_Temp ) ) ) &&
    // INCHI✔️❌:                          CTableCreate( Ct_Temp, num_atoms, pCD + iBase ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         CtFullCopy( Ct_Temp, Ct_Base );
    // INCHI✔️❌:                         /* since Ct_Base does not have Ct_Base->iso_sort_key we
    // INCHI✔️❌:                            have to fill out Ct_Temp->iso_sort_key separately */
    // INCHI✔️❌:                         if (iso_sort_keyBase)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             for (i = 0; i < num_at_tg; i++)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 Ct_Temp->iso_sort_key[nCanonRankBase[i] - 1] = iso_sort_keyBase[i];
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             Ct_Temp->len_iso_sort_key = num_at_tg;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             Ct_Temp->len_iso_sort_key = 0;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (iso_exchg_atnos)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 Ct_Temp->iso_exchg_atnos[nCanonRankBase[i] - 1] = iso_exchg_atnos[i];
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             Ct_Temp->len_iso_exchg_atnos = num_at_tg;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             Ct_Temp->len_iso_exchg_atnos = 0;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         goto exit_error_alloc;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     Ct_BaseIso = Ct_Temp;
    // INCHI✔️❌:                     Ct_Temp = NULL;
    // INCHI✔️❌:                     memcpy(nSymmRankBaseIso, nSymmRankBase, num_at_tg * sizeof(nSymmRankBaseIso[0]));
    // INCHI✔️❌:                     memcpy(nCanonRankBaseIso, nCanonRankBase, num_at_tg * sizeof(nCanonRankBaseIso[0]));
    // INCHI✔️❌:                     memcpy(nAtomNumberCanonBaseIso, nAtomNumberCanonBase, num_at_tg * sizeof(nAtomNumberCanonBaseIso[0]));
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* in case of non-tautomeric component this is the final result */
    // INCHI✔️❌:                 /* i = CtFullCompare( Ct_BaseIso, Ct_Temp, num_at_tg, 0, 0 );*/
    // INCHI✔️❌:
    // INCHI✔️❌:                 t_group_info->bIgnoreIsotopic = 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /**********************************************************************************/
    // INCHI✔️❌:     /* get "fixed H" canonical numbering, connection table, and equivalence partition */
    // INCHI✔️❌:     /**********************************************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:     if (bReqTaut && bReqNonTaut)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         maxlenNumHfixed = num_atoms + 1;
    // INCHI✔️❌:         nSymmRankFixH = (AT_RANK *) inchi_calloc( num_max, sizeof( nSymmRankFixH[0] ) );
    // INCHI✔️❌:         nCanonRankFixH = (AT_RANK *) inchi_calloc( num_max, sizeof( nCanonRankFixH[0] ) );
    // INCHI✔️❌:         nAtomNumberCanonFixH = (AT_NUMB *) inchi_calloc( num_max, sizeof( nAtomNumberCanonFixH[0] ) );
    // INCHI✔️❌:         NumHfixed = (NUM_H *) inchi_calloc( maxlenNumHfixed, sizeof( NumHfixed[0] ) );
    // INCHI✔️❌:         if (!NumHfixed || !nSymmRankFixH || !nCanonRankFixH || !nAtomNumberCanonFixH)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_error_alloc;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* fixed and non-tautomeric H different in taut and non-taut structures */
    // INCHI✔️❌:
    // INCHI✔️❌:             if (at_base[i].endpoint)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 NumHfixed[i] = at_other[i].num_H ? at_other[i].num_H + BASE_H_NUMBER : EMPTY_H_NUMBER;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else if (at_other[i].num_H != at_base[i].num_H)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 NumHfixed[i] = (NUM_H) at_other[i].num_H - (NUM_H) at_base[i].num_H + BASE_H_NUMBER;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 NumHfixed[i] = EMPTY_H_NUMBER;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* pointers */
    // INCHI✔️❌:         pCD[iOther].LinearCT = NULL; /* LinearCT; */
    // INCHI✔️❌:         pCD[iOther].NumH = numHNoTautH;
    // INCHI✔️❌:         pCD[iOther].NumHfixed = NumHfixed;/* variables - unchanged */
    // INCHI✔️❌:         pCD[iOther].iso_sort_key = NULL;
    // INCHI✔️❌:         pCD[iOther].iso_exchg_atnos = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:         pCD[iOther].iso_sort_key_Hfixed = NULL;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         pCD[iOther].ulTimeOutTime = pBCN->ulTimeOutTime;
    // INCHI✔️❌:         pCD[iOther].nMaxLenLinearCT = s[iOther].nLenCTAtOnly + 1;
    // INCHI✔️❌:         pCD[iOther].maxlenNumH = maxlenNumHNoTautH;
    // INCHI✔️❌:         pCD[iOther].maxlenNumHfixed = maxlenNumHfixed;
    // INCHI✔️❌:         /* return values & input/output */
    // INCHI✔️❌:         pCD[iOther].nLenLinearCT = s[iOther].nLenCTAtOnly;
    // INCHI✔️❌:         pCD[iOther].nLenCTAtOnly = s[iOther].nLenCTAtOnly;
    // INCHI✔️❌:         pCD[iOther].lenNumH = num_atoms; /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         pCD[iOther].lenNumHfixed = num_atoms;
    // INCHI✔️❌:         pCD[iOther].len_iso_sort_key = 0;
    // INCHI✔️❌:         pCD[iOther].maxlen_iso_sort_key = 0;
    // INCHI✔️❌:         pCD[iOther].len_iso_exchg_atnos = 0;
    // INCHI✔️❌:         pCD[iOther].maxlen_iso_exchg_atnos = 0;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:         pCD[iOther].len_iso_sort_key_Hfixed = 0;
    // INCHI✔️❌:         pCD[iOther].maxlen_iso_sort_key_Hfixed = 0;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         pCD[iOther].nAuxRank = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_AUX_RANKING == 1 )
    // INCHI✔️❌:         if (!nRankAux)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nRankAux = (AT_RANK *) inchi_calloc( num_max, sizeof( nRankAux[0] ) );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (!nAtomNumberAux)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nAtomNumberAux = (AT_NUMB *) inchi_calloc( num_max, sizeof( nAtomNumberAux[0] ) );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (!pAtomInvariantAux)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pAtomInvariantAux = (ATOM_INVARIANT2 *) inchi_malloc( num_max * sizeof( pAtomInvariantAux[0] ) );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (!nRankAux || !nAtomNumberAux || !pAtomInvariantAux)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_error_alloc;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* refine no-H partition according to non-taut H distribution */
    // INCHI✔️❌:         memset( pAtomInvariantAux, 0, num_max * sizeof( pAtomInvariantAux[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pAtomInvariantAux[i].val[0] = nSymmRankBase[i];
    // INCHI✔️❌:             pAtomInvariantAux[i].val[1] = NumHfixed[i]; /* additional differentiation */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* initial ranks for t-group(s) */
    // INCHI✔️❌:         nNumCurrRanks = SetInitialRanks2( num_atoms, pAtomInvariantAux, nRankAux, nAtomNumberAux, pCG );
    // INCHI✔️❌:
    // INCHI✔️❌:         /* make equitable, digraph procedure */
    // INCHI✔️❌:         nNumCurrRanks = DifferentiateRanks2( pCG, num_atoms, NeighList[TAUT_NON],
    // INCHI✔️❌:                                             nNumCurrRanks, nRankAux,
    // INCHI✔️❌:                                             nTempRank, nAtomNumberAux, &lCount, 0 /* 0 means use qsort */ ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:         /* to accelerate do not call CanonGraph() to find really equivalent atoms */
    // INCHI✔️❌:         pCD[iOther].nAuxRank = nRankAux;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         ret = CanonGraph07( ic, pCG, num_atoms, num_atoms, num_max, 0,
    // INCHI✔️❌:                             NeighList[TAUT_NON], (Partition *) pBCN->pRankStack,
    // INCHI✔️❌:                             nSymmRankFixH, nCanonRankFixH, nAtomNumberCanonFixH,
    // INCHI✔️❌:                             pCD + iOther, pCC, &Ct_NoTautH, &Ct_FixH, LargeMolecules );
    // INCHI✔️❌:         if (ret < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_error;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /*******************************************************************************************/
    // INCHI✔️❌:         /* get "fixed H" isotopic canonical numbering, connection table, and equivalence partition */
    // INCHI✔️❌:         /*******************************************************************************************/
    // INCHI✔️❌:         iflag = (s[iBase].num_isotopic_atoms && !s[iBase].bIgnoreIsotopic) ||
    // INCHI✔️❌:             (s[iBase].bHasIsotopicTautGroups && !bTautIgnoreIsotopic); /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         if (bFixIsoFixedH) /* #if ( FIX_ISO_FIXEDH_BUG == 1 )  */
    // INCHI✔️❌:              /* fix bug when iso H was removed as a proton and fixed-H isotopic layer is missing -  2008-09-24 DT*/
    // INCHI✔️❌:         {
    // INCHI✔️❌:             iflag = iflag || (s[iOther].num_isotopic_atoms && !s[iOther].bIgnoreIsotopic); /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (iflag)
    // INCHI✔️❌:         {
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:             maxlen_iso_sort_key_Hfixed =
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             maxlen_iso_sort_key_NoTautH = num_atoms + 1;
    // INCHI✔️❌:             nSymmRankFixHIso = (AT_RANK *) inchi_calloc( num_max, sizeof( nSymmRankFixHIso[0] ) );
    // INCHI✔️❌:             nCanonRankFixHIso = (AT_RANK *) inchi_calloc( num_max, sizeof( nCanonRankFixHIso[0] ) );
    // INCHI✔️❌:             nAtomNumberCanonFixHIso = (AT_NUMB *) inchi_calloc( num_max, sizeof( nAtomNumberCanonFixHIso[0] ) );
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:             iso_sort_key_Hfixed = (AT_ISO_SORT_KEY *) inchi_calloc( maxlen_iso_sort_key_Hfixed, sizeof( iso_sort_key_Hfixed[0] ) );
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:             iso_sort_key_NoTautH = (AT_ISO_SORT_KEY *) inchi_calloc( maxlen_iso_sort_key_NoTautH, sizeof( iso_sort_key_NoTautH[0] ) );
    // INCHI✔️❌:
    // INCHI✔️❌:             if (!nSymmRankFixHIso || !nCanonRankFixHIso || !nAtomNumberCanonFixHIso ||
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:                  !iso_sort_key_Hfixed ||
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:                  !iso_sort_key_NoTautH)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto exit_error_alloc;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:             /* fill out isotopic non-tautomeric keys */
    // INCHI✔️❌:             for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at_base[i].endpoint)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     iso_sort_key = make_iso_sort_key( at_base[i].iso_atw_diff, 0, 0, 0 );
    // INCHI✔️❌:                     iso_sort_key2 = make_iso_sort_key( 0, at_base[i].num_iso_H[0], at_base[i].num_iso_H[1], at_base[i].num_iso_H[2] );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     iso_sort_key = make_iso_sort_key( at_base[i].iso_atw_diff, at_base[i].num_iso_H[0], at_base[i].num_iso_H[1], at_base[i].num_iso_H[2] );
    // INCHI✔️❌:                     iso_sort_key2 = 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (iso_sort_key)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     iso_sort_key_NoTautH[i] = iso_sort_key;
    // INCHI✔️❌:                     num_iso_NoTautH++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     iso_sort_key_NoTautH[i] = EMPTY_ISO_SORT_KEY;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (iso_sort_key2)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     num_iso_Hfixed++;
    // INCHI✔️❌:                     iso_sort_key_Hfixed[i] = iso_sort_key2;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     iso_sort_key_Hfixed[i] = EMPTY_ISO_SORT_KEY;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌: #else
    // INCHI✔️❌:             /* fill out isotopic non-tautomeric keys */
    // INCHI✔️❌:             for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (bFixIsoFixedH) /* #if ( FIX_ISO_FIXEDH_BUG == 1 )  */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* fix bug when iso H was removed as a proton and fixed-H isotopic layer is missing -  2008-09-24 DT*/
    // INCHI✔️❌:                     if (at_other)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         iso_sort_key = make_iso_sort_key( at_other[i].iso_atw_diff, at_other[i].num_iso_H[0], at_other[i].num_iso_H[1], at_other[i].num_iso_H[2] );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         iso_sort_key = make_iso_sort_key( at_base[i].iso_atw_diff, at_base[i].num_iso_H[0], at_base[i].num_iso_H[1], at_base[i].num_iso_H[2] );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                     iso_sort_key = make_iso_sort_key( at_base[i].iso_atw_diff, at_base[i].num_iso_H[0], at_base[i].num_iso_H[1], at_base[i].num_iso_H[2] );
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (iso_sort_key)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     iso_sort_key_NoTautH[i] = iso_sort_key;
    // INCHI✔️❌:                     num_iso_NoTautH++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     iso_sort_key_NoTautH[i] = EMPTY_ISO_SORT_KEY;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:             /* pointers */
    // INCHI✔️❌:             pCD[iOther].LinearCT = NULL; /* LinearCT; */
    // INCHI✔️❌:             pCD[iOther].NumH = numHNoTautH;
    // INCHI✔️❌:             pCD[iOther].NumHfixed = NumHfixed;/* variables - unchanged */
    // INCHI✔️❌:             pCD[iOther].iso_sort_key = iso_sort_key_NoTautH;
    // INCHI✔️❌:             pCD[iOther].iso_exchg_atnos = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:             pCD[iOther].iso_sort_key_Hfixed = iso_sort_key_Hfixed;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:             pCD[iOther].ulTimeOutTime = pBCN->ulTimeOutTime;
    // INCHI✔️❌:             pCD[iOther].nMaxLenLinearCT = s[iOther].nLenCTAtOnly + 1;
    // INCHI✔️❌:             pCD[iOther].maxlenNumH = maxlenNumHNoTautH;
    // INCHI✔️❌:             pCD[iOther].maxlenNumHfixed = maxlenNumHfixed;
    // INCHI✔️❌:             /* return values & input/output */
    // INCHI✔️❌:             pCD[iOther].nLenLinearCT = s[iOther].nLenCTAtOnly;
    // INCHI✔️❌:             pCD[iOther].nLenCTAtOnly = s[iOther].nLenCTAtOnly;
    // INCHI✔️❌:             pCD[iOther].lenNumH = num_atoms; /* djb-rwth: removing redundant code */
    // INCHI✔️❌:             pCD[iOther].lenNumHfixed = num_atoms;
    // INCHI✔️❌:             pCD[iOther].len_iso_sort_key = len_iso_sort_key_NoTautH = num_atoms;
    // INCHI✔️❌:             pCD[iOther].maxlen_iso_sort_key = maxlen_iso_sort_key_NoTautH;
    // INCHI✔️❌:             pCD[iOther].len_iso_exchg_atnos = 0;
    // INCHI✔️❌:             pCD[iOther].maxlen_iso_exchg_atnos = 0;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:             pCD[iOther].len_iso_sort_key_Hfixed = len_iso_sort_key_Hfixed = num_atoms;
    // INCHI✔️❌:             pCD[iOther].maxlen_iso_sort_key_Hfixed = maxlen_iso_sort_key_Hfixed;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:             pCD[iOther].nAuxRank = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:             if (num_iso_Hfixed || num_iso_NoTautH)
    // INCHI✔️❌: #else
    // INCHI✔️❌:             if (num_iso_NoTautH)
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* check whether we need NoTautH cononicalization */
    // INCHI✔️❌:                 memset( nTempRank, 0, num_max * sizeof( nTempRank[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:                 for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (nTempRank[nSymmRankFixH[i] - 1] < i)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         nTempRank[nSymmRankFixH[i] - 1] = i; /* greatest class representative */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:                     if (iso_sort_key_Hfixed[i] != iso_sort_key_Hfixed[nTempRank[nSymmRankFixH[i] - 1]])
    // INCHI✔️❌:                         break;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     if (iso_sort_key_NoTautH[i] != iso_sort_key_NoTautH[nTempRank[nSymmRankFixH[i] - 1]])
    // INCHI✔️❌:                         break; /* atoms so far found to be equivalent have different isotopic shifts; the canonicalization is needed */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 i = num_atoms; /* should not happen */
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* i = 0; *//* debug: force to call the canonicalization */
    // INCHI✔️❌:
    // INCHI✔️❌:             if (i < num_atoms)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pCD[iOther].nCanonFlags |= CANON_FLAG_ISO_FIXED_H_DIFF;
    // INCHI✔️❌:                 /* we need canonicalization */
    // INCHI✔️❌:                 /* get aux canonical ranking of the structure with isotopic non-tautomeric H */
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_AUX_RANKING == 1 )
    // INCHI✔️❌:                 /* refine fixed-taut-H partition according to the isotopic distribution */
    // INCHI✔️❌:                 memset( pAtomInvariantAux, 0, num_max * sizeof( pAtomInvariantAux[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:                 for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     pAtomInvariantAux[i].val[0] = nSymmRankFixH[i];
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:                     iso_sort_key = 0;
    // INCHI✔️❌:                     if (iso_sort_key_NoTautH[i] != EMPTY_ISO_SORT_KEY)
    // INCHI✔️❌:                         iso_sort_key |= iso_sort_key_NoTautH[i];
    // INCHI✔️❌:                     if (iso_sort_key_Hfixed[i] != EMPTY_ISO_SORT_KEY)
    // INCHI✔️❌:                         iso_sort_key |= iso_sort_key_Hfixed[i];
    // INCHI✔️❌:                     if (!iso_sort_key)
    // INCHI✔️❌:                         iso_sort_key = EMPTY_ISO_SORT_KEY;
    // INCHI✔️❌: #else
    // INCHI✔️❌:                     iso_sort_key = iso_sort_key_NoTautH[i];
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:                     pAtomInvariantAux[i].iso_sort_key = iso_sort_key; /* additional differentiation */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* initial ranks for non-taut H isotopic distribution */
    // INCHI✔️❌:                 nNumCurrRanks = SetInitialRanks2( num_atoms, pAtomInvariantAux, nRankAux, nAtomNumberAux, pCG );
    // INCHI✔️❌:                 /* make equitable, digraph procedure */
    // INCHI✔️❌:                 nNumCurrRanks = DifferentiateRanks2( pCG, num_atoms, NeighList[TAUT_NON],
    // INCHI✔️❌:                                                     nNumCurrRanks, nRankAux,
    // INCHI✔️❌:                                                     nTempRank, nAtomNumberAux, &lCount, 0 /* 0 means use qsort */ ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* to accelerate do not call CanonGraph() to find really equivalent atoms */
    // INCHI✔️❌:                 pCD[iOther].nAuxRank = nRankAux;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:                 ret = CanonGraph08( ic, pCG, num_atoms, num_atoms, num_max, 1 /* digraph?? was 0 */,
    // INCHI✔️❌:                                     NeighList[TAUT_NON], (Partition *) pBCN->pRankStack,
    // INCHI✔️❌:                                     nSymmRankFixHIso, nCanonRankFixHIso, nAtomNumberCanonFixHIso,
    // INCHI✔️❌:                                     pCD + iOther, pCC, &Ct_FixH, &Ct_FixHIso, LargeMolecules );
    // INCHI✔️❌:                 if (ret < 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     goto exit_error;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* in case of non-tautomeric structure the final results are in:
    // INCHI✔️❌:
    // INCHI✔️❌:                            nSymmRankFixHIso
    // INCHI✔️❌:                            nCanonRankFixHIso
    // INCHI✔️❌:                            nAtomNumberCanonFixHIso
    // INCHI✔️❌:                            Ct_FixHIso
    // INCHI✔️❌:                            iso_sort_keyBase     ([0..num_atoms] original isotopic atom positions)
    // INCHI✔️❌:                            iso_sort_key_Hfixed  (original fixed tautomeric H distribution)
    // INCHI✔️❌:                 */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* copy the results of the previous (no taut H) canonicalization */
    // INCHI✔️❌:                 /* in this case numHNoTautH[] is not needed for the next canonicalization(s) */
    // INCHI✔️❌:                 if (( Ct_Temp = (ConTable *) inchi_calloc( 1, sizeof( *Ct_Temp ) ) ) &&
    // INCHI✔️❌:                      CTableCreate( Ct_Temp, num_atoms, pCD + iOther ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     CtFullCopy( Ct_Temp, Ct_FixH );
    // INCHI✔️❌:                     /* since Ct_FixH does not have Ct_FixH->iso_sort_key and Ct_FixH->iso_sort_key_Hfixed we
    // INCHI✔️❌:                        have to fill out Ct_Temp->iso_sort_key and Ct_Temp->iso_sort_key_Hfixed separately */
    // INCHI✔️❌:                     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         Ct_Temp->iso_sort_key[nCanonRankFixH[i] - 1] = iso_sort_key_NoTautH[i];
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:                         Ct_Temp->iso_sort_key_Hfixed[nCanonRankFixH[i] - 1] = iso_sort_key_Hfixed[i];
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     Ct_Temp->len_iso_sort_key = num_atoms;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:                     Ct_Temp->len_iso_sort_key_Hfixed = num_atoms;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:                     /*Ct_Temp->lenNumH = num_atoms;*/
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     goto exit_error_alloc;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 Ct_FixHIso = Ct_Temp;
    // INCHI✔️❌:                 Ct_Temp = NULL;
    // INCHI✔️❌:                 memcpy(nSymmRankFixHIso, nSymmRankFixH, num_atoms * sizeof(nSymmRankFixHIso[0]));
    // INCHI✔️❌:                 memcpy(nCanonRankFixHIso, nCanonRankFixH, num_atoms * sizeof(nCanonRankFixHIso[0]));
    // INCHI✔️❌:                 memcpy(nAtomNumberCanonFixHIso, nAtomNumberCanonFixH, num_atoms * sizeof(nAtomNumberCanonFixHIso[0]));
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* in case of non-tautomeric component this is the final result */
    // INCHI✔️❌:             /* i = CtFullCompare( Ct_NoTautHIso, Ct_Temp, num_atoms, 0, 0 );*/
    // INCHI✔️❌:         }
    // INCHI✔️❌:     } /* "fixed H" canonical numbering */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* consistency check: compare canonical connection tables, H-atoms, isotopic H & taut groups */
    // INCHI✔️❌:     ret = 0;
    // INCHI✔️❌:     ret |= ( Ct_NoH->lenCt != Ct_NoTautH->lenCt ) || memcmp( Ct_NoH->Ctbl, Ct_NoTautH->Ctbl, Ct_NoH->lenCt * sizeof( Ct_NoH->Ctbl[0] ) );
    // INCHI✔️❌:     if (bReqTaut)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Ct_FixH)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret |= ( Ct_NoTautH->lenCt != Ct_FixH->lenCt ) || memcmp( Ct_NoTautH->Ctbl, Ct_FixH->Ctbl, Ct_NoTautH->lenCt * sizeof( Ct_NoTautH->Ctbl[0] ) );
    // INCHI✔️❌:             ret |= ( Ct_NoTautH->lenNumH != Ct_FixH->lenNumH ) || memcmp( Ct_NoTautH->NumH, Ct_FixH->NumH, Ct_NoTautH->lenNumH * sizeof( Ct_Base->NumH[0] ) );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         ret |= ( Ct_NoTautH->lenCt > Ct_Base->lenCt ) || memcmp( Ct_NoTautH->Ctbl, Ct_Base->Ctbl, Ct_NoTautH->lenCt * sizeof( Ct_NoTautH->Ctbl[0] ) );
    // INCHI✔️❌:         ret |= ( Ct_NoTautH->lenNumH > Ct_Base->lenNumH ) || memcmp( Ct_NoTautH->NumH, Ct_Base->NumH, Ct_NoTautH->lenNumH * sizeof( Ct_Base->NumH[0] ) );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* isotopic canonicalization */
    // INCHI✔️❌:     if (Ct_NoTautHIso)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret |= ( Ct_NoH->lenCt != Ct_NoTautHIso->lenCt ) || memcmp( Ct_NoH->Ctbl, Ct_NoTautHIso->Ctbl, Ct_NoH->lenCt * sizeof( Ct_NoH->Ctbl[0] ) );
    // INCHI✔️❌:         ret |= ( Ct_NoTautH->lenNumH != Ct_NoTautHIso->lenNumH ) || memcmp( Ct_NoTautH->NumH, Ct_NoTautHIso->NumH, Ct_NoTautH->lenNumH * sizeof( Ct_Base->NumH[0] ) );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else if (Ct_BaseIso)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret |= ( Ct_BaseIso->lenCt != Ct_Base->lenCt ) || memcmp( Ct_BaseIso->Ctbl, Ct_Base->Ctbl, Ct_BaseIso->lenCt * sizeof( Ct_BaseIso->Ctbl[0] ) );
    // INCHI✔️❌:         ret |= ( Ct_BaseIso->lenNumH != Ct_Base->lenNumH ) || memcmp( Ct_BaseIso->NumH, Ct_Base->NumH, Ct_BaseIso->lenNumH * sizeof( Ct_BaseIso->NumH[0] ) );
    // INCHI✔️❌:         if (Ct_FixHIso)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret |= ( Ct_FixHIso->lenCt > Ct_BaseIso->lenCt ) || memcmp( Ct_FixHIso->Ctbl, Ct_BaseIso->Ctbl, Ct_FixHIso->lenCt * sizeof( Ct_FixHIso->Ctbl[0] ) );
    // INCHI✔️❌:             ret |= ( Ct_FixHIso->lenNumH > Ct_BaseIso->lenNumH ) || memcmp( Ct_FixHIso->NumH, Ct_BaseIso->NumH, Ct_FixHIso->lenNumH * sizeof( Ct_BaseIso->NumH[0] ) );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ret)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_error;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (bReqTaut)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* restore save "process isotopic" mark; temporarily set it to NO */
    // INCHI✔️❌:         t_group_info->bIgnoreIsotopic = bTautIgnoreIsotopic;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* output the canonicalization results */
    // INCHI✔️❌:     pBCN->num_max = num_max;
    // INCHI✔️❌:     pBCN->num_at_tg = num_at_tg;
    // INCHI✔️❌:     pBCN->num_atoms = num_atoms;
    // INCHI✔️❌:
    // INCHI✔️❌:     pBCN->ftcn[TAUT_NON].NeighList = NeighList[TAUT_NON]; NeighList[TAUT_NON] = NULL;
    // INCHI✔️❌:     pBCN->ftcn[TAUT_YES].NeighList = NeighList[TAUT_YES]; NeighList[TAUT_YES] = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (bReqTaut)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* tautomeric results */
    // INCHI✔️❌:
    // INCHI✔️❌:         /* base tautomeric structure, iBase = TAUT_YES */
    // INCHI✔️❌:
    // INCHI✔️❌:         pBCN->ftcn[TAUT_YES].num_at_tg = num_at_tg;
    // INCHI✔️❌:         pBCN->ftcn[TAUT_YES].num_atoms = num_atoms;
    // INCHI✔️❌:
    // INCHI✔️❌:         pBCN->ftcn[TAUT_YES].LinearCt = Ct_Base->Ctbl;            Ct_Base->Ctbl = NULL;
    // INCHI✔️❌:         pBCN->ftcn[TAUT_YES].nLenLinearCtAtOnly = s[iBase].nLenCTAtOnly;
    // INCHI✔️❌:         pBCN->ftcn[TAUT_YES].nMaxLenLinearCt = s[iBase].nLenCT + 1;
    // INCHI✔️❌:         pBCN->ftcn[TAUT_YES].nLenLinearCt = s[iBase].nLenCT;
    // INCHI✔️❌:
    // INCHI✔️❌:         pBCN->ftcn[TAUT_YES].PartitionCt.Rank = nCanonRankBase;        nCanonRankBase = NULL;
    // INCHI✔️❌:         pBCN->ftcn[TAUT_YES].PartitionCt.AtNumber = nAtomNumberCanonBase;  nAtomNumberCanonBase = NULL;
    // INCHI✔️❌:         pBCN->ftcn[TAUT_YES].nSymmRankCt = nSymmRankBase;         nSymmRankBase = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:         pBCN->ftcn[TAUT_YES].nNumHOrig = numH;                 numH = NULL;
    // INCHI✔️❌:         pBCN->ftcn[TAUT_YES].nNumH = Ct_Base->NumH;        Ct_Base->NumH = NULL;
    // INCHI✔️❌:         pBCN->ftcn[TAUT_YES].nLenNumH = inchi_min( maxlenNumH, Ct_Base->maxlenNumH );
    // INCHI✔️❌:
    // INCHI✔️❌:         /* fixed H structure: exists only if the structure is tautomeric */
    // INCHI✔️❌:         pBCN->ftcn[TAUT_YES].nNumHOrigFixH = NULL;
    // INCHI✔️❌:         pBCN->ftcn[TAUT_YES].nNumHFixH = NULL;
    // INCHI✔️❌:         pBCN->ftcn[TAUT_YES].nLenNumHFixH = 0;
    // INCHI✔️❌:         pBCN->ftcn[TAUT_YES].nCanonFlags |= pCD[iBase].nCanonFlags;
    // INCHI✔️❌:
    // INCHI✔️❌:         CleanNumH( pBCN->ftcn[TAUT_YES].nNumHOrig, pBCN->ftcn[TAUT_YES].nLenNumH );
    // INCHI✔️❌:         CleanNumH( pBCN->ftcn[TAUT_YES].nNumH, pBCN->ftcn[TAUT_YES].nLenNumH );
    // INCHI✔️❌:         CleanCt( pBCN->ftcn[TAUT_YES].LinearCt, pBCN->ftcn[TAUT_YES].nLenLinearCt );
    // INCHI✔️❌:
    // INCHI✔️❌:         /* isotopic canonicalization */
    // INCHI✔️❌:         if (Ct_BaseIso)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pBCN->ftcn[TAUT_YES].PartitionCtIso.Rank = nCanonRankBaseIso;           nCanonRankBaseIso = NULL;
    // INCHI✔️❌:             pBCN->ftcn[TAUT_YES].PartitionCtIso.AtNumber = nAtomNumberCanonBaseIso;     nAtomNumberCanonBaseIso = NULL;
    // INCHI✔️❌:             pBCN->ftcn[TAUT_YES].nSymmRankCtIso = nSymmRankBaseIso;            nSymmRankBaseIso = NULL;
    // INCHI✔️❌:             pBCN->ftcn[TAUT_YES].iso_sort_keys = Ct_BaseIso->iso_sort_key;    Ct_BaseIso->iso_sort_key = NULL;
    // INCHI✔️❌:             pBCN->ftcn[TAUT_YES].iso_sort_keysOrig = iso_sort_keyBase;            iso_sort_keyBase = NULL;
    // INCHI✔️❌:             pBCN->ftcn[TAUT_YES].len_iso_sort_keys = len_iso_sort_keyBase;
    // INCHI✔️❌:             pBCN->ftcn[TAUT_YES].iso_exchg_atnos = Ct_BaseIso->iso_exchg_atnos; Ct_BaseIso->iso_exchg_atnos = NULL;
    // INCHI✔️❌:             pBCN->ftcn[TAUT_YES].iso_exchg_atnosOrig = iso_exchg_atnos;             iso_exchg_atnos = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:             CleanIsoSortKeys( pBCN->ftcn[TAUT_YES].iso_sort_keys, pBCN->ftcn[TAUT_YES].len_iso_sort_keys );
    // INCHI✔️❌:             CleanIsoSortKeys( pBCN->ftcn[TAUT_YES].iso_sort_keysOrig, pBCN->ftcn[TAUT_YES].len_iso_sort_keys );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     } /* tautomeric results */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (bReqNonTaut)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* non-tautomeric results */
    // INCHI✔️❌:
    // INCHI✔️❌:         /* TAUT_NON if tautomeric + non-tautomeric or special non-taut request
    // INCHI✔️❌:            TAUT_YES if the structure happened to be non-tautomeric while user requested tautomeric processing
    // INCHI✔️❌:            In both cases the correct index is iOther. TAUT_NON replaced with iOther 4-2-2004 */
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!bReqTaut)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* rearrange the results for a non-tautomeric structure */
    // INCHI✔️❌:             nSymmRankFixH = nSymmRankNoTautH;           nSymmRankNoTautH = NULL;
    // INCHI✔️❌:             nCanonRankFixH = nCanonRankNoTautH;          nCanonRankNoTautH = NULL;
    // INCHI✔️❌:             nAtomNumberCanonFixH = nAtomNumberCanonNoTautH;    nAtomNumberCanonNoTautH = NULL;
    // INCHI✔️❌:             Ct_FixH = Ct_NoTautH;                 Ct_NoTautH = NULL;
    // INCHI✔️❌:             /* isotopic canonicalization */
    // INCHI✔️❌:             nSymmRankFixHIso = nSymmRankNoTautHIso;        nSymmRankNoTautHIso = NULL;
    // INCHI✔️❌:             nCanonRankFixHIso = nCanonRankNoTautHIso;       nCanonRankNoTautHIso = NULL;
    // INCHI✔️❌:             nAtomNumberCanonFixHIso = nAtomNumberCanonNoTautHIso; nAtomNumberCanonNoTautHIso = NULL;
    // INCHI✔️❌:             Ct_FixHIso = Ct_NoTautHIso;              Ct_NoTautHIso = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:             if (iOther == TAUT_YES && pBCN->ftcn[TAUT_NON].NeighList && !pBCN->ftcn[TAUT_YES].NeighList)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* here only non-taut results go to pBCN->ftcn[TAUT_YES]
    // INCHI✔️❌:                    Since non-taut NeighList is always in pBCN->ftcn[TAUT_NON].NeighList, move it to
    // INCHI✔️❌:                    pBCN->ftcn[TAUT_YES].NeighList. 2004-04-02.
    // INCHI✔️❌:                 */
    // INCHI✔️❌:                 pBCN->ftcn[TAUT_YES].NeighList = pBCN->ftcn[TAUT_NON].NeighList;
    // INCHI✔️❌:                 pBCN->ftcn[TAUT_NON].NeighList = NULL;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         pBCN->ftcn[iOther].num_at_tg = num_atoms;
    // INCHI✔️❌:         pBCN->ftcn[iOther].num_atoms = num_atoms;
    // INCHI✔️❌:
    // INCHI✔️❌:         pBCN->ftcn[iOther].LinearCt = Ct_FixH->Ctbl;
    // INCHI✔️❌:         Ct_FixH->Ctbl = NULL;
    // INCHI✔️❌:         pBCN->ftcn[iOther].nLenLinearCtAtOnly = s[iOther].nLenCTAtOnly;
    // INCHI✔️❌:         pBCN->ftcn[iOther].nMaxLenLinearCt = s[iOther].nLenCTAtOnly + 1;
    // INCHI✔️❌:         pBCN->ftcn[iOther].nLenLinearCt = s[iOther].nLenCTAtOnly;
    // INCHI✔️❌:
    // INCHI✔️❌:         pBCN->ftcn[iOther].PartitionCt.Rank = nCanonRankFixH;
    // INCHI✔️❌:         nCanonRankFixH = NULL;
    // INCHI✔️❌:         pBCN->ftcn[iOther].PartitionCt.AtNumber = nAtomNumberCanonFixH;
    // INCHI✔️❌:         nAtomNumberCanonFixH = NULL;
    // INCHI✔️❌:         pBCN->ftcn[iOther].nSymmRankCt = nSymmRankFixH;
    // INCHI✔️❌:         nSymmRankFixH = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:         pBCN->ftcn[iOther].nNumHOrig = numHNoTautH;
    // INCHI✔️❌:         numHNoTautH = NULL;
    // INCHI✔️❌:         pBCN->ftcn[iOther].nNumH = Ct_FixH->NumH;
    // INCHI✔️❌:         Ct_FixH->NumH = NULL;
    // INCHI✔️❌:         pBCN->ftcn[iOther].nLenNumH = inchi_min( maxlenNumHNoTautH, Ct_FixH->maxlenNumH );
    // INCHI✔️❌:
    // INCHI✔️❌:         /* fixed H structure: exists only if the structure is tautomeric */
    // INCHI✔️❌:         pBCN->ftcn[iOther].nNumHOrigFixH = NumHfixed;
    // INCHI✔️❌:         NumHfixed = NULL;
    // INCHI✔️❌:         pBCN->ftcn[iOther].nNumHFixH = Ct_FixH->NumHfixed;
    // INCHI✔️❌:         Ct_FixH->NumHfixed = NULL;
    // INCHI✔️❌:         pBCN->ftcn[iOther].nLenNumHFixH = num_atoms;
    // INCHI✔️❌:         pBCN->ftcn[iOther].nCanonFlags |= pCD[iOther].nCanonFlags;
    // INCHI✔️❌:
    // INCHI✔️❌:         /* original H */
    // INCHI✔️❌:         CleanNumH( pBCN->ftcn[iOther].nNumHOrig, pBCN->ftcn[iOther].nLenNumH );
    // INCHI✔️❌:         CleanNumH( pBCN->ftcn[iOther].nNumHOrigFixH, pBCN->ftcn[iOther].nLenNumH );
    // INCHI✔️❌:         /* canonical H positions */
    // INCHI✔️❌:         CleanNumH( pBCN->ftcn[iOther].nNumH, pBCN->ftcn[iOther].nLenNumH );
    // INCHI✔️❌:         CleanNumH( pBCN->ftcn[iOther].nNumHFixH, pBCN->ftcn[iOther].nLenNumH );
    // INCHI✔️❌:
    // INCHI✔️❌:         /* connection table */
    // INCHI✔️❌:         CleanCt( pBCN->ftcn[iOther].LinearCt, pBCN->ftcn[iOther].nLenLinearCt );
    // INCHI✔️❌:
    // INCHI✔️❌:        /* isotopic canonicalization */
    // INCHI✔️❌:         if (Ct_FixHIso)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pBCN->ftcn[iOther].PartitionCtIso.Rank = nCanonRankFixHIso;        nCanonRankFixHIso = NULL;
    // INCHI✔️❌:             pBCN->ftcn[iOther].PartitionCtIso.AtNumber = nAtomNumberCanonFixHIso;  nAtomNumberCanonFixHIso = NULL;
    // INCHI✔️❌:             pBCN->ftcn[iOther].nSymmRankCtIso = nSymmRankFixHIso;         nSymmRankFixHIso = NULL;
    // INCHI✔️❌:             pBCN->ftcn[iOther].iso_sort_keys = Ct_FixHIso->iso_sort_key; Ct_FixHIso->iso_sort_key = NULL;
    // INCHI✔️❌:             pBCN->ftcn[iOther].iso_sort_keysOrig = iso_sort_key_NoTautH;     iso_sort_key_NoTautH = NULL;
    // INCHI✔️❌:             pBCN->ftcn[iOther].len_iso_sort_keys = len_iso_sort_key_NoTautH;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:             MergeCleanIsoSortKeys( pBCN->ftcn[iOther].iso_sort_keys, Ct_FixHIso->iso_sort_key_Hfixed, pBCN->ftcn[iOther].len_iso_sort_keys );
    // INCHI✔️❌:             MergeCleanIsoSortKeys( pBCN->ftcn[iOther].iso_sort_keysOrig, iso_sort_key_Hfixed, pBCN->ftcn[iOther].len_iso_sort_keys );
    // INCHI✔️❌: #else
    // INCHI✔️❌:             CleanIsoSortKeys( pBCN->ftcn[iOther].iso_sort_keys, pBCN->ftcn[iOther].len_iso_sort_keys );
    // INCHI✔️❌:             CleanIsoSortKeys( pBCN->ftcn[iOther].iso_sort_keysOrig, pBCN->ftcn[iOther].len_iso_sort_keys );
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }  /* non-tautomeric results */
    // INCHI✔️❌:
    // INCHI✔️❌:     goto exit_function;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: exit_error_alloc:
    // INCHI✔️❌:     ret = CT_OUT_OF_RAM;
    // INCHI✔️❌:     goto exit_function;
    // INCHI✔️❌:
    // INCHI✔️❌: exit_error:
    // INCHI✔️❌:     if (!RETURNED_ERROR( ret ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = CT_CANON_ERR;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     goto exit_function;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:
    // INCHI✔️❌: #define FREE_CONTABLE( X) if (X) {CTableFree( X);inchi_free( X);}
    // INCHI✔️❌: #define FREE_ARRAY( X) if (X) inchi_free( X);
    // INCHI✔️❌:
    // INCHI✔️❌:     FreeNeighList( NeighList[TAUT_NON] );
    // INCHI✔️❌:     FreeNeighList( NeighList[TAUT_YES] );
    // INCHI✔️❌:
    // INCHI✔️❌:     FREE_CONTABLE( Ct_NoH )
    // INCHI✔️❌:         FREE_CONTABLE( Ct_NoTautH )
    // INCHI✔️❌:         FREE_CONTABLE( Ct_Base )
    // INCHI✔️❌:         FREE_CONTABLE( Ct_FixH )
    // INCHI✔️❌:         FREE_CONTABLE( Ct_Temp )
    // INCHI✔️❌:         /* isotopic canonicalization */
    // INCHI✔️❌:         FREE_CONTABLE( Ct_NoTautHIso )
    // INCHI✔️❌:         FREE_CONTABLE( Ct_BaseIso )
    // INCHI✔️❌:         FREE_CONTABLE( Ct_FixHIso )
    // INCHI✔️❌:
    // INCHI✔️❌:         /* free the first two pointers from pBCN->pRankStack */
    // INCHI✔️❌:         FREE_ARRAY( nRank )
    // INCHI✔️❌:         FREE_ARRAY( nAtomNumber )
    // INCHI✔️❌:
    // INCHI✔️❌:         if (pBCN->pRankStack)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pBCN->pRankStack[0] =
    // INCHI✔️❌:                 pBCN->pRankStack[1] = NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_AUX_RANKING == 1 )
    // INCHI✔️❌:     FREE_ARRAY( nRankAux )
    // INCHI✔️❌:         FREE_ARRAY( nAtomNumberAux )
    // INCHI✔️❌:         FREE_ARRAY( pAtomInvariantAux )
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         FREE_ARRAY( pAtomInvariant )
    // INCHI✔️❌:
    // INCHI✔️❌:         FREE_ARRAY( nCanonRankNoH )
    // INCHI✔️❌:         FREE_ARRAY( nAtomNumberCanonNoH )
    // INCHI✔️❌:         FREE_ARRAY( nSymmRankNoH )
    // INCHI✔️❌:
    // INCHI✔️❌:         FREE_ARRAY( nSymmRankNoTautH )
    // INCHI✔️❌:         FREE_ARRAY( nCanonRankNoTautH )
    // INCHI✔️❌:         FREE_ARRAY( nAtomNumberCanonNoTautH )
    // INCHI✔️❌:         FREE_ARRAY( numHNoTautH )
    // INCHI✔️❌:
    // INCHI✔️❌:         FREE_ARRAY( nSymmRankBase )
    // INCHI✔️❌:         FREE_ARRAY( nCanonRankBase )
    // INCHI✔️❌:         FREE_ARRAY( nAtomNumberCanonBase )
    // INCHI✔️❌:         FREE_ARRAY( numH )
    // INCHI✔️❌:
    // INCHI✔️❌:         FREE_ARRAY( nSymmRankFixH )
    // INCHI✔️❌:         FREE_ARRAY( nCanonRankFixH )
    // INCHI✔️❌:         FREE_ARRAY( nAtomNumberCanonFixH )
    // INCHI✔️❌:         FREE_ARRAY( NumHfixed )
    // INCHI✔️❌:
    // INCHI✔️❌:         /* isotopic canonicalization */
    // INCHI✔️❌:
    // INCHI✔️❌:         FREE_ARRAY( nSymmRankNoTautHIso )
    // INCHI✔️❌:         FREE_ARRAY( nCanonRankNoTautHIso )
    // INCHI✔️❌:         FREE_ARRAY( nAtomNumberCanonNoTautHIso )
    // INCHI✔️❌:         FREE_ARRAY( iso_sort_key_NoTautH )
    // INCHI✔️❌:
    // INCHI✔️❌:         FREE_ARRAY( nSymmRankBaseIso )
    // INCHI✔️❌:         FREE_ARRAY( nCanonRankBaseIso )
    // INCHI✔️❌:         FREE_ARRAY( nAtomNumberCanonBaseIso )
    // INCHI✔️❌:         FREE_ARRAY( iso_sort_keyBase )
    // INCHI✔️❌:         FREE_ARRAY( iso_exchg_atnos )
    // INCHI✔️❌:
    // INCHI✔️❌:         FREE_ARRAY( nSymmRankFixHIso )
    // INCHI✔️❌:         FREE_ARRAY( nCanonRankFixHIso )
    // INCHI✔️❌:         FREE_ARRAY( nAtomNumberCanonFixHIso )
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    // INCHI✔️❌:         FREE_ARRAY( iso_sort_key_Hfixed )
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         FREE_ARRAY( nTempRank )
    // INCHI✔️❌:
    // INCHI✔️❌: #undef FREE_CONTABLE
    // INCHI✔️❌: #undef FREE_ARRAY
    // INCHI✔️❌:
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: GetBaseCanonRanking
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetBaseCanonRanking
    // INCHI✔️❌: #define COMPILE_ANSI_ONLY
    // INCHI✔️❌: #define TARGET_API_LIB
    // INCHI✔️❌: GCC/Linux branch
    // INCHI✔️❌: #define USE_AUX_RANKING 1
    // INCHI✔️❌: #define USE_ISO_SORT_KEY_HFIXED 0
    // END INCHI ACTIVE MACRO CONFIGURATION: GetBaseCanonRanking

    // CTableCreate does not require nextAtRank, iso_sort_key, or
    // iso_exchg_atnos for success, although CanonGraph later dereferences each
    // requested field. Allocation failure at those ordinals has no defined C
    // result. Rust returns CT_OUT_OF_RAM before dereference and runs the source
    // cleanup order. The behavioral marker covers the source-defined state
    // space; it does not claim a return-value match for undefined C execution.

    let mut ret = 0_i32;
    let mut i_base = TAUT_NON as usize;
    let mut i_other = TAUT_NON as usize;
    let mut b_req_non_taut = 0_i32;
    let mut b_req_taut = 0_i32;
    let mut b_taut_ignore_isotopic = 0_i32;
    let mut b_use_iso_aux_base = [false; TAUT_NUM as usize];
    let num_max = num_at_tg;
    let mut l_count = 0_i64;
    let mut canon_data = std::array::from_fn::<_, 3, _>(|_| CANON_DATA::default());
    let mut canon_counts = CANON_COUNTS::default();

    let mut p_atom_invariant = SourceMutPointer::<ATOM_INVARIANT2>::null();
    let mut n_rank = SourceMutPointer::<AT_RANK>::null();
    let mut n_atom_number = SourceMutPointer::<AT_NUMB>::null();
    let mut n_temp_rank = SourceMutPointer::<AT_RANK>::null();
    let mut n_rank_aux = SourceMutPointer::<AT_RANK>::null();
    let mut n_atom_number_aux = SourceMutPointer::<AT_NUMB>::null();
    let mut p_atom_invariant_aux = SourceMutPointer::<ATOM_INVARIANT2>::null();
    let mut n_symm_rank_no_h = SourceMutPointer::<AT_RANK>::null();
    let mut n_canon_rank_no_h = SourceMutPointer::<AT_RANK>::null();
    let mut n_atom_number_canon_no_h = SourceMutPointer::<AT_NUMB>::null();
    let mut n_symm_rank_no_taut_h = SourceMutPointer::<AT_RANK>::null();
    let mut n_canon_rank_no_taut_h = SourceMutPointer::<AT_RANK>::null();
    let mut n_atom_number_canon_no_taut_h = SourceMutPointer::<AT_NUMB>::null();
    let mut num_h_no_taut_h = SourceMutPointer::<NUM_H>::null();
    let mut len_num_h_no_taut_h = 0_i32;
    let mut max_len_num_h_no_taut_h = 0_i32;
    let mut n_symm_rank_no_taut_h_iso = SourceMutPointer::<AT_RANK>::null();
    let mut n_canon_rank_no_taut_h_iso = SourceMutPointer::<AT_RANK>::null();
    let mut n_atom_number_canon_no_taut_h_iso = SourceMutPointer::<AT_NUMB>::null();
    let mut iso_sort_key_no_taut_h = SourceMutPointer::<AT_ISO_SORT_KEY>::null();
    let mut len_iso_sort_key_no_taut_h = 0_i32;
    let mut max_len_iso_sort_key_no_taut_h = 0_i32;
    let mut n_symm_rank_base = SourceMutPointer::<AT_RANK>::null();
    let mut n_canon_rank_base = SourceMutPointer::<AT_RANK>::null();
    let mut n_atom_number_canon_base = SourceMutPointer::<AT_NUMB>::null();
    let mut num_h = SourceMutPointer::<NUM_H>::null();
    let mut len_num_h = num_atoms;
    let mut max_len_num_h = 0_i32;
    let mut n_symm_rank_base_iso = SourceMutPointer::<AT_RANK>::null();
    let mut n_canon_rank_base_iso = SourceMutPointer::<AT_RANK>::null();
    let mut n_atom_number_canon_base_iso = SourceMutPointer::<AT_NUMB>::null();
    let mut iso_sort_key_base = SourceMutPointer::<AT_ISO_SORT_KEY>::null();
    let mut len_iso_sort_key_base = 0_i32;
    let mut max_len_iso_sort_key_base = 0_i32;
    let mut iso_exchg_atnos = SourceMutPointer::<i8>::null();
    let mut len_iso_exchg_atnos = 0_i32;
    let mut max_len_iso_exchg_atnos = 0_i32;
    let mut n_symm_rank_fix_h = SourceMutPointer::<AT_RANK>::null();
    let mut n_canon_rank_fix_h = SourceMutPointer::<AT_RANK>::null();
    let mut n_atom_number_canon_fix_h = SourceMutPointer::<AT_NUMB>::null();
    let mut num_h_fixed = SourceMutPointer::<NUM_H>::null();
    let mut max_len_num_h_fixed = 0_i32;
    let mut n_symm_rank_fix_h_iso = SourceMutPointer::<AT_RANK>::null();
    let mut n_canon_rank_fix_h_iso = SourceMutPointer::<AT_RANK>::null();
    let mut n_atom_number_canon_fix_h_iso = SourceMutPointer::<AT_NUMB>::null();
    let mut neigh_list = [SourceMutPointer::<NEIGH_LIST>::null(); TAUT_NUM as usize];
    let mut ct_no_h = SourceMutPointer::<ConTable>::null();
    let mut ct_no_taut_h = SourceMutPointer::<ConTable>::null();
    let mut ct_no_taut_h_iso = SourceMutPointer::<ConTable>::null();
    let mut ct_base = SourceMutPointer::<ConTable>::null();
    let mut ct_base_iso = SourceMutPointer::<ConTable>::null();
    let mut ct_fix_h = SourceMutPointer::<ConTable>::null();
    let mut ct_fix_h_iso = SourceMutPointer::<ConTable>::null();
    let mut ct_temp = SourceMutPointer::<ConTable>::null();

    macro_rules! alloc_or_exit {
        ($ty:ty, $count:expr, $source_size:expr, $label:lifetime) => {{
            let count = match u64::try_from($count) {
                Ok(count) => count,
                Err(_) => {
                    ret = CT_OUT_OF_RAM;
                    break $label;
                }
            };
            match inchi_calloc::<$ty>(heap, count, $source_size) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => {
                    ret = CT_OUT_OF_RAM;
                    break $label;
                }
                Err(error) => return Err(error),
            }
        }};
    }

    'source: {
        let taut_yes = TAUT_YES as usize;
        let taut_non = TAUT_NON as usize;
        let has_taut_input = !at[taut_yes].is_null() && sizes[taut_yes].nLenCT != 0;
        let taut_info_is_active = t_group_info.as_deref().is_some_and(|info| {
            sizes[taut_yes].nLenLinearCTTautomer > 0
                && ((!info.t_group.is_null() && info.num_t_groups > 0)
                    || (info.tni.bNormalizationFlags & u64::from(FLAG_NORM_CONSIDER_TAUT)) != 0
                    || (info.nNumIsotopicEndpoints > 1
                        && (info.bTautFlagsDone
                            & u64::from(
                                TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE,
                            ))
                            != 0))
        });

        if has_taut_input && taut_info_is_active {
            i_base = taut_yes;
            b_req_taut = 1;
            let info = t_group_info
                .as_deref()
                .ok_or(SourceHeapError::NullPointer)?;
            if !at[taut_non].is_null() && sizes[taut_non].nLenCT != 0 {
                i_other = taut_non;
                b_req_non_taut = 1;
            } else {
                i_other = i_base;
            }
            let use_iso_aux = sizes[i_base].nLenIsotopicEndpoints > 1
                && (info.bTautFlagsDone
                    & u64::from(TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE))
                    != 0;
            b_use_iso_aux_base[i_base] = use_iso_aux;
        } else if !at[taut_non].is_null() && sizes[taut_non].nLenCT != 0 {
            i_base = taut_non;
            i_other = i_base;
            b_req_non_taut = 1;
            num_at_tg = num_atoms;
        } else if has_taut_input {
            i_base = taut_yes;
            i_other = i_base;
            b_req_non_taut = 1;
            num_at_tg = num_atoms;
            let use_iso_aux = sizes[i_base].nLenIsotopicEndpoints > 1;
            b_use_iso_aux_base[i_base] = use_iso_aux;
        } else {
            ret = CT_UNKNOWN_ERR;
            break 'source;
        }

        if b_req_taut != 0 {
            let info = t_group_info
                .as_deref_mut()
                .ok_or(SourceHeapError::NullPointer)?;
            b_taut_ignore_isotopic = info.bIgnoreIsotopic;
            info.bIgnoreIsotopic = 1;
        }

        let at_base = at[i_base];
        p_atom_invariant = alloc_or_exit!(ATOM_INVARIANT2, num_max, 40, 'source);
        n_symm_rank_no_h = alloc_or_exit!(AT_RANK, num_max, 2, 'source);
        n_canon_rank_no_h = alloc_or_exit!(AT_RANK, num_max, 2, 'source);
        n_atom_number_canon_no_h = alloc_or_exit!(AT_NUMB, num_max, 2, 'source);
        n_rank = alloc_or_exit!(AT_RANK, num_max, 2, 'source);
        n_atom_number = alloc_or_exit!(AT_NUMB, num_max, 2, 'source);
        n_temp_rank = alloc_or_exit!(AT_RANK, num_max, 2, 'source);
        n_rank_aux = alloc_or_exit!(AT_RANK, num_max, 2, 'source);
        n_atom_number_aux = alloc_or_exit!(AT_NUMB, num_max, 2, 'source);
        p_atom_invariant_aux = alloc_or_exit!(ATOM_INVARIANT2, num_max, 40, 'source);

        if b_req_taut != 0 {
            neigh_list[taut_yes] = CreateNeighList(
                heap,
                num_atoms,
                num_at_tg,
                at_base.as_const(),
                0,
                t_group_info.as_deref(),
            )?;
            if neigh_list[taut_yes].is_null() {
                ret = CT_OUT_OF_RAM;
                break 'source;
            }
            neigh_list[taut_non] =
                CreateNeighList(heap, num_atoms, num_atoms, at_base.as_const(), 0, None)?;
        } else {
            neigh_list[taut_non] =
                CreateNeighList(heap, num_atoms, num_atoms, at_base.as_const(), 0, None)?;
        }
        if neigh_list[taut_non].is_null() {
            ret = CT_OUT_OF_RAM;
            break 'source;
        }

        pBCN.nMaxLenRankStack = 0;
        pBCN.num_max = num_max;
        pBCN.num_at_tg = num_at_tg;
        pBCN.num_atoms = num_atoms;
        pBCN.ulTimeOutTime = ulTimeOutTime;

        FillOutAtomInvariant2(
            heap,
            at_base.as_const(),
            num_atoms,
            num_atoms,
            p_atom_invariant,
            1,
            0,
            0,
            0,
            0,
            None,
        )?;
        let mut n_num_curr_ranks = SetInitialRanks2(
            heap,
            num_atoms,
            p_atom_invariant,
            n_rank,
            n_atom_number,
            pCG,
        )?;
        n_num_curr_ranks = DifferentiateRanks2(
            heap,
            pCG,
            num_atoms,
            neigh_list[taut_non],
            n_num_curr_ranks,
            n_rank,
            n_temp_rank,
            n_atom_number,
            &mut l_count,
            0,
        )?;

        let n_max_len_rank_stack = 2_i32
            .wrapping_mul(num_at_tg.wrapping_sub(n_num_curr_ranks))
            .wrapping_add(8);
        pBCN.pRankStack = alloc_or_exit!(
            SourceMutPointer<AT_RANK>,
            n_max_len_rank_stack,
            SOURCE_SIZEOF_POINTER,
            'source
        );
        pBCN.nMaxLenRankStack = n_max_len_rank_stack;
        source_set(heap, pBCN.pRankStack, 0, n_rank)?;
        source_set(heap, pBCN.pRankStack, 1, n_atom_number)?;
        let partition_count = usize::try_from(n_max_len_rank_stack / 2)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let mut partitions = vec![Partition::default(); partition_count];
        partitions[0] = Partition {
            Rank: n_rank,
            AtNumber: n_atom_number,
        };

        let data = &mut canon_data[i_other];
        data.ulTimeOutTime = pBCN.ulTimeOutTime;
        data.nMaxLenLinearCT = sizes[i_other].nLenCTAtOnly.wrapping_add(1);
        data.nLenLinearCT = sizes[i_other].nLenCTAtOnly;
        data.nLenCTAtOnly = sizes[i_other].nLenCTAtOnly;
        ret = CanonGraph(
            heap,
            ic,
            pCG,
            num_atoms,
            num_atoms,
            num_max,
            0,
            neigh_list[taut_non],
            &mut partitions,
            n_symm_rank_no_h,
            n_canon_rank_no_h,
            n_atom_number_canon_no_h,
            data,
            &mut canon_counts,
            SourceMutPointer::null(),
            Some(&mut ct_no_h),
            LargeMolecules,
            user_action,
            console_quit,
            clock_result,
        )?;
        sync_partition_stack(heap, pBCN.pRankStack, &partitions)?;
        if ret < 0 {
            if !(CT_ERR_MIN..=CT_ERR_MAX).contains(&ret) {
                ret = CT_CANON_ERR;
            }
            break 'source;
        }

        let mut changed = 0_i32;
        n_num_curr_ranks = FixCanonEquivalenceInfo(
            heap,
            pCG,
            num_atoms,
            n_symm_rank_no_h,
            n_rank,
            n_temp_rank,
            n_atom_number,
            Some(&mut changed),
        )?;
        if changed & 3 != 0 {
            if !ct_no_h.is_null() {
                let mut table = source_clone(heap, ct_no_h, 0)?;
                CTableFree(heap, Some(&mut table))?;
                inchi_free(heap, ct_no_h)?;
                ct_no_h = SourceMutPointer::null();
            }
            canon_data[i_other].nCanonFlags |= CANON_FLAG_NO_H_RECANON as i32;
            ret = CanonGraph(
                heap,
                ic,
                pCG,
                num_atoms,
                num_atoms,
                num_max,
                0,
                neigh_list[taut_non],
                &mut partitions,
                n_symm_rank_no_h,
                n_canon_rank_no_h,
                n_atom_number_canon_no_h,
                &mut canon_data[i_other],
                &mut canon_counts,
                SourceMutPointer::null(),
                Some(&mut ct_no_h),
                LargeMolecules,
                user_action,
                console_quit,
                clock_result,
            )?;
            sync_partition_stack(heap, pBCN.pRankStack, &partitions)?;
            if ret < 0 {
                if !(CT_ERR_MIN..=CT_ERR_MAX).contains(&ret) {
                    ret = CT_CANON_ERR;
                }
                break 'source;
            }
        }

        max_len_num_h_no_taut_h = num_atoms.wrapping_add(1);
        n_symm_rank_no_taut_h = alloc_or_exit!(AT_RANK, num_max, 2, 'source);
        n_canon_rank_no_taut_h = alloc_or_exit!(AT_RANK, num_max, 2, 'source);
        n_atom_number_canon_no_taut_h = alloc_or_exit!(AT_NUMB, num_max, 2, 'source);
        num_h_no_taut_h = alloc_or_exit!(NUM_H, max_len_num_h_no_taut_h, 2, 'source);
        let atom_count =
            usize::try_from(num_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        for index in 0..atom_count {
            let atom = heap
                .slice(at_base.as_const())?
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let value = if atom.endpoint == 0 && atom.num_H != 0 {
                i32::from(atom.num_H).wrapping_add(BASE_H_NUMBER as i32) as NUM_H
            } else {
                EMPTY_H_NUMBER as NUM_H
            };
            source_set(heap, num_h_no_taut_h, index as i32, value)?;
        }

        {
            let data = &mut canon_data[i_other];
            data.LinearCT = SourceMutPointer::null();
            data.NumH = num_h_no_taut_h;
            data.NumHfixed = SourceMutPointer::null();
            data.iso_sort_key = SourceMutPointer::null();
            data.iso_exchg_atnos = SourceMutPointer::null();
            data.ulTimeOutTime = pBCN.ulTimeOutTime;
            data.nMaxLenLinearCT = sizes[i_other].nLenCTAtOnly.wrapping_add(1);
            data.maxlenNumH = max_len_num_h_no_taut_h;
            data.nLenLinearCT = sizes[i_other].nLenCTAtOnly;
            data.nLenCTAtOnly = sizes[i_other].nLenCTAtOnly;
            len_num_h_no_taut_h = num_atoms;
            data.lenNumH = len_num_h_no_taut_h;
            data.lenNumHfixed = 0;
            data.len_iso_sort_key = 0;
            data.maxlen_iso_sort_key = 0;
            data.len_iso_exchg_atnos = 0;
            data.maxlen_iso_exchg_atnos = 0;
            data.nAuxRank = SourceMutPointer::null();
        }

        for value in heap
            .slice_mut(n_temp_rank)?
            .get_mut(..usize::try_from(num_max).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
        {
            *value = 0;
        }
        for index in 0..atom_count {
            let symmetry = source_get(heap, n_symm_rank_no_h, index as i32)?;
            let class_index = usize::from(symmetry.wrapping_sub(1));
            let representatives = heap.slice_mut(n_temp_rank)?;
            if representatives
                .get(class_index)
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                < index as AT_RANK
            {
                representatives[class_index] = index as AT_RANK;
            }
        }
        let mut needs_no_taut_h_canon = false;
        for index in 0..atom_count {
            let symmetry = source_get(heap, n_symm_rank_no_h, index as i32)?;
            let representative = source_get(heap, n_temp_rank, i32::from(symmetry) - 1)?;
            if source_get(heap, num_h_no_taut_h, index as i32)?
                != source_get(heap, num_h_no_taut_h, i32::from(representative))?
            {
                canon_data[i_other].nCanonFlags |= CANON_FLAG_NO_TAUT_H_DIFF as i32;
                needs_no_taut_h_canon = true;
                break;
            }
        }

        if needs_no_taut_h_canon {
            for value in heap
                .slice_mut(p_atom_invariant_aux)?
                .get_mut(..atom_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
            {
                *value = ATOM_INVARIANT2::default();
            }
            for index in 0..atom_count {
                let symmetry = source_get(heap, n_symm_rank_no_h, index as i32)?;
                let hydrogen = source_get(heap, num_h_no_taut_h, index as i32)? as AT_NUMB;
                let invariant = heap
                    .slice_mut(p_atom_invariant_aux)?
                    .get_mut(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                invariant.val[0] = symmetry;
                invariant.val[1] = hydrogen;
            }
            n_num_curr_ranks = SetInitialRanks2(
                heap,
                num_atoms,
                p_atom_invariant_aux,
                n_rank_aux,
                n_atom_number_aux,
                pCG,
            )?;
            n_num_curr_ranks = DifferentiateRanks2(
                heap,
                pCG,
                num_atoms,
                neigh_list[taut_non],
                n_num_curr_ranks,
                n_rank_aux,
                n_temp_rank,
                n_atom_number_aux,
                &mut l_count,
                0,
            )?;
            canon_data[i_other].nAuxRank = n_rank_aux;
            ret = CanonGraph(
                heap,
                ic,
                pCG,
                num_atoms,
                num_atoms,
                num_max,
                1,
                neigh_list[taut_non],
                &mut partitions,
                n_symm_rank_no_taut_h,
                n_canon_rank_no_taut_h,
                n_atom_number_canon_no_taut_h,
                &mut canon_data[i_other],
                &mut canon_counts,
                ct_no_h,
                Some(&mut ct_no_taut_h),
                LargeMolecules,
                user_action,
                console_quit,
                clock_result,
            )?;
            sync_partition_stack(heap, pBCN.pRankStack, &partitions)?;
            if ret < 0 {
                if !(CT_ERR_MIN..=CT_ERR_MAX).contains(&ret) {
                    ret = CT_CANON_ERR;
                }
                break 'source;
            }
        } else {
            ct_temp = alloc_or_exit!(ConTable, 1, 104, 'source);
            let mut temporary = ConTable::default();
            if CTableCreate(heap, &mut temporary, num_atoms, &canon_data[i_other])? == 0
                || !con_table_is_usable(&temporary, &canon_data[i_other])
            {
                ret = CT_OUT_OF_RAM;
                heap.slice_mut(ct_temp)?[0] = temporary;
                break 'source;
            }
            let no_h = source_clone(heap, ct_no_h, 0)?;
            CtFullCopy(heap, &mut temporary, &no_h)?;
            for index in 0..atom_count {
                let rank = source_get(heap, n_canon_rank_no_h, index as i32)?;
                source_set(
                    heap,
                    temporary.NumH,
                    i32::from(rank) - 1,
                    source_get(heap, num_h_no_taut_h, index as i32)?,
                )?;
            }
            temporary.lenNumH = num_atoms;
            heap.slice_mut(ct_temp)?[0] = temporary;
            ct_no_taut_h = ct_temp;
            ct_temp = SourceMutPointer::null();
            for index in 0..atom_count {
                source_set(
                    heap,
                    n_symm_rank_no_taut_h,
                    index as i32,
                    source_get(heap, n_symm_rank_no_h, index as i32)?,
                )?;
                source_set(
                    heap,
                    n_canon_rank_no_taut_h,
                    index as i32,
                    source_get(heap, n_canon_rank_no_h, index as i32)?,
                )?;
                source_set(
                    heap,
                    n_atom_number_canon_no_taut_h,
                    index as i32,
                    source_get(heap, n_atom_number_canon_no_h, index as i32)?,
                )?;
            }
        }

        if sizes[i_other].num_isotopic_atoms != 0
            && sizes[i_other].bIgnoreIsotopic == 0
            && b_req_taut == 0
            && b_req_non_taut != 0
        {
            max_len_iso_sort_key_no_taut_h = num_atoms.wrapping_add(1);
            n_symm_rank_no_taut_h_iso = alloc_or_exit!(AT_RANK, num_max, 2, 'source);
            n_canon_rank_no_taut_h_iso = alloc_or_exit!(AT_RANK, num_max, 2, 'source);
            n_atom_number_canon_no_taut_h_iso = alloc_or_exit!(AT_NUMB, num_max, 2, 'source);
            iso_sort_key_no_taut_h = alloc_or_exit!(
                AT_ISO_SORT_KEY,
                max_len_iso_sort_key_no_taut_h,
                8,
                'source
            );

            let mut num_iso_no_taut_h = 0_i32;
            for index in 0..atom_count {
                let atom = heap
                    .slice(at_base.as_const())?
                    .get(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let key = if atom.endpoint != 0 {
                    make_iso_sort_key(i32::from(atom.iso_atw_diff), 0, 0, 0)
                } else {
                    make_iso_sort_key(
                        i32::from(atom.iso_atw_diff),
                        i32::from(atom.num_iso_H[0]),
                        i32::from(atom.num_iso_H[1]),
                        i32::from(atom.num_iso_H[2]),
                    )
                };
                if key != 0 {
                    num_iso_no_taut_h = num_iso_no_taut_h.wrapping_add(1);
                    source_set(heap, iso_sort_key_no_taut_h, index as i32, key)?;
                } else {
                    source_set(
                        heap,
                        iso_sort_key_no_taut_h,
                        index as i32,
                        AT_ISO_SORT_KEY::MAX,
                    )?;
                }
            }

            {
                let data = &mut canon_data[i_other];
                data.LinearCT = SourceMutPointer::null();
                data.NumH = num_h_no_taut_h;
                data.NumHfixed = SourceMutPointer::null();
                data.iso_sort_key = iso_sort_key_no_taut_h;
                data.iso_exchg_atnos = SourceMutPointer::null();
                data.ulTimeOutTime = pBCN.ulTimeOutTime;
                data.nMaxLenLinearCT = sizes[i_other].nLenCTAtOnly.wrapping_add(1);
                data.maxlenNumH = max_len_num_h_no_taut_h;
                data.nLenLinearCT = sizes[i_other].nLenCTAtOnly;
                data.nLenCTAtOnly = sizes[i_other].nLenCTAtOnly;
                data.lenNumH = len_num_h_no_taut_h;
                data.lenNumHfixed = 0;
                len_iso_sort_key_no_taut_h = num_atoms;
                data.len_iso_sort_key = len_iso_sort_key_no_taut_h;
                data.maxlen_iso_sort_key = max_len_iso_sort_key_no_taut_h;
                data.len_iso_exchg_atnos = 0;
                data.maxlen_iso_exchg_atnos = 0;
                data.nAuxRank = SourceMutPointer::null();
            }

            let mut needs_isotope_canon = false;
            if num_iso_no_taut_h != 0 {
                for value in heap
                    .slice_mut(n_temp_rank)?
                    .get_mut(
                        ..usize::try_from(num_max)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                {
                    *value = 0;
                }
                for index in 0..atom_count {
                    let symmetry = source_get(heap, n_symm_rank_no_taut_h, index as i32)?;
                    let class_index = i32::from(symmetry) - 1;
                    if source_get(heap, n_temp_rank, class_index)? < index as AT_RANK {
                        source_set(heap, n_temp_rank, class_index, index as AT_RANK)?;
                    }
                }
                for index in 0..atom_count {
                    let symmetry = source_get(heap, n_symm_rank_no_taut_h, index as i32)?;
                    let representative = source_get(heap, n_temp_rank, i32::from(symmetry) - 1)?;
                    if source_get(heap, iso_sort_key_no_taut_h, index as i32)?
                        != source_get(heap, iso_sort_key_no_taut_h, i32::from(representative))?
                    {
                        canon_data[i_other].nCanonFlags |= CANON_FLAG_ISO_ONLY_NON_TAUT_DIFF as i32;
                        needs_isotope_canon = true;
                        break;
                    }
                }
            }

            if needs_isotope_canon {
                for value in heap
                    .slice_mut(p_atom_invariant_aux)?
                    .get_mut(..atom_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                {
                    *value = ATOM_INVARIANT2::default();
                }
                for index in 0..atom_count {
                    let symmetry = source_get(heap, n_symm_rank_no_taut_h, index as i32)?;
                    let isotope = source_get(heap, iso_sort_key_no_taut_h, index as i32)?;
                    let invariant = heap
                        .slice_mut(p_atom_invariant_aux)?
                        .get_mut(index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    invariant.val[0] = symmetry;
                    invariant.iso_sort_key = isotope;
                }
                n_num_curr_ranks = SetInitialRanks2(
                    heap,
                    num_atoms,
                    p_atom_invariant_aux,
                    n_rank_aux,
                    n_atom_number_aux,
                    pCG,
                )?;
                n_num_curr_ranks = DifferentiateRanks2(
                    heap,
                    pCG,
                    num_atoms,
                    neigh_list[taut_non],
                    n_num_curr_ranks,
                    n_rank_aux,
                    n_temp_rank,
                    n_atom_number_aux,
                    &mut l_count,
                    0,
                )?;
                canon_data[i_other].nAuxRank = n_rank_aux;
                ret = CanonGraph(
                    heap,
                    ic,
                    pCG,
                    num_atoms,
                    num_atoms,
                    num_max,
                    1,
                    neigh_list[taut_non],
                    &mut partitions,
                    n_symm_rank_no_taut_h_iso,
                    n_canon_rank_no_taut_h_iso,
                    n_atom_number_canon_no_taut_h_iso,
                    &mut canon_data[i_other],
                    &mut canon_counts,
                    ct_no_taut_h,
                    Some(&mut ct_no_taut_h_iso),
                    LargeMolecules,
                    user_action,
                    console_quit,
                    clock_result,
                )?;
                sync_partition_stack(heap, pBCN.pRankStack, &partitions)?;
                if ret < 0 {
                    if !(CT_ERR_MIN..=CT_ERR_MAX).contains(&ret) {
                        ret = CT_CANON_ERR;
                    }
                    break 'source;
                }
            } else {
                ct_temp = alloc_or_exit!(ConTable, 1, 104, 'source);
                let mut temporary = ConTable::default();
                if CTableCreate(heap, &mut temporary, num_atoms, &canon_data[i_other])? == 0
                    || !con_table_is_usable(&temporary, &canon_data[i_other])
                {
                    ret = CT_OUT_OF_RAM;
                    heap.slice_mut(ct_temp)?[0] = temporary;
                    break 'source;
                }
                let no_taut_h = source_clone(heap, ct_no_taut_h, 0)?;
                CtFullCopy(heap, &mut temporary, &no_taut_h)?;
                for index in 0..atom_count {
                    let rank = source_get(heap, n_canon_rank_no_taut_h, index as i32)?;
                    source_set(
                        heap,
                        temporary.iso_sort_key,
                        i32::from(rank) - 1,
                        source_get(heap, iso_sort_key_no_taut_h, index as i32)?,
                    )?;
                }
                temporary.len_iso_sort_key = num_atoms;
                heap.slice_mut(ct_temp)?[0] = temporary;
                ct_no_taut_h_iso = ct_temp;
                ct_temp = SourceMutPointer::null();
                for index in 0..atom_count {
                    source_set(
                        heap,
                        n_symm_rank_no_taut_h_iso,
                        index as i32,
                        source_get(heap, n_symm_rank_no_taut_h, index as i32)?,
                    )?;
                    source_set(
                        heap,
                        n_canon_rank_no_taut_h_iso,
                        index as i32,
                        source_get(heap, n_canon_rank_no_taut_h, index as i32)?,
                    )?;
                    source_set(
                        heap,
                        n_atom_number_canon_no_taut_h_iso,
                        index as i32,
                        source_get(heap, n_atom_number_canon_no_taut_h, index as i32)?,
                    )?;
                }
            }
        }

        if b_req_taut != 0 {
            max_len_num_h = num_atoms
                .wrapping_add(
                    (T_NUM_NO_ISOTOPIC as i32).wrapping_mul(num_at_tg.wrapping_sub(num_atoms)),
                )
                .wrapping_add(1);
            n_symm_rank_base = alloc_or_exit!(AT_RANK, num_max, 2, 'source);
            n_canon_rank_base = alloc_or_exit!(AT_RANK, num_max, 2, 'source);
            n_atom_number_canon_base = alloc_or_exit!(AT_NUMB, num_max, 2, 'source);
            num_h = alloc_or_exit!(NUM_H, max_len_num_h, 2, 'source);

            for index in 0..atom_count {
                let atom = heap
                    .slice(at_base.as_const())?
                    .get(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let value = if atom.endpoint == 0 && atom.num_H != 0 {
                    i32::from(atom.num_H).wrapping_add(BASE_H_NUMBER as i32) as NUM_H
                } else {
                    EMPTY_H_NUMBER as NUM_H
                };
                source_set(heap, num_h, index as i32, value)?;
            }
            let info = t_group_info
                .as_deref()
                .ok_or(SourceHeapError::NullPointer)?;
            let mut next_h = num_atoms;
            let mut index = num_atoms;
            while index < num_at_tg {
                let group = source_clone(heap, info.t_group, index.wrapping_sub(num_atoms))?;
                for value in group.num.iter().take(T_NUM_NO_ISOTOPIC as usize) {
                    source_set(
                        heap,
                        num_h,
                        next_h,
                        if *value != 0 {
                            i32::from(*value).wrapping_add(BASE_H_NUMBER as i32) as NUM_H
                        } else {
                            EMPTY_H_NUMBER as NUM_H
                        },
                    )?;
                    next_h = next_h.wrapping_add(1);
                }
                index = index.wrapping_add(1);
            }
            len_num_h = next_h;

            {
                let data = &mut canon_data[i_base];
                data.LinearCT = SourceMutPointer::null();
                data.NumH = num_h;
                data.NumHfixed = SourceMutPointer::null();
                data.iso_sort_key = SourceMutPointer::null();
                data.iso_exchg_atnos = SourceMutPointer::null();
                data.ulTimeOutTime = pBCN.ulTimeOutTime;
                data.nMaxLenLinearCT = sizes[i_base].nLenCT.wrapping_add(1);
                data.maxlenNumH = max_len_num_h;
                data.nLenLinearCT = sizes[i_base].nLenCT;
                data.nLenCTAtOnly = sizes[i_base].nLenCTAtOnly;
                data.lenNumH = len_num_h;
                data.lenNumHfixed = 0;
                data.len_iso_sort_key = 0;
                data.maxlen_iso_sort_key = 0;
                data.len_iso_exchg_atnos = 0;
                data.maxlen_iso_exchg_atnos = 0;
                data.nAuxRank = SourceMutPointer::null();
            }

            FillOutAtomInvariant2(
                heap,
                at_base.as_const(),
                num_atoms,
                num_at_tg,
                p_atom_invariant,
                1,
                0,
                0,
                1,
                1,
                t_group_info.as_deref(),
            )?;
            for index in 0..atom_count {
                let rank = source_get(heap, partitions[0].Rank, index as i32)?;
                heap.slice_mut(p_atom_invariant)?[index].val[0] = rank;
            }
            n_num_curr_ranks = SetInitialRanks2(
                heap,
                num_at_tg,
                p_atom_invariant,
                n_rank,
                n_atom_number,
                pCG,
            )?;
            n_num_curr_ranks = DifferentiateRanks4(
                heap,
                pCG,
                num_at_tg,
                neigh_list[taut_yes],
                n_num_curr_ranks,
                partitions[0].Rank,
                n_temp_rank,
                partitions[0].AtNumber,
                num_atoms as AT_RANK,
                &mut l_count,
            )?;

            for value in heap
                .slice_mut(p_atom_invariant_aux)?
                .get_mut(
                    ..usize::try_from(num_max).map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )
                .ok_or(SourceHeapError::PointerOutOfBounds)?
            {
                *value = ATOM_INVARIANT2::default();
            }
            for index in 0..atom_count {
                let symmetry = source_get(heap, n_symm_rank_no_taut_h, index as i32)?;
                let hydrogen = source_get(heap, num_h, index as i32)? as AT_NUMB;
                let invariant = &mut heap.slice_mut(p_atom_invariant_aux)?[index];
                invariant.val[0] = symmetry;
                invariant.val[1] = hydrogen;
            }
            let total_count =
                usize::try_from(num_at_tg).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            for index in atom_count..total_count {
                heap.slice_mut(p_atom_invariant_aux)?[index].val[0] =
                    source_get(heap, n_rank, index as i32)?;
            }
            n_num_curr_ranks = SetInitialRanks2(
                heap,
                num_at_tg,
                p_atom_invariant_aux,
                n_rank_aux,
                n_atom_number_aux,
                pCG,
            )?;
            n_num_curr_ranks = DifferentiateRanks4(
                heap,
                pCG,
                num_at_tg,
                neigh_list[taut_yes],
                n_num_curr_ranks,
                n_rank_aux,
                n_temp_rank,
                n_atom_number_aux,
                num_atoms as AT_RANK,
                &mut l_count,
            )?;
            canon_data[i_base].nAuxRank = n_rank_aux;
            ret = CanonGraph(
                heap,
                ic,
                pCG,
                num_atoms,
                num_at_tg,
                num_max,
                1,
                neigh_list[taut_yes],
                &mut partitions,
                n_symm_rank_base,
                n_canon_rank_base,
                n_atom_number_canon_base,
                &mut canon_data[i_base],
                &mut canon_counts,
                ct_no_taut_h,
                Some(&mut ct_base),
                LargeMolecules,
                user_action,
                console_quit,
                clock_result,
            )?;
            sync_partition_stack(heap, pBCN.pRankStack, &partitions)?;
            if ret < 0 {
                if !(CT_ERR_MIN..=CT_ERR_MAX).contains(&ret) {
                    ret = CT_CANON_ERR;
                }
                break 'source;
            }

            let needs_taut_isotope = (sizes[i_base].num_isotopic_atoms != 0
                && sizes[i_base].bIgnoreIsotopic == 0)
                || (sizes[i_base].bHasIsotopicTautGroups != 0 && b_taut_ignore_isotopic == 0)
                || (b_use_iso_aux_base[i_base] && b_taut_ignore_isotopic == 0);
            if needs_taut_isotope {
                t_group_info
                    .as_deref_mut()
                    .ok_or(SourceHeapError::NullPointer)?
                    .bIgnoreIsotopic = b_taut_ignore_isotopic;
                n_symm_rank_base_iso = alloc_or_exit!(AT_RANK, num_max, 2, 'source);
                n_canon_rank_base_iso = alloc_or_exit!(AT_RANK, num_max, 2, 'source);
                n_atom_number_canon_base_iso = alloc_or_exit!(AT_NUMB, num_max, 2, 'source);
                if b_use_iso_aux_base[i_base] {
                    max_len_iso_exchg_atnos = num_max.wrapping_add(1);
                    iso_exchg_atnos = alloc_or_exit!(i8, max_len_iso_exchg_atnos, 1, 'source);
                }
                max_len_iso_sort_key_base = num_max.wrapping_add(1);
                iso_sort_key_base = alloc_or_exit!(
                    AT_ISO_SORT_KEY,
                    max_len_iso_sort_key_base,
                    8,
                    'source
                );

                let mut num_iso_no_taut_h = 0_i32;
                let mut num_iso_no_aux_base = 0_i32;
                if !iso_exchg_atnos.is_null() {
                    len_iso_exchg_atnos = num_at_tg;
                }
                for index in 0..atom_count {
                    let atom = heap.slice(at_base.as_const())?[index].clone();
                    let key = if atom.endpoint != 0
                        || (!iso_exchg_atnos.is_null()
                            && (atom.cFlags & AT_FLAG_ISO_H_POINT as i8) != 0)
                    {
                        if !iso_exchg_atnos.is_null() && atom.endpoint == 0 {
                            num_iso_no_aux_base = num_iso_no_aux_base.wrapping_add(1);
                        }
                        make_iso_sort_key(i32::from(atom.iso_atw_diff), 0, 0, 0)
                    } else {
                        if !iso_exchg_atnos.is_null() {
                            source_set(heap, iso_exchg_atnos, index as i32, 1)?;
                        }
                        make_iso_sort_key(
                            i32::from(atom.iso_atw_diff),
                            i32::from(atom.num_iso_H[0]),
                            i32::from(atom.num_iso_H[1]),
                            i32::from(atom.num_iso_H[2]),
                        )
                    };
                    if key != 0 {
                        num_iso_no_taut_h = num_iso_no_taut_h.wrapping_add(1);
                        source_set(heap, iso_sort_key_base, index as i32, key)?;
                    } else {
                        source_set(heap, iso_sort_key_base, index as i32, AT_ISO_SORT_KEY::MAX)?;
                    }
                }

                if !iso_exchg_atnos.is_null() {
                    let info = t_group_info
                        .as_deref()
                        .ok_or(SourceHeapError::NullPointer)?;
                    if num_iso_no_aux_base
                        != i32::from(source_get(heap, info.nIsotopicEndpointAtomNumber, 0)?)
                    {
                        ret = CT_ISOCOUNT_ERR;
                        break 'source;
                    }
                    let mut endpoint = 1_i32;
                    while endpoint <= num_iso_no_aux_base {
                        let atom_index = i32::from(source_get(
                            heap,
                            info.nIsotopicEndpointAtomNumber,
                            endpoint,
                        )?);
                        let atom = source_clone(heap, at_base, atom_index)?;
                        if atom.endpoint != 0 || (atom.cFlags & AT_FLAG_ISO_H_POINT as i8) == 0 {
                            ret = CT_ISOCOUNT_ERR;
                            break 'source;
                        }
                        endpoint = endpoint.wrapping_add(1);
                    }
                }

                let mut num_iso_base = 0_i32;
                if !iso_exchg_atnos.is_null() {
                    for index in num_atoms..num_at_tg {
                        source_set(heap, iso_sort_key_base, index, AT_ISO_SORT_KEY::MAX)?;
                    }
                } else {
                    let info = t_group_info
                        .as_deref()
                        .ok_or(SourceHeapError::NullPointer)?;
                    for index in num_atoms..num_at_tg {
                        let group =
                            source_clone(heap, info.t_group, index.wrapping_sub(num_atoms))?;
                        if group.iWeight != 0 {
                            num_iso_base = num_iso_base.wrapping_add(1);
                            source_set(heap, iso_sort_key_base, index, group.iWeight)?;
                        } else {
                            source_set(heap, iso_sort_key_base, index, AT_ISO_SORT_KEY::MAX)?;
                        }
                    }
                }
                if num_iso_no_aux_base == 0 && !iso_exchg_atnos.is_null() {
                    inchi_free(heap, iso_exchg_atnos)?;
                    iso_exchg_atnos = SourceMutPointer::null();
                    len_iso_exchg_atnos = 0;
                    max_len_iso_exchg_atnos = 0;
                }
                if num_iso_no_taut_h == 0 && num_iso_base == 0 {
                    inchi_free(heap, iso_sort_key_base)?;
                    iso_sort_key_base = SourceMutPointer::null();
                    max_len_iso_sort_key_base = 0;
                } else {
                    len_iso_sort_key_base = num_at_tg;
                }

                if iso_exchg_atnos.is_null() && iso_sort_key_base.is_null() {
                    inchi_free(heap, n_symm_rank_base_iso)?;
                    inchi_free(heap, n_canon_rank_base_iso)?;
                    inchi_free(heap, n_atom_number_canon_base_iso)?;
                    n_symm_rank_base_iso = SourceMutPointer::null();
                    n_canon_rank_base_iso = SourceMutPointer::null();
                    n_atom_number_canon_base_iso = SourceMutPointer::null();
                } else {
                    {
                        let data = &mut canon_data[i_base];
                        data.LinearCT = SourceMutPointer::null();
                        data.NumH = num_h;
                        data.NumHfixed = SourceMutPointer::null();
                        data.iso_sort_key = iso_sort_key_base;
                        data.iso_exchg_atnos = iso_exchg_atnos;
                        data.ulTimeOutTime = pBCN.ulTimeOutTime;
                        data.nMaxLenLinearCT = sizes[i_base].nLenCT.wrapping_add(1);
                        data.maxlenNumH = max_len_num_h;
                        data.nLenLinearCT = sizes[i_base].nLenCT;
                        data.nLenCTAtOnly = sizes[i_base].nLenCTAtOnly;
                        data.lenNumH = len_num_h;
                        data.lenNumHfixed = 0;
                        data.len_iso_sort_key = len_iso_sort_key_base;
                        data.maxlen_iso_sort_key = max_len_iso_sort_key_base;
                        data.len_iso_exchg_atnos = len_iso_exchg_atnos;
                        data.maxlen_iso_exchg_atnos = max_len_iso_exchg_atnos;
                        data.nAuxRank = SourceMutPointer::null();
                    }
                    let mut needs_iso_canon = false;
                    if num_iso_no_taut_h != 0 || num_iso_base != 0 || num_iso_no_aux_base != 0 {
                        for value in heap
                            .slice_mut(n_temp_rank)?
                            .get_mut(
                                ..usize::try_from(num_max)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                        {
                            *value = 0;
                        }
                        for index in 0..total_count {
                            let rank = source_get(heap, n_symm_rank_base, index as i32)?;
                            if source_get(heap, n_temp_rank, i32::from(rank) - 1)?
                                < index as AT_RANK
                            {
                                source_set(
                                    heap,
                                    n_temp_rank,
                                    i32::from(rank) - 1,
                                    index as AT_RANK,
                                )?;
                            }
                        }
                        for index in 0..total_count {
                            let rank = source_get(heap, n_symm_rank_base, index as i32)?;
                            let representative =
                                source_get(heap, n_temp_rank, i32::from(rank) - 1)?;
                            let isotope_differs = !iso_sort_key_base.is_null()
                                && source_get(heap, iso_sort_key_base, index as i32)?
                                    != source_get(
                                        heap,
                                        iso_sort_key_base,
                                        i32::from(representative),
                                    )?;
                            let exchange_differs = !iso_exchg_atnos.is_null()
                                && source_get(heap, iso_exchg_atnos, index as i32)?
                                    != source_get(
                                        heap,
                                        iso_exchg_atnos,
                                        i32::from(representative),
                                    )?;
                            if isotope_differs || exchange_differs {
                                canon_data[i_base].nCanonFlags |= CANON_FLAG_ISO_TAUT_DIFF as i32;
                                needs_iso_canon = true;
                                break;
                            }
                        }
                    }

                    if needs_iso_canon {
                        for value in heap
                            .slice_mut(p_atom_invariant_aux)?
                            .get_mut(
                                ..usize::try_from(num_max)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                        {
                            *value = ATOM_INVARIANT2::default();
                        }
                        for index in 0..total_count {
                            let symmetry = source_get(heap, n_symm_rank_base, index as i32)?;
                            let isotope = if iso_sort_key_base.is_null() {
                                0
                            } else {
                                source_get(heap, iso_sort_key_base, index as i32)?
                            };
                            let auxiliary = if iso_exchg_atnos.is_null() {
                                0
                            } else {
                                source_get(heap, iso_exchg_atnos, index as i32)?
                            };
                            let invariant = &mut heap.slice_mut(p_atom_invariant_aux)?[index];
                            invariant.val[0] = symmetry;
                            invariant.iso_sort_key = isotope;
                            invariant.iso_aux_key = auxiliary;
                        }
                        n_num_curr_ranks = SetInitialRanks2(
                            heap,
                            num_at_tg,
                            p_atom_invariant_aux,
                            n_rank_aux,
                            n_atom_number_aux,
                            pCG,
                        )?;
                        n_num_curr_ranks = DifferentiateRanks2(
                            heap,
                            pCG,
                            num_at_tg,
                            neigh_list[taut_yes],
                            n_num_curr_ranks,
                            n_rank_aux,
                            n_temp_rank,
                            n_atom_number_aux,
                            &mut l_count,
                            0,
                        )?;
                        canon_data[i_base].nAuxRank = n_rank_aux;
                        ret = CanonGraph(
                            heap,
                            ic,
                            pCG,
                            num_atoms,
                            num_at_tg,
                            num_max,
                            1,
                            neigh_list[taut_yes],
                            &mut partitions,
                            n_symm_rank_base_iso,
                            n_canon_rank_base_iso,
                            n_atom_number_canon_base_iso,
                            &mut canon_data[i_base],
                            &mut canon_counts,
                            ct_base,
                            Some(&mut ct_base_iso),
                            LargeMolecules,
                            user_action,
                            console_quit,
                            clock_result,
                        )?;
                        sync_partition_stack(heap, pBCN.pRankStack, &partitions)?;
                        if ret < 0 {
                            if !(CT_ERR_MIN..=CT_ERR_MAX).contains(&ret) {
                                ret = CT_CANON_ERR;
                            }
                            break 'source;
                        }
                    } else {
                        ct_temp = alloc_or_exit!(ConTable, 1, 104, 'source);
                        let mut temporary = ConTable::default();
                        if CTableCreate(heap, &mut temporary, num_atoms, &canon_data[i_base])? == 0
                            || !con_table_is_usable(&temporary, &canon_data[i_base])
                        {
                            ret = CT_OUT_OF_RAM;
                            heap.slice_mut(ct_temp)?[0] = temporary;
                            break 'source;
                        }
                        let base = source_clone(heap, ct_base, 0)?;
                        CtFullCopy(heap, &mut temporary, &base)?;
                        if !iso_sort_key_base.is_null() {
                            for index in 0..total_count {
                                let rank = source_get(heap, n_canon_rank_base, index as i32)?;
                                source_set(
                                    heap,
                                    temporary.iso_sort_key,
                                    i32::from(rank) - 1,
                                    source_get(heap, iso_sort_key_base, index as i32)?,
                                )?;
                            }
                            temporary.len_iso_sort_key = num_at_tg;
                        } else {
                            temporary.len_iso_sort_key = 0;
                        }
                        if !iso_exchg_atnos.is_null() {
                            for index in 0..atom_count {
                                let rank = source_get(heap, n_canon_rank_base, index as i32)?;
                                source_set(
                                    heap,
                                    temporary.iso_exchg_atnos,
                                    i32::from(rank) - 1,
                                    source_get(heap, iso_exchg_atnos, index as i32)?,
                                )?;
                            }
                            temporary.len_iso_exchg_atnos = num_at_tg;
                        } else {
                            temporary.len_iso_exchg_atnos = 0;
                        }
                        heap.slice_mut(ct_temp)?[0] = temporary;
                        ct_base_iso = ct_temp;
                        ct_temp = SourceMutPointer::null();
                        for index in 0..total_count {
                            source_set(
                                heap,
                                n_symm_rank_base_iso,
                                index as i32,
                                source_get(heap, n_symm_rank_base, index as i32)?,
                            )?;
                            source_set(
                                heap,
                                n_canon_rank_base_iso,
                                index as i32,
                                source_get(heap, n_canon_rank_base, index as i32)?,
                            )?;
                            source_set(
                                heap,
                                n_atom_number_canon_base_iso,
                                index as i32,
                                source_get(heap, n_atom_number_canon_base, index as i32)?,
                            )?;
                        }
                    }
                    t_group_info
                        .as_deref_mut()
                        .ok_or(SourceHeapError::NullPointer)?
                        .bIgnoreIsotopic = 1;
                }
            }
        }

        if b_req_taut != 0 && b_req_non_taut != 0 {
            max_len_num_h_fixed = num_atoms.wrapping_add(1);
            n_symm_rank_fix_h = alloc_or_exit!(AT_RANK, num_max, 2, 'source);
            n_canon_rank_fix_h = alloc_or_exit!(AT_RANK, num_max, 2, 'source);
            n_atom_number_canon_fix_h = alloc_or_exit!(AT_NUMB, num_max, 2, 'source);
            num_h_fixed = alloc_or_exit!(NUM_H, max_len_num_h_fixed, 2, 'source);
            let at_other = at[i_other];
            for index in 0..atom_count {
                let base_atom = heap.slice(at_base.as_const())?[index].clone();
                let other_atom = heap.slice(at_other.as_const())?[index].clone();
                let value = if base_atom.endpoint != 0 {
                    if other_atom.num_H != 0 {
                        i32::from(other_atom.num_H).wrapping_add(BASE_H_NUMBER as i32) as NUM_H
                    } else {
                        EMPTY_H_NUMBER as NUM_H
                    }
                } else if other_atom.num_H != base_atom.num_H {
                    i32::from(other_atom.num_H)
                        .wrapping_sub(i32::from(base_atom.num_H))
                        .wrapping_add(BASE_H_NUMBER as i32) as NUM_H
                } else {
                    EMPTY_H_NUMBER as NUM_H
                };
                source_set(heap, num_h_fixed, index as i32, value)?;
            }
            {
                let data = &mut canon_data[i_other];
                data.LinearCT = SourceMutPointer::null();
                data.NumH = num_h_no_taut_h;
                data.NumHfixed = num_h_fixed;
                data.iso_sort_key = SourceMutPointer::null();
                data.iso_exchg_atnos = SourceMutPointer::null();
                data.ulTimeOutTime = pBCN.ulTimeOutTime;
                data.nMaxLenLinearCT = sizes[i_other].nLenCTAtOnly.wrapping_add(1);
                data.maxlenNumH = max_len_num_h_no_taut_h;
                data.maxlenNumHfixed = max_len_num_h_fixed;
                data.nLenLinearCT = sizes[i_other].nLenCTAtOnly;
                data.nLenCTAtOnly = sizes[i_other].nLenCTAtOnly;
                data.lenNumH = num_atoms;
                data.lenNumHfixed = num_atoms;
                data.len_iso_sort_key = 0;
                data.maxlen_iso_sort_key = 0;
                data.len_iso_exchg_atnos = 0;
                data.maxlen_iso_exchg_atnos = 0;
                data.nAuxRank = SourceMutPointer::null();
            }
            for value in heap
                .slice_mut(p_atom_invariant_aux)?
                .get_mut(
                    ..usize::try_from(num_max).map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )
                .ok_or(SourceHeapError::PointerOutOfBounds)?
            {
                *value = ATOM_INVARIANT2::default();
            }
            for index in 0..atom_count {
                let symmetry = source_get(heap, n_symm_rank_base, index as i32)?;
                let hydrogen = source_get(heap, num_h_fixed, index as i32)? as AT_NUMB;
                let invariant = &mut heap.slice_mut(p_atom_invariant_aux)?[index];
                invariant.val[0] = symmetry;
                invariant.val[1] = hydrogen;
            }
            n_num_curr_ranks = SetInitialRanks2(
                heap,
                num_atoms,
                p_atom_invariant_aux,
                n_rank_aux,
                n_atom_number_aux,
                pCG,
            )?;
            n_num_curr_ranks = DifferentiateRanks2(
                heap,
                pCG,
                num_atoms,
                neigh_list[taut_non],
                n_num_curr_ranks,
                n_rank_aux,
                n_temp_rank,
                n_atom_number_aux,
                &mut l_count,
                0,
            )?;
            canon_data[i_other].nAuxRank = n_rank_aux;
            ret = CanonGraph(
                heap,
                ic,
                pCG,
                num_atoms,
                num_atoms,
                num_max,
                0,
                neigh_list[taut_non],
                &mut partitions,
                n_symm_rank_fix_h,
                n_canon_rank_fix_h,
                n_atom_number_canon_fix_h,
                &mut canon_data[i_other],
                &mut canon_counts,
                ct_no_taut_h,
                Some(&mut ct_fix_h),
                LargeMolecules,
                user_action,
                console_quit,
                clock_result,
            )?;
            sync_partition_stack(heap, pBCN.pRankStack, &partitions)?;
            if ret < 0 {
                if !(CT_ERR_MIN..=CT_ERR_MAX).contains(&ret) {
                    ret = CT_CANON_ERR;
                }
                break 'source;
            }

            let mut isotopic_fixed_h = (sizes[i_base].num_isotopic_atoms != 0
                && sizes[i_base].bIgnoreIsotopic == 0)
                || (sizes[i_base].bHasIsotopicTautGroups != 0 && b_taut_ignore_isotopic == 0);
            if bFixIsoFixedH != 0 {
                isotopic_fixed_h |=
                    sizes[i_other].num_isotopic_atoms != 0 && sizes[i_other].bIgnoreIsotopic == 0;
            }
            if isotopic_fixed_h {
                max_len_iso_sort_key_no_taut_h = num_atoms.wrapping_add(1);
                n_symm_rank_fix_h_iso = alloc_or_exit!(AT_RANK, num_max, 2, 'source);
                n_canon_rank_fix_h_iso = alloc_or_exit!(AT_RANK, num_max, 2, 'source);
                n_atom_number_canon_fix_h_iso = alloc_or_exit!(AT_NUMB, num_max, 2, 'source);
                iso_sort_key_no_taut_h = alloc_or_exit!(
                    AT_ISO_SORT_KEY,
                    max_len_iso_sort_key_no_taut_h,
                    8,
                    'source
                );
                let mut num_iso_no_taut_h = 0_i32;
                for index in 0..atom_count {
                    let atom = if bFixIsoFixedH != 0 {
                        heap.slice(at_other.as_const())?[index].clone()
                    } else {
                        heap.slice(at_base.as_const())?[index].clone()
                    };
                    let key = make_iso_sort_key(
                        i32::from(atom.iso_atw_diff),
                        i32::from(atom.num_iso_H[0]),
                        i32::from(atom.num_iso_H[1]),
                        i32::from(atom.num_iso_H[2]),
                    );
                    if key != 0 {
                        num_iso_no_taut_h = num_iso_no_taut_h.wrapping_add(1);
                        source_set(heap, iso_sort_key_no_taut_h, index as i32, key)?;
                    } else {
                        source_set(
                            heap,
                            iso_sort_key_no_taut_h,
                            index as i32,
                            AT_ISO_SORT_KEY::MAX,
                        )?;
                    }
                }
                len_iso_sort_key_no_taut_h = num_atoms;
                {
                    let data = &mut canon_data[i_other];
                    data.LinearCT = SourceMutPointer::null();
                    data.NumH = num_h_no_taut_h;
                    data.NumHfixed = num_h_fixed;
                    data.iso_sort_key = iso_sort_key_no_taut_h;
                    data.iso_exchg_atnos = SourceMutPointer::null();
                    data.ulTimeOutTime = pBCN.ulTimeOutTime;
                    data.nMaxLenLinearCT = sizes[i_other].nLenCTAtOnly.wrapping_add(1);
                    data.maxlenNumH = max_len_num_h_no_taut_h;
                    data.maxlenNumHfixed = max_len_num_h_fixed;
                    data.nLenLinearCT = sizes[i_other].nLenCTAtOnly;
                    data.nLenCTAtOnly = sizes[i_other].nLenCTAtOnly;
                    data.lenNumH = num_atoms;
                    data.lenNumHfixed = num_atoms;
                    data.len_iso_sort_key = len_iso_sort_key_no_taut_h;
                    data.maxlen_iso_sort_key = max_len_iso_sort_key_no_taut_h;
                    data.len_iso_exchg_atnos = 0;
                    data.maxlen_iso_exchg_atnos = 0;
                    data.nAuxRank = SourceMutPointer::null();
                }
                let mut needs_iso_canon = false;
                if num_iso_no_taut_h != 0 {
                    for value in heap
                        .slice_mut(n_temp_rank)?
                        .get_mut(
                            ..usize::try_from(num_max)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                    {
                        *value = 0;
                    }
                    for index in 0..atom_count {
                        let rank = source_get(heap, n_symm_rank_fix_h, index as i32)?;
                        if source_get(heap, n_temp_rank, i32::from(rank) - 1)? < index as AT_RANK {
                            source_set(heap, n_temp_rank, i32::from(rank) - 1, index as AT_RANK)?;
                        }
                    }
                    for index in 0..atom_count {
                        let rank = source_get(heap, n_symm_rank_fix_h, index as i32)?;
                        let representative = source_get(heap, n_temp_rank, i32::from(rank) - 1)?;
                        if source_get(heap, iso_sort_key_no_taut_h, index as i32)?
                            != source_get(heap, iso_sort_key_no_taut_h, i32::from(representative))?
                        {
                            needs_iso_canon = true;
                            break;
                        }
                    }
                }
                if needs_iso_canon {
                    canon_data[i_other].nCanonFlags |= CANON_FLAG_ISO_FIXED_H_DIFF as i32;
                    for value in heap
                        .slice_mut(p_atom_invariant_aux)?
                        .get_mut(
                            ..usize::try_from(num_max)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                    {
                        *value = ATOM_INVARIANT2::default();
                    }
                    for index in 0..atom_count {
                        let symmetry = source_get(heap, n_symm_rank_fix_h, index as i32)?;
                        let isotope = source_get(heap, iso_sort_key_no_taut_h, index as i32)?;
                        let invariant = &mut heap.slice_mut(p_atom_invariant_aux)?[index];
                        invariant.val[0] = symmetry;
                        invariant.iso_sort_key = isotope;
                    }
                    n_num_curr_ranks = SetInitialRanks2(
                        heap,
                        num_atoms,
                        p_atom_invariant_aux,
                        n_rank_aux,
                        n_atom_number_aux,
                        pCG,
                    )?;
                    n_num_curr_ranks = DifferentiateRanks2(
                        heap,
                        pCG,
                        num_atoms,
                        neigh_list[taut_non],
                        n_num_curr_ranks,
                        n_rank_aux,
                        n_temp_rank,
                        n_atom_number_aux,
                        &mut l_count,
                        0,
                    )?;
                    canon_data[i_other].nAuxRank = n_rank_aux;
                    ret = CanonGraph(
                        heap,
                        ic,
                        pCG,
                        num_atoms,
                        num_atoms,
                        num_max,
                        1,
                        neigh_list[taut_non],
                        &mut partitions,
                        n_symm_rank_fix_h_iso,
                        n_canon_rank_fix_h_iso,
                        n_atom_number_canon_fix_h_iso,
                        &mut canon_data[i_other],
                        &mut canon_counts,
                        ct_fix_h,
                        Some(&mut ct_fix_h_iso),
                        LargeMolecules,
                        user_action,
                        console_quit,
                        clock_result,
                    )?;
                    sync_partition_stack(heap, pBCN.pRankStack, &partitions)?;
                    if ret < 0 {
                        if !(CT_ERR_MIN..=CT_ERR_MAX).contains(&ret) {
                            ret = CT_CANON_ERR;
                        }
                        break 'source;
                    }
                } else {
                    ct_temp = alloc_or_exit!(ConTable, 1, 104, 'source);
                    let mut temporary = ConTable::default();
                    if CTableCreate(heap, &mut temporary, num_atoms, &canon_data[i_other])? == 0
                        || !con_table_is_usable(&temporary, &canon_data[i_other])
                    {
                        ret = CT_OUT_OF_RAM;
                        heap.slice_mut(ct_temp)?[0] = temporary;
                        break 'source;
                    }
                    let fixed = source_clone(heap, ct_fix_h, 0)?;
                    CtFullCopy(heap, &mut temporary, &fixed)?;
                    for index in 0..atom_count {
                        let rank = source_get(heap, n_canon_rank_fix_h, index as i32)?;
                        source_set(
                            heap,
                            temporary.iso_sort_key,
                            i32::from(rank) - 1,
                            source_get(heap, iso_sort_key_no_taut_h, index as i32)?,
                        )?;
                    }
                    temporary.len_iso_sort_key = num_atoms;
                    heap.slice_mut(ct_temp)?[0] = temporary;
                    ct_fix_h_iso = ct_temp;
                    ct_temp = SourceMutPointer::null();
                    for index in 0..atom_count {
                        source_set(
                            heap,
                            n_symm_rank_fix_h_iso,
                            index as i32,
                            source_get(heap, n_symm_rank_fix_h, index as i32)?,
                        )?;
                        source_set(
                            heap,
                            n_canon_rank_fix_h_iso,
                            index as i32,
                            source_get(heap, n_canon_rank_fix_h, index as i32)?,
                        )?;
                        source_set(
                            heap,
                            n_atom_number_canon_fix_h_iso,
                            index as i32,
                            source_get(heap, n_atom_number_canon_fix_h, index as i32)?,
                        )?;
                    }
                }
            }
        }

        ret = 0;
        let no_h = source_clone(heap, ct_no_h, 0)?;
        let no_taut_h = source_clone(heap, ct_no_taut_h, 0)?;
        let mut inconsistent = no_h.lenCt != no_taut_h.lenCt
            || !source_prefix_equal(heap, no_h.Ctbl, no_taut_h.Ctbl, no_h.lenCt)?;
        if b_req_taut != 0 {
            if !ct_fix_h.is_null() {
                let fixed = source_clone(heap, ct_fix_h, 0)?;
                inconsistent |= no_taut_h.lenCt != fixed.lenCt
                    || !source_prefix_equal(heap, no_taut_h.Ctbl, fixed.Ctbl, no_taut_h.lenCt)?;
                inconsistent |= no_taut_h.lenNumH != fixed.lenNumH
                    || !source_prefix_equal(heap, no_taut_h.NumH, fixed.NumH, no_taut_h.lenNumH)?;
            }
            let base = source_clone(heap, ct_base, 0)?;
            inconsistent |= no_taut_h.lenCt > base.lenCt
                || !source_prefix_equal(heap, no_taut_h.Ctbl, base.Ctbl, no_taut_h.lenCt)?;
            inconsistent |= no_taut_h.lenNumH > base.lenNumH
                || !source_prefix_equal(heap, no_taut_h.NumH, base.NumH, no_taut_h.lenNumH)?;
        }
        if !ct_no_taut_h_iso.is_null() {
            let isotopic = source_clone(heap, ct_no_taut_h_iso, 0)?;
            inconsistent |= no_h.lenCt != isotopic.lenCt
                || !source_prefix_equal(heap, no_h.Ctbl, isotopic.Ctbl, no_h.lenCt)?;
            inconsistent |= no_taut_h.lenNumH != isotopic.lenNumH
                || !source_prefix_equal(heap, no_taut_h.NumH, isotopic.NumH, no_taut_h.lenNumH)?;
        } else if !ct_base_iso.is_null() {
            let base = source_clone(heap, ct_base, 0)?;
            let base_iso = source_clone(heap, ct_base_iso, 0)?;
            inconsistent |= base_iso.lenCt != base.lenCt
                || !source_prefix_equal(heap, base_iso.Ctbl, base.Ctbl, base_iso.lenCt)?;
            inconsistent |= base_iso.lenNumH != base.lenNumH
                || !source_prefix_equal(heap, base_iso.NumH, base.NumH, base_iso.lenNumH)?;
            if !ct_fix_h_iso.is_null() {
                let fixed_iso = source_clone(heap, ct_fix_h_iso, 0)?;
                inconsistent |= fixed_iso.lenCt > base_iso.lenCt
                    || !source_prefix_equal(heap, fixed_iso.Ctbl, base_iso.Ctbl, fixed_iso.lenCt)?;
                inconsistent |= fixed_iso.lenNumH > base_iso.lenNumH
                    || !source_prefix_equal(
                        heap,
                        fixed_iso.NumH,
                        base_iso.NumH,
                        fixed_iso.lenNumH,
                    )?;
            }
        }
        if inconsistent {
            ret = CT_CANON_ERR;
            break 'source;
        }

        if b_req_taut != 0 {
            t_group_info
                .as_deref_mut()
                .ok_or(SourceHeapError::NullPointer)?
                .bIgnoreIsotopic = b_taut_ignore_isotopic;
        }
        pBCN.num_max = num_max;
        pBCN.num_at_tg = num_at_tg;
        pBCN.num_atoms = num_atoms;
        pBCN.ftcn[taut_non].NeighList = neigh_list[taut_non];
        neigh_list[taut_non] = SourceMutPointer::null();
        pBCN.ftcn[taut_yes].NeighList = neigh_list[taut_yes];
        neigh_list[taut_yes] = SourceMutPointer::null();

        if b_req_taut != 0 {
            let base = heap
                .slice_mut(ct_base)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let output = &mut pBCN.ftcn[taut_yes];
            output.num_at_tg = num_at_tg;
            output.num_atoms = num_atoms;
            output.LinearCt = base.Ctbl;
            base.Ctbl = SourceMutPointer::null();
            output.nLenLinearCtAtOnly = sizes[i_base].nLenCTAtOnly;
            output.nMaxLenLinearCt = sizes[i_base].nLenCT.wrapping_add(1);
            output.nLenLinearCt = sizes[i_base].nLenCT;
            output.PartitionCt.Rank = n_canon_rank_base;
            n_canon_rank_base = SourceMutPointer::null();
            output.PartitionCt.AtNumber = n_atom_number_canon_base;
            n_atom_number_canon_base = SourceMutPointer::null();
            output.nSymmRankCt = n_symm_rank_base;
            n_symm_rank_base = SourceMutPointer::null();
            output.nNumHOrig = num_h;
            num_h = SourceMutPointer::null();
            output.nNumH = base.NumH;
            base.NumH = SourceMutPointer::null();
            output.nLenNumH = max_len_num_h.min(base.maxlenNumH);
            output.nNumHOrigFixH = SourceMutPointer::null();
            output.nNumHFixH = SourceMutPointer::null();
            output.nLenNumHFixH = 0;
            output.nCanonFlags |= canon_data[i_base].nCanonFlags;
            CleanNumH(heap, output.nNumHOrig, output.nLenNumH)?;
            CleanNumH(heap, output.nNumH, output.nLenNumH)?;
            CleanCt(heap, output.LinearCt, output.nLenLinearCt)?;
            if !ct_base_iso.is_null() {
                let base_iso = heap
                    .slice_mut(ct_base_iso)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                output.PartitionCtIso.Rank = n_canon_rank_base_iso;
                n_canon_rank_base_iso = SourceMutPointer::null();
                output.PartitionCtIso.AtNumber = n_atom_number_canon_base_iso;
                n_atom_number_canon_base_iso = SourceMutPointer::null();
                output.nSymmRankCtIso = n_symm_rank_base_iso;
                n_symm_rank_base_iso = SourceMutPointer::null();
                output.iso_sort_keys = base_iso.iso_sort_key;
                base_iso.iso_sort_key = SourceMutPointer::null();
                output.iso_sort_keysOrig = iso_sort_key_base;
                iso_sort_key_base = SourceMutPointer::null();
                output.len_iso_sort_keys = len_iso_sort_key_base;
                output.iso_exchg_atnos = base_iso.iso_exchg_atnos;
                base_iso.iso_exchg_atnos = SourceMutPointer::null();
                output.iso_exchg_atnosOrig = iso_exchg_atnos;
                iso_exchg_atnos = SourceMutPointer::null();
                CleanIsoSortKeys(heap, output.iso_sort_keys, output.len_iso_sort_keys)?;
                CleanIsoSortKeys(heap, output.iso_sort_keysOrig, output.len_iso_sort_keys)?;
            }
        }

        if b_req_non_taut != 0 {
            if b_req_taut == 0 {
                n_symm_rank_fix_h = n_symm_rank_no_taut_h;
                n_symm_rank_no_taut_h = SourceMutPointer::null();
                n_canon_rank_fix_h = n_canon_rank_no_taut_h;
                n_canon_rank_no_taut_h = SourceMutPointer::null();
                n_atom_number_canon_fix_h = n_atom_number_canon_no_taut_h;
                n_atom_number_canon_no_taut_h = SourceMutPointer::null();
                ct_fix_h = ct_no_taut_h;
                ct_no_taut_h = SourceMutPointer::null();
                n_symm_rank_fix_h_iso = n_symm_rank_no_taut_h_iso;
                n_symm_rank_no_taut_h_iso = SourceMutPointer::null();
                n_canon_rank_fix_h_iso = n_canon_rank_no_taut_h_iso;
                n_canon_rank_no_taut_h_iso = SourceMutPointer::null();
                n_atom_number_canon_fix_h_iso = n_atom_number_canon_no_taut_h_iso;
                n_atom_number_canon_no_taut_h_iso = SourceMutPointer::null();
                ct_fix_h_iso = ct_no_taut_h_iso;
                ct_no_taut_h_iso = SourceMutPointer::null();
                if i_other == taut_yes
                    && !pBCN.ftcn[taut_non].NeighList.is_null()
                    && pBCN.ftcn[taut_yes].NeighList.is_null()
                {
                    pBCN.ftcn[taut_yes].NeighList = pBCN.ftcn[taut_non].NeighList;
                    pBCN.ftcn[taut_non].NeighList = SourceMutPointer::null();
                }
            }
            let fixed = heap
                .slice_mut(ct_fix_h)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let output = &mut pBCN.ftcn[i_other];
            output.num_at_tg = num_atoms;
            output.num_atoms = num_atoms;
            output.LinearCt = fixed.Ctbl;
            fixed.Ctbl = SourceMutPointer::null();
            output.nLenLinearCtAtOnly = sizes[i_other].nLenCTAtOnly;
            output.nMaxLenLinearCt = sizes[i_other].nLenCTAtOnly.wrapping_add(1);
            output.nLenLinearCt = sizes[i_other].nLenCTAtOnly;
            output.PartitionCt.Rank = n_canon_rank_fix_h;
            n_canon_rank_fix_h = SourceMutPointer::null();
            output.PartitionCt.AtNumber = n_atom_number_canon_fix_h;
            n_atom_number_canon_fix_h = SourceMutPointer::null();
            output.nSymmRankCt = n_symm_rank_fix_h;
            n_symm_rank_fix_h = SourceMutPointer::null();
            output.nNumHOrig = num_h_no_taut_h;
            num_h_no_taut_h = SourceMutPointer::null();
            output.nNumH = fixed.NumH;
            fixed.NumH = SourceMutPointer::null();
            output.nLenNumH = max_len_num_h_no_taut_h.min(fixed.maxlenNumH);
            output.nNumHOrigFixH = num_h_fixed;
            num_h_fixed = SourceMutPointer::null();
            output.nNumHFixH = fixed.NumHfixed;
            fixed.NumHfixed = SourceMutPointer::null();
            output.nLenNumHFixH = num_atoms;
            output.nCanonFlags |= canon_data[i_other].nCanonFlags;
            CleanNumH(heap, output.nNumHOrig, output.nLenNumH)?;
            CleanNumH(heap, output.nNumHOrigFixH, output.nLenNumH)?;
            CleanNumH(heap, output.nNumH, output.nLenNumH)?;
            CleanNumH(heap, output.nNumHFixH, output.nLenNumH)?;
            CleanCt(heap, output.LinearCt, output.nLenLinearCt)?;
            if !ct_fix_h_iso.is_null() {
                let fixed_iso = heap
                    .slice_mut(ct_fix_h_iso)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                output.PartitionCtIso.Rank = n_canon_rank_fix_h_iso;
                n_canon_rank_fix_h_iso = SourceMutPointer::null();
                output.PartitionCtIso.AtNumber = n_atom_number_canon_fix_h_iso;
                n_atom_number_canon_fix_h_iso = SourceMutPointer::null();
                output.nSymmRankCtIso = n_symm_rank_fix_h_iso;
                n_symm_rank_fix_h_iso = SourceMutPointer::null();
                output.iso_sort_keys = fixed_iso.iso_sort_key;
                fixed_iso.iso_sort_key = SourceMutPointer::null();
                output.iso_sort_keysOrig = iso_sort_key_no_taut_h;
                iso_sort_key_no_taut_h = SourceMutPointer::null();
                output.len_iso_sort_keys = len_iso_sort_key_no_taut_h;
                CleanIsoSortKeys(heap, output.iso_sort_keys, output.len_iso_sort_keys)?;
                CleanIsoSortKeys(heap, output.iso_sort_keysOrig, output.len_iso_sort_keys)?;
            }
        }
    }

    FreeNeighList(heap, neigh_list[TAUT_NON as usize])?;
    FreeNeighList(heap, neigh_list[TAUT_YES as usize])?;
    for table in [
        ct_no_h,
        ct_no_taut_h,
        ct_base,
        ct_fix_h,
        ct_temp,
        ct_no_taut_h_iso,
        ct_base_iso,
        ct_fix_h_iso,
    ] {
        free_con_table(heap, table)?;
    }
    for pointer in [n_rank, n_atom_number] {
        if !pointer.is_null() {
            inchi_free(heap, pointer)?;
        }
    }
    if !pBCN.pRankStack.is_null() {
        source_set(heap, pBCN.pRankStack, 0, SourceMutPointer::null())?;
        source_set(heap, pBCN.pRankStack, 1, SourceMutPointer::null())?;
    }
    for pointer in [n_rank_aux, n_atom_number_aux] {
        if !pointer.is_null() {
            inchi_free(heap, pointer)?;
        }
    }
    if !p_atom_invariant_aux.is_null() {
        inchi_free(heap, p_atom_invariant_aux)?;
    }
    if !p_atom_invariant.is_null() {
        inchi_free(heap, p_atom_invariant)?;
    }
    for pointer in [
        n_canon_rank_no_h,
        n_atom_number_canon_no_h,
        n_symm_rank_no_h,
        n_symm_rank_no_taut_h,
        n_canon_rank_no_taut_h,
        n_atom_number_canon_no_taut_h,
        n_symm_rank_base,
        n_canon_rank_base,
        n_atom_number_canon_base,
        n_symm_rank_fix_h,
        n_canon_rank_fix_h,
        n_atom_number_canon_fix_h,
        n_symm_rank_no_taut_h_iso,
        n_canon_rank_no_taut_h_iso,
        n_atom_number_canon_no_taut_h_iso,
        n_symm_rank_base_iso,
        n_canon_rank_base_iso,
        n_atom_number_canon_base_iso,
        n_symm_rank_fix_h_iso,
        n_canon_rank_fix_h_iso,
        n_atom_number_canon_fix_h_iso,
        n_temp_rank,
    ] {
        if !pointer.is_null() {
            inchi_free(heap, pointer)?;
        }
    }
    for pointer in [num_h_no_taut_h, num_h, num_h_fixed] {
        if !pointer.is_null() {
            inchi_free(heap, pointer)?;
        }
    }
    for pointer in [iso_sort_key_no_taut_h, iso_sort_key_base] {
        if !pointer.is_null() {
            inchi_free(heap, pointer)?;
        }
    }
    if !iso_exchg_atnos.is_null() {
        inchi_free(heap, iso_exchg_atnos)?;
    }

    Ok(ret)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_port__ichican2__getbasecanonranking__line_5105() {
        fn atom(symbol: &[u8], isotope_difference: i8) -> sp_ATOM {
            let mut atom = sp_ATOM {
                iso_atw_diff: isotope_difference,
                ..sp_ATOM::default()
            };
            for (target, source) in atom.elname.iter_mut().zip(symbol.iter().copied()) {
                *target = source as i8;
            }
            atom
        }

        fn run_taut_aux_allocation_case(
            failure_ordinal: Option<u64>,
        ) -> (Result<i32, SourceHeapError>, u64, usize) {
            let mut heap = SourceHeap::default();
            let mut endpoint = atom(b"N", 0);
            endpoint.endpoint = 1;
            let mut exchange = atom(b"N", 0);
            exchange.cFlags = AT_FLAG_ISO_H_POINT as i8;
            let atoms = heap
                .allocate_model_storage(vec![endpoint, exchange])
                .unwrap();
            let groups = heap
                .allocate_model_storage(vec![T_GROUP {
                    num: [1, 0, 0, 0, 0],
                    nGroupNumber: 1,
                    nNumEndpoints: 1,
                    nFirstEndpointAtNoPos: 0,
                    ..T_GROUP::default()
                }])
                .unwrap();
            let endpoints = heap.allocate_model_storage(vec![0_u16]).unwrap();
            let isotopic_endpoints = heap.allocate_model_storage(vec![1_u16, 1]).unwrap();
            let mut info = T_GROUP_INFO {
                t_group: groups,
                nEndpointAtomNumber: endpoints,
                nIsotopicEndpointAtomNumber: isotopic_endpoints,
                nNumIsotopicEndpoints: 2,
                num_t_groups: 1,
                max_num_t_groups: 1,
                bTautFlagsDone: u64::from(TG_FLAG_FOUND_ISOTOPIC_H_DONE),
                ..T_GROUP_INFO::default()
            };
            let mut sizes = std::array::from_fn(|_| ATOM_SIZES::default());
            sizes[TAUT_YES as usize] = ATOM_SIZES {
                nLenCT: 4,
                nLenCTAtOnly: 2,
                nLenLinearCTTautomer: 1,
                nLenIsotopicEndpoints: 2,
                ..ATOM_SIZES::default()
            };
            let mut bcn = BCN::default();
            let mut globals = CANON_GLOBALS::default();
            if let Some(ordinal) = failure_ordinal {
                heap.fail_after_allocations(ordinal);
            }
            let result = GetBaseCanonRanking(
                &mut heap,
                &mut INCHI_CLOCK::default(),
                2,
                3,
                [SourceMutPointer::null(), atoms],
                Some(&mut info),
                &sizes,
                &mut bcn,
                SourceMutPointer::null(),
                &mut globals,
                1,
                0,
                None,
                None,
                0,
            );
            let allocation_calls = heap.source_allocation_calls();
            DeAllocBCN(&mut heap, Some(&mut bcn)).unwrap();
            SetBitFree(&mut heap, &mut globals).unwrap();
            inchi_free(&mut heap, atoms).unwrap();
            inchi_free(&mut heap, groups).unwrap();
            inchi_free(&mut heap, endpoints).unwrap();
            inchi_free(&mut heap, isotopic_endpoints).unwrap();
            (result, allocation_calls, heap.live_allocation_count())
        }

        fn run_fixed_iso_allocation_case(
            failure_ordinal: Option<u64>,
        ) -> (Result<i32, SourceHeapError>, u64, usize) {
            let mut heap = SourceHeap::default();
            let mut base_first = atom(b"N", 0);
            base_first.endpoint = 1;
            base_first.num_H = 1;
            let base_second = base_first.clone();
            let mut other_first = atom(b"N", 0);
            other_first.num_H = 1;
            let mut other_second = other_first.clone();
            other_second.iso_atw_diff = 1;
            let base_atoms = heap
                .allocate_model_storage(vec![base_first, base_second])
                .unwrap();
            let other_atoms = heap
                .allocate_model_storage(vec![other_first, other_second])
                .unwrap();
            let groups = heap
                .allocate_model_storage(vec![T_GROUP {
                    num: [2, 0, 0, 0, 0],
                    nGroupNumber: 1,
                    nNumEndpoints: 2,
                    nFirstEndpointAtNoPos: 0,
                    ..T_GROUP::default()
                }])
                .unwrap();
            let endpoints = heap.allocate_model_storage(vec![0_u16, 1]).unwrap();
            let mut info = T_GROUP_INFO {
                t_group: groups,
                nEndpointAtomNumber: endpoints,
                num_t_groups: 1,
                max_num_t_groups: 1,
                ..T_GROUP_INFO::default()
            };
            let mut sizes = std::array::from_fn(|_| ATOM_SIZES::default());
            sizes[TAUT_YES as usize] = ATOM_SIZES {
                nLenCT: 5,
                nLenCTAtOnly: 2,
                nLenLinearCTTautomer: 1,
                ..ATOM_SIZES::default()
            };
            sizes[TAUT_NON as usize] = ATOM_SIZES {
                num_isotopic_atoms: 1,
                nLenCT: 2,
                nLenCTAtOnly: 2,
                ..ATOM_SIZES::default()
            };
            let mut bcn = BCN::default();
            let mut globals = CANON_GLOBALS::default();
            if let Some(ordinal) = failure_ordinal {
                heap.fail_after_allocations(ordinal);
            }
            let result = GetBaseCanonRanking(
                &mut heap,
                &mut INCHI_CLOCK::default(),
                2,
                3,
                [other_atoms, base_atoms],
                Some(&mut info),
                &sizes,
                &mut bcn,
                SourceMutPointer::null(),
                &mut globals,
                1,
                0,
                None,
                None,
                0,
            );
            let allocation_calls = heap.source_allocation_calls();
            DeAllocBCN(&mut heap, Some(&mut bcn)).unwrap();
            SetBitFree(&mut heap, &mut globals).unwrap();
            inchi_free(&mut heap, base_atoms).unwrap();
            inchi_free(&mut heap, other_atoms).unwrap();
            inchi_free(&mut heap, groups).unwrap();
            inchi_free(&mut heap, endpoints).unwrap();
            (result, allocation_calls, heap.live_allocation_count())
        }

        fn run_single_atom(isotopic: bool) -> (SourceHeap, BCN, SourceMutPointer<sp_ATOM>, i32) {
            let mut heap = SourceHeap::default();
            let atoms = heap
                .allocate_model_storage(vec![atom(b"C", i8::from(isotopic))])
                .unwrap();
            let mut sizes = std::array::from_fn(|_| ATOM_SIZES::default());
            sizes[TAUT_NON as usize] = ATOM_SIZES {
                num_isotopic_atoms: i32::from(isotopic),
                nLenCT: 1,
                nLenCTAtOnly: 1,
                ..ATOM_SIZES::default()
            };
            let mut bcn = BCN::default();
            let mut globals = CANON_GLOBALS::default();
            let result = GetBaseCanonRanking(
                &mut heap,
                &mut INCHI_CLOCK::default(),
                1,
                1,
                [atoms, SourceMutPointer::null()],
                None,
                &sizes,
                &mut bcn,
                SourceMutPointer::null(),
                &mut globals,
                1,
                0,
                None,
                None,
                0,
            )
            .unwrap();
            SetBitFree(&mut heap, &mut globals).unwrap();
            (heap, bcn, atoms, result)
        }

        let mut empty_heap = SourceHeap::default();
        let mut empty_bcn = BCN::default();
        assert_eq!(
            GetBaseCanonRanking(
                &mut empty_heap,
                &mut INCHI_CLOCK::default(),
                0,
                0,
                [SourceMutPointer::null(); TAUT_NUM as usize],
                None,
                &std::array::from_fn(|_| ATOM_SIZES::default()),
                &mut empty_bcn,
                SourceMutPointer::null(),
                &mut CANON_GLOBALS::default(),
                1,
                0,
                None,
                None,
                0,
            ),
            Ok(CT_UNKNOWN_ERR)
        );
        assert_eq!(empty_heap.live_allocation_count(), 0);

        let mut non_isotopic_allocation_calls = 0_u64;
        for isotopic in [false, true] {
            let (mut heap, mut bcn, atoms, result) = run_single_atom(isotopic);
            assert_eq!(result, 0, "isotopic={isotopic}");
            assert_eq!((bcn.num_atoms, bcn.num_at_tg, bcn.num_max), (1, 1, 1));
            let output = &bcn.ftcn[TAUT_NON as usize];
            assert_eq!((output.num_atoms, output.num_at_tg), (1, 1));
            assert_eq!(heap.slice(output.LinearCt.as_const()).unwrap()[0], 1);
            assert_eq!(
                heap.slice(output.PartitionCt.Rank.as_const()).unwrap(),
                &[1]
            );
            assert_eq!(
                heap.slice(output.PartitionCt.AtNumber.as_const()).unwrap(),
                &[0]
            );
            assert_eq!(heap.slice(output.nSymmRankCt.as_const()).unwrap(), &[1]);
            assert_eq!(heap.slice(output.nNumHOrig.as_const()).unwrap()[0], 0);
            assert_eq!(heap.slice(output.nNumH.as_const()).unwrap()[0], 0);
            if isotopic {
                assert_eq!(
                    heap.slice(output.PartitionCtIso.Rank.as_const()).unwrap(),
                    &[1]
                );
                assert_eq!(output.len_iso_sort_keys, 1);
                assert_eq!(
                    heap.slice(output.iso_sort_keys.as_const()).unwrap()[0],
                    32768
                );
                assert_eq!(
                    heap.slice(output.iso_sort_keysOrig.as_const()).unwrap()[0],
                    32768
                );
            } else {
                assert!(output.PartitionCtIso.Rank.is_null());
                assert!(output.iso_sort_keys.is_null());
            }
            DeAllocBCN(&mut heap, Some(&mut bcn)).unwrap();
            inchi_free(&mut heap, atoms).unwrap();
            if !isotopic {
                non_isotopic_allocation_calls = heap.source_allocation_calls();
            }
            assert_eq!(heap.live_allocation_count(), 0, "isotopic={isotopic}");
        }

        assert!(non_isotopic_allocation_calls > 0);
        for allocation_ordinal in 0..non_isotopic_allocation_calls {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(vec![atom(b"C", 0)]).unwrap();
            let mut sizes = std::array::from_fn(|_| ATOM_SIZES::default());
            sizes[TAUT_NON as usize] = ATOM_SIZES {
                nLenCT: 1,
                nLenCTAtOnly: 1,
                ..ATOM_SIZES::default()
            };
            let mut bcn = BCN::default();
            let mut globals = CANON_GLOBALS::default();
            heap.fail_after_allocations(allocation_ordinal);
            let result = GetBaseCanonRanking(
                &mut heap,
                &mut INCHI_CLOCK::default(),
                1,
                1,
                [atoms, SourceMutPointer::null()],
                None,
                &sizes,
                &mut bcn,
                SourceMutPointer::null(),
                &mut globals,
                1,
                0,
                None,
                None,
                0,
            );
            match result {
                Ok(code) if code < 0 => {}
                other => panic!(
                    "allocation ordinal {allocation_ordinal}, calls {}: {other:?}",
                    heap.source_allocation_calls()
                ),
            }
            DeAllocBCN(&mut heap, Some(&mut bcn)).unwrap();
            SetBitFree(&mut heap, &mut globals).unwrap();
            inchi_free(&mut heap, atoms).unwrap();
            assert_eq!(
                heap.live_allocation_count(),
                0,
                "allocation ordinal {allocation_ordinal}"
            );
        }

        let mut differentiated_heap = SourceHeap::default();
        let first = atom(b"C", 0);
        let mut second = atom(b"C", 1);
        second.num_H = 1;
        let differentiated_atoms = differentiated_heap
            .allocate_model_storage(vec![first, second])
            .unwrap();
        let mut differentiated_sizes = std::array::from_fn(|_| ATOM_SIZES::default());
        differentiated_sizes[TAUT_NON as usize] = ATOM_SIZES {
            num_isotopic_atoms: 1,
            nLenCT: 2,
            nLenCTAtOnly: 2,
            ..ATOM_SIZES::default()
        };
        let mut differentiated_bcn = BCN::default();
        let mut differentiated_globals = CANON_GLOBALS::default();
        assert_eq!(
            GetBaseCanonRanking(
                &mut differentiated_heap,
                &mut INCHI_CLOCK::default(),
                2,
                2,
                [differentiated_atoms, SourceMutPointer::null()],
                None,
                &differentiated_sizes,
                &mut differentiated_bcn,
                SourceMutPointer::null(),
                &mut differentiated_globals,
                1,
                0,
                None,
                None,
                0,
            ),
            Ok(0)
        );
        let differentiated_output = &differentiated_bcn.ftcn[TAUT_NON as usize];
        assert_ne!(
            differentiated_output.nCanonFlags & CANON_FLAG_NO_TAUT_H_DIFF as i32,
            0
        );
        assert_eq!(
            differentiated_output.nCanonFlags & CANON_FLAG_ISO_ONLY_NON_TAUT_DIFF as i32,
            0
        );
        assert_eq!(
            differentiated_heap
                .slice(differentiated_output.PartitionCt.Rank.as_const())
                .unwrap(),
            &[2, 1]
        );
        assert_eq!(
            &differentiated_heap
                .slice(differentiated_output.nNumHOrig.as_const())
                .unwrap()[..2],
            &[0, 1]
        );
        assert_eq!(differentiated_output.len_iso_sort_keys, 2);
        assert_eq!(
            &differentiated_heap
                .slice(differentiated_output.iso_sort_keysOrig.as_const())
                .unwrap()[..2],
            &[0, 32768]
        );
        DeAllocBCN(&mut differentiated_heap, Some(&mut differentiated_bcn)).unwrap();
        SetBitFree(&mut differentiated_heap, &mut differentiated_globals).unwrap();
        inchi_free(&mut differentiated_heap, differentiated_atoms).unwrap();
        assert_eq!(differentiated_heap.live_allocation_count(), 0);

        let mut isotope_split_heap = SourceHeap::default();
        let isotope_split_atoms = isotope_split_heap
            .allocate_model_storage(vec![atom(b"C", 0), atom(b"C", 1)])
            .unwrap();
        let mut isotope_split_sizes = std::array::from_fn(|_| ATOM_SIZES::default());
        isotope_split_sizes[TAUT_NON as usize] = ATOM_SIZES {
            num_isotopic_atoms: 1,
            nLenCT: 2,
            nLenCTAtOnly: 2,
            ..ATOM_SIZES::default()
        };
        let mut isotope_split_bcn = BCN::default();
        let mut isotope_split_globals = CANON_GLOBALS::default();
        assert_eq!(
            GetBaseCanonRanking(
                &mut isotope_split_heap,
                &mut INCHI_CLOCK::default(),
                2,
                2,
                [isotope_split_atoms, SourceMutPointer::null()],
                None,
                &isotope_split_sizes,
                &mut isotope_split_bcn,
                SourceMutPointer::null(),
                &mut isotope_split_globals,
                1,
                0,
                None,
                None,
                0,
            ),
            Ok(0)
        );
        let isotope_split_output = &isotope_split_bcn.ftcn[TAUT_NON as usize];
        assert_eq!(
            isotope_split_output.nCanonFlags & CANON_FLAG_NO_TAUT_H_DIFF as i32,
            0
        );
        assert_ne!(
            isotope_split_output.nCanonFlags & CANON_FLAG_ISO_ONLY_NON_TAUT_DIFF as i32,
            0
        );
        assert_eq!(isotope_split_output.len_iso_sort_keys, 2);
        DeAllocBCN(&mut isotope_split_heap, Some(&mut isotope_split_bcn)).unwrap();
        SetBitFree(&mut isotope_split_heap, &mut isotope_split_globals).unwrap();
        inchi_free(&mut isotope_split_heap, isotope_split_atoms).unwrap();
        assert_eq!(isotope_split_heap.live_allocation_count(), 0);

        fn quit_action() -> u32 {
            USER_ACTION_QUIT
        }
        for timed_out in [false, true] {
            let mut error_heap = SourceHeap::default();
            let error_atoms = error_heap
                .allocate_model_storage(vec![atom(b"C", 0), atom(b"C", 0)])
                .unwrap();
            let timeout = if timed_out {
                error_heap
                    .allocate_model_storage(vec![inchiTime { clockTime: 0 }])
                    .unwrap()
            } else {
                SourceMutPointer::null()
            };
            let mut error_sizes = std::array::from_fn(|_| ATOM_SIZES::default());
            error_sizes[TAUT_NON as usize] = ATOM_SIZES {
                nLenCT: 2,
                nLenCTAtOnly: 2,
                ..ATOM_SIZES::default()
            };
            let mut error_bcn = BCN::default();
            let mut error_globals = CANON_GLOBALS::default();
            assert_eq!(
                GetBaseCanonRanking(
                    &mut error_heap,
                    &mut INCHI_CLOCK::default(),
                    2,
                    2,
                    [error_atoms, SourceMutPointer::null()],
                    None,
                    &error_sizes,
                    &mut error_bcn,
                    timeout,
                    &mut error_globals,
                    1,
                    0,
                    if timed_out { None } else { Some(quit_action) },
                    None,
                    i64::from(timed_out),
                ),
                Ok(if timed_out {
                    CT_TIMEOUT_ERR
                } else {
                    CT_USER_QUIT_ERR
                })
            );
            DeAllocBCN(&mut error_heap, Some(&mut error_bcn)).unwrap();
            SetBitFree(&mut error_heap, &mut error_globals).unwrap();
            inchi_free(&mut error_heap, error_atoms).unwrap();
            if !timeout.is_null() {
                inchi_free(&mut error_heap, timeout).unwrap();
            }
            assert_eq!(
                error_heap.live_allocation_count(),
                0,
                "timed_out={timed_out}"
            );
        }

        let mut taut_heap = SourceHeap::default();
        let mut endpoint_atom = atom(b"N", 0);
        endpoint_atom.endpoint = 1;
        endpoint_atom.num_H = 1;
        let taut_atoms = taut_heap
            .allocate_model_storage(vec![endpoint_atom])
            .unwrap();
        let groups = taut_heap
            .allocate_model_storage(vec![T_GROUP {
                num: [1, 0, 0, 0, 0],
                nGroupNumber: 1,
                nNumEndpoints: 1,
                nFirstEndpointAtNoPos: 0,
                ..T_GROUP::default()
            }])
            .unwrap();
        let endpoint_numbers = taut_heap.allocate_model_storage(vec![0_u16]).unwrap();
        let mut taut_info = T_GROUP_INFO {
            t_group: groups,
            nEndpointAtomNumber: endpoint_numbers,
            num_t_groups: 1,
            max_num_t_groups: 1,
            bIgnoreIsotopic: 0,
            ..T_GROUP_INFO::default()
        };
        let mut taut_sizes = std::array::from_fn(|_| ATOM_SIZES::default());
        taut_sizes[TAUT_YES as usize] = ATOM_SIZES {
            nLenCT: 3,
            nLenCTAtOnly: 1,
            nLenLinearCTTautomer: 1,
            ..ATOM_SIZES::default()
        };
        let mut taut_bcn = BCN::default();
        let mut taut_globals = CANON_GLOBALS::default();
        assert_eq!(
            GetBaseCanonRanking(
                &mut taut_heap,
                &mut INCHI_CLOCK::default(),
                1,
                2,
                [SourceMutPointer::null(), taut_atoms],
                Some(&mut taut_info),
                &taut_sizes,
                &mut taut_bcn,
                SourceMutPointer::null(),
                &mut taut_globals,
                1,
                0,
                None,
                None,
                0,
            ),
            Ok(0)
        );
        assert_eq!(taut_info.bIgnoreIsotopic, 0);
        let taut_output = &taut_bcn.ftcn[TAUT_YES as usize];
        assert_eq!((taut_output.num_atoms, taut_output.num_at_tg), (1, 2));
        assert_eq!(taut_output.nLenLinearCt, 3);
        assert_eq!(
            &taut_heap.slice(taut_output.LinearCt.as_const()).unwrap()[..3],
            &[1, 2, 1]
        );
        assert_eq!(
            taut_heap
                .slice(taut_output.PartitionCt.Rank.as_const())
                .unwrap(),
            &[1, 2]
        );
        assert_eq!(
            &taut_heap.slice(taut_output.nNumHOrig.as_const()).unwrap()[..3],
            &[0, 1, 0]
        );
        DeAllocBCN(&mut taut_heap, Some(&mut taut_bcn)).unwrap();
        SetBitFree(&mut taut_heap, &mut taut_globals).unwrap();
        inchi_free(&mut taut_heap, taut_atoms).unwrap();
        inchi_free(&mut taut_heap, groups).unwrap();
        inchi_free(&mut taut_heap, endpoint_numbers).unwrap();
        assert_eq!(taut_heap.live_allocation_count(), 0);

        let mut taut_isotope_heap = SourceHeap::default();
        let mut taut_isotope_first = atom(b"N", 0);
        taut_isotope_first.endpoint = 1;
        let mut taut_isotope_second = atom(b"N", 1);
        taut_isotope_second.endpoint = 1;
        let taut_isotope_atoms = taut_isotope_heap
            .allocate_model_storage(vec![taut_isotope_first, taut_isotope_second])
            .unwrap();
        let taut_isotope_groups = taut_isotope_heap
            .allocate_model_storage(vec![T_GROUP {
                num: [1, 0, 0, 0, 0],
                nGroupNumber: 1,
                nNumEndpoints: 2,
                nFirstEndpointAtNoPos: 0,
                ..T_GROUP::default()
            }])
            .unwrap();
        let taut_isotope_endpoints = taut_isotope_heap
            .allocate_model_storage(vec![0_u16, 1])
            .unwrap();
        let mut taut_isotope_info = T_GROUP_INFO {
            t_group: taut_isotope_groups,
            nEndpointAtomNumber: taut_isotope_endpoints,
            num_t_groups: 1,
            max_num_t_groups: 1,
            ..T_GROUP_INFO::default()
        };
        let mut taut_isotope_sizes = std::array::from_fn(|_| ATOM_SIZES::default());
        taut_isotope_sizes[TAUT_YES as usize] = ATOM_SIZES {
            num_isotopic_atoms: 1,
            nLenCT: 5,
            nLenCTAtOnly: 2,
            nLenLinearCTTautomer: 1,
            ..ATOM_SIZES::default()
        };
        let mut taut_isotope_bcn = BCN::default();
        let mut taut_isotope_globals = CANON_GLOBALS::default();
        assert_eq!(
            GetBaseCanonRanking(
                &mut taut_isotope_heap,
                &mut INCHI_CLOCK::default(),
                2,
                3,
                [SourceMutPointer::null(), taut_isotope_atoms],
                Some(&mut taut_isotope_info),
                &taut_isotope_sizes,
                &mut taut_isotope_bcn,
                SourceMutPointer::null(),
                &mut taut_isotope_globals,
                1,
                0,
                None,
                None,
                0,
            ),
            Ok(0)
        );
        let taut_isotope_output = &taut_isotope_bcn.ftcn[TAUT_YES as usize];
        assert_ne!(
            taut_isotope_output.nCanonFlags & CANON_FLAG_ISO_TAUT_DIFF as i32,
            0
        );
        assert_eq!(taut_isotope_output.len_iso_sort_keys, 3);
        assert_eq!(
            &taut_isotope_heap
                .slice(taut_isotope_output.iso_sort_keysOrig.as_const())
                .unwrap()[..3],
            &[0, 32768, 0]
        );
        assert_eq!(
            taut_isotope_heap
                .slice(taut_isotope_output.PartitionCtIso.Rank.as_const())
                .unwrap()
                .len(),
            3
        );
        DeAllocBCN(&mut taut_isotope_heap, Some(&mut taut_isotope_bcn)).unwrap();
        SetBitFree(&mut taut_isotope_heap, &mut taut_isotope_globals).unwrap();
        inchi_free(&mut taut_isotope_heap, taut_isotope_atoms).unwrap();
        inchi_free(&mut taut_isotope_heap, taut_isotope_groups).unwrap();
        inchi_free(&mut taut_isotope_heap, taut_isotope_endpoints).unwrap();
        assert_eq!(taut_isotope_heap.live_allocation_count(), 0);

        let mut fixed_heap = SourceHeap::default();
        let mut fixed_base_atom = atom(b"N", 0);
        fixed_base_atom.endpoint = 1;
        fixed_base_atom.num_H = 1;
        let mut fixed_other_atom = atom(b"N", 1);
        fixed_other_atom.num_H = 2;
        let fixed_base_atoms = fixed_heap
            .allocate_model_storage(vec![fixed_base_atom])
            .unwrap();
        let fixed_other_atoms = fixed_heap
            .allocate_model_storage(vec![fixed_other_atom])
            .unwrap();
        let fixed_groups = fixed_heap
            .allocate_model_storage(vec![T_GROUP {
                num: [1, 0, 0, 0, 0],
                nGroupNumber: 1,
                nNumEndpoints: 1,
                nFirstEndpointAtNoPos: 0,
                ..T_GROUP::default()
            }])
            .unwrap();
        let fixed_endpoint_numbers = fixed_heap.allocate_model_storage(vec![0_u16]).unwrap();
        let mut fixed_info = T_GROUP_INFO {
            t_group: fixed_groups,
            nEndpointAtomNumber: fixed_endpoint_numbers,
            num_t_groups: 1,
            max_num_t_groups: 1,
            ..T_GROUP_INFO::default()
        };
        let mut fixed_sizes = std::array::from_fn(|_| ATOM_SIZES::default());
        fixed_sizes[TAUT_YES as usize] = ATOM_SIZES {
            nLenCT: 3,
            nLenCTAtOnly: 1,
            nLenLinearCTTautomer: 1,
            ..ATOM_SIZES::default()
        };
        fixed_sizes[TAUT_NON as usize] = ATOM_SIZES {
            num_isotopic_atoms: 1,
            nLenCT: 1,
            nLenCTAtOnly: 1,
            ..ATOM_SIZES::default()
        };
        let mut fixed_bcn = BCN::default();
        let mut fixed_globals = CANON_GLOBALS::default();
        assert_eq!(
            GetBaseCanonRanking(
                &mut fixed_heap,
                &mut INCHI_CLOCK::default(),
                1,
                2,
                [fixed_other_atoms, fixed_base_atoms],
                Some(&mut fixed_info),
                &fixed_sizes,
                &mut fixed_bcn,
                SourceMutPointer::null(),
                &mut fixed_globals,
                1,
                0,
                None,
                None,
                0,
            ),
            Ok(0)
        );
        let fixed_taut_output = &fixed_bcn.ftcn[TAUT_YES as usize];
        assert_eq!(
            (fixed_taut_output.num_atoms, fixed_taut_output.num_at_tg),
            (1, 2)
        );
        let fixed_non_output = &fixed_bcn.ftcn[TAUT_NON as usize];
        assert_eq!(
            (fixed_non_output.num_atoms, fixed_non_output.num_at_tg),
            (1, 1)
        );
        assert_eq!(fixed_non_output.nLenNumHFixH, 1);
        assert_eq!(
            fixed_heap
                .slice(fixed_non_output.nNumHOrigFixH.as_const())
                .unwrap()[0],
            2
        );
        assert_eq!(
            fixed_heap
                .slice(fixed_non_output.nNumHFixH.as_const())
                .unwrap()[0],
            2
        );
        assert_eq!(fixed_non_output.len_iso_sort_keys, 1);
        assert_eq!(
            fixed_heap
                .slice(fixed_non_output.iso_sort_keys.as_const())
                .unwrap()[0],
            32768
        );
        assert_eq!(
            fixed_heap
                .slice(fixed_non_output.iso_sort_keysOrig.as_const())
                .unwrap()[0],
            32768
        );
        DeAllocBCN(&mut fixed_heap, Some(&mut fixed_bcn)).unwrap();
        SetBitFree(&mut fixed_heap, &mut fixed_globals).unwrap();
        inchi_free(&mut fixed_heap, fixed_base_atoms).unwrap();
        inchi_free(&mut fixed_heap, fixed_other_atoms).unwrap();
        inchi_free(&mut fixed_heap, fixed_groups).unwrap();
        inchi_free(&mut fixed_heap, fixed_endpoint_numbers).unwrap();
        assert_eq!(fixed_heap.live_allocation_count(), 0);

        let mut fixed_iso_heap = SourceHeap::default();
        let mut fixed_iso_base_first = atom(b"N", 0);
        fixed_iso_base_first.endpoint = 1;
        fixed_iso_base_first.num_H = 1;
        let fixed_iso_base_second = fixed_iso_base_first.clone();
        let mut fixed_iso_other_first = atom(b"N", 0);
        fixed_iso_other_first.num_H = 1;
        let mut fixed_iso_other_second = fixed_iso_other_first.clone();
        fixed_iso_other_second.iso_atw_diff = 1;
        let fixed_iso_base_atoms = fixed_iso_heap
            .allocate_model_storage(vec![fixed_iso_base_first, fixed_iso_base_second])
            .unwrap();
        let fixed_iso_other_atoms = fixed_iso_heap
            .allocate_model_storage(vec![fixed_iso_other_first, fixed_iso_other_second])
            .unwrap();
        let fixed_iso_groups = fixed_iso_heap
            .allocate_model_storage(vec![T_GROUP {
                num: [2, 0, 0, 0, 0],
                nGroupNumber: 1,
                nNumEndpoints: 2,
                nFirstEndpointAtNoPos: 0,
                ..T_GROUP::default()
            }])
            .unwrap();
        let fixed_iso_endpoint_numbers = fixed_iso_heap
            .allocate_model_storage(vec![0_u16, 1])
            .unwrap();
        let mut fixed_iso_info = T_GROUP_INFO {
            t_group: fixed_iso_groups,
            nEndpointAtomNumber: fixed_iso_endpoint_numbers,
            num_t_groups: 1,
            max_num_t_groups: 1,
            ..T_GROUP_INFO::default()
        };
        let mut fixed_iso_sizes = std::array::from_fn(|_| ATOM_SIZES::default());
        fixed_iso_sizes[TAUT_YES as usize] = ATOM_SIZES {
            nLenCT: 5,
            nLenCTAtOnly: 2,
            nLenLinearCTTautomer: 1,
            ..ATOM_SIZES::default()
        };
        fixed_iso_sizes[TAUT_NON as usize] = ATOM_SIZES {
            num_isotopic_atoms: 1,
            nLenCT: 2,
            nLenCTAtOnly: 2,
            ..ATOM_SIZES::default()
        };
        let mut fixed_iso_bcn = BCN::default();
        let mut fixed_iso_globals = CANON_GLOBALS::default();
        assert_eq!(
            GetBaseCanonRanking(
                &mut fixed_iso_heap,
                &mut INCHI_CLOCK::default(),
                2,
                3,
                [fixed_iso_other_atoms, fixed_iso_base_atoms],
                Some(&mut fixed_iso_info),
                &fixed_iso_sizes,
                &mut fixed_iso_bcn,
                SourceMutPointer::null(),
                &mut fixed_iso_globals,
                1,
                0,
                None,
                None,
                0,
            ),
            Ok(0)
        );
        let fixed_iso_output = &fixed_iso_bcn.ftcn[TAUT_NON as usize];
        assert_ne!(
            fixed_iso_output.nCanonFlags & CANON_FLAG_ISO_FIXED_H_DIFF as i32,
            0
        );
        assert_eq!(fixed_iso_output.len_iso_sort_keys, 2);
        assert_eq!(
            &fixed_iso_heap
                .slice(fixed_iso_output.iso_sort_keysOrig.as_const())
                .unwrap()[..2],
            &[0, 32768]
        );
        let fixed_iso_ranks = fixed_iso_heap
            .slice(fixed_iso_output.PartitionCtIso.Rank.as_const())
            .unwrap();
        assert_eq!(fixed_iso_ranks.len(), 3);
        assert_eq!(&fixed_iso_ranks[..2], &[2, 1]);
        DeAllocBCN(&mut fixed_iso_heap, Some(&mut fixed_iso_bcn)).unwrap();
        SetBitFree(&mut fixed_iso_heap, &mut fixed_iso_globals).unwrap();
        inchi_free(&mut fixed_iso_heap, fixed_iso_base_atoms).unwrap();
        inchi_free(&mut fixed_iso_heap, fixed_iso_other_atoms).unwrap();
        inchi_free(&mut fixed_iso_heap, fixed_iso_groups).unwrap();
        inchi_free(&mut fixed_iso_heap, fixed_iso_endpoint_numbers).unwrap();
        assert_eq!(fixed_iso_heap.live_allocation_count(), 0);

        for isotopic_endpoint_numbers in [vec![1_u16, 1], vec![0_u16], vec![1_u16, 0]] {
            let mut aux_heap = SourceHeap::default();
            let mut aux_endpoint = atom(b"N", 0);
            aux_endpoint.endpoint = 1;
            let mut aux_exchange = atom(b"N", 0);
            aux_exchange.cFlags = AT_FLAG_ISO_H_POINT as i8;
            let aux_atoms = aux_heap
                .allocate_model_storage(vec![aux_endpoint, aux_exchange])
                .unwrap();
            let aux_groups = aux_heap
                .allocate_model_storage(vec![T_GROUP {
                    num: [1, 0, 0, 0, 0],
                    nGroupNumber: 1,
                    nNumEndpoints: 1,
                    nFirstEndpointAtNoPos: 0,
                    ..T_GROUP::default()
                }])
                .unwrap();
            let aux_endpoints = aux_heap.allocate_model_storage(vec![0_u16]).unwrap();
            let aux_isotopic_endpoints = aux_heap
                .allocate_model_storage(isotopic_endpoint_numbers.clone())
                .unwrap();
            let mut aux_info = T_GROUP_INFO {
                t_group: aux_groups,
                nEndpointAtomNumber: aux_endpoints,
                nIsotopicEndpointAtomNumber: aux_isotopic_endpoints,
                nNumIsotopicEndpoints: isotopic_endpoint_numbers.len() as i32,
                num_t_groups: 1,
                max_num_t_groups: 1,
                bTautFlagsDone: u64::from(TG_FLAG_FOUND_ISOTOPIC_H_DONE),
                ..T_GROUP_INFO::default()
            };
            let mut aux_sizes = std::array::from_fn(|_| ATOM_SIZES::default());
            aux_sizes[TAUT_YES as usize] = ATOM_SIZES {
                nLenCT: 4,
                nLenCTAtOnly: 2,
                nLenLinearCTTautomer: 1,
                nLenIsotopicEndpoints: 2,
                ..ATOM_SIZES::default()
            };
            let mut aux_bcn = BCN::default();
            aux_bcn.ftcn[TAUT_YES as usize].nCanonFlags = 0x4000_0000;
            let mut aux_globals = CANON_GLOBALS::default();
            let result = GetBaseCanonRanking(
                &mut aux_heap,
                &mut INCHI_CLOCK::default(),
                2,
                3,
                [SourceMutPointer::null(), aux_atoms],
                Some(&mut aux_info),
                &aux_sizes,
                &mut aux_bcn,
                SourceMutPointer::null(),
                &mut aux_globals,
                1,
                0,
                None,
                None,
                0,
            );
            let valid = isotopic_endpoint_numbers == [1, 1];
            assert_eq!(result, Ok(if valid { 0 } else { CT_ISOCOUNT_ERR }));
            assert_eq!(aux_info.bIgnoreIsotopic, 0);
            if valid {
                let output = &aux_bcn.ftcn[TAUT_YES as usize];
                assert_eq!(output.len_iso_sort_keys, 0);
                assert!(output.iso_sort_keys.is_null());
                assert!(output.iso_sort_keysOrig.is_null());
                assert_eq!(
                    &aux_heap
                        .slice(output.iso_exchg_atnosOrig.as_const())
                        .unwrap()[..3],
                    &[0, 0, 0]
                );
                assert_eq!(
                    &aux_heap.slice(output.iso_exchg_atnos.as_const()).unwrap()[..3],
                    &[0, 0, 0]
                );
            } else {
                let output = &aux_bcn.ftcn[TAUT_YES as usize];
                assert_eq!(output.nCanonFlags, 0x4000_0000);
                assert!(output.LinearCt.is_null());
                assert!(output.PartitionCt.Rank.is_null());
                assert!(output.NeighList.is_null());
            }
            DeAllocBCN(&mut aux_heap, Some(&mut aux_bcn)).unwrap();
            SetBitFree(&mut aux_heap, &mut aux_globals).unwrap();
            inchi_free(&mut aux_heap, aux_atoms).unwrap();
            inchi_free(&mut aux_heap, aux_groups).unwrap();
            inchi_free(&mut aux_heap, aux_endpoints).unwrap();
            inchi_free(&mut aux_heap, aux_isotopic_endpoints).unwrap();
            assert_eq!(
                aux_heap.live_allocation_count(),
                0,
                "isotopic_endpoint_numbers={isotopic_endpoint_numbers:?}"
            );
        }

        for (name, run_case) in [
            (
                "taut_aux",
                run_taut_aux_allocation_case
                    as fn(Option<u64>) -> (Result<i32, SourceHeapError>, u64, usize),
            ),
            (
                "fixed_iso",
                run_fixed_iso_allocation_case
                    as fn(Option<u64>) -> (Result<i32, SourceHeapError>, u64, usize),
            ),
        ] {
            let (success, allocation_calls, success_live) = run_case(None);
            assert_eq!(success, Ok(0), "{name} success");
            assert_eq!(success_live, 0, "{name} success cleanup");
            for ordinal in 0..allocation_calls {
                let (result, calls, live) = run_case(Some(ordinal));
                assert!(
                    matches!(result, Ok(code) if code < 0),
                    "{name} allocation ordinal {ordinal}, calls {calls}: {result:?}"
                );
                assert_eq!(live, 0, "{name} allocation ordinal {ordinal}");
            }
        }

        let mut inactive_taut_heap = SourceHeap::default();
        let inactive_taut_atoms = inactive_taut_heap
            .allocate_model_storage(vec![atom(b"C", 0)])
            .unwrap();
        let mut inactive_taut_sizes = std::array::from_fn(|_| ATOM_SIZES::default());
        inactive_taut_sizes[TAUT_YES as usize] = ATOM_SIZES {
            nLenCT: 1,
            nLenCTAtOnly: 1,
            ..ATOM_SIZES::default()
        };
        let mut inactive_taut_bcn = BCN::default();
        let mut inactive_taut_globals = CANON_GLOBALS::default();
        assert_eq!(
            GetBaseCanonRanking(
                &mut inactive_taut_heap,
                &mut INCHI_CLOCK::default(),
                1,
                1,
                [SourceMutPointer::null(), inactive_taut_atoms],
                None,
                &inactive_taut_sizes,
                &mut inactive_taut_bcn,
                SourceMutPointer::null(),
                &mut inactive_taut_globals,
                1,
                0,
                None,
                None,
                0,
            ),
            Ok(0)
        );
        assert!(
            inactive_taut_bcn.ftcn[TAUT_NON as usize]
                .NeighList
                .is_null()
        );
        let inactive_taut_output = &inactive_taut_bcn.ftcn[TAUT_YES as usize];
        assert!(!inactive_taut_output.NeighList.is_null());
        assert_eq!(
            (
                inactive_taut_output.num_atoms,
                inactive_taut_output.num_at_tg
            ),
            (1, 1)
        );
        assert_eq!(
            inactive_taut_heap
                .slice(inactive_taut_output.LinearCt.as_const())
                .unwrap()[0],
            1
        );
        DeAllocBCN(&mut inactive_taut_heap, Some(&mut inactive_taut_bcn)).unwrap();
        SetBitFree(&mut inactive_taut_heap, &mut inactive_taut_globals).unwrap();
        inchi_free(&mut inactive_taut_heap, inactive_taut_atoms).unwrap();
        assert_eq!(inactive_taut_heap.live_allocation_count(), 0);
    }

    #[test]
    fn source_port__ichican2__setinitialranks2__line_4674() {
        fn invariant(first_value: AT_NUMB) -> ATOM_INVARIANT2 {
            let mut invariant = ATOM_INVARIANT2::default();
            invariant.val[0] = first_value;
            invariant
        }

        let mut empty_heap = SourceHeap::default();
        let mut empty_globals = CANON_GLOBALS::default();
        assert_eq!(
            SetInitialRanks2(
                &mut empty_heap,
                0,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut empty_globals,
            ),
            Ok(1)
        );
        assert!(empty_globals.m_pAtomInvariant2ForSort.is_null());

        let mut one_heap = SourceHeap::default();
        let one_invariant = one_heap.allocate_model_storage(vec![invariant(7)]).unwrap();
        let one_rank = one_heap.allocate_model_storage(vec![99_u16]).unwrap();
        let one_number = one_heap.allocate_model_storage(vec![99_u16]).unwrap();
        let mut one_globals = CANON_GLOBALS::default();
        assert_eq!(
            SetInitialRanks2(
                &mut one_heap,
                1,
                one_invariant,
                one_rank,
                one_number,
                &mut one_globals,
            ),
            Ok(1)
        );
        assert_eq!(one_heap.slice(one_number.as_const()).unwrap(), &[0]);
        assert_eq!(one_heap.slice(one_rank.as_const()).unwrap(), &[1]);
        assert_eq!(
            one_globals.m_pAtomInvariant2ForSort,
            one_invariant.as_const()
        );

        let mut heap = SourceHeap::default();
        let invariants = heap
            .allocate_model_storage(vec![invariant(5), invariant(1), invariant(1), invariant(3)])
            .unwrap();
        let ranks = heap.allocate_model_storage(vec![99_u16; 4]).unwrap();
        let atom_numbers = heap.allocate_model_storage(vec![99_u16; 4]).unwrap();
        let mut globals = CANON_GLOBALS::default();
        assert_eq!(
            SetInitialRanks2(&mut heap, 4, invariants, ranks, atom_numbers, &mut globals,),
            Ok(3)
        );
        assert_eq!(heap.slice(atom_numbers.as_const()).unwrap(), &[1, 2, 3, 0]);
        assert_eq!(heap.slice(ranks.as_const()).unwrap(), &[4, 2, 2, 3]);
        assert_eq!(globals.m_pAtomInvariant2ForSort, invariants.as_const());

        let mut second_loop_heap = SourceHeap::default();
        let mut second_loop_invariants = vec![ATOM_INVARIANT2::default(); 2];
        second_loop_invariants[0].val[7..].copy_from_slice(&[10, 20, 30]);
        second_loop_invariants[1].val[7..].copy_from_slice(&[40, 20, 60]);
        let second_loop_invariants = second_loop_heap
            .allocate_model_storage(second_loop_invariants)
            .unwrap();
        let second_loop_ranks = second_loop_heap
            .allocate_model_storage(vec![99_u16; 2])
            .unwrap();
        let second_loop_numbers = second_loop_heap
            .allocate_model_storage(vec![99_u16; 2])
            .unwrap();
        let mut second_loop_globals = CANON_GLOBALS::default();
        // Official InChI continues past unequal slot 7, then returns zero at
        // equal slot 8. SetInitialRanks2 must therefore keep one rank.
        assert_eq!(
            SetInitialRanks2(
                &mut second_loop_heap,
                2,
                second_loop_invariants,
                second_loop_ranks,
                second_loop_numbers,
                &mut second_loop_globals,
            ),
            Ok(1)
        );
        assert_eq!(
            second_loop_heap
                .slice(second_loop_numbers.as_const())
                .unwrap(),
            &[0, 1]
        );
        assert_eq!(
            second_loop_heap
                .slice(second_loop_ranks.as_const())
                .unwrap(),
            &[2, 2]
        );

        let mut error_heap = SourceHeap::default();
        let short_invariants = error_heap
            .allocate_model_storage(vec![invariant(1)])
            .unwrap();
        let untouched_ranks = error_heap.allocate_model_storage(vec![88_u16, 88]).unwrap();
        let initialized_numbers = error_heap.allocate_model_storage(vec![77_u16, 77]).unwrap();
        let mut error_globals = CANON_GLOBALS::default();
        assert_eq!(
            SetInitialRanks2(
                &mut error_heap,
                2,
                short_invariants,
                untouched_ranks,
                initialized_numbers,
                &mut error_globals,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            error_heap.slice(initialized_numbers.as_const()).unwrap(),
            &[0, 1]
        );
        assert_eq!(
            error_heap.slice(untouched_ranks.as_const()).unwrap(),
            &[88, 88]
        );
        assert_eq!(
            SetInitialRanks2(
                &mut error_heap,
                -1,
                short_invariants,
                untouched_ranks,
                initialized_numbers,
                &mut error_globals,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichican2__filloutatominvariant2__line_4719() {
        fn atom(symbol: &[u8], valence: i8, hydrogens: i8, isotope: i64) -> sp_ATOM {
            let mut atom = sp_ATOM {
                valence,
                num_H: hydrogens,
                iso_sort_key: isotope,
                ..sp_ATOM::default()
            };
            for (target, source) in atom.elname.iter_mut().zip(symbol.iter().copied()) {
                *target = source as i8;
            }
            atom
        }

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![
                atom(b"O", 2, 1, 10),
                atom(b"C", 4, 0, 20),
                atom(b"H", 0, 0, 30),
                atom(b"Br", 1, 0, 40),
                atom(b"Cl", 1, 0, 50),
                atom(b"D", 0, 0, 60),
            ])
            .unwrap();
        let invariants = heap
            .allocate_model_storage(vec![ATOM_INVARIANT2::default(); 6])
            .unwrap();
        FillOutAtomInvariant2(
            &mut heap,
            atoms.as_const(),
            6,
            6,
            invariants,
            0,
            1,
            1,
            0,
            0,
            None,
        )
        .unwrap();
        let values = heap.slice(invariants.as_const()).unwrap();
        assert_eq!(
            values.iter().map(|value| value.val).collect::<Vec<_>>(),
            vec![
                [4, 2, 1, 0, 0, 0, 0, 0, 0, 0],
                [1, 4, 0, 0, 0, 0, 0, 0, 0, 0],
                [5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [2, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                [3, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                [5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            ]
        );
        assert_eq!(
            values
                .iter()
                .map(|value| value.iso_sort_key)
                .collect::<Vec<_>>(),
            vec![10, 20, 30, 40, 50, 60]
        );

        let mut taut_heap = SourceHeap::default();
        let mut endpoint = atom(b"N", 3, 3, 11);
        endpoint.endpoint = 1;
        let taut_atoms = taut_heap
            .allocate_model_storage(vec![endpoint, atom(b"O", 2, 1, 22)])
            .unwrap();
        let groups = taut_heap
            .allocate_model_storage(vec![
                T_GROUP {
                    num: [4, 5, 6, 7, 8],
                    nNumEndpoints: 2,
                    ..T_GROUP::default()
                },
                T_GROUP {
                    num: [14, 15, 16, 17, 18],
                    nNumEndpoints: 3,
                    ..T_GROUP::default()
                },
            ])
            .unwrap();
        let info = T_GROUP_INFO {
            t_group: groups,
            num_t_groups: 2,
            ..T_GROUP_INFO::default()
        };
        let taut_invariants = taut_heap
            .allocate_model_storage(vec![ATOM_INVARIANT2::default(); 4])
            .unwrap();
        FillOutAtomInvariant2(
            &mut taut_heap,
            taut_atoms.as_const(),
            2,
            4,
            taut_invariants,
            0,
            1,
            1,
            0,
            0,
            Some(&info),
        )
        .unwrap();
        let values = taut_heap.slice(taut_invariants.as_const()).unwrap();
        assert_eq!(values[0].val, [1, 3, 0, 2, 4, 5, 3, 6, 7, 8]);
        assert_eq!(values[0].iso_sort_key, 11);
        assert_eq!(values[1].val, [2, 2, 1, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(values[1].iso_sort_key, 22);
        assert_eq!(values[2].val, [3, 21, 0, 2, 0, 0, 0, 6, 7, 8]);
        assert_eq!(values[3].val, [3, 21, 0, 3, 0, 0, 0, 16, 17, 18]);

        let digraph_invariants = taut_heap
            .allocate_model_storage(vec![ATOM_INVARIANT2::default(); 4])
            .unwrap();
        let mut isotopes_ignored = info.clone();
        isotopes_ignored.bIgnoreIsotopic = 1;
        FillOutAtomInvariant2(
            &mut taut_heap,
            taut_atoms.as_const(),
            2,
            4,
            digraph_invariants,
            1,
            1,
            1,
            1,
            0,
            Some(&isotopes_ignored),
        )
        .unwrap();
        let values = taut_heap.slice(digraph_invariants.as_const()).unwrap();
        assert_eq!(values[0].val, [1, 3, 0, 0, 0, 0, 3, 0, 0, 0]);
        assert_eq!(values[0].iso_sort_key, 0);
        assert_eq!(values[2].val, [3, 21, 0, 2, 0, 0, 0, 4, 5, 0]);
        assert_eq!(values[3].val, [3, 21, 0, 3, 0, 0, 0, 14, 15, 0]);

        let groups_only = taut_heap
            .allocate_model_storage(vec![
                ATOM_INVARIANT2 {
                    val: [99; 10],
                    iso_sort_key: 99,
                    iso_aux_key: 99,
                };
                4
            ])
            .unwrap();
        FillOutAtomInvariant2(
            &mut taut_heap,
            taut_atoms.as_const(),
            2,
            4,
            groups_only,
            0,
            0,
            0,
            0,
            1,
            Some(&info),
        )
        .unwrap();
        let values = taut_heap.slice(groups_only.as_const()).unwrap();
        assert_eq!(values[0], ATOM_INVARIANT2::default());
        assert_eq!(values[1], ATOM_INVARIANT2::default());
        assert_eq!(values[2].val, [4, 21, 0, 2, 0, 0, 0, 6, 7, 8]);
        assert_eq!(values[3].val, [4, 21, 0, 3, 0, 0, 0, 16, 17, 18]);

        let no_groups = taut_heap
            .allocate_model_storage(vec![
                ATOM_INVARIANT2 {
                    val: [77; 10],
                    ..ATOM_INVARIANT2::default()
                };
                4
            ])
            .unwrap();
        FillOutAtomInvariant2(
            &mut taut_heap,
            taut_atoms.as_const(),
            2,
            4,
            no_groups,
            0,
            0,
            0,
            0,
            0,
            None,
        )
        .unwrap();
        let values = taut_heap.slice(no_groups.as_const()).unwrap();
        assert_eq!(values[0].val, [1, 3, 0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(values[2], ATOM_INVARIANT2::default());
        assert_eq!(values[3], ATOM_INVARIANT2::default());

        assert_eq!(
            FillOutAtomInvariant2(
                &mut taut_heap,
                taut_atoms.as_const(),
                3,
                2,
                no_groups,
                0,
                0,
                0,
                0,
                0,
                None,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            FillOutAtomInvariant2(
                &mut taut_heap,
                taut_atoms.as_const(),
                -1,
                0,
                no_groups,
                0,
                0,
                0,
                0,
                0,
                None,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichican2__cleannumh__line_4881() {
        let mut heap = SourceHeap::default();
        let values = heap
            .allocate_model_storage(vec![
                EMPTY_H_NUMBER as NUM_H,
                BASE_H_NUMBER as NUM_H,
                BASE_H_NUMBER as NUM_H + 7,
                0,
                NUM_H::MIN,
                NUM_H::MAX,
                91,
            ])
            .unwrap();
        assert_eq!(CleanNumH(&mut heap, values, 6), Ok(()));
        assert_eq!(
            heap.slice(values.as_const()).unwrap(),
            &[0, 0, 7, -16383, 16385, 16384, 91]
        );

        assert_eq!(CleanNumH(&mut heap, values, 0), Ok(()));
        assert_eq!(CleanNumH(&mut heap, values, -1), Ok(()));
        assert_eq!(
            CleanNumH(&mut heap, SourceMutPointer::null(), i32::MAX),
            Ok(())
        );
        assert_eq!(
            CleanNumH(&mut heap, values, 8),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(values.as_const()).unwrap()[6], 91);
    }

    #[test]
    fn source_port__ichican2__cleanct__line_4902() {
        let mut heap = SourceHeap::default();
        let sentinel = heap
            .allocate_model_storage(vec![11_u16, 22, EMPTY_CT as AT_RANK, 44])
            .unwrap();
        assert_eq!(CleanCt(&mut heap, sentinel, 2), Ok(1));
        assert_eq!(heap.slice(sentinel.as_const()).unwrap(), &[11, 22, 0, 44]);
        assert_eq!(CleanCt(&mut heap, sentinel, 2), Ok(1));
        assert_eq!(CleanCt(&mut heap, sentinel, 1), Ok(0));
        assert_eq!(
            CleanCt(&mut heap, SourceMutPointer::null(), i32::MIN),
            Ok(0)
        );
        assert_eq!(
            CleanCt(&mut heap, sentinel, -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            CleanCt(&mut heap, sentinel, 4),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichican2__cleanisosortkeys__line_4915() {
        let mut heap = SourceHeap::default();
        let values = heap
            .allocate_model_storage(vec![AT_ISO_SORT_KEY::MAX, 0, -1, AT_ISO_SORT_KEY::MIN, 77])
            .unwrap();
        assert_eq!(CleanIsoSortKeys(&mut heap, values, 4), Ok(()));
        assert_eq!(
            heap.slice(values.as_const()).unwrap(),
            &[0, 0, -1, AT_ISO_SORT_KEY::MIN, 77]
        );
        assert_eq!(CleanIsoSortKeys(&mut heap, values, 0), Ok(()));
        assert_eq!(CleanIsoSortKeys(&mut heap, values, -1), Ok(()));
        assert_eq!(
            CleanIsoSortKeys(&mut heap, SourceMutPointer::null(), i32::MAX),
            Ok(())
        );
        assert_eq!(
            CleanIsoSortKeys(&mut heap, values, 6),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(values.as_const()).unwrap()[4], 77);
    }

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
    fn source_port__ichican2__nodesetfromvertices__line_1052() {
        let mut heap = SourceHeap::default();
        let masks = heap
            .allocate((0..16).map(|bit| 1_u16 << bit).collect())
            .unwrap();
        let globals = CANON_GLOBALS {
            m_bBit: masks,
            m_num_bit: 16,
            ..CANON_GLOBALS::default()
        };
        let row0 = heap.allocate(vec![0x1111_u16; 3]).unwrap();
        let row1 = heap.allocate(vec![0x2222_u16; 3]).unwrap();
        let rows = heap.allocate(vec![row0, row1]).unwrap();
        let set = NodeSet {
            bitword: rows,
            num_set: 2,
            len_set: 3,
        };
        let vertices = heap.allocate(vec![1_u16, 16, 17, 32, 33]).unwrap();

        assert_eq!(
            NodeSetFromVertices(&mut heap, &globals, &set, 2, vertices, 5),
            Ok(())
        );
        assert_eq!(heap.slice(row0.as_const()).unwrap(), &[0x1111; 3]);
        assert_eq!(
            heap.slice(row1.as_const()).unwrap(),
            &[0x8001, 0x8001, 0x0001]
        );

        let duplicate = heap.allocate(vec![2_u16, 2]).unwrap();
        assert_eq!(
            NodeSetFromVertices(&mut heap, &globals, &set, 1, duplicate, 2),
            Ok(())
        );
        assert_eq!(heap.slice(row0.as_const()).unwrap(), &[0x0002, 0, 0]);

        assert_eq!(
            NodeSetFromVertices(&mut heap, &globals, &set, 1, duplicate, 0),
            Ok(())
        );
        assert_eq!(heap.slice(row0.as_const()).unwrap(), &[0; 3]);
        heap.slice_mut(row0).unwrap().fill(0x3333);
        assert_eq!(
            NodeSetFromVertices(&mut heap, &globals, &set, 1, duplicate, -1),
            Ok(())
        );
        assert_eq!(heap.slice(row0.as_const()).unwrap(), &[0; 3]);

        heap.slice_mut(row0).unwrap().fill(0x4444);
        assert_eq!(
            NodeSetFromVertices(&mut heap, &globals, &set, 0, duplicate, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(row0.as_const()).unwrap(), &[0x4444; 3]);

        let partly_invalid = heap.allocate(vec![1_u16, 49]).unwrap();
        assert_eq!(
            NodeSetFromVertices(&mut heap, &globals, &set, 1, partly_invalid, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(row0.as_const()).unwrap(), &[1, 0, 0]);

        let zero_vertex = heap.allocate(vec![0_u16]).unwrap();
        assert_eq!(
            NodeSetFromVertices(&mut heap, &globals, &set, 1, zero_vertex, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(row0.as_const()).unwrap(), &[0; 3]);
    }

    #[test]
    fn source_port__ichican2__allnodesareinset__line_1071() {
        let mut heap = SourceHeap::default();
        let nodes0 = heap.allocate(vec![0x0001_u16, 0x0080, 0x4000]).unwrap();
        let nodes1 = heap.allocate(vec![0x0003_u16, 0x0100, 0x8000]).unwrap();
        let node_rows = heap.allocate(vec![nodes0, nodes1]).unwrap();
        let cur_nodes = NodeSet {
            bitword: node_rows,
            num_set: 2,
            len_set: 3,
        };
        let set0 = heap.allocate(vec![0x0001_u16, 0x0081, 0x4000]).unwrap();
        let set1 = heap.allocate(vec![0x0001_u16, 0x0100, 0x7fff]).unwrap();
        let set_rows = heap.allocate(vec![set0, set1]).unwrap();
        let set = NodeSet {
            bitword: set_rows,
            num_set: 2,
            len_set: 99,
        };

        assert_eq!(AllNodesAreInSet(&heap, &cur_nodes, 1, &set, 1), Ok(1));
        assert_eq!(AllNodesAreInSet(&heap, &cur_nodes, 2, &set, 2), Ok(0));

        heap.slice_mut(set1).unwrap()[0] = 0x0003;
        assert_eq!(AllNodesAreInSet(&heap, &cur_nodes, 2, &set, 2), Ok(0));
        heap.slice_mut(set1).unwrap()[2] = 0xffff;
        assert_eq!(AllNodesAreInSet(&heap, &cur_nodes, 2, &set, 2), Ok(1));

        let short_node = heap.allocate(vec![0x0002_u16, 0]).unwrap();
        let short_set = heap.allocate(vec![0x0001_u16]).unwrap();
        let short_node_rows = heap.allocate(vec![short_node]).unwrap();
        let short_set_rows = heap.allocate(vec![short_set]).unwrap();
        let short_cur = NodeSet {
            bitword: short_node_rows,
            num_set: 1,
            len_set: 2,
        };
        let short_container = NodeSet {
            bitword: short_set_rows,
            num_set: 1,
            len_set: 1,
        };
        assert_eq!(
            AllNodesAreInSet(&heap, &short_cur, 1, &short_container, 1),
            Ok(0)
        );
        heap.slice_mut(short_set).unwrap()[0] = 0x0002;
        assert_eq!(
            AllNodesAreInSet(&heap, &short_cur, 1, &short_container, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let negative = NodeSet {
            bitword: node_rows,
            num_set: 2,
            len_set: -1,
        };
        assert_eq!(AllNodesAreInSet(&heap, &negative, 1, &set, 1), Ok(1));
        assert_eq!(
            AllNodesAreInSet(&heap, &cur_nodes, 0, &set, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            AllNodesAreInSet(&heap, &cur_nodes, 1, &set, 3),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichican2__partitiongetmcrandfixset__line_1098() {
        let mut heap = SourceHeap::default();
        let masks = heap
            .allocate((0..16).map(|bit| 1_u16 << bit).collect())
            .unwrap();
        let globals = CANON_GLOBALS {
            m_bBit: masks,
            m_num_bit: 16,
            ..CANON_GLOBALS::default()
        };
        let mark = !rank_mask_bit();
        let ranks = heap
            .allocate(vec![4 | mark, 5, 4, 1 | mark, 7, 4 | mark, 7 | mark])
            .unwrap();
        let at_numbers = heap.allocate(vec![3_u16, 5, 2, 0, 1, 6, 4]).unwrap();
        let partition = Partition {
            Rank: ranks,
            AtNumber: at_numbers,
        };
        let mcr0 = heap.allocate(vec![0xaaaa_u16]).unwrap();
        let mcr1 = heap.allocate(vec![0xbbbb_u16]).unwrap();
        let mcr_rows = heap.allocate(vec![mcr0, mcr1]).unwrap();
        let mcr = NodeSet {
            bitword: mcr_rows,
            num_set: 2,
            len_set: 1,
        };
        let fix0 = heap.allocate(vec![0xcccc_u16]).unwrap();
        let fix1 = heap.allocate(vec![0xdddd_u16]).unwrap();
        let fix_rows = heap.allocate(vec![fix0, fix1]).unwrap();
        let fix = NodeSet {
            bitword: fix_rows,
            num_set: 2,
            len_set: 1,
        };

        assert_eq!(
            PartitionGetMcrAndFixSet(&mut heap, &globals, &partition, &mcr, &fix, 7, 2),
            Ok(())
        );
        assert_eq!(heap.slice(mcr0.as_const()).unwrap(), &[0xaaaa]);
        assert_eq!(heap.slice(fix0.as_const()).unwrap(), &[0xcccc]);
        assert_eq!(heap.slice(mcr1.as_const()).unwrap(), &[0x001b]);
        assert_eq!(heap.slice(fix1.as_const()).unwrap(), &[0x000a]);

        heap.slice_mut(mcr1).unwrap()[0] = 0xffff;
        heap.slice_mut(fix1).unwrap()[0] = 0xffff;
        assert_eq!(
            PartitionGetMcrAndFixSet(&mut heap, &globals, &partition, &mcr, &fix, 0, 2),
            Ok(())
        );
        assert_eq!(heap.slice(mcr1.as_const()).unwrap(), &[0]);
        assert_eq!(heap.slice(fix1.as_const()).unwrap(), &[0]);

        heap.slice_mut(mcr1).unwrap()[0] = 0xffff;
        heap.slice_mut(fix1).unwrap()[0] = 0xffff;
        assert_eq!(
            PartitionGetMcrAndFixSet(&mut heap, &globals, &partition, &mcr, &fix, -1, 2),
            Ok(())
        );
        assert_eq!(heap.slice(mcr1.as_const()).unwrap(), &[0]);
        assert_eq!(heap.slice(fix1.as_const()).unwrap(), &[0]);

        let short_fix = heap.allocate(Vec::<u16>::new()).unwrap();
        let short_fix_rows = heap.allocate(vec![short_fix]).unwrap();
        let short_fix_set = NodeSet {
            bitword: short_fix_rows,
            num_set: 1,
            len_set: 0,
        };
        heap.slice_mut(mcr0).unwrap()[0] = 0xffff;
        assert_eq!(
            PartitionGetMcrAndFixSet(&mut heap, &globals, &partition, &mcr, &short_fix_set, 7, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(mcr0.as_const()).unwrap(), &[0]);

        assert_eq!(
            PartitionGetMcrAndFixSet(&mut heap, &globals, &partition, &mcr, &fix, 7, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichican2__partitiongettransposition__line_1285() {
        let mut heap = SourceHeap::default();
        let from_atoms = heap.allocate(vec![2_u16, 0, 3, 1]).unwrap();
        let to_atoms = heap.allocate(vec![1_u16, 3, 0, 2]).unwrap();
        let from = Partition {
            AtNumber: from_atoms,
            ..Partition::default()
        };
        let to = Partition {
            AtNumber: to_atoms,
            ..Partition::default()
        };
        let output = heap.allocate(vec![9_u16; 4]).unwrap();
        let gamma = Transposition { nAtNumb: output };

        assert_eq!(
            PartitionGetTransposition(&mut heap, &from, &to, 4, &gamma),
            Ok(())
        );
        assert_eq!(heap.slice(output.as_const()).unwrap(), &[3, 2, 1, 0]);

        heap.slice_mut(output).unwrap().fill(8);
        assert_eq!(
            PartitionGetTransposition(&mut heap, &from, &to, 0, &gamma),
            Ok(())
        );
        assert_eq!(heap.slice(output.as_const()).unwrap(), &[8; 4]);
        assert_eq!(
            PartitionGetTransposition(&mut heap, &from, &to, -1, &gamma),
            Ok(())
        );
        assert_eq!(heap.slice(output.as_const()).unwrap(), &[8; 4]);

        let bad_from_atoms = heap.allocate(vec![1_u16, 4]).unwrap();
        let bad_to_atoms = heap.allocate(vec![7_u16, 6]).unwrap();
        let bad_from = Partition {
            AtNumber: bad_from_atoms,
            ..Partition::default()
        };
        let bad_to = Partition {
            AtNumber: bad_to_atoms,
            ..Partition::default()
        };
        let short_output = heap.allocate(vec![5_u16; 3]).unwrap();
        let short_gamma = Transposition {
            nAtNumb: short_output,
        };
        assert_eq!(
            PartitionGetTransposition(&mut heap, &bad_from, &bad_to, 2, &short_gamma),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(short_output.as_const()).unwrap(), &[5, 7, 5]);
    }

    #[test]
    fn source_port__ichican2__ngetmcr2__line_1305() {
        let mut heap = SourceHeap::default();

        let self_root = heap.allocate(vec![0_u16, 1, 2]).unwrap();
        assert_eq!(nGetMcr2(&mut heap, self_root, 2), Ok(2));
        assert_eq!(heap.slice(self_root.as_const()).unwrap(), &[0, 1, 2]);

        let chain = heap.allocate(vec![0_u16, 1, 1, 2, 3, 4]).unwrap();
        assert_eq!(nGetMcr2(&mut heap, chain, 5), Ok(1));
        assert_eq!(heap.slice(chain.as_const()).unwrap(), &[0, 1, 1, 1, 1, 1]);

        let partly_compressed = heap.allocate(vec![0_u16, 0, 1, 1, 3, 3]).unwrap();
        assert_eq!(nGetMcr2(&mut heap, partly_compressed, 5), Ok(0));
        assert_eq!(
            heap.slice(partly_compressed.as_const()).unwrap(),
            &[0, 0, 1, 0, 3, 0]
        );

        let short = heap.allocate(vec![0_u16, 1]).unwrap();
        assert_eq!(
            nGetMcr2(&mut heap, short, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(short.as_const()).unwrap(), &[0, 1]);
    }

    #[test]
    fn source_port__ichican2__njoin2mcrs2__line_1342() {
        let mut heap = SourceHeap::default();

        let lower_first = heap.allocate(vec![0_u16, 0, 2, 2, 4]).unwrap();
        assert_eq!(nJoin2Mcrs2(&mut heap, lower_first, 1, 3), Ok(1));
        assert_eq!(
            heap.slice(lower_first.as_const()).unwrap(),
            &[0, 0, 0, 2, 4]
        );

        let lower_second = heap.allocate(vec![0_u16, 1, 1, 3, 3]).unwrap();
        assert_eq!(nJoin2Mcrs2(&mut heap, lower_second, 4, 2), Ok(1));
        assert_eq!(
            heap.slice(lower_second.as_const()).unwrap(),
            &[0, 1, 1, 1, 3]
        );

        let same = heap.allocate(vec![0_u16, 0, 1, 2, 3]).unwrap();
        assert_eq!(nJoin2Mcrs2(&mut heap, same, 4, 2), Ok(0));
        assert_eq!(heap.slice(same.as_const()).unwrap(), &[0, 0, 0, 0, 0]);

        let second_error = heap.allocate(vec![0_u16, 0, 1, 2]).unwrap();
        assert_eq!(
            nJoin2Mcrs2(&mut heap, second_error, 3, 4),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(second_error.as_const()).unwrap(), &[0, 0, 0, 0]);
    }

    #[test]
    fn source_port__ichican2__getunorderedpartitionmcrnode__line_1372() {
        let mut heap = SourceHeap::default();
        let chain = heap.allocate(vec![0_u16, 0, 1, 2, 3, 4]).unwrap();
        let partition = UnorderedPartition { equ2: chain };
        assert_eq!(
            GetUnorderedPartitionMcrNode(&mut heap, &partition, 6),
            Ok(1)
        );
        assert_eq!(heap.slice(chain.as_const()).unwrap(), &[0, 0, 0, 0, 0, 0]);

        let short = heap.allocate(vec![0_u16, 1]).unwrap();
        let short_partition = UnorderedPartition { equ2: short };
        assert_eq!(
            GetUnorderedPartitionMcrNode(&mut heap, &short_partition, 3),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut boundary_values = vec![0_u16; usize::from(u16::MAX) + 1];
        boundary_values[usize::from(u16::MAX)] = u16::MAX;
        let boundary = heap.allocate(boundary_values).unwrap();
        let boundary_partition = UnorderedPartition { equ2: boundary };
        assert_eq!(
            GetUnorderedPartitionMcrNode(&mut heap, &boundary_partition, 0),
            Ok(0)
        );
        assert_eq!(
            heap.slice(boundary.as_const()).unwrap()[usize::from(u16::MAX)],
            u16::MAX
        );
    }

    #[test]
    fn source_port__ichican2__unorderedpartitionjoin__line_1385() {
        let mut heap = SourceHeap::default();
        let p1_values = heap.allocate(vec![0_u16, 0, 1, 3, 3, 4]).unwrap();
        let p2_values = heap.allocate(vec![0_u16, 1, 2, 3, 4, 5]).unwrap();
        let p1 = UnorderedPartition { equ2: p1_values };
        let p2 = UnorderedPartition { equ2: p2_values };
        assert_eq!(UnorderedPartitionJoin(&mut heap, &p1, &p2, 6), Ok(4));
        assert_eq!(
            heap.slice(p2_values.as_const()).unwrap(),
            &[0, 0, 0, 3, 3, 3]
        );

        let same_root_p1_values = heap.allocate(vec![0_u16, 0, 1]).unwrap();
        let same_root_p2_values = heap.allocate(vec![0_u16, 0, 1]).unwrap();
        let same_root_p1 = UnorderedPartition {
            equ2: same_root_p1_values,
        };
        let same_root_p2 = UnorderedPartition {
            equ2: same_root_p2_values,
        };
        assert_eq!(
            UnorderedPartitionJoin(&mut heap, &same_root_p1, &same_root_p2, 3),
            Ok(0)
        );
        assert_eq!(
            heap.slice(same_root_p2_values.as_const()).unwrap(),
            &[0, 0, 0]
        );

        let self_only = heap.allocate(vec![0_u16]).unwrap();
        let self_partition = UnorderedPartition { equ2: self_only };
        let null_partition = UnorderedPartition::default();
        assert_eq!(
            UnorderedPartitionJoin(&mut heap, &self_partition, &null_partition, 1),
            Ok(0)
        );
        assert_eq!(
            UnorderedPartitionJoin(&mut heap, &self_partition, &null_partition, 0),
            Ok(0)
        );
        assert_eq!(
            UnorderedPartitionJoin(&mut heap, &self_partition, &null_partition, -1),
            Ok(0)
        );

        let partial_p1_values = heap.allocate(vec![0_u16, 0, 5]).unwrap();
        let partial_p2_values = heap.allocate(vec![0_u16, 1, 2]).unwrap();
        let partial_p1 = UnorderedPartition {
            equ2: partial_p1_values,
        };
        let partial_p2 = UnorderedPartition {
            equ2: partial_p2_values,
        };
        assert_eq!(
            UnorderedPartitionJoin(&mut heap, &partial_p1, &partial_p2, 3),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            heap.slice(partial_p2_values.as_const()).unwrap(),
            &[0, 0, 2]
        );
    }

    #[test]
    fn source_port__ichican2__partitionsatisfieslemma_2_25__line_1407() {
        fn partition(
            heap: &mut SourceHeap,
            ranks: Vec<AT_RANK>,
            mark_indices: &[usize],
        ) -> Partition {
            let mark = !rank_mask_bit();
            let mut ranks = ranks;
            for &index in mark_indices {
                ranks[index] |= mark;
            }
            let atoms = heap
                .allocate((0..ranks.len() as AT_NUMB).collect())
                .unwrap();
            Partition {
                Rank: heap.allocate(ranks).unwrap(),
                AtNumber: atoms,
            }
        }

        let mut heap = SourceHeap::default();
        let first_condition = partition(&mut heap, vec![5, 5, 5, 5, 5], &[0, 4]);
        assert_eq!(
            PartitionSatisfiesLemma_2_25(&heap, &first_condition, 5),
            Ok(1)
        );

        let second_condition =
            partition(&mut heap, vec![2, 2, 4, 4, 6, 6, 8, 8, 10, 10], &[1, 4, 8]);
        assert_eq!(
            PartitionSatisfiesLemma_2_25(&heap, &second_condition, 10),
            Ok(1)
        );

        let third_condition = partition(
            &mut heap,
            vec![3, 3, 3, 5, 5, 7, 7, 9, 9, 11, 11],
            &[2, 6, 10],
        );
        assert_eq!(
            PartitionSatisfiesLemma_2_25(&heap, &third_condition, 11),
            Ok(1)
        );

        let false_case = partition(
            &mut heap,
            vec![4, 4, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12],
            &[3, 7, 11],
        );
        assert_eq!(PartitionSatisfiesLemma_2_25(&heap, &false_case, 12), Ok(0));

        assert_eq!(
            PartitionSatisfiesLemma_2_25(&heap, &Partition::default(), 0),
            Ok(1)
        );
        assert_eq!(
            PartitionSatisfiesLemma_2_25(&heap, &Partition::default(), -1),
            Ok(1)
        );

        let short_atoms = heap.allocate(vec![0_u16]).unwrap();
        let short_ranks = heap.allocate(vec![1_u16]).unwrap();
        let short = Partition {
            Rank: short_ranks,
            AtNumber: short_atoms,
        };
        assert_eq!(
            PartitionSatisfiesLemma_2_25(&heap, &short, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichican2__partitioncopy__line_1445() {
        let mut heap = SourceHeap::default();
        let mark = !rank_mask_bit();
        let from_atoms = heap.allocate(vec![3_u16, 1, 2, 0, 4]).unwrap();
        let from_ranks = heap
            .allocate(vec![5_u16 | mark, 4, 3 | mark, 2, 1 | mark])
            .unwrap();
        let from = Partition {
            AtNumber: from_atoms,
            Rank: from_ranks,
        };
        let to_atoms = heap.allocate(vec![90_u16, 91, 92, 93, 94, 95]).unwrap();
        let to_ranks = heap.allocate(vec![80_u16, 81, 82, 83, 84, 85]).unwrap();
        let to = Partition {
            AtNumber: to_atoms,
            Rank: to_ranks,
        };

        assert_eq!(PartitionCopy(&mut heap, &to, &from, 4), Ok(()));
        assert_eq!(
            heap.slice(to_atoms.as_const()).unwrap(),
            &[3, 1, 2, 0, 94, 95]
        );
        assert_eq!(
            heap.slice(to_ranks.as_const()).unwrap(),
            &[5, 4, 3, 2, 84, 85]
        );
        assert_eq!(heap.slice(from_atoms.as_const()).unwrap(), &[3, 1, 2, 0, 4]);
        assert_eq!(
            heap.slice(from_ranks.as_const()).unwrap(),
            &[5 | mark, 4, 3 | mark, 2, 1 | mark]
        );
        assert!(!heap.has_proven_index_bound(to_atoms.as_const(), 4, 4));

        heap.record_index_bound(from_atoms, 4, 4).unwrap();
        assert_eq!(PartitionCopy(&mut heap, &to, &from, 4), Ok(()));
        assert!(heap.has_proven_index_bound(to_atoms.as_const(), 4, 4));
        assert!(!heap.has_proven_index_bound(to_atoms.as_const(), 5, 5));

        heap.slice_mut(to_atoms).unwrap().fill(70);
        heap.slice_mut(to_ranks).unwrap().fill(60);
        assert_eq!(PartitionCopy(&mut heap, &to, &from, 0), Ok(()));
        assert_eq!(heap.slice(to_atoms.as_const()).unwrap(), &[70; 6]);
        assert_eq!(heap.slice(to_ranks.as_const()).unwrap(), &[60; 6]);
        assert_eq!(PartitionCopy(&mut heap, &to, &from, -1), Ok(()));
        assert_eq!(heap.slice(to_atoms.as_const()).unwrap(), &[70; 6]);
        assert_eq!(heap.slice(to_ranks.as_const()).unwrap(), &[60; 6]);

        let partial_from_atoms = heap.allocate(vec![7_u16, 6, 5]).unwrap();
        let partial_from_ranks = heap.allocate(vec![3_u16 | mark, 2]).unwrap();
        let partial_from = Partition {
            AtNumber: partial_from_atoms,
            Rank: partial_from_ranks,
        };
        let partial_to_atoms = heap.allocate(vec![11_u16, 12, 13]).unwrap();
        let partial_to_ranks = heap.allocate(vec![21_u16, 22, 23]).unwrap();
        let partial_to = Partition {
            AtNumber: partial_to_atoms,
            Rank: partial_to_ranks,
        };
        assert_eq!(
            PartitionCopy(&mut heap, &partial_to, &partial_from, 3),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(partial_to_atoms.as_const()).unwrap(), &[7, 6, 5]);
        assert_eq!(
            heap.slice(partial_to_ranks.as_const()).unwrap(),
            &[21, 22, 23]
        );
    }

    #[test]
    fn source_port__ichican2__partitioncolorvertex__line_1464() {
        fn fixture(
            heap: &mut SourceHeap,
            atoms: Vec<AT_NUMB>,
            ranks: Vec<AT_RANK>,
        ) -> (SourceMutPointer<NEIGH_LIST>, Vec<Partition>) {
            let atom_count = atoms.len();
            let at_number = heap.allocate(atoms).unwrap();
            let rank = heap.allocate(ranks).unwrap();
            let mut lists = Vec::with_capacity(atom_count);
            for _ in 0..atom_count {
                lists.push(heap.allocate(vec![0_u16]).unwrap());
            }
            let graph = heap.allocate(lists).unwrap();
            (
                graph,
                vec![
                    Partition {
                        Rank: rank,
                        AtNumber: at_number,
                    },
                    Partition::default(),
                    Partition::default(),
                ],
            )
        }

        fn contiguous_empty_graph(
            heap: &mut SourceHeap,
            atom_count: usize,
        ) -> SourceMutPointer<NEIGH_LIST> {
            let storage = heap
                .allocate_model_storage(vec![0_u16; atom_count])
                .unwrap();
            let rows = (0..atom_count)
                .map(|index| storage.offset(index as i64).unwrap())
                .collect();
            let graph = heap.allocate_model_storage(rows).unwrap();
            heap.record_contiguous_neighbor_layout(graph, storage, atom_count, atom_count)
                .unwrap();
            graph
        }

        for failed_allocation in 0..4 {
            let mut heap = SourceHeap::default();
            let source_atoms = heap.allocate(vec![0_u16]).unwrap();
            let source_ranks = heap.allocate(vec![1_u16]).unwrap();
            let mut partitions = vec![
                Partition {
                    Rank: source_ranks,
                    AtNumber: source_atoms,
                },
                Partition::default(),
                Partition::default(),
            ];
            heap.fail_after_allocations(failed_allocation);
            assert_eq!(
                PartitionColorVertex(
                    &mut heap,
                    &mut CANON_GLOBALS::default(),
                    SourceMutPointer::null(),
                    &mut partitions,
                    1,
                    1,
                    1,
                    1,
                    0,
                    1,
                ),
                Ok(CT_OUT_OF_RAM)
            );
            assert_eq!(
                heap.source_allocation_calls(),
                if failed_allocation < 2 { 2 } else { 4 }
            );
            assert_eq!(partitions[1].AtNumber.is_null(), failed_allocation == 0);
            assert_eq!(partitions[1].Rank.is_null(), failed_allocation == 1);
            assert_eq!(partitions[2].AtNumber.is_null(), failed_allocation <= 2);
            assert_eq!(partitions[2].Rank.is_null(), failed_allocation != 2);
        }

        let mut negative_heap = SourceHeap::default();
        let mut negative_partitions = vec![Partition::default(); 3];
        assert_eq!(
            PartitionColorVertex(
                &mut negative_heap,
                &mut CANON_GLOBALS::default(),
                SourceMutPointer::null(),
                &mut negative_partitions,
                1,
                0,
                0,
                -1,
                0,
                0,
            ),
            Ok(CT_OUT_OF_RAM)
        );
        assert!(negative_partitions[1].AtNumber.is_null());
        assert!(negative_partitions[1].Rank.is_null());

        let mut existing_heap = SourceHeap::default();
        let source_atoms = existing_heap.allocate(vec![0_u16, 1, 2, 3]).unwrap();
        let source_ranks = existing_heap.allocate(vec![4_u16; 4]).unwrap();
        let first_atoms = existing_heap.allocate(vec![90_u16; 6]).unwrap();
        let first_ranks = existing_heap.allocate(vec![80_u16; 6]).unwrap();
        let second_atoms = existing_heap.allocate(vec![70_u16; 6]).unwrap();
        let second_ranks = existing_heap.allocate(vec![60_u16; 6]).unwrap();
        let mut existing = vec![
            Partition {
                Rank: source_ranks,
                AtNumber: source_atoms,
            },
            Partition {
                Rank: first_ranks,
                AtNumber: first_atoms,
            },
            Partition {
                Rank: second_ranks,
                AtNumber: second_atoms,
            },
        ];
        existing_heap.trace_source_allocations();
        assert_eq!(
            PartitionColorVertex(
                &mut existing_heap,
                &mut CANON_GLOBALS::default(),
                SourceMutPointer::null(),
                &mut existing,
                5,
                4,
                4,
                6,
                0,
                1,
            ),
            Ok(CT_CANON_ERR)
        );
        assert_eq!(existing_heap.source_allocation_calls(), 0);
        assert_eq!(existing[1].AtNumber, first_atoms);
        assert_eq!(existing[1].Rank, first_ranks);
        assert_eq!(existing[2].AtNumber, second_atoms);
        assert_eq!(existing[2].Rank, second_ranks);
        assert_eq!(
            existing_heap.slice(first_atoms.as_const()).unwrap(),
            &[0, 1, 2, 3, 90, 90]
        );
        assert_eq!(
            existing_heap.slice(first_ranks.as_const()).unwrap(),
            &[4, 4, 4, 4, 80, 80]
        );
        assert_eq!(
            existing_heap.slice(second_atoms.as_const()).unwrap(),
            &[70; 6]
        );
        assert_eq!(
            existing_heap.slice(second_ranks.as_const()).unwrap(),
            &[60; 6]
        );

        let mut zero_vertex_heap = SourceHeap::default();
        let (zero_graph, mut zero_vertex) =
            fixture(&mut zero_vertex_heap, vec![0, 1, 2, 3], vec![4, 4, 4, 4]);
        assert_eq!(
            PartitionColorVertex(
                &mut zero_vertex_heap,
                &mut CANON_GLOBALS::default(),
                zero_graph,
                &mut zero_vertex,
                0,
                4,
                4,
                4,
                0,
                1,
            ),
            Ok(CT_CANON_ERR)
        );

        let mut absent_heap = SourceHeap::default();
        let (absent_graph, mut absent) =
            fixture(&mut absent_heap, vec![0, 2, 3, 3], vec![4, 4, 4, 4]);
        assert_eq!(
            PartitionColorVertex(
                &mut absent_heap,
                &mut CANON_GLOBALS::default(),
                absent_graph,
                &mut absent,
                2,
                4,
                4,
                4,
                0,
                1,
            ),
            Ok(CT_CANON_ERR)
        );
        assert_eq!(
            absent_heap.slice(absent[1].AtNumber.as_const()).unwrap(),
            &[0, 2, 3, 3]
        );
        assert_eq!(
            absent_heap.slice(absent[1].Rank.as_const()).unwrap(),
            &[4, 4, 4, 4]
        );

        let mut start_heap = SourceHeap::default();
        let start_atoms = start_heap.allocate(vec![0_u16, 1, 2, 3, 91, 92]).unwrap();
        let start_ranks = start_heap.allocate(vec![4_u16, 4, 4, 4, 81, 82]).unwrap();
        let start_graph = contiguous_empty_graph(&mut start_heap, 4);
        let mut start = vec![
            Partition {
                AtNumber: start_atoms,
                Rank: start_ranks,
            },
            Partition::default(),
            Partition::default(),
        ];
        start_heap
            .record_index_bound(start[0].AtNumber, 4, 4)
            .unwrap();
        let mut start_globals = CANON_GLOBALS::default();
        assert_eq!(
            PartitionColorVertex(
                &mut start_heap,
                &mut start_globals,
                start_graph,
                &mut start,
                2,
                4,
                4,
                6,
                0,
                1,
            ),
            Ok(2)
        );
        assert_eq!(
            start_heap.slice(start[1].AtNumber.as_const()).unwrap(),
            &[1, 0, 2, 3, 0, 0]
        );
        assert_eq!(
            start_heap.slice(start[1].Rank.as_const()).unwrap(),
            &[4, 1, 4, 4, 0, 0]
        );
        assert_eq!(
            start_heap.slice(start[2].Rank.as_const()).unwrap(),
            &[4, 1, 4, 4, 0, 0]
        );
        assert!(start_heap.has_proven_index_bound(start[1].AtNumber.as_const(), 4, 4));

        let mut middle_heap = SourceHeap::default();
        let (middle_graph, mut middle) =
            fixture(&mut middle_heap, vec![0, 1, 2, 3], vec![1, 4, 4, 4]);
        assert_eq!(
            PartitionColorVertex(
                &mut middle_heap,
                &mut CANON_GLOBALS::default(),
                middle_graph,
                &mut middle,
                3,
                4,
                4,
                4,
                0,
                2,
            ),
            Ok(3)
        );
        assert_eq!(
            middle_heap.slice(middle[1].AtNumber.as_const()).unwrap(),
            &[0, 2, 1, 3]
        );
        assert_eq!(
            middle_heap.slice(middle[1].Rank.as_const()).unwrap(),
            &[1, 4, 2, 4]
        );
        assert_eq!(
            middle_heap.slice(middle[2].Rank.as_const()).unwrap(),
            &[1, 4, 2, 4]
        );

        let mut digraph_heap = SourceHeap::default();
        let (digraph_graph, mut digraph) =
            fixture(&mut digraph_heap, vec![0, 1, 2, 3], vec![4, 4, 4, 4]);
        let mut digraph_globals = CANON_GLOBALS::default();
        assert_eq!(
            PartitionColorVertex(
                &mut digraph_heap,
                &mut digraph_globals,
                digraph_graph,
                &mut digraph,
                2,
                2,
                4,
                4,
                1,
                1,
            ),
            Ok(2)
        );
        assert_eq!(digraph_globals.m_nMaxAtNeighRankForSort, 2);
        assert_eq!(
            digraph_heap.slice(digraph[2].Rank.as_const()).unwrap(),
            &[4, 1, 4, 4]
        );

        let mut partial_heap = SourceHeap::default();
        let partial_atoms = partial_heap.allocate(vec![0_u16, 1, 2]).unwrap();
        let partial_ranks = partial_heap.allocate(vec![3_u16; 3]).unwrap();
        let mut partial = vec![
            Partition {
                Rank: partial_ranks,
                AtNumber: partial_atoms,
            },
            Partition::default(),
            Partition::default(),
        ];
        assert_eq!(
            PartitionColorVertex(
                &mut partial_heap,
                &mut CANON_GLOBALS::default(),
                SourceMutPointer::null(),
                &mut partial,
                1,
                3,
                3,
                3,
                0,
                1,
            ),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(
            partial_heap.slice(partial[1].AtNumber.as_const()).unwrap(),
            &[0, 1, 2]
        );
        assert_eq!(
            partial_heap.slice(partial[1].Rank.as_const()).unwrap(),
            &[1, 3, 3]
        );
    }

    #[test]
    fn source_port__ichican2__cellgetminnode__line_1586() {
        RANK_MARK_BIT.store(0x8000, Ordering::Relaxed);
        RANK_MASK_BIT.store(0x7fff, Ordering::Relaxed);
        let infinity = INCHI_CANON_INFINITY as Node;

        let empty = Cell {
            first: i32::from(infinity),
            next: 0,
            prev: 91,
        };
        assert_eq!(
            CellGetMinNode(
                &SourceHeap::default(),
                &Partition::default(),
                &empty,
                0,
                None,
            ),
            Ok(infinity)
        );

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate(vec![3_u16, 1, 2, 0]).unwrap();
        let ranks = heap
            .allocate(vec![4_u16, 4 | rank_mark_bit(), 4, 4])
            .unwrap();
        let partition = Partition {
            Rank: ranks,
            AtNumber: atoms,
        };
        let whole = Cell {
            first: 0,
            next: 4,
            prev: -7,
        };
        assert_eq!(CellGetMinNode(&heap, &partition, &whole, 0, None), Ok(1));
        assert_eq!(CellGetMinNode(&heap, &partition, &whole, 2, None), Ok(3));
        assert_eq!(CellGetMinNode(&heap, &partition, &whole, 3, None), Ok(4));
        assert_eq!(
            CellGetMinNode(&heap, &partition, &whole, 4, None),
            Ok(infinity)
        );

        let null_aux = CANON_DATA::default();
        assert_eq!(
            CellGetMinNode(&heap, &partition, &whole, 2, Some(&null_aux)),
            Ok(3)
        );
        heap.slice_mut(ranks).unwrap()[3] |= rank_mark_bit();
        assert_eq!(
            CellGetMinNode(&heap, &partition, &whole, 4, None),
            Ok(infinity)
        );

        let wrapped_first = Cell {
            first: 65_536,
            next: 1,
            prev: 0,
        };
        assert_eq!(
            CellGetMinNode(&heap, &partition, &wrapped_first, 0, None),
            Ok(infinity)
        );
        let negative_first = Cell {
            first: -1,
            next: 4,
            prev: 0,
        };
        assert_eq!(
            CellGetMinNode(&heap, &partition, &negative_first, 0, None),
            Ok(infinity)
        );

        let mut aux_heap = SourceHeap::default();
        let aux_atoms = aux_heap.allocate(vec![4_u16, 1, 3, 0, 2]).unwrap();
        let aux_ranks = aux_heap.allocate(vec![5_u16; 5]).unwrap();
        let aux_values = aux_heap.allocate(vec![5_u16, 2, 1, 2, 1]).unwrap();
        let aux_partition = Partition {
            Rank: aux_ranks,
            AtNumber: aux_atoms,
        };
        let aux_cell = Cell {
            first: 0,
            next: 5,
            prev: 44,
        };
        let aux_data = CANON_DATA {
            nAuxRank: aux_values,
            ..CANON_DATA::default()
        };
        assert_eq!(
            CellGetMinNode(&aux_heap, &aux_partition, &aux_cell, 0, Some(&aux_data)),
            Ok(3)
        );
        assert_eq!(
            CellGetMinNode(&aux_heap, &aux_partition, &aux_cell, 3, Some(&aux_data)),
            Ok(5)
        );
        assert_eq!(
            CellGetMinNode(&aux_heap, &aux_partition, &aux_cell, 5, Some(&aux_data)),
            Ok(2)
        );

        aux_heap.slice_mut(aux_ranks).unwrap()[4] |= rank_mark_bit();
        assert_eq!(
            CellGetMinNode(&aux_heap, &aux_partition, &aux_cell, 3, Some(&aux_data)),
            Ok(2)
        );
        aux_heap.slice_mut(aux_ranks).unwrap().fill(rank_mark_bit());
        assert_eq!(
            CellGetMinNode(&aux_heap, &aux_partition, &aux_cell, 0, Some(&aux_data)),
            Ok(infinity)
        );

        let short_cell = Cell {
            first: 0,
            next: 2,
            prev: 0,
        };
        assert_eq!(
            CellGetMinNode(
                &aux_heap,
                &Partition {
                    Rank: SourceMutPointer::null(),
                    AtNumber: aux_atoms,
                },
                &short_cell,
                0,
                None,
            ),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichican2__cellgetnumberofnodes__line_1748() {
        RANK_MARK_BIT.store(0x8000, Ordering::Relaxed);
        RANK_MASK_BIT.store(0x7fff, Ordering::Relaxed);

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate(vec![4_u16, 1, 3, 0, 2]).unwrap();
        let ranks = heap
            .allocate(vec![5_u16, 5 | rank_mark_bit(), 5, 5 | rank_mark_bit(), 5])
            .unwrap();
        let partition = Partition {
            Rank: ranks,
            AtNumber: atoms,
        };
        assert_eq!(
            CellGetNumberOfNodes(
                &heap,
                &partition,
                &Cell {
                    first: 0,
                    next: 5,
                    prev: 99,
                },
            ),
            Ok(3)
        );
        assert_eq!(
            CellGetNumberOfNodes(
                &heap,
                &partition,
                &Cell {
                    first: 1,
                    next: 4,
                    prev: 0,
                },
            ),
            Ok(1)
        );
        assert_eq!(
            CellGetNumberOfNodes(
                &SourceHeap::default(),
                &Partition::default(),
                &Cell {
                    first: 7,
                    next: 7,
                    prev: -1,
                },
            ),
            Ok(0)
        );
        assert_eq!(
            CellGetNumberOfNodes(
                &heap,
                &Partition {
                    Rank: SourceMutPointer::null(),
                    AtNumber: atoms,
                },
                &Cell {
                    first: 0,
                    next: 1,
                    prev: 0,
                },
            ),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(
            CellGetNumberOfNodes(
                &heap,
                &partition,
                &Cell {
                    first: -1,
                    next: 0,
                    prev: 0,
                },
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichican2__cellintersectwithset__line_1768() {
        RANK_MARK_BIT.store(0x8000, Ordering::Relaxed);
        RANK_MASK_BIT.store(0x7fff, Ordering::Relaxed);

        let mut heap = SourceHeap::default();
        let masks = heap
            .allocate(vec![
                0x0001_u16, 0x0002, 0x0004, 0x0008, 0x0010, 0x0020, 0x0040, 0x0080, 0x0100, 0x0200,
                0x0400, 0x0800, 0x1000, 0x2000, 0x4000, 0x8000,
            ])
            .unwrap();
        let globals = CANON_GLOBALS {
            m_bBit: masks,
            m_num_bit: 16,
            ..CANON_GLOBALS::default()
        };
        let set0 = heap.allocate(vec![0b0000_1010_u16]).unwrap();
        let set1 = heap.allocate(vec![0b0000_0100_u16]).unwrap();
        let rows = heap.allocate(vec![set0, set1]).unwrap();
        let mcr = NodeSet {
            bitword: rows,
            num_set: 2,
            len_set: 1,
        };
        let atoms = heap.allocate(vec![1_u16, 2, 3, 4, 2, 0]).unwrap();
        let ranks = heap
            .allocate(vec![7_u16, 7, 7 | rank_mark_bit(), 7, 7, 7])
            .unwrap();
        let partition = Partition {
            Rank: ranks,
            AtNumber: atoms,
        };

        assert_eq!(
            CellIntersectWithSet(
                &mut heap,
                &globals,
                &partition,
                &Cell {
                    first: 0,
                    next: 6,
                    prev: 99,
                },
                &mcr,
                1,
            ),
            Ok(2)
        );
        assert_eq!(
            heap.slice(ranks.as_const()).unwrap(),
            &[
                7 | rank_mark_bit(),
                7,
                7 | rank_mark_bit(),
                7,
                7 | rank_mark_bit(),
                7
            ]
        );

        heap.slice_mut(ranks).unwrap().fill(7);
        assert_eq!(
            CellIntersectWithSet(
                &mut heap,
                &globals,
                &partition,
                &Cell {
                    first: 1,
                    next: 5,
                    prev: 0,
                },
                &mcr,
                2,
            ),
            Ok(2)
        );
        assert_eq!(
            heap.slice(ranks.as_const()).unwrap(),
            &[7, 7, 7, 7 | rank_mark_bit(), 7 | rank_mark_bit(), 7]
        );

        assert_eq!(
            CellIntersectWithSet(
                &mut SourceHeap::default(),
                &CANON_GLOBALS::default(),
                &Partition::default(),
                &Cell {
                    first: 4,
                    next: 4,
                    prev: -1,
                },
                &NodeSet::default(),
                1,
            ),
            Err(SourceHeapError::NullPointer)
        );

        heap.slice_mut(ranks).unwrap().fill(7);
        let invalid_atoms = heap.allocate(vec![0_u16, 16]).unwrap();
        let invalid_partition = Partition {
            Rank: ranks,
            AtNumber: invalid_atoms,
        };
        assert_eq!(
            CellIntersectWithSet(
                &mut heap,
                &globals,
                &invalid_partition,
                &Cell {
                    first: 0,
                    next: 2,
                    prev: 0,
                },
                &mcr,
                1,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            heap.slice(ranks.as_const()).unwrap()[0],
            7 | rank_mark_bit()
        );
    }

    #[test]
    fn source_port__ichican2__ctpartclear__line_1798() {
        let mut heap = SourceHeap::default();
        let table = heap.allocate(vec![10_u16, 11, 12, 13, 14, 15, 16]).unwrap();
        let positions = heap.allocate(vec![2_u16, 4, 6]).unwrap();
        let mut ct = ConTable {
            Ctbl: table,
            nextCtblPos: positions,
            lenCt: 7,
            lenPos: 91,
            lenNumH: 92,
            ..ConTable::default()
        };

        assert_eq!(CtPartClear(&mut heap, &mut ct, 2), Ok(()));
        assert_eq!(
            heap.slice(table.as_const()).unwrap(),
            &[10, 11, 12, 13, 0, 0, 0]
        );
        assert_eq!(ct.lenCt, 4);
        assert_eq!(ct.lenPos, 2);
        assert_eq!(ct.lenNumH, 92);

        heap.slice_mut(table)
            .unwrap()
            .copy_from_slice(&[20, 21, 22, 23, 24, 25, 26]);
        ct.lenCt = 7;
        assert_eq!(CtPartClear(&mut heap, &mut ct, 1), Ok(()));
        assert_eq!(heap.slice(table.as_const()).unwrap(), &[0; 7]);
        assert_eq!(ct.lenCt, 0);
        assert_eq!(ct.lenPos, 1);

        heap.slice_mut(table)
            .unwrap()
            .copy_from_slice(&[30, 31, 32, 33, 34, 35, 36]);
        ct.lenCt = 4;
        assert_eq!(CtPartClear(&mut heap, &mut ct, 3), Ok(()));
        assert_eq!(
            heap.slice(table.as_const()).unwrap(),
            &[30, 31, 32, 33, 34, 35, 36]
        );
        assert_eq!(ct.lenCt, 6);
        assert_eq!(ct.lenPos, 3);

        ct.Ctbl = SourceMutPointer::null();
        ct.lenCt = 3;
        assert_eq!(CtPartClear(&mut heap, &mut ct, 3), Ok(()));
        assert_eq!(ct.lenCt, 6);
        assert_eq!(ct.lenPos, 3);

        ct.nextCtblPos = SourceMutPointer::null();
        ct.lenCt = 5;
        ct.lenPos = 77;
        assert_eq!(
            CtPartClear(&mut heap, &mut ct, 2),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(ct.lenCt, 5);
        assert_eq!(ct.lenPos, 77);
    }

    #[test]
    fn source_port__ichican2__insertions_sort_neighlist_at_numbers2__line_1819() {
        RANK_MARK_BIT.store(0x8000, Ordering::Relaxed);
        RANK_MASK_BIT.store(0x7fff, Ordering::Relaxed);

        let mut heap = SourceHeap::default();
        let ranks = heap
            .allocate(vec![30_u16, 10 | rank_mark_bit(), 20, 10, 40, 5])
            .unwrap();
        let list = heap.allocate(vec![6_u16, 0, 1, 2, 3, 4, 5]).unwrap();
        assert_eq!(
            insertions_sort_NeighList_AT_NUMBERS2(&mut heap, list, ranks, u16::MAX),
            Ok(())
        );
        assert_eq!(heap.slice(list.as_const()).unwrap(), &[6, 5, 1, 3, 2, 0, 4]);

        let gated = heap.allocate(vec![5_u16, 0, 4, 2, 5, 1]).unwrap();
        assert_eq!(
            insertions_sort_NeighList_AT_NUMBERS2(&mut heap, gated, ranks, 25),
            Ok(())
        );
        assert_eq!(heap.slice(gated.as_const()).unwrap(), &[5, 5, 1, 2, 0, 4]);

        let empty = heap.allocate(vec![0_u16]).unwrap();
        let singleton = heap.allocate(vec![1_u16, 4]).unwrap();
        assert_eq!(
            insertions_sort_NeighList_AT_NUMBERS2(&mut heap, empty, SourceMutPointer::null(), 1,),
            Ok(())
        );
        assert_eq!(
            insertions_sort_NeighList_AT_NUMBERS2(
                &mut heap,
                singleton,
                SourceMutPointer::null(),
                1,
            ),
            Ok(())
        );
        assert_eq!(heap.slice(singleton.as_const()).unwrap(), &[1, 4]);

        let short = heap.allocate(vec![3_u16, 4, 2]).unwrap();
        assert_eq!(
            insertions_sort_NeighList_AT_NUMBERS2(&mut heap, short, ranks, u16::MAX),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(short.as_const()).unwrap(), &[3, 2, 4]);
    }

    #[test]
    fn source_port__ichican2__ctpartfill__line_1849() {
        RANK_MARK_BIT.store(0x8000, Ordering::Relaxed);
        RANK_MASK_BIT.store(0x7fff, Ordering::Relaxed);

        let mut heap = SourceHeap::default();
        let ranks = heap.allocate(vec![2_u16, 3, 1, 4, 6]).unwrap();
        let order = heap.allocate(vec![2_u16, 0, 1, 3, 4]).unwrap();
        let partition = Partition {
            Rank: ranks,
            AtNumber: order,
        };
        let graph0 = heap.allocate(vec![1_u16, 2]).unwrap();
        let graph1 = heap.allocate(vec![3_u16, 0, 4, 2]).unwrap();
        let graph2 = heap.allocate(vec![0_u16]).unwrap();
        let graph3 = heap.allocate(vec![2_u16, 1, 0]).unwrap();
        let graph4 = heap.allocate(vec![0_u16]).unwrap();
        let graph = heap
            .allocate(vec![graph0, graph1, graph2, graph3, graph4])
            .unwrap();

        let source_h = heap.allocate(vec![100_i16, 101, 102, 130, 131]).unwrap();
        let source_fixed_h = heap.allocate(vec![200_i16, 201, 202, 230, 231]).unwrap();
        let source_iso = heap.allocate(vec![300_i64, 301, 302, 303, 304]).unwrap();
        let source_exchange = heap.allocate(vec![10_i8, 11, 12, 13, 14]).unwrap();
        let canon = CANON_DATA {
            NumH: source_h,
            NumHfixed: source_fixed_h,
            iso_sort_key: source_iso,
            iso_exchg_atnos: source_exchange,
            ..CANON_DATA::default()
        };

        let table = heap.allocate(vec![999_u16; 12]).unwrap();
        let positions = heap.allocate(vec![88_u16; 2]).unwrap();
        let next_ranks = heap.allocate(vec![77_u16; 2]).unwrap();
        let target_h = heap.allocate(vec![-1_i16; 7]).unwrap();
        let target_fixed_h = heap.allocate(vec![-2_i16; 7]).unwrap();
        let target_iso = heap.allocate(vec![-3_i64; 7]).unwrap();
        let target_exchange = heap.allocate(vec![-4_i8; 7]).unwrap();
        let mut ct = ConTable {
            Ctbl: table,
            nextCtblPos: positions,
            nextAtRank: next_ranks,
            NumH: target_h,
            lenNumH: 90,
            NumHfixed: target_fixed_h,
            iso_sort_key: target_iso,
            len_iso_sort_key: 91,
            iso_exchg_atnos: target_exchange,
            len_iso_exchg_atnos: 92,
            lenCt: 93,
            lenPos: 94,
            ..ConTable::default()
        };

        assert_eq!(
            CtPartFill(
                &mut heap,
                graph,
                &canon,
                &partition,
                Some(&mut ct),
                1,
                3,
                5,
                5,
            ),
            Ok(())
        );
        assert_eq!(
            heap.slice(table.as_const()).unwrap(),
            &[1, 2, 1, 3, 1, 2, 4, 2, 3, 999, 999, 999]
        );
        assert_eq!(heap.slice(graph1.as_const()).unwrap(), &[3, 2, 0, 4]);
        assert_eq!(heap.slice(graph3.as_const()).unwrap(), &[2, 0, 1]);
        assert_eq!(
            heap.slice(target_h.as_const()).unwrap(),
            &[102, 100, 101, 130, 131, -1, -1]
        );
        assert_eq!(
            heap.slice(target_fixed_h.as_const()).unwrap(),
            &[202, 200, 201, -2, -2, -2, -2]
        );
        assert_eq!(
            heap.slice(target_iso.as_const()).unwrap(),
            &[302, 300, 301, 303, -3, -3, -3]
        );
        assert_eq!(
            heap.slice(target_exchange.as_const()).unwrap(),
            &[12, 10, 11, 13, -4, -4, -4]
        );
        assert_eq!(ct.lenCt, 9);
        assert_eq!(ct.lenNumH, 5);
        assert_eq!(ct.len_iso_sort_key, 4);
        assert_eq!(ct.len_iso_exchg_atnos, 4);
        assert_eq!(ct.lenPos, 1);
        assert_eq!(heap.slice(positions.as_const()).unwrap(), &[9, 88]);
        assert_eq!(heap.slice(next_ranks.as_const()).unwrap(), &[5, 77]);

        let continuation_positions = heap.allocate(vec![9_u16, 66]).unwrap();
        let continuation_ranks = heap.allocate(vec![5_u16, 55]).unwrap();
        let mut continuation = ConTable {
            Ctbl: table,
            nextCtblPos: continuation_positions,
            nextAtRank: continuation_ranks,
            lenNumH: 7,
            len_iso_sort_key: 8,
            len_iso_exchg_atnos: 9,
            ..ConTable::default()
        };
        assert_eq!(
            CtPartFill(
                &mut heap,
                graph,
                &CANON_DATA::default(),
                &partition,
                Some(&mut continuation),
                2,
                3,
                4,
                5,
            ),
            Ok(())
        );
        assert_eq!(continuation.lenCt, 9);
        assert_eq!(continuation.lenNumH, 0);
        assert_eq!(continuation.len_iso_sort_key, 0);
        assert_eq!(continuation.len_iso_exchg_atnos, 0);
        assert_eq!(continuation.lenPos, 2);
        assert_eq!(
            heap.slice(continuation_positions.as_const()).unwrap(),
            &[9, 9]
        );
        assert_eq!(heap.slice(continuation_ranks.as_const()).unwrap(), &[5, 6]);

        assert_eq!(
            CtPartFill(
                &mut SourceHeap::default(),
                SourceMutPointer::null(),
                &CANON_DATA::default(),
                &Partition::default(),
                None,
                1,
                0,
                0,
                0,
            ),
            Ok(())
        );

        let negative_positions = heap.allocate(vec![0_u16]).unwrap();
        let negative_ranks = heap.allocate(vec![0_u16]).unwrap();
        let mut negative = ConTable {
            nextCtblPos: negative_positions,
            nextAtRank: negative_ranks,
            lenCt: 41,
            lenPos: 42,
            ..ConTable::default()
        };
        assert_eq!(
            CtPartFill(
                &mut heap,
                SourceMutPointer::null(),
                &CANON_DATA::default(),
                &Partition::default(),
                Some(&mut negative),
                2,
                0,
                0,
                0,
            ),
            Ok(())
        );
        assert_eq!(negative.lenCt, 41);
        assert_eq!(negative.lenPos, 42);

        let mut limited = ConTable {
            lenCt: 51,
            lenPos: 52,
            ..ConTable::default()
        };
        assert_eq!(
            CtPartFill(
                &mut heap,
                SourceMutPointer::null(),
                &CANON_DATA::default(),
                &Partition::default(),
                Some(&mut limited),
                1,
                0,
                0,
                0,
            ),
            Ok(())
        );
        assert_eq!(limited.lenCt, 51);
        assert_eq!(limited.lenPos, 52);

        let short_table = heap.allocate(vec![700_u16]).unwrap();
        let short_positions = heap.allocate(vec![0_u16]).unwrap();
        let short_ranks = heap.allocate(vec![0_u16]).unwrap();
        let mut partial = ConTable {
            Ctbl: short_table,
            nextCtblPos: short_positions,
            nextAtRank: short_ranks,
            ..ConTable::default()
        };
        assert_eq!(
            CtPartFill(
                &mut heap,
                graph,
                &CANON_DATA::default(),
                &partition,
                Some(&mut partial),
                1,
                3,
                5,
                5,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(short_table.as_const()).unwrap(), &[1]);
        assert_eq!(partial.lenCt, 0);
        assert_eq!(partial.lenPos, 0);
    }

    #[test]
    fn source_port__ichican2__ctpartinchi_canon_infinity__line_2065() {
        let mut heap = SourceHeap::default();
        let table = heap.allocate(vec![11_u16, 12, 13, 14, 15, 16]).unwrap();
        let positions = heap.allocate(vec![2_u16, 4, 5]).unwrap();
        let cmp = heap.allocate(vec![9_i8, 8, 7, 6]).unwrap();
        let ct = ConTable {
            Ctbl: table,
            nextCtblPos: positions,
            lenCt: 6,
            lenPos: 3,
            ..ConTable::default()
        };

        assert_eq!(
            CtPartINCHI_CANON_INFINITY(&mut heap, &ct, SourceMutPointer::null(), 1),
            Ok(())
        );
        assert_eq!(
            heap.slice(table.as_const()).unwrap(),
            &[0, 12, 13, 14, 15, 16]
        );
        assert_eq!(ct.lenCt, 6);
        assert_eq!(ct.lenPos, 3);

        assert_eq!(CtPartINCHI_CANON_INFINITY(&mut heap, &ct, cmp, 3), Ok(()));
        assert_eq!(heap.slice(cmp.as_const()).unwrap(), &[0, 0, 7, 6]);
        assert_eq!(
            heap.slice(table.as_const()).unwrap(),
            &[0, 12, 13, 14, 0, 16]
        );

        heap.slice_mut(table).unwrap()[3] = 0;
        heap.slice_mut(table).unwrap()[4] = 55;
        heap.slice_mut(cmp).unwrap().fill(5);
        assert_eq!(CtPartINCHI_CANON_INFINITY(&mut heap, &ct, cmp, 3), Ok(()));
        assert_eq!(heap.slice(cmp.as_const()).unwrap(), &[0, 0, 5, 5]);
        assert_eq!(heap.slice(table.as_const()).unwrap()[4], 55);

        let short_cmp = heap.allocate(vec![4_i8]).unwrap();
        heap.slice_mut(table).unwrap()[4] = 66;
        assert_eq!(
            CtPartINCHI_CANON_INFINITY(&mut heap, &ct, short_cmp, 3),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(short_cmp.as_const()).unwrap(), &[4]);
        assert_eq!(heap.slice(table.as_const()).unwrap()[4], 66);

        let bad_positions = ConTable {
            Ctbl: table,
            nextCtblPos: SourceMutPointer::null(),
            ..ConTable::default()
        };
        assert_eq!(
            CtPartINCHI_CANON_INFINITY(&mut heap, &bad_positions, cmp, 2),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichican2__ctfullcompare__line_2623() {
        fn make_ct(
            heap: &mut SourceHeap,
            ctbl: &[u16],
            h: &[i16],
            fixed_h: &[i16],
            iso: &[i64],
            atoms_only: i32,
        ) -> ConTable {
            let len_ct = ctbl.len() as i32;
            let atom_len = h.len() as i32;
            ConTable {
                Ctbl: heap.allocate(ctbl.to_vec()).unwrap(),
                lenCt: len_ct,
                nLenCTAtOnly: atoms_only,
                maxlenCt: len_ct,
                maxVert: atom_len,
                lenPos: 1,
                nextCtblPos: heap.allocate(vec![len_ct as u16]).unwrap(),
                nextAtRank: heap.allocate(vec![(atom_len + 1) as u16]).unwrap(),
                NumH: heap.allocate(h.to_vec()).unwrap(),
                lenNumH: atom_len,
                NumHfixed: heap.allocate(fixed_h.to_vec()).unwrap(),
                iso_sort_key: heap.allocate(iso.to_vec()).unwrap(),
                len_iso_sort_key: atom_len,
                ..ConTable::default()
            }
        }

        let mut heap = SourceHeap::default();
        let base = make_ct(&mut heap, &[1, 2, 3, 4], &[10, 11], &[20, 21], &[30, 31], 2);
        let equal = make_ct(&mut heap, &[1, 2, 3, 4], &[10, 11], &[20, 21], &[30, 31], 2);
        assert_eq!(CtFullCompare(&mut heap, &base, &equal, 0, 1), Ok(0));

        let shorter = make_ct(&mut heap, &[1, 2, 3], &[10, 11], &[20, 21], &[30, 31], 2);
        assert_eq!(CtFullCompare(&mut heap, &base, &shorter, 0, 1), Ok(1));

        let layer0 = make_ct(&mut heap, &[2, 2, 3, 4], &[10, 11], &[20, 21], &[30, 31], 2);
        assert_eq!(CtFullCompare(&mut heap, &base, &layer0, 0, 1), Ok(1));

        let layer1 = make_ct(&mut heap, &[1, 2, 3, 4], &[11, 11], &[20, 21], &[30, 31], 2);
        assert_eq!(CtFullCompare(&mut heap, &base, &layer1, 0, 1), Ok(2));

        let layer2 = make_ct(&mut heap, &[1, 2, 4, 4], &[10, 11], &[20, 21], &[30, 31], 2);
        assert_eq!(CtFullCompare(&mut heap, &base, &layer2, 0, 1), Ok(3));
        assert_eq!(CtFullCompare(&mut heap, &base, &layer2, 0, 0), Ok(1));

        let layer3 = make_ct(&mut heap, &[1, 2, 3, 4], &[10, 11], &[21, 21], &[30, 31], 2);
        assert_eq!(CtFullCompare(&mut heap, &base, &layer3, 0, 1), Ok(4));

        let layer4 = make_ct(&mut heap, &[1, 2, 3, 4], &[10, 11], &[20, 21], &[31, 31], 2);
        assert_eq!(CtFullCompare(&mut heap, &base, &layer4, 0, 1), Ok(5));

        let common_left = ConTable {
            Ctbl: heap.allocate(vec![1_u16, 2, 3, 9, 0]).unwrap(),
            lenCt: 4,
            nLenCTAtOnly: 2,
            maxVert: 2,
            lenPos: 1,
            nextCtblPos: heap.allocate(vec![4_u16]).unwrap(),
            nextAtRank: heap.allocate(vec![3_u16]).unwrap(),
            ..ConTable::default()
        };
        let common_right = ConTable {
            Ctbl: heap.allocate(vec![1_u16, 2, 3, 8, 0]).unwrap(),
            lenCt: 4,
            nLenCTAtOnly: 2,
            maxVert: 2,
            lenPos: 1,
            nextCtblPos: heap.allocate(vec![4_u16]).unwrap(),
            nextAtRank: heap.allocate(vec![3_u16]).unwrap(),
            ..ConTable::default()
        };
        assert_eq!(
            CtFullCompare(&mut heap, &common_left, &common_right, 1, 1),
            Ok(0)
        );

        let invalid = ConTable {
            lenPos: 1,
            nextCtblPos: SourceMutPointer::null(),
            ..ConTable::default()
        };
        assert_eq!(
            CtFullCompare(&mut heap, &invalid, &equal, 0, 1),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichican2__ctfullcomparelayers__line_2863() {
        let mut layers: KLeastLayers = std::array::from_fn(|_| kLeast::default());
        assert_eq!(CtFullCompareLayers(&layers), Ok(0));

        layers[4] = kLeast { k: 7, i: 91 };
        assert_eq!(CtFullCompareLayers(&layers), Ok(5));

        layers[4].k = -7;
        assert_eq!(CtFullCompareLayers(&layers), Ok(-5));

        layers[1] = kLeast { k: 1, i: -1 };
        assert_eq!(CtFullCompareLayers(&layers), Ok(2));
    }

    #[test]
    fn source_port__ichican2__ctcomparelayersgetfirstdiff__line_2880() {
        let mut layers: KLeastLayers = std::array::from_fn(|_| kLeast::default());
        let (mut layer, mut item, mut k) = (91, 92, 93);

        assert_eq!(
            CtCompareLayersGetFirstDiff(None, 7, &mut layer, &mut item, &mut k,),
            Ok(-1)
        );
        assert_eq!((layer, item, k), (91, 92, 93));

        assert_eq!(
            CtCompareLayersGetFirstDiff(Some(&layers), 7, &mut layer, &mut item, &mut k,),
            Ok(0)
        );
        assert_eq!((layer, item, k), (7, -1, 0));

        assert_eq!(
            CtCompareLayersGetFirstDiff(Some(&layers), 0, &mut layer, &mut item, &mut k,),
            Ok(0)
        );
        assert_eq!((layer, item, k), (INCHI_CANON_INFINITY as i32, -1, 0));

        layers[8] = kLeast { k: -11, i: 17 };
        assert_eq!(
            CtCompareLayersGetFirstDiff(Some(&layers), 4, &mut layer, &mut item, &mut k,),
            Ok(1)
        );
        assert_eq!((layer, item, k), (8, 17, -11));

        layers[2] = kLeast { k: 5, i: 6 };
        assert_eq!(
            CtCompareLayersGetFirstDiff(Some(&layers), 4, &mut layer, &mut item, &mut k,),
            Ok(1)
        );
        assert_eq!((layer, item, k), (2, 6, 5));
    }

    #[test]
    fn source_port__ichican2__ctpartcomparelayers__line_2929() {
        let mut layers: KLeastLayers = std::array::from_fn(|_| kLeast::default());
        assert_eq!(CtPartCompareLayers(Some(&layers), 99, 0), Ok(0));
        assert_eq!(CtPartCompareLayers(Some(&layers), 99, 7), Ok(0));
        assert_eq!(CtPartCompareLayers(None, 99, 0), Ok(0));

        layers[4] = kLeast { k: 8, i: 3 };
        assert_eq!(CtPartCompareLayers(Some(&layers), 4, 0), Ok(5));
        assert_eq!(CtPartCompareLayers(Some(&layers), 3, 0), Ok(0));

        layers[4].k = -8;
        assert_eq!(CtPartCompareLayers(Some(&layers), 4, 0), Ok(-5));
    }

    #[test]
    fn source_port__ichican2__updatecomparelayers__line_2947() {
        let mut layers: KLeastLayers = std::array::from_fn(|_| kLeast::default());
        layers[0] = kLeast { k: 2, i: 20 };
        layers[1] = kLeast { k: 3, i: 30 };
        layers[2] = kLeast { k: -4, i: 40 };

        assert_eq!(UpdateCompareLayers(Some(&mut layers), 3), Ok(()));
        assert_eq!(layers[0], kLeast { k: 2, i: 20 });
        assert_eq!(layers[1], kLeast::default());
        assert_eq!(layers[2], kLeast::default());

        layers[3] = kLeast { k: 1, i: 9 };
        assert_eq!(UpdateCompareLayers(Some(&mut layers), 0), Ok(()));
        assert!(layers.iter().all(|value| *value == kLeast::default()));

        assert_eq!(UpdateCompareLayers(None, 1), Ok(()));
    }

    #[test]
    fn source_port__ichican2__transpositiongetmcrandfixsetandunorderedpartition__line_3125() {
        let mut heap = SourceHeap::default();
        let mut globals = CANON_GLOBALS::default();
        assert_eq!(SetBitCreate(&mut heap, &mut globals), Ok(1));

        let mcr_row = heap.allocate(vec![u16::MAX]).unwrap();
        let fix_row = heap.allocate(vec![u16::MAX]).unwrap();
        let mcr_rows = heap.allocate(vec![mcr_row]).unwrap();
        let fix_rows = heap.allocate(vec![fix_row]).unwrap();
        let mcr_set = NodeSet {
            bitword: mcr_rows,
            num_set: 1,
            len_set: 1,
        };
        let fix_set = NodeSet {
            bitword: fix_rows,
            num_set: 1,
            len_set: 1,
        };
        let permutation = heap.allocate(vec![0_u16, 2, 1, 4, 5, 3]).unwrap();
        let mut gamma = Transposition {
            nAtNumb: permutation,
        };
        let equ2 = heap.allocate(vec![99_u16; 6]).unwrap();
        let mut partition = UnorderedPartition { equ2 };

        assert_eq!(
            TranspositionGetMcrAndFixSetAndUnorderedPartition(
                &mut heap,
                &globals,
                &mut gamma,
                &mcr_set,
                &fix_set,
                6,
                1,
                &mut partition,
            ),
            Ok(())
        );
        assert_eq!(heap.slice(mcr_row.as_const()).unwrap(), &[0b001011]);
        assert_eq!(heap.slice(fix_row.as_const()).unwrap(), &[0b000001]);
        assert_eq!(heap.slice(equ2.as_const()).unwrap(), &[0, 1, 1, 3, 3, 3]);
        assert_eq!(
            heap.slice(permutation.as_const()).unwrap(),
            &[0, 2, 1, 4, 5, 3]
        );

        heap.slice_mut(permutation)
            .unwrap()
            .copy_from_slice(&[0, 1, 2, 3, 4, 5]);
        assert_eq!(
            TranspositionGetMcrAndFixSetAndUnorderedPartition(
                &mut heap,
                &globals,
                &mut gamma,
                &mcr_set,
                &fix_set,
                6,
                1,
                &mut partition,
            ),
            Ok(())
        );
        assert_eq!(heap.slice(mcr_row.as_const()).unwrap(), &[0b111111]);
        assert_eq!(heap.slice(fix_row.as_const()).unwrap(), &[0b111111]);
        assert_eq!(heap.slice(equ2.as_const()).unwrap(), &[0, 1, 2, 3, 4, 5]);

        heap.slice_mut(mcr_row).unwrap()[0] = u16::MAX;
        heap.slice_mut(fix_row).unwrap()[0] = u16::MAX;
        assert_eq!(
            TranspositionGetMcrAndFixSetAndUnorderedPartition(
                &mut heap,
                &CANON_GLOBALS::default(),
                &mut Transposition::default(),
                &mcr_set,
                &fix_set,
                0,
                1,
                &mut UnorderedPartition::default(),
            ),
            Ok(())
        );
        assert_eq!(heap.slice(mcr_row.as_const()).unwrap(), &[0]);
        assert_eq!(heap.slice(fix_row.as_const()).unwrap(), &[0]);

        let bad_set = NodeSet {
            bitword: SourceMutPointer::null(),
            num_set: 1,
            len_set: 1,
        };
        assert_eq!(
            TranspositionGetMcrAndFixSetAndUnorderedPartition(
                &mut heap,
                &globals,
                &mut gamma,
                &bad_set,
                &fix_set,
                6,
                1,
                &mut partition,
            ),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichican2__getoneadditionallayer__line_3566() {
        let mut heap = SourceHeap::default();
        let num_h = heap.allocate(vec![1_i16]).unwrap();
        let fixed_h = heap.allocate(vec![2_i16]).unwrap();
        let iso = heap.allocate(vec![3_i64]).unwrap();

        assert_eq!(GetOneAdditionalLayer(None, None), 0);
        assert_eq!(GetOneAdditionalLayer(Some(&CANON_DATA::default()), None), 0);
        assert_eq!(GetOneAdditionalLayer(None, Some(&ConTable::default())), 0);

        let data_l1 = CANON_DATA {
            NumH: num_h,
            ..CANON_DATA::default()
        };
        assert_eq!(
            GetOneAdditionalLayer(Some(&data_l1), Some(&ConTable::default())),
            1
        );

        let data_l2 = CANON_DATA {
            nLenCTAtOnly: 1,
            nLenLinearCT: 2,
            ..CANON_DATA::default()
        };
        let fixed_l2 = ConTable {
            nLenCTAtOnly: 3,
            lenCt: 3,
            ..ConTable::default()
        };
        assert_eq!(GetOneAdditionalLayer(Some(&data_l2), Some(&fixed_l2)), 2);

        let data_l3 = CANON_DATA {
            NumHfixed: fixed_h,
            ..CANON_DATA::default()
        };
        assert_eq!(
            GetOneAdditionalLayer(Some(&data_l3), Some(&ConTable::default())),
            3
        );

        let data_l4 = CANON_DATA {
            iso_sort_key: iso,
            ..CANON_DATA::default()
        };
        assert_eq!(
            GetOneAdditionalLayer(Some(&data_l4), Some(&ConTable::default())),
            4
        );

        let multiple = CANON_DATA {
            NumH: num_h,
            NumHfixed: fixed_h,
            ..CANON_DATA::default()
        };
        assert_eq!(
            GetOneAdditionalLayer(Some(&multiple), Some(&ConTable::default())),
            0
        );

        let no_missing = ConTable {
            NumH: num_h,
            NumHfixed: fixed_h,
            iso_sort_key: iso,
            ..ConTable::default()
        };
        let all_data = CANON_DATA {
            NumH: num_h,
            NumHfixed: fixed_h,
            iso_sort_key: iso,
            ..CANON_DATA::default()
        };
        assert_eq!(GetOneAdditionalLayer(Some(&all_data), Some(&no_missing)), 0);
    }

    #[test]
    fn source_port__ichican2__transpositioncreate__line_694() {
        let mut heap = SourceHeap::default();
        let old = heap.allocate(vec![7_u16, 8]).unwrap();
        let mut transposition = Transposition { nAtNumb: old };

        assert_eq!(TranspositionCreate(&mut heap, &mut transposition, 4), Ok(1));
        assert_ne!(transposition.nAtNumb, old);
        assert_eq!(
            heap.slice(transposition.nAtNumb.as_const()).unwrap(),
            &[0, 0, 0, 0]
        );
        assert_eq!(heap.slice(old.as_const()).unwrap(), &[7, 8]);

        let successful = transposition.nAtNumb;
        heap.fail_after_allocations(0);
        assert_eq!(TranspositionCreate(&mut heap, &mut transposition, 3), Ok(0));
        assert!(transposition.nAtNumb.is_null());
        assert_eq!(heap.slice(successful.as_const()).unwrap(), &[0, 0, 0, 0]);

        let mut negative_heap = SourceHeap::default();
        let mut negative = Transposition::default();
        assert_eq!(
            TranspositionCreate(&mut negative_heap, &mut negative, -1),
            Ok(0)
        );
        assert!(negative.nAtNumb.is_null());
    }

    #[test]
    fn source_port__ichican2__transpositionfree__line_707() {
        let mut heap = SourceHeap::default();
        let allocation = heap.allocate(vec![1_u16, 2, 3]).unwrap();
        let mut transposition = Transposition {
            nAtNumb: allocation,
        };

        assert_eq!(
            TranspositionFree(&mut heap, Some(&mut transposition)),
            Ok(())
        );
        assert!(transposition.nAtNumb.is_null());
        assert!(heap.slice(allocation.as_const()).is_err());

        assert_eq!(
            TranspositionFree(&mut heap, Some(&mut transposition)),
            Ok(())
        );
        assert_eq!(TranspositionFree(&mut heap, None), Ok(()));
    }

    #[test]
    fn source_port__ichican2__ctablecreate__line_769() {
        fn full_data(heap: &mut SourceHeap) -> CANON_DATA {
            CANON_DATA {
                nMaxLenLinearCT: 4,
                nLenCTAtOnly: 3,
                NumH: heap.allocate(vec![1_i16]).unwrap(),
                maxlenNumH: 2,
                NumHfixed: heap.allocate(vec![1_i16]).unwrap(),
                maxlenNumHfixed: 1,
                maxlen_iso_sort_key: 2,
                iso_exchg_atnos: heap.allocate(vec![1_i8]).unwrap(),
                maxlen_iso_exchg_atnos: 4,
                ..CANON_DATA::default()
            }
        }

        let mut heap = SourceHeap::default();
        let data = full_data(&mut heap);
        let old = heap.allocate(vec![91_u16]).unwrap();
        let mut table = ConTable {
            Ctbl: old,
            lenCt: 88,
            ..ConTable::default()
        };
        assert_eq!(CTableCreate(&mut heap, &mut table, 3, &data), Ok(1));
        assert_eq!(table.maxVert, 3);
        assert_eq!(table.maxPos, 4);
        assert_eq!(table.nLenCTAtOnly, 3);
        assert_eq!(table.maxlenCt, 5);
        assert_eq!(table.maxlenNumH, 3);
        assert_eq!(table.maxlen_iso_sort_key, 3);
        assert_eq!(table.maxlen_iso_exchg_atnos, 3);
        assert_eq!((table.lenCt, table.lenNumH, table.lenPos), (0, 0, 0));
        assert_eq!(heap.slice(table.Ctbl.as_const()).unwrap(), &[0; 5]);
        assert_eq!(heap.slice(table.nextCtblPos.as_const()).unwrap(), &[0; 4]);
        assert_eq!(heap.slice(table.nextAtRank.as_const()).unwrap(), &[0; 4]);
        assert_eq!(heap.slice(table.NumH.as_const()).unwrap(), &[0; 3]);
        assert_eq!(heap.slice(table.NumHfixed.as_const()).unwrap(), &[0; 2]);
        assert_eq!(heap.slice(table.iso_sort_key.as_const()).unwrap(), &[0; 3]);
        assert_eq!(
            heap.slice(table.iso_exchg_atnos.as_const()).unwrap(),
            &[0; 5]
        );
        assert_eq!(heap.slice(old.as_const()).unwrap(), &[91]);

        let expected = [0, 0, 1, 0, 0, 1, 1, 1];
        for (failure_ordinal, expected_return) in expected.into_iter().enumerate() {
            let mut failure_heap = SourceHeap::default();
            let failure_data = full_data(&mut failure_heap);
            failure_heap.fail_after_allocations(failure_ordinal as u64);
            let mut partial = ConTable::default();
            assert_eq!(
                CTableCreate(&mut failure_heap, &mut partial, 3, &failure_data),
                Ok(expected_return),
                "allocation ordinal {failure_ordinal}"
            );
            assert_eq!(
                failure_heap.source_allocation_calls(),
                7,
                "allocation ordinal {failure_ordinal}"
            );
            if failure_ordinal >= 1 {
                assert!(!partial.Ctbl.is_null());
                assert_eq!(
                    failure_heap.slice(partial.Ctbl.as_const()).unwrap(),
                    &[0; 5]
                );
            }
            if failure_ordinal == 2 {
                assert!(!partial.nextCtblPos.is_null());
                assert!(partial.nextAtRank.is_null());
                assert!(!partial.NumH.is_null());
                assert!(!partial.NumHfixed.is_null());
            }
        }

        let mut omission_heap = SourceHeap::default();
        omission_heap.fail_after_allocations(2);
        let mut omission = ConTable::default();
        assert_eq!(
            CTableCreate(
                &mut omission_heap,
                &mut omission,
                2,
                &CANON_DATA {
                    nMaxLenLinearCT: 1,
                    ..CANON_DATA::default()
                },
            ),
            Ok(1)
        );
        assert!(omission.nextAtRank.is_null());
    }

    #[test]
    fn source_port__ichican2__ctablefree__line_851() {
        let mut heap = SourceHeap::default();
        let mut table = ConTable {
            Ctbl: heap.allocate(vec![1_u16]).unwrap(),
            nextCtblPos: heap.allocate(vec![2_u16]).unwrap(),
            nextAtRank: heap.allocate(vec![3_u16]).unwrap(),
            NumH: heap.allocate(vec![4_i16]).unwrap(),
            NumHfixed: heap.allocate(vec![5_i16]).unwrap(),
            iso_sort_key: heap.allocate(vec![6_i64]).unwrap(),
            iso_exchg_atnos: heap.allocate(vec![7_i8]).unwrap(),
            lenCt: 9,
            maxVert: 10,
            ..ConTable::default()
        };
        let pointers = (
            table.Ctbl,
            table.nextCtblPos,
            table.nextAtRank,
            table.NumH,
            table.NumHfixed,
            table.iso_sort_key,
            table.iso_exchg_atnos,
        );

        assert_eq!(CTableFree(&mut heap, Some(&mut table)), Ok(()));
        assert_eq!(table, ConTable::default());
        assert!(heap.slice(pointers.0.as_const()).is_err());
        assert!(heap.slice(pointers.1.as_const()).is_err());
        assert!(heap.slice(pointers.2.as_const()).is_err());
        assert!(heap.slice(pointers.3.as_const()).is_err());
        assert!(heap.slice(pointers.4.as_const()).is_err());
        assert!(heap.slice(pointers.5.as_const()).is_err());
        assert!(heap.slice(pointers.6.as_const()).is_err());

        assert_eq!(CTableFree(&mut heap, Some(&mut table)), Ok(()));
        assert_eq!(CTableFree(&mut heap, None), Ok(()));

        let partial_pointer = heap.allocate(vec![11_i64]).unwrap();
        let mut partial = ConTable {
            iso_sort_key: partial_pointer,
            len_iso_sort_key: 1,
            ..ConTable::default()
        };
        assert_eq!(CTableFree(&mut heap, Some(&mut partial)), Ok(()));
        assert_eq!(partial, ConTable::default());
        assert!(heap.slice(partial_pointer.as_const()).is_err());
    }

    #[test]
    fn source_port__ichican2__unorderedpartitioncreate__line_904() {
        let mut heap = SourceHeap::default();
        let old = heap.allocate(vec![8_u16, 9]).unwrap();
        let mut partition = UnorderedPartition { equ2: old };
        assert_eq!(
            UnorderedPartitionCreate(&mut heap, &mut partition, 4),
            Ok(1)
        );
        assert_ne!(partition.equ2, old);
        assert_eq!(heap.slice(partition.equ2.as_const()).unwrap(), &[0; 4]);
        assert_eq!(heap.slice(old.as_const()).unwrap(), &[8, 9]);

        let allocated = partition.equ2;
        heap.fail_after_allocations(0);
        assert_eq!(
            UnorderedPartitionCreate(&mut heap, &mut partition, 3),
            Ok(0)
        );
        assert!(partition.equ2.is_null());
        assert_eq!(heap.slice(allocated.as_const()).unwrap(), &[0; 4]);

        let mut negative_heap = SourceHeap::default();
        let mut negative = UnorderedPartition::default();
        assert_eq!(
            UnorderedPartitionCreate(&mut negative_heap, &mut negative, -1),
            Ok(0)
        );
        assert!(negative.equ2.is_null());

        let mut zero_heap = SourceHeap::default();
        let mut zero = UnorderedPartition::default();
        assert_eq!(
            UnorderedPartitionCreate(&mut zero_heap, &mut zero, 0),
            Ok(1)
        );
        assert_eq!(
            zero_heap.slice(zero.equ2.as_const()).unwrap(),
            &[] as &[AT_NUMB]
        );
    }

    #[test]
    fn source_port__ichican2__unorderedpartitionfree__line_919() {
        let mut heap = SourceHeap::default();
        let allocation = heap.allocate(vec![1_u16, 2, 3]).unwrap();
        let mut partition = UnorderedPartition { equ2: allocation };
        assert_eq!(UnorderedPartitionFree(&mut heap, &mut partition), Ok(()));
        assert!(partition.equ2.is_null());
        assert!(heap.slice(allocation.as_const()).is_err());
        assert_eq!(UnorderedPartitionFree(&mut heap, &mut partition), Ok(()));

        let mut empty = UnorderedPartition::default();
        assert_eq!(UnorderedPartitionFree(&mut heap, &mut empty), Ok(()));
        assert!(empty.equ2.is_null());
    }

    #[test]
    fn source_port__ichican2__unorderedpartitionmakediscrete__line_933() {
        let mut heap = SourceHeap::default();
        let values = heap.allocate(vec![99_u16; 5]).unwrap();
        let partition = UnorderedPartition { equ2: values };
        assert_eq!(
            UnorderedPartitionMakeDiscrete(&mut heap, &partition, 5),
            Ok(())
        );
        assert_eq!(heap.slice(values.as_const()).unwrap(), &[0, 1, 2, 3, 4]);

        heap.slice_mut(values).unwrap().fill(77);
        assert_eq!(
            UnorderedPartitionMakeDiscrete(&mut heap, &partition, 0),
            Ok(())
        );
        assert_eq!(heap.slice(values.as_const()).unwrap(), &[77; 5]);
        assert_eq!(
            UnorderedPartitionMakeDiscrete(&mut heap, &partition, -1),
            Ok(())
        );
        assert_eq!(heap.slice(values.as_const()).unwrap(), &[77; 5]);

        let short = heap.allocate(vec![88_u16; 2]).unwrap();
        assert_eq!(
            UnorderedPartitionMakeDiscrete(&mut heap, &UnorderedPartition { equ2: short }, 3,),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(short.as_const()).unwrap(), &[0, 1]);

        let boundary = heap.allocate(vec![0_u16; 65_537]).unwrap();
        assert_eq!(
            UnorderedPartitionMakeDiscrete(
                &mut heap,
                &UnorderedPartition { equ2: boundary },
                65_537,
            ),
            Ok(())
        );
        assert_eq!(heap.slice(boundary.as_const()).unwrap()[65_535], u16::MAX);
        assert_eq!(heap.slice(boundary.as_const()).unwrap()[65_536], 0);
    }

    #[test]
    fn source_port__ichican2__partitioncreate__line_946() {
        let mut heap = SourceHeap::default();
        let old_atoms = heap.allocate(vec![8_u16]).unwrap();
        let old_ranks = heap.allocate(vec![9_u16]).unwrap();
        let mut partition = Partition {
            AtNumber: old_atoms,
            Rank: old_ranks,
        };
        assert_eq!(PartitionCreate(&mut heap, &mut partition, 4), Ok(1));
        assert_eq!(heap.slice(partition.AtNumber.as_const()).unwrap(), &[0; 4]);
        assert_eq!(heap.slice(partition.Rank.as_const()).unwrap(), &[0; 4]);
        assert_eq!(heap.slice(old_atoms.as_const()).unwrap(), &[8]);
        assert_eq!(heap.slice(old_ranks.as_const()).unwrap(), &[9]);

        for failure_ordinal in 0..=2_u64 {
            let mut failure_heap = SourceHeap::default();
            failure_heap.fail_after_allocations(failure_ordinal);
            let mut partial = Partition::default();
            assert_eq!(
                PartitionCreate(&mut failure_heap, &mut partial, 3),
                Ok(if failure_ordinal == 2 { 1 } else { 0 })
            );
            assert_eq!(failure_heap.source_allocation_calls(), 2);
            if failure_ordinal == 0 {
                assert!(partial.AtNumber.is_null());
                assert!(!partial.Rank.is_null());
            } else if failure_ordinal == 1 {
                assert!(!partial.AtNumber.is_null());
                assert!(partial.Rank.is_null());
            }
        }

        let mut negative_heap = SourceHeap::default();
        let mut negative = Partition::default();
        assert_eq!(
            PartitionCreate(&mut negative_heap, &mut negative, -1),
            Ok(0)
        );
        assert!(negative.AtNumber.is_null());
        assert!(negative.Rank.is_null());

        let mut zero_heap = SourceHeap::default();
        let mut zero = Partition::default();
        assert_eq!(PartitionCreate(&mut zero_heap, &mut zero, 0), Ok(1));
        assert_eq!(
            zero_heap.slice(zero.AtNumber.as_const()).unwrap(),
            &[] as &[AT_NUMB]
        );
        assert_eq!(
            zero_heap.slice(zero.Rank.as_const()).unwrap(),
            &[] as &[AT_RANK]
        );
    }

    #[test]
    fn source_port__ichican2__partitionisdiscrete__line_978() {
        let mut heap = SourceHeap::default();
        let atoms = heap.allocate(vec![2_u16, 0, 1]).unwrap();
        let ranks = heap.allocate(vec![2_u16, 3, 1]).unwrap();
        let partition = Partition {
            AtNumber: atoms,
            Rank: ranks,
        };
        assert_eq!(PartitionIsDiscrete(&heap, &partition, 3), Ok(1));

        heap.slice_mut(ranks).unwrap()[0] |= rank_mark_bit();
        assert_eq!(PartitionIsDiscrete(&heap, &partition, 3), Ok(1));

        heap.slice_mut(ranks).unwrap()[1] = 2;
        assert_eq!(PartitionIsDiscrete(&heap, &partition, 3), Ok(0));
        heap.slice_mut(ranks).unwrap()[1] = 3;
        heap.slice_mut(ranks).unwrap()[2] = 9;
        assert_eq!(PartitionIsDiscrete(&heap, &partition, 3), Ok(0));

        assert_eq!(PartitionIsDiscrete(&heap, &Partition::default(), 0), Ok(1));
        assert_eq!(PartitionIsDiscrete(&heap, &Partition::default(), -1), Ok(1));
        assert_eq!(
            PartitionIsDiscrete(&heap, &Partition::default(), 1),
            Err(SourceHeapError::NullPointer)
        );

        let bad_atoms = heap.allocate(vec![7_u16]).unwrap();
        assert_eq!(
            PartitionIsDiscrete(
                &heap,
                &Partition {
                    AtNumber: bad_atoms,
                    Rank: ranks,
                },
                1,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichican2__partitiongetfirstcell__line_998() {
        let mut heap = SourceHeap::default();
        let mut globals = CANON_GLOBALS::default();
        assert_eq!(SetBitCreate(&mut heap, &mut globals), Ok(1));
        let atoms = heap.allocate(vec![0_u16, 1, 2, 3, 4, 5]).unwrap();
        let ranks = heap.allocate(vec![1_u16, 2, 4, 4, 6, 6]).unwrap();
        let partition = Partition {
            AtNumber: atoms,
            Rank: ranks,
        };
        let mut cells = vec![
            Cell {
                first: 90,
                next: 91,
                prev: 77,
            },
            Cell {
                first: 92,
                next: 93,
                prev: 78,
            },
        ];
        assert_eq!(
            PartitionGetFirstCell(&heap, &partition, &mut cells, 1, 6),
            Ok(2)
        );
        assert_eq!(
            cells[0],
            Cell {
                first: 2,
                next: 4,
                prev: 77
            }
        );
        assert_eq!(
            PartitionGetFirstCell(&heap, &partition, &mut cells, 2, 6),
            Ok(2)
        );
        assert_eq!(
            cells[1],
            Cell {
                first: 4,
                next: 6,
                prev: 78
            }
        );

        heap.slice_mut(ranks).unwrap()[2] |= rank_mark_bit();
        heap.slice_mut(ranks).unwrap()[3] |= rank_mark_bit();
        assert_eq!(
            PartitionGetFirstCell(&heap, &partition, &mut cells, 1, 6),
            Ok(2)
        );

        let discrete_ranks = heap.allocate(vec![1_u16, 2, 3]).unwrap();
        let discrete_atoms = heap.allocate(vec![0_u16, 1, 2]).unwrap();
        let mut no_cell = [Cell {
            first: 10,
            next: 11,
            prev: 12,
        }];
        assert_eq!(
            PartitionGetFirstCell(
                &heap,
                &Partition {
                    AtNumber: discrete_atoms,
                    Rank: discrete_ranks,
                },
                &mut no_cell,
                1,
                3,
            ),
            Ok(0)
        );
        assert_eq!(
            no_cell[0],
            Cell {
                first: INCHI_CANON_INFINITY as i32,
                next: 0,
                prev: 12,
            }
        );

        assert_eq!(
            PartitionGetFirstCell(&heap, &partition, &mut [], 1, 6),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let bad_atoms = heap.allocate(vec![0_u16, 9]).unwrap();
        let bad_ranks = heap.allocate(vec![2_u16]).unwrap();
        let mut partial = [Cell {
            first: 40,
            next: 55,
            prev: 60,
        }];
        assert_eq!(
            PartitionGetFirstCell(
                &heap,
                &Partition {
                    AtNumber: bad_atoms,
                    Rank: bad_ranks,
                },
                &mut partial,
                1,
                2,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            partial[0],
            Cell {
                first: 0,
                next: 55,
                prev: 60
            }
        );
    }

    #[test]
    fn source_port__ichican2__canongraph__line_3662() {
        fn fixture(
            heap: &mut SourceHeap,
            ranks: &[AT_RANK],
        ) -> (
            SourceMutPointer<NEIGH_LIST>,
            Vec<Partition>,
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_NUMB>,
            CANON_DATA,
        ) {
            let count = ranks.len();
            let rank = heap.allocate(ranks.to_vec()).unwrap();
            let at_number = heap
                .allocate((0..count).map(|index| index as AT_NUMB).collect())
                .unwrap();
            let mut graph_rows = Vec::with_capacity(count);
            for _ in 0..count {
                graph_rows.push(heap.allocate(vec![0_u16]).unwrap());
            }
            let graph = heap.allocate(graph_rows).unwrap();
            let mut partitions = vec![Partition::default(); count + 3];
            partitions[0] = Partition {
                Rank: rank,
                AtNumber: at_number,
            };
            let symmetry = heap.allocate(vec![77_u16; count]).unwrap();
            let canonical = heap.allocate(vec![88_u16; count]).unwrap();
            let atom_order = heap.allocate(vec![99_u16; count]).unwrap();
            let linear = heap.allocate(vec![66_u16; count * 2 + 2]).unwrap();
            (
                graph,
                partitions,
                symmetry,
                canonical,
                atom_order,
                CANON_DATA {
                    LinearCT: linear,
                    nMaxLenLinearCT: (count * 2 + 1) as i32,
                    nLenCTAtOnly: count as i32,
                    ..CANON_DATA::default()
                },
            )
        }

        let mut heap = SourceHeap::default();
        let (graph, mut pi, symmetry, canonical, atom_order, mut data) = fixture(&mut heap, &[1]);
        let mut globals = CANON_GLOBALS::default();
        let mut counts = CANON_COUNTS::default();
        let mut rho_output = SourceMutPointer::null();
        assert_eq!(
            CanonGraph(
                &mut heap,
                &mut INCHI_CLOCK::default(),
                &mut globals,
                1,
                1,
                1,
                0,
                graph,
                &mut pi,
                symmetry,
                canonical,
                atom_order,
                &mut data,
                &mut counts,
                SourceMutPointer::null(),
                Some(&mut rho_output),
                0,
                None,
                None,
                0,
            ),
            Ok(0)
        );
        assert_eq!(heap.slice(symmetry.as_const()).unwrap(), &[1]);
        assert_eq!(heap.slice(canonical.as_const()).unwrap(), &[1]);
        assert_eq!(heap.slice(atom_order.as_const()).unwrap(), &[0]);
        assert_eq!(data.nLenLinearCT, 0);
        assert_eq!(heap.slice(data.LinearCT.as_const()).unwrap()[0], 66);
        assert_eq!(counts.lNumTotCT, 1);
        assert_eq!(counts.dGroupSize, 1.0);
        assert!(!rho_output.is_null());
        let fixed = source_clone(&heap, rho_output, 0).unwrap();
        assert_eq!(fixed.lenCt, 1);
        assert_eq!(&heap.slice(fixed.Ctbl.as_const()).unwrap()[..2], &[1, 0]);

        let (graph2, mut pi2, symmetry2, canonical2, atom_order2, mut data2) =
            fixture(&mut heap, &[1]);
        let mut fixed_counts = CANON_COUNTS::default();
        assert_eq!(
            CanonGraph(
                &mut heap,
                &mut INCHI_CLOCK::default(),
                &mut globals,
                1,
                1,
                1,
                0,
                graph2,
                &mut pi2,
                symmetry2,
                canonical2,
                atom_order2,
                &mut data2,
                &mut fixed_counts,
                rho_output,
                None,
                0,
                None,
                None,
                0,
            ),
            Ok(0)
        );
        assert_eq!(heap.slice(symmetry2.as_const()).unwrap(), &[1]);

        fn quit() -> u32 {
            USER_ACTION_QUIT
        }
        let mut quit_heap = SourceHeap::default();
        let (quit_graph, mut quit_pi, quit_symmetry, quit_canonical, quit_atoms, mut quit_data) =
            fixture(&mut quit_heap, &[2, 2]);
        assert_eq!(
            CanonGraph(
                &mut quit_heap,
                &mut INCHI_CLOCK::default(),
                &mut CANON_GLOBALS::default(),
                2,
                2,
                2,
                0,
                quit_graph,
                &mut quit_pi,
                quit_symmetry,
                quit_canonical,
                quit_atoms,
                &mut quit_data,
                &mut CANON_COUNTS::default(),
                SourceMutPointer::null(),
                None,
                0,
                Some(quit),
                None,
                0,
            ),
            Ok(CT_USER_QUIT_ERR)
        );

        let mut timeout_heap = SourceHeap::default();
        let (
            timeout_graph,
            mut timeout_pi,
            timeout_symmetry,
            timeout_canonical,
            timeout_atoms,
            mut timeout_data,
        ) = fixture(&mut timeout_heap, &[2, 2]);
        timeout_data.ulTimeOutTime = timeout_heap
            .allocate(vec![inchiTime { clockTime: 0 }])
            .unwrap();
        assert_eq!(
            CanonGraph(
                &mut timeout_heap,
                &mut INCHI_CLOCK::default(),
                &mut CANON_GLOBALS::default(),
                2,
                2,
                2,
                0,
                timeout_graph,
                &mut timeout_pi,
                timeout_symmetry,
                timeout_canonical,
                timeout_atoms,
                &mut timeout_data,
                &mut CANON_COUNTS::default(),
                SourceMutPointer::null(),
                None,
                0,
                None,
                None,
                1,
            ),
            Ok(CT_TIMEOUT_ERR)
        );

        let mut search_heap = SourceHeap::default();
        let (
            search_graph,
            mut search_pi,
            search_symmetry,
            search_canonical,
            search_atoms,
            mut search_data,
        ) = fixture(&mut search_heap, &[2, 2]);
        let mut search_counts = CANON_COUNTS::default();
        assert_eq!(
            CanonGraph(
                &mut search_heap,
                &mut INCHI_CLOCK::default(),
                &mut CANON_GLOBALS::default(),
                2,
                2,
                2,
                0,
                search_graph,
                &mut search_pi,
                search_symmetry,
                search_canonical,
                search_atoms,
                &mut search_data,
                &mut search_counts,
                SourceMutPointer::null(),
                None,
                0,
                None,
                None,
                0,
            ),
            Ok(0)
        );
        assert_eq!(
            search_heap.slice(search_symmetry.as_const()).unwrap(),
            &[1, 1]
        );
        assert_eq!(search_counts.dGroupSize, 2.0);
        assert_eq!(search_counts.lNumGenerators, 1);

        fn console_quit() -> i32 {
            1
        }
        let mut console_heap = SourceHeap::default();
        let (
            console_graph,
            mut console_pi,
            console_symmetry,
            console_canonical,
            console_atoms,
            mut console_data,
        ) = fixture(&mut console_heap, &[2, 2]);
        assert_eq!(
            CanonGraph(
                &mut console_heap,
                &mut INCHI_CLOCK::default(),
                &mut CANON_GLOBALS::default(),
                2,
                2,
                2,
                0,
                console_graph,
                &mut console_pi,
                console_symmetry,
                console_canonical,
                console_atoms,
                &mut console_data,
                &mut CANON_COUNTS::default(),
                SourceMutPointer::null(),
                None,
                0,
                None,
                Some(console_quit),
                0,
            ),
            Ok(CT_USER_QUIT_ERR)
        );

        let mut mismatch_heap = SourceHeap::default();
        let (
            mismatch_graph,
            mut mismatch_pi,
            mismatch_symmetry,
            mismatch_canonical,
            mismatch_atoms,
            mut mismatch_data,
        ) = fixture(&mut mismatch_heap, &[1]);
        let mismatch_fixed = mismatch_heap
            .allocate(vec![ConTable {
                nLenCTAtOnly: mismatch_data.nLenCTAtOnly + 1,
                ..ConTable::default()
            }])
            .unwrap();
        assert_eq!(
            CanonGraph(
                &mut mismatch_heap,
                &mut INCHI_CLOCK::default(),
                &mut CANON_GLOBALS::default(),
                1,
                1,
                1,
                0,
                mismatch_graph,
                &mut mismatch_pi,
                mismatch_symmetry,
                mismatch_canonical,
                mismatch_atoms,
                &mut mismatch_data,
                &mut CANON_COUNTS::default(),
                mismatch_fixed,
                None,
                0,
                None,
                None,
                0,
            ),
            Ok(-2)
        );

        for failure_ordinal in 0_u64..27 {
            let mut failure_heap = SourceHeap::default();
            let (
                failure_graph,
                mut failure_pi,
                failure_symmetry,
                failure_canonical,
                failure_atoms,
                mut failure_data,
            ) = fixture(&mut failure_heap, &[1]);
            let baseline = failure_heap.live_allocation_count();
            failure_heap.fail_after_allocations(failure_ordinal);
            let mut failure_globals = CANON_GLOBALS::default();
            let result = CanonGraph(
                &mut failure_heap,
                &mut INCHI_CLOCK::default(),
                &mut failure_globals,
                1,
                1,
                1,
                0,
                failure_graph,
                &mut failure_pi,
                failure_symmetry,
                failure_canonical,
                failure_atoms,
                &mut failure_data,
                &mut CANON_COUNTS::default(),
                SourceMutPointer::null(),
                None,
                0,
                None,
                None,
                0,
            );
            let expected = match failure_ordinal {
                0 => Ok(-1),
                _ => Ok(CT_OUT_OF_RAM),
            };
            assert_eq!(result, expected, "allocation ordinal {failure_ordinal}");
            assert_eq!(
                failure_heap.live_allocation_count(),
                baseline + usize::from(!failure_globals.m_bBit.is_null()),
                "cleanup allocation ordinal {failure_ordinal}"
            );
            assert_eq!(
                failure_heap.source_allocation_calls(),
                match failure_ordinal {
                    0 => 1,
                    12 => 24,
                    16 | 18 | 20 => 26,
                    _ => 27,
                },
                "allocation call count ordinal {failure_ordinal}"
            );
        }

        let mut owned_fixed = source_clone(&heap, rho_output, 0).unwrap();
        CTableFree(&mut heap, Some(&mut owned_fixed)).unwrap();
        heap.free(rho_output).unwrap();
    }

    #[test]
    fn source_port__ichican2__ctpartcompare__line_2123() {
        fn make_ct(
            heap: &mut SourceHeap,
            ctbl: &[u16],
            h: &[i16],
            fixed_h: &[i16],
            iso: &[i64],
            exchange: &[i8],
            atoms_only: i32,
        ) -> ConTable {
            let len_ct = ctbl.len() as i32;
            let atom_len = h.len() as i32;
            ConTable {
                Ctbl: heap.allocate(ctbl.to_vec()).unwrap(),
                lenCt: len_ct,
                nLenCTAtOnly: atoms_only,
                maxlenCt: len_ct,
                maxVert: atom_len,
                lenPos: 1,
                nextCtblPos: heap.allocate(vec![len_ct as u16]).unwrap(),
                nextAtRank: heap.allocate(vec![(atom_len + 1) as u16]).unwrap(),
                NumH: heap.allocate(h.to_vec()).unwrap(),
                lenNumH: atom_len,
                NumHfixed: heap.allocate(fixed_h.to_vec()).unwrap(),
                iso_sort_key: heap.allocate(iso.to_vec()).unwrap(),
                len_iso_sort_key: atom_len,
                iso_exchg_atnos: heap.allocate(exchange.to_vec()).unwrap(),
                len_iso_exchg_atnos: atom_len,
                ..ConTable::default()
            }
        }

        let mut heap = SourceHeap::default();
        let base = make_ct(
            &mut heap,
            &[1, 2, 3, 4],
            &[10, 11],
            &[20, 21],
            &[30, 31],
            &[1, 2],
            2,
        );
        let equal = make_ct(
            &mut heap,
            &[1, 2, 3, 4],
            &[10, 11],
            &[20, 21],
            &[30, 31],
            &[1, 2],
            2,
        );
        assert_eq!(
            CtPartCompare(
                &mut heap,
                &base,
                &equal,
                SourceMutPointer::null(),
                None,
                1,
                0,
                1,
            ),
            Ok(0)
        );

        let layer0 = make_ct(
            &mut heap,
            &[2, 2, 3, 4],
            &[10, 11],
            &[20, 21],
            &[30, 31],
            &[1, 2],
            2,
        );
        assert_eq!(
            CtPartCompare(
                &mut heap,
                &base,
                &layer0,
                SourceMutPointer::null(),
                None,
                1,
                0,
                1,
            ),
            Ok(1)
        );

        let layer1 = make_ct(
            &mut heap,
            &[1, 2, 3, 4],
            &[11, 11],
            &[20, 21],
            &[30, 31],
            &[1, 2],
            2,
        );
        assert_eq!(
            CtPartCompare(
                &mut heap,
                &base,
                &layer1,
                SourceMutPointer::null(),
                None,
                1,
                0,
                1,
            ),
            Ok(2)
        );

        let layer2 = make_ct(
            &mut heap,
            &[1, 2, 4, 4],
            &[10, 11],
            &[20, 21],
            &[30, 31],
            &[1, 2],
            2,
        );
        assert_eq!(
            CtPartCompare(
                &mut heap,
                &base,
                &layer2,
                SourceMutPointer::null(),
                None,
                1,
                0,
                1,
            ),
            Ok(3)
        );
        assert_eq!(
            CtPartCompare(
                &mut heap,
                &base,
                &layer2,
                SourceMutPointer::null(),
                None,
                1,
                0,
                0,
            ),
            Ok(1)
        );

        let layer3 = make_ct(
            &mut heap,
            &[1, 2, 3, 4],
            &[10, 11],
            &[21, 21],
            &[30, 31],
            &[1, 2],
            2,
        );
        assert_eq!(
            CtPartCompare(
                &mut heap,
                &base,
                &layer3,
                SourceMutPointer::null(),
                None,
                1,
                0,
                1,
            ),
            Ok(4)
        );

        let layer4 = make_ct(
            &mut heap,
            &[1, 2, 3, 4],
            &[10, 11],
            &[20, 21],
            &[31, 31],
            &[1, 2],
            2,
        );
        assert_eq!(
            CtPartCompare(
                &mut heap,
                &base,
                &layer4,
                SourceMutPointer::null(),
                None,
                1,
                0,
                1,
            ),
            Ok(5)
        );
        let exchange4 = make_ct(
            &mut heap,
            &[1, 2, 3, 4],
            &[10, 11],
            &[20, 21],
            &[30, 31],
            &[2, 2],
            2,
        );
        assert_eq!(
            CtPartCompare(
                &mut heap,
                &base,
                &exchange4,
                SourceMutPointer::null(),
                None,
                1,
                0,
                1,
            ),
            Ok(5)
        );

        let mut least: KLeastLayers = std::array::from_fn(|_| kLeast::default());
        assert_eq!(
            CtPartCompare(
                &mut heap,
                &base,
                &layer1,
                SourceMutPointer::null(),
                Some(&mut least),
                1,
                0,
                1,
            ),
            Ok(0)
        );
        assert_eq!(least[1], kLeast { k: 1, i: 0 });

        let cmp = heap.allocate(vec![1_i8, 0]).unwrap();
        assert_eq!(
            CtPartCompare(
                &mut heap,
                &ConTable::default(),
                &ConTable::default(),
                cmp,
                None,
                2,
                0,
                1,
            ),
            Ok(1)
        );
        assert_eq!(heap.slice(cmp.as_const()).unwrap(), &[1, 1]);

        let common_short = make_ct(
            &mut heap,
            &[1, 2, 3],
            &[10, 11],
            &[20, 21],
            &[30, 31],
            &[1, 2],
            2,
        );
        assert_eq!(
            CtPartCompare(
                &mut heap,
                &base,
                &common_short,
                SourceMutPointer::null(),
                None,
                1,
                1,
                1,
            ),
            Ok(0)
        );
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

    #[test]
    fn source_port__ichican2__ctpartcopy__line_2965() {
        let mut heap = SourceHeap::default();
        let source_ctbl = heap
            .allocate(vec![10_u16, 11, 12, 13, 14, 15, 16, 17])
            .unwrap();
        let source_positions = heap.allocate(vec![2_u16, 5, 7]).unwrap();
        let source_ranks = heap.allocate(vec![3_u16, 6, 9]).unwrap();
        let source_h = heap
            .allocate(vec![100_i16, 101, 102, 103, 104, 105, 106])
            .unwrap();
        let source_fixed_h = heap
            .allocate(vec![200_i16, 201, 202, 203, 204, 205, 206])
            .unwrap();
        let source_iso = heap
            .allocate(vec![300_i64, 301, 302, 303, 304, 305, 306])
            .unwrap();
        let source_exchange = heap.allocate(vec![10_i8, 11, 12, 13, 14, 15, 16]).unwrap();
        let source = ConTable {
            Ctbl: source_ctbl,
            nextCtblPos: source_positions,
            nextAtRank: source_ranks,
            NumH: source_h,
            lenNumH: 7,
            maxVert: 4,
            NumHfixed: source_fixed_h,
            iso_sort_key: source_iso,
            iso_exchg_atnos: source_exchange,
            lenPos: 3,
            ..ConTable::default()
        };

        let target_ctbl = heap.allocate(vec![999_u16; 8]).unwrap();
        let target_positions = heap.allocate(vec![77_u16; 3]).unwrap();
        let target_ranks = heap.allocate(vec![88_u16; 3]).unwrap();
        let target_h = heap.allocate(vec![-1_i16; 8]).unwrap();
        let target_fixed_h = heap.allocate(vec![-2_i16; 8]).unwrap();
        let target_iso = heap.allocate(vec![-3_i64; 8]).unwrap();
        let target_exchange = heap.allocate(vec![-4_i8; 8]).unwrap();
        let mut target = ConTable {
            Ctbl: target_ctbl,
            nextCtblPos: target_positions,
            nextAtRank: target_ranks,
            NumH: target_h,
            lenNumH: 91,
            NumHfixed: target_fixed_h,
            iso_sort_key: target_iso,
            len_iso_sort_key: 92,
            iso_exchg_atnos: target_exchange,
            len_iso_exchg_atnos: 93,
            lenPos: 94,
            ..ConTable::default()
        };

        assert_eq!(CtPartCopy(&mut heap, &mut target, &source, 1), Ok(()));
        assert_eq!(
            heap.slice(target_ctbl.as_const()).unwrap(),
            &[10, 11, 999, 999, 999, 999, 999, 999]
        );
        assert_eq!(
            heap.slice(target_positions.as_const()).unwrap(),
            &[2, 77, 77]
        );
        assert_eq!(heap.slice(target_ranks.as_const()).unwrap(), &[3, 88, 88]);
        assert_eq!(
            heap.slice(target_h.as_const()).unwrap(),
            &[100, 101, -1, -1, -1, -1, -1, -1]
        );
        assert_eq!(
            heap.slice(target_fixed_h.as_const()).unwrap(),
            &[200, 201, -2, -2, -2, -2, -2, -2]
        );
        assert_eq!(
            heap.slice(target_iso.as_const()).unwrap(),
            &[300, 301, -3, -3, -3, -3, -3, -3]
        );
        assert_eq!(
            heap.slice(target_exchange.as_const()).unwrap(),
            &[10, 11, -4, -4, -4, -4, -4, -4]
        );
        assert_eq!(target.lenCt, 2);
        assert_eq!(target.lenNumH, 2);
        assert_eq!(target.len_iso_sort_key, 2);
        assert_eq!(target.len_iso_exchg_atnos, 2);
        assert_eq!(target.lenPos, 1);

        heap.slice_mut(target_positions).unwrap()[0] = 1;
        heap.slice_mut(target_ranks).unwrap()[0] = 2;
        assert_eq!(CtPartCopy(&mut heap, &mut target, &source, 2), Ok(()));
        assert_eq!(
            heap.slice(target_ctbl.as_const()).unwrap(),
            &[10, 12, 13, 14, 999, 999, 999, 999]
        );
        assert_eq!(
            heap.slice(target_positions.as_const()).unwrap(),
            &[1, 4, 77]
        );
        assert_eq!(heap.slice(target_ranks.as_const()).unwrap(), &[2, 6, 88]);
        assert_eq!(
            heap.slice(target_h.as_const()).unwrap(),
            &[100, 102, 103, 104, 105, 106, -1, -1]
        );
        assert_eq!(
            heap.slice(target_fixed_h.as_const()).unwrap(),
            &[200, 202, 203, 204, -2, -2, -2, -2]
        );
        assert_eq!(
            heap.slice(target_iso.as_const()).unwrap(),
            &[300, 302, 303, 304, -3, -3, -3, -3]
        );
        assert_eq!(
            heap.slice(target_exchange.as_const()).unwrap(),
            &[10, 12, 13, 14, -4, -4, -4, -4]
        );
        assert_eq!(target.lenCt, 4);
        assert_eq!(target.lenNumH, 6);
        assert_eq!(target.len_iso_sort_key, 4);
        assert_eq!(target.len_iso_exchg_atnos, 4);
        assert_eq!(target.lenPos, 2);

        let empty_source_ctbl = heap.allocate(vec![7_u16]).unwrap();
        let empty_source_positions = heap.allocate(vec![0_u16]).unwrap();
        let empty_source_ranks = heap.allocate(vec![1_u16]).unwrap();
        let empty_source = ConTable {
            Ctbl: empty_source_ctbl,
            nextCtblPos: empty_source_positions,
            nextAtRank: empty_source_ranks,
            ..ConTable::default()
        };
        let empty_target_ctbl = heap.allocate(vec![8_u16]).unwrap();
        let empty_target_positions = heap.allocate(vec![9_u16]).unwrap();
        let empty_target_ranks = heap.allocate(vec![10_u16]).unwrap();
        let mut empty_target = ConTable {
            Ctbl: empty_target_ctbl,
            nextCtblPos: empty_target_positions,
            nextAtRank: empty_target_ranks,
            lenNumH: 51,
            len_iso_sort_key: 52,
            len_iso_exchg_atnos: 53,
            lenPos: 54,
            ..ConTable::default()
        };
        assert_eq!(
            CtPartCopy(&mut heap, &mut empty_target, &empty_source, 1),
            Ok(())
        );
        assert_eq!(heap.slice(empty_target_ctbl.as_const()).unwrap(), &[8]);
        assert_eq!(heap.slice(empty_target_positions.as_const()).unwrap(), &[0]);
        assert_eq!(heap.slice(empty_target_ranks.as_const()).unwrap(), &[1]);
        assert_eq!(empty_target.lenCt, 0);
        assert_eq!(empty_target.lenNumH, 51);
        assert_eq!(empty_target.len_iso_sort_key, 52);
        assert_eq!(empty_target.len_iso_exchg_atnos, 53);
        assert_eq!(empty_target.lenPos, 1);
    }

    #[test]
    fn source_port__ichican2__ctfullcopy__line_3113() {
        let mut heap = SourceHeap::default();
        let source_ctbl = heap.allocate(vec![10_u16, 11, 12, 13, 14]).unwrap();
        let source_positions = heap.allocate(vec![2_u16, 5]).unwrap();
        let source_ranks = heap.allocate(vec![3_u16, 6]).unwrap();
        let source_h = heap.allocate(vec![20_i16, 21, 22, 23, 24]).unwrap();
        let source = ConTable {
            Ctbl: source_ctbl,
            nextCtblPos: source_positions,
            nextAtRank: source_ranks,
            NumH: source_h,
            lenNumH: 5,
            maxVert: 8,
            lenPos: 2,
            ..ConTable::default()
        };
        let target_ctbl = heap.allocate(vec![99_u16; 6]).unwrap();
        let target_positions = heap.allocate(vec![98_u16; 2]).unwrap();
        let target_ranks = heap.allocate(vec![97_u16; 2]).unwrap();
        let target_h = heap.allocate(vec![-1_i16; 6]).unwrap();
        let mut target = ConTable {
            Ctbl: target_ctbl,
            nextCtblPos: target_positions,
            nextAtRank: target_ranks,
            NumH: target_h,
            lenCt: 91,
            lenNumH: 92,
            lenPos: 93,
            ..ConTable::default()
        };

        assert_eq!(CtFullCopy(&mut heap, &mut target, &source), Ok(()));
        assert_eq!(
            heap.slice(target_ctbl.as_const()).unwrap(),
            &[10, 11, 12, 13, 14, 99]
        );
        assert_eq!(heap.slice(target_positions.as_const()).unwrap(), &[2, 5]);
        assert_eq!(heap.slice(target_ranks.as_const()).unwrap(), &[3, 6]);
        assert_eq!(
            heap.slice(target_h.as_const()).unwrap(),
            &[20, 21, 22, 23, 24, -1]
        );
        assert_eq!(target.lenCt, 5);
        assert_eq!(target.lenNumH, 5);
        assert_eq!(target.lenPos, 2);

        let empty = ConTable::default();
        let before = target.clone();
        assert_eq!(CtFullCopy(&mut heap, &mut target, &empty), Ok(()));
        assert_eq!(target, before);
    }

    #[test]
    fn source_port__ichican2__cellmakeempty__line_1040() {
        let original = [
            Cell {
                first: 10,
                next: 11,
                prev: 12,
            },
            Cell {
                first: 20,
                next: 21,
                prev: 22,
            },
            Cell {
                first: 30,
                next: 31,
                prev: 32,
            },
        ];

        let mut first = original.clone();
        assert_eq!(CellMakeEmpty(&mut first, 1), Ok(()));
        assert_eq!(first[0].first, INCHI_CANON_INFINITY as i32);
        assert_eq!(first[0].next, 0);
        assert_eq!(first[0].prev, -1);
        assert_eq!(&first[1..], &original[1..]);

        let mut middle = original.clone();
        assert_eq!(CellMakeEmpty(&mut middle, 2), Ok(()));
        assert_eq!(middle[0], original[0]);
        assert_eq!(middle[1].first, INCHI_CANON_INFINITY as i32);
        assert_eq!(middle[1].next, 0);
        assert_eq!(middle[1].prev, -1);
        assert_eq!(middle[2], original[2]);

        let mut zero = original.clone();
        assert_eq!(
            CellMakeEmpty(&mut zero, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(zero, original);

        let mut past_end = original.clone();
        assert_eq!(
            CellMakeEmpty(&mut past_end, 4),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(past_end, original);
    }

    #[test]
    fn source_port__ichican2__partitionfree__line_959() {
        let mut heap = SourceHeap::default();

        assert_eq!(PartitionFree(&mut heap, None), Ok(()));

        let mut empty = Partition::default();
        assert_eq!(PartitionFree(&mut heap, Some(&mut empty)), Ok(()));
        assert!(empty.AtNumber.is_null());
        assert!(empty.Rank.is_null());

        let at_number_only = heap.allocate(vec![2_u16, 0, 1]).unwrap();
        let mut with_at_number = Partition {
            AtNumber: at_number_only,
            ..Partition::default()
        };
        assert_eq!(PartitionFree(&mut heap, Some(&mut with_at_number)), Ok(()));
        assert!(with_at_number.AtNumber.is_null());
        assert!(with_at_number.Rank.is_null());
        assert_eq!(
            heap.slice(at_number_only.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let rank_only = heap.allocate(vec![1_u16, 2, 3]).unwrap();
        let mut with_rank = Partition {
            Rank: rank_only,
            ..Partition::default()
        };
        assert_eq!(PartitionFree(&mut heap, Some(&mut with_rank)), Ok(()));
        assert!(with_rank.AtNumber.is_null());
        assert!(with_rank.Rank.is_null());
        assert_eq!(
            heap.slice(rank_only.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let at_number = heap.allocate(vec![1_u16, 0]).unwrap();
        let rank = heap.allocate(vec![2_u16, 1]).unwrap();
        let mut both = Partition {
            AtNumber: at_number,
            Rank: rank,
        };
        assert_eq!(PartitionFree(&mut heap, Some(&mut both)), Ok(()));
        assert!(both.AtNumber.is_null());
        assert!(both.Rank.is_null());
        assert_eq!(
            heap.slice(at_number.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(rank.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(PartitionFree(&mut heap, Some(&mut both)), Ok(()));
    }

    #[test]
    fn source_port__ichican2__deallocbcn__line_4960() {
        fn assert_freed<T: 'static>(heap: &SourceHeap, pointer: SourceMutPointer<T>) {
            assert!(matches!(
                heap.slice(pointer.as_const()),
                Err(SourceHeapError::MissingAllocation)
            ));
        }

        let mut heap = SourceHeap::default();
        assert_eq!(DeAllocBCN(&mut heap, None), Ok(()));
        assert_eq!(DeAllocBCN(&mut heap, Some(&mut BCN::default())), Ok(()));

        let rank0 = heap.allocate(vec![1_u16, 2]).unwrap();
        let rank1 = heap.allocate(vec![3_u16, 4]).unwrap();
        let rank2 = heap.allocate(vec![5_u16, 6]).unwrap();
        let rank_stack = heap.allocate(vec![rank0, rank1, rank2]).unwrap();
        let timeout = heap.allocate(vec![Default::default()]).unwrap();
        let mut bcn = BCN {
            pRankStack: rank_stack,
            nMaxLenRankStack: 3,
            num_max: 41,
            num_at_tg: 42,
            num_atoms: 43,
            ulTimeOutTime: timeout,
            ..BCN::default()
        };

        let mut neighbor_storage = [SourceMutPointer::null(); TAUT_NUM as usize];
        let mut partition_ct_rank = [SourceMutPointer::null(); TAUT_NUM as usize];
        let mut partition_ct_atoms = [SourceMutPointer::null(); TAUT_NUM as usize];
        let mut partition_iso_rank = [SourceMutPointer::null(); TAUT_NUM as usize];
        let mut partition_iso_atoms = [SourceMutPointer::null(); TAUT_NUM as usize];

        for k in 0..TAUT_NUM as usize {
            let neighbor_atoms = heap.allocate(vec![1_u16, k as u16]).unwrap();
            neighbor_storage[k] = neighbor_atoms;
            bcn.ftcn[k].NeighList = heap.allocate(vec![neighbor_atoms]).unwrap();
            bcn.ftcn[k].LinearCt = heap.allocate(vec![10_u16 + k as u16]).unwrap();
            partition_ct_rank[k] = heap.allocate(vec![20_u16 + k as u16]).unwrap();
            partition_ct_atoms[k] = heap.allocate(vec![30_u16 + k as u16]).unwrap();
            bcn.ftcn[k].PartitionCt = Partition {
                Rank: partition_ct_rank[k],
                AtNumber: partition_ct_atoms[k],
            };
            bcn.ftcn[k].nSymmRankCt = heap.allocate(vec![40_u16 + k as u16]).unwrap();
            bcn.ftcn[k].nNumHOrig = heap.allocate(vec![50_i16 + k as i16]).unwrap();
            bcn.ftcn[k].nNumH = heap.allocate(vec![60_i16 + k as i16]).unwrap();
            bcn.ftcn[k].nNumHOrigFixH = heap.allocate(vec![70_i16 + k as i16]).unwrap();
            bcn.ftcn[k].nNumHFixH = heap.allocate(vec![80_i16 + k as i16]).unwrap();
            partition_iso_rank[k] = heap.allocate(vec![90_u16 + k as u16]).unwrap();
            partition_iso_atoms[k] = heap.allocate(vec![100_u16 + k as u16]).unwrap();
            bcn.ftcn[k].PartitionCtIso = Partition {
                Rank: partition_iso_rank[k],
                AtNumber: partition_iso_atoms[k],
            };
            bcn.ftcn[k].nSymmRankCtIso = heap.allocate(vec![110_u16 + k as u16]).unwrap();
            bcn.ftcn[k].iso_sort_keys = heap.allocate(vec![k as i64]).unwrap();
            bcn.ftcn[k].iso_sort_keysOrig = heap.allocate(vec![k as i64 + 1]).unwrap();
            bcn.ftcn[k].iso_exchg_atnos = heap.allocate(vec![1_i8 + k as i8]).unwrap();
            bcn.ftcn[k].iso_exchg_atnosOrig = heap.allocate(vec![3_i8 + k as i8]).unwrap();
            bcn.ftcn[k].num_atoms = 100 + k as i32;
        }

        assert_eq!(DeAllocBCN(&mut heap, Some(&mut bcn)), Ok(()));

        assert_freed(&heap, rank0);
        assert_freed(&heap, rank1);
        assert_freed(&heap, rank2);
        assert_freed(&heap, rank_stack);
        assert_eq!(bcn.pRankStack, rank_stack);
        assert_eq!(bcn.nMaxLenRankStack, 3);
        assert_eq!(bcn.num_max, 41);
        assert_eq!(bcn.num_at_tg, 42);
        assert_eq!(bcn.num_atoms, 43);
        assert_eq!(bcn.ulTimeOutTime, timeout);
        assert_eq!(heap.slice(timeout.as_const()).unwrap().len(), 1);

        for k in 0..TAUT_NUM as usize {
            assert_freed(&heap, neighbor_storage[k]);
            assert_freed(&heap, bcn.ftcn[k].NeighList);
            assert_freed(&heap, bcn.ftcn[k].LinearCt);
            assert_freed(&heap, partition_ct_rank[k]);
            assert_freed(&heap, partition_ct_atoms[k]);
            assert!(bcn.ftcn[k].PartitionCt.Rank.is_null());
            assert!(bcn.ftcn[k].PartitionCt.AtNumber.is_null());
            assert_freed(&heap, bcn.ftcn[k].nSymmRankCt);
            assert_freed(&heap, bcn.ftcn[k].nNumHOrig);
            assert_freed(&heap, bcn.ftcn[k].nNumH);
            assert_freed(&heap, bcn.ftcn[k].nNumHOrigFixH);
            assert_freed(&heap, bcn.ftcn[k].nNumHFixH);
            assert_freed(&heap, partition_iso_rank[k]);
            assert_freed(&heap, partition_iso_atoms[k]);
            assert!(bcn.ftcn[k].PartitionCtIso.Rank.is_null());
            assert!(bcn.ftcn[k].PartitionCtIso.AtNumber.is_null());
            assert_freed(&heap, bcn.ftcn[k].nSymmRankCtIso);
            assert_freed(&heap, bcn.ftcn[k].iso_sort_keys);
            assert_freed(&heap, bcn.ftcn[k].iso_sort_keysOrig);
            assert_freed(&heap, bcn.ftcn[k].iso_exchg_atnos);
            assert_freed(&heap, bcn.ftcn[k].iso_exchg_atnosOrig);
            assert_eq!(bcn.ftcn[k].num_atoms, 100 + k as i32);
        }
    }
}
