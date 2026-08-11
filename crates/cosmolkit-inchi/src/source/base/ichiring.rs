use crate::source::base::util::{inchi_calloc, inchi_free};
use crate::source_types::{
    AT_RANK, MAX_ATOMS, QUEUE, S_CHAR, SourceHeap, SourceHeapError, SourceMutPointer,
    StableSourceConstSlice, StableSourceSlice, inp_ATOM, qInt,
};

struct RingSearchWorkspace {
    atom: StableSourceConstSlice<inp_ATOM>,
    queue: StableSourceSlice<QUEUE>,
    values: StableSourceSlice<qInt>,
    atom_level: StableSourceSlice<AT_RANK>,
    source: StableSourceSlice<S_CHAR>,
}

impl RingSearchWorkspace {
    fn new(
        heap: &mut SourceHeap,
        atom: SourceMutPointer<inp_ATOM>,
        q: SourceMutPointer<QUEUE>,
        n_atom_level: SourceMutPointer<AT_RANK>,
        source: SourceMutPointer<S_CHAR>,
    ) -> Result<Self, SourceHeapError> {
        let queue_value = heap
            .slice(q.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let writable_ids = [
            q.allocation_identity(),
            queue_value.Val.allocation_identity(),
            n_atom_level.allocation_identity(),
            source.allocation_identity(),
        ];
        if writable_ids.iter().any(Option::is_none)
            || writable_ids
                .iter()
                .enumerate()
                .any(|(index, id)| writable_ids[..index].contains(id))
            || writable_ids.contains(&atom.allocation_identity())
        {
            return Err(SourceHeapError::PointerAllocationMismatch);
        }

        let total = usize::try_from(queue_value.nTotLength)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let first =
            usize::try_from(queue_value.nFirst).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let length = usize::try_from(queue_value.nLength)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if total == 0 || first >= total || length > total {
            return Err(SourceHeapError::PointerOutOfBounds);
        }

        // SAFETY: the five source arrays have distinct allocation identities,
        // and GetMinRingSize neither frees nor resizes them. The source graph
        // contract makes every active neighbor and queued atom an index into
        // atom[], nAtomLevel[], and cSource[].
        let atom = unsafe { heap.stable_slice(atom.as_const())? };
        let queue = unsafe { heap.stable_slice_mut(q)? };
        let values = unsafe { heap.stable_slice_mut(queue_value.Val)? };
        let atom_level = unsafe { heap.stable_slice_mut(n_atom_level)? };
        let source = unsafe { heap.stable_slice_mut(source)? };
        if values.len() < total || atom.len() < atom_level.len() || source.len() < atom_level.len()
        {
            return Err(SourceHeapError::PointerOutOfBounds);
        }

        Ok(Self {
            atom,
            queue,
            values,
            atom_level,
            source,
        })
    }

    #[inline(always)]
    unsafe fn queue(&self) -> &QUEUE {
        // SAFETY: construction requires the QUEUE allocation to contain one item.
        unsafe { self.queue.get_unchecked(0) }
    }

    #[inline(always)]
    unsafe fn queue_mut(&mut self) -> &mut QUEUE {
        // SAFETY: construction requires the QUEUE allocation to contain one item.
        unsafe { self.queue.get_unchecked_mut(0) }
    }

    #[inline(always)]
    unsafe fn queue_get(&mut self, value: &mut qInt) -> i32 {
        let queue = unsafe { self.queue() };
        if queue.nLength <= 0 {
            return -1;
        }
        let first = queue.nFirst as usize;
        let total = queue.nTotLength;
        // SAFETY: the source QUEUE invariants validated at construction and
        // preserved below keep nFirst within the allocated circular buffer.
        *value = *unsafe { self.values.get_unchecked(first) };
        let queue = unsafe { self.queue_mut() };
        queue.nFirst = if queue.nFirst == total.wrapping_sub(1) {
            0
        } else {
            queue.nFirst.wrapping_add(1)
        };
        queue.nLength = queue.nLength.wrapping_sub(1);
        queue.nLength
    }

    #[inline(always)]
    unsafe fn queue_add(&mut self, value: qInt) -> i32 {
        let queue = unsafe { self.queue() };
        if queue.nLength >= queue.nTotLength {
            return -1;
        }
        let destination = queue.nFirst.wrapping_add(queue.nLength) % queue.nTotLength;
        let length = queue.nLength.wrapping_add(1);
        // SAFETY: the source QUEUE invariants make the modulo result a valid
        // circular-buffer index.
        *unsafe { self.values.get_unchecked_mut(destination as usize) } = value;
        unsafe { self.queue_mut() }.nLength = length;
        length
    }
}

fn get_min_ring_size_with_workspace(
    workspace: &mut RingSearchWorkspace,
    n_max_ring_size: AT_RANK,
) -> i32 {
    let mut n_min_ring_size = (MAX_ATOMS + 1) as AT_RANK;
    let mut at_no: qInt = 0;

    loop {
        // SAFETY: RingSearchWorkspace validates and exclusively owns the queue.
        let q_len = unsafe { workspace.queue() }.nLength;
        if q_len == 0 {
            break;
        }
        for _i in 0..q_len {
            // SAFETY: the workspace maintains the source QUEUE invariants.
            if unsafe { workspace.queue_get(&mut at_no) } >= 0 {
                let iat_no = usize::from(at_no);
                // SAFETY: queued atom numbers obey the source graph-array contract.
                let n_cur_level =
                    unsafe { workspace.atom_level.get_unchecked(iat_no) }.wrapping_add(1);
                if 2 * i32::from(n_cur_level) > i32::from(n_max_ring_size) + 4 {
                    if u32::from(n_min_ring_size) < MAX_ATOMS + 1 {
                        return if n_min_ring_size >= n_max_ring_size {
                            0
                        } else {
                            i32::from(n_min_ring_size)
                        };
                    }
                    return 0;
                }

                let current_source = *unsafe { workspace.source.get_unchecked(iat_no) };
                // SAFETY: iat_no is a queued source atom index. Each read is
                // kept as a scalar so queue mutation does not extend a borrow
                // across the independent atom allocation.
                let valence = i32::from(unsafe { workspace.atom.get_unchecked(iat_no) }.valence);
                for j in 0..valence {
                    let next = unsafe { workspace.atom.get_unchecked(iat_no) }.neighbor[j as usize]
                        as qInt;
                    let inext = usize::from(next);
                    // SAFETY: active inp_ATOM neighbors obey the source graph contract.
                    let next_level = *unsafe { workspace.atom_level.get_unchecked(inext) };
                    if next_level == 0 {
                        // SAFETY: the workspace maintains the source QUEUE invariants.
                        if unsafe { workspace.queue_add(next) } >= 0 {
                            *unsafe { workspace.atom_level.get_unchecked_mut(inext) } = n_cur_level;
                            *unsafe { workspace.source.get_unchecked_mut(inext) } = current_source;
                        } else {
                            return -1;
                        }
                    } else {
                        let next_source = *unsafe { workspace.source.get_unchecked(inext) };
                        if i32::from(next_level) + 1 >= i32::from(n_cur_level)
                            && next_source != current_source
                        {
                            if next_source == -1 {
                                return -1;
                            }
                            let n_ring_size = next_level.wrapping_add(n_cur_level).wrapping_sub(2);
                            if n_ring_size < n_min_ring_size {
                                n_min_ring_size = n_ring_size;
                            }
                        }
                    }
                }
            } else {
                return -1;
            }
        }
    }

    if u32::from(n_min_ring_size) < MAX_ATOMS + 1 {
        return if n_min_ring_size >= n_max_ring_size {
            0
        } else {
            i32::from(n_min_ring_size)
        };
    }
    0
}

#[allow(non_snake_case)]
pub(crate) fn QueueCreate(
    heap: &mut SourceHeap,
    nTotLength: i32,
    nSize: i32,
) -> Result<SourceMutPointer<QUEUE>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:67 QueueCreate
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap access adds overhead.
    /*
    QUEUE *QueueCreate( int nTotLength, int nSize )
    {
        QUEUE     *q = NULL;
        QINT_TYPE *Val = NULL;
        if (nTotLength < 1 || nSize != ( int )sizeof( QINT_TYPE ) ||
             !( q = (QUEUE     *) inchi_calloc( 1, sizeof( QUEUE ) ) ) ||
             !( Val = (QINT_TYPE *) inchi_calloc( nTotLength, nSize ) ))
        {
            if (q) inchi_free( q );
            return NULL;
        }
        q->Val = Val;
        /* q->nSize      = nSize; */
        q->nTotLength = nTotLength;

        return q;
    }
    */
    // END INCHI C FUNCTION: QueueCreate
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: QueueCreate
    // INCHI✔️❌: #define QUEUE_QINT 1
    // INCHI✔️❌: typedef AT_RANK qInt;
    // INCHI✔️❌: typedef qInt QINT_TYPE;
    // END INCHI ACTIVE MACRO CONFIGURATION: QueueCreate

    if nTotLength < 1 || nSize != 2 {
        return Ok(SourceMutPointer::null());
    }

    let queue = match inchi_calloc::<QUEUE>(heap, 1, 24) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return Ok(SourceMutPointer::null()),
        Err(error) => return Err(error),
    };
    let values = match inchi_calloc::<qInt>(heap, nTotLength as u64, nSize as u64) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => {
            inchi_free(heap, queue)?;
            return Ok(SourceMutPointer::null());
        }
        Err(error) => {
            inchi_free(heap, queue)?;
            return Err(error);
        }
    };

    let queue_value = heap
        .slice_mut(queue)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    queue_value.Val = values;
    queue_value.nTotLength = nTotLength;
    Ok(queue)
}

#[allow(non_snake_case)]
pub(crate) fn QueueAdd(
    heap: &mut SourceHeap,
    q: SourceMutPointer<QUEUE>,
    Val: SourceMutPointer<qInt>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:87 QueueAdd
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap access adds overhead.
    /*
    int QueueAdd( QUEUE *q, QINT_TYPE *Val )
    {
        if (q && Val && q->nLength < q->nTotLength)
        {
            q->Val[( q->nFirst + q->nLength ) % q->nTotLength] = *Val;
            q->nLength++;
            return q->nLength;
        }

        return -1;
    }
    */
    // END INCHI C FUNCTION: QueueAdd

    if q.is_null() || Val.is_null() {
        return Ok(-1);
    }
    let queue = heap
        .slice(q.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    if queue.nLength >= queue.nTotLength {
        return Ok(-1);
    }
    if queue.nTotLength == 0 {
        return Err(SourceHeapError::SourceIntegerOverflow);
    }
    let destination = queue.nFirst.wrapping_add(queue.nLength) % queue.nTotLength;
    let destination =
        usize::try_from(destination).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let value = *heap
        .slice(Val.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    *heap
        .slice_mut(queue.Val)?
        .get_mut(destination)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = value;
    let queue_mut = heap
        .slice_mut(q)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    queue_mut.nLength = queue.nLength.wrapping_add(1);
    Ok(queue_mut.nLength)
}

#[allow(non_snake_case)]
pub(crate) fn QueueGet(
    heap: &mut SourceHeap,
    q: SourceMutPointer<QUEUE>,
    Val: SourceMutPointer<qInt>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:101 QueueGet
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap access adds overhead.
    /*
    int QueueGet( QUEUE *q, QINT_TYPE *Val )
    {
        if (q && Val && q->nLength > 0)
        {
            *Val = q->Val[q->nFirst];
            /* new: do not allow to overwrite the retrieved value */
            q->nFirst = ( q->nFirst == q->nTotLength - 1 ) ? 0 : q->nFirst + 1;
            q->nLength--;
            /* -- old --
            if ( -- q->nLength ) {
                q->nFirst = (q->nFirst == q->nTotLength - 1)? 0 : q->nFirst + 1;
            }
            */
            return q->nLength;
        }

        return -1;
    }
    */
    // END INCHI C FUNCTION: QueueGet

    if q.is_null() || Val.is_null() {
        return Ok(-1);
    }
    let queue = heap
        .slice(q.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    if queue.nLength <= 0 {
        return Ok(-1);
    }
    let first = usize::try_from(queue.nFirst).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let value = *heap
        .slice(queue.Val.as_const())?
        .get(first)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    *heap
        .slice_mut(Val)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = value;
    let queue_mut = heap
        .slice_mut(q)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    queue_mut.nFirst = if queue.nFirst == queue.nTotLength.wrapping_sub(1) {
        0
    } else {
        queue.nFirst.wrapping_add(1)
    };
    queue_mut.nLength = queue.nLength.wrapping_sub(1);
    Ok(queue_mut.nLength)
}

#[allow(non_snake_case)]
pub(crate) fn QueueGetAny(
    heap: &mut SourceHeap,
    q: SourceMutPointer<QUEUE>,
    Val: SourceMutPointer<qInt>,
    ord: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:122 QueueGetAny
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap access adds overhead.
    /*
    int QueueGetAny( QUEUE *q, QINT_TYPE *Val, int ord )
    {
        if (0 <= ord && ord < q->nTotLength)
        {
            *Val = q->Val[ord];
            return  1; /*  success */
        }
        else
        {
            return -1; /*  error */
        }
    }
    */
    // END INCHI C FUNCTION: QueueGetAny

    let queue = heap
        .slice(q.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    if 0 <= ord && ord < queue.nTotLength {
        let value = *heap
            .slice(queue.Val.as_const())?
            .get(usize::try_from(ord).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        *heap
            .slice_mut(Val)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)? = value;
        Ok(1)
    } else {
        Ok(-1)
    }
}

#[allow(non_snake_case)]
pub(crate) fn QueueReinit(
    heap: &mut SourceHeap,
    q: SourceMutPointer<QUEUE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:216 QueueReinit
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap access adds overhead.
    /*
    int QueueReinit( QUEUE *q )
    {
        if (q)
        {
            q->nFirst = 0;
            q->nLength = 0;
            /* memset( q->Val, 0, q->nTotLength*sizeof(q->Val[0])); */ /*  for debug only */
            return q->nTotLength;
        }

        return -1;
    }
    */
    // END INCHI C FUNCTION: QueueReinit

    if q.is_null() {
        return Ok(-1);
    }
    let queue = heap
        .slice_mut(q)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    queue.nFirst = 0;
    queue.nLength = 0;
    Ok(queue.nTotLength)
}

#[allow(non_snake_case)]
pub(crate) fn QueueLength(
    heap: &SourceHeap,
    q: SourceMutPointer<QUEUE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:231 QueueLength
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap access adds overhead.
    /*
    int QueueLength( QUEUE *q )
    {
        if (q)
        {
            return q->nLength;
        }
        else
        {
            return 0;
        }
    }
    */
    // END INCHI C FUNCTION: QueueLength

    if q.is_null() {
        return Ok(0);
    }
    Ok(heap
        .slice(q.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .nLength)
}

#[allow(non_snake_case)]
pub(crate) fn QueueWrittenLength(
    heap: &SourceHeap,
    q: SourceMutPointer<QUEUE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:245 QueueWrittenLength
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap access adds overhead.
    /*
    int QueueWrittenLength( QUEUE *q )
    {
        if (q)
        {
            int len = q->nFirst + q->nLength;
            return ( len > q->nTotLength ) ? q->nTotLength : len;
        }
        else
        {
            return 0;
        }
    }
    */
    // END INCHI C FUNCTION: QueueWrittenLength

    if q.is_null() {
        return Ok(0);
    }
    let queue = heap
        .slice(q.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let len = queue.nFirst.wrapping_add(queue.nLength);
    Ok(if len > queue.nTotLength {
        queue.nTotLength
    } else {
        len
    })
}

#[allow(non_snake_case)]
pub(crate) fn GetMinRingSize(
    heap: &mut SourceHeap,
    atom: SourceMutPointer<inp_ATOM>,
    q: SourceMutPointer<QUEUE>,
    nAtomLevel: SourceMutPointer<AT_RANK>,
    cSource: SourceMutPointer<S_CHAR>,
    nMaxRingSize: AT_RANK,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:262 GetMinRingSize
    // INCHI✔️✔️: int GetMinRingSize( inp_ATOM* atom,
    // INCHI✔️✔️:                     QUEUE *q,
    // INCHI✔️✔️:                     AT_RANK *nAtomLevel,
    // INCHI✔️✔️:                     S_CHAR *cSource,
    // INCHI✔️✔️:                     AT_RANK nMaxRingSize )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int qLen, i, j;
    // INCHI✔️✔️:     AT_RANK nCurLevel, nRingSize, nMinRingSize = MAX_ATOMS + 1;
    // INCHI✔️✔️:     qInt at_no, next;
    // INCHI✔️✔️:     int  iat_no, inext;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     while ((qLen = QueueLength( q ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /*  traverse the next level (next outer ring) */
    // INCHI✔️✔️:         for (i = 0; i < qLen; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (0 <= QueueGet( q, &at_no ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 iat_no = (int) at_no;
    // INCHI✔️✔️:                 nCurLevel = nAtomLevel[iat_no] + 1;
    // INCHI✔️✔️:                 if (2 * nCurLevel > nMaxRingSize + 4)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /*  2*nCurLevel = nRingSize + 3 + k, k = 0 or 1  */
    // INCHI✔️✔️:                     if (nMinRingSize < MAX_ATOMS + 1)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         return ( nMinRingSize >= nMaxRingSize ) ? 0 : nMinRingSize;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     return 0; /*  min. ring size > nMaxRingSize */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 for (j = 0; j < atom[iat_no].valence; j++)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     next = (qInt) atom[iat_no].neighbor[j];
    // INCHI✔️✔️:                     inext = (int) next;
    // INCHI✔️✔️:                     if (!nAtomLevel[inext])
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         /*  the at_no neighbor has not been traversed yet. Add it to the queue */
    // INCHI✔️✔️:                         if (0 <= QueueAdd( q, &next ))
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             nAtomLevel[inext] = nCurLevel;
    // INCHI✔️✔️:                             cSource[inext] = cSource[iat_no]; /*  keep the path number */
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         else
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             return -1; /*  error */
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     else
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         if (nAtomLevel[inext] + 1 >= nCurLevel &&
    // INCHI✔️✔️:                                 cSource[inext] != cSource[iat_no]
    // INCHI✔️✔️:                                 /*  && cSource[(int)next] != -1 */
    // INCHI✔️✔️:                               )
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             /*  found a ring closure */
    // INCHI✔️✔️:                             /*  debug */
    // INCHI✔️✔️:                             if (cSource[inext] == -1)
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 return -1;  /*  error */
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                             if (( nRingSize = nAtomLevel[inext] + nCurLevel - 2 ) < nMinRingSize)
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 nMinRingSize = nRingSize;
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                             /* return (nRingSize >= nMaxRingSize)? 0 : nRingSize; */
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 return -1; /*  error */
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (nMinRingSize < MAX_ATOMS + 1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return ( nMinRingSize >= nMaxRingSize ) ? 0 : nMinRingSize;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: GetMinRingSize

    if q.is_null() {
        return Ok(0);
    }
    let initial_length = heap
        .slice(q.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .nLength;
    if initial_length == 0 {
        return Ok(0);
    }
    let mut workspace = RingSearchWorkspace::new(heap, atom, q, nAtomLevel, cSource)?;
    Ok(get_min_ring_size_with_workspace(
        &mut workspace,
        nMaxRingSize,
    ))
}

#[allow(non_snake_case)]
pub(crate) fn is_bond_in_Nmax_memb_ring(
    heap: &mut SourceHeap,
    atom: SourceMutPointer<inp_ATOM>,
    at_no: i32,
    neigh_ord: i32,
    q: SourceMutPointer<QUEUE>,
    nAtomLevel: SourceMutPointer<AT_RANK>,
    cSource: SourceMutPointer<S_CHAR>,
    nMaxRingSize: AT_RANK,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:362 is_bond_in_Nmax_memb_ring
    // INCHI✔️✔️: int is_bond_in_Nmax_memb_ring( inp_ATOM* atom,
    // INCHI✔️✔️:                                int at_no,
    // INCHI✔️✔️:                                int neigh_ord,
    // INCHI✔️✔️:                                QUEUE *q,
    // INCHI✔️✔️:                                AT_RANK *nAtomLevel,
    // INCHI✔️✔️:                                S_CHAR *cSource,
    // INCHI✔️✔️:                                AT_RANK nMaxRingSize )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int  nMinRingSize = -1, i;
    // INCHI✔️✔️:     qInt n;
    // INCHI✔️✔️:     int  nTotLen;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (nMaxRingSize < 3)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     QueueReinit( q );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*  mark the starting atom */
    // INCHI✔️✔️:     nAtomLevel[at_no] = 1;
    // INCHI✔️✔️:     cSource[at_no] = -1;
    // INCHI✔️✔️:     /*  add neighbors */
    // INCHI✔️✔️:     for (i = 0; i < atom[at_no].valence; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         n = (qInt) atom[at_no].neighbor[i];
    // INCHI✔️✔️:         nAtomLevel[(int) n] = 2;
    // INCHI✔️✔️:         cSource[(int) n] = 1 + ( i == neigh_ord );
    // INCHI✔️✔️:         QueueAdd( q, &n );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     nMinRingSize = GetMinRingSize( atom, q, nAtomLevel, cSource, nMaxRingSize );
    // INCHI✔️✔️:     /*  cleanup */
    // INCHI✔️✔️:     nTotLen = QueueWrittenLength( q );
    // INCHI✔️✔️:     for (i = 0; i < nTotLen; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (0 < QueueGetAny( q, &n, i ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             nAtomLevel[(int) n] = 0;
    // INCHI✔️✔️:             cSource[(int) n] = 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     nAtomLevel[at_no] = 0;
    // INCHI✔️✔️:     cSource[at_no] = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*
    // INCHI✔️✔️:     if ( nAtomLevel )
    // INCHI✔️✔️:         inchi_free ( nAtomLevel );
    // INCHI✔️✔️:     if ( cSource )
    // INCHI✔️✔️:         inchi_free ( cSource );
    // INCHI✔️✔️:     QueueDelete( q );
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return nMinRingSize;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: is_bond_in_Nmax_memb_ring

    if nMaxRingSize < 3 {
        return Ok(0);
    }

    let queue = heap
        .slice_mut(q)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    queue.nFirst = 0;
    queue.nLength = 0;

    let mut workspace = RingSearchWorkspace::new(heap, atom, q, nAtomLevel, cSource)?;
    let start = usize::try_from(at_no).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if start >= workspace.atom.len()
        || start >= workspace.atom_level.len()
        || start >= workspace.source.len()
    {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    // SAFETY: start was checked against all three source arrays above.
    *unsafe { workspace.atom_level.get_unchecked_mut(start) } = 1;
    *unsafe { workspace.source.get_unchecked_mut(start) } = -1;

    // SAFETY: start was checked against atom[] above.
    let valence = i32::from(unsafe { workspace.atom.get_unchecked(start) }.valence);
    for i in 0..valence {
        let n = unsafe { workspace.atom.get_unchecked(start) }.neighbor[i as usize] as qInt;
        let neighbor = usize::from(n);
        // SAFETY: active inp_ATOM neighbors obey the source graph contract.
        *unsafe { workspace.atom_level.get_unchecked_mut(neighbor) } = 2;
        *unsafe { workspace.source.get_unchecked_mut(neighbor) } = 1 + i8::from(i == neigh_ord);
        // SAFETY: the workspace maintains the source QUEUE invariants. The C
        // caller intentionally ignores QueueAdd's return value here.
        let _ = unsafe { workspace.queue_add(n) };
    }

    let n_min_ring_size = get_min_ring_size_with_workspace(&mut workspace, nMaxRingSize);
    // SAFETY: the workspace owns a validated QUEUE.
    let queue = unsafe { workspace.queue() };
    let written = queue.nFirst.wrapping_add(queue.nLength);
    let n_tot_len = if written > queue.nTotLength {
        queue.nTotLength
    } else {
        written
    };
    for i in 0..n_tot_len {
        // SAFETY: QueueWrittenLength bounds i by nTotLength, and construction
        // validates the backing circular-buffer capacity.
        let n = usize::from(*unsafe { workspace.values.get_unchecked(i as usize) });
        *unsafe { workspace.atom_level.get_unchecked_mut(n) } = 0;
        *unsafe { workspace.source.get_unchecked_mut(n) } = 0;
    }
    *unsafe { workspace.atom_level.get_unchecked_mut(start) } = 0;
    *unsafe { workspace.source.get_unchecked_mut(start) } = 0;

    Ok(n_min_ring_size)
}

pub(crate) fn is_atom_in_3memb_ring(
    heap: &SourceHeap,
    atom: SourceMutPointer<inp_ATOM>,
    at_no: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:420 is_atom_in_3memb_ring
    // INCHI✔️✔️: int is_atom_in_3memb_ring( inp_ATOM* atom, int at_no )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     AT_NUMB   neigh_neigh;
    // INCHI✔️✔️:     int       i, j, k, val, val_neigh, neigh;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (atom[at_no].nNumAtInRingSystem < 3)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i = 0, val = atom[at_no].valence; i < val; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         neigh = (int) atom[at_no].neighbor[i];
    // INCHI✔️✔️:         if (atom[at_no].nRingSystem != atom[neigh].nRingSystem)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             continue;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         for (j = 0, val_neigh = atom[neigh].valence; j < val_neigh; j++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             neigh_neigh = atom[neigh].neighbor[j];
    // INCHI✔️✔️:             if ((int) neigh_neigh == at_no)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 continue;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             for (k = 0; k < val; k++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (atom[at_no].neighbor[k] == neigh_neigh)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     return 1;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: is_atom_in_3memb_ring

    let start_index = usize::try_from(at_no).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    // SAFETY: this function is read-only and does not free or resize atom[].
    let atoms = unsafe { heap.stable_slice(atom.as_const())? };
    let start = atoms
        .get(start_index)
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if start.nNumAtInRingSystem < 3 {
        return Ok(0);
    }
    let valence = i32::from(start.valence);
    for i in 0..valence {
        let neighbor_index = usize::from(start.neighbor[i as usize]);
        let neighbor = atoms
            .get(neighbor_index)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if start.nRingSystem != neighbor.nRingSystem {
            continue;
        }
        for j in 0..i32::from(neighbor.valence) {
            let neighbor_neighbor = neighbor.neighbor[j as usize];
            if i32::from(neighbor_neighbor) == at_no {
                continue;
            }
            for k in 0..valence {
                if start.neighbor[k as usize] == neighbor_neighbor {
                    return Ok(1);
                }
            }
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn QueueDelete(
    heap: &mut SourceHeap,
    q: SourceMutPointer<QUEUE>,
) -> Result<SourceMutPointer<QUEUE>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:203 QueueDelete
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap access adds overhead.
    /*
    QUEUE *QueueDelete( QUEUE *q )
    {
        if (q)
        {
            if (q->Val) inchi_free( q->Val );
            inchi_free( q );
        }

        return NULL;
    }
    */
    // END INCHI C FUNCTION: QueueDelete

    if q.is_null() {
        return Ok(SourceMutPointer::null());
    }

    let value = heap
        .slice(q.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .Val;
    let value_result = inchi_free(heap, value);
    let queue_result = inchi_free(heap, q);
    value_result?;
    queue_result?;
    Ok(SourceMutPointer::null())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ring_atom(neighbors: &[u16]) -> inp_ATOM {
        let mut atom = inp_ATOM::default();
        atom.valence = neighbors.len() as i8;
        atom.neighbor[..neighbors.len()].copy_from_slice(neighbors);
        atom
    }

    fn ring_queue(
        heap: &mut SourceHeap,
        values: Vec<u16>,
        total: i32,
        first: i32,
        length: i32,
    ) -> (SourceMutPointer<u16>, SourceMutPointer<QUEUE>) {
        let values = heap.allocate_model_storage(values).unwrap();
        let queue = heap
            .allocate_model_storage(vec![QUEUE {
                Val: values,
                nTotLength: total,
                nFirst: first,
                nLength: length,
            }])
            .unwrap();
        (values, queue)
    }

    #[test]
    fn source_port__ichiring__queuecreate__line_67() {
        for (length, size) in [(0, 2), (-1, 2), (1, 0), (1, 1), (1, 3)] {
            let mut heap = SourceHeap::default();
            assert_eq!(
                QueueCreate(&mut heap, length, size),
                Ok(SourceMutPointer::null())
            );
            assert_eq!(heap.live_allocation_count(), 0);
            assert_eq!(heap.source_allocation_calls(), 0);
        }

        let mut heap = SourceHeap::default();
        heap.fail_after_allocations(0);
        assert_eq!(QueueCreate(&mut heap, 3, 2), Ok(SourceMutPointer::null()));
        assert_eq!(heap.live_allocation_count(), 0);
        assert_eq!(heap.source_allocation_calls(), 1);

        let mut heap = SourceHeap::default();
        heap.fail_after_allocations(1);
        assert_eq!(QueueCreate(&mut heap, 3, 2), Ok(SourceMutPointer::null()));
        assert_eq!(heap.live_allocation_count(), 0);
        assert_eq!(heap.source_allocation_calls(), 2);

        let mut heap = SourceHeap::default();
        let queue = QueueCreate(&mut heap, 3, 2).unwrap();
        assert!(!queue.is_null());
        assert_eq!(heap.live_allocation_count(), 2);
        let queue_value = heap.slice(queue.as_const()).unwrap()[0].clone();
        assert_eq!(queue_value.nTotLength, 3);
        assert_eq!(queue_value.nFirst, 0);
        assert_eq!(queue_value.nLength, 0);
        assert_eq!(heap.slice(queue_value.Val.as_const()).unwrap(), &[0, 0, 0]);
        assert_eq!(QueueDelete(&mut heap, queue), Ok(SourceMutPointer::null()));
        assert_eq!(heap.live_allocation_count(), 0);
    }

    #[test]
    fn source_port__ichiring__queueadd__line_87() {
        let mut heap = SourceHeap::default();
        let input = heap.allocate_model_storage(vec![7_u16]).unwrap();
        assert_eq!(QueueAdd(&mut heap, SourceMutPointer::null(), input), Ok(-1));

        let values = heap.allocate_model_storage(vec![10_u16, 20, 30]).unwrap();
        let queue = heap
            .allocate_model_storage(vec![QUEUE {
                Val: values,
                nTotLength: 3,
                nFirst: 2,
                nLength: 0,
            }])
            .unwrap();
        assert_eq!(QueueAdd(&mut heap, queue, input), Ok(1));
        assert_eq!(heap.slice(values.as_const()).unwrap(), &[10, 20, 7]);

        heap.slice_mut(input).unwrap()[0] = 11;
        assert_eq!(QueueAdd(&mut heap, queue, input), Ok(2));
        assert_eq!(heap.slice(values.as_const()).unwrap(), &[11, 20, 7]);

        assert_eq!(QueueAdd(&mut heap, queue, values.offset(1).unwrap()), Ok(3));
        assert_eq!(heap.slice(values.as_const()).unwrap(), &[11, 20, 7]);
        assert_eq!(heap.slice(queue.as_const()).unwrap()[0].nLength, 3);

        heap.slice_mut(input).unwrap()[0] = 99;
        assert_eq!(QueueAdd(&mut heap, queue, input), Ok(-1));
        assert_eq!(QueueAdd(&mut heap, queue, SourceMutPointer::null()), Ok(-1));
        assert_eq!(heap.slice(values.as_const()).unwrap(), &[11, 20, 7]);
        assert_eq!(heap.slice(queue.as_const()).unwrap()[0].nLength, 3);

        heap.slice_mut(queue).unwrap()[0].nLength = -1;
        heap.slice_mut(queue).unwrap()[0].nTotLength = 0;
        assert_eq!(
            QueueAdd(&mut heap, queue, input),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
    }

    #[test]
    fn source_port__ichiring__queueget__line_101() {
        let mut heap = SourceHeap::default();
        let output = heap.allocate_model_storage(vec![99_u16]).unwrap();
        assert_eq!(
            QueueGet(&mut heap, SourceMutPointer::null(), output),
            Ok(-1)
        );

        let values = heap.allocate_model_storage(vec![11_u16, 22, 33]).unwrap();
        let queue = heap
            .allocate_model_storage(vec![QUEUE {
                Val: values,
                nTotLength: 3,
                nFirst: 2,
                nLength: 2,
            }])
            .unwrap();
        assert_eq!(QueueGet(&mut heap, queue, output), Ok(1));
        assert_eq!(heap.slice(output.as_const()).unwrap()[0], 33);
        assert_eq!(heap.slice(queue.as_const()).unwrap()[0].nFirst, 0);
        assert_eq!(heap.slice(queue.as_const()).unwrap()[0].nLength, 1);

        assert_eq!(QueueGet(&mut heap, queue, output), Ok(0));
        assert_eq!(heap.slice(output.as_const()).unwrap()[0], 11);
        assert_eq!(heap.slice(queue.as_const()).unwrap()[0].nFirst, 1);
        assert_eq!(heap.slice(queue.as_const()).unwrap()[0].nLength, 0);

        heap.slice_mut(output).unwrap()[0] = 77;
        assert_eq!(QueueGet(&mut heap, queue, output), Ok(-1));
        assert_eq!(heap.slice(output.as_const()).unwrap()[0], 77);
        assert_eq!(heap.slice(queue.as_const()).unwrap()[0].nFirst, 1);
        assert_eq!(QueueGet(&mut heap, queue, SourceMutPointer::null()), Ok(-1));

        heap.slice_mut(queue).unwrap()[0].nLength = 1;
        heap.slice_mut(queue).unwrap()[0].nFirst = -1;
        assert_eq!(
            QueueGet(&mut heap, queue, output),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(output.as_const()).unwrap()[0], 77);
    }

    #[test]
    fn source_port__ichiring__queuegetany__line_122() {
        let mut heap = SourceHeap::default();
        let values = heap
            .allocate_model_storage(vec![7_u16, 0, u16::MAX])
            .unwrap();
        let queue = heap
            .allocate_model_storage(vec![QUEUE {
                Val: values,
                nTotLength: 3,
                nFirst: 2,
                nLength: 0,
            }])
            .unwrap();
        let output = heap.allocate_model_storage(vec![99_u16]).unwrap();

        for (ord, expected) in [(0, 7_u16), (1, 0), (2, u16::MAX)] {
            heap.slice_mut(output).unwrap()[0] = 99;
            assert_eq!(QueueGetAny(&mut heap, queue, output, ord), Ok(1));
            assert_eq!(heap.slice(output.as_const()).unwrap()[0], expected);
        }
        for ord in [i32::MIN, -1, 3, i32::MAX] {
            heap.slice_mut(output).unwrap()[0] = 99;
            assert_eq!(
                QueueGetAny(&mut heap, queue, SourceMutPointer::null(), ord),
                Ok(-1)
            );
            assert_eq!(heap.slice(output.as_const()).unwrap()[0], 99);
        }
        assert_eq!(
            QueueGetAny(&mut heap, queue, SourceMutPointer::null(), 0),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(
            QueueGetAny(&mut heap, SourceMutPointer::null(), output, 0),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichiring__queuereinit__line_216() {
        let mut heap = SourceHeap::default();
        assert_eq!(QueueReinit(&mut heap, SourceMutPointer::null()), Ok(-1));

        let values = heap.allocate_model_storage(vec![7_u16, 11, 13]).unwrap();
        let queue = heap
            .allocate_model_storage(vec![QUEUE {
                Val: values,
                nTotLength: 3,
                nFirst: 2,
                nLength: 3,
            }])
            .unwrap();
        assert_eq!(QueueReinit(&mut heap, queue), Ok(3));
        assert_eq!(heap.slice(values.as_const()).unwrap(), &[7, 11, 13]);
        assert_eq!(heap.slice(queue.as_const()).unwrap()[0].nFirst, 0);
        assert_eq!(heap.slice(queue.as_const()).unwrap()[0].nLength, 0);

        heap.slice_mut(queue).unwrap()[0].nTotLength = i32::MIN;
        heap.slice_mut(queue).unwrap()[0].nFirst = i32::MAX;
        heap.slice_mut(queue).unwrap()[0].nLength = -7;
        assert_eq!(QueueReinit(&mut heap, queue), Ok(i32::MIN));
        assert_eq!(heap.slice(queue.as_const()).unwrap()[0].nFirst, 0);
        assert_eq!(heap.slice(queue.as_const()).unwrap()[0].nLength, 0);
    }

    #[test]
    fn source_port__ichiring__queuelength__line_231() {
        let mut heap = SourceHeap::default();
        assert_eq!(QueueLength(&heap, SourceMutPointer::null()), Ok(0));
        let queue = heap.allocate_model_storage(vec![QUEUE::default()]).unwrap();
        for length in [i32::MIN, -1, 0, 1, i32::MAX] {
            heap.slice_mut(queue).unwrap()[0].nLength = length;
            assert_eq!(QueueLength(&heap, queue), Ok(length));
        }
    }

    #[test]
    fn source_port__ichiring__queuewrittenlength__line_245() {
        let mut heap = SourceHeap::default();
        assert_eq!(QueueWrittenLength(&heap, SourceMutPointer::null()), Ok(0));
        let queue = heap
            .allocate_model_storage(vec![QUEUE {
                nTotLength: 7,
                nFirst: 2,
                nLength: 3,
                ..QUEUE::default()
            }])
            .unwrap();

        for (first, length, total, expected) in [
            (0, 0, 7, 0),
            (2, 3, 7, 5),
            (4, 3, 7, 7),
            (6, 3, 7, 7),
            (-3, 1, 7, -2),
            (i32::MAX, 0, i32::MAX, i32::MAX),
            (i32::MAX, 1, i32::MAX, i32::MIN),
            (i32::MIN, -1, i32::MAX, i32::MAX),
        ] {
            let value = &mut heap.slice_mut(queue).unwrap()[0];
            value.nFirst = first;
            value.nLength = length;
            value.nTotLength = total;
            assert_eq!(QueueWrittenLength(&heap, queue), Ok(expected));
        }
    }

    #[test]
    fn source_port__ichiring__getminringsize__line_262() {
        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(vec![ring_atom(&[])]).unwrap();
        let levels = heap.allocate_model_storage(vec![7_u16]).unwrap();
        let sources = heap.allocate_model_storage(vec![3_i8]).unwrap();
        let (_values, queue) = ring_queue(&mut heap, vec![0], 1, 0, 0);
        let allocations = heap.live_allocation_count();
        assert_eq!(
            GetMinRingSize(&mut heap, atoms, queue, levels, sources, 10),
            Ok(0)
        );
        assert_eq!(heap.slice(levels.as_const()).unwrap(), &[7]);
        assert_eq!(heap.slice(sources.as_const()).unwrap(), &[3]);
        assert_eq!(heap.live_allocation_count(), allocations);

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![ring_atom(&[1]), ring_atom(&[2]), ring_atom(&[])])
            .unwrap();
        let levels = heap.allocate_model_storage(vec![1_u16, 0, 0]).unwrap();
        let sources = heap.allocate_model_storage(vec![4_i8, 0, 0]).unwrap();
        let (values, queue) = ring_queue(&mut heap, vec![0, 91, 92], 3, 0, 1);
        assert_eq!(
            GetMinRingSize(&mut heap, atoms, queue, levels, sources, 10),
            Ok(0)
        );
        assert_eq!(heap.slice(levels.as_const()).unwrap(), &[1, 2, 3]);
        assert_eq!(heap.slice(sources.as_const()).unwrap(), &[4, 4, 4]);
        assert_eq!(heap.slice(values.as_const()).unwrap(), &[0, 1, 2]);
        assert_eq!(heap.slice(queue.as_const()).unwrap()[0].nFirst, 0);
        assert_eq!(heap.slice(queue.as_const()).unwrap()[0].nLength, 0);

        for (maximum, expected) in [(10_u16, 2), (2_u16, 0)] {
            let mut heap = SourceHeap::default();
            let atoms = heap
                .allocate_model_storage(vec![ring_atom(&[2]), ring_atom(&[2]), ring_atom(&[0, 1])])
                .unwrap();
            let levels = heap.allocate_model_storage(vec![1_u16, 1, 0]).unwrap();
            let sources = heap.allocate_model_storage(vec![1_i8, 2, 0]).unwrap();
            let (_values, queue) = ring_queue(&mut heap, vec![0, 1, 99], 3, 0, 2);
            assert_eq!(
                GetMinRingSize(&mut heap, atoms, queue, levels, sources, maximum),
                Ok(expected)
            );
            assert_eq!(heap.slice(levels.as_const()).unwrap(), &[1, 1, 2]);
            assert_eq!(heap.slice(sources.as_const()).unwrap(), &[1, 2, 1]);
        }

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![ring_atom(&[2]), ring_atom(&[]), ring_atom(&[])])
            .unwrap();
        let levels = heap.allocate_model_storage(vec![2_u16, 7, 2]).unwrap();
        let sources = heap.allocate_model_storage(vec![1_i8, 1, 2]).unwrap();
        let (_values, queue) = ring_queue(&mut heap, vec![0, 1], 2, 0, 2);
        assert_eq!(
            GetMinRingSize(&mut heap, atoms, queue, levels, sources, 10),
            Ok(3)
        );
        assert_eq!(heap.slice(queue.as_const()).unwrap()[0].nLength, 0);

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![ring_atom(&[1, 2]), ring_atom(&[]), ring_atom(&[])])
            .unwrap();
        let levels = heap.allocate_model_storage(vec![1_u16, 0, 0]).unwrap();
        let sources = heap.allocate_model_storage(vec![5_i8, 0, 0]).unwrap();
        let (values, queue) = ring_queue(&mut heap, vec![0], 1, 0, 1);
        assert_eq!(
            GetMinRingSize(&mut heap, atoms, queue, levels, sources, 10),
            Ok(-1)
        );
        assert_eq!(heap.slice(levels.as_const()).unwrap(), &[1, 2, 0]);
        assert_eq!(heap.slice(sources.as_const()).unwrap(), &[5, 5, 0]);
        assert_eq!(heap.slice(values.as_const()).unwrap(), &[1]);
        assert_eq!(heap.slice(queue.as_const()).unwrap()[0].nLength, 1);

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![ring_atom(&[1]), ring_atom(&[])])
            .unwrap();
        let levels = heap.allocate_model_storage(vec![1_u16, 1]).unwrap();
        let sources = heap.allocate_model_storage(vec![1_i8, -1]).unwrap();
        let (_values, queue) = ring_queue(&mut heap, vec![0], 1, 0, 1);
        assert_eq!(
            GetMinRingSize(&mut heap, atoms, queue, levels, sources, 10),
            Ok(-1)
        );

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(vec![ring_atom(&[])]).unwrap();
        let levels = heap.allocate_model_storage(vec![u16::MAX]).unwrap();
        let sources = heap.allocate_model_storage(vec![i8::MIN]).unwrap();
        let (_values, queue) = ring_queue(&mut heap, vec![0], 1, 0, 1);
        assert_eq!(
            GetMinRingSize(&mut heap, atoms, queue, levels, sources, u16::MAX),
            Ok(0)
        );

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(vec![ring_atom(&[])]).unwrap();
        let levels = heap.allocate_model_storage(vec![1_u16]).unwrap();
        let sources = heap.allocate_model_storage(vec![1_i8]).unwrap();
        let (_values, queue) = ring_queue(&mut heap, vec![0], 1, 1, 1);
        let allocations = heap.live_allocation_count();
        assert_eq!(
            GetMinRingSize(&mut heap, atoms, queue, levels, sources, 10),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.live_allocation_count(), allocations);
    }

    #[test]
    fn get_min_ring_size_accepts_official_atom_tail_layout() {
        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![
                ring_atom(&[2]),
                ring_atom(&[2]),
                ring_atom(&[0, 1]),
                ring_atom(&[]),
            ])
            .unwrap();
        let levels = heap.allocate_model_storage(vec![1_u16, 1, 0]).unwrap();
        let sources = heap.allocate_model_storage(vec![1_i8, 2, 0]).unwrap();
        let (_values, queue) = ring_queue(&mut heap, vec![0, 1, 99, 99], 4, 0, 2);

        assert_eq!(
            GetMinRingSize(&mut heap, atoms, queue, levels, sources, 10),
            Ok(2)
        );
        assert_eq!(heap.slice(levels.as_const()).unwrap(), &[1, 1, 2]);
        assert_eq!(heap.slice(sources.as_const()).unwrap(), &[1, 2, 1]);
    }

    #[test]
    fn source_port__ichiring__is_bond_in_nmax_memb_ring__line_362() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            is_bond_in_Nmax_memb_ring(
                &mut heap,
                SourceMutPointer::null(),
                i32::MIN,
                i32::MAX,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                2,
            ),
            Ok(0)
        );
        assert_eq!(heap.live_allocation_count(), 0);

        for (maximum, expected) in [(4_u16, 3), (3_u16, 0)] {
            let mut heap = SourceHeap::default();
            let atoms = heap
                .allocate_model_storage(vec![
                    ring_atom(&[1, 2]),
                    ring_atom(&[0, 2]),
                    ring_atom(&[0, 1]),
                ])
                .unwrap();
            let levels = heap.allocate_model_storage(vec![8_u16, 9, 10]).unwrap();
            let sources = heap.allocate_model_storage(vec![8_i8, 9, 10]).unwrap();
            let (values, queue) = ring_queue(&mut heap, vec![91, 92, 93], 3, 2, 1);
            let allocations = heap.live_allocation_count();
            assert_eq!(
                is_bond_in_Nmax_memb_ring(&mut heap, atoms, 0, 0, queue, levels, sources, maximum,),
                Ok(expected)
            );
            assert_eq!(heap.slice(levels.as_const()).unwrap(), &[0, 0, 0]);
            assert_eq!(heap.slice(sources.as_const()).unwrap(), &[0, 0, 0]);
            assert_eq!(heap.slice(values.as_const()).unwrap(), &[1, 2, 93]);
            assert_eq!(heap.live_allocation_count(), allocations);
        }

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![ring_atom(&[1, 2]), ring_atom(&[0]), ring_atom(&[0])])
            .unwrap();
        let levels = heap.allocate_model_storage(vec![0_u16; 3]).unwrap();
        let sources = heap.allocate_model_storage(vec![0_i8; 3]).unwrap();
        let (_values, queue) = ring_queue(&mut heap, vec![0, 0, 0], 3, 0, 0);
        assert_eq!(
            is_bond_in_Nmax_memb_ring(&mut heap, atoms, 0, 1, queue, levels, sources, 8),
            Ok(0)
        );
        assert_eq!(heap.slice(levels.as_const()).unwrap(), &[0, 0, 0]);
        assert_eq!(heap.slice(sources.as_const()).unwrap(), &[0, 0, 0]);

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![ring_atom(&[1, 2]), ring_atom(&[]), ring_atom(&[])])
            .unwrap();
        let levels = heap.allocate_model_storage(vec![0_u16; 3]).unwrap();
        let sources = heap.allocate_model_storage(vec![0_i8; 3]).unwrap();
        let (values, queue) = ring_queue(&mut heap, vec![77], 1, 0, 0);
        assert_eq!(
            is_bond_in_Nmax_memb_ring(&mut heap, atoms, 0, 0, queue, levels, sources, 8),
            Ok(0)
        );
        assert_eq!(heap.slice(values.as_const()).unwrap(), &[1]);
        assert_eq!(heap.slice(levels.as_const()).unwrap(), &[0, 2, 2]);
        assert_eq!(heap.slice(sources.as_const()).unwrap(), &[0, 2, 1]);
    }

    #[test]
    fn source_port__ichiring__is_atom_in_3memb_ring__line_420() {
        fn evaluate(atoms: Vec<inp_ATOM>, at_no: i32) -> Result<i32, SourceHeapError> {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            is_atom_in_3memb_ring(&heap, atoms, at_no)
        }

        let mut too_small = inp_ATOM::default();
        too_small.nNumAtInRingSystem = 2;
        too_small.valence = 1;
        too_small.neighbor[0] = u16::MAX;
        assert_eq!(evaluate(vec![too_small], 0), Ok(0));

        let mut triangle = vec![ring_atom(&[1, 2]), ring_atom(&[0, 2]), ring_atom(&[0, 1])];
        for atom in &mut triangle {
            atom.nRingSystem = 7;
            atom.nNumAtInRingSystem = 3;
        }
        assert_eq!(evaluate(triangle.clone(), 0), Ok(1));
        assert_eq!(evaluate(triangle.clone(), 1), Ok(1));
        assert_eq!(evaluate(triangle.clone(), 2), Ok(1));

        let mut open = triangle.clone();
        open[1] = ring_atom(&[0]);
        open[1].nRingSystem = 7;
        open[1].nNumAtInRingSystem = 3;
        open[2] = ring_atom(&[0]);
        open[2].nRingSystem = 7;
        open[2].nNumAtInRingSystem = 3;
        assert_eq!(evaluate(open, 0), Ok(0));

        let mut split_system = triangle;
        split_system[1].nRingSystem = 8;
        split_system[2].nRingSystem = 9;
        assert_eq!(evaluate(split_system, 0), Ok(0));

        let mut dangling = inp_ATOM::default();
        dangling.nNumAtInRingSystem = 3;
        dangling.valence = 1;
        dangling.neighbor[0] = u16::MAX;
        assert_eq!(
            evaluate(vec![dangling], 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            evaluate(vec![inp_ATOM::default()], -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            evaluate(vec![inp_ATOM::default()], 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichiring__queuedelete__line_203() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            QueueDelete(&mut heap, SourceMutPointer::null()),
            Ok(SourceMutPointer::null())
        );
        assert_eq!(heap.live_allocation_count(), 0);

        let queue = heap.allocate_model_storage(vec![QUEUE::default()]).unwrap();
        assert_eq!(heap.live_allocation_count(), 1);
        assert_eq!(QueueDelete(&mut heap, queue), Ok(SourceMutPointer::null()));
        assert_eq!(heap.live_allocation_count(), 0);

        let values = heap.allocate_model_storage(vec![7_u16, 11]).unwrap();
        let queue = heap
            .allocate_model_storage(vec![QUEUE {
                Val: values,
                nTotLength: 2,
                nFirst: 0,
                nLength: 2,
            }])
            .unwrap();
        assert_eq!(heap.live_allocation_count(), 2);
        assert_eq!(QueueDelete(&mut heap, queue), Ok(SourceMutPointer::null()));
        assert_eq!(heap.live_allocation_count(), 0);
    }
}
