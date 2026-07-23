use crate::source::base::util::{inchi_calloc, inchi_free};
use crate::source_types::{
    AT_NUMB, AT_RANK, NEIGH_LIST, S_CHAR, SourceConstPointer, SourceHeap, SourceHeapError,
    SourceMutPointer,
};

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
pub(crate) fn iisort(
    heap: &mut SourceHeap,
    list: SourceMutPointer<i32>,
    number: i32,
) -> Result<SourceMutPointer<i32>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:1014 iisort
    // INCHI❌❌: int * iisort( int *list, int num )
    // INCHI❌❌: {
    // INCHI❌❌:     int i;
    // INCHI❌❌:     for (i = 1; i < num; i++)
    // INCHI❌❌:     {
    // INCHI❌❌:         int tmp = list[i];
    // INCHI❌❌:         int j = i - 1;
    // INCHI❌❌:         while (j >= 0 && list[j] > tmp)
    // INCHI❌❌:         {
    // INCHI❌❌:             list[j + 1] = list[j];
    // INCHI❌❌:             j--;
    // INCHI❌❌:         }
    // INCHI❌❌:         list[j + 1] = tmp;
    // INCHI❌❌:     }
    // INCHI❌❌:
    // INCHI❌❌:     return list;
    // INCHI❌❌: }
    // END INCHI C FUNCTION: iisort

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

        let null = SourceMutPointer::null();
        assert_eq!(iisort(&mut heap, null, i32::MIN), Ok(null));
        assert_eq!(iisort(&mut heap, null, 0), Ok(null));
        assert_eq!(iisort(&mut heap, null, 1), Ok(null));
        assert_eq!(
            iisort(&mut heap, descending, 6),
            Err(SourceHeapError::PointerOutOfBounds)
        );
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
