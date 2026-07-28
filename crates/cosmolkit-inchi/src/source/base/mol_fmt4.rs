use crate::source::base::ichi_bns::nBondsValenceInpAt;
use crate::source::base::ichi_io::inchi_ios_print_nodisplay;
use crate::source::base::util::{
    get_atomic_mass_from_elnum, inchi_calloc, inchi_free, is_in_the_ilist,
    needed_unusual_el_valence,
};
use crate::source_types::{
    _IS_ERROR, FILE, INCHI_IOSTREAM, INT_ARRAY, MOL_FMT_M_CONN_EU, MOL_FMT_M_CONN_HH,
    MOL_FMT_M_CONN_HT, MOL_FMT_M_SST_ALT, MOL_FMT_M_SST_BLK, MOL_FMT_M_SST_RAN, OAD_PolymerUnit,
    ORIG_ATOM_DATA, RADICAL_DOUBLET, RADICAL_SINGLET, RADICAL_TRIPLET, SourceConstPointer,
    SourceFormatArgument, SourceHeap, SourceHeapError, SourceMutPointer, SourceVaList, inp_ATOM,
};

fn mol_fmt4_print(
    heap: &mut SourceHeap,
    stream: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    format: &[i8],
    arguments: Vec<SourceFormatArgument>,
) -> Result<i32, SourceHeapError> {
    let format_pointer = heap.allocate_model_storage(format.to_vec())?;
    let result = inchi_ios_print_nodisplay(
        heap,
        stream,
        stdout,
        format_pointer.as_const(),
        &SourceVaList {
            arguments,
            position: 0,
        },
    );
    heap.free(format_pointer)?;
    result
}

#[allow(non_snake_case)]
pub(crate) fn IntArray_Alloc(
    heap: &mut SourceHeap,
    items: &mut INT_ARRAY,
    nitems: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:510 IntArray_Alloc
    // INCHI✔️❌: int IntArray_Alloc(INT_ARRAY *items, int nitems)
    // INCHI✔️❌: {
    // INCHI✔️❌:     if ((items->item = (int *)inchi_calloc(nitems, sizeof(int)))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         items->increment = items->allocated = nitems;
    // INCHI✔️❌:         items->used = 0;
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return -1;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: IntArray_Alloc
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: IntArray_Alloc
    // INCHI✔️❌: #define inchi_calloc calloc
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: IntArray_Alloc

    match inchi_calloc::<i32>(heap, nitems as u64, 4) {
        Ok(pointer) => {
            items.item = pointer;
            items.allocated = nitems;
            items.increment = nitems;
            items.used = 0;
            Ok(0)
        }
        Err(SourceHeapError::AllocationFailed)
        | Err(SourceHeapError::AllocationElementCountOutOfRange)
        | Err(SourceHeapError::AllocationSizeOverflow) => {
            items.item = SourceMutPointer::null();
            Ok(-1)
        }
        Err(error) => Err(error),
    }
}

#[allow(non_snake_case)]
pub(crate) fn IntArray_ReAlloc(
    heap: &mut SourceHeap,
    items: Option<&mut INT_ARRAY>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:525 IntArray_ReAlloc
    // INCHI✔️❌: int IntArray_ReAlloc(INT_ARRAY *items)
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (items)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (items->item && items->allocated > 0 && items->increment > 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             void *p = items->item;
    // INCHI✔️❌:             if ((items->item =
    // INCHI✔️❌:                      (int *)inchi_calloc((long long)items->allocated + (long long)items->increment, sizeof(items->item[0])))) /* djb-rwth: cast operators added; addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 memcpy(items->item, p, items->used * sizeof(items->item[0]));
    // INCHI✔️❌:                 inchi_free(p);
    // INCHI✔️❌:                 items->allocated += items->increment;
    // INCHI✔️❌:                 return 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_free(p);
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return -1;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: IntArray_ReAlloc
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: IntArray_ReAlloc
    // INCHI✔️❌: #define inchi_calloc calloc
    // INCHI✔️❌: #define inchi_free(X) do { if (X) free(X); } while (0)
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: IntArray_ReAlloc

    let Some(items) = items else {
        return Ok(-1);
    };
    if items.item.is_null() || items.allocated <= 0 || items.increment <= 0 {
        return Ok(-1);
    }

    let old = items.item;
    let new_count = i64::from(items.allocated) + i64::from(items.increment);
    let replacement = inchi_calloc::<i32>(heap, new_count as u64, 4);
    let replacement = match replacement {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed)
        | Err(SourceHeapError::AllocationElementCountOutOfRange)
        | Err(SourceHeapError::AllocationSizeOverflow) => {
            items.item = SourceMutPointer::null();
            inchi_free(heap, old)?;
            return Ok(-1);
        }
        Err(error) => return Err(error),
    };
    items.item = replacement;

    let used = usize::try_from(items.used).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let old_values = heap
        .slice(old.as_const())?
        .get(..used)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    heap.slice_mut(items.item)?
        .get_mut(..used)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .copy_from_slice(&old_values);
    inchi_free(heap, old)?;
    items.allocated = items
        .allocated
        .checked_add(items.increment)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn IntArray_Append(
    heap: &mut SourceHeap,
    items: Option<&mut INT_ARRAY>,
    new_item: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:550 IntArray_Append
    // INCHI✔️❌: int IntArray_Append(INT_ARRAY *items, int new_item)
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (items)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (items->used + 1 > items->allocated)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* need to expand buffer */
    // INCHI✔️❌:             if (IntArray_ReAlloc(items))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         items->item[items->used++] = new_item;
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return -1;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: IntArray_Append
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: IntArray_Append
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: IntArray_Append

    let Some(items) = items else {
        return Ok(-1);
    };
    let next_used = items
        .used
        .checked_add(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    if next_used > items.allocated && IntArray_ReAlloc(heap, Some(items))? != 0 {
        return Ok(-1);
    }
    let index = usize::try_from(items.used).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    *heap
        .slice_mut(items.item)?
        .get_mut(index)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = new_item;
    items.used = next_used;
    Ok(0)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn IntArray_AppendIfAbsent(
    heap: &mut SourceHeap,
    items: &mut INT_ARRAY,
    new_item: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:572 IntArray_AppendIfAbsent
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
int IntArray_AppendIfAbsent(INT_ARRAY *items, int new_item)
{
    if (!is_in_the_ilist(items->item, new_item, items->used))
    {
        return IntArray_Append(items, new_item);
    }
    return 0;
}
    */
    // END INCHI C FUNCTION: IntArray_AppendIfAbsent
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: IntArray_AppendIfAbsent
    // INCHI✔️✔️: COMPILE_ANSI_ONLY and TARGET_API_LIB do not alter this production helper.
    // INCHI✔️✔️: The Rust lookup is the already-ported linear is_in_the_ilist implementation.
    // INCHI✔️✔️: A missing item delegates directly to the already-ported IntArray_Append.
    // END INCHI ACTIVE MACRO CONFIGURATION: IntArray_AppendIfAbsent

    let item_slice = if items.used == 0 {
        None
    } else {
        Some(heap.slice(items.item.as_const())?)
    };
    if is_in_the_ilist(item_slice, new_item, items.used)?.is_none() {
        return IntArray_Append(heap, Some(items), new_item);
    }
    Ok(0)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn IntArray_DebugPrint(items: SourceMutPointer<INT_ARRAY>) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:582 IntArray_DebugPrint
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
void IntArray_DebugPrint(INT_ARRAY *items)
{
    if (items)
    {
        int i;
        if (items->used > 0)
        {
            for (i = 0; i < items->used - 1; i++)
            {
                ITRACE_("%-d, ", items->item[i]);
            }
            ITRACE_("%-d\n", items->item[items->used - 1]);
        }
        else
        {
            ; /*ITRACE_( "[None]\n");*/
        }
    }
}
    */
    // END INCHI C FUNCTION: IntArray_DebugPrint
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: IntArray_DebugPrint
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux, so ITRACE_ is 0 && _inchi_trace.
    // INCHI✔️✔️: Short-circuiting prevents evaluation of every item access and produces no output.
    // INCHI✔️✔️: The optimized active production behavior is therefore an allocation-free O(1) no-op.
    // END INCHI ACTIVE MACRO CONFIGURATION: IntArray_DebugPrint

    let _ = items;
}

#[allow(non_snake_case)]
pub(crate) fn IntArray_Reset(items: &mut INT_ARRAY) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:603 IntArray_Reset
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
    void IntArray_Reset(INT_ARRAY *items)
    {
        items->used = 0;
        return;
    }
        */
    // END INCHI C FUNCTION: IntArray_Reset
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: IntArray_Reset
    // INCHI✔️✔️: COMPILE_ANSI_ONLY and TARGET_API_LIB do not alter this production helper.
    // INCHI✔️✔️: The allocation pointer, capacity, increment, and element bytes remain untouched.
    // INCHI✔️✔️: Rust performs the same single field write with no allocation.
    // END INCHI ACTIVE MACRO CONFIGURATION: IntArray_Reset

    items.used = 0;
}

#[allow(non_snake_case)]
pub(crate) fn IntArray_Free(
    heap: &mut SourceHeap,
    items: Option<&mut INT_ARRAY>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:610 IntArray_Free
    // INCHI✔️❌: void IntArray_Free(INT_ARRAY *items)
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (items)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (items->item)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free(items->item);
    // INCHI✔️❌:         }
    // INCHI✔️❌:         memset(items, 0, sizeof(*items)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: IntArray_Free
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: IntArray_Free
    // INCHI✔️❌: #define inchi_free(X) do { if (X) free(X); } while (0)
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: IntArray_Free

    let Some(items) = items else {
        return Ok(());
    };
    if !items.item.is_null() {
        inchi_free(heap, items.item)?;
    }
    *items = INT_ARRAY::default();
    Ok(())
}

fn mol_fmt4_c_text(
    heap: &SourceHeap,
    pointer: SourceConstPointer<i8>,
) -> Result<Vec<u8>, SourceHeapError> {
    if pointer.is_null() {
        return Ok(Vec::new());
    }
    let bytes = heap.slice(pointer)?;
    let length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    Ok(bytes[..length].iter().map(|byte| *byte as u8).collect())
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OrigAtData_WriteToSDfile(
    heap: &mut SourceHeap,
    input: &ORIG_ATOM_DATA,
    mut stream: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    name: SourceConstPointer<i8>,
    comment: SourceConstPointer<i8>,
    chiral_flag: i32,
    atoms_dt: i32,
    label: SourceConstPointer<i8>,
    value: SourceConstPointer<i8>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:783 OrigAtData_WriteToSDfile
    // INCHI✔️❌: int OrigAtData_WriteToSDfile(const ORIG_ATOM_DATA *inp_at_data,
    // INCHI✔️❌:                              INCHI_IOSTREAM *fcb,
    // INCHI✔️❌:                              const char *name,
    // INCHI✔️❌:                              const char *comment,
    // INCHI✔️❌:                              int bChiralFlag,
    // INCHI✔️❌:                              int bAtomsDT,
    // INCHI✔️❌:                              const char *szLabel,
    // INCHI✔️❌:                              const char *szValue)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int num_bonds = 0, nNumAddLines = 0, nNumIsoLines = 0, nNumChargeLines = 0,
    // INCHI✔️❌:         nNumRadicalLines = 0, nNumAliasLines = 0, ret = 0;
    // INCHI✔️❌:     INT_ARRAY written_bond_ends;
    // INCHI✔️❌:
    // INCHI✔️❌:     OrigAtData_WriteToSDfileHeaderAndCountThings((ORIG_ATOM_DATA *)inp_at_data,
    // INCHI✔️❌:                                                  fcb, name, comment,
    // INCHI✔️❌:                                                  bChiralFlag, bAtomsDT,
    // INCHI✔️❌:                                                  szLabel, szValue,
    // INCHI✔️❌:                                                  &nNumAliasLines,
    // INCHI✔️❌:                                                  &nNumChargeLines,
    // INCHI✔️❌:                                                  &nNumRadicalLines,
    // INCHI✔️❌:                                                  &nNumIsoLines,
    // INCHI✔️❌:                                                  &nNumAddLines,
    // INCHI✔️❌:                                                  &num_bonds);
    // INCHI✔️❌:
    // INCHI✔️❌:     if (IntArray_Alloc(&written_bond_ends, num_bonds ? num_bonds : 255))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = _IS_ERROR;
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     OrigAtData_WriteToSDfileAtomsBlock(inp_at_data, fcb, name, comment,
    // INCHI✔️❌:                                        bAtomsDT, szLabel, szValue);
    // INCHI✔️❌:
    // INCHI✔️❌:     OrigAtData_WriteToSDfileBondsBlock(inp_at_data, fcb, name, comment,
    // INCHI✔️❌:                                        szLabel, szValue, &written_bond_ends);
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nNumAddLines)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         OrigAtData_WriteToSDfileAdditionalLines(inp_at_data, fcb, name, comment,
    // INCHI✔️❌:                                                 bAtomsDT, szLabel, szValue,
    // INCHI✔️❌:                                                 nNumAliasLines, nNumChargeLines,
    // INCHI✔️❌:                                                 nNumRadicalLines, nNumIsoLines,
    // INCHI✔️❌:                                                 &written_bond_ends);
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Add field with label/ID if applicable and mark the end of record */
    // INCHI✔️❌:     if (szValue && szValue[0])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (szLabel && szLabel[0])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_print_nodisplay(fcb, "> <%s>\n", szLabel);
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_print_nodisplay(fcb, "> <ID>\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_ios_print_nodisplay(fcb, " %s\n\n", szValue);
    // INCHI✔️❌:     }
    // INCHI✔️❌:     inchi_ios_print_nodisplay(fcb, "$$$$\n");
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     IntArray_Free(&written_bond_ends);
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OrigAtData_WriteToSDfile
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_WriteToSDfile
    // INCHI✔️❌: #define _IS_ERROR 2
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_WriteToSDfile

    let mut num_bonds = 0_i32;
    let mut num_add_lines = 0_i32;
    let mut num_iso_lines = 0_i32;
    let mut num_charge_lines = 0_i32;
    let mut num_radical_lines = 0_i32;
    let mut num_alias_lines = 0_i32;
    let _ = OrigAtData_WriteToSDfileHeaderAndCountThings(
        heap,
        input,
        stream.as_deref_mut(),
        stdout,
        name,
        comment,
        chiral_flag,
        atoms_dt,
        label,
        value,
        &mut num_alias_lines,
        &mut num_charge_lines,
        &mut num_radical_lines,
        &mut num_iso_lines,
        &mut num_add_lines,
        &mut num_bonds,
    )?;

    let mut written_bond_ends = INT_ARRAY::default();
    if IntArray_Alloc(
        heap,
        &mut written_bond_ends,
        if num_bonds != 0 { num_bonds } else { 255 },
    )? != 0
    {
        IntArray_Free(heap, Some(&mut written_bond_ends))?;
        return Ok(_IS_ERROR as i32);
    }

    let operation = (|| -> Result<i32, SourceHeapError> {
        let _ = OrigAtData_WriteToSDfileAtomsBlock(
            heap,
            input,
            stream.as_deref_mut(),
            stdout,
            name,
            comment,
            atoms_dt,
            label,
            value,
        )?;
        let _ = OrigAtData_WriteToSDfileBondsBlock(
            heap,
            input,
            stream.as_deref_mut(),
            stdout,
            name,
            comment,
            label,
            value,
            Some(&mut written_bond_ends),
        )?;
        if num_add_lines != 0 {
            let _ = OrigAtData_WriteToSDfileAdditionalLines(
                heap,
                Some(input),
                stream.as_deref_mut(),
                stdout,
                name,
                comment,
                atoms_dt,
                label,
                value,
                num_alias_lines,
                num_charge_lines,
                num_radical_lines,
                num_iso_lines,
                Some(&written_bond_ends),
            )?;
        }

        let value_bytes = mol_fmt4_c_text(heap, value)?;
        if !value_bytes.is_empty() {
            let label_bytes = mol_fmt4_c_text(heap, label)?;
            if !label_bytes.is_empty() {
                let format = (*b"> <%s>\n\0").map(|byte| byte as i8);
                mol_fmt4_print_dynamic_bytes(
                    heap,
                    stream.as_deref_mut(),
                    stdout,
                    &format,
                    vec![],
                    &label_bytes,
                )?;
            } else {
                let format = (*b"> <ID>\n\0").map(|byte| byte as i8);
                let _ = mol_fmt4_print(heap, stream.as_deref_mut(), stdout, &format, vec![])?;
            }
            let format = (*b" %s\n\n\0").map(|byte| byte as i8);
            mol_fmt4_print_dynamic_bytes(
                heap,
                stream.as_deref_mut(),
                stdout,
                &format,
                vec![],
                &value_bytes,
            )?;
        }
        let format = (*b"$$$$\n\0").map(|byte| byte as i8);
        let _ = mol_fmt4_print(heap, stream.as_deref_mut(), stdout, &format, vec![])?;
        Ok(0)
    })();
    let cleanup = IntArray_Free(heap, Some(&mut written_bond_ends));
    match (operation, cleanup) {
        (Ok(result), Ok(())) => Ok(result),
        (Err(error), Ok(())) | (Ok(_), Err(error)) => Err(error),
        (Err(error), Err(_cleanup_error)) => Err(error),
    }
}

fn mol_fmt4_header_text(
    heap: &SourceHeap,
    pointer: SourceConstPointer<i8>,
) -> Result<Vec<u8>, SourceHeapError> {
    if pointer.is_null() {
        return Ok(Vec::new());
    }
    let bytes = heap.slice(pointer)?;
    let Some(first) = bytes.first() else {
        return Err(SourceHeapError::PointerOutOfBounds);
    };
    if *first == 0 {
        return Ok(Vec::new());
    }
    Ok(bytes
        .iter()
        .take(80)
        .take_while(|byte| **byte != 0)
        .map(|byte| *byte as u8)
        .collect())
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OrigAtData_WriteToSDfileHeaderAndCountThings(
    heap: &mut SourceHeap,
    input: &ORIG_ATOM_DATA,
    mut stream: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    name: SourceConstPointer<i8>,
    comment: SourceConstPointer<i8>,
    chiral_flag: i32,
    atoms_dt: i32,
    _label: SourceConstPointer<i8>,
    _value: SourceConstPointer<i8>,
    num_alias_lines: &mut i32,
    num_charge_lines: &mut i32,
    num_radical_lines: &mut i32,
    num_iso_lines: &mut i32,
    num_add_lines: &mut i32,
    num_bonds: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:852 OrigAtData_WriteToSDfileHeaderAndCountThings
    // INCHI✔️❌: int OrigAtData_WriteToSDfileHeaderAndCountThings(const ORIG_ATOM_DATA *inp_at_data,
    // INCHI✔️❌:                                                  INCHI_IOSTREAM *fcb,
    // INCHI✔️❌:                                                  const char *name,
    // INCHI✔️❌:                                                  const char *comment,
    // INCHI✔️❌:                                                  int bChiralFlag,
    // INCHI✔️❌:                                                  int bAtomsDT,
    // INCHI✔️❌:                                                  const char *szLabel,
    // INCHI✔️❌:                                                  const char *szValue,
    // INCHI✔️❌:                                                  int *nNumAliasLines,
    // INCHI✔️❌:                                                  int *nNumChargeLines,
    // INCHI✔️❌:                                                  int *nNumRadicalLines,
    // INCHI✔️❌:                                                  int *nNumIsoLines,
    // INCHI✔️❌:                                                  int *nNumAddLines,
    // INCHI✔️❌:                                                  int *num_bonds)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, ret = 0;
    // INCHI✔️❌:     int bAtomNeedsAlias,
    // INCHI✔️❌:         nNumNecessaryIsoLines = 0,
    // INCHI✔️❌:         nNumNecessaryChgLines = 0,
    // INCHI✔️❌:         nNumNecessaryRadLines = 0;
    // INCHI✔️❌:     int num_atoms = inp_at_data->num_inp_atoms;
    // INCHI✔️❌:     int bV2000 = SDF_OUTPUT_V2000;
    // INCHI✔️❌:     const inp_ATOM *at = inp_at_data->at;
    // INCHI✔️❌:
    // INCHI✔️❌:     {
    // INCHI✔️❌:         char strLocName[82];
    // INCHI✔️❌:         memset(strLocName, 0, sizeof(strLocName)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         if (name && *name)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             strncpy(strLocName, name, 80);
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_ios_print_nodisplay(fcb, "%s\n", strLocName);
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /**********************************************************************/
    // INCHI✔️❌:     /**                                                                  **/
    // INCHI✔️❌:     /** Important: Atoms with alias cannot have charge, radical          **/
    // INCHI✔️❌:     /**            isotope differences are allowed                       **/
    // INCHI✔️❌:     /**                                                                  **/
    // INCHI✔️❌:     /**            Atoms with alias cannot be abnormal.                  **/
    // INCHI✔️❌:     /**                                                                  **/
    // INCHI✔️❌:     /** Abnormal atoms are atoms which need M  CHG, M RAD, M  ISO        **/
    // INCHI✔️❌:     /**                                                                  **/
    // INCHI✔️❌:     /** Output aliased atoms if they have implicit D or T                **/
    // INCHI✔️❌:     /**                                                                  **/
    // INCHI✔️❌:     /**********************************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:     /*                                    F10.5     F12.5       I6
    // INCHI✔️❌:                      IIPPPPPPPPMMDDYYHHmmddSSssssssssssEEEEEEEEEEEERRRRRR
    // INCHI✔️❌:     inchi_ios_eprint( fcb,"NISTTRANHP09089809272D 1   1.0         0.0    %6ld\n", lEpa);*/
    // INCHI✔️❌:     /*^^^
    // INCHI✔️❌:     inchi_ios_print_nodisplay( fcb,"  %s v%s SDfile Output                       \n", INCHI_NAME, INCHI_VERSION);
    // INCHI✔️❌:
    // INCHI✔️❌:     Changed 01/10/2009 to conform CTFile specification (by Symyx request)*/
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_ios_print_nodisplay(fcb,
    // INCHI✔️❌:                               /*   IIPPPPPPPPMMDDYYHHmmddSSssssssssssEEEEEEEEEEEERRRRRR*/
    // INCHI✔️❌:                               "  InChIV10                                     \n");
    // INCHI✔️❌:     /*y_fprintf(fcb, "  -CPSS-  1213981200n\n");*/
    // INCHI✔️❌:
    // INCHI✔️❌:     {
    // INCHI✔️❌:         char strLocName[82];
    // INCHI✔️❌:
    // INCHI✔️❌:         memset(strLocName, 0, sizeof(strLocName)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         if (comment && *comment)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             strncpy(strLocName, comment, 80);
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_ios_print_nodisplay(fcb, "%s\n", strLocName);
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     *num_bonds = 0;
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         (*num_bonds) += at[i].valence;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     (*num_bonds) /= 2;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*find if we need "M  CHG" and "M  RAD"*/
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((bAtomNeedsAlias = ALIASED_AT(i))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* has isotopic implicit D or T; ignoring pure 1H */
    // INCHI✔️❌:             (*nNumAliasLines) += 2 * bAtomNeedsAlias;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* abnormal means atom needs CHG, RAD, or ISO entry */
    // INCHI✔️❌:             /* nNumAddLines    += ABNORMAL_AT(i); */
    // INCHI✔️❌:             /* nNumIso         += ( 0 == strcmp( at[i].elname, "D" ) || ( 0 == strcmp( at[i].elname, "T" ) || at[i].iso_atw_diff ) ); */
    // INCHI✔️❌:             /* nNumAddIso      += at[i].iso_atw_diff && (at[i].iso_atw_diff == 1 || at[i].iso_atw_diff < -3 || at[i].iso_atw_diff > 5 ); */
    // INCHI✔️❌:             nNumNecessaryIsoLines += ABNORMAL_ISO(i);
    // INCHI✔️❌:             nNumNecessaryChgLines += ABNORMAL_CHG(i);
    // INCHI✔️❌:             nNumNecessaryRadLines += ABNORMAL_RAD(i);
    // INCHI✔️❌:             (*nNumIsoLines) += ANY_ISO(i, bAtomsDT);
    // INCHI✔️❌:             (*nNumChargeLines) += ANY_CHG(i);
    // INCHI✔️❌:             (*nNumRadicalLines) += ANY_RAD(i);
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     *nNumChargeLines = (*nNumChargeLines + 7) / 8;
    // INCHI✔️❌:     *nNumRadicalLines = (*nNumRadicalLines + 7) / 8;
    // INCHI✔️❌:     *nNumIsoLines = (*nNumIsoLines + 7) / 8;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!bV2000)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (!nNumNecessaryRadLines && !nNumNecessaryChgLines)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *nNumRadicalLines = 0;
    // INCHI✔️❌:             *nNumChargeLines = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (!nNumNecessaryIsoLines)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *nNumIsoLines = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* recalculate number of added lines */
    // INCHI✔️❌:     *nNumAddLines = *nNumChargeLines + *nNumRadicalLines + *nNumIsoLines + *nNumAliasLines; /* 1 for M  END*/
    // INCHI✔️❌:
    // INCHI✔️❌:     if (*nNumAddLines || bV2000)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *nNumAddLines += 1; /* add 1 for "M  END" line*/
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*                         aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv                                      */
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_ios_print_nodisplay(fcb, "%3d%3d  0  0%3d  0  0  0  0  0%3d%s\n",
    // INCHI✔️❌:                               num_atoms, *num_bonds, bChiralFlag ? 1 : 0, *nNumAddLines, *nNumAddLines ? " V2000" : "");
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OrigAtData_WriteToSDfileHeaderAndCountThings
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_WriteToSDfileHeaderAndCountThings
    // INCHI✔️❌: #define SDF_OUTPUT_V2000 1
    // INCHI✔️❌: #define ALIASED_AT(i) (0 < NUM_ISO_H(at, i))
    // INCHI✔️❌: #define NUM_ISO_H(AT, N) (AT[N].num_iso_H[0] + AT[N].num_iso_H[1] + AT[N].num_iso_H[2])
    // INCHI✔️❌: #define IS_DEUTERIUM(i) (!strcmp(at[i].elname, "D") || at[i].iso_atw_diff == 2 && !strcmp(at[i].elname, "H"))
    // INCHI✔️❌: #define IS_TRITIUM(i) (!strcmp(at[i].elname, "T") || at[i].iso_atw_diff == 3 && !strcmp(at[i].elname, "H"))
    // INCHI✔️❌: #define ABNORMAL_ISO(i) (at[i].iso_atw_diff == 1 || at[i].iso_atw_diff < -3 || at[i].iso_atw_diff > 5)
    // INCHI✔️❌: #define ABNORMAL_CHG(i) (abs(at[i].charge) > 3)
    // INCHI✔️❌: #define ABNORMAL_RAD(i) (RADICAL_SINGLET <= at[i].radical && at[i].radical <= RADICAL_TRIPLET)
    // INCHI✔️❌: #define ANY_ISO(i, X) ((X) ? (at[i].iso_atw_diff && !IS_DEUTERIUM(i) && !IS_TRITIUM(i)) : (at[i].iso_atw_diff || IS_DEUTERIUM(i) || IS_TRITIUM(i)))
    // INCHI✔️❌: #define ANY_CHG(i) (0 != at[i].charge)
    // INCHI✔️❌: #define ANY_RAD(i) (RADICAL_SINGLET <= at[i].radical && at[i].radical <= RADICAL_TRIPLET)
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_WriteToSDfileHeaderAndCountThings

    let name = mol_fmt4_header_text(heap, name)?;
    let format = (*b"%s\n\0").map(|byte| byte as i8);
    mol_fmt4_print_dynamic_bytes(heap, stream.as_deref_mut(), stdout, &format, vec![], &name)?;
    let format = (*b"  InChIV10                                     \n\0").map(|byte| byte as i8);
    let _ = mol_fmt4_print(heap, stream.as_deref_mut(), stdout, &format, vec![])?;
    let comment = mol_fmt4_header_text(heap, comment)?;
    let format = (*b"%s\n\0").map(|byte| byte as i8);
    mol_fmt4_print_dynamic_bytes(
        heap,
        stream.as_deref_mut(),
        stdout,
        &format,
        vec![],
        &comment,
    )?;

    let atom_count =
        usize::try_from(input.num_inp_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atoms = if atom_count == 0 {
        &[][..]
    } else {
        heap.slice(input.at.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
    };
    *num_bonds = 0;
    for atom in atoms {
        *num_bonds = num_bonds.wrapping_add(i32::from(atom.valence));
    }
    *num_bonds /= 2;

    let mut necessary_iso_lines = 0_i32;
    let mut necessary_charge_lines = 0_i32;
    let mut necessary_radical_lines = 0_i32;
    for atom in atoms {
        let alias = atom
            .num_iso_H
            .iter()
            .map(|value| i32::from(*value))
            .sum::<i32>()
            > 0;
        if alias {
            *num_alias_lines = num_alias_lines.wrapping_add(2);
            continue;
        }
        let atom_name = mol_fmt4_atom_name(atom)?;
        let deuterium = atom_name == "D" || (atom.iso_atw_diff == 2 && atom_name == "H");
        let tritium = atom_name == "T" || (atom.iso_atw_diff == 3 && atom_name == "H");
        necessary_iso_lines = necessary_iso_lines.wrapping_add(i32::from(
            atom.iso_atw_diff == 1 || atom.iso_atw_diff < -3 || atom.iso_atw_diff > 5,
        ));
        necessary_charge_lines =
            necessary_charge_lines.wrapping_add(i32::from(i32::from(atom.charge).abs() > 3));
        let valid_radical =
            (RADICAL_SINGLET as i32..=RADICAL_TRIPLET as i32).contains(&i32::from(atom.radical));
        necessary_radical_lines = necessary_radical_lines.wrapping_add(i32::from(valid_radical));
        let any_iso = if atoms_dt != 0 {
            atom.iso_atw_diff != 0 && !deuterium && !tritium
        } else {
            atom.iso_atw_diff != 0 || deuterium || tritium
        };
        *num_iso_lines = num_iso_lines.wrapping_add(i32::from(any_iso));
        *num_charge_lines = num_charge_lines.wrapping_add(i32::from(atom.charge != 0));
        *num_radical_lines = num_radical_lines.wrapping_add(i32::from(valid_radical));
    }

    *num_charge_lines = num_charge_lines.wrapping_add(7) / 8;
    *num_radical_lines = num_radical_lines.wrapping_add(7) / 8;
    *num_iso_lines = num_iso_lines.wrapping_add(7) / 8;

    // SDF_OUTPUT_V2000 is 1 in the selected production configuration.
    let _inactive_v2000_counters = (
        necessary_iso_lines,
        necessary_charge_lines,
        necessary_radical_lines,
    );
    *num_add_lines = num_charge_lines
        .wrapping_add(*num_radical_lines)
        .wrapping_add(*num_iso_lines)
        .wrapping_add(*num_alias_lines);
    *num_add_lines = num_add_lines.wrapping_add(1);

    let format = (*b"%3d%3d  0  0%3d  0  0  0  0  0%3d%s\n\0").map(|byte| byte as i8);
    mol_fmt4_print_dynamic_bytes(
        heap,
        stream.as_deref_mut(),
        stdout,
        &format,
        vec![
            SourceFormatArgument::Signed(i64::from(input.num_inp_atoms)),
            SourceFormatArgument::Signed(i64::from(*num_bonds)),
            SourceFormatArgument::Signed(i64::from(i32::from(chiral_flag != 0))),
            SourceFormatArgument::Signed(i64::from(*num_add_lines)),
        ],
        b" V2000",
    )?;
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OrigAtData_WriteToSDfileAtomsBlock(
    heap: &mut SourceHeap,
    input: &ORIG_ATOM_DATA,
    mut stream: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    _name: SourceConstPointer<i8>,
    _comment: SourceConstPointer<i8>,
    atoms_dt: i32,
    _label: SourceConstPointer<i8>,
    _value: SourceConstPointer<i8>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:989 OrigAtData_WriteToSDfileAtomsBlock
    // INCHI✔️❌: int OrigAtData_WriteToSDfileAtomsBlock(const ORIG_ATOM_DATA *inp_at_data,
    // INCHI✔️❌:                                         INCHI_IOSTREAM       *fcb,
    // INCHI✔️❌:                                         const char           *name,
    // INCHI✔️❌:                                         const char           *comment,
    // INCHI✔️❌:                                         int                  bAtomsDT,
    // INCHI✔️❌:                                         const char           *szLabel,
    // INCHI✔️❌:                                         const char           *szValue)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, ret = 0;
    // INCHI✔️❌:     int bAtomNeedsAlias;
    // INCHI✔️❌:     /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     int num_atoms = inp_at_data->num_inp_atoms;
    // INCHI✔️❌:     const inp_ATOM *at = inp_at_data->at;
    // INCHI✔️❌:     double x, y, z;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         char elname[ATOM_EL_LEN] = "\0\0\0\0\0";
    // INCHI✔️❌:         int iso = 0;
    // INCHI✔️❌:         int charge = 0;
    // INCHI✔️❌:         int valence = 0;
    // INCHI✔️❌:         int nIsotopeH = IS_DEUTERIUM( i ) ? 1 : IS_TRITIUM( i ) ? 2 : 0;
    // INCHI✔️❌:         int bonds_val;
    // INCHI✔️❌:         bAtomNeedsAlias = ALIASED_AT( i );
    // INCHI✔️❌:         memset( elname, 0, sizeof( elname ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:         if (bAtomNeedsAlias)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* alias */
    // INCHI✔️❌:             strcpy(elname, "C");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* isotope*/
    // INCHI✔️❌:             if (nIsotopeH)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 strcpy(elname, bAtomsDT ? (nIsotopeH == 1 ? "D" : "T") : "H");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 strncpy(elname, at[i].elname, sizeof(elname) - 1);
    // INCHI✔️❌:                 elname[sizeof(elname) - 1] = 0; /* adding zero termination after strncpy */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (!ABNORMAL_CHG(i) && !ANY_RAD(i))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* charge*/
    // INCHI✔️❌:                 /* Only atoms without alias can be here*/
    // INCHI✔️❌:                 switch (at[i].charge)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                 case 3:
    // INCHI✔️❌:                     charge = 1;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case 2:
    // INCHI✔️❌:                     charge = 2;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case 1:
    // INCHI✔️❌:                     charge = 3;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case -1:
    // INCHI✔️❌:                     charge = 5;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case -2:
    // INCHI✔️❌:                     charge = 6;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case -3:
    // INCHI✔️❌:                     charge = 7;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case 0:
    // INCHI✔️❌:                     charge = 0;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 default:
    // INCHI✔️❌:                     break; /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* radical*/
    // INCHI✔️❌:             if (ANY_RAD(i) && !ANY_CHG(i))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at[i].radical == RADICAL_DOUBLET)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     charge = 4;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* allow isotopic shift for aliased atoms */
    // INCHI✔️❌:         if (NORMAL_ISO(i, bAtomsDT))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             iso = at[i].iso_atw_diff > 0 ? at[i].iso_atw_diff - 1 :
    // INCHI✔️❌:                     at[i].iso_atw_diff < 0 ? at[i].iso_atw_diff :
    // INCHI✔️❌:                         nIsotopeH ? nIsotopeH : 0; /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         x = at[i].x;
    // INCHI✔️❌:         y = at[i].y;
    // INCHI✔️❌:         z = at[i].z;
    // INCHI✔️❌:
    // INCHI✔️❌:         /* valence -- set only if needed */
    // INCHI✔️❌:         bonds_val = nBondsValenceInpAt(at + i, NULL, NULL);
    // INCHI✔️❌:         valence = needed_unusual_el_valence(at[i].el_number, at[i].charge, at[i].radical,
    // INCHI✔️❌:                                             at[i].chem_bonds_valence, bonds_val, NUMH(at, i), at[i].valence);
    // INCHI✔️❌:
    // INCHI✔️❌:         if (valence < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             valence = 15; /* means no bonds nor H */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Convert "Zz" to "*" element symbol */
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!strcmp(elname, "Zz") || !strcmp(elname, "Zy"))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             strcpy(elname, "*");
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /*inchi_ios_eprint(fcb,"%10.4f%10.4f%10.4f %-3.3s%2d%3d  0     0  0  0  0  0  0\n",*/
    // INCHI✔️❌:         /*    (float)at[i].x, (float)(-at[i].y), fzero, at[i].elname, iso, charge);*/
    // INCHI✔️❌:         /*              xxxxxxyyyyyyzzzzzz aaa____ddcccsssnnnbbbvvvrrriiimmmeee  */
    // INCHI✔️❌:         inchi_ios_print_nodisplay(fcb, "%10.4f%10.4f%10.4f %-3.3s%2d%3d  0     0%3d  0  0  0  0\n",
    // INCHI✔️❌:                                   x, y, z, elname, (int)iso, (int)charge, valence /* at[i].special*/);
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Reflect image against x-axis;                                    */
    // INCHI✔️❌:         /* when transforming MOLfile back to STDATA in mol_to_stdata(...),  */
    // INCHI✔️❌:         /* make one more reflection to restore original orientation.        */
    // INCHI✔️❌:         /* Reason: in MS Search y-axis is directed from top to bottom,      */
    // INCHI✔️❌:         /*         while in MOLfile y-axis goes from bottom to top.         */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OrigAtData_WriteToSDfileAtomsBlock
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_WriteToSDfileAtomsBlock
    // INCHI✔️❌: #define ATOM_EL_LEN 6
    // INCHI✔️❌: #define ALIASED_AT(i) (0 < NUM_ISO_H(at, i))
    // INCHI✔️❌: #define NUM_ISO_H(AT, N) (AT[N].num_iso_H[0] + AT[N].num_iso_H[1] + AT[N].num_iso_H[2])
    // INCHI✔️❌: #define NUMH(AT, N) (AT[N].num_H + NUM_ISO_H(AT, N))
    // INCHI✔️❌: #define NORMAL_ISO(i, X) (ANY_ISO(i, X) && !ABNORMAL_ISO(i))
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_WriteToSDfileAtomsBlock

    let atom_count =
        usize::try_from(input.num_inp_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atoms = if atom_count == 0 {
        Vec::new()
    } else {
        heap.slice(input.at.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };
    for atom in &atoms {
        let atom_name = mol_fmt4_atom_name(atom)?;
        let deuterium = atom_name == "D" || (atom.iso_atw_diff == 2 && atom_name == "H");
        let tritium = atom_name == "T" || (atom.iso_atw_diff == 3 && atom_name == "H");
        let isotope_hydrogen = if deuterium {
            1
        } else if tritium {
            2
        } else {
            0
        };
        let alias = atom
            .num_iso_H
            .iter()
            .map(|value| i32::from(*value))
            .sum::<i32>()
            > 0;
        let mut element = if alias {
            b"C".to_vec()
        } else if isotope_hydrogen != 0 {
            if atoms_dt != 0 {
                if isotope_hydrogen == 1 {
                    b"D".to_vec()
                } else {
                    b"T".to_vec()
                }
            } else {
                b"H".to_vec()
            }
        } else {
            atom.elname
                .iter()
                .take(5)
                .take_while(|byte| **byte != 0)
                .map(|byte| *byte as u8)
                .collect()
        };
        let valid_radical =
            (RADICAL_SINGLET as i32..=RADICAL_TRIPLET as i32).contains(&i32::from(atom.radical));
        let mut molfile_charge = 0_i32;
        if !alias {
            if i32::from(atom.charge).abs() <= 3 && !valid_radical {
                molfile_charge = match atom.charge {
                    3 => 1,
                    2 => 2,
                    1 => 3,
                    -1 => 5,
                    -2 => 6,
                    -3 => 7,
                    _ => 0,
                };
            }
            if valid_radical && atom.charge == 0 && atom.radical == RADICAL_DOUBLET as i8 {
                molfile_charge = 4;
            }
        }
        let any_iso = if atoms_dt != 0 {
            atom.iso_atw_diff != 0 && !deuterium && !tritium
        } else {
            atom.iso_atw_diff != 0 || deuterium || tritium
        };
        let abnormal_iso =
            atom.iso_atw_diff == 1 || atom.iso_atw_diff < -3 || atom.iso_atw_diff > 5;
        let isotope = if any_iso && !abnormal_iso {
            if atom.iso_atw_diff > 0 {
                i32::from(atom.iso_atw_diff) - 1
            } else if atom.iso_atw_diff < 0 {
                i32::from(atom.iso_atw_diff)
            } else {
                isotope_hydrogen
            }
        } else {
            0
        };
        let actual_bonds_valence = nBondsValenceInpAt(atom, None, None);
        let total_hydrogens = i32::from(atom.num_H)
            .wrapping_add(i32::from(atom.num_iso_H[0]))
            .wrapping_add(i32::from(atom.num_iso_H[1]))
            .wrapping_add(i32::from(atom.num_iso_H[2]));
        let mut valence = needed_unusual_el_valence(
            i32::from(atom.el_number),
            i32::from(atom.charge),
            i32::from(atom.radical),
            i32::from(atom.chem_bonds_valence),
            actual_bonds_valence,
            total_hydrogens,
            i32::from(atom.valence),
        )?;
        if valence < 0 {
            valence = 15;
        }
        if element == b"Zz" || element == b"Zy" {
            element.clear();
            element.push(b'*');
        }

        let element_pointer = heap.allocate_model_storage(
            element
                .iter()
                .copied()
                .chain(std::iter::once(0))
                .map(|byte| byte as i8)
                .collect(),
        )?;
        let format = (*b"%10.4f%10.4f%10.4f %-3.3s%2d%3d  0     0%3d  0  0  0  0\n\0")
            .map(|byte| byte as i8);
        let result = mol_fmt4_print(
            heap,
            stream.as_deref_mut(),
            stdout,
            &format,
            vec![
                SourceFormatArgument::Float(atom.x),
                SourceFormatArgument::Float(atom.y),
                SourceFormatArgument::Float(atom.z),
                SourceFormatArgument::Bytes(element_pointer.as_const()),
                SourceFormatArgument::Signed(i64::from(isotope)),
                SourceFormatArgument::Signed(i64::from(molfile_charge)),
                SourceFormatArgument::Signed(i64::from(valence)),
            ],
        );
        heap.free(element_pointer)?;
        let _ = result?;
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OrigAtData_WriteToSDfileBondsBlock(
    heap: &mut SourceHeap,
    input: &ORIG_ATOM_DATA,
    mut stream: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    _name: SourceConstPointer<i8>,
    _comment: SourceConstPointer<i8>,
    _label: SourceConstPointer<i8>,
    _value: SourceConstPointer<i8>,
    mut written_bond_ends: Option<&mut INT_ARRAY>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:1122 OrigAtData_WriteToSDfileBondsBlock
    // INCHI✔️❌: int OrigAtData_WriteToSDfileBondsBlock( const ORIG_ATOM_DATA *inp_at_data,
    // INCHI✔️❌:                                         INCHI_IOSTREAM       *fcb,
    // INCHI✔️❌:                                         const char           *name,
    // INCHI✔️❌:                                         const char           *comment,
    // INCHI✔️❌:                                         const char           *szLabel,
    // INCHI✔️❌:                                         const char           *szValue,
    // INCHI✔️❌:                                         INT_ARRAY            *written_bond_ends )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, k, ret = 0;
    // INCHI✔️❌:     int num_atoms = inp_at_data->num_inp_atoms;
    // INCHI✔️❌:     const inp_ATOM *at = inp_at_data->at;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* bonds*/
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (j = 0; j < at[i].valence; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i < at[i].neighbor[j])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 unsigned a1, a2;
    // INCHI✔️❌:                 if ((k = at[i].bond_stereo[j])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* bond stereo */
    // INCHI✔️❌:                     if (k < 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* transposition */
    // INCHI✔️❌:                         a1 = (unsigned)(at[i].neighbor[j] + 1);
    // INCHI✔️❌:                         a2 = (unsigned)(i + 1);
    // INCHI✔️❌:                         inchi_ios_print_nodisplay(fcb, "%3u%3u%3u%3u  0  0  0\n",
    // INCHI✔️❌:                                                   a1, a2, (unsigned)(at[i].bond_type[j]), (unsigned)abs(k));
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* no transposition*/
    // INCHI✔️❌:                         a1 = (unsigned)(i + 1);
    // INCHI✔️❌:                         a2 = (unsigned)(at[i].neighbor[j] + 1);
    // INCHI✔️❌:                         inchi_ios_print_nodisplay(fcb, "%3u%3u%3u%3u  0  0  0\n",
    // INCHI✔️❌:                                                   a1, a2, (unsigned)(at[i].bond_type[j]), (unsigned)abs(k));
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     a1 = (unsigned)(i + 1);
    // INCHI✔️❌:                     a2 = (unsigned)(at[i].neighbor[j] + 1);
    // INCHI✔️❌:                     inchi_ios_print_nodisplay(fcb, "%3u%3u%3u  0  0  0  0\n",
    // INCHI✔️❌:                                               a1, a2, (unsigned)(at[i].bond_type[j]));
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 IntArray_Append(written_bond_ends, a1);
    // INCHI✔️❌:                 IntArray_Append(written_bond_ends, a2);
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OrigAtData_WriteToSDfileBondsBlock
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_WriteToSDfileBondsBlock
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_WriteToSDfileBondsBlock

    let atom_count =
        usize::try_from(input.num_inp_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atoms = if atom_count == 0 {
        Vec::new()
    } else {
        heap.slice(input.at.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };
    for (atom_index, atom) in atoms.iter().enumerate() {
        let mut bond_index = 0_i32;
        while bond_index < i32::from(atom.valence) {
            let order =
                usize::try_from(bond_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let neighbor = i32::from(atom.neighbor[order]);
            if atom_index as i32 >= neighbor {
                bond_index += 1;
                continue;
            }
            let stereo = i32::from(atom.bond_stereo[order]);
            let (first, second) = if stereo < 0 {
                (neighbor.wrapping_add(1) as u32, atom_index as u32 + 1)
            } else {
                (atom_index as u32 + 1, neighbor.wrapping_add(1) as u32)
            };
            let bond_type = u32::from(atom.bond_type[order]);
            if stereo != 0 {
                let format = (*b"%3u%3u%3u%3u  0  0  0\n\0").map(|byte| byte as i8);
                let _ = mol_fmt4_print(
                    heap,
                    stream.as_deref_mut(),
                    stdout,
                    &format,
                    vec![
                        SourceFormatArgument::Unsigned(u64::from(first)),
                        SourceFormatArgument::Unsigned(u64::from(second)),
                        SourceFormatArgument::Unsigned(u64::from(bond_type)),
                        SourceFormatArgument::Unsigned(u64::from(stereo.unsigned_abs())),
                    ],
                )?;
            } else {
                let format = (*b"%3u%3u%3u  0  0  0  0\n\0").map(|byte| byte as i8);
                let _ = mol_fmt4_print(
                    heap,
                    stream.as_deref_mut(),
                    stdout,
                    &format,
                    vec![
                        SourceFormatArgument::Unsigned(u64::from(first)),
                        SourceFormatArgument::Unsigned(u64::from(second)),
                        SourceFormatArgument::Unsigned(u64::from(bond_type)),
                    ],
                )?;
            }
            let _ = IntArray_Append(heap, written_bond_ends.as_deref_mut(), first as i32)?;
            let _ = IntArray_Append(heap, written_bond_ends.as_deref_mut(), second as i32)?;
            bond_index += 1;
        }
    }
    Ok(0)
}

fn mol_fmt4_print_symbol(
    heap: &mut SourceHeap,
    stream: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    format: &[i8],
    first: i32,
    symbol: &[u8],
) -> Result<(), SourceHeapError> {
    let symbol = heap.allocate_model_storage(
        symbol
            .iter()
            .copied()
            .chain(std::iter::once(0))
            .map(|byte| byte as i8)
            .collect(),
    )?;
    let result = mol_fmt4_print(
        heap,
        stream,
        stdout,
        format,
        vec![
            SourceFormatArgument::Signed(i64::from(first)),
            SourceFormatArgument::Bytes(symbol.as_const()),
        ],
    );
    heap.free(symbol)?;
    let _ = result?;
    Ok(())
}

fn mol_fmt4_print_dynamic_bytes(
    heap: &mut SourceHeap,
    stream: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    format: &[i8],
    mut arguments: Vec<SourceFormatArgument>,
    bytes: &[u8],
) -> Result<(), SourceHeapError> {
    let bytes = heap.allocate_model_storage(
        bytes
            .iter()
            .copied()
            .chain(std::iter::once(0))
            .map(|byte| byte as i8)
            .collect(),
    )?;
    arguments.push(SourceFormatArgument::Bytes(bytes.as_const()));
    let result = mol_fmt4_print(heap, stream, stdout, format, arguments);
    heap.free(bytes)?;
    let _ = result?;
    Ok(())
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OrigAtData_WriteToSDfilePolymerData(
    heap: &mut SourceHeap,
    input: &ORIG_ATOM_DATA,
    mut stream: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    _name: SourceConstPointer<i8>,
    _comment: SourceConstPointer<i8>,
    _label: SourceConstPointer<i8>,
    _value: SourceConstPointer<i8>,
    written_bond_ends: Option<&INT_ARRAY>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:1390 OrigAtData_WriteToSDfilePolymerData
    // INCHI✔️❌: int OrigAtData_WriteToSDfilePolymerData( const ORIG_ATOM_DATA *inp_at_data,
    // INCHI✔️❌:                                          INCHI_IOSTREAM       *fcb,
    // INCHI✔️❌:                                          const char           * name,
    // INCHI✔️❌:                                          const char           *comment,
    // INCHI✔️❌:                                          const char           *szLabel,
    // INCHI✔️❌:                                          const char           *szValue,
    // INCHI✔️❌:                                          INT_ARRAY            *written_bond_ends )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int j, k, ju, jj, jprev, ret = 0;
    // INCHI✔️❌:     const char *sty[] = {"NON", "SRU", "MON", "COP", "MOD", "CRO", "MER"};
    // INCHI✔️❌:     const char *sst[] = {"NON", "ALT", "RAN", "BLO"};
    // INCHI✔️❌:     const char *con[] = {"NON", "HT", "HH", "EU"};
    // INCHI✔️❌:     OAD_PolymerUnit *u = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* STY */
    // INCHI✔️❌:     jj = 0;
    // INCHI✔️❌:     jprev = -1;
    // INCHI✔️❌:     for (j = 0; j < inp_at_data->polymer->n; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         u = inp_at_data->polymer->units[j];
    // INCHI✔️❌:         if (u->type > 0 && u->type <= 6)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             jj++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (jj == 8 || j == inp_at_data->polymer->n - 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_print_nodisplay(fcb, "M  STY%3d", jj % 8 ? jj % 8 : 8);
    // INCHI✔️❌:             for (k = jprev + 1; k <= j; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 u = inp_at_data->polymer->units[k];
    // INCHI✔️❌:                 if (u->type > 0 && u->type <= 6)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     inchi_ios_print_nodisplay(fcb, " %3d %3s", u->id, sty[u->type]);
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_ios_print_nodisplay(fcb, "\n");
    // INCHI✔️❌:             jj = 0;
    // INCHI✔️❌:             jprev = j;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* SLB */
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     jprev = -1;
    // INCHI✔️❌:     for (j = 0; j < inp_at_data->polymer->n; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         if (j == 8 || j == inp_at_data->polymer->n - 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             jj = j + 1;
    // INCHI✔️❌:             inchi_ios_print_nodisplay(fcb, "M  SLB%3d", jj % 8 ? jj % 8 : 8);
    // INCHI✔️❌:             for (k = jprev + 1; k < jj; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 u = inp_at_data->polymer->units[k];
    // INCHI✔️❌:                 inchi_ios_print_nodisplay(fcb, " %3d %3d", u->id, u->label);
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_ios_print_nodisplay(fcb, "\n");
    // INCHI✔️❌:             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:             jprev = j;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* SST */
    // INCHI✔️❌:     jj = 0;
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     for (j = 0; j < inp_at_data->polymer->n; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         u = inp_at_data->polymer->units[j];
    // INCHI✔️❌:         if (u->subtype == MOL_FMT_M_SST_ALT || u->subtype == MOL_FMT_M_SST_RAN || u->subtype == MOL_FMT_M_SST_BLK)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             jj++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (jj)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         jj = 0;
    // INCHI✔️❌:         jprev = -1;
    // INCHI✔️❌:         for (j = 0; j < inp_at_data->polymer->n; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             u = inp_at_data->polymer->units[j];
    // INCHI✔️❌:             if (u->subtype == MOL_FMT_M_SST_ALT || u->subtype == MOL_FMT_M_SST_RAN || u->subtype == MOL_FMT_M_SST_BLK)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 jj++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (jj == 8 || j == inp_at_data->polymer->n - 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_print_nodisplay(fcb, "M  SST%3d", jj % 8 ? jj % 8 : 8);
    // INCHI✔️❌:                 for (k = jprev + 1; k <= j; k++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     u = inp_at_data->polymer->units[k];
    // INCHI✔️❌:                     if (u->subtype == MOL_FMT_M_SST_ALT || u->subtype == MOL_FMT_M_SST_RAN || u->subtype == MOL_FMT_M_SST_BLK)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         inchi_ios_print_nodisplay(fcb, " %3d %3s", u->id, sst[u->subtype]);
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 inchi_ios_print_nodisplay(fcb, "\n");
    // INCHI✔️❌:                 jj = 0;
    // INCHI✔️❌:                 jprev = j;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* SCN */
    // INCHI✔️❌:     jj = 0;
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     for (j = 0; j < inp_at_data->polymer->n; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         u = inp_at_data->polymer->units[j];
    // INCHI✔️❌:         if (u->conn == MOL_FMT_M_CONN_HT || u->conn == MOL_FMT_M_CONN_HH || u->conn == MOL_FMT_M_CONN_EU)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             jj++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (jj)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         jj = 0;
    // INCHI✔️❌:         jprev = -1;
    // INCHI✔️❌:         for (j = 0; j < inp_at_data->polymer->n; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             u = inp_at_data->polymer->units[j];
    // INCHI✔️❌:             if (u->conn == MOL_FMT_M_CONN_HT || u->conn == MOL_FMT_M_CONN_HH || u->conn == MOL_FMT_M_CONN_EU)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 jj++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (jj == 8 || j == inp_at_data->polymer->n - 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_print_nodisplay(fcb, "M  SCN%3d", jj % 8 ? jj % 8 : 8);
    // INCHI✔️❌:                 for (k = jprev + 1; k <= j; k++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     u = inp_at_data->polymer->units[k];
    // INCHI✔️❌:                     if (u->conn == MOL_FMT_M_CONN_HT || u->conn == MOL_FMT_M_CONN_HH || u->conn == MOL_FMT_M_CONN_EU)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         inchi_ios_print_nodisplay(fcb, " %3d %3s", u->id, con[u->conn]);
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 inchi_ios_print_nodisplay(fcb, "\n");
    // INCHI✔️❌:                 jj = 0;
    // INCHI✔️❌:                 jprev = j;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* SAL */
    // INCHI✔️❌:     for (ju = 0; ju < inp_at_data->polymer->n; ju++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         u = inp_at_data->polymer->units[ju];
    // INCHI✔️❌:         jj = 0;
    // INCHI✔️❌:         jprev = -1;
    // INCHI✔️❌:         for (j = 0; j < u->na; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             jj++;
    // INCHI✔️❌:             if (jj == 15 || j == u->na - 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_print_nodisplay(fcb, "M  SAL %3d%3d", u->id, jj % 15 ? jj % 15 : 15);
    // INCHI✔️❌:                 for (k = jprev + 1; k <= j; k++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     inchi_ios_print_nodisplay(fcb, " %3d", u->alist[k]);
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 inchi_ios_print_nodisplay(fcb, "\n");
    // INCHI✔️❌:                 jj = 0;
    // INCHI✔️❌:                 jprev = j;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* SBL */
    // INCHI✔️❌:     for (ju = 0; ju < inp_at_data->polymer->n; ju++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         u = inp_at_data->polymer->units[ju];
    // INCHI✔️❌:         jj = 0;
    // INCHI✔️❌:         jprev = -1;
    // INCHI✔️❌:         for (j = 0; j < u->nb; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             jj++;
    // INCHI✔️❌:             if (jj == 15 || j == u->nb - 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_print_nodisplay(fcb, "M  SBL %3d%3d", u->id, jj % 15 ? jj % 15 : 15);
    // INCHI✔️❌:
    // INCHI✔️❌:                 for (k = jprev + 1; k <= j; k++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     int a1, a2, e1, e2, wb, bond_num = 0;
    // INCHI✔️❌:                     a1 = u->blist[2 * k];
    // INCHI✔️❌:                     a2 = u->blist[2 * k + 1];
    // INCHI✔️❌:                     for (wb = 0; wb < written_bond_ends->used / 2; wb++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         e1 = written_bond_ends->item[2 * wb];
    // INCHI✔️❌:                         e2 = written_bond_ends->item[2 * wb + 1];
    // INCHI✔️❌:                         if ((a1 == e1 && a2 == e2) || (a2 == e1 && a1 == e2))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             bond_num = wb + 1;
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (bond_num)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         inchi_ios_print_nodisplay(fcb, " %3d", bond_num);
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 inchi_ios_print_nodisplay(fcb, "\n");
    // INCHI✔️❌:                 jj = 0;
    // INCHI✔️❌:                 jprev = j;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* SDI */
    // INCHI✔️❌:     for (j = 0; j < inp_at_data->polymer->n; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* better than nothing */
    // INCHI✔️❌:         float xmin, xmax, ymin, ymax;
    // INCHI✔️❌:         xmin = ymin = -1.0 * ((float)j + 1.0); /* djb-rwth: cast operator added, 1->1.0 */
    // INCHI✔️❌:         xmax = ymax = +1.0 * ((float)j + 1.0); /* djb-rwth: cast operator added, 1->1.0 */
    // INCHI✔️❌:         u = inp_at_data->polymer->units[j];
    // INCHI✔️❌:         /* u->xbr1[0], x1, y1, x2, y2 u->xbr1[1], u->xbr1[2], u->xbr1[3] */
    // INCHI✔️❌:         inchi_ios_print_nodisplay(fcb, "M  SDI %3d%3d%10.4f%10.4f%10.4f%10.4f\n", u->id, 4, xmin, ymin, xmin, ymax);
    // INCHI✔️❌:         /* u->xbr2[0], u->xbr2[1], u->xbr2[2], u->xbr2[3] */
    // INCHI✔️❌:         inchi_ios_print_nodisplay(fcb, "M  SDI %3d%3d%10.4f%10.4f%10.4f%10.4f\n", u->id, 4, xmax, ymax, xmax, ymin);
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OrigAtData_WriteToSDfilePolymerData
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_WriteToSDfilePolymerData
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_WriteToSDfilePolymerData

    let polymer = heap
        .slice(input.polymer.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    let count = usize::try_from(polymer.n).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let unit_pointers = if count == 0 {
        Vec::new()
    } else {
        heap.slice(polymer.units.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };
    let mut units = Vec::<OAD_PolymerUnit>::with_capacity(count);
    for pointer in unit_pointers {
        units.push(
            heap.slice(pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone(),
        );
    }

    let sty = [
        b"NON".as_slice(),
        b"SRU",
        b"MON",
        b"COP",
        b"MOD",
        b"CRO",
        b"MER",
    ];
    let sst = [b"NON".as_slice(), b"ALT", b"RAN", b"BLO"];
    let con = [b"NON".as_slice(), b"HT", b"HH", b"EU"];
    let signed = |value: i32| SourceFormatArgument::Signed(i64::from(value));

    let mut matched = 0_i32;
    let mut previous = -1_i32;
    for j in 0..polymer.n {
        let unit = &units[j as usize];
        if unit.type_ > 0 && unit.type_ <= 6 {
            matched += 1;
        }
        if matched == 8 || j == polymer.n - 1 {
            let shown = if matched % 8 != 0 { matched % 8 } else { 8 };
            let format = (*b"M  STY%3d\0").map(|byte| byte as i8);
            let _ = mol_fmt4_print(
                heap,
                stream.as_deref_mut(),
                stdout,
                &format,
                vec![signed(shown)],
            )?;
            for k in previous + 1..=j {
                let unit = &units[k as usize];
                if unit.type_ > 0 && unit.type_ <= 6 {
                    let format = (*b" %3d %3s\0").map(|byte| byte as i8);
                    mol_fmt4_print_symbol(
                        heap,
                        stream.as_deref_mut(),
                        stdout,
                        &format,
                        unit.id,
                        sty[unit.type_ as usize],
                    )?;
                }
            }
            let format = [b'\n' as i8, 0];
            let _ = mol_fmt4_print(heap, stream.as_deref_mut(), stdout, &format, vec![])?;
            matched = 0;
            previous = j;
        }
    }

    previous = -1;
    for j in 0..polymer.n {
        if j == 8 || j == polymer.n - 1 {
            let end = j + 1;
            let shown = if end % 8 != 0 { end % 8 } else { 8 };
            let format = (*b"M  SLB%3d\0").map(|byte| byte as i8);
            let _ = mol_fmt4_print(
                heap,
                stream.as_deref_mut(),
                stdout,
                &format,
                vec![signed(shown)],
            )?;
            for k in previous + 1..end {
                let unit = &units[k as usize];
                let format = (*b" %3d %3d\0").map(|byte| byte as i8);
                let _ = mol_fmt4_print(
                    heap,
                    stream.as_deref_mut(),
                    stdout,
                    &format,
                    vec![signed(unit.id), signed(unit.label)],
                )?;
            }
            let format = [b'\n' as i8, 0];
            let _ = mol_fmt4_print(heap, stream.as_deref_mut(), stdout, &format, vec![])?;
            previous = j;
        }
    }

    for (field, names, valid) in [
        (
            0_i32,
            sst.as_slice(),
            [
                MOL_FMT_M_SST_ALT as i32,
                MOL_FMT_M_SST_RAN as i32,
                MOL_FMT_M_SST_BLK as i32,
            ],
        ),
        (
            1_i32,
            con.as_slice(),
            [
                MOL_FMT_M_CONN_HT as i32,
                MOL_FMT_M_CONN_HH as i32,
                MOL_FMT_M_CONN_EU as i32,
            ],
        ),
    ] {
        let field_value =
            |unit: &OAD_PolymerUnit| if field == 0 { unit.subtype } else { unit.conn };
        if !units.iter().any(|unit| valid.contains(&field_value(unit))) {
            continue;
        }
        matched = 0;
        previous = -1;
        for j in 0..polymer.n {
            let value = field_value(&units[j as usize]);
            if valid.contains(&value) {
                matched += 1;
            }
            if matched == 8 || j == polymer.n - 1 {
                let shown = if matched % 8 != 0 { matched % 8 } else { 8 };
                let format = if field == 0 {
                    (*b"M  SST%3d\0").map(|byte| byte as i8)
                } else {
                    (*b"M  SCN%3d\0").map(|byte| byte as i8)
                };
                let _ = mol_fmt4_print(
                    heap,
                    stream.as_deref_mut(),
                    stdout,
                    &format,
                    vec![signed(shown)],
                )?;
                for k in previous + 1..=j {
                    let unit = &units[k as usize];
                    let value = field_value(unit);
                    if valid.contains(&value) {
                        let format = (*b" %3d %3s\0").map(|byte| byte as i8);
                        mol_fmt4_print_symbol(
                            heap,
                            stream.as_deref_mut(),
                            stdout,
                            &format,
                            unit.id,
                            names[value as usize],
                        )?;
                    }
                }
                let format = [b'\n' as i8, 0];
                let _ = mol_fmt4_print(heap, stream.as_deref_mut(), stdout, &format, vec![])?;
                matched = 0;
                previous = j;
            }
        }
    }

    for unit in &units {
        matched = 0;
        previous = -1;
        for j in 0..unit.na {
            matched += 1;
            if matched == 15 || j == unit.na - 1 {
                let shown = if matched % 15 != 0 { matched % 15 } else { 15 };
                let format = (*b"M  SAL %3d%3d\0").map(|byte| byte as i8);
                let _ = mol_fmt4_print(
                    heap,
                    stream.as_deref_mut(),
                    stdout,
                    &format,
                    vec![signed(unit.id), signed(shown)],
                )?;
                let atoms = heap.slice(unit.alist.as_const())?.to_vec();
                for k in previous + 1..=j {
                    let value = *atoms
                        .get(k as usize)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let format = (*b" %3d\0").map(|byte| byte as i8);
                    let _ = mol_fmt4_print(
                        heap,
                        stream.as_deref_mut(),
                        stdout,
                        &format,
                        vec![signed(value)],
                    )?;
                }
                let format = [b'\n' as i8, 0];
                let _ = mol_fmt4_print(heap, stream.as_deref_mut(), stdout, &format, vec![])?;
                matched = 0;
                previous = j;
            }
        }
    }

    let written = if units.iter().any(|unit| unit.nb > 0) {
        let written = written_bond_ends.ok_or(SourceHeapError::NullPointer)?;
        let pair_count = (written.used / 2).max(0) as usize;
        heap.slice(written.item.as_const())?
            .get(..pair_count.saturating_mul(2))
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    } else {
        Vec::new()
    };
    for unit in &units {
        matched = 0;
        previous = -1;
        for j in 0..unit.nb {
            matched += 1;
            if matched == 15 || j == unit.nb - 1 {
                let shown = if matched % 15 != 0 { matched % 15 } else { 15 };
                let format = (*b"M  SBL %3d%3d\0").map(|byte| byte as i8);
                let _ = mol_fmt4_print(
                    heap,
                    stream.as_deref_mut(),
                    stdout,
                    &format,
                    vec![signed(unit.id), signed(shown)],
                )?;
                let bonds = heap.slice(unit.blist.as_const())?.to_vec();
                for k in previous + 1..=j {
                    let offset = (2 * k) as usize;
                    let first = *bonds
                        .get(offset)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let second = *bonds
                        .get(offset + 1)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let mut bond_number = 0_i32;
                    for (index, endpoints) in written.chunks_exact(2).enumerate() {
                        if (first == endpoints[0] && second == endpoints[1])
                            || (second == endpoints[0] && first == endpoints[1])
                        {
                            bond_number = index as i32 + 1;
                            break;
                        }
                    }
                    if bond_number != 0 {
                        let format = (*b" %3d\0").map(|byte| byte as i8);
                        let _ = mol_fmt4_print(
                            heap,
                            stream.as_deref_mut(),
                            stdout,
                            &format,
                            vec![signed(bond_number)],
                        )?;
                    }
                }
                let format = [b'\n' as i8, 0];
                let _ = mol_fmt4_print(heap, stream.as_deref_mut(), stdout, &format, vec![])?;
                matched = 0;
                previous = j;
            }
        }
    }

    for (j, unit) in units.iter().enumerate() {
        let negative = -1.0_f32 * (j as f32 + 1.0_f32);
        let positive = 1.0_f32 * (j as f32 + 1.0_f32);
        let format = (*b"M  SDI %3d%3d%10.4f%10.4f%10.4f%10.4f\n\0").map(|byte| byte as i8);
        let _ = mol_fmt4_print(
            heap,
            stream.as_deref_mut(),
            stdout,
            &format,
            vec![
                signed(unit.id),
                signed(4),
                SourceFormatArgument::Float(f64::from(negative)),
                SourceFormatArgument::Float(f64::from(negative)),
                SourceFormatArgument::Float(f64::from(negative)),
                SourceFormatArgument::Float(f64::from(positive)),
            ],
        )?;
        let _ = mol_fmt4_print(
            heap,
            stream.as_deref_mut(),
            stdout,
            &format,
            vec![
                signed(unit.id),
                signed(4),
                SourceFormatArgument::Float(f64::from(positive)),
                SourceFormatArgument::Float(f64::from(positive)),
                SourceFormatArgument::Float(f64::from(positive)),
                SourceFormatArgument::Float(f64::from(negative)),
            ],
        )?;
    }
    Ok(0)
}

fn mol_fmt4_atom_name(atom: &inp_ATOM) -> Result<String, SourceHeapError> {
    let end = atom
        .elname
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    String::from_utf8(atom.elname[..end].iter().map(|byte| *byte as u8).collect())
        .map_err(|_| SourceHeapError::InvalidSourceTextEncoding)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OrigAtData_WriteToSDfileAdditionalLines(
    heap: &mut SourceHeap,
    input: Option<&ORIG_ATOM_DATA>,
    mut stream: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    name: SourceConstPointer<i8>,
    comment: SourceConstPointer<i8>,
    atoms_dt: i32,
    label: SourceConstPointer<i8>,
    value: SourceConstPointer<i8>,
    num_alias_lines: i32,
    num_charge_lines: i32,
    num_radical_lines: i32,
    num_iso_lines: i32,
    written_bond_ends: Option<&INT_ARRAY>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:1182 OrigAtData_WriteToSDfileAdditionalLines
    // INCHI✔️❌: int OrigAtData_WriteToSDfileAdditionalLines( const ORIG_ATOM_DATA *inp_at_data,
    // INCHI✔️❌:                                              INCHI_IOSTREAM       *fcb,
    // INCHI✔️❌:                                              const char           *name,
    // INCHI✔️❌:                                              const char           *comment,
    // INCHI✔️❌:                                              int                  bAtomsDT,
    // INCHI✔️❌:                                              const char           *szLabel,
    // INCHI✔️❌:                                              const char           *szValue,
    // INCHI✔️❌:                                              int                  nNumAliasLines,
    // INCHI✔️❌:                                              int                  nNumChargeLines,
    // INCHI✔️❌:                                              int                  nNumRadicalLines,
    // INCHI✔️❌:                                              int                  nNumIsoLines,
    // INCHI✔️❌:                                              INT_ARRAY            *written_bond_ends )
    // INCHI✔️❌: {
    // INCHI✔️❌:     char str_m[66], entry[25];
    // INCHI✔️❌:     int i, num_m, k, j, ret = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (inp_at_data) /* djb-rwth: fixing a NULL pointer dereference */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int num_atoms = inp_at_data->num_inp_atoms;
    // INCHI✔️❌:         int is_polymer = inp_at_data && inp_at_data->polymer && inp_at_data->polymer->n > 0 && inp_at_data->valid_polymer;
    // INCHI✔️❌:         const inp_ATOM *at = inp_at_data->at;
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Aliases. 5-3-99 DCh.*/
    // INCHI✔️❌:         if (nNumAliasLines)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num_m = 0;
    // INCHI✔️❌:             for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (ALIASED_AT(i))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     int len;
    // INCHI✔️❌:                     inchi_ios_print_nodisplay(fcb, "A  %d\n", i + 1);
    // INCHI✔️❌:                     num_m++;
    // INCHI✔️❌:                     len = sprintf(str_m, "%s", at[i].elname);
    // INCHI✔️❌:                     /* add isotopic H to the alias */
    // INCHI✔️❌:                     for (k = 0; k < NUM_H_ISOTOPES; k++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         int num_H = at[i].num_iso_H[k] + (k ? 0 : at[i].num_H);
    // INCHI✔️❌:                         if (num_H)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             len += sprintf(str_m + len, "%s", k == 0 ? "H" : k == 1 ? "D" : k == 2 ? "T" : "?");
    // INCHI✔️❌:                             if (num_H != 1)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 len += sprintf(str_m + len, "%d", num_H);
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* Add charge to the Alias */
    // INCHI✔️❌:                     if (at[i].charge)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len += sprintf(str_m + len, "%s", at[i].charge > 0 ? "+" : "-");
    // INCHI✔️❌:                         if (1 < (j = abs(at[i].charge)))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             len += sprintf(str_m + len, "%d", j);
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* Add radical to the Alias */
    // INCHI✔️❌:                     if (at[i].radical == RADICAL_SINGLET)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len += sprintf(str_m + len, "%s", ":");
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else if (at[i].radical == RADICAL_DOUBLET)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len += sprintf(str_m + len, "%s", "^");
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else if (at[i].radical == RADICAL_TRIPLET)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len += sprintf(str_m + len, "%s", "^^");
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     inchi_ios_print_nodisplay(fcb, "%s\n", str_m);
    // INCHI✔️❌:                     num_m++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (num_m != nNumAliasLines)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* error in lines counting*/
    // INCHI✔️❌:                 ret++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* charges*/
    // INCHI✔️❌:         str_m[0] = 0;
    // INCHI✔️❌:         num_m = 0;
    // INCHI✔️❌:         if (nNumChargeLines)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at[i].charge && !ALIASED_AT(i))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     sprintf(entry, " %3d %3d", i + 1, (int)at[i].charge);
    // INCHI✔️❌:                     strcat(str_m, entry);
    // INCHI✔️❌:                     num_m++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if ((i == num_atoms - 1 && num_m) || num_m == 8) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     inchi_ios_print_nodisplay(fcb, "M  CHG%3d%s\n", num_m, str_m);
    // INCHI✔️❌:                     str_m[0] = 0;
    // INCHI✔️❌:                     num_m = 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* radicals*/
    // INCHI✔️❌:         str_m[0] = 0;
    // INCHI✔️❌:         num_m = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nNumRadicalLines)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (at[i].radical && !ALIASED_AT(i))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     int radical = (at[i].radical == RADICAL_SINGLET ||
    // INCHI✔️❌:                                    at[i].radical == RADICAL_DOUBLET ||
    // INCHI✔️❌:                                    at[i].radical == RADICAL_TRIPLET) ? at[i].radical : 0;
    // INCHI✔️❌:                     if (radical)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         sprintf(entry, " %3d %3d", i + 1, radical);
    // INCHI✔️❌:                         strcat(str_m, entry);
    // INCHI✔️❌:                         num_m++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if ((i == num_atoms - 1 && num_m) || num_m == 8) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     inchi_ios_print_nodisplay(fcb, "M  RAD%3d%s\n", num_m, str_m);
    // INCHI✔️❌:                     str_m[0] = 0;
    // INCHI✔️❌:                     num_m = 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* isotopes*/
    // INCHI✔️❌:         str_m[0] = 0;
    // INCHI✔️❌:         num_m = 0;
    // INCHI✔️❌:         if (nNumIsoLines)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int el_num, iso;
    // INCHI✔️❌:             for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /*
    // INCHI✔️❌:                 if ( 0 == strcmp( at[i].elname, "D" ) ) {
    // INCHI✔️❌:                     sprintf( entry, " %3d %3d", i+1, 2 );
    // INCHI✔️❌:                     strcat( str_m, entry );
    // INCHI✔️❌:                     num_m ++;
    // INCHI✔️❌:                 } else
    // INCHI✔️❌:                 if ( 0 == strcmp( at[i].elname, "T" ) ) {
    // INCHI✔️❌:                     sprintf( entry, " %3d %3d", i+1, 3 );
    // INCHI✔️❌:                     strcat( str_m, entry );
    // INCHI✔️❌:                     num_m ++;
    // INCHI✔️❌:                 } else
    // INCHI✔️❌:                 if ( k = at[i].iso_atw_diff ) {
    // INCHI✔️❌:                     int mw = get_atomic_mass_from_elnum( at[i].el_number );
    // INCHI✔️❌:                     mw += (k > 0)? k-1 : k;
    // INCHI✔️❌:                     sprintf( entry, " %3d %3d", i+1, mw );
    // INCHI✔️❌:                     strcat( str_m, entry );
    // INCHI✔️❌:                     num_m ++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 */
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (ANY_ISO(i, bAtomsDT) && !ALIASED_AT(i))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (IS_DEUTERIUM(i))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         iso = 1;
    // INCHI✔️❌:                         el_num = 1;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else if (IS_TRITIUM(i))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         iso = 2;
    // INCHI✔️❌:                         el_num = 1;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         iso = at[i].iso_atw_diff > 0 ? at[i].iso_atw_diff - 1 : at[i].iso_atw_diff;
    // INCHI✔️❌:                         el_num = at[i].el_number;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     iso += get_atomic_mass_from_elnum(el_num);
    // INCHI✔️❌:                     sprintf(entry, " %3d %3d", i + 1, iso);
    // INCHI✔️❌:                     strcat(str_m, entry);
    // INCHI✔️❌:                     num_m++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 if ((i == num_atoms - 1 && num_m) || num_m == 8) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     inchi_ios_print_nodisplay(fcb, "M  ISO%3d%s\n", num_m, str_m);
    // INCHI✔️❌:                     str_m[0] = 0;
    // INCHI✔️❌:                     num_m = 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (is_polymer)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             OrigAtData_WriteToSDfilePolymerData(inp_at_data, fcb, name, comment,
    // INCHI✔️❌:                                                 szLabel, szValue, written_bond_ends);
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         inchi_ios_print_nodisplay(fcb, "M  END\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OrigAtData_WriteToSDfileAdditionalLines
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_WriteToSDfileAdditionalLines
    // INCHI✔️❌: #define ALIASED_AT(i) (0 < NUM_ISO_H(at, i))
    // INCHI✔️❌: #define NUM_ISO_H(AT, N) (AT[N].num_iso_H[0] + AT[N].num_iso_H[1] + AT[N].num_iso_H[2])
    // INCHI✔️❌: #define IS_DEUTERIUM(i) (!strcmp(at[i].elname, "D") || at[i].iso_atw_diff == 2 && !strcmp(at[i].elname, "H"))
    // INCHI✔️❌: #define IS_TRITIUM(i) (!strcmp(at[i].elname, "T") || at[i].iso_atw_diff == 3 && !strcmp(at[i].elname, "H"))
    // INCHI✔️❌: #define ANY_ISO(i, X) ((X) ? (at[i].iso_atw_diff && !IS_DEUTERIUM(i) && !IS_TRITIUM(i)) : (at[i].iso_atw_diff || IS_DEUTERIUM(i) || IS_TRITIUM(i)))
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_WriteToSDfileAdditionalLines

    let Some(input) = input else {
        return Ok(0);
    };
    let atom_count =
        usize::try_from(input.num_inp_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atoms = if atom_count == 0 {
        Vec::new()
    } else {
        heap.slice(input.at.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };
    let is_alias = |atom: &inp_ATOM| {
        atom.num_iso_H
            .iter()
            .map(|value| i32::from(*value))
            .sum::<i32>()
            > 0
    };
    let is_deuterium =
        |atom: &inp_ATOM, name: &str| name == "D" || (atom.iso_atw_diff == 2 && name == "H");
    let is_tritium =
        |atom: &inp_ATOM, name: &str| name == "T" || (atom.iso_atw_diff == 3 && name == "H");
    let signed = |number: i32| SourceFormatArgument::Signed(i64::from(number));
    let mut ret = 0_i32;

    if num_alias_lines != 0 {
        let mut lines = 0_i32;
        for (index, atom) in atoms.iter().enumerate() {
            if !is_alias(atom) {
                continue;
            }
            let format = (*b"A  %d\n\0").map(|byte| byte as i8);
            let _ = mol_fmt4_print(
                heap,
                stream.as_deref_mut(),
                stdout,
                &format,
                vec![signed(index as i32 + 1)],
            )?;
            lines += 1;
            let mut alias = mol_fmt4_atom_name(atom)?;
            for isotope in 0..3 {
                let hydrogen_count = i32::from(atom.num_iso_H[isotope])
                    + if isotope == 0 {
                        i32::from(atom.num_H)
                    } else {
                        0
                    };
                if hydrogen_count != 0 {
                    alias.push_str(["H", "D", "T"][isotope]);
                    if hydrogen_count != 1 {
                        alias.push_str(&hydrogen_count.to_string());
                    }
                }
            }
            if atom.charge != 0 {
                alias.push(if atom.charge > 0 { '+' } else { '-' });
                let magnitude = i32::from(atom.charge).abs();
                if magnitude > 1 {
                    alias.push_str(&magnitude.to_string());
                }
            }
            alias.push_str(match i32::from(atom.radical) {
                value if value == RADICAL_SINGLET as i32 => ":",
                value if value == RADICAL_DOUBLET as i32 => "^",
                value if value == RADICAL_TRIPLET as i32 => "^^",
                _ => "",
            });
            let format = (*b"%s\n\0").map(|byte| byte as i8);
            mol_fmt4_print_dynamic_bytes(
                heap,
                stream.as_deref_mut(),
                stdout,
                &format,
                vec![],
                alias.as_bytes(),
            )?;
            lines += 1;
        }
        if lines != num_alias_lines {
            ret += 1;
        }
    }

    for (enabled, tag, select) in [
        (num_charge_lines, b"M  CHG".as_slice(), 0_i32),
        (num_radical_lines, b"M  RAD".as_slice(), 1_i32),
    ] {
        if enabled == 0 {
            continue;
        }
        let mut entries = Vec::<u8>::new();
        let mut count = 0_i32;
        for (index, atom) in atoms.iter().enumerate() {
            let value = if select == 0 {
                i32::from(atom.charge)
            } else if matches!(
                i32::from(atom.radical),
                value if value == RADICAL_SINGLET as i32
                    || value == RADICAL_DOUBLET as i32
                    || value == RADICAL_TRIPLET as i32
            ) {
                i32::from(atom.radical)
            } else {
                0
            };
            if value != 0 && !is_alias(atom) {
                entries.extend(format!(" {:3} {:3}", index + 1, value).bytes());
                count += 1;
            }
            if ((index + 1 == atoms.len()) && count != 0) || count == 8 {
                let mut format = tag.to_vec();
                format.extend_from_slice(b"%3d%s\n\0");
                let format = format
                    .into_iter()
                    .map(|byte| byte as i8)
                    .collect::<Vec<_>>();
                mol_fmt4_print_dynamic_bytes(
                    heap,
                    stream.as_deref_mut(),
                    stdout,
                    &format,
                    vec![signed(count)],
                    &entries,
                )?;
                entries.clear();
                count = 0;
            }
        }
    }

    if num_iso_lines != 0 {
        let mut entries = Vec::<u8>::new();
        let mut count = 0_i32;
        for (index, atom) in atoms.iter().enumerate() {
            let name = mol_fmt4_atom_name(atom)?;
            let deuterium = is_deuterium(atom, &name);
            let tritium = is_tritium(atom, &name);
            let any_iso = if atoms_dt != 0 {
                atom.iso_atw_diff != 0 && !deuterium && !tritium
            } else {
                atom.iso_atw_diff != 0 || deuterium || tritium
            };
            if any_iso && !is_alias(atom) {
                let (difference, element_number) = if deuterium {
                    (1_i32, 1_i32)
                } else if tritium {
                    (2_i32, 1_i32)
                } else {
                    (
                        if atom.iso_atw_diff > 0 {
                            i32::from(atom.iso_atw_diff) - 1
                        } else {
                            i32::from(atom.iso_atw_diff)
                        },
                        i32::from(atom.el_number),
                    )
                };
                let isotope = difference + get_atomic_mass_from_elnum(element_number)?;
                entries.extend(format!(" {:3} {:3}", index + 1, isotope).bytes());
                count += 1;
            }
            if ((index + 1 == atoms.len()) && count != 0) || count == 8 {
                let format = (*b"M  ISO%3d%s\n\0").map(|byte| byte as i8);
                mol_fmt4_print_dynamic_bytes(
                    heap,
                    stream.as_deref_mut(),
                    stdout,
                    &format,
                    vec![signed(count)],
                    &entries,
                )?;
                entries.clear();
                count = 0;
            }
        }
    }

    let is_polymer = if input.polymer.is_null() {
        false
    } else {
        heap.slice(input.polymer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .n
            > 0
            && input.valid_polymer != 0
    };
    if is_polymer {
        let _ = OrigAtData_WriteToSDfilePolymerData(
            heap,
            input,
            stream.as_deref_mut(),
            stdout,
            name,
            comment,
            label,
            value,
            written_bond_ends,
        )?;
    }
    let format = (*b"M  END\n\0").map(|byte| byte as i8);
    let _ = mol_fmt4_print(heap, stream.as_deref_mut(), stdout, &format, vec![])?;
    Ok(ret)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{INCHI_IOS_STRING, INCHI_IOS_TYPE_STRING, OAD_Polymer, inp_ATOM};

    fn string_stream() -> INCHI_IOSTREAM {
        INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            s: INCHI_IOS_STRING::default(),
            ..INCHI_IOSTREAM::default()
        }
    }

    fn stream_bytes(heap: &SourceHeap, stream: &INCHI_IOSTREAM) -> Vec<u8> {
        if stream.s.pStr.is_null() {
            return Vec::new();
        }
        heap.slice(stream.s.pStr.as_const()).unwrap()[..stream.s.nUsedLength as usize]
            .iter()
            .map(|byte| *byte as u8)
            .collect()
    }

    fn c_text(heap: &mut SourceHeap, bytes: &[u8]) -> SourceConstPointer<i8> {
        heap.allocate_model_storage(
            bytes
                .iter()
                .copied()
                .chain(std::iter::once(0))
                .map(|byte| byte as i8)
                .collect(),
        )
        .unwrap()
        .as_const()
    }

    fn named_atom(name: &[u8]) -> inp_ATOM {
        let mut atom = inp_ATOM::default();
        for (target, source) in atom.elname.iter_mut().zip(name.iter().copied()) {
            *target = source as i8;
        }
        atom
    }

    fn bonded_atom(name: &[u8], element_number: u8, bonds: usize) -> inp_ATOM {
        let mut atom = named_atom(name);
        atom.el_number = element_number;
        atom.valence = bonds as i8;
        atom.chem_bonds_valence = bonds as i8;
        atom.bond_type[..bonds].fill(1);
        atom
    }

    fn expected_atom_line(
        atom: &inp_ATOM,
        symbol: &str,
        isotope: i32,
        charge: i32,
        valence: i32,
    ) -> String {
        format!(
            "{:10.4}{:10.4}{:10.4} {:<3.3}{:2}{:3}  0     0{:3}  0  0  0  0\n",
            atom.x, atom.y, atom.z, symbol, isotope, charge, valence
        )
    }

    #[test]
    fn source_port__mol_fmt4__origatdata_writetosdfileatomsblock__line_989() {
        let mut heap = SourceHeap::default();
        let mut atoms = Vec::new();
        let mut expected = String::new();

        for (source_charge, molfile_charge) in [
            (3, 1),
            (2, 2),
            (1, 3),
            (0, 0),
            (-1, 5),
            (-2, 6),
            (-3, 7),
            (4, 0),
        ] {
            let mut atom = bonded_atom(b"C", 1, 1);
            atom.charge = source_charge;
            atom.x = f64::from(source_charge);
            expected.push_str(&expected_atom_line(&atom, "C", 0, molfile_charge, 0));
            atoms.push(atom);
        }

        for (radical, molfile_charge, valence) in [
            (RADICAL_SINGLET as i8, 0, 0),
            (RADICAL_DOUBLET as i8, 4, 1),
            (RADICAL_TRIPLET as i8, 0, 1),
        ] {
            let mut atom = bonded_atom(b"C", 1, 1);
            atom.radical = radical;
            expected.push_str(&expected_atom_line(&atom, "C", 0, molfile_charge, valence));
            atoms.push(atom);
        }
        let mut charged_radical = bonded_atom(b"C", 1, 1);
        charged_radical.charge = 1;
        charged_radical.radical = RADICAL_DOUBLET as i8;
        expected.push_str(&expected_atom_line(&charged_radical, "C", 0, 0, 0));
        atoms.push(charged_radical);

        let mut alias = bonded_atom(b"O", 6, 3);
        alias.num_iso_H = [1, 0, 0];
        alias.charge = 3;
        alias.radical = RADICAL_DOUBLET as i8;
        alias.iso_atw_diff = 2;
        expected.push_str(&expected_atom_line(&alias, "C", 1, 0, 4));
        atoms.push(alias);

        for (difference, isotope) in [(-3, -3), (-4, 0), (5, 4), (6, 0), (1, 0)] {
            let mut atom = bonded_atom(b"C", 6, 4);
            atom.iso_atw_diff = difference;
            expected.push_str(&expected_atom_line(&atom, "C", isotope, 0, 0));
            atoms.push(atom);
        }

        for name in [b"Zy".as_slice(), b"Zz".as_slice()] {
            let atom = bonded_atom(name, 1, 1);
            expected.push_str(&expected_atom_line(&atom, "*", 0, 0, 0));
            atoms.push(atom);
        }
        let long = bonded_atom(b"ABCDE", 1, 1);
        expected.push_str(&expected_atom_line(&long, "ABC", 0, 0, 0));
        atoms.push(long);

        let atom_pointer = heap.allocate_model_storage(atoms).unwrap();
        let input = ORIG_ATOM_DATA {
            at: atom_pointer,
            num_inp_atoms: heap.slice(atom_pointer.as_const()).unwrap().len() as i32,
            ..ORIG_ATOM_DATA::default()
        };
        let mut stream = string_stream();
        assert_eq!(
            OrigAtData_WriteToSDfileAtomsBlock(
                &mut heap,
                &input,
                Some(&mut stream),
                SourceMutPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
            ),
            Ok(0)
        );
        let actual = String::from_utf8(stream_bytes(&heap, &stream)).unwrap();
        let actual_lines = actual.lines().collect::<Vec<_>>();
        let expected_lines = expected.lines().collect::<Vec<_>>();
        assert_eq!(actual_lines.len(), expected_lines.len());
        for (index, (actual, expected)) in actual_lines.iter().zip(&expected_lines).enumerate() {
            assert_eq!(actual, expected, "atom line {index}");
        }

        let mut deuterium = bonded_atom(b"H", 1, 1);
        deuterium.iso_atw_diff = 2;
        let tritium = bonded_atom(b"T", 1, 1);
        let dt_pointer = heap
            .allocate_model_storage(vec![deuterium.clone(), tritium.clone()])
            .unwrap();
        let dt_input = ORIG_ATOM_DATA {
            at: dt_pointer,
            num_inp_atoms: 2,
            ..ORIG_ATOM_DATA::default()
        };
        for (atoms_dt, symbols, isotopes) in [(0, ["H", "H"], [1, 2]), (1, ["D", "T"], [0, 0])] {
            let mut stream = string_stream();
            assert_eq!(
                OrigAtData_WriteToSDfileAtomsBlock(
                    &mut heap,
                    &dt_input,
                    Some(&mut stream),
                    SourceMutPointer::null(),
                    SourceConstPointer::null(),
                    SourceConstPointer::null(),
                    atoms_dt,
                    SourceConstPointer::null(),
                    SourceConstPointer::null(),
                ),
                Ok(0)
            );
            let expected = expected_atom_line(&deuterium, symbols[0], isotopes[0], 0, 0)
                + &expected_atom_line(&tritium, symbols[1], isotopes[1], 0, 0);
            assert_eq!(stream_bytes(&heap, &stream), expected.as_bytes());
        }

        let empty = ORIG_ATOM_DATA::default();
        let mut empty_stream = string_stream();
        assert_eq!(
            OrigAtData_WriteToSDfileAtomsBlock(
                &mut heap,
                &empty,
                Some(&mut empty_stream),
                SourceMutPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(stream_bytes(&heap, &empty_stream), b"");
    }

    #[test]
    fn source_port__mol_fmt4__origatdata_writetosdfile__line_783() {
        let mut heap = SourceHeap::default();
        let name = c_text(&mut heap, b"Record");
        let comment = c_text(&mut heap, b"Comment");
        let label = c_text(&mut heap, b"TAG");
        let value = c_text(&mut heap, b"value");
        let mut atom = bonded_atom(b"C", 6, 0);
        atom.num_H = 4;
        let atoms = heap.allocate_model_storage(vec![atom.clone()]).unwrap();
        let input = ORIG_ATOM_DATA {
            at: atoms,
            num_inp_atoms: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let mut stream = string_stream();
        assert_eq!(
            OrigAtData_WriteToSDfile(
                &mut heap,
                &input,
                Some(&mut stream),
                SourceMutPointer::null(),
                name,
                comment,
                1,
                0,
                label,
                value,
            ),
            Ok(0)
        );
        let expected = format!(
            "Record\n  InChIV10                                     \nComment\n  1  0  0  0  1  0  0  0  0  0  1 V2000\n{}M  END\n> <TAG>\n value\n\n$$$$\n",
            expected_atom_line(&atom, "C", 0, 0, 0)
        );
        assert_eq!(stream_bytes(&heap, &stream), expected.as_bytes());

        let mut first = bonded_atom(b"H", 1, 1);
        first.neighbor[0] = 1;
        let mut second = bonded_atom(b"H", 1, 1);
        second.neighbor[0] = 0;
        let pair = heap
            .allocate_model_storage(vec![first.clone(), second.clone()])
            .unwrap();
        let bonded = ORIG_ATOM_DATA {
            at: pair,
            num_inp_atoms: 2,
            ..ORIG_ATOM_DATA::default()
        };
        let empty_label = c_text(&mut heap, b"");
        let id_value = c_text(&mut heap, b"42");
        let mut bonded_stream = string_stream();
        assert_eq!(
            OrigAtData_WriteToSDfile(
                &mut heap,
                &bonded,
                Some(&mut bonded_stream),
                SourceMutPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                0,
                0,
                empty_label,
                id_value,
            ),
            Ok(0)
        );
        let expected = format!(
            "\n  InChIV10                                     \n\n  2  1  0  0  0  0  0  0  0  0  1 V2000\n{}{}  1  2  1  0  0  0  0\nM  END\n> <ID>\n 42\n\n$$$$\n",
            expected_atom_line(&first, "H", 0, 0, 0),
            expected_atom_line(&second, "H", 0, 0, 0),
        );
        assert_eq!(stream_bytes(&heap, &bonded_stream), expected.as_bytes());

        let mut saw_int_array_failure = false;
        for allocation_limit in 0..32 {
            let mut failing_heap = SourceHeap::default();
            failing_heap.fail_after_allocations(allocation_limit);
            let mut failing_stream = string_stream();
            let result = OrigAtData_WriteToSDfile(
                &mut failing_heap,
                &ORIG_ATOM_DATA::default(),
                Some(&mut failing_stream),
                SourceMutPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                0,
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
            );
            if result == Ok(_IS_ERROR as i32) {
                saw_int_array_failure = true;
                assert!(!stream_bytes(&failing_heap, &failing_stream).ends_with(b"$$$$\n"));
            }
        }
        assert!(saw_int_array_failure);
    }

    #[test]
    fn source_port__mol_fmt4__origatdata_writetosdfileheaderandcountthings__line_852() {
        let mut heap = SourceHeap::default();
        let name = c_text(&mut heap, &vec![b'N'; 81]);
        let comment = c_text(&mut heap, &vec![b'C'; 80]);
        let mut atoms = Vec::new();

        let mut alias = named_atom(b"C");
        alias.valence = 1;
        alias.num_iso_H = [1, 0, 0];
        alias.charge = 7;
        alias.radical = RADICAL_TRIPLET as i8;
        alias.iso_atw_diff = 6;
        atoms.push(alias);

        for isotope in [1_i8, -4, 6, -3, 5, -1, 2, 3] {
            let mut atom = named_atom(b"C");
            atom.valence = 1;
            atom.charge = 1;
            atom.radical = RADICAL_SINGLET as i8;
            atom.iso_atw_diff = isotope;
            atoms.push(atom);
        }
        let atom_pointer = heap.allocate_model_storage(atoms).unwrap();
        let input = ORIG_ATOM_DATA {
            at: atom_pointer,
            num_inp_atoms: 9,
            ..ORIG_ATOM_DATA::default()
        };
        let mut stream = string_stream();
        let (mut aliases, mut charges, mut radicals, mut isotopes, mut additions, mut bonds) =
            (5, 1, 1, 1, -99, -99);
        assert_eq!(
            OrigAtData_WriteToSDfileHeaderAndCountThings(
                &mut heap,
                &input,
                Some(&mut stream),
                SourceMutPointer::null(),
                name,
                comment,
                -7,
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                &mut aliases,
                &mut charges,
                &mut radicals,
                &mut isotopes,
                &mut additions,
                &mut bonds,
            ),
            Ok(0)
        );
        assert_eq!((aliases, charges, radicals, isotopes), (7, 2, 2, 2));
        assert_eq!((additions, bonds), (14, 4));
        let mut expected = vec![b'N'; 80];
        expected.extend_from_slice(b"\n  InChIV10                                     \n");
        expected.extend(std::iter::repeat_n(b'C', 80));
        expected.extend_from_slice(b"\n  9  4  0  0  1  0  0  0  0  0 14 V2000\n");
        assert_eq!(stream_bytes(&heap, &stream), expected);

        let dt_atoms = heap
            .allocate_model_storage(vec![named_atom(b"D"), named_atom(b"T")])
            .unwrap();
        let dt_input = ORIG_ATOM_DATA {
            at: dt_atoms,
            num_inp_atoms: 2,
            ..ORIG_ATOM_DATA::default()
        };
        for (atoms_dt, expected_iso, expected_add) in [(0, 1, 2), (1, 0, 1)] {
            let mut stream = string_stream();
            let (mut aliases, mut charges, mut radicals, mut isotopes, mut additions, mut bonds) =
                (0, 0, 0, 0, 0, 0);
            assert_eq!(
                OrigAtData_WriteToSDfileHeaderAndCountThings(
                    &mut heap,
                    &dt_input,
                    Some(&mut stream),
                    SourceMutPointer::null(),
                    SourceConstPointer::null(),
                    SourceConstPointer::null(),
                    0,
                    atoms_dt,
                    SourceConstPointer::null(),
                    SourceConstPointer::null(),
                    &mut aliases,
                    &mut charges,
                    &mut radicals,
                    &mut isotopes,
                    &mut additions,
                    &mut bonds,
                ),
                Ok(0)
            );
            assert_eq!(
                (aliases, charges, radicals, isotopes, additions, bonds),
                (0, 0, 0, expected_iso, expected_add, 0)
            );
            assert_eq!(
                stream_bytes(&heap, &stream),
                format!(
                    "\n  InChIV10                                     \n\n  2  0  0  0  0  0  0  0  0  0{:3} V2000\n",
                    expected_add
                )
                .into_bytes()
            );
        }
    }

    #[test]
    fn source_port__mol_fmt4__intarray_alloc__line_510() {
        let mut heap = SourceHeap::default();
        let stale = heap.allocate_model_storage(vec![9_i32]).unwrap();
        let mut items = INT_ARRAY {
            item: stale,
            allocated: 11,
            used: 7,
            increment: 5,
        };
        assert_eq!(IntArray_Alloc(&mut heap, &mut items, 3), Ok(0));
        assert_eq!((items.allocated, items.used, items.increment), (3, 0, 3));
        assert_eq!(heap.slice(items.item.as_const()).unwrap(), &[0, 0, 0]);
        assert_eq!(heap.slice(stale.as_const()).unwrap(), &[9]);

        let mut zero = INT_ARRAY {
            allocated: -1,
            used: -2,
            increment: -3,
            ..INT_ARRAY::default()
        };
        assert_eq!(IntArray_Alloc(&mut heap, &mut zero, 0), Ok(0));
        assert!(!zero.item.is_null());
        assert_eq!((zero.allocated, zero.used, zero.increment), (0, 0, 0));
        assert_eq!(heap.slice(zero.item.as_const()).unwrap(), &[] as &[i32]);

        let mut failure_heap = SourceHeap::default();
        let old = failure_heap.allocate_model_storage(vec![1_i32, 2]).unwrap();
        let mut failure = INT_ARRAY {
            item: old,
            allocated: 8,
            used: 6,
            increment: 4,
        };
        failure_heap.fail_after_allocations(0);
        assert_eq!(IntArray_Alloc(&mut failure_heap, &mut failure, 2), Ok(-1));
        assert!(failure.item.is_null());
        assert_eq!(
            (failure.allocated, failure.used, failure.increment),
            (8, 6, 4)
        );
        assert_eq!(failure_heap.slice(old.as_const()).unwrap(), &[1, 2]);

        let old = heap.allocate_model_storage(vec![3_i32]).unwrap();
        let mut negative = INT_ARRAY {
            item: old,
            allocated: 1,
            used: 1,
            increment: 1,
        };
        assert_eq!(IntArray_Alloc(&mut heap, &mut negative, -1), Ok(-1));
        assert!(negative.item.is_null());
        assert_eq!(
            (negative.allocated, negative.used, negative.increment),
            (1, 1, 1)
        );
        assert_eq!(heap.slice(old.as_const()).unwrap(), &[3]);
    }

    #[test]
    fn source_port__mol_fmt4__intarray_realloc__line_525() {
        let mut heap = SourceHeap::default();
        assert_eq!(IntArray_ReAlloc(&mut heap, None), Ok(-1));

        for mut invalid in [
            INT_ARRAY::default(),
            INT_ARRAY {
                item: heap.allocate_model_storage(vec![7_i32]).unwrap(),
                allocated: 0,
                increment: 1,
                used: 1,
            },
            INT_ARRAY {
                item: heap.allocate_model_storage(vec![8_i32]).unwrap(),
                allocated: 1,
                increment: 0,
                used: 1,
            },
        ] {
            let before = invalid.clone();
            assert_eq!(IntArray_ReAlloc(&mut heap, Some(&mut invalid)), Ok(-1));
            assert_eq!(invalid, before);
        }

        let old = heap
            .allocate_model_storage(vec![11_i32, -22, 33, 44])
            .unwrap();
        let mut items = INT_ARRAY {
            item: old,
            allocated: 4,
            increment: 3,
            used: 3,
        };
        assert_eq!(IntArray_ReAlloc(&mut heap, Some(&mut items)), Ok(0));
        assert_eq!((items.allocated, items.increment, items.used), (7, 3, 3));
        assert_ne!(items.item, old);
        assert_eq!(
            heap.slice(items.item.as_const()).unwrap(),
            &[11, -22, 33, 0, 0, 0, 0]
        );
        assert_eq!(
            heap.slice(old.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let mut failure_heap = SourceHeap::default();
        let failure_old = failure_heap.allocate_model_storage(vec![5_i32, 6]).unwrap();
        let mut failure = INT_ARRAY {
            item: failure_old,
            allocated: 2,
            increment: 2,
            used: 2,
        };
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            IntArray_ReAlloc(&mut failure_heap, Some(&mut failure)),
            Ok(-1)
        );
        assert!(failure.item.is_null());
        assert_eq!(
            (failure.allocated, failure.increment, failure.used),
            (2, 2, 2)
        );
        assert_eq!(
            failure_heap.slice(failure_old.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__mol_fmt4__intarray_append__line_550() {
        let mut heap = SourceHeap::default();
        assert_eq!(IntArray_Append(&mut heap, None, 7), Ok(-1));

        let storage = heap
            .allocate_model_storage(vec![10_i32, 20, 91, 92])
            .unwrap();
        let mut items = INT_ARRAY {
            item: storage,
            allocated: 4,
            used: 2,
            increment: 2,
        };
        assert_eq!(IntArray_Append(&mut heap, Some(&mut items), -30), Ok(0));
        assert_eq!((items.allocated, items.used), (4, 3));
        assert_eq!(
            heap.slice(items.item.as_const()).unwrap(),
            &[10, 20, -30, 92]
        );

        items.used = 4;
        assert_eq!(IntArray_Append(&mut heap, Some(&mut items), 50), Ok(0));
        assert_eq!((items.allocated, items.used), (6, 5));
        assert_eq!(
            heap.slice(items.item.as_const()).unwrap(),
            &[10, 20, -30, 92, 50, 0]
        );

        let mut failure_heap = SourceHeap::default();
        let old = failure_heap.allocate_model_storage(vec![1_i32, 2]).unwrap();
        let mut failure = INT_ARRAY {
            item: old,
            allocated: 2,
            used: 2,
            increment: 2,
        };
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            IntArray_Append(&mut failure_heap, Some(&mut failure), 3),
            Ok(-1)
        );
        assert!(failure.item.is_null());
        assert_eq!(
            (failure.allocated, failure.used, failure.increment),
            (2, 2, 2)
        );
        assert_eq!(
            failure_heap.slice(old.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let boundary_storage = heap.allocate_model_storage(vec![0_i32]).unwrap();
        let mut boundary = INT_ARRAY {
            item: boundary_storage,
            allocated: i32::MAX,
            used: i32::MAX,
            increment: 1,
        };
        assert_eq!(
            IntArray_Append(&mut heap, Some(&mut boundary), 1),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(boundary.used, i32::MAX);
    }

    #[test]
    fn source_port__mol_fmt4__intarray_appendifabsent__line_572() {
        let mut heap = SourceHeap::default();
        let storage = heap
            .allocate_model_storage(vec![i32::MIN, 7, i32::MAX, 91, 92])
            .unwrap();
        let mut items = INT_ARRAY {
            item: storage,
            allocated: 5,
            used: 3,
            increment: 2,
        };

        for duplicate in [i32::MIN, 7, i32::MAX] {
            let before = heap.slice(items.item.as_const()).unwrap().to_vec();
            let pointer = items.item;
            assert_eq!(
                IntArray_AppendIfAbsent(&mut heap, &mut items, duplicate),
                Ok(0)
            );
            assert_eq!(items.item, pointer);
            assert_eq!((items.allocated, items.used, items.increment), (5, 3, 2));
            assert_eq!(heap.slice(items.item.as_const()).unwrap(), before);
        }

        assert_eq!(IntArray_AppendIfAbsent(&mut heap, &mut items, -11), Ok(0));
        assert_eq!((items.allocated, items.used), (5, 4));
        assert_eq!(
            heap.slice(items.item.as_const()).unwrap(),
            &[i32::MIN, 7, i32::MAX, -11, 92]
        );

        items.used = 5;
        let old = items.item;
        assert_eq!(IntArray_AppendIfAbsent(&mut heap, &mut items, 123), Ok(0));
        assert_ne!(items.item, old);
        assert_eq!((items.allocated, items.used, items.increment), (7, 6, 2));
        assert_eq!(
            heap.slice(items.item.as_const()).unwrap(),
            &[i32::MIN, 7, i32::MAX, -11, 92, 123, 0]
        );
        assert_eq!(
            heap.slice(old.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let empty_storage = heap.allocate_model_storage(vec![44_i32]).unwrap();
        let mut empty = INT_ARRAY {
            item: empty_storage,
            allocated: 1,
            used: 0,
            increment: 1,
        };
        assert_eq!(
            IntArray_AppendIfAbsent(&mut heap, &mut empty, i32::MIN),
            Ok(0)
        );
        assert_eq!(empty.used, 1);
        assert_eq!(heap.slice(empty.item.as_const()).unwrap(), &[i32::MIN]);

        let mut failure_heap = SourceHeap::default();
        let failure_old = failure_heap.allocate_model_storage(vec![1_i32, 2]).unwrap();
        let mut failure = INT_ARRAY {
            item: failure_old,
            allocated: 2,
            used: 2,
            increment: 2,
        };
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            IntArray_AppendIfAbsent(&mut failure_heap, &mut failure, 3),
            Ok(-1)
        );
        assert!(failure.item.is_null());
        assert_eq!(
            (failure.allocated, failure.used, failure.increment),
            (2, 2, 2)
        );
        assert_eq!(
            failure_heap.slice(failure_old.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let invalid_storage = heap.allocate_model_storage(vec![5_i32]).unwrap();
        let mut negative = INT_ARRAY {
            item: invalid_storage,
            allocated: 1,
            used: -1,
            increment: 1,
        };
        assert_eq!(
            IntArray_AppendIfAbsent(&mut heap, &mut negative, 5),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(negative.used, -1);

        let mut out_of_bounds = INT_ARRAY {
            item: invalid_storage,
            allocated: 2,
            used: 2,
            increment: 1,
        };
        assert_eq!(
            IntArray_AppendIfAbsent(&mut heap, &mut out_of_bounds, 8),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!((out_of_bounds.allocated, out_of_bounds.used), (2, 2));
    }

    #[test]
    fn source_port__mol_fmt4__intarray_debugprint__line_582() {
        IntArray_DebugPrint(SourceMutPointer::null());

        let mut heap = SourceHeap::default();
        for used in [i32::MIN, -1, 0, 1, 2, i32::MAX] {
            let items = heap
                .allocate_model_storage(vec![INT_ARRAY {
                    item: SourceMutPointer::null(),
                    allocated: i32::MIN,
                    used,
                    increment: i32::MAX,
                }])
                .unwrap();
            IntArray_DebugPrint(items);
            let items = &heap.slice(items.as_const()).unwrap()[0];
            assert!(items.item.is_null());
            assert_eq!(items.used, used);
        }

        let item = heap
            .allocate_model_storage(vec![i32::MIN, 0, i32::MAX])
            .unwrap();
        let items = heap
            .allocate_model_storage(vec![INT_ARRAY {
                item,
                allocated: 3,
                used: 3,
                increment: 3,
            }])
            .unwrap();
        IntArray_DebugPrint(items);
        assert_eq!(
            heap.slice(item.as_const()).unwrap(),
            &[i32::MIN, 0, i32::MAX]
        );
    }

    #[test]
    fn source_port__mol_fmt4__intarray_reset__line_603() {
        let mut heap = SourceHeap::default();
        let pointer = heap
            .allocate_model_storage(vec![11_i32, -7, i32::MAX, i32::MIN])
            .unwrap();
        let mut items = INT_ARRAY {
            item: pointer,
            allocated: 4,
            used: 3,
            increment: -9,
        };
        IntArray_Reset(&mut items);
        assert_eq!(items.item, pointer);
        assert_eq!((items.allocated, items.used, items.increment), (4, 0, -9));
        assert_eq!(
            heap.slice(pointer.as_const()).unwrap(),
            &[11, -7, i32::MAX, i32::MIN]
        );

        items.used = i32::MIN;
        IntArray_Reset(&mut items);
        assert_eq!(items.used, 0);
        assert_eq!(items.item, pointer);

        let mut null_storage = INT_ARRAY {
            item: SourceMutPointer::null(),
            allocated: i32::MIN,
            used: i32::MAX,
            increment: i32::MAX,
        };
        IntArray_Reset(&mut null_storage);
        assert_eq!(
            null_storage,
            INT_ARRAY {
                item: SourceMutPointer::null(),
                allocated: i32::MIN,
                used: 0,
                increment: i32::MAX,
            }
        );
    }

    #[test]
    fn source_port__mol_fmt4__intarray_free__line_610() {
        let mut heap = SourceHeap::default();
        assert_eq!(IntArray_Free(&mut heap, None), Ok(()));

        let mut no_item = INT_ARRAY {
            item: SourceMutPointer::null(),
            allocated: 9,
            used: -4,
            increment: 3,
        };
        assert_eq!(IntArray_Free(&mut heap, Some(&mut no_item)), Ok(()));
        assert_eq!(no_item, INT_ARRAY::default());

        let storage = heap.allocate_model_storage(vec![11_i32, 22, 33]).unwrap();
        let mut items = INT_ARRAY {
            item: storage,
            allocated: 3,
            used: 2,
            increment: 3,
        };
        assert_eq!(IntArray_Free(&mut heap, Some(&mut items)), Ok(()));
        assert_eq!(items, INT_ARRAY::default());
        assert_eq!(
            heap.slice(storage.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__mol_fmt4__origatdata_writetosdfilebondsblock__line_1122() {
        let mut heap = SourceHeap::default();
        let mut atoms = vec![inp_ATOM::default(); 3];
        atoms[0].valence = 2;
        atoms[0].neighbor[..2].copy_from_slice(&[1, 2]);
        atoms[0].bond_type[..2].copy_from_slice(&[1, 2]);
        atoms[0].bond_stereo[..2].copy_from_slice(&[0, -6]);
        atoms[1].valence = 2;
        atoms[1].neighbor[..2].copy_from_slice(&[0, 2]);
        atoms[1].bond_type[..2].copy_from_slice(&[1, 3]);
        atoms[1].bond_stereo[..2].copy_from_slice(&[0, 1]);
        atoms[2].valence = 2;
        atoms[2].neighbor[..2].copy_from_slice(&[0, 1]);
        atoms[2].bond_type[..2].copy_from_slice(&[2, 3]);
        atoms[2].bond_stereo[..2].copy_from_slice(&[6, 1]);
        let input = ORIG_ATOM_DATA {
            at: heap.allocate_model_storage(atoms).unwrap(),
            num_inp_atoms: 3,
            num_inp_bonds: 3,
            ..ORIG_ATOM_DATA::default()
        };
        let written_storage = heap.allocate_model_storage(vec![0_i32; 6]).unwrap();
        let mut written = INT_ARRAY {
            item: written_storage,
            allocated: 6,
            used: 0,
            increment: 6,
        };
        let mut stream = string_stream();
        assert_eq!(
            OrigAtData_WriteToSDfileBondsBlock(
                &mut heap,
                &input,
                Some(&mut stream),
                SourceMutPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                Some(&mut written),
            ),
            Ok(0)
        );
        assert_eq!(
            stream_bytes(&heap, &stream),
            b"  1  2  1  0  0  0  0\n  3  1  2  6  0  0  0\n  2  3  3  1  0  0  0\n"
        );
        assert_eq!(written.used, 6);
        assert_eq!(
            heap.slice(written.item.as_const()).unwrap(),
            &[1, 2, 3, 1, 2, 3]
        );

        let null_written_storage = heap.allocate_model_storage(vec![0_i32; 6]).unwrap();
        let mut null_written = INT_ARRAY {
            item: null_written_storage,
            allocated: 6,
            used: 0,
            increment: 6,
        };
        assert_eq!(
            OrigAtData_WriteToSDfileBondsBlock(
                &mut heap,
                &input,
                None,
                SourceMutPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                Some(&mut null_written),
            ),
            Ok(0)
        );
        assert_eq!(null_written.used, 6);
        assert_eq!(
            heap.slice(null_written.item.as_const()).unwrap(),
            &[1, 2, 3, 1, 2, 3]
        );
    }

    #[test]
    fn source_port__mol_fmt4__origatdata_writetosdfileadditionallines__line_1182() {
        fn set_name(atom: &mut inp_ATOM, name: &str) {
            atom.elname.fill(0);
            for (target, source) in atom.elname.iter_mut().zip(name.bytes()) {
                *target = source as i8;
            }
        }

        let mut heap = SourceHeap::default();
        let mut null_stream = string_stream();
        assert_eq!(
            OrigAtData_WriteToSDfileAdditionalLines(
                &mut heap,
                None,
                Some(&mut null_stream),
                SourceMutPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                1,
                1,
                1,
                1,
                None,
            ),
            Ok(0)
        );
        assert_eq!(stream_bytes(&heap, &null_stream), b"");

        let mut atoms = vec![inp_ATOM::default(); 12];
        set_name(&mut atoms[0], "C");
        atoms[0].num_H = 1;
        atoms[0].num_iso_H = [1, 0, 0];
        atoms[0].radical = RADICAL_SINGLET as i8;
        set_name(&mut atoms[1], "N");
        atoms[1].num_iso_H = [0, 1, 0];
        atoms[1].radical = RADICAL_DOUBLET as i8;
        set_name(&mut atoms[2], "O");
        atoms[2].num_iso_H = [0, 0, 1];
        atoms[2].charge = -2;
        atoms[2].radical = RADICAL_TRIPLET as i8;

        for (offset, atom) in atoms[3..].iter_mut().enumerate() {
            set_name(atom, "C");
            atom.el_number = 6;
            atom.charge = if offset == 8 { -4 } else { 1 };
            atom.radical = if offset < 8 {
                (offset % 3 + 1) as i8
            } else {
                4
            };
        }
        set_name(&mut atoms[3], "D");
        atoms[3].el_number = 1;
        set_name(&mut atoms[4], "T");
        atoms[4].el_number = 1;
        set_name(&mut atoms[5], "H");
        atoms[5].el_number = 1;
        atoms[5].iso_atw_diff = 2;
        set_name(&mut atoms[6], "H");
        atoms[6].el_number = 1;
        atoms[6].iso_atw_diff = 3;
        atoms[7].iso_atw_diff = 2;
        atoms[8].iso_atw_diff = -1;

        let atom_pointer = heap.allocate_model_storage(atoms).unwrap();
        let input = ORIG_ATOM_DATA {
            at: atom_pointer,
            num_inp_atoms: 12,
            ..ORIG_ATOM_DATA::default()
        };
        let mut stream = string_stream();
        assert_eq!(
            OrigAtData_WriteToSDfileAdditionalLines(
                &mut heap,
                Some(&input),
                Some(&mut stream),
                SourceMutPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                6,
                1,
                1,
                1,
                None,
            ),
            Ok(0)
        );
        assert_eq!(
            stream_bytes(&heap, &stream),
            concat!(
                "A  1\n",
                "CH2:\n",
                "A  2\n",
                "ND^\n",
                "A  3\n",
                "OT-2^^\n",
                "M  CHG  8   4   1   5   1   6   1   7   1   8   1   9   1  10   1  11   1\n",
                "M  CHG  1  12  -4\n",
                "M  RAD  8   4   1   5   2   6   3   7   1   8   2   9   3  10   1  11   2\n",
                "M  ISO  6   4   2   5   3   6   2   7   3   8  13   9  11\n",
                "M  END\n",
            )
            .as_bytes()
        );

        let mut mismatch_stream = string_stream();
        assert_eq!(
            OrigAtData_WriteToSDfileAdditionalLines(
                &mut heap,
                Some(&input),
                Some(&mut mismatch_stream),
                SourceMutPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                1,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                5,
                0,
                0,
                1,
                None,
            ),
            Ok(1)
        );
        assert_eq!(
            stream_bytes(&heap, &mismatch_stream),
            concat!(
                "A  1\n",
                "CH2:\n",
                "A  2\n",
                "ND^\n",
                "A  3\n",
                "OT-2^^\n",
                "M  ISO  2   8  13   9  11\n",
                "M  END\n",
            )
            .as_bytes()
        );

        let polymer_unit = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                id: 2,
                type_: 1,
                label: 3,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let polymer_units = heap.allocate_model_storage(vec![polymer_unit]).unwrap();
        let polymer = heap
            .allocate_model_storage(vec![OAD_Polymer {
                units: polymer_units,
                n: 1,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        let polymer_input = ORIG_ATOM_DATA {
            polymer,
            valid_polymer: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let mut polymer_stream = string_stream();
        assert_eq!(
            OrigAtData_WriteToSDfileAdditionalLines(
                &mut heap,
                Some(&polymer_input),
                Some(&mut polymer_stream),
                SourceMutPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                0,
                0,
                0,
                0,
                None,
            ),
            Ok(0)
        );
        assert_eq!(
            stream_bytes(&heap, &polymer_stream),
            concat!(
                "M  STY  1   2 SRU\n",
                "M  SLB  1   2   3\n",
                "M  SDI   2  4   -1.0000   -1.0000   -1.0000    1.0000\n",
                "M  SDI   2  4    1.0000    1.0000    1.0000   -1.0000\n",
                "M  END\n",
            )
            .as_bytes()
        );
    }

    #[test]
    fn source_port__mol_fmt4__origatdata_writetosdfilepolymerdata__line_1390() {
        let mut heap = SourceHeap::default();
        let empty_polymer = heap
            .allocate_model_storage(vec![OAD_Polymer::default()])
            .unwrap();
        let empty_input = ORIG_ATOM_DATA {
            polymer: empty_polymer,
            ..ORIG_ATOM_DATA::default()
        };
        let mut empty_stream = string_stream();
        assert_eq!(
            OrigAtData_WriteToSDfilePolymerData(
                &mut heap,
                &empty_input,
                Some(&mut empty_stream),
                SourceMutPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                None,
            ),
            Ok(0)
        );
        assert_eq!(stream_bytes(&heap, &empty_stream), b"");

        let atom_list = heap
            .allocate_model_storage((1_i32..=16).collect::<Vec<_>>())
            .unwrap();
        let mut bond_values = Vec::new();
        let mut written_values = Vec::new();
        for index in 1_i32..=16 {
            bond_values.extend([index, index + 1]);
            if index <= 15 {
                written_values.extend([index + 1, index]);
            }
        }
        let bond_list = heap.allocate_model_storage(bond_values).unwrap();
        let written_storage = heap.allocate_model_storage(written_values.clone()).unwrap();
        let written = INT_ARRAY {
            item: written_storage,
            allocated: written_values.len() as i32,
            used: written_values.len() as i32,
            increment: 0,
        };
        let unit = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                id: 7,
                type_: 6,
                subtype: MOL_FMT_M_SST_BLK as i32,
                conn: MOL_FMT_M_CONN_EU as i32,
                label: 9,
                na: 16,
                nb: 16,
                alist: atom_list,
                blist: bond_list,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let units = heap.allocate_model_storage(vec![unit]).unwrap();
        let polymer = heap
            .allocate_model_storage(vec![OAD_Polymer {
                units,
                n: 1,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        let input = ORIG_ATOM_DATA {
            polymer,
            ..ORIG_ATOM_DATA::default()
        };
        let mut stream = string_stream();
        assert_eq!(
            OrigAtData_WriteToSDfilePolymerData(
                &mut heap,
                &input,
                Some(&mut stream),
                SourceMutPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                Some(&written),
            ),
            Ok(0)
        );
        assert_eq!(
            stream_bytes(&heap, &stream),
            concat!(
                "M  STY  1   7 MER\n",
                "M  SLB  1   7   9\n",
                "M  SST  1   7 BLO\n",
                "M  SCN  1   7  EU\n",
                "M  SAL   7 15   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15\n",
                "M  SAL   7  1  16\n",
                "M  SBL   7 15   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15\n",
                "M  SBL   7  1\n",
                "M  SDI   7  4   -1.0000   -1.0000   -1.0000    1.0000\n",
                "M  SDI   7  4    1.0000    1.0000    1.0000   -1.0000\n",
            )
            .as_bytes()
        );

        let mut chunk_units = Vec::new();
        for index in 0_i32..10 {
            chunk_units.push(
                heap.allocate_model_storage(vec![OAD_PolymerUnit {
                    id: index + 1,
                    type_: if index < 8 { index % 6 + 1 } else { 0 },
                    label: 20 + index,
                    ..OAD_PolymerUnit::default()
                }])
                .unwrap(),
            );
        }
        let chunk_units = heap.allocate_model_storage(chunk_units).unwrap();
        let chunk_polymer = heap
            .allocate_model_storage(vec![OAD_Polymer {
                units: chunk_units,
                n: 10,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        let chunk_input = ORIG_ATOM_DATA {
            polymer: chunk_polymer,
            ..ORIG_ATOM_DATA::default()
        };
        let mut chunk_stream = string_stream();
        assert_eq!(
            OrigAtData_WriteToSDfilePolymerData(
                &mut heap,
                &chunk_input,
                Some(&mut chunk_stream),
                SourceMutPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                None,
            ),
            Ok(0)
        );
        let bytes = stream_bytes(&heap, &chunk_stream);
        let text = std::str::from_utf8(&bytes).unwrap();
        let lines = text.lines().collect::<Vec<_>>();
        assert_eq!(lines.len(), 24);
        assert_eq!(
            &lines[..4],
            &[
                "M  STY  8   1 SRU   2 MON   3 COP   4 MOD   5 CRO   6 MER   7 SRU   8 MON",
                "M  STY  8",
                "M  SLB  1   1  20   2  21   3  22   4  23   5  24   6  25   7  26   8  27   9  28",
                "M  SLB  2  10  29",
            ]
        );

        assert_eq!(
            OrigAtData_WriteToSDfilePolymerData(
                &mut heap,
                &chunk_input,
                None,
                SourceMutPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                None,
            ),
            Ok(0)
        );
    }
}
