use crate::source_types::{
    AT_NUMB, BOND_TYPE_ALTERN, BOND_TYPE_MASK, BOND_TYPE_TRIPLE, NUM_H_ISOTOPES, RADICAL_DOUBLET,
    RADICAL_SINGLET, RADICAL_TRIPLET, S_CHAR, SourceConstPointer, SourceHeap, SourceHeapError,
    SourceMutPointer, inp_ATOM,
    local_util::{ERR_ELEM, MAX_ATOM_CHARGE, MAX_NUM_VALENCES, MIN_ATOM_CHARGE, NEUTRAL_STATE},
};

pub(crate) fn inchi_malloc(
    heap: &mut SourceHeap,
    byte_count: u64,
) -> Result<SourceMutPointer<i8>, SourceHeapError> {
    // BEGIN INCHI ACTIVE MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mode.h:1171 inchi_malloc
    // INCHI✔️❌: #ifndef inchi_malloc
    // INCHI✔️❌: #define inchi_malloc   malloc
    // INCHI✔️❌: #endif
    // END INCHI ACTIVE MACRO: inchi_malloc

    let byte_count = usize::try_from(byte_count)
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    let mut bytes = Vec::new();
    bytes
        .try_reserve_exact(byte_count)
        .map_err(|_| SourceHeapError::AllocationFailed)?;
    bytes.resize(byte_count, 0_i8);
    heap.allocate(bytes)
}

pub(crate) fn inchi_calloc<T: Default + 'static>(
    heap: &mut SourceHeap,
    count: u64,
    source_element_size: u64,
) -> Result<SourceMutPointer<T>, SourceHeapError> {
    // BEGIN INCHI ACTIVE MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mode.h:1174 inchi_calloc
    // INCHI✔️❌: #ifndef inchi_calloc
    // INCHI✔️❌: #define inchi_calloc   calloc
    // INCHI✔️❌: #endif
    // END INCHI ACTIVE MACRO: inchi_calloc

    count
        .checked_mul(source_element_size)
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    let count =
        usize::try_from(count).map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| SourceHeapError::AllocationFailed)?;
    values.resize_with(count, T::default);
    heap.allocate(values)
}

pub(crate) fn inchi_realloc<T: Clone + Default + 'static>(
    heap: &mut SourceHeap,
    pointer: SourceMutPointer<T>,
    count: u64,
) -> Result<SourceMutPointer<T>, SourceHeapError> {
    // BEGIN INCHI ACTIVE MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mode.h:1177 inchi_realloc
    // INCHI✔️❌: #ifndef inchi_realloc
    // INCHI✔️❌: #define inchi_realloc   realloc
    // INCHI✔️❌: #endif
    // END INCHI ACTIVE MACRO: inchi_realloc

    let count =
        usize::try_from(count).map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    let old = if pointer.is_null() {
        Vec::new()
    } else {
        heap.slice(pointer.as_const())?.to_vec()
    };
    let mut replacement = Vec::new();
    replacement
        .try_reserve_exact(count)
        .map_err(|_| SourceHeapError::AllocationFailed)?;
    replacement.resize_with(count, T::default);
    let copied = old.len().min(count);
    replacement[..copied].clone_from_slice(&old[..copied]);
    let replacement = match heap.allocate(replacement) {
        Ok(replacement) => replacement,
        Err(SourceHeapError::AllocationFailed) => return Ok(SourceMutPointer::null()),
        Err(error) => return Err(error),
    };
    inchi_free(heap, pointer)?;
    Ok(replacement)
}

pub(crate) fn inchi_free<T: 'static>(
    heap: &mut SourceHeap,
    pointer: SourceMutPointer<T>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI ACTIVE MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mode.h:1180 inchi_free
    // INCHI✔️❌: #ifndef inchi_free
    // INCHI✔️❌: #define inchi_free(X)  do{ if(X) free(X); }while(0)
    // INCHI✔️❌: #endif
    // END INCHI ACTIVE MACRO: inchi_free

    heap.free(pointer)
}

#[rustfmt::skip]
pub(crate) fn normalize_string(
    heap: &mut SourceHeap,
    name: SourceMutPointer<i8>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1589 normalize_string
    // INCHI✔️❌: complete source frame follows verbatim; SourceHeap checked access adds allocation-map overhead.
    /*
int normalize_string( char* name )
{
    int i, len, n;

    len = (int) strlen( name );

    for (i = 0, n = 0; i < len; i++)
    {
        if (isspace( UCINT name[i] ) /*|| !isprint( UCINT name[i] )*/)
        {
            name[i] = ' '; /* exterminate tabs !!! */
            n++;
        }
        else
        {
            if (n > 0)
            {
                memmove((void*)&name[i - n], (void*)&name[i], (long long)len - (long long)i + 1); /* djb-rwth: cast operators added */
                i -= n;
                len -= n;
            }
            n = -1;
        }
    }
    if (n == len) /* empty line */
    {
        name[len = 0] = '\0';
    }
    else if (++n && n <= len)
    {
        len -= n;
        name[len] = '\0';
    }

    return len;
}
    */
    // END INCHI C FUNCTION: normalize_string
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: normalize_string
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; UCINT casts each signed char through unsigned char.
    // INCHI✔️❌: The selected C locale makes isspace true exactly for HT, LF, VT, FF, CR, and space.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: normalize_string

    let bytes = heap.slice_mut(name)?;
    let source_length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let mut len = i32::try_from(source_length)
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let mut i = 0_i32;
    let mut n = 0_i32;

    while i < len {
        let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if source_is_ascii_space(bytes[index]) {
            bytes[index] = b' ' as i8;
            n = n.wrapping_add(1);
        } else {
            if n > 0 {
                let source_start = usize::try_from(i)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let destination_start = usize::try_from(i.wrapping_sub(n))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let source_end = usize::try_from(len.wrapping_add(1))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                bytes.copy_within(source_start..source_end, destination_start);
                i = i.wrapping_sub(n);
                len = len.wrapping_sub(n);
            }
            n = -1;
        }
        i = i.wrapping_add(1);
    }
    if n == len {
        len = 0;
        bytes[0] = 0;
    } else {
        n = n.wrapping_add(1);
        if n != 0 && n <= len {
            len = len.wrapping_sub(n);
            bytes[usize::try_from(len).map_err(|_| SourceHeapError::PointerOutOfBounds)?] = 0;
        }
    }

    Ok(len)
}

#[rustfmt::skip]
pub(crate) fn is_matching_any_delim(
    character: i8,
    delimiters: Option<&[i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1710 is_matching_any_delim
    // INCHI✔️✔️: complete source frame follows verbatim; Rust performs the same allocation-free linear scan.
    /*
int is_matching_any_delim( char c, char* delims )
{
    int ic = UCINT c;
    while (*delims)
    {
        if (ic == *delims)
        {
            return 1;
        }
        delims++;
    }
    return 0;
}
    */
    // END INCHI C FUNCTION: is_matching_any_delim
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: is_matching_any_delim
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux uses signed char and UCINT casts only c through unsigned char.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: is_matching_any_delim

    let delimiters = delimiters.ok_or(SourceHeapError::NullPointer)?;
    let character = i32::from(character as u8);
    for delimiter in delimiters {
        if *delimiter == 0 {
            return Ok(0);
        }
        if character == i32::from(*delimiter) {
            return Ok(1);
        }
    }
    Err(SourceHeapError::MissingNulTerminator)
}

#[rustfmt::skip]
pub(crate) fn read_upto_delim(
    heap: &mut SourceHeap,
    pstring: &mut SourceMutPointer<i8>,
    field: SourceMutPointer<i8>,
    maxlen: i32,
    delimiters: Option<&[i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1658 read_upto_delim
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap accesses add allocation-map overhead.
    /*
int read_upto_delim( char **pstring, char *field, int maxlen, char* delims )
{
    int i, n;
    char *p = *pstring;

    if (!p)
    {
        return -1;
    }

    /* skip leading spaces */
    for (i = 0; p[i] && isspace( UCINT p[i] ); i++)
    {
        ;
    }
    p += i;

    /* read up to next delim or eol */
    n = 0;
    while (p[n] && !is_matching_any_delim( p[n], delims ))
    {
        n++;
    }

    if (n + 1 > maxlen)
    {
        return -1;
    }

    mystrncpy( field, p, n + 1 );
    field[n + 1] = '\0';

    if (!p[n])
    {
        /* reached EOL */
        *pstring = NULL;
    }
    else
    {
        /* advance reading pos */
        *pstring = *pstring + i + n;
    }

    return n;
}
    */
    // END INCHI C FUNCTION: read_upto_delim
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: read_upto_delim
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux signed-char and C-locale isspace behavior.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: read_upto_delim

    let original = *pstring;
    if original.is_null() {
        return Ok(-1);
    }
    let input = heap.slice(original.as_const())?;
    let mut leading = 0_usize;
    loop {
        let byte = *input
            .get(leading)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        if byte == 0 || !source_is_ascii_space(byte) {
            break;
        }
        leading += 1;
    }
    let mut copied = 0_usize;
    let reached_end = loop {
        let byte = *input
            .get(leading + copied)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        if byte == 0 {
            break true;
        }
        if is_matching_any_delim(byte, delimiters)? != 0 {
            break false;
        }
        copied += 1;
    };
    let copied_i32 = i32::try_from(copied).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if copied_i32.wrapping_add(1) > maxlen {
        return Ok(-1);
    }
    let source = original.offset(
        i64::try_from(leading).map_err(|_| SourceHeapError::PointerOffsetOverflow)?,
    )?;
    mystrncpy(
        heap,
        field,
        source.as_const(),
        u32::try_from(copied_i32.wrapping_add(1))
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
    )?;
    *heap
        .slice_mut(field)?
        .get_mut(copied + 1)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
    if reached_end {
        *pstring = SourceMutPointer::null();
    } else {
        *pstring = original.offset(
            i64::try_from(leading + copied)
                .map_err(|_| SourceHeapError::PointerOffsetOverflow)?,
        )?;
    }
    Ok(copied_i32)
}

#[rustfmt::skip]
pub(crate) fn remove_trailing_spaces(
    heap: &mut SourceHeap,
    pointer: SourceMutPointer<i8>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1728 remove_trailing_spaces
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap access adds allocation-map overhead.
    /*
void remove_trailing_spaces( char* p )
{
    int   len;
    for (len = (int) strlen( p ) - 1; len >= 0 && isspace( UCINT p[len] ); len--)
    {
        ;
    }
    p[++len] = '\0';
}
    */
    // END INCHI C FUNCTION: remove_trailing_spaces
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: remove_trailing_spaces
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux signed-char and C-locale isspace behavior.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: remove_trailing_spaces

    let bytes = heap.slice_mut(pointer)?;
    let length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    i32::try_from(length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let mut retained = length;
    while retained != 0 && source_is_ascii_space(bytes[retained - 1]) {
        retained -= 1;
    }
    bytes[retained] = 0;
    Ok(())
}

pub(crate) fn get_element_chemical_symbol(
    mut atomic_number: i32,
    element: &mut [i8],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:289 get_element_chemical_symbol
    // INCHI✔️❌: int get_element_chemical_symbol( int nAtNum, char *szElement )
    // INCHI✔️❌: {
    // INCHI✔️❌:     nAtNum -= 1;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (0 < nAtNum)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nAtNum += 2; /*  bypass D, T */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (0 <= nAtNum && nAtNum < nElDataLen)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* valid element symbol found */
    // INCHI✔️❌:         strcpy(szElement, ElData[nAtNum].szElName);
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* not found */
    // INCHI✔️❌:     strcpy(szElement, "??");
    // INCHI✔️❌:     return -1;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: get_element_chemical_symbol
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: get_element_chemical_symbol
    // INCHI✔️❌: const int nElDataLen = sizeof( ElData ) / sizeof( ElData[0] ) - 1;
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: get_element_chemical_symbol

    atomic_number = atomic_number.wrapping_sub(1);
    if atomic_number > 0 {
        atomic_number = atomic_number.wrapping_add(2);
    }
    let (symbol, result) = usize::try_from(atomic_number)
        .ok()
        .filter(|index| *index < EL_DATA_SYMBOLS.len())
        .and_then(|index| EL_DATA_SYMBOLS.get(index).copied())
        .map_or((b"??".as_slice(), -1), |symbol| (symbol, 0));
    let required = symbol
        .len()
        .checked_add(1)
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    let destination = element
        .get_mut(..required)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for (target, source) in destination.iter_mut().zip(symbol.iter().copied()) {
        *target = source as i8;
    }
    destination[symbol.len()] = 0;
    Ok(result)
}

pub(crate) fn get_element_or_pseudoelement_symbol(
    mut atomic_number: i32,
    element: &mut [i8],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:316 get_element_or_pseudoelement_symbol
    // INCHI✔️❌: int get_element_or_pseudoelement_symbol( int nAtNum,
    // INCHI✔️❌:                                          char *szElement )
    // INCHI✔️❌: {
    // INCHI✔️❌:     nAtNum -= 1;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (0 < nAtNum)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nAtNum += 2; /*  bypass D, T */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (0 <= nAtNum && nAtNum < nElDataLen)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* valid element symbol found */
    // INCHI✔️❌:         strcpy(szElement, ElData[nAtNum].szElName);
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!strcmp( szElement, "Zy" ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             strcpy(szElement, "Zz");
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* not found */
    // INCHI✔️❌:     strcpy(szElement, "??");
    // INCHI✔️❌:
    // INCHI✔️❌:     return -1;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: get_element_or_pseudoelement_symbol

    atomic_number = atomic_number.wrapping_sub(1);
    if atomic_number > 0 {
        atomic_number = atomic_number.wrapping_add(2);
    }
    let (symbol, result) = usize::try_from(atomic_number)
        .ok()
        .and_then(|index| EL_DATA_SYMBOLS.get(index).copied())
        .map_or((b"??".as_slice(), -1), |symbol| {
            if symbol == b"Zy" {
                (b"Zz".as_slice(), 0)
            } else {
                (symbol, 0)
            }
        });
    let required = symbol
        .len()
        .checked_add(1)
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    let destination = element
        .get_mut(..required)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for (target, source) in destination.iter_mut().zip(symbol.iter().copied()) {
        *target = source as i8;
    }
    destination[symbol.len()] = 0;
    Ok(result)
}

#[allow(non_snake_case)]
pub(crate) fn remove_one_lf(
    heap: &mut SourceHeap,
    pointer: SourceMutPointer<i8>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1740 remove_one_lf
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    void remove_one_lf( char* p )
    {
        size_t len;
        if (p && 0 < ( len = strlen( p ) ) && p[len - 1] == '\n')
        {
            p[len - 1] = '\0';
            if (len >= 2 && p[len - 2] == '\r')
            {
                p[len - 2] = '\0';
            }
        }
    }
        */
    // END INCHI C FUNCTION: remove_one_lf
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: remove_one_lf
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; size_t is 64-bit.
    // INCHI✔️❌: strlen is reproduced by the first NUL search; SourceHeap reports unterminated modeled storage explicitly.
    // INCHI✔️❌: The algorithm remains one linear scan with no allocation, but SourceHeap allocation-map lookup is absent in C.
    // END INCHI ACTIVE MACRO CONFIGURATION: remove_one_lf

    if pointer.is_null() {
        return Ok(());
    }
    let bytes = heap.slice_mut(pointer)?;
    let length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    if length != 0 && bytes[length - 1] == b'\n' as i8 {
        bytes[length - 1] = 0;
        if length >= 2 && bytes[length - 2] == b'\r' as i8 {
            bytes[length - 2] = 0;
        }
    }
    Ok(())
}

pub(crate) fn lrtrim(
    heap: &mut SourceHeap,
    pointer: SourceMutPointer<i8>,
    output_length: Option<&mut i32>,
) -> Result<SourceMutPointer<i8>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1804 lrtrim
    // INCHI✔️❌: char* lrtrim( char *p, int* nLen )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, len = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (p && ( len = (int) strlen( p ) ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = 0; i < len && is_ascii( p[i] ) && isspace( p[i] ); i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             len -= i; /* djb-rwth: variable has to be decreased before memmove */
    // INCHI✔️❌:             (memmove)(p, p + i, ((long long)len + 1)); /* djb-rwth: now cast operator can be added */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         for (; 0 < len && is_ascii( p[len - 1] ) && isspace( p[len - 1] ); len--)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         p[len] = '\0';
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nLen)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *nLen = len;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return p;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: lrtrim

    let mut length = 0_usize;
    if !pointer.is_null() {
        let bytes = heap.slice_mut(pointer)?;
        length = bytes
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        i32::try_from(length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if length != 0 {
            let leading = bytes[..length]
                .iter()
                .take_while(|byte| source_is_ascii_space(**byte))
                .count();
            if leading != 0 {
                length -= leading;
                bytes.copy_within(leading..=leading + length, 0);
            }
            while length != 0 && source_is_ascii_space(bytes[length - 1]) {
                length -= 1;
            }
            bytes[length] = 0;
        }
    }

    if let Some(output_length) = output_length {
        *output_length =
            i32::try_from(length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    }
    Ok(pointer)
}

pub(crate) fn extract_inchi_substring(
    heap: &mut SourceHeap,
    output: &mut SourceMutPointer<i8>,
    string: SourceConstPointer<i8>,
    string_length: u64,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1860 extract_inchi_substring
    // INCHI✔️❌: void extract_inchi_substring( char ** buf, const char *str, size_t slen )
    // INCHI✔️❌: {
    // INCHI✔️❌:     size_t i;
    // INCHI✔️❌:     const char *p;
    // INCHI✔️❌:     char* bufp;
    // INCHI✔️❌:     char pp;
    // INCHI✔️❌:
    // INCHI✔️❌:     bufp = *buf;
    // INCHI✔️❌:     *buf = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (str == NULL)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (strlen( str ) < 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     p = strstr( str, "InChI=" );
    // INCHI✔️❌:     if (NULL == p)
    // INCHI✔️❌:         return;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < slen; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         pp = p[i];
    // INCHI✔️❌:
    // INCHI✔️❌:         if (pp >= 'A' && pp <= 'Z')   continue;
    // INCHI✔️❌:         if (pp >= 'a' && pp <= 'z')   continue;
    // INCHI✔️❌:         if (pp >= '0' && pp <= '9')   continue;
    // INCHI✔️❌:         switch (pp)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             case '(':
    // INCHI✔️❌:             case ')':
    // INCHI✔️❌:             case '*':
    // INCHI✔️❌:             case '+':
    // INCHI✔️❌:             case ',':
    // INCHI✔️❌:             case '-':
    // INCHI✔️❌:             case '.':
    // INCHI✔️❌:             case '/':
    // INCHI✔️❌:             case ';':
    // INCHI✔️❌:             case '=':
    // INCHI✔️❌:             case '?':
    // INCHI✔️❌:             case '@':    continue;
    // INCHI✔️❌:
    // INCHI✔️❌:             default:    break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         break;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     *buf = (char*) inchi_calloc( i + 1, sizeof( char ) );
    // INCHI✔️❌:     memcpy(*buf, p, i);
    // INCHI✔️❌:     if (*buf)
    // INCHI✔️❌:         (*buf)[i] = '\0';
    // INCHI✔️❌:
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: extract_inchi_substring

    let _old_output = *output;
    *output = SourceMutPointer::null();
    if string.is_null() {
        return Ok(());
    }
    let bytes = heap.slice(string)?;
    let nul = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    if nul == 0 {
        return Ok(());
    }
    let Some(start) = bytes[..nul].windows(6).position(|window| {
        window
            == [
                b'I' as i8, b'n' as i8, b'C' as i8, b'h' as i8, b'I' as i8, b'=' as i8,
            ]
    }) else {
        return Ok(());
    };
    let limit = usize::try_from(string_length)
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    let mut length = 0_usize;
    while length < limit {
        let byte = *bytes
            .get(start + length)
            .ok_or(SourceHeapError::PointerOutOfBounds)? as u8;
        if byte.is_ascii_alphanumeric()
            || matches!(
                byte,
                b'(' | b')' | b'*' | b'+' | b',' | b'-' | b'.' | b'/' | b';' | b'=' | b'?' | b'@'
            )
        {
            length += 1;
        } else {
            break;
        }
    }
    let allocation_count = u64::try_from(length)
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
        .checked_add(1)
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    let extracted = inchi_calloc::<i8>(heap, allocation_count, 1)?;
    let source = heap.slice(string)?[start..start + length].to_vec();
    heap.slice_mut(extracted)?[..length].copy_from_slice(&source);
    heap.slice_mut(extracted)?[length] = 0;
    *output = extracted;
    Ok(())
}

#[rustfmt::skip]
pub(crate) fn extract_auxinfo_substring(
    heap: &mut SourceHeap,
    output: &mut SourceMutPointer<i8>,
    string: SourceConstPointer<i8>,
    string_length: u64,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1921 extract_auxinfo_substring
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
void extract_auxinfo_substring( char ** buf, const char *str, size_t slen )
{
    size_t i;
    const char *p;
    char* bufp;
    char pp;

    bufp = *buf;
    *buf = NULL;

    if (str == NULL)
    {
        return;
    }
    if (strlen( str ) < 1)
    {
        return;
    }

    p = strstr( str, "AuxInfo=" );
    if (NULL == p)
    {
        return;
    }

    for (i = 0; i < slen; i++)
    {
        pp = p[i];
        if (isspace( UCINT pp ))    break;
    }

    *buf = (char*) inchi_calloc( i + 1, sizeof( char ) );
    memcpy(*buf, p, i);
    if (*buf)
        (*buf)[i] = '\0';

    return;
}
    */
    // END INCHI C FUNCTION: extract_auxinfo_substring
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: extract_auxinfo_substring
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; UCINT casts through unsigned char.
    // INCHI✔️❌: The selected C locale classifies only ASCII space, tab, LF, VT, FF, and CR here.
    // INCHI✔️❌: Rust clones the source prefix before allocation to satisfy SourceHeap aliasing rules.
    // END INCHI ACTIVE MACRO CONFIGURATION: extract_auxinfo_substring

    let _old_output = *output;
    *output = SourceMutPointer::null();
    if string.is_null() {
        return Ok(());
    }
    let bytes = heap.slice(string)?;
    let nul = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    if nul == 0 {
        return Ok(());
    }
    let Some(start) = bytes[..nul].windows(8).position(|window| {
        window
            == [
                b'A' as i8,
                b'u' as i8,
                b'x' as i8,
                b'I' as i8,
                b'n' as i8,
                b'f' as i8,
                b'o' as i8,
                b'=' as i8,
            ]
    }) else {
        return Ok(());
    };
    let limit = usize::try_from(string_length)
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    let mut length = 0_usize;
    while length < limit {
        let byte = *bytes
            .get(start + length)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if source_is_ascii_space(byte) {
            break;
        }
        length += 1;
    }
    let source = bytes
        .get(start..start + length)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    let allocation_count = u64::try_from(length)
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
        .checked_add(1)
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    let extracted = inchi_calloc::<i8>(heap, allocation_count, 1)?;
    heap.slice_mut(extracted)?[..length].copy_from_slice(&source);
    heap.slice_mut(extracted)?[length] = 0;
    *output = extracted;
    Ok(())
}

pub(crate) fn mystrncpy(
    heap: &mut SourceHeap,
    target: SourceMutPointer<i8>,
    source: SourceConstPointer<i8>,
    maximum_length: u32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1760 mystrncpy
    // INCHI✔️❌: int mystrncpy( char *target, const char *source, unsigned maxlen )
    // INCHI✔️❌: {
    // INCHI✔️❌:     const char *p;
    // INCHI✔️❌:     unsigned len, source_len;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (target == NULL || maxlen == 0 || source == NULL)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* giallu: PR #163 */
    // INCHI✔️❌:     /* Find actual source length first to limit memchr search */
    // INCHI✔️❌:     source_len = (unsigned)strlen(source);
    // INCHI✔️❌:
    // INCHI✔️❌:     if (source_len < maxlen)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Source is shorter than maxlen, use actual source length */
    // INCHI✔️❌:         len = source_len;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else if ((p = (const char*)memchr(source, 0, maxlen))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* maxlen does not include the found zero termination */
    // INCHI✔️❌:         len = (int) ( p - source );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  reduced length does not include one more byte for zero termination */
    // INCHI✔️❌:         len = maxlen - 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (len)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         memmove(target, source, len);
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     memset(target + len, 0, maxlen - len); /*  zero termination */ /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:     return 1;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: mystrncpy

    if target.is_null() || maximum_length == 0 || source.is_null() {
        return Ok(0);
    }

    let copied = {
        let source_bytes = heap.slice(source)?;
        let source_length = source_bytes
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        u32::try_from(source_length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let maximum_length = usize::try_from(maximum_length)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        let copied_length = if source_length < maximum_length {
            source_length
        } else if source_bytes[..maximum_length].iter().any(|byte| *byte == 0) {
            source_bytes[..maximum_length]
                .iter()
                .position(|byte| *byte == 0)
                .expect("presence checked")
        } else {
            maximum_length - 1
        };
        let mut copied = Vec::new();
        copied
            .try_reserve_exact(copied_length + 1)
            .map_err(|_| SourceHeapError::AllocationFailed)?;
        copied.extend_from_slice(&source_bytes[..copied_length]);
        copied.push(0);
        copied
    };

    let target_bytes = heap.slice_mut(target)?;
    mystrncpy_slice(Some(target_bytes), Some(&copied), maximum_length)
}

// Storage adapter for the closed `mystrncpy` behavior when the C target is an
// inline array rather than a separately modeled source pointer.
pub(crate) fn mystrncpy_slice(
    target: Option<&mut [i8]>,
    source: Option<&[i8]>,
    maximum_length: u32,
) -> Result<i32, SourceHeapError> {
    let (Some(target), Some(source)) = (target, source) else {
        return Ok(0);
    };
    if maximum_length == 0 {
        return Ok(0);
    }
    let maximum_length = usize::try_from(maximum_length)
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    let source_length = source
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    u32::try_from(source_length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let copied_length = source_length.min(maximum_length - 1);
    let target_prefix = target
        .get_mut(..maximum_length)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    target_prefix[..copied_length].copy_from_slice(&source[..copied_length]);
    target_prefix[copied_length..].fill(0);
    Ok(1)
}

fn source_is_ascii_space(value: i8) -> bool {
    matches!(value as u8, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

fn source_is_ascii_lower(value: i8) -> Result<bool, SourceHeapError> {
    let value = value as u8;
    if !value.is_ascii() {
        return Err(SourceHeapError::InvalidSourceTextEncoding);
    }
    Ok(value.is_ascii_lowercase())
}

// Configured `ElData[].nAtMass` from util.c:103-273. The inactive
// `INCHI_ZFRAG` rows are excluded; the terminal sentinel is not part of
// `nElDataLen`.
const EL_DATA_ATOMIC_MASSES: [i32; 122] = [
    1, 2, 3, 4, 7, 9, 11, 12, 14, 16, 19, 20, 23, 24, 27, 28, 31, 32, 35, 40, 39, 40, 45, 48, 51,
    52, 55, 56, 59, 59, 64, 65, 70, 73, 75, 79, 80, 84, 85, 88, 89, 91, 93, 96, 98, 101, 103, 106,
    108, 112, 115, 119, 122, 128, 127, 131, 133, 137, 139, 140, 141, 144, 145, 150, 152, 157, 159,
    163, 165, 167, 169, 173, 175, 178, 181, 184, 186, 190, 192, 195, 197, 201, 204, 207, 209, 209,
    210, 222, 223, 226, 227, 232, 231, 238, 237, 244, 243, 247, 247, 251, 252, 257, 258, 259, 260,
    261, 270, 269, 270, 270, 278, 281, 281, 285, 278, 289, 289, 293, 297, 294, 0, 0,
];

// Configured `ElData[].nType & IS_METAL` projection from util.c:103-273.
// The inactive `INCHI_ZFRAG` rows and terminal sentinel are excluded.
const EL_DATA_IS_METAL: [bool; 122] = [
    false, false, false, false, true, true, false, false, false, false, false, false, true, true,
    true, false, false, false, false, false, true, true, true, true, true, true, true, true, true,
    true, true, true, true, false, false, false, false, false, true, true, true, true, true, true,
    true, true, true, true, true, true, true, true, true, false, false, false, true, true, true,
    true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    true, true, true, true, true, true, true, true, true, true, true, false, false, true, true,
    true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    true, true, true, true, true, true, true, true, true, true, true, true, true, true, false,
    false,
];

// Configured `ElData[].nType` projection from util.c:103-273. `METAL2` is
// distinct from `METAL` even though both satisfy the metal bit test.
const EL_DATA_TYPES: [i32; 122] = {
    let mut values = [0_i32; 122];
    let mut index = 0;
    while index < values.len() {
        if EL_DATA_IS_METAL[index] {
            values[index] = 1;
        }
        index += 1;
    }
    let metal2 = [
        26_usize, 27, 28, 29, 51, 59, 60, 63, 64, 66, 69, 70, 75, 76, 77, 78, 79, 81, 82, 83, 85,
        91, 92, 93, 94, 95, 96,
    ];
    index = 0;
    while index < metal2.len() {
        values[metal2[index]] = 3;
        index += 1;
    }
    values
};

// Configured `ElData[].bSkipAddingH` projection from util.c:103-273.
// The true runs are exactly Sc-Zn, Y-Cd, and La-Zz in the active table.
const EL_DATA_SKIP_ADDING_H: [bool; 122] = {
    let mut values = [false; 122];
    let mut index = 22;
    while index <= 31 {
        values[index] = true;
        index += 1;
    }
    index = 40;
    while index <= 49 {
        values[index] = true;
        index += 1;
    }
    index = 58;
    while index <= 121 {
        values[index] = true;
        index += 1;
    }
    values
};

const EL_DATA_SYMBOLS: [&[u8]; 122] = [
    b"H", b"D", b"T", b"He", b"Li", b"Be", b"B", b"C", b"N", b"O", b"F", b"Ne", b"Na", b"Mg",
    b"Al", b"Si", b"P", b"S", b"Cl", b"Ar", b"K", b"Ca", b"Sc", b"Ti", b"V", b"Cr", b"Mn", b"Fe",
    b"Co", b"Ni", b"Cu", b"Zn", b"Ga", b"Ge", b"As", b"Se", b"Br", b"Kr", b"Rb", b"Sr", b"Y",
    b"Zr", b"Nb", b"Mo", b"Tc", b"Ru", b"Rh", b"Pd", b"Ag", b"Cd", b"In", b"Sn", b"Sb", b"Te",
    b"I", b"Xe", b"Cs", b"Ba", b"La", b"Ce", b"Pr", b"Nd", b"Pm", b"Sm", b"Eu", b"Gd", b"Tb",
    b"Dy", b"Ho", b"Er", b"Tm", b"Yb", b"Lu", b"Hf", b"Ta", b"W", b"Re", b"Os", b"Ir", b"Pt",
    b"Au", b"Hg", b"Tl", b"Pb", b"Bi", b"Po", b"At", b"Rn", b"Fr", b"Ra", b"Ac", b"Th", b"Pa",
    b"U", b"Np", b"Pu", b"Am", b"Cm", b"Bk", b"Cf", b"Es", b"Fm", b"Md", b"No", b"Lr", b"Rf",
    b"Db", b"Sg", b"Bh", b"Hs", b"Mt", b"Ds", b"Rg", b"Cn", b"Nh", b"Fl", b"Mc", b"Lv", b"Ts",
    b"Og", b"Zy", b"Zz",
];

// Configured `ElData[].cValence` projection from util.c:103-273. The inactive
// `INCHI_ZFRAG` rows and terminal sentinel are excluded.
const EL_DATA_VALENCES: [[[i8; 5]; 5]; 122] = [
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // H
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // D
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // T
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // He
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Li
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Be
    [
        [3, 0, 0, 0, 0],
        [4, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [2, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
    ], // B
    [
        [2, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [4, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [2, 0, 0, 0, 0],
    ], // C
    [
        [1, 0, 0, 0, 0],
        [2, 0, 0, 0, 0],
        [3, 5, 0, 0, 0],
        [4, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
    ], // N
    [
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [2, 0, 0, 0, 0],
        [3, 5, 0, 0, 0],
        [4, 0, 0, 0, 0],
    ], // O
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [2, 0, 0, 0, 0],
        [3, 5, 0, 0, 0],
    ], // F
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Ne
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Na
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Mg
    [
        [3, 5, 0, 0, 0],
        [4, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [2, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
    ], // Al
    [
        [2, 0, 0, 0, 0],
        [3, 5, 0, 0, 0],
        [4, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [2, 0, 0, 0, 0],
    ], // Si
    [
        [1, 3, 5, 7, 0],
        [2, 4, 6, 0, 0],
        [3, 5, 0, 0, 0],
        [4, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
    ], // P
    [
        [0, 0, 0, 0, 0],
        [1, 3, 5, 7, 0],
        [2, 4, 6, 0, 0],
        [3, 5, 0, 0, 0],
        [4, 0, 0, 0, 0],
    ], // S
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 3, 5, 7, 0],
        [2, 4, 6, 0, 0],
        [3, 5, 0, 0, 0],
    ], // Cl
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Ar
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // K
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Ca
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Sc
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 4, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Ti
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 3, 4, 5, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // V
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 3, 6, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Cr
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 3, 4, 6, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Mn
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 3, 4, 6, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Fe
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 3, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Co
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 3, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Ni
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 2, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Cu
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Zn
    [
        [3, 5, 0, 0, 0],
        [4, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
    ], // Ga
    [
        [2, 4, 6, 0, 0],
        [3, 5, 0, 0, 0],
        [4, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Ge
    [
        [1, 3, 5, 7, 0],
        [2, 4, 6, 0, 0],
        [3, 5, 0, 0, 0],
        [4, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
    ], // As
    [
        [0, 0, 0, 0, 0],
        [1, 3, 5, 7, 0],
        [2, 4, 6, 0, 0],
        [3, 5, 0, 0, 0],
        [4, 0, 0, 0, 0],
    ], // Se
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 3, 5, 7, 0],
        [2, 4, 6, 0, 0],
        [3, 5, 0, 0, 0],
    ], // Br
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Kr
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Rb
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Sr
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Y
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [4, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Zr
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 5, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Nb
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 4, 5, 6, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Mo
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [7, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Tc
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 3, 4, 6, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Ru
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 3, 4, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Rh
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 4, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Pd
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Ag
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Cd
    [
        [3, 5, 0, 0, 0],
        [2, 4, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
    ], // In
    [
        [2, 4, 6, 0, 0],
        [3, 5, 0, 0, 0],
        [2, 4, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Sn
    [
        [1, 3, 5, 7, 0],
        [2, 4, 6, 0, 0],
        [3, 5, 0, 0, 0],
        [2, 4, 0, 0, 0],
        [3, 0, 0, 0, 0],
    ], // Sb
    [
        [0, 0, 0, 0, 0],
        [1, 3, 5, 7, 0],
        [2, 4, 6, 0, 0],
        [3, 5, 0, 0, 0],
        [2, 4, 0, 0, 0],
    ], // Te
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 3, 5, 7, 0],
        [2, 4, 6, 0, 0],
        [3, 5, 0, 0, 0],
    ], // I
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Xe
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Cs
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Ba
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // La
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 4, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Ce
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 4, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Pr
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Nd
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Pm
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 3, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Sm
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 3, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Eu
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Gd
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 4, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Tb
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Dy
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Ho
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Er
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 3, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Tm
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 3, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Yb
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Lu
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [4, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Hf
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [5, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Ta
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 4, 5, 6, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // W
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 4, 6, 7, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Re
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 3, 4, 6, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Os
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 3, 4, 6, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Ir
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 4, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Pt
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 3, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Au
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 2, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Hg
    [
        [3, 5, 0, 0, 0],
        [2, 4, 0, 0, 0],
        [1, 3, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Tl
    [
        [2, 4, 6, 0, 0],
        [3, 5, 0, 0, 0],
        [2, 4, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Pb
    [
        [1, 3, 5, 7, 0],
        [2, 4, 6, 0, 0],
        [3, 5, 0, 0, 0],
        [2, 4, 0, 0, 0],
        [3, 0, 0, 0, 0],
    ], // Bi
    [
        [0, 0, 0, 0, 0],
        [1, 3, 5, 7, 0],
        [2, 4, 6, 0, 0],
        [3, 5, 0, 0, 0],
        [2, 4, 0, 0, 0],
    ], // Po
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 3, 5, 7, 0],
        [2, 4, 6, 0, 0],
        [3, 5, 0, 0, 0],
    ], // At
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Rn
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Fr
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Ra
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Ac
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 4, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Th
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 4, 5, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Pa
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 4, 5, 6, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // U
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 4, 5, 6, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Np
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 4, 5, 6, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Pu
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 4, 5, 6, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Am
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Cm
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 4, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Bk
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Cf
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Es
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Fm
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Md
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // No
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Lr
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Rf
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Db
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Sg
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Bh
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Hs
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Mt
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Ds
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Rg
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Cn
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Nh
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Fl
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Mc
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Lv
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Ts
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Og
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Zy
    [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ], // Zz
];

pub(crate) fn el_number_in_internal_ref_table(
    element_name: Option<&[i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:347 el_number_in_internal_ref_table
    // INCHI✔❌: int el_number_in_internal_ref_table( const char* elname )
    // INCHI✔❌: {
    // INCHI✔❌:     int i;
    // INCHI✔❌:     const char *p;
    // INCHI✔❌:
    // INCHI✔❌:     for (i = 0; ( p = ElData[i].szElName )[0] && strcmp( p, elname ); i++)
    // INCHI✔❌:     {
    // INCHI✔❌:         ;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return p[0] ? i : ERR_ELEM;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: el_number_in_internal_ref_table

    let element_name = element_name.ok_or(SourceHeapError::NullPointer)?;
    let length = element_name
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let element_name = &element_name[..length];

    for (index, symbol) in EL_DATA_SYMBOLS.iter().enumerate() {
        if symbol.len() == element_name.len()
            && symbol
                .iter()
                .zip(element_name)
                .all(|(left, right)| *left == *right as u8)
        {
            return i32::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow);
        }
    }
    Ok(ERR_ELEM)
}

pub(crate) fn get_periodic_table_number(
    element_name: Option<&[i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:364 get_periodic_table_number
    // INCHI✔❌: int get_periodic_table_number( const char* elname )
    // INCHI✔❌: {
    // INCHI✔❌:     int num;
    // INCHI✔❌:
    // INCHI✔❌:     if (elname == NULL)
    // INCHI✔❌:     {
    // INCHI✔❌:         return ERR_ELEM;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (strlen(elname) == 0)
    // INCHI✔❌:     {
    // INCHI✔❌:         return ERR_ELEM;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* the single letter (common) elements */
    // INCHI✔❌:     if (!elname[1])
    // INCHI✔❌:     {
    // INCHI✔❌:         switch (elname[0])
    // INCHI✔❌:         {
    // INCHI✔❌:             case 'H':
    // INCHI✔❌:                 return EL_NUMBER_H;
    // INCHI✔❌:                 break;
    // INCHI✔❌:             case 'B':
    // INCHI✔❌:                 return EL_NUMBER_B;
    // INCHI✔❌:                 break;
    // INCHI✔❌:             case 'C':
    // INCHI✔❌:                 return EL_NUMBER_C;
    // INCHI✔❌:                 break;
    // INCHI✔❌:             case 'N':
    // INCHI✔❌:                 return EL_NUMBER_N;
    // INCHI✔❌:                 break;
    // INCHI✔❌:             case 'O':
    // INCHI✔❌:                 return EL_NUMBER_O;
    // INCHI✔❌:                 break;
    // INCHI✔❌:             case 'P':
    // INCHI✔❌:                 return EL_NUMBER_P;
    // INCHI✔❌:                 break;
    // INCHI✔❌:             case 'S':
    // INCHI✔❌:                 return EL_NUMBER_S;
    // INCHI✔❌:                 break;
    // INCHI✔❌:             case 'F':
    // INCHI✔❌:                 return EL_NUMBER_F;
    // INCHI✔❌:                 break;
    // INCHI✔❌:             case 'I':
    // INCHI✔❌:                 return EL_NUMBER_I;
    // INCHI✔❌:                 break;
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     num = el_number_in_internal_ref_table( elname );
    // INCHI✔❌:
    // INCHI✔❌:     if (num < ERR_ELEM)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* account for D,T in internal table (but not Mendeleev's table) */
    // INCHI✔❌:         num = inchi_max( 1, num - 1 );
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return num;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: get_periodic_table_number
    // BEGIN INCHI ACTIVE HEADER MACROS: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.h:46 EL_NUMBER_*
    // INCHI✔❌: #define EL_NUMBER_H  ((U_CHAR)1)
    // INCHI✔❌: #define EL_NUMBER_B  ((U_CHAR)5)
    // INCHI✔❌: #define EL_NUMBER_C  ((U_CHAR)6)
    // INCHI✔❌: #define EL_NUMBER_N  ((U_CHAR)7)
    // INCHI✔❌: #define EL_NUMBER_O  ((U_CHAR)8)
    // INCHI✔❌: #define EL_NUMBER_F  ((U_CHAR)9)
    // INCHI✔❌: #define EL_NUMBER_P  ((U_CHAR)15)
    // INCHI✔❌: #define EL_NUMBER_S  ((U_CHAR)16)
    // INCHI✔❌: #define EL_NUMBER_I  ((U_CHAR)53)
    // END INCHI ACTIVE HEADER MACROS: EL_NUMBER_*

    let Some(element_name) = element_name else {
        return Ok(ERR_ELEM);
    };
    let length = element_name
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    if length == 0 {
        return Ok(ERR_ELEM);
    }

    if length == 1 {
        let periodic_number = match element_name[0] as u8 {
            b'H' => Some(1),
            b'B' => Some(5),
            b'C' => Some(6),
            b'N' => Some(7),
            b'O' => Some(8),
            b'P' => Some(15),
            b'S' => Some(16),
            b'F' => Some(9),
            b'I' => Some(53),
            _ => None,
        };
        if let Some(periodic_number) = periodic_number {
            return Ok(periodic_number);
        }
    }

    let mut number = el_number_in_internal_ref_table(Some(element_name))?;
    if number < ERR_ELEM {
        number = 1.max(
            number
                .checked_sub(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?,
        );
    }
    Ok(number)
}

pub(crate) fn if_skip_add_H(periodic_number: i32) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:428 if_skip_add_H
    // INCHI✔️❌: int if_skip_add_H( int nPeriodicNum )
    // INCHI✔️❌: /* was called if_skip_add_H(, renamed to avoid confusion with other procedures   */
    // INCHI✔️❌: {
    // INCHI✔️❌:     return
    // INCHI✔️❌:         ElData[nPeriodicNum > 1 ? nPeriodicNum + 1 : 0].bSkipAddingH;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: if_skip_add_H
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: if_skip_add_H
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; INCHI_ZFRAG inactive
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: if_skip_add_H

    let table_index = if periodic_number > 1 {
        periodic_number
            .checked_add(1)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?
    } else {
        0
    };
    let table_index =
        usize::try_from(table_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    EL_DATA_SKIP_ADDING_H
        .get(table_index)
        .copied()
        .map(i32::from)
        .ok_or(SourceHeapError::PointerOutOfBounds)
}

pub(crate) fn get_el_valence(
    periodic_number: i32,
    charge: i32,
    valence_number: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:439 get_el_valence
    // INCHI✔❌: int get_el_valence( int nPeriodicNum, int charge, int val_num )
    // INCHI✔❌: {
    // INCHI✔❌:     if (charge < MIN_ATOM_CHARGE || charge > MAX_ATOM_CHARGE || val_num >= MAX_NUM_VALENCES)
    // INCHI✔❌:     {
    // INCHI✔❌:         return 0;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return
    // INCHI✔❌:         ElData[nPeriodicNum > 1 ? nPeriodicNum + 1 : 0].cValence[NEUTRAL_STATE + charge][val_num];
    // INCHI✔❌: }
    // END INCHI C FUNCTION: get_el_valence

    if charge < MIN_ATOM_CHARGE
        || charge > MAX_ATOM_CHARGE as i32
        || valence_number >= MAX_NUM_VALENCES as i32
    {
        return Ok(0);
    }
    let charge_index = NEUTRAL_STATE as i32 + charge;
    let charge_index =
        usize::try_from(charge_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let valence_number =
        usize::try_from(valence_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let table_index = if periodic_number > 1 {
        periodic_number
            .checked_add(1)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?
    } else {
        0
    };
    let table_index =
        usize::try_from(table_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;

    if table_index == EL_DATA_VALENCES.len() {
        return Ok(0);
    }
    let value = EL_DATA_VALENCES
        .get(table_index)
        .and_then(|charges| charges.get(charge_index))
        .and_then(|valences| valences.get(valence_number))
        .copied()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    Ok(i32::from(value))
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn get_unusual_el_valence(
    periodic_number: i32,
    charge: i32,
    radical: i32,
    bonds_valence: i32,
    hydrogen_count: i32,
    bond_count: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:454 get_unusual_el_valence
    // INCHI✔️❌: complete behavior; checked Rust table access adds modeled-pointer validation.
    /*
    int get_unusual_el_valence( int nPeriodicNum,
                                int charge,
                                int radical,
                                int bonds_valence,
                                int num_H,
                                int num_bonds )
    {
        int i, num_found, chem_valence, rad_adj, known_chem_valence, exact_found;

        if (!num_bonds && !num_H)
        {
            return 0;
        }

        if (charge < MIN_ATOM_CHARGE || charge > MAX_ATOM_CHARGE)
        {
            if (bonds_valence == num_bonds)
            {
                return 0; /* all single bonds */
            }
            return bonds_valence;
        }

        if (!get_el_valence( nPeriodicNum, charge, 0 ) && bonds_valence == num_bonds)
        {
            return 0;
        }

        chem_valence = bonds_valence + num_H;
        rad_adj = 0;
        num_found = 0;
        exact_found = 0;

        /* Take into account a radical */
        if (radical == RADICAL_DOUBLET)
        {
            rad_adj = 1;
        }
        else if (radical == RADICAL_TRIPLET)
        {
            rad_adj = 2;
        }

        for (i = 0; i < MAX_NUM_VALENCES; i++)
        {
            if (0 < ( known_chem_valence = get_el_valence( nPeriodicNum, charge, i ) - rad_adj ) &&
                 num_bonds <= known_chem_valence && known_chem_valence <= chem_valence)
            {
                num_found++;
                if (known_chem_valence == chem_valence)
                {
                    exact_found = 1;
                    break;
                }
            }
        }

        return ( exact_found && 1 == num_found ) ? 0 : chem_valence;
    }
    */
    // END INCHI C FUNCTION: get_unusual_el_valence
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: get_unusual_el_valence
    // INCHI✔️❌: #define MIN_ATOM_CHARGE (-2)
    // INCHI✔️❌: #define MAX_ATOM_CHARGE 2
    // INCHI✔️❌: #define MAX_NUM_VALENCES 5
    // INCHI✔️❌: #define RADICAL_DOUBLET 2
    // INCHI✔️❌: #define RADICAL_TRIPLET 3
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: get_unusual_el_valence

    if bond_count == 0 && hydrogen_count == 0 {
        return Ok(0);
    }
    if charge < MIN_ATOM_CHARGE || charge > MAX_ATOM_CHARGE as i32 {
        return Ok(if bonds_valence == bond_count {
            0
        } else {
            bonds_valence
        });
    }
    if get_el_valence(periodic_number, charge, 0)? == 0 && bonds_valence == bond_count {
        return Ok(0);
    }

    let chemical_valence = bonds_valence
        .checked_add(hydrogen_count)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let radical_adjustment = if radical == RADICAL_DOUBLET as i32 {
        1
    } else if radical == RADICAL_TRIPLET as i32 {
        2
    } else {
        0
    };
    let mut found_count = 0_i32;
    let mut exact_found = false;
    for index in 0..MAX_NUM_VALENCES as i32 {
        let known_valence = get_el_valence(periodic_number, charge, index)? - radical_adjustment;
        if known_valence > 0 && bond_count <= known_valence && known_valence <= chemical_valence {
            found_count = found_count.wrapping_add(1);
            if known_valence == chemical_valence {
                exact_found = true;
                break;
            }
        }
    }

    if exact_found && found_count == 1 {
        Ok(0)
    } else {
        Ok(chemical_valence)
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn needed_unusual_el_valence(
    periodic_number: i32,
    charge: i32,
    radical: i32,
    bonds_valence: i32,
    actual_bonds_valence: i32,
    hydrogen_count: i32,
    bond_count: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:518 needed_unusual_el_valence
    // INCHI✔️❌: int needed_unusual_el_valence( int nPeriodicNum,
    // INCHI✔️❌:                                int charge,
    // INCHI✔️❌:                                int radical,
    // INCHI✔️❌:                                int bonds_valence,
    // INCHI✔️❌:                                int actual_bonds_valence,
    // INCHI✔️❌:                                int num_H, int
    // INCHI✔️❌:                                num_bonds )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int chem_valence, num_H_expected; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:     char szElement[4];
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( !num_bonds && !num_H )
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (num_bonds && get_element_chemical_symbol( nPeriodicNum, szElement ) != -1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         num_H_expected = get_num_H( szElement, 0, NULL, charge, radical, actual_bonds_valence, 0, 0, 0, 0 );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         num_H_expected = num_H;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     chem_valence = bonds_valence + num_H;
    // INCHI✔️❌:
    // INCHI❌❌: #if ( (BUILD_WITH_ENG_OPTIONS==1) && (SDF_OUTPUT_HETERO_VALENCE==1) )
    // INCHI❌❌:     if ((nPeriodicNum == 1 && chem_valence != 1) /* H */ || (nPeriodicNum == 6 && chem_valence != 4) /* C */ ||
    // INCHI❌❌:         (nPeriodicNum != 1 && nPeriodicNum != 6) || charge || radical) /* djb-rwth: addressing LLVM warning */
    // INCHI❌❌:     {
    // INCHI❌❌:         return chem_valence ? chem_valence : -1;
    // INCHI❌❌:     }
    // INCHI❌❌:     else
    // INCHI❌❌:     {
    // INCHI❌❌:         return 0;
    // INCHI❌❌:     }
    // INCHI✔️❌: #else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int i, num_found, num_found_known, rad_adj, known_chem_valence, exact_found;
    // INCHI✔️❌:
    // INCHI✔️❌:         if (charge < MIN_ATOM_CHARGE || charge > MAX_ATOM_CHARGE ||
    // INCHI✔️❌:              !get_el_valence( nPeriodicNum, charge, 0 ) ||
    // INCHI✔️❌:              if_skip_add_H( nPeriodicNum ) || bonds_valence != actual_bonds_valence ||
    // INCHI✔️❌:              num_H_expected != num_H)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!num_H && !num_H_expected && bonds_valence == actual_bonds_valence)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 0; /* no H */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return chem_valence; /* needs to add H-atoms */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* take into account radical */
    // INCHI✔️❌:         if (radical == RADICAL_DOUBLET)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             rad_adj = 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (radical == RADICAL_TRIPLET)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             rad_adj = 2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             rad_adj = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         num_found_known = 0;
    // INCHI✔️❌:         num_found = 0;
    // INCHI✔️❌:         exact_found = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = 0; i < MAX_NUM_VALENCES; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (0 < ( known_chem_valence = get_el_valence( nPeriodicNum, charge, i ) ) &&
    // INCHI✔️❌:                  bonds_valence <= ( known_chem_valence -= rad_adj ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* found known valence that fits without H */
    // INCHI✔️❌:                 num_found_known++;
    // INCHI✔️❌:                 if (known_chem_valence <= chem_valence)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* known valence is large enough to accommodate (implicit) H */
    // INCHI✔️❌:                     num_found++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (known_chem_valence == chem_valence)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     exact_found = 1;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         return ( exact_found && 1 == num_found && 1 == num_found_known )
    // INCHI✔️❌:             ? 0
    // INCHI✔️❌:             : chem_valence ? chem_valence : -1;    /* needs zero */
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: needed_unusual_el_valence
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: needed_unusual_el_valence
    // INCHI✔️❌: #define BUILD_WITH_ENG_OPTIONS 0
    // INCHI✔️❌: #define MIN_ATOM_CHARGE (-2)
    // INCHI✔️❌: #define MAX_ATOM_CHARGE 2
    // INCHI✔️❌: #define MAX_NUM_VALENCES 5
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: needed_unusual_el_valence

    let mut element = [0_i8; 4];
    let expected_hydrogen_count =
        if bond_count != 0 && get_element_chemical_symbol(periodic_number, &mut element)? != -1 {
            get_num_H(
                Some(&element),
                0,
                None,
                charge,
                radical,
                actual_bonds_valence,
                0,
                0,
                0,
                0,
            )?
        } else {
            hydrogen_count
        };
    let chemical_valence = bonds_valence.wrapping_add(hydrogen_count);

    if charge < MIN_ATOM_CHARGE
        || charge > MAX_ATOM_CHARGE as i32
        || get_el_valence(periodic_number, charge, 0)? == 0
        || if_skip_add_H(periodic_number)? != 0
        || bonds_valence != actual_bonds_valence
        || expected_hydrogen_count != hydrogen_count
    {
        if hydrogen_count == 0
            && expected_hydrogen_count == 0
            && bonds_valence == actual_bonds_valence
        {
            return Ok(0);
        }
        return Ok(chemical_valence);
    }

    let radical_adjustment = if radical == RADICAL_DOUBLET as i32 {
        1
    } else if radical == RADICAL_TRIPLET as i32 {
        2
    } else {
        0
    };
    let mut known_count = 0_i32;
    let mut fitting_count = 0_i32;
    let mut exact_found = false;
    for index in 0..MAX_NUM_VALENCES as i32 {
        let mut known_valence = get_el_valence(periodic_number, charge, index)?;
        if known_valence > 0 {
            known_valence = known_valence.wrapping_sub(radical_adjustment);
            if bonds_valence <= known_valence {
                known_count = known_count.wrapping_add(1);
                if known_valence <= chemical_valence {
                    fitting_count = fitting_count.wrapping_add(1);
                }
                if known_valence == chemical_valence {
                    exact_found = true;
                    break;
                }
            }
        }
    }
    if exact_found && fitting_count == 1 && known_count == 1 {
        Ok(0)
    } else if chemical_valence != 0 {
        Ok(chemical_valence)
    } else {
        Ok(-1)
    }
}

pub(crate) fn detect_unusual_el_valence(
    periodic_number: i32,
    charge: i32,
    radical: i32,
    bonds_valence: i32,
    hydrogen_count: i32,
    bond_count: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:620 detect_unusual_el_valence
    // INCHI✔❌: int detect_unusual_el_valence( int nPeriodicNum,
    // INCHI✔❌:                                int charge,
    // INCHI✔❌:                                int radical,
    // INCHI✔❌:                                int bonds_valence,
    // INCHI✔❌:                                int num_H,
    // INCHI✔❌:                                int num_bonds )
    // INCHI✔❌: {
    // INCHI✔❌:     int i, chem_valence, rad_adj, known_chem_valence;
    // INCHI✔❌:
    // INCHI✔❌:     if (!num_bonds && !num_H)
    // INCHI✔❌:     {
    // INCHI✔❌:         return 0;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (charge < MIN_ATOM_CHARGE || charge > MAX_ATOM_CHARGE)
    // INCHI✔❌:     {
    // INCHI✔❌:         if (bonds_valence == num_bonds)
    // INCHI✔❌:         {
    // INCHI✔❌:             return 0; /* all single bonds */
    // INCHI✔❌:         }
    // INCHI✔❌:         return bonds_valence;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (!get_el_valence( nPeriodicNum, charge, 0 ) && bonds_valence == num_bonds)
    // INCHI✔❌:     {
    // INCHI✔❌:         return 0;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     chem_valence = bonds_valence + num_H;
    // INCHI✔❌:     rad_adj = 0;
    // INCHI✔❌:
    // INCHI✔❌:     /* take into account radical */
    // INCHI✔❌:     if (radical == RADICAL_DOUBLET)
    // INCHI✔❌:     {
    // INCHI✔❌:         rad_adj = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (radical == RADICAL_TRIPLET || radical == RADICAL_SINGLET)
    // INCHI✔❌:     {
    // INCHI✔❌:         rad_adj = 2;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     for (i = 0; i < MAX_NUM_VALENCES; i++)
    // INCHI✔❌:     {
    // INCHI✔❌:         if (0 < ( known_chem_valence = get_el_valence( nPeriodicNum, charge, i ) - rad_adj ))
    // INCHI✔❌:         {
    // INCHI✔❌:             if (known_chem_valence == chem_valence)
    // INCHI✔❌:             {
    // INCHI✔❌:                 return 0;
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return chem_valence;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: detect_unusual_el_valence

    if bond_count == 0 && hydrogen_count == 0 {
        return Ok(0);
    }
    if charge < MIN_ATOM_CHARGE || charge > MAX_ATOM_CHARGE as i32 {
        return Ok(if bonds_valence == bond_count {
            0
        } else {
            bonds_valence
        });
    }
    if get_el_valence(periodic_number, charge, 0)? == 0 && bonds_valence == bond_count {
        return Ok(0);
    }
    let chemical_valence = bonds_valence
        .checked_add(hydrogen_count)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let radical_adjustment = match radical {
        2 => 1,
        1 | 3 => 2,
        _ => 0,
    };
    for index in 0..MAX_NUM_VALENCES as i32 {
        let known = get_el_valence(periodic_number, charge, index)? - radical_adjustment;
        if known > 0 && known == chemical_valence {
            return Ok(0);
        }
    }
    Ok(chemical_valence)
}

fn source_strtol_base_10(bytes: &[i8], start: usize) -> Result<(i32, usize), SourceHeapError> {
    let mut cursor = start;
    while bytes
        .get(cursor)
        .is_some_and(|byte| source_is_ascii_space(*byte))
    {
        cursor = cursor
            .checked_add(1)
            .ok_or(SourceHeapError::PointerOffsetOverflow)?;
    }
    let mut negative = false;
    if let Some(sign) = bytes.get(cursor).map(|byte| *byte as u8) {
        if sign == b'+' || sign == b'-' {
            negative = sign == b'-';
            cursor = cursor
                .checked_add(1)
                .ok_or(SourceHeapError::PointerOffsetOverflow)?;
        }
    }
    let first_digit = cursor;
    let limit = if negative {
        i64::MAX as u128 + 1
    } else {
        i64::MAX as u128
    };
    let mut magnitude = 0_u128;
    while let Some(digit) = bytes
        .get(cursor)
        .map(|byte| *byte as u8)
        .filter(u8::is_ascii_digit)
    {
        magnitude = magnitude
            .saturating_mul(10)
            .saturating_add(u128::from(digit - b'0'))
            .min(limit);
        cursor = cursor
            .checked_add(1)
            .ok_or(SourceHeapError::PointerOffsetOverflow)?;
    }
    if cursor == first_digit {
        return Ok((0, start));
    }
    let value = if negative {
        if magnitude == i64::MAX as u128 + 1 {
            i64::MIN
        } else {
            -(magnitude as i64)
        }
    } else {
        magnitude as i64
    };
    Ok((value as i32, cursor))
}

pub(crate) fn extract_charges_and_radicals(
    element_name: Option<&mut [i8]>,
    radical: Option<&mut i32>,
    charge: Option<&mut i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:700 extract_charges_and_radicals
    // INCHI✔❌: int extract_charges_and_radicals( char *elname, int *pnRadical, int *pnCharge )
    // INCHI✔❌: {
    // INCHI✔❌:     char *q, *r, *p;
    // INCHI✔❌:     int  nCharge = 0, nRad = 0, charge_len = 0, k, nVal, nSign, nLastSign = 1; /* djb-rwth: removing redundant variables */
    // INCHI✔❌:
    // INCHI✔❌:     p = elname;
    // INCHI✔❌:
    // INCHI✔❌:     /*  extract radicals & charges */
    // INCHI✔❌:     while ((q = strpbrk(p, "+-^"))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:     {
    // INCHI✔❌:         switch (*q)
    // INCHI✔❌:         {
    // INCHI✔❌:             case '+':
    // INCHI✔❌:             case '-':
    // INCHI✔❌:                 for (k = 0, nVal = 0; ( nSign = ( '+' == q[k] ) ) || ( nSign = -( '-' == q[k] ) ); k++)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     nVal += ( nLastSign = nSign );
    // INCHI✔❌:                     charge_len++;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 if ((nSign = (int)strtol(q + k, &r, 10))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:                 {
    // INCHI✔❌:                     /*  fixed 12-5-2001 */
    // INCHI✔❌:                     nVal += nLastSign * ( nSign - 1 );
    // INCHI✔❌:                 }
    // INCHI✔❌:                 charge_len = (int) ( r - q );
    // INCHI✔❌:                 nCharge += nVal;
    // INCHI✔❌:                 break;
    // INCHI✔❌:             /* case '.': */ /*  singlet '.' may be confused with '.' in formulas like CaO.H2O */
    // INCHI✔❌:             case '^':
    // INCHI✔❌:                 nRad = 1; /* doublet here is 1. See below */
    // INCHI✔❌:                 charge_len = 1;
    // INCHI✔❌:                 for (k = 1; q[0] == q[k]; k++)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     nRad++;
    // INCHI✔❌:                     charge_len++;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 break;
    // INCHI✔❌:         }
    // INCHI✔❌:         memmove(q, q + charge_len, strlen(q + charge_len) + 1);
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* djb-rwth: removing redundant code */
    // INCHI✔❌:
    // INCHI✔❌:     /*  radical */
    // INCHI✔❌:     if (( q = strrchr( p, ':' ) ) && !q[1])
    // INCHI✔❌:     {
    // INCHI✔❌:         nRad = RADICAL_SINGLET;
    // INCHI✔❌:         q[0] = '\0';
    // INCHI✔❌:         /* djb-rwth: removing redundant code */
    // INCHI✔❌:     }
    // INCHI✔❌:     else
    // INCHI✔❌:     {
    // INCHI✔❌:         while (( q = strrchr( p, '.' ) ) && !q[1])
    // INCHI✔❌:         {
    // INCHI✔❌:             nRad++;
    // INCHI✔❌:             q[0] = '\0';
    // INCHI✔❌:             /* djb-rwth: removing redundant code */
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         nRad = nRad == 1 ? RADICAL_DOUBLET :
    // INCHI✔❌:             nRad == 2 ? RADICAL_TRIPLET : 0;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     *pnRadical = nRad;
    // INCHI✔❌:     *pnCharge = nCharge;
    // INCHI✔❌:
    // INCHI✔❌:     return ( nRad || nCharge );
    // INCHI✔❌: }
    // END INCHI C FUNCTION: extract_charges_and_radicals

    let bytes = element_name.ok_or(SourceHeapError::NullPointer)?;
    let radical = radical.ok_or(SourceHeapError::NullPointer)?;
    let charge = charge.ok_or(SourceHeapError::NullPointer)?;
    let mut length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let mut total_charge = 0_i32;
    let mut radical_count = 0_i32;
    let mut charge_length = 0_i32;

    while let Some(token) = bytes[..length]
        .iter()
        .position(|byte| matches!(*byte as u8, b'+' | b'-' | b'^'))
    {
        if matches!(bytes[token] as u8, b'+' | b'-') {
            let mut offset = 0_usize;
            let mut value = 0_i32;
            let mut last_sign = 1_i32;
            while let Some(sign) = bytes
                .get(
                    token
                        .checked_add(offset)
                        .ok_or(SourceHeapError::PointerOffsetOverflow)?,
                )
                .map(|byte| *byte as u8)
            {
                let sign = match sign {
                    b'+' => 1,
                    b'-' => -1,
                    _ => break,
                };
                last_sign = sign;
                value = value
                    .checked_add(sign)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                charge_length = charge_length
                    .checked_add(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                offset = offset
                    .checked_add(1)
                    .ok_or(SourceHeapError::PointerOffsetOverflow)?;
            }
            let number_start = token
                .checked_add(offset)
                .ok_or(SourceHeapError::PointerOffsetOverflow)?;
            let (parsed, end) = source_strtol_base_10(&bytes[..length], number_start)?;
            if parsed != 0 {
                let adjustment = last_sign
                    .checked_mul(
                        parsed
                            .checked_sub(1)
                            .ok_or(SourceHeapError::SourceIntegerOverflow)?,
                    )
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                value = value
                    .checked_add(adjustment)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            }
            charge_length =
                i32::try_from(end - token).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            total_charge = total_charge
                .checked_add(value)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        } else {
            radical_count = 1;
            charge_length = 1;
            loop {
                let run_length = usize::try_from(charge_length)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let next = token
                    .checked_add(run_length)
                    .ok_or(SourceHeapError::PointerOffsetOverflow)?;
                if next >= length || bytes[token] != bytes[next] {
                    break;
                }
                radical_count = radical_count
                    .checked_add(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                charge_length = charge_length
                    .checked_add(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            }
        }
        let removed =
            usize::try_from(charge_length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let source = token
            .checked_add(removed)
            .ok_or(SourceHeapError::PointerOffsetOverflow)?;
        if source > length {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        bytes.copy_within(source..=length, token);
        length = length
            .checked_sub(removed)
            .ok_or(SourceHeapError::PointerOffsetOverflow)?;
    }

    if let Some(position) = bytes[..length].iter().rposition(|byte| *byte as u8 == b':')
        && position
            .checked_add(1)
            .ok_or(SourceHeapError::PointerOffsetOverflow)?
            == length
    {
        radical_count = RADICAL_SINGLET as i32;
        bytes[position] = 0;
    } else {
        while length > 0 && bytes[length - 1] as u8 == b'.' {
            radical_count = radical_count
                .checked_add(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            length -= 1;
            bytes[length] = 0;
        }
        radical_count = if radical_count == 1 {
            RADICAL_DOUBLET as i32
        } else if radical_count == 2 {
            RADICAL_TRIPLET as i32
        } else {
            0
        };
    }
    *radical = radical_count;
    *charge = total_charge;
    Ok(i32::from(radical_count != 0 || total_charge != 0))
}

pub(crate) fn extract_h_atoms(
    element_name: Option<&mut [i8]>,
    isotope_hydrogens: Option<&mut [S_CHAR]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:774 extract_H_atoms
    // INCHI✔❌: int extract_H_atoms( char *elname, S_CHAR num_iso_H[] )
    // INCHI✔❌: {
    // INCHI✔❌:     int i, len, c, k, num_H, val;
    // INCHI✔❌:     char *q;
    // INCHI✔❌:     char elname1 = '\0';
    // INCHI✔❌:
    // INCHI✔❌:     i = 0;
    // INCHI✔❌:     num_H = 0;
    // INCHI✔❌:     len = (int) strlen( elname );
    // INCHI✔❌:     c = UCINT elname[0];
    // INCHI✔❌:
    // INCHI✔❌:     if (len > 1)
    // INCHI✔❌:     {
    // INCHI✔❌:         elname1 = elname[1];
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     while (i < len)
    // INCHI✔❌:     {
    // INCHI✔❌:         switch (c)
    // INCHI✔❌:         {
    // INCHI✔❌:             case 'H':
    // INCHI✔❌:                 k = 0;
    // INCHI✔❌:                 break;
    // INCHI✔❌:             case 'D':
    // INCHI✔❌:                 k = 1;
    // INCHI✔❌:                 break;
    // INCHI✔❌:             case 'T':
    // INCHI✔❌:                 k = 2;
    // INCHI✔❌:                 break;
    // INCHI✔❌:             default:
    // INCHI✔❌:                 k = -1;
    // INCHI✔❌:                 break;
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         q = elname + i + 1; /*  pointer to the next to elname[i] character */
    // INCHI✔❌:         c = UCINT q[0];
    // INCHI✔❌:
    // INCHI✔❌:         if (k >= 0 && !islower( c ))
    // INCHI✔❌:         {
    // INCHI✔❌:             /*  found a hydrogen */
    // INCHI✔❌:             if (isdigit( c ))
    // INCHI✔❌:             {
    // INCHI✔❌:                 val = (int) strtol( q, &q, 10 );
    // INCHI✔❌:                 /*  q = pointer to the next to number of hydrogen atom(s) character */
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 val = 1;
    // INCHI✔❌:             }
    // INCHI✔❌:             if (k)
    // INCHI✔❌:             {
    // INCHI✔❌:                 num_iso_H[k] += val;
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 num_H += val;
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             /*  remove the hydrogen atom from the string */
    // INCHI✔❌:             len -= (int) ( q - elname ) - i;
    // INCHI✔❌:             memmove(elname + i, q, (long long)len + 1); /* djb-rwth: cast operator added */
    // INCHI✔❌:             /*  c =  UCINT elname[i]; */
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:         {
    // INCHI✔❌:             i++;
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         c = UCINT elname[i]; /*  moved here 11-04-2002 */
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     len = (int) strlen( elname );
    // INCHI✔❌:     if (len == 2)
    // INCHI✔❌:     {
    // INCHI✔❌:         if (elname[1] != elname1)
    // INCHI✔❌:             /* Error, incorrect 2nd char of elname appears after 'subtracting' {H,D,T}  */
    // INCHI✔❌:             /* See a bug reported to inchi-discuss by A. Dalke for alias atom "pH4d"    */
    // INCHI✔❌:             /*^^^ 2017-01-06                                                            */
    // INCHI✔❌:             elname[1] = '?';
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return num_H;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: extract_H_atoms

    let bytes = element_name.ok_or(SourceHeapError::NullPointer)?;
    let isotope_hydrogens = isotope_hydrogens.ok_or(SourceHeapError::NullPointer)?;
    if isotope_hydrogens.len() < 3 {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let mut length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    i32::try_from(length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let original_second = if length > 1 { bytes[1] } else { 0 };
    let mut index = 0_usize;
    let mut hydrogen_count = 0_i32;
    let mut current = *bytes.first().ok_or(SourceHeapError::PointerOutOfBounds)?;

    while index < length {
        let isotope_index = match current as u8 {
            b'H' => Some(0_usize),
            b'D' => Some(1_usize),
            b'T' => Some(2_usize),
            _ => None,
        };
        let mut next = index
            .checked_add(1)
            .ok_or(SourceHeapError::PointerOffsetOverflow)?;
        let next_character = *bytes.get(next).ok_or(SourceHeapError::PointerOutOfBounds)?;

        if let Some(isotope_index) = isotope_index
            && !source_is_ascii_lower(next_character)?
        {
            let value = if (next_character as u8).is_ascii_digit() {
                let (value, end) = source_strtol_base_10(&bytes[..length], next)?;
                next = end;
                value
            } else {
                1
            };
            if isotope_index == 0 {
                hydrogen_count = hydrogen_count
                    .checked_add(value)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            } else {
                let accumulated = i32::from(isotope_hydrogens[isotope_index])
                    .checked_add(value)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                isotope_hydrogens[isotope_index] = accumulated as S_CHAR;
            }

            let consumed = next
                .checked_sub(index)
                .ok_or(SourceHeapError::PointerOffsetOverflow)?;
            length = length
                .checked_sub(consumed)
                .ok_or(SourceHeapError::PointerOffsetOverflow)?;
            let copied_length = length
                .checked_add(1)
                .ok_or(SourceHeapError::PointerOffsetOverflow)?;
            let source_end = next
                .checked_add(copied_length)
                .ok_or(SourceHeapError::PointerOffsetOverflow)?;
            let target_end = index
                .checked_add(copied_length)
                .ok_or(SourceHeapError::PointerOffsetOverflow)?;
            if source_end > bytes.len() || target_end > bytes.len() {
                return Err(SourceHeapError::PointerOutOfBounds);
            }
            bytes.copy_within(next..source_end, index);
        } else {
            index = index
                .checked_add(1)
                .ok_or(SourceHeapError::PointerOffsetOverflow)?;
        }
        current = *bytes
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
    }

    let final_length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    i32::try_from(final_length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if final_length == 2 && bytes[1] != original_second {
        bytes[1] = b'?' as i8;
    }
    Ok(hydrogen_count)
}

pub(crate) fn is_el_a_metal(periodic_number: i32) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:688 is_el_a_metal
    // INCHI✔❌: int is_el_a_metal( int nPeriodicNum )
    // INCHI✔❌: {
    // INCHI✔❌:     return 0 != ( ElData[nPeriodicNum + 1].nType & IS_METAL );
    // INCHI✔❌: }
    // END INCHI C FUNCTION: is_el_a_metal

    let table_index = periodic_number
        .checked_add(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let table_index =
        usize::try_from(table_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if table_index == EL_DATA_IS_METAL.len() {
        return Ok(0);
    }
    EL_DATA_IS_METAL
        .get(table_index)
        .copied()
        .map(i32::from)
        .ok_or(SourceHeapError::PointerOutOfBounds)
}

pub(crate) fn get_el_type(periodic_number: i32) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:679 get_el_type
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: get_el_type
    // INCHI✔️✔️: int get_el_type( int nPeriodicNum )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     return ElData[nPeriodicNum + 1].nType;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: get_el_type
    // END INCHI C FUNCTION: get_el_type
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: get_el_type
    // INCHI✔️✔️: #define METAL  1
    // INCHI✔️✔️: #define METAL2 3
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: get_el_type

    let table_index = periodic_number
        .checked_add(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let table_index =
        usize::try_from(table_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if table_index == EL_DATA_TYPES.len() {
        return Ok(0);
    }
    EL_DATA_TYPES
        .get(table_index)
        .copied()
        .ok_or(SourceHeapError::PointerOutOfBounds)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn get_num_H(
    element_name: Option<&[i8]>,
    input_hydrogens: i32,
    input_isotope_hydrogens: Option<&[S_CHAR]>,
    charge: i32,
    radical: i32,
    chemical_bond_valence: i32,
    atom_input_valence: i32,
    aliased: i32,
    do_not_add_hydrogen: i32,
    has_metal_neighbor: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:862 get_num_H
    // INCHI✔❌: int get_num_H( const char* elname,
    // INCHI✔❌:                 int inp_num_H,
    // INCHI✔❌:                 S_CHAR inp_num_iso_H[],
    // INCHI✔❌:                 int charge,
    // INCHI✔❌:                 int radical,
    // INCHI✔❌:                 int chem_bonds_valence,
    // INCHI✔❌:                 int atom_input_valence,
    // INCHI✔❌:                 int bAliased,
    // INCHI✔❌:                 int bDoNotAddH,
    // INCHI✔❌:                 int bHasMetalNeighbor )
    // INCHI✔❌: {
    // INCHI✔❌:     int val, i, el_number, num_H = 0, num_iso_H;
    // INCHI✔❌:     static int intl_el_number_N = 0, intl_el_number_S=0, intl_el_number_O=0, intl_el_number_C=0;
    // INCHI✔❌:
    // INCHI✔❌:     if (!intl_el_number_N)
    // INCHI✔❌:     {
    // INCHI✔❌:         intl_el_number_N = el_number_in_internal_ref_table( "N" );
    // INCHI✔❌:     }
    // INCHI✔❌:     if (!intl_el_number_S)
    // INCHI✔❌:     {
    // INCHI✔❌:         intl_el_number_S = el_number_in_internal_ref_table( "S" );
    // INCHI✔❌:     }
    // INCHI✔❌:     if (!intl_el_number_O)
    // INCHI✔❌:     {
    // INCHI✔❌:         intl_el_number_O = el_number_in_internal_ref_table( "O" );
    // INCHI✔❌:     }
    // INCHI✔❌:     if (!intl_el_number_C)
    // INCHI✔❌:     {
    // INCHI✔❌:         intl_el_number_C = el_number_in_internal_ref_table( "C" );
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌:     /*  atom_input_valence (cValence) cannot be specified in case of */
    // INCHI✔❌:     /*  aliased MOLFile atom with known inp_num_H or inp_num_iso_H[] */
    // INCHI✔❌:
    // INCHI✔❌:     if (bAliased)
    // INCHI✔❌:     {
    // INCHI✔❌:         num_H = inp_num_H;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (atom_input_valence && ( atom_input_valence != 15 || chem_bonds_valence ))
    // INCHI✔❌:     {
    // INCHI✔❌:         num_H = inchi_max( 0, atom_input_valence - chem_bonds_valence );
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (atom_input_valence == 15 && !chem_bonds_valence)
    // INCHI✔❌:     {
    // INCHI✔❌:         num_H = 0;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (MIN_ATOM_CHARGE <= charge &&
    // INCHI✔❌:               MAX_ATOM_CHARGE >= charge &&
    // INCHI✔❌:               ERR_ELEM != ( el_number = el_number_in_internal_ref_table( elname ) ) &&
    // INCHI✔❌:               !ElData[el_number].bSkipAddingH && !bDoNotAddH)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* add hydrogen atoms according to standard element valence */
    // INCHI✔❌:         if (radical && radical != RADICAL_SINGLET)
    // INCHI✔❌:         {
    // INCHI✔❌:             if ((val = ElData[el_number].cValence[NEUTRAL_STATE + charge][0])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:             {
    // INCHI✔❌:                 val -= ( radical == RADICAL_DOUBLET ) ? 1
    // INCHI✔❌:                     : ( radical == RADICAL_SINGLET || radical == RADICAL_TRIPLET ) ? 2 : val;
    // INCHI✔❌:                 /* if unknown radical then do not add H */
    // INCHI✔❌:                 num_H = inchi_max( 0, val - chem_bonds_valence );
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:         {
    // INCHI✔❌:             /* find the smallest valence that is greater than the sum of the chemical bond valences */
    // INCHI✔❌:             for (i = 0;
    // INCHI✔❌:                  ( val = ElData[el_number].cValence[NEUTRAL_STATE + charge][i] ) &&
    // INCHI✔❌:                  val < chem_bonds_valence;
    // INCHI✔❌:                  i++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ;
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             /* special case: do not add H to N(IV), S(III), S+(II), S-(II) */ /* S ions added 2004-05-10 */
    // INCHI✔❌:             if (el_number == intl_el_number_N && !charge && !radical && val == 5)
    // INCHI✔❌:             {
    // INCHI✔❌:                 val = 3;
    // INCHI✔❌:             }
    // INCHI✔❌:             /*else if ( el_number == el_number_N && !charge && !radical && val == 3 &&
    // INCHI✔❌:                  chem_bonds_valence == 2 && bHasMetalNeighbor )
    // INCHI✔❌:               {
    // INCHI✔❌:               val = 2;
    // INCHI✔❌:               }
    // INCHI✔❌:             */
    // INCHI✔❌:             else if (el_number == intl_el_number_S && !charge && !radical && val == 4 && chem_bonds_valence == 3)
    // INCHI✔❌:             {
    // INCHI✔❌:                 val = 3;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (bHasMetalNeighbor && el_number != intl_el_number_C && val > 0)
    // INCHI✔❌:             {
    // INCHI✔❌:                 val--;
    // INCHI✔❌:             }
    // INCHI✔❌:             /*
    // INCHI✔❌:             if ( (el_number == el_number_S || el_number == el_number_O) &&
    // INCHI✔❌:                  abs(charge)==1 && !radical && val == 3 && chem_bonds_valence == 2 && bHasMetalNeighbor )
    // INCHI✔❌:               {
    // INCHI✔❌:               val = 2;
    // INCHI✔❌:               }
    // INCHI✔❌:             else
    // INCHI✔❌:             */
    // INCHI✔❌:
    // INCHI✔❌:             num_H = inchi_max( 0, val - chem_bonds_valence );
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         num_iso_H = 0;
    // INCHI✔❌:         if (inp_num_iso_H)
    // INCHI✔❌:         {
    // INCHI✔❌:             for (i = 0; i < NUM_H_ISOTOPES; i++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 num_iso_H += inp_num_iso_H[i];
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         /*  should not happen because atom here is not aliased */
    // INCHI✔❌:         if (num_iso_H)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (num_H >= num_iso_H)
    // INCHI✔❌:             {
    // INCHI✔❌:                 num_H -= num_iso_H;
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 num_H = inp_num_H; /*  as requested in the alias */
    // INCHI✔❌:                 /* num_H = (num_iso_H - num_H) % 2; */ /*  keep unchanged parity of the total number of H atoms */
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         /*  should not happen because atom here is not aliased */
    // INCHI✔❌:         if (inp_num_H > num_H)
    // INCHI✔❌:         {
    // INCHI✔❌:             num_H = inp_num_H;  /*  as requested in the alias */
    // INCHI✔❌:             /* num_H = inp_num_H + (inp_num_H - num_H)%2; */ /*  keep unchanged parity of the number of non-isotopic H atoms */
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:     else
    // INCHI✔❌:     {
    // INCHI✔❌:         num_H = inp_num_H;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return num_H;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: get_num_H

    let nitrogen = el_number_in_internal_ref_table(Some(&[b'N' as i8, 0]))?;
    let sulfur = el_number_in_internal_ref_table(Some(&[b'S' as i8, 0]))?;
    let _oxygen = el_number_in_internal_ref_table(Some(&[b'O' as i8, 0]))?;
    let carbon = el_number_in_internal_ref_table(Some(&[b'C' as i8, 0]))?;

    if aliased != 0 {
        return Ok(input_hydrogens);
    }
    if atom_input_valence != 0 && (atom_input_valence != 15 || chemical_bond_valence != 0) {
        return Ok(0.max(atom_input_valence - chemical_bond_valence));
    }
    if atom_input_valence == 15 && chemical_bond_valence == 0 {
        return Ok(0);
    }
    if !(MIN_ATOM_CHARGE <= charge && charge <= MAX_ATOM_CHARGE as i32) {
        return Ok(input_hydrogens);
    }
    let element = el_number_in_internal_ref_table(element_name)?;
    if element == ERR_ELEM || EL_DATA_SKIP_ADDING_H[element as usize] || do_not_add_hydrogen != 0 {
        return Ok(input_hydrogens);
    }

    let charge_index = (NEUTRAL_STATE as i32 + charge) as usize;
    let mut hydrogen_count = 0;
    if radical != 0 && radical != RADICAL_SINGLET as i32 {
        let mut valence = EL_DATA_VALENCES[element as usize][charge_index][0] as i32;
        if valence != 0 {
            valence -= if radical == RADICAL_DOUBLET as i32 {
                1
            } else if radical == RADICAL_SINGLET as i32 || radical == RADICAL_TRIPLET as i32 {
                2
            } else {
                valence
            };
            hydrogen_count = 0.max(valence - chemical_bond_valence);
        }
    } else {
        let valences = &EL_DATA_VALENCES[element as usize][charge_index];
        let mut index = 0;
        let mut valence = valences[index] as i32;
        while valence != 0 && valence < chemical_bond_valence {
            index += 1;
            valence = valences[index] as i32;
        }
        if element == nitrogen && charge == 0 && radical == 0 && valence == 5 {
            valence = 3;
        } else if element == sulfur
            && charge == 0
            && radical == 0
            && valence == 4
            && chemical_bond_valence == 3
        {
            valence = 3;
        } else if has_metal_neighbor != 0 && element != carbon && valence > 0 {
            valence -= 1;
        }
        hydrogen_count = 0.max(valence - chemical_bond_valence);
    }

    let mut isotope_hydrogens = 0_i32;
    if let Some(input) = input_isotope_hydrogens {
        for index in 0..NUM_H_ISOTOPES as usize {
            isotope_hydrogens = isotope_hydrogens.wrapping_add(input[index] as i32);
        }
    }
    if isotope_hydrogens != 0 {
        if hydrogen_count >= isotope_hydrogens {
            hydrogen_count -= isotope_hydrogens;
        } else {
            hydrogen_count = input_hydrogens;
        }
    }
    if input_hydrogens > hydrogen_count {
        hydrogen_count = input_hydrogens;
    }
    Ok(hydrogen_count)
}

pub(crate) fn num_of_H(atoms: &[inp_ATOM], iat: i32) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1129 num_of_H
    // INCHI✔️✔️: int num_of_H( inp_ATOM *at, int iat )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     static int el_number_H = (int)EL_NUMBER_H;
    // INCHI✔️✔️:     int    i, n, num_explicit_H = 0;
    // INCHI✔️✔️:     inp_ATOM *a = at + iat;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i = 0; i < a->valence; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         n = a->neighbor[i];
    // INCHI✔️✔️:         num_explicit_H += ( 1 == at[n].valence && el_number_H == at[n].el_number );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return num_explicit_H + NUMH( at, iat );
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: num_of_H
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: num_of_H
    // INCHI✔️✔️: #define EL_NUMBER_H  ((U_CHAR)1)
    // INCHI✔️✔️: #define NUM_ISO_H(AT, N) (AT[N].num_iso_H[0] + AT[N].num_iso_H[1] + AT[N].num_iso_H[2])
    // INCHI✔️✔️: #define NUMH(AT, N) (AT[N].num_H + NUM_ISO_H(AT, N))
    // END INCHI ACTIVE MACRO CONFIGURATION: num_of_H

    let atom_index = usize::try_from(iat).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut explicit_h = 0_i32;
    for index in 0..i32::from(atom.valence) {
        let neighbor_index =
            usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let neighbor_atom_index = usize::from(
            *atom
                .neighbor
                .get(neighbor_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let neighbor_atom = atoms
            .get(neighbor_atom_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        explicit_h += i32::from(neighbor_atom.valence == 1 && neighbor_atom.el_number == 1);
    }
    Ok(explicit_h
        + i32::from(atom.num_H)
        + i32::from(atom.num_iso_H[0])
        + i32::from(atom.num_iso_H[1])
        + i32::from(atom.num_iso_H[2]))
}

pub(crate) fn ion_el_group(el: i32) -> u8 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1151 ion_el_group
    // INCHI✔️✔️: U_CHAR ion_el_group( int el )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     switch ( el )
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         case EL_NUMBER_C: /* fallthrough */
    // INCHI❌❌: #if ( FIX_REM_ION_PAIRS_Si_BUG == 1 )
    // INCHI❌❌:         case EL_NUMBER_SI:
    // INCHI❌❌: #endif
    // INCHI✔️✔️:             return EL_NUMBER_C;
    // INCHI✔️✔️:         case EL_NUMBER_N: /* fallthrough */
    // INCHI✔️✔️:         case EL_NUMBER_P:
    // INCHI✔️✔️:         case EL_NUMBER_AS:
    // INCHI✔️✔️:         case EL_NUMBER_SB:
    // INCHI✔️✔️:             return EL_NUMBER_N;
    // INCHI✔️✔️:         case EL_NUMBER_O: /* fallthrough */
    // INCHI✔️✔️:         case EL_NUMBER_S:
    // INCHI✔️✔️:         case EL_NUMBER_SE:
    // INCHI✔️✔️:         case EL_NUMBER_TE:
    // INCHI✔️✔️:             return EL_NUMBER_O;
    // INCHI✔️✔️:         default:
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: ion_el_group
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: ion_el_group
    // INCHI✔️✔️: #define EL_NUMBER_C  ((U_CHAR)6)
    // INCHI✔️✔️: #define EL_NUMBER_N  ((U_CHAR)7)
    // INCHI✔️✔️: #define EL_NUMBER_O  ((U_CHAR)8)
    // INCHI✔️✔️: #define EL_NUMBER_P  ((U_CHAR)15)
    // INCHI✔️✔️: #define EL_NUMBER_S  ((U_CHAR)16)
    // INCHI✔️✔️: #define EL_NUMBER_AS ((U_CHAR)33)
    // INCHI✔️✔️: #define EL_NUMBER_SE ((U_CHAR)34)
    // INCHI✔️✔️: #define EL_NUMBER_SB ((U_CHAR)51)
    // INCHI✔️✔️: #define EL_NUMBER_TE ((U_CHAR)52)
    // INCHI✔️✔️: #define FIX_REM_ION_PAIRS_Si_BUG 0
    // END INCHI ACTIVE MACRO CONFIGURATION: ion_el_group

    match el {
        6 => 6,
        7 | 15 | 33 | 51 => 7,
        8 | 16 | 34 | 52 => 8,
        _ => 0,
    }
}

pub(crate) fn has_other_ion_neigh(
    atoms: &[inp_ATOM],
    atom_number: i32,
    ion_neighbor_number: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1175 has_other_ion_neigh
    // INCHI✔️✔️: int has_other_ion_neigh( inp_ATOM *at,
    // INCHI✔️✔️:                          int iat,
    // INCHI✔️✔️:                          int iat_ion_neigh)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int charge = at[iat_ion_neigh].charge;
    // INCHI✔️✔️:     int i, neigh;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i = 0; i < at[iat].valence; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         neigh = at[iat].neighbor[i];
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (neigh != iat_ion_neigh && at[neigh].charge == charge &&
    // INCHI✔️✔️:             ion_el_group( at[neigh].el_number ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return 1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: has_other_ion_neigh
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: has_other_ion_neigh
    // INCHI✔️✔️: typedef unsigned short AT_NUMB;
    // INCHI✔️✔️: typedef signed char S_CHAR;
    // INCHI✔️✔️: U_CHAR ion_el_group( int el );
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: has_other_ion_neigh

    let ion_neighbor = atoms
        .get(
            usize::try_from(ion_neighbor_number)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
        )
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let charge = ion_neighbor.charge;
    let atom = atoms
        .get(usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for ordinal in 0..i32::from(atom.valence) {
        let neighbor_number = i32::from(
            *atom
                .neighbor
                .get(usize::try_from(ordinal).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let neighbor = atoms
            .get(
                usize::try_from(neighbor_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
            )
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if neighbor_number != ion_neighbor_number
            && neighbor.charge == charge
            && ion_el_group(i32::from(neighbor.el_number)) != 0
        {
            return Ok(1);
        }
    }
    Ok(0)
}

pub(crate) fn has_other_ion_in_sphere_2(
    atoms: &mut [inp_ATOM],
    atom_number: i32,
    ion_neighbor_number: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1201 has_other_ion_in_sphere_2
    // INCHI✔️✔️: int has_other_ion_in_sphere_2( inp_ATOM *at, int iat,
    // INCHI✔️✔️:                                int iat_ion_neigh )
    // INCHI✔️✔️: {
    // INCHI✔️✔️: #define MAXQ 16
    // INCHI✔️✔️:     AT_NUMB q[MAXQ];
    // INCHI✔️✔️:     int lenq = 0, lenq2, dist = 0, i = 0, iq, neigh, j, nRet = 0;
    // INCHI✔️✔️:     q[lenq++] = iat;
    // INCHI✔️✔️:     at[iat].cFlags = 1;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     iq = 0;
    // INCHI✔️✔️:     dist = 1;
    // INCHI✔️✔️:     /* use at->cFlags as an indicator */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     while (dist <= 2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         for (lenq2 = lenq; iq < lenq2; iq++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             i = q[iq];
    // INCHI✔️✔️:
    // INCHI✔️✔️:             for (j = 0; j < at[i].valence; j++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 neigh = at[i].neighbor[j];
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 if (!at[neigh].cFlags &&
    // INCHI✔️✔️:                      at[neigh].valence <= 3 &&
    // INCHI✔️✔️:                      ion_el_group( at[neigh].el_number ))
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     q[lenq++] = neigh;
    // INCHI✔️✔️:                     at[neigh].cFlags = 1;
    // INCHI✔️✔️:                     if (neigh != iat_ion_neigh &&
    // INCHI✔️✔️:                          at[iat_ion_neigh].charge == at[neigh].charge)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         nRet++;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         dist++;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (iq = 0; iq < lenq; iq++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         i = q[iq];
    // INCHI✔️✔️:         at[i].cFlags = 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return nRet;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: has_other_ion_in_sphere_2
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: has_other_ion_in_sphere_2
    // INCHI✔️✔️: #define MAXQ 16
    // INCHI✔️✔️: typedef unsigned short AT_NUMB;
    // INCHI✔️✔️: typedef signed char S_CHAR;
    // INCHI✔️✔️: U_CHAR ion_el_group( int el );
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: has_other_ion_in_sphere_2

    const MAXQ: usize = 16;
    let center_index =
        usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    atoms
        .get(center_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;

    let mut queue = [0_u16; MAXQ];
    let mut queue_length = 0_usize;
    queue[queue_length] = atom_number as u16;
    queue_length += 1;
    atoms[center_index].cFlags = 1;

    let mut queue_index = 0_usize;
    let mut distance = 1_i32;
    let mut result = 0_i32;
    while distance <= 2 {
        let layer_end = queue_length;
        while queue_index < layer_end {
            let current_index = usize::from(queue[queue_index]);
            let valence = i32::from(
                atoms
                    .get(current_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .valence,
            );
            for ordinal in 0..valence {
                let neighbor_index = usize::from(
                    *atoms
                        .get(current_index)
                        .and_then(|atom| atom.neighbor.get(usize::try_from(ordinal).ok()?))
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let neighbor = atoms
                    .get(neighbor_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if neighbor.cFlags == 0
                    && neighbor.valence <= 3
                    && ion_el_group(i32::from(neighbor.el_number)) != 0
                {
                    let queue_slot = queue
                        .get_mut(queue_length)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    *queue_slot = neighbor_index as u16;
                    queue_length += 1;
                    atoms[neighbor_index].cFlags = 1;
                    if neighbor_index as i32 != ion_neighbor_number {
                        let ion_charge = atoms
                            .get(
                                usize::try_from(ion_neighbor_number)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .charge;
                        if ion_charge == atoms[neighbor_index].charge {
                            result += 1;
                        }
                    }
                }
            }
            queue_index += 1;
        }
        distance += 1;
    }

    for &queued_atom in &queue[..queue_length] {
        atoms[usize::from(queued_atom)].cFlags = 0;
    }
    Ok(result)
}

pub(crate) fn get_endpoint_valence(el_number: u8) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1509 get_endpoint_valence
    // INCHI✔️✔️: int get_endpoint_valence( U_CHAR el_number )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     switch (el_number)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         case EL_NUMBER_O:  /* fallthrough */
    // INCHI✔️✔️:         case EL_NUMBER_S:
    // INCHI✔️✔️:         case EL_NUMBER_SE:
    // INCHI✔️✔️:         case EL_NUMBER_TE:
    // INCHI✔️✔️:             return 2;
    // INCHI✔️✔️:         case EL_NUMBER_N:
    // INCHI✔️✔️:             return 3;
    // INCHI✔️✔️:         default:
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: get_endpoint_valence
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: get_endpoint_valence
    // INCHI✔️✔️: #define EL_NUMBER_N  ((U_CHAR)7)
    // INCHI✔️✔️: #define EL_NUMBER_O  ((U_CHAR)8)
    // INCHI✔️✔️: #define EL_NUMBER_S  ((U_CHAR)16)
    // INCHI✔️✔️: #define EL_NUMBER_SE ((U_CHAR)34)
    // INCHI✔️✔️: #define EL_NUMBER_TE ((U_CHAR)52)
    // END INCHI ACTIVE MACRO CONFIGURATION: get_endpoint_valence

    match el_number {
        8 | 16 | 34 | 52 => 2,
        7 => 3,
        _ => 0,
    }
}

#[allow(non_snake_case)]
pub(crate) fn get_endpoint_valence_KET(el_number: u8) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1530 get_endpoint_valence_KET
    // INCHI✔️✔️: int get_endpoint_valence_KET( U_CHAR el_number )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     switch (el_number)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         case EL_NUMBER_C:
    // INCHI✔️✔️:             return 4;
    // INCHI✔️✔️:         case EL_NUMBER_O:
    // INCHI✔️✔️:             return 2;
    // INCHI✔️✔️:         default:
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: get_endpoint_valence_KET
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: get_endpoint_valence_KET
    // INCHI✔️✔️: #define KETO_ENOL_TAUT             1 /* include keto-enol tautomerism */
    // INCHI✔️✔️: #define EL_NUMBER_C  ((U_CHAR)6)
    // INCHI✔️✔️: #define EL_NUMBER_O  ((U_CHAR)8)
    // END INCHI ACTIVE MACRO CONFIGURATION: get_endpoint_valence_KET

    match el_number {
        6 => 4,
        8 => 2,
        _ => 0,
    }
}

pub(crate) fn get_atomic_mass_from_elnum(atomic_number: i32) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1007 get_atomic_mass_from_elnum
    // INCHI✔❌: int get_atomic_mass_from_elnum( int nAtNum )
    // INCHI✔❌: {
    // INCHI✔❌:     nAtNum -= 1;
    // INCHI✔❌:
    // INCHI✔❌:     if (0 < nAtNum)
    // INCHI✔❌:     {
    // INCHI✔❌:         nAtNum += 2; /*  bypass D, T */
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (0 <= nAtNum && nAtNum < nElDataLen)
    // INCHI✔❌:     {
    // INCHI✔❌:         return (int) ElData[nAtNum].nAtMass;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return 0;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: get_atomic_mass_from_elnum

    let mut table_index = atomic_number
        .checked_sub(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    if table_index > 0 {
        table_index = table_index
            .checked_add(2)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    }
    if table_index >= 0 && table_index < EL_DATA_ATOMIC_MASSES.len() as i32 {
        return Ok(EL_DATA_ATOMIC_MASSES[table_index as usize]);
    }
    Ok(0)
}

#[rustfmt::skip]
pub(crate) fn get_atomic_mass(element_name: Option<&[i8]>) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1040 get_atomic_mass
    // INCHI✔️❌: complete source frame follows verbatim; Rust slice validation and Result propagation add overhead.
    /*
int get_atomic_mass( const char *elname )
{
    int el_number, atw;
    if (ERR_ELEM != ( el_number = el_number_in_internal_ref_table( elname ) ))
    {
        atw = ElData[el_number].nAtMass;
    }
    else
    {
        atw = 0;
    }

    return atw;
}
    */
    // END INCHI C FUNCTION: get_atomic_mass
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: get_atomic_mass
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; inactive INCHI_ZFRAG ElData rows are excluded.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: get_atomic_mass

    let element_number = el_number_in_internal_ref_table(element_name)?;
    if element_number == ERR_ELEM {
        return Ok(0);
    }
    EL_DATA_ATOMIC_MASSES
        .get(usize::try_from(element_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .copied()
        .ok_or(SourceHeapError::PointerOutOfBounds)
}

pub(crate) fn is_in_the_list(
    path_atom: Option<&[AT_NUMB]>,
    next_atom: AT_NUMB,
    path_length: i32,
) -> Result<Option<usize>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1059 is_in_the_list
    // INCHI✔❌: AT_NUMB *is_in_the_list( AT_NUMB *pathAtom, AT_NUMB nNextAtom, int nPathLen )
    // INCHI✔❌: {
    // INCHI✔❌:     for (; nPathLen && *pathAtom != nNextAtom; nPathLen--, pathAtom++)
    // INCHI✔❌:     {
    // INCHI✔❌:         ;
    // INCHI✔❌:     }
    // INCHI✔❌:     return nPathLen ? pathAtom : NULL;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: is_in_the_list

    if path_length < 0 {
        return Err(SourceHeapError::SourceIntegerOverflow);
    }
    if path_length == 0 {
        return Ok(None);
    }
    let path_length =
        usize::try_from(path_length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let path = path_atom.ok_or(SourceHeapError::NullPointer)?;
    for index in 0..path_length {
        let candidate = path.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        if *candidate == next_atom {
            return Ok(Some(index));
        }
    }
    Ok(None)
}

pub(crate) fn is_in_the_ilist(
    path_atom: Option<&[i32]>,
    next_atom: i32,
    path_length: i32,
) -> Result<Option<usize>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1072 is_in_the_ilist
    // INCHI✔️✔️: int *is_in_the_ilist( int *pathAtom, int nNextAtom, int nPathLen )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     for (; nPathLen && *pathAtom != nNextAtom; nPathLen--, pathAtom++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     return nPathLen ? pathAtom : NULL;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: is_in_the_ilist

    if path_length < 0 {
        return Err(SourceHeapError::SourceIntegerOverflow);
    }
    if path_length == 0 {
        return Ok(None);
    }
    let path_length =
        usize::try_from(path_length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let path = path_atom.ok_or(SourceHeapError::NullPointer)?;
    for index in 0..path_length {
        let candidate = path.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
        if *candidate == next_atom {
            return Ok(Some(index));
        }
    }
    Ok(None)
}

pub(crate) fn is_ilist_inside(
    integer_list: Option<&[i32]>,
    list_length: i32,
    embedding_list: Option<&[i32]>,
    embedding_length: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1085 is_ilist_inside
    // INCHI✔️✔️: int is_ilist_inside( int *ilist, int nlist, int *ilist2, int nlist2 )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int k;
    // INCHI✔️✔️:     for (k = 0; k < nlist; k++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (!is_in_the_ilist( ilist2, ilist[k], nlist2 ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     return 1;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: is_ilist_inside
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: is_ilist_inside
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; no conditional source branches.
    // INCHI✔️✔️: The direct is_in_the_ilist callee is invoked once per source item until failure.
    // END INCHI ACTIVE MACRO CONFIGURATION: is_ilist_inside

    if list_length <= 0 {
        return Ok(1);
    }
    let list_length =
        usize::try_from(list_length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let integer_list = integer_list.ok_or(SourceHeapError::NullPointer)?;
    for index in 0..list_length {
        let value = *integer_list
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if is_in_the_ilist(embedding_list, value, embedding_length)?.is_none() {
            return Ok(0);
        }
    }
    Ok(1)
}

pub(crate) fn n_bonds_val_to_metal(
    atoms: Option<&[inp_ATOM]>,
    atom_index: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1100 nBondsValToMetal
    // INCHI✔❌: int nBondsValToMetal( inp_ATOM* at, int iat )
    // INCHI✔❌: {
    // INCHI✔❌:     int i, bond_type, nVal2Metal = 0; /* djb-rwth: removing redundant variables */
    // INCHI✔❌:     inp_ATOM* a = at + iat;
    // INCHI✔❌:
    // INCHI✔❌:     for (i = 0; i < a->valence; i++)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* djb-rwth: removing redundant code */
    // INCHI✔❌:
    // INCHI✔❌:         if (is_el_a_metal( at[(int) a->neighbor[i]].el_number ))
    // INCHI✔❌:         {
    // INCHI✔❌:             bond_type = a->bond_type[i];
    // INCHI✔❌:
    // INCHI✔❌:             if (bond_type <= BOND_TYPE_TRIPLE)
    // INCHI✔❌:             {
    // INCHI✔❌:                 nVal2Metal += bond_type;
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 return -1;  /* bond to metal order is not well defined */
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return nVal2Metal;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: nBondsValToMetal

    let atoms = atoms.ok_or(SourceHeapError::NullPointer)?;
    let atom_index =
        usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut valence_to_metal = 0_i32;

    for index in 0..i32::from(atom.valence) {
        let index = usize::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let neighbor_index = *atom
            .neighbor
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let neighbor = atoms
            .get(usize::from(neighbor_index))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if is_el_a_metal(i32::from(neighbor.el_number))? != 0 {
            let bond_type = i32::from(
                *atom
                    .bond_type
                    .get(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            if bond_type <= BOND_TYPE_TRIPLE as i32 {
                valence_to_metal = valence_to_metal
                    .checked_add(bond_type)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            } else {
                return Ok(-1);
            }
        }
    }

    Ok(valence_to_metal)
}

pub(crate) fn n_no_metal_num_bonds(
    atoms: Option<&[inp_ATOM]>,
    atom_index: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1253 nNoMetalNumBonds
    // INCHI✔️✔️: int nNoMetalNumBonds( inp_ATOM *at, int at_no )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     inp_ATOM *a = at + at_no;
    // INCHI✔️✔️:     int num_H = NUMH( a, 0 );
    // INCHI✔️✔️:     int std_chem_bonds_valence = get_el_valence( a->el_number, a->charge, 0 );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (a->chem_bonds_valence + num_H > std_chem_bonds_valence)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         int valence_to_metal = 0;
    // INCHI✔️✔️:         int num_bonds_to_metal = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         for (i = 0; i < a->valence; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (is_el_a_metal( at[(int) a->neighbor[i]].el_number ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (( a->bond_type[i] & BOND_TYPE_MASK ) >= BOND_TYPE_ALTERN)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     return a->valence; /* fall back */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 num_bonds_to_metal++;
    // INCHI✔️✔️:                 valence_to_metal += ( a->bond_type[i] & BOND_TYPE_MASK );
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (a->chem_bonds_valence + num_H - valence_to_metal == std_chem_bonds_valence)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* removing bonds to metal produces standard valence */
    // INCHI✔️✔️:             return a->valence - num_bonds_to_metal;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️: #if ( S_VI_O_PLUS_METAL_FIX_BOND == 1 )
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (1 == a->charge && 2 == get_endpoint_valence( a->el_number ) &&
    // INCHI✔️✔️:                 a->chem_bonds_valence + num_H == std_chem_bonds_valence)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             int valence_to_metal = 0;
    // INCHI✔️✔️:             int num_bonds_to_metal = 0;
    // INCHI✔️✔️:             for (i = 0; i < a->valence; i++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (is_el_a_metal( at[(int) a->neighbor[i]].el_number ))
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     if (( a->bond_type[i] & BOND_TYPE_MASK ) >= BOND_TYPE_ALTERN)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         return a->valence; /* fall back */
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     num_bonds_to_metal++;
    // INCHI✔️✔️:                     valence_to_metal += ( a->bond_type[i] & BOND_TYPE_MASK );
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (1 == valence_to_metal)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 /* removing bonds to metal produces standard valence */
    // INCHI✔️✔️:                 return a->valence - num_bonds_to_metal;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return a->valence;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: nNoMetalNumBonds
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: nNoMetalNumBonds
    // INCHI✔️✔️: #define BOND_TYPE_ALTERN 4
    // INCHI✔️✔️: #define BOND_TYPE_MASK 15
    // INCHI✔️✔️: #define S_VI_O_PLUS_METAL_FIX_BOND 1
    // END INCHI ACTIVE MACRO CONFIGURATION: nNoMetalNumBonds

    let atoms = atoms.ok_or(SourceHeapError::NullPointer)?;
    let atom_index =
        usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let num_h = i32::from(atom.num_H)
        + i32::from(atom.num_iso_H[0])
        + i32::from(atom.num_iso_H[1])
        + i32::from(atom.num_iso_H[2]);
    let std_chem_bonds_valence =
        get_el_valence(i32::from(atom.el_number), i32::from(atom.charge), 0)?;

    if i32::from(atom.chem_bonds_valence) + num_h > std_chem_bonds_valence {
        let mut valence_to_metal = 0_i32;
        let mut num_bonds_to_metal = 0_i32;
        for index in 0..i32::from(atom.valence) {
            let index =
                usize::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let neighbor_index = *atom
                .neighbor
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let neighbor = atoms
                .get(usize::from(neighbor_index))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if is_el_a_metal(i32::from(neighbor.el_number))? != 0 {
                let bond_type = i32::from(
                    *atom
                        .bond_type
                        .get(index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                if (bond_type & BOND_TYPE_MASK as i32) >= BOND_TYPE_ALTERN as i32 {
                    return Ok(i32::from(atom.valence));
                }
                num_bonds_to_metal += 1;
                valence_to_metal += bond_type & BOND_TYPE_MASK as i32;
            }
        }
        if i32::from(atom.chem_bonds_valence) + num_h - valence_to_metal == std_chem_bonds_valence {
            return Ok(i32::from(atom.valence) - num_bonds_to_metal);
        }
    } else if i32::from(atom.charge) == 1
        && get_endpoint_valence(atom.el_number) == 2
        && i32::from(atom.chem_bonds_valence) + num_h == std_chem_bonds_valence
    {
        let mut valence_to_metal = 0_i32;
        let mut num_bonds_to_metal = 0_i32;
        for index in 0..i32::from(atom.valence) {
            let index =
                usize::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let neighbor_index = *atom
                .neighbor
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let neighbor = atoms
                .get(usize::from(neighbor_index))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if is_el_a_metal(i32::from(neighbor.el_number))? != 0 {
                let bond_type = i32::from(
                    *atom
                        .bond_type
                        .get(index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                if (bond_type & BOND_TYPE_MASK as i32) >= BOND_TYPE_ALTERN as i32 {
                    return Ok(i32::from(atom.valence));
                }
                num_bonds_to_metal += 1;
                valence_to_metal += bond_type & BOND_TYPE_MASK as i32;
            }
        }
        if valence_to_metal == 1 {
            return Ok(i32::from(atom.valence) - num_bonds_to_metal);
        }
    }

    Ok(i32::from(atom.valence))
}

pub(crate) fn n_no_metal_bonds_valence(
    atoms: Option<&[inp_ATOM]>,
    atom_index: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1320 nNoMetalBondsValence
    // INCHI✔️✔️: int nNoMetalBondsValence( inp_ATOM *at, int at_no )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     inp_ATOM *a = at + at_no;
    // INCHI✔️✔️:     int num_H = NUMH( a, 0 );
    // INCHI✔️✔️:     int std_chem_bonds_valence = get_el_valence( a->el_number, a->charge, 0 );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (a->chem_bonds_valence + num_H > std_chem_bonds_valence)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         int valence_to_metal = 0;
    // INCHI✔️✔️:         /*int num_bonds_to_metal = 0;*/
    // INCHI✔️✔️:
    // INCHI✔️✔️:         for (i = 0; i < a->valence; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (is_el_a_metal( at[(int) a->neighbor[i]].el_number ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (( a->bond_type[i] & BOND_TYPE_MASK ) >= BOND_TYPE_ALTERN)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     return a->valence; /* fall back */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 /* num_bonds_to_metal ++;*/
    // INCHI✔️✔️:                 valence_to_metal += ( a->bond_type[i] & BOND_TYPE_MASK );
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (a->chem_bonds_valence + num_H - valence_to_metal == std_chem_bonds_valence)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* removing bonds to metal produces standard valence */
    // INCHI✔️✔️:             return a->chem_bonds_valence - valence_to_metal;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️: #if ( S_VI_O_PLUS_METAL_FIX_BOND == 1 )
    // INCHI✔️✔️:     else if (1 == a->charge && 2 == get_endpoint_valence( a->el_number ) &&
    // INCHI✔️✔️:          a->chem_bonds_valence + num_H == std_chem_bonds_valence)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:
    // INCHI✔️✔️:         int valence_to_metal = 0;
    // INCHI✔️✔️:         /* int num_bonds_to_metal = 0;*/
    // INCHI✔️✔️:
    // INCHI✔️✔️:         for (i = 0; i < a->valence; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (is_el_a_metal( at[(int) a->neighbor[i]].el_number ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (( a->bond_type[i] & BOND_TYPE_MASK ) >= BOND_TYPE_ALTERN)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     return a->valence; /* fall back */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 /* num_bonds_to_metal ++;*/
    // INCHI✔️✔️:                 valence_to_metal += ( a->bond_type[i] & BOND_TYPE_MASK );
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (1 == valence_to_metal)
    // INCHI✔️✔️:         {/* removing bonds to metal produces standard valence */
    // INCHI✔️✔️:             return a->chem_bonds_valence - valence_to_metal;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return a->chem_bonds_valence;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: nNoMetalBondsValence
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: nNoMetalBondsValence
    // INCHI✔️✔️: #define BOND_TYPE_ALTERN 4
    // INCHI✔️✔️: #define BOND_TYPE_MASK 15
    // INCHI✔️✔️: #define S_VI_O_PLUS_METAL_FIX_BOND 1
    // END INCHI ACTIVE MACRO CONFIGURATION: nNoMetalBondsValence

    let atoms = atoms.ok_or(SourceHeapError::NullPointer)?;
    let atom_index =
        usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let num_h = i32::from(atom.num_H)
        + i32::from(atom.num_iso_H[0])
        + i32::from(atom.num_iso_H[1])
        + i32::from(atom.num_iso_H[2]);
    let atom_chem_valence = i32::from(atom.chem_bonds_valence);
    let std_chem_bonds_valence =
        get_el_valence(i32::from(atom.el_number), i32::from(atom.charge), 0)?;

    if atom_chem_valence + num_h > std_chem_bonds_valence {
        let mut valence_to_metal = 0_i32;
        for index in 0..i32::from(atom.valence) {
            let index =
                usize::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let neighbor_index = *atom
                .neighbor
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let neighbor = atoms
                .get(usize::from(neighbor_index))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if is_el_a_metal(i32::from(neighbor.el_number))? != 0 {
                let bond_type = i32::from(
                    *atom
                        .bond_type
                        .get(index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                if (bond_type & BOND_TYPE_MASK as i32) >= BOND_TYPE_ALTERN as i32 {
                    return Ok(i32::from(atom.valence));
                }
                valence_to_metal += bond_type & BOND_TYPE_MASK as i32;
            }
        }
        if atom_chem_valence + num_h - valence_to_metal == std_chem_bonds_valence {
            return Ok(atom_chem_valence - valence_to_metal);
        }
    } else if i32::from(atom.charge) == 1
        && get_endpoint_valence(atom.el_number) == 2
        && atom_chem_valence + num_h == std_chem_bonds_valence
    {
        let mut valence_to_metal = 0_i32;
        for index in 0..i32::from(atom.valence) {
            let index =
                usize::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let neighbor_index = *atom
                .neighbor
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let neighbor = atoms
                .get(usize::from(neighbor_index))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if is_el_a_metal(i32::from(neighbor.el_number))? != 0 {
                let bond_type = i32::from(
                    *atom
                        .bond_type
                        .get(index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                if (bond_type & BOND_TYPE_MASK as i32) >= BOND_TYPE_ALTERN as i32 {
                    return Ok(i32::from(atom.valence));
                }
                valence_to_metal += bond_type & BOND_TYPE_MASK as i32;
            }
        }
        if valence_to_metal == 1 {
            return Ok(atom_chem_valence - valence_to_metal);
        }
    }

    Ok(atom_chem_valence)
}

pub(crate) fn n_no_metal_neigh_index(
    atoms: Option<&[inp_ATOM]>,
    atom_index: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1386 nNoMetalNeighIndex
    // INCHI✔️✔️: int nNoMetalNeighIndex( inp_ATOM *at, int at_no )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     inp_ATOM *a = at + at_no;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i = 0; i < a->valence; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (!is_el_a_metal( at[(int) a->neighbor[i]].el_number ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return i;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return -1;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: nNoMetalNeighIndex

    let atoms = atoms.ok_or(SourceHeapError::NullPointer)?;
    let atom_index =
        usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for index in 0..i32::from(atom.valence) {
        let index_usize =
            usize::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let neighbor_index = *atom
            .neighbor
            .get(index_usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let neighbor = atoms
            .get(usize::from(neighbor_index))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if is_el_a_metal(i32::from(neighbor.el_number))? == 0 {
            return Ok(index);
        }
    }
    Ok(-1)
}

pub(crate) fn n_no_metal_other_neigh_index(
    atoms: Option<&[inp_ATOM]>,
    atom_index: i32,
    current_neigh: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1405 nNoMetalOtherNeighIndex
    // INCHI✔️✔️: int nNoMetalOtherNeighIndex( inp_ATOM *at, int at_no, int cur_neigh )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, neigh;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     inp_ATOM *a = at + at_no;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i = 0; i < a->valence; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         neigh = (int) a->neighbor[i];
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (neigh != cur_neigh && !is_el_a_metal( at[neigh].el_number ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return i;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return -1;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: nNoMetalOtherNeighIndex

    let atoms = atoms.ok_or(SourceHeapError::NullPointer)?;
    let atom_index =
        usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for index in 0..i32::from(atom.valence) {
        let index_usize =
            usize::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let neigh = i32::from(
            *atom
                .neighbor
                .get(index_usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        if neigh != current_neigh {
            let neighbor = atoms
                .get(usize::try_from(neigh).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if is_el_a_metal(i32::from(neighbor.el_number))? == 0 {
                return Ok(index);
            }
        }
    }
    Ok(-1)
}

pub(crate) fn n_no_metal_other_neigh_index2(
    atoms: Option<&[inp_ATOM]>,
    atom_index: i32,
    current_neigh: i32,
    current_neigh2: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1426 nNoMetalOtherNeighIndex2
    // INCHI✔️✔️: int nNoMetalOtherNeighIndex2( inp_ATOM *at,
    // INCHI✔️✔️:                               int at_no,
    // INCHI✔️✔️:                               int cur_neigh,
    // INCHI✔️✔️:                               int cur_neigh2 )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, neigh;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     inp_ATOM *a = at + at_no;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     for (i = 0; i < a->valence; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         neigh = (int) a->neighbor[i];
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (neigh != cur_neigh && neigh != cur_neigh2 && !is_el_a_metal( at[neigh].el_number ))
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return i;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return -1;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: nNoMetalOtherNeighIndex2

    let atoms = atoms.ok_or(SourceHeapError::NullPointer)?;
    let atom_index =
        usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = atoms
        .get(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for index in 0..i32::from(atom.valence) {
        let index_usize =
            usize::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let neigh = i32::from(
            *atom
                .neighbor
                .get(index_usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        if neigh != current_neigh && neigh != current_neigh2 {
            let neighbor = atoms
                .get(usize::try_from(neigh).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if is_el_a_metal(i32::from(neighbor.el_number))? == 0 {
                return Ok(index);
            }
        }
    }
    Ok(-1)
}

pub(crate) fn inchi_memicmp(
    heap: &SourceHeap,
    first: SourceConstPointer<i8>,
    second: SourceConstPointer<i8>,
    mut length: u64,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1971 inchi_memicmp
    // INCHI✔️❌: int inchi_memicmp( const void * p1, const void * p2, size_t length )
    // INCHI✔️❌: {
    // INCHI✔️❌:     const U_CHAR *s1 = (const U_CHAR*) p1;
    // INCHI✔️❌:     const U_CHAR *s2 = (const U_CHAR*) p2;
    // INCHI✔️❌:     while (length--)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (*s1 == *s2 ||
    // INCHI✔️❌:               __MYTOLOWER( (int) *s1 ) == __MYTOLOWER( (int) *s2 ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             s1++;
    // INCHI✔️❌:             s2++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return
    // INCHI✔️❌:                 __MYTOLOWER( (int) *s1 ) - __MYTOLOWER( (int) *s2 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_memicmp

    let lower = |byte: u8| {
        if byte.is_ascii_uppercase() {
            byte + (b'a' - b'A')
        } else {
            byte
        }
    };
    let mut index = 0_usize;
    while length != 0 {
        length -= 1;
        let left = *heap
            .slice(first)?
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? as u8;
        let right = *heap
            .slice(second)?
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? as u8;
        if left != right && lower(left) != lower(right) {
            return Ok(i32::from(lower(left)) - i32::from(lower(right)));
        }
        index = index
            .checked_add(1)
            .ok_or(SourceHeapError::PointerOffsetOverflow)?;
    }
    Ok(0)
}

pub(crate) fn inchi_stricmp(
    heap: &SourceHeap,
    first: SourceConstPointer<i8>,
    second: SourceConstPointer<i8>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1995 inchi_stricmp
    // INCHI✔️❌: int inchi_stricmp( const char *s1, const char *s2 )
    // INCHI✔️❌: {
    // INCHI✔️❌:     while (*s1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (*s1 == *s2 ||
    // INCHI✔️❌:               __MYTOLOWER( (int) *s1 ) == __MYTOLOWER( (int) *s2 ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             s1++;
    // INCHI✔️❌:             s2++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return
    // INCHI✔️❌:                 __MYTOLOWER( (int) *s1 ) - __MYTOLOWER( (int) *s2 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (*s2)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_stricmp

    let lower = |byte: i8| {
        let value = i32::from(byte);
        if value >= i32::from(b'A') && value <= i32::from(b'Z') {
            value - i32::from(b'A') + i32::from(b'a')
        } else {
            value
        }
    };
    let first_bytes = heap.slice(first)?;
    let second_bytes = heap.slice(second)?;
    let mut index = 0_usize;
    loop {
        let left = *first_bytes
            .get(index)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        if left == 0 {
            let right = *second_bytes
                .get(index)
                .ok_or(SourceHeapError::MissingNulTerminator)?;
            return Ok(if right != 0 { -1 } else { 0 });
        }
        let right = *second_bytes
            .get(index)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        if left != right && lower(left) != lower(right) {
            return Ok(lower(left) - lower(right));
        }
        index = index
            .checked_add(1)
            .ok_or(SourceHeapError::PointerOffsetOverflow)?;
    }
}

pub(crate) fn inchi__strdup(
    heap: &mut SourceHeap,
    string: SourceConstPointer<i8>,
) -> Result<SourceMutPointer<i8>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:2035 inchi__strdup
    // INCHI✔❌: char *inchi__strdup( const char *string )
    // INCHI✔❌: {
    // INCHI✔❌:     char *p = NULL;
    // INCHI✔❌:     if (string)
    // INCHI✔❌:     {
    // INCHI✔❌:         size_t length = strlen( string );
    // INCHI✔❌:         p = (char *) inchi_malloc( length + 1 );
    // INCHI✔❌:         if (p)
    // INCHI✔❌:         {
    // INCHI✔❌:             strcpy(p, string);
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return p;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: inchi__strdup

    if string.is_null() {
        return Ok(SourceMutPointer::null());
    }
    let copied = {
        let source = heap.slice(string)?;
        let length = source
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        source[..=length].to_vec()
    };
    let byte_count = u64::try_from(copied.len())
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    let duplicate = inchi_malloc(heap, byte_count)?;
    heap.slice_mut(duplicate)?.copy_from_slice(&copied);
    Ok(duplicate)
}

#[rustfmt::skip]
pub(crate) fn dotify_non_printable_chars(
    heap: &mut SourceHeap,
    line: SourceMutPointer<i8>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1630 dotify_non_printable_chars
    // INCHI✔️❌: complete source frame follows verbatim; checked source-heap access adds overhead.
    /*
int dotify_non_printable_chars( char *line )
{
    int i, c, num = 0;

    if (line)
    {
        for (i = 0; (c = UCINT line[i]); i++) /* djb-rwth: addressing LLVM warning */
        {
            /* assuming ASCII charset */
            if (c < ' ' || c >= 0x7F)
            {
                line[i] = '.';
                num++;
            }
        }
    }

    return num;
}
    */
    // END INCHI C FUNCTION: dotify_non_printable_chars
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: dotify_non_printable_chars
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; UCINT casts each signed char through unsigned char.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: dotify_non_printable_chars

    if line.is_null() {
        return Ok(0);
    }
    let bytes = heap.slice_mut(line)?;
    let mut replacements = 0_i32;
    for byte in bytes {
        let value = *byte as u8;
        if value == 0 {
            return Ok(replacements);
        }
        if value < b' ' || value >= 0x7f {
            *byte = b'.' as i8;
            replacements = replacements.wrapping_add(1);
        }
    }
    Err(SourceHeapError::MissingNulTerminator)
}

#[cfg(test)]
#[allow(non_snake_case)]
mod tests {
    use super::*;
    use crate::source_types::{BOND_DOUBLE, BOND_SINGLE};

    #[test]
    fn source_port__util__dotify_non_printable_chars__line_1630() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            dotify_non_printable_chars(&mut heap, SourceMutPointer::null()),
            Ok(0)
        );

        let empty = heap.allocate(vec![0_i8, 7]).unwrap();
        assert_eq!(dotify_non_printable_chars(&mut heap, empty), Ok(0));
        assert_eq!(heap.slice(empty.as_const()).unwrap(), &[0, 7]);

        let line = heap
            .allocate(vec![1_i8, 31, 32, 126, 127, -1, b'A' as i8, 0, 2])
            .unwrap();
        assert_eq!(dotify_non_printable_chars(&mut heap, line), Ok(4));
        assert_eq!(
            heap.slice(line.as_const()).unwrap(),
            &[
                b'.' as i8, b'.' as i8, 32, 126, b'.' as i8, b'.' as i8, b'A' as i8, 0, 2
            ]
        );

        let unterminated = heap.allocate(vec![b'A' as i8]).unwrap();
        assert_eq!(
            dotify_non_printable_chars(&mut heap, unterminated),
            Err(SourceHeapError::MissingNulTerminator)
        );
    }

    #[test]
    fn source_port__util__get_element_or_pseudoelement_symbol__line_316() {
        let expected = [
            "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
            "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
            "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
            "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr",
            "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
            "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
            "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
            "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc",
            "Lv", "Ts", "Og", "Zz", "Zz",
        ];
        assert_eq!(expected.len(), 120);
        for (index, symbol) in expected.iter().enumerate() {
            let mut output = [b'X' as i8; 5];
            assert_eq!(
                get_element_or_pseudoelement_symbol(index as i32 + 1, &mut output),
                Ok(0),
                "atomic number {}",
                index + 1
            );
            let bytes = symbol.as_bytes();
            assert_eq!(
                &output[..bytes.len()],
                &bytes.iter().map(|byte| *byte as i8).collect::<Vec<_>>()
            );
            assert_eq!(output[bytes.len()], 0);
            assert!(
                output[bytes.len() + 1..]
                    .iter()
                    .all(|byte| *byte == b'X' as i8)
            );
        }

        for atomic_number in [i32::MIN, -100, -1, 0, 121, 122, i32::MAX] {
            let mut output = [b'X' as i8; 5];
            assert_eq!(
                get_element_or_pseudoelement_symbol(atomic_number, &mut output),
                Ok(-1),
                "{atomic_number}"
            );
            assert_eq!(output, [b'?' as i8, b'?' as i8, 0, b'X' as i8, b'X' as i8]);
        }

        let mut short = [b'X' as i8; 2];
        assert_eq!(
            get_element_or_pseudoelement_symbol(2, &mut short),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(short, [b'X' as i8; 2]);
    }
    use crate::test_support::{allocate_source_fixture, assert_source_slice_eq};

    #[test]
    fn source_port__util__inchi_memicmp__line_1971() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            inchi_memicmp(
                &heap,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                0,
            ),
            Ok(0)
        );
        let left = allocate_source_fixture(
            &mut heap,
            vec![b'A' as i8, b'z' as i8, 0x80_u8 as i8, b'Q' as i8],
        );
        let equal = allocate_source_fixture(
            &mut heap,
            vec![b'a' as i8, b'Z' as i8, 0x80_u8 as i8, b'q' as i8],
        );
        assert_eq!(
            inchi_memicmp(&heap, left.as_const(), equal.as_const(), 4),
            Ok(0)
        );
        let low = allocate_source_fixture(&mut heap, vec![b'a' as i8, 0]);
        let high = allocate_source_fixture(&mut heap, vec![0xff_u8 as i8, 0]);
        assert_eq!(
            inchi_memicmp(&heap, low.as_const(), high.as_const(), 1),
            Ok(97 - 255)
        );
        assert_eq!(
            inchi_memicmp(&heap, high.as_const(), low.as_const(), 1),
            Ok(255 - 97)
        );
        assert_eq!(
            inchi_memicmp(&heap, left.as_const(), equal.as_const(), 5),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            inchi_memicmp(&heap, SourceConstPointer::null(), equal.as_const(), 1),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__util__inchi_stricmp__line_1995() {
        let mut heap = SourceHeap::default();
        let empty = allocate_source_fixture(&mut heap, vec![0_i8]);
        let a = allocate_source_fixture(&mut heap, vec![b'A' as i8, 0]);
        let lower_a = allocate_source_fixture(&mut heap, vec![b'a' as i8, 0]);
        let ab = allocate_source_fixture(&mut heap, vec![b'a' as i8, b'B' as i8, 0]);
        assert_eq!(
            inchi_stricmp(&heap, empty.as_const(), empty.as_const()),
            Ok(0)
        );
        assert_eq!(
            inchi_stricmp(&heap, a.as_const(), lower_a.as_const()),
            Ok(0)
        );
        assert_eq!(inchi_stricmp(&heap, a.as_const(), ab.as_const()), Ok(-1));
        assert_eq!(inchi_stricmp(&heap, ab.as_const(), a.as_const()), Ok(98));

        let negative = allocate_source_fixture(&mut heap, vec![0xff_u8 as i8, 0]);
        assert_eq!(
            inchi_stricmp(&heap, negative.as_const(), a.as_const()),
            Ok(-1 - 97)
        );
        assert_eq!(
            inchi_stricmp(&heap, SourceConstPointer::null(), a.as_const()),
            Err(SourceHeapError::NullPointer)
        );
        let unterminated = allocate_source_fixture(&mut heap, vec![b'A' as i8]);
        assert_eq!(
            inchi_stricmp(&heap, unterminated.as_const(), lower_a.as_const()),
            Err(SourceHeapError::MissingNulTerminator)
        );
    }

    #[test]
    fn source_port__util__inchi_strdup__line_2035() {
        let mut heap = SourceHeap::default();
        heap.trace_source_allocations();
        assert_eq!(
            inchi__strdup(&mut heap, SourceConstPointer::null()),
            Ok(SourceMutPointer::null())
        );
        assert_eq!(heap.source_allocation_calls(), 0);

        let empty = allocate_source_fixture(&mut heap, vec![0_i8, 77]);
        heap.trace_source_allocations();
        let empty_copy = inchi__strdup(&mut heap, empty.as_const()).unwrap();
        assert_eq!(heap.source_allocation_calls(), 1);
        assert_eq!(heap.slice(empty_copy.as_const()).unwrap(), &[0]);

        let source = allocate_source_fixture(
            &mut heap,
            vec![b'A' as i8, 0xff_u8 as i8, b'Z' as i8, 0, 88, 99],
        );
        heap.trace_source_allocations();
        let duplicate = inchi__strdup(&mut heap, source.as_const()).unwrap();
        assert_eq!(heap.source_allocation_calls(), 1);
        assert_ne!(duplicate, source);
        assert_eq!(
            heap.slice(duplicate.as_const()).unwrap(),
            &[b'A' as i8, 0xff_u8 as i8, b'Z' as i8, 0]
        );
        heap.slice_mut(duplicate).unwrap()[0] = b'Q' as i8;
        assert_eq!(heap.slice(source.as_const()).unwrap()[0], b'A' as i8);

        heap.fail_after_allocations(0);
        assert_eq!(
            inchi__strdup(&mut heap, source.as_const()),
            Err(SourceHeapError::AllocationFailed)
        );

        let unterminated = allocate_source_fixture(&mut heap, vec![b'x' as i8]);
        heap.trace_source_allocations();
        assert_eq!(
            inchi__strdup(&mut heap, unterminated.as_const()),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(heap.source_allocation_calls(), 0);
    }

    #[test]
    fn source_port__util__is_in_the_list__line_1059() {
        assert_eq!(is_in_the_list(None, 7, 0), Ok(None));

        let path = [2_u16, 7, 7, 9];
        assert_eq!(is_in_the_list(Some(&path), 7, 4), Ok(Some(1)));
        assert_eq!(is_in_the_list(Some(&path), 9, 2), Ok(None));
        assert_eq!(is_in_the_list(Some(&path), 9, 4), Ok(Some(3)));
        assert_eq!(
            is_in_the_list(Some(&path), 1, -1),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(
            is_in_the_list(Some(&path), 1, 5),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            is_in_the_list(None, 1, 1),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__util__is_in_the_ilist__line_1072() {
        assert_eq!(is_in_the_ilist(None, 7, 0), Ok(None));

        let path = [i32::MIN, 7, 7, i32::MAX];
        assert_eq!(is_in_the_ilist(Some(&path), 7, 4), Ok(Some(1)));
        assert_eq!(is_in_the_ilist(Some(&path), i32::MAX, 3), Ok(None));
        assert_eq!(is_in_the_ilist(Some(&path), i32::MAX, 4), Ok(Some(3)));
        assert_eq!(
            is_in_the_ilist(Some(&path), 1, -1),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(
            is_in_the_ilist(Some(&path), 1, 5),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            is_in_the_ilist(None, 1, 1),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__util__is_ilist_inside__line_1085() {
        assert_eq!(is_ilist_inside(None, 0, None, 0), Ok(1));
        assert_eq!(is_ilist_inside(None, -1, None, i32::MIN), Ok(1));

        let list = [7, i32::MIN, 7, i32::MAX];
        let embedding = [i32::MAX, 7, 9, i32::MIN];
        assert_eq!(is_ilist_inside(Some(&list), 4, Some(&embedding), 4), Ok(1));
        assert_eq!(is_ilist_inside(Some(&list), 3, Some(&embedding), 3), Ok(0));
        assert_eq!(is_ilist_inside(Some(&[7, 7, 7]), 3, Some(&[7]), 1), Ok(1));
        assert_eq!(is_ilist_inside(Some(&[7]), 1, Some(&embedding), 2), Ok(1));
        assert_eq!(is_ilist_inside(Some(&[7]), 1, None, 0), Ok(0));
        assert_eq!(
            is_ilist_inside(None, 1, Some(&embedding), 4),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(
            is_ilist_inside(Some(&list), 5, Some(&embedding), 4),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(is_ilist_inside(Some(&[7]), 1, Some(&embedding), 5), Ok(1));
        assert_eq!(
            is_ilist_inside(Some(&[8]), 1, Some(&embedding), 5),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            is_ilist_inside(Some(&[7]), 1, Some(&embedding), -1),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(
            is_ilist_inside(Some(&[8, 7]), 2, Some(&embedding), 1),
            Ok(0)
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__util__is_in_the_list__line_1059() {
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
            .arg("--list-records")
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
        let mut record_count = 0;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["operation"], "is_in_the_list");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let path = official["input"]["path"]
                .as_array()
                .expect("path must be an array")
                .iter()
                .map(|value| value.as_u64().unwrap() as AT_NUMB)
                .collect::<Vec<_>>();
            let target = official["input"]["target"].as_u64().unwrap() as AT_NUMB;
            let path_length = official["input"]["path_length"].as_i64().unwrap() as i32;
            let rust = if case_id == "null-zero" {
                is_in_the_list(None, target, path_length)
            } else {
                is_in_the_list(Some(&path), target, path_length)
            }
            .unwrap_or_else(|error| panic!("{case_id}: {error:?}"));
            let expected_offset = official["output"]["offset"].as_i64().unwrap();
            let expected = usize::try_from(expected_offset).ok();
            assert_eq!(rust, expected, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 8_021);
    }

    #[test]
    fn source_port__util__get_element_chemical_symbol__line_289() {
        for (atomic_number, expected_result, expected_symbol) in [
            (-1, -1, b"??".as_slice()),
            (0, -1, b"??".as_slice()),
            (1, 0, b"H".as_slice()),
            (2, 0, b"He".as_slice()),
            (5, 0, b"B".as_slice()),
            (6, 0, b"C".as_slice()),
            (7, 0, b"N".as_slice()),
            (8, 0, b"O".as_slice()),
            (17, 0, b"Cl".as_slice()),
            (53, 0, b"I".as_slice()),
            (118, 0, b"Og".as_slice()),
            (119, 0, b"Zy".as_slice()),
            (120, 0, b"Zz".as_slice()),
            (121, -1, b"??".as_slice()),
            (i32::MAX - 2, -1, b"??".as_slice()),
        ] {
            let mut output = [b'!' as i8; 8];
            assert_eq!(
                get_element_chemical_symbol(atomic_number, &mut output),
                Ok(expected_result),
                "atomic number {atomic_number}"
            );
            let expected = expected_symbol
                .iter()
                .copied()
                .chain(std::iter::once(0))
                .map(|byte| byte as i8)
                .collect::<Vec<_>>();
            assert_eq!(&output[..expected.len()], expected.as_slice());
            assert!(
                output[expected.len()..]
                    .iter()
                    .all(|byte| *byte == b'!' as i8)
            );
        }

        let mut short = [0_i8; 2];
        assert_eq!(
            get_element_chemical_symbol(2, &mut short),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__util__el_number_in_internal_ref_table__line_347() {
        assert_eq!(
            el_number_in_internal_ref_table(None),
            Err(SourceHeapError::NullPointer)
        );

        // Fixed expected indices are emitted by the independent official C
        // oracle; they are not derived from EL_DATA_SYMBOLS.
        for (symbol, expected) in [
            (b"H\0".as_slice(), 0),
            (b"D\0".as_slice(), 1),
            (b"T\0".as_slice(), 2),
            (b"He\0".as_slice(), 3),
            (b"C\0".as_slice(), 7),
            (b"Og\0".as_slice(), 119),
            (b"Zy\0".as_slice(), 120),
            (b"Zz\0".as_slice(), 121),
        ] {
            let c_string: Vec<i8> = symbol.iter().map(|byte| *byte as i8).collect();
            assert_eq!(
                el_number_in_internal_ref_table(Some(&c_string)),
                Ok(expected),
                "symbol {}",
                String::from_utf8_lossy(&symbol[..symbol.len() - 1])
            );
        }

        for missing in ["", "h", "Carbon", "C1", "Zu", "Zv", "Zw", "Zx", "??"] {
            let mut c_string: Vec<i8> = missing.bytes().map(|byte| byte as i8).collect();
            c_string.push(0);
            assert_eq!(
                el_number_in_internal_ref_table(Some(&c_string)),
                Ok(ERR_ELEM),
                "missing symbol {missing}"
            );
        }

        assert_eq!(
            el_number_in_internal_ref_table(Some(&[b'C' as i8, 0, b'l' as i8, 0])),
            Ok(7)
        );
        assert_eq!(
            el_number_in_internal_ref_table(Some(&[b'C' as i8])),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(
            el_number_in_internal_ref_table(Some(&[-1_i8, 0])),
            Ok(ERR_ELEM)
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__util__el_number_in_internal_ref_table__line_347() {
        use std::collections::BTreeMap;
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
            .arg("--element-lookup-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );

        let records = String::from_utf8(oracle.stdout)
            .expect("official C oracle output must be UTF-8")
            .lines()
            .map(|line| serde_json::from_str::<Value>(line).expect("oracle record must be JSON"))
            .collect::<Vec<_>>();
        assert_eq!(records.len(), 718);

        let mut discovered = BTreeMap::new();
        for official in records {
            assert_eq!(official["operation"], "el_number_in_internal_ref_table");
            assert_eq!(official["output"]["n_el_data_len"], 122);
            assert_eq!(official["output"]["err_elem"], 255);
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let symbol = official["input"]["symbol_bytes"]
                .as_array()
                .expect("symbol bytes must be an array")
                .iter()
                .map(|value| value.as_u64().expect("symbol byte must be u8") as u8)
                .collect::<Vec<_>>();
            let rust_symbol = symbol.iter().map(|byte| *byte as i8).collect::<Vec<_>>();
            let rust = el_number_in_internal_ref_table(Some(&rust_symbol))
                .unwrap_or_else(|error| panic!("{case_id}: {error:?}"));
            let expected = official["output"]["result"]
                .as_i64()
                .and_then(|value| i32::try_from(value).ok())
                .expect("result must be i32");
            assert_eq!(rust, expected, "{case_id}");
            if case_id.starts_with("candidate-") && expected != ERR_ELEM {
                assert!(
                    discovered.insert(expected, symbol).is_none(),
                    "duplicate official table index {expected}"
                );
            }
        }
        assert_eq!(discovered.len(), 122);
        assert_eq!(
            discovered.keys().copied().collect::<Vec<_>>(),
            (0..122).collect::<Vec<_>>()
        );
    }

    #[test]
    fn source_port__util__get_periodic_table_number__line_364() {
        assert_eq!(get_periodic_table_number(None), Ok(ERR_ELEM));
        assert_eq!(get_periodic_table_number(Some(&[0])), Ok(ERR_ELEM));

        for (symbol, expected) in [
            (b"H\0".as_slice(), 1),
            (b"D\0".as_slice(), 1),
            (b"T\0".as_slice(), 1),
            (b"He\0".as_slice(), 2),
            (b"B\0".as_slice(), 5),
            (b"C\0".as_slice(), 6),
            (b"N\0".as_slice(), 7),
            (b"O\0".as_slice(), 8),
            (b"P\0".as_slice(), 15),
            (b"S\0".as_slice(), 16),
            (b"F\0".as_slice(), 9),
            (b"I\0".as_slice(), 53),
            (b"Og\0".as_slice(), 118),
            (b"Zy\0".as_slice(), 119),
            (b"Zz\0".as_slice(), 120),
        ] {
            let c_string: Vec<i8> = symbol.iter().map(|byte| *byte as i8).collect();
            assert_eq!(get_periodic_table_number(Some(&c_string)), Ok(expected));
        }

        for missing in ["X", "h", "Carbon", "Zu", "??"] {
            let mut c_string: Vec<i8> = missing.bytes().map(|byte| byte as i8).collect();
            c_string.push(0);
            assert_eq!(
                get_periodic_table_number(Some(&c_string)),
                Ok(ERR_ELEM),
                "missing symbol {missing}"
            );
        }
        assert_eq!(
            get_periodic_table_number(Some(&[b'C' as i8, 0, b'l' as i8, 0])),
            Ok(6)
        );
        assert_eq!(
            get_periodic_table_number(Some(&[b'C' as i8])),
            Err(SourceHeapError::MissingNulTerminator)
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__util__get_periodic_table_number__line_364() {
        use std::collections::BTreeSet;
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
            .arg("--periodic-lookup-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let records = String::from_utf8(oracle.stdout)
            .expect("official C oracle output must be UTF-8")
            .lines()
            .map(|line| serde_json::from_str::<Value>(line).expect("oracle record must be JSON"))
            .collect::<Vec<_>>();
        assert_eq!(records.len(), 719);

        let mut valid_candidate_count = 0;
        let mut periodic_numbers = BTreeSet::new();
        for official in records {
            assert_eq!(official["operation"], "get_periodic_table_number");
            assert_eq!(official["output"]["err_elem"], 255);
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let is_null = official["input"]["is_null"]
                .as_bool()
                .expect("is_null must be bool");
            let rust = if is_null {
                get_periodic_table_number(None).unwrap()
            } else {
                let symbol = official["input"]["symbol_bytes"]
                    .as_array()
                    .expect("symbol bytes must be an array")
                    .iter()
                    .map(|value| value.as_u64().expect("symbol byte must be u8") as i8)
                    .collect::<Vec<_>>();
                get_periodic_table_number(Some(&symbol))
                    .unwrap_or_else(|error| panic!("{case_id}: {error:?}"))
            };
            let expected = official["output"]["result"]
                .as_i64()
                .and_then(|value| i32::try_from(value).ok())
                .expect("result must be i32");
            assert_eq!(rust, expected, "{case_id}");
            if case_id.starts_with("candidate-") && expected != ERR_ELEM {
                valid_candidate_count += 1;
                periodic_numbers.insert(expected);
            }
        }
        assert_eq!(valid_candidate_count, 122);
        assert_eq!(periodic_numbers, (1..=120).collect::<BTreeSet<_>>());
    }

    #[test]
    fn source_port__util__if_skip_add_h__line_428() {
        for periodic_number in -4_i32..=120 {
            let expected = i32::from(
                (21..=30).contains(&periodic_number)
                    || (39..=48).contains(&periodic_number)
                    || (57..=120).contains(&periodic_number),
            );
            assert_eq!(
                if_skip_add_H(periodic_number),
                Ok(expected),
                "periodic number {periodic_number}"
            );
        }
        assert_eq!(if_skip_add_H(121), Err(SourceHeapError::PointerOutOfBounds));
        assert_eq!(
            if_skip_add_H(i32::MAX),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(if_skip_add_H(i32::MIN), Ok(0));
    }

    #[test]
    fn source_port__util__get_el_valence__line_439() {
        // Fixed values come from the independent official C table oracle.
        for (periodic_number, charge, slot, expected) in [
            (2, 0, 0, 0),
            (6, -2, 0, 2),
            (6, 0, 0, 4),
            (7, 0, 0, 3),
            (7, 0, 1, 5),
            (8, 1, 0, 3),
            (8, 1, 1, 5),
            (15, 0, 0, 3),
            (15, 0, 1, 5),
            (118, 0, 0, 1),
            (119, 0, 0, 1),
            (120, 0, 0, 1),
        ] {
            assert_eq!(
                get_el_valence(periodic_number, charge, slot),
                Ok(expected),
                "periodic={periodic_number} charge={charge} slot={slot}"
            );
        }

        for periodic_number in [i32::MIN, -1, 0, 1] {
            assert_eq!(get_el_valence(periodic_number, 0, 0), Ok(1));
            assert_eq!(get_el_valence(periodic_number, 0, 1), Ok(0));
        }
        for charge in MIN_ATOM_CHARGE..=MAX_ATOM_CHARGE as i32 {
            for valence_number in 0..MAX_NUM_VALENCES as i32 {
                assert_eq!(get_el_valence(121, charge, valence_number), Ok(0));
            }
        }

        assert_eq!(get_el_valence(6, -3, 0), Ok(0));
        assert_eq!(get_el_valence(6, 3, 0), Ok(0));
        assert_eq!(get_el_valence(i32::MAX, 3, -1), Ok(0));
        assert_eq!(get_el_valence(6, 0, 5), Ok(0));
        assert_eq!(get_el_valence(6, 0, i32::MAX), Ok(0));
        assert_eq!(
            get_el_valence(6, 0, -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            get_el_valence(122, 0, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            get_el_valence(i32::MAX, 0, 0),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__util__get_el_valence__line_439() {
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
            .arg("--el-valence-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let records = String::from_utf8(oracle.stdout)
            .expect("official C oracle output must be UTF-8")
            .lines()
            .map(|line| serde_json::from_str::<Value>(line).expect("oracle record must be JSON"))
            .collect::<Vec<_>>();
        assert_eq!(records.len(), 5208);

        for official in records {
            assert_eq!(official["operation"], "get_el_valence");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let input = &official["input"];
            let periodic_number = input["periodic_number"].as_i64().unwrap() as i32;
            let charge = input["charge"].as_i64().unwrap() as i32;
            let slot = input["slot"].as_i64().unwrap() as i32;
            let rust = get_el_valence(periodic_number, charge, slot)
                .unwrap_or_else(|error| panic!("{case_id}: {error:?}"));
            let expected = official["output"]["result"].as_i64().unwrap() as i32;
            assert_eq!(rust, expected, "{case_id}");
        }
    }

    #[test]
    fn source_port__util__get_unusual_el_valence__line_454() {
        for (case, input, expected) in [
            ("empty", (i32::MAX, 99, 99, i32::MAX, 0, 0), 0),
            ("charge-high-single", (6, 3, 0, 2, 1, 2), 0),
            ("charge-low-nonsingle", (6, -3, 0, 2, 1, 1), 2),
            ("unknown-element-single", (121, 0, 0, 2, 1, 2), 0),
            ("carbon-exact", (6, 0, 0, 3, 1, 3), 0),
            ("carbon-nonexact", (6, 0, 0, 2, 1, 2), 3),
            ("nitrogen-first-exact", (7, 0, 0, 1, 2, 1), 0),
            ("nitrogen-multiple-fit", (7, 0, 0, 1, 4, 1), 5),
            (
                "doublet-adjustment",
                (6, 0, RADICAL_DOUBLET as i32, 2, 1, 2),
                0,
            ),
            (
                "triplet-adjustment",
                (6, 0, RADICAL_TRIPLET as i32, 1, 1, 1),
                0,
            ),
            (
                "singlet-is-not-adjusted",
                (6, 0, RADICAL_SINGLET as i32, 1, 1, 1),
                2,
            ),
            ("unknown-radical-is-not-adjusted", (6, 0, 99, 1, 1, 1), 2),
            ("negative-chemical-valence", (6, 0, 0, -3, 1, 1), -2),
        ] {
            let (periodic, charge, radical, bonds, hydrogens, bond_count) = input;
            assert_eq!(
                get_unusual_el_valence(periodic, charge, radical, bonds, hydrogens, bond_count,),
                Ok(expected),
                "{case}"
            );
        }

        assert_eq!(
            get_unusual_el_valence(6, 0, 0, i32::MAX, 1, 1),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(
            get_unusual_el_valence(6, 0, 0, i32::MIN, -1, 1),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(
            get_unusual_el_valence(i32::MAX, 0, 0, 1, 1, 0),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
    }

    #[test]
    fn source_port__util__needed_unusual_el_valence__line_518() {
        for (case, input, expected) in [
            ("carbon-exact", (6, 0, 0, 4, 4, 0, 1), 0),
            ("carbon-implicit-h", (6, 0, 0, 3, 3, 1, 1), 0),
            ("no-bonds-zero-needs-zero", (6, 0, 0, 0, 0, 0, 0), -1),
            ("helium-no-h", (2, 0, 0, 0, 0, 0, 1), 0),
            ("invalid-charge", (6, 3, 0, 2, 2, 1, 1), 3),
            ("skip-h-no-h", (21, 0, 0, 3, 3, 0, 1), 0),
            ("skip-h-with-h", (21, 0, 0, 3, 3, 1, 1), 4),
            ("bond-valence-mismatch", (6, 0, 0, 3, 4, 1, 1), 4),
            ("expected-h-mismatch", (6, 0, 0, 4, 4, 1, 1), 5),
            ("nitrogen-exact-first", (7, 0, 0, 2, 2, 1, 1), 0),
            ("sulfur-multiple-known", (16, 0, 0, 3, 3, 0, 1), 3),
            (
                "doublet-adjustment",
                (6, 0, RADICAL_DOUBLET as i32, 3, 3, 0, 1),
                0,
            ),
            (
                "triplet-adjustment",
                (6, 0, RADICAL_TRIPLET as i32, 2, 2, 0, 1),
                0,
            ),
            ("unknown-radical-zero", (6, 0, 99, 0, 0, 0, 1), -1),
            ("invalid-element-no-h", (121, 0, 0, 0, 0, 0, 1), 0),
        ] {
            let (periodic, charge, radical, bonds, actual, hydrogens, bond_count) = input;
            assert_eq!(
                needed_unusual_el_valence(
                    periodic, charge, radical, bonds, actual, hydrogens, bond_count,
                ),
                Ok(expected),
                "{case}"
            );
        }

        assert_eq!(
            needed_unusual_el_valence(i32::MAX, 0, 0, 0, 0, 0, 1),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
    }

    #[test]
    fn source_port__util__detect_unusual_el_valence__line_620() {
        assert_eq!(
            detect_unusual_el_valence(i32::MAX, 99, 99, i32::MAX, 0, 0),
            Ok(0)
        );
        assert_eq!(detect_unusual_el_valence(6, 3, 0, 2, 2, 2), Ok(0));
        assert_eq!(detect_unusual_el_valence(6, 3, 0, 2, 1, 1), Ok(2));
        assert_eq!(detect_unusual_el_valence(2, 0, 0, 2, 0, 2), Ok(0));
        assert_eq!(detect_unusual_el_valence(6, 0, 0, 3, 1, 3), Ok(0));
        assert_eq!(detect_unusual_el_valence(6, 0, 0, 3, 0, 3), Ok(3));
        assert_eq!(detect_unusual_el_valence(6, 0, 2, 3, 0, 3), Ok(0));
        assert_eq!(detect_unusual_el_valence(6, 0, 1, 2, 0, 1), Ok(0));
        assert_eq!(detect_unusual_el_valence(6, 0, 3, 2, 0, 1), Ok(0));
        assert_eq!(detect_unusual_el_valence(6, 0, -1, 0, 4, 1), Ok(0));
        assert_eq!(detect_unusual_el_valence(6, 0, 4, 0, 4, 1), Ok(0));
        assert_eq!(detect_unusual_el_valence(7, 0, 0, 0, 5, 1), Ok(0));
        assert_eq!(detect_unusual_el_valence(2, 0, 3, 0, 1, 1), Ok(1));
        assert_eq!(detect_unusual_el_valence(121, 0, 0, 0, 7, 1), Ok(7));
        assert_eq!(detect_unusual_el_valence(i32::MIN, 0, 0, 0, 1, 1), Ok(0));
        for (bonds_valence, hydrogen_count, expected) in [
            (i32::MIN, 0, i32::MIN),
            (i32::MAX, 0, i32::MAX),
            (i32::MIN, i32::MAX, -1),
            (i32::MAX, i32::MIN, -1),
        ] {
            assert_eq!(
                detect_unusual_el_valence(6, 0, 0, bonds_valence, hydrogen_count, 1,),
                Ok(expected)
            );
        }
        assert_eq!(
            detect_unusual_el_valence(6, 0, 0, i32::MAX, 1, 1),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(
            detect_unusual_el_valence(6, 0, 0, i32::MIN, -1, 1),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__util__detect_unusual_el_valence__line_620() {
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
            .arg("--detect-unusual-valence-records")
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
        let mut record_count = 0;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["operation"], "detect_unusual_el_valence");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let input = &official["input"];
            let read_i32 = |field: &str| {
                input[field]
                    .as_i64()
                    .and_then(|value| i32::try_from(value).ok())
                    .unwrap_or_else(|| panic!("{case_id}: {field} must be i32"))
            };
            let rust = detect_unusual_el_valence(
                read_i32("periodic_number"),
                read_i32("charge"),
                read_i32("radical"),
                read_i32("bonds_valence"),
                read_i32("hydrogen_count"),
                read_i32("bond_count"),
            )
            .unwrap_or_else(|error| panic!("{case_id}: {error:?}"));
            let expected = official["output"]["result"]
                .as_i64()
                .and_then(|value| i32::try_from(value).ok())
                .expect("result must be i32");
            assert_eq!(rust, expected, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 131_956);
    }

    #[test]
    fn source_port__util__extract_charges_and_radicals__line_700() {
        fn assert_case(
            input: &str,
            expected_name: &str,
            expected_radical: i32,
            expected_charge: i32,
            expected_status: i32,
        ) {
            let mut bytes: Vec<i8> = input.bytes().map(|byte| byte as i8).collect();
            bytes.extend([0, 0x55, 0x66]);
            let mut radical = -91;
            let mut charge = -92;

            assert_eq!(
                extract_charges_and_radicals(
                    Some(&mut bytes),
                    Some(&mut radical),
                    Some(&mut charge),
                ),
                Ok(expected_status),
                "input {input:?}"
            );
            let nul = bytes.iter().position(|byte| *byte == 0).unwrap();
            let actual_name: Vec<u8> = bytes[..nul].iter().map(|byte| *byte as u8).collect();
            assert_eq!(actual_name, expected_name.as_bytes(), "input {input:?}");
            assert_eq!(radical, expected_radical, "input {input:?}");
            assert_eq!(charge, expected_charge, "input {input:?}");
        }

        for (input, name, radical, charge, status) in [
            ("", "", 0, 0, 0),
            ("C", "C", 0, 0, 0),
            ("C+", "C", 0, 1, 1),
            ("C-", "C", 0, -1, 1),
            ("C++", "C", 0, 2, 1),
            ("C---", "C", 0, -3, 1),
            ("C+0", "C", 0, 1, 1),
            ("C-0", "C", 0, -1, 1),
            ("C+2", "C", 0, 2, 1),
            ("C-2", "C", 0, -2, 1),
            ("C+-", "C", 0, 0, 0),
            ("C-+", "C", 0, 0, 0),
            ("C+-2", "C", 0, -1, 1),
            ("C-+2", "C", 0, 1, 1),
            ("C+2H-3", "CH", 0, -1, 1),
            ("C+2H+3", "CH", 0, 5, 1),
            ("C^", "C", RADICAL_DOUBLET as i32, 0, 1),
            ("C^^", "C", RADICAL_TRIPLET as i32, 0, 1),
            ("C^^^", "C", 0, 0, 0),
            ("C^x^^", "Cx", RADICAL_TRIPLET as i32, 0, 1),
            ("C:", "C", RADICAL_SINGLET as i32, 0, 1),
            ("C.", "C", RADICAL_DOUBLET as i32, 0, 1),
            ("C..", "C", RADICAL_TRIPLET as i32, 0, 1),
            ("C...", "C", 0, 0, 0),
            ("C.O", "C.O", 0, 0, 0),
            ("C:O", "C:O", 0, 0, 0),
            ("C^.", "C", RADICAL_TRIPLET as i32, 0, 1),
            ("C^^.", "C", 0, 0, 0),
            ("C^^:", "C", RADICAL_SINGLET as i32, 0, 1),
            ("Fe+ 2", "Fe", 0, 2, 1),
            ("Fe+\t-2", "Fe", 0, -2, 1),
            ("Fe+abc", "Feabc", 0, 1, 1),
            ("Fe+  abc", "Fe  abc", 0, 1, 1),
        ] {
            assert_case(input, name, radical, charge, status);
        }

        let parse = |input: &str| {
            let bytes: Vec<i8> = input.bytes().map(|byte| byte as i8).collect();
            source_strtol_base_10(&bytes, 0)
        };
        assert_eq!(parse("2147483647x"), Ok((i32::MAX, 10)));
        assert_eq!(parse("2147483648x"), Ok((i32::MIN, 10)));
        assert_eq!(parse("4294967295x"), Ok((-1, 10)));
        assert_eq!(parse("4294967296x"), Ok((0, 10)));
        assert_eq!(parse("9223372036854775807x"), Ok((-1, 19)));
        assert_eq!(parse("9223372036854775808x"), Ok((-1, 19)));
        assert_eq!(parse("-9223372036854775808x"), Ok((0, 20)));
        assert_eq!(parse("-9223372036854775809x"), Ok((0, 20)));
        assert_eq!(parse("-2147483649x"), Ok((i32::MAX, 11)));
        assert_eq!(parse(" \t+12x"), Ok((12, 5)));
        assert_eq!(parse(" \t+x"), Ok((0, 0)));

        let mut exact = b"C+2H-3\0\x55\x66".map(|byte| byte as i8);
        let mut radical = -1;
        let mut charge = -1;
        assert_eq!(
            extract_charges_and_radicals(Some(&mut exact), Some(&mut radical), Some(&mut charge),),
            Ok(1)
        );
        assert_eq!(
            exact,
            [
                b'C' as i8, b'H' as i8, 0, b'3' as i8, 0, b'3' as i8, 0, 0x55, 0x66,
            ]
        );
        assert_eq!((radical, charge), (0, -1));

        let mut unchanged = b"C+2\0\x55".map(|byte| byte as i8);
        let original = unchanged;
        let mut radical = 17;
        let mut charge = 19;
        assert_eq!(
            extract_charges_and_radicals(None, Some(&mut radical), Some(&mut charge)),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!((radical, charge), (17, 19));
        assert_eq!(
            extract_charges_and_radicals(Some(&mut unchanged), None, Some(&mut charge)),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(unchanged, original);
        assert_eq!((radical, charge), (17, 19));
        assert_eq!(
            extract_charges_and_radicals(Some(&mut unchanged), Some(&mut radical), None),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(unchanged, original);
        assert_eq!((radical, charge), (17, 19));

        let mut unterminated = b"C+2".map(|byte| byte as i8);
        assert_eq!(
            extract_charges_and_radicals(
                Some(&mut unterminated),
                Some(&mut radical),
                Some(&mut charge),
            ),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(unterminated, b"C+2".map(|byte| byte as i8));
        assert_eq!((radical, charge), (17, 19));

        let mut local_overflow = b"X++2147483647\0".map(|byte| byte as i8);
        assert_eq!(
            extract_charges_and_radicals(
                Some(&mut local_overflow),
                Some(&mut radical),
                Some(&mut charge),
            ),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(local_overflow, b"X++2147483647\0".map(|byte| byte as i8));
        assert_eq!((radical, charge), (17, 19));

        let mut narrowed_overflow = b"X+-2147483648\0".map(|byte| byte as i8);
        assert_eq!(
            extract_charges_and_radicals(
                Some(&mut narrowed_overflow),
                Some(&mut radical),
                Some(&mut charge),
            ),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!((radical, charge), (17, 19));

        let mut accumulated_overflow = b"X+2147483647Y+\0".map(|byte| byte as i8);
        assert_eq!(
            extract_charges_and_radicals(
                Some(&mut accumulated_overflow),
                Some(&mut radical),
                Some(&mut charge),
            ),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(
            &accumulated_overflow[..4],
            &[b'X' as i8, b'Y' as i8, b'+' as i8, 0]
        );
        assert_eq!((radical, charge), (17, 19));
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__util__extract_charges_and_radicals__line_700() {
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
            .arg("--extract-charge-records")
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
        let mut record_count = 0;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["operation"], "extract_charges_and_radicals");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let text = official["input"]["text"]
                .as_str()
                .expect("input text must be text");
            let mut buffer = vec![0x55_i8; 128];
            for (target, source) in buffer.iter_mut().zip(text.bytes()) {
                *target = source as i8;
            }
            buffer[text.len()] = 0;
            let mut radical = -1901;
            let mut charge = -1902;
            let status = extract_charges_and_radicals(
                Some(&mut buffer),
                Some(&mut radical),
                Some(&mut charge),
            )
            .unwrap_or_else(|error| panic!("{case_id}: {error:?}"));
            let expected = &official["output"];
            assert_eq!(
                status,
                expected["status"].as_i64().unwrap() as i32,
                "{case_id}"
            );
            assert_eq!(
                radical,
                expected["radical"].as_i64().unwrap() as i32,
                "{case_id}"
            );
            assert_eq!(
                charge,
                expected["charge"].as_i64().unwrap() as i32,
                "{case_id}"
            );
            let nul_offset = buffer
                .iter()
                .position(|byte| *byte == 0)
                .expect("result must remain NUL terminated");
            assert_eq!(
                nul_offset,
                expected["nul_offset"].as_u64().unwrap() as usize,
                "{case_id}"
            );
            let expected_buffer = expected["buffer"]
                .as_array()
                .expect("buffer must be an array")
                .iter()
                .map(|value| value.as_i64().unwrap() as i8)
                .collect::<Vec<_>>();
            assert_eq!(buffer, expected_buffer, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 1_060);
    }

    #[test]
    fn source_port__util__extract_h_atoms__line_774() {
        fn assert_case(
            input: &[u8],
            initial_isotopes: [S_CHAR; 3],
            expected_name: &[u8],
            expected_hydrogens: i32,
            expected_isotopes: [S_CHAR; 3],
        ) {
            let mut bytes: Vec<i8> = input.iter().map(|byte| *byte as i8).collect();
            bytes.push(0);
            bytes.resize(bytes.len() + 32, 0x55);
            let mut isotopes = initial_isotopes;

            assert_eq!(
                extract_h_atoms(Some(&mut bytes), Some(&mut isotopes)),
                Ok(expected_hydrogens),
                "input {:?}",
                String::from_utf8_lossy(input)
            );
            let nul = bytes.iter().position(|byte| *byte == 0).unwrap();
            let actual_name: Vec<u8> = bytes[..nul].iter().map(|byte| *byte as u8).collect();
            assert_eq!(
                actual_name,
                expected_name,
                "input {:?}",
                String::from_utf8_lossy(input)
            );
            assert_eq!(isotopes, expected_isotopes);
        }

        for (input, name, hydrogens, isotopes) in [
            (b"".as_slice(), b"".as_slice(), 0, [0, 0, 0]),
            (b"C".as_slice(), b"C".as_slice(), 0, [0, 0, 0]),
            (b"H".as_slice(), b"".as_slice(), 1, [0, 0, 0]),
            (b"H2".as_slice(), b"".as_slice(), 2, [0, 0, 0]),
            (b"H0".as_slice(), b"".as_slice(), 0, [0, 0, 0]),
            (b"H01".as_slice(), b"".as_slice(), 1, [0, 0, 0]),
            (b"D".as_slice(), b"".as_slice(), 0, [0, 1, 0]),
            (b"D2".as_slice(), b"".as_slice(), 0, [0, 2, 0]),
            (b"T".as_slice(), b"".as_slice(), 0, [0, 0, 1]),
            (b"T3".as_slice(), b"".as_slice(), 0, [0, 0, 3]),
            (b"H2D3T4".as_slice(), b"".as_slice(), 2, [0, 3, 4]),
            (b"HH".as_slice(), b"".as_slice(), 2, [0, 0, 0]),
            (b"DH".as_slice(), b"".as_slice(), 1, [0, 1, 0]),
            (b"HC".as_slice(), b"C".as_slice(), 1, [0, 0, 0]),
            (b"CH3".as_slice(), b"C".as_slice(), 3, [0, 0, 0]),
            (b"HaH2".as_slice(), b"Ha".as_slice(), 2, [0, 0, 0]),
            (b"H2a".as_slice(), b"a".as_slice(), 2, [0, 0, 0]),
            (b"He".as_slice(), b"He".as_slice(), 0, [0, 0, 0]),
            (b"Dy".as_slice(), b"Dy".as_slice(), 0, [0, 0, 0]),
            (b"Th".as_slice(), b"Th".as_slice(), 0, [0, 0, 0]),
            (b"pH4d".as_slice(), b"p?".as_slice(), 4, [0, 0, 0]),
            (b"CH3x".as_slice(), b"C?".as_slice(), 3, [0, 0, 0]),
            (b"CD2x".as_slice(), b"C?".as_slice(), 0, [0, 2, 0]),
            (b"CXD2".as_slice(), b"CX".as_slice(), 0, [0, 2, 0]),
        ] {
            assert_case(input, [0, 0, 0], name, hydrogens, isotopes);
        }

        assert_case(b"H2D3T4", [9, 10, 11], b"", 2, [9, 13, 15]);
        assert_case(b"D", [0, S_CHAR::MAX, 0], b"", 0, [0, S_CHAR::MIN, 0]);
        assert_case(b"D255", [0, 0, 0], b"", 0, [0, -1, 0]);
        assert_case(b"D2147483647", [0, 0, 0], b"", 0, [0, -1, 0]);
        assert_case(b"D2147483648", [0, 0, 0], b"", 0, [0, 0, 0]);
        assert_case(b"H2147483648", [0, 0, 0], b"", i32::MIN, [0, 0, 0]);
        assert_case(b"H9223372036854775808", [0, 0, 0], b"", -1, [0, 0, 0]);

        let mut non_ascii = [b'H' as i8, 0xe0_u8 as i8, 0, 0x55, 0x66];
        let non_ascii_original = non_ascii;
        let mut non_ascii_isotopes = [3, 5, 7];
        assert_eq!(
            extract_h_atoms(Some(&mut non_ascii), Some(&mut non_ascii_isotopes)),
            Err(SourceHeapError::InvalidSourceTextEncoding)
        );
        assert_eq!(non_ascii, non_ascii_original);
        assert_eq!(non_ascii_isotopes, [3, 5, 7]);

        let mut exact = b"CH3X\0AB".map(|byte| byte as i8);
        let mut isotopes = [0; 3];
        assert_eq!(
            extract_h_atoms(Some(&mut exact), Some(&mut isotopes)),
            Ok(3)
        );
        assert_eq!(
            exact,
            [
                b'C' as i8, b'?' as i8, 0, b'A' as i8, 0, b'A' as i8, b'B' as i8
            ]
        );
        assert_eq!(isotopes, [0, 0, 0]);

        let mut unchanged = b"H2\0XYZ".map(|byte| byte as i8);
        let original = unchanged;
        let mut isotopes = [3, 5, 7];
        assert_eq!(
            extract_h_atoms(None, Some(&mut isotopes)),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(isotopes, [3, 5, 7]);
        assert_eq!(
            extract_h_atoms(Some(&mut unchanged), None),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(unchanged, original);
        assert_eq!(isotopes, [3, 5, 7]);

        let mut too_short = [0_i8; 2];
        assert_eq!(
            extract_h_atoms(Some(&mut unchanged), Some(&mut too_short)),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(unchanged, original);

        let mut unterminated = b"H2".map(|byte| byte as i8);
        assert_eq!(
            extract_h_atoms(Some(&mut unterminated), Some(&mut isotopes)),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(unterminated, b"H2".map(|byte| byte as i8));
        assert_eq!(isotopes, [3, 5, 7]);

        let mut source_bounds = b"ABCDH\0".map(|byte| byte as i8);
        assert_eq!(
            extract_h_atoms(Some(&mut source_bounds), Some(&mut isotopes)),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(source_bounds, b"ABCDH\0".map(|byte| byte as i8));
        assert_eq!(isotopes, [3, 6, 7]);

        let mut isotope_overflow = b"D2147483647\0".map(|byte| byte as i8);
        let mut maximum_isotope = [0, S_CHAR::MAX, 0];
        assert_eq!(
            extract_h_atoms(Some(&mut isotope_overflow), Some(&mut maximum_isotope)),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(isotope_overflow, b"D2147483647\0".map(|byte| byte as i8));
        assert_eq!(maximum_isotope, [0, S_CHAR::MAX, 0]);

        let mut count_overflow: Vec<i8> =
            b"H2147483647H\0".iter().map(|byte| *byte as i8).collect();
        count_overflow.resize(count_overflow.len() + 16, 0x55);
        let mut zero_isotopes = [0; 3];
        assert_eq!(
            extract_h_atoms(Some(&mut count_overflow), Some(&mut zero_isotopes)),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(&count_overflow[..2], &[b'H' as i8, 0]);
        assert_eq!(zero_isotopes, [0, 0, 0]);
    }

    #[test]
    fn source_port__util__get_atomic_mass_from_elnum__line_1007() {
        // Fixed values come from the independent official C table oracle.
        for (atomic_number, expected) in [
            (1, 1),
            (2, 4),
            (6, 12),
            (26, 56),
            (53, 127),
            (92, 238),
            (118, 294),
            (119, 0),
            (120, 0),
        ] {
            assert_eq!(
                get_atomic_mass_from_elnum(atomic_number),
                Ok(expected),
                "atomic number {atomic_number}"
            );
        }
        assert_eq!(get_atomic_mass_from_elnum(-1), Ok(0));
        assert_eq!(get_atomic_mass_from_elnum(0), Ok(0));
        assert_eq!(get_atomic_mass_from_elnum(121), Ok(0));
        assert_eq!(get_atomic_mass_from_elnum(i32::MIN + 1), Ok(0));
        assert_eq!(get_atomic_mass_from_elnum(i32::MAX - 1), Ok(0));
        assert_eq!(
            get_atomic_mass_from_elnum(i32::MIN),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(
            get_atomic_mass_from_elnum(i32::MAX),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
    }

    #[test]
    fn source_port__util__get_atomic_mass__line_1040() {
        for (symbol, expected) in [
            (&b"H\0"[..], 1),
            (&b"D\0"[..], 2),
            (&b"T\0"[..], 3),
            (&b"C\0"[..], 12),
            (&b"Fe\0"[..], 56),
            (&b"Og\0"[..], 294),
            (&b"Zy\0"[..], 0),
            (&b"Zz\0"[..], 0),
            (&b"\0"[..], 0),
            (&b"c\0"[..], 0),
            (&b"FE\0"[..], 0),
            (&b"Qq\0"[..], 0),
        ] {
            let input: Vec<i8> = symbol.iter().map(|byte| *byte as i8).collect();
            assert_eq!(get_atomic_mass(Some(&input)), Ok(expected), "{symbol:?}");
        }
        assert_eq!(get_atomic_mass(None), Err(SourceHeapError::NullPointer));
        assert_eq!(
            get_atomic_mass(Some(&[b'C' as i8])),
            Err(SourceHeapError::MissingNulTerminator)
        );
    }

    #[test]
    fn source_port__util__get_num_h__line_862() {
        let c = [b'C' as i8, 0];
        let n = [b'N' as i8, 0];
        let o = [b'O' as i8, 0];
        let s = [b'S' as i8, 0];
        let sc = [b'S' as i8, b'c' as i8, 0];
        let unknown = [b'Q' as i8, 0];
        let call = |element: &[i8],
                    input_h,
                    isotopes: Option<&[i8]>,
                    charge,
                    radical,
                    bonds,
                    input_valence,
                    aliased,
                    no_h,
                    metal| {
            get_num_H(
                Some(element),
                input_h,
                isotopes,
                charge,
                radical,
                bonds,
                input_valence,
                aliased,
                no_h,
                metal,
            )
            .unwrap()
        };

        assert_eq!(
            call(&unknown, i32::MIN, None, 99, 99, 99, 99, 1, 1, 1),
            i32::MIN
        );
        assert_eq!(call(&c, 0, None, 0, 0, 1, 4, 0, 0, 0), 3);
        assert_eq!(call(&c, 0, None, 0, 0, 5, 4, 0, 0, 0), 0);
        assert_eq!(call(&c, 7, None, 0, 0, 0, 15, 0, 0, 0), 0);
        assert_eq!(call(&c, 7, None, MIN_ATOM_CHARGE - 1, 0, 0, 0, 0, 0, 0), 7);
        assert_eq!(
            call(&c, 7, None, MAX_ATOM_CHARGE as i32 + 1, 0, 0, 0, 0, 0, 0,),
            7
        );
        assert_eq!(call(&unknown, 7, None, 0, 0, 0, 0, 0, 0, 0), 7);
        assert_eq!(call(&sc, 7, None, 0, 0, 0, 0, 0, 0, 0), 7);
        assert_eq!(call(&c, 7, None, 0, 0, 0, 0, 0, 1, 0), 7);

        assert_eq!(
            call(&c, 0, None, 0, RADICAL_DOUBLET as i32, 1, 0, 0, 0, 0),
            2
        );
        assert_eq!(
            call(&c, 0, None, 0, RADICAL_TRIPLET as i32, 1, 0, 0, 0, 0),
            1
        );
        assert_eq!(call(&c, 0, None, 0, 99, 1, 0, 0, 0, 0), 0);
        assert_eq!(
            call(&c, 0, None, 0, RADICAL_SINGLET as i32, 1, 0, 0, 0, 0),
            3
        );

        assert_eq!(call(&n, 0, None, 0, 0, 4, 0, 0, 0, 0), 0);
        assert_eq!(call(&s, 0, None, 0, 0, 3, 0, 0, 0, 0), 0);
        assert_eq!(call(&o, 0, None, 0, 0, 1, 0, 0, 0, 0), 1);
        assert_eq!(call(&o, 0, None, 0, 0, 1, 0, 0, 0, 1), 0);
        assert_eq!(call(&c, 0, None, 0, 0, 1, 0, 0, 0, 1), 3);

        assert_eq!(call(&c, 0, Some(&[1, 1, 0]), 0, 0, 0, 0, 0, 0, 0), 2);
        assert_eq!(call(&c, 6, Some(&[3, 3, 3]), 0, 0, 0, 0, 0, 0, 0), 6);
        assert_eq!(call(&c, 5, None, 0, 0, 3, 0, 0, 0, 0), 5);
        assert_eq!(call(&c, 0, Some(&[-1, 0, 0]), 0, 0, 3, 0, 0, 0, 0), 2);
    }

    #[test]
    fn source_port__util__num_of_h__line_1129() {
        let mut heap = SourceHeap::default();
        let mut center = inp_ATOM {
            valence: 3,
            num_H: 2,
            num_iso_H: [1, -1, 3],
            ..inp_ATOM::default()
        };
        center.neighbor[0] = 1;
        center.neighbor[1] = 2;
        center.neighbor[2] = 3;

        let mut h0 = inp_ATOM {
            valence: 1,
            el_number: 1,
            num_H: 1,
            ..inp_ATOM::default()
        };
        h0.neighbor[0] = 0;

        let mut h1 = inp_ATOM {
            valence: 1,
            el_number: 1,
            ..inp_ATOM::default()
        };
        h1.neighbor[0] = 0;

        let mut carbon = inp_ATOM {
            valence: 1,
            el_number: 6,
            ..inp_ATOM::default()
        };
        carbon.neighbor[0] = 0;

        let atoms = heap.allocate(vec![center, h0, h1, carbon]).unwrap();
        let atoms = heap.slice(atoms.as_const()).unwrap();
        assert_eq!(num_of_H(atoms, 0), Ok(7));
        assert_eq!(num_of_H(atoms, 1), Ok(1));
        assert_eq!(
            num_of_H(atoms, -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(num_of_H(atoms, 4), Err(SourceHeapError::PointerOutOfBounds));

        let mut negative_heap = SourceHeap::default();
        let negative = inp_ATOM {
            valence: -1,
            num_H: -4,
            num_iso_H: [-1, 0, 1],
            ..inp_ATOM::default()
        };
        let negative_atoms = negative_heap.allocate(vec![negative]).unwrap();
        assert_eq!(
            num_of_H(negative_heap.slice(negative_atoms.as_const()).unwrap(), 0),
            Ok(-4)
        );
    }

    #[test]
    fn source_port__util__ion_el_group__line_1151() {
        assert_eq!(ion_el_group(6), 6);
        assert_eq!(ion_el_group(7), 7);
        assert_eq!(ion_el_group(15), 7);
        assert_eq!(ion_el_group(33), 7);
        assert_eq!(ion_el_group(51), 7);
        assert_eq!(ion_el_group(8), 8);
        assert_eq!(ion_el_group(16), 8);
        assert_eq!(ion_el_group(34), 8);
        assert_eq!(ion_el_group(52), 8);
        assert_eq!(ion_el_group(14), 0);
        assert_eq!(ion_el_group(1), 0);
        assert_eq!(ion_el_group(-1), 0);
        assert_eq!(ion_el_group(0), 0);
    }

    #[test]
    fn source_port__util__has_other_ion_neigh__line_1175() {
        let mut atoms = vec![inp_ATOM::default(); 5];
        atoms[0].valence = 4;
        atoms[0].neighbor[..4].copy_from_slice(&[1, 2, 3, 4]);
        atoms[1].charge = -1;
        atoms[1].el_number = 8;
        atoms[2].charge = 1;
        atoms[2].el_number = 8;
        atoms[3].charge = -1;
        atoms[3].el_number = 14;
        atoms[4].charge = -1;
        atoms[4].el_number = 16;

        assert_eq!(has_other_ion_neigh(&atoms, 0, 1), Ok(1));
        atoms[4].charge = 0;
        assert_eq!(has_other_ion_neigh(&atoms, 0, 1), Ok(0));
        atoms[4].charge = -1;
        atoms[4].el_number = 14;
        assert_eq!(has_other_ion_neigh(&atoms, 0, 1), Ok(0));

        atoms[0].valence = 0;
        assert_eq!(has_other_ion_neigh(&atoms, 0, 1), Ok(0));
        atoms[0].valence = i8::MIN;
        assert_eq!(has_other_ion_neigh(&atoms, 0, 1), Ok(0));

        assert_eq!(
            has_other_ion_neigh(&atoms, i32::MAX, -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            has_other_ion_neigh(&atoms, -1, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut early_return = vec![inp_ATOM::default(); 3];
        early_return[0].valence = 2;
        early_return[0].neighbor[0] = 2;
        early_return[0].neighbor[1] = u16::MAX;
        early_return[1].charge = i8::MAX;
        early_return[2].charge = i8::MAX;
        early_return[2].el_number = 51;
        assert_eq!(has_other_ion_neigh(&early_return, 0, 1), Ok(1));

        early_return[2].charge = i8::MIN;
        assert_eq!(
            has_other_ion_neigh(&early_return, 0, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__util__has_other_ion_in_sphere_2__line_1201() {
        let mut atoms = vec![inp_ATOM::default(); 8];
        atoms[0].valence = 4;
        atoms[0].neighbor[..4].copy_from_slice(&[1, 2, 5, 6]);
        atoms[0].cFlags = 77;
        atoms[1].valence = 2;
        atoms[1].neighbor[..2].copy_from_slice(&[0, 3]);
        atoms[1].el_number = 8;
        atoms[1].charge = -1;
        atoms[2].valence = 2;
        atoms[2].neighbor[..2].copy_from_slice(&[0, 4]);
        atoms[2].el_number = 7;
        atoms[2].charge = -1;
        atoms[3].valence = 1;
        atoms[3].neighbor[0] = 7;
        atoms[3].el_number = 16;
        atoms[3].charge = -1;
        atoms[4].valence = 1;
        atoms[4].el_number = 14;
        atoms[4].charge = -1;
        atoms[5].valence = 4;
        atoms[5].el_number = 8;
        atoms[5].charge = -1;
        atoms[6].valence = 1;
        atoms[6].el_number = 8;
        atoms[6].charge = -1;
        atoms[6].cFlags = 9;
        atoms[7].el_number = 8;
        atoms[7].charge = -1;

        assert_eq!(has_other_ion_in_sphere_2(&mut atoms, 0, 1), Ok(2));
        assert_eq!(atoms[0].cFlags, 0);
        assert_eq!(atoms[1].cFlags, 0);
        assert_eq!(atoms[2].cFlags, 0);
        assert_eq!(atoms[3].cFlags, 0);
        assert_eq!(atoms[4].cFlags, 0);
        assert_eq!(atoms[5].cFlags, 0);
        assert_eq!(atoms[6].cFlags, 9);
        assert_eq!(atoms[7].cFlags, 0);

        atoms[0].valence = i8::MIN;
        atoms[0].cFlags = -7;
        assert_eq!(has_other_ion_in_sphere_2(&mut atoms, 0, i32::MAX), Ok(0));
        assert_eq!(atoms[0].cFlags, 0);

        let mut full_queue = vec![inp_ATOM::default(); 16];
        full_queue[0].valence = 15;
        for neighbor in 1..16 {
            full_queue[0].neighbor[neighbor - 1] = neighbor as u16;
            full_queue[neighbor].el_number = 8;
            full_queue[neighbor].charge = i8::MAX;
        }
        assert_eq!(has_other_ion_in_sphere_2(&mut full_queue, 0, 1), Ok(14));
        assert!(full_queue.iter().all(|atom| atom.cFlags == 0));

        let mut overflowing_queue = vec![inp_ATOM::default(); 17];
        overflowing_queue[0].valence = 16;
        for neighbor in 1..17 {
            overflowing_queue[0].neighbor[neighbor - 1] = neighbor as u16;
            overflowing_queue[neighbor].el_number = 8;
        }
        assert_eq!(
            has_other_ion_in_sphere_2(&mut overflowing_queue, 0, 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert!(overflowing_queue[..16].iter().all(|atom| atom.cFlags == 1));
        assert_eq!(overflowing_queue[16].cFlags, 0);

        let mut empty = [];
        assert_eq!(
            has_other_ion_in_sphere_2(&mut empty, -1, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__util__get_endpoint_valence__line_1509() {
        assert_eq!(get_endpoint_valence(8), 2);
        assert_eq!(get_endpoint_valence(16), 2);
        assert_eq!(get_endpoint_valence(34), 2);
        assert_eq!(get_endpoint_valence(52), 2);
        assert_eq!(get_endpoint_valence(7), 3);
        assert_eq!(get_endpoint_valence(6), 0);
        assert_eq!(get_endpoint_valence(15), 0);
        assert_eq!(get_endpoint_valence(0), 0);
        assert_eq!(get_endpoint_valence(255), 0);
    }

    #[test]
    fn source_port__util__get_endpoint_valence_ket__line_1530() {
        assert_eq!(get_endpoint_valence_KET(6), 4);
        assert_eq!(get_endpoint_valence_KET(8), 2);
        assert_eq!(get_endpoint_valence_KET(7), 0);
        assert_eq!(get_endpoint_valence_KET(16), 0);
        assert_eq!(get_endpoint_valence_KET(0), 0);
        assert_eq!(get_endpoint_valence_KET(255), 0);
    }

    #[test]
    fn source_port__util__nnometalnumbonds__line_1253() {
        let mut atoms = vec![
            inp_ATOM {
                valence: 2,
                chem_bonds_valence: 3,
                el_number: 8,
                charge: 0,
                num_H: 0,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                valence: 1,
                el_number: 11,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                valence: 1,
                el_number: 6,
                ..inp_ATOM::default()
            },
        ];
        atoms[0].neighbor[0] = 1;
        atoms[0].neighbor[1] = 2;
        atoms[0].bond_type[0] = BOND_SINGLE as u8;
        atoms[0].bond_type[1] = BOND_DOUBLE as u8;
        assert_eq!(n_no_metal_num_bonds(Some(&atoms), 0), Ok(1));

        let mut fallback = atoms.clone();
        fallback[0].bond_type[0] = BOND_TYPE_ALTERN as u8;
        assert_eq!(n_no_metal_num_bonds(Some(&fallback), 0), Ok(2));

        let mut charge_one = vec![
            inp_ATOM {
                valence: 2,
                chem_bonds_valence: 3,
                el_number: 8,
                charge: 1,
                num_H: 0,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                valence: 1,
                el_number: 11,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                valence: 1,
                el_number: 6,
                ..inp_ATOM::default()
            },
        ];
        charge_one[0].neighbor[0] = 1;
        charge_one[0].neighbor[1] = 2;
        charge_one[0].bond_type[0] = BOND_SINGLE as u8;
        charge_one[0].bond_type[1] = BOND_DOUBLE as u8;
        assert_eq!(n_no_metal_num_bonds(Some(&charge_one), 0), Ok(1));

        assert_eq!(
            n_no_metal_num_bonds(None, 0),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(
            n_no_metal_num_bonds(Some(&atoms), 3),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__util__nnometalbondsvalence__line_1320() {
        let mut atoms = vec![
            inp_ATOM {
                valence: 2,
                chem_bonds_valence: 3,
                el_number: 8,
                charge: 0,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                valence: 1,
                el_number: 11,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                valence: 1,
                el_number: 6,
                ..inp_ATOM::default()
            },
        ];
        atoms[0].neighbor[0] = 1;
        atoms[0].neighbor[1] = 2;
        atoms[0].bond_type[0] = BOND_SINGLE as u8;
        atoms[0].bond_type[1] = BOND_DOUBLE as u8;
        assert_eq!(n_no_metal_bonds_valence(Some(&atoms), 0), Ok(2));

        let mut fallback = atoms.clone();
        fallback[0].bond_type[0] = BOND_TYPE_ALTERN as u8;
        assert_eq!(n_no_metal_bonds_valence(Some(&fallback), 0), Ok(2));

        let mut charge_one = atoms.clone();
        charge_one[0].charge = 1;
        assert_eq!(n_no_metal_bonds_valence(Some(&charge_one), 0), Ok(2));

        let mut no_metal = atoms.clone();
        no_metal[1].el_number = 6;
        assert_eq!(n_no_metal_bonds_valence(Some(&no_metal), 0), Ok(3));

        assert_eq!(
            n_no_metal_bonds_valence(None, 0),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(
            n_no_metal_bonds_valence(Some(&atoms), -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__util__nnometalneighindex__line_1386() {
        let mut atoms = vec![
            inp_ATOM {
                valence: 3,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 11,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 12,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 6,
                ..inp_ATOM::default()
            },
        ];
        atoms[0].neighbor[0] = 1;
        atoms[0].neighbor[1] = 2;
        atoms[0].neighbor[2] = 3;
        assert_eq!(n_no_metal_neigh_index(Some(&atoms), 0), Ok(2));

        atoms[3].el_number = 13;
        assert_eq!(n_no_metal_neigh_index(Some(&atoms), 0), Ok(-1));

        let mut negative = atoms.clone();
        negative[0].valence = -1;
        assert_eq!(n_no_metal_neigh_index(Some(&negative), 0), Ok(-1));

        assert_eq!(
            n_no_metal_neigh_index(None, 0),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(
            n_no_metal_neigh_index(Some(&atoms), -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__util__nnometalotherneighindex__line_1405() {
        let mut atoms = vec![
            inp_ATOM {
                valence: 3,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 11,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 6,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 7,
                ..inp_ATOM::default()
            },
        ];
        atoms[0].neighbor[0] = 1;
        atoms[0].neighbor[1] = 2;
        atoms[0].neighbor[2] = 3;
        assert_eq!(n_no_metal_other_neigh_index(Some(&atoms), 0, 2), Ok(2));
        assert_eq!(n_no_metal_other_neigh_index(Some(&atoms), 0, -1), Ok(1));

        atoms[3].el_number = 12;
        assert_eq!(n_no_metal_other_neigh_index(Some(&atoms), 0, 2), Ok(-1));

        let mut skipped_invalid = atoms.clone();
        skipped_invalid[0].neighbor[1] = 99;
        assert_eq!(
            n_no_metal_other_neigh_index(Some(&skipped_invalid), 0, 99),
            Ok(-1)
        );
        assert_eq!(
            n_no_metal_other_neigh_index(Some(&skipped_invalid), 0, -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__util__nnometalotherneighindex2__line_1426() {
        let mut atoms = vec![
            inp_ATOM {
                valence: 4,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 11,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 12,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 6,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 7,
                ..inp_ATOM::default()
            },
        ];
        atoms[0].neighbor[0] = 1;
        atoms[0].neighbor[1] = 2;
        atoms[0].neighbor[2] = 3;
        atoms[0].neighbor[3] = 4;
        assert_eq!(n_no_metal_other_neigh_index2(Some(&atoms), 0, 1, 2), Ok(2));
        assert_eq!(n_no_metal_other_neigh_index2(Some(&atoms), 0, 1, -1), Ok(2));

        atoms[3].el_number = 12;
        atoms[4].el_number = 7;
        assert_eq!(n_no_metal_other_neigh_index2(Some(&atoms), 0, 1, 2), Ok(3));
        assert_eq!(n_no_metal_other_neigh_index2(Some(&atoms), 0, 1, 3), Ok(3));

        let mut skipped_invalid = atoms.clone();
        skipped_invalid[0].neighbor[2] = 99;
        assert_eq!(
            n_no_metal_other_neigh_index2(Some(&skipped_invalid), 0, 1, 99),
            Ok(3)
        );
        assert_eq!(
            n_no_metal_other_neigh_index2(Some(&skipped_invalid), 0, -1, -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut negative = atoms.clone();
        negative[0].valence = -1;
        assert_eq!(
            n_no_metal_other_neigh_index2(Some(&negative), 0, 1, 2),
            Ok(-1)
        );
        assert_eq!(
            n_no_metal_other_neigh_index2(None, 0, 1, 2),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(
            n_no_metal_other_neigh_index2(Some(&atoms), -1, 1, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__util__get_atomic_mass_from_elnum__line_1007() {
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
            .arg("--atomic-mass-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let records = String::from_utf8(oracle.stdout)
            .expect("official C oracle output must be UTF-8")
            .lines()
            .map(|line| serde_json::from_str::<Value>(line).expect("oracle record must be JSON"))
            .collect::<Vec<_>>();
        assert_eq!(records.len(), 143);
        for official in records {
            assert_eq!(official["operation"], "get_atomic_mass_from_elnum");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let atomic_number = official["input"]["atomic_number"].as_i64().unwrap() as i32;
            let rust = get_atomic_mass_from_elnum(atomic_number)
                .unwrap_or_else(|error| panic!("{case_id}: {error:?}"));
            let expected = official["output"]["result"].as_i64().unwrap() as i32;
            assert_eq!(rust, expected, "{case_id}");
        }
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__util__extract_h_atoms__line_774() {
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
            .arg("--extract-hydrogen-records")
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
        let mut record_count = 0;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["operation"], "extract_H_atoms");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let text = official["input"]["text"]
                .as_str()
                .expect("input text must be text");
            let mut buffer = vec![0x55_i8; 128];
            for (target, source) in buffer.iter_mut().zip(text.bytes()) {
                *target = source as i8;
            }
            buffer[text.len()] = 0;
            let mut isotopes = official["input"]["initial_isotopes"]
                .as_array()
                .expect("initial isotopes must be an array")
                .iter()
                .map(|value| value.as_i64().unwrap() as S_CHAR)
                .collect::<Vec<_>>();
            let hydrogen_count = extract_h_atoms(Some(&mut buffer), Some(&mut isotopes))
                .unwrap_or_else(|error| panic!("{case_id}: {error:?}"));
            let expected = &official["output"];
            assert_eq!(
                hydrogen_count,
                expected["hydrogen_count"].as_i64().unwrap() as i32,
                "{case_id}"
            );
            let expected_isotopes = expected["isotopes"]
                .as_array()
                .expect("isotopes must be an array")
                .iter()
                .map(|value| value.as_i64().unwrap() as S_CHAR)
                .collect::<Vec<_>>();
            assert_eq!(isotopes, expected_isotopes, "{case_id}");
            let nul_offset = buffer
                .iter()
                .position(|byte| *byte == 0)
                .expect("result must remain NUL terminated");
            assert_eq!(
                nul_offset,
                expected["nul_offset"].as_u64().unwrap() as usize,
                "{case_id}"
            );
            let expected_buffer = expected["buffer"]
                .as_array()
                .expect("buffer must be an array")
                .iter()
                .map(|value| value.as_i64().unwrap() as i8)
                .collect::<Vec<_>>();
            assert_eq!(buffer, expected_buffer, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 56);
    }

    #[test]
    fn source_port__util__is_el_a_metal__line_688() {
        for periodic_number in -1..=120 {
            let expected = i32::from(matches!(
                periodic_number,
                3..=4 | 11..=13 | 19..=31 | 37..=51 | 55..=84 | 87..=118
            ));
            assert_eq!(
                is_el_a_metal(periodic_number),
                Ok(expected),
                "periodic number {periodic_number}"
            );
        }

        assert_eq!(is_el_a_metal(-2), Err(SourceHeapError::PointerOutOfBounds));
        assert_eq!(is_el_a_metal(121), Ok(0));
        assert_eq!(
            is_el_a_metal(i32::MIN),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            is_el_a_metal(i32::MAX),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
    }

    #[test]
    fn source_port__util__get_el_type__line_679() {
        const METAL2_PERIODIC_NUMBERS: [i32; 27] = [
            25, 26, 27, 28, 50, 58, 59, 62, 63, 65, 68, 69, 74, 75, 76, 77, 78, 80, 81, 82, 84, 90,
            91, 92, 93, 94, 95,
        ];
        for periodic_number in -1..=120 {
            let expected = if METAL2_PERIODIC_NUMBERS.contains(&periodic_number) {
                3
            } else if matches!(
                periodic_number,
                3..=4 | 11..=13 | 19..=31 | 37..=51 | 55..=84 | 87..=118
            ) {
                1
            } else {
                0
            };
            assert_eq!(
                get_el_type(periodic_number),
                Ok(expected),
                "periodic number {periodic_number}"
            );
        }
        assert_eq!(get_el_type(-2), Err(SourceHeapError::PointerOutOfBounds));
        assert_eq!(get_el_type(121), Ok(0));
        assert_eq!(
            get_el_type(i32::MIN),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            get_el_type(i32::MAX),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__util__is_el_a_metal__line_688() {
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
            .arg("--metal-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let records = String::from_utf8(oracle.stdout)
            .expect("official C oracle output must be UTF-8")
            .lines()
            .map(|line| serde_json::from_str::<Value>(line).expect("oracle record must be JSON"))
            .collect::<Vec<_>>();
        assert_eq!(records.len(), 123);
        for official in records {
            assert_eq!(official["operation"], "is_el_a_metal");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let periodic_number = official["input"]["periodic_number"].as_i64().unwrap() as i32;
            let rust = is_el_a_metal(periodic_number)
                .unwrap_or_else(|error| panic!("{case_id}: {error:?}"));
            let expected = official["output"]["result"].as_i64().unwrap() as i32;
            assert_eq!(rust, expected, "{case_id}");
        }
    }

    #[test]
    fn source_port__util__nbondsvaltometal__line_1100() {
        assert_eq!(
            n_bonds_val_to_metal(None, 0),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(
            n_bonds_val_to_metal(Some(&[]), 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut negative_valence = inp_ATOM::default();
        negative_valence.valence = -1;
        assert_eq!(n_bonds_val_to_metal(Some(&[negative_valence]), 0), Ok(0));

        let mut center = inp_ATOM::default();
        center.valence = 5;
        center.neighbor[..5].copy_from_slice(&[1, 2, 3, 4, 5]);
        center.bond_type[..5].copy_from_slice(&[1, 255, 2, 0, 3]);
        let mut lithium = inp_ATOM::default();
        lithium.el_number = 3;
        let mut carbon = inp_ATOM::default();
        carbon.el_number = 6;
        let mut iron = inp_ATOM::default();
        iron.el_number = 26;
        let mut oganesson = inp_ATOM::default();
        oganesson.el_number = 118;
        let mut oxygen = inp_ATOM::default();
        oxygen.el_number = 8;
        let atoms = [center, lithium, carbon, iron, oganesson, oxygen];
        assert_eq!(n_bonds_val_to_metal(Some(&atoms), 0), Ok(3));

        let mut all_orders = inp_ATOM::default();
        all_orders.valence = 4;
        all_orders.neighbor[..4].copy_from_slice(&[1, 1, 1, 1]);
        all_orders.bond_type[..4].copy_from_slice(&[0, 1, 2, 3]);
        assert_eq!(
            n_bonds_val_to_metal(Some(&[all_orders, atoms[1].clone()]), 0),
            Ok(6)
        );

        for undefined_order in [4_u8, u8::MAX] {
            let mut undefined = inp_ATOM::default();
            undefined.valence = 1;
            undefined.neighbor[0] = 1;
            undefined.bond_type[0] = undefined_order;
            assert_eq!(
                n_bonds_val_to_metal(Some(&[undefined, atoms[3].clone()]), 0),
                Ok(-1)
            );
        }

        let mut ignored_undefined = inp_ATOM::default();
        ignored_undefined.valence = 1;
        ignored_undefined.neighbor[0] = 1;
        ignored_undefined.bond_type[0] = u8::MAX;
        assert_eq!(
            n_bonds_val_to_metal(Some(&[ignored_undefined, atoms[2].clone()]), 0),
            Ok(0)
        );

        let mut invalid_neighbor = inp_ATOM::default();
        invalid_neighbor.valence = 1;
        invalid_neighbor.neighbor[0] = 2;
        assert_eq!(
            n_bonds_val_to_metal(Some(&[invalid_neighbor]), 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut invalid_element = inp_ATOM::default();
        invalid_element.valence = 1;
        invalid_element.neighbor[0] = 1;
        let mut element_121 = inp_ATOM::default();
        element_121.el_number = 121;
        assert_eq!(
            n_bonds_val_to_metal(Some(&[invalid_element, element_121]), 0),
            Ok(0)
        );

        let mut excessive_valence = inp_ATOM::default();
        excessive_valence.valence = 21;
        let hydrogen = inp_ATOM {
            el_number: 1,
            ..inp_ATOM::default()
        };
        assert_eq!(
            n_bonds_val_to_metal(Some(&[excessive_valence, hydrogen]), 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            n_bonds_val_to_metal(Some(&atoms), -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            n_bonds_val_to_metal(Some(&atoms), i32::MAX),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__util__nbondsvaltometal__line_1100() {
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
            .arg("--bonds-to-metal-records")
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
        let mut record_count = 0;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["operation"], "nBondsValToMetal");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let valence = official["input"]["valence"].as_i64().unwrap() as i8;
            let elements = official["input"]["elements"]
                .as_array()
                .expect("elements must be an array");
            let bond_types = official["input"]["bond_types"]
                .as_array()
                .expect("bond types must be an array");
            assert_eq!(elements.len(), bond_types.len(), "{case_id}");
            let mut atoms = vec![inp_ATOM::default(); elements.len() + 1];
            atoms[0].valence = valence;
            for index in 0..elements.len() {
                atoms[0].neighbor[index] = (index + 1) as AT_NUMB;
                atoms[0].bond_type[index] = bond_types[index].as_u64().unwrap() as u8;
                atoms[index + 1].el_number = elements[index].as_u64().unwrap() as u8;
            }
            let rust = n_bonds_val_to_metal(Some(&atoms), 0)
                .unwrap_or_else(|error| panic!("{case_id}: {error:?}"));
            let expected = official["output"]["result"].as_i64().unwrap() as i32;
            assert_eq!(rust, expected, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 31_239);
    }

    #[test]
    fn source_port__util__remove_one_lf__line_1740() {
        fn run_case(input: &[i8], expected: &[i8]) {
            let mut heap = SourceHeap::default();
            let pointer = allocate_source_fixture(&mut heap, input.to_vec());
            assert_eq!(remove_one_lf(&mut heap, pointer), Ok(()));
            assert_eq!(heap.slice(pointer.as_const()).unwrap(), expected);
        }

        let mut heap = SourceHeap::default();
        assert_eq!(remove_one_lf(&mut heap, SourceMutPointer::null()), Ok(()));

        run_case(&[0, 91], &[0, 91]);
        run_case(
            &[b'a' as i8, b'\n' as i8, b'b' as i8, 0, 91],
            &[b'a' as i8, b'\n' as i8, b'b' as i8, 0, 91],
        );
        run_case(&[b'a' as i8, b'\n' as i8, 0, 91], &[b'a' as i8, 0, 0, 91]);
        run_case(
            &[b'a' as i8, b'\r' as i8, b'\n' as i8, 0, 91],
            &[b'a' as i8, 0, 0, 0, 91],
        );
        run_case(
            &[b'a' as i8, b'\n' as i8, b'\n' as i8, 0, 91],
            &[b'a' as i8, b'\n' as i8, 0, 0, 91],
        );
        run_case(
            &[b'a' as i8, b'\r' as i8, 0, 91],
            &[b'a' as i8, b'\r' as i8, 0, 91],
        );

        let unterminated = allocate_source_fixture(&mut heap, vec![b'a' as i8, b'\n' as i8]);
        assert_eq!(
            remove_one_lf(&mut heap, unterminated),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(
            heap.slice(unterminated.as_const()).unwrap(),
            &[b'a' as i8, b'\n' as i8]
        );
    }

    #[test]
    fn source_port__util__mystrncpy__line_1760() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            mystrncpy(
                &mut heap,
                SourceMutPointer::null(),
                SourceConstPointer::null(),
                5
            ),
            Ok(0)
        );

        let source =
            allocate_source_fixture(&mut heap, vec![b'a' as i8, b'b' as i8, b'c' as i8, 0]);
        let target = allocate_source_fixture(&mut heap, vec![99_i8; 6]);
        assert_eq!(mystrncpy(&mut heap, target, source.as_const(), 6), Ok(1));
        assert_eq!(
            heap.slice(target.as_const()).unwrap(),
            &[b'a' as i8, b'b' as i8, b'c' as i8, 0, 0, 0]
        );
        assert_eq!(inchi_free(&mut heap, source), Ok(()));
        assert_eq!(inchi_free(&mut heap, target), Ok(()));

        let source = allocate_source_fixture(
            &mut heap,
            b"abcdef"
                .iter()
                .map(|byte| *byte as i8)
                .chain([0])
                .collect(),
        );
        let target = allocate_source_fixture(&mut heap, vec![99_i8; 4]);
        assert_eq!(mystrncpy(&mut heap, target, source.as_const(), 4), Ok(1));
        assert_eq!(
            heap.slice(target.as_const()).unwrap(),
            &[b'a' as i8, b'b' as i8, b'c' as i8, 0]
        );
        assert_eq!(inchi_free(&mut heap, source), Ok(()));
        assert_eq!(inchi_free(&mut heap, target), Ok(()));

        let overlap = allocate_source_fixture(
            &mut heap,
            b"abcdef"
                .iter()
                .map(|byte| *byte as i8)
                .chain([0, 99, 99])
                .collect(),
        );
        assert_eq!(
            mystrncpy(&mut heap, overlap.offset(2).unwrap(), overlap.as_const(), 5),
            Ok(1)
        );
        assert_eq!(
            heap.slice(overlap.as_const()).unwrap(),
            &[
                b'a' as i8, b'b' as i8, b'a' as i8, b'b' as i8, b'c' as i8, b'd' as i8, 0, 99, 99,
            ]
        );
        assert_eq!(inchi_free(&mut heap, overlap), Ok(()));

        let unterminated = allocate_source_fixture(&mut heap, vec![b'x' as i8]);
        let target = allocate_source_fixture(&mut heap, vec![99_i8; 2]);
        assert_eq!(
            mystrncpy(&mut heap, target, unterminated.as_const(), 2),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(heap.slice(target.as_const()).unwrap(), &[99, 99]);
        assert_eq!(inchi_free(&mut heap, unterminated), Ok(()));
        assert_eq!(inchi_free(&mut heap, target), Ok(()));
    }

    #[test]
    fn source_port__util__inchi_malloc__line_1552() {
        let mut heap = SourceHeap::default();

        let bytes = inchi_malloc(&mut heap, 7).unwrap();
        assert!(!bytes.is_null());
        assert_eq!(heap.slice(bytes.as_const()).unwrap().len(), 7);
        assert_eq!(inchi_free(&mut heap, bytes), Ok(()));
    }

    #[test]
    fn source_port__util__extract_inchi_substring__line_1860() {
        let mut heap = SourceHeap::default();
        let old = allocate_source_fixture(&mut heap, vec![b'o' as i8, 0]);
        let input = allocate_source_fixture(
            &mut heap,
            b"prefix InChI=1S/C2H6/c1-2/h1-2H3?@(),*+.;= trailing"
                .iter()
                .map(|byte| *byte as i8)
                .chain(std::iter::once(0))
                .collect(),
        );
        let mut output = old;
        assert_eq!(
            extract_inchi_substring(&mut heap, &mut output, input.as_const(), 128),
            Ok(())
        );
        let extracted = heap.slice(output.as_const()).unwrap();
        let length = extracted.iter().position(|byte| *byte == 0).unwrap();
        assert_eq!(
            extracted[..length]
                .iter()
                .map(|byte| *byte as u8)
                .collect::<Vec<_>>(),
            b"InChI=1S/C2H6/c1-2/h1-2H3?@(),*+.;="
        );
        assert_eq!(heap.slice(old.as_const()).unwrap(), &[b'o' as i8, 0]);

        let mut truncated = SourceMutPointer::null();
        assert_eq!(
            extract_inchi_substring(&mut heap, &mut truncated, input.as_const(), 4),
            Ok(())
        );
        assert_eq!(
            heap.slice(truncated.as_const()).unwrap(),
            &[b'I' as i8, b'n' as i8, b'C' as i8, b'h' as i8, 0]
        );

        let no_match = allocate_source_fixture(
            &mut heap,
            vec![
                b'i' as i8, b'n' as i8, b'c' as i8, b'h' as i8, b'i' as i8, b'=' as i8, 0,
            ],
        );
        let mut unchanged = old;
        assert_eq!(
            extract_inchi_substring(&mut heap, &mut unchanged, no_match.as_const(), 6),
            Ok(())
        );
        assert!(unchanged.is_null());
        let mut null_string = old;
        assert_eq!(
            extract_inchi_substring(&mut heap, &mut null_string, SourceConstPointer::null(), 0,),
            Ok(())
        );
        assert!(null_string.is_null());

        heap.fail_after_allocations(0);
        let mut allocation_failure = old;
        assert_eq!(
            extract_inchi_substring(&mut heap, &mut allocation_failure, input.as_const(), 128,),
            Err(SourceHeapError::AllocationFailed)
        );
        assert!(allocation_failure.is_null());
    }

    #[test]
    fn source_port__util__extract_auxinfo_substring__line_1921() {
        fn allocated(
            heap: &mut SourceHeap,
            bytes: impl IntoIterator<Item = i8>,
        ) -> SourceMutPointer<i8> {
            allocate_source_fixture(heap, bytes.into_iter().collect())
        }

        for whitespace in [b' ', b'\t', b'\n', 0x0b, 0x0c, b'\r'] {
            let mut heap = SourceHeap::default();
            let old = allocated(&mut heap, [b'o' as i8, 0]);
            let input = allocated(
                &mut heap,
                b"prefix:AuxInfo=1/0/N:1"
                    .iter()
                    .copied()
                    .chain([whitespace, b'X', 0])
                    .map(|byte| byte as i8),
            );
            let mut output = old;
            assert_eq!(
                extract_auxinfo_substring(&mut heap, &mut output, input.as_const(), 128),
                Ok(())
            );
            assert_eq!(
                heap.slice(output.as_const()).unwrap(),
                b"AuxInfo=1/0/N:1\0"
                    .iter()
                    .map(|byte| *byte as i8)
                    .collect::<Vec<_>>()
            );
            assert_eq!(heap.slice(old.as_const()).unwrap(), &[b'o' as i8, 0]);
        }

        let mut heap = SourceHeap::default();
        let input = allocated(
            &mut heap,
            b"xxAuxInfo=1/1/N:1\xA0Z \0".iter().map(|byte| *byte as i8),
        );
        let mut truncated = SourceMutPointer::null();
        assert_eq!(
            extract_auxinfo_substring(&mut heap, &mut truncated, input.as_const(), 4),
            Ok(())
        );
        assert_eq!(
            heap.slice(truncated.as_const()).unwrap(),
            &[b'A' as i8, b'u' as i8, b'x' as i8, b'I' as i8, 0]
        );
        let mut full = SourceMutPointer::null();
        assert_eq!(
            extract_auxinfo_substring(&mut heap, &mut full, input.as_const(), 128),
            Ok(())
        );
        assert_eq!(
            heap.slice(full.as_const()).unwrap(),
            b"AuxInfo=1/1/N:1\xA0Z\0"
                .iter()
                .map(|byte| *byte as i8)
                .collect::<Vec<_>>()
        );

        let mut zero_length = SourceMutPointer::null();
        assert_eq!(
            extract_auxinfo_substring(&mut heap, &mut zero_length, input.as_const(), 0),
            Ok(())
        );
        assert_eq!(heap.slice(zero_length.as_const()).unwrap(), &[0]);

        let embedded_nul = allocated(&mut heap, b"AuxInfo=A\0Q \0".iter().map(|byte| *byte as i8));
        let mut through_nul = SourceMutPointer::null();
        assert_eq!(
            extract_auxinfo_substring(&mut heap, &mut through_nul, embedded_nul.as_const(), 32),
            Ok(())
        );
        assert_eq!(
            heap.slice(through_nul.as_const()).unwrap(),
            b"AuxInfo=A\0Q\0"
                .iter()
                .map(|byte| *byte as i8)
                .collect::<Vec<_>>()
        );

        let old = allocated(
            &mut heap,
            [b'l' as i8, b'e' as i8, b'a' as i8, b'k' as i8, 0],
        );
        let empty = allocated(&mut heap, [0]);
        let mut output = old;
        assert_eq!(
            extract_auxinfo_substring(&mut heap, &mut output, empty.as_const(), 7),
            Ok(())
        );
        assert!(output.is_null());
        assert!(heap.slice(old.as_const()).is_ok());

        let no_match = allocated(&mut heap, b"auxinfo=1\0".iter().map(|byte| *byte as i8));
        output = old;
        assert_eq!(
            extract_auxinfo_substring(&mut heap, &mut output, no_match.as_const(), 11),
            Ok(())
        );
        assert!(output.is_null());
        output = old;
        assert_eq!(
            extract_auxinfo_substring(&mut heap, &mut output, SourceConstPointer::null(), u64::MAX,),
            Ok(())
        );
        assert!(output.is_null());

        let unterminated = allocated(&mut heap, b"AuxInfo=1".iter().map(|byte| *byte as i8));
        output = old;
        assert_eq!(
            extract_auxinfo_substring(&mut heap, &mut output, unterminated.as_const(), 2),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert!(output.is_null());

        let bounded = allocated(&mut heap, b"AuxInfo=1\0".iter().map(|byte| *byte as i8));
        output = old;
        assert_eq!(
            extract_auxinfo_substring(&mut heap, &mut output, bounded.as_const(), 12),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert!(output.is_null());

        heap.fail_after_allocations(0);
        output = old;
        assert_eq!(
            extract_auxinfo_substring(&mut heap, &mut output, input.as_const(), 4),
            Err(SourceHeapError::AllocationFailed)
        );
        assert!(output.is_null());
    }

    #[test]
    fn source_port__util__lrtrim__line_1804() {
        let mut heap = SourceHeap::default();

        let mut null_length = 9;
        assert_eq!(
            lrtrim(&mut heap, SourceMutPointer::null(), Some(&mut null_length)),
            Ok(SourceMutPointer::null())
        );
        assert_eq!(null_length, 0);

        let text = allocate_source_fixture(
            &mut heap,
            vec![
                b' ' as i8,
                b'\t' as i8,
                b'A' as i8,
                b'b' as i8,
                b' ' as i8,
                b'c' as i8,
                b'\r' as i8,
                b'\n' as i8,
                0,
                77,
            ],
        );
        let mut length = -1;
        assert_eq!(lrtrim(&mut heap, text, Some(&mut length)), Ok(text));
        assert_eq!(length, 4);
        assert_eq!(
            &heap.slice(text.as_const()).unwrap()[..5],
            &[b'A' as i8, b'b' as i8, b' ' as i8, b'c' as i8, 0]
        );
        assert_eq!(inchi_free(&mut heap, text), Ok(()));

        let all_space =
            allocate_source_fixture(&mut heap, vec![b' ' as i8, b'\t' as i8, b'\n' as i8, 0]);
        let mut length = -1;
        assert_eq!(
            lrtrim(&mut heap, all_space, Some(&mut length)),
            Ok(all_space)
        );
        assert_eq!(length, 0);
        assert_eq!(heap.slice(all_space.as_const()).unwrap()[0], 0);
        assert_eq!(inchi_free(&mut heap, all_space), Ok(()));

        let non_ascii =
            allocate_source_fixture(&mut heap, vec![b' ' as i8, 0xa0_u8 as i8, b' ' as i8, 0]);
        assert_eq!(lrtrim(&mut heap, non_ascii, None), Ok(non_ascii));
        assert_eq!(
            &heap.slice(non_ascii.as_const()).unwrap()[..2],
            &[0xa0_u8 as i8, 0]
        );
        assert_eq!(inchi_free(&mut heap, non_ascii), Ok(()));

        let unterminated = allocate_source_fixture(&mut heap, vec![b'x' as i8]);
        let mut unchanged_length = 17;
        assert_eq!(
            lrtrim(&mut heap, unterminated, Some(&mut unchanged_length)),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(unchanged_length, 17);
        assert_eq!(inchi_free(&mut heap, unterminated), Ok(()));
    }

    #[test]
    fn source_port__util__normalize_string__line_1589() {
        let mut heap = SourceHeap::default();

        let mixed = allocate_source_fixture(
            &mut heap,
            vec![
                b' ' as i8,
                b' ' as i8,
                b'A' as i8,
                b'\t' as i8,
                b'\t' as i8,
                b'B' as i8,
                b'\r' as i8,
                b'\n' as i8,
                0,
                b'X' as i8,
                b'Y' as i8,
            ],
        );
        assert_eq!(normalize_string(&mut heap, mixed), Ok(3));
        assert_eq!(
            heap.slice(mixed.as_const()).unwrap(),
            &[
                b'A' as i8,
                b' ' as i8,
                b'B' as i8,
                0,
                b' ' as i8,
                0,
                0,
                b'\n' as i8,
                0,
                b'X' as i8,
                b'Y' as i8,
            ]
        );

        let all_space = allocate_source_fixture(
            &mut heap,
            vec![
                b'\t' as i8,
                b'\n' as i8,
                0x0b,
                0x0c,
                b'\r' as i8,
                b' ' as i8,
                0,
                77,
            ],
        );
        assert_eq!(normalize_string(&mut heap, all_space), Ok(0));
        assert_eq!(
            heap.slice(all_space.as_const()).unwrap(),
            &[
                0, b' ' as i8, b' ' as i8, b' ' as i8, b' ' as i8, b' ' as i8, 0, 77
            ]
        );

        let empty = allocate_source_fixture(&mut heap, vec![0_i8, 91]);
        assert_eq!(normalize_string(&mut heap, empty), Ok(0));
        assert_eq!(heap.slice(empty.as_const()).unwrap(), &[0, 91]);

        let non_ascii = allocate_source_fixture(
            &mut heap,
            vec![
                b' ' as i8,
                0xa0_u8 as i8,
                b' ' as i8,
                b'C' as i8,
                b' ' as i8,
                0,
            ],
        );
        assert_eq!(normalize_string(&mut heap, non_ascii), Ok(3));
        assert_eq!(
            &heap.slice(non_ascii.as_const()).unwrap()[..4],
            &[0xa0_u8 as i8, b' ' as i8, b'C' as i8, 0]
        );

        let interior = allocate_source_fixture(
            &mut heap,
            vec![99_i8, b' ' as i8, b'D' as i8, b' ' as i8, 0, 88],
        );
        let interior_start = interior.offset(1).unwrap();
        assert_eq!(normalize_string(&mut heap, interior_start), Ok(1));
        assert_eq!(
            heap.slice(interior.as_const()).unwrap(),
            &[99, b'D' as i8, 0, 0, 0, 88]
        );

        assert_eq!(
            normalize_string(&mut heap, SourceMutPointer::null()),
            Err(SourceHeapError::NullPointer)
        );
        let unterminated = allocate_source_fixture(&mut heap, vec![b' ' as i8, b'Q' as i8]);
        assert_eq!(
            normalize_string(&mut heap, unterminated),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(
            heap.slice(unterminated.as_const()).unwrap(),
            &[b' ' as i8, b'Q' as i8]
        );
    }

    #[test]
    fn source_port__util__is_matching_any_delim__line_1710() {
        let whitespace = [
            b' ' as i8,
            b'\t' as i8,
            b'\n' as i8,
            0x0b,
            0x0c,
            b'\r' as i8,
            0,
        ];
        for character in &whitespace[..6] {
            assert_eq!(is_matching_any_delim(*character, Some(&whitespace)), Ok(1));
        }
        assert_eq!(is_matching_any_delim(b'X' as i8, Some(&whitespace)), Ok(0));
        assert_eq!(is_matching_any_delim(0, Some(&whitespace)), Ok(0));
        assert_eq!(
            is_matching_any_delim(b',' as i8, Some(&[b';' as i8, b',' as i8, b',' as i8, 0])),
            Ok(1)
        );
        assert_eq!(
            is_matching_any_delim(0xff_u8 as i8, Some(&[0xff_u8 as i8, 0])),
            Ok(0)
        );
        assert_eq!(
            is_matching_any_delim(b'A' as i8, None),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(
            is_matching_any_delim(b'A' as i8, Some(&[b'B' as i8, b'C' as i8])),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(
            is_matching_any_delim(b'A' as i8, Some(&[b'A' as i8])),
            Ok(1)
        );
    }

    #[test]
    fn source_port__util__read_upto_delim__line_1658() {
        let mut heap = SourceHeap::default();
        let delimiters = [b',' as i8, b';' as i8, 0];
        let input = allocate_source_fixture(
            &mut heap,
            b" \talpha,beta\0".iter().map(|byte| *byte as i8).collect(),
        );
        let field = allocate_source_fixture(&mut heap, vec![0x55_i8; 16]);
        let mut cursor = input;
        assert_eq!(
            read_upto_delim(&mut heap, &mut cursor, field, 16, Some(&delimiters)),
            Ok(5)
        );
        assert_eq!(cursor, input.offset(7).unwrap());
        assert_eq!(
            &heap.slice(field.as_const()).unwrap()[..8],
            &[
                b'a' as i8, b'l' as i8, b'p' as i8, b'h' as i8, b'a' as i8, 0, 0, 0x55
            ]
        );

        heap.slice_mut(field).unwrap().fill(0x44);
        assert_eq!(
            read_upto_delim(&mut heap, &mut cursor, field, 1, Some(&delimiters)),
            Ok(0)
        );
        assert_eq!(cursor, input.offset(7).unwrap());
        assert_eq!(&heap.slice(field.as_const()).unwrap()[..3], &[0, 0, 0x44]);

        cursor = cursor.offset(1).unwrap();
        heap.slice_mut(field).unwrap().fill(0x33);
        assert_eq!(
            read_upto_delim(&mut heap, &mut cursor, field, 5, Some(&delimiters)),
            Ok(4)
        );
        assert!(cursor.is_null());
        assert_eq!(
            &heap.slice(field.as_const()).unwrap()[..7],
            &[b'b' as i8, b'e' as i8, b't' as i8, b'a' as i8, 0, 0, 0x33]
        );

        let too_long = allocate_source_fixture(
            &mut heap,
            b"  abc,\0".iter().map(|byte| *byte as i8).collect(),
        );
        cursor = too_long;
        heap.slice_mut(field).unwrap().fill(0x22);
        assert_eq!(
            read_upto_delim(&mut heap, &mut cursor, field, 3, Some(&delimiters)),
            Ok(-1)
        );
        assert_eq!(cursor, too_long);
        assert!(
            heap.slice(field.as_const())
                .unwrap()
                .iter()
                .all(|byte| *byte == 0x22)
        );
        assert_eq!(
            read_upto_delim(&mut heap, &mut cursor, field, -1, Some(&delimiters)),
            Ok(-1)
        );
        assert_eq!(cursor, too_long);

        let spaces =
            allocate_source_fixture(&mut heap, vec![b' ' as i8, b'\t' as i8, b'\n' as i8, 0]);
        cursor = spaces;
        assert_eq!(
            read_upto_delim(&mut heap, &mut cursor, field, 1, None),
            Ok(0)
        );
        assert!(cursor.is_null());

        cursor = input.offset(2).unwrap();
        assert_eq!(
            read_upto_delim(&mut heap, &mut cursor, field, 16, None),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(cursor, input.offset(2).unwrap());

        cursor = SourceMutPointer::null();
        assert_eq!(
            read_upto_delim(&mut heap, &mut cursor, field, 16, Some(&delimiters)),
            Ok(-1)
        );
        let empty = allocate_source_fixture(&mut heap, vec![0_i8]);
        cursor = empty;
        assert_eq!(
            read_upto_delim(
                &mut heap,
                &mut cursor,
                SourceMutPointer::null(),
                1,
                Some(&delimiters)
            ),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(cursor, empty);
    }

    #[test]
    fn source_port__util__remove_trailing_spaces__line_1728() {
        let mut heap = SourceHeap::default();
        let mixed = allocate_source_fixture(
            &mut heap,
            vec![
                b'A' as i8,
                b' ' as i8,
                b'B' as i8,
                b'\t' as i8,
                b'\n' as i8,
                0x0b,
                0x0c,
                b'\r' as i8,
                b' ' as i8,
                0,
                77,
            ],
        );
        assert_eq!(remove_trailing_spaces(&mut heap, mixed), Ok(()));
        assert_eq!(
            heap.slice(mixed.as_const()).unwrap(),
            &[
                b'A' as i8,
                b' ' as i8,
                b'B' as i8,
                0,
                b'\n' as i8,
                0x0b,
                0x0c,
                b'\r' as i8,
                b' ' as i8,
                0,
                77,
            ]
        );
        let all_space =
            allocate_source_fixture(&mut heap, vec![b' ' as i8, b'\t' as i8, b'\n' as i8, 0, 88]);
        assert_eq!(remove_trailing_spaces(&mut heap, all_space), Ok(()));
        assert_eq!(
            heap.slice(all_space.as_const()).unwrap(),
            &[0, b'\t' as i8, b'\n' as i8, 0, 88]
        );
        let non_ascii =
            allocate_source_fixture(&mut heap, vec![b'X' as i8, 0xa0_u8 as i8, b' ' as i8, 0]);
        assert_eq!(remove_trailing_spaces(&mut heap, non_ascii), Ok(()));
        assert_eq!(
            heap.slice(non_ascii.as_const()).unwrap(),
            &[b'X' as i8, 0xa0_u8 as i8, 0, 0]
        );
        let empty = allocate_source_fixture(&mut heap, vec![0_i8, 99]);
        assert_eq!(remove_trailing_spaces(&mut heap, empty), Ok(()));
        assert_eq!(heap.slice(empty.as_const()).unwrap(), &[0, 99]);
        assert_eq!(
            remove_trailing_spaces(&mut heap, SourceMutPointer::null()),
            Err(SourceHeapError::NullPointer)
        );
        let unterminated = allocate_source_fixture(&mut heap, vec![b'X' as i8, b' ' as i8]);
        assert_eq!(
            remove_trailing_spaces(&mut heap, unterminated),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(
            heap.slice(unterminated.as_const()).unwrap(),
            &[b'X' as i8, b' ' as i8]
        );
    }

    #[test]
    fn source_port__util__inchi_calloc__line_1561() {
        let mut heap = SourceHeap::default();

        let integers = inchi_calloc::<i32>(&mut heap, 4, 4).unwrap();
        assert_source_slice_eq("zeroed-i32-allocation", &heap, integers.as_const(), &[0; 4]);
        assert_eq!(inchi_free(&mut heap, integers), Ok(()));

        assert_eq!(
            inchi_calloc::<u8>(&mut heap, u64::MAX, 2),
            Err(SourceHeapError::AllocationSizeOverflow)
        );
    }

    #[test]
    fn source_port__util__inchi_free__line_1570() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            inchi_free(&mut heap, SourceMutPointer::<u8>::null()),
            Ok(())
        );

        let allocation = allocate_source_fixture(&mut heap, vec![3_u8, 1, 4]);
        let alias = allocation.as_const();
        assert_source_slice_eq("before-free", &heap, alias, &[3, 1, 4]);
        assert_eq!(inchi_free(&mut heap, allocation), Ok(()));
        assert_eq!(heap.slice(alias), Err(SourceHeapError::MissingAllocation));

        let allocation = allocate_source_fixture(&mut heap, vec![10_i32, 20]);
        let interior = allocation.offset(1).unwrap();
        assert_eq!(
            inchi_free(&mut heap, interior),
            Err(SourceHeapError::FreeOfInteriorPointer)
        );
        assert_source_slice_eq(
            "rejected-interior-free-preserves-allocation",
            &heap,
            allocation.as_const(),
            &[10, 20],
        );
        assert_eq!(inchi_free(&mut heap, allocation), Ok(()));
    }
}
