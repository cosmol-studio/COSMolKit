use crate::source::base::ichi_io::{
    inchi_ios_close, inchi_ios_init, inchi_strbuf_addline, inchi_strbuf_close, inchi_strbuf_reset, source_vformat,
};
use crate::source::base::ichierr::AddErrorMessage;
use crate::source::base::ichiprt2::{inchi_strtod, inchi_strtol};
use crate::source::base::mol_fmt4::{NumLists_Alloc, NumLists_Append, NumLists_Free};
use crate::source::base::util::{
    dotify_non_printable_chars, get_atomic_mass, inchi_calloc, inchi_free, mystrncpy, mystrncpy_slice, read_upto_delim,
    remove_one_lf, remove_trailing_spaces,
};
use crate::source_types::{
    ATOM_EL_LEN, INCHI_IOS_STRING, INCHI_IOS_TYPE_STRING, INCHI_IOSTREAM, MOL_FMT_CHAR_INT_DATA, MOL_FMT_CTAB,
    MOL_FMT_DOUBLE_DATA, MOL_FMT_FLOAT_DATA, MOL_FMT_INT_DATA, MOL_FMT_LONG_INT_DATA, MOL_FMT_SHORT_INT_DATA,
    MOL_FMT_STRING_DATA, MOL_FMT_V3000_MAXFIELDLEN, MOL_FMT_V3000_STEABS, MOL_FMT_V3000_STENON, MOL_FMT_V3000_STERAC,
    MOL_FMT_V3000_STEREL, MOL_FMT_v3000, NUM_LISTS, SD_FMT_END_OF_DATA, SourceConstPointer, SourceFormatArgument,
    SourceHeap, SourceHeapError, SourceMutPointer, SourceVaList,
};

fn mol_fmt3_c_string_eq(value: &[i8], expected: &[u8]) -> bool {
    let length = value.iter().position(|byte| *byte == 0).unwrap_or(value.len());
    length == expected.len()
        && value[..length]
            .iter()
            .zip(expected)
            .all(|(left, right)| *left as u8 == *right)
}

fn mol_fmt3_add_message(target: Option<&mut [i8]>, message: &[u8]) -> Result<i32, SourceHeapError> {
    let terminated = message
        .iter()
        .copied()
        .chain(std::iter::once(0))
        .map(|byte| byte as i8)
        .collect::<Vec<_>>();
    AddErrorMessage(target, Some(&terminated))
}

fn mol_fmt3_treat_error(
    err: &mut i32,
    new_err: i32,
    target: Option<&mut [i8]>,
    message: &[u8],
) -> Result<(), SourceHeapError> {
    if *err == 0 && new_err != 0 {
        *err = new_err;
    }
    let _ = mol_fmt3_add_message(target, message)?;
    Ok(())
}

pub(crate) enum MolfileV3000FieldData<'a> {
    String(&'a mut [i8]),
    Char(&'a mut i8),
    Short(&'a mut i16),
    Long(&'a mut i64),
    Int(&'a mut i32),
    Double(&'a mut f64),
    Float(&'a mut f32),
    None,
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileV3000ReadField(
    heap: &mut SourceHeap,
    data: MolfileV3000FieldData<'_>,
    data_type: i32,
    line_ptr: &mut SourceMutPointer<i8>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:257 MolfileV3000ReadField
    // INCHI✔️❌: complete source frame follows verbatim; modeled stack storage and checked typed output access add overhead.
    /*
int MolfileV3000ReadField(void *data,
                          int data_type,
                          char **line_ptr)
{
    int nread = 0;
    char field[MOL_FMT_V3000_MAXFIELDLEN];
    const int max_field_len = sizeof(field);
    long ldata = 0L;
    double ddata = 0.0;
    char *p_end;

    memset(field, 0, max_field_len); /* djb-rwth: memset_s C11/Annex K variant? */

    nread = read_upto_delim(line_ptr, field, max_field_len, " \t\n\v\f\r");

    switch (data_type)
    {
    case MOL_FMT_STRING_DATA:
    {
        if (nread && (nread <= ATOM_EL_LEN)) /* djb-rwth: fixing GHI #133 -- updated 28/09/2025 */
        {
            mystrncpy((char *)data, field, nread + 1);
        }
        else
        {
            ((char *)data)[0] = '\0';
        }
    }
    break;

    case MOL_FMT_CHAR_INT_DATA:
    case MOL_FMT_SHORT_INT_DATA:
    case MOL_FMT_LONG_INT_DATA:
    case MOL_FMT_INT_DATA:
    {
        /* assume that field ends at first non-digit */
        ldata = strtol(field, &p_end, 10);

        if (p_end == field)
        {
            nread = 0;
        }

        if (data_type == MOL_FMT_LONG_INT_DATA)
        {
            if (LONG_MIN < ldata && ldata < LONG_MAX)
            {
                *(long *)data = (long)ldata;
            }
            else
            {
                *(long *)data = 0L;
                nread = -1;
            }
        }
        else if (data_type == MOL_FMT_INT_DATA)
        {
            if (INT_MIN <= ldata && ldata <= INT_MAX) /* djb-rwth: addressing coverity ID #499496/499553 -- ldata check seems to be necessary */
            {
                *(int *)data = (int)ldata;
            }
            else
            {
                *(int *)data = (int)0;
                nread = -1;
            }
        }

        else if (data_type == MOL_FMT_CHAR_INT_DATA)
        {
            if (SCHAR_MIN <= ldata && ldata <= SCHAR_MAX)
            {
                *(S_CHAR *)data = (S_CHAR)ldata;
            }
            else
            {
                *(S_CHAR *)data = (S_CHAR)0;
                nread = -1;
            }
        }

        else if (data_type == MOL_FMT_SHORT_INT_DATA)
        {
            if (SHRT_MIN <= ldata && ldata <= SHRT_MAX)
            {
                *(S_SHORT *)data = (S_SHORT)ldata;
            }
            else
            {
                *(S_SHORT *)data = (S_SHORT)0;
                nread = -1;
            }
        }
        else
        {
            nread = -1;
        }
    }
    break; /* INT's */

    case MOL_FMT_DOUBLE_DATA:
    case MOL_FMT_FLOAT_DATA:
    {
        /* assume that field ends at first non-digit */
        ddata = strtod(field, &p_end);

        if (p_end == field)
        {
            nread = 0;
        }

        if (data_type == MOL_FMT_DOUBLE_DATA)
        {
            if (ddata != HUGE_VAL && /*ldata*/ ddata != -HUGE_VAL)
            {
                *(double *)data = ddata;
            }
            else
            {
                *(double *)data = 0.0;
                nread = -1;
            }
        }
        else if (data_type == MOL_FMT_FLOAT_DATA)
        {
            if (fabs(ddata) <= (double)FLT_MIN)
            {
                *(float *)data = 0.0;
            }
            else if (fabs(ddata) >= (double)FLT_MAX)
            {
                *(float *)data = 0.0;
                nread = -1;
            }
        }
        else
        {
            *(float *)data = (float)ddata; /* djb-rwth: addressing coverity ID #499519 -- probably never reached */
        }
    }
    break; /* REAL's */

    default:
        nread = -1;
    }

    return nread;
}
    */
    // END INCHI C FUNCTION: MolfileV3000ReadField
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadField
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; long is signed 64-bit, int 32-bit, short 16-bit, char 8-bit.
    // INCHI✔️❌: ATOM_EL_LEN is 6 and MOL_FMT_V3000_MAXFIELDLEN is 4096.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadField

    let field = heap.allocate_model_storage(vec![0_i8; MOL_FMT_V3000_MAXFIELDLEN as usize])?;
    let result = (|| -> Result<i32, SourceHeapError> {
        let delimiters = [b' ' as i8, b'\t' as i8, b'\n' as i8, 0x0b, 0x0c, b'\r' as i8, 0];
        let mut nread = read_upto_delim(
            heap,
            line_ptr,
            field,
            MOL_FMT_V3000_MAXFIELDLEN as i32,
            Some(&delimiters),
        )?;
        if data_type == i32::from(MOL_FMT_STRING_DATA) {
            let MolfileV3000FieldData::String(output) = data else {
                return Err(SourceHeapError::UnsupportedSourceBehavior);
            };
            if nread != 0 && nread <= ATOM_EL_LEN as i32 {
                let source = heap.slice(field.as_const())?.to_vec();
                let _ = mystrncpy_slice(Some(output), Some(&source), nread.wrapping_add(1) as u32)?;
            } else {
                *output.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
            }
            return Ok(nread);
        }

        let integer = data_type == i32::from(MOL_FMT_CHAR_INT_DATA)
            || data_type == i32::from(MOL_FMT_SHORT_INT_DATA)
            || data_type == i32::from(MOL_FMT_LONG_INT_DATA)
            || data_type == i32::from(MOL_FMT_INT_DATA);
        if integer {
            let mut end = SourceConstPointer::null();
            let value = inchi_strtol(heap, field.as_const(), Some(&mut end), 10)?;
            if end == field.as_const() {
                nread = 0;
            }
            match (data_type, data) {
                (kind, MolfileV3000FieldData::Long(output)) if kind == i32::from(MOL_FMT_LONG_INT_DATA) => {
                    if value > i64::MIN && value < i64::MAX { *output = value; } else { *output = 0; nread = -1; }
                }
                (kind, MolfileV3000FieldData::Int(output)) if kind == i32::from(MOL_FMT_INT_DATA) => {
                    if let Ok(value) = i32::try_from(value) { *output = value; } else { *output = 0; nread = -1; }
                }
                (kind, MolfileV3000FieldData::Char(output)) if kind == i32::from(MOL_FMT_CHAR_INT_DATA) => {
                    if let Ok(value) = i8::try_from(value) { *output = value; } else { *output = 0; nread = -1; }
                }
                (kind, MolfileV3000FieldData::Short(output)) if kind == i32::from(MOL_FMT_SHORT_INT_DATA) => {
                    if let Ok(value) = i16::try_from(value) { *output = value; } else { *output = 0; nread = -1; }
                }
                _ => nread = -1,
            }
            return Ok(nread);
        }

        if data_type == i32::from(MOL_FMT_DOUBLE_DATA) || data_type == i32::from(MOL_FMT_FLOAT_DATA) {
            let mut end = SourceConstPointer::null();
            let value = inchi_strtod(heap, field.as_const(), Some(&mut end))?;
            if end == field.as_const() {
                nread = 0;
            }
            match (data_type, data) {
                (kind, MolfileV3000FieldData::Double(output)) if kind == i32::from(MOL_FMT_DOUBLE_DATA) => {
                    if value.is_infinite() { *output = 0.0; nread = -1; } else { *output = value; }
                }
                (kind, MolfileV3000FieldData::Float(output)) if kind == i32::from(MOL_FMT_FLOAT_DATA) => {
                    let magnitude = value.abs();
                    if magnitude <= f64::from(f32::MIN_POSITIVE) { *output = 0.0; }
                    else if magnitude >= f64::from(f32::MAX) { *output = 0.0; nread = -1; }
                    // The active source has no assignment for an ordinary finite float.
                }
                _ => nread = -1,
            }
            return Ok(nread);
        }
        let _ = data;
        Ok(-1)
    })();
    let cleanup = heap.free(field);
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn get_V3000_input_line_to_strbuf(
    heap: &mut SourceHeap,
    buf: &mut INCHI_IOS_STRING,
    mut inp_stream: Option<&mut INCHI_IOSTREAM>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1889 get_V3000_input_line_to_strbuf
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap memmove and delegated buffer operations add overhead.
    /*
int get_V3000_input_line_to_strbuf(INCHI_IOS_STRING *buf,
                                   INCHI_IOSTREAM *inp_stream)
{
    const int prefix_len = 7; /* "M  V30 " */
    int old_used, crlf2lf = 1, preserve_lf = 0;

    inchi_strbuf_reset(buf);

    old_used = buf->nUsedLength;
    while (1)
    {
        inchi_strbuf_addline(buf, inp_stream, crlf2lf, preserve_lf);

        if (buf->nUsedLength - old_used < 8)
        {
            return -1;
        }
        if (strncmp(buf->pStr + old_used, "M  V30 ", prefix_len))
        {
            return -1;
        }

        memmove((void *)(buf->pStr + old_used), (void *)(buf->pStr + old_used + prefix_len), (long long)buf->nUsedLength - (long long)old_used - (long long)prefix_len + 1); /* djb-rwth: cast operators added */ /* ricrogz: fixing memory overflow error */
        buf->nUsedLength -= prefix_len;

        if (buf->pStr[buf->nUsedLength - 1] != '-')
        {
            break;
        }
        buf->pStr[--buf->nUsedLength] = '\0';

        old_used = buf->nUsedLength;
    }

    remove_trailing_spaces(buf->pStr);
    buf->nUsedLength = strlen(buf->pStr);

    return buf->nUsedLength;
}
    */
    // END INCHI C FUNCTION: get_V3000_input_line_to_strbuf
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: get_V3000_input_line_to_strbuf
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; prefix length is exactly seven bytes.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: get_V3000_input_line_to_strbuf

    inchi_strbuf_reset(heap, Some(&mut *buf))?;
    let mut old_used = buf.nUsedLength;
    loop {
        let _ = inchi_strbuf_addline(heap, &mut *buf, inp_stream.as_deref_mut(), 1, 0)?;
        if buf.nUsedLength.wrapping_sub(old_used) < 8 {
            return Ok(-1);
        }
        let old = usize::try_from(old_used).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let used = usize::try_from(buf.nUsedLength).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let bytes = heap.slice(buf.pStr.as_const())?;
        if bytes.get(old..old + 7) != Some(&b"M  V30 ".map(|byte| byte as i8)) {
            return Ok(-1);
        }
        heap.slice_mut(buf.pStr)?.copy_within(old + 7..=used, old);
        buf.nUsedLength = buf.nUsedLength.wrapping_sub(7);
        let last = usize::try_from(buf.nUsedLength.wrapping_sub(1))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if heap.slice(buf.pStr.as_const())?[last] != b'-' as i8 {
            break;
        }
        buf.nUsedLength = buf.nUsedLength.wrapping_sub(1);
        let end = usize::try_from(buf.nUsedLength).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        heap.slice_mut(buf.pStr)?[end] = 0;
        old_used = buf.nUsedLength;
    }
    remove_trailing_spaces(heap, buf.pStr)?;
    let length = heap
        .slice(buf.pStr.as_const())?
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    buf.nUsedLength = i32::try_from(length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    Ok(buf.nUsedLength)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileV3000ReadKeyword(
    heap: &mut SourceHeap,
    key: SourceMutPointer<i8>,
    line_ptr: &mut SourceMutPointer<i8>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:410 MolfileV3000ReadKeyword
    // INCHI✔️❌: complete source frame follows verbatim; modeled stack field and checked SourceHeap access add overhead.
    /*
int MolfileV3000ReadKeyword(char *key,
                            char **line_ptr)
{
    int nread = 0;
    char field[MOL_FMT_V3000_MAXFIELDLEN];
    const int max_field_len = sizeof(field);

    memset(field, 0, max_field_len); /* djb-rwth: memset_s C11/Annex K variant? */

    nread = read_upto_delim(line_ptr, field, max_field_len, "= \t\n\v\f\r");

    if (nread)
    {
        mystrncpy(key, field, nread + 1);
        /* consume '=' sign if present */
        if (*line_ptr)
        {
            if (*line_ptr[0] == '=')
            {
                *line_ptr = *line_ptr + 1;
            }
        }
    }
    else
    {
        key[0] = '\0';
    }

    return nread;
}
    */
    // END INCHI C FUNCTION: MolfileV3000ReadKeyword
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadKeyword
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; MOL_FMT_V3000_MAXFIELDLEN is 4096.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadKeyword

    let field = heap.allocate_model_storage(vec![0_i8; MOL_FMT_V3000_MAXFIELDLEN as usize])?;
    let result = (|| -> Result<i32, SourceHeapError> {
        let delimiters = b"= \t\n\x0b\x0c\r\0".map(|byte| byte as i8);
        let nread = read_upto_delim(
            heap,
            line_ptr,
            field,
            MOL_FMT_V3000_MAXFIELDLEN as i32,
            Some(&delimiters),
        )?;
        if nread != 0 {
            let _ = mystrncpy(
                heap,
                key,
                field.as_const(),
                nread.wrapping_add(1) as u32,
            )?;
            if !line_ptr.is_null() && heap.slice(line_ptr.as_const())?[0] == b'=' as i8 {
                *line_ptr = line_ptr.offset(1)?;
            }
        } else {
            *heap
                .slice_mut(key)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
        }
        Ok(nread)
    })();
    let cleanup = heap.free(field);
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileV3000ReadCTABBeginAndCountsLine(
    heap: &mut SourceHeap,
    ctab: &mut MOL_FMT_CTAB,
    mut inp_file: Option<&mut INCHI_IOSTREAM>,
    mut p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:444 MolfileV3000ReadCTABBeginAndCountsLine
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap stream and descriptor access add overhead.
    /*
int MolfileV3000ReadCTABBeginAndCountsLine(MOL_FMT_CTAB *ctab,
                                           INCHI_IOSTREAM *inp_file,
                                           char *pStrErr)
{
    char field[MOL_FMT_V3000_MAXFIELDLEN];
    int err = 0, len; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    int failed = 0;

    int nc;
    char *p = NULL, *line = NULL;
    INCHI_IOSTREAM tmpin;
    INCHI_IOS_STRING *pin = &tmpin.s;
    inchi_ios_init(&tmpin, INCHI_IOS_TYPE_STRING, NULL);

    field[0] = '\0'; /* djb-rwth: adding zero termination */

    /* Check for proper start */

    /*p = inchi_fgetsLf_V3000( line, inp_file );*/
    inchi_strbuf_reset(pin);
    nc = get_V3000_input_line_to_strbuf(pin, inp_file);
    if (nc < 1)
    {
        p = NULL;
    }
    else
    {
        p = line = pin->pStr;
    }
    if (!p || strcmp(p, "BEGIN CTAB"))
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Error: No V3000 CTab start marker");
    }
    remove_one_lf(line);

    /* Reset all previosly read data from quasi-counts line            */
    /* (which contains only single meaningful value, 'V3000' marker    */
    ctab->n_atoms = -1;
    ctab->n_bonds = -1;
    ctab->chiral_flag = -1;
    ctab->n_stext_entries = -1;
    /* Relax stricthness of V3000 conformance: */
    /* Do not check if '999' supplied, just use this. */
    ctab->n_property_lines = 999;

    /* Read counts line */
    /*p = inchi_fgetsLf_V3000( line, inp_file );*/
    inchi_strbuf_reset(pin);
    nc = get_V3000_input_line_to_strbuf(pin, inp_file);
    if (nc < 1)
    {
        p = NULL;
    }
    else
    {
        p = line = pin->pStr;
    }
    if (!p)
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Cannot read V3000 counts line");
    }
    remove_one_lf(line);

    /* Parse counts line */
    len = MolfileV3000ReadField(field, MOL_FMT_STRING_DATA, &p); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    if (strcmp(field, "COUNTS"))
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Cannot read V3000 counts line");
    }
    failed = 0;
    if (0 > MolfileV3000ReadField(&ctab->n_atoms, MOL_FMT_INT_DATA, &p))
    {
        failed = 2;
    }
    else if (0 > MolfileV3000ReadField(&ctab->n_bonds, MOL_FMT_INT_DATA, &p))
    {
        failed = 1;
    }
    else if (0 > MolfileV3000ReadField(&ctab->v3000->n_sgroups, MOL_FMT_INT_DATA, &p))
    {
        failed = 1;
    }
    else if (0 > MolfileV3000ReadField(&ctab->v3000->n_3d_constraints, MOL_FMT_INT_DATA, &p))
    {
        failed = 1;
    }
    else if (0 > MolfileV3000ReadField(&ctab->chiral_flag, MOL_FMT_CHAR_INT_DATA, &p))
    {
        failed = 1;
    }

    if (failed)
    {
        err = 3;
        if (failed == 2)
        {
            TREAT_ERR(err, 3, "Number of atoms too large. V3000 counts line:");
        }
        else
        {
            /* too long input file line or other value min-max range mismatch */
            TREAT_ERR(err, 3, "Cannot interpret V3000 counts line:");
        }
        dotify_non_printable_chars(line);
        AddErrorMessage(pStrErr, line);
        goto err_fin;
    }

err_fin:
    inchi_strbuf_close(pin);

    return err;
}
    */
    // END INCHI C FUNCTION: MolfileV3000ReadCTABBeginAndCountsLine
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadCTABBeginAndCountsLine
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; char is signed 8-bit, short is signed 16-bit, and int is signed 32-bit.
    // INCHI✔️❌: TREAT_ERR_AND_FIN jumps to string-buffer cleanup; TREAT_ERR preserves an existing error and appends through AddErrorMessage.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadCTABBeginAndCountsLine

    let mut field = [0_i8; MOL_FMT_V3000_MAXFIELDLEN as usize];
    let mut err = 0_i32;
    let mut tmpin = INCHI_IOSTREAM::default();
    inchi_ios_init(Some(&mut tmpin), INCHI_IOS_TYPE_STRING as i32, SourceMutPointer::null())?;
    field[0] = 0;
    let result = (|| -> Result<i32, SourceHeapError> {
        inchi_strbuf_reset(heap, Some(&mut tmpin.s))?;
        let nc = get_V3000_input_line_to_strbuf(heap, &mut tmpin.s, inp_file.as_deref_mut())?;
        let p = if nc < 1 { SourceMutPointer::null() } else { tmpin.s.pStr };
        let line = p;
        if p.is_null() || !mol_fmt3_c_string_eq(heap.slice(p.as_const())?, b"BEGIN CTAB") {
            mol_fmt3_treat_error(&mut err, 1, p_str_err.as_deref_mut(), b"Error: No V3000 CTab start marker")?;
            return Ok(err);
        }
        remove_one_lf(heap, line)?;

        ctab.n_atoms = -1;
        ctab.n_bonds = -1;
        ctab.chiral_flag = -1;
        ctab.n_stext_entries = -1;
        ctab.n_property_lines = 999;

        inchi_strbuf_reset(heap, Some(&mut tmpin.s))?;
        let nc = get_V3000_input_line_to_strbuf(heap, &mut tmpin.s, inp_file.as_deref_mut())?;
        let mut p = if nc < 1 { SourceMutPointer::null() } else { tmpin.s.pStr };
        let line = p;
        if p.is_null() {
            mol_fmt3_treat_error(&mut err, 1, p_str_err.as_deref_mut(), b"Cannot read V3000 counts line")?;
            return Ok(err);
        }
        remove_one_lf(heap, line)?;

        let _len = MolfileV3000ReadField(heap, MolfileV3000FieldData::String(&mut field), i32::from(MOL_FMT_STRING_DATA), &mut p)?;
        if !mol_fmt3_c_string_eq(&field, b"COUNTS") {
            mol_fmt3_treat_error(&mut err, 1, p_str_err.as_deref_mut(), b"Cannot read V3000 counts line")?;
            return Ok(err);
        }

        let mut failed = 0_i32;
        let mut value = ctab.n_atoms;
        let read = MolfileV3000ReadField(heap, MolfileV3000FieldData::Int(&mut value), i32::from(MOL_FMT_INT_DATA), &mut p)?;
        ctab.n_atoms = value;
        if read < 0 {
            failed = 2;
        } else {
            value = ctab.n_bonds;
            let read = MolfileV3000ReadField(heap, MolfileV3000FieldData::Int(&mut value), i32::from(MOL_FMT_INT_DATA), &mut p)?;
            ctab.n_bonds = value;
            if read < 0 {
                failed = 1;
            } else {
                let current = heap.slice(ctab.v3000.as_const())?.first().ok_or(SourceHeapError::PointerOutOfBounds)?.n_sgroups;
                value = current;
                let read = MolfileV3000ReadField(heap, MolfileV3000FieldData::Int(&mut value), i32::from(MOL_FMT_INT_DATA), &mut p)?;
                heap.slice_mut(ctab.v3000)?.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)?.n_sgroups = value;
                if read < 0 {
                    failed = 1;
                } else {
                    let current = heap.slice(ctab.v3000.as_const())?.first().ok_or(SourceHeapError::PointerOutOfBounds)?.n_3d_constraints;
                    value = current;
                    let read = MolfileV3000ReadField(heap, MolfileV3000FieldData::Int(&mut value), i32::from(MOL_FMT_INT_DATA), &mut p)?;
                    heap.slice_mut(ctab.v3000)?.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)?.n_3d_constraints = value;
                    if read < 0 {
                        failed = 1;
                    } else {
                        let mut chiral = ctab.chiral_flag;
                        let read = MolfileV3000ReadField(heap, MolfileV3000FieldData::Char(&mut chiral), i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)?;
                        ctab.chiral_flag = chiral;
                        if read < 0 { failed = 1; }
                    }
                }
            }
        }

        if failed != 0 {
            err = 3;
            if failed == 2 {
                mol_fmt3_treat_error(&mut err, 3, p_str_err.as_deref_mut(), b"Number of atoms too large. V3000 counts line:")?;
            } else {
                mol_fmt3_treat_error(&mut err, 3, p_str_err.as_deref_mut(), b"Cannot interpret V3000 counts line:")?;
            }
            let _ = dotify_non_printable_chars(heap, line)?;
            let source_line = heap.slice(line.as_const())?.to_vec();
            let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(&source_line))?;
        }
        Ok(err)
    })();
    let cleanup = inchi_strbuf_close(heap, Some(&mut tmpin.s));
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileV3000ReadHapticBond(
    heap: &mut SourceHeap,
    _ctab: Option<&MOL_FMT_CTAB>,
    line_ptr: &mut SourceMutPointer<i8>,
    num_list: &mut SourceMutPointer<i32>,
    _p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1732 MolfileV3000ReadHapticBond
    // INCHI✔️❌: complete source frame follows verbatim; modeled stack field and checked heap traversal add overhead.
    /*
int MolfileV3000ReadHapticBond(MOL_FMT_CTAB *ctab,
                               char **line_ptr,
                               int **num_list,
                               char *pStrErr)
{
    int nread = 0;
    char field[MOL_FMT_V3000_MAXFIELDLEN];
    const int max_field_len = sizeof(field);
    char *p_end;
    int i, nnum = 0;

    *num_list = NULL;

    memset(field, 0, max_field_len); /* djb-rwth: memset_s C11/Annex K variant? */

    nread = read_upto_delim(line_ptr, field, max_field_len, "1234567890 \t\n\v\f\r"); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    if (strcmp(field, "("))
    {
        return -1;
    }

    nread = read_upto_delim(line_ptr, field, max_field_len, " \t\n\v\f\r"); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    nnum = strtol(field, &p_end, 10);

    if (p_end == field)
    {
        return -1; /* paranoia */
    }
    if (nnum < 0)
    {
        return -1;
    }

    *num_list = (int *)inchi_calloc((long long)nnum + 3, sizeof(int)); /* djb-rwth: cast operator added */

    if (!*num_list)
    {
        nread = -1;
        goto ret;
    }

    (*num_list)[0] = -1; /* will be bond type, to be filled by caller */
    (*num_list)[1] = -1; /* will be atom number, to be filled by caller */
    (*num_list)[2] = nnum;

    for (i = 3; i < nnum + 3; i++)
    {
        if (0 > MolfileV3000ReadField(&((*num_list)[i]), MOL_FMT_INT_DATA, line_ptr))
        {
            nread = -1;
            goto ret;
        }
    }

    /* ')' should have been consumed  by strtol */

    /* check for ATTACH=ALL */

    nread = read_upto_delim(line_ptr, field, max_field_len, " \t\n\v\f\r");
    if (nread > 0)
    {
        if (strcmp(field, "ATTACH=ALL"))
        {
            nread = -1;
            goto ret;
        }
    }

ret:
    if (nread < 0)
    {
        if (*num_list)
        {
            inchi_free(*num_list);
            *num_list = NULL;
        }
    }

    return nread;
}
    */
    // END INCHI C FUNCTION: MolfileV3000ReadHapticBond
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadHapticBond
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; inchi_calloc is calloc and int is signed 32-bit.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadHapticBond

    *num_list = SourceMutPointer::null();
    let field = heap.allocate_model_storage(vec![0_i8; MOL_FMT_V3000_MAXFIELDLEN as usize])?;
    let result = (|| -> Result<i32, SourceHeapError> {
        let opener_delimiters = b"1234567890 \t\n\x0b\x0c\r\0".map(|byte| byte as i8);
        let whitespace = b" \t\n\x0b\x0c\r\0".map(|byte| byte as i8);
        let mut nread = read_upto_delim(
            heap,
            line_ptr,
            field,
            MOL_FMT_V3000_MAXFIELDLEN as i32,
            Some(&opener_delimiters),
        )?;
        if !mol_fmt3_c_string_eq(heap.slice(field.as_const())?, b"(") {
            return Ok(-1);
        }
        nread = read_upto_delim(
            heap,
            line_ptr,
            field,
            MOL_FMT_V3000_MAXFIELDLEN as i32,
            Some(&whitespace),
        )?;
        let mut end = SourceConstPointer::null();
        let nnum = inchi_strtol(heap, field.as_const(), Some(&mut end), 10)? as i32;
        if end == field.as_const() || nnum < 0 {
            return Ok(-1);
        }
        let allocation_count = i64::from(nnum).wrapping_add(3);
        let list = match inchi_calloc::<i32>(heap, allocation_count as u64, 4) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed)
            | Err(SourceHeapError::AllocationElementCountOutOfRange)
            | Err(SourceHeapError::AllocationSizeOverflow) => return Ok(-1),
            Err(error) => return Err(error),
        };
        *num_list = list;
        {
            let values = heap.slice_mut(list)?;
            values[0] = -1;
            values[1] = -1;
            values[2] = nnum;
        }
        for i in 3..nnum.wrapping_add(3) {
            let mut value = 0_i32;
            if MolfileV3000ReadField(
                heap,
                MolfileV3000FieldData::Int(&mut value),
                i32::from(MOL_FMT_INT_DATA),
                line_ptr,
            )? < 0
            {
                nread = -1;
                break;
            }
            let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            heap.slice_mut(list)?[index] = value;
        }
        if nread >= 0 {
            nread = read_upto_delim(
                heap,
                line_ptr,
                field,
                MOL_FMT_V3000_MAXFIELDLEN as i32,
                Some(&whitespace),
            )?;
            if nread > 0 && !mol_fmt3_c_string_eq(heap.slice(field.as_const())?, b"ATTACH=ALL") {
                nread = -1;
            }
        }
        if nread < 0 {
            inchi_free(heap, list)?;
            *num_list = SourceMutPointer::null();
        }
        Ok(nread)
    })();
    let cleanup = heap.free(field);
    if result.is_err() && !num_list.is_null() {
        let _ = inchi_free(heap, *num_list);
        *num_list = SourceMutPointer::null();
    }
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileV3000ReadStereoCollection(
    heap: &mut SourceHeap,
    _ctab: Option<&MOL_FMT_CTAB>,
    line_ptr: &mut SourceMutPointer<i8>,
    num_list: &mut SourceMutPointer<i32>,
    _p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1817 MolfileV3000ReadStereoCollection
    // INCHI✔️❌: complete source frame follows verbatim; modeled stack field and checked SourceHeap access add overhead.
    /*
int MolfileV3000ReadStereoCollection(MOL_FMT_CTAB *ctab,
                                     char **line_ptr,
                                     int **num_list,
                                     char *pStrErr)
{
    int nread = 0;
    char field[MOL_FMT_V3000_MAXFIELDLEN];
    const int max_field_len = sizeof(field);
    char *p_end;
    int i, nnum = 0;

    *num_list = NULL;

    memset(field, 0, max_field_len); /* djb-rwth: memset_s C11/Annex K variant? */

    nread = read_upto_delim(line_ptr, field, max_field_len, "1234567890 \t\n\v\f\r"); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    if (strcmp(field, "("))
    {
        return -1;
    }

    nread = read_upto_delim(line_ptr, field, max_field_len, " \t\n\v\f\r");

    nnum = strtol(field, &p_end, 10);

    if (p_end == field)
    {
        return -1; /* paranoia */
    }
    if (nnum < 0)
    {
        return -1;
    }

    *num_list = (int *)inchi_calloc((long long)nnum + 3, sizeof(int)); /* djb-rwth: cast operator added */

    if (!*num_list)
    {
        nread = -1;
        goto ret;
    }

    (*num_list)[0] = -1; /* reserved, may be filled by caller */
    (*num_list)[1] = nnum;

    for (i = 2; i < nnum + 2; i++)
    {
        if (0 > MolfileV3000ReadField(&((*num_list)[i]), MOL_FMT_INT_DATA, line_ptr))
        {
            nread = -1;
            goto ret;
        }
    }

    /* ')' should have been consumed  by strtol */

ret:
    if (nread < 0)
    {
        if (*num_list)
        {
            inchi_free(*num_list);
            *num_list = NULL;
        }
    }

    return nread;
}
    */
    // END INCHI C FUNCTION: MolfileV3000ReadStereoCollection
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadStereoCollection
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; inchi_calloc is calloc and int is signed 32-bit.
    // INCHI✔️❌: ctab and pStrErr are unused; the final allocated int remains zero because allocation is count + 3 while only count + 2 entries are initialized.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadStereoCollection

    *num_list = SourceMutPointer::null();
    let field = heap.allocate_model_storage(vec![0_i8; MOL_FMT_V3000_MAXFIELDLEN as usize])?;
    let result = (|| -> Result<i32, SourceHeapError> {
        let opener_delimiters = b"1234567890 \t\n\x0b\x0c\r\0".map(|byte| byte as i8);
        let whitespace = b" \t\n\x0b\x0c\r\0".map(|byte| byte as i8);
        let _ = read_upto_delim(heap, line_ptr, field, MOL_FMT_V3000_MAXFIELDLEN as i32, Some(&opener_delimiters))?;
        if !mol_fmt3_c_string_eq(heap.slice(field.as_const())?, b"(") { return Ok(-1); }
        let mut nread = read_upto_delim(heap, line_ptr, field, MOL_FMT_V3000_MAXFIELDLEN as i32, Some(&whitespace))?;
        let mut end = SourceConstPointer::null();
        let nnum = inchi_strtol(heap, field.as_const(), Some(&mut end), 10)? as i32;
        if end == field.as_const() || nnum < 0 { return Ok(-1); }
        let allocation_count = i64::from(nnum).wrapping_add(3);
        let list = match inchi_calloc::<i32>(heap, allocation_count as u64, 4) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed)
            | Err(SourceHeapError::AllocationElementCountOutOfRange)
            | Err(SourceHeapError::AllocationSizeOverflow) => return Ok(-1),
            Err(error) => return Err(error),
        };
        *num_list = list;
        heap.slice_mut(list)?[0] = -1;
        heap.slice_mut(list)?[1] = nnum;
        let mut i = 2_i32;
        while i < nnum.wrapping_add(2) {
            let mut value = 0_i32;
            if MolfileV3000ReadField(heap, MolfileV3000FieldData::Int(&mut value), i32::from(MOL_FMT_INT_DATA), line_ptr)? < 0 {
                nread = -1;
                break;
            }
            let position = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            heap.slice_mut(list)?[position] = value;
            i = i.wrapping_add(1);
        }
        if nread < 0 {
            inchi_free(heap, list)?;
            *num_list = SourceMutPointer::null();
        }
        Ok(nread)
    })();
    let cleanup = heap.free(field);
    if result.is_err() && !num_list.is_null() {
        let _ = inchi_free(heap, *num_list);
        *num_list = SourceMutPointer::null();
    }
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileV3000ReadAtomsBlock(
    heap: &mut SourceHeap,
    ctab: &mut MOL_FMT_CTAB,
    mut inp_file: Option<&mut INCHI_IOSTREAM>,
    mut err: i32,
    mut p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:862 MolfileV3000ReadAtomsBlock
    // INCHI✔️❌: complete source frame follows verbatim; modeled stack storage, checked SourceHeap access, and owned formatter output add overhead.
    /*
int MolfileV3000ReadAtomsBlock(MOL_FMT_CTAB *ctab,
                               INCHI_IOSTREAM *inp_file,
                               int err,
                               char *pStrErr)
{
    int i;
    /* djb-rwth: removing redundant variables */
    char field[MOL_FMT_V3000_MAXFIELDLEN];
    int nc, failed = 0;
    char *p = NULL, *line = NULL;
    INCHI_IOSTREAM tmpin;
    INCHI_IOS_STRING *pin = &tmpin.s;
    inchi_ios_init(&tmpin, INCHI_IOS_TYPE_STRING, NULL);

    /* Check for proper start */
    /*p = inchi_fgetsLf_V3000( line, inp_file );*/

    nc = get_V3000_input_line_to_strbuf(pin, inp_file);

    if (nc < 1)
    {
        p = NULL;
    }
    else
    {
        p = line = pin->pStr;
    }
    if (!p || strcmp(p, "BEGIN ATOM"))
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Error: No V3000 Atom block start marker");
    }
    remove_one_lf(line);

    ctab->v3000->n_non_star_atoms = 0;
    for (i = 0; i < ctab->n_atoms; i++)
    {
        int ii = -1;

        /*p = inchi_fgetsLf_V3000( line, inp_file );*/
        inchi_strbuf_reset(pin);

        nc = get_V3000_input_line_to_strbuf(pin, inp_file);

        if (nc < 1)
        {
            p = NULL;
        }
        else
        {
            p = line = pin->pStr;
        }
        if (!p)
        {
            if (!err)
            {
                TREAT_ERR(err, 2, "Cannot read V3000 atom block line");
            }
            break;
        }
        remove_one_lf(line);

        if (err)
        {
            if (!strcmp(line, SD_FMT_END_OF_DATA))
            {
                err = -abs(err);
                break;
            }
            continue; /* bypass the rest of the Atom block */
        }

        if (ctab->atoms)
        {
            int index, aamap; /* not used actually, just read them */
            int len;
            char symbol[6]; /* TODO: treat possibly long V3000 atom names */
            double fx = 0.0, fy = 0.0, fz = 0.0;
#ifdef GHI100_FIX
#if (SPRINTF_FLAG == 2)
            char *fxs, *fys, *fzs;
            int rfxs, rfys, rfzs;

            fxs = (char *)inchi_malloc((10 + 3) * sizeof(double));
            fys = (char *)inchi_malloc((10 + 3) * sizeof(double));
            fzs = (char *)inchi_malloc((10 + 3) * sizeof(double));

            if (fxs || fys || fzs)
            {
                failed = 1;
            }
#endif
#endif
            symbol[0] = '\0'; /* djb-rwth: adding zero termination */

            /* Read positional parameters */
            failed = 0;

            if (0 > MolfileV3000ReadField(&index, MOL_FMT_INT_DATA, &p))
            {
                failed = 1;
            }
            else if (0 > MolfileV3000ReadField(&symbol, MOL_FMT_STRING_DATA, &p))
            {
                failed = 1;
            }
            else if (0 > MolfileV3000ReadField(&fx, MOL_FMT_DOUBLE_DATA, &p))
            {
                failed = 1;
            }
            else if (0 > MolfileV3000ReadField(&fy, MOL_FMT_DOUBLE_DATA, &p))
            {
                failed = 1;
            }
            else if (0 > MolfileV3000ReadField(&fz, MOL_FMT_DOUBLE_DATA, &p))
            {
                failed = 1;
            }
            else if (0 > MolfileV3000ReadField(&aamap, MOL_FMT_INT_DATA, &p))
            {
                failed = 1;
            }

            if (failed)
            {

                err = 4;
                TREAT_ERR(err, 4, "Cannot interpret V3000 atom block line:"); /* djb-rwth: addressing coverity ID #499547 -- TREAT_ERR properly used */
                dotify_non_printable_chars(line);
                AddErrorMessage(pStrErr, line);

                if (!strcmp(line, SD_FMT_END_OF_DATA))
                {
                    err = -abs(err);
                    break;
                }
                continue; /* can't interpret a first half of atom block line */
            }

            /* emulate V2000 coordinates substring */
            if (ctab->coords)
            {
                char szcoords[40];
#ifdef GHI100_FIX
#if (SPRINTF_FLAG == 2)
                if (fxs)
                {
                    rfxs = dbl2int(fxs, 10, -1, 'g', fx);
                }
                if (fys)
                {
                    rfys = dbl2int(fys, 10, -1, 'g', fy);
                }
                if (fzs)
                {
                    rfzs = dbl2int(fzs, 10, -1, 'g', fz);
                }
                if ((rfxs >= 0) && (rfys >= 0) && (rfzs >= 0))
                {
                    sprintf(szcoords, "%s%s%s", fxs, fys, fzs);
                }
                else
                {
                    failed = 1;
                }
                inchi_free(fxs);
                inchi_free(fys);
                inchi_free(fzs);
#elif (SPRINTF_FLAG == 2)
                stbsp_sprintf(szcoords, "%10g%10g%10g", fx, fy, fz);
#else
                sprintf(szcoords, "%10g%10g%10g", fx, fy, fz);
#endif
#endif
                sprintf(szcoords, "%10g%10g%10g", fx, fy, fz);
                strcpy(ctab->coords[i], szcoords);
            }

            if (!strcmp(symbol, "*"))
            {
                /* ignore star atoms but save index info */
                ctab->v3000->atom_index_orig[i] = index;
                ctab->v3000->atom_index_fin[i] = -1;
                ctab->v3000->n_star_atoms++;
                continue;
            }

            ctab->v3000->n_non_star_atoms++;
            ctab->v3000->atom_index_orig[i] = index;
            ctab->v3000->atom_index_fin[i] = ctab->v3000->n_non_star_atoms;
            ii = ctab->v3000->n_non_star_atoms - 1;

            mystrncpy(ctab->atoms[ii].symbol, symbol, sizeof(ctab->atoms[ii].symbol));
            if (2 == strlen(ctab->atoms[ii].symbol) && isupper(UCINT ctab->atoms[ii].symbol[1]))
            {
                ctab->atoms[ii].symbol[1] = (char)tolower(UCINT ctab->atoms[ii].symbol[1]); /* 5-4-99 DCh*/
            }
            ctab->atoms[ii].fx = fx;
            ctab->atoms[ii].fy = fy;
            ctab->atoms[ii].fz = fz;

            /* Read key-val pairs if any */
            while (p && (len = MolfileV3000ReadKeyword(field, &p)) > 0) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
            {

                int itmp;
                char ctmp;
                char stmp[MOL_FMT_V3000_MAXFIELDLEN];

                failed = 0;
                if (!strcmp(field, "CHG"))
                {
                    if (0 > MolfileV3000ReadField(&ctab->atoms[ii].charge, MOL_FMT_CHAR_INT_DATA, &p))
                    {
                        failed = 1;
                    }
                }
                else if (!strcmp(field, "RAD"))
                {
                    if (0 > MolfileV3000ReadField(&ctab->atoms[ii].radical, MOL_FMT_CHAR_INT_DATA, &p))
                    {
                        failed = 1;
                    }
                }
                else if (!strcmp(field, "CFG"))
                {
                    if (0 > MolfileV3000ReadField(&ctab->atoms[ii].stereo_parity, MOL_FMT_CHAR_INT_DATA, &p))
                    {
                        failed = 1;
                    }
                }

                else if (!strcmp(field, "MASS"))
                {
                    /*
                        Default = natural abundance
                        A specified value indicates the absolute
                        atomic weight of the designated atom.
                    */
                    S_SHORT iso_mass;
                    if (0 > MolfileV3000ReadField(&iso_mass, MOL_FMT_SHORT_INT_DATA, &p))
                    {
                        failed = 1;
                        TREAT_ERR(err, 0, "Isotopic data not recognized:");
                        AddErrorMessage(pStrErr, line);
                        /* ignore isotopic error for now */
                    }
                    else
                    {
                        /*  What we read is an absolute isotopic mass, by V3000 spec.
                            Adjust this to old convention for further processing:
                            set 'ctab->atoms[ii].mass_difference' to 127
                            if isotopic mass is the same as element mass
                            in Periodic Table (rounded avg by all isotopes), 'atw'
                            delta otherwise, the value of difference 'delta' = ( isotopic mass - 'atw')
                        */
                        int atw, delta;
                        atw = get_atomic_mass(ctab->atoms[ii].symbol);
                        delta = (int)iso_mass - atw;
                        ctab->atoms[ii].mass_difference = (char)(delta ? delta : ZERO_ATW_DIFF);
                    }
                }

                else if (!strcmp(field, "VAL"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        failed = 1;
                    }
                    else
                    {
                        /* adjust to old convention: was 15 for zero, now -1 for zero */
                        if (itmp == -1)
                        {
                            ctmp = 15;
                        }
                        else
                        {
                            ctmp = (char)itmp;
                        }
                        ctab->atoms[ii].valence = ctmp;
                    }
                }
                else if (!strcmp(field, "HCOUNT"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        ; /* skip query-related stuff */
                    }
                }
                else if (!strcmp(field, "STBOX"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        ; /* skip for now */
                    }
                }
                else if (!strcmp(field, "INVRET") || !strcmp(field, "EXACHG"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        ; /* skip reaction-related stuff */
                    }
                }
                else if (!strcmp(field, "SUBST") || !strcmp(field, "UNSAT") || !strcmp(field, "RBCNT"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        ; /* skip query-related stuff */
                    }
                }
                else if (!strcmp(field, "ATTCHPT"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        ;
                    }
                }
                else if (!strcmp(field, "RGROUPS"))
                {
                    if (0 > MolfileV3000ReadField(&stmp, MOL_FMT_STRING_DATA, &p))
                    {
                        ;
                    }
                }
                else if (!strcmp(field, "ATTCHORD"))
                {
                    if (0 > MolfileV3000ReadField(&stmp, MOL_FMT_STRING_DATA, &p))
                    {
                        ;
                    }
                }
                else if (!strcmp(field, "CLASS"))
                {
                    if (0 > MolfileV3000ReadField(&stmp, MOL_FMT_STRING_DATA, &p))
                    {
                        ;
                    }
                }
                else if (!strcmp(field, "SEQID"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        ;
                    }
                }

                if (failed)
                {
                    err = 4;
                    TREAT_ERR(err, 4, "Cannot interpret V3000 atom block key-value pair");
                    dotify_non_printable_chars(line);
                    AddErrorMessage(pStrErr, line);

                    if (!strcmp(line, SD_FMT_END_OF_DATA))
                    {
                        err = -abs(err);
                        break;
                    }
                    continue;
                }
            }
        } /* if ( NULL != ctab->atoms )  */
    } /* for ( i = 0; i < ctab->n_atoms; i++ )  */

    if (ctab->v3000->n_star_atoms)
    {
        AddErrorMessage(pStrErr, "V3000 star atoms ignored");
        ctab->n_atoms = ctab->v3000->n_non_star_atoms;
    }

    /* Check for proper finish */

    /*p = inchi_fgetsLf_V3000( line, inp_file );*/
    inchi_strbuf_reset(pin);

    nc = get_V3000_input_line_to_strbuf(pin, inp_file);

    if (nc < 1)
    {
        p = NULL;
    }
    else
    {
        p = line = pin->pStr;
    }
    if (!p || strcmp(p, "END ATOM"))
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Error: No V3000 Atom block end marker");
    }
    remove_one_lf(line);

err_fin:
    inchi_strbuf_close(pin);

    return err;
}
    */
    // END INCHI C FUNCTION: MolfileV3000ReadAtomsBlock
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadAtomsBlock
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; GHI100_FIX is undefined, so the sole active coordinate formatter is sprintf("%10g%10g%10g", fx, fy, fz).
    // INCHI✔️❌: SD_FMT_END_OF_DATA is "$$$$"; plain char is signed; short is signed 16-bit; int is signed 32-bit; ZERO_ATW_DIFF is 127.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadAtomsBlock

    let mut tmpin = INCHI_IOSTREAM::default();
    inchi_ios_init(Some(&mut tmpin), INCHI_IOS_TYPE_STRING as i32, SourceMutPointer::null())?;
    let field = heap.allocate_model_storage(vec![0_i8; MOL_FMT_V3000_MAXFIELDLEN as usize])?;
    let coordinate_format = heap.allocate_model_storage(
        b"%10g%10g%10g\0".iter().map(|byte| *byte as i8).collect::<Vec<_>>(),
    )?;
    let result = (|| -> Result<i32, SourceHeapError> {
        let nc = get_V3000_input_line_to_strbuf(heap, &mut tmpin.s, inp_file.as_deref_mut())?;
        let mut p = if nc < 1 { SourceMutPointer::null() } else { tmpin.s.pStr };
        let mut line = p;
        if p.is_null() || !mol_fmt3_c_string_eq(heap.slice(p.as_const())?, b"BEGIN ATOM") {
            mol_fmt3_treat_error(&mut err, 1, p_str_err.as_deref_mut(), b"Error: No V3000 Atom block start marker")?;
            return Ok(err);
        }
        remove_one_lf(heap, line)?;
        heap.slice_mut(ctab.v3000)?[0].n_non_star_atoms = 0;

        let mut i = 0_i32;
        while i < ctab.n_atoms {
            inchi_strbuf_reset(heap, Some(&mut tmpin.s))?;
            let nc = get_V3000_input_line_to_strbuf(heap, &mut tmpin.s, inp_file.as_deref_mut())?;
            p = if nc < 1 { SourceMutPointer::null() } else { tmpin.s.pStr };
            line = p;
            if p.is_null() {
                if err == 0 {
                    mol_fmt3_treat_error(&mut err, 2, p_str_err.as_deref_mut(), b"Cannot read V3000 atom block line")?;
                }
                break;
            }
            remove_one_lf(heap, line)?;
            if err != 0 {
                if mol_fmt3_c_string_eq(heap.slice(line.as_const())?, &SD_FMT_END_OF_DATA[..SD_FMT_END_OF_DATA.len() - 1]) {
                    err = err.wrapping_abs().wrapping_neg();
                    break;
                }
                i = i.wrapping_add(1);
                continue;
            }

            if !ctab.atoms.is_null() {
                let mut index = 0_i32;
                let mut aamap = 0_i32;
                let mut symbol = [0_i8; 6];
                let mut fx = 0.0_f64;
                let mut fy = 0.0_f64;
                let mut fz = 0.0_f64;
                let mut failed = 0_i32;
                if MolfileV3000ReadField(heap, MolfileV3000FieldData::Int(&mut index), i32::from(MOL_FMT_INT_DATA), &mut p)? < 0 {
                    failed = 1;
                } else if MolfileV3000ReadField(heap, MolfileV3000FieldData::String(&mut symbol), i32::from(MOL_FMT_STRING_DATA), &mut p)? < 0 {
                    failed = 1;
                } else if MolfileV3000ReadField(heap, MolfileV3000FieldData::Double(&mut fx), i32::from(MOL_FMT_DOUBLE_DATA), &mut p)? < 0 {
                    failed = 1;
                } else if MolfileV3000ReadField(heap, MolfileV3000FieldData::Double(&mut fy), i32::from(MOL_FMT_DOUBLE_DATA), &mut p)? < 0 {
                    failed = 1;
                } else if MolfileV3000ReadField(heap, MolfileV3000FieldData::Double(&mut fz), i32::from(MOL_FMT_DOUBLE_DATA), &mut p)? < 0 {
                    failed = 1;
                } else if MolfileV3000ReadField(heap, MolfileV3000FieldData::Int(&mut aamap), i32::from(MOL_FMT_INT_DATA), &mut p)? < 0 {
                    failed = 1;
                }
                if failed != 0 {
                    err = 4;
                    mol_fmt3_treat_error(&mut err, 4, p_str_err.as_deref_mut(), b"Cannot interpret V3000 atom block line:")?;
                    let _ = dotify_non_printable_chars(heap, line)?;
                    let source_line = heap.slice(line.as_const())?.to_vec();
                    let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(&source_line))?;
                    if mol_fmt3_c_string_eq(heap.slice(line.as_const())?, &SD_FMT_END_OF_DATA[..SD_FMT_END_OF_DATA.len() - 1]) {
                        err = err.wrapping_abs().wrapping_neg();
                        break;
                    }
                    i = i.wrapping_add(1);
                    continue;
                }

                if !ctab.coords.is_null() {
                    let mut arguments = SourceVaList {
                        arguments: vec![
                            SourceFormatArgument::Float(fx),
                            SourceFormatArgument::Float(fy),
                            SourceFormatArgument::Float(fz),
                        ],
                        position: 0,
                    };
                    let formatted = source_vformat(heap, coordinate_format.as_const(), &mut arguments)?;
                    let destination = heap.slice_mut(ctab.coords)?
                        .get_mut(usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let required = formatted.len().checked_add(1).ok_or(SourceHeapError::AllocationSizeOverflow)?;
                    if destination.len() < required {
                        return Err(SourceHeapError::PointerOutOfBounds);
                    }
                    for (target, source) in destination.iter_mut().zip(&formatted) {
                        *target = *source as i8;
                    }
                    destination[formatted.len()] = 0;
                }

                let atom_position = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if mol_fmt3_c_string_eq(&symbol, b"*") {
                    let v3000 = heap.slice(ctab.v3000.as_const())?.first().cloned().ok_or(SourceHeapError::PointerOutOfBounds)?;
                    heap.slice_mut(v3000.atom_index_orig)?[atom_position] = index;
                    heap.slice_mut(v3000.atom_index_fin)?[atom_position] = -1;
                    heap.slice_mut(ctab.v3000)?[0].n_star_atoms = v3000.n_star_atoms.wrapping_add(1);
                    i = i.wrapping_add(1);
                    continue;
                }

                let ii = {
                    let v3000 = heap.slice_mut(ctab.v3000)?.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)?;
                    v3000.n_non_star_atoms = v3000.n_non_star_atoms.wrapping_add(1);
                    v3000.n_non_star_atoms
                };
                let v3000 = heap.slice(ctab.v3000.as_const())?.first().cloned().ok_or(SourceHeapError::PointerOutOfBounds)?;
                heap.slice_mut(v3000.atom_index_orig)?[atom_position] = index;
                heap.slice_mut(v3000.atom_index_fin)?[atom_position] = ii;
                let ii = usize::try_from(ii.wrapping_sub(1)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                {
                    let atom = heap.slice_mut(ctab.atoms)?.get_mut(ii).ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let _ = mystrncpy_slice(Some(&mut atom.symbol), Some(&symbol), 6)?;
                    let symbol_length = atom.symbol.iter().position(|byte| *byte == 0).unwrap_or(atom.symbol.len());
                    if symbol_length == 2 && (atom.symbol[1] as u8).is_ascii_uppercase() {
                        atom.symbol[1] = (atom.symbol[1] as u8).to_ascii_lowercase() as i8;
                    }
                    atom.fx = fx;
                    atom.fy = fy;
                    atom.fz = fz;
                }

                while !p.is_null() && MolfileV3000ReadKeyword(heap, field, &mut p)? > 0 {
                    let key = heap.slice(field.as_const())?.to_vec();
                    let mut keyword_failed = 0_i32;
                    if mol_fmt3_c_string_eq(&key, b"CHG") {
                        let mut value = heap.slice(ctab.atoms.as_const())?[ii].charge;
                        let read = MolfileV3000ReadField(heap, MolfileV3000FieldData::Char(&mut value), i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)?;
                        heap.slice_mut(ctab.atoms)?[ii].charge = value;
                        if read < 0 { keyword_failed = 1; }
                    } else if mol_fmt3_c_string_eq(&key, b"RAD") {
                        let mut value = heap.slice(ctab.atoms.as_const())?[ii].radical;
                        let read = MolfileV3000ReadField(heap, MolfileV3000FieldData::Char(&mut value), i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)?;
                        heap.slice_mut(ctab.atoms)?[ii].radical = value;
                        if read < 0 { keyword_failed = 1; }
                    } else if mol_fmt3_c_string_eq(&key, b"CFG") {
                        let mut value = heap.slice(ctab.atoms.as_const())?[ii].stereo_parity;
                        let read = MolfileV3000ReadField(heap, MolfileV3000FieldData::Char(&mut value), i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)?;
                        heap.slice_mut(ctab.atoms)?[ii].stereo_parity = value;
                        if read < 0 { keyword_failed = 1; }
                    } else if mol_fmt3_c_string_eq(&key, b"MASS") {
                        let mut iso_mass = 0_i16;
                        if MolfileV3000ReadField(heap, MolfileV3000FieldData::Short(&mut iso_mass), i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0 {
                            keyword_failed = 1;
                            mol_fmt3_treat_error(&mut err, 0, p_str_err.as_deref_mut(), b"Isotopic data not recognized:")?;
                            let source_line = heap.slice(line.as_const())?.to_vec();
                            let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(&source_line))?;
                        } else {
                            let atom_symbol = heap.slice(ctab.atoms.as_const())?[ii].symbol;
                            let atw = get_atomic_mass(Some(&atom_symbol))?;
                            let delta = i32::from(iso_mass).wrapping_sub(atw);
                            heap.slice_mut(ctab.atoms)?[ii].mass_difference = if delta != 0 { delta as i8 } else { 127_i8 };
                        }
                    } else if mol_fmt3_c_string_eq(&key, b"VAL") {
                        let mut value = 0_i32;
                        if MolfileV3000ReadField(heap, MolfileV3000FieldData::Int(&mut value), i32::from(MOL_FMT_INT_DATA), &mut p)? < 0 {
                            keyword_failed = 1;
                        } else {
                            heap.slice_mut(ctab.atoms)?[ii].valence = if value == -1 { 15 } else { value as i8 };
                        }
                    } else if mol_fmt3_c_string_eq(&key, b"HCOUNT")
                        || mol_fmt3_c_string_eq(&key, b"STBOX")
                        || mol_fmt3_c_string_eq(&key, b"INVRET")
                        || mol_fmt3_c_string_eq(&key, b"EXACHG")
                        || mol_fmt3_c_string_eq(&key, b"SUBST")
                        || mol_fmt3_c_string_eq(&key, b"UNSAT")
                        || mol_fmt3_c_string_eq(&key, b"RBCNT")
                        || mol_fmt3_c_string_eq(&key, b"ATTCHPT")
                        || mol_fmt3_c_string_eq(&key, b"SEQID")
                    {
                        let mut ignored = 0_i32;
                        let _ = MolfileV3000ReadField(heap, MolfileV3000FieldData::Int(&mut ignored), i32::from(MOL_FMT_INT_DATA), &mut p)?;
                    } else if mol_fmt3_c_string_eq(&key, b"RGROUPS")
                        || mol_fmt3_c_string_eq(&key, b"ATTCHORD")
                        || mol_fmt3_c_string_eq(&key, b"CLASS")
                    {
                        let mut ignored = [0_i8; MOL_FMT_V3000_MAXFIELDLEN as usize];
                        let _ = MolfileV3000ReadField(heap, MolfileV3000FieldData::String(&mut ignored), i32::from(MOL_FMT_STRING_DATA), &mut p)?;
                    }

                    if keyword_failed != 0 {
                        err = 4;
                        mol_fmt3_treat_error(&mut err, 4, p_str_err.as_deref_mut(), b"Cannot interpret V3000 atom block key-value pair")?;
                        let _ = dotify_non_printable_chars(heap, line)?;
                        let source_line = heap.slice(line.as_const())?.to_vec();
                        let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(&source_line))?;
                        if mol_fmt3_c_string_eq(heap.slice(line.as_const())?, &SD_FMT_END_OF_DATA[..SD_FMT_END_OF_DATA.len() - 1]) {
                            err = err.wrapping_abs().wrapping_neg();
                            break;
                        }
                        continue;
                    }
                }
                let _ = aamap;
            }
            i = i.wrapping_add(1);
        }

        let v3000 = heap.slice(ctab.v3000.as_const())?.first().cloned().ok_or(SourceHeapError::PointerOutOfBounds)?;
        if v3000.n_star_atoms != 0 {
            let _ = mol_fmt3_add_message(p_str_err.as_deref_mut(), b"V3000 star atoms ignored")?;
            ctab.n_atoms = v3000.n_non_star_atoms;
        }
        inchi_strbuf_reset(heap, Some(&mut tmpin.s))?;
        let nc = get_V3000_input_line_to_strbuf(heap, &mut tmpin.s, inp_file.as_deref_mut())?;
        p = if nc < 1 { SourceMutPointer::null() } else { tmpin.s.pStr };
        line = p;
        if p.is_null() || !mol_fmt3_c_string_eq(heap.slice(p.as_const())?, b"END ATOM") {
            mol_fmt3_treat_error(&mut err, 1, p_str_err.as_deref_mut(), b"Error: No V3000 Atom block end marker")?;
            return Ok(err);
        }
        remove_one_lf(heap, line)?;
        Ok(err)
    })();
    let field_cleanup = heap.free(field);
    let format_cleanup = heap.free(coordinate_format);
    let stream_cleanup = inchi_strbuf_close(heap, Some(&mut tmpin.s));
    match (result, field_cleanup, format_cleanup, stream_cleanup) {
        (Err(error), _, _, _)
        | (Ok(_), Err(error), _, _)
        | (Ok(_), Ok(()), Err(error), _)
        | (Ok(_), Ok(()), Ok(()), Err(error)) => Err(error),
        (Ok(value), Ok(()), Ok(()), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileV3000ReadBondsBlock(
    heap: &mut SourceHeap,
    ctab: &mut MOL_FMT_CTAB,
    mut inp_file: Option<&mut INCHI_IOSTREAM>,
    mut err: i32,
    mut p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1262 MolfileV3000ReadBondsBlock
    // INCHI✔️❌: complete source frame follows verbatim; modeled stack storage and checked SourceHeap access add overhead.
    /*
int MolfileV3000ReadBondsBlock(MOL_FMT_CTAB *ctab,
                               INCHI_IOSTREAM *inp_file,
                               int err,
                               char *pStrErr)
{
    int i;
    char field[MOL_FMT_V3000_MAXFIELDLEN];
    int nc;
    char *p = NULL, *line = NULL;
    INCHI_IOSTREAM tmpin;
    INCHI_IOS_STRING *pin = &tmpin.s;

    if (!ctab->n_bonds)
    {
        return 0;
    }
    inchi_ios_init(&tmpin, INCHI_IOS_TYPE_STRING, NULL);

    /* Check for proper start */
    /*p = inchi_fgetsLf_V3000( line, inp_file );*/

    nc = get_V3000_input_line_to_strbuf(pin, inp_file);

    if (nc < 1)
    {
        p = NULL;
    }
    else
    {
        p = line = pin->pStr;
    }
    if (!p || strcmp(p, "BEGIN BOND"))
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Error: No V3000 Bond block start marker");
    }
    remove_one_lf(line);

    ctab->v3000->n_haptic_bonds = 0;
    ctab->v3000->n_non_haptic_bonds = 0;

    for (i = 0; i < ctab->n_bonds; i++)
    {
        int is_haptic = 0;

        /*p = inchi_fgetsLf_V3000( line, inp_file );*/
        inchi_strbuf_reset(pin);

        nc = get_V3000_input_line_to_strbuf(pin, inp_file);

        if (nc < 1)
        {
            p = NULL;
        }
        else
        {
            p = line = pin->pStr;
        }
        if (!p)
        {
            if (!err)
            {
                TREAT_ERR(err, 2, "Cannot read V3000 bond block line"); /* djb-rwth: addressing coverity ID #499565 -- TREAT_ERR properly used */
            }
            break;
        }
        remove_one_lf(line);

        if (err)
        {
            if (!strcmp(line, SD_FMT_END_OF_DATA))
            {
                err = -abs(err);
                break;
            }
            continue;
        }

        if (ctab->bonds)
        {
            int index, n_orig_at, len;
            short int atnum1 = -1, atnum2 = -1;
            char bond_type = 0, stereo = 0;
            int failed = 0;
            /* djb-rwth: removing redundant variables */

            n_orig_at = ctab->v3000->n_non_star_atoms + ctab->v3000->n_star_atoms;

            /* read positional parameters */
            if (0 > MolfileV3000ReadField(&index, MOL_FMT_INT_DATA, &p))
            {
                failed = 1;
            }
            else if (0 > MolfileV3000ReadField(&bond_type, MOL_FMT_CHAR_INT_DATA, &p))
            {
                failed = 1;
            }
            else if (0 > MolfileV3000ReadField(&atnum1, MOL_FMT_SHORT_INT_DATA, &p))
            {
                failed = 1;
            }
            else if (0 > MolfileV3000ReadField(&atnum2, MOL_FMT_SHORT_INT_DATA, &p))
            {
                failed = 1;
            }

            atnum1 = get_actual_atom_number(atnum1, n_orig_at,
                                            ctab->v3000->atom_index_orig,
                                            ctab->v3000->atom_index_fin);

            atnum2 = get_actual_atom_number(atnum2, n_orig_at,
                                            ctab->v3000->atom_index_orig,
                                            ctab->v3000->atom_index_fin);

            if ((atnum1 < 0) && (atnum2 < 0))
            {
                failed = 1;
            }
            /* djb-rwth: removing redundant code */

            if (failed)
            {

                if (!err)
                {
                    /* can't interpret bonds block line */
                    TREAT_ERR(err, 4, "Cannot interpret V3000 bond block line:");
                    dotify_non_printable_chars(line);
                    AddErrorMessage(pStrErr, line);
                }
                if (!strcmp(line, SD_FMT_END_OF_DATA))
                {
                    err = -abs(err);
                    break;
                }
            }

            /* TODO: treat new bond types  9 10 */
            /* read key-val pairs if any */
            while (p && (len = MolfileV3000ReadKeyword(field, &p)) > 0) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
            {

                int itmp;
                char stmp[MOL_FMT_V3000_MAXFIELDLEN];
                failed = 0;

                if (!strcmp(field, "CFG"))
                {
                    if (0 > MolfileV3000ReadField(&stereo, MOL_FMT_CHAR_INT_DATA, &p))
                    {
                        failed = 1;
                    }
                    else
                    {
                        /*    adjust stereo to old convention for wedges which was:
                                0 = not stereo, 1 = Up,  4 = Either, 6 = Down
                            now:
                                0 = none (default), 1 = up, 2 = either, 3 = down
                        */
                        if (stereo == 2)
                        {
                            stereo = 4;
                        }
                        else if (stereo == 3)
                        {
                            stereo = 6;
                        }
                    }
                }
                else if (!strcmp(field, "TOPO"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        ; /* skip query-related stuff */
                    }
                }
                else if (!strcmp(field, "RXCTR"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        ; /* skip reaction-related stuff */
                    }
                }
                else if (!strcmp(field, "STBOX"))
                {
                    if (0 > MolfileV3000ReadField(&itmp, MOL_FMT_INT_DATA, &p))
                    {
                        ; /* skip for now */
                    }
                }
                else if (!strcmp(field, "ENDPTS"))
                {
                    int res, *num_list = NULL;
                    if (0 > MolfileV3000ReadHapticBond(ctab, &p, &num_list, pStrErr))
                    {
                        failed = 1;
                    }
                    else if (!num_list)
                    {
                        failed = 1;
                    }
                    else
                    {
                        int existent_atom = atnum1;
                        if (existent_atom < 0)
                        {
                            existent_atom = atnum2;
                        }
                        if (existent_atom < 0) /* should not be here */
                        {
                            failed = 1;
                        }
                        else
                        {
                            int k, nnum;
                            nnum = num_list[2];
                            num_list[1] = existent_atom;
                            for (k = 3; k < nnum; k++)
                            {
                                num_list[k] = get_actual_atom_number(num_list[k],
                                                                     n_orig_at,
                                                                     ctab->v3000->atom_index_orig,
                                                                     ctab->v3000->atom_index_fin);
                            }
                            res = NumLists_Append(ctab->v3000->haptic_bonds, num_list);
                            if (res < 0)
                            {
                                failed = 1;
                            }
                            else
                            {
                                is_haptic = 1;
                            }
                        }
                    }
                    /* djb-rwth: addressing coverity ID #499489 -- false positive as num_atoms allocated in MolfileV3000ReadHapticBond and returns a value in this block */
                }
                else if (!strcmp(field, "DISP"))
                {
                    if (0 > MolfileV3000ReadField(&stmp, MOL_FMT_STRING_DATA, &p))
                    {
                        ;
                    }
                }
                else if (!strcmp(field, "ATTACH"))
                {
                    if (0 > MolfileV3000ReadField(&stmp, MOL_FMT_STRING_DATA, &p))
                    {
                        ;
                    }
                }

                if (failed)
                {
                    if (!err)
                    {
                        /* can't interpret bonds block line */
                        TREAT_ERR(err, 4, "Cannot interpret V3000 bond block line:");
                        dotify_non_printable_chars(line);
                        AddErrorMessage(pStrErr, line);
                    }
                    if (!strcmp(line, SD_FMT_END_OF_DATA))
                    {
                        err = -abs(err);
                        break;
                    }
                }
            } /* while ( p && (len=MolfileV3000ReadKeyword(field, &p)) > 0 ) */

            if (is_haptic)
            {
                int ii = ctab->v3000->n_haptic_bonds;
                ctab->v3000->haptic_bonds->lists[ii][0] = bond_type;
                ctab->v3000->n_haptic_bonds++;
                continue;
            }
            else
            {
                int ii = ctab->v3000->n_non_haptic_bonds;
                ctab->bonds[ii].atnum1 = atnum1;
                ctab->bonds[ii].atnum2 = atnum2;
                ctab->bonds[ii].bond_type = bond_type;
                ctab->bonds[ii].bond_stereo = stereo;
                ctab->v3000->n_non_haptic_bonds++;
            }
        } /* if ctab->bonds */
    } /* for ( i = 0; i < ctab->n_bonds; i++ )  */

    if (ctab->v3000->n_haptic_bonds)
    {
        AddErrorMessage(pStrErr, "V3000 haptic bonds read/stored but ignored");
        ctab->n_bonds = ctab->v3000->n_non_haptic_bonds;
    }

    /* Check for proper finish */
    /*p = inchi_fgetsLf_V3000( line, inp_file );*/
    inchi_strbuf_reset(pin);

    nc = get_V3000_input_line_to_strbuf(pin, inp_file);

    if (nc < 1)
    {
        p = NULL;
    }
    else
    {
        p = line = pin->pStr;
    }
    if (!p || strcmp(p, "END BOND"))
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Error: No V3000 Bond block end marker");
    }
    remove_one_lf(line);

err_fin:

    inchi_ios_close(&tmpin); /* ricrogz: fixing memory leak */
    return err;
}
    */
    // END INCHI C FUNCTION: MolfileV3000ReadBondsBlock
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadBondsBlock
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; TREAT_ERR preserves a pre-existing error and always invokes AddErrorMessage.
    // INCHI✔️❌: SD_FMT_END_OF_DATA is "$$$$"; char is signed; short is signed 16-bit; int is signed 32-bit.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadBondsBlock

    if ctab.n_bonds == 0 { return Ok(0); }
    let mut tmpin = INCHI_IOSTREAM::default();
    inchi_ios_init(Some(&mut tmpin), INCHI_IOS_TYPE_STRING as i32, SourceMutPointer::null())?;
    let field = heap.allocate_model_storage(vec![0_i8; MOL_FMT_V3000_MAXFIELDLEN as usize])?;
    let result = (|| -> Result<i32, SourceHeapError> {
        let nc = get_V3000_input_line_to_strbuf(heap, &mut tmpin.s, inp_file.as_deref_mut())?;
        let mut p = if nc < 1 { SourceMutPointer::null() } else { tmpin.s.pStr };
        let mut line = p;
        if p.is_null() || !mol_fmt3_c_string_eq(heap.slice(p.as_const())?, b"BEGIN BOND") {
            mol_fmt3_treat_error(&mut err, 1, p_str_err.as_deref_mut(), b"Error: No V3000 Bond block start marker")?;
            return Ok(err);
        }
        remove_one_lf(heap, line)?;
        {
            let v3000 = heap.slice_mut(ctab.v3000)?.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)?;
            v3000.n_haptic_bonds = 0;
            v3000.n_non_haptic_bonds = 0;
        }

        let mut i = 0_i32;
        while i < ctab.n_bonds {
            let mut is_haptic = 0_i32;
            inchi_strbuf_reset(heap, Some(&mut tmpin.s))?;
            let nc = get_V3000_input_line_to_strbuf(heap, &mut tmpin.s, inp_file.as_deref_mut())?;
            p = if nc < 1 { SourceMutPointer::null() } else { tmpin.s.pStr };
            line = p;
            if p.is_null() {
                if err == 0 {
                    mol_fmt3_treat_error(&mut err, 2, p_str_err.as_deref_mut(), b"Cannot read V3000 bond block line")?;
                }
                break;
            }
            remove_one_lf(heap, line)?;
            if err != 0 {
                if mol_fmt3_c_string_eq(heap.slice(line.as_const())?, &SD_FMT_END_OF_DATA[..SD_FMT_END_OF_DATA.len() - 1]) {
                    err = err.wrapping_abs().wrapping_neg();
                    break;
                }
                i = i.wrapping_add(1);
                continue;
            }

            if !ctab.bonds.is_null() {
                let mut index = 0_i32;
                let mut atnum1 = -1_i16;
                let mut atnum2 = -1_i16;
                let mut bond_type = 0_i8;
                let mut stereo = 0_i8;
                let mut failed = 0_i32;
                let v3000 = heap.slice(ctab.v3000.as_const())?.first().cloned().ok_or(SourceHeapError::PointerOutOfBounds)?;
                let n_orig_at = v3000.n_non_star_atoms.wrapping_add(v3000.n_star_atoms);

                if MolfileV3000ReadField(heap, MolfileV3000FieldData::Int(&mut index), i32::from(MOL_FMT_INT_DATA), &mut p)? < 0 {
                    failed = 1;
                } else if MolfileV3000ReadField(heap, MolfileV3000FieldData::Char(&mut bond_type), i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)? < 0 {
                    failed = 1;
                } else if MolfileV3000ReadField(heap, MolfileV3000FieldData::Short(&mut atnum1), i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0 {
                    failed = 1;
                } else if MolfileV3000ReadField(heap, MolfileV3000FieldData::Short(&mut atnum2), i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0 {
                    failed = 1;
                }
                let orig = heap.slice(v3000.atom_index_orig.as_const())?.to_vec();
                let fin = heap.slice(v3000.atom_index_fin.as_const())?.to_vec();
                atnum1 = get_actual_atom_number(i32::from(atnum1), n_orig_at, Some(&orig), Some(&fin))? as i16;
                atnum2 = get_actual_atom_number(i32::from(atnum2), n_orig_at, Some(&orig), Some(&fin))? as i16;
                if atnum1 < 0 && atnum2 < 0 { failed = 1; }
                if failed != 0 {
                    if err == 0 {
                        mol_fmt3_treat_error(&mut err, 4, p_str_err.as_deref_mut(), b"Cannot interpret V3000 bond block line:")?;
                        let _ = dotify_non_printable_chars(heap, line)?;
                        let source_line = heap.slice(line.as_const())?.to_vec();
                        let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(&source_line))?;
                    }
                    if mol_fmt3_c_string_eq(heap.slice(line.as_const())?, &SD_FMT_END_OF_DATA[..SD_FMT_END_OF_DATA.len() - 1]) {
                        err = err.wrapping_abs().wrapping_neg();
                        break;
                    }
                }

                while !p.is_null() && MolfileV3000ReadKeyword(heap, field, &mut p)? > 0 {
                    let key = heap.slice(field.as_const())?;
                    let mut keyword_failed = 0_i32;
                    if mol_fmt3_c_string_eq(key, b"CFG") {
                        if MolfileV3000ReadField(heap, MolfileV3000FieldData::Char(&mut stereo), i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)? < 0 {
                            keyword_failed = 1;
                        } else if stereo == 2 { stereo = 4; } else if stereo == 3 { stereo = 6; }
                    } else if mol_fmt3_c_string_eq(key, b"TOPO") || mol_fmt3_c_string_eq(key, b"RXCTR") || mol_fmt3_c_string_eq(key, b"STBOX") {
                        let mut ignored = 0_i32;
                        let _ = MolfileV3000ReadField(heap, MolfileV3000FieldData::Int(&mut ignored), i32::from(MOL_FMT_INT_DATA), &mut p)?;
                    } else if mol_fmt3_c_string_eq(key, b"ENDPTS") {
                        let mut num_list = SourceMutPointer::null();
                        if MolfileV3000ReadHapticBond(heap, Some(&*ctab), &mut p, &mut num_list, p_str_err.as_deref_mut())? < 0 || num_list.is_null() {
                            keyword_failed = 1;
                        } else {
                            let mut existent_atom = atnum1;
                            if existent_atom < 0 { existent_atom = atnum2; }
                            if existent_atom < 0 {
                                keyword_failed = 1;
                            } else {
                                let nnum = heap.slice(num_list.as_const())?[2];
                                heap.slice_mut(num_list)?[1] = i32::from(existent_atom);
                                let mut k = 3_i32;
                                while k < nnum {
                                    let position = usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                    let atom = heap.slice(num_list.as_const())?[position];
                                    let actual = get_actual_atom_number(atom, n_orig_at, Some(&orig), Some(&fin))?;
                                    heap.slice_mut(num_list)?[position] = actual;
                                    k = k.wrapping_add(1);
                                }
                                let haptic = heap.slice(ctab.v3000.as_const())?.first().ok_or(SourceHeapError::PointerOutOfBounds)?.haptic_bonds;
                                let mut lists = heap.slice(haptic.as_const())?.first().cloned().ok_or(SourceHeapError::PointerOutOfBounds)?;
                                let append = NumLists_Append(heap, Some(&mut lists), num_list)?;
                                *heap.slice_mut(haptic)?.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)? = lists;
                                if append < 0 { keyword_failed = 1; } else { is_haptic = 1; }
                            }
                        }
                    } else if mol_fmt3_c_string_eq(key, b"DISP") || mol_fmt3_c_string_eq(key, b"ATTACH") {
                        let mut ignored = [0_i8; MOL_FMT_V3000_MAXFIELDLEN as usize];
                        let _ = MolfileV3000ReadField(heap, MolfileV3000FieldData::String(&mut ignored), i32::from(MOL_FMT_STRING_DATA), &mut p)?;
                    }
                    if keyword_failed != 0 {
                        if err == 0 {
                            mol_fmt3_treat_error(&mut err, 4, p_str_err.as_deref_mut(), b"Cannot interpret V3000 bond block line:")?;
                            let _ = dotify_non_printable_chars(heap, line)?;
                            let source_line = heap.slice(line.as_const())?.to_vec();
                            let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(&source_line))?;
                        }
                        if mol_fmt3_c_string_eq(heap.slice(line.as_const())?, &SD_FMT_END_OF_DATA[..SD_FMT_END_OF_DATA.len() - 1]) {
                            err = err.wrapping_abs().wrapping_neg();
                            break;
                        }
                    }
                }

                if is_haptic != 0 {
                    let v3000 = heap.slice(ctab.v3000.as_const())?.first().cloned().ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let lists = heap.slice(v3000.haptic_bonds.as_const())?.first().cloned().ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let ii = usize::try_from(v3000.n_haptic_bonds).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let num_list = *heap.slice(lists.lists.as_const())?.get(ii).ok_or(SourceHeapError::PointerOutOfBounds)?;
                    heap.slice_mut(num_list)?[0] = i32::from(bond_type);
                    heap.slice_mut(ctab.v3000)?[0].n_haptic_bonds = v3000.n_haptic_bonds.wrapping_add(1);
                    i = i.wrapping_add(1);
                    continue;
                }
                let non_haptic = heap.slice(ctab.v3000.as_const())?.first().ok_or(SourceHeapError::PointerOutOfBounds)?.n_non_haptic_bonds;
                let ii = usize::try_from(non_haptic).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let bond = heap.slice_mut(ctab.bonds)?.get_mut(ii).ok_or(SourceHeapError::PointerOutOfBounds)?;
                bond.atnum1 = atnum1;
                bond.atnum2 = atnum2;
                bond.bond_type = bond_type;
                bond.bond_stereo = stereo;
                heap.slice_mut(ctab.v3000)?[0].n_non_haptic_bonds = non_haptic.wrapping_add(1);
                let _ = index;
            }
            i = i.wrapping_add(1);
        }

        let v3000 = heap.slice(ctab.v3000.as_const())?.first().cloned().ok_or(SourceHeapError::PointerOutOfBounds)?;
        if v3000.n_haptic_bonds != 0 {
            let _ = mol_fmt3_add_message(p_str_err.as_deref_mut(), b"V3000 haptic bonds read/stored but ignored")?;
            ctab.n_bonds = v3000.n_non_haptic_bonds;
        }
        inchi_strbuf_reset(heap, Some(&mut tmpin.s))?;
        let nc = get_V3000_input_line_to_strbuf(heap, &mut tmpin.s, inp_file.as_deref_mut())?;
        p = if nc < 1 { SourceMutPointer::null() } else { tmpin.s.pStr };
        line = p;
        if p.is_null() || !mol_fmt3_c_string_eq(heap.slice(p.as_const())?, b"END BOND") {
            mol_fmt3_treat_error(&mut err, 1, p_str_err.as_deref_mut(), b"Error: No V3000 Bond block end marker")?;
            return Ok(err);
        }
        remove_one_lf(heap, line)?;
        Ok(err)
    })();
    let field_cleanup = heap.free(field);
    let stream_cleanup = inchi_ios_close(heap, Some(&mut tmpin));
    match (result, field_cleanup, stream_cleanup) {
        (Err(error), _, _) | (Ok(_), Err(error), _) | (Ok(_), Ok(()), Err(error)) => Err(error),
        (Ok(value), Ok(()), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileV3000ReadSGroup(
    heap: &mut SourceHeap,
    _ctab: Option<&mut MOL_FMT_CTAB>,
    mut inp_file: Option<&mut INCHI_IOSTREAM>,
    _err: i32,
    _p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:561 MolfileV3000ReadSGroup
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap stream buffering adds overhead.
    /*
int MolfileV3000ReadSGroup(MOL_FMT_CTAB *ctab,
                           INCHI_IOSTREAM *inp_file,
                           int err,
                           char *pStrErr)
{
    int nc;
    char *p = NULL, *line = NULL;
    INCHI_IOSTREAM tmpin;
    INCHI_IOS_STRING *pin = &tmpin.s;

    inchi_ios_init(&tmpin, INCHI_IOS_TYPE_STRING, NULL);
    /*p = inchi_fgetsLf_V3000( line, inp_file );*/

    while (1)
    {
        nc = get_V3000_input_line_to_strbuf(pin, inp_file);

        if (nc < 1)
        {
            p = NULL;
        }
        else
        {
            p = line = pin->pStr;
            remove_one_lf(line);
        }
        if (p && !strcmp(p, "END SGROUP"))
        {
            inchi_ios_close(&tmpin); /* ricrogz: fixing memory leak */
            return 0;
        }
    }

    /* if (  !p  || strcmp(p, "END SGROUP") ) */
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Error: No V3000 SGroup end marker");
    }

err_fin:

    inchi_ios_close(&tmpin); /* ricrogz: fixing memory leak */
    return err;
}
    */
    // END INCHI C FUNCTION: MolfileV3000ReadSGroup
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadSGroup
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; the post-loop TREAT_ERR block is unreachable.
    // INCHI✔️❌: EOF does not terminate the active source loop; terminating modeled inputs must contain an exact END SGROUP marker.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadSGroup

    let mut tmpin = INCHI_IOSTREAM::default();
    inchi_ios_init(Some(&mut tmpin), INCHI_IOS_TYPE_STRING as i32, SourceMutPointer::null())?;
    loop {
        let nc = get_V3000_input_line_to_strbuf(heap, &mut tmpin.s, inp_file.as_deref_mut())?;
        if nc >= 1 {
            let line = tmpin.s.pStr;
            remove_one_lf(heap, line)?;
            if mol_fmt3_c_string_eq(heap.slice(line.as_const())?, b"END SGROUP") {
                inchi_ios_close(heap, Some(&mut tmpin))?;
                return Ok(0);
            }
        }
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileV3000Read3DBlock(
    heap: &mut SourceHeap,
    _ctab: Option<&mut MOL_FMT_CTAB>,
    mut inp_file: Option<&mut INCHI_IOSTREAM>,
    mut err: i32,
    mut p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:608 MolfileV3000Read3DBlock
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap stream buffering adds overhead.
    /*
int MolfileV3000Read3DBlock(MOL_FMT_CTAB *ctab,
                            INCHI_IOSTREAM *inp_file,
                            int err,
                            char *pStrErr)
{
    int nc;
    char *p = NULL, *line = NULL;
    INCHI_IOSTREAM tmpin;
    INCHI_IOS_STRING *pin = &tmpin.s;
    inchi_ios_init(&tmpin, INCHI_IOS_TYPE_STRING, NULL);
    /*p = inchi_fgetsLf_V3000( line, inp_file );*/

    nc = get_V3000_input_line_to_strbuf(pin, inp_file);

    if (nc < 1)
    {
        p = NULL;
    }
    else
    {
        p = line = pin->pStr;
    }
    remove_one_lf(line);

    if (!p || strcmp(p, "END OBJ3D"))
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Error: No V3000 3DBlock end marker");
    }
    goto err_fin;

err_fin:

    inchi_ios_close(&tmpin); /* ricrogz: fixing memory leak */
    return err;
}
    */
    // END INCHI C FUNCTION: MolfileV3000Read3DBlock
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000Read3DBlock
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; TREAT_ERR preserves a pre-existing nonzero error and always appends the message.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000Read3DBlock

    let mut tmpin = INCHI_IOSTREAM::default();
    inchi_ios_init(Some(&mut tmpin), INCHI_IOS_TYPE_STRING as i32, SourceMutPointer::null())?;
    let result = (|| -> Result<i32, SourceHeapError> {
        let nc = get_V3000_input_line_to_strbuf(heap, &mut tmpin.s, inp_file.as_deref_mut())?;
        let p = if nc < 1 { SourceMutPointer::null() } else { tmpin.s.pStr };
        remove_one_lf(heap, p)?;
        if p.is_null() || !mol_fmt3_c_string_eq(heap.slice(p.as_const())?, b"END OBJ3D") {
            mol_fmt3_treat_error(&mut err, 1, p_str_err.as_deref_mut(), b"Error: No V3000 3DBlock end marker")?;
        }
        Ok(err)
    })();
    let cleanup = inchi_ios_close(heap, Some(&mut tmpin));
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileV3000ReadCollections(
    heap: &mut SourceHeap,
    ctab: &mut MOL_FMT_CTAB,
    mut inp_file: Option<&mut INCHI_IOSTREAM>,
    mut err: i32,
    mut p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:647 MolfileV3000ReadCollections
    // INCHI✔️❌: complete source frame follows verbatim; modeled stack storage, descriptor copies, and checked SourceHeap access add overhead.
    /*
int MolfileV3000ReadCollections(MOL_FMT_CTAB *ctab,
                                INCHI_IOSTREAM *inp_file,
                                int err,
                                char *pStrErr)
{
    char field[MOL_FMT_V3000_MAXFIELDLEN];
    const int max_field_len = sizeof(field);
    int nread, len, n_coll = 0;
    int failed = 0;
    int nc;
    char *p = NULL, *line = NULL;
    INCHI_IOSTREAM tmpin;
    INCHI_IOS_STRING *pin = &tmpin.s;

    inchi_ios_init(&tmpin, INCHI_IOS_TYPE_STRING, NULL);
    /*p = inchi_fgetsLf_V3000( line, inp_file );*/

    nc = get_V3000_input_line_to_strbuf(pin, inp_file);

    if (nc < 1)
    {
        p = NULL;
    }
    else
    {
        p = line = pin->pStr;
    }
    remove_one_lf(line);

    while (p && strcmp(p, "END COLLECTION"))
    {
        int stereo_kind = MOL_FMT_V3000_STENON;
        /* stereo collection of interest */
        NUM_LISTS *ste_coll = NULL;

        nread = read_upto_delim(&p, field, max_field_len, "/");
        if (nread < 6)
        {
            failed = 1;
            break;
        }
        if (strcmp(field, "MDLV30"))
        {
            failed = 1;
            break;
        }

        nread = read_upto_delim(&p, field, max_field_len, "1234567890 \t\n\v\f\r"); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        if (!strcmp(field, "/STEABS"))
        {
            n_coll = 1;
            stereo_kind = MOL_FMT_V3000_STEABS;
            ste_coll = ctab->v3000->steabs;
        }
        else if (!strcmp(field, "/STEREL"))
        {
            /* get number of collection */
            if (0 > MolfileV3000ReadField(&n_coll, MOL_FMT_CHAR_INT_DATA, &p))
            {
                failed = 1;
                break;
            }
            stereo_kind = MOL_FMT_V3000_STEREL;
            ste_coll = ctab->v3000->sterel;
        }
        else if (!strcmp(field, "/STERAC"))
        {
            /* get number of collection */
            if (0 > MolfileV3000ReadField(&n_coll, MOL_FMT_CHAR_INT_DATA, &p))
            {
                failed = 1;
                break;
            }
            stereo_kind = MOL_FMT_V3000_STERAC;
            ste_coll = ctab->v3000->sterac;
        }
        else
        {
            ;
        }

        if (stereo_kind != MOL_FMT_V3000_STENON)
        /* currently skip non-stereo collections */
        {
            /* consume atoms= */
            if ((len = MolfileV3000ReadKeyword(field, &p) > 0)) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
            {
                if (!strcmp(field, "ATOMS"))
                {
                    int res, *num_list = NULL;

                    if (0 > MolfileV3000ReadStereoCollection(ctab, &p, &num_list, pStrErr))
                    {
                        failed = 1;
                    }
                    else if (!num_list)
                    {
                        failed = 1;
                    }
                    else
                    {
                        int k, nnum;
                        num_list[0] = n_coll;
                        nnum = num_list[1];
                        for (k = 2; k < nnum; k++)
                        {
                            num_list[k] =
                                get_actual_atom_number(num_list[k],
                                                       ctab->v3000->n_non_star_atoms + ctab->v3000->n_star_atoms,
                                                       ctab->v3000->atom_index_orig,
                                                       ctab->v3000->atom_index_fin);
                        }
                        res = NumLists_Append(ste_coll, num_list);
                        if (res < 0)
                        {
                            failed = 1;
                        }
                        else
                        {
                            if (stereo_kind == MOL_FMT_V3000_STEABS)
                            {
                                ctab->v3000->n_steabs++;
                                ctab->v3000->n_collections++;
                            }
                            else if (stereo_kind == MOL_FMT_V3000_STEREL)
                            {
                                ctab->v3000->n_sterel++;
                                ctab->v3000->n_collections++;
                            }
                            else if (stereo_kind == MOL_FMT_V3000_STERAC)
                            {
                                ctab->v3000->n_sterac++;
                                ctab->v3000->n_collections++;
                            }
                        }
                    }
                }
            }
            else
            {
                failed = 1;
                break;
            }
        }

        /*next_line:*/
        /*p = inchi_fgetsLf_V3000( line, inp_file );*/
        inchi_strbuf_reset(pin);

        nc = get_V3000_input_line_to_strbuf(pin, inp_file);

        if (nc < 1)
        {
            p = NULL;
        }
        else
        {
            p = line = pin->pStr;
            remove_one_lf(line);
        }
    }

    if (failed)
    {
        /*p = inchi_fgetsLf_V3000( line, inp_file );*/
        inchi_strbuf_reset(pin);
        line = NULL; /* Reset line pointer since buffer was freed */

        nc = get_V3000_input_line_to_strbuf(pin, inp_file);

        if (nc < 1)
        {
            p = NULL;
        }
        else
        {
            p = line = pin->pStr;
            remove_one_lf(line);
        }
    }

    if (!p)
    {
        failed = 1;
    }

    if (failed)
    {
        err = 7;
        TREAT_ERR(err, 7, "Cannot interpret V3000 collection line(s)"); /* djb-rwth: addressing coverity ID #499531 -- TREAT_ERR properly used */
        if (line)
        {
            dotify_non_printable_chars(line);
            AddErrorMessage(pStrErr, line);
        }
        goto err_fin;
    }

    /* Error: No V3000 Collection end marker */
    if (ctab->v3000->n_steabs ||
        ctab->v3000->n_sterel ||
        ctab->v3000->n_sterac)
    {
        AddErrorMessage(pStrErr, "V3000 enhanced stereo read/stored but ignored");
    }

err_fin:

    inchi_ios_close(&tmpin); /* ricrogz: fixing memory leak */
    return err;
}
    */
    // END INCHI C FUNCTION: MolfileV3000ReadCollections
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadCollections
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux little-endian LP64; MOL_FMT_CHAR_INT_DATA writes only the low byte of n_coll.
    // INCHI✔️❌: The source atom-remapping loop is k = 2; k < nnum, and a failed line causes one additional logical line to be consumed before error reporting.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadCollections

    let field = heap.allocate_model_storage(vec![0_i8; MOL_FMT_V3000_MAXFIELDLEN as usize])?;
    let mut tmpin = INCHI_IOSTREAM::default();
    inchi_ios_init(Some(&mut tmpin), INCHI_IOS_TYPE_STRING as i32, SourceMutPointer::null())?;
    let result = (|| -> Result<i32, SourceHeapError> {
        let mut n_coll = 0_i32;
        let mut failed = 0_i32;
        let nc = get_V3000_input_line_to_strbuf(heap, &mut tmpin.s, inp_file.as_deref_mut())?;
        let mut p = if nc < 1 { SourceMutPointer::null() } else { tmpin.s.pStr };
        let mut line = p;
        remove_one_lf(heap, line)?;
        while !p.is_null() && !mol_fmt3_c_string_eq(heap.slice(p.as_const())?, b"END COLLECTION") {
            let mut stereo_kind = MOL_FMT_V3000_STENON;
            let mut ste_coll = SourceMutPointer::null();
            let slash = [b'/' as i8, 0];
            if read_upto_delim(heap, &mut p, field, MOL_FMT_V3000_MAXFIELDLEN as i32, Some(&slash))? < 6
                || !mol_fmt3_c_string_eq(heap.slice(field.as_const())?, b"MDLV30")
            {
                failed = 1;
                break;
            }
            let delimiters = b"1234567890 \t\n\x0b\x0c\r\0".map(|byte| byte as i8);
            let _ = read_upto_delim(heap, &mut p, field, MOL_FMT_V3000_MAXFIELDLEN as i32, Some(&delimiters))?;
            let v3000 = heap.slice(ctab.v3000.as_const())?.first().cloned().ok_or(SourceHeapError::PointerOutOfBounds)?;
            if mol_fmt3_c_string_eq(heap.slice(field.as_const())?, b"/STEABS") {
                n_coll = 1;
                stereo_kind = MOL_FMT_V3000_STEABS as i32;
                ste_coll = v3000.steabs;
            } else if mol_fmt3_c_string_eq(heap.slice(field.as_const())?, b"/STEREL")
                || mol_fmt3_c_string_eq(heap.slice(field.as_const())?, b"/STERAC")
            {
                let is_relative = mol_fmt3_c_string_eq(heap.slice(field.as_const())?, b"/STEREL");
                let mut low_byte = n_coll as i8;
                if MolfileV3000ReadField(heap, MolfileV3000FieldData::Char(&mut low_byte), i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)? < 0 {
                    failed = 1;
                    break;
                }
                n_coll = (n_coll & !0xff) | i32::from(low_byte as u8);
                if is_relative {
                    stereo_kind = MOL_FMT_V3000_STEREL as i32;
                    ste_coll = v3000.sterel;
                } else {
                    stereo_kind = MOL_FMT_V3000_STERAC as i32;
                    ste_coll = v3000.sterac;
                }
            }

            if stereo_kind != MOL_FMT_V3000_STENON {
                if MolfileV3000ReadKeyword(heap, field, &mut p)? > 0 {
                    if mol_fmt3_c_string_eq(heap.slice(field.as_const())?, b"ATOMS") {
                        let mut num_list = SourceMutPointer::null();
                        if MolfileV3000ReadStereoCollection(heap, Some(&*ctab), &mut p, &mut num_list, p_str_err.as_deref_mut())? < 0 || num_list.is_null() {
                            failed = 1;
                        } else {
                            heap.slice_mut(num_list)?[0] = n_coll;
                            let nnum = heap.slice(num_list.as_const())?[1];
                            let orig = heap.slice(v3000.atom_index_orig.as_const())?.to_vec();
                            let fin = heap.slice(v3000.atom_index_fin.as_const())?.to_vec();
                            let n_orig = v3000.n_non_star_atoms.wrapping_add(v3000.n_star_atoms);
                            let mut k = 2_i32;
                            while k < nnum {
                                let position = usize::try_from(k).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                let atom = heap.slice(num_list.as_const())?[position];
                                heap.slice_mut(num_list)?[position] = get_actual_atom_number(atom, n_orig, Some(&orig), Some(&fin))?;
                                k = k.wrapping_add(1);
                            }
                            let mut lists = heap.slice(ste_coll.as_const())?.first().cloned().ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let append = NumLists_Append(heap, Some(&mut lists), num_list)?;
                            *heap.slice_mut(ste_coll)?.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)? = lists;
                            if append < 0 {
                                failed = 1;
                            } else {
                                let current = heap.slice(ctab.v3000.as_const())?.first().cloned().ok_or(SourceHeapError::PointerOutOfBounds)?;
                                let target = heap.slice_mut(ctab.v3000)?.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)?;
                                if stereo_kind == MOL_FMT_V3000_STEABS as i32 { target.n_steabs = current.n_steabs.wrapping_add(1); }
                                else if stereo_kind == MOL_FMT_V3000_STEREL as i32 { target.n_sterel = current.n_sterel.wrapping_add(1); }
                                else if stereo_kind == MOL_FMT_V3000_STERAC as i32 { target.n_sterac = current.n_sterac.wrapping_add(1); }
                                target.n_collections = current.n_collections.wrapping_add(1);
                            }
                        }
                    }
                } else {
                    failed = 1;
                    break;
                }
            }
            inchi_strbuf_reset(heap, Some(&mut tmpin.s))?;
            let nc = get_V3000_input_line_to_strbuf(heap, &mut tmpin.s, inp_file.as_deref_mut())?;
            p = if nc < 1 { SourceMutPointer::null() } else { tmpin.s.pStr };
            line = p;
            if !line.is_null() { remove_one_lf(heap, line)?; }
        }
        if failed != 0 {
            inchi_strbuf_reset(heap, Some(&mut tmpin.s))?;
            line = SourceMutPointer::null();
            let nc = get_V3000_input_line_to_strbuf(heap, &mut tmpin.s, inp_file.as_deref_mut())?;
            p = if nc < 1 { SourceMutPointer::null() } else { tmpin.s.pStr };
            line = p;
            if !line.is_null() { remove_one_lf(heap, line)?; }
        }
        if p.is_null() { failed = 1; }
        if failed != 0 {
            err = 7;
            mol_fmt3_treat_error(&mut err, 7, p_str_err.as_deref_mut(), b"Cannot interpret V3000 collection line(s)")?;
            if !line.is_null() {
                let _ = dotify_non_printable_chars(heap, line)?;
                let source_line = heap.slice(line.as_const())?.to_vec();
                let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(&source_line))?;
            }
            return Ok(err);
        }
        let v3000 = heap.slice(ctab.v3000.as_const())?.first().ok_or(SourceHeapError::PointerOutOfBounds)?;
        if v3000.n_steabs != 0 || v3000.n_sterel != 0 || v3000.n_sterac != 0 {
            let _ = mol_fmt3_add_message(p_str_err.as_deref_mut(), b"V3000 enhanced stereo read/stored but ignored")?;
        }
        Ok(err)
    })();
    let field_cleanup = heap.free(field);
    let stream_cleanup = inchi_ios_close(heap, Some(&mut tmpin));
    match (result, field_cleanup, stream_cleanup) {
        (Err(error), _, _) | (Ok(_), Err(error), _) | (Ok(_), Ok(()), Err(error)) => Err(error),
        (Ok(value), Ok(()), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
pub(crate) fn get_actual_atom_number(
    index: i32,
    n: i32,
    orig: Option<&[i32]>,
    fin: Option<&[i32]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1585 get_actual_atom_number
    // INCHI✔️✔️: complete source frame follows verbatim; Rust performs the same linear scan without allocation.
    /*
int get_actual_atom_number(int index, int n, int *orig, int *fin)
{
    int i;
    for (i = 0; i < n; i++)
    {
        if (orig[i] == index)
        {
            return fin[i];
        }
    }

    return -1;
}
    */
    // END INCHI C FUNCTION: get_actual_atom_number
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: get_actual_atom_number
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; int is signed 32-bit.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: get_actual_atom_number

    if n <= 0 {
        return Ok(-1);
    }
    let count = usize::try_from(n).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let orig = orig.ok_or(SourceHeapError::NullPointer)?;
    let fin = fin.ok_or(SourceHeapError::NullPointer)?;
    let orig = orig.get(..count).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let fin = fin.get(..count).ok_or(SourceHeapError::PointerOutOfBounds)?;
    for i in 0..count {
        if orig[i] == index {
            return Ok(fin[i]);
        }
    }
    Ok(-1)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileV3000ReadTailOfCTAB(
    heap: &mut SourceHeap,
    ctab: &mut MOL_FMT_CTAB,
    mut inp_file: Option<&mut INCHI_IOSTREAM>,
    mut err: i32,
    mut p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1602 MolfileV3000ReadTailOfCTAB
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap stream and descriptor access add overhead.
    /*
int MolfileV3000ReadTailOfCTAB(MOL_FMT_CTAB *ctab,
                               INCHI_IOSTREAM *inp_file,
                               int err,
                               char *pStrErr)
{
    int retcode = err;
    int nc;
    char *p = NULL, *line = NULL;
    INCHI_IOSTREAM tmpin;
    INCHI_IOS_STRING *pin = &tmpin.s;
    inchi_ios_init(&tmpin, INCHI_IOS_TYPE_STRING, NULL);

    /*p = inchi_fgetsLf_V3000( line, inp_file );*/

    nc = get_V3000_input_line_to_strbuf(pin, inp_file);

    if (nc < 1)
    {
        p = NULL;
    }
    else
    {
        p = line = pin->pStr;
    }
    remove_one_lf(line);

    if (p && !strcmp(p, "BEGIN SGROUP"))
    {
        retcode = MolfileV3000ReadSGroup(ctab, inp_file, retcode, pStrErr);
        if (retcode)
        {
            retcode += 70;
            TREAT_ERR_AND_FIN(retcode, 1, err_fin, pStrErr);
        }
        /*p = inchi_fgetsLf_V3000( line, inp_file );*/
        inchi_strbuf_reset(pin);
        nc = get_V3000_input_line_to_strbuf(pin, inp_file);
        if (nc < 1)
        {
            p = NULL;
        }
        else
        {
            p = line = pin->pStr;
            remove_one_lf(line);
        }
    }

    if (p && !strcmp(p, "BEGIN OBJ3D"))
    {
        retcode = MolfileV3000Read3DBlock(ctab, inp_file, retcode, pStrErr);
        if (retcode)
        {
            retcode += 70;
            TREAT_ERR_AND_FIN(retcode, 1, err_fin, pStrErr);
        }
        /*p = inchi_fgetsLf_V3000( line, inp_file );*/
        inchi_strbuf_reset(pin);

        nc = get_V3000_input_line_to_strbuf(pin, inp_file);
        if (nc < 1)
        {
            p = NULL;
        }
        else
        {
            p = line = pin->pStr;
            remove_one_lf(line);
        }
    }

    while (p && !strcmp(p, "LINKNODE"))
    {
        /* skip for now */
        /*p = inchi_fgetsLf_V3000( line, inp_file );*/
        inchi_strbuf_reset(pin);

        nc = get_V3000_input_line_to_strbuf(pin, inp_file);
        if (nc < 1)
        {
            p = NULL;
        }
        else
        {
            p = line = pin->pStr;
            remove_one_lf(line);
        }
    }

    /* Collections */
    while (p && !strcmp(p, "BEGIN COLLECTION"))
    {
        retcode = MolfileV3000ReadCollections(ctab, inp_file, retcode, pStrErr);
        if (retcode)
        {
            retcode += 70;
            TREAT_ERR_AND_FIN(retcode, 1, err_fin, pStrErr);
        }
        /*p = inchi_fgetsLf_V3000( line, inp_file );*/
        inchi_strbuf_reset(pin);

        nc = get_V3000_input_line_to_strbuf(pin, inp_file);

        if (nc < 1)
        {
            p = NULL;
        }
        else
        {
            p = line = pin->pStr;
            remove_one_lf(line);
        }
    }

    if (!p || strcmp(p, "END CTAB"))
    {
        TREAT_ERR_AND_FIN(err, 1, err_fin, "Error: No V3000 CTAB end marker");
    }

    remove_one_lf(line);

err_fin:
    inchi_strbuf_close(pin);

    return err;
}
    */
    // END INCHI C FUNCTION: MolfileV3000ReadTailOfCTAB
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadTailOfCTAB
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; TREAT_ERR_AND_FIN preserves nonzero errors, appends through AddErrorMessage, and jumps to cleanup.
    // INCHI✔️❌: A child error changes only retcode; the active source still returns the original err parameter after closing the temporary string buffer.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000ReadTailOfCTAB

    let mut retcode = err;
    let mut tmpin = INCHI_IOSTREAM::default();
    inchi_ios_init(Some(&mut tmpin), INCHI_IOS_TYPE_STRING as i32, SourceMutPointer::null())?;
    let result = (|| -> Result<i32, SourceHeapError> {
        let nc = get_V3000_input_line_to_strbuf(heap, &mut tmpin.s, inp_file.as_deref_mut())?;
        let mut p = if nc < 1 { SourceMutPointer::null() } else { tmpin.s.pStr };
        let mut line = p;
        remove_one_lf(heap, line)?;

        if !p.is_null() && mol_fmt3_c_string_eq(heap.slice(p.as_const())?, b"BEGIN SGROUP") {
            retcode = MolfileV3000ReadSGroup(heap, Some(&mut *ctab), inp_file.as_deref_mut(), retcode, p_str_err.as_deref_mut())?;
            if retcode != 0 {
                retcode = retcode.wrapping_add(70);
                let message = p_str_err.as_deref().map(<[i8]>::to_vec);
                let _ = AddErrorMessage(p_str_err.as_deref_mut(), message.as_deref())?;
                return Ok(err);
            }
            inchi_strbuf_reset(heap, Some(&mut tmpin.s))?;
            let nc = get_V3000_input_line_to_strbuf(heap, &mut tmpin.s, inp_file.as_deref_mut())?;
            if nc < 1 {
                p = SourceMutPointer::null();
            } else {
                p = tmpin.s.pStr;
                line = p;
                remove_one_lf(heap, line)?;
            }
        }

        if !p.is_null() && mol_fmt3_c_string_eq(heap.slice(p.as_const())?, b"BEGIN OBJ3D") {
            retcode = MolfileV3000Read3DBlock(heap, Some(&mut *ctab), inp_file.as_deref_mut(), retcode, p_str_err.as_deref_mut())?;
            if retcode != 0 {
                retcode = retcode.wrapping_add(70);
                let message = p_str_err.as_deref().map(<[i8]>::to_vec);
                let _ = AddErrorMessage(p_str_err.as_deref_mut(), message.as_deref())?;
                return Ok(err);
            }
            inchi_strbuf_reset(heap, Some(&mut tmpin.s))?;
            let nc = get_V3000_input_line_to_strbuf(heap, &mut tmpin.s, inp_file.as_deref_mut())?;
            if nc < 1 {
                p = SourceMutPointer::null();
            } else {
                p = tmpin.s.pStr;
                line = p;
                remove_one_lf(heap, line)?;
            }
        }

        while !p.is_null() && mol_fmt3_c_string_eq(heap.slice(p.as_const())?, b"LINKNODE") {
            inchi_strbuf_reset(heap, Some(&mut tmpin.s))?;
            let nc = get_V3000_input_line_to_strbuf(heap, &mut tmpin.s, inp_file.as_deref_mut())?;
            if nc < 1 {
                p = SourceMutPointer::null();
            } else {
                p = tmpin.s.pStr;
                line = p;
                remove_one_lf(heap, line)?;
            }
        }

        while !p.is_null() && mol_fmt3_c_string_eq(heap.slice(p.as_const())?, b"BEGIN COLLECTION") {
            retcode = MolfileV3000ReadCollections(heap, ctab, inp_file.as_deref_mut(), retcode, p_str_err.as_deref_mut())?;
            if retcode != 0 {
                retcode = retcode.wrapping_add(70);
                let message = p_str_err.as_deref().map(<[i8]>::to_vec);
                let _ = AddErrorMessage(p_str_err.as_deref_mut(), message.as_deref())?;
                return Ok(err);
            }
            inchi_strbuf_reset(heap, Some(&mut tmpin.s))?;
            let nc = get_V3000_input_line_to_strbuf(heap, &mut tmpin.s, inp_file.as_deref_mut())?;
            if nc < 1 {
                p = SourceMutPointer::null();
            } else {
                p = tmpin.s.pStr;
                line = p;
                remove_one_lf(heap, line)?;
            }
        }

        if p.is_null() || !mol_fmt3_c_string_eq(heap.slice(p.as_const())?, b"END CTAB") {
            mol_fmt3_treat_error(&mut err, 1, p_str_err.as_deref_mut(), b"Error: No V3000 CTAB end marker")?;
            return Ok(err);
        }
        remove_one_lf(heap, line)?;
        Ok(err)
    })();
    let cleanup = inchi_strbuf_close(heap, Some(&mut tmpin.s));
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileV3000Init(
    heap: &mut SourceHeap,
    ctab: &mut MOL_FMT_CTAB,
    mut p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:68 MolfileV3000Init
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap descriptor access adds overhead.
    /*
int MolfileV3000Init(MOL_FMT_CTAB *ctab,
                     char *pStrErr)
{
    int ret = 0;
    int i;

    /* STAR ATOMS */
    ctab->v3000->n_star_atoms = 0;
    ctab->v3000->n_non_star_atoms = 0;

    if (ctab->n_atoms)
    {
        ctab->v3000->atom_index_orig = (int *)inchi_calloc(ctab->n_atoms, sizeof(int));
        ctab->v3000->atom_index_fin = (int *)inchi_calloc(ctab->n_atoms, sizeof(int));
        if (ctab->v3000->atom_index_orig && ctab->v3000->atom_index_fin) /* djb-rwth: fixing a NULL pointer dereference */
        {
            for (i = 0; i < ctab->n_atoms; i++) /* protective */
            {
                ctab->v3000->atom_index_orig[i] = -1;
                ctab->v3000->atom_index_fin[i] = -1;
            }
        }
    }
    else
    {
        ctab->v3000->atom_index_orig = NULL;
        ctab->v3000->atom_index_fin = NULL;
    }

    /* HAPTIC BONDS */
    ctab->v3000->n_haptic_bonds = 0;
    ctab->v3000->haptic_bonds = (NUM_LISTS *)inchi_calloc(1, sizeof(NUM_LISTS));
    if (!ctab->v3000->haptic_bonds)
    {
        AddErrorMessage(pStrErr, "Out of RAM");
        return -1;
    }
    ret = NumLists_Alloc(ctab->v3000->haptic_bonds, 8);
    if (ret < 0)
    {
        ctab->v3000->haptic_bonds = NULL;
        AddErrorMessage(pStrErr, "Out of RAM");
        return -1;
    }

    /* STEABS */
    ctab->v3000->n_steabs = 0;
    ctab->v3000->steabs = (NUM_LISTS *)inchi_calloc(1, sizeof(NUM_LISTS));
    if (!ctab->v3000->steabs)
    {
        AddErrorMessage(pStrErr, "Out of RAM");
        return -1;
    }
    ret = NumLists_Alloc(ctab->v3000->steabs, 1);
    if (ret < 0)
    {
        ctab->v3000->steabs = NULL;
        AddErrorMessage(pStrErr, "Out of RAM");
        return -1;
    }
    /* STEREL */
    ctab->v3000->n_sterel = 0;
    ctab->v3000->sterel = (NUM_LISTS *)inchi_calloc(1, sizeof(NUM_LISTS));
    if (!ctab->v3000->sterel)
    {
        AddErrorMessage(pStrErr, "Out of RAM");
        return -1;
    }
    ret = NumLists_Alloc(ctab->v3000->sterel, 4);
    if (ret < 0)
    {
        ctab->v3000->sterel = NULL;
        AddErrorMessage(pStrErr, "Out of RAM");
        return -1;
    }
    /* STERAC */
    ctab->v3000->n_sterac = 0;
    ctab->v3000->sterac = (NUM_LISTS *)inchi_calloc(1, sizeof(NUM_LISTS));
    if (!ctab->v3000->sterac)
    {
        AddErrorMessage(pStrErr, "Out of RAM");
        return -1;
    }
    ret = NumLists_Alloc(ctab->v3000->sterac, 4);
    if (ret < 0)
    {
        ctab->v3000->sterac = NULL;
        AddErrorMessage(pStrErr, "Out of RAM");
        return -1;
    }

    return ret;
}
    */
    // END INCHI C FUNCTION: MolfileV3000Init
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000Init
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; int is signed 32-bit, sizeof(int) is four, and sizeof(NUM_LISTS) is 24.
    // INCHI✔️❌: inchi_calloc resolves to libc calloc; atom-index allocation failure is not reported, while list allocation failure appends Out of RAM and does not roll back earlier allocations.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileV3000Init

    macro_rules! calloc_or_null {
        ($type:ty, $count:expr, $size:expr) => {
            match inchi_calloc::<$type>(heap, $count, $size) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed)
                | Err(SourceHeapError::AllocationElementCountOutOfRange)
                | Err(SourceHeapError::AllocationSizeOverflow) => SourceMutPointer::null(),
                Err(error) => return Err(error),
            }
        };
    }

    heap.slice_mut(ctab.v3000)?[0].n_star_atoms = 0;
    heap.slice_mut(ctab.v3000)?[0].n_non_star_atoms = 0;
    if ctab.n_atoms != 0 {
        let count = ctab.n_atoms as u64;
        let orig = calloc_or_null!(i32, count, 4);
        heap.slice_mut(ctab.v3000)?[0].atom_index_orig = orig;
        let fin = calloc_or_null!(i32, count, 4);
        heap.slice_mut(ctab.v3000)?[0].atom_index_fin = fin;
        if !orig.is_null() && !fin.is_null() {
            let mut i = 0_i32;
            while i < ctab.n_atoms {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                heap.slice_mut(orig)?[index] = -1;
                heap.slice_mut(fin)?[index] = -1;
                i = i.wrapping_add(1);
            }
        }
    } else {
        heap.slice_mut(ctab.v3000)?[0].atom_index_orig = SourceMutPointer::null();
        heap.slice_mut(ctab.v3000)?[0].atom_index_fin = SourceMutPointer::null();
    }

    let stages = [
        (0_u8, 8_i32),
        (1_u8, 1_i32),
        (2_u8, 4_i32),
        (3_u8, 4_i32),
    ];
    let mut ret = 0_i32;
    for (stage, capacity) in stages {
        {
            let v3000 = &mut heap.slice_mut(ctab.v3000)?[0];
            match stage {
                0 => v3000.n_haptic_bonds = 0,
                1 => v3000.n_steabs = 0,
                2 => v3000.n_sterel = 0,
                _ => v3000.n_sterac = 0,
            }
        }
        let descriptor = calloc_or_null!(NUM_LISTS, 1, 24);
        {
            let v3000 = &mut heap.slice_mut(ctab.v3000)?[0];
            match stage {
                0 => v3000.haptic_bonds = descriptor,
                1 => v3000.steabs = descriptor,
                2 => v3000.sterel = descriptor,
                _ => v3000.sterac = descriptor,
            }
        }
        if descriptor.is_null() {
            let _ = mol_fmt3_add_message(p_str_err.as_deref_mut(), b"Out of RAM")?;
            return Ok(-1);
        }
        let mut value = heap.slice(descriptor.as_const())?[0].clone();
        ret = NumLists_Alloc(heap, Some(&mut value), capacity)?;
        heap.slice_mut(descriptor)?[0] = value;
        if ret < 0 {
            let v3000 = &mut heap.slice_mut(ctab.v3000)?[0];
            match stage {
                0 => v3000.haptic_bonds = SourceMutPointer::null(),
                1 => v3000.steabs = SourceMutPointer::null(),
                2 => v3000.sterel = SourceMutPointer::null(),
                _ => v3000.sterac = SourceMutPointer::null(),
            }
            let _ = mol_fmt3_add_message(p_str_err.as_deref_mut(), b"Out of RAM")?;
            return Ok(-1);
        }
    }
    Ok(ret)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn DeleteMolfileV3000Info(
    heap: &mut SourceHeap,
    mut v3000: SourceMutPointer<MOL_FMT_v3000>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:165 DeleteMolfileV3000Info
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap object access adds overhead.
    /*
int DeleteMolfileV3000Info(MOL_FMT_v3000 *v3000)
{
    if (v3000)
    {

        if (v3000->atom_index_orig)
        {
            inchi_free(v3000->atom_index_orig);
        }

        if (v3000->atom_index_fin)
        {
            inchi_free(v3000->atom_index_fin);
        }

        if (v3000->haptic_bonds)
        {
            NumLists_Free(v3000->haptic_bonds);
            free(v3000->haptic_bonds);
        }

        if (v3000->steabs)
        {
            NumLists_Free(v3000->steabs);
            free(v3000->steabs);
        }

        if (v3000->sterel)
        {
            NumLists_Free(v3000->sterel);
            free(v3000->sterel);
        }

        if (v3000->sterac)
        {
            NumLists_Free(v3000->sterac);
            free(v3000->sterac);
        }

        inchi_free(v3000);
        v3000 = NULL;
    }

    return 0;
}
    */
    // END INCHI C FUNCTION: DeleteMolfileV3000Info
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: DeleteMolfileV3000Info
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; inchi_free and free both resolve to libc free behavior.
    // INCHI✔️❌: The final local NULL assignment does not modify the caller because the source pointer parameter is passed by value.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: DeleteMolfileV3000Info

    if v3000.is_null() {
        return Ok(0);
    }
    let value = heap
        .slice(v3000.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if !value.atom_index_orig.is_null() {
        inchi_free(heap, value.atom_index_orig)?;
    }
    if !value.atom_index_fin.is_null() {
        inchi_free(heap, value.atom_index_fin)?;
    }
    for pointer in [value.haptic_bonds, value.steabs, value.sterel, value.sterac] {
        if !pointer.is_null() {
            let mut lists = heap
                .slice(pointer.as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            NumLists_Free(heap, Some(&mut lists))?;
            inchi_free(heap, pointer)?;
        }
    }
    inchi_free(heap, v3000)?;
    v3000 = SourceMutPointer::null();
    let _ = v3000;
    Ok(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source::base::ichi_io::inchi_strbuf_init;
    use crate::source_types::{INCHI_IOS_TYPE_STRING, MOL_COORD, MOL_FMT_ATOM, MOL_FMT_BOND};

    fn input(heap: &mut SourceHeap, value: &str) -> SourceMutPointer<i8> {
        heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain(std::iter::once(0)).collect())
            .unwrap()
    }

    fn stream(heap: &mut SourceHeap, value: &str) -> INCHI_IOSTREAM {
        let data = input(heap, value);
        INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            s: INCHI_IOS_STRING {
                pStr: data,
                nAllocatedLength: value.len() as i32 + 1,
                nUsedLength: value.len() as i32,
                nPtr: 0,
            },
            ..INCHI_IOSTREAM::default()
        }
    }

    #[test]
    fn source_port__mol_fmt3__get_v3000_input_line_to_strbuf__line_1889() {
        fn output(heap: &SourceHeap, buffer: &INCHI_IOS_STRING) -> Vec<u8> {
            heap.slice(buffer.pStr.as_const()).unwrap()[..buffer.nUsedLength as usize]
                .iter()
                .map(|byte| *byte as u8)
                .collect()
        }

        let mut heap = SourceHeap::default();
        let mut buffer = INCHI_IOS_STRING::default();
        inchi_strbuf_init(&mut heap, &mut buffer, 8, 8).unwrap();
        heap.slice_mut(buffer.pStr).unwrap()[..4].copy_from_slice(&[b'o' as i8, b'l' as i8, b'd' as i8, 0]);
        buffer.nUsedLength = 3;

        let mut single = stream(&mut heap, "M  V30 BEGIN CTAB\nrest\n");
        assert_eq!(
            get_V3000_input_line_to_strbuf(&mut heap, &mut buffer, Some(&mut single)),
            Ok(10)
        );
        assert_eq!(output(&heap, &buffer), b"BEGIN CTAB");
        assert_eq!(single.s.nPtr, "M  V30 BEGIN CTAB\n".len() as i32);

        let mut continued = stream(&mut heap, "M  V30 ABC-\r\nM  V30 DEF   \r\nfollowing\n");
        assert_eq!(
            get_V3000_input_line_to_strbuf(&mut heap, &mut buffer, Some(&mut continued)),
            Ok(6)
        );
        assert_eq!(output(&heap, &buffer), b"ABCDEF");
        assert_eq!(continued.s.nPtr, "M  V30 ABC-\r\nM  V30 DEF   \r\n".len() as i32);

        for malformed in ["", "M  V30 \n", "X  V30 VALUE\n"] {
            let mut malformed_stream = stream(&mut heap, malformed);
            assert_eq!(
                get_V3000_input_line_to_strbuf(&mut heap, &mut buffer, Some(&mut malformed_stream)),
                Ok(-1),
                "{malformed:?}"
            );
        }
        assert_eq!(
            get_V3000_input_line_to_strbuf(&mut heap, &mut buffer, None),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__mol_fmt3__molfilev3000readkeyword__line_410() {
        let mut heap = SourceHeap::default();
        let key = heap.allocate_model_storage(vec![0x55_i8; 16]).unwrap();

        let mut cursor = input(&mut heap, "ATOMS=VALUE");
        assert_eq!(MolfileV3000ReadKeyword(&mut heap, key, &mut cursor), Ok(5));
        assert_eq!(
            &heap.slice(key.as_const()).unwrap()[..7],
            &[b'A' as i8, b'T' as i8, b'O' as i8, b'M' as i8, b'S' as i8, 0, 0x55]
        );
        assert_eq!(heap.slice(cursor.as_const()).unwrap()[0], b'V' as i8);

        heap.slice_mut(key).unwrap().fill(0x44);
        cursor = input(&mut heap, "KEY VALUE");
        assert_eq!(MolfileV3000ReadKeyword(&mut heap, key, &mut cursor), Ok(3));
        assert_eq!(
            &heap.slice(key.as_const()).unwrap()[..4],
            &[b'K' as i8, b'E' as i8, b'Y' as i8, 0]
        );
        assert_eq!(heap.slice(cursor.as_const()).unwrap()[0], b' ' as i8);

        heap.slice_mut(key).unwrap().fill(0x33);
        cursor = input(&mut heap, "=VALUE");
        let equals = cursor;
        assert_eq!(MolfileV3000ReadKeyword(&mut heap, key, &mut cursor), Ok(0));
        assert_eq!(heap.slice(key.as_const()).unwrap()[..2], [0, 0x33]);
        assert_eq!(cursor, equals);

        cursor = input(&mut heap, "");
        heap.slice_mut(key).unwrap()[0] = 0x22;
        assert_eq!(MolfileV3000ReadKeyword(&mut heap, key, &mut cursor), Ok(0));
        assert_eq!(heap.slice(key.as_const()).unwrap()[0], 0);
        assert!(cursor.is_null());

        cursor = SourceMutPointer::null();
        heap.slice_mut(key).unwrap()[0] = 0x11;
        assert_eq!(MolfileV3000ReadKeyword(&mut heap, key, &mut cursor), Ok(-1));
        assert_eq!(heap.slice(key.as_const()).unwrap()[0], 0x11);
        assert!(cursor.is_null());
    }

    #[test]
    fn source_port__mol_fmt3__molfilev3000readctabbeginandcountsline__line_444() {
        fn text(buffer: &[i8]) -> String {
            let end = buffer.iter().position(|byte| *byte == 0).unwrap();
            String::from_utf8(buffer[..end].iter().map(|byte| *byte as u8).collect()).unwrap()
        }

        fn ctab(heap: &mut SourceHeap) -> MOL_FMT_CTAB {
            let v3000 = heap
                .allocate_model_storage(vec![MOL_FMT_v3000 {
                    n_sgroups: 77,
                    n_3d_constraints: 88,
                    ..MOL_FMT_v3000::default()
                }])
                .unwrap();
            MOL_FMT_CTAB {
                n_atoms: 91,
                n_bonds: 92,
                chiral_flag: 3,
                n_stext_entries: 94,
                n_property_lines: 95,
                v3000,
                ..MOL_FMT_CTAB::default()
            }
        }

        let mut heap = SourceHeap::default();
        let mut valid_ctab = ctab(&mut heap);
        let consumed = "M  V30 BEGIN CTAB\r\nM  V30 COUNTS 12 13 2 3 1\r\n";
        let mut valid = stream(&mut heap, &format!("{consumed}M  V30 NEXT\n"));
        let mut errors = [0_i8; 256];
        assert_eq!(
            MolfileV3000ReadCTABBeginAndCountsLine(&mut heap, &mut valid_ctab, Some(&mut valid), Some(&mut errors),),
            Ok(0)
        );
        assert_eq!(valid.s.nPtr, consumed.len() as i32);
        assert_eq!(
            (
                valid_ctab.n_atoms,
                valid_ctab.n_bonds,
                valid_ctab.chiral_flag,
                valid_ctab.n_stext_entries,
                valid_ctab.n_property_lines,
            ),
            (12, 13, 1, -1, 999)
        );
        let valid_v3000 = &heap.slice(valid_ctab.v3000.as_const()).unwrap()[0];
        assert_eq!((valid_v3000.n_sgroups, valid_v3000.n_3d_constraints), (2, 3));
        assert_eq!(errors[0], 0);

        let mut wrong_ctab = ctab(&mut heap);
        let wrong_before = wrong_ctab.clone();
        let mut wrong = stream(&mut heap, "M  V30 COUNTS 1 2 3 4 0\n");
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadCTABBeginAndCountsLine(&mut heap, &mut wrong_ctab, Some(&mut wrong), Some(&mut errors),),
            Ok(1)
        );
        assert_eq!(wrong_ctab, wrong_before);
        assert_eq!(text(&errors), "Error: No V3000 CTab start marker");

        let mut eof_ctab = ctab(&mut heap);
        let mut eof = stream(&mut heap, "M  V30 BEGIN CTAB\n");
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadCTABBeginAndCountsLine(&mut heap, &mut eof_ctab, Some(&mut eof), Some(&mut errors),),
            Ok(1)
        );
        assert_eq!(
            (
                eof_ctab.n_atoms,
                eof_ctab.n_bonds,
                eof_ctab.chiral_flag,
                eof_ctab.n_stext_entries,
                eof_ctab.n_property_lines,
            ),
            (-1, -1, -1, -1, 999)
        );
        let eof_v3000 = &heap.slice(eof_ctab.v3000.as_const()).unwrap()[0];
        assert_eq!((eof_v3000.n_sgroups, eof_v3000.n_3d_constraints), (77, 88));
        assert_eq!(text(&errors), "Cannot read V3000 counts line");

        let mut keyword_ctab = ctab(&mut heap);
        let mut keyword = stream(&mut heap, "M  V30 BEGIN CTAB\nM  V30 COUNT 1 2 3 4 0\n");
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadCTABBeginAndCountsLine(&mut heap, &mut keyword_ctab, Some(&mut keyword), Some(&mut errors),),
            Ok(1)
        );
        assert_eq!(keyword_ctab.n_atoms, -1);
        assert_eq!(text(&errors), "Cannot read V3000 counts line");

        let mut atom_overflow_ctab = ctab(&mut heap);
        let overflow_line = "COUNTS 2147483648 2 3 4 1\x01";
        let mut atom_overflow = stream(&mut heap, &format!("M  V30 BEGIN CTAB\nM  V30 {overflow_line}\n"));
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadCTABBeginAndCountsLine(
                &mut heap,
                &mut atom_overflow_ctab,
                Some(&mut atom_overflow),
                Some(&mut errors),
            ),
            Ok(3)
        );
        assert_eq!(atom_overflow_ctab.n_atoms, 0);
        assert_eq!(atom_overflow_ctab.n_bonds, -1);
        assert_eq!(
            text(&errors),
            "Number of atoms too large. V3000 counts line: COUNTS 2147483648 2 3 4 1."
        );

        let mut later_overflow_ctab = ctab(&mut heap);
        let mut later_overflow = stream(&mut heap, "M  V30 BEGIN CTAB\nM  V30 COUNTS 5 6 7 8 128\n");
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadCTABBeginAndCountsLine(
                &mut heap,
                &mut later_overflow_ctab,
                Some(&mut later_overflow),
                Some(&mut errors),
            ),
            Ok(3)
        );
        assert_eq!(
            (
                later_overflow_ctab.n_atoms,
                later_overflow_ctab.n_bonds,
                later_overflow_ctab.chiral_flag,
            ),
            (5, 6, 0)
        );
        let later_v3000 = &heap.slice(later_overflow_ctab.v3000.as_const()).unwrap()[0];
        assert_eq!((later_v3000.n_sgroups, later_v3000.n_3d_constraints), (7, 8));
        assert_eq!(text(&errors), "Cannot interpret V3000 counts line: COUNTS 5 6 7 8 128");

        let mut nonnumeric_ctab = ctab(&mut heap);
        let mut nonnumeric = stream(&mut heap, "M  V30 BEGIN CTAB\nM  V30 COUNTS x y z q r\n");
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadCTABBeginAndCountsLine(
                &mut heap,
                &mut nonnumeric_ctab,
                Some(&mut nonnumeric),
                Some(&mut errors),
            ),
            Ok(0)
        );
        assert_eq!(
            (
                nonnumeric_ctab.n_atoms,
                nonnumeric_ctab.n_bonds,
                nonnumeric_ctab.chiral_flag,
            ),
            (0, 0, 0)
        );
        let nonnumeric_v3000 = &heap.slice(nonnumeric_ctab.v3000.as_const()).unwrap()[0];
        assert_eq!((nonnumeric_v3000.n_sgroups, nonnumeric_v3000.n_3d_constraints), (0, 0));
        assert_eq!(errors[0], 0);
    }

    #[test]
    fn source_port__mol_fmt3__molfilev3000readhapticbond__line_1732() {
        let mut heap = SourceHeap::default();
        let mut output = heap.allocate_model_storage(vec![77_i32]).unwrap();
        let stale = output;
        let mut cursor = input(&mut heap, "(3 4 5 6) ATTACH=ALL");
        assert_eq!(
            MolfileV3000ReadHapticBond(&mut heap, None, &mut cursor, &mut output, None),
            Ok(10)
        );
        assert_ne!(output, stale);
        assert_eq!(heap.slice(output.as_const()).unwrap(), &[-1, -1, 3, 4, 5, 6]);
        assert!(cursor.is_null());
        inchi_free(&mut heap, output).unwrap();

        cursor = input(&mut heap, "(1 9) ");
        output = stale;
        assert_eq!(
            MolfileV3000ReadHapticBond(&mut heap, None, &mut cursor, &mut output, None),
            Ok(0)
        );
        assert_eq!(heap.slice(output.as_const()).unwrap(), &[-1, -1, 1, 9]);
        assert!(cursor.is_null());
        inchi_free(&mut heap, output).unwrap();

        for malformed in ["X(1 9) ", "( x ", "( -1 ", "(1 999999999999) ", "(1 9) ATTACH=ANY"] {
            cursor = input(&mut heap, malformed);
            output = stale;
            assert_eq!(
                MolfileV3000ReadHapticBond(&mut heap, None, &mut cursor, &mut output, None),
                Ok(-1),
                "{malformed}"
            );
            assert!(output.is_null(), "{malformed}");
        }

        cursor = input(&mut heap, "(1 9) ");
        output = stale;
        heap.fail_after_allocations(0);
        assert_eq!(
            MolfileV3000ReadHapticBond(&mut heap, None, &mut cursor, &mut output, None),
            Ok(-1)
        );
        assert!(output.is_null());
    }

    #[test]
    fn source_port__mol_fmt3__molfilev3000readstereocollection__line_1817() {
        let mut heap = SourceHeap::default();
        let mut output = heap.allocate_model_storage(vec![77_i32]).unwrap();
        let stale = output;
        let mut cursor = input(&mut heap, "(3 10 20 30) trailing");
        assert_eq!(
            MolfileV3000ReadStereoCollection(&mut heap, None, &mut cursor, &mut output, None,),
            Ok(1)
        );
        assert_ne!(output, stale);
        assert_eq!(heap.slice(output.as_const()).unwrap(), &[-1, 3, 10, 20, 30, 0]);
        assert_eq!(heap.slice(cursor.as_const()).unwrap()[0], b' ' as i8);
        inchi_free(&mut heap, output).unwrap();

        cursor = input(&mut heap, "(0) ");
        output = stale;
        assert_eq!(
            MolfileV3000ReadStereoCollection(&mut heap, None, &mut cursor, &mut output, None,),
            Ok(2)
        );
        assert_eq!(heap.slice(output.as_const()).unwrap(), &[-1, 0, 0]);
        inchi_free(&mut heap, output).unwrap();

        for malformed in ["X(1 9)", "( x", "( -1", "(2 1 999999999999)"] {
            cursor = input(&mut heap, malformed);
            output = stale;
            assert_eq!(
                MolfileV3000ReadStereoCollection(&mut heap, None, &mut cursor, &mut output, None,),
                Ok(-1),
                "{malformed}"
            );
            assert!(output.is_null(), "{malformed}");
        }

        cursor = input(&mut heap, "(2 1 bad)");
        output = stale;
        assert_eq!(
            MolfileV3000ReadStereoCollection(&mut heap, None, &mut cursor, &mut output, None,),
            Ok(1)
        );
        assert_eq!(heap.slice(output.as_const()).unwrap(), &[-1, 2, 1, 0, 0]);
        inchi_free(&mut heap, output).unwrap();

        cursor = input(&mut heap, "(2 1)");
        output = stale;
        assert_eq!(
            MolfileV3000ReadStereoCollection(&mut heap, None, &mut cursor, &mut output, None,),
            Ok(1)
        );
        assert_eq!(heap.slice(output.as_const()).unwrap(), &[-1, 2, 1, 0, 0]);
        inchi_free(&mut heap, output).unwrap();

        cursor = input(&mut heap, "(1 9)");
        output = stale;
        heap.fail_after_allocations(0);
        assert_eq!(
            MolfileV3000ReadStereoCollection(&mut heap, None, &mut cursor, &mut output, None,),
            Ok(-1)
        );
        assert!(output.is_null());
    }

    #[test]
    fn source_port__mol_fmt3__molfilev3000readbondsblock__line_1262() {
        fn error_text(buffer: &[i8]) -> String {
            let end = buffer.iter().position(|byte| *byte == 0).unwrap();
            String::from_utf8(buffer[..end].iter().map(|byte| *byte as u8).collect()).unwrap()
        }

        fn ctab(heap: &mut SourceHeap, bond_count: i32) -> MOL_FMT_CTAB {
            let orig = heap.allocate_model_storage(vec![1_i32, 2, 3]).unwrap();
            let fin = heap.allocate_model_storage(vec![1_i32, 2, -1]).unwrap();
            let list_array = heap.allocate_model_storage(vec![SourceMutPointer::null(); 2]).unwrap();
            let haptic = heap
                .allocate_model_storage(vec![NUM_LISTS {
                    lists: list_array,
                    allocated: 2,
                    used: 0,
                    increment: 2,
                }])
                .unwrap();
            let v3000 = heap
                .allocate_model_storage(vec![MOL_FMT_v3000 {
                    n_non_star_atoms: 2,
                    n_star_atoms: 1,
                    atom_index_orig: orig,
                    atom_index_fin: fin,
                    haptic_bonds: haptic,
                    ..MOL_FMT_v3000::default()
                }])
                .unwrap();
            let bonds = heap
                .allocate_model_storage(vec![MOL_FMT_BOND::default(); bond_count.max(0) as usize])
                .unwrap();
            MOL_FMT_CTAB {
                n_bonds: bond_count,
                bonds,
                v3000,
                ..MOL_FMT_CTAB::default()
            }
        }

        let mut heap = SourceHeap::default();
        let mut empty = MOL_FMT_CTAB::default();
        assert_eq!(MolfileV3000ReadBondsBlock(&mut heap, &mut empty, None, 17, None), Ok(0));

        let mut structure = ctab(&mut heap, 3);
        let mut input_stream = stream(
            &mut heap,
            "M  V30 BEGIN BOND\n\
             M  V30 1 1 1 2 CFG=2 TOPO=1 RXCTR=0 STBOX=0 DISP=X ATTACH=ALL\n\
             M  V30 2 2 1 3 ENDPTS=(2 2 3) ATTACH=ALL\n\
             M  V30 3 3 2 1 CFG=3\n\
             M  V30 END BOND\n",
        );
        let mut errors = [0_i8; 256];
        assert_eq!(
            MolfileV3000ReadBondsBlock(&mut heap, &mut structure, Some(&mut input_stream), 0, Some(&mut errors),),
            Ok(0)
        );
        assert_eq!(structure.n_bonds, 2);
        let v3000 = heap.slice(structure.v3000.as_const()).unwrap()[0].clone();
        assert_eq!((v3000.n_non_haptic_bonds, v3000.n_haptic_bonds), (2, 1));
        assert_eq!(
            &heap.slice(structure.bonds.as_const()).unwrap()[..2],
            &[
                MOL_FMT_BOND {
                    atnum1: 1,
                    atnum2: 2,
                    bond_type: 1,
                    bond_stereo: 4,
                },
                MOL_FMT_BOND {
                    atnum1: 2,
                    atnum2: 1,
                    bond_type: 3,
                    bond_stereo: 6,
                },
            ]
        );
        let lists = heap.slice(v3000.haptic_bonds.as_const()).unwrap()[0].clone();
        assert_eq!(lists.used, 1);
        let haptic = heap.slice(lists.lists.as_const()).unwrap()[0];
        assert_eq!(heap.slice(haptic.as_const()).unwrap(), &[2, 1, 2, 2, 3]);
        assert_eq!(error_text(&errors), "V3000 haptic bonds read/stored but ignored");

        let mut bad_start = ctab(&mut heap, 1);
        let mut bad_start_stream = stream(&mut heap, "M  V30 BEGIN ATOM\n");
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadBondsBlock(
                &mut heap,
                &mut bad_start,
                Some(&mut bad_start_stream),
                0,
                Some(&mut errors),
            ),
            Ok(1)
        );
        assert_eq!(error_text(&errors), "Error: No V3000 Bond block start marker");

        let mut preexisting = ctab(&mut heap, 1);
        let mut end_data = stream(&mut heap, "M  V30 BEGIN BOND\nM  V30 $$$$\nM  V30 END BOND\n");
        assert_eq!(
            MolfileV3000ReadBondsBlock(&mut heap, &mut preexisting, Some(&mut end_data), 7, None,),
            Ok(-7)
        );

        let mut malformed = ctab(&mut heap, 1);
        let mut malformed_stream = stream(&mut heap, "M  V30 BEGIN BOND\nM  V30 bad\nM  V30 END BOND\n");
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadBondsBlock(
                &mut heap,
                &mut malformed,
                Some(&mut malformed_stream),
                0,
                Some(&mut errors),
            ),
            Ok(4)
        );
        let malformed_v3000 = &heap.slice(malformed.v3000.as_const()).unwrap()[0];
        assert_eq!(malformed_v3000.n_non_haptic_bonds, 1);
        assert_eq!(
            heap.slice(malformed.bonds.as_const()).unwrap()[0],
            MOL_FMT_BOND {
                atnum1: -1,
                atnum2: -1,
                bond_type: 0,
                bond_stereo: 0,
            }
        );
        assert_eq!(error_text(&errors), "Cannot interpret V3000 bond block line: bad");
    }

    #[test]
    fn source_port__mol_fmt3__molfilev3000readsgroup__line_561() {
        let mut heap = SourceHeap::default();
        let mut ctab = MOL_FMT_CTAB {
            n_atoms: 7,
            n_bonds: 3,
            ..MOL_FMT_CTAB::default()
        };
        let original = ctab.clone();
        let consumed = "M  V30 1 SUP 0 ATOMS=(1 1)\r\n\
                        M  V30 END SGROUP X\n\
                        M  V30 END SGROUP\n";
        let mut input_stream = stream(&mut heap, &format!("{consumed}M  V30 FOLLOWING\n"));
        let mut errors = [0x55_i8; 256];
        assert_eq!(
            MolfileV3000ReadSGroup(
                &mut heap,
                Some(&mut ctab),
                Some(&mut input_stream),
                73,
                Some(&mut errors),
            ),
            Ok(0)
        );
        assert_eq!(ctab, original);
        assert_eq!(errors, [0x55_i8; 256]);
        assert_eq!(input_stream.s.nPtr, consumed.len() as i32);
        assert_eq!(
            &heap.slice(input_stream.s.pStr.as_const()).unwrap()
                [input_stream.s.nPtr as usize..input_stream.s.nPtr as usize + 7],
            &b"M  V30 ".map(|byte| byte as i8)
        );
    }

    #[test]
    fn source_port__mol_fmt3__molfilev3000read3dblock__line_608() {
        fn text(buffer: &[i8]) -> String {
            let end = buffer.iter().position(|byte| *byte == 0).unwrap();
            String::from_utf8(buffer[..end].iter().map(|byte| *byte as u8).collect()).unwrap()
        }

        let mut heap = SourceHeap::default();
        let mut ctab = MOL_FMT_CTAB {
            n_atoms: 4,
            ..MOL_FMT_CTAB::default()
        };
        let original = ctab.clone();
        let mut valid = stream(&mut heap, "M  V30 END OBJ3D\r\nM  V30 NEXT\n");
        let mut errors = [0_i8; 256];
        assert_eq!(
            MolfileV3000Read3DBlock(&mut heap, Some(&mut ctab), Some(&mut valid), 9, Some(&mut errors),),
            Ok(9)
        );
        assert_eq!(ctab, original);
        assert_eq!(errors[0], 0);
        assert_eq!(valid.s.nPtr, "M  V30 END OBJ3D\r\n".len() as i32);

        let mut wrong = stream(&mut heap, "M  V30 BEGIN OBJ3D\n");
        assert_eq!(
            MolfileV3000Read3DBlock(&mut heap, None, Some(&mut wrong), 0, Some(&mut errors),),
            Ok(1)
        );
        assert_eq!(text(&errors), "Error: No V3000 3DBlock end marker");

        errors.fill(0);
        let mut wrong_existing = stream(&mut heap, "M  V30 WRONG\n");
        assert_eq!(
            MolfileV3000Read3DBlock(&mut heap, None, Some(&mut wrong_existing), 27, Some(&mut errors),),
            Ok(27)
        );
        assert_eq!(text(&errors), "Error: No V3000 3DBlock end marker");

        errors.fill(0);
        let mut eof = stream(&mut heap, "");
        assert_eq!(
            MolfileV3000Read3DBlock(&mut heap, None, Some(&mut eof), 0, Some(&mut errors),),
            Ok(1)
        );
        assert_eq!(text(&errors), "Error: No V3000 3DBlock end marker");
    }

    #[test]
    fn source_port__mol_fmt3__molfilev3000readcollections__line_647() {
        fn text(buffer: &[i8]) -> String {
            let end = buffer.iter().position(|byte| *byte == 0).unwrap();
            String::from_utf8(buffer[..end].iter().map(|byte| *byte as u8).collect()).unwrap()
        }

        fn descriptor(heap: &mut SourceHeap) -> SourceMutPointer<NUM_LISTS> {
            let array = heap.allocate_model_storage(vec![SourceMutPointer::null(); 2]).unwrap();
            heap.allocate_model_storage(vec![NUM_LISTS {
                lists: array,
                allocated: 2,
                used: 0,
                increment: 2,
            }])
            .unwrap()
        }

        fn collection_ctab(heap: &mut SourceHeap) -> MOL_FMT_CTAB {
            let orig = heap.allocate_model_storage(vec![10_i32, 20, 30]).unwrap();
            let fin = heap.allocate_model_storage(vec![1_i32, 2, 3]).unwrap();
            let steabs = descriptor(heap);
            let sterel = descriptor(heap);
            let sterac = descriptor(heap);
            let v3000 = heap
                .allocate_model_storage(vec![MOL_FMT_v3000 {
                    n_non_star_atoms: 3,
                    atom_index_orig: orig,
                    atom_index_fin: fin,
                    steabs,
                    sterel,
                    sterac,
                    ..MOL_FMT_v3000::default()
                }])
                .unwrap();
            MOL_FMT_CTAB {
                v3000,
                ..MOL_FMT_CTAB::default()
            }
        }

        fn only_list(heap: &SourceHeap, descriptor: SourceMutPointer<NUM_LISTS>) -> Vec<i32> {
            let lists = heap.slice(descriptor.as_const()).unwrap()[0].clone();
            let pointer = heap.slice(lists.lists.as_const()).unwrap()[0];
            heap.slice(pointer.as_const()).unwrap().to_vec()
        }

        let mut heap = SourceHeap::default();
        let mut ctab = collection_ctab(&mut heap);
        let consumed = "M  V30 MDLV30/STEABS ATOMS=(3 10 20 30)\n\
                        M  V30 MDLV30/STEREL2 ATOMS=(2 20 30)\n\
                        M  V30 MDLV30/STERAC3 ATOMS=(1 30)\n\
                        M  V30 MDLV30/UNKNOWN DATA=ignored\n\
                        M  V30 END COLLECTION\n";
        let mut input_stream = stream(&mut heap, &format!("{consumed}M  V30 NEXT\n"));
        let mut errors = [0_i8; 256];
        assert_eq!(
            MolfileV3000ReadCollections(&mut heap, &mut ctab, Some(&mut input_stream), 0, Some(&mut errors),),
            Ok(0)
        );
        assert_eq!(input_stream.s.nPtr, consumed.len() as i32);
        let v3000 = heap.slice(ctab.v3000.as_const()).unwrap()[0].clone();
        assert_eq!(
            (v3000.n_collections, v3000.n_steabs, v3000.n_sterel, v3000.n_sterac,),
            (3, 1, 1, 1)
        );
        assert_eq!(only_list(&heap, v3000.steabs), vec![1, 3, 1, 20, 30, 0]);
        assert_eq!(only_list(&heap, v3000.sterel), vec![2, 2, 20, 30, 0]);
        assert_eq!(only_list(&heap, v3000.sterac), vec![3, 1, 30, 0]);
        assert_eq!(text(&errors), "V3000 enhanced stereo read/stored but ignored");

        let mut failed_ctab = collection_ctab(&mut heap);
        let failure_consumed = "M  V30 BAD\nM  V30 AFTER\x01\n";
        let mut failed_stream = stream(&mut heap, &format!("{failure_consumed}M  V30 UNREAD\n"));
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadCollections(
                &mut heap,
                &mut failed_ctab,
                Some(&mut failed_stream),
                99,
                Some(&mut errors),
            ),
            Ok(7)
        );
        assert_eq!(failed_stream.s.nPtr, failure_consumed.len() as i32);
        assert_eq!(text(&errors), "Cannot interpret V3000 collection line(s); AFTER.");
        let failed_v3000 = &heap.slice(failed_ctab.v3000.as_const()).unwrap()[0];
        assert_eq!(failed_v3000.n_collections, 0);
    }

    #[test]
    fn source_port__mol_fmt3__molfilev3000readtailofctab__line_1602() {
        fn text(buffer: &[i8]) -> String {
            let end = buffer.iter().position(|byte| *byte == 0).unwrap();
            String::from_utf8(buffer[..end].iter().map(|byte| *byte as u8).collect()).unwrap()
        }

        fn ctab(heap: &mut SourceHeap) -> MOL_FMT_CTAB {
            let v3000 = heap.allocate_model_storage(vec![MOL_FMT_v3000::default()]).unwrap();
            MOL_FMT_CTAB {
                v3000,
                ..MOL_FMT_CTAB::default()
            }
        }

        let mut heap = SourceHeap::default();
        let mut structure = ctab(&mut heap);
        let consumed = "M  V30 BEGIN SGROUP\n\
                        M  V30 ignored\n\
                        M  V30 END SGROUP\n\
                        M  V30 BEGIN OBJ3D\n\
                        M  V30 END OBJ3D\n\
                        M  V30 LINKNODE\n\
                        M  V30 LINKNODE\n\
                        M  V30 BEGIN COLLECTION\n\
                        M  V30 END COLLECTION\n\
                        M  V30 END CTAB\n";
        let mut full = stream(&mut heap, &format!("{consumed}M  V30 FOLLOWING\n"));
        let mut errors = [0_i8; 256];
        assert_eq!(
            MolfileV3000ReadTailOfCTAB(&mut heap, &mut structure, Some(&mut full), 0, Some(&mut errors),),
            Ok(0)
        );
        assert_eq!(full.s.nPtr, consumed.len() as i32);
        assert_eq!(errors[0], 0);

        let mut direct_structure = ctab(&mut heap);
        let mut direct = stream(&mut heap, "M  V30 END CTAB\r\nM  V30 NEXT\n");
        assert_eq!(
            MolfileV3000ReadTailOfCTAB(
                &mut heap,
                &mut direct_structure,
                Some(&mut direct),
                23,
                Some(&mut errors),
            ),
            Ok(23)
        );
        assert_eq!(direct.s.nPtr, "M  V30 END CTAB\r\n".len() as i32);

        let mut missing_structure = ctab(&mut heap);
        let mut missing = stream(&mut heap, "M  V30 WRONG\nM  V30 UNREAD\n");
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadTailOfCTAB(
                &mut heap,
                &mut missing_structure,
                Some(&mut missing),
                0,
                Some(&mut errors),
            ),
            Ok(1)
        );
        assert_eq!(text(&errors), "Error: No V3000 CTAB end marker");
        assert_eq!(missing.s.nPtr, "M  V30 WRONG\n".len() as i32);

        let mut eof_structure = ctab(&mut heap);
        let mut eof = stream(&mut heap, "");
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadTailOfCTAB(&mut heap, &mut eof_structure, Some(&mut eof), 31, Some(&mut errors),),
            Ok(31)
        );
        assert_eq!(text(&errors), "Error: No V3000 CTAB end marker");

        let mut child_structure = ctab(&mut heap);
        let mut child_error = stream(&mut heap, "M  V30 BEGIN OBJ3D\nM  V30 WRONG\nM  V30 END CTAB\n");
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadTailOfCTAB(
                &mut heap,
                &mut child_structure,
                Some(&mut child_error),
                0,
                Some(&mut errors),
            ),
            Ok(0)
        );
        assert_eq!(text(&errors), "Error: No V3000 3DBlock end marker");
        assert_eq!(child_error.s.nPtr, "M  V30 BEGIN OBJ3D\nM  V30 WRONG\n".len() as i32);

        let mut existing_structure = ctab(&mut heap);
        let mut existing = stream(&mut heap, "M  V30 BEGIN OBJ3D\nM  V30 END OBJ3D\nM  V30 END CTAB\n");
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadTailOfCTAB(
                &mut heap,
                &mut existing_structure,
                Some(&mut existing),
                9,
                Some(&mut errors),
            ),
            Ok(9)
        );
        assert_eq!(errors[0], 0);
        assert_eq!(existing.s.nPtr, "M  V30 BEGIN OBJ3D\nM  V30 END OBJ3D\n".len() as i32);
    }

    #[test]
    fn source_port__mol_fmt3__molfilev3000init__line_68() {
        fn text(buffer: &[i8]) -> String {
            let end = buffer.iter().position(|byte| *byte == 0).unwrap();
            String::from_utf8(buffer[..end].iter().map(|byte| *byte as u8).collect()).unwrap()
        }

        fn ctab(heap: &mut SourceHeap, atom_count: i32) -> MOL_FMT_CTAB {
            let v3000 = heap
                .allocate_model_storage(vec![MOL_FMT_v3000 {
                    n_star_atoms: 11,
                    n_non_star_atoms: 12,
                    n_haptic_bonds: 13,
                    n_steabs: 14,
                    n_sterel: 15,
                    n_sterac: 16,
                    ..MOL_FMT_v3000::default()
                }])
                .unwrap();
            MOL_FMT_CTAB {
                n_atoms: atom_count,
                v3000,
                ..MOL_FMT_CTAB::default()
            }
        }

        fn descriptor(heap: &SourceHeap, pointer: SourceMutPointer<NUM_LISTS>) -> NUM_LISTS {
            heap.slice(pointer.as_const()).unwrap()[0].clone()
        }

        let mut heap = SourceHeap::default();
        let old_index = heap.allocate_model_storage(vec![44_i32]).unwrap();
        let mut valid = ctab(&mut heap, 3);
        heap.slice_mut(valid.v3000).unwrap()[0].atom_index_orig = old_index;
        let mut errors = [0_i8; 256];
        assert_eq!(MolfileV3000Init(&mut heap, &mut valid, Some(&mut errors)), Ok(0));
        let initialized = heap.slice(valid.v3000.as_const()).unwrap()[0].clone();
        assert_eq!((initialized.n_star_atoms, initialized.n_non_star_atoms), (0, 0));
        assert_eq!(
            heap.slice(initialized.atom_index_orig.as_const()).unwrap(),
            &[-1, -1, -1]
        );
        assert_eq!(
            heap.slice(initialized.atom_index_fin.as_const()).unwrap(),
            &[-1, -1, -1]
        );
        assert_eq!(heap.slice(old_index.as_const()).unwrap(), &[44]);
        assert_eq!(
            (
                initialized.n_haptic_bonds,
                initialized.n_steabs,
                initialized.n_sterel,
                initialized.n_sterac,
            ),
            (0, 0, 0, 0)
        );
        for (pointer, capacity) in [
            (initialized.haptic_bonds, 8),
            (initialized.steabs, 1),
            (initialized.sterel, 4),
            (initialized.sterac, 4),
        ] {
            let lists = descriptor(&heap, pointer);
            assert_eq!((lists.allocated, lists.used, lists.increment), (capacity, 0, capacity));
            assert_eq!(heap.slice(lists.lists.as_const()).unwrap().len(), capacity as usize);
            assert!(
                heap.slice(lists.lists.as_const())
                    .unwrap()
                    .iter()
                    .all(|pointer| pointer.is_null())
            );
        }
        assert_eq!(errors[0], 0);

        let mut zero = ctab(&mut heap, 0);
        assert_eq!(MolfileV3000Init(&mut heap, &mut zero, None), Ok(0));
        let zero_v3000 = &heap.slice(zero.v3000.as_const()).unwrap()[0];
        assert!(zero_v3000.atom_index_orig.is_null());
        assert!(zero_v3000.atom_index_fin.is_null());

        let mut negative = ctab(&mut heap, -1);
        assert_eq!(MolfileV3000Init(&mut heap, &mut negative, None), Ok(0));
        let negative_v3000 = &heap.slice(negative.v3000.as_const()).unwrap()[0];
        assert!(negative_v3000.atom_index_orig.is_null());
        assert!(negative_v3000.atom_index_fin.is_null());

        let mut atom_failure_heap = SourceHeap::default();
        let mut atom_failure = ctab(&mut atom_failure_heap, 2);
        atom_failure_heap.fail_after_allocations(0);
        assert_eq!(MolfileV3000Init(&mut atom_failure_heap, &mut atom_failure, None), Ok(0));
        let atom_failure_v3000 = &atom_failure_heap.slice(atom_failure.v3000.as_const()).unwrap()[0];
        assert!(atom_failure_v3000.atom_index_orig.is_null());
        assert_eq!(
            atom_failure_heap
                .slice(atom_failure_v3000.atom_index_fin.as_const())
                .unwrap(),
            &[0, 0]
        );

        for failure_ordinal in 0_u64..8 {
            let mut failure_heap = SourceHeap::default();
            let mut failure = ctab(&mut failure_heap, 0);
            let mut failure_errors = [0_i8; 256];
            failure_heap.fail_after_allocations(failure_ordinal);
            assert_eq!(
                MolfileV3000Init(&mut failure_heap, &mut failure, Some(&mut failure_errors),),
                Ok(-1),
                "allocation ordinal {failure_ordinal}"
            );
            assert_eq!(text(&failure_errors), "Out of RAM");
            assert_eq!(
                failure_heap.source_allocation_calls(),
                failure_ordinal + 1,
                "allocation ordinal {failure_ordinal}"
            );
        }
    }

    #[test]
    fn source_port__mol_fmt3__molfilev3000readfield__line_257() {
        let mut heap = SourceHeap::default();

        let mut cursor = input(&mut heap, " C rest");
        let mut string = [0x55_i8; 8];
        assert_eq!(
            MolfileV3000ReadField(
                &mut heap,
                MolfileV3000FieldData::String(&mut string),
                i32::from(MOL_FMT_STRING_DATA),
                &mut cursor,
            ),
            Ok(1)
        );
        assert_eq!(&string[..3], &[b'C' as i8, 0, 0x55]);
        assert_eq!(heap.slice(cursor.as_const()).unwrap()[0], b' ' as i8);

        cursor = input(&mut heap, "ABCDEF ");
        string.fill(0x44);
        assert_eq!(
            MolfileV3000ReadField(
                &mut heap,
                MolfileV3000FieldData::String(&mut string),
                i32::from(MOL_FMT_STRING_DATA),
                &mut cursor,
            ),
            Ok(6)
        );
        assert_eq!(&string[..7], b"ABCDEF\0".map(|byte| byte as i8).as_slice());
        cursor = input(&mut heap, "ABCDEFG ");
        string.fill(0x33);
        assert_eq!(
            MolfileV3000ReadField(
                &mut heap,
                MolfileV3000FieldData::String(&mut string),
                i32::from(MOL_FMT_STRING_DATA),
                &mut cursor,
            ),
            Ok(7)
        );
        assert_eq!(&string[..3], &[0, 0x33, 0x33]);

        let mut char_value = 9_i8;
        cursor = input(&mut heap, "127 ");
        assert_eq!(
            MolfileV3000ReadField(
                &mut heap,
                MolfileV3000FieldData::Char(&mut char_value),
                i32::from(MOL_FMT_CHAR_INT_DATA),
                &mut cursor
            ),
            Ok(3)
        );
        assert_eq!(char_value, 127);
        cursor = input(&mut heap, "128 ");
        assert_eq!(
            MolfileV3000ReadField(
                &mut heap,
                MolfileV3000FieldData::Char(&mut char_value),
                i32::from(MOL_FMT_CHAR_INT_DATA),
                &mut cursor
            ),
            Ok(-1)
        );
        assert_eq!(char_value, 0);

        let mut short_value = 0_i16;
        cursor = input(&mut heap, "-32768 ");
        assert_eq!(
            MolfileV3000ReadField(
                &mut heap,
                MolfileV3000FieldData::Short(&mut short_value),
                i32::from(MOL_FMT_SHORT_INT_DATA),
                &mut cursor
            ),
            Ok(6)
        );
        assert_eq!(short_value, i16::MIN);

        let mut int_value = 0_i32;
        cursor = input(&mut heap, "2147483647 ");
        assert_eq!(
            MolfileV3000ReadField(
                &mut heap,
                MolfileV3000FieldData::Int(&mut int_value),
                i32::from(MOL_FMT_INT_DATA),
                &mut cursor
            ),
            Ok(10)
        );
        assert_eq!(int_value, i32::MAX);
        cursor = input(&mut heap, "12x ");
        assert_eq!(
            MolfileV3000ReadField(
                &mut heap,
                MolfileV3000FieldData::Int(&mut int_value),
                i32::from(MOL_FMT_INT_DATA),
                &mut cursor
            ),
            Ok(3)
        );
        assert_eq!(int_value, 12);

        let mut long_value = 7_i64;
        cursor = input(&mut heap, "9223372036854775807 ");
        assert_eq!(
            MolfileV3000ReadField(
                &mut heap,
                MolfileV3000FieldData::Long(&mut long_value),
                i32::from(MOL_FMT_LONG_INT_DATA),
                &mut cursor
            ),
            Ok(-1)
        );
        assert_eq!(long_value, 0);
        cursor = input(&mut heap, "9223372036854775806 ");
        assert_eq!(
            MolfileV3000ReadField(
                &mut heap,
                MolfileV3000FieldData::Long(&mut long_value),
                i32::from(MOL_FMT_LONG_INT_DATA),
                &mut cursor
            ),
            Ok(19)
        );
        assert_eq!(long_value, i64::MAX - 1);

        let mut double_value = 9.0_f64;
        cursor = input(&mut heap, "1.25 ");
        assert_eq!(
            MolfileV3000ReadField(
                &mut heap,
                MolfileV3000FieldData::Double(&mut double_value),
                i32::from(MOL_FMT_DOUBLE_DATA),
                &mut cursor
            ),
            Ok(4)
        );
        assert_eq!(double_value, 1.25);
        cursor = input(&mut heap, "1e9999 ");
        assert_eq!(
            MolfileV3000ReadField(
                &mut heap,
                MolfileV3000FieldData::Double(&mut double_value),
                i32::from(MOL_FMT_DOUBLE_DATA),
                &mut cursor
            ),
            Ok(-1)
        );
        assert_eq!(double_value, 0.0);

        let mut float_value = 7.25_f32;
        cursor = input(&mut heap, "1.5 ");
        assert_eq!(
            MolfileV3000ReadField(
                &mut heap,
                MolfileV3000FieldData::Float(&mut float_value),
                i32::from(MOL_FMT_FLOAT_DATA),
                &mut cursor
            ),
            Ok(3)
        );
        assert_eq!(float_value, 7.25);
        cursor = input(&mut heap, "1e-50 ");
        assert_eq!(
            MolfileV3000ReadField(
                &mut heap,
                MolfileV3000FieldData::Float(&mut float_value),
                i32::from(MOL_FMT_FLOAT_DATA),
                &mut cursor
            ),
            Ok(5)
        );
        assert_eq!(float_value, 0.0);
        cursor = input(&mut heap, "1e40 ");
        assert_eq!(
            MolfileV3000ReadField(
                &mut heap,
                MolfileV3000FieldData::Float(&mut float_value),
                i32::from(MOL_FMT_FLOAT_DATA),
                &mut cursor
            ),
            Ok(-1)
        );
        assert_eq!(float_value, 0.0);

        cursor = input(&mut heap, "");
        int_value = 17;
        assert_eq!(
            MolfileV3000ReadField(
                &mut heap,
                MolfileV3000FieldData::Int(&mut int_value),
                i32::from(MOL_FMT_INT_DATA),
                &mut cursor
            ),
            Ok(0)
        );
        assert_eq!(int_value, 0);
        assert!(cursor.is_null());
        cursor = input(&mut heap, "token ");
        assert_eq!(
            MolfileV3000ReadField(&mut heap, MolfileV3000FieldData::None, -77, &mut cursor),
            Ok(-1)
        );
    }

    #[test]
    fn source_port__mol_fmt3__get_actual_atom_number__line_1585() {
        let orig = [10, -7, 10, i32::MAX];
        let fin = [1, -1, 3, i32::MIN];
        assert_eq!(get_actual_atom_number(10, 4, Some(&orig), Some(&fin)), Ok(1));
        assert_eq!(get_actual_atom_number(-7, 4, Some(&orig), Some(&fin)), Ok(-1));
        assert_eq!(
            get_actual_atom_number(i32::MAX, 4, Some(&orig), Some(&fin)),
            Ok(i32::MIN)
        );
        assert_eq!(get_actual_atom_number(10, 1, Some(&orig), Some(&fin)), Ok(1));
        assert_eq!(get_actual_atom_number(10, 0, None, None), Ok(-1));
        assert_eq!(get_actual_atom_number(10, i32::MIN, None, None), Ok(-1));
        assert_eq!(
            get_actual_atom_number(10, 1, None, Some(&fin)),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(
            get_actual_atom_number(10, 5, Some(&orig), Some(&fin)),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            get_actual_atom_number(10, 4, Some(&orig), Some(&fin[..3])),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__mol_fmt3__deletemolfilev3000info__line_165() {
        let mut heap = SourceHeap::default();
        assert_eq!(DeleteMolfileV3000Info(&mut heap, SourceMutPointer::null()), Ok(0));

        let atom_index_orig = heap.allocate_model_storage(vec![3_i32, 1]).unwrap();
        let atom_index_fin = heap.allocate_model_storage(vec![1_i32, 3]).unwrap();
        let mut child_lists = Vec::new();
        let mut list_arrays = Vec::new();
        let mut descriptors = Vec::new();
        for value in 10_i32..14 {
            let child = heap.allocate_model_storage(vec![value]).unwrap();
            let array = heap.allocate_model_storage(vec![child]).unwrap();
            let descriptor = heap
                .allocate_model_storage(vec![NUM_LISTS {
                    lists: array,
                    allocated: 1,
                    used: 1,
                    increment: 1,
                }])
                .unwrap();
            child_lists.push(child);
            list_arrays.push(array);
            descriptors.push(descriptor);
        }
        let object = heap
            .allocate_model_storage(vec![MOL_FMT_v3000 {
                atom_index_orig,
                atom_index_fin,
                haptic_bonds: descriptors[0],
                steabs: descriptors[1],
                sterel: descriptors[2],
                sterac: descriptors[3],
                ..MOL_FMT_v3000::default()
            }])
            .unwrap();
        assert_eq!(DeleteMolfileV3000Info(&mut heap, object), Ok(0));
        for pointer in child_lists {
            assert_eq!(heap.slice(pointer.as_const()), Err(SourceHeapError::MissingAllocation));
        }
        for pointer in list_arrays {
            assert_eq!(heap.slice(pointer.as_const()), Err(SourceHeapError::MissingAllocation));
        }
        for pointer in descriptors {
            assert_eq!(heap.slice(pointer.as_const()), Err(SourceHeapError::MissingAllocation));
        }
        assert_eq!(
            heap.slice(atom_index_orig.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(atom_index_fin.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(heap.slice(object.as_const()), Err(SourceHeapError::MissingAllocation));

        let empty = heap.allocate_model_storage(vec![MOL_FMT_v3000::default()]).unwrap();
        assert_eq!(DeleteMolfileV3000Info(&mut heap, empty), Ok(0));
        assert_eq!(heap.slice(empty.as_const()), Err(SourceHeapError::MissingAllocation));
    }

    #[test]
    fn source_port__mol_fmt3__molfilev3000readatomsblock__line_862() {
        fn text(buffer: &[i8]) -> String {
            let end = buffer.iter().position(|byte| *byte == 0).unwrap();
            String::from_utf8(buffer[..end].iter().map(|byte| *byte as u8).collect()).unwrap()
        }

        fn ctab(heap: &mut SourceHeap, atom_count: i32, with_atoms: bool, with_coords: bool) -> MOL_FMT_CTAB {
            let count = atom_count.max(0) as usize;
            let atom_index_orig = heap.allocate_model_storage(vec![-9_i32; count]).unwrap();
            let atom_index_fin = heap.allocate_model_storage(vec![-8_i32; count]).unwrap();
            let v3000 = heap
                .allocate_model_storage(vec![MOL_FMT_v3000 {
                    atom_index_orig,
                    atom_index_fin,
                    ..MOL_FMT_v3000::default()
                }])
                .unwrap();
            let atoms = if with_atoms {
                heap.allocate_model_storage(vec![MOL_FMT_ATOM::default(); count])
                    .unwrap()
            } else {
                SourceMutPointer::null()
            };
            let coords = if with_coords {
                heap.allocate_model_storage(vec![[0x55_i8; 32] as MOL_COORD; count])
                    .unwrap()
            } else {
                SourceMutPointer::null()
            };
            MOL_FMT_CTAB {
                n_atoms: atom_count,
                atoms,
                coords,
                v3000,
                ..MOL_FMT_CTAB::default()
            }
        }

        let mut heap = SourceHeap::default();
        let mut structure = ctab(&mut heap, 3, true, true);
        let consumed = "M  V30 BEGIN ATOM\n\
                        M  V30 1 CL 1 -2.5 0 0 CHG=-2 RAD=3 CFG=2 MASS=35 VAL=-1 HCOUNT=x STBOX=x INVRET=x EXACHG=x SUBST=x UNSAT=x RBCNT=x ATTCHPT=x RGROUPS=X ATTCHORD=Y CLASS=Z SEQID=x\n\
                        M  V30 2 * 4 5 6 0\n\
                        M  V30 3 C 0.00001 100000 -0 7 MASS=13 VAL=258 UNKNOWN=9\n\
                        M  V30 END ATOM\n";
        let mut valid = stream(&mut heap, &format!("{consumed}M  V30 NEXT\n"));
        let mut errors = [0_i8; 1024];
        assert_eq!(
            MolfileV3000ReadAtomsBlock(&mut heap, &mut structure, Some(&mut valid), 0, Some(&mut errors),),
            Ok(0)
        );
        assert_eq!(valid.s.nPtr, consumed.len() as i32);
        assert_eq!(structure.n_atoms, 2);
        let v3000 = heap.slice(structure.v3000.as_const()).unwrap()[0].clone();
        assert_eq!((v3000.n_non_star_atoms, v3000.n_star_atoms), (2, 1));
        assert_eq!(heap.slice(v3000.atom_index_orig.as_const()).unwrap(), &[1, 2, 3]);
        assert_eq!(heap.slice(v3000.atom_index_fin.as_const()).unwrap(), &[1, -1, 2]);
        let atoms = heap.slice(structure.atoms.as_const()).unwrap();
        assert_eq!(&atoms[0].symbol, &[b'C' as i8, b'l' as i8, 0, 0, 0, 0]);
        assert_eq!((atoms[0].fx, atoms[0].fy, atoms[0].fz), (1.0, -2.5, 0.0));
        assert_eq!(
            (
                atoms[0].charge,
                atoms[0].radical,
                atoms[0].stereo_parity,
                atoms[0].mass_difference,
                atoms[0].valence,
            ),
            (-2, 3, 2, 127, 15)
        );
        assert_eq!(&atoms[1].symbol, &[b'C' as i8, 0, 0, 0, 0, 0]);
        assert_eq!((atoms[1].fx, atoms[1].fy, atoms[1].fz), (0.00001, 100000.0, -0.0));
        assert_eq!((atoms[1].mass_difference, atoms[1].valence), (1, 2));
        let coords = heap.slice(structure.coords.as_const()).unwrap();
        assert_eq!(
            &coords[0][..31],
            b"         1      -2.5         0\0".map(|byte| byte as i8).as_slice()
        );
        assert_eq!(coords[0][31], 0x55);
        assert_eq!(
            &coords[1][..31],
            b"         4         5         6\0".map(|byte| byte as i8).as_slice()
        );
        assert_eq!(
            &coords[2][..31],
            b"     1e-05    100000        -0\0".map(|byte| byte as i8).as_slice()
        );
        assert_eq!(text(&errors), "V3000 star atoms ignored");

        let mut no_atoms = ctab(&mut heap, 1, false, false);
        let mut ignored = stream(&mut heap, "M  V30 BEGIN ATOM\nM  V30 unparsed input\nM  V30 END ATOM\n");
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadAtomsBlock(&mut heap, &mut no_atoms, Some(&mut ignored), 0, Some(&mut errors),),
            Ok(0)
        );
        assert_eq!(errors[0], 0);

        let mut bad_start = ctab(&mut heap, 1, true, false);
        let mut wrong = stream(&mut heap, "M  V30 BEGIN BOND\n");
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadAtomsBlock(&mut heap, &mut bad_start, Some(&mut wrong), 0, Some(&mut errors),),
            Ok(1)
        );
        assert_eq!(text(&errors), "Error: No V3000 Atom block start marker");

        let mut positional = ctab(&mut heap, 1, true, false);
        let mut overflow = stream(
            &mut heap,
            "M  V30 BEGIN ATOM\nM  V30 2147483648 C 0 0 0 0\nM  V30 END ATOM\n",
        );
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadAtomsBlock(&mut heap, &mut positional, Some(&mut overflow), 0, Some(&mut errors),),
            Ok(4)
        );
        assert_eq!(
            text(&errors),
            "Cannot interpret V3000 atom block line: 2147483648 C 0 0 0 0"
        );
        assert_eq!(heap.slice(positional.v3000.as_const()).unwrap()[0].n_non_star_atoms, 0);

        let mut keyword = ctab(&mut heap, 1, true, false);
        let mut narrowing = stream(
            &mut heap,
            "M  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0 CHG=128 RAD=2\nM  V30 END ATOM\n",
        );
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadAtomsBlock(&mut heap, &mut keyword, Some(&mut narrowing), 0, Some(&mut errors),),
            Ok(4)
        );
        let keyword_atom = &heap.slice(keyword.atoms.as_const()).unwrap()[0];
        assert_eq!((keyword_atom.charge, keyword_atom.radical), (0, 2));
        assert_eq!(
            text(&errors),
            "Cannot interpret V3000 atom block key-value pair; 1 C 0 0 0 0 CHG=128 RAD=2"
        );

        let mut mass = ctab(&mut heap, 1, true, false);
        let mut mass_overflow = stream(
            &mut heap,
            "M  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0 MASS=32768\nM  V30 END ATOM\n",
        );
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadAtomsBlock(&mut heap, &mut mass, Some(&mut mass_overflow), 0, Some(&mut errors),),
            Ok(4)
        );
        assert_eq!(
            text(&errors),
            "Isotopic data not recognized: 1 C 0 0 0 0 MASS=32768; Cannot interpret V3000 atom block key-value pair"
        );

        let mut preexisting = ctab(&mut heap, 2, true, false);
        let mut bypass = stream(
            &mut heap,
            "M  V30 BEGIN ATOM\nM  V30 skipped\nM  V30 $$$$\nM  V30 END ATOM\n",
        );
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadAtomsBlock(&mut heap, &mut preexisting, Some(&mut bypass), 9, Some(&mut errors),),
            Ok(-9)
        );
        assert_eq!(errors[0], 0);

        let mut missing_line = ctab(&mut heap, 1, true, false);
        let mut eof = stream(&mut heap, "M  V30 BEGIN ATOM\n");
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadAtomsBlock(&mut heap, &mut missing_line, Some(&mut eof), 0, Some(&mut errors),),
            Ok(2)
        );
        assert_eq!(
            text(&errors),
            "Cannot read V3000 atom block line; Error: No V3000 Atom block end marker"
        );

        let mut missing_end = ctab(&mut heap, 1, true, false);
        let mut wrong_end = stream(&mut heap, "M  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 WRONG\n");
        errors.fill(0);
        assert_eq!(
            MolfileV3000ReadAtomsBlock(&mut heap, &mut missing_end, Some(&mut wrong_end), 0, Some(&mut errors),),
            Ok(1)
        );
        assert_eq!(text(&errors), "Error: No V3000 Atom block end marker");
    }
}
