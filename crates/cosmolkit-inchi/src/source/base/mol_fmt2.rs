use crate::source::base::ichi_io::{inchi_fgetsLf_bytes, source_decimal_i64};
use crate::source::base::ichierr::AddErrorMessage;
use crate::source::base::ichierr::add_error_message_bytes;
use crate::source::base::ichiparm::source_strtod_bytes_with_end;
use crate::source::base::ichiprt2::{inchi_strtod, inchi_strtol, inchi_strtol_bytes};
use crate::source::base::mol_fmt3::DeleteMolfileV3000Info;
use crate::source::base::mol_fmt4::MolFmtSgroups_Free;
use crate::source::base::util::{inchi_free, inchi_memicmp, lrtrim_bytes, mystrncpy_within};
use crate::source_types::{
    FILE, INCHI_IOS_TYPE_FILE, INCHI_IOSTREAM, MOL_FMT_CHAR_INT_DATA, MOL_FMT_DATA,
    MOL_FMT_DOUBLE_DATA, MOL_FMT_FLOAT_DATA, MOL_FMT_HEADER_BLOCK, MOL_FMT_INPLINELEN,
    MOL_FMT_JUMP_TO_RIGHT, MOL_FMT_LONG_INT_DATA, MOL_FMT_MAX_VALUE_LEN, MOL_FMT_SHORT_INT_DATA,
    MOL_FMT_STRING_DATA, SourceConstPointer, SourceHeap, SourceHeapError, SourceMutPointer,
};

pub(crate) enum MolfileFieldData<'a> {
    None,
    String(&'a mut [i8]),
    Char(&'a mut i8),
    Short(&'a mut i16),
    ShortIntoInt(&'a mut i32),
    Long(&'a mut i64),
    Double(&'a mut f64),
    Float(&'a mut f32),
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileStrnread(
    heap: &mut SourceHeap,
    dest: SourceMutPointer<i8>,
    source: SourceConstPointer<i8>,
    len: i32,
    first_space: &mut SourceMutPointer<i8>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:65 MolfileStrnread
    // INCHI✔️❌: complete source frame follows verbatim; checked source-heap access is materially slower than raw pointers.
    /*
int MolfileStrnread(char *dest,
                    char *source,
                    int len,
                    char **first_space)
{
    /* required len >= 0; dest must have at least len+1 bytes */

    int i, c;

    if (len > 0)
    {
        strncpy(dest, source, len);
    }
    dest[len] = '\0';

    len = (len > 0) ? (int)strlen(dest) : 0;

    for (i = ( len - 1 ); i >= 0 && 0 != ( c = source[i] ) && isspace( UCINT c ); i--);

    *first_space = dest + ((long long)i + 1); /* first blank or zero terminating byte in dest */ /* djb-rwth: cast operator added */

    return len; /* number of actually processed bytes excluding zero terminator */
}
    */
    // END INCHI C FUNCTION: MolfileStrnread
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileStrnread
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; UCINT casts through unsigned char.
    // INCHI✔️❌: The selected C locale classifies only ASCII space, tab, LF, VT, FF, and CR.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileStrnread

    let (actual_length, first_space_offset) = molfile_strnread_destination(
        heap, MolfileStringDestination::Stored(dest), MolfileFieldSource::Stored(source), len,
    )?;
    *first_space = dest.offset(
        i64::try_from(first_space_offset).map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
    )?;
    i32::try_from(actual_length).map_err(|_| SourceHeapError::SourceIntegerOverflow)
}

#[derive(Clone, Copy)]
enum MolfileFieldSource<'a> {
    Stored(SourceConstPointer<i8>),
    Borrowed(&'a [i8], usize),
}

impl MolfileFieldSource<'_> {
    fn bytes<'a>(&'a self, heap: &'a SourceHeap) -> Result<&'a [i8], SourceHeapError> {
        // INCHI✔️❗: Reached source view only, no copy/allocation/admission pass.
        match self {
            Self::Stored(pointer) => heap.slice(*pointer),
            Self::Borrowed(bytes, start) => bytes
                .get(*start..)
                .ok_or(SourceHeapError::PointerOutOfBounds),
        }
    }

    fn offset(self, amount: i64) -> Result<Self, SourceHeapError> {
        // Only the source nonnegative field prefix is projected here.
        match self {
            Self::Stored(pointer) => Ok(Self::Stored(pointer.offset(amount)?)),
            Self::Borrowed(bytes, start) => {
                let amount =
                    usize::try_from(amount).map_err(|_| SourceHeapError::PointerOffsetOverflow)?;
                Ok(Self::Borrowed(
                    bytes,
                    start
                        .checked_add(amount)
                        .ok_or(SourceHeapError::PointerOffsetOverflow)?,
                ))
            }
        }
    }
}

enum MolfileStringDestination<'a> {
    Stored(SourceMutPointer<i8>),
    Borrowed(&'a mut [i8]),
}

impl MolfileStringDestination<'_> {
    fn bytes_mut<'a>(
        &'a mut self,
        heap: &'a mut SourceHeap,
    ) -> Result<&'a mut [i8], SourceHeapError> {
        // INCHI source mol_fmt2.c:79:
        //     dest[len] = '\0';
        // INCHI✔️❗: Checked concrete destination view, O(1), no allocation.
        match self {
            Self::Stored(pointer) => heap.slice_mut(*pointer),
            Self::Borrowed(bytes) => Ok(bytes),
        }
    }
}

fn molfile_strnread_destination(
    heap: &mut SourceHeap,
    mut destination: MolfileStringDestination<'_>,
    source: MolfileFieldSource<'_>,
    len: i32,
) -> Result<(usize, usize), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:65 MolfileStrnread
    // INCHI✔️❌: complete source frame follows verbatim; checked source-heap access is materially slower than raw pointers.
    /*
    int MolfileStrnread(char *dest,
                        char *source,
                        int len,
                        char **first_space)
    {
        /* required len >= 0; dest must have at least len+1 bytes */

        int i, c;

        if (len > 0)
        {
            strncpy(dest, source, len);
        }
        dest[len] = '\0';

        len = (len > 0) ? (int)strlen(dest) : 0;

        for (i = ( len - 1 ); i >= 0 && 0 != ( c = source[i] ) && isspace( UCINT c ); i--);

        *first_space = dest + ((long long)i + 1); /* first blank or zero terminating byte in dest */ /* djb-rwth: cast operator added */

        return len; /* number of actually processed bytes excluding zero terminator */
    }
        */
    // END INCHI C FUNCTION: MolfileStrnread
    // INCHI✔️❗: One complete bounded strncpy/padding/trim kernel for stored
    // and borrowed destinations. Bounded first-NUL discovery supplies the same
    // resulting strlen; bulk copy/padding need no snapshot or arena temporary.
    // Re-read source after writes for native trailing-space/alias semantics.
    // O(len + reached trim), O(1) scratch, no allocation; checked copy/arena
    // constants remain OPEN INCHI-PERF-MOLFILE-STRING.
    let requested = usize::try_from(len).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if destination.bytes_mut(heap)?.len() <= requested {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let actual_length = if requested > 0 {
        let (available, nul) = {
            let bytes = source.bytes(heap)?;
            let available = requested.min(bytes.len());
            (
                available,
                memchr::memchr(0, bytemuck::cast_slice(&bytes[..available])),
            )
        };
        let copied = nul.map_or(available, |index| index + 1);
        match &mut destination {
            MolfileStringDestination::Stored(pointer) => match source {
                MolfileFieldSource::Stored(source) => {
                    heap.copy_nonoverlapping_prefix(*pointer, source, copied, true)?;
                }
                MolfileFieldSource::Borrowed(bytes, start) => {
                    heap.slice_mut(*pointer)?[..copied]
                        .copy_from_slice(&bytes[start..start + copied]);
                }
            },
            MolfileStringDestination::Borrowed(bytes) => {
                bytes[..copied].copy_from_slice(&source.bytes(heap)?[..copied]);
            }
        }
        // Preserve the previously written prefix on short unterminated input.
        if nul.is_none() && available < requested {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let actual = nul.unwrap_or(requested);
        destination.bytes_mut(heap)?[actual..=requested].fill(0);
        actual
    } else {
        destination.bytes_mut(heap)?[0] = 0;
        0
    };
    let mut first_space_offset = actual_length;
    if first_space_offset > 0 {
        let source_bytes = source.bytes(heap)?;
        while first_space_offset > 0 {
            let byte = *source_bytes
                .get(first_space_offset - 1)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if byte == 0 || !matches!(byte as u8, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r') {
                break;
            }
            first_space_offset -= 1;
        }
    }
    Ok((actual_length, first_space_offset))
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileReadField(
    heap: &mut SourceHeap,
    data: MolfileFieldData<'_>,
    field_len: i32,
    data_type: i32,
    line_ptr: &mut SourceMutPointer<i8>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:105 MolfileReadField
    // INCHI✔️❌: complete source frame follows verbatim; typed void-pointer replacement and checked source-heap access add overhead.
    /*
int MolfileReadField(void *data,
                     int field_len,
                     int data_type,
                     char **line_ptr)
{
    char *p = *line_ptr, *q, *p_end;
    int  i, c, len, ret = 1;
    long ldata;
    double ddata;

    int DEFINITE_LENGTH_FIELD = 0;
    int FIELD_ENDS_AT_FIRST_NON_DIGIT = 0;
    int TOO_LONG_FIELD = 0;

    if (field_len > MOL_FMT_MAX_VALUE_LEN)
    {
        TOO_LONG_FIELD = 1;
    }
    else if (field_len <= 0)
    {
        FIELD_ENDS_AT_FIRST_NON_DIGIT = 1;
    }
    else
    {
        DEFINITE_LENGTH_FIELD = 1;
    }

    switch (data_type)
    {
    case MOL_FMT_STRING_DATA:
        /* pass by all leading spaces */
        for (i = 0;
             i < field_len && 0 != (c = p[i]) && isspace(UCINT c);
             i++)
        {
            ;
        }

        len = MolfileStrnread((char *)data, &p[i], field_len - i, &q);

            ret = ( q - (char*) data );/* actual data length */
            *q = '\0'; /* add zero termination to data if it is not there yet*/
            *line_ptr += ( (long long)len + (long long)i ); /* ptr to the 1st byte of the next input field or to zero termination */ /* djb-rwth: cast operators added */
        break;

    case MOL_FMT_CHAR_INT_DATA:
    case MOL_FMT_SHORT_INT_DATA:
    case MOL_FMT_LONG_INT_DATA:
    {
        char str[MOL_FMT_MAX_VALUE_LEN + 1];
        ldata = 0L;
        if (TOO_LONG_FIELD)
        {
            ret = -1;
        }
        else if (DEFINITE_LENGTH_FIELD)
        {
            /* fixed length */
            *line_ptr += (len = MolfileStrnread(str, p, field_len, &q));

            *q = '\0';
            if (!len || !(q - str))
            {
                ret = 0; /* empty string */
            }
            else
            {
                if ((ldata = strtol(str, &p_end, 10), p_end != q))
                {
                    ret = -1; /* wrong data: incompletely interpreted */
                }
            }
        }
        else if (FIELD_ENDS_AT_FIRST_NON_DIGIT)
        {
            /* free format: field_len <= 0 */
            ldata = strtol(p, &p_end, 10);
            *line_ptr += (len = p_end - p);
            if (len == 0)
            {
                ret = 0;
            }
        }
        else
        {
            /* should not come here */
            ret = -1;
        }

        switch (data_type)
        {
        case MOL_FMT_CHAR_INT_DATA:
            if (SCHAR_MIN <= ldata && ldata <= SCHAR_MAX)
            {
                /* from || to &&: 11-19-96 */
                *(S_CHAR *)data = (S_CHAR)ldata;
            }
            else
            {
                *(S_CHAR *)data = (S_CHAR)0;
                ret = -1;
            }
            break;
        case MOL_FMT_SHORT_INT_DATA:
            if (SHRT_MIN <= ldata && ldata <= SHRT_MAX)
            {
                *(S_SHORT *)data = (S_SHORT)ldata;
            }
            else
            {
                *(S_SHORT *)data = (S_SHORT)0;
                ret = -1;
            }
            break;
        case MOL_FMT_LONG_INT_DATA:
            if (LONG_MIN < ldata && ldata < LONG_MAX)
            {
                *(long *)data = (long)ldata;
            }
            else
            {
                *(long *)data = 0L;
                ret = -1;
            }
            break;
        default:
            ret = -1;
        }
    } /* MOL_FMT_CHAR_INT_DATA... */
    break;

    case MOL_FMT_DOUBLE_DATA:
    case MOL_FMT_FLOAT_DATA:
    {
        char str[MOL_FMT_MAX_VALUE_LEN + 1];

        if (TOO_LONG_FIELD)
        {
            ret = -1;
            ddata = 0.0;
        }
        else if (DEFINITE_LENGTH_FIELD)
        {
            *line_ptr += (len = MolfileStrnread(str, p, field_len, &q));
            *q = '\0';
            if (!len || !(q - str))
            {
                /* empty string */
                ddata = 0.0;
                ret = 0;
            }
            else if ((ddata = strtod(str, &p_end), p_end != q))
            {
                /* wrong data */
                ret = -1;
            }
        }
        else if (FIELD_ENDS_AT_FIRST_NON_DIGIT)
        {
            /* free format */
            ddata = strtod(p, &p_end);
            *line_ptr += (len = p_end - p);
            if (len == 0)
            {
                ret = 0;
            }
        }
        else
        {
            /* should not come here */
            ret = -1; /* djb-rwth: addressing coverity ID #499478 -- see the original comment above */
        }

        switch (data_type)
        {

        case MOL_FMT_DOUBLE_DATA:
            if (ddata != HUGE_VAL && /*ldata*/ ddata != -HUGE_VAL)
            { /* replaced ldata with ddata 6-30-98 DCh */
                *(double *)data = ddata;
            }
            else
            {
                *(double *)data = 0.0;
                ret = -1;
            }
            break;

        case MOL_FMT_FLOAT_DATA:
            if (fabs(ddata) <= (double)FLT_MIN)
            {
                *(float *)data = 0.0;
            }
            else if (fabs(ddata) >= (double)FLT_MAX)
            {
                *(float *)data = 0.0;
                ret = -1;
            }
            else
            {
                *(float *)data = (float)ddata;
            }
            break;
        }
    } /* MOL_FMT_DOUBLE_DATA... */
    break;

    case MOL_FMT_JUMP_TO_RIGHT:
    {

        for (i = 0; i < field_len && p[i]; i++)
            ;

        *line_ptr += i;
        ret = i;
    }
    break;

    default:
        ret = -1;
    }

    return ret;
}
    */
    // END INCHI C FUNCTION: MolfileReadField
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileReadField
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64 uses signed 8-bit S_CHAR, 16-bit S_SHORT, and 64-bit long.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileReadField

    // INCHI✔️❗: Stored cursor keeps its original publication/error points;
    // one shared source kernel, no argument copy/allocation or parser duplicate.
    let p = MolfileFieldSource::Stored(line_ptr.as_const());
    molfile_read_field_source(heap, data, field_len, data_type, p, |amount| {
        *line_ptr = line_ptr.offset(amount)?;
        Ok(())
    })
}

pub(crate) fn molfile_read_field_bytes(
    heap: &mut SourceHeap,
    data: MolfileFieldData<'_>,
    field_len: i32,
    data_type: i32,
    bytes: &[i8],
    position: &mut usize,
) -> Result<i32, SourceHeapError> {
    // INCHI source mol_fmt1.c485..488 passes its local line cursor to this same
    // source field algorithm. O(1) carrier/no arena registration or allocation;
    // borrowed input admission occurs only at reached reads, checked cost OPEN.
    let source = MolfileFieldSource::Borrowed(bytes, *position);
    molfile_read_field_source(heap, data, field_len, data_type, source, |amount| {
        let amount = usize::try_from(amount).map_err(|_| SourceHeapError::PointerOffsetOverflow)?;
        *position = position
            .checked_add(amount)
            .ok_or(SourceHeapError::PointerOffsetOverflow)?;
        Ok(())
    })
}

fn molfile_read_field_source(
    heap: &mut SourceHeap,
    data: MolfileFieldData<'_>,
    field_len: i32,
    data_type: i32,
    p: MolfileFieldSource<'_>,
    mut advance: impl FnMut(i64) -> Result<(), SourceHeapError>,
) -> Result<i32, SourceHeapError> {
    // VERBATIM SOURCE FRAGMENT: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:105..328 C393 shared field kernel
    // int MolfileReadField(void *data,
    //                      int field_len,
    //                      int data_type,
    //                      char **line_ptr)
    // {
    //     char *p = *line_ptr, *q, *p_end;
    //     int  i, c, len, ret = 1;
    //     long ldata;
    //     double ddata;
    //
    //     int DEFINITE_LENGTH_FIELD = 0;
    //     int FIELD_ENDS_AT_FIRST_NON_DIGIT = 0;
    //     int TOO_LONG_FIELD = 0;
    //
    //     if (field_len > MOL_FMT_MAX_VALUE_LEN)
    //     {
    //         TOO_LONG_FIELD = 1;
    //     }
    //     else if (field_len <= 0)
    //     {
    //         FIELD_ENDS_AT_FIRST_NON_DIGIT = 1;
    //     }
    //     else
    //     {
    //         DEFINITE_LENGTH_FIELD = 1;
    //     }
    //
    //     switch (data_type)
    //     {
    //     case MOL_FMT_STRING_DATA:
    //         /* pass by all leading spaces */
    //         for (i = 0;
    //              i < field_len && 0 != (c = p[i]) && isspace(UCINT c);
    //              i++)
    //         {
    //             ;
    //         }
    //
    //         len = MolfileStrnread((char *)data, &p[i], field_len - i, &q);
    //
    //             ret = ( q - (char*) data );/* actual data length */
    //             *q = '\0'; /* add zero termination to data if it is not there yet*/
    //             *line_ptr += ( (long long)len + (long long)i ); /* ptr to the 1st byte of the next input field or to zero termination */ /* djb-rwth: cast operators added */
    //         break;
    //
    //     case MOL_FMT_CHAR_INT_DATA:
    //     case MOL_FMT_SHORT_INT_DATA:
    //     case MOL_FMT_LONG_INT_DATA:
    //     {
    //         char str[MOL_FMT_MAX_VALUE_LEN + 1];
    //         ldata = 0L;
    //         if (TOO_LONG_FIELD)
    //         {
    //             ret = -1;
    //         }
    //         else if (DEFINITE_LENGTH_FIELD)
    //         {
    //             /* fixed length */
    //             *line_ptr += (len = MolfileStrnread(str, p, field_len, &q));
    //
    //             *q = '\0';
    //             if (!len || !(q - str))
    //             {
    //                 ret = 0; /* empty string */
    //             }
    //             else
    //             {
    //                 if ((ldata = strtol(str, &p_end, 10), p_end != q))
    //                 {
    //                     ret = -1; /* wrong data: incompletely interpreted */
    //                 }
    //             }
    //         }
    //         else if (FIELD_ENDS_AT_FIRST_NON_DIGIT)
    //         {
    //             /* free format: field_len <= 0 */
    //             ldata = strtol(p, &p_end, 10);
    //             *line_ptr += (len = p_end - p);
    //             if (len == 0)
    //             {
    //                 ret = 0;
    //             }
    //         }
    //         else
    //         {
    //             /* should not come here */
    //             ret = -1;
    //         }
    //
    //         switch (data_type)
    //         {
    //         case MOL_FMT_CHAR_INT_DATA:
    //             if (SCHAR_MIN <= ldata && ldata <= SCHAR_MAX)
    //             {
    //                 /* from || to &&: 11-19-96 */
    //                 *(S_CHAR *)data = (S_CHAR)ldata;
    //             }
    //             else
    //             {
    //                 *(S_CHAR *)data = (S_CHAR)0;
    //                 ret = -1;
    //             }
    //             break;
    //         case MOL_FMT_SHORT_INT_DATA:
    //             if (SHRT_MIN <= ldata && ldata <= SHRT_MAX)
    //             {
    //                 *(S_SHORT *)data = (S_SHORT)ldata;
    //             }
    //             else
    //             {
    //                 *(S_SHORT *)data = (S_SHORT)0;
    //                 ret = -1;
    //             }
    //             break;
    //         case MOL_FMT_LONG_INT_DATA:
    //             if (LONG_MIN < ldata && ldata < LONG_MAX)
    //             {
    //                 *(long *)data = (long)ldata;
    //             }
    //             else
    //             {
    //                 *(long *)data = 0L;
    //                 ret = -1;
    //             }
    //             break;
    //         default:
    //             ret = -1;
    //         }
    //     } /* MOL_FMT_CHAR_INT_DATA... */
    //     break;
    //
    //     case MOL_FMT_DOUBLE_DATA:
    //     case MOL_FMT_FLOAT_DATA:
    //     {
    //         char str[MOL_FMT_MAX_VALUE_LEN + 1];
    //
    //         if (TOO_LONG_FIELD)
    //         {
    //             ret = -1;
    //             ddata = 0.0;
    //         }
    //         else if (DEFINITE_LENGTH_FIELD)
    //         {
    //             *line_ptr += (len = MolfileStrnread(str, p, field_len, &q));
    //             *q = '\0';
    //             if (!len || !(q - str))
    //             {
    //                 /* empty string */
    //                 ddata = 0.0;
    //                 ret = 0;
    //             }
    //             else if ((ddata = strtod(str, &p_end), p_end != q))
    //             {
    //                 /* wrong data */
    //                 ret = -1;
    //             }
    //         }
    //         else if (FIELD_ENDS_AT_FIRST_NON_DIGIT)
    //         {
    //             /* free format */
    //             ddata = strtod(p, &p_end);
    //             *line_ptr += (len = p_end - p);
    //             if (len == 0)
    //             {
    //                 ret = 0;
    //             }
    //         }
    //         else
    //         {
    //             /* should not come here */
    //             ret = -1; /* djb-rwth: addressing coverity ID #499478 -- see the original comment above */
    //         }
    //
    //         switch (data_type)
    //         {
    //
    //         case MOL_FMT_DOUBLE_DATA:
    //             if (ddata != HUGE_VAL && /*ldata*/ ddata != -HUGE_VAL)
    //             { /* replaced ldata with ddata 6-30-98 DCh */
    //                 *(double *)data = ddata;
    //             }
    //             else
    //             {
    //                 *(double *)data = 0.0;
    //                 ret = -1;
    //             }
    //             break;
    //
    //         case MOL_FMT_FLOAT_DATA:
    //             if (fabs(ddata) <= (double)FLT_MIN)
    //             {
    //                 *(float *)data = 0.0;
    //             }
    //             else if (fabs(ddata) >= (double)FLT_MAX)
    //             {
    //                 *(float *)data = 0.0;
    //                 ret = -1;
    //             }
    //             else
    //             {
    //                 *(float *)data = (float)ddata;
    //             }
    //             break;
    //         }
    //     } /* MOL_FMT_DOUBLE_DATA... */
    //     break;
    //
    //     case MOL_FMT_JUMP_TO_RIGHT:
    //     {
    //
    //         for (i = 0; i < field_len && p[i]; i++)
    //             ;
    //
    //         *line_ptr += i;
    //         ret = i;
    //     }
    //     break;
    //
    //     default:
    //         ret = -1;
    //     }
    //
    //     return ret;
    // }
    // END VERBATIM SOURCE FRAGMENT: C393 shared field kernel
    // INCHI✔️❗: Sole native field/control/type/errno algorithm for canonical
    // pointer and local borrowed source. Advance callback runs at the original
    // publication point before conversion; partial effects survive Result errors.
    // Stored copy alias and numeric wrappers unchanged, local bytes use their
    // same canonical conversion cores. O(field)/fixed33-byte scratch, no new
    // allocation/clone/preflight/algorithm; whole-call costs remain OPEN.
    let too_long = field_len > MOL_FMT_MAX_VALUE_LEN as i32;
    let definite = field_len > 0 && !too_long;

    if data_type == i32::from(MOL_FMT_STRING_DATA) {
        let MolfileFieldData::String(output) = data else {
            return Err(SourceHeapError::UnsupportedSourceBehavior);
        };
        let width =
            usize::try_from(field_len).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        // INCHI source mol_fmt2.c:136..141:
        //         for (i = 0;
        //              i < field_len && 0 != (c = p[i]) && isspace(UCINT c);
        //              i++)
        //         {
        //             ;
        //         }
        // INCHI✔️❗: Admit only reached leading bytes; destination capacity
        // applies to field_len-i, as in the complete shared bounded reader.
        // Borrow output directly, preserving its untouched tail/partial writes.
        // No Vec, arena allocation or whole writeback; whole-call cost OPEN.
        let mut skipped = 0_usize;
        if width > 0 {
            let input = p.bytes(heap)?;
            while skipped < width {
                let byte = *input
                    .get(skipped)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if byte == 0 || !matches!(byte as u8, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r') {
                    break;
                }
                skipped += 1;
            }
        }
        let (read, q) = molfile_strnread_destination(
            heap,
            MolfileStringDestination::Borrowed(&mut *output),
            p.offset(skipped as i64)?,
            field_len - skipped as i32,
        )?;
        output[q] = 0;
        advance(read as i64 + skipped as i64)?;
        return i32::try_from(q).map_err(|_| SourceHeapError::SourceIntegerOverflow);
    }

    let integer = data_type == i32::from(MOL_FMT_CHAR_INT_DATA)
        || data_type == i32::from(MOL_FMT_SHORT_INT_DATA)
        || data_type == i32::from(MOL_FMT_LONG_INT_DATA);
    if integer {
        let mut ret = 1;
        let mut value = 0_i64;
        if too_long {
            ret = -1;
        } else if definite {
            // INCHI source mol_fmt2.c:154..177: native fixed stack field.
            // INCHI✔️❗: Native33byte local storage and borrowed canonical
            // parser; cursor advances before conversion. No Vec, allocation,
            // arena registration/free or extra failure. Parser guards remain.
            let mut temporary = [0_i8; MOL_FMT_MAX_VALUE_LEN as usize + 1];
            let (read, q) = molfile_strnread_destination(
                heap,
                MolfileStringDestination::Borrowed(&mut temporary),
                p,
                field_len,
            )?;
            advance(read as i64)?;
            temporary[q] = 0;
            if read == 0 || q == 0 {
                ret = 0;
            } else {
                let (parsed, end, overflowed) = inchi_strtol_bytes(&temporary, 10)?;
                value = parsed;
                if overflowed {
                    heap.set_source_errno(34);
                }
                if end != q {
                    ret = -1;
                }
            }
        } else {
            let read = match p {
                MolfileFieldSource::Stored(pointer) => {
                    let mut end = SourceConstPointer::null();
                    value = inchi_strtol(heap, pointer, Some(&mut end), 10)?;
                    let read = i64::try_from(
                        heap.slice(pointer)?
                            .len()
                            .checked_sub(heap.slice(end)?.len())
                            .ok_or(SourceHeapError::PointerDifferenceOverflow)?,
                    )
                    .map_err(|_| SourceHeapError::PointerDifferenceOverflow)?;
                    read
                }
                MolfileFieldSource::Borrowed(_, _) => {
                    let (parsed, read, overflowed) = inchi_strtol_bytes(p.bytes(heap)?, 10)?;
                    value = parsed;
                    if overflowed {
                        heap.set_source_errno(34);
                    }
                    i64::try_from(read).map_err(|_| SourceHeapError::PointerDifferenceOverflow)?
                }
            };
            advance(read)?;
            if read == 0 {
                ret = 0;
            }
        }
        return match (data_type, data) {
            (x, MolfileFieldData::Char(out)) if x == i32::from(MOL_FMT_CHAR_INT_DATA) => {
                if (-128..=127).contains(&value) {
                    *out = value as i8;
                } else {
                    *out = 0;
                    ret = -1;
                }
                Ok(ret)
            }
            (x, MolfileFieldData::Short(out)) if x == i32::from(MOL_FMT_SHORT_INT_DATA) => {
                if (-32768..=32767).contains(&value) {
                    *out = value as i16;
                } else {
                    *out = 0;
                    ret = -1;
                }
                Ok(ret)
            }
            (x, MolfileFieldData::ShortIntoInt(out)) if x == i32::from(MOL_FMT_SHORT_INT_DATA) => {
                let narrowed = if (-32768..=32767).contains(&value) {
                    value as i16
                } else {
                    ret = -1;
                    0
                };
                let mut bytes = out.to_ne_bytes();
                bytes[..2].copy_from_slice(&narrowed.to_ne_bytes());
                *out = i32::from_ne_bytes(bytes);
                Ok(ret)
            }
            (x, MolfileFieldData::Long(out)) if x == i32::from(MOL_FMT_LONG_INT_DATA) => {
                if i64::MIN < value && value < i64::MAX {
                    *out = value;
                } else {
                    *out = 0;
                    ret = -1;
                }
                Ok(ret)
            }
            _ => Err(SourceHeapError::UnsupportedSourceBehavior),
        };
    }

    let floating =
        data_type == i32::from(MOL_FMT_DOUBLE_DATA) || data_type == i32::from(MOL_FMT_FLOAT_DATA);
    if floating {
        let mut ret = 1;
        let value;
        if too_long {
            ret = -1;
            value = 0.0;
        } else if definite {
            // INCHI source mol_fmt2.c:238..256: native fixed stack field.
            // INCHI✔️❗: Native33byte local storage and borrowed canonical
            // parser; cursor advances before conversion. No Vec, allocation,
            // arena registration/free or extra failure. Parser guards remain.
            let mut temporary = [0_i8; MOL_FMT_MAX_VALUE_LEN as usize + 1];
            let (read, q) = molfile_strnread_destination(
                heap,
                MolfileStringDestination::Borrowed(&mut temporary),
                p,
                field_len,
            )?;
            advance(read as i64)?;
            temporary[q] = 0;
            if read == 0 || q == 0 {
                value = 0.0;
                ret = 0;
            } else {
                let mut source_errno = heap.source_errno();
                let (parsed, end) = source_strtod_bytes_with_end(&temporary, &mut source_errno)?;
                heap.set_source_errno(source_errno);
                value = parsed;
                if end != q {
                    ret = -1;
                }
            }
        } else {
            let read = match p {
                MolfileFieldSource::Stored(pointer) => {
                    let mut end = SourceConstPointer::null();
                    value = inchi_strtod(heap, pointer, Some(&mut end))?;
                    let read = i64::try_from(
                        heap.slice(pointer)?
                            .len()
                            .checked_sub(heap.slice(end)?.len())
                            .ok_or(SourceHeapError::PointerDifferenceOverflow)?,
                    )
                    .map_err(|_| SourceHeapError::PointerDifferenceOverflow)?;
                    read
                }
                MolfileFieldSource::Borrowed(_, _) => {
                    let mut source_errno = heap.source_errno();
                    let (parsed, read) =
                        source_strtod_bytes_with_end(p.bytes(heap)?, &mut source_errno)?;
                    heap.set_source_errno(source_errno);
                    value = parsed;
                    i64::try_from(read).map_err(|_| SourceHeapError::PointerDifferenceOverflow)?
                }
            };
            advance(read)?;
            if read == 0 {
                ret = 0;
            }
        }
        return match (data_type, data) {
            (x, MolfileFieldData::Double(out)) if x == i32::from(MOL_FMT_DOUBLE_DATA) => {
                if value.is_infinite() {
                    *out = 0.0;
                    ret = -1;
                } else {
                    *out = value;
                }
                Ok(ret)
            }
            (x, MolfileFieldData::Float(out)) if x == i32::from(MOL_FMT_FLOAT_DATA) => {
                let magnitude = value.abs();
                if magnitude <= f64::from(f32::MIN_POSITIVE) {
                    *out = 0.0;
                } else if magnitude >= f64::from(f32::MAX) {
                    *out = 0.0;
                    ret = -1;
                } else {
                    *out = value as f32;
                }
                Ok(ret)
            }
            _ => Err(SourceHeapError::UnsupportedSourceBehavior),
        };
    }

    if data_type == i32::from(MOL_FMT_JUMP_TO_RIGHT) {
        // INCHI source mol_fmt2.c:315..318:
        //         for (i = 0; i < field_len && p[i]; i++)
        //             ;
        //
        //         *line_ptr += i;
        // INCHI✔️❗: Nonpositive width does not admit unused source; only
        // reached bytes are checked. O(reached)/O(1), no allocation.
        let mut count = 0;
        if field_len > 0 {
            let input = p.bytes(heap)?;
            while count < field_len {
                if *input
                    .get(count as usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    == 0
                {
                    break;
                }
                count += 1;
            }
        }
        advance(i64::from(count))?;
        return Ok(count);
    }
    Ok(-1)
}

#[allow(non_snake_case)]
pub(crate) fn MolfileExtractStrucNum(
    heap: &mut SourceHeap,
    pHdr: Option<&MOL_FMT_HEADER_BLOCK>,
) -> Result<i64, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:333 MolfileExtractStrucNum
    // INCHI✔️❌: complete source frame follows verbatim; SourceHeap pointer modeling requires one temporary checked byte buffer.
    /*
    long MolfileExtractStrucNum(MOL_FMT_HEADER_BLOCK *pHdr)
    {
        static char sStruct[] = "Structure #";
        static char sINCHI[] = INCHI_NAME;
        long   lMolfileNumber = 0;
        char   *p, *q = NULL;

        if (pHdr)
        {
            if (!inchi_memicmp(pHdr->molname, sStruct, sizeof(sStruct) - 1))
            {
                p = pHdr->molname + sizeof(sStruct) - 1;
                lMolfileNumber = strtol(p, &q, 10);
                p = pHdr->line2;
                if (!q || *q ||
                    inchi_memicmp(p, sINCHI, sizeof(sINCHI) - 1) ||
                    !strstr(p + sizeof(sINCHI) - 1, "SDfile Output"))
                {
                    lMolfileNumber = 0;
                }
            }
        }

        return lMolfileNumber;
    }
        */
    // END INCHI C FUNCTION: MolfileExtractStrucNum
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileExtractStrucNum
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; long is signed 64-bit and INCHI_NAME is "InChI".
    // INCHI✔️❌: inchi_memicmp and inchi_strtol are completed source ports; strstr is reproduced as a case-sensitive contiguous byte search after the five-byte INCHI_NAME prefix.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileExtractStrucNum

    let Some(header) = pHdr else {
        return Ok(0);
    };

    // INCHI source mol_fmt2.c:342..349:
    //         if (!inchi_memicmp(pHdr->molname, sStruct, sizeof(sStruct) - 1))
    //         {
    //             p = pHdr->molname + sizeof(sStruct) - 1;
    //             lMolfileNumber = strtol(p, &q, 10);
    //             p = pHdr->line2;
    //             if (!q || *q ||
    //                 inchi_memicmp(p, sINCHI, sizeof(sINCHI) - 1) ||
    //                 !strstr(p + sizeof(sINCHI) - 1, "SDfile Output"))
    // INCHI✔️❗: Borrow the two independently terminated header fields;
    // mol_fmt1.c485..488/534..537 and mol_fmt2.c145 insert their NULs.
    // Never concatenate fields with artificial source-literal terminators.
    // Malformed reached field ends retain MissingNulTerminator. No preflight
    // for an unreached field. Conversion publishes errno before later gates.
    // O(reached numeric token + reached search prefix), O(1) scratch, no
    // allocation/clone/full-field copy. Checked/libc constants remain OPEN.
    if crate::source::base::util::inchi_memicmp_bytes(
        &header.molname,
        bytemuck::cast_slice(b"Structure #"),
        11,
    )? != 0
    {
        return Ok(0);
    }

    let numeric = &header.molname[11..];
    let (molfile_number, end, overflowed) = inchi_strtol_bytes(numeric, 10)?;
    if overflowed {
        heap.set_source_errno(34);
    }
    if *numeric
        .get(end)
        .ok_or(SourceHeapError::MissingNulTerminator)?
        != 0
        || crate::source::base::util::inchi_memicmp_bytes(
            &header.line2,
            bytemuck::cast_slice(b"InChI"),
            5,
        )? != 0
    {
        return Ok(0);
    }

    // Native strstr stops at a match or the first NUL. The fixed13-byte
    // needle bounds each candidate comparison; skip noncandidate bytes with
    // the existing allocation-free memchr2 rather than a full strlen scan.
    let bytes: &[u8] = bytemuck::cast_slice(&header.line2[5..]);
    let needle = b"SDfile Output";
    let mut position = 0_usize;
    loop {
        let offset = memchr::memchr2(0, needle[0], &bytes[position..])
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        position += offset;
        if bytes[position] == 0 {
            return Ok(0);
        }
        let mut matched = true;
        for (index, expected) in needle.iter().enumerate().skip(1) {
            let actual = *bytes
                .get(position + index)
                .ok_or(SourceHeapError::MissingNulTerminator)?;
            if actual == 0 {
                return Ok(0);
            }
            if actual != *expected {
                matched = false;
                break;
            }
        }
        if matched {
            return Ok(molfile_number);
        }
        position += 1;
    }
}

#[allow(non_snake_case)]
pub(crate) fn MolfileHasNoChemStruc(mfdata: Option<&MOL_FMT_DATA>) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:362 MolfileHasNoChemStruc
    // INCHI✔️✔️: int MolfileHasNoChemStruc(MOL_FMT_DATA *mfdata)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     if (!mfdata || !mfdata->ctab.atoms)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (mfdata->ctab.n_atoms <= 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (0 < mfdata->ctab.n_bonds && !mfdata->ctab.bonds)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: MolfileHasNoChemStruc

    let Some(mfdata) = mfdata else {
        return 1;
    };
    if mfdata.ctab.atoms.is_null() {
        return 1;
    }
    if mfdata.ctab.n_atoms <= 0 {
        return 1;
    }
    if mfdata.ctab.n_bonds > 0 && mfdata.ctab.bonds.is_null() {
        return 1;
    }
    0
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MolfileGetXYZDimAndNormFactors(
    heap: &SourceHeap,
    mfdata: Option<&MOL_FMT_DATA>,
    find_norm_factors: i32,
    x0: &mut f64,
    y0: &mut f64,
    z0: &mut f64,
    xmin: &mut f64,
    ymin: &mut f64,
    zmin: &mut f64,
    scaler: &mut f64,
    err: &mut i32,
    mut p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:474 MolfileGetXYZDimAndNormFactors
    // INCHI✔️❌: complete source frame follows verbatim; checked heap access and warning materialization add overhead.
    /*
int MolfileGetXYZDimAndNormFactors(MOL_FMT_DATA *mfdata,
                                   int find_norm_factors,
                                   double *x0,
                                   double *y0,
                                   double *z0,
                                   double *xmin,
                                   double *ymin,
                                   double *zmin,
                                   double *scaler,
                                   int *err,
                                   char *pStrErr)

{
    int i;
    int num_dimensions = 0, num_atoms, num_bonds;
    double max_x = -1.0e32, max_y = -1.0e32, max_z = -1.0e32;
    double min_x = 1.0e32, min_y = 1.0e32, min_z = 1.0e32;
    double macheps = 1.0e-10, small_coeff = 0.00001;
    double x_coeff, y_coeff, z_coeff, coeff = 1.0, average_bond_length;

    *x0 = MIN_STDATA_X_COORD;
    *y0 = MIN_STDATA_Y_COORD;
    *z0 = MIN_STDATA_Z_COORD;
    *xmin = *ymin = *zmin = 0.0;
    *scaler = coeff;

    if (MolfileHasNoChemStruc(mfdata))
    {
        goto exit_function;
    }

    num_atoms = mfdata->ctab.n_atoms;
    for (i = 0; i < num_atoms; i++)
    {
        max_x = inchi_max(mfdata->ctab.atoms[i].fx, max_x);
        min_x = inchi_min(mfdata->ctab.atoms[i].fx, min_x);
        max_y = inchi_max(mfdata->ctab.atoms[i].fy, max_y);
        min_y = inchi_min(mfdata->ctab.atoms[i].fy, min_y);
        max_z = inchi_max(mfdata->ctab.atoms[i].fz, max_z);
        min_z = inchi_min(mfdata->ctab.atoms[i].fz, min_z);
    }

    num_bonds = 0;
    average_bond_length = 0.0;
    for (i = 0; i < mfdata->ctab.n_bonds; i++)
    {
        double dx, dy, dz;
        int a1 = mfdata->ctab.bonds[i].atnum1 - 1;
        int a2 = mfdata->ctab.bonds[i].atnum2 - 1;

        if (a1 < 0 || a1 >= num_atoms ||
            a2 < 0 || a2 >= num_atoms ||
            a1 == a2)
        {
            *err |= 1; /*  bond for invalid atom number(s); ignored */
            TREAT_ERR(*err, 0, "Bond to nonexistent atom");
            continue;
        }

        dx = mfdata->ctab.atoms[a1].fx - mfdata->ctab.atoms[a2].fx;
        dy = mfdata->ctab.atoms[a1].fy - mfdata->ctab.atoms[a2].fy;
        dz = mfdata->ctab.atoms[a1].fz - mfdata->ctab.atoms[a2].fz;

        average_bond_length += sqrt(dx*dx + dy*dy + dz*dz);
        num_bonds++;
    }

    if (max_x - min_x <= small_coeff * (fabs(max_x) + fabs(min_x)))
    {
        x_coeff = 0.0;
    }
    else
    {
        x_coeff = (MAX_STDATA_X_COORD - MIN_STDATA_X_COORD) / (max_x - min_x);
    }

    if (max_y - min_y <= small_coeff * (fabs(max_y) + fabs(min_y)))
    {
        y_coeff = 0.0;
    }
    else
    {
        y_coeff = (MAX_STDATA_Y_COORD - MIN_STDATA_Y_COORD) / (max_y - min_y);
    }

    if (max_z - min_z <= small_coeff * (fabs(max_z) + fabs(min_z)))
    {
        z_coeff = 0.0;
    }
    else
    {
        z_coeff = (MAX_STDATA_Z_COORD - MIN_STDATA_Z_COORD) / (max_z - min_z);
    }

    num_dimensions = ((x_coeff > macheps || y_coeff > macheps) && fabs(z_coeff) < macheps)
                        ? 2
                        : ( fabs( z_coeff ) > macheps ) ? 3 : 0;


    if (!find_norm_factors)
    {
        goto exit_function;
    }

    /* Find normalization parameters */
    switch (num_dimensions)
    {
    case 0:
        coeff = 0.0;
        break;

    case 2:
        /* choose the smallest stretching coefficient */
        if (x_coeff > macheps && y_coeff > macheps)
        {
            coeff = inchi_min(x_coeff, y_coeff);
        }
        else if (x_coeff > macheps)
        {
            coeff = x_coeff;
        }
        else if (y_coeff > macheps)
        {
            coeff = y_coeff;
        }
        else
        {
            coeff = 1.0;
        }
        break;

    case 3:
        /* choose the smallest stretching coefficient */
        if (x_coeff > macheps && y_coeff > macheps)
        {
            coeff = inchi_min(x_coeff, y_coeff);
            coeff = inchi_min(coeff, z_coeff);
        }
        else if (x_coeff > macheps)
        {
            coeff = inchi_min(x_coeff, z_coeff);
        }
        else if (y_coeff > macheps)
        {
            coeff = inchi_min(y_coeff, z_coeff);
        }
        else
        {
            coeff = z_coeff;
        }
        break;

    default:
        coeff = 0.0;
    }

    if (num_bonds > 0)
    {

        average_bond_length /= (double)num_bonds;
        if (average_bond_length * coeff > MAX_STDATA_AVE_BOND_LENGTH)
        {
            coeff = MAX_STDATA_AVE_BOND_LENGTH / average_bond_length; /* avoid too long bonds */
        }
        else if (average_bond_length * coeff < macheps)
        {
            coeff = 1.0; /* all lengths are of zero length */
        }
        else if (average_bond_length * coeff < MIN_STDATA_AVE_BOND_LENGTH)
        {
            coeff = MIN_STDATA_AVE_BOND_LENGTH / average_bond_length; /* avoid too short bonds */
        }
    }

exit_function:;

    *x0 = min_x;
    *y0 = min_y;
    *z0 = min_z;
    *xmin = MIN_STDATA_X_COORD;
    *ymin = MIN_STDATA_Y_COORD;
    *zmin = MIN_STDATA_Z_COORD;
    *scaler = coeff;

    return num_dimensions;
}
    */
    // END INCHI C FUNCTION: MolfileGetXYZDimAndNormFactors
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileGetXYZDimAndNormFactors
    // INCHI✔️❌: inchi_max(a,b)=(((a)>(b))?(a):(b)); inchi_min(a,b)=(((a)<(b))?(a):(b)).
    // INCHI✔️❌: MIN_STDATA_*_COORD=0.0; MAX_STDATA_*_COORD=256.0.
    // INCHI✔️❌: MAX_STDATA_AVE_BOND_LENGTH=20.0; MIN_STDATA_AVE_BOND_LENGTH=10.0.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileGetXYZDimAndNormFactors

    const MIN_COORD: f64 = 0.0;
    const MAX_COORD: f64 = 256.0;
    const MAX_AVERAGE_BOND_LENGTH: f64 = 20.0;
    const MIN_AVERAGE_BOND_LENGTH: f64 = 10.0;
    const MACHINE_EPSILON: f64 = 1.0e-10;
    const SMALL_COEFFICIENT: f64 = 0.00001;

    let source_max = |left: f64, right: f64| if left > right { left } else { right };
    let source_min = |left: f64, right: f64| if left < right { left } else { right };

    *x0 = MIN_COORD;
    *y0 = MIN_COORD;
    *z0 = MIN_COORD;
    *xmin = 0.0;
    *ymin = 0.0;
    *zmin = 0.0;
    let mut coefficient = 1.0;
    *scaler = coefficient;

    let mut max_x = -1.0e32;
    let mut max_y = -1.0e32;
    let mut max_z = -1.0e32;
    let mut min_x = 1.0e32;
    let mut min_y = 1.0e32;
    let mut min_z = 1.0e32;
    let mut num_dimensions = 0;

    if let Some(mfdata) = mfdata.filter(|data| MolfileHasNoChemStruc(Some(data)) == 0) {
        let atom_count = usize::try_from(mfdata.ctab.n_atoms)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let atoms = heap
            .slice(mfdata.ctab.atoms.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for atom in atoms {
            max_x = source_max(atom.fx, max_x);
            min_x = source_min(atom.fx, min_x);
            max_y = source_max(atom.fy, max_y);
            min_y = source_min(atom.fy, min_y);
            max_z = source_max(atom.fz, max_z);
            min_z = source_min(atom.fz, min_z);
        }

        let source_bond_count = if mfdata.ctab.n_bonds > 0 {
            usize::try_from(mfdata.ctab.n_bonds)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
        } else {
            0
        };
        let bonds = if source_bond_count == 0 {
            &[][..]
        } else {
            heap.slice(mfdata.ctab.bonds.as_const())?
                .get(..source_bond_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        };
        let mut valid_bonds = 0_i32;
        let mut average_bond_length = 0.0;
        for bond in bonds {
            let a1 = i32::from(bond.atnum1).wrapping_sub(1);
            let a2 = i32::from(bond.atnum2).wrapping_sub(1);
            if a1 < 0
                || a1 >= mfdata.ctab.n_atoms
                || a2 < 0
                || a2 >= mfdata.ctab.n_atoms
                || a1 == a2
            {
                *err |= 1;
                add_error_message_bytes(p_str_err.as_deref_mut(), b"Bond to nonexistent atom")?;
                continue;
            }
            let first = &atoms[a1 as usize];
            let second = &atoms[a2 as usize];
            let dx = first.fx - second.fx;
            let dy = first.fy - second.fy;
            let dz = first.fz - second.fz;
            average_bond_length += (dx * dx + dy * dy + dz * dz).sqrt();
            valid_bonds = valid_bonds.wrapping_add(1);
        }

        let x_difference = max_x - min_x;
        let x_coefficient = if x_difference
            <= SMALL_COEFFICIENT * (max_x.abs() + min_x.abs())
        {
            0.0
        } else {
            (MAX_COORD - MIN_COORD) / x_difference
        };
        let y_difference = max_y - min_y;
        let y_coefficient = if y_difference
            <= SMALL_COEFFICIENT * (max_y.abs() + min_y.abs())
        {
            0.0
        } else {
            (MAX_COORD - MIN_COORD) / y_difference
        };
        let z_difference = max_z - min_z;
        let z_coefficient = if z_difference
            <= SMALL_COEFFICIENT * (max_z.abs() + min_z.abs())
        {
            0.0
        } else {
            (MAX_COORD - MIN_COORD) / z_difference
        };

        num_dimensions = if (x_coefficient > MACHINE_EPSILON
            || y_coefficient > MACHINE_EPSILON)
            && z_coefficient.abs() < MACHINE_EPSILON
        {
            2
        } else if z_coefficient.abs() > MACHINE_EPSILON {
            3
        } else {
            0
        };

        if find_norm_factors != 0 {
            coefficient = match num_dimensions {
                0 => 0.0,
                2 => {
                    if x_coefficient > MACHINE_EPSILON && y_coefficient > MACHINE_EPSILON {
                        source_min(x_coefficient, y_coefficient)
                    } else if x_coefficient > MACHINE_EPSILON {
                        x_coefficient
                    } else if y_coefficient > MACHINE_EPSILON {
                        y_coefficient
                    } else {
                        1.0
                    }
                }
                3 => {
                    if x_coefficient > MACHINE_EPSILON && y_coefficient > MACHINE_EPSILON {
                        source_min(source_min(x_coefficient, y_coefficient), z_coefficient)
                    } else if x_coefficient > MACHINE_EPSILON {
                        source_min(x_coefficient, z_coefficient)
                    } else if y_coefficient > MACHINE_EPSILON {
                        source_min(y_coefficient, z_coefficient)
                    } else {
                        z_coefficient
                    }
                }
                _ => 0.0,
            };

            if valid_bonds > 0 {
                average_bond_length /= f64::from(valid_bonds);
                if average_bond_length * coefficient > MAX_AVERAGE_BOND_LENGTH {
                    coefficient = MAX_AVERAGE_BOND_LENGTH / average_bond_length;
                } else if average_bond_length * coefficient < MACHINE_EPSILON {
                    coefficient = 1.0;
                } else if average_bond_length * coefficient < MIN_AVERAGE_BOND_LENGTH {
                    coefficient = MIN_AVERAGE_BOND_LENGTH / average_bond_length;
                }
            }
        }
    }

    *x0 = min_x;
    *y0 = min_y;
    *z0 = min_z;
    *xmin = MIN_COORD;
    *ymin = MIN_COORD;
    *zmin = MIN_COORD;
    *scaler = coefficient;
    Ok(num_dimensions)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn FreeMolfileData(
    heap: &mut SourceHeap,
    mfdata: SourceMutPointer<MOL_FMT_DATA>,
) -> Result<SourceMutPointer<MOL_FMT_DATA>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:664 FreeMolfileData
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap ownership traversal adds overhead.
    /*
MOL_FMT_DATA *FreeMolfileData(MOL_FMT_DATA *mfdata)
{
    if (mfdata)
    {

        if (mfdata->ctab.atoms)
        {
            inchi_free(mfdata->ctab.atoms);
        }

        if (mfdata->ctab.bonds)
        {
            inchi_free(mfdata->ctab.bonds);
        }

        if (mfdata->ctab.coords)
        {
            inchi_free(mfdata->ctab.coords);
        }

        /*if ( 0!=mfdata->ctab.sgroups.used )*/
        MolFmtSgroups_Free(&(mfdata->ctab.sgroups));

        if (mfdata->ctab.v3000)
        {
            DeleteMolfileV3000Info(mfdata->ctab.v3000);
        }

        inchi_free(mfdata);
        mfdata = NULL;
    }

    return mfdata;
}
    */
    // END INCHI C FUNCTION: FreeMolfileData
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: FreeMolfileData
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; inchi_free is the active checked-free macro behavior.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: FreeMolfileData

    if mfdata.is_null() {
        return Ok(SourceMutPointer::null());
    }
    let value = heap
        .slice(mfdata.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if !value.ctab.atoms.is_null() {
        inchi_free(heap, value.ctab.atoms)?;
    }
    if !value.ctab.bonds.is_null() {
        inchi_free(heap, value.ctab.bonds)?;
    }
    if !value.ctab.coords.is_null() {
        inchi_free(heap, value.ctab.coords)?;
    }
    let mut sgroups = value.ctab.sgroups;
    MolFmtSgroups_Free(heap, Some(&mut sgroups))?;
    if !value.ctab.v3000.is_null() {
        let _ = DeleteMolfileV3000Info(heap, value.ctab.v3000)?;
    }
    inchi_free(heap, mfdata)?;
    Ok(SourceMutPointer::null())
}

fn source_fseek(
    heap: &mut SourceHeap,
    file: SourceMutPointer<FILE>,
    position: i64,
) -> Result<i32, SourceHeapError> {
    if file.is_null() || position < 0 {
        return Ok(1);
    }
    let file = heap
        .slice_mut(file)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    // INCHI source mol_fmt2.c417: fseek(infile, fPtrStart, SEEK_SET)
    // SUPPLEMENTAL SOURCE: musl/src/stdio/fseek.c25..31
    // 	if (f->seek(f, off, whence) < 0) return -1;
    //
    // 	/* If seek succeeded, file is seekable and we discard read buffer. */
    // 	f->rpos = f->rend = 0;
    // 	f->flags &= ~F_EOF;
    //
    // 	return 0;
    // INCHI✔️❗: This owned-byte backend is seekable, with no pending native
    // write buffer. Retain its error indicator; only EOF is cleared. Forced
    // read/write errors remain in their respective IO operations. O(1), no
    // allocation/clone/scan; INCHI-PERF-STICKY-ERROR-SEEK-TELL costs OPEN.
    file.position = position as u64;
    file.eof = false;
    Ok(0)
}

fn source_ftell(heap: &SourceHeap, file: SourceMutPointer<FILE>) -> Result<i64, SourceHeapError> {
    if file.is_null() {
        return Ok(-1);
    }
    let file = heap
        .slice(file.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    // INCHI source mol_fmt2.c422: ftell(infile)
    // SUPPLEMENTAL SOURCE: musl/src/stdio/ftell.c12..17
    // 	/* Adjust for data in buffer. */
    // 	if (f->rend)
    // 		pos += f->rpos - f->rend;
    // 	else if (f->wbase)
    // 		pos += f->wpos - f->wbase;
    // 	return pos;
    // INCHI✔️❗: Position belongs to the existing seekable owned-byte stream;
    // its retained IO error indicator does not invalidate that position.
    // O(1), no allocation/clone/scan; checked conversion and typed errors
    // unchanged. INCHI-PERF-STICKY-ERROR-SEEK-TELL whole costs remain OPEN.
    i64::try_from(file.position).map_err(|_| SourceHeapError::SourceIntegerOverflow)
}

fn source_fputs(
    heap: &mut SourceHeap,
    string: SourceMutPointer<i8>,
    file: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    if file.is_null() {
        return Ok(-1);
    }
    let input = heap.slice(string.as_const())?;
    let length = input
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    // Preserve source-byte admission before FILE type/offset/extent errors.
    // A valid byte input and FILE record cannot share one typed allocation.
    let state = heap
        .slice_mut(file)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if length == 0 {
        return Ok(0);
    }
    if !source_fputs_can_write(state) {
        return Ok(-1);
    }
    heap.with_slice_mut_and_heap(file, |file, heap| {
        let input = &heap.slice(string.as_const())?[..length];
        source_fputs_content(&mut file[0], input)
    })
}

fn source_fputs_bytes(
    heap: &mut SourceHeap,
    string: &[i8],
    file: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    if file.is_null() {
        return Ok(-1);
    }
    let length = string
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let file = heap
        .slice_mut(file)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if length == 0 {
        return Ok(0);
    }
    if !source_fputs_can_write(file) {
        return Ok(-1);
    }
    source_fputs_content(file, &string[..length])
}

fn source_fputs_can_write(file: &mut FILE) -> bool {
    // INCHI source mol_fmt2.c449:             fputs(line, outfile);
    // SUPPLEMENTAL SOURCE: musl/src/stdio/__stdio_write.c21..23
    // 		if (cnt < 0) {
    // 			f->wpos = f->wbase = f->wend = 0;
    // 			f->flags |= F_ERR;
    // INCHI✔️❗: A forced backend failure sets the independent sticky flag;
    // earlier F_ERR alone never suppresses a subsequent successful write.
    // Empty fputs returns0 before this backend attempt (fwrite count0).
    // O(1), no allocation/scan/clone; FILE-ERROR-STATE costs remain OPEN.
    if file.error {
        file.error_indicator = true;
    }
    !file.error
}

fn source_fputs_content(file: &mut FILE, input: &[i8]) -> Result<i32, SourceHeapError> {
    // INCHI source mol_fmt2.c444:             fputs(line, outfile);
    // BEGIN INCHI SUPPLEMENTAL LIBC: third_party/musl/src/stdio/fputs.c
    // int fputs(const char *restrict s, FILE *restrict f)
    // {
    // 	size_t l = strlen(s);
    // 	return (fwrite(s, 1, l, f)==l) - 1;
    // }
    // END INCHI SUPPLEMENTAL LIBC: fputs
    // INCHI✔️❗: Borrow bytes without a Vec snapshot; signed/u8 projection is
    // bit-preserving. The caller distinguishes forced failure from the
    // independent sticky indicator. O(length), O(1) scratch; output growth
    // is unchanged. No configured-libc certification is claimed. Safe
    // arena/whole-call costs remain OPEN INCHI-PERF-MOL-COPY-STACK.
    let input = bytemuck::cast_slice::<i8, u8>(input);
    let position =
        usize::try_from(file.position).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if file.bytes.len() < position {
        file.bytes.resize(position, 0);
    }
    let overwrite = input.len().min(file.bytes.len().saturating_sub(position));
    file.bytes[position..position + overwrite].copy_from_slice(&input[..overwrite]);
    file.bytes.extend_from_slice(&input[overwrite..]);
    file.position = file
        .position
        .checked_add(
            u64::try_from(input.len()).map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
        )
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn MolfileSaveCopy(
    heap: &mut SourceHeap,
    input: Option<&mut INCHI_IOSTREAM>,
    start: i64,
    end: i64,
    output: SourceMutPointer<FILE>,
    number: i64,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:385 MolfileSaveCopy
    // INCHI✔️❌: int MolfileSaveCopy(INCHI_IOSTREAM *inp_file,
    // INCHI✔️❌:                     long fPtrStart,
    // INCHI✔️❌:                     long fPtrEnd,
    // INCHI✔️❌:                     FILE *outfile,
    // INCHI✔️❌:                     long num)
    // INCHI✔️❌: {
    // INCHI✔️❌:     char line[MOL_FMT_INPLINELEN], *p;
    // INCHI✔️❌:     long fPtr;
    // INCHI✔️❌:     int ret = 1;
    // INCHI✔️❌:     char szNumber[32];
    // INCHI✔️❌:
    // INCHI✔️❌:     if (inp_file->type == INCHI_IOS_TYPE_FILE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         FILE *infile = inp_file->f;
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!infile)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!outfile)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (fPtrStart < 0L && fPtrEnd <= fPtrStart)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (0 != fseek(infile, fPtrStart, SEEK_SET))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         while (fPtrEnd > (fPtr = ftell(infile)) && fPtr >= 0L
    // INCHI✔️❌:                 && inchi_fgetsLf(line, sizeof(line) - 1, inp_file))
    // INCHI✔️❌:         {
    // INCHI✔️❌:
    // INCHI✔️❌:             line[sizeof(line) - 1] = '\0'; /*  unnecessary extra precaution */
    // INCHI✔️❌:
    // INCHI✔️❌:             if (fPtr == fPtrStart && num)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 int len;
    // INCHI✔️❌:                 lrtrim(line, &len);
    // INCHI✔️❌:                 len = sprintf(szNumber, "#%ld%s", num, len ? "/" : "");
    // INCHI✔️❌:                 mystrncpy(line + len, line, sizeof(line) - len - 1);
    // INCHI✔️❌:                 memcpy(line, szNumber, len);
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (!strchr(line, '\n'))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 p = line + strlen(line);
    // INCHI✔️❌:                 p[0] = '\n';
    // INCHI✔️❌:                 p[1] = '\0';
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             fputs(line, outfile);
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         ret = fseek(infile, fPtrEnd, SEEK_SET);
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else if (inp_file->type == INCHI_IOS_TYPE_STRING)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MolfileSaveCopy

    // INCHI✔️❗: Complete fixed-array path shares canonical gets/trim/strncpy
    // kernels. Decimal prefix uses bounded stack bytes; overlap is memmove
    // followed by source zero padding. No heap scratch or output snapshot.
    // O(reached line bytes), fixed stack; whole-call safe costs remain OPEN
    // INCHI-PERF-MOL-COPY-STACK. Independent FILE backends remain guarded.
    let input = input.ok_or(SourceHeapError::NullPointer)?;
    let mut result = 1;
    if input.type_ == INCHI_IOS_TYPE_FILE as i32 {
        if input.f.is_null() || output.is_null() || (start < 0 && end <= start) {
            return Ok(1);
        }
        if source_fseek(heap, input.f, start)? != 0 {
            return Ok(1);
        }

        // Native line[204]/szNumber[32] have stack lifetime. No scratch
        // allocation/registration/failure; bounded safe zeroing cost OPEN.
        let mut line = [0_i8; MOL_FMT_INPLINELEN as usize];
        loop {
            let position = source_ftell(heap, input.f)?;
            if end <= position || position < 0 {
                break;
            }
            if !inchi_fgetsLf_bytes(heap, &mut line, MOL_FMT_INPLINELEN as i32 - 1, Some(input))? {
                break;
            }
            line[MOL_FMT_INPLINELEN as usize - 1] = 0;
            if position == start && number != 0 {
                let trimmed_length = lrtrim_bytes(&mut line)?;
                let mut digits = [0_u8; 21];
                let number_bytes = source_decimal_i64(number, Some(b'#'), &mut digits);
                let mut prefix = [0_i8; 32];
                let mut prefix_length = number_bytes.len();
                prefix[..prefix_length].copy_from_slice(bytemuck::cast_slice(number_bytes));
                if trimmed_length != 0 {
                    prefix[prefix_length] = b'/' as i8;
                    prefix_length += 1;
                }
                mystrncpy_within(
                    &mut line,
                    prefix_length,
                    0,
                    u32::try_from(MOL_FMT_INPLINELEN as usize - prefix_length - 1)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                )?;
                line[..prefix_length].copy_from_slice(&prefix[..prefix_length]);
            }
            // strchr determines LF/NUL once; strlen is reached only if no LF.
            let length = line
                .iter()
                .position(|byte| *byte == 0 || *byte == b'\n' as i8)
                .ok_or(SourceHeapError::MissingNulTerminator)?;
            if line[length] != b'\n' as i8 {
                if length + 1 >= line.len() {
                    return Err(SourceHeapError::PointerOutOfBounds);
                }
                line[length] = b'\n' as i8;
                line[length + 1] = 0;
            }
            source_fputs_bytes(heap, &line, output)?;
        }
        result = source_fseek(heap, input.f, end)?;
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{
        INCHI_IOS_TYPE_STRING, INT_ARRAY, MOL_COORD, MOL_FMT_ATOM, MOL_FMT_BOND, MOL_FMT_CTAB,
        MOL_FMT_SGROUP, MOL_FMT_SGROUPS, MOL_FMT_v3000, SourceFile,
    };

    #[test]
    fn source_port__mol_fmt2__freemolfiledata__line_664() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            FreeMolfileData(&mut heap, SourceMutPointer::null()),
            Ok(SourceMutPointer::null())
        );

        let atoms = heap
            .allocate_model_storage(vec![MOL_FMT_ATOM::default(); 2])
            .unwrap();
        let bonds = heap
            .allocate_model_storage(vec![MOL_FMT_BOND::default()])
            .unwrap();
        let coords = heap
            .allocate_model_storage(vec![MOL_COORD::default(); 2])
            .unwrap();
        let atom_list = heap.allocate_model_storage(vec![1_i32, 2]).unwrap();
        let bond_list = heap.allocate_model_storage(vec![1_i32]).unwrap();
        let group = heap
            .allocate_model_storage(vec![MOL_FMT_SGROUP {
                alist: INT_ARRAY {
                    item: atom_list,
                    allocated: 2,
                    used: 2,
                    increment: 2,
                },
                blist: INT_ARRAY {
                    item: bond_list,
                    allocated: 1,
                    used: 1,
                    increment: 1,
                },
                ..MOL_FMT_SGROUP::default()
            }])
            .unwrap();
        let group_array = heap.allocate_model_storage(vec![group]).unwrap();
        let v3000 = heap
            .allocate_model_storage(vec![MOL_FMT_v3000::default()])
            .unwrap();
        let data = heap
            .allocate_model_storage(vec![MOL_FMT_DATA {
                ctab: MOL_FMT_CTAB {
                    n_atoms: 2,
                    n_bonds: 1,
                    atoms,
                    bonds,
                    coords,
                    sgroups: MOL_FMT_SGROUPS {
                        group: group_array,
                        allocated: 1,
                        used: 1,
                        increment: 1,
                    },
                    v3000,
                    ..MOL_FMT_CTAB::default()
                },
                ..MOL_FMT_DATA::default()
            }])
            .unwrap();
        assert_eq!(
            FreeMolfileData(&mut heap, data),
            Ok(SourceMutPointer::null())
        );
        assert_eq!(
            heap.slice(atoms.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(bonds.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(coords.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(atom_list.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(bond_list.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(group.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(group_array.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(v3000.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(data.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__mol_fmt2__molfilestrnread__line_65() {
        fn bytes(values: &[u8]) -> Vec<i8> {
            values.iter().map(|byte| *byte as i8).collect()
        }

        let mut heap = SourceHeap::default();
        let source = heap.allocate(bytes(b"unused\0")).unwrap();
        let destination = heap.allocate(bytes(b"ABCD")).unwrap();
        let mut first_space = SourceMutPointer::null();
        assert_eq!(
            MolfileStrnread(
                &mut heap,
                destination.offset(1).unwrap(),
                source.as_const(),
                0,
                &mut first_space,
            ),
            Ok(0)
        );
        assert_eq!(first_space, destination.offset(1).unwrap());
        assert_eq!(
            heap.slice(destination.as_const()).unwrap(),
            &bytes(b"A\0CD")
        );

        let short_source = heap.allocate(bytes(b"xy\0ignored")).unwrap();
        let short_destination = heap.allocate(vec![0x55_i8; 8]).unwrap();
        assert_eq!(
            MolfileStrnread(
                &mut heap,
                short_destination,
                short_source.as_const(),
                5,
                &mut first_space,
            ),
            Ok(2)
        );
        assert_eq!(first_space, short_destination.offset(2).unwrap());
        assert_eq!(
            heap.slice(short_destination.as_const()).unwrap(),
            &[b'x' as i8, b'y' as i8, 0, 0, 0, 0, 0x55, 0x55]
        );

        let trailing_source = heap.allocate(bytes(b"ab \t\r\n\x0b\x0cTAIL")).unwrap();
        let trailing_destination = heap.allocate(vec![0x33_i8; 12]).unwrap();
        assert_eq!(
            MolfileStrnread(
                &mut heap,
                trailing_destination,
                trailing_source.as_const(),
                8,
                &mut first_space,
            ),
            Ok(8)
        );
        assert_eq!(first_space, trailing_destination.offset(2).unwrap());
        assert_eq!(
            &heap.slice(trailing_destination.as_const()).unwrap()[..10],
            &[
                b'a' as i8,
                b'b' as i8,
                b' ' as i8,
                b'\t' as i8,
                b'\r' as i8,
                b'\n' as i8,
                0x0b,
                0x0c,
                0,
                0x33
            ]
        );

        let interior_source = heap.allocate(bytes(b"a b\0")).unwrap();
        let interior_destination = heap.allocate(vec![0x22_i8; 5]).unwrap();
        assert_eq!(
            MolfileStrnread(
                &mut heap,
                interior_destination,
                interior_source.as_const(),
                3,
                &mut first_space,
            ),
            Ok(3)
        );
        assert_eq!(first_space, interior_destination.offset(3).unwrap());

        let all_space_source = heap.allocate(bytes(b" \t\n\r\x0b\x0c")).unwrap();
        let all_space_destination = heap.allocate(vec![0_i8; 7]).unwrap();
        assert_eq!(
            MolfileStrnread(
                &mut heap,
                all_space_destination,
                all_space_source.as_const(),
                6,
                &mut first_space,
            ),
            Ok(6)
        );
        assert_eq!(first_space, all_space_destination);

        let high_source = heap.allocate(vec![-1_i8, b' ' as i8, 0]).unwrap();
        let high_destination = heap.allocate(vec![0_i8; 4]).unwrap();
        assert_eq!(
            MolfileStrnread(
                &mut heap,
                high_destination,
                high_source.as_const(),
                2,
                &mut first_space,
            ),
            Ok(2)
        );
        assert_eq!(first_space, high_destination.offset(1).unwrap());
        assert_eq!(
            heap.slice(high_destination.as_const()).unwrap(),
            &[-1_i8, b' ' as i8, 0, 0]
        );

        assert_eq!(
            MolfileStrnread(
                &mut heap,
                high_destination,
                high_source.as_const(),
                -1,
                &mut first_space,
            ),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
    }

    #[test]
    fn source_port__mol_fmt2__molfilereadfield__line_105() {
        fn line(heap: &mut SourceHeap, text: &str) -> SourceMutPointer<i8> {
            heap.allocate(
                text.as_bytes()
                    .iter()
                    .map(|byte| *byte as i8)
                    .chain(std::iter::once(0))
                    .collect(),
            )
            .unwrap()
        }

        let mut heap = SourceHeap::default();
        let base = line(&mut heap, "  AB  X");
        let mut cursor = base;
        let mut string = [0x55_i8; 8];
        assert_eq!(
            MolfileReadField(
                &mut heap,
                MolfileFieldData::String(&mut string),
                6,
                i32::from(MOL_FMT_STRING_DATA),
                &mut cursor,
            ),
            Ok(2)
        );
        assert_eq!(cursor, base.offset(6).unwrap());
        assert_eq!(&string[..5], &[b'A' as i8, b'B' as i8, 0, b' ' as i8, 0]);

        let fixed = line(&mut heap, " 12 12x 128");
        cursor = fixed;
        let mut char_value = 9_i8;
        assert_eq!(
            MolfileReadField(
                &mut heap,
                MolfileFieldData::Char(&mut char_value),
                4,
                i32::from(MOL_FMT_CHAR_INT_DATA),
                &mut cursor
            ),
            Ok(1)
        );
        assert_eq!((char_value, cursor), (12, fixed.offset(4).unwrap()));
        assert_eq!(
            MolfileReadField(
                &mut heap,
                MolfileFieldData::Char(&mut char_value),
                3,
                i32::from(MOL_FMT_CHAR_INT_DATA),
                &mut cursor
            ),
            Ok(-1)
        );
        assert_eq!((char_value, cursor), (12, fixed.offset(7).unwrap()));
        cursor = fixed.offset(8).unwrap();
        assert_eq!(
            MolfileReadField(
                &mut heap,
                MolfileFieldData::Char(&mut char_value),
                3,
                i32::from(MOL_FMT_CHAR_INT_DATA),
                &mut cursor
            ),
            Ok(-1)
        );
        assert_eq!(char_value, 0);

        let free = line(&mut heap, "-32768x -9223372036854775808!");
        cursor = free;
        let mut short_value = 0_i16;
        assert_eq!(
            MolfileReadField(
                &mut heap,
                MolfileFieldData::Short(&mut short_value),
                0,
                i32::from(MOL_FMT_SHORT_INT_DATA),
                &mut cursor
            ),
            Ok(1)
        );
        assert_eq!((short_value, cursor), (-32768, free.offset(6).unwrap()));
        cursor = free.offset(8).unwrap();
        let mut long_value = 7_i64;
        assert_eq!(
            MolfileReadField(
                &mut heap,
                MolfileFieldData::Long(&mut long_value),
                -1,
                i32::from(MOL_FMT_LONG_INT_DATA),
                &mut cursor
            ),
            Ok(-1)
        );
        assert_eq!(long_value, 0);
        assert_eq!(cursor, free.offset(28).unwrap());

        let too_long = line(&mut heap, "9");
        cursor = too_long;
        short_value = 4;
        assert_eq!(
            MolfileReadField(
                &mut heap,
                MolfileFieldData::Short(&mut short_value),
                33,
                i32::from(MOL_FMT_SHORT_INT_DATA),
                &mut cursor
            ),
            Ok(-1)
        );
        assert_eq!((short_value, cursor), (0, too_long));

        let floating = line(&mut heap, "nan inf 1e-50 1e40 nope");
        cursor = floating;
        let mut double_value = 0.0;
        assert_eq!(
            MolfileReadField(
                &mut heap,
                MolfileFieldData::Double(&mut double_value),
                0,
                i32::from(MOL_FMT_DOUBLE_DATA),
                &mut cursor
            ),
            Ok(1)
        );
        assert!(double_value.is_nan());
        cursor = floating.offset(4).unwrap();
        assert_eq!(
            MolfileReadField(
                &mut heap,
                MolfileFieldData::Double(&mut double_value),
                0,
                i32::from(MOL_FMT_DOUBLE_DATA),
                &mut cursor
            ),
            Ok(-1)
        );
        assert_eq!(double_value, 0.0);
        cursor = floating.offset(8).unwrap();
        let mut float_value = 3.0_f32;
        assert_eq!(
            MolfileReadField(
                &mut heap,
                MolfileFieldData::Float(&mut float_value),
                0,
                i32::from(MOL_FMT_FLOAT_DATA),
                &mut cursor
            ),
            Ok(1)
        );
        assert_eq!(float_value, 0.0);
        cursor = floating.offset(14).unwrap();
        assert_eq!(
            MolfileReadField(
                &mut heap,
                MolfileFieldData::Float(&mut float_value),
                0,
                i32::from(MOL_FMT_FLOAT_DATA),
                &mut cursor
            ),
            Ok(-1)
        );
        assert_eq!(float_value, 0.0);
        cursor = floating.offset(19).unwrap();
        assert_eq!(
            MolfileReadField(
                &mut heap,
                MolfileFieldData::Double(&mut double_value),
                0,
                i32::from(MOL_FMT_DOUBLE_DATA),
                &mut cursor
            ),
            Ok(0)
        );
        assert_eq!((double_value, cursor), (0.0, floating.offset(19).unwrap()));

        let jump = line(&mut heap, "ab");
        cursor = jump;
        assert_eq!(
            MolfileReadField(
                &mut heap,
                MolfileFieldData::None,
                8,
                i32::from(MOL_FMT_JUMP_TO_RIGHT),
                &mut cursor
            ),
            Ok(2)
        );
        assert_eq!(cursor, jump.offset(2).unwrap());
        cursor = jump;
        assert_eq!(
            MolfileReadField(
                &mut heap,
                MolfileFieldData::None,
                -1,
                i32::from(MOL_FMT_JUMP_TO_RIGHT),
                &mut cursor
            ),
            Ok(0)
        );
        assert_eq!(
            MolfileReadField(&mut heap, MolfileFieldData::None, 1, -7, &mut cursor),
            Ok(-1)
        );
        assert_eq!(cursor, jump);
    }

    #[test]
    fn source_port__mol_fmt2__molfileextractstrucnum__line_333() {
        fn fixed<const N: usize>(text: &str) -> [i8; N] {
            assert!(text.len() < N);
            let mut result = [0_i8; N];
            for (target, source) in result.iter_mut().zip(text.bytes()) {
                *target = source as i8;
            }
            result
        }

        fn header(name: &str, line2: &str) -> MOL_FMT_HEADER_BLOCK {
            MOL_FMT_HEADER_BLOCK {
                molname: fixed(name),
                line2: fixed(line2),
                ..MOL_FMT_HEADER_BLOCK::default()
            }
        }

        let mut heap = SourceHeap::default();
        assert_eq!(MolfileExtractStrucNum(&mut heap, None), Ok(0));

        for (name, line2, expected) in [
            ("Structure #22", "InChI v1.07 SDfile Output", 22_i64),
            ("sTrUcTuRe #-22", "iNcHi SDfile Output", -22),
            ("Structure #  +22", "InChI--SDfile Output--", 22),
            ("Structure #0", "InChI SDfile Output", 0),
            (
                "Structure #9223372036854775807",
                "InChI SDfile Output",
                i64::MAX,
            ),
            (
                "Structure #-9223372036854775808",
                "InChI SDfile Output",
                i64::MIN,
            ),
            (
                "Structure #9223372036854775808",
                "InChI SDfile Output",
                i64::MAX,
            ),
            (
                "Structure #-9223372036854775809",
                "InChI SDfile Output",
                i64::MIN,
            ),
            ("structure #17x", "InChI SDfile Output", 0),
            ("Structure #17 ", "InChI SDfile Output", 0),
            ("Structure #", "InChI SDfile Output", 0),
            ("Structure 17", "InChI SDfile Output", 0),
            ("Xtructure #17", "InChI SDfile Output", 0),
            ("Structure #17", "XnChI SDfile Output", 0),
            ("Structure #17", "InChI sdfile Output", 0),
            ("Structure #17", "InChI SDfile  Output", 0),
            ("Structure #17", "InChISDfile Output", 17),
        ] {
            let input = header(name, line2);
            assert_eq!(
                MolfileExtractStrucNum(&mut heap, Some(&input)),
                Ok(expected),
                "name={name:?}, line2={line2:?}"
            );
            assert_eq!(heap.live_allocations_of::<i8>(), 0);
        }
        assert_eq!(heap.source_errno(), 34);
    }

    #[test]
    fn source_port__mol_fmt2__molfilehasnochemstruc__line_362() {
        assert_eq!(MolfileHasNoChemStruc(None), 1);

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate(vec![MOL_FMT_ATOM::default()]).unwrap();
        let bonds = heap.allocate(vec![MOL_FMT_BOND::default()]).unwrap();

        for atom_count in [i32::MIN, -1, 0] {
            let input = MOL_FMT_DATA {
                ctab: MOL_FMT_CTAB {
                    n_atoms: atom_count,
                    atoms,
                    ..MOL_FMT_CTAB::default()
                },
                ..MOL_FMT_DATA::default()
            };
            assert_eq!(MolfileHasNoChemStruc(Some(&input)), 1);
        }

        let no_atoms = MOL_FMT_DATA {
            ctab: MOL_FMT_CTAB {
                n_atoms: i32::MAX,
                n_bonds: i32::MAX,
                bonds,
                ..MOL_FMT_CTAB::default()
            },
            ..MOL_FMT_DATA::default()
        };
        assert_eq!(MolfileHasNoChemStruc(Some(&no_atoms)), 1);

        for bond_count in [1, i32::MAX] {
            let input = MOL_FMT_DATA {
                ctab: MOL_FMT_CTAB {
                    n_atoms: 1,
                    n_bonds: bond_count,
                    atoms,
                    ..MOL_FMT_CTAB::default()
                },
                ..MOL_FMT_DATA::default()
            };
            assert_eq!(MolfileHasNoChemStruc(Some(&input)), 1);
        }

        for bond_count in [i32::MIN, -1, 0] {
            let input = MOL_FMT_DATA {
                ctab: MOL_FMT_CTAB {
                    n_atoms: 1,
                    n_bonds: bond_count,
                    atoms,
                    ..MOL_FMT_CTAB::default()
                },
                ..MOL_FMT_DATA::default()
            };
            assert_eq!(MolfileHasNoChemStruc(Some(&input)), 0);
        }

        let complete = MOL_FMT_DATA {
            ctab: MOL_FMT_CTAB {
                n_atoms: 1,
                n_bonds: 1,
                atoms,
                bonds,
                ..MOL_FMT_CTAB::default()
            },
            ..MOL_FMT_DATA::default()
        };
        assert_eq!(MolfileHasNoChemStruc(Some(&complete)), 0);
    }

    #[test]
    fn source_port__mol_fmt2__molfilegetxyzdimandnormfactors__line_474() {
        fn fixture(
            heap: &mut SourceHeap,
            coordinates: &[(f64, f64, f64)],
            bonds: &[(i16, i16)],
        ) -> MOL_FMT_DATA {
            let atoms = heap
                .allocate(
                    coordinates
                        .iter()
                        .map(|&(fx, fy, fz)| MOL_FMT_ATOM {
                            fx,
                            fy,
                            fz,
                            ..MOL_FMT_ATOM::default()
                        })
                        .collect(),
                )
                .unwrap();
            let bond_pointer = if bonds.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate(
                    bonds
                        .iter()
                        .map(|&(atnum1, atnum2)| MOL_FMT_BOND {
                            atnum1,
                            atnum2,
                            ..MOL_FMT_BOND::default()
                        })
                        .collect(),
                )
                .unwrap()
            };
            MOL_FMT_DATA {
                ctab: MOL_FMT_CTAB {
                    n_atoms: coordinates.len() as i32,
                    n_bonds: bonds.len() as i32,
                    atoms,
                    bonds: bond_pointer,
                    ..MOL_FMT_CTAB::default()
                },
                ..MOL_FMT_DATA::default()
            }
        }

        fn evaluate(
            heap: &SourceHeap,
            input: Option<&MOL_FMT_DATA>,
            find_norm_factors: i32,
            err: &mut i32,
            errors: Option<&mut [i8]>,
        ) -> (i32, f64, f64, f64, f64, f64, f64, f64) {
            let (mut x0, mut y0, mut z0) = (91.0, 92.0, 93.0);
            let (mut xmin, mut ymin, mut zmin) = (94.0, 95.0, 96.0);
            let mut scaler = 97.0;
            let dimensions = MolfileGetXYZDimAndNormFactors(
                heap,
                input,
                find_norm_factors,
                &mut x0,
                &mut y0,
                &mut z0,
                &mut xmin,
                &mut ymin,
                &mut zmin,
                &mut scaler,
                err,
                errors,
            )
            .unwrap();
            (dimensions, x0, y0, z0, xmin, ymin, zmin, scaler)
        }

        let empty_heap = SourceHeap::default();
        let mut err = 37;
        assert_eq!(
            evaluate(&empty_heap, None, 1, &mut err, None),
            (0, 1.0e32, 1.0e32, 1.0e32, 0.0, 0.0, 0.0, 1.0)
        );
        assert_eq!(err, 37);

        let mut heap = SourceHeap::default();
        let point = fixture(&mut heap, &[(5.0, -7.0, 11.0)], &[]);
        let mut err = 0;
        assert_eq!(
            evaluate(&heap, Some(&point), 1, &mut err, None),
            (0, 5.0, -7.0, 11.0, 0.0, 0.0, 0.0, 0.0)
        );

        let planar = fixture(&mut heap, &[(0.0, 0.0, 0.0), (1.0, 2.0, 0.0)], &[]);
        assert_eq!(
            evaluate(&heap, Some(&planar), 0, &mut err, None),
            (2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)
        );
        assert_eq!(evaluate(&heap, Some(&planar), 1, &mut err, None).7, 128.0);

        let spatial = fixture(&mut heap, &[(0.0, 0.0, 0.0), (3.0, 4.0, 12.0)], &[(1, 2)]);
        let spatial_result = evaluate(&heap, Some(&spatial), 1, &mut err, None);
        assert_eq!(spatial_result.0, 3);
        assert_eq!(spatial_result.7, 20.0 / 13.0);

        let short = fixture(
            &mut heap,
            &[(0.0, 0.0, 0.0), (1000.0, 0.0, 0.0), (0.1, 0.0, 0.0)],
            &[(1, 3)],
        );
        assert_eq!(evaluate(&heap, Some(&short), 1, &mut err, None).7, 100.0);

        let zero_length = fixture(&mut heap, &[(2.0, 3.0, 4.0), (2.0, 3.0, 4.0)], &[(1, 2)]);
        assert_eq!(
            evaluate(&heap, Some(&zero_length), 1, &mut err, None).7,
            1.0
        );

        let nan_after_finite = fixture(
            &mut heap,
            &[(1.0, 2.0, 3.0), (f64::NAN, f64::NAN, f64::NAN)],
            &[],
        );
        let nan_result = evaluate(&heap, Some(&nan_after_finite), 1, &mut err, None);
        assert_eq!(
            (nan_result.0, nan_result.1, nan_result.2, nan_result.3),
            (0, 1.0, 2.0, 3.0)
        );

        let invalid = fixture(&mut heap, &[(0.0, 0.0, 0.0)], &[(0, 1), (1, 1)]);
        let mut errors = [0_i8; 256];
        let invalid_result = evaluate(&heap, Some(&invalid), 1, &mut err, Some(&mut errors));
        assert_eq!(invalid_result.0, 0);
        assert_eq!(err, 1);
        let error_length = errors.iter().position(|byte| *byte == 0).unwrap();
        assert_eq!(
            errors[..error_length]
                .iter()
                .map(|byte| *byte as u8)
                .collect::<Vec<_>>(),
            b"Bond to nonexistent atom"
        );
    }

    #[test]
    fn source_port__mol_fmt2__molfilesavecopy__line_385() {
        let mut heap = SourceHeap::default();
        let input_file = heap
            .allocate(vec![SourceFile {
                bytes: b"skip\n  name  \nline without lf".to_vec(),
                ..SourceFile::default()
            }])
            .unwrap();
        let output_file = heap
            .allocate(vec![SourceFile {
                bytes: b"old".to_vec(),
                position: 1,
                ..SourceFile::default()
            }])
            .unwrap();
        let mut stream = INCHI_IOSTREAM {
            f: input_file,
            type_: INCHI_IOS_TYPE_FILE as i32,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(
            MolfileSaveCopy(&mut heap, Some(&mut stream), 5, 30, output_file, 7),
            Ok(0)
        );
        let output = &heap.slice(output_file.as_const()).unwrap()[0];
        assert_eq!(output.bytes, b"o#7/name\nline without lf\n".to_vec());
        assert_eq!(output.position, 25);
        assert_eq!(heap.slice(input_file.as_const()).unwrap()[0].position, 30);

        heap.slice_mut(input_file).unwrap()[0].position = 20;
        heap.slice_mut(output_file).unwrap()[0] = SourceFile::default();
        assert_eq!(
            MolfileSaveCopy(&mut heap, Some(&mut stream), 0, 5, output_file, 0),
            Ok(0)
        );
        assert_eq!(
            heap.slice(output_file.as_const()).unwrap()[0].bytes,
            b"skip\n".to_vec()
        );

        assert_eq!(
            MolfileSaveCopy(&mut heap, Some(&mut stream), -2, -3, output_file, 0),
            Ok(1)
        );
        let mut string_stream = INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(
            MolfileSaveCopy(&mut heap, Some(&mut string_stream), 0, 1, output_file, 0),
            Ok(1)
        );
        let mut null_file_stream = INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_FILE as i32,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(
            MolfileSaveCopy(&mut heap, Some(&mut null_file_stream), 0, 1, output_file, 0),
            Ok(1)
        );
        assert_eq!(
            MolfileSaveCopy(
                &mut heap,
                Some(&mut stream),
                0,
                1,
                SourceMutPointer::null(),
                0
            ),
            Ok(1)
        );
        assert_eq!(
            MolfileSaveCopy(&mut heap, None, 0, 1, output_file, 0),
            Err(SourceHeapError::NullPointer)
        );
    }
}
