use crate::source::base::util::{inchi_calloc, inchi_free, lrtrim};
use crate::source_types::{
    FILE, INCHI_IOS_STRING, INCHI_IOS_TYPE_FILE, INCHI_IOS_TYPE_STRING, INCHI_IOSTREAM,
    INCHI_STRBUF_INITIAL_SIZE, INCHI_STRBUF_SIZE_INCREMENT, SourceConstPointer,
    SourceFormatArgument, SourceHeap, SourceHeapError, SourceMutPointer, SourceVaList,
    local_ichi_io::INCHI_ADD_STR_LEN,
};

const SOURCE_EOF: i32 = -1;

pub(crate) fn source_decimal_i64(value: i64, prefix: Option<u8>, buffer: &mut [u8; 21]) -> &[u8] {
    // INCHI✔️❗: native unpadded signed-decimal source anchors:
    // readinch.c:166 / 306
    // sprintf(szTemp, "#%d", i0D + 1);
    // readinch.c:432
    // sprintf( szErrCode, "%d", k );
    // mol_fmt4.c:1225 / 1236
    // len += sprintf(str_m + len, "%d", num_H);
    // len += sprintf(str_m + len, "%d", j);
    // The configured C int/long arguments are widened by callers. Maximum
    // signed magnitude has 19 decimal digits; sign and optional '#' fit in21.
    // unsigned_abs also handles i64::MIN. Only the initialized suffix escapes.
    // O(digits) time, fixed21byte scratch, no String/allocation/UTF8 scan;
    // checked indexing/division/whole-call cost DECIMAL-MESSAGES remains OPEN.
    let mut magnitude = value.unsigned_abs();
    let mut start = buffer.len();
    loop {
        start -= 1;
        buffer[start] = b'0' + (magnitude % 10) as u8;
        magnitude /= 10;
        if magnitude == 0 {
            break;
        }
    }
    if value < 0 {
        start -= 1;
        buffer[start] = b'-';
    }
    if let Some(prefix) = prefix {
        start -= 1;
        buffer[start] = prefix;
    }
    &buffer[start..]
}

pub(crate) fn inchi_strbuf_init(
    heap: &mut SourceHeap,
    buf: &mut INCHI_IOS_STRING,
    mut start_size: i32,
    mut incr_size: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1370 inchi_strbuf_init
    // INCHI✔❌: int inchi_strbuf_init(INCHI_IOS_STRING* buf, int start_size, int incr_size)
    // INCHI✔❌: {
    // INCHI✔❌:     char* new_str = NULL;
    // INCHI✔❌:     memset(buf, 0, sizeof(*buf)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:
    // INCHI✔❌:     if (start_size <= 0)
    // INCHI✔❌:     {
    // INCHI✔❌:         start_size = INCHI_STRBUF_INITIAL_SIZE;
    // INCHI✔❌:     }
    // INCHI✔❌:     if (incr_size <= 0)
    // INCHI✔❌:     {
    // INCHI✔❌:         incr_size = INCHI_STRBUF_SIZE_INCREMENT;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     new_str = (char*)inchi_calloc(start_size, sizeof(char));
    // INCHI✔❌:
    // INCHI✔❌:     if (!new_str)
    // INCHI✔❌:     {
    // INCHI✔❌:         return -1;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     buf->pStr = new_str;
    // INCHI✔❌:     buf->nAllocatedLength = start_size;
    // INCHI✔❌:     buf->nPtr = incr_size;
    // INCHI✔❌:
    // INCHI✔❌:     return start_size;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: inchi_strbuf_init

    *buf = INCHI_IOS_STRING::default();
    if start_size <= 0 {
        start_size = INCHI_STRBUF_INITIAL_SIZE as i32;
    }
    if incr_size <= 0 {
        incr_size = INCHI_STRBUF_SIZE_INCREMENT as i32;
    }
    let new_string = match inchi_calloc::<i8>(heap, start_size as u64, 1) {
        Ok(pointer) => pointer,
        Err(_) => return Ok(-1),
    };
    buf.pStr = new_string;
    buf.nAllocatedLength = start_size;
    buf.nPtr = incr_size;
    Ok(start_size)
}

pub(crate) fn inchi_strbuf_reset(
    heap: &mut SourceHeap,
    buf: Option<&mut INCHI_IOS_STRING>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1403 inchi_strbuf_reset
    // INCHI✔❌: void inchi_strbuf_reset(INCHI_IOS_STRING* buf)
    // INCHI✔❌: {
    // INCHI✔❌:     if (!buf)
    // INCHI✔❌:     {
    // INCHI✔❌:         return;
    // INCHI✔❌:     }
    // INCHI✔❌:     if (buf->pStr)
    // INCHI✔❌:     {
    // INCHI✔❌:         buf->pStr[0] = '\0';
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     buf->nUsedLength = buf->nPtr = 0;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: inchi_strbuf_reset

    let Some(buf) = buf else {
        return Ok(());
    };
    if !buf.pStr.is_null() {
        heap.slice_mut(buf.pStr)?[0] = 0;
    }
    buf.nPtr = 0;
    buf.nUsedLength = 0;
    Ok(())
}

pub(crate) fn inchi_strbuf_close(
    heap: &mut SourceHeap,
    buf: Option<&mut INCHI_IOS_STRING>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1422 inchi_strbuf_close
    // INCHI✔❌: void inchi_strbuf_close(INCHI_IOS_STRING* buf)
    // INCHI✔❌: {
    // INCHI✔❌:     if (!buf)
    // INCHI✔❌:     {
    // INCHI✔❌:         return;
    // INCHI✔❌:     }
    // INCHI✔❌:     if (buf->pStr)
    // INCHI✔❌:     {
    // INCHI✔❌:         inchi_free(buf->pStr);
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     memset(buf, 0, sizeof(*buf)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌: }
    // END INCHI C FUNCTION: inchi_strbuf_close

    let Some(buf) = buf else {
        return Ok(());
    };
    if !buf.pStr.is_null() {
        inchi_free(heap, buf.pStr)?;
    }
    *buf = INCHI_IOS_STRING::default();
    Ok(())
}

#[rustfmt::skip]
pub(crate) fn inchi_strbuf_create_copy(
    heap: &mut SourceHeap,
    destination: &mut INCHI_IOS_STRING,
    source: &INCHI_IOS_STRING,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1438 inchi_strbuf_create_copy
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int inchi_strbuf_create_copy(INCHI_IOS_STRING* buf2, INCHI_IOS_STRING* buf)
{
    char* new_str = NULL;

    new_str = (char*)inchi_calloc(buf->nAllocatedLength, sizeof(char));
    buf2->pStr = new_str;
    if (!new_str)
    {
        return -1;
    }
    buf2->nAllocatedLength = buf->nAllocatedLength;
    buf2->nUsedLength = buf->nUsedLength;
    buf2->nPtr = buf->nPtr;

    return 0;
}
    */
    // END INCHI C FUNCTION: inchi_strbuf_create_copy
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: inchi_strbuf_create_copy
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; inchi_calloc is the active calloc macro.
    // INCHI✔️❌: The source does not copy source bytes and does not free a pre-existing destination buffer.
    // INCHI✔️❌: SourceHeap allocation bookkeeping is materially slower than the source calloc assignment.
    // END INCHI ACTIVE MACRO CONFIGURATION: inchi_strbuf_create_copy

    let new_string = match inchi_calloc::<i8>(heap, source.nAllocatedLength as u64, 1) {
        Ok(pointer) => pointer,
        Err(
            SourceHeapError::AllocationFailed
            | SourceHeapError::AllocationSizeOverflow
            | SourceHeapError::AllocationElementCountOutOfRange,
        ) => SourceMutPointer::null(),
        Err(error) => return Err(error),
    };
    destination.pStr = new_string;
    if new_string.is_null() {
        return Ok(-1);
    }
    destination.nAllocatedLength = source.nAllocatedLength;
    destination.nUsedLength = source.nUsedLength;
    destination.nPtr = source.nPtr;
    Ok(0)
}

pub(crate) fn inchi_strbuf_update(
    heap: &mut SourceHeap,
    buf: Option<&mut INCHI_IOS_STRING>,
    new_addition_size: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1459 inchi_strbuf_update
    // INCHI✔️❌: int inchi_strbuf_update(INCHI_IOS_STRING* buf, int new_addition_size)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int requsted_len;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!buf)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (new_addition_size <= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return buf->nAllocatedLength;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     requsted_len = buf->nUsedLength + new_addition_size;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (requsted_len >= buf->nAllocatedLength)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Expand */
    // INCHI✔️❌:         int  nAddLength = inchi_max(buf->nPtr, new_addition_size);
    // INCHI✔️❌:         /* buf->nPtr stores size increment for this buffer */
    // INCHI✔️❌:         char* new_str =
    // INCHI✔️❌:             (char*)inchi_calloc((long long)buf->nAllocatedLength + (long long)nAddLength,
    // INCHI✔️❌:                 sizeof(new_str[0])); /* djb-rwth: cast operators added */
    // INCHI✔️❌:         if (!new_str)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return -1; /* failed */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (buf->pStr)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (buf->nUsedLength > 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 memcpy(new_str, buf->pStr, sizeof(new_str[0]) * buf->nUsedLength);
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_free(buf->pStr);
    // INCHI✔️❌:         }
    // INCHI✔️❌:         buf->pStr = new_str;
    // INCHI✔️❌:         buf->nAllocatedLength += nAddLength;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return buf->nAllocatedLength;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_strbuf_update

    let Some(buf) = buf else {
        return Ok(-1);
    };
    if new_addition_size <= 0 {
        return Ok(buf.nAllocatedLength);
    }
    let requested_length = buf.nUsedLength.wrapping_add(new_addition_size);
    if requested_length >= buf.nAllocatedLength {
        let addition_length = buf.nPtr.max(new_addition_size);
        let allocation_length = i64::from(buf.nAllocatedLength) + i64::from(addition_length);
        let new_string = match inchi_calloc::<i8>(heap, allocation_length as u64, 1) {
            Ok(pointer) => pointer,
            Err(_) => return Ok(-1),
        };
        if !buf.pStr.is_null() {
            if buf.nUsedLength > 0 {
                let used = usize::try_from(buf.nUsedLength)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                // Source memcpy reads the old live allocation directly. The
                // freshly allocated destination cannot alias that allocation.
                heap.slice(buf.pStr.as_const())?
                    .get(..used)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                heap.with_slice_mut_and_heap(new_string, |destination, heap| {
                    let source = heap.slice(buf.pStr.as_const())?;
                    destination
                        .get_mut(..used)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .copy_from_slice(&source[..used]);
                    Ok(())
                })?;
            }
            inchi_free(heap, buf.pStr)?;
        }
        buf.pStr = new_string;
        buf.nAllocatedLength = buf.nAllocatedLength.wrapping_add(addition_length);
    }
    Ok(buf.nAllocatedLength)
}

pub(crate) fn inchi_strbuf_printf(
    heap: &mut SourceHeap,
    buf: Option<&mut INCHI_IOS_STRING>,
    format: SourceConstPointer<i8>,
    arguments: &SourceVaList,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1507 inchi_strbuf_printf
    // INCHI✔️❌: int inchi_strbuf_printf(INCHI_IOS_STRING* buf, const char* lpszFormat, ...)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0, max_len;
    // INCHI✔️❌:     va_list argList;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!buf)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     my_va_start(argList, lpszFormat);
    // INCHI✔️❌:     max_len = GetMaxPrintfLength(lpszFormat, argList);
    // INCHI✔️❌:     va_end(argList);
    // INCHI✔️❌:     if (max_len < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_strbuf_update(buf, max_len);
    // INCHI✔️❌:
    // INCHI✔️❌:     my_va_start(argList, lpszFormat);
    // INCHI✔️❌:     ret = vsprintf(buf->pStr + buf->nUsedLength, lpszFormat, argList);
    // INCHI✔️❌:     va_end(argList);
    // INCHI✔️❌:     if (ret >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         buf->nUsedLength += ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_strbuf_printf

    inchi_strbuf_printf_arguments(
        heap,
        buf,
        format,
        &SourcePrintArguments {
            arguments: &arguments.arguments,
            position: arguments.position,
        },
    )
}

pub(crate) fn inchi_strbuf_printf_arguments(
    heap: &mut SourceHeap,
    buf: Option<&mut INCHI_IOS_STRING>,
    format: SourceConstPointer<i8>,
    arguments: &SourcePrintArguments<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1507 inchi_strbuf_printf
    // INCHI✔️❌: int inchi_strbuf_printf(INCHI_IOS_STRING* buf, const char* lpszFormat, ...)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0, max_len;
    // INCHI✔️❌:     va_list argList;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!buf)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     my_va_start(argList, lpszFormat);
    // INCHI✔️❌:     max_len = GetMaxPrintfLength(lpszFormat, argList);
    // INCHI✔️❌:     va_end(argList);
    // INCHI✔️❌:     if (max_len < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_strbuf_update(buf, max_len);
    // INCHI✔️❌:
    // INCHI✔️❌:     my_va_start(argList, lpszFormat);
    // INCHI✔️❌:     ret = vsprintf(buf->pStr + buf->nUsedLength, lpszFormat, argList);
    // INCHI✔️❌:     va_end(argList);
    // INCHI✔️❌:     if (ret >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         buf->nUsedLength += ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_strbuf_printf

    inchi_strbuf_printf_format_arguments(heap, buf, SourceFormat::Pointer(format), arguments)
}

#[rustfmt::skip]
pub(crate) fn inchi_strbuf_addline(
    heap: &mut SourceHeap,
    buf: &mut INCHI_IOS_STRING,
    inp_stream: Option<&mut INCHI_IOSTREAM>,
    crlf2lf: i32,
    preserve_lf: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1635 inchi_strbuf_addline
    // INCHI✔️❌: complete source frame follows verbatim; modeled format literal and typed varargs add overhead.
    /*
int inchi_strbuf_addline(INCHI_IOS_STRING* buf,
    INCHI_IOSTREAM* inp_stream,
    int crlf2lf,
    int preserve_lf)
{
    int c;

    while (1)
    {
        c = inchi_ios_str_getc(inp_stream);
        if (c == EOF)
        {
            return -1;
        }
        inchi_strbuf_printf(buf, "%c", c);
        if (c == '\n')
        {
            break;
        }
    }
    if (crlf2lf)
    {
        if (buf->nUsedLength > 2)
        {
            if (buf->pStr[buf->nUsedLength - 2] == '\r')
            {
                buf->pStr[buf->nUsedLength - 2] = '\n';
                buf->pStr[--buf->nUsedLength] = '\0';
            }
        }
    }
    if (!preserve_lf)
    {
        buf->pStr[--buf->nUsedLength] = '\0';
    }

    return buf->nUsedLength;
}
    */
    // END INCHI C FUNCTION: inchi_strbuf_addline
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: inchi_strbuf_addline
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; EOF is -1 and char is signed.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: inchi_strbuf_addline

    let inp_stream = inp_stream.ok_or(SourceHeapError::NullPointer)?;
    let format = crate::source::base::ichi_io::SourceFormat::Bytes(&[b'%' as i8, b'c' as i8, 0]);
    let result = (|| -> Result<i32, SourceHeapError> {
        loop {
            let character = inchi_ios_str_getc(heap, inp_stream)?;
            if character == SOURCE_EOF {
                return Ok(-1);
            }
            let arguments = SourcePrintArguments {
                arguments: &[SourceFormatArgument::Byte(character as i8)],
                position: 0,
            };
            let _ = crate::source::base::ichi_io::inchi_strbuf_printf_format_arguments(heap, Some(&mut *buf), format, &arguments)?;
            if character == i32::from(b'\n') {
                break;
            }
        }
        if crlf2lf != 0 && buf.nUsedLength > 2 {
            let cr_index = usize::try_from(buf.nUsedLength.wrapping_sub(2))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            if heap.slice(buf.pStr.as_const())?[cr_index] == b'\r' as i8 {
                heap.slice_mut(buf.pStr)?[cr_index] = b'\n' as i8;
                buf.nUsedLength = buf.nUsedLength.wrapping_sub(1);
                let end = usize::try_from(buf.nUsedLength)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                heap.slice_mut(buf.pStr)?[end] = 0;
            }
        }
        if preserve_lf == 0 {
            buf.nUsedLength = buf.nUsedLength.wrapping_sub(1);
            let end = usize::try_from(buf.nUsedLength)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            heap.slice_mut(buf.pStr)?[end] = 0;
        }
        Ok(buf.nUsedLength)
    })();
    result
}

// Concrete private printf-format storage projection; no owned format buffer.
// Bytes retain explicit first-NUL admission; Text models the old byte copy
// followed by one NUL, including any earlier embedded NUL. Pointer wrappers
// preserve their original typed-heap admission and state/error order.
// INCHI✔️❗: Format storage is O(1), eliminating the introduced O(F) copy,
// allocation, arena retention and free. Original estimate/render admission
// scans remain; dispatch/checked access and whole-call costs remain OPEN.
#[derive(Clone, Copy)]
pub(crate) enum SourceFormat<'a> {
    Pointer(SourceConstPointer<i8>),
    Bytes(&'a [i8]),
    Text(&'a str),
}

pub(crate) fn inchi_strbuf_printf_format_arguments(
    heap: &mut SourceHeap,
    buf: Option<&mut INCHI_IOS_STRING>,
    format: SourceFormat<'_>,
    arguments: &SourcePrintArguments<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1507 inchi_strbuf_printf
    // INCHI✔️❌: int inchi_strbuf_printf(INCHI_IOS_STRING* buf, const char* lpszFormat, ...)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0, max_len;
    // INCHI✔️❌:     va_list argList;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!buf)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     my_va_start(argList, lpszFormat);
    // INCHI✔️❌:     max_len = GetMaxPrintfLength(lpszFormat, argList);
    // INCHI✔️❌:     va_end(argList);
    // INCHI✔️❌:     if (max_len < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_strbuf_update(buf, max_len);
    // INCHI✔️❌:
    // INCHI✔️❌:     my_va_start(argList, lpszFormat);
    // INCHI✔️❌:     ret = vsprintf(buf->pStr + buf->nUsedLength, lpszFormat, argList);
    // INCHI✔️❌:     va_end(argList);
    // INCHI✔️❌:     if (ret >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         buf->nUsedLength += ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_strbuf_printf
    // INCHI✔️❗: Same canonical formatter and two independent argument cursors.
    inchi_strbuf_printf_arguments_core(heap, buf, format, arguments, None)
}

pub(crate) fn inchi_strbuf_printf_string(
    heap: &mut SourceHeap,
    buf: Option<&mut INCHI_IOS_STRING>,
    value: SourceFormat<'_>,
) -> Result<i32, SourceHeapError> {
    // INCHI source ichiprt2.c:457 / 2101:
    //         n = inchi_strbuf_printf( buf, "%s", szValue );
    //             inchi_strbuf_printf( strbuf, "%s", szValue );
    // INCHI✔️❗: Fixed %s adapter to the same canonical two-pass printer.
    // Borrow local scratch directly; no arena payload/owned argument buffer.
    // Existing update/format/error guards stay independent; dispatch cost OPEN.
    inchi_strbuf_printf_arguments_core(
        heap,
        buf,
        SourceFormat::Text("%s"),
        &SourcePrintArguments {
            arguments: &[],
            position: 0,
        },
        Some(value),
    )
}

pub(crate) fn inchi_strbuf_printf_string_format(
    heap: &mut SourceHeap,
    buf: Option<&mut INCHI_IOS_STRING>,
    value: SourceFormat<'_>,
) -> Result<i32, SourceHeapError> {
    // INCHI source ichiprt2.c:536:
    //     n = inchi_strbuf_printf( buf, "%-s", szTailingDelim );
    // INCHI✔️❗: One borrowed string argument, only the fixed source %-s call.
    // Both native va_start passes use the existing estimator/renderer core.
    // O(1) carrier storage; whole-call backend and dispatch costs remain OPEN.
    inchi_strbuf_printf_arguments_core(
        heap,
        buf,
        SourceFormat::Text("%-s"),
        &SourcePrintArguments {
            arguments: &[],
            position: 0,
        },
        Some(value),
    )
}

pub(crate) fn inchi_strbuf_printf_polymer_label(
    heap: &mut SourceHeap,
    buf: Option<&mut INCHI_IOS_STRING>,
    value: SourceFormat<'_>,
) -> Result<i32, SourceHeapError> {
    // INCHI source ichiprt1.c:4959:
    //     inchi_strbuf_printf( strbuf, "%-s-", u->smt[0] ? u->smt : "n" );
    // INCHI✔️❗: This fixed source call uses the existing two-pass formatter,
    // including its single borrowed string cursor and source growth order.
    // O(1) carrier, no arena payload; formatter/whole-call cost stays OPEN.
    inchi_strbuf_printf_arguments_core(
        heap,
        buf,
        SourceFormat::Text("%-s-"),
        &SourcePrintArguments {
            arguments: &[],
            position: 0,
        },
        Some(value),
    )
}

fn inchi_strbuf_printf_arguments_core<'a>(
    heap: &mut SourceHeap,
    buf: Option<&mut INCHI_IOS_STRING>,
    format: SourceFormat<'a>,
    arguments: &SourcePrintArguments<'a>,
    borrowed_string: Option<SourceFormat<'a>>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1507 inchi_strbuf_printf
    // INCHI✔️❌: int inchi_strbuf_printf(INCHI_IOS_STRING* buf, const char* lpszFormat, ...)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0, max_len;
    // INCHI✔️❌:     va_list argList;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!buf)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     my_va_start(argList, lpszFormat);
    // INCHI✔️❌:     max_len = GetMaxPrintfLength(lpszFormat, argList);
    // INCHI✔️❌:     va_end(argList);
    // INCHI✔️❌:     if (max_len < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_strbuf_update(buf, max_len);
    // INCHI✔️❌:
    // INCHI✔️❌:     my_va_start(argList, lpszFormat);
    // INCHI✔️❌:     ret = vsprintf(buf->pStr + buf->nUsedLength, lpszFormat, argList);
    // INCHI✔️❌:     va_end(argList);
    // INCHI✔️❌:     if (ret >= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         buf->nUsedLength += ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_strbuf_printf
    // INCHI✔️❗: Argument admission borrows a fixed slice and preserves the
    // initial cursor for each native va_start pass, without allocation or
    // cloning the list. Whole formatter/backend costs remain independently open.

    let Some(buf) = buf else {
        return Ok(-1);
    };
    let mut length_arguments_position = arguments.position;
    let mut length_arguments = FormatCursor {
        borrowed_string: borrowed_string
            .map(|value| BorrowedPrintStrings::One(0, BorrowedFormatString::Format(value))),
        arguments: arguments.arguments,
        position: &mut length_arguments_position,
    };
    let mut literal_format = false;
    let maximum_length = get_max_printf_length_format_cursor(
        heap,
        format,
        &mut length_arguments,
        Some(&mut literal_format),
    )?;
    if maximum_length < 0 {
        return Ok(0);
    }
    let updated = inchi_strbuf_update(heap, Some(&mut *buf), maximum_length)?;
    let direct_string = if literal_format {
        Some((format, Some(maximum_length as usize)))
    } else if let (Some(value), SourceFormat::Text("%s" | "%-s")) = (borrowed_string, format) {
        Some((value, None))
    } else {
        None
    };
    if let Some((value, literal_length)) = direct_string {
        // INCHI source ichi_io.c:1528:
        //     ret = vsprintf(buf->pStr + buf->nUsedLength, lpszFormat, argList);
        // INCHI✔️❗: Literal formats and the two declared borrowed-string
        // consumers emit only the first-NUL prefix: no conversions or %n.
        // Literal length comes from the existing estimator traversal; there
        // is no second format preflight or NUL scan. All other percent-bearing
        // formats, including %%, retain the original general renderer.
        // Keep the native estimate/update before the second string admission.
        // Copy directly into the actual destination, with no output-sized
        // temporary, allocation, string clone or second formatting algorithm.
        // Native ignores update failure; continue only if the complete output
        // including NUL fits the surviving backing. Otherwise preserve the
        // modeled AllocationFailed before any write (native would overflow).
        // O(string length) time/O(1) scratch. Checked heap access, carrier
        // dispatch and whole-call cost remain OPEN in performance_pending.md.
        let (string, implicit_nul) = format_string_bytes(heap, value)?;
        let length = if let Some(length) = literal_length {
            string
                .get(..length)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            length
        } else {
            memchr::memchr(0, bytemuck::cast_slice(string))
                .or_else(|| implicit_nul.then_some(string.len()))
                .ok_or(SourceHeapError::MissingNulTerminator)?
        };
        let rendered_length =
            i32::try_from(length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let fit_error = |error| {
            if updated < 0 {
                SourceHeapError::AllocationFailed
            } else {
                error
            }
        };
        let used = usize::try_from(buf.nUsedLength)
            .map_err(|_| fit_error(SourceHeapError::PointerOutOfBounds))?;
        let end = used
            .checked_add(length)
            .ok_or_else(|| fit_error(SourceHeapError::AllocationSizeOverflow))?;
        if updated < 0 {
            heap.slice(buf.pStr.as_const())
                .map_err(fit_error)?
                .get(end)
                .ok_or(SourceHeapError::AllocationFailed)?;
        }
        match value {
            SourceFormat::Pointer(source) => {
                let destination = buf.pStr.offset(
                    i64::try_from(used).map_err(|_| SourceHeapError::PointerOffsetOverflow)?,
                )?;
                // Disjoint ranges within one allocation are native-defined.
                // Reuse the existing safe same-allocation copy; this also
                // retains the previous modeled snapshot's overlap behavior.
                if destination.allocation_identity() == source.as_mut().allocation_identity() {
                    heap.copy_prefix_within_allocation(destination, source, length)?;
                } else {
                    heap.copy_nonoverlapping_prefix(destination, source, length, true)?;
                }
            }
            SourceFormat::Bytes(bytes) => {
                heap.slice_mut(buf.pStr)?
                    .get_mut(used..end)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .copy_from_slice(&bytes[..length]);
            }
            SourceFormat::Text(text) => {
                heap.slice_mut(buf.pStr)?
                    .get_mut(used..end)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .copy_from_slice(bytemuck::cast_slice(&text.as_bytes()[..length]));
            }
        }
        *heap
            .slice_mut(buf.pStr)?
            .get_mut(end)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
        buf.nUsedLength = buf.nUsedLength.wrapping_add(rendered_length);
        return Ok(rendered_length);
    }
    // INCHI source ichi_io.c:1525..1528:
    //     inchi_strbuf_update(buf, max_len);
    //
    //     my_va_start(argList, lpszFormat);
    //     ret = vsprintf(buf->pStr + buf->nUsedLength, lpszFormat, argList);
    // INCHI✔️❌: Failed conservative growth does not reject a short result
    // fitting old backing. Render once, preserving reached %n and cursor order;
    // never replay effects. On failed growth, admit the entire result and NUL
    // before buffer publication; no-fit/invalid backing retains AllocationFailed
    // instead of reproducing C overflow. Existing failed-allocation error
    // priority retained for Rust-only invalid formatting. Extra full-output
    // staging remains deferred in INCHI-PERF-PRINTF-FULL-OUTPUT.
    let fit_error = |error| {
        if updated < 0 {
            SourceHeapError::AllocationFailed
        } else {
            error
        }
    };
    let mut render_arguments_position = arguments.position;
    let mut render_arguments = FormatCursor {
        borrowed_string: borrowed_string
            .map(|value| BorrowedPrintStrings::One(0, BorrowedFormatString::Format(value))),
        arguments: arguments.arguments,
        position: &mut render_arguments_position,
    };
    let rendered =
        source_vformat_format_cursor(heap, format, &mut render_arguments).map_err(fit_error)?;
    let rendered_length = i32::try_from(rendered.len())
        .map_err(|_| fit_error(SourceHeapError::SourceIntegerOverflow))?;
    let used = usize::try_from(buf.nUsedLength)
        .map_err(|_| fit_error(SourceHeapError::PointerOutOfBounds))?;
    let end = used
        .checked_add(rendered.len())
        .ok_or_else(|| fit_error(SourceHeapError::AllocationSizeOverflow))?;
    if updated < 0 {
        heap.slice(buf.pStr.as_const())
            .map_err(fit_error)?
            .get(end)
            .ok_or(SourceHeapError::AllocationFailed)?;
    }
    let destination = heap.slice_mut(buf.pStr)?;
    destination
        .get_mut(used..end)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .copy_from_slice(bytemuck::cast_slice::<u8, i8>(&rendered));
    *destination
        .get_mut(end)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
    buf.nUsedLength = buf.nUsedLength.wrapping_add(rendered_length);
    Ok(rendered_length)
}

pub(crate) fn inchi_print_nodisplay_format_arguments(
    heap: &mut SourceHeap,
    file: Option<SourceMutPointer<FILE>>,
    stdout: SourceMutPointer<FILE>,
    format: SourceFormat<'_>,
    arguments: &SourcePrintArguments<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:867 inchi_print_nodisplay
    // INCHI✔️❌: int inchi_print_nodisplay(FILE* f, const char* lpszFormat, ...)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0;
    // INCHI✔️❌:     va_list argList;
    // INCHI✔️❌:     FILE* fi;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (f)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         fi = f;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         fi = stdout;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     my_va_start(argList, lpszFormat);
    // INCHI✔️❌:     ret = vfprintf(fi, lpszFormat, argList);
    // INCHI✔️❌:     va_end(argList);
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_print_nodisplay
    // INCHI✔️❗: Argument admission borrows a fixed slice and preserves the
    // initial cursor for each native va_start pass, without allocation or
    // cloning the list. Whole formatter/backend costs remain independently open.

    let destination = file.unwrap_or(stdout);
    let mut arguments_position = arguments.position;
    let mut arguments = FormatCursor {
        borrowed_string: None,
        arguments: arguments.arguments,
        position: &mut arguments_position,
    };
    inchi_vfprintf_format_cursor(heap, destination, format, &mut arguments)
}

fn source_stream_literal(
    heap: &mut SourceHeap,
    stream: &mut INCHI_IOS_STRING,
    format: SourceFormat<'_>,
    length: i32,
) -> Result<i32, SourceHeapError> {
    // INCHI✔️❗:             ret = vsprintf(ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList); /* djb-rwth: not using vsnprintf due to variable length argument */
    // INCHI✔️❗:             va_end(argList);
    // INCHI✔️❗:             if (ret >= 0)
    // INCHI✔️❗:             {
    // INCHI✔️❗:                 ios->s.nUsedLength += ret;
    // INCHI✔️❗:             }
    // INCHI✔️❗: Pure literals use the existing estimator's exact length;
    // its same traversal admitted first NUL (implicit only for Text). Re-admit
    // live storage after native growth, without another scan or output Vec.
    // Check the complete destination including NUL before any prefix write,
    // retaining this stream owner's checked-error state and update order.
    // O(length) time/O(1) scratch; arena/dispatch whole-call cost remains OPEN.
    let count = usize::try_from(length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    format_string_bytes(heap, format)?
        .0
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let used =
        usize::try_from(stream.nUsedLength).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let end = used
        .checked_add(count)
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    let required = end
        .checked_add(1)
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    heap.slice_mut(stream.pStr)?
        .get(..required)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    match format {
        SourceFormat::Pointer(source) => {
            let destination = stream
                .pStr
                .offset(i64::try_from(used).map_err(|_| SourceHeapError::PointerOffsetOverflow)?)?;
            if destination.allocation_identity() == source.as_mut().allocation_identity() {
                heap.copy_prefix_within_allocation(destination, source, count)?;
            } else {
                heap.copy_nonoverlapping_prefix(destination, source, count, true)?;
            }
        }
        SourceFormat::Bytes(bytes) => {
            heap.slice_mut(stream.pStr)?[used..end].copy_from_slice(&bytes[..count]);
        }
        SourceFormat::Text(text) => {
            heap.slice_mut(stream.pStr)?[used..end]
                .copy_from_slice(bytemuck::cast_slice(&text.as_bytes()[..count]));
        }
    }
    heap.slice_mut(stream.pStr)?[end] = 0;
    stream.nUsedLength = stream
        .nUsedLength
        .checked_add(length)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    Ok(length)
}

pub(crate) fn inchi_ios_print_format_arguments(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    format: SourceFormat<'_>,
    arguments: &SourcePrintArguments<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:477 inchi_ios_print
    // INCHI✔❌: int inchi_ios_print(INCHI_IOSTREAM* ios, const char* lpszFormat, ...)
    // INCHI✔❌: {
    // INCHI✔❌:     int ret = 0, ret2 = 0;
    // INCHI✔❌:     va_list argList;
    // INCHI✔❌:
    // INCHI✔❌:     if (!ios)
    // INCHI✔❌:     {
    // INCHI✔❌:         return -1;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (ios->type == INCHI_IOS_TYPE_STRING)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* output to string buffer */
    // INCHI✔❌:         int max_len;
    // INCHI✔❌:         my_va_start(argList, lpszFormat);
    // INCHI✔❌:         max_len = GetMaxPrintfLength(lpszFormat, argList);
    // INCHI✔❌:         va_end(argList);
    // INCHI✔❌:         if (max_len >= 0)
    // INCHI✔❌:         {
    // INCHI✔❌:             /* djb-rwth: fixing oss-fuzz issue #30152 */
    // INCHI✔❌:             int  nAddLength = inchi_max(INCHI_ADD_STR_LEN, max_len);
    // INCHI✔❌:             long long new_str_len = (long long)ios->s.nAllocatedLength + (long long)nAddLength;
    // INCHI✔❌:             if (ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len)
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* enlarge output string */
    // INCHI✔❌:                 char* new_str = (char*)inchi_calloc(new_str_len, sizeof(char)); /* djb-rwth: cast operators added */
    // INCHI✔❌:                 if (new_str)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     if (ios->s.pStr)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         if (ios->s.nUsedLength > 0)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             memcpy(new_str, ios->s.pStr, sizeof(new_str[0]) * ios->s.nUsedLength);
    // INCHI✔❌:                         }
    // INCHI✔❌:                         inchi_free(ios->s.pStr);
    // INCHI✔❌:                     }
    // INCHI✔❌:                     ios->s.pStr = new_str;
    // INCHI✔❌:                     ios->s.nAllocatedLength += nAddLength;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                     return -1; /* failed */
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:             /* output */
    // INCHI✔❌:             my_va_start(argList, lpszFormat);
    // INCHI✔❌:             ret = vsprintf(ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList); /* djb-rwth: not using vsnprintf due to variable length argument */
    // INCHI✔❌:             va_end(argList);
    // INCHI✔❌:             if (ret >= 0)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ios->s.nUsedLength += ret;
    // INCHI✔❌:             }
    // INCHI✔❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌: #if 0
    // INCHI✔❌:             if (FWPRINT)
    // INCHI✔❌:             {
    // INCHI✔❌:                 my_va_start(argList, lpszFormat);
    // INCHI✔❌:                 FWPRINT(lpszFormat, argList);
    // INCHI✔❌:                 va_end(argList);
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌: #endif
    // INCHI✔❌:             return ret;
    // INCHI✔❌:         }
    // INCHI✔❌:         return -1;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     else if (ios->type == INCHI_IOS_TYPE_FILE)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* output to file */
    // INCHI✔❌:         if (ios->f)
    // INCHI✔❌:         {
    // INCHI✔❌:             my_va_start(argList, lpszFormat);
    // INCHI✔❌:             ret = vfprintf(ios->f, lpszFormat, argList);
    // INCHI✔❌:             va_end(argList);
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:         {
    // INCHI✔❌:             my_va_start(argList, lpszFormat);
    // INCHI✔❌:             ret2 = vfprintf(stdout, lpszFormat, argList);
    // INCHI✔❌:             va_end(argList);
    // INCHI✔❌:         }
    // INCHI✔❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌:         if (FWPRINT)
    // INCHI✔❌:         {
    // INCHI✔❌:             my_va_start(argList, lpszFormat);
    // INCHI✔❌:             FWPRINT(lpszFormat, argList);
    // INCHI✔❌:             va_end(argList);
    // INCHI✔❌:         }
    // INCHI✔❌: #endif
    // INCHI✔❌:         return ret ? ret : ret2;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* no output */
    // INCHI✔❌:     return 0;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: inchi_ios_print
    // INCHI✔️❗: Argument admission borrows a fixed slice and preserves the
    // initial cursor for each native va_start pass, without allocation or
    // cloning the list. Whole formatter/backend costs remain independently open.

    source_print_arguments(heap, ios, stdout, format, arguments, None)
}

pub(crate) fn inchi_ios_print_format_strings<'a>(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    format: SourceFormat<'a>,
    arguments: &SourcePrintArguments<'a>,
    strings: &'a [Option<SourceFormat<'a>>],
) -> Result<i32, SourceHeapError> {
    // INCHI source ichiread.c11977:
    //     inchi_ios_print(is, "%-s%-s\n", edited_s, s2 ? s2 : "");
    // INCHI✔️❗: borrow concrete source strings through the existing cursor,
    // including its implicit static NUL. One canonical printer; no staging
    // allocation/copy/free. Argument/renderer/backend costs remain OPEN.
    source_print_arguments(
        heap,
        ios,
        stdout,
        format,
        arguments,
        Some(BorrowedPrintStrings::Indexed(strings)),
    )
}

fn source_print_arguments<'a>(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    format: SourceFormat<'a>,
    arguments: &SourcePrintArguments<'a>,
    borrowed_string: Option<BorrowedPrintStrings<'a>>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:477 inchi_ios_print
    // INCHI✔❌: int inchi_ios_print(INCHI_IOSTREAM* ios, const char* lpszFormat, ...)
    // INCHI✔❌: {
    // INCHI✔❌:     int ret = 0, ret2 = 0;
    // INCHI✔❌:     va_list argList;
    // INCHI✔❌:
    // INCHI✔❌:     if (!ios)
    // INCHI✔❌:     {
    // INCHI✔❌:         return -1;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (ios->type == INCHI_IOS_TYPE_STRING)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* output to string buffer */
    // INCHI✔❌:         int max_len;
    // INCHI✔❌:         my_va_start(argList, lpszFormat);
    // INCHI✔❌:         max_len = GetMaxPrintfLength(lpszFormat, argList);
    // INCHI✔❌:         va_end(argList);
    // INCHI✔❌:         if (max_len >= 0)
    // INCHI✔❌:         {
    // INCHI✔❌:             /* djb-rwth: fixing oss-fuzz issue #30152 */
    // INCHI✔❌:             int  nAddLength = inchi_max(INCHI_ADD_STR_LEN, max_len);
    // INCHI✔❌:             long long new_str_len = (long long)ios->s.nAllocatedLength + (long long)nAddLength;
    // INCHI✔❌:             if (ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len)
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* enlarge output string */
    // INCHI✔❌:                 char* new_str = (char*)inchi_calloc(new_str_len, sizeof(char)); /* djb-rwth: cast operators added */
    // INCHI✔❌:                 if (new_str)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     if (ios->s.pStr)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         if (ios->s.nUsedLength > 0)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             memcpy(new_str, ios->s.pStr, sizeof(new_str[0]) * ios->s.nUsedLength);
    // INCHI✔❌:                         }
    // INCHI✔❌:                         inchi_free(ios->s.pStr);
    // INCHI✔❌:                     }
    // INCHI✔❌:                     ios->s.pStr = new_str;
    // INCHI✔❌:                     ios->s.nAllocatedLength += nAddLength;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                     return -1; /* failed */
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:             /* output */
    // INCHI✔❌:             my_va_start(argList, lpszFormat);
    // INCHI✔❌:             ret = vsprintf(ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList); /* djb-rwth: not using vsnprintf due to variable length argument */
    // INCHI✔❌:             va_end(argList);
    // INCHI✔❌:             if (ret >= 0)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ios->s.nUsedLength += ret;
    // INCHI✔❌:             }
    // INCHI✔❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌: #if 0
    // INCHI✔❌:             if (FWPRINT)
    // INCHI✔❌:             {
    // INCHI✔❌:                 my_va_start(argList, lpszFormat);
    // INCHI✔❌:                 FWPRINT(lpszFormat, argList);
    // INCHI✔❌:                 va_end(argList);
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌: #endif
    // INCHI✔❌:             return ret;
    // INCHI✔❌:         }
    // INCHI✔❌:         return -1;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     else if (ios->type == INCHI_IOS_TYPE_FILE)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* output to file */
    // INCHI✔❌:         if (ios->f)
    // INCHI✔❌:         {
    // INCHI✔❌:             my_va_start(argList, lpszFormat);
    // INCHI✔❌:             ret = vfprintf(ios->f, lpszFormat, argList);
    // INCHI✔❌:             va_end(argList);
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:         {
    // INCHI✔❌:             my_va_start(argList, lpszFormat);
    // INCHI✔❌:             ret2 = vfprintf(stdout, lpszFormat, argList);
    // INCHI✔❌:             va_end(argList);
    // INCHI✔❌:         }
    // INCHI✔❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌:         if (FWPRINT)
    // INCHI✔❌:         {
    // INCHI✔❌:             my_va_start(argList, lpszFormat);
    // INCHI✔❌:             FWPRINT(lpszFormat, argList);
    // INCHI✔❌:             va_end(argList);
    // INCHI✔❌:         }
    // INCHI✔❌: #endif
    // INCHI✔❌:         return ret ? ret : ret2;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* no output */
    // INCHI✔❌:     return 0;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: inchi_ios_print
    // INCHI✔️❗: Argument admission borrows a fixed slice and preserves the
    // initial cursor for each native va_start pass, without allocation or
    // cloning the list. Whole formatter/backend costs remain independently open.

    let Some(ios) = ios else {
        return Ok(-1);
    };
    if ios.type_ == INCHI_IOS_TYPE_STRING as i32 {
        let mut estimate_arguments_position = arguments.position;
        let mut estimate_arguments = FormatCursor {
            borrowed_string,
            arguments: arguments.arguments,
            position: &mut estimate_arguments_position,
        };
        let mut literal_format = false;
        let maximum_length = get_max_printf_length_format_cursor(
            heap,
            format,
            &mut estimate_arguments,
            Some(&mut literal_format),
        )?;
        if maximum_length < 0 {
            return Ok(-1);
        }
        // The source estimate counts printable bytes, not the trailing NUL
        // written by vsprintf. Reserve that slot explicitly in safe storage;
        // retain the source growth quantum and formatted output unchanged.
        let required_length = maximum_length
            .checked_add(1)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        let additional_length = (INCHI_ADD_STR_LEN as i32).max(required_length);
        let new_string_length = i64::from(ios.s.nAllocatedLength)
            .checked_add(i64::from(additional_length))
            .ok_or(SourceHeapError::AllocationSizeOverflow)?;
        if ios.s.nAllocatedLength - ios.s.nUsedLength <= maximum_length {
            let new_string = match inchi_calloc::<i8>(heap, new_string_length as u64, 1) {
                Ok(pointer) => pointer,
                Err(_) => return Ok(-1),
            };
            if !ios.s.pStr.is_null() {
                if ios.s.nUsedLength > 0 {
                    let used = usize::try_from(ios.s.nUsedLength)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    // Fresh destination and old live storage are disjoint;
                    // copy directly without an intermediate prefix allocation.
                    heap.slice(ios.s.pStr.as_const())?
                        .get(..used)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    heap.with_slice_mut_and_heap(new_string, |destination, heap| {
                        let source = heap.slice(ios.s.pStr.as_const())?;
                        destination[..used].copy_from_slice(&source[..used]);
                        Ok(())
                    })?;
                }
                inchi_free(heap, ios.s.pStr)?;
            }
            ios.s.pStr = new_string;
            ios.s.nAllocatedLength = ios
                .s
                .nAllocatedLength
                .checked_add(additional_length)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        }
        if literal_format {
            return source_stream_literal(heap, &mut ios.s, format, maximum_length);
        }
        let mut output_arguments_position = arguments.position;
        let mut output_arguments = FormatCursor {
            borrowed_string,
            arguments: arguments.arguments,
            position: &mut output_arguments_position,
        };
        let output = source_vformat_format_cursor(heap, format, &mut output_arguments)?;
        let used =
            usize::try_from(ios.s.nUsedLength).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let destination = heap.slice_mut(ios.s.pStr)?;
        let end = used
            .checked_add(output.len())
            .ok_or(SourceHeapError::AllocationSizeOverflow)?;
        if end >= destination.len() {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        for (destination, source) in destination[used..end].iter_mut().zip(&output) {
            *destination = *source as i8;
        }
        destination[end] = 0;
        let result =
            i32::try_from(output.len()).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        ios.s.nUsedLength = ios
            .s
            .nUsedLength
            .checked_add(result)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        return Ok(result);
    }
    if ios.type_ == INCHI_IOS_TYPE_FILE as i32 {
        let destination = if ios.f.is_null() { stdout } else { ios.f };
        let mut arguments_position = arguments.position;
        let mut arguments = FormatCursor {
            borrowed_string,
            arguments: arguments.arguments,
            position: &mut arguments_position,
        };
        return inchi_vfprintf_format_cursor(heap, destination, format, &mut arguments);
    }
    Ok(0)
}

pub(crate) fn inchi_ios_print_nodisplay_format_arguments(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    format: SourceFormat<'_>,
    arguments: &SourcePrintArguments<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:604 inchi_ios_print_nodisplay
    // INCHI✔️❌: int inchi_ios_print_nodisplay( INCHI_IOSTREAM * ios,
    // INCHI✔️❌:                                const char* lpszFormat,
    // INCHI✔️❌:                                ... )
    // INCHI✔️❌: {
    // INCHI✔️❌:     va_list argList;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!ios)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ios->type == INCHI_IOS_TYPE_STRING)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* output to string buffer */
    // INCHI✔️❌:         int ret = 0, max_len;
    // INCHI✔️❌:         my_va_start( argList, lpszFormat );
    // INCHI✔️❌:         max_len = GetMaxPrintfLength( lpszFormat, argList );
    // INCHI✔️❌:         va_end( argList );
    // INCHI✔️❌:         if (max_len >= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* djb-rwth: fixing oss-fuzz issue #30163 */
    // INCHI✔️❌:             int  nAddLength = inchi_max(INCHI_ADD_STR_LEN, max_len);
    // INCHI✔️❌:             long long new_str_len = (long long)ios->s.nAllocatedLength + (long long)nAddLength;
    // INCHI✔️❌:             if (ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* enlarge output string */
    // INCHI✔️❌:                 char* new_str = (char*)inchi_calloc(new_str_len, sizeof(new_str[0])); /* djb-rwth: cast operators added */
    // INCHI✔️❌:                 if (new_str)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (ios->s.pStr)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (ios->s.nUsedLength > 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             memcpy( new_str, ios->s.pStr, sizeof( new_str[0] )*ios->s.nUsedLength );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         inchi_free( ios->s.pStr );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     ios->s.pStr = new_str;
    // INCHI✔️❌:                     ios->s.nAllocatedLength += nAddLength;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return -1; /* failed */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* output */
    // INCHI✔️❌:             /* djb-rwth: fixing oss-fuzz issue #67676 */
    // INCHI✔️❌:             my_va_start(argList, lpszFormat);
    // INCHI✔️❌:             ret = vsprintf(ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList); /* djb-rwth: not using vsnprintf due to variable length argument; fixing GHI #71 */
    // INCHI✔️❌:             va_end(argList);
    // INCHI✔️❌:             if (ret >= 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ios->s.nUsedLength += ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     else if (ios->type == INCHI_IOS_TYPE_FILE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         my_va_start( argList, lpszFormat );
    // INCHI✔️❌:         inchi_print_nodisplay( ios->f, lpszFormat, argList );
    // INCHI✔️❌:         va_end( argList );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* no output */
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_ios_print_nodisplay
    // INCHI✔️❗: Argument admission borrows a fixed slice and preserves the
    // initial cursor for each native va_start pass, without allocation or
    // cloning the list. Whole formatter/backend costs remain independently open.

    source_nodisplay_arguments(heap, ios, stdout, format, arguments, None)
}

pub(crate) fn inchi_ios_print_nodisplay_content_argument<'a>(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    format: SourceFormat<'a>,
    arguments: &SourcePrintArguments<'a>,
    slot: usize,
    content: &'a [u8],
) -> Result<i32, SourceHeapError> {
    // Native mol_fmt4.c882/921/978/831/837/1103/1221/1250/1367 passes borrowed
    // string arguments directly. One concrete indexed content slot, no arena
    // registration or owned argument vector; canonical printer stays shared.
    source_nodisplay_arguments(
        heap,
        ios,
        stdout,
        format,
        arguments,
        Some(BorrowedPrintStrings::One(
            slot,
            BorrowedFormatString::Content(content),
        )),
    )
}

pub(crate) fn inchi_ios_print_nodisplay_format_strings<'a>(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    format: SourceFormat<'a>,
    arguments: &SourcePrintArguments<'a>,
    strings: &'a [Option<SourceFormat<'a>>],
) -> Result<i32, SourceHeapError> {
    // INCHI✔️❗: int inchi_ios_print_nodisplay( INCHI_IOSTREAM * ios,
    // INCHI✔️❗:                                const char* lpszFormat,
    // INCHI✔️❗:                                ... )
    // INCHI✔️❗: {
    // INCHI✔️❗:     va_list argList;
    // INCHI✔️❗:
    // INCHI✔️❗:     if (!ios)
    // INCHI✔️❗:     {
    // INCHI✔️❗:         return -1;
    // INCHI✔️❗:     }
    // INCHI✔️❗:
    // INCHI✔️❗:     if (ios->type == INCHI_IOS_TYPE_STRING)
    // INCHI✔️❗:     {
    // INCHI✔️❗:         /* output to string buffer */
    // INCHI✔️❗:         int ret = 0, max_len;
    // INCHI✔️❗:         my_va_start( argList, lpszFormat );
    // INCHI✔️❗:         max_len = GetMaxPrintfLength( lpszFormat, argList );
    // INCHI✔️❗:         va_end( argList );
    // INCHI✔️❗:         if (max_len >= 0)
    // INCHI✔️❗:         {
    // INCHI✔️❗:             /* djb-rwth: fixing oss-fuzz issue #30163 */
    // INCHI✔️❗:             int  nAddLength = inchi_max(INCHI_ADD_STR_LEN, max_len);
    // INCHI✔️❗:             long long new_str_len = (long long)ios->s.nAllocatedLength + (long long)nAddLength;
    // INCHI✔️❗:             if (ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len)
    // INCHI✔️❗:             {
    // INCHI✔️❗:                 /* enlarge output string */
    // INCHI✔️❗:                 char* new_str = (char*)inchi_calloc(new_str_len, sizeof(new_str[0])); /* djb-rwth: cast operators added */
    // INCHI✔️❗:                 if (new_str)
    // INCHI✔️❗:                 {
    // INCHI✔️❗:                     if (ios->s.pStr)
    // INCHI✔️❗:                     {
    // INCHI✔️❗:                         if (ios->s.nUsedLength > 0)
    // INCHI✔️❗:                         {
    // INCHI✔️❗:                             memcpy( new_str, ios->s.pStr, sizeof( new_str[0] )*ios->s.nUsedLength );
    // INCHI✔️❗:                         }
    // INCHI✔️❗:                         inchi_free( ios->s.pStr );
    // INCHI✔️❗:                     }
    // INCHI✔️❗:                     ios->s.pStr = new_str;
    // INCHI✔️❗:                     ios->s.nAllocatedLength += nAddLength;
    // INCHI✔️❗:                 }
    // INCHI✔️❗:                 else
    // INCHI✔️❗:                 {
    // INCHI✔️❗:                     return -1; /* failed */
    // INCHI✔️❗:                 }
    // INCHI✔️❗:             }
    // INCHI✔️❗:             /* output */
    // INCHI✔️❗:             /* djb-rwth: fixing oss-fuzz issue #67676 */
    // INCHI✔️❗:             my_va_start(argList, lpszFormat);
    // INCHI✔️❗:             ret = vsprintf(ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList); /* djb-rwth: not using vsnprintf due to variable length argument; fixing GHI #71 */
    // INCHI✔️❗:             va_end(argList);
    // INCHI✔️❗:             if (ret >= 0)
    // INCHI✔️❗:             {
    // INCHI✔️❗:                 ios->s.nUsedLength += ret;
    // INCHI✔️❗:             }
    // INCHI✔️❗:             return ret;
    // INCHI✔️❗:         }
    // INCHI✔️❗:         return -1;
    // INCHI✔️❗:     }
    // INCHI✔️❗:
    // INCHI✔️❗:     else if (ios->type == INCHI_IOS_TYPE_FILE)
    // INCHI✔️❗:     {
    // INCHI✔️❗:         my_va_start( argList, lpszFormat );
    // INCHI✔️❗:         inchi_print_nodisplay( ios->f, lpszFormat, argList );
    // INCHI✔️❗:         va_end( argList );
    // INCHI✔️❗:     }
    // INCHI✔️❗:
    // INCHI✔️❗:     /* no output */
    // INCHI✔️❗:     return 0;
    // INCHI✔️❗: }
    // INCHI✔️❗: Concrete ichiparm variadic strings use the existing cursor
    // and complete canonical printer. No formatter duplication, argument
    // allocation, string copy or eager NUL/error admission. Indexed dispatch
    // and ABI-relative whole-call cost remain OPEN, never a performance pass.
    source_nodisplay_arguments(
        heap,
        ios,
        stdout,
        format,
        arguments,
        Some(BorrowedPrintStrings::Indexed(strings)),
    )
}

fn source_nodisplay_arguments<'a>(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    format: SourceFormat<'a>,
    arguments: &SourcePrintArguments<'a>,
    borrowed_string: Option<BorrowedPrintStrings<'a>>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:604 inchi_ios_print_nodisplay
    // INCHI✔️❌: int inchi_ios_print_nodisplay( INCHI_IOSTREAM * ios,
    // INCHI✔️❌:                                const char* lpszFormat,
    // INCHI✔️❌:                                ... )
    // INCHI✔️❌: {
    // INCHI✔️❌:     va_list argList;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!ios)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ios->type == INCHI_IOS_TYPE_STRING)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* output to string buffer */
    // INCHI✔️❌:         int ret = 0, max_len;
    // INCHI✔️❌:         my_va_start( argList, lpszFormat );
    // INCHI✔️❌:         max_len = GetMaxPrintfLength( lpszFormat, argList );
    // INCHI✔️❌:         va_end( argList );
    // INCHI✔️❌:         if (max_len >= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* djb-rwth: fixing oss-fuzz issue #30163 */
    // INCHI✔️❌:             int  nAddLength = inchi_max(INCHI_ADD_STR_LEN, max_len);
    // INCHI✔️❌:             long long new_str_len = (long long)ios->s.nAllocatedLength + (long long)nAddLength;
    // INCHI✔️❌:             if (ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* enlarge output string */
    // INCHI✔️❌:                 char* new_str = (char*)inchi_calloc(new_str_len, sizeof(new_str[0])); /* djb-rwth: cast operators added */
    // INCHI✔️❌:                 if (new_str)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (ios->s.pStr)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (ios->s.nUsedLength > 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             memcpy( new_str, ios->s.pStr, sizeof( new_str[0] )*ios->s.nUsedLength );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         inchi_free( ios->s.pStr );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     ios->s.pStr = new_str;
    // INCHI✔️❌:                     ios->s.nAllocatedLength += nAddLength;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return -1; /* failed */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* output */
    // INCHI✔️❌:             /* djb-rwth: fixing oss-fuzz issue #67676 */
    // INCHI✔️❌:             my_va_start(argList, lpszFormat);
    // INCHI✔️❌:             ret = vsprintf(ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList); /* djb-rwth: not using vsnprintf due to variable length argument; fixing GHI #71 */
    // INCHI✔️❌:             va_end(argList);
    // INCHI✔️❌:             if (ret >= 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ios->s.nUsedLength += ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     else if (ios->type == INCHI_IOS_TYPE_FILE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         my_va_start( argList, lpszFormat );
    // INCHI✔️❌:         inchi_print_nodisplay( ios->f, lpszFormat, argList );
    // INCHI✔️❌:         va_end( argList );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* no output */
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_ios_print_nodisplay
    // INCHI✔️❗: Argument admission borrows a fixed slice and preserves the
    // initial cursor for each native va_start pass, without allocation or
    // cloning the list. Whole formatter/backend costs remain independently open.

    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:867 inchi_print_nodisplay
    // INCHI✔️❌: int inchi_print_nodisplay(FILE* f, const char* lpszFormat, ...)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0;
    // INCHI✔️❌:     va_list argList;
    // INCHI✔️❌:     FILE* fi;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (f)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         fi = f;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         fi = stdout;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     my_va_start(argList, lpszFormat);
    // INCHI✔️❌:     ret = vfprintf(fi, lpszFormat, argList);
    // INCHI✔️❌:     va_end(argList);
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_print_nodisplay
    // INCHI✔️❗: Shared full STRING/FILE printer, one borrowed string slot.
    // Independent va_start positions, type/error advancement, allocation and
    // publication retained. FILE keeps stdout fallback and discarded return.
    // No argument copy/allocation; renderer guards remain. Indexed dispatch
    // and whole-call cost OPEN INCHI-PERF-PRINTF-CONTENT-ARGUMENT.
    let Some(ios) = ios else {
        return Ok(-1);
    };
    if ios.type_ == INCHI_IOS_TYPE_STRING as i32 {
        let mut estimate_arguments_position = arguments.position;
        let mut estimate_arguments = FormatCursor {
            borrowed_string,
            arguments: arguments.arguments,
            position: &mut estimate_arguments_position,
        };
        let mut literal_format = false;
        let maximum_length = get_max_printf_length_format_cursor(
            heap,
            format,
            &mut estimate_arguments,
            Some(&mut literal_format),
        )?;
        if maximum_length < 0 {
            return Ok(-1);
        }
        // vsprintf also writes the terminating NUL, outside max_len's count.
        let required_length = maximum_length
            .checked_add(1)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        let additional_length = (INCHI_ADD_STR_LEN as i32).max(required_length);
        let new_string_length = i64::from(ios.s.nAllocatedLength)
            .checked_add(i64::from(additional_length))
            .ok_or(SourceHeapError::AllocationSizeOverflow)?;
        if ios.s.nAllocatedLength - ios.s.nUsedLength <= maximum_length {
            let new_string = match inchi_calloc::<i8>(heap, new_string_length as u64, 1) {
                Ok(pointer) => pointer,
                Err(_) => return Ok(-1),
            };
            if !ios.s.pStr.is_null() {
                if ios.s.nUsedLength > 0 {
                    let used = usize::try_from(ios.s.nUsedLength)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    // Fresh destination and old live storage are disjoint;
                    // copy directly without an intermediate prefix allocation.
                    let _ = &heap.slice(ios.s.pStr.as_const())?[..used];
                    heap.with_slice_mut_and_heap(new_string, |destination, heap| {
                        let source = heap.slice(ios.s.pStr.as_const())?;
                        destination[..used].copy_from_slice(&source[..used]);
                        Ok(())
                    })?;
                }
                inchi_free(heap, ios.s.pStr)?;
            }
            ios.s.pStr = new_string;
            ios.s.nAllocatedLength = ios
                .s
                .nAllocatedLength
                .checked_add(additional_length)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        }
        if literal_format {
            return source_stream_literal(heap, &mut ios.s, format, maximum_length);
        }
        let mut output_arguments_position = arguments.position;
        let mut output_arguments = FormatCursor {
            borrowed_string,
            arguments: arguments.arguments,
            position: &mut output_arguments_position,
        };
        let output = source_vformat_format_cursor(heap, format, &mut output_arguments)?;
        let used =
            usize::try_from(ios.s.nUsedLength).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let destination = heap.slice_mut(ios.s.pStr)?;
        let end = used
            .checked_add(output.len())
            .ok_or(SourceHeapError::AllocationSizeOverflow)?;
        if end >= destination.len() {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        for (destination, source) in destination[used..end].iter_mut().zip(&output) {
            *destination = *source as i8;
        }
        destination[end] = 0;
        let result =
            i32::try_from(output.len()).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        ios.s.nUsedLength = ios
            .s
            .nUsedLength
            .checked_add(result)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        return Ok(result);
    }
    if ios.type_ == INCHI_IOS_TYPE_FILE as i32 {
        let destination = if ios.f.is_null() { stdout } else { ios.f };
        let mut position = arguments.position;
        let mut cursor = FormatCursor {
            borrowed_string,
            arguments: arguments.arguments,
            position: &mut position,
        };
        inchi_vfprintf_format_cursor(heap, destination, format, &mut cursor)?;
    }
    Ok(0)
}

pub(crate) fn inchi_ios_eprint_format_arguments(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    format: SourceFormat<'_>,
    arguments: &SourcePrintArguments<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:708 inchi_ios_eprint
    // INCHI✔️❌: int inchi_ios_eprint(INCHI_IOSTREAM* ios, const char* lpszFormat, ...)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0, ret2 = 0;
    // INCHI✔️❌:     va_list argList;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!ios)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ios->type == INCHI_IOS_TYPE_STRING)
    // INCHI✔️❌:         /* was #if ( defined(TARGET_API_LIB) || defined(INCHI_STANDALONE_EXE) ) */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* output to string buffer */
    // INCHI✔️❌:         int max_len, nAddLength = 0;
    // INCHI✔️❌:         char* new_str = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:         my_va_start(argList, lpszFormat);
    // INCHI✔️❌:         max_len = GetMaxPrintfLength(lpszFormat, argList);
    // INCHI✔️❌:         va_end(argList);
    // INCHI✔️❌:
    // INCHI✔️❌:         if (max_len >= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* enlarge output string */
    // INCHI✔️❌:                 nAddLength = inchi_max(INCHI_ADD_STR_LEN, max_len);
    // INCHI✔️❌:                 new_str = (char*)inchi_calloc((long long)ios->s.nAllocatedLength + (long long)nAddLength, sizeof(new_str[0])); /* djb-rwth: cast operators added */
    // INCHI✔️❌:                 if (new_str)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (ios->s.pStr)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (ios->s.nUsedLength > 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             memcpy(new_str, ios->s.pStr, sizeof(new_str[0]) * ios->s.nUsedLength);
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         inchi_free(ios->s.pStr);
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     ios->s.pStr = new_str;
    // INCHI✔️❌:                     ios->s.nAllocatedLength += nAddLength;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return -1; /* failed */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* output */
    // INCHI✔️❌:             my_va_start(argList, lpszFormat);
    // INCHI✔️❌:             ret = vsprintf(ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList);
    // INCHI✔️❌:             va_end(argList);
    // INCHI✔️❌:             if (ret >= 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ios->s.nUsedLength += ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     else if (ios->type == INCHI_IOS_TYPE_FILE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (ios->f)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* output to plain file */
    // INCHI✔️❌:             my_va_start(argList, lpszFormat);
    // INCHI✔️❌:             ret = inchi_vfprintf(ios->f, lpszFormat, argList);
    // INCHI✔️❌:             va_end(argList);
    // INCHI✔️❌:             /*  No output to stderr from within dll or GUI program */
    // INCHI✔️❌:             return ret ? ret : ret2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* no output */
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_ios_eprint
    // INCHI✔️❗: Argument admission borrows a fixed slice and preserves the
    // initial cursor for each native va_start pass, without allocation or
    // cloning the list. Whole formatter/backend costs remain independently open.

    source_eprint_arguments(heap, ios, format, arguments, None)
}

pub(crate) fn inchi_ios_eprint_content(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    bytes: &[u8],
) -> Result<i32, SourceHeapError> {
    // INCHI source runichi2.c370/732/735/768/790/813: inchi_ios_eprint.
    // Borrow existing message bytes directly, preserving first-NUL semantics.
    // O(1) carrier, no staging allocation/copy/free; independent renderer and
    // full printer costs remain OPEN INCHI-PERF-EPRINT-BORROWED-CONTENT.
    source_eprint_arguments(
        heap,
        ios,
        SourceFormat::Bytes(&[b'%' as i8, b's' as i8, 0]),
        &SourcePrintArguments {
            arguments: &[],
            position: 0,
        },
        Some(BorrowedPrintStrings::One(
            0,
            BorrowedFormatString::Content(bytes),
        )),
    )
}

pub(crate) fn inchi_ios_eprint_format_strings<'a>(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    format: SourceFormat<'a>,
    arguments: &SourcePrintArguments<'a>,
    strings: &'a [Option<SourceFormat<'a>>],
) -> Result<i32, SourceHeapError> {
    // INCHI source runichi2.c:732:
    //             inchi_ios_eprint( log_file, "%s inp structure #%ld: End of file.%s%s%s%s    \n", sd->pStrErrStruct, *num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
    // INCHI✔️❗: Borrow concrete indexed strings through the same complete
    // eprint kernel/cursor. No owned argument/message buffer or arena staging.
    // Explicit byte arrays retain required-NUL errors only when reached;
    // Text keeps native static literals. Whole-call cost remains OPEN.
    source_eprint_arguments(
        heap,
        ios,
        format,
        arguments,
        Some(BorrowedPrintStrings::Indexed(strings)),
    )
}

fn source_eprint_arguments<'a>(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    format: SourceFormat<'a>,
    arguments: &SourcePrintArguments<'a>,
    borrowed_string: Option<BorrowedPrintStrings<'a>>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:708 inchi_ios_eprint
    // INCHI✔️❌: int inchi_ios_eprint(INCHI_IOSTREAM* ios, const char* lpszFormat, ...)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0, ret2 = 0;
    // INCHI✔️❌:     va_list argList;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!ios)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ios->type == INCHI_IOS_TYPE_STRING)
    // INCHI✔️❌:         /* was #if ( defined(TARGET_API_LIB) || defined(INCHI_STANDALONE_EXE) ) */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* output to string buffer */
    // INCHI✔️❌:         int max_len, nAddLength = 0;
    // INCHI✔️❌:         char* new_str = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:         my_va_start(argList, lpszFormat);
    // INCHI✔️❌:         max_len = GetMaxPrintfLength(lpszFormat, argList);
    // INCHI✔️❌:         va_end(argList);
    // INCHI✔️❌:
    // INCHI✔️❌:         if (max_len >= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* enlarge output string */
    // INCHI✔️❌:                 nAddLength = inchi_max(INCHI_ADD_STR_LEN, max_len);
    // INCHI✔️❌:                 new_str = (char*)inchi_calloc((long long)ios->s.nAllocatedLength + (long long)nAddLength, sizeof(new_str[0])); /* djb-rwth: cast operators added */
    // INCHI✔️❌:                 if (new_str)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (ios->s.pStr)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (ios->s.nUsedLength > 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             memcpy(new_str, ios->s.pStr, sizeof(new_str[0]) * ios->s.nUsedLength);
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         inchi_free(ios->s.pStr);
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     ios->s.pStr = new_str;
    // INCHI✔️❌:                     ios->s.nAllocatedLength += nAddLength;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return -1; /* failed */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* output */
    // INCHI✔️❌:             my_va_start(argList, lpszFormat);
    // INCHI✔️❌:             ret = vsprintf(ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList);
    // INCHI✔️❌:             va_end(argList);
    // INCHI✔️❌:             if (ret >= 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ios->s.nUsedLength += ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     else if (ios->type == INCHI_IOS_TYPE_FILE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (ios->f)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* output to plain file */
    // INCHI✔️❌:             my_va_start(argList, lpszFormat);
    // INCHI✔️❌:             ret = inchi_vfprintf(ios->f, lpszFormat, argList);
    // INCHI✔️❌:             va_end(argList);
    // INCHI✔️❌:             /*  No output to stderr from within dll or GUI program */
    // INCHI✔️❌:             return ret ? ret : ret2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* no output */
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_ios_eprint
    // INCHI✔️❗: Argument admission borrows a fixed slice and preserves the
    // initial cursor for each native va_start pass, without allocation or
    // cloning the list. Whole formatter/backend costs remain independently open.

    // INCHI✔️❗: One full canonical eprint kernel; original cursor starts,
    // STRING publication/failure and FILE null/return preserved. O(1) carrier,
    // no argument allocation/clone; whole-call cost deferred.
    let Some(ios) = ios else {
        return Ok(-1);
    };
    if ios.type_ == INCHI_IOS_TYPE_STRING as i32 {
        let mut estimate_arguments_position = arguments.position;
        let mut estimate_arguments = FormatCursor {
            borrowed_string,
            arguments: arguments.arguments,
            position: &mut estimate_arguments_position,
        };
        let mut literal_format = false;
        let maximum_length = get_max_printf_length_format_cursor(
            heap,
            format,
            &mut estimate_arguments,
            Some(&mut literal_format),
        )?;
        if maximum_length < 0 {
            return Ok(-1);
        }
        if ios.s.nAllocatedLength - ios.s.nUsedLength <= maximum_length {
            // vsprintf also writes the terminating NUL, outside max_len's count.
            let required_length = maximum_length
                .checked_add(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            let additional_length = (INCHI_ADD_STR_LEN as i32).max(required_length);
            let allocation_length = i64::from(ios.s.nAllocatedLength)
                .checked_add(i64::from(additional_length))
                .ok_or(SourceHeapError::AllocationSizeOverflow)?;
            let new_string = match inchi_calloc::<i8>(heap, allocation_length as u64, 1) {
                Ok(pointer) => pointer,
                Err(_) => return Ok(-1),
            };
            if !ios.s.pStr.is_null() {
                if ios.s.nUsedLength > 0 {
                    let used = usize::try_from(ios.s.nUsedLength)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let old = heap.slice(ios.s.pStr.as_const())?;
                    if old.len() < used {
                        return Err(SourceHeapError::PointerOutOfBounds);
                    }
                    // Fresh destination and old live storage are disjoint;
                    // copy directly without an intermediate prefix allocation.
                    heap.with_slice_mut_and_heap(new_string, |destination, heap| {
                        let source = heap.slice(ios.s.pStr.as_const())?;
                        destination[..used].copy_from_slice(&source[..used]);
                        Ok(())
                    })?;
                }
                inchi_free(heap, ios.s.pStr)?;
            }
            ios.s.pStr = new_string;
            ios.s.nAllocatedLength = ios
                .s
                .nAllocatedLength
                .checked_add(additional_length)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        }

        if literal_format {
            return source_stream_literal(heap, &mut ios.s, format, maximum_length);
        }
        let mut output_arguments_position = arguments.position;
        let mut output_arguments = FormatCursor {
            borrowed_string,
            arguments: arguments.arguments,
            position: &mut output_arguments_position,
        };
        let output = source_vformat_format_cursor(heap, format, &mut output_arguments)?;
        let used =
            usize::try_from(ios.s.nUsedLength).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let required = used
            .checked_add(output.len())
            .and_then(|value| value.checked_add(1))
            .ok_or(SourceHeapError::AllocationSizeOverflow)?;
        let destination = heap.slice_mut(ios.s.pStr)?;
        if destination.len() < required {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        // u8/i8 have identical byte representation; borrow the rendered bytes
        // rather than allocating and filling a second output vector.
        destination[used..used + output.len()]
            .copy_from_slice(bytemuck::cast_slice::<u8, i8>(&output));
        destination[used + output.len()] = 0;
        let result =
            i32::try_from(output.len()).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        ios.s.nUsedLength = ios
            .s
            .nUsedLength
            .checked_add(result)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        return Ok(result);
    }

    if ios.type_ == INCHI_IOS_TYPE_FILE as i32 && !ios.f.is_null() {
        let mut output_arguments_position = arguments.position;
        let mut output_arguments = FormatCursor {
            borrowed_string,
            arguments: arguments.arguments,
            position: &mut output_arguments_position,
        };
        return inchi_vfprintf_format_cursor(heap, ios.f, format, &mut output_arguments);
    }
    Ok(0)
}

pub(crate) fn source_vformat_format_arguments(
    heap: &mut SourceHeap,
    format: SourceFormat<'_>,
    arguments: &SourcePrintArguments<'_>,
) -> Result<Vec<u8>, SourceHeapError> {
    // INCHI source mol_fmt3.c:1032:
    //                 sprintf(szcoords, "%10g%10g%10g", fx, fy, fz);
    // INCHI✔️❗: Borrow fixed arguments with the same cursor start; no argument
    // allocation. Renderer/checked-access costs and defects remain independent.
    let mut position = arguments.position;
    let mut cursor = FormatCursor {
        borrowed_string: None,
        arguments: arguments.arguments,
        position: &mut position,
    };
    source_vformat_format_cursor(heap, format, &mut cursor)
}

pub(crate) fn source_vformat_format_strings<'a>(
    heap: &mut SourceHeap,
    format: SourceFormat<'a>,
    arguments: &SourcePrintArguments<'a>,
    strings: &'a [Option<SourceFormat<'a>>],
) -> Result<Vec<u8>, SourceHeapError> {
    // VERBATIM SOURCE FRAGMENT: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:1344..1346 C384 title sprintf
    //                 sprintf(szTitle, "%sInput Structure #%ld.%s%s%s%s%s",
    //                     bStructurePreprocessed ? "Preprocessed " : "",
    //                     num_inp, SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue), iINChI ? " (Reconnected)" : "");
    // END VERBATIM SOURCE FRAGMENT: C384 title sprintf
    // INCHI✔️❌: the existing borrowed-string cursor projects source varargs;
    // no new renderer/argument allocation/clone/scan. Keep initial position.
    // O(result bytes) output Vec/publication copy remains known native overhead
    // deferred in INCHI-PERF-PRINTF-FULL-OUTPUT; no no-worse claim.
    let mut position = arguments.position;
    let mut cursor = FormatCursor {
        borrowed_string: Some(BorrowedPrintStrings::Indexed(strings)),
        arguments: arguments.arguments,
        position: &mut position,
    };
    source_vformat_format_cursor(heap, format, &mut cursor)
}

fn source_vformat_bytes(
    heap: &mut SourceHeap,
    format: &[i8],
    arguments: &mut FormatCursor<'_>,
) -> Result<Vec<u8>, SourceHeapError> {
    // INCHI source ichi_io.c:1527..1528:
    //     my_va_start(argList, lpszFormat);
    //     ret = vsprintf(buf->pStr + buf->nUsedLength, lpszFormat, argList);
    // INCHI❗❌: One existing renderer for admitted bytes; full-output Vec adds
    // O(result bytes) scratch and a publication copy versus direct vsprintf.
    // HUMAN defers this cost in INCHI-PERF-PRINTF-FULL-OUTPUT; removing0136
    // is not cost repair/equivalence. Independent float/backend guards remain.
    // BEGIN SUPPLEMENTAL LIBC BYTE CONVERSIONS: third_party/musl/src/stdio/vfprintf.c591..604
    // INCHI✔️❗: 		narrow_c:
    // INCHI✔️❗: 		case 'c':
    // INCHI✔️❗: 			*(a=z-(p=1))=arg.i;
    // INCHI✔️❗: 			fl &= ~ZERO_PAD;
    // INCHI✔️❗: 			break;
    // INCHI✔️❗: 		case 'm':
    // INCHI✔️❗: 			if (1) a = strerror(errno); else
    // INCHI✔️❗: 		case 's':
    // INCHI✔️❗: 			a = arg.p ? arg.p : "(null)";
    // INCHI✔️❗: 			z = a + strnlen(a, p<0 ? INT_MAX : p);
    // INCHI✔️❗: 			if (p<0 && *z) goto overflow;
    // INCHI✔️❗: 			p = z-a;
    // INCHI✔️❗: 			fl &= ~ZERO_PAD;
    // INCHI✔️❗: 			break;
    // END SUPPLEMENTAL LIBC BYTE CONVERSIONS
    let mut output = Vec::with_capacity(format.len());
    let mut position = 0_usize;
    while position < format.len() {
        if format[position] as u8 != b'%' {
            output.push(format[position] as u8);
            position += 1;
            continue;
        }
        position += 1;
        if position < format.len() && format[position] as u8 == b'%' {
            output.push(b'%');
            position += 1;
            continue;
        }

        let mut alternate = false;
        let mut left = false;
        let mut plus = false;
        let mut space = false;
        let mut zero = false;
        while position < format.len() {
            match format[position] as u8 {
                b'#' => alternate = true,
                b'-' => left = true,
                b'+' => plus = true,
                b' ' => space = true,
                b'0' => zero = true,
                _ => break,
            }
            position += 1;
        }
        let mut width = 0_i32;
        if position < format.len() && format[position] as u8 == b'*' {
            width = next_format_int(arguments)?;
            position += 1;
            if width < 0 {
                left = true;
                width = width
                    .checked_abs()
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            }
        } else {
            while position < format.len() && (format[position] as u8).is_ascii_digit() {
                width = width
                    .checked_mul(10)
                    .and_then(|value| value.checked_add(i32::from(format[position] - b'0' as i8)))
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                position += 1;
            }
        }
        let mut precision = None;
        if position < format.len() && format[position] as u8 == b'.' {
            position += 1;
            let mut value = 0_i32;
            if position < format.len() && format[position] as u8 == b'*' {
                value = next_format_int(arguments)?;
                position += 1;
                if value >= 0 {
                    precision = Some(value as usize);
                }
            } else {
                while position < format.len() && (format[position] as u8).is_ascii_digit() {
                    value = value
                        .checked_mul(10)
                        .and_then(|current| {
                            current.checked_add(i32::from(format[position] - b'0' as i8))
                        })
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                    position += 1;
                }
                precision = Some(value as usize);
            }
        }
        let length =
            if position < format.len() && matches!(format[position] as u8, b'h' | b'l' | b'L') {
                let value = format[position] as u8;
                position += 1;
                if position < format.len() && format[position] as u8 == value && value != b'L' {
                    position += 1;
                    if value == b'h' { b'H' } else { b'q' }
                } else {
                    value
                }
            } else {
                0
            };
        let specifier = *format
            .get(position)
            .ok_or(SourceHeapError::UnsupportedSourceBehavior)? as u8;
        position += 1;

        let rendered = match specifier {
            b'c' => {
                let value = match next_format_argument(arguments)? {
                    SourceFormatArgument::Signed(value) => value as u8,
                    SourceFormatArgument::Byte(value) => value as u8,
                    _ => return Err(SourceHeapError::AllocationTypeMismatch),
                };
                // Narrow C char is one byte, including high-bit values and NUL.
                append_printf_bytes(&mut output, &[value], width, left)?;
                continue;
            }
            b's' => {
                let value = next_format_string(arguments)?;
                let (bytes, implicit_nul) = argument_string_bytes(heap, value)?;
                // strnlen stops at precision: a fully readable bounded prefix
                // needs no NUL beyond it. Implicit NUL belongs only to borrowed
                // content; real pointer strings preserve reached missing-NUL errors.
                let limit = precision.unwrap_or(usize::MAX).min(bytes.len());
                let selected = &bytes[..limit];
                let count = match memchr::memchr(0, bytemuck::cast_slice(selected)) {
                    Some(count) => count,
                    None if precision.is_some_and(|bound| bound <= bytes.len()) => limit,
                    None if implicit_nul => bytes.len(),
                    None => return Err(SourceHeapError::MissingNulTerminator),
                };
                append_printf_bytes(
                    &mut output,
                    bytemuck::cast_slice(&bytes[..count]),
                    width,
                    left,
                )?;
                continue;
            }
            b'd' | b'i' => {
                let raw = match next_format_argument(arguments)? {
                    SourceFormatArgument::Signed(value) => value,
                    _ => return Err(SourceHeapError::AllocationTypeMismatch),
                };
                let value = match length {
                    b'h' => i64::from(raw as i16),
                    b'H' => i64::from(raw as i8),
                    0 => i64::from(raw as i32),
                    _ => raw,
                };
                let mut scratch = [0_u8; 21];
                let signed = source_decimal_i64(value, None, &mut scratch);
                let digits = if value < 0 { &signed[1..] } else { signed };
                let digits = if precision == Some(0) && value == 0 {
                    &[]
                } else {
                    digits
                };
                let prefix: &[u8] = if value < 0 {
                    b"-"
                } else if plus {
                    b"+"
                } else if space {
                    b" "
                } else {
                    b""
                };
                append_printf_numeric(
                    &mut output,
                    prefix,
                    digits,
                    precision.unwrap_or(0),
                    width,
                    left,
                    zero && precision.is_none(),
                )?;
                continue;
            }
            b'u' | b'o' | b'x' | b'X' => {
                let raw = match next_format_argument(arguments)? {
                    SourceFormatArgument::Unsigned(value) => value,
                    _ => return Err(SourceHeapError::AllocationTypeMismatch),
                };
                let value = match length {
                    b'h' => u64::from(raw as u16),
                    b'H' => u64::from(raw as u8),
                    0 => u64::from(raw as u32),
                    _ => raw,
                };
                let mut scratch = [0_u8; 22];
                let digits = source_printf_unsigned(value, specifier, &mut scratch);
                let digits = if precision == Some(0) && value == 0 {
                    &[]
                } else {
                    digits
                };
                let mut minimum_digits = precision.unwrap_or(0);
                let prefix: &[u8] = if alternate && value != 0 && specifier == b'x' {
                    b"0x"
                } else if alternate && value != 0 && specifier == b'X' {
                    b"0X"
                } else {
                    b""
                };
                if alternate && specifier == b'o' && digits.first() != Some(&b'0') {
                    minimum_digits = minimum_digits.max(digits.len() + 1);
                }
                append_printf_numeric(
                    &mut output,
                    prefix,
                    digits,
                    minimum_digits,
                    width,
                    left,
                    zero && precision.is_none(),
                )?;
                continue;
            }
            b'e' | b'E' | b'f' | b'F' | b'g' | b'G' => {
                let value = match next_format_argument(arguments)? {
                    SourceFormatArgument::Float(value) => value,
                    _ => return Err(SourceHeapError::AllocationTypeMismatch),
                };
                let negative = value.is_sign_negative();
                // Supplemental musl vfprintf.c204 masks ZERO_PAD for special
                // floats. Signs remain independent, including signed NaNs.
                if !value.is_finite() {
                    zero = false;
                }
                let mut value = c_printf_float(value.abs(), specifier, precision, alternate);
                if negative {
                    value.insert(0, '-');
                } else if plus {
                    value.insert(0, '+');
                } else if space {
                    value.insert(0, ' ');
                }
                value
            }
            b'p' => {
                // Native ichi_io.c1289..1293 consumes a void*; c1528
                // delegates rendering to vsprintf. HUMAN 2026-10-17 selects
                // the actual Rust object's address for these diagnostics,
                // not another process's address or C object layout.
                // INCHI✔️✔️: Read address bits only, never dereference. NULL
                // remains (nil); nonNULL uses lowercase 0x hex. Existing
                // padding emits directly with O(1) stack scratch, no item
                // String, address map, allocation, scan or chemistry change.
                let pointer = match next_format_argument(arguments)? {
                    SourceFormatArgument::Pointer(pointer) => pointer,
                    _ => return Err(SourceHeapError::AllocationTypeMismatch),
                };
                if pointer.is_null() {
                    append_printf_bytes(&mut output, b"(nil)", width, left)?;
                } else {
                    let mut scratch = [0_u8; 22];
                    let digits = source_printf_unsigned(pointer.addr() as u64, b'x', &mut scratch);
                    append_printf_numeric(
                        &mut output,
                        b"0x",
                        digits,
                        precision.unwrap_or(1),
                        width,
                        left,
                        zero && precision.is_none(),
                    )?;
                }
                continue;
            }
            b'n' => {
                let pointer = match next_format_argument(arguments)? {
                    SourceFormatArgument::MutSigned(pointer) => pointer,
                    _ => return Err(SourceHeapError::AllocationTypeMismatch),
                };
                let count = i32::try_from(output.len())
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                *heap
                    .slice_mut(pointer)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = count;
                String::new()
            }
            _ => return Err(SourceHeapError::UnsupportedSourceBehavior),
        };
        // Remaining float/pointer/count paths preserve their existing item
        // semantics, but padding goes directly into the shared byte output.
        let bytes = rendered.as_bytes();
        let prefix_length = if bytes
            .first()
            .is_some_and(|byte| matches!(*byte, b'+' | b'-' | b' '))
        {
            1
        } else if bytes.starts_with(b"0x") || bytes.starts_with(b"0X") {
            2
        } else {
            0
        };
        append_printf_numeric(
            &mut output,
            &bytes[..prefix_length],
            &bytes[prefix_length..],
            0,
            width,
            left,
            zero,
        )?;
    }
    Ok(output)
}

fn source_printf_unsigned(value: u64, specifier: u8, scratch: &mut [u8; 22]) -> &[u8] {
    // BEGIN SUPPLEMENTAL LIBC INTEGER DIGITS: third_party/musl/src/stdio/vfprintf.c153..172
    // INCHI✔️❗: static char *fmt_x(uintmax_t x, char *s, int lower)
    // INCHI✔️❗: {
    // INCHI✔️❗: 	for (; x; x>>=4) *--s = xdigits[(x&15)]|lower;
    // INCHI✔️❗: 	return s;
    // INCHI✔️❗: }
    // INCHI✔️❗:
    // INCHI✔️❗: static char *fmt_o(uintmax_t x, char *s)
    // INCHI✔️❗: {
    // INCHI✔️❗: 	for (; x; x>>=3) *--s = '0' + (x&7);
    // INCHI✔️❗: 	return s;
    // INCHI✔️❗: }
    // INCHI✔️❗:
    // INCHI✔️❗: static char *fmt_u(uintmax_t x, char *s)
    // INCHI✔️❗: {
    // INCHI✔️❗: 	unsigned long y;
    // INCHI✔️❗: 	for (   ; x>ULONG_MAX; x/=10) *--s = '0' + x%10;
    // INCHI✔️❗: 	for (y=x;           y; y/=10) *--s = '0' + y%10;
    // INCHI✔️❗: 	return s;
    // INCHI✔️❗: }
    // INCHI✔️❗:
    // END SUPPLEMENTAL LIBC INTEGER DIGITS: third_party/musl/src/stdio/vfprintf.c153..172
    // INCHI✔️❗: Same backward unsigned radix conversion, full u64 range.
    // At most22octal/20decimal/16hex bytes; only initialized suffix escapes.
    // O(digits), fixed22bytes/no allocation; checked dispatch cost deferred.
    let mut value = value;
    let mut start = scratch.len();
    loop {
        let digit = match specifier {
            b'o' => {
                let digit = (value & 7) as u8;
                value >>= 3;
                digit
            }
            b'x' | b'X' => {
                let digit = (value & 15) as u8;
                value >>= 4;
                digit
            }
            _ => {
                let digit = (value % 10) as u8;
                value /= 10;
                digit
            }
        };
        start -= 1;
        scratch[start] = if digit < 10 {
            b'0' + digit
        } else if specifier == b'X' {
            b'A' + digit - 10
        } else {
            b'a' + digit - 10
        };
        if value == 0 {
            break;
        }
    }
    &scratch[start..]
}

#[allow(clippy::too_many_arguments)]
fn append_printf_numeric(
    output: &mut Vec<u8>,
    prefix: &[u8],
    digits: &[u8],
    minimum_digits: usize,
    width: i32,
    left: bool,
    zero: bool,
) -> Result<(), SourceHeapError> {
    // BEGIN SUPPLEMENTAL LIBC NUMBER PADDING: third_party/musl/src/stdio/vfprintf.c631..645
    // INCHI✔️❗:
    // INCHI✔️❗: 		if (p < z-a) p = z-a;
    // INCHI✔️❗: 		if (p > INT_MAX-pl) goto overflow;
    // INCHI✔️❗: 		if (w < pl+p) w = pl+p;
    // INCHI✔️❗: 		if (w > INT_MAX-cnt) goto overflow;
    // INCHI✔️❗:
    // INCHI✔️❗: 		pad(f, ' ', w, pl+p, fl);
    // INCHI✔️❗: 		out(f, prefix, pl);
    // INCHI✔️❗: 		pad(f, '0', w, pl+p, fl^ZERO_PAD);
    // INCHI✔️❗: 		pad(f, '0', p, z-a, 0);
    // INCHI✔️❗: 		out(f, a, z-a);
    // INCHI✔️❗: 		pad(f, ' ', w, pl+p, fl^LEFT_ADJ);
    // INCHI✔️❗:
    // INCHI✔️❗: 		l = w;
    // INCHI✔️❗: 	}
    // END SUPPLEMENTAL LIBC NUMBER PADDING: third_party/musl/src/stdio/vfprintf.c631..645
    // INCHI✔️❗: Prefix, precision zeros, digits and field padding are emitted
    // directly into the existing byte output, no item/padding allocation or
    // insertion copy. O(emitted bytes), O(1) scratch; full staging0136 cost is
    // deferred, independent float conversions remain guarded, whole-call OPEN.
    let precision_zeros = minimum_digits.saturating_sub(digits.len());
    let length = prefix
        .len()
        .checked_add(digits.len())
        .and_then(|length| length.checked_add(precision_zeros))
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    let padding = (width.max(0) as usize).saturating_sub(length);
    let end = output
        .len()
        .checked_add(length)
        .and_then(|length| length.checked_add(padding))
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    if !left && !zero {
        output.resize(output.len() + padding, b' ');
    }
    output.extend_from_slice(prefix);
    if !left && zero {
        output.resize(output.len() + padding, b'0');
    }
    output.resize(output.len() + precision_zeros, b'0');
    output.extend_from_slice(digits);
    if left {
        output.resize(end, b' ');
    }
    Ok(())
}

fn append_printf_bytes(
    output: &mut Vec<u8>,
    bytes: &[u8],
    width: i32,
    left: bool,
) -> Result<(), SourceHeapError> {
    // BEGIN SUPPLEMENTAL LIBC BYTE PADDING: third_party/musl/src/stdio/vfprintf.c638..643
    // INCHI✔️❗: 		out(f, prefix, pl);
    // INCHI✔️❗: 		pad(f, '0', w, pl+p, fl^ZERO_PAD);
    // INCHI✔️❗: 		pad(f, '0', p, z-a, 0);
    // INCHI✔️❗: 		out(f, a, z-a);
    // INCHI✔️❗: 		pad(f, ' ', w, pl+p, fl^LEFT_ADJ);
    // INCHI✔️❗:
    // END SUPPLEMENTAL LIBC BYTE PADDING
    // INCHI✔️❗: c/s clear ZERO_PAD; emit raw bytes and spaces directly into
    // the existing output. No item String/UTF8 conversion/padding Vec. O(width
    // or selected bytes), O(1) scratch; full-output staging0136 cost is deferred.
    let padding = (width.max(0) as usize).saturating_sub(bytes.len());
    let final_length = output
        .len()
        .checked_add(bytes.len())
        .and_then(|length| length.checked_add(padding))
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    if !left {
        output.resize(output.len() + padding, b' ');
    }
    output.extend_from_slice(bytes);
    if left {
        output.resize(final_length, b' ');
    }
    Ok(())
}

fn source_vformat_format_cursor(
    heap: &mut SourceHeap,
    format: SourceFormat<'_>,
    arguments: &mut FormatCursor<'_>,
) -> Result<Vec<u8>, SourceHeapError> {
    match format {
        SourceFormat::Pointer(pointer) => source_vformat_cursor(heap, pointer, arguments),
        SourceFormat::Bytes(bytes) => {
            let length = memchr::memchr(0, bytemuck::cast_slice(bytes))
                .ok_or(SourceHeapError::MissingNulTerminator)?;
            source_vformat_bytes(heap, &bytes[..length], arguments)
        }
        SourceFormat::Text(text) => {
            let bytes = text.as_bytes();
            let length = memchr::memchr(0, bytes).unwrap_or(bytes.len());
            source_vformat_bytes(heap, bytemuck::cast_slice(&bytes[..length]), arguments)
        }
    }
}

fn get_max_printf_length_format_cursor(
    heap: &SourceHeap,
    format: SourceFormat<'_>,
    arguments: &mut FormatCursor<'_>,
    literal_format: Option<&mut bool>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1065 GetMaxPrintfLength
    // INCHI✔️❌: int GetMaxPrintfLength(const char* lpszFormat, va_list argList)
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*ASSERT(AfxIsValidString(lpszFormat, FALSE));*/
    // INCHI✔️❌:     const char* lpsz;
    // INCHI✔️❌:     int nMaxLen, nWidth, nPrecision, nModifier, nItemLen;
    // INCHI✔️❌:
    // INCHI✔️❌:     nMaxLen = 0;
    // INCHI✔️❌:     /* make a guess at the maximum length of the resulting string */
    // INCHI✔️❌:     for (lpsz = lpszFormat; *lpsz; lpsz++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* moved from below for C syntax reason - 2024-09-01 DT */
    // INCHI✔️❌:         /* djb-rwth: return values needed for va_arg; djb-rwth: ignoring LLVM warnings: function returning value */
    // INCHI✔️❌:         int ivarg;
    // INCHI✔️❌:         double dvarg;
    // INCHI✔️❌:         void* ivvarg;
    // INCHI✔️❌:         int* ipvarg;
    // INCHI✔️❌:
    // INCHI✔️❌:         /* handle '%' character, but watch out for '%%' */
    // INCHI✔️❌:         if (*lpsz != '%' || *(++lpsz) == '%')
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nMaxLen += 1;
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         nItemLen = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         /*  handle '%' character with format */
    // INCHI✔️❌:         nWidth = 0;
    // INCHI✔️❌:         for (; *lpsz; lpsz++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* check for valid flags */
    // INCHI✔️❌:             if (*lpsz == '#')
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nMaxLen += 2;   /* for '0x' */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else if (*lpsz == '*')
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nWidth = va_arg(argList, int);
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else if (*lpsz == '-' || *lpsz == '+' || *lpsz == '0'
    // INCHI✔️❌:                 || *lpsz == ' ')
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else /* hit non-flag character */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* get width and skip it */
    // INCHI✔️❌:         if (nWidth == 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* width indicated by */
    // INCHI✔️❌:             nWidth = atoi(lpsz);
    // INCHI✔️❌:             for (; *lpsz && isdigit(*lpsz); lpsz++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /*ASSERT(nWidth >= 0);*/
    // INCHI✔️❌:         if (nWidth < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_error; /* instead of exception */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         nPrecision = 0;
    // INCHI✔️❌:         if (*lpsz == '.')
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* skip past '.' separator (width.precision)*/
    // INCHI✔️❌:             lpsz++;
    // INCHI✔️❌:
    // INCHI✔️❌:             /* get precision and skip it*/
    // INCHI✔️❌:             if (*lpsz == '*')
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nPrecision = va_arg(argList, int);
    // INCHI✔️❌:                 lpsz++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nPrecision = atoi(lpsz);
    // INCHI✔️❌:                 for (; *lpsz && isdigit(*lpsz); lpsz++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (nPrecision < 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto exit_error; /* instead of exception */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* should be on type modifier or specifier */
    // INCHI✔️❌:         nModifier = 0;
    // INCHI✔️❌:         switch (*lpsz)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* modifiers that affect size */
    // INCHI✔️❌:         case 'h':
    // INCHI✔️❌:             switch (lpsz[1])
    // INCHI✔️❌:             {
    // INCHI✔️❌:             case 'd':
    // INCHI✔️❌:             case 'i':
    // INCHI✔️❌:             case 'o':
    // INCHI✔️❌:             case 'x':
    // INCHI✔️❌:             case 'X':
    // INCHI✔️❌:             case 'u':
    // INCHI✔️❌:                 /* short unsigned, short double, etc. -- added to the original MS example */
    // INCHI✔️❌:                 /* ignore the fact that these modifiers do affect size */
    // INCHI✔️❌:                 lpsz++;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             default:
    // INCHI✔️❌:                 nModifier = FORCE_ANSI;
    // INCHI✔️❌:                 lpsz++;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         case 'l':
    // INCHI✔️❌:             switch (lpsz[1])
    // INCHI✔️❌:             {
    // INCHI✔️❌:             case 'd':
    // INCHI✔️❌:             case 'i':
    // INCHI✔️❌:             case 'o':
    // INCHI✔️❌:             case 'x':
    // INCHI✔️❌:             case 'X':
    // INCHI✔️❌:             case 'u':
    // INCHI✔️❌:             case 'f': /* long float -- post ANSI C */
    // INCHI✔️❌:                 /* long unsigned, long double, etc. -- added to the original MS example */
    // INCHI✔️❌:                 /* ignore the fact that these modifiers do affect size */
    // INCHI✔️❌:                 lpsz++;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             default:
    // INCHI✔️❌:                 /*
    // INCHI✔️❌:                 nModifier = FORCE_UNICODE;
    // INCHI✔️❌:                 lpsz ++;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:                 */
    // INCHI✔️❌:                 goto exit_error;  /* no UNICODE, please */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             break;
    // INCHI✔️❌:             /* modifiers that do not affect size */
    // INCHI✔️❌:         case 'F':
    // INCHI✔️❌:         case 'N':
    // INCHI✔️❌:         case 'L':
    // INCHI✔️❌:             lpsz++;
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* now should be on specifier */
    // INCHI✔️❌:
    // INCHI✔️❌:         switch (*lpsz | nModifier)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* single characters*/
    // INCHI✔️❌:         case 'c':
    // INCHI✔️❌:         case 'C':
    // INCHI✔️❌:             nItemLen = 2;
    // INCHI✔️❌:             ivarg = va_arg(argList, int); /* djb-rwth: int return value; ignoring LLVM warning */
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         case 'c' | FORCE_ANSI:
    // INCHI✔️❌:         case 'C' | FORCE_ANSI:
    // INCHI✔️❌:             nItemLen = 2;
    // INCHI✔️❌:             ivarg = va_arg(argList, int); /* djb-rwth: int return value; ignoring LLVM warning */
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         case 'c' | FORCE_UNICODE:
    // INCHI✔️❌:         case 'C' | FORCE_UNICODE:
    // INCHI✔️❌:             goto exit_error;  /* no UNICODE, please */
    // INCHI✔️❌:             /*
    // INCHI✔️❌:             nItemLen = 2;
    // INCHI✔️❌:             va_arg(argList, wchar_t);
    // INCHI✔️❌:             break;
    // INCHI✔️❌:             */
    // INCHI✔️❌:
    // INCHI✔️❌:             /* strings*/
    // INCHI✔️❌:         case 's':
    // INCHI✔️❌:         case 'S':
    // INCHI✔️❌:             nItemLen = (int)strlen(va_arg(argList, char*));
    // INCHI✔️❌:             nItemLen = inchi_max(1, nItemLen);
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         case 's' | FORCE_ANSI:
    // INCHI✔️❌:         case 'S' | FORCE_ANSI:
    // INCHI✔️❌:             nItemLen = (int)strlen(va_arg(argList, char*));
    // INCHI✔️❌:             nItemLen = inchi_max(1, nItemLen);
    // INCHI✔️❌:             break;
    // INCHI✔️❌:
    // INCHI✔️❌:         case 's' | FORCE_UNICODE:
    // INCHI✔️❌:         case 'S' | FORCE_UNICODE:
    // INCHI✔️❌:             goto exit_error;  /* no UNICODE, please */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* adjust nItemLen for strings */
    // INCHI✔️❌:         if (nItemLen != 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nItemLen = inchi_max(nItemLen, nWidth);
    // INCHI✔️❌:             if (nPrecision != 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nItemLen = inchi_min(nItemLen, nPrecision);
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             switch (*lpsz)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* integers */
    // INCHI✔️❌:             case 'd':
    // INCHI✔️❌:             case 'i':
    // INCHI✔️❌:             case 'u':
    // INCHI✔️❌:             case 'x':
    // INCHI✔️❌:             case 'X':
    // INCHI✔️❌:             case 'o':
    // INCHI✔️❌:                 ivarg = va_arg(argList, int); /* djb-rwth: int return value; ignoring LLVM warning */
    // INCHI✔️❌:                 nItemLen = 32;
    // INCHI✔️❌:                 nItemLen = inchi_max(nItemLen, nWidth + nPrecision);
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:
    // INCHI✔️❌:             case 'e':
    // INCHI✔️❌:             case 'f':
    // INCHI✔️❌:             case 'g':
    // INCHI✔️❌:             case 'G':
    // INCHI✔️❌:                 dvarg = va_arg(argList, double); /* djb-rwth: double return value; ignoring LLVM warning */
    // INCHI✔️❌:                 nItemLen = 32;
    // INCHI✔️❌:                 nItemLen = inchi_max(nItemLen, nWidth + nPrecision);
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:
    // INCHI✔️❌:             case 'p':
    // INCHI✔️❌:                 ivvarg = va_arg(argList, void*); /* djb-rwth: void* return value; ignoring LLVM warning */
    // INCHI✔️❌:                 nItemLen = 32;
    // INCHI✔️❌:                 nItemLen = inchi_max(nItemLen, nWidth + nPrecision);
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* no output */
    // INCHI✔️❌:             case 'n':
    // INCHI✔️❌:                 ipvarg = va_arg(argList, int*); /* djb-rwth: int* return value; ignoring LLVM warning */
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:
    // INCHI✔️❌:             default:
    // INCHI✔️❌:                 /*ASSERT(FALSE);*/  /* unknown formatting option*/
    // INCHI✔️❌:                 goto exit_error; /* instead of exception */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* adjust nMaxLen for output nItemLen */
    // INCHI✔️❌:         nMaxLen += nItemLen;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return nMaxLen;
    // INCHI✔️❌:
    // INCHI✔️❌: exit_error:
    // INCHI✔️❌:     return -1; /* wrong format */
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: GetMaxPrintfLength

    const FORCE_ANSI: i32 = 0x10000;
    let bytes = match format {
        SourceFormat::Pointer(pointer) => heap.slice(pointer)?,
        SourceFormat::Bytes(bytes) => bytes,
        SourceFormat::Text(text) => bytemuck::cast_slice(text.as_bytes()),
    };
    // INCHI source ichi_io.c:1073..1074:
    //     for (lpsz = lpszFormat; *lpsz; lpsz++)
    //     {
    // INCHI✔️❗: Read only reached format bytes, including source early-error
    // exits. No full suffix preflight, allocation or extra traversal. Missing
    // terminators remain structured errors when reached; Text retains exactly
    // its previous appended zero. Safe access/dispatch whole-call cost is OPEN.
    let byte_at = |position: usize| match bytes.get(position) {
        Some(byte) => Ok(*byte),
        None if matches!(format, SourceFormat::Text(_)) && position == bytes.len() => Ok(0),
        None => Err(SourceHeapError::MissingNulTerminator),
    };
    let mut position = 0_usize;
    let mut maximum_length = 0_i32;
    // INCHI source ichi_io.c:1082..1085:
    //         /* handle '%' character, but watch out for '%%' */
    //         if (*lpsz != '%' || *(++lpsz) == '%')
    //         {
    //             nMaxLen += 1;
    // INCHI✔️❗: Retain one bit from this existing traversal for the shared
    // direct literal consumer. No extra scan/allocation or argument access.
    // False on any reached percent, even an escape; all general formatting
    // semantics and early errors stay unchanged. Fixed dispatch cost OPEN.
    let mut has_percent = false;
    while byte_at(position)? != 0 {
        if byte_at(position)? != b'%' as i8 {
            maximum_length = maximum_length
                .checked_add(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            position += 1;
            continue;
        }
        has_percent = true;
        position += 1;
        if byte_at(position)? == b'%' as i8 {
            maximum_length = maximum_length
                .checked_add(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            position += 1;
            continue;
        }

        let mut item_length = 0_i32;
        let mut width = 0_i32;
        while byte_at(position)? != 0 {
            match byte_at(position)? as u8 {
                b'#' => {
                    maximum_length = maximum_length
                        .checked_add(2)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                }
                b'*' => width = next_format_int(arguments)?,
                b'-' | b'+' | b'0' | b' ' => {}
                _ => break,
            }
            position += 1;
        }
        if width == 0 {
            while (byte_at(position)? as u8).is_ascii_digit() {
                let digit = i32::from(byte_at(position)? - b'0' as i8);
                width = width
                    .checked_mul(10)
                    .and_then(|value| value.checked_add(digit))
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                position += 1;
            }
        }
        if width < 0 {
            return Ok(-1);
        }

        let mut precision = 0_i32;
        if byte_at(position)? == b'.' as i8 {
            position += 1;
            if byte_at(position)? == b'*' as i8 {
                precision = next_format_int(arguments)?;
                position += 1;
            } else {
                while (byte_at(position)? as u8).is_ascii_digit() {
                    let digit = i32::from(byte_at(position)? - b'0' as i8);
                    precision = precision
                        .checked_mul(10)
                        .and_then(|value| value.checked_add(digit))
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                    position += 1;
                }
            }
            if precision < 0 {
                return Ok(-1);
            }
        }

        let mut modifier = 0_i32;
        match byte_at(position)? as u8 {
            b'h' => {
                position += 1;
                if !matches!(
                    byte_at(position)? as u8,
                    b'd' | b'i' | b'o' | b'x' | b'X' | b'u'
                ) {
                    modifier = FORCE_ANSI;
                }
            }
            b'l' => {
                position += 1;
                if !matches!(
                    byte_at(position)? as u8,
                    b'd' | b'i' | b'o' | b'x' | b'X' | b'u' | b'f'
                ) {
                    return Ok(-1);
                }
            }
            b'F' | b'N' | b'L' => position += 1,
            _ => {}
        }

        let specifier = byte_at(position)? as u8;
        match (i32::from(specifier) | modifier) as u32 {
            value
                if value == b'c' as u32
                    || value == b'C' as u32
                    || value == (i32::from(b'c') | FORCE_ANSI) as u32
                    || value == (i32::from(b'C') | FORCE_ANSI) as u32 =>
            {
                item_length = 2;
                next_format_int(arguments)?;
            }
            value
                if value == b's' as u32
                    || value == b'S' as u32
                    || value == (i32::from(b's') | FORCE_ANSI) as u32
                    || value == (i32::from(b'S') | FORCE_ANSI) as u32 =>
            {
                let value = next_format_string(arguments)?;
                let (string, implicit_nul) = argument_string_bytes(heap, value)?;
                let length = string
                    .iter()
                    .position(|byte| *byte == 0)
                    .or_else(|| implicit_nul.then_some(string.len()))
                    .ok_or(SourceHeapError::MissingNulTerminator)?;
                // INCHI source ichi_io.c:1237..1238:
                //             nItemLen = (int)strlen(va_arg(argList, char*));
                //             nItemLen = inchi_max(1, nItemLen);
                // INCHI✔️❗: GCC/Linux LP64 signed conversion retains the low
                // 32 bits before max(1), including lengths above INT_MAX.
                // O(1) narrowing; source-required strlen remains one traversal.
                item_length = (length as i32).max(1);
            }
            _ => {}
        }

        if item_length != 0 {
            item_length = item_length.max(width);
            if precision != 0 {
                item_length = item_length.min(precision);
            }
        } else {
            match specifier {
                b'd' | b'i' | b'u' | b'x' | b'X' | b'o' => {
                    // INCHI source ichi_io.c:1277..1279:
                    //                 ivarg = va_arg(argList, int); /* djb-rwth: int return value; ignoring LLVM warning */
                    //                 nItemLen = 32;
                    //                 nItemLen = inchi_max(nItemLen, nWidth + nPrecision);
                    // INCHI✔️❗: Consume the integer slot once; its value is
                    // discarded even for source long modifiers. No narrowing
                    // admission, allocation or scan; tagged ABI cost stays OPEN.
                    match next_format_argument(arguments)? {
                        SourceFormatArgument::Signed(_)
                        | SourceFormatArgument::Unsigned(_)
                        | SourceFormatArgument::Byte(_) => {}
                        _ => return Err(SourceHeapError::AllocationTypeMismatch),
                    }
                    item_length = 32_i32.max(
                        width
                            .checked_add(precision)
                            .ok_or(SourceHeapError::SourceIntegerOverflow)?,
                    );
                }
                b'e' | b'f' | b'g' | b'G' => {
                    match next_format_argument(arguments)? {
                        SourceFormatArgument::Float(_) => {}
                        _ => return Err(SourceHeapError::AllocationTypeMismatch),
                    }
                    item_length = 32_i32.max(
                        width
                            .checked_add(precision)
                            .ok_or(SourceHeapError::SourceIntegerOverflow)?,
                    );
                }
                b'p' => {
                    match next_format_argument(arguments)? {
                        SourceFormatArgument::Pointer(_) => {}
                        _ => return Err(SourceHeapError::AllocationTypeMismatch),
                    }
                    item_length = 32_i32.max(
                        width
                            .checked_add(precision)
                            .ok_or(SourceHeapError::SourceIntegerOverflow)?,
                    );
                }
                b'n' => match next_format_argument(arguments)? {
                    // INCHI source ichi_io.c:1297..1299:
                    //             case 'n':
                    //                 ipvarg = va_arg(argList, int*); /* djb-rwth: int* return value; ignoring LLVM warning */
                    //                 break;
                    // INCHI✔️❗: Consume the same typed pointer once, without
                    // dereferencing or writing. The renderer accepts MutSigned.
                    SourceFormatArgument::Pointer(_) | SourceFormatArgument::MutSigned(_) => {}
                    _ => return Err(SourceHeapError::AllocationTypeMismatch),
                },
                _ => return Ok(-1),
            }
        }
        maximum_length = maximum_length
            .checked_add(item_length)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        position += 1;
    }
    if let Some(literal_format) = literal_format {
        *literal_format = !has_percent;
    }
    Ok(maximum_length)
}

fn inchi_vfprintf_format_cursor(
    heap: &mut SourceHeap,
    file: SourceMutPointer<FILE>,
    format: SourceFormat<'_>,
    arguments: &mut FormatCursor<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:828 inchi_vfprintf
    // INCHI✔️❌: int inchi_vfprintf(FILE* f, const char* lpszFormat, va_list argList)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0;
    // INCHI✔️❌:     if (lpszFormat && lpszFormat[0]) /* djb-rwth: condition added as lpszFormat == 0 may lead to undefined ret value */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (f == stderr && '\r' == lpszFormat[strlen(lpszFormat) - 1])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = vfprintf(f, lpszFormat, argList);
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = vfprintf(f, lpszFormat, argList);
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_vfprintf

    let first = match format {
        SourceFormat::Pointer(pointer) => {
            if pointer.is_null() {
                return Ok(0);
            }
            heap.slice(pointer)?.first().copied()
        }
        SourceFormat::Bytes(bytes) => bytes.first().copied(),
        SourceFormat::Text(text) => Some(text.as_bytes().first().copied().unwrap_or(0) as i8),
    };
    if first == Some(0) {
        return Ok(0);
    }
    let output = source_vformat_format_cursor(heap, format, arguments)?;
    let file = heap
        .slice_mut(file)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if file.error {
        // musl __stdio_write.c: f->flags |= F_ERR;
        // Sticky state records failure; it does not cause later failures.
        file.error_indicator = true;
        return Ok(-1);
    }
    let position =
        usize::try_from(file.position).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let end = position
        .checked_add(output.len())
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    if file.bytes.len() < end {
        file.bytes.resize(end, 0);
    }
    file.bytes[position..end].copy_from_slice(&output);
    file.position = u64::try_from(end).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    i32::try_from(output.len()).map_err(|_| SourceHeapError::SourceIntegerOverflow)
}

pub(crate) fn inchi_strbuf_printf_format(
    heap: &mut SourceHeap,
    buf: Option<&mut INCHI_IOS_STRING>,
    format: SourceFormat<'_>,
    arguments: &SourceVaList,
) -> Result<i32, SourceHeapError> {
    // INCHI✔️❗: Private borrowed format projection; argument order and initial
    // cursor are unchanged. The canonical core retains the full native frame.
    inchi_strbuf_printf_format_arguments(
        heap,
        buf,
        format,
        &SourcePrintArguments {
            arguments: &arguments.arguments,
            position: arguments.position,
        },
    )
}

pub(crate) fn inchi_ios_print_nodisplay_format(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    format: SourceFormat<'_>,
    arguments: &SourceVaList,
) -> Result<i32, SourceHeapError> {
    // INCHI✔️❗: Private borrowed format projection; argument order and initial
    // cursor are unchanged. The canonical core retains the full native frame.
    inchi_ios_print_nodisplay_format_arguments(
        heap,
        ios,
        stdout,
        format,
        &SourcePrintArguments {
            arguments: &arguments.arguments,
            position: arguments.position,
        },
    )
}

pub(crate) fn inchi_ios_eprint_format(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    format: SourceFormat<'_>,
    arguments: &SourceVaList,
) -> Result<i32, SourceHeapError> {
    // INCHI✔️❗: Private borrowed format projection; argument order and initial
    // cursor are unchanged. The canonical core retains the full native frame.
    inchi_ios_eprint_format_arguments(
        heap,
        ios,
        format,
        &SourcePrintArguments {
            arguments: &arguments.arguments,
            position: arguments.position,
        },
    )
}

// Source va_start/va_arg cursor state; the argument storage is borrowed, never copied.
// A reference to position preserves partial advancement on errors and unwinding.
/// Borrowed arguments for one private printf call; no owned argument buffer.
#[derive(Clone, Copy)]
pub(crate) struct SourcePrintArguments<'a> {
    pub(crate) arguments: &'a [SourceFormatArgument],
    pub(crate) position: u64,
}

#[derive(Clone, Copy)]
enum BorrowedFormatString<'a> {
    Format(SourceFormat<'a>),
    Content(&'a [u8]),
}

#[derive(Clone, Copy)]
enum BorrowedPrintStrings<'a> {
    One(usize, BorrowedFormatString<'a>),
    Indexed(&'a [Option<SourceFormat<'a>>]),
}

fn borrowed_print_string<'a>(
    strings: Option<BorrowedPrintStrings<'a>>,
    index: usize,
) -> Option<BorrowedFormatString<'a>> {
    // INCHI source ichi_io.c:1237:
    //             nItemLen = (int)strlen(va_arg(argList, char*));
    // INCHI✔️❗: Concrete source argument slots, direct O(1) selection. Keep
    // native pointer or explicit/implicit NUL views; no scan/allocation/clone.
    match strings? {
        BorrowedPrintStrings::One(slot, value) => (slot == index).then_some(value),
        BorrowedPrintStrings::Indexed(values) => values
            .get(index)
            .copied()
            .flatten()
            .map(BorrowedFormatString::Format),
    }
}

struct FormatCursor<'a> {
    borrowed_string: Option<BorrowedPrintStrings<'a>>,
    arguments: &'a [SourceFormatArgument],
    position: &'a mut u64,
}

fn next_format_argument(
    arguments: &mut FormatCursor<'_>,
) -> Result<SourceFormatArgument, SourceHeapError> {
    // INCHI source ichi_io.c:1277:
    //                 ivarg = va_arg(argList, int); /* djb-rwth: int return value; ignoring LLVM warning */
    // INCHI✔️❗: O(1) tagged selection/advance without allocation; checked
    // bounds, tag footprint and whole ABI-relative cost remain deferred.
    let index =
        usize::try_from(*arguments.position).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if borrowed_print_string(arguments.borrowed_string, index).is_some() {
        *arguments.position = arguments
            .position
            .checked_add(1)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        return Err(SourceHeapError::AllocationTypeMismatch);
    }
    let argument = arguments
        .arguments
        .get(index)
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    *arguments.position = arguments
        .position
        .checked_add(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    Ok(argument)
}

fn next_format_string<'a>(
    arguments: &mut FormatCursor<'a>,
) -> Result<BorrowedFormatString<'a>, SourceHeapError> {
    // INCHI source ichi_io.c:1237:
    //             nItemLen = (int)strlen(va_arg(argList, char*));
    // INCHI✔️❗: Borrowed %-s argument uses the same reached selection/advance.
    // Pointer callers keep their existing tag errors after advancement.
    if arguments.borrowed_string.is_some() {
        let index = usize::try_from(*arguments.position)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if let Some(value) = borrowed_print_string(arguments.borrowed_string, index) {
            *arguments.position = arguments
                .position
                .checked_add(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            return Ok(value);
        }
    }
    match next_format_argument(arguments)? {
        SourceFormatArgument::Bytes(pointer) => {
            Ok(BorrowedFormatString::Format(SourceFormat::Pointer(pointer)))
        }
        _ => Err(SourceHeapError::AllocationTypeMismatch),
    }
}

fn argument_string_bytes<'a>(
    heap: &'a SourceHeap,
    value: BorrowedFormatString<'a>,
) -> Result<(&'a [i8], bool), SourceHeapError> {
    // INCHI source ichi_io.c1237: nItemLen = (int)strlen(va_arg(argList, char*));
    // INCHI✔️❗: Concrete borrowed content replaces the caller's temporary
    // copy plus explicit trailing zero. Embedded NUL wins; no UTF8 admission,
    // allocation, copy or extra scan. O(1) projection; whole-call cost OPEN.
    match value {
        BorrowedFormatString::Format(format) => format_string_bytes(heap, format),
        BorrowedFormatString::Content(bytes) => Ok((bytemuck::cast_slice(bytes), true)),
    }
}

fn format_string_bytes<'a>(
    heap: &'a SourceHeap,
    value: SourceFormat<'a>,
) -> Result<(&'a [i8], bool), SourceHeapError> {
    // INCHI source ichiprt2.c:536:
    //     n = inchi_strbuf_printf( buf, "%-s", szTailingDelim );
    // INCHI✔️❗: Text replaces the former byte copy plus explicit NUL.
    // First embedded NUL and trailing bytes stay unchanged; no scan or copy.
    match value {
        SourceFormat::Pointer(pointer) => Ok((heap.slice(pointer)?, false)),
        SourceFormat::Bytes(bytes) => Ok((bytes, false)),
        SourceFormat::Text(text) => Ok((bytemuck::cast_slice(text.as_bytes()), true)),
    }
}

fn next_format_int(arguments: &mut FormatCursor<'_>) -> Result<i32, SourceHeapError> {
    // INCHI source ichi_io.c:1103:
    //                 nWidth = va_arg(argList, int);
    // INCHI✔️❗: Preserve typed/range errors after the same cursor advance.
    match next_format_argument(arguments)? {
        SourceFormatArgument::Signed(value) => {
            i32::try_from(value).map_err(|_| SourceHeapError::SourceIntegerOverflow)
        }
        SourceFormatArgument::Unsigned(value) => Ok(value as u32 as i32),
        SourceFormatArgument::Byte(value) => Ok(i32::from(value)),
        _ => Err(SourceHeapError::AllocationTypeMismatch),
    }
}

fn c_printf_exponential(value: f64, precision: usize, upper: bool) -> String {
    // BEGIN SUPPLEMENTAL LIBC FLOAT: third_party/musl/src/stdio/vfprintf.c370..381
    // INCHI❗❗: 	} else {
    // INCHI❗❗: 		estr=fmt_u(e<0 ? -e : e, ebuf);
    // INCHI❗❗: 		while(ebuf-estr<2) *--estr='0';
    // INCHI❗❗: 		*--estr = (e<0 ? '-' : '+');
    // INCHI❗❗: 		*--estr = t;
    // INCHI❗❗: 		if (ebuf-estr > INT_MAX-l) return -1;
    // INCHI❗❗: 		l += ebuf-estr;
    // INCHI❗❗: 	}
    // INCHI❗❗:
    // INCHI❗❗: 	if (l > INT_MAX-pl) return -1;
    // INCHI❗❗: 	pad(f, ' ', w, pl+l, fl);
    // INCHI❗❗: 	out(f, prefix, pl);
    // END SUPPLEMENTAL LIBC FLOAT
    // INCHI❗❗: Retain existing Rust rounding, normalize its exponent in place.
    // One output String, no second mantissa String or numeric re-conversion.
    let mut result = format!("{value:.precision$e}");
    let exponent = take_printf_exponent(&mut result);
    push_printf_exponent(&mut result, exponent, upper);
    result
}

fn take_printf_exponent(value: &mut String) -> i32 {
    // INCHI source ichi_io.c:1528:
    //     ret = vsprintf(buf->pStr + buf->nUsedLength, lpszFormat, argList);
    // INCHI❗❗: Decode the already-rounded Rust scientific exponent only;
    // finite binary64 gives at most three decimal exponent digits. O(1), no
    // allocation. This is not parsing the mantissa back into floating point.
    let marker = value.rfind('e').expect("scientific format has exponent");
    let bytes = &value.as_bytes()[marker + 1..];
    let negative = bytes.first() == Some(&b'-');
    let start = usize::from(matches!(bytes.first(), Some(b'-' | b'+')));
    let mut exponent = 0_i32;
    for &digit in &bytes[start..] {
        debug_assert!(digit.is_ascii_digit());
        exponent = exponent * 10 + i32::from(digit - b'0');
    }
    value.truncate(marker);
    if negative { -exponent } else { exponent }
}

fn push_printf_exponent(value: &mut String, exponent: i32, upper: bool) {
    // BEGIN SUPPLEMENTAL LIBC FLOAT: third_party/musl/src/stdio/vfprintf.c371..380
    // INCHI❗❗: 		estr=fmt_u(e<0 ? -e : e, ebuf);
    // INCHI❗❗: 		while(ebuf-estr<2) *--estr='0';
    // INCHI❗❗: 		*--estr = (e<0 ? '-' : '+');
    // INCHI❗❗: 		*--estr = t;
    // INCHI❗❗: 		if (ebuf-estr > INT_MAX-l) return -1;
    // INCHI❗❗: 		l += ebuf-estr;
    // INCHI❗❗: 	}
    // INCHI❗❗:
    // INCHI❗❗: 	if (l > INT_MAX-pl) return -1;
    // INCHI❗❗: 	pad(f, ' ', w, pl+l, fl);
    // END SUPPLEMENTAL LIBC FLOAT
    // INCHI✔️❗: Explicit sign and at least two exponent digits, fixed stack
    // scratch appended to the existing output. No owned suffix or format!.
    value.push(if upper { 'E' } else { 'e' });
    value.push(if exponent < 0 { '-' } else { '+' });
    let mut scratch = [0_u8; 21];
    let digits = source_decimal_i64(i64::from(exponent.abs()), None, &mut scratch);
    if digits.len() < 2 {
        value.push('0');
    }
    value.push_str(std::str::from_utf8(digits).expect("decimal digits are ASCII"));
}

fn c_printf_float(value: f64, specifier: u8, precision: Option<usize>, alternate: bool) -> String {
    // BEGIN SUPPLEMENTAL LIBC FLOAT: third_party/musl/src/stdio/vfprintf.c200..208
    // INCHI❗❗: 	if (!isfinite(y)) {
    // INCHI❗❗: 		char *s = (t&32)?"inf":"INF";
    // INCHI❗❗: 		if (y!=y) s=(t&32)?"nan":"NAN";
    // INCHI❗❗: 		pad(f, ' ', w, 3+pl, fl&~ZERO_PAD);
    // INCHI❗❗: 		out(f, prefix, pl);
    // INCHI❗❗: 		out(f, s, 3);
    // INCHI❗❗: 		pad(f, ' ', w, 3+pl, fl^LEFT_ADJ);
    // INCHI❗❗: 		return MAX(w, 3+pl);
    // INCHI❗❗: 	}
    // END SUPPLEMENTAL LIBC FLOAT
    // BEGIN SUPPLEMENTAL LIBC FLOAT: third_party/musl/src/stdio/vfprintf.c307..367
    // INCHI❗❗: 	else e=0;
    // INCHI❗❗:
    // INCHI❗❗: 	/* Perform rounding: j is precision after the radix (possibly neg) */
    // INCHI❗❗: 	j = p - ((t|32)!='f')*e - ((t|32)=='g' && p);
    // INCHI❗❗: 	if (j < 9*(z-r-1)) {
    // INCHI❗❗: 		uint32_t x;
    // INCHI❗❗: 		/* We avoid C's broken division of negative numbers */
    // INCHI❗❗: 		d = r + 1 + ((j+9*LDBL_MAX_EXP)/9 - LDBL_MAX_EXP);
    // INCHI❗❗: 		j += 9*LDBL_MAX_EXP;
    // INCHI❗❗: 		j %= 9;
    // INCHI❗❗: 		for (i=10, j++; j<9; i*=10, j++);
    // INCHI❗❗: 		x = *d % i;
    // INCHI❗❗: 		/* Are there any significant digits past j? */
    // INCHI❗❗: 		if (x || d+1!=z) {
    // INCHI❗❗: 			long double round = 2/LDBL_EPSILON;
    // INCHI❗❗: 			long double small;
    // INCHI❗❗: 			if ((*d/i & 1) || (i==1000000000 && d>a && (d[-1]&1)))
    // INCHI❗❗: 				round += 2;
    // INCHI❗❗: 			if (x<i/2) small=0x0.8p0;
    // INCHI❗❗: 			else if (x==i/2 && d+1==z) small=0x1.0p0;
    // INCHI❗❗: 			else small=0x1.8p0;
    // INCHI❗❗: 			if (pl && *prefix=='-') round*=-1, small*=-1;
    // INCHI❗❗: 			*d -= x;
    // INCHI❗❗: 			/* Decide whether to round by probing round+small */
    // INCHI❗❗: 			if (round+small != round) {
    // INCHI❗❗: 				*d = *d + i;
    // INCHI❗❗: 				while (*d > 999999999) {
    // INCHI❗❗: 					*d--=0;
    // INCHI❗❗: 					if (d<a) *--a=0;
    // INCHI❗❗: 					(*d)++;
    // INCHI❗❗: 				}
    // INCHI❗❗: 				for (i=10, e=9*(r-a); *a>=i; i*=10, e++);
    // INCHI❗❗: 			}
    // INCHI❗❗: 		}
    // INCHI❗❗: 		if (z>d+1) z=d+1;
    // INCHI❗❗: 	}
    // INCHI❗❗: 	for (; z>a && !z[-1]; z--);
    // INCHI❗❗:
    // INCHI❗❗: 	if ((t|32)=='g') {
    // INCHI❗❗: 		if (!p) p++;
    // INCHI❗❗: 		if (p>e && e>=-4) {
    // INCHI❗❗: 			t--;
    // INCHI❗❗: 			p-=e+1;
    // INCHI❗❗: 		} else {
    // INCHI❗❗: 			t-=2;
    // INCHI❗❗: 			p--;
    // INCHI❗❗: 		}
    // INCHI❗❗: 		if (!(fl&ALT_FORM)) {
    // INCHI❗❗: 			/* Count trailing zeros in last place */
    // INCHI❗❗: 			if (z>a && z[-1]) for (i=10, j=0; z[-1]%i==0; i*=10, j++);
    // INCHI❗❗: 			else j=9;
    // INCHI❗❗: 			if ((t|32)=='f')
    // INCHI❗❗: 				p = MIN(p,MAX(0,9*(z-r-1)-j));
    // INCHI❗❗: 			else
    // INCHI❗❗: 				p = MIN(p,MAX(0,9*(z-r-1)+e-j));
    // INCHI❗❗: 		}
    // INCHI❗❗: 	}
    // INCHI❗❗: 	if (p > INT_MAX-1-(p || (fl&ALT_FORM)))
    // INCHI❗❗: 		return -1;
    // INCHI❗❗: 	l = 1 + p + (p || (fl&ALT_FORM));
    // INCHI❗❗: 	if ((t|32)=='f') {
    // END SUPPLEMENTAL LIBC FLOAT
    // INCHI❗❌: General style uses the rounded exponent and same rounded
    // digits, never log10 or a second floating conversion. One materialized
    // item plus shifts/copy remains extra work versus native direct emission;
    // INCHI-PERF-PRINTF-FLOAT deferred, not performance equivalence. Existing
    // Rust rounding retained; configured libc/fenv comparison remains OPEN.
    let upper = matches!(specifier, b'E' | b'F' | b'G');
    if value.is_nan() {
        return if upper { "NAN" } else { "nan" }.to_owned();
    }
    if value.is_infinite() {
        return if upper { "INF" } else { "inf" }.to_owned();
    }
    match specifier {
        b'f' | b'F' => {
            let precision = precision.unwrap_or(6);
            let mut result = format!("{value:.precision$}");
            if alternate && precision == 0 {
                result.push('.');
            }
            result
        }
        b'e' | b'E' => {
            let precision = precision.unwrap_or(6);
            let mut result = c_printf_exponential(value, precision, upper);
            if alternate && precision == 0 {
                result.insert(1, '.');
            }
            result
        }
        b'g' | b'G' => {
            let significant = precision.unwrap_or(6).max(1);
            // Binary64=M*2^E with M<10^16 and E>=-1074. For E<0,
            // M*5^(-E)<10^768; E>=0 needs at most309decimal digits.
            // Above768significant digits the representation is exact: extra
            // zeros would only be trimmed. Alternate requires those zeros.
            let digits = if alternate {
                significant
            } else {
                significant.min(768)
            };
            let fractional = digits - 1;
            let mut result = format!("{value:.fractional$e}");
            let exponent = take_printf_exponent(&mut result);
            if fractional != 0 {
                result.remove(1);
            }
            if !alternate {
                while result.len() > 1 && result.ends_with('0') {
                    result.pop();
                }
            }
            if exponent >= -4 && (exponent < 0 || (exponent as usize) < significant) {
                let point = exponent + 1;
                if point <= 0 {
                    let prefix = match point {
                        -3 => "0.000",
                        -2 => "0.00",
                        -1 => "0.0",
                        0 => "0.",
                        _ => unreachable!("fixed general exponent is at least -4"),
                    };
                    result.insert_str(0, prefix);
                } else {
                    let point = point as usize;
                    while result.len() < point {
                        result.push('0');
                    }
                    if point < result.len() {
                        result.insert(point, '.');
                    } else if alternate {
                        result.push('.');
                    }
                }
            } else {
                if result.len() > 1 || alternate {
                    result.insert(1, '.');
                }
                push_printf_exponent(&mut result, exponent, upper);
            }
            result
        }
        _ => unreachable!(),
    }
}

pub(crate) fn source_vformat(
    heap: &mut SourceHeap,
    format: SourceConstPointer<i8>,
    arguments: &mut SourceVaList,
) -> Result<Vec<u8>, SourceHeapError> {
    let mut cursor = FormatCursor {
        borrowed_string: None,
        arguments: &arguments.arguments,
        position: &mut arguments.position,
    };
    source_vformat_cursor(heap, format, &mut cursor)
}

pub(crate) fn source_vformat_arguments(
    heap: &mut SourceHeap,
    format: SourceConstPointer<i8>,
    arguments: &SourcePrintArguments<'_>,
) -> Result<Vec<u8>, SourceHeapError> {
    // INCHI source mol_fmt3.c:1032:
    //                 sprintf(szcoords, "%10g%10g%10g", fx, fy, fz);
    // INCHI✔️❗: Borrow fixed arguments with the same cursor start; no argument
    // allocation. Renderer/checked-access costs and defects remain independent.
    source_vformat_format_arguments(heap, SourceFormat::Pointer(format), arguments)
}

fn source_vformat_cursor(
    heap: &mut SourceHeap,
    format: SourceConstPointer<i8>,
    arguments: &mut FormatCursor<'_>,
) -> Result<Vec<u8>, SourceHeapError> {
    let format_values = heap.slice(format)?;
    let format_length = memchr::memchr(0, bytemuck::cast_slice(format_values))
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    // SAFETY: formatting never frees or resizes the format allocation. `%n`
    // writes through an i32 allocation and therefore cannot overlap this i8
    // view under SourceHeap's allocation type checks.
    let format_view = unsafe { heap.stable_slice(format)? };
    let format = format_view.prefix(format_length)?;
    source_vformat_bytes(heap, format, arguments)
}

pub(crate) fn inchi_vfprintf(
    heap: &mut SourceHeap,
    file: SourceMutPointer<FILE>,
    format: SourceConstPointer<i8>,
    arguments: &mut SourceVaList,
) -> Result<i32, SourceHeapError> {
    let mut cursor = FormatCursor {
        borrowed_string: None,
        arguments: &arguments.arguments,
        position: &mut arguments.position,
    };
    inchi_vfprintf_cursor(heap, file, format, &mut cursor)
}

fn inchi_vfprintf_cursor(
    heap: &mut SourceHeap,
    file: SourceMutPointer<FILE>,
    format: SourceConstPointer<i8>,
    arguments: &mut FormatCursor<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:828 inchi_vfprintf
    // INCHI✔️❌: int inchi_vfprintf(FILE* f, const char* lpszFormat, va_list argList)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0;
    // INCHI✔️❌:     if (lpszFormat && lpszFormat[0]) /* djb-rwth: condition added as lpszFormat == 0 may lead to undefined ret value */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (f == stderr && '\r' == lpszFormat[strlen(lpszFormat) - 1])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = vfprintf(f, lpszFormat, argList);
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = vfprintf(f, lpszFormat, argList);
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_vfprintf

    inchi_vfprintf_format_cursor(heap, file, SourceFormat::Pointer(format), arguments)
}

pub(crate) fn inchi_print_nodisplay(
    heap: &mut SourceHeap,
    file: Option<SourceMutPointer<FILE>>,
    stdout: SourceMutPointer<FILE>,
    format: SourceConstPointer<i8>,
    arguments: &SourceVaList,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:867 inchi_print_nodisplay
    // INCHI✔️❌: int inchi_print_nodisplay(FILE* f, const char* lpszFormat, ...)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0;
    // INCHI✔️❌:     va_list argList;
    // INCHI✔️❌:     FILE* fi;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (f)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         fi = f;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         fi = stdout;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     my_va_start(argList, lpszFormat);
    // INCHI✔️❌:     ret = vfprintf(fi, lpszFormat, argList);
    // INCHI✔️❌:     va_end(argList);
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_print_nodisplay

    inchi_print_nodisplay_arguments(
        heap,
        file,
        stdout,
        format,
        &SourcePrintArguments {
            arguments: &arguments.arguments,
            position: arguments.position,
        },
    )
}

pub(crate) fn inchi_print_nodisplay_arguments(
    heap: &mut SourceHeap,
    file: Option<SourceMutPointer<FILE>>,
    stdout: SourceMutPointer<FILE>,
    format: SourceConstPointer<i8>,
    arguments: &SourcePrintArguments<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:867 inchi_print_nodisplay
    // INCHI✔️❌: int inchi_print_nodisplay(FILE* f, const char* lpszFormat, ...)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0;
    // INCHI✔️❌:     va_list argList;
    // INCHI✔️❌:     FILE* fi;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (f)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         fi = f;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         fi = stdout;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     my_va_start(argList, lpszFormat);
    // INCHI✔️❌:     ret = vfprintf(fi, lpszFormat, argList);
    // INCHI✔️❌:     va_end(argList);
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_print_nodisplay

    inchi_print_nodisplay_format_arguments(
        heap,
        file,
        stdout,
        SourceFormat::Pointer(format),
        arguments,
    )
}

pub(crate) fn inchi_ios_print(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    format: SourceConstPointer<i8>,
    arguments: &SourceVaList,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:477 inchi_ios_print
    // INCHI✔❌: int inchi_ios_print(INCHI_IOSTREAM* ios, const char* lpszFormat, ...)
    // INCHI✔❌: {
    // INCHI✔❌:     int ret = 0, ret2 = 0;
    // INCHI✔❌:     va_list argList;
    // INCHI✔❌:
    // INCHI✔❌:     if (!ios)
    // INCHI✔❌:     {
    // INCHI✔❌:         return -1;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (ios->type == INCHI_IOS_TYPE_STRING)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* output to string buffer */
    // INCHI✔❌:         int max_len;
    // INCHI✔❌:         my_va_start(argList, lpszFormat);
    // INCHI✔❌:         max_len = GetMaxPrintfLength(lpszFormat, argList);
    // INCHI✔❌:         va_end(argList);
    // INCHI✔❌:         if (max_len >= 0)
    // INCHI✔❌:         {
    // INCHI✔❌:             /* djb-rwth: fixing oss-fuzz issue #30152 */
    // INCHI✔❌:             int  nAddLength = inchi_max(INCHI_ADD_STR_LEN, max_len);
    // INCHI✔❌:             long long new_str_len = (long long)ios->s.nAllocatedLength + (long long)nAddLength;
    // INCHI✔❌:             if (ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len)
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* enlarge output string */
    // INCHI✔❌:                 char* new_str = (char*)inchi_calloc(new_str_len, sizeof(char)); /* djb-rwth: cast operators added */
    // INCHI✔❌:                 if (new_str)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     if (ios->s.pStr)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         if (ios->s.nUsedLength > 0)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             memcpy(new_str, ios->s.pStr, sizeof(new_str[0]) * ios->s.nUsedLength);
    // INCHI✔❌:                         }
    // INCHI✔❌:                         inchi_free(ios->s.pStr);
    // INCHI✔❌:                     }
    // INCHI✔❌:                     ios->s.pStr = new_str;
    // INCHI✔❌:                     ios->s.nAllocatedLength += nAddLength;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                     return -1; /* failed */
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:             /* output */
    // INCHI✔❌:             my_va_start(argList, lpszFormat);
    // INCHI✔❌:             ret = vsprintf(ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList); /* djb-rwth: not using vsnprintf due to variable length argument */
    // INCHI✔❌:             va_end(argList);
    // INCHI✔❌:             if (ret >= 0)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ios->s.nUsedLength += ret;
    // INCHI✔❌:             }
    // INCHI✔❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌: #if 0
    // INCHI✔❌:             if (FWPRINT)
    // INCHI✔❌:             {
    // INCHI✔❌:                 my_va_start(argList, lpszFormat);
    // INCHI✔❌:                 FWPRINT(lpszFormat, argList);
    // INCHI✔❌:                 va_end(argList);
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌: #endif
    // INCHI✔❌:             return ret;
    // INCHI✔❌:         }
    // INCHI✔❌:         return -1;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     else if (ios->type == INCHI_IOS_TYPE_FILE)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* output to file */
    // INCHI✔❌:         if (ios->f)
    // INCHI✔❌:         {
    // INCHI✔❌:             my_va_start(argList, lpszFormat);
    // INCHI✔❌:             ret = vfprintf(ios->f, lpszFormat, argList);
    // INCHI✔❌:             va_end(argList);
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:         {
    // INCHI✔❌:             my_va_start(argList, lpszFormat);
    // INCHI✔❌:             ret2 = vfprintf(stdout, lpszFormat, argList);
    // INCHI✔❌:             va_end(argList);
    // INCHI✔❌:         }
    // INCHI✔❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌:         if (FWPRINT)
    // INCHI✔❌:         {
    // INCHI✔❌:             my_va_start(argList, lpszFormat);
    // INCHI✔❌:             FWPRINT(lpszFormat, argList);
    // INCHI✔❌:             va_end(argList);
    // INCHI✔❌:         }
    // INCHI✔❌: #endif
    // INCHI✔❌:         return ret ? ret : ret2;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* no output */
    // INCHI✔❌:     return 0;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: inchi_ios_print

    inchi_ios_print_arguments(
        heap,
        ios,
        stdout,
        format,
        &SourcePrintArguments {
            arguments: &arguments.arguments,
            position: arguments.position,
        },
    )
}

pub(crate) fn inchi_ios_print_arguments(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    format: SourceConstPointer<i8>,
    arguments: &SourcePrintArguments<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:477 inchi_ios_print
    // INCHI✔❌: int inchi_ios_print(INCHI_IOSTREAM* ios, const char* lpszFormat, ...)
    // INCHI✔❌: {
    // INCHI✔❌:     int ret = 0, ret2 = 0;
    // INCHI✔❌:     va_list argList;
    // INCHI✔❌:
    // INCHI✔❌:     if (!ios)
    // INCHI✔❌:     {
    // INCHI✔❌:         return -1;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (ios->type == INCHI_IOS_TYPE_STRING)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* output to string buffer */
    // INCHI✔❌:         int max_len;
    // INCHI✔❌:         my_va_start(argList, lpszFormat);
    // INCHI✔❌:         max_len = GetMaxPrintfLength(lpszFormat, argList);
    // INCHI✔❌:         va_end(argList);
    // INCHI✔❌:         if (max_len >= 0)
    // INCHI✔❌:         {
    // INCHI✔❌:             /* djb-rwth: fixing oss-fuzz issue #30152 */
    // INCHI✔❌:             int  nAddLength = inchi_max(INCHI_ADD_STR_LEN, max_len);
    // INCHI✔❌:             long long new_str_len = (long long)ios->s.nAllocatedLength + (long long)nAddLength;
    // INCHI✔❌:             if (ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len)
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* enlarge output string */
    // INCHI✔❌:                 char* new_str = (char*)inchi_calloc(new_str_len, sizeof(char)); /* djb-rwth: cast operators added */
    // INCHI✔❌:                 if (new_str)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     if (ios->s.pStr)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         if (ios->s.nUsedLength > 0)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             memcpy(new_str, ios->s.pStr, sizeof(new_str[0]) * ios->s.nUsedLength);
    // INCHI✔❌:                         }
    // INCHI✔❌:                         inchi_free(ios->s.pStr);
    // INCHI✔❌:                     }
    // INCHI✔❌:                     ios->s.pStr = new_str;
    // INCHI✔❌:                     ios->s.nAllocatedLength += nAddLength;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                     return -1; /* failed */
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:             /* output */
    // INCHI✔❌:             my_va_start(argList, lpszFormat);
    // INCHI✔❌:             ret = vsprintf(ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList); /* djb-rwth: not using vsnprintf due to variable length argument */
    // INCHI✔❌:             va_end(argList);
    // INCHI✔❌:             if (ret >= 0)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ios->s.nUsedLength += ret;
    // INCHI✔❌:             }
    // INCHI✔❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌: #if 0
    // INCHI✔❌:             if (FWPRINT)
    // INCHI✔❌:             {
    // INCHI✔❌:                 my_va_start(argList, lpszFormat);
    // INCHI✔❌:                 FWPRINT(lpszFormat, argList);
    // INCHI✔❌:                 va_end(argList);
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌: #endif
    // INCHI✔❌:             return ret;
    // INCHI✔❌:         }
    // INCHI✔❌:         return -1;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     else if (ios->type == INCHI_IOS_TYPE_FILE)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* output to file */
    // INCHI✔❌:         if (ios->f)
    // INCHI✔❌:         {
    // INCHI✔❌:             my_va_start(argList, lpszFormat);
    // INCHI✔❌:             ret = vfprintf(ios->f, lpszFormat, argList);
    // INCHI✔❌:             va_end(argList);
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:         {
    // INCHI✔❌:             my_va_start(argList, lpszFormat);
    // INCHI✔❌:             ret2 = vfprintf(stdout, lpszFormat, argList);
    // INCHI✔❌:             va_end(argList);
    // INCHI✔❌:         }
    // INCHI✔❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌:         if (FWPRINT)
    // INCHI✔❌:         {
    // INCHI✔❌:             my_va_start(argList, lpszFormat);
    // INCHI✔❌:             FWPRINT(lpszFormat, argList);
    // INCHI✔❌:             va_end(argList);
    // INCHI✔❌:         }
    // INCHI✔❌: #endif
    // INCHI✔❌:         return ret ? ret : ret2;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* no output */
    // INCHI✔❌:     return 0;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: inchi_ios_print

    inchi_ios_print_format_arguments(heap, ios, stdout, SourceFormat::Pointer(format), arguments)
}

pub(crate) fn inchi_ios_print_nodisplay(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    format: SourceConstPointer<i8>,
    arguments: &SourceVaList,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:604 inchi_ios_print_nodisplay
    // INCHI✔️❌: int inchi_ios_print_nodisplay( INCHI_IOSTREAM * ios,
    // INCHI✔️❌:                                const char* lpszFormat,
    // INCHI✔️❌:                                ... )
    // INCHI✔️❌: {
    // INCHI✔️❌:     va_list argList;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!ios)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ios->type == INCHI_IOS_TYPE_STRING)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* output to string buffer */
    // INCHI✔️❌:         int ret = 0, max_len;
    // INCHI✔️❌:         my_va_start( argList, lpszFormat );
    // INCHI✔️❌:         max_len = GetMaxPrintfLength( lpszFormat, argList );
    // INCHI✔️❌:         va_end( argList );
    // INCHI✔️❌:         if (max_len >= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* djb-rwth: fixing oss-fuzz issue #30163 */
    // INCHI✔️❌:             int  nAddLength = inchi_max(INCHI_ADD_STR_LEN, max_len);
    // INCHI✔️❌:             long long new_str_len = (long long)ios->s.nAllocatedLength + (long long)nAddLength;
    // INCHI✔️❌:             if (ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* enlarge output string */
    // INCHI✔️❌:                 char* new_str = (char*)inchi_calloc(new_str_len, sizeof(new_str[0])); /* djb-rwth: cast operators added */
    // INCHI✔️❌:                 if (new_str)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (ios->s.pStr)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (ios->s.nUsedLength > 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             memcpy( new_str, ios->s.pStr, sizeof( new_str[0] )*ios->s.nUsedLength );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         inchi_free( ios->s.pStr );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     ios->s.pStr = new_str;
    // INCHI✔️❌:                     ios->s.nAllocatedLength += nAddLength;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return -1; /* failed */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* output */
    // INCHI✔️❌:             /* djb-rwth: fixing oss-fuzz issue #67676 */
    // INCHI✔️❌:             my_va_start(argList, lpszFormat);
    // INCHI✔️❌:             ret = vsprintf(ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList); /* djb-rwth: not using vsnprintf due to variable length argument; fixing GHI #71 */
    // INCHI✔️❌:             va_end(argList);
    // INCHI✔️❌:             if (ret >= 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ios->s.nUsedLength += ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     else if (ios->type == INCHI_IOS_TYPE_FILE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         my_va_start( argList, lpszFormat );
    // INCHI✔️❌:         inchi_print_nodisplay( ios->f, lpszFormat, argList );
    // INCHI✔️❌:         va_end( argList );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* no output */
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_ios_print_nodisplay

    inchi_ios_print_nodisplay_arguments(
        heap,
        ios,
        stdout,
        format,
        &SourcePrintArguments {
            arguments: &arguments.arguments,
            position: arguments.position,
        },
    )
}

pub(crate) fn inchi_ios_print_nodisplay_arguments(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    format: SourceConstPointer<i8>,
    arguments: &SourcePrintArguments<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:604 inchi_ios_print_nodisplay
    // INCHI✔️❌: int inchi_ios_print_nodisplay( INCHI_IOSTREAM * ios,
    // INCHI✔️❌:                                const char* lpszFormat,
    // INCHI✔️❌:                                ... )
    // INCHI✔️❌: {
    // INCHI✔️❌:     va_list argList;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!ios)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ios->type == INCHI_IOS_TYPE_STRING)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* output to string buffer */
    // INCHI✔️❌:         int ret = 0, max_len;
    // INCHI✔️❌:         my_va_start( argList, lpszFormat );
    // INCHI✔️❌:         max_len = GetMaxPrintfLength( lpszFormat, argList );
    // INCHI✔️❌:         va_end( argList );
    // INCHI✔️❌:         if (max_len >= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* djb-rwth: fixing oss-fuzz issue #30163 */
    // INCHI✔️❌:             int  nAddLength = inchi_max(INCHI_ADD_STR_LEN, max_len);
    // INCHI✔️❌:             long long new_str_len = (long long)ios->s.nAllocatedLength + (long long)nAddLength;
    // INCHI✔️❌:             if (ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* enlarge output string */
    // INCHI✔️❌:                 char* new_str = (char*)inchi_calloc(new_str_len, sizeof(new_str[0])); /* djb-rwth: cast operators added */
    // INCHI✔️❌:                 if (new_str)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (ios->s.pStr)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (ios->s.nUsedLength > 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             memcpy( new_str, ios->s.pStr, sizeof( new_str[0] )*ios->s.nUsedLength );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         inchi_free( ios->s.pStr );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     ios->s.pStr = new_str;
    // INCHI✔️❌:                     ios->s.nAllocatedLength += nAddLength;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return -1; /* failed */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* output */
    // INCHI✔️❌:             /* djb-rwth: fixing oss-fuzz issue #67676 */
    // INCHI✔️❌:             my_va_start(argList, lpszFormat);
    // INCHI✔️❌:             ret = vsprintf(ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList); /* djb-rwth: not using vsnprintf due to variable length argument; fixing GHI #71 */
    // INCHI✔️❌:             va_end(argList);
    // INCHI✔️❌:             if (ret >= 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ios->s.nUsedLength += ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     else if (ios->type == INCHI_IOS_TYPE_FILE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         my_va_start( argList, lpszFormat );
    // INCHI✔️❌:         inchi_print_nodisplay( ios->f, lpszFormat, argList );
    // INCHI✔️❌:         va_end( argList );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* no output */
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_ios_print_nodisplay

    inchi_ios_print_nodisplay_format_arguments(
        heap,
        ios,
        stdout,
        SourceFormat::Pointer(format),
        arguments,
    )
}

pub(crate) fn inchi_ios_eprint(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    format: SourceConstPointer<i8>,
    arguments: &SourceVaList,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:708 inchi_ios_eprint
    // INCHI✔️❌: int inchi_ios_eprint(INCHI_IOSTREAM* ios, const char* lpszFormat, ...)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0, ret2 = 0;
    // INCHI✔️❌:     va_list argList;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!ios)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ios->type == INCHI_IOS_TYPE_STRING)
    // INCHI✔️❌:         /* was #if ( defined(TARGET_API_LIB) || defined(INCHI_STANDALONE_EXE) ) */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* output to string buffer */
    // INCHI✔️❌:         int max_len, nAddLength = 0;
    // INCHI✔️❌:         char* new_str = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:         my_va_start(argList, lpszFormat);
    // INCHI✔️❌:         max_len = GetMaxPrintfLength(lpszFormat, argList);
    // INCHI✔️❌:         va_end(argList);
    // INCHI✔️❌:
    // INCHI✔️❌:         if (max_len >= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* enlarge output string */
    // INCHI✔️❌:                 nAddLength = inchi_max(INCHI_ADD_STR_LEN, max_len);
    // INCHI✔️❌:                 new_str = (char*)inchi_calloc((long long)ios->s.nAllocatedLength + (long long)nAddLength, sizeof(new_str[0])); /* djb-rwth: cast operators added */
    // INCHI✔️❌:                 if (new_str)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (ios->s.pStr)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (ios->s.nUsedLength > 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             memcpy(new_str, ios->s.pStr, sizeof(new_str[0]) * ios->s.nUsedLength);
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         inchi_free(ios->s.pStr);
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     ios->s.pStr = new_str;
    // INCHI✔️❌:                     ios->s.nAllocatedLength += nAddLength;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return -1; /* failed */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* output */
    // INCHI✔️❌:             my_va_start(argList, lpszFormat);
    // INCHI✔️❌:             ret = vsprintf(ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList);
    // INCHI✔️❌:             va_end(argList);
    // INCHI✔️❌:             if (ret >= 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ios->s.nUsedLength += ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     else if (ios->type == INCHI_IOS_TYPE_FILE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (ios->f)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* output to plain file */
    // INCHI✔️❌:             my_va_start(argList, lpszFormat);
    // INCHI✔️❌:             ret = inchi_vfprintf(ios->f, lpszFormat, argList);
    // INCHI✔️❌:             va_end(argList);
    // INCHI✔️❌:             /*  No output to stderr from within dll or GUI program */
    // INCHI✔️❌:             return ret ? ret : ret2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* no output */
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_ios_eprint

    inchi_ios_eprint_arguments(
        heap,
        ios,
        format,
        &SourcePrintArguments {
            arguments: &arguments.arguments,
            position: arguments.position,
        },
    )
}

pub(crate) fn inchi_ios_eprint_arguments(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
    format: SourceConstPointer<i8>,
    arguments: &SourcePrintArguments<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:708 inchi_ios_eprint
    // INCHI✔️❌: int inchi_ios_eprint(INCHI_IOSTREAM* ios, const char* lpszFormat, ...)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0, ret2 = 0;
    // INCHI✔️❌:     va_list argList;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!ios)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ios->type == INCHI_IOS_TYPE_STRING)
    // INCHI✔️❌:         /* was #if ( defined(TARGET_API_LIB) || defined(INCHI_STANDALONE_EXE) ) */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* output to string buffer */
    // INCHI✔️❌:         int max_len, nAddLength = 0;
    // INCHI✔️❌:         char* new_str = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:         my_va_start(argList, lpszFormat);
    // INCHI✔️❌:         max_len = GetMaxPrintfLength(lpszFormat, argList);
    // INCHI✔️❌:         va_end(argList);
    // INCHI✔️❌:
    // INCHI✔️❌:         if (max_len >= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* enlarge output string */
    // INCHI✔️❌:                 nAddLength = inchi_max(INCHI_ADD_STR_LEN, max_len);
    // INCHI✔️❌:                 new_str = (char*)inchi_calloc((long long)ios->s.nAllocatedLength + (long long)nAddLength, sizeof(new_str[0])); /* djb-rwth: cast operators added */
    // INCHI✔️❌:                 if (new_str)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (ios->s.pStr)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (ios->s.nUsedLength > 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             memcpy(new_str, ios->s.pStr, sizeof(new_str[0]) * ios->s.nUsedLength);
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         inchi_free(ios->s.pStr);
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     ios->s.pStr = new_str;
    // INCHI✔️❌:                     ios->s.nAllocatedLength += nAddLength;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return -1; /* failed */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* output */
    // INCHI✔️❌:             my_va_start(argList, lpszFormat);
    // INCHI✔️❌:             ret = vsprintf(ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList);
    // INCHI✔️❌:             va_end(argList);
    // INCHI✔️❌:             if (ret >= 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ios->s.nUsedLength += ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     else if (ios->type == INCHI_IOS_TYPE_FILE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (ios->f)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* output to plain file */
    // INCHI✔️❌:             my_va_start(argList, lpszFormat);
    // INCHI✔️❌:             ret = inchi_vfprintf(ios->f, lpszFormat, argList);
    // INCHI✔️❌:             va_end(argList);
    // INCHI✔️❌:             /*  No output to stderr from within dll or GUI program */
    // INCHI✔️❌:             return ret ? ret : ret2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* no output */
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_ios_eprint

    inchi_ios_eprint_format_arguments(heap, ios, SourceFormat::Pointer(format), arguments)
}

// Concrete output storage for the two existing gets kernels. String reads
// keep their per-byte arena lookup, including alias-visible prior writes.
enum SourceGetsDestination<'a> {
    Pointer(SourceMutPointer<i8>),
    Bytes(&'a mut [i8]),
}

impl SourceGetsDestination<'_> {
    fn bytes_mut<'a>(
        &'a mut self,
        heap: &'a mut SourceHeap,
    ) -> Result<&'a mut [i8], SourceHeapError> {
        match self {
            Self::Pointer(pointer) => heap.slice_mut(*pointer),
            Self::Bytes(bytes) => Ok(bytes),
        }
    }

    fn put(
        &mut self,
        heap: &mut SourceHeap,
        index: usize,
        value: i8,
    ) -> Result<(), SourceHeapError> {
        let bytes = match self {
            Self::Pointer(pointer) => heap.slice_mut(*pointer)?,
            Self::Bytes(bytes) => bytes,
        };
        *bytes
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = value;
        Ok(())
    }

    fn with_bytes<R>(
        &mut self,
        heap: &mut SourceHeap,
        operation: impl FnOnce(&mut [i8], &SourceHeap) -> Result<R, SourceHeapError>,
    ) -> Result<R, SourceHeapError> {
        match self {
            Self::Pointer(pointer) => {
                // slice_mut admits type before signed offset; extraction's
                // internal order differs for malformed typed pointers.
                heap.slice_mut(*pointer)?;
                heap.with_slice_mut_and_heap(*pointer, operation)
            }
            Self::Bytes(bytes) => operation(bytes, heap),
        }
    }
}

pub(crate) fn inchi_sgets(
    heap: &mut SourceHeap,
    output: SourceMutPointer<i8>,
    capacity: i32,
    ios: Option<&mut INCHI_IOSTREAM>,
) -> Result<SourceMutPointer<i8>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1322 inchi_sgets
    // INCHI✔️❌: char* inchi_sgets(char* s, int n, INCHI_IOSTREAM* ios)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int c = 0;
    // INCHI✔️❌:     char* p;
    // INCHI✔️❌:     char* inp;
    // INCHI✔️❌:
    // INCHI✔️❌:     inp = ios->s.pStr + ios->s.nPtr;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (n <= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (NULL == inp)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* like EOF */
    // INCHI✔️❌:         return NULL; /* djb-rwth: addressing coverity ID #499480 -- inp can be NULL */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     p = s;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( *inp == '\0' )
    // INCHI✔️❌:         s = '\0';
    // INCHI✔️❌:     else
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     while (--n > 0 && (c = *inp++))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ios->s.nPtr++;
    // INCHI✔️❌:         if ((*p++ = c) == '\n')
    // INCHI✔️❌:         {
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     *p = '\0';
    // INCHI✔️❌:
    // INCHI✔️❌:     /* printf("\n*** {%-s}",s); */
    // INCHI✔️❌:
    // INCHI✔️❌:     return (c == '\0' && p == s)
    // INCHI✔️❌:         ? NULL /* like EOF reached */
    // INCHI✔️❌:         : s;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_sgets

    let mut destination = SourceGetsDestination::Pointer(output);
    if source_sgets_destination(heap, &mut destination, capacity, ios)? {
        Ok(output)
    } else {
        Ok(SourceMutPointer::null())
    }
}

fn source_sgets_destination(
    heap: &mut SourceHeap,
    destination: &mut SourceGetsDestination<'_>,
    mut capacity: i32,
    ios: Option<&mut INCHI_IOSTREAM>,
) -> Result<bool, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1322 inchi_sgets
    // INCHI✔️❌: char* inchi_sgets(char* s, int n, INCHI_IOSTREAM* ios)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int c = 0;
    // INCHI✔️❌:     char* p;
    // INCHI✔️❌:     char* inp;
    // INCHI✔️❌:
    // INCHI✔️❌:     inp = ios->s.pStr + ios->s.nPtr;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (n <= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (NULL == inp)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* like EOF */
    // INCHI✔️❌:         return NULL; /* djb-rwth: addressing coverity ID #499480 -- inp can be NULL */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     p = s;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( *inp == '\0' )
    // INCHI✔️❌:         s = '\0';
    // INCHI✔️❌:     else
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     while (--n > 0 && (c = *inp++))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ios->s.nPtr++;
    // INCHI✔️❌:         if ((*p++ = c) == '\n')
    // INCHI✔️❌:         {
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     *p = '\0';
    // INCHI✔️❌:
    // INCHI✔️❌:     /* printf("\n*** {%-s}",s); */
    // INCHI✔️❌:
    // INCHI✔️❌:     return (c == '\0' && p == s)
    // INCHI✔️❌:         ? NULL /* like EOF reached */
    // INCHI✔️❌:         : s;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_sgets

    // INCHI✔️❗: One canonical bytewise kernel for arena and fixed stack output.
    // Preserve read, nPtr increment, output store and terminator order; no
    // snapshot hides writes through aliased string input/output pointers.
    // O(reached bytes), O(1) state, no allocation/clone. Safe storage dispatch
    // and whole-call cost remain OPEN INCHI-PERF-GETS-STACK-DISCARD.
    let ios = ios.ok_or(SourceHeapError::NullPointer)?;
    if capacity <= 0 {
        return Ok(false);
    }
    if ios.s.pStr.is_null() {
        return Ok(false);
    }
    let input_offset =
        usize::try_from(ios.s.nPtr).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mut input_index = input_offset;
    let mut output_index = 0_usize;
    let mut character = 0_i8;
    loop {
        capacity -= 1;
        if capacity <= 0 {
            break;
        }
        character = *heap
            .slice(ios.s.pStr.as_const())?
            .get(input_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        input_index += 1;
        if character == 0 {
            break;
        }
        ios.s.nPtr = ios
            .s
            .nPtr
            .checked_add(1)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        destination.put(heap, output_index, character)?;
        output_index += 1;
        if character == b'\n' as i8 {
            break;
        }
    }
    destination.put(heap, output_index, 0)?;
    if character == 0 && output_index == 0 {
        Ok(false)
    } else {
        Ok(true)
    }
}

fn source_file_fgets(
    heap: &mut SourceHeap,
    output: SourceMutPointer<i8>,
    capacity: i32,
    file: SourceMutPointer<FILE>,
) -> Result<SourceMutPointer<i8>, SourceHeapError> {
    // INCHI source ichi_io.c:1006,1012:
    //         if (NULL != (p = fgets(line, line_len, finp)) &&
    //             while (NULL != fgets(temp, sizeof(temp), finp) && NULL == strchr(temp, '\n'))
    // BEGIN INCHI SUPPLEMENTAL LIBC: third_party/musl/src/stdio/fgets.c (configured native identity unresolved)
    /*
    char *fgets(char *restrict s, int n, FILE *restrict f)
    {
        char *p = s;
        unsigned char *z;
        size_t k;
        int c;

        FLOCK(f);

        if (n<=1) {
            f->mode |= f->mode-1;
            FUNLOCK(f);
            if (n<1) return 0;
            *s = 0;
            return s;
        }
        n--;

        while (n) {
            if (f->rpos != f->rend) {
                z = memchr(f->rpos, '\n', f->rend - f->rpos);
                k = z ? z - f->rpos + 1 : f->rend - f->rpos;
                k = MIN(k, n);
                memcpy(p, f->rpos, k);
                f->rpos += k;
                p += k;
                n -= k;
                if (z || !n) break;
            }
            if ((c = getc_unlocked(f)) < 0) {
                if (p==s || !feof(f)) s = 0;
                break;
            }
            n--;
            if ((*p++ = c) == '\n') break;
        }
        if (s) *p = 0;

        FUNLOCK(f);

        return s;
    }
        */
    // END INCHI SUPPLEMENTAL LIBC: fgets
    // INCHI source third_party/musl/src/stdio/__stdio_read.c:
    // 	if (cnt <= 0) {
    // 		f->flags |= cnt ? F_ERR : F_EOF;
    // 		return 0;
    // 	}
    // INCHI source third_party/musl/src/stdio/__toread.c:
    // 	return (f->flags & F_EOF) ? EOF : 0;
    // INCHI✔️❗: Reproduce bounded prefix/terminator and EOF probe semantics
    // for the existing owned-byte SourceFile model, including no-read capacity1.
    // SourceFile.error retains its existing modeled forced-read-failure meaning.
    // Supplemental libc is source evidence, not configured-native certification.
    // Cost: O(copied prefix), O(1) scratch, no whole-file clone/allocation.
    // Safe disjoint arena extraction/checks versus native buffering constants
    // and WHOLE-call cost remain independently OPEN C057.
    let mut destination = SourceGetsDestination::Pointer(output);
    if source_file_fgets_destination(heap, &mut destination, capacity, file)? {
        Ok(output)
    } else {
        Ok(SourceMutPointer::null())
    }
}

fn source_file_fgets_destination(
    heap: &mut SourceHeap,
    destination: &mut SourceGetsDestination<'_>,
    capacity: i32,
    file: SourceMutPointer<FILE>,
) -> Result<bool, SourceHeapError> {
    // INCHI source ichi_io.c:1006,1012:
    //         if (NULL != (p = fgets(line, line_len, finp)) &&
    //             while (NULL != fgets(temp, sizeof(temp), finp) && NULL == strchr(temp, '\n'))
    // BEGIN INCHI SUPPLEMENTAL LIBC: third_party/musl/src/stdio/fgets.c (configured native identity unresolved)
    /*
    char *fgets(char *restrict s, int n, FILE *restrict f)
    {
        char *p = s;
        unsigned char *z;
        size_t k;
        int c;

        FLOCK(f);

        if (n<=1) {
            f->mode |= f->mode-1;
            FUNLOCK(f);
            if (n<1) return 0;
            *s = 0;
            return s;
        }
        n--;

        while (n) {
            if (f->rpos != f->rend) {
                z = memchr(f->rpos, '\n', f->rend - f->rpos);
                k = z ? z - f->rpos + 1 : f->rend - f->rpos;
                k = MIN(k, n);
                memcpy(p, f->rpos, k);
                f->rpos += k;
                p += k;
                n -= k;
                if (z || !n) break;
            }
            if ((c = getc_unlocked(f)) < 0) {
                if (p==s || !feof(f)) s = 0;
                break;
            }
            n--;
            if ((*p++ = c) == '\n') break;
        }
        if (s) *p = 0;

        FUNLOCK(f);

        return s;
    }
        */
    // END INCHI SUPPLEMENTAL LIBC: fgets
    // INCHI source third_party/musl/src/stdio/__stdio_read.c:
    // 	if (cnt <= 0) {
    // 		f->flags |= cnt ? F_ERR : F_EOF;
    // 		return 0;
    // 	}
    // INCHI source third_party/musl/src/stdio/__toread.c:
    // 	return (f->flags & F_EOF) ? EOF : 0;
    // INCHI✔️❗: Reproduce bounded prefix/terminator and EOF probe semantics
    // for the existing owned-byte SourceFile model, including no-read capacity1.
    // SourceFile.error retains its existing modeled forced-read-failure meaning.
    // Supplemental libc is source evidence, not configured-native certification.
    // Cost: O(copied prefix), O(1) scratch, no whole-file clone/allocation.
    // Safe disjoint arena extraction/checks versus native buffering constants
    // and WHOLE-call cost remain independently OPEN C057.
    // INCHI✔️❗: Share FILE state/EOF/copy semantics for arena and fixed stack
    // output. Destination admission precedes the second FILE access and all
    // position/EOF publication. O(copied prefix), O(1) state, no allocation.
    // Whole-call constants remain OPEN INCHI-PERF-GETS-STACK-DISCARD.
    if capacity <= 0 {
        return Ok(false);
    }
    let (start, end, reached_eof) = {
        let state = heap
            .slice(file.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if capacity == 1 {
            destination.put(heap, 0, 0)?;
            return Ok(true);
        }
        let start =
            usize::try_from(state.position).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if state.eof {
            return Ok(false);
        }
        if state.error {
            // musl __stdio_read.c: f->flags |= cnt ? F_ERR : F_EOF;
            // Capacity1 above performs no backend read or sticky update.
            heap.slice_mut(file)?[0].error_indicator = true;
            return Ok(false);
        }
        if start >= state.bytes.len() {
            heap.slice_mut(file)?[0].eof = true;
            return Ok(false);
        }
        let limit =
            usize::try_from(capacity - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let mut end = start;
        let mut newline = false;
        while end - start < limit && end < state.bytes.len() {
            let byte = state.bytes[end];
            end += 1;
            if byte == b'\n' {
                newline = true;
                break;
            }
        }
        let reached_eof = !newline && end - start < limit && end == state.bytes.len();
        (start, end, reached_eof)
    };
    let count = end - start;
    // Preserve the original destination type/offset/extent admission order.
    // A valid byte destination cannot alias the already-admitted FILE record.
    destination.with_bytes(heap, |destination, heap| {
        if destination.len() < count + 1 {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let state = heap
            .slice(file.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        destination[..count]
            .copy_from_slice(bytemuck::cast_slice::<u8, i8>(&state.bytes[start..end]));
        destination[count] = 0;
        Ok(())
    })?;
    let state = &mut heap.slice_mut(file)?[0];
    state.position = u64::try_from(end).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if reached_eof {
        state.eof = true;
    }
    Ok(true)
}

fn source_line_ending_or_nul(bytes: &[i8], ending: i8) -> Result<usize, SourceHeapError> {
    // INCHI source ichi_io.c:1039:
    //         if ((q = strchr(line, '\r'))) /* djb-rwth: addressing LLVM warning */
    // INCHI source ichi_io.c:405:
    //         p = strchr(szLine, '\n');
    // INCHI✔️❗: Inspect the determining prefix only. Reached missing bytes
    // remain MissingNulTerminator; current callers provide terminated lines.
    // O(first ending/NUL position), O(1) space, no allocation or extra scan.
    // Safe iterator/Result versus libc whole-call constants remain OPEN C055.
    bytes
        .iter()
        .position(|byte| *byte == 0 || *byte == ending)
        .ok_or(SourceHeapError::MissingNulTerminator)
}

fn source_line_too_long(bytes: &[i8], content_limit: i32) -> Result<i32, SourceHeapError> {
    // INCHI source ichi_io.c:405..406 (also getsTab/getsTab1):
    //         p = strchr(szLine, '\n');
    //         *bTooLongLine = (!p && ((int)strlen(szLine)) == len - 2);
    // INCHI✔️❗: With LF, do not evaluate strlen or narrow its result. Without
    // LF, the determining NUL position already supplies strlen exactly.
    // All callers wrote NUL at nonnegative i32 len-1 before admission, so
    // the reached length fits i32. Existing typed overflow checks stay intact.
    // O(first LF/NUL), O(1) space; no second scan, allocation or clone.
    // Independent libc/check/whole-call cost remains OPEN C055.
    let position = source_line_ending_or_nul(bytes, b'\n' as i8)?;
    if bytes[position] == b'\n' as i8 {
        Ok(0)
    } else {
        let length = i32::try_from(position).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        Ok(i32::from(length == content_limit))
    }
}

#[allow(non_snake_case)]
pub(crate) fn inchi_fgetsLf(
    heap: &mut SourceHeap,
    line: SourceMutPointer<i8>,
    line_length: i32,
    input: Option<&mut INCHI_IOSTREAM>,
) -> Result<SourceMutPointer<i8>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:996 inchi_fgetsLf
    // INCHI✔️❌: char* inchi_fgetsLf(char* line, int line_len, INCHI_IOSTREAM* inp_stream)
    // INCHI✔️❌: {
    // INCHI✔️❌:     char* p = NULL, * q;
    // INCHI✔️❌:     FILE* finp = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (inp_stream->type == INCHI_IOS_TYPE_FILE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Read from file */
    // INCHI✔️❌:         finp = inp_stream->f;
    // INCHI✔️❌:         memset(line, 0, line_len); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         if (NULL != (p = fgets(line, line_len, finp)) &&
    // INCHI✔️❌:             NULL == strchr(p, '\n'))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             char temp[64];
    // INCHI✔️❌:             /* bypass up to '\n' or up to end of file whichever comes first */
    // INCHI✔️❌:             while (NULL != fgets(temp, sizeof(temp), finp) && NULL == strchr(temp, '\n'))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else if (inp_stream->type == INCHI_IOS_TYPE_STRING)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Read from supplied string representing Molfile */
    // INCHI✔️❌:         memset(line, 0, line_len); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         if (NULL != (p = inchi_sgets(line, line_len, inp_stream)) &&
    // INCHI✔️❌:             NULL == strchr(p, '\n'))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             char temp[64];
    // INCHI✔️❌:             /* bypass up to '\n' or up to end of file whichever comes first */
    // INCHI✔️❌:             while (NULL != inchi_sgets(temp, sizeof(temp), inp_stream) && NULL == strchr(temp, '\n'))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (p)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((q = strchr(line, '\r'))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  fix CR CR LF line terminator. */
    // INCHI✔️❌:             q[0] = '\n';
    // INCHI✔️❌:             q[1] = '\0';
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return p;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_fgetsLf

    let mut destination = SourceGetsDestination::Pointer(line);
    if source_fgets_lf_destination(heap, &mut destination, line_length, input)? {
        Ok(line)
    } else {
        Ok(SourceMutPointer::null())
    }
}

#[allow(non_snake_case)]
pub(crate) fn inchi_fgetsLf_bytes(
    heap: &mut SourceHeap,
    line: &mut [i8],
    line_length: i32,
    input: Option<&mut INCHI_IOSTREAM>,
) -> Result<bool, SourceHeapError> {
    // Native mol_fmt2.c419..420 passes its local line array to inchi_fgetsLf.
    // Same complete kernel as pointer callers; no arena registration/copy.
    source_fgets_lf_destination(
        heap,
        &mut SourceGetsDestination::Bytes(line),
        line_length,
        input,
    )
}

fn source_fgets_lf_destination(
    heap: &mut SourceHeap,
    line: &mut SourceGetsDestination<'_>,
    line_length: i32,
    input: Option<&mut INCHI_IOSTREAM>,
) -> Result<bool, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:996 inchi_fgetsLf
    // INCHI✔️❌: char* inchi_fgetsLf(char* line, int line_len, INCHI_IOSTREAM* inp_stream)
    // INCHI✔️❌: {
    // INCHI✔️❌:     char* p = NULL, * q;
    // INCHI✔️❌:     FILE* finp = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (inp_stream->type == INCHI_IOS_TYPE_FILE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Read from file */
    // INCHI✔️❌:         finp = inp_stream->f;
    // INCHI✔️❌:         memset(line, 0, line_len); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         if (NULL != (p = fgets(line, line_len, finp)) &&
    // INCHI✔️❌:             NULL == strchr(p, '\n'))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             char temp[64];
    // INCHI✔️❌:             /* bypass up to '\n' or up to end of file whichever comes first */
    // INCHI✔️❌:             while (NULL != fgets(temp, sizeof(temp), finp) && NULL == strchr(temp, '\n'))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else if (inp_stream->type == INCHI_IOS_TYPE_STRING)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Read from supplied string representing Molfile */
    // INCHI✔️❌:         memset(line, 0, line_len); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:         if (NULL != (p = inchi_sgets(line, line_len, inp_stream)) &&
    // INCHI✔️❌:             NULL == strchr(p, '\n'))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             char temp[64];
    // INCHI✔️❌:             /* bypass up to '\n' or up to end of file whichever comes first */
    // INCHI✔️❌:             while (NULL != inchi_sgets(temp, sizeof(temp), inp_stream) && NULL == strchr(temp, '\n'))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (p)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((q = strchr(line, '\r'))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  fix CR CR LF line terminator. */
    // INCHI✔️❌:             q[0] = '\n';
    // INCHI✔️❌:             q[1] = '\0';
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return p;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_fgetsLf

    // INCHI✔️❗: Sole full line/discard/CR kernel for source pointer and local
    // array storage. Pointer type/offset admission precedes extraction; string
    // alias reads/writes stay bytewise. O(reached bytes), fixed64byte scratch,
    // no allocation/clone; whole-call cost OPEN INCHI-PERF-MOL-COPY-STACK.
    let input = input.ok_or(SourceHeapError::NullPointer)?;
    let mut result = false;
    if input.type_ == INCHI_IOS_TYPE_FILE as i32 {
        if line_length > 0 {
            let length =
                usize::try_from(line_length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let bytes = line.bytes_mut(heap)?;
            if bytes.len() < length {
                return Err(SourceHeapError::PointerOutOfBounds);
            }
            bytes[..length].fill(0);
        }
        result = source_file_fgets_destination(heap, line, line_length, input.f)?;
        if result
            && !line
                .bytes_mut(heap)?
                .iter()
                .take_while(|byte| **byte != 0)
                .any(|byte| *byte == b'\n' as i8)
        {
            // Native temp[64] has stack lifetime. Safe initialization adds
            // bounded zeroing; whole-call cost is deferred, not certified.
            let mut temporary = [0_i8; 64];
            loop {
                let mut destination = SourceGetsDestination::Bytes(&mut temporary);
                let chunk = source_file_fgets_destination(heap, &mut destination, 64, input.f)?;
                if !chunk
                    || temporary
                        .iter()
                        .take_while(|byte| **byte != 0)
                        .any(|byte| *byte == b'\n' as i8)
                {
                    break;
                }
            }
        }
    } else if input.type_ == INCHI_IOS_TYPE_STRING as i32 {
        if line_length > 0 {
            let length =
                usize::try_from(line_length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let bytes = line.bytes_mut(heap)?;
            if bytes.len() < length {
                return Err(SourceHeapError::PointerOutOfBounds);
            }
            bytes[..length].fill(0);
        }
        result = source_sgets_destination(heap, line, line_length, Some(input))?;
        if result
            && !line
                .bytes_mut(heap)?
                .iter()
                .take_while(|byte| **byte != 0)
                .any(|byte| *byte == b'\n' as i8)
        {
            // Native temp[64] has stack lifetime. Safe initialization adds
            // bounded zeroing; whole-call cost is deferred, not certified.
            let mut temporary = [0_i8; 64];
            loop {
                let mut destination = SourceGetsDestination::Bytes(&mut temporary);
                let chunk = source_sgets_destination(heap, &mut destination, 64, Some(input))?;
                if !chunk
                    || temporary
                        .iter()
                        .take_while(|byte| **byte != 0)
                        .any(|byte| *byte == b'\n' as i8)
                {
                    break;
                }
            }
        }
    }

    if result {
        let bytes = line.bytes_mut(heap)?;
        // INCHI✔️❗: Stop at the first CR or NUL, as source strchr does.
        // Produced successful lines are NUL-terminated; no whole-line preflight.
        let carriage_return = source_line_ending_or_nul(bytes, b'\r' as i8)?;
        if bytes[carriage_return] == b'\r' as i8 {
            bytes[carriage_return] = b'\n' as i8;
            *bytes
                .get_mut(carriage_return + 1)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
        }
    }
    Ok(result)
}

#[allow(non_snake_case)]
pub(crate) fn GetMaxPrintfLength(
    heap: &SourceHeap,
    format: SourceConstPointer<i8>,
    arguments: &mut SourceVaList,
) -> Result<i32, SourceHeapError> {
    let mut cursor = FormatCursor {
        borrowed_string: None,
        arguments: &arguments.arguments,
        position: &mut arguments.position,
    };
    get_max_printf_length_cursor(heap, format, &mut cursor)
}

fn get_max_printf_length_cursor(
    heap: &SourceHeap,
    format: SourceConstPointer<i8>,
    arguments: &mut FormatCursor<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1065 GetMaxPrintfLength
    // INCHI✔️❌: int GetMaxPrintfLength(const char* lpszFormat, va_list argList)
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*ASSERT(AfxIsValidString(lpszFormat, FALSE));*/
    // INCHI✔️❌:     const char* lpsz;
    // INCHI✔️❌:     int nMaxLen, nWidth, nPrecision, nModifier, nItemLen;
    // INCHI✔️❌:
    // INCHI✔️❌:     nMaxLen = 0;
    // INCHI✔️❌:     /* make a guess at the maximum length of the resulting string */
    // INCHI✔️❌:     for (lpsz = lpszFormat; *lpsz; lpsz++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* moved from below for C syntax reason - 2024-09-01 DT */
    // INCHI✔️❌:         /* djb-rwth: return values needed for va_arg; djb-rwth: ignoring LLVM warnings: function returning value */
    // INCHI✔️❌:         int ivarg;
    // INCHI✔️❌:         double dvarg;
    // INCHI✔️❌:         void* ivvarg;
    // INCHI✔️❌:         int* ipvarg;
    // INCHI✔️❌:
    // INCHI✔️❌:         /* handle '%' character, but watch out for '%%' */
    // INCHI✔️❌:         if (*lpsz != '%' || *(++lpsz) == '%')
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nMaxLen += 1;
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         nItemLen = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         /*  handle '%' character with format */
    // INCHI✔️❌:         nWidth = 0;
    // INCHI✔️❌:         for (; *lpsz; lpsz++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* check for valid flags */
    // INCHI✔️❌:             if (*lpsz == '#')
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nMaxLen += 2;   /* for '0x' */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else if (*lpsz == '*')
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nWidth = va_arg(argList, int);
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else if (*lpsz == '-' || *lpsz == '+' || *lpsz == '0'
    // INCHI✔️❌:                 || *lpsz == ' ')
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else /* hit non-flag character */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* get width and skip it */
    // INCHI✔️❌:         if (nWidth == 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* width indicated by */
    // INCHI✔️❌:             nWidth = atoi(lpsz);
    // INCHI✔️❌:             for (; *lpsz && isdigit(*lpsz); lpsz++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /*ASSERT(nWidth >= 0);*/
    // INCHI✔️❌:         if (nWidth < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_error; /* instead of exception */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         nPrecision = 0;
    // INCHI✔️❌:         if (*lpsz == '.')
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* skip past '.' separator (width.precision)*/
    // INCHI✔️❌:             lpsz++;
    // INCHI✔️❌:
    // INCHI✔️❌:             /* get precision and skip it*/
    // INCHI✔️❌:             if (*lpsz == '*')
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nPrecision = va_arg(argList, int);
    // INCHI✔️❌:                 lpsz++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nPrecision = atoi(lpsz);
    // INCHI✔️❌:                 for (; *lpsz && isdigit(*lpsz); lpsz++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (nPrecision < 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto exit_error; /* instead of exception */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* should be on type modifier or specifier */
    // INCHI✔️❌:         nModifier = 0;
    // INCHI✔️❌:         switch (*lpsz)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* modifiers that affect size */
    // INCHI✔️❌:         case 'h':
    // INCHI✔️❌:             switch (lpsz[1])
    // INCHI✔️❌:             {
    // INCHI✔️❌:             case 'd':
    // INCHI✔️❌:             case 'i':
    // INCHI✔️❌:             case 'o':
    // INCHI✔️❌:             case 'x':
    // INCHI✔️❌:             case 'X':
    // INCHI✔️❌:             case 'u':
    // INCHI✔️❌:                 /* short unsigned, short double, etc. -- added to the original MS example */
    // INCHI✔️❌:                 /* ignore the fact that these modifiers do affect size */
    // INCHI✔️❌:                 lpsz++;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             default:
    // INCHI✔️❌:                 nModifier = FORCE_ANSI;
    // INCHI✔️❌:                 lpsz++;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         case 'l':
    // INCHI✔️❌:             switch (lpsz[1])
    // INCHI✔️❌:             {
    // INCHI✔️❌:             case 'd':
    // INCHI✔️❌:             case 'i':
    // INCHI✔️❌:             case 'o':
    // INCHI✔️❌:             case 'x':
    // INCHI✔️❌:             case 'X':
    // INCHI✔️❌:             case 'u':
    // INCHI✔️❌:             case 'f': /* long float -- post ANSI C */
    // INCHI✔️❌:                 /* long unsigned, long double, etc. -- added to the original MS example */
    // INCHI✔️❌:                 /* ignore the fact that these modifiers do affect size */
    // INCHI✔️❌:                 lpsz++;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             default:
    // INCHI✔️❌:                 /*
    // INCHI✔️❌:                 nModifier = FORCE_UNICODE;
    // INCHI✔️❌:                 lpsz ++;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:                 */
    // INCHI✔️❌:                 goto exit_error;  /* no UNICODE, please */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             break;
    // INCHI✔️❌:             /* modifiers that do not affect size */
    // INCHI✔️❌:         case 'F':
    // INCHI✔️❌:         case 'N':
    // INCHI✔️❌:         case 'L':
    // INCHI✔️❌:             lpsz++;
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* now should be on specifier */
    // INCHI✔️❌:
    // INCHI✔️❌:         switch (*lpsz | nModifier)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* single characters*/
    // INCHI✔️❌:         case 'c':
    // INCHI✔️❌:         case 'C':
    // INCHI✔️❌:             nItemLen = 2;
    // INCHI✔️❌:             ivarg = va_arg(argList, int); /* djb-rwth: int return value; ignoring LLVM warning */
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         case 'c' | FORCE_ANSI:
    // INCHI✔️❌:         case 'C' | FORCE_ANSI:
    // INCHI✔️❌:             nItemLen = 2;
    // INCHI✔️❌:             ivarg = va_arg(argList, int); /* djb-rwth: int return value; ignoring LLVM warning */
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         case 'c' | FORCE_UNICODE:
    // INCHI✔️❌:         case 'C' | FORCE_UNICODE:
    // INCHI✔️❌:             goto exit_error;  /* no UNICODE, please */
    // INCHI✔️❌:             /*
    // INCHI✔️❌:             nItemLen = 2;
    // INCHI✔️❌:             va_arg(argList, wchar_t);
    // INCHI✔️❌:             break;
    // INCHI✔️❌:             */
    // INCHI✔️❌:
    // INCHI✔️❌:             /* strings*/
    // INCHI✔️❌:         case 's':
    // INCHI✔️❌:         case 'S':
    // INCHI✔️❌:             nItemLen = (int)strlen(va_arg(argList, char*));
    // INCHI✔️❌:             nItemLen = inchi_max(1, nItemLen);
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         case 's' | FORCE_ANSI:
    // INCHI✔️❌:         case 'S' | FORCE_ANSI:
    // INCHI✔️❌:             nItemLen = (int)strlen(va_arg(argList, char*));
    // INCHI✔️❌:             nItemLen = inchi_max(1, nItemLen);
    // INCHI✔️❌:             break;
    // INCHI✔️❌:
    // INCHI✔️❌:         case 's' | FORCE_UNICODE:
    // INCHI✔️❌:         case 'S' | FORCE_UNICODE:
    // INCHI✔️❌:             goto exit_error;  /* no UNICODE, please */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* adjust nItemLen for strings */
    // INCHI✔️❌:         if (nItemLen != 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nItemLen = inchi_max(nItemLen, nWidth);
    // INCHI✔️❌:             if (nPrecision != 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nItemLen = inchi_min(nItemLen, nPrecision);
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             switch (*lpsz)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* integers */
    // INCHI✔️❌:             case 'd':
    // INCHI✔️❌:             case 'i':
    // INCHI✔️❌:             case 'u':
    // INCHI✔️❌:             case 'x':
    // INCHI✔️❌:             case 'X':
    // INCHI✔️❌:             case 'o':
    // INCHI✔️❌:                 ivarg = va_arg(argList, int); /* djb-rwth: int return value; ignoring LLVM warning */
    // INCHI✔️❌:                 nItemLen = 32;
    // INCHI✔️❌:                 nItemLen = inchi_max(nItemLen, nWidth + nPrecision);
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:
    // INCHI✔️❌:             case 'e':
    // INCHI✔️❌:             case 'f':
    // INCHI✔️❌:             case 'g':
    // INCHI✔️❌:             case 'G':
    // INCHI✔️❌:                 dvarg = va_arg(argList, double); /* djb-rwth: double return value; ignoring LLVM warning */
    // INCHI✔️❌:                 nItemLen = 32;
    // INCHI✔️❌:                 nItemLen = inchi_max(nItemLen, nWidth + nPrecision);
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:
    // INCHI✔️❌:             case 'p':
    // INCHI✔️❌:                 ivvarg = va_arg(argList, void*); /* djb-rwth: void* return value; ignoring LLVM warning */
    // INCHI✔️❌:                 nItemLen = 32;
    // INCHI✔️❌:                 nItemLen = inchi_max(nItemLen, nWidth + nPrecision);
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* no output */
    // INCHI✔️❌:             case 'n':
    // INCHI✔️❌:                 ipvarg = va_arg(argList, int*); /* djb-rwth: int* return value; ignoring LLVM warning */
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:
    // INCHI✔️❌:             default:
    // INCHI✔️❌:                 /*ASSERT(FALSE);*/  /* unknown formatting option*/
    // INCHI✔️❌:                 goto exit_error; /* instead of exception */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* adjust nMaxLen for output nItemLen */
    // INCHI✔️❌:         nMaxLen += nItemLen;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return nMaxLen;
    // INCHI✔️❌:
    // INCHI✔️❌: exit_error:
    // INCHI✔️❌:     return -1; /* wrong format */
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: GetMaxPrintfLength

    get_max_printf_length_format_cursor(heap, SourceFormat::Pointer(format), arguments, None)
}

pub(crate) fn inchi_ios_init(
    ios: Option<&mut INCHI_IOSTREAM>,
    io_type: i32,
    file: SourceMutPointer<FILE>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:85 inchi_ios_init
    // INCHI✔️✔️: void inchi_ios_init(INCHI_IOSTREAM* ios, int io_type, FILE* f)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     memset(ios, 0, sizeof(*ios)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:     switch (io_type)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:     case INCHI_IOS_TYPE_FILE:
    // INCHI✔️✔️:         ios->type = INCHI_IOS_TYPE_FILE;
    // INCHI✔️✔️:         break;
    // INCHI✔️✔️:     case INCHI_IOS_TYPE_STRING:
    // INCHI✔️✔️:     default:
    // INCHI✔️✔️:         ios->type = INCHI_IOS_TYPE_STRING;
    // INCHI✔️✔️:         break;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     ios->f = f;
    // INCHI✔️✔️:     return;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: inchi_ios_init

    let ios = ios.ok_or(SourceHeapError::NullPointer)?;
    *ios = INCHI_IOSTREAM::default();
    ios.type_ = if io_type == INCHI_IOS_TYPE_FILE as i32 {
        INCHI_IOS_TYPE_FILE as i32
    } else {
        INCHI_IOS_TYPE_STRING as i32
    };
    ios.f = file;
    Ok(())
}

pub(crate) fn inchi_ios_close(
    heap: &mut SourceHeap,
    ios: Option<&mut INCHI_IOSTREAM>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:229 inchi_ios_close
    // INCHI✔️❌: void inchi_ios_close(INCHI_IOSTREAM* ios)
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (NULL == ios)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (ios->s.pStr)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free(ios->s.pStr);
    // INCHI✔️❌:     }
    // INCHI✔️❌:     ios->s.pStr = NULL;
    // INCHI✔️❌:     ios->s.nUsedLength = ios->s.nAllocatedLength = ios->s.nPtr = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (NULL != ios->f && stdout != ios->f && stderr != ios->f && stdin != ios->f)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         fclose(ios->f);
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_ios_close

    let Some(ios) = ios else {
        return Ok(());
    };
    if !ios.s.pStr.is_null() {
        inchi_free(heap, ios.s.pStr)?;
    }
    ios.s.pStr = SourceMutPointer::null();
    ios.s.nUsedLength = 0;
    ios.s.nAllocatedLength = 0;
    ios.s.nPtr = 0;
    if !ios.f.is_null() {
        let is_standard = heap
            .slice(ios.f.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .is_standard_stream;
        if !is_standard {
            heap.free(ios.f)?;
        }
    }
    Ok(())
}

pub(crate) fn inchi_ios_reset(
    heap: &mut SourceHeap,
    ios: &mut INCHI_IOSTREAM,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:255 inchi_ios_reset
    // INCHI✔️❌: void inchi_ios_reset(INCHI_IOSTREAM* ios)
    // INCHI✔️❌: {
    // INCHI✔️❌:     ios->s.pStr = NULL;
    // INCHI✔️❌:     ios->s.nUsedLength = ios->s.nAllocatedLength = ios->s.nPtr = 0;
    // INCHI✔️❌:     if (NULL != ios->f && stdout != ios->f && stderr != ios->f && stdin != ios->f)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         fclose(ios->f);
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_ios_reset

    ios.s.pStr = SourceMutPointer::null();
    ios.s.nUsedLength = 0;
    ios.s.nAllocatedLength = 0;
    ios.s.nPtr = 0;
    if !ios.f.is_null() {
        let is_standard = heap
            .slice(ios.f.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .is_standard_stream;
        if !is_standard {
            heap.free(ios.f)?;
        }
    }
    Ok(())
}

pub(crate) fn inchi_ios_str_getc(
    heap: &mut SourceHeap,
    ios: &mut INCHI_IOSTREAM,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:294 inchi_ios_str_getc
    // INCHI✔️❌: int inchi_ios_str_getc(INCHI_IOSTREAM* ios)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int c;
    // INCHI✔️❌:     if (ios->type == INCHI_IOS_TYPE_STRING)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (ios->s.nPtr < ios->s.nUsedLength)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return (int)ios->s.pStr[ios->s.nPtr++];
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return EOF;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     else if (ios->type == INCHI_IOS_TYPE_FILE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         c = fgetc(ios->f);
    // INCHI✔️❌:         if (ferror(ios->f))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             c = EOF;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return c;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* error */
    // INCHI✔️❌:     return EOF;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_ios_str_getc

    if ios.type_ == INCHI_IOS_TYPE_STRING as i32 {
        if ios.s.nPtr < ios.s.nUsedLength {
            let index =
                usize::try_from(ios.s.nPtr).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let value = *heap
                .slice(ios.s.pStr.as_const())?
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            ios.s.nPtr += 1;
            return Ok(i32::from(value));
        }
        return Ok(SOURCE_EOF);
    }

    if ios.type_ == INCHI_IOS_TYPE_FILE as i32 {
        let file = heap
            .slice_mut(ios.f)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        // Native ichi_io.c308..314 reads first, then inspects sticky ferror.
        // Unlike getInChIChar, it returns EOF even after a successful read if
        // an earlier error indicator remains set; the consumed byte stays.
        let value = file.read_byte()?;
        return Ok(if file.error_indicator {
            SOURCE_EOF
        } else {
            value
        });
    }

    Ok(SOURCE_EOF)
}

pub(crate) fn inchi_ios_str_gets(
    heap: &mut SourceHeap,
    line: SourceMutPointer<i8>,
    len: i32,
    ios: &mut INCHI_IOSTREAM,
) -> Result<SourceMutPointer<i8>, SourceHeapError> {
    let mut target = SourceLineBuffer::Pointer(line);
    let position = target.start();
    let found = source_read_line_raw(heap, &mut target, position, len, ios, false)?;
    Ok(if found {
        line
    } else {
        SourceMutPointer::null()
    })
}

pub(crate) fn inchi_ios_str_getsTab(
    heap: &mut SourceHeap,
    line: SourceMutPointer<i8>,
    len: i32,
    ios: &mut INCHI_IOSTREAM,
) -> Result<SourceMutPointer<i8>, SourceHeapError> {
    let mut target = SourceLineBuffer::Pointer(line);
    let position = target.start();
    let found = source_read_line_raw(heap, &mut target, position, len, ios, true)?;
    Ok(if found {
        line
    } else {
        SourceMutPointer::null()
    })
}

pub(crate) fn inchi_ios_gets(
    heap: &mut SourceHeap,
    line: SourceMutPointer<i8>,
    len: i32,
    ios: &mut INCHI_IOSTREAM,
    too_long_line: &mut i32,
) -> Result<i32, SourceHeapError> {
    let mut target = SourceLineBuffer::Pointer(line);
    let position = target.start();
    source_read_line(
        heap,
        &mut target,
        position,
        len,
        ios,
        too_long_line,
        false,
        true,
    )
}

pub(crate) fn inchi_ios_getsTab(
    heap: &mut SourceHeap,
    line: SourceMutPointer<i8>,
    len: i32,
    ios: &mut INCHI_IOSTREAM,
    too_long_line: &mut i32,
) -> Result<i32, SourceHeapError> {
    let mut target = SourceLineBuffer::Pointer(line);
    let position = target.start();
    source_read_line(
        heap,
        &mut target,
        position,
        len,
        ios,
        too_long_line,
        true,
        true,
    )
}

pub(crate) fn inchi_ios_getsTab1(
    heap: &mut SourceHeap,
    line: SourceMutPointer<i8>,
    len: i32,
    ios: &mut INCHI_IOSTREAM,
    too_long_line: &mut i32,
) -> Result<i32, SourceHeapError> {
    let mut target = SourceLineBuffer::Pointer(line);
    let position = target.start();
    source_read_line(
        heap,
        &mut target,
        position,
        len,
        ios,
        too_long_line,
        true,
        false,
    )
}

// Temporary source line projection: real arena pointers or real fixed arrays.
// No source identity is fabricated for a borrowed array.
pub(crate) enum SourceLineBuffer<'a> {
    Pointer(SourceMutPointer<i8>),
    Borrowed(&'a mut [i8]),
}

#[derive(Clone, Copy)]
pub(crate) enum SourceLineCursor {
    Stored(SourceMutPointer<i8>),
    Relative(i64),
}

impl SourceLineCursor {
    pub(crate) fn is_null(self) -> bool {
        matches!(self, Self::Stored(p) if p.is_null())
    }

    pub(crate) fn offset(self, amount: i64) -> Result<Self, SourceHeapError> {
        match self {
            Self::Stored(p) => p.offset(amount).map(Self::Stored),
            Self::Relative(p) => p
                .checked_add(amount)
                .map(Self::Relative)
                .ok_or(SourceHeapError::PointerOffsetOverflow),
        }
    }

    pub(crate) fn stored(self) -> Result<SourceMutPointer<i8>, SourceHeapError> {
        match self {
            Self::Stored(p) => Ok(p),
            Self::Relative(_) => Err(SourceHeapError::PointerAllocationMismatch),
        }
    }

    pub(crate) fn relative(self) -> Result<i64, SourceHeapError> {
        match self {
            Self::Relative(p) => Ok(p),
            Self::Stored(_) => Err(SourceHeapError::PointerAllocationMismatch),
        }
    }
}

impl SourceLineBuffer<'_> {
    pub(crate) fn start(&self) -> SourceLineCursor {
        match self {
            Self::Pointer(p) => SourceLineCursor::Stored(*p),
            Self::Borrowed(_) => SourceLineCursor::Relative(0),
        }
    }

    pub(crate) fn distance(&self, cursor: SourceLineCursor) -> Result<i64, SourceHeapError> {
        match (self, cursor) {
            (Self::Pointer(base), SourceLineCursor::Stored(p)) => p.difference(*base),
            (Self::Borrowed(_), SourceLineCursor::Relative(p)) => Ok(p),
            _ => Err(SourceHeapError::PointerAllocationMismatch),
        }
    }

    pub(crate) fn bytes<'a>(
        &'a self,
        heap: &'a SourceHeap,
        cursor: SourceLineCursor,
    ) -> Result<&'a [i8], SourceHeapError> {
        // INCHI✔️❗: FindToken source445 searches p before using szLine.
        // Stored cursors resolve directly, without admitting an unused base.
        // Borrowed cursors select the actual array, no heap ID or snapshot.
        // O(1), checked transport cost stays OPEN under C542.
        match (self, cursor) {
            (Self::Pointer(_), SourceLineCursor::Stored(p)) => heap.slice(p.as_const()),
            (Self::Borrowed(bytes), SourceLineCursor::Relative(p)) => {
                let p = usize::try_from(p).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                bytes.get(p..).ok_or(SourceHeapError::PointerOutOfBounds)
            }
            _ => Err(SourceHeapError::PointerAllocationMismatch),
        }
    }

    pub(crate) fn bytes_mut<'a>(
        &'a mut self,
        heap: &'a mut SourceHeap,
        cursor: SourceLineCursor,
    ) -> Result<&'a mut [i8], SourceHeapError> {
        match (self, cursor) {
            (Self::Pointer(_), SourceLineCursor::Stored(p)) => heap.slice_mut(p),
            (Self::Borrowed(bytes), SourceLineCursor::Relative(p)) => {
                let p = usize::try_from(p).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                bytes
                    .get_mut(p..)
                    .ok_or(SourceHeapError::PointerOutOfBounds)
            }
            _ => Err(SourceHeapError::PointerAllocationMismatch),
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn source_read_line_raw(
    heap: &mut SourceHeap,
    target: &mut SourceLineBuffer<'_>,
    position: SourceLineCursor,
    len: i32,
    ios: &mut INCHI_IOSTREAM,
    tabs: bool,
) -> Result<bool, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:324 inchi_ios_str_gets
    // INCHI✔️❌: char* inchi_ios_str_gets(char* szLine, int len, INCHI_IOSTREAM* f)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int  length = 0, c = 0;
    // INCHI✔️❌:     if (--len < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     while (length < len && EOF != (c = inchi_ios_str_getc(f)))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         szLine[length++] = (char)c;
    // INCHI✔️❌:         if (c == '\n')
    // INCHI✔️❌:         {
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (!length && EOF == c)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     szLine[length] = '\0';
    // INCHI✔️❌:
    // INCHI✔️❌:     return szLine;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_ios_str_gets
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:354 inchi_ios_str_getsTab
    // INCHI✔️❌: char* inchi_ios_str_getsTab(char* szLine, int len, INCHI_IOSTREAM* f)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int  length = 0, c = 0;
    // INCHI✔️❌:     if (--len < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     while (length < len && EOF != (c = inchi_ios_str_getc(f)))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (c == '\t')
    // INCHI✔️❌:         {
    // INCHI✔️❌:             c = '\n';
    // INCHI✔️❌:         }
    // INCHI✔️❌:         szLine[length++] = (char)c;
    // INCHI✔️❌:         if (c == '\n')
    // INCHI✔️❌:         {
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (!length && EOF == c)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     szLine[length] = '\0';
    // INCHI✔️❌:
    // INCHI✔️❌:     return szLine;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_ios_str_getsTab
    // INCHI✔️❗: One reached-byte loop for the real pointer/array target;
    // preserve len/EOF/getc-before-write/tab/NUL order, no buffer or clone.
    // INCHI-PERF-AUXINFO-FIXED-LINES: checked/dispatch costs stay OPEN.

    let len = len
        .checked_sub(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    if len < 0 {
        return Ok(false);
    }

    let mut length = 0_i32;
    let mut c = 0_i32;
    while length < len {
        c = inchi_ios_str_getc(heap, ios)?;
        if c == SOURCE_EOF {
            break;
        }
        if tabs && c == i32::from(b'\t') {
            c = i32::from(b'\n');
        }
        let index = usize::try_from(length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        *target
            .bytes_mut(heap, position)?
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = c as i8;
        length += 1;
        if c == i32::from(b'\n') {
            break;
        }
    }
    if length == 0 && c == SOURCE_EOF {
        return Ok(false);
    }
    let index = usize::try_from(length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    *target
        .bytes_mut(heap, position)?
        .get_mut(index)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
    Ok(true)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn source_read_line(
    heap: &mut SourceHeap,
    target: &mut SourceLineBuffer<'_>,
    position: SourceLineCursor,
    len: i32,
    ios: &mut INCHI_IOSTREAM,
    too_long_line: &mut i32,
    tabs: bool,
    repeat_empty: bool,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:386 inchi_ios_gets
    // INCHI✔️❌: int inchi_ios_gets(char* szLine,
    // INCHI✔️❌:     int len,
    // INCHI✔️❌:     INCHI_IOSTREAM* f,
    // INCHI✔️❌:     int* bTooLongLine)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int  length;
    // INCHI✔️❌:     char* p;
    // INCHI✔️❌:     do
    // INCHI✔️❌:     {
    // INCHI✔️❌:         p = inchi_ios_str_gets(szLine, len - 1, f);
    // INCHI✔️❌:         if (!p)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *bTooLongLine = 0;
    // INCHI✔️❌:             return -1; /* end of file or cannot read */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         szLine[len - 1] = '\0';
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         *bTooLongLine = !strchr( szLine, '\n' );
    // INCHI✔️❌:         */
    // INCHI✔️❌:         p = strchr(szLine, '\n');
    // INCHI✔️❌:         *bTooLongLine = (!p && ((int)strlen(szLine)) == len - 2);
    // INCHI✔️❌:         lrtrim(szLine, &length);
    // INCHI✔️❌:     } while (!length);
    // INCHI✔️❌:
    // INCHI✔️❌:     return length;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_ios_gets
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:420 inchi_ios_getsTab
    // INCHI✔️❌: int inchi_ios_getsTab(char* szLine,
    // INCHI✔️❌:     int len,
    // INCHI✔️❌:     INCHI_IOSTREAM* f,
    // INCHI✔️❌:     int* bTooLongLine)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int  length;
    // INCHI✔️❌:     char* p;
    // INCHI✔️❌:     do
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         p = inchi_ios_str_getsTab(szLine, len - 1, f);
    // INCHI✔️❌:         if (!p)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *bTooLongLine = 0;
    // INCHI✔️❌:             return -1; /* end of file or cannot read */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         szLine[len - 1] = '\0';
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         *bTooLongLine = !strchr( szLine, '\n' );
    // INCHI✔️❌:         */
    // INCHI✔️❌:         p = strchr(szLine, '\n');
    // INCHI✔️❌:         *bTooLongLine = (!p && ((int)strlen(szLine)) == len - 2);
    // INCHI✔️❌:         lrtrim(szLine, &length);
    // INCHI✔️❌:
    // INCHI✔️❌:     } while (!length);
    // INCHI✔️❌:
    // INCHI✔️❌:     return length;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_ios_getsTab
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:451 inchi_ios_getsTab1
    // INCHI✔️❌: int inchi_ios_getsTab1(char* szLine,
    // INCHI✔️❌:     int len,
    // INCHI✔️❌:     INCHI_IOSTREAM* f,
    // INCHI✔️❌:     int* bTooLongLine)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int  length;
    // INCHI✔️❌:     char* p;
    // INCHI✔️❌:
    // INCHI✔️❌:     p = inchi_ios_str_getsTab(szLine, len - 1, f);
    // INCHI✔️❌:     if (!p)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *bTooLongLine = 0;
    // INCHI✔️❌:         return -1; /* end of file or cannot read */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     szLine[len - 1] = '\0';
    // INCHI✔️❌:     p = strchr(szLine, '\n');
    // INCHI✔️❌:     *bTooLongLine = (!p && ((int)strlen(szLine)) == len - 2);
    // INCHI✔️❌:     lrtrim(szLine, &length);
    // INCHI✔️❌:
    // INCHI✔️❌:     return length;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: inchi_ios_getsTab1
    // INCHI✔️❗: Preserve source raw-read/EOF/terminator/too-long/trim order.
    // Only source gets/getsTab retry empty lines; getsTab1 returns immediately.
    // INCHI-PERF-AUXINFO-FIXED-LINES: no allocation/copy/extra scan, costs OPEN.
    let helper_len = len
        .checked_sub(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let content_limit = len
        .checked_sub(2)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    loop {
        if !source_read_line_raw(heap, target, position, helper_len, ios, tabs)? {
            *too_long_line = 0;
            return Ok(-1);
        }
        let terminator_index =
            usize::try_from(helper_len).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        *target
            .bytes_mut(heap, position)?
            .get_mut(terminator_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
        *too_long_line = source_line_too_long(target.bytes(heap, position)?, content_limit)?;
        let length = crate::source::base::util::lrtrim_bytes(target.bytes_mut(heap, position)?)?;
        if !repeat_empty || length != 0 {
            return Ok(length);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{
        INCHI_IOS_STRING, SourceFile, SourceFormatArgument, SourceMutPointer, SourceVaList,
        SourceVoid,
    };
    use crate::test_support::allocate_source_fixture;

    #[test]
    fn inchi_string_printers_reserve_nul_at_growth_boundaries() {
        let block = INCHI_ADD_STR_LEN as usize;
        for printer in 0..3 {
            for prefix in ["", "prefix:"] {
                for length in [block - 1, block, block + 1, 2 * block] {
                    let mut heap = SourceHeap::default();
                    let stdout = heap.allocate(vec![SourceFile::default()]).unwrap();
                    let format = source_format(&mut heap, "%s");
                    let payload = "x".repeat(length);
                    let text = source_format(&mut heap, &payload);
                    let arguments = SourceVaList {
                        arguments: vec![SourceFormatArgument::Bytes(text.as_const())],
                        position: 0,
                    };
                    let mut stream = INCHI_IOSTREAM {
                        type_: INCHI_IOS_TYPE_STRING as i32,
                        ..Default::default()
                    };
                    if !prefix.is_empty() {
                        stream.s.pStr = source_format(&mut heap, prefix);
                        stream.s.nUsedLength = prefix.len() as i32;
                        stream.s.nAllocatedLength = prefix.len() as i32 + 1;
                    }
                    let result = match printer {
                        0 => inchi_ios_print(
                            &mut heap,
                            Some(&mut stream),
                            stdout,
                            format.as_const(),
                            &arguments,
                        ),
                        1 => inchi_ios_print_nodisplay(
                            &mut heap,
                            Some(&mut stream),
                            stdout,
                            format.as_const(),
                            &arguments,
                        ),
                        _ => inchi_ios_eprint(
                            &mut heap,
                            Some(&mut stream),
                            format.as_const(),
                            &arguments,
                        ),
                    };
                    assert_eq!(
                        result,
                        Ok(length as i32),
                        "printer={printer}, prefix={prefix}, length={length}"
                    );
                    let expected = format!("{prefix}{payload}\0");
                    let bytes = heap.slice(stream.s.pStr.as_const()).unwrap();
                    assert_eq!(stream.s.nUsedLength as usize, expected.len() - 1);
                    assert_eq!(stream.s.nAllocatedLength as usize, bytes.len());
                    assert_eq!(
                        bytes[..expected.len()]
                            .iter()
                            .map(|x| *x as u8)
                            .collect::<Vec<_>>(),
                        expected.as_bytes()
                    );
                    assert_eq!(stream.s.nPtr, 0);
                }
            }
        }
    }

    #[test]
    fn source_port__ichi_io__inchi_strbuf_addline__line_1635() {
        fn stream(heap: &mut SourceHeap, value: &[u8]) -> INCHI_IOSTREAM {
            let data = heap
                .allocate_model_storage(value.iter().map(|byte| *byte as i8).collect())
                .unwrap();
            INCHI_IOSTREAM {
                type_: INCHI_IOS_TYPE_STRING as i32,
                s: INCHI_IOS_STRING {
                    pStr: data,
                    nAllocatedLength: value.len() as i32,
                    nUsedLength: value.len() as i32,
                    nPtr: 0,
                },
                ..INCHI_IOSTREAM::default()
            }
        }
        fn bytes(heap: &SourceHeap, buffer: &INCHI_IOS_STRING) -> Vec<u8> {
            heap.slice(buffer.pStr.as_const()).unwrap()[..buffer.nUsedLength as usize]
                .iter()
                .map(|byte| *byte as u8)
                .collect()
        }

        let mut heap = SourceHeap::default();
        let mut buffer = INCHI_IOS_STRING::default();
        assert!(inchi_strbuf_init(&mut heap, &mut buffer, 4, 4).unwrap() > 0);
        heap.slice_mut(buffer.pStr).unwrap()[..4]
            .copy_from_slice(&[b'p' as i8, b'r' as i8, b'e' as i8, 0]);
        buffer.nUsedLength = 3;

        let mut crlf = stream(&mut heap, b"a\r\nrest");
        assert_eq!(
            inchi_strbuf_addline(&mut heap, &mut buffer, Some(&mut crlf), 1, 1),
            Ok(5)
        );
        assert_eq!(bytes(&heap, &buffer), b"prea\n");
        assert_eq!(crlf.s.nPtr, 3);

        let mut no_lf = stream(&mut heap, b"b\ntail");
        assert_eq!(
            inchi_strbuf_addline(&mut heap, &mut buffer, Some(&mut no_lf), 0, 0),
            Ok(6)
        );
        assert_eq!(bytes(&heap, &buffer), b"prea\nb");

        let mut eof = stream(&mut heap, b"xy");
        assert_eq!(
            inchi_strbuf_addline(&mut heap, &mut buffer, Some(&mut eof), 1, 0),
            Ok(-1)
        );
        assert_eq!(bytes(&heap, &buffer), b"prea\nbxy");

        let mut empty_buffer = INCHI_IOS_STRING::default();
        inchi_strbuf_init(&mut heap, &mut empty_buffer, 2, 2).unwrap();
        let mut empty_crlf = stream(&mut heap, b"\r\n");
        assert_eq!(
            inchi_strbuf_addline(&mut heap, &mut empty_buffer, Some(&mut empty_crlf), 1, 1),
            Ok(2)
        );
        assert_eq!(bytes(&heap, &empty_buffer), b"\r\n");

        let mut bare_lf = stream(&mut heap, b"\n");
        assert_eq!(
            inchi_strbuf_addline(&mut heap, &mut empty_buffer, Some(&mut bare_lf), 0, 0),
            Ok(2)
        );
        assert_eq!(bytes(&heap, &empty_buffer), b"\r\n");
        assert_eq!(
            inchi_strbuf_addline(&mut heap, &mut empty_buffer, None, 0, 0),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichi_io__inchi_strbuf_init__line_1370() {
        let mut heap = SourceHeap::default();
        let old = heap.allocate_model_storage(vec![b'X' as i8, 0]).unwrap();
        let mut explicit = INCHI_IOS_STRING {
            pStr: old,
            nAllocatedLength: 91,
            nUsedLength: 73,
            nPtr: 55,
        };
        assert_eq!(inchi_strbuf_init(&mut heap, &mut explicit, 4, 7), Ok(4));
        assert_ne!(explicit.pStr, old);
        assert_eq!(explicit.nAllocatedLength, 4);
        assert_eq!(explicit.nUsedLength, 0);
        assert_eq!(explicit.nPtr, 7);
        assert_eq!(heap.slice(explicit.pStr.as_const()).unwrap(), &[0; 4]);
        assert_eq!(heap.slice(old.as_const()).unwrap(), &[b'X' as i8, 0]);

        let mut defaults = INCHI_IOS_STRING::default();
        assert_eq!(
            inchi_strbuf_init(&mut heap, &mut defaults, 0, -1),
            Ok(INCHI_STRBUF_INITIAL_SIZE as i32)
        );
        assert_eq!(defaults.nAllocatedLength, INCHI_STRBUF_INITIAL_SIZE as i32);
        assert_eq!(defaults.nUsedLength, 0);
        assert_eq!(defaults.nPtr, INCHI_STRBUF_SIZE_INCREMENT as i32);
        assert_eq!(
            heap.slice(defaults.pStr.as_const()).unwrap(),
            vec![0; INCHI_STRBUF_INITIAL_SIZE as usize]
        );

        let mut failing_heap = SourceHeap::default();
        let failed_old = failing_heap
            .allocate_model_storage(vec![b'Y' as i8, 0])
            .unwrap();
        let mut failed = INCHI_IOS_STRING {
            pStr: failed_old,
            nAllocatedLength: 9,
            nUsedLength: 8,
            nPtr: 7,
        };
        failing_heap.fail_after_allocations(0);
        assert_eq!(
            inchi_strbuf_init(&mut failing_heap, &mut failed, -3, 5),
            Ok(-1)
        );
        assert_eq!(failed, INCHI_IOS_STRING::default());
        assert_eq!(
            failing_heap.slice(failed_old.as_const()).unwrap(),
            &[b'Y' as i8, 0]
        );
    }

    #[test]
    fn source_port__ichi_io__inchi_strbuf_close__line_1422() {
        let mut heap = SourceHeap::default();
        assert_eq!(inchi_strbuf_close(&mut heap, None), Ok(()));

        let mut empty = INCHI_IOS_STRING {
            nAllocatedLength: 11,
            nUsedLength: 7,
            nPtr: 3,
            ..INCHI_IOS_STRING::default()
        };
        assert_eq!(inchi_strbuf_close(&mut heap, Some(&mut empty)), Ok(()));
        assert_eq!(empty, INCHI_IOS_STRING::default());

        let allocation = heap.allocate(vec![b'A' as i8, 0]).unwrap();
        let mut populated = INCHI_IOS_STRING {
            pStr: allocation,
            nAllocatedLength: 2,
            nUsedLength: 1,
            nPtr: 19,
        };
        assert_eq!(inchi_strbuf_close(&mut heap, Some(&mut populated)), Ok(()));
        assert_eq!(populated, INCHI_IOS_STRING::default());
        assert_eq!(
            heap.slice(allocation.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__ichi_io__inchi_strbuf_create_copy__line_1438() {
        let mut heap = SourceHeap::default();
        let source_string = heap
            .allocate_model_storage(vec![b'A' as i8, b'B' as i8, 0, 91, 92, 93])
            .unwrap();
        let source = INCHI_IOS_STRING {
            pStr: source_string,
            nAllocatedLength: 6,
            nUsedLength: 2,
            nPtr: 17,
        };
        let old_destination = heap.allocate_model_storage(vec![b'X' as i8, 0]).unwrap();
        let mut destination = INCHI_IOS_STRING {
            pStr: old_destination,
            nAllocatedLength: 71,
            nUsedLength: 72,
            nPtr: 73,
        };
        assert_eq!(
            inchi_strbuf_create_copy(&mut heap, &mut destination, &source),
            Ok(0)
        );
        assert_ne!(destination.pStr, source.pStr);
        assert_ne!(destination.pStr, old_destination);
        assert_eq!(
            (
                destination.nAllocatedLength,
                destination.nUsedLength,
                destination.nPtr,
            ),
            (6, 2, 17)
        );
        assert_eq!(heap.slice(destination.pStr.as_const()).unwrap(), &[0; 6]);
        assert_eq!(
            heap.slice(source_string.as_const()).unwrap(),
            &[b'A' as i8, b'B' as i8, 0, 91, 92, 93]
        );
        assert_eq!(
            heap.slice(old_destination.as_const()).unwrap(),
            &[b'X' as i8, 0]
        );

        let source_without_storage = INCHI_IOS_STRING {
            pStr: SourceMutPointer::null(),
            nAllocatedLength: 3,
            nUsedLength: i32::MIN,
            nPtr: i32::MAX,
        };
        let mut destination_without_storage = INCHI_IOS_STRING::default();
        assert_eq!(
            inchi_strbuf_create_copy(
                &mut heap,
                &mut destination_without_storage,
                &source_without_storage,
            ),
            Ok(0)
        );
        assert_eq!(
            heap.slice(destination_without_storage.pStr.as_const())
                .unwrap(),
            &[0; 3]
        );
        assert_eq!(destination_without_storage.nUsedLength, i32::MIN);
        assert_eq!(destination_without_storage.nPtr, i32::MAX);

        let mut zero_destination = INCHI_IOS_STRING {
            nAllocatedLength: 7,
            nUsedLength: 8,
            nPtr: 9,
            ..INCHI_IOS_STRING::default()
        };
        assert_eq!(
            inchi_strbuf_create_copy(
                &mut heap,
                &mut zero_destination,
                &INCHI_IOS_STRING::default(),
            ),
            Ok(0)
        );
        assert!(!zero_destination.pStr.is_null());
        assert!(
            heap.slice(zero_destination.pStr.as_const())
                .unwrap()
                .is_empty()
        );
        assert_eq!(
            zero_destination,
            INCHI_IOS_STRING {
                pStr: zero_destination.pStr,
                ..INCHI_IOS_STRING::default()
            }
        );

        let mut failing_heap = SourceHeap::default();
        let failed_old = failing_heap
            .allocate_model_storage(vec![b'Y' as i8, 0])
            .unwrap();
        let mut failed = INCHI_IOS_STRING {
            pStr: failed_old,
            nAllocatedLength: 31,
            nUsedLength: 32,
            nPtr: 33,
        };
        let failed_source = INCHI_IOS_STRING {
            nAllocatedLength: 4,
            nUsedLength: 2,
            nPtr: 1,
            ..INCHI_IOS_STRING::default()
        };
        failing_heap.fail_after_allocations(0);
        assert_eq!(
            inchi_strbuf_create_copy(&mut failing_heap, &mut failed, &failed_source),
            Ok(-1)
        );
        assert!(failed.pStr.is_null());
        assert_eq!(
            (failed.nAllocatedLength, failed.nUsedLength, failed.nPtr),
            (31, 32, 33)
        );
        assert_eq!(
            failing_heap.slice(failed_old.as_const()).unwrap(),
            &[b'Y' as i8, 0]
        );

        let negative_old = heap.allocate_model_storage(vec![7_i8]).unwrap();
        let mut negative_destination = INCHI_IOS_STRING {
            pStr: negative_old,
            nAllocatedLength: 41,
            nUsedLength: 42,
            nPtr: 43,
        };
        let negative_source = INCHI_IOS_STRING {
            nAllocatedLength: -1,
            nUsedLength: 51,
            nPtr: 52,
            ..INCHI_IOS_STRING::default()
        };
        assert_eq!(
            inchi_strbuf_create_copy(&mut heap, &mut negative_destination, &negative_source),
            Ok(-1)
        );
        assert!(negative_destination.pStr.is_null());
        assert_eq!(
            (
                negative_destination.nAllocatedLength,
                negative_destination.nUsedLength,
                negative_destination.nPtr,
            ),
            (41, 42, 43)
        );
        assert_eq!(heap.slice(negative_old.as_const()).unwrap(), &[7]);
    }

    #[test]
    fn source_port__ichi_io__inchi_strbuf_reset__line_1403() {
        let mut heap = SourceHeap::default();
        assert_eq!(inchi_strbuf_reset(&mut heap, None), Ok(()));

        let mut empty = INCHI_IOS_STRING {
            nAllocatedLength: 17,
            nUsedLength: 9,
            nPtr: 4,
            ..INCHI_IOS_STRING::default()
        };
        assert_eq!(inchi_strbuf_reset(&mut heap, Some(&mut empty)), Ok(()));
        assert!(empty.pStr.is_null());
        assert_eq!(empty.nAllocatedLength, 17);
        assert_eq!(empty.nUsedLength, 0);
        assert_eq!(empty.nPtr, 0);

        let allocation = heap.allocate(vec![65_i8, 66, 67, 0, 99]).unwrap();
        let mut populated = INCHI_IOS_STRING {
            pStr: allocation,
            nAllocatedLength: 5,
            nUsedLength: 3,
            nPtr: 2,
        };
        assert_eq!(inchi_strbuf_reset(&mut heap, Some(&mut populated)), Ok(()));
        assert_eq!(populated.pStr, allocation);
        assert_eq!(populated.nAllocatedLength, 5);
        assert_eq!(populated.nUsedLength, 0);
        assert_eq!(populated.nPtr, 0);
        assert_eq!(
            heap.slice(allocation.as_const()).unwrap(),
            &[0, 66, 67, 0, 99]
        );
    }

    #[test]
    fn source_port__ichi_io__inchi_strbuf_update__line_1459() {
        let mut heap = SourceHeap::default();
        assert_eq!(inchi_strbuf_update(&mut heap, None, 8), Ok(-1));

        let pointer = heap.allocate(vec![1_i8, 2, 3, 0, 77, 88]).unwrap();
        let mut buffer = INCHI_IOS_STRING {
            pStr: pointer,
            nAllocatedLength: 6,
            nUsedLength: 3,
            nPtr: 4,
        };
        assert_eq!(inchi_strbuf_update(&mut heap, Some(&mut buffer), 0), Ok(6));
        assert_eq!(inchi_strbuf_update(&mut heap, Some(&mut buffer), -7), Ok(6));
        assert_eq!(buffer.pStr, pointer);
        assert_eq!(inchi_strbuf_update(&mut heap, Some(&mut buffer), 2), Ok(6));
        assert_eq!(buffer.pStr, pointer);
        assert_eq!(
            heap.slice(pointer.as_const()).unwrap(),
            &[1, 2, 3, 0, 77, 88]
        );

        assert_eq!(inchi_strbuf_update(&mut heap, Some(&mut buffer), 3), Ok(10));
        assert_ne!(buffer.pStr, pointer);
        assert_eq!(buffer.nAllocatedLength, 10);
        assert_eq!(buffer.nUsedLength, 3);
        assert_eq!(buffer.nPtr, 4);
        assert_eq!(
            &heap.slice(buffer.pStr.as_const()).unwrap()[..3],
            &[1, 2, 3]
        );
        assert_eq!(&heap.slice(buffer.pStr.as_const()).unwrap()[3..], &[0; 7]);
        assert_eq!(
            heap.slice(pointer.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let second_pointer = buffer.pStr;
        assert_eq!(inchi_strbuf_update(&mut heap, Some(&mut buffer), 9), Ok(19));
        assert_eq!(buffer.nAllocatedLength, 19);
        assert_ne!(buffer.pStr, second_pointer);
        assert_eq!(
            &heap.slice(buffer.pStr.as_const()).unwrap()[..3],
            &[1, 2, 3]
        );

        let mut failing_heap = SourceHeap::default();
        let old = failing_heap.allocate(vec![9_i8, 8, 0, 7]).unwrap();
        let mut failed = INCHI_IOS_STRING {
            pStr: old,
            nAllocatedLength: 4,
            nUsedLength: 2,
            nPtr: 8,
        };
        let before = failed.clone();
        failing_heap.fail_after_allocations(0);
        assert_eq!(
            inchi_strbuf_update(&mut failing_heap, Some(&mut failed), 2),
            Ok(-1)
        );
        assert_eq!(failed, before);
        assert_eq!(failing_heap.slice(old.as_const()).unwrap(), &[9, 8, 0, 7]);

        let mut null_buffer = INCHI_IOS_STRING {
            nAllocatedLength: 0,
            nUsedLength: 0,
            nPtr: 5,
            ..INCHI_IOS_STRING::default()
        };
        assert_eq!(
            inchi_strbuf_update(&mut heap, Some(&mut null_buffer), 1),
            Ok(5)
        );
        assert_eq!(heap.slice(null_buffer.pStr.as_const()).unwrap(), &[0; 5]);
    }

    #[test]
    fn source_port__ichi_io__inchi_strbuf_printf__line_1507() {
        let mut heap = SourceHeap::default();
        let format = source_format(&mut heap, "%s:%+d");
        let value = source_format(&mut heap, "xy");
        let arguments = SourceVaList {
            arguments: vec![
                SourceFormatArgument::Bytes(value.as_const()),
                SourceFormatArgument::Signed(7),
            ],
            ..SourceVaList::default()
        };
        assert_eq!(
            inchi_strbuf_printf(&mut heap, None, format.as_const(), &arguments),
            Ok(-1)
        );

        let mut buffer = INCHI_IOS_STRING::default();
        assert_eq!(inchi_strbuf_init(&mut heap, &mut buffer, 4, 3), Ok(4));
        assert_eq!(
            inchi_strbuf_printf(&mut heap, Some(&mut buffer), format.as_const(), &arguments,),
            Ok(5)
        );
        assert_eq!(buffer.nUsedLength, 5);
        assert_eq!(
            &heap.slice(buffer.pStr.as_const()).unwrap()[..6],
            &[
                b'x' as i8, b'y' as i8, b':' as i8, b'+' as i8, b'7' as i8, 0
            ]
        );

        let suffix = source_format(&mut heap, "/%02u");
        let suffix_arguments = SourceVaList {
            arguments: vec![SourceFormatArgument::Unsigned(3)],
            ..SourceVaList::default()
        };
        assert_eq!(
            inchi_strbuf_printf(
                &mut heap,
                Some(&mut buffer),
                suffix.as_const(),
                &suffix_arguments,
            ),
            Ok(3)
        );
        assert_eq!(buffer.nUsedLength, 8);
        assert_eq!(
            &heap.slice(buffer.pStr.as_const()).unwrap()[..9],
            &[120, 121, 58, 43, 55, 47, 48, 51, 0]
        );

        let invalid = source_format(&mut heap, "%*d");
        let invalid_arguments = SourceVaList {
            arguments: vec![
                SourceFormatArgument::Signed(-1),
                SourceFormatArgument::Signed(5),
            ],
            ..SourceVaList::default()
        };
        let before = buffer.clone();
        assert_eq!(
            inchi_strbuf_printf(
                &mut heap,
                Some(&mut buffer),
                invalid.as_const(),
                &invalid_arguments,
            ),
            Ok(0)
        );
        assert_eq!(buffer, before);

        let mut failing_heap = SourceHeap::default();
        let old = failing_heap.allocate(vec![b'A' as i8, 0]).unwrap();
        let format = source_format(&mut failing_heap, "%s");
        let value = source_format(&mut failing_heap, "long");
        let arguments = SourceVaList {
            arguments: vec![SourceFormatArgument::Bytes(value.as_const())],
            ..SourceVaList::default()
        };
        let mut failed = INCHI_IOS_STRING {
            pStr: old,
            nAllocatedLength: 2,
            nUsedLength: 1,
            nPtr: 2,
        };
        let before = failed.clone();
        failing_heap.fail_after_allocations(0);
        assert_eq!(
            inchi_strbuf_printf(
                &mut failing_heap,
                Some(&mut failed),
                format.as_const(),
                &arguments,
            ),
            Err(SourceHeapError::AllocationFailed)
        );
        assert_eq!(failed, before);
        assert_eq!(
            failing_heap.slice(old.as_const()).unwrap(),
            &[b'A' as i8, 0]
        );
    }

    #[test]
    fn source_port__ichi_io__inchi_ios_print__line_477() {
        let mut heap = SourceHeap::default();
        let stdout = heap.allocate(vec![SourceFile::default()]).unwrap();
        let format = source_format(&mut heap, "%s:%04d");
        let value = source_format(&mut heap, "x");
        let arguments = SourceVaList {
            arguments: vec![
                SourceFormatArgument::Bytes(value.as_const()),
                SourceFormatArgument::Signed(7),
            ],
            position: 0,
        };
        assert_eq!(
            inchi_ios_print(&mut heap, None, stdout, format.as_const(), &arguments),
            Ok(-1)
        );

        let mut string = INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(
            inchi_ios_print(
                &mut heap,
                Some(&mut string),
                stdout,
                format.as_const(),
                &arguments,
            ),
            Ok(6)
        );
        let first = string.s.pStr;
        assert_eq!(string.s.nUsedLength, 6);
        assert_eq!(string.s.nAllocatedLength, INCHI_ADD_STR_LEN as i32);
        assert_eq!(
            &heap.slice(first.as_const()).unwrap()[..7],
            b"x:0007\0".map(|byte| byte as i8).as_slice()
        );
        assert_eq!(
            inchi_ios_print(
                &mut heap,
                Some(&mut string),
                stdout,
                format.as_const(),
                &arguments,
            ),
            Ok(6)
        );
        assert_eq!(string.s.pStr, first);
        assert_eq!(string.s.nUsedLength, 12);
        assert_eq!(
            &heap.slice(first.as_const()).unwrap()[..13],
            b"x:0007x:0007\0".map(|byte| byte as i8).as_slice()
        );

        let old = heap.allocate(vec![b'A' as i8, 0]).unwrap();
        let mut exact = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: old,
                nAllocatedLength: 2,
                nUsedLength: 1,
                nPtr: 31,
            },
            type_: INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        let literal = source_format(&mut heap, "B");
        assert_eq!(
            inchi_ios_print(
                &mut heap,
                Some(&mut exact),
                stdout,
                literal.as_const(),
                &SourceVaList::default(),
            ),
            Ok(1)
        );
        assert_ne!(exact.s.pStr, old);
        assert_eq!(exact.s.nPtr, 31);
        assert_eq!(
            heap.slice(old.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            &heap.slice(exact.s.pStr.as_const()).unwrap()[..3],
            &[65, 66, 0]
        );

        let file = heap.allocate(vec![SourceFile::default()]).unwrap();
        let mut file_ios = INCHI_IOSTREAM {
            f: file,
            type_: INCHI_IOS_TYPE_FILE as i32,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(
            inchi_ios_print(
                &mut heap,
                Some(&mut file_ios),
                stdout,
                literal.as_const(),
                &SourceVaList::default(),
            ),
            Ok(1)
        );
        assert_eq!(heap.slice(file.as_const()).unwrap()[0].bytes, b"B");
        file_ios.f = SourceMutPointer::null();
        assert_eq!(
            inchi_ios_print(
                &mut heap,
                Some(&mut file_ios),
                stdout,
                literal.as_const(),
                &SourceVaList::default(),
            ),
            Ok(1)
        );
        assert_eq!(heap.slice(stdout.as_const()).unwrap()[0].bytes, b"B");
        heap.slice_mut(file).unwrap()[0].error = true;
        file_ios.f = file;
        assert_eq!(
            inchi_ios_print(
                &mut heap,
                Some(&mut file_ios),
                stdout,
                literal.as_const(),
                &SourceVaList::default(),
            ),
            Ok(-1)
        );

        file_ios.type_ = 99;
        assert_eq!(
            inchi_ios_print(
                &mut heap,
                Some(&mut file_ios),
                stdout,
                literal.as_const(),
                &SourceVaList::default(),
            ),
            Ok(0)
        );
        let invalid = source_format(&mut heap, "%q");
        exact.type_ = INCHI_IOS_TYPE_STRING as i32;
        let before = exact.clone();
        assert_eq!(
            inchi_ios_print(
                &mut heap,
                Some(&mut exact),
                stdout,
                invalid.as_const(),
                &SourceVaList::default(),
            ),
            Ok(-1)
        );
        assert_eq!(exact, before);

        let mut failing_heap = SourceHeap::default();
        let failing_format = failing_heap
            .allocate_model_storage(vec![b'Z' as i8, 0])
            .unwrap();
        let failing_stdout = failing_heap
            .allocate_model_storage(vec![SourceFile::default()])
            .unwrap();
        let mut failing = INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        failing_heap.fail_after_allocations(0);
        assert_eq!(
            inchi_ios_print(
                &mut failing_heap,
                Some(&mut failing),
                failing_stdout,
                failing_format.as_const(),
                &SourceVaList::default(),
            ),
            Ok(-1)
        );
        assert_eq!(
            failing,
            INCHI_IOSTREAM {
                type_: INCHI_IOS_TYPE_STRING as i32,
                ..INCHI_IOSTREAM::default()
            }
        );
    }

    fn source_format(heap: &mut SourceHeap, value: &str) -> SourceMutPointer<i8> {
        allocate_source_fixture(
            heap,
            value.bytes().map(|byte| byte as i8).chain([0]).collect(),
        )
    }

    fn get_max_length(
        heap: &mut SourceHeap,
        format: &str,
        arguments: Vec<SourceFormatArgument>,
    ) -> (i32, u64) {
        let format_pointer = source_format(heap, format);
        let mut arguments = SourceVaList {
            arguments,
            position: 0,
        };
        let result = GetMaxPrintfLength(heap, format_pointer.as_const(), &mut arguments).unwrap();
        inchi_free(heap, format_pointer).unwrap();
        (result, arguments.position)
    }

    #[test]
    fn source_port__ichi_io__getmaxprintflength__line_1065() {
        let mut heap = SourceHeap::default();
        assert_eq!(get_max_length(&mut heap, "ab%%cd", vec![]), (5, 0));
        assert_eq!(
            get_max_length(&mut heap, "%#08.3d", vec![SourceFormatArgument::Signed(7)],),
            (34, 1)
        );
        assert_eq!(
            get_max_length(
                &mut heap,
                "%*d",
                vec![
                    SourceFormatArgument::Signed(40),
                    SourceFormatArgument::Signed(7),
                ],
            ),
            (40, 2)
        );
        assert_eq!(
            get_max_length(
                &mut heap,
                "%*d",
                vec![
                    SourceFormatArgument::Signed(-1),
                    SourceFormatArgument::Signed(7),
                ],
            ),
            (-1, 1)
        );
        assert_eq!(
            get_max_length(
                &mut heap,
                "%.*d",
                vec![
                    SourceFormatArgument::Signed(-1),
                    SourceFormatArgument::Signed(7),
                ],
            ),
            (-1, 1)
        );

        let string = source_format(&mut heap, "abc");
        let empty = source_format(&mut heap, "");
        assert_eq!(
            get_max_length(
                &mut heap,
                "%s%5.2s%s",
                vec![
                    SourceFormatArgument::Bytes(string.as_const()),
                    SourceFormatArgument::Bytes(string.as_const()),
                    SourceFormatArgument::Bytes(empty.as_const()),
                ],
            ),
            (6, 3)
        );
        assert_eq!(
            get_max_length(
                &mut heap,
                "%c%5.1C",
                vec![
                    SourceFormatArgument::Signed(65),
                    SourceFormatArgument::Signed(66),
                ],
            ),
            (3, 2)
        );
        assert_eq!(
            get_max_length(
                &mut heap,
                "%d%i%u%x%X%o",
                vec![SourceFormatArgument::Signed(1); 6],
            ),
            (192, 6)
        );
        assert_eq!(
            get_max_length(
                &mut heap,
                "%e%f%g%G",
                vec![SourceFormatArgument::Float(1.0); 4],
            ),
            (128, 4)
        );
        assert_eq!(
            get_max_length(
                &mut heap,
                "%p%n",
                vec![
                    SourceFormatArgument::Pointer(std::ptr::null()),
                    SourceFormatArgument::Pointer(std::ptr::null()),
                ],
            ),
            (32, 2)
        );
        assert_eq!(
            get_max_length(
                &mut heap,
                "%hd%ld%lf%hs%hc%Ff%Nd%Ld",
                vec![
                    SourceFormatArgument::Signed(1),
                    SourceFormatArgument::Signed(2),
                    SourceFormatArgument::Float(3.0),
                    SourceFormatArgument::Bytes(string.as_const()),
                    SourceFormatArgument::Signed(65),
                    SourceFormatArgument::Float(4.0),
                    SourceFormatArgument::Signed(5),
                    SourceFormatArgument::Signed(6),
                ],
            ),
            (197, 8)
        );
        assert_eq!(get_max_length(&mut heap, "%q", vec![]), (-1, 0));
        assert_eq!(get_max_length(&mut heap, "%ls", vec![]), (-1, 0));
        assert_eq!(get_max_length(&mut heap, "%", vec![]), (-1, 0));
        inchi_free(&mut heap, string).unwrap();
        inchi_free(&mut heap, empty).unwrap();
    }

    #[test]
    fn source_port__ichi_io__inchi_vfprintf__line_828() {
        let mut heap = SourceHeap::default();
        let mut no_arguments = SourceVaList::default();
        assert_eq!(
            inchi_vfprintf(
                &mut heap,
                SourceMutPointer::null(),
                SourceConstPointer::null(),
                &mut no_arguments,
            ),
            Ok(0)
        );
        let empty_format = source_format(&mut heap, "");
        assert_eq!(
            inchi_vfprintf(
                &mut heap,
                SourceMutPointer::null(),
                empty_format.as_const(),
                &mut no_arguments,
            ),
            Ok(0)
        );

        let string = source_format(&mut heap, "abcdef");
        let format = source_format(
            &mut heap,
            "A:%+06d U:%#x S:%-5.3s F:%.2f E:%.1e G:%.3g %% %nZ\r",
        );
        let count = heap.allocate(vec![-1_i32]).unwrap();
        let file = heap
            .allocate(vec![SourceFile {
                bytes: b"xxxxxxxx".to_vec(),
                position: 2,
                ..SourceFile::default()
            }])
            .unwrap();
        let mut arguments = SourceVaList {
            arguments: vec![
                SourceFormatArgument::Signed(-12),
                SourceFormatArgument::Unsigned(42),
                SourceFormatArgument::Bytes(string.as_const()),
                SourceFormatArgument::Float(1.25),
                SourceFormatArgument::Float(12.0),
                SourceFormatArgument::Float(1234.0),
                SourceFormatArgument::MutSigned(count),
            ],
            position: 0,
        };
        let expected = "A:-00012 U:0x2a S:abc   F:1.25 E:1.2e+01 G:1.23e+03 % Z\r";
        assert_eq!(
            inchi_vfprintf(&mut heap, file, format.as_const(), &mut arguments),
            Ok(expected.len() as i32)
        );
        assert_eq!(arguments.position, 7);
        assert_eq!(
            heap.slice(count.as_const()).unwrap()[0],
            expected.find('Z').unwrap() as i32
        );
        let file_value = &heap.slice(file.as_const()).unwrap()[0];
        assert_eq!(&file_value.bytes[..2], b"xx");
        assert_eq!(&file_value.bytes[2..], expected.as_bytes());
        assert_eq!(file_value.position, (2 + expected.len()) as u64);

        let special_format = source_format(&mut heap, "%p %#o %.0u %hhd %lld %+f %F");
        let special_file = heap.allocate(vec![SourceFile::default()]).unwrap();
        let mut special_arguments = SourceVaList {
            arguments: vec![
                SourceFormatArgument::Pointer(std::ptr::null()),
                SourceFormatArgument::Unsigned(8),
                SourceFormatArgument::Unsigned(0),
                SourceFormatArgument::Signed(255),
                SourceFormatArgument::Signed(i64::MIN),
                SourceFormatArgument::Float(f64::INFINITY),
                SourceFormatArgument::Float(f64::NAN),
            ],
            position: 0,
        };
        let special_expected = "(nil) 010  -1 -9223372036854775808 +inf NAN";
        assert_eq!(
            inchi_vfprintf(
                &mut heap,
                special_file,
                special_format.as_const(),
                &mut special_arguments,
            ),
            Ok(special_expected.len() as i32)
        );
        assert_eq!(
            heap.slice(special_file.as_const()).unwrap()[0].bytes,
            special_expected.as_bytes()
        );

        let opaque = heap.allocate(vec![SourceVoid]).unwrap();
        let pointer_format = source_format(&mut heap, "%p");
        let pointer_file = heap.allocate(vec![SourceFile::default()]).unwrap();
        let rust_pointer = heap.slice(opaque.as_const()).unwrap().as_ptr();
        let pointer_expected = format!("{rust_pointer:p}");
        let mut pointer_arguments = SourceVaList {
            arguments: vec![SourceFormatArgument::Pointer(rust_pointer)],
            position: 0,
        };
        assert_eq!(
            inchi_vfprintf(
                &mut heap,
                pointer_file,
                pointer_format.as_const(),
                &mut pointer_arguments,
            ),
            Ok(pointer_expected.len() as i32)
        );

        assert_eq!(
            heap.slice(pointer_file.as_const()).unwrap()[0].bytes,
            pointer_expected.as_bytes()
        );

        let error_file = heap
            .allocate(vec![SourceFile {
                bytes: b"unchanged".to_vec(),
                error: true,
                ..SourceFile::default()
            }])
            .unwrap();
        let literal = source_format(&mut heap, "literal");
        assert_eq!(
            inchi_vfprintf(
                &mut heap,
                error_file,
                literal.as_const(),
                &mut SourceVaList::default(),
            ),
            Ok(-1)
        );
        assert_eq!(
            heap.slice(error_file.as_const()).unwrap()[0].bytes,
            b"unchanged"
        );
    }

    #[test]
    fn source_port__ichi_io__inchi_print_nodisplay__line_867() {
        let mut heap = SourceHeap::default();
        let format = source_format(&mut heap, "[%s] %d");
        let value = source_format(&mut heap, "hidden");
        let arguments = SourceVaList {
            arguments: vec![
                SourceFormatArgument::Bytes(value.as_const()),
                SourceFormatArgument::Signed(-17),
            ],
            position: 0,
        };
        let file = heap.allocate(vec![SourceFile::default()]).unwrap();
        let stdout = heap.allocate(vec![SourceFile::default()]).unwrap();

        assert_eq!(
            inchi_print_nodisplay(&mut heap, Some(file), stdout, format.as_const(), &arguments,),
            Ok(12)
        );
        assert_eq!(
            heap.slice(file.as_const()).unwrap()[0].bytes,
            b"[hidden] -17"
        );
        assert!(heap.slice(stdout.as_const()).unwrap()[0].bytes.is_empty());

        assert_eq!(
            inchi_print_nodisplay(&mut heap, None, stdout, format.as_const(), &arguments),
            Ok(12)
        );
        assert_eq!(
            heap.slice(stdout.as_const()).unwrap()[0].bytes,
            b"[hidden] -17"
        );

        heap.slice_mut(file).unwrap()[0].error = true;
        assert_eq!(
            inchi_print_nodisplay(&mut heap, Some(file), stdout, format.as_const(), &arguments,),
            Ok(-1)
        );
        assert_eq!(
            heap.slice(file.as_const()).unwrap()[0].bytes,
            b"[hidden] -17"
        );
        assert_eq!(arguments.position, 0);
    }

    #[test]
    fn source_port__ichi_io__inchi_ios_print_nodisplay__line_604() {
        let mut heap = SourceHeap::default();
        let format = source_format(&mut heap, "%s:%04d");
        let value = source_format(&mut heap, "x");
        let arguments = SourceVaList {
            arguments: vec![
                SourceFormatArgument::Bytes(value.as_const()),
                SourceFormatArgument::Signed(7),
            ],
            position: 0,
        };
        let stdout = heap.allocate(vec![SourceFile::default()]).unwrap();
        let mut string = INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(
            inchi_ios_print_nodisplay(
                &mut heap,
                Some(&mut string),
                stdout,
                format.as_const(),
                &arguments,
            ),
            Ok(6)
        );
        assert_eq!(string.s.nUsedLength, 6);
        assert_eq!(string.s.nAllocatedLength, INCHI_ADD_STR_LEN as i32);
        assert_eq!(
            &heap.slice(string.s.pStr.as_const()).unwrap()[..7],
            &[
                b'x' as i8, b':' as i8, b'0' as i8, b'0' as i8, b'0' as i8, b'7' as i8, 0
            ]
        );
        assert_eq!(
            inchi_ios_print_nodisplay(
                &mut heap,
                Some(&mut string),
                stdout,
                format.as_const(),
                &arguments,
            ),
            Ok(6)
        );
        assert_eq!(string.s.nUsedLength, 12);
        assert_eq!(
            &heap.slice(string.s.pStr.as_const()).unwrap()[..13],
            b"x:0007x:0007\0".map(|byte| byte as i8).as_slice()
        );

        let file = heap.allocate(vec![SourceFile::default()]).unwrap();
        let mut file_stream = INCHI_IOSTREAM {
            f: file,
            type_: INCHI_IOS_TYPE_FILE as i32,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(
            inchi_ios_print_nodisplay(
                &mut heap,
                Some(&mut file_stream),
                stdout,
                format.as_const(),
                &arguments,
            ),
            Ok(0)
        );
        assert_eq!(heap.slice(file.as_const()).unwrap()[0].bytes, b"x:0007");

        let mut null_file_stream = INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_FILE as i32,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(
            inchi_ios_print_nodisplay(
                &mut heap,
                Some(&mut null_file_stream),
                stdout,
                format.as_const(),
                &arguments,
            ),
            Ok(0)
        );
        assert_eq!(heap.slice(stdout.as_const()).unwrap()[0].bytes, b"x:0007");

        let mut none_stream = INCHI_IOSTREAM::default();
        assert_eq!(
            inchi_ios_print_nodisplay(
                &mut heap,
                Some(&mut none_stream),
                stdout,
                format.as_const(),
                &arguments,
            ),
            Ok(0)
        );
        assert_eq!(
            inchi_ios_print_nodisplay(&mut heap, None, stdout, format.as_const(), &arguments),
            Ok(-1)
        );

        let mut failing_heap = SourceHeap::default();
        let failing_format = source_format(&mut failing_heap, "abc");
        let failing_stdout = failing_heap
            .allocate_model_storage(vec![SourceFile::default()])
            .unwrap();
        let mut failing_stream = INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        failing_heap.fail_after_allocations(0);
        assert_eq!(
            inchi_ios_print_nodisplay(
                &mut failing_heap,
                Some(&mut failing_stream),
                failing_stdout,
                failing_format.as_const(),
                &SourceVaList::default(),
            ),
            Ok(-1)
        );
        assert_eq!(
            failing_stream,
            INCHI_IOSTREAM {
                type_: INCHI_IOS_TYPE_STRING as i32,
                ..INCHI_IOSTREAM::default()
            }
        );
    }

    #[test]
    fn source_port__ichi_io__inchi_ios_eprint__line_708() {
        let mut heap = SourceHeap::default();
        let format = source_format(&mut heap, "value=%d");
        let arguments = SourceVaList {
            arguments: vec![SourceFormatArgument::Signed(17)],
            position: 0,
        };
        assert_eq!(
            inchi_ios_eprint(&mut heap, None, format.as_const(), &arguments),
            Ok(-1)
        );

        let mut string_ios = INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(
            inchi_ios_eprint(
                &mut heap,
                Some(&mut string_ios),
                format.as_const(),
                &arguments,
            ),
            Ok(8)
        );
        assert_eq!(arguments.position, 0);
        assert_eq!(string_ios.s.nUsedLength, 8);
        assert_eq!(string_ios.s.nAllocatedLength, INCHI_ADD_STR_LEN as i32);
        let first_allocation = string_ios.s.pStr;
        assert_eq!(
            &heap.slice(first_allocation.as_const()).unwrap()[..9],
            &[
                b'v' as i8, b'a' as i8, b'l' as i8, b'u' as i8, b'e' as i8, b'=' as i8, b'1' as i8,
                b'7' as i8, 0,
            ]
        );
        let suffix = source_format(&mut heap, "/%s");
        let text = source_format(&mut heap, "ok");
        let suffix_arguments = SourceVaList {
            arguments: vec![SourceFormatArgument::Bytes(text.as_const())],
            position: 0,
        };
        assert_eq!(
            inchi_ios_eprint(
                &mut heap,
                Some(&mut string_ios),
                suffix.as_const(),
                &suffix_arguments,
            ),
            Ok(3)
        );
        assert_eq!(string_ios.s.pStr, first_allocation);
        assert_eq!(string_ios.s.nUsedLength, 11);
        assert_eq!(
            &heap.slice(first_allocation.as_const()).unwrap()[..12],
            &[
                b'v' as i8, b'a' as i8, b'l' as i8, b'u' as i8, b'e' as i8, b'=' as i8, b'1' as i8,
                b'7' as i8, b'/' as i8, b'o' as i8, b'k' as i8, 0,
            ]
        );

        let exact_old = heap.allocate(vec![b'A' as i8, 0]).unwrap();
        let mut exact = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: exact_old,
                nAllocatedLength: 2,
                nUsedLength: 1,
                nPtr: 0,
            },
            type_: INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        let literal = source_format(&mut heap, "B");
        assert_eq!(
            inchi_ios_eprint(
                &mut heap,
                Some(&mut exact),
                literal.as_const(),
                &SourceVaList::default(),
            ),
            Ok(1)
        );
        assert_ne!(exact.s.pStr, exact_old);
        assert_eq!(
            heap.slice(exact_old.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            &heap.slice(exact.s.pStr.as_const()).unwrap()[..3],
            &[65, 66, 0]
        );

        let invalid = source_format(&mut heap, "%q");
        let exact_before = exact.clone();
        assert_eq!(
            inchi_ios_eprint(
                &mut heap,
                Some(&mut exact),
                invalid.as_const(),
                &SourceVaList::default(),
            ),
            Ok(-1)
        );
        assert_eq!(exact, exact_before);

        let allocation_old = heap.allocate(vec![b'X' as i8, 0]).unwrap();
        let mut allocation_failure = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: allocation_old,
                nAllocatedLength: 2,
                nUsedLength: 1,
                nPtr: 0,
            },
            type_: INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        let allocation_before = allocation_failure.clone();
        heap.fail_after_allocations(0);
        assert_eq!(
            inchi_ios_eprint(
                &mut heap,
                Some(&mut allocation_failure),
                literal.as_const(),
                &SourceVaList::default(),
            ),
            Ok(-1)
        );
        assert_eq!(allocation_failure, allocation_before);
        assert!(heap.slice(allocation_old.as_const()).is_ok());

        let file = heap.allocate(vec![SourceFile::default()]).unwrap();
        let mut file_ios = INCHI_IOSTREAM {
            f: file,
            type_: INCHI_IOS_TYPE_FILE as i32,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(
            inchi_ios_eprint(
                &mut heap,
                Some(&mut file_ios),
                format.as_const(),
                &arguments,
            ),
            Ok(8)
        );
        assert_eq!(heap.slice(file.as_const()).unwrap()[0].bytes, b"value=17");

        file_ios.f = SourceMutPointer::null();
        assert_eq!(
            inchi_ios_eprint(
                &mut heap,
                Some(&mut file_ios),
                format.as_const(),
                &arguments,
            ),
            Ok(0)
        );
        file_ios.type_ = 99;
        assert_eq!(
            inchi_ios_eprint(
                &mut heap,
                Some(&mut file_ios),
                format.as_const(),
                &arguments,
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichi_io__inchi_sgets__line_1322() {
        let mut heap = SourceHeap::default();
        let input = source_format(&mut heap, "ab\ncd");
        let output = heap.allocate(vec![99_i8; 8]).unwrap();
        let mut ios = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: input,
                nAllocatedLength: 6,
                nUsedLength: 5,
                nPtr: 0,
            },
            type_: INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(
            inchi_sgets(&mut heap, output, 8, Some(&mut ios)),
            Ok(output)
        );
        assert_eq!(
            &heap.slice(output.as_const()).unwrap()[..4],
            &[b'a' as i8, b'b' as i8, b'\n' as i8, 0]
        );
        assert_eq!(ios.s.nPtr, 3);
        assert_eq!(
            inchi_sgets(&mut heap, output, 3, Some(&mut ios)),
            Ok(output)
        );
        assert_eq!(
            &heap.slice(output.as_const()).unwrap()[..3],
            &[b'c' as i8, b'd' as i8, 0]
        );
        assert_eq!(ios.s.nPtr, 5);
        let before_eof = heap.slice(output.as_const()).unwrap().to_vec();
        assert!(
            inchi_sgets(&mut heap, output, 8, Some(&mut ios))
                .unwrap()
                .is_null()
        );
        assert_eq!(ios.s.nPtr, 5);
        assert_eq!(heap.slice(output.as_const()).unwrap()[0], 0);
        assert_ne!(heap.slice(output.as_const()).unwrap(), before_eof);

        heap.slice_mut(output).unwrap()[0] = 77;
        assert!(
            inchi_sgets(
                &mut heap,
                output,
                1,
                Some(&mut INCHI_IOSTREAM {
                    s: INCHI_IOS_STRING {
                        pStr: input,
                        nPtr: 0,
                        ..INCHI_IOS_STRING::default()
                    },
                    ..INCHI_IOSTREAM::default()
                })
            )
            .unwrap()
            .is_null()
        );
        assert_eq!(heap.slice(output.as_const()).unwrap()[0], 0);

        heap.slice_mut(output).unwrap()[0] = 88;
        assert!(
            inchi_sgets(&mut heap, output, 0, Some(&mut ios))
                .unwrap()
                .is_null()
        );
        assert_eq!(heap.slice(output.as_const()).unwrap()[0], 88);
        let mut null_input = INCHI_IOSTREAM::default();
        assert!(
            inchi_sgets(&mut heap, output, 8, Some(&mut null_input))
                .unwrap()
                .is_null()
        );
        assert_eq!(heap.slice(output.as_const()).unwrap()[0], 88);
        assert_eq!(
            inchi_sgets(&mut heap, output, 8, None),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichi_io__inchi_fgetslf__line_996() {
        let mut heap = SourceHeap::default();
        let line = heap.allocate(vec![99_i8; 8]).unwrap();

        let file = heap
            .allocate(vec![SourceFile {
                bytes: b"abcdefghi\nnext\r\r\nlast\n".to_vec(),
                ..SourceFile::default()
            }])
            .unwrap();
        let mut file_stream = INCHI_IOSTREAM {
            f: file,
            type_: INCHI_IOS_TYPE_FILE as i32,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(
            inchi_fgetsLf(&mut heap, line, 5, Some(&mut file_stream)),
            Ok(line)
        );
        assert_eq!(
            &heap.slice(line.as_const()).unwrap()[..5],
            &[b'a' as i8, b'b' as i8, b'c' as i8, b'd' as i8, 0]
        );
        assert_eq!(heap.slice(file.as_const()).unwrap()[0].position, 10);
        assert_eq!(
            inchi_fgetsLf(&mut heap, line, 8, Some(&mut file_stream)),
            Ok(line)
        );
        assert_eq!(
            &heap.slice(line.as_const()).unwrap()[..6],
            &[
                b'n' as i8,
                b'e' as i8,
                b'x' as i8,
                b't' as i8,
                b'\n' as i8,
                0
            ]
        );
        assert_eq!(heap.slice(file.as_const()).unwrap()[0].position, 17);
        assert_eq!(
            inchi_fgetsLf(&mut heap, line, 8, Some(&mut file_stream)),
            Ok(line)
        );
        assert_eq!(
            &heap.slice(line.as_const()).unwrap()[..6],
            &[
                b'l' as i8,
                b'a' as i8,
                b's' as i8,
                b't' as i8,
                b'\n' as i8,
                0
            ]
        );
        assert!(
            inchi_fgetsLf(&mut heap, line, 8, Some(&mut file_stream))
                .unwrap()
                .is_null()
        );
        assert_eq!(heap.slice(line.as_const()).unwrap()[..8], [0; 8]);

        let string = source_format(&mut heap, "123456789\nxy\r\r\nz\n");
        let mut string_stream = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: string,
                nAllocatedLength: 19,
                nUsedLength: 18,
                nPtr: 0,
            },
            type_: INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(
            inchi_fgetsLf(&mut heap, line, 5, Some(&mut string_stream)),
            Ok(line)
        );
        assert_eq!(
            &heap.slice(line.as_const()).unwrap()[..5],
            &[b'1' as i8, b'2' as i8, b'3' as i8, b'4' as i8, 0]
        );
        assert_eq!(string_stream.s.nPtr, 10);
        assert_eq!(
            inchi_fgetsLf(&mut heap, line, 8, Some(&mut string_stream)),
            Ok(line)
        );
        assert_eq!(
            &heap.slice(line.as_const()).unwrap()[..4],
            &[b'x' as i8, b'y' as i8, b'\n' as i8, 0]
        );
        assert_eq!(string_stream.s.nPtr, 15);

        heap.slice_mut(line).unwrap().fill(77);
        let mut unknown = INCHI_IOSTREAM {
            type_: 99,
            ..INCHI_IOSTREAM::default()
        };
        assert!(
            inchi_fgetsLf(&mut heap, line, 8, Some(&mut unknown))
                .unwrap()
                .is_null()
        );
        assert_eq!(heap.slice(line.as_const()).unwrap(), &[77; 8]);
        assert_eq!(
            inchi_fgetsLf(&mut heap, line, 8, None),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichi_io__inchi_ios_close__line_229() {
        let mut heap = SourceHeap::default();
        assert_eq!(inchi_ios_close(&mut heap, None), Ok(()));

        let string = allocate_source_fixture(&mut heap, vec![b'x' as i8, 0]);
        let file = allocate_source_fixture(
            &mut heap,
            vec![SourceFile {
                bytes: b"file".to_vec(),
                ..SourceFile::default()
            }],
        );
        let mut ios = INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            f: file,
            s: INCHI_IOS_STRING {
                pStr: string,
                nAllocatedLength: 2,
                nUsedLength: 1,
                nPtr: 1,
            },
        };
        assert_eq!(inchi_ios_close(&mut heap, Some(&mut ios)), Ok(()));
        assert!(ios.s.pStr.is_null());
        assert_eq!(
            (ios.s.nAllocatedLength, ios.s.nUsedLength, ios.s.nPtr),
            (0, 0, 0)
        );
        assert_eq!(ios.f, file);
        assert_eq!(
            heap.slice(string.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(file.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );

        let standard = allocate_source_fixture(
            &mut heap,
            vec![SourceFile {
                is_standard_stream: true,
                ..SourceFile::default()
            }],
        );
        let mut standard_ios = INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_FILE as i32,
            f: standard,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(inchi_ios_close(&mut heap, Some(&mut standard_ios)), Ok(()));
        assert_eq!(heap.slice(standard.as_const()).unwrap().len(), 1);

        let mut empty = INCHI_IOSTREAM::default();
        assert_eq!(inchi_ios_close(&mut heap, Some(&mut empty)), Ok(()));
    }

    #[test]
    fn source_port__ichi_io__inchi_ios_reset__line_255() {
        let mut heap = SourceHeap::default();
        let string = allocate_source_fixture(&mut heap, vec![b'k' as i8, 0]);
        let file = allocate_source_fixture(&mut heap, vec![SourceFile::default()]);
        let mut ios = INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            f: file,
            s: INCHI_IOS_STRING {
                pStr: string,
                nAllocatedLength: 9,
                nUsedLength: 7,
                nPtr: 5,
            },
        };
        assert_eq!(inchi_ios_reset(&mut heap, &mut ios), Ok(()));
        assert!(ios.s.pStr.is_null());
        assert_eq!(
            (ios.s.nAllocatedLength, ios.s.nUsedLength, ios.s.nPtr),
            (0, 0, 0)
        );
        assert_eq!(heap.slice(string.as_const()).unwrap(), &[b'k' as i8, 0]);
        assert_eq!(
            heap.slice(file.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(ios.f, file);

        let standard = allocate_source_fixture(
            &mut heap,
            vec![SourceFile {
                is_standard_stream: true,
                ..SourceFile::default()
            }],
        );
        let mut standard_ios = INCHI_IOSTREAM {
            f: standard,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(inchi_ios_reset(&mut heap, &mut standard_ios), Ok(()));
        assert_eq!(heap.slice(standard.as_const()).unwrap().len(), 1);

        let mut empty = INCHI_IOSTREAM::default();
        assert_eq!(inchi_ios_reset(&mut heap, &mut empty), Ok(()));
    }

    #[test]
    fn source_port__ichi_io__inchi_ios_init__line_85() {
        let mut heap = SourceHeap::default();
        let old_string = allocate_source_fixture(&mut heap, vec![b'o' as i8, 0]);
        let old_file = allocate_source_fixture(
            &mut heap,
            vec![SourceFile {
                bytes: b"old".to_vec(),
                ..SourceFile::default()
            }],
        );
        let new_file = allocate_source_fixture(
            &mut heap,
            vec![SourceFile {
                bytes: b"new".to_vec(),
                ..SourceFile::default()
            }],
        );
        let mut ios = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: old_string,
                nAllocatedLength: 2,
                nUsedLength: 1,
                nPtr: 1,
            },
            f: old_file,
            type_: 97,
        };

        assert_eq!(
            inchi_ios_init(Some(&mut ios), INCHI_IOS_TYPE_FILE as i32, new_file),
            Ok(())
        );
        assert_eq!(ios.s, INCHI_IOS_STRING::default());
        assert_eq!(ios.type_, INCHI_IOS_TYPE_FILE as i32);
        assert_eq!(ios.f, new_file);
        assert_eq!(heap.slice(old_string.as_const()).unwrap(), &[b'o' as i8, 0]);
        assert_eq!(heap.slice(old_file.as_const()).unwrap()[0].bytes, b"old");

        for io_type in [INCHI_IOS_TYPE_STRING as i32, 0, -1, i32::MIN, i32::MAX] {
            ios.s = INCHI_IOS_STRING {
                pStr: old_string,
                nAllocatedLength: 11,
                nUsedLength: 7,
                nPtr: 3,
            };
            ios.type_ = INCHI_IOS_TYPE_FILE as i32;
            assert_eq!(
                inchi_ios_init(Some(&mut ios), io_type, SourceMutPointer::null()),
                Ok(())
            );
            assert_eq!(
                ios,
                INCHI_IOSTREAM {
                    type_: INCHI_IOS_TYPE_STRING as i32,
                    ..INCHI_IOSTREAM::default()
                }
            );
        }

        assert_eq!(
            inchi_ios_init(None, INCHI_IOS_TYPE_FILE as i32, new_file),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(heap.slice(new_file.as_const()).unwrap()[0].bytes, b"new");

        assert_eq!(inchi_free(&mut heap, old_string), Ok(()));
        assert_eq!(inchi_free(&mut heap, old_file), Ok(()));
        assert_eq!(inchi_free(&mut heap, new_file), Ok(()));
    }

    #[test]
    fn source_port__ichi_io__inchi_ios_getstab1__line_451() {
        let mut heap = SourceHeap::default();
        let input = allocate_source_fixture(
            &mut heap,
            b"\t  value \t".iter().map(|byte| *byte as i8).collect(),
        );
        let mut ios = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: input,
                nAllocatedLength: 10,
                nUsedLength: 10,
                nPtr: 0,
            },
            f: SourceMutPointer::null(),
            type_: INCHI_IOS_TYPE_STRING as i32,
        };
        let line = allocate_source_fixture(&mut heap, vec![99_i8; 12]);
        let mut too_long = -1;

        assert_eq!(
            inchi_ios_getsTab1(&mut heap, line, 12, &mut ios, &mut too_long),
            Ok(0)
        );
        assert_eq!(heap.slice(line.as_const()).unwrap()[0], 0);
        assert_eq!(too_long, 0);
        assert_eq!(ios.s.nPtr, 1);

        assert_eq!(
            inchi_ios_getsTab1(&mut heap, line, 12, &mut ios, &mut too_long),
            Ok(5)
        );
        assert_eq!(
            &heap.slice(line.as_const()).unwrap()[..6],
            &[
                b'v' as i8, b'a' as i8, b'l' as i8, b'u' as i8, b'e' as i8, 0
            ]
        );
        assert_eq!(too_long, 0);
        assert_eq!(ios.s.nPtr, 10);

        assert_eq!(inchi_free(&mut heap, input), Ok(()));
        assert_eq!(inchi_free(&mut heap, line), Ok(()));
    }

    #[test]
    fn source_port__ichi_io__inchi_ios_getstab__line_420() {
        let mut heap = SourceHeap::default();
        let input = allocate_source_fixture(
            &mut heap,
            b"\t  value \t".iter().map(|byte| *byte as i8).collect(),
        );
        let mut ios = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: input,
                nAllocatedLength: 10,
                nUsedLength: 10,
                nPtr: 0,
            },
            f: SourceMutPointer::null(),
            type_: INCHI_IOS_TYPE_STRING as i32,
        };
        let line = allocate_source_fixture(&mut heap, vec![99_i8; 12]);
        let mut too_long = -1;

        assert_eq!(
            inchi_ios_getsTab(&mut heap, line, 12, &mut ios, &mut too_long),
            Ok(5)
        );
        assert_eq!(
            &heap.slice(line.as_const()).unwrap()[..6],
            &[
                b'v' as i8, b'a' as i8, b'l' as i8, b'u' as i8, b'e' as i8, 0
            ]
        );
        assert_eq!(too_long, 0);
        assert_eq!(ios.s.nPtr, 10);

        assert_eq!(inchi_free(&mut heap, input), Ok(()));
        assert_eq!(inchi_free(&mut heap, line), Ok(()));
    }

    #[test]
    fn source_port__ichi_io__inchi_ios_str_getstab__line_354() {
        let mut heap = SourceHeap::default();
        let input = allocate_source_fixture(
            &mut heap,
            b"ab\tcd\n".iter().map(|byte| *byte as i8).collect(),
        );
        let mut ios = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: input,
                nAllocatedLength: 6,
                nUsedLength: 6,
                nPtr: 0,
            },
            f: SourceMutPointer::null(),
            type_: INCHI_IOS_TYPE_STRING as i32,
        };
        let line = allocate_source_fixture(&mut heap, vec![99_i8; 8]);

        assert_eq!(
            inchi_ios_str_getsTab(&mut heap, line, 8, &mut ios),
            Ok(line)
        );
        assert_eq!(
            &heap.slice(line.as_const()).unwrap()[..4],
            &[b'a' as i8, b'b' as i8, b'\n' as i8, 0]
        );
        assert_eq!(ios.s.nPtr, 3);

        assert_eq!(
            inchi_ios_str_getsTab(&mut heap, line, 8, &mut ios),
            Ok(line)
        );
        assert_eq!(
            &heap.slice(line.as_const()).unwrap()[..4],
            &[b'c' as i8, b'd' as i8, b'\n' as i8, 0]
        );
        assert_eq!(ios.s.nPtr, 6);
        let before_eof = heap.slice(line.as_const()).unwrap().to_vec();
        assert_eq!(
            inchi_ios_str_getsTab(&mut heap, line, 8, &mut ios),
            Ok(SourceMutPointer::null())
        );
        assert_eq!(heap.slice(line.as_const()).unwrap(), before_eof);

        assert_eq!(inchi_free(&mut heap, input), Ok(()));
        assert_eq!(inchi_free(&mut heap, line), Ok(()));
    }

    #[test]
    fn source_port__ichi_io__inchi_ios_gets__line_386() {
        let mut heap = SourceHeap::default();
        let input = allocate_source_fixture(
            &mut heap,
            b" \n  abc  \n".iter().map(|byte| *byte as i8).collect(),
        );
        let mut ios = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: input,
                nAllocatedLength: 10,
                nUsedLength: 10,
                nPtr: 0,
            },
            f: SourceMutPointer::null(),
            type_: INCHI_IOS_TYPE_STRING as i32,
        };
        let line = allocate_source_fixture(&mut heap, vec![99_i8; 10]);
        let mut too_long = -1;
        assert_eq!(
            inchi_ios_gets(&mut heap, line, 10, &mut ios, &mut too_long),
            Ok(3)
        );
        assert_eq!(
            &heap.slice(line.as_const()).unwrap()[..4],
            &[b'a' as i8, b'b' as i8, b'c' as i8, 0]
        );
        assert_eq!(too_long, 0);
        assert_eq!(ios.s.nPtr, 10);
        assert_eq!(inchi_free(&mut heap, input), Ok(()));
        assert_eq!(inchi_free(&mut heap, line), Ok(()));

        let input = allocate_source_fixture(
            &mut heap,
            b"abcdef".iter().map(|byte| *byte as i8).collect(),
        );
        let mut ios = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: input,
                nAllocatedLength: 6,
                nUsedLength: 6,
                nPtr: 0,
            },
            f: SourceMutPointer::null(),
            type_: INCHI_IOS_TYPE_STRING as i32,
        };
        let line = allocate_source_fixture(&mut heap, vec![0_i8; 5]);
        let mut too_long = 0;
        assert_eq!(
            inchi_ios_gets(&mut heap, line, 5, &mut ios, &mut too_long),
            Ok(3)
        );
        assert_eq!(
            &heap.slice(line.as_const()).unwrap()[..4],
            &[b'a' as i8, b'b' as i8, b'c' as i8, 0]
        );
        assert_eq!(too_long, 1);
        assert_eq!(ios.s.nPtr, 3);
        assert_eq!(inchi_free(&mut heap, input), Ok(()));
        assert_eq!(inchi_free(&mut heap, line), Ok(()));

        let empty = allocate_source_fixture(&mut heap, Vec::<i8>::new());
        let mut empty_ios = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: empty,
                nAllocatedLength: 0,
                nUsedLength: 0,
                nPtr: 0,
            },
            f: SourceMutPointer::null(),
            type_: INCHI_IOS_TYPE_STRING as i32,
        };
        let line = allocate_source_fixture(&mut heap, vec![99_i8; 5]);
        let mut too_long = 7;
        assert_eq!(
            inchi_ios_gets(&mut heap, line, 5, &mut empty_ios, &mut too_long),
            Ok(-1)
        );
        assert_eq!(too_long, 0);
        assert_eq!(heap.slice(line.as_const()).unwrap(), &[99_i8; 5]);
        assert_eq!(inchi_free(&mut heap, empty), Ok(()));
        assert_eq!(inchi_free(&mut heap, line), Ok(()));
    }

    #[test]
    fn source_port__ichi_io__inchi_ios_str_gets__line_324() {
        let mut heap = SourceHeap::default();
        let input = allocate_source_fixture(
            &mut heap,
            vec![b'a' as i8, b'b' as i8, b'\n' as i8, b'z' as i8],
        );
        let mut ios = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: input,
                nAllocatedLength: 4,
                nUsedLength: 4,
                nPtr: 0,
            },
            f: SourceMutPointer::null(),
            type_: INCHI_IOS_TYPE_STRING as i32,
        };
        let line = allocate_source_fixture(&mut heap, vec![99_i8; 8]);
        assert_eq!(inchi_ios_str_gets(&mut heap, line, 8, &mut ios), Ok(line));
        assert_eq!(
            heap.slice(line.as_const()).unwrap(),
            &[b'a' as i8, b'b' as i8, b'\n' as i8, 0, 99, 99, 99, 99]
        );
        assert_eq!(ios.s.nPtr, 3);

        assert_eq!(inchi_ios_str_gets(&mut heap, line, 2, &mut ios), Ok(line));
        assert_eq!(heap.slice(line.as_const()).unwrap()[..2], [b'z' as i8, 0]);
        assert_eq!(ios.s.nPtr, 4);
        let before_eof = heap.slice(line.as_const()).unwrap().to_vec();
        assert_eq!(
            inchi_ios_str_gets(&mut heap, line, 8, &mut ios),
            Ok(SourceMutPointer::null())
        );
        assert_eq!(heap.slice(line.as_const()).unwrap(), before_eof);
        assert_eq!(inchi_free(&mut heap, input), Ok(()));
        assert_eq!(inchi_free(&mut heap, line), Ok(()));

        let one_byte_input = allocate_source_fixture(&mut heap, vec![b'Q' as i8]);
        let mut one_byte_ios = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: one_byte_input,
                nAllocatedLength: 1,
                nUsedLength: 1,
                nPtr: 0,
            },
            f: SourceMutPointer::null(),
            type_: INCHI_IOS_TYPE_STRING as i32,
        };
        let one_byte_line = allocate_source_fixture(&mut heap, vec![99_i8]);
        assert_eq!(
            inchi_ios_str_gets(&mut heap, one_byte_line, 1, &mut one_byte_ios),
            Ok(one_byte_line)
        );
        assert_eq!(heap.slice(one_byte_line.as_const()).unwrap(), &[0]);
        assert_eq!(one_byte_ios.s.nPtr, 0);
        heap.slice_mut(one_byte_line).unwrap()[0] = 99;
        assert_eq!(
            inchi_ios_str_gets(&mut heap, one_byte_line, 0, &mut one_byte_ios),
            Ok(SourceMutPointer::null())
        );
        assert_eq!(heap.slice(one_byte_line.as_const()).unwrap(), &[99]);
        assert_eq!(one_byte_ios.s.nPtr, 0);
        assert_eq!(
            inchi_ios_str_gets(&mut heap, one_byte_line, i32::MIN, &mut one_byte_ios),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(inchi_free(&mut heap, one_byte_input), Ok(()));
        assert_eq!(inchi_free(&mut heap, one_byte_line), Ok(()));

        let file = allocate_source_fixture(
            &mut heap,
            vec![SourceFile {
                bytes: vec![255],
                position: 0,
                error: false,
                error_indicator: false,
                eof: false,
                is_standard_stream: false,
            }],
        );
        let mut file_ios = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING::default(),
            f: file,
            type_: INCHI_IOS_TYPE_FILE as i32,
        };
        let file_line = allocate_source_fixture(&mut heap, vec![0_i8; 2]);
        assert_eq!(
            inchi_ios_str_gets(&mut heap, file_line, 2, &mut file_ios),
            Ok(file_line)
        );
        assert_eq!(heap.slice(file_line.as_const()).unwrap(), &[-1, 0]);
        assert_eq!(inchi_free(&mut heap, file), Ok(()));
        assert_eq!(inchi_free(&mut heap, file_line), Ok(()));
    }

    #[test]
    fn source_port__ichi_io__inchi_ios_str_getc__line_294() {
        let mut heap = SourceHeap::default();
        let string = allocate_source_fixture(&mut heap, vec![b'A' as i8, 0, -1]);
        let mut string_ios = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: string,
                nAllocatedLength: 3,
                nUsedLength: 3,
                nPtr: 0,
            },
            f: SourceMutPointer::null(),
            type_: INCHI_IOS_TYPE_STRING as i32,
        };
        assert_eq!(inchi_ios_str_getc(&mut heap, &mut string_ios), Ok(65));
        assert_eq!(string_ios.s.nPtr, 1);
        assert_eq!(inchi_ios_str_getc(&mut heap, &mut string_ios), Ok(0));
        assert_eq!(string_ios.s.nPtr, 2);
        assert_eq!(
            inchi_ios_str_getc(&mut heap, &mut string_ios),
            Ok(SOURCE_EOF)
        );
        assert_eq!(string_ios.s.nPtr, 3);
        assert_eq!(
            inchi_ios_str_getc(&mut heap, &mut string_ios),
            Ok(SOURCE_EOF)
        );
        assert_eq!(string_ios.s.nPtr, 3);
        assert_eq!(inchi_free(&mut heap, string), Ok(()));

        let file = allocate_source_fixture(
            &mut heap,
            vec![SourceFile {
                bytes: vec![0, 255],
                position: 0,
                error: false,
                error_indicator: false,
                eof: false,
                is_standard_stream: false,
            }],
        );
        let mut file_ios = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING::default(),
            f: file,
            type_: INCHI_IOS_TYPE_FILE as i32,
        };
        assert_eq!(inchi_ios_str_getc(&mut heap, &mut file_ios), Ok(0));
        assert_eq!(inchi_ios_str_getc(&mut heap, &mut file_ios), Ok(255));
        assert_eq!(inchi_ios_str_getc(&mut heap, &mut file_ios), Ok(SOURCE_EOF));
        let file_state = &heap.slice(file.as_const()).unwrap()[0];
        assert_eq!(file_state.position, 2);
        assert!(file_state.eof);
        assert!(!file_state.error);
        assert_eq!(inchi_free(&mut heap, file), Ok(()));

        let mut invalid_ios = INCHI_IOSTREAM::default();
        assert_eq!(
            inchi_ios_str_getc(&mut heap, &mut invalid_ios),
            Ok(SOURCE_EOF)
        );
    }
}
