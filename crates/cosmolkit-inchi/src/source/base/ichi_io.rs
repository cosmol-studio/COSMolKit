use crate::source::base::util::{inchi_calloc, inchi_free, lrtrim};
use crate::source_types::{
    FILE, INCHI_IOS_STRING, INCHI_IOS_TYPE_FILE, INCHI_IOS_TYPE_STRING, INCHI_IOSTREAM,
    INCHI_STRBUF_INITIAL_SIZE, INCHI_STRBUF_SIZE_INCREMENT, SourceConstPointer,
    SourceFormatArgument, SourceHeap, SourceHeapError, SourceMutPointer, SourceVaList,
    local_ichi_io::INCHI_ADD_STR_LEN,
};

const SOURCE_EOF: i32 = -1;

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
                let old = heap
                    .slice(buf.pStr.as_const())?
                    .get(..used)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                heap.slice_mut(new_string)?
                    .get_mut(..used)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .copy_from_slice(&old);
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

    let Some(buf) = buf else {
        return Ok(-1);
    };
    let mut length_arguments = arguments.clone();
    let maximum_length = GetMaxPrintfLength(heap, format, &mut length_arguments)?;
    if maximum_length < 0 {
        return Ok(0);
    }
    if inchi_strbuf_update(heap, Some(&mut *buf), maximum_length)? < 0 {
        return Err(SourceHeapError::AllocationFailed);
    }
    let mut render_arguments = arguments.clone();
    let rendered = source_vformat(heap, format, &mut render_arguments)?;
    let rendered_length =
        i32::try_from(rendered.len()).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let used = usize::try_from(buf.nUsedLength).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let end = used
        .checked_add(rendered.len())
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    let destination = heap.slice_mut(buf.pStr)?;
    destination
        .get_mut(used..end)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .copy_from_slice(&rendered.iter().map(|byte| *byte as i8).collect::<Vec<_>>());
    *destination
        .get_mut(end)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
    buf.nUsedLength = buf.nUsedLength.wrapping_add(rendered_length);
    Ok(rendered_length)
}

fn next_format_argument(
    arguments: &mut SourceVaList,
) -> Result<SourceFormatArgument, SourceHeapError> {
    let index =
        usize::try_from(arguments.position).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let argument = arguments
        .arguments
        .get(index)
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    arguments.position = arguments
        .position
        .checked_add(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    Ok(argument)
}

fn next_format_int(arguments: &mut SourceVaList) -> Result<i32, SourceHeapError> {
    match next_format_argument(arguments)? {
        SourceFormatArgument::Signed(value) => {
            i32::try_from(value).map_err(|_| SourceHeapError::SourceIntegerOverflow)
        }
        SourceFormatArgument::Unsigned(value) => Ok(value as u32 as i32),
        SourceFormatArgument::Byte(value) => Ok(i32::from(value)),
        _ => Err(SourceHeapError::AllocationTypeMismatch),
    }
}

fn apply_printf_width(mut value: String, width: i32, left: bool, zero: bool) -> String {
    let width = width.max(0) as usize;
    if value.len() >= width {
        return value;
    }
    let count = width - value.len();
    if left {
        value.extend(std::iter::repeat_n(' ', count));
        return value;
    }
    if zero {
        let prefix_length =
            if value.starts_with('+') || value.starts_with('-') || value.starts_with(' ') {
                1
            } else if value.starts_with("0x") || value.starts_with("0X") {
                2
            } else {
                0
            };
        value.insert_str(prefix_length, &"0".repeat(count));
        return value;
    }
    format!("{}{}", " ".repeat(count), value)
}

fn c_printf_exponential(value: f64, precision: usize, upper: bool) -> String {
    let rust = format!("{value:.precision$e}");
    let (mantissa, exponent) = rust
        .split_once('e')
        .expect("scientific format has exponent");
    let exponent = exponent
        .parse::<i32>()
        .expect("scientific exponent is decimal");
    let marker = if upper { 'E' } else { 'e' };
    format!("{mantissa}{marker}{exponent:+03}")
}

fn trim_general_float(mut value: String, alternate: bool) -> String {
    if alternate {
        return value;
    }
    let exponent_index = value.find(['e', 'E']).unwrap_or(value.len());
    let exponent = value.split_off(exponent_index);
    while value.ends_with('0') && value.contains('.') {
        value.pop();
    }
    if value.ends_with('.') {
        value.pop();
    }
    value.push_str(&exponent);
    value
}

fn c_printf_float(value: f64, specifier: u8, precision: Option<usize>, alternate: bool) -> String {
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
            let mut significant = precision.unwrap_or(6);
            if significant == 0 {
                significant = 1;
            }
            let exponent = if value == 0.0 {
                0
            } else {
                value.abs().log10().floor() as i32
            };
            if exponent < -4 || exponent >= significant as i32 {
                trim_general_float(
                    c_printf_exponential(value, significant.saturating_sub(1), upper),
                    alternate,
                )
            } else {
                let decimals = (significant as i32 - exponent - 1).max(0) as usize;
                trim_general_float(format!("{value:.decimals$}"), alternate)
            }
        }
        _ => unreachable!(),
    }
}

fn source_vformat(
    heap: &mut SourceHeap,
    format: SourceConstPointer<i8>,
    arguments: &mut SourceVaList,
) -> Result<Vec<u8>, SourceHeapError> {
    let format = heap.slice(format)?;
    let format_length = format
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let format = format[..format_length].to_vec();
    let mut output = Vec::new();
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

        let mut rendered = match specifier {
            b'c' => {
                let value = match next_format_argument(arguments)? {
                    SourceFormatArgument::Signed(value) => value as u8,
                    SourceFormatArgument::Byte(value) => value as u8,
                    _ => return Err(SourceHeapError::AllocationTypeMismatch),
                };
                String::from(value as char)
            }
            b's' => {
                let pointer = match next_format_argument(arguments)? {
                    SourceFormatArgument::Bytes(pointer) => pointer,
                    _ => return Err(SourceHeapError::AllocationTypeMismatch),
                };
                let bytes = heap.slice(pointer)?;
                let length = bytes
                    .iter()
                    .position(|byte| *byte == 0)
                    .ok_or(SourceHeapError::MissingNulTerminator)?;
                let length = precision.map_or(length, |precision| precision.min(length));
                String::from_utf8(bytes[..length].iter().map(|byte| *byte as u8).collect())
                    .map_err(|_| SourceHeapError::InvalidSourceTextEncoding)?
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
                let negative = value < 0;
                let mut digits = value.unsigned_abs().to_string();
                if precision == Some(0) && value == 0 {
                    digits.clear();
                }
                if let Some(precision) = precision {
                    if digits.len() < precision {
                        digits.insert_str(0, &"0".repeat(precision - digits.len()));
                    }
                }
                if negative {
                    digits.insert(0, '-');
                } else if plus {
                    digits.insert(0, '+');
                } else if space {
                    digits.insert(0, ' ');
                }
                digits
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
                let mut digits = match specifier {
                    b'u' => value.to_string(),
                    b'o' => format!("{value:o}"),
                    b'x' => format!("{value:x}"),
                    b'X' => format!("{value:X}"),
                    _ => unreachable!(),
                };
                if precision == Some(0) && value == 0 {
                    digits.clear();
                }
                if let Some(precision) = precision {
                    if digits.len() < precision {
                        digits.insert_str(0, &"0".repeat(precision - digits.len()));
                    }
                }
                if alternate && specifier == b'o' && !digits.starts_with('0') {
                    digits.insert(0, '0');
                } else if alternate && value != 0 && specifier == b'x' {
                    digits.insert_str(0, "0x");
                } else if alternate && value != 0 && specifier == b'X' {
                    digits.insert_str(0, "0X");
                }
                digits
            }
            b'e' | b'E' | b'f' | b'F' | b'g' | b'G' => {
                let value = match next_format_argument(arguments)? {
                    SourceFormatArgument::Float(value) => value,
                    _ => return Err(SourceHeapError::AllocationTypeMismatch),
                };
                let negative = value.is_sign_negative();
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
            b'p' => match next_format_argument(arguments)? {
                SourceFormatArgument::Pointer(pointer) if pointer.is_null() => "(nil)".to_owned(),
                SourceFormatArgument::Pointer(_) => {
                    return Err(SourceHeapError::UnsupportedSourceBehavior);
                }
                _ => return Err(SourceHeapError::AllocationTypeMismatch),
            },
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
        if precision.is_some() && matches!(specifier, b'd' | b'i' | b'u' | b'o' | b'x' | b'X') {
            zero = false;
        }
        rendered = apply_printf_width(rendered, width, left, zero && !left);
        output.extend_from_slice(rendered.as_bytes());
    }
    Ok(output)
}

pub(crate) fn inchi_vfprintf(
    heap: &mut SourceHeap,
    file: SourceMutPointer<FILE>,
    format: SourceConstPointer<i8>,
    arguments: &mut SourceVaList,
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

    if format.is_null() {
        return Ok(0);
    }
    let format_bytes = heap.slice(format)?;
    if format_bytes.first().copied() == Some(0) {
        return Ok(0);
    }
    let output = source_vformat(heap, format, arguments)?;
    let file = heap
        .slice_mut(file)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if file.error {
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

    let destination = file.unwrap_or(stdout);
    let mut arguments = arguments.clone();
    inchi_vfprintf(heap, destination, format, &mut arguments)
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

    let Some(ios) = ios else {
        return Ok(-1);
    };
    if ios.type_ == INCHI_IOS_TYPE_STRING as i32 {
        let mut estimate_arguments = arguments.clone();
        let maximum_length = GetMaxPrintfLength(heap, format, &mut estimate_arguments)?;
        if maximum_length < 0 {
            return Ok(-1);
        }
        let additional_length = (INCHI_ADD_STR_LEN as i32).max(maximum_length);
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
                    let bytes = heap
                        .slice(ios.s.pStr.as_const())?
                        .get(..used)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .to_vec();
                    heap.slice_mut(new_string)?[..used].copy_from_slice(&bytes);
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
        let mut output_arguments = arguments.clone();
        let output = source_vformat(heap, format, &mut output_arguments)?;
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
        let mut arguments = arguments.clone();
        return inchi_vfprintf(heap, destination, format, &mut arguments);
    }
    Ok(0)
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

    let Some(ios) = ios else {
        return Ok(-1);
    };
    if ios.type_ == INCHI_IOS_TYPE_STRING as i32 {
        let mut estimate_arguments = arguments.clone();
        let maximum_length = GetMaxPrintfLength(heap, format, &mut estimate_arguments)?;
        if maximum_length < 0 {
            return Ok(-1);
        }
        let additional_length = (INCHI_ADD_STR_LEN as i32).max(maximum_length);
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
                    let bytes = heap.slice(ios.s.pStr.as_const())?[..used].to_vec();
                    heap.slice_mut(new_string)?[..used].copy_from_slice(&bytes);
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
        let mut output_arguments = arguments.clone();
        let output = source_vformat(heap, format, &mut output_arguments)?;
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
        let file = (!ios.f.is_null()).then_some(ios.f);
        inchi_print_nodisplay(heap, file, stdout, format, arguments)?;
    }
    Ok(0)
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

    let Some(ios) = ios else {
        return Ok(-1);
    };
    if ios.type_ == INCHI_IOS_TYPE_STRING as i32 {
        let mut estimate_arguments = arguments.clone();
        let maximum_length = GetMaxPrintfLength(heap, format, &mut estimate_arguments)?;
        if maximum_length < 0 {
            return Ok(-1);
        }
        if ios.s.nAllocatedLength - ios.s.nUsedLength <= maximum_length {
            let additional_length = (INCHI_ADD_STR_LEN as i32).max(maximum_length);
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
                    let bytes = old[..used].to_vec();
                    heap.slice_mut(new_string)?[..used].copy_from_slice(&bytes);
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

        let mut output_arguments = arguments.clone();
        let output = source_vformat(heap, format, &mut output_arguments)?;
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
        destination[used..used + output.len()]
            .copy_from_slice(&output.iter().map(|byte| *byte as i8).collect::<Vec<_>>());
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
        let mut output_arguments = arguments.clone();
        return inchi_vfprintf(heap, ios.f, format, &mut output_arguments);
    }
    Ok(0)
}

pub(crate) fn inchi_sgets(
    heap: &mut SourceHeap,
    output: SourceMutPointer<i8>,
    mut capacity: i32,
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

    let ios = ios.ok_or(SourceHeapError::NullPointer)?;
    if capacity <= 0 {
        return Ok(SourceMutPointer::null());
    }
    if ios.s.pStr.is_null() {
        return Ok(SourceMutPointer::null());
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
        *heap
            .slice_mut(output)?
            .get_mut(output_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = character;
        output_index += 1;
        if character == b'\n' as i8 {
            break;
        }
    }
    *heap
        .slice_mut(output)?
        .get_mut(output_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
    if character == 0 && output_index == 0 {
        Ok(SourceMutPointer::null())
    } else {
        Ok(output)
    }
}

fn source_file_fgets(
    heap: &mut SourceHeap,
    output: SourceMutPointer<i8>,
    capacity: i32,
    file: SourceMutPointer<FILE>,
) -> Result<SourceMutPointer<i8>, SourceHeapError> {
    if capacity <= 0 {
        return Ok(SourceMutPointer::null());
    }
    let (bytes, start, has_error) = {
        let file = heap
            .slice(file.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        (
            file.bytes.clone(),
            usize::try_from(file.position).map_err(|_| SourceHeapError::PointerOutOfBounds)?,
            file.error,
        )
    };
    if has_error {
        return Ok(SourceMutPointer::null());
    }
    if start >= bytes.len() {
        heap.slice_mut(file)?[0].eof = true;
        return Ok(SourceMutPointer::null());
    }
    let limit = usize::try_from(capacity - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mut end = start;
    while end - start < limit && end < bytes.len() {
        let byte = bytes[end];
        end += 1;
        if byte == b'\n' {
            break;
        }
    }
    let destination = heap.slice_mut(output)?;
    if destination.len() < end - start + 1 {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    for (destination, source) in destination.iter_mut().zip(&bytes[start..end]) {
        *destination = *source as i8;
    }
    destination[end - start] = 0;
    heap.slice_mut(file)?[0].position =
        u64::try_from(end).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    Ok(output)
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

    let input = input.ok_or(SourceHeapError::NullPointer)?;
    let mut result = SourceMutPointer::null();
    if input.type_ == INCHI_IOS_TYPE_FILE as i32 {
        if line_length > 0 {
            let length =
                usize::try_from(line_length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let destination = heap.slice_mut(line)?;
            if destination.len() < length {
                return Err(SourceHeapError::PointerOutOfBounds);
            }
            destination[..length].fill(0);
        }
        result = source_file_fgets(heap, line, line_length, input.f)?;
        if !result.is_null()
            && !heap
                .slice(result.as_const())?
                .iter()
                .take_while(|byte| **byte != 0)
                .any(|byte| *byte == b'\n' as i8)
        {
            let temporary = heap.allocate_model_storage(vec![0_i8; 64])?;
            loop {
                let chunk = source_file_fgets(heap, temporary, 64, input.f)?;
                if chunk.is_null()
                    || heap
                        .slice(chunk.as_const())?
                        .iter()
                        .take_while(|byte| **byte != 0)
                        .any(|byte| *byte == b'\n' as i8)
                {
                    break;
                }
            }
            heap.free(temporary)?;
        }
    } else if input.type_ == INCHI_IOS_TYPE_STRING as i32 {
        if line_length > 0 {
            let length =
                usize::try_from(line_length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let destination = heap.slice_mut(line)?;
            if destination.len() < length {
                return Err(SourceHeapError::PointerOutOfBounds);
            }
            destination[..length].fill(0);
        }
        result = inchi_sgets(heap, line, line_length, Some(input))?;
        if !result.is_null()
            && !heap
                .slice(result.as_const())?
                .iter()
                .take_while(|byte| **byte != 0)
                .any(|byte| *byte == b'\n' as i8)
        {
            let temporary = heap.allocate_model_storage(vec![0_i8; 64])?;
            loop {
                let chunk = inchi_sgets(heap, temporary, 64, Some(input))?;
                if chunk.is_null()
                    || heap
                        .slice(chunk.as_const())?
                        .iter()
                        .take_while(|byte| **byte != 0)
                        .any(|byte| *byte == b'\n' as i8)
                {
                    break;
                }
            }
            heap.free(temporary)?;
        }
    }

    if !result.is_null() {
        let bytes = heap.slice_mut(line)?;
        let length = bytes
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        if let Some(carriage_return) = bytes[..length].iter().position(|byte| *byte == b'\r' as i8)
        {
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
    let bytes = heap.slice(format)?;
    let format_length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let bytes = &bytes[..=format_length];
    let mut position = 0_usize;
    let mut maximum_length = 0_i32;
    while bytes[position] != 0 {
        if bytes[position] != b'%' as i8 {
            maximum_length = maximum_length
                .checked_add(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            position += 1;
            continue;
        }
        position += 1;
        if bytes[position] == b'%' as i8 {
            maximum_length = maximum_length
                .checked_add(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            position += 1;
            continue;
        }

        let mut item_length = 0_i32;
        let mut width = 0_i32;
        while bytes[position] != 0 {
            match bytes[position] as u8 {
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
            while (bytes[position] as u8).is_ascii_digit() {
                width = width
                    .checked_mul(10)
                    .and_then(|value| value.checked_add(i32::from(bytes[position] - b'0' as i8)))
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                position += 1;
            }
        }
        if width < 0 {
            return Ok(-1);
        }

        let mut precision = 0_i32;
        if bytes[position] == b'.' as i8 {
            position += 1;
            if bytes[position] == b'*' as i8 {
                precision = next_format_int(arguments)?;
                position += 1;
            } else {
                while (bytes[position] as u8).is_ascii_digit() {
                    precision = precision
                        .checked_mul(10)
                        .and_then(|value| {
                            value.checked_add(i32::from(bytes[position] - b'0' as i8))
                        })
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                    position += 1;
                }
            }
            if precision < 0 {
                return Ok(-1);
            }
        }

        let mut modifier = 0_i32;
        match bytes[position] as u8 {
            b'h' => {
                position += 1;
                if !matches!(
                    bytes[position] as u8,
                    b'd' | b'i' | b'o' | b'x' | b'X' | b'u'
                ) {
                    modifier = FORCE_ANSI;
                }
            }
            b'l' => {
                position += 1;
                if !matches!(
                    bytes[position] as u8,
                    b'd' | b'i' | b'o' | b'x' | b'X' | b'u' | b'f'
                ) {
                    return Ok(-1);
                }
            }
            b'F' | b'N' | b'L' => position += 1,
            _ => {}
        }

        let specifier = bytes[position] as u8;
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
                let pointer = match next_format_argument(arguments)? {
                    SourceFormatArgument::Bytes(pointer) => pointer,
                    _ => return Err(SourceHeapError::AllocationTypeMismatch),
                };
                let string = heap.slice(pointer)?;
                let length = string
                    .iter()
                    .position(|byte| *byte == 0)
                    .ok_or(SourceHeapError::MissingNulTerminator)?;
                item_length = i32::try_from(length)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                    .max(1);
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
                    next_format_int(arguments)?;
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
                    SourceFormatArgument::Pointer(_) => {}
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
    Ok(maximum_length)
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
        let index =
            usize::try_from(file.position).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let value = if let Some(value) = file.bytes.get(index).copied() {
            file.position += 1;
            i32::from(value)
        } else {
            file.eof = true;
            SOURCE_EOF
        };
        return Ok(if file.error { SOURCE_EOF } else { value });
    }

    Ok(SOURCE_EOF)
}

pub(crate) fn inchi_ios_str_gets(
    heap: &mut SourceHeap,
    line: SourceMutPointer<i8>,
    len: i32,
    ios: &mut INCHI_IOSTREAM,
) -> Result<SourceMutPointer<i8>, SourceHeapError> {
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

    let len = len
        .checked_sub(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    if len < 0 {
        return Ok(SourceMutPointer::null());
    }

    let mut length = 0_i32;
    let mut c = 0_i32;
    while length < len {
        c = inchi_ios_str_getc(heap, ios)?;
        if c == SOURCE_EOF {
            break;
        }
        let index = usize::try_from(length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        *heap
            .slice_mut(line)?
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = c as i8;
        length += 1;
        if c == i32::from(b'\n') {
            break;
        }
    }
    if length == 0 && c == SOURCE_EOF {
        return Ok(SourceMutPointer::null());
    }
    let index = usize::try_from(length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    *heap
        .slice_mut(line)?
        .get_mut(index)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
    Ok(line)
}

pub(crate) fn inchi_ios_str_getsTab(
    heap: &mut SourceHeap,
    line: SourceMutPointer<i8>,
    len: i32,
    ios: &mut INCHI_IOSTREAM,
) -> Result<SourceMutPointer<i8>, SourceHeapError> {
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

    let len = len
        .checked_sub(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    if len < 0 {
        return Ok(SourceMutPointer::null());
    }

    let mut length = 0_i32;
    let mut c = 0_i32;
    while length < len {
        c = inchi_ios_str_getc(heap, ios)?;
        if c == SOURCE_EOF {
            break;
        }
        if c == i32::from(b'\t') {
            c = i32::from(b'\n');
        }
        let index = usize::try_from(length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        *heap
            .slice_mut(line)?
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = c as i8;
        length += 1;
        if c == i32::from(b'\n') {
            break;
        }
    }
    if length == 0 && c == SOURCE_EOF {
        return Ok(SourceMutPointer::null());
    }
    let index = usize::try_from(length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    *heap
        .slice_mut(line)?
        .get_mut(index)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
    Ok(line)
}

pub(crate) fn inchi_ios_gets(
    heap: &mut SourceHeap,
    line: SourceMutPointer<i8>,
    len: i32,
    ios: &mut INCHI_IOSTREAM,
    too_long_line: &mut i32,
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

    let helper_len = len
        .checked_sub(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let content_limit = len
        .checked_sub(2)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    loop {
        let result = inchi_ios_str_gets(heap, line, helper_len, ios)?;
        if result.is_null() {
            *too_long_line = 0;
            return Ok(-1);
        }

        let terminator_index =
            usize::try_from(helper_len).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        *heap
            .slice_mut(line)?
            .get_mut(terminator_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
        let (has_newline, string_length) = {
            let bytes = heap.slice(line.as_const())?;
            let string_length = bytes
                .iter()
                .position(|byte| *byte == 0)
                .ok_or(SourceHeapError::MissingNulTerminator)?;
            (
                bytes[..string_length]
                    .iter()
                    .any(|byte| *byte == b'\n' as i8),
                i32::try_from(string_length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
            )
        };
        *too_long_line = i32::from(!has_newline && string_length == content_limit);

        let mut length = 0;
        lrtrim(heap, line, Some(&mut length))?;
        if length != 0 {
            return Ok(length);
        }
    }
}

pub(crate) fn inchi_ios_getsTab(
    heap: &mut SourceHeap,
    line: SourceMutPointer<i8>,
    len: i32,
    ios: &mut INCHI_IOSTREAM,
    too_long_line: &mut i32,
) -> Result<i32, SourceHeapError> {
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

    let helper_len = len
        .checked_sub(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let content_limit = len
        .checked_sub(2)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    loop {
        let result = inchi_ios_str_getsTab(heap, line, helper_len, ios)?;
        if result.is_null() {
            *too_long_line = 0;
            return Ok(-1);
        }

        let terminator_index =
            usize::try_from(helper_len).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        *heap
            .slice_mut(line)?
            .get_mut(terminator_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
        let (has_newline, string_length) = {
            let bytes = heap.slice(line.as_const())?;
            let string_length = bytes
                .iter()
                .position(|byte| *byte == 0)
                .ok_or(SourceHeapError::MissingNulTerminator)?;
            (
                bytes[..string_length]
                    .iter()
                    .any(|byte| *byte == b'\n' as i8),
                i32::try_from(string_length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
            )
        };
        *too_long_line = i32::from(!has_newline && string_length == content_limit);

        let mut length = 0;
        lrtrim(heap, line, Some(&mut length))?;
        if length != 0 {
            return Ok(length);
        }
    }
}

pub(crate) fn inchi_ios_getsTab1(
    heap: &mut SourceHeap,
    line: SourceMutPointer<i8>,
    len: i32,
    ios: &mut INCHI_IOSTREAM,
    too_long_line: &mut i32,
) -> Result<i32, SourceHeapError> {
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

    let helper_len = len
        .checked_sub(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let content_limit = len
        .checked_sub(2)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let result = inchi_ios_str_getsTab(heap, line, helper_len, ios)?;
    if result.is_null() {
        *too_long_line = 0;
        return Ok(-1);
    }

    let terminator_index =
        usize::try_from(helper_len).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    *heap
        .slice_mut(line)?
        .get_mut(terminator_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
    let (has_newline, string_length) = {
        let bytes = heap.slice(line.as_const())?;
        let string_length = bytes
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        (
            bytes[..string_length]
                .iter()
                .any(|byte| *byte == b'\n' as i8),
            i32::try_from(string_length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
        )
    };
    *too_long_line = i32::from(!has_newline && string_length == content_limit);

    let mut length = 0;
    lrtrim(heap, line, Some(&mut length))?;
    Ok(length)
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
                    SourceFormatArgument::Pointer(SourceConstPointer::<SourceVoid>::null()),
                    SourceFormatArgument::Pointer(SourceConstPointer::<SourceVoid>::null()),
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
                SourceFormatArgument::Pointer(SourceConstPointer::<SourceVoid>::null()),
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
        let mut pointer_arguments = SourceVaList {
            arguments: vec![SourceFormatArgument::Pointer(opaque.as_const())],
            position: 0,
        };
        assert_eq!(
            inchi_vfprintf(
                &mut heap,
                pointer_file,
                pointer_format.as_const(),
                &mut pointer_arguments,
            ),
            Err(SourceHeapError::UnsupportedSourceBehavior)
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
