use crate::source_types::{SourceConstPointer, SourceHeap, SourceMutPointer};
use crate::source_types::{
    BNS_ALTBOND_ERR, BNS_RADICAL_ERR, BNS_TIMEOUT, CT_ATOMCOUNT_ERR, CT_CALC_STEREO_ERR,
    CT_CANON_ERR, CT_ISO_H_ERR, CT_ISOCOUNT_ERR, CT_ISOTAUCOUNT_ERR, CT_LEN_MISMATCH,
    CT_MAPCOUNT_ERR, CT_OUT_OF_RAM, CT_OVERFLOW, CT_RANKING_ERR, CT_REMOVE_STEREO_ERR,
    CT_STEREO_CANON_ERR, CT_STEREOBOND_ERROR, CT_STEREOCOUNT_ERR, CT_TAUCOUNT_ERR, CT_TIMEOUT_ERR,
    CT_UNKNOWN_ERR, CT_USER_QUIT_ERR, CT_WRONG_FORMULA, STR_ERR_LEN, SourceHeapError,
};

// A temporary text destination: actual legacy pointer, local array or field.
// No allocation identity is inferred from a borrowed slice.
pub(crate) enum SourceTextTarget<'a> {
    Pointer(Option<SourceMutPointer<i8>>),
    Borrowed(&'a mut [i8]),
}

impl SourceTextTarget<'_> {
    pub(crate) fn has_target(&self) -> bool {
        // API Some(NULL) is distinct from None: keep its reached heap error.
        !matches!(self, Self::Pointer(None))
    }

    pub(crate) fn clear(&mut self, heap: &mut SourceHeap) -> Result<(), SourceHeapError> {
        // Native ichilnct.c231..234 and inchi_dll_b.c685/690/756/760/1898/1902/1961:
        // clear the actual destination at the reached source point. Preserve
        // legacy Some(NULL) heap admission; only the merge entry normalizes it.
        // INCHI✔️❗: O(1), no buffer/scan/registration; C513 costs OPEN.
        let bytes = match self {
            Self::Pointer(None) => return Ok(()),
            Self::Pointer(Some(pointer)) => heap.slice_mut(*pointer)?,
            Self::Borrowed(bytes) => &mut **bytes,
        };
        *bytes.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
        Ok(())
    }

    pub(crate) fn first_byte(&self, heap: &SourceHeap) -> Result<Option<i8>, SourceHeapError> {
        // Native ichilnct.c398: pStrErr && !pStrErr[0]. Keep the legacy
        // empty-slice observation (None), with no added extent admission.
        // INCHI✔️❗: O(1) read only; checked/whole-call costs OPEN.
        let bytes = match self {
            Self::Pointer(None) => return Ok(None),
            Self::Pointer(Some(pointer)) => heap.slice(pointer.as_const())?,
            Self::Borrowed(bytes) => &**bytes,
        };
        Ok(bytes.first().copied())
    }

    pub(crate) fn copy_from_pointer(
        &mut self, heap: &mut SourceHeap, source: SourceConstPointer<i8>, maximum_length: u32,
    ) -> Result<i32, SourceHeapError> {
        // Native util.c1760..1798 is the canonical mystrncpy implementation.
        // Pointer mode delegates its COMPLETE overlap/NULL/admission path;
        // the actual independent local arrays use its existing slice kernel.
        // INCHI✔️❗: same required strlen/copy/padding, no snapshot/extra scan;
        // source pointer includes its actual NUL. C513 whole-call costs OPEN.
        match self {
            Self::Pointer(None) => Ok(0),
            Self::Pointer(Some(pointer)) =>
                crate::source::base::util::mystrncpy(heap, *pointer, source, maximum_length),
            Self::Borrowed(bytes) => {
                if source.is_null() || maximum_length == 0 {
                    return Ok(0);
                }
                crate::source::base::util::mystrncpy_slice(
                    Some(&mut **bytes), Some(heap.slice(source)?), maximum_length)
            }
        }
    }

    pub(crate) fn copy_from_bytes(
        &mut self, heap: &mut SourceHeap, source: &[i8], maximum_length: u32,
    ) -> Result<i32, SourceHeapError> {
        // Native inchi_dll_b.c697/708/1925/1937 passes actual global-line
        // bytes to mystrncpy; no arena model/identity is required for them.
        // Preserve original pointer-target NULL/zero gates and reuse the ONE
        // mystrncpy_slice length/copy/padding algorithm. C542 costs stay OPEN.
        match self {
            Self::Pointer(None) => Ok(0),
            Self::Pointer(Some(p)) if p.is_null() || maximum_length == 0 => Ok(0),
            Self::Pointer(Some(p)) => crate::source::base::util::mystrncpy_slice(
                Some(heap.slice_mut(*p)?), Some(source), maximum_length,
            ),
            Self::Borrowed(bytes) => crate::source::base::util::mystrncpy_slice(
                Some(&mut **bytes), Some(source), maximum_length,
            ),
        }
    }

    pub(crate) fn trim(&mut self, heap: &mut SourceHeap) -> Result<i32, SourceHeapError> {
        // Native util.c1804..1835: preserve ONE canonical in-place lrtrim
        // algorithm; pointer mode retains its NULL/admission/overlap behavior.
        // INCHI✔️❗: original scan/copy only, no allocation; C513 costs OPEN.
        match self {
            Self::Pointer(None) => Ok(0),
            Self::Pointer(Some(pointer)) => {
                let mut length = 0;
                crate::source::base::util::lrtrim(heap, *pointer, Some(&mut length))?;
                Ok(length)
            }
            Self::Borrowed(bytes) => crate::source::base::util::lrtrim_bytes(&mut **bytes),
        }
    }

    pub(crate) fn append_from_distinct_pointer(
        &mut self, heap: &mut SourceHeap, source: SourceConstPointer<i8>, length: usize,
    ) -> Result<(), SourceHeapError> {
        // Native inchi_dll_b.c1961..1965:
        //                     pStrErr[0] = '\0';
        //                     c = *q;
        //                     *q = '\0';
        //                     TREAT_ERR( *err, 0, p );
        //                     *q = c;
        // INCHI✔️❗: Caller owns a fresh, distinct parser line allocation.
        // Preserve destination-before-source admission; same canonical byte
        // append supplies the implicit quote-end NUL. No line write/copy,
        // guessed origin or general overlap waiver. O(1) projection, C524 OPEN.
        match self {
            Self::Pointer(None) => Ok(()),
            Self::Pointer(Some(pointer)) => heap.with_slice_mut_and_heap(*pointer, |output, heap| {
                let message = heap.slice(source)?.get(..length)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                add_error_message_bytes(Some(output), bytemuck::cast_slice(message))?;
                Ok(())
            }),
            Self::Borrowed(output) => {
                let message = heap.slice(source)?.get(..length)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                add_error_message_bytes(Some(&mut **output), bytemuck::cast_slice(message))?;
                Ok(())
            }
        }
    }

    pub(crate) fn append_bytes(
        &mut self, heap: &mut SourceHeap, message: &[u8],
    ) -> Result<(), SourceHeapError> {
        // Native ichierr.c106..157 is implemented by the ONE canonical
        // add_error_message_bytes kernel below. Only project its destination;
        // preserve API heap admission before message gates for Some(pointer).
        // INCHI✔️❗: O(1) projection, no buffer/clone/scan/arena registration;
        // shared append and whole-call costs remain C512 OPEN.
        let output = match self {
            Self::Pointer(None) => return Ok(()),
            Self::Pointer(Some(pointer)) => heap.slice_mut(*pointer)?,
            Self::Borrowed(bytes) => &mut **bytes,
        };
        add_error_message_bytes(Some(output), message)?;
        Ok(())
    }
}

fn c_string_length(value: &[i8]) -> Result<usize, SourceHeapError> {
    value
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)
}

/// A private message projection; implicit NUL replaces an appended temporary.
#[derive(Clone, Copy)]
struct ErrorMessage<'a> {
    bytes: &'a [i8],
    implicit_final_nul: bool,
}

impl ErrorMessage<'_> {
    fn first(self) -> Result<i8, SourceHeapError> {
        self.bytes.first().copied()
            .or_else(|| self.implicit_final_nul.then_some(0))
            .ok_or(SourceHeapError::MissingNulTerminator)
    }

    fn length(self) -> Result<usize, SourceHeapError> {
        self.bytes.iter().position(|byte| *byte == 0)
            .or_else(|| self.implicit_final_nul.then_some(self.bytes.len()))
            .ok_or(SourceHeapError::MissingNulTerminator)
    }
}

fn find_bytes(haystack: &[i8], needle: &[i8]) -> Option<usize> {
    // INCHI source ichierr.c:146:
    //     if (strstr( all_messages, "..." ))
    // INCHI✔️❗: These slices already end at their known first NUL. Borrow
    // their byte representation for the existing linear, allocation-free search.
    memchr::memmem::find(bytemuck::cast_slice(haystack), bytemuck::cast_slice(needle))
}

fn find_message_in_c_string(
    previous: &[i8],
    message: ErrorMessage<'_>,
) -> Result<(Option<usize>, usize), SourceHeapError> {
    // INCHI source ichierr.c:164:
    //     char *p = strstr( prev_messages, new_message );
    // INCHI✔️❗: Preserve the first C-string match without a full previous
    // string preflight. Existing memchr Finder borrows the needle; no heap
    // storage. Windows overlap by M-1 and advance by at least M, so total
    // search/terminator work is O(N+M), with constant scratch. Chunk startup,
    // safe bounds and configured libc-relative whole-call costs remain open.
    let first = message.first()?;
    if first == 0 {
        return Ok((Some(0), 0));
    }
    let previous_bytes: &[u8] = bytemuck::cast_slice(previous);
    let Some(first_position) = memchr::memchr2(0, first as u8, previous_bytes) else {
        return Err(SourceHeapError::MissingNulTerminator);
    };
    if previous_bytes[first_position] == 0 {
        return Ok((None, 0));
    }
    let message_length = message.length()?;
    let message: &[u8] = bytemuck::cast_slice(&message.bytes[..message_length]);
    let previous: &[u8] = bytemuck::cast_slice(previous);
    let finder = memchr::memmem::Finder::new(message);
    // This is a search window, never an input-size cap or fallback domain.
    let window_length = message_length.max(256);
    let mut scanned = first_position;
    loop {
        if scanned == previous.len() {
            return Err(SourceHeapError::MissingNulTerminator);
        }
        let next = scanned.saturating_add(window_length).min(previous.len());
        let terminator = memchr::memchr(0, &previous[scanned..next])
            .map(|offset| scanned + offset);
        let end = terminator.unwrap_or(next);
        let begin = scanned.saturating_sub(message_length - 1).max(first_position);
        if let Some(offset) = finder.find(&previous[begin..end]) {
            return Ok((Some(begin + offset), message_length));
        }
        if terminator.is_some() {
            return Ok((None, message_length));
        }
        scanned = next;
    }
}

pub(crate) fn find_c_substring(previous: &[i8], needle: &[u8]) -> Result<Option<usize>, SourceHeapError> {
    // INCHI source ichiprt1.c:5066,5070:
    //     if (!strstr(s, "/z"))
    //     s0 = strstr(s, "InChI=1B/");
    // INCHI✔️❗: One canonical first-match/NUL search, borrowing a literal
    // needle with implicit NUL. O(N+M), constant scratch, no allocation;
    // libc-relative whole-call costs remain OPEN with the existing kernel.
    find_message_in_c_string(previous, ErrorMessage {
        bytes: bytemuck::cast_slice(needle), implicit_final_nul: true,
    }).map(|(position, _)| position)
}

fn error_byte_at(value: &[i8], index: i64) -> Result<i8, SourceHeapError> {
    let index = usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    value.get(index).copied().ok_or(SourceHeapError::PointerOutOfBounds)
}

fn append_error_part(
    destination: &mut [i8],
    end: &mut usize,
    part: &[i8],
) -> Result<(), SourceHeapError> {
    // INCHI source ichierr.c:137,139,141:
    //                 strcat(all_messages, ";");
    //             strcat(all_messages, " ");
    //         strcat(all_messages, new_message);
    // INCHI✔️❗: Each source strcat uses the actual NUL position, independently
    // of narrowed room arithmetic, and installs its own trailing NUL. Safe
    // extent errors are checked only for the reached append; no allocation.
    let next = end.checked_add(part.len()).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let required = next.checked_add(1).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let output = destination.get_mut(*end..required).ok_or(SourceHeapError::PointerOutOfBounds)?;
    output[..part.len()].copy_from_slice(part);
    output[part.len()] = 0;
    *end = next;
    Ok(())
}

pub(crate) fn err_msg_text(error_code: i32, scratch: &mut [u8; 64]) -> &str {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.c:54 ErrMsg
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: ErrMsg
    // INCHI✔️❌: const char *ErrMsg( int nErrorCode )
    // INCHI✔️❌: {
    // INCHI✔️❌:     const char *p;
    // INCHI✔️❌:     static char szErrMsg[64];
    // INCHI✔️❌:     switch (nErrorCode)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         case 0:                      p = "";                      break;
    // INCHI✔️❌:         case CT_OVERFLOW:            p = "ARRAY OVERFLOW";        break;
    // INCHI✔️❌:         case CT_LEN_MISMATCH:        p = "LENGTH_MISMATCH";       break;
    // INCHI✔️❌:         case CT_OUT_OF_RAM:          p = "Out of RAM";            break;
    // INCHI✔️❌:         case CT_RANKING_ERR:         p = "RANKING_ERR";           break;
    // INCHI✔️❌:         case CT_ISOCOUNT_ERR:        p = "ISOCOUNT_ERR";          break;
    // INCHI✔️❌:         case CT_TAUCOUNT_ERR:        p = "TAUCOUNT_ERR";          break;
    // INCHI✔️❌:         case CT_ISOTAUCOUNT_ERR:     p = "ISOTAUCOUNT_ERR";       break;
    // INCHI✔️❌:         case CT_MAPCOUNT_ERR:        p = "MAPCOUNT_ERR";          break;
    // INCHI✔️❌:         case CT_TIMEOUT_ERR:         p = "Time limit exceeded";   break;
    // INCHI✔️❌:         case CT_ISO_H_ERR:           p = "ISO_H_ERR";             break;
    // INCHI✔️❌:         case CT_STEREOCOUNT_ERR:     p = "STEREOCOUNT_ERR";       break;
    // INCHI✔️❌:         case CT_ATOMCOUNT_ERR:       p = "ATOMCOUNT_ERR";         break;
    // INCHI✔️❌:         case CT_STEREOBOND_ERROR:    p = "STEREOBOND_ERR";        break;
    // INCHI✔️❌:         case CT_USER_QUIT_ERR:       p = "User requested termination"; break;
    // INCHI✔️❌:         case CT_REMOVE_STEREO_ERR:   p = "REMOVE_STEREO_ERR";     break;
    // INCHI✔️❌:         case CT_CALC_STEREO_ERR:     p = "CALC_STEREO_ERR";       break;
    // INCHI✔️❌:         case CT_STEREO_CANON_ERR:    p = "STEREO_CANON_ERR";      break;
    // INCHI✔️❌:         case CT_CANON_ERR:           p = "CANON_ERR";             break;
    // INCHI✔️❌:         case CT_WRONG_FORMULA:       p = "Wrong or missing chemical formula";  break;
    // INCHI✔️❌:         /*case CT_CANON_ERR2:          p = "CT_CANON_ERR2";         break;*/
    // INCHI✔️❌:         case CT_UNKNOWN_ERR:         p = "UNKNOWN_ERR";           break;
    // INCHI✔️❌:         case BNS_RADICAL_ERR:        p = "Cannot process free radical center"; break;
    // INCHI✔️❌:         case BNS_ALTBOND_ERR:        p = "Cannot process aromatic bonds";      break;
    // INCHI✔️❌:         /* v. 1.05 */
    // INCHI✔️❌:         case BNS_TIMEOUT:             p = "Structure normalization timeout";      break;
    // INCHI✔️❌:
    // INCHI✔️❌:         default:
    // INCHI✔️❌:             if (nErrorCode > CT_UNKNOWN_ERR)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 sprintf(szErrMsg, "No description(%d)", nErrorCode);
    // INCHI✔️❌:                 p = szErrMsg;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 sprintf(szErrMsg, "UNKNOWN_ERR(%d)", CT_UNKNOWN_ERR - nErrorCode);
    // INCHI✔️❌:                 p = szErrMsg;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             break;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return p;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: ErrMsg
    // END INCHI C FUNCTION: ErrMsg
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: ErrMsg
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️❌: Owned Rust text adds allocation versus C literals/static szErrMsg storage.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: ErrMsg
    // INCHI✔️❗: One canonical switch; known messages borrow static literals.
    // Defaults use source-sized64byte scratch and shared exact signed decimal,
    // no format!/String/Vec staging. O(1) bounded time/space, zero allocation.
    // Checked UTF8 on constructed ASCII and scratch initialization cost OPEN.
    match error_code {
        0 => "",
        CT_OVERFLOW => "ARRAY OVERFLOW",
        CT_LEN_MISMATCH => "LENGTH_MISMATCH",
        CT_OUT_OF_RAM => "Out of RAM",
        CT_RANKING_ERR => "RANKING_ERR",
        CT_ISOCOUNT_ERR => "ISOCOUNT_ERR",
        CT_TAUCOUNT_ERR => "TAUCOUNT_ERR",
        CT_ISOTAUCOUNT_ERR => "ISOTAUCOUNT_ERR",
        CT_MAPCOUNT_ERR => "MAPCOUNT_ERR",
        CT_TIMEOUT_ERR => "Time limit exceeded",
        CT_ISO_H_ERR => "ISO_H_ERR",
        CT_STEREOCOUNT_ERR => "STEREOCOUNT_ERR",
        CT_ATOMCOUNT_ERR => "ATOMCOUNT_ERR",
        CT_STEREOBOND_ERROR => "STEREOBOND_ERR",
        CT_USER_QUIT_ERR => "User requested termination",
        CT_REMOVE_STEREO_ERR => "REMOVE_STEREO_ERR",
        CT_CALC_STEREO_ERR => "CALC_STEREO_ERR",
        CT_STEREO_CANON_ERR => "STEREO_CANON_ERR",
        CT_CANON_ERR => "CANON_ERR",
        CT_WRONG_FORMULA => "Wrong or missing chemical formula",
        CT_UNKNOWN_ERR => "UNKNOWN_ERR",
        BNS_RADICAL_ERR => "Cannot process free radical center",
        BNS_ALTBOND_ERR => "Cannot process aromatic bonds",
        BNS_TIMEOUT => "Structure normalization timeout",
        code => {
            let (prefix, number) = if code > CT_UNKNOWN_ERR {
                (b"No description(".as_slice(), code)
            } else {
                (b"UNKNOWN_ERR(".as_slice(), CT_UNKNOWN_ERR.wrapping_sub(code))
            };
            let mut decimal = [0_u8; 21];
            let digits = crate::source::base::ichi_io::source_decimal_i64(
                i64::from(number), None, &mut decimal,
            );
            let end = prefix.len() + digits.len();
            // i32 needs at most11 decimal bytes; prefix15 +11 +')'+NUL <64.
            scratch[..prefix.len()].copy_from_slice(prefix);
            scratch[prefix.len()..end].copy_from_slice(digits);
            scratch[end] = b')';
            scratch[end + 1] = 0;
            std::str::from_utf8(&scratch[..end + 1])
                .expect("source error literals and decimal digits are ASCII")
        }
    }
}

#[allow(non_snake_case)]
pub(crate) fn ErrMsg(error_code: i32) -> String {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.c:54 ErrMsg
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: ErrMsg
    // INCHI✔️❌: const char *ErrMsg( int nErrorCode )
    // INCHI✔️❌: {
    // INCHI✔️❌:     const char *p;
    // INCHI✔️❌:     static char szErrMsg[64];
    // INCHI✔️❌:     switch (nErrorCode)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         case 0:                      p = "";                      break;
    // INCHI✔️❌:         case CT_OVERFLOW:            p = "ARRAY OVERFLOW";        break;
    // INCHI✔️❌:         case CT_LEN_MISMATCH:        p = "LENGTH_MISMATCH";       break;
    // INCHI✔️❌:         case CT_OUT_OF_RAM:          p = "Out of RAM";            break;
    // INCHI✔️❌:         case CT_RANKING_ERR:         p = "RANKING_ERR";           break;
    // INCHI✔️❌:         case CT_ISOCOUNT_ERR:        p = "ISOCOUNT_ERR";          break;
    // INCHI✔️❌:         case CT_TAUCOUNT_ERR:        p = "TAUCOUNT_ERR";          break;
    // INCHI✔️❌:         case CT_ISOTAUCOUNT_ERR:     p = "ISOTAUCOUNT_ERR";       break;
    // INCHI✔️❌:         case CT_MAPCOUNT_ERR:        p = "MAPCOUNT_ERR";          break;
    // INCHI✔️❌:         case CT_TIMEOUT_ERR:         p = "Time limit exceeded";   break;
    // INCHI✔️❌:         case CT_ISO_H_ERR:           p = "ISO_H_ERR";             break;
    // INCHI✔️❌:         case CT_STEREOCOUNT_ERR:     p = "STEREOCOUNT_ERR";       break;
    // INCHI✔️❌:         case CT_ATOMCOUNT_ERR:       p = "ATOMCOUNT_ERR";         break;
    // INCHI✔️❌:         case CT_STEREOBOND_ERROR:    p = "STEREOBOND_ERR";        break;
    // INCHI✔️❌:         case CT_USER_QUIT_ERR:       p = "User requested termination"; break;
    // INCHI✔️❌:         case CT_REMOVE_STEREO_ERR:   p = "REMOVE_STEREO_ERR";     break;
    // INCHI✔️❌:         case CT_CALC_STEREO_ERR:     p = "CALC_STEREO_ERR";       break;
    // INCHI✔️❌:         case CT_STEREO_CANON_ERR:    p = "STEREO_CANON_ERR";      break;
    // INCHI✔️❌:         case CT_CANON_ERR:           p = "CANON_ERR";             break;
    // INCHI✔️❌:         case CT_WRONG_FORMULA:       p = "Wrong or missing chemical formula";  break;
    // INCHI✔️❌:         /*case CT_CANON_ERR2:          p = "CT_CANON_ERR2";         break;*/
    // INCHI✔️❌:         case CT_UNKNOWN_ERR:         p = "UNKNOWN_ERR";           break;
    // INCHI✔️❌:         case BNS_RADICAL_ERR:        p = "Cannot process free radical center"; break;
    // INCHI✔️❌:         case BNS_ALTBOND_ERR:        p = "Cannot process aromatic bonds";      break;
    // INCHI✔️❌:         /* v. 1.05 */
    // INCHI✔️❌:         case BNS_TIMEOUT:             p = "Structure normalization timeout";      break;
    // INCHI✔️❌:
    // INCHI✔️❌:         default:
    // INCHI✔️❌:             if (nErrorCode > CT_UNKNOWN_ERR)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 sprintf(szErrMsg, "No description(%d)", nErrorCode);
    // INCHI✔️❌:                 p = szErrMsg;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 sprintf(szErrMsg, "UNKNOWN_ERR(%d)", CT_UNKNOWN_ERR - nErrorCode);
    // INCHI✔️❌:                 p = szErrMsg;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             break;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return p;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: ErrMsg
    // END INCHI C FUNCTION: ErrMsg
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: ErrMsg
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️❌: Owned Rust text adds allocation versus C literals/static szErrMsg storage.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: ErrMsg
    // INCHI✔️❌: Original owned String signature retained. Its one output
    // allocation/copy versus native static pointer is explicitly deferred in
    // INCHI-PERF-ERRMSG-OWNED; not a no-worse-cost or allocator acceptance.
    // Shared borrowed callers avoid this cost without changing existing API.
    let mut scratch = [0_u8; 64];
    err_msg_text(error_code, &mut scratch).to_owned()
}

pub(crate) fn already_have_this_message(
    prev_messages: &[i8],
    new_message: &[i8],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.c:160 already_have_this_message
    // INCHI✔️❗: int already_have_this_message( char *prev_messages, const char *new_message )
    // INCHI✔️❗: {
    // INCHI✔️❗:     int have = 0;
    // INCHI✔️❗:
    // INCHI✔️❗:     char *p = strstr( prev_messages, new_message );
    // INCHI✔️❗:
    // INCHI✔️❗:     if (p)
    // INCHI✔️❗:     {
    // INCHI✔️❗:         have = ( p == prev_messages || (*( p - 1 ) == ' ' && ( *( p - 2 ) == ';' || *( p - 2 ) == ':' )) ); /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❗:         if (have)
    // INCHI✔️❗:         {
    // INCHI✔️❗:             int len_prev = (int) strlen( prev_messages );
    // INCHI✔️❗:             int len = (int) strlen( new_message );
    // INCHI✔️❗:             have = ( p + len == prev_messages + len_prev || (p[len] == ';' && p[len + 1] == ' ') || (p[len - 1] == ':' && p[len] == ' ') ); /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❗:         }
    // INCHI✔️❗:     }
    // INCHI✔️❗:
    // INCHI✔️❗:     return have;
    // INCHI✔️❗: }
    // END INCHI C FUNCTION: already_have_this_message

    already_have_message(prev_messages, ErrorMessage {
        bytes: new_message,
        implicit_final_nul: false,
    })
}

fn already_have_message(
    prev_messages: &[i8],
    new_message: ErrorMessage<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.c:160 already_have_this_message
    // INCHI✔️❗: int already_have_this_message( char *prev_messages, const char *new_message )
    // INCHI✔️❗: {
    // INCHI✔️❗:     int have = 0;
    // INCHI✔️❗:
    // INCHI✔️❗:     char *p = strstr( prev_messages, new_message );
    // INCHI✔️❗:
    // INCHI✔️❗:     if (p)
    // INCHI✔️❗:     {
    // INCHI✔️❗:         have = ( p == prev_messages || (*( p - 1 ) == ' ' && ( *( p - 2 ) == ';' || *( p - 2 ) == ':' )) ); /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❗:         if (have)
    // INCHI✔️❗:         {
    // INCHI✔️❗:             int len_prev = (int) strlen( prev_messages );
    // INCHI✔️❗:             int len = (int) strlen( new_message );
    // INCHI✔️❗:             have = ( p + len == prev_messages + len_prev || (p[len] == ';' && p[len + 1] == ' ') || (p[len - 1] == ':' && p[len] == ' ') ); /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❗:         }
    // INCHI✔️❗:     }
    // INCHI✔️❗:
    // INCHI✔️❗:     return have;
    // INCHI✔️❗: }
    // END INCHI C FUNCTION: already_have_this_message

    let (position, message_length) = find_message_in_c_string(prev_messages, new_message)?;
    let Some(position) = position else {
        return Ok(0);
    };
    let starts_at_message_boundary = position == 0
        || (position >= 2
            && prev_messages[position - 1] == b' ' as i8
            && matches!(prev_messages[position - 2] as u8, b';' | b':'));
    if !starts_at_message_boundary {
        return Ok(0);
    }

    // GCC/LP64 casts strlen to 32-bit int; pointer offsets use those narrowed
    // values, not the physical string length. Keep negative offsets checked.
    let previous_length = c_string_length(prev_messages)? as i32;
    let message_length = message_length as i32;
    // A Rust slice's byte length is at most isize::MAX; source-defined pointer
    // offsets therefore fit i64 here. Reject only an out-of-range sum.
    let end = (position as i64)
        .checked_add(i64::from(message_length))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if end == i64::from(previous_length) {
        return Ok(1);
    }
    let at_end = error_byte_at(prev_messages, end)?;
    if at_end == b';' as i8 && error_byte_at(prev_messages, end + 1)? == b' ' as i8 {
        return Ok(1);
    }
    // The native empty-needle p[-1] case is undefined when p is the base;
    // retain the previous safe non-match for that unmodeled source state.
    if end == 0 && message_length == 0 {
        return Ok(0);
    }
    Ok(i32::from(
        error_byte_at(prev_messages, end - 1)? == b':' as i8 && at_end == b' ' as i8,
    ))
}

pub(crate) fn AddErrorMessage(
    all_messages: Option<&mut [i8]>,
    new_message: Option<&[i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.c:106 AddErrorMessage
    // INCHI✔️❗: int AddErrorMessage( char *all_messages, const char *new_message )
    // INCHI✔️❗: {
    // INCHI✔️❗:     int len_all, len;
    // INCHI✔️❗:
    // INCHI✔️❗:     if (!all_messages)
    // INCHI✔️❗:     {
    // INCHI✔️❗:         return 0;
    // INCHI✔️❗:     }
    // INCHI✔️❗:     if (!new_message)
    // INCHI✔️❗:     {
    // INCHI✔️❗:         return 0;
    // INCHI✔️❗:     }
    // INCHI✔️❗:     if (!new_message[0])
    // INCHI✔️❗:     {
    // INCHI✔️❗:         return 0;
    // INCHI✔️❗:     }
    // INCHI✔️❗:     if (already_have_this_message( all_messages, new_message ))
    // INCHI✔️❗:     {
    // INCHI✔️❗:         return 1;
    // INCHI✔️❗:     }
    // INCHI✔️❗:
    // INCHI✔️❗:     len_all = (int) strlen( all_messages );
    // INCHI✔️❗:     len = (int) strlen( new_message );
    // INCHI✔️❗:
    // INCHI✔️❗:     if (len_all + len + 2 * ( len_all > 0 ) < STR_ERR_LEN)
    // INCHI✔️❗:     {
    // INCHI✔️❗:         /* enough room... add message and return */
    // INCHI✔️❗:         if (len_all > 0)
    // INCHI✔️❗:         {
    // INCHI✔️❗:             if (all_messages[len_all - 1] != ':')
    // INCHI✔️❗:             {
    // INCHI✔️❗:                 strcat(all_messages, ";");
    // INCHI✔️❗:             }
    // INCHI✔️❗:             strcat(all_messages, " ");
    // INCHI✔️❗:         }
    // INCHI✔️❗:         strcat(all_messages, new_message);
    // INCHI✔️❗:         return 1;
    // INCHI✔️❗:     }
    // INCHI✔️❗:
    // INCHI✔️❗:     /*  not enough room... add no-room marker if not yet added */
    // INCHI✔️❗:     if (strstr( all_messages, "..." ))
    // INCHI✔️❗:     {
    // INCHI✔️❗:         return 0;
    // INCHI✔️❗:     }
    // INCHI✔️❗:     if (len_all + 3 < STR_ERR_LEN)
    // INCHI✔️❗:     {
    // INCHI✔️❗:         strcat(all_messages, "...");
    // INCHI✔️❗:     }
    // INCHI✔️❗:
    // INCHI✔️❗:     return 0;
    // INCHI✔️❗: }
    // END INCHI C FUNCTION: AddErrorMessage

    add_error_message(all_messages, new_message.map(|bytes| ErrorMessage {
        bytes,
        implicit_final_nul: false,
    }))
}

/// Source AddErrorMessage with the same message as both pointer arguments.
pub(crate) fn add_error_message_to_self(
    message: Option<&[i8]>,
) -> Result<i32, SourceHeapError> {
    // INCHI source ichierr.c:106..157; mol_fmt3.c TREAT_ERR_AND_FIN self-message:
    // int AddErrorMessage( char *all_messages, const char *new_message )
    // {
    //     int len_all, len;
    // 
    //     if (!all_messages)
    //     {
    //         return 0;
    //     }
    //     if (!new_message)
    //     {
    //         return 0;
    //     }
    //     if (!new_message[0])
    //     {
    //         return 0;
    //     }
    //     if (already_have_this_message( all_messages, new_message ))
    //     {
    //         return 1;
    //     }
    // 
    //     len_all = (int) strlen( all_messages );
    //     len = (int) strlen( new_message );
    // 
    //     if (len_all + len + 2 * ( len_all > 0 ) < STR_ERR_LEN)
    //     {
    //         /* enough room... add message and return */
    //         if (len_all > 0)
    //         {
    //             if (all_messages[len_all - 1] != ':')
    //             {
    //                 strcat(all_messages, ";");
    //             }
    //             strcat(all_messages, " ");
    //         }
    //         strcat(all_messages, new_message);
    //         return 1;
    //     }
    // 
    //     /*  not enough room... add no-room marker if not yet added */
    //     if (strstr( all_messages, "..." ))
    //     {
    //         return 0;
    //     }
    //     if (len_all + 3 < STR_ERR_LEN)
    //     {
    //         strcat(all_messages, "...");
    //     }
    // 
    //     return 0;
    // }
    // 
    // INCHI✔️❗: A nonempty string matches itself at zero, with equal narrowed
    // lengths, so the canonical duplicate reader returns before any append.
    // Preserve its missing-NUL/extent errors; no guessed successful return.
    // O(L)/O(1), no snapshot/allocation. Shared reader and whole-call costs open.
    let Some(message) = message else {
        return Ok(0);
    };
    let view = ErrorMessage {
        bytes: message,
        implicit_final_nul: false,
    };
    if view.first()? == 0 {
        return Ok(0);
    }
    already_have_message(message, view)
}

pub(crate) fn add_error_message_bytes(
    all_messages: Option<&mut [i8]>,
    new_message: &[u8],
) -> Result<i32, SourceHeapError> {
    // INCHI source ichierr.h:173:
    //         if ( !(err) && (new_err) ) { (err) = (new_err);} AddErrorMessage(pStrErr, (msg))
    // INCHI✔️❗: Borrow literal/detached message bytes. The helper's appended
    // final zero is implicit, while an earlier zero still ends the message.
    // Shared native operation/search bodies retain every gate and write;
    // no temporary buffer. View/checked/whole-call costs remain open.
    add_error_message(all_messages, Some(ErrorMessage {
        bytes: bytemuck::cast_slice(new_message),
        implicit_final_nul: true,
    }))
}

fn add_error_message(
    all_messages: Option<&mut [i8]>,
    new_message: Option<ErrorMessage<'_>>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.c:106 AddErrorMessage
    // INCHI✔️❗: int AddErrorMessage( char *all_messages, const char *new_message )
    // INCHI✔️❗: {
    // INCHI✔️❗:     int len_all, len;
    // INCHI✔️❗:
    // INCHI✔️❗:     if (!all_messages)
    // INCHI✔️❗:     {
    // INCHI✔️❗:         return 0;
    // INCHI✔️❗:     }
    // INCHI✔️❗:     if (!new_message)
    // INCHI✔️❗:     {
    // INCHI✔️❗:         return 0;
    // INCHI✔️❗:     }
    // INCHI✔️❗:     if (!new_message[0])
    // INCHI✔️❗:     {
    // INCHI✔️❗:         return 0;
    // INCHI✔️❗:     }
    // INCHI✔️❗:     if (already_have_this_message( all_messages, new_message ))
    // INCHI✔️❗:     {
    // INCHI✔️❗:         return 1;
    // INCHI✔️❗:     }
    // INCHI✔️❗:
    // INCHI✔️❗:     len_all = (int) strlen( all_messages );
    // INCHI✔️❗:     len = (int) strlen( new_message );
    // INCHI✔️❗:
    // INCHI✔️❗:     if (len_all + len + 2 * ( len_all > 0 ) < STR_ERR_LEN)
    // INCHI✔️❗:     {
    // INCHI✔️❗:         /* enough room... add message and return */
    // INCHI✔️❗:         if (len_all > 0)
    // INCHI✔️❗:         {
    // INCHI✔️❗:             if (all_messages[len_all - 1] != ':')
    // INCHI✔️❗:             {
    // INCHI✔️❗:                 strcat(all_messages, ";");
    // INCHI✔️❗:             }
    // INCHI✔️❗:             strcat(all_messages, " ");
    // INCHI✔️❗:         }
    // INCHI✔️❗:         strcat(all_messages, new_message);
    // INCHI✔️❗:         return 1;
    // INCHI✔️❗:     }
    // INCHI✔️❗:
    // INCHI✔️❗:     /*  not enough room... add no-room marker if not yet added */
    // INCHI✔️❗:     if (strstr( all_messages, "..." ))
    // INCHI✔️❗:     {
    // INCHI✔️❗:         return 0;
    // INCHI✔️❗:     }
    // INCHI✔️❗:     if (len_all + 3 < STR_ERR_LEN)
    // INCHI✔️❗:     {
    // INCHI✔️❗:         strcat(all_messages, "...");
    // INCHI✔️❗:     }
    // INCHI✔️❗:
    // INCHI✔️❗:     return 0;
    // INCHI✔️❗: }
    // END INCHI C FUNCTION: AddErrorMessage

    let Some(all_messages) = all_messages else {
        return Ok(0);
    };
    let Some(new_message) = new_message else {
        return Ok(0);
    };
    if new_message.first()? == 0 {
        return Ok(0);
    }
    if already_have_message(all_messages, new_message)? != 0 {
        return Ok(1);
    }

    let all_length = c_string_length(all_messages)?;
    let message_length = new_message.length()?;
    let native_all_length = all_length as i32;
    let native_message_length = message_length as i32;
    // Native signed addition overflow is outside defined source behavior;
    // reject that state structurally, rather than admitting a heuristic room.
    let room_sum = native_all_length
        .checked_add(native_message_length)
        .and_then(|sum| sum.checked_add(2 * i32::from(native_all_length > 0)))
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    if room_sum < STR_ERR_LEN as i32 {
        let mut output = all_length;
        if native_all_length > 0 {
            if all_messages[native_all_length as usize - 1] != b':' as i8 {
                append_error_part(all_messages, &mut output, &[b';' as i8])?;
            }
            append_error_part(all_messages, &mut output, &[b' ' as i8])?;
        }
        // strcat copies the entire real string even when int(strlen) wraps.
        append_error_part(all_messages, &mut output, &new_message.bytes[..message_length])?;
        return Ok(1);
    }
    if find_bytes(&all_messages[..all_length], &[b'.' as i8; 3]).is_some() {
        return Ok(0);
    }
    if native_all_length.checked_add(3).ok_or(SourceHeapError::SourceIntegerOverflow)?
        < STR_ERR_LEN as i32
    {
        let mut output = all_length;
        append_error_part(all_messages, &mut output, &[b'.' as i8; 3])?;
    }
    Ok(0)
}

#[cfg(test)]
#[allow(non_snake_case)]
mod tests {
    use super::*;

    fn c_bytes(value: &str) -> Vec<i8> {
        value.bytes().map(|byte| byte as i8).chain([0]).collect()
    }

    fn error_buffer(value: &str) -> [i8; STR_ERR_LEN as usize] {
        let mut result = [0; STR_ERR_LEN as usize];
        let value = c_bytes(value);
        result[..value.len()].copy_from_slice(&value);
        result
    }

    fn buffer_text(value: &[i8]) -> String {
        let length = c_string_length(value).unwrap();
        String::from_utf8(value[..length].iter().map(|byte| *byte as u8).collect()).unwrap()
    }

    #[test]
    fn source_port__ichierr__errmsg__line_54() {
        for (code, expected) in [
            (0, ""),
            (CT_OVERFLOW, "ARRAY OVERFLOW"),
            (CT_LEN_MISMATCH, "LENGTH_MISMATCH"),
            (CT_OUT_OF_RAM, "Out of RAM"),
            (CT_RANKING_ERR, "RANKING_ERR"),
            (CT_ISOCOUNT_ERR, "ISOCOUNT_ERR"),
            (CT_TAUCOUNT_ERR, "TAUCOUNT_ERR"),
            (CT_ISOTAUCOUNT_ERR, "ISOTAUCOUNT_ERR"),
            (CT_MAPCOUNT_ERR, "MAPCOUNT_ERR"),
            (CT_TIMEOUT_ERR, "Time limit exceeded"),
            (CT_ISO_H_ERR, "ISO_H_ERR"),
            (CT_STEREOCOUNT_ERR, "STEREOCOUNT_ERR"),
            (CT_ATOMCOUNT_ERR, "ATOMCOUNT_ERR"),
            (CT_STEREOBOND_ERROR, "STEREOBOND_ERR"),
            (CT_USER_QUIT_ERR, "User requested termination"),
            (CT_REMOVE_STEREO_ERR, "REMOVE_STEREO_ERR"),
            (CT_CALC_STEREO_ERR, "CALC_STEREO_ERR"),
            (CT_STEREO_CANON_ERR, "STEREO_CANON_ERR"),
            (CT_CANON_ERR, "CANON_ERR"),
            (CT_WRONG_FORMULA, "Wrong or missing chemical formula"),
            (CT_UNKNOWN_ERR, "UNKNOWN_ERR"),
            (BNS_RADICAL_ERR, "Cannot process free radical center"),
            (BNS_ALTBOND_ERR, "Cannot process aromatic bonds"),
            (BNS_TIMEOUT, "Structure normalization timeout"),
        ] {
            assert_eq!(ErrMsg(code), expected, "error code {code}");
        }

        assert_eq!(ErrMsg(-29_999), "No description(-29999)");
        assert_eq!(ErrMsg(1), "No description(1)");
        assert_eq!(ErrMsg(i32::MAX), "No description(2147483647)");
        assert_eq!(ErrMsg(CT_UNKNOWN_ERR - 1), "UNKNOWN_ERR(1)");
        assert_eq!(ErrMsg(i32::MIN), "UNKNOWN_ERR(2147453629)");
    }

    #[test]
    fn source_port__ichierr__already_have_this_message__line_160() {
        assert_eq!(
            already_have_this_message(&c_bytes("alpha; beta"), &c_bytes("alpha")),
            Ok(1)
        );
        assert_eq!(
            already_have_this_message(&c_bytes("alpha; beta"), &c_bytes("beta")),
            Ok(1)
        );
        assert_eq!(
            already_have_this_message(&c_bytes("prefix alphabet"), &c_bytes("alpha")),
            Ok(0)
        );
        assert_eq!(
            already_have_this_message(&c_bytes("head: detail; tail"), &c_bytes("detail")),
            Ok(1)
        );
        assert_eq!(
            already_have_this_message(&c_bytes("head: detail text"), &c_bytes("head:")),
            Ok(1)
        );
        assert_eq!(
            already_have_this_message(&c_bytes("alphabet; alpha"), &c_bytes("alpha")),
            Ok(0),
            "the C function checks only the first strstr match"
        );
        assert_eq!(
            already_have_this_message(&[b'a' as i8], &c_bytes("a")),
            Err(SourceHeapError::MissingNulTerminator)
        );
    }

    #[test]
    fn source_port__ichierr__adderrormessage__line_106() {
        let message = c_bytes("first");
        assert_eq!(AddErrorMessage(None, Some(&message)), Ok(0));
        let mut messages = error_buffer("");
        assert_eq!(AddErrorMessage(Some(&mut messages), None), Ok(0));
        assert_eq!(
            AddErrorMessage(Some(&mut messages), Some(&c_bytes(""))),
            Ok(0)
        );

        assert_eq!(AddErrorMessage(Some(&mut messages), Some(&message)), Ok(1));
        assert_eq!(buffer_text(&messages), "first");
        assert_eq!(AddErrorMessage(Some(&mut messages), Some(&message)), Ok(1));
        assert_eq!(buffer_text(&messages), "first");
        assert_eq!(
            AddErrorMessage(Some(&mut messages), Some(&c_bytes("second"))),
            Ok(1)
        );
        assert_eq!(buffer_text(&messages), "first; second");

        let mut colon = error_buffer("header:");
        assert_eq!(
            AddErrorMessage(Some(&mut colon), Some(&c_bytes("detail"))),
            Ok(1)
        );
        assert_eq!(buffer_text(&colon), "header: detail");

        let mut no_room = error_buffer(&"x".repeat(253));
        assert_eq!(
            AddErrorMessage(Some(&mut no_room), Some(&c_bytes("more"))),
            Ok(0)
        );
        assert_eq!(buffer_text(&no_room), "x".repeat(253));

        let mut marker_room = error_buffer(&"x".repeat(252));
        assert_eq!(
            AddErrorMessage(Some(&mut marker_room), Some(&c_bytes("more"))),
            Ok(0)
        );
        assert_eq!(buffer_text(&marker_room), format!("{}...", "x".repeat(252)));
        assert_eq!(
            AddErrorMessage(Some(&mut marker_room), Some(&c_bytes("again"))),
            Ok(0)
        );
        assert_eq!(buffer_text(&marker_room), format!("{}...", "x".repeat(252)));

        let mut short = [0_i8; 8];
        assert_eq!(
            AddErrorMessage(Some(&mut short), Some(&message)),
            // ichierr.c130..142: 0+5<256; strcat writes six bytes into eight.
            Ok(1)
        );
    }
}
