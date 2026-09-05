use crate::source_types::{
    BNS_ALTBOND_ERR, BNS_RADICAL_ERR, BNS_TIMEOUT, CT_ATOMCOUNT_ERR, CT_CALC_STEREO_ERR, CT_CANON_ERR, CT_ISO_H_ERR,
    CT_ISOCOUNT_ERR, CT_ISOTAUCOUNT_ERR, CT_LEN_MISMATCH, CT_MAPCOUNT_ERR, CT_OUT_OF_RAM, CT_OVERFLOW, CT_RANKING_ERR,
    CT_REMOVE_STEREO_ERR, CT_STEREO_CANON_ERR, CT_STEREOBOND_ERROR, CT_STEREOCOUNT_ERR, CT_TAUCOUNT_ERR,
    CT_TIMEOUT_ERR, CT_UNKNOWN_ERR, CT_USER_QUIT_ERR, CT_WRONG_FORMULA, STR_ERR_LEN, SourceHeapError,
};

fn c_string_length(value: &[i8]) -> Result<usize, SourceHeapError> {
    value
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)
}

fn find_bytes(haystack: &[i8], needle: &[i8]) -> Option<usize> {
    if needle.is_empty() {
        return Some(0);
    }
    haystack.windows(needle.len()).position(|window| window == needle)
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
    match error_code {
        0 => "".to_owned(),
        CT_OVERFLOW => "ARRAY OVERFLOW".to_owned(),
        CT_LEN_MISMATCH => "LENGTH_MISMATCH".to_owned(),
        CT_OUT_OF_RAM => "Out of RAM".to_owned(),
        CT_RANKING_ERR => "RANKING_ERR".to_owned(),
        CT_ISOCOUNT_ERR => "ISOCOUNT_ERR".to_owned(),
        CT_TAUCOUNT_ERR => "TAUCOUNT_ERR".to_owned(),
        CT_ISOTAUCOUNT_ERR => "ISOTAUCOUNT_ERR".to_owned(),
        CT_MAPCOUNT_ERR => "MAPCOUNT_ERR".to_owned(),
        CT_TIMEOUT_ERR => "Time limit exceeded".to_owned(),
        CT_ISO_H_ERR => "ISO_H_ERR".to_owned(),
        CT_STEREOCOUNT_ERR => "STEREOCOUNT_ERR".to_owned(),
        CT_ATOMCOUNT_ERR => "ATOMCOUNT_ERR".to_owned(),
        CT_STEREOBOND_ERROR => "STEREOBOND_ERR".to_owned(),
        CT_USER_QUIT_ERR => "User requested termination".to_owned(),
        CT_REMOVE_STEREO_ERR => "REMOVE_STEREO_ERR".to_owned(),
        CT_CALC_STEREO_ERR => "CALC_STEREO_ERR".to_owned(),
        CT_STEREO_CANON_ERR => "STEREO_CANON_ERR".to_owned(),
        CT_CANON_ERR => "CANON_ERR".to_owned(),
        CT_WRONG_FORMULA => "Wrong or missing chemical formula".to_owned(),
        CT_UNKNOWN_ERR => "UNKNOWN_ERR".to_owned(),
        BNS_RADICAL_ERR => "Cannot process free radical center".to_owned(),
        BNS_ALTBOND_ERR => "Cannot process aromatic bonds".to_owned(),
        BNS_TIMEOUT => "Structure normalization timeout".to_owned(),
        code if code > CT_UNKNOWN_ERR => format!("No description({code})"),
        code => format!("UNKNOWN_ERR({})", CT_UNKNOWN_ERR.wrapping_sub(code)),
    }
}

pub(crate) fn already_have_this_message(prev_messages: &[i8], new_message: &[i8]) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.c:160 already_have_this_message
    // INCHI✔️✔️: int already_have_this_message( char *prev_messages, const char *new_message )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int have = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     char *p = strstr( prev_messages, new_message );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (p)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         have = ( p == prev_messages || (*( p - 1 ) == ' ' && ( *( p - 2 ) == ';' || *( p - 2 ) == ':' )) ); /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         if (have)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             int len_prev = (int) strlen( prev_messages );
    // INCHI✔️✔️:             int len = (int) strlen( new_message );
    // INCHI✔️✔️:             have = ( p + len == prev_messages + len_prev || (p[len] == ';' && p[len + 1] == ' ') || (p[len - 1] == ':' && p[len] == ' ') ); /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return have;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: already_have_this_message

    let previous_length = c_string_length(prev_messages)?;
    let message_length = c_string_length(new_message)?;
    let previous = &prev_messages[..previous_length];
    let message = &new_message[..message_length];
    let Some(position) = find_bytes(previous, message) else {
        return Ok(0);
    };

    let starts_at_message_boundary = position == 0
        || (position >= 2
            && previous[position - 1] == b' ' as i8
            && matches!(previous[position - 2] as u8, b';' | b':'));
    if !starts_at_message_boundary {
        return Ok(0);
    }

    let end = position + message_length;
    let ends_at_message_boundary = end == previous_length
        || (end + 1 < previous_length && previous[end] == b';' as i8 && previous[end + 1] == b' ' as i8)
        || (message_length > 0
            && message[message_length - 1] == b':' as i8
            && end < previous_length
            && previous[end] == b' ' as i8);
    Ok(i32::from(ends_at_message_boundary))
}

pub(crate) fn AddErrorMessage(
    all_messages: Option<&mut [i8]>,
    new_message: Option<&[i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.c:106 AddErrorMessage
    // INCHI✔️✔️: int AddErrorMessage( char *all_messages, const char *new_message )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int len_all, len;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (!all_messages)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (!new_message)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (!new_message[0])
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (already_have_this_message( all_messages, new_message ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     len_all = (int) strlen( all_messages );
    // INCHI✔️✔️:     len = (int) strlen( new_message );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (len_all + len + 2 * ( len_all > 0 ) < STR_ERR_LEN)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* enough room... add message and return */
    // INCHI✔️✔️:         if (len_all > 0)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (all_messages[len_all - 1] != ':')
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 strcat(all_messages, ";");
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             strcat(all_messages, " ");
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         strcat(all_messages, new_message);
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*  not enough room... add no-room marker if not yet added */
    // INCHI✔️✔️:     if (strstr( all_messages, "..." ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (len_all + 3 < STR_ERR_LEN)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         strcat(all_messages, "...");
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: AddErrorMessage

    let Some(all_messages) = all_messages else {
        return Ok(0);
    };
    let Some(new_message) = new_message else {
        return Ok(0);
    };
    let message_length = c_string_length(new_message)?;
    if message_length == 0 {
        return Ok(0);
    }
    if all_messages.len() < STR_ERR_LEN as usize {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    if already_have_this_message(all_messages, new_message)? != 0 {
        return Ok(1);
    }

    let all_length = c_string_length(all_messages)?;
    if all_length + message_length + 2 * usize::from(all_length > 0) < STR_ERR_LEN as usize {
        let mut output = all_length;
        if all_length > 0 {
            if all_messages[all_length - 1] != b':' as i8 {
                all_messages[output] = b';' as i8;
                output += 1;
            }
            all_messages[output] = b' ' as i8;
            output += 1;
        }
        all_messages[output..output + message_length].copy_from_slice(&new_message[..message_length]);
        all_messages[output + message_length] = 0;
        return Ok(1);
    }

    if find_bytes(&all_messages[..all_length], &[b'.' as i8; 3]).is_some() {
        return Ok(0);
    }
    if all_length + 3 < STR_ERR_LEN as usize {
        all_messages[all_length..all_length + 3].fill(b'.' as i8);
        all_messages[all_length + 3] = 0;
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
        assert_eq!(AddErrorMessage(Some(&mut messages), Some(&c_bytes(""))), Ok(0));

        assert_eq!(AddErrorMessage(Some(&mut messages), Some(&message)), Ok(1));
        assert_eq!(buffer_text(&messages), "first");
        assert_eq!(AddErrorMessage(Some(&mut messages), Some(&message)), Ok(1));
        assert_eq!(buffer_text(&messages), "first");
        assert_eq!(AddErrorMessage(Some(&mut messages), Some(&c_bytes("second"))), Ok(1));
        assert_eq!(buffer_text(&messages), "first; second");

        let mut colon = error_buffer("header:");
        assert_eq!(AddErrorMessage(Some(&mut colon), Some(&c_bytes("detail"))), Ok(1));
        assert_eq!(buffer_text(&colon), "header: detail");

        let mut no_room = error_buffer(&"x".repeat(253));
        assert_eq!(AddErrorMessage(Some(&mut no_room), Some(&c_bytes("more"))), Ok(0));
        assert_eq!(buffer_text(&no_room), "x".repeat(253));

        let mut marker_room = error_buffer(&"x".repeat(252));
        assert_eq!(AddErrorMessage(Some(&mut marker_room), Some(&c_bytes("more"))), Ok(0));
        assert_eq!(buffer_text(&marker_room), format!("{}...", "x".repeat(252)));
        assert_eq!(AddErrorMessage(Some(&mut marker_room), Some(&c_bytes("again"))), Ok(0));
        assert_eq!(buffer_text(&marker_room), format!("{}...", "x".repeat(252)));

        let mut short = [0_i8; 8];
        assert_eq!(
            AddErrorMessage(Some(&mut short), Some(&message)),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }
}
