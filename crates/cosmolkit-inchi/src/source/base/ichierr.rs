use crate::source_types::{STR_ERR_LEN, SourceHeapError};

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
    haystack
        .windows(needle.len())
        .position(|window| window == needle)
}

pub(crate) fn already_have_this_message(
    prev_messages: &[i8],
    new_message: &[i8],
) -> Result<i32, SourceHeapError> {
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
        || (end + 1 < previous_length
            && previous[end] == b';' as i8
            && previous[end + 1] == b' ' as i8)
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
        all_messages[output..output + message_length]
            .copy_from_slice(&new_message[..message_length]);
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
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }
}
