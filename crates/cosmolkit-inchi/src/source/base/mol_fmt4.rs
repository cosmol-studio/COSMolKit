use crate::source::base::ichi_bns::nBondsValenceInpAt;
use crate::source::base::ichi_io::{inchi_fgetsLf, inchi_ios_print_nodisplay};
use crate::source::base::ichierr::AddErrorMessage;
use crate::source::base::util::{
    dotify_non_printable_chars, get_atomic_mass_from_elnum, inchi_calloc, inchi_free,
    inchi_memicmp, is_in_the_ilist, mystrncpy, needed_unusual_el_valence, normalize_string,
    remove_trailing_spaces,
};
use crate::source_types::local_mol_fmt4::{
    SD_FMT_END_OF_DATA_BLOCK, SD_FMT_END_OF_DATA_ITEM, SDF_DATA_HEADER, SDF_DATA_HEADER_CAS,
    SDF_DATA_HEADER_COMMENT, SDF_DATA_HEADER_NAME, SDF_DATA_HEADER_USER, SDF_DATA_LINE,
    SDF_EMPTY_LINE, SDF_START,
};
use crate::source_types::{
    _IS_ERROR, FILE, INCHI_IOSTREAM, INT_ARRAY, MAX_SDF_VALUE, MOL_FMT_INPLINELEN,
    MOL_FMT_M_CONN_EU, MOL_FMT_M_CONN_HH, MOL_FMT_M_CONN_HT, MOL_FMT_M_SST_ALT, MOL_FMT_M_SST_BLK,
    MOL_FMT_M_SST_RAN, MOL_FMT_MAXLINELEN, MOL_FMT_SGROUP, MOL_FMT_SGROUPS, NUM_LISTS,
    OAD_PolymerUnit, ORIG_ATOM_DATA, RADICAL_DOUBLET, RADICAL_SINGLET, RADICAL_TRIPLET,
    SD_FMT_END_OF_DATA, SourceConstPointer, SourceFormatArgument, SourceHeap, SourceHeapError,
    SourceMutPointer, SourceVaList, inp_ATOM,
};

fn mol_fmt4_add_ascii_message(
    target: Option<&mut [i8]>,
    message: &[u8],
) -> Result<i32, SourceHeapError> {
    let terminated = message
        .iter()
        .map(|byte| *byte as i8)
        .chain(std::iter::once(0))
        .collect::<Vec<_>>();
    AddErrorMessage(target, Some(&terminated))
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn SDFileIdentifyLabel(
    heap: &mut SourceHeap,
    inp_line: SourceConstPointer<i8>,
    sdf_label: SourceConstPointer<i8>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:346 SDFileIdentifyLabel
    // INCHI✔️❌: complete source frame follows verbatim; SourceHeap stack-array and literal-pointer modeling adds allocation-map overhead.
    /*
int SDFileIdentifyLabel(char *inp_line, const char *pSdfLabel)
{
    char line[MOL_FMT_MAXLINELEN];
    char *p, *q;
    int i, j, len, tmp1 = 0, tmp2 = 0, cnd = 0;

    if ((p = strchr(inp_line, '<')) &&
        (q = strchr(p, '>')) &&
        (len = q - p - 1) > 0 && len < (int)sizeof(line))
    {
        memcpy(line, p + 1, len);
        line[len] = '\0';

        for (i = 0; (i < len) && (cnd == 0); i++)
        {
            if (isspace(UCINT line[i]))
            {
                tmp1 = i;
                cnd = 1;
            }
        }

        cnd = 0;

        for (j = len - 1; (j >= tmp1) && (cnd == 0); j--)
        {
            if (isspace(UCINT line[j]))
            {
                tmp2 = j;
                cnd = 1;
            }
        }

        len = tmp2 - tmp1 + 1;
        p = line + tmp1;

        if (pSdfLabel && pSdfLabel[0] && len == (int)strlen(pSdfLabel) && !inchi_memicmp(p, pSdfLabel, len))
        {
            return SDF_DATA_HEADER_USER;
        }

        if (len == sizeof(sdf_data_hdr_name) - 1 && !inchi_memicmp(p, sdf_data_hdr_name, len))
        {
            return SDF_DATA_HEADER_NAME;
        }

        if (len == sizeof(sdf_data_hdr_comm) - 1 && !inchi_memicmp(p, sdf_data_hdr_comm, len))
        {
            return SDF_DATA_HEADER_COMMENT;
        }

        if (!inchi_memicmp(p, "CAS", 3))
        {
            return SDF_DATA_HEADER_CAS;
        }
    }

    return SDF_DATA_HEADER;
}
    */
    // END INCHI C FUNCTION: SDFileIdentifyLabel
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: SDFileIdentifyLabel
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; UCINT casts to unsigned char and the selected C locale has ASCII isspace behavior.
    // INCHI✔️❌: strchr, memcpy, strlen, and completed inchi_memicmp are active; the source tmp1/tmp2 initialization and fixed three-byte CAS comparison are preserved literally.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: SDFileIdentifyLabel

    let line_capacity = usize::try_from(MOL_FMT_MAXLINELEN)
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let line = heap.allocate_model_storage(vec![0_i8; line_capacity])?;
    let name = heap.allocate_model_storage(b"NAME\0".iter().map(|byte| *byte as i8).collect())?;
    let comment = heap.allocate_model_storage(b"COMMENT\0".iter().map(|byte| *byte as i8).collect())?;
    let cas = heap.allocate_model_storage(b"CAS\0".iter().map(|byte| *byte as i8).collect())?;
    let result = (|| -> Result<i32, SourceHeapError> {
        let input = heap.slice(inp_line)?;
        let input_length = input.iter().position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        let input = &input[..input_length];
        let Some(open) = input.iter().position(|byte| *byte == b'<' as i8) else {
            return Ok(SDF_DATA_HEADER as i32);
        };
        let Some(relative_close) = input[open..].iter().position(|byte| *byte == b'>' as i8) else {
            return Ok(SDF_DATA_HEADER as i32);
        };
        let close = open.wrapping_add(relative_close);
        let original_length = close.wrapping_sub(open).wrapping_sub(1);
        if original_length == 0 || original_length >= line_capacity {
            return Ok(SDF_DATA_HEADER as i32);
        }
        let copied_input = input[open + 1..close].to_vec();
        heap.slice_mut(line)?[..original_length]
            .copy_from_slice(&copied_input);
        heap.slice_mut(line)?[original_length] = 0;

        let copied = heap.slice(line.as_const())?;
        let mut tmp1 = 0_usize;
        if let Some(index) = copied[..original_length].iter().position(|byte| {
            matches!(*byte as u8, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
        }) {
            tmp1 = index;
        }
        let mut tmp2 = 0_usize;
        if let Some(index) = copied[tmp1..original_length].iter().rposition(|byte| {
            matches!(*byte as u8, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
        }) {
            tmp2 = tmp1.wrapping_add(index);
        }
        let length = tmp2.wrapping_sub(tmp1).wrapping_add(1);
        let candidate = line.as_const().offset(i64::try_from(tmp1)
            .map_err(|_| SourceHeapError::PointerOffsetOverflow)?)?;

        if !sdf_label.is_null() {
            let label = heap.slice(sdf_label)?;
            let label_length = label.iter().position(|byte| *byte == 0)
                .ok_or(SourceHeapError::MissingNulTerminator)?;
            if label.first().copied().unwrap_or(0) != 0
                && length == label_length
                && inchi_memicmp(heap, candidate, sdf_label, length as u64)? == 0
            {
                return Ok(SDF_DATA_HEADER_USER as i32);
            }
        }
        if length == 4 && inchi_memicmp(heap, candidate, name.as_const(), length as u64)? == 0 {
            return Ok(SDF_DATA_HEADER_NAME as i32);
        }
        if length == 7 && inchi_memicmp(heap, candidate, comment.as_const(), length as u64)? == 0 {
            return Ok(SDF_DATA_HEADER_COMMENT as i32);
        }
        if inchi_memicmp(heap, candidate, cas.as_const(), 3)? == 0 {
            return Ok(SDF_DATA_HEADER_CAS as i32);
        }
        Ok(SDF_DATA_HEADER as i32)
    })();
    let cleanup = [line, name, comment, cas]
        .into_iter()
        .try_for_each(|pointer| heap.free(pointer));
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn SDFileExtractCASNo(
    heap: &mut SourceHeap,
    line: SourceMutPointer<i8>,
) -> Result<u64, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:407 SDFileExtractCASNo
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap allocation-map access adds constant overhead.
    /*
unsigned long SDFileExtractCASNo(char *line)
{
    int i, j;

    i = line[0] == '-' ? 1 : 0;

    for (j = i; line[i]; i++)
    {
        if (isdigit(UCINT line[i]))
        {
            line[j++] = line[i];
        }
        else if (line[i] != '-')
        {
            break;
        }
    }

    line[j] = '\0';

    return strtoul(line, NULL, 10);
}
    */
    // END INCHI C FUNCTION: SDFileExtractCASNo
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: SDFileExtractCASNo
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; unsigned long is 64 bits and UCINT casts to unsigned char.
    // INCHI✔️❌: The selected C locale recognizes only ASCII decimal digits; libc strtoul base 10 returns ULONG_MAX on range error and applies a leading minus by unsigned negation when the magnitude is representable.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: SDFileExtractCASNo

    let bytes = heap.slice_mut(line)?;
    let nul = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let mut i = usize::from(bytes.first().copied() == Some(b'-' as i8));
    let mut j = i;
    while i < nul {
        let byte = bytes[i] as u8;
        if byte.is_ascii_digit() {
            bytes[j] = bytes[i];
            j = j.wrapping_add(1);
        } else if byte != b'-' {
            break;
        }
        i = i.wrapping_add(1);
    }
    bytes[j] = 0;

    let negative = j != 0 && bytes[0] == b'-' as i8;
    let first_digit = usize::from(negative);
    let mut value = 0_u64;
    let mut overflow = false;
    for byte in &bytes[first_digit..j] {
        let digit = (*byte as u8).wrapping_sub(b'0') as u64;
        if !overflow {
            match value.checked_mul(10).and_then(|value| value.checked_add(digit)) {
                Some(parsed) => value = parsed,
                None => overflow = true,
            }
        }
    }
    if overflow {
        Ok(u64::MAX)
    } else if negative {
        Ok(0_u64.wrapping_sub(value))
    } else {
        Ok(value)
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn SDFileSkipExtraData(
    heap: &mut SourceHeap,
    mut inp_file: Option<&mut INCHI_IOSTREAM>,
    mut CAS_num: Option<&mut u64>,
    comment: SourceMutPointer<i8>,
    lcomment: i32,
    name: SourceMutPointer<i8>,
    lname: i32,
    prev_err: i32,
    pSdfLabel: SourceConstPointer<i8>,
    pSdfValue: SourceMutPointer<i8>,
    mut pStrErr: Option<&mut [i8]>,
    bNoWarnings: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:161 SDFileSkipExtraData
    // INCHI✔️❌: complete source frame follows verbatim; modeled stack buffers and checked SourceHeap pointer access add allocation-map overhead.
    /*
int SDFileSkipExtraData(INCHI_IOSTREAM *inp_file,
                        unsigned long *CAS_num,
                        char *comment,
                        int lcomment,
                        char *name,
                        int lname,
                        int prev_err,
                        const char *pSdfLabel,
                        char *pSdfValue,
                        char *pStrErr,
                        int bNoWarnings)
{
    char *p = NULL;
    char line[MOL_FMT_INPLINELEN];
    const int line_len = sizeof( line );
    int n_blank_lines = 0, n_lines = 0;
    int current_state = SDF_START;
    int err = 0;
    int wait_for_CAS = 0;
    int CAS_num_is_user = 0;
    int wait_for_name = name && lname > 0 && !name[0];
    int wait_for_comment = comment && lcomment > 0 && !comment[0];
    int wait_for_user = pSdfLabel && pSdfLabel[0] && pSdfValue;

    if (CAS_num != NULL)
    {
        wait_for_CAS = 1;
        *CAS_num = 0LU;
        CAS_num_is_user = (wait_for_user && !inchi_memicmp(pSdfLabel, "CAS", 3));
    }

    while (!err &&
           current_state != SD_FMT_END_OF_DATA_BLOCK &&
           NULL != (p = inchi_fgetsLf(line, line_len, inp_file)))
    {
        if (!n_lines && !memcmp(line, "M  END", 6))
        {
            /*  allow subtle errors */
            continue;
        }

        n_lines++;
        remove_trailing_spaces(line);

        if (line[MOL_FMT_MAXLINELEN])
        {
            if (current_state != SDF_DATA_HEADER &&
                current_state != SDF_DATA_LINE &&
                current_state != SDF_DATA_HEADER_NAME &&
                current_state != SDF_DATA_HEADER_USER &&
                current_state != SDF_DATA_HEADER_COMMENT)
            {
                line[MOL_FMT_MAXLINELEN] = '\0';
                if (!prev_err)
                {
                    TREAT_ERR(err, 0, "Too long SData line truncated");
                }
            }
            else
            {
                /* allow long lines in SDF data. 9-29-00 DCh */
                line[MOL_FMT_MAXLINELEN] = '\0';
            }
        }

        n_blank_lines += (*line == '\0');

        switch (current_state)
        {
            case SDF_START:
            case SD_FMT_END_OF_DATA_ITEM:
            case SDF_EMPTY_LINE: /* Added 9-25-97 DCh */

                if (!strcmp(line, SD_FMT_END_OF_DATA))
                {
                    current_state = SD_FMT_END_OF_DATA_BLOCK;
                }
                else if ('>' == *line)
                {
                    current_state = (wait_for_name || wait_for_comment || wait_for_CAS || wait_for_user) ? SDFileIdentifyLabel( line, pSdfLabel ) : SDF_DATA_HEADER;
                }
                else if (*line == '\0')
                {
                    /* Added 9-25-97 DCh */
                    /* Relax the strictness: Allow more than 1 empty line. */
                    current_state = SDF_EMPTY_LINE;
                }
                else if (!prev_err)
                {
                    TREAT_ERR(err, 3, "Unexpected SData header line:"); /* djb-rwth: addressing coverity ID #499557 -- TREAT_ERR properly used */
                    dotify_non_printable_chars(line);
                    AddErrorMessage(pStrErr, line);
                    /* unexpected contents of data header line */
                }
                else
                {
                    err = 3;
                }
                break;

            case SDF_DATA_HEADER_NAME:

                if (wait_for_name && 0 < normalize_string(line))
                {
                    wait_for_name = 0;
                    mystrncpy(name, line, lname);
                }
                goto got_data_line;

            case SDF_DATA_HEADER_COMMENT:

                if (wait_for_comment && 0 < normalize_string(line))
                {
                    wait_for_comment = 0;
                    mystrncpy(comment, line, lcomment);
                }
                goto got_data_line;

            case SDF_DATA_HEADER_USER:

                if (wait_for_user && 0 < normalize_string(line))
                {
                    wait_for_user = 0;
                    mystrncpy(pSdfValue, line, MAX_SDF_VALUE + 1);

                    if (CAS_num_is_user && wait_for_CAS)
                    {
                        *CAS_num = SDFileExtractCASNo(line);
                        wait_for_CAS = (0LU == *CAS_num);
                    }
                }
                goto got_data_line;

            case SDF_DATA_HEADER_CAS:

                if (wait_for_CAS && 0 < normalize_string(line))
                {
                    *CAS_num = SDFileExtractCASNo(line);
                    wait_for_CAS = (0LU == *CAS_num);
                }
                goto got_data_line;

            case SDF_DATA_HEADER:
            case SDF_DATA_LINE:

            got_data_line:
                current_state = *line ? SDF_DATA_LINE : SD_FMT_END_OF_DATA_ITEM;
                break;
        }
    }

    if (!err && SD_FMT_END_OF_DATA_BLOCK != current_state && NULL == p)
    {
        ; /* err = 4; */ /* unexpected end of file: missing $$$$ */
    }

    else if (err && (n_blank_lines == n_lines && *line == '\0'))
    {
        /* empty lines -- do not know when this can happen */
        err = 5;
    }

    if (err && err != 5 && current_state != SD_FMT_END_OF_DATA_BLOCK && p)
    {
        /*  bypass up to $$$$ */
        while ((p = inchi_fgetsLf(line, line_len, inp_file)) &&
               memcmp(line, SD_FMT_END_OF_DATA, 4))
        {
            ;
        }
        if (p)
        {
            /*  arrived to $$$$; non-fatal */
            err = 9;
            if (!bNoWarnings)
            {
                WarningMessage(pStrErr, "Bypassing to next structure");
            }
        }
    }

    return err;
}
    */
    // END INCHI C FUNCTION: SDFileSkipExtraData
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: SDFileSkipExtraData
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; unsigned long is u64, line_len is 204, and MOL_FMT_MAXLINELEN is 200.
    // INCHI✔️❌: WarningMessage expands to completed AddErrorMessage; TREAT_ERR preserves a nonzero error, assigns only a nonzero new error, and always calls AddErrorMessage.
    // INCHI✔️❌: All direct callees are completed source ports: inchi_fgetsLf, remove_trailing_spaces, SDFileIdentifyLabel, normalize_string, mystrncpy, SDFileExtractCASNo, dotify_non_printable_chars, inchi_memicmp, and AddErrorMessage.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: SDFileSkipExtraData

    let line = heap.allocate_model_storage(vec![0_i8; MOL_FMT_INPLINELEN as usize])?;
    let cas_literal = heap.allocate_model_storage(b"CAS\0".iter().map(|byte| *byte as i8).collect())?;
    let result = (|| -> Result<i32, SourceHeapError> {
        let mut p = SourceMutPointer::null();
        let mut n_blank_lines = 0_i32;
        let mut n_lines = 0_i32;
        let mut current_state = SDF_START as i32;
        let mut err = 0_i32;
        let mut wait_for_CAS = false;
        let mut CAS_num_is_user = false;
        let mut wait_for_name = !name.is_null()
            && lname > 0
            && heap.slice(name.as_const())?.first().copied() == Some(0);
        let mut wait_for_comment = !comment.is_null()
            && lcomment > 0
            && heap.slice(comment.as_const())?.first().copied() == Some(0);
        let mut wait_for_user = !pSdfLabel.is_null()
            && heap.slice(pSdfLabel)?.first().copied().unwrap_or(0) != 0
            && !pSdfValue.is_null();

        if let Some(value) = CAS_num.as_deref_mut() {
            wait_for_CAS = true;
            *value = 0;
            CAS_num_is_user = wait_for_user
                && inchi_memicmp(heap, pSdfLabel, cas_literal.as_const(), 3)? == 0;
        }

        while err == 0 && current_state != SD_FMT_END_OF_DATA_BLOCK as i32 {
            p = inchi_fgetsLf(
                heap,
                line,
                MOL_FMT_INPLINELEN as i32,
                inp_file.as_deref_mut(),
            )?;
            if p.is_null() {
                break;
            }
            if n_lines == 0
                && heap.slice(line.as_const())?[..6]
                    == [b'M' as i8, b' ' as i8, b' ' as i8, b'E' as i8, b'N' as i8, b'D' as i8]
            {
                continue;
            }

            n_lines = n_lines.wrapping_add(1);
            remove_trailing_spaces(heap, line)?;

            if heap.slice(line.as_const())?[MOL_FMT_MAXLINELEN as usize] != 0 {
                if current_state != SDF_DATA_HEADER as i32
                    && current_state != SDF_DATA_LINE as i32
                    && current_state != SDF_DATA_HEADER_NAME as i32
                    && current_state != SDF_DATA_HEADER_USER as i32
                    && current_state != SDF_DATA_HEADER_COMMENT as i32
                {
                    heap.slice_mut(line)?[MOL_FMT_MAXLINELEN as usize] = 0;
                    if prev_err == 0 {
                        let _ = mol_fmt4_add_ascii_message(
                            pStrErr.as_deref_mut(),
                            b"Too long SData line truncated",
                        )?;
                    }
                } else {
                    heap.slice_mut(line)?[MOL_FMT_MAXLINELEN as usize] = 0;
                }
            }

            let line_is_empty = heap.slice(line.as_const())?.first().copied() == Some(0);
            n_blank_lines = n_blank_lines.wrapping_add(i32::from(line_is_empty));

            match current_state {
                state if state == SDF_START as i32
                    || state == SD_FMT_END_OF_DATA_ITEM as i32
                    || state == SDF_EMPTY_LINE as i32 =>
                {
                    let bytes = heap.slice(line.as_const())?;
                    let length = bytes
                        .iter()
                        .position(|byte| *byte == 0)
                        .ok_or(SourceHeapError::MissingNulTerminator)?;
                    let is_end = length == SD_FMT_END_OF_DATA.len() - 1
                        && bytes[..length]
                            .iter()
                            .zip(&SD_FMT_END_OF_DATA[..length])
                            .all(|(left, right)| *left as u8 == *right);
                    if is_end {
                        current_state = SD_FMT_END_OF_DATA_BLOCK as i32;
                    } else if bytes.first().copied() == Some(b'>' as i8) {
                        current_state = if wait_for_name
                            || wait_for_comment
                            || wait_for_CAS
                            || wait_for_user
                        {
                            SDFileIdentifyLabel(heap, line.as_const(), pSdfLabel)?
                        } else {
                            SDF_DATA_HEADER as i32
                        };
                    } else if line_is_empty {
                        current_state = SDF_EMPTY_LINE as i32;
                    } else if prev_err == 0 {
                        if err == 0 {
                            err = 3;
                        }
                        let _ = mol_fmt4_add_ascii_message(
                            pStrErr.as_deref_mut(),
                            b"Unexpected SData header line:",
                        )?;
                        let _ = dotify_non_printable_chars(heap, line)?;
                        let message = heap.slice(line.as_const())?;
                        let _ = AddErrorMessage(pStrErr.as_deref_mut(), Some(message))?;
                    } else {
                        err = 3;
                    }
                }
                state if state == SDF_DATA_HEADER_NAME as i32 => {
                    if wait_for_name && normalize_string(heap, line)? > 0 {
                        wait_for_name = false;
                        let _ = mystrncpy(heap, name, line.as_const(), lname as u32)?;
                    }
                    current_state = if heap.slice(line.as_const())?.first().copied() != Some(0) {
                        SDF_DATA_LINE as i32
                    } else {
                        SD_FMT_END_OF_DATA_ITEM as i32
                    };
                }
                state if state == SDF_DATA_HEADER_COMMENT as i32 => {
                    if wait_for_comment && normalize_string(heap, line)? > 0 {
                        wait_for_comment = false;
                        let _ = mystrncpy(heap, comment, line.as_const(), lcomment as u32)?;
                    }
                    current_state = if heap.slice(line.as_const())?.first().copied() != Some(0) {
                        SDF_DATA_LINE as i32
                    } else {
                        SD_FMT_END_OF_DATA_ITEM as i32
                    };
                }
                state if state == SDF_DATA_HEADER_USER as i32 => {
                    if wait_for_user && normalize_string(heap, line)? > 0 {
                        wait_for_user = false;
                        let _ = mystrncpy(
                            heap,
                            pSdfValue,
                            line.as_const(),
                            MAX_SDF_VALUE + 1,
                        )?;
                        if CAS_num_is_user && wait_for_CAS {
                            if let Some(value) = CAS_num.as_deref_mut() {
                                *value = SDFileExtractCASNo(heap, line)?;
                                wait_for_CAS = *value == 0;
                            }
                        }
                    }
                    current_state = if heap.slice(line.as_const())?.first().copied() != Some(0) {
                        SDF_DATA_LINE as i32
                    } else {
                        SD_FMT_END_OF_DATA_ITEM as i32
                    };
                }
                state if state == SDF_DATA_HEADER_CAS as i32 => {
                    if wait_for_CAS && normalize_string(heap, line)? > 0 {
                        if let Some(value) = CAS_num.as_deref_mut() {
                            *value = SDFileExtractCASNo(heap, line)?;
                            wait_for_CAS = *value == 0;
                        }
                    }
                    current_state = if heap.slice(line.as_const())?.first().copied() != Some(0) {
                        SDF_DATA_LINE as i32
                    } else {
                        SD_FMT_END_OF_DATA_ITEM as i32
                    };
                }
                state if state == SDF_DATA_HEADER as i32 || state == SDF_DATA_LINE as i32 => {
                    current_state = if heap.slice(line.as_const())?.first().copied() != Some(0) {
                        SDF_DATA_LINE as i32
                    } else {
                        SD_FMT_END_OF_DATA_ITEM as i32
                    };
                }
                _ => {}
            }
        }

        if err == 0 && current_state != SD_FMT_END_OF_DATA_BLOCK as i32 && p.is_null() {
        } else if err != 0
            && n_blank_lines == n_lines
            && heap.slice(line.as_const())?.first().copied() == Some(0)
        {
            err = 5;
        }

        if err != 0
            && err != 5
            && current_state != SD_FMT_END_OF_DATA_BLOCK as i32
            && !p.is_null()
        {
            loop {
                p = inchi_fgetsLf(
                    heap,
                    line,
                    MOL_FMT_INPLINELEN as i32,
                    inp_file.as_deref_mut(),
                )?;
                if p.is_null()
                    || heap.slice(line.as_const())?[..4]
                        == [b'$' as i8, b'$' as i8, b'$' as i8, b'$' as i8]
                {
                    break;
                }
            }
            if !p.is_null() {
                err = 9;
                if bNoWarnings == 0 {
                    let _ = mol_fmt4_add_ascii_message(
                        pStrErr.as_deref_mut(),
                        b"Bypassing to next structure",
                    )?;
                }
            }
        }
        Ok(err)
    })();
    let cleanup = heap.free(line).and_then(|_| heap.free(cas_literal));
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn NumLists_Alloc(
    heap: &mut SourceHeap,
    num_lists: Option<&mut NUM_LISTS>,
    nlists: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:433 NumLists_Alloc
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap allocation-map access adds overhead.
    /*
int NumLists_Alloc(NUM_LISTS *num_lists, int nlists)
{
    if (num_lists)
    {
        if ((num_lists->lists = (int **)inchi_calloc(nlists, sizeof(int *)))) /* djb-rwth: addressing LLVM warning */
        {
            num_lists->increment =
                num_lists->allocated = nlists;
            return 0; /*  ok */
        }
    }

    return -1; /*  error */
}
    */
    // END INCHI C FUNCTION: NumLists_Alloc
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: NumLists_Alloc
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; inchi_calloc resolves to libc calloc and pointer size is eight bytes.
    // INCHI✔️❌: Assignment occurs before the allocation-result test, so failure stores NULL; success does not reset used or free a prior lists pointer.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: NumLists_Alloc

    let Some(num_lists) = num_lists else { return Ok(-1); };
    let replacement = match inchi_calloc::<SourceMutPointer<i32>>(heap, nlists as u64, 8) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed)
        | Err(SourceHeapError::AllocationElementCountOutOfRange)
        | Err(SourceHeapError::AllocationSizeOverflow) => {
            num_lists.lists = SourceMutPointer::null();
            return Ok(-1);
        }
        Err(error) => return Err(error),
    };
    num_lists.lists = replacement;
    num_lists.allocated = nlists;
    num_lists.increment = nlists;
    Ok(0)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn NumLists_ReAlloc(
    heap: &mut SourceHeap,
    num_lists: Option<&mut NUM_LISTS>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:449 NumLists_ReAlloc
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap copy adds overhead.
    /*
int NumLists_ReAlloc(NUM_LISTS *num_lists)
{
    if (num_lists)
    {
        if (num_lists->lists && num_lists->allocated > 0 && num_lists->increment > 0)
        {
            void *p = num_lists->lists;
            if ((num_lists->lists =
                     (int **)inchi_calloc((long long)num_lists->allocated + (long long)num_lists->increment, sizeof(int *)))) /* djb-rwth: cast operators added; addressing LLVM warning */
            {
                memcpy(num_lists->lists, p, num_lists->used * sizeof(num_lists->lists[0]));
                inchi_free(p);
                num_lists->allocated += num_lists->increment;
                return 0; /*  ok */
            }
        }
    }

    return -1; /*  error */
}
    */
    // END INCHI C FUNCTION: NumLists_ReAlloc
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: NumLists_ReAlloc
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; pointer size is eight bytes.
    // INCHI✔️❌: Failed calloc leaves num_lists->lists NULL and does not free the saved old array.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: NumLists_ReAlloc

    let Some(num_lists) = num_lists else { return Ok(-1); };
    if num_lists.lists.is_null() || num_lists.allocated <= 0 || num_lists.increment <= 0 {
        return Ok(-1);
    }
    let old = num_lists.lists;
    let count = i64::from(num_lists.allocated) + i64::from(num_lists.increment);
    let replacement = match inchi_calloc::<SourceMutPointer<i32>>(heap, count as u64, 8) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed)
        | Err(SourceHeapError::AllocationElementCountOutOfRange)
        | Err(SourceHeapError::AllocationSizeOverflow) => {
            num_lists.lists = SourceMutPointer::null();
            return Ok(-1);
        }
        Err(error) => return Err(error),
    };
    num_lists.lists = replacement;
    let used = usize::try_from(num_lists.used).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let copied = heap
        .slice(old.as_const())?
        .get(..used)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    heap.slice_mut(replacement)?[..used].copy_from_slice(&copied);
    inchi_free(heap, old)?;
    num_lists.allocated = num_lists.allocated.wrapping_add(num_lists.increment);
    Ok(0)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn NumLists_Append(
    heap: &mut SourceHeap,
    num_lists: Option<&mut NUM_LISTS>,
    list: SourceMutPointer<i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:471 NumLists_Append
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap indexing adds overhead.
    /*
int NumLists_Append(NUM_LISTS *num_lists, int *list)
{
    if (num_lists)
    {
        if (num_lists->used + 1 > num_lists->allocated)
        {
            /* need to expand buffer */
            if (NumLists_ReAlloc(num_lists))
            {
                return -1; /*  error */
            }
        }
        num_lists->lists[num_lists->used++] = list;
        return 0;
    }

    return -1;
}
    */
    // END INCHI C FUNCTION: NumLists_Append
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: NumLists_Append
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; NumLists_ReAlloc is the completed source callee.
    // INCHI✔️❌: A NULL list value is stored; only a NULL NUM_LISTS descriptor is rejected.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: NumLists_Append

    let Some(num_lists) = num_lists else { return Ok(-1); };
    if num_lists.used.wrapping_add(1) > num_lists.allocated
        && NumLists_ReAlloc(heap, Some(num_lists))? != 0
    {
        return Ok(-1);
    }
    let used = usize::try_from(num_lists.used).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    *heap
        .slice_mut(num_lists.lists)?
        .get_mut(used)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = list;
    num_lists.used = num_lists.used.wrapping_add(1);
    Ok(0)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn NumLists_Free(
    heap: &mut SourceHeap,
    num_lists: Option<&mut NUM_LISTS>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:491 NumLists_Free
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap frees add allocation-map overhead.
    /*
void NumLists_Free(NUM_LISTS *num_lists)
{
    if (num_lists)
    {
        int i;
        for (i = 0; i < num_lists->used; i++)
            inchi_free(num_lists->lists[i]); /* djb-rwth: unresolved issue -- revision required? -- false positive as this function just does the clean-up job */
        inchi_free(num_lists->lists);
        memset(num_lists, 0, sizeof(*num_lists)); /* djb-rwth: memset_s C11/Annex K variant? */
    }
}
    */
    // END INCHI C FUNCTION: NumLists_Free
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: NumLists_Free
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; inchi_free(X) checks X before free(X).
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: NumLists_Free

    let Some(num_lists) = num_lists else {
        return Ok(());
    };
    let used = num_lists.used.max(0) as usize;
    if used != 0 {
        let pointers = heap
            .slice(num_lists.lists.as_const())?
            .get(..used)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec();
        for pointer in pointers {
            inchi_free(heap, pointer)?;
        }
    }
    inchi_free(heap, num_lists.lists)?;
    *num_lists = NUM_LISTS::default();
    Ok(())
}

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

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolFmtSgroups_Alloc(
    heap: &mut SourceHeap,
    sgroups: Option<&mut MOL_FMT_SGROUPS>,
    nsgroups: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:675 MolFmtSgroups_Alloc
    // INCHI✔️❌: complete source frame follows verbatim; checked source-heap allocation adds overhead.
    /*
int MolFmtSgroups_Alloc(MOL_FMT_SGROUPS *sgroups, int nsgroups)
{
    if (sgroups)
    {
        if (NULL != (sgroups->group = (MOL_FMT_SGROUP **)inchi_calloc(nsgroups, sizeof(MOL_FMT_SGROUP *))))
        {
            /* ITRACE_( "\nAllocated sgroups->group at %-p \n", sgroups->group ); */
            sgroups->increment = sgroups->allocated = nsgroups;
            return 0;
        }
    }

    return -1;
}
    */
    // END INCHI C FUNCTION: MolFmtSgroups_Alloc
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolFmtSgroups_Alloc
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; inchi_calloc is calloc and a pointer is 8 bytes.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolFmtSgroups_Alloc

    let Some(sgroups) = sgroups else {
        return Ok(-1);
    };
    match inchi_calloc::<SourceMutPointer<MOL_FMT_SGROUP>>(heap, nsgroups as u64, 8) {
        Ok(group) => {
            sgroups.group = group;
            sgroups.allocated = nsgroups;
            sgroups.increment = nsgroups;
            Ok(0)
        }
        Err(SourceHeapError::AllocationFailed)
        | Err(SourceHeapError::AllocationElementCountOutOfRange)
        | Err(SourceHeapError::AllocationSizeOverflow) => {
            sgroups.group = SourceMutPointer::null();
            Ok(-1)
        }
        Err(error) => Err(error),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolFmtSgroups_ReAlloc(
    heap: &mut SourceHeap,
    sgroups: Option<&mut MOL_FMT_SGROUPS>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:693 MolFmtSgroups_ReAlloc
    // INCHI✔️❌: complete source frame follows verbatim; checked source-heap copying adds overhead.
    /*
int MolFmtSgroups_ReAlloc(MOL_FMT_SGROUPS *sgroups)
{
    if (sgroups)
    {
        if (sgroups->group && sgroups->allocated > 0 && sgroups->increment > 0)
        {
            void *p = sgroups->group;
            if ((sgroups->group = (MOL_FMT_SGROUP **)inchi_calloc( (long long)sgroups->allocated + (long long)sgroups->increment,
                sizeof( sgroups->group[0] ) ))) /* djb-rwth: cast operators added; addressing LLVM warning */
            {
                memcpy(sgroups->group, p, sgroups->used * sizeof(sgroups->group[0]));
                inchi_free(p);
                sgroups->allocated += sgroups->increment;
                return 0; /*  ok */
            }
        }
    }

    return -1;
}
    */
    // END INCHI C FUNCTION: MolFmtSgroups_ReAlloc
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolFmtSgroups_ReAlloc
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; inchi_calloc is calloc and pointer size is 8.
    // INCHI✔️❌: On calloc failure the assignment leaves sgroups->group NULL and the saved old allocation is not freed.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolFmtSgroups_ReAlloc

    let Some(sgroups) = sgroups else {
        return Ok(-1);
    };
    if sgroups.group.is_null() || sgroups.allocated <= 0 || sgroups.increment <= 0 {
        return Ok(-1);
    }

    let old_group = sgroups.group;
    let count = i64::from(sgroups.allocated) + i64::from(sgroups.increment);
    let replacement = match inchi_calloc::<SourceMutPointer<MOL_FMT_SGROUP>>(
        heap,
        count as u64,
        8,
    ) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed)
        | Err(SourceHeapError::AllocationElementCountOutOfRange)
        | Err(SourceHeapError::AllocationSizeOverflow) => {
            sgroups.group = SourceMutPointer::null();
            return Ok(-1);
        }
        Err(error) => return Err(error),
    };
    sgroups.group = replacement;

    let used = usize::try_from(sgroups.used).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let copied = heap
        .slice(old_group.as_const())?
        .get(..used)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    heap.slice_mut(replacement)?
        .get_mut(..used)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .copy_from_slice(&copied);
    inchi_free(heap, old_group)?;
    sgroups.allocated = sgroups.allocated.wrapping_add(sgroups.increment);
    Ok(0)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolFmtSgroups_Append(
    heap: &mut SourceHeap,
    sgroups: Option<&mut MOL_FMT_SGROUPS>,
    id: i32,
    type_: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:715 MolFmtSgroups_Append
    // INCHI✔️❌: complete source frame follows verbatim; checked source-heap pointer writes add overhead.
    /*
int MolFmtSgroups_Append(MOL_FMT_SGROUPS *sgroups, int id, int type)
{
    if (sgroups)
    {
        /* Make new Sgroup */
        MOL_FMT_SGROUP *sgroup = NULL;
        if (0 != MolFmtSgroup_Create(&sgroup, id, type))
        {
            return -1;
        }
        /* Add new created Sgroup to Sgroups */
        if (sgroups->used + 1 > sgroups->allocated)
        {
            /* expand buffer */
            if (MolFmtSgroups_ReAlloc(sgroups))
            {
                MolFmtSgroup_Free(sgroup); /* djb-rwth: avoiding memory leak */
                return -1; /*  no RAM */
            }
        }
        sgroups->group[sgroups->used++] = sgroup;

        /*
        {
        int num = sgroups->used-1;
        printf("\nCreated/added Sgroup: id=%-d ( num in Sgroups=%-d ) of type=%-d \n", sgroups->group[num]->id, num, sgroups->group[num]->type );
        }
        */

        return 0;
    }

    return -1;
}
    */
    // END INCHI C FUNCTION: MolFmtSgroups_Append
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolFmtSgroups_Append
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; the printf diagnostic is inside a source comment and inactive.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolFmtSgroups_Append

    let Some(sgroups) = sgroups else {
        return Ok(-1);
    };
    let mut sgroup = SourceMutPointer::null();
    if MolFmtSgroup_Create(heap, &mut sgroup, id, type_)? != 0 {
        return Ok(-1);
    }
    if sgroups.used.wrapping_add(1) > sgroups.allocated
        && MolFmtSgroups_ReAlloc(heap, Some(sgroups))? != 0
    {
        MolFmtSgroup_Free(heap, sgroup)?;
        return Ok(-1);
    }

    let index = usize::try_from(sgroups.used).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    *heap
        .slice_mut(sgroups.group)?
        .get_mut(index)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = sgroup;
    sgroups.used = sgroups.used.wrapping_add(1);
    Ok(0)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolFmtSgroups_Free(
    heap: &mut SourceHeap,
    sgroups: Option<&mut MOL_FMT_SGROUPS>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:751 MolFmtSgroups_Free
    // INCHI✔️❌: complete source frame follows verbatim; checked source-heap reads add overhead.
    /*
void MolFmtSgroups_Free(MOL_FMT_SGROUPS *sgroups)
{
    if (sgroups)
    {
        int i;
        for (i = 0; i < sgroups->used; i++)
        {
            MolFmtSgroup_Free(sgroups->group[i]);
        }

        /* ITRACE_( "\nAbout to free sgroups->group at %-p\n", sgroups->group ); */
        inchi_free(sgroups->group);

        memset(sgroups, 0, sizeof(MOL_FMT_SGROUPS)); /* djb-rwth: memset_s C11/Annex K variant? */
    }
}
    */
    // END INCHI C FUNCTION: MolFmtSgroups_Free
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolFmtSgroups_Free
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; ITRACE_ is inactive and inchi_free is the null-checking free macro.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolFmtSgroups_Free

    let Some(sgroups) = sgroups else {
        return Ok(());
    };
    let used = usize::try_from(sgroups.used).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let children = if used == 0 {
        Vec::new()
    } else {
        heap.slice(sgroups.group.as_const())?
            .get(..used)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };
    for child in children {
        MolFmtSgroup_Free(heap, child)?;
    }
    inchi_free(heap, sgroups.group)?;
    *sgroups = MOL_FMT_SGROUPS::default();
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolFmtSgroups_GetIndexBySgroupId(
    heap: &SourceHeap,
    id: i32,
    sgroups: &MOL_FMT_SGROUPS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:769 MolFmtSgroups_GetIndexBySgroupId
    // INCHI✔️❌: complete source frame follows verbatim; checked source-heap allocation-map lookups add overhead.
    /*
int MolFmtSgroups_GetIndexBySgroupId(int id, MOL_FMT_SGROUPS *sgroups)
{
    int i;
    for (i = 0; i < sgroups->used; i++)
    {
        if (sgroups->group[i]->id == id)
        {
            return i;
        }
    }
    return -1;
}
    */
    // END INCHI C FUNCTION: MolFmtSgroups_GetIndexBySgroupId
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolFmtSgroups_GetIndexBySgroupId
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; no conditional branch alters this helper.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolFmtSgroups_GetIndexBySgroupId

    if sgroups.used <= 0 {
        return Ok(-1);
    }
    let used = usize::try_from(sgroups.used).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let groups = heap
        .slice(sgroups.group.as_const())?
        .get(..used)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for (index, group) in groups.iter().enumerate() {
        let value = heap
            .slice(group.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if value.id == id {
            return i32::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds);
        }
    }
    Ok(-1)
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

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolFmtSgroup_Create(
    heap: &mut SourceHeap,
    sgroup: &mut SourceMutPointer<MOL_FMT_SGROUP>,
    id: i32,
    type_: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:634 MolFmtSgroup_Create
    // INCHI✔️❌: complete source frame follows verbatim; checked source-heap object access adds overhead.
    /*
int MolFmtSgroup_Create(MOL_FMT_SGROUP **sgroup, int id, int type)
{
    *sgroup = (MOL_FMT_SGROUP *)inchi_calloc(1, sizeof(MOL_FMT_SGROUP));
    if (*sgroup)
    {
        if (IntArray_Alloc(&((*sgroup)->alist), 8) ||
            IntArray_Alloc(&((*sgroup)->blist), 8))
        {
            MolFmtSgroup_Free(*sgroup);
            return -1;
        }
        (*sgroup)->id = id;
        (*sgroup)->type = type;

        (*sgroup)->subtype = 0;
        (*sgroup)->conn = 0;
        (*sgroup)->label = 0;

        return 0;
    }
    return -1;
}
    */
    // END INCHI C FUNCTION: MolFmtSgroup_Create
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolFmtSgroup_Create
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; inchi_calloc is calloc.
    // INCHI✔️❌: MOL_FMT_SGROUP is 216 bytes and calloc zero-initializes every field before the two IntArray_Alloc calls.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolFmtSgroup_Create

    let allocated = match inchi_calloc::<MOL_FMT_SGROUP>(heap, 1, 216) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed)
        | Err(SourceHeapError::AllocationElementCountOutOfRange)
        | Err(SourceHeapError::AllocationSizeOverflow) => {
            *sgroup = SourceMutPointer::null();
            return Ok(-1);
        }
        Err(error) => return Err(error),
    };
    *sgroup = allocated;

    let mut value = MOL_FMT_SGROUP::default();
    if IntArray_Alloc(heap, &mut value.alist, 8)? != 0 {
        heap.slice_mut(allocated)?[0] = value;
        MolFmtSgroup_Free(heap, allocated)?;
        return Ok(-1);
    }
    heap.slice_mut(allocated)?[0] = value.clone();
    if IntArray_Alloc(heap, &mut value.blist, 8)? != 0 {
        heap.slice_mut(allocated)?[0] = value;
        MolFmtSgroup_Free(heap, allocated)?;
        return Ok(-1);
    }

    value.id = id;
    value.type_ = type_;
    value.subtype = 0;
    value.conn = 0;
    value.label = 0;
    heap.slice_mut(allocated)?[0] = value;
    Ok(0)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolFmtSgroup_Free(
    heap: &mut SourceHeap,
    sgroup: SourceMutPointer<MOL_FMT_SGROUP>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:658 MolFmtSgroup_Free
    // INCHI✔️❌: complete source frame follows verbatim; checked source-heap object extraction adds overhead.
    /*
void MolFmtSgroup_Free(MOL_FMT_SGROUP *sgroup)
{
    if (sgroup)
    {
        IntArray_Free(&(sgroup->alist));
        IntArray_Free(&(sgroup->blist));
        inchi_free(sgroup);
    }
}
    */
    // END INCHI C FUNCTION: MolFmtSgroup_Free
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolFmtSgroup_Free
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; inchi_free is the mode.h null-checking free macro.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolFmtSgroup_Free

    if sgroup.is_null() {
        return Ok(());
    }
    let mut value = heap
        .slice(sgroup.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    IntArray_Free(heap, Some(&mut value.alist))?;
    IntArray_Free(heap, Some(&mut value.blist))?;
    inchi_free(heap, sgroup)
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
    use crate::source_types::{
        INCHI_IOS_STRING, INCHI_IOS_TYPE_FILE, INCHI_IOS_TYPE_STRING, OAD_Polymer, STR_ERR_LEN,
        SourceFile, inp_ATOM,
    };

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

    fn input_stream(heap: &mut SourceHeap, bytes: &[u8]) -> INCHI_IOSTREAM {
        let data = heap
            .allocate_model_storage(
                bytes
                    .iter()
                    .copied()
                    .chain(std::iter::once(0))
                    .map(|byte| byte as i8)
                    .collect(),
            )
            .unwrap();
        INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            s: INCHI_IOS_STRING {
                pStr: data,
                nAllocatedLength: bytes.len() as i32 + 1,
                nUsedLength: bytes.len() as i32,
                nPtr: 0,
            },
            ..INCHI_IOSTREAM::default()
        }
    }

    fn mutable_c_buffer(heap: &mut SourceHeap, capacity: usize) -> SourceMutPointer<i8> {
        heap.allocate_model_storage(vec![0_i8; capacity]).unwrap()
    }

    fn c_buffer_bytes(heap: &SourceHeap, pointer: SourceConstPointer<i8>) -> Vec<u8> {
        let bytes = heap.slice(pointer).unwrap();
        let length = bytes.iter().position(|byte| *byte == 0).unwrap();
        bytes[..length].iter().map(|byte| *byte as u8).collect()
    }

    fn error_bytes(errors: &[i8]) -> Vec<u8> {
        let length = errors.iter().position(|byte| *byte == 0).unwrap();
        errors[..length].iter().map(|byte| *byte as u8).collect()
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
    fn source_port__mol_fmt4__sdfileidentifylabel__line_346() {
        let mut heap = SourceHeap::default();

        for (line, expected) in [
            (b">  <CAS>\n".as_slice(), SDF_DATA_HEADER_CAS as i32),
            (b">  <cas-number>\n".as_slice(), SDF_DATA_HEADER_CAS as i32),
            (b">  <NAME>\n".as_slice(), SDF_DATA_HEADER as i32),
            (b">  <COMMENT>\n".as_slice(), SDF_DATA_HEADER as i32),
            (b">  <FIELD>\n".as_slice(), SDF_DATA_HEADER as i32),
            (b"no angle brackets".as_slice(), SDF_DATA_HEADER as i32),
            (b">  <>".as_slice(), SDF_DATA_HEADER as i32),
            (b">  <unterminated".as_slice(), SDF_DATA_HEADER as i32),
        ] {
            let input = c_text(&mut heap, line);
            assert_eq!(
                SDFileIdentifyLabel(&mut heap, input, SourceConstPointer::null()),
                Ok(expected),
                "{}",
                String::from_utf8_lossy(line)
            );
        }

        let one_character = c_text(&mut heap, b">  <x>\n");
        let one_character_label = c_text(&mut heap, b"X");
        assert_eq!(
            SDFileIdentifyLabel(&mut heap, one_character, one_character_label),
            Ok(SDF_DATA_HEADER_USER as i32)
        );

        let spaced = c_text(&mut heap, b">  < FIELD >\n");
        let spaced_label = c_text(&mut heap, b" FIELD ");
        assert_eq!(
            SDFileIdentifyLabel(&mut heap, spaced, spaced_label),
            Ok(SDF_DATA_HEADER_USER as i32)
        );

        let ordinary = c_text(&mut heap, b">  <FIELD>\n");
        let ordinary_label = c_text(&mut heap, b"FIELD");
        assert_eq!(
            SDFileIdentifyLabel(&mut heap, ordinary, ordinary_label),
            Ok(SDF_DATA_HEADER as i32)
        );

        let empty_label = c_text(&mut heap, b"");
        assert_eq!(
            SDFileIdentifyLabel(&mut heap, ordinary, empty_label),
            Ok(SDF_DATA_HEADER as i32)
        );

        let long_inner = vec![b'A'; MOL_FMT_MAXLINELEN as usize];
        let mut long_line = b">  <".to_vec();
        long_line.extend_from_slice(&long_inner);
        long_line.push(b'>');
        let long = c_text(&mut heap, &long_line);
        assert_eq!(
            SDFileIdentifyLabel(&mut heap, long, SourceConstPointer::null()),
            Ok(SDF_DATA_HEADER as i32)
        );
        assert_eq!(heap.live_source_allocation_count(), 0);
    }

    #[test]
    fn source_port__mol_fmt4__sdfileextractcasno__line_407() {
        fn check(input: &[u8], expected_value: u64, expected_bytes: &[u8]) {
            let mut heap = SourceHeap::default();
            let storage = heap
                .allocate_model_storage(input.iter().map(|byte| *byte as i8).collect())
                .unwrap();
            assert_eq!(
                SDFileExtractCASNo(&mut heap, storage),
                Ok(expected_value),
                "{}",
                String::from_utf8_lossy(input)
            );
            assert_eq!(
                heap.slice(storage.as_const())
                    .unwrap()
                    .iter()
                    .map(|byte| *byte as u8)
                    .collect::<Vec<_>>(),
                expected_bytes
            );
            heap.free(storage).unwrap();
            assert_eq!(heap.live_source_allocation_count(), 0);
        }

        check(b"\0", 0, b"\0");
        check(b"-\0", 0, b"-\0");
        check(b"123-45\0", 12_345, b"12345\0\0");
        check(b"-123-45\0", 0_u64.wrapping_sub(12_345), b"-12345\0\0");
        check(b"--12\0", 0_u64.wrapping_sub(12), b"-12\0\0");
        check(b"12A34\0\x55\x66", 12, b"12\034\0\x55\x66");
        check(b"A12\0\x55", 0, b"\012\0\x55");
        check(b"+12\0", 0, b"\012\0");
        check(b" 12\0", 0, b"\012\0");
        check(b"1-2--3\0", 123, b"123\0-3\0");
        check(
            b"18446744073709551615\0",
            u64::MAX,
            b"18446744073709551615\0",
        );
        check(
            b"18446744073709551616\0",
            u64::MAX,
            b"18446744073709551616\0",
        );
        check(b"-18446744073709551615\0", 1, b"-18446744073709551615\0");
        check(
            b"-18446744073709551616\0",
            u64::MAX,
            b"-18446744073709551616\0",
        );
        check(b"\xff12\0", 0, b"\012\0");

        let mut heap = SourceHeap::default();
        let unterminated = heap
            .allocate_model_storage(vec![b'1' as i8, b'2' as i8])
            .unwrap();
        assert_eq!(
            SDFileExtractCASNo(&mut heap, unterminated),
            Err(SourceHeapError::MissingNulTerminator)
        );
    }

    #[test]
    fn source_port__mol_fmt4__sdfileskipextradata__line_161() {
        let mut heap = SourceHeap::default();

        let mut direct_end = input_stream(&mut heap, b"$$$$\nNEXT\n");
        let mut cas = 77_u64;
        let name = mutable_c_buffer(&mut heap, 16);
        let comment = mutable_c_buffer(&mut heap, 16);
        let mut errors = [0_i8; STR_ERR_LEN as usize];
        assert_eq!(
            SDFileSkipExtraData(
                &mut heap,
                Some(&mut direct_end),
                Some(&mut cas),
                comment,
                16,
                name,
                16,
                0,
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                Some(&mut errors),
                0,
            ),
            Ok(0)
        );
        assert_eq!(cas, 0);
        assert_eq!(direct_end.s.nPtr, 5);
        assert_eq!(c_buffer_bytes(&heap, name.as_const()), b"");
        assert_eq!(c_buffer_bytes(&heap, comment.as_const()), b"");
        assert_eq!(error_bytes(&errors), b"");

        let mut cas_input = input_stream(
            &mut heap,
            b"M  END\nM  END trailing text\n>  <CAS-number>\n  12-34-5  \n\n$$$$\nNEXT\n",
        );
        cas = 9;
        assert_eq!(
            SDFileSkipExtraData(
                &mut heap,
                Some(&mut cas_input),
                Some(&mut cas),
                SourceMutPointer::null(),
                0,
                SourceMutPointer::null(),
                0,
                0,
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                Some(&mut errors),
                0,
            ),
            Ok(0)
        );
        assert_eq!(cas, 12_345);
        assert_eq!(
            &heap.slice(cas_input.s.pStr.as_const()).unwrap()[cas_input.s.nPtr as usize..][..5],
            &[b'N' as i8, b'E' as i8, b'X' as i8, b'T' as i8, b'\n' as i8]
        );

        let label = c_text(&mut heap, b"x");
        let value = mutable_c_buffer(&mut heap, MAX_SDF_VALUE as usize + 1);
        let mut user_input = input_stream(&mut heap, b">  <X>\n  Alpha\t\tBeta  \n\n$$$$\n");
        assert_eq!(
            SDFileSkipExtraData(
                &mut heap,
                Some(&mut user_input),
                None,
                SourceMutPointer::null(),
                0,
                SourceMutPointer::null(),
                0,
                0,
                label,
                value,
                Some(&mut errors),
                0,
            ),
            Ok(0)
        );
        assert_eq!(c_buffer_bytes(&heap, value.as_const()), b"Alpha Beta");

        let mut eof = input_stream(&mut heap, b">  <FIELD>\ndata\n");
        assert_eq!(
            SDFileSkipExtraData(
                &mut heap,
                Some(&mut eof),
                None,
                SourceMutPointer::null(),
                0,
                SourceMutPointer::null(),
                0,
                0,
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                Some(&mut errors),
                0,
            ),
            Ok(0)
        );
        assert_eq!(eof.s.nPtr, eof.s.nUsedLength);

        errors.fill(0);
        let mut malformed = input_stream(&mut heap, b"bad\x01header\nignored\n$$$$\nNEXT\n");
        assert_eq!(
            SDFileSkipExtraData(
                &mut heap,
                Some(&mut malformed),
                None,
                SourceMutPointer::null(),
                0,
                SourceMutPointer::null(),
                0,
                0,
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                Some(&mut errors),
                0,
            ),
            Ok(9)
        );
        assert_eq!(
            error_bytes(&errors),
            b"Unexpected SData header line: bad.header; Bypassing to next structure"
        );
        assert_eq!(
            &heap.slice(malformed.s.pStr.as_const()).unwrap()[malformed.s.nPtr as usize..][..5],
            &[b'N' as i8, b'E' as i8, b'X' as i8, b'T' as i8, b'\n' as i8]
        );

        errors.fill(0);
        errors[..4].copy_from_slice(&[b'o' as i8, b'l' as i8, b'd' as i8, 0]);
        let mut suppressed = input_stream(&mut heap, b"bad\n$$$$\n");
        assert_eq!(
            SDFileSkipExtraData(
                &mut heap,
                Some(&mut suppressed),
                None,
                SourceMutPointer::null(),
                0,
                SourceMutPointer::null(),
                0,
                1,
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                Some(&mut errors),
                1,
            ),
            Ok(9)
        );
        assert_eq!(error_bytes(&errors), b"old");

        errors.fill(0);
        let mut malformed_eof = input_stream(&mut heap, b"bad\nignored\n");
        assert_eq!(
            SDFileSkipExtraData(
                &mut heap,
                Some(&mut malformed_eof),
                None,
                SourceMutPointer::null(),
                0,
                SourceMutPointer::null(),
                0,
                0,
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                Some(&mut errors),
                0,
            ),
            Ok(3)
        );
        assert_eq!(error_bytes(&errors), b"Unexpected SData header line: bad");

        errors.fill(0);
        let mut long_header_data = vec![b'A'; MOL_FMT_MAXLINELEN as usize + 20];
        long_header_data.extend_from_slice(b"\n$$$$\n");
        let mut long_header = input_stream(&mut heap, &long_header_data);
        assert_eq!(
            SDFileSkipExtraData(
                &mut heap,
                Some(&mut long_header),
                None,
                SourceMutPointer::null(),
                0,
                SourceMutPointer::null(),
                0,
                0,
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                Some(&mut errors),
                0,
            ),
            Ok(9)
        );
        let long_errors = error_bytes(&errors);
        assert!(
            long_errors
                .starts_with(b"Too long SData line truncated; Unexpected SData header line:")
        );
        assert!(long_errors.ends_with(b"Bypassing to next structure"));

        errors.fill(0);
        let mut long_data_input = b">  <FIELD>\n".to_vec();
        long_data_input.extend(std::iter::repeat_n(b'A', MOL_FMT_MAXLINELEN as usize + 20));
        long_data_input.extend_from_slice(b"\n\n$$$$\n");
        let mut long_data = input_stream(&mut heap, &long_data_input);
        assert_eq!(
            SDFileSkipExtraData(
                &mut heap,
                Some(&mut long_data),
                None,
                SourceMutPointer::null(),
                0,
                SourceMutPointer::null(),
                0,
                0,
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                Some(&mut errors),
                0,
            ),
            Ok(0)
        );
        assert_eq!(error_bytes(&errors), b"");

        let file = heap
            .allocate(vec![SourceFile {
                bytes: b"$$$$\nNEXT\n".to_vec(),
                ..SourceFile::default()
            }])
            .unwrap();
        let mut file_stream = INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_FILE as i32,
            f: file,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(
            SDFileSkipExtraData(
                &mut heap,
                Some(&mut file_stream),
                None,
                SourceMutPointer::null(),
                0,
                SourceMutPointer::null(),
                0,
                0,
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                None,
                0,
            ),
            Ok(0)
        );
        assert_eq!(heap.slice(file.as_const()).unwrap()[0].position, 5);
        heap.free(file).unwrap();
    }

    #[test]
    fn source_port__mol_fmt4__numlists_alloc__line_433() {
        let mut heap = SourceHeap::default();
        assert_eq!(NumLists_Alloc(&mut heap, None, 3), Ok(-1));

        let old = heap
            .allocate_model_storage(vec![SourceMutPointer::<i32>::null()])
            .unwrap();
        let mut lists = NUM_LISTS {
            lists: old,
            allocated: 11,
            used: 7,
            increment: 13,
        };
        assert_eq!(NumLists_Alloc(&mut heap, Some(&mut lists), 3), Ok(0));
        assert_ne!(lists.lists, old);
        assert_eq!((lists.allocated, lists.used, lists.increment), (3, 7, 3));
        assert_eq!(
            heap.slice(lists.lists.as_const()).unwrap(),
            &[
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ]
        );
        assert_eq!(
            heap.slice(old.as_const()).unwrap(),
            &[SourceMutPointer::null()]
        );

        let mut zero = NUM_LISTS {
            used: 9,
            ..NUM_LISTS::default()
        };
        assert_eq!(NumLists_Alloc(&mut heap, Some(&mut zero), 0), Ok(0));
        assert!(!zero.lists.is_null());
        assert_eq!((zero.allocated, zero.used, zero.increment), (0, 9, 0));
        assert!(heap.slice(zero.lists.as_const()).unwrap().is_empty());

        let stale = heap
            .allocate_model_storage(vec![SourceMutPointer::<i32>::null()])
            .unwrap();
        let mut negative = NUM_LISTS {
            lists: stale,
            allocated: 4,
            used: 2,
            increment: 6,
        };
        assert_eq!(NumLists_Alloc(&mut heap, Some(&mut negative), -1), Ok(-1));
        assert!(negative.lists.is_null());
        assert_eq!(
            (negative.allocated, negative.used, negative.increment),
            (4, 2, 6)
        );
        assert_eq!(
            heap.slice(stale.as_const()).unwrap(),
            &[SourceMutPointer::null()]
        );

        let leaked = heap
            .allocate_model_storage(vec![SourceMutPointer::<i32>::null()])
            .unwrap();
        let mut failing = NUM_LISTS {
            lists: leaked,
            allocated: 8,
            used: 5,
            increment: 9,
        };
        heap.fail_after_allocations(0);
        assert_eq!(NumLists_Alloc(&mut heap, Some(&mut failing), 2), Ok(-1));
        assert!(failing.lists.is_null());
        assert_eq!(
            (failing.allocated, failing.used, failing.increment),
            (8, 5, 9)
        );
        assert_eq!(
            heap.slice(leaked.as_const()).unwrap(),
            &[SourceMutPointer::null()]
        );
    }

    #[test]
    fn source_port__mol_fmt4__numlists_realloc__line_449() {
        let mut heap = SourceHeap::default();
        assert_eq!(NumLists_ReAlloc(&mut heap, None), Ok(-1));
        let mut invalid = NUM_LISTS::default();
        assert_eq!(NumLists_ReAlloc(&mut heap, Some(&mut invalid)), Ok(-1));

        let first = heap.allocate_model_storage(vec![1_i32]).unwrap();
        let second = heap.allocate_model_storage(vec![2_i32]).unwrap();
        let old = heap.allocate_model_storage(vec![first, second]).unwrap();
        let mut lists = NUM_LISTS {
            lists: old,
            allocated: 2,
            used: 2,
            increment: 2,
        };
        assert_eq!(NumLists_ReAlloc(&mut heap, Some(&mut lists)), Ok(0));
        assert_eq!((lists.allocated, lists.used, lists.increment), (4, 2, 2));
        assert_ne!(lists.lists, old);
        assert_eq!(
            heap.slice(lists.lists.as_const()).unwrap(),
            &[
                first,
                second,
                SourceMutPointer::null(),
                SourceMutPointer::null()
            ]
        );
        assert_eq!(
            heap.slice(old.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let leaked = heap.allocate_model_storage(vec![first]).unwrap();
        let mut failing = NUM_LISTS {
            lists: leaked,
            allocated: 1,
            used: 1,
            increment: 1,
        };
        heap.fail_after_allocations(0);
        assert_eq!(NumLists_ReAlloc(&mut heap, Some(&mut failing)), Ok(-1));
        assert!(failing.lists.is_null());
        assert_eq!(heap.slice(leaked.as_const()).unwrap(), &[first]);
        assert_eq!(
            (failing.allocated, failing.used, failing.increment),
            (1, 1, 1)
        );
    }

    #[test]
    fn source_port__mol_fmt4__numlists_append__line_471() {
        let mut heap = SourceHeap::default();
        let first = heap.allocate_model_storage(vec![11_i32]).unwrap();
        let second = heap.allocate_model_storage(vec![22_i32]).unwrap();
        assert_eq!(NumLists_Append(&mut heap, None, first), Ok(-1));

        let array = heap
            .allocate_model_storage(vec![SourceMutPointer::null(), SourceMutPointer::null()])
            .unwrap();
        let mut lists = NUM_LISTS {
            lists: array,
            allocated: 2,
            used: 0,
            increment: 2,
        };
        assert_eq!(NumLists_Append(&mut heap, Some(&mut lists), first), Ok(0));
        assert_eq!((lists.used, lists.allocated), (1, 2));
        assert_eq!(heap.slice(lists.lists.as_const()).unwrap()[0], first);
        assert_eq!(
            NumLists_Append(&mut heap, Some(&mut lists), SourceMutPointer::null()),
            Ok(0)
        );
        assert_eq!((lists.used, lists.allocated), (2, 2));
        assert!(heap.slice(lists.lists.as_const()).unwrap()[1].is_null());

        assert_eq!(NumLists_Append(&mut heap, Some(&mut lists), second), Ok(0));
        assert_eq!((lists.used, lists.allocated, lists.increment), (3, 4, 2));
        assert_eq!(
            &heap.slice(lists.lists.as_const()).unwrap()[..3],
            &[first, SourceMutPointer::null(), second]
        );
        assert_eq!(
            heap.slice(array.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let old = heap.allocate_model_storage(vec![first]).unwrap();
        let mut failing = NUM_LISTS {
            lists: old,
            allocated: 1,
            used: 1,
            increment: 1,
        };
        heap.fail_after_allocations(0);
        assert_eq!(
            NumLists_Append(&mut heap, Some(&mut failing), second),
            Ok(-1)
        );
        assert!(failing.lists.is_null());
        assert_eq!(
            (failing.used, failing.allocated, failing.increment),
            (1, 1, 1)
        );
        assert_eq!(heap.slice(old.as_const()).unwrap(), &[first]);
    }

    #[test]
    fn source_port__mol_fmt4__numlists_free__line_491() {
        let mut heap = SourceHeap::default();
        assert_eq!(NumLists_Free(&mut heap, None), Ok(()));

        let first = heap.allocate_model_storage(vec![1_i32, 2]).unwrap();
        let second = heap.allocate_model_storage(vec![3_i32]).unwrap();
        let unused = heap.allocate_model_storage(vec![4_i32, 5]).unwrap();
        let list_array = heap
            .allocate_model_storage(vec![first, second, unused])
            .unwrap();
        let mut lists = NUM_LISTS {
            lists: list_array,
            allocated: 3,
            used: 2,
            increment: 7,
        };
        assert_eq!(NumLists_Free(&mut heap, Some(&mut lists)), Ok(()));
        assert_eq!(lists, NUM_LISTS::default());
        assert_eq!(
            heap.slice(first.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(second.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(list_array.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(heap.slice(unused.as_const()).unwrap(), &[4, 5]);

        let preserved = heap.allocate_model_storage(vec![9_i32]).unwrap();
        let empty_array = heap.allocate_model_storage(vec![preserved]).unwrap();
        let mut negative = NUM_LISTS {
            lists: empty_array,
            allocated: 1,
            used: -1,
            increment: i32::MIN,
        };
        assert_eq!(NumLists_Free(&mut heap, Some(&mut negative)), Ok(()));
        assert_eq!(negative, NUM_LISTS::default());
        assert_eq!(heap.slice(preserved.as_const()).unwrap(), &[9]);
        assert_eq!(
            heap.slice(empty_array.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let mut all_zero = NUM_LISTS::default();
        assert_eq!(NumLists_Free(&mut heap, Some(&mut all_zero)), Ok(()));
        assert_eq!(all_zero, NUM_LISTS::default());
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
    fn source_port__mol_fmt4__molfmtsgroups_alloc__line_675() {
        let mut heap = SourceHeap::default();
        assert_eq!(MolFmtSgroups_Alloc(&mut heap, None, 3), Ok(-1));

        let stale = heap
            .allocate_model_storage(vec![SourceMutPointer::<MOL_FMT_SGROUP>::null()])
            .unwrap();
        let mut groups = MOL_FMT_SGROUPS {
            group: stale,
            allocated: 9,
            used: 4,
            increment: 7,
        };
        assert_eq!(MolFmtSgroups_Alloc(&mut heap, Some(&mut groups), 3), Ok(0));
        assert_eq!((groups.allocated, groups.used, groups.increment), (3, 4, 3));
        assert_eq!(
            heap.slice(groups.group.as_const()).unwrap(),
            &[
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null()
            ]
        );
        assert_eq!(heap.slice(stale.as_const()).unwrap().len(), 1);

        let mut zero = MOL_FMT_SGROUPS {
            allocated: -1,
            used: -2,
            increment: -3,
            ..MOL_FMT_SGROUPS::default()
        };
        assert_eq!(MolFmtSgroups_Alloc(&mut heap, Some(&mut zero), 0), Ok(0));
        assert!(!zero.group.is_null());
        assert_eq!((zero.allocated, zero.used, zero.increment), (0, -2, 0));
        assert_eq!(heap.slice(zero.group.as_const()).unwrap(), &[]);

        let mut failure_heap = SourceHeap::default();
        let old = failure_heap
            .allocate_model_storage(vec![SourceMutPointer::<MOL_FMT_SGROUP>::null()])
            .unwrap();
        let mut failure = MOL_FMT_SGROUPS {
            group: old,
            allocated: 8,
            used: 6,
            increment: 4,
        };
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            MolFmtSgroups_Alloc(&mut failure_heap, Some(&mut failure), 2),
            Ok(-1)
        );
        assert!(failure.group.is_null());
        assert_eq!(
            (failure.allocated, failure.used, failure.increment),
            (8, 6, 4)
        );
        assert_eq!(failure_heap.slice(old.as_const()).unwrap().len(), 1);

        let old = heap
            .allocate_model_storage(vec![SourceMutPointer::<MOL_FMT_SGROUP>::null()])
            .unwrap();
        let mut negative = MOL_FMT_SGROUPS {
            group: old,
            allocated: 1,
            used: 1,
            increment: 1,
        };
        assert_eq!(
            MolFmtSgroups_Alloc(&mut heap, Some(&mut negative), -1),
            Ok(-1)
        );
        assert!(negative.group.is_null());
        assert_eq!(
            (negative.allocated, negative.used, negative.increment),
            (1, 1, 1)
        );
        assert_eq!(heap.slice(old.as_const()).unwrap().len(), 1);
    }

    #[test]
    fn source_port__mol_fmt4__molfmtsgroups_realloc__line_693() {
        let mut heap = SourceHeap::default();
        assert_eq!(MolFmtSgroups_ReAlloc(&mut heap, None), Ok(-1));

        for mut invalid in [
            MOL_FMT_SGROUPS::default(),
            MOL_FMT_SGROUPS {
                group: heap
                    .allocate_model_storage(vec![SourceMutPointer::null()])
                    .unwrap(),
                allocated: 0,
                used: 0,
                increment: 1,
            },
            MOL_FMT_SGROUPS {
                group: heap
                    .allocate_model_storage(vec![SourceMutPointer::null()])
                    .unwrap(),
                allocated: 1,
                used: 0,
                increment: 0,
            },
        ] {
            let before = invalid.clone();
            assert_eq!(MolFmtSgroups_ReAlloc(&mut heap, Some(&mut invalid)), Ok(-1));
            assert_eq!(invalid, before);
        }

        let first = heap
            .allocate_model_storage(vec![MOL_FMT_SGROUP::default()])
            .unwrap();
        let second = heap
            .allocate_model_storage(vec![MOL_FMT_SGROUP::default()])
            .unwrap();
        let old = heap
            .allocate(vec![first, second, SourceMutPointer::null()])
            .unwrap();
        let mut groups = MOL_FMT_SGROUPS {
            group: old,
            allocated: 3,
            used: 2,
            increment: 2,
        };
        assert_eq!(MolFmtSgroups_ReAlloc(&mut heap, Some(&mut groups)), Ok(0));
        assert_eq!((groups.allocated, groups.used, groups.increment), (5, 2, 2));
        assert_eq!(
            heap.slice(groups.group.as_const()).unwrap(),
            &[
                first,
                second,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null()
            ]
        );
        assert_eq!(
            heap.slice(old.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let empty_old = heap
            .allocate(Vec::<SourceMutPointer<MOL_FMT_SGROUP>>::new())
            .unwrap();
        let mut empty = MOL_FMT_SGROUPS {
            group: empty_old,
            allocated: 1,
            used: 0,
            increment: 1,
        };
        assert_eq!(MolFmtSgroups_ReAlloc(&mut heap, Some(&mut empty)), Ok(0));
        assert_eq!(empty.allocated, 2);
        assert_eq!(
            heap.slice(empty.group.as_const()).unwrap(),
            &[SourceMutPointer::null(), SourceMutPointer::null()]
        );

        let mut failure_heap = SourceHeap::default();
        let failure_old = failure_heap
            .allocate(vec![SourceMutPointer::<MOL_FMT_SGROUP>::null()])
            .unwrap();
        let mut failure = MOL_FMT_SGROUPS {
            group: failure_old,
            allocated: 1,
            used: 1,
            increment: 2,
        };
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            MolFmtSgroups_ReAlloc(&mut failure_heap, Some(&mut failure)),
            Ok(-1)
        );
        assert!(failure.group.is_null());
        assert_eq!(
            (failure.allocated, failure.used, failure.increment),
            (1, 1, 2)
        );
        assert_eq!(
            failure_heap.slice(failure_old.as_const()).unwrap(),
            &[SourceMutPointer::null()]
        );
        assert_eq!(failure_heap.live_source_allocation_count(), 1);
    }

    #[test]
    fn source_port__mol_fmt4__molfmtsgroups_append__line_715() {
        let mut heap = SourceHeap::default();
        assert_eq!(MolFmtSgroups_Append(&mut heap, None, 1, 2), Ok(-1));
        assert_eq!(heap.live_source_allocation_count(), 0);

        let array = heap
            .allocate(vec![SourceMutPointer::<MOL_FMT_SGROUP>::null(); 2])
            .unwrap();
        let mut groups = MOL_FMT_SGROUPS {
            group: array,
            allocated: 2,
            used: 0,
            increment: 2,
        };
        assert_eq!(
            MolFmtSgroups_Append(&mut heap, Some(&mut groups), 17, 19),
            Ok(0)
        );
        assert_eq!((groups.allocated, groups.used, groups.increment), (2, 1, 2));
        let appended = heap.slice(groups.group.as_const()).unwrap()[0];
        let value = &heap.slice(appended.as_const()).unwrap()[0];
        assert_eq!((value.id, value.type_), (17, 19));
        assert_eq!(heap.live_source_allocation_count(), 4);
        assert_eq!(MolFmtSgroup_Free(&mut heap, appended), Ok(()));
        assert_eq!(inchi_free(&mut heap, groups.group), Ok(()));

        let full_array = heap
            .allocate(vec![SourceMutPointer::<MOL_FMT_SGROUP>::null()])
            .unwrap();
        let mut full = MOL_FMT_SGROUPS {
            group: full_array,
            allocated: 1,
            used: 1,
            increment: 2,
        };
        assert_eq!(
            MolFmtSgroups_Append(&mut heap, Some(&mut full), -7, 31),
            Ok(0)
        );
        assert_eq!((full.allocated, full.used, full.increment), (3, 2, 2));
        let expanded = heap.slice(full.group.as_const()).unwrap();
        assert!(expanded[0].is_null());
        let expanded_group = expanded[1];
        assert!(expanded[2].is_null());
        assert_eq!(
            (
                heap.slice(expanded_group.as_const()).unwrap()[0].id,
                heap.slice(expanded_group.as_const()).unwrap()[0].type_
            ),
            (-7, 31)
        );
        assert_eq!(
            heap.slice(full_array.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(MolFmtSgroup_Free(&mut heap, expanded_group), Ok(()));
        assert_eq!(inchi_free(&mut heap, full.group), Ok(()));
        assert_eq!(heap.live_source_allocation_count(), 0);

        for successful_allocations in 0..=2 {
            let mut failure_heap = SourceHeap::default();
            let old = failure_heap
                .allocate(vec![SourceMutPointer::<MOL_FMT_SGROUP>::null()])
                .unwrap();
            let mut failure = MOL_FMT_SGROUPS {
                group: old,
                allocated: 1,
                used: 0,
                increment: 1,
            };
            failure_heap.fail_after_allocations(successful_allocations);
            assert_eq!(
                MolFmtSgroups_Append(&mut failure_heap, Some(&mut failure), 3, 5),
                Ok(-1)
            );
            assert_eq!(failure.group, old);
            assert_eq!(
                (failure.allocated, failure.used, failure.increment),
                (1, 0, 1)
            );
            assert_eq!(failure_heap.live_source_allocation_count(), 1);
        }

        let mut realloc_failure_heap = SourceHeap::default();
        let realloc_old = realloc_failure_heap
            .allocate(vec![SourceMutPointer::<MOL_FMT_SGROUP>::null()])
            .unwrap();
        let mut realloc_failure = MOL_FMT_SGROUPS {
            group: realloc_old,
            allocated: 1,
            used: 1,
            increment: 1,
        };
        realloc_failure_heap.fail_after_allocations(3);
        assert_eq!(
            MolFmtSgroups_Append(
                &mut realloc_failure_heap,
                Some(&mut realloc_failure),
                23,
                29,
            ),
            Ok(-1)
        );
        assert!(realloc_failure.group.is_null());
        assert_eq!(
            (
                realloc_failure.allocated,
                realloc_failure.used,
                realloc_failure.increment
            ),
            (1, 1, 1)
        );
        assert_eq!(
            realloc_failure_heap.slice(realloc_old.as_const()).unwrap(),
            &[SourceMutPointer::null()]
        );
        assert_eq!(realloc_failure_heap.live_source_allocation_count(), 1);

        let mut invalid_growth_heap = SourceHeap::default();
        let invalid_old = invalid_growth_heap
            .allocate(vec![SourceMutPointer::<MOL_FMT_SGROUP>::null()])
            .unwrap();
        let mut invalid_growth = MOL_FMT_SGROUPS {
            group: invalid_old,
            allocated: 1,
            used: 1,
            increment: 0,
        };
        assert_eq!(
            MolFmtSgroups_Append(&mut invalid_growth_heap, Some(&mut invalid_growth), 2, 3),
            Ok(-1)
        );
        assert_eq!(invalid_growth.group, invalid_old);
        assert_eq!(
            (
                invalid_growth.allocated,
                invalid_growth.used,
                invalid_growth.increment
            ),
            (1, 1, 0)
        );
        assert_eq!(invalid_growth_heap.live_source_allocation_count(), 1);
    }

    #[test]
    fn source_port__mol_fmt4__molfmtsgroups_free__line_751() {
        let mut heap = SourceHeap::default();
        assert_eq!(MolFmtSgroups_Free(&mut heap, None), Ok(()));

        let mut first = SourceMutPointer::null();
        let mut second = SourceMutPointer::null();
        assert_eq!(MolFmtSgroup_Create(&mut heap, &mut first, 1, 2), Ok(0));
        assert_eq!(MolFmtSgroup_Create(&mut heap, &mut second, 3, 4), Ok(0));
        let first_alist = heap.slice(first.as_const()).unwrap()[0].alist.item;
        let first_blist = heap.slice(first.as_const()).unwrap()[0].blist.item;
        let second_alist = heap.slice(second.as_const()).unwrap()[0].alist.item;
        let second_blist = heap.slice(second.as_const()).unwrap()[0].blist.item;
        let array = heap
            .allocate(vec![first, second, SourceMutPointer::null()])
            .unwrap();
        let mut groups = MOL_FMT_SGROUPS {
            group: array,
            allocated: 3,
            used: 2,
            increment: 3,
        };
        assert_eq!(heap.live_source_allocation_count(), 7);
        assert_eq!(MolFmtSgroups_Free(&mut heap, Some(&mut groups)), Ok(()));
        assert_eq!(groups, MOL_FMT_SGROUPS::default());
        assert_eq!(heap.live_source_allocation_count(), 0);
        for pointer_result in [
            heap.slice(first.as_const()).map(|_| ()),
            heap.slice(second.as_const()).map(|_| ()),
            heap.slice(first_alist.as_const()).map(|_| ()),
            heap.slice(first_blist.as_const()).map(|_| ()),
            heap.slice(second_alist.as_const()).map(|_| ()),
            heap.slice(second_blist.as_const()).map(|_| ()),
            heap.slice(array.as_const()).map(|_| ()),
        ] {
            assert_eq!(pointer_result, Err(SourceHeapError::MissingAllocation));
        }

        let empty_array = heap
            .allocate(Vec::<SourceMutPointer<MOL_FMT_SGROUP>>::new())
            .unwrap();
        let mut empty = MOL_FMT_SGROUPS {
            group: empty_array,
            allocated: 0,
            used: 0,
            increment: 0,
        };
        assert_eq!(MolFmtSgroups_Free(&mut heap, Some(&mut empty)), Ok(()));
        assert_eq!(empty, MOL_FMT_SGROUPS::default());
        assert_eq!(heap.live_source_allocation_count(), 0);

        let mut used_child = SourceMutPointer::null();
        let mut unused_child = SourceMutPointer::null();
        assert_eq!(MolFmtSgroup_Create(&mut heap, &mut used_child, 5, 6), Ok(0));
        assert_eq!(
            MolFmtSgroup_Create(&mut heap, &mut unused_child, 7, 8),
            Ok(0)
        );
        let partial_array = heap.allocate(vec![used_child, unused_child]).unwrap();
        let mut partial = MOL_FMT_SGROUPS {
            group: partial_array,
            allocated: 2,
            used: 1,
            increment: 2,
        };
        assert_eq!(MolFmtSgroups_Free(&mut heap, Some(&mut partial)), Ok(()));
        assert_eq!(partial, MOL_FMT_SGROUPS::default());
        assert_eq!(
            heap.slice(used_child.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(heap.slice(unused_child.as_const()).unwrap()[0].id, 7);
        assert_eq!(heap.live_source_allocation_count(), 3);
        assert_eq!(MolFmtSgroup_Free(&mut heap, unused_child), Ok(()));
        assert_eq!(heap.live_source_allocation_count(), 0);

        let mut null_group = MOL_FMT_SGROUPS {
            group: SourceMutPointer::null(),
            allocated: 9,
            used: 0,
            increment: 4,
        };
        assert_eq!(MolFmtSgroups_Free(&mut heap, Some(&mut null_group)), Ok(()));
        assert_eq!(null_group, MOL_FMT_SGROUPS::default());
    }

    #[test]
    fn source_port__mol_fmt4__molfmtsgroups_getindexbysgroupid__line_769() {
        let mut heap = SourceHeap::default();
        let empty = MOL_FMT_SGROUPS::default();
        assert_eq!(MolFmtSgroups_GetIndexBySgroupId(&heap, 1, &empty), Ok(-1));
        let negative_used = MOL_FMT_SGROUPS {
            used: -1,
            ..MOL_FMT_SGROUPS::default()
        };
        assert_eq!(
            MolFmtSgroups_GetIndexBySgroupId(&heap, 1, &negative_used),
            Ok(-1)
        );

        let first = heap
            .allocate_model_storage(vec![MOL_FMT_SGROUP {
                id: 17,
                ..MOL_FMT_SGROUP::default()
            }])
            .unwrap();
        let second = heap
            .allocate_model_storage(vec![MOL_FMT_SGROUP {
                id: -3,
                ..MOL_FMT_SGROUP::default()
            }])
            .unwrap();
        let duplicate = heap
            .allocate_model_storage(vec![MOL_FMT_SGROUP {
                id: 17,
                ..MOL_FMT_SGROUP::default()
            }])
            .unwrap();
        let outside_used = heap
            .allocate_model_storage(vec![MOL_FMT_SGROUP {
                id: 99,
                ..MOL_FMT_SGROUP::default()
            }])
            .unwrap();
        let array = heap
            .allocate_model_storage(vec![first, second, duplicate, outside_used])
            .unwrap();
        let groups = MOL_FMT_SGROUPS {
            group: array,
            allocated: 4,
            used: 3,
            increment: 4,
        };
        assert_eq!(MolFmtSgroups_GetIndexBySgroupId(&heap, 17, &groups), Ok(0));
        assert_eq!(MolFmtSgroups_GetIndexBySgroupId(&heap, -3, &groups), Ok(1));
        assert_eq!(MolFmtSgroups_GetIndexBySgroupId(&heap, 99, &groups), Ok(-1));
        assert_eq!(
            MolFmtSgroups_GetIndexBySgroupId(&heap, 123, &groups),
            Ok(-1)
        );
    }

    #[test]
    fn source_port__mol_fmt4__molfmtsgroup_free__line_658() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            MolFmtSgroup_Free(&mut heap, SourceMutPointer::null()),
            Ok(())
        );

        let alist = heap.allocate(vec![1_i32, 2, 3]).unwrap();
        let blist = heap.allocate(vec![4_i32, 5]).unwrap();
        let group = heap
            .allocate(vec![MOL_FMT_SGROUP {
                id: 17,
                alist: INT_ARRAY {
                    item: alist,
                    allocated: 3,
                    used: 3,
                    increment: 3,
                },
                blist: INT_ARRAY {
                    item: blist,
                    allocated: 2,
                    used: 2,
                    increment: 2,
                },
                ..MOL_FMT_SGROUP::default()
            }])
            .unwrap();
        assert_eq!(heap.live_source_allocation_count(), 3);
        assert_eq!(MolFmtSgroup_Free(&mut heap, group), Ok(()));
        assert_eq!(heap.live_source_allocation_count(), 0);
        assert_eq!(
            heap.slice(alist.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(blist.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(group.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let group_with_null_lists = heap
            .allocate(vec![MOL_FMT_SGROUP {
                id: 18,
                ..MOL_FMT_SGROUP::default()
            }])
            .unwrap();
        assert_eq!(MolFmtSgroup_Free(&mut heap, group_with_null_lists), Ok(()));
        assert_eq!(heap.live_source_allocation_count(), 0);
    }

    #[test]
    fn source_port__mol_fmt4__molfmtsgroup_create__line_634() {
        let mut heap = SourceHeap::default();
        let stale = heap
            .allocate_model_storage(vec![MOL_FMT_SGROUP::default()])
            .unwrap();
        let mut group = stale;
        assert_eq!(MolFmtSgroup_Create(&mut heap, &mut group, -17, 23), Ok(0));
        assert_ne!(group, stale);
        let created = &heap.slice(group.as_const()).unwrap()[0];
        assert_eq!((created.id, created.type_), (-17, 23));
        assert_eq!((created.subtype, created.conn, created.label), (0, 0, 0));
        assert_eq!(created.xbr1, [0.0; 4]);
        assert_eq!(created.xbr2, [0.0; 4]);
        assert_eq!(created.smt, [0; 80]);
        assert_eq!(
            (
                created.alist.allocated,
                created.alist.used,
                created.alist.increment
            ),
            (8, 0, 8)
        );
        assert_eq!(
            (
                created.blist.allocated,
                created.blist.used,
                created.blist.increment
            ),
            (8, 0, 8)
        );
        assert_eq!(heap.slice(created.alist.item.as_const()).unwrap(), &[0; 8]);
        assert_eq!(heap.slice(created.blist.item.as_const()).unwrap(), &[0; 8]);
        assert_eq!(heap.live_source_allocation_count(), 3);
        assert_eq!(MolFmtSgroup_Free(&mut heap, group), Ok(()));
        assert_eq!(heap.live_source_allocation_count(), 0);

        for successful_allocations in 0..=2 {
            let mut failure_heap = SourceHeap::default();
            let old = failure_heap
                .allocate_model_storage(vec![MOL_FMT_SGROUP::default()])
                .unwrap();
            let mut output = old;
            failure_heap.fail_after_allocations(successful_allocations);
            assert_eq!(
                MolFmtSgroup_Create(&mut failure_heap, &mut output, 7, 11),
                Ok(-1),
                "failure after {successful_allocations} successful allocations"
            );
            if successful_allocations == 0 {
                assert!(output.is_null());
            } else {
                assert!(!output.is_null());
                assert_eq!(
                    failure_heap.slice(output.as_const()),
                    Err(SourceHeapError::MissingAllocation)
                );
            }
            assert_eq!(failure_heap.live_source_allocation_count(), 0);
            assert_eq!(failure_heap.slice(old.as_const()).unwrap().len(), 1);
        }
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
