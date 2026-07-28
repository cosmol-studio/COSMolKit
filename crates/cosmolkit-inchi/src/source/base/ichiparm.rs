use crate::source::base::ichi_io::{inchi_ios_eprint, inchi_ios_print_nodisplay};
use crate::source::base::util::{
    inchi__strdup, inchi_free, inchi_malloc, inchi_memicmp, inchi_stricmp, lrtrim, mystrncpy,
};
use crate::source_types::*;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct CommonOptionsByParg {
    pub(crate) ver1_default_mode: INCHI_MODE,
    pub(crate) mode: i32,
    pub(crate) inchi_output_options: i32,
    pub(crate) inchi_output_options2: i32,
    pub(crate) std_format: i32,
    pub(crate) hash_key: i32,
    pub(crate) hash_xtra1: i32,
    pub(crate) hash_xtra2: i32,
    pub(crate) fix_sp3_bug: i32,
    pub(crate) fix_fb2: i32,
    pub(crate) add_phosphine_stereo: i32,
    pub(crate) add_arsine_stereo: i32,
    pub(crate) no_struct_labels: i32,
    pub(crate) pointed_edge_stereo: i32,
    pub(crate) do_not_add_h: i32,
    pub(crate) forced_chiral_flag: i32,
    pub(crate) reconnect_coord: i32,
    pub(crate) keto_enol_taut: i32,
    pub(crate) taut_15_non_ring: i32,
    pub(crate) pt_06_00_taut: i32,
    pub(crate) pt_13_00_taut: i32,
    pub(crate) pt_16_00_taut: i32,
    pub(crate) pt_18_00_taut: i32,
    pub(crate) pt_22_00_taut: i32,
    pub(crate) pt_39_00_taut: i32,
    pub(crate) loose_tsa_check: i32,
    pub(crate) large_molecules: i32,
    pub(crate) polymers: i32,
    pub(crate) fold_polymer_sru: i32,
    pub(crate) frame_shift_scheme: i32,
    pub(crate) stereo_at_zz: i32,
    pub(crate) np_zz: i32,
    pub(crate) no_warnings: i32,
    pub(crate) merge_hash: i32,
    pub(crate) hide_inchi: i32,
}

fn option_literal(
    heap: &mut SourceHeap,
    argument: SourceConstPointer<i8>,
    literal: &str,
) -> Result<bool, SourceHeapError> {
    let literal = heap.allocate_model_storage(
        literal
            .bytes()
            .chain(std::iter::once(0))
            .map(|byte| byte as i8)
            .collect(),
    )?;
    let result = inchi_stricmp(heap, argument, literal.as_const());
    heap.free(literal)?;
    result.map(|ordering| ordering == 0)
}

fn option_prefix(
    heap: &mut SourceHeap,
    argument: SourceConstPointer<i8>,
    literal: &str,
) -> Result<bool, SourceHeapError> {
    let literal_pointer = heap.allocate_model_storage(
        literal
            .bytes()
            .chain(std::iter::once(0))
            .map(|byte| byte as i8)
            .collect(),
    )?;
    let available = heap.slice(argument)?.len();
    let padded_argument = if available < literal.len() {
        let mut bytes = heap.slice(argument)?.to_vec();
        bytes.resize(literal.len(), 0);
        Some(heap.allocate_model_storage(bytes)?)
    } else {
        None
    };
    let compared_argument = padded_argument
        .map(SourceMutPointer::as_const)
        .unwrap_or(argument);
    let result = inchi_memicmp(
        heap,
        compared_argument,
        literal_pointer.as_const(),
        literal.len() as u64,
    );
    heap.free(literal_pointer)?;
    if let Some(padded_argument) = padded_argument {
        heap.free(padded_argument)?;
    }
    result.map(|ordering| ordering == 0)
}

fn source_strtol_long(
    heap: &mut SourceHeap,
    pointer: SourceConstPointer<i8>,
) -> Result<i64, SourceHeapError> {
    source_strtol_long_with_end(heap, pointer).map(|(value, _)| value)
}

fn source_strtol_long_with_end(
    heap: &mut SourceHeap,
    pointer: SourceConstPointer<i8>,
) -> Result<(i64, usize), SourceHeapError> {
    let bytes = heap.slice(pointer)?.to_vec();
    let nul = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let mut index = 0_usize;
    while index < nul
        && matches!(
            bytes[index] as u8,
            b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r'
        )
    {
        index += 1;
    }
    let negative = if index < nul && matches!(bytes[index] as u8, b'+' | b'-') {
        let negative = bytes[index] as u8 == b'-';
        index += 1;
        negative
    } else {
        false
    };
    let digit_start = index;
    let mut value = 0_i64;
    let mut overflowed = false;
    while index < nul && (bytes[index] as u8).is_ascii_digit() {
        let digit = i64::from((bytes[index] as u8) - b'0');
        let next = if negative {
            value
                .checked_mul(10)
                .and_then(|value| value.checked_sub(digit))
        } else {
            value
                .checked_mul(10)
                .and_then(|value| value.checked_add(digit))
        };
        value = next.unwrap_or_else(|| {
            overflowed = true;
            if negative { i64::MIN } else { i64::MAX }
        });
        index += 1;
    }
    if overflowed {
        heap.set_source_errno(34);
    }
    if index == digit_start {
        Ok((0, 0))
    } else {
        Ok((value, index))
    }
}

pub(crate) fn source_strtod_with_end(
    heap: &mut SourceHeap,
    pointer: SourceConstPointer<i8>,
) -> Result<(f64, usize), SourceHeapError> {
    let bytes = heap.slice(pointer)?.to_vec();
    let nul = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let mut index = 0_usize;
    while index < nul
        && matches!(
            bytes[index] as u8,
            b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r'
        )
    {
        index += 1;
    }
    let token_start = index;
    if index < nul && matches!(bytes[index] as u8, b'+' | b'-') {
        index += 1;
    }
    let number_start = index;
    let remaining = &bytes[number_start..nul];
    let lower = |byte: i8| (byte as u8).to_ascii_lowercase();
    if remaining.len() >= 3
        && lower(remaining[0]) == b'i'
        && lower(remaining[1]) == b'n'
        && lower(remaining[2]) == b'f'
    {
        index += 3;
        let tail = &bytes[index..nul];
        if tail.len() >= 5
            && tail[..5]
                .iter()
                .map(|byte| lower(*byte))
                .eq(b"inity".iter().copied())
        {
            index += 5;
        }
        let negative = bytes[token_start] as u8 == b'-';
        return Ok((
            if negative {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            },
            index,
        ));
    }
    if remaining.len() >= 3
        && lower(remaining[0]) == b'n'
        && lower(remaining[1]) == b'a'
        && lower(remaining[2]) == b'n'
    {
        index += 3;
        let negative = bytes[token_start] as u8 == b'-';
        let mut payload = 0_u64;
        if bytes.get(index).map(|byte| *byte as u8) == Some(b'(')
            && let Some(close) = bytes[index + 1..nul]
                .iter()
                .position(|byte| *byte as u8 == b')')
        {
            let payload_bytes = &bytes[index + 1..index + 1 + close];
            if payload_bytes.iter().all(|byte| {
                let byte = *byte as u8;
                byte.is_ascii_alphanumeric() || byte == b'_'
            }) {
                index += close + 2;

                let (base, digit_start) = if payload_bytes.len() > 2
                    && payload_bytes[0] as u8 == b'0'
                    && matches!(payload_bytes[1] as u8, b'x' | b'X')
                {
                    (16_u64, 2_usize)
                } else if payload_bytes.first().map(|byte| *byte as u8) == Some(b'0') {
                    (8_u64, 0_usize)
                } else {
                    (10_u64, 0_usize)
                };
                let mut parsed = 0_u64;
                let mut saw_digit = false;
                let mut valid_number = digit_start < payload_bytes.len();
                for &byte in &payload_bytes[digit_start..] {
                    let digit = match byte as u8 {
                        b'0'..=b'9' => u64::from(byte as u8 - b'0'),
                        b'a'..=b'f' if base == 16 => u64::from(byte as u8 - b'a') + 10,
                        b'A'..=b'F' if base == 16 => u64::from(byte as u8 - b'A') + 10,
                        _ => {
                            valid_number = false;
                            break;
                        }
                    };
                    if digit >= base {
                        valid_number = false;
                        break;
                    }
                    saw_digit = true;
                    parsed = parsed
                        .checked_mul(base)
                        .and_then(|value| value.checked_add(digit))
                        .unwrap_or(u64::MAX);
                }
                if valid_number && saw_digit {
                    payload = parsed;
                }
            }
        }
        let bits =
            u64::from(negative) << 63 | 0x7ff8_0000_0000_0000 | (payload & 0x000f_ffff_ffff_ffff);
        return Ok((f64::from_bits(bits), index));
    }

    if bytes
        .get(index)
        .map(|byte| (*byte as u8).to_ascii_lowercase())
        == Some(b'0')
        && bytes
            .get(index + 1)
            .map(|byte| (*byte as u8).to_ascii_lowercase())
            == Some(b'x')
    {
        index += 2;
        let mut digits = 0_usize;
        let mut has_nonzero_digit = false;
        while index < nul && (bytes[index] as u8).is_ascii_hexdigit() {
            has_nonzero_digit |= bytes[index] as u8 != b'0';
            digits += 1;
            index += 1;
        }
        if index < nul && bytes[index] as u8 == b'.' {
            index += 1;
            while index < nul && (bytes[index] as u8).is_ascii_hexdigit() {
                has_nonzero_digit |= bytes[index] as u8 != b'0';
                digits += 1;
                index += 1;
            }
        }
        if digits == 0 {
            return Ok((0.0, 0));
        }
        let mut has_complete_exponent = false;
        if index < nul && matches!(bytes[index] as u8, b'p' | b'P') {
            let exponent = index;
            index += 1;
            if index < nul && matches!(bytes[index] as u8, b'+' | b'-') {
                index += 1;
            }
            let exponent_digits = index;
            while index < nul && (bytes[index] as u8).is_ascii_digit() {
                index += 1;
            }
            if index == exponent_digits {
                index = exponent;
            } else {
                has_complete_exponent = true;
            }
        }
        let mut token: Vec<u8> = bytes[token_start..index]
            .iter()
            .map(|byte| *byte as u8)
            .collect();
        if !has_complete_exponent {
            token.extend_from_slice(b"p0");
        }
        let token =
            std::str::from_utf8(&token).map_err(|_| SourceHeapError::InvalidSourceTextEncoding)?;
        let value = hexf_parse::parse_hexf64(token, false)
            .map_err(|_| SourceHeapError::UnsupportedSourceBehavior)?;
        if value.is_infinite()
            || (has_nonzero_digit && (value == 0.0 || value.abs() < f64::MIN_POSITIVE))
        {
            heap.set_source_errno(34);
        }
        return Ok((value, index));
    }

    let mut has_nonzero_digit = false;
    let mut digits = 0_usize;
    while index < nul && (bytes[index] as u8).is_ascii_digit() {
        has_nonzero_digit |= bytes[index] as u8 != b'0';
        digits += 1;
        index += 1;
    }
    if index < nul && bytes[index] as u8 == b'.' {
        index += 1;
        while index < nul && (bytes[index] as u8).is_ascii_digit() {
            has_nonzero_digit |= bytes[index] as u8 != b'0';
            digits += 1;
            index += 1;
        }
    }
    if digits == 0 {
        return Ok((0.0, 0));
    }
    if index < nul && matches!(bytes[index] as u8, b'e' | b'E') {
        let exponent = index;
        index += 1;
        if index < nul && matches!(bytes[index] as u8, b'+' | b'-') {
            index += 1;
        }
        let exponent_digits = index;
        while index < nul && (bytes[index] as u8).is_ascii_digit() {
            index += 1;
        }
        if index == exponent_digits {
            index = exponent;
        }
    }
    let token: Vec<u8> = bytes[token_start..index]
        .iter()
        .map(|byte| *byte as u8)
        .collect();
    let token =
        std::str::from_utf8(&token).map_err(|_| SourceHeapError::InvalidSourceTextEncoding)?;
    let value = token
        .parse::<f64>()
        .map_err(|_| SourceHeapError::UnsupportedSourceBehavior)?;
    if value.is_infinite()
        || (has_nonzero_digit && (value == 0.0 || value.abs() < f64::MIN_POSITIVE))
    {
        heap.set_source_errno(34);
    }
    Ok((value, index))
}

#[allow(non_snake_case)]
pub(crate) fn set_common_options_by_parg(
    heap: &mut SourceHeap,
    pArg: SourceConstPointer<i8>,
    developer_options: i32,
    ip: &mut INPUT_PARMS,
    options: &mut CommonOptionsByParg,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c:131 set_common_options_by_parg
    // INCHI✔❌: int set_common_options_by_parg(const char* pArg,
    // INCHI✔❌:     int  developer_options,
    // INCHI✔❌:     INPUT_PARMS* ip,
    // INCHI✔❌:     INCHI_MODE* pbVer1DefaultMode,
    // INCHI✔❌:     int* pnMode,
    // INCHI✔❌:     int* pbINChIOutputOptions,
    // INCHI✔❌:     int* pbINChIOutputOptions2,
    // INCHI✔❌:     int* pbStdFormat,
    // INCHI✔❌:     int* pbHashKey,
    // INCHI✔❌:     int* pbHashXtra1,
    // INCHI✔❌:     int* pbHashXtra2,
    // INCHI✔❌:     int* pbFixSp3bug,
    // INCHI✔❌:     int* pbFixFB2,
    // INCHI✔❌:     int* pbAddPhosphineStereo,
    // INCHI✔❌:     int* pbAddArsineStereo,
    // INCHI✔❌:     int* pbNoStructLabels,
    // INCHI✔❌:     int* pbPointedEdgeStereo,
    // INCHI✔❌:     int* pbDoNotAddH,
    // INCHI✔❌:     int* pbForcedChiralFlag,
    // INCHI✔❌:     int* pbReconnectCoord,
    // INCHI✔❌:     int* pbKetoEnolTaut,
    // INCHI✔❌:     int* pb15TautNonRing,
    // INCHI✔❌:     int* pbPT_06_00_Taut,
    // INCHI✔❌:     int* pbPT_13_00_Taut,
    // INCHI✔❌:     int* pbPT_16_00_Taut,
    // INCHI✔❌:     int* pbPT_18_00_Taut,
    // INCHI✔❌:     int* pbPT_22_00_Taut,
    // INCHI✔❌:     int* pbPT_39_00_Taut,
    // INCHI✔❌:     int* pbLooseTSACheck,
    // INCHI✔❌:     int* pbLargeMolecules,
    // INCHI✔❌:     int* pbPolymers,
    // INCHI✔❌:     int* pbFoldPolymerSRU,
    // INCHI✔❌:     int* pbFrameShiftScheme,
    // INCHI✔❌:     int* pbStereoAtZz,
    // INCHI✔❌:     int* pbNPZz,
    // INCHI✔❌:     int* pbNoWarnings,
    // INCHI✔❌:     int* pbMergeHash,
    // INCHI✔❌:     int* pbHideInChI
    // INCHI✔❌: )
    // INCHI✔❌: {
    // INCHI✔❌:     int got = 0;
    // INCHI✔❌:
    // INCHI✔❌:     /* Input */
    // INCHI✔❌:     if (!inchi_stricmp(pArg, "INPAUX"))
    // INCHI✔❌:     {
    // INCHI✔❌:         if (INPUT_NONE == ip->nInputType)
    // INCHI✔❌:         {
    // INCHI✔❌:             ip->nInputType = INPUT_INCHI_PLAIN;
    // INCHI✔❌:         }
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_memicmp(pArg, "START:", 6))
    // INCHI✔❌:     {
    // INCHI✔❌:         ip->first_struct_number = strtol(pArg + 6, NULL, 10);
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_memicmp(pArg, "END:", 4))
    // INCHI✔❌:     {
    // INCHI✔❌:         ip->last_struct_number = strtol(pArg + 4, NULL, 10);
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_memicmp(pArg, "RECORD:", 7))
    // INCHI✔❌:     {
    // INCHI✔❌:         long num = strtol(pArg + 7, NULL, 10);
    // INCHI✔❌:         /* djb-rwth: removing redundant code */
    // INCHI✔❌:         ip->first_struct_number = num;
    // INCHI✔❌:         ip->last_struct_number = num;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "NOLABELS"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbNoStructLabels = 1;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "SAVEOPT"))
    // INCHI✔❌:     {
    // INCHI✔❌:         (*pbINChIOutputOptions) |= INCHI_OUT_SAVEOPT;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     /* Generation */
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "AUXNONE"))
    // INCHI✔❌:     {
    // INCHI✔❌:         /* no aux. info */
    // INCHI✔❌:         (*pbINChIOutputOptions) |= INCHI_OUT_NO_AUX_INFO;  /* no aux info */
    // INCHI✔❌:         (*pbINChIOutputOptions) &= ~INCHI_OUT_SHORT_AUX_INFO;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "MISMATCHISERROR"))
    // INCHI✔❌:     {
    // INCHI✔❌:         /* Consider InChI conversion "problem/mismatch" as error */
    // INCHI✔❌:         (*pbINChIOutputOptions2) |= INCHI_OUT_MISMATCH_AS_ERROR;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "OUTERRINCHI"))
    // INCHI✔❌:     {
    // INCHI✔❌:         /* Signify InChI error generation on InChI strings output, not only report to log file */
    // INCHI✔❌:         (*pbINChIOutputOptions2) |= INCHI_OUT_INCHI_GEN_ERROR;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     /* InChIKey/InChI hash */
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "Key"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbHashKey = 1;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "XHash1"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbHashXtra1 = 1;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "XHash2"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbHashXtra2 = 1;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     /* All modes (std and non-std InChI) structure perception options */
    // INCHI✔❌:     /* These options DO NOT TURN OFF Std flag                         */
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "SNON"))
    // INCHI✔❌:     {
    // INCHI✔❌:         (*pbVer1DefaultMode) &= ~REQ_MODE_STEREO; /* no stereo */
    // INCHI✔❌:         (*pnMode) &= ~(REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO);
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "NEWPSOFF"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbPointedEdgeStereo = 0;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "DONOTADDH"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbDoNotAddH = 1;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "LooseTSACheck"))
    // INCHI✔❌:     {
    // INCHI✔❌:         (*pbLooseTSACheck) = 1;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌: #ifndef USE_STDINCHI_API
    // INCHI✔❌:     /* These options DO TURN OFF Std flag   */
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "SREL"))
    // INCHI✔❌:     {
    // INCHI✔❌:         if ((*pnMode) & REQ_MODE_RACEMIC_STEREO)
    // INCHI✔❌:         {
    // INCHI✔❌:             (*pnMode) ^= REQ_MODE_RACEMIC_STEREO;
    // INCHI✔❌:         }
    // INCHI✔❌:         if ((*pnMode) & REQ_MODE_CHIR_FLG_STEREO)
    // INCHI✔❌:         {
    // INCHI✔❌:             (*pnMode) ^= REQ_MODE_CHIR_FLG_STEREO;
    // INCHI✔❌:         }
    // INCHI✔❌:         (*pnMode) |= REQ_MODE_RELATIVE_STEREO;
    // INCHI✔❌:         (*pnMode) |= REQ_MODE_STEREO;
    // INCHI✔❌:         *pbStdFormat = 0;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "SRAC"))
    // INCHI✔❌:     {
    // INCHI✔❌:         if ((*pnMode) & REQ_MODE_RELATIVE_STEREO)
    // INCHI✔❌:         {
    // INCHI✔❌:             (*pnMode) ^= REQ_MODE_RELATIVE_STEREO;
    // INCHI✔❌:         }
    // INCHI✔❌:         if ((*pnMode) & REQ_MODE_CHIR_FLG_STEREO)
    // INCHI✔❌:         {
    // INCHI✔❌:             (*pnMode) ^= REQ_MODE_CHIR_FLG_STEREO;
    // INCHI✔❌:         }
    // INCHI✔❌:         (*pnMode) |= REQ_MODE_RACEMIC_STEREO;
    // INCHI✔❌:         (*pnMode) |= REQ_MODE_STEREO;
    // INCHI✔❌:         *pbStdFormat = 0;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "SUCF"))
    // INCHI✔❌:     {
    // INCHI✔❌:         if ((*pnMode) & REQ_MODE_RELATIVE_STEREO)
    // INCHI✔❌:         {
    // INCHI✔❌:             (*pnMode) ^= REQ_MODE_RELATIVE_STEREO;
    // INCHI✔❌:         }
    // INCHI✔❌:         if ((*pnMode) & REQ_MODE_RACEMIC_STEREO)
    // INCHI✔❌:         {
    // INCHI✔❌:             (*pnMode) ^= REQ_MODE_RACEMIC_STEREO;
    // INCHI✔❌:         }
    // INCHI✔❌:         (*pnMode) |= REQ_MODE_CHIR_FLG_STEREO; /* stereo defined by the Chiral flag */
    // INCHI✔❌:         (*pnMode) |= REQ_MODE_STEREO;
    // INCHI✔❌:         *pbStdFormat = 0;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "ChiralFlagON"))
    // INCHI✔❌:     {
    // INCHI✔❌:         /* used only with /SUCF */
    // INCHI✔❌:         /* NB: do not toggle off bStdFormat! (if necessary SUCF will do) */
    // INCHI✔❌:         (*pbForcedChiralFlag) &= ~FLAG_SET_INP_AT_NONCHIRAL;
    // INCHI✔❌:         (*pbForcedChiralFlag) |= FLAG_SET_INP_AT_CHIRAL;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "ChiralFlagOFF"))
    // INCHI✔❌:     {
    // INCHI✔❌:         /* used only with /SUCF */
    // INCHI✔❌:         /* NB: do not toggle off bStdFormat! (if necessary SUCF will do) */
    // INCHI✔❌:         (*pbForcedChiralFlag) &= ~FLAG_SET_INP_AT_CHIRAL;
    // INCHI✔❌:         (*pbForcedChiralFlag) |= FLAG_SET_INP_AT_NONCHIRAL;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /*--- Non-std InChI creation options ---*/
    // INCHI✔❌:     /* These options DO TURN OFF Std flag       */
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "SUU"))
    // INCHI✔❌:     {
    // INCHI✔❌:         /* include omitted undef/unknown stereo */
    // INCHI✔❌:         (*pbVer1DefaultMode) &= ~(REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU);
    // INCHI✔❌:         *pbStdFormat = 0;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "SLUUD"))
    // INCHI✔❌:     {
    // INCHI✔❌:         /* make labels for unknown and undefined stereo different */
    // INCHI✔❌:         (*pbVer1DefaultMode) |= REQ_MODE_DIFF_UU_STEREO;
    // INCHI✔❌:         *pbStdFormat = 0;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     /* FixedH */
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "FIXEDH"))
    // INCHI✔❌:     {
    // INCHI✔❌:         (*pbVer1DefaultMode) |= REQ_MODE_BASIC;  /* tautomeric */
    // INCHI✔❌:         *pbStdFormat = 0;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     /* RecMet */
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "RECMET"))
    // INCHI✔❌:     {
    // INCHI✔❌:         /* reconnect metals */
    // INCHI✔❌:         *pbReconnectCoord = 1;
    // INCHI✔❌:         *pbStdFormat = 0;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌: #if ( KETO_ENOL_TAUT == 1 )
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "KET"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbKetoEnolTaut = 1;
    // INCHI✔❌:         *pbStdFormat = 0;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌: #endif
    // INCHI✔❌: #if ( TAUT_15_NON_RING == 1 )
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "15T"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pb15TautNonRing = 1;
    // INCHI✔❌:         *pbStdFormat = 0;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌: #if ( TAUT_PT_22_00  == 1 )
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "PT_22_00"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbPT_22_00_Taut = 1;
    // INCHI✔❌:         *pbStdFormat = 0;
    // INCHI✔❌:     }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌: #if ( TAUT_PT_16_00  == 1 )
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "PT_16_00"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbPT_16_00_Taut = 1;
    // INCHI✔❌:         *pbStdFormat = 0;
    // INCHI✔❌:     }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌: #if ( TAUT_PT_06_00  == 1 )
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "PT_06_00"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbPT_06_00_Taut = 1;
    // INCHI✔❌:         *pbStdFormat = 0;
    // INCHI✔❌:     }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌: #if ( TAUT_PT_39_00  == 1 )
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "PT_39_00"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbPT_39_00_Taut = 1;
    // INCHI✔❌:         *pbStdFormat = 0;
    // INCHI✔❌:     }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌: #if ( TAUT_PT_13_00  == 1 )
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "PT_13_00"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbPT_13_00_Taut = 1;
    // INCHI✔❌:         *pbStdFormat = 0;
    // INCHI✔❌:     }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌: #if ( TAUT_PT_18_00  == 1 )
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "PT_18_00"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbPT_18_00_Taut = 1;
    // INCHI✔❌:         *pbStdFormat = 0;
    // INCHI✔❌:     }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "LargeMolecules"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbLargeMolecules = 1;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "Polymers"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbPolymers = POLYMERS_MODERN;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "Polymers105"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbPolymers = POLYMERS_LEGACY;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "NPZz"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbNPZz = 1;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "NoWarnings"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbNoWarnings = 1;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "MergeHash"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbMergeHash = 1;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "NoInChI") || !inchi_stricmp(pArg, "HideInChI"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbHideInChI = 1;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "FoldCRU")) /* v. 1.06 */
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbFoldPolymerSRU = 1;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "FoldSRU")) /* v. 1.06 */
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbFoldPolymerSRU = 1;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_memicmp(pArg, "FrameShift:", 11))
    // INCHI✔❌:     {
    // INCHI✔❌:         int k;
    // INCHI✔❌:         char wrd[256];
    // INCHI✔❌:         k = 0;
    // INCHI✔❌:         mystrncpy(wrd, pArg + 11, 256);
    // INCHI✔❌:         lrtrim(wrd, &k);
    // INCHI✔❌:         if (k)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (!inchi_stricmp(wrd, "None"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 *pbFrameShiftScheme = FSS_NONE;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(wrd, "Cyclize"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 *pbFrameShiftScheme = FSS_STARS_CYCLED;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(wrd, "MoveStars"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 *pbFrameShiftScheme = FSS_STARS_OPENED;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(wrd, "MoveBrackets"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 *pbFrameShiftScheme = FSS_STARS_ENDS_OPENED;
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:         {
    // INCHI✔❌:             *pbFrameShiftScheme = FSS_STARS_CYCLED;
    // INCHI✔❌:         }
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "NoFrameShift"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbFrameShiftScheme = FSS_NONE;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "NoEdits"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbFoldPolymerSRU = 0;
    // INCHI✔❌:         *pbFrameShiftScheme = FSS_NONE;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!inchi_stricmp(pArg, "SATZZ"))
    // INCHI✔❌:     {
    // INCHI✔❌:         *pbStereoAtZz = 1;
    // INCHI✔❌:         got = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌: #endif /* ifndef USE_STDINCHI_API */
    // INCHI✔❌:
    // INCHI✔❌:     if (!got && developer_options)
    // INCHI✔❌:     {
    // INCHI✔❌:
    // INCHI✔❌:         if (!inchi_stricmp(pArg, "PGO"))
    // INCHI✔❌:         {
    // INCHI✔❌:             /* PGO : extract all good MOLfiles into the problem file */
    // INCHI✔❌:             ip->bSaveAllGoodStructsAsProblem = 1;
    // INCHI✔❌:             got = 1;
    // INCHI✔❌:         }
    // INCHI✔❌: #if ( ALLOW_SUBSTRUCTURE_FILTERING== 1 )
    // INCHI✔❌:         else if (!inchi_stricmp(pArg, "FilterSS"))
    // INCHI✔❌:         {
    // INCHI✔❌:             ip->bFilterSS = 1;
    // INCHI✔❌:             got = 1;
    // INCHI✔❌:         }
    // INCHI✔❌:         else if (!inchi_stricmp(pArg, "InvFilterSS"))
    // INCHI✔❌:         {
    // INCHI✔❌:             ip->bFilterSS = -1;
    // INCHI✔❌:             got = 1;
    // INCHI✔❌:         }
    // INCHI✔❌: #endif
    // INCHI✔❌:         /* Options below DO TURN OFF Std flag   */
    // INCHI✔❌:         if (!inchi_stricmp(pArg, "FNUDOFF"))
    // INCHI✔❌:         {
    // INCHI✔❌:             ip->bFixNonUniformDraw = 0;
    // INCHI✔❌:             *pbStdFormat = 0;
    // INCHI✔❌:             got = 1;
    // INCHI✔❌:         }
    // INCHI✔❌:         else if (!inchi_stricmp(pArg, "FixSp3bugOFF"))
    // INCHI✔❌:         {
    // INCHI✔❌:             *pbFixSp3bug = 0;
    // INCHI✔❌:             *pbStdFormat = 0;
    // INCHI✔❌:             got = 1;
    // INCHI✔❌:         }
    // INCHI✔❌:         else if (!inchi_stricmp(pArg, "FBOFF"))
    // INCHI✔❌:         {
    // INCHI✔❌:             *pbFixSp3bug = 0;
    // INCHI✔❌:             *pbStdFormat = 0;
    // INCHI✔❌:             got = 1;
    // INCHI✔❌:         }
    // INCHI✔❌:         else if (!inchi_stricmp(pArg, "FB2OFF"))
    // INCHI✔❌:         {
    // INCHI✔❌:             *pbFixFB2 = 0;
    // INCHI✔❌:             *pbStdFormat = 0;
    // INCHI✔❌:             got = 1;
    // INCHI✔❌:         }
    // INCHI✔❌:         else if (!inchi_stricmp(pArg, "SPXYZOFF"))
    // INCHI✔❌:         {
    // INCHI✔❌:             *pbAddPhosphineStereo = 0;
    // INCHI✔❌:             *pbStdFormat = 0;
    // INCHI✔❌:             got = 1;
    // INCHI✔❌:         }
    // INCHI✔❌:         else if (!inchi_stricmp(pArg, "SASXYZOFF"))
    // INCHI✔❌:         {
    // INCHI✔❌:             *pbAddArsineStereo = 0;
    // INCHI✔❌:             *pbStdFormat = 0;
    // INCHI✔❌:             got = 1;
    // INCHI✔❌:         }
    // INCHI✔❌:         else if (!inchi_stricmp(pArg, "Polymers105+"))
    // INCHI✔❌:         {
    // INCHI✔❌:             *pbPolymers = POLYMERS_LEGACY_PLUS;
    // INCHI✔❌:             *pbStdFormat = 0;
    // INCHI✔❌:             got = 1;
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return got;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: set_common_options_by_parg

    let mut got = 0;
    if option_literal(heap, pArg, "INPAUX")? {
        if ip.nInputType == tagInputType_INPUT_NONE {
            ip.nInputType = tagInputType_INPUT_INCHI_PLAIN;
        }
        got = 1;
    } else if option_prefix(heap, pArg, "START:")? {
        ip.first_struct_number = source_strtol_long(heap, pArg.offset(6)?)?;
        got = 1;
    } else if option_prefix(heap, pArg, "END:")? {
        ip.last_struct_number = source_strtol_long(heap, pArg.offset(4)?)?;
        got = 1;
    } else if option_prefix(heap, pArg, "RECORD:")? {
        let number = source_strtol_long(heap, pArg.offset(7)?)?;
        ip.first_struct_number = number;
        ip.last_struct_number = number;
        got = 1;
    } else if option_literal(heap, pArg, "NOLABELS")? {
        options.no_struct_labels = 1;
        got = 1;
    } else if option_literal(heap, pArg, "SAVEOPT")? {
        options.inchi_output_options |= INCHI_OUT_SAVEOPT as i32;
        got = 1;
    } else if option_literal(heap, pArg, "AUXNONE")? {
        options.inchi_output_options |= INCHI_OUT_NO_AUX_INFO as i32;
        options.inchi_output_options &= !(INCHI_OUT_SHORT_AUX_INFO as i32);
        got = 1;
    } else if option_literal(heap, pArg, "MISMATCHISERROR")? {
        options.inchi_output_options2 |= INCHI_OUT_MISMATCH_AS_ERROR as i32;
        got = 1;
    } else if option_literal(heap, pArg, "OUTERRINCHI")? {
        options.inchi_output_options2 |= INCHI_OUT_INCHI_GEN_ERROR as i32;
        got = 1;
    } else if option_literal(heap, pArg, "Key")? {
        options.hash_key = 1;
        got = 1;
    } else if option_literal(heap, pArg, "XHash1")? {
        options.hash_xtra1 = 1;
        got = 1;
    } else if option_literal(heap, pArg, "XHash2")? {
        options.hash_xtra2 = 1;
        got = 1;
    } else if option_literal(heap, pArg, "SNON")? {
        options.ver1_default_mode &= !(REQ_MODE_STEREO as INCHI_MODE);
        options.mode &= !((REQ_MODE_RACEMIC_STEREO
            | REQ_MODE_RELATIVE_STEREO
            | REQ_MODE_CHIR_FLG_STEREO) as i32);
        got = 1;
    } else if option_literal(heap, pArg, "NEWPSOFF")? {
        options.pointed_edge_stereo = 0;
        got = 1;
    } else if option_literal(heap, pArg, "DONOTADDH")? {
        options.do_not_add_h = 1;
        got = 1;
    } else if option_literal(heap, pArg, "LooseTSACheck")? {
        options.loose_tsa_check = 1;
        got = 1;
    } else if option_literal(heap, pArg, "SREL")? {
        options.mode &= !(REQ_MODE_RACEMIC_STEREO as i32);
        options.mode &= !(REQ_MODE_CHIR_FLG_STEREO as i32);
        options.mode |= (REQ_MODE_RELATIVE_STEREO | REQ_MODE_STEREO) as i32;
        options.std_format = 0;
        got = 1;
    } else if option_literal(heap, pArg, "SRAC")? {
        options.mode &= !(REQ_MODE_RELATIVE_STEREO as i32);
        options.mode &= !(REQ_MODE_CHIR_FLG_STEREO as i32);
        options.mode |= (REQ_MODE_RACEMIC_STEREO | REQ_MODE_STEREO) as i32;
        options.std_format = 0;
        got = 1;
    } else if option_literal(heap, pArg, "SUCF")? {
        options.mode &= !(REQ_MODE_RELATIVE_STEREO as i32);
        options.mode &= !(REQ_MODE_RACEMIC_STEREO as i32);
        options.mode |= (REQ_MODE_CHIR_FLG_STEREO | REQ_MODE_STEREO) as i32;
        options.std_format = 0;
        got = 1;
    } else if option_literal(heap, pArg, "ChiralFlagON")? {
        options.forced_chiral_flag &= !(FLAG_SET_INP_AT_NONCHIRAL as i32);
        options.forced_chiral_flag |= FLAG_SET_INP_AT_CHIRAL as i32;
        got = 1;
    } else if option_literal(heap, pArg, "ChiralFlagOFF")? {
        options.forced_chiral_flag &= !(FLAG_SET_INP_AT_CHIRAL as i32);
        options.forced_chiral_flag |= FLAG_SET_INP_AT_NONCHIRAL as i32;
        got = 1;
    } else if option_literal(heap, pArg, "SUU")? {
        options.ver1_default_mode &=
            !((REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU) as INCHI_MODE);
        options.std_format = 0;
        got = 1;
    } else if option_literal(heap, pArg, "SLUUD")? {
        options.ver1_default_mode |= REQ_MODE_DIFF_UU_STEREO as INCHI_MODE;
        options.std_format = 0;
        got = 1;
    } else if option_literal(heap, pArg, "FIXEDH")? {
        options.ver1_default_mode |= REQ_MODE_BASIC as INCHI_MODE;
        options.std_format = 0;
        got = 1;
    } else if option_literal(heap, pArg, "RECMET")? {
        options.reconnect_coord = 1;
        options.std_format = 0;
        got = 1;
    } else if option_literal(heap, pArg, "KET")? {
        options.keto_enol_taut = 1;
        options.std_format = 0;
        got = 1;
    } else if option_literal(heap, pArg, "15T")? {
        options.taut_15_non_ring = 1;
        options.std_format = 0;
        got = 1;
    } else if option_literal(heap, pArg, "PT_22_00")? {
        options.pt_22_00_taut = 1;
        options.std_format = 0;
    } else if option_literal(heap, pArg, "PT_16_00")? {
        options.pt_16_00_taut = 1;
        options.std_format = 0;
    } else if option_literal(heap, pArg, "PT_06_00")? {
        options.pt_06_00_taut = 1;
        options.std_format = 0;
    } else if option_literal(heap, pArg, "PT_39_00")? {
        options.pt_39_00_taut = 1;
        options.std_format = 0;
    } else if option_literal(heap, pArg, "PT_13_00")? {
        options.pt_13_00_taut = 1;
        options.std_format = 0;
    } else if option_literal(heap, pArg, "PT_18_00")? {
        options.pt_18_00_taut = 1;
        options.std_format = 0;
    } else if option_literal(heap, pArg, "LargeMolecules")? {
        options.large_molecules = 1;
        got = 1;
    } else if option_literal(heap, pArg, "Polymers")? {
        options.polymers = POLYMERS_MODERN as i32;
        got = 1;
    } else if option_literal(heap, pArg, "Polymers105")? {
        options.polymers = POLYMERS_LEGACY as i32;
        got = 1;
    } else if option_literal(heap, pArg, "NPZz")? {
        options.np_zz = 1;
        got = 1;
    } else if option_literal(heap, pArg, "NoWarnings")? {
        options.no_warnings = 1;
        got = 1;
    } else if option_literal(heap, pArg, "MergeHash")? {
        options.merge_hash = 1;
        got = 1;
    } else if option_literal(heap, pArg, "NoInChI")? || option_literal(heap, pArg, "HideInChI")? {
        options.hide_inchi = 1;
        got = 1;
    } else if option_literal(heap, pArg, "FoldCRU")? || option_literal(heap, pArg, "FoldSRU")? {
        options.fold_polymer_sru = 1;
        got = 1;
    } else if option_prefix(heap, pArg, "FrameShift:")? {
        let word = heap.allocate_model_storage(vec![0_i8; 256])?;
        mystrncpy(heap, word, pArg.offset(11)?, 256)?;
        let mut length = 0;
        lrtrim(heap, word, Some(&mut length))?;
        if length != 0 {
            if option_literal(heap, word.as_const(), "None")? {
                options.frame_shift_scheme = tagFrameShifScheme_FSS_NONE as i32;
            } else if option_literal(heap, word.as_const(), "Cyclize")? {
                options.frame_shift_scheme = tagFrameShifScheme_FSS_STARS_CYCLED as i32;
            } else if option_literal(heap, word.as_const(), "MoveStars")? {
                options.frame_shift_scheme = tagFrameShifScheme_FSS_STARS_OPENED as i32;
            } else if option_literal(heap, word.as_const(), "MoveBrackets")? {
                options.frame_shift_scheme = tagFrameShifScheme_FSS_STARS_ENDS_OPENED as i32;
            }
        } else {
            options.frame_shift_scheme = tagFrameShifScheme_FSS_STARS_CYCLED as i32;
        }
        heap.free(word)?;
        got = 1;
    } else if option_literal(heap, pArg, "NoFrameShift")? {
        options.frame_shift_scheme = tagFrameShifScheme_FSS_NONE as i32;
        got = 1;
    } else if option_literal(heap, pArg, "NoEdits")? {
        options.fold_polymer_sru = 0;
        options.frame_shift_scheme = tagFrameShifScheme_FSS_NONE as i32;
        got = 1;
    } else if option_literal(heap, pArg, "SATZZ")? {
        options.stereo_at_zz = 1;
        got = 1;
    }

    if got == 0 && developer_options != 0 {
        if option_literal(heap, pArg, "PGO")? {
            ip.bSaveAllGoodStructsAsProblem = 1;
            got = 1;
        } else if option_literal(heap, pArg, "FilterSS")? {
            ip.bFilterSS = 1;
            got = 1;
        } else if option_literal(heap, pArg, "InvFilterSS")? {
            ip.bFilterSS = -1;
            got = 1;
        }
        if option_literal(heap, pArg, "FNUDOFF")? {
            ip.bFixNonUniformDraw = 0;
            options.std_format = 0;
            got = 1;
        } else if option_literal(heap, pArg, "FixSp3bugOFF")?
            || option_literal(heap, pArg, "FBOFF")?
        {
            options.fix_sp3_bug = 0;
            options.std_format = 0;
            got = 1;
        } else if option_literal(heap, pArg, "FB2OFF")? {
            options.fix_fb2 = 0;
            options.std_format = 0;
            got = 1;
        } else if option_literal(heap, pArg, "SPXYZOFF")? {
            options.add_phosphine_stereo = 0;
            options.std_format = 0;
            got = 1;
        } else if option_literal(heap, pArg, "SASXYZOFF")? {
            options.add_arsine_stereo = 0;
            options.std_format = 0;
            got = 1;
        } else if option_literal(heap, pArg, "Polymers105+")? {
            options.polymers = POLYMERS_LEGACY_PLUS as i32;
            options.std_format = 0;
            got = 1;
        }
    }
    Ok(got)
}

fn source_c_string(
    heap: &SourceHeap,
    pointer: SourceConstPointer<i8>,
) -> Result<Vec<i8>, SourceHeapError> {
    let bytes = heap.slice(pointer)?;
    let nul = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    Ok(bytes[..=nul].to_vec())
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn ReadCommandLineParms(
    heap: &mut SourceHeap,
    argc: i32,
    argv: &[SourceConstPointer<i8>],
    ip: &mut INPUT_PARMS,
    szSdfDataValue: SourceMutPointer<i8>,
    ulDisplTime: &mut u64,
    _bReleaseVersion: i32,
    mut log_file: Option<&mut INCHI_IOSTREAM>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c:602 ReadCommandLineParms
    // INCHI✔❌: int ReadCommandLineParms(int argc,
    // INCHI✔❌:     const char* argv[],
    // INCHI✔❌:     INPUT_PARMS* ip,
    // INCHI✔❌:     char* szSdfDataValue,
    // INCHI✔❌:     unsigned long* ulDisplTime,
    // INCHI✔❌:     int bReleaseVersion,
    // INCHI✔❌:     INCHI_IOSTREAM* log_file)
    // INCHI✔❌: {
    // INCHI✔❌: #if (BUILD_WITH_ENG_OPTIONS==1)
    // INCHI✔❌:     const int developer_options = 1;
    // INCHI✔❌: #else
    // INCHI✔❌:     const int developer_options = 0;
    // INCHI✔❌: #endif
    // INCHI✔❌:     const char* q;
    // INCHI✔❌:     const char* ext[MAX_NUM_PATHS + 1];
    // INCHI✔❌:     const char* pArg;
    // INCHI✔❌:     char szNameSuffix[32], szOutputPath[512];
    // INCHI✔❌:
    // INCHI✔❌:     unsigned long ul;
    // INCHI✔❌:
    // INCHI✔❌:     double t = 0;
    // INCHI✔❌:
    // INCHI✔❌:     int bVer1Options = 1;
    // INCHI✔❌:     int nMode = 0;
    // INCHI✔❌:     int nReleaseMode = nMode | (REQ_MODE_BASIC | REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_STEREO);
    // INCHI✔❌:     INCHI_MODE bVer1DefaultMode = VER103_DEFAULT_MODE;
    // INCHI✔❌:
    // INCHI✔❌:     int bNameSuffix;
    // INCHI✔❌:     int bOutputPath;
    // INCHI✔❌:     int bMergeAllInputStructures;
    // INCHI✔❌:     int bForcedChiralFlag = 0;
    // INCHI✔❌:
    // INCHI✔❌:     int bDisconnectSalts = (DISCONNECT_SALTS == 1);
    // INCHI✔❌:     int bDoNotAddH = 0;
    // INCHI✔❌:     int bRecognizedOption;
    // INCHI✔❌:     int bTgFlagVariableProtons = 1;
    // INCHI✔❌:     int bTgFlagHardAddRenProtons = 1;
    // INCHI✔❌:     int bReconnectCoord = (RECONNECT_METALS == 1);
    // INCHI✔❌:     int bDisconnectCoord = (DISCONNECT_METALS == 1);
    // INCHI✔❌:     int bDisconnectCoordChkVal = (CHECK_METAL_VALENCE == 1);
    // INCHI✔❌:     int bMovePositiveCharges = (MOVE_CHARGES == 1);
    // INCHI✔❌:     int bAcidTautomerism = (DISCONNECT_SALTS == 1) ? (TEST_REMOVE_S_ATOMS == 1 ? 2 : 1) : 0;
    // INCHI✔❌:     int bUnchargedAcidTaut = (CHARGED_SALTS_ONLY == 0);
    // INCHI✔❌:     int bMergeSaltTGroups = (DISCONNECT_SALTS == 1);
    // INCHI✔❌: #if ( MIN_SB_RING_SIZE > 0 )
    // INCHI✔❌:     int nMinDbRinSize = MIN_SB_RING_SIZE, mdbr = 0;
    // INCHI✔❌: #endif
    // INCHI✔❌: #ifdef STEREO_WEDGE_ONLY
    // INCHI✔❌:     int bPointedEdgeStereo = STEREO_WEDGE_ONLY; /* NEWPS TG_FLAG_POINTED_EDGE_STEREO */
    // INCHI✔❌: #endif
    // INCHI✔❌: #if ( FIX_ADJ_RAD == 1 )
    // INCHI✔❌:     int bFixAdjacentRad = 0;
    // INCHI✔❌: #endif
    // INCHI✔❌:     int bAddPhosphineStereo = 1;
    // INCHI✔❌:     int bAddArsineStereo = 1;
    // INCHI✔❌:     int bFixSp3bug = 1;
    // INCHI✔❌:     int bFixFB2 = 1;
    // INCHI✔❌:     int bKetoEnolTaut = 0;
    // INCHI✔❌:     int b15TautNonRing = 0;
    // INCHI✔❌:     int bPT_22_00_Taut = 0;
    // INCHI✔❌:     int bPT_16_00_Taut = 0;
    // INCHI✔❌:     int bPT_06_00_Taut = 0;
    // INCHI✔❌:     int bPT_39_00_Taut = 0;
    // INCHI✔❌:     int bPT_13_00_Taut = 0;
    // INCHI✔❌:     int bPT_18_00_Taut = 0;
    // INCHI✔❌:     int bStdFormat = 1;
    // INCHI✔❌:     int bHashKey = 0;
    // INCHI✔❌:     int bHashXtra1 = 0;
    // INCHI✔❌:     int bHashXtra2 = 0;
    // INCHI✔❌:     int bLargeMolecules = 0;
    // INCHI✔❌:     int bPolymers = POLYMERS_NO;
    // INCHI✔❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌:     int bFoldPolymerSRU = 0;
    // INCHI✔❌:     int bFrameShiftScheme = FSS_STARS_CYCLED;
    // INCHI✔❌: #else
    // INCHI✔❌:     int bFoldPolymerSRU = 0;
    // INCHI✔❌:     int bFrameShiftScheme = FSS_STARS_CYCLED;
    // INCHI✔❌: #endif
    // INCHI✔❌:     int bLooseTSACheck = 0;
    // INCHI✔❌:     int bStereoAtZz = 0;
    // INCHI✔❌:     int bNPZz = 0;
    // INCHI✔❌:     int bNoWarnings = 0;
    // INCHI✔❌:     int bMergeHash = 0;
    // INCHI✔❌:     int bHideInChI = 0;
    // INCHI✔❌:
    // INCHI✔❌:     int bOutputStyle = INCHI_OUT_PLAIN_TEXT;
    // INCHI✔❌:     int bDisplay = 0;
    // INCHI✔❌:     int bNoStructLabels = 0;
    // INCHI✔❌:     int bOutputMolfileOnly = 0;
    // INCHI✔❌:     int bOutputMolfileDT = 0;
    // INCHI✔❌:     int bOutputMolfileSplit = 0;
    // INCHI✔❌:     int bDisplayCompositeResults = 0; /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔❌:     int nFontSize = -9; /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔❌:     int bINChIOutputOptions2 = 0;
    // INCHI✔❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌:     int is_gui = 1;
    // INCHI✔❌:     int bINChIOutputOptions = INCHI_OUT_EMBED_REC; /* embed reconnected & output full aux info */
    // INCHI✔❌:     int bCompareComponents = CMP_COMPONENTS;
    // INCHI✔❌: #else
    // INCHI✔❌:     int is_gui = 0;
    // INCHI✔❌:     int bINChIOutputOptions = ((EMBED_REC_METALS_INCHI == 1) ? INCHI_OUT_EMBED_REC : 0);
    // INCHI✔❌:     int bCompareComponents = 0;
    // INCHI✔❌: #endif
    // INCHI✔❌: #if ( READ_INCHI_STRING == 1 )
    // INCHI✔❌:     int bDisplayIfRestoreWarnings = 0;
    // INCHI✔❌: #endif
    // INCHI✔❌: #if ( BUILD_WITH_AMI == 1 ) && ( OUTPUT_FILE_EXT == 1 )
    // INCHI✔❌:     int numOutNameExt;
    // INCHI✔❌:     char szOutNameExt[3][128];
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌:     int i, k, c, got;
    // INCHI✔❌:     int timeout_set_warning = 0;
    // INCHI✔❌:     int timeout_set_error = 0;
    // INCHI✔❌:
    // INCHI✔❌:     ext[0] = ".mol";
    // INCHI✔❌:     ext[1] = bVer1Options ? ".txt" : ".ich";
    // INCHI✔❌:     ext[2] = ".log";
    // INCHI✔❌:     ext[3] = ".prb";
    // INCHI✔❌:     ext[4] = "";
    // INCHI✔❌:
    // INCHI✔❌:     /*  Init table of parameters */
    // INCHI✔❌:     memset(ip, 0, sizeof(*ip)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:
    // INCHI✔❌:     /* Default are StdInChI generation options */
    // INCHI✔❌:     bVer1DefaultMode &= ~REQ_MODE_BASIC; /* "FIXEDH - OFF" */
    // INCHI✔❌:     bReconnectCoord = 0;                /* "RECMET - OFF" */
    // INCHI✔❌:     bPointedEdgeStereo = 1;                /* "NEWPS"        */
    // INCHI✔❌:     ip->bFixNonUniformDraw = 1;                /* "FNUD"         */
    // INCHI✔❌:
    // INCHI✔❌: #ifndef COMPILE_ANSI_ONLY
    // INCHI✔❌:     strcpy(ip->tdp.ReqShownFoundTxt[ilSHOWN], "Shown");
    // INCHI✔❌:     ip->dp.sdp.tdp = &ip->tdp;
    // INCHI✔❌:     ip->dp.pdp = &ip->pdp;
    // INCHI✔❌: #endif
    // INCHI✔❌:     memset(szNameSuffix, 0, sizeof(szNameSuffix)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:     bNameSuffix = 0;
    // INCHI✔❌:     memset(szOutputPath, 0, sizeof(szOutputPath)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:     bOutputPath = 0;
    // INCHI✔❌: #if( OUTPUT_FILE_EXT == 1 )
    // INCHI✔❌:     memset(szOutNameExt, 0, sizeof(szOutNameExt)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:     numOutNameExt = 0;
    // INCHI✔❌: #endif
    // INCHI✔❌:     bMergeAllInputStructures = 0;
    // INCHI✔❌: #ifdef TARGET_API_LIB
    // INCHI✔❌:     ip->msec_MaxTime = 0;      /*  milliseconds, default in libinchi: unlimited */
    // INCHI✔❌: #else
    // INCHI✔❌:     ip->msec_MaxTime = 60000;  /*  milliseconds, default: 60 sec */
    // INCHI✔❌: #endif
    // INCHI✔❌:     * ulDisplTime = 0;
    // INCHI✔❌:
    // INCHI✔❌:     if (bReleaseVersion)
    // INCHI✔❌:     {
    // INCHI✔❌:         /*  normal */
    // INCHI✔❌:         ip->bAbcNumbers = 0;
    // INCHI✔❌:         ip->bCtPredecessors = 0;
    // INCHI✔❌:     }
    // INCHI✔❌:     else
    // INCHI✔❌:     {
    // INCHI✔❌:         nReleaseMode = 0;
    // INCHI✔❌:     }
    // INCHI✔❌:     if (bVer1Options)
    // INCHI✔❌:     {
    // INCHI✔❌:         bNameSuffix = 1;
    // INCHI✔❌:         szNameSuffix[0] = '\0';
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌: #if ( ALLOW_SUBSTRUCTURE_FILTERING== 1 )
    // INCHI✔❌:     ip->bFilterSS = 0;
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌:     /* Analyze command line switches */
    // INCHI✔❌:     for (i = 1; i < argc; i++)
    // INCHI✔❌:     {
    // INCHI✔❌:
    // INCHI✔❌:         if (is_gui && INCHI_OPTION_PREFX == argv[i][0] && INCHI_OPTION_PREFX != argv[i][1])
    // INCHI✔❌:         {
    // INCHI✔❌:             /* Parsing TARGET_LIB_FOR_WINCHI GUI options (and v. 0.9xx Beta as well)  */
    // INCHI✔❌:             pArg = argv[i] + 1;
    // INCHI✔❌:             got = set_common_options_by_parg(pArg, developer_options, ip, &bVer1DefaultMode, &nMode,
    // INCHI✔❌:                 &bINChIOutputOptions, &bINChIOutputOptions2,
    // INCHI✔❌:                 &bStdFormat, &bHashKey, &bHashXtra1, &bHashXtra2,
    // INCHI✔❌:                 &bFixSp3bug, &bFixFB2,
    // INCHI✔❌:                 &bAddPhosphineStereo, &bAddArsineStereo,
    // INCHI✔❌:                 &bNoStructLabels, &bPointedEdgeStereo,
    // INCHI✔❌:                 &bDoNotAddH, &bForcedChiralFlag, &bReconnectCoord,
    // INCHI✔❌:                 &bKetoEnolTaut,
    // INCHI✔❌:                 &b15TautNonRing,
    // INCHI✔❌:                 &bPT_06_00_Taut, &bPT_13_00_Taut, &bPT_16_00_Taut,
    // INCHI✔❌:                 &bPT_18_00_Taut, &bPT_22_00_Taut, &bPT_39_00_Taut,
    // INCHI✔❌:                 &bLooseTSACheck,
    // INCHI✔❌:                 &bLargeMolecules, &bPolymers,
    // INCHI✔❌:                 &bFoldPolymerSRU, &bFrameShiftScheme,
    // INCHI✔❌:                 &bStereoAtZz, &bNPZz,
    // INCHI✔❌:                 &bNoWarnings, &bMergeHash, &bHideInChI);
    // INCHI✔❌:             if (got)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (INPUT_NONE == ip->nInputType &&
    // INCHI✔❌:                 (!inchi_memicmp(pArg, "SDF", 3)) &&
    // INCHI✔❌:                 (pArg[3] == ':'))
    // INCHI✔❌:             {
    // INCHI✔❌:                 k = 0;
    // INCHI✔❌:                 mystrncpy(ip->szSdfDataHeader, pArg + 4, MAX_SDF_HEADER + 1);
    // INCHI✔❌:                 lrtrim(ip->szSdfDataHeader, &k);
    // INCHI✔❌:                 if (k)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     ip->pSdfLabel = ip->szSdfDataHeader;
    // INCHI✔❌:                     ip->pSdfValue = szSdfDataValue;
    // INCHI✔❌:                     ip->nInputType = INPUT_SDFILE;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                     ip->pSdfLabel = NULL;
    // INCHI✔❌:                     ip->pSdfValue = NULL;
    // INCHI✔❌:                     ip->nInputType = INPUT_MOLFILE;
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (INPUT_NONE == ip->nInputType && !inchi_stricmp(pArg, "MOL"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->nInputType = INPUT_MOLFILE;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (INPUT_NONE == ip->nInputType && !inchi_stricmp(pArg, "SDF"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->nInputType = INPUT_MOLFILE;
    // INCHI✔❌:             }
    // INCHI✔❌:             /*--- Output options ---*/
    // INCHI✔❌: #if ( !defined(TARGET_API_LIB) && !defined(TARGET_LIB_FOR_WINCHI) )
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "Tabbed") || !inchi_stricmp(pArg, "Tab"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bOutputStyle |= INCHI_OUT_TABBED_OUTPUT;
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌:             /* Removed condition
    // INCHI✔❌:                 #if ( defined(BUILD_WITH_ENG_OPTIONS) || defined(TARGET_LIB_FOR_WINCHI) )
    // INCHI✔❌:                 which is always true, as we already are under condition
    // INCHI✔❌:                 'if ( is_gui && ...)' and is_gui==1 means TARGET_LIB_FOR_WINCHI
    // INCHI✔❌:             */
    // INCHI✔❌:             /* #if ( defined(BUILD_WITH_ENG_OPTIONS) || defined(TARGET_LIB_FOR_WINCHI) )  */
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "SDFID"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bGetSdfileId = 1;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "PLAIN"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bOutputStyle |= INCHI_OUT_PLAIN_TEXT;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "ANNPLAIN"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bOutputStyle |= INCHI_OUT_PLAIN_TEXT_COMMENTS;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "AUXINFO:", 8) && isdigit(UCINT pArg[8]))
    // INCHI✔❌:             {
    // INCHI✔❌:                 k = strtol(pArg + 8, NULL, 10);
    // INCHI✔❌:                 if (k == 0)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     bINChIOutputOptions |= INCHI_OUT_NO_AUX_INFO;  /* no aux info */
    // INCHI✔❌:                     bINChIOutputOptions &= ~INCHI_OUT_SHORT_AUX_INFO;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else if (k == 1)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     bINChIOutputOptions &= ~(INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO); /* include full aux info */
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else if (k == 2)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     bINChIOutputOptions &= ~INCHI_OUT_NO_AUX_INFO; /* include short aux info */
    // INCHI✔❌:                     bINChIOutputOptions |= INCHI_OUT_SHORT_AUX_INFO;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                     bINChIOutputOptions = k;  /* override everything */
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "MERGE"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bMergeAllInputStructures = 1;
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "PGO"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bSaveAllGoodStructsAsProblem = 1;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "DCR"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bDisplayCompositeResults = 1;
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "DSB"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 nMode |= REQ_MODE_NO_ALT_SBONDS;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "NOHDR"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bNoStructLabels = 1;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "NoVarH"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bTgFlagVariableProtons = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             /*--- (hidden) Old structure-perception and InChI creation options ---*/
    // INCHI✔❌:             /*--- (engineering) Old structure-perception and InChI creation options ---*/
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "NOUUSB"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 nMode |= REQ_MODE_SB_IGN_ALL_UU;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "NOUUSC"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 nMode |= REQ_MODE_SC_IGN_ALL_UU;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌: #if ( FIX_ADJ_RAD == 1 )
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "FixRad"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bFixAdjacentRad = 1;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌: #if ( RENUMBER_ATOMS_AND_RECALC_V106 == 1 )
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "TestRenum") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bRenumber = 1;
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌: #if ( UNDERIVATIZE == 1 )
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "DoDRV"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bUnderivatize = 1;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌: #if( UNDERIVATIZE_REPORT == 1 )
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "DoDrvReport"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bUnderivatize = 3;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌: #endif
    // INCHI✔❌: #if ( RING2CHAIN == 1 )
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "DoR2C"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bRing2Chain = 1;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌: #if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "DoneOnly"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bIgnoreUnchanged = 1;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "NoADP"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bTgFlagHardAddRenProtons = 0;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "DISCONSALT:", 11))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bDisconnectSalts = (0 != strtol(pArg + 11, NULL, 10));
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "DISCONMETAL:", 12))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bDisconnectCoord = (0 != strtol(pArg + 12, NULL, 10));
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "RECONMETAL:", 11))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bReconnectCoord = (0 != strtol(pArg + 11, NULL, 10));
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "DISCONMETALCHKVAL:", 18))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bDisconnectCoordChkVal = (0 != strtol(pArg + 18, NULL, 10));
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "MOVEPOS:", 8))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bMovePositiveCharges = (0 != strtol(pArg + 8, NULL, 10));
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "MERGESALTTG:", 12))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bMergeSaltTGroups = (0 != strtol(pArg + 12, NULL, 10));
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "UNCHARGEDACIDS:", 15))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bUnchargedAcidTaut = (0 != strtol(pArg + 15, NULL, 16));
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "ACIDTAUT:", 9))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bAcidTautomerism = c = (int)strtol(pArg + 9, NULL, 10);
    // INCHI✔❌:                 if (0 <= c && c <= 2)  bAcidTautomerism = c;
    // INCHI✔❌:                 /*else bNotRecognized = 2*bReleaseVersion;*/
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             /*--- (hidden) Old output and other options ---*/
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "O:", 2))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bNameSuffix = 1;
    // INCHI✔❌:                 strncpy(szNameSuffix, pArg + 2, sizeof(szNameSuffix) - 1);
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "OP:", 3))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bOutputPath = 1;
    // INCHI✔❌:                 strncpy(szOutputPath, pArg + 3, sizeof(szOutputPath) - 1);
    // INCHI✔❌:                 }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "ALT"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bAbcNumbers = 1;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "SCT"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bCtPredecessors = 1;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "CMP"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bCompareComponents = CMP_COMPONENTS;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "CMPNONISO"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bCompareComponents = CMP_COMPONENTS | CMP_COMPONENTS_NONISO;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "PW"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bSaveWarningStructsAsProblem = 1;
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "RSB:", 4) && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 mdbr = (int)strtol(pArg + 4, NULL, 10);
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "DISCONSALT:", 11) && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bDisconnectSalts = (0 != strtol(pArg + 11, NULL, 10));
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "DISCONMETAL:", 12) && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bDisconnectCoord = (0 != strtol(pArg + 12, NULL, 10));
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "RECONMETAL:", 11) && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bReconnectCoord = (0 != strtol(pArg + 11, NULL, 10));
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "DISCONMETALCHKVAL:", 18) && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bDisconnectCoordChkVal = (0 != strtol(pArg + 18, NULL, 10));
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "MOVEPOS:", 8) && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bMovePositiveCharges = (0 != strtol(pArg + 8, NULL, 10));
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "MERGESALTTG:", 12) && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bMergeSaltTGroups = (0 != strtol(pArg + 12, NULL, 10));
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "UNCHARGEDACIDS:", 15) && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bUnchargedAcidTaut = (0 != strtol(pArg + 15, NULL, 16));;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "ACIDTAUT:", 9) && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bAcidTautomerism = c = (int)strtol(pArg + 9, NULL, 10);
    // INCHI✔❌:                 if (0 <= c && c <= 2)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     bAcidTautomerism = c;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 /*else bNotRecognized = 2*bReleaseVersion;*/
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 /*for ( k = 0; c=pArg[k]; k ++ )*/
    // INCHI✔❌:                 k = 0;
    // INCHI✔❌:                 c = pArg[k]; /* prohibit multiple option concatenations, strict syntax check 2008-11-05 DT  */
    // INCHI✔❌:                 {
    // INCHI✔❌:                     c = toupper(c);
    // INCHI✔❌:                     switch (c)
    // INCHI✔❌:                     {
    // INCHI✔❌:                     case 'D':
    // INCHI✔❌:                         bDisplay |= 1;
    // INCHI✔❌:                         if ((pArg[k + 1] == 'C' || pArg[k + 1] == 'c') && !pArg[k + 2])
    // INCHI✔❌:                         {
    // INCHI✔❌:                             bDisplay |= 1;
    // INCHI✔❌:                             k++;
    // INCHI✔❌:                             ip->bDisplayEachComponentINChI = 1;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         else if (!pArg[k + 1])
    // INCHI✔❌:                         {
    // INCHI✔❌:                             bDisplay |= 1;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         break;
    // INCHI✔❌:                     case 'W':
    // INCHI✔❌:                         if (pArg[k + 1] == 'D')
    // INCHI✔❌:                         {
    // INCHI✔❌:                             /* restore Display Time functionality */
    // INCHI✔❌:                             c = 'D';
    // INCHI✔❌:                             k++;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         t = strtod(pArg + k + 1, (char**)&q); /*  cast deliberately discards 'const' qualifier */
    // INCHI✔❌:                         if ((q > pArg + k + 1 && errno == ERANGE) || t < 0.0 || t * 1000.0 >(double)ULONG_MAX) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:                         {
    // INCHI✔❌:                             ul = 0;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         else
    // INCHI✔❌:                         {
    // INCHI✔❌:                             ul = (unsigned long)(t * 1000.0);
    // INCHI✔❌:                         }
    // INCHI✔❌:                         if ( /*q > pArg+k &&*/ !*q)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             k = q - pArg - 1; /* k will be incremented by the for() cycle */
    // INCHI✔❌:                             switch (c)
    // INCHI✔❌:                             {
    // INCHI✔❌:                             case 'D':
    // INCHI✔❌:                                 *ulDisplTime = ul;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 'W':
    // INCHI✔❌:                                 ip->msec_MaxTime = ul;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             }
    // INCHI✔❌:                         }
    // INCHI✔❌:                         break;
    // INCHI✔❌:                     case 'F':
    // INCHI✔❌:                         c = (int)strtol(pArg + k + 1, (char**)&q, 10); /*  cast deliberately discards 'const' qualifier */
    // INCHI✔❌:                         if (q > pArg + k && !*q)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             k = q - pArg - 1;
    // INCHI✔❌:                             if (abs(c) > 5)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 nFontSize = -c;  /* font size 5 or less is too small */
    // INCHI✔❌:                             }
    // INCHI✔❌:                         }
    // INCHI✔❌:                         break;
    // INCHI✔❌:                     default:
    // INCHI✔❌:                         if (!pArg[k + 1])
    // INCHI✔❌:                         {
    // INCHI✔❌:                             switch (c)
    // INCHI✔❌:                             {
    // INCHI✔❌:                             case 'B':
    // INCHI✔❌:                                 nMode |= REQ_MODE_BASIC;
    // INCHI✔❌:                                 nReleaseMode = 0;
    // INCHI✔❌:                                 /*bNotRecognized = bReleaseVersion;*/
    // INCHI✔❌:                                 bStdFormat = 0;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 'T':
    // INCHI✔❌:                                 nMode |= REQ_MODE_TAUT;
    // INCHI✔❌:                                 nReleaseMode = 0;
    // INCHI✔❌:                                 /*bNotRecognized = bReleaseVersion;*/
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 'I':
    // INCHI✔❌:                                 nMode |= REQ_MODE_ISO;
    // INCHI✔❌:                                 nReleaseMode = 0;
    // INCHI✔❌:                                 /*bNotRecognized = bReleaseVersion;*/
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 'N':
    // INCHI✔❌:                                 nMode |= REQ_MODE_NON_ISO;
    // INCHI✔❌:                                 nReleaseMode = 0;
    // INCHI✔❌:                                 bStdFormat = 0;
    // INCHI✔❌:                                 /*bNotRecognized = bReleaseVersion;*/
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 'S':
    // INCHI✔❌:                                 nMode |= REQ_MODE_STEREO;
    // INCHI✔❌:                                 nReleaseMode = 0;
    // INCHI✔❌:                                 /*bNotRecognized = bReleaseVersion;*/
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 'E':
    // INCHI✔❌:                                 if (nReleaseMode & REQ_MODE_STEREO)
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     nReleaseMode ^= REQ_MODE_STEREO;
    // INCHI✔❌:                                     bStdFormat = 0;
    // INCHI✔❌:                                 }
    // INCHI✔❌:                                 break;
    // INCHI✔❌:
    // INCHI✔❌: #ifndef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌:                             default:
    // INCHI✔❌:                                 inchi_ios_eprint(log_file, "Unrecognized optionQ1: \"%c\".\n", c);
    // INCHI✔❌:
    // INCHI✔❌: #endif
    // INCHI✔❌:                             }
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌: #ifndef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌:                         else
    // INCHI✔❌:                         {
    // INCHI✔❌:                             inchi_ios_eprint(log_file, "Unrecognized optionQ2: \"%c\".\n", c);
    // INCHI✔❌:                         }
    // INCHI✔❌: #endif
    // INCHI✔❌:                     }
    // INCHI✔❌:                     /*
    // INCHI✔❌:                         if ( bNotRecognized && bNotRecognized == bReleaseVersion ) {
    // INCHI✔❌:                         inchi_ios_eprint(stderr, "Unrecognized option: \"%c\".\n", c);
    // INCHI✔❌:                         bNotRecognized = 0;
    // INCHI✔❌:                         }
    // INCHI✔❌:                     */
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             /*
    // INCHI✔❌:                 if ( bNotRecognized && bNotRecognized == 2*bReleaseVersion )
    // INCHI✔❌:                 {
    // INCHI✔❌:                 inchi_ios_eprint(stderr, "Unrecognized option: \"%s\".\n", argv[i]);
    // INCHI✔❌:                 bNotRecognized = 0;
    // INCHI✔❌:                 }
    // INCHI✔❌:             */
    // INCHI✔❌:
    // INCHI✔❌:         } /* eof Parsing TARGET_LIB_FOR_WINCHI GUI options (and v. 0.9xx Beta as well)  */
    // INCHI✔❌:
    // INCHI✔❌:         else if ((bVer1Options & 1) && INCHI_OPTION_PREFX == argv[i][0] && argv[i][1])
    // INCHI✔❌:         {
    // INCHI✔❌:             /* Parsing stand-alone executable/libinchi options */
    // INCHI✔❌:
    // INCHI✔❌:             pArg = argv[i] + 1;
    // INCHI✔❌:
    // INCHI✔❌:             bRecognizedOption = 2;
    // INCHI✔❌:             bVer1Options += 2;
    // INCHI✔❌:             /* always on: REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_STEREO */
    // INCHI✔❌:
    // INCHI✔❌:             got = set_common_options_by_parg(pArg, developer_options, ip, &bVer1DefaultMode, &nMode,
    // INCHI✔❌:                 &bINChIOutputOptions, &bINChIOutputOptions2,
    // INCHI✔❌:                 &bStdFormat, &bHashKey, &bHashXtra1, &bHashXtra2,
    // INCHI✔❌:                 &bFixSp3bug, &bFixFB2,
    // INCHI✔❌:                 &bAddPhosphineStereo, &bAddArsineStereo,
    // INCHI✔❌:                 &bNoStructLabels, &bPointedEdgeStereo,
    // INCHI✔❌:                 &bDoNotAddH, &bForcedChiralFlag, &bReconnectCoord,
    // INCHI✔❌:                 &bKetoEnolTaut,
    // INCHI✔❌:                 &b15TautNonRing,
    // INCHI✔❌:                 &bPT_06_00_Taut, &bPT_13_00_Taut, &bPT_16_00_Taut,
    // INCHI✔❌:                 &bPT_18_00_Taut, &bPT_22_00_Taut, &bPT_39_00_Taut,
    // INCHI✔❌:                 &bLooseTSACheck,
    // INCHI✔❌:                 &bLargeMolecules, &bPolymers,
    // INCHI✔❌:                 &bFoldPolymerSRU, &bFrameShiftScheme,
    // INCHI✔❌:                 &bStereoAtZz, &bNPZz,
    // INCHI✔❌:                 &bNoWarnings, &bMergeHash, &bHideInChI);
    // INCHI✔❌:
    // INCHI✔❌:             if (got)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ;
    // INCHI✔❌:             }
    // INCHI✔❌:             /* Input */
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "STDIO"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bNameSuffix = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if ( /* INPUT_NONE == ip->nInputType &&*/
    // INCHI✔❌:                 !inchi_memicmp(pArg, "SDF:", 4))
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* SDfile label */
    // INCHI✔❌:                 k = 0;
    // INCHI✔❌:                 mystrncpy(ip->szSdfDataHeader, pArg + 4, MAX_SDF_HEADER + 1);
    // INCHI✔❌:                 lrtrim(ip->szSdfDataHeader, &k);
    // INCHI✔❌:                 if (k)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     ip->pSdfLabel = ip->szSdfDataHeader;
    // INCHI✔❌:                     ip->pSdfValue = szSdfDataValue;
    // INCHI✔❌:                     if (INPUT_NONE == ip->nInputType)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         ip->nInputType = INPUT_SDFILE;
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                     ip->pSdfLabel = NULL;
    // INCHI✔❌:                     ip->pSdfValue = NULL;
    // INCHI✔❌:                     if (INPUT_NONE == ip->nInputType)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         ip->nInputType = INPUT_MOLFILE;
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "RSB:", 4) && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 mdbr = (int)strtol(pArg + 4, NULL, 10);
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             /* Output */
    // INCHI✔❌: #if ( !defined(TARGET_API_LIB) && !defined(TARGET_LIB_FOR_WINCHI) )
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "Tabbed") || !inchi_stricmp(pArg, "Tab"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bOutputStyle |= INCHI_OUT_TABBED_OUTPUT;
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "OUTPUTSDF"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* output SDfile */
    // INCHI✔❌:                 bOutputMolfileOnly = 1;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "SdfAtomsDT"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* output isotopes H as D and T in SDfile */
    // INCHI✔❌:                 bOutputMolfileDT = 1;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "D"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* display the structures */
    // INCHI✔❌:                 bDisplay |= 1;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "F", 1) && (c = (int)strtol(pArg + 1, (char**)&q, 10), q > pArg + 1))
    // INCHI✔❌:             {
    // INCHI✔❌:                 nFontSize = -c;                      /* struct. display font size */
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "EQU"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bCompareComponents = CMP_COMPONENTS;
    // INCHI✔❌:             }
    // INCHI✔❌: #if( OUTPUT_FILE_EXT == 1 )
    // INCHI✔❌:             else if (pArg[0] == '.' && numOutNameExt < (int)(sizeof(szOutNameExt) / sizeof(szOutNameExt[0])))
    // INCHI✔❌:             {
    // INCHI✔❌:                 strncpy(szOutNameExt[numOutNameExt], pArg, sizeof(szOutNameExt[0]) - 1);
    // INCHI✔❌:                 numOutNameExt++;
    // INCHI✔❌:                 if (ip->path[numOutNameExt])
    // INCHI✔❌:                 {
    // INCHI✔❌:                     ip->path[numOutNameExt] = ""; /*strcpy( ip->path[numOutNameExt], "");*/
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌:             /* djb-rwth: avoiding Error 98 for empty .mol files -- GH issue #25, thanks to @wijnand1 */
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "WarnOnEmptyStructure"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bAllowEmptyStructure = 1;
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             /* Generation options */
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "W", 1))
    // INCHI✔❌:             {
    // INCHI✔❌:                 long timeout_value;
    // INCHI✔❌:                 const char c1 = *(pArg + 1);
    // INCHI✔❌:                 if (c1 && (c1 == 'M' || c1 == 'm'))
    // INCHI✔❌:                 {
    // INCHI✔❌:                     /* "WMnumber", milliseconds */
    // INCHI✔❌:                     timeout_value = strtol(pArg + 2, (char**)&q, 10);
    // INCHI✔❌:                     if (timeout_value && q > pArg + 2 && *q == '\0')
    // INCHI✔❌:                     {
    // INCHI✔❌:                         if (errno == ERANGE || timeout_value < 0.0 || timeout_value>LONG_MAX) /* djb-rwth: addressing coverity ID #499550 -- the condition takes into account all possible overflows/errors */
    // INCHI✔❌:                         {
    // INCHI✔❌:                             timeout_value = 0;
    // INCHI✔❌:                             timeout_set_warning = 1;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         timeout_set_error = 0;
    // INCHI✔❌:                     }
    // INCHI✔❌:                     else
    // INCHI✔❌:                     {
    // INCHI✔❌:                         timeout_set_error = 1;
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                     /* expect "Wnumber", seconds */
    // INCHI✔❌:                     t = strtod(pArg + 1, (char**)&q);
    // INCHI✔❌:                     if (t && q > pArg + 1)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         if (*q != '\0')
    // INCHI✔❌:                         {
    // INCHI✔❌:                             timeout_set_warning = 1;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         if (errno == ERANGE || t < 0.0 || t * 1000.0 >(double)LONG_MAX)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             timeout_value = 0;
    // INCHI✔❌:                             timeout_set_warning = 1;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         else
    // INCHI✔❌:                         {
    // INCHI✔❌:                             timeout_value = (long)(t * 1000.0);  /* max. time per structure */
    // INCHI✔❌:                         }
    // INCHI✔❌:                         timeout_set_error = 0;
    // INCHI✔❌:                     }
    // INCHI✔❌:                     else
    // INCHI✔❌:                     {
    // INCHI✔❌:                         timeout_set_error = 1;
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:                 if (timeout_set_error == 0)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     ip->msec_MaxTime = timeout_value;
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             /*--- Conversion modes ---*/
    // INCHI✔❌: #if ( READ_INCHI_STRING == 1 )
    // INCHI✔❌:
    // INCHI✔❌: /*#if (BUILD_WITH_ENG_OPTIONS==1)*/
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "InChI2InChI"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* Read InChI Identifiers and output InChI Identifiers */
    // INCHI✔❌:                 ip->nInputType = INPUT_INCHI;
    // INCHI✔❌:                 ip->bReadInChIOptions |= READ_INCHI_OUTPUT_INCHI;
    // INCHI✔❌:                 ip->bReadInChIOptions &= ~READ_INCHI_TO_STRUCTURE;
    // INCHI✔❌:             }
    // INCHI✔❌:             /*#endif*/
    // INCHI✔❌:
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "InChI2Struct"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* Split InChI Identifiers into components */
    // INCHI✔❌:                 ip->bReadInChIOptions |= READ_INCHI_TO_STRUCTURE;
    // INCHI✔❌:                 ip->bReadInChIOptions &= ~READ_INCHI_OUTPUT_INCHI;
    // INCHI✔❌:                 ip->nInputType = INPUT_INCHI;
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "KeepBalanceP") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* When spliting InChI Identifiers into components: */
    // INCHI✔❌:                 /* If MobileH output then add p to each component;  */
    // INCHI✔❌:                 /* Otherwise add one more component containing balance */
    // INCHI✔❌:                 /* of protons and exchangeable isotopic H */
    // INCHI✔❌:                 ip->bReadInChIOptions |= READ_INCHI_KEEP_BALANCE_P;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌:             /*--- (engineering) Undo bug/draw fixes options ---*/
    // INCHI✔❌:
    // INCHI✔❌:             /* (developer_options) Old structure-perception and InChI creation options */
    // INCHI✔❌: #if ( FIX_ADJ_RAD == 1 )
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "FixRad") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bFixAdjacentRad = 1;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌: #if ( RENUMBER_ATOMS_AND_RECALC_V106 == 1 )
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "TestRenum") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bRenumber = 1;
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌: #if ( UNDERIVATIZE == 1 )
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "DoDRV") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bUnderivatize = 1;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌: #if( UNDERIVATIZE_REPORT == 1 )
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "DoDrvReport"))
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bUnderivatize = 3;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌: #endif
    // INCHI✔❌: #if ( RING2CHAIN == 1 )
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "DoR2C") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bRing2Chain = 1;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌: #if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "DoneOnly") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bIgnoreUnchanged = 1;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "MOVEPOS:", 8) && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bMovePositiveCharges = (0 != strtol(pArg + 8, NULL, 10));
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "NoADP") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bTgFlagHardAddRenProtons = 0;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             /* Tautomer perception off */
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "EXACT") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bVer1DefaultMode |= REQ_MODE_BASIC;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "ONLYRECSALT") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* do not disconnect salts */
    // INCHI✔❌:                 bDisconnectSalts = 0;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if ((!inchi_stricmp(pArg, "ONLYEXACT") || !inchi_stricmp(pArg, "ONLYFIXEDH")) && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bVer1DefaultMode |= REQ_MODE_BASIC;
    // INCHI✔❌:                 bVer1DefaultMode &= ~REQ_MODE_TAUT;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "ONLYNONISO") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bVer1DefaultMode |= REQ_MODE_NON_ISO;
    // INCHI✔❌:                 bVer1DefaultMode &= ~REQ_MODE_ISO;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "TAUT") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bVer1DefaultMode &= ~REQ_MODE_BASIC;
    // INCHI✔❌:                 bVer1DefaultMode |= REQ_MODE_TAUT;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "ONLYRECMET") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* do not disconnect metals */
    // INCHI✔❌:                 bDisconnectCoord = 0;
    // INCHI✔❌:                 bStdFormat = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             /*--- (hidden) Old output and other options ---*/
    // INCHI✔❌:
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "SdfSplit") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* Split single Molfiles into disconnected components */
    // INCHI✔❌:                 bOutputMolfileSplit = 1;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "DCR") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bDisplayCompositeResults = 1;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if ((!inchi_stricmp(pArg, "AUXFULL") || !inchi_stricmp(pArg, "AUXMAX")) && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* full aux info */
    // INCHI✔❌:                 bINChIOutputOptions &= ~(INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO); /* include short aux info */
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "AUXMIN") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* minimal aux info */
    // INCHI✔❌:                 bINChIOutputOptions &= ~INCHI_OUT_NO_AUX_INFO; /* include short aux info */
    // INCHI✔❌:                 bINChIOutputOptions |= INCHI_OUT_SHORT_AUX_INFO;
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌: #if ( READ_INCHI_STRING == 1 )
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "DDSRC") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bDisplayIfRestoreWarnings = 1;  /* InChI->Structure debugging: Display Debug Structure Restore Components */
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "NoVarH") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bTgFlagVariableProtons = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "FULL") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bVer1DefaultMode = VER103_DEFAULT_MODE;
    // INCHI✔❌:                 nMode = 0;
    // INCHI✔❌:                 bReconnectCoord = 1;            /* full output */
    // INCHI✔❌:                 bINChIOutputOptions = ((EMBED_REC_METALS_INCHI == 1) ? INCHI_OUT_EMBED_REC : 0) | INCHI_OUT_SHORT_AUX_INFO;
    // INCHI✔❌:                 ip->bCtPredecessors = 0;
    // INCHI✔❌:                 ip->bAbcNumbers = 0;
    // INCHI✔❌:                 bOutputStyle |= INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "MIN") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bVer1DefaultMode = VER103_DEFAULT_MODE;
    // INCHI✔❌:                 nMode = 0;
    // INCHI✔❌:                 bReconnectCoord = 1;            /* minimal output */
    // INCHI✔❌:                 bINChIOutputOptions = ((EMBED_REC_METALS_INCHI == 1) ? INCHI_OUT_EMBED_REC : 0) | INCHI_OUT_NO_AUX_INFO;            /* minimal compressed output */
    // INCHI✔❌:                 ip->bCtPredecessors = 1;
    // INCHI✔❌:                 ip->bAbcNumbers = 1;
    // INCHI✔❌:                 bOutputStyle |= INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "COMPRESS") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bAbcNumbers = 1;
    // INCHI✔❌:                 ip->bCtPredecessors = 1;             /* compressed output */
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌: #if ( READ_INCHI_STRING == 1 )
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "InChI2InChI")) /*&& developer_options)*/
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* Read InChI Identifiers and output InChI Identifiers */
    // INCHI✔❌:                 ip->nInputType = INPUT_INCHI;
    // INCHI✔❌:                 ip->bReadInChIOptions |= READ_INCHI_OUTPUT_INCHI;
    // INCHI✔❌:                 ip->bReadInChIOptions &= ~READ_INCHI_TO_STRUCTURE;
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "SplitInChI") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* Split InChI Identifiers into components */
    // INCHI✔❌:                 ip->bReadInChIOptions |= READ_INCHI_SPLIT_OUTPUT;
    // INCHI✔❌:             }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "MOLFILENUMBER") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bGetMolfileNumber |= 1;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "OutputPLAIN") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bOutputStyle |= INCHI_OUT_PLAIN_TEXT;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "OutputANNPLAIN") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bOutputStyle |= INCHI_OUT_PLAIN_TEXT_COMMENTS;
    // INCHI✔❌:                 bOutputStyle |= INCHI_OUT_WINCHI_WINDOW; /* debug */
    // INCHI✔❌:             }
    // INCHI✔❌:             else if ((!inchi_stricmp(pArg, "ONLYEXACT") || !inchi_stricmp(pArg, "ONLYFIXEDH")) && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bVer1DefaultMode |= REQ_MODE_BASIC;
    // INCHI✔❌:                 bVer1DefaultMode &= ~REQ_MODE_TAUT;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "ONLYNONISO") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bVer1DefaultMode |= REQ_MODE_NON_ISO;
    // INCHI✔❌:                 bVer1DefaultMode &= ~REQ_MODE_ISO;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "TAUT") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bVer1DefaultMode &= ~REQ_MODE_BASIC;
    // INCHI✔❌:                 bVer1DefaultMode |= REQ_MODE_TAUT;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "ONLYRECMET") && developer_options)
    // INCHI✔❌:             {  /* do not disconnect metals */
    // INCHI✔❌:                 bDisconnectCoord = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "ONLYRECSALT") && developer_options)
    // INCHI✔❌:             {  /* do not disconnect salts */
    // INCHI✔❌:                 bDisconnectSalts = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "MOVEPOS:", 8) && developer_options)
    // INCHI✔❌:             {   /* added -- 2010-03-01 DT */
    // INCHI✔❌:                 bMovePositiveCharges = (0 != strtol(pArg + 8, NULL, 10));
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "RSB:", 4) && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 mdbr = (int)strtol(pArg + 4, NULL, 10);
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "EQU") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bCompareComponents = CMP_COMPONENTS;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(pArg, "EQUNONISO") && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bCompareComponents = CMP_COMPONENTS | CMP_COMPONENTS_NONISO;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_memicmp(pArg, "OP:", 3) && developer_options)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bOutputPath = 1;
    // INCHI✔❌:                 strncpy(szOutputPath, pArg + 3, sizeof(szOutputPath) - 1);
    // INCHI✔❌:             }
    // INCHI✔❌:             /* eof developer_options */
    // INCHI✔❌:
    // INCHI✔❌:             /* Display unrecognized option */
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 bRecognizedOption = 0;
    // INCHI✔❌:
    // INCHI✔❌: #ifndef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌:                 inchi_ios_eprint(log_file, "Unrecognized optionQ3: \"%s\".\n", pArg);
    // INCHI✔❌: #endif
    // INCHI✔❌:             }
    // INCHI✔❌:             bVer1Options |= bRecognizedOption;
    // INCHI✔❌:
    // INCHI✔❌:         } /* eof Parsing stand-alone executable/libinchi options */
    // INCHI✔❌:
    // INCHI✔❌:         else if (ip->num_paths < MAX_NUM_PATHS)
    // INCHI✔❌:         {
    // INCHI✔❌:             char* sz;
    // INCHI✔❌: #if( ALLOW_EMPTY_PATHS == 1 )
    // INCHI✔❌:             if (argv[i])
    // INCHI✔❌: #else
    // INCHI✔❌:             if (argv[i] && argv[i][0])
    // INCHI✔❌: #endif
    // INCHI✔❌:             {
    // INCHI✔❌:                 if ((sz = (char*)inchi_malloc((strlen(argv[i]) + 1) * sizeof(sz[0])))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:                 {
    // INCHI✔❌:                     strcpy(sz, argv[i]);
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 ip->path[ip->num_paths++] = sz;
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:     } /*  eof parsing argv loop */
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌:     /* Print messages and set controil variables according to just parsed options */
    // INCHI✔❌:
    // INCHI✔❌:     /* Timeout option(s) */
    // INCHI✔❌:     if (timeout_set_warning)
    // INCHI✔❌:     {
    // INCHI✔❌:         inchi_ios_eprint(log_file, "Warning: timeout value may have been modified (truncated?) due to number formatting issues;\n");
    // INCHI✔❌:     }
    // INCHI✔❌:     if (timeout_set_error)
    // INCHI✔❌:     {
    // INCHI✔❌:         inchi_ios_eprint(log_file, "Warning: specified timeout value was ignored due to invalid number format, using the default;\n");
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* InChIKey option(s) */
    // INCHI✔❌:     if (bHashKey != 0)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* Suppress InChIKey calculation if:                */
    // INCHI✔❌:         /* compressed output OR Inchi2Struct OR Inchi2Inchi */
    // INCHI✔❌:         if ((ip->bAbcNumbers == 1) && (ip->bCtPredecessors == 1))
    // INCHI✔❌:         {
    // INCHI✔❌:             bHashKey = 0;
    // INCHI✔❌: #ifndef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌:             inchi_ios_eprint(log_file, "Terminating: generation of InChIKey is not available with 'Compress' option\n");
    // INCHI✔❌:             return -1;
    // INCHI✔❌: #endif
    // INCHI✔❌:         }
    // INCHI✔❌:         if (ip->nInputType == INPUT_INCHI)
    // INCHI✔❌:         {
    // INCHI✔❌:             bHashKey = 0;
    // INCHI✔❌: #ifndef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌:             inchi_ios_eprint(log_file, "Terminating: generation of InChIKey is not available in InChI conversion mode\n");
    // INCHI✔❌:             return -1;
    // INCHI✔❌: #endif
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:         {
    // INCHI✔❌:             if (bOutputMolfileOnly == 1)
    // INCHI✔❌:             {
    // INCHI✔❌:                 bHashKey = 0;
    // INCHI✔❌: #ifndef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌:                 inchi_ios_eprint(log_file, "Terminating: generation of InChIKey is not available with 'OutputSDF' option\n");
    // INCHI✔❌:                 return -1;
    // INCHI✔❌: #endif
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (bNameSuffix || bOutputPath)
    // INCHI✔❌:     {
    // INCHI✔❌:         const char szNUL[] = "NUL";
    // INCHI✔❌:         /* fix for AMD processor: use const char[] instead of just "NUL" constant 2008-11-5 DT */
    // INCHI✔❌:         const char* p = NULL;
    // INCHI✔❌:         char* p_prev = NULL; /* djb-rwth: avoiding use of memory after it is freed */
    // INCHI✔❌:         char* r = NULL;
    // INCHI✔❌:         char* sz;
    // INCHI✔❌:         int len;
    // INCHI✔❌:         /*  find the 1st path */
    // INCHI✔❌:         for (i = 0; i < MAX_NUM_PATHS; i++)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (!p && ip->path[i] && ip->path[i][0])
    // INCHI✔❌:             {
    // INCHI✔❌:                 p = ip->path[i];
    // INCHI✔❌:                 break;
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:         /* fix output path */
    // INCHI✔❌:         if (bOutputPath && szOutputPath[0] && p)
    // INCHI✔❌:         {
    // INCHI✔❌:             /* remove last slash */
    // INCHI✔❌:             len = (int)strlen(szOutputPath);
    // INCHI✔❌:             if (len > 0 && szOutputPath[len - 1] != INCHI_PATH_DELIM)
    // INCHI✔❌:             {
    // INCHI✔❌:                 szOutputPath[len++] = INCHI_PATH_DELIM;
    // INCHI✔❌:                 szOutputPath[len] = '\0';
    // INCHI✔❌:             }
    // INCHI✔❌:             if (len > 0 && (r = (char*)strrchr(p, INCHI_PATH_DELIM)) && r[1])
    // INCHI✔❌:             {
    // INCHI✔❌:                 strcat(szOutputPath, r + 1);
    // INCHI✔❌:                 p = szOutputPath;
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌: /* djb-rwth: copying the value of p */
    // INCHI✔❌: #ifdef _WIN32
    // INCHI✔❌:         p_prev = _strdup(p);
    // INCHI✔❌: #else
    // INCHI✔❌:         p_prev = inchi__strdup(p);
    // INCHI✔❌: #endif
    // INCHI✔❌:         /*  add missing paths */
    // INCHI✔❌:         /* djb-rwth: this whole block had to be rewritten to avoid use of memory after it is freed */
    // INCHI✔❌:         for (i = 0; p_prev && i < MAX_NUM_PATHS; i++)
    // INCHI✔❌:         {
    // INCHI✔❌:             /* fix for AMD processor: changed order 2008-11-5 DT */
    // INCHI✔❌:             if (!ip->path[i] || !ip->path[i][0])
    // INCHI✔❌:             {
    // INCHI✔❌: #if ( BUILD_WITH_AMI == 1 ) && ( OUTPUT_FILE_EXT == 1 )
    // INCHI✔❌:                 char* pLastExt = (i && numOutNameExt >= i) ? strrchr((char*)p_prev, '.') : 0;
    // INCHI✔❌:                 char* pLastSlash = (i && numOutNameExt >= i) ? strrchr((char*)p_prev, INCHI_PATH_DELIM) : 0;
    // INCHI✔❌:                 if (pLastExt && pLastSlash && pLastSlash > pLastExt)
    // INCHI✔❌:                     pLastExt = NULL;
    // INCHI✔❌: #else
    // INCHI✔❌:                 char* pLastExt = NULL;
    // INCHI✔❌: #endif
    // INCHI✔❌:                 len = (int)strlen(p_prev) + strlen(szNameSuffix) + strlen(ext[i]);
    // INCHI✔❌:                 if ((sz = (char*)inchi_malloc(((long long)len + 1) * sizeof(sz[0])))) /* djb-rwth: cast operators added; addressing and ignoring LLVM warnings */
    // INCHI✔❌:                 {
    // INCHI✔❌:                     strcpy(sz, p_prev); /* djb-rwth: fix for use of memory after being freed */
    // INCHI✔❌: #if ( BUILD_WITH_AMI == 1 ) && ( OUTPUT_FILE_EXT == 1 )
    // INCHI✔❌:                     if (pLastExt)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         strcpy(sz + (pLastExt - p), szOutNameExt[i - 1]);
    // INCHI✔❌:                     }
    // INCHI✔❌: #endif
    // INCHI✔❌:                     strcat(sz, szNameSuffix);
    // INCHI✔❌:                     if (!pLastExt)
    // INCHI✔❌:                         strcat(sz, ext[i]);
    // INCHI✔❌:                     ip->num_paths++;
    // INCHI✔❌:                     if (ip->path[i])
    // INCHI✔❌:                     {
    // INCHI✔❌:                         inchi_free((char*)ip->path[i]); /* eliminate memory leak 2013-12-18 DCh */
    // INCHI✔❌:                     }
    // INCHI✔❌:                     ip->path[i] = sz;
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!inchi_stricmp(ip->path[i], szNUL))
    // INCHI✔❌:             {
    // INCHI✔❌:                 inchi_free((char*)ip->path[i]); /* cast deliberately const qualifier */
    // INCHI✔❌:                 ip->path[i] = NULL;
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:         free(p_prev); /* djb-rwth: freeing memory reserved for auxiliary variable */
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* inchi2inchi and inchi2struct option(s) */
    // INCHI✔❌: #if ( READ_INCHI_STRING == 1 )
    // INCHI✔❌:     if (INPUT_INCHI == ip->nInputType)
    // INCHI✔❌:     {
    // INCHI✔❌:         bCompareComponents = 0;
    // INCHI✔❌:         /*bDisplayCompositeResults = 0;*/
    // INCHI✔❌: #if ( I2S_MODIFY_OUTPUT == 1 )
    // INCHI✔❌:         if (!(ip->bReadInChIOptions & READ_INCHI_TO_STRUCTURE))
    // INCHI✔❌: #endif
    // INCHI✔❌:         {
    // INCHI✔❌:             bOutputMolfileOnly = 0;
    // INCHI✔❌:             bINChIOutputOptions |= INCHI_OUT_NO_AUX_INFO;
    // INCHI✔❌:             bINChIOutputOptions &= ~INCHI_OUT_SHORT_AUX_INFO;
    // INCHI✔❌:             bINChIOutputOptions &= ~INCHI_OUT_ONLY_AUX_INFO;
    // INCHI✔❌:             /* bNoStructLabels   = 1; */
    // INCHI✔❌:         }
    // INCHI✔❌:         ip->bDisplayIfRestoreWarnings = bDisplayIfRestoreWarnings;
    // INCHI✔❌:         if (!(bINChIOutputOptions &
    // INCHI✔❌:             (INCHI_OUT_SDFILE_ONLY |       /* not in bINChIOutputOptions yet */
    // INCHI✔❌:                 INCHI_OUT_PLAIN_TEXT |        /* not in bINChIOutputOptions yet */
    // INCHI✔❌:                 INCHI_OUT_PLAIN_TEXT_COMMENTS /* not in bINChIOutputOptions yet */
    // INCHI✔❌:                 )
    // INCHI✔❌:             )
    // INCHI✔❌: #if ( I2S_MODIFY_OUTPUT == 1 )
    // INCHI✔❌:             && !bOutputMolfileOnly
    // INCHI✔❌:             && !(bOutputStyle & (INCHI_OUT_PLAIN_TEXT))
    // INCHI✔❌: #endif
    // INCHI✔❌:             )
    // INCHI✔❌:         {
    // INCHI✔❌:             bINChIOutputOptions |= INCHI_OUT_PLAIN_TEXT;
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌:     if (bVer1Options)
    // INCHI✔❌:     {
    // INCHI✔❌:         nMode |= bVer1DefaultMode;
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (bReleaseVersion)
    // INCHI✔❌:     {
    // INCHI✔❌:         nMode |= nReleaseMode;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌: #if ( defined(COMPILE_ANSI_ONLY) || defined(TARGET_LIB_FOR_WINCHI) )
    // INCHI✔❌:     if (bCompareComponents && !(bDisplay & 1))
    // INCHI✔❌:     {
    // INCHI✔❌:         bCompareComponents = 0;
    // INCHI✔❌:     }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌:     /* Save original options */
    // INCHI✔❌:     /* nOrigMode = nMode; */
    // INCHI✔❌: #ifndef COMPILE_ANSI_ONLY
    // INCHI✔❌:     ip->dp.sdp.nFontSize = nFontSize;
    // INCHI✔❌:     ip->dp.sdp.ulDisplTime = *ulDisplTime;
    // INCHI✔❌:     ip->bDisplay = bDisplay;
    // INCHI✔❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌:     ip->bDisplayCompositeResults = bDisplay;
    // INCHI✔❌: #else
    // INCHI✔❌:     ip->bDisplayCompositeResults = bDisplayCompositeResults;
    // INCHI✔❌: #endif
    // INCHI✔❌: #else
    // INCHI✔❌:     ip->bDisplayEachComponentINChI = 0;
    // INCHI✔❌:     bCompareComponents = 0;
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌:     ip->bMergeAllInputStructures = bMergeAllInputStructures;
    // INCHI✔❌:     ip->bDoNotAddH = bDoNotAddH;
    // INCHI✔❌:
    // INCHI✔❌:     /*  Set default options */
    // INCHI✔❌:     if (!nMode || nMode == REQ_MODE_STEREO)
    // INCHI✔❌:     {
    // INCHI✔❌:         /*  requested all output */
    // INCHI✔❌:         nMode |= (REQ_MODE_BASIC | REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_NON_ISO | REQ_MODE_STEREO);
    // INCHI✔❌:     }
    // INCHI✔❌:     else
    // INCHI✔❌:     {
    // INCHI✔❌:         if (!(nMode & (REQ_MODE_BASIC | REQ_MODE_TAUT)))
    // INCHI✔❌:         {
    // INCHI✔❌:             nMode |= (REQ_MODE_BASIC | REQ_MODE_TAUT);
    // INCHI✔❌:         }
    // INCHI✔❌:         if ((nMode & REQ_MODE_STEREO) && !(nMode & (REQ_MODE_ISO | REQ_MODE_NON_ISO)))
    // INCHI✔❌:         {
    // INCHI✔❌:             nMode |= (REQ_MODE_ISO | REQ_MODE_NON_ISO);
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:     /*  if the user requested isotopic then unconditionally add non-isotopic output. */
    // INCHI✔❌:     if (nMode & REQ_MODE_ISO)
    // INCHI✔❌:     {
    // INCHI✔❌:         nMode |= REQ_MODE_NON_ISO;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌: #if ( MIN_SB_RING_SIZE > 0 )
    // INCHI✔❌:     if (mdbr)
    // INCHI✔❌:     {
    // INCHI✔❌:         nMinDbRinSize = mdbr;
    // INCHI✔❌:     }
    // INCHI✔❌:     nMode |= (nMinDbRinSize << REQ_MODE_MIN_SB_RING_SHFT) & REQ_MODE_MIN_SB_RING_MASK;
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌:     /* Input file */
    // INCHI✔❌:     if (ip->nInputType == INPUT_NONE && ip->num_paths > 0)
    // INCHI✔❌:     {
    // INCHI✔❌:         ip->nInputType = INPUT_MOLFILE; /*  default */
    // INCHI✔❌:     }
    // INCHI✔❌:     ip->nMode = nMode;
    // INCHI✔❌:     /* Compare components */
    // INCHI✔❌:     if ((bCompareComponents & CMP_COMPONENTS) && (nMode & REQ_MODE_BASIC))
    // INCHI✔❌:     {
    // INCHI✔❌:         bCompareComponents |= CMP_COMPONENTS_NONTAUT; /* compare non-tautomeric */
    // INCHI✔❌:     }
    // INCHI✔❌:     ip->bCompareComponents = bCompareComponents;
    // INCHI✔❌:     /* Output */
    // INCHI✔❌:     ip->bINChIOutputOptions = bINChIOutputOptions | (bOutputMolfileOnly ? INCHI_OUT_SDFILE_ONLY : 0);
    // INCHI✔❌:     if (bOutputMolfileOnly)
    // INCHI✔❌:     {
    // INCHI✔❌:         bOutputStyle &= ~(INCHI_OUT_PLAIN_TEXT |
    // INCHI✔❌:             INCHI_OUT_PLAIN_TEXT_COMMENTS |
    // INCHI✔❌:             INCHI_OUT_TABBED_OUTPUT);
    // INCHI✔❌: #if ( SDF_OUTPUT_DT == 1 )
    // INCHI✔❌:         ip->bINChIOutputOptions |= bOutputMolfileDT ? INCHI_OUT_SDFILE_ATOMS_DT : 0;
    // INCHI✔❌:         ip->bINChIOutputOptions |= bOutputMolfileSplit ? INCHI_OUT_SDFILE_SPLIT : 0;
    // INCHI✔❌: #endif
    // INCHI✔❌:     }
    // INCHI✔❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌:     if (!(bDisplay & 1))
    // INCHI✔❌:     {
    // INCHI✔❌:         bOutputStyle &= ~(INCHI_OUT_PLAIN_TEXT_COMMENTS); /* do not ouput comments in wINChI text file results */
    // INCHI✔❌:     }
    // INCHI✔❌:     else
    // INCHI✔❌:     {
    // INCHI✔❌:         bOutputStyle |= INCHI_OUT_WINCHI_WINDOW;
    // INCHI✔❌:     }
    // INCHI✔❌: #endif
    // INCHI✔❌:     ip->bINChIOutputOptions |= bOutputStyle;
    // INCHI✔❌:     ip->bNoStructLabels = bNoStructLabels;
    // INCHI✔❌:
    // INCHI✔❌:     /* Processing options */
    // INCHI✔❌:     if (bForcedChiralFlag)
    // INCHI✔❌:     {
    // INCHI✔❌:         ip->bChiralFlag = bForcedChiralFlag;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* Tautomeric/salts options */
    // INCHI✔❌:     ip->bTautFlags = 0;
    // INCHI✔❌:     ip->bTautFlagsDone = 0;
    // INCHI✔❌:
    // INCHI✔❌:     /* Find regular tautomerism */
    // INCHI✔❌:     ip->bTautFlags |= TG_FLAG_TEST_TAUT__ATOMS;
    // INCHI✔❌:
    // INCHI✔❌:     /* Disconnect salts */
    // INCHI✔❌:     ip->bTautFlags |= bDisconnectSalts ? TG_FLAG_DISCONNECT_SALTS : 0;
    // INCHI✔❌:
    // INCHI✔❌:     /* If possible, find long-range H/(-) taut. on =C-OH, >C=O    */
    // INCHI✔❌:     ip->bTautFlags |= bAcidTautomerism ? TG_FLAG_TEST_TAUT__SALTS : 0;
    // INCHI✔❌:
    // INCHI✔❌:     /* Allow long-range movement of N(+), P(+) charges           */
    // INCHI✔❌:     ip->bTautFlags |= bMovePositiveCharges ? TG_FLAG_MOVE_POS_CHARGES : 0;
    // INCHI✔❌:
    // INCHI✔❌:     /* Multi-attachement long-range H/(-) taut. on =C-OH, >C=O   */
    // INCHI✔❌:     ip->bTautFlags |= (bAcidTautomerism > 1) ? TG_FLAG_TEST_TAUT2_SALTS : 0;
    // INCHI✔❌:
    // INCHI✔❌:     /* (Debug) allow to find long-range H-only tautomerism on =C-OH, >C=O */
    // INCHI✔❌:     ip->bTautFlags |= (bUnchargedAcidTaut == 1) ? TG_FLAG_ALLOW_NO_NEGTV_O : 0;
    // INCHI✔❌:
    // INCHI✔❌:     /* Merge =C-OH and >C=O containing t-groups and other =C-OH groups */
    // INCHI✔❌:     ip->bTautFlags |= bMergeSaltTGroups ? TG_FLAG_MERGE_TAUT_SALTS : 0;
    // INCHI✔❌:     ip->bTautFlags |= bDisconnectCoord ? TG_FLAG_DISCONNECT_COORD : 0;
    // INCHI✔❌:     ip->bTautFlags |= (bDisconnectCoord && bReconnectCoord) ? TG_FLAG_RECONNECT_COORD : 0;
    // INCHI✔❌:     ip->bTautFlags |= bDisconnectCoordChkVal ? TG_FLAG_CHECK_VALENCE_COORD : 0;
    // INCHI✔❌:     ip->bTautFlags |= bTgFlagVariableProtons ? TG_FLAG_VARIABLE_PROTONS : 0;
    // INCHI✔❌:     ip->bTautFlags |= bTgFlagHardAddRenProtons ? TG_FLAG_HARD_ADD_REM_PROTONS : 0;
    // INCHI✔❌:     ip->bTautFlags |= bKetoEnolTaut ? TG_FLAG_KETO_ENOL_TAUT : 0;
    // INCHI✔❌:     ip->bTautFlags |= b15TautNonRing ? TG_FLAG_1_5_TAUT : 0;
    // INCHI✔❌:
    // INCHI✔❌:     /*^^^ IPl 2020-04-02 added to forcefully enable new tauto rules (in test_ixa, etc.) */
    // INCHI✔❌:     /*^^^ IPl 2020-10-26 removed set to 1
    // INCHI✔❌:     bPT_22_00_Taut = 1;
    // INCHI✔❌:     bPT_16_00_Taut = 1;
    // INCHI✔❌:     bPT_06_00_Taut = 1;
    // INCHI✔❌:     bPT_39_00_Taut = 1;
    // INCHI✔❌:     */
    // INCHI✔❌:
    // INCHI✔❌:     /*^^^ IPl 2020-04-02 */
    // INCHI✔❌:     ip->bTautFlags |= bPT_22_00_Taut ? TG_FLAG_PT_22_00 : 0;
    // INCHI✔❌:     ip->bTautFlags |= bPT_16_00_Taut ? TG_FLAG_PT_16_00 : 0;
    // INCHI✔❌:     ip->bTautFlags |= bPT_06_00_Taut ? TG_FLAG_PT_06_00 : 0;
    // INCHI✔❌:     ip->bTautFlags |= bPT_39_00_Taut ? TG_FLAG_PT_39_00 : 0;
    // INCHI✔❌:     ip->bTautFlags |= bPT_13_00_Taut ? TG_FLAG_PT_13_00 : 0;
    // INCHI✔❌:     ip->bTautFlags |= bPT_18_00_Taut ? TG_FLAG_PT_18_00 : 0;
    // INCHI✔❌:
    // INCHI✔❌: #ifdef STEREO_WEDGE_ONLY
    // INCHI✔❌:     ip->bTautFlags |= bPointedEdgeStereo ? TG_FLAG_POINTED_EDGE_STEREO : 0;
    // INCHI✔❌: #endif
    // INCHI✔❌: #if ( FIX_ADJ_RAD == 1 )
    // INCHI✔❌:     ip->bTautFlags |= bFixAdjacentRad ? TG_FLAG_FIX_ADJ_RADICALS : 0;
    // INCHI✔❌: #endif
    // INCHI✔❌:     ip->bTautFlags |= bAddPhosphineStereo ? TG_FLAG_PHOSPHINE_STEREO : 0;
    // INCHI✔❌:     ip->bTautFlags |= bAddArsineStereo ? TG_FLAG_ARSINE_STEREO : 0;
    // INCHI✔❌:     ip->bTautFlags |= bFixSp3bug ? TG_FLAG_FIX_SP3_BUG : 0;
    // INCHI✔❌:
    // INCHI✔❌:     /* Bug fixes */
    // INCHI✔❌:     if (bFixFB2)
    // INCHI✔❌:     {
    // INCHI✔❌: #if ( FIX_ISO_FIXEDH_BUG == 1 )
    // INCHI✔❌:         ip->bTautFlags |= TG_FLAG_FIX_ISO_FIXEDH_BUG; /* accomodate FIX_ISO_FIXEDH_BUG */
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌: #if ( FIX_TERM_H_CHRG_BUG == 1 )
    // INCHI✔❌:         ip->bTautFlags |= TG_FLAG_FIX_TERM_H_CHRG_BUG; /* accomodate FIX_TERM_H_CHRG_BUG */
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌: #if ( FIX_TRANSPOSITION_CHARGE_BUG == 1 )
    // INCHI✔❌:         ip->bINChIOutputOptions |= INCHI_OUT_FIX_TRANSPOSITION_CHARGE_BUG;
    // INCHI✔❌: #endif
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (!ip->nInputType)
    // INCHI✔❌:     {
    // INCHI✔❌:         ip->nInputType = INPUT_MOLFILE;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* Check if /SNon requested turn OFF SUU/SLUUD */
    // INCHI✔❌:     if (!(ip->nMode & REQ_MODE_STEREO))
    // INCHI✔❌:     {
    // INCHI✔❌:         ip->nMode &= ~REQ_MODE_DIFF_UU_STEREO;
    // INCHI✔❌:         ip->nMode &= ~(REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU);
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* Standard InChI ? */
    // INCHI✔❌:     if (bStdFormat)
    // INCHI✔❌:     {
    // INCHI✔❌:         ip->bINChIOutputOptions |= INCHI_OUT_STDINCHI;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* InChIKey ? */
    // INCHI✔❌:     if (!bHashKey)
    // INCHI✔❌:     {
    // INCHI✔❌:         ip->bCalcInChIHash = INCHIHASH_NONE;
    // INCHI✔❌:     }
    // INCHI✔❌:     else
    // INCHI✔❌:     {
    // INCHI✔❌:         ip->bCalcInChIHash = INCHIHASH_KEY;
    // INCHI✔❌:     }
    // INCHI✔❌:     /* Extension(s) to hash (in non-std mode only) ? */
    // INCHI✔❌:     if (!bHashKey)
    // INCHI✔❌:     {
    // INCHI✔❌:         if ((bHashXtra1 != 0) || (bHashXtra2 != 0))
    // INCHI✔❌:         {
    // INCHI✔❌:             inchi_ios_eprint(log_file, "Hash extension(s) not generated: InChIKey not requested");
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:     else
    // INCHI✔❌:     {
    // INCHI✔❌:         if (bHashXtra1)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (bHashXtra2)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bCalcInChIHash = INCHIHASH_KEY_XTRA1_XTRA2;
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bCalcInChIHash = INCHIHASH_KEY_XTRA1;
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:         else if (bHashXtra2)
    // INCHI✔❌:         {
    // INCHI✔❌:             ip->bCalcInChIHash = INCHIHASH_KEY_XTRA2;
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     ip->bLargeMolecules = bLargeMolecules;
    // INCHI✔❌:     ip->bLooseTSACheck = bLooseTSACheck;
    // INCHI✔❌:
    // INCHI✔❌:     ip->bNPZz = bNPZz;
    // INCHI✔❌:     ip->bStereoAtZz = bStereoAtZz;  /*STEREO_AT_ZZ;*/
    // INCHI✔❌:
    // INCHI✔❌:     ip->bNoWarnings = bNoWarnings;
    // INCHI✔❌:     ip->bMergeHash = bMergeHash;
    // INCHI✔❌:     ip->bHideInChI = bHideInChI;
    // INCHI✔❌:
    // INCHI✔❌:     ip->bPolymers = bPolymers;
    // INCHI✔❌:     ip->bFoldPolymerSRU = bFoldPolymerSRU;
    // INCHI✔❌:     ip->bFrameShiftScheme = bFrameShiftScheme;
    // INCHI✔❌:
    // INCHI✔❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔❌:     ip->bLargeMolecules = 1;
    // INCHI✔❌:     /*ip->bPolymers = POLYMERS_MODERN; */
    // INCHI✔❌:     /*ip->bNPZz = 1;*/
    // INCHI✔❌:     /*ip->bFrameShiftScheme = FSS_NONE;*/
    // INCHI✔❌:     /*ip->bFoldPolymerSRU = 1;*/
    // INCHI✔❌:
    // INCHI✔❌: #if ( UNDERIVATIZE == 1 )
    // INCHI✔❌:     ip->bUnderivatize = 1;
    // INCHI✔❌:     if (ip->bUnderivatize)
    // INCHI✔❌:     {
    // INCHI✔❌:         ip->bINChIOutputOptions &= ~INCHI_OUT_STDINCHI;
    // INCHI✔❌:     }
    // INCHI✔❌: #endif
    // INCHI✔❌:     ip->bNoWarnings = 0;
    // INCHI✔❌:     ip->bMergeHash = 0;
    // INCHI✔❌:     ip->bHideInChI = 0;
    // INCHI✔❌:
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌:     ip->bINChIOutputOptions2 = bINChIOutputOptions2;
    // INCHI✔❌:
    // INCHI✔❌:     return 0;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: ReadCommandLineParms

    let mut options = CommonOptionsByParg {
        ver1_default_mode: (REQ_MODE_TAUT
            | REQ_MODE_ISO
            | REQ_MODE_STEREO
            | REQ_MODE_SB_IGN_ALL_UU
            | REQ_MODE_SC_IGN_ALL_UU) as INCHI_MODE,
        mode: 0,
        inchi_output_options: INCHI_OUT_EMBED_REC as i32,
        inchi_output_options2: 0,
        std_format: 1,
        fix_sp3_bug: 1,
        fix_fb2: 1,
        add_phosphine_stereo: 1,
        add_arsine_stereo: 1,
        pointed_edge_stereo: 1,
        polymers: POLYMERS_NO as i32,
        frame_shift_scheme: tagFrameShifScheme_FSS_STARS_CYCLED as i32,
        ..CommonOptionsByParg::default()
    };
    options.ver1_default_mode &= !(REQ_MODE_BASIC as INCHI_MODE);
    let mut b_ver1_options = 1_i32;
    let mut b_name_suffix = 1_i32;
    let b_merge_all_input_structures = 0_i32;
    let mut b_output_style = INCHI_OUT_PLAIN_TEXT as i32;
    let mut b_output_molfile_only = 0_i32;
    let mut b_output_molfile_dt = 0_i32;
    let b_output_molfile_split = 0_i32;
    let mut timeout_set_warning = 0_i32;
    let mut timeout_set_error = 0_i32;

    *ip = INPUT_PARMS::default();
    ip.bFixNonUniformDraw = 1;
    ip.msec_MaxTime = 0;
    *ulDisplTime = 0;

    let argc = if argc <= 1 {
        1_usize
    } else {
        usize::try_from(argc).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
    };
    if argc > argv.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    for argument_index in 1..argc {
        let argument = argv[argument_index];
        let argument_bytes = heap.slice(argument)?;
        let argument_nul = argument_bytes
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        if argument_nul >= 2
            && argument_bytes[0] as u8 == INCHI_OPTION_PREFX
            && argument_bytes[1] != 0
        {
            let p_arg = argument.offset(1)?;
            b_ver1_options = b_ver1_options
                .checked_add(2)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            let got = set_common_options_by_parg(heap, p_arg, 0, ip, &mut options)?;
            if got != 0 {
                continue;
            }
            if option_literal(heap, p_arg, "STDIO")? {
                b_name_suffix = 0;
            } else if option_prefix(heap, p_arg, "SDF:")? {
                let header = heap.allocate_model_storage(vec![0_i8; 65])?;
                mystrncpy(heap, header, p_arg.offset(4)?, 65)?;
                let mut length = 0_i32;
                lrtrim(heap, header, Some(&mut length))?;
                ip.szSdfDataHeader
                    .copy_from_slice(heap.slice(header.as_const())?);
                if length != 0 {
                    ip.pSdfLabel = header;
                    ip.pSdfValue = szSdfDataValue;
                    if ip.nInputType == tagInputType_INPUT_NONE {
                        ip.nInputType = tagInputType_INPUT_SDFILE;
                    }
                } else {
                    inchi_free(heap, header)?;
                    ip.pSdfLabel = SourceMutPointer::null();
                    ip.pSdfValue = SourceMutPointer::null();
                    if ip.nInputType == tagInputType_INPUT_NONE {
                        ip.nInputType = tagInputType_INPUT_MOLFILE;
                    }
                }
            } else if option_literal(heap, p_arg, "OUTPUTSDF")? {
                b_output_molfile_only = 1;
            } else if option_literal(heap, p_arg, "SdfAtomsDT")? {
                b_output_molfile_dt = 1;
            } else if option_literal(heap, p_arg, "D")? {
                // Display state is discarded by the active COMPILE_ANSI_ONLY branch.
            } else if option_prefix(heap, p_arg, "F")? {
                let (_, end) = source_strtol_long_with_end(heap, p_arg.offset(1)?)?;
                if end == 0 {
                    let bytes = source_c_string(heap, p_arg)?;
                    eprint_call(
                        heap,
                        log_file.as_deref_mut(),
                        "Unrecognized optionQ3: \"%s\".\n",
                        vec![PrintArgument::Bytes(bytes)],
                    )?;
                }
            } else if option_literal(heap, p_arg, "EQU")? {
                // bCompareComponents is unconditionally cleared by COMPILE_ANSI_ONLY.
            } else if option_literal(heap, p_arg, "WarnOnEmptyStructure")? {
                ip.bAllowEmptyStructure = 1;
            } else if option_prefix(heap, p_arg, "W")? {
                let bytes = heap.slice(p_arg)?;
                let c1 = *bytes.get(1).ok_or(SourceHeapError::PointerOutOfBounds)? as u8;
                let timeout_value;
                if matches!(c1, b'M' | b'm') {
                    let (parsed, end) = source_strtol_long_with_end(heap, p_arg.offset(2)?)?;
                    let tail = heap.slice(p_arg.offset(2)?)?;
                    if parsed != 0 && end > 0 && tail.get(end) == Some(&0) {
                        timeout_value = if heap.source_errno() == 34 || parsed < 0 {
                            timeout_set_warning = 1;
                            0
                        } else {
                            parsed
                        };
                        timeout_set_error = 0;
                    } else {
                        timeout_value = 0;
                        timeout_set_error = 1;
                    }
                } else {
                    let (parsed, end) = source_strtod_with_end(heap, p_arg.offset(1)?)?;
                    let tail = heap.slice(p_arg.offset(1)?)?;
                    if parsed != 0.0 && end > 0 {
                        if tail.get(end) != Some(&0) {
                            timeout_set_warning = 1;
                        }
                        timeout_value = if heap.source_errno() == 34
                            || parsed < 0.0
                            || parsed * 1000.0 > i64::MAX as f64
                        {
                            timeout_set_warning = 1;
                            0
                        } else if parsed.is_nan() {
                            return Err(SourceHeapError::UnsupportedSourceBehavior);
                        } else {
                            (parsed * 1000.0) as i64
                        };
                        timeout_set_error = 0;
                    } else {
                        timeout_value = 0;
                        timeout_set_error = 1;
                    }
                }
                if timeout_set_error == 0 {
                    ip.msec_MaxTime = timeout_value;
                }
            } else if option_literal(heap, p_arg, "InChI2InChI")? {
                ip.nInputType = tagInputType_INPUT_INCHI;
                ip.bReadInChIOptions |= READ_INCHI_OUTPUT_INCHI as i32;
                ip.bReadInChIOptions &= !(READ_INCHI_TO_STRUCTURE as i32);
            } else if option_literal(heap, p_arg, "InChI2Struct")? {
                ip.bReadInChIOptions |= READ_INCHI_TO_STRUCTURE as i32;
                ip.bReadInChIOptions &= !(READ_INCHI_OUTPUT_INCHI as i32);
                ip.nInputType = tagInputType_INPUT_INCHI;
            } else if option_literal(heap, p_arg, "DoDrvReport")? {
                ip.bUnderivatize = 3;
                options.std_format = 0;
            } else {
                let bytes = source_c_string(heap, p_arg)?;
                eprint_call(
                    heap,
                    log_file.as_deref_mut(),
                    "Unrecognized optionQ3: \"%s\".\n",
                    vec![PrintArgument::Bytes(bytes)],
                )?;
            }
        } else if ip.num_paths < MAX_NUM_PATHS as i32 && argument_nul != 0 {
            let copied = argument_bytes[..=argument_nul].to_vec();
            let allocation = match inchi_malloc(
                heap,
                u64::try_from(copied.len())
                    .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?,
            ) {
                Ok(allocation) => {
                    heap.slice_mut(allocation)?.copy_from_slice(&copied);
                    allocation.as_const()
                }
                Err(SourceHeapError::AllocationFailed) => SourceConstPointer::null(),
                Err(error) => return Err(error),
            };
            let path_index = usize::try_from(ip.num_paths)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            ip.path[path_index] = allocation;
            ip.num_paths += 1;
        }
    }

    if timeout_set_warning != 0 {
        eprint_call(
            heap,
            log_file.as_deref_mut(),
            "Warning: timeout value may have been modified (truncated?) due to number formatting issues;\n",
            Vec::new(),
        )?;
    }
    if timeout_set_error != 0 {
        eprint_call(
            heap,
            log_file.as_deref_mut(),
            "Warning: specified timeout value was ignored due to invalid number format, using the default;\n",
            Vec::new(),
        )?;
    }

    if options.hash_key != 0 {
        if ip.bAbcNumbers == 1 && ip.bCtPredecessors == 1 {
            eprint_call(
                heap,
                log_file.as_deref_mut(),
                "Terminating: generation of InChIKey is not available with 'Compress' option\n",
                Vec::new(),
            )?;
            return Ok(-1);
        }
        if ip.nInputType == tagInputType_INPUT_INCHI {
            eprint_call(
                heap,
                log_file.as_deref_mut(),
                "Terminating: generation of InChIKey is not available in InChI conversion mode\n",
                Vec::new(),
            )?;
            return Ok(-1);
        }
        if b_output_molfile_only == 1 {
            eprint_call(
                heap,
                log_file.as_deref_mut(),
                "Terminating: generation of InChIKey is not available with 'OutputSDF' option\n",
                Vec::new(),
            )?;
            return Ok(-1);
        }
    }

    if b_name_suffix != 0 {
        let mut first_path = SourceConstPointer::null();
        for path in ip.path {
            if !path.is_null() && heap.slice(path)?.first() != Some(&0) {
                first_path = path;
                break;
            }
        }
        let previous = match inchi__strdup(heap, first_path) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        let extensions = [".mol", ".txt", ".log", ".prb"];
        for (index, extension) in extensions.iter().enumerate() {
            if previous.is_null() {
                break;
            }
            let path_is_empty =
                ip.path[index].is_null() || heap.slice(ip.path[index])?.first() == Some(&0);
            if path_is_empty {
                let previous_bytes = source_c_string(heap, previous.as_const())?;
                let mut generated = previous_bytes[..previous_bytes.len() - 1].to_vec();
                generated.extend(extension.bytes().map(|byte| byte as i8));
                generated.push(0);
                let allocation = match inchi_malloc(
                    heap,
                    u64::try_from(generated.len())
                        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?,
                ) {
                    Ok(allocation) => allocation,
                    Err(SourceHeapError::AllocationFailed) => continue,
                    Err(error) => return Err(error),
                };
                heap.slice_mut(allocation)?.copy_from_slice(&generated);
                ip.num_paths += 1;
                if !ip.path[index].is_null() {
                    inchi_free(heap, ip.path[index].as_mut())?;
                }
                ip.path[index] = allocation.as_const();
            } else if option_literal(heap, ip.path[index], "NUL")? {
                inchi_free(heap, ip.path[index].as_mut())?;
                ip.path[index] = SourceConstPointer::null();
            }
        }
        inchi_free(heap, previous)?;
    }

    if ip.nInputType == tagInputType_INPUT_INCHI {
        b_output_molfile_only = 0;
        options.inchi_output_options |= INCHI_OUT_NO_AUX_INFO as i32;
        options.inchi_output_options &= !(INCHI_OUT_SHORT_AUX_INFO as i32);
        options.inchi_output_options &= !(INCHI_OUT_ONLY_AUX_INFO as i32);
        if options.inchi_output_options
            & (INCHI_OUT_SDFILE_ONLY | INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS) as i32
            == 0
        {
            options.inchi_output_options |= INCHI_OUT_PLAIN_TEXT as i32;
        }
    }

    if b_ver1_options != 0 {
        options.mode |= options.ver1_default_mode as i32;
    }
    ip.bDisplayEachComponentINChI = 0;
    ip.bCompareComponents = 0;
    ip.bMergeAllInputStructures = b_merge_all_input_structures;
    ip.bDoNotAddH = options.do_not_add_h;

    if options.mode == 0 || options.mode == REQ_MODE_STEREO as i32 {
        options.mode |=
            (REQ_MODE_BASIC | REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_NON_ISO | REQ_MODE_STEREO)
                as i32;
    } else {
        if options.mode & (REQ_MODE_BASIC | REQ_MODE_TAUT) as i32 == 0 {
            options.mode |= (REQ_MODE_BASIC | REQ_MODE_TAUT) as i32;
        }
        if options.mode & REQ_MODE_STEREO as i32 != 0
            && options.mode & (REQ_MODE_ISO | REQ_MODE_NON_ISO) as i32 == 0
        {
            options.mode |= (REQ_MODE_ISO | REQ_MODE_NON_ISO) as i32;
        }
    }
    if options.mode & REQ_MODE_ISO as i32 != 0 {
        options.mode |= REQ_MODE_NON_ISO as i32;
    }
    options.mode |=
        ((MIN_SB_RING_SIZE as i32) << REQ_MODE_MIN_SB_RING_SHFT) & REQ_MODE_MIN_SB_RING_MASK as i32;
    if ip.nInputType == tagInputType_INPUT_NONE && ip.num_paths > 0 {
        ip.nInputType = tagInputType_INPUT_MOLFILE;
    }
    ip.nMode = options.mode as INCHI_MODE;
    ip.bINChIOutputOptions = options.inchi_output_options
        | if b_output_molfile_only != 0 {
            INCHI_OUT_SDFILE_ONLY as i32
        } else {
            0
        };
    if b_output_molfile_only != 0 {
        b_output_style &= !((INCHI_OUT_PLAIN_TEXT
            | INCHI_OUT_PLAIN_TEXT_COMMENTS
            | INCHI_OUT_TABBED_OUTPUT) as i32);
        if b_output_molfile_dt != 0 {
            ip.bINChIOutputOptions |= INCHI_OUT_SDFILE_ATOMS_DT as i32;
        }
        if b_output_molfile_split != 0 {
            ip.bINChIOutputOptions |= INCHI_OUT_SDFILE_SPLIT as i32;
        }
    }
    ip.bINChIOutputOptions |= b_output_style;
    ip.bNoStructLabels = options.no_struct_labels;
    if options.forced_chiral_flag != 0 {
        ip.bChiralFlag = options.forced_chiral_flag;
    }

    ip.bTautFlags = (TG_FLAG_TEST_TAUT__ATOMS
        | TG_FLAG_DISCONNECT_SALTS
        | TG_FLAG_TEST_TAUT__SALTS
        | TG_FLAG_MOVE_POS_CHARGES
        | TG_FLAG_TEST_TAUT2_SALTS
        | TG_FLAG_MERGE_TAUT_SALTS
        | TG_FLAG_DISCONNECT_COORD
        | TG_FLAG_VARIABLE_PROTONS
        | TG_FLAG_HARD_ADD_REM_PROTONS) as INCHI_MODE;
    if options.reconnect_coord != 0 {
        ip.bTautFlags |= TG_FLAG_RECONNECT_COORD as INCHI_MODE;
    }
    if options.keto_enol_taut != 0 {
        ip.bTautFlags |= TG_FLAG_KETO_ENOL_TAUT as INCHI_MODE;
    }
    if options.taut_15_non_ring != 0 {
        ip.bTautFlags |= TG_FLAG_1_5_TAUT as INCHI_MODE;
    }
    for (enabled, flag) in [
        (options.pt_22_00_taut, TG_FLAG_PT_22_00),
        (options.pt_16_00_taut, TG_FLAG_PT_16_00),
        (options.pt_06_00_taut, TG_FLAG_PT_06_00),
        (options.pt_39_00_taut, TG_FLAG_PT_39_00),
        (options.pt_13_00_taut, TG_FLAG_PT_13_00),
        (options.pt_18_00_taut, TG_FLAG_PT_18_00),
        (options.pointed_edge_stereo, TG_FLAG_POINTED_EDGE_STEREO),
        (options.add_phosphine_stereo, TG_FLAG_PHOSPHINE_STEREO),
        (options.add_arsine_stereo, TG_FLAG_ARSINE_STEREO),
        (options.fix_sp3_bug, TG_FLAG_FIX_SP3_BUG),
    ] {
        if enabled != 0 {
            ip.bTautFlags |= flag as INCHI_MODE;
        }
    }
    ip.bTautFlagsDone = 0;
    if options.fix_fb2 != 0 {
        ip.bTautFlags |= (TG_FLAG_FIX_ISO_FIXEDH_BUG | TG_FLAG_FIX_TERM_H_CHRG_BUG) as INCHI_MODE;
        ip.bINChIOutputOptions |= INCHI_OUT_FIX_TRANSPOSITION_CHARGE_BUG as i32;
    }
    if ip.nInputType == tagInputType_INPUT_NONE {
        ip.nInputType = tagInputType_INPUT_MOLFILE;
    }
    if ip.nMode & REQ_MODE_STEREO as INCHI_MODE == 0 {
        ip.nMode &= !(REQ_MODE_DIFF_UU_STEREO as INCHI_MODE);
        ip.nMode &= !((REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU) as INCHI_MODE);
    }
    if options.std_format != 0 {
        ip.bINChIOutputOptions |= INCHI_OUT_STDINCHI as i32;
    }
    ip.bCalcInChIHash = if options.hash_key == 0 {
        tagInChIHashCalc_INCHIHASH_NONE as i32
    } else if options.hash_xtra1 != 0 && options.hash_xtra2 != 0 {
        tagInChIHashCalc_INCHIHASH_KEY_XTRA1_XTRA2 as i32
    } else if options.hash_xtra1 != 0 {
        tagInChIHashCalc_INCHIHASH_KEY_XTRA1 as i32
    } else if options.hash_xtra2 != 0 {
        tagInChIHashCalc_INCHIHASH_KEY_XTRA2 as i32
    } else {
        tagInChIHashCalc_INCHIHASH_KEY as i32
    };
    if options.hash_key == 0 && (options.hash_xtra1 != 0 || options.hash_xtra2 != 0) {
        eprint_call(
            heap,
            log_file.as_deref_mut(),
            "Hash extension(s) not generated: InChIKey not requested",
            Vec::new(),
        )?;
    }
    ip.bLargeMolecules = options.large_molecules;
    ip.bLooseTSACheck = options.loose_tsa_check;
    ip.bNPZz = options.np_zz;
    ip.bStereoAtZz = options.stereo_at_zz;
    ip.bNoWarnings = options.no_warnings;
    ip.bMergeHash = options.merge_hash;
    ip.bHideInChI = options.hide_inchi;
    ip.bPolymers = options.polymers;
    ip.bFoldPolymerSRU = options.fold_polymer_sru;
    ip.bFrameShiftScheme = options.frame_shift_scheme;
    ip.bINChIOutputOptions2 = options.inchi_output_options2;
    Ok(0)
}

enum PrintArgument {
    Signed(i64),
    Text(&'static str),
    Bytes(Vec<i8>),
}

fn eprint_call(
    heap: &mut SourceHeap,
    log_file: Option<&mut INCHI_IOSTREAM>,
    format: &str,
    arguments: Vec<PrintArgument>,
) -> Result<i32, SourceHeapError> {
    let Some(log_file) = log_file else {
        return Ok(-1);
    };
    let format_pointer = heap.allocate_model_storage(
        format
            .bytes()
            .chain(std::iter::once(0))
            .map(|byte| byte as i8)
            .collect(),
    )?;
    let mut text_pointers = Vec::new();
    let mut source_arguments = Vec::with_capacity(arguments.len());
    for argument in arguments {
        match argument {
            PrintArgument::Signed(value) => {
                source_arguments.push(SourceFormatArgument::Signed(value));
            }
            PrintArgument::Text(value) => {
                let pointer = heap.allocate_model_storage(
                    value
                        .bytes()
                        .chain(std::iter::once(0))
                        .map(|byte| byte as i8)
                        .collect(),
                )?;
                text_pointers.push(pointer);
                source_arguments.push(SourceFormatArgument::Bytes(pointer.as_const()));
            }
            PrintArgument::Bytes(value) => {
                let pointer = heap.allocate_model_storage(value)?;
                text_pointers.push(pointer);
                source_arguments.push(SourceFormatArgument::Bytes(pointer.as_const()));
            }
        }
    }
    let result = inchi_ios_eprint(
        heap,
        Some(log_file),
        format_pointer.as_const(),
        &SourceVaList {
            arguments: source_arguments,
            position: 0,
        },
    );
    heap.free(format_pointer)?;
    for pointer in text_pointers {
        heap.free(pointer)?;
    }
    result
}

fn nodisplay_call(
    heap: &mut SourceHeap,
    stream: &mut INCHI_IOSTREAM,
    stdout: SourceMutPointer<FILE>,
    format: &str,
    arguments: Vec<PrintArgument>,
) -> Result<i32, SourceHeapError> {
    let format_pointer = heap.allocate_model_storage(
        format
            .bytes()
            .chain(std::iter::once(0))
            .map(|byte| byte as i8)
            .collect(),
    )?;
    let mut text_pointers = Vec::new();
    let mut source_arguments = Vec::with_capacity(arguments.len());
    for argument in arguments {
        match argument {
            PrintArgument::Signed(value) => {
                source_arguments.push(SourceFormatArgument::Signed(value));
            }
            PrintArgument::Text(value) => {
                let pointer = heap.allocate_model_storage(
                    value
                        .bytes()
                        .chain(std::iter::once(0))
                        .map(|byte| byte as i8)
                        .collect(),
                )?;
                text_pointers.push(pointer);
                source_arguments.push(SourceFormatArgument::Bytes(pointer.as_const()));
            }
            PrintArgument::Bytes(value) => {
                let pointer = heap.allocate_model_storage(value)?;
                text_pointers.push(pointer);
                source_arguments.push(SourceFormatArgument::Bytes(pointer.as_const()));
            }
        }
    }
    let result = inchi_ios_print_nodisplay(
        heap,
        Some(stream),
        stdout,
        format_pointer.as_const(),
        &SourceVaList {
            arguments: source_arguments,
            position: 0,
        },
    );
    heap.free(format_pointer)?;
    for pointer in text_pointers {
        heap.free(pointer)?;
    }
    result
}

fn c_text(value: &str) -> Vec<i8> {
    value
        .bytes()
        .chain(std::iter::once(0))
        .map(|byte| byte as i8)
        .collect()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct InchiBuildMetadata<'a> {
    pub(crate) compiler: &'a str,
    pub(crate) date: &'a str,
    pub(crate) time: &'a str,
}

fn pointer_starts_with_nonzero(
    heap: &SourceHeap,
    pointer: SourceConstPointer<i8>,
) -> Result<bool, SourceHeapError> {
    if pointer.is_null() {
        return Ok(false);
    }
    Ok(heap.slice(pointer)?.first().is_some_and(|byte| *byte != 0))
}

pub(crate) fn PrintInputParms(
    heap: &mut SourceHeap,
    mut log_file: Option<&mut INCHI_IOSTREAM>,
    ip: &INPUT_PARMS,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c:2130 PrintInputParms
    // INCHI✔️❌: int PrintInputParms(INCHI_IOSTREAM* log_file,
    // INCHI✔️❌:     INPUT_PARMS* ip)
    // INCHI✔️❌: {
    // INCHI✔️❌:
    // INCHI✔️❌: #if (BUILD_WITH_ENG_OPTIONS==1)
    // INCHI✔️❌:     const int developer_options = 1;
    // INCHI✔️❌: #else
    // INCHI✔️❌:     const int developer_options = 0;
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔️❌:     int bInChI2Struct = 0; /* winchi-1 can not convert InChI to structure */
    // INCHI✔️❌: #else
    // INCHI✔️❌:     int bInChI2Struct = (ip->bReadInChIOptions & READ_INCHI_TO_STRUCTURE) && ip->nInputType == INPUT_INCHI;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_MODE nMode = ip->nMode;
    // INCHI✔️❌:     int k;
    // INCHI✔️❌:     int bStdFormat = 1;
    // INCHI✔️❌:     int first = 1;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!(ip->bINChIOutputOptions & INCHI_OUT_STDINCHI))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         bStdFormat = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* Some stereo */
    // INCHI✔️❌:     if (!(nMode & REQ_MODE_STEREO))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Using specific structure perception features:\n");
    // INCHI✔️❌:         first = 0;
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "  Stereo OFF\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (!(TG_FLAG_POINTED_EDGE_STEREO & ip->bTautFlags))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (first)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_eprint(log_file, "Using specific structure perception features:\n");
    // INCHI✔️❌:                 first = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Both ends of wedge point to stereocenters\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (ip->bDoNotAddH)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (first)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "Using specific structure perception features:\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "  Do not add H\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( RENUMBER_ATOMS_AND_RECALC_V106 == 1 )
    // INCHI✔️❌:     if (ip->bRenumber == 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "\nGenerate InChI upon random atom renumbering\n\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( UNDERIVATIZE == 1 )
    // INCHI✔️❌:     if (ip->bUnderivatize == 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "\nConvert input structure to derivative precursor before InChI calculation\n\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else if (ip->bUnderivatize == 3)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "\nOutputs derivative information for the input structure\n\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  Generation/conversion indicator */
    // INCHI✔️❌:     if (bStdFormat)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (!(ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY) && !bInChI2Struct)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "Generating standard InChI\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( !defined(TARGET_API_LIB) && !defined(TARGET_LIB_FOR_WINCHI) && !defined(TARGET_EXE_USING_API) )
    // INCHI✔️❌:         /* effective only in command line program InChI or stdInChI */
    // INCHI✔️❌:         else if (bInChI2Struct)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "Converting InChI(s) to structure(s) in %s\n",
    // INCHI✔️❌:                 (ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY) ?
    // INCHI✔️❌:                 "MOL format" : "aux. info format");
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Generating non-standard InChI with the options: \n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* SDfile output */
    // INCHI✔️❌:     if (ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file,
    // INCHI✔️❌:             "Output SDfile only without stereochemical information and atom coordinates%s\n",
    // INCHI✔️❌:             (ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ATOMS_DT) ?
    // INCHI✔️❌:             "\n(write H isotopes as D, T)" : "");
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Fixed/Mobile H */
    // INCHI✔️❌:     if (!bStdFormat)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((nMode & (REQ_MODE_BASIC | REQ_MODE_TAUT)) == (REQ_MODE_BASIC | REQ_MODE_TAUT))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Mobile H Perception OFF (include FixedH layer)\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if ((nMode & (REQ_MODE_BASIC | REQ_MODE_TAUT)) == (REQ_MODE_TAUT))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Mobile H Perception ON  (omit FixedH layer)\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if ((nMode & (REQ_MODE_BASIC | REQ_MODE_TAUT)) == (REQ_MODE_BASIC))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Mobile H ignored\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Undefined Mobile H mode\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if ((ip->bTautFlags & TG_FLAG_VARIABLE_PROTONS))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!(ip->bTautFlags & TG_FLAG_HARD_ADD_REM_PROTONS))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_eprint(log_file, "  Disabled Aggressive (De)protonation\n");
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( FIND_RING_SYSTEMS != 1 )
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "  %s5-, 6-, 7-memb. ring taut. ignored\n", i ? "; " : "");
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         /* RecMet */
    // INCHI✔️❌:         if (ip->bTautFlags & TG_FLAG_DISCONNECT_COORD)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (ip->bTautFlags & TG_FLAG_RECONNECT_COORD)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_eprint(log_file, "  Include bonds to metals\n");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_eprint(log_file, "  Do not reconnect metals (omit RecMet layer)\n");
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Do not disconnect metals\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Isotopic - always ON, output disabled. 09-17-2009*/
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         if ( nMode & REQ_MODE_ISO )
    // INCHI✔️❌:         inchi_ios_eprint( log_file, "  Isotopic ON\n");
    // INCHI✔️❌:         else if ( nMode & REQ_MODE_NON_ISO )
    // INCHI✔️❌:         inchi_ios_eprint( log_file, "  Isotopic OFF\n");
    // INCHI✔️❌:         */
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( FIX_ADJ_RAD == 1 )
    // INCHI✔️❌:         if (ip->bTautFlags & TG_FLAG_FIX_ADJ_RADICALS)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "Fix Adjacent Radicals\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Stereo */
    // INCHI✔️❌:         if (nMode & REQ_MODE_STEREO)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  %s%s%s%sStereo ON\n",
    // INCHI✔️❌:                 (nMode & REQ_MODE_NOEQ_STEREO) ? "Slow " : "",
    // INCHI✔️❌:                 (nMode & REQ_MODE_REDNDNT_STEREO) ? "Redund. " : "",
    // INCHI✔️❌:                 (nMode & REQ_MODE_NO_ALT_SBONDS) ? "No AltBond " : "",
    // INCHI✔️❌:
    // INCHI✔️❌:                 (nMode & REQ_MODE_RACEMIC_STEREO) ? "Racemic " :
    // INCHI✔️❌:                 (nMode & REQ_MODE_RELATIVE_STEREO) ? "Relative " :
    // INCHI✔️❌:                 (nMode & REQ_MODE_CHIR_FLG_STEREO) ? "Chiral Flag " : "Absolute ");
    // INCHI✔️❌:
    // INCHI✔️❌:             if (0 == (nMode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_eprint(log_file, "  Include undefined/unknown stereogenic centers and bonds\n");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else if (REQ_MODE_SC_IGN_ALL_UU == (nMode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_eprint(log_file, "  Omit undefined/unknown stereogenic centers\n");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else if (REQ_MODE_SB_IGN_ALL_UU == (nMode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_eprint(log_file, "  Omit undefined/unknown stereogenic bonds\n");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /*case REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU*/
    // INCHI✔️❌:                 inchi_ios_eprint(log_file, "  Omit undefined/unknown stereogenic centers and bonds\n");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (0 != (nMode & REQ_MODE_DIFF_UU_STEREO))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_eprint(log_file, "  Make labels for unknown and undefined stereo different\n");
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( defined(MIN_SB_RING_SIZE) && MIN_SB_RING_SIZE > 0 )
    // INCHI✔️❌:             k = (ip->nMode & REQ_MODE_MIN_SB_RING_MASK) >> REQ_MODE_MIN_SB_RING_SHFT;
    // INCHI✔️❌:             if (bRELEASE_VERSION != 1 || k != MIN_SB_RING_SIZE)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (k >= 3)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     inchi_ios_eprint(log_file, "  Min. stereobond ring size: %d\n", k);
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     inchi_ios_eprint(log_file, "  Min. stereobond ring size: NONE\n");
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         }   /* Stereo */
    // INCHI✔️❌:     }   /* !bStdFormat */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!bStdFormat)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (TG_FLAG_KETO_ENOL_TAUT & ip->bTautFlags)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Account for keto-enol tautomerism\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Do not account for keto-enol tautomerism\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (TG_FLAG_1_5_TAUT & ip->bTautFlags)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Account for 1,5-tautomerism\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Do not account for 1,5-tautomerism\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (TG_FLAG_PT_22_00 & ip->bTautFlags)
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Account for PT_22_00 tautomerism\n");
    // INCHI✔️❌:         else
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Do not account for PT_22_00 tautomerism\n");
    // INCHI✔️❌:         if (TG_FLAG_PT_16_00 & ip->bTautFlags)
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Account for PT_16_00 tautomerism\n");
    // INCHI✔️❌:         else
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Do not account for PT_16_00 tautomerism\n");
    // INCHI✔️❌:         if (TG_FLAG_PT_06_00 & ip->bTautFlags)
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Account for PT_06_00 tautomerism\n");
    // INCHI✔️❌:         else
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Do not account for PT_06_00 tautomerism\n");
    // INCHI✔️❌:         if (TG_FLAG_PT_39_00 & ip->bTautFlags)
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Account for PT_39_00 tautomerism\n");
    // INCHI✔️❌:         else
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Do not account for PT_39_00 tautomerism\n");
    // INCHI✔️❌:         if (TG_FLAG_PT_13_00 & ip->bTautFlags)
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Account for PT_13_00 tautomerism\n");
    // INCHI✔️❌:         else
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Do not account for PT_13_00 tautomerism\n");
    // INCHI✔️❌:         if (TG_FLAG_PT_18_00 & ip->bTautFlags)
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Account for PT_18_00 tautomerism\n");
    // INCHI✔️❌:         else
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  Do not account for PT_18_00 tautomerism\n");
    // INCHI✔️❌:
    // INCHI✔️❌:         if (developer_options)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (TG_FLAG_PHOSPHINE_STEREO & ip->bTautFlags)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_eprint(log_file, "  Include phosphine stereochemistry\n");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_eprint(log_file, "  Do not include phosphine stereochemistry\n");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (TG_FLAG_ARSINE_STEREO & ip->bTautFlags)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_eprint(log_file, "  Include arsine stereochemistry\n");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_eprint(log_file, "  Do not include arsine stereochemistry\n");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (!(TG_FLAG_FIX_SP3_BUG & ip->bTautFlags))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_eprint(log_file, "  Turned OFF fix of bug leading to missing or undefined sp3 parity\n");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (!(TG_FLAG_FIX_ISO_FIXEDH_BUG & ip->bTautFlags))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_eprint(log_file, "  Turned OFF bug-fixes found after v.1.02b release\n");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (!(ip->bFixNonUniformDraw))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_eprint(log_file, "  Turned OFF fixes of non-uniform drawing issues\n");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (!(TG_FLAG_MOVE_POS_CHARGES & ip->bTautFlags))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_eprint(log_file, "  MovePos turned OFF\n");
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:     } /* !bStdFormat */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ip->bCalcInChIHash != INCHIHASH_NONE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (bStdFormat)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "Generating standard InChIKey\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "Generating InChIKey\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (ip->bCalcInChIHash == INCHIHASH_KEY_XTRA1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "Generating hash extension (1st block)\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (ip->bCalcInChIHash == INCHIHASH_KEY_XTRA2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "Generating hash extension (2nd block)\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (ip->bCalcInChIHash == INCHIHASH_KEY_XTRA1_XTRA2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "Generating hash extension (two blocks)\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ip->bINChIOutputOptions & INCHI_OUT_SAVEOPT)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Saving InChI creation options");
    // INCHI✔️❌:         if (bStdFormat)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, " suppressed for standard InChI");
    // INCHI✔️❌:             /* NB: actual suppression takes place on InChI serialization */
    // INCHI✔️❌:             /* (as on e.g. Inchi2Inchi conversion it may appear that we create non-std */
    // INCHI✔️❌:             /*  InChI instead of standard one) */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ip->bAllowEmptyStructure)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Issue warning on empty structure\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Input */
    // INCHI✔️❌:     if (ip->nInputType)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Input format: %s",
    // INCHI✔️❌:             ip->nInputType == INPUT_MOLFILE ? "MOLfile" :
    // INCHI✔️❌:             ip->nInputType == INPUT_SDFILE ? "SDfile" :
    // INCHI✔️❌: #if ( READ_INCHI_STRING == 1 )
    // INCHI✔️❌:             ip->nInputType == INPUT_INCHI ? "InChI (plain identifier)" :
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             ip->nInputType == INPUT_INCHI_PLAIN ? "InChI AuxInfo (plain)" : "Unknown");
    // INCHI✔️❌:         if ((ip->nInputType == INPUT_MOLFILE || ip->nInputType == INPUT_SDFILE) &&
    // INCHI✔️❌:             ip->bGetMolfileNumber)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "  (attempting to read Molfile number)");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ip->szSdfDataHeader[0] && ip->nInputType != INPUT_SDFILE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "  SDfile data header: \"%s\"\n", ip->szSdfDataHeader);
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Output */
    // INCHI✔️❌:     inchi_ios_eprint(log_file, "Output format: %s%s\n",
    // INCHI✔️❌:         (ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT) ? "Plain text" :
    // INCHI✔️❌:
    // INCHI✔️❌:         ((ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY) && bInChI2Struct) ? "SDfile only (without stereochemical info and atom coordinates)" :
    // INCHI✔️❌:         ((ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY) && !bInChI2Struct) ? "SDfile only" : "Unknown",
    // INCHI✔️❌:
    // INCHI✔️❌:         ((ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT) &&
    // INCHI✔️❌:             (ip->bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT)) ? ", tabbed" : "");
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( bRELEASE_VERSION == 1 )
    // INCHI✔️❌:     if (ip->bCtPredecessors || ip->bAbcNumbers)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (ip->bCtPredecessors && ip->bAbcNumbers)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "Representation: Compressed\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "Connection table: %s, %s\n",
    // INCHI✔️❌:                 ip->bCtPredecessors ? "Predecessor_numbers(closures)" : "Canon_numbers(branching, ring closures)",
    // INCHI✔️❌:                 ip->bAbcNumbers ? "Shorter alternative" : "Numerical");
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #else
    // INCHI✔️❌:     if ((bRELEASE_VERSION != 1) || ip->bCtPredecessors || ip->bAbcNumbers)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Connection table: %s, %s\n",
    // INCHI✔️❌:             ip->bCtPredecessors ? "Predecessor_numbers(closures)" : "Canon_numbers(branching, ring closures)",
    // INCHI✔️❌:             ip->bAbcNumbers ? "Shorter alternative" : "Numerical");
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Representation: Numerical");
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ip->bNoWarnings)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Warnings suppressed\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ip->bHideInChI)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Printing InChI string itself suppressed\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (ip->bMergeHash)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "InChIKey combined with extra hash(es)\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!(ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (ip->bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "Aux. info suppressed\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (ip->bINChIOutputOptions & INCHI_OUT_SHORT_AUX_INFO)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "Minimal Aux. info\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "Full Aux. info\n");
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (ip->first_struct_number > 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Skipping %ld structure%s\n", ip->first_struct_number - 1, ip->first_struct_number == 2 ? "" : "s");
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (ip->last_struct_number > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Terminate after structure #%ld\n", ip->last_struct_number);
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (ip->bSaveWarningStructsAsProblem && ip->path[3] && ip->path[3][0])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Saving warning structures into the problem file\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (ip->bSaveAllGoodStructsAsProblem && ip->path[3] && ip->path[3][0])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Saving only all good structures into the problem file\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ip->bINChIOutputOptions2 & INCHI_OUT_INCHI_GEN_ERROR)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Print empty InChI if generation fails\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (ip->bINChIOutputOptions2 & INCHI_OUT_MISMATCH_AS_ERROR)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Consider problem/mismatch on InChI conversion as error\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (ip->msec_MaxTime)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Timeout per structure: %ld msec\n", ip->msec_MaxTime);
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         unsigned long seconds = ip->msec_MaxTime / 1000;
    // INCHI✔️❌:         unsigned long milliseconds = (ip->msec_MaxTime%1000);
    // INCHI✔️❌:         inchi_ios_eprint( log_file, "Timeout per structure: %lu/*.%03lu sec\n", seconds, milliseconds); -- djb-rwth: ignoring LLVM warning
    // INCHI✔️❌:         inchi_ios_eprint( log_file, "Timeout per structure: %lu sec\n", seconds );
    // INCHI✔️❌:         */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "No timeout\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ip->bLooseTSACheck)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Relax criteria of ambiguous drawing for in-ring stereo centers\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int maxna = NORMALLY_ALLOWED_INP_MAX_ATOMS;
    // INCHI✔️❌:         if (ip->bLargeMolecules)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "Experimental mode: ");
    // INCHI✔️❌:             maxna = MAX_ATOMS;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Up to %d atoms per structure\n", maxna);
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ip->bPolymers != POLYMERS_NO)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Experimental mode: Treating polymers");
    // INCHI✔️❌:
    // INCHI✔️❌:         if (ip->bPolymers == POLYMERS_MODERN)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (ip->bPolymers == POLYMERS_LEGACY)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, " (v. 1.05 legacy mode)");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (ip->bPolymers == POLYMERS_LEGACY_PLUS)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, " (v. 1.05 legacy mode with senior link placed at start)");
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (ip->bFoldPolymerSRU)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(log_file, "; CRU folding enabled");
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     inchi_ios_eprint(log_file, "\n");
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ip->bNPZz == 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Allowing non-polymer Zz pseudo atoms\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (ip->bStereoAtZz == 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Allowing stereo at atoms connected to Zz\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  Report debug modes */
    // INCHI✔️❌: #if ( bRELEASE_VERSION != 1 )
    // INCHI✔️❌:     inchi_ios_eprint(log_file, "Release version = NO\n");
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #if ( TRACE_MEMORY_LEAKS == 1 && defined(_DEBUG) )
    // INCHI✔️❌:     inchi_ios_eprint(log_file, "Tracing memory leaks (SLOW)\n");
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: #if (BUILD_WITH_ENG_OPTIONS==1)
    // INCHI✔️❌:     inchi_ios_eprint(log_file, "! Working in engineering mode\n");
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     inchi_ios_eprint(log_file, "\n");
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( bRELEASE_VERSION != 1 )
    // INCHI✔️❌: #if ( FIND_RING_SYSTEMS == 1 )
    // INCHI✔️❌:     inchi_ios_eprint(log_file, "Find ring systems=Y\nTautomers:\n");
    // INCHI✔️❌:     inchi_ios_eprint(log_file, " 4-pyridinol=%s\n", TAUT_4PYRIDINOL_RINGS == 1 ? "Y" : "N");
    // INCHI✔️❌:     inchi_ios_eprint(log_file, " pyrazole=%s\n", TAUT_PYRAZOLE_RINGS == 1 ? "Y" : "N");
    // INCHI✔️❌:     inchi_ios_eprint(log_file, " tropolone=%s\n", TAUT_TROPOLONE_7 == 1 ? "Y" : "N");
    // INCHI✔️❌:     inchi_ios_eprint(log_file, " tropolone-5=%s\n", TAUT_TROPOLONE_5 == 1 ? "Y" : "N");
    // INCHI✔️❌:     inchi_ios_eprint(log_file, "Only chain attachments to tautomeric rings=%s\n", TAUT_RINGS_ATTACH_CHAIN == 1 ? "Y" : "N");
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ip->bGetSdfileId)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint(log_file, "Extracting SDfile IDs\n");
    // INCHI✔️❌:     }
    // INCHI✔️❌:     inchi_ios_eprint(log_file, "\nDbg: MOVE_CHARGES=%d\n",
    // INCHI✔️❌:         0 != (ip->bTautFlags & TG_FLAG_MOVE_POS_CHARGES));
    // INCHI✔️❌:     inchi_ios_eprint(log_file, "     REPLACE_ALT_WITH_TAUT=%d; NEUTRALIZE_ENDPOINTS=%d; BNS_PROTECT_FROM_TAUT=%d\n",
    // INCHI✔️❌:         REPLACE_ALT_WITH_TAUT, NEUTRALIZE_ENDPOINTS, BNS_PROTECT_FROM_TAUT);
    // INCHI✔️❌:     inchi_ios_eprint(log_file, "     DISCONNECT_SALTS=%d;   TEST_TAUT_SALTS=%d;    TEST_TAUT2_SALTS=%d\n",
    // INCHI✔️❌:         0 != (ip->bTautFlags & TG_FLAG_DISCONNECT_SALTS),
    // INCHI✔️❌:         0 != (ip->bTautFlags & TG_FLAG_TEST_TAUT__SALTS),
    // INCHI✔️❌:         0 != (ip->bTautFlags & TG_FLAG_TEST_TAUT2_SALTS));
    // INCHI✔️❌:     inchi_ios_eprint(log_file, "     CHARGED_ACID_TAUT_ONLY=%d MERGE_TAUT_SALTS=%d\n",
    // INCHI✔️❌:         0 == (ip->bTautFlags & TG_FLAG_ALLOW_NO_NEGTV_O),
    // INCHI✔️❌:         0 != (ip->bTautFlags & TG_FLAG_MERGE_TAUT_SALTS));
    // INCHI✔️❌:     inchi_ios_eprint(log_file, "     DISCONNECT_COORD=%d\n",
    // INCHI✔️❌:         0 != (ip->bTautFlags & TG_FLAG_DISCONNECT_COORD));
    // INCHI✔️❌: #endif /* ( bRELEASE_VERSION != 1 ) */
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: PrintInputParms

    macro_rules! print {
        ($format:expr $(, $argument:expr)* $(,)?) => {{
            eprint_call(
                heap,
                log_file.as_deref_mut(),
                $format,
                vec![$($argument),*],
            )?;
        }};
    }

    let b_inchi_to_struct = ip.bReadInChIOptions & READ_INCHI_TO_STRUCTURE as i32 != 0
        && ip.nInputType == tagInputType_INPUT_INCHI;
    let n_mode = ip.nMode;
    let mut standard_format = true;
    let mut first = true;

    if ip.bINChIOutputOptions & INCHI_OUT_STDINCHI as i32 == 0 {
        standard_format = false;
    }
    if n_mode & REQ_MODE_STEREO as u64 == 0 {
        print!("Using specific structure perception features:\n");
        first = false;
        print!("  Stereo OFF\n");
    } else if ip.bTautFlags & TG_FLAG_POINTED_EDGE_STEREO as u64 == 0 {
        if first {
            print!("Using specific structure perception features:\n");
            first = false;
        }
        print!("  Both ends of wedge point to stereocenters\n");
    }
    if ip.bDoNotAddH != 0 {
        if first {
            print!("Using specific structure perception features:\n");
        }
        print!("  Do not add H\n");
    }
    if ip.bRenumber == 1 {
        print!("\nGenerate InChI upon random atom renumbering\n\n");
    }
    if ip.bUnderivatize == 1 {
        print!("\nConvert input structure to derivative precursor before InChI calculation\n\n");
    } else if ip.bUnderivatize == 3 {
        print!("\nOutputs derivative information for the input structure\n\n");
    }

    if standard_format {
        if ip.bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY as i32 == 0 && !b_inchi_to_struct {
            print!("Generating standard InChI\n");
        }
    } else {
        print!("Generating non-standard InChI with the options: \n");
    }
    if ip.bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY as i32 != 0 {
        print!(
            "Output SDfile only without stereochemical information and atom coordinates%s\n",
            PrintArgument::Text(
                if ip.bINChIOutputOptions & INCHI_OUT_SDFILE_ATOMS_DT as i32 != 0 {
                    "\n(write H isotopes as D, T)"
                } else {
                    ""
                }
            )
        );
    }

    if !standard_format {
        match n_mode & (REQ_MODE_BASIC | REQ_MODE_TAUT) as u64 {
            value if value == (REQ_MODE_BASIC | REQ_MODE_TAUT) as u64 => {
                print!("  Mobile H Perception OFF (include FixedH layer)\n");
            }
            value if value == REQ_MODE_TAUT as u64 => {
                print!("  Mobile H Perception ON  (omit FixedH layer)\n");
            }
            value if value == REQ_MODE_BASIC as u64 => print!("  Mobile H ignored\n"),
            _ => print!("  Undefined Mobile H mode\n"),
        }
        if ip.bTautFlags & TG_FLAG_VARIABLE_PROTONS as u64 != 0
            && ip.bTautFlags & TG_FLAG_HARD_ADD_REM_PROTONS as u64 == 0
        {
            print!("  Disabled Aggressive (De)protonation\n");
        }
        if ip.bTautFlags & TG_FLAG_DISCONNECT_COORD as u64 != 0 {
            if ip.bTautFlags & TG_FLAG_RECONNECT_COORD as u64 != 0 {
                print!("  Include bonds to metals\n");
            } else {
                print!("  Do not reconnect metals (omit RecMet layer)\n");
            }
        } else {
            print!("  Do not disconnect metals\n");
        }
        if n_mode & REQ_MODE_STEREO as u64 != 0 {
            let stereo_type = if n_mode & REQ_MODE_RACEMIC_STEREO as u64 != 0 {
                "Racemic "
            } else if n_mode & REQ_MODE_RELATIVE_STEREO as u64 != 0 {
                "Relative "
            } else if n_mode & REQ_MODE_CHIR_FLG_STEREO as u64 != 0 {
                "Chiral Flag "
            } else {
                "Absolute "
            };
            print!(
                "  %s%s%s%sStereo ON\n",
                PrintArgument::Text(if n_mode & REQ_MODE_NOEQ_STEREO as u64 != 0 {
                    "Slow "
                } else {
                    ""
                }),
                PrintArgument::Text(if n_mode & REQ_MODE_REDNDNT_STEREO as u64 != 0 {
                    "Redund. "
                } else {
                    ""
                }),
                PrintArgument::Text(if n_mode & REQ_MODE_NO_ALT_SBONDS as u64 != 0 {
                    "No AltBond "
                } else {
                    ""
                }),
                PrintArgument::Text(stereo_type)
            );
            let undefined_mode = n_mode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU) as u64;
            match undefined_mode {
                0 => print!("  Include undefined/unknown stereogenic centers and bonds\n"),
                value if value == REQ_MODE_SC_IGN_ALL_UU as u64 => {
                    print!("  Omit undefined/unknown stereogenic centers\n");
                }
                value if value == REQ_MODE_SB_IGN_ALL_UU as u64 => {
                    print!("  Omit undefined/unknown stereogenic bonds\n");
                }
                _ => print!("  Omit undefined/unknown stereogenic centers and bonds\n"),
            }
            if n_mode & REQ_MODE_DIFF_UU_STEREO as u64 != 0 {
                print!("  Make labels for unknown and undefined stereo different\n");
            }
            let ring_size =
                ((ip.nMode & REQ_MODE_MIN_SB_RING_MASK as u64) >> REQ_MODE_MIN_SB_RING_SHFT) as i32;
            if ring_size != MIN_SB_RING_SIZE as i32 {
                if ring_size >= 3 {
                    print!(
                        "  Min. stereobond ring size: %d\n",
                        PrintArgument::Signed(i64::from(ring_size))
                    );
                } else {
                    print!("  Min. stereobond ring size: NONE\n");
                }
            }
        }
    }

    if !standard_format {
        for (flag, enabled, disabled) in [
            (
                TG_FLAG_KETO_ENOL_TAUT,
                "  Account for keto-enol tautomerism\n",
                "  Do not account for keto-enol tautomerism\n",
            ),
            (
                TG_FLAG_1_5_TAUT,
                "  Account for 1,5-tautomerism\n",
                "  Do not account for 1,5-tautomerism\n",
            ),
            (
                TG_FLAG_PT_22_00,
                "  Account for PT_22_00 tautomerism\n",
                "  Do not account for PT_22_00 tautomerism\n",
            ),
            (
                TG_FLAG_PT_16_00,
                "  Account for PT_16_00 tautomerism\n",
                "  Do not account for PT_16_00 tautomerism\n",
            ),
            (
                TG_FLAG_PT_06_00,
                "  Account for PT_06_00 tautomerism\n",
                "  Do not account for PT_06_00 tautomerism\n",
            ),
            (
                TG_FLAG_PT_39_00,
                "  Account for PT_39_00 tautomerism\n",
                "  Do not account for PT_39_00 tautomerism\n",
            ),
            (
                TG_FLAG_PT_13_00,
                "  Account for PT_13_00 tautomerism\n",
                "  Do not account for PT_13_00 tautomerism\n",
            ),
            (
                TG_FLAG_PT_18_00,
                "  Account for PT_18_00 tautomerism\n",
                "  Do not account for PT_18_00 tautomerism\n",
            ),
        ] {
            print!(if ip.bTautFlags & flag as u64 != 0 {
                enabled
            } else {
                disabled
            });
        }
    }

    if ip.bCalcInChIHash != tagInChIHashCalc_INCHIHASH_NONE as i32 {
        print!(if standard_format {
            "Generating standard InChIKey\n"
        } else {
            "Generating InChIKey\n"
        });
        match ip.bCalcInChIHash as u32 {
            tagInChIHashCalc_INCHIHASH_KEY_XTRA1 => {
                print!("Generating hash extension (1st block)\n")
            }
            tagInChIHashCalc_INCHIHASH_KEY_XTRA2 => {
                print!("Generating hash extension (2nd block)\n")
            }
            tagInChIHashCalc_INCHIHASH_KEY_XTRA1_XTRA2 => {
                print!("Generating hash extension (two blocks)\n")
            }
            _ => {}
        }
    }
    if ip.bINChIOutputOptions & INCHI_OUT_SAVEOPT as i32 != 0 {
        print!("Saving InChI creation options");
        if standard_format {
            print!(" suppressed for standard InChI");
        }
        print!("\n");
    }
    if ip.bAllowEmptyStructure != 0 {
        print!("Issue warning on empty structure\n");
    }

    if ip.nInputType != tagInputType_INPUT_NONE {
        let input_name = match ip.nInputType {
            tagInputType_INPUT_MOLFILE => "MOLfile",
            tagInputType_INPUT_SDFILE => "SDfile",
            tagInputType_INPUT_INCHI => "InChI (plain identifier)",
            tagInputType_INPUT_INCHI_PLAIN => "InChI AuxInfo (plain)",
            _ => "Unknown",
        };
        print!("Input format: %s", PrintArgument::Text(input_name));
        if matches!(
            ip.nInputType,
            tagInputType_INPUT_MOLFILE | tagInputType_INPUT_SDFILE
        ) && ip.bGetMolfileNumber != 0
        {
            print!("  (attempting to read Molfile number)");
        }
        print!("\n");
    }
    if ip.szSdfDataHeader[0] != 0 && ip.nInputType != tagInputType_INPUT_SDFILE {
        print!(
            "  SDfile data header: \"%s\"\n",
            PrintArgument::Bytes(ip.szSdfDataHeader.to_vec())
        );
    }

    let output_name = if ip.bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT as i32 != 0 {
        "Plain text"
    } else if ip.bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY as i32 != 0 && b_inchi_to_struct {
        "SDfile only (without stereochemical info and atom coordinates)"
    } else if ip.bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY as i32 != 0 {
        "SDfile only"
    } else {
        "Unknown"
    };
    let tabbed = if ip.bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT as i32 != 0
        && ip.bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT as i32 != 0
    {
        ", tabbed"
    } else {
        ""
    };
    print!(
        "Output format: %s%s\n",
        PrintArgument::Text(output_name),
        PrintArgument::Text(tabbed)
    );

    if ip.bCtPredecessors != 0 || ip.bAbcNumbers != 0 {
        if ip.bCtPredecessors != 0 && ip.bAbcNumbers != 0 {
            print!("Representation: Compressed\n");
        } else {
            print!(
                "Connection table: %s, %s\n",
                PrintArgument::Text(if ip.bCtPredecessors != 0 {
                    "Predecessor_numbers(closures)"
                } else {
                    "Canon_numbers(branching, ring closures)"
                }),
                PrintArgument::Text(if ip.bAbcNumbers != 0 {
                    "Shorter alternative"
                } else {
                    "Numerical"
                })
            );
        }
    }
    if ip.bNoWarnings != 0 {
        print!("Warnings suppressed\n");
    }
    if ip.bHideInChI != 0 {
        print!("Printing InChI string itself suppressed\n");
    }
    if ip.bMergeHash != 0 {
        print!("InChIKey combined with extra hash(es)\n");
    }

    if ip.bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY as i32 == 0 {
        if ip.bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO as i32 != 0 {
            print!("Aux. info suppressed\n");
        } else if ip.bINChIOutputOptions & INCHI_OUT_SHORT_AUX_INFO as i32 != 0 {
            print!("Minimal Aux. info\n");
        } else {
            print!("Full Aux. info\n");
        }
    }
    if ip.first_struct_number > 1 {
        print!(
            "Skipping %ld structure%s\n",
            PrintArgument::Signed(ip.first_struct_number.wrapping_sub(1)),
            PrintArgument::Text(if ip.first_struct_number == 2 { "" } else { "s" })
        );
    }
    if ip.last_struct_number > 0 {
        print!(
            "Terminate after structure #%ld\n",
            PrintArgument::Signed(ip.last_struct_number)
        );
    }
    let problem_path_present = pointer_starts_with_nonzero(heap, ip.path[3])?;
    if ip.bSaveWarningStructsAsProblem != 0 && problem_path_present {
        print!("Saving warning structures into the problem file\n");
    }
    if ip.bSaveAllGoodStructsAsProblem != 0 && problem_path_present {
        print!("Saving only all good structures into the problem file\n");
    }
    if ip.bINChIOutputOptions2 & INCHI_OUT_INCHI_GEN_ERROR as i32 != 0 {
        print!("Print empty InChI if generation fails\n");
    }
    if ip.bINChIOutputOptions2 & INCHI_OUT_MISMATCH_AS_ERROR as i32 != 0 {
        print!("Consider problem/mismatch on InChI conversion as error\n");
    }
    if ip.msec_MaxTime != 0 {
        print!(
            "Timeout per structure: %ld msec\n",
            PrintArgument::Signed(ip.msec_MaxTime)
        );
    } else {
        print!("No timeout\n");
    }
    if ip.bLooseTSACheck != 0 {
        print!("Relax criteria of ambiguous drawing for in-ring stereo centers\n");
    }
    let max_atoms = if ip.bLargeMolecules != 0 {
        print!("Experimental mode: ");
        MAX_ATOMS
    } else {
        NORMALLY_ALLOWED_INP_MAX_ATOMS
    };
    print!(
        "Up to %d atoms per structure\n",
        PrintArgument::Signed(i64::from(max_atoms))
    );

    if ip.bPolymers != POLYMERS_NO as i32 {
        print!("Experimental mode: Treating polymers");
        if ip.bPolymers == POLYMERS_LEGACY as i32 {
            print!(" (v. 1.05 legacy mode)");
        } else if ip.bPolymers == POLYMERS_LEGACY_PLUS as i32 {
            print!(" (v. 1.05 legacy mode with senior link placed at start)");
        }
        if ip.bFoldPolymerSRU != 0 {
            print!("; CRU folding enabled");
        }
    }
    print!("\n");
    if ip.bNPZz == 1 {
        print!("Allowing non-polymer Zz pseudo atoms\n");
    }
    if ip.bStereoAtZz == 1 {
        print!("Allowing stereo at atoms connected to Zz\n");
    }
    print!("\n");
    Ok(0)
}

pub(crate) fn HelpCommandLineParms(
    heap: &mut SourceHeap,
    stream: Option<&mut INCHI_IOSTREAM>,
    stdout: SourceMutPointer<FILE>,
    build: InchiBuildMetadata<'_>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c:2700 HelpCommandLineParms
    // INCHI✔️❌: void HelpCommandLineParms(INCHI_IOSTREAM* f)
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (!f)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( bRELEASE_VERSION == 1 )
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f,
    // INCHI✔️❌: #ifdef TARGET_EXE_USING_API
    // INCHI✔️❌:         "%s %-s\n%-s Build (%-s%-s) of %s %-s %-s\n\nUsage:\ninchi_main inputFile [outputFile [logFile [problemFile]]] [%coption[ %coption...]]\n",
    // INCHI✔️❌:         APP_DESCRIPTION, INCHI_SRC_REV,
    // INCHI✔️❌:         INCHI_BUILD_PLATFORM, INCHI_BUILD_COMPILER, INCHI_BUILD_DEBUG, __DATE__, __TIME__,
    // INCHI✔️❌:         RELEASE_IS_FINAL ? "" : " *** pre-release, for evaluation only ***",
    // INCHI✔️❌:         INCHI_OPTION_PREFX, INCHI_OPTION_PREFX);
    // INCHI✔️❌: #else
    // INCHI✔️❌:         "%s %-s\n%-s Build (%-s%-s) of %s %-s %-s\n\nUsage:\ninchi-1 inputFile [outputFile [logFile [problemFile]]] [%coption[ %coption...]]\n",
    // INCHI✔️❌:         APP_DESCRIPTION, INCHI_SRC_REV,
    // INCHI✔️❌:         INCHI_BUILD_PLATFORM, INCHI_BUILD_COMPILER, INCHI_BUILD_DEBUG, __DATE__, __TIME__,
    // INCHI✔️❌:         RELEASE_IS_FINAL ? "" : " *** pre-release, for evaluation only ***",
    // INCHI✔️❌:         INCHI_OPTION_PREFX, INCHI_OPTION_PREFX );
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( BUILD_WITH_AMI == 1 )
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f,
    // INCHI✔️❌:         "inchi-1 inputFiles... %cAMI [%coption[ %coption...]]\n",
    // INCHI✔️❌:         INCHI_OPTION_PREFX, INCHI_OPTION_PREFX, INCHI_OPTION_PREFX);
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "\nOptions:\n");
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "\nInput\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  STDIO       Use standard input/output streams\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  InpAux      Input structures in %s default aux. info format\n              (for use with STDIO)\n", INCHI_NAME);
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  SDF:DataHeader Read from the input SDfile the ID under this DataHeader\n");
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  START:n     Start at n-th input structure\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  END:n       Stop after n-th input structure\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  RECORD:n    Treat only n-th input structure\n");
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( BUILD_WITH_AMI == 1 )
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  AMI         Allow multiple input files (wildcards supported)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  AMIOutStd   Write output to stdout (in AMI mode)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  AMILogStd   Write log to stderr (in AMI mode)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  AMIPrbNone  Suppress creation of problem files (in AMI mode)\n");
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "Output\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  NoLabels    Omit structure number, DataHeader and ID from %s output\n", INCHI_NAME);
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  NoWarnings  Suppress all warning messages\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  AuxNone     Omit auxiliary information\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  SaveOpt     Save custom InChI creation options (non-standard InChI)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  Tabbed      Separate structure number, %s, and AuxInfo with tabs\n", INCHI_NAME);
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  MergeHash   Combine InChIKey with extra hash(es) if present\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  NoInChI     Do not print InChI string itself\n");
    // INCHI✔️❌: #ifndef TARGET_EXE_USING_API
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  OutErrInChI On fail, print empty InChI (default: nothing)\n");
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #if ( defined(_WIN32) && !defined(COMPILE_ANSI_ONLY) && !defined(TARGET_API_LIB) ) /* djb-rwth: check if this is working on GCC for Windows */
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  D           Display the structures\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  EQU         Display sets of identical components\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  Fnumber     Set display Font size in number of points\n");
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  OutputSDF   Convert %s created with default aux. info to SDfile\n", INCHI_NAME);
    // INCHI✔️❌: #if ( SDF_OUTPUT_DT == 1 )
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  SdfAtomsDT  Output Hydrogen Isotopes to SDfile as Atoms D and T\n");
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "Structure perception\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  SNon        Exclude stereo (default: include absolute stereo)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  NEWPSOFF    Both ends of wedge point to stereocenters (default: a narrow end)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  LooseTSACheck   Relax criteria of ambiguous drawing for in-ring tetrahedral stereo\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  DoNotAddH   All H are explicit (default: add H according to usual valences)\n");
    // INCHI✔️❌: #ifndef USE_STDINCHI_API
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "Stereo perception modifiers (non-standard InChI)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  SRel        Relative stereo\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  SRac        Racemic stereo\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  SUCF        Use Chiral Flag: On means Absolute stereo, Off - Relative\n");
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "Customizing InChI creation (non-standard InChI)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  SUU         Always include omitted unknown/undefined stereo\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  SLUUD       Make labels for unknown and undefined stereo different\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  RecMet      Include reconnected metals results\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  FixedH      Include Fixed H layer\n");
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  KET         Consider keto-enol tautomerism (experimental)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  15T         Consider 1,5-tautomerism (experimental)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  PT_06_00    Consider 1,3 heteroatom shift (experimental)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  PT_13_00    Consider keten-ynol exchange (experimental)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  PT_16_00    Consider nitroso-oxime tautomerism (experimental)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  PT_18_00    Consider cyanic/iso-cyanic acids (experimental)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  PT_22_00    Consider imine/imine tautomerism (experimental)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  PT_39_00    Consider nitrone/azoxy or Behrend rearrangement (experimental)\n");
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "Generation\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  Wnumber     Set time-out per structure in seconds; W0 means unlimited\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  WMnumber    Set time-out per structure in milliseconds (int); WM0 means unlimited\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  LargeMolecules Treat molecules up to 32766 atoms (experimental)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  WarnOnEmptyStructure Warn and produce empty %s for empty structure\n", INCHI_NAME);
    // INCHI✔️❌:     /*inchi_ios_print_nodisplay( f, "  MismatchIsError Treat problem/mismatch on inchi2struct conversion as error\n");*/
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  Polymers    Allow processing of polymers (experimental)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  Polymers105 Allow processing of polymers (experimental, legacy mode of v. 1.05)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  FoldCRU     Fold polymer CRU if inner repeats occur\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  NoFrameShift Disable polymer CRU frame shift\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  NoEdits     Disable polymer CRU frame shift and folding\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  NPZz        Allow non-polymer-related Zz atoms (pseudo element placeholders)\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  SAtZz       Allow stereo at atoms connected to Zz(default: disabled)\n");
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  Key         Generate InChIKey\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  XHash1      Generate hash extension (to 256 bits) for 1st block of InChIKey\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  XHash2      Generate hash extension (to 256 bits) for 2nd block of InChIKey\n");
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "Conversion\n");
    // INCHI✔️❌: #ifdef TARGET_EXE_USING_API
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  InChI2Struct Test mode: Mol/SDfile -> %s -> Structure -> (%s+AuxInfo)\n", INCHI_NAME, INCHI_NAME);
    // INCHI✔️❌: #else
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  InChI2Struct Convert InChI string(s) to structure(s) in InChI aux.info format\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  InChI2InChI  Convert  Convert %s string(s) into %s string(s)\n", INCHI_NAME, INCHI_NAME);
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: #if (BUILD_WITH_ENG_OPTIONS==1)
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "Engineering/hidden\n");
    // INCHI✔️❌: #ifdef TARGET_EXE_USING_API
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  InChI2InChI  Test mode: Mol/SDfile -> %s -> %s\n", INCHI_NAME, INCHI_NAME);
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     /*inchi_ios_print_nodisplay( f, "  Compress    Compressed output\n"); */
    // INCHI✔️❌:     /*inchi_ios_print_nodisplay( f, "    FULL        Standard set of options for Full Verbose Output\n");*/
    // INCHI✔️❌:     /*inchi_ios_print_nodisplay( f, "    MIN         Standard set of options for Minimal Concise Output\n");*/
    // INCHI✔️❌:
    // INCHI✔️❌: #if ALLOW_SUBSTRUCTURE_FILTERING==1
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  FilterSS    Select input SDF records using (hard-coded) substructure filter\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  InvFilterSS Invert match for (hard-coded) substructure filter\n");
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  Compress    Compressed output\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  MERGE       Use bMergeAllInputStructures\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  PGO         Use bSaveAllGoodStructsAsProblem\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  DCR         Use bDisplayCompositeResults\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  DSB         Use REQ_MODE_NO_ALT_SBONDS \n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  NOHDR       Use bNoStructLabels\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  NoVarH      Set bTgFlagVariableProtons=0\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  NOUUSB      Use REQ_MODE_SB_IGN_ALL_UU\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  NOUUSC      Use REQ_MODE_SC_IGN_ALL_UU\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  FixRad      Set bFixAdjacentRad\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  TestRenum   Generate InChI upon random atom renumbering\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  DoDRV       Set bUnderivatize=1\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  DoDrvReport Set bUnderivatize=3\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  DoR2C       Set bRing2Chain\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  DoneOnly    Set bIgnoreUnchanged\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  NoADP       Set bTgFlagHardAddRenProtons=0\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  MOVEPOS:0|1 Set bMovePositiveCharges\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  RSB:n       Set nMinDbRinSize\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  DISCONSALT:0|1     Set bDisconnectSalts\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  DISCONMETAL:0|1    Set bDisconnectCoord\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  DISCONMETALCHKVAL:0|1 Set bDisconnectCoordChkVal \n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  RECONMETAL:0|1     Set bReconnectCoord\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  MERGESALTTG:0|1    Set bMergeSaltTGroups\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  UNCHARGEDACIDS:0|1 Set bUnchargedAcidTaut \n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  ACIDTAUT:0|1|2     Set bAcidTautomerism\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  AUXINFO:0|1|2      Set AuxInfo print options\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  KeepBalanceP...  \n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  SDFID       ...\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  PLAINP      ....\n");
    // INCHI✔️❌:     inchi_ios_print_nodisplay(f, "  ANNPLAIN    ....\n");
    // INCHI✔️❌:
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: #endif
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: HelpCommandLineParms

    let Some(stream) = stream else {
        return Ok(());
    };
    macro_rules! print {
        ($format:expr $(, $argument:expr)* $(,)?) => {{
            nodisplay_call(heap, stream, stdout, $format, vec![$($argument),*])?;
        }};
    }

    print!(
        "%s %-s\n%-s Build (%-s%-s) of %s %-s %-s\n\nUsage:\ninchi-1 inputFile [outputFile [logFile [problemFile]]] [%coption[ %coption...]]\n",
        PrintArgument::Text("InChI version 1, Software 1.07.5 (API Library)"),
        PrintArgument::Text(""),
        PrintArgument::Text("Linux 64-bit"),
        PrintArgument::Bytes(c_text(build.compiler)),
        PrintArgument::Text(""),
        PrintArgument::Bytes(c_text(build.date)),
        PrintArgument::Bytes(c_text(build.time)),
        PrintArgument::Text(""),
        PrintArgument::Signed(i64::from(INCHI_OPTION_PREFX)),
        PrintArgument::Signed(i64::from(INCHI_OPTION_PREFX)),
    );

    for line in [
        "\nOptions:\n",
        "\nInput\n",
        "  STDIO       Use standard input/output streams\n",
    ] {
        print!(line);
    }
    print!(
        "  InpAux      Input structures in %s default aux. info format\n              (for use with STDIO)\n",
        PrintArgument::Text("InChI")
    );
    for line in [
        "  SDF:DataHeader Read from the input SDfile the ID under this DataHeader\n",
        "  START:n     Start at n-th input structure\n",
        "  END:n       Stop after n-th input structure\n",
        "  RECORD:n    Treat only n-th input structure\n",
        "Output\n",
    ] {
        print!(line);
    }
    print!(
        "  NoLabels    Omit structure number, DataHeader and ID from %s output\n",
        PrintArgument::Text("InChI")
    );
    for line in [
        "  NoWarnings  Suppress all warning messages\n",
        "  AuxNone     Omit auxiliary information\n",
        "  SaveOpt     Save custom InChI creation options (non-standard InChI)\n",
    ] {
        print!(line);
    }
    print!(
        "  Tabbed      Separate structure number, %s, and AuxInfo with tabs\n",
        PrintArgument::Text("InChI")
    );
    for line in [
        "  MergeHash   Combine InChIKey with extra hash(es) if present\n",
        "  NoInChI     Do not print InChI string itself\n",
        "  OutErrInChI On fail, print empty InChI (default: nothing)\n",
    ] {
        print!(line);
    }
    print!(
        "  OutputSDF   Convert %s created with default aux. info to SDfile\n",
        PrintArgument::Text("InChI")
    );
    for line in [
        "  SdfAtomsDT  Output Hydrogen Isotopes to SDfile as Atoms D and T\n",
        "Structure perception\n",
        "  SNon        Exclude stereo (default: include absolute stereo)\n",
        "  NEWPSOFF    Both ends of wedge point to stereocenters (default: a narrow end)\n",
        "  LooseTSACheck   Relax criteria of ambiguous drawing for in-ring tetrahedral stereo\n",
        "  DoNotAddH   All H are explicit (default: add H according to usual valences)\n",
        "Stereo perception modifiers (non-standard InChI)\n",
        "  SRel        Relative stereo\n",
        "  SRac        Racemic stereo\n",
        "  SUCF        Use Chiral Flag: On means Absolute stereo, Off - Relative\n",
        "Customizing InChI creation (non-standard InChI)\n",
        "  SUU         Always include omitted unknown/undefined stereo\n",
        "  SLUUD       Make labels for unknown and undefined stereo different\n",
        "  RecMet      Include reconnected metals results\n",
        "  FixedH      Include Fixed H layer\n",
        "  KET         Consider keto-enol tautomerism (experimental)\n",
        "  15T         Consider 1,5-tautomerism (experimental)\n",
        "  PT_06_00    Consider 1,3 heteroatom shift (experimental)\n",
        "  PT_13_00    Consider keten-ynol exchange (experimental)\n",
        "  PT_16_00    Consider nitroso-oxime tautomerism (experimental)\n",
        "  PT_18_00    Consider cyanic/iso-cyanic acids (experimental)\n",
        "  PT_22_00    Consider imine/imine tautomerism (experimental)\n",
        "  PT_39_00    Consider nitrone/azoxy or Behrend rearrangement (experimental)\n",
        "Generation\n",
        "  Wnumber     Set time-out per structure in seconds; W0 means unlimited\n",
        "  WMnumber    Set time-out per structure in milliseconds (int); WM0 means unlimited\n",
        "  LargeMolecules Treat molecules up to 32766 atoms (experimental)\n",
    ] {
        print!(line);
    }
    print!(
        "  WarnOnEmptyStructure Warn and produce empty %s for empty structure\n",
        PrintArgument::Text("InChI")
    );
    for line in [
        "  Polymers    Allow processing of polymers (experimental)\n",
        "  Polymers105 Allow processing of polymers (experimental, legacy mode of v. 1.05)\n",
        "  FoldCRU     Fold polymer CRU if inner repeats occur\n",
        "  NoFrameShift Disable polymer CRU frame shift\n",
        "  NoEdits     Disable polymer CRU frame shift and folding\n",
        "  NPZz        Allow non-polymer-related Zz atoms (pseudo element placeholders)\n",
        "  SAtZz       Allow stereo at atoms connected to Zz(default: disabled)\n",
        "  Key         Generate InChIKey\n",
        "  XHash1      Generate hash extension (to 256 bits) for 1st block of InChIKey\n",
        "  XHash2      Generate hash extension (to 256 bits) for 2nd block of InChIKey\n",
        "Conversion\n",
        "  InChI2Struct Convert InChI string(s) to structure(s) in InChI aux.info format\n",
    ] {
        print!(line);
    }
    print!(
        "  InChI2InChI  Convert  Convert %s string(s) into %s string(s)\n",
        PrintArgument::Text("InChI"),
        PrintArgument::Text("InChI")
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_common_option(
        heap: &mut SourceHeap,
        argument: &str,
        developer_options: i32,
        ip: &mut INPUT_PARMS,
        options: &mut CommonOptionsByParg,
    ) -> Result<i32, SourceHeapError> {
        let pointer = heap.allocate_model_storage(
            argument
                .bytes()
                .chain(std::iter::once(0))
                .map(|byte| byte as i8)
                .collect(),
        )?;
        let result =
            set_common_options_by_parg(heap, pointer.as_const(), developer_options, ip, options);
        heap.free(pointer)?;
        result
    }

    fn allocate_argv(heap: &mut SourceHeap, arguments: &[&str]) -> Vec<SourceConstPointer<i8>> {
        arguments
            .iter()
            .map(|argument| {
                heap.allocate_model_storage(
                    argument
                        .bytes()
                        .chain(std::iter::once(0))
                        .map(|byte| byte as i8)
                        .collect(),
                )
                .unwrap()
                .as_const()
            })
            .collect()
    }

    fn path_text(heap: &SourceHeap, pointer: SourceConstPointer<i8>) -> Option<String> {
        if pointer.is_null() {
            return None;
        }
        let bytes = heap.slice(pointer).unwrap();
        let nul = bytes.iter().position(|byte| *byte == 0).unwrap();
        Some(String::from_utf8(bytes[..nul].iter().map(|byte| *byte as u8).collect()).unwrap())
    }

    #[test]
    fn source_port__ichiparm__readcommandlineparms__line_602() {
        let default_taut_flags = (TG_FLAG_TEST_TAUT__ATOMS
            | TG_FLAG_DISCONNECT_SALTS
            | TG_FLAG_TEST_TAUT__SALTS
            | TG_FLAG_MOVE_POS_CHARGES
            | TG_FLAG_TEST_TAUT2_SALTS
            | TG_FLAG_MERGE_TAUT_SALTS
            | TG_FLAG_DISCONNECT_COORD
            | TG_FLAG_VARIABLE_PROTONS
            | TG_FLAG_HARD_ADD_REM_PROTONS
            | TG_FLAG_POINTED_EDGE_STEREO
            | TG_FLAG_PHOSPHINE_STEREO
            | TG_FLAG_ARSINE_STEREO
            | TG_FLAG_FIX_SP3_BUG
            | TG_FLAG_FIX_ISO_FIXEDH_BUG
            | TG_FLAG_FIX_TERM_H_CHRG_BUG) as INCHI_MODE;
        let default_mode = (REQ_MODE_TAUT
            | REQ_MODE_ISO
            | REQ_MODE_NON_ISO
            | REQ_MODE_STEREO
            | REQ_MODE_SB_IGN_ALL_UU
            | REQ_MODE_SC_IGN_ALL_UU
            | (MIN_SB_RING_SIZE << REQ_MODE_MIN_SB_RING_SHFT))
            as INCHI_MODE;

        let mut heap = SourceHeap::default();
        let argv = allocate_argv(&mut heap, &["inchi"]);
        let old_path = heap
            .allocate_model_storage(vec![b'o' as i8, b'l' as i8, b'd' as i8, 0])
            .unwrap();
        let mut ip = INPUT_PARMS {
            path: [
                old_path.as_const(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
            ],
            num_paths: 1,
            bNoWarnings: 9,
            ..INPUT_PARMS::default()
        };
        let mut display_time = 99_u64;
        assert_eq!(
            ReadCommandLineParms(
                &mut heap,
                argv.len() as i32,
                &argv,
                &mut ip,
                SourceMutPointer::null(),
                &mut display_time,
                1,
                None,
            ),
            Ok(0)
        );
        assert_eq!(display_time, 0);
        assert_eq!(ip.nInputType, tagInputType_INPUT_MOLFILE);
        assert_eq!(ip.nMode, default_mode);
        assert_eq!(ip.bTautFlags, default_taut_flags);
        assert_eq!(ip.bTautFlagsDone, 0);
        assert_eq!(
            ip.bINChIOutputOptions,
            (INCHI_OUT_EMBED_REC
                | INCHI_OUT_PLAIN_TEXT
                | INCHI_OUT_FIX_TRANSPOSITION_CHARGE_BUG
                | INCHI_OUT_STDINCHI) as i32
        );
        assert_eq!(ip.num_paths, 0);
        assert_eq!(ip.bFixNonUniformDraw, 1);

        let argv = allocate_argv(&mut heap, &["inchi", "sample"]);
        let mut paths = INPUT_PARMS::default();
        assert_eq!(
            ReadCommandLineParms(
                &mut heap,
                argv.len() as i32,
                &argv,
                &mut paths,
                SourceMutPointer::null(),
                &mut display_time,
                1,
                None,
            ),
            Ok(0)
        );
        assert_eq!(paths.num_paths, 4);
        assert_eq!(
            paths.path.map(|pointer| path_text(&heap, pointer)),
            [
                Some("sample".to_owned()),
                Some("sample.txt".to_owned()),
                Some("sample.log".to_owned()),
                Some("sample.prb".to_owned()),
            ]
        );

        let argv = allocate_argv(&mut heap, &["inchi", "NUL"]);
        let mut nul_paths = INPUT_PARMS::default();
        assert_eq!(
            ReadCommandLineParms(
                &mut heap,
                argv.len() as i32,
                &argv,
                &mut nul_paths,
                SourceMutPointer::null(),
                &mut display_time,
                1,
                None,
            ),
            Ok(0)
        );
        assert!(nul_paths.path[0].is_null());
        assert_eq!(
            path_text(&heap, nul_paths.path[1]).as_deref(),
            Some("NUL.txt")
        );
        assert_eq!(nul_paths.num_paths, 4);

        let sdf_value = heap
            .allocate_model_storage(vec![b'v' as i8, b'a' as i8, b'l' as i8, 0])
            .unwrap();
        let argv = allocate_argv(
            &mut heap,
            &[
                "inchi",
                "-SDF:  ID  ",
                "-WarnOnEmptyStructure",
                "-DoDrvReport",
                "-OUTPUTSDF",
                "-SdfAtomsDT",
            ],
        );
        let mut sdf = INPUT_PARMS::default();
        assert_eq!(
            ReadCommandLineParms(
                &mut heap,
                argv.len() as i32,
                &argv,
                &mut sdf,
                sdf_value,
                &mut display_time,
                0,
                None,
            ),
            Ok(0)
        );
        assert_eq!(sdf.nInputType, tagInputType_INPUT_SDFILE);
        assert_eq!(&sdf.szSdfDataHeader[..3], &[b'I' as i8, b'D' as i8, 0]);
        assert_eq!(
            source_c_string(&heap, sdf.pSdfLabel.as_const()).unwrap(),
            vec![b'I' as i8, b'D' as i8, 0]
        );
        assert_eq!(sdf.pSdfValue, sdf_value);
        assert_eq!(sdf.bAllowEmptyStructure, 1);
        assert_eq!(sdf.bUnderivatize, 3);
        assert_ne!(sdf.bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY as i32, 0);
        assert_ne!(
            sdf.bINChIOutputOptions & INCHI_OUT_SDFILE_ATOMS_DT as i32,
            0
        );
        assert_eq!(sdf.bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT as i32, 0);
        assert_eq!(sdf.bINChIOutputOptions & INCHI_OUT_STDINCHI as i32, 0);

        let argv = allocate_argv(
            &mut heap,
            &[
                "inchi",
                "-D",
                "-F12",
                "-EQU",
                "-WM250",
                "-W1.25tail",
                "-WM0",
                "-unknown",
            ],
        );
        let mut log = string_stream();
        let mut timeout = INPUT_PARMS::default();
        heap.set_source_errno(0);
        assert_eq!(
            ReadCommandLineParms(
                &mut heap,
                argv.len() as i32,
                &argv,
                &mut timeout,
                SourceMutPointer::null(),
                &mut display_time,
                1,
                Some(&mut log),
            ),
            Ok(0)
        );
        assert_eq!(timeout.msec_MaxTime, 1250);
        let messages = output(&heap, &log);
        assert!(messages.contains("may have been modified"));
        assert!(messages.contains("specified timeout value was ignored"));
        assert!(messages.contains("Unrecognized optionQ3: \"unknown\"."));

        let argv = allocate_argv(&mut heap, &["inchi", "-W0x1.8p1"]);
        let mut hexadecimal_timeout = INPUT_PARMS::default();
        heap.set_source_errno(0);
        assert_eq!(
            ReadCommandLineParms(
                &mut heap,
                argv.len() as i32,
                &argv,
                &mut hexadecimal_timeout,
                SourceMutPointer::null(),
                &mut display_time,
                1,
                None,
            ),
            Ok(0)
        );
        assert_eq!(hexadecimal_timeout.msec_MaxTime, 3000);

        let argv = allocate_argv(&mut heap, &["inchi", "-W2"]);
        let mut stale_errno_timeout = INPUT_PARMS::default();
        heap.set_source_errno(34);
        assert_eq!(
            ReadCommandLineParms(
                &mut heap,
                argv.len() as i32,
                &argv,
                &mut stale_errno_timeout,
                SourceMutPointer::null(),
                &mut display_time,
                1,
                None,
            ),
            Ok(0)
        );
        assert_eq!(stale_errno_timeout.msec_MaxTime, 0);

        let argv = allocate_argv(
            &mut heap,
            &["inchi", "-InChI2Struct", "-InChI2InChI", "-XHash1"],
        );
        let mut conversion_log = string_stream();
        let mut conversion = INPUT_PARMS::default();
        heap.set_source_errno(0);
        assert_eq!(
            ReadCommandLineParms(
                &mut heap,
                argv.len() as i32,
                &argv,
                &mut conversion,
                SourceMutPointer::null(),
                &mut display_time,
                1,
                Some(&mut conversion_log),
            ),
            Ok(0)
        );
        assert_eq!(conversion.nInputType, tagInputType_INPUT_INCHI);
        assert_ne!(
            conversion.bReadInChIOptions & READ_INCHI_OUTPUT_INCHI as i32,
            0
        );
        assert_eq!(
            conversion.bReadInChIOptions & READ_INCHI_TO_STRUCTURE as i32,
            0
        );
        assert_ne!(
            conversion.bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO as i32,
            0
        );
        assert!(output(&heap, &conversion_log).contains("InChIKey not requested"));

        for arguments in [
            ["inchi", "-Key", "-InChI2InChI"],
            ["inchi", "-Key", "-OUTPUTSDF"],
        ] {
            let argv = allocate_argv(&mut heap, &arguments);
            let mut early_log = string_stream();
            let mut early = INPUT_PARMS::default();
            assert_eq!(
                ReadCommandLineParms(
                    &mut heap,
                    argv.len() as i32,
                    &argv,
                    &mut early,
                    SourceMutPointer::null(),
                    &mut display_time,
                    1,
                    Some(&mut early_log),
                ),
                Ok(-1),
                "{arguments:?}"
            );
            assert!(output(&heap, &early_log).contains("Terminating:"));
        }

        let argv = allocate_argv(
            &mut heap,
            &[
                "inchi",
                "-PT_22_00",
                "-KET",
                "-RECMET",
                "-Key",
                "-XHash1",
                "-XHash2",
            ],
        );
        let mut common_log = string_stream();
        let mut common = INPUT_PARMS::default();
        assert_eq!(
            ReadCommandLineParms(
                &mut heap,
                argv.len() as i32,
                &argv,
                &mut common,
                SourceMutPointer::null(),
                &mut display_time,
                1,
                Some(&mut common_log),
            ),
            Ok(0)
        );
        assert_ne!(common.bTautFlags & TG_FLAG_PT_22_00 as u64, 0);
        assert_ne!(common.bTautFlags & TG_FLAG_KETO_ENOL_TAUT as u64, 0);
        assert_ne!(common.bTautFlags & TG_FLAG_RECONNECT_COORD as u64, 0);
        assert_eq!(
            common.bCalcInChIHash,
            tagInChIHashCalc_INCHIHASH_KEY_XTRA1_XTRA2 as i32
        );
        assert!(output(&heap, &common_log).contains("PT_22_00"));

        let argv = allocate_argv(&mut heap, &["inchi", "lost"]);
        let mut allocation_failure = INPUT_PARMS::default();
        heap.fail_after_allocations(0);
        assert_eq!(
            ReadCommandLineParms(
                &mut heap,
                argv.len() as i32,
                &argv,
                &mut allocation_failure,
                SourceMutPointer::null(),
                &mut display_time,
                1,
                None,
            ),
            Ok(0)
        );
        assert_eq!(allocation_failure.num_paths, 1);
        assert!(allocation_failure.path[0].is_null());

        let argv = allocate_argv(&mut heap, &["inchi"]);
        let mut bounds = INPUT_PARMS::default();
        assert_eq!(
            ReadCommandLineParms(
                &mut heap,
                2,
                &argv,
                &mut bounds,
                SourceMutPointer::null(),
                &mut display_time,
                1,
                None,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let unterminated = heap
            .allocate_model_storage(vec![b'i' as i8, b'n' as i8])
            .unwrap();
        assert_eq!(
            ReadCommandLineParms(
                &mut heap,
                2,
                &[argv[0], unterminated.as_const()],
                &mut bounds,
                SourceMutPointer::null(),
                &mut display_time,
                1,
                None,
            ),
            Err(SourceHeapError::MissingNulTerminator)
        );
    }

    #[test]
    fn source_port__ichiparm__set_common_options_by_parg__line_131() {
        let mut heap = SourceHeap::default();
        let mut ip = INPUT_PARMS {
            bFixNonUniformDraw: 7,
            ..INPUT_PARMS::default()
        };
        let mut options = CommonOptionsByParg {
            ver1_default_mode: u64::MAX,
            mode: -1,
            inchi_output_options: INCHI_OUT_SHORT_AUX_INFO as i32,
            inchi_output_options2: 0,
            std_format: 7,
            fix_sp3_bug: 7,
            fix_fb2: 7,
            add_phosphine_stereo: 7,
            add_arsine_stereo: 7,
            pointed_edge_stereo: 7,
            forced_chiral_flag: (FLAG_SET_INP_AT_CHIRAL | FLAG_SET_INP_AT_NONCHIRAL) as i32,
            fold_polymer_sru: 7,
            frame_shift_scheme: 7,
            ..CommonOptionsByParg::default()
        };

        assert_eq!(
            run_common_option(&mut heap, "iNpAuX", 0, &mut ip, &mut options),
            Ok(1)
        );
        assert_eq!(ip.nInputType, tagInputType_INPUT_INCHI_PLAIN);
        ip.nInputType = tagInputType_INPUT_MOLFILE;
        assert_eq!(
            run_common_option(&mut heap, "INPAUX", 0, &mut ip, &mut options),
            Ok(1)
        );
        assert_eq!(ip.nInputType, tagInputType_INPUT_MOLFILE);

        for argument in [
            "START:9223372036854775808",
            "END:-9223372036854775809",
            "RECORD:not-a-number",
        ] {
            assert_eq!(
                run_common_option(&mut heap, argument, 0, &mut ip, &mut options),
                Ok(1),
                "{argument}"
            );
        }
        assert_eq!(ip.first_struct_number, 0);
        assert_eq!(ip.last_struct_number, 0);
        assert_eq!(
            run_common_option(
                &mut heap,
                "START:9223372036854775808",
                0,
                &mut ip,
                &mut options
            ),
            Ok(1)
        );
        assert_eq!(ip.first_struct_number, i64::MAX);
        assert_eq!(
            run_common_option(
                &mut heap,
                "END:-9223372036854775809",
                0,
                &mut ip,
                &mut options
            ),
            Ok(1)
        );
        assert_eq!(ip.last_struct_number, i64::MIN);
        assert_eq!(
            run_common_option(&mut heap, "RECORD:+42tail", 0, &mut ip, &mut options),
            Ok(1)
        );
        assert_eq!((ip.first_struct_number, ip.last_struct_number), (42, 42));

        for argument in [
            "NOLABELS",
            "SAVEOPT",
            "AUXNONE",
            "MISMATCHISERROR",
            "OUTERRINCHI",
            "Key",
            "XHash1",
            "XHash2",
            "SNON",
            "NEWPSOFF",
            "DONOTADDH",
            "LooseTSACheck",
            "SREL",
            "SRAC",
            "SUCF",
            "ChiralFlagON",
            "ChiralFlagOFF",
            "SUU",
            "SLUUD",
            "FIXEDH",
            "RECMET",
            "KET",
            "15T",
            "LargeMolecules",
            "Polymers",
            "Polymers105",
            "NPZz",
            "NoWarnings",
            "MergeHash",
            "NoInChI",
            "HideInChI",
            "FoldCRU",
            "FoldSRU",
            "SATZZ",
        ] {
            assert_eq!(
                run_common_option(&mut heap, argument, 0, &mut ip, &mut options),
                Ok(1),
                "{argument}"
            );
        }
        assert_eq!(options.no_struct_labels, 1);
        assert_eq!(
            options.inchi_output_options,
            (INCHI_OUT_SAVEOPT | INCHI_OUT_NO_AUX_INFO) as i32
        );
        assert_eq!(
            options.inchi_output_options2,
            (INCHI_OUT_MISMATCH_AS_ERROR | INCHI_OUT_INCHI_GEN_ERROR) as i32
        );
        assert_eq!(
            (options.hash_key, options.hash_xtra1, options.hash_xtra2),
            (1, 1, 1)
        );
        assert_eq!(options.pointed_edge_stereo, 0);
        assert_eq!(options.do_not_add_h, 1);
        assert_eq!(
            options.mode & (REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO) as i32,
            0
        );
        assert_ne!(options.mode & REQ_MODE_CHIR_FLG_STEREO as i32, 0);
        assert_eq!(
            options.forced_chiral_flag
                & (FLAG_SET_INP_AT_CHIRAL | FLAG_SET_INP_AT_NONCHIRAL) as i32,
            FLAG_SET_INP_AT_NONCHIRAL as i32
        );
        assert_eq!(options.reconnect_coord, 1);
        assert_eq!((options.keto_enol_taut, options.taut_15_non_ring), (1, 1));
        assert_eq!((options.large_molecules, options.loose_tsa_check), (1, 1));
        assert_eq!(options.polymers, POLYMERS_LEGACY as i32);
        assert_eq!((options.np_zz, options.stereo_at_zz), (1, 1));
        assert_eq!(
            (options.no_warnings, options.merge_hash, options.hide_inchi),
            (1, 1, 1)
        );
        assert_eq!(options.fold_polymer_sru, 1);

        for (argument, field) in [
            ("PT_06_00", 6),
            ("PT_13_00", 13),
            ("PT_16_00", 16),
            ("PT_18_00", 18),
            ("PT_22_00", 22),
            ("PT_39_00", 39),
        ] {
            assert_eq!(
                run_common_option(&mut heap, argument, 0, &mut ip, &mut options),
                Ok(0),
                "source omits got=1 for PT_{field:02}_00"
            );
        }
        assert_eq!(
            (
                options.pt_06_00_taut,
                options.pt_13_00_taut,
                options.pt_16_00_taut,
                options.pt_18_00_taut,
                options.pt_22_00_taut,
                options.pt_39_00_taut,
            ),
            (1, 1, 1, 1, 1, 1)
        );

        for (argument, expected) in [
            ("FrameShift:None", tagFrameShifScheme_FSS_NONE as i32),
            (
                "FrameShift:Cyclize",
                tagFrameShifScheme_FSS_STARS_CYCLED as i32,
            ),
            (
                "FrameShift:MoveStars",
                tagFrameShifScheme_FSS_STARS_OPENED as i32,
            ),
            (
                "FrameShift:MoveBrackets",
                tagFrameShifScheme_FSS_STARS_ENDS_OPENED as i32,
            ),
            (
                "FrameShift:  MoveStars \t",
                tagFrameShifScheme_FSS_STARS_OPENED as i32,
            ),
            ("FrameShift:   ", tagFrameShifScheme_FSS_STARS_CYCLED as i32),
        ] {
            options.frame_shift_scheme = 99;
            assert_eq!(
                run_common_option(&mut heap, argument, 0, &mut ip, &mut options),
                Ok(1)
            );
            assert_eq!(options.frame_shift_scheme, expected, "{argument}");
        }
        options.frame_shift_scheme = 99;
        assert_eq!(
            run_common_option(&mut heap, "FrameShift:unknown", 0, &mut ip, &mut options),
            Ok(1)
        );
        assert_eq!(options.frame_shift_scheme, 99);
        assert_eq!(
            run_common_option(&mut heap, "NoFrameShift", 0, &mut ip, &mut options),
            Ok(1)
        );
        assert_eq!(
            options.frame_shift_scheme,
            tagFrameShifScheme_FSS_NONE as i32
        );
        options.fold_polymer_sru = 9;
        options.frame_shift_scheme = 9;
        assert_eq!(
            run_common_option(&mut heap, "NoEdits", 0, &mut ip, &mut options),
            Ok(1)
        );
        assert_eq!(options.fold_polymer_sru, 0);
        assert_eq!(
            options.frame_shift_scheme,
            tagFrameShifScheme_FSS_NONE as i32
        );

        assert_eq!(
            run_common_option(&mut heap, "PGO", 0, &mut ip, &mut options),
            Ok(0)
        );
        for argument in [
            "PGO",
            "FilterSS",
            "InvFilterSS",
            "FNUDOFF",
            "FixSp3bugOFF",
            "FBOFF",
            "FB2OFF",
            "SPXYZOFF",
            "SASXYZOFF",
            "Polymers105+",
        ] {
            assert_eq!(
                run_common_option(&mut heap, argument, 1, &mut ip, &mut options),
                Ok(1),
                "{argument}"
            );
        }
        assert_eq!(ip.bSaveAllGoodStructsAsProblem, 1);
        assert_eq!(ip.bFilterSS, -1);
        assert_eq!(ip.bFixNonUniformDraw, 0);
        assert_eq!(options.fix_sp3_bug, 0);
        assert_eq!(options.fix_fb2, 0);
        assert_eq!(options.add_phosphine_stereo, 0);
        assert_eq!(options.add_arsine_stereo, 0);
        assert_eq!(options.polymers, POLYMERS_LEGACY_PLUS as i32);
        assert_eq!(options.std_format, 0);

        let before_ip = ip.clone();
        let before_options = options.clone();
        assert_eq!(
            run_common_option(&mut heap, "not-an-option", 1, &mut ip, &mut options),
            Ok(0)
        );
        assert_eq!(ip, before_ip);
        assert_eq!(options, before_options);
    }

    fn string_stream() -> INCHI_IOSTREAM {
        INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        }
    }

    fn output(heap: &SourceHeap, stream: &INCHI_IOSTREAM) -> String {
        let bytes = heap.slice(stream.s.pStr.as_const()).unwrap();
        String::from_utf8(
            bytes[..stream.s.nUsedLength as usize]
                .iter()
                .map(|byte| *byte as u8)
                .collect(),
        )
        .unwrap()
    }

    #[test]
    fn source_port__ichiparm__printinputparms__line_2130() {
        let mut heap = SourceHeap::default();
        let mut stream = string_stream();
        let standard = INPUT_PARMS {
            nMode: REQ_MODE_DEFAULT as u64,
            bTautFlags: TG_FLAG_POINTED_EDGE_STEREO as u64,
            bINChIOutputOptions: (INCHI_OUT_STDINCHI | INCHI_OUT_PLAIN_TEXT) as i32,
            ..INPUT_PARMS::default()
        };
        assert_eq!(
            PrintInputParms(&mut heap, Some(&mut stream), &standard),
            Ok(0)
        );
        assert_eq!(
            output(&heap, &stream),
            "Generating standard InChI\nOutput format: Plain text\nFull Aux. info\nNo timeout\nUp to 1024 atoms per structure\n\n\n"
        );

        let path = heap.allocate_model_storage(vec![b'p' as i8, 0]).unwrap();
        let mut header = [0_i8; 65];
        header[..7].copy_from_slice(&[
            b'H' as i8, b'E' as i8, b'A' as i8, b'D' as i8, b'E' as i8, b'R' as i8, 0,
        ]);
        let mut rich = INPUT_PARMS {
            szSdfDataHeader: header,
            path: [
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                path.as_const(),
            ],
            first_struct_number: 3,
            last_struct_number: 9,
            nInputType: tagInputType_INPUT_MOLFILE,
            nMode: (REQ_MODE_BASIC
                | REQ_MODE_STEREO
                | REQ_MODE_NOEQ_STEREO
                | REQ_MODE_REDNDNT_STEREO
                | REQ_MODE_NO_ALT_SBONDS
                | REQ_MODE_RACEMIC_STEREO
                | REQ_MODE_SC_IGN_ALL_UU
                | REQ_MODE_DIFF_UU_STEREO
                | (3 << REQ_MODE_MIN_SB_RING_SHFT)) as u64,
            bAbcNumbers: 1,
            bINChIOutputOptions: (INCHI_OUT_PLAIN_TEXT
                | INCHI_OUT_TABBED_OUTPUT
                | INCHI_OUT_SAVEOPT
                | INCHI_OUT_NO_AUX_INFO) as i32,
            bINChIOutputOptions2: (INCHI_OUT_INCHI_GEN_ERROR | INCHI_OUT_MISMATCH_AS_ERROR) as i32,
            bCtPredecessors: 1,
            msec_MaxTime: -7,
            bSaveWarningStructsAsProblem: 1,
            bSaveAllGoodStructsAsProblem: 1,
            bGetMolfileNumber: 1,
            bDoNotAddH: 1,
            bCalcInChIHash: tagInChIHashCalc_INCHIHASH_KEY_XTRA1_XTRA2 as i32,
            bAllowEmptyStructure: 1,
            bLargeMolecules: 1,
            bLooseTSACheck: 1,
            bPolymers: POLYMERS_LEGACY_PLUS as i32,
            bFoldPolymerSRU: 1,
            bNPZz: 1,
            bStereoAtZz: 1,
            bMergeHash: 1,
            bNoWarnings: 1,
            bHideInChI: 1,
            bTautFlags: (TG_FLAG_VARIABLE_PROTONS
                | TG_FLAG_DISCONNECT_COORD
                | TG_FLAG_RECONNECT_COORD
                | TG_FLAG_KETO_ENOL_TAUT
                | TG_FLAG_1_5_TAUT
                | TG_FLAG_PT_22_00
                | TG_FLAG_PT_16_00
                | TG_FLAG_PT_06_00
                | TG_FLAG_PT_39_00
                | TG_FLAG_PT_13_00
                | TG_FLAG_PT_18_00) as u64,
            bRenumber: 1,
            bUnderivatize: 1,
            ..INPUT_PARMS::default()
        };
        let mut rich_stream = string_stream();
        assert_eq!(
            PrintInputParms(&mut heap, Some(&mut rich_stream), &rich),
            Ok(0)
        );
        let rich_output = output(&heap, &rich_stream);
        for expected in [
            "Both ends of wedge point to stereocenters",
            "Generate InChI upon random atom renumbering",
            "Convert input structure to derivative precursor",
            "Slow Redund. No AltBond Racemic Stereo ON",
            "Min. stereobond ring size: 3",
            "Account for PT_18_00 tautomerism",
            "Generating hash extension (two blocks)",
            "SDfile data header: \"HEADER\"",
            "Output format: Plain text, tabbed",
            "Representation: Compressed",
            "Skipping 2 structures",
            "Timeout per structure: -7 msec",
            "Up to 32766 atoms per structure",
            "v. 1.05 legacy mode with senior link placed at start",
            "Allowing stereo at atoms connected to Zz",
        ] {
            assert!(rich_output.contains(expected), "missing {expected:?}");
        }

        rich.bUnderivatize = 3;
        rich.nMode = 0;
        rich.bTautFlags = 0;
        rich.bINChIOutputOptions = (INCHI_OUT_SDFILE_ONLY | INCHI_OUT_SDFILE_ATOMS_DT) as i32;
        rich.bReadInChIOptions = READ_INCHI_TO_STRUCTURE as i32;
        rich.nInputType = tagInputType_INPUT_INCHI;
        rich.bPolymers = POLYMERS_LEGACY as i32;
        let mut alternate_stream = string_stream();
        assert_eq!(
            PrintInputParms(&mut heap, Some(&mut alternate_stream), &rich),
            Ok(0)
        );
        let alternate = output(&heap, &alternate_stream);
        for expected in [
            "Stereo OFF",
            "Outputs derivative information",
            "Undefined Mobile H mode",
            "Do not disconnect metals",
            "Do not account for keto-enol tautomerism",
            "Output format: SDfile only (without stereochemical info and atom coordinates)",
            "(v. 1.05 legacy mode)",
        ] {
            assert!(alternate.contains(expected), "missing {expected:?}");
        }
        assert_eq!(
            PrintInputParms(&mut heap, None, &INPUT_PARMS::default()),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichiparm__helpcommandlineparms__line_2700() {
        let mut heap = SourceHeap::default();
        let stdout = heap.allocate(vec![SourceFile::default()]).unwrap();
        let build = InchiBuildMetadata {
            compiler: "gcc oracle",
            date: "Jul 20 2026",
            time: "15:30:00",
        };
        assert_eq!(HelpCommandLineParms(&mut heap, None, stdout, build), Ok(()));
        assert!(heap.slice(stdout.as_const()).unwrap()[0].bytes.is_empty());

        let mut stream = string_stream();
        assert_eq!(
            HelpCommandLineParms(&mut heap, Some(&mut stream), stdout, build),
            Ok(())
        );
        let expected = [
            "InChI version 1, Software 1.07.5 (API Library) \nLinux 64-bit Build (gcc oracle) of Jul 20 2026 15:30:00 \n\nUsage:\ninchi-1 inputFile [outputFile [logFile [problemFile]]] [-option[ -option...]]\n",
            "\nOptions:\n",
            "\nInput\n",
            "  STDIO       Use standard input/output streams\n",
            "  InpAux      Input structures in InChI default aux. info format\n              (for use with STDIO)\n",
            "  SDF:DataHeader Read from the input SDfile the ID under this DataHeader\n",
            "  START:n     Start at n-th input structure\n",
            "  END:n       Stop after n-th input structure\n",
            "  RECORD:n    Treat only n-th input structure\n",
            "Output\n",
            "  NoLabels    Omit structure number, DataHeader and ID from InChI output\n",
            "  NoWarnings  Suppress all warning messages\n",
            "  AuxNone     Omit auxiliary information\n",
            "  SaveOpt     Save custom InChI creation options (non-standard InChI)\n",
            "  Tabbed      Separate structure number, InChI, and AuxInfo with tabs\n",
            "  MergeHash   Combine InChIKey with extra hash(es) if present\n",
            "  NoInChI     Do not print InChI string itself\n",
            "  OutErrInChI On fail, print empty InChI (default: nothing)\n",
            "  OutputSDF   Convert InChI created with default aux. info to SDfile\n",
            "  SdfAtomsDT  Output Hydrogen Isotopes to SDfile as Atoms D and T\n",
            "Structure perception\n",
            "  SNon        Exclude stereo (default: include absolute stereo)\n",
            "  NEWPSOFF    Both ends of wedge point to stereocenters (default: a narrow end)\n",
            "  LooseTSACheck   Relax criteria of ambiguous drawing for in-ring tetrahedral stereo\n",
            "  DoNotAddH   All H are explicit (default: add H according to usual valences)\n",
            "Stereo perception modifiers (non-standard InChI)\n",
            "  SRel        Relative stereo\n",
            "  SRac        Racemic stereo\n",
            "  SUCF        Use Chiral Flag: On means Absolute stereo, Off - Relative\n",
            "Customizing InChI creation (non-standard InChI)\n",
            "  SUU         Always include omitted unknown/undefined stereo\n",
            "  SLUUD       Make labels for unknown and undefined stereo different\n",
            "  RecMet      Include reconnected metals results\n",
            "  FixedH      Include Fixed H layer\n",
            "  KET         Consider keto-enol tautomerism (experimental)\n",
            "  15T         Consider 1,5-tautomerism (experimental)\n",
            "  PT_06_00    Consider 1,3 heteroatom shift (experimental)\n",
            "  PT_13_00    Consider keten-ynol exchange (experimental)\n",
            "  PT_16_00    Consider nitroso-oxime tautomerism (experimental)\n",
            "  PT_18_00    Consider cyanic/iso-cyanic acids (experimental)\n",
            "  PT_22_00    Consider imine/imine tautomerism (experimental)\n",
            "  PT_39_00    Consider nitrone/azoxy or Behrend rearrangement (experimental)\n",
            "Generation\n",
            "  Wnumber     Set time-out per structure in seconds; W0 means unlimited\n",
            "  WMnumber    Set time-out per structure in milliseconds (int); WM0 means unlimited\n",
            "  LargeMolecules Treat molecules up to 32766 atoms (experimental)\n",
            "  WarnOnEmptyStructure Warn and produce empty InChI for empty structure\n",
            "  Polymers    Allow processing of polymers (experimental)\n",
            "  Polymers105 Allow processing of polymers (experimental, legacy mode of v. 1.05)\n",
            "  FoldCRU     Fold polymer CRU if inner repeats occur\n",
            "  NoFrameShift Disable polymer CRU frame shift\n",
            "  NoEdits     Disable polymer CRU frame shift and folding\n",
            "  NPZz        Allow non-polymer-related Zz atoms (pseudo element placeholders)\n",
            "  SAtZz       Allow stereo at atoms connected to Zz(default: disabled)\n",
            "  Key         Generate InChIKey\n",
            "  XHash1      Generate hash extension (to 256 bits) for 1st block of InChIKey\n",
            "  XHash2      Generate hash extension (to 256 bits) for 2nd block of InChIKey\n",
            "Conversion\n",
            "  InChI2Struct Convert InChI string(s) to structure(s) in InChI aux.info format\n",
            "  InChI2InChI  Convert  Convert InChI string(s) into InChI string(s)\n",
        ]
        .concat();
        assert_eq!(output(&heap, &stream), expected);
        assert!(!expected.contains("AMIOutStd"));
        assert!(!expected.contains("Engineering/hidden"));
        assert!(!expected.contains("Display the structures"));
    }
}
