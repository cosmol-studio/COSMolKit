use crate::source::base::ikey_base26::{
    base26_dublet_for_bits_28_to_36, base26_dublet_for_bits_56_to_64, base26_triplet_1,
    base26_triplet_2, base26_triplet_3, base26_triplet_4, get_xtra_hash_major_hex,
    get_xtra_hash_minor_hex,
};
use crate::source::base::sha2::sha2_csum;
use crate::source::base::util::{extract_inchi_substring, inchi_calloc, inchi_free};
use crate::source_types::{
    INCHIKEY_EMPTY_INPUT, INCHIKEY_INVALID_INCHI, INCHIKEY_INVALID_INCHI_PREFIX,
    INCHIKEY_INVALID_STD_INCHI, INCHIKEY_NOT_ENOUGH_MEMORY, INCHIKEY_OK, SourceConstPointer,
    SourceHeap, SourceHeapError, SourceMutPointer,
};

fn source_strtol_decimal(bytes: &[u8], start: usize) -> i32 {
    let mut index = start;
    let negative = match bytes.get(index).copied() {
        Some(b'-') => {
            index += 1;
            true
        }
        Some(b'+') => {
            index += 1;
            false
        }
        _ => false,
    };
    let mut saw_digit = false;
    let mut value = 0_u64;
    let limit = if negative {
        (i64::MAX as u64) + 1
    } else {
        i64::MAX as u64
    };
    while let Some(digit) = bytes
        .get(index)
        .copied()
        .filter(u8::is_ascii_digit)
        .map(|byte| u64::from(byte - b'0'))
    {
        saw_digit = true;
        value = value
            .checked_mul(10)
            .and_then(|current| current.checked_add(digit))
            .unwrap_or(limit)
            .min(limit);
        index += 1;
    }
    if !saw_digit {
        return 0;
    }
    let signed = if negative {
        if value == (i64::MAX as u64) + 1 {
            i64::MIN
        } else {
            -(value as i64)
        }
    } else {
        value as i64
    };
    signed as i32
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn GetINCHIKeyFromINCHI(
    heap: &mut SourceHeap,
    szINCHISource: SourceConstPointer<i8>,
    xtra1: i32,
    xtra2: i32,
    szINCHIKey: SourceMutPointer<i8>,
    szXtra1: SourceMutPointer<i8>,
    szXtra2: SourceMutPointer<i8>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c:113 GetINCHIKeyFromINCHI
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
        EXPIMP_TEMPLATE INCHI_API
        int INCHI_DECL GetINCHIKeyFromINCHI( const char* szINCHISource,
                                             const int xtra1,
                                             const int xtra2,
                                             char* szINCHIKey,
                                             char* szXtra1,
                                             char* szXtra2 )
        {
            int ret = INCHIKEY_OK;
            int ret1 = INCHIKEY_OK;
            int cn;
            size_t slen, i, j, jproto = 0, ncp, pos_slash1 = 0;
            char *str = NULL, *smajor = NULL, *sminor = NULL,
                *sproto = NULL, *stmp = NULL, tmp[MINOUTLENGTH];
            unsigned char
                digest_major[32], digest_minor[32];

            char flagstd = 'S', /* standard key */
                flagnonstd = 'N', /* non-standard key */
                flagexptl = 'B', /* experimental ('beta') key */
                flagver = 'A', /* InChI v. 1 */
                flagproto = 'N'; /* no [de]protonization , by default */
            int  nprotons;
            /*
            Protonization encoding:
            N 0
            O +1 P +2 Q +3 R +4 S +5 T +6 U +7 V +8 W +9 X +10 Y +11 Z +12
            M -1 L-2 K -3 J -4 I -5 H -6 G -7 F -8 E -9 D -10 C -11 B -12
            A < -12 or > +12
            */
            static const char *pplus = "OPQRSTUVWXYZ";
            static const char *pminus = "MLKJIHGFEDCB";
            int is_stdinchi = 0;    /* 0 -  non-standard,
                                       1    standard
                                      -1    experimental ('beta') */



            if (NULL != szXtra1) /* Software version 1.06 added check to fix bug with NULL szXtra, thanks to WDI */
            {
                szXtra1[0] = '\0';
            }
            if (NULL != szXtra2)
            {
                szXtra2[0] = '\0';
            }

            /* Check if input is a valid InChI string */

            /* .. non-empty */
            if (szINCHISource == NULL)
            {
                return INCHIKEY_EMPTY_INPUT;
            }

            slen = strlen( szINCHISource );

            /* .. has valid prefix */
            if (slen < LEN_INCHI_STRING_PREFIX + 3)
            {
                return INCHIKEY_INVALID_INCHI_PREFIX;
            }
            if (memcmp( szINCHISource, INCHI_STRING_PREFIX, LEN_INCHI_STRING_PREFIX ))
            {
                return INCHIKEY_INVALID_INCHI_PREFIX;
            }

            /* .. has InChI version 1 */
            /* if (!isdigit(szINCHISource[LEN_INCHI_STRING_PREFIX]) )  */
            if (szINCHISource[LEN_INCHI_STRING_PREFIX] != '1')
                return INCHIKEY_INVALID_INCHI_PREFIX;

            /* .. optionally has a 'standard' flag character */
            pos_slash1 = LEN_INCHI_STRING_PREFIX + 1;
            if (szINCHISource[pos_slash1] == 'S')
            {
                /* Standard InChI ==> standard InChIKey */
                is_stdinchi = 1;
                pos_slash1++;
            }
            else if (szINCHISource[pos_slash1] == 'B')
            {
                /* v. 1.05 Experimental ('beta') InChI ==> corresponding InChIKey */
                is_stdinchi = -1;
                pos_slash1++;
            }

            /* .. has trailing slash in the right place */
            if (szINCHISource[pos_slash1] != '/')
            {
                return INCHIKEY_INVALID_INCHI_PREFIX;
            }

            /* .. the rest of source string contains at least one a..Z0.9 or slash */
            if (!isalnum( szINCHISource[pos_slash1 + 1] ))
            {
                int valid_empty_inchi = szINCHISource[pos_slash1 + 1] == '/';
                int allowed_and_valid_err_inchi = 0;
                if (!valid_empty_inchi)
                {
                    allowed_and_valid_err_inchi = szINCHISource[pos_slash1 + 1] == '?';
                    if (!allowed_and_valid_err_inchi)
                    {
                        return INCHIKEY_INVALID_INCHI;
                    }
                }
            }

            /* Ok. Will use a local copy of the source. */

            extract_inchi_substring( &str, szINCHISource, slen );
            if (NULL == str)
            {
                ret = INCHIKEY_NOT_ENOUGH_MEMORY;
                goto fin;
            }
            slen = strlen( str );


            /* Make buffers. */
            smajor = (char*) inchi_calloc( slen + 1, sizeof( char ) );
            if (NULL == smajor)
            {
                ret = INCHIKEY_NOT_ENOUGH_MEMORY; goto fin;
            }
            sminor = (char*) inchi_calloc( 2 * slen + 2, sizeof( char ) ); /* we may double the length ... */
            if (NULL == sminor)
            {
                ret = INCHIKEY_NOT_ENOUGH_MEMORY; goto fin;
            }
            stmp = (char*) inchi_calloc( slen + 1, sizeof( char ) );
            if (NULL == stmp)
            {
                ret = INCHIKEY_NOT_ENOUGH_MEMORY; goto fin;
            }
            sproto = (char*) inchi_calloc( slen + 1, sizeof( char ) );
            if (NULL == sproto)
            {
                ret = INCHIKEY_NOT_ENOUGH_MEMORY; goto fin;
            }

            szINCHIKey[0] = '\0';

            /* Extract the major block */
            smajor[0] = '\0';
            for (j = pos_slash1 + 1; j < slen - 1; j++)
            {
                if (str[j] == '/')
                {
                    cn = str[j + 1];
                    switch (cn)
                    {
                        /* anything allowed from a major part */
                        case 'c':
                        case 'h':
                        case 'q':   continue;

                        /* "/p"; protons now go to to special string, not to minor hash */
                        case 'p':
                            jproto = j;
                            continue;

                /* "/f",  "/r" : may not occur in stdInChI */
                        case 'f':
                        case 'r':
                            if (is_stdinchi==1)
                            {
                                ret = INCHIKEY_INVALID_STD_INCHI;
                                goto fin;
                            }
                            break;

                /* anything allowed from a minor part */
                        default:
                            break;
                    }
                    break;
                }
            }
            j++;
            if (j == slen)
            {
                j++;
            }
            else
            {
                j--;
            }

            if (jproto)
            {
                ncp = jproto - pos_slash1 - 1;
            }
            else
            {
                ncp = j - pos_slash1 - 1;
            }


            /* Trim 'InChI=1[S]/' */
            memcpy(smajor, str + pos_slash1 + 1, ncp * sizeof(str[0]));
            smajor[ncp] = '\0';


            /* Treat protonization */
            if (jproto)
            {
                /* 2009-01-07 fix bug/typo: assigned incorrect length to the protonation segment of
                   source string ( was sproto[ncp]='\0'; should be sproto[lenproto]='\0'; )  */
                int lenproto = j - (int) jproto;
                if (lenproto < 3)
                {
                    /* empty "/p", should not occur */
                    ret = INCHIKEY_INVALID_INCHI;
                    goto fin;
                }

                memcpy(sproto, str + pos_slash1 + ncp + 1, lenproto * sizeof(str[0]));   
                sproto[lenproto] = '\0';

                nprotons = strtol( sproto + 2, NULL, 10 );

                if (nprotons > 0)
                {
                    if (nprotons > 12)
                    {
                        flagproto = 'A';
                    }
                    else
                    {
                        flagproto = pplus[nprotons - 1];
                    }
                }
                else if (nprotons < 0)
                {
                    if (nprotons < -12)
                    {
                        flagproto = 'A';
                    }
                    else
                    {
                        flagproto = pminus[-nprotons - 1];
                    }
                }
                else
                {
                    /* should never occur */
                    ret = INCHIKEY_INVALID_STD_INCHI;
                    goto fin;
                }
            }

            /* Extract the minor block. */

            if (j != slen + 1)    /* check that something exists at right.*/
            {
                ncp = slen - j;
                memcpy(sminor, str + j, (ncp) * sizeof(str[0]));
                sminor[ncp] = '\0';
            }
            else
            {
                sminor[0] = '\0';
            }


#if INCHIKEY_DEBUG
            ITRACE_( "Source:  {%-s}\n", str );
            ITRACE_( "SMajor:  {%-s}\n", smajor );
            ITRACE_( "SMinor:  {%-s}\n", sminor );
            ITRACE_( "SProto:  {%-s}\n", sproto );
#endif

            /* Compute and compose the InChIKey string. */

            /* Major hash sub-string. */
            for (i = 0; i < 32; i++)
            {
                digest_major[i] = 0;
            }

            sha2_csum( (unsigned char *) smajor, (int) strlen( smajor ), digest_major );

            sprintf(tmp, "%-.3s%-.3s%-.3s%-.3s%-.2s",
                base26_triplet_1(digest_major), base26_triplet_2(digest_major),
                base26_triplet_3(digest_major), base26_triplet_4(digest_major),
                base26_dublet_for_bits_56_to_64(digest_major));
            strcat(szINCHIKey, tmp);
#if (INCHIKEY_DEBUG>1)
            fprint_digest( stderr, "Major hash, full SHA-256", digest_major );
#endif


            /* Minor hash sub-string. */
            for (i = 0; i < 32; i++)
            {
                digest_minor[i] = 0;
            }
            slen = strlen( sminor );
            if (( slen > 0 ) && ( slen < 255 ))
            {
                strcpy(stmp, sminor);
                strcpy(sminor + slen, stmp);
            }

            sha2_csum( (unsigned char *) sminor, (int) strlen( sminor ), digest_minor );

#if (INCHIKEY_DEBUG>1)
            fprint_digest( stderr, "Minor hash, full SHA-256", digest_minor );
#endif

            strcat(szINCHIKey, "-");
            sprintf(tmp, "%-.3s%-.3s%-.2s",
                base26_triplet_1(digest_minor),
                base26_triplet_2(digest_minor),
                base26_dublet_for_bits_28_to_36(digest_minor));
            strcat(szINCHIKey, tmp);
            /* Append a standard/non-standard flag */
            slen = strlen( szINCHIKey );
            if (is_stdinchi == 1)
            {
                szINCHIKey[slen] = flagstd;
            }
            else if (is_stdinchi == -1)
            {
                szINCHIKey[slen] = flagexptl;
            }
            else
            {
                szINCHIKey[slen] = flagnonstd;
            }

            /* Append InChI v.1 flag */
            szINCHIKey[slen + 1] = flagver;

            /* Append dash  */
            szINCHIKey[slen + 2] = '-';

            /* Append protonization flag */
            szINCHIKey[slen + 3] = flagproto;
            szINCHIKey[slen + 4] = '\0';


#if INCHIKEY_DEBUG
            ITRACE_( "szINCHIKey:  {%-s}\n", szINCHIKey );
#endif

            /* Hash extensions */
            if (xtra1 && szXtra1)
            {
                get_xtra_hash_major_hex( digest_major, szXtra1 );
#if INCHIKEY_DEBUG
                fprintf( stderr, "XHash1=%-s\n", szXtra1 );
                fprintf( stderr, "j=%-d\n", j );
#endif
            }
            if (xtra2 && szXtra2)
            {
                get_xtra_hash_minor_hex( digest_minor, szXtra2 );
#if INCHIKEY_DEBUG
                fprintf( stderr, "XHash2=%-s\n", szXtra2 );
                fprintf( stderr, "j=%-d\n", j );
#endif
            }


        fin:
            /* djb-rwth: fixing GH issue #87 (thanks to Ricardo Rodriguez) and oss-fuzz issue #66746 */
            if (NULL != str)
            {
                inchi_free( str );
            }
            if (NULL != smajor)
            {
                inchi_free( smajor );
            }
            if (NULL != sminor)
            {
                inchi_free( sminor );
                sminor = NULL;
            }
            if (NULL != stmp)
            {
                inchi_free( stmp );
            }
            if (NULL != sproto)
            {
                inchi_free( sproto );
            }
            if (( ret == INCHIKEY_OK ) && ( ret1 != INCHIKEY_OK ))
            {
                ret = ret1;
            }

            return ret;
        }
            */
    // END INCHI C FUNCTION: GetINCHIKeyFromINCHI
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: GetINCHIKeyFromINCHI
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; INCHIKEY_DEBUG=0.
    // INCHI✔️❌: inchi_calloc and inchi_free are active allocation macros; libc strlen,
    // INCHI✔️❌: memcmp, memcpy, isalnum, strtol, sprintf, strcat, and strcpy behavior is active.
    // INCHI✔️❌: SourceHeap checked access and byte-buffer modeling add known overhead.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: GetINCHIKeyFromINCHI

    if !szXtra1.is_null() {
        heap.slice_mut(szXtra1)?[0] = 0;
    }
    if !szXtra2.is_null() {
        heap.slice_mut(szXtra2)?[0] = 0;
    }
    if szINCHISource.is_null() {
        return Ok(INCHIKEY_EMPTY_INPUT as i32);
    }

    let source_length = heap
        .slice(szINCHISource)?
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let source = heap.slice(szINCHISource)?;
    if source_length < 9
        || source[..6]
            != [
                b'I' as i8, b'n' as i8, b'C' as i8, b'h' as i8, b'I' as i8, b'=' as i8,
            ]
        || source[6] != b'1' as i8
    {
        return Ok(INCHIKEY_INVALID_INCHI_PREFIX as i32);
    }

    let mut pos_slash1 = 7_usize;
    let is_stdinchi = if source[pos_slash1] == b'S' as i8 {
        pos_slash1 += 1;
        1_i32
    } else if source[pos_slash1] == b'B' as i8 {
        pos_slash1 += 1;
        -1_i32
    } else {
        0_i32
    };
    if source[pos_slash1] != b'/' as i8 {
        return Ok(INCHIKEY_INVALID_INCHI_PREFIX as i32);
    }
    let first_body_byte = source[pos_slash1 + 1] as u8;
    if !first_body_byte.is_ascii_alphanumeric()
        && first_body_byte != b'/'
        && first_body_byte != b'?'
    {
        return Ok(INCHIKEY_INVALID_INCHI as i32);
    }

    let mut str_pointer = SourceMutPointer::null();
    let mut smajor = SourceMutPointer::null();
    let mut sminor = SourceMutPointer::null();
    let mut stmp = SourceMutPointer::null();
    let mut sproto = SourceMutPointer::null();
    let mut ret = INCHIKEY_OK as i32;

    match extract_inchi_substring(heap, &mut str_pointer, szINCHISource, source_length as u64) {
        Ok(()) => {}
        Err(SourceHeapError::AllocationFailed) => ret = INCHIKEY_NOT_ENOUGH_MEMORY as i32,
        Err(error) => return Err(error),
    }
    if ret == INCHIKEY_OK as i32 && str_pointer.is_null() {
        ret = INCHIKEY_NOT_ENOUGH_MEMORY as i32;
    }
    let mut slen = 0_usize;
    if ret == INCHIKEY_OK as i32 {
        slen = heap
            .slice(str_pointer.as_const())?
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        let allocate = |result: Result<SourceMutPointer<u8>, SourceHeapError>| match result {
            Ok(pointer) => Ok(pointer),
            Err(SourceHeapError::AllocationFailed) => Ok(SourceMutPointer::null()),
            Err(error) => Err(error),
        };
        smajor = allocate(inchi_calloc(heap, (slen + 1) as u64, 1))?;
        if smajor.is_null() {
            ret = INCHIKEY_NOT_ENOUGH_MEMORY as i32;
        }
        if ret == INCHIKEY_OK as i32 {
            sminor = allocate(inchi_calloc(heap, (2 * slen + 2) as u64, 1))?;
            if sminor.is_null() {
                ret = INCHIKEY_NOT_ENOUGH_MEMORY as i32;
            }
        }
        if ret == INCHIKEY_OK as i32 {
            stmp = allocate(inchi_calloc(heap, (slen + 1) as u64, 1))?;
            if stmp.is_null() {
                ret = INCHIKEY_NOT_ENOUGH_MEMORY as i32;
            }
        }
        if ret == INCHIKEY_OK as i32 {
            sproto = allocate(inchi_calloc(heap, (slen + 1) as u64, 1))?;
            if sproto.is_null() {
                ret = INCHIKEY_NOT_ENOUGH_MEMORY as i32;
            }
        }
    }

    if ret == INCHIKEY_OK as i32 {
        heap.slice_mut(szINCHIKey)?[0] = 0;
        let mut jproto = 0_usize;
        let mut j = pos_slash1 + 1;
        while j < slen - 1 {
            if heap.slice(str_pointer.as_const())?[j] == b'/' as i8 {
                match heap.slice(str_pointer.as_const())?[j + 1] as u8 {
                    b'c' | b'h' | b'q' => {
                        j += 1;
                        continue;
                    }
                    b'p' => {
                        jproto = j;
                        j += 1;
                        continue;
                    }
                    b'f' | b'r' if is_stdinchi == 1 => {
                        ret = INCHIKEY_INVALID_STD_INCHI as i32;
                        break;
                    }
                    _ => break,
                }
            }
            j += 1;
        }
        if ret == INCHIKEY_OK as i32 {
            j += 1;
            if j == slen {
                j += 1;
            } else {
                j -= 1;
            }
            let major_length = if jproto != 0 {
                jproto - pos_slash1 - 1
            } else {
                j - pos_slash1 - 1
            };
            {
                for index in 0..major_length {
                    let byte = heap.slice(str_pointer.as_const())?[pos_slash1 + 1 + index] as u8;
                    heap.slice_mut(smajor)?[index] = byte;
                }
                heap.slice_mut(smajor)?[major_length] = 0;
            }

            let mut flagproto = b'N';
            if jproto != 0 {
                let proto_length = j - jproto;
                if proto_length < 3 {
                    ret = INCHIKEY_INVALID_INCHI as i32;
                } else {
                    {
                        for index in 0..proto_length {
                            let byte = heap.slice(str_pointer.as_const())?[jproto + index] as u8;
                            heap.slice_mut(sproto)?[index] = byte;
                        }
                        heap.slice_mut(sproto)?[proto_length] = 0;
                    }
                    let nprotons =
                        source_strtol_decimal(&heap.slice(sproto.as_const())?[..proto_length], 2);
                    flagproto = match nprotons {
                        1..=12 => b"OPQRSTUVWXYZ"[(nprotons - 1) as usize],
                        -12..=-1 => b"MLKJIHGFEDCB"[(-nprotons - 1) as usize],
                        value if value > 12 || value < -12 => b'A',
                        _ => {
                            ret = INCHIKEY_INVALID_STD_INCHI as i32;
                            b'N'
                        }
                    };
                }
            }

            if ret == INCHIKEY_OK as i32 {
                let minor_length = if j != slen + 1 { slen - j } else { 0 };
                {
                    for index in 0..minor_length {
                        let byte = heap.slice(str_pointer.as_const())?[j + index] as u8;
                        heap.slice_mut(sminor)?[index] = byte;
                    }
                    heap.slice_mut(sminor)?[minor_length] = 0;
                }

                let major_length = heap
                    .slice(smajor.as_const())?
                    .iter()
                    .position(|byte| *byte == 0)
                    .ok_or(SourceHeapError::MissingNulTerminator)?;
                let mut digest_major = [0_u8; 32];
                sha2_csum(
                    &heap.slice(smajor.as_const())?[..major_length],
                    major_length as i32,
                    &mut digest_major,
                );

                let mut key = [0_u8; 28];
                let mut key_index = 0_usize;
                for encoded in [
                    base26_triplet_1(&digest_major).as_slice(),
                    base26_triplet_2(&digest_major).as_slice(),
                    base26_triplet_3(&digest_major).as_slice(),
                    base26_triplet_4(&digest_major).as_slice(),
                    base26_dublet_for_bits_56_to_64(&digest_major).as_slice(),
                ] {
                    key[key_index..key_index + encoded.len()].copy_from_slice(encoded);
                    key_index += encoded.len();
                }

                let minor_length = heap
                    .slice(sminor.as_const())?
                    .iter()
                    .position(|byte| *byte == 0)
                    .ok_or(SourceHeapError::MissingNulTerminator)?;
                if minor_length > 0 && minor_length < 255 {
                    for index in 0..=minor_length {
                        let byte = heap.slice(sminor.as_const())?[index];
                        heap.slice_mut(stmp)?[index] = byte;
                    }
                    for index in 0..=minor_length {
                        let byte = heap.slice(stmp.as_const())?[index];
                        heap.slice_mut(sminor)?[minor_length + index] = byte;
                    }
                }
                let hashed_minor_length = heap
                    .slice(sminor.as_const())?
                    .iter()
                    .position(|byte| *byte == 0)
                    .ok_or(SourceHeapError::MissingNulTerminator)?;
                let mut digest_minor = [0_u8; 32];
                sha2_csum(
                    &heap.slice(sminor.as_const())?[..hashed_minor_length],
                    hashed_minor_length as i32,
                    &mut digest_minor,
                );

                key[key_index] = b'-';
                key_index += 1;
                for encoded in [
                    base26_triplet_1(&digest_minor).as_slice(),
                    base26_triplet_2(&digest_minor).as_slice(),
                    base26_dublet_for_bits_28_to_36(&digest_minor).as_slice(),
                ] {
                    key[key_index..key_index + encoded.len()].copy_from_slice(encoded);
                    key_index += encoded.len();
                }
                key[key_index] = match is_stdinchi {
                    1 => b'S',
                    -1 => b'B',
                    _ => b'N',
                };
                key[key_index + 1] = b'A';
                key[key_index + 2] = b'-';
                key[key_index + 3] = flagproto;
                key[key_index + 4] = 0;
                for (destination, source) in heap.slice_mut(szINCHIKey)?.iter_mut().zip(key.iter())
                {
                    *destination = *source as i8;
                }

                if xtra1 != 0 && !szXtra1.is_null() {
                    let mut xtra = [0_u8; 49];
                    get_xtra_hash_major_hex(&digest_major, &mut xtra);
                    for (destination, source) in
                        heap.slice_mut(szXtra1)?.iter_mut().zip(xtra.iter())
                    {
                        *destination = *source as i8;
                    }
                }
                if xtra2 != 0 && !szXtra2.is_null() {
                    let mut xtra = [0_u8; 57];
                    get_xtra_hash_minor_hex(&digest_minor, &mut xtra);
                    for (destination, source) in
                        heap.slice_mut(szXtra2)?.iter_mut().zip(xtra.iter())
                    {
                        *destination = *source as i8;
                    }
                }
            }
        }
    }

    if !str_pointer.is_null() {
        inchi_free(heap, str_pointer)?;
    }
    if !smajor.is_null() {
        inchi_free(heap, smajor)?;
    }
    if !sminor.is_null() {
        inchi_free(heap, sminor)?;
    }
    if !stmp.is_null() {
        inchi_free(heap, stmp)?;
    }
    if !sproto.is_null() {
        inchi_free(heap, sproto)?;
    }
    Ok(ret)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    fn allocate_bytes(heap: &mut SourceHeap, bytes: &[u8]) -> SourceMutPointer<i8> {
        heap.allocate_model_storage(bytes.iter().map(|byte| *byte as i8).collect())
            .expect("fixture allocation")
    }

    fn c_string(heap: &SourceHeap, pointer: SourceConstPointer<i8>) -> Vec<u8> {
        heap.slice(pointer)
            .expect("fixture pointer")
            .iter()
            .take_while(|byte| **byte != 0)
            .map(|byte| *byte as u8)
            .collect()
    }

    #[test]
    fn source_port__ikey_dll__getinchikeyfrominchi__line_113() {
        let cases = [
            ("InChI=1S/CH4/h1H4\0", b'S', b'N'),
            ("InChI=1/CH4/h1H4\0", b'N', b'N'),
            ("InChI=1B/CH4/h1H4\0", b'B', b'N'),
            ("InChI=1/C2H6/c1-2/h1-2H3/q+1\0", b'N', b'N'),
            ("InChI=1/CH4/t1-/m0/s1\0", b'N', b'N'),
            ("InChI=1/CH4/fh1H4\0", b'N', b'N'),
            ("InChI=1/CH4/rC/h1H4\0", b'N', b'N'),
            ("InChI=1S//\0", b'S', b'N'),
            ("InChI=1S/?\0", b'S', b'N'),
            ("InChI=1S/CH4/h1H4 trailing bytes\0", b'S', b'N'),
            ("InChI=1S/CH4/h1H4/p+1\0", b'S', b'O'),
            ("InChI=1S/CH4/h1H4/p+12\0", b'S', b'Z'),
            ("InChI=1S/CH4/h1H4/p+13\0", b'S', b'A'),
            ("InChI=1S/CH4/h1H4/p-1\0", b'S', b'M'),
            ("InChI=1S/CH4/h1H4/p-12\0", b'S', b'B'),
            ("InChI=1S/CH4/h1H4/p-13\0", b'S', b'A'),
            ("InChI=1S/CH4/h1H4/p+1junk/t1-\0", b'S', b'O'),
            ("InChI=1S/CH4/h1H4/p9223372036854775808\0", b'S', b'M'),
        ];
        for (text, expected_standardness, expected_proton) in cases {
            let mut heap = SourceHeap::default();
            let input = allocate_bytes(&mut heap, text.as_bytes());
            let key = allocate_bytes(&mut heap, &[0xa5; 64]);
            let xtra1 = allocate_bytes(&mut heap, &[0xa5; 64]);
            let xtra2 = allocate_bytes(&mut heap, &[0xa5; 64]);
            let status = GetINCHIKeyFromINCHI(&mut heap, input.as_const(), 1, 1, key, xtra1, xtra2)
                .expect("source execution");
            assert_eq!(status, INCHIKEY_OK as i32, "{text}");
            let key_bytes = c_string(&heap, key.as_const());
            assert_eq!(key_bytes.len(), 27, "{text}");
            assert_eq!(key_bytes[23], expected_standardness, "{text}");
            assert_eq!(key_bytes[24], b'A', "{text}");
            assert_eq!(key_bytes[25], b'-', "{text}");
            assert_eq!(key_bytes[26], expected_proton, "{text}");
            assert_eq!(c_string(&heap, xtra1.as_const()).len(), 48, "{text}");
            assert_eq!(c_string(&heap, xtra2.as_const()).len(), 56, "{text}");
            assert_eq!(heap.slice(key.as_const()).expect("key")[28], 0xa5_u8 as i8);
            assert_eq!(
                heap.slice(xtra1.as_const()).expect("xtra1")[49],
                0xa5_u8 as i8
            );
            assert_eq!(
                heap.slice(xtra2.as_const()).expect("xtra2")[57],
                0xa5_u8 as i8
            );
            assert_eq!(heap.live_source_allocation_count(), 0, "{text}");
        }

        let mut heap = SourceHeap::default();
        let input = allocate_bytes(&mut heap, b"InChI=1S/CH4/h1H4\0");
        let key = allocate_bytes(&mut heap, &[0xa5; 64]);
        let xtra1 = allocate_bytes(&mut heap, &[0xa5; 64]);
        let xtra2 = allocate_bytes(&mut heap, &[0xa5; 64]);
        assert_eq!(
            GetINCHIKeyFromINCHI(&mut heap, input.as_const(), 0, 0, key, xtra1, xtra2)
                .expect("source execution"),
            INCHIKEY_OK as i32
        );
        assert_eq!(
            c_string(&heap, key.as_const()),
            b"VNWKTOKETHGBQD-UHFFFAOYSA-N"
        );
        assert_eq!(heap.slice(xtra1.as_const()).expect("xtra1")[0], 0);
        assert_eq!(
            heap.slice(xtra1.as_const()).expect("xtra1")[1],
            0xa5_u8 as i8
        );
        assert_eq!(heap.slice(xtra2.as_const()).expect("xtra2")[0], 0);
        assert_eq!(
            heap.slice(xtra2.as_const()).expect("xtra2")[1],
            0xa5_u8 as i8
        );
        assert_eq!(heap.live_source_allocation_count(), 0);

        for minor_length in [254_usize, 255] {
            let text = format!("InChI=1/C/{}\0", "x".repeat(minor_length - 1));
            let mut heap = SourceHeap::default();
            let input = allocate_bytes(&mut heap, text.as_bytes());
            let key = allocate_bytes(&mut heap, &[0xa5; 64]);
            assert_eq!(
                GetINCHIKeyFromINCHI(
                    &mut heap,
                    input.as_const(),
                    1,
                    1,
                    key,
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                )
                .expect("source execution"),
                INCHIKEY_OK as i32,
                "minor length {minor_length}"
            );
            assert_eq!(c_string(&heap, key.as_const()).len(), 27);
            assert_eq!(heap.live_source_allocation_count(), 0);
        }

        let mut heap = SourceHeap::default();
        let key = allocate_bytes(&mut heap, &[0xa5; 64]);
        let xtra1 = allocate_bytes(&mut heap, &[0xa5; 64]);
        let xtra2 = allocate_bytes(&mut heap, &[0xa5; 64]);
        assert_eq!(
            GetINCHIKeyFromINCHI(
                &mut heap,
                SourceConstPointer::null(),
                1,
                1,
                key,
                xtra1,
                xtra2,
            )
            .expect("source execution"),
            INCHIKEY_EMPTY_INPUT as i32
        );
        assert_eq!(heap.slice(key.as_const()).expect("key")[0], 0xa5_u8 as i8);
        assert_eq!(heap.slice(xtra1.as_const()).expect("xtra1")[0], 0);
        assert_eq!(heap.slice(xtra2.as_const()).expect("xtra2")[0], 0);

        for (text, expected) in [
            ("InChI=\0", INCHIKEY_INVALID_INCHI_PREFIX),
            ("inchi=1S/CH4\0", INCHIKEY_INVALID_INCHI_PREFIX),
            ("InChI=2S/CH4\0", INCHIKEY_INVALID_INCHI_PREFIX),
            ("InChI=1S:CH4\0", INCHIKEY_INVALID_INCHI_PREFIX),
            ("InChI=1S/\0", INCHIKEY_INVALID_INCHI),
            ("InChI=1S/!bad\0", INCHIKEY_INVALID_INCHI),
            ("InChI=1S/CH4/fh1\0", INCHIKEY_INVALID_STD_INCHI),
            ("InChI=1S/CH4/rC\0", INCHIKEY_INVALID_STD_INCHI),
            ("InChI=1S/CH4/p0\0", INCHIKEY_INVALID_STD_INCHI),
            (
                "InChI=1S/CH4/p-9223372036854775809\0",
                INCHIKEY_INVALID_STD_INCHI,
            ),
            ("InChI=1S/CH4/p\0", INCHIKEY_INVALID_STD_INCHI),
            ("InChI=1S/CH4/p/t1-\0", INCHIKEY_INVALID_INCHI),
        ] {
            let mut heap = SourceHeap::default();
            let input = allocate_bytes(&mut heap, text.as_bytes());
            let key = allocate_bytes(&mut heap, &[0xa5; 64]);
            assert_eq!(
                GetINCHIKeyFromINCHI(
                    &mut heap,
                    input.as_const(),
                    0,
                    0,
                    key,
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                )
                .expect("source execution"),
                expected as i32,
                "{text}"
            );
            assert_eq!(
                heap.slice(key.as_const()).expect("key")[0],
                if matches!(
                    expected,
                    INCHIKEY_INVALID_STD_INCHI | INCHIKEY_INVALID_INCHI
                ) && text.starts_with("InChI=1S/CH4/")
                {
                    0
                } else {
                    0xa5_u8 as i8
                },
                "{text}"
            );
            assert_eq!(heap.live_source_allocation_count(), 0, "{text}");
        }

        for successful_allocations in 0..5 {
            let mut heap = SourceHeap::default();
            let input = allocate_bytes(&mut heap, b"InChI=1S/CH4/h1H4\0");
            let key = allocate_bytes(&mut heap, &[0xa5; 64]);
            let xtra = allocate_bytes(&mut heap, &[0xa5; 64]);
            heap.fail_after_allocations(successful_allocations);
            assert_eq!(
                GetINCHIKeyFromINCHI(
                    &mut heap,
                    input.as_const(),
                    1,
                    0,
                    key,
                    xtra,
                    SourceMutPointer::null(),
                )
                .expect("source execution"),
                INCHIKEY_NOT_ENOUGH_MEMORY as i32
            );
            assert_eq!(heap.live_source_allocation_count(), 0);
            assert_eq!(heap.slice(xtra.as_const()).expect("xtra")[0], 0);
            assert_eq!(heap.slice(key.as_const()).expect("key")[0], 0xa5_u8 as i8);
        }
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    pub(crate) fn official_c_oracle__getinchikeyfrominchi__exact() {
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
            .arg("--inchi-key-records")
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
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "get_inchi_key_from_inchi");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_bytes = |value: &Value, expected_length: usize| {
                let values = value
                    .as_array()
                    .unwrap_or_else(|| panic!("{case_id}: bytes must be an array"));
                assert_eq!(values.len(), expected_length, "{case_id}");
                values
                    .iter()
                    .map(|value| {
                        u8::try_from(value.as_u64().expect("byte must be unsigned"))
                            .unwrap_or_else(|_| panic!("{case_id}: byte exceeds u8"))
                    })
                    .collect::<Vec<_>>()
            };

            let input_before = parse_bytes(&official["input"]["bytes"], 640);
            let key_before = parse_bytes(&official["input"]["key"], 64);
            let xtra1_before = parse_bytes(&official["input"]["xtra1_bytes"], 64);
            let xtra2_before = parse_bytes(&official["input"]["xtra2_bytes"], 64);
            let mut heap = SourceHeap::default();
            let input = allocate_bytes(&mut heap, &input_before);
            let key = allocate_bytes(&mut heap, &key_before);
            let xtra1_pointer = allocate_bytes(&mut heap, &xtra1_before);
            let xtra2_pointer = allocate_bytes(&mut heap, &xtra2_before);
            let allocation_failure_ordinal = official["input"]["allocation_failure_ordinal"]
                .as_u64()
                .expect("allocation_failure_ordinal must be unsigned");
            if allocation_failure_ordinal == 0 {
                heap.trace_source_allocations();
            } else {
                heap.fail_after_allocations(allocation_failure_ordinal - 1);
            }
            let input_pointer_null = official["input"]["input_pointer_null"]
                .as_bool()
                .expect("input_pointer_null must be boolean");
            let xtra1_pointer_null = official["input"]["xtra1_pointer_null"]
                .as_bool()
                .expect("xtra1_pointer_null must be boolean");
            let xtra2_pointer_null = official["input"]["xtra2_pointer_null"]
                .as_bool()
                .expect("xtra2_pointer_null must be boolean");
            let status = GetINCHIKeyFromINCHI(
                &mut heap,
                if input_pointer_null {
                    SourceConstPointer::null()
                } else {
                    input.as_const()
                },
                i32::try_from(
                    official["input"]["xtra1"]
                        .as_i64()
                        .expect("xtra1 must be signed"),
                )
                .expect("xtra1 exceeds i32"),
                i32::try_from(
                    official["input"]["xtra2"]
                        .as_i64()
                        .expect("xtra2 must be signed"),
                )
                .expect("xtra2 exceeds i32"),
                key,
                if xtra1_pointer_null {
                    SourceMutPointer::null()
                } else {
                    xtra1_pointer
                },
                if xtra2_pointer_null {
                    SourceMutPointer::null()
                } else {
                    xtra2_pointer
                },
            )
            .unwrap_or_else(|error| panic!("{case_id}: source execution failed: {error:?}"));

            assert_eq!(
                i64::from(status),
                official["output"]["status"]
                    .as_i64()
                    .expect("status must be signed"),
                "{case_id}"
            );
            let rust_input = heap
                .slice(input.as_const())
                .expect("input storage")
                .iter()
                .map(|byte| *byte as u8)
                .collect::<Vec<_>>();
            assert_eq!(
                rust_input,
                parse_bytes(&official["output"]["bytes"], 640),
                "{case_id}"
            );
            assert_eq!(rust_input, input_before, "{case_id}");
            assert_eq!(official["output"]["input_unchanged"], true, "{case_id}");
            for (pointer, field) in [
                (key.as_const(), "key"),
                (xtra1_pointer.as_const(), "xtra1_bytes"),
                (xtra2_pointer.as_const(), "xtra2_bytes"),
            ] {
                let rust_bytes = heap
                    .slice(pointer)
                    .unwrap_or_else(|error| panic!("{case_id}: {field}: {error:?}"))
                    .iter()
                    .map(|byte| *byte as u8)
                    .collect::<Vec<_>>();
                assert_eq!(
                    rust_bytes,
                    parse_bytes(&official["output"][field], 64),
                    "{case_id}: {field}"
                );
            }

            let official_allocation_calls = official["output"]["allocation_calls"]
                .as_u64()
                .expect("allocation_calls must be unsigned");
            assert_eq!(
                heap.source_allocation_calls(),
                official_allocation_calls,
                "{case_id}"
            );
            let failure_observed = u64::from(
                allocation_failure_ordinal > 0
                    && official_allocation_calls >= allocation_failure_ordinal,
            );
            assert_eq!(
                official["output"]["deferred_frees"]
                    .as_u64()
                    .expect("deferred_frees must be unsigned"),
                official_allocation_calls - failure_observed,
                "{case_id}"
            );
            assert_eq!(official["output"]["live_allocations"], 0, "{case_id}");
            assert_eq!(heap.live_source_allocation_count(), 0, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 38);
    }
}
