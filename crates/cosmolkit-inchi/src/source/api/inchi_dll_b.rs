use crate::source::api::inchi_dll::parse_options_string;
use crate::source::base::ichi_io::{
    inchi_ios_close, inchi_ios_eprint, inchi_ios_gets, inchi_ios_getsTab, inchi_ios_getsTab1,
    inchi_ios_init, inchi_strbuf_close, inchi_strbuf_init,
};
use crate::source::base::ichican2::SetBitFree;
use crate::source::base::ichierr::AddErrorMessage;
use crate::source::base::ichiparm::ReadCommandLineParms;
use crate::source::base::mol2atom::FreeOrigAtData;
use crate::source::base::readinch::{FindToken, FreeInchi_Stereo0D, LoadLine};
use crate::source::base::runichi4::FreeAllINChIArrays;
use crate::source::base::util::{inchi_calloc, inchi_free, inchi_malloc, lrtrim, mystrncpy};
use crate::source_types::{
    AB_MAX_WELL_DEFINED_PARITY, AB_MIN_WELL_DEFINED_PARITY, CANON_GLOBALS, CHAR_MASK,
    FLAG_INP_AT_CHIRAL, FLAG_INP_AT_NONCHIRAL, INCHI_INP_EOF_ERR, INCHI_INP_EOF_RET,
    INCHI_INP_ERROR_ERR, INCHI_INP_ERROR_RET, INCHI_INP_FATAL_ERR, INCHI_INP_FATAL_RET,
    INCHI_IOS_STRING, INCHI_IOS_TYPE_STRING, INCHI_IOSTREAM, INCHI_LINE_ADD, INCHI_LINE_LEN,
    INCHI_MAX_NUM_ARG, INCHI_MODE, INCHI_NUM, INCHI_OUT_NO_AUX_INFO, INCHI_STRBUF_INITIAL_SIZE,
    INCHI_STRBUF_SIZE_INCREMENT, INPUT_PARMS, INPUT_TYPE, MAX_INPUT_BOND_TYPE, MAX_NUM_PATHS,
    MAX_SDF_HEADER, MAX_SDF_VALUE, MIN_INPUT_BOND_TYPE, MOL_COORD, MOL2INCHI_BAD_COMMAND_LINE,
    MOL2INCHI_NO_RAM, NO_ATOM, ORIG_ATOM_DATA, PINChI_Aux2, PINChI2, S_SHORT, SB_PARITY_MASK,
    SB_PARITY_SHFT, STRUCT_DATA, SourceArgvPointer, SourceConstPointer, SourceHeap,
    SourceHeapError, SourceMutPointer, SourceVaList, bRELEASE_VERSION, inchi_Atom, inchi_Input,
    inchi_Stereo0D, tagINCHIBondStereo2D_INCHI_BOND_STEREO_DOUBLE_EITHER,
    tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1DOWN,
    tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1EITHER,
    tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1UP,
    tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2DOWN,
    tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2EITHER,
    tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2UP, tagINCHIBondType_INCHI_BOND_TYPE_ALTERN,
    tagINCHIBondType_INCHI_BOND_TYPE_DOUBLE, tagINCHIBondType_INCHI_BOND_TYPE_SINGLE,
    tagINCHIBondType_INCHI_BOND_TYPE_TRIPLE, tagINCHIStereoParity0D_INCHI_PARITY_EVEN,
    tagINCHIStereoParity0D_INCHI_PARITY_ODD, tagINCHIStereoParity0D_INCHI_PARITY_UNDEFINED,
    tagINCHIStereoParity0D_INCHI_PARITY_UNKNOWN, tagINCHIStereoType0D_INCHI_StereoType_Allene,
    tagINCHIStereoType0D_INCHI_StereoType_DoubleBond, tagINCHIStereoType0D_INCHI_StereoType_None,
    tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral, tagInputType_INPUT_INCHI_PLAIN,
    tagInputType_INPUT_INCHI_XML, tagInputType_INPUT_MOLFILE,
};

// `inchi_api.h` fields under the pinned LP64 ABI: three doubles, three
// 20-element source arrays, and the trailing scalar fields total 120 bytes.
const SOURCE_SIZEOF_INCHI_ATOM: u64 = 3 * 8 + 20 * 2 + 20 + 20 + 6 + 2 + 4 + 2 + 1 + 1;

fn source_strtol_decimal(bytes: &[i8], start: usize) -> (i64, usize) {
    let mut cursor = start;
    while bytes
        .get(cursor)
        .is_some_and(|byte| matches!(*byte as u8, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r'))
    {
        cursor += 1;
    }

    let negative = match bytes.get(cursor).map(|byte| *byte as u8) {
        Some(b'-') => {
            cursor += 1;
            true
        }
        Some(b'+') => {
            cursor += 1;
            false
        }
        _ => false,
    };
    let digit_start = cursor;
    let limit = if negative {
        i64::MAX as u64 + 1
    } else {
        i64::MAX as u64
    };
    let mut magnitude = 0_u64;
    let mut overflowed = false;
    while let Some(digit) = bytes
        .get(cursor)
        .map(|byte| *byte as u8)
        .filter(u8::is_ascii_digit)
    {
        let digit = u64::from(digit - b'0');
        match magnitude
            .checked_mul(10)
            .and_then(|value| value.checked_add(digit))
        {
            Some(value) if value <= limit => magnitude = value,
            _ => overflowed = true,
        }
        cursor += 1;
    }
    if cursor == digit_start {
        return (0, start);
    }
    if overflowed {
        return (if negative { i64::MIN } else { i64::MAX }, cursor);
    }
    let value = if negative {
        if magnitude == i64::MAX as u64 + 1 {
            i64::MIN
        } else {
            -(magnitude as i64)
        }
    } else {
        magnitude as i64
    };
    (value, cursor)
}

fn source_strtod_decimal(bytes: &[i8], start: usize) -> (f64, usize) {
    let mut cursor = start;
    while bytes
        .get(cursor)
        .is_some_and(|byte| matches!(*byte as u8, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r'))
    {
        cursor += 1;
    }
    let number_start = cursor;
    if matches!(bytes.get(cursor).map(|byte| *byte as u8), Some(b'+' | b'-')) {
        cursor += 1;
    }
    let integer_start = cursor;
    while bytes
        .get(cursor)
        .is_some_and(|byte| (*byte as u8).is_ascii_digit())
    {
        cursor += 1;
    }
    let mut has_digits = cursor != integer_start;
    if bytes.get(cursor).map(|byte| *byte as u8) == Some(b'.') {
        cursor += 1;
        let fraction_start = cursor;
        while bytes
            .get(cursor)
            .is_some_and(|byte| (*byte as u8).is_ascii_digit())
        {
            cursor += 1;
        }
        has_digits |= cursor != fraction_start;
    }
    if !has_digits {
        return (0.0, start);
    }
    if matches!(bytes.get(cursor).map(|byte| *byte as u8), Some(b'e' | b'E')) {
        let exponent_mark = cursor;
        cursor += 1;
        if matches!(bytes.get(cursor).map(|byte| *byte as u8), Some(b'+' | b'-')) {
            cursor += 1;
        }
        let exponent_start = cursor;
        while bytes
            .get(cursor)
            .is_some_and(|byte| (*byte as u8).is_ascii_digit())
        {
            cursor += 1;
        }
        if cursor == exponent_start {
            cursor = exponent_mark;
        }
    }
    let text: String = bytes[number_start..cursor]
        .iter()
        .map(|byte| *byte as u8 as char)
        .collect();
    (text.parse::<f64>().unwrap_or(0.0), cursor)
}

fn source_find_bytes(haystack: &[i8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() {
        return Some(0);
    }
    haystack.windows(needle.len()).position(|window| {
        window
            .iter()
            .zip(needle)
            .all(|(left, right)| *left as u8 == *right)
    })
}

fn source_starts_with(haystack: &[i8], needle: &[u8]) -> bool {
    haystack.len() >= needle.len()
        && haystack[..needle.len()]
            .iter()
            .zip(needle)
            .all(|(left, right)| *left as u8 == *right)
}

fn source_c_string(bytes: &[i8]) -> Result<&[i8], SourceHeapError> {
    let length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    Ok(&bytes[..length])
}

pub(crate) fn add_source_error(
    heap: &mut SourceHeap,
    error_buffer: Option<SourceMutPointer<i8>>,
    message: &str,
) -> Result<(), SourceHeapError> {
    let Some(error_buffer) = error_buffer else {
        return Ok(());
    };
    let mut source: Vec<i8> = message.bytes().map(|byte| byte as i8).collect();
    source.push(0);
    AddErrorMessage(Some(heap.slice_mut(error_buffer)?), Some(&source))?;
    Ok(())
}

pub(crate) fn FreeInchi_Atom(
    heap: &mut SourceHeap,
    atom_slot: Option<&mut SourceMutPointer<inchi_Atom>>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:106 FreeInchi_Atom
    // INCHI✔️❌: void FreeInchi_Atom( inchi_Atom **at )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (at && *at)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( *at );
    // INCHI✔️❌:         *at = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: FreeInchi_Atom

    let Some(atom_slot) = atom_slot else {
        return Ok(());
    };
    if atom_slot.is_null() {
        return Ok(());
    }
    inchi_free(heap, *atom_slot)?;
    *atom_slot = SourceMutPointer::null();
    Ok(())
}

pub(crate) fn CreateInchiAtom(
    heap: &mut SourceHeap,
    num_atoms: i32,
) -> Result<SourceMutPointer<inchi_Atom>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:117 CreateInchiAtom
    // INCHI✔️❌: inchi_Atom *CreateInchiAtom( int num_atoms )
    // INCHI✔️❌: {
    // INCHI✔️❌:     inchi_Atom *p = (inchi_Atom*) inchi_calloc( num_atoms, sizeof( inchi_Atom ) );
    // INCHI✔️❌:     return p;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CreateInchiAtom

    inchi_calloc(heap, num_atoms as u64, SOURCE_SIZEOF_INCHI_ATOM)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn PrepareToMakeINCHI(
    heap: &mut SourceHeap,
    sd: &mut STRUCT_DATA,
    ip: &mut INPUT_PARMS,
    orig_inp_data: &mut ORIG_ATOM_DATA,
    prep_inp_data: &mut [ORIG_ATOM_DATA; INCHI_NUM as usize],
    p_inchi: &mut [SourceMutPointer<PINChI2>; INCHI_NUM as usize],
    p_inchi_aux: &mut [SourceMutPointer<PINChI_Aux2>; INCHI_NUM as usize],
    pout: &mut INCHI_IOSTREAM,
    plog: &mut INCHI_IOSTREAM,
    pprb: &mut INCHI_IOSTREAM,
    inp_file: &mut INCHI_IOSTREAM,
    moltext: SourceConstPointer<i8>,
    options: SourceConstPointer<i8>,
    strbuf: &mut INCHI_IOS_STRING,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:285 PrepareToMakeINCHI
    // INCHI✔️❌: complete source frame follows verbatim; SourceHeap stack-buffer and checked pointer modeling add overhead.
    /*
int PrepareToMakeINCHI( STRUCT_DATA *sd,
                        INPUT_PARMS *ip,
                        ORIG_ATOM_DATA *orig_inp_data,
                        ORIG_ATOM_DATA *prep_inp_data,
                        PINChI2 *pINChI[INCHI_NUM],
                        PINChI_Aux2 *pINChI_Aux[INCHI_NUM],
                        INCHI_IOSTREAM *pout,
                        INCHI_IOSTREAM *plog,
                        INCHI_IOSTREAM *pprb,
                        INCHI_IOSTREAM *inp_file,
                        const char *moltext,
                        char *options,
                        INCHI_IOS_STRING *strbuf )
{
    int retcode = 0;
    unsigned long  ulDisplTime = 0;
    int bReleaseVersion = bRELEASE_VERSION;
    char szSdfDataValue[MAX_SDF_VALUE + 1];


    const char *quasi_argv[INCHI_MAX_NUM_ARG + 1];
    int   quasi_argc;
    char *quasi_options = NULL;

    if (options)
    {
        quasi_options = (char*) inchi_malloc( strlen( options ) + 1 );
    }
    if (quasi_options)
    {
        strcpy( quasi_options, options );
        quasi_argc = parse_options_string( quasi_options, quasi_argv, INCHI_MAX_NUM_ARG );
    }
    else
    {
        quasi_argc = 1;
        quasi_argv[0] = "";
        quasi_argv[1] = NULL;
    }

    /* I/O streams */
    inchi_ios_init( pout, INCHI_IOS_TYPE_STRING, NULL );
    inchi_ios_init( plog, INCHI_IOS_TYPE_STRING, NULL );
    inchi_ios_init( pprb, INCHI_IOS_TYPE_STRING, NULL );


    /* input ( string of Molfile )*/
    inchi_ios_init( inp_file, INCHI_IOS_TYPE_STRING, NULL );
    inp_file->s.pStr = (char *) moltext;
    inp_file->s.nPtr = 0;
    inp_file->s.nUsedLength = strlen( moltext ) + 1;
    inp_file->f = NULL;

    memset( szSdfDataValue, 0, sizeof( szSdfDataValue ) ); /* djb-rwth: memset_s C11/Annex K variant? */


    /* data structs */
    memset( sd, 0, sizeof( *sd ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( ip, 0, sizeof( *ip ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    memset( orig_inp_data, 0, sizeof( *orig_inp_data ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( prep_inp_data, 0, 2 * sizeof( *prep_inp_data ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    pINChI[0] = pINChI[1] = NULL;
    pINChI_Aux[0] = pINChI_Aux[1] = NULL;

    /* Parse command line */
    if (0 > ReadCommandLineParms( quasi_argc,
        quasi_argv,
        ip,
        szSdfDataValue,
        &ulDisplTime,
        bReleaseVersion,
        plog ))
    {
        return  MOL2INCHI_BAD_COMMAND_LINE;
    }

    ip->nInputType = INPUT_MOLFILE;
    ip->bNoStructLabels = 1;
    ip->pSdfLabel = NULL;
    ip->pSdfValue = NULL;
    /* ip->bINChIOutputOptions |= INCHI_OUT_NO_AUX_INFO; */

    /* Supply expandable string buffer */
    if (0 >= inchi_strbuf_init( strbuf, INCHI_STRBUF_INITIAL_SIZE, INCHI_STRBUF_SIZE_INCREMENT ))
    {
        if (plog && plog->s.pStr)
        {
            inchi_ios_eprint( plog, "Cannot allocate output string buffer. Terminating\n" );
        }
        retcode = MOL2INCHI_NO_RAM;
        goto ret;
    }

ret:
    if (quasi_options)
    {
        inchi_free( quasi_options );
    }

    return retcode;
}
    */
    // END INCHI C FUNCTION: PrepareToMakeINCHI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: PrepareToMakeINCHI
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; bRELEASE_VERSION == 1.
    // INCHI✔️❌: inchi_malloc/free are the active libc allocation macros; INCHI_NUM == 2.
    // INCHI✔️❌: INPUT_MOLFILE resolves to tagInputType_INPUT_MOLFILE.
    // END INCHI ACTIVE MACRO CONFIGURATION: PrepareToMakeINCHI

    let moltext_length = heap
        .slice(moltext)?
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;

    let mut quasi_options = SourceMutPointer::null();
    let mut parsed_arguments = vec![SourceArgvPointer::Null; INCHI_MAX_NUM_ARG as usize + 1];
    let quasi_argc;
    if !options.is_null() {
        let option_length = heap
            .slice(options)?
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        quasi_options = match inchi_malloc(heap, (option_length + 1) as u64) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        if !quasi_options.is_null() {
            let copied = heap.slice(options)?[..=option_length].to_vec();
            heap.slice_mut(quasi_options)?[..=option_length].copy_from_slice(&copied);
        }
    }
    if !quasi_options.is_null() {
        quasi_argc = parse_options_string(
            heap,
            quasi_options,
            &mut parsed_arguments,
            INCHI_MAX_NUM_ARG as i32,
        )?;
    } else {
        quasi_argc = 1;
        parsed_arguments[0] = SourceArgvPointer::EmptyLiteral;
        parsed_arguments[1] = SourceArgvPointer::Null;
    }

    inchi_ios_init(
        Some(pout),
        INCHI_IOS_TYPE_STRING as i32,
        SourceMutPointer::null(),
    )?;
    inchi_ios_init(
        Some(plog),
        INCHI_IOS_TYPE_STRING as i32,
        SourceMutPointer::null(),
    )?;
    inchi_ios_init(
        Some(pprb),
        INCHI_IOS_TYPE_STRING as i32,
        SourceMutPointer::null(),
    )?;
    inchi_ios_init(
        Some(inp_file),
        INCHI_IOS_TYPE_STRING as i32,
        SourceMutPointer::null(),
    )?;
    inp_file.s.pStr = moltext.as_mut();
    inp_file.s.nPtr = 0;
    inp_file.s.nUsedLength = moltext_length.wrapping_add(1) as i32;
    inp_file.f = SourceMutPointer::null();

    *sd = STRUCT_DATA::default();
    *ip = INPUT_PARMS::default();
    *orig_inp_data = ORIG_ATOM_DATA::default();
    *prep_inp_data = std::array::from_fn(|_| ORIG_ATOM_DATA::default());
    *p_inchi = std::array::from_fn(|_| SourceMutPointer::null());
    *p_inchi_aux = std::array::from_fn(|_| SourceMutPointer::null());

    let sdf_value = heap.allocate_model_storage(vec![0_i8; MAX_SDF_VALUE as usize + 1])?;
    let empty_argument = heap.allocate_model_storage(vec![0_i8])?;
    let mut argv = Vec::with_capacity(quasi_argc as usize);
    for argument in parsed_arguments.iter().take(quasi_argc as usize) {
        argv.push(match *argument {
            SourceArgvPointer::EmptyLiteral => empty_argument.as_const(),
            SourceArgvPointer::Command(pointer) => pointer.as_const(),
            SourceArgvPointer::Null => SourceConstPointer::null(),
        });
    }
    let mut display_time = 0_u64;
    let command_result = ReadCommandLineParms(
        heap,
        quasi_argc,
        &argv,
        ip,
        sdf_value,
        &mut display_time,
        bRELEASE_VERSION as i32,
        Some(plog),
    );
    heap.free(empty_argument)?;
    heap.free(sdf_value)?;
    let command_result = command_result?;
    if command_result < 0 {
        return Ok(MOL2INCHI_BAD_COMMAND_LINE as i32);
    }

    ip.nInputType = tagInputType_INPUT_MOLFILE;
    ip.bNoStructLabels = 1;
    ip.pSdfLabel = SourceMutPointer::null();
    ip.pSdfValue = SourceMutPointer::null();

    let mut retcode = 0_i32;
    if inchi_strbuf_init(
        heap,
        strbuf,
        INCHI_STRBUF_INITIAL_SIZE as i32,
        INCHI_STRBUF_SIZE_INCREMENT as i32,
    )? <= 0
    {
        if !plog.s.pStr.is_null() {
            let format = heap.allocate_model_storage(
                b"Cannot allocate output string buffer. Terminating\n\0"
                    .iter()
                    .map(|byte| *byte as i8)
                    .collect(),
            )?;
            let print_result = inchi_ios_eprint(
                heap,
                Some(plog),
                format.as_const(),
                &SourceVaList::default(),
            );
            heap.free(format)?;
            print_result?;
        }
        retcode = MOL2INCHI_NO_RAM as i32;
    }

    if !quasi_options.is_null() {
        inchi_free(heap, quasi_options)?;
    }
    Ok(retcode)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn PostMakeINCHICleanup(
    heap: &mut SourceHeap,
    canonical_globals: &mut CANON_GLOBALS,
    sd: &mut STRUCT_DATA,
    ip: &mut INPUT_PARMS,
    orig_inp_data: &mut ORIG_ATOM_DATA,
    prep_inp_data: &mut [ORIG_ATOM_DATA; INCHI_NUM as usize],
    p_inchi: &mut [SourceMutPointer<PINChI2>; INCHI_NUM as usize],
    p_inchi_aux: &mut [SourceMutPointer<PINChI_Aux2>; INCHI_NUM as usize],
    pout: &mut INCHI_IOSTREAM,
    _plog: &mut INCHI_IOSTREAM,
    pprb: &mut INCHI_IOSTREAM,
    _inp_file: &mut INCHI_IOSTREAM,
    _moltext: SourceConstPointer<i8>,
    strbuf: &mut INCHI_IOS_STRING,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:391 PostMakeINCHICleanup
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap frees add overhead.
    /*
int PostMakeINCHICleanup( struct tagCANON_GLOBALS *pCG,
                          STRUCT_DATA *sd,
                          INPUT_PARMS *ip,
                          ORIG_ATOM_DATA *orig_inp_data,
                          ORIG_ATOM_DATA *prep_inp_data,
                          PINChI2 *pINChI[INCHI_NUM],
                          PINChI_Aux2 *pINChI_Aux[INCHI_NUM],
                          INCHI_IOSTREAM *pout,
                          INCHI_IOSTREAM *plog,
                          INCHI_IOSTREAM *pprb,
                          INCHI_IOSTREAM *inp_file,
                          const char *moltext,
                          INCHI_IOS_STRING *strbuf )
{
    int retcode = 0;
    int i;


    /* Free structure data */

    /*  Free INChI memory */
    FreeAllINChIArrays( pINChI, pINChI_Aux, sd->num_components );

    FreeOrigAtData( orig_inp_data );
    FreeOrigAtData( prep_inp_data );
    FreeOrigAtData( prep_inp_data + 1 );

    inchi_ios_close( pout );
    inchi_ios_close( pprb );

    inchi_strbuf_close( strbuf );

    for (i = 0; i < MAX_NUM_PATHS; i++)
    {
        if (ip->path[i])
        {
            inchi_free( (void*) ip->path[i] );
            /*  cast deliberately discards 'const' qualifier */
            ip->path[i] = NULL;
        }
    }
    SetBitFree( pCG );

    return retcode;
}
    */
    // END INCHI C FUNCTION: PostMakeINCHICleanup
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: PostMakeINCHICleanup
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; INCHI_NUM == 2 and MAX_NUM_PATHS == 4.
    // INCHI✔️❌: plog, inp_file, and moltext are intentionally not closed or freed by this function.
    // END INCHI ACTIVE MACRO CONFIGURATION: PostMakeINCHICleanup

    FreeAllINChIArrays(heap, p_inchi, p_inchi_aux, &mut sd.num_components)?;
    FreeOrigAtData(heap, Some(orig_inp_data))?;
    FreeOrigAtData(heap, Some(&mut prep_inp_data[0]))?;
    FreeOrigAtData(heap, Some(&mut prep_inp_data[1]))?;
    inchi_ios_close(heap, Some(pout))?;
    inchi_ios_close(heap, Some(pprb))?;
    inchi_strbuf_close(heap, Some(strbuf))?;
    for path in ip.path.iter_mut().take(MAX_NUM_PATHS as usize) {
        if !path.is_null() {
            inchi_free(heap, path.as_mut())?;
            *path = SourceConstPointer::null();
        }
    }
    SetBitFree(heap, canonical_globals)?;
    Ok(0)
}

pub(crate) fn FreeInchi_Input(
    heap: &mut SourceHeap,
    input: &mut inchi_Input,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:439 FreeInchi_Input
    // INCHI✔️❌: void FreeInchi_Input( inchi_Input *inp_at_data )
    // INCHI✔️❌: {
    // INCHI✔️❌:     FreeInchi_Atom( &inp_at_data->atom );
    // INCHI✔️❌:     FreeInchi_Stereo0D( &inp_at_data->stereo0D );
    // INCHI✔️❌:     memset( inp_at_data, 0, sizeof( *inp_at_data ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: FreeInchi_Input

    FreeInchi_Atom(heap, Some(&mut input.atom))?;
    FreeInchi_Stereo0D(heap, Some(&mut input.stereo0D))?;
    *input = inchi_Input::default();
    Ok(())
}

pub(crate) fn is_in_the_slist(
    path_atom: Option<&[S_SHORT]>,
    next_atom: S_SHORT,
    mut path_length: i32,
) -> Result<Option<usize>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:531 is_in_the_slist
    // INCHI✔️❌: S_SHORT *is_in_the_slist( S_SHORT *pathAtom, S_SHORT nNextAtom, int nPathLen )
    // INCHI✔️❌: {
    // INCHI✔️❌:     for (; nPathLen && *pathAtom != nNextAtom; nPathLen--, pathAtom++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return nPathLen ? pathAtom : NULL;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: is_in_the_slist

    if path_length == 0 {
        return Ok(None);
    }
    let path_atom = path_atom.ok_or(SourceHeapError::NullPointer)?;
    let mut position = 0_usize;
    while path_length != 0 {
        let current = *path_atom
            .get(position)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if current == next_atom {
            return Ok(Some(position));
        }
        path_length = path_length
            .checked_sub(1)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        position = position
            .checked_add(1)
            .ok_or(SourceHeapError::PointerOffsetOverflow)?;
    }
    Ok(None)
}

pub(crate) fn is_element_a_metal(element: &[i8]) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:543 is_element_a_metal
    // INCHI✔️✔️: int is_element_a_metal( char szEl[] )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     static const char szMetals[] = "K;V;Y;W;U;"
    // INCHI✔️✔️:         "Li;Be;Na;Mg;Al;Ca;Sc;Ti;Cr;Mn;Fe;Co;Ni;Cu;Zn;Ga;Rb;Sr;Zr;"
    // INCHI✔️✔️:         "Nb;Mo;Tc;Ru;Rh;Pd;Ag;Cd;In;Sn;Sb;Cs;Ba;La;Ce;Pr;Nd;Pm;Sm;"
    // INCHI✔️✔️:         "Eu;Gd;Tb;Dy;Ho;Er;Tm;Yb;Lu;Hf;Ta;Re;Os;Ir;Pt;Au;Hg;Tl;Pb;"
    // INCHI✔️✔️:         "Bi;Po;Fr;Ra;Ac;Th;Pa;Np;Pu;Am;Cm;Bk;Cf;Es;Fm;Md;No;Lr;Rf;";
    // INCHI✔️✔️:     const int len = (int) strlen( szEl );
    // INCHI✔️✔️:     const char *p;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (0 < len && len <= 2 &&
    // INCHI✔️✔️:          isalpha( UCINT szEl[0] ) && isupper( szEl[0] ) &&
    // INCHI✔️✔️:          ( p = strstr( szMetals, szEl ) ) && p[len] == ';')
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:
    // INCHI✔️✔️:         return 1; /*return AtType_Metal;*/
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: is_element_a_metal

    const METALS: &[u8] = b"K;V;Y;W;U;\
Li;Be;Na;Mg;Al;Ca;Sc;Ti;Cr;Mn;Fe;Co;Ni;Cu;Zn;Ga;Rb;Sr;Zr;\
Nb;Mo;Tc;Ru;Rh;Pd;Ag;Cd;In;Sn;Sb;Cs;Ba;La;Ce;Pr;Nd;Pm;Sm;\
Eu;Gd;Tb;Dy;Ho;Er;Tm;Yb;Lu;Hf;Ta;Re;Os;Ir;Pt;Au;Hg;Tl;Pb;\
Bi;Po;Fr;Ra;Ac;Th;Pa;Np;Pu;Am;Cm;Bk;Cf;Es;Fm;Md;No;Lr;Rf;";

    let length = element
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    if length == 0 || length > 2 {
        return Ok(0);
    }

    let first = element[0] as u8;
    if !first.is_ascii_alphabetic() || !first.is_ascii_uppercase() {
        return Ok(0);
    }

    let position = METALS.windows(length).position(|candidate| {
        candidate
            .iter()
            .zip(&element[..length])
            .all(|(candidate, needle)| *candidate == *needle as u8)
    });
    Ok(i32::from(position.is_some_and(|position| {
        METALS.get(position + length) == Some(&b';')
    })))
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn InchiToInchiAtom(
    heap: &mut SourceHeap,
    input: Option<&mut INCHI_IOSTREAM>,
    mut stereo0d_slot: Option<&mut SourceMutPointer<inchi_Stereo0D>>,
    mut stereo0d_count: Option<&mut i32>,
    do_not_add_h: i32,
    _ab_parity_unknown: i32,
    input_type: INPUT_TYPE,
    mut atom_slot: Option<&mut SourceMutPointer<inchi_Atom>>,
    maximum_atom_count: i32,
    dimensions: Option<&mut i32>,
    bond_count: Option<&mut i32>,
    sdf_label: Option<SourceMutPointer<i8>>,
    sdf_value: Option<SourceMutPointer<i8>>,
    mut id: Option<&mut i64>,
    mut input_atom_flags: Option<&mut INCHI_MODE>,
    error: Option<&mut i32>,
    error_text: Option<SourceMutPointer<i8>>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:582 InchiToInchiAtom
    // INCHI✔❌: int InchiToInchiAtom( INCHI_IOSTREAM *inp_file,
    // INCHI✔❌:                       inchi_Stereo0D **stereo0D,
    // INCHI✔❌:                       int *num_stereo0D,
    // INCHI✔❌:                       int bDoNotAddH,
    // INCHI✔❌:                       int vABParityUnknown,
    // INCHI✔❌:                       INPUT_TYPE nInputType,
    // INCHI✔❌:                       inchi_Atom **at,
    // INCHI✔❌:                       int max_num_at,
    // INCHI✔❌:                       int *num_dimensions,
    // INCHI✔❌:                       int *num_bonds,
    // INCHI✔❌:                       char *pSdfLabel,
    // INCHI✔❌:                       char *pSdfValue,
    // INCHI✔❌:                       long *Id,
    // INCHI✔❌:                       INCHI_MODE *pInpAtomFlags,
    // INCHI✔❌:                       int *err,
    // INCHI✔❌:                       char *pStrErr )
    // INCHI✔❌: {
    // INCHI✔❌:     int      num_atoms = 0, bFindNext = 0, len = 0, bHeaderRead, bItemIsOver, bErrorMsg, bRestoreInfo;; /* djb-rwth: removing redundant variables; initialising variables -- updated 28/09/2025 */
    // INCHI✔❌:     int      bFatal = 0, num_struct = 0;
    // INCHI✔❌:     int      i, k, k2, res, bond_type, bond_stereo1, bond_stereo2, bond_char, neigh, bond_parity, bond_parityNM;
    // INCHI✔❌:     int      bTooLongLine, res2, bTooLongLine2, pos, hlen, hk;
    // INCHI✔❌:     long     longID;
    // INCHI✔❌:     char     *p, *q, *s, parity;
    // INCHI✔❌:     int      b2D = 0, b3D = 0, b23D, nNumBonds = 0, bNonZeroXYZ, bNonMetal;
    // INCHI✔❌:     int      len_stereo0D = 0, max_len_stereo0D = 0;
    // INCHI✔❌:     inchi_Stereo0D  *atom_stereo0D = NULL;
    // INCHI✔❌:     inchi_Atom      *atom = NULL;
    // INCHI✔❌:     MOL_COORD       *pszCoord = NULL;
    // INCHI✔❌:     INCHI_MODE InpAtomFlags = 0; /* 0 or FLAG_INP_AT_NONCHIRAL or FLAG_INP_AT_CHIRAL */
    // INCHI✔❌:     static const char szIsoH[] = "hdt";
    // INCHI✔❌:     /* plain tags */
    // INCHI✔❌:     static const char sStructHdrPln[] = "Structure:";
    // INCHI✔❌:     static const char sStructHdrPlnNoLblVal[] = " is missing";
    // INCHI✔❌:     static char sStructHdrPlnAuxStart[64] = ""; /*"$1.1Beta/";*/
    // INCHI✔❌:     static int  lenStructHdrPlnAuxStart = 0;
    // INCHI✔❌:     static const char sStructHdrPlnRevAt[] = "/rA:";
    // INCHI✔❌:     static const char sStructHdrPlnRevBn[] = "/rB:";
    // INCHI✔❌:     static const char sStructHdrPlnRevXYZ[] = "/rC:";
    // INCHI✔❌:     const  char *sToken;
    // INCHI✔❌:     int  lToken;
    // INCHI✔❌:
    // INCHI✔❌:     if (!lenStructHdrPlnAuxStart)
    // INCHI✔❌:     {
    // INCHI✔❌:         lenStructHdrPlnAuxStart = sprintf( sStructHdrPlnAuxStart, "AuxInfo=" );
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (at)
    // INCHI✔❌:     {
    // INCHI✔❌:
    // INCHI✔❌:         if (*at && max_num_at)
    // INCHI✔❌:         {
    // INCHI✔❌:             memset( *at, 0, max_num_at * sizeof( **at ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         if (stereo0D && num_stereo0D)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (*stereo0D && *num_stereo0D)
    // INCHI✔❌:             {
    // INCHI✔❌:                 max_len_stereo0D = *num_stereo0D;
    // INCHI✔❌:                 memset( *stereo0D, 0, max_len_stereo0D * sizeof( **stereo0D ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 max_len_stereo0D = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:     else  /* if ( at )  */
    // INCHI✔❌:     {
    // INCHI✔❌:         bFindNext = 1;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     bHeaderRead = bErrorMsg = bRestoreInfo = 0;
    // INCHI✔❌:     *num_dimensions = *num_bonds = 0;
    // INCHI✔❌:
    // INCHI✔❌:     /*************************************************************/
    // INCHI✔❌:     /*   extract reversibility info from plain text INChI format */
    // INCHI✔❌:     /*************************************************************/
    // INCHI✔❌:
    // INCHI✔❌:     if (nInputType == INPUT_INCHI_PLAIN)
    // INCHI✔❌:     {
    // INCHI✔❌:
    // INCHI✔❌:         bHeaderRead = hk = 0; /* djb-rwth: ignoring LLVM warning: value used */
    // INCHI✔❌:
    // INCHI✔❌:         while (0 < ( res = inchi_ios_getsTab(szLine_i2ia, sizeof(szLine_i2ia) - 1, inp_file, &bTooLongLine ) ))
    // INCHI✔❌:         {
    // INCHI✔❌:
    // INCHI✔❌:             /********************* find and interpret structure header ************/
    // INCHI✔❌:             if (!bTooLongLine &&
    // INCHI✔❌:                 ( hlen = sizeof( sStructHdrPln ) - 1, !memcmp(szLine_i2ia, sStructHdrPln, hlen ) ))
    // INCHI✔❌:             {
    // INCHI✔❌:                 p = szLine_i2ia + hlen;
    // INCHI✔❌:                 longID = 0;
    // INCHI✔❌:                 num_atoms = 0;
    // INCHI✔❌:
    // INCHI✔❌:                 /* structure number */
    // INCHI✔❌:                 longID = strtol( p, &q, 10 );
    // INCHI✔❌:                 if (q && q[0] == '.' && q[1] == ' ')
    // INCHI✔❌:                     p = q + 2;
    // INCHI✔❌:                 p = p + strspn( p, " \n\r" );
    // INCHI✔❌:
    // INCHI✔❌:                 if (pSdfLabel)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     pSdfLabel[0] = '\0';
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 if (pSdfValue)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     pSdfValue[0] = '\0';
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌:                 if (*p)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     /* has label name */
    // INCHI✔❌:
    // INCHI✔❌:                     /*p ++;*/
    // INCHI✔❌:                     if ((q = strchr( p, '=' ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:                     {
    // INCHI✔❌:
    // INCHI✔❌:                         /* '=' separates label name from the value */
    // INCHI✔❌:                         len = inchi_min( q - p + 1, MAX_SDF_HEADER - 1 );
    // INCHI✔❌:
    // INCHI✔❌:                         if (pSdfLabel)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             mystrncpy( pSdfLabel, p, len );
    // INCHI✔❌:                             lrtrim( pSdfLabel, &len );
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         p = q + 1;
    // INCHI✔❌:                         q = p + (int) strlen( p );
    // INCHI✔❌:
    // INCHI✔❌:                         if (q - p > 0)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             len = inchi_min( q - p + 1, MAX_SDF_VALUE - 1 );
    // INCHI✔❌:                             if (pSdfValue)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 mystrncpy( pSdfValue, p, len );
    // INCHI✔❌:                             }
    // INCHI✔❌:                             p = q; /* djb-rwth: ignoring LLVM warning: value used */
    // INCHI✔❌:                         }
    // INCHI✔❌:                     }
    // INCHI✔❌:                     else if ((q = strstr( p, sStructHdrPlnNoLblVal ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:                     {
    // INCHI✔❌:                         len = inchi_min( q - p + 1, MAX_SDF_HEADER - 1 );
    // INCHI✔❌:                         if (pSdfLabel)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             mystrncpy( pSdfLabel, p, len );
    // INCHI✔❌:                         }
    // INCHI✔❌:                         p = q + 1; /* djb-rwth: ignoring LLVM warning: value used */
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 if (Id)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     *Id = longID;
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 bHeaderRead = 1;
    // INCHI✔❌:                 bErrorMsg = bRestoreInfo = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (!memcmp(szLine_i2ia, sStructHdrPlnAuxStart, lenStructHdrPlnAuxStart ))
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* found the header of the AuxInfo, read AuxInfo head of the line */
    // INCHI✔❌:
    // INCHI✔❌:                 if (!bHeaderRead)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     longID = 0;
    // INCHI✔❌:                     if (Id)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         *Id = longID;
    // INCHI✔❌:                     }
    // INCHI✔❌:                     if (pSdfLabel)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         pSdfLabel[0] = '\0';
    // INCHI✔❌:                     }
    // INCHI✔❌:                     if (pSdfValue)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         pSdfValue[0] = '\0';
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 bHeaderRead = 0; /* djb-rwth: ignoring LLVM warning: value used */
    // INCHI✔❌:
    // INCHI✔❌:                 /* check for empty "AuxInfo=ver//" */
    // INCHI✔❌:
    // INCHI✔❌:                 p = strchr(szLine_i2ia + lenStructHdrPlnAuxStart, '/' );
    // INCHI✔❌:
    // INCHI✔❌:                 if (p && p[1] == '/' && ( !p[2] || '\n' == p[2] ))
    // INCHI✔❌:                 {
    // INCHI✔❌:                     goto bypass_end_of_INChI_plain;
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 /***************** search for atoms block (plain) **********************/
    // INCHI✔❌:
    // INCHI✔❌:                 p = szLine_i2ia;
    // INCHI✔❌:                 sToken = sStructHdrPlnRevAt;
    // INCHI✔❌:                 lToken = sizeof( sStructHdrPlnRevAt ) - 1;
    // INCHI✔❌:
    // INCHI✔❌:                 /* search for sToken in the line; load next segments of the line if sToken has not found */
    // INCHI✔❌:
    // INCHI✔❌:                 p = FindToken( inp_file, &bTooLongLine, sToken, lToken,
    // INCHI✔❌:                     szLine_i2ia, sizeof(szLine_i2ia), p, &res );
    // INCHI✔❌:
    // INCHI✔❌:                 if (!p)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                     num_atoms = INCHI_INP_ERROR_RET;
    // INCHI✔❌:                     TREAT_ERR( *err, 0, "Missing atom data" );
    // INCHI✔❌:                     goto bypass_end_of_INChI_plain;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                     /* atoms block started */
    // INCHI✔❌:
    // INCHI✔❌:                     i = 0;
    // INCHI✔❌:                     res2 = bTooLongLine2 = -1; /* djb-rwth: ignoring LLVM warning: value used */
    // INCHI✔❌:                     bItemIsOver = ( s = strchr( p, '/' ) ) || !bTooLongLine;
    // INCHI✔❌:
    // INCHI✔❌:                     while (1)
    // INCHI✔❌:                     {
    // INCHI✔❌:
    // INCHI✔❌:                         p = LoadLine( inp_file, &bTooLongLine, &bItemIsOver, &s,
    // INCHI✔❌:                             szLine_i2ia, sizeof(szLine_i2ia), INCHI_LINE_ADD, p, &res );
    // INCHI✔❌:
    // INCHI✔❌:                         if (!i)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             /* allocate atom */
    // INCHI✔❌:                             num_atoms = strtol( p, &q, 10 );
    // INCHI✔❌:
    // INCHI✔❌:                             if (!num_atoms || !q || !*q)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 num_atoms = 0; /* no atom data */
    // INCHI✔❌:                                 goto bypass_end_of_INChI_plain;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             p = q;
    // INCHI✔❌:
    // INCHI✔❌:                             /* Molfile chirality flag */
    // INCHI✔❌:                             switch (*p)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 case 'c':
    // INCHI✔❌:                                     InpAtomFlags |= FLAG_INP_AT_CHIRAL;
    // INCHI✔❌:                                     p++;
    // INCHI✔❌:                                     break;
    // INCHI✔❌:                                 case 'n':
    // INCHI✔❌:                                     InpAtomFlags |= FLAG_INP_AT_NONCHIRAL;
    // INCHI✔❌:                                     p++;
    // INCHI✔❌:                                     break;
    // INCHI✔❌:                             }
    // INCHI✔❌:
    // INCHI✔❌:                             if (at && *at)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 if (num_atoms > max_num_at)
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     inchi_free( *at );
    // INCHI✔❌:                                     *at = NULL;
    // INCHI✔❌:                                 }
    // INCHI✔❌:                                 else
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     memset( *at, 0, max_num_at * sizeof( **at ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:                                     atom = *at;
    // INCHI✔❌:                                 }
    // INCHI✔❌:                             }
    // INCHI✔❌:
    // INCHI✔❌:                             if (!at || !*at)
    // INCHI✔❌:                             {
    // INCHI✔❌:
    // INCHI✔❌:                                 atom = CreateInchiAtom( num_atoms + 1 );
    // INCHI✔❌:
    // INCHI✔❌:                                 if (!atom)
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     num_atoms = INCHI_INP_FATAL_RET; /* was -1; error */
    // INCHI✔❌:                                     *err = INCHI_INP_FATAL_ERR;
    // INCHI✔❌:                                     TREAT_ERR( *err, 0, "Out of RAM" );
    // INCHI✔❌:                                     goto bypass_end_of_INChI_plain;
    // INCHI✔❌:                                 }
    // INCHI✔❌:                             }
    // INCHI✔❌:
    // INCHI✔❌:                             if (stereo0D && *stereo0D)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 if (num_atoms > max_len_stereo0D)
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     FreeInchi_Stereo0D( stereo0D );
    // INCHI✔❌:                                 }
    // INCHI✔❌:                                 else
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     memset( *stereo0D, 0, max_len_stereo0D * sizeof( **stereo0D ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:                                     atom_stereo0D = *stereo0D;
    // INCHI✔❌:                                 }
    // INCHI✔❌:                             }
    // INCHI✔❌:
    // INCHI✔❌:                             if (!stereo0D || !*stereo0D)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 max_len_stereo0D = num_atoms + 1;
    // INCHI✔❌:
    // INCHI✔❌:                                 atom_stereo0D = CreateInchi_Stereo0D( max_len_stereo0D );
    // INCHI✔❌:
    // INCHI✔❌:                                 if (!atom_stereo0D)
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     num_atoms = INCHI_INP_FATAL_RET; /* fatal error: cannot allocate */
    // INCHI✔❌:                                     *err = INCHI_INP_FATAL_ERR;
    // INCHI✔❌:                                     TREAT_ERR( *err, 0, "Out of RAM" );
    // INCHI✔❌:                                     /* djb-rwth: avoiding memory leak */
    // INCHI✔❌:                                     inchi_free(atom);
    // INCHI✔❌:                                     atom = NULL;
    // INCHI✔❌:                                     goto bypass_end_of_INChI_plain;
    // INCHI✔❌:                                 }
    // INCHI✔❌:                             }
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         /* element, first char */
    // INCHI✔❌:                         if (!isalpha( UCINT *p ) || !isupper( UCINT *p ) || i >= num_atoms)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             break; /* end of atoms block */
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         atom[i].elname[0] = *p++;
    // INCHI✔❌:
    // INCHI✔❌:                         /* element, second char */
    // INCHI✔❌:                         if (isalpha( UCINT *p ) && islower( UCINT *p ))
    // INCHI✔❌:                         {
    // INCHI✔❌:                             atom[i].elname[1] = *p++;
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         /* bonds' valence + number of non-isotopic H */
    // INCHI✔❌:                         if (isdigit( UCINT *p ))
    // INCHI✔❌:                         {
    // INCHI✔❌:                             AT_BONDS_VAL( atom, i ) = (char) strtol( p, &q, 10 );
    // INCHI✔❌:                             if (!AT_BONDS_VAL( atom, i ))
    // INCHI✔❌:                                 AT_BONDS_VAL( atom, i ) = ISOLATED_ATOM; /* same convention as in MOLfile, found zero bonds valence */
    // INCHI✔❌:                             p = q;
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         /* charge */
    // INCHI✔❌:                         atom[i].charge = ( *p == '+' ) ? 1 : ( *p == '-' ) ? -1 : 0;
    // INCHI✔❌:                         if (atom[i].charge)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             p++;
    // INCHI✔❌:                             if (isdigit( UCINT *p ))
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 atom[i].charge *= (S_CHAR) ( strtol( p, &q, 10 ) & CHAR_MASK );
    // INCHI✔❌:                                 p = q;
    // INCHI✔❌:                             }
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         /* radical */
    // INCHI✔❌:                         if (*p == '.')
    // INCHI✔❌:                         {
    // INCHI✔❌:                             p++;
    // INCHI✔❌:                             if (isdigit( UCINT *p ))
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 atom[i].radical = (S_CHAR) strtol( p, &q, 10 );
    // INCHI✔❌:                                 p = q;
    // INCHI✔❌:                             }
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         /* isotopic mass */
    // INCHI✔❌:                         if (*p == 'i')
    // INCHI✔❌:                         {
    // INCHI✔❌:                             p++;
    // INCHI✔❌:                             if (isdigit( UCINT *p ))
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 int mw = strtol( p, &q, 10 );
    // INCHI✔❌:                                 p = q;
    // INCHI✔❌:
    // INCHI✔❌:                                 atom[i].isotopic_mass = mw;
    // INCHI✔❌:                             }
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         /* parity */
    // INCHI✔❌:                         switch (*p)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             case 'o':
    // INCHI✔❌:                                 parity = INCHI_PARITY_ODD;
    // INCHI✔❌:                                 p++;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 'e':
    // INCHI✔❌:                                 parity = INCHI_PARITY_EVEN;
    // INCHI✔❌:                                 p++;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 'u':
    // INCHI✔❌:                                 parity = INCHI_PARITY_UNKNOWN;
    // INCHI✔❌:                                 p++;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case '?':
    // INCHI✔❌:                                 parity = INCHI_PARITY_UNDEFINED;
    // INCHI✔❌:                                 p++;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             default:
    // INCHI✔❌:                                 parity = 0;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         if (parity)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             atom_stereo0D[len_stereo0D].central_atom = i;
    // INCHI✔❌:                             atom_stereo0D[len_stereo0D].parity = parity;
    // INCHI✔❌:                             atom_stereo0D[len_stereo0D].type = INCHI_StereoType_Tetrahedral;
    // INCHI✔❌:                             len_stereo0D++;
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         /* isotopic h, d, t */
    // INCHI✔❌:                         for (k = 0; k < NUM_H_ISOTOPES; k++)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             if (*p == szIsoH[k])
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 NUM_ISO_Hk( atom, i, k ) = 1;
    // INCHI✔❌:                                 p++;
    // INCHI✔❌:                                 if (isdigit( UCINT *p ))
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     NUM_ISO_Hk( atom, i, k ) = (char) strtol( p, &q, 10 );
    // INCHI✔❌:                                     p = q;
    // INCHI✔❌:                                 }
    // INCHI✔❌:                             }
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         i++;
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌:                     if (!bItemIsOver || i != num_atoms || (s && p != s)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:                     {
    // INCHI✔❌:                         num_atoms = INCHI_INP_ERROR_RET; /* error */
    // INCHI✔❌:                         *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                         TREAT_ERR( *err, 0, "Wrong number of atoms" );
    // INCHI✔❌:                         /* djb-rwth: avoiding memory leak */
    // INCHI✔❌:                         inchi_free(atom);
    // INCHI✔❌:                         atom = NULL;
    // INCHI✔❌:                         goto bypass_end_of_INChI_plain;
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 /***************** search for bonds block (plain) and read it *****************/
    // INCHI✔❌:
    // INCHI✔❌:                 /*p = szLine;*/
    // INCHI✔❌:                 sToken = sStructHdrPlnRevBn;
    // INCHI✔❌:                 lToken = sizeof( sStructHdrPlnRevBn ) - 1;
    // INCHI✔❌:
    // INCHI✔❌:                 /* search for sToken in the line; load next segments of the line if sToken has not found */
    // INCHI✔❌:
    // INCHI✔❌:                 p = FindToken( inp_file, &bTooLongLine, sToken, lToken,
    // INCHI✔❌:                     szLine_i2ia, sizeof(szLine_i2ia), p, &res );
    // INCHI✔❌:
    // INCHI✔❌:                 if (!p)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     num_atoms = INCHI_INP_ERROR_RET; /* error */
    // INCHI✔❌:                     *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                     TREAT_ERR( *err, 0, "Missing bonds data" );
    // INCHI✔❌:                     goto bypass_end_of_INChI_plain;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                     /* bonds block started */
    // INCHI✔❌:
    // INCHI✔❌:                     i = 1;
    // INCHI✔❌:
    // INCHI✔❌:                     res2 = bTooLongLine2 = -1; /* djb-rwth: ignoring LLVM warning: value used */
    // INCHI✔❌:
    // INCHI✔❌:                     bItemIsOver = ( s = strchr( p, '/' ) ) || !bTooLongLine;
    // INCHI✔❌:
    // INCHI✔❌:                     if (1 == num_atoms)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         /* needed because the next '/' may be still out of szLine */
    // INCHI✔❌:
    // INCHI✔❌:                         p = LoadLine( inp_file, &bTooLongLine, &bItemIsOver, &s,
    // INCHI✔❌:                             szLine_i2ia, sizeof(szLine_i2ia), INCHI_LINE_ADD, p, &res );
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     while (i < num_atoms)
    // INCHI✔❌:                     {
    // INCHI✔❌:
    // INCHI✔❌:                         p = LoadLine( inp_file, &bTooLongLine, &bItemIsOver, &s,
    // INCHI✔❌:                             szLine_i2ia, sizeof(szLine_i2ia), INCHI_LINE_ADD, p, &res );
    // INCHI✔❌:
    // INCHI✔❌:                         if (i >= num_atoms || (s && p >= s)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:                         {
    // INCHI✔❌:                             break; /* end of bonds (plain) */
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         /* bond, first char */
    // INCHI✔❌:                         if (*p == ';')
    // INCHI✔❌:                         {
    // INCHI✔❌:                             p++;
    // INCHI✔❌:                             i++;
    // INCHI✔❌:                             continue;
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         if (!isalpha( UCINT *p ))
    // INCHI✔❌:                         {
    // INCHI✔❌:                             num_atoms = INCHI_INP_ERROR_RET; /* error */
    // INCHI✔❌:                             *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                             TREAT_ERR( *err, 0, "Wrong bonds data" );
    // INCHI✔❌:                             goto bypass_end_of_INChI_plain;
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         bond_char = *p++;
    // INCHI✔❌:
    // INCHI✔❌:                         /* bond parity */
    // INCHI✔❌:                         switch (*p)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             case '-':
    // INCHI✔❌:                                 bond_parity = INCHI_PARITY_ODD;
    // INCHI✔❌:                                 p++;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case '+':
    // INCHI✔❌:                                 bond_parity = INCHI_PARITY_EVEN;
    // INCHI✔❌:                                 p++;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 'u':
    // INCHI✔❌:                                 bond_parity = INCHI_PARITY_UNKNOWN;
    // INCHI✔❌:                                 p++;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case '?':
    // INCHI✔❌:                                 bond_parity = INCHI_PARITY_UNDEFINED;
    // INCHI✔❌:                                 p++;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             default:
    // INCHI✔❌:                                 bond_parity = 0;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         if (bond_parity)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             switch (*p)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 case '-':
    // INCHI✔❌:                                     bond_parityNM = INCHI_PARITY_ODD;
    // INCHI✔❌:                                     p++;
    // INCHI✔❌:                                     break;
    // INCHI✔❌:                                 case '+':
    // INCHI✔❌:                                     bond_parityNM = INCHI_PARITY_EVEN;
    // INCHI✔❌:                                     p++;
    // INCHI✔❌:                                     break;
    // INCHI✔❌:                                 case 'u':
    // INCHI✔❌:                                     bond_parityNM = INCHI_PARITY_UNKNOWN;
    // INCHI✔❌:                                     p++;
    // INCHI✔❌:                                     break;
    // INCHI✔❌:                                 case '?':
    // INCHI✔❌:                                     bond_parityNM = INCHI_PARITY_UNDEFINED;
    // INCHI✔❌:                                     p++;
    // INCHI✔❌:                                     break;
    // INCHI✔❌:                                 default:
    // INCHI✔❌:                                     bond_parityNM = 0;
    // INCHI✔❌:                                     break;
    // INCHI✔❌:                             }
    // INCHI✔❌:                         }
    // INCHI✔❌:                         else
    // INCHI✔❌:                         {
    // INCHI✔❌:                             bond_parityNM = 0;
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         /* neighbor of the current atom */
    // INCHI✔❌:                         if (!isdigit( UCINT *p ))
    // INCHI✔❌:                         {
    // INCHI✔❌:                             num_atoms = INCHI_INP_ERROR_RET; /* error */
    // INCHI✔❌:                             *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                             TREAT_ERR( *err, 0, "Wrong bonds data" );
    // INCHI✔❌:                             goto bypass_end_of_INChI_plain;
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         neigh = (int) strtol( p, &q, 10 ) - 1;
    // INCHI✔❌:
    // INCHI✔❌:                         if (i >= num_atoms || neigh >= num_atoms)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             num_atoms = INCHI_INP_ERROR_RET; /* error */
    // INCHI✔❌:                             *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                             TREAT_ERR( *err, 0, "Bond to nonexistent atom" );
    // INCHI✔❌:                             goto bypass_end_of_INChI_plain;
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         p = q;
    // INCHI✔❌:                         bond_stereo1 = bond_stereo2 = 0;
    // INCHI✔❌:
    // INCHI✔❌:                         /* bond type & 2D stereo */
    // INCHI✔❌:                         switch (bond_char)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             case 'v':
    // INCHI✔❌:                                 bond_type = INCHI_BOND_TYPE_SINGLE;
    // INCHI✔❌:                                 bond_stereo1 = INCHI_BOND_STEREO_SINGLE_1EITHER;
    // INCHI✔❌:                                 bond_stereo2 = INCHI_BOND_STEREO_SINGLE_2EITHER;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 'V':
    // INCHI✔❌:                                 bond_type = INCHI_BOND_TYPE_SINGLE;
    // INCHI✔❌:                                 bond_stereo1 = INCHI_BOND_STEREO_SINGLE_2EITHER;
    // INCHI✔❌:                                 bond_stereo2 = INCHI_BOND_STEREO_SINGLE_1EITHER;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 'w':
    // INCHI✔❌:                                 bond_type = INCHI_BOND_TYPE_DOUBLE;
    // INCHI✔❌:                                 bond_stereo1 =
    // INCHI✔❌:                                     bond_stereo2 = INCHI_BOND_STEREO_DOUBLE_EITHER;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 's':
    // INCHI✔❌:                                 bond_type = INCHI_BOND_TYPE_SINGLE;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 'd':
    // INCHI✔❌:                                 bond_type = INCHI_BOND_TYPE_DOUBLE;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 't':
    // INCHI✔❌:                                 bond_type = INCHI_BOND_TYPE_TRIPLE;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 'a':
    // INCHI✔❌:                                 bond_type = INCHI_BOND_TYPE_ALTERN;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 'p':
    // INCHI✔❌:                                 bond_type = INCHI_BOND_TYPE_SINGLE;
    // INCHI✔❌:                                 bond_stereo1 = INCHI_BOND_STEREO_SINGLE_1UP;
    // INCHI✔❌:                                 bond_stereo2 = INCHI_BOND_STEREO_SINGLE_2UP;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 'P':
    // INCHI✔❌:                                 bond_type = INCHI_BOND_TYPE_SINGLE;
    // INCHI✔❌:                                 bond_stereo1 = INCHI_BOND_STEREO_SINGLE_2UP;
    // INCHI✔❌:                                 bond_stereo2 = INCHI_BOND_STEREO_SINGLE_1UP;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 'n':
    // INCHI✔❌:                                 bond_type = INCHI_BOND_TYPE_SINGLE;
    // INCHI✔❌:                                 bond_stereo1 = INCHI_BOND_STEREO_SINGLE_1DOWN;
    // INCHI✔❌:                                 bond_stereo2 = INCHI_BOND_STEREO_SINGLE_2DOWN;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 'N':
    // INCHI✔❌:                                 bond_type = INCHI_BOND_TYPE_SINGLE;
    // INCHI✔❌:                                 bond_stereo1 = INCHI_BOND_STEREO_SINGLE_2DOWN;
    // INCHI✔❌:                                 bond_stereo2 = INCHI_BOND_STEREO_SINGLE_1DOWN;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             default:
    // INCHI✔❌:                                 num_atoms = INCHI_INP_ERROR_RET; /* error */
    // INCHI✔❌:                                 *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                                 TREAT_ERR( *err, 0, "Wrong bond type" );
    // INCHI✔❌:                                 goto bypass_end_of_INChI_plain;
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         k = AT_NUM_BONDS( atom[i] )++;
    // INCHI✔❌:
    // INCHI✔❌:                         atom[i].bond_type[k] = bond_type;
    // INCHI✔❌:                         atom[i].bond_stereo[k] = bond_stereo1;
    // INCHI✔❌:                         atom[i].neighbor[k] = (ATOM_NUMBER) neigh;
    // INCHI✔❌:
    // INCHI✔❌:                         k2 = AT_NUM_BONDS( atom[neigh] )++;
    // INCHI✔❌:                         atom[neigh].bond_type[k2] = bond_type;
    // INCHI✔❌:                         atom[neigh].bond_stereo[k2] = bond_stereo2;
    // INCHI✔❌:                         atom[neigh].neighbor[k2] = (ATOM_NUMBER) i;
    // INCHI✔❌:
    // INCHI✔❌:                         bond_parity |= ( bond_parityNM << SB_PARITY_SHFT );
    // INCHI✔❌:
    // INCHI✔❌:                         if (bond_parity)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             if (max_len_stereo0D <= len_stereo0D)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 /* realloc atom_Stereo0D */
    // INCHI✔❌:
    // INCHI✔❌:                                 inchi_Stereo0D *new_atom_stereo0D = CreateInchi_Stereo0D( max_len_stereo0D + num_atoms );
    // INCHI✔❌:
    // INCHI✔❌:                                 if (!new_atom_stereo0D)
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     num_atoms = INCHI_INP_FATAL_RET; /* fatal error: cannot allocate */
    // INCHI✔❌:                                     *err = INCHI_INP_FATAL_ERR;
    // INCHI✔❌:                                     TREAT_ERR( *err, 0, "Out of RAM" );
    // INCHI✔❌:                                     goto bypass_end_of_INChI_plain;
    // INCHI✔❌:                                 }
    // INCHI✔❌:
    // INCHI✔❌:                                 memcpy(new_atom_stereo0D, atom_stereo0D, len_stereo0D * sizeof(*atom_stereo0D));
    // INCHI✔❌:                                 FreeInchi_Stereo0D( &atom_stereo0D );
    // INCHI✔❌:                                 atom_stereo0D = new_atom_stereo0D;
    // INCHI✔❌:                                 max_len_stereo0D += num_atoms;
    // INCHI✔❌:                             }
    // INCHI✔❌:
    // INCHI✔❌:                             /* (a) i may be allene endpoint and     neigh = allene middle point or
    // INCHI✔❌:                                 (b) i may be allene middle point and neigh = allene endpoint
    // INCHI✔❌:                                 !!!!! CURRENTLY ONLY (b) IS ALLOWED !!!!!
    // INCHI✔❌:                             */
    // INCHI✔❌:
    // INCHI✔❌:                             atom_stereo0D[len_stereo0D].neighbor[1] = neigh; /* neigh < i */
    // INCHI✔❌:                             atom_stereo0D[len_stereo0D].neighbor[2] = i;
    // INCHI✔❌:                             atom_stereo0D[len_stereo0D].parity = bond_parity;
    // INCHI✔❌:                             atom_stereo0D[len_stereo0D].type = INCHI_StereoType_DoubleBond; /* incl allenes & cumulenes */
    // INCHI✔❌:                             len_stereo0D++;
    // INCHI✔❌:                         }
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     if (!bItemIsOver || i != num_atoms || (s && p != s)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:                     {
    // INCHI✔❌:                         num_atoms = INCHI_INP_ERROR_RET; /* error */
    // INCHI✔❌:                         *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                         TREAT_ERR( *err, 0, "Wrong number of bonds" );
    // INCHI✔❌:                         goto bypass_end_of_INChI_plain;
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 /***************** search for coordinates block (plain) **********************/
    // INCHI✔❌:                 /*p = szLine;*/
    // INCHI✔❌:
    // INCHI✔❌:                 sToken = sStructHdrPlnRevXYZ;
    // INCHI✔❌:                 lToken = sizeof( sStructHdrPlnRevXYZ ) - 1;
    // INCHI✔❌:
    // INCHI✔❌:                 /* search for sToken in the line; load next segments of the line if sToken has not found */
    // INCHI✔❌:
    // INCHI✔❌:                 p = FindToken( inp_file, &bTooLongLine, sToken, lToken,
    // INCHI✔❌:                     szLine_i2ia, sizeof(szLine_i2ia), p, &res );
    // INCHI✔❌:
    // INCHI✔❌:                 if (!p)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     num_atoms = INCHI_INP_ERROR_RET; /* error */
    // INCHI✔❌:                     *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                     TREAT_ERR( *err, 0, "Missing atom coordinates data" );
    // INCHI✔❌:                     goto bypass_end_of_INChI_plain;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                     /* coordinates block started */
    // INCHI✔❌:                     if ((pszCoord = (MOL_COORD*) inchi_malloc( inchi_max( num_atoms, 1 ) * sizeof( MOL_COORD ) ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:                     {
    // INCHI✔❌:                         memset( pszCoord, ' ', inchi_max( num_atoms, 1 ) * sizeof( MOL_COORD ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:                     }
    // INCHI✔❌:                     else
    // INCHI✔❌:                     {
    // INCHI✔❌:                         num_atoms = INCHI_INP_FATAL_RET; /* allocation error */
    // INCHI✔❌:                         *err = INCHI_INP_FATAL_ERR;
    // INCHI✔❌:                         TREAT_ERR( *err, 0, "Out of RAM" );
    // INCHI✔❌:                         goto bypass_end_of_INChI_plain;
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     i = 0;
    // INCHI✔❌:                     res2 = bTooLongLine2 = -1; /* djb-rwth: ignoring LLVM warning: value used */
    // INCHI✔❌:                     bItemIsOver = ( s = strchr( p, '/' ) ) || !bTooLongLine;
    // INCHI✔❌:
    // INCHI✔❌:                     while (i < num_atoms)
    // INCHI✔❌:                     {
    // INCHI✔❌:
    // INCHI✔❌:                         p = LoadLine( inp_file, &bTooLongLine, &bItemIsOver, &s,
    // INCHI✔❌:                             szLine_i2ia, sizeof(szLine_i2ia), INCHI_LINE_ADD, p, &res );
    // INCHI✔❌:
    // INCHI✔❌:                         if (i >= num_atoms || (s && p >= s)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:                         {
    // INCHI✔❌:                             break; /* end of bonds (plain) */
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         /* coord, first char */
    // INCHI✔❌:                         if (*p == ';')
    // INCHI✔❌:                         {
    // INCHI✔❌:                             for (k = 0; k < NUM_COORD; k++)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 pszCoord[i][LEN_COORD*k + 4] = '0';
    // INCHI✔❌:                             }
    // INCHI✔❌:                             p++;
    // INCHI✔❌:                             i++;
    // INCHI✔❌:                             continue;
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         for (k = 0; k < 3; k++)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             double xyz;
    // INCHI✔❌:                             bNonZeroXYZ = 0;
    // INCHI✔❌:                             if (*p == ';')
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 pszCoord[i][LEN_COORD*k + 4] = '0';
    // INCHI✔❌:                                 xyz = 0.0;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             else
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 if (*p == ',')
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     /* empty */
    // INCHI✔❌:                                     pszCoord[i][LEN_COORD*k + 4] = '0';
    // INCHI✔❌:                                     xyz = 0.0;
    // INCHI✔❌:                                     p++;
    // INCHI✔❌:                                 }
    // INCHI✔❌:                                 else
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     xyz = strtod( p, &q );
    // INCHI✔❌:                                     bNonZeroXYZ = fabs( xyz ) > MIN_BOND_LENGTH;
    // INCHI✔❌:                                     if (q != NULL)
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         memcpy( pszCoord[i] + LEN_COORD*(long long)k, p, q - p ); /* djb-rwth: cast operator added */
    // INCHI✔❌:                                         if (*q == ',')
    // INCHI✔❌:                                             q++;
    // INCHI✔❌:                                         p = q;
    // INCHI✔❌:                                     }
    // INCHI✔❌:                                     else
    // INCHI✔❌:                                         pszCoord[i][LEN_COORD*k + 4] = '0';
    // INCHI✔❌:                                 }
    // INCHI✔❌:                             }
    // INCHI✔❌:
    // INCHI✔❌:                             switch (k)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 case 0:
    // INCHI✔❌:                                     atom[i].x = xyz;
    // INCHI✔❌:                                     b2D |= bNonZeroXYZ;
    // INCHI✔❌:                                     break;
    // INCHI✔❌:                                 case 1:
    // INCHI✔❌:                                     atom[i].y = xyz;
    // INCHI✔❌:                                     b2D |= bNonZeroXYZ;
    // INCHI✔❌:                                     break;
    // INCHI✔❌:                                 case 2:
    // INCHI✔❌:                                     b3D |= bNonZeroXYZ;
    // INCHI✔❌:                                     atom[i].z = xyz;
    // INCHI✔❌:                                     break;
    // INCHI✔❌:                             }
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         if (*p == ';')
    // INCHI✔❌:                         {
    // INCHI✔❌:                             p++; /* end of this triple of coordinates */
    // INCHI✔❌:                             i++;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         else
    // INCHI✔❌:                         {
    // INCHI✔❌:                             num_atoms = INCHI_INP_ERROR_RET; /* error in input data: atoms, bonds & coord must be present together */
    // INCHI✔❌:                             *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                             TREAT_ERR( *err, 0, "Wrong atom coordinates data" );
    // INCHI✔❌:                             goto bypass_end_of_INChI_plain;
    // INCHI✔❌:                         }
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     if (!bItemIsOver || (s && p != s) || i != num_atoms) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:                     {
    // INCHI✔❌:                         num_atoms = INCHI_INP_ERROR_RET; /* error */
    // INCHI✔❌:                         *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                         TREAT_ERR( *err, 0, "Wrong number of coordinates" );
    // INCHI✔❌:                         goto bypass_end_of_INChI_plain;
    // INCHI✔❌:                     }
    // INCHI✔❌:                 } /* end of coordinates */
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌:                 /* set special valences and implicit H (xml) */
    // INCHI✔❌:
    // INCHI✔❌:                 b23D = b2D | b3D;
    // INCHI✔❌:                 b2D = b3D = 0; /* djb-rwth: ignoring LLVM warning: values used */
    // INCHI✔❌:                 if (at)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     if (!*at)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         int a1, a2, n1, n2, valence;
    // INCHI✔❌:                         int chem_bonds_valence;
    // INCHI✔❌:                         int    nX = 0, nY = 0, nZ = 0, nXYZ;
    // INCHI✔❌:                         *at = atom;
    // INCHI✔❌:
    // INCHI✔❌:                         /* special valences */
    // INCHI✔❌:
    // INCHI✔❌:                         for (bNonMetal = 0; bNonMetal < 1; bNonMetal++)
    // INCHI✔❌:                         {
    // INCHI✔❌:
    // INCHI✔❌:                             for (a1 = 0; a1 < num_atoms; a1++)
    // INCHI✔❌:                             {
    // INCHI✔❌:
    // INCHI✔❌:                                 int num_bond_type[MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE + 1];
    // INCHI✔❌:
    // INCHI✔❌:                                 memset( num_bond_type, 0, sizeof( num_bond_type ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:
    // INCHI✔❌:                                 valence = AT_BONDS_VAL( atom, a1 ); /*  save atom valence if available */
    // INCHI✔❌:                                 AT_BONDS_VAL( atom, a1 ) = 0;
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌:                                 nX = nY = nZ = 0;
    // INCHI✔❌:
    // INCHI✔❌:                                 for (n1 = 0; n1 < AT_NUM_BONDS( atom[a1] ); n1++)
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     bond_type = atom[a1].bond_type[n1] - MIN_INPUT_BOND_TYPE;
    // INCHI✔❌:                                     if (bond_type < 0 || bond_type > MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE)
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         bond_type = 0;
    // INCHI✔❌:                                         TREAT_ERR( *err, 0, "Unknown bond type in InChI aux assigned as a single bond" );
    // INCHI✔❌:                                     }
    // INCHI✔❌:
    // INCHI✔❌:                                     num_bond_type[bond_type] ++;
    // INCHI✔❌:                                     nNumBonds++;
    // INCHI✔❌:                                     if (b23D)
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         neigh = atom[a1].neighbor[n1];
    // INCHI✔❌:                                         nX |= ( fabs( atom[a1].x - atom[neigh].x ) > MIN_BOND_LENGTH );
    // INCHI✔❌:                                         nY |= ( fabs( atom[a1].y - atom[neigh].y ) > MIN_BOND_LENGTH );
    // INCHI✔❌:                                         nZ |= ( fabs( atom[a1].z - atom[neigh].z ) > MIN_BOND_LENGTH );
    // INCHI✔❌:                                     }
    // INCHI✔❌:                                 }
    // INCHI✔❌:
    // INCHI✔❌:                                 chem_bonds_valence = 0;
    // INCHI✔❌:                                 for (n1 = 0; MIN_INPUT_BOND_TYPE + n1 <= 3 && MIN_INPUT_BOND_TYPE + n1 <= MAX_INPUT_BOND_TYPE; n1++)
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     chem_bonds_valence += ( MIN_INPUT_BOND_TYPE + n1 ) * num_bond_type[n1];
    // INCHI✔❌:                                 }
    // INCHI✔❌:
    // INCHI✔❌:                                 if (MIN_INPUT_BOND_TYPE <= INCHI_BOND_TYPE_ALTERN && INCHI_BOND_TYPE_ALTERN <= MAX_INPUT_BOND_TYPE &&
    // INCHI✔❌:                                     ( n2 = num_bond_type[INCHI_BOND_TYPE_ALTERN - MIN_INPUT_BOND_TYPE] ))
    // INCHI✔❌:                                 {
    // INCHI✔❌:
    // INCHI✔❌:                                     /* accept input aromatic bonds for now */
    // INCHI✔❌:
    // INCHI✔❌:                                     switch (n2)
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         case 2:
    // INCHI✔❌:                                             chem_bonds_valence += 3;  /* =A- */
    // INCHI✔❌:                                             break;
    // INCHI✔❌:
    // INCHI✔❌:                                         case 3:
    // INCHI✔❌:                                             chem_bonds_valence += 4;  /* =A< */
    // INCHI✔❌:                                             break;
    // INCHI✔❌:
    // INCHI✔❌:                                         default:
    // INCHI✔❌:                                             /*  if 1 or >= 4 aromatic bonds then replace such bonds with single bonds */
    // INCHI✔❌:                                             for (n1 = 0; n1 < AT_NUM_BONDS( atom[a1] ); n1++)
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 if (atom[a1].bond_type[n1] == INCHI_BOND_TYPE_ALTERN)
    // INCHI✔❌:                                                 {
    // INCHI✔❌:                                                     ATOM_NUMBER *p1;
    // INCHI✔❌:                                                     a2 = atom[a1].neighbor[n1];
    // INCHI✔❌:                                                     p1 = IN_NEIGH_LIST( atom[a2].neighbor, (ATOM_NUMBER) a1, AT_NUM_BONDS( atom[a2] ) );
    // INCHI✔❌:                                                     if (p1)
    // INCHI✔❌:                                                     {
    // INCHI✔❌:                                                         atom[a1].bond_type[n1] =
    // INCHI✔❌:                                                             atom[a2].bond_type[p1 - atom[a2].neighbor] = INCHI_BOND_TYPE_SINGLE;
    // INCHI✔❌:                                                     }
    // INCHI✔❌:                                                     else
    // INCHI✔❌:                                                     {
    // INCHI✔❌:                                                         *err = -2;  /*  Program error */
    // INCHI✔❌:                                                         TREAT_ERR( *err, 0, "Program error interpreting InChI aux" );
    // INCHI✔❌:                                                         num_atoms = INCHI_INP_FATAL_RET;
    // INCHI✔❌:                                                         goto bypass_end_of_INChI_plain; /*  no structure */
    // INCHI✔❌:                                                     }
    // INCHI✔❌:                                                 }
    // INCHI✔❌:                                             }
    // INCHI✔❌:
    // INCHI✔❌:                                             chem_bonds_valence += n2;
    // INCHI✔❌:                                             *err |= 32; /*  Unrecognized aromatic bond(s) replaced with single */
    // INCHI✔❌:                                             TREAT_ERR( *err, 0, "Atom has 1 or more than 3 aromatic bonds" );
    // INCHI✔❌:                                             break;
    // INCHI✔❌:                                     }
    // INCHI✔❌:                                 }
    // INCHI✔❌:                                 /********************************
    // INCHI✔❌:                                     *
    // INCHI✔❌:                                     *  Set number of hydrogen atoms
    // INCHI✔❌:                                     */
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     int num_iso_H;
    // INCHI✔❌:                                     num_iso_H = atom[a1].num_iso_H[1] + atom[a1].num_iso_H[2] + atom[a1].num_iso_H[3];
    // INCHI✔❌:                                     if (valence == ISOLATED_ATOM)
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         atom[a1].num_iso_H[0] = 0;
    // INCHI✔❌:                                     }
    // INCHI✔❌:                                     else
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         if (valence && valence >= chem_bonds_valence)
    // INCHI✔❌:                                         {
    // INCHI✔❌:                                             atom[a1].num_iso_H[0] = valence - chem_bonds_valence;
    // INCHI✔❌:                                         }
    // INCHI✔❌:                                         else
    // INCHI✔❌:                                         {
    // INCHI✔❌:                                             if (valence || bDoNotAddH)
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 atom[a1].num_iso_H[0] = 0;
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                             else
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 if (!bDoNotAddH)
    // INCHI✔❌:                                                 {
    // INCHI✔❌:                                                     atom[a1].num_iso_H[0] = -1; /* auto add H */
    // INCHI✔❌:                                                 }
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                         }
    // INCHI✔❌:                                     }
    // INCHI✔❌:                                 }
    // INCHI✔❌:                             }
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         nNumBonds /= 2;
    // INCHI✔❌:
    // INCHI✔❌:                         if (b23D && nNumBonds)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             nXYZ = nX + nY + nZ;
    // INCHI✔❌:                             b2D = ( nXYZ > 0 );
    // INCHI✔❌:                             b3D = ( nXYZ == 3 );
    // INCHI✔❌:                             *num_dimensions = b3D ? 3 : b2D ? 2 : 0;
    // INCHI✔❌:                             *num_bonds = nNumBonds;
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         /*======= 0D parities =================================*/
    // INCHI✔❌:                         if (len_stereo0D > 0 && atom_stereo0D && stereo0D)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             *stereo0D = atom_stereo0D;
    // INCHI✔❌:                             *num_stereo0D = len_stereo0D;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         else
    // INCHI✔❌:                         {
    // INCHI✔❌:                             FreeInchi_Stereo0D( &atom_stereo0D );
    // INCHI✔❌:                             *num_stereo0D = len_stereo0D = 0;
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         for (i = 0; i < len_stereo0D; i++)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             ATOM_NUMBER *p1, *p2;
    // INCHI✔❌:                             int     sb_ord_from_a1 = -1, sb_ord_from_a2 = -1, bEnd1 = 0, bEnd2 = 0;
    // INCHI✔❌:
    // INCHI✔❌:                             switch (atom_stereo0D[i].type)
    // INCHI✔❌:                             {
    // INCHI✔❌:
    // INCHI✔❌:                                 case INCHI_StereoType_Tetrahedral:
    // INCHI✔❌:                                     a1 = atom_stereo0D[i].central_atom;
    // INCHI✔❌:                                     if (atom_stereo0D[i].parity && ( AT_NUM_BONDS( atom[a1] ) == 3 || AT_NUM_BONDS( atom[a1] ) == 4 ))
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         int ii, kk = 0;
    // INCHI✔❌:                                         if (AT_NUM_BONDS( atom[a1] ) == 3)
    // INCHI✔❌:                                             atom_stereo0D[i].neighbor[kk++] = a1;
    // INCHI✔❌:                                         for (ii = 0; ii < AT_NUM_BONDS( atom[a1] ); ii++)
    // INCHI✔❌:                                             atom_stereo0D[i].neighbor[kk++] = atom[a1].neighbor[ii];
    // INCHI✔❌:                                     }
    // INCHI✔❌:
    // INCHI✔❌:                                     break;
    // INCHI✔❌:
    // INCHI✔❌:                                 case INCHI_StereoType_DoubleBond:
    // INCHI✔❌: #define MAX_CHAIN_LEN 20
    // INCHI✔❌:                                     a1 = atom_stereo0D[i].neighbor[1];
    // INCHI✔❌:                                     a2 = atom_stereo0D[i].neighbor[2];
    // INCHI✔❌:                                     p1 = IN_NEIGH_LIST( atom[a1].neighbor, (ATOM_NUMBER) a2, AT_NUM_BONDS( atom[a1] ) );
    // INCHI✔❌:                                     p2 = IN_NEIGH_LIST( atom[a2].neighbor, (ATOM_NUMBER) a1, AT_NUM_BONDS( atom[a2] ) );
    // INCHI✔❌:                                     if (!p1 || !p2)
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         atom_stereo0D[i].type = INCHI_StereoType_None;
    // INCHI✔❌:                                         atom_stereo0D[i].central_atom = NO_ATOM;
    // INCHI✔❌:                                         atom_stereo0D[i].neighbor[0] =
    // INCHI✔❌:                                             atom_stereo0D[i].neighbor[3] = -1;
    // INCHI✔❌:                                         *err |= 64; /* Error in cumulene stereo */
    // INCHI✔❌:                                         TREAT_ERR( *err, 0, "0D stereobond not recognized" );
    // INCHI✔❌:                                         break;
    // INCHI✔❌:                                     }
    // INCHI✔❌:
    // INCHI✔❌:                                     /* streobond, allene, or cumulene */
    // INCHI✔❌:
    // INCHI✔❌:                                     sb_ord_from_a1 = p1 - atom[a1].neighbor;
    // INCHI✔❌:                                     sb_ord_from_a2 = p2 - atom[a2].neighbor;
    // INCHI✔❌:
    // INCHI✔❌:                                     if (AT_NUM_BONDS( atom[a1] ) == 2 &&
    // INCHI✔❌:                                             atom[a1].bond_type[0] + atom[a1].bond_type[1] == 2 * INCHI_BOND_TYPE_DOUBLE &&
    // INCHI✔❌:                                             0 == inchi_NUMH2( atom, a1 ) &&
    // INCHI✔❌:                                             ( AT_NUM_BONDS( atom[a2] ) != 2 ||
    // INCHI✔❌:                                                 atom[a2].bond_type[0] + atom[a2].bond_type[1] != 2 * INCHI_BOND_TYPE_DOUBLE ))
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         bEnd2 = 1; /* a2 is the end-atom, a1 is middle atom */
    // INCHI✔❌:                                     }
    // INCHI✔❌:
    // INCHI✔❌:                                     if (AT_NUM_BONDS( atom[a2] ) == 2 &&
    // INCHI✔❌:                                             atom[a2].bond_type[0] + atom[a2].bond_type[1] == 2 * INCHI_BOND_TYPE_DOUBLE &&
    // INCHI✔❌:                                             0 == inchi_NUMH2( atom, a2 ) &&
    // INCHI✔❌:                                             ( AT_NUM_BONDS( atom[a1] ) != 2 ||
    // INCHI✔❌:                                                 atom[a1].bond_type[0] + atom[a1].bond_type[1] != 2 * INCHI_BOND_TYPE_DOUBLE ))
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         bEnd1 = 1; /* a1 is the end-atom, a2 is middle atom */
    // INCHI✔❌:                                     }
    // INCHI✔❌:
    // INCHI✔❌:                                     if (bEnd2 + bEnd1 == 1)
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         /* allene or cumulene */
    // INCHI✔❌:
    // INCHI✔❌:                                         ATOM_NUMBER  chain[MAX_CHAIN_LEN + 1], prev, cur, next;
    // INCHI✔❌:
    // INCHI✔❌:                                         if (bEnd2 && !bEnd1)
    // INCHI✔❌:                                         {
    // INCHI✔❌:                                             cur = a1;
    // INCHI✔❌:                                             a1 = a2;
    // INCHI✔❌:                                             a2 = cur;
    // INCHI✔❌:                                             sb_ord_from_a1 = sb_ord_from_a2;
    // INCHI✔❌:                                         }
    // INCHI✔❌:
    // INCHI✔❌:                                         sb_ord_from_a2 = -1;
    // INCHI✔❌:                                         cur = a1;
    // INCHI✔❌:                                         next = a2;
    // INCHI✔❌:                                         len = 0;
    // INCHI✔❌:                                         chain[len++] = cur;
    // INCHI✔❌:                                         chain[len++] = next;
    // INCHI✔❌:
    // INCHI✔❌:                                         while (len < MAX_CHAIN_LEN)
    // INCHI✔❌:                                         {
    // INCHI✔❌:                                             /* arbitrary very high upper limit to prevent infinite loop */
    // INCHI✔❌:
    // INCHI✔❌:                                             prev = cur;
    // INCHI✔❌:                                             cur = next;
    // INCHI✔❌:                                                 /* follow double bond path && avoid going back */
    // INCHI✔❌:                                             if (AT_NUM_BONDS( atom[cur] ) == 2 &&
    // INCHI✔❌:                                                     atom[cur].bond_type[0] + atom[cur].bond_type[1] == 2 * INCHI_BOND_TYPE_DOUBLE &&
    // INCHI✔❌:                                                     0 == inchi_NUMH2( atom, cur ))
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 next = atom[cur].neighbor[atom[cur].neighbor[0] == prev];
    // INCHI✔❌:                                                 chain[len++] = next;
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                             else
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 break;
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                         }
    // INCHI✔❌:                                         if (len > 2 &&
    // INCHI✔❌:                                             ( p2 = IN_NEIGH_LIST( atom[cur].neighbor, (ATOM_NUMBER) prev, AT_NUM_BONDS( atom[cur] ) ) ))
    // INCHI✔❌:                                         {
    // INCHI✔❌:                                             sb_ord_from_a2 = p2 - atom[cur].neighbor;
    // INCHI✔❌:                                             a2 = cur;
    // INCHI✔❌:                                             /* by design we need to pick up the first non-stereo-bond-neighbor as "sn"-atom */
    // INCHI✔❌:                                             atom_stereo0D[i].neighbor[0] = atom[a1].neighbor[sb_ord_from_a1 == 0];
    // INCHI✔❌:                                             atom_stereo0D[i].neighbor[1] = a1;
    // INCHI✔❌:                                             atom_stereo0D[i].neighbor[2] = a2;
    // INCHI✔❌:                                             atom_stereo0D[i].neighbor[3] = atom[a2].neighbor[sb_ord_from_a2 == 0];
    // INCHI✔❌:
    // INCHI✔❌:                                             if (len % 2)
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 atom_stereo0D[i].central_atom = chain[len / 2];
    // INCHI✔❌:                                                 atom_stereo0D[i].type = INCHI_StereoType_Allene;
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                             else
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 atom_stereo0D[i].central_atom = NO_ATOM;
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                         }
    // INCHI✔❌:                                         else
    // INCHI✔❌:                                         {
    // INCHI✔❌:                                             /* error */
    // INCHI✔❌:                                             atom_stereo0D[i].type = INCHI_StereoType_None;
    // INCHI✔❌:                                             atom_stereo0D[i].central_atom = NO_ATOM;
    // INCHI✔❌:                                             atom_stereo0D[i].neighbor[0] =
    // INCHI✔❌:                                                 atom_stereo0D[i].neighbor[3] = -1;
    // INCHI✔❌:                                             *err |= 64; /* Error in cumulene stereo */
    // INCHI✔❌:                                             TREAT_ERR( *err, 0, "Cumulene stereo not recognized (0D)" );
    // INCHI✔❌:                                         }
    // INCHI✔❌: #undef MAX_CHAIN_LEN
    // INCHI✔❌:                                     }
    // INCHI✔❌:                                     else
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         /****** a normal possibly stereogenic bond -- not an allene or cumulene *******/
    // INCHI✔❌:                                         /* by design we need to pick up the first non-stereo-bond-neighbor as "sn"-atom */
    // INCHI✔❌:                                         sb_ord_from_a1 = p1 - atom[a1].neighbor;
    // INCHI✔❌:                                         sb_ord_from_a2 = p2 - atom[a2].neighbor;
    // INCHI✔❌:                                         atom_stereo0D[i].neighbor[0] = atom[a1].neighbor[p1 == atom[a1].neighbor];
    // INCHI✔❌:                                         atom_stereo0D[i].neighbor[3] = atom[a2].neighbor[p2 == atom[a2].neighbor];
    // INCHI✔❌:                                         atom_stereo0D[i].central_atom = NO_ATOM;
    // INCHI✔❌:                                     }
    // INCHI✔❌:
    // INCHI✔❌:                                     if (atom_stereo0D[i].type != INCHI_StereoType_None &&
    // INCHI✔❌:                                             sb_ord_from_a1 >= 0 && sb_ord_from_a2 >= 0 &&
    // INCHI✔❌:                                             ATOM_PARITY_WELL_DEF( SB_PARITY_2( atom_stereo0D[i].parity ) ))
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         /* Detected well-defined disconnected stereo
    // INCHI✔❌:                                             * locate first non-metal neighbors */
    // INCHI✔❌:
    // INCHI✔❌:                                         int    a, n, j, /* k,*/ sb_ord, cur_neigh, min_neigh;
    // INCHI✔❌:
    // INCHI✔❌:                                         for (k = 0; k < 2; k++)
    // INCHI✔❌:                                         {
    // INCHI✔❌:                                             a = k ? atom_stereo0D[i].neighbor[2] : atom_stereo0D[i].neighbor[1];
    // INCHI✔❌:                                             sb_ord = k ? sb_ord_from_a2 : sb_ord_from_a1;
    // INCHI✔❌:                                             min_neigh = num_atoms;
    // INCHI✔❌:                                             for (n = j = 0; j < AT_NUM_BONDS( atom[a] ); j++)
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 cur_neigh = atom[a].neighbor[j];
    // INCHI✔❌:                                                 if (j != sb_ord && !IS_METAL_ATOM( atom, cur_neigh ))
    // INCHI✔❌:                                                 {
    // INCHI✔❌:                                                     min_neigh = inchi_min( cur_neigh, min_neigh );
    // INCHI✔❌:                                                 }
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                             if (min_neigh < num_atoms)
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 atom_stereo0D[i].neighbor[k ? 3 : 0] = min_neigh;
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                             else
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 TREAT_ERR( *err, 0, "Cannot find non-metal stereobond neighor (0D)" );
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                         }
    // INCHI✔❌:                                     }
    // INCHI✔❌:
    // INCHI✔❌:                                     break;
    // INCHI✔❌:                             }
    // INCHI✔❌:                         }
    // INCHI✔❌:                         /* end of 0D parities extraction */
    // INCHI✔❌: /*exit_cycle:;*/
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     if (pInpAtomFlags)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         /* save chirality flag */
    // INCHI✔❌:                         *pInpAtomFlags |= InpAtomFlags;
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else if (atom)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     inchi_free( atom );
    // INCHI✔❌:                     atom = NULL;
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 if (pszCoord)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     inchi_free( pszCoord );
    // INCHI✔❌:                     pszCoord = NULL;
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 goto bypass_end_of_INChI_plain;
    // INCHI✔❌:                 /*return num_atoms;*/
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         if (atom_stereo0D)
    // INCHI✔❌:         {
    // INCHI✔❌:             FreeInchi_Stereo0D( &atom_stereo0D );
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌:         /* end of structure reading cycle */
    // INCHI✔❌:
    // INCHI✔❌:         if (res <= 0)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (*err == INCHI_INP_ERROR_ERR)
    // INCHI✔❌:                 return num_atoms;
    // INCHI✔❌:
    // INCHI✔❌:             *err = INCHI_INP_EOF_ERR;
    // INCHI✔❌:             return INCHI_INP_EOF_RET; /* no more data */
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌:     bypass_end_of_INChI_plain:
    // INCHI✔❌:
    // INCHI✔❌:             /* cleanup */
    // INCHI✔❌:         if (num_atoms == INCHI_INP_ERROR_RET || num_atoms == INCHI_INP_FATAL_RET)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (atom_stereo0D) /* djb-rwth: fixing coverity ID #500395 */
    // INCHI✔❌:             {
    // INCHI✔❌:                 if (stereo0D && *stereo0D == atom_stereo0D)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     *stereo0D = NULL;
    // INCHI✔❌:                     *num_stereo0D = 0;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 FreeInchi_Stereo0D(&atom_stereo0D);
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             if (atom) /* djb-rwth: fixing coverity ID #499615 */
    // INCHI✔❌:             {
    // INCHI✔❌:                 inchi_free(atom);
    // INCHI✔❌:                 atom = NULL;
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             if (pszCoord) /* djb-rwth: fixing coverity ID #500397 */
    // INCHI✔❌:             {
    // INCHI✔❌:                 inchi_free(pszCoord);
    // INCHI✔❌:                 pszCoord = NULL;
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         while (bTooLongLine &&
    // INCHI✔❌:                 0 < inchi_ios_getsTab1(szLine_i2ia, sizeof(szLine_i2ia) - 1, inp_file, &bTooLongLine ))
    // INCHI✔❌:         {
    // INCHI✔❌:             ;
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌:         /* cleanup
    // INCHI✔❌:         if (at && !*at)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (atom)
    // INCHI✔❌:             {
    // INCHI✔❌:                 inchi_free( atom );
    // INCHI✔❌:
    // INCHI✔❌:             }
    // INCHI✔❌:             if (pszCoord)
    // INCHI✔❌:             {
    // INCHI✔❌:                 inchi_free( pszCoord );
    // INCHI✔❌:                 pszCoord = NULL;
    // INCHI✔❌:             }
    // INCHI✔❌:         }*/
    // INCHI✔❌:
    // INCHI✔❌:         return num_atoms;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /***********************************************************/
    // INCHI✔❌:     /*   extract reversibility info from xml text INChI format */
    // INCHI✔❌:     /*                                                         */
    // INCHI✔❌:     /*   OBSOLETE CODE because InChI output in XML             */
    // INCHI✔❌:     /*      does not exist anymore. Unsupported.               */
    // INCHI✔❌:     /*                                                         */
    // INCHI✔❌:     /***********************************************************/
    // INCHI✔❌:
    // INCHI✔❌:     if (nInputType == INPUT_INCHI_XML)
    // INCHI✔❌:     {
    // INCHI✔❌:
    // INCHI✔❌:         /* xml tags */
    // INCHI✔❌:
    // INCHI✔❌:         static const char sStructHdrXml[] = "<structure";
    // INCHI✔❌:         static const char sStructHdrXmlEnd[] = "</structure";
    // INCHI✔❌:         static const char sStructHdrXmlNumber[] = "number=\"";
    // INCHI✔❌:         static const char sStructHdrXmlIdName[] = "id.name=\"";
    // INCHI✔❌:         static const char sStructHdrXmlIdValue[] = "id.value=\"";
    // INCHI✔❌:         static const char sStructMsgXmlErr[] = "<message type=\"error (no InChI)\" value=\"";
    // INCHI✔❌:         static const char sStructMsgXmlErrFatal[] = "<message type=\"fatal (aborted)\" value=\"";
    // INCHI✔❌:         static const char sStructRevXmlRevHdr[] = "<reversibility>";
    // INCHI✔❌:         static const char sStructRevXmlRevAt[] = "<atoms>";
    // INCHI✔❌:         static const char sStructRevXmlRevAtEnd[] = "</atoms>";
    // INCHI✔❌:         static const char sStructRevXmlRevBn[] = "<bonds>";
    // INCHI✔❌:         static const char sStructRevXmlRevBnEnd[] = "</bonds>";
    // INCHI✔❌:         static const char sStructRevXmlRevXYZ[] = "<xyz>";
    // INCHI✔❌:         static const char sStructRevXmlRevXYZEnd[] = "</xyz>";
    // INCHI✔❌:         static const char sStructAuxXml[] = "<identifier.auxiliary-info";
    // INCHI✔❌:         static const char sStructAuxXmlEnd[] = "</identifier.auxiliary-info";
    // INCHI✔❌:         int         bInTheAuxInfo = 0;
    // INCHI✔❌:
    // INCHI✔❌:         while (0 < ( res = inchi_ios_gets(szLine_i2ia, sizeof(szLine_i2ia) - 1, inp_file, &bTooLongLine ) ))
    // INCHI✔❌:         {
    // INCHI✔❌:             /********************* find and interpret structure header ************/
    // INCHI✔❌:
    // INCHI✔❌:             if (!memcmp(szLine_i2ia, sStructHdrXml, sizeof( sStructHdrXml ) - 1 ))
    // INCHI✔❌:             {
    // INCHI✔❌:                 num_struct = 1;
    // INCHI✔❌:                 p = szLine_i2ia + sizeof( sStructHdrXml ) - 1;
    // INCHI✔❌:                 longID = 0;
    // INCHI✔❌:                 num_atoms = 0;
    // INCHI✔❌:                 /* structure number */
    // INCHI✔❌:                 if ((q = strstr( p, sStructHdrXmlNumber ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:                 {
    // INCHI✔❌:                     p = q + sizeof( sStructHdrXmlNumber ) - 1;
    // INCHI✔❌:                     longID = strtol( p, &q, 10 );
    // INCHI✔❌:                     if (q && *q == '\"')
    // INCHI✔❌:                     {
    // INCHI✔❌:                         p = q + 1;
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:                 if (pSdfLabel)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     pSdfLabel[0] = '\0';
    // INCHI✔❌:                 }
    // INCHI✔❌:                 if (pSdfValue)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     pSdfValue[0] = '\0';
    // INCHI✔❌:                 }
    // INCHI✔❌:                 /* pSdfLabel */
    // INCHI✔❌:                 if ((q = strstr( p, sStructHdrXmlIdName ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:                 {
    // INCHI✔❌:                     p = q + sizeof( sStructHdrXmlIdName ) - 1;
    // INCHI✔❌:                     q = strchr( p, '\"' );
    // INCHI✔❌:                     if (q)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         len = inchi_min( q - p + 1, MAX_SDF_HEADER - 1 );
    // INCHI✔❌:                         if (pSdfLabel)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             mystrncpy( pSdfLabel, p, len );
    // INCHI✔❌:                         }
    // INCHI✔❌:                         p = q + 1;
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:                 /* pSdfValue */
    // INCHI✔❌:                 if ((q = strstr( p, sStructHdrXmlIdValue ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:                 {
    // INCHI✔❌:                     p = q + sizeof( sStructHdrXmlIdValue ) - 1;
    // INCHI✔❌:                     q = strchr( p, '\"' );
    // INCHI✔❌:                     if (q)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         len = inchi_min( q - p + 1, MAX_SDF_VALUE - 1 );
    // INCHI✔❌:                         if (pSdfValue)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             mystrncpy( pSdfValue, p, len );
    // INCHI✔❌:                         }
    // INCHI✔❌:                         p = q + 1; /* djb-rwth: ignoring LLVM warning: value used */
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:                 if (Id)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     *Id = longID;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 bHeaderRead = 1;
    // INCHI✔❌:                 bErrorMsg = bRestoreInfo = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if ( (bHeaderRead && !memcmp(szLine_i2ia, sStructMsgXmlErr, sizeof(sStructMsgXmlErr) - 1) ) ||
    // INCHI✔❌:                     (bHeaderRead && !memcmp(szLine_i2ia, sStructMsgXmlErrFatal, sizeof(sStructMsgXmlErrFatal) - 1) ) ) /* djb-rwth: fixed incorrectly written operators */
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* djb-rwth: remaining block from the above condition */
    // INCHI✔❌:                 if (bHeaderRead && !memcmp(szLine_i2ia, sStructMsgXmlErr, sizeof(sStructMsgXmlErr) - 1))
    // INCHI✔❌:                 {
    // INCHI✔❌:                     bFatal = 0;
    // INCHI✔❌:                     len = sizeof(sStructMsgXmlErr) - 1;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 if (bHeaderRead && !memcmp(szLine_i2ia, sStructMsgXmlErrFatal, sizeof(sStructMsgXmlErrFatal) - 1))
    // INCHI✔❌:                 {
    // INCHI✔❌:                     bFatal = 1;
    // INCHI✔❌:                     len = sizeof(sStructMsgXmlErrFatal) - 1;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 p = szLine_i2ia + len;
    // INCHI✔❌:                 q = strchr( p, '\"' );
    // INCHI✔❌:                 if (q && !bFindNext)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     int c;
    // INCHI✔❌:                     bErrorMsg = 1;
    // INCHI✔❌:                     pStrErr[0] = '\0';
    // INCHI✔❌:                     c = *q;
    // INCHI✔❌:                     *q = '\0';
    // INCHI✔❌:                     TREAT_ERR( *err, 0, p );
    // INCHI✔❌:                     *q = c;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 *err = bFatal ? INCHI_INP_FATAL_ERR : INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                 num_atoms = bFatal ? INCHI_INP_FATAL_RET : INCHI_INP_ERROR_RET;
    // INCHI✔❌:                 goto bypass_end_of_INChI;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (bHeaderRead && !memcmp(szLine_i2ia, sStructAuxXml, sizeof( sStructAuxXml ) - 1 ))
    // INCHI✔❌:             {
    // INCHI✔❌:                 bInTheAuxInfo = 1;
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (bHeaderRead && !memcmp(szLine_i2ia, sStructAuxXmlEnd, sizeof( sStructAuxXmlEnd ) - 1 ))
    // INCHI✔❌:             {
    // INCHI✔❌:                 *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                 num_atoms = INCHI_INP_ERROR_RET;
    // INCHI✔❌:                 TREAT_ERR( *err, 0, "Missing reversibility info" );
    // INCHI✔❌:                 goto bypass_end_of_INChI; /* reversibility info not found */
    // INCHI✔❌:             }
    // INCHI✔❌:             else if (bHeaderRead && bInTheAuxInfo && !memcmp(szLine_i2ia, sStructRevXmlRevHdr, sizeof( sStructRevXmlRevHdr ) - 1 ))
    // INCHI✔❌:             {
    // INCHI✔❌:                 /***********************  atoms xml ***************************/
    // INCHI✔❌:                 num_struct = 1;
    // INCHI✔❌:                 res = inchi_ios_gets(szLine_i2ia, sizeof(szLine_i2ia) - 1, inp_file, &bTooLongLine );
    // INCHI✔❌:                 if (res <= 0)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     num_atoms = INCHI_INP_EOF_RET; /* no data, probably end of file */
    // INCHI✔❌:                     *err = INCHI_INP_EOF_ERR;
    // INCHI✔❌:                     goto bypass_end_of_INChI;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 if (memcmp(szLine_i2ia, sStructRevXmlRevAt, sizeof( sStructRevXmlRevAt ) - 1 ))
    // INCHI✔❌:                 {
    // INCHI✔❌:                     bHeaderRead = 0; /* invalid reversibility info; look for another header */
    // INCHI✔❌:                     continue;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 /* read (the head of) the atoms line */
    // INCHI✔❌:                 res = inchi_ios_gets(szLine_i2ia, sizeof(szLine_i2ia) - 1, inp_file, &bTooLongLine );
    // INCHI✔❌:                 if (res <= 0)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     num_atoms = INCHI_INP_EOF_RET; /* no data */
    // INCHI✔❌:                     *err = INCHI_INP_EOF_ERR;
    // INCHI✔❌:                     goto bypass_end_of_INChI;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 p = szLine_i2ia;
    // INCHI✔❌:                 num_atoms = strtol( p, &q, 10 );
    // INCHI✔❌:                 if (!num_atoms || !q || !*q)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     num_atoms = INCHI_INP_EOF_RET; /* no atom data */
    // INCHI✔❌:                     *err = INCHI_INP_EOF_ERR;
    // INCHI✔❌:                     goto bypass_end_of_INChI;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 p = q;
    // INCHI✔❌:                 /* Molfile chirality flag */
    // INCHI✔❌:                 switch (*p)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     case 'c':
    // INCHI✔❌:                         InpAtomFlags |= FLAG_INP_AT_CHIRAL;
    // INCHI✔❌:                         p++;
    // INCHI✔❌:                         break;
    // INCHI✔❌:                     case 'n':
    // INCHI✔❌:                         InpAtomFlags |= FLAG_INP_AT_NONCHIRAL;
    // INCHI✔❌:                         p++;
    // INCHI✔❌:                         break;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 if (at && *at)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     if (num_atoms > max_num_at)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         inchi_free( *at );
    // INCHI✔❌:                         *at = NULL;
    // INCHI✔❌:                     }
    // INCHI✔❌:                     else
    // INCHI✔❌:                     {
    // INCHI✔❌:                         memset( *at, 0, max_num_at * sizeof( **at ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:                         atom = *at;
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:                 if (!at || !*at)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     atom = CreateInchiAtom( num_atoms + 1 );
    // INCHI✔❌:                     if (!atom)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         num_atoms = INCHI_INP_FATAL_RET; /* fatal error: cannot allocate */
    // INCHI✔❌:                         *err = INCHI_INP_FATAL_ERR;
    // INCHI✔❌:                         TREAT_ERR( *err, 0, "Out of RAM" );
    // INCHI✔❌:                         goto bypass_end_of_INChI;
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:                 if (stereo0D && *stereo0D)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     if (num_atoms > max_len_stereo0D)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         FreeInchi_Stereo0D( stereo0D );
    // INCHI✔❌:                     }
    // INCHI✔❌:                     else
    // INCHI✔❌:                     {
    // INCHI✔❌:                         memset( *stereo0D, 0, max_len_stereo0D * sizeof( **stereo0D ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:                         atom_stereo0D = *stereo0D;
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:                 if (!stereo0D || !*stereo0D)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     max_len_stereo0D = num_atoms + 1;
    // INCHI✔❌:                     atom_stereo0D = CreateInchi_Stereo0D( max_len_stereo0D );
    // INCHI✔❌:                     if (!atom_stereo0D)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         num_atoms = INCHI_INP_FATAL_RET; /* fatal error: cannot allocate */
    // INCHI✔❌:                         *err = INCHI_INP_FATAL_ERR;
    // INCHI✔❌:                         TREAT_ERR( *err, 0, "Out of RAM" );
    // INCHI✔❌:                         /* djb-rwth: avoiding memory leak */
    // INCHI✔❌:                         inchi_free(atom);
    // INCHI✔❌:                         atom = NULL;
    // INCHI✔❌:                         goto bypass_end_of_INChI;
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 i = 0;
    // INCHI✔❌:                 bItemIsOver = 0;
    // INCHI✔❌:                 res2 = bTooLongLine2 = -1; /* djb-rwth: ignoring LLVM warning: value used */
    // INCHI✔❌:
    // INCHI✔❌:                 /* read all atoms xml */
    // INCHI✔❌:                 while (i < num_atoms)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     pos = p - szLine_i2ia;
    // INCHI✔❌:                     if (!bItemIsOver && ( int )sizeof(szLine_i2ia) - res + pos > ( int )sizeof(szNextLine_i2ia))
    // INCHI✔❌:                     {
    // INCHI✔❌:                         /* load next line if possible */
    // INCHI✔❌:                         res2 = inchi_ios_gets( szNextLine_i2ia, sizeof(szNextLine_i2ia) - 1, inp_file, &bTooLongLine2 );
    // INCHI✔❌:                         if (res2 > 0 && memcmp( szNextLine_i2ia, sStructRevXmlRevAtEnd, sizeof( sStructRevXmlRevAtEnd ) - 1 ))
    // INCHI✔❌:                         {
    // INCHI✔❌:                             if (pos)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 res -= pos;  /* number of chars left to process in szLine */
    // INCHI✔❌:                                 memmove(szLine_i2ia, p, res * sizeof(szLine_i2ia[0] ) ); /* move them to the start of the line */
    // INCHI✔❌:                             }
    // INCHI✔❌:                             memcpy(szLine_i2ia + res, szNextLine_i2ia, ( (long long)res2 + 1 ) * sizeof(szNextLine_i2ia[0] ) ); /* djb-rwth: cast operator added */
    // INCHI✔❌:                             res += res2;
    // INCHI✔❌:                             szLine_i2ia[res] = '\0';
    // INCHI✔❌:                             bTooLongLine = bTooLongLine2;
    // INCHI✔❌:                             p = szLine_i2ia;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         else
    // INCHI✔❌:                         {
    // INCHI✔❌:                             bItemIsOver = 1;
    // INCHI✔❌:                         }
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     /* element, first char */
    // INCHI✔❌:                     if (!isalpha( UCINT *p ) || !isupper( UCINT *p ) || i >= num_atoms)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         bHeaderRead = 0; /* wrong atom data */ /* djb-rwth: ignoring LLVM warning: value used */
    // INCHI✔❌:                         num_atoms = INCHI_INP_ERROR_RET; /* was 0, error */
    // INCHI✔❌:                         *err = INCHI_INP_ERROR_ERR;     /* 40 */
    // INCHI✔❌:                         TREAT_ERR( *err, 0, "Wrong atoms data" );
    // INCHI✔❌:                         /* djb-rwth: avoiding memory leak */
    // INCHI✔❌:                         inchi_free(atom);
    // INCHI✔❌:                         atom = NULL;
    // INCHI✔❌:                         goto bypass_end_of_INChI;
    // INCHI✔❌:                     }
    // INCHI✔❌:                     atom[i].elname[0] = *p++;
    // INCHI✔❌:                     /* element, second char */
    // INCHI✔❌:                     if (isalpha( UCINT *p ) && islower( UCINT *p ))
    // INCHI✔❌:                     {
    // INCHI✔❌:                         atom[i].elname[1] = *p++;
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     /* bonds' valence */
    // INCHI✔❌:                     if (isdigit( UCINT *p ))
    // INCHI✔❌:                     {
    // INCHI✔❌:                         AT_BONDS_VAL( atom, i ) = (char) strtol( p, &q, 10 );
    // INCHI✔❌:                         if (!AT_BONDS_VAL( atom, i ))
    // INCHI✔❌:                             AT_BONDS_VAL( atom, i ) = ISOLATED_ATOM; /* same convention as in MOLfile, found zero bonds valence */
    // INCHI✔❌:                         p = q;
    // INCHI✔❌:                     }
    // INCHI✔❌:                     /* charge */
    // INCHI✔❌:                     atom[i].charge = ( *p == '+' ) ? 1 : ( *p == '-' ) ? -1 : 0;
    // INCHI✔❌:                     if (atom[i].charge)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         p++;
    // INCHI✔❌:                         if (isdigit( UCINT *p ))
    // INCHI✔❌:                         {
    // INCHI✔❌:                             atom[i].charge *= (S_CHAR) ( strtol( p, &q, 10 ) & CHAR_MASK );
    // INCHI✔❌:                             p = q;
    // INCHI✔❌:                         }
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     /* radical */
    // INCHI✔❌:                     if (*p == '.')
    // INCHI✔❌:                     {
    // INCHI✔❌:                         p++;
    // INCHI✔❌:                         if (isdigit( UCINT *p ))
    // INCHI✔❌:                         {
    // INCHI✔❌:                             atom[i].radical = (S_CHAR) strtol( p, &q, 10 );
    // INCHI✔❌:                             p = q;
    // INCHI✔❌:                         }
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     /* isotopic mass */
    // INCHI✔❌:                     if (*p == 'i')
    // INCHI✔❌:                     {
    // INCHI✔❌:                         p++;
    // INCHI✔❌:                         if (isdigit( UCINT *p ))
    // INCHI✔❌:                         {
    // INCHI✔❌:                             int mw = strtol( p, &q, 10 );
    // INCHI✔❌:                             p = q;
    // INCHI✔❌:
    // INCHI✔❌:                             atom[i].isotopic_mass = mw;
    // INCHI✔❌:                         }
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     /* parity */
    // INCHI✔❌:                     switch (*p)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         case 'o':
    // INCHI✔❌:                             parity = INCHI_PARITY_ODD;
    // INCHI✔❌:                             p++;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         case 'e':
    // INCHI✔❌:                             parity = INCHI_PARITY_EVEN;
    // INCHI✔❌:                             p++;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         case 'u':
    // INCHI✔❌:                             parity = INCHI_PARITY_UNKNOWN;
    // INCHI✔❌:                             p++;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         case '?':
    // INCHI✔❌:                             parity = INCHI_PARITY_UNDEFINED;
    // INCHI✔❌:                             p++;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         default:
    // INCHI✔❌:                             parity = 0;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                     }
    // INCHI✔❌:                     if (parity)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         atom_stereo0D[len_stereo0D].central_atom = i;
    // INCHI✔❌:                         atom_stereo0D[len_stereo0D].parity = parity;
    // INCHI✔❌:                         atom_stereo0D[len_stereo0D].type = INCHI_StereoType_Tetrahedral;
    // INCHI✔❌:                         len_stereo0D++;
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     /* isotopic h, d, t */
    // INCHI✔❌:                     for (k = 0; k < NUM_H_ISOTOPES; k++)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         if (*p == szIsoH[k])
    // INCHI✔❌:                         {
    // INCHI✔❌:                             NUM_ISO_Hk( atom, i, k ) = 1;
    // INCHI✔❌:                             p++;
    // INCHI✔❌:                             if (isdigit( UCINT *p ))
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 NUM_ISO_Hk( atom, i, k ) = (char) strtol( p, &q, 10 );
    // INCHI✔❌:                                 p = q;
    // INCHI✔❌:                             }
    // INCHI✔❌:                         }
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     i++;
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 if (!bItemIsOver || p - szLine_i2ia != res || i != num_atoms)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     num_atoms = INCHI_INP_ERROR_RET; /* error */
    // INCHI✔❌:                     *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                     TREAT_ERR( *err, 0, "Wrong number of atoms" );
    // INCHI✔❌:                     /* djb-rwth: avoiding memory leak */
    // INCHI✔❌:                     inchi_free(atom);
    // INCHI✔❌:                     atom = NULL;
    // INCHI✔❌:                     goto bypass_end_of_INChI;
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 /********************** bonds xml ****************************/
    // INCHI✔❌:
    // INCHI✔❌:                 res = inchi_ios_gets( szLine_i2ia, sizeof(szLine_i2ia) - 1, inp_file, &bTooLongLine );
    // INCHI✔❌:                 if (res <= 0)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     num_atoms = 0; /* no data */
    // INCHI✔❌:                     goto bypass_end_of_INChI;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 if (memcmp( szLine_i2ia, sStructRevXmlRevBn, sizeof( sStructRevXmlRevBn ) - 1 ))
    // INCHI✔❌:                 {
    // INCHI✔❌:                     bHeaderRead = 0; /* invalid reversibility info; look for another header */
    // INCHI✔❌:                     continue;
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 /* read (the head of) the xml bonds line */
    // INCHI✔❌:
    // INCHI✔❌:                 res = inchi_ios_gets( szLine_i2ia, sizeof(szLine_i2ia) - 1, inp_file, &bTooLongLine );
    // INCHI✔❌:
    // INCHI✔❌:                 if (res <= 0)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     num_atoms = INCHI_INP_ERROR_RET; /* was 0; error: no data -- eof? */
    // INCHI✔❌:                     *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                     goto bypass_end_of_INChI;
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 i = 1;
    // INCHI✔❌:                 bItemIsOver = 0;
    // INCHI✔❌:                 res2 = bTooLongLine2 = -1; /* djb-rwth: ignoring LLVM warning: value used */
    // INCHI✔❌:                 p = szLine_i2ia;
    // INCHI✔❌:
    // INCHI✔❌:                 if (!memcmp( szLine_i2ia, sStructRevXmlRevBnEnd, sizeof( sStructRevXmlRevBnEnd ) - 1 ))
    // INCHI✔❌:                 {
    // INCHI✔❌:                     /* empty bonds section */
    // INCHI✔❌:                     res = 0;
    // INCHI✔❌:                     bItemIsOver = 1;
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 /* read all bonds (xml), starting from atom 1 (not 0) */
    // INCHI✔❌:
    // INCHI✔❌:                 while (i < num_atoms)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     pos = p - szLine_i2ia;
    // INCHI✔❌:                     if (!bItemIsOver &&
    // INCHI✔❌:                         ( int )sizeof(szLine_i2ia) - res + pos > ( int )sizeof(szNextLine_i2ia))
    // INCHI✔❌:                     {
    // INCHI✔❌:                         /* load next line if possible */
    // INCHI✔❌:
    // INCHI✔❌:                         res2 = inchi_ios_gets( szNextLine_i2ia, sizeof(szNextLine_i2ia) - 1, inp_file, &bTooLongLine2 );
    // INCHI✔❌:
    // INCHI✔❌:                         if (res2 > 0 && memcmp( szNextLine_i2ia, sStructRevXmlRevBnEnd, sizeof( sStructRevXmlRevBnEnd ) - 1 ))
    // INCHI✔❌:                         {
    // INCHI✔❌:                             if (pos)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 res -= pos;  /* number of chars left to process in szLine */
    // INCHI✔❌:                                 memmove(szLine_i2ia, p, res * sizeof(szLine_i2ia[0] ) ); /* move them to the start of the line */
    // INCHI✔❌:                             }
    // INCHI✔❌:                             memcpy( szLine_i2ia + res, szNextLine_i2ia, ( (long long)res2 + 1 ) * sizeof(szNextLine_i2ia[0] ) ); /* djb-rwth: cast operator added */
    // INCHI✔❌:                             res += res2;
    // INCHI✔❌:                             szLine_i2ia[res] = '\0';
    // INCHI✔❌:                             bTooLongLine = bTooLongLine2;
    // INCHI✔❌:                             p = szLine_i2ia;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         else
    // INCHI✔❌:                             bItemIsOver = 1;
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     if (i >= num_atoms)
    // INCHI✔❌:                         break;
    // INCHI✔❌:
    // INCHI✔❌:                     /* bond, first char */
    // INCHI✔❌:
    // INCHI✔❌:                     if (*p == ';')
    // INCHI✔❌:                     {
    // INCHI✔❌:                         p++;
    // INCHI✔❌:                         i++;
    // INCHI✔❌:                         continue;
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     if (!isalpha( UCINT *p ))
    // INCHI✔❌:                     {
    // INCHI✔❌:                         num_atoms = INCHI_INP_ERROR_RET; /* error in input data */
    // INCHI✔❌:                         *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                         TREAT_ERR( *err, 0, "Wrong bonds data" );
    // INCHI✔❌:                         goto bypass_end_of_INChI;
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     bond_char = *p++;
    // INCHI✔❌:
    // INCHI✔❌:                     /* bond parity */
    // INCHI✔❌:
    // INCHI✔❌:                     switch (*p)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         case '-':
    // INCHI✔❌:                             bond_parity = INCHI_PARITY_ODD;
    // INCHI✔❌:                             p++;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         case '+':
    // INCHI✔❌:                             bond_parity = INCHI_PARITY_EVEN;
    // INCHI✔❌:                             p++;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         case 'u':
    // INCHI✔❌:                             bond_parity = INCHI_PARITY_UNKNOWN;
    // INCHI✔❌:                             p++;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         case '?':
    // INCHI✔❌:                             bond_parity = INCHI_PARITY_UNDEFINED;
    // INCHI✔❌:                             p++;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         default:
    // INCHI✔❌:                             bond_parity = 0;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     if (bond_parity)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         switch (*p)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             case '-':
    // INCHI✔❌:                                 bond_parityNM = INCHI_PARITY_ODD;
    // INCHI✔❌:                                 p++;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case '+':
    // INCHI✔❌:                                 bond_parityNM = INCHI_PARITY_EVEN;
    // INCHI✔❌:                                 p++;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case 'u':
    // INCHI✔❌:                                 bond_parityNM = INCHI_PARITY_UNKNOWN;
    // INCHI✔❌:                                 p++;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             case '?':
    // INCHI✔❌:                                 bond_parityNM = INCHI_PARITY_UNDEFINED;
    // INCHI✔❌:                                 p++;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                             default:
    // INCHI✔❌:                                 bond_parityNM = 0;
    // INCHI✔❌:                                 break;
    // INCHI✔❌:                         }
    // INCHI✔❌:                     }
    // INCHI✔❌:                     else
    // INCHI✔❌:                     {
    // INCHI✔❌:                         bond_parityNM = 0;
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     /* neighbor of the current atom */
    // INCHI✔❌:
    // INCHI✔❌:                     if (!isdigit( UCINT *p ))
    // INCHI✔❌:                     {
    // INCHI✔❌:                         num_atoms = INCHI_INP_ERROR_RET; /* error in input data */
    // INCHI✔❌:                         *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                         TREAT_ERR( *err, 0, "Wrong bonds data" );
    // INCHI✔❌:                         goto bypass_end_of_INChI;
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     neigh = (int) strtol( p, &q, 10 ) - 1;
    // INCHI✔❌:
    // INCHI✔❌:                     if (i >= num_atoms || neigh >= num_atoms)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         num_atoms = INCHI_INP_ERROR_RET; /* error in input data */
    // INCHI✔❌:                         *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                         TREAT_ERR( *err, 0, "Bond to nonexistent atom" );
    // INCHI✔❌:                         goto bypass_end_of_INChI;
    // INCHI✔❌:                     }
    // INCHI✔❌:                     p = q;
    // INCHI✔❌:                     bond_stereo1 = bond_stereo2 = 0;
    // INCHI✔❌:
    // INCHI✔❌:                     /* bond type & 2D stereo */
    // INCHI✔❌:                     switch (bond_char)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         case 'v':
    // INCHI✔❌:                             bond_type = INCHI_BOND_TYPE_SINGLE;
    // INCHI✔❌:                             bond_stereo1 = INCHI_BOND_STEREO_SINGLE_1EITHER;
    // INCHI✔❌:                             bond_stereo2 = INCHI_BOND_STEREO_SINGLE_2EITHER;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         case 'V':
    // INCHI✔❌:                             bond_type = INCHI_BOND_TYPE_SINGLE;
    // INCHI✔❌:                             bond_stereo1 = INCHI_BOND_STEREO_SINGLE_2EITHER;
    // INCHI✔❌:                             bond_stereo2 = INCHI_BOND_STEREO_SINGLE_1EITHER;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         case 'w':
    // INCHI✔❌:                             bond_type = INCHI_BOND_TYPE_DOUBLE;
    // INCHI✔❌:                             bond_stereo1 =
    // INCHI✔❌:                                 bond_stereo2 = INCHI_BOND_STEREO_DOUBLE_EITHER;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         case 's':
    // INCHI✔❌:                             bond_type = INCHI_BOND_TYPE_SINGLE;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         case 'd':
    // INCHI✔❌:                             bond_type = INCHI_BOND_TYPE_DOUBLE;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         case 't':
    // INCHI✔❌:                             bond_type = INCHI_BOND_TYPE_TRIPLE;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         case 'a':
    // INCHI✔❌:                             bond_type = INCHI_BOND_TYPE_ALTERN;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         case 'p':
    // INCHI✔❌:                             bond_type = INCHI_BOND_TYPE_SINGLE;
    // INCHI✔❌:                             bond_stereo1 = INCHI_BOND_STEREO_SINGLE_1UP;
    // INCHI✔❌:                             bond_stereo2 = INCHI_BOND_STEREO_SINGLE_2UP;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         case 'P':
    // INCHI✔❌:                             bond_type = INCHI_BOND_TYPE_SINGLE;
    // INCHI✔❌:                             bond_stereo1 = INCHI_BOND_STEREO_SINGLE_2UP;
    // INCHI✔❌:                             bond_stereo2 = INCHI_BOND_STEREO_SINGLE_1UP;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         case 'n':
    // INCHI✔❌:                             bond_type = INCHI_BOND_TYPE_SINGLE;
    // INCHI✔❌:                             bond_stereo1 = INCHI_BOND_STEREO_SINGLE_1DOWN;
    // INCHI✔❌:                             bond_stereo2 = INCHI_BOND_STEREO_SINGLE_2DOWN;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         case 'N':
    // INCHI✔❌:                             bond_type = INCHI_BOND_TYPE_SINGLE;
    // INCHI✔❌:                             bond_stereo1 = INCHI_BOND_STEREO_SINGLE_2DOWN;
    // INCHI✔❌:                             bond_stereo2 = INCHI_BOND_STEREO_SINGLE_1DOWN;
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         default:
    // INCHI✔❌:                             num_atoms = INCHI_INP_ERROR_RET; /* error */
    // INCHI✔❌:                             *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                             TREAT_ERR( *err, 0, "Wrong bond type" );
    // INCHI✔❌:                             goto bypass_end_of_INChI;
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     k = AT_NUM_BONDS( atom[i] )++;
    // INCHI✔❌:                     atom[i].bond_type[k] = bond_type;
    // INCHI✔❌:                     atom[i].bond_stereo[k] = bond_stereo1;
    // INCHI✔❌:                     atom[i].neighbor[k] = (ATOM_NUMBER) neigh;
    // INCHI✔❌:
    // INCHI✔❌:                     k2 = AT_NUM_BONDS( atom[neigh] )++;
    // INCHI✔❌:                     atom[neigh].bond_type[k2] = bond_type;
    // INCHI✔❌:                     atom[neigh].bond_stereo[k2] = bond_stereo2;
    // INCHI✔❌:                     atom[neigh].neighbor[k2] = (ATOM_NUMBER) i;
    // INCHI✔❌:
    // INCHI✔❌:                     bond_parity |= ( bond_parityNM << SB_PARITY_SHFT );
    // INCHI✔❌:
    // INCHI✔❌:                     if (bond_parity)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         if (max_len_stereo0D <= len_stereo0D)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             /* realloc atom_Stereo0D */
    // INCHI✔❌:
    // INCHI✔❌:                             inchi_Stereo0D *new_atom_stereo0D = CreateInchi_Stereo0D( max_len_stereo0D + num_atoms );
    // INCHI✔❌:
    // INCHI✔❌:                             if (!new_atom_stereo0D)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 num_atoms = INCHI_INP_FATAL_RET; /* fatal error: cannot allocate */
    // INCHI✔❌:                                 *err = INCHI_INP_FATAL_ERR;
    // INCHI✔❌:                                 TREAT_ERR( *err, 0, "Out of RAM" );
    // INCHI✔❌:                                 goto bypass_end_of_INChI;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             memcpy( new_atom_stereo0D, atom_stereo0D, len_stereo0D * sizeof( *atom_stereo0D ) );
    // INCHI✔❌:                             FreeInchi_Stereo0D( &atom_stereo0D );
    // INCHI✔❌:                             atom_stereo0D = new_atom_stereo0D;
    // INCHI✔❌:                             max_len_stereo0D += num_atoms;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         /* (a) i may be allene endpoint and     neigh = allene middle point or
    // INCHI✔❌:                             (b) i may be allene middle point and neigh = allene endpoint
    // INCHI✔❌:                             !!!!! CURRENTLY ONLY (b) IS ALLOWED !!!!!
    // INCHI✔❌:                         */
    // INCHI✔❌:                         atom_stereo0D[len_stereo0D].neighbor[1] = neigh; /* neigh < i */
    // INCHI✔❌:                         atom_stereo0D[len_stereo0D].neighbor[2] = i;
    // INCHI✔❌:                         atom_stereo0D[len_stereo0D].parity = bond_parity;
    // INCHI✔❌:                         atom_stereo0D[len_stereo0D].type = INCHI_StereoType_DoubleBond; /* incl allenes & cumulenes */
    // INCHI✔❌:                         len_stereo0D++;
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 if (!bItemIsOver || p - szLine_i2ia != res || i != num_atoms)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     num_atoms = INCHI_INP_ERROR_RET; /* error in input data */
    // INCHI✔❌:                     *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                     TREAT_ERR( *err, 0, "Wrong number of bonds" );
    // INCHI✔❌:                     goto bypass_end_of_INChI;
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 /********************** coordinates xml ****************************/
    // INCHI✔❌:
    // INCHI✔❌:                 pszCoord = (MOL_COORD*) inchi_malloc( inchi_max( num_atoms, 1 ) * sizeof( MOL_COORD ) );
    // INCHI✔❌:
    // INCHI✔❌:                 if (pszCoord)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     memset( pszCoord, ' ', inchi_max( num_atoms, 1 ) * sizeof( MOL_COORD ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:                     res = inchi_ios_gets( szLine_i2ia, sizeof(szLine_i2ia) - 1, inp_file, &bTooLongLine );
    // INCHI✔❌:
    // INCHI✔❌:                     if (res <= 0 ||
    // INCHI✔❌:                             /* compare the header */
    // INCHI✔❌:                             memcmp( szLine_i2ia, sStructRevXmlRevXYZ, sizeof( sStructRevXmlRevXYZ ) - 1 ) ||
    // INCHI✔❌:                             /* read (the head of) the coordinates (xml) line */
    // INCHI✔❌:                             0 >= ( res = inchi_ios_gets( szLine_i2ia, sizeof(szLine_i2ia) - 1, inp_file, &bTooLongLine ) ))
    // INCHI✔❌:                     {
    // INCHI✔❌:                         num_atoms = INCHI_INP_ERROR_RET; /* error in input data: atoms, bonds & coord must be present together */
    // INCHI✔❌:                         *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                         TREAT_ERR( *err, 0, "Missing atom coordinates data" );
    // INCHI✔❌:                         goto bypass_end_of_INChI;
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     i = 0;
    // INCHI✔❌:                     bItemIsOver = 0;
    // INCHI✔❌:                     res2 = bTooLongLine2 = -1; /* djb-rwth: ignoring LLVM warning: value used */
    // INCHI✔❌:                     p = szLine_i2ia;
    // INCHI✔❌:                     if (!memcmp( szLine_i2ia, sStructRevXmlRevXYZEnd, sizeof( sStructRevXmlRevXYZEnd ) - 1 ))
    // INCHI✔❌:                     {
    // INCHI✔❌:                         /* empty bonds section */
    // INCHI✔❌:                         res = 0;
    // INCHI✔❌:                         bItemIsOver = 1;
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     /* read all coordinates (xml), starting from atom 1 (not 0) */
    // INCHI✔❌:
    // INCHI✔❌:                     while (i < num_atoms)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         pos = p - szLine_i2ia;
    // INCHI✔❌:
    // INCHI✔❌:                         if (!bItemIsOver &&
    // INCHI✔❌:                             ( int )sizeof(szLine_i2ia) - res + pos > ( int )sizeof(szNextLine_i2ia))
    // INCHI✔❌:                         {
    // INCHI✔❌:
    // INCHI✔❌:                             /* load next line if possible */
    // INCHI✔❌:
    // INCHI✔❌:                             res2 = inchi_ios_gets( szNextLine_i2ia, sizeof(szNextLine_i2ia) - 1, inp_file, &bTooLongLine2 );
    // INCHI✔❌:
    // INCHI✔❌:                             if (res2 > 0 && memcmp( szNextLine_i2ia, sStructRevXmlRevXYZEnd, sizeof( sStructRevXmlRevXYZEnd ) - 1 ))
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 if (pos)
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     res -= pos;  /* number of chars left to process in szLine */
    // INCHI✔❌:                                     memmove(szLine_i2ia, p, res * sizeof(szLine_i2ia[0] ) ); /* move them to the start of the line */
    // INCHI✔❌:                                 }
    // INCHI✔❌:                                 memcpy(szLine_i2ia + res, szNextLine_i2ia, ( (long long)res2 + 1 ) * sizeof( szNextLine_i2ia[0] ) ); /* djb-rwth: cast operator added */
    // INCHI✔❌:                                 res += res2;
    // INCHI✔❌:                                 szLine_i2ia[res] = '\0';
    // INCHI✔❌:                                 bTooLongLine = bTooLongLine2;
    // INCHI✔❌:                                 p = szLine_i2ia;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             else
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 bItemIsOver = 1;
    // INCHI✔❌:                             }
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         /* coord, first char */
    // INCHI✔❌:
    // INCHI✔❌:                         if (*p == ';')
    // INCHI✔❌:                         {
    // INCHI✔❌:                             for (k = 0; k < NUM_COORD; k++)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 pszCoord[i][LEN_COORD*k + 4] = '0';
    // INCHI✔❌:                             }
    // INCHI✔❌:                             p++;
    // INCHI✔❌:                             i++;
    // INCHI✔❌:                             continue;
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         for (k = 0; k < 3; k++)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             double xyz;
    // INCHI✔❌:                             bNonZeroXYZ = 0;
    // INCHI✔❌:                             if (*p == ';')
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 pszCoord[i][LEN_COORD*k + 4] = '0';
    // INCHI✔❌:                                 xyz = 0.0;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             else if (*p == ',')
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 /* empty */
    // INCHI✔❌:                                 pszCoord[i][LEN_COORD*k + 4] = '0';
    // INCHI✔❌:                                 xyz = 0.0;
    // INCHI✔❌:                                 p++;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             else
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 xyz = strtod( p, &q );
    // INCHI✔❌:                                 bNonZeroXYZ = fabs( xyz ) > MIN_BOND_LENGTH;
    // INCHI✔❌:                                 if (q != NULL)
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     memcpy( pszCoord[i] + LEN_COORD*(long long)k, p, q - p ); /* djb-rwth: cast operator added */
    // INCHI✔❌:                                     if (*q == ',')
    // INCHI✔❌:                                         q++;
    // INCHI✔❌:                                     p = q;
    // INCHI✔❌:                                 }
    // INCHI✔❌:                                 else
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     pszCoord[i][LEN_COORD*k + 4] = '0';
    // INCHI✔❌:                                 }
    // INCHI✔❌:                             }
    // INCHI✔❌:
    // INCHI✔❌:                             switch (k)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 case 0:
    // INCHI✔❌:                                     atom[i].x = xyz;
    // INCHI✔❌:                                     b2D |= bNonZeroXYZ;
    // INCHI✔❌:                                     break;
    // INCHI✔❌:                                 case 1:
    // INCHI✔❌:                                     atom[i].y = xyz;
    // INCHI✔❌:                                     b2D |= bNonZeroXYZ;
    // INCHI✔❌:                                     break;
    // INCHI✔❌:                                 case 2:
    // INCHI✔❌:                                     b3D |= bNonZeroXYZ;
    // INCHI✔❌:                                     atom[i].z = xyz;
    // INCHI✔❌:                                     break;
    // INCHI✔❌:                             }
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         if (*p == ';')
    // INCHI✔❌:                         {
    // INCHI✔❌:                             p++;  /* end of this triple of coordinates */
    // INCHI✔❌:                             i++;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         else
    // INCHI✔❌:                         {
    // INCHI✔❌:                             num_atoms = INCHI_INP_ERROR_RET; /* error in input data: atoms, bonds & coord must be present together */
    // INCHI✔❌:                             *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                             TREAT_ERR( *err, 0, "Wrong atom coordinates data" );
    // INCHI✔❌:                             goto bypass_end_of_INChI;
    // INCHI✔❌:                         }
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     if (!bItemIsOver || p - szLine_i2ia != res || i != num_atoms)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         num_atoms = INCHI_INP_ERROR_RET; /* error in input data: atoms, bonds & coord must be present together */
    // INCHI✔❌:                         *err = INCHI_INP_ERROR_ERR;
    // INCHI✔❌:                         TREAT_ERR( *err, 0, "Wrong number of coordinates" );
    // INCHI✔❌:                         goto bypass_end_of_INChI;
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                     /* allocation failed */
    // INCHI✔❌:                     num_atoms = INCHI_INP_FATAL_RET;
    // INCHI✔❌:                     *err = INCHI_INP_FATAL_ERR;
    // INCHI✔❌:                     TREAT_ERR( *err, 0, "Out of RAM" );
    // INCHI✔❌:                     goto bypass_end_of_INChI;
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 /* set special valences and implicit H (xml) */
    // INCHI✔❌:
    // INCHI✔❌:                 b23D = b2D | b3D;
    // INCHI✔❌:                 b2D = b3D = 0; /* djb-rwth: ignoring LLVM warning: values used */
    // INCHI✔❌:
    // INCHI✔❌:                 if (at)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     if (!*at)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         int a1, a2, n1, n2, valence;
    // INCHI✔❌:                         int chem_bonds_valence;
    // INCHI✔❌:                         int    nX = 0, nY = 0, nZ = 0, nXYZ;
    // INCHI✔❌:                         *at = atom;
    // INCHI✔❌:
    // INCHI✔❌:                         /* special valences */
    // INCHI✔❌:                         for (bNonMetal = 0; bNonMetal < 1 /*2*/; bNonMetal++)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             for (a1 = 0; a1 < num_atoms; a1++)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 int num_bond_type[MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE + 1];
    // INCHI✔❌:
    // INCHI✔❌:                                 memset( num_bond_type, 0, sizeof( num_bond_type ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:
    // INCHI✔❌:                                 valence = AT_BONDS_VAL( atom, a1 ); /*  save atom valence if available */
    // INCHI✔❌:                                 AT_BONDS_VAL( atom, a1 ) = 0;
    // INCHI✔❌:
    // INCHI✔❌:                                 nX = nY = nZ = 0;
    // INCHI✔❌:
    // INCHI✔❌:                                 for (n1 = 0; n1 < AT_NUM_BONDS( atom[a1] ); n1++)
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     bond_type = atom[a1].bond_type[n1] - MIN_INPUT_BOND_TYPE;
    // INCHI✔❌:                                     if (bond_type < 0 || bond_type > MAX_INPUT_BOND_TYPE - MIN_INPUT_BOND_TYPE)
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         bond_type = 0; /* cannot happen */
    // INCHI✔❌:                                         TREAT_ERR( *err, 0, "Unknown bond type in InChI aux assigned as a single bond" );
    // INCHI✔❌:                                     }
    // INCHI✔❌:
    // INCHI✔❌:                                     num_bond_type[bond_type] ++;
    // INCHI✔❌:                                     nNumBonds++;
    // INCHI✔❌:
    // INCHI✔❌:                                     if (b23D)
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         neigh = atom[a1].neighbor[n1];
    // INCHI✔❌:                                         nX |= ( fabs( atom[a1].x - atom[neigh].x ) > MIN_BOND_LENGTH );
    // INCHI✔❌:                                         nY |= ( fabs( atom[a1].y - atom[neigh].y ) > MIN_BOND_LENGTH );
    // INCHI✔❌:                                         nZ |= ( fabs( atom[a1].z - atom[neigh].z ) > MIN_BOND_LENGTH );
    // INCHI✔❌:                                     }
    // INCHI✔❌:                                 }
    // INCHI✔❌:
    // INCHI✔❌:                                 chem_bonds_valence = 0;
    // INCHI✔❌:                                 for (n1 = 0; MIN_INPUT_BOND_TYPE + n1 <= 3 && MIN_INPUT_BOND_TYPE + n1 <= MAX_INPUT_BOND_TYPE; n1++)
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     chem_bonds_valence += ( MIN_INPUT_BOND_TYPE + n1 ) * num_bond_type[n1];
    // INCHI✔❌:                                 }
    // INCHI✔❌:
    // INCHI✔❌:                                 if (MIN_INPUT_BOND_TYPE <= INCHI_BOND_TYPE_ALTERN &&
    // INCHI✔❌:                                         INCHI_BOND_TYPE_ALTERN <= MAX_INPUT_BOND_TYPE &&
    // INCHI✔❌:                                         ( n2 = num_bond_type[INCHI_BOND_TYPE_ALTERN - MIN_INPUT_BOND_TYPE] ))
    // INCHI✔❌:                                 {
    // INCHI✔❌:
    // INCHI✔❌:                                     /* accept input aromatic bonds for now */
    // INCHI✔❌:
    // INCHI✔❌:                                     switch (n2)
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         case 2:
    // INCHI✔❌:                                             chem_bonds_valence += 3;  /* =A- */
    // INCHI✔❌:                                             break;
    // INCHI✔❌:                                         case 3:
    // INCHI✔❌:                                             chem_bonds_valence += 4;  /* =A< */
    // INCHI✔❌:                                             break;
    // INCHI✔❌:                                         default:
    // INCHI✔❌:                                             /*  if 1 or >= 4 aromatic bonds then replace such bonds with single bonds */
    // INCHI✔❌:                                             for (n1 = 0; n1 < AT_NUM_BONDS( atom[a1] ); n1++)
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 if (atom[a1].bond_type[n1] == INCHI_BOND_TYPE_ALTERN)
    // INCHI✔❌:                                                 {
    // INCHI✔❌:                                                     ATOM_NUMBER *p1;
    // INCHI✔❌:                                                     a2 = atom[a1].neighbor[n1];
    // INCHI✔❌:                                                     p1 = IN_NEIGH_LIST( atom[a2].neighbor, (ATOM_NUMBER) a1, AT_NUM_BONDS( atom[a2] ) );
    // INCHI✔❌:                                                     if (p1)
    // INCHI✔❌:                                                     {
    // INCHI✔❌:                                                         atom[a1].bond_type[n1] =
    // INCHI✔❌:                                                             atom[a2].bond_type[p1 - atom[a2].neighbor] = INCHI_BOND_TYPE_SINGLE;
    // INCHI✔❌:                                                     }
    // INCHI✔❌:                                                     else
    // INCHI✔❌:                                                     {
    // INCHI✔❌:                                                         *err = -2;  /*  Program error */
    // INCHI✔❌:                                                         TREAT_ERR( *err, 0, "Program error interpreting InChI aux" );
    // INCHI✔❌:                                                         num_atoms = 0;
    // INCHI✔❌:                                                         goto bypass_end_of_INChI; /*  no structure */
    // INCHI✔❌:                                                     }
    // INCHI✔❌:                                                 }
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                             chem_bonds_valence += n2;
    // INCHI✔❌:                                             *err |= 32; /*  Unrecognized aromatic bond(s) replaced with single */
    // INCHI✔❌:                                             TREAT_ERR( *err, 0, "Atom has 1 or more than 3 aromatic bonds" );
    // INCHI✔❌:                                             break;
    // INCHI✔❌:                                     }
    // INCHI✔❌:                                 }
    // INCHI✔❌:
    // INCHI✔❌:                                 /*************************************************************************************
    // INCHI✔❌:                                     *
    // INCHI✔❌:                                     *  Set number of hydrogen atoms
    // INCHI✔❌:                                     */
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     int num_iso_H;
    // INCHI✔❌:                                     num_iso_H = atom[a1].num_iso_H[1] + atom[a1].num_iso_H[2] + atom[a1].num_iso_H[3];
    // INCHI✔❌:                                     if (valence == ISOLATED_ATOM)
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         atom[a1].num_iso_H[0] = 0;
    // INCHI✔❌:                                     }
    // INCHI✔❌:                                     else
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         if (valence && valence >= chem_bonds_valence)
    // INCHI✔❌:                                         {
    // INCHI✔❌:                                             atom[a1].num_iso_H[0] = valence - chem_bonds_valence;
    // INCHI✔❌:                                         }
    // INCHI✔❌:                                         else
    // INCHI✔❌:                                         {
    // INCHI✔❌:                                             if (valence || bDoNotAddH)
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 atom[a1].num_iso_H[0] = 0;
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                             else
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 if (!bDoNotAddH)
    // INCHI✔❌:                                                 {
    // INCHI✔❌:                                                     atom[a1].num_iso_H[0] = -1; /* auto add H */
    // INCHI✔❌:                                                 }
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                         }
    // INCHI✔❌:                                     }
    // INCHI✔❌:                                 }
    // INCHI✔❌:                             }
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         nNumBonds /= 2;
    // INCHI✔❌:
    // INCHI✔❌:                         if (b23D && nNumBonds)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             nXYZ = nX + nY + nZ;
    // INCHI✔❌:                             b2D = ( nXYZ > 0 );
    // INCHI✔❌:                             b3D = ( nXYZ == 3 );
    // INCHI✔❌:                             *num_dimensions = b3D ? 3 : b2D ? 2 : 0;
    // INCHI✔❌:                             *num_bonds = nNumBonds;
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         /*======= 0D parities =================================*/
    // INCHI✔❌:                         if (len_stereo0D > 0 && atom_stereo0D && stereo0D)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             *stereo0D = atom_stereo0D;
    // INCHI✔❌:                             *num_stereo0D = len_stereo0D;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         else
    // INCHI✔❌:                         {
    // INCHI✔❌:                             FreeInchi_Stereo0D( &atom_stereo0D );
    // INCHI✔❌:                             *num_stereo0D = len_stereo0D = 0;
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         for (i = 0; i < len_stereo0D; i++)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             ATOM_NUMBER *p1, *p2;
    // INCHI✔❌:                             int     sb_ord_from_a1 = -1, sb_ord_from_a2 = -1, bEnd1 = 0, bEnd2 = 0;
    // INCHI✔❌:                             switch (atom_stereo0D[i].type)
    // INCHI✔❌:                             {
    // INCHI✔❌:
    // INCHI✔❌:                                 case INCHI_StereoType_Tetrahedral:
    // INCHI✔❌:                                     a1 = atom_stereo0D[i].central_atom;
    // INCHI✔❌:                                     if (atom_stereo0D[i].parity &&
    // INCHI✔❌:                                         ( AT_NUM_BONDS( atom[a1] ) == 3 || AT_NUM_BONDS( atom[a1] ) == 4 ))
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         int ii, kk = 0;
    // INCHI✔❌:                                         if (AT_NUM_BONDS( atom[a1] ) == 3)
    // INCHI✔❌:                                         {
    // INCHI✔❌:                                             atom_stereo0D[i].neighbor[kk++] = a1;
    // INCHI✔❌:                                         }
    // INCHI✔❌:                                         for (ii = 0; ii < AT_NUM_BONDS( atom[a1] ); ii++)
    // INCHI✔❌:                                         {
    // INCHI✔❌:                                             atom_stereo0D[i].neighbor[kk++] = atom[a1].neighbor[ii];
    // INCHI✔❌:                                         }
    // INCHI✔❌:                                     }
    // INCHI✔❌:
    // INCHI✔❌:                                     break;
    // INCHI✔❌:
    // INCHI✔❌:                                 case INCHI_StereoType_DoubleBond:
    // INCHI✔❌:
    // INCHI✔❌: #define MAX_CHAIN_LEN 20
    // INCHI✔❌:
    // INCHI✔❌:                                     a1 = atom_stereo0D[i].neighbor[1];
    // INCHI✔❌:                                     a2 = atom_stereo0D[i].neighbor[2];
    // INCHI✔❌:                                     p1 = IN_NEIGH_LIST( atom[a1].neighbor, (ATOM_NUMBER) a2, AT_NUM_BONDS( atom[a1] ) );
    // INCHI✔❌:                                     p2 = IN_NEIGH_LIST( atom[a2].neighbor, (ATOM_NUMBER) a1, AT_NUM_BONDS( atom[a2] ) );
    // INCHI✔❌:
    // INCHI✔❌:                                     if (!p1 || !p2)
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         atom_stereo0D[i].type = INCHI_StereoType_None;
    // INCHI✔❌:                                         atom_stereo0D[i].central_atom = NO_ATOM;
    // INCHI✔❌:                                         atom_stereo0D[i].neighbor[0] =
    // INCHI✔❌:                                             atom_stereo0D[i].neighbor[3] = -1;
    // INCHI✔❌:                                         *err |= 64; /* Error in cumulene stereo */
    // INCHI✔❌:                                         TREAT_ERR( *err, 0, "0D stereobond not recognized" );
    // INCHI✔❌:                                         break;
    // INCHI✔❌:                                     }
    // INCHI✔❌:
    // INCHI✔❌:                                     /* streobond, allene, or cumulene */
    // INCHI✔❌:
    // INCHI✔❌:                                     sb_ord_from_a1 = p1 - atom[a1].neighbor;
    // INCHI✔❌:                                     sb_ord_from_a2 = p2 - atom[a2].neighbor;
    // INCHI✔❌:
    // INCHI✔❌:                                     if (AT_NUM_BONDS( atom[a1] ) == 2 &&
    // INCHI✔❌:                                             atom[a1].bond_type[0] + atom[a1].bond_type[1] == 2 * INCHI_BOND_TYPE_DOUBLE &&
    // INCHI✔❌:                                             0 == inchi_NUMH2( atom, a1 ) &&
    // INCHI✔❌:                                             ( AT_NUM_BONDS( atom[a2] ) != 2 ||
    // INCHI✔❌:                                                 atom[a2].bond_type[0] + atom[a2].bond_type[1] != 2 * INCHI_BOND_TYPE_DOUBLE ))
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         bEnd2 = 1; /* a2 is the end-atom, a1 is middle atom */
    // INCHI✔❌:                                     }
    // INCHI✔❌:                                     if (AT_NUM_BONDS( atom[a2] ) == 2 &&
    // INCHI✔❌:                                             atom[a2].bond_type[0] + atom[a2].bond_type[1] == 2 * INCHI_BOND_TYPE_DOUBLE &&
    // INCHI✔❌:                                             0 == inchi_NUMH2( atom, a2 ) &&
    // INCHI✔❌:                                             ( AT_NUM_BONDS( atom[a1] ) != 2 ||
    // INCHI✔❌:                                                 atom[a1].bond_type[0] + atom[a1].bond_type[1] != 2 * INCHI_BOND_TYPE_DOUBLE ))
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         bEnd1 = 1; /* a1 is the end-atom, a2 is middle atom */
    // INCHI✔❌:                                     }
    // INCHI✔❌:
    // INCHI✔❌:                                     if (bEnd2 + bEnd1 == 1)
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         /* allene or cumulene */
    // INCHI✔❌:                                         ATOM_NUMBER  chain[MAX_CHAIN_LEN + 1], prev, cur, next;
    // INCHI✔❌:                                         if (bEnd2 && !bEnd1)
    // INCHI✔❌:                                         {
    // INCHI✔❌:                                             cur = a1;
    // INCHI✔❌:                                             a1 = a2;
    // INCHI✔❌:                                             a2 = cur;
    // INCHI✔❌:                                             sb_ord_from_a1 = sb_ord_from_a2;
    // INCHI✔❌:                                         }
    // INCHI✔❌:
    // INCHI✔❌:                                         sb_ord_from_a2 = -1;
    // INCHI✔❌:                                         cur = a1;
    // INCHI✔❌:                                         next = a2;
    // INCHI✔❌:                                         len = 0;
    // INCHI✔❌:                                         chain[len++] = cur;
    // INCHI✔❌:                                         chain[len++] = next;
    // INCHI✔❌:
    // INCHI✔❌:                                         while (len < MAX_CHAIN_LEN)
    // INCHI✔❌:                                         {
    // INCHI✔❌:                                             /* arbitrary very high upper limit to prevent infinite loop */
    // INCHI✔❌:
    // INCHI✔❌:                                             prev = cur;
    // INCHI✔❌:                                             cur = next;
    // INCHI✔❌:                                                 /* follow double bond path && avoid going back */
    // INCHI✔❌:                                             if (AT_NUM_BONDS( atom[cur] ) == 2 &&
    // INCHI✔❌:                                                     atom[cur].bond_type[0] + atom[cur].bond_type[1] == 2 * INCHI_BOND_TYPE_DOUBLE &&
    // INCHI✔❌:                                                     0 == inchi_NUMH2( atom, cur ))
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 next = atom[cur].neighbor[atom[cur].neighbor[0] == prev];
    // INCHI✔❌:                                                 chain[len++] = next;
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                             else
    // INCHI✔❌:                                                 break;
    // INCHI✔❌:                                         }
    // INCHI✔❌:
    // INCHI✔❌:                                         if (len > 2 &&
    // INCHI✔❌:                                             ( p2 = IN_NEIGH_LIST( atom[cur].neighbor, (ATOM_NUMBER) prev, AT_NUM_BONDS( atom[cur] ) ) ))
    // INCHI✔❌:                                         {
    // INCHI✔❌:                                             sb_ord_from_a2 = p2 - atom[cur].neighbor;
    // INCHI✔❌:                                             a2 = cur;
    // INCHI✔❌:
    // INCHI✔❌:                                             /* by design we need to pick up the first non-stereo-bond-neighbor as "sn"-atom */
    // INCHI✔❌:
    // INCHI✔❌:                                             atom_stereo0D[i].neighbor[0] = atom[a1].neighbor[sb_ord_from_a1 == 0];
    // INCHI✔❌:                                             atom_stereo0D[i].neighbor[1] = a1;
    // INCHI✔❌:                                             atom_stereo0D[i].neighbor[2] = a2;
    // INCHI✔❌:                                             atom_stereo0D[i].neighbor[3] = atom[a2].neighbor[sb_ord_from_a2 == 0];
    // INCHI✔❌:
    // INCHI✔❌:                                             if (len % 2)
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 atom_stereo0D[i].central_atom = chain[len / 2];
    // INCHI✔❌:                                                 atom_stereo0D[i].type = INCHI_StereoType_Allene;
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                             else
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 atom_stereo0D[i].central_atom = NO_ATOM;
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                         }
    // INCHI✔❌:                                         else
    // INCHI✔❌:                                         {
    // INCHI✔❌:                                             /* error */
    // INCHI✔❌:                                             atom_stereo0D[i].type = INCHI_StereoType_None;
    // INCHI✔❌:                                             atom_stereo0D[i].central_atom = NO_ATOM;
    // INCHI✔❌:                                             atom_stereo0D[i].neighbor[0] =
    // INCHI✔❌:                                                 atom_stereo0D[i].neighbor[3] = -1;
    // INCHI✔❌:                                             *err |= 64; /* Error in cumulene stereo */
    // INCHI✔❌:                                             TREAT_ERR( *err, 0, "Cumulene stereo not recognized (0D)" );
    // INCHI✔❌:                                         }
    // INCHI✔❌:
    // INCHI✔❌: #undef MAX_CHAIN_LEN
    // INCHI✔❌:                                     }
    // INCHI✔❌:                                     else
    // INCHI✔❌:                                     {
    // INCHI✔❌:
    // INCHI✔❌:                                         /****** a normal possibly stereogenic bond -- not an allene or cumulene *******/
    // INCHI✔❌:                                         /* by design we need to pick up the first non-stereo-bond-neighbor as "sn"-atom */
    // INCHI✔❌:
    // INCHI✔❌:                                         sb_ord_from_a1 = p1 - atom[a1].neighbor;
    // INCHI✔❌:                                         sb_ord_from_a2 = p2 - atom[a2].neighbor;
    // INCHI✔❌:                                         atom_stereo0D[i].neighbor[0] = atom[a1].neighbor[p1 == atom[a1].neighbor];
    // INCHI✔❌:                                         atom_stereo0D[i].neighbor[3] = atom[a2].neighbor[p2 == atom[a2].neighbor];
    // INCHI✔❌:                                         atom_stereo0D[i].central_atom = NO_ATOM;
    // INCHI✔❌:                                     }
    // INCHI✔❌:
    // INCHI✔❌:                                     if (atom_stereo0D[i].type != INCHI_StereoType_None &&
    // INCHI✔❌:                                             sb_ord_from_a1 >= 0 && sb_ord_from_a2 >= 0 &&
    // INCHI✔❌:                                             ATOM_PARITY_WELL_DEF( SB_PARITY_2( atom_stereo0D[i].parity ) ))
    // INCHI✔❌:                                     {
    // INCHI✔❌:
    // INCHI✔❌:                                         /* Detected well-defined disconnected stereo
    // INCHI✔❌:                                             * locate first non-metal neighbors */
    // INCHI✔❌:
    // INCHI✔❌:                                         int    a, n, j, /* k,*/ sb_ord, cur_neigh, min_neigh; /* djb-rwth: removing redundant variables */
    // INCHI✔❌:                                         for (k = 0; k < 2; k++)
    // INCHI✔❌:                                         {
    // INCHI✔❌:                                             a = k ? atom_stereo0D[i].neighbor[2] : atom_stereo0D[i].neighbor[1];
    // INCHI✔❌:                                             sb_ord = k ? sb_ord_from_a2 : sb_ord_from_a1;
    // INCHI✔❌:                                             min_neigh = num_atoms;
    // INCHI✔❌:                                             for (n = j = 0; j < AT_NUM_BONDS( atom[a] ); j++) /* djb-rwth: removing redundant code */
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 cur_neigh = atom[a].neighbor[j];
    // INCHI✔❌:                                                 if (j != sb_ord && !IS_METAL_ATOM( atom, cur_neigh ))
    // INCHI✔❌:                                                 {
    // INCHI✔❌:                                                     min_neigh = inchi_min( cur_neigh, min_neigh );
    // INCHI✔❌:                                                 }
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                             if (min_neigh < num_atoms)
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 atom_stereo0D[i].neighbor[k ? 3 : 0] = min_neigh;
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                             else
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 TREAT_ERR( *err, 0, "Cannot find non-metal stereobond neighor (0D)" );
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                         }
    // INCHI✔❌:                                     }
    // INCHI✔❌:
    // INCHI✔❌:                                     break;
    // INCHI✔❌:                             }
    // INCHI✔❌:                         }
    // INCHI✔❌:
    // INCHI✔❌:                         /* end of 0D parities extraction */
    // INCHI✔❌: /*exit_cycle:;*/
    // INCHI✔❌:                     }
    // INCHI✔❌:
    // INCHI✔❌:                     if (pInpAtomFlags)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         /* save chirality flag */
    // INCHI✔❌:                         *pInpAtomFlags |= InpAtomFlags;
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else if (atom)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     inchi_free( atom );
    // INCHI✔❌:                     atom = NULL;
    // INCHI✔❌:                 }
    // INCHI✔❌:
    // INCHI✔❌:                 if (pszCoord)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     inchi_free( pszCoord );
    // INCHI✔❌:                 }
    // INCHI✔❌:                 goto bypass_end_of_INChI;
    // INCHI✔❌:                 /*return num_atoms;*/
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         if (atom_stereo0D)
    // INCHI✔❌:             FreeInchi_Stereo0D(&atom_stereo0D);
    // INCHI✔❌:
    // INCHI✔❌:         /* end of struct. reading cycle, code never used? */
    // INCHI✔❌:
    // INCHI✔❌:         if (res <= 0)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (*err == INCHI_INP_ERROR_ERR)
    // INCHI✔❌:             {
    // INCHI✔❌:                 return num_atoms;
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             *err = INCHI_INP_EOF_ERR;
    // INCHI✔❌:             return INCHI_INP_EOF_RET; /* no more data */
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:     bypass_end_of_INChI:
    // INCHI✔❌:
    // INCHI✔❌:         /* cleanup */
    // INCHI✔❌:         if (num_atoms == INCHI_INP_ERROR_RET && atom_stereo0D)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (stereo0D && *stereo0D == atom_stereo0D)
    // INCHI✔❌:             {
    // INCHI✔❌:                 *stereo0D = NULL;
    // INCHI✔❌:                 *num_stereo0D = 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             FreeInchi_Stereo0D( &atom_stereo0D );
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         if (!memcmp( szLine_i2ia, sStructHdrXmlEnd, sizeof( sStructHdrXmlEnd ) - 1 ))
    // INCHI✔❌:         {
    // INCHI✔❌:             num_struct--;
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         if (!memcmp( szLine_i2ia, sStructHdrXml, sizeof( sStructHdrXml ) - 1 ))
    // INCHI✔❌:         {
    // INCHI✔❌:             num_struct++;
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         while (num_struct > 0 && 0 < inchi_ios_gets( szLine_i2ia, sizeof(szLine_i2ia) - 1, inp_file, &bTooLongLine ))
    // INCHI✔❌:         {
    // INCHI✔❌:             if (!memcmp( szLine_i2ia, sStructHdrXmlEnd, sizeof( sStructHdrXmlEnd ) - 1 ))
    // INCHI✔❌:             {
    // INCHI✔❌:                 num_struct--;
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 if (!memcmp( szLine_i2ia, sStructHdrXml, sizeof( sStructHdrXml ) - 1 ))
    // INCHI✔❌:                 {
    // INCHI✔❌:                     num_struct++;
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:         return num_atoms;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* djb-rwth: avoiding memory leak */
    // INCHI✔❌:     inchi_free(atom);
    // INCHI✔❌:     return num_atoms;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: InchiToInchiAtom

    let input = input.ok_or(SourceHeapError::NullPointer)?;
    let dimensions = dimensions.ok_or(SourceHeapError::NullPointer)?;
    let bond_count = bond_count.ok_or(SourceHeapError::NullPointer)?;
    let error = error.ok_or(SourceHeapError::NullPointer)?;
    *dimensions = 0;
    *bond_count = 0;

    if maximum_atom_count < 0 {
        return Err(SourceHeapError::SourceIntegerOverflow);
    }
    if let Some(slot) = atom_slot.as_deref_mut()
        && !slot.is_null()
        && maximum_atom_count != 0
    {
        let maximum = usize::try_from(maximum_atom_count)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        heap.slice_mut(*slot)?
            .get_mut(..maximum)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .fill(inchi_Atom::default());
    }

    let mut initial_stereo_capacity = 0_i32;
    if atom_slot.is_some()
        && let (Some(slot), Some(count)) =
            (stereo0d_slot.as_deref_mut(), stereo0d_count.as_deref_mut())
    {
        if !slot.is_null() && *count != 0 {
            if *count < 0 {
                return Err(SourceHeapError::SourceIntegerOverflow);
            }
            initial_stereo_capacity = *count;
            let capacity =
                usize::try_from(*count).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            heap.slice_mut(*slot)?
                .get_mut(..capacity)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .fill(inchi_Stereo0D::default());
        }
    }

    let line = heap.allocate_model_storage(vec![0_i8; INCHI_LINE_LEN as usize])?;
    let next_line = heap.allocate_model_storage(vec![0_i8; INCHI_LINE_ADD as usize])?;
    let atom_token =
        heap.allocate_model_storage(b"/rA:\0".iter().map(|byte| *byte as i8).collect::<Vec<_>>())?;

    let operation = (|| -> Result<i32, SourceHeapError> {
        let mut working_atom = SourceMutPointer::<inchi_Atom>::null();
        let mut working_atom_is_caller = false;
        let mut working_stereo = SourceMutPointer::<inchi_Stereo0D>::null();
        let mut working_coordinates = SourceMutPointer::<MOL_COORD>::null();
        let mut working_stereo_capacity = 0_i32;
        let mut xml_structure_depth = 0_i32;
        let mut too_long = 0_i32;

        macro_rules! drain_xml_structure {
            () => {{
                while xml_structure_depth > 0 {
                    let result = inchi_ios_gets(
                        heap,
                        line,
                        INCHI_LINE_LEN as i32 - 1,
                        input,
                        &mut too_long,
                    )?;
                    if result <= 0 {
                        break;
                    }
                    let content = source_c_string(heap.slice(line.as_const())?)?;
                    if source_starts_with(content, b"</structure") {
                        xml_structure_depth -= 1;
                    } else if source_starts_with(content, b"<structure") {
                        xml_structure_depth += 1;
                    }
                }
            }};
        }

        macro_rules! fail_input {
            ($status:expr, $code:expr, $message:expr) => {{
                let status = $status as i32;
                *error = $code as i32;
                add_source_error(heap, error_text, $message)?;
                let source_cleans_work = (input_type == tagInputType_INPUT_INCHI_PLAIN
                    && (status == INCHI_INP_ERROR_RET || status == INCHI_INP_FATAL_RET as i32))
                    || (input_type == tagInputType_INPUT_INCHI_XML
                        && status == INCHI_INP_ERROR_RET);
                if source_cleans_work {
                    if !working_stereo.is_null() {
                        if let Some(slot) = stereo0d_slot.as_deref_mut()
                            && *slot == working_stereo
                        {
                            *slot = SourceMutPointer::null();
                            if let Some(count) = stereo0d_count.as_deref_mut() {
                                *count = 0;
                            }
                        }
                        FreeInchi_Stereo0D(heap, Some(&mut working_stereo))?;
                    }
                    // Preserve the freed caller allocation as a readable snapshot of C's
                    // dangling output pointer; SourceHeap cannot expose use-after-free.
                    if !working_atom.is_null() && !working_atom_is_caller {
                        inchi_free(heap, working_atom)?;
                        working_atom = SourceMutPointer::null();
                    }
                    if !working_coordinates.is_null() {
                        inchi_free(heap, working_coordinates)?;
                        working_coordinates = SourceMutPointer::null();
                    }
                }
                if input_type == tagInputType_INPUT_INCHI_XML {
                    drain_xml_structure!();
                }
                return Ok(status);
            }};
        }

        let find_next_only = atom_slot.is_none();
        let mut xml_last_atom_count = 0_i32;

        'structure_attempt: loop {
            let mut atom_section = Vec::<i8>::new();
            let mut bond_section = Vec::<i8>::new();
            let mut coordinate_section = Vec::<i8>::new();
            let mut coordinate_requires_exact_end = false;
            let mut plain_coordinate_section_ended = false;
            let mut found_structure = false;
            let mut xml_sections_pending = false;
            let mut xml_coordinates_pending = false;
            let mut plain_after_atom_section = None::<Vec<i8>>;
            let mut plain_coordinates_missing = false;
            let mut plain_spanned_buffers = false;

            if input_type == tagInputType_INPUT_INCHI_PLAIN {
                let mut header_read = false;
                loop {
                    let mut result = inchi_ios_getsTab(
                        heap,
                        line,
                        INCHI_LINE_LEN as i32 - 1,
                        input,
                        &mut too_long,
                    )?;
                    if result <= 0 {
                        if *error == INCHI_INP_ERROR_ERR as i32 {
                            return Ok(0);
                        }
                        *error = INCHI_INP_EOF_ERR as i32;
                        return Ok(INCHI_INP_EOF_RET as i32);
                    }

                    let current = source_c_string(heap.slice(line.as_const())?)?.to_vec();
                    if too_long == 0 && source_starts_with(&current, b"Structure:") {
                        let mut cursor = b"Structure:".len();
                        let (long_id, end) = source_strtol_decimal(&current, cursor);
                        cursor = if current.get(end) == Some(&(b'.' as i8))
                            && current.get(end + 1) == Some(&(b' ' as i8))
                        {
                            end + 2
                        } else {
                            cursor
                        };
                        while current
                            .get(cursor)
                            .is_some_and(|byte| matches!(*byte as u8, b' ' | b'\n' | b'\r'))
                        {
                            cursor += 1;
                        }
                        if let Some(label) = sdf_label {
                            *heap
                                .slice_mut(label)?
                                .first_mut()
                                .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
                        }
                        if let Some(value) = sdf_value {
                            *heap
                                .slice_mut(value)?
                                .first_mut()
                                .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
                        }
                        if cursor < current.len() {
                            if let Some(equal_relative) = current[cursor..]
                                .iter()
                                .position(|byte| *byte == b'=' as i8)
                            {
                                let equal = cursor + equal_relative;
                                let mut length =
                                    (equal - cursor + 1).min(MAX_SDF_HEADER as usize - 1);
                                if let Some(label) = sdf_label {
                                    mystrncpy(
                                        heap,
                                        label,
                                        line.offset(cursor as i64)?.as_const(),
                                        length as u32,
                                    )?;
                                    let mut trimmed_length = i32::try_from(length)
                                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                                    lrtrim(heap, label, Some(&mut trimmed_length))?;
                                    length = usize::try_from(trimmed_length)
                                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                                    let _ = length;
                                }
                                let value_start = equal + 1;
                                if value_start < current.len() {
                                    let length = (current.len() - value_start + 1)
                                        .min(MAX_SDF_VALUE as usize - 1);
                                    if let Some(value) = sdf_value {
                                        mystrncpy(
                                            heap,
                                            value,
                                            line.offset(value_start as i64)?.as_const(),
                                            length as u32,
                                        )?;
                                    }
                                }
                            } else if let Some(missing_relative) =
                                source_find_bytes(&current[cursor..], b" is missing")
                            {
                                let missing = cursor + missing_relative;
                                let length =
                                    (missing - cursor + 1).min(MAX_SDF_HEADER as usize - 1);
                                if let Some(label) = sdf_label {
                                    mystrncpy(
                                        heap,
                                        label,
                                        line.offset(cursor as i64)?.as_const(),
                                        length as u32,
                                    )?;
                                }
                            }
                        }
                        if let Some(id) = id.as_deref_mut() {
                            *id = long_id;
                        }
                        header_read = true;
                        continue;
                    }

                    if !source_starts_with(&current, b"AuxInfo=") {
                        continue;
                    }
                    if !header_read {
                        if let Some(id) = id.as_deref_mut() {
                            *id = 0;
                        }
                        if let Some(label) = sdf_label {
                            *heap
                                .slice_mut(label)?
                                .first_mut()
                                .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
                        }
                        if let Some(value) = sdf_value {
                            *heap
                                .slice_mut(value)?
                                .first_mut()
                                .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
                        }
                    }
                    if let Some(slash) = current[b"AuxInfo=".len()..]
                        .iter()
                        .position(|byte| *byte == b'/' as i8)
                        .map(|position| position + b"AuxInfo=".len())
                        && current.get(slash + 1) == Some(&(b'/' as i8))
                        && current
                            .get(slash + 2)
                            .is_none_or(|byte| *byte == b'\n' as i8)
                    {
                        return Ok(0);
                    }

                    let mut cursor = FindToken(
                        heap,
                        input,
                        &mut too_long,
                        atom_token.as_const(),
                        4,
                        line,
                        INCHI_LINE_LEN as i32,
                        line,
                        &mut result,
                    )?;
                    if cursor.is_null() {
                        fail_input!(
                            INCHI_INP_ERROR_RET,
                            INCHI_INP_ERROR_ERR,
                            "Missing atom data"
                        );
                    }

                    let slash = {
                        let bytes = source_c_string(heap.slice(cursor.as_const())?)?;
                        bytes
                            .iter()
                            .position(|byte| *byte == b'/' as i8)
                            .map(|position| cursor.offset(position as i64))
                            .transpose()?
                            .unwrap_or(SourceMutPointer::null())
                    };
                    let mut slash = slash;
                    let mut item_is_over = i32::from(!slash.is_null() || too_long == 0);
                    cursor = LoadLine(
                        heap,
                        input,
                        &mut too_long,
                        &mut item_is_over,
                        &mut slash,
                        line,
                        INCHI_LINE_LEN as i32,
                        INCHI_LINE_ADD as i32,
                        cursor,
                        &mut result,
                    )?;
                    atom_section
                        .extend_from_slice(source_c_string(heap.slice(cursor.as_const())?)?);

                    plain_spanned_buffers = too_long != 0;
                    while too_long != 0 {
                        let loaded = inchi_ios_getsTab1(
                            heap,
                            next_line,
                            INCHI_LINE_ADD as i32 - 1,
                            input,
                            &mut too_long,
                        )?;
                        if loaded < 0 {
                            break;
                        }
                        atom_section
                            .extend_from_slice(source_c_string(heap.slice(next_line.as_const())?)?);
                    }

                    if let Some(bond_marker) = source_find_bytes(&atom_section, b"/rB:") {
                        plain_after_atom_section = Some(atom_section[bond_marker + 4..].to_vec());
                        atom_section.truncate(bond_marker);
                    } else if let Some(section_end) =
                        atom_section.iter().position(|byte| *byte == b'/' as i8)
                    {
                        atom_section.truncate(section_end);
                    }
                    found_structure = true;
                    break;
                }
            } else if input_type == tagInputType_INPUT_INCHI_XML {
                let mut header_read = false;
                let mut in_aux_info = false;
                loop {
                    let result = inchi_ios_gets(
                        heap,
                        line,
                        INCHI_LINE_LEN as i32 - 1,
                        input,
                        &mut too_long,
                    )?;
                    if result <= 0 {
                        if !working_stereo.is_null() {
                            let stereo_is_exposed = stereo0d_slot
                                .as_deref()
                                .is_some_and(|slot| *slot == working_stereo);
                            if stereo_is_exposed {
                                // C frees this allocation while leaving the caller's output pointer
                                // dangling. Keep the SourceHeap snapshot readable for exact observation.
                                working_stereo = SourceMutPointer::null();
                            } else {
                                FreeInchi_Stereo0D(heap, Some(&mut working_stereo))?;
                            }
                        }
                        if *error == INCHI_INP_ERROR_ERR as i32 {
                            return Ok(xml_last_atom_count);
                        }
                        *error = INCHI_INP_EOF_ERR as i32;
                        return Ok(INCHI_INP_EOF_RET as i32);
                    }
                    let current = source_c_string(heap.slice(line.as_const())?)?.to_vec();
                    if source_starts_with(&current, b"<structure") {
                        xml_last_atom_count = 0;
                        xml_structure_depth = 1;
                        if let Some(label) = sdf_label {
                            *heap
                                .slice_mut(label)?
                                .first_mut()
                                .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
                        }
                        if let Some(value) = sdf_value {
                            *heap
                                .slice_mut(value)?
                                .first_mut()
                                .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
                        }
                        let mut cursor = b"<structure".len();
                        let mut long_id = 0_i64;
                        if let Some(relative) = source_find_bytes(&current[cursor..], b"number=\"")
                        {
                            let start = cursor + relative + b"number=\"".len();
                            let (value, end) = source_strtol_decimal(&current, start);
                            long_id = value;
                            cursor = end + usize::from(current.get(end) == Some(&(b'\"' as i8)));
                        }
                        if let Some(relative) = source_find_bytes(&current[cursor..], b"id.name=\"")
                        {
                            let start = cursor + relative + b"id.name=\"".len();
                            if let Some(end_relative) = current[start..]
                                .iter()
                                .position(|byte| *byte == b'\"' as i8)
                            {
                                let end = start + end_relative;
                                let length = (end - start + 1).min(MAX_SDF_HEADER as usize - 1);
                                if let Some(label) = sdf_label {
                                    mystrncpy(
                                        heap,
                                        label,
                                        line.offset(start as i64)?.as_const(),
                                        length as u32,
                                    )?;
                                }
                                cursor = end + 1;
                            }
                        }
                        if let Some(relative) =
                            source_find_bytes(&current[cursor..], b"id.value=\"")
                        {
                            let start = cursor + relative + b"id.value=\"".len();
                            if let Some(end_relative) = current[start..]
                                .iter()
                                .position(|byte| *byte == b'\"' as i8)
                            {
                                let end = start + end_relative;
                                let length = (end - start + 1).min(MAX_SDF_VALUE as usize - 1);
                                if let Some(value) = sdf_value {
                                    mystrncpy(
                                        heap,
                                        value,
                                        line.offset(start as i64)?.as_const(),
                                        length as u32,
                                    )?;
                                }
                            }
                        }
                        if let Some(id) = id.as_deref_mut() {
                            *id = long_id;
                        }
                        header_read = true;
                        in_aux_info = false;
                        continue;
                    }

                    let error_prefix = b"<message type=\"error (no InChI)\" value=\"";
                    let fatal_prefix = b"<message type=\"fatal (aborted)\" value=\"";
                    if header_read
                        && (source_starts_with(&current, error_prefix)
                            || source_starts_with(&current, fatal_prefix))
                    {
                        let fatal = source_starts_with(&current, fatal_prefix);
                        let prefix_length = if fatal {
                            fatal_prefix.len()
                        } else {
                            error_prefix.len()
                        };
                        if !find_next_only
                            && let Some(end) = current[prefix_length..]
                                .iter()
                                .position(|byte| *byte == b'\"' as i8)
                        {
                            if let Some(error_text) = error_text {
                                *heap
                                    .slice_mut(error_text)?
                                    .first_mut()
                                    .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
                            }
                            let message: String = current[prefix_length..prefix_length + end]
                                .iter()
                                .map(|byte| *byte as u8 as char)
                                .collect();
                            add_source_error(heap, error_text, &message)?;
                        }
                        *error = if fatal {
                            INCHI_INP_FATAL_ERR as i32
                        } else {
                            INCHI_INP_ERROR_ERR as i32
                        };
                        drain_xml_structure!();
                        return Ok(if fatal {
                            INCHI_INP_FATAL_RET as i32
                        } else {
                            INCHI_INP_ERROR_RET
                        });
                    }
                    if header_read && source_starts_with(&current, b"<identifier.auxiliary-info") {
                        in_aux_info = true;
                        continue;
                    }
                    if header_read && source_starts_with(&current, b"</identifier.auxiliary-info") {
                        fail_input!(
                            INCHI_INP_ERROR_RET,
                            INCHI_INP_ERROR_ERR,
                            "Missing reversibility info"
                        );
                    }
                    if !(header_read
                        && in_aux_info
                        && source_starts_with(&current, b"<reversibility>"))
                    {
                        continue;
                    }

                    let result = inchi_ios_gets(
                        heap,
                        line,
                        INCHI_LINE_LEN as i32 - 1,
                        input,
                        &mut too_long,
                    )?;
                    if result <= 0 {
                        *error = INCHI_INP_EOF_ERR as i32;
                        return Ok(INCHI_INP_EOF_RET as i32);
                    }
                    if !source_starts_with(
                        source_c_string(heap.slice(line.as_const())?)?,
                        b"<atoms>",
                    ) {
                        header_read = false;
                        continue;
                    }
                    let mut read_atom_data = false;
                    loop {
                        let result = inchi_ios_gets(
                            heap,
                            line,
                            INCHI_LINE_LEN as i32 - 1,
                            input,
                            &mut too_long,
                        )?;
                        if result <= 0 {
                            if !read_atom_data {
                                *error = INCHI_INP_EOF_ERR as i32;
                                return Ok(INCHI_INP_EOF_RET as i32);
                            }
                            break;
                        }
                        read_atom_data = true;
                        let content = source_c_string(heap.slice(line.as_const())?)?;
                        if source_starts_with(content, b"</atoms>") {
                            break;
                        }
                        atom_section.extend_from_slice(content);
                    }
                    found_structure = true;
                    xml_sections_pending = true;
                    break;
                }
            }

            if !found_structure {
                inchi_free(heap, SourceMutPointer::<inchi_Atom>::null())?;
                return Ok(0);
            }

            let (source_atom_count, mut atom_cursor) = source_strtol_decimal(&atom_section, 0);
            let atom_count_i32 = source_atom_count as i32;
            if input_type == tagInputType_INPUT_INCHI_XML {
                xml_last_atom_count = atom_count_i32;
            }
            if atom_count_i32 == 0 || atom_cursor >= atom_section.len() {
                if input_type == tagInputType_INPUT_INCHI_XML {
                    *error = INCHI_INP_EOF_ERR as i32;
                }
                return Ok(0);
            }
            if atom_count_i32 < 0 {
                fail_input!(
                    INCHI_INP_ERROR_RET,
                    INCHI_INP_ERROR_ERR,
                    "Wrong number of atoms"
                );
            }
            let atom_count = usize::try_from(atom_count_i32)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let reuse_caller_atoms = atom_slot
                .as_deref()
                .is_some_and(|slot| !slot.is_null() && atom_count_i32 <= maximum_atom_count);

            if atom_slot.as_deref().is_none_or(|slot| slot.is_null()) {
                // A failed XML attempt leaves C's prior local allocation unreachable and a
                // subsequent attempt overwrites the local pointer with a fresh allocation.
                working_atom = SourceMutPointer::null();
                working_atom_is_caller = false;
            }

            if let Some(slot) = atom_slot.as_deref_mut()
                && !slot.is_null()
            {
                if atom_count_i32 > maximum_atom_count {
                    inchi_free(heap, *slot)?;
                    *slot = SourceMutPointer::null();
                    working_atom = SourceMutPointer::null();
                    working_atom_is_caller = false;
                } else {
                    let maximum = usize::try_from(maximum_atom_count)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    heap.slice_mut(*slot)?
                        .get_mut(..maximum)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .fill(inchi_Atom::default());
                    working_atom = *slot;
                    working_atom_is_caller = true;
                }
            }
            if working_atom.is_null() {
                working_atom = match CreateInchiAtom(heap, atom_count_i32 + 1) {
                    Ok(pointer) => pointer,
                    Err(SourceHeapError::AllocationFailed) => {
                        fail_input!(INCHI_INP_FATAL_RET, INCHI_INP_FATAL_ERR, "Out of RAM")
                    }
                    Err(error) => return Err(error),
                };
            }

            let source_initial_stereo_capacity = if atom_slot.is_some() {
                initial_stereo_capacity
            } else {
                0
            };
            if stereo0d_slot.as_deref().is_none_or(|slot| slot.is_null()) {
                working_stereo = SourceMutPointer::null();
                working_stereo_capacity = 0;
            }
            if let Some(slot) = stereo0d_slot.as_deref_mut()
                && !slot.is_null()
            {
                if atom_count_i32 > source_initial_stereo_capacity {
                    FreeInchi_Stereo0D(heap, Some(slot))?;
                    working_stereo = SourceMutPointer::null();
                    working_stereo_capacity = 0;
                } else {
                    let capacity = usize::try_from(source_initial_stereo_capacity)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    heap.slice_mut(*slot)?
                        .get_mut(..capacity)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .fill(inchi_Stereo0D::default());
                    working_stereo = *slot;
                    working_stereo_capacity = source_initial_stereo_capacity;
                }
            }
            if working_stereo.is_null() {
                working_stereo_capacity = atom_count_i32
                    .checked_add(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                working_stereo = match crate::source::base::readinch::CreateInchi_Stereo0D(
                    heap,
                    working_stereo_capacity,
                ) {
                    Ok(pointer) => pointer,
                    Err(SourceHeapError::AllocationFailed) => {
                        if !working_atom_is_caller {
                            inchi_free(heap, working_atom)?;
                            working_atom = SourceMutPointer::null();
                        }
                        fail_input!(INCHI_INP_FATAL_RET, INCHI_INP_FATAL_ERR, "Out of RAM")
                    }
                    Err(error) => return Err(error),
                };
            }

            let mut parsed_atoms = Vec::new();
            if parsed_atoms.try_reserve_exact(atom_count + 1).is_err() {
                fail_input!(
                    INCHI_INP_FATAL_RET as i32,
                    INCHI_INP_FATAL_ERR,
                    "Out of RAM"
                );
            }
            parsed_atoms.resize_with(atom_count + 1, inchi_Atom::default);
            let mut parsed_stereo = Vec::<inchi_Stereo0D>::new();
            if parsed_stereo.try_reserve_exact(atom_count + 1).is_err() {
                fail_input!(
                    INCHI_INP_FATAL_RET as i32,
                    INCHI_INP_FATAL_ERR,
                    "Out of RAM"
                );
            }
            let mut parsed_flags = 0_u64;
            match atom_section.get(atom_cursor).map(|byte| *byte as u8) {
                Some(b'c') => {
                    parsed_flags |= u64::from(FLAG_INP_AT_CHIRAL);
                    atom_cursor += 1;
                }
                Some(b'n') => {
                    parsed_flags |= u64::from(FLAG_INP_AT_NONCHIRAL);
                    atom_cursor += 1;
                }
                _ => {}
            }

            for atom_index in 0..atom_count {
                let Some(first) = atom_section.get(atom_cursor).map(|byte| *byte as u8) else {
                    fail_input!(
                        INCHI_INP_ERROR_RET,
                        INCHI_INP_ERROR_ERR,
                        if input_type == tagInputType_INPUT_INCHI_XML {
                            "Wrong atoms data"
                        } else {
                            "Wrong number of atoms"
                        }
                    );
                };
                if !first.is_ascii_alphabetic() || !first.is_ascii_uppercase() {
                    fail_input!(
                        INCHI_INP_ERROR_RET,
                        INCHI_INP_ERROR_ERR,
                        if input_type == tagInputType_INPUT_INCHI_XML {
                            "Wrong atoms data"
                        } else {
                            "Wrong number of atoms"
                        }
                    );
                }
                parsed_atoms[atom_index].elname[0] = first as i8;
                atom_cursor += 1;
                if atom_section
                    .get(atom_cursor)
                    .map(|byte| *byte as u8)
                    .is_some_and(|byte| byte.is_ascii_alphabetic() && byte.is_ascii_lowercase())
                {
                    parsed_atoms[atom_index].elname[1] = atom_section[atom_cursor];
                    atom_cursor += 1;
                }
                if atom_section
                    .get(atom_cursor)
                    .is_some_and(|byte| (*byte as u8).is_ascii_digit())
                {
                    let (value, end) = source_strtol_decimal(&atom_section, atom_cursor);
                    parsed_atoms[atom_index].num_iso_H[0] = value as i8;
                    if parsed_atoms[atom_index].num_iso_H[0] == 0 {
                        parsed_atoms[atom_index].num_iso_H[0] = -15;
                    }
                    atom_cursor = end;
                }
                let charge_sign = match atom_section.get(atom_cursor).map(|byte| *byte as u8) {
                    Some(b'+') => 1_i32,
                    Some(b'-') => -1_i32,
                    _ => 0,
                };
                parsed_atoms[atom_index].charge = charge_sign as i8;
                if charge_sign != 0 {
                    atom_cursor += 1;
                    if atom_section
                        .get(atom_cursor)
                        .is_some_and(|byte| (*byte as u8).is_ascii_digit())
                    {
                        let (value, end) = source_strtol_decimal(&atom_section, atom_cursor);
                        let magnitude = (value & i64::from(CHAR_MASK)) as i8;
                        parsed_atoms[atom_index].charge =
                            (charge_sign * i32::from(magnitude)) as i8;
                        atom_cursor = end;
                    }
                }
                if atom_section.get(atom_cursor) == Some(&(b'.' as i8)) {
                    atom_cursor += 1;
                    if atom_section
                        .get(atom_cursor)
                        .is_some_and(|byte| (*byte as u8).is_ascii_digit())
                    {
                        let (value, end) = source_strtol_decimal(&atom_section, atom_cursor);
                        parsed_atoms[atom_index].radical = value as i8;
                        atom_cursor = end;
                    }
                }
                if atom_section.get(atom_cursor) == Some(&(b'i' as i8)) {
                    atom_cursor += 1;
                    if atom_section
                        .get(atom_cursor)
                        .is_some_and(|byte| (*byte as u8).is_ascii_digit())
                    {
                        let (value, end) = source_strtol_decimal(&atom_section, atom_cursor);
                        parsed_atoms[atom_index].isotopic_mass = value as i16;
                        atom_cursor = end;
                    }
                }
                let parity = match atom_section.get(atom_cursor).map(|byte| *byte as u8) {
                    Some(b'o') => tagINCHIStereoParity0D_INCHI_PARITY_ODD as i8,
                    Some(b'e') => tagINCHIStereoParity0D_INCHI_PARITY_EVEN as i8,
                    Some(b'u') => tagINCHIStereoParity0D_INCHI_PARITY_UNKNOWN as i8,
                    Some(b'?') => tagINCHIStereoParity0D_INCHI_PARITY_UNDEFINED as i8,
                    _ => 0,
                };
                if parity != 0 {
                    atom_cursor += 1;
                    parsed_stereo.push(inchi_Stereo0D {
                        central_atom: atom_index as i16,
                        parity,
                        type_: tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral as i8,
                        ..inchi_Stereo0D::default()
                    });
                    let stereo_index = parsed_stereo.len() - 1;
                    heap.slice_mut(working_stereo)?[stereo_index] =
                        parsed_stereo[stereo_index].clone();
                }
                for (isotope, marker) in b"hdt".iter().enumerate() {
                    if atom_section.get(atom_cursor).map(|byte| *byte as u8) == Some(*marker) {
                        parsed_atoms[atom_index].num_iso_H[isotope + 1] = 1;
                        atom_cursor += 1;
                        if atom_section
                            .get(atom_cursor)
                            .is_some_and(|byte| (*byte as u8).is_ascii_digit())
                        {
                            let (value, end) = source_strtol_decimal(&atom_section, atom_cursor);
                            parsed_atoms[atom_index].num_iso_H[isotope + 1] = value as i8;
                            atom_cursor = end;
                        }
                    }
                }
                heap.slice_mut(working_atom)?[atom_index] = parsed_atoms[atom_index].clone();
            }
            if atom_cursor != atom_section.len() {
                fail_input!(
                    INCHI_INP_ERROR_RET,
                    INCHI_INP_ERROR_ERR,
                    "Wrong number of atoms"
                );
            }

            if input_type == tagInputType_INPUT_INCHI_PLAIN {
                let Some(remaining) = plain_after_atom_section.take() else {
                    fail_input!(
                        INCHI_INP_ERROR_RET,
                        INCHI_INP_ERROR_ERR,
                        "Missing bonds data"
                    );
                };
                if let Some(coordinate_marker) = source_find_bytes(&remaining, b"/rC:") {
                    bond_section.extend_from_slice(&remaining[..coordinate_marker]);
                    let coordinates = &remaining[coordinate_marker + 4..];
                    let coordinate_end = coordinates.iter().position(|byte| *byte == b'/' as i8);
                    plain_coordinate_section_ended = coordinate_end.is_some();
                    coordinate_requires_exact_end =
                        coordinate_end.is_some() || plain_spanned_buffers;
                    let coordinate_end = coordinate_end.unwrap_or(coordinates.len());
                    coordinate_section.extend_from_slice(&coordinates[..coordinate_end]);
                } else {
                    let bond_end = remaining
                        .iter()
                        .position(|byte| *byte == b'/' as i8)
                        .unwrap_or(remaining.len());
                    bond_section.extend_from_slice(&remaining[..bond_end]);
                    plain_coordinates_missing = true;
                }
            }

            if xml_sections_pending {
                let result =
                    inchi_ios_gets(heap, line, INCHI_LINE_LEN as i32 - 1, input, &mut too_long)?;
                if result <= 0 {
                    return Ok(0);
                }
                if !source_starts_with(source_c_string(heap.slice(line.as_const())?)?, b"<bonds>") {
                    continue 'structure_attempt;
                }
                let mut read_bond_data = false;
                loop {
                    let result = inchi_ios_gets(
                        heap,
                        line,
                        INCHI_LINE_LEN as i32 - 1,
                        input,
                        &mut too_long,
                    )?;
                    if result <= 0 {
                        if !read_bond_data {
                            fail_input!(INCHI_INP_ERROR_RET, INCHI_INP_ERROR_ERR, "");
                        }
                        break;
                    }
                    read_bond_data = true;
                    let content = source_c_string(heap.slice(line.as_const())?)?;
                    if source_starts_with(content, b"</bonds>") {
                        break;
                    }
                    bond_section.extend_from_slice(content);
                }
                xml_coordinates_pending = true;
            }

            let mut current_atom = 1_usize;
            let mut cursor = 0_usize;
            while current_atom < atom_count {
                if cursor >= bond_section.len() {
                    if input_type == tagInputType_INPUT_INCHI_XML {
                        fail_input!(INCHI_INP_ERROR_RET, INCHI_INP_ERROR_ERR, "Wrong bonds data");
                    }
                    break;
                }
                if bond_section.get(cursor) == Some(&(b';' as i8)) {
                    cursor += 1;
                    current_atom += 1;
                    continue;
                }
                let Some(bond_character) = bond_section.get(cursor).map(|byte| *byte as u8) else {
                    fail_input!(INCHI_INP_ERROR_RET, INCHI_INP_ERROR_ERR, "Wrong bonds data");
                };
                if !bond_character.is_ascii_alphabetic() {
                    fail_input!(INCHI_INP_ERROR_RET, INCHI_INP_ERROR_ERR, "Wrong bonds data");
                }
                cursor += 1;
                let bond_parity = match bond_section.get(cursor).map(|byte| *byte as u8) {
                    Some(b'-') => tagINCHIStereoParity0D_INCHI_PARITY_ODD as i32,
                    Some(b'+') => tagINCHIStereoParity0D_INCHI_PARITY_EVEN as i32,
                    Some(b'u') => tagINCHIStereoParity0D_INCHI_PARITY_UNKNOWN as i32,
                    Some(b'?') => tagINCHIStereoParity0D_INCHI_PARITY_UNDEFINED as i32,
                    _ => 0,
                };
                if bond_parity != 0 {
                    cursor += 1;
                }
                let bond_parity_nonmetal = if bond_parity != 0 {
                    let parity = match bond_section.get(cursor).map(|byte| *byte as u8) {
                        Some(b'-') => tagINCHIStereoParity0D_INCHI_PARITY_ODD as i32,
                        Some(b'+') => tagINCHIStereoParity0D_INCHI_PARITY_EVEN as i32,
                        Some(b'u') => tagINCHIStereoParity0D_INCHI_PARITY_UNKNOWN as i32,
                        Some(b'?') => tagINCHIStereoParity0D_INCHI_PARITY_UNDEFINED as i32,
                        _ => 0,
                    };
                    if parity != 0 {
                        cursor += 1;
                    }
                    parity
                } else {
                    0
                };
                if !bond_section
                    .get(cursor)
                    .is_some_and(|byte| (*byte as u8).is_ascii_digit())
                {
                    fail_input!(INCHI_INP_ERROR_RET, INCHI_INP_ERROR_ERR, "Wrong bonds data");
                }
                let (neighbor_number, end) = source_strtol_decimal(&bond_section, cursor);
                cursor = end;
                let neighbor = (neighbor_number as i32).wrapping_sub(1);
                if neighbor < 0 || neighbor as usize >= atom_count {
                    fail_input!(
                        INCHI_INP_ERROR_RET,
                        INCHI_INP_ERROR_ERR,
                        "Bond to nonexistent atom"
                    );
                }
                let neighbor = neighbor as usize;
                let (bond_type, stereo1, stereo2) = match bond_character {
                    b'v' => (
                        tagINCHIBondType_INCHI_BOND_TYPE_SINGLE as i8,
                        tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1EITHER as i8,
                        tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2EITHER as i8,
                    ),
                    b'V' => (
                        tagINCHIBondType_INCHI_BOND_TYPE_SINGLE as i8,
                        tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2EITHER as i8,
                        tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1EITHER as i8,
                    ),
                    b'w' => (
                        tagINCHIBondType_INCHI_BOND_TYPE_DOUBLE as i8,
                        tagINCHIBondStereo2D_INCHI_BOND_STEREO_DOUBLE_EITHER as i8,
                        tagINCHIBondStereo2D_INCHI_BOND_STEREO_DOUBLE_EITHER as i8,
                    ),
                    b's' => (tagINCHIBondType_INCHI_BOND_TYPE_SINGLE as i8, 0, 0),
                    b'd' => (tagINCHIBondType_INCHI_BOND_TYPE_DOUBLE as i8, 0, 0),
                    b't' => (tagINCHIBondType_INCHI_BOND_TYPE_TRIPLE as i8, 0, 0),
                    b'a' => (tagINCHIBondType_INCHI_BOND_TYPE_ALTERN as i8, 0, 0),
                    b'p' => (
                        tagINCHIBondType_INCHI_BOND_TYPE_SINGLE as i8,
                        tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1UP as i8,
                        tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2UP as i8,
                    ),
                    b'P' => (
                        tagINCHIBondType_INCHI_BOND_TYPE_SINGLE as i8,
                        tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2UP as i8,
                        tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1UP as i8,
                    ),
                    b'n' => (
                        tagINCHIBondType_INCHI_BOND_TYPE_SINGLE as i8,
                        tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1DOWN as i8,
                        tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2DOWN as i8,
                    ),
                    b'N' => (
                        tagINCHIBondType_INCHI_BOND_TYPE_SINGLE as i8,
                        tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2DOWN as i8,
                        tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1DOWN as i8,
                    ),
                    _ => fail_input!(INCHI_INP_ERROR_RET, INCHI_INP_ERROR_ERR, "Wrong bond type"),
                };
                let current_order = usize::try_from(parsed_atoms[current_atom].num_bonds)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                if current_order >= parsed_atoms[current_atom].neighbor.len() {
                    fail_input!(INCHI_INP_ERROR_RET, INCHI_INP_ERROR_ERR, "Wrong bonds data");
                }
                parsed_atoms[current_atom].neighbor[current_order] = neighbor as i16;
                parsed_atoms[current_atom].bond_type[current_order] = bond_type;
                parsed_atoms[current_atom].bond_stereo[current_order] = stereo1;
                parsed_atoms[current_atom].num_bonds += 1;
                let neighbor_order = usize::try_from(parsed_atoms[neighbor].num_bonds)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                if neighbor_order >= parsed_atoms[neighbor].neighbor.len() {
                    fail_input!(INCHI_INP_ERROR_RET, INCHI_INP_ERROR_ERR, "Wrong bonds data");
                }
                parsed_atoms[neighbor].neighbor[neighbor_order] = current_atom as i16;
                parsed_atoms[neighbor].bond_type[neighbor_order] = bond_type;
                parsed_atoms[neighbor].bond_stereo[neighbor_order] = stereo2;
                parsed_atoms[neighbor].num_bonds += 1;
                heap.slice_mut(working_atom)?[current_atom] = parsed_atoms[current_atom].clone();
                if neighbor != current_atom {
                    heap.slice_mut(working_atom)?[neighbor] = parsed_atoms[neighbor].clone();
                }
                let combined_parity = bond_parity | (bond_parity_nonmetal << SB_PARITY_SHFT);
                if combined_parity != 0 {
                    if working_stereo_capacity <= parsed_stereo.len() as i32 {
                        let new_capacity = working_stereo_capacity
                            .checked_add(atom_count_i32)
                            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                        let replacement = match crate::source::base::readinch::CreateInchi_Stereo0D(
                            heap,
                            new_capacity,
                        ) {
                            Ok(pointer) => pointer,
                            Err(SourceHeapError::AllocationFailed) => {
                                fail_input!(INCHI_INP_FATAL_RET, INCHI_INP_FATAL_ERR, "Out of RAM")
                            }
                            Err(error) => return Err(error),
                        };
                        let existing = parsed_stereo.len();
                        if existing != 0 {
                            let prior = heap.slice(working_stereo.as_const())?[..existing].to_vec();
                            heap.slice_mut(replacement)?[..existing].clone_from_slice(&prior);
                        }
                        FreeInchi_Stereo0D(heap, Some(&mut working_stereo))?;
                        working_stereo = replacement;
                        working_stereo_capacity = new_capacity;
                    }
                    parsed_stereo.push(inchi_Stereo0D {
                        neighbor: [0, neighbor as i16, current_atom as i16, 0],
                        parity: combined_parity as i8,
                        type_: tagINCHIStereoType0D_INCHI_StereoType_DoubleBond as i8,
                        ..inchi_Stereo0D::default()
                    });
                    let stereo_index = parsed_stereo.len() - 1;
                    heap.slice_mut(working_stereo)?[stereo_index] =
                        parsed_stereo[stereo_index].clone();
                }
            }
            if current_atom != atom_count || cursor != bond_section.len() {
                fail_input!(
                    INCHI_INP_ERROR_RET,
                    INCHI_INP_ERROR_ERR,
                    "Wrong number of bonds"
                );
            }

            if plain_coordinates_missing {
                fail_input!(
                    INCHI_INP_ERROR_RET,
                    INCHI_INP_ERROR_ERR,
                    "Missing atom coordinates data"
                );
            }

            if xml_coordinates_pending {
                let result =
                    inchi_ios_gets(heap, line, INCHI_LINE_LEN as i32 - 1, input, &mut too_long)?;
                if result <= 0
                    || !source_starts_with(source_c_string(heap.slice(line.as_const())?)?, b"<xyz>")
                {
                    fail_input!(
                        INCHI_INP_ERROR_RET,
                        INCHI_INP_ERROR_ERR,
                        "Missing atom coordinates data"
                    );
                }
                let mut read_coordinate_data = false;
                loop {
                    let result = inchi_ios_gets(
                        heap,
                        line,
                        INCHI_LINE_LEN as i32 - 1,
                        input,
                        &mut too_long,
                    )?;
                    if result <= 0 {
                        if !read_coordinate_data {
                            fail_input!(
                                INCHI_INP_ERROR_RET,
                                INCHI_INP_ERROR_ERR,
                                "Missing atom coordinates data"
                            );
                        }
                        break;
                    }
                    read_coordinate_data = true;
                    let content = source_c_string(heap.slice(line.as_const())?)?;
                    if source_starts_with(content, b"</xyz>") {
                        break;
                    }
                    coordinate_section.extend_from_slice(content);
                }
                coordinate_requires_exact_end = true;
            }

            working_coordinates = match heap.allocate(vec![[b' ' as i8; 32]; atom_count.max(1)]) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => {
                    fail_input!(INCHI_INP_FATAL_RET, INCHI_INP_FATAL_ERR, "Out of RAM")
                }
                Err(error) => return Err(error),
            };

            let mut has_xy = false;
            let mut has_z = false;
            cursor = 0;
            let mut parsed_coordinate_count = 0_usize;
            for (atom_index, atom) in parsed_atoms.iter_mut().take(atom_count).enumerate() {
                if plain_coordinate_section_ended && cursor >= coordinate_section.len() {
                    break;
                }
                if coordinate_section.get(cursor) == Some(&(b';' as i8)) {
                    cursor += 1;
                    parsed_coordinate_count += 1;
                    continue;
                }
                for coordinate in 0..3 {
                    let (value, end) = if coordinate_section.get(cursor) == Some(&(b';' as i8)) {
                        (0.0, cursor)
                    } else if coordinate_section.get(cursor) == Some(&(b',' as i8)) {
                        (0.0, cursor + 1)
                    } else {
                        source_strtod_decimal(&coordinate_section, cursor)
                    };
                    cursor = end;
                    let nonzero = value.abs() > 1.0e-6;
                    match coordinate {
                        0 => {
                            atom.x = value;
                            has_xy |= nonzero;
                        }
                        1 => {
                            atom.y = value;
                            has_xy |= nonzero;
                        }
                        _ => {
                            atom.z = value;
                            has_z |= nonzero;
                        }
                    }
                    if coordinate_section.get(cursor) == Some(&(b',' as i8)) {
                        cursor += 1;
                    }
                }
                heap.slice_mut(working_atom)?[atom_index] = atom.clone();
                if coordinate_section.get(cursor) != Some(&(b';' as i8)) {
                    fail_input!(
                        INCHI_INP_ERROR_RET,
                        INCHI_INP_ERROR_ERR,
                        "Wrong atom coordinates data"
                    );
                }
                cursor += 1;
                parsed_coordinate_count += 1;
            }
            if (plain_coordinate_section_ended && parsed_coordinate_count != atom_count)
                || (coordinate_requires_exact_end && cursor != coordinate_section.len())
            {
                fail_input!(
                    INCHI_INP_ERROR_RET,
                    INCHI_INP_ERROR_ERR,
                    "Wrong number of coordinates"
                );
            }

            if !reuse_caller_atoms && !find_next_only {
                let coordinates_present = has_xy || has_z;
                let mut counted_bonds = 0_i32;
                let mut last_x = 0_i32;
                let mut last_y = 0_i32;
                let mut last_z = 0_i32;
                for atom_index in 0..atom_count {
                    let mut bond_type_counts = [0_i32; 4];
                    let valence = i32::from(parsed_atoms[atom_index].num_iso_H[0]);
                    parsed_atoms[atom_index].num_iso_H[0] = 0;
                    last_x = 0;
                    last_y = 0;
                    last_z = 0;
                    let local_bond_count = usize::try_from(parsed_atoms[atom_index].num_bonds)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    for order in 0..local_bond_count {
                        let mut index = i32::from(parsed_atoms[atom_index].bond_type[order])
                            - MIN_INPUT_BOND_TYPE as i32;
                        if index < 0
                            || index > MAX_INPUT_BOND_TYPE as i32 - MIN_INPUT_BOND_TYPE as i32
                        {
                            index = 0;
                            add_source_error(
                                heap,
                                error_text,
                                "Unknown bond type in InChI aux assigned as a single bond",
                            )?;
                        }
                        bond_type_counts[index as usize] += 1;
                        counted_bonds += 1;
                        if coordinates_present {
                            let neighbor =
                                usize::try_from(parsed_atoms[atom_index].neighbor[order])
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            last_x |= i32::from(
                                (parsed_atoms[atom_index].x - parsed_atoms[neighbor].x).abs()
                                    > 1.0e-6,
                            );
                            last_y |= i32::from(
                                (parsed_atoms[atom_index].y - parsed_atoms[neighbor].y).abs()
                                    > 1.0e-6,
                            );
                            last_z |= i32::from(
                                (parsed_atoms[atom_index].z - parsed_atoms[neighbor].z).abs()
                                    > 1.0e-6,
                            );
                        }
                    }
                    let mut chemical_valence = 0_i32;
                    for bond_type in MIN_INPUT_BOND_TYPE..=3.min(MAX_INPUT_BOND_TYPE) {
                        chemical_valence += bond_type as i32
                            * bond_type_counts[(bond_type - MIN_INPUT_BOND_TYPE) as usize];
                    }
                    let aromatic_index =
                        tagINCHIBondType_INCHI_BOND_TYPE_ALTERN as i32 - MIN_INPUT_BOND_TYPE as i32;
                    if aromatic_index >= 0 && aromatic_index < bond_type_counts.len() as i32 {
                        let aromatic_count = bond_type_counts[aromatic_index as usize];
                        match aromatic_count {
                            0 => {}
                            2 => chemical_valence += 3,
                            3 => chemical_valence += 4,
                            count => {
                                for order in 0..local_bond_count {
                                    if parsed_atoms[atom_index].bond_type[order]
                                        == tagINCHIBondType_INCHI_BOND_TYPE_ALTERN as i8
                                    {
                                        let neighbor = usize::try_from(
                                            parsed_atoms[atom_index].neighbor[order],
                                        )
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                        let neighbor_count =
                                            parsed_atoms[neighbor].num_bonds as i32;
                                        let reciprocal = is_in_the_slist(
                                            Some(&parsed_atoms[neighbor].neighbor),
                                            atom_index as i16,
                                            neighbor_count,
                                        )?;
                                        let Some(reciprocal) = reciprocal else {
                                            *error = -2;
                                            add_source_error(
                                                heap,
                                                error_text,
                                                "Program error interpreting InChI aux",
                                            )?;
                                            return Ok(INCHI_INP_FATAL_RET as i32);
                                        };
                                        parsed_atoms[atom_index].bond_type[order] =
                                            tagINCHIBondType_INCHI_BOND_TYPE_SINGLE as i8;
                                        parsed_atoms[neighbor].bond_type[reciprocal] =
                                            tagINCHIBondType_INCHI_BOND_TYPE_SINGLE as i8;
                                    }
                                }
                                chemical_valence += count;
                                *error |= 32;
                                add_source_error(
                                    heap,
                                    error_text,
                                    "Atom has 1 or more than 3 aromatic bonds",
                                )?;
                            }
                        }
                    }
                    parsed_atoms[atom_index].num_iso_H[0] = if valence == -15 {
                        0
                    } else if valence != 0 && valence >= chemical_valence {
                        (valence - chemical_valence) as i8
                    } else if valence != 0 || do_not_add_h != 0 {
                        0
                    } else {
                        -1
                    };
                }
                counted_bonds /= 2;
                if coordinates_present && counted_bonds != 0 {
                    let coordinate_axes = last_x + last_y + last_z;
                    *dimensions = if coordinate_axes == 3 {
                        3
                    } else if coordinate_axes > 0 {
                        2
                    } else {
                        0
                    };
                    *bond_count = counted_bonds;
                }

                for stereo in &mut parsed_stereo {
                    match u32::try_from(stereo.type_).unwrap_or(0) {
                        tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral => {
                            let center = usize::try_from(stereo.central_atom)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let neighbors = parsed_atoms[center].num_bonds as usize;
                            if stereo.parity != 0 && matches!(neighbors, 3 | 4) {
                                let mut output = 0;
                                if neighbors == 3 {
                                    stereo.neighbor[output] = center as i16;
                                    output += 1;
                                }
                                for neighbor in parsed_atoms[center].neighbor.iter().take(neighbors)
                                {
                                    stereo.neighbor[output] = *neighbor;
                                    output += 1;
                                }
                            }
                        }
                        tagINCHIStereoType0D_INCHI_StereoType_DoubleBond => {
                            let mut first = usize::try_from(stereo.neighbor[1])
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let mut second = usize::try_from(stereo.neighbor[2])
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let first_order = is_in_the_slist(
                                Some(&parsed_atoms[first].neighbor),
                                second as i16,
                                parsed_atoms[first].num_bonds as i32,
                            )?;
                            let second_order = is_in_the_slist(
                                Some(&parsed_atoms[second].neighbor),
                                first as i16,
                                parsed_atoms[second].num_bonds as i32,
                            )?;
                            let (Some(mut first_order), Some(mut second_order)) =
                                (first_order, second_order)
                            else {
                                stereo.type_ = tagINCHIStereoType0D_INCHI_StereoType_None as i8;
                                stereo.central_atom = NO_ATOM as i16;
                                stereo.neighbor[0] = -1;
                                stereo.neighbor[3] = -1;
                                *error |= 64;
                                add_source_error(heap, error_text, "0D stereobond not recognized")?;
                                continue;
                            };
                            let inchi_num_h2 = |atom: &inchi_Atom| {
                                i32::from(atom.num_iso_H[0].max(0))
                                    + atom.num_iso_H[1..]
                                        .iter()
                                        .map(|value| i32::from(*value))
                                        .sum::<i32>()
                            };
                            let first_middle = parsed_atoms[first].num_bonds == 2
                                && i32::from(parsed_atoms[first].bond_type[0])
                                    + i32::from(parsed_atoms[first].bond_type[1])
                                    == 2 * tagINCHIBondType_INCHI_BOND_TYPE_DOUBLE as i32
                                && inchi_num_h2(&parsed_atoms[first]) == 0
                                && (parsed_atoms[second].num_bonds != 2
                                    || i32::from(parsed_atoms[second].bond_type[0])
                                        + i32::from(parsed_atoms[second].bond_type[1])
                                        != 2 * tagINCHIBondType_INCHI_BOND_TYPE_DOUBLE as i32);
                            let second_middle = parsed_atoms[second].num_bonds == 2
                                && i32::from(parsed_atoms[second].bond_type[0])
                                    + i32::from(parsed_atoms[second].bond_type[1])
                                    == 2 * tagINCHIBondType_INCHI_BOND_TYPE_DOUBLE as i32
                                && inchi_num_h2(&parsed_atoms[second]) == 0
                                && (parsed_atoms[first].num_bonds != 2
                                    || i32::from(parsed_atoms[first].bond_type[0])
                                        + i32::from(parsed_atoms[first].bond_type[1])
                                        != 2 * tagINCHIBondType_INCHI_BOND_TYPE_DOUBLE as i32);
                            if first_middle ^ second_middle {
                                if first_middle {
                                    std::mem::swap(&mut first, &mut second);
                                    first_order = second_order;
                                }
                                let mut chain = [0_i16; 21];
                                let mut length = 0_usize;
                                let mut current = first;
                                let mut next = second;
                                chain[length] = current as i16;
                                length += 1;
                                chain[length] = next as i16;
                                length += 1;
                                let mut previous = current;
                                while length < 20 {
                                    previous = current;
                                    current = next;
                                    if parsed_atoms[current].num_bonds == 2
                                        && i32::from(parsed_atoms[current].bond_type[0])
                                            + i32::from(parsed_atoms[current].bond_type[1])
                                            == 2 * tagINCHIBondType_INCHI_BOND_TYPE_DOUBLE as i32
                                        && inchi_num_h2(&parsed_atoms[current]) == 0
                                    {
                                        let order = usize::from(
                                            parsed_atoms[current].neighbor[0] == previous as i16,
                                        );
                                        next = parsed_atoms[current].neighbor[order] as usize;
                                        chain[length] = next as i16;
                                        length += 1;
                                    } else {
                                        break;
                                    }
                                }
                                let final_order = if length > 2 {
                                    is_in_the_slist(
                                        Some(&parsed_atoms[current].neighbor),
                                        previous as i16,
                                        parsed_atoms[current].num_bonds as i32,
                                    )?
                                } else {
                                    None
                                };
                                if let Some(order) = final_order {
                                    second_order = order;
                                    second = current;
                                    stereo.neighbor[0] =
                                        parsed_atoms[first].neighbor[usize::from(first_order == 0)];
                                    stereo.neighbor[1] = first as i16;
                                    stereo.neighbor[2] = second as i16;
                                    stereo.neighbor[3] = parsed_atoms[second].neighbor
                                        [usize::from(second_order == 0)];
                                    if length % 2 == 1 {
                                        stereo.central_atom = chain[length / 2];
                                        stereo.type_ =
                                            tagINCHIStereoType0D_INCHI_StereoType_Allene as i8;
                                    } else {
                                        stereo.central_atom = NO_ATOM as i16;
                                    }
                                } else {
                                    stereo.type_ = tagINCHIStereoType0D_INCHI_StereoType_None as i8;
                                    stereo.central_atom = NO_ATOM as i16;
                                    stereo.neighbor[0] = -1;
                                    stereo.neighbor[3] = -1;
                                    *error |= 64;
                                    add_source_error(
                                        heap,
                                        error_text,
                                        "Cumulene stereo not recognized (0D)",
                                    )?;
                                }
                            } else {
                                stereo.neighbor[0] =
                                    parsed_atoms[first].neighbor[usize::from(first_order == 0)];
                                stereo.neighbor[3] =
                                    parsed_atoms[second].neighbor[usize::from(second_order == 0)];
                                stereo.central_atom = NO_ATOM as i16;
                            }

                            let second_parity = (i32::from(stereo.parity) >> SB_PARITY_SHFT)
                                & SB_PARITY_MASK as i32;
                            if stereo.type_ != tagINCHIStereoType0D_INCHI_StereoType_None as i8
                                && (AB_MIN_WELL_DEFINED_PARITY as i32
                                    ..=AB_MAX_WELL_DEFINED_PARITY as i32)
                                    .contains(&second_parity)
                            {
                                for endpoint in 0..2 {
                                    let atom_index = if endpoint == 0 {
                                        stereo.neighbor[1]
                                    } else {
                                        stereo.neighbor[2]
                                    } as usize;
                                    let stereo_order = if endpoint == 0 {
                                        first_order
                                    } else {
                                        second_order
                                    };
                                    let mut minimum = atom_count;
                                    for order in 0..parsed_atoms[atom_index].num_bonds as usize {
                                        let neighbor =
                                            parsed_atoms[atom_index].neighbor[order] as usize;
                                        if order != stereo_order
                                            && is_element_a_metal(&parsed_atoms[neighbor].elname)?
                                                == 0
                                        {
                                            minimum = minimum.min(neighbor);
                                        }
                                    }
                                    if minimum < atom_count {
                                        stereo.neighbor[if endpoint == 0 { 0 } else { 3 }] =
                                            minimum as i16;
                                    } else {
                                        add_source_error(
                                            heap,
                                            error_text,
                                            "Cannot find non-metal stereobond neighor (0D)",
                                        )?;
                                    }
                                }
                            }
                        }
                        _ => {}
                    }
                }

                if let Some(flags) = input_atom_flags.as_deref_mut() {
                    *flags |= parsed_flags;
                }
            }

            if let Some(slot) = atom_slot.as_deref_mut() {
                let target = heap.slice_mut(working_atom)?;
                let target = target
                    .get_mut(..parsed_atoms.len())
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                target.clone_from_slice(&parsed_atoms);
                if slot.is_null() {
                    *slot = working_atom;
                }
            }

            if !find_next_only {
                if reuse_caller_atoms {
                    if let Some(slot) = stereo0d_slot.as_deref_mut()
                        && !slot.is_null()
                        && *slot == working_stereo
                        && !parsed_stereo.is_empty()
                    {
                        heap.slice_mut(working_stereo)?
                            .get_mut(..parsed_stereo.len())
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .clone_from_slice(&parsed_stereo);
                    }
                } else {
                    if let Some(count) = stereo0d_count.as_deref_mut() {
                        *count = 0;
                    }
                    if !parsed_stereo.is_empty()
                        && let Some(slot) = stereo0d_slot.as_deref_mut()
                    {
                        heap.slice_mut(working_stereo)?
                            .get_mut(..parsed_stereo.len())
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .clone_from_slice(&parsed_stereo);
                        *slot = working_stereo;
                        if let Some(count) = stereo0d_count.as_deref_mut() {
                            *count = parsed_stereo.len() as i32;
                        }
                    } else if !working_stereo.is_null() {
                        FreeInchi_Stereo0D(heap, Some(&mut working_stereo))?;
                    }
                }
            }

            if atom_slot.is_none() && !working_atom.is_null() {
                inchi_free(heap, working_atom)?;
                working_atom = SourceMutPointer::null();
            }
            if !working_coordinates.is_null() {
                inchi_free(heap, working_coordinates)?;
                working_coordinates = SourceMutPointer::null();
            }
            if input_type == tagInputType_INPUT_INCHI_XML {
                drain_xml_structure!();
            }

            return Ok(atom_count_i32);
        }
    })();

    let line_cleanup = inchi_free(heap, line);
    let next_line_cleanup = inchi_free(heap, next_line);
    let token_cleanup = inchi_free(heap, atom_token);
    let cleanup = line_cleanup
        .and(next_line_cleanup)
        .and(token_cleanup)
        .map(|_| ());
    match operation {
        Ok(status) => cleanup.map(|_| status),
        Err(error) => Err(error),
    }
}

#[cfg(test)]
#[allow(non_snake_case)]
mod tests {
    use super::*;
    use crate::source_types::{
        INCHI_IOS_STRING, INCHI_IOS_TYPE_FILE, INCHI_IOS_TYPE_STRING, STR_ERR_LEN, SourceFile,
        bitWord, inchi_Stereo0D, inp_ATOM,
    };
    use crate::test_support::allocate_source_fixture;
    use serde_json::{Value, json};
    use std::path::Path;
    use std::process::Command;

    #[test]
    fn source_port__inchi_dll_b__createinchiatom__line_117() {
        let mut heap = SourceHeap::default();

        let atoms = CreateInchiAtom(&mut heap, 3).unwrap();
        assert_eq!(
            heap.slice(atoms.as_const()).unwrap(),
            &[
                inchi_Atom::default(),
                inchi_Atom::default(),
                inchi_Atom::default(),
            ]
        );
        assert_eq!(inchi_free(&mut heap, atoms), Ok(()));

        assert_eq!(
            CreateInchiAtom(&mut heap, -1),
            Err(SourceHeapError::AllocationSizeOverflow)
        );
    }

    #[test]
    fn source_port__inchi_dll_b__preparetomakeinchi__line_285() {
        fn call_prepare(
            heap: &mut SourceHeap,
            moltext: SourceConstPointer<i8>,
            options: SourceConstPointer<i8>,
        ) -> (
            Result<i32, SourceHeapError>,
            STRUCT_DATA,
            INPUT_PARMS,
            ORIG_ATOM_DATA,
            [ORIG_ATOM_DATA; INCHI_NUM as usize],
            [SourceMutPointer<PINChI2>; INCHI_NUM as usize],
            [SourceMutPointer<PINChI_Aux2>; INCHI_NUM as usize],
            [INCHI_IOSTREAM; 4],
            INCHI_IOS_STRING,
        ) {
            let mut sd = STRUCT_DATA {
                nErrorCode: 91,
                ..STRUCT_DATA::default()
            };
            let mut ip = INPUT_PARMS::default();
            ip.nInputType = 92;
            let mut original = ORIG_ATOM_DATA {
                num_inp_atoms: 93,
                ..ORIG_ATOM_DATA::default()
            };
            let mut prepared = std::array::from_fn(|index| ORIG_ATOM_DATA {
                num_inp_atoms: 94 + index as i32,
                ..ORIG_ATOM_DATA::default()
            });
            let sentinel = moltext.as_mut().cast::<PINChI2>();
            let aux_sentinel = moltext.as_mut().cast::<PINChI_Aux2>();
            let mut inchi = [sentinel; INCHI_NUM as usize];
            let mut aux = [aux_sentinel; INCHI_NUM as usize];
            let mut streams: [INCHI_IOSTREAM; 4] = std::array::from_fn(|_| {
                let mut stream = INCHI_IOSTREAM::default();
                stream.type_ = 95;
                stream
            });
            let mut strbuf = INCHI_IOS_STRING::default();
            let (pout, remainder) = streams.split_at_mut(1);
            let (plog, remainder) = remainder.split_at_mut(1);
            let (pprb, inp_file) = remainder.split_at_mut(1);
            let result = PrepareToMakeINCHI(
                heap,
                &mut sd,
                &mut ip,
                &mut original,
                &mut prepared,
                &mut inchi,
                &mut aux,
                &mut pout[0],
                &mut plog[0],
                &mut pprb[0],
                &mut inp_file[0],
                moltext,
                options,
                &mut strbuf,
            );
            (
                result, sd, ip, original, prepared, inchi, aux, streams, strbuf,
            )
        }

        let mut heap = SourceHeap::default();
        let moltext = allocate_source_fixture(
            &mut heap,
            b"empty molfile text\n\0"
                .iter()
                .map(|byte| *byte as i8)
                .collect(),
        );
        let options = allocate_source_fixture(
            &mut heap,
            b"-AuxNone\0".iter().map(|byte| *byte as i8).collect(),
        );
        let (result, sd, ip, original, prepared, inchi, aux, streams, strbuf) =
            call_prepare(&mut heap, moltext.as_const(), options.as_const());
        assert_eq!(result, Ok(0));
        assert_eq!(sd, STRUCT_DATA::default());
        assert_eq!(original, ORIG_ATOM_DATA::default());
        assert_eq!(prepared, std::array::from_fn(|_| ORIG_ATOM_DATA::default()));
        assert_eq!(inchi, [SourceMutPointer::null(); INCHI_NUM as usize]);
        assert_eq!(aux, [SourceMutPointer::null(); INCHI_NUM as usize]);
        assert_eq!(ip.nInputType, tagInputType_INPUT_MOLFILE);
        assert_eq!(ip.bNoStructLabels, 1);
        assert!(ip.pSdfLabel.is_null());
        assert!(ip.pSdfValue.is_null());
        assert_ne!(ip.bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO as i32, 0);
        assert_eq!(streams[0].type_, INCHI_IOS_TYPE_STRING as i32);
        assert_eq!(streams[1].type_, INCHI_IOS_TYPE_STRING as i32);
        assert_eq!(streams[2].type_, INCHI_IOS_TYPE_STRING as i32);
        assert_eq!(streams[3].type_, INCHI_IOS_TYPE_STRING as i32);
        assert_eq!(streams[3].s.pStr, moltext);
        assert_eq!(streams[3].s.nPtr, 0);
        assert_eq!(streams[3].s.nUsedLength, 20);
        assert!(streams[3].f.is_null());
        assert!(!strbuf.pStr.is_null());
        assert_eq!(strbuf.nAllocatedLength, INCHI_STRBUF_INITIAL_SIZE as i32);
        assert_eq!(strbuf.nPtr, INCHI_STRBUF_SIZE_INCREMENT as i32);
        inchi_free(&mut heap, strbuf.pStr).unwrap();

        let mut fallback_heap = SourceHeap::default();
        let fallback_moltext = allocate_source_fixture(
            &mut fallback_heap,
            b"x\0".iter().map(|byte| *byte as i8).collect(),
        );
        let fallback_options = allocate_source_fixture(
            &mut fallback_heap,
            b"-AuxNone\0".iter().map(|byte| *byte as i8).collect(),
        );
        fallback_heap.fail_after_allocations(0);
        let (fallback_result, _, fallback_ip, _, _, _, _, _, fallback_strbuf) = call_prepare(
            &mut fallback_heap,
            fallback_moltext.as_const(),
            fallback_options.as_const(),
        );
        assert_eq!(fallback_result, Ok(0));
        assert_eq!(fallback_heap.source_allocation_calls(), 2);
        assert_eq!(
            fallback_ip.bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO as i32,
            0
        );
        inchi_free(&mut fallback_heap, fallback_strbuf.pStr).unwrap();

        let mut no_ram_heap = SourceHeap::default();
        let no_ram_moltext = allocate_source_fixture(
            &mut no_ram_heap,
            b"x\0".iter().map(|byte| *byte as i8).collect(),
        );
        no_ram_heap.fail_after_allocations(0);
        let (no_ram_result, _, _, _, _, _, _, _, no_ram_strbuf) = call_prepare(
            &mut no_ram_heap,
            no_ram_moltext.as_const(),
            SourceConstPointer::null(),
        );
        assert_eq!(no_ram_result, Ok(MOL2INCHI_NO_RAM as i32));
        assert!(no_ram_strbuf.pStr.is_null());

        let mut bad_heap = SourceHeap::default();
        let bad_moltext = allocate_source_fixture(
            &mut bad_heap,
            b"x\0".iter().map(|byte| *byte as i8).collect(),
        );
        let bad_options = allocate_source_fixture(
            &mut bad_heap,
            b"-Key -OUTPUTSDF\0"
                .iter()
                .map(|byte| *byte as i8)
                .collect(),
        );
        bad_heap.trace_source_allocations();
        let (bad_result, _, _, _, _, _, _, _, bad_strbuf) = call_prepare(
            &mut bad_heap,
            bad_moltext.as_const(),
            bad_options.as_const(),
        );
        assert_eq!(bad_result, Ok(MOL2INCHI_BAD_COMMAND_LINE as i32));
        assert!(bad_strbuf.pStr.is_null());
        assert!(bad_heap.live_source_allocation_count() >= 1);
    }

    #[test]
    fn source_port__inchi_dll_b__postmakeinchicleanup__line_391() {
        let mut heap = SourceHeap::default();
        let canonical_bits = allocate_source_fixture(&mut heap, vec![0 as bitWord]);
        let mut canonical_globals = CANON_GLOBALS {
            m_bBit: canonical_bits,
            ..CANON_GLOBALS::default()
        };

        let inchi_rows = allocate_source_fixture(&mut heap, vec![[SourceMutPointer::null(); 2]]);
        let aux_rows = allocate_source_fixture(&mut heap, vec![[SourceMutPointer::null(); 2]]);
        let mut p_inchi = [inchi_rows, SourceMutPointer::null()];
        let mut p_inchi_aux = [aux_rows, SourceMutPointer::null()];
        let mut sd = STRUCT_DATA {
            num_components: [1, 0],
            ..STRUCT_DATA::default()
        };

        let original_atoms = allocate_source_fixture(&mut heap, vec![inp_ATOM::default()]);
        let prepared_atoms = allocate_source_fixture(&mut heap, vec![inp_ATOM::default()]);
        let reconnected_atoms = allocate_source_fixture(&mut heap, vec![inp_ATOM::default()]);
        let mut original = ORIG_ATOM_DATA {
            at: original_atoms,
            num_inp_atoms: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let mut prepared = [
            ORIG_ATOM_DATA {
                at: prepared_atoms,
                num_inp_atoms: 1,
                ..ORIG_ATOM_DATA::default()
            },
            ORIG_ATOM_DATA {
                at: reconnected_atoms,
                num_inp_atoms: 1,
                ..ORIG_ATOM_DATA::default()
            },
        ];

        let output_text = allocate_source_fixture(&mut heap, vec![b'o' as i8, 0]);
        let log_text = allocate_source_fixture(&mut heap, vec![b'l' as i8, 0]);
        let problem_text = allocate_source_fixture(&mut heap, vec![b'p' as i8, 0]);
        let input_text = allocate_source_fixture(&mut heap, vec![b'i' as i8, 0]);
        let string_buffer = allocate_source_fixture(&mut heap, vec![b's' as i8, 0]);
        let moltext = allocate_source_fixture(&mut heap, vec![b'm' as i8, 0]);
        let make_stream = |pointer| INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            s: INCHI_IOS_STRING {
                pStr: pointer,
                nAllocatedLength: 2,
                nUsedLength: 1,
                nPtr: 0,
            },
            ..INCHI_IOSTREAM::default()
        };
        let mut pout = make_stream(output_text);
        let mut plog = make_stream(log_text);
        let mut pprb = make_stream(problem_text);
        let mut inp_file = make_stream(input_text);
        let mut strbuf = INCHI_IOS_STRING {
            pStr: string_buffer,
            nAllocatedLength: 2,
            nUsedLength: 1,
            nPtr: 8,
        };

        let path0 = allocate_source_fixture(&mut heap, vec![b'a' as i8, 0]);
        let path2 = allocate_source_fixture(&mut heap, vec![b'b' as i8, 0]);
        let mut ip = INPUT_PARMS::default();
        ip.path[0] = path0.as_const();
        ip.path[2] = path2.as_const();

        assert_eq!(
            PostMakeINCHICleanup(
                &mut heap,
                &mut canonical_globals,
                &mut sd,
                &mut ip,
                &mut original,
                &mut prepared,
                &mut p_inchi,
                &mut p_inchi_aux,
                &mut pout,
                &mut plog,
                &mut pprb,
                &mut inp_file,
                moltext.as_const(),
                &mut strbuf,
            ),
            Ok(0)
        );

        for pointer in [
            canonical_bits.cast::<i8>(),
            inchi_rows.cast::<i8>(),
            aux_rows.cast::<i8>(),
            original_atoms.cast::<i8>(),
            prepared_atoms.cast::<i8>(),
            reconnected_atoms.cast::<i8>(),
            output_text,
            problem_text,
            string_buffer,
            path0,
            path2,
        ] {
            assert_eq!(
                heap.slice(pointer.as_const()).map(|_| ()),
                Err(SourceHeapError::MissingAllocation)
            );
        }
        assert!(heap.slice(log_text.as_const()).is_ok());
        assert!(heap.slice(input_text.as_const()).is_ok());
        assert!(heap.slice(moltext.as_const()).is_ok());
        assert_eq!(sd.num_components, [0, 0]);
        assert_eq!(p_inchi, [SourceMutPointer::null(); INCHI_NUM as usize]);
        assert_eq!(p_inchi_aux, [SourceMutPointer::null(); INCHI_NUM as usize]);
        assert!(original.at.is_null());
        assert!(prepared[0].at.is_null());
        assert!(prepared[1].at.is_null());
        assert_eq!(pout.s, INCHI_IOS_STRING::default());
        assert!(pout.f.is_null());
        assert_eq!(pout.type_, INCHI_IOS_TYPE_STRING as i32);
        assert_eq!(pprb.s, INCHI_IOS_STRING::default());
        assert!(pprb.f.is_null());
        assert_eq!(pprb.type_, INCHI_IOS_TYPE_STRING as i32);
        assert_eq!(plog.s.pStr, log_text);
        assert_eq!(inp_file.s.pStr, input_text);
        assert_eq!(strbuf, INCHI_IOS_STRING::default());
        assert!(ip.path.iter().all(|path| path.is_null()));
        assert!(canonical_globals.m_bBit.is_null());
    }

    #[test]
    fn source_port__inchi_dll_b__freeinchi_input__line_439() {
        let mut heap = SourceHeap::default();
        let atom = allocate_source_fixture(&mut heap, vec![inchi_Atom::default(); 2]);
        let atom_alias = atom.as_const();
        let stereo = allocate_source_fixture(&mut heap, vec![inchi_Stereo0D::default()]);
        let stereo_alias = stereo.as_const();
        let options = allocate_source_fixture(&mut heap, vec![b'-' as i8, b'S' as i8, 0]);
        let mut input = inchi_Input {
            atom,
            stereo0D: stereo,
            szOptions: options,
            num_atoms: 2,
            num_stereo0D: 1,
        };

        assert_eq!(FreeInchi_Input(&mut heap, &mut input), Ok(()));
        assert_eq!(input, inchi_Input::default());
        assert_eq!(
            heap.slice(atom_alias),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(stereo_alias),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(options.as_const()).unwrap(),
            &[b'-' as i8, b'S' as i8, 0]
        );
        assert_eq!(inchi_free(&mut heap, options), Ok(()));

        let atom = allocate_source_fixture(&mut heap, vec![inchi_Atom::default()]);
        let atom_alias = atom.as_const();
        let stereo = allocate_source_fixture(&mut heap, vec![inchi_Stereo0D::default(); 2]);
        let mut partial_failure_input = inchi_Input {
            atom,
            stereo0D: stereo.offset(1).unwrap(),
            szOptions: SourceMutPointer::null(),
            num_atoms: 1,
            num_stereo0D: 2,
        };
        assert_eq!(
            FreeInchi_Input(&mut heap, &mut partial_failure_input),
            Err(SourceHeapError::FreeOfInteriorPointer)
        );
        assert!(partial_failure_input.atom.is_null());
        assert_eq!(partial_failure_input.stereo0D, stereo.offset(1).unwrap());
        assert_eq!(partial_failure_input.num_atoms, 1);
        assert_eq!(partial_failure_input.num_stereo0D, 2);
        assert_eq!(
            heap.slice(atom_alias),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(heap.slice(stereo.as_const()).unwrap().len(), 2);
        assert_eq!(inchi_free(&mut heap, stereo), Ok(()));
    }

    #[test]
    fn source_port__inchi_dll_b__freeinchi_atom__line_106() {
        let mut heap = SourceHeap::default();
        assert_eq!(FreeInchi_Atom(&mut heap, None), Ok(()));

        let mut null_slot = SourceMutPointer::null();
        assert_eq!(FreeInchi_Atom(&mut heap, Some(&mut null_slot)), Ok(()));
        assert!(null_slot.is_null());

        let mut atom_slot = allocate_source_fixture(
            &mut heap,
            vec![inchi_Atom::default(), inchi_Atom::default()],
        );
        let alias = atom_slot.as_const();
        assert_eq!(FreeInchi_Atom(&mut heap, Some(&mut atom_slot)), Ok(()));
        assert!(atom_slot.is_null());
        assert_eq!(heap.slice(alias), Err(SourceHeapError::MissingAllocation));

        let allocation = allocate_source_fixture(&mut heap, vec![inchi_Atom::default(); 2]);
        let mut interior_slot = allocation.offset(1).unwrap();
        assert_eq!(
            FreeInchi_Atom(&mut heap, Some(&mut interior_slot)),
            Err(SourceHeapError::FreeOfInteriorPointer)
        );
        assert_eq!(interior_slot, allocation.offset(1).unwrap());
        assert_eq!(heap.slice(allocation.as_const()).unwrap().len(), 2);
        assert_eq!(inchi_free(&mut heap, allocation), Ok(()));
    }

    #[test]
    fn source_port__inchi_dll_b__is_in_the_slist__line_531() {
        let path = [7_i16, 11, 13, 17];

        assert_eq!(is_in_the_slist(Some(&path), 7, 4), Ok(Some(0)));
        assert_eq!(is_in_the_slist(Some(&path), 13, 4), Ok(Some(2)));
        assert_eq!(is_in_the_slist(Some(&path), 17, 4), Ok(Some(3)));
        assert_eq!(is_in_the_slist(Some(&path), 19, 4), Ok(None));

        assert_eq!(
            is_in_the_slist(None, 7, 0),
            Ok(None),
            "the zero-length C branch does not dereference pathAtom"
        );
        assert_eq!(
            is_in_the_slist(None, 7, 1),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(
            is_in_the_slist(Some(&path[1..]), 11, -1),
            Ok(Some(0)),
            "a negative C int remains true and may match before decrement"
        );
        assert_eq!(
            is_in_the_slist(Some(&path[3..]), 19, 2),
            Err(SourceHeapError::PointerOutOfBounds),
            "the Rust source model reports the C out-of-bounds precondition violation"
        );
    }

    #[test]
    fn source_port__inchi_dll_b__is_element_a_metal__line_543() {
        fn c_element(symbol: &str) -> Vec<i8> {
            symbol.bytes().map(|byte| byte as i8).chain([0]).collect()
        }

        let source_metals = [
            "K", "V", "Y", "W", "U", "Li", "Be", "Na", "Mg", "Al", "Ca", "Sc", "Ti", "Cr", "Mn",
            "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
            "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
            "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "Re", "Os", "Ir",
            "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra", "Ac", "Th", "Pa", "Np", "Pu",
            "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf",
        ];
        for symbol in source_metals {
            assert_eq!(is_element_a_metal(&c_element(symbol)), Ok(1), "{symbol}");
        }

        for symbol in [
            "H", "B", "C", "N", "O", "F", "Ne", "Si", "P", "S", "Cl", "Ar", "As", "Se", "Br", "Kr",
            "Te", "I", "Xe", "At", "Rn", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
        ] {
            assert_eq!(is_element_a_metal(&c_element(symbol)), Ok(0), "{symbol}");
        }

        for symbol in ["", "k", "li", "LI", "L", "A", "K;", "Lithium"] {
            assert_eq!(is_element_a_metal(&c_element(symbol)), Ok(0), "{symbol}");
        }
        assert_eq!(
            is_element_a_metal(&[b'F' as i8, b'e' as i8]),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(is_element_a_metal(&[-1, 0]), Ok(0));
    }

    #[derive(Debug)]
    struct InchiToAtomFixtureOutput {
        status: i32,
        error: i32,
        error_text: String,
        atoms: Vec<inchi_Atom>,
        stereo: Vec<inchi_Stereo0D>,
        dimensions: i32,
        bonds: i32,
        label: String,
        value: String,
        id: i64,
        flags: INCHI_MODE,
        stereo_count: i32,
        atom_pointer_state: &'static str,
        stereo_pointer_state: &'static str,
        input_position: i64,
    }

    #[derive(Clone, Copy)]
    struct InchiToAtomFixtureConfig {
        do_not_add_h: i32,
        caller_atom_capacity: i32,
        caller_stereo_capacity: i32,
        allocation_failure_ordinal: i32,
        omit_atom_output: bool,
        omit_stereo_output: bool,
        omit_label_output: bool,
        omit_value_output: bool,
        omit_id_output: bool,
        omit_flags_output: bool,
        pass_zero_atom_capacity: bool,
        pass_zero_stereo_count: bool,
        initial_error_code: i32,
    }

    impl Default for InchiToAtomFixtureConfig {
        fn default() -> Self {
            Self {
                do_not_add_h: 0,
                caller_atom_capacity: 0,
                caller_stereo_capacity: 0,
                allocation_failure_ordinal: 0,
                omit_atom_output: false,
                omit_stereo_output: false,
                omit_label_output: false,
                omit_value_output: false,
                omit_id_output: false,
                omit_flags_output: false,
                pass_zero_atom_capacity: false,
                pass_zero_stereo_count: false,
                initial_error_code: 0,
            }
        }
    }

    fn test_c_text(bytes: &[i8]) -> String {
        let bytes = source_c_string(bytes).unwrap();
        String::from_utf8(bytes.iter().map(|byte| *byte as u8).collect()).unwrap()
    }

    fn run_inchi_to_atom_fixture(
        input_text: &[u8],
        input_type: INPUT_TYPE,
        file_mode: bool,
    ) -> InchiToAtomFixtureOutput {
        run_inchi_to_atom_fixture_with_config(
            input_text,
            input_type,
            file_mode,
            InchiToAtomFixtureConfig::default(),
        )
    }

    fn run_inchi_to_atom_fixture_with_config(
        input_text: &[u8],
        input_type: INPUT_TYPE,
        file_mode: bool,
        config: InchiToAtomFixtureConfig,
    ) -> InchiToAtomFixtureOutput {
        let mut heap = SourceHeap::default();
        let mut string_allocation = SourceMutPointer::null();
        let mut file_allocation = SourceMutPointer::null();
        let mut stream = if file_mode {
            file_allocation = allocate_source_fixture(
                &mut heap,
                vec![SourceFile {
                    bytes: input_text.to_vec(),
                    ..SourceFile::default()
                }],
            );
            INCHI_IOSTREAM {
                f: file_allocation,
                type_: INCHI_IOS_TYPE_FILE as i32,
                ..INCHI_IOSTREAM::default()
            }
        } else {
            string_allocation = allocate_source_fixture(
                &mut heap,
                input_text.iter().map(|byte| *byte as i8).collect(),
            );
            INCHI_IOSTREAM {
                s: INCHI_IOS_STRING {
                    pStr: string_allocation,
                    nAllocatedLength: input_text.len() as i32,
                    nUsedLength: input_text.len() as i32,
                    nPtr: 0,
                },
                type_: INCHI_IOS_TYPE_STRING as i32,
                ..INCHI_IOSTREAM::default()
            }
        };
        let label = allocate_source_fixture(&mut heap, vec![0_i8; MAX_SDF_HEADER as usize]);
        let value = allocate_source_fixture(&mut heap, vec![0_i8; MAX_SDF_VALUE as usize]);
        let error_text = allocate_source_fixture(&mut heap, vec![0_i8; STR_ERR_LEN as usize]);
        let mut atoms = if !config.omit_atom_output && config.caller_atom_capacity > 0 {
            allocate_source_fixture(
                &mut heap,
                vec![inchi_Atom::default(); config.caller_atom_capacity as usize],
            )
        } else {
            SourceMutPointer::null()
        };
        let original_atoms = atoms;
        let mut stereo = if !config.omit_stereo_output && config.caller_stereo_capacity > 0 {
            allocate_source_fixture(
                &mut heap,
                vec![inchi_Stereo0D::default(); config.caller_stereo_capacity as usize],
            )
        } else {
            SourceMutPointer::null()
        };
        let original_stereo = stereo;
        let mut stereo_count = if config.pass_zero_stereo_count {
            0
        } else {
            config.caller_stereo_capacity
        };
        let mut dimensions = -1;
        let mut bonds = -1;
        let mut id = -1;
        let mut flags = 0;
        let mut error = config.initial_error_code;

        if config.allocation_failure_ordinal > 0 {
            heap.fail_after_allocations(config.allocation_failure_ordinal as u64 - 1);
        }

        let status = InchiToInchiAtom(
            &mut heap,
            Some(&mut stream),
            (!config.omit_stereo_output).then_some(&mut stereo),
            Some(&mut stereo_count),
            config.do_not_add_h,
            0,
            input_type,
            (!config.omit_atom_output).then_some(&mut atoms),
            if config.pass_zero_atom_capacity {
                0
            } else {
                config.caller_atom_capacity
            },
            Some(&mut dimensions),
            Some(&mut bonds),
            (!config.omit_label_output).then_some(label),
            (!config.omit_value_output).then_some(value),
            (!config.omit_id_output).then_some(&mut id),
            (!config.omit_flags_output).then_some(&mut flags),
            Some(&mut error),
            Some(error_text),
        )
        .unwrap();

        let atom_storage_count =
            if atoms == original_atoms && config.caller_atom_capacity > status.max(0) {
                config.caller_atom_capacity
            } else {
                status.max(0)
            };
        let parsed_atoms = if atom_storage_count > 0 && !atoms.is_null() {
            heap.slice(atoms.as_const()).unwrap()[..atom_storage_count as usize].to_vec()
        } else {
            Vec::new()
        };
        let stereo_storage_count =
            if stereo == original_stereo && config.caller_stereo_capacity > stereo_count.max(0) {
                config.caller_stereo_capacity
            } else {
                stereo_count.max(0)
            };
        let parsed_stereo = if stereo_storage_count > 0 && !stereo.is_null() {
            heap.slice(stereo.as_const()).unwrap()[..stereo_storage_count as usize].to_vec()
        } else {
            Vec::new()
        };
        let atom_pointer_state = if atoms.is_null() {
            "null"
        } else if !original_atoms.is_null() && atoms == original_atoms {
            "reused"
        } else {
            "allocated"
        };
        let stereo_pointer_state = if stereo.is_null() {
            "null"
        } else if !original_stereo.is_null() && stereo == original_stereo {
            "reused"
        } else {
            "allocated"
        };
        let input_position = if file_mode {
            heap.slice(file_allocation.as_const()).unwrap()[0].position as i64
        } else {
            i64::from(stream.s.nPtr)
        };
        let output = InchiToAtomFixtureOutput {
            status,
            error,
            error_text: test_c_text(heap.slice(error_text.as_const()).unwrap()),
            atoms: parsed_atoms,
            stereo: parsed_stereo,
            dimensions,
            bonds,
            label: test_c_text(heap.slice(label.as_const()).unwrap()),
            value: test_c_text(heap.slice(value.as_const()).unwrap()),
            id,
            flags,
            stereo_count,
            atom_pointer_state,
            stereo_pointer_state,
            input_position,
        };

        if !atoms.is_null() {
            inchi_free(&mut heap, atoms).unwrap();
        }
        if !stereo.is_null() {
            inchi_free(&mut heap, stereo).unwrap();
        }
        inchi_free(&mut heap, label).unwrap();
        inchi_free(&mut heap, value).unwrap();
        inchi_free(&mut heap, error_text).unwrap();
        if file_mode {
            inchi_free(&mut heap, file_allocation).unwrap();
        } else {
            inchi_free(&mut heap, string_allocation).unwrap();
        }
        output
    }

    #[test]
    fn source_port__inchi_dll_b__inchitoinchiatom__line_582() {
        let mut heap = SourceHeap::default();
        let mut dimensions = 0;
        let mut bonds = 0;
        let mut error = 0;
        assert_eq!(
            InchiToInchiAtom(
                &mut heap,
                None,
                None,
                None,
                0,
                0,
                tagInputType_INPUT_INCHI_PLAIN,
                None,
                0,
                Some(&mut dimensions),
                Some(&mut bonds),
                None,
                None,
                None,
                None,
                Some(&mut error),
                None,
            ),
            Err(SourceHeapError::NullPointer)
        );

        let labeled = b"Structure: 42. LABEL =value\n\
AuxInfo=1/0/N:2/rA:2cCO/rB:s1;/rC:0,0,0;1,1,1;\n";
        let string_output =
            run_inchi_to_atom_fixture(labeled, tagInputType_INPUT_INCHI_PLAIN, false);
        let file_output = run_inchi_to_atom_fixture(labeled, tagInputType_INPUT_INCHI_PLAIN, true);
        for output in [&string_output, &file_output] {
            assert_eq!(output.status, 2);
            assert_eq!(output.error, 0);
            assert_eq!(output.error_text, "");
            assert_eq!(output.dimensions, 3);
            assert_eq!(output.bonds, 1);
            assert_eq!(output.label, "LABEL");
            assert_eq!(output.value, "value");
            assert_eq!(output.id, 42);
            assert_eq!(output.flags, u64::from(FLAG_INP_AT_CHIRAL));
            assert_eq!(&output.atoms[0].elname[..2], &[b'C' as i8, 0]);
            assert_eq!(&output.atoms[1].elname[..2], &[b'O' as i8, 0]);
            assert_eq!(output.atoms[0].neighbor[0], 1);
            assert_eq!(output.atoms[1].neighbor[0], 0);
            assert_eq!(output.atoms[0].bond_type[0], 1);
            assert_eq!(output.atoms[0].num_iso_H[0], -1);
            assert_eq!(output.atoms[1].num_iso_H[0], -1);
        }

        let undefined = run_inchi_to_atom_fixture(
            b"AuxInfo=1/rA:2CO/rB:s1;/rC:;;\n",
            tagInputType_INPUT_INCHI_PLAIN,
            false,
        );
        assert_eq!((undefined.dimensions, undefined.bonds), (0, 0));
        let two_dimensional = run_inchi_to_atom_fixture(
            b"AuxInfo=1/rA:2CO/rB:s1;/rC:0,0,0;1,0,0;\n",
            tagInputType_INPUT_INCHI_PLAIN,
            false,
        );
        assert_eq!((two_dimensional.dimensions, two_dimensional.bonds), (2, 1));

        let integer_boundaries = run_inchi_to_atom_fixture(
            b"AuxInfo=1/rA:1C4+255.128i32768oh255d2t3/rB:/rC:;\n",
            tagInputType_INPUT_INCHI_PLAIN,
            false,
        );
        assert_eq!(integer_boundaries.status, 1);
        assert_eq!(integer_boundaries.atoms[0].charge, -1);
        assert_eq!(integer_boundaries.atoms[0].radical, -128);
        assert_eq!(integer_boundaries.atoms[0].isotopic_mass, -32768);
        assert_eq!(integer_boundaries.atoms[0].num_iso_H, [4, -1, 2, 3]);
        assert_eq!(integer_boundaries.stereo.len(), 1);
        assert_eq!(
            integer_boundaries.stereo[0].type_,
            tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral as i8
        );
        assert_eq!(
            integer_boundaries.stereo[0].parity,
            tagINCHIStereoParity0D_INCHI_PARITY_ODD as i8
        );

        let bond_cases = [
            (b'v', 1, 4, -4),
            (b'V', 1, -4, 4),
            (b'w', 2, 3, 3),
            (b's', 1, 0, 0),
            (b'd', 2, 0, 0),
            (b't', 3, 0, 0),
            (b'p', 1, 1, -1),
            (b'P', 1, -1, 1),
            (b'n', 1, 6, -6),
            (b'N', 1, -6, 6),
        ];
        for (bond_character, expected_type, stereo1, stereo2) in bond_cases {
            let text = format!(
                "AuxInfo=1/rA:2CO/rB:{}1;/rC:0,0,0;1,0,0;\n",
                bond_character as char
            );
            let output =
                run_inchi_to_atom_fixture(text.as_bytes(), tagInputType_INPUT_INCHI_PLAIN, false);
            assert_eq!(output.status, 2, "{}", bond_character as char);
            assert_eq!(output.atoms[1].bond_type[0], expected_type);
            assert_eq!(output.atoms[1].bond_stereo[0], stereo1);
            assert_eq!(output.atoms[0].bond_stereo[0], stereo2);
        }
        let aromatic_warning = run_inchi_to_atom_fixture(
            b"AuxInfo=1/rA:2CO/rB:a1;/rC:0,0,0;1,0,0;\n",
            tagInputType_INPUT_INCHI_PLAIN,
            false,
        );
        assert_eq!(aromatic_warning.error, 32);
        assert_eq!(
            aromatic_warning.error_text,
            "Atom has 1 or more than 3 aromatic bonds"
        );
        assert_eq!(aromatic_warning.atoms[0].bond_type[0], 1);
        assert_eq!(aromatic_warning.atoms[1].bond_type[0], 1);

        let metal_stereo = run_inchi_to_atom_fixture(
            b"AuxInfo=1/rA:5FeOCCC/rB:;s1s2;d-+3;s4;/rC:;;;;;\n",
            tagInputType_INPUT_INCHI_PLAIN,
            false,
        );
        assert_eq!(metal_stereo.status, 5);
        assert_eq!(metal_stereo.stereo.len(), 1);
        assert_eq!(metal_stereo.stereo[0].parity, 17);
        assert_eq!(metal_stereo.stereo[0].neighbor, [1, 2, 3, 4]);
        assert_eq!(metal_stereo.error_text, "");

        let xml = b"<structure number=\"7\" id.name=\"xml-label\" id.value=\"xml-value\">\n\
<identifier.auxiliary-info>\n\
<reversibility>\n\
<atoms>\n2nCO\n</atoms>\n\
<bonds>\ns1;\n</bonds>\n\
<xyz>\n0,0,0;1,1,1;\n</xyz>\n";
        let xml_output = run_inchi_to_atom_fixture(xml, tagInputType_INPUT_INCHI_XML, false);
        assert_eq!(xml_output.status, 2);
        assert_eq!(xml_output.id, 7);
        assert_eq!(xml_output.label, "xml-label");
        assert_eq!(xml_output.value, "xml-value");
        assert_eq!(xml_output.flags, u64::from(FLAG_INP_AT_NONCHIRAL));
        assert_eq!((xml_output.dimensions, xml_output.bonds), (3, 1));

        let xml_error = run_inchi_to_atom_fixture(
            b"<structure number=\"8\">\n\
<message type=\"error (no InChI)\" value=\"source error\">\n",
            tagInputType_INPUT_INCHI_XML,
            false,
        );
        assert_eq!(xml_error.status, INCHI_INP_ERROR_RET);
        assert_eq!(xml_error.error, INCHI_INP_ERROR_ERR as i32);
        assert_eq!(xml_error.error_text, "source error");
        let xml_fatal = run_inchi_to_atom_fixture(
            b"<structure number=\"9\">\n\
<message type=\"fatal (aborted)\" value=\"source fatal\">\n",
            tagInputType_INPUT_INCHI_XML,
            false,
        );
        assert_eq!(xml_fatal.status, INCHI_INP_FATAL_RET as i32);
        assert_eq!(xml_fatal.error, INCHI_INP_FATAL_ERR as i32);
        assert_eq!(xml_fatal.error_text, "source fatal");

        for (text, expected_message) in [
            (
                b"AuxInfo=1/no-reversibility\n".as_slice(),
                "Missing atom data",
            ),
            (
                b"AuxInfo=1/rA:2CO/rB:x/rC:;;\n".as_slice(),
                "Wrong bonds data",
            ),
            (
                b"AuxInfo=1/rA:2CO/rB:s1;/rC:0,0x;1,0,0;\n".as_slice(),
                "Wrong atom coordinates data",
            ),
        ] {
            let malformed = run_inchi_to_atom_fixture(text, tagInputType_INPUT_INCHI_PLAIN, false);
            assert_eq!(malformed.status, INCHI_INP_ERROR_RET);
            assert_eq!(malformed.error, INCHI_INP_ERROR_ERR as i32);
            assert_eq!(malformed.error_text, expected_message);
            assert!(malformed.atoms.is_empty());
            assert!(malformed.stereo.is_empty());
        }
        let eof = run_inchi_to_atom_fixture(b"", tagInputType_INPUT_INCHI_PLAIN, false);
        assert_eq!(eof.status, INCHI_INP_EOF_RET as i32);
        assert_eq!(eof.error, INCHI_INP_EOF_ERR as i32);

        let mut token_in_final_segment = b"AuxInfo=".to_vec();
        token_in_final_segment.extend(std::iter::repeat_n(b'X', INCHI_LINE_LEN as usize + 32));
        token_in_final_segment.extend_from_slice(b"/rA:1C/rB:/rC:;\n");
        let final_segment = run_inchi_to_atom_fixture(
            &token_in_final_segment,
            tagInputType_INPUT_INCHI_PLAIN,
            false,
        );
        assert_eq!(final_segment.status, INCHI_INP_ERROR_RET);
        assert_eq!(final_segment.error, INCHI_INP_ERROR_ERR as i32);
        assert_eq!(final_segment.error_text, "Missing atom data");

        let mut long_input = b"AuxInfo=".to_vec();
        long_input.extend(std::iter::repeat_n(b'X', INCHI_LINE_LEN as usize + 32));
        long_input.extend_from_slice(b"/rA:1C/rB:/rC:;/");
        long_input.extend(std::iter::repeat_n(b'X', INCHI_LINE_LEN as usize));
        long_input.push(b'\n');
        let cross_buffer =
            run_inchi_to_atom_fixture(&long_input, tagInputType_INPUT_INCHI_PLAIN, false);
        assert_eq!(cross_buffer.status, 1);
        assert_eq!(&cross_buffer.atoms[0].elname[..2], &[b'C' as i8, 0]);

        let mut heap = SourceHeap::default();
        let input_bytes = b"AuxInfo=1/rA:2CO/rB:s1;/rC:0,0,0;1,0,0;\n";
        let input_allocation = allocate_source_fixture(
            &mut heap,
            input_bytes.iter().map(|byte| *byte as i8).collect(),
        );
        let mut stream = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: input_allocation,
                nAllocatedLength: input_bytes.len() as i32,
                nUsedLength: input_bytes.len() as i32,
                nPtr: 0,
            },
            type_: INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        let original_atoms = allocate_source_fixture(&mut heap, vec![inchi_Atom::default(); 3]);
        let mut atoms = original_atoms;
        let mut stereo = SourceMutPointer::null();
        let mut stereo_count = 0;
        let mut dimensions = 0;
        let mut bonds = 0;
        let mut error = 0;
        assert_eq!(
            InchiToInchiAtom(
                &mut heap,
                Some(&mut stream),
                Some(&mut stereo),
                Some(&mut stereo_count),
                0,
                0,
                tagInputType_INPUT_INCHI_PLAIN,
                Some(&mut atoms),
                3,
                Some(&mut dimensions),
                Some(&mut bonds),
                None,
                None,
                None,
                None,
                Some(&mut error),
                None,
            ),
            Ok(2)
        );
        assert_eq!(
            atoms, original_atoms,
            "sufficient caller capacity is reused"
        );
        assert_eq!(heap.slice(atoms.as_const()).unwrap().len(), 3);
        inchi_free(&mut heap, atoms).unwrap();
        inchi_free(&mut heap, input_allocation).unwrap();

        for allocations_before_failure in [0_u64, 1] {
            let mut heap = SourceHeap::default();
            let text = b"AuxInfo=1/rA:1Coe/rB:/rC:;\n";
            let input_allocation =
                allocate_source_fixture(&mut heap, text.iter().map(|byte| *byte as i8).collect());
            let error_buffer = allocate_source_fixture(&mut heap, vec![0_i8; STR_ERR_LEN as usize]);
            let mut stream = INCHI_IOSTREAM {
                s: INCHI_IOS_STRING {
                    pStr: input_allocation,
                    nAllocatedLength: text.len() as i32,
                    nUsedLength: text.len() as i32,
                    nPtr: 0,
                },
                type_: INCHI_IOS_TYPE_STRING as i32,
                ..INCHI_IOSTREAM::default()
            };
            let mut atoms = SourceMutPointer::null();
            let mut stereo = SourceMutPointer::null();
            let mut stereo_count = 0;
            let mut dimensions = 0;
            let mut bonds = 0;
            let mut error = 0;
            heap.fail_after_allocations(allocations_before_failure);
            assert_eq!(
                InchiToInchiAtom(
                    &mut heap,
                    Some(&mut stream),
                    Some(&mut stereo),
                    Some(&mut stereo_count),
                    0,
                    0,
                    tagInputType_INPUT_INCHI_PLAIN,
                    Some(&mut atoms),
                    0,
                    Some(&mut dimensions),
                    Some(&mut bonds),
                    None,
                    None,
                    None,
                    None,
                    Some(&mut error),
                    Some(error_buffer),
                ),
                Ok(INCHI_INP_FATAL_RET as i32)
            );
            assert_eq!(error, INCHI_INP_FATAL_ERR as i32);
            assert_eq!(
                test_c_text(heap.slice(error_buffer.as_const()).unwrap()),
                "Out of RAM"
            );
            assert!(atoms.is_null(), "partially allocated atoms must be cleaned");
            assert!(
                stereo.is_null(),
                "partially allocated stereo must be cleaned"
            );
            assert_eq!(stereo_count, 0);
            inchi_free(&mut heap, input_allocation).unwrap();
            inchi_free(&mut heap, error_buffer).unwrap();
        }

        for (marker, expected_parity) in [
            ('o', tagINCHIStereoParity0D_INCHI_PARITY_ODD as i8),
            ('e', tagINCHIStereoParity0D_INCHI_PARITY_EVEN as i8),
            ('u', tagINCHIStereoParity0D_INCHI_PARITY_UNKNOWN as i8),
            ('?', tagINCHIStereoParity0D_INCHI_PARITY_UNDEFINED as i8),
        ] {
            let text = format!("AuxInfo=1/rA:4CCCCl{marker}/rB:;;s1s2s3;/rC:;;;;\n");
            let output =
                run_inchi_to_atom_fixture(text.as_bytes(), tagInputType_INPUT_INCHI_PLAIN, false);
            assert_eq!(output.status, 4, "tetrahedral parity {marker}");
            assert_eq!(output.stereo.len(), 1, "tetrahedral parity {marker}");
            assert_eq!(output.stereo[0].parity, expected_parity);
            assert_eq!(output.stereo[0].central_atom, 3);
            assert_eq!(output.stereo[0].neighbor, [3, 0, 1, 2]);
        }

        let parity_tokens = [
            ('-', tagINCHIStereoParity0D_INCHI_PARITY_ODD as i32),
            ('+', tagINCHIStereoParity0D_INCHI_PARITY_EVEN as i32),
            ('u', tagINCHIStereoParity0D_INCHI_PARITY_UNKNOWN as i32),
            ('?', tagINCHIStereoParity0D_INCHI_PARITY_UNDEFINED as i32),
        ];
        for (first_token, first_parity) in parity_tokens {
            for (second_token, second_parity) in parity_tokens {
                let text = format!(
                    "AuxInfo=1/rA:5FeOCCC/rB:;s1s2;d{first_token}{second_token}3;s4;/rC:;;;;;\n"
                );
                let output = run_inchi_to_atom_fixture(
                    text.as_bytes(),
                    tagInputType_INPUT_INCHI_PLAIN,
                    false,
                );
                assert_eq!(output.status, 5, "bond parity {first_token}{second_token}");
                assert_eq!(output.stereo.len(), 1);
                assert_eq!(
                    i32::from(output.stereo[0].parity),
                    first_parity | (second_parity << SB_PARITY_SHFT)
                );
            }
        }

        let do_not_add_h = run_inchi_to_atom_fixture_with_config(
            b"AuxInfo=1/rA:1C/rB:/rC:;\n",
            tagInputType_INPUT_INCHI_PLAIN,
            false,
            InchiToAtomFixtureConfig {
                do_not_add_h: 1,
                ..InchiToAtomFixtureConfig::default()
            },
        );
        assert_eq!(do_not_add_h.status, 1);
        assert_eq!(do_not_add_h.atoms[0].num_iso_H[0], 0);

        let caller_capacity = run_inchi_to_atom_fixture_with_config(
            b"AuxInfo=1/rA:2CCle/rB:s1;/rC:0,0,0;1,0,0;\n",
            tagInputType_INPUT_INCHI_PLAIN,
            false,
            InchiToAtomFixtureConfig {
                caller_atom_capacity: 3,
                caller_stereo_capacity: 3,
                ..InchiToAtomFixtureConfig::default()
            },
        );
        assert_eq!(caller_capacity.status, 2);
        assert_eq!(caller_capacity.atom_pointer_state, "reused");
        assert_eq!(caller_capacity.stereo_pointer_state, "reused");
        assert_eq!(caller_capacity.atoms.len(), 3);
        assert_eq!(caller_capacity.atoms[0].num_iso_H[0], 0);
        assert_eq!(caller_capacity.atoms[1].num_iso_H[0], 0);
        assert_eq!(caller_capacity.atoms[2], inchi_Atom::default());
        assert_eq!((caller_capacity.dimensions, caller_capacity.bonds), (0, 0));
        assert_eq!(caller_capacity.flags, 0);
        assert_eq!(caller_capacity.stereo_count, 3);
        assert_eq!(caller_capacity.stereo.len(), 3);
        assert_eq!(caller_capacity.stereo[1], inchi_Stereo0D::default());
        assert_eq!(caller_capacity.stereo[2], inchi_Stereo0D::default());

        for (text, expected_message) in [
            (b"AuxInfo=1/rA:1C/rC:;\n".as_slice(), "Missing bonds data"),
            (
                b"AuxInfo=1/rA:1C/rB:/\n".as_slice(),
                "Missing atom coordinates data",
            ),
            (
                b"AuxInfo=1/rA:2C/rB:;/rC:;;\n".as_slice(),
                "Wrong number of atoms",
            ),
            (
                b"AuxInfo=1/rA:2CO/rB:q1;/rC:;;\n".as_slice(),
                "Wrong bond type",
            ),
            (
                b"AuxInfo=1/rA:2CO/rB:s3;/rC:;;\n".as_slice(),
                "Bond to nonexistent atom",
            ),
            (
                b"AuxInfo=1/rA:3CCO/rB:s1/rC:;;;\n".as_slice(),
                "Wrong number of bonds",
            ),
            (
                b"AuxInfo=1/rA:2CO/rB:s1;/rC:;;;/tail\n".as_slice(),
                "Wrong number of coordinates",
            ),
        ] {
            let output = run_inchi_to_atom_fixture(text, tagInputType_INPUT_INCHI_PLAIN, false);
            assert_eq!(output.status, INCHI_INP_ERROR_RET, "{expected_message}");
            assert_eq!(output.error, INCHI_INP_ERROR_ERR as i32);
            assert_eq!(output.error_text, expected_message);
        }

        let trailing_coordinate_data = run_inchi_to_atom_fixture(
            b"AuxInfo=1/rA:2CO/rB:s1;/rC:;;;\n",
            tagInputType_INPUT_INCHI_PLAIN,
            false,
        );
        assert_eq!(trailing_coordinate_data.status, 2);
        assert_eq!(trailing_coordinate_data.error, 0);

        let xml_missing_reversibility = run_inchi_to_atom_fixture(
            b"<structure number=\"10\">\n\
<identifier.auxiliary-info>\n\
</identifier.auxiliary-info>\n",
            tagInputType_INPUT_INCHI_XML,
            false,
        );
        assert_eq!(xml_missing_reversibility.status, INCHI_INP_ERROR_RET);
        assert_eq!(
            xml_missing_reversibility.error_text,
            "Missing reversibility info"
        );
    }

    fn official_c_oracle_record(
        case_id: &str,
        text: &[u8],
        input_type: INPUT_TYPE,
        file_mode: bool,
        config: InchiToAtomFixtureConfig,
    ) -> Value {
        let output = run_inchi_to_atom_fixture_with_config(text, input_type, file_mode, config);
        let atom_count = output.status.max(0) as usize;
        let atoms = output
            .atoms
            .iter()
            .map(|atom| {
                json!({
                    "x_bits": atom.x.to_bits(),
                    "y_bits": atom.y.to_bits(),
                    "z_bits": atom.z.to_bits(),
                    "neighbor": atom.neighbor,
                    "bond_type": atom.bond_type,
                    "bond_stereo": atom.bond_stereo,
                    "elname_bytes": atom.elname,
                    "num_bonds": atom.num_bonds,
                    "num_iso_h": atom.num_iso_H,
                    "isotopic_mass": atom.isotopic_mass,
                    "radical": atom.radical,
                    "charge": atom.charge,
                })
            })
            .collect::<Vec<_>>();
        let bond_fields = output
            .atoms
            .iter()
            .take(atom_count)
            .enumerate()
            .flat_map(|(atom_index, atom)| {
                let num_bonds = usize::try_from(atom.num_bonds)
                    .expect("official oracle atom bond count must be nonnegative");
                (0..num_bonds).map(move |slot| {
                    json!({
                        "atom_index": atom_index,
                        "slot": slot,
                        "neighbor": atom.neighbor[slot],
                        "bond_type": atom.bond_type[slot],
                        "bond_stereo": atom.bond_stereo[slot],
                    })
                })
            })
            .collect::<Vec<_>>();
        let stereo0d = output
            .stereo
            .iter()
            .map(|stereo| {
                json!({
                    "central_atom": stereo.central_atom,
                    "neighbor": stereo.neighbor,
                    "type": stereo.type_,
                    "parity": stereo.parity,
                })
            })
            .collect::<Vec<_>>();
        json!({
            "schema_version": "cosmolkit-inchi-official-c-v1",
            "oracle": {
                "project": "IUPAC InChI",
                "tag": "v1.07.5",
                "commit": "11a87982bb518f57ac013f0b258c283655e1ea1d",
                "tree_sha256": "4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd",
                "api_version": "1.07.5",
            },
            "case_id": case_id,
            "operation": "inchi_to_inchi_atom",
            "input": {
                "mode": if file_mode { "file" } else { "string" },
                "input_type": input_type,
                "text": String::from_utf8(text.to_vec()).unwrap(),
                "do_not_add_h": config.do_not_add_h,
                "caller_atom_capacity": config.caller_atom_capacity,
                "caller_stereo_capacity": config.caller_stereo_capacity,
                "allocation_failure_ordinal": config.allocation_failure_ordinal,
                "omit_atom_output": config.omit_atom_output as i32,
                "omit_stereo_output": config.omit_stereo_output as i32,
                "omit_label_output": config.omit_label_output as i32,
                "omit_value_output": config.omit_value_output as i32,
                "omit_id_output": config.omit_id_output as i32,
                "omit_flags_output": config.omit_flags_output as i32,
                "pass_zero_atom_capacity": config.pass_zero_atom_capacity as i32,
                "pass_zero_stereo_count": config.pass_zero_stereo_count as i32,
                "initial_error_code": config.initial_error_code,
            },
            "output": {
                "status": output.status,
                "error_code": output.error,
                "error_text": output.error_text,
                "atom_count": atom_count,
                "bond_count": output.bonds,
                "coordinate_dimensions": output.dimensions,
                "atoms": atoms,
                "bond_fields": bond_fields,
                "stereo0d": stereo0d,
                "sdf_label": (!config.omit_label_output).then_some(output.label),
                "sdf_value": (!config.omit_value_output).then_some(output.value),
                "id": (!config.omit_id_output).then_some(output.id),
                "atom_flags": (!config.omit_flags_output).then_some(output.flags),
                "stereo_count": output.stereo_count,
                "atom_pointer_state": output.atom_pointer_state,
                "stereo_pointer_state": output.stereo_pointer_state,
                "input_position": output.input_position,
            },
        })
    }

    #[test]
    fn official_c_oracle_inchi_to_inchi_atom_exact() {
        let expected_case_ids = [
            "plain-labeled-string",
            "plain-labeled-file",
            "plain-integer-boundaries",
            "plain-metal-stereo",
            "plain-aromatic-warning",
            "xml-labeled",
            "xml-error",
            "plain-malformed-coordinates",
            "plain-eof",
            "xml-fatal",
            "plain-missing-atom",
            "plain-wrong-bonds",
            "plain-undefined-coordinates",
            "plain-2d-coordinates",
            "plain-bond-v",
            "plain-bond-V",
            "plain-bond-w",
            "plain-bond-t",
            "plain-bond-p",
            "plain-bond-P",
            "plain-bond-n",
            "plain-bond-N",
            "plain-tetrahedral-odd",
            "plain-tetrahedral-even",
            "plain-tetrahedral-unknown",
            "plain-tetrahedral-undefined",
            "plain-parity-odd-odd",
            "plain-parity-odd-even",
            "plain-parity-odd-unknown",
            "plain-parity-odd-undefined",
            "plain-parity-even-odd",
            "plain-parity-even-even",
            "plain-parity-even-unknown",
            "plain-parity-even-undefined",
            "plain-parity-unknown-odd",
            "plain-parity-unknown-even",
            "plain-parity-unknown-unknown",
            "plain-parity-unknown-undefined",
            "plain-parity-undefined-odd",
            "plain-parity-undefined-even",
            "plain-parity-undefined-unknown",
            "plain-parity-undefined-undefined",
            "plain-do-not-add-h",
            "plain-caller-capacity",
            "plain-allocation-failure-atom",
            "plain-allocation-failure-stereo",
            "plain-cross-buffer",
            "plain-token-final-segment",
            "plain-empty-reversibility",
            "plain-nonempty-double-slash",
            "plain-missing-bonds",
            "plain-missing-coordinates",
            "plain-bonds-no-section-terminator",
            "plain-single-atom-trailing-bond-data",
            "plain-wrong-atom-count",
            "plain-wrong-bond-type",
            "plain-nonexistent-bond",
            "plain-wrong-bond-count",
            "plain-wrong-coordinate-count",
            "plain-too-few-coordinates-section-terminator",
            "plain-trailing-coordinate-data",
            "plain-long-integer-narrowing",
            "xml-missing-reversibility",
            "xml-reversibility-without-auxinfo",
            "plain-find-next-only",
            "plain-no-stereo-output",
            "plain-error-no-stereo-output",
            "plain-no-metadata-outputs",
            "plain-caller-atom-too-small",
            "plain-caller-stereo-too-small",
            "plain-header-no-dot-space",
            "plain-header-is-missing",
            "plain-nonchiral-isolated",
            "plain-empty-atom-quantities",
            "plain-aromatic-two",
            "plain-aromatic-three",
            "plain-aromatic-four",
            "plain-allene",
            "plain-cumulene",
            "xml-eof-after-reversibility",
            "xml-wrong-atoms-tag",
            "xml-eof-after-atoms-tag",
            "xml-no-structure-number",
            "xml-caller-capacities-too-small",
            "plain-zero-atom-count",
            "plain-bond-nonalphabetic",
            "plain-stereo-reallocation",
            "plain-allocation-failure-stereo-reallocation",
            "plain-allocation-failure-coordinate",
            "plain-partial-coordinate-components",
            "plain-empty-coordinate-components",
            "xml-rich-atom",
            "xml-all-bond-codes",
            "xml-tetrahedral",
            "xml-allene",
            "xml-aromatic-two",
            "xml-aromatic-three",
            "xml-aromatic-four",
            "xml-caller-capacities-reused",
            "xml-caller-capacities-reused-error",
            "xml-no-stereo-output",
            "xml-metal-stereo",
            "xml-find-next-only",
            "xml-allocation-failure-atom",
            "xml-allocation-failure-stereo",
            "xml-allocation-failure-coordinate",
            "xml-zero-atom-count",
            "xml-wrong-atom-count",
            "xml-bond-nonalphabetic",
            "xml-missing-bond-neighbor",
            "xml-nonexistent-bond",
            "xml-wrong-bond-count",
            "xml-missing-coordinate-tag",
            "xml-malformed-coordinate",
            "xml-nested-error-cleanup",
            "xml-success-structure-cleanup",
            "plain-header-empty-label-negative-charge",
            "plain-header-no-label",
            "xml-atom-parity-even",
            "xml-atom-parity-unknown",
            "xml-atom-parity-undefined",
            "xml-tetrahedral-three-neighbors-valid",
            "xml-tetrahedral-four-neighbors-valid",
            "xml-bond-parity-even-odd",
            "xml-bond-parity-unknown-undefined",
            "xml-wrong-bond-type",
            "xml-empty-bonds-direct-end",
            "xml-empty-coordinates-direct-end",
            "xml-coordinate-leading-commas",
            "xml-stereo-reallocation",
            "xml-allocation-failure-stereo-reallocation",
            "xml-eof-before-bonds-header",
            "xml-eof-after-bonds-header",
            "xml-atom-trailing-data",
            "xml-bond-trailing-data",
            "xml-coordinate-trailing-data",
            "xml-cross-line-sections",
            "xml-bond-parity-undefined-unknown",
            "xml-explicit-valence",
            "plain-caller-error-cleanup",
            "plain-reverse-allene",
            "xml-even-cumulene",
            "xml-wrong-bonds-header-eof",
            "xml-wrong-bonds-header-next-structure",
            "xml-structure-end-as-atom-data",
            "xml-nested-structure-as-atom-data",
            "plain-caller-zero-passed-capacities",
            "plain-no-header-no-metadata-outputs",
            "plain-empty-structure-header",
            "plain-empty-label-value",
            "plain-missing-label-no-output",
            "plain-skip-unrecognized-line",
            "plain-empty-reversibility-no-newline",
            "plain-positive-atom-count-no-payload",
            "plain-extra-atom",
            "plain-secondary-parity-default",
            "plain-tetrahedral-four-neighbors-valid",
            "xml-no-metadata-outputs",
            "xml-number-not-terminated-at-endptr",
            "xml-id-name-unterminated",
            "xml-id-value-unterminated",
            "xml-error-find-next-only",
            "xml-error-unterminated",
            "xml-atom-optional-numeric-suffixes",
            "xml-secondary-parity-default",
            "xml-wrong-bonds-header-initial-error",
            "plain-auxinfo-no-slash",
            "plain-mixed-aromatic-single",
            "plain-cumulene-max-chain",
            "plain-nonzero-identical-bond-coordinates",
            "xml-negative-charge",
            "xml-cumulene-max-chain",
            "xml-eof-before-coordinate-header",
            "xml-eof-after-coordinate-header",
            "xml-caller-capacity-reused-valid",
            "xml-nonzero-identical-bond-coordinates",
            "xml-cross-buffer-all-sections",
            "xml-eof-mid-atoms",
            "xml-eof-mid-bonds",
            "xml-eof-mid-coordinates",
            "xml-long-open-atoms",
            "xml-long-open-bonds",
            "xml-long-open-coordinates",
            "plain-long-open-atoms",
            "plain-long-open-bonds",
            "plain-long-open-coordinates",
            "plain-eof-initial-error",
            "plain-stereo-simple-endpoints",
            "plain-stereo-degree-two-single-double",
            "plain-stereo-middle-isotope-h",
            "plain-stereo-middle-isotope-h-reverse",
            "plain-stereo-later-interior-isotope-h",
            "plain-stereo-central-cumulene",
            "xml-stereo-simple-endpoints",
            "xml-stereo-degree-two-single-double",
            "xml-stereo-middle-isotope-h",
            "xml-stereo-middle-isotope-h-reverse",
            "xml-stereo-later-interior-isotope-h",
            "xml-stereo-central-cumulene",
            "xml-success-file-mode",
            "invalid-input-type",
        ];

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/official_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--inchi-to-inchi-atom-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let official_records = String::from_utf8(oracle.stdout)
            .expect("official C oracle output must be UTF-8")
            .lines()
            .map(|line| serde_json::from_str::<Value>(line).expect("oracle record must be JSON"))
            .collect::<Vec<_>>();
        let official_case_ids = official_records
            .iter()
            .map(|record| record["case_id"].as_str().expect("case_id must be text"))
            .collect::<Vec<_>>();
        assert_eq!(official_case_ids, expected_case_ids);

        for official in &official_records {
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let input = &official["input"];
            let text = input["text"].as_str().expect("input text must be text");
            let input_type = input["input_type"]
                .as_i64()
                .expect("input_type must be signed integer")
                as INPUT_TYPE;
            let file_mode = input["mode"].as_str() == Some("file");
            let config = InchiToAtomFixtureConfig {
                do_not_add_h: input["do_not_add_h"]
                    .as_i64()
                    .expect("do_not_add_h must be integer") as i32,
                caller_atom_capacity: input["caller_atom_capacity"]
                    .as_i64()
                    .expect("caller_atom_capacity must be integer")
                    as i32,
                caller_stereo_capacity: input["caller_stereo_capacity"]
                    .as_i64()
                    .expect("caller_stereo_capacity must be integer")
                    as i32,
                allocation_failure_ordinal: input["allocation_failure_ordinal"]
                    .as_i64()
                    .expect("allocation_failure_ordinal must be integer")
                    as i32,
                omit_atom_output: input["omit_atom_output"]
                    .as_i64()
                    .expect("omit_atom_output must be integer")
                    != 0,
                omit_stereo_output: input["omit_stereo_output"]
                    .as_i64()
                    .expect("omit_stereo_output must be integer")
                    != 0,
                omit_label_output: input["omit_label_output"]
                    .as_i64()
                    .expect("omit_label_output must be integer")
                    != 0,
                omit_value_output: input["omit_value_output"]
                    .as_i64()
                    .expect("omit_value_output must be integer")
                    != 0,
                omit_id_output: input["omit_id_output"]
                    .as_i64()
                    .expect("omit_id_output must be integer")
                    != 0,
                omit_flags_output: input["omit_flags_output"]
                    .as_i64()
                    .expect("omit_flags_output must be integer")
                    != 0,
                pass_zero_atom_capacity: input["pass_zero_atom_capacity"]
                    .as_i64()
                    .expect("pass_zero_atom_capacity must be integer")
                    != 0,
                pass_zero_stereo_count: input["pass_zero_stereo_count"]
                    .as_i64()
                    .expect("pass_zero_stereo_count must be integer")
                    != 0,
                initial_error_code: input["initial_error_code"]
                    .as_i64()
                    .expect("initial_error_code must be integer")
                    as i32,
            };
            let rust =
                official_c_oracle_record(case_id, text.as_bytes(), input_type, file_mode, config);
            assert_eq!(official, &rust, "exact field mismatch for {case_id}");
        }
    }
}
