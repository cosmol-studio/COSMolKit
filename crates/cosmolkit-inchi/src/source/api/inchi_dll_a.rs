use crate::source::api::inchi_dll::{
    ExtractOneStructure, input_erroneously_contains_pseudoatoms, parse_options_string,
};
use crate::source::api::inchi_dll_a2::{CanonOneStructureINChI, NormOneStructureINChI, make_norm_atoms_from_inp_atoms};
use crate::source::base::ichi_io::{
    inchi_ios_close, inchi_ios_init, inchi_strbuf_close, inchi_strbuf_init, inchi_strbuf_reset,
};
use crate::source::base::ichican2::SetBitFree;
use crate::source::base::ichierr::AddErrorMessage;
use crate::source::base::ichinorm::FreeInpAtomData;
use crate::source::base::ichiparm::{HelpCommandLineParms, InchiBuildMetadata, PrintInputParms, ReadCommandLineParms};
use crate::source::base::ichiprt1::{OrigStruct_FillOut, OrigStruct_Free};
use crate::source::base::ichitaut::free_t_group_info;
use crate::source::base::mol_fmt4::OrigAtData_WriteToSDfile;
use crate::source::base::mol2atom::{FreeCompAtomData, FreeOrigAtData};
use crate::source::base::runichi4::{FreeAllINChIArrays, SortAndPrintINChI, TreatCreateINChIWarning, bIsStructChiral};
use crate::source::base::util::{inchi_free, inchi_malloc, inchi_stricmp};
use crate::source_types::{
    _IS_EOF, _IS_ERROR, _IS_FATAL, _IS_OKAY, _IS_SKIP, _IS_UNKNOWN, _IS_WARNING, CANON_GLOBALS, FILE,
    FLAG_INP_AT_CHIRAL, FLAG_NORM_CONSIDER_TAUT, FLAG_PROTON_CHARGE_CANCEL, INCHI_BAS, INCHI_CLOCK,
    INCHI_IOS_TYPE_FILE, INCHI_IOS_TYPE_STRING, INCHI_MAX_NUM_ARG, INCHI_NUM, INCHI_OPTION_PREFX,
    INCHI_OUT_NO_AUX_INFO, INCHI_OUT_PLAIN_TEXT, INCHI_OUT_SAVEOPT, INCHI_OUT_SDFILE_ATOMS_DT, INCHI_OUT_SDFILE_ONLY,
    INCHI_OUT_SHORT_AUX_INFO, INCHI_OUT_STDINCHI, INCHI_REC, INCHI_STRBUF_INITIAL_SIZE, INCHI_STRBUF_SIZE_INCREMENT,
    INCHIGEN_CONTROL, INCHIGEN_DATA, INCHIGEN_HANDLE, INP_ATOM_DATA, INPUT_PARMS, MAX_NUM_PATHS, MAX_SDF_VALUE,
    ORIG_ATOM_DATA, ORIG_STRUCT, REQ_MODE_BASIC, REQ_MODE_CHIR_FLG_STEREO, REQ_MODE_DIFF_UU_STEREO,
    REQ_MODE_RACEMIC_STEREO, REQ_MODE_RELATIVE_STEREO, REQ_MODE_SB_IGN_ALL_UU, REQ_MODE_SC_IGN_ALL_UU, REQ_MODE_STEREO,
    SAVE_OPT_15T, SAVE_OPT_FIXEDH, SAVE_OPT_KET, SAVE_OPT_RECMET, SAVE_OPT_SLUUD, SAVE_OPT_SUU, STRUCT_DATA,
    SourceArgvPointer, SourceConstPointer, SourceHeap, SourceHeapError, SourceMutPointer, TAUT_NON, TAUT_NUM, TAUT_YES,
    TG_FLAG_1_5_TAUT, TG_FLAG_DISCONNECT_COORD_DONE, TG_FLAG_KETO_ENOL_TAUT, TG_FLAG_RECONNECT_COORD, bRELEASE_VERSION,
    inchi_Input, inchi_InputEx, inchi_Output, tagRetValGetINCHI_inchi_Ret_EOF, tagRetValGetINCHI_inchi_Ret_ERROR,
    tagRetValGetINCHI_inchi_Ret_FATAL, tagRetValGetINCHI_inchi_Ret_OKAY, tagRetValGetINCHI_inchi_Ret_SKIP,
    tagRetValGetINCHI_inchi_Ret_UNKNOWN, tagRetValGetINCHI_inchi_Ret_WARNING,
};

fn free_inp_atom_data_array(
    heap: &mut SourceHeap,
    pointer: SourceMutPointer<INP_ATOM_DATA>,
    count: i32,
) -> Result<(), SourceHeapError> {
    for index in 0..count {
        let item_pointer = pointer.offset(i64::from(index))?;
        let mut item = heap
            .slice(item_pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        FreeInpAtomData(heap, &mut item)?;
        *heap
            .slice_mut(item_pointer)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)? = item;
    }
    inchi_free(heap, pointer)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn STDINCHIGEN_Setup(
    heap: &mut SourceHeap,
    h_gen: INCHIGEN_HANDLE,
    mut p_gen_data: Option<&mut INCHIGEN_DATA>,
    p_inp: Option<&inchi_Input>,
    stdout: SourceMutPointer<FILE>,
    build: InchiBuildMetadata<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:185 STDINCHIGEN_Setup
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int INCHI_DECL STDINCHIGEN_Setup( INCHIGEN_HANDLE _HGen,
                                 INCHIGEN_DATA * pGenData,
                                 inchi_Input * pInp )
{
    INCHIGEN_CONTROL *HGen = (INCHIGEN_CONTROL *) _HGen;
    INPUT_PARMS *ip = &( HGen->InpParms );
    STRUCT_DATA *sd = &( HGen->StructData );
    int retcode = inchi_Ret_OKAY;

    retcode = INCHIGEN_Setup( _HGen, pGenData, pInp );

    /* Ensure standardness */
    if (ip->bINChIOutputOptions & INCHI_OUT_SAVEOPT)
    {
        ip->bINChIOutputOptions &= ~INCHI_OUT_SAVEOPT;
        retcode = _IS_WARNING;
    }
    if (0 != ( ip->bTautFlags & TG_FLAG_RECONNECT_COORD ))
    {
        ip->bTautFlags &= ~TG_FLAG_RECONNECT_COORD;
        retcode = _IS_WARNING;
    }
    if (0 != ( ip->nMode & REQ_MODE_BASIC ))
    {
        ip->nMode &= ~REQ_MODE_BASIC;
        retcode = _IS_WARNING;
    }
    if (0 != ( ip->nMode & REQ_MODE_RELATIVE_STEREO ))
    {
        ip->nMode &= ~( REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO );
        retcode = _IS_WARNING;
    }
    if (0 != ( ip->nMode & REQ_MODE_RACEMIC_STEREO ))
    {
        ip->nMode &= ~( REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO );
        retcode = _IS_WARNING;
    }
    if (0 != ( ip->nMode & REQ_MODE_CHIR_FLG_STEREO ))
    {
        ip->nMode &= ~( REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO );
        retcode = _IS_WARNING;
    }
    if (0 != ( ip->nMode & REQ_MODE_DIFF_UU_STEREO ))
    {
        ip->nMode &= ~REQ_MODE_DIFF_UU_STEREO;
        retcode = _IS_WARNING;
    }
    if (0 == ( ip->nMode & ( REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU ) ))
    {
        ip->nMode |= REQ_MODE_SB_IGN_ALL_UU;
        ip->nMode |= REQ_MODE_SC_IGN_ALL_UU;
        retcode = _IS_WARNING;
    }
    if (0 != ( ip->bTautFlags & TG_FLAG_KETO_ENOL_TAUT ))
    {
        ip->bTautFlags &= ~TG_FLAG_KETO_ENOL_TAUT;
        retcode = _IS_WARNING;
    }
    if (0 != ( ip->bTautFlags & TG_FLAG_1_5_TAUT ))
    {
        ip->bTautFlags &= ~TG_FLAG_1_5_TAUT;
        retcode = _IS_WARNING;
    }

    /* And anyway... */
    ip->bINChIOutputOptions |= INCHI_OUT_STDINCHI;
    ip->bINChIOutputOptions &= ~INCHI_OUT_SAVEOPT;

    strcpy( pGenData->pStrErrStruct, sd->pStrErrStruct );

    return retcode;
}
    */
    // END INCHI C FUNCTION: STDINCHIGEN_Setup
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: STDINCHIGEN_Setup
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; no function-local conditional compilation.
    // INCHI✔️❌: SourceHeap checks and cloned generator state are materially slower than native direct field access.
    // END INCHI ACTIVE MACRO CONFIGURATION: STDINCHIGEN_Setup

    let mut retcode = INCHIGEN_Setup(
        heap,
        h_gen,
        p_gen_data.as_deref_mut(),
        p_inp,
        stdout,
        build,
    )?;
    let control_pointer = h_gen.cast::<INCHIGEN_CONTROL>();
    let mut control = heap
        .slice(control_pointer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let input_parameters = &mut control.InpParms;
    if input_parameters.bINChIOutputOptions & INCHI_OUT_SAVEOPT as i32 != 0 {
        input_parameters.bINChIOutputOptions &= !(INCHI_OUT_SAVEOPT as i32);
        retcode = _IS_WARNING as i32;
    }
    if input_parameters.bTautFlags & u64::from(TG_FLAG_RECONNECT_COORD) != 0 {
        input_parameters.bTautFlags &= !u64::from(TG_FLAG_RECONNECT_COORD);
        retcode = _IS_WARNING as i32;
    }
    if input_parameters.nMode & u64::from(REQ_MODE_BASIC) != 0 {
        input_parameters.nMode &= !u64::from(REQ_MODE_BASIC);
        retcode = _IS_WARNING as i32;
    }
    let stereo_modes = u64::from(
        REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO,
    );
    if input_parameters.nMode & u64::from(REQ_MODE_RELATIVE_STEREO) != 0 {
        input_parameters.nMode &= !stereo_modes;
        retcode = _IS_WARNING as i32;
    }
    if input_parameters.nMode & u64::from(REQ_MODE_RACEMIC_STEREO) != 0 {
        input_parameters.nMode &= !stereo_modes;
        retcode = _IS_WARNING as i32;
    }
    if input_parameters.nMode & u64::from(REQ_MODE_CHIR_FLG_STEREO) != 0 {
        input_parameters.nMode &= !stereo_modes;
        retcode = _IS_WARNING as i32;
    }
    if input_parameters.nMode & u64::from(REQ_MODE_DIFF_UU_STEREO) != 0 {
        input_parameters.nMode &= !u64::from(REQ_MODE_DIFF_UU_STEREO);
        retcode = _IS_WARNING as i32;
    }
    if input_parameters.nMode & u64::from(REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU) == 0 {
        input_parameters.nMode |=
            u64::from(REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU);
        retcode = _IS_WARNING as i32;
    }
    if input_parameters.bTautFlags & u64::from(TG_FLAG_KETO_ENOL_TAUT) != 0 {
        input_parameters.bTautFlags &= !u64::from(TG_FLAG_KETO_ENOL_TAUT);
        retcode = _IS_WARNING as i32;
    }
    if input_parameters.bTautFlags & u64::from(TG_FLAG_1_5_TAUT) != 0 {
        input_parameters.bTautFlags &= !u64::from(TG_FLAG_1_5_TAUT);
        retcode = _IS_WARNING as i32;
    }
    input_parameters.bINChIOutputOptions |= INCHI_OUT_STDINCHI as i32;
    input_parameters.bINChIOutputOptions &= !(INCHI_OUT_SAVEOPT as i32);

    if let Some(data) = p_gen_data {
        let nul = control
            .StructData
            .pStrErrStruct
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        data.pStrErrStruct[..=nul]
            .copy_from_slice(&control.StructData.pStrErrStruct[..=nul]);
    } else {
        *heap
            .slice_mut(control_pointer)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)? = control;
        return Err(SourceHeapError::NullPointer);
    }
    *heap
        .slice_mut(control_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = control;
    Ok(retcode)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn INCHIGEN_Setup(
    heap: &mut SourceHeap,
    h_gen: INCHIGEN_HANDLE,
    mut p_gen_data: Option<&mut INCHIGEN_DATA>,
    p_inp: Option<&inchi_Input>,
    stdout: SourceMutPointer<FILE>,
    build: InchiBuildMetadata<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:260 INCHIGEN_Setup
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int INCHI_DECL INCHIGEN_Setup( INCHIGEN_HANDLE _HGen,
                              INCHIGEN_DATA * pGenData,
                              inchi_Input * pInp )
{
    int retcode = inchi_Ret_OKAY;

    INCHIGEN_CONTROL *HGen = (INCHIGEN_CONTROL *) _HGen;

    ORIG_ATOM_DATA *orig_inp_data = &( HGen->OrigInpData );
    STRUCT_DATA *sd = &( HGen->StructData );
    INPUT_PARMS *ip = &( HGen->InpParms );
    INCHI_IOSTREAM *log_file = HGen->inchi_file + 1;
    INCHI_IOSTREAM prbstr, *prb_file = &prbstr;

    const char *argv[INCHI_MAX_NUM_ARG + 1];
    int   argc;
    char *szOptions = NULL;
    char szSdfDataValue[MAX_SDF_VALUE + 1];
    int bReleaseVersion = bRELEASE_VERSION;
    unsigned long  ulDisplTime = 0;    /*  infinite, milliseconds */
    int p;
    inchi_InputEx inpEx;
    inchi_InputEx *pInpEx = &inpEx;

    /* No '*' or 'Zz' elements are allowed in the input . */
    if ( input_erroneously_contains_pseudoatoms( pInp, NULL) )
    {
        AddErrorMessage(sd->pStrErrStruct, "Pseudoatoms are not supported in current API mode");
        sd->nStructReadError = 99;
        sd->nErrorType = _IS_ERROR;
        retcode = _IS_ERROR;
        goto ret;
    }

    pInpEx->atom = pInp->atom;
    pInpEx->num_atoms = pInp->num_atoms;
    pInpEx->num_stereo0D = pInp->num_stereo0D;
    pInpEx->stereo0D = pInp->stereo0D;
    pInpEx->szOptions = pInp->szOptions;
    pInpEx->polymer = NULL; /* v. 1.05 ext not supported in modularized API */
    pInpEx->v3000 = NULL; /* v. 1.05 ext not supported in modularized API */


    /* Allocate/init */
    if (!pGenData)
    {
        retcode = _IS_ERROR;
        goto ret;
    }
    memset( pGenData, 0, sizeof( *pGenData ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    /* Parse 'command-line' options and fill internal INPUT_PARMS structure */
    if (pInp && pInp->szOptions)
    {
        szOptions = (char*) inchi_malloc( strlen( pInp->szOptions ) + 1 );
        if (!szOptions)
            return _IS_FATAL;   /* Not enough memory.... */
        else
        {
            /* Parse. */
            strcpy( szOptions, pInp->szOptions );
            argc = parse_options_string( szOptions, argv, INCHI_MAX_NUM_ARG );
        }
    }
    else
    {
        /* Got NULL options string or NULL 'pInp', will use defaults. */
        argc = 1;
        argv[0] = "";
        argv[1] = NULL;
    }


    if ((argc == 1
#ifdef TARGET_API_LIB
        && ( !pInp || pInp->num_atoms <= 0 || !pInp->atom ))
#endif
        || (argc == 2 && ( argv[1][0] == INCHI_OPTION_PREFX ) &&
        ( !strcmp( argv[1] + 1, "?" ) || !inchi_stricmp( argv[1] + 1, "help" ) ))) /* djb-rwth: addressing LLVM warnings */
    {

        HelpCommandLineParms( log_file );
        memset( log_file, 0, sizeof( *log_file ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        return _IS_EOF;
    }

    memset( szSdfDataValue, 0, sizeof( szSdfDataValue ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    /* Decrypt command line    */
    retcode = ReadCommandLineParms( argc, argv, ip, szSdfDataValue, &ulDisplTime, bReleaseVersion, log_file );

    if (szOptions)
    {
        inchi_free( szOptions );
    }

    /* INChI DLL specific */
    ip->bNoStructLabels = 1;

    if (0 > retcode)
    {
        goto ret;
    }

    if (ip->bNoStructLabels)
    {
        ip->pSdfLabel = NULL;
        ip->pSdfValue = NULL;
    }
    else
    {
        if (ip->nInputType == INPUT_INCHI_XML ||
          ip->nInputType == INPUT_INCHI_PLAIN ||
          ip->nInputType == INPUT_CMLFILE)
        {
            /* the input may contain both the header and the label of the structure */
            if (!ip->pSdfLabel)
            {
                ip->pSdfLabel = ip->szSdfDataHeader;
            }
            if (!ip->pSdfValue)
            {
                ip->pSdfValue = szSdfDataValue;
            }
        }
    }

    if (retcode != inchi_Ret_OKAY) goto ret;

    PrintInputParms( log_file, ip );

    /* Extract the structure */
    retcode = ExtractOneStructure( sd, ip, HGen->szTitle, pInpEx, log_file,
                                   HGen->inchi_file, /* out_file */
                                   prb_file, orig_inp_data,
                                   &( HGen->num_inp ) );

ret:switch (retcode)
{
    case _IS_OKAY: retcode = inchi_Ret_OKAY; HGen->init_passed = 1; break;    /* Success; break; no errors or warnings */

    case _IS_ERROR: ( HGen->num_err )++;  retcode = inchi_Ret_ERROR; break;
                                                            /* Error: no INChI has been created */
    case _IS_FATAL: ( HGen->num_err )++;  retcode = inchi_Ret_FATAL; break;
                                                            /* Severe error: no INChI has been created
                                                            (typically; break; memory allocation failed) */
    case _IS_SKIP: retcode = inchi_Ret_SKIP; break;   /* not used in INChI dll */
    case _IS_EOF: retcode = inchi_Ret_EOF; break;   /* no structural data has been provided */
    case _IS_WARNING: retcode = inchi_Ret_WARNING; HGen->init_passed = 1; break;    /* Success; break; warning(s) issued */
    case _IS_UNKNOWN:
    default: retcode = inchi_Ret_UNKNOWN; break;   /* Unlnown program error */
}

    if (NULL!=pGenData)
    {
        strcpy( pGenData->pStrErrStruct, sd->pStrErrStruct );
        for (p = 0; p < INCHI_NUM; p++)
        {
            pGenData->num_components[p] = sd->num_components[p];
        }
    }

    return retcode;
}
    */
    // END INCHI C FUNCTION: INCHIGEN_Setup
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: INCHIGEN_Setup
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; the TARGET_API_LIB no-input term is active.
    // INCHI✔️❌: SourceHeap checks, cloned generator state, and temporary argument vectors are materially slower than native stack arrays and direct pointers.
    // END INCHI ACTIVE MACRO CONFIGURATION: INCHIGEN_Setup

    let control_pointer = h_gen.cast::<INCHIGEN_CONTROL>();
    let mut control = heap
        .slice(control_pointer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut retcode = tagRetValGetINCHI_inchi_Ret_OKAY;

    'source: {
        let input = p_inp.ok_or(SourceHeapError::NullPointer)?;
        if input_erroneously_contains_pseudoatoms(heap, Some(input), None)? != 0 {
            let message = b"Pseudoatoms are not supported in current API mode\0"
                .iter()
                .map(|byte| *byte as i8)
                .collect::<Vec<_>>();
            AddErrorMessage(
                Some(&mut control.StructData.pStrErrStruct),
                Some(&message),
            )?;
            control.StructData.nStructReadError = 99;
            control.StructData.nErrorType = _IS_ERROR as i32;
            retcode = _IS_ERROR as i32;
            break 'source;
        }

        let extended_input = inchi_InputEx {
            atom: input.atom,
            num_atoms: input.num_atoms,
            num_stereo0D: input.num_stereo0D,
            stereo0D: input.stereo0D,
            szOptions: input.szOptions,
            polymer: SourceMutPointer::null(),
            v3000: SourceMutPointer::null(),
        };

        let Some(data) = p_gen_data.as_deref_mut() else {
            retcode = _IS_ERROR as i32;
            break 'source;
        };
        *data = INCHIGEN_DATA::default();

        let mut parsed_arguments = vec![SourceArgvPointer::Null; INCHI_MAX_NUM_ARG as usize + 1];
        let empty_argument = heap.allocate_model_storage(vec![0_i8])?;
        let mut options = SourceMutPointer::null();
        let argument_count = if input.szOptions.is_null() {
            parsed_arguments[0] = SourceArgvPointer::EmptyLiteral;
            1
        } else {
            let source = heap.slice(input.szOptions.as_const())?;
            let length = source
                .iter()
                .position(|byte| *byte == 0)
                .ok_or(SourceHeapError::MissingNulTerminator)?;
            let source = source[..=length].to_vec();
            options = match inchi_malloc(heap, (length + 1) as u64) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => {
                    heap.free(empty_argument)?;
                    *heap
                        .slice_mut(control_pointer)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = control;
                    return Ok(_IS_FATAL as i32);
                }
                Err(error) => return Err(error),
            };
            heap.slice_mut(options)?[..=length].copy_from_slice(&source);
            parse_options_string(
                heap,
                options,
                &mut parsed_arguments,
                INCHI_MAX_NUM_ARG as i32,
            )?
        };
        let mut argv = Vec::with_capacity(argument_count as usize);
        for argument in parsed_arguments.iter().copied().take(argument_count as usize) {
            argv.push(match argument {
                SourceArgvPointer::EmptyLiteral => empty_argument.as_const(),
                SourceArgvPointer::Command(pointer) => pointer.as_const(),
                SourceArgvPointer::Null => return Err(SourceHeapError::NullPointer),
            });
        }

        let no_input = input.num_atoms <= 0 || input.atom.is_null();
        let help_requested = if argument_count == 2 {
            let argument = argv[1];
            let bytes = heap.slice(argument)?;
            if bytes.first().copied() == Some(INCHI_OPTION_PREFX as i8) {
                let suffix = argument.offset(1)?;
                let suffix_bytes = heap.slice(suffix)?;
                let is_question = suffix_bytes.starts_with(&[b'?' as i8, 0]);
                let help = heap.allocate_model_storage(
                    b"help\0".iter().map(|byte| *byte as i8).collect(),
                )?;
                let is_help = inchi_stricmp(heap, suffix, help.as_const())? == 0;
                heap.free(help)?;
                is_question || is_help
            } else {
                false
            }
        } else {
            false
        };
        if (argument_count == 1 && no_input) || help_requested {
            HelpCommandLineParms(heap, Some(&mut control.inchi_file[1]), stdout, build)?;
            control.inchi_file[1] = Default::default();
            heap.free(empty_argument)?;
            *heap
                .slice_mut(control_pointer)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)? = control;
            return Ok(_IS_EOF);
        }

        let sdf_value = heap.allocate_model_storage(vec![0_i8; MAX_SDF_VALUE as usize + 1])?;
        let mut display_time = 0_u64;
        retcode = ReadCommandLineParms(
            heap,
            argument_count,
            &argv,
            &mut control.InpParms,
            sdf_value,
            &mut display_time,
            bRELEASE_VERSION as i32,
            Some(&mut control.inchi_file[1]),
        )?;
        heap.free(sdf_value)?;
        if !options.is_null() {
            inchi_free(heap, options)?;
        }
        heap.free(empty_argument)?;
        control.InpParms.bNoStructLabels = 1;
        if retcode < 0 {
            break 'source;
        }
        control.InpParms.pSdfLabel = SourceMutPointer::null();
        control.InpParms.pSdfValue = SourceMutPointer::null();
        if retcode != tagRetValGetINCHI_inchi_Ret_OKAY {
            break 'source;
        }

        PrintInputParms(heap, Some(&mut control.inchi_file[1]), &control.InpParms)?;
        let mut problem_stream = Default::default();
        let (output_streams, remaining_streams) = control.inchi_file.split_at_mut(1);
        retcode = ExtractOneStructure(
            heap,
            &mut control.StructData,
            &mut control.InpParms,
            None,
            Some(&extended_input),
            Some(&mut remaining_streams[0]),
            Some(&mut output_streams[0]),
            Some(&mut problem_stream),
            &mut control.OrigInpData,
            &mut control.num_inp,
        )?;
    }

    retcode = match retcode {
        value if value == _IS_OKAY as i32 => {
            control.init_passed = 1;
            tagRetValGetINCHI_inchi_Ret_OKAY
        }
        value if value == _IS_ERROR as i32 => {
            control.num_err = control.num_err.wrapping_add(1);
            tagRetValGetINCHI_inchi_Ret_ERROR
        }
        value if value == _IS_FATAL as i32 => {
            control.num_err = control.num_err.wrapping_add(1);
            tagRetValGetINCHI_inchi_Ret_FATAL
        }
        _IS_SKIP => tagRetValGetINCHI_inchi_Ret_SKIP,
        _IS_EOF => tagRetValGetINCHI_inchi_Ret_EOF,
        value if value == _IS_WARNING as i32 => {
            control.init_passed = 1;
            tagRetValGetINCHI_inchi_Ret_WARNING
        }
        value if value == _IS_UNKNOWN as i32 => tagRetValGetINCHI_inchi_Ret_UNKNOWN,
        _ => tagRetValGetINCHI_inchi_Ret_UNKNOWN,
    };

    if let Some(data) = p_gen_data {
        let nul = control
            .StructData
            .pStrErrStruct
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        data.pStrErrStruct[..=nul]
            .copy_from_slice(&control.StructData.pStrErrStruct[..=nul]);
        for index in 0..INCHI_NUM as usize {
            data.num_components[index] = control.StructData.num_components[index];
        }
    }
    *heap
        .slice_mut(control_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = control;
    Ok(retcode)
}

#[allow(non_snake_case)]
pub(crate) fn STDINCHIGEN_DoNormalization(
    heap: &mut SourceHeap,
    h_gen: INCHIGEN_HANDLE,
    p_gen_data: Option<&mut INCHIGEN_DATA>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:435 STDINCHIGEN_DoNormalization
    // INCHI✔️❌: int INCHI_DECL STDINCHIGEN_DoNormalization( INCHIGEN_HANDLE HGen, INCHIGEN_DATA * pGenData )
    // INCHI✔️❌: {
    // INCHI✔️❌:     return INCHIGEN_DoNormalization( HGen, pGenData );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: STDINCHIGEN_DoNormalization
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: STDINCHIGEN_DoNormalization
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; INCHI_DECL has no behavioral body.
    // INCHI✔️❌: Performance inherits the known SourceHeap overhead of INCHIGEN_DoNormalization.
    // END INCHI ACTIVE MACRO CONFIGURATION: STDINCHIGEN_DoNormalization

    INCHIGEN_DoNormalization(heap, h_gen, p_gen_data)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn INCHIGEN_DoNormalization(
    heap: &mut SourceHeap,
    h_gen: INCHIGEN_HANDLE,
    p_gen_data: Option<&mut INCHIGEN_DATA>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:442 INCHIGEN_DoNormalization
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int INCHI_DECL INCHIGEN_DoNormalization( INCHIGEN_HANDLE _HGen, INCHIGEN_DATA *pGenData )
{
    int nRet = 0, nRet1 = 0;
    /* int maxINChI=0; */

    INCHIGEN_CONTROL * HGen = (INCHIGEN_CONTROL *) _HGen;
    INPUT_PARMS *ip = &( HGen->InpParms );
    STRUCT_DATA *sd = &( HGen->StructData );
    NORM_CANON_FLAGS *pncFlags = &( HGen->ncFlags );
    INCHI_IOSTREAM *out_file = HGen->inchi_file;
    INCHI_IOSTREAM inpstr, *inp_file = &inpstr;
    ORIG_ATOM_DATA *orig_inp_data = &( HGen->OrigInpData );
    ORIG_STRUCT      *pOrigStruct = NULL;
    INCHI_CLOCK ic;
    CANON_GLOBALS CG;
    int k;

    memset( &CG, 0, sizeof( CG ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( &ic, 0, sizeof( ic ) ); /* djb-rwth: memset_s C11/Annex K variant? */

#if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )
    int ret1 = 0, ret2 = 0; /* djb-rwth: ignoring LLVM warning: variables used */
#endif

/* Set debug output */
#if (TRACE_MEMORY_LEAKS == 1)

    _CrtSetDbgFlag( _CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF );
/* for execution outside the VC++ debugger uncomment one of the following two */
#ifdef MY_REPORT_FILE
    _CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_WARN, MY_REPORT_FILE );
    _CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_ERROR, MY_REPORT_FILE );
    _CrtSetReportMode( _CRT_ASSERT, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_ASSERT, MY_REPORT_FILE );
#else
    _CrtSetReportMode( _CRT_WARN | _CRT_ERROR, _CRTDBG_MODE_DEBUG );
#endif
#if ( !defined(__STDC__) || __STDC__ != 1 )
    /* turn on floating point exceptions */
    {
        /* Get the default control word. */
        int cw = _controlfp( 0, 0 );
        /* Set the exception masks OFF, turn exceptions on. */
        /*cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_INEXACT|EM_ZERODIVIDE|EM_DENORMAL);*/
        cw &= ~( EM_OVERFLOW | EM_UNDERFLOW | EM_ZERODIVIDE | EM_DENORMAL );
        /* Set the control word. */
        _controlfp( cw, MCW_EM );
    }
#endif /* ( !defined(__STDC__) || __STDC__ != 1 ) */

#endif /* (TRACE_MEMORY_LEAKS == 1) */

    if (HGen->init_passed == 0)
    {
        AddErrorMessage( sd->pStrErrStruct, "InChI generator not initialized" );
        sd->nStructReadError = 99;
        sd->nErrorType = _IS_ERROR;
        nRet = _IS_ERROR;
        goto exit_function;
    }

    inchi_ios_init( inp_file, INCHI_IOS_TYPE_FILE, NULL );

    sd->bUserQuitComponent = 0;
    sd->bUserQuitComponentDisplay = 0;
    memset( HGen->composite_norm_data, 0, sizeof( HGen->composite_norm_data ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( pncFlags, 0, sizeof( *pncFlags ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    /* for testing only */
#if( REMOVE_ION_PAIRS_ORIG_STRU == 1 )
    fix_odd_things( orig_inp_data->num_inp_atoms, orig_inp_data->at, 0 );
#endif

#if 0 /*** FOR NOW, DISABLE UNDERIVATIZATION IN  MODULAR INTERFACE LIBINCHI ***/
#if( UNDERIVATIZE == 1 )  /***** post v.1 feature *****/
    /*if (ip->bUnderivatize && 0 > ( ret2 = underivatize( orig_inp_data ) ))*/

    if (ip->bUnderivatize
         && (ret2 = OAD_Edit_Underivatize( &ic, &CG, orig_inp_data, ( ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY ), ip->bUnderivatize & 2, ip->pSdfValue ) )
        )

    {
        long num_inp2 = HGen->num_inp;
        AddErrorMessage( sd->pStrErrStruct, "Underivatization error" );
        sd->nStructReadError = 99;
        sd->nErrorType = _IS_ERROR;
        nRet = _IS_ERROR;
        TreatReadTheStructureErrors( sd, ip, LOG_MASK_ALL, inp_file, log_file, out_file, prb_file,
                                        prep_inp_data, &num_inp2, HGen->pStr, PSTR_BUFFER_SIZE );
        goto exit_function; /* output only if derivatives found */
    }
#endif /* UNDERIVATIZE == 1 */
#if( RING2CHAIN == 1 )  /***** post v.1 feature *****/
    if (ip->bRing2Chain && 0 > ( ret1 = Ring2Chain( orig_inp_data ) ))
    {
        long num_inp2 = HGen->num_inp;
        AddErrorMessage( sd->pStrErrStruct, "Ring to chain error" );
        sd->nStructReadError = 99;
        sd->nErrorType = _IS_ERROR;
        nRet = _IS_ERROR;
        TreatReadTheStructureErrors( sd, ip, LOG_MASK_ALL, inp_file, log_file, out_file, prb_file,
                                        prep_inp_data, &num_inp2, HGen->pStr, PSTR_BUFFER_SIZE );
        goto exit_function; /* output only if derivatives found */
    }
#endif /* RING2CHAIN == 1 */
#if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )  /***** post v.1 feature *****/
    if (ip->bIngnoreUnchanged && !ret1 && !ret2)
    {
        goto exit_function; /* output only if derivatives or ring/chain found */
    }
#endif /* RING2CHAIN == 1 || UNDERIVATIZE == 1 */
#endif /*** FOR NOW, DISABLE UNDERIVATIZATION IN  MODULAR INTERFACE LIBINCHI ***/


    /***** output MOLfile ***************/
    if (ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY)
    {
        char szNumber[32];
        int ret1a = 0, ret2a = 0; /* for derivatives and ring-chain */ /* djb-rwth: ignoring LLVM warning: variables used to store function return values */
        ret1a = sprintf( szNumber, "Structure #%ld", HGen->num_inp );
        ret2a = OrigAtData_WriteToSDfile( orig_inp_data, out_file, szNumber, NULL,
            ( sd->bChiralFlag & FLAG_INP_AT_CHIRAL ) ? 1 : 0,
            ( ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ATOMS_DT ) ? 1 : 0, ip->pSdfLabel, ip->pSdfValue );
        goto exit_function;
    }

    /******* create full reversibility information **************/
    if (!( ip->bINChIOutputOptions & ( INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO ) ))
    {
        pOrigStruct = &( HGen->OrigStruct );
        memset( pOrigStruct, 0, sizeof( *pOrigStruct ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        if (OrigStruct_FillOut( &CG, orig_inp_data, pOrigStruct, sd ))
        {
            AddErrorMessage( sd->pStrErrStruct, "Cannot interpret reversibility information" );
            sd->nStructReadError = 99;
            sd->nErrorType = _IS_ERROR;
            nRet = _IS_ERROR;
        }
    }

    sd->bUserQuit = 0;
    if (sd->bUserQuit)
    {
        goto exit_function;
    }

    /* Normalize the whole disconnected or original structure */
    if (nRet != _IS_FATAL && nRet != _IS_ERROR)
    {
        nRet1 = NormOneStructureINChI( &CG, &ic, pGenData, HGen, INCHI_BAS, inp_file );
        nRet = inchi_max( nRet, nRet1 );
    }

    /*
    if ( nRet != _IS_FATAL && nRet != _IS_ERROR )
        maxINChI = 1;
    */

    if (nRet != _IS_FATAL && nRet != _IS_ERROR &&
        ( sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE ) &&
        ( ip->bTautFlags               & TG_FLAG_RECONNECT_COORD ))
    {
        /* Normalize  the whole reconnected structure */
        nRet1 = NormOneStructureINChI( &CG, &ic, pGenData, HGen, INCHI_REC, inp_file );
        nRet = inchi_max( nRet, nRet1 );
        /*
        if ( nRet != _IS_FATAL && nRet != _IS_ERROR )
                maxINChI = 2;
        */
    }

exit_function:

    if (nRet != _IS_FATAL && nRet != _IS_ERROR)
    {
        HGen->norm_passed = 1;
    }

    for (k = 0; k < INCHI_NUM; k++)
    {
        pGenData->num_components[k] = sd->num_components[k];
    }

    /* Emit normalization warnings */
    if (nRet != _IS_FATAL && nRet != _IS_ERROR)
    {
        int ics, istruct, itaut, nc[2];
        int warn_prot = 0, warn_neutr = 0;
        INP_ATOM_DATA *inp_norm_data[TAUT_NUM]; /*  = { &InpNormAtData, &InpNormTautData }; */
        nc[0] = pGenData->num_components[0];
        nc[1] = pGenData->num_components[1];
        for (istruct = 0; istruct < 2; istruct++)
        {
            if (nc[istruct] > 0)
            {
                for (ics = 0; ics < nc[istruct]; ics++)
                {
                    inp_norm_data[0] = &( HGen->InpNormAtData[istruct][ics] );
                    inp_norm_data[1] = &( HGen->InpNormTautData[istruct][ics] );
                    for (itaut = 0; itaut < 2; itaut++)
                    {
                        if (NULL != inp_norm_data[itaut])
                        {
                            if (inp_norm_data[itaut]->bTautomeric)
                            {
                                if (inp_norm_data[itaut]->bNormalizationFlags & ( FLAG_NORM_CONSIDER_TAUT &~FLAG_PROTON_CHARGE_CANCEL ))
                                {
                                    if (warn_prot == 0)
                                    {
                                        warn_prot++;
                                        WarningMessage( sd->pStrErrStruct, "Proton(s) added/removed" );
                                    }
                                }
                                if (inp_norm_data[itaut]->bNormalizationFlags & FLAG_PROTON_CHARGE_CANCEL)
                                {
                                    if (warn_neutr == 0)
                                    {
                                        warn_neutr++;
                                        WarningMessage( sd->pStrErrStruct, "Charges neutralized" );
                                    }
                                }
                            }
                        }
                    } /* itaut */
                }
            }
        }
    }

    strcpy( pGenData->pStrErrStruct, sd->pStrErrStruct );
    make_norm_atoms_from_inp_atoms( pGenData, HGen );

    return nRet;
}
    */
    // END INCHI C FUNCTION: INCHIGEN_DoNormalization
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: INCHIGEN_DoNormalization
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64.
    // INCHI✔️❌: TRACE_MEMORY_LEAKS == 0 and REMOVE_ION_PAIRS_ORIG_STRU == 0.
    // INCHI✔️❌: The outer #if 0 excludes all UNDERIVATIZE and RING2CHAIN code.
    // INCHI✔️❌: WarningMessage is the active AddErrorMessage macro alias.
    // INCHI✔️❌: Typed SourceHeap state cloning and shallow-copy conversion retain known overhead.
    // END INCHI ACTIVE MACRO CONFIGURATION: INCHIGEN_DoNormalization

    let control_pointer = h_gen.cast::<INCHIGEN_CONTROL>();
    let mut control = heap
        .slice(control_pointer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let generation_data = p_gen_data.ok_or(SourceHeapError::NullPointer)?;
    let clock = heap.allocate_model_storage(vec![INCHI_CLOCK::default()])?;
    let mut canonical_globals = CANON_GLOBALS::default();

    let operation = (|| -> Result<i32, SourceHeapError> {
        let mut result = 0_i32;
        let mut component_result;
        let mut input_stream = Default::default();

        'normalization: {
            if control.init_passed == 0 {
                let message = b"InChI generator not initialized\0".map(|byte| byte as i8);
                AddErrorMessage(Some(&mut control.StructData.pStrErrStruct), Some(&message))?;
                control.StructData.nStructReadError = 99;
                control.StructData.nErrorType = _IS_ERROR as i32;
                result = _IS_ERROR as i32;
                break 'normalization;
            }

            inchi_ios_init(
                Some(&mut input_stream),
                INCHI_IOS_TYPE_FILE as i32,
                SourceMutPointer::null(),
            )?;
            control.StructData.bUserQuitComponent = 0;
            control.StructData.bUserQuitComponentDisplay = 0;
            control.composite_norm_data =
                std::array::from_fn(|_| std::array::from_fn(|_| Default::default()));
            control.ncFlags = Default::default();

            if control.InpParms.bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY as i32 != 0 {
                let title = format!("Structure #{}", control.num_inp);
                if title.len() >= 32 {
                    return Err(SourceHeapError::PointerOutOfBounds);
                }
                let title = heap.allocate_model_storage(
                    title
                        .bytes()
                        .chain(std::iter::once(0))
                        .map(|byte| byte as i8)
                        .collect(),
                )?;
                let chiral = i32::from(
                    control.StructData.bChiralFlag & FLAG_INP_AT_CHIRAL as i32 != 0,
                );
                let atoms_dt = i32::from(
                    control.InpParms.bINChIOutputOptions & INCHI_OUT_SDFILE_ATOMS_DT as i32 != 0,
                );
                let label = control.InpParms.pSdfLabel.as_const();
                let value = control.InpParms.pSdfValue.as_const();
                let write_result = OrigAtData_WriteToSDfile(
                    heap,
                    &control.OrigInpData,
                    Some(&mut control.inchi_file[0]),
                    SourceMutPointer::null(),
                    title.as_const(),
                    SourceConstPointer::null(),
                    chiral,
                    atoms_dt,
                    label,
                    value,
                );
                heap.free(title)?;
                let _ = write_result?;
                break 'normalization;
            }

            if control.InpParms.bINChIOutputOptions
                & (INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO) as i32
                == 0
            {
                control.OrigStruct = ORIG_STRUCT::default();
                if OrigStruct_FillOut(
                    heap,
                    &mut canonical_globals,
                    &mut control.OrigInpData,
                    &mut control.OrigStruct,
                    &mut control.StructData,
                )? != 0
                {
                    let message =
                        b"Cannot interpret reversibility information\0".map(|byte| byte as i8);
                    AddErrorMessage(
                        Some(&mut control.StructData.pStrErrStruct),
                        Some(&message),
                    )?;
                    control.StructData.nStructReadError = 99;
                    control.StructData.nErrorType = _IS_ERROR as i32;
                    result = _IS_ERROR as i32;
                }
            }

            control.StructData.bUserQuit = 0;
            if control.StructData.bUserQuit != 0 {
                break 'normalization;
            }

            if result != _IS_FATAL as i32 && result != _IS_ERROR as i32 {
                component_result = NormOneStructureINChI(
                    heap,
                    &mut canonical_globals,
                    clock,
                    generation_data,
                    &mut control,
                    INCHI_BAS as i32,
                    Some(&mut input_stream),
                    0,
                )?;
                result = result.max(component_result);
            }

            if result != _IS_FATAL as i32
                && result != _IS_ERROR as i32
                && control.StructData.bTautFlagsDone[INCHI_BAS as usize]
                    & u64::from(TG_FLAG_DISCONNECT_COORD_DONE)
                    != 0
                && control.InpParms.bTautFlags & u64::from(TG_FLAG_RECONNECT_COORD) != 0
            {
                component_result = NormOneStructureINChI(
                    heap,
                    &mut canonical_globals,
                    clock,
                    generation_data,
                    &mut control,
                    INCHI_REC as i32,
                    Some(&mut input_stream),
                    0,
                )?;
                result = result.max(component_result);
            }
        }

        if result != _IS_FATAL as i32 && result != _IS_ERROR as i32 {
            control.norm_passed = 1;
        }
        for k in 0..INCHI_NUM as usize {
            generation_data.num_components[k] = control.StructData.num_components[k];
        }

        if result != _IS_FATAL as i32 && result != _IS_ERROR as i32 {
            let mut warned_protons = false;
            let mut warned_neutralization = false;
            let proton_mask =
                u64::from(FLAG_NORM_CONSIDER_TAUT & !FLAG_PROTON_CHARGE_CANCEL);
            for structure in 0..INCHI_NUM as usize {
                let count = usize::try_from(generation_data.num_components[structure])
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                for component in 0..count {
                    for source in [
                        control.InpNormAtData[structure],
                        control.InpNormTautData[structure],
                    ] {
                        let normalized = heap
                            .slice(source.as_const())?
                            .get(component)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if normalized.bTautomeric != 0 {
                            if normalized.bNormalizationFlags & proton_mask != 0
                                && !warned_protons
                            {
                                warned_protons = true;
                                let message =
                                    b"Proton(s) added/removed\0".map(|byte| byte as i8);
                                AddErrorMessage(
                                    Some(&mut control.StructData.pStrErrStruct),
                                    Some(&message),
                                )?;
                            }
                            if normalized.bNormalizationFlags
                                & u64::from(FLAG_PROTON_CHARGE_CANCEL)
                                != 0
                                && !warned_neutralization
                            {
                                warned_neutralization = true;
                                let message = b"Charges neutralized\0".map(|byte| byte as i8);
                                AddErrorMessage(
                                    Some(&mut control.StructData.pStrErrStruct),
                                    Some(&message),
                                )?;
                            }
                        }
                    }
                }
            }
        }

        let message_length = control
            .StructData
            .pStrErrStruct
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        generation_data.pStrErrStruct[..=message_length]
            .copy_from_slice(&control.StructData.pStrErrStruct[..=message_length]);
        make_norm_atoms_from_inp_atoms(heap, generation_data, &control)?;

        Ok(result)
    })();

    let free_result = heap.free(clock);
    *heap
        .slice_mut(control_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = control;
    let result = operation?;
    free_result?;
    Ok(result)
}

#[allow(non_snake_case)]
pub(crate) fn STDINCHIGEN_DoCanonicalization(
    heap: &mut SourceHeap,
    h_gen: INCHIGEN_HANDLE,
    p_gen_data: Option<&mut INCHIGEN_DATA>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:689 STDINCHIGEN_DoCanonicalization
    // INCHI✔️❌: int INCHI_DECL STDINCHIGEN_DoCanonicalization
    // INCHI✔️❌: ( INCHIGEN_HANDLE HGen, INCHIGEN_DATA * pGenData )
    // INCHI✔️❌: {
    // INCHI✔️❌:     return INCHIGEN_DoCanonicalization( HGen, pGenData );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: STDINCHIGEN_DoCanonicalization
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: STDINCHIGEN_DoCanonicalization
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; INCHI_DECL has no behavioral body.
    // INCHI✔️❌: Performance inherits the known SourceHeap overhead of INCHIGEN_DoCanonicalization.
    // END INCHI ACTIVE MACRO CONFIGURATION: STDINCHIGEN_DoCanonicalization

    INCHIGEN_DoCanonicalization(heap, h_gen, p_gen_data)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn INCHIGEN_DoCanonicalization(
    heap: &mut SourceHeap,
    h_gen: INCHIGEN_HANDLE,
    p_gen_data: Option<&mut INCHIGEN_DATA>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:697 INCHIGEN_DoCanonicalization
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap handle access and cloned control state add overhead.
    /*
int INCHI_DECL INCHIGEN_DoCanonicalization
( INCHIGEN_HANDLE _HGen, INCHIGEN_DATA *pGenData )
{
    int nRet = 0, nRet1 /*, maxINChI=0*/;
    INCHIGEN_CONTROL * HGen = (INCHIGEN_CONTROL *) _HGen;


    STRUCT_DATA *sd = &( HGen->StructData );
    INPUT_PARMS *ip = &( HGen->InpParms );
    INCHI_IOSTREAM *out_file = HGen->inchi_file, *log_file = HGen->inchi_file + 1;
    INCHI_IOSTREAM prbstr, *prb_file = &prbstr;
    INCHI_IOSTREAM inpstr, *inp_file = &inpstr;

    ORIG_ATOM_DATA *prep_inp_data = &( HGen->PrepInpData[0] );

    INCHI_CLOCK ic;
    CANON_GLOBALS CG;

    int k;

    memset( &ic, 0, sizeof( ic ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( &CG, 0, sizeof( CG ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    /* Set debug output */
#if (TRACE_MEMORY_LEAKS == 1)

    _CrtSetDbgFlag( _CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF );
    /* for execution outside the VC++ debugger uncomment one of the following two */
#ifdef MY_REPORT_FILE
    _CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_WARN, MY_REPORT_FILE );
    _CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_ERROR, MY_REPORT_FILE );
    _CrtSetReportMode( _CRT_ASSERT, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_ASSERT, MY_REPORT_FILE );
#else
    _CrtSetReportMode( _CRT_WARN | _CRT_ERROR, _CRTDBG_MODE_DEBUG );
#endif
#if ( !defined(__STDC__) || __STDC__ != 1 )
    /* turn on floating point exceptions */
    {
        /* Get the default control word. */
        int cw = _controlfp( 0, 0 );
        /* Set the exception masks OFF, turn exceptions on. */
        /*cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_INEXACT|EM_ZERODIVIDE|EM_DENORMAL);*/
        cw &= ~( EM_OVERFLOW | EM_UNDERFLOW | EM_ZERODIVIDE | EM_DENORMAL );
        /* Set the control word. */
        _controlfp( cw, MCW_EM );
    }
#endif /* ( !defined(__STDC__) || __STDC__ != 1 ) */

#endif /* (TRACE_MEMORY_LEAKS == 1) */



    if (HGen->norm_passed == 0)
    {
        AddErrorMessage( sd->pStrErrStruct, "Got non-normalized structure" );
        sd->nStructReadError = 99;
        sd->nErrorType = _IS_ERROR;
        nRet = _IS_ERROR;
        goto exit_function;
    }


    inchi_ios_init( inp_file, INCHI_IOS_TYPE_FILE, NULL );
    inchi_ios_init( prb_file, INCHI_IOS_TYPE_FILE, NULL );

    sd->bUserQuit = 0;
    if (sd->bUserQuit)
    {
        goto exit_function;
    }

   /* create INChI for each connected component of the structure and optionally display them */
   /* output INChI for the whole structure */



    /* create INChI for each connected component of the structure and optionally display them */
    /* create INChI for the whole disconnected or original structure */
    if (nRet != _IS_FATAL && nRet != _IS_ERROR)
    {
        nRet1 = CanonOneStructureINChI( &CG, &ic, HGen, INCHI_BAS, inp_file );
        nRet = inchi_max( nRet, nRet1 );
    }

    /*
    if ( nRet != _IS_FATAL && nRet != _IS_ERROR )
        maxINChI = 1;
    */

    if (nRet != _IS_FATAL && nRet != _IS_ERROR &&
        ( sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE ) &&
        ( ip->bTautFlags               & TG_FLAG_RECONNECT_COORD ))
    {
        /* create INChI for the whole reconnected structure */
        nRet1 = CanonOneStructureINChI( &CG, &ic, HGen, INCHI_REC, inp_file );
        nRet = inchi_max( nRet, nRet1 );
        /*
        if ( nRet != _IS_FATAL && nRet != _IS_ERROR )
                maxINChI = 2;
        */
    }

    if (nRet != _IS_FATAL && nRet != _IS_ERROR)
    {
        if (( sd->bChiralFlag & FLAG_INP_AT_CHIRAL ) &&
            ( ip->nMode & REQ_MODE_STEREO ) &&
              !( ip->nMode & ( REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO ) ) &&
              !bIsStructChiral( HGen->pINChI, sd->num_components ))
        {
            WarningMessage( sd->pStrErrStruct, "Not chiral" );
        }

        /*************************************/
        /*       Output err/warn messages    */
        /*************************************/
        if ( /*!sd->nErrorCode &&*/ !sd->bUserQuitComponent && !sd->bUserQuit)
        {
            /*  if successful then returns 0, otherwise returns _IS_FATAL */
            /*  extract the structure if requested */
            nRet1 = TreatCreateINChIWarning( sd, ip, prep_inp_data, HGen->num_inp,
                                 inp_file, log_file, out_file, prb_file );
            nRet = inchi_max( nRet, nRet1 );
        }
    }

    switch (nRet)
    {
        case _IS_SKIP: nRet = inchi_Ret_SKIP; break; /* not used in INChI dll */
        case _IS_EOF: nRet = inchi_Ret_EOF; break; /* no structural data has been provided */
        case _IS_OKAY: nRet = inchi_Ret_OKAY; HGen->canon_passed = 1; break;
                                                    /* Success; break; no errors or warnings */
        case _IS_WARNING: nRet = inchi_Ret_WARNING; HGen->canon_passed = 1; break;
                                                    /* Success; break; warning(s) issued */
        case _IS_ERROR: nRet = inchi_Ret_ERROR; break; /* Error: no INChI has been created */
        case _IS_FATAL: nRet = inchi_Ret_FATAL; break; /* Severe error: no INChI has been created (typically; break; memory allocation failed) */
        case _IS_UNKNOWN:
        default: nRet = inchi_Ret_UNKNOWN; break; /* Unknown program error */
    }
exit_function:

    strcpy( pGenData->pStrErrStruct, sd->pStrErrStruct );
    for (k = 0; k < INCHI_NUM; k++)
    {
        pGenData->num_components[k] = sd->num_components[k];
    }

    return nRet;
} /* INCHIGEN_DoCanonicalization */
    */
    // END INCHI C FUNCTION: INCHIGEN_DoCanonicalization
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: INCHIGEN_DoCanonicalization
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; TRACE_MEMORY_LEAKS == 0.
    // INCHI✔️❌: WarningMessage is the active AddErrorMessage macro alias.
    // END INCHI ACTIVE MACRO CONFIGURATION: INCHIGEN_DoCanonicalization

    let control_pointer = h_gen.cast::<INCHIGEN_CONTROL>();
    let mut control = heap
        .slice(control_pointer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let generation_data = p_gen_data.ok_or(SourceHeapError::NullPointer)?;
    let clock = heap.allocate_model_storage(vec![INCHI_CLOCK::default()])?;
    let mut canonical_globals = CANON_GLOBALS::default();

    let operation = (|| -> Result<i32, SourceHeapError> {
        let mut result = 0_i32;
        let mut input_stream = crate::source_types::INCHI_IOSTREAM::default();
        let mut problem_file = crate::source_types::INCHI_IOSTREAM::default();

        'canonicalization: {
            if control.norm_passed == 0 {
                let message = b"Got non-normalized structure\0".map(|byte| byte as i8);
                AddErrorMessage(Some(&mut control.StructData.pStrErrStruct), Some(&message))?;
                control.StructData.nStructReadError = 99;
                control.StructData.nErrorType = _IS_ERROR as i32;
                result = _IS_ERROR as i32;
                break 'canonicalization;
            }

            inchi_ios_init(
                Some(&mut input_stream),
                INCHI_IOS_TYPE_FILE as i32,
                SourceMutPointer::null(),
            )?;
            inchi_ios_init(
                Some(&mut problem_file),
                INCHI_IOS_TYPE_FILE as i32,
                SourceMutPointer::null(),
            )?;
            control.StructData.bUserQuit = 0;
            if control.StructData.bUserQuit != 0 {
                break 'canonicalization;
            }

            if result != _IS_FATAL as i32 && result != _IS_ERROR as i32 {
                let component_result = CanonOneStructureINChI(
                    heap,
                    &mut canonical_globals,
                    clock,
                    &mut control,
                    INCHI_BAS as i32,
                    Some(&mut input_stream),
                    0,
                )?;
                result = result.max(component_result);
            }

            if result != _IS_FATAL as i32
                && result != _IS_ERROR as i32
                && control.StructData.bTautFlagsDone[INCHI_BAS as usize]
                    & u64::from(TG_FLAG_DISCONNECT_COORD_DONE)
                    != 0
                && control.InpParms.bTautFlags & u64::from(TG_FLAG_RECONNECT_COORD) != 0
            {
                let component_result = CanonOneStructureINChI(
                    heap,
                    &mut canonical_globals,
                    clock,
                    &mut control,
                    INCHI_REC as i32,
                    Some(&mut input_stream),
                    0,
                )?;
                result = result.max(component_result);
            }

            if result != _IS_FATAL as i32 && result != _IS_ERROR as i32 {
                if control.StructData.bChiralFlag & FLAG_INP_AT_CHIRAL as i32 != 0
                    && control.InpParms.nMode & u64::from(REQ_MODE_STEREO) != 0
                    && control.InpParms.nMode
                        & u64::from(REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO)
                        == 0
                    && bIsStructChiral(
                        heap,
                        control.pINChI,
                        control.StructData.num_components,
                    )? == 0
                {
                    let message = b"Not chiral\0".map(|byte| byte as i8);
                    AddErrorMessage(
                        Some(&mut control.StructData.pStrErrStruct),
                        Some(&message),
                    )?;
                }

                if control.StructData.bUserQuitComponent == 0
                    && control.StructData.bUserQuit == 0
                {
                    let input_parameters = control.InpParms.clone();
                    let prepared = control.PrepInpData[0].clone();
                    let input_number = control.num_inp;
                    let (out_file, remainder) = control.inchi_file.split_at_mut(1);
                    let warning_result = TreatCreateINChIWarning(
                        heap,
                        &mut control.StructData,
                        &input_parameters,
                        &prepared,
                        input_number,
                        Some(&mut input_stream),
                        Some(&mut remainder[0]),
                        Some(&mut out_file[0]),
                        Some(&mut problem_file),
                    )?;
                    result = result.max(warning_result);
                }
            }
        }

        result = match result {
            value if value == _IS_SKIP as i32 => tagRetValGetINCHI_inchi_Ret_SKIP,
            value if value == _IS_EOF as i32 => tagRetValGetINCHI_inchi_Ret_EOF,
            value if value == _IS_OKAY as i32 => {
                control.canon_passed = 1;
                tagRetValGetINCHI_inchi_Ret_OKAY
            }
            value if value == _IS_WARNING as i32 => {
                control.canon_passed = 1;
                tagRetValGetINCHI_inchi_Ret_WARNING
            }
            value if value == _IS_ERROR as i32 => tagRetValGetINCHI_inchi_Ret_ERROR,
            value if value == _IS_FATAL as i32 => tagRetValGetINCHI_inchi_Ret_FATAL,
            _ => tagRetValGetINCHI_inchi_Ret_UNKNOWN,
        };

        let message_length = control
            .StructData
            .pStrErrStruct
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        generation_data.pStrErrStruct[..=message_length]
            .copy_from_slice(&control.StructData.pStrErrStruct[..=message_length]);
        for k in 0..INCHI_NUM as usize {
            generation_data.num_components[k] = control.StructData.num_components[k];
        }
        Ok(result)
    })();

    let free_result = heap.free(clock);
    *heap
        .slice_mut(control_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = control;
    let result = operation?;
    free_result?;
    Ok(result)
}

#[allow(non_snake_case)]
pub(crate) fn STDINCHIGEN_DoSerialization(
    heap: &mut SourceHeap,
    h_gen: INCHIGEN_HANDLE,
    p_gen_data: Option<&mut INCHIGEN_DATA>,
    p_results: Option<&mut inchi_Output>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:859 STDINCHIGEN_DoSerialization
    // INCHI✔️❌: int INCHI_DECL STDINCHIGEN_DoSerialization( INCHIGEN_HANDLE HGen,
    // INCHI✔️❌:                                                                      INCHIGEN_DATA * pGenData,
    // INCHI✔️❌:                                                                      inchi_Output * pResults )
    // INCHI✔️❌: {
    // INCHI✔️❌:     return INCHIGEN_DoSerialization( HGen, pGenData, pResults );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: STDINCHIGEN_DoSerialization
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: STDINCHIGEN_DoSerialization
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; INCHI_DECL has no behavioral body.
    // INCHI✔️❌: Performance inherits the known SourceHeap overhead of INCHIGEN_DoSerialization.
    // END INCHI ACTIVE MACRO CONFIGURATION: STDINCHIGEN_DoSerialization

    INCHIGEN_DoSerialization(heap, h_gen, p_gen_data, p_results)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn INCHIGEN_DoSerialization(
    heap: &mut SourceHeap,
    h_gen: INCHIGEN_HANDLE,
    p_gen_data: Option<&mut INCHIGEN_DATA>,
    p_results: Option<&mut inchi_Output>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:869 INCHIGEN_DoSerialization
    // INCHI✔️❌: complete source frame follows verbatim; SourceHeap stack-object and checked pointer modeling add overhead.
    /*
int INCHI_DECL INCHIGEN_DoSerialization(INCHIGEN_HANDLE _HGen,
    INCHIGEN_DATA* pGenData,
    inchi_Output* pResults)
{
    int nRet = 0, nRet1 = 0, i, k;



    INCHIGEN_CONTROL* HGen = (INCHIGEN_CONTROL*)_HGen;

    INPUT_PARMS* ip = &(HGen->InpParms);
    INCHI_IOSTREAM* out_file = HGen->inchi_file, * log_file = HGen->inchi_file + 1;
    INCHI_IOSTREAM inpstr, * inp_file = &inpstr;
    INCHI_IOSTREAM prbstr, * prb_file = &prbstr;

    STRUCT_DATA* sd = &(HGen->StructData);
    NORM_CANON_FLAGS* pncFlags = &(HGen->ncFlags);
    ORIG_ATOM_DATA* orig_inp_data = &(HGen->OrigInpData);
    ORIG_ATOM_DATA* prep_inp_data = &(HGen->PrepInpData[0]);
    ORIG_STRUCT* pOrigStruct = &(HGen->OrigStruct);
    int bSortPrintINChIFlags = 0;
    unsigned char save_opt_bits = 0;
    int retcode = 0;

    CANON_GLOBALS CG;
    memset(&CG, 0, sizeof(CG)); /* djb-rwth: memset_s C11/Annex K variant? */

    /* Post-1.02b - added initialization of pResults to 0; thanks to David Foss */
    memset(pResults, 0, sizeof(*pResults)); /* djb-rwth: memset_s C11/Annex K variant? */
    pResults->szLog = log_file->s.pStr;
    inchi_ios_init(inp_file, INCHI_IOS_TYPE_FILE, NULL);
    inchi_ios_init(prb_file, INCHI_IOS_TYPE_FILE, NULL);

    /* Set debug output */

#if (TRACE_MEMORY_LEAKS == 1)

    _CrtSetDbgFlag(_CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF);

    /* for execution outside the VC++ debugger uncomment one of the following two */
#ifdef MY_REPORT_FILE
    _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_WARN, MY_REPORT_FILE);
    _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_ERROR, MY_REPORT_FILE);
    _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_ASSERT, MY_REPORT_FILE);
#else
    _CrtSetReportMode(_CRT_WARN | _CRT_ERROR, _CRTDBG_MODE_DEBUG);
#endif
#if ( !defined(__STDC__) || __STDC__ != 1 )
    /* turn on floating point exceptions */
    {
        /* Get the default control word. */
        int cw = _controlfp(0, 0);
        /* Set the exception masks OFF, turn exceptions on. */
        /*cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_INEXACT|EM_ZERODIVIDE|EM_DENORMAL);*/
        cw &= ~(EM_OVERFLOW | EM_UNDERFLOW | EM_ZERODIVIDE | EM_DENORMAL);
        /* Set the control word. */
        _controlfp(cw, MCW_EM);
    }
#endif /* ( !defined(__STDC__) || __STDC__ != 1 ) */

#endif /* (TRACE_MEMORY_LEAKS == 1) */


    /*****************************/


    if (HGen->canon_passed == 0)
    {
        AddErrorMessage(sd->pStrErrStruct, "Got non-canonicalized structure");
        sd->nStructReadError = 99;
        sd->nErrorType = _IS_ERROR;
        retcode = _IS_ERROR;
        goto frees;
    }

    /************************************************/
    /*  sort and print INChI for the whole structure */
    /************************************************/

    /* Prepare SaveOpt bits */
    if (ip->bINChIOutputOptions & INCHI_OUT_SAVEOPT)
    {
        if (0 != (ip->bTautFlags & TG_FLAG_RECONNECT_COORD))
        {
            save_opt_bits |= SAVE_OPT_RECMET;
        }
        if (0 != (ip->nMode & REQ_MODE_BASIC))
        {
            save_opt_bits |= SAVE_OPT_FIXEDH;
        }
        if (0 != (ip->nMode & REQ_MODE_DIFF_UU_STEREO))
        {
            save_opt_bits |= SAVE_OPT_SLUUD;
        }
        if (0 == (ip->nMode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)))
        {
            save_opt_bits |= SAVE_OPT_SUU;
        }
        if (0 != (ip->bTautFlags & TG_FLAG_KETO_ENOL_TAUT))
        {
            save_opt_bits |= SAVE_OPT_KET;
        }
        if (0 != (ip->bTautFlags & TG_FLAG_1_5_TAUT))
        {
            save_opt_bits |= SAVE_OPT_15T;
        }
    }

    nRet = SortAndPrintINChI(&CG, out_file, &(HGen->strbuf_container),
        log_file, ip, orig_inp_data, prep_inp_data,
        HGen->composite_norm_data, pOrigStruct,
        sd->num_components, sd->num_non_taut, sd->num_taut,
        sd->bTautFlags, sd->bTautFlagsDone,
        pncFlags, HGen->num_inp,
        HGen->pINChI, HGen->pINChI_Aux,
        &bSortPrintINChIFlags, save_opt_bits);

    if (nRet != _IS_FATAL && nRet != _IS_ERROR)
    {
        /* Special mode: extract all good MOLfiles into the problem file
        * Do not extract any MOLfile that could not be processed (option /PGO)
        */
        if (prb_file->f && 0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd && ip->bSaveAllGoodStructsAsProblem)
        {
            MolfileSaveCopy(inp_file, sd->fPtrStart, sd->fPtrEnd, prb_file->f, 0);
        }
#if( /*bRELEASE_VERSION != 1 &&*/ EXTR_FLAGS == EXTR_TRANSPOSITION_EXAMPLES && EXTR_MASK == EXTR_FLAGS )
        else
            if (prb_file->f && (bSortPrintINChIFlags &
                (FLAG_SORT_PRINT_TRANSPOS_BAS | FLAG_SORT_PRINT_TRANSPOS_REC))
                )
            {
                MolfileSaveCopy(inp_file, sd->fPtrStart, sd->fPtrEnd, prb_file->f, 0);
            }
#endif
    }

    for (i = 0; i < INCHI_NUM; i++)
    {
        for (k = 0; k < TAUT_NUM + 1; k++)
        {
            FreeCompAtomData(&(HGen->composite_norm_data[i][k]));
        }
    }

    /*****************************/

    /* Prepare output message(s). */

    /* Error/warning. */
    if (sd->pStrErrStruct[0])
    {
        if (pGenData && (pResults->szMessage = (char*)inchi_malloc(strlen(sd->pStrErrStruct) + 1)))
        {
            strcpy(pResults->szMessage, sd->pStrErrStruct);
        }
    }

    /* InChI, AuxInfo  (go to  pResults->szInChI, pResults->szAuxInfo) */
    if (out_file->s.pStr && out_file->s.nUsedLength > 0 && pGenData)
    {
        char* p;
        pResults->szInChI = out_file->s.pStr;
        pResults->szAuxInfo = NULL;
        if (!(INCHI_OUT_SDFILE_ONLY & ip->bINChIOutputOptions)) /* do not remove last LF from SDF output - 2008-12-23 DT */
        {
            for (p = strchr(pResults->szInChI, '\n'); p; p = strchr(p + 1, '\n'))
            {
                if (!memcmp(p, "\nAuxInfo", 8))
                {
                    *p = '\0';            /* remove LF after INChI */
                    pResults->szAuxInfo = p + 1; /* save pointer to AuxInfo */
                }
                else
                {
                    if (pResults->szAuxInfo || !p[1])
                    {
                        /* remove LF after aux info or from the last char */
                        *p = '\0';
                        break;
                    }
                }
            }
        }
        out_file->s.pStr = NULL;
    }

    /* Log message. */
    if (log_file->s.pStr && log_file->s.nUsedLength > 0)
    {
        while (log_file->s.nUsedLength && '\n' == log_file->s.pStr[log_file->s.nUsedLength - 1])
            log_file->s.pStr[--log_file->s.nUsedLength] = '\0'; /* remove last LF */
        if (pGenData)
        {
            pResults->szLog = log_file->s.pStr;
            log_file->s.pStr = NULL;
        }
    }

    if (out_file->s.pStr)
    {
        inchi_free(out_file->s.pStr);
        out_file->s.pStr = NULL;
    }
    if (log_file->s.pStr)
    {
        inchi_free(log_file->s.pStr);
        log_file->s.pStr = NULL;
    }

    HGen->ulTotalProcessingTime += sd->ulStructTime;
    nRet = inchi_max(nRet, nRet1);

    switch (nRet)
    {
    case _IS_FATAL:
    case _IS_ERROR:     HGen->num_err++;
    }

frees:
    /* Free all. */
    for (i = 0; i < MAX_NUM_PATHS; i++)
    {
        if (ip->path[i])
        {
            inchi_free((char*)ip->path[i]); /*  cast deliberately discards 'const' qualifier */
            ip->path[i] = NULL;
        }
    }
    SetBitFree(&CG);


    if (pGenData) /* djb-rwth: fixing a NULL pointer dereference */
    {
        strcpy(pGenData->pStrErrStruct, sd->pStrErrStruct);
        for (k = 0; k < INCHI_NUM; k++)
        {
            pGenData->num_components[k] = sd->num_components[k];
        }
    }

    return retcode;
}
    */
    // END INCHI C FUNCTION: INCHIGEN_DoSerialization
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: INCHIGEN_DoSerialization
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; TRACE_MEMORY_LEAKS == 0.
    // INCHI✔️❌: EXTR_FLAGS == 0, EXTR_MASK == 0, and EXTR_TRANSPOSITION_EXAMPLES != 0, so the transposition copy branch is inactive.
    // INCHI✔️❌: prb_file is a newly initialized file stream with a null FILE pointer, so the active PGO copy condition is unreachable in this function.
    // END INCHI ACTIVE MACRO CONFIGURATION: INCHIGEN_DoSerialization

    let control_pointer = h_gen.cast::<INCHIGEN_CONTROL>();
    let mut control = heap
        .slice(control_pointer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let results = p_results.ok_or(SourceHeapError::NullPointer)?;
    *results = inchi_Output::default();
    results.szLog = control.inchi_file[1].s.pStr;

    let mut input_stream = crate::source_types::INCHI_IOSTREAM::default();
    let mut problem_stream = crate::source_types::INCHI_IOSTREAM::default();
    inchi_ios_init(
        Some(&mut input_stream),
        INCHI_IOS_TYPE_FILE as i32,
        SourceMutPointer::null(),
    )?;
    inchi_ios_init(
        Some(&mut problem_stream),
        INCHI_IOS_TYPE_FILE as i32,
        SourceMutPointer::null(),
    )?;

    let canonical_pointer = heap.allocate_model_storage(vec![CANON_GLOBALS::default()])?;
    let original_pointer = heap.allocate_model_storage(vec![control.OrigInpData.clone()])?;
    let prepared_pointer = heap.allocate_model_storage(vec![control.PrepInpData[0].clone()])?;
    let original_structure_pointer = heap.allocate_model_storage(vec![control.OrigStruct.clone()])?;
    let sort_flags_pointer = heap.allocate_model_storage(vec![0_i32])?;

    let operation = (|| -> Result<i32, SourceHeapError> {
        let mut retcode = 0_i32;
        let mut n_ret = 0_i32;

        'serialization: {
            if control.canon_passed == 0 {
                let message = b"Got non-canonicalized structure\0".map(|byte| byte as i8);
                AddErrorMessage(Some(&mut control.StructData.pStrErrStruct), Some(&message))?;
                control.StructData.nStructReadError = 99;
                control.StructData.nErrorType = _IS_ERROR as i32;
                retcode = _IS_ERROR as i32;
                break 'serialization;
            }

            let mut save_option_bits = 0_u8;
            if control.InpParms.bINChIOutputOptions & INCHI_OUT_SAVEOPT as i32 != 0 {
                if control.InpParms.bTautFlags & u64::from(TG_FLAG_RECONNECT_COORD) != 0 {
                    save_option_bits |= SAVE_OPT_RECMET as u8;
                }
                if control.InpParms.nMode & u64::from(REQ_MODE_BASIC) != 0 {
                    save_option_bits |= SAVE_OPT_FIXEDH as u8;
                }
                if control.InpParms.nMode & u64::from(REQ_MODE_DIFF_UU_STEREO) != 0 {
                    save_option_bits |= SAVE_OPT_SLUUD as u8;
                }
                if control.InpParms.nMode
                    & u64::from(REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)
                    == 0
                {
                    save_option_bits |= SAVE_OPT_SUU as u8;
                }
                if control.InpParms.bTautFlags & u64::from(TG_FLAG_KETO_ENOL_TAUT) != 0 {
                    save_option_bits |= SAVE_OPT_KET as u8;
                }
                if control.InpParms.bTautFlags & u64::from(TG_FLAG_1_5_TAUT) != 0 {
                    save_option_bits |= SAVE_OPT_15T as u8;
                }
            }

            let input_parameters = control.InpParms.clone();
            let normalization_flags = control.ncFlags.clone();
            let component_counts = control.StructData.num_components;
            let non_tautomeric_counts = control.StructData.num_non_taut;
            let tautomeric_counts = control.StructData.num_taut;
            let mut taut_flags = control.StructData.bTautFlags;
            let mut taut_flags_done = control.StructData.bTautFlagsDone;
            let input_number = control.num_inp;
            let inchi_components = control.pINChI;
            let aux_components = control.pINChI_Aux;
            let (out_file, remaining_files) = control.inchi_file.split_at_mut(1);
            let log_file = &mut remaining_files[0];
            n_ret = SortAndPrintINChI(
                heap,
                canonical_pointer,
                &mut out_file[0],
                Some(&mut control.strbuf_container),
                log_file,
                &input_parameters,
                original_pointer.as_const(),
                prepared_pointer.as_const(),
                Some(&mut control.composite_norm_data),
                original_structure_pointer.as_const(),
                &component_counts,
                &non_tautomeric_counts,
                &tautomeric_counts,
                &mut taut_flags,
                &mut taut_flags_done,
                &normalization_flags,
                input_number,
                inchi_components,
                aux_components,
                sort_flags_pointer,
                save_option_bits,
                SourceMutPointer::null(),
            )?;
            control.StructData.bTautFlags = taut_flags;
            control.StructData.bTautFlagsDone = taut_flags_done;

            for domain in &mut control.composite_norm_data {
                for representation in domain.iter_mut().take(TAUT_NUM as usize + 1) {
                    FreeCompAtomData(heap, representation)?;
                }
            }

            if control.StructData.pStrErrStruct[0] != 0 {
                if p_gen_data.is_some() {
                    let message_length = control
                        .StructData
                        .pStrErrStruct
                        .iter()
                        .position(|byte| *byte == 0)
                        .ok_or(SourceHeapError::MissingNulTerminator)?;
                    let message_pointer = match inchi_malloc(heap, (message_length + 1) as u64) {
                        Ok(pointer) => pointer,
                        Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                        Err(error) => return Err(error),
                    };
                    results.szMessage = message_pointer;
                    if !message_pointer.is_null() {
                        heap.slice_mut(message_pointer)?[..=message_length]
                            .copy_from_slice(&control.StructData.pStrErrStruct[..=message_length]);
                    }
                }
            }

            let output_pointer = control.inchi_file[0].s.pStr;
            if !output_pointer.is_null()
                && control.inchi_file[0].s.nUsedLength > 0
                && p_gen_data.is_some()
            {
                results.szInChI = output_pointer;
                results.szAuxInfo = SourceMutPointer::null();
                if control.InpParms.bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY as i32 == 0 {
                    let allocation = heap.slice(output_pointer.as_const())?;
                    let nul = allocation
                        .iter()
                        .position(|byte| *byte == 0)
                        .ok_or(SourceHeapError::MissingNulTerminator)?;
                    let mut search = 0_usize;
                    loop {
                        let next = heap.slice(output_pointer.as_const())?[search..nul]
                            .iter()
                            .position(|byte| *byte == b'\n' as i8)
                            .map(|offset| search + offset);
                        let Some(position) = next else {
                            break;
                        };
                        let is_aux = nul - position >= 8
                            && heap.slice(output_pointer.as_const())?[position..position + 8]
                                == [
                                    b'\n' as i8,
                                    b'A' as i8,
                                    b'u' as i8,
                                    b'x' as i8,
                                    b'I' as i8,
                                    b'n' as i8,
                                    b'f' as i8,
                                    b'o' as i8,
                                ];
                        if is_aux {
                            heap.slice_mut(output_pointer)?[position] = 0;
                            results.szAuxInfo = output_pointer.offset((position + 1) as i64)?;
                            search = position + 1;
                        } else {
                            let next_is_nul = heap
                                .slice(output_pointer.as_const())?
                                .get(position + 1)
                                == Some(&0);
                            if !results.szAuxInfo.is_null() || next_is_nul {
                                heap.slice_mut(output_pointer)?[position] = 0;
                                break;
                            }
                            search = position + 1;
                        }
                    }
                }
                control.inchi_file[0].s.pStr = SourceMutPointer::null();
            }

            if !control.inchi_file[1].s.pStr.is_null()
                && control.inchi_file[1].s.nUsedLength > 0
            {
                while control.inchi_file[1].s.nUsedLength > 0 {
                    let last = usize::try_from(control.inchi_file[1].s.nUsedLength - 1)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    if heap.slice(control.inchi_file[1].s.pStr.as_const())?[last]
                        != b'\n' as i8
                    {
                        break;
                    }
                    control.inchi_file[1].s.nUsedLength -= 1;
                    heap.slice_mut(control.inchi_file[1].s.pStr)?[last] = 0;
                }
                if p_gen_data.is_some() {
                    results.szLog = control.inchi_file[1].s.pStr;
                    control.inchi_file[1].s.pStr = SourceMutPointer::null();
                }
            }

            if !control.inchi_file[0].s.pStr.is_null() {
                inchi_free(heap, control.inchi_file[0].s.pStr)?;
                control.inchi_file[0].s.pStr = SourceMutPointer::null();
            }
            if !control.inchi_file[1].s.pStr.is_null() {
                inchi_free(heap, control.inchi_file[1].s.pStr)?;
                control.inchi_file[1].s.pStr = SourceMutPointer::null();
            }

            control.ulTotalProcessingTime = control
                .ulTotalProcessingTime
                .wrapping_add(control.StructData.ulStructTime);
            n_ret = n_ret.max(0);
            if n_ret == _IS_FATAL as i32 || n_ret == _IS_ERROR as i32 {
                control.num_err = control.num_err.wrapping_add(1);
            }
        }

        for path in control
            .InpParms
            .path
            .iter_mut()
            .take(MAX_NUM_PATHS as usize)
        {
            if !path.is_null() {
                inchi_free(heap, path.as_mut())?;
                *path = SourceConstPointer::null();
            }
        }
        let mut canonical_globals = heap.slice(canonical_pointer.as_const())?[0].clone();
        SetBitFree(heap, &mut canonical_globals)?;
        heap.slice_mut(canonical_pointer)?[0] = canonical_globals;

        if let Some(generation_data) = p_gen_data {
            let message_length = control
                .StructData
                .pStrErrStruct
                .iter()
                .position(|byte| *byte == 0)
                .ok_or(SourceHeapError::MissingNulTerminator)?;
            generation_data.pStrErrStruct[..=message_length]
                .copy_from_slice(&control.StructData.pStrErrStruct[..=message_length]);
            for k in 0..INCHI_NUM as usize {
                generation_data.num_components[k] = control.StructData.num_components[k];
            }
        }

        Ok(retcode)
    })();

    let cleanup = [
        heap.free(sort_flags_pointer),
        heap.free(original_structure_pointer),
        heap.free(prepared_pointer),
        heap.free(original_pointer),
        heap.free(canonical_pointer),
    ];
    *heap
        .slice_mut(control_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = control;
    let result = operation?;
    for cleanup_result in cleanup {
        cleanup_result?;
    }
    Ok(result)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn STDINCHIGEN_Create(
    heap: &mut SourceHeap,
) -> Result<INCHIGEN_HANDLE, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:122 STDINCHIGEN_Create
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
INCHIGEN_HANDLE INCHI_DECL STDINCHIGEN_Create( void )
{
    return INCHIGEN_Create( );
}
    */
    // END INCHI C FUNCTION: STDINCHIGEN_Create
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: STDINCHIGEN_Create
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; no function-local conditional compilation.
    // INCHI✔️❌: Performance inherits the SourceHeap allocation overhead of INCHIGEN_Create.
    // END INCHI ACTIVE MACRO CONFIGURATION: STDINCHIGEN_Create

    INCHIGEN_Create(heap)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn INCHIGEN_Create(
    heap: &mut SourceHeap,
) -> Result<INCHIGEN_HANDLE, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:130 INCHIGEN_Create
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
INCHIGEN_HANDLE INCHI_DECL INCHIGEN_Create( void )
{
    INCHIGEN_CONTROL * HGen = NULL;


    HGen = (INCHIGEN_CONTROL *) inchi_malloc( sizeof( INCHIGEN_CONTROL ) );

    if (!HGen)
    {
        return (INCHIGEN_HANDLE) NULL;
    }

    memset( HGen, 0, sizeof( INCHIGEN_CONTROL ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    /* Set/init aliases */
    memset( &( HGen->InpParms ), 0, sizeof( INPUT_PARMS ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( &( HGen->StructData ), 0, sizeof( STRUCT_DATA ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    HGen->ulTotalProcessingTime = 0;
    HGen->num_err = 0;
    HGen->num_inp = 0;
    HGen->szTitle[0] = '\0';


    /* Initialize output streams as string buffers */
    inchi_ios_init( &( HGen->inchi_file[0] ), INCHI_IOS_TYPE_STRING, NULL );
    inchi_ios_init( &( HGen->inchi_file[1] ), INCHI_IOS_TYPE_STRING, NULL );
    inchi_ios_init( &( HGen->inchi_file[2] ), INCHI_IOS_TYPE_STRING, NULL );


    memset( &( HGen->OrigInpData ), 0, sizeof( HGen->OrigInpData ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( &( HGen->PrepInpData[0] ), 0, 2 * sizeof( HGen->PrepInpData[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    memset( HGen->pINChI, 0, sizeof( HGen->pINChI ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( HGen->pINChI_Aux, 0, sizeof( HGen->pINChI_Aux ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    /* Supply expandable string buffer */
    if (0 >= inchi_strbuf_init( &( HGen->strbuf_container ), INCHI_STRBUF_INITIAL_SIZE, INCHI_STRBUF_SIZE_INCREMENT ))
    {
        inchi_free( HGen );
        return (INCHIGEN_HANDLE) NULL;
    }

    return (INCHIGEN_HANDLE) HGen;
}
    */
    // END INCHI C FUNCTION: INCHIGEN_Create
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: INCHIGEN_Create
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; inchi_malloc is the active libc malloc macro.
    // INCHI✔️❌: Typed SourceHeap allocation and map bookkeeping are materially slower than native malloc and direct field stores.
    // END INCHI ACTIVE MACRO CONFIGURATION: INCHIGEN_Create

    let control_pointer = match heap.allocate(vec![INCHIGEN_CONTROL::default()]) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return Ok(SourceMutPointer::null()),
        Err(error) => return Err(error),
    };
    let mut control = INCHIGEN_CONTROL::default();
    for stream in &mut control.inchi_file {
        inchi_ios_init(
            Some(stream),
            INCHI_IOS_TYPE_STRING as i32,
            SourceMutPointer::null(),
        )?;
    }
    if inchi_strbuf_init(
        heap,
        &mut control.strbuf_container,
        INCHI_STRBUF_INITIAL_SIZE as i32,
        INCHI_STRBUF_SIZE_INCREMENT as i32,
    )? <= 0
    {
        inchi_free(heap, control_pointer)?;
        return Ok(SourceMutPointer::null());
    }
    *heap
        .slice_mut(control_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = control;
    Ok(control_pointer.cast())
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn STDINCHIGEN_Destroy(
    heap: &mut SourceHeap,
    h_gen: INCHIGEN_HANDLE,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1308 STDINCHIGEN_Destroy
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
void INCHI_DECL STDINCHIGEN_Destroy( INCHIGEN_HANDLE HGen )
{
    INCHIGEN_Destroy( HGen );
}
    */
    // END INCHI C FUNCTION: STDINCHIGEN_Destroy
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: STDINCHIGEN_Destroy
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; no function-local conditional compilation.
    // INCHI✔️❌: Performance inherits the SourceHeap and cloning overhead of INCHIGEN_Destroy.
    // END INCHI ACTIVE MACRO CONFIGURATION: STDINCHIGEN_Destroy

    INCHIGEN_Destroy(heap, h_gen)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn INCHIGEN_Destroy(
    heap: &mut SourceHeap,
    h_gen: INCHIGEN_HANDLE,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1315 INCHIGEN_Destroy
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
void INCHI_DECL INCHIGEN_Destroy( INCHIGEN_HANDLE _HGen )
{
    INCHIGEN_CONTROL * HGen = (INCHIGEN_CONTROL *) _HGen;

    if (NULL != HGen)
    {
        inchi_strbuf_close( &( HGen->strbuf_container ) );

        inchi_ios_close( &( HGen->inchi_file[0] ) );
        inchi_ios_close( &( HGen->inchi_file[1] ) );
        inchi_ios_close( &( HGen->inchi_file[2] ) );

        inchi_free( HGen );
    }
}
    */
    // END INCHI C FUNCTION: INCHIGEN_Destroy
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: INCHIGEN_Destroy
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; no function-local conditional compilation.
    // INCHI✔️❌: SourceHeap map access, control cloning, and checked frees are materially slower than native pointer operations.
    // END INCHI ACTIVE MACRO CONFIGURATION: INCHIGEN_Destroy

    if h_gen.is_null() {
        return Ok(());
    }
    let control_pointer = h_gen.cast::<INCHIGEN_CONTROL>();
    let mut control = heap
        .slice(control_pointer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    inchi_strbuf_close(heap, Some(&mut control.strbuf_container))?;
    for stream in &mut control.inchi_file {
        inchi_ios_close(heap, Some(stream))?;
    }
    inchi_free(heap, control_pointer)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn STDINCHIGEN_Reset(
    heap: &mut SourceHeap,
    h_gen: INCHIGEN_HANDLE,
    p_gen_data: Option<&mut INCHIGEN_DATA>,
    p_results: Option<&mut inchi_Output>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1125 STDINCHIGEN_Reset
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
void INCHI_DECL STDINCHIGEN_Reset( INCHIGEN_HANDLE HGen,
                               INCHIGEN_DATA * pGenData,
                               inchi_Output * pResults )
{
    INCHIGEN_Reset( HGen, pGenData, pResults );
}
    */
    // END INCHI C FUNCTION: STDINCHIGEN_Reset
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: STDINCHIGEN_Reset
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; no function-local conditional compilation.
    // INCHI✔️❌: Performance inherits the SourceHeap and cloning overhead of INCHIGEN_Reset.
    // END INCHI ACTIVE MACRO CONFIGURATION: STDINCHIGEN_Reset

    INCHIGEN_Reset(heap, h_gen, p_gen_data, p_results)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn INCHIGEN_Reset(
    heap: &mut SourceHeap,
    h_gen: INCHIGEN_HANDLE,
    p_gen_data: Option<&mut INCHIGEN_DATA>,
    p_results: Option<&mut inchi_Output>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1134 INCHIGEN_Reset
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
void INCHI_DECL INCHIGEN_Reset( INCHIGEN_HANDLE _HGen,
                                 INCHIGEN_DATA * pGenData,
                                 inchi_Output * pResults )
{
    int i, k, nc;
    INCHIGEN_CONTROL * HGen = (INCHIGEN_CONTROL *) _HGen;

    if (pResults->szInChI)
    {
        inchi_free( pResults->szInChI );
    }
    if (pResults->szLog)
    {
        inchi_free( pResults->szLog );
    }
    if (pResults->szMessage)
    {
        inchi_free( pResults->szMessage );
    }

    /* Free all data associated with components of disconn/conn structures */
    if (NULL != HGen)
    {

        /* Re-initialize output streams/string buffers */
        inchi_ios_close( &( HGen->inchi_file[0] ) );
        inchi_ios_close( &( HGen->inchi_file[1] ) );
        inchi_ios_close( &( HGen->inchi_file[2] ) );
        inchi_ios_init( &( HGen->inchi_file[0] ), INCHI_IOS_TYPE_STRING, NULL );
        inchi_ios_init( &( HGen->inchi_file[1] ), INCHI_IOS_TYPE_STRING, NULL );
        inchi_ios_init( &( HGen->inchi_file[2] ), INCHI_IOS_TYPE_STRING, NULL );

        inchi_strbuf_reset( &( HGen->strbuf_container ) );

        for (i = 0; i < MAX_NUM_PATHS; i++)
        {
            if (HGen->InpParms.path[i])
            {
                inchi_free( (char*) HGen->InpParms.path[i] ); /*  cast deliberately discards 'const' qualifier */
                HGen->InpParms.path[i] = NULL;
            }
        }
        memset( &( HGen->InpParms ), 0, sizeof( INPUT_PARMS ) ); /* djb-rwth: memset_s C11/Annex K variant? */

        FreeOrigAtData( &( HGen->OrigInpData ) );
        memset( &( HGen->OrigInpData ), 0, sizeof( HGen->OrigInpData ) ); /* djb-rwth: memset_s C11/Annex K variant? */


        FreeOrigAtData( &( HGen->PrepInpData[0] ) );
        FreeOrigAtData( &( HGen->PrepInpData[1] ) );
        memset( &( HGen->PrepInpData[0] ), 0, 2 * sizeof( HGen->PrepInpData[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */

        OrigStruct_Free( &( HGen->OrigStruct ) );
        memset( &( HGen->OrigStruct ), 0, sizeof( HGen->OrigStruct ) ); /* djb-rwth: memset_s C11/Annex K variant? */

        for (i = 0; i < INCHI_NUM; i++)
        {
            for (k = 0; k < TAUT_NUM + 1; k++)
            {
                FreeCompAtomData( &( HGen->composite_norm_data[i][k] ) );
            }
        }

        for (k = 0; k < INCHI_NUM; k++)
        {
            nc = HGen->StructData.num_components[k];

            if (HGen->InpCurAtData[k])
            {
                for (i = 0; i < nc; i++)
                {
                    FreeInpAtomData( &( HGen->InpCurAtData[k][i] ) );
                }
                inchi_free( HGen->InpCurAtData[k] );
                HGen->InpCurAtData[k] = NULL;
            }

            if (HGen->cti[k])
            {
                if (( HGen->cti[k] )->at[TAUT_YES])
                {
                    inchi_free( ( HGen->cti[k] )->at[TAUT_YES] );
                    ( HGen->cti[k] )->at[TAUT_YES] = NULL;
                }

                if (( HGen->cti[k] )->at[TAUT_NON])
                {
                    inchi_free( ( HGen->cti[k] )->at[TAUT_NON] );
                    ( HGen->cti[k] )->at[TAUT_NON] = NULL;
                }

                if (&( ( HGen->cti[k] )->vt_group_info )) /* djb-rwth: &((HGen->cti[k])->vt_group_info) will always evaluate to true?  */
                {
                    free_t_group_info( &( ( HGen->cti[k] )->vt_group_info ) );
                }

                if (&( ( HGen->cti[k] )->vt_group_info_orig )) /* djb-rwth: &((HGen->cti[k])->vt_group_info_orig) will always evaluate to true?  */
                {
                    free_t_group_info( &( ( HGen->cti[k] )->vt_group_info_orig ) );
                }

                inchi_free( HGen->cti[k] );
                HGen->cti[k] = NULL;
            }
        }

        for (k = 0; k < INCHI_NUM; k++)
        {
            nc = HGen->StructData.num_components[k];
            if (HGen->InpNormAtData[k])
            {
                for (i = 0; i < nc; i++)
                {
                    FreeInpAtomData( &( HGen->InpNormAtData[k][i] ) );
                }
                inchi_free( HGen->InpNormAtData[k] );
                HGen->InpNormAtData[k] = NULL;
            }

            if (HGen->InpNormTautData[k])
            {
                for (i = 0; i < nc; i++)
                {
                    FreeInpAtomData( &( HGen->InpNormTautData[k][i] ) );
                }
                inchi_free( HGen->InpNormTautData[k] );
                HGen->InpNormTautData[k] = NULL;
            }

            if (pGenData->NormAtomsTaut[k])
            {
                /*
                for ( i = 0;  i < nc; i ++ )
                    FreeInpAtomData( &(pGenData->NormAtomsTaut[k][i]) );
                */
                inchi_free( pGenData->NormAtomsTaut[k] );
                pGenData->NormAtomsTaut[k] = NULL;
            }
            if (pGenData->NormAtomsNontaut[k])
            {
                /* for ( i = 0;  i < nc; i ++ )
                    FreeInpAtomData( &(pGenData->NormAtomsNontaut[k][i]) );
                */
                inchi_free( pGenData->NormAtomsNontaut[k] );
                pGenData->NormAtomsNontaut[k] = NULL;
            }
        }

        /*  free INChI memory */
        FreeAllINChIArrays( HGen->pINChI, HGen->pINChI_Aux, HGen->StructData.num_components );
        memset( HGen->pINChI, 0, sizeof( HGen->pINChI ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        memset( HGen->pINChI_Aux, 0, sizeof( HGen->pINChI_Aux ) ); /* djb-rwth: memset_s C11/Annex K variant? */

        HGen->szTitle[0] = '\0';
    }

    if (HGen) /* djb-rwth: fixing a NULL pointer dereference */
        memset( &( HGen->StructData ), 0, sizeof( STRUCT_DATA ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( pResults, 0, sizeof( *pResults ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( pGenData, 0, sizeof( *pGenData ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    return;
}
    */
    // END INCHI C FUNCTION: INCHIGEN_Reset
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: INCHIGEN_Reset
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; no function-local conditional compilation.
    // INCHI✔️❌: SourceHeap map access, checked frees, and control-record cloning are materially slower than native pointer operations.
    // END INCHI ACTIVE MACRO CONFIGURATION: INCHIGEN_Reset

    let p_gen_data = p_gen_data.ok_or(SourceHeapError::NullPointer)?;
    let p_results = p_results.ok_or(SourceHeapError::NullPointer)?;

    if !p_results.szInChI.is_null() {
        inchi_free(heap, p_results.szInChI)?;
    }
    if !p_results.szLog.is_null() {
        inchi_free(heap, p_results.szLog)?;
    }
    if !p_results.szMessage.is_null() {
        inchi_free(heap, p_results.szMessage)?;
    }

    if h_gen.is_null() {
        *p_results = inchi_Output::default();
        *p_gen_data = INCHIGEN_DATA::default();
        return Ok(());
    }

    let control_pointer = h_gen.cast::<INCHIGEN_CONTROL>();
    let mut control = heap
        .slice(control_pointer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;

    let operation = (|| -> Result<(), SourceHeapError> {
        for stream in &mut control.inchi_file {
            inchi_ios_close(heap, Some(stream))?;
        }
        for stream in &mut control.inchi_file {
            inchi_ios_init(
                Some(stream),
                INCHI_IOS_TYPE_STRING as i32,
                SourceMutPointer::null(),
            )?;
        }

        inchi_strbuf_reset(heap, Some(&mut control.strbuf_container))?;

        for path in &mut control.InpParms.path {
            if !path.is_null() {
                inchi_free(heap, path.as_mut())?;
                *path = SourceConstPointer::null();
            }
        }
        control.InpParms = INPUT_PARMS::default();

        FreeOrigAtData(heap, Some(&mut control.OrigInpData))?;
        control.OrigInpData = ORIG_ATOM_DATA::default();
        for prepared in &mut control.PrepInpData {
            FreeOrigAtData(heap, Some(prepared))?;
        }
        control.PrepInpData = std::array::from_fn(|_| ORIG_ATOM_DATA::default());

        OrigStruct_Free(heap, Some(&mut control.OrigStruct))?;
        control.OrigStruct = ORIG_STRUCT::default();

        for representation in &mut control.composite_norm_data {
            for tautomer in representation {
                FreeCompAtomData(heap, tautomer)?;
            }
        }

        for representation in 0..INCHI_NUM as usize {
            let component_count = control.StructData.num_components[representation];
            let current = control.InpCurAtData[representation];
            if !current.is_null() {
                free_inp_atom_data_array(heap, current, component_count)?;
                control.InpCurAtData[representation] = SourceMutPointer::null();
            }

            let cti_pointer = control.cti[representation];
            if !cti_pointer.is_null() {
                let mut cti = heap
                    .slice(cti_pointer.as_const())?
                    .first()
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if !cti.at[TAUT_YES as usize].is_null() {
                    inchi_free(heap, cti.at[TAUT_YES as usize])?;
                    cti.at[TAUT_YES as usize] = SourceMutPointer::null();
                }
                if !cti.at[TAUT_NON as usize].is_null() {
                    inchi_free(heap, cti.at[TAUT_NON as usize])?;
                    cti.at[TAUT_NON as usize] = SourceMutPointer::null();
                }
                free_t_group_info(heap, Some(&mut cti.vt_group_info))?;
                free_t_group_info(heap, Some(&mut cti.vt_group_info_orig))?;
                *heap
                    .slice_mut(cti_pointer)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = cti;
                inchi_free(heap, cti_pointer)?;
                control.cti[representation] = SourceMutPointer::null();
            }
        }

        for representation in 0..INCHI_NUM as usize {
            let component_count = control.StructData.num_components[representation];
            let normalized = control.InpNormAtData[representation];
            if !normalized.is_null() {
                free_inp_atom_data_array(heap, normalized, component_count)?;
                control.InpNormAtData[representation] = SourceMutPointer::null();
            }

            let tautomeric = control.InpNormTautData[representation];
            if !tautomeric.is_null() {
                free_inp_atom_data_array(heap, tautomeric, component_count)?;
                control.InpNormTautData[representation] = SourceMutPointer::null();
            }

            if !p_gen_data.NormAtomsTaut[representation].is_null() {
                inchi_free(heap, p_gen_data.NormAtomsTaut[representation])?;
                p_gen_data.NormAtomsTaut[representation] = SourceMutPointer::null();
            }
            if !p_gen_data.NormAtomsNontaut[representation].is_null() {
                inchi_free(heap, p_gen_data.NormAtomsNontaut[representation])?;
                p_gen_data.NormAtomsNontaut[representation] = SourceMutPointer::null();
            }
        }

        FreeAllINChIArrays(
            heap,
            &mut control.pINChI,
            &mut control.pINChI_Aux,
            &mut control.StructData.num_components,
        )?;
        control.pINChI = std::array::from_fn(|_| SourceMutPointer::null());
        control.pINChI_Aux = std::array::from_fn(|_| SourceMutPointer::null());
        control.szTitle[0] = 0;
        control.StructData = STRUCT_DATA::default();
        Ok(())
    })();

    *heap
        .slice_mut(control_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = control;
    operation?;
    *p_results = inchi_Output::default();
    *p_gen_data = INCHIGEN_DATA::default();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{
        INCHI_IOS_STRING, NORM_ATOMS, PINChI_Aux2, PINChI2, SourceVoid, inchi_Atom, inp_ATOM, sp_ATOM,
    };

    fn assert_freed<T: 'static>(heap: &SourceHeap, pointer: SourceMutPointer<T>) {
        assert_eq!(
            heap.slice(pointer.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    fn assert_live<T: 'static>(heap: &SourceHeap, pointer: SourceMutPointer<T>) {
        assert!(heap.slice(pointer.as_const()).is_ok());
    }

    fn setup_build() -> InchiBuildMetadata<'static> {
        InchiBuildMetadata {
            compiler: "gcc",
            date: "Jan 01 1970",
            time: "00:00:00",
        }
    }

    fn setup_atom(element: &[u8]) -> inchi_Atom {
        let mut atom = inchi_Atom::default();
        for (target, source) in atom.elname.iter_mut().zip(element) {
            *target = *source as i8;
        }
        atom
    }

    fn setup_input(heap: &mut SourceHeap, atom: Option<inchi_Atom>, options: Option<&str>) -> inchi_Input {
        let atom = atom
            .map(|atom| heap.allocate_model_storage(vec![atom]).unwrap())
            .unwrap_or_else(SourceMutPointer::null);
        let options = options
            .map(|options| {
                heap.allocate_model_storage(
                    options
                        .bytes()
                        .chain(std::iter::once(0))
                        .map(|byte| byte as i8)
                        .collect(),
                )
                .unwrap()
            })
            .unwrap_or_else(SourceMutPointer::null);
        inchi_Input {
            atom,
            szOptions: options,
            num_atoms: i16::from(!atom.is_null()),
            ..inchi_Input::default()
        }
    }

    #[test]
    fn source_port__inchi_dll_a__inchigen_setup__line_260() {
        let stdout = SourceMutPointer::null();

        let mut null_heap = SourceHeap::default();
        let null_input = inchi_Input::default();
        let mut null_data = INCHIGEN_DATA::default();
        assert_eq!(
            INCHIGEN_Setup(
                &mut null_heap,
                SourceMutPointer::null(),
                Some(&mut null_data),
                Some(&null_input),
                stdout,
                setup_build(),
            ),
            Err(SourceHeapError::NullPointer)
        );

        let mut pseudo_heap = SourceHeap::default();
        let pseudo_handle = INCHIGEN_Create(&mut pseudo_heap).unwrap();
        let pseudo_input = setup_input(&mut pseudo_heap, Some(setup_atom(b"Zz\0")), None);
        let retained = pseudo_heap.allocate_model_storage(vec![19_i8]).unwrap();
        let mut pseudo_data = INCHIGEN_DATA {
            pStrErrStruct: [77; 256],
            num_components: [8, 9],
            NormAtomsNontaut: [retained.cast(), SourceMutPointer::null()],
            ..INCHIGEN_DATA::default()
        };
        assert_eq!(
            INCHIGEN_Setup(
                &mut pseudo_heap,
                pseudo_handle,
                Some(&mut pseudo_data),
                Some(&pseudo_input),
                stdout,
                setup_build(),
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_ERROR)
        );
        let pseudo_message = b"Pseudoatoms are not supported in current API mode\0";
        assert_eq!(
            &pseudo_data.pStrErrStruct[..pseudo_message.len()],
            &pseudo_message.iter().map(|byte| *byte as i8).collect::<Vec<_>>()
        );
        assert_eq!(pseudo_data.pStrErrStruct[pseudo_message.len()], 77);
        assert_eq!(pseudo_data.num_components, [0, 0]);
        assert_eq!(pseudo_data.NormAtomsNontaut[0], retained.cast());
        let pseudo_control = pseudo_heap
            .slice(pseudo_handle.cast::<INCHIGEN_CONTROL>().as_const())
            .unwrap()[0]
            .clone();
        assert_eq!(pseudo_control.num_err, 1);
        assert_eq!(pseudo_control.init_passed, 0);
        assert_eq!(pseudo_control.StructData.nStructReadError, 99);
        assert_eq!(pseudo_control.StructData.nErrorType, _IS_ERROR as i32);

        let mut missing_data_heap = SourceHeap::default();
        let missing_data_handle = INCHIGEN_Create(&mut missing_data_heap).unwrap();
        let missing_data_input = setup_input(&mut missing_data_heap, Some(setup_atom(b"C\0")), None);
        assert_eq!(
            INCHIGEN_Setup(
                &mut missing_data_heap,
                missing_data_handle,
                None,
                Some(&missing_data_input),
                stdout,
                setup_build(),
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_ERROR)
        );
        let missing_data_control = missing_data_heap
            .slice(missing_data_handle.cast::<INCHIGEN_CONTROL>().as_const())
            .unwrap()[0]
            .clone();
        assert_eq!(missing_data_control.num_err, 1);
        assert_eq!(missing_data_control.init_passed, 0);

        let mut allocation_heap = SourceHeap::default();
        let allocation_handle = INCHIGEN_Create(&mut allocation_heap).unwrap();
        let allocation_input = setup_input(&mut allocation_heap, Some(setup_atom(b"C\0")), Some("-AuxNone"));
        let mut allocation_data = INCHIGEN_DATA {
            pStrErrStruct: [41; 256],
            num_components: [3, 4],
            ..INCHIGEN_DATA::default()
        };
        allocation_heap.fail_after_allocations(0);
        assert_eq!(
            INCHIGEN_Setup(
                &mut allocation_heap,
                allocation_handle,
                Some(&mut allocation_data),
                Some(&allocation_input),
                stdout,
                setup_build(),
            ),
            Ok(_IS_FATAL as i32)
        );
        assert_eq!(allocation_heap.source_allocation_calls(), 1);
        assert_eq!(allocation_data, INCHIGEN_DATA::default());
        let allocation_control = allocation_heap
            .slice(allocation_handle.cast::<INCHIGEN_CONTROL>().as_const())
            .unwrap()[0]
            .clone();
        assert_eq!(allocation_control.num_err, 0);
        assert_eq!(allocation_control.InpParms, INPUT_PARMS::default());

        let mut no_input_heap = SourceHeap::default();
        let no_input_handle = INCHIGEN_Create(&mut no_input_heap).unwrap();
        let no_input = inchi_Input::default();
        let mut no_input_data = INCHIGEN_DATA {
            pStrErrStruct: [22; 256],
            num_components: [5, 6],
            ..INCHIGEN_DATA::default()
        };
        no_input_heap.trace_source_allocations();
        assert_eq!(
            INCHIGEN_Setup(
                &mut no_input_heap,
                no_input_handle,
                Some(&mut no_input_data),
                Some(&no_input),
                stdout,
                setup_build(),
            ),
            Ok(_IS_EOF)
        );
        assert_eq!(no_input_data, INCHIGEN_DATA::default());
        let no_input_control = no_input_heap
            .slice(no_input_handle.cast::<INCHIGEN_CONTROL>().as_const())
            .unwrap()[0]
            .clone();
        assert_eq!(no_input_control.inchi_file[1], Default::default());
        assert_eq!(no_input_control.num_err, 0);
        assert_eq!(no_input_control.init_passed, 0);
        assert!(no_input_heap.live_source_allocation_count() > 2);

        let mut help_heap = SourceHeap::default();
        let help_handle = INCHIGEN_Create(&mut help_heap).unwrap();
        let help_input = setup_input(&mut help_heap, Some(setup_atom(b"C\0")), Some("-HeLp"));
        let mut help_data = INCHIGEN_DATA {
            pStrErrStruct: [33; 256],
            ..INCHIGEN_DATA::default()
        };
        help_heap.trace_source_allocations();
        assert_eq!(
            INCHIGEN_Setup(
                &mut help_heap,
                help_handle,
                Some(&mut help_data),
                Some(&help_input),
                stdout,
                setup_build(),
            ),
            Ok(_IS_EOF)
        );
        assert_eq!(help_data, INCHIGEN_DATA::default());
        assert!(help_heap.source_allocation_calls() >= 2);
        let help_live_before_destroy = help_heap.live_source_allocation_count();
        assert!(help_live_before_destroy >= 4);
        INCHIGEN_Destroy(&mut help_heap, help_handle).unwrap();
        assert_eq!(help_heap.live_source_allocation_count(), help_live_before_destroy - 2);

        let mut empty_heap = SourceHeap::default();
        let empty_handle = INCHIGEN_Create(&mut empty_heap).unwrap();
        let empty_input = setup_input(&mut empty_heap, None, Some("-AuxNone"));
        let mut empty_data = INCHIGEN_DATA::default();
        assert_eq!(
            INCHIGEN_Setup(
                &mut empty_heap,
                empty_handle,
                Some(&mut empty_data),
                Some(&empty_input),
                stdout,
                setup_build(),
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_ERROR)
        );
        assert_eq!(empty_data.pStrErrStruct[0], b'E' as i8);
        let empty_control = empty_heap
            .slice(empty_handle.cast::<INCHIGEN_CONTROL>().as_const())
            .unwrap()[0]
            .clone();
        assert_eq!(empty_control.num_err, 1);
        assert_eq!(empty_control.init_passed, 0);
        assert_eq!(empty_control.OrigInpData.num_inp_atoms, 0);

        let mut success_heap = SourceHeap::default();
        let success_handle = INCHIGEN_Create(&mut success_heap).unwrap();
        let success_input = setup_input(&mut success_heap, Some(setup_atom(b"C\0")), None);
        let mut success_data = INCHIGEN_DATA {
            pStrErrStruct: [55; 256],
            num_components: [7, 8],
            ..INCHIGEN_DATA::default()
        };
        assert_eq!(
            INCHIGEN_Setup(
                &mut success_heap,
                success_handle,
                Some(&mut success_data),
                Some(&success_input),
                stdout,
                setup_build(),
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY)
        );
        let success_control = success_heap
            .slice(success_handle.cast::<INCHIGEN_CONTROL>().as_const())
            .unwrap()[0]
            .clone();
        assert_eq!(success_control.init_passed, 1);
        assert_eq!(success_control.num_err, 0);
        assert_eq!(success_control.num_inp, 1);
        assert_eq!(success_control.OrigInpData.num_inp_atoms, 1);
        assert_eq!(success_data.num_components, success_control.StructData.num_components);
        assert_eq!(success_data.pStrErrStruct[0], 0);
        assert!(success_control.InpParms.bNoStructLabels != 0);
        assert!(success_control.InpParms.pSdfLabel.is_null());
        assert!(success_control.InpParms.pSdfValue.is_null());
    }

    #[test]
    fn source_port__inchi_dll_a__stdinchigen_setup__line_185() {
        let stdout = SourceMutPointer::null();
        let ignore_unknown = u64::from(REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU);
        let cases = [
            (INCHI_OUT_SAVEOPT as i32, 0_u64, ignore_unknown),
            (0, u64::from(TG_FLAG_RECONNECT_COORD), ignore_unknown),
            (0, 0, ignore_unknown | u64::from(REQ_MODE_BASIC)),
            (0, 0, ignore_unknown | u64::from(REQ_MODE_RELATIVE_STEREO)),
            (0, 0, ignore_unknown | u64::from(REQ_MODE_RACEMIC_STEREO)),
            (0, 0, ignore_unknown | u64::from(REQ_MODE_CHIR_FLG_STEREO)),
            (0, 0, ignore_unknown | u64::from(REQ_MODE_DIFF_UU_STEREO)),
            (0, 0, 0),
            (0, u64::from(TG_FLAG_KETO_ENOL_TAUT), ignore_unknown),
            (0, u64::from(TG_FLAG_1_5_TAUT), ignore_unknown),
        ];
        for (output_options, taut_flags, mode) in cases {
            let mut heap = SourceHeap::default();
            let handle = INCHIGEN_Create(&mut heap).unwrap();
            let pointer = handle.cast::<INCHIGEN_CONTROL>();
            let control = &mut heap.slice_mut(pointer).unwrap()[0];
            control.InpParms.bINChIOutputOptions = output_options;
            control.InpParms.bTautFlags = taut_flags;
            control.InpParms.nMode = mode;
            let input = setup_input(&mut heap, Some(setup_atom(b"Zz\0")), None);
            let mut data = INCHIGEN_DATA {
                pStrErrStruct: [61; 256],
                ..INCHIGEN_DATA::default()
            };
            assert_eq!(
                STDINCHIGEN_Setup(&mut heap, handle, Some(&mut data), Some(&input), stdout, setup_build(),),
                Ok(_IS_WARNING as i32),
                "output={output_options} taut={taut_flags} mode={mode}"
            );
            let control = heap.slice(pointer.as_const()).unwrap()[0].clone();
            assert_eq!(control.InpParms.bINChIOutputOptions & INCHI_OUT_SAVEOPT as i32, 0);
            assert_ne!(control.InpParms.bINChIOutputOptions & INCHI_OUT_STDINCHI as i32, 0);
            assert_eq!(
                control.InpParms.bTautFlags
                    & u64::from(TG_FLAG_RECONNECT_COORD | TG_FLAG_KETO_ENOL_TAUT | TG_FLAG_1_5_TAUT),
                0
            );
            assert_eq!(
                control.InpParms.nMode
                    & u64::from(
                        REQ_MODE_BASIC
                            | REQ_MODE_RELATIVE_STEREO
                            | REQ_MODE_RACEMIC_STEREO
                            | REQ_MODE_CHIR_FLG_STEREO
                            | REQ_MODE_DIFF_UU_STEREO
                    ),
                0
            );
            assert_eq!(control.InpParms.nMode & ignore_unknown, ignore_unknown);
            assert_eq!(control.num_err, 1);
            assert_eq!(control.init_passed, 0);
            let message = b"Pseudoatoms are not supported in current API mode\0";
            assert_eq!(
                &data.pStrErrStruct[..message.len()],
                &message.iter().map(|byte| *byte as i8).collect::<Vec<_>>()
            );
            assert_eq!(data.pStrErrStruct[message.len()], 61);
        }

        let mut unchanged_heap = SourceHeap::default();
        let unchanged_handle = INCHIGEN_Create(&mut unchanged_heap).unwrap();
        let unchanged_pointer = unchanged_handle.cast::<INCHIGEN_CONTROL>();
        unchanged_heap.slice_mut(unchanged_pointer).unwrap()[0].InpParms.nMode = ignore_unknown;
        let unchanged_input = setup_input(&mut unchanged_heap, Some(setup_atom(b"Zz\0")), None);
        let mut unchanged_data = INCHIGEN_DATA::default();
        assert_eq!(
            STDINCHIGEN_Setup(
                &mut unchanged_heap,
                unchanged_handle,
                Some(&mut unchanged_data),
                Some(&unchanged_input),
                stdout,
                setup_build(),
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_ERROR)
        );
        let unchanged = &unchanged_heap.slice(unchanged_pointer.as_const()).unwrap()[0];
        assert_eq!(unchanged.InpParms.nMode, ignore_unknown);
        assert_eq!(unchanged.InpParms.bINChIOutputOptions, INCHI_OUT_STDINCHI as i32);

        let mut missing_heap = SourceHeap::default();
        let missing_handle = INCHIGEN_Create(&mut missing_heap).unwrap();
        let missing_pointer = missing_handle.cast::<INCHIGEN_CONTROL>();
        missing_heap.slice_mut(missing_pointer).unwrap()[0]
            .InpParms
            .bINChIOutputOptions = INCHI_OUT_SAVEOPT as i32;
        let missing_input = setup_input(&mut missing_heap, Some(setup_atom(b"C\0")), None);
        assert_eq!(
            STDINCHIGEN_Setup(
                &mut missing_heap,
                missing_handle,
                None,
                Some(&missing_input),
                stdout,
                setup_build(),
            ),
            Err(SourceHeapError::NullPointer)
        );
        let missing = &missing_heap.slice(missing_pointer.as_const()).unwrap()[0];
        assert_eq!(missing.InpParms.bINChIOutputOptions & INCHI_OUT_SAVEOPT as i32, 0);
        assert_ne!(missing.InpParms.bINChIOutputOptions & INCHI_OUT_STDINCHI as i32, 0);

        let mut success_heap = SourceHeap::default();
        let success_handle = INCHIGEN_Create(&mut success_heap).unwrap();
        let success_input = setup_input(&mut success_heap, Some(setup_atom(b"C\0")), None);
        let mut success_data = INCHIGEN_DATA::default();
        assert_eq!(
            STDINCHIGEN_Setup(
                &mut success_heap,
                success_handle,
                Some(&mut success_data),
                Some(&success_input),
                stdout,
                setup_build(),
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY)
        );
        let success = &success_heap
            .slice(success_handle.cast::<INCHIGEN_CONTROL>().as_const())
            .unwrap()[0];
        assert_eq!(success.init_passed, 1);
        assert_ne!(success.InpParms.bINChIOutputOptions & INCHI_OUT_STDINCHI as i32, 0);
    }

    #[test]
    fn source_port__inchi_dll_a__stdinchigen_donormalization__line_435() {
        let mut null_heap = SourceHeap::default();
        let mut null_data = INCHIGEN_DATA::default();
        assert_eq!(
            STDINCHIGEN_DoNormalization(&mut null_heap, SourceMutPointer::null(), Some(&mut null_data),),
            Err(SourceHeapError::NullPointer)
        );

        let mut heap = SourceHeap::default();
        let handle = INCHIGEN_Create(&mut heap).unwrap();
        let mut data = INCHIGEN_DATA::default();
        assert_eq!(
            STDINCHIGEN_DoNormalization(&mut heap, handle, Some(&mut data)),
            Ok(_IS_ERROR as i32)
        );
        let control = &heap.slice(handle.cast::<INCHIGEN_CONTROL>().as_const()).unwrap()[0];
        assert_eq!(control.StructData.nStructReadError, 99);
        let message = b"InChI generator not initialized\0";
        assert_eq!(
            &data.pStrErrStruct[..message.len()],
            &message.iter().map(|byte| *byte as i8).collect::<Vec<_>>()
        );
    }

    #[test]
    fn source_port__inchi_dll_a__inchigen_donormalization__line_442() {
        let mut null_heap = SourceHeap::default();
        let mut null_data = INCHIGEN_DATA::default();
        assert_eq!(
            INCHIGEN_DoNormalization(&mut null_heap, SourceMutPointer::null(), Some(&mut null_data),),
            Err(SourceHeapError::NullPointer)
        );

        let mut uninitialized_heap = SourceHeap::default();
        let uninitialized_handle = INCHIGEN_Create(&mut uninitialized_heap).unwrap();
        let mut uninitialized_data = INCHIGEN_DATA {
            pStrErrStruct: [77; 256],
            ..INCHIGEN_DATA::default()
        };
        assert_eq!(
            INCHIGEN_DoNormalization(
                &mut uninitialized_heap,
                uninitialized_handle,
                Some(&mut uninitialized_data),
            ),
            Ok(_IS_ERROR as i32)
        );
        let uninitialized = &uninitialized_heap
            .slice(uninitialized_handle.cast::<INCHIGEN_CONTROL>().as_const())
            .unwrap()[0];
        assert_eq!(uninitialized.StructData.nStructReadError, 99);
        assert_eq!(uninitialized.StructData.nErrorType, _IS_ERROR as i32);
        assert_eq!(uninitialized.norm_passed, 0);
        let init_message = b"InChI generator not initialized\0";
        assert_eq!(
            &uninitialized_data.pStrErrStruct[..init_message.len()],
            &init_message.iter().map(|byte| *byte as i8).collect::<Vec<_>>()
        );
        assert_eq!(uninitialized_data.pStrErrStruct[init_message.len()], 77);

        let mut normalization_heap = SourceHeap::default();
        let atoms = normalization_heap
            .allocate_model_storage(vec![inp_ATOM {
                elname: [b'C' as i8, 0, 0, 0, 0, 0],
                el_number: 6,
                orig_at_number: 1,
                component: 1,
                ..inp_ATOM::default()
            }])
            .unwrap();
        let base_lengths = normalization_heap.allocate_model_storage(vec![1_u16]).unwrap();
        let reconnected_lengths = normalization_heap.allocate_model_storage(vec![1_u16]).unwrap();
        let original = ORIG_ATOM_DATA {
            at: atoms,
            num_inp_atoms: 1,
            num_components: 1,
            bDisconnectCoord: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let prepared_base = ORIG_ATOM_DATA {
            nCurAtLen: base_lengths,
            ..original.clone()
        };
        let prepared_reconnected = ORIG_ATOM_DATA {
            nCurAtLen: reconnected_lengths,
            ..original.clone()
        };
        let mut normalization_control = INCHIGEN_CONTROL::default();
        normalization_control.init_passed = 1;
        normalization_control.InpParms.nMode = u64::from(REQ_MODE_BASIC);
        normalization_control.InpParms.bINChIOutputOptions = INCHI_OUT_NO_AUX_INFO as i32;
        normalization_control.InpParms.bTautFlags =
            u64::from(crate::source_types::TG_FLAG_H_ALREADY_REMOVED | TG_FLAG_RECONNECT_COORD);
        normalization_control.OrigInpData = original;
        normalization_control.PrepInpData = [prepared_base, prepared_reconnected];
        normalization_control.StructData.bTautFlagsDone[INCHI_BAS as usize] = u64::from(TG_FLAG_DISCONNECT_COORD_DONE);
        let normalization_pointer = normalization_heap
            .allocate_model_storage(vec![normalization_control])
            .unwrap();
        let normalization_handle = normalization_pointer.cast::<SourceVoid>();
        let mut normalization_data = INCHIGEN_DATA::default();
        assert_eq!(
            INCHIGEN_DoNormalization(
                &mut normalization_heap,
                normalization_handle,
                Some(&mut normalization_data),
            ),
            Ok(0)
        );
        let normalized = &normalization_heap.slice(normalization_pointer.as_const()).unwrap()[0];
        assert_eq!(normalized.norm_passed, 1);
        assert_eq!(normalization_data.num_components, [1, 1]);
        for representation in 0..INCHI_NUM as usize {
            assert!(!normalization_data.NormAtomsNontaut[representation].is_null());
            assert!(!normalization_data.NormAtomsTaut[representation].is_null());
            let public = &normalization_heap
                .slice(normalization_data.NormAtomsNontaut[representation].as_const())
                .unwrap()[0];
            let internal = &normalization_heap
                .slice(normalized.InpNormAtData[representation].as_const())
                .unwrap()[0];
            assert_eq!(public.at, internal.at.cast());
            assert_eq!(public.num_at, internal.num_at);
            assert_eq!(public.bNormalizationFlags, internal.bNormalizationFlags);
        }

        let mut warning_heap = SourceHeap::default();
        let warning_handle = INCHIGEN_Create(&mut warning_heap).unwrap();
        let warning_pointer = warning_handle.cast::<INCHIGEN_CONTROL>();
        let warning_atoms = warning_heap
            .allocate_model_storage(vec![inp_ATOM {
                elname: [b'C' as i8, 0, 0, 0, 0, 0],
                el_number: 6,
                orig_at_number: 1,
                ..inp_ATOM::default()
            }])
            .unwrap();
        let warning_entry = INP_ATOM_DATA {
            num_at: 31,
            bTautomeric: 1,
            bNormalizationFlags: u64::from(
                (FLAG_NORM_CONSIDER_TAUT & !FLAG_PROTON_CHARGE_CANCEL) | FLAG_PROTON_CHARGE_CANCEL,
            ),
            ..INP_ATOM_DATA::default()
        };
        let warning_non = warning_heap
            .allocate_model_storage(vec![warning_entry.clone(); 2])
            .unwrap();
        let warning_taut = warning_heap
            .allocate_model_storage(vec![warning_entry.clone(); 2])
            .unwrap();
        let warning_non_output = warning_heap
            .allocate_model_storage(vec![NORM_ATOMS::default(); 2])
            .unwrap();
        let warning_taut_output = warning_heap
            .allocate_model_storage(vec![NORM_ATOMS::default(); 2])
            .unwrap();
        {
            let warning_control = &mut warning_heap.slice_mut(warning_pointer).unwrap()[0];
            warning_control.init_passed = 1;
            warning_control.InpParms.bINChIOutputOptions = INCHI_OUT_SDFILE_ONLY as i32;
            warning_control.OrigInpData = ORIG_ATOM_DATA {
                at: warning_atoms,
                num_inp_atoms: 1,
                num_components: 1,
                ..ORIG_ATOM_DATA::default()
            };
            warning_control.StructData.num_components = [2, 0];
            warning_control.InpNormAtData[0] = warning_non;
            warning_control.InpNormTautData[0] = warning_taut;
        }
        let mut warning_data = INCHIGEN_DATA {
            pStrErrStruct: [93; 256],
            NormAtomsNontaut: [warning_non_output, SourceMutPointer::null()],
            NormAtomsTaut: [warning_taut_output, SourceMutPointer::null()],
            ..INCHIGEN_DATA::default()
        };
        assert_eq!(
            INCHIGEN_DoNormalization(&mut warning_heap, warning_handle, Some(&mut warning_data),),
            Ok(0)
        );
        let warning_control = &warning_heap.slice(warning_pointer.as_const()).unwrap()[0];
        assert_eq!(warning_control.norm_passed, 1);
        assert!(!warning_control.inchi_file[0].s.pStr.is_null());
        let warning_message = b"Proton(s) added/removed; Charges neutralized\0";
        assert_eq!(
            &warning_data.pStrErrStruct[..warning_message.len()],
            &warning_message.iter().map(|byte| *byte as i8).collect::<Vec<_>>()
        );
        assert_eq!(warning_data.pStrErrStruct[warning_message.len()], 93);
        assert_eq!(warning_data.num_components, [2, 0]);
        assert_eq!(
            warning_heap.slice(warning_data.NormAtomsNontaut[0].as_const()).unwrap()[0].num_at,
            31
        );
        assert_eq!(
            warning_heap.slice(warning_data.NormAtomsTaut[0].as_const()).unwrap()[1].bNormalizationFlags,
            warning_entry.bNormalizationFlags
        );
    }

    #[test]
    fn source_port__inchi_dll_a__stdinchigen_docanonicalization__line_689() {
        let mut null_heap = SourceHeap::default();
        let mut null_data = INCHIGEN_DATA::default();
        assert_eq!(
            STDINCHIGEN_DoCanonicalization(&mut null_heap, SourceMutPointer::null(), Some(&mut null_data),),
            Err(SourceHeapError::NullPointer)
        );

        let mut wrapper_heap = SourceHeap::default();
        let wrapper_handle = INCHIGEN_Create(&mut wrapper_heap).unwrap();
        let mut wrapper_data = INCHIGEN_DATA {
            pStrErrStruct: [77; 256],
            ..INCHIGEN_DATA::default()
        };
        let wrapper_result = STDINCHIGEN_DoCanonicalization(&mut wrapper_heap, wrapper_handle, Some(&mut wrapper_data));

        let mut direct_heap = SourceHeap::default();
        let direct_handle = INCHIGEN_Create(&mut direct_heap).unwrap();
        let mut direct_data = INCHIGEN_DATA {
            pStrErrStruct: [77; 256],
            ..INCHIGEN_DATA::default()
        };
        let direct_result = INCHIGEN_DoCanonicalization(&mut direct_heap, direct_handle, Some(&mut direct_data));

        assert_eq!(wrapper_result, direct_result);
        assert_eq!(wrapper_data, direct_data);
        let wrapper_control = &wrapper_heap
            .slice(wrapper_handle.cast::<INCHIGEN_CONTROL>().as_const())
            .unwrap()[0];
        let direct_control = &direct_heap
            .slice(direct_handle.cast::<INCHIGEN_CONTROL>().as_const())
            .unwrap()[0];
        assert_eq!(wrapper_control.StructData, direct_control.StructData);
        assert_eq!(wrapper_control.canon_passed, direct_control.canon_passed);
    }

    #[test]
    fn source_port__inchi_dll_a__inchigen_docanonicalization__line_697() {
        let mut unnormalized_heap = SourceHeap::default();
        let unnormalized = INCHIGEN_Create(&mut unnormalized_heap).unwrap();
        let mut unnormalized_data = INCHIGEN_DATA {
            pStrErrStruct: [77; 256],
            ..INCHIGEN_DATA::default()
        };
        assert_eq!(
            INCHIGEN_DoCanonicalization(&mut unnormalized_heap, unnormalized, Some(&mut unnormalized_data),),
            Ok(tagRetValGetINCHI_inchi_Ret_ERROR)
        );
        let message = b"Got non-normalized structure\0";
        assert_eq!(
            &unnormalized_data.pStrErrStruct[..message.len()],
            &message.iter().map(|byte| *byte as i8).collect::<Vec<_>>()
        );
        assert_eq!(unnormalized_data.pStrErrStruct[message.len()], 77);

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![inp_ATOM {
                elname: [b'C' as i8, 0, 0, 0, 0, 0],
                el_number: 6,
                orig_at_number: 1,
                component: 1,
                ..inp_ATOM::default()
            }])
            .unwrap();
        let lengths = heap.allocate_model_storage(vec![1_u16]).unwrap();
        let original = ORIG_ATOM_DATA {
            at: atoms,
            num_inp_atoms: 1,
            num_components: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let mut control = INCHIGEN_CONTROL::default();
        control.init_passed = 1;
        control.InpParms.nMode = u64::from(REQ_MODE_BASIC | crate::source_types::REQ_MODE_NON_ISO);
        control.InpParms.bINChIOutputOptions = INCHI_OUT_NO_AUX_INFO as i32;
        control.InpParms.bTautFlags = u64::from(crate::source_types::TG_FLAG_H_ALREADY_REMOVED);
        control.OrigInpData = original.clone();
        control.PrepInpData[INCHI_BAS as usize] = ORIG_ATOM_DATA {
            nCurAtLen: lengths,
            ..original
        };
        let control_pointer = heap.allocate_model_storage(vec![control]).unwrap();
        let handle = control_pointer.cast::<SourceVoid>();
        let mut data = INCHIGEN_DATA::default();
        assert_eq!(
            INCHIGEN_DoNormalization(&mut heap, handle, Some(&mut data)),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(
            INCHIGEN_DoCanonicalization(&mut heap, handle, Some(&mut data)),
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY)
        );
        let canonicalized = &heap.slice(control_pointer.as_const()).unwrap()[0];
        assert_eq!(canonicalized.canon_passed, 1);
        assert_eq!(data.num_components, [1, 0]);
        let row = heap.slice(canonicalized.pINChI[INCHI_BAS as usize].as_const()).unwrap()[0];
        let inchi = &heap.slice(row[TAUT_NON as usize].as_const()).unwrap()[0];
        assert_eq!((inchi.nErrorCode, inchi.nNumberOfAtoms), (0, 1));
        assert_eq!(heap.slice(inchi.nAtom.as_const()).unwrap(), &[6]);
    }

    #[test]
    fn source_port__inchi_dll_a__stdinchigen_doserialization__line_859() {
        let mut wrapper_heap = SourceHeap::default();
        let wrapper_handle = INCHIGEN_Create(&mut wrapper_heap).unwrap();
        let mut wrapper_data = INCHIGEN_DATA {
            pStrErrStruct: [71; 256],
            num_components: [3, 5],
            ..INCHIGEN_DATA::default()
        };
        let mut wrapper_results = inchi_Output::default();
        let wrapper_result = STDINCHIGEN_DoSerialization(
            &mut wrapper_heap,
            wrapper_handle,
            Some(&mut wrapper_data),
            Some(&mut wrapper_results),
        );

        let mut direct_heap = SourceHeap::default();
        let direct_handle = INCHIGEN_Create(&mut direct_heap).unwrap();
        let mut direct_data = INCHIGEN_DATA {
            pStrErrStruct: [71; 256],
            num_components: [3, 5],
            ..INCHIGEN_DATA::default()
        };
        let mut direct_results = inchi_Output::default();
        let direct_result = INCHIGEN_DoSerialization(
            &mut direct_heap,
            direct_handle,
            Some(&mut direct_data),
            Some(&mut direct_results),
        );

        assert_eq!(wrapper_result, direct_result);
        assert_eq!(wrapper_data, direct_data);
        assert_eq!(wrapper_results, direct_results);
        let wrapper_control = &wrapper_heap
            .slice(wrapper_handle.cast::<INCHIGEN_CONTROL>().as_const())
            .unwrap()[0];
        let direct_control = &direct_heap
            .slice(direct_handle.cast::<INCHIGEN_CONTROL>().as_const())
            .unwrap()[0];
        assert_eq!(wrapper_control, direct_control);
    }

    #[test]
    fn source_port__inchi_dll_a__inchigen_doserialization__line_869() {
        let mut rejected_heap = SourceHeap::default();
        let rejected_handle = INCHIGEN_Create(&mut rejected_heap).unwrap();
        let retained_log = rejected_heap
            .allocate_model_storage(vec![b'l' as i8, b'o' as i8, b'g' as i8, 0])
            .unwrap();
        let released_path = rejected_heap.allocate_model_storage(vec![b'p' as i8, 0]).unwrap();
        {
            let rejected_control = &mut rejected_heap
                .slice_mut(rejected_handle.cast::<INCHIGEN_CONTROL>())
                .unwrap()[0];
            rejected_control.inchi_file[1].s.pStr = retained_log;
            rejected_control.inchi_file[1].s.nUsedLength = 3;
            rejected_control.inchi_file[1].s.nAllocatedLength = 4;
            rejected_control.InpParms.path[0] = released_path.as_const();
        }
        let mut rejected_data = INCHIGEN_DATA {
            pStrErrStruct: [77; 256],
            num_components: [3, 4],
            ..INCHIGEN_DATA::default()
        };
        let sentinel = rejected_heap.allocate_model_storage(vec![b's' as i8, 0]).unwrap();
        let mut rejected_results = inchi_Output {
            szInChI: sentinel,
            szAuxInfo: sentinel,
            szMessage: sentinel,
            szLog: sentinel,
        };
        assert_eq!(
            INCHIGEN_DoSerialization(
                &mut rejected_heap,
                rejected_handle,
                Some(&mut rejected_data),
                Some(&mut rejected_results),
            ),
            Ok(_IS_ERROR as i32)
        );
        let rejection = b"Got non-canonicalized structure\0";
        assert_eq!(
            &rejected_data.pStrErrStruct[..rejection.len()],
            &rejection.iter().map(|byte| *byte as i8).collect::<Vec<_>>()
        );
        assert_eq!(rejected_results.szInChI, SourceMutPointer::null());
        assert_eq!(rejected_results.szAuxInfo, SourceMutPointer::null());
        assert_eq!(rejected_results.szMessage, SourceMutPointer::null());
        assert_eq!(rejected_results.szLog, retained_log);
        assert_freed(&rejected_heap, released_path);
        assert_live(&rejected_heap, retained_log);
        assert_live(&rejected_heap, sentinel);
        let rejected_control = &rejected_heap
            .slice(rejected_handle.cast::<INCHIGEN_CONTROL>().as_const())
            .unwrap()[0];
        assert_eq!(rejected_control.StructData.nStructReadError, 99);
        assert_eq!(rejected_control.StructData.nErrorType, _IS_ERROR as i32);
        assert_eq!(rejected_control.InpParms.path[0], SourceConstPointer::null());

        let mut heap = SourceHeap::default();
        let handle = INCHIGEN_Create(&mut heap).unwrap();
        let mut methane_carbon = setup_atom(b"C\0");
        methane_carbon.num_iso_H[0] = -1;
        let input = setup_input(&mut heap, Some(methane_carbon), Some("-AuxNone"));
        let mut data = INCHIGEN_DATA::default();
        assert_eq!(
            INCHIGEN_Setup(
                &mut heap,
                handle,
                Some(&mut data),
                Some(&input),
                SourceMutPointer::null(),
                setup_build(),
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY)
        );
        {
            let control = &mut heap.slice_mut(handle.cast::<INCHIGEN_CONTROL>()).unwrap()[0];
            control.InpParms.bINChIOutputOptions |= INCHI_OUT_PLAIN_TEXT as i32;
        }
        assert_eq!(
            INCHIGEN_DoNormalization(&mut heap, handle, Some(&mut data)),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(
            INCHIGEN_DoCanonicalization(&mut heap, handle, Some(&mut data)),
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY)
        );
        let mut results = inchi_Output::default();
        assert_eq!(
            INCHIGEN_DoSerialization(&mut heap, handle, Some(&mut data), Some(&mut results),),
            Ok(0)
        );
        assert!(!results.szInChI.is_null());
        assert!(results.szAuxInfo.is_null());
        let inchi_length = heap
            .slice(results.szInChI.as_const())
            .unwrap()
            .iter()
            .position(|byte| *byte == 0)
            .unwrap();
        let inchi = heap.slice(results.szInChI.as_const()).unwrap()[..inchi_length]
            .iter()
            .map(|byte| *byte as u8)
            .collect::<Vec<_>>();
        assert_eq!(inchi, b"InChI=1S/CH4/h1H4");
        let serialized_control = &heap.slice(handle.cast::<INCHIGEN_CONTROL>().as_const()).unwrap()[0];
        assert!(serialized_control.inchi_file[0].s.pStr.is_null());
        assert_eq!(serialized_control.StructData.num_components, data.num_components);
    }

    #[test]
    fn source_port__inchi_dll_a__inchigen_reset__line_1134() {
        let mut null_heap = SourceHeap::default();
        assert_eq!(
            INCHIGEN_Reset(&mut null_heap, SourceMutPointer::null(), None, None,),
            Err(SourceHeapError::NullPointer)
        );

        let inchi = null_heap.allocate_model_storage(vec![1_i8, 0]).unwrap();
        let aux = null_heap.allocate_model_storage(vec![2_i8, 0]).unwrap();
        let log = null_heap.allocate_model_storage(vec![3_i8, 0]).unwrap();
        let message = null_heap.allocate_model_storage(vec![4_i8, 0]).unwrap();
        let normalized = null_heap.allocate_model_storage(vec![NORM_ATOMS::default()]).unwrap();
        let mut results = inchi_Output {
            szInChI: inchi,
            szAuxInfo: aux,
            szMessage: message,
            szLog: log,
        };
        let mut data = INCHIGEN_DATA {
            pStrErrStruct: [7; 256],
            num_components: [8, 9],
            NormAtomsNontaut: [normalized, SourceMutPointer::null()],
            NormAtomsTaut: [SourceMutPointer::null(), normalized],
        };
        assert_eq!(
            INCHIGEN_Reset(
                &mut null_heap,
                SourceMutPointer::null(),
                Some(&mut data),
                Some(&mut results),
            ),
            Ok(())
        );
        assert_eq!(results, inchi_Output::default());
        assert_eq!(data, INCHIGEN_DATA::default());
        assert_freed(&null_heap, inchi);
        assert_freed(&null_heap, log);
        assert_freed(&null_heap, message);
        assert_live(&null_heap, aux);
        assert_live(&null_heap, normalized);

        let mut heap = SourceHeap::default();
        let mut control = INCHIGEN_CONTROL::default();
        control.init_passed = 17;
        control.norm_passed = 18;
        control.canon_passed = 19;
        control.ulTotalProcessingTime = 23;
        control.num_err = 29;
        control.num_inp = 31;
        control.szTitle[0] = b'A' as i8;
        control.szTitle[1] = b'B' as i8;

        let stream_buffers: [SourceMutPointer<i8>; 3] =
            std::array::from_fn(|index| heap.allocate_model_storage(vec![index as i8 + 1, 0]).unwrap());
        for (stream, buffer) in control.inchi_file.iter_mut().zip(stream_buffers) {
            stream.type_ = 99;
            stream.s = INCHI_IOS_STRING {
                pStr: buffer,
                nUsedLength: 1,
                nAllocatedLength: 2,
                nPtr: 7,
            };
        }
        let retained_buffer = heap.allocate_model_storage(vec![b'x' as i8, b'y' as i8, 0]).unwrap();
        control.strbuf_container = INCHI_IOS_STRING {
            pStr: retained_buffer,
            nUsedLength: 2,
            nAllocatedLength: 3,
            nPtr: 11,
        };

        let path = heap.allocate_model_storage(vec![b'p' as i8, 0]).unwrap();
        control.InpParms.path[0] = path.as_const();
        control.InpParms.num_paths = 1;
        control.InpParms.nMode = u64::MAX;

        let original_atom = heap.allocate_model_storage(vec![inp_ATOM::default()]).unwrap();
        control.OrigInpData.at = original_atom;
        control.OrigInpData.num_inp_atoms = 1;
        let prepared_atoms: [SourceMutPointer<inp_ATOM>; 2] =
            std::array::from_fn(|_| heap.allocate_model_storage(vec![inp_ATOM::default()]).unwrap());
        for (prepared, atoms) in control.PrepInpData.iter_mut().zip(prepared_atoms) {
            prepared.at = atoms;
            prepared.num_inp_atoms = 1;
        }
        let original_texts: [SourceMutPointer<i8>; 3] =
            std::array::from_fn(|index| heap.allocate_model_storage(vec![index as i8 + 40, 0]).unwrap());
        control.OrigStruct.szAtoms = original_texts[0];
        control.OrigStruct.szBonds = original_texts[1];
        control.OrigStruct.szCoord = original_texts[2];
        control.OrigStruct.num_atoms = 1;

        let composite_atoms = heap.allocate_model_storage(vec![inp_ATOM::default()]).unwrap();
        let composite_offsets = heap.allocate_model_storage(vec![1_u16]).unwrap();
        control.composite_norm_data[0][0].at = composite_atoms;
        control.composite_norm_data[0][0].nOffsetAtAndH = composite_offsets;
        control.composite_norm_data[0][0].num_at = 1;

        fn atom_data(
            heap: &mut SourceHeap,
        ) -> (
            crate::source_types::INP_ATOM_DATA,
            SourceMutPointer<inp_ATOM>,
            SourceMutPointer<inp_ATOM>,
        ) {
            let at = heap.allocate_model_storage(vec![inp_ATOM::default()]).unwrap();
            let fixed = heap.allocate_model_storage(vec![inp_ATOM::default()]).unwrap();
            (
                crate::source_types::INP_ATOM_DATA {
                    at,
                    at_fixed_bonds: fixed,
                    num_at: 1,
                    ..crate::source_types::INP_ATOM_DATA::default()
                },
                at,
                fixed,
            )
        }

        let (current_data, current_at, current_fixed) = atom_data(&mut heap);
        let current = heap.allocate_model_storage(vec![current_data]).unwrap();
        control.InpCurAtData[0] = current;
        let (norm_data, norm_at, norm_fixed) = atom_data(&mut heap);
        let norm = heap.allocate_model_storage(vec![norm_data]).unwrap();
        control.InpNormAtData[0] = norm;
        let (taut_data, taut_at, taut_fixed) = atom_data(&mut heap);
        let taut = heap.allocate_model_storage(vec![taut_data]).unwrap();
        control.InpNormTautData[0] = taut;

        let cti_taut = heap.allocate_model_storage(vec![sp_ATOM::default()]).unwrap();
        let cti_nontaut = heap.allocate_model_storage(vec![sp_ATOM::default()]).unwrap();
        let cti = heap
            .allocate_model_storage(vec![crate::source_types::COMPONENT_TREAT_INFO {
                at: [cti_nontaut, cti_taut],
                ..crate::source_types::COMPONENT_TREAT_INFO::default()
            }])
            .unwrap();
        control.cti[0] = cti;

        let inchi_array = heap.allocate_model_storage(vec![PINChI2::default()]).unwrap();
        let aux_array = heap.allocate_model_storage(vec![PINChI_Aux2::default()]).unwrap();
        let leaked_zero_count_inchi = heap.allocate_model_storage(vec![PINChI2::default()]).unwrap();
        control.pINChI[0] = inchi_array;
        control.pINChI_Aux[0] = aux_array;
        control.pINChI[1] = leaked_zero_count_inchi;
        control.StructData.num_components = [1, 0];

        let control_pointer = heap.allocate_model_storage(vec![control]).unwrap();
        let handle: SourceMutPointer<SourceVoid> = control_pointer.cast();
        let normal_taut = heap.allocate_model_storage(vec![NORM_ATOMS::default()]).unwrap();
        let normal_nontaut = heap.allocate_model_storage(vec![NORM_ATOMS::default()]).unwrap();
        let mut data = INCHIGEN_DATA {
            pStrErrStruct: [33; 256],
            num_components: [1, 0],
            NormAtomsNontaut: [normal_nontaut, SourceMutPointer::null()],
            NormAtomsTaut: [normal_taut, SourceMutPointer::null()],
        };
        let result_inchi = heap.allocate_model_storage(vec![51_i8, 0]).unwrap();
        let result_aux = heap.allocate_model_storage(vec![52_i8, 0]).unwrap();
        let result_log = heap.allocate_model_storage(vec![53_i8, 0]).unwrap();
        let result_message = heap.allocate_model_storage(vec![54_i8, 0]).unwrap();
        let mut results = inchi_Output {
            szInChI: result_inchi,
            szAuxInfo: result_aux,
            szMessage: result_message,
            szLog: result_log,
        };

        assert_eq!(
            INCHIGEN_Reset(&mut heap, handle, Some(&mut data), Some(&mut results)),
            Ok(())
        );
        assert_eq!(results, inchi_Output::default());
        assert_eq!(data, INCHIGEN_DATA::default());

        for pointer in stream_buffers {
            assert_freed(&heap, pointer);
        }
        for pointer in [
            path,
            original_texts[0],
            original_texts[1],
            original_texts[2],
            result_inchi,
            result_log,
            result_message,
        ] {
            assert_freed(&heap, pointer);
        }
        for pointer in [
            original_atom,
            prepared_atoms[0],
            prepared_atoms[1],
            composite_atoms,
            current_at,
            current_fixed,
            norm_at,
            norm_fixed,
            taut_at,
            taut_fixed,
        ] {
            assert_freed(&heap, pointer);
        }
        assert_freed(&heap, composite_offsets);
        for pointer in [current, norm, taut] {
            assert_freed(&heap, pointer);
        }
        for pointer in [cti_taut, cti_nontaut] {
            assert_freed(&heap, pointer);
        }
        assert_freed(&heap, cti);
        assert_freed(&heap, inchi_array);
        assert_freed(&heap, aux_array);
        assert_freed(&heap, normal_taut);
        assert_freed(&heap, normal_nontaut);

        assert_live(&heap, retained_buffer);
        assert_eq!(heap.slice(retained_buffer.as_const()).unwrap()[0], 0);
        assert_live(&heap, leaked_zero_count_inchi);
        assert_live(&heap, result_aux);

        let reset = heap.slice(control_pointer.as_const()).unwrap()[0].clone();
        assert_eq!(reset.init_passed, 17);
        assert_eq!(reset.norm_passed, 18);
        assert_eq!(reset.canon_passed, 19);
        assert_eq!(reset.ulTotalProcessingTime, 23);
        assert_eq!(reset.num_err, 29);
        assert_eq!(reset.num_inp, 31);
        assert_eq!(reset.InpParms, INPUT_PARMS::default());
        assert_eq!(reset.OrigInpData, ORIG_ATOM_DATA::default());
        assert_eq!(reset.PrepInpData, std::array::from_fn(|_| ORIG_ATOM_DATA::default()));
        assert_eq!(reset.OrigStruct, ORIG_STRUCT::default());
        assert!(
            reset
                .composite_norm_data
                .iter()
                .flatten()
                .all(|data| *data == crate::source_types::COMP_ATOM_DATA::default())
        );
        assert!(reset.InpCurAtData.iter().all(|pointer| pointer.is_null()));
        assert!(reset.InpNormAtData.iter().all(|pointer| pointer.is_null()));
        assert!(reset.InpNormTautData.iter().all(|pointer| pointer.is_null()));
        assert!(reset.cti.iter().all(|pointer| pointer.is_null()));
        assert!(reset.pINChI.iter().all(|pointer| pointer.is_null()));
        assert!(reset.pINChI_Aux.iter().all(|pointer| pointer.is_null()));
        assert_eq!(reset.StructData, STRUCT_DATA::default());
        assert_eq!(reset.szTitle[0], 0);
        assert_eq!(reset.szTitle[1], b'B' as i8);
        assert_eq!(reset.strbuf_container.pStr, retained_buffer);
        assert_eq!(reset.strbuf_container.nUsedLength, 0);
        assert_eq!(reset.strbuf_container.nPtr, 0);
        assert_eq!(reset.strbuf_container.nAllocatedLength, 3);
        assert!(reset.inchi_file.iter().all(|stream| {
            stream.type_ == INCHI_IOS_TYPE_STRING as i32
                && stream.s == INCHI_IOS_STRING::default()
                && stream.f.is_null()
        }));
    }

    #[test]
    fn source_port__inchi_dll_a__inchigen_create__line_130() {
        for successful_allocations in [0, 1] {
            let mut heap = SourceHeap::default();
            heap.fail_after_allocations(successful_allocations);
            assert_eq!(INCHIGEN_Create(&mut heap), Ok(SourceMutPointer::null()));
            assert_eq!(heap.source_allocation_calls(), successful_allocations + 1);
            assert_eq!(heap.live_source_allocation_count(), 0);
        }

        let mut heap = SourceHeap::default();
        heap.trace_source_allocations();
        let handle = INCHIGEN_Create(&mut heap).unwrap();
        assert!(!handle.is_null());
        assert_eq!(heap.source_allocation_calls(), 2);
        assert_eq!(heap.live_source_allocation_count(), 2);
        assert_eq!(heap.live_source_allocations_of::<INCHIGEN_CONTROL>(), 1);
        assert_eq!(heap.live_source_allocations_of::<i8>(), 1);

        let control_pointer = handle.cast::<INCHIGEN_CONTROL>();
        let control = heap.slice(control_pointer.as_const()).unwrap()[0].clone();
        assert_eq!(control.InpParms, INPUT_PARMS::default());
        assert_eq!(control.StructData, STRUCT_DATA::default());
        assert_eq!(control.OrigInpData, ORIG_ATOM_DATA::default());
        assert_eq!(control.PrepInpData, std::array::from_fn(|_| ORIG_ATOM_DATA::default()));
        assert!(control.pINChI.iter().all(|pointer| pointer.is_null()));
        assert!(control.pINChI_Aux.iter().all(|pointer| pointer.is_null()));
        assert_eq!(control.ulTotalProcessingTime, 0);
        assert_eq!(control.num_err, 0);
        assert_eq!(control.num_inp, 0);
        assert_eq!(control.szTitle[0], 0);
        assert!(control.inchi_file.iter().all(|stream| {
            stream.type_ == INCHI_IOS_TYPE_STRING as i32
                && stream.s == INCHI_IOS_STRING::default()
                && stream.f.is_null()
        }));
        assert!(!control.strbuf_container.pStr.is_null());
        assert_eq!(
            control.strbuf_container.nAllocatedLength,
            INCHI_STRBUF_INITIAL_SIZE as i32
        );
        assert_eq!(control.strbuf_container.nPtr, INCHI_STRBUF_SIZE_INCREMENT as i32);
        assert_eq!(control.strbuf_container.nUsedLength, 0);
        assert!(
            heap.slice(control.strbuf_container.pStr.as_const())
                .unwrap()
                .iter()
                .all(|byte| *byte == 0)
        );

        inchi_free(&mut heap, control.strbuf_container.pStr).unwrap();
        inchi_free(&mut heap, control_pointer).unwrap();
        assert_eq!(heap.live_source_allocation_count(), 0);
    }

    #[test]
    fn source_port__inchi_dll_a__inchigen_destroy__line_1315() {
        let mut null_heap = SourceHeap::default();
        assert_eq!(INCHIGEN_Destroy(&mut null_heap, SourceMutPointer::null()), Ok(()));

        let mut created_heap = SourceHeap::default();
        created_heap.trace_source_allocations();
        let created = INCHIGEN_Create(&mut created_heap).unwrap();
        assert_eq!(created_heap.live_source_allocation_count(), 2);
        assert_eq!(INCHIGEN_Destroy(&mut created_heap, created), Ok(()));
        assert_eq!(created_heap.live_source_allocation_count(), 0);

        let mut heap = SourceHeap::default();
        let strbuf = heap.allocate_model_storage(vec![1_i8, 0]).unwrap();
        let stream_buffers: [SourceMutPointer<i8>; 3] =
            std::array::from_fn(|index| heap.allocate_model_storage(vec![index as i8 + 2, 0]).unwrap());
        let retained_path = heap.allocate_model_storage(vec![b'p' as i8, 0]).unwrap();
        let mut control = INCHIGEN_CONTROL::default();
        control.strbuf_container = INCHI_IOS_STRING {
            pStr: strbuf,
            nUsedLength: 1,
            nAllocatedLength: 2,
            nPtr: 4,
        };
        for (stream, buffer) in control.inchi_file.iter_mut().zip(stream_buffers) {
            stream.s = INCHI_IOS_STRING {
                pStr: buffer,
                nUsedLength: 1,
                nAllocatedLength: 2,
                nPtr: 4,
            };
        }
        control.InpParms.path[0] = retained_path.as_const();
        let control_pointer = heap.allocate_model_storage(vec![control]).unwrap();
        let handle: SourceMutPointer<SourceVoid> = control_pointer.cast();
        assert_eq!(INCHIGEN_Destroy(&mut heap, handle), Ok(()));
        assert_freed(&heap, strbuf);
        for buffer in stream_buffers {
            assert_freed(&heap, buffer);
        }
        assert_freed(&heap, control_pointer);
        assert_live(&heap, retained_path);
    }

    #[test]
    fn source_port__inchi_dll_a__stdinchigen_destroy__line_1308() {
        let mut null_heap = SourceHeap::default();
        assert_eq!(STDINCHIGEN_Destroy(&mut null_heap, SourceMutPointer::null()), Ok(()));
        assert_eq!(null_heap.live_source_allocation_count(), 0);

        let mut heap = SourceHeap::default();
        heap.trace_source_allocations();
        let handle = STDINCHIGEN_Create(&mut heap).unwrap();
        assert!(!handle.is_null());
        assert_eq!(heap.live_source_allocation_count(), 2);
        assert_eq!(heap.live_source_allocations_of::<INCHIGEN_CONTROL>(), 1);
        assert_eq!(heap.live_source_allocations_of::<i8>(), 1);
        assert_eq!(STDINCHIGEN_Destroy(&mut heap, handle), Ok(()));
        assert_eq!(heap.live_source_allocation_count(), 0);
    }

    #[test]
    fn source_port__inchi_dll_a__stdinchigen_create__line_122() {
        for successful_allocations in [0, 1] {
            let mut heap = SourceHeap::default();
            heap.fail_after_allocations(successful_allocations);
            assert_eq!(STDINCHIGEN_Create(&mut heap), Ok(SourceMutPointer::null()));
            assert_eq!(heap.source_allocation_calls(), successful_allocations + 1);
            assert_eq!(heap.live_source_allocation_count(), 0);
        }

        let mut heap = SourceHeap::default();
        heap.trace_source_allocations();
        let handle = STDINCHIGEN_Create(&mut heap).unwrap();
        assert!(!handle.is_null());
        assert_eq!(heap.source_allocation_calls(), 2);
        assert_eq!(heap.live_source_allocation_count(), 2);
        let control_pointer = handle.cast::<INCHIGEN_CONTROL>();
        let buffer = heap.slice(control_pointer.as_const()).unwrap()[0].strbuf_container.pStr;
        assert!(!buffer.is_null());
        inchi_free(&mut heap, buffer).unwrap();
        inchi_free(&mut heap, control_pointer).unwrap();
        assert_eq!(heap.live_source_allocation_count(), 0);
    }

    #[test]
    fn source_port__inchi_dll_a__stdinchigen_reset__line_1125() {
        let mut heap = SourceHeap::default();
        let inchi = heap.allocate_model_storage(vec![1_i8, 0]).unwrap();
        let log = heap.allocate_model_storage(vec![2_i8, 0]).unwrap();
        let message = heap.allocate_model_storage(vec![3_i8, 0]).unwrap();
        let mut results = inchi_Output {
            szInChI: inchi,
            szLog: log,
            szMessage: message,
            ..inchi_Output::default()
        };
        let mut data = INCHIGEN_DATA {
            pStrErrStruct: [9; 256],
            num_components: [7, 8],
            ..INCHIGEN_DATA::default()
        };
        assert_eq!(
            STDINCHIGEN_Reset(&mut heap, SourceMutPointer::null(), Some(&mut data), Some(&mut results),),
            Ok(())
        );
        assert_eq!(results, inchi_Output::default());
        assert_eq!(data, INCHIGEN_DATA::default());
        assert_freed(&heap, inchi);
        assert_freed(&heap, log);
        assert_freed(&heap, message);

        let mut untouched_data = INCHIGEN_DATA {
            num_components: [4, 5],
            ..INCHIGEN_DATA::default()
        };
        assert_eq!(
            STDINCHIGEN_Reset(&mut heap, SourceMutPointer::null(), Some(&mut untouched_data), None,),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(untouched_data.num_components, [4, 5]);
    }
}
