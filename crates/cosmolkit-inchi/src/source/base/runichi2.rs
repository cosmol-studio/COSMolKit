use crate::source::base::ichi_io::inchi_ios_eprint;
use crate::source::base::ichierr::AddErrorMessage;
use crate::source::base::ichimake::GetInpStructErrorType;
use crate::source::base::ichinorm::{OAD_Edit_Underivatize, Ring2Chain};
use crate::source::base::mol_fmt2::MolfileSaveCopy;
use crate::source_types::{
    _IS_EOF, _IS_ERROR, _IS_FATAL, _IS_OKAY, _IS_SKIP, _IS_WARNING, CANON_GLOBALS, INCHI_CLOCK,
    INCHI_IOSTREAM, INCHI_OUT_SDFILE_ONLY, INPUT_PARMS, LOG_MASK_ALL, LOG_MASK_ERR, LOG_MASK_FATAL,
    LOG_MASK_WARN, ORIG_ATOM_DATA, STRUCT_DATA, SourceFormatArgument, SourceHeap, SourceHeapError,
    SourceMutPointer, SourceVaList, clock_t,
};

fn source_c_bytes(
    heap: &SourceHeap,
    pointer: SourceMutPointer<i8>,
) -> Result<Vec<u8>, SourceHeapError> {
    if pointer.is_null() {
        return Ok(Vec::new());
    }
    let bytes = heap.slice(pointer.as_const())?;
    let length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    Ok(bytes[..length].iter().map(|byte| *byte as u8).collect())
}

fn source_array_c_bytes(bytes: &[i8]) -> Result<Vec<u8>, SourceHeapError> {
    let length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    Ok(bytes[..length].iter().map(|byte| *byte as u8).collect())
}

fn sdf_label_value(heap: &SourceHeap, input: &INPUT_PARMS) -> Result<Vec<u8>, SourceHeapError> {
    let label = source_c_bytes(heap, input.pSdfLabel)?;
    let value = source_c_bytes(heap, input.pSdfValue)?;
    let mut output = Vec::new();
    if !label.is_empty() {
        output.push(b' ');
        output.extend_from_slice(&label);
        if !value.is_empty() {
            output.push(b'=');
            output.extend_from_slice(&value);
        } else {
            output.push(b' ');
            output.extend_from_slice(b"is missing");
        }
    } else if !value.is_empty() {
        output.extend_from_slice(&value);
    }
    Ok(output)
}

fn eprint_bytes(
    heap: &mut SourceHeap,
    stream: Option<&mut INCHI_IOSTREAM>,
    bytes: &[u8],
) -> Result<i32, SourceHeapError> {
    let format = heap.allocate_model_storage(vec![b'%' as i8, b's' as i8, 0])?;
    let mut string = bytes.iter().map(|byte| *byte as i8).collect::<Vec<_>>();
    string.push(0);
    let string = heap.allocate_model_storage(string)?;
    let result = inchi_ios_eprint(
        heap,
        stream,
        format.as_const(),
        &SourceVaList {
            arguments: vec![SourceFormatArgument::Bytes(string.as_const())],
            position: 0,
        },
    );
    heap.free(string)?;
    heap.free(format)?;
    result
}

fn append_i64(output: &mut Vec<u8>, value: i64) {
    output.extend_from_slice(value.to_string().as_bytes());
}

fn append_i32(output: &mut Vec<u8>, value: i32) {
    output.extend_from_slice(value.to_string().as_bytes());
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn TreatErrorsInReadTheStructure(
    heap: &mut SourceHeap,
    structure: &mut STRUCT_DATA,
    input_parameters: &INPUT_PARMS,
    log_mask: i32,
    mut input_file: Option<&mut INCHI_IOSTREAM>,
    mut log_file: Option<&mut INCHI_IOSTREAM>,
    _output_file: Option<&mut INCHI_IOSTREAM>,
    mut problem_file: Option<&mut INCHI_IOSTREAM>,
    original_input: &ORIG_ATOM_DATA,
    input_number: &mut i64,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:716 TreatErrorsInReadTheStructure
    // INCHI✔️❌: int TreatErrorsInReadTheStructure( STRUCT_DATA      *sd,
    // INCHI✔️❌:                                    INPUT_PARMS      *ip,
    // INCHI✔️❌:                                    int              nLogMask,
    // INCHI✔️❌:                                    INCHI_IOSTREAM   *inp_file,
    // INCHI✔️❌:                                    INCHI_IOSTREAM   *log_file,
    // INCHI✔️❌:                                    INCHI_IOSTREAM   *out_file,
    // INCHI✔️❌:                                    INCHI_IOSTREAM   *prb_file,
    // INCHI✔️❌:                                    ORIG_ATOM_DATA   *orig_inp_data,
    // INCHI✔️❌:                                    long             *num_inp )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int nRet = _IS_OKAY;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (10 < sd->nStructReadError && sd->nStructReadError < 20)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  End of file */
    // INCHI✔️❌:         if (sd->pStrErrStruct[0])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint( log_file, "%s inp structure #%ld: End of file.%s%s%s%s    \n", sd->pStrErrStruct, *num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         inchi_ios_eprint( log_file, "End of file detected after structure #%ld.   \n", *num_inp - 1 );
    // INCHI✔️❌:
    // INCHI✔️❌:         nRet = _IS_EOF;
    // INCHI✔️❌:         goto exit_function; /*  end of file */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*(*num_inp) ++;*/
    // INCHI✔️❌:
    // INCHI✔️❌:     if (*num_inp < ip->first_struct_number)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  Skip the structure */
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( !defined(TARGET_API_LIB) && !defined(TARGET_EXE_STANDALONE) )
    // INCHI✔️❌: /* #ifndef TARGET_API_LIB */
    // INCHI✔️❌:         if (log_file->f != stderr)
    // INCHI✔️❌:             inchi_fprintf( stderr, "\rSkipping structure #%ld.%s%s%s%s...", *num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         nRet = sd->nErrorType = _IS_SKIP;
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     sd->nErrorType = GetInpStructErrorType( ip, sd->nStructReadError,
    // INCHI✔️❌:                                             sd->pStrErrStruct,
    // INCHI✔️❌:                                             orig_inp_data->num_inp_atoms );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (sd->nErrorType == _IS_FATAL)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  Fatal error */
    // INCHI✔️❌:         if (nLogMask & LOG_MASK_FATAL)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint( log_file, "Fatal Error %d (aborted; %s) inp structure #%ld.%s%s%s%s\n",
    // INCHI✔️❌:                      sd->nStructReadError, sd->pStrErrStruct, *num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( bRELEASE_VERSION == 1 || EXTR_FLAGS == 0 )
    // INCHI✔️❌:         /* djb-rwth: fixing oss-fuzz issue #27902 */
    // INCHI✔️❌:         if (prb_file && prb_file->f && 0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd && !ip->bSaveAllGoodStructsAsProblem)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             MolfileSaveCopy( inp_file, sd->fPtrStart, sd->fPtrEnd, prb_file->f, *num_inp );
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         /* goto exit_function; */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (sd->nErrorType == _IS_ERROR)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  Non-fatal errors: do not produce INChI */
    // INCHI✔️❌:
    // INCHI✔️❌:         /*  70 => too many atoms */
    // INCHI✔️❌:         if (nLogMask & LOG_MASK_ERR)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint( log_file, "Error %d (no %s; %s) inp structure #%ld.%s%s%s%s\n",
    // INCHI✔️❌:                  sd->nStructReadError, ( ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY ) ? "Molfile" : INCHI_NAME,
    // INCHI✔️❌:                  sd->pStrErrStruct, *num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( bRELEASE_VERSION == 1 || EXTR_FLAGS == 0 )
    // INCHI✔️❌:         if (prb_file && prb_file->f && 0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd && !ip->bSaveAllGoodStructsAsProblem)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             MolfileSaveCopy( inp_file, sd->fPtrStart, sd->fPtrEnd, prb_file->f, *num_inp ); /* djb-rwth: addressing coverity ID #499477 -- return values handled properly */
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (sd->nErrorType == _IS_WARNING)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  Warnings: try to produce INChI */
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nLogMask & LOG_MASK_WARN)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint( log_file, "Warning: (%s) inp structure #%ld.%s%s%s%s\n",
    // INCHI✔️❌:                sd->pStrErrStruct, *num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔️❌:     if (( ip->bINChIOutputOptions & INCHI_OUT_WINCHI_WINDOW ) &&
    // INCHI✔️❌:         ( ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (sd->nErrorType != _IS_OKAY && sd->nErrorType != _IS_WARNING)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sd->nErrorType =
    // INCHI✔️❌:                 ProcessStructError( out_file, log_file, sd->pStrErrStruct, sd->nErrorType, *num_inp, ip );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     if (nRet <= _IS_OKAY && sd->nErrorType > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nRet = sd->nErrorType;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return nRet;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: TreatErrorsInReadTheStructure

    let mut result = _IS_OKAY as i32;
    let error_text = source_array_c_bytes(&structure.pStrErrStruct)?;
    let suffix = sdf_label_value(heap, input_parameters)?;
    if structure.nStructReadError > 10 && structure.nStructReadError < 20 {
        if !error_text.is_empty() {
            let mut message = error_text.clone();
            message.extend_from_slice(b" inp structure #");
            append_i64(&mut message, *input_number);
            message.extend_from_slice(b": End of file.");
            message.extend_from_slice(&suffix);
            message.extend_from_slice(b"    \n");
            eprint_bytes(heap, log_file.as_deref_mut(), &message)?;
        }
        let mut message = b"End of file detected after structure #".to_vec();
        append_i64(&mut message, input_number.wrapping_sub(1));
        message.extend_from_slice(b".   \n");
        eprint_bytes(heap, log_file.as_deref_mut(), &message)?;
        result = _IS_EOF;
    } else if *input_number < input_parameters.first_struct_number {
        structure.nErrorType = _IS_SKIP;
        result = _IS_SKIP;
    } else {
        structure.nErrorType = GetInpStructErrorType(
            Some(input_parameters),
            structure.nStructReadError,
            Some(&structure.pStrErrStruct),
            original_input.num_inp_atoms,
        )?;
        if structure.nErrorType == _IS_FATAL as i32 {
            if log_mask & LOG_MASK_FATAL as i32 != 0 {
                let mut message = b"Fatal Error ".to_vec();
                append_i32(&mut message, structure.nStructReadError);
                message.extend_from_slice(b" (aborted; ");
                message.extend_from_slice(&error_text);
                message.extend_from_slice(b") inp structure #");
                append_i64(&mut message, *input_number);
                message.push(b'.');
                message.extend_from_slice(&suffix);
                message.push(b'\n');
                eprint_bytes(heap, log_file.as_deref_mut(), &message)?;
            }
            if let Some(problem) = problem_file.as_deref_mut()
                && input_file.is_some()
                && !problem.f.is_null()
                && structure.fPtrStart >= 0
                && structure.fPtrStart < structure.fPtrEnd
                && input_parameters.bSaveAllGoodStructsAsProblem == 0
            {
                MolfileSaveCopy(
                    heap,
                    input_file.as_deref_mut(),
                    structure.fPtrStart,
                    structure.fPtrEnd,
                    problem.f,
                    *input_number,
                )?;
            }
        }
        if structure.nErrorType == _IS_ERROR as i32 {
            if log_mask & LOG_MASK_ERR as i32 != 0 {
                let mut message = b"Error ".to_vec();
                append_i32(&mut message, structure.nStructReadError);
                message.extend_from_slice(b" (no ");
                message.extend_from_slice(
                    if input_parameters.bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY as i32 != 0 {
                        b"Molfile"
                    } else {
                        b"InChI"
                    },
                );
                message.extend_from_slice(b"; ");
                message.extend_from_slice(&error_text);
                message.extend_from_slice(b") inp structure #");
                append_i64(&mut message, *input_number);
                message.push(b'.');
                message.extend_from_slice(&suffix);
                message.push(b'\n');
                eprint_bytes(heap, log_file.as_deref_mut(), &message)?;
            }
            if let Some(problem) = problem_file.as_deref_mut()
                && input_file.is_some()
                && !problem.f.is_null()
                && structure.fPtrStart >= 0
                && structure.fPtrStart < structure.fPtrEnd
                && input_parameters.bSaveAllGoodStructsAsProblem == 0
            {
                MolfileSaveCopy(
                    heap,
                    input_file.as_deref_mut(),
                    structure.fPtrStart,
                    structure.fPtrEnd,
                    problem.f,
                    *input_number,
                )?;
            }
        }
        if structure.nErrorType == _IS_WARNING as i32 && log_mask & LOG_MASK_WARN as i32 != 0 {
            let mut message = b"Warning: (".to_vec();
            message.extend_from_slice(&error_text);
            message.extend_from_slice(b") inp structure #");
            append_i64(&mut message, *input_number);
            message.push(b'.');
            message.extend_from_slice(&suffix);
            message.push(b'\n');
            eprint_bytes(heap, log_file.as_deref_mut(), &message)?;
        }
    }

    if result <= _IS_OKAY as i32 && structure.nErrorType > 0 {
        result = structure.nErrorType;
    }
    Ok(result)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn DoOneStructureEarlyPreprocessing(
    heap: &mut SourceHeap,
    clock: SourceMutPointer<INCHI_CLOCK>,
    canon_globals: &mut CANON_GLOBALS,
    input_number: i64,
    structure: &mut STRUCT_DATA,
    input_parameters: &INPUT_PARMS,
    mut input_file: Option<&mut INCHI_IOSTREAM>,
    mut log_file: Option<&mut INCHI_IOSTREAM>,
    mut output_file: Option<&mut INCHI_IOSTREAM>,
    mut problem_file: Option<&mut INCHI_IOSTREAM>,
    original_input: &mut ORIG_ATOM_DATA,
    prepared_input: &ORIG_ATOM_DATA,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:493 DoOneStructureEarlyPreprocessing
    // INCHI✔️❌: int DoOneStructureEarlyPreprocessing( INCHI_CLOCK *ic,
    // INCHI✔️❌:                                       CANON_GLOBALS *pCG,
    // INCHI✔️❌:                                       long num_inp,
    // INCHI✔️❌:                                       STRUCT_DATA *sd,
    // INCHI✔️❌:                                       INPUT_PARMS *ip,
    // INCHI✔️❌:                                       INCHI_IOSTREAM *inp_file,
    // INCHI✔️❌:                                       INCHI_IOSTREAM *log_file,
    // INCHI✔️❌:                                       INCHI_IOSTREAM *out_file,
    // INCHI✔️❌:                                       INCHI_IOSTREAM *prb_file,
    // INCHI✔️❌:                                       ORIG_ATOM_DATA *orig_inp_data,
    // INCHI✔️❌:                                       ORIG_ATOM_DATA *prep_inp_data )
    // INCHI✔️❌: {
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )
    // INCHI✔️❌:     int ret1 = 0, ret2 = 0; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( REMOVE_ION_PAIRS_ORIG_STRU == 1 )
    // INCHI✔️❌:     fix_odd_things( orig_inp_data->num_inp_atoms, orig_inp_data->at, 0, ip->bFixNonUniformDraw );
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( UNDERIVATIZE == 1 )
    // INCHI✔️❌:     if (ip->bUnderivatize)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (0 > ( ret2 = OAD_Edit_Underivatize( ic, pCG, orig_inp_data, ( ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY ), ip->bUnderivatize & 2, ip->pSdfValue ) ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             long num_inp2 = num_inp;
    // INCHI✔️❌:             AddErrorMessage( sd->pStrErrStruct, "Underivatization error" );
    // INCHI✔️❌:             sd->nStructReadError = 99;
    // INCHI✔️❌:             sd->nErrorType = _IS_ERROR;
    // INCHI✔️❌:             TreatErrorsInReadTheStructure( sd, ip, LOG_MASK_ALL, inp_file, log_file, out_file, prb_file,
    // INCHI✔️❌:                                         prep_inp_data, &num_inp2 );
    // INCHI✔️❌:             return _IS_ERROR; /* output only if derivatives found */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (0 < ret2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!ip->bNoWarnings)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 WarningMessage( sd->pStrErrStruct, "Input structure underivatized" );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif /* UNDERIVATIZE == 1 */
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( RING2CHAIN == 1 )
    // INCHI✔️❌:     if (ip->bRing2Chain && 0 > ( ret1 = Ring2Chain( ic, pCG, orig_inp_data ) ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         long num_inp2 = num_inp;
    // INCHI✔️❌:         AddErrorMessage( sd->pStrErrStruct, "Ring to chain error" );
    // INCHI✔️❌:         sd->nStructReadError = 99;
    // INCHI✔️❌:         sd->nErrorType = _IS_ERROR;
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         TreatErrorsInReadTheStructure( sd, ip, LOG_MASK_ALL,
    // INCHI✔️❌:                                        inp_file, log_file,
    // INCHI✔️❌:                                        out_file, prb_file,
    // INCHI✔️❌:                                        prep_inp_data, &num_inp2 );
    // INCHI✔️❌:         return _IS_ERROR; /* output only if derivatives found */
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif /* RING2CHAIN == 1 */
    // INCHI✔️❌: #if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )  /***** post v.1 feature *****/
    // INCHI✔️❌:     if (ip->bIgnoreUnchanged && !ret1 && !ret2)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return _IS_SKIP; /* output only if derivatives or ring/chain found */
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif /* RING2CHAIN == 1 || UNDERIVATIZE == 1 */
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: DoOneStructureEarlyPreprocessing
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: DoOneStructureEarlyPreprocessing
    // INCHI✔️❌: #define UNDERIVATIZE               1 /* split to possible underivatized fragments */
    // INCHI✔️❌: #define RING2CHAIN                 1 /* open rings R-C(-OH)-O-R => R-C(=O) OH-R   */
    // INCHI✔️❌: #define REMOVE_ION_PAIRS_ORIG_STRU 0
    // INCHI✔️❌: #define WarningMessage AddErrorMessage
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: DoOneStructureEarlyPreprocessing

    let mut ring_result = 0_i32;
    let mut underivatize_result = 0_i32;

    if input_parameters.bUnderivatize != 0 {
        underivatize_result = OAD_Edit_Underivatize(
            heap,
            clock,
            canon_globals,
            original_input,
            input_parameters.bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY as i32,
            input_parameters.bUnderivatize & 2,
            input_parameters.pSdfValue,
            clock_result,
        )?;
        if underivatize_result < 0 {
            let mut copied_input_number = input_number;
            let message = (*b"Underivatization error\0").map(|byte| byte as i8);
            AddErrorMessage(Some(&mut structure.pStrErrStruct), Some(&message))?;
            structure.nStructReadError = 99;
            structure.nErrorType = _IS_ERROR as i32;
            TreatErrorsInReadTheStructure(
                heap,
                structure,
                input_parameters,
                LOG_MASK_ALL as i32,
                input_file.as_deref_mut(),
                log_file.as_deref_mut(),
                output_file.as_deref_mut(),
                problem_file.as_deref_mut(),
                prepared_input,
                &mut copied_input_number,
            )?;
            return Ok(_IS_ERROR as i32);
        }
        if underivatize_result > 0 && input_parameters.bNoWarnings == 0 {
            let message = (*b"Input structure underivatized\0").map(|byte| byte as i8);
            AddErrorMessage(Some(&mut structure.pStrErrStruct), Some(&message))?;
        }
    }

    if input_parameters.bRing2Chain != 0 {
        ring_result = Ring2Chain(heap, clock, canon_globals, original_input, clock_result)?;
        if ring_result < 0 {
            let mut copied_input_number = input_number;
            let message = (*b"Ring to chain error\0").map(|byte| byte as i8);
            AddErrorMessage(Some(&mut structure.pStrErrStruct), Some(&message))?;
            structure.nStructReadError = 99;
            structure.nErrorType = _IS_ERROR as i32;
            TreatErrorsInReadTheStructure(
                heap,
                structure,
                input_parameters,
                LOG_MASK_ALL as i32,
                input_file.as_deref_mut(),
                log_file.as_deref_mut(),
                output_file.as_deref_mut(),
                problem_file.as_deref_mut(),
                prepared_input,
                &mut copied_input_number,
            )?;
            return Ok(_IS_ERROR as i32);
        }
    }

    if input_parameters.bIgnoreUnchanged != 0 && ring_result == 0 && underivatize_result == 0 {
        return Ok(_IS_SKIP);
    }
    Ok(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{
        BOND_DOUBLE, BOND_SINGLE, INCHI_IOS_STRING, INCHI_IOS_TYPE_FILE, SourceFile, inp_ATOM,
    };

    const EL_NUMBER_C: u8 = 6;
    const EL_NUMBER_O: u8 = 8;

    fn c_string(heap: &mut SourceHeap, value: &[u8]) -> SourceMutPointer<i8> {
        let mut bytes = value.iter().map(|byte| *byte as i8).collect::<Vec<_>>();
        bytes.push(0);
        heap.allocate(bytes).unwrap()
    }

    fn log_stream() -> INCHI_IOSTREAM {
        INCHI_IOSTREAM {
            type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
            s: INCHI_IOS_STRING::default(),
            ..INCHI_IOSTREAM::default()
        }
    }

    fn log_bytes(heap: &SourceHeap, stream: &INCHI_IOSTREAM) -> Vec<u8> {
        if stream.s.pStr.is_null() {
            return Vec::new();
        }
        heap.slice(stream.s.pStr.as_const()).unwrap()[..stream.s.nUsedLength as usize]
            .iter()
            .map(|byte| *byte as u8)
            .collect()
    }

    fn early_original(heap: &mut SourceHeap, atoms: Vec<inp_ATOM>) -> ORIG_ATOM_DATA {
        let num_inp_bonds = atoms
            .iter()
            .map(|atom| i32::from(atom.valence))
            .sum::<i32>()
            / 2;
        ORIG_ATOM_DATA {
            num_inp_atoms: atoms.len() as i32,
            num_inp_bonds,
            at: heap.allocate_model_storage(atoms).unwrap(),
            ..ORIG_ATOM_DATA::default()
        }
    }

    fn early_carbon() -> inp_ATOM {
        let mut atom = inp_ATOM {
            el_number: EL_NUMBER_C,
            num_H: 4,
            orig_at_number: 1,
            ..inp_ATOM::default()
        };
        atom.elname[0] = b'C' as i8;
        atom
    }

    fn connect(atoms: &mut [inp_ATOM], first: usize, second: usize) {
        let first_order = atoms[first].valence as usize;
        atoms[first].neighbor[first_order] = second as u16;
        atoms[first].bond_type[first_order] = BOND_SINGLE as u8;
        atoms[first].valence += 1;
        atoms[first].chem_bonds_valence += 1;
        let second_order = atoms[second].valence as usize;
        atoms[second].neighbor[second_order] = first as u16;
        atoms[second].bond_type[second_order] = BOND_SINGLE as u8;
        atoms[second].valence += 1;
        atoms[second].chem_bonds_valence += 1;
    }

    fn early_ring() -> Vec<inp_ATOM> {
        let mut atoms = vec![inp_ATOM::default(); 8];
        for (index, atom) in atoms.iter_mut().enumerate() {
            atom.orig_at_number = index as u16 + 1;
        }
        for index in [0, 3, 4, 6, 7] {
            atoms[index].el_number = EL_NUMBER_C;
            atoms[index].elname[0] = b'C' as i8;
        }
        for index in [1, 2, 5] {
            atoms[index].el_number = EL_NUMBER_O;
            atoms[index].elname[0] = b'O' as i8;
        }
        atoms[1].num_H = 1;
        atoms[3].num_H = 1;
        atoms[4].num_H = 3;
        atoms[5].num_H = 1;
        atoms[6].num_H = 2;
        atoms[7].num_H = 2;
        connect(&mut atoms, 0, 1);
        connect(&mut atoms, 0, 2);
        connect(&mut atoms, 2, 6);
        connect(&mut atoms, 6, 7);
        connect(&mut atoms, 7, 3);
        connect(&mut atoms, 3, 0);
        connect(&mut atoms, 0, 4);
        connect(&mut atoms, 3, 5);
        atoms
    }

    #[test]
    fn source_port__runichi__doonestructureearlypreprocessing__line_493() {
        let run = |parameters: INPUT_PARMS, atoms: Vec<inp_ATOM>| {
            let mut heap = SourceHeap::default();
            let mut original = early_original(&mut heap, atoms);
            let prepared = original.clone();
            let clock = heap
                .allocate_model_storage(vec![INCHI_CLOCK::default()])
                .unwrap();
            let mut structure = STRUCT_DATA::default();
            let result = DoOneStructureEarlyPreprocessing(
                &mut heap,
                clock,
                &mut CANON_GLOBALS::default(),
                7,
                &mut structure,
                &parameters,
                None,
                None,
                None,
                None,
                &mut original,
                &prepared,
                0,
            );
            (heap, original, structure, result)
        };

        let (_, unchanged, structure, result) = run(INPUT_PARMS::default(), vec![early_carbon()]);
        assert_eq!(result, Ok(0));
        assert_eq!(unchanged.num_inp_atoms, 1);
        assert_eq!(structure, STRUCT_DATA::default());

        let (_, _, structure, result) = run(
            INPUT_PARMS {
                bIgnoreUnchanged: 1,
                ..INPUT_PARMS::default()
            },
            vec![early_carbon()],
        );
        assert_eq!(result, Ok(_IS_SKIP));
        assert_eq!(structure, STRUCT_DATA::default());

        let (_, unchanged, structure, result) = run(
            INPUT_PARMS {
                bUnderivatize: 1,
                bRing2Chain: 1,
                ..INPUT_PARMS::default()
            },
            vec![early_carbon()],
        );
        assert_eq!(result, Ok(0));
        assert_eq!(unchanged.num_inp_atoms, 1);
        assert_eq!(structure.pStrErrStruct[0], 0);

        let (heap, converted, structure, result) = run(
            INPUT_PARMS {
                bRing2Chain: 1,
                bIgnoreUnchanged: 1,
                ..INPUT_PARMS::default()
            },
            early_ring(),
        );
        assert_eq!(result, Ok(0));
        assert_eq!(structure.pStrErrStruct[0], 0);
        let atoms = &heap.slice(converted.at.as_const()).unwrap()[..8];
        let carbon = atoms
            .iter()
            .position(|atom| atom.orig_at_number == 1)
            .unwrap();
        let hydroxyl = atoms
            .iter()
            .position(|atom| atom.orig_at_number == 2)
            .unwrap();
        let opened = atoms
            .iter()
            .position(|atom| atom.orig_at_number == 3)
            .unwrap();
        let carbon_to_hydroxyl = (0..atoms[carbon].valence as usize)
            .find(|order| atoms[carbon].neighbor[*order] as usize == hydroxyl)
            .unwrap();
        assert_eq!(
            atoms[carbon].bond_type[carbon_to_hydroxyl],
            BOND_DOUBLE as u8
        );
        assert!(
            !(0..atoms[carbon].valence as usize)
                .any(|order| atoms[carbon].neighbor[order] as usize == opened)
        );
    }

    #[test]
    fn source_port__runichi2__treaterrorsinreadthestructure__line_716() {
        let mut heap = SourceHeap::default();
        let label = c_string(&mut heap, b"ID");
        let value = c_string(&mut heap, b"42");
        let input_parameters = INPUT_PARMS {
            pSdfLabel: label,
            pSdfValue: value,
            first_struct_number: 1,
            ..INPUT_PARMS::default()
        };
        let input_file_pointer = heap
            .allocate(vec![SourceFile {
                bytes: b"name\nbody\n".to_vec(),
                ..SourceFile::default()
            }])
            .unwrap();
        let mut input_file = INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_FILE as i32,
            f: input_file_pointer,
            ..INCHI_IOSTREAM::default()
        };
        let mut log = log_stream();
        let mut eof = STRUCT_DATA {
            nStructReadError: 11,
            ..STRUCT_DATA::default()
        };
        eof.pStrErrStruct[..4].copy_from_slice(&[b's' as i8, b't' as i8, b'o' as i8, b'p' as i8]);
        let mut number = 3;
        assert_eq!(
            TreatErrorsInReadTheStructure(
                &mut heap,
                &mut eof,
                &input_parameters,
                7,
                Some(&mut input_file),
                Some(&mut log),
                None,
                None,
                &ORIG_ATOM_DATA::default(),
                &mut number,
            ),
            Ok(_IS_EOF)
        );
        assert_eq!(
            log_bytes(&heap, &log),
            b"stop inp structure #3: End of file. ID=42    \nEnd of file detected after structure #2.   \n"
        );

        let mut skipped = STRUCT_DATA::default();
        number = 0;
        assert_eq!(
            TreatErrorsInReadTheStructure(
                &mut heap,
                &mut skipped,
                &input_parameters,
                7,
                Some(&mut input_file),
                None,
                None,
                None,
                &ORIG_ATOM_DATA::default(),
                &mut number,
            ),
            Ok(_IS_SKIP)
        );
        assert_eq!(skipped.nErrorType, _IS_SKIP);

        let problem_pointer = heap.allocate(vec![SourceFile::default()]).unwrap();
        let mut problem = INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_FILE as i32,
            f: problem_pointer,
            ..INCHI_IOSTREAM::default()
        };
        let mut fatal = STRUCT_DATA {
            nStructReadError: 5,
            fPtrStart: 0,
            fPtrEnd: 5,
            ..STRUCT_DATA::default()
        };
        fatal.pStrErrStruct[..3].copy_from_slice(&[b'b' as i8, b'a' as i8, b'd' as i8]);
        let mut fatal_log = log_stream();
        number = 8;
        assert_eq!(
            TreatErrorsInReadTheStructure(
                &mut heap,
                &mut fatal,
                &input_parameters,
                LOG_MASK_FATAL as i32,
                Some(&mut input_file),
                Some(&mut fatal_log),
                None,
                Some(&mut problem),
                &ORIG_ATOM_DATA {
                    num_inp_atoms: 1,
                    ..ORIG_ATOM_DATA::default()
                },
                &mut number,
            ),
            Ok(_IS_FATAL as i32)
        );
        assert_eq!(fatal.nErrorType, _IS_FATAL as i32);
        assert_eq!(
            log_bytes(&heap, &fatal_log),
            b"Fatal Error 5 (aborted; bad) inp structure #8. ID=42\n"
        );
        assert_eq!(
            heap.slice(problem_pointer.as_const()).unwrap()[0].bytes,
            b"#8/name\n".to_vec()
        );

        let mut error = STRUCT_DATA {
            nStructReadError: 9,
            ..STRUCT_DATA::default()
        };
        error.pStrErrStruct[..4].copy_from_slice(&[b'n' as i8, b'o' as i8, b'p' as i8, b'e' as i8]);
        let mut error_log = log_stream();
        let sdf_parameters = INPUT_PARMS {
            bINChIOutputOptions: INCHI_OUT_SDFILE_ONLY as i32,
            ..input_parameters.clone()
        };
        assert_eq!(
            TreatErrorsInReadTheStructure(
                &mut heap,
                &mut error,
                &sdf_parameters,
                LOG_MASK_ERR as i32,
                Some(&mut input_file),
                Some(&mut error_log),
                None,
                None,
                &ORIG_ATOM_DATA::default(),
                &mut number,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(
            log_bytes(&heap, &error_log),
            b"Error 9 (no Molfile; nope) inp structure #8. ID=42\n"
        );

        let mut warning = STRUCT_DATA::default();
        warning.pStrErrStruct[0] = b'w' as i8;
        let mut warning_log = log_stream();
        assert_eq!(
            TreatErrorsInReadTheStructure(
                &mut heap,
                &mut warning,
                &input_parameters,
                LOG_MASK_WARN as i32,
                Some(&mut input_file),
                Some(&mut warning_log),
                None,
                None,
                &ORIG_ATOM_DATA {
                    num_inp_atoms: 1,
                    ..ORIG_ATOM_DATA::default()
                },
                &mut number,
            ),
            Ok(_IS_WARNING as i32)
        );
        assert_eq!(
            log_bytes(&heap, &warning_log),
            b"Warning: (w) inp structure #8. ID=42\n"
        );

        let mut okay = STRUCT_DATA::default();
        assert_eq!(
            TreatErrorsInReadTheStructure(
                &mut heap,
                &mut okay,
                &input_parameters,
                0,
                Some(&mut input_file),
                None,
                None,
                None,
                &ORIG_ATOM_DATA {
                    num_inp_atoms: 1,
                    ..ORIG_ATOM_DATA::default()
                },
                &mut number,
            ),
            Ok(_IS_OKAY as i32)
        );
    }
}
