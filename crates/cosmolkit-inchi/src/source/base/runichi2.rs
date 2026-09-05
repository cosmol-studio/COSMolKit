use crate::source::base::ichi_io::{
    inchi_ios_close, inchi_ios_eprint, inchi_ios_init, inchi_strbuf_close, inchi_strbuf_create_copy, inchi_strbuf_init,
    inchi_strbuf_printf,
};
use crate::source::base::ichicano::{InchiTimeElapsed, InchiTimeGet};
use crate::source::base::ichierr::AddErrorMessage;
use crate::source::base::ichimake::GetInpStructErrorType;
use crate::source::base::ichinorm::{OAD_Edit_Underivatize, Ring2Chain};
use crate::source::base::ichiprt2::inchi_strtol;
use crate::source::base::ichiread::{
    extract_all_backbone_bonds_from_inchi_string, extract_stereo_info_from_inchi_string,
};
use crate::source::base::ichisort::iisort;
use crate::source::base::ichister::{cross_prod3, dot_prod3};
use crate::source::base::mol_fmt2::MolfileSaveCopy;
use crate::source::base::mol_fmt4::{
    IntArray_Alloc, IntArray_Append, IntArray_AppendIfAbsent, IntArray_DebugPrint, IntArray_Free,
};
use crate::source::base::mol2atom::{CreateInpAtomData, FreeOrigAtData};
use crate::source::base::runichi3::{
    Inp_Atom_GetBondType, OAD_CollectBackboneBonds, OAD_CollectReachableAtoms, OAD_PolymerUnit_DebugTrace,
    OAD_PolymerUnit_FindEndsAndCaps, OAD_ValidatePolymerAndPseudoElementData, OrigAtData_Duplicate,
};
use crate::source::base::runichi4::FreeAllINChIArrays;
use crate::source::base::strutil::ExtractConnectedComponent;
use crate::source::base::util::{inchi_calloc, inchi_free, inchi_malloc, is_in_the_ilist};
use crate::source_types::local_runichi2::{DiylFrag, ModSCenterInfo};
use crate::source_types::{
    _IS_EOF, _IS_ERROR, _IS_FATAL, _IS_OKAY, _IS_SKIP, _IS_WARNING, BOND_TYPE_SINGLE, CANON_GLOBALS, CT_ATOMCOUNT_ERR,
    CT_UNKNOWN_ERR, INCHI_CLOCK, INCHI_IOS_STRING, INCHI_IOS_TYPE_STRING, INCHI_IOSTREAM, INCHI_NUM,
    INCHI_OUT_SDFILE_ONLY, INCHI_STRBUF_INITIAL_SIZE, INCHI_STRBUF_SIZE_INCREMENT, INP_ATOM_DATA, INPUT_PARMS,
    INT_ARRAY, LOG_MASK_ALL, LOG_MASK_ERR, LOG_MASK_FATAL, LOG_MASK_WARN, MAX_ATOMS, OAD_StructureEdits,
    ORIG_ATOM_DATA, PINChI_Aux2, PINChI2, POLYMERS_MODERN, POSEContext, STEREO_SNGL_DOWN, STEREO_SNGL_UP, STR_ERR_LEN,
    STRUCT_DATA, SourceConstPointer, SourceFormatArgument, SourceHeap, SourceHeapError, SourceMutPointer, SourceVaList,
    clock_t, inchiTime, tagINCHIStereoParity0D_INCHI_PARITY_EVEN as INCHI_PARITY_EVEN,
    tagINCHIStereoParity0D_INCHI_PARITY_ODD as INCHI_PARITY_ODD,
};

pub(crate) fn source_c_bytes(heap: &SourceHeap, pointer: SourceMutPointer<i8>) -> Result<Vec<u8>, SourceHeapError> {
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

pub(crate) fn source_array_c_bytes(bytes: &[i8]) -> Result<Vec<u8>, SourceHeapError> {
    let length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    Ok(bytes[..length].iter().map(|byte| *byte as u8).collect())
}

pub(crate) fn sdf_label_value(heap: &SourceHeap, input: &INPUT_PARMS) -> Result<Vec<u8>, SourceHeapError> {
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

pub(crate) fn eprint_bytes(
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
pub(crate) fn GetOneComponent(
    heap: &mut SourceHeap,
    clock: &mut INCHI_CLOCK,
    structure: &mut STRUCT_DATA,
    input_parameters: &INPUT_PARMS,
    mut log_file: Option<&mut INCHI_IOSTREAM>,
    _output_file: Option<&mut INCHI_IOSTREAM>,
    current_input: &mut INP_ATOM_DATA,
    original_input: &ORIG_ATOM_DATA,
    component_index: i32,
    input_number: i64,
    start_clock_result: clock_t,
    end_clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:345 GetOneComponent
    // INCHI✔️❌: complete source frame follows verbatim; SourceHeap borrowing requires temporary atom copies.
    /*
    int GetOneComponent( INCHI_CLOCK        *ic,
                         STRUCT_DATA        *sd,
                         INPUT_PARMS        *ip,
                         INCHI_IOSTREAM     *log_file,
                         INCHI_IOSTREAM     *out_file,
                         INP_ATOM_DATA      *inp_cur_data,
                         ORIG_ATOM_DATA     *orig_inp_data,
                         int                i,
                         long               num_inp )
    {
        inchiTime ulTStart;

        InchiTimeGet( &ulTStart );

        CreateInpAtomData( inp_cur_data, orig_inp_data->nCurAtLen[i], 0 );

        inp_cur_data->num_at = ExtractConnectedComponent( orig_inp_data->at, orig_inp_data->num_inp_atoms, i + 1, inp_cur_data->at );

        sd->ulStructTime += InchiTimeElapsed( ic, &ulTStart );


        /*  Error processing */
        if (inp_cur_data->num_at <= 0 || orig_inp_data->nCurAtLen[i] != inp_cur_data->num_at)
        {
            /*  Log error message */
            AddErrorMessage( sd->pStrErrStruct, "Cannot extract Component" );
            inchi_ios_eprint( log_file, "%s #%d structure #%ld.%s%s%s%s\n", sd->pStrErrStruct, i + 1, num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
            sd->nErrorCode = inp_cur_data->num_at < 0 ? inp_cur_data->num_at : ( orig_inp_data->nCurAtLen[i] != inp_cur_data->num_at ) ? CT_ATOMCOUNT_ERR : CT_UNKNOWN_ERR;
            /* num_err ++; */
            sd->nErrorType = _IS_ERROR;
    #ifdef TARGET_LIB_FOR_WINCHI
            if (( ip->bINChIOutputOptions & INCHI_OUT_WINCHI_WINDOW ) &&
                ( ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT ))
            {
                sd->nErrorType = ProcessStructError( out_file,
                                                     log_file,
                                                     sd->pStrErrStruct,
                                                     sd->nErrorType,
                                                     num_inp,
                                                     ip );
            }
    #endif
        }

        return sd->nErrorType;
    }
    */
    // END INCHI C FUNCTION: GetOneComponent
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: GetOneComponent
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; TARGET_LIB_FOR_WINCHI inactive.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: GetOneComponent

    let component_offset = usize::try_from(component_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let expected_atoms = i32::from(
        *heap
            .slice(original_input.nCurAtLen.as_const())?
            .get(component_offset)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let mut start = inchiTime::default();
    InchiTimeGet(&mut start, start_clock_result);

    let _ = CreateInpAtomData(heap, current_input, expected_atoms, 0)?;
    let source_atoms = heap
        .slice(original_input.at.as_const())?
        .get(..usize::try_from(original_input.num_inp_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    let mut component_atoms = if current_input.at.is_null() {
        Vec::new()
    } else {
        heap.slice(current_input.at.as_const())?.to_vec()
    };
    current_input.num_at = ExtractConnectedComponent(
        heap,
        &source_atoms,
        original_input.num_inp_atoms,
        component_index.wrapping_add(1),
        &mut component_atoms,
    )?;
    if !current_input.at.is_null() {
        heap.slice_mut(current_input.at)?.clone_from_slice(&component_atoms);
    }

    let elapsed = InchiTimeElapsed(clock, Some(&start), end_clock_result);
    structure.ulStructTime = structure.ulStructTime.wrapping_add(elapsed as u64);

    if current_input.num_at <= 0 || expected_atoms != current_input.num_at {
        let message = (*b"Cannot extract Component\0").map(|byte| byte as i8);
        AddErrorMessage(Some(&mut structure.pStrErrStruct), Some(&message))?;

        let mut output = source_array_c_bytes(&structure.pStrErrStruct)?;
        output.extend_from_slice(b" #");
        append_i32(&mut output, component_index.wrapping_add(1));
        output.extend_from_slice(b" structure #");
        append_i64(&mut output, input_number);
        output.push(b'.');
        output.extend_from_slice(&sdf_label_value(heap, input_parameters)?);
        output.push(b'\n');
        eprint_bytes(heap, log_file.as_deref_mut(), &output)?;

        structure.nErrorCode = if current_input.num_at < 0 {
            current_input.num_at
        } else if expected_atoms != current_input.num_at {
            CT_ATOMCOUNT_ERR
        } else {
            CT_UNKNOWN_ERR
        };
        structure.nErrorType = _IS_ERROR as i32;
    }

    Ok(structure.nErrorType)
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

#[allow(non_snake_case)]
pub(crate) const fn bIsSameBond(a1: i32, a2: i32, b1: i32, b2: i32) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1293 bIsSameBond
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: bIsSameBond
    // INCHI✔️✔️: int bIsSameBond(int a1, int a2, int b1, int b2)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     if (a1 == b1 && a2 == b2) return 1;
    // INCHI✔️✔️:     if (a1 == b2 && a2 == b1) return -1;
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: bIsSameBond
    // END INCHI C FUNCTION: bIsSameBond
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: bIsSameBond
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; int is signed 32-bit.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: bIsSameBond

    if a1 == b1 && a2 == b2 {
        1
    } else if a1 == b2 && a2 == b1 {
        -1
    } else {
        0
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn GetFrameShiftInfoFrom105PlusInChI(
    heap: &mut SourceHeap,
    sinchi: SourceConstPointer<i8>,
    frame_shift_info: &mut [i32],
    max_crossing: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1305 GetFrameShiftInfoFrom105PlusInChI
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
static int GetFrameShiftInfoFrom105PlusInChI(char *sinchi,
    int *frame_shift_info,
    int max_crossing)
{
    int k, c = 0, j, aindex = 0, iunit = 0;
    const char *p, *q;

    p = strstr(sinchi, "/z");  /* must always be there */

                               /* each frame_shift_info triple(iunit, iunit_a1, iunit_a2) contains,
                               for each eligible frame-shiftable unit,
                               iunit - unit_no
                               iunit_a1, iunit_a2 - atom numbers for the senior bkbond
                               note that iunit_a1 is more senior then iunit_a2
                               */

                               /* eligible unit has Z-layer pattern
                               "range-of-numbers(number1,number2,nimbers...)"
                               >=2 bkbonds in CRU; relink may be necessary to shift frame or swap bkbond atoms?
                               OPTIONALLY DO NOT DO SWAP NOW?
                               senior bkbond and right atoms order is (number1,number2)
                               "range-of-numbers(number1.number2)"
                               1-bkbond CRU; relink may still be necessary to swap bkbond atoms so that
                               more senior atom is connected to lesser-numbered Zz
                               OPTIONALLY DO NOT DO SWAP NOW?
                               senior atoms order in bkbond is (number1,number2)
                               "range-of-numbers(number)"  1-atom CRU
                               relink is not applicable, skip it
                               */

    while (p)
    {
        int num[2] = { -1,-1 };
        p = strstr(p + 2, "(");
        if (!p)
        {
            break;
        }
        p++;
        q = p;
        j = 0;
        while ((k = (int)inchi_strtol(p, &q, 10)) && j < 2)
        {
            num[j] = k;
            j++;
            c = UCINT *q;
            if (j==1 && c == '-') /* do not consider pattern "(cap-end, cap-end)" */
            {
                goto find_next_unit;
            }
            else if (c != ')')
            {
                p = q + 1;
            }
            else
            {
                goto find_next_unit;
            }
        }
        if (j < 2)
        {
            goto find_next_unit;
        }
        frame_shift_info[3 * aindex] = iunit;
        frame_shift_info[3 * aindex + 1] = num[0];
        frame_shift_info[3 * aindex + 2] = num[1];
        aindex++;
        if (aindex >= max_crossing)
        {
            break;
        }

    find_next_unit:
        p = strstr(p, ";");
        iunit++;
    }

    return aindex;
}
    */
    // END INCHI C FUNCTION: GetFrameShiftInfoFrom105PlusInChI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetFrameShiftInfoFrom105PlusInChI
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; UCINT casts the inspected char through unsigned char.
    // INCHI✔️✔️: strstr is the libc implementation and inchi_strtol dispatches to LP64 libc strtol through its completed source port.
    // INCHI✔️✔️: The left operand of && invokes inchi_strtol once more after j reaches 2; its value and end pointer are intentionally discarded.
    // INCHI✔️✔️: SourceHeap pointer traversal has the same linear scans and no additional allocation or cloning.
    // END INCHI ACTIVE MACRO CONFIGURATION: GetFrameShiftInfoFrom105PlusInChI

    fn find(
        heap: &SourceHeap,
        haystack: SourceConstPointer<i8>,
        needle: &[u8],
    ) -> Result<Option<SourceConstPointer<i8>>, SourceHeapError> {
        let bytes = heap.slice(haystack)?;
        let nul = bytes
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        let Some(offset) = bytes[..nul]
            .windows(needle.len())
            .position(|window| {
                window
                    .iter()
                    .zip(needle)
                    .all(|(actual, expected)| *actual as u8 == *expected)
            })
        else {
            return Ok(None);
        };
        Ok(Some(haystack.offset(
            i64::try_from(offset).map_err(|_| SourceHeapError::PointerOffsetOverflow)?,
        )?))
    }

    let mut p = find(heap, sinchi, b"/z")?;
    let mut aindex = 0_i32;
    let mut iunit = 0_i32;
    while let Some(current) = p {
        let search_start = current.offset(2)?;
        let Some(open_parenthesis) = find(heap, search_start, b"(")? else {
            break;
        };
        let mut current_number = open_parenthesis.offset(1)?;
        let mut q = current_number;
        let mut j = 0_i32;
        let mut numbers = [-1_i32; 2];
        let mut find_next_unit = false;

        loop {
            let k = inchi_strtol(heap, current_number, Some(&mut q), 10)? as i32;
            if k == 0 || j >= 2 {
                break;
            }
            numbers[usize::try_from(j).map_err(|_| SourceHeapError::SourceIntegerOverflow)?] = k;
            j = j.wrapping_add(1);
            let c = *heap
                .slice(q)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)? as u8;
            if j == 1 && c == b'-' {
                find_next_unit = true;
                break;
            } else if c != b')' {
                current_number = q.offset(1)?;
            } else {
                find_next_unit = true;
                break;
            }
        }

        if !find_next_unit && j >= 2 {
            let base = usize::try_from(aindex)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                .checked_mul(3)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            *frame_shift_info
                .get_mut(base)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = iunit;
            *frame_shift_info
                .get_mut(base + 1)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = numbers[0];
            *frame_shift_info
                .get_mut(base + 2)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = numbers[1];
            aindex = aindex.wrapping_add(1);
            if aindex >= max_crossing {
                break;
            }
        }

        p = find(heap, current_number, b";")?;
        iunit = iunit.wrapping_add(1);
    }

    Ok(aindex)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_PrepareFrameShiftEdits(
    heap: &mut SourceHeap,
    original_atom_data: &ORIG_ATOM_DATA,
    sinchi: SourceConstPointer<i8>,
    saux: SourceConstPointer<i8>,
    edits: &OAD_StructureEdits,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2523 OAD_Polymer_PrepareFrameShiftEdits
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int  OAD_Polymer_PrepareFrameShiftEdits( ORIG_ATOM_DATA *orig_at_data,
                                         char *sinchi,
                                         char *saux,
                                         OAD_StructureEdits *ed)
{
    int ret = _IS_OKAY;
    int *orig = NULL, *frame_shift_info = NULL;
    int n_frame_shifts, j;
    ModSCenterInfo *scinfo = NULL;		/* 4 elements; [0]th for old_end1, [1] old_end2, [2] end1, [3] end2	*/
    
    OAD_Polymer *p = orig_at_data->polymer;
    int nu = orig_at_data->polymer->n;
    int nat = orig_at_data->num_inp_atoms;
    
    /* Extract cano_nums-->orig_nums mapping for InChI AuxInfo Main Layer */
    orig = (int *)inchi_calloc((long long)nat + 1, sizeof(int)); /* djb-rwth: cast operator added */
    if (!orig)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    ret = extract_orig_nums_from_auxinfo_string(saux, orig);
    if (ret != _IS_OKAY && ret != _IS_WARNING)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }

    scinfo = (ModSCenterInfo *)inchi_calloc(4, sizeof(scinfo[0]));
    if (!scinfo)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    

    /* Parse InChI and extract, for each 'bistar' CRU, the senior bkbond (to frame-shift brackets to its ends) */
    frame_shift_info = (int *)inchi_calloc(3 * ((long long)nu + 1), sizeof(int)); /* djb-rwth: cast operator added */
    if (!frame_shift_info)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    n_frame_shifts = GetFrameShiftInfoFrom105PlusInChI(sinchi, frame_shift_info, nu);
    /* translate atom numbers to orig numbers */
    for (j = 0; j < n_frame_shifts; j++)
    {
        frame_shift_info[3 * j + 1] = orig[frame_shift_info[3 * j + 1]];
        frame_shift_info[3 * j + 2] = orig[frame_shift_info[3 * j + 2]];
    }

    /* Collect OAD edits */
    for (j = 0; j < n_frame_shifts; j++)
    {
        OAD_PolymerUnit *u = NULL;
        int k, iu = -1; /*int iu = frame_shift_info[3 * j];*/
        int end1, cap1, cap1_is_star, end2, cap2, cap2_is_star, old_end1, old_end2, err, fail = 0;

        end1 = frame_shift_info[3 * j + 1];
        end2 = frame_shift_info[3 * j + 2];
        
        /* Find the unit to edit (== that unit whose alist contains the new end atoms) */
        for (k = 0; k < p->n; k++)
        {
            int ak, present=0;
            
            if (NULL == p->units[k]->blist || p->units[k]->nb < 2 )
            {
                /* No crossing bonds in the unit */
                continue;
            }
            /* Find the unit to edit (== that unit whose backbone contains the new end atoms)
            for (bk = 0; bk < p->units[k]->nbkbonds; bk++ )
            {
                if ( bIsSameBond(end1, end2, p->units[k]->bkbonds[bk][0], p->units[k]->bkbonds[bk][1] ) )
                {
                    iu = k;
                    break;
                }
            }*/
            for (ak = 0; ak < p->units[k]->na; ak++ )
            {
                if (p->units[k]->alist[ak] == end1 || p->units[k]->alist[ak] == end2)
                {
                    present++;
                }
                if (present==2)
                {
                    iu = k;
                    break;
                }
            }
        }
        if (iu < 0)
        {
            /* Unit to edit unexpectedly not found, that's an error */
            ret = _IS_ERROR;
            goto exit_function;
        }

        u = p->units[iu];

        OAD_PolymerUnit_FindEndsAndCaps(u, orig_at_data,
                                        &old_end1, &cap1, &cap1_is_star,
                                        &old_end2, &cap2, &cap2_is_star,
                                        &err, NULL);

        if (!err && cap1_is_star && cap2_is_star && end1 && end2 && cap1 && cap2)
        {
            /* find old CRU ends */
            if (cap1 == u->blist[0])		old_end1 = u->blist[1];
            else if (cap1 == u->blist[1])	old_end1 = u->blist[0];
            else if (cap1 == u->blist[2])	old_end1 = u->blist[3];
            else if (cap1 == u->blist[3])	old_end1 = u->blist[2];
            else /* something wrong */
                continue;
            if (cap2 == u->blist[0])		old_end2 = u->blist[1];
            else if (cap2 == u->blist[1])	old_end2 = u->blist[0];
            else if (cap2 == u->blist[2])	old_end2 = u->blist[3];
            else if (cap2 == u->blist[3])	old_end2 = u->blist[2];
            else /* something wrong */
                continue;

            if (!old_end1 || !old_end2 || old_end1 == old_end2)
            {
                continue;
            }
            if (bIsSameBond(old_end1, cap1, end1, cap1) && bIsSameBond(old_end2, cap2, end2, cap2))
            {
                continue; /* ignore swaps for now */
            }

            /*	If applicable, collect bonds to modify */
            
            /* Check if atoms involved in modifications are stereocenters (needs additional care) */
            ModSCenter_Init(&scinfo[0], orig_at_data->at, old_end1 - 1);
            ModSCenter_Init(&scinfo[1], orig_at_data->at, old_end2 - 1);
            ModSCenter_Init(&scinfo[2], orig_at_data->at, end1 - 1);
            ModSCenter_Init(&scinfo[3], orig_at_data->at, end2 - 1);

            /* djb-rwth: removing redundant code */
            if (!bIsSameBond(old_end1, cap1, end1, cap1))
            {
                /* Modify bond: (old_end1-cap1) --> (end1-cap1) */
                fail = 0;
                fail += IntArray_Append(ed->mod_bond, old_end1);
                fail += IntArray_Append(ed->mod_bond, cap1);
                fail += IntArray_Append(ed->mod_bond, end1);
                fail += IntArray_Append(ed->mod_bond, cap1);
                if (fail)
                {
                    ret = _IS_ERROR;
                    goto exit_function;
                }
                ModSCenter_DelFrom(&scinfo[0], cap1 - 1);
                ModSCenter_AddTo(&scinfo[2], cap1-1 );
            }
            if (!bIsSameBond(old_end2, cap2, end2, cap2))
            {
                /* Modify bond: (old_end2-cap2) --> (end2-cap2) */
                fail = 0;
                fail += IntArray_Append(ed->mod_bond, old_end2);
                fail += IntArray_Append(ed->mod_bond, cap2);
                fail += IntArray_Append(ed->mod_bond, end2);
                fail += IntArray_Append(ed->mod_bond, cap2);
                if (fail)
                {
                    ret = _IS_ERROR;
                    goto exit_function;
                }
                ModSCenter_DelFrom(&scinfo[1], cap2 - 1);
                ModSCenter_AddTo(&scinfo[3], cap2 - 1);
            }
            /* Modify bond: (end1-end2) --> (old_end1-old_end2) */
            fail = 0;
            fail += IntArray_Append(ed->mod_bond, end1);
            fail += IntArray_Append(ed->mod_bond, end2);
            fail += IntArray_Append(ed->mod_bond, old_end1);
            fail += IntArray_Append(ed->mod_bond, old_end2);
            if (fail)
            {
                ret = _IS_ERROR;
                goto exit_function;
            }
            ModSCenter_DelFrom(&scinfo[2], end2 - 1);
            ModSCenter_DelFrom(&scinfo[3], end1 - 1);
            ModSCenter_AddTo(&scinfo[0], old_end2 - 1);
            ModSCenter_AddTo(&scinfo[1], old_end1 - 1);

        }

        /* djb-rwth: n_flip and ModSCenter_IsChanged function completely redundant? -- discussion required */
        if (orig_at_data->num_dimensions)
        {
            /* Check if we must flip stereocenter configuration */
            /* (ignore errrors signaled by returning -1)		*/
            int n_flip = 0;
            if (0 < ModSCenter_IsChanged(&scinfo[0], orig_at_data->at))
            {
                n_flip++;
            }
            if (0 < ModSCenter_IsChanged(&scinfo[1], orig_at_data->at))
            {
                n_flip++;
            }
            if (0 < ModSCenter_IsChanged(&scinfo[2], orig_at_data->at))
            {
                n_flip++;
            }
            if (0 < ModSCenter_IsChanged(&scinfo[3], orig_at_data->at))
            {
                n_flip++;
            }
            n_flip = 1;
        }
    }

exit_function:
    if (orig)
    {
        inchi_free(orig);
    }
    if (frame_shift_info)
    {
        inchi_free(frame_shift_info);
    }
    if (scinfo)
    {
        inchi_free(scinfo);
    }

    return ret;
}
    */
    // END INCHI C FUNCTION: OAD_Polymer_PrepareFrameShiftEdits
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OAD_Polymer_PrepareFrameShiftEdits
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; sizeof(int)=4 and sizeof(ModSCenterInfo)=172.
    // INCHI✔️❌: inchi_calloc/free resolve to the active libc macros; every direct function callee is an independently completed source port.
    // INCHI✔️❌: Allocation return mapping, all aggregate append attempts, unit selection, one-based atom translation, stereo calls, and cleanup order are preserved.
    // INCHI✔️❌: SourceHeap pointer-map access and atom/unit snapshots add material overhead beyond native pointer indexing.
    // END INCHI ACTIVE MACRO CONFIGURATION: OAD_Polymer_PrepareFrameShiftEdits

    fn is_source_allocation_error(error: SourceHeapError) -> bool {
        matches!(
            error,
            SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange
                | SourceHeapError::AllocationFailed
        )
    }

    fn append_modified_bond(
        heap: &mut SourceHeap,
        edits: &OAD_StructureEdits,
        value: i32,
    ) -> Result<i32, SourceHeapError> {
        if edits.mod_bond.is_null() {
            return IntArray_Append(heap, None, value);
        }
        heap.with_slice_mut_and_heap_mut(edits.mod_bond, |arrays, heap| {
            IntArray_Append(heap, arrays.first_mut(), value)
        })
    }

    let mut orig = SourceMutPointer::null();
    let mut frame_shift_info = SourceMutPointer::null();
    let mut stereo_info = SourceMutPointer::null();

    let operation = (|| -> Result<i32, SourceHeapError> {
        let polymer = heap
            .slice(original_atom_data.polymer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let unit_count = polymer.n;
        let atom_count = original_atom_data.num_inp_atoms;

        let orig_count = i64::from(atom_count).wrapping_add(1) as u64;
        orig = match inchi_calloc::<i32>(heap, orig_count, 4) {
            Ok(pointer) => pointer,
            Err(error) if is_source_allocation_error(error) => return Ok(_IS_ERROR as i32),
            Err(error) => return Err(error),
        };
        let result = heap.with_slice_mut_and_heap_mut(orig, |values, heap| {
            extract_orig_nums_from_auxinfo_string(heap, saux, values)
        })?;
        if result != _IS_OKAY as i32 && result != _IS_WARNING as i32 {
            return Ok(_IS_ERROR as i32);
        }

        stereo_info = match inchi_calloc::<ModSCenterInfo>(heap, 4, 172) {
            Ok(pointer) => pointer,
            Err(error) if is_source_allocation_error(error) => return Ok(_IS_ERROR as i32),
            Err(error) => return Err(error),
        };

        let frame_count = i64::from(unit_count)
            .wrapping_add(1)
            .wrapping_mul(3) as u64;
        frame_shift_info = match inchi_calloc::<i32>(heap, frame_count, 4) {
            Ok(pointer) => pointer,
            Err(error) if is_source_allocation_error(error) => return Ok(_IS_ERROR as i32),
            Err(error) => return Err(error),
        };
        let frame_shift_count =
            heap.with_slice_mut_and_heap_mut(frame_shift_info, |values, heap| {
                GetFrameShiftInfoFrom105PlusInChI(heap, sinchi, values, unit_count)
            })?;

        for frame_index in 0..frame_shift_count {
            let base = usize::try_from(frame_index)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                .checked_mul(3)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            for offset in [1_usize, 2] {
                let canonical = *heap
                    .slice(frame_shift_info.as_const())?
                    .get(base + offset)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let original = *heap
                    .slice(orig.as_const())?
                    .get(
                        usize::try_from(canonical)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                *heap
                    .slice_mut(frame_shift_info)?
                    .get_mut(base + offset)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = original;
            }
        }

        let atoms = if frame_shift_count > 0 {
            let count = usize::try_from(atom_count)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            heap.slice(original_atom_data.at.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec()
        } else {
            Vec::new()
        };

        for frame_index in 0..frame_shift_count {
            let base = usize::try_from(frame_index)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                .checked_mul(3)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            let end1 = *heap
                .slice(frame_shift_info.as_const())?
                .get(base + 1)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let end2 = *heap
                .slice(frame_shift_info.as_const())?
                .get(base + 2)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;

            let mut unit_index = -1_i32;
            for candidate_index in 0..polymer.n {
                let candidate_pointer = *heap
                    .slice(polymer.units.as_const())?
                    .get(
                        usize::try_from(candidate_index)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let candidate = heap
                    .slice(candidate_pointer.as_const())?
                    .first()
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if candidate.blist.is_null() || candidate.nb < 2 {
                    continue;
                }
                let mut present = 0_i32;
                for atom_index in 0..candidate.na {
                    let atom = *heap
                        .slice(candidate.alist.as_const())?
                        .get(
                            usize::try_from(atom_index)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if atom == end1 || atom == end2 {
                        present = present.wrapping_add(1);
                    }
                    if present == 2 {
                        unit_index = candidate_index;
                        break;
                    }
                }
            }
            if unit_index < 0 {
                return Ok(_IS_ERROR as i32);
            }

            let unit_pointer = *heap
                .slice(polymer.units.as_const())?
                .get(
                    usize::try_from(unit_index)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mut unit = heap
                .slice(unit_pointer.as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mut old_end1 = 0_i32;
            let mut cap1 = 0_i32;
            let mut cap1_is_star = 0_i32;
            let mut old_end2 = 0_i32;
            let mut cap2 = 0_i32;
            let mut cap2_is_star = 0_i32;
            let mut error = 0_i32;
            OAD_PolymerUnit_FindEndsAndCaps(
                heap,
                &mut unit,
                original_atom_data,
                &mut old_end1,
                &mut cap1,
                &mut cap1_is_star,
                &mut old_end2,
                &mut cap2,
                &mut cap2_is_star,
                &mut error,
                None,
            )?;
            heap.slice_mut(unit_pointer)?[0] = unit.clone();

            if error == 0
                && cap1_is_star != 0
                && cap2_is_star != 0
                && end1 != 0
                && end2 != 0
                && cap1 != 0
                && cap2 != 0
            {
                let crossing = heap.slice(unit.blist.as_const())?;
                old_end1 = if cap1
                    == *crossing
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                {
                    *crossing
                        .get(1)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                } else if cap1
                    == *crossing
                        .get(1)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                {
                    *crossing
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                } else if cap1
                    == *crossing
                        .get(2)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                {
                    *crossing
                        .get(3)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                } else if cap1
                    == *crossing
                        .get(3)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                {
                    *crossing
                        .get(2)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                } else {
                    continue;
                };
                old_end2 = if cap2
                    == *crossing
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                {
                    *crossing
                        .get(1)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                } else if cap2
                    == *crossing
                        .get(1)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                {
                    *crossing
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                } else if cap2
                    == *crossing
                        .get(2)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                {
                    *crossing
                        .get(3)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                } else if cap2
                    == *crossing
                        .get(3)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                {
                    *crossing
                        .get(2)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                } else {
                    continue;
                };

                if old_end1 == 0 || old_end2 == 0 || old_end1 == old_end2 {
                    continue;
                }
                if bIsSameBond(old_end1, cap1, end1, cap1) != 0
                    && bIsSameBond(old_end2, cap2, end2, cap2) != 0
                {
                    continue;
                }

                let mut centers = heap
                    .slice(stereo_info.as_const())?
                    .get(..4)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                ModSCenter_Init(&mut centers[0], &atoms, old_end1.wrapping_sub(1))?;
                ModSCenter_Init(&mut centers[1], &atoms, old_end2.wrapping_sub(1))?;
                ModSCenter_Init(&mut centers[2], &atoms, end1.wrapping_sub(1))?;
                ModSCenter_Init(&mut centers[3], &atoms, end2.wrapping_sub(1))?;
                heap.slice_mut(stereo_info)?
                    .get_mut(..4)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone_from_slice(&centers);

                if bIsSameBond(old_end1, cap1, end1, cap1) == 0 {
                    let mut fail = 0_i32;
                    for value in [old_end1, cap1, end1, cap1] {
                        fail = fail.wrapping_add(append_modified_bond(heap, edits, value)?);
                    }
                    if fail != 0 {
                        return Ok(_IS_ERROR as i32);
                    }
                    ModSCenter_DelFrom(&mut centers[0], cap1.wrapping_sub(1))?;
                    ModSCenter_AddTo(&mut centers[2], cap1.wrapping_sub(1))?;
                }
                if bIsSameBond(old_end2, cap2, end2, cap2) == 0 {
                    let mut fail = 0_i32;
                    for value in [old_end2, cap2, end2, cap2] {
                        fail = fail.wrapping_add(append_modified_bond(heap, edits, value)?);
                    }
                    if fail != 0 {
                        return Ok(_IS_ERROR as i32);
                    }
                    ModSCenter_DelFrom(&mut centers[1], cap2.wrapping_sub(1))?;
                    ModSCenter_AddTo(&mut centers[3], cap2.wrapping_sub(1))?;
                }

                let mut fail = 0_i32;
                for value in [end1, end2, old_end1, old_end2] {
                    fail = fail.wrapping_add(append_modified_bond(heap, edits, value)?);
                }
                if fail != 0 {
                    return Ok(_IS_ERROR as i32);
                }
                ModSCenter_DelFrom(&mut centers[2], end2.wrapping_sub(1))?;
                ModSCenter_DelFrom(&mut centers[3], end1.wrapping_sub(1))?;
                ModSCenter_AddTo(&mut centers[0], old_end2.wrapping_sub(1))?;
                ModSCenter_AddTo(&mut centers[1], old_end1.wrapping_sub(1))?;
                heap.slice_mut(stereo_info)?
                    .get_mut(..4)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone_from_slice(&centers);
            }

            if original_atom_data.num_dimensions != 0 {
                let mut flip_count = 0_i32;
                for center_index in 0..4 {
                    let mut center = heap
                        .slice(stereo_info.as_const())?
                        .get(center_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if ModSCenter_IsChanged(heap, &mut center, &atoms)? > 0 {
                        flip_count = flip_count.wrapping_add(1);
                    }
                    heap.slice_mut(stereo_info)?[center_index] = center;
                }
                flip_count = 1;
                let _ = flip_count;
            }
        }

        Ok(result)
    })();

    let mut cleanup_error = None;
    if !orig.is_null() {
        if let Err(error) = inchi_free(heap, orig) {
            cleanup_error = Some(error);
        }
    }
    if !frame_shift_info.is_null() {
        if let Err(error) = inchi_free(heap, frame_shift_info)
            && cleanup_error.is_none()
        {
            cleanup_error = Some(error);
        }
    }
    if !stereo_info.is_null() {
        if let Err(error) = inchi_free(heap, stereo_info)
            && cleanup_error.is_none()
        {
            cleanup_error = Some(error);
        }
    }

    if let Some(error) = cleanup_error {
        Err(error)
    } else {
        operation
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn ModSCenter_Init(
    stereo_center: &mut ModSCenterInfo,
    atoms: &[crate::source_types::inp_ATOM],
    atom_index: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2760 ModSCenter_Init
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
void ModSCenter_Init(ModSCenterInfo *scinfo, inp_ATOM *at, int iatom)
{
    int i;
    scinfo->num = iatom;
    scinfo->valence = at[iatom].valence;
    scinfo->n_stereo = NDefStereoBonds(at, iatom, 1); /* , bOnlyPointedEndMatters=1 */
    for (i = 0; i < scinfo->valence; i++)
    {
        scinfo->nbr[i] = at[iatom].neighbor[i];
        scinfo->new_nbr[i] = scinfo->nbr[i];
    }

    return;
}
    */
    // END INCHI C FUNCTION: ModSCenter_Init
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: ModSCenter_Init
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; S_CHAR valence and AT_NUMB neighbors promote to int on assignment.
    // INCHI✔️✔️: NDefStereoBonds is the independently completed source port; all array accesses are direct and allocation-free.
    // END INCHI ACTIVE MACRO CONFIGURATION: ModSCenter_Init

    stereo_center.num = atom_index;
    let atom = atoms
        .get(usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    stereo_center.valence = i32::from(atom.valence);
    stereo_center.n_stereo = NDefStereoBonds(atoms, atom_index, 1)?;
    for index in 0..stereo_center.valence {
        let index = usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let neighbor = i32::from(
            *atom
                .neighbor
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        *stereo_center
            .nbr
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = neighbor;
        *stereo_center
            .new_nbr
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = stereo_center.nbr[index];
    }
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn NDefStereoBonds(
    atoms: &[crate::source_types::inp_ATOM],
    atom_index: i32,
    only_pointed_end_matters: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2775 NDefStereoBonds
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
int NDefStereoBonds(inp_ATOM *at, int iatom, int bOnlyPointedEndMatters)
{
    int i, n_stereo = 0;
    int stereo_value, stereo_type;
    for (i = 0; i < at[iatom].valence; i++)
    {
        stereo_value = at[iatom].bond_stereo[i];
        if (bOnlyPointedEndMatters)
        {
            /* establish the stereo considering only the pointed end of stereo bond */
            stereo_type = stereo_value;
        }
        else
        {
            stereo_type = abs(stereo_value);
        }
        if (stereo_type == STEREO_SNGL_UP || stereo_type == STEREO_SNGL_DOWN)
        {
            n_stereo++;
        }
    }
    return n_stereo;
}
    */
    // END INCHI C FUNCTION: NDefStereoBonds
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: NDefStereoBonds
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; S_CHAR promotes to signed int before abs and comparisons.
    // INCHI✔️✔️: STEREO_SNGL_UP=1 and STEREO_SNGL_DOWN=6; every fixed-array access remains direct and allocation-free.
    // END INCHI ACTIVE MACRO CONFIGURATION: NDefStereoBonds

    let atom = atoms
        .get(usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut stereo_count = 0_i32;
    for index in 0..i32::from(atom.valence) {
        let stereo_value = i32::from(
            *atom
                .bond_stereo
                .get(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let stereo_type = if only_pointed_end_matters != 0 {
            stereo_value
        } else {
            stereo_value.abs()
        };
        if stereo_type == STEREO_SNGL_UP as i32 || stereo_type == STEREO_SNGL_DOWN as i32 {
            stereo_count = stereo_count.wrapping_add(1);
        }
    }
    Ok(stereo_count)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn ModSCenter_AddTo(
    stereo_center: &mut ModSCenterInfo,
    atom_to_add: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2803 ModSCenter_AddTo
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
void ModSCenter_AddTo(ModSCenterInfo *scinfo, int iadd)
{
    if (!is_in_the_ilist(scinfo->new_nbr, iadd, scinfo->valence))
    {
        scinfo->new_nbr[scinfo->valence] = iadd;
        scinfo->valence++;
    }
    return;
}
    */
    // END INCHI C FUNCTION: ModSCenter_AddTo
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: ModSCenter_AddTo
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; is_in_the_ilist is the independently completed source port.
    // INCHI✔️✔️: Search and append are linear/direct, allocation-free, and preserve duplicate and tail-slot mutation order.
    // END INCHI ACTIVE MACRO CONFIGURATION: ModSCenter_AddTo

    if is_in_the_ilist(
        Some(&stereo_center.new_nbr),
        atom_to_add,
        stereo_center.valence,
    )?
    .is_none()
    {
        let index = usize::try_from(stereo_center.valence)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        *stereo_center
            .new_nbr
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = atom_to_add;
        stereo_center.valence = stereo_center.valence.wrapping_add(1);
    }
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn ModSCenter_DelFrom(
    stereo_center: &mut ModSCenterInfo,
    atom_to_delete: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2815 ModSCenter_DelFrom
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
void ModSCenter_DelFrom(ModSCenterInfo *scinfo, int idel)
{
    int i, j;
    for (i = 0; i < scinfo->valence; i++)
    {
        if (scinfo->nbr[i]==idel )
        {
            for (j=i+1; j < scinfo->valence; j++)
            {
                scinfo->new_nbr[j-1] = scinfo->new_nbr[j];
            }
            scinfo->valence--;
            return;
        }
    }
    return;
}
    */
    // END INCHI C FUNCTION: ModSCenter_DelFrom
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: ModSCenter_DelFrom
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; signed int loop and decrement semantics are preserved.
    // INCHI✔️✔️: Search uses nbr while shifting uses new_nbr; direct indexed mutation is allocation-free and leaves the old tail unchanged.
    // END INCHI ACTIVE MACRO CONFIGURATION: ModSCenter_DelFrom

    for index in 0..stereo_center.valence {
        let index = usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let neighbor = *stereo_center
            .nbr
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if neighbor == atom_to_delete {
            let mut shifted_index = index
                .checked_add(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            let valence = usize::try_from(stereo_center.valence)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            while shifted_index < valence {
                let shifted = *stereo_center
                    .new_nbr
                    .get(shifted_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                *stereo_center
                    .new_nbr
                    .get_mut(shifted_index - 1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = shifted;
                shifted_index += 1;
            }
            stereo_center.valence = stereo_center.valence.wrapping_sub(1);
            return Ok(());
        }
    }
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn ModSCenter_IsChanged(
    heap: &mut SourceHeap,
    stereo_center: &mut ModSCenterInfo,
    atoms: &[crate::source_types::inp_ATOM],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2836 ModSCenter_IsChanged
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int ModSCenter_IsChanged(ModSCenterInfo *scinfo, inp_ATOM *at)
{
    int i, ns, base1=-1, base2=-1, new_base2=-1, n_changed=0;
    double a[3], b[3], new_b[3], z[3], new_z[3], zz; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    if (scinfo->n_stereo < 1)
    {
        return 0;
    }
    if (scinfo->valence != at[scinfo->num].valence )
    {
        return -1; /* something went wrong */
    }
    iisort(scinfo->nbr, scinfo->valence);
    iisort(scinfo->new_nbr, scinfo->valence);
    /* Find the kept stereo base atom */
    for (i = 0; i < at[scinfo->num].valence; i++)
    {
        if ( is_in_the_ilist(scinfo->nbr, scinfo->new_nbr[i], scinfo->valence) )
        {
            ns = NDefStereoBonds(at, scinfo->new_nbr[i], 0); /* bOnlyPointedEndMatters=0 */
            if (ns==0)
            {
                base1 = scinfo->new_nbr[i];
                break;
            }
        }
    }
    if (base1==-1)
    {
        return -1; /* something went wrong */
    }
    /* Find the newly appeared stereo base atom */
    for (i = 0; i < at[scinfo->num].valence; i++)
    {
        /*!!! TUT NADO NE TAK
         base2 tot, kogo net v new_nbr
         new_base2 - tot, kogo net v  nbr
        */
        if ( !is_in_the_ilist(scinfo->nbr, scinfo->new_nbr[i], scinfo->valence))
        {
            ns = NDefStereoBonds(at, scinfo->nbr[i], 0);
            if (ns == 0)
            {
                new_base2 = scinfo->new_nbr[i];
                base2 = scinfo->nbr[i];
                n_changed++;
            }
        }
    }
    if (n_changed > 1 || new_base2 == -1 || base2 == -1)
    {
        return -1; /* something went wrong */
    }
    a[0] = at[base1].x - at[scinfo->num].x; a[1] = at[base1].y - at[scinfo->num].y; a[2] = at[base1].z - at[scinfo->num].z;
    b[0] = at[base2].x - at[scinfo->num].x; b[1] = at[base2].y - at[scinfo->num].y; b[2] = at[base2].z - at[scinfo->num].z;
    new_b[0] = at[new_base2].x - at[scinfo->num].x; new_b[1] = at[new_base2].y - at[scinfo->num].y; new_b[2] = at[new_base2].z - at[scinfo->num].z;

    cross_prod3(a, b, z);
    cross_prod3(a, new_b, new_z);
    zz = dot_prod3(z, new_z); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    return -1;
}
    */
    // END INCHI C FUNCTION: ModSCenter_IsChanged
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: ModSCenter_IsChanged
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; all direct callees are independently completed source ports.
    // INCHI✔️❌: Both in-place iisort side effects, integer comparisons, atom-index selection, floating evaluation order, and unconditional final -1 are preserved.
    // INCHI✔️❌: SourceHeap cannot address embedded struct arrays, so each iisort call uses and frees one model-storage mirror, adding material overhead without a source allocation call.
    // END INCHI ACTIVE MACRO CONFIGURATION: ModSCenter_IsChanged

    if stereo_center.n_stereo < 1 {
        return Ok(0);
    }
    let center_index = usize::try_from(stereo_center.num)
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let center_atom = atoms
        .get(center_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if stereo_center.valence != i32::from(center_atom.valence) {
        return Ok(-1);
    }

    let neighbors = heap.allocate_model_storage(stereo_center.nbr.to_vec())?;
    let sort_neighbors = iisort(heap, neighbors, stereo_center.valence);
    match sort_neighbors {
        Ok(_) => stereo_center
            .nbr
            .copy_from_slice(heap.slice(neighbors.as_const())?),
        Err(error) => {
            heap.free(neighbors)?;
            return Err(error);
        }
    }
    heap.free(neighbors)?;

    let new_neighbors = heap.allocate_model_storage(stereo_center.new_nbr.to_vec())?;
    let sort_new_neighbors = iisort(heap, new_neighbors, stereo_center.valence);
    match sort_new_neighbors {
        Ok(_) => stereo_center
            .new_nbr
            .copy_from_slice(heap.slice(new_neighbors.as_const())?),
        Err(error) => {
            heap.free(new_neighbors)?;
            return Err(error);
        }
    }
    heap.free(new_neighbors)?;

    let mut base1 = -1_i32;
    for index in 0..i32::from(center_atom.valence) {
        let index = usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let candidate = *stereo_center
            .new_nbr
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if is_in_the_ilist(
            Some(&stereo_center.nbr),
            candidate,
            stereo_center.valence,
        )?
        .is_some()
            && NDefStereoBonds(atoms, candidate, 0)? == 0
        {
            base1 = candidate;
            break;
        }
    }
    if base1 == -1 {
        return Ok(-1);
    }

    let mut base2 = -1_i32;
    let mut new_base2 = -1_i32;
    let mut changed = 0_i32;
    for index in 0..i32::from(center_atom.valence) {
        let index = usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let candidate = *stereo_center
            .new_nbr
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if is_in_the_ilist(
            Some(&stereo_center.nbr),
            candidate,
            stereo_center.valence,
        )?
        .is_none()
        {
            let old_candidate = *stereo_center
                .nbr
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if NDefStereoBonds(atoms, old_candidate, 0)? == 0 {
                new_base2 = candidate;
                base2 = old_candidate;
                changed = changed.wrapping_add(1);
            }
        }
    }
    if changed > 1 || new_base2 == -1 || base2 == -1 {
        return Ok(-1);
    }

    let base1_atom = atoms
        .get(usize::try_from(base1).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let base2_atom = atoms
        .get(usize::try_from(base2).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let new_base2_atom = atoms
        .get(usize::try_from(new_base2).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let a = [
        base1_atom.x - center_atom.x,
        base1_atom.y - center_atom.y,
        base1_atom.z - center_atom.z,
    ];
    let b = [
        base2_atom.x - center_atom.x,
        base2_atom.y - center_atom.y,
        base2_atom.z - center_atom.z,
    ];
    let new_b = [
        new_base2_atom.x - center_atom.x,
        new_base2_atom.y - center_atom.y,
        new_base2_atom.z - center_atom.z,
    ];
    let mut z = [0.0_f64; 3];
    let mut new_z = [0.0_f64; 3];
    cross_prod3(&a, &b, &mut z);
    cross_prod3(&a, &new_b, &mut new_z);
    let _zz = dot_prod3(&z, &new_z);

    Ok(-1)
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

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn POSEContext_Init(
    heap: &mut SourceHeap,
    context: &mut POSEContext,
    structure_data: Option<&STRUCT_DATA>,
    input_parameters: Option<&INPUT_PARMS>,
    title: &[i8],
    inchi_components: Option<&[SourceMutPointer<PINChI2>; INCHI_NUM as usize]>,
    aux_components: Option<&[SourceMutPointer<PINChI_Aux2>; INCHI_NUM as usize]>,
    input_file: SourceMutPointer<INCHI_IOSTREAM>,
    _log_file: SourceMutPointer<INCHI_IOSTREAM>,
    _output_file: SourceMutPointer<INCHI_IOSTREAM>,
    _problem_file: SourceMutPointer<INCHI_IOSTREAM>,
    original_input: Option<&ORIG_ATOM_DATA>,
    prepared_input: Option<&[ORIG_ATOM_DATA]>,
    input_number: i64,
    string_buffer: Option<&INCHI_IOS_STRING>,
    save_option_bits: u8,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1503 POSEContext_Init
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int  POSEContext_Init(POSEContext *context,
                      STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle,
                      PINChI2 *pINChI2[INCHI_NUM], PINChI_Aux2 *pINChI_Aux2[INCHI_NUM],
                      INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file,
                      INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file,
                      ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data,
                      long num_inp, INCHI_IOS_STRING *strbuf, unsigned char save_opt_bits)
{
    char *sz = NULL;
    int ret = _IS_OKAY, res = 0, i;

    memset(context, 0, sizeof(*context)); /* djb-rwth: memset_s C11/Annex K variant? */

    if (!sd)
    {
        memset(&context->sd, 0, sizeof(context->sd)); /* djb-rwth: memset_s C11/Annex K variant? */
    }
    else
    {
        memcpy(&context->sd, sd, sizeof(context->sd));
    }

    if (!ip)
    {
        memset(&context->ip, 0, sizeof(context->ip)); /* djb-rwth: memset_s C11/Annex K variant? */
    }
    else
    {
        memcpy(&context->ip, ip, sizeof(context->ip));
        for (i = 0; i < MAX_NUM_PATHS; i++)
        {
            if (ip->path[i])
            {
                sz = (char*)inchi_malloc((strlen(ip->path[i]) + 1) * sizeof(sz[0]));
                if (!sz)
                {
                    ret = _IS_ERROR;
                    goto exit_function;
                }
                strcpy(sz, context->ip.path[i]);
                context->ip.path[i] = sz;
            }
        }
    }

    if (strlen(szTitle))
    {
        strcpy(context->szTitle, szTitle);
    }
    else
    {
        context->szTitle[0] = '\0';
    }

    /* pINChI2, pINChI_Aux2: We do not fill/allocate elements of these structures   */
    /* assuming that NULL's are there. If not just raise an error.                  */

    context->pINChI2[0] = context->pINChI2[1] = NULL;
    if (pINChI2 && (pINChI2[0] || pINChI2[1])) /* djb-rwth: condition corrected */
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    context->pINChI_Aux2[0] = context->pINChI_Aux2[1] = NULL; 
    if (pINChI_Aux2 && (pINChI_Aux2[0] || pINChI_Aux2[1])) /* djb-rwth: condition corrected */
    {
        ret = _IS_ERROR;
        goto exit_function;
    }

    context->out_file = context->inchi_file;
    context->log_file = context->inchi_file + 1;
    context->prb_file = context->inchi_file + 2;
    /* Initialize internal for this function output streams as string buffers */
    inchi_ios_init(context->out_file, INCHI_IOS_TYPE_STRING, NULL);
    inchi_ios_init(context->log_file, INCHI_IOS_TYPE_STRING, NULL);
    inchi_ios_init(context->prb_file, INCHI_IOS_TYPE_STRING, NULL);
    context->inp_file = NULL;
    if (inp_file)
    {
        context->inp_file = inp_file;
    }

    context->orig_inp_data = &context->OrigAtData;
    context->prep_inp_data = context->PrepAtData;

    if (orig_inp_data)
    {
        memset(context->orig_inp_data, 0, sizeof(*context->orig_inp_data)); /* djb-rwth: memset_s C11/Annex K variant? */
        res = OrigAtData_Duplicate(context->orig_inp_data, orig_inp_data);
        if (res)
        {
            ret = _IS_ERROR;
            goto exit_function;
        }
    }

    if (prep_inp_data)
    {
        memset(context->prep_inp_data, 0, 2 * sizeof(*context->prep_inp_data)); /* djb-rwth: memset_s C11/Annex K variant? */
        res = OrigAtData_Duplicate(context->prep_inp_data, prep_inp_data);
        if (res)
        {
            ret = _IS_ERROR;
            goto exit_function;
        }
    }

    /* num_inp, strbuf, save_opt_bits */
    context->num_inp = num_inp;
    context->save_opt_bits = save_opt_bits;
    context->strbuf = &context->temp_string_container;
    if (strbuf)
    {
        res = inchi_strbuf_create_copy(context->strbuf, strbuf);
    }
    else
    {
        res = inchi_strbuf_init(context->strbuf, INCHI_STRBUF_INITIAL_SIZE, INCHI_STRBUF_SIZE_INCREMENT);
    }
    if (res == -1)
    {
        ret = _IS_FATAL;
        goto exit_function;
    }

exit_function:

    return ret;
}
    */
    // END INCHI C FUNCTION: POSEContext_Init
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: POSEContext_Init
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; MAX_NUM_PATHS=4 and INCHI_NUM=2.
    // INCHI✔️❌: inchi_malloc, inchi_ios_init, OrigAtData_Duplicate, inchi_strbuf_create_copy,
    // INCHI✔️❌: and inchi_strbuf_init dispatch to their independently ported active definitions.
    // INCHI✔️❌: Model-storage mirrors preserve C self-pointer aliasing but add allocation-map overhead.
    // END INCHI ACTIVE MACRO CONFIGURATION: POSEContext_Init

    *context = POSEContext::default();
    if let Some(structure_data) = structure_data {
        context.sd = structure_data.clone();
    }
    if let Some(input_parameters) = input_parameters {
        context.ip = input_parameters.clone();
        for index in 0..context.ip.path.len() {
            let source_path = input_parameters.path[index];
            if source_path.is_null() {
                continue;
            }
            let source = heap.slice(source_path)?;
            let length = source
                .iter()
                .position(|byte| *byte == 0)
                .ok_or(SourceHeapError::MissingNulTerminator)?;
            let source = source[..=length].to_vec();
            let byte_count = u64::try_from(length)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                .checked_add(1)
                .ok_or(SourceHeapError::AllocationSizeOverflow)?;
            let copied = match inchi_malloc(heap, byte_count) {
                Ok(pointer) => pointer,
                Err(
                    SourceHeapError::AllocationFailed
                    | SourceHeapError::AllocationSizeOverflow
                    | SourceHeapError::AllocationElementCountOutOfRange,
                ) => return Ok(_IS_ERROR as i32),
                Err(error) => return Err(error),
            };
            heap.slice_mut(copied)?[..=length].copy_from_slice(&source);
            context.ip.path[index] = copied.as_const();
        }
    }

    let title_length = title
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    if title_length >= context.szTitle.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    if title_length != 0 {
        context.szTitle[..=title_length].copy_from_slice(&title[..=title_length]);
    } else {
        context.szTitle[0] = 0;
    }

    context.pINChI2 = [SourceMutPointer::null(); INCHI_NUM as usize];
    if inchi_components.is_some_and(|values| values.iter().any(|pointer| !pointer.is_null())) {
        return Ok(_IS_ERROR as i32);
    }
    context.pINChI_Aux2 = [SourceMutPointer::null(); INCHI_NUM as usize];
    if aux_components.is_some_and(|values| values.iter().any(|pointer| !pointer.is_null())) {
        return Ok(_IS_ERROR as i32);
    }

    let stream_storage = heap.allocate_model_storage(vec![INCHI_IOSTREAM::default(); 3])?;
    context.out_file = stream_storage;
    context.log_file = stream_storage.offset(1)?;
    context.prb_file = stream_storage.offset(2)?;
    for index in 0..3_i64 {
        let stream = stream_storage.offset(index)?;
        inchi_ios_init(
            heap.slice_mut(stream)?.first_mut(),
            INCHI_IOS_TYPE_STRING as i32,
            SourceMutPointer::null(),
        )?;
    }
    context.inchi_file.clone_from_slice(&heap.slice(stream_storage.as_const())?[..3]);
    context.inp_file = input_file;

    let original_storage = heap.allocate_model_storage(vec![ORIG_ATOM_DATA::default()])?;
    context.orig_inp_data = original_storage;
    if let Some(original_input) = original_input {
        let duplicate_result = heap.with_slice_mut_and_heap_mut(original_storage, |values, heap| {
            let destination = values
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            OrigAtData_Duplicate(heap, destination, original_input)
        })?;
        context.OrigAtData = heap.slice(original_storage.as_const())?[0].clone();
        if duplicate_result != 0 {
            return Ok(_IS_ERROR as i32);
        }
    }

    let prepared_storage = heap.allocate_model_storage(vec![ORIG_ATOM_DATA::default(); 2])?;
    context.prep_inp_data = prepared_storage;
    if let Some(prepared_input) = prepared_input {
        let source = prepared_input
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let duplicate_result = heap.with_slice_mut_and_heap_mut(prepared_storage, |values, heap| {
            let destination = values
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            OrigAtData_Duplicate(heap, destination, source)
        })?;
        context.PrepAtData.clone_from_slice(&heap.slice(prepared_storage.as_const())?[..2]);
        if duplicate_result != 0 {
            return Ok(_IS_ERROR as i32);
        }
    }

    context.num_inp = input_number;
    context.save_opt_bits = save_option_bits;
    let string_storage = heap.allocate_model_storage(vec![INCHI_IOS_STRING::default()])?;
    context.strbuf = string_storage;
    let string_result = heap.with_slice_mut_and_heap_mut(string_storage, |values, heap| {
        let destination = values
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if let Some(string_buffer) = string_buffer {
            inchi_strbuf_create_copy(heap, destination, string_buffer)
        } else {
            inchi_strbuf_init(
                heap,
                destination,
                INCHI_STRBUF_INITIAL_SIZE as i32,
                INCHI_STRBUF_SIZE_INCREMENT as i32,
            )
        }
    })?;
    context.temp_string_container = heap.slice(string_storage.as_const())?[0].clone();
    if string_result == -1 {
        return Ok(_IS_FATAL as i32);
    }

    Ok(_IS_OKAY as i32)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn POSEContext_Free(
    heap: &mut SourceHeap,
    context: &mut POSEContext,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1636 POSEContext_Free
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
void POSEContext_Free(POSEContext *context)
{
    int i;
    for (i = 0; i < MAX_NUM_PATHS; i++)
    {
        if (context->ip.path[i])
        {
            inchi_free((void*)context->ip.path[i]);
            /*  cast deliberately discards 'const' qualifier */
            context->ip.path[i] = NULL;
        }
    }
    FreeAllINChIArrays(context->pINChI2, context->pINChI_Aux2, context->sd.num_components);
    if (context->inp_file)
    {
        ;
    }
    else
    {
        ;
    }
    inchi_ios_close(context->out_file);
    inchi_ios_close(context->log_file);
    inchi_ios_close(context->prb_file);
    FreeOrigAtData(context->orig_inp_data);
    FreeOrigAtData(context->prep_inp_data);
    FreeOrigAtData( context->prep_inp_data+1);
    context->num_inp = 0;
    context->save_opt_bits = 0;
    inchi_strbuf_close(context->strbuf);

    return;
}
    */
    // END INCHI C FUNCTION: POSEContext_Free
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: POSEContext_Free
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; MAX_NUM_PATHS=4 and INCHI_NUM=2.
    // INCHI✔️❌: inchi_free is the GCC/Linux macro expansion to free; all cleanup callees dispatch
    // INCHI✔️❌: to their independently ported active definitions. Model-storage self-pointer mirrors
    // INCHI✔️❌: add allocation-map and snapshot synchronization overhead absent from the C object.
    // END INCHI ACTIVE MACRO CONFIGURATION: POSEContext_Free

    for path in &mut context.ip.path {
        if !path.is_null() {
            inchi_free(heap, path.as_mut())?;
            *path = SourceConstPointer::null();
        }
    }

    FreeAllINChIArrays(
        heap,
        &mut context.pINChI2,
        &mut context.pINChI_Aux2,
        &mut context.sd.num_components,
    )?;

    if !context.out_file.is_null() {
        heap.with_slice_mut_and_heap_mut(context.out_file, |streams, heap| {
            let streams = streams
                .get_mut(..3)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            for stream in streams.iter_mut() {
                inchi_ios_close(heap, Some(stream))?;
            }
            Ok(())
        })?;
        context
            .inchi_file
            .clone_from_slice(&heap.slice(context.out_file.as_const())?[..3]);
    }

    if !context.orig_inp_data.is_null() {
        heap.with_slice_mut_and_heap_mut(context.orig_inp_data, |values, heap| {
            FreeOrigAtData(heap, values.first_mut())
        })?;
        context.OrigAtData = heap.slice(context.orig_inp_data.as_const())?[0].clone();
    }

    if !context.prep_inp_data.is_null() {
        heap.with_slice_mut_and_heap_mut(context.prep_inp_data, |values, heap| {
            let values = values
                .get_mut(..2)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let (first, second) = values.split_at_mut(1);
            FreeOrigAtData(heap, first.first_mut())?;
            FreeOrigAtData(heap, second.first_mut())
        })?;
        context
            .PrepAtData
            .clone_from_slice(&heap.slice(context.prep_inp_data.as_const())?[..2]);
    }

    context.num_inp = 0;
    context.save_opt_bits = 0;

    if !context.strbuf.is_null() {
        heap.with_slice_mut_and_heap_mut(context.strbuf, |values, heap| {
            inchi_strbuf_close(heap, values.first_mut())
        })?;
        context.temp_string_container = heap.slice(context.strbuf.as_const())?[0].clone();
    }

    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn extract_orig_nums_from_auxinfo_string(
    heap: &mut SourceHeap,
    saux: SourceConstPointer<i8>,
    orig: &mut [i32],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1389 extract_orig_nums_from_auxinfo_string
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
int extract_orig_nums_from_auxinfo_string(char *saux, int *orig)
{
    int res = _IS_OKAY;
    int k, c = 0, cano_num = 1 /*0*/;
    const char *p, *q;

    p = strstr(saux, "/N:");  /* must always be there */
    if (!p || !p[3] || !isdigit(UCINT p[3]))
    {
        res = _IS_ERROR;
        goto exit_function;
    }

    p += 3;
    q = p;

    while ((k = inchi_strtol(p, &q, 10))) /* djb-rwth: addressing LLVM warning */
    {
        orig[cano_num++] = k/* - 1*/; /* 1-based numbers */
        if ((c = UCINT *q) && c != '/') /* djb-rwth: addressing LLVM warning */
        {
            p = q + 1;
        }
        else
        {
            break;
        }
    }

exit_function:

    return res;
}
    */
    // END INCHI C FUNCTION: extract_orig_nums_from_auxinfo_string
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: extract_orig_nums_from_auxinfo_string
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; C locale.
    // INCHI✔️✔️: isdigit receives UCINT and accepts only ASCII decimal digits here.
    // INCHI✔️✔️: LP64 long values narrow to GCC int two's-complement low bits before loop testing.
    // END INCHI ACTIVE MACRO CONFIGURATION: extract_orig_nums_from_auxinfo_string

    let bytes = heap.slice(saux)?;
    let nul = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let header = bytes[..nul]
        .windows(3)
        .position(|window| window == [b'/' as i8, b'N' as i8, b':' as i8]);
    let Some(header) = header else {
        return Ok(_IS_ERROR as i32);
    };
    let first = header
        .checked_add(3)
        .ok_or(SourceHeapError::PointerOffsetOverflow)?;
    if first >= nul || !(bytes[first] as u8).is_ascii_digit() {
        return Ok(_IS_ERROR as i32);
    }

    let mut p = saux.offset(
        i64::try_from(first).map_err(|_| SourceHeapError::PointerOffsetOverflow)?,
    )?;
    let mut canonical_number = 1_i32;
    loop {
        let mut q = p;
        let k = inchi_strtol(heap, p, Some(&mut q), 10)? as i32;
        if k == 0 {
            break;
        }
        let index = usize::try_from(canonical_number)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        *orig
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = k;
        canonical_number = canonical_number.wrapping_add(1);

        let c = *heap
            .slice(q)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)? as u8;
        if c != 0 && c != b'/' {
            p = q.offset(1)?;
        } else {
            break;
        }
    }
    Ok(crate::source_types::_IS_OKAY as i32)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn extract_nonstereo_eq_classes_from_auxinfo_string(
    heap: &mut SourceHeap,
    saux: SourceConstPointer<i8>,
    nat: i32,
    orig: &[i32],
    nclasses: &mut i32,
    eclass: &mut [i32],
    eclass_by_origs: &mut [i32],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1427 extract_nonstereo_eq_classes_from_auxinfo_string
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
int extract_nonstereo_eq_classes_from_auxinfo_string( char *saux,
                                                      int nat,
                                                      int *orig,
                                                      int *nclasses,
                                                      int *eclass,
                                                      int *eclass_by_origs)
{
    int res = _IS_OKAY;
    int k, c = 0, cano_num = 1, orig_num = 1;
    const char *p, *q;

    /* Note that all atom and class numbers here are 1-based */

    *nclasses = 0;
    memset(eclass, -1, ((long long)nat+1) * sizeof(int)); /* djb-rwth: cast operator added; memset_s C11/Annex K variant? */
    memset(eclass_by_origs, -1, ((long long)nat+1) * sizeof(int)); /* djb-rwth: cast operator added; memset_s C11/Annex K variant? */

    p = strstr(saux, "/E:");
    if (!p)
    {
        /* No "/E" means that all atoms are different  */
        return res;
    }

    p += 3;
    q = p;
    while ((k = (AT_NUMB)inchi_strtol(p + 1, &q, 10))) /* djb-rwth: addressing LLVM warning */
    {
        c = UCINT *q;
        if (c == '/')
        {
            break;
        }
        else if (c == ',' || c == ')')
        {
            eclass[k] = *nclasses;
            if (c == ')')
            {
                (*nclasses)++;
                q++;
                c = UCINT *q;
                if (c == '/')
                    break;
                else
                    ;
            }
            p = q;
        }
        else
        {
            return _IS_ERROR;
        }
    }
    /* NB: cano, origs start from 0 */
    for (cano_num = 1; cano_num <= nat; cano_num++)
    {
        if (eclass[cano_num] == -1) /* the atom is unique, add one more eq class for him */
        {
            (*nclasses)++;
            eclass[cano_num] = *nclasses;
        }
    }

    for (cano_num = 1; cano_num <= nat; cano_num++)
    {
        orig_num = orig[cano_num];  /* NB: cano, origs start from 0 */
        eclass_by_origs[orig_num] = eclass[cano_num];
    }

    return res;
}
    */
    // END INCHI C FUNCTION: extract_nonstereo_eq_classes_from_auxinfo_string
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: extract_nonstereo_eq_classes_from_auxinfo_string
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; C locale.
    // INCHI✔️✔️: AT_NUMB is unsigned short, so each LP64 long is narrowed modulo 2^16.
    // INCHI✔️✔️: UCINT reads delimiters as unsigned char; memset writes nat+1 int slots in call order.
    // END INCHI ACTIVE MACRO CONFIGURATION: extract_nonstereo_eq_classes_from_auxinfo_string

    *nclasses = 0;
    let count = i64::from(nat)
        .checked_add(1)
        .and_then(|value| usize::try_from(value).ok())
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    eclass
        .get_mut(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .fill(-1);
    eclass_by_origs
        .get_mut(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .fill(-1);

    let header = {
        let bytes = heap.slice(saux)?;
        let nul = bytes
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        bytes[..nul]
            .windows(3)
            .position(|window| window == [b'/' as i8, b'E' as i8, b':' as i8])
    };
    let Some(header) = header else {
        return Ok(_IS_OKAY as i32);
    };

    let start = header
        .checked_add(3)
        .ok_or(SourceHeapError::PointerOffsetOverflow)?;
    let mut p = saux.offset(
        i64::try_from(start).map_err(|_| SourceHeapError::PointerOffsetOverflow)?,
    )?;
    loop {
        let parse_start = p.offset(1)?;
        let mut q = p;
        let k = i32::from(inchi_strtol(heap, parse_start, Some(&mut q), 10)? as u16);
        if k == 0 {
            break;
        }
        let mut c = *heap
            .slice(q)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)? as u8;
        if c == b'/' {
            break;
        } else if c == b',' || c == b')' {
            *eclass
                .get_mut(usize::try_from(k).map_err(|_| SourceHeapError::SourceIntegerOverflow)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = *nclasses;
            if c == b')' {
                *nclasses = nclasses.wrapping_add(1);
                q = q.offset(1)?;
                c = *heap
                    .slice(q)?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)? as u8;
                if c == b'/' {
                    break;
                }
            }
            p = q;
        } else {
            return Ok(_IS_ERROR as i32);
        }
    }

    for canonical_number in 1..=nat {
        let index = usize::try_from(canonical_number)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let class = eclass
            .get_mut(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if *class == -1 {
            *nclasses = nclasses.wrapping_add(1);
            *class = *nclasses;
        }
    }

    for canonical_number in 1..=nat {
        let canonical_index = usize::try_from(canonical_number)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let original_number = *orig
            .get(canonical_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let original_index = usize::try_from(original_number)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        *eclass_by_origs
            .get_mut(original_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = eclass[canonical_index];
    }

    Ok(_IS_OKAY as i32)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn DiylFrag_New(
    heap: &mut SourceHeap,
    na: i32,
    end1: i32,
    end2: i32,
    s: SourceConstPointer<i8>,
) -> Result<SourceMutPointer<DiylFrag>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2007 DiylFrag_New
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
DiylFrag* DiylFrag_New(int na, int end1, int end2, char *s)
{
    int err = 0;

    DiylFrag *pfrag = NULL;

    pfrag = (DiylFrag *)inchi_calloc(1, sizeof(DiylFrag));
    if (NULL == pfrag)
    {
        err = 1;
        goto exit_function;
    }

    pfrag->na = na; 
    pfrag->end1 = end1;
    pfrag->end2 = end2;
    pfrag->alist = NULL;
    pfrag->xclist = NULL;

    if (na > 0 )
    {
        pfrag->alist = (int *)inchi_calloc(na, sizeof(int));
        pfrag->xclist = (int *)inchi_calloc(na, sizeof(int));
        if (!pfrag->alist || !pfrag->xclist)
        {
            err = 2;
            goto exit_function;
        }
    }

    inchi_strbuf_printf(&pfrag->sig, "%-s", s);

exit_function:
    if (err)
    {
        DiylFrag_Free(pfrag);
        inchi_free(pfrag); /* djb-rwth: addressing coverity ID #499507 */
        return NULL;
    }
    return pfrag;
}
    */
    // END INCHI C FUNCTION: DiylFrag_New
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: DiylFrag_New
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; sizeof(DiylFrag)=56.
    // INCHI✔️❌: Both atom arrays are allocated before either null result is tested; cleanup order is DiylFrag_Free then owner free.
    // INCHI✔️❌: The sole active caller passes an empty signature; nonempty text exposes the source formatter's one-byte NUL overrun.
    // INCHI✔️❌: The source ignores the signature formatter return; a formatter allocation failure reaches undefined null-write behavior.
    // INCHI✔️❌: The safe pointer model reports that undefined formatter path as a structured error after complete cleanup.
    // INCHI✔️❌: SourceHeap allocation-map operations add overhead absent from direct C allocations.
    // END INCHI ACTIVE MACRO CONFIGURATION: DiylFrag_New

    let pfrag: SourceMutPointer<DiylFrag> = match inchi_calloc(heap, 1, 56) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return Ok(SourceMutPointer::null()),
        Err(error) => return Err(error),
    };
    {
        let fragment = heap
            .slice_mut(pfrag)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        fragment.na = na;
        fragment.end1 = end1;
        fragment.end2 = end2;
        fragment.alist = SourceMutPointer::null();
        fragment.xclist = SourceMutPointer::null();
    }

    if na > 0 {
        let alist_result = inchi_calloc(heap, na as u64, 4);
        let xclist_result = inchi_calloc(heap, na as u64, 4);
        let unexpected_error = match (&alist_result, &xclist_result) {
            (Err(error), _) if *error != SourceHeapError::AllocationFailed => Some(*error),
            (_, Err(error)) if *error != SourceHeapError::AllocationFailed => Some(*error),
            _ => None,
        };
        if let Some(error) = unexpected_error {
            if let Ok(pointer) = &alist_result {
                inchi_free(heap, *pointer)?;
            }
            if let Ok(pointer) = &xclist_result {
                inchi_free(heap, *pointer)?;
            }
            heap.with_slice_mut_and_heap_mut(pfrag, |fragments, heap| {
                DiylFrag_Free(heap, fragments.first_mut())
            })?;
            inchi_free(heap, pfrag)?;
            return Err(error);
        }
        let alist = match alist_result {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        let xclist = match xclist_result {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        {
            let fragment = heap
                .slice_mut(pfrag)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            fragment.alist = alist;
            fragment.xclist = xclist;
        }
        if alist.is_null() || xclist.is_null() {
            heap.with_slice_mut_and_heap_mut(pfrag, |fragments, heap| {
                DiylFrag_Free(heap, fragments.first_mut())
            })?;
            inchi_free(heap, pfrag)?;
            return Ok(SourceMutPointer::null());
        }
    }

    let format = heap.allocate_model_storage(vec![b'%' as i8, b'-' as i8, b's' as i8, 0])?;
    let print_result = heap.with_slice_mut_and_heap_mut(pfrag, |fragments, heap| {
        let fragment = fragments
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        inchi_strbuf_printf(
            heap,
            Some(&mut fragment.sig),
            format.as_const(),
            &SourceVaList {
                arguments: vec![SourceFormatArgument::Bytes(s)],
                position: 0,
            },
        )
    });
    heap.free(format)?;
    if let Err(error) = print_result {
        heap.with_slice_mut_and_heap_mut(pfrag, |fragments, heap| {
            DiylFrag_Free(heap, fragments.first_mut())
        })?;
        inchi_free(heap, pfrag)?;
        return Err(error);
    }

    Ok(pfrag)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn DiylFrag_Free(
    heap: &mut SourceHeap,
    pfrag: Option<&mut DiylFrag>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2049 DiylFrag_Free
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
void DiylFrag_Free(DiylFrag *pfrag)
{
    if (!pfrag)
    {
        return;
    }
    if (pfrag->alist)
    {
        inchi_free(pfrag->alist);
        pfrag->alist = NULL;
    }
    if (pfrag->xclist)
    {
        inchi_free(pfrag->xclist);
        pfrag->xclist = NULL;
    }
    inchi_strbuf_close(&pfrag->sig);
    return;
}
    */
    // END INCHI C FUNCTION: DiylFrag_Free
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: DiylFrag_Free
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64.
    // INCHI✔️❌: inchi_free is the GCC/Linux free macro; inchi_strbuf_close frees pStr and zeroes sig.
    // INCHI✔️❌: SourceHeap allocation-map lookups add overhead absent from direct C frees.
    // END INCHI ACTIVE MACRO CONFIGURATION: DiylFrag_Free

    let Some(pfrag) = pfrag else {
        return Ok(());
    };
    if !pfrag.alist.is_null() {
        inchi_free(heap, pfrag.alist)?;
        pfrag.alist = SourceMutPointer::null();
    }
    if !pfrag.xclist.is_null() {
        inchi_free(heap, pfrag.xclist)?;
        pfrag.xclist = SourceMutPointer::null();
    }
    inchi_strbuf_close(heap, Some(&mut pfrag.sig))?;
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn DiylFrag_MakeSignature(
    heap: &mut SourceHeap,
    pfrag: &mut DiylFrag,
    nxc: i32,
    xc: &[i32],
    cnt: &mut [i32],
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2069 DiylFrag_MakeSignature
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
void DiylFrag_MakeSignature(DiylFrag *pfrag, 
                            int nxc,            /* n xclasses (molecule-wide)       */
                            int *xc,            /* xclasses (molecule-wide)         */
                            int *cnt )          /* temp storage: counts of xclasses */
{
    int i, k, nxc_frag; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    
    inchi_strbuf_printf(&pfrag->sig, "%-d,%-d,%-d{", pfrag->na, xc[pfrag->end1], xc[pfrag->end2]);
    for (i = 0; i < pfrag->na; i++)
    {
        pfrag->xclist[i] = xc[pfrag->alist[i]];
    }  
    nxc_frag = count_colors_in_sequence(pfrag->xclist, pfrag->na, nxc+1, cnt); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    for (k = 0; k < nxc; k++)
    {
        if (cnt[k] > 0)
        {
            /* (xclass:cnt)*/
            inchi_strbuf_printf(&pfrag->sig, "(%-d:%-d)", k, cnt[k]);
        }
    }

    inchi_strbuf_printf(&pfrag->sig, "}");

    return;
}
    */
    // END INCHI C FUNCTION: DiylFrag_MakeSignature
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: DiylFrag_MakeSignature
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; sizeof(int)=4.
    // INCHI✔️❌: The completed inchi_strbuf_printf and count_colors_in_sequence callees preserve formatting, append, clear, and wrapping semantics.
    // INCHI✔️❌: The output loop intentionally excludes class nxc even though the counting call accepts it through maxcol=nxc+1.
    // INCHI✔️❌: SourceHeap format-literal allocations and allocation-map accesses are materially slower than C static literals and direct pointers.
    // INCHI✔️❌: Invalid pointers and the source formatter's unchecked allocation-failure writes are reported as structured errors after the same defined mutation prefix.
    // END INCHI ACTIVE MACRO CONFIGURATION: DiylFrag_MakeSignature

    let first_end = usize::try_from(pfrag.end1)
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let second_end = usize::try_from(pfrag.end2)
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let first_color = *xc
        .get(first_end)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second_color = *xc
        .get(second_end)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;

    let header_format = heap.allocate_model_storage(
        b"%-d,%-d,%-d{\0"
            .iter()
            .map(|byte| *byte as i8)
            .collect(),
    )?;
    let header_result = inchi_strbuf_printf(
        heap,
        Some(&mut pfrag.sig),
        header_format.as_const(),
        &SourceVaList {
            arguments: vec![
                SourceFormatArgument::Signed(i64::from(pfrag.na)),
                SourceFormatArgument::Signed(i64::from(first_color)),
                SourceFormatArgument::Signed(i64::from(second_color)),
            ],
            position: 0,
        },
    );
    heap.free(header_format)?;
    header_result?;

    if pfrag.na > 0 {
        let atom_count = usize::try_from(pfrag.na)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let alist_pointer = pfrag.alist;
        heap.with_slice_mut_and_heap(pfrag.xclist, |xclist, heap| {
            let alist = heap.slice(alist_pointer.as_const())?;
            for index in 0..atom_count {
                let atom = *alist
                    .get(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let atom = usize::try_from(atom)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                *xclist
                    .get_mut(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = *xc
                    .get(atom)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
            }
            Ok(())
        })?;
    }

    let maximum_colors = nxc
        .checked_add(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let _number_of_fragment_colors = if pfrag.na > 0 {
        count_colors_in_sequence(
            heap.slice(pfrag.xclist.as_const())?,
            pfrag.na,
            maximum_colors,
            cnt,
        )?
    } else {
        count_colors_in_sequence(&[], pfrag.na, maximum_colors, cnt)?
    };

    for color in 0..nxc {
        let color_index = usize::try_from(color)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let count = *cnt
            .get(color_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if count > 0 {
            let count_format = heap.allocate_model_storage(
                b"(%-d:%-d)\0"
                    .iter()
                    .map(|byte| *byte as i8)
                    .collect(),
            )?;
            let count_result = inchi_strbuf_printf(
                heap,
                Some(&mut pfrag.sig),
                count_format.as_const(),
                &SourceVaList {
                    arguments: vec![
                        SourceFormatArgument::Signed(i64::from(color)),
                        SourceFormatArgument::Signed(i64::from(count)),
                    ],
                    position: 0,
                },
            );
            heap.free(count_format)?;
            count_result?;
        }
    }

    let close_format = heap.allocate_model_storage(vec![b'}' as i8, 0])?;
    let close_result = inchi_strbuf_printf(
        heap,
        Some(&mut pfrag.sig),
        close_format.as_const(),
        &SourceVaList::default(),
    );
    heap.free(close_format)?;
    close_result?;
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn DiylFrag_Diff(
    heap: &SourceHeap,
    pfrag1: &DiylFrag,
    pfrag2: &DiylFrag,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2098 DiylFrag_Diff
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int DiylFrag_Diff(DiylFrag *pfrag1, DiylFrag *pfrag2)
{
    if (pfrag1->na != pfrag2->na)
    {
        return 1;
    }
    if (pfrag1->nb != pfrag2->nb)
    {
        return 1;
    }
    if (pfrag1->sig.nUsedLength && pfrag2->sig.nUsedLength)
    {
        int cmp = strcmp(pfrag1->sig.pStr, pfrag2->sig.pStr);
        return cmp;
    }

    return 0;
}
    */
    // END INCHI C FUNCTION: DiylFrag_Diff
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: DiylFrag_Diff
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; strcmp compares unsigned char values.
    // INCHI✔️❌: nUsedLength is only a truth test; either zero suppresses every pStr read, while negative values still select strcmp.
    // INCHI✔️❌: The direct byte loop preserves strcmp O(n) time and O(1) space, but SourceHeap BTreeMap pointer resolution adds logarithmic lookup overhead.
    // INCHI✔️❌: Non-string pointer states that make source strcmp undefined return structured pointer or missing-NUL errors.
    // END INCHI ACTIVE MACRO CONFIGURATION: DiylFrag_Diff

    if pfrag1.na != pfrag2.na {
        return Ok(1);
    }
    if pfrag1.nb != pfrag2.nb {
        return Ok(1);
    }
    if pfrag1.sig.nUsedLength != 0 && pfrag2.sig.nUsedLength != 0 {
        let first = heap.slice(pfrag1.sig.pStr.as_const())?;
        let second = heap.slice(pfrag2.sig.pStr.as_const())?;
        let mut index = 0_usize;
        loop {
            let first_byte = *first
                .get(index)
                .ok_or(SourceHeapError::MissingNulTerminator)? as u8;
            let second_byte = *second
                .get(index)
                .ok_or(SourceHeapError::MissingNulTerminator)? as u8;
            if first_byte != second_byte {
                return Ok(i32::from(first_byte) - i32::from(second_byte));
            }
            if first_byte == 0 {
                return Ok(0);
            }
            index = index
                .checked_add(1)
                .ok_or(SourceHeapError::PointerOffsetOverflow)?;
        }
    }
    Ok(0)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn DiylFrag_DebugTrace(
    pfrag: Option<&DiylFrag>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2119 DiylFrag_DebugTrace
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
void DiylFrag_DebugTrace(DiylFrag *pfrag)
{
    int k, na;

    if (!pfrag)
    {
        return;
    }
    
    ITRACE_("DiylFrag @ %-p ", pfrag);
    na = pfrag->na;
    ITRACE_("\n\t%-d atoms. List of atoms and their xclasses : { ", na);
    for (k = 0; k < na - 1; k++)
    {
        ITRACE_(" %-d(%-d), ", pfrag->alist[k], pfrag->xclist[k]);
    }
    ITRACE_(" %-d(%-d) }\n", pfrag->alist[na - 1], pfrag->xclist[na - 1]);

    ITRACE_("\tend1 = %-d, end2 = %-d, nb = %-d\n", pfrag->end1, pfrag->end2, pfrag->nb);
    
    ITRACE_("\tSignature = '%-s'\n", pfrag->sig.pStr);

    return;
}
    */
    // END INCHI C FUNCTION: DiylFrag_DebugTrace
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: DiylFrag_DebugTrace
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux, so ichi_io.h defines ITRACE_ as 0 && _inchi_trace.
    // INCHI✔️✔️: Every trace argument is unevaluated; alist, xclist, end fields, nb, and sig.pStr are never read.
    // INCHI✔️✔️: The null check, na load, na-1 arithmetic, empty O(max(na,0)) loop, and return remain active with no allocation or output.
    // END INCHI ACTIVE MACRO CONFIGURATION: DiylFrag_DebugTrace

    let Some(pfrag) = pfrag else {
        return Ok(());
    };
    let na = pfrag.na;
    let loop_limit = na
        .checked_sub(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    for _k in 0..loop_limit {
        // ITRACE_ is short-circuited, so the source loop body evaluates nothing.
    }
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OAD_Polymer_PrepareFoldCRUEdits(
    heap: &mut SourceHeap,
    original_atom_data: &mut ORIG_ATOM_DATA,
    sinchi_noedits: SourceConstPointer<i8>,
    _saux_noedits: SourceConstPointer<i8>,
    sinchi: SourceConstPointer<i8>,
    saux: SourceConstPointer<i8>,
    edits: &OAD_StructureEdits,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1838 OAD_Polymer_PrepareFoldCRUEdits
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int  OAD_Polymer_PrepareFoldCRUEdits( ORIG_ATOM_DATA *orig_at_data,
                                      char *sinchi_noedits, 
                                      char *saux_noedits,
                                      char *sinchi,
                                      char *saux,
                                      OAD_StructureEdits *ed)
{
    int ret = _IS_OKAY;
    int i, j, k;
    int err;
    char pStrErr[STR_ERR_LEN];
    int *orig = NULL;
    int nat = orig_at_data->num_inp_atoms;
    int neclasses = 0;		/* No of constitutional equivalence classses for the atoms		*/
    int nxclasses = 0;      /* No of extended (stereo-aware) atom classses == 3*neclasses   */
    
    int *all_bkb_orig = NULL, n_all_bkb_orig = 0;
    OAD_Polymer *p = orig_at_data->polymer;
    int nu = orig_at_data->polymer->n;

    /* Extract cano_nums-->orig_nums mapping from AuxInfo AuxInfo Main Layer */
    orig = (int*)inchi_calloc((long long)nat + 1, sizeof(int)); /* djb-rwth: cast operator added */
    if (!orig)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    ret = extract_orig_nums_from_auxinfo_string(saux, orig);
    if (ret != _IS_OKAY && ret != _IS_WARNING)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    /* Extract non-stereo eq. classes data from AuxInfo */
    ret = extract_nonstereo_eq_classes_from_auxinfo_string(saux, nat, orig, &neclasses, ec_cano_opp, ec_opp);
    if (ret != _IS_OKAY && ret != _IS_WARNING)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    if (neclasses == 0)
    {
        goto exit_function;
    }
    /* Extract stereocenter data from InChI */

    /*ret = extract_stereo_info_from_inchi_string(sinchi, nat, orig, at_stereo_mark_orig);*/
    ret = extract_stereo_info_from_inchi_string(sinchi_noedits, nat, orig, at_stereo_mark_orig_opp);
    if (ret != _IS_OKAY && ret != _IS_WARNING)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    /* Make extended stereo-aware atom classes */
    nxclasses = neclasses * 3;
    for (i = 1; i <= nat; i++)  /* orig # */
    {
        int atom_class = ec_opp[i];

        if (at_stereo_mark_orig_opp[i] == INCHI_PARITY_ODD)
        {
            atom_class += neclasses;
        }
        else if (at_stereo_mark_orig_opp[i] == INCHI_PARITY_EVEN)
        {
            atom_class += 2 * neclasses;
        }
        xc_opp[i] = atom_class;
    }
    /* Extract all backbone bonds, in all units, from InChI (z layer).
        NB: we assume that units are not 'inter-crossing' so
        any particular bkbond belongs to some unique CRU.
    */
    all_bkb_orig = (int*)inchi_calloc(2 * ((long long)orig_at_data->num_inp_bonds + 1), sizeof(int)); /* djb-rwth: cast operator added */
    if (!all_bkb_orig)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    memset(all_bkb_orig, 0, ((long long)orig_at_data->num_inp_bonds + 1) * sizeof(int)); /* djb-rwth: cast operator added; memset_s C11/Annex K variant? */
    ret = extract_all_backbone_bonds_from_inchi_string(sinchi, &n_all_bkb_orig, orig, all_bkb_orig);
    if (ret != _IS_OKAY && ret != _IS_WARNING)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    /* just for case, remove those bkbonds which are not single (alternate may be here) */
    for (k = n_all_bkb_orig - 1; k >= 0; k--)
    {
        int orig1 = all_bkb_orig[2 * k];
        int orig2 = all_bkb_orig[2 * k + 1];
        int bond_type = Inp_Atom_GetBondType(orig_at_data->at, orig1 - 1, orig2 - 1);
        if (bond_type > BOND_TYPE_SINGLE) /* not == intentionally, to keep -1 ("no bond") */
        {
            /* remove k-th bond and shift others to start of list */
            int kk;
            for (kk = k; kk < n_all_bkb_orig; kk++)
            {
                all_bkb_orig[2 * kk] = all_bkb_orig[2 * (kk + 1)];
                all_bkb_orig[2 * kk + 1] = all_bkb_orig[2 * (kk + 1) + 1];
            }
            all_bkb_orig[2 * n_all_bkb_orig] = 0;
            all_bkb_orig[2 * n_all_bkb_orig + 1] = 0;
            n_all_bkb_orig--;
        }
    }

    err = OAD_ValidatePolymerAndPseudoElementData(orig_at_data,
        POLYMERS_MODERN,
        1, /* ip->bNPZz,*/
        pStrErr,
        0 /*ip->bNoWarnings*/);
    if (err)
    {
        goto exit_function;
    }

    /* For each unit analyze a possibility of folding (i.e., removal of excess in-CRU repeats) */
    for (j = 0; j < nu; j++)
    {
        OAD_PolymerUnit* u = p->units[j];

        if (u->na < 2)
        {
            goto nextj;
        }
        if (u->nb < 2)
        {
            goto nextj;
        }
        /* this is only for bi-star CRU's */
        if (!u->cap1_is_undef)
        {
            goto nextj;
        }
        if (!u->cap2_is_undef)
        {
            goto nextj;
        }

        err = analyze_CRU_folding(orig_at_data, j,
            n_all_bkb_orig, all_bkb_orig,
            nxclasses, xc_opp,
            ed);
        if (err)
        {
            ret = inchi_max(_IS_WARNING, err);
            goto nextj;
        }

    nextj:;
    }

exit_function:
    if (orig)
    {
        inchi_free(orig);
    }
    if (all_bkb_orig)
    {
        inchi_free(all_bkb_orig);
    }

    return ret;
}
    */
    // END INCHI C FUNCTION: OAD_Polymer_PrepareFoldCRUEdits
    // BEGIN INCHI ACTIVE GLOBAL STORAGE: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1827
    // INCHI✔️❌: int ec_opp[MAX_ATOMS], ec_cano_opp[MAX_ATOMS], at_stereo_mark_orig_opp[MAX_ATOMS], xc_opp[MAX_ATOMS];
    // END INCHI ACTIVE GLOBAL STORAGE: OAD_Polymer_PrepareFoldCRUEdits
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OAD_Polymer_PrepareFoldCRUEdits
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; sizeof(int)=4 and sizeof(long long)=8.
    // INCHI✔️❌: `saux_noedits` is an active formal parameter but is never evaluated by the function body.
    // INCHI✔️❌: All direct callees are completed source ports; `inchi_max` selects the greater signed int.
    // INCHI✔️❌: Rust-local mirrors of the four process-global arrays avoid shared mutation but add allocation and initialization work absent from each C call.
    // INCHI✔️❌: SourceHeap pointer-map access and temporary clones add further material overhead beyond native pointer indexing.
    // END INCHI ACTIVE MACRO CONFIGURATION: OAD_Polymer_PrepareFoldCRUEdits

    fn is_allocation_error(error: SourceHeapError) -> bool {
        matches!(
            error,
            SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange
                | SourceHeapError::AllocationFailed
        )
    }
    fn global_array() -> Result<Vec<i32>, SourceHeapError> {
        let mut values = Vec::new();
        values
            .try_reserve_exact(MAX_ATOMS as usize)
            .map_err(|_| SourceHeapError::AllocationFailed)?;
        values.resize(MAX_ATOMS as usize, 0);
        Ok(values)
    }
    fn accepted(ret: i32) -> bool {
        ret == _IS_OKAY as i32 || ret == _IS_WARNING as i32
    }
    fn pair_offset(index: i32) -> Result<usize, SourceHeapError> {
        usize::try_from(index)
            .ok()
            .and_then(|value| value.checked_mul(2))
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }

    let polymer = heap
        .slice(original_atom_data.polymer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let number_of_units = polymer.n;
    let number_of_atoms = original_atom_data.num_inp_atoms;

    let mut ec_by_original = global_array()?;
    let mut ec_by_canonical = global_array()?;
    let mut stereo_by_original = global_array()?;
    let mut xclass_by_original = global_array()?;
    let mut original_numbers = SourceMutPointer::null();
    let mut all_backbone_bonds = SourceMutPointer::null();

    let execution = (|| -> Result<i32, SourceHeapError> {
        let original_count = i64::from(number_of_atoms).wrapping_add(1) as u64;
        original_numbers = match inchi_calloc::<i32>(heap, original_count, 4) {
            Ok(pointer) => pointer,
            Err(error) if is_allocation_error(error) => return Ok(_IS_ERROR as i32),
            Err(error) => return Err(error),
        };

        let mut ret = heap.with_slice_mut_and_heap_mut(original_numbers, |values, heap| {
            extract_orig_nums_from_auxinfo_string(heap, saux, values)
        })?;
        if !accepted(ret) {
            return Ok(_IS_ERROR as i32);
        }

        let original_map = heap.slice(original_numbers.as_const())?.to_vec();
        let mut number_of_classes = 0_i32;
        ret = extract_nonstereo_eq_classes_from_auxinfo_string(
            heap,
            saux,
            number_of_atoms,
            &original_map,
            &mut number_of_classes,
            &mut ec_by_canonical,
            &mut ec_by_original,
        )?;
        if !accepted(ret) {
            return Ok(_IS_ERROR as i32);
        }
        if number_of_classes == 0 {
            return Ok(ret);
        }

        ret = extract_stereo_info_from_inchi_string(
            heap,
            sinchi_noedits,
            number_of_atoms,
            &original_map,
            &mut stereo_by_original,
        )?;
        if !accepted(ret) {
            return Ok(_IS_ERROR as i32);
        }

        let number_of_xclasses = number_of_classes
            .checked_mul(3)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        if number_of_atoms > 0 {
            for atom_number in 1..=number_of_atoms {
                let index = usize::try_from(atom_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let mut atom_class = *ec_by_original
                    .get(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let parity = *stereo_by_original
                    .get(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if parity == INCHI_PARITY_ODD as i32 {
                    atom_class = atom_class.wrapping_add(number_of_classes);
                } else if parity == INCHI_PARITY_EVEN as i32 {
                    atom_class = atom_class.wrapping_add(number_of_classes.wrapping_mul(2));
                }
                *xclass_by_original
                    .get_mut(index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = atom_class;
            }
        }

        let backbone_slots = i64::from(original_atom_data.num_inp_bonds)
            .wrapping_add(1)
            .wrapping_mul(2) as u64;
        all_backbone_bonds = match inchi_calloc::<i32>(heap, backbone_slots, 4) {
            Ok(pointer) => pointer,
            Err(error) if is_allocation_error(error) => return Ok(_IS_ERROR as i32),
            Err(error) => return Err(error),
        };

        let mut number_of_backbone_bonds = 0_i32;
        ret = heap.with_slice_mut_and_heap_mut(all_backbone_bonds, |values, heap| {
            extract_all_backbone_bonds_from_inchi_string(
                heap,
                sinchi,
                &mut number_of_backbone_bonds,
                &original_map,
                values,
            )
        })?;
        if !accepted(ret) {
            return Ok(_IS_ERROR as i32);
        }

        let mut k = number_of_backbone_bonds.wrapping_sub(1);
        while k >= 0 {
            let offset = pair_offset(k)?;
            let (original1, original2) = {
                let bonds = heap.slice(all_backbone_bonds.as_const())?;
                (
                    *bonds.get(offset).ok_or(SourceHeapError::PointerOutOfBounds)?,
                    *bonds
                        .get(offset + 1)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                )
            };
            let first = original1
                .checked_sub(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            let second = original2
                .checked_sub(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            let bond_type = Inp_Atom_GetBondType(
                heap,
                original_atom_data.at.as_const(),
                first,
                second,
            )?;
            if bond_type > BOND_TYPE_SINGLE as i32 {
                let mut kk = k;
                while kk < number_of_backbone_bonds {
                    let target = pair_offset(kk)?;
                    let source = target
                        .checked_add(2)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let source_pair = {
                        let bonds = heap.slice(all_backbone_bonds.as_const())?;
                        (
                            *bonds.get(source).ok_or(SourceHeapError::PointerOutOfBounds)?,
                            *bonds
                                .get(source + 1)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?,
                        )
                    };
                    let bonds = heap.slice_mut(all_backbone_bonds)?;
                    *bonds.get_mut(target).ok_or(SourceHeapError::PointerOutOfBounds)? =
                        source_pair.0;
                    *bonds
                        .get_mut(target + 1)
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = source_pair.1;
                    kk = kk.wrapping_add(1);
                }
                let end = pair_offset(number_of_backbone_bonds)?;
                let bonds = heap.slice_mut(all_backbone_bonds)?;
                *bonds.get_mut(end).ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
                *bonds
                    .get_mut(end + 1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
                number_of_backbone_bonds = number_of_backbone_bonds.wrapping_sub(1);
            }
            k = k.wrapping_sub(1);
        }

        let mut error_text = [0_i8; STR_ERR_LEN as usize];
        let validation_error = OAD_ValidatePolymerAndPseudoElementData(
            heap,
            original_atom_data,
            POLYMERS_MODERN as i32,
            1,
            Some(&mut error_text),
            0,
        )?;
        if validation_error != 0 {
            return Ok(ret);
        }

        if number_of_units > 0 {
            for unit_index in 0..number_of_units {
                let unit_index_usize = usize::try_from(unit_index)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let unit_pointer = *heap
                    .slice(polymer.units.as_const())?
                    .get(unit_index_usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let unit = heap
                    .slice(unit_pointer.as_const())?
                    .first()
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if unit.na < 2
                    || unit.nb < 2
                    || unit.cap1_is_undef == 0
                    || unit.cap2_is_undef == 0
                {
                    continue;
                }

                let backbone_count = usize::try_from(number_of_backbone_bonds)
                    .ok()
                    .and_then(|value| value.checked_mul(2))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let backbone_values = heap
                    .slice(all_backbone_bonds.as_const())?
                    .get(..backbone_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                let error = analyze_CRU_folding(
                    heap,
                    original_atom_data,
                    unit_index,
                    number_of_backbone_bonds,
                    &backbone_values,
                    number_of_xclasses,
                    &xclass_by_original,
                    edits,
                )?;
                if error != 0 {
                    ret = (_IS_WARNING as i32).max(error);
                }
            }
        }
        Ok(ret)
    })();

    let free_original = if original_numbers.is_null() {
        Ok(())
    } else {
        inchi_free(heap, original_numbers)
    };
    let free_backbone = if all_backbone_bonds.is_null() {
        Ok(())
    } else {
        inchi_free(heap, all_backbone_bonds)
    };
    match execution {
        Err(error) => Err(error),
        Ok(ret) => {
            free_original?;
            free_backbone?;
            Ok(ret)
        }
    }
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn analyze_CRU_folding(
    heap: &mut SourceHeap,
    original_atom_data: &ORIG_ATOM_DATA,
    unit_index: i32,
    all_backbone_bond_count: i32,
    all_backbone_bonds: &[i32],
    xclass_count: i32,
    xclasses: &[i32],
    edits: &OAD_StructureEdits,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2146 analyze_CRU_folding
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int analyze_CRU_folding(ORIG_ATOM_DATA *orig_at_data,
                        int iunit,
                        int n_all_bkb,
                        int *all_bkb,
                        int nxclasses, 
                        int *xc,
                        OAD_StructureEdits *ed)
{
    int ret = _IS_OKAY;
    int err, i, j, k, m, fail, a1, a2;
    int n_cuts = 0, n_frags = 0; 
    int n_frags_in_repeating_subunit = 0;
    int n_fold, n_frag_classes = 0;
    int subunit_last_atom, next_subunit_first_atom = 0;
    int *cut = NULL;        /* [ bkbond1at1, bkbond1at2,  bkbond2at1,bkbond2at2, ... ] 
                               these are (atoms of) backbone bonds which are non-cyclic and non-multiple ('breakable')      */
    DiylFrag **frag=NULL;   /* frag is divalent fragment surrounded by 'cut' bonds, so it may be a repeating CRU sub-unit   */
    int *frag_class=NULL;   /* fragments are classified, by their signatures, to produce unique labelling; 
                            if the two fragments have the same class, they have the same signature and whence are equivalent */
    int *frag_xc_counts = NULL; /* counts of xclass atoms in CRU, order of class numbers    */
    char pStrErr[STR_ERR_LEN];

    OAD_PolymerUnit *u = orig_at_data->polymer->units[iunit];
    ITRACE_("\n\n%-s\t\t%-s:%-d", "analyze_CRU_folding()", __FILE__,__LINE__);

    pStrErr[0] = '\0'; /* djb-rwth: fixing coverity ID #499611; pStrErr is a dummy parameter in this function and is never used */

    /* Reserve space for frag-specific xclass counts */
    frag_xc_counts = (int *)inchi_calloc((long long)nxclasses + 1, sizeof(int)); /* djb-rwth: cast operator added */
    if (!frag_xc_counts)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }


    /* Prepare list of cuts - backbone lying on the way from cap1 to cap2 */
    cut = (int *)inchi_calloc(2 * (long long)n_all_bkb, sizeof(int)); /* djb-rwth: cast operator added */
    if (!cut)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }

    OAD_PolymerUnit_DebugTrace(u);
    OAD_CollectBackboneBonds(orig_at_data,
                            u->na, u->alist,
                            u->end_atom1, u->end_atom2,
                            &(u->nbkbonds), u->bkbonds,
                            &err, pStrErr);
    if (err)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    OAD_PolymerUnit_DebugTrace(u);
    if (u->nbkbonds < 1)
    {
        goto exit_function;
    }

    /* Make 'cut' list from the bonds which are both in all_bkb and u->bkb
       (all_bkb eliminates bonds with order >1 and cyclic ones, 
       but may contain artificial cyclizing bond) 
    */
    for (i = 0; i <u->nbkbonds; i++)
    {
        a1 = u->bkbonds[i][0];
        a2 = u->bkbonds[i][1];
        for (j = 0; j < n_all_bkb; j++)
        {
            if (bIsSameBond(a1, a2, all_bkb[2 * j], all_bkb[2 * j + 1]))
            {
                cut[2 * n_cuts] = a1; /* djb-rwth: buffer overrun implicitly avoided in loop condition */
                cut[2 * n_cuts + 1] = a2;
                n_cuts++;
                break;
            }
        }
    }
    if (n_cuts < 1)
    {
        /* no valid sub-units is available */
        goto exit_function;
    }

    /* Collect fragments */
    n_frags = n_cuts + 1;
    frag = (DiylFrag**) inchi_calloc(n_frags, sizeof(DiylFrag *));
    if (!frag)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    frag_class = (int *) inchi_calloc(n_frags, sizeof(int));
    if (!frag_class)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    n_frag_classes = 0;
    for (i = 0; i < n_frags; i++)
    {
        /* Create fragment */
        int forbidden[4], novel=1;
        DiylFrag *pfrag = NULL; 

        /* Calculate and store signature of the fragment */
        /* 
            end_atom1...cut[i-1])---frag[i]---cut[i]---...end_atom2
        */
        if (i == 0)
        {
            a1 = u->end_atom1;
            forbidden[0] = u->cap1;
        }
        else
        {
            a1 = cut[2 * (i-1) + 1];            /* near end of prev cut */
            forbidden[0] = cut[2 * (i - 1) ];   /* far  end of prev cut */
        }
        forbidden[1] = a1;
        if (i==n_frags-1)
        {
            a2 = u->end_atom2;
            forbidden[2] = u->cap2;
        }
        else
        {
            a2 = cut[2 * i];                    /* near end of next cut */
            forbidden[2] = cut[2 * i + 1];      /* far  end of next cut */
        }
        forbidden[3] = a2;

        pfrag = DiylFrag_New(u->na, a1, a2, "");
        if (!pfrag)
        {
            ret = _IS_ERROR;
            goto exit_function;
        }
        frag[i] = pfrag;

        ret = OAD_CollectReachableAtoms(orig_at_data, a1, 2, forbidden,
                                  &pfrag->na, pfrag->alist, &err, pStrErr);
        if (ret==_IS_ERROR)
        {
            goto exit_function;
        }

        DiylFrag_MakeSignature(pfrag, nxclasses, xc, frag_xc_counts); 

        novel = 1;
        for (j = 0; j < i; j++)
        {
            if ( !DiylFrag_Diff(frag[i], frag[j]) )
            {
                frag_class[i] = frag_class[j];
                novel = 0;
                break;
            }
        }
        if (novel)
        {
            frag_class[i] = n_frag_classes++;
        }

        ITRACE_("\nCANDIDATE CRU SUBUNIT %-d/%-d (CLASS #%-d)\t", i+1, n_frags, frag_class[i]);
        DiylFrag_DebugTrace(pfrag);
    }

    if (n_frag_classes == n_frags)
    {
        /* All classes are distinct ==> no repeats, folding is impossible, skip the CRU */
        goto exit_function;
    }
        
    n_frags_in_repeating_subunit = len_repeating_subsequence(frag_class, NULL, n_frags);
    if (0 == n_frags_in_repeating_subunit)
    {
        /* valid repeating pattern not found */
        goto exit_function;
    }
    n_fold = n_frags / n_frags_in_repeating_subunit;
    if (1==n_fold || (0!=n_frags%n_frags_in_repeating_subunit) )    
    {
        /* valid repeating pattern not found */
        goto exit_function;
    }
    ITRACE_("\n");

    /* {1...2}---{5...6}---{8...9} */

    ITRACE_("\n* Found %-d times foldable unit of %-d fragments\n* First repeating sub-unit formed by %-d-fragment backbone : ",
            n_fold, n_frags, n_frags_in_repeating_subunit);
    
    for (k = 0; k < n_frags_in_repeating_subunit && n_frags_in_repeating_subunit < n_frags && frag[k]; k++) /* djb-rwth: fixing a NULL pointer dereference and buffer overflow */
    {
        if (frag[k]->end1 == frag[k]->end2)
        {
            ITRACE_("-{%-d}-", frag[k]->end1, frag[k]->end2);
        }
        else
        {
            ITRACE_("-{%-d...%-d}-", frag[k]->end1, frag[k]->end2);
        }
    }

    ITRACE_("\n");
    ITRACE_("* Backbone pattern for %-d fragments that may be removed :  ", n_frags - n_frags_in_repeating_subunit);
    for (k = n_frags_in_repeating_subunit; k < n_frags; k++)
    {
        if (frag[k]->end1 == frag[k]->end2)
        {
            ITRACE_("-{%-d}-", frag[k]->end1, frag[k]->end2);
        }
        else
        {
            ITRACE_("-{%-d...%-d}-", frag[k]->end1, frag[k]->end2);
        }
    }
    ITRACE_("\n");
    
    /* Folding is possible, prepare the edits 
        Keep the least in-CRU repeating subunit 
            { frag[0] ... frag[n_frags_in_repeating_subunit-1] }
        and remove 
            { frag[n_frags_in_repeating_subunit]...frag[n_frags-1] } and all side chain attached to that  
    
        NB: which bond is modified and which is broke is important for applying these edits further!			
    */

    /* Break bond from the subunit to the next fragment and replace an original 
       bond to "right" cap with bond from the subunit "right" atom  
    */

    /*djb-rwth: the whole block had to be rewritten to fix NULL pointer dereference */
    if (n_frags_in_repeating_subunit < n_frags && frag[n_frags_in_repeating_subunit] && frag[n_frags_in_repeating_subunit - 1]) /* djb-rwth: fixing a NULL pointer dereference and buffer overflow */
    {
        subunit_last_atom        = frag[n_frags_in_repeating_subunit - 1]->end2;
        next_subunit_first_atom  = frag[n_frags_in_repeating_subunit]->end1;

        fail = 0;
        fail += IntArray_Append(ed->del_bond, subunit_last_atom);
        fail += IntArray_Append(ed->del_bond, next_subunit_first_atom);

        fail += IntArray_Append(ed->mod_bond, u->end_atom2);
        fail += IntArray_Append(ed->mod_bond, u->cap2);
        fail += IntArray_Append(ed->mod_bond, subunit_last_atom);
        fail += IntArray_Append(ed->mod_bond, u->cap2);

        if (fail)
        {
            ret = _IS_ERROR;
            goto exit_function;
        }
    }
            
    /*	Now collect all backbone atoms to be deleted (we will then delete the
    associated side chains also, but no need to reveal them at the moment)	*/

    for (k = n_frags_in_repeating_subunit; k < n_frags; k++)
    {
        if (frag[k]) /* djb-rwth: fixing a NULL pointer dereference */
        {
            for (m = 0; m < frag[k]->na; m++)
            {
                fail = IntArray_AppendIfAbsent(ed->del_atom, frag[k]->alist[m]);
                if (fail)
                {
                    ret = _IS_ERROR;
                    goto exit_function;
                }
            }
        }
    }
    /* Care on atom coordinates: as bond to cap2 changes, 
       we use coordinates of next_subunit_first_atom for cap2 
    */
    fail = 0;
    fail += IntArray_Append(ed->mod_coord, next_subunit_first_atom);
    fail += IntArray_Append(ed->mod_coord, u->cap2);
    if (fail)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }

exit_function:
    if (cut)
    {
        inchi_free(cut);
    }
    if (frag)
    {
        for (i = 0; i < n_frags; i++)
        {
            DiylFrag_Free(frag[i]);
            inchi_free(frag[i]);
        }
        inchi_free(frag);
    }
    if (frag_class)
    {
        inchi_free(frag_class);
    }
    if (frag_xc_counts)
    {
        inchi_free(frag_xc_counts);
    }

    return ret;
}
    */
    // END INCHI C FUNCTION: analyze_CRU_folding
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: analyze_CRU_folding
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; sizeof(int)=4 and sizeof(void*)=8.
    // INCHI✔️❌: ITRACE_ is 0 && _inchi_trace; trace arguments are unevaluated, while surrounding loop tests and end1/end2 branches remain active.
    // INCHI✔️❌: All direct callees are completed source ports; allocation failures map to _IS_ERROR and every owned temporary follows source cleanup order.
    // INCHI✔️❌: Rust preserves partial edit-array mutations and all six append attempts before the first aggregate fail check.
    // INCHI✔️❌: SourceHeap pointer-map access, fragment clones, and mirror vectors add material overhead beyond native pointer indexing.
    // END INCHI ACTIVE MACRO CONFIGURATION: analyze_CRU_folding

    fn is_allocation_error(error: SourceHeapError) -> bool {
        matches!(
            error,
            SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange
                | SourceHeapError::AllocationFailed
        )
    }
    fn pointed_value<T: Clone + 'static>(
        heap: &SourceHeap,
        pointer: SourceMutPointer<T>,
    ) -> Result<T, SourceHeapError> {
        heap.slice(pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }
    fn append(
        heap: &mut SourceHeap,
        pointer: SourceMutPointer<INT_ARRAY>,
        value: i32,
    ) -> Result<i32, SourceHeapError> {
        let mut array = pointed_value(heap, pointer)?;
        let result = IntArray_Append(heap, Some(&mut array), value);
        heap.slice_mut(pointer)?[0] = array;
        result
    }
    fn append_if_absent(
        heap: &mut SourceHeap,
        pointer: SourceMutPointer<INT_ARRAY>,
        value: i32,
    ) -> Result<i32, SourceHeapError> {
        let mut array = pointed_value(heap, pointer)?;
        let result = IntArray_AppendIfAbsent(heap, &mut array, value);
        heap.slice_mut(pointer)?[0] = array;
        result
    }

    let polymer = pointed_value(heap, original_atom_data.polymer)?;
    let unit_offset = usize::try_from(unit_index)
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let unit_pointer = *heap
        .slice(polymer.units.as_const())?
        .get(unit_offset)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut unit = pointed_value(heap, unit_pointer)?;

    let mut cut = SourceMutPointer::null();
    let mut fragments: SourceMutPointer<SourceMutPointer<DiylFrag>> = SourceMutPointer::null();
    let mut fragment_classes = SourceMutPointer::null();
    let mut fragment_xclass_counts = SourceMutPointer::null();
    let mut empty_literal = SourceMutPointer::null();
    let mut fragment_count = 0_i32;

    let operation = (|| -> Result<i32, SourceHeapError> {
        let xclass_slots = i64::from(xclass_count).wrapping_add(1) as u64;
        fragment_xclass_counts = match inchi_calloc::<i32>(heap, xclass_slots, 4) {
            Ok(pointer) => pointer,
            Err(error) if is_allocation_error(error) => return Ok(_IS_ERROR as i32),
            Err(error) => return Err(error),
        };

        let cut_slots = i64::from(all_backbone_bond_count).wrapping_mul(2) as u64;
        cut = match inchi_calloc::<i32>(heap, cut_slots, 4) {
            Ok(pointer) => pointer,
            Err(error) if is_allocation_error(error) => return Ok(_IS_ERROR as i32),
            Err(error) => return Err(error),
        };

        OAD_PolymerUnit_DebugTrace(Some(&unit));
        let atom_count = usize::try_from(unit.na.max(0))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let atom_list = heap
            .slice(unit.alist.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec();
        let mut local_error = 0_i32;
        let mut error_text = [0_i8; STR_ERR_LEN as usize];
        let collect_result = OAD_CollectBackboneBonds(
            heap,
            original_atom_data,
            unit.na,
            &atom_list,
            unit.end_atom1,
            unit.end_atom2,
            &mut unit.nbkbonds,
            unit.bkbonds,
            &mut local_error,
            Some(&mut error_text),
        );
        heap.slice_mut(unit_pointer)?[0] = unit.clone();
        collect_result?;
        if local_error != 0 {
            return Ok(_IS_ERROR as i32);
        }
        OAD_PolymerUnit_DebugTrace(Some(&unit));
        if unit.nbkbonds < 1 {
            return Ok(_IS_OKAY as i32);
        }

        let mut cut_count = 0_i32;
        for backbone_index in 0..unit.nbkbonds {
            let row_pointer = *heap
                .slice(unit.bkbonds.as_const())?
                .get(usize::try_from(backbone_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let row = heap.slice(row_pointer.as_const())?;
            let atom1 = *row.first().ok_or(SourceHeapError::PointerOutOfBounds)?;
            let atom2 = *row.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?;
            for all_index in 0..all_backbone_bond_count.max(0) {
                let offset = usize::try_from(all_index)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                    .checked_mul(2)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                let all_atom1 = *all_backbone_bonds
                    .get(offset)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let all_atom2 = *all_backbone_bonds
                    .get(offset + 1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if bIsSameBond(atom1, atom2, all_atom1, all_atom2) != 0 {
                    let output = usize::try_from(cut_count)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                        .checked_mul(2)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                    *heap
                        .slice_mut(cut)?
                        .get_mut(output)
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = atom1;
                    *heap
                        .slice_mut(cut)?
                        .get_mut(output + 1)
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = atom2;
                    cut_count = cut_count
                        .checked_add(1)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                    break;
                }
            }
        }
        if cut_count < 1 {
            return Ok(_IS_OKAY as i32);
        }

        fragment_count = cut_count
            .checked_add(1)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        fragments = match inchi_calloc::<SourceMutPointer<DiylFrag>>(
            heap,
            fragment_count as u64,
            8,
        ) {
            Ok(pointer) => pointer,
            Err(error) if is_allocation_error(error) => return Ok(_IS_ERROR as i32),
            Err(error) => return Err(error),
        };
        fragment_classes = match inchi_calloc::<i32>(heap, fragment_count as u64, 4) {
            Ok(pointer) => pointer,
            Err(error) if is_allocation_error(error) => return Ok(_IS_ERROR as i32),
            Err(error) => return Err(error),
        };

        empty_literal = heap.allocate_model_storage(vec![0_i8])?;
        let mut fragment_class_count = 0_i32;
        let count_mirror_length = usize::try_from(xclass_slots)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        let mut count_mirror = Vec::new();
        count_mirror
            .try_reserve_exact(count_mirror_length)
            .map_err(|_| SourceHeapError::AllocationFailed)?;
        count_mirror.resize(count_mirror_length, 0_i32);
        for fragment_index in 0..fragment_count {
            let mut forbidden = [0_i32; 4];
            let atom1;
            if fragment_index == 0 {
                atom1 = unit.end_atom1;
                forbidden[0] = unit.cap1;
            } else {
                let previous = usize::try_from(fragment_index.wrapping_sub(1))
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                    .checked_mul(2)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                atom1 = heap.slice(cut.as_const())?[previous + 1];
                forbidden[0] = heap.slice(cut.as_const())?[previous];
            }
            forbidden[1] = atom1;
            let atom2;
            if fragment_index == fragment_count.wrapping_sub(1) {
                atom2 = unit.end_atom2;
                forbidden[2] = unit.cap2;
            } else {
                let next = usize::try_from(fragment_index)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                    .checked_mul(2)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                atom2 = heap.slice(cut.as_const())?[next];
                forbidden[2] = heap.slice(cut.as_const())?[next + 1];
            }
            forbidden[3] = atom2;

            let fragment_pointer =
                DiylFrag_New(heap, unit.na, atom1, atom2, empty_literal.as_const())?;
            if fragment_pointer.is_null() {
                return Ok(_IS_ERROR as i32);
            }
            heap.slice_mut(fragments)?
                [usize::try_from(fragment_index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?] = fragment_pointer;

            let forbidden_pointer = heap.allocate_model_storage(forbidden.to_vec())?;
            let mut fragment = pointed_value(heap, fragment_pointer)?;
            let reachable_result = OAD_CollectReachableAtoms(
                heap,
                original_atom_data,
                atom1,
                2,
                forbidden_pointer,
                &mut fragment.na,
                fragment.alist,
                &mut local_error,
                Some(&mut error_text),
            );
            heap.slice_mut(fragment_pointer)?[0] = fragment.clone();
            heap.free(forbidden_pointer)?;
            let reachable_status = reachable_result?;
            if reachable_status == _IS_ERROR as i32 {
                return Ok(_IS_ERROR as i32);
            }

            let signature_result = heap.with_slice_mut_and_heap_mut(
                fragment_pointer,
                |values, heap| {
                    let fragment = values
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    DiylFrag_MakeSignature(
                        heap,
                        fragment,
                        xclass_count,
                        xclasses,
                        &mut count_mirror,
                    )
                },
            );
            signature_result?;

            let mut novel = true;
            let current = pointed_value(heap, fragment_pointer)?;
            for previous_index in 0..fragment_index {
                let previous_pointer = heap.slice(fragments.as_const())?
                    [usize::try_from(previous_index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?];
                let previous = pointed_value(heap, previous_pointer)?;
                if DiylFrag_Diff(heap, &current, &previous)? == 0 {
                    let class = heap.slice(fragment_classes.as_const())?
                        [usize::try_from(previous_index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?];
                    heap.slice_mut(fragment_classes)?
                        [usize::try_from(fragment_index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?] = class;
                    novel = false;
                    break;
                }
            }
            if novel {
                heap.slice_mut(fragment_classes)?
                    [usize::try_from(fragment_index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?] = fragment_class_count;
                fragment_class_count = fragment_class_count
                    .checked_add(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            }
            DiylFrag_DebugTrace(Some(&current))?;
        }
        if fragment_class_count == fragment_count {
            return Ok(_IS_OKAY as i32);
        }
        let classes = heap.slice(fragment_classes.as_const())?;
        let repeating_fragment_count = len_repeating_subsequence(
            Some(classes),
            None,
            fragment_count,
        )?;
        if repeating_fragment_count == 0 {
            return Ok(_IS_OKAY as i32);
        }
        let fold_count = fragment_count / repeating_fragment_count;
        if fold_count == 1 || fragment_count % repeating_fragment_count != 0 {
            return Ok(_IS_OKAY as i32);
        }

        for index in 0..repeating_fragment_count {
            if repeating_fragment_count >= fragment_count {
                break;
            }
            let pointer = heap.slice(fragments.as_const())?
                [usize::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?];
            if pointer.is_null() {
                break;
            }
            let fragment = pointed_value(heap, pointer)?;
            let _same_end = fragment.end1 == fragment.end2;
        }
        for index in repeating_fragment_count..fragment_count {
            let pointer = heap.slice(fragments.as_const())?
                [usize::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?];
            let fragment = pointed_value(heap, pointer)?;
            let _same_end = fragment.end1 == fragment.end2;
        }

        let next_index = usize::try_from(repeating_fragment_count)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let previous_index = usize::try_from(repeating_fragment_count.wrapping_sub(1))
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let next_pointer = *heap
            .slice(fragments.as_const())?
            .get(next_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let previous_pointer = *heap
            .slice(fragments.as_const())?
            .get(previous_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut next_subunit_first_atom = 0_i32;
        if repeating_fragment_count < fragment_count
            && !next_pointer.is_null()
            && !previous_pointer.is_null()
        {
            let previous = pointed_value(heap, previous_pointer)?;
            let next = pointed_value(heap, next_pointer)?;
            let subunit_last_atom = previous.end2;
            next_subunit_first_atom = next.end1;
            let mut fail = 0_i32;
            fail = fail.wrapping_add(append(heap, edits.del_bond, subunit_last_atom)?);
            fail = fail.wrapping_add(append(heap, edits.del_bond, next_subunit_first_atom)?);
            fail = fail.wrapping_add(append(heap, edits.mod_bond, unit.end_atom2)?);
            fail = fail.wrapping_add(append(heap, edits.mod_bond, unit.cap2)?);
            fail = fail.wrapping_add(append(heap, edits.mod_bond, subunit_last_atom)?);
            fail = fail.wrapping_add(append(heap, edits.mod_bond, unit.cap2)?);
            if fail != 0 {
                return Ok(_IS_ERROR as i32);
            }
        }

        for fragment_index in repeating_fragment_count..fragment_count {
            let pointer = heap.slice(fragments.as_const())?
                [usize::try_from(fragment_index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?];
            if !pointer.is_null() {
                let fragment = pointed_value(heap, pointer)?;
                for atom_index in 0..fragment.na {
                    let atom = *heap
                        .slice(fragment.alist.as_const())?
                        .get(usize::try_from(atom_index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if append_if_absent(heap, edits.del_atom, atom)? != 0 {
                        return Ok(_IS_ERROR as i32);
                    }
                }
            }
        }
        let mut fail = 0_i32;
        fail = fail.wrapping_add(append(heap, edits.mod_coord, next_subunit_first_atom)?);
        fail = fail.wrapping_add(append(heap, edits.mod_coord, unit.cap2)?);
        if fail != 0 {
            return Ok(_IS_ERROR as i32);
        }
        Ok(_IS_OKAY as i32)
    })();

    let empty_cleanup = heap.free(empty_literal);
    let cut_cleanup = inchi_free(heap, cut);
    let mut fragment_cleanup = Ok(());
    if !fragments.is_null() {
        let pointers = match heap.slice(fragments.as_const()) {
            Ok(values) => values
                .get(..usize::try_from(fragment_count.max(0)).unwrap_or(0))
                .unwrap_or(values)
                .to_vec(),
            Err(error) => {
                fragment_cleanup = Err(error);
                Vec::new()
            }
        };
        for pointer in pointers {
            if pointer.is_null() {
                continue;
            }
            let close = heap.with_slice_mut_and_heap_mut(pointer, |values, heap| {
                DiylFrag_Free(heap, values.first_mut())
            });
            if fragment_cleanup.is_ok() {
                fragment_cleanup = close;
            }
            let free = inchi_free(heap, pointer);
            if fragment_cleanup.is_ok() {
                fragment_cleanup = free;
            }
        }
        let free = inchi_free(heap, fragments);
        if fragment_cleanup.is_ok() {
            fragment_cleanup = free;
        }
    }
    let class_cleanup = inchi_free(heap, fragment_classes);
    let count_cleanup = inchi_free(heap, fragment_xclass_counts);

    match operation {
        Err(error) => Err(error),
        Ok(status) => {
            empty_cleanup?;
            cut_cleanup?;
            fragment_cleanup?;
            class_cleanup?;
            count_cleanup?;
            Ok(status)
        }
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn count_colors_in_sequence(
    color: &[i32],
    n: i32,
    maxcol: i32,
    counts: &mut [i32],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2463 count_colors_in_sequence
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
int count_colors_in_sequence( int *color, int n, int maxcol, int *counts)
{
    int i, ncol=0;
    memset(counts, 0, maxcol * sizeof(int)); /* djb-rwth: memset_s C11/Annex K variant? */
    for (i = 0; i<n; i++) 
    {
        int colori = color[i];
        if (colori < 0) /* removed orig atom (H D etc.) */
        {
            continue;
        }
        if (0==counts[colori ])
        {
            ncol++;
        }
        counts[ color[i] ]++;
    }
    return ncol;
}
    */
    // END INCHI C FUNCTION: count_colors_in_sequence
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: count_colors_in_sequence
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; sizeof(int)=4.
    // INCHI✔️✔️: Direct indexed loops preserve O(maxcol+n) time, O(1) extra space, mutation order, and cache locality.
    // END INCHI ACTIVE MACRO CONFIGURATION: count_colors_in_sequence

    let maximum_colors = usize::try_from(maxcol)
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    counts
        .get_mut(..maximum_colors)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .fill(0);

    let mut number_of_colors = 0_i32;
    for index in 0..n {
        let index = usize::try_from(index)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let color_index = *color
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if color_index < 0 {
            continue;
        }
        let color_index = usize::try_from(color_index)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let count = counts
            .get_mut(color_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if *count == 0 {
            number_of_colors = number_of_colors.wrapping_add(1);
        }
        *count = count.wrapping_add(1);
    }
    Ok(number_of_colors)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn len_repeating_subsequence(
    color: Option<&[i32]>,
    color2: Option<&[i32]>,
    n: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2489 len_repeating_subsequence
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
int len_repeating_subsequence(int *color, int *color2, int n)
{
    int m, k;

    if (n < 2 || !color)
    {
        return 0;
    }

    for (m = 0; m < (n + 1) / 2; m++)
    {
        for (k = m + 1; k < n; k++)
        {
            if (color[k] != color[k - m - 1]) 
            { 
                goto nextm; 
            }
            if (color2 && color2[k] != color2[k - m - 1])
            {
                goto nextm;
            }
        }
        return (m + 1);
nextm:	;
    }

    return 0;
}
    */
    // END INCHI C FUNCTION: len_repeating_subsequence
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: len_repeating_subsequence
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; sizeof(int)=4.
    // INCHI✔️✔️: Candidate periods, nested-loop order, partial-final-repeat acceptance, optional color2 short-circuit, and O(n^2) worst-case cost are unchanged.
    // INCHI✔️✔️: Source signed-overflow and invalid-pointer states return structured errors without changing any defined result.
    // END INCHI ACTIVE MACRO CONFIGURATION: len_repeating_subsequence

    if n < 2 || color.is_none() {
        return Ok(0);
    }
    let color = color.expect("checked above");
    let candidate_count = n
        .checked_add(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?
        / 2;
    for m in 0..candidate_count {
        let mut matches = true;
        let period = m
            .checked_add(1)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        for k in period..n {
            let current = usize::try_from(k)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let previous = usize::try_from(k.wrapping_sub(period))
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            if color
                .get(current)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                != color
                    .get(previous)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
            {
                matches = false;
                break;
            }
            if let Some(color2) = color2
                && color2
                    .get(current)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    != color2
                        .get(previous)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
            {
                matches = false;
                break;
            }
        }
        if matches {
            return Ok(period);
        }
    }
    Ok(0)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn OAD_StructureEdits_Clear(
    heap: &mut SourceHeap,
    edits: &mut OAD_StructureEdits,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1769 OAD_StructureEdits_Clear
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
void OAD_StructureEdits_Clear(OAD_StructureEdits *ed)
{
    if (ed->del_atom)
    {
        IntArray_Free(ed->del_atom);
        inchi_free(ed->del_atom);
        ed->del_atom = NULL;
    }
    if (ed->del_bond)
    {
        IntArray_Free(ed->del_bond);
        inchi_free(ed->del_bond);
        ed->del_bond = NULL;
    }
    if (ed->mod_bond)
    {
        IntArray_Free(ed->mod_bond);
        inchi_free(ed->mod_bond);
        ed->mod_bond = NULL;
    }
    if (ed->new_bond)
    {
        IntArray_Free(ed->new_bond);
        inchi_free(ed->new_bond);
        ed->new_bond = NULL;
    }
    if (ed->mod_coord)
    {
        IntArray_Free(ed->mod_coord);
        inchi_free(ed->mod_coord);
        ed->mod_coord = NULL;
    }

    return;
}
    */
    // END INCHI C FUNCTION: OAD_StructureEdits_Clear
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OAD_StructureEdits_Clear
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64.
    // INCHI✔️❌: inchi_free is the GCC/Linux macro expansion to free non-null pointers.
    // INCHI✔️❌: SourceHeap allocation-map lookups add overhead absent from direct C frees.
    // END INCHI ACTIVE MACRO CONFIGURATION: OAD_StructureEdits_Clear

    if !edits.del_atom.is_null() {
        heap.with_slice_mut_and_heap_mut(edits.del_atom, |arrays, heap| {
            IntArray_Free(heap, arrays.first_mut())
        })?;
        inchi_free(heap, edits.del_atom)?;
        edits.del_atom = SourceMutPointer::null();
    }
    if !edits.del_bond.is_null() {
        heap.with_slice_mut_and_heap_mut(edits.del_bond, |arrays, heap| {
            IntArray_Free(heap, arrays.first_mut())
        })?;
        inchi_free(heap, edits.del_bond)?;
        edits.del_bond = SourceMutPointer::null();
    }
    if !edits.mod_bond.is_null() {
        heap.with_slice_mut_and_heap_mut(edits.mod_bond, |arrays, heap| {
            IntArray_Free(heap, arrays.first_mut())
        })?;
        inchi_free(heap, edits.mod_bond)?;
        edits.mod_bond = SourceMutPointer::null();
    }
    if !edits.new_bond.is_null() {
        heap.with_slice_mut_and_heap_mut(edits.new_bond, |arrays, heap| {
            IntArray_Free(heap, arrays.first_mut())
        })?;
        inchi_free(heap, edits.new_bond)?;
        edits.new_bond = SourceMutPointer::null();
    }
    if !edits.mod_coord.is_null() {
        heap.with_slice_mut_and_heap_mut(edits.mod_coord, |arrays, heap| {
            IntArray_Free(heap, arrays.first_mut())
        })?;
        inchi_free(heap, edits.mod_coord)?;
        edits.mod_coord = SourceMutPointer::null();
    }

    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn OAD_StructureEdits_Init(
    heap: &mut SourceHeap,
    edits: &mut OAD_StructureEdits,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1735 OAD_StructureEdits_Init
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int OAD_StructureEdits_Init(OAD_StructureEdits *ed)
{
    ed->del_side_chains = 0; /* by default, do not delete */

    ed->del_atom = (INT_ARRAY *)inchi_calloc(1, sizeof(INT_ARRAY));
    if (!ed->del_atom)						goto exitf;
    if (0 != IntArray_Alloc(ed->del_atom, 2))	goto exitf;

    ed->del_bond = (INT_ARRAY *)inchi_calloc(1, sizeof(INT_ARRAY));
    if (!ed->del_bond)						goto exitf;
    if (0 != IntArray_Alloc(ed->del_bond, 2))	goto exitf;

    ed->new_bond = (INT_ARRAY *)inchi_calloc(1, sizeof(INT_ARRAY));
    if (!ed->new_bond)						goto exitf;
    if (0 != IntArray_Alloc(ed->new_bond, 2))	goto exitf;

    ed->mod_bond = (INT_ARRAY *)inchi_calloc(1, sizeof(INT_ARRAY));
    if (!ed->mod_bond)						goto exitf;
    if (0 != IntArray_Alloc(ed->mod_bond, 12))	goto exitf;

    ed->mod_coord = (INT_ARRAY *)inchi_calloc(1, sizeof(INT_ARRAY));
    if (!ed->mod_coord)						goto exitf;
    if (0 != IntArray_Alloc(ed->mod_coord, 4))	goto exitf;


    return 0;

exitf:
    OAD_StructureEdits_Clear(ed);
    return _IS_ERROR;
}
    */
    // END INCHI C FUNCTION: OAD_StructureEdits_Init
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OAD_StructureEdits_Init
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; sizeof(INT_ARRAY)=24.
    // INCHI✔️❌: inchi_calloc is the active calloc macro and allocation failures become null.
    // INCHI✔️❌: SourceHeap allocation-map operations add overhead absent from direct C allocation.
    // END INCHI ACTIVE MACRO CONFIGURATION: OAD_StructureEdits_Init

    edits.del_side_chains = 0;
    let execution = (|| -> Result<i32, SourceHeapError> {
        edits.del_atom = match inchi_calloc(heap, 1, std::mem::size_of::<crate::source_types::INT_ARRAY>() as u64) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        if edits.del_atom.is_null() {
            return Ok(_IS_ERROR as i32);
        }
        if heap.with_slice_mut_and_heap_mut(edits.del_atom, |arrays, heap| {
            let array = arrays.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)?;
            IntArray_Alloc(heap, array, 2)
        })? != 0 {
            return Ok(_IS_ERROR as i32);
        }

        edits.del_bond = match inchi_calloc(heap, 1, std::mem::size_of::<crate::source_types::INT_ARRAY>() as u64) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        if edits.del_bond.is_null() {
            return Ok(_IS_ERROR as i32);
        }
        if heap.with_slice_mut_and_heap_mut(edits.del_bond, |arrays, heap| {
            let array = arrays.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)?;
            IntArray_Alloc(heap, array, 2)
        })? != 0 {
            return Ok(_IS_ERROR as i32);
        }

        edits.new_bond = match inchi_calloc(heap, 1, std::mem::size_of::<crate::source_types::INT_ARRAY>() as u64) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        if edits.new_bond.is_null() {
            return Ok(_IS_ERROR as i32);
        }
        if heap.with_slice_mut_and_heap_mut(edits.new_bond, |arrays, heap| {
            let array = arrays.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)?;
            IntArray_Alloc(heap, array, 2)
        })? != 0 {
            return Ok(_IS_ERROR as i32);
        }

        edits.mod_bond = match inchi_calloc(heap, 1, std::mem::size_of::<crate::source_types::INT_ARRAY>() as u64) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        if edits.mod_bond.is_null() {
            return Ok(_IS_ERROR as i32);
        }
        if heap.with_slice_mut_and_heap_mut(edits.mod_bond, |arrays, heap| {
            let array = arrays.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)?;
            IntArray_Alloc(heap, array, 12)
        })? != 0 {
            return Ok(_IS_ERROR as i32);
        }

        edits.mod_coord = match inchi_calloc(heap, 1, std::mem::size_of::<crate::source_types::INT_ARRAY>() as u64) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        if edits.mod_coord.is_null() {
            return Ok(_IS_ERROR as i32);
        }
        if heap.with_slice_mut_and_heap_mut(edits.mod_coord, |arrays, heap| {
            let array = arrays.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)?;
            IntArray_Alloc(heap, array, 4)
        })? != 0 {
            return Ok(_IS_ERROR as i32);
        }

        Ok(0)
    })();

    match execution {
        Ok(0) => Ok(0),
        Ok(_) => {
            OAD_StructureEdits_Clear(heap, edits)?;
            Ok(_IS_ERROR as i32)
        }
        Err(error) => {
            let _ = OAD_StructureEdits_Clear(heap, edits);
            Err(error)
        }
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn OAD_StructureEdits_DebugPrint(edits: &OAD_StructureEdits) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1807 OAD_StructureEdits_DebugPrint
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
void OAD_StructureEdits_DebugPrint(OAD_StructureEdits *ed)
{
    ITRACE_("\n*****************************\nOAD_StructureEdits @ %-p\n*****************************", ed); 
    ITRACE_("\nDel_side_chains :\t%-d\n", ed->del_side_chains); 
    ITRACE_("Del_atom:\t%-s", ed->del_atom->used ? "" : "(empty)\n");
    IntArray_DebugPrint(ed->del_atom);
    ITRACE_("Del_bond:\t%-s", ed->del_bond->used ? "" : "(empty)\n");
    IntArray_DebugPrint(ed->del_bond);
    ITRACE_("New_bond:\t%-s", ed->new_bond->used ? "" : "(empty)\n");
    IntArray_DebugPrint(ed->new_bond);
    ITRACE_("Mod_bond:\t%-s", ed->mod_bond->used ? "" : "(empty)\n");
    IntArray_DebugPrint(ed->mod_bond);
    ITRACE_("Mod_coord:\t%-s", ed->mod_coord->used ? "" : "(empty)\n");
    IntArray_DebugPrint(ed->mod_coord);
    
}
    */
    // END INCHI C FUNCTION: OAD_StructureEdits_DebugPrint
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OAD_StructureEdits_DebugPrint
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux, so ITRACE_ is 0 && _inchi_trace.
    // INCHI✔️✔️: ITRACE_ arguments are not evaluated; five active calls reach the O(1) no-op callee.
    // INCHI✔️✔️: Pointer values and call order are preserved without allocation or output.
    // END INCHI ACTIVE MACRO CONFIGURATION: OAD_StructureEdits_DebugPrint

    IntArray_DebugPrint(edits.del_atom);
    IntArray_DebugPrint(edits.del_bond);
    IntArray_DebugPrint(edits.new_bond);
    IntArray_DebugPrint(edits.mod_bond);
    IntArray_DebugPrint(edits.mod_coord);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{
        BOND_DOUBLE, BOND_SINGLE, CT_OUT_OF_RAM, INCHI_IOS_STRING, INCHI_IOS_TYPE_FILE, INChI, INChI_Aux, INT_ARRAY,
        OAD_Polymer, OAD_PolymerUnit, SourceFile, inp_ATOM,
    };

    #[test]
    fn source_port__runichi2__extract_orig_nums_from_auxinfo_string__line_1389() {
        let mut heap = SourceHeap::default();
        let input = c_string(&mut heap, b"AuxInfo=1/0/N:3,1;2/E:(1,2)");
        let mut orig = [91_i32; 6];
        assert_eq!(
            extract_orig_nums_from_auxinfo_string(&mut heap, input.as_const(), &mut orig),
            Ok(0)
        );
        assert_eq!(orig, [91, 3, 1, 2, 91, 91]);

        let repeated = c_string(&mut heap, b"/N:1,,2");
        let mut repeated_orig = [8_i32; 4];
        assert_eq!(
            extract_orig_nums_from_auxinfo_string(&mut heap, repeated.as_const(), &mut repeated_orig,),
            Ok(0)
        );
        assert_eq!(repeated_orig, [8, 1, 8, 8]);

        let zero = c_string(&mut heap, b"/N:0,7");
        let mut zero_orig = [6_i32; 3];
        assert_eq!(
            extract_orig_nums_from_auxinfo_string(&mut heap, zero.as_const(), &mut zero_orig),
            Ok(0)
        );
        assert_eq!(zero_orig, [6, 6, 6]);

        let narrowing = c_string(&mut heap, b"/N:2147483648,4294967296,7");
        let mut narrowed = [5_i32; 5];
        assert_eq!(
            extract_orig_nums_from_auxinfo_string(&mut heap, narrowing.as_const(), &mut narrowed,),
            Ok(0)
        );
        assert_eq!(narrowed, [5, i32::MIN, 5, 5, 5]);

        let overflow = c_string(&mut heap, b"/N:9223372036854775808/");
        let mut overflow_orig = [4_i32; 3];
        heap.set_source_errno(0);
        assert_eq!(
            extract_orig_nums_from_auxinfo_string(&mut heap, overflow.as_const(), &mut overflow_orig,),
            Ok(0)
        );
        assert_eq!(overflow_orig, [4, -1, 4]);
        assert_eq!(heap.source_errno(), 34);

        for invalid in [b"" as &[u8], b"/E:1", b"/N:", b"/N:-1", b"/N:A"] {
            let input = c_string(&mut heap, invalid);
            let mut untouched = [3_i32; 3];
            assert_eq!(
                extract_orig_nums_from_auxinfo_string(&mut heap, input.as_const(), &mut untouched,),
                Ok(_IS_ERROR as i32),
                "invalid={invalid:?}"
            );
            assert_eq!(untouched, [3, 3, 3]);
        }

        let high_byte = heap
            .allocate_model_storage(vec![b'/' as i8, b'N' as i8, b':' as i8, -1, 0])
            .unwrap();
        assert_eq!(
            extract_orig_nums_from_auxinfo_string(&mut heap, high_byte.as_const(), &mut [2_i32; 2],),
            Ok(_IS_ERROR as i32)
        );

        let short_input = c_string(&mut heap, b"/N:1,2");
        let mut short = [9_i32; 2];
        assert_eq!(
            extract_orig_nums_from_auxinfo_string(&mut heap, short_input.as_const(), &mut short,),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(short, [9, 1]);
    }

    #[test]
    fn source_port__runichi2__extract_nonstereo_eq_classes_from_auxinfo_string__line_1427() {
        let mut heap = SourceHeap::default();
        let grouped = c_string(&mut heap, b"AuxInfo=1/0/N:4,2,1,3/E:(1,3)(2,4)/rA:");
        let orig = [99, 4, 2, 1, 3];
        let mut nclasses = -7;
        let mut eclass = [77; 7];
        let mut by_orig = [88; 7];
        assert_eq!(
            extract_nonstereo_eq_classes_from_auxinfo_string(
                &mut heap,
                grouped.as_const(),
                4,
                &orig,
                &mut nclasses,
                &mut eclass,
                &mut by_orig,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(nclasses, 2);
        assert_eq!(eclass, [-1, 0, 1, 0, 1, 77, 77]);
        assert_eq!(by_orig, [-1, 0, 1, 1, 0, 88, 88]);

        let absent = c_string(&mut heap, b"AuxInfo=1/0/N:2,1/rA:");
        let mut absent_count = 9;
        let mut absent_classes = [5; 5];
        let mut absent_by_orig = [6; 5];
        assert_eq!(
            extract_nonstereo_eq_classes_from_auxinfo_string(
                &mut heap,
                absent.as_const(),
                2,
                &[0, 2, 1],
                &mut absent_count,
                &mut absent_classes,
                &mut absent_by_orig,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(absent_count, 0);
        assert_eq!(absent_classes, [-1, -1, -1, 5, 5]);
        assert_eq!(absent_by_orig, [-1, -1, -1, 6, 6]);

        let empty = c_string(&mut heap, b"/E:/x");
        let mut empty_count = -1;
        let mut empty_classes = [4; 4];
        let mut empty_by_orig = [3; 4];
        assert_eq!(
            extract_nonstereo_eq_classes_from_auxinfo_string(
                &mut heap,
                empty.as_const(),
                2,
                &[0, 2, 1],
                &mut empty_count,
                &mut empty_classes,
                &mut empty_by_orig,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(empty_count, 2);
        assert_eq!(empty_classes, [-1, 1, 2, 4]);
        assert_eq!(empty_by_orig, [-1, 2, 1, 3]);

        let slash = c_string(&mut heap, b"/E:(1/");
        let mut slash_count = -1;
        let mut slash_classes = [4; 4];
        let mut slash_by_orig = [3; 4];
        assert_eq!(
            extract_nonstereo_eq_classes_from_auxinfo_string(
                &mut heap,
                slash.as_const(),
                2,
                &[0, 1, 2],
                &mut slash_count,
                &mut slash_classes,
                &mut slash_by_orig,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(slash_count, 2);
        assert_eq!(slash_classes, [-1, 1, 2, 4]);
        assert_eq!(slash_by_orig, [-1, 1, 2, 3]);

        let unique = c_string(&mut heap, b"/E:(1,3)/x");
        let mut unique_count = -1;
        let mut unique_classes = [0; 5];
        let mut unique_by_orig = [0; 5];
        assert_eq!(
            extract_nonstereo_eq_classes_from_auxinfo_string(
                &mut heap,
                unique.as_const(),
                4,
                &[0, 1, 2, 3, 4],
                &mut unique_count,
                &mut unique_classes,
                &mut unique_by_orig,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(unique_count, 3);
        assert_eq!(unique_classes, [-1, 0, 2, 0, 3]);
        assert_eq!(unique_by_orig, unique_classes);

        let malformed = c_string(&mut heap, b"/E:(1,2X/");
        let mut malformed_count = 7;
        let mut malformed_classes = [4; 4];
        let mut malformed_by_orig = [3; 4];
        assert_eq!(
            extract_nonstereo_eq_classes_from_auxinfo_string(
                &mut heap,
                malformed.as_const(),
                2,
                &[0, 1, 2],
                &mut malformed_count,
                &mut malformed_classes,
                &mut malformed_by_orig,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(malformed_count, 0);
        assert_eq!(malformed_classes, [-1, 0, -1, 4]);
        assert_eq!(malformed_by_orig, [-1, -1, -1, 3]);

        let narrowed = c_string(&mut heap, b"/E:(-1,X");
        let mut narrowed_count = 8;
        let mut narrowed_classes = vec![7; 65_536];
        let mut narrowed_by_orig = [6; 2];
        assert_eq!(
            extract_nonstereo_eq_classes_from_auxinfo_string(
                &mut heap,
                narrowed.as_const(),
                1,
                &[0, 1],
                &mut narrowed_count,
                &mut narrowed_classes,
                &mut narrowed_by_orig,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(narrowed_count, 1);
        assert_eq!(narrowed_classes[0], -1);
        assert_eq!(narrowed_classes[1], 1);
        assert_eq!(narrowed_classes[65_535], 0);
        assert_eq!(narrowed_classes[65_534], 7);
        assert_eq!(narrowed_by_orig, [-1, 1]);

        let narrowed_zero = c_string(&mut heap, b"/E:(65536,1)/x");
        let mut zero_count = 8;
        let mut zero_classes = [7; 3];
        let mut zero_by_orig = [6; 3];
        assert_eq!(
            extract_nonstereo_eq_classes_from_auxinfo_string(
                &mut heap,
                narrowed_zero.as_const(),
                2,
                &[0, 1, 2],
                &mut zero_count,
                &mut zero_classes,
                &mut zero_by_orig,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(zero_count, 2);
        assert_eq!(zero_classes, [-1, 1, 2]);
        assert_eq!(zero_by_orig, [-1, 1, 2]);

        let overflow = c_string(&mut heap, b"/E:(9223372036854775808,X");
        let mut overflow_count = 8;
        let mut overflow_classes = vec![7; 65_536];
        let mut overflow_by_orig = [6; 2];
        heap.set_source_errno(0);
        assert_eq!(
            extract_nonstereo_eq_classes_from_auxinfo_string(
                &mut heap,
                overflow.as_const(),
                1,
                &[0, 1],
                &mut overflow_count,
                &mut overflow_classes,
                &mut overflow_by_orig,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(overflow_classes[65_535], 0);
        assert_eq!(heap.source_errno(), 34);

        let bad_index = c_string(&mut heap, b"/E:(3,1)/x");
        let mut bad_index_count = 8;
        let mut bad_index_classes = [7; 3];
        let mut bad_index_by_orig = [6; 3];
        assert_eq!(
            extract_nonstereo_eq_classes_from_auxinfo_string(
                &mut heap,
                bad_index.as_const(),
                2,
                &[0, 1, 2],
                &mut bad_index_count,
                &mut bad_index_classes,
                &mut bad_index_by_orig,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(bad_index_count, 0);
        assert_eq!(bad_index_classes, [-1, -1, -1]);
        assert_eq!(bad_index_by_orig, [-1, -1, -1]);

        let bad_orig = c_string(&mut heap, b"/E:/x");
        let mut bad_orig_count = 8;
        let mut bad_orig_classes = [7; 3];
        let mut bad_orig_by_orig = [6; 3];
        assert_eq!(
            extract_nonstereo_eq_classes_from_auxinfo_string(
                &mut heap,
                bad_orig.as_const(),
                2,
                &[0, 1, 9],
                &mut bad_orig_count,
                &mut bad_orig_classes,
                &mut bad_orig_by_orig,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(bad_orig_count, 2);
        assert_eq!(bad_orig_classes, [-1, 1, 2]);
        assert_eq!(bad_orig_by_orig, [-1, 1, -1]);

        let no_nul = heap
            .allocate_model_storage(b"no AuxInfo terminator".iter().map(|byte| *byte as i8).collect())
            .unwrap();
        let mut no_nul_count = 8;
        let mut no_nul_classes = [7; 2];
        let mut no_nul_by_orig = [6; 2];
        assert_eq!(
            extract_nonstereo_eq_classes_from_auxinfo_string(
                &mut heap,
                no_nul.as_const(),
                1,
                &[0, 1],
                &mut no_nul_count,
                &mut no_nul_classes,
                &mut no_nul_by_orig,
            ),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(no_nul_count, 0);
        assert_eq!(no_nul_classes, [-1, -1]);
        assert_eq!(no_nul_by_orig, [-1, -1]);

        let mut negative_count = 8;
        assert_eq!(
            extract_nonstereo_eq_classes_from_auxinfo_string(
                &mut heap,
                absent.as_const(),
                -1,
                &[],
                &mut negative_count,
                &mut [],
                &mut [],
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(negative_count, 0);

        let mut invalid_negative_count = 8;
        let mut invalid_negative_classes = [7];
        let mut invalid_negative_by_orig = [6];
        assert_eq!(
            extract_nonstereo_eq_classes_from_auxinfo_string(
                &mut heap,
                absent.as_const(),
                -2,
                &[],
                &mut invalid_negative_count,
                &mut invalid_negative_classes,
                &mut invalid_negative_by_orig,
            ),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(invalid_negative_count, 0);
        assert_eq!(invalid_negative_classes, [7]);
        assert_eq!(invalid_negative_by_orig, [6]);

        let mut short_count = 8;
        let mut short_classes = [7; 2];
        let mut short_by_orig = [6; 3];
        assert_eq!(
            extract_nonstereo_eq_classes_from_auxinfo_string(
                &mut heap,
                absent.as_const(),
                2,
                &[0, 1, 2],
                &mut short_count,
                &mut short_classes,
                &mut short_by_orig,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(short_count, 0);
        assert_eq!(short_classes, [7, 7]);
        assert_eq!(short_by_orig, [6, 6, 6]);

        let mut second_short_count = 8;
        let mut second_short_classes = [7; 3];
        let mut second_short_by_orig = [6; 2];
        assert_eq!(
            extract_nonstereo_eq_classes_from_auxinfo_string(
                &mut heap,
                absent.as_const(),
                2,
                &[0, 1, 2],
                &mut second_short_count,
                &mut second_short_classes,
                &mut second_short_by_orig,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(second_short_count, 0);
        assert_eq!(second_short_classes, [-1, -1, -1]);
        assert_eq!(second_short_by_orig, [6, 6]);
    }

    #[test]
    fn source_port__runichi2__diylfrag_free__line_2049() {
        let mut heap = SourceHeap::default();
        assert_eq!(DiylFrag_Free(&mut heap, None), Ok(()));
        assert_eq!(heap.live_allocation_count(), 0);

        let alist = heap.allocate_model_storage(vec![4_i32, 8, 15]).unwrap();
        let xclist = heap.allocate_model_storage(vec![16_i32, 23, 42]).unwrap();
        let signature = heap
            .allocate_model_storage(vec![b's' as i8, b'i' as i8, b'g' as i8, 0])
            .unwrap();
        let mut fragment = DiylFrag {
            na: 3,
            nb: 2,
            end1: 11,
            end2: 19,
            alist,
            xclist,
            sig: INCHI_IOS_STRING {
                pStr: signature,
                nAllocatedLength: 4,
                nUsedLength: 3,
                nPtr: 2,
            },
        };
        assert_eq!(heap.live_allocation_count(), 3);
        assert_eq!(DiylFrag_Free(&mut heap, Some(&mut fragment)), Ok(()));
        assert_eq!(fragment.alist, SourceMutPointer::null());
        assert_eq!(fragment.xclist, SourceMutPointer::null());
        assert_eq!(fragment.sig, INCHI_IOS_STRING::default());
        assert_eq!((fragment.na, fragment.nb, fragment.end1, fragment.end2), (3, 2, 11, 19));
        assert!(heap.slice(alist.as_const()).is_err());
        assert!(heap.slice(xclist.as_const()).is_err());
        assert!(heap.slice(signature.as_const()).is_err());
        assert_eq!(heap.live_allocation_count(), 0);

        assert_eq!(DiylFrag_Free(&mut heap, Some(&mut fragment)), Ok(()));
        assert_eq!(
            fragment,
            DiylFrag {
                na: 3,
                nb: 2,
                end1: 11,
                end2: 19,
                ..DiylFrag::default()
            }
        );

        let only_xclist = heap.allocate_model_storage(vec![7_i32]).unwrap();
        let mut xclist_fragment = DiylFrag {
            na: i32::MIN,
            nb: i32::MAX,
            end1: -1,
            end2: -2,
            xclist: only_xclist,
            sig: INCHI_IOS_STRING {
                pStr: SourceMutPointer::null(),
                nAllocatedLength: 91,
                nUsedLength: 92,
                nPtr: 93,
            },
            ..DiylFrag::default()
        };
        assert_eq!(DiylFrag_Free(&mut heap, Some(&mut xclist_fragment)), Ok(()));
        assert_eq!(xclist_fragment.alist, SourceMutPointer::null());
        assert_eq!(xclist_fragment.xclist, SourceMutPointer::null());
        assert_eq!(xclist_fragment.sig, INCHI_IOS_STRING::default());
        assert_eq!(
            (
                xclist_fragment.na,
                xclist_fragment.nb,
                xclist_fragment.end1,
                xclist_fragment.end2,
            ),
            (i32::MIN, i32::MAX, -1, -2)
        );
        assert!(heap.slice(only_xclist.as_const()).is_err());

        let only_alist = heap.allocate_model_storage(vec![9_i32]).unwrap();
        let mut alist_fragment = DiylFrag {
            alist: only_alist,
            ..DiylFrag::default()
        };
        assert_eq!(DiylFrag_Free(&mut heap, Some(&mut alist_fragment)), Ok(()));
        assert_eq!(alist_fragment, DiylFrag::default());
        assert!(heap.slice(only_alist.as_const()).is_err());
        assert_eq!(heap.live_allocation_count(), 0);
    }

    #[test]
    fn source_port__runichi2__diylfrag_new__line_2007() {
        fn release(heap: &mut SourceHeap, pointer: SourceMutPointer<DiylFrag>) {
            heap.with_slice_mut_and_heap_mut(pointer, |fragments, heap| DiylFrag_Free(heap, fragments.first_mut()))
                .unwrap();
            inchi_free(heap, pointer).unwrap();
        }

        let mut heap = SourceHeap::default();
        let signature = c_string(&mut heap, b"");
        heap.trace_source_allocations();
        let fragment_pointer = DiylFrag_New(&mut heap, 3, -7, 19, signature.as_const()).unwrap();
        assert!(!fragment_pointer.is_null());
        assert_eq!(heap.source_allocation_calls(), 4);
        let fragment = heap.slice(fragment_pointer.as_const()).unwrap()[0].clone();
        assert_eq!((fragment.na, fragment.nb, fragment.end1, fragment.end2), (3, 0, -7, 19));
        assert!(!fragment.alist.is_null());
        assert!(!fragment.xclist.is_null());
        assert_eq!(heap.slice(fragment.alist.as_const()).unwrap(), &[0, 0, 0]);
        assert_eq!(heap.slice(fragment.xclist.as_const()).unwrap(), &[0, 0, 0]);
        assert_eq!(fragment.sig.nUsedLength, 0);
        assert_eq!(fragment.sig.nPtr, 0);
        assert_eq!(fragment.sig.nAllocatedLength, 1);
        assert_eq!(heap.slice(fragment.sig.pStr.as_const()).unwrap(), &[0]);
        release(&mut heap, fragment_pointer);
        assert_eq!(heap.live_allocation_count(), 1);
        heap.free(signature).unwrap();
        assert_eq!(heap.live_allocation_count(), 0);

        for na in [0, -1, i32::MIN] {
            let mut branch_heap = SourceHeap::default();
            let empty = c_string(&mut branch_heap, b"");
            branch_heap.trace_source_allocations();
            let pointer = DiylFrag_New(&mut branch_heap, na, i32::MIN, i32::MAX, empty.as_const()).unwrap();
            assert!(!pointer.is_null(), "na={na}");
            assert_eq!(branch_heap.source_allocation_calls(), 2, "na={na}");
            let fragment = &branch_heap.slice(pointer.as_const()).unwrap()[0];
            assert_eq!(fragment.na, na);
            assert_eq!(fragment.nb, 0);
            assert_eq!(fragment.end1, i32::MIN);
            assert_eq!(fragment.end2, i32::MAX);
            assert!(fragment.alist.is_null());
            assert!(fragment.xclist.is_null());
            assert_eq!(fragment.sig.nUsedLength, 0);
            assert_eq!(branch_heap.slice(fragment.sig.pStr.as_const()).unwrap()[0], 0);
            release(&mut branch_heap, pointer);
            assert_eq!(branch_heap.live_allocation_count(), 1);
        }

        for successful_allocations in 0..3 {
            let mut failure_heap = SourceHeap::default();
            let value = c_string(&mut failure_heap, b"x");
            failure_heap.fail_after_allocations(successful_allocations);
            let pointer = DiylFrag_New(&mut failure_heap, 2, 1, 2, value.as_const()).unwrap();
            assert!(pointer.is_null(), "successful_allocations={successful_allocations}");
            assert_eq!(
                failure_heap.live_allocation_count(),
                1,
                "successful_allocations={successful_allocations}"
            );
            assert_eq!(
                failure_heap.source_allocation_calls(),
                if successful_allocations == 0 { 1 } else { 3 },
                "successful_allocations={successful_allocations}"
            );
        }

        let mut signature_failure_heap = SourceHeap::default();
        let value = c_string(&mut signature_failure_heap, b"x");
        signature_failure_heap.fail_after_allocations(3);
        assert_eq!(
            DiylFrag_New(&mut signature_failure_heap, 2, 1, 2, value.as_const()),
            Err(SourceHeapError::AllocationFailed)
        );
        assert_eq!(signature_failure_heap.source_allocation_calls(), 4);
        assert_eq!(signature_failure_heap.live_allocation_count(), 1);

        let mut zero_signature_failure_heap = SourceHeap::default();
        let value = c_string(&mut zero_signature_failure_heap, b"x");
        zero_signature_failure_heap.fail_after_allocations(1);
        assert_eq!(
            DiylFrag_New(&mut zero_signature_failure_heap, 0, 1, 2, value.as_const(),),
            Err(SourceHeapError::AllocationFailed)
        );
        assert_eq!(zero_signature_failure_heap.source_allocation_calls(), 2);
        assert_eq!(zero_signature_failure_heap.live_allocation_count(), 1);

        let mut nonempty_heap = SourceHeap::default();
        let nonempty = c_string(&mut nonempty_heap, b"fragment");
        assert_eq!(
            DiylFrag_New(&mut nonempty_heap, 2, 1, 2, nonempty.as_const()),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(nonempty_heap.live_allocation_count(), 1);

        let mut invalid_heap = SourceHeap::default();
        let missing_nul = invalid_heap
            .allocate_model_storage(vec![b'b' as i8, b'a' as i8, b'd' as i8])
            .unwrap();
        assert_eq!(
            DiylFrag_New(&mut invalid_heap, 2, 1, 2, missing_nul.as_const()),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(invalid_heap.live_allocation_count(), 1);
        assert_eq!(
            DiylFrag_New(&mut invalid_heap, 2, 1, 2, SourceConstPointer::null(),),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(invalid_heap.live_allocation_count(), 1);
    }

    #[test]
    fn source_port__runichi2__analyze_cru_folding__line_2146() {
        #[derive(Clone)]
        struct Fixture {
            data: ORIG_ATOM_DATA,
            unit: SourceMutPointer<OAD_PolymerUnit>,
            edits: OAD_StructureEdits,
        }

        fn int_array(heap: &mut SourceHeap, capacity: i32, increment: i32) -> SourceMutPointer<INT_ARRAY> {
            let item = heap
                .allocate_model_storage(vec![-99_i32; capacity.max(0) as usize])
                .unwrap();
            heap.allocate_model_storage(vec![INT_ARRAY {
                item,
                allocated: capacity,
                used: 0,
                increment,
            }])
            .unwrap()
        }

        fn fixture(heap: &mut SourceHeap, unit_atom_count: usize, edit_capacity: i32) -> Fixture {
            let total_atoms = unit_atom_count + 2;
            let mut atoms = vec![inp_ATOM::default(); total_atoms];
            for index in 0..total_atoms {
                atoms[index].orig_at_number = (index + 1) as _;
                let mut degree = 0_usize;
                if index != 0 {
                    atoms[index].neighbor[degree] = (index - 1) as _;
                    atoms[index].bond_type[degree] = 1;
                    degree += 1;
                }
                if index + 1 < total_atoms {
                    atoms[index].neighbor[degree] = (index + 1) as _;
                    atoms[index].bond_type[degree] = 1;
                    degree += 1;
                }
                atoms[index].valence = degree as i8;
            }
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let atom_list = heap
                .allocate_model_storage((2..=(unit_atom_count + 1) as i32).collect::<Vec<_>>())
                .unwrap();
            let backbone_rows = (0..unit_atom_count + 2)
                .map(|_| heap.allocate_model_storage(vec![-1_i32, -1]).unwrap())
                .collect::<Vec<_>>();
            let backbone_bonds = heap.allocate_model_storage(backbone_rows).unwrap();
            let unit = heap
                .allocate_model_storage(vec![OAD_PolymerUnit {
                    na: unit_atom_count as i32,
                    nb: 2,
                    cap1: 1,
                    end_atom1: 2,
                    end_atom2: (unit_atom_count + 1) as i32,
                    cap2: (unit_atom_count + 2) as i32,
                    alist: atom_list,
                    maxbkbonds: (unit_atom_count + 2) as i32,
                    nbkbonds: -77,
                    bkbonds: backbone_bonds,
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
            let edits = OAD_StructureEdits {
                del_atom: int_array(heap, edit_capacity, 8),
                del_bond: int_array(heap, edit_capacity, 8),
                mod_bond: int_array(heap, edit_capacity, 8),
                mod_coord: int_array(heap, edit_capacity, 8),
                del_side_chains: 73,
                ..OAD_StructureEdits::default()
            };
            Fixture {
                data: ORIG_ATOM_DATA {
                    at: atoms,
                    num_inp_atoms: total_atoms as i32,
                    num_inp_bonds: total_atoms.saturating_sub(1) as i32,
                    polymer,
                    ..ORIG_ATOM_DATA::default()
                },
                unit,
                edits,
            }
        }

        fn values(heap: &SourceHeap, pointer: SourceMutPointer<INT_ARRAY>) -> Vec<i32> {
            let array = &heap.slice(pointer.as_const()).unwrap()[0];
            heap.slice(array.item.as_const()).unwrap()[..array.used as usize].to_vec()
        }

        fn prime_edit_arrays(heap: &mut SourceHeap, edits: &OAD_StructureEdits) {
            for pointer in [edits.del_atom, edits.del_bond, edits.mod_bond, edits.mod_coord] {
                heap.slice_mut(pointer).unwrap()[0].used = 1;
            }
        }

        fn all_unit_bonds(unit_atom_count: usize) -> Vec<i32> {
            (2..=(unit_atom_count as i32))
                .flat_map(|atom| [atom, atom + 1])
                .collect()
        }

        let mut heap = SourceHeap::default();
        let normal = fixture(&mut heap, 6, 16);
        let baseline = heap.live_allocation_count();
        let xclasses = vec![0_i32; 9];
        heap.trace_source_allocations();
        assert_eq!(
            analyze_CRU_folding(&mut heap, &normal.data, 0, 1, &[4, 5], 1, &xclasses, &normal.edits,),
            Ok(_IS_OKAY as i32)
        );
        let allocations_before_edits = heap.source_allocation_calls();
        assert_eq!(allocations_before_edits, 59);
        assert_eq!(heap.live_allocation_count(), baseline);
        assert_eq!(heap.slice(normal.unit.as_const()).unwrap()[0].nbkbonds, 5);
        assert_eq!(values(&heap, normal.edits.del_bond), vec![4, 5]);
        assert_eq!(values(&heap, normal.edits.mod_bond), vec![7, 8, 4, 8]);
        assert_eq!(values(&heap, normal.edits.del_atom), vec![5, 6, 7]);
        assert_eq!(values(&heap, normal.edits.mod_coord), vec![5, 8]);
        assert_eq!(normal.edits.del_side_chains, 73);

        let mut single_atom_fragment_heap = SourceHeap::default();
        let single_atom_fragment = fixture(&mut single_atom_fragment_heap, 6, 16);
        let single_atom_bonds = all_unit_bonds(6);
        assert_eq!(
            analyze_CRU_folding(
                &mut single_atom_fragment_heap,
                &single_atom_fragment.data,
                0,
                5,
                &single_atom_bonds,
                1,
                &xclasses,
                &single_atom_fragment.edits,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(
            values(&single_atom_fragment_heap, single_atom_fragment.edits.del_bond),
            vec![2, 3]
        );
        assert_eq!(
            values(&single_atom_fragment_heap, single_atom_fragment.edits.del_atom),
            vec![3, 4, 5, 6, 7]
        );
        assert_eq!(
            values(&single_atom_fragment_heap, single_atom_fragment.edits.mod_coord),
            vec![3, 8]
        );

        let mut no_backbone_heap = SourceHeap::default();
        let no_backbone = fixture(&mut no_backbone_heap, 1, 16);
        assert_eq!(
            analyze_CRU_folding(
                &mut no_backbone_heap,
                &no_backbone.data,
                0,
                0,
                &[],
                1,
                &[0_i32; 4],
                &no_backbone.edits,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(
            no_backbone_heap.slice(no_backbone.unit.as_const()).unwrap()[0].nbkbonds,
            0
        );
        assert!(values(&no_backbone_heap, no_backbone.edits.del_atom).is_empty());

        let mut no_cut_heap = SourceHeap::default();
        let no_cut = fixture(&mut no_cut_heap, 6, 16);
        let no_cut_baseline = no_cut_heap.live_allocation_count();
        assert_eq!(
            analyze_CRU_folding(
                &mut no_cut_heap,
                &no_cut.data,
                0,
                1,
                &[1, 2],
                1,
                &xclasses,
                &no_cut.edits,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(no_cut_heap.live_allocation_count(), no_cut_baseline);
        assert!(values(&no_cut_heap, no_cut.edits.del_bond).is_empty());

        let mut distinct_heap = SourceHeap::default();
        let distinct = fixture(&mut distinct_heap, 6, 16);
        let mut distinct_classes = vec![0_i32; 9];
        for atom in 2..=7 {
            distinct_classes[atom] = atom as i32;
        }
        assert_eq!(
            analyze_CRU_folding(
                &mut distinct_heap,
                &distinct.data,
                0,
                1,
                &[4, 5],
                8,
                &distinct_classes,
                &distinct.edits,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert!(values(&distinct_heap, distinct.edits.del_atom).is_empty());

        for (classes, expected_branch) in [
            ([0, 1, 0, 2, 0], "no repeating subsequence"),
            ([0, 1, 0, 1, 0], "non-divisible repeat"),
            ([0, 1, 2, 0, 1], "single fold"),
        ] {
            let mut branch_heap = SourceHeap::default();
            let branch = fixture(&mut branch_heap, 5, 16);
            let bonds = all_unit_bonds(5);
            let mut branch_classes = vec![0_i32; 8];
            branch_classes[2..7].copy_from_slice(&classes);
            assert_eq!(
                analyze_CRU_folding(
                    &mut branch_heap,
                    &branch.data,
                    0,
                    4,
                    &bonds,
                    3,
                    &branch_classes,
                    &branch.edits,
                ),
                Ok(_IS_OKAY as i32),
                "{expected_branch}"
            );
            assert!(
                values(&branch_heap, branch.edits.del_atom).is_empty(),
                "{expected_branch}"
            );
        }

        for successful_allocations in [0_u64, 1, 15, 16, 17, 21] {
            let mut failure_heap = SourceHeap::default();
            let failure = fixture(&mut failure_heap, 6, 16);
            let failure_baseline = failure_heap.live_allocation_count();
            failure_heap.fail_after_allocations(successful_allocations);
            assert_eq!(
                analyze_CRU_folding(
                    &mut failure_heap,
                    &failure.data,
                    0,
                    1,
                    &[4, 5],
                    1,
                    &xclasses,
                    &failure.edits,
                ),
                Ok(_IS_ERROR as i32),
                "allocation ordinal {successful_allocations}"
            );
            assert_eq!(
                failure_heap.live_allocation_count(),
                failure_baseline,
                "allocation ordinal {successful_allocations}"
            );
            assert!(values(&failure_heap, failure.edits.del_bond).is_empty());
        }

        let mut first_append_heap = SourceHeap::default();
        let first_append = fixture(&mut first_append_heap, 6, 1);
        prime_edit_arrays(&mut first_append_heap, &first_append.edits);
        let first_append_baseline = first_append_heap.live_allocation_count();
        first_append_heap.fail_after_allocations(allocations_before_edits);
        assert_eq!(
            analyze_CRU_folding(
                &mut first_append_heap,
                &first_append.data,
                0,
                1,
                &[4, 5],
                1,
                &xclasses,
                &first_append.edits,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(first_append_heap.live_allocation_count(), first_append_baseline - 1);
        let failed_del_bond = &first_append_heap.slice(first_append.edits.del_bond.as_const()).unwrap()[0];
        assert!(failed_del_bond.item.is_null());
        assert_eq!((failed_del_bond.allocated, failed_del_bond.used), (1, 1));
        assert_eq!(
            values(&first_append_heap, first_append.edits.mod_bond),
            vec![-99, 7, 8, 4, 8]
        );
        assert_eq!(values(&first_append_heap, first_append.edits.del_atom), vec![-99]);
        assert_eq!(values(&first_append_heap, first_append.edits.mod_coord), vec![-99]);

        let mut atom_append_heap = SourceHeap::default();
        let atom_append = fixture(&mut atom_append_heap, 6, 1);
        prime_edit_arrays(&mut atom_append_heap, &atom_append.edits);
        let atom_append_baseline = atom_append_heap.live_allocation_count();
        atom_append_heap.fail_after_allocations(allocations_before_edits + 2);
        assert_eq!(
            analyze_CRU_folding(
                &mut atom_append_heap,
                &atom_append.data,
                0,
                1,
                &[4, 5],
                1,
                &xclasses,
                &atom_append.edits,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(atom_append_heap.live_allocation_count(), atom_append_baseline - 1);
        assert_eq!(values(&atom_append_heap, atom_append.edits.del_bond), vec![-99, 4, 5]);
        assert_eq!(
            values(&atom_append_heap, atom_append.edits.mod_bond),
            vec![-99, 7, 8, 4, 8]
        );
        let failed_del_atom = &atom_append_heap.slice(atom_append.edits.del_atom.as_const()).unwrap()[0];
        assert!(failed_del_atom.item.is_null());
        assert_eq!((failed_del_atom.allocated, failed_del_atom.used), (1, 1));
        assert_eq!(values(&atom_append_heap, atom_append.edits.mod_coord), vec![-99]);

        let mut coord_append_heap = SourceHeap::default();
        let coord_append = fixture(&mut coord_append_heap, 6, 1);
        prime_edit_arrays(&mut coord_append_heap, &coord_append.edits);
        let coord_append_baseline = coord_append_heap.live_allocation_count();
        coord_append_heap.fail_after_allocations(allocations_before_edits + 3);
        assert_eq!(
            analyze_CRU_folding(
                &mut coord_append_heap,
                &coord_append.data,
                0,
                1,
                &[4, 5],
                1,
                &xclasses,
                &coord_append.edits,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(coord_append_heap.live_allocation_count(), coord_append_baseline - 1);
        assert_eq!(
            values(&coord_append_heap, coord_append.edits.del_atom),
            vec![-99, 5, 6, 7]
        );
        let failed_coord = &coord_append_heap
            .slice(coord_append.edits.mod_coord.as_const())
            .unwrap()[0];
        assert!(failed_coord.item.is_null());
        assert_eq!((failed_coord.allocated, failed_coord.used), (1, 1));

        let mut invalid_heap = SourceHeap::default();
        let invalid = fixture(&mut invalid_heap, 6, 16);
        assert_eq!(
            analyze_CRU_folding(
                &mut invalid_heap,
                &invalid.data,
                -1,
                1,
                &[4, 5],
                1,
                &xclasses,
                &invalid.edits,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            analyze_CRU_folding(
                &mut invalid_heap,
                &invalid.data,
                0,
                1,
                &[4],
                1,
                &xclasses,
                &invalid.edits,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__runichi2__oad_polymer_preparefoldcruedits__line_1838() {
        const POLYMER_INCHI: &[u8] = b"InChI=1B/C4H4N4.2Zz/c1-5-2-7-4-8-3-6-1;;/h1-4H;;/z101-1-8(9,10-8,3,1,6,2,5,2,7,3,6,1,5,4,7,4,8)/b5-1-,5-2+,6-1+,6-3-,7-2+,7-4+,8-3+,8-4+;;";
        const AUX_MAPPING: &[u8] = b"AuxInfo=1/0/N:2,3,5,1,6,8,7,4/E:(1,2,3,4,5,6,7,8)/";
        const LACTIC_ODD: &[u8] = b"InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/t2-/m0/s1";
        const LACTIC_EVEN: &[u8] = b"InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/t2-/m1/s1";

        fn named_atom(name: &[u8]) -> inp_ATOM {
            let mut atom = inp_ATOM::default();
            for (target, source) in atom.elname.iter_mut().zip(name.iter().copied()) {
                *target = source as i8;
            }
            atom
        }

        fn int_array(heap: &mut SourceHeap) -> SourceMutPointer<INT_ARRAY> {
            let item = heap.allocate_model_storage(vec![-91_i32; 16]).unwrap();
            heap.allocate_model_storage(vec![INT_ARRAY {
                item,
                allocated: 16,
                used: 0,
                increment: 8,
            }])
            .unwrap()
        }

        fn edits(heap: &mut SourceHeap) -> OAD_StructureEdits {
            OAD_StructureEdits {
                del_atom: int_array(heap),
                del_bond: int_array(heap),
                mod_bond: int_array(heap),
                mod_coord: int_array(heap),
                ..OAD_StructureEdits::default()
            }
        }

        fn values(heap: &SourceHeap, pointer: SourceMutPointer<INT_ARRAY>) -> Vec<i32> {
            let array = &heap.slice(pointer.as_const()).unwrap()[0];
            heap.slice(array.item.as_const()).unwrap()[..array.used as usize].to_vec()
        }

        fn empty_polymer_data(heap: &mut SourceHeap, atom_count: usize, bond_count: i32) -> ORIG_ATOM_DATA {
            let atoms = heap
                .allocate_model_storage(vec![inp_ATOM::default(); atom_count])
                .unwrap();
            let polymer = heap.allocate_model_storage(vec![OAD_Polymer::default()]).unwrap();
            ORIG_ATOM_DATA {
                at: atoms,
                num_inp_atoms: atom_count as i32,
                num_inp_bonds: bond_count,
                polymer,
                ..ORIG_ATOM_DATA::default()
            }
        }

        fn chain_fixture(
            heap: &mut SourceHeap,
            cut_bond_type: u8,
            cap_names: (&[u8], &[u8]),
            unit_atom_count: i32,
            crossing_bond_count: i32,
        ) -> (ORIG_ATOM_DATA, SourceMutPointer<OAD_PolymerUnit>, OAD_StructureEdits) {
            let mut atoms = (0..8)
                .map(|index| {
                    if index == 0 {
                        named_atom(cap_names.0)
                    } else if index == 7 {
                        named_atom(cap_names.1)
                    } else {
                        named_atom(b"C\0")
                    }
                })
                .collect::<Vec<_>>();
            for index in 0..atoms.len() {
                let mut degree = 0_usize;
                if index != 0 {
                    atoms[index].neighbor[degree] = (index - 1) as _;
                    atoms[index].bond_type[degree] = if index == 4 { cut_bond_type } else { 1 };
                    degree += 1;
                }
                if index + 1 < atoms.len() {
                    atoms[index].neighbor[degree] = (index + 1) as _;
                    atoms[index].bond_type[degree] = if index == 3 { cut_bond_type } else { 1 };
                    degree += 1;
                }
                atoms[index].valence = degree as i8;
                atoms[index].chem_bonds_valence = atoms[index].bond_type[..degree]
                    .iter()
                    .fold(0_i8, |sum, value| sum.wrapping_add(*value as i8));
                atoms[index].orig_at_number = (index + 1) as _;
            }
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let atom_list = heap.allocate_model_storage((2..=7).collect::<Vec<i32>>()).unwrap();
            let crossing_bonds = heap.allocate_model_storage(vec![1_i32, 2, 8, 7]).unwrap();
            let unit = heap
                .allocate_model_storage(vec![OAD_PolymerUnit {
                    type_: crate::source_types::POLYMER_STY_SRU as i32,
                    subtype: crate::source_types::POLYMER_SST_NON as i32,
                    conn: crate::source_types::POLYMER_CONN_HT as i32,
                    na: unit_atom_count,
                    nb: crossing_bond_count,
                    alist: atom_list,
                    blist: crossing_bonds,
                    ..OAD_PolymerUnit::default()
                }])
                .unwrap();
            let units = heap.allocate_model_storage(vec![unit]).unwrap();
            let polymer = heap
                .allocate_model_storage(vec![OAD_Polymer {
                    n: 1,
                    units,
                    treat: POLYMERS_MODERN as i32,
                    ..OAD_Polymer::default()
                }])
                .unwrap();
            (
                ORIG_ATOM_DATA {
                    at: atoms,
                    num_inp_atoms: 8,
                    num_inp_bonds: 8,
                    polymer,
                    ..ORIG_ATOM_DATA::default()
                },
                unit,
                edits(heap),
            )
        }

        fn invoke(
            heap: &mut SourceHeap,
            data: &mut ORIG_ATOM_DATA,
            no_edits: &[u8],
            inchi: &[u8],
            aux: &[u8],
            edits: &OAD_StructureEdits,
        ) -> Result<i32, SourceHeapError> {
            let no_edits = c_string(heap, no_edits);
            let inchi = c_string(heap, inchi);
            let aux = c_string(heap, aux);
            OAD_Polymer_PrepareFoldCRUEdits(
                heap,
                data,
                no_edits.as_const(),
                SourceConstPointer::null(),
                inchi.as_const(),
                aux.as_const(),
                edits,
            )
        }

        let mut early_heap = SourceHeap::default();
        let mut early = empty_polymer_data(&mut early_heap, 1, 0);
        let early_edits = edits(&mut early_heap);
        let aux_without_classes = c_string(&mut early_heap, b"AuxInfo=1/0/N:1");
        let early_baseline = early_heap.live_allocation_count();
        assert_eq!(
            OAD_Polymer_PrepareFoldCRUEdits(
                &mut early_heap,
                &mut early,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                aux_without_classes.as_const(),
                &early_edits,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(early_heap.live_allocation_count(), early_baseline);

        let malformed_aux = c_string(&mut early_heap, b"/N:");
        assert_eq!(
            OAD_Polymer_PrepareFoldCRUEdits(
                &mut early_heap,
                &mut early,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                malformed_aux.as_const(),
                &early_edits,
            ),
            Ok(_IS_ERROR as i32)
        );
        let malformed_classes = c_string(&mut early_heap, b"/N:1/E:(1X/");
        assert_eq!(
            OAD_Polymer_PrepareFoldCRUEdits(
                &mut early_heap,
                &mut early,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                malformed_classes.as_const(),
                &early_edits,
            ),
            Ok(_IS_ERROR as i32)
        );

        let mut allocation_heap = SourceHeap::default();
        let mut allocation_data = empty_polymer_data(&mut allocation_heap, 1, 0);
        let allocation_edits = edits(&mut allocation_heap);
        let allocation_aux = c_string(&mut allocation_heap, b"AuxInfo=1/0/N:1");
        let allocation_baseline = allocation_heap.live_allocation_count();
        allocation_heap.fail_after_allocations(0);
        assert_eq!(
            OAD_Polymer_PrepareFoldCRUEdits(
                &mut allocation_heap,
                &mut allocation_data,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                allocation_aux.as_const(),
                &allocation_edits,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(allocation_heap.live_allocation_count(), allocation_baseline);

        let mut stereo_error_heap = SourceHeap::default();
        let mut stereo_error_data = empty_polymer_data(&mut stereo_error_heap, 1, 0);
        let stereo_error_edits = edits(&mut stereo_error_heap);
        assert_eq!(
            invoke(
                &mut stereo_error_heap,
                &mut stereo_error_data,
                b"InChI=1S/CH4/i\n",
                b"",
                b"AuxInfo=1/0/N:1/E:(1)/",
                &stereo_error_edits,
            ),
            Ok(_IS_ERROR as i32)
        );

        for stereochemical_input in [LACTIC_ODD, LACTIC_EVEN] {
            let mut traversal_heap = SourceHeap::default();
            let mut traversal = empty_polymer_data(&mut traversal_heap, 8, 8);
            let traversal_edits = edits(&mut traversal_heap);
            let baseline = traversal_heap.live_allocation_count();
            assert_eq!(
                invoke(
                    &mut traversal_heap,
                    &mut traversal,
                    stereochemical_input,
                    POLYMER_INCHI,
                    AUX_MAPPING,
                    &traversal_edits,
                ),
                Ok(_IS_OKAY as i32)
            );
            assert_eq!(traversal_heap.live_allocation_count(), baseline + 3);
            assert!(values(&traversal_heap, traversal_edits.del_atom).is_empty());
        }

        let mut backbone_error_heap = SourceHeap::default();
        let mut backbone_error = empty_polymer_data(&mut backbone_error_heap, 8, 8);
        let backbone_error_edits = edits(&mut backbone_error_heap);
        assert_eq!(
            invoke(
                &mut backbone_error_heap,
                &mut backbone_error,
                LACTIC_ODD,
                b"InChI=1/C32767/z()",
                AUX_MAPPING,
                &backbone_error_edits,
            ),
            Ok(_IS_ERROR as i32)
        );

        let mut validation_heap = SourceHeap::default();
        let mut validation = empty_polymer_data(&mut validation_heap, 8, 8);
        validation_heap.slice_mut(validation.polymer).unwrap()[0].n = -1;
        let validation_edits = edits(&mut validation_heap);
        assert_eq!(
            invoke(
                &mut validation_heap,
                &mut validation,
                LACTIC_ODD,
                POLYMER_INCHI,
                AUX_MAPPING,
                &validation_edits,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(validation.valid_polymer, 0);

        let mut fold_heap = SourceHeap::default();
        let (mut fold, fold_unit, fold_edits) = chain_fixture(&mut fold_heap, 1, (b"Zz\0", b"Zz\0"), 6, 2);
        assert_eq!(
            invoke(
                &mut fold_heap,
                &mut fold,
                POLYMER_INCHI,
                POLYMER_INCHI,
                AUX_MAPPING,
                &fold_edits,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(fold.valid_polymer, 1);
        assert_eq!(fold_heap.slice(fold_unit.as_const()).unwrap()[0].nbkbonds, 5);
        assert_eq!(values(&fold_heap, fold_edits.del_bond), vec![4, 5]);
        assert_eq!(values(&fold_heap, fold_edits.mod_bond), vec![7, 8, 4, 8]);
        assert_eq!(values(&fold_heap, fold_edits.del_atom), vec![5, 6, 7]);
        assert_eq!(values(&fold_heap, fold_edits.mod_coord), vec![5, 8]);

        let mut multiple_heap = SourceHeap::default();
        let (mut multiple, multiple_unit, multiple_edits) =
            chain_fixture(&mut multiple_heap, 2, (b"Zz\0", b"Zz\0"), 6, 2);
        assert_eq!(
            invoke(
                &mut multiple_heap,
                &mut multiple,
                POLYMER_INCHI,
                POLYMER_INCHI,
                AUX_MAPPING,
                &multiple_edits,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(multiple_heap.slice(multiple_unit.as_const()).unwrap()[0].nbkbonds, 5);
        assert!(values(&multiple_heap, multiple_edits.del_bond).is_empty());
        assert!(values(&multiple_heap, multiple_edits.del_atom).is_empty());

        for (unit_atoms, crossing_bonds, cap_names) in [
            (1, 2, (b"Zz\0".as_slice(), b"Zz\0".as_slice())),
            (6, 0, (b"Zz\0".as_slice(), b"Zz\0".as_slice())),
            (6, 2, (b"C\0".as_slice(), b"Zz\0".as_slice())),
            (6, 2, (b"Zz\0".as_slice(), b"C\0".as_slice())),
        ] {
            let mut skip_heap = SourceHeap::default();
            let (mut skip, _, skip_edits) = chain_fixture(&mut skip_heap, 1, cap_names, unit_atoms, crossing_bonds);
            assert_eq!(
                invoke(
                    &mut skip_heap,
                    &mut skip,
                    POLYMER_INCHI,
                    POLYMER_INCHI,
                    AUX_MAPPING,
                    &skip_edits,
                ),
                Ok(_IS_OKAY as i32)
            );
            assert!(values(&skip_heap, skip_edits.del_atom).is_empty());
            assert!(values(&skip_heap, skip_edits.del_bond).is_empty());
        }
    }

    #[test]
    fn source_port__runichi2__count_colors_in_sequence__line_2463() {
        let mut counts = [91_i32; 7];
        assert_eq!(count_colors_in_sequence(&[2, 0, 2, -1, 3, 0], 6, 5, &mut counts), Ok(3));
        assert_eq!(counts, [2, 0, 2, 1, 0, 91, 91]);

        let mut negative_counts = [7_i32; 4];
        assert_eq!(
            count_colors_in_sequence(&[-1, i32::MIN, -7], 3, 4, &mut negative_counts),
            Ok(0)
        );
        assert_eq!(negative_counts, [0, 0, 0, 0]);

        for n in [0, -1, i32::MIN] {
            let mut empty_counts = [5_i32; 3];
            assert_eq!(count_colors_in_sequence(&[], n, 2, &mut empty_counts), Ok(0), "n={n}");
            assert_eq!(empty_counts, [0, 0, 5], "n={n}");
        }

        let mut zero_counts = [4_i32; 2];
        assert_eq!(count_colors_in_sequence(&[], 0, 0, &mut zero_counts), Ok(0));
        assert_eq!(zero_counts, [4, 4]);

        let mut zero_limit_counts = [4_i32; 2];
        assert_eq!(count_colors_in_sequence(&[0], 1, 0, &mut zero_limit_counts), Ok(0));
        assert_eq!(zero_limit_counts, [5, 4]);

        let mut negative_max_counts = [3_i32; 2];
        assert_eq!(
            count_colors_in_sequence(&[], 0, -1, &mut negative_max_counts),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(negative_max_counts, [3, 3]);

        let mut short_counts = [3_i32; 2];
        assert_eq!(
            count_colors_in_sequence(&[], 0, 3, &mut short_counts),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(short_counts, [3, 3]);

        let mut short_color_counts = [9_i32; 4];
        assert_eq!(
            count_colors_in_sequence(&[1, 1], 3, 4, &mut short_color_counts),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(short_color_counts, [0, 2, 0, 0]);

        let mut out_of_range_counts = [9_i32; 4];
        assert_eq!(
            count_colors_in_sequence(&[1, 4, 2], 3, 4, &mut out_of_range_counts),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(out_of_range_counts, [0, 1, 0, 0]);

        let mut boundary_counts = [5_i32; 4];
        assert_eq!(count_colors_in_sequence(&[3, 3, 3], 3, 4, &mut boundary_counts), Ok(1));
        assert_eq!(boundary_counts, [0, 0, 0, 3]);
    }

    #[test]
    fn source_port__runichi2__len_repeating_subsequence__line_2489() {
        for n in [i32::MIN, -1, 0, 1] {
            assert_eq!(len_repeating_subsequence(None, None, n), Ok(0), "n={n}");
            assert_eq!(len_repeating_subsequence(Some(&[]), Some(&[]), n), Ok(0), "n={n}");
        }
        assert_eq!(len_repeating_subsequence(None, None, 2), Ok(0));
        assert_eq!(len_repeating_subsequence(None, None, i32::MAX), Ok(0));
        assert_eq!(
            len_repeating_subsequence(Some(&[]), None, i32::MAX),
            Err(SourceHeapError::SourceIntegerOverflow)
        );

        assert_eq!(len_repeating_subsequence(Some(&[7, 7]), None, 2), Ok(1));
        assert_eq!(len_repeating_subsequence(Some(&[1, 2]), None, 2), Ok(0));
        assert_eq!(len_repeating_subsequence(Some(&[1, 2, 1, 2]), None, 4), Ok(2));
        assert_eq!(len_repeating_subsequence(Some(&[1, 2, 1, 2, 1]), None, 5), Ok(2));
        assert_eq!(len_repeating_subsequence(Some(&[1, 2, 3, 1, 2, 3]), None, 6), Ok(3));
        assert_eq!(len_repeating_subsequence(Some(&[1, 2, 3, 1, 2]), None, 5), Ok(3));
        assert_eq!(len_repeating_subsequence(Some(&[1, 2, 3, 4]), None, 4), Ok(0));

        assert_eq!(
            len_repeating_subsequence(Some(&[9, 9, 9, 9]), Some(&[1, 2, 1, 2]), 4),
            Ok(2)
        );
        assert_eq!(len_repeating_subsequence(Some(&[9, 9, 9]), Some(&[4, 4, 4]), 3), Ok(1));
        assert_eq!(len_repeating_subsequence(Some(&[1, 2]), Some(&[]), 2), Ok(0));
        assert_eq!(
            len_repeating_subsequence(Some(&[1, 1]), Some(&[]), 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            len_repeating_subsequence(Some(&[1]), None, 3),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            len_repeating_subsequence(Some(&[1, 2, 1]), Some(&[5, 6]), 3),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__runichi2__diylfrag_makesignature__line_2069() {
        fn new_fragment(
            heap: &mut SourceHeap,
            na: i32,
            end1: i32,
            end2: i32,
            atoms: &[i32],
        ) -> SourceMutPointer<DiylFrag> {
            let empty = c_string(heap, b"");
            let fragment = DiylFrag_New(heap, na, end1, end2, empty.as_const()).unwrap();
            heap.free(empty).unwrap();
            assert!(!fragment.is_null());
            if na > 0 {
                let alist = heap.slice(fragment.as_const()).unwrap()[0].alist;
                heap.slice_mut(alist).unwrap()[..atoms.len()].copy_from_slice(atoms);
            }
            fragment
        }

        fn make_signature(
            heap: &mut SourceHeap,
            fragment: SourceMutPointer<DiylFrag>,
            nxc: i32,
            xc: &[i32],
            counts: &mut [i32],
        ) -> Result<(), SourceHeapError> {
            heap.with_slice_mut_and_heap_mut(fragment, |fragments, heap| {
                DiylFrag_MakeSignature(
                    heap,
                    fragments.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)?,
                    nxc,
                    xc,
                    counts,
                )
            })
        }

        fn signature(heap: &SourceHeap, fragment: SourceMutPointer<DiylFrag>) -> Vec<u8> {
            let fragment = &heap.slice(fragment.as_const()).unwrap()[0];
            heap.slice(fragment.sig.pStr.as_const()).unwrap()[..fragment.sig.nUsedLength as usize]
                .iter()
                .map(|byte| *byte as u8)
                .collect()
        }

        fn release(heap: &mut SourceHeap, fragment: SourceMutPointer<DiylFrag>) {
            heap.with_slice_mut_and_heap_mut(fragment, |fragments, heap| DiylFrag_Free(heap, fragments.first_mut()))
                .unwrap();
            heap.free(fragment).unwrap();
        }

        let mut heap = SourceHeap::default();
        let fragment = new_fragment(&mut heap, 5, 0, 5, &[0, 1, 2, 4, 5]);
        let mut counts = [91_i32; 5];
        heap.trace_source_allocations();
        assert_eq!(
            make_signature(&mut heap, fragment, 4, &[3, 1, -1, 0, 4, 2], &mut counts,),
            Ok(())
        );
        assert_eq!(counts, [0, 1, 1, 1, 1]);
        assert_eq!(signature(&heap, fragment), b"5,3,2{(1:1)(2:1)(3:1)}");
        assert_eq!(heap.source_allocation_calls(), 1);
        let copied = heap.slice(fragment.as_const()).unwrap()[0].clone();
        assert_eq!(heap.slice(copied.xclist.as_const()).unwrap(), &[3, 1, -1, 4, 2]);

        counts.fill(77);
        assert_eq!(
            make_signature(&mut heap, fragment, 4, &[3, 1, -1, 0, 4, 2], &mut counts,),
            Ok(())
        );
        assert_eq!(counts, [0, 1, 1, 1, 1]);
        assert_eq!(
            signature(&heap, fragment),
            b"5,3,2{(1:1)(2:1)(3:1)}5,3,2{(1:1)(2:1)(3:1)}"
        );
        release(&mut heap, fragment);
        assert_eq!(heap.live_allocation_count(), 0);

        let empty_fragment = new_fragment(&mut heap, 0, 0, 1, &[]);
        let mut empty_counts = [12_i32; 1];
        assert_eq!(
            make_signature(&mut heap, empty_fragment, 0, &[i32::MIN, i32::MAX], &mut empty_counts,),
            Ok(())
        );
        assert_eq!(empty_counts, [0]);
        assert_eq!(signature(&heap, empty_fragment), b"0,-2147483648,2147483647{}");
        release(&mut heap, empty_fragment);

        let excluded_fragment = new_fragment(&mut heap, 1, 0, 0, &[0]);
        let mut excluded_counts = [8_i32; 1];
        assert_eq!(
            make_signature(&mut heap, excluded_fragment, 0, &[0], &mut excluded_counts,),
            Ok(())
        );
        assert_eq!(excluded_counts, [1]);
        assert_eq!(signature(&heap, excluded_fragment), b"1,0,0{}");
        release(&mut heap, excluded_fragment);

        let endpoint_fragment = new_fragment(&mut heap, 1, -1, 0, &[0]);
        let mut endpoint_counts = [6_i32; 2];
        assert_eq!(
            make_signature(&mut heap, endpoint_fragment, 1, &[4], &mut endpoint_counts,),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(endpoint_counts, [6, 6]);
        assert_eq!(signature(&heap, endpoint_fragment), b"");
        release(&mut heap, endpoint_fragment);

        let bad_atom_fragment = new_fragment(&mut heap, 3, 0, 1, &[0, 9, 1]);
        let mut bad_atom_counts = [6_i32; 3];
        assert_eq!(
            make_signature(&mut heap, bad_atom_fragment, 2, &[4, 5], &mut bad_atom_counts,),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(bad_atom_counts, [6, 6, 6]);
        assert_eq!(signature(&heap, bad_atom_fragment), b"3,4,5{");
        let copied = heap.slice(bad_atom_fragment.as_const()).unwrap()[0].clone();
        assert_eq!(heap.slice(copied.xclist.as_const()).unwrap(), &[4, 0, 0]);
        release(&mut heap, bad_atom_fragment);

        let overflow_fragment = new_fragment(&mut heap, 0, 0, 0, &[]);
        let mut overflow_counts = [3_i32; 1];
        assert_eq!(
            make_signature(&mut heap, overflow_fragment, i32::MAX, &[7], &mut overflow_counts,),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(overflow_counts, [3]);
        assert_eq!(signature(&heap, overflow_fragment), b"0,7,7{");
        release(&mut heap, overflow_fragment);

        let short_counts_fragment = new_fragment(&mut heap, 1, 0, 0, &[0]);
        let mut short_counts = [3_i32; 2];
        assert_eq!(
            make_signature(&mut heap, short_counts_fragment, 2, &[1], &mut short_counts,),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(short_counts, [3, 3]);
        assert_eq!(signature(&heap, short_counts_fragment), b"1,1,1{");
        release(&mut heap, short_counts_fragment);

        let first_allocation_failure = new_fragment(&mut heap, 0, 0, 0, &[]);
        let mut failure_counts = [3_i32; 1];
        heap.fail_after_allocations(0);
        assert_eq!(
            make_signature(&mut heap, first_allocation_failure, 0, &[7], &mut failure_counts,),
            Err(SourceHeapError::AllocationFailed)
        );
        assert_eq!(heap.source_allocation_calls(), 1);
        assert_eq!(failure_counts, [3]);
        assert_eq!(signature(&heap, first_allocation_failure), b"");
        release(&mut heap, first_allocation_failure);

        let later_allocation_failure = new_fragment(&mut heap, 8, 0, 7, &[0, 1, 2, 3, 4, 5, 6, 7]);
        let mut later_failure_counts = [3_i32; 9];
        heap.fail_after_allocations(1);
        assert_eq!(
            make_signature(
                &mut heap,
                later_allocation_failure,
                8,
                &[0, 1, 2, 3, 4, 5, 6, 7],
                &mut later_failure_counts,
            ),
            Err(SourceHeapError::AllocationFailed)
        );
        assert_eq!(heap.source_allocation_calls(), 2);
        assert_eq!(later_failure_counts, [1, 1, 1, 1, 1, 1, 1, 1, 0]);
        assert_eq!(
            signature(&heap, later_allocation_failure),
            b"8,0,7{(0:1)(1:1)(2:1)(3:1)(4:1)(5:1)"
        );
        release(&mut heap, later_allocation_failure);
        assert_eq!(heap.live_allocation_count(), 0);
    }

    #[test]
    fn source_port__runichi2__diylfrag_diff__line_2098() {
        fn fragment(pointer: SourceMutPointer<i8>, used: i32, na: i32, nb: i32) -> DiylFrag {
            DiylFrag {
                na,
                nb,
                sig: INCHI_IOS_STRING {
                    pStr: pointer,
                    nAllocatedLength: 91,
                    nUsedLength: used,
                    nPtr: 92,
                },
                ..DiylFrag::default()
            }
        }

        let mut heap = SourceHeap::default();
        let null = SourceMutPointer::null();
        assert_eq!(
            DiylFrag_Diff(&heap, &fragment(null, 1, i32::MIN, 7), &fragment(null, 1, i32::MAX, 8),),
            Ok(1)
        );
        assert_eq!(
            DiylFrag_Diff(&heap, &fragment(null, 1, 4, i32::MIN), &fragment(null, 1, 4, i32::MAX),),
            Ok(1)
        );

        for (first_used, second_used) in [(0, 0), (0, 1), (1, 0), (0, -1), (-1, 0)] {
            assert_eq!(
                DiylFrag_Diff(
                    &heap,
                    &fragment(null, first_used, 4, 3),
                    &fragment(null, second_used, 4, 3),
                ),
                Ok(0),
                "first_used={first_used}, second_used={second_used}"
            );
        }

        let equal_left = c_string(&mut heap, b"same");
        let equal_right = c_string(&mut heap, b"same");
        assert_eq!(
            DiylFrag_Diff(&heap, &fragment(equal_left, 1, 4, 3), &fragment(equal_right, 99, 4, 3),),
            Ok(0)
        );

        let lower = c_string(&mut heap, b"az");
        let upper = c_string(&mut heap, b"aZ");
        assert_eq!(
            DiylFrag_Diff(&heap, &fragment(lower, 1, 4, 3), &fragment(upper, 1, 4, 3),),
            Ok(32)
        );
        assert_eq!(
            DiylFrag_Diff(&heap, &fragment(upper, -1, 4, 3), &fragment(lower, i32::MIN, 4, 3),),
            Ok(-32)
        );

        let prefix = c_string(&mut heap, b"a");
        let longer = c_string(&mut heap, b"aa");
        assert_eq!(
            DiylFrag_Diff(&heap, &fragment(prefix, 77, 4, 3), &fragment(longer, 1, 4, 3),),
            Ok(-97)
        );
        assert_eq!(
            DiylFrag_Diff(&heap, &fragment(longer, 1, 4, 3), &fragment(prefix, 77, 4, 3),),
            Ok(97)
        );

        let high = heap.allocate_model_storage(vec![-1_i8, 0]).unwrap();
        let low = heap.allocate_model_storage(vec![1_i8, 0]).unwrap();
        assert_eq!(
            DiylFrag_Diff(&heap, &fragment(high, 1, 4, 3), &fragment(low, 1, 4, 3),),
            Ok(254)
        );

        assert_eq!(
            DiylFrag_Diff(&heap, &fragment(null, 1, 4, 3), &fragment(equal_right, 1, 4, 3),),
            Err(SourceHeapError::NullPointer)
        );
        let unterminated_left = heap.allocate_model_storage(vec![b'x' as i8]).unwrap();
        let unterminated_right = heap.allocate_model_storage(vec![b'x' as i8]).unwrap();
        assert_eq!(
            DiylFrag_Diff(
                &heap,
                &fragment(unterminated_left, 1, 4, 3),
                &fragment(unterminated_right, 1, 4, 3),
            ),
            Err(SourceHeapError::MissingNulTerminator)
        );
    }

    #[test]
    fn source_port__runichi2__diylfrag_debugtrace__line_2119() {
        assert_eq!(DiylFrag_DebugTrace(None), Ok(()));

        let dangling = SourceMutPointer::<i32>::null();
        let dangling_string = SourceMutPointer::<i8>::null();
        for na in [-7, -1, 0, 1, 2, 5, 32_767] {
            let fragment = DiylFrag {
                na,
                nb: i32::MIN,
                end1: i32::MIN,
                end2: i32::MAX,
                alist: dangling,
                xclist: dangling,
                sig: INCHI_IOS_STRING {
                    pStr: dangling_string,
                    nAllocatedLength: i32::MIN,
                    nUsedLength: i32::MAX,
                    nPtr: i32::MIN,
                },
            };
            let before = fragment.clone();
            assert_eq!(DiylFrag_DebugTrace(Some(&fragment)), Ok(()), "na={na}");
            assert_eq!(fragment, before, "na={na}");
        }

        let minimum = DiylFrag {
            na: i32::MIN,
            alist: dangling,
            xclist: dangling,
            sig: INCHI_IOS_STRING {
                pStr: dangling_string,
                ..INCHI_IOS_STRING::default()
            },
            ..DiylFrag::default()
        };
        assert_eq!(
            DiylFrag_DebugTrace(Some(&minimum)),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
    }

    #[test]
    fn source_port__runichi2__oad_structureedits_clear__line_1769() {
        fn allocated_array(
            heap: &mut SourceHeap,
            values: &[i32],
        ) -> (SourceMutPointer<INT_ARRAY>, SourceMutPointer<i32>) {
            let item = if values.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(values.to_vec()).unwrap()
            };
            let array = heap
                .allocate_model_storage(vec![INT_ARRAY {
                    item,
                    used: values.len() as i32,
                    allocated: values.len() as i32,
                    increment: 7,
                }])
                .unwrap();
            (array, item)
        }

        let mut heap = SourceHeap::default();
        let (del_atom, del_atom_item) = allocated_array(&mut heap, &[1, 2]);
        let (del_bond, del_bond_item) = allocated_array(&mut heap, &[]);
        let (new_bond, new_bond_item) = allocated_array(&mut heap, &[3, 4]);
        let (mod_bond, mod_bond_item) = allocated_array(&mut heap, &[5; 12]);
        let (mod_coord, mod_coord_item) = allocated_array(&mut heap, &[6; 4]);
        let mut edits = OAD_StructureEdits {
            del_atom,
            del_bond,
            new_bond,
            mod_bond,
            mod_coord,
            del_side_chains: i32::MIN,
        };

        assert_eq!(OAD_StructureEdits_Clear(&mut heap, &mut edits), Ok(()));
        assert_eq!(
            edits,
            OAD_StructureEdits {
                del_side_chains: i32::MIN,
                ..OAD_StructureEdits::default()
            }
        );
        for pointer in [del_atom, del_bond, new_bond, mod_bond, mod_coord] {
            assert!(heap.slice(pointer.as_const()).is_err());
        }
        for pointer in [del_atom_item, new_bond_item, mod_bond_item, mod_coord_item] {
            assert!(heap.slice(pointer.as_const()).is_err());
        }
        assert!(del_bond_item.is_null());

        assert_eq!(OAD_StructureEdits_Clear(&mut heap, &mut edits), Ok(()));
        assert_eq!(edits.del_side_chains, i32::MIN);

        let mut sparse_heap = SourceHeap::default();
        let (only_new_bond, only_new_bond_item) = allocated_array(&mut sparse_heap, &[9]);
        let mut sparse = OAD_StructureEdits {
            new_bond: only_new_bond,
            del_side_chains: i32::MAX,
            ..OAD_StructureEdits::default()
        };
        assert_eq!(OAD_StructureEdits_Clear(&mut sparse_heap, &mut sparse), Ok(()));
        assert!(sparse.new_bond.is_null());
        assert_eq!(sparse.del_side_chains, i32::MAX);
        assert!(sparse_heap.slice(only_new_bond.as_const()).is_err());
        assert!(sparse_heap.slice(only_new_bond_item.as_const()).is_err());
    }

    #[test]
    fn source_port__runichi2__oad_structureedits_init__line_1735() {
        fn assert_array(heap: &SourceHeap, pointer: SourceMutPointer<INT_ARRAY>, capacity: i32) {
            let array = &heap.slice(pointer.as_const()).unwrap()[0];
            assert_eq!(array.used, 0);
            assert_eq!(array.allocated, capacity);
            assert_eq!(array.increment, capacity);
            assert_eq!(
                heap.slice(array.item.as_const()).unwrap(),
                vec![0_i32; capacity as usize]
            );
        }

        let mut heap = SourceHeap::default();
        heap.trace_source_allocations();
        let mut edits = OAD_StructureEdits {
            del_side_chains: i32::MAX,
            ..OAD_StructureEdits::default()
        };
        assert_eq!(OAD_StructureEdits_Init(&mut heap, &mut edits), Ok(0));
        assert_eq!(heap.source_allocation_calls(), 10);
        assert_eq!(heap.live_allocation_count(), 10);
        assert_eq!(edits.del_side_chains, 0);
        assert_array(&heap, edits.del_atom, 2);
        assert_array(&heap, edits.del_bond, 2);
        assert_array(&heap, edits.new_bond, 2);
        assert_array(&heap, edits.mod_bond, 12);
        assert_array(&heap, edits.mod_coord, 4);
        assert_eq!(OAD_StructureEdits_Clear(&mut heap, &mut edits), Ok(()));
        assert_eq!(heap.live_allocation_count(), 0);

        for successful_allocations in 0..10 {
            let mut failure_heap = SourceHeap::default();
            failure_heap.fail_after_allocations(successful_allocations);
            let mut failed = OAD_StructureEdits {
                del_side_chains: -1,
                ..OAD_StructureEdits::default()
            };
            assert_eq!(
                OAD_StructureEdits_Init(&mut failure_heap, &mut failed),
                Ok(_IS_ERROR as i32),
                "successful_allocations={successful_allocations}"
            );
            assert_eq!(failed, OAD_StructureEdits::default());
            assert_eq!(
                failure_heap.live_allocation_count(),
                0,
                "successful_allocations={successful_allocations}"
            );
        }
    }

    #[test]
    fn source_port__runichi2__oad_structureedits_debugprint__line_1807() {
        OAD_StructureEdits_DebugPrint(&OAD_StructureEdits::default());

        let mut heap = SourceHeap::default();
        let pointers: [SourceMutPointer<INT_ARRAY>; 5] = std::array::from_fn(|index| {
            heap.allocate_model_storage(vec![INT_ARRAY {
                item: SourceMutPointer::null(),
                allocated: index as i32,
                used: if index == 4 { i32::MAX } else { index as i32 },
                increment: -(index as i32),
            }])
            .unwrap()
        });
        let edits = OAD_StructureEdits {
            del_atom: pointers[0],
            del_bond: pointers[1],
            new_bond: pointers[2],
            mod_bond: pointers[3],
            mod_coord: pointers[4],
            del_side_chains: i32::MIN,
        };
        let before = pointers.map(|pointer| heap.slice(pointer.as_const()).unwrap()[0].clone());
        OAD_StructureEdits_DebugPrint(&edits);
        assert_eq!(edits.del_side_chains, i32::MIN);
        for (pointer, expected) in pointers.into_iter().zip(before) {
            assert_eq!(heap.slice(pointer.as_const()).unwrap()[0], expected);
        }
    }

    #[test]
    fn source_port__runichi2__posecontext_free__line_1636() {
        let mut heap = SourceHeap::default();
        let external_input = heap.allocate_model_storage(vec![INCHI_IOSTREAM::default()]).unwrap();
        let path0 = c_string(&mut heap, b"first");
        let path3 = c_string(&mut heap, b"");
        let parameters = INPUT_PARMS {
            path: [
                path0.as_const(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                path3.as_const(),
            ],
            ..INPUT_PARMS::default()
        };
        let original_atoms = heap.allocate_model_storage(vec![early_carbon()]).unwrap();
        let original = ORIG_ATOM_DATA {
            at: original_atoms,
            num_inp_atoms: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let prepared_atoms = heap.allocate_model_storage(vec![early_carbon()]).unwrap();
        let prepared = [
            ORIG_ATOM_DATA {
                at: prepared_atoms,
                num_inp_atoms: 1,
                ..ORIG_ATOM_DATA::default()
            },
            ORIG_ATOM_DATA::default(),
        ];
        let source_string_allocation = heap.allocate_model_storage(vec![b'X' as i8, 0, 9, 9]).unwrap();
        let source_string = INCHI_IOS_STRING {
            pStr: source_string_allocation,
            nAllocatedLength: 4,
            nUsedLength: 1,
            nPtr: 2,
        };
        let mut context = POSEContext::default();
        assert_eq!(
            POSEContext_Init(
                &mut heap,
                &mut context,
                Some(&STRUCT_DATA::default()),
                Some(&parameters),
                &[b'T' as i8, 0],
                None,
                None,
                external_input,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                Some(&original),
                Some(&prepared),
                i64::MAX,
                Some(&source_string),
                u8::MAX,
            ),
            Ok(_IS_OKAY as i32)
        );

        let copied_paths = [context.ip.path[0], context.ip.path[3]];
        assert_ne!(copied_paths[0], path0.as_const());
        assert_ne!(copied_paths[1], path3.as_const());

        let inchi = heap.allocate(vec![INChI::default()]).unwrap();
        let inchi_rows = heap.allocate(vec![[inchi, SourceMutPointer::null()]]).unwrap();
        let aux = heap.allocate(vec![INChI_Aux::default()]).unwrap();
        let aux_rows = heap.allocate(vec![[aux, SourceMutPointer::null()]]).unwrap();
        let retained_zero_inchi = heap
            .allocate(vec![[SourceMutPointer::null(), SourceMutPointer::null()]])
            .unwrap();
        let retained_zero_aux = heap
            .allocate(vec![[SourceMutPointer::null(), SourceMutPointer::null()]])
            .unwrap();
        context.pINChI2 = [inchi_rows, retained_zero_inchi];
        context.pINChI_Aux2 = [aux_rows, retained_zero_aux];
        context.sd.num_components = [1, 0];

        let mut stream_strings = Vec::new();
        let mut stream_files = Vec::new();
        for index in 0..3_i64 {
            let string = heap.allocate(vec![b'A' as i8 + index as i8, 0]).unwrap();
            let file = heap
                .allocate(vec![SourceFile {
                    is_standard_stream: false,
                    ..SourceFile::default()
                }])
                .unwrap();
            let stream_pointer = context.out_file.offset(index).unwrap();
            let stream = heap.slice_mut(stream_pointer).unwrap().first_mut().unwrap();
            stream.s = INCHI_IOS_STRING {
                pStr: string,
                nAllocatedLength: 2,
                nUsedLength: 1,
                nPtr: 1,
            };
            stream.f = file;
            stream_strings.push(string);
            stream_files.push(file);
        }
        let copied_original_atoms = heap.slice(context.orig_inp_data.as_const()).unwrap()[0].at;
        let copied_prepared_atoms = heap.slice(context.prep_inp_data.as_const()).unwrap()[0].at;
        let copied_string = heap.slice(context.strbuf.as_const()).unwrap()[0].pStr;

        assert_eq!(POSEContext_Free(&mut heap, &mut context), Ok(()));

        assert_eq!(context.ip.path, [SourceConstPointer::null(); 4]);
        for path in copied_paths {
            assert_eq!(heap.slice(path), Err(SourceHeapError::MissingAllocation));
        }
        assert_eq!(context.sd.num_components, [0, 0]);
        assert_eq!(context.pINChI2[0], SourceMutPointer::null());
        assert_eq!(context.pINChI_Aux2[0], SourceMutPointer::null());
        assert_eq!(context.pINChI2[1], retained_zero_inchi);
        assert_eq!(context.pINChI_Aux2[1], retained_zero_aux);
        for released in [
            inchi.cast::<u8>(),
            aux.cast::<u8>(),
            inchi_rows.cast::<u8>(),
            aux_rows.cast::<u8>(),
        ] {
            assert_eq!(heap.slice(released.as_const()), Err(SourceHeapError::MissingAllocation));
        }
        assert!(heap.slice(retained_zero_inchi.as_const()).is_ok());
        assert!(heap.slice(retained_zero_aux.as_const()).is_ok());

        let closed_streams = &heap.slice(context.out_file.as_const()).unwrap()[..3];
        for (index, stream) in closed_streams.iter().enumerate() {
            assert_eq!(stream.s, INCHI_IOS_STRING::default());
            assert_eq!(stream.f, stream_files[index]);
            assert_eq!(
                heap.slice(stream_strings[index].as_const()),
                Err(SourceHeapError::MissingAllocation)
            );
            assert_eq!(
                heap.slice(stream_files[index].as_const()),
                Err(SourceHeapError::MissingAllocation)
            );
        }
        assert_eq!(context.inchi_file, closed_streams);
        assert_eq!(context.inp_file, external_input);
        assert!(heap.slice(external_input.as_const()).is_ok());

        assert_eq!(context.OrigAtData, ORIG_ATOM_DATA::default());
        assert_eq!(context.PrepAtData, std::array::from_fn(|_| ORIG_ATOM_DATA::default()));
        assert_eq!(
            heap.slice(copied_original_atoms.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(copied_prepared_atoms.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(context.num_inp, 0);
        assert_eq!(context.save_opt_bits, 0);
        assert_eq!(context.temp_string_container, INCHI_IOS_STRING::default());
        assert_eq!(
            heap.slice(copied_string.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(source_string_allocation.as_const()).unwrap(),
            &[b'X' as i8, 0, 9, 9]
        );

        let mut empty = POSEContext {
            num_inp: i64::MIN,
            save_opt_bits: 17,
            ..POSEContext::default()
        };
        assert_eq!(POSEContext_Free(&mut SourceHeap::default(), &mut empty), Ok(()));
        assert_eq!(empty.num_inp, 0);
        assert_eq!(empty.save_opt_bits, 0);

        let mut partial_heap = SourceHeap::default();
        let aliased_path = c_string(&mut partial_heap, b"aliased-on-failure");
        let partial_parameters = INPUT_PARMS {
            path: [
                aliased_path.as_const(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
            ],
            ..INPUT_PARMS::default()
        };
        partial_heap.fail_after_allocations(0);
        let mut partial = POSEContext::default();
        assert_eq!(
            POSEContext_Init(
                &mut partial_heap,
                &mut partial,
                None,
                Some(&partial_parameters),
                &[0],
                None,
                None,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                None,
                None,
                5,
                None,
                6,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(partial.ip.path[0], aliased_path.as_const());
        assert_eq!(POSEContext_Free(&mut partial_heap, &mut partial), Ok(()));
        assert!(partial.ip.path[0].is_null());
        assert_eq!(
            partial_heap.slice(aliased_path.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let mut string_failure_heap = SourceHeap::default();
        string_failure_heap.fail_after_allocations(0);
        let mut string_failure = POSEContext::default();
        assert_eq!(
            POSEContext_Init(
                &mut string_failure_heap,
                &mut string_failure,
                None,
                None,
                &[0],
                None,
                None,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                None,
                None,
                -7,
                None,
                8,
            ),
            Ok(_IS_FATAL as i32)
        );
        assert_eq!(POSEContext_Free(&mut string_failure_heap, &mut string_failure), Ok(()));
        assert_eq!(string_failure.num_inp, 0);
        assert_eq!(string_failure.save_opt_bits, 0);
        assert_eq!(string_failure.temp_string_container, INCHI_IOS_STRING::default());
    }

    #[test]
    fn source_port__runichi2__posecontext_init__line_1503() {
        fn invoke(
            heap: &mut SourceHeap,
            context: &mut POSEContext,
            structure: Option<&STRUCT_DATA>,
            parameters: Option<&INPUT_PARMS>,
            title: &[i8],
            inchi: Option<&[SourceMutPointer<PINChI2>; INCHI_NUM as usize]>,
            aux: Option<&[SourceMutPointer<PINChI_Aux2>; INCHI_NUM as usize]>,
            input: SourceMutPointer<INCHI_IOSTREAM>,
            original: Option<&ORIG_ATOM_DATA>,
            prepared: Option<&[ORIG_ATOM_DATA]>,
            string: Option<&INCHI_IOS_STRING>,
        ) -> Result<i32, SourceHeapError> {
            POSEContext_Init(
                heap,
                context,
                structure,
                parameters,
                title,
                inchi,
                aux,
                input,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                original,
                prepared,
                i64::MIN,
                string,
                0xa5,
            )
        }

        let mut heap = SourceHeap::default();
        let external_input = heap.allocate_model_storage(vec![INCHI_IOSTREAM::default()]).unwrap();
        let mut empty_context = POSEContext {
            num_inp: 77,
            save_opt_bits: 88,
            ..POSEContext::default()
        };
        assert_eq!(
            invoke(
                &mut heap,
                &mut empty_context,
                None,
                None,
                &[0],
                None,
                None,
                external_input,
                None,
                None,
                None,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(empty_context.sd, STRUCT_DATA::default());
        assert_eq!(empty_context.ip, INPUT_PARMS::default());
        assert_eq!(empty_context.szTitle, [0; 575]);
        assert_eq!(empty_context.inp_file, external_input);
        assert_eq!(empty_context.num_inp, i64::MIN);
        assert_eq!(empty_context.save_opt_bits, 0xa5);
        assert!(!empty_context.out_file.is_null());
        assert_eq!(empty_context.log_file.difference(empty_context.out_file), Ok(1));
        assert_eq!(empty_context.prb_file.difference(empty_context.out_file), Ok(2));
        for stream in heap.slice(empty_context.out_file.as_const()).unwrap().iter().take(3) {
            assert_eq!(stream.type_, INCHI_IOS_TYPE_STRING as i32);
            assert!(stream.f.is_null());
            assert_eq!(stream.s, INCHI_IOS_STRING::default());
        }
        assert_eq!(
            heap.slice(empty_context.orig_inp_data.as_const()).unwrap()[0],
            ORIG_ATOM_DATA::default()
        );
        assert_eq!(
            &heap.slice(empty_context.prep_inp_data.as_const()).unwrap()[..2],
            &[ORIG_ATOM_DATA::default(), ORIG_ATOM_DATA::default()]
        );
        let initialized_string = &heap.slice(empty_context.strbuf.as_const()).unwrap()[0];
        assert_eq!(initialized_string.nAllocatedLength, INCHI_STRBUF_INITIAL_SIZE as i32);
        assert_eq!(initialized_string.nPtr, INCHI_STRBUF_SIZE_INCREMENT as i32);
        assert_eq!(initialized_string.nUsedLength, 0);
        assert_eq!(
            heap.slice(initialized_string.pStr.as_const()).unwrap(),
            vec![0; INCHI_STRBUF_INITIAL_SIZE as usize]
        );

        let path0 = c_string(&mut heap, b"alpha");
        let path2 = c_string(&mut heap, b"");
        let structure = STRUCT_DATA {
            nErrorCode: -93,
            num_components: [2, 3],
            ..STRUCT_DATA::default()
        };
        let parameters = INPUT_PARMS {
            path: [
                path0.as_const(),
                SourceConstPointer::null(),
                path2.as_const(),
                SourceConstPointer::null(),
            ],
            num_paths: 3,
            bNoWarnings: 41,
            ..INPUT_PARMS::default()
        };
        let source_atoms = heap.allocate_model_storage(vec![early_carbon()]).unwrap();
        let original = ORIG_ATOM_DATA {
            at: source_atoms,
            num_inp_atoms: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let prepared_atoms = heap.allocate_model_storage(vec![early_carbon()]).unwrap();
        let prepared = [
            ORIG_ATOM_DATA {
                at: prepared_atoms,
                num_inp_atoms: 1,
                ..ORIG_ATOM_DATA::default()
            },
            ORIG_ATOM_DATA {
                num_inp_atoms: 99,
                ..ORIG_ATOM_DATA::default()
            },
        ];
        let source_string_pointer = heap.allocate_model_storage(vec![b'Q' as i8, 0, 7, 8]).unwrap();
        let source_string = INCHI_IOS_STRING {
            pStr: source_string_pointer,
            nAllocatedLength: 4,
            nUsedLength: 1,
            nPtr: 13,
        };
        let mut rich = POSEContext::default();
        assert_eq!(
            invoke(
                &mut heap,
                &mut rich,
                Some(&structure),
                Some(&parameters),
                &[b'T' as i8, b'i' as i8, b't' as i8, b'l' as i8, b'e' as i8, 0],
                Some(&[SourceMutPointer::null(); INCHI_NUM as usize]),
                Some(&[SourceMutPointer::null(); INCHI_NUM as usize]),
                external_input,
                Some(&original),
                Some(&prepared),
                Some(&source_string),
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(rich.sd, structure);
        assert_eq!(
            &rich.szTitle[..6],
            &[b'T' as i8, b'i' as i8, b't' as i8, b'l' as i8, b'e' as i8, 0]
        );
        assert_ne!(rich.ip.path[0], path0.as_const());
        assert_ne!(rich.ip.path[2], path2.as_const());
        assert_eq!(source_c_bytes(&heap, rich.ip.path[0].as_mut()).unwrap(), b"alpha");
        assert_eq!(source_c_bytes(&heap, rich.ip.path[2].as_mut()).unwrap(), b"");
        let copied_original = &heap.slice(rich.orig_inp_data.as_const()).unwrap()[0];
        assert_ne!(copied_original.at, source_atoms);
        assert_eq!(copied_original.num_inp_atoms, 1);
        assert_eq!(heap.slice(copied_original.at.as_const()).unwrap()[0], early_carbon());
        let copied_prepared = heap.slice(rich.prep_inp_data.as_const()).unwrap();
        assert_eq!(copied_prepared[0].num_inp_atoms, 1);
        assert_ne!(copied_prepared[0].at, prepared_atoms);
        assert_eq!(copied_prepared[1], ORIG_ATOM_DATA::default());
        let copied_string = &heap.slice(rich.strbuf.as_const()).unwrap()[0];
        assert_eq!(
            (
                copied_string.nAllocatedLength,
                copied_string.nUsedLength,
                copied_string.nPtr,
            ),
            (4, 1, 13)
        );
        assert_eq!(heap.slice(copied_string.pStr.as_const()).unwrap(), &[0; 4]);
        assert_eq!(
            heap.slice(source_string_pointer.as_const()).unwrap(),
            &[b'Q' as i8, 0, 7, 8]
        );

        for reject_aux in [false, true] {
            let mut reject_heap = SourceHeap::default();
            let nonnull = reject_heap
                .allocate_model_storage(vec![if reject_aux {
                    PINChI_Aux2::default()
                } else {
                    PINChI_Aux2::default()
                }])
                .unwrap();
            let mut rejected = POSEContext::default();
            let inchi_values = if reject_aux {
                [SourceMutPointer::null(); INCHI_NUM as usize]
            } else {
                [nonnull.cast::<PINChI2>(), SourceMutPointer::null()]
            };
            let aux_values = if reject_aux {
                [nonnull, SourceMutPointer::null()]
            } else {
                [SourceMutPointer::null(); INCHI_NUM as usize]
            };
            assert_eq!(
                invoke(
                    &mut reject_heap,
                    &mut rejected,
                    None,
                    None,
                    &[0],
                    Some(&inchi_values),
                    Some(&aux_values),
                    SourceMutPointer::null(),
                    None,
                    None,
                    None,
                ),
                Ok(_IS_ERROR as i32)
            );
            assert!(rejected.out_file.is_null());
        }

        let mut path_failure_heap = SourceHeap::default();
        let source_path = c_string(&mut path_failure_heap, b"owned-later");
        let path_parameters = INPUT_PARMS {
            path: [
                source_path.as_const(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
            ],
            ..INPUT_PARMS::default()
        };
        path_failure_heap.fail_after_allocations(0);
        let mut path_failure = POSEContext::default();
        assert_eq!(
            invoke(
                &mut path_failure_heap,
                &mut path_failure,
                None,
                Some(&path_parameters),
                &[0],
                None,
                None,
                SourceMutPointer::null(),
                None,
                None,
                None,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(path_failure.ip.path[0], source_path.as_const());
        assert!(path_failure.out_file.is_null());

        let mut duplicate_failure_heap = SourceHeap::default();
        duplicate_failure_heap.fail_after_allocations(0);
        let mut duplicate_failure = POSEContext::default();
        assert_eq!(
            invoke(
                &mut duplicate_failure_heap,
                &mut duplicate_failure,
                None,
                None,
                &[0],
                None,
                None,
                SourceMutPointer::null(),
                Some(&ORIG_ATOM_DATA::default()),
                None,
                None,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert!(!duplicate_failure.orig_inp_data.is_null());
        assert!(duplicate_failure.prep_inp_data.is_null());

        let mut string_failure_heap = SourceHeap::default();
        string_failure_heap.fail_after_allocations(0);
        let mut string_failure = POSEContext::default();
        assert_eq!(
            invoke(
                &mut string_failure_heap,
                &mut string_failure,
                None,
                None,
                &[0],
                None,
                None,
                SourceMutPointer::null(),
                None,
                None,
                None,
            ),
            Ok(_IS_FATAL as i32)
        );
        assert!(!string_failure.strbuf.is_null());
        assert!(
            string_failure_heap.slice(string_failure.strbuf.as_const()).unwrap()[0]
                .pStr
                .is_null()
        );

        let mut malformed_context = POSEContext::default();
        assert_eq!(
            invoke(
                &mut SourceHeap::default(),
                &mut malformed_context,
                None,
                None,
                &[b'X' as i8],
                None,
                None,
                SourceMutPointer::null(),
                None,
                None,
                None,
            ),
            Err(SourceHeapError::MissingNulTerminator)
        );
        let oversized_title = vec![b'X' as i8; 576];
        let mut oversized_context = POSEContext::default();
        assert_eq!(
            invoke(
                &mut SourceHeap::default(),
                &mut oversized_context,
                None,
                None,
                &oversized_title,
                None,
                None,
                SourceMutPointer::null(),
                None,
                None,
                None,
            ),
            Err(SourceHeapError::MissingNulTerminator)
        );
        let mut terminated_oversized = vec![b'X' as i8; 576];
        terminated_oversized[575] = 0;
        assert_eq!(
            invoke(
                &mut SourceHeap::default(),
                &mut oversized_context,
                None,
                None,
                &terminated_oversized,
                None,
                None,
                SourceMutPointer::null(),
                None,
                None,
                None,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__runichi2__bissamebond__line_1293() {
        assert_eq!(bIsSameBond(1, 2, 1, 2), 1);
        assert_eq!(bIsSameBond(1, 2, 2, 1), -1);
        assert_eq!(bIsSameBond(1, 2, 1, 3), 0);
        assert_eq!(bIsSameBond(7, 7, 7, 7), 1);
        assert_eq!(bIsSameBond(i32::MIN, i32::MAX, i32::MIN, i32::MAX), 1);
        assert_eq!(bIsSameBond(i32::MIN, i32::MAX, i32::MAX, i32::MIN), -1);
        assert_eq!(bIsSameBond(i32::MIN, i32::MIN, i32::MAX, i32::MAX), 0);
    }

    #[test]
    fn source_port__runichi2__getframeshiftinfofrom105plusinchi__line_1305() {
        let mut heap = SourceHeap::default();

        for input in [b"" as &[u8], b"InChI=1S/CH4/h1H4", b"/x(1,2,3)"] {
            let input = c_string(&mut heap, input);
            let mut output = [91_i32; 6];
            assert_eq!(
                GetFrameShiftInfoFrom105PlusInChI(&mut heap, input.as_const(), &mut output, 2,),
                Ok(0)
            );
            assert_eq!(output, [91; 6]);
        }

        let all_branches = c_string(
            &mut heap,
            b"InChI=1B/x/zA(1,2,3);B(4.5.6);C(7);D(8-9,10);E(0,11,12);F(-1,+2,9);G(2147483648,2147483647,1);H(4294967296,5,1)",
        );
        let mut output = [91_i32; 15];
        assert_eq!(
            GetFrameShiftInfoFrom105PlusInChI(&mut heap, all_branches.as_const(), &mut output, 5,),
            Ok(4)
        );
        assert_eq!(output, [0, 1, 2, 1, 4, 5, 5, -1, 2, 6, i32::MIN, i32::MAX, 91, 91, 91,]);

        let closed_second = c_string(&mut heap, b"/zA(1,2);B(3,4,5)");
        let mut closed_output = [73_i32; 6];
        assert_eq!(
            GetFrameShiftInfoFrom105PlusInChI(&mut heap, closed_second.as_const(), &mut closed_output, 2,),
            Ok(1)
        );
        assert_eq!(closed_output, [1, 3, 4, 73, 73, 73]);

        let overflow_third = c_string(&mut heap, b"/zA(1,2,9223372036854775808)");
        let mut overflow_output = [17_i32; 3];
        heap.set_source_errno(0);
        assert_eq!(
            GetFrameShiftInfoFrom105PlusInChI(&mut heap, overflow_third.as_const(), &mut overflow_output, 1,),
            Ok(1)
        );
        assert_eq!(overflow_output, [0, 1, 2]);
        assert_eq!(heap.source_errno(), 34);

        for maximum in [1, 0, -1] {
            let input = c_string(&mut heap, b"/zA(10,11,12);B(20,21,22)");
            let mut limited = [63_i32; 6];
            assert_eq!(
                GetFrameShiftInfoFrom105PlusInChI(&mut heap, input.as_const(), &mut limited, maximum,),
                Ok(1),
                "maximum={maximum}"
            );
            assert_eq!(limited, [0, 10, 11, 63, 63, 63]);
        }

        let interior_allocation = c_string(&mut heap, b"junk/zA(+3,-4,5)");
        let interior = interior_allocation.as_const().offset(4).unwrap();
        let mut interior_output = [29_i32; 3];
        assert_eq!(
            GetFrameShiftInfoFrom105PlusInChI(&mut heap, interior, &mut interior_output, 1),
            Ok(1)
        );
        assert_eq!(interior_output, [0, 3, -4]);

        let partial_input = c_string(&mut heap, b"/zA(7,8,9)");
        let mut one_slot = [55_i32; 1];
        assert_eq!(
            GetFrameShiftInfoFrom105PlusInChI(&mut heap, partial_input.as_const(), &mut one_slot, 1,),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(one_slot, [0]);
        let mut two_slots = [55_i32; 2];
        assert_eq!(
            GetFrameShiftInfoFrom105PlusInChI(&mut heap, partial_input.as_const(), &mut two_slots, 1,),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(two_slots, [0, 7]);

        let unterminated = heap.allocate_model_storage(vec![b'/' as i8, b'z' as i8]).unwrap();
        assert_eq!(
            GetFrameShiftInfoFrom105PlusInChI(&mut heap, unterminated.as_const(), &mut [0_i32; 3], 1,),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(
            GetFrameShiftInfoFrom105PlusInChI(&mut heap, SourceConstPointer::null(), &mut [0_i32; 3], 1,),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__runichi2__oad_polymer_prepareframeshiftedits__line_2523() {
        const IDENTITY_AUX: &[u8] = b"AuxInfo=1/0/N:1,2,3,4,5,6,7,8";
        const SHIFT_3_6: &[u8] = b"InChI=1B/C6H12.2Zz/zA(3,6,1)";

        #[derive(Clone)]
        struct Fixture {
            data: ORIG_ATOM_DATA,
            units: Vec<SourceMutPointer<OAD_PolymerUnit>>,
            edits: OAD_StructureEdits,
        }

        fn named_atom(name: &[u8]) -> inp_ATOM {
            let mut atom = inp_ATOM::default();
            for (target, source) in atom.elname.iter_mut().zip(name.iter().copied()) {
                *target = source as i8;
            }
            atom
        }

        fn fixture(
            heap: &mut SourceHeap,
            unit_specs: Vec<(Vec<i32>, Option<Vec<i32>>, i32)>,
            star_caps: bool,
            edit_capacity: i32,
        ) -> Fixture {
            let mut atoms = vec![named_atom(b"C"); 8];
            if star_caps {
                atoms[0] = named_atom(b"Zz");
                atoms[7] = named_atom(b"Zz");
            }
            for (index, atom) in atoms.iter_mut().enumerate() {
                atom.orig_at_number = (index + 1) as u16;
            }
            let atoms = heap.allocate_model_storage(atoms).unwrap();

            let mut units = Vec::new();
            for (atom_list, crossing_bonds, crossing_count) in unit_specs {
                let atom_count = atom_list.len() as i32;
                let atom_list = heap.allocate_model_storage(atom_list).unwrap();
                let crossing_bonds = crossing_bonds
                    .map(|values| heap.allocate_model_storage(values).unwrap())
                    .unwrap_or_else(SourceMutPointer::null);
                units.push(
                    heap.allocate_model_storage(vec![OAD_PolymerUnit {
                        na: atom_count,
                        nb: crossing_count,
                        alist: atom_list,
                        blist: crossing_bonds,
                        ..OAD_PolymerUnit::default()
                    }])
                    .unwrap(),
                );
            }
            let unit_table = heap.allocate_model_storage(units.clone()).unwrap();
            let polymer = heap
                .allocate_model_storage(vec![OAD_Polymer {
                    units: unit_table,
                    n: units.len() as i32,
                    ..OAD_Polymer::default()
                }])
                .unwrap();

            let edit_item = heap
                .allocate_model_storage(vec![-91_i32; edit_capacity as usize])
                .unwrap();
            let mod_bond = heap
                .allocate_model_storage(vec![INT_ARRAY {
                    item: edit_item,
                    allocated: edit_capacity,
                    used: 0,
                    increment: edit_capacity.max(1),
                }])
                .unwrap();
            Fixture {
                data: ORIG_ATOM_DATA {
                    at: atoms,
                    num_inp_atoms: 8,
                    polymer,
                    ..ORIG_ATOM_DATA::default()
                },
                units,
                edits: OAD_StructureEdits {
                    mod_bond,
                    ..OAD_StructureEdits::default()
                },
            }
        }

        fn invoke(heap: &mut SourceHeap, fixture: &Fixture, inchi: &[u8], aux: &[u8]) -> Result<i32, SourceHeapError> {
            let inchi = c_string(heap, inchi);
            let aux = c_string(heap, aux);
            heap.trace_source_allocations();
            OAD_Polymer_PrepareFrameShiftEdits(heap, &fixture.data, inchi.as_const(), aux.as_const(), &fixture.edits)
        }

        fn edit_values(heap: &SourceHeap, edits: &OAD_StructureEdits) -> Vec<i32> {
            let array = &heap.slice(edits.mod_bond.as_const()).unwrap()[0];
            heap.slice(array.item.as_const()).unwrap()[..array.used as usize].to_vec()
        }

        let standard_unit = || vec![(vec![2, 3, 4, 5, 6, 7], Some(vec![1, 2, 8, 7]), 2)];

        let mut empty_heap = SourceHeap::default();
        let empty = fixture(&mut empty_heap, standard_unit(), true, 16);
        let empty_baseline = empty_heap.live_allocation_count();
        assert_eq!(
            invoke(&mut empty_heap, &empty, b"InChI=1B/C6H12.2Zz", IDENTITY_AUX),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(empty_heap.source_allocation_calls(), 3);
        assert_eq!(empty_heap.live_allocation_count(), empty_baseline + 2);
        assert!(edit_values(&empty_heap, &empty.edits).is_empty());

        let mut malformed_heap = SourceHeap::default();
        let malformed = fixture(&mut malformed_heap, standard_unit(), true, 16);
        let malformed_inchi = c_string(&mut malformed_heap, SHIFT_3_6);
        let malformed_aux = c_string(&mut malformed_heap, b"/N:");
        let malformed_baseline = malformed_heap.live_allocation_count();
        malformed_heap.trace_source_allocations();
        assert_eq!(
            OAD_Polymer_PrepareFrameShiftEdits(
                &mut malformed_heap,
                &malformed.data,
                malformed_inchi.as_const(),
                malformed_aux.as_const(),
                &malformed.edits,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(malformed_heap.source_allocation_calls(), 1);
        assert_eq!(malformed_heap.live_allocation_count(), malformed_baseline);

        for successful_allocations in 0..3 {
            let mut failure_heap = SourceHeap::default();
            let failed = fixture(&mut failure_heap, standard_unit(), true, 16);
            let inchi = c_string(&mut failure_heap, SHIFT_3_6);
            let aux = c_string(&mut failure_heap, IDENTITY_AUX);
            let baseline = failure_heap.live_allocation_count();
            failure_heap.fail_after_allocations(successful_allocations);
            assert_eq!(
                OAD_Polymer_PrepareFrameShiftEdits(
                    &mut failure_heap,
                    &failed.data,
                    inchi.as_const(),
                    aux.as_const(),
                    &failed.edits,
                ),
                Ok(_IS_ERROR as i32),
                "successful_allocations={successful_allocations}"
            );
            assert_eq!(failure_heap.source_allocation_calls(), successful_allocations + 1);
            assert_eq!(failure_heap.live_allocation_count(), baseline);
            assert!(edit_values(&failure_heap, &failed.edits).is_empty());
        }

        let mut skipped_heap = SourceHeap::default();
        let skipped = fixture(
            &mut skipped_heap,
            vec![(vec![3, 6], None, 2), (vec![3, 6], Some(vec![1, 2, 8, 7]), 1)],
            true,
            16,
        );
        assert_eq!(
            invoke(&mut skipped_heap, &skipped, SHIFT_3_6, IDENTITY_AUX),
            Ok(_IS_ERROR as i32)
        );
        assert!(edit_values(&skipped_heap, &skipped.edits).is_empty());

        let mut missing_heap = SourceHeap::default();
        let missing = fixture(
            &mut missing_heap,
            vec![(vec![2, 7], Some(vec![1, 2, 8, 7]), 2)],
            true,
            16,
        );
        assert_eq!(
            invoke(&mut missing_heap, &missing, SHIFT_3_6, IDENTITY_AUX),
            Ok(_IS_ERROR as i32)
        );

        let mut find_error_heap = SourceHeap::default();
        let find_error = fixture(
            &mut find_error_heap,
            vec![(vec![3, 6, 7], Some(vec![3, 6, 8, 7]), 2)],
            true,
            16,
        );
        assert_eq!(
            invoke(&mut find_error_heap, &find_error, SHIFT_3_6, IDENTITY_AUX,),
            Ok(_IS_OKAY as i32)
        );
        assert!(edit_values(&find_error_heap, &find_error.edits).is_empty());

        let mut nonstar_heap = SourceHeap::default();
        let nonstar = fixture(&mut nonstar_heap, standard_unit(), false, 16);
        assert_eq!(
            invoke(&mut nonstar_heap, &nonstar, SHIFT_3_6, IDENTITY_AUX),
            Ok(_IS_OKAY as i32)
        );
        assert!(edit_values(&nonstar_heap, &nonstar.edits).is_empty());

        let mut unchanged_heap = SourceHeap::default();
        let unchanged = fixture(&mut unchanged_heap, standard_unit(), true, 16);
        assert_eq!(
            invoke(
                &mut unchanged_heap,
                &unchanged,
                b"InChI=1B/C6H12.2Zz/zA(2,7,1)",
                IDENTITY_AUX,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert!(edit_values(&unchanged_heap, &unchanged.edits).is_empty());

        let mut same_end_heap = SourceHeap::default();
        let same_end = fixture(
            &mut same_end_heap,
            vec![(vec![2, 3, 6], Some(vec![1, 2, 8, 2]), 2)],
            true,
            16,
        );
        assert_eq!(
            invoke(&mut same_end_heap, &same_end, SHIFT_3_6, IDENTITY_AUX),
            Ok(_IS_OKAY as i32)
        );
        assert!(edit_values(&same_end_heap, &same_end.edits).is_empty());

        for (crossing, expected) in [
            (vec![1, 3, 8, 7], vec![7, 8, 6, 8, 3, 6, 3, 7]),
            (vec![1, 2, 8, 6], vec![2, 1, 3, 1, 3, 6, 2, 6]),
        ] {
            let mut one_side_heap = SourceHeap::default();
            let one_side = fixture(
                &mut one_side_heap,
                vec![(vec![2, 3, 4, 5, 6, 7], Some(crossing), 2)],
                true,
                16,
            );
            assert_eq!(
                invoke(&mut one_side_heap, &one_side, SHIFT_3_6, IDENTITY_AUX),
                Ok(_IS_OKAY as i32)
            );
            assert_eq!(edit_values(&one_side_heap, &one_side.edits), expected);
        }

        let mut selected_heap = SourceHeap::default();
        let mut selected = fixture(
            &mut selected_heap,
            vec![
                (vec![2, 3, 4, 5, 6, 7], Some(vec![1, 2, 8, 7]), 2),
                (vec![3, 4, 5, 6], Some(vec![1, 4, 8, 5]), 2),
            ],
            true,
            16,
        );
        selected.data.num_dimensions = 3;
        let selected_baseline = selected_heap.live_allocation_count();
        assert_eq!(
            invoke(&mut selected_heap, &selected, SHIFT_3_6, IDENTITY_AUX),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(selected_heap.source_allocation_calls(), 3);
        assert_eq!(selected_heap.live_allocation_count(), selected_baseline + 2);
        assert_eq!(
            edit_values(&selected_heap, &selected.edits),
            [4, 1, 3, 1, 5, 8, 6, 8, 3, 6, 4, 5]
        );
        assert_eq!(selected_heap.slice(selected.units[0].as_const()).unwrap()[0].cap1, 0);
        let selected_unit = &selected_heap.slice(selected.units[1].as_const()).unwrap()[0];
        assert_eq!(
            (
                selected_unit.cap1,
                selected_unit.end_atom1,
                selected_unit.end_atom2,
                selected_unit.cap2,
            ),
            (1, 4, 5, 8)
        );

        let mut append_failure_heap = SourceHeap::default();
        let append_failure = fixture(&mut append_failure_heap, standard_unit(), true, 4);
        let inchi = c_string(&mut append_failure_heap, SHIFT_3_6);
        let aux = c_string(&mut append_failure_heap, IDENTITY_AUX);
        let append_baseline = append_failure_heap.live_allocation_count();
        append_failure_heap.fail_after_allocations(3);
        assert_eq!(
            OAD_Polymer_PrepareFrameShiftEdits(
                &mut append_failure_heap,
                &append_failure.data,
                inchi.as_const(),
                aux.as_const(),
                &append_failure.edits,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(append_failure_heap.source_allocation_calls(), 4);
        assert_eq!(append_failure_heap.live_allocation_count(), append_baseline - 1);
        let failed_array = &append_failure_heap
            .slice(append_failure.edits.mod_bond.as_const())
            .unwrap()[0];
        assert!(failed_array.item.is_null());
        assert_eq!(
            (failed_array.allocated, failed_array.used, failed_array.increment),
            (4, 4, 4)
        );
    }

    #[test]
    fn source_port__runichi2__modscenter_init__line_2760() {
        fn sentinel() -> ModSCenterInfo {
            ModSCenterInfo {
                num: 91,
                valence: 92,
                n_stereo: 93,
                nbr: [94; 20],
                new_nbr: [95; 20],
            }
        }

        let mut selected = inp_ATOM {
            valence: 4,
            ..inp_ATOM::default()
        };
        selected.neighbor[..4].copy_from_slice(&[0, u16::MAX, 123, 7]);
        selected.bond_stereo[..4].copy_from_slice(&[1, -1, 6, -6]);
        let atoms = [inp_ATOM::default(), selected];
        let mut center = sentinel();
        assert_eq!(ModSCenter_Init(&mut center, &atoms, 1), Ok(()));
        assert_eq!(center.num, 1);
        assert_eq!(center.valence, 4);
        assert_eq!(center.n_stereo, 2);
        assert_eq!(&center.nbr[..4], &[0, i32::from(u16::MAX), 123, 7]);
        assert_eq!(&center.new_nbr[..4], &[0, i32::from(u16::MAX), 123, 7]);
        assert_eq!(&center.nbr[4..], &[94; 16]);
        assert_eq!(&center.new_nbr[4..], &[95; 16]);

        for valence in [-1_i8, 0] {
            let atom = inp_ATOM {
                valence,
                neighbor: [u16::MAX; 20],
                bond_stereo: [1; 20],
                ..inp_ATOM::default()
            };
            let mut empty = sentinel();
            assert_eq!(ModSCenter_Init(&mut empty, &[atom], 0), Ok(()));
            assert_eq!(empty.num, 0);
            assert_eq!(empty.valence, i32::from(valence));
            assert_eq!(empty.n_stereo, 0);
            assert_eq!(empty.nbr, [94; 20]);
            assert_eq!(empty.new_nbr, [95; 20]);
        }

        let mut full_atom = inp_ATOM {
            valence: 20,
            ..inp_ATOM::default()
        };
        for index in 0..20 {
            full_atom.neighbor[index] = (u16::MAX as usize - index) as u16;
            full_atom.bond_stereo[index] = if index % 2 == 0 { 1 } else { 6 };
        }
        let mut full = sentinel();
        assert_eq!(ModSCenter_Init(&mut full, &[full_atom.clone()], 0), Ok(()));
        assert_eq!(full.valence, 20);
        assert_eq!(full.n_stereo, 20);
        for index in 0..20 {
            assert_eq!(full.nbr[index], i32::from(full_atom.neighbor[index]));
            assert_eq!(full.new_nbr[index], full.nbr[index]);
        }

        for invalid_index in [-1, 1] {
            let mut invalid = sentinel();
            assert_eq!(
                ModSCenter_Init(&mut invalid, &[inp_ATOM::default()], invalid_index),
                Err(SourceHeapError::PointerOutOfBounds)
            );
            assert_eq!(invalid.num, invalid_index);
            assert_eq!(invalid.valence, 92);
            assert_eq!(invalid.n_stereo, 93);
            assert_eq!(invalid.nbr, [94; 20]);
            assert_eq!(invalid.new_nbr, [95; 20]);
        }

        let invalid_valence = inp_ATOM {
            valence: 21,
            bond_stereo: [1; 20],
            ..inp_ATOM::default()
        };
        let mut invalid = sentinel();
        assert_eq!(
            ModSCenter_Init(&mut invalid, &[invalid_valence], 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(invalid.num, 0);
        assert_eq!(invalid.valence, 21);
        assert_eq!(invalid.n_stereo, 93);
        assert_eq!(invalid.nbr, [94; 20]);
        assert_eq!(invalid.new_nbr, [95; 20]);
    }

    #[test]
    fn source_port__runichi2__ndefstereobonds__line_2775() {
        for stereo in i8::MIN..=i8::MAX {
            let mut atom = inp_ATOM {
                valence: 1,
                ..inp_ATOM::default()
            };
            atom.bond_stereo[0] = stereo;
            let atoms = [atom];
            assert_eq!(
                NDefStereoBonds(&atoms, 0, 1),
                Ok(i32::from(matches!(stereo, 1 | 6))),
                "pointed stereo={stereo}"
            );
            assert_eq!(
                NDefStereoBonds(&atoms, 0, 0),
                Ok(i32::from(matches!(stereo, -6 | -1 | 1 | 6))),
                "absolute stereo={stereo}"
            );
        }

        let mut all_slots = inp_ATOM {
            valence: 20,
            ..inp_ATOM::default()
        };
        all_slots.bond_stereo = [
            1, -1, 6, -6, 0, 2, -2, 3, -3, 4, -4, 5, -5, 7, -7, 127, -127, -128, 1, 6,
        ];
        let original = all_slots.clone();
        let all_atoms = [all_slots];
        assert_eq!(NDefStereoBonds(&all_atoms, 0, i32::MIN), Ok(4));
        assert_eq!(NDefStereoBonds(&all_atoms, 0, i32::MAX), Ok(4));
        assert_eq!(NDefStereoBonds(&all_atoms, 0, 0), Ok(6));
        assert_eq!(all_atoms[0], original);

        let mut zero_valence = inp_ATOM {
            valence: 0,
            ..inp_ATOM::default()
        };
        zero_valence.bond_stereo.fill(1);
        assert_eq!(NDefStereoBonds(&[zero_valence.clone()], 0, 0), Ok(0));
        zero_valence.valence = -1;
        assert_eq!(NDefStereoBonds(&[zero_valence], 0, 0), Ok(0));

        let invalid_valence = inp_ATOM {
            valence: 21,
            bond_stereo: [1; 20],
            ..inp_ATOM::default()
        };
        assert_eq!(
            NDefStereoBonds(&[invalid_valence], 0, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            NDefStereoBonds(&[inp_ATOM::default()], -1, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            NDefStereoBonds(&[inp_ATOM::default()], 1, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(NDefStereoBonds(&[], 0, 0), Err(SourceHeapError::PointerOutOfBounds));
    }

    #[test]
    fn source_port__runichi2__modscenter_addto__line_2803() {
        let mut empty = ModSCenterInfo {
            valence: 0,
            nbr: [91; 20],
            new_nbr: [92; 20],
            ..ModSCenterInfo::default()
        };
        assert_eq!(ModSCenter_AddTo(&mut empty, i32::MIN), Ok(()));
        assert_eq!(empty.valence, 1);
        assert_eq!(empty.new_nbr[0], i32::MIN);
        assert_eq!(&empty.new_nbr[1..], &[92; 19]);
        assert_eq!(empty.nbr, [91; 20]);

        let mut center = ModSCenterInfo {
            num: 17,
            valence: 2,
            n_stereo: 19,
            nbr: [71; 20],
            new_nbr: [73; 20],
        };
        center.new_nbr[..2].copy_from_slice(&[4, 5]);
        let duplicate = center.clone();
        assert_eq!(ModSCenter_AddTo(&mut center, 4), Ok(()));
        assert_eq!(center, duplicate);
        assert_eq!(ModSCenter_AddTo(&mut center, i32::MAX), Ok(()));
        assert_eq!(center.num, 17);
        assert_eq!(center.valence, 3);
        assert_eq!(center.n_stereo, 19);
        assert_eq!(&center.new_nbr[..3], &[4, 5, i32::MAX]);
        assert_eq!(&center.new_nbr[3..], &[73; 17]);
        assert_eq!(center.nbr, [71; 20]);

        let mut full = ModSCenterInfo {
            valence: 20,
            new_nbr: std::array::from_fn(|index| index as i32),
            ..ModSCenterInfo::default()
        };
        let full_before = full.clone();
        assert_eq!(ModSCenter_AddTo(&mut full, 19), Ok(()));
        assert_eq!(full, full_before);
        assert_eq!(
            ModSCenter_AddTo(&mut full, 20),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(full, full_before);

        for (invalid_valence, expected_error) in [
            (-1, SourceHeapError::SourceIntegerOverflow),
            (21, SourceHeapError::PointerOutOfBounds),
            (i32::MAX, SourceHeapError::PointerOutOfBounds),
        ] {
            let mut invalid = ModSCenterInfo {
                valence: invalid_valence,
                new_nbr: [3; 20],
                ..ModSCenterInfo::default()
            };
            let before = invalid.clone();
            assert_eq!(ModSCenter_AddTo(&mut invalid, 4), Err(expected_error));
            assert_eq!(invalid, before);
        }

        let mut searches_new_neighbors = ModSCenterInfo {
            valence: 1,
            nbr: [8; 20],
            new_nbr: [9; 20],
            ..ModSCenterInfo::default()
        };
        assert_eq!(ModSCenter_AddTo(&mut searches_new_neighbors, 8), Ok(()));
        assert_eq!(searches_new_neighbors.valence, 2);
        assert_eq!(&searches_new_neighbors.new_nbr[..2], &[9, 8]);
    }

    #[test]
    fn source_port__runichi2__modscenter_delfrom__line_2815() {
        fn center(neighbors: &[i32], new_neighbors: &[i32]) -> ModSCenterInfo {
            let mut result = ModSCenterInfo {
                num: 41,
                valence: neighbors.len() as i32,
                n_stereo: 43,
                nbr: [91; 20],
                new_nbr: [92; 20],
            };
            result.nbr[..neighbors.len()].copy_from_slice(neighbors);
            result.new_nbr[..new_neighbors.len()].copy_from_slice(new_neighbors);
            result
        }

        let mut middle = center(&[1, 2, 3, 4], &[10, 20, 30, 40]);
        let original_neighbors = middle.nbr;
        assert_eq!(ModSCenter_DelFrom(&mut middle, 2), Ok(()));
        assert_eq!(middle.num, 41);
        assert_eq!(middle.valence, 3);
        assert_eq!(middle.n_stereo, 43);
        assert_eq!(middle.nbr, original_neighbors);
        assert_eq!(&middle.new_nbr[..4], &[10, 30, 40, 40]);
        assert_eq!(&middle.new_nbr[4..], &[92; 16]);

        let mut first = center(&[1, 2, 3], &[10, 20, 30]);
        assert_eq!(ModSCenter_DelFrom(&mut first, 1), Ok(()));
        assert_eq!(first.valence, 2);
        assert_eq!(&first.new_nbr[..3], &[20, 30, 30]);

        let mut last = center(&[1, 2, 3], &[10, 20, 30]);
        assert_eq!(ModSCenter_DelFrom(&mut last, 3), Ok(()));
        assert_eq!(last.valence, 2);
        assert_eq!(&last.new_nbr[..3], &[10, 20, 30]);

        let mut duplicate = center(&[5, 5, 6], &[50, 51, 60]);
        assert_eq!(ModSCenter_DelFrom(&mut duplicate, 5), Ok(()));
        assert_eq!(duplicate.valence, 2);
        assert_eq!(&duplicate.new_nbr[..3], &[51, 60, 60]);

        let mut absent = center(&[1, 2], &[7, 8]);
        let absent_before = absent.clone();
        assert_eq!(ModSCenter_DelFrom(&mut absent, 7), Ok(()));
        assert_eq!(absent, absent_before);

        for valence in [-1, 0] {
            let mut empty = center(&[], &[]);
            empty.valence = valence;
            let before = empty.clone();
            assert_eq!(ModSCenter_DelFrom(&mut empty, 91), Ok(()));
            assert_eq!(empty, before);
        }

        let mut full_neighbors = [0_i32; 20];
        let mut full_new_neighbors = [0_i32; 20];
        for index in 0..20 {
            full_neighbors[index] = index as i32;
            full_new_neighbors[index] = 100 + index as i32;
        }
        let mut full = center(&full_neighbors, &full_new_neighbors);
        assert_eq!(ModSCenter_DelFrom(&mut full, 0), Ok(()));
        assert_eq!(full.valence, 19);
        for index in 0..19 {
            assert_eq!(full.new_nbr[index], 101 + index as i32);
        }
        assert_eq!(full.new_nbr[19], 119);

        let mut oversized_absent = center(&[1; 20], &[2; 20]);
        oversized_absent.valence = 21;
        let absent_before = oversized_absent.clone();
        assert_eq!(
            ModSCenter_DelFrom(&mut oversized_absent, 3),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(oversized_absent, absent_before);

        let mut partial = center(&[1; 20], &[100; 20]);
        partial.nbr[18] = 9;
        partial.new_nbr[19] = 119;
        partial.valence = 21;
        assert_eq!(
            ModSCenter_DelFrom(&mut partial, 9),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(partial.valence, 21);
        assert_eq!(partial.new_nbr[18], 119);
        assert_eq!(partial.new_nbr[19], 119);
    }

    #[test]
    fn source_port__runichi2__modscenter_ischanged__line_2836() {
        fn center(old: &[i32], new: &[i32]) -> ModSCenterInfo {
            let mut result = ModSCenterInfo {
                num: 0,
                valence: old.len() as i32,
                n_stereo: 1,
                nbr: [91; 20],
                new_nbr: [92; 20],
            };
            result.nbr[..old.len()].copy_from_slice(old);
            result.new_nbr[..new.len()].copy_from_slice(new);
            result
        }

        fn plain_atom(x: f64, y: f64, z: f64) -> inp_ATOM {
            inp_ATOM {
                x,
                y,
                z,
                ..inp_ATOM::default()
            }
        }

        fn stereo_atom() -> inp_ATOM {
            let mut atom = inp_ATOM {
                valence: 1,
                ..inp_ATOM::default()
            };
            atom.bond_stereo[0] = 1;
            atom
        }

        let mut early_heap = SourceHeap::default();
        let mut no_stereo = center(&[3, 1, 2], &[4, 1, 2]);
        no_stereo.n_stereo = 0;
        no_stereo.num = i32::MAX;
        let no_stereo_before = no_stereo.clone();
        assert_eq!(ModSCenter_IsChanged(&mut early_heap, &mut no_stereo, &[]), Ok(0));
        assert_eq!(no_stereo, no_stereo_before);
        assert_eq!(early_heap.live_allocation_count(), 0);

        let mut mismatch = center(&[3, 1, 2], &[4, 1, 2]);
        let mismatch_before = mismatch.clone();
        let mismatch_atoms = [inp_ATOM {
            valence: 2,
            ..inp_ATOM::default()
        }];
        assert_eq!(
            ModSCenter_IsChanged(&mut early_heap, &mut mismatch, &mismatch_atoms),
            Ok(-1)
        );
        assert_eq!(mismatch, mismatch_before);
        assert_eq!(early_heap.live_allocation_count(), 0);

        let mut atoms = vec![plain_atom(0.0, 0.0, 0.0); 7];
        atoms[0].valence = 3;
        atoms[1] = plain_atom(1.0, 0.0, 0.0);
        atoms[2] = plain_atom(0.0, 1.0, 0.0);
        atoms[3] = plain_atom(0.0, 0.0, 1.0);
        atoms[4] = plain_atom(-1.0, 0.0, f64::NAN);

        let mut full = center(&[3, 1, 2], &[4, 2, 1]);
        assert_eq!(ModSCenter_IsChanged(&mut early_heap, &mut full, &atoms), Ok(-1));
        assert_eq!(&full.nbr[..3], &[1, 2, 3]);
        assert_eq!(&full.new_nbr[..3], &[1, 2, 4]);
        assert_eq!(&full.nbr[3..], &[91; 17]);
        assert_eq!(&full.new_nbr[3..], &[92; 17]);
        assert_eq!(early_heap.live_allocation_count(), 0);

        let mut no_common = center(&[3, 1, 2], &[6, 4, 5]);
        assert_eq!(ModSCenter_IsChanged(&mut early_heap, &mut no_common, &atoms), Ok(-1));
        assert_eq!(&no_common.nbr[..3], &[1, 2, 3]);
        assert_eq!(&no_common.new_nbr[..3], &[4, 5, 6]);

        let mut identical = center(&[3, 1, 2], &[3, 2, 1]);
        assert_eq!(ModSCenter_IsChanged(&mut early_heap, &mut identical, &atoms), Ok(-1));
        assert_eq!(&identical.nbr[..3], &[1, 2, 3]);
        assert_eq!(&identical.new_nbr[..3], &[1, 2, 3]);

        let mut multiple = center(&[3, 1, 2], &[5, 1, 4]);
        assert_eq!(ModSCenter_IsChanged(&mut early_heap, &mut multiple, &atoms), Ok(-1));
        assert_eq!(&multiple.nbr[..3], &[1, 2, 3]);
        assert_eq!(&multiple.new_nbr[..3], &[1, 4, 5]);

        let mut stereo_base_atoms = atoms.clone();
        stereo_base_atoms[1] = stereo_atom();
        stereo_base_atoms[2] = stereo_atom();
        stereo_base_atoms[3] = stereo_atom();
        let mut no_nonstereo_base = center(&[3, 1, 2], &[3, 2, 1]);
        assert_eq!(
            ModSCenter_IsChanged(&mut early_heap, &mut no_nonstereo_base, &stereo_base_atoms,),
            Ok(-1)
        );
        assert_eq!(&no_nonstereo_base.nbr[..3], &[1, 2, 3]);

        let mut removed_stereo_atoms = atoms.clone();
        removed_stereo_atoms[3] = stereo_atom();
        let mut removed_stereo = center(&[3, 1, 2], &[4, 2, 1]);
        assert_eq!(
            ModSCenter_IsChanged(&mut early_heap, &mut removed_stereo, &removed_stereo_atoms,),
            Ok(-1)
        );
        assert_eq!(&removed_stereo.nbr[..3], &[1, 2, 3]);
        assert_eq!(&removed_stereo.new_nbr[..3], &[1, 2, 4]);

        let mut invalid_atoms = atoms.clone();
        invalid_atoms[0].valence = 1;
        let mut invalid_index = center(&[99], &[99]);
        assert_eq!(
            ModSCenter_IsChanged(&mut early_heap, &mut invalid_index, &invalid_atoms),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(invalid_index.nbr[0], 99);
        assert_eq!(invalid_index.new_nbr[0], 99);
        assert_eq!(early_heap.live_allocation_count(), 0);

        let mut allocation_injection_heap = SourceHeap::default();
        let mut allocation_injection = center(&[3, 1, 2], &[4, 2, 1]);
        allocation_injection_heap.fail_after_allocations(0);
        assert_eq!(
            ModSCenter_IsChanged(&mut allocation_injection_heap, &mut allocation_injection, &atoms,),
            Ok(-1)
        );
        assert_eq!(&allocation_injection.nbr[..3], &[1, 2, 3]);
        assert_eq!(&allocation_injection.new_nbr[..3], &[1, 2, 4]);
        assert_eq!(allocation_injection_heap.source_allocation_calls(), 0);
        assert_eq!(allocation_injection_heap.live_allocation_count(), 0);

        let mut invalid_center = center(&[1], &[1]);
        invalid_center.num = -1;
        assert_eq!(
            ModSCenter_IsChanged(&mut SourceHeap::default(), &mut invalid_center, &atoms),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

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
        let num_inp_bonds = atoms.iter().map(|atom| i32::from(atom.valence)).sum::<i32>() / 2;
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

    fn component_original(heap: &mut SourceHeap, lengths: Vec<u16>) -> ORIG_ATOM_DATA {
        let mut atoms = vec![inp_ATOM::default(); 3];
        atoms[0].component = 1;
        atoms[0].orig_at_number = 1;
        atoms[1].component = 2;
        atoms[1].orig_at_number = 2;
        atoms[1].valence = 1;
        atoms[1].neighbor[0] = 2;
        atoms[2].component = 2;
        atoms[2].orig_at_number = 3;
        atoms[2].valence = 1;
        atoms[2].neighbor[0] = 1;
        ORIG_ATOM_DATA {
            at: heap.allocate_model_storage(atoms).unwrap(),
            num_inp_atoms: 3,
            num_inp_bonds: 1,
            num_components: 2,
            nCurAtLen: heap.allocate_model_storage(lengths).unwrap(),
            ..ORIG_ATOM_DATA::default()
        }
    }

    #[test]
    fn source_port__runichi2__getonecomponent__line_345() {
        let mut heap = SourceHeap::default();
        let original = component_original(&mut heap, vec![1, 2, 0]);
        let mut clock = INCHI_CLOCK::default();
        let mut structure = STRUCT_DATA {
            ulStructTime: 5,
            ..STRUCT_DATA::default()
        };
        let mut current = INP_ATOM_DATA::default();
        assert_eq!(
            GetOneComponent(
                &mut heap,
                &mut clock,
                &mut structure,
                &INPUT_PARMS::default(),
                None,
                None,
                &mut current,
                &original,
                1,
                12,
                100,
                2_100,
            ),
            Ok(0)
        );
        assert_eq!(current.num_at, 2);
        let extracted = heap.slice(current.at.as_const()).unwrap();
        assert_eq!(extracted[0].orig_at_number, 2);
        assert_eq!(extracted[0].orig_compt_at_numb, 1);
        assert_eq!(extracted[0].neighbor[0], 1);
        assert_eq!(extracted[1].orig_at_number, 3);
        assert_eq!(extracted[1].orig_compt_at_numb, 2);
        assert_eq!(extracted[1].neighbor[0], 0);
        assert_eq!(structure.ulStructTime, 7);
        assert_eq!(structure.nErrorType, 0);

        let mut mismatch_heap = SourceHeap::default();
        let mismatch_original = component_original(&mut mismatch_heap, vec![2, 2, 0]);
        let label = c_string(&mut mismatch_heap, b"ID");
        let value = c_string(&mut mismatch_heap, b"42");
        let parameters = INPUT_PARMS {
            pSdfLabel: label,
            pSdfValue: value,
            ..INPUT_PARMS::default()
        };
        let mut mismatch_log = log_stream();
        let mut mismatch_structure = STRUCT_DATA::default();
        let mut mismatch_current = INP_ATOM_DATA::default();
        assert_eq!(
            GetOneComponent(
                &mut mismatch_heap,
                &mut INCHI_CLOCK::default(),
                &mut mismatch_structure,
                &parameters,
                Some(&mut mismatch_log),
                None,
                &mut mismatch_current,
                &mismatch_original,
                0,
                7,
                0,
                0,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(mismatch_current.num_at, 1);
        assert_eq!(mismatch_structure.nErrorCode, CT_ATOMCOUNT_ERR);
        assert_eq!(mismatch_structure.nErrorType, _IS_ERROR as i32);
        assert_eq!(
            log_bytes(&mismatch_heap, &mismatch_log),
            b"Cannot extract Component #1 structure #7. ID=42\n"
        );

        let mut empty_structure = STRUCT_DATA::default();
        let mut empty_current = INP_ATOM_DATA::default();
        let mut empty_log = log_stream();
        assert_eq!(
            GetOneComponent(
                &mut heap,
                &mut INCHI_CLOCK::default(),
                &mut empty_structure,
                &INPUT_PARMS::default(),
                Some(&mut empty_log),
                None,
                &mut empty_current,
                &original,
                2,
                8,
                0,
                0,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(empty_current.num_at, 0);
        assert_eq!(empty_structure.nErrorCode, CT_UNKNOWN_ERR);
        assert_eq!(
            log_bytes(&heap, &empty_log),
            b"Cannot extract Component #3 structure #8.\n"
        );

        let mut failure_heap = SourceHeap::default();
        let failure_original = component_original(&mut failure_heap, vec![1, 2]);
        failure_heap.fail_after_allocations(1);
        let mut failure_structure = STRUCT_DATA::default();
        let mut failure_current = INP_ATOM_DATA::default();
        assert_eq!(
            GetOneComponent(
                &mut failure_heap,
                &mut INCHI_CLOCK::default(),
                &mut failure_structure,
                &INPUT_PARMS::default(),
                None,
                None,
                &mut failure_current,
                &failure_original,
                0,
                9,
                0,
                0,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(failure_current.num_at, CT_OUT_OF_RAM);
        assert_eq!(failure_structure.nErrorCode, CT_OUT_OF_RAM);
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
            let clock = heap.allocate_model_storage(vec![INCHI_CLOCK::default()]).unwrap();
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
        let carbon = atoms.iter().position(|atom| atom.orig_at_number == 1).unwrap();
        let hydroxyl = atoms.iter().position(|atom| atom.orig_at_number == 2).unwrap();
        let opened = atoms.iter().position(|atom| atom.orig_at_number == 3).unwrap();
        let carbon_to_hydroxyl = (0..atoms[carbon].valence as usize)
            .find(|order| atoms[carbon].neighbor[*order] as usize == hydroxyl)
            .unwrap();
        assert_eq!(atoms[carbon].bond_type[carbon_to_hydroxyl], BOND_DOUBLE as u8);
        assert!(!(0..atoms[carbon].valence as usize).any(|order| atoms[carbon].neighbor[order] as usize == opened));
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
