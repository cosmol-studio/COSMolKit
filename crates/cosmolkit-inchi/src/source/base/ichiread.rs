use crate::source::base::ichi_io::{inchi_ios_eprint, inchi_strbuf_close, inchi_strbuf_init};
use crate::source::base::ichicano::{InchiTimeElapsed, InchiTimeGet};
use crate::source::base::ichimak2::AllocateAndFillHillFormula;
use crate::source::base::ichimake::CompareReversedINChI;
use crate::source::base::ichirvr7::FreeInpInChI;
use crate::source::base::runichi4::{FreeAllINChIArrays, SortAndPrintINChI};
use crate::source::base::util::{inchi_calloc, inchi_free, is_el_a_metal};
use crate::source_types::{
    _IS_ERROR, _IS_FATAL, AB_MAX_WELL_DEFINED_PARITY, AB_MIN_WELL_DEFINED_PARITY, AT_NUMB,
    CANON_GLOBALS, FILE, FLAG_SORT_PRINT_ReChI_PREFIX, INCHI_BAS, INCHI_CLOCK, INCHI_IOS_STRING,
    INCHI_IOSTREAM, INCHI_MODE, INCHI_NUM, INCHI_OUT_INCHI_GEN_ERROR, INCHI_OUT_SAVEOPT,
    INCHI_OUT_SDFILE_ONLY, INCHI_OUT_STDINCHI, INCHI_REC, INCHI_STRBUF_INITIAL_SIZE,
    INCHI_STRBUF_SIZE_INCREMENT, INChI, INChI_Aux, INPUT_PARMS, InpInChI, NORM_CANON_FLAGS,
    PINChI_Aux2, PINChI2, READ_INCHI_SPLIT_OUTPUT, READ_INCHI_TO_STRUCTURE, REQ_MODE_BASIC,
    REQ_MODE_ISO_STEREO, REQ_MODE_RACEMIC_STEREO, REQ_MODE_RELATIVE_STEREO, REQ_MODE_SB_IGN_ALL_UU,
    REQ_MODE_SC_IGN_ALL_UU, REQ_MODE_STEREO, RI_ERR_ALLOC, RI_ERR_EOF, RI_ERR_EOL, RI_ERR_PROGR,
    RI_ERR_SYNTAX, SAVE_OPT_15T, SAVE_OPT_FIXEDH, SAVE_OPT_KET, SAVE_OPT_PT_06_00,
    SAVE_OPT_PT_13_00, SAVE_OPT_PT_16_00, SAVE_OPT_PT_18_00, SAVE_OPT_PT_22_00, SAVE_OPT_PT_39_00,
    SAVE_OPT_RECMET, SAVE_OPT_SLUUD, SAVE_OPT_SUU, STRUCT_DATA, SourceConstPointer,
    SourceFormatArgument, SourceHeap, SourceHeapError, SourceMutPointer, SourceVaList,
    StrFromINChI, T_GROUP_HDR_LEN, TAUT_NON, TAUT_NUM, TAUT_YES, TG_FLAG_ARSINE_STEREO,
    TG_FLAG_DISCONNECT_COORD, TG_FLAG_DISCONNECT_COORD_DONE, TG_FLAG_PHOSPHINE_STEREO,
    TG_FLAG_RECONNECT_COORD, clock_t, inchiTime, local_ichiread::MODE_PIXH,
    local_ichiread::tagModeProtonIsoExchgH_MODE_PIXH_ADD_A_PIXH_COMPONENT as MODE_PIXH_ADD_A_PIXH_COMPONENT,
    local_ichiread::tagModeProtonIsoExchgH_MODE_PIXH_ADD_TO_EACH as MODE_PIXH_ADD_TO_EACH,
    local_ichiread::tagModeProtonIsoExchgH_MODE_PIXH_ADD_TO_FIRST as MODE_PIXH_ADD_TO_FIRST,
    tagInputType_INPUT_INCHI,
};

#[allow(non_snake_case)]
pub(crate) fn SetHillFormFromInChI(
    heap: &mut SourceHeap,
    one_input: &mut InpInChI,
) -> Result<i32, SourceHeapError> {
    // Active libinchi configuration: FIX_DALKE_BUGS == 1 includes this function.
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:583 SetHillFormFromInChI
    // INCHI✔️❌: int SetHillFormFromInChI(InpInChI* OneInput)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int iINChI, iTaut, iComp, num_diff;
    // INCHI✔️❌:     INChI* pINChI;
    // INCHI✔️❌:     char* szHillFormulaOld;
    // INCHI✔️❌:     for (iINChI = 0, num_diff = 0; iINChI < INCHI_NUM; iINChI++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (iTaut = TAUT_NON; iTaut < TAUT_NUM; iTaut++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (iComp = 0; iComp < OneInput->nNumComponents[iINChI][iTaut]; iComp++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pINChI = &OneInput->pInpInChI[iINChI][iTaut][iComp];
    // INCHI✔️❌:                 if (!pINChI->nNumberOfAtoms || pINChI->bDeleted || !pINChI->szHillFormula || !pINChI->szHillFormula[0])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 szHillFormulaOld = pINChI->szHillFormula;
    // INCHI✔️❌:                 pINChI->szHillFormula = AllocateAndFillHillFormula(pINChI);
    // INCHI✔️❌:                 num_diff += !pINChI->szHillFormula || !pINChI->szHillFormula[0] || strcmp(pINChI->szHillFormula, szHillFormulaOld);
    // INCHI✔️❌:                 inchi_free(szHillFormulaOld);
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return num_diff;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: SetHillFormFromInChI

    let mut num_differences = 0_i32;
    for representation in 0..INCHI_NUM as usize {
        for tautomer in TAUT_NON as usize..TAUT_NUM as usize {
            let component_count = one_input.nNumComponents[representation][tautomer];
            let component_count = if component_count > 0 {
                usize::try_from(component_count)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
            } else {
                0
            };
            let components = one_input.pInpInChI[representation][tautomer];
            for component in 0..component_count {
                let snapshot = heap
                    .slice(components.as_const())?
                    .get(component)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if snapshot.nNumberOfAtoms == 0
                    || snapshot.bDeleted != 0
                    || snapshot.szHillFormula.is_null()
                {
                    continue;
                }
                let old_formula = snapshot.szHillFormula;
                let old_bytes = heap.slice(old_formula.as_const())?;
                if old_bytes
                    .first()
                    .copied()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    == 0
                {
                    continue;
                }

                let new_formula = AllocateAndFillHillFormula(heap, &snapshot)?;
                heap.slice_mut(components)?
                    .get_mut(component)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .szHillFormula = new_formula;

                let differs = if new_formula.is_null() {
                    true
                } else {
                    let new_bytes = heap.slice(new_formula.as_const())?;
                    if new_bytes
                        .first()
                        .copied()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        == 0
                    {
                        true
                    } else {
                        let new_end = new_bytes
                            .iter()
                            .position(|byte| *byte == 0)
                            .ok_or(SourceHeapError::MissingNulTerminator)?;
                        let old_bytes = heap.slice(old_formula.as_const())?;
                        let old_end = old_bytes
                            .iter()
                            .position(|byte| *byte == 0)
                            .ok_or(SourceHeapError::MissingNulTerminator)?;
                        new_bytes[..new_end] != old_bytes[..old_end]
                    }
                };
                num_differences = num_differences.wrapping_add(i32::from(differs));
                inchi_free(heap, old_formula)?;
            }
        }
    }
    Ok(num_differences)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn ConvertInChI2InChI(
    heap: &mut SourceHeap,
    input_parameters: &mut INPUT_PARMS,
    one_input: &mut InpInChI,
    out: &mut INCHI_IOSTREAM,
    log: &mut INCHI_IOSTREAM,
    structure_data: &STRUCT_DATA,
    num_components: &mut [i32; INCHI_NUM as usize],
    proton_modes: &[MODE_PIXH; INCHI_NUM as usize],
    current_header: &mut SourceMutPointer<i8>,
    input_number: i64,
    _num_errors: &mut i64,
    save_option_bits: u8,
    start_time: &mut inchiTime,
    processing_time: &mut i64,
    clock: &mut INCHI_CLOCK,
    canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    stdout: SourceMutPointer<FILE>,
    start_clock_result: clock_t,
    end_clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // Active libinchi configuration: TARGET_API_LIB excludes the standalone
    // InChIKey block and the non-library error-reporting block.
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:11003 ConvertInChI2InChI
    // INCHI✔️❌: int ConvertInChI2InChI(INPUT_PARMS* ip,
    // INCHI✔️❌:     InpInChI* pOneInput,
    // INCHI✔️❌:     INCHI_IOSTREAM* pOut,
    // INCHI✔️❌:     INCHI_IOSTREAM* pLog,
    // INCHI✔️❌:     STRUCT_DATA* sd,
    // INCHI✔️❌:     int            num_components[INCHI_NUM],
    // INCHI✔️❌:     MODE_PIXH      nModeProtonIsoExchgH[INCHI_NUM],
    // INCHI✔️❌:     char** pszCurHdr,
    // INCHI✔️❌:     long           num_inp,
    // INCHI✔️❌:     long* num_errors,
    // INCHI✔️❌:     unsigned char  save_opt_bits,
    // INCHI✔️❌:     inchiTime* pulTStart,
    // INCHI✔️❌:     long* ulProcessingTime,
    // INCHI✔️❌:     struct         tagINCHI_CLOCK* ic,
    // INCHI✔️❌:     struct         tagCANON_GLOBALS* pCG)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret, tmp;
    // INCHI✔️❌:
    // INCHI✔️❌:     InchiTimeGet(pulTStart);
    // INCHI✔️❌:
    // INCHI✔️❌:     tmp = ip->bNoStructLabels;
    // INCHI✔️❌:     ip->bNoStructLabels = 1;
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:         ip->pSdfValue = NULL;
    // INCHI✔️❌:     ip->pSdfLabel = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( FIX_DALKE_BUGS == 1 )
    // INCHI✔️❌:     SetHillFormFromInChI(pOneInput);
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = OutputInChIAsRequested(pCG, pOut, pLog, ip, sd,
    // INCHI✔️❌:         pOneInput, num_components,
    // INCHI✔️❌:         nModeProtonIsoExchgH,
    // INCHI✔️❌:         num_inp, save_opt_bits);
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( !defined(TARGET_API_LIB) && defined(TARGET_EXE_STANDALONE) )
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Calculate InChIKey if requested */
    // INCHI✔️❌:     /* However, do not calculate/write it if this function is called from within dll */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         char ik_string[256];    /* Resulting InChIKey string */
    // INCHI✔️❌:         int ik_ret = 0;           /* InChIKey-calc result code */
    // INCHI✔️❌:         int xhash1, xhash2;
    // INCHI✔️❌:         char szXtra1[65], szXtra2[65];
    // INCHI✔️❌:
    // INCHI✔️❌:         inchi_ios_flush2(pLog, stderr);
    // INCHI✔️❌:
    // INCHI✔️❌:         /* post-1.02b addition - correctly treat tabbed output with InChIKey */
    // INCHI✔️❌:         if (ip->bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (ip->bCalcInChIHash != INCHIHASH_NONE)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (pOut->s.pStr)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (pOut->s.nUsedLength > 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (pOut->s.pStr[pOut->s.nUsedLength - 1] == '\n')
    // INCHI✔️❌:                         {    /* replace LF with TAB */
    // INCHI✔️❌:                             pOut->s.pStr[pOut->s.nUsedLength - 1] = '\t';
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (ip->bCalcInChIHash == INCHIHASH_NONE)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* inchi_ios_flush(pOut); */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             char* buf = NULL;
    // INCHI✔️❌:             size_t slen = pOut->s.nUsedLength;
    // INCHI✔️❌:             extract_inchi_substring(&buf, pOut->s.pStr, slen);
    // INCHI✔️❌:
    // INCHI✔️❌:             if (NULL != buf)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 xhash1 = xhash2 = 0;
    // INCHI✔️❌:                 if ((ip->bCalcInChIHash == INCHIHASH_KEY_XTRA1) ||
    // INCHI✔️❌:                     (ip->bCalcInChIHash == INCHIHASH_KEY_XTRA1_XTRA2))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     xhash1 = 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if ((ip->bCalcInChIHash == INCHIHASH_KEY_XTRA2) ||
    // INCHI✔️❌:                     (ip->bCalcInChIHash == INCHIHASH_KEY_XTRA1_XTRA2))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     xhash2 = 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 ik_ret = GetINCHIKeyFromINCHI(buf,
    // INCHI✔️❌:                     xhash1,
    // INCHI✔️❌:                     xhash2,
    // INCHI✔️❌:                     ik_string,
    // INCHI✔️❌:                     szXtra1,
    // INCHI✔️❌:                     szXtra2);
    // INCHI✔️❌:                 inchi_free(buf);
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ik_ret = INCHIKEY_NOT_ENOUGH_MEMORY;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:             if (ik_ret == INCHIKEY_OK)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_print(pOut, "InChIKey=%-s\n", ik_string);
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_print(pLog, "Warning (Could not compute InChIKey: ", num_inp);
    // INCHI✔️❌:                 switch (ik_ret)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                 case INCHIKEY_UNKNOWN_ERROR:
    // INCHI✔️❌:                     inchi_ios_print(pLog, "unresolved error)");
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case INCHIKEY_EMPTY_INPUT:
    // INCHI✔️❌:                     inchi_ios_print(pLog, "got an empty string)");
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case INCHIKEY_INVALID_INCHI_PREFIX:
    // INCHI✔️❌:                 case INCHIKEY_INVALID_INCHI:
    // INCHI✔️❌:                 case INCHIKEY_INVALID_STD_INCHI:
    // INCHI✔️❌:                     inchi_ios_print(pLog, "no valid InChI string found)");
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case INCHIKEY_NOT_ENOUGH_MEMORY:
    // INCHI✔️❌:                     inchi_ios_print(pLog, "not enough memory to treat the string)");
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 default:inchi_ios_print(pLog, "internal program error)");
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 inchi_ios_print(pLog, " structure #%-lu.\n", num_inp);
    // INCHI✔️❌:                 if (ip->bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     inchi_ios_print(pOut, "\n");
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             } /* if (ip->bCalcInChIHash!=INCHIHASH_NONE) */
    // INCHI✔️❌:
    // INCHI✔️❌:             inchi_ios_flush(pOut);
    // INCHI✔️❌:             inchi_ios_flush2(pLog, stderr);
    // INCHI✔️❌:         }
    // INCHI✔️❌:     } /* Calculate InChIKey if requested */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     ip->bNoStructLabels = tmp;
    // INCHI✔️❌:
    // INCHI✔️❌: #ifndef TARGET_API_LIB
    // INCHI✔️❌:     if (ret < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         if (*pszCurHdr && (*pszCurHdr)[0])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(pLog, "Error %d creating InChI string %s\n", ret, *pszCurHdr);
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_eprint(pLog, "Error %d creating InChI string, Structure %ld\n", ret, num_inp);
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (ip->bINChIOutputOptions2 & INCHI_OUT_INCHI_GEN_ERROR)
    // INCHI✔️❌:         {/* inchi_ios_eprint( pOut, "InChICreationError!\n"); *//* emit err string */
    // INCHI✔️❌:             if (ip->bINChIOutputOptions & INCHI_OUT_STDINCHI)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_eprint(pOut, "InChI=1S//\n");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_eprint(pOut, "InChI=1//\n");
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         (*num_errors)++;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( !defined(TARGET_API_LIB) && !defined(TARGET_EXE_STANDALONE) )
    // INCHI✔️❌:     else
    // INCHI✔️❌:         if (*pszCurHdr && (*pszCurHdr)[0])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_fprintf(stderr, "%s\r", *pszCurHdr);
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (*pszCurHdr)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free(*pszCurHdr);
    // INCHI✔️❌:         *pszCurHdr = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:         * ulProcessingTime += InchiTimeElapsed(ic, pulTStart);
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: ConvertInChI2InChI

    InchiTimeGet(start_time, start_clock_result);
    let old_no_structure_labels = input_parameters.bNoStructLabels;
    input_parameters.bNoStructLabels = 1;
    input_parameters.pSdfValue = SourceMutPointer::null();
    input_parameters.pSdfLabel = SourceMutPointer::null();

    let result = (|| {
        SetHillFormFromInChI(heap, one_input)?;
        OutputInChIAsRequested(
            heap,
            canonical_globals,
            out,
            log,
            input_parameters,
            structure_data,
            one_input,
            num_components,
            proton_modes,
            input_number,
            save_option_bits,
            stdout,
        )
    })();

    input_parameters.bNoStructLabels = old_no_structure_labels;
    if !current_header.is_null() {
        inchi_free(heap, *current_header)?;
        *current_header = SourceMutPointer::null();
    }
    *processing_time =
        processing_time.wrapping_add(InchiTimeElapsed(clock, Some(start_time), end_clock_result));
    result
}

#[allow(non_snake_case)]
pub(crate) fn GetNumNeighborsFromInchi(
    heap: &SourceHeap,
    p_inchi: Option<&INChI>,
    mut n_at_number: AT_NUMB,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1646 GetNumNeighborsFromInchi
    // INCHI✔️❌: int GetNumNeighborsFromInchi(INChI* pInChI, AT_NUMB nAtNumber)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, n_vertex, n_neigh, nNumNeigh, bTautAtom, nNumH, nTotNumNeigh, num_atoms;
    // INCHI✔️❌:     AT_NUMB  taut_at_number;
    // INCHI✔️❌:     nAtNumber -= 1;
    // INCHI✔️❌:     nNumNeigh = 0; /* number of bonds */
    // INCHI✔️❌:     bTautAtom = 0; /* 1 if atom belongs to a Mobile-H group */
    // INCHI✔️❌:     nNumH = 0; /* number of terminal neighbors H */
    // INCHI✔️❌:     num_atoms = 0; /* djb-rwth: initialisation with pInChI below */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (pInChI) /* djb-rwth: fixing a NULL pointer dereference */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         num_atoms = pInChI->nNumberOfAtoms;
    // INCHI✔️❌:         /* from RestoreAtomConnectionsSetStereo() */
    // INCHI✔️❌:         /* Connection table structure:
    // INCHI✔️❌:         Vert(1) [, Neigh(11), Neigh(12),...], Vert(2) [, Neigh(2,1), Neigh(2,2),...] ...
    // INCHI✔️❌:         where Neigh(i,1) < Neigh(i,2) <... < Vert(i);
    // INCHI✔️❌:         Vert(i) < Vert(i+1)
    // INCHI✔️❌:         */
    // INCHI✔️❌:         for (i = 1, n_vertex = pInChI->nConnTable[0] - 1; i < pInChI->lenConnTable; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((n_neigh = pInChI->nConnTable[i] - 1) < n_vertex)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /*  vertex - neighbor connection */
    // INCHI✔️❌:                 nNumNeigh += (nAtNumber == n_vertex || nAtNumber == n_neigh);
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {/* n_neigh is the next vertex */
    // INCHI✔️❌:                 if ((n_vertex = n_neigh) >= num_atoms)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return  RI_ERR_PROGR;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* is atom tautomeric, from GetTgroupInfoFromInChI() */
    // INCHI✔️❌:     if (pInChI && pInChI->lenTautomer > 1 && pInChI->nTautomer && pInChI->nTautomer[0] > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int itg, len_tg;
    // INCHI✔️❌:         int tot_len_tg = pInChI->lenTautomer - T_GROUP_HDR_LEN * pInChI->nTautomer[0] - 1; /* number of endpoints */
    // INCHI✔️❌:         j = 1; /* index in pInChI->nTautomer[] */
    // INCHI✔️❌:         i = 0; /* index in ti->nEndpointAtomNumber[] */
    // INCHI✔️❌:         for (itg = 0; itg < pInChI->nTautomer[0]; itg++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             len_tg = pInChI->nTautomer[j]; /* t-group length not including pInChI->nTautomer[j] */
    // INCHI✔️❌:             j += T_GROUP_HDR_LEN;   /* skip t-group header */
    // INCHI✔️❌:             len_tg -= T_GROUP_HDR_LEN - 1;
    // INCHI✔️❌:             for (; 0 < len_tg--; j++, i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 taut_at_number = pInChI->nTautomer[j] - 1; /* Mobile-H group atom number */
    // INCHI✔️❌:                 bTautAtom += (taut_at_number == nAtNumber);
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i != tot_len_tg)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return RI_ERR_PROGR;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* count hydrogen neighbors */
    // INCHI✔️❌:     if (pInChI && pInChI->nNum_H) /* djb-rwth: condition added for fixing a NULL pointer dereference */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nNumH = pInChI->nNum_H[nAtNumber];
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* conclusion: if not tautomeric then return positive number, otherwise add 1000 */
    // INCHI✔️❌:     nTotNumNeigh = nNumNeigh + nNumH;
    // INCHI✔️❌:     if (bTautAtom)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nTotNumNeigh += 1000;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return nTotNumNeigh;
    // INCHI✔️❌:
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: GetNumNeighborsFromInchi

    n_at_number = n_at_number.wrapping_sub(1);
    let mut number_of_neighbors = 0_i32;
    let mut tautomeric_atom = 0_i32;
    let mut number_of_hydrogens = 0_i32;

    if let Some(inchi) = p_inchi {
        let connection_table = heap.slice(inchi.nConnTable.as_const())?;
        let mut vertex = i32::from(
            connection_table
                .first()
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        ) - 1;
        let connection_length = if inchi.lenConnTable > 1 {
            usize::try_from(inchi.lenConnTable)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
        } else {
            1
        };
        for index in 1..connection_length {
            let neighbor = i32::from(
                connection_table
                    .get(index)
                    .copied()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            ) - 1;
            if neighbor < vertex {
                number_of_neighbors = number_of_neighbors.wrapping_add(i32::from(
                    i32::from(n_at_number) == vertex || i32::from(n_at_number) == neighbor,
                ));
            } else {
                vertex = neighbor;
                if vertex >= inchi.nNumberOfAtoms {
                    return Ok(RI_ERR_PROGR);
                }
            }
        }
    }

    if let Some(inchi) = p_inchi.filter(|value| value.lenTautomer > 1 && !value.nTautomer.is_null())
    {
        let tautomer = heap.slice(inchi.nTautomer.as_const())?;
        let group_count = tautomer
            .first()
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if group_count > 0 {
            let total_endpoint_count = inchi
                .lenTautomer
                .wrapping_sub((T_GROUP_HDR_LEN as i32).wrapping_mul(i32::from(group_count)))
                .wrapping_sub(1);
            let mut index = 1_usize;
            let mut endpoint_count = 0_i32;
            for _ in 0..group_count {
                let mut group_length = i32::from(
                    tautomer
                        .get(index)
                        .copied()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                index = index.wrapping_add(T_GROUP_HDR_LEN as usize);
                group_length = group_length.wrapping_sub(T_GROUP_HDR_LEN as i32 - 1);
                while group_length > 0 {
                    group_length = group_length.wrapping_sub(1);
                    let tautomer_atom = tautomer
                        .get(index)
                        .copied()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .wrapping_sub(1);
                    tautomeric_atom =
                        tautomeric_atom.wrapping_add(i32::from(tautomer_atom == n_at_number));
                    index = index.wrapping_add(1);
                    endpoint_count = endpoint_count.wrapping_add(1);
                }
            }
            if endpoint_count != total_endpoint_count {
                return Ok(RI_ERR_PROGR);
            }
        }
    }

    if let Some(inchi) = p_inchi.filter(|value| !value.nNum_H.is_null()) {
        number_of_hydrogens = i32::from(
            *heap
                .slice(inchi.nNum_H.as_const())?
                .get(usize::from(n_at_number))
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
    }
    let mut total = number_of_neighbors.wrapping_add(number_of_hydrogens);
    if tautomeric_atom != 0 {
        total = total.wrapping_add(1000);
    }
    Ok(total)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn CountStereoTypes(
    heap: &SourceHeap,
    p_inchi: &INChI,
    num_known_sb: &mut i32,
    num_known_sc: &mut i32,
    num_unk_und_sb: &mut i32,
    num_unk_und_sc: &mut i32,
    num_sc_piii: &mut i32,
    num_sc_asiii: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1723 CountStereoTypes
    // INCHI✔️❌: int CountStereoTypes(INChI* pInChI,
    // INCHI✔️❌:     int* num_known_SB,
    // INCHI✔️❌:     int* num_known_SC,
    // INCHI✔️❌:     int* num_unk_und_SB,
    // INCHI✔️❌:     int* num_unk_und_SC,
    // INCHI✔️❌:     int* num_SC_PIII,
    // INCHI✔️❌:     int* num_SC_AsIII)
    // INCHI✔️❌: {
    // INCHI✔️❌:     INChI_Stereo* Stereo;
    // INCHI✔️❌:     int           i, ret;
    // INCHI✔️❌:     AT_NUMB       nAtNumber;
    // INCHI✔️❌:     U_CHAR        el_number;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!pInChI->nNumberOfAtoms || pInChI->bDeleted)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0; /* no InChI */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     Stereo = (pInChI->StereoIsotopic &&
    // INCHI✔️❌:         (pInChI->StereoIsotopic->nNumberOfStereoBonds +
    // INCHI✔️❌:             pInChI->StereoIsotopic->nNumberOfStereoCenters)) ? pInChI->StereoIsotopic :
    // INCHI✔️❌:         (pInChI->Stereo &&
    // INCHI✔️❌:             (pInChI->Stereo->nNumberOfStereoBonds +
    // INCHI✔️❌:                 pInChI->Stereo->nNumberOfStereoCenters)) ? pInChI->Stereo : NULL;
    // INCHI✔️❌:     if (!Stereo)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 1; /* No Stereo */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* count SB and cumulenes */
    // INCHI✔️❌:     for (i = 0; i < Stereo->nNumberOfStereoBonds; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (ATOM_PARITY_WELL_DEF(Stereo->b_parity[i]))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             (*num_known_SB)++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             (*num_unk_und_SB)++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* count SC and allenes */
    // INCHI✔️❌:     for (i = 0; i < Stereo->nNumberOfStereoCenters; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (!(nAtNumber = Stereo->nNumber[i]) || nAtNumber > pInChI->nNumberOfAtoms)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return RI_ERR_PROGR; /* wrong data, should never happen */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (ATOM_PARITY_WELL_DEF(Stereo->t_parity[i]))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             (*num_known_SC)++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             (*num_unk_und_SC)++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         el_number = pInChI->nAtom[nAtNumber - 1];
    // INCHI✔️❌:         if (el_number != EL_NUMBER_P && el_number != EL_NUMBER_AS)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         ret = GetNumNeighborsFromInchi(pInChI, nAtNumber);
    // INCHI✔️❌:         if (ret < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (3 == ret)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *num_SC_PIII += (EL_NUMBER_P == el_number);
    // INCHI✔️❌:             *num_SC_AsIII += (EL_NUMBER_AS == el_number);
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 2; /* Has Stereo */
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CountStereoTypes
    // BEGIN ACTIVE INCHI MACROS: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/extr_ct.h:294 and util.h:53,57
    // INCHI✔️❌: #define ATOM_PARITY_WELL_DEF(X)     (AB_MIN_WELL_DEFINED_PARITY <= (X) && (X) <= AB_MAX_WELL_DEFINED_PARITY)
    // INCHI✔️❌: #define EL_NUMBER_P  ((U_CHAR)15)
    // INCHI✔️❌: #define EL_NUMBER_AS ((U_CHAR)33)
    // END ACTIVE INCHI MACROS

    const EL_NUMBER_P: u8 = 15;
    const EL_NUMBER_AS: u8 = 33;
    let parity_is_well_defined = |parity: i8| {
        AB_MIN_WELL_DEFINED_PARITY as i32 <= i32::from(parity)
            && i32::from(parity) <= AB_MAX_WELL_DEFINED_PARITY as i32
    };

    if p_inchi.nNumberOfAtoms == 0 || p_inchi.bDeleted != 0 {
        return Ok(0);
    }
    let isotopic = if p_inchi.StereoIsotopic.is_null() {
        None
    } else {
        Some(
            heap.slice(p_inchi.StereoIsotopic.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
    };
    let regular = if p_inchi.Stereo.is_null() {
        None
    } else {
        Some(
            heap.slice(p_inchi.Stereo.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
    };
    let stereo = isotopic
        .filter(|value| {
            value
                .nNumberOfStereoBonds
                .wrapping_add(value.nNumberOfStereoCenters)
                != 0
        })
        .or_else(|| {
            regular.filter(|value| {
                value
                    .nNumberOfStereoBonds
                    .wrapping_add(value.nNumberOfStereoCenters)
                    != 0
            })
        });
    let Some(stereo) = stereo else {
        return Ok(1);
    };

    let stereo_bond_count = usize::try_from(stereo.nNumberOfStereoBonds.max(0))
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if stereo_bond_count > 0 {
        let bond_parities = heap.slice(stereo.b_parity.as_const())?;
        for index in 0..stereo_bond_count {
            let parity = *bond_parities
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if parity_is_well_defined(parity) {
                *num_known_sb = num_known_sb.wrapping_add(1);
            } else {
                *num_unk_und_sb = num_unk_und_sb.wrapping_add(1);
            }
        }
    }

    let stereo_center_count = usize::try_from(stereo.nNumberOfStereoCenters.max(0))
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if stereo_center_count > 0 {
        let atom_numbers = heap.slice(stereo.nNumber.as_const())?;
        let center_parities = heap.slice(stereo.t_parity.as_const())?;
        for index in 0..stereo_center_count {
            let atom_number = *atom_numbers
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if atom_number == 0 || i32::from(atom_number) > p_inchi.nNumberOfAtoms {
                return Ok(RI_ERR_PROGR);
            }
            let parity = *center_parities
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if parity_is_well_defined(parity) {
                *num_known_sc = num_known_sc.wrapping_add(1);
            } else {
                *num_unk_und_sc = num_unk_und_sc.wrapping_add(1);
            }
            let element = *heap
                .slice(p_inchi.nAtom.as_const())?
                .get(usize::from(atom_number.wrapping_sub(1)))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if element != EL_NUMBER_P && element != EL_NUMBER_AS {
                continue;
            }
            let result = GetNumNeighborsFromInchi(heap, Some(p_inchi), atom_number)?;
            if result < 0 {
                return Ok(result);
            }
            if result == 3 {
                *num_sc_piii = num_sc_piii.wrapping_add(i32::from(element == EL_NUMBER_P));
                *num_sc_asiii = num_sc_asiii.wrapping_add(i32::from(element == EL_NUMBER_AS));
            }
        }
    }
    Ok(2)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn DetectInpInchiCreationOptions(
    heap: &SourceHeap,
    one_input: &InpInChI,
    has_reconnected: &mut i32,
    has_metal: &mut i32,
    has_fixed_h: &mut i32,
    mode_flags_stereo: &mut i32,
    taut_flags_stereo: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1880 DetectInpInchiCreationOptions
    // INCHI✔️❌: int DetectInpInchiCreationOptions(InpInChI* pOneInput,
    // INCHI✔️❌:     int* bHasReconnected,
    // INCHI✔️❌:     int* bHasMetal,
    // INCHI✔️❌:     int* bHasFixedH,
    // INCHI✔️❌:     int* nModeFlagsStereo,
    // INCHI✔️❌:     int* bTautFlagsStereo)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0, bHasStereo;
    // INCHI✔️❌:     int nModeFlagsValue = 0, bTautFlagsValue; /* stereo flags */
    // INCHI✔️❌:     int iInChI, iMobileH, bIso, k, max_components, num_components;
    // INCHI✔️❌:     INChI* pInChI;
    // INCHI✔️❌:     int num_known_SB /*Stereo Bonds & Cumulenes >C==C==C==C< */;
    // INCHI✔️❌:     int num_known_SC /* Stereo Centers & Allenes >C=C=C< */;
    // INCHI✔️❌:     int num_unk_und_SB, num_unk_und_SC;
    // INCHI✔️❌:     int num_SC_PIII, num_SC_AsIII; /* has Phosphine or Arsine stereo center(s) */
    // INCHI✔️❌:
    // INCHI✔️❌:     *bHasReconnected = *bHasFixedH = *nModeFlagsStereo = *bTautFlagsStereo = 0;
    // INCHI✔️❌:     nModeFlagsValue = bTautFlagsValue = bHasStereo = 0;
    // INCHI✔️❌:     num_known_SB = num_known_SC = num_unk_und_SB = num_unk_und_SC = num_SC_PIII = num_SC_AsIII = 0;
    // INCHI✔️❌:     *bHasMetal = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (iInChI = 0; iInChI < INCHI_NUM; iInChI++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (iMobileH = 0; iMobileH < TAUT_NUM; iMobileH++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (bIso = 1; !nModeFlagsValue && 0 <= bIso; bIso--)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 switch (pOneInput->s[iInChI][iMobileH][bIso])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                 case 1: /* SABS */
    // INCHI✔️❌:                     nModeFlagsValue |= REQ_MODE_STEREO | REQ_MODE_ISO_STEREO;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case 2:
    // INCHI✔️❌:                     nModeFlagsValue |= REQ_MODE_STEREO | REQ_MODE_ISO_STEREO | REQ_MODE_RELATIVE_STEREO;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case 3:
    // INCHI✔️❌:                     nModeFlagsValue |= REQ_MODE_STEREO | REQ_MODE_ISO_STEREO | REQ_MODE_RACEMIC_STEREO;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             max_components = pOneInput->pInpInChI[iInChI][iMobileH] ?
    // INCHI✔️❌:                 pOneInput->nNumComponents[iInChI][iMobileH] : 0;
    // INCHI✔️❌:
    // INCHI✔️❌:             for (k = num_components = 0; k < max_components; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pInChI = pOneInput->pInpInChI[iInChI][iMobileH] + k;
    // INCHI✔️❌:                 ret = CountStereoTypes(pInChI,
    // INCHI✔️❌:                     &num_known_SB, &num_known_SC,
    // INCHI✔️❌:                     &num_unk_und_SB, &num_unk_und_SC,
    // INCHI✔️❌:                     &num_SC_PIII, &num_SC_AsIII);
    // INCHI✔️❌:                 if (ret < 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return ret; /* error */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 bHasStereo += (ret == 2);
    // INCHI✔️❌:                 if ((ret > 0))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* ret == 0 => Empty InChI, 1=> No Stereo, 2=> Has Stereo */
    // INCHI✔️❌:                     num_components++;
    // INCHI✔️❌:                     *bHasReconnected |= (iInChI == INCHI_REC);
    // INCHI✔️❌:                     *bHasFixedH |= (iMobileH == TAUT_NON);
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 *bHasMetal |= bInChIHasReconnectedMetal(pInChI);
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if ((nModeFlagsValue & REQ_MODE_RELATIVE_STEREO) && (nModeFlagsValue & REQ_MODE_RACEMIC_STEREO))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return RI_ERR_SYNTAX;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (bHasStereo && !nModeFlagsValue) /* REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU*/
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* inversion does not change the stereo or no stereo at all */
    // INCHI✔️❌:         nModeFlagsValue = REQ_MODE_STEREO | REQ_MODE_ISO_STEREO; /* Abs */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!num_known_SB && num_unk_und_SB)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ; /* full SUU option or SB part of it */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nModeFlagsValue |= REQ_MODE_SB_IGN_ALL_UU; /* ignore Unknown/Undefind SB if no well-defined SB exist */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!num_known_SC && num_unk_und_SC)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ; /* full SUU option or SB part of it */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         nModeFlagsValue |= REQ_MODE_SC_IGN_ALL_UU; /* ignore Unknown/Undefind SC if no well-defined SB exist */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* Phosphine and Arsine Stereo */
    // INCHI✔️❌:     if (num_SC_PIII)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         bTautFlagsValue |= TG_FLAG_PHOSPHINE_STEREO;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* Phosphine and Arsine Stereo */
    // INCHI✔️❌:     if (num_SC_AsIII)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         bTautFlagsValue |= TG_FLAG_ARSINE_STEREO;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     *nModeFlagsStereo = nModeFlagsValue;
    // INCHI✔️❌:     *bTautFlagsStereo = bTautFlagsValue;
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: DetectInpInchiCreationOptions

    *has_reconnected = 0;
    *has_fixed_h = 0;
    *mode_flags_stereo = 0;
    *taut_flags_stereo = 0;
    *has_metal = 0;
    let mut mode_flags_value = 0_i32;
    let mut taut_flags_value = 0_i32;
    let mut has_stereo = 0_i32;
    let mut num_known_sb = 0_i32;
    let mut num_known_sc = 0_i32;
    let mut num_unk_und_sb = 0_i32;
    let mut num_unk_und_sc = 0_i32;
    let mut num_sc_piii = 0_i32;
    let mut num_sc_asiii = 0_i32;

    for representation in 0..INCHI_NUM as usize {
        for mobile_h in 0..TAUT_NUM as usize {
            let mut isotopic = 1_i32;
            while mode_flags_value == 0 && isotopic >= 0 {
                match one_input.s[representation][mobile_h][isotopic as usize] {
                    1 => {
                        mode_flags_value |= (REQ_MODE_STEREO | REQ_MODE_ISO_STEREO) as i32;
                    }
                    2 => {
                        mode_flags_value |=
                            (REQ_MODE_STEREO | REQ_MODE_ISO_STEREO | REQ_MODE_RELATIVE_STEREO)
                                as i32;
                    }
                    3 => {
                        mode_flags_value |=
                            (REQ_MODE_STEREO | REQ_MODE_ISO_STEREO | REQ_MODE_RACEMIC_STEREO)
                                as i32;
                    }
                    _ => {}
                }
                isotopic = isotopic.wrapping_sub(1);
            }

            let components = one_input.pInpInChI[representation][mobile_h];
            let max_components = if components.is_null() {
                0
            } else {
                one_input.nNumComponents[representation][mobile_h]
            };
            let mut component = 0_i32;
            let mut num_components = 0_i32;
            while component < max_components {
                let inchi = heap
                    .slice(components.as_const())?
                    .get(
                        usize::try_from(component)
                            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let result = CountStereoTypes(
                    heap,
                    inchi,
                    &mut num_known_sb,
                    &mut num_known_sc,
                    &mut num_unk_und_sb,
                    &mut num_unk_und_sc,
                    &mut num_sc_piii,
                    &mut num_sc_asiii,
                )?;
                if result < 0 {
                    return Ok(result);
                }
                has_stereo = has_stereo.wrapping_add(i32::from(result == 2));
                if result > 0 {
                    num_components = num_components.wrapping_add(1);
                    *has_reconnected |= i32::from(representation == INCHI_REC as usize);
                    *has_fixed_h |= i32::from(mobile_h == TAUT_NON as usize);
                }
                *has_metal |= bInChIHasReconnectedMetal(heap, Some(inchi))?;
                component = component.wrapping_add(1);
            }
            let _ = num_components;
        }
    }

    if mode_flags_value & REQ_MODE_RELATIVE_STEREO as i32 != 0
        && mode_flags_value & REQ_MODE_RACEMIC_STEREO as i32 != 0
    {
        return Ok(RI_ERR_SYNTAX);
    }
    if has_stereo != 0 && mode_flags_value == 0 {
        mode_flags_value = (REQ_MODE_STEREO | REQ_MODE_ISO_STEREO) as i32;
    }
    if num_known_sb != 0 || num_unk_und_sb == 0 {
        mode_flags_value |= REQ_MODE_SB_IGN_ALL_UU as i32;
    }
    if num_known_sc != 0 || num_unk_und_sc == 0 {
        mode_flags_value |= REQ_MODE_SC_IGN_ALL_UU as i32;
    }
    if num_sc_piii != 0 {
        taut_flags_value |= TG_FLAG_PHOSPHINE_STEREO as i32;
    }
    if num_sc_asiii != 0 {
        taut_flags_value |= TG_FLAG_ARSINE_STEREO as i32;
    }
    *mode_flags_stereo = mode_flags_value;
    *taut_flags_stereo = taut_flags_value;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn bRevInchiComponentExists(
    heap: &SourceHeap,
    structure: Option<&StrFromINChI>,
    representation: i32,
    mobile_h: i32,
    component: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1838 bRevInchiComponentExists
    // INCHI✔️❌: int bRevInchiComponentExists(StrFromINChI* pStruct,
    // INCHI✔️❌:     int           iInChI,
    // INCHI✔️❌:     int           bMobileH,
    // INCHI✔️❌:     int           k)
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (!pStruct || /*!pStruct->at2 ||*/ !pStruct->num_atoms ||
    // INCHI✔️❌:         (INCHI_BAS != iInChI && iInChI != INCHI_REC) ||
    // INCHI✔️❌:         (TAUT_NON != bMobileH && TAUT_YES != bMobileH) || k < 0) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return (k < pStruct->RevInChI.num_components[iInChI] &&
    // INCHI✔️❌:         pStruct->RevInChI.pINChI[iInChI] &&
    // INCHI✔️❌:         pStruct->RevInChI.pINChI[iInChI][k][bMobileH] &&
    // INCHI✔️❌:         pStruct->RevInChI.pINChI[iInChI][k][bMobileH]->nNumberOfAtoms > 0 &&
    // INCHI✔️❌:         !pStruct->RevInChI.pINChI[iInChI][k][bMobileH]->bDeleted);
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: bRevInchiComponentExists

    let Some(structure) = structure else {
        return Ok(0);
    };
    if structure.num_atoms == 0
        || (representation != INCHI_BAS as i32 && representation != INCHI_REC as i32)
        || (mobile_h != TAUT_NON as i32 && mobile_h != TAUT_YES as i32)
        || component < 0
    {
        return Ok(0);
    }
    let representation = representation as usize;
    if component >= structure.RevInChI.num_components[representation] {
        return Ok(0);
    }
    let components = structure.RevInChI.pINChI[representation];
    if components.is_null() {
        return Ok(0);
    }
    let inchi = heap
        .slice(components.as_const())?
        .get(usize::try_from(component).map_err(|_| SourceHeapError::SourceIntegerOverflow)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?[mobile_h as usize];
    if inchi.is_null() {
        return Ok(0);
    }
    let inchi = heap
        .slice(inchi.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    Ok(i32::from(inchi.nNumberOfAtoms > 0 && inchi.bDeleted == 0))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputInChIAsRequested(
    heap: &mut SourceHeap,
    canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    out: &mut INCHI_IOSTREAM,
    log: &mut INCHI_IOSTREAM,
    input_parameters: &INPUT_PARMS,
    structure_data: &STRUCT_DATA,
    one_input: &mut InpInChI,
    num_components: &mut [i32; INCHI_NUM as usize],
    proton_modes: &[MODE_PIXH; INCHI_NUM as usize],
    input_number: i64,
    save_option_bits: u8,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1213 OutputInChIAsRequested
    // INCHI✔️❌: int OutputInChIAsRequested(struct tagCANON_GLOBALS* pCG,
    // INCHI✔️❌:     INCHI_IOSTREAM* pOut,
    // INCHI✔️❌:     INCHI_IOSTREAM* pLog,
    // INCHI✔️❌:     ICHICONST INPUT_PARMS* ip_inp,
    // INCHI✔️❌:     STRUCT_DATA* sd_inp,
    // INCHI✔️❌:     InpInChI* OneInput,
    // INCHI✔️❌:     int                     num_components[INCHI_NUM],
    // INCHI✔️❌:     MODE_PIXH               nModeProtonIsoExchgH[INCHI_NUM],
    // INCHI✔️❌:     long                    num_inp,
    // INCHI✔️❌:     unsigned char           save_opt_bits)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int      j, k, k1, k2, ret2 = 0, iINChI, iINChI1; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     PINChI2* pINChI[INCHI_NUM], * newPTR1;
    // INCHI✔️❌:     PINChI_Aux2* pINChI_Aux[INCHI_NUM], * newPTR2;
    // INCHI✔️❌:     int bReqNonTaut;
    // INCHI✔️❌:     int bHasSomeReconnected;
    // INCHI✔️❌:
    // INCHI✔️❌:     INPUT_PARMS ip_local;
    // INCHI✔️❌:     STRUCT_DATA sd_local;
    // INCHI✔️❌:     INPUT_PARMS* ip = &ip_local;
    // INCHI✔️❌:     STRUCT_DATA* sd = &sd_local;
    // INCHI✔️❌:     NORM_CANON_FLAGS ncFlags;
    // INCHI✔️❌:     NORM_CANON_FLAGS* pncFlags = &ncFlags;
    // INCHI✔️❌:     int  nRet1, bSortPrintINChIFlags;
    // INCHI✔️❌:     int  bReqSplitOutputInChI;
    // INCHI✔️❌:     int  nNumOutputComponents;
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_IOS_STRING temp_string_container;
    // INCHI✔️❌:     INCHI_IOS_STRING* strbuf = &temp_string_container;
    // INCHI✔️❌:     memset(strbuf, 0, sizeof(*strbuf)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (0 >= inchi_strbuf_init(strbuf, INCHI_STRBUF_INITIAL_SIZE, INCHI_STRBUF_SIZE_INCREMENT))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret2 = RI_ERR_ALLOC;
    // INCHI✔️❌:         goto exit_error;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nRet1 = 0;
    // INCHI✔️❌:     k1 = k2 = 0;
    // INCHI✔️❌:     memset(pncFlags, 0, sizeof(*pncFlags)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset(pINChI, 0, sizeof(pINChI)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset(pINChI_Aux, 0, sizeof(pINChI_Aux)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:
    // INCHI✔️❌:     *ip = *ip_inp;
    // INCHI✔️❌:     *sd = *sd_inp;
    // INCHI✔️❌:     bHasSomeReconnected = 0;
    // INCHI✔️❌:     bSortPrintINChIFlags = 0;
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     bReqNonTaut = (0 != (ip->nMode & REQ_MODE_BASIC));
    // INCHI✔️❌:     bReqSplitOutputInChI = (0 != (ip->bReadInChIOptions & READ_INCHI_SPLIT_OUTPUT));
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:         if (num_components[INCHI_BAS])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* djb-rwth: MYREALLOC2( PINChI2, PINChI_Aux2, pINChI[INCHI_BAS], pINChI_Aux[INCHI_BAS], num_components[INCHI_BAS], (long long)num_components[INCHI_BAS], k1 ); has been replaced and the whole block rewritten to address memory leaks and reading from freed memory locations */
    // INCHI✔️❌:
    // INCHI✔️❌:             do {
    // INCHI✔️❌:                 if ((num_components[INCHI_BAS]) <= ((long long)num_components[INCHI_BAS]))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     newPTR1 = (PINChI2*)inchi_calloc(((long long)num_components[INCHI_BAS]) + 1, sizeof(PINChI2));
    // INCHI✔️❌:                     newPTR2 = (PINChI_Aux2*)inchi_calloc(((long long)num_components[INCHI_BAS]) + 1, sizeof(PINChI_Aux2));
    // INCHI✔️❌:                     if (newPTR1 && newPTR2) {
    // INCHI✔️❌:                         if ((pINChI[INCHI_BAS]) && (num_components[INCHI_BAS]) > 0)
    // INCHI✔️❌:                             memcpy(newPTR1, (pINChI[INCHI_BAS]), (num_components[INCHI_BAS]) * sizeof(PINChI2));
    // INCHI✔️❌:                         if ((pINChI_Aux[INCHI_BAS]) && (num_components[INCHI_BAS]) > 0)
    // INCHI✔️❌:                             memcpy(newPTR2, (pINChI_Aux[INCHI_BAS]), (num_components[INCHI_BAS]) * sizeof(PINChI_Aux2));
    // INCHI✔️❌:                         if (pINChI[INCHI_BAS])
    // INCHI✔️❌:                             inchi_free(pINChI[INCHI_BAS]);
    // INCHI✔️❌:                         if (pINChI_Aux[INCHI_BAS])
    // INCHI✔️❌:                             inchi_free(pINChI_Aux[INCHI_BAS]);
    // INCHI✔️❌:                         pINChI[INCHI_BAS] = newPTR1;
    // INCHI✔️❌:                         pINChI_Aux[INCHI_BAS] = newPTR2;
    // INCHI✔️❌:                         num_components[INCHI_BAS] = (long long)num_components[INCHI_BAS];
    // INCHI✔️❌:                         k1 = 0;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         inchi_free(newPTR1);
    // INCHI✔️❌:                         inchi_free(newPTR2);
    // INCHI✔️❌:                         k1 = 1;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else { k1 = 0; }
    // INCHI✔️❌:             } while (0);
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (num_components[INCHI_REC])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* djb-rwth: MYREALLOC2( PINChI2, PINChI_Aux2, pINChI[INCHI_REC], pINChI_Aux[INCHI_REC], num_components[INCHI_REC], (long long)num_components[INCHI_REC], k2 ); has been replaced and the whole block rewritten to address memory leaks and reading from freed memory locations */
    // INCHI✔️❌:
    // INCHI✔️❌:         do {
    // INCHI✔️❌:             if ((num_components[INCHI_REC]) <= ((long long)num_components[INCHI_REC]))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 newPTR1 = (PINChI2*)inchi_calloc(((long long)num_components[INCHI_REC]) + 1, sizeof(PINChI2));
    // INCHI✔️❌:                 newPTR2 = (PINChI_Aux2*)inchi_calloc(((long long)num_components[INCHI_REC]) + 1, sizeof(PINChI_Aux2));
    // INCHI✔️❌:                 if (newPTR1 && newPTR2) {
    // INCHI✔️❌:                     if ((pINChI[INCHI_REC]) && (num_components[INCHI_REC]) > 0)
    // INCHI✔️❌:                         memcpy(newPTR1, (pINChI[INCHI_REC]), (num_components[INCHI_REC]) * sizeof(PINChI2));
    // INCHI✔️❌:                     if ((pINChI_Aux[INCHI_REC]) && (num_components[INCHI_REC]) > 0)
    // INCHI✔️❌:                         memcpy(newPTR2, (pINChI_Aux[INCHI_REC]), (num_components[INCHI_REC]) * sizeof(PINChI_Aux2));
    // INCHI✔️❌:                     if (pINChI[INCHI_REC])
    // INCHI✔️❌:                         inchi_free(pINChI[INCHI_REC]);
    // INCHI✔️❌:                     if (pINChI_Aux[INCHI_REC])
    // INCHI✔️❌:                         inchi_free(pINChI_Aux[INCHI_REC]);
    // INCHI✔️❌:                     pINChI[INCHI_REC] = newPTR1;
    // INCHI✔️❌:                     pINChI_Aux[INCHI_REC] = newPTR2;
    // INCHI✔️❌:                     num_components[INCHI_REC] = (long long)num_components[INCHI_REC];
    // INCHI✔️❌:                     k2 = 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     inchi_free(newPTR1);
    // INCHI✔️❌:                     inchi_free(newPTR2);
    // INCHI✔️❌:                     k2 = 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else { k2 = 0; }
    // INCHI✔️❌:         } while (0);
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:         if (k1 || k2 /*|| !pStr*/)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* djb-rwth: avoiding memory leak */
    // INCHI✔️❌:             free(pINChI[INCHI_BAS]);
    // INCHI✔️❌:             free(pINChI_Aux[INCHI_BAS]);
    // INCHI✔️❌:             free(pINChI[INCHI_REC]);
    // INCHI✔️❌:             free(pINChI_Aux[INCHI_REC]);
    // INCHI✔️❌:             ret2 = RI_ERR_ALLOC;
    // INCHI✔️❌:             goto exit_error;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (num_components[INCHI_REC] &&
    // INCHI✔️❌:         (ip->bTautFlags & TG_FLAG_RECONNECT_COORD) &&
    // INCHI✔️❌:         (ip->bTautFlags & TG_FLAG_DISCONNECT_COORD))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sd->bTautFlagsDone[0] |= TG_FLAG_DISCONNECT_COORD_DONE;
    // INCHI✔️❌:         bHasSomeReconnected = 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (iINChI = 0; iINChI < INCHI_NUM; iINChI++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (j = 0; j < TAUT_NUM; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (bReqNonTaut || (j != TAUT_NON && OneInput->pInpInChI[iINChI][j])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (k = 0; k < num_components[iINChI]; k++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* allocate InChI & AuxInfo */
    // INCHI✔️❌:                     if (!(pINChI[iINChI][k][j] = (INChI*)inchi_calloc(1, sizeof(INChI))))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* djb-rwth: avoiding memory leak */
    // INCHI✔️❌:                         free(pINChI[INCHI_BAS]);
    // INCHI✔️❌:                         free(pINChI_Aux[INCHI_BAS]);
    // INCHI✔️❌:                         free(pINChI[INCHI_REC]);
    // INCHI✔️❌:                         free(pINChI_Aux[INCHI_REC]);
    // INCHI✔️❌:                         ret2 = RI_ERR_ALLOC;
    // INCHI✔️❌:                         goto exit_error;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (!(pINChI_Aux[iINChI][k][j] = (INChI_Aux*)inchi_calloc(1, sizeof(INChI_Aux))))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* djb-rwth: avoiding memory leak */
    // INCHI✔️❌:                         free(pINChI[INCHI_BAS]);
    // INCHI✔️❌:                         free(pINChI_Aux[INCHI_BAS]);
    // INCHI✔️❌:                         free(pINChI[INCHI_REC]);
    // INCHI✔️❌:                         free(pINChI_Aux[INCHI_REC]);
    // INCHI✔️❌:                         ret2 = RI_ERR_ALLOC;
    // INCHI✔️❌:                         goto exit_error;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /* copy InChI & AuxInfo */
    // INCHI✔️❌:                     if (k < OneInput->nNumComponents[iINChI][j])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:
    // INCHI✔️❌:                         /* copy InChI */
    // INCHI✔️❌:                         *pINChI[iINChI][k][j] = OneInput->pInpInChI[iINChI][j][k];
    // INCHI✔️❌:                         memset(&OneInput->pInpInChI[iINChI][j][k], 0, sizeof(OneInput->pInpInChI[iINChI][j][k])); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:                         INCHI_HEAPCHK
    // INCHI✔️❌:                             /* take care of protons in AuxInfo */
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (nModeProtonIsoExchgH[iINChI] == MODE_PIXH_ADD_TO_EACH && j == TAUT_YES)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 pINChI_Aux[iINChI][k][j]->nNumRemovedProtons =
    // INCHI✔️❌:                                     OneInput->nNumProtons[iINChI][j].pNumProtons[k].nNumRemovedProtons;
    // INCHI✔️❌:                                 for (k1 = 0; k1 < NUM_H_ISOTOPES; k1++)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     pINChI_Aux[iINChI][k][j]->nNumRemovedIsotopicH[k1] =
    // INCHI✔️❌:                                         OneInput->nNumProtons[iINChI][j].pNumProtons[k].nNumRemovedIsotopicH[k1];
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 INCHI_HEAPCHK
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else if ((!k && nModeProtonIsoExchgH[iINChI] == MODE_PIXH_ADD_TO_FIRST) ||
    // INCHI✔️❌:                                 (k + 1 == OneInput->nNumComponents[iINChI][j] &&
    // INCHI✔️❌:                                     nModeProtonIsoExchgH[iINChI] == MODE_PIXH_ADD_A_PIXH_COMPONENT)) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* add protons and exchangeable isotopic H to the first component's AuxInfo */
    // INCHI✔️❌:                                 pINChI_Aux[iINChI][k][j]->nNumRemovedProtons = OneInput->nNumProtons[iINChI][j].nNumRemovedProtons;
    // INCHI✔️❌:                                 for (k1 = 0; k1 < NUM_H_ISOTOPES; k1++)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     pINChI_Aux[iINChI][k][j]->nNumRemovedIsotopicH[k1] =
    // INCHI✔️❌:                                         OneInput->nNumProtons[iINChI][j].nNumRemovedIsotopicH[k1];
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 INCHI_HEAPCHK
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 pINChI_Aux[iINChI][k][j]->bDeleted = pINChI[iINChI][k][j]->bDeleted;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (j == TAUT_YES && pINChI[iINChI][k][j] && pINChI[iINChI][k][j]->nNumberOfAtoms &&
    // INCHI✔️❌:                             !pINChI[iINChI][k][j]->nNum_H_fixed)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* serializer crashes if it is not allocated */
    // INCHI✔️❌:                             pINChI[iINChI][k][j]->nNum_H_fixed = (S_CHAR*)inchi_calloc((long long)pINChI[iINChI][k][j]->nNumberOfAtoms + 1, sizeof(pINChI[0][0][0]->nNum_H_fixed[0])); /* djb-rwth: cast operator added */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (j == TAUT_YES && k < OneInput->nNumComponents[iINChI][TAUT_NON] &&
    // INCHI✔️❌:                             pINChI[iINChI][k][j] && pINChI[iINChI][k][j]->nNumberOfAtoms &&
    // INCHI✔️❌:                             pINChI[iINChI][k][TAUT_NON] && pINChI[iINChI][k][TAUT_NON]->nNumberOfAtoms &&
    // INCHI✔️❌:                             !CompareReversedINChI(pINChI[iINChI][k][j], pINChI[iINChI][k][TAUT_NON], NULL, NULL))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             pINChI[iINChI][k][TAUT_NON]->nNumberOfAtoms = 0; /* eliminate non-taut equal to taut */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* extra component, usually it is a Mobile H component */
    // INCHI✔️❌:                         /* corresponding to a free proton component in Fixed H */
    // INCHI✔️❌:                         pINChI[iINChI][k][j]->bDeleted = 1;
    // INCHI✔️❌:                         pINChI_Aux[iINChI][k][j]->bDeleted = 1;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 } /* k */
    // INCHI✔️❌:             } /* if ( bReqNonTaut || j != TAUT_NON && OneInput->pInpInChI[iINChI][j] )  */
    // INCHI✔️❌:
    // INCHI✔️❌:             if (OneInput->pInpInChI[iINChI][j])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 INCHI_HEAPCHK
    // INCHI✔️❌:                     inchi_free(OneInput->pInpInChI[iINChI][j]);
    // INCHI✔️❌:                 OneInput->pInpInChI[iINChI][j] = NULL;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         } /* j */
    // INCHI✔️❌:     } /* iINChI */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (bReqSplitOutputInChI)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (bHasSomeReconnected)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             iINChI1 = INCHI_REC; /* only reconnected */
    // INCHI✔️❌:             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:             sd->num_components[INCHI_BAS] = sd->num_components[INCHI_REC];
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             iINChI1 = 0;         /* only disconnected */
    // INCHI✔️❌:             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         sd->num_components[INCHI_REC] = 0;  /* treat reconnected as connected */
    // INCHI✔️❌:         nNumOutputComponents = sd->num_components[INCHI_BAS];
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         iINChI1 = 0;
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         nNumOutputComponents = 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (k1 = 0, k2 = (bReqSplitOutputInChI ? k1 + 1 : nNumOutputComponents); k1 < k2 && k1 < nNumOutputComponents; k1 = k2, k2++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         if (bReqSplitOutputInChI)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sd->num_components[INCHI_BAS] = 1;
    // INCHI✔️❌:             sd->num_components[INCHI_REC] = 0;
    // INCHI✔️❌:             /* additional data */
    // INCHI✔️❌:             sd->num_non_taut[INCHI_BAS] =
    // INCHI✔️❌:                 sd->num_taut[INCHI_BAS] =
    // INCHI✔️❌:                 sd->num_non_taut[INCHI_REC] =
    // INCHI✔️❌:                 sd->num_taut[INCHI_REC] = 0;
    // INCHI✔️❌:             iINChI = iINChI1;
    // INCHI✔️❌:             for (j = 0; j < TAUT_NUM && sd->num_components[iINChI]; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (k = k1; k < k2; k++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /*  find where the current processed structure is located */
    // INCHI✔️❌:                     int cur_is_in_non_taut = (pINChI[iINChI][k][TAUT_NON] && pINChI[iINChI][k][TAUT_NON]->nNumberOfAtoms > 0);
    // INCHI✔️❌:                     int cur_is_in_taut = (pINChI[iINChI][k][TAUT_YES] && pINChI[iINChI][k][TAUT_YES]->nNumberOfAtoms > 0);
    // INCHI✔️❌:                     int cur_is_non_taut = (cur_is_in_non_taut && 0 == pINChI[iINChI][k][TAUT_NON]->lenTautomer) ||
    // INCHI✔️❌:                         (cur_is_in_taut && 0 == pINChI[iINChI][k][TAUT_YES]->lenTautomer); /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️❌:                     int cur_is_taut = cur_is_in_taut && 0 < pINChI[iINChI][k][TAUT_YES]->lenTautomer;
    // INCHI✔️❌:                     if (cur_is_non_taut + cur_is_taut)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /*  count tautomeric and non-tautomeric components of the structures */
    // INCHI✔️❌:                         /*
    // INCHI✔️❌:                         int j1 = cur_is_in_non_taut? TAUT_NON:TAUT_YES;
    // INCHI✔️❌:                         int j2 = cur_is_in_taut?     TAUT_YES:TAUT_NON;
    // INCHI✔️❌:                         */
    // INCHI✔️❌:                         sd->num_non_taut[INCHI_BAS] += cur_is_non_taut;
    // INCHI✔️❌:                         sd->num_taut[INCHI_BAS] += cur_is_taut;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             INCHI_HEAPCHK
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sd->num_components[INCHI_BAS] = inchi_max(OneInput->nNumComponents[INCHI_BAS][TAUT_YES],
    // INCHI✔️❌:                 OneInput->nNumComponents[INCHI_BAS][TAUT_NON]);
    // INCHI✔️❌:             sd->num_components[INCHI_REC] = inchi_max(OneInput->nNumComponents[INCHI_REC][TAUT_YES],
    // INCHI✔️❌:                 OneInput->nNumComponents[INCHI_REC][TAUT_NON]);
    // INCHI✔️❌:             /* additional data needed for SortAndPrintINChI() */
    // INCHI✔️❌:             for (iINChI = 0; iINChI < INCHI_NUM; iINChI++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 sd->num_non_taut[iINChI] =
    // INCHI✔️❌:                     sd->num_taut[iINChI] = 0;
    // INCHI✔️❌:                 for (j = 0; j < TAUT_NUM && sd->num_components[iINChI]; j++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     for (k = k1; k < k2; k++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /*  find where the current processed structure is located */
    // INCHI✔️❌:                         int cur_is_in_non_taut = (pINChI[iINChI][k][TAUT_NON] && pINChI[iINChI][k][TAUT_NON]->nNumberOfAtoms > 0);
    // INCHI✔️❌:                         int cur_is_in_taut = (pINChI[iINChI][k][TAUT_YES] && pINChI[iINChI][k][TAUT_YES]->nNumberOfAtoms > 0);
    // INCHI✔️❌:                         int cur_is_non_taut = (cur_is_in_non_taut && 0 == pINChI[iINChI][k][TAUT_NON]->lenTautomer) ||
    // INCHI✔️❌:                             (cur_is_in_taut && 0 == pINChI[iINChI][k][TAUT_YES]->lenTautomer); /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️❌:                         int cur_is_taut = cur_is_in_taut && 0 < pINChI[iINChI][k][TAUT_YES]->lenTautomer;
    // INCHI✔️❌:                         if (cur_is_non_taut + cur_is_taut)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /*  count tautomeric and non-tautomeric components of the structures */
    // INCHI✔️❌:                             /*
    // INCHI✔️❌:                             int j1 = cur_is_in_non_taut? TAUT_NON:TAUT_YES;
    // INCHI✔️❌:                             int j2 = cur_is_in_taut?     TAUT_YES:TAUT_NON;
    // INCHI✔️❌:                             */
    // INCHI✔️❌:                             sd->num_non_taut[iINChI] += cur_is_non_taut;
    // INCHI✔️❌:                             sd->num_taut[iINChI] += cur_is_taut;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             INCHI_HEAPCHK
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (bReqSplitOutputInChI)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* output components one by one (for splitting input InChI into components) */
    // INCHI✔️❌:             PINChI2* pInChI_2[INCHI_NUM];
    // INCHI✔️❌:             PINChI_Aux2* pInChI_Aux_2[INCHI_NUM];
    // INCHI✔️❌:             INChI* pInChI_1[1][2];
    // INCHI✔️❌:             INChI_Aux* pInChI_Aux_1[1][2];
    // INCHI✔️❌:             memset(pInChI_2, 0, sizeof(pInChI_2)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:             memset(pInChI_Aux_2, 0, sizeof(pInChI_Aux_2)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:             for (j = 0; j < TAUT_NUM; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pInChI_1[0][j] = pINChI[iINChI1][k1][j];
    // INCHI✔️❌:                 pInChI_Aux_1[0][j] = pINChI_Aux[iINChI1][k1][j];
    // INCHI✔️❌:             }
    // INCHI✔️❌:             pInChI_2[INCHI_BAS] = pInChI_1;
    // INCHI✔️❌:             pInChI_Aux_2[INCHI_BAS] = pInChI_Aux_1;
    // INCHI✔️❌:             /* make sure purely reconnected InChI is marked as ReChI, not InChI */
    // INCHI✔️❌:             if (bHasSomeReconnected &&
    // INCHI✔️❌:                 (bInChIHasReconnectedMetal(pInChI_1[0][TAUT_YES]) ||
    // INCHI✔️❌:                     bInChIHasReconnectedMetal(pInChI_1[0][TAUT_NON])))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 bSortPrintINChIFlags = FLAG_SORT_PRINT_ReChI_PREFIX;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 bSortPrintINChIFlags = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             INCHI_HEAPCHK
    // INCHI✔️❌:                 nRet1 = SortAndPrintINChI(pCG, pOut, strbuf, pLog, ip,
    // INCHI✔️❌:                     NULL /*orig_inp_data*/,
    // INCHI✔️❌:                     NULL  /*prep_inp_data*/,
    // INCHI✔️❌:                     NULL /*composite_norm_data*/,
    // INCHI✔️❌:                     NULL /*pOrigStruct*/,
    // INCHI✔️❌:                     sd->num_components, sd->num_non_taut,
    // INCHI✔️❌:                     sd->num_taut, sd->bTautFlags,
    // INCHI✔️❌:                     sd->bTautFlagsDone, pncFlags, num_inp,
    // INCHI✔️❌:                     pInChI_2, pInChI_Aux_2,
    // INCHI✔️❌:                     &bSortPrintINChIFlags, save_opt_bits);
    // INCHI✔️❌:             INCHI_HEAPCHK
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:                 bSortPrintINChIFlags = 0;
    // INCHI✔️❌:             nRet1 = SortAndPrintINChI(pCG, pOut, strbuf, pLog, ip,
    // INCHI✔️❌:                 NULL /*orig_inp_data*/, NULL  /*prep_inp_data*/,
    // INCHI✔️❌:                 NULL /*composite_norm_data*/, NULL /*pOrigStruct*/,
    // INCHI✔️❌:                 sd->num_components, sd->num_non_taut, sd->num_taut,
    // INCHI✔️❌:                 sd->bTautFlags, sd->bTautFlagsDone, pncFlags, num_inp,
    // INCHI✔️❌:                 pINChI, pINChI_Aux, &bSortPrintINChIFlags,
    // INCHI✔️❌:                 save_opt_bits);
    // INCHI✔️❌:             INCHI_HEAPCHK
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (nRet1 == _IS_FATAL || nRet1 == _IS_ERROR)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:         FreeAllINChIArrays(pINChI, pINChI_Aux, num_components);
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:         for (iINChI = 0; iINChI < INCHI_NUM; iINChI++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (j = 0; j < TAUT_NUM; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (OneInput->nNumProtons[iINChI][j].pNumProtons)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     inchi_free(OneInput->nNumProtons[iINChI][j].pNumProtons);
    // INCHI✔️❌:                     OneInput->nNumProtons[iINChI][j].pNumProtons = NULL;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_HEAPCHK
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nRet1 == _IS_FATAL || nRet1 == _IS_ERROR)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret2 = RI_ERR_PROGR;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_error:
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_strbuf_close(strbuf);
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret2;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OutputInChIAsRequested

    let mut string_buffer = INCHI_IOS_STRING::default();
    let initialization = inchi_strbuf_init(
        heap,
        &mut string_buffer,
        INCHI_STRBUF_INITIAL_SIZE as i32,
        INCHI_STRBUF_SIZE_INCREMENT as i32,
    )?;
    if initialization <= 0 {
        inchi_strbuf_close(heap, Some(&mut string_buffer))?;
        return Ok(RI_ERR_ALLOC);
    }

    let execution = (|| -> Result<i32, SourceHeapError> {
        let mut ip = input_parameters.clone();
        let mut sd = structure_data.clone();
        let normalization_flags = NORM_CANON_FLAGS::default();
        let mut p_inchi = [SourceMutPointer::<PINChI2>::null(); INCHI_NUM as usize];
        let mut p_inchi_aux = [SourceMutPointer::<PINChI_Aux2>::null(); INCHI_NUM as usize];
        let requested_non_tautomeric = ip.nMode & REQ_MODE_BASIC as INCHI_MODE != 0;
        let split_output = ip.bReadInChIOptions & READ_INCHI_SPLIT_OUTPUT as i32 != 0;

        let mut allocation_failed = false;
        for representation in 0..INCHI_NUM as usize {
            if num_components[representation] == 0 {
                continue;
            }
            let allocation_count = i64::from(num_components[representation]).wrapping_add(1);
            let count = match usize::try_from(allocation_count) {
                Ok(count) => count,
                Err(_) => {
                    allocation_failed = true;
                    continue;
                }
            };
            let new_inchi = inchi_calloc::<PINChI2>(heap, count as u64, 1);
            let new_aux = inchi_calloc::<PINChI_Aux2>(heap, count as u64, 1);
            match (new_inchi, new_aux) {
                (Ok(inchi), Ok(aux)) => {
                    p_inchi[representation] = inchi;
                    p_inchi_aux[representation] = aux;
                }
                (inchi, aux) => {
                    if let Ok(pointer) = inchi {
                        inchi_free(heap, pointer)?;
                    }
                    if let Ok(pointer) = aux {
                        inchi_free(heap, pointer)?;
                    }
                    allocation_failed = true;
                }
            }
        }
        if allocation_failed {
            for representation in 0..INCHI_NUM as usize {
                if !p_inchi[representation].is_null() {
                    inchi_free(heap, p_inchi[representation])?;
                }
                if !p_inchi_aux[representation].is_null() {
                    inchi_free(heap, p_inchi_aux[representation])?;
                }
            }
            return Ok(RI_ERR_ALLOC);
        }

        let has_reconnected = num_components[INCHI_REC as usize] != 0
            && ip.bTautFlags & TG_FLAG_RECONNECT_COORD as INCHI_MODE != 0
            && ip.bTautFlags & TG_FLAG_DISCONNECT_COORD as INCHI_MODE != 0;
        if has_reconnected {
            sd.bTautFlagsDone[INCHI_BAS as usize] |= TG_FLAG_DISCONNECT_COORD_DONE as INCHI_MODE;
        }

        for representation in 0..INCHI_NUM as usize {
            for tautomer in 0..TAUT_NUM as usize {
                let source_pointer = one_input.pInpInChI[representation][tautomer];
                if requested_non_tautomeric
                    || (tautomer != TAUT_NON as usize && !source_pointer.is_null())
                {
                    let count = if num_components[representation] > 0 {
                        usize::try_from(num_components[representation])
                            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                    } else {
                        0
                    };
                    for component in 0..count {
                        let inchi_owner = match inchi_calloc::<INChI>(heap, 1, 1) {
                            Ok(pointer) => pointer,
                            Err(SourceHeapError::AllocationFailed) => {
                                for index in 0..INCHI_NUM as usize {
                                    if !p_inchi[index].is_null() {
                                        inchi_free(heap, p_inchi[index])?;
                                    }
                                    if !p_inchi_aux[index].is_null() {
                                        inchi_free(heap, p_inchi_aux[index])?;
                                    }
                                }
                                return Ok(RI_ERR_ALLOC);
                            }
                            Err(error) => return Err(error),
                        };
                        heap.slice_mut(p_inchi[representation])?[component][tautomer] = inchi_owner;

                        let aux_owner = match inchi_calloc::<INChI_Aux>(heap, 1, 1) {
                            Ok(pointer) => pointer,
                            Err(SourceHeapError::AllocationFailed) => {
                                for index in 0..INCHI_NUM as usize {
                                    if !p_inchi[index].is_null() {
                                        inchi_free(heap, p_inchi[index])?;
                                    }
                                    if !p_inchi_aux[index].is_null() {
                                        inchi_free(heap, p_inchi_aux[index])?;
                                    }
                                }
                                return Ok(RI_ERR_ALLOC);
                            }
                            Err(error) => return Err(error),
                        };
                        heap.slice_mut(p_inchi_aux[representation])?[component][tautomer] =
                            aux_owner;

                        if component
                            < usize::try_from(
                                one_input.nNumComponents[representation][tautomer].max(0),
                            )
                            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                        {
                            let moved = heap
                                .slice(source_pointer.as_const())?
                                .get(component)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            heap.slice_mut(inchi_owner)?[0] = moved;
                            heap.slice_mut(source_pointer)?[component] = INChI::default();

                            let moved_inchi = heap.slice(inchi_owner.as_const())?[0].clone();
                            let mut aux = INChI_Aux::default();
                            if proton_modes[representation] == MODE_PIXH_ADD_TO_EACH
                                && tautomer == TAUT_YES as usize
                            {
                                let proton_pointer =
                                    one_input.nNumProtons[representation][tautomer].pNumProtons;
                                let protons = heap
                                    .slice(proton_pointer.as_const())?
                                    .get(component)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                aux.nNumRemovedProtons = protons.nNumRemovedProtons;
                                aux.nNumRemovedIsotopicH = protons.nNumRemovedIsotopicH;
                            } else if (component == 0
                                && proton_modes[representation] == MODE_PIXH_ADD_TO_FIRST)
                                || (component.wrapping_add(1)
                                    == usize::try_from(
                                        one_input.nNumComponents[representation][tautomer].max(0),
                                    )
                                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                                    && proton_modes[representation]
                                        == MODE_PIXH_ADD_A_PIXH_COMPONENT)
                            {
                                aux.nNumRemovedProtons = one_input.nNumProtons[representation]
                                    [tautomer]
                                    .nNumRemovedProtons;
                                aux.nNumRemovedIsotopicH = one_input.nNumProtons[representation]
                                    [tautomer]
                                    .nNumRemovedIsotopicH;
                            } else {
                                aux.bDeleted = moved_inchi.bDeleted;
                            }
                            heap.slice_mut(aux_owner)?[0] = aux;

                            if tautomer == TAUT_YES as usize
                                && moved_inchi.nNumberOfAtoms != 0
                                && moved_inchi.nNum_H_fixed.is_null()
                            {
                                let fixed_count =
                                    i64::from(moved_inchi.nNumberOfAtoms).wrapping_add(1);
                                if let Ok(fixed_count) = u64::try_from(fixed_count) {
                                    if let Ok(pointer) = inchi_calloc::<i8>(heap, fixed_count, 1) {
                                        heap.slice_mut(inchi_owner)?[0].nNum_H_fixed = pointer;
                                    }
                                }
                            }

                            if tautomer == TAUT_YES as usize
                                && component
                                    < usize::try_from(
                                        one_input.nNumComponents[representation][TAUT_NON as usize]
                                            .max(0),
                                    )
                                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                            {
                                let non_taut_owner = heap
                                    .slice(p_inchi[representation].as_const())?[component]
                                    [TAUT_NON as usize];
                                let mobile = heap.slice(inchi_owner.as_const())?[0].clone();
                                let non_taut = if non_taut_owner.is_null() {
                                    None
                                } else {
                                    Some(heap.slice(non_taut_owner.as_const())?[0].clone())
                                };
                                if mobile.nNumberOfAtoms != 0
                                    && non_taut
                                        .as_ref()
                                        .is_some_and(|value| value.nNumberOfAtoms != 0)
                                    && CompareReversedINChI(
                                        heap,
                                        Some(&mobile),
                                        non_taut.as_ref(),
                                        None,
                                        None,
                                    )? == 0
                                {
                                    heap.slice_mut(non_taut_owner)?[0].nNumberOfAtoms = 0;
                                }
                            }
                        } else {
                            heap.slice_mut(inchi_owner)?[0].bDeleted = 1;
                            heap.slice_mut(aux_owner)?[0].bDeleted = 1;
                        }
                    }
                }

                if !source_pointer.is_null() {
                    inchi_free(heap, source_pointer)?;
                    one_input.pInpInChI[representation][tautomer] = SourceMutPointer::null();
                }
            }
        }

        let (selected_representation, output_components) = if split_output {
            let selected = if has_reconnected {
                sd.num_components[INCHI_BAS as usize] = sd.num_components[INCHI_REC as usize];
                INCHI_REC as usize
            } else {
                INCHI_BAS as usize
            };
            sd.num_components[INCHI_REC as usize] = 0;
            (selected, sd.num_components[INCHI_BAS as usize])
        } else {
            (INCHI_BAS as usize, 1)
        };

        let mut output_result = 0_i32;
        let mut first = 0_i32;
        let mut end = if split_output { 1 } else { output_components };
        while first < end && first < output_components {
            if split_output {
                sd.num_components = [1, 0];
                sd.num_non_taut = [0, 0];
                sd.num_taut = [0, 0];
                let component =
                    usize::try_from(first).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                if sd.num_components[selected_representation] != 0 {
                    for _tautomer in 0..TAUT_NUM as usize {
                        let row = heap.slice(p_inchi[selected_representation].as_const())?;
                        let non_taut_owner = row[component][TAUT_NON as usize];
                        let taut_owner = row[component][TAUT_YES as usize];
                        let non_taut = if non_taut_owner.is_null() {
                            None
                        } else {
                            Some(heap.slice(non_taut_owner.as_const())?[0].clone())
                        };
                        let taut = if taut_owner.is_null() {
                            None
                        } else {
                            Some(heap.slice(taut_owner.as_const())?[0].clone())
                        };
                        let is_in_non_taut = non_taut
                            .as_ref()
                            .is_some_and(|value| value.nNumberOfAtoms > 0);
                        let is_in_taut =
                            taut.as_ref().is_some_and(|value| value.nNumberOfAtoms > 0);
                        let is_non_taut = is_in_non_taut
                            && non_taut
                                .as_ref()
                                .is_some_and(|value| value.lenTautomer == 0)
                            || is_in_taut
                                && taut.as_ref().is_some_and(|value| value.lenTautomer == 0);
                        let is_taut =
                            is_in_taut && taut.as_ref().is_some_and(|value| value.lenTautomer > 0);
                        if is_non_taut || is_taut {
                            sd.num_non_taut[INCHI_BAS as usize] = sd.num_non_taut
                                [INCHI_BAS as usize]
                                .wrapping_add(i32::from(is_non_taut));
                            sd.num_taut[INCHI_BAS as usize] =
                                sd.num_taut[INCHI_BAS as usize].wrapping_add(i32::from(is_taut));
                        }
                    }
                }
            } else {
                for representation in 0..INCHI_NUM as usize {
                    sd.num_components[representation] = one_input.nNumComponents[representation]
                        [TAUT_YES as usize]
                        .max(one_input.nNumComponents[representation][TAUT_NON as usize]);
                    sd.num_non_taut[representation] = 0;
                    sd.num_taut[representation] = 0;
                    if sd.num_components[representation] == 0 {
                        continue;
                    }
                    let component = usize::try_from(first)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    for _tautomer in 0..TAUT_NUM as usize {
                        let row = heap.slice(p_inchi[representation].as_const())?;
                        let non_taut_owner = row[component][TAUT_NON as usize];
                        let taut_owner = row[component][TAUT_YES as usize];
                        let non_taut = if non_taut_owner.is_null() {
                            None
                        } else {
                            Some(heap.slice(non_taut_owner.as_const())?[0].clone())
                        };
                        let taut = if taut_owner.is_null() {
                            None
                        } else {
                            Some(heap.slice(taut_owner.as_const())?[0].clone())
                        };
                        let is_in_non_taut = non_taut
                            .as_ref()
                            .is_some_and(|value| value.nNumberOfAtoms > 0);
                        let is_in_taut =
                            taut.as_ref().is_some_and(|value| value.nNumberOfAtoms > 0);
                        let is_non_taut = is_in_non_taut
                            && non_taut
                                .as_ref()
                                .is_some_and(|value| value.lenTautomer == 0)
                            || is_in_taut
                                && taut.as_ref().is_some_and(|value| value.lenTautomer == 0);
                        let is_taut =
                            is_in_taut && taut.as_ref().is_some_and(|value| value.lenTautomer > 0);
                        if is_non_taut || is_taut {
                            sd.num_non_taut[representation] = sd.num_non_taut[representation]
                                .wrapping_add(i32::from(is_non_taut));
                            sd.num_taut[representation] =
                                sd.num_taut[representation].wrapping_add(i32::from(is_taut));
                        }
                    }
                }
            }

            let flags = heap.allocate_model_storage(vec![0_i32])?;
            if split_output {
                let component =
                    usize::try_from(first).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let row = heap.slice(p_inchi[selected_representation].as_const())?[component];
                let aux_row =
                    heap.slice(p_inchi_aux[selected_representation].as_const())?[component];
                let split_rows = heap.allocate_model_storage(vec![row])?;
                let split_aux_rows = heap.allocate_model_storage(vec![aux_row])?;
                let has_reconnected_metal = if has_reconnected {
                    let mobile = if row[TAUT_YES as usize].is_null() {
                        None
                    } else {
                        Some(heap.slice(row[TAUT_YES as usize].as_const())?[0].clone())
                    };
                    let fixed = if row[TAUT_NON as usize].is_null() {
                        None
                    } else {
                        Some(heap.slice(row[TAUT_NON as usize].as_const())?[0].clone())
                    };
                    bInChIHasReconnectedMetal(heap, mobile.as_ref())? != 0
                        || bInChIHasReconnectedMetal(heap, fixed.as_ref())? != 0
                } else {
                    false
                };
                heap.slice_mut(flags)?[0] = if has_reconnected_metal {
                    FLAG_SORT_PRINT_ReChI_PREFIX as i32
                } else {
                    0
                };
                output_result = SortAndPrintINChI(
                    heap,
                    canonical_globals,
                    out,
                    Some(&mut string_buffer),
                    log,
                    &ip,
                    SourceConstPointer::null(),
                    SourceConstPointer::null(),
                    None,
                    SourceConstPointer::null(),
                    &sd.num_components,
                    &sd.num_non_taut,
                    &sd.num_taut,
                    &mut sd.bTautFlags,
                    &mut sd.bTautFlagsDone,
                    &normalization_flags,
                    input_number,
                    [split_rows, SourceMutPointer::null()],
                    [split_aux_rows, SourceMutPointer::null()],
                    flags,
                    save_option_bits,
                    stdout,
                )?;
                heap.free(split_rows)?;
                heap.free(split_aux_rows)?;
            } else {
                output_result = SortAndPrintINChI(
                    heap,
                    canonical_globals,
                    out,
                    Some(&mut string_buffer),
                    log,
                    &ip,
                    SourceConstPointer::null(),
                    SourceConstPointer::null(),
                    None,
                    SourceConstPointer::null(),
                    &sd.num_components,
                    &sd.num_non_taut,
                    &sd.num_taut,
                    &mut sd.bTautFlags,
                    &mut sd.bTautFlagsDone,
                    &normalization_flags,
                    input_number,
                    p_inchi,
                    p_inchi_aux,
                    flags,
                    save_option_bits,
                    stdout,
                )?;
            }
            heap.free(flags)?;
            if output_result == _IS_FATAL as i32 || output_result == _IS_ERROR as i32 {
                break;
            }
            first = end;
            end = end.wrapping_add(1);
        }

        FreeAllINChIArrays(heap, &mut p_inchi, &mut p_inchi_aux, num_components)?;
        for representation in 0..INCHI_NUM as usize {
            for tautomer in 0..TAUT_NUM as usize {
                let pointer = one_input.nNumProtons[representation][tautomer].pNumProtons;
                if !pointer.is_null() {
                    inchi_free(heap, pointer)?;
                    one_input.nNumProtons[representation][tautomer].pNumProtons =
                        SourceMutPointer::null();
                }
            }
        }
        Ok(
            if output_result == _IS_FATAL as i32 || output_result == _IS_ERROR as i32 {
                RI_ERR_PROGR
            } else {
                0
            },
        )
    })();

    let close_result = inchi_strbuf_close(heap, Some(&mut string_buffer));
    let result = execution?;
    close_result?;
    Ok(result)
}

#[allow(non_snake_case)]
pub(crate) fn bInChIHasReconnectedMetal(
    heap: &SourceHeap,
    inchi: Option<&INChI>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1993 bInChIHasReconnectedMetal
    // INCHI✔❌: int bInChIHasReconnectedMetal(INChI* pInChI)
    // INCHI✔❌: {
    // INCHI✔❌:     int i;
    // INCHI✔❌:     if (pInChI && !pInChI->bDeleted && pInChI->nNumberOfAtoms && pInChI->nAtom)
    // INCHI✔❌:     {
    // INCHI✔❌:         for (i = 0; i < pInChI->nNumberOfAtoms; i++)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (is_el_a_metal((int)pInChI->nAtom[i]))
    // INCHI✔❌:             {
    // INCHI✔❌:                 if (pInChI->nNumberOfAtoms > 1 || (pInChI->nNum_H && pInChI->nNum_H[0])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:                 {
    // INCHI✔❌:                     return 1;
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return 0;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: bInChIHasReconnectedMetal

    let Some(inchi) = inchi else {
        return Ok(0);
    };
    if inchi.bDeleted != 0 || inchi.nNumberOfAtoms == 0 || inchi.nAtom.is_null() {
        return Ok(0);
    }
    let count = if inchi.nNumberOfAtoms > 0 {
        usize::try_from(inchi.nNumberOfAtoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
    } else {
        0
    };
    let atoms = heap
        .slice(inchi.nAtom.as_const())?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for &atom in atoms {
        if is_el_a_metal(i32::from(atom))? != 0 {
            if inchi.nNumberOfAtoms > 1 {
                return Ok(1);
            }
            if !inchi.nNum_H.is_null() && heap.slice(inchi.nNum_H.as_const())?[0] != 0 {
                return Ok(1);
            }
        }
    }
    Ok(0)
}

// INCHI source: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:380-438 irErrMsg
const IR_ERR_MSG: [&str; 37] = [
    "MOBILE_H_FORMULA",
    "MOBILE_H_CONNECTIONS",
    "MOBILE_H",
    "MOBILE_H_CHARGE",
    "MOBILE_H_PROTONS",
    "MOBILE_H_SP2",
    "MOBILE_H_SP3",
    "MOBILE_H_SP3_/m",
    "MOBILE_H_SP3_/s",
    "MOBILE_H_ISO_LAYER_FORK",
    "MOBILE_H_ISO_ATOMS",
    "MOBILE_H_ISO_EXCH_H",
    "MOBILE_H_ISO_SP2",
    "MOBILE_H_ISO_SP3",
    "MOBILE_H_ISO_SP3_/m",
    "MOBILE_H_ISO_SP3_/s",
    "FIXED_H_LAYER_FORK",
    "FIXED_H_FORMULA",
    "FIXED_H",
    "FIXED_H_CHARGE",
    "FIXED_H_SP2",
    "FIXED_H_SP3",
    "FIXED_H_SP3_/m",
    "FIXED_H_SP3_/s",
    "FIXED_H_PERMUTATION",
    "FIXED_H_ISO_LAYER_FORK",
    "FIXED_H_ISO_ATOMS",
    "FIXED_H_ISO_LAYER",
    "FIXED_H_ISO_SP2",
    "FIXED_H_ISO_SP3",
    "FIXED_H_ISO_SP3_m",
    "FIXED_H_ISO_SP3_s",
    "FIXED_H_ISO_PERMUTATION",
    "RECONNECTED_LAYER_FORK",
    "RECONNECTED_FORMULA",
    "MATERIAL_BALANCE",
    "POLYMER_LAYER",
];

#[allow(non_snake_case)]
pub(crate) fn getInchiStateReadErr(
    heap: &mut SourceHeap,
    mut stat: i32,
    sz_msg: SourceMutPointer<i8>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:536 getInchiStateReadErr
    // INCHI✔❌: void getInchiStateReadErr(int stat, char* szMsg)
    // INCHI✔❌: /* const char *getInchiStateReadErr(int stat) */
    // INCHI✔❌: {
    // INCHI✔❌:     int i, bRecMet = 0;
    // INCHI✔❌:
    // INCHI✔❌:     if (stat >= IST_HAPPENED_IN_RECMET)
    // INCHI✔❌:     {
    // INCHI✔❌:         bRecMet = 1;
    // INCHI✔❌:         stat -= IST_HAPPENED_IN_RECMET;
    // INCHI✔❌:     }
    // INCHI✔❌:     for (i = 0; 0 <= irErrMsg[i].stat && stat != irErrMsg[i].stat; i++)
    // INCHI✔❌:     {
    // INCHI✔❌:         ;
    // INCHI✔❌:     }
    // INCHI✔❌:     sprintf(szMsg,
    // INCHI✔❌: #if ( FIX_DALKE_BUGS == 1 )
    // INCHI✔❌:         "%s%.100s",
    // INCHI✔❌: #else
    // INCHI✔❌:         "%s%s",
    // INCHI✔❌: #endif
    // INCHI✔❌:         irErrMsg[i].msg, bRecMet ? ", Reconnected layer" : "");
    // INCHI✔❌:
    // INCHI✔❌: }
    // END INCHI C FUNCTION: getInchiStateReadErr

    let reconnected = stat >= 100;
    if reconnected {
        stat -= 100;
    }
    let message = usize::try_from(stat)
        .ok()
        .and_then(|index| IR_ERR_MSG.get(index))
        .copied()
        .unwrap_or("Unknown Error");
    let suffix = if reconnected {
        ", Reconnected layer"
    } else {
        ""
    };
    let required = message
        .len()
        .checked_add(suffix.len())
        .and_then(|length| length.checked_add(1))
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    let destination = heap.slice_mut(sz_msg)?;
    if destination.len() < required {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let mut position = 0;
    for byte in message.bytes().chain(suffix.bytes()) {
        destination[position] = byte as i8;
        position += 1;
    }
    destination[position] = 0;
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn getInchiErrName(n_err: i32) -> &'static [u8] {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:562 getInchiErrName
    // INCHI✔❌: const char* getInchiErrName(int nErr)
    // INCHI✔❌: {
    // INCHI✔❌:     switch (nErr)
    // INCHI✔❌:     {
    // INCHI✔❌:     case RI_ERR_ALLOC:
    // INCHI✔❌:         return "Allocation failed";
    // INCHI✔❌:     case RI_ERR_PROGR:
    // INCHI✔❌:         return "Program error";
    // INCHI✔❌:     case RI_ERR_SYNTAX:
    // INCHI✔❌:         return "Syntax error";
    // INCHI✔❌:     case RI_ERR_EOL:
    // INCHI✔❌:         return "End of line";
    // INCHI✔❌:     }
    // INCHI✔❌:     return "Unknown error";
    // INCHI✔❌: }
    // END INCHI C FUNCTION: getInchiErrName

    match n_err {
        RI_ERR_ALLOC => b"Allocation failed\0",
        RI_ERR_PROGR => b"Program error\0",
        RI_ERR_SYNTAX => b"Syntax error\0",
        RI_ERR_EOL => b"End of line\0",
        _ => b"Unknown error\0",
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn TreatErrorsInReadInChIString(
    heap: &mut SourceHeap,
    n_read_status: i32,
    n_err: i32,
    p_state: i32,
    ip: &INPUT_PARMS,
    mut p_out: Option<&mut INCHI_IOSTREAM>,
    mut p_log: Option<&mut INCHI_IOSTREAM>,
    num_inp: &mut i64,
    num_errors: &mut i64,
    num_processed: &mut i64,
    pstr_hdr: &mut SourceMutPointer<i8>,
    psz_cur_hdr: &mut SourceMutPointer<i8>,
    p_one_input: &mut InpInChI,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10901 TreatErrorsInReadInChIString
    // INCHI✔❌: void TreatErrorsInReadInChIString(int            nReadStatus,
    // INCHI✔❌:     int            nErr,
    // INCHI✔❌:     int            pState,
    // INCHI✔❌:     INPUT_PARMS* ip,
    // INCHI✔❌:     INCHI_IOSTREAM* pOut,
    // INCHI✔❌:     INCHI_IOSTREAM* pLog,
    // INCHI✔❌:     long* num_inp,
    // INCHI✔❌:     long* num_errors,
    // INCHI✔❌:     long* num_processed,
    // INCHI✔❌:     char** pstrHdr,
    // INCHI✔❌:     char** pszCurHdr,
    // INCHI✔❌:     InpInChI* pOneInput)
    // INCHI✔❌: {
    // INCHI✔❌:     int bInChI2Struct = (ip->bReadInChIOptions & READ_INCHI_TO_STRUCTURE) && ip->nInputType == INPUT_INCHI;
    // INCHI✔❌:
    // INCHI✔❌:     /* InChI could not be read */
    // INCHI✔❌:     if (nReadStatus == RI_ERR_EOF && nErr == 0 && pState == 0) /* && !(*pstrHdr) )  */
    // INCHI✔❌:     {
    // INCHI✔❌:         /*if ( !(*pstrHdr) ) */
    // INCHI✔❌:         ;/*inchi_ios_eprint( pLog, "\nEnd of file detected after structure %ld.    \n", *num_inp );*/
    // INCHI✔❌:     }
    // INCHI✔❌:     else
    // INCHI✔❌:     {
    // INCHI✔❌:         /* Output InChI parsing error message */
    // INCHI✔❌:         char szHdrSimulation[128];
    // INCHI✔❌:         char szMsg2[1024];
    // INCHI✔❌:         (*num_inp)++;
    // INCHI✔❌:         sprintf(szHdrSimulation, "Structure: %ld", *num_inp);
    // INCHI✔❌:         getInchiStateReadErr(pState, szMsg2);
    // INCHI✔❌:
    // INCHI✔❌: #ifdef TARGET_EXE_STANDALONE
    // INCHI✔❌:         if (pOneInput->polymer &&
    // INCHI✔❌:             bInChI2Struct &&
    // INCHI✔❌:             !(ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY))
    // INCHI✔❌:         {
    // INCHI✔❌:             inchi_ios_eprint(pLog, "%s Skipping polymer InChI (only conversion to Molfile is available, use OutputSDF option)\n",
    // INCHI✔❌:                 *pstrHdr ? *pstrHdr : szHdrSimulation);
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌: #endif
    // INCHI✔❌:         {
    // INCHI✔❌:             if (!bInChI2Struct &&
    // INCHI✔❌:                 (pState == IST_MOBILE_H_POLYMER && !ip->bPolymers))
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* TO DO: implement InChI2InChI for polymers in a way similar to InChI2Struct
    // INCHI✔❌:                 thru an external (to ReadWriteInchi) loop                                */
    // INCHI✔❌:                 inchi_ios_eprint(pLog, "%s Skipping polymer InChI for conversion of InChI to InChI\n",
    // INCHI✔❌:                     *pstrHdr ? *pstrHdr : szHdrSimulation);
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 inchi_ios_eprint(pLog, "\n%s %s (%d) in %s (%d)\n",
    // INCHI✔❌:                     *pstrHdr ? *pstrHdr : szHdrSimulation,
    // INCHI✔❌:                     getInchiErrName(nErr), nErr,
    // INCHI✔❌:                     szMsg2, pState);
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         if (ip->bINChIOutputOptions2 & INCHI_OUT_INCHI_GEN_ERROR)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (!(ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY))
    // INCHI✔❌:             {
    // INCHI✔❌:                 inchi_ios_eprint(pOut, "%s\n", *pstrHdr ? *pstrHdr : szHdrSimulation);
    // INCHI✔❌:                 if (ip->bINChIOutputOptions & INCHI_OUT_STDINCHI)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     inchi_ios_eprint(pOut, "InChI=1S//\n");
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                     inchi_ios_eprint(pOut, "InChI=1//\n");
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         if (0 != (ip->bReadInChIOptions & READ_INCHI_TO_STRUCTURE))
    // INCHI✔❌:
    // INCHI✔❌:             (*num_errors)++;
    // INCHI✔❌:         (*num_processed)++;
    // INCHI✔❌:     }
    // INCHI✔❌:     if (*pstrHdr)
    // INCHI✔❌:     {
    // INCHI✔❌:         inchi_free(*pstrHdr);
    // INCHI✔❌:         *pstrHdr = NULL;
    // INCHI✔❌:     }
    // INCHI✔❌:     if (*pszCurHdr)
    // INCHI✔❌:     {
    // INCHI✔❌:         inchi_free(*pszCurHdr);
    // INCHI✔❌:         *pszCurHdr = NULL;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     FreeInpInChI(pOneInput);
    // INCHI✔❌:
    // INCHI✔❌:     return;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: TreatErrorsInReadInChIString

    let b_inchi_to_struct = ip.bReadInChIOptions & READ_INCHI_TO_STRUCTURE as i32 != 0
        && ip.nInputType == tagInputType_INPUT_INCHI;
    if !(n_read_status == RI_ERR_EOF as i32 && n_err == 0 && p_state == 0) {
        *num_inp = num_inp
            .checked_add(1)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        let simulation = format!("Structure: {num_inp}");
        let simulation = heap.allocate_model_storage(
            simulation
                .bytes()
                .chain(std::iter::once(0))
                .map(|byte| byte as i8)
                .collect(),
        )?;
        let state_message = heap.allocate_model_storage(vec![0_i8; 1024])?;
        getInchiStateReadErr(heap, p_state, state_message)?;
        let header = if pstr_hdr.is_null() {
            simulation.as_const()
        } else {
            pstr_hdr.as_const()
        };

        if !b_inchi_to_struct && p_state == 36 && ip.bPolymers == 0 {
            let _ = eprint_call(
                heap,
                p_log.as_deref_mut(),
                "%s Skipping polymer InChI for conversion of InChI to InChI\n",
                vec![SourceFormatArgument::Bytes(header)],
            )?;
        } else {
            let error_name = heap.allocate_model_storage(
                getInchiErrName(n_err)
                    .iter()
                    .map(|byte| *byte as i8)
                    .collect(),
            )?;
            let _ = eprint_call(
                heap,
                p_log.as_deref_mut(),
                "\n%s %s (%d) in %s (%d)\n",
                vec![
                    SourceFormatArgument::Bytes(header),
                    SourceFormatArgument::Bytes(error_name.as_const()),
                    SourceFormatArgument::Signed(i64::from(n_err)),
                    SourceFormatArgument::Bytes(state_message.as_const()),
                    SourceFormatArgument::Signed(i64::from(p_state)),
                ],
            )?;
            heap.free(error_name)?;
        }

        if ip.bINChIOutputOptions2 & INCHI_OUT_INCHI_GEN_ERROR as i32 != 0
            && ip.bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY as i32 == 0
        {
            let _ = eprint_call(
                heap,
                p_out.as_deref_mut(),
                "%s\n",
                vec![SourceFormatArgument::Bytes(header)],
            )?;
            let empty_inchi = if ip.bINChIOutputOptions & INCHI_OUT_STDINCHI as i32 != 0 {
                "InChI=1S//\n"
            } else {
                "InChI=1//\n"
            };
            let _ = eprint_call(heap, p_out.as_deref_mut(), empty_inchi, vec![])?;
        }
        if ip.bReadInChIOptions & READ_INCHI_TO_STRUCTURE as i32 != 0 {
            *num_errors = num_errors
                .checked_add(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        }
        *num_processed = num_processed
            .checked_add(1)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        heap.free(state_message)?;
        heap.free(simulation)?;
    }
    if !pstr_hdr.is_null() {
        inchi_free(heap, *pstr_hdr)?;
        *pstr_hdr = SourceMutPointer::null();
    }
    if !psz_cur_hdr.is_null() {
        inchi_free(heap, *psz_cur_hdr)?;
        *psz_cur_hdr = SourceMutPointer::null();
    }
    FreeInpInChI(heap, p_one_input)?;
    Ok(())
}

fn eprint_call(
    heap: &mut SourceHeap,
    log: Option<&mut INCHI_IOSTREAM>,
    format: &str,
    arguments: Vec<SourceFormatArgument>,
) -> Result<i32, SourceHeapError> {
    let format = heap.allocate_model_storage(
        format
            .bytes()
            .chain(std::iter::once(0))
            .map(|byte| byte as i8)
            .collect(),
    )?;
    let result = inchi_ios_eprint(
        heap,
        log,
        format.as_const(),
        &SourceVaList {
            arguments,
            position: 0,
        },
    );
    heap.free(format)?;
    result
}

#[allow(non_snake_case)]
pub(crate) fn PrepareSaveOptBits(
    heap: &mut SourceHeap,
    ip: &mut INPUT_PARMS,
    mut log: Option<&mut INCHI_IOSTREAM>,
    num_inp: i64,
    sz_cur_hdr: SourceConstPointer<i8>,
    input_has_save_opt: i32,
    input_save_opt_bits: u8,
    save_opt_bits: &mut u8,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10771 PrepareSaveOptBits
    // INCHI✔❌: void PrepareSaveOptBits(INPUT_PARMS* ip,
    // INCHI✔❌:     INCHI_IOSTREAM* pLog,
    // INCHI✔❌:     const long     num_inp,
    // INCHI✔❌:     const char* szCurHdr,
    // INCHI✔❌:     int            input_has_save_opt,
    // INCHI✔❌:     unsigned char  input_save_opt_bits,
    // INCHI✔❌:     unsigned char* save_opt_bits)
    // INCHI✔❌: {
    // INCHI✔❌:
    // INCHI✔❌:     if (!input_has_save_opt)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* Does not allow to create SaveOpt if the source lacks appendix */
    // INCHI✔❌:         ip->bINChIOutputOptions &= ~INCHI_OUT_SAVEOPT;
    // INCHI✔❌:         if (szCurHdr && szCurHdr[0])
    // INCHI✔❌:         {
    // INCHI✔❌:             inchi_ios_eprint(pLog,
    // INCHI✔❌:                 "Warning: ignore SaveOpt request for SaveOpt-less input, %s\n",
    // INCHI✔❌:                 szCurHdr);
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:         {
    // INCHI✔❌:             inchi_ios_eprint(pLog,
    // INCHI✔❌:                 "Warning: ignore SaveOpt request for SaveOpt-less input, Structure %ld\n",
    // INCHI✔❌:                 num_inp);
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:     else
    // INCHI✔❌:     {
    // INCHI✔❌:         /* Analyze existing and prepare new SaveOpt appendix */
    // INCHI✔❌:         /* djb-rwth: addressing coverity ID #499490 -- these are variable initialisers setting values to 0 */
    // INCHI✔❌:         int input_save_opt_has_recmet = input_save_opt_bits & SAVE_OPT_RECMET;
    // INCHI✔❌:         int input_save_opt_has_fixedh = input_save_opt_bits & SAVE_OPT_FIXEDH;
    // INCHI✔❌:         int input_save_opt_has_suu = input_save_opt_bits & SAVE_OPT_SUU;
    // INCHI✔❌:         int input_save_opt_has_sluud = input_save_opt_bits & SAVE_OPT_SLUUD;
    // INCHI✔❌:         int input_save_opt_has_ket = input_save_opt_bits & SAVE_OPT_KET;
    // INCHI✔❌:         int input_save_opt_has_15t = input_save_opt_bits & SAVE_OPT_15T;
    // INCHI✔❌:         int input_save_opt_has_pt_22_00 = input_save_opt_bits & SAVE_OPT_PT_22_00;
    // INCHI✔❌:         int input_save_opt_has_pt_16_00 = input_save_opt_bits & SAVE_OPT_PT_16_00;
    // INCHI✔❌:         int input_save_opt_has_pt_06_00 = input_save_opt_bits & SAVE_OPT_PT_06_00;
    // INCHI✔❌:         int input_save_opt_has_pt_39_00 = input_save_opt_bits & SAVE_OPT_PT_39_00;
    // INCHI✔❌:         int input_save_opt_has_pt_13_00 = input_save_opt_bits & SAVE_OPT_PT_13_00;
    // INCHI✔❌:         int input_save_opt_has_pt_18_00 = input_save_opt_bits & SAVE_OPT_PT_18_00;
    // INCHI✔❌:
    // INCHI✔❌:         if (0 != (ip->bTautFlags & TG_FLAG_RECONNECT_COORD))
    // INCHI✔❌:         {
    // INCHI✔❌:             /* RecMet requested */
    // INCHI✔❌:             if (input_save_opt_has_recmet)
    // INCHI✔❌:             {
    // INCHI✔❌:                 *save_opt_bits |= SAVE_OPT_RECMET;
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->bTautFlags &= ~TG_FLAG_RECONNECT_COORD;
    // INCHI✔❌:                 if (szCurHdr && szCurHdr[0])
    // INCHI✔❌:                 {
    // INCHI✔❌:                     inchi_ios_eprint(pLog, "Warning: input created w/o RecMet - ignoring RecMet request, %s\n", szCurHdr);
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                     inchi_ios_eprint(pLog, "Warning: input created w/o RecMet - ignoring RecMet request, Structure %ld\n", num_inp);
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         if (0 != (ip->nMode & REQ_MODE_BASIC))
    // INCHI✔❌:         {
    // INCHI✔❌:             /* FixedH requested */
    // INCHI✔❌:             if (input_save_opt_has_fixedh)
    // INCHI✔❌:             {
    // INCHI✔❌:                 *save_opt_bits |= SAVE_OPT_FIXEDH;
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 ip->nMode &= ~REQ_MODE_BASIC;
    // INCHI✔❌:                 if (szCurHdr && szCurHdr[0])
    // INCHI✔❌:                 {
    // INCHI✔❌:                     inchi_ios_eprint(pLog, "Warning: input created w/o FixedH - ignoring FixedH request, %s\n", szCurHdr);
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                     inchi_ios_eprint(pLog, "Warning: input created w/o FixedH - ignoring FixedH request, Structure %ld\n", num_inp);
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         /* Copy from source SaveOpt those bits which we do not touch    */
    // INCHI✔❌:         /* while converting InChI:         SUU SLUUD KET 15T            */
    // INCHI✔❌:         if (input_save_opt_has_suu)
    // INCHI✔❌:         {
    // INCHI✔❌:             *save_opt_bits |= SAVE_OPT_SUU;
    // INCHI✔❌:         }
    // INCHI✔❌:         if (input_save_opt_has_sluud)
    // INCHI✔❌:         {
    // INCHI✔❌:             *save_opt_bits |= SAVE_OPT_SLUUD;
    // INCHI✔❌:         }
    // INCHI✔❌:         if (input_save_opt_has_ket)
    // INCHI✔❌:         {
    // INCHI✔❌:             *save_opt_bits |= SAVE_OPT_KET;
    // INCHI✔❌:         }
    // INCHI✔❌:         if (input_save_opt_has_15t)
    // INCHI✔❌:         {
    // INCHI✔❌:             *save_opt_bits |= SAVE_OPT_15T;
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         if (input_save_opt_has_pt_22_00)
    // INCHI✔❌:             *save_opt_bits |= SAVE_OPT_PT_22_00;
    // INCHI✔❌:         if (input_save_opt_has_pt_16_00)
    // INCHI✔❌:             *save_opt_bits |= SAVE_OPT_PT_16_00;
    // INCHI✔❌:         if (input_save_opt_has_pt_06_00)
    // INCHI✔❌:             *save_opt_bits |= SAVE_OPT_PT_06_00;
    // INCHI✔❌:         if (input_save_opt_has_pt_39_00)
    // INCHI✔❌:             *save_opt_bits |= SAVE_OPT_PT_39_00;
    // INCHI✔❌:         if (input_save_opt_has_pt_13_00)
    // INCHI✔❌:             *save_opt_bits |= SAVE_OPT_PT_13_00;
    // INCHI✔❌:         if (input_save_opt_has_pt_18_00)
    // INCHI✔❌:             *save_opt_bits |= SAVE_OPT_PT_18_00;
    // INCHI✔❌:
    // INCHI✔❌:         /* Check if /SNon requested and turn OFF stereo bits if so */
    // INCHI✔❌:         if (!(ip->nMode & REQ_MODE_STEREO))
    // INCHI✔❌:         {
    // INCHI✔❌:             *save_opt_bits &= ~SAVE_OPT_SUU;
    // INCHI✔❌:             *save_opt_bits &= ~SAVE_OPT_SLUUD;
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: PrepareSaveOptBits

    let header_is_present = if sz_cur_hdr.is_null() {
        false
    } else {
        heap.slice(sz_cur_hdr)?
            .first()
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            != 0
    };
    let mut warn =
        |format_with_header: &str, format_with_number: &str| -> Result<(), SourceHeapError> {
            if header_is_present {
                let _ = eprint_call(
                    heap,
                    log.as_deref_mut(),
                    format_with_header,
                    vec![SourceFormatArgument::Bytes(sz_cur_hdr)],
                )?;
            } else {
                let _ = eprint_call(
                    heap,
                    log.as_deref_mut(),
                    format_with_number,
                    vec![SourceFormatArgument::Signed(num_inp)],
                )?;
            }
            Ok(())
        };

    if input_has_save_opt == 0 {
        ip.bINChIOutputOptions &= !(INCHI_OUT_SAVEOPT as i32);
        warn(
            "Warning: ignore SaveOpt request for SaveOpt-less input, %s\n",
            "Warning: ignore SaveOpt request for SaveOpt-less input, Structure %ld\n",
        )?;
        return Ok(());
    }

    let has = |mask: u32| u32::from(input_save_opt_bits) & mask != 0;
    if ip.bTautFlags & TG_FLAG_RECONNECT_COORD as INCHI_MODE != 0 {
        if has(SAVE_OPT_RECMET) {
            *save_opt_bits |= SAVE_OPT_RECMET as u8;
        } else {
            ip.bTautFlags &= !(TG_FLAG_RECONNECT_COORD as INCHI_MODE);
            warn(
                "Warning: input created w/o RecMet - ignoring RecMet request, %s\n",
                "Warning: input created w/o RecMet - ignoring RecMet request, Structure %ld\n",
            )?;
        }
    }
    if ip.nMode & REQ_MODE_BASIC as INCHI_MODE != 0 {
        if has(SAVE_OPT_FIXEDH) {
            *save_opt_bits |= SAVE_OPT_FIXEDH as u8;
        } else {
            ip.nMode &= !(REQ_MODE_BASIC as INCHI_MODE);
            warn(
                "Warning: input created w/o FixedH - ignoring FixedH request, %s\n",
                "Warning: input created w/o FixedH - ignoring FixedH request, Structure %ld\n",
            )?;
        }
    }

    for mask in [
        SAVE_OPT_SUU,
        SAVE_OPT_SLUUD,
        SAVE_OPT_KET,
        SAVE_OPT_15T,
        SAVE_OPT_PT_22_00,
        SAVE_OPT_PT_16_00,
        SAVE_OPT_PT_06_00,
        SAVE_OPT_PT_39_00,
        SAVE_OPT_PT_13_00,
        SAVE_OPT_PT_18_00,
    ] {
        if has(mask) {
            *save_opt_bits |= mask as u8;
        }
    }
    if ip.nMode & REQ_MODE_STEREO as INCHI_MODE == 0 {
        *save_opt_bits &= !(SAVE_OPT_SUU as u8);
        *save_opt_bits &= !(SAVE_OPT_SLUUD as u8);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{
        COMPONENT_REM_PROTONS, INCHI_IOS_TYPE_STRING,
        local_ichiread::tagModeProtonIsoExchgH_MODE_PIXH_UNDEFINED as MODE_PIXH_UNDEFINED,
    };

    fn output_stream() -> INCHI_IOSTREAM {
        INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        }
    }

    #[test]
    fn source_port__ichiread__sethillformfrominchi__line_583() {
        fn text(heap: &mut SourceHeap, value: &[u8]) -> SourceMutPointer<i8> {
            heap.allocate_model_storage(
                value
                    .iter()
                    .copied()
                    .map(|byte| byte as i8)
                    .chain(std::iter::once(0))
                    .collect(),
            )
            .unwrap()
        }

        fn component(
            heap: &mut SourceHeap,
            atoms: &[u8],
            hydrogens: &[i8],
            formula: SourceMutPointer<i8>,
            deleted: i32,
        ) -> INChI {
            INChI {
                nNumberOfAtoms: atoms.len() as i32,
                szHillFormula: formula,
                nAtom: if atoms.is_empty() {
                    SourceMutPointer::null()
                } else {
                    heap.allocate_model_storage(atoms.to_vec()).unwrap()
                },
                nNum_H: if hydrogens.is_empty() {
                    SourceMutPointer::null()
                } else {
                    heap.allocate_model_storage(hydrogens.to_vec()).unwrap()
                },
                bDeleted: deleted,
                ..INChI::default()
            }
        }

        let mut heap = SourceHeap::default();
        let equal_old = text(&mut heap, b"C");
        let different_old = text(&mut heap, b"incorrect");
        let deleted_old = text(&mut heap, b"kept-deleted");
        let equal = component(&mut heap, &[6], &[0], equal_old, 0);
        let deleted = component(&mut heap, &[8], &[0], deleted_old, 1);
        let different = component(&mut heap, &[8], &[0], different_old, 0);
        let fixed_components = heap.allocate_model_storage(vec![equal, deleted]).unwrap();
        let mobile_components = heap.allocate_model_storage(vec![different]).unwrap();
        let mut input = InpInChI::default();
        input.pInpInChI[INCHI_BAS as usize][TAUT_NON as usize] = fixed_components;
        input.nNumComponents[INCHI_BAS as usize][TAUT_NON as usize] = 2;
        input.pInpInChI[INCHI_REC as usize][TAUT_YES as usize] = mobile_components;
        input.nNumComponents[INCHI_REC as usize][TAUT_YES as usize] = 1;

        assert_eq!(SetHillFormFromInChI(&mut heap, &mut input), Ok(1));
        let equal_new = heap.slice(fixed_components.as_const()).unwrap()[0].szHillFormula;
        let different_new = heap.slice(mobile_components.as_const()).unwrap()[0].szHillFormula;
        assert_ne!(equal_new, equal_old);
        assert_ne!(different_new, different_old);
        assert_eq!(heap.slice(equal_new.as_const()).unwrap(), &[b'C' as i8, 0]);
        assert_eq!(
            heap.slice(different_new.as_const()).unwrap(),
            &[b'O' as i8, 0]
        );
        assert_eq!(
            heap.slice(equal_old.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(different_old.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(fixed_components.as_const()).unwrap()[1].szHillFormula,
            deleted_old
        );
        assert!(heap.slice(deleted_old.as_const()).is_ok());

        let mut skip_heap = SourceHeap::default();
        let zero_old = text(&mut skip_heap, b"zero-atoms");
        let deleted_old = text(&mut skip_heap, b"deleted");
        let empty_old = text(&mut skip_heap, b"");
        let zero = component(&mut skip_heap, &[], &[], zero_old, 0);
        let deleted = component(&mut skip_heap, &[8], &[0], deleted_old, 1);
        let null_formula = component(&mut skip_heap, &[8], &[0], SourceMutPointer::null(), 0);
        let empty = component(&mut skip_heap, &[8], &[0], empty_old, 0);
        let skips = skip_heap
            .allocate_model_storage(vec![zero, deleted, null_formula, empty])
            .unwrap();
        let mut skip_input = InpInChI::default();
        skip_input.pInpInChI[0][0] = skips;
        skip_input.nNumComponents[0][0] = 4;
        skip_heap.trace_source_allocations();
        assert_eq!(SetHillFormFromInChI(&mut skip_heap, &mut skip_input), Ok(0));
        assert_eq!(skip_heap.source_allocation_calls(), 0);
        assert_eq!(
            skip_heap.slice(skips.as_const()).unwrap()[0].szHillFormula,
            zero_old
        );
        assert_eq!(
            skip_heap.slice(skips.as_const()).unwrap()[1].szHillFormula,
            deleted_old
        );
        assert!(
            skip_heap.slice(skips.as_const()).unwrap()[2]
                .szHillFormula
                .is_null()
        );
        assert_eq!(
            skip_heap.slice(skips.as_const()).unwrap()[3].szHillFormula,
            empty_old
        );

        let mut failure_heap = SourceHeap::default();
        let old = text(&mut failure_heap, b"O");
        let value = component(&mut failure_heap, &[8], &[0], old, 0);
        let components = failure_heap.allocate_model_storage(vec![value]).unwrap();
        let mut failure_input = InpInChI::default();
        failure_input.pInpInChI[0][0] = components;
        failure_input.nNumComponents[0][0] = 1;
        failure_heap.trace_source_allocations();
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            SetHillFormFromInChI(&mut failure_heap, &mut failure_input),
            Ok(1)
        );
        assert_eq!(failure_heap.source_allocation_calls(), 1);
        assert!(
            failure_heap.slice(components.as_const()).unwrap()[0]
                .szHillFormula
                .is_null()
        );
        assert_eq!(
            failure_heap.slice(old.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__ichiread__convertinchi2inchi__line_11003() {
        fn run(
            heap: &mut SourceHeap,
            input_parameters: &mut INPUT_PARMS,
            current_header: &mut SourceMutPointer<i8>,
            processing_time: &mut i64,
        ) -> Result<i32, SourceHeapError> {
            ConvertInChI2InChI(
                heap,
                input_parameters,
                &mut InpInChI::default(),
                &mut output_stream(),
                &mut output_stream(),
                &STRUCT_DATA::default(),
                &mut [0, 0],
                &[MODE_PIXH_UNDEFINED, MODE_PIXH_UNDEFINED],
                current_header,
                17,
                &mut 23,
                0xa5,
                &mut inchiTime { clockTime: -1 },
                processing_time,
                &mut INCHI_CLOCK::default(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                1_000,
                6_000,
            )
        }

        let mut heap = SourceHeap::default();
        let sdf_value = heap.allocate_model_storage(vec![11_i8, 0]).unwrap();
        let sdf_label = heap.allocate_model_storage(vec![12_i8, 0]).unwrap();
        let header = heap.allocate_model_storage(vec![13_i8, 0]).unwrap();
        let mut current_header = header;
        let mut processing_time = 7_i64;
        let mut input_parameters = INPUT_PARMS {
            bNoStructLabels: 91,
            pSdfValue: sdf_value,
            pSdfLabel: sdf_label,
            ..INPUT_PARMS::default()
        };
        assert_eq!(
            run(
                &mut heap,
                &mut input_parameters,
                &mut current_header,
                &mut processing_time,
            ),
            Ok(0)
        );
        assert_eq!(input_parameters.bNoStructLabels, 91);
        assert!(input_parameters.pSdfValue.is_null());
        assert!(input_parameters.pSdfLabel.is_null());
        assert!(current_header.is_null());
        assert_eq!(processing_time, 12);
        assert_eq!(
            heap.slice(header.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert!(heap.slice(sdf_value.as_const()).is_ok());
        assert!(heap.slice(sdf_label.as_const()).is_ok());

        let mut failure_heap = SourceHeap::default();
        let header = failure_heap
            .allocate_model_storage(vec![b'H' as i8, 0])
            .unwrap();
        let mut current_header = header;
        let mut processing_time = -8_i64;
        let mut input_parameters = INPUT_PARMS {
            bNoStructLabels: -37,
            ..INPUT_PARMS::default()
        };
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            run(
                &mut failure_heap,
                &mut input_parameters,
                &mut current_header,
                &mut processing_time,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(input_parameters.bNoStructLabels, -37);
        assert!(current_header.is_null());
        assert_eq!(processing_time, -3);
        assert_eq!(
            failure_heap.slice(header.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__ichiread__getnumneighborsfrominchi__line_1646() {
        let mut heap = SourceHeap::default();
        assert_eq!(GetNumNeighborsFromInchi(&heap, None, 0), Ok(0));

        let connection_table = heap
            .allocate_model_storage(vec![1_u16, 2, 1, 3, 2])
            .unwrap();
        let hydrogens = heap.allocate_model_storage(vec![2_i8, -1, 0]).unwrap();
        let mut inchi = INChI {
            nNumberOfAtoms: 3,
            lenConnTable: 5,
            nConnTable: connection_table,
            nNum_H: hydrogens,
            ..INChI::default()
        };
        assert_eq!(GetNumNeighborsFromInchi(&heap, Some(&inchi), 1), Ok(3));
        assert_eq!(GetNumNeighborsFromInchi(&heap, Some(&inchi), 2), Ok(1));
        assert_eq!(GetNumNeighborsFromInchi(&heap, Some(&inchi), 3), Ok(1));

        let tautomer = heap
            .allocate_model_storage(vec![1_u16, 4, 1, 0, 1, 3])
            .unwrap();
        inchi.lenTautomer = 6;
        inchi.nTautomer = tautomer;
        assert_eq!(GetNumNeighborsFromInchi(&heap, Some(&inchi), 1), Ok(1003));
        assert_eq!(GetNumNeighborsFromInchi(&heap, Some(&inchi), 2), Ok(1));
        assert_eq!(GetNumNeighborsFromInchi(&heap, Some(&inchi), 3), Ok(1001));

        let one_vertex = heap.allocate_model_storage(vec![1_u16]).unwrap();
        let no_hydrogens = INChI {
            nNumberOfAtoms: 1,
            lenConnTable: 1,
            nConnTable: one_vertex,
            ..INChI::default()
        };
        assert_eq!(
            GetNumNeighborsFromInchi(&heap, Some(&no_hydrogens), 1),
            Ok(0)
        );

        let invalid_vertex = heap.allocate_model_storage(vec![1_u16, 3]).unwrap();
        let invalid_connection = INChI {
            nNumberOfAtoms: 2,
            lenConnTable: 2,
            nConnTable: invalid_vertex,
            ..INChI::default()
        };
        assert_eq!(
            GetNumNeighborsFromInchi(&heap, Some(&invalid_connection), 1),
            Ok(RI_ERR_PROGR)
        );

        let malformed_tautomer = heap
            .allocate_model_storage(vec![1_u16, 3, 1, 0, 1, 0])
            .unwrap();
        let malformed = INChI {
            nNumberOfAtoms: 1,
            lenConnTable: 1,
            nConnTable: one_vertex,
            lenTautomer: 6,
            nTautomer: malformed_tautomer,
            ..INChI::default()
        };
        assert_eq!(
            GetNumNeighborsFromInchi(&heap, Some(&malformed), 1),
            Ok(RI_ERR_PROGR)
        );
    }

    #[test]
    fn source_port__ichiread__countstereotypes__line_1723() {
        fn count(
            heap: &SourceHeap,
            inchi: &INChI,
            counters: &mut [i32; 6],
        ) -> Result<i32, SourceHeapError> {
            let [known_sb, known_sc, unknown_sb, unknown_sc, piii, asiii] = counters;
            CountStereoTypes(
                heap, inchi, known_sb, known_sc, unknown_sb, unknown_sc, piii, asiii,
            )
        }

        let mut heap = SourceHeap::default();
        let mut counters = [7, 11, 13, 17, 19, 23];
        assert_eq!(count(&heap, &INChI::default(), &mut counters), Ok(0));
        assert_eq!(counters, [7, 11, 13, 17, 19, 23]);
        assert_eq!(
            count(
                &heap,
                &INChI {
                    nNumberOfAtoms: 1,
                    bDeleted: 1,
                    ..INChI::default()
                },
                &mut counters,
            ),
            Ok(0)
        );

        let no_stereo = INChI {
            nNumberOfAtoms: 1,
            ..INChI::default()
        };
        assert_eq!(count(&heap, &no_stereo, &mut counters), Ok(1));

        let zero_isotopic = heap
            .allocate_model_storage(vec![crate::source_types::INChI_Stereo::default()])
            .unwrap();
        let atom_numbers = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
        let center_parities = heap.allocate_model_storage(vec![1_i8, 3]).unwrap();
        let bond_parities = heap.allocate_model_storage(vec![2_i8, 0]).unwrap();
        let regular = heap
            .allocate_model_storage(vec![crate::source_types::INChI_Stereo {
                nNumberOfStereoCenters: 2,
                nNumber: atom_numbers,
                t_parity: center_parities,
                nNumberOfStereoBonds: 2,
                b_parity: bond_parities,
                ..crate::source_types::INChI_Stereo::default()
            }])
            .unwrap();
        let elements = heap.allocate_model_storage(vec![15_u8, 33]).unwrap();
        let hydrogens = heap.allocate_model_storage(vec![3_i8, 3]).unwrap();
        let connection = heap.allocate_model_storage(vec![1_u16]).unwrap();
        let inchi = INChI {
            nNumberOfAtoms: 2,
            nAtom: elements,
            lenConnTable: 1,
            nConnTable: connection,
            nNum_H: hydrogens,
            Stereo: regular,
            StereoIsotopic: zero_isotopic,
            ..INChI::default()
        };
        let mut counters = [0; 6];
        assert_eq!(count(&heap, &inchi, &mut counters), Ok(2));
        assert_eq!(counters, [1, 1, 1, 1, 1, 1]);

        let isotopic_numbers = heap.allocate_model_storage(vec![1_u16]).unwrap();
        let isotopic_parity = heap.allocate_model_storage(vec![2_i8]).unwrap();
        let isotopic = heap
            .allocate_model_storage(vec![crate::source_types::INChI_Stereo {
                nNumberOfStereoCenters: 1,
                nNumber: isotopic_numbers,
                t_parity: isotopic_parity,
                ..crate::source_types::INChI_Stereo::default()
            }])
            .unwrap();
        let carbon = heap.allocate_model_storage(vec![6_u8, 33]).unwrap();
        let isotopic_inchi = INChI {
            StereoIsotopic: isotopic,
            nAtom: carbon,
            ..inchi.clone()
        };
        let mut counters = [0; 6];
        assert_eq!(count(&heap, &isotopic_inchi, &mut counters), Ok(2));
        assert_eq!(counters, [0, 1, 0, 0, 0, 0]);

        let invalid_number = heap.allocate_model_storage(vec![0_u16]).unwrap();
        let invalid_bond_parity = heap.allocate_model_storage(vec![1_i8]).unwrap();
        let invalid = heap
            .allocate_model_storage(vec![crate::source_types::INChI_Stereo {
                nNumberOfStereoCenters: 1,
                nNumber: invalid_number,
                t_parity: isotopic_parity,
                nNumberOfStereoBonds: 1,
                b_parity: invalid_bond_parity,
                ..crate::source_types::INChI_Stereo::default()
            }])
            .unwrap();
        let invalid_inchi = INChI {
            Stereo: invalid,
            StereoIsotopic: SourceMutPointer::null(),
            ..inchi.clone()
        };
        let mut counters = [0; 6];
        assert_eq!(
            count(&heap, &invalid_inchi, &mut counters),
            Ok(RI_ERR_PROGR)
        );
        assert_eq!(counters, [1, 0, 0, 0, 0, 0]);

        let phosphorus_number = heap.allocate_model_storage(vec![1_u16]).unwrap();
        let phosphorus_stereo = heap
            .allocate_model_storage(vec![crate::source_types::INChI_Stereo {
                nNumberOfStereoCenters: 1,
                nNumber: phosphorus_number,
                t_parity: isotopic_parity,
                ..crate::source_types::INChI_Stereo::default()
            }])
            .unwrap();
        let bad_connection = heap.allocate_model_storage(vec![1_u16, 3]).unwrap();
        let helper_error = INChI {
            nNumberOfAtoms: 2,
            nAtom: elements,
            lenConnTable: 2,
            nConnTable: bad_connection,
            nNum_H: hydrogens,
            Stereo: phosphorus_stereo,
            StereoIsotopic: SourceMutPointer::null(),
            ..INChI::default()
        };
        let mut counters = [0; 6];
        assert_eq!(count(&heap, &helper_error, &mut counters), Ok(RI_ERR_PROGR));
        assert_eq!(counters, [0, 1, 0, 0, 0, 0]);
    }

    #[test]
    fn source_port__ichiread__detectinpinchicreationoptions__line_1880() {
        fn detect(
            heap: &SourceHeap,
            input: &InpInChI,
            output: &mut [i32; 5],
        ) -> Result<i32, SourceHeapError> {
            let [reconnected, metal, fixed_h, mode, taut] = output;
            DetectInpInchiCreationOptions(heap, input, reconnected, metal, fixed_h, mode, taut)
        }

        let heap = SourceHeap::default();
        let mut output = [1; 5];
        assert_eq!(detect(&heap, &InpInChI::default(), &mut output), Ok(0));
        assert_eq!(
            output,
            [
                0,
                0,
                0,
                (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU) as i32,
                0
            ]
        );

        for (source_mode, expected) in [
            (1, REQ_MODE_STEREO | REQ_MODE_ISO_STEREO),
            (
                2,
                REQ_MODE_STEREO | REQ_MODE_ISO_STEREO | REQ_MODE_RELATIVE_STEREO,
            ),
            (
                3,
                REQ_MODE_STEREO | REQ_MODE_ISO_STEREO | REQ_MODE_RACEMIC_STEREO,
            ),
        ] {
            let mut input = InpInChI::default();
            input.s[0][0][1] = source_mode;
            input.s[0][0][0] = if source_mode == 2 { 3 } else { 2 };
            let mut output = [0; 5];
            assert_eq!(detect(&heap, &input, &mut output), Ok(0));
            assert_eq!(
                output[3],
                (expected | REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU) as i32
            );
        }
        let mut nonisotopic_mode = InpInChI::default();
        nonisotopic_mode.s[0][0][0] = 1;
        let mut output = [0; 5];
        assert_eq!(detect(&heap, &nonisotopic_mode, &mut output), Ok(0));
        assert_eq!(
            output[3],
            (REQ_MODE_STEREO
                | REQ_MODE_ISO_STEREO
                | REQ_MODE_SB_IGN_ALL_UU
                | REQ_MODE_SC_IGN_ALL_UU) as i32
        );

        let mut heap = SourceHeap::default();
        let no_stereo = heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                ..INChI::default()
            }])
            .unwrap();
        let mut component_input = InpInChI::default();
        component_input.pInpInChI[INCHI_REC as usize][TAUT_NON as usize] = no_stereo;
        component_input.nNumComponents[INCHI_REC as usize][TAUT_NON as usize] = 1;
        let mut output = [0; 5];
        assert_eq!(detect(&heap, &component_input, &mut output), Ok(0));
        assert_eq!(output[0], 1);
        assert_eq!(output[2], 1);

        let atom_numbers = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
        let center_parities = heap.allocate_model_storage(vec![1_i8, 3]).unwrap();
        let bond_parities = heap.allocate_model_storage(vec![0_i8]).unwrap();
        let stereo = heap
            .allocate_model_storage(vec![crate::source_types::INChI_Stereo {
                nNumberOfStereoCenters: 2,
                nNumber: atom_numbers,
                t_parity: center_parities,
                nNumberOfStereoBonds: 1,
                b_parity: bond_parities,
                ..crate::source_types::INChI_Stereo::default()
            }])
            .unwrap();
        let elements = heap.allocate_model_storage(vec![15_u8, 33]).unwrap();
        let hydrogens = heap.allocate_model_storage(vec![3_i8, 3]).unwrap();
        let connection = heap.allocate_model_storage(vec![1_u16]).unwrap();
        let stereo_component = heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 2,
                nAtom: elements,
                lenConnTable: 1,
                nConnTable: connection,
                nNum_H: hydrogens,
                Stereo: stereo,
                ..INChI::default()
            }])
            .unwrap();
        let mut stereo_input = InpInChI::default();
        stereo_input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize] = stereo_component;
        stereo_input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 1;
        let mut output = [0; 5];
        assert_eq!(detect(&heap, &stereo_input, &mut output), Ok(0));
        assert_eq!(output[0], 0);
        assert_eq!(output[2], 0);
        assert_eq!(
            output[3],
            (REQ_MODE_STEREO | REQ_MODE_ISO_STEREO | REQ_MODE_SC_IGN_ALL_UU) as i32
        );
        assert_eq!(
            output[4],
            (TG_FLAG_PHOSPHINE_STEREO | TG_FLAG_ARSINE_STEREO) as i32
        );

        let sodium = heap.allocate_model_storage(vec![11_u8]).unwrap();
        let sodium_h = heap.allocate_model_storage(vec![1_i8]).unwrap();
        let metal_component = heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                nAtom: sodium,
                nNum_H: sodium_h,
                ..INChI::default()
            }])
            .unwrap();
        let mut metal_input = InpInChI::default();
        metal_input.pInpInChI[0][TAUT_YES as usize] = metal_component;
        metal_input.nNumComponents[0][TAUT_YES as usize] = 1;
        let mut output = [0; 5];
        assert_eq!(detect(&heap, &metal_input, &mut output), Ok(0));
        assert_eq!(output[1], 1);

        let invalid_number = heap.allocate_model_storage(vec![0_u16]).unwrap();
        let invalid_parity = heap.allocate_model_storage(vec![1_i8]).unwrap();
        let invalid_stereo = heap
            .allocate_model_storage(vec![crate::source_types::INChI_Stereo {
                nNumberOfStereoCenters: 1,
                nNumber: invalid_number,
                t_parity: invalid_parity,
                ..crate::source_types::INChI_Stereo::default()
            }])
            .unwrap();
        let invalid_component = heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                Stereo: invalid_stereo,
                ..INChI::default()
            }])
            .unwrap();
        let mut invalid_input = InpInChI::default();
        invalid_input.pInpInChI[0][0] = invalid_component;
        invalid_input.nNumComponents[0][0] = 1;
        let mut output = [9; 5];
        assert_eq!(detect(&heap, &invalid_input, &mut output), Ok(RI_ERR_PROGR));
        assert_eq!(output, [0; 5]);
    }

    #[test]
    fn source_port__ichiread__brevinchicomponentexists__line_1838() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            bRevInchiComponentExists(&heap, None, INCHI_BAS as i32, TAUT_NON as i32, 0),
            Ok(0)
        );
        let empty = StrFromINChI::default();
        for (representation, mobile_h, component) in [
            (INCHI_BAS as i32, TAUT_NON as i32, 0),
            (-1, TAUT_NON as i32, 0),
            (INCHI_NUM as i32, TAUT_NON as i32, 0),
            (INCHI_BAS as i32, -1, 0),
            (INCHI_BAS as i32, TAUT_NUM as i32, 0),
            (INCHI_BAS as i32, TAUT_NON as i32, -1),
        ] {
            assert_eq!(
                bRevInchiComponentExists(&heap, Some(&empty), representation, mobile_h, component,),
                Ok(0)
            );
        }

        let mut null_table = StrFromINChI {
            num_atoms: 1,
            ..StrFromINChI::default()
        };
        null_table.RevInChI.num_components[INCHI_BAS as usize] = 1;
        assert_eq!(
            bRevInchiComponentExists(
                &heap,
                Some(&null_table),
                INCHI_BAS as i32,
                TAUT_NON as i32,
                0,
            ),
            Ok(0)
        );

        let zero = heap.allocate_model_storage(vec![INChI::default()]).unwrap();
        let deleted = heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                bDeleted: 1,
                ..INChI::default()
            }])
            .unwrap();
        let live = heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                ..INChI::default()
            }])
            .unwrap();
        let components = heap
            .allocate_model_storage(vec![
                [zero, SourceMutPointer::null()],
                [SourceMutPointer::null(), deleted],
                [live, live],
            ])
            .unwrap();
        let mut structure = StrFromINChI {
            num_atoms: 3,
            ..StrFromINChI::default()
        };
        structure.RevInChI.pINChI[INCHI_BAS as usize] = components;
        structure.RevInChI.num_components[INCHI_BAS as usize] = 3;
        for (mobile_h, component, expected) in [
            (TAUT_NON as i32, 0, 0),
            (TAUT_YES as i32, 0, 0),
            (TAUT_YES as i32, 1, 0),
            (TAUT_NON as i32, 2, 1),
            (TAUT_YES as i32, 2, 1),
            (TAUT_NON as i32, 3, 0),
        ] {
            assert_eq!(
                bRevInchiComponentExists(
                    &heap,
                    Some(&structure),
                    INCHI_BAS as i32,
                    mobile_h,
                    component,
                ),
                Ok(expected)
            );
        }
    }

    fn output_requested(
        heap: &mut SourceHeap,
        ip: &INPUT_PARMS,
        sd: &STRUCT_DATA,
        input: &mut InpInChI,
        counts: &mut [i32; INCHI_NUM as usize],
        modes: &[MODE_PIXH; INCHI_NUM as usize],
    ) -> Result<i32, SourceHeapError> {
        OutputInChIAsRequested(
            heap,
            SourceMutPointer::null(),
            &mut output_stream(),
            &mut output_stream(),
            ip,
            sd,
            input,
            counts,
            modes,
            i64::MAX,
            u8::MAX,
            SourceMutPointer::null(),
        )
    }

    #[test]
    fn source_port__ichiread__outputinchiasrequested__line_1213() {
        let mut initialization_failure = SourceHeap::default();
        initialization_failure.fail_after_allocations(0);
        assert_eq!(
            output_requested(
                &mut initialization_failure,
                &INPUT_PARMS::default(),
                &STRUCT_DATA::default(),
                &mut InpInChI::default(),
                &mut [0, 0],
                &[MODE_PIXH_UNDEFINED, MODE_PIXH_UNDEFINED],
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(initialization_failure.source_allocation_calls(), 1);

        for failure_after in 1..=6 {
            let mut heap = SourceHeap::default();
            let source = heap.allocate_model_storage(vec![INChI::default()]).unwrap();
            let proton = heap
                .allocate_model_storage(vec![COMPONENT_REM_PROTONS {
                    nNumRemovedProtons: 7,
                    nNumRemovedIsotopicH: [11, 13, 17],
                }])
                .unwrap();
            let mut input = InpInChI::default();
            input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize] = source;
            input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 1;
            input.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize].pNumProtons = proton;
            let mut counts = [1, 0];
            heap.fail_after_allocations(failure_after);
            let result = output_requested(
                &mut heap,
                &INPUT_PARMS::default(),
                &STRUCT_DATA {
                    num_components: [1, 0],
                    ..STRUCT_DATA::default()
                },
                &mut input,
                &mut counts,
                &[MODE_PIXH_ADD_TO_EACH, MODE_PIXH_UNDEFINED],
            );
            if failure_after <= 4 {
                assert_eq!(result, Ok(RI_ERR_ALLOC), "failure position {failure_after}");
                assert_eq!(counts, [1, 0]);
                assert_eq!(input.pInpInChI[0][TAUT_YES as usize], source);
                assert!(heap.slice(source.as_const()).is_ok());
                assert!(heap.slice(proton.as_const()).is_ok());
            } else {
                assert_eq!(result, Ok(0), "failure position {failure_after}");
                assert_eq!(counts, [0, 0]);
                assert!(input.pInpInChI[0][TAUT_YES as usize].is_null());
                assert!(
                    input.nNumProtons[0][TAUT_YES as usize]
                        .pNumProtons
                        .is_null()
                );
                assert_eq!(
                    heap.slice(source.as_const()),
                    Err(SourceHeapError::MissingAllocation)
                );
                assert_eq!(
                    heap.slice(proton.as_const()),
                    Err(SourceHeapError::MissingAllocation)
                );
            }
        }

        for mode in [
            MODE_PIXH_UNDEFINED,
            MODE_PIXH_ADD_TO_FIRST,
            MODE_PIXH_ADD_TO_EACH,
            MODE_PIXH_ADD_A_PIXH_COMPONENT,
        ] {
            let mut heap = SourceHeap::default();
            let source = heap
                .allocate_model_storage(vec![INChI {
                    bDeleted: i32::from(mode == MODE_PIXH_UNDEFINED),
                    ..INChI::default()
                }])
                .unwrap();
            let per_component = heap
                .allocate_model_storage(vec![COMPONENT_REM_PROTONS {
                    nNumRemovedProtons: i16::MIN,
                    nNumRemovedIsotopicH: [i16::MIN, 0, i16::MAX],
                }])
                .unwrap();
            let mut input = InpInChI::default();
            input.pInpInChI[0][TAUT_YES as usize] = source;
            input.nNumComponents[0][TAUT_YES as usize] = 1;
            input.nNumProtons[0][TAUT_YES as usize].nNumRemovedProtons = i16::MAX;
            input.nNumProtons[0][TAUT_YES as usize].nNumRemovedIsotopicH = [i16::MAX, -1, i16::MIN];
            input.nNumProtons[0][TAUT_YES as usize].pNumProtons = per_component;
            let mut counts = [2, 0];
            assert_eq!(
                output_requested(
                    &mut heap,
                    &INPUT_PARMS::default(),
                    &STRUCT_DATA {
                        num_components: [2, 0],
                        ..STRUCT_DATA::default()
                    },
                    &mut input,
                    &mut counts,
                    &[mode, MODE_PIXH_UNDEFINED],
                ),
                Ok(0),
                "proton mode {mode}"
            );
            assert_eq!(counts, [0, 0]);
            assert!(input.pInpInChI[0][TAUT_YES as usize].is_null());
            assert!(
                input.nNumProtons[0][TAUT_YES as usize]
                    .pNumProtons
                    .is_null()
            );
            assert_eq!(
                heap.slice(source.as_const()),
                Err(SourceHeapError::MissingAllocation)
            );
            assert_eq!(
                heap.slice(per_component.as_const()),
                Err(SourceHeapError::MissingAllocation)
            );
        }

        let mut equal_heap = SourceHeap::default();
        let make_equal = |heap: &mut SourceHeap, tautomer_length: i32| {
            let formula = heap.allocate_model_storage(vec![b'C' as i8, 0]).unwrap();
            let atom = heap.allocate_model_storage(vec![6_u8]).unwrap();
            let hydrogen = heap.allocate_model_storage(vec![0_i8]).unwrap();
            (
                INChI {
                    nNumberOfAtoms: 1,
                    szHillFormula: formula,
                    nAtom: atom,
                    nNum_H: hydrogen,
                    lenTautomer: tautomer_length,
                    ..INChI::default()
                },
                (formula, atom, hydrogen),
            )
        };
        let (fixed, fixed_allocations) = make_equal(&mut equal_heap, 0);
        let (mobile, mobile_allocations) = make_equal(&mut equal_heap, 0);
        let fixed_source = equal_heap.allocate_model_storage(vec![fixed]).unwrap();
        let mobile_source = equal_heap.allocate_model_storage(vec![mobile]).unwrap();
        let mut equal_input = InpInChI::default();
        equal_input.pInpInChI[0] = [fixed_source, mobile_source];
        equal_input.nNumComponents[0] = [1, 1];
        let mut equal_counts = [1, 0];
        assert_eq!(
            output_requested(
                &mut equal_heap,
                &INPUT_PARMS {
                    nMode: REQ_MODE_BASIC as INCHI_MODE,
                    ..INPUT_PARMS::default()
                },
                &STRUCT_DATA {
                    num_components: [1, 0],
                    ..STRUCT_DATA::default()
                },
                &mut equal_input,
                &mut equal_counts,
                &[MODE_PIXH_UNDEFINED, MODE_PIXH_UNDEFINED],
            ),
            Ok(0)
        );
        assert_eq!(equal_counts, [0, 0]);
        for pointer in [fixed_source.as_const(), mobile_source.as_const()] {
            assert_eq!(
                equal_heap.slice(pointer),
                Err(SourceHeapError::MissingAllocation)
            );
        }
        for (formula, atom, hydrogen) in [fixed_allocations, mobile_allocations] {
            assert_eq!(
                equal_heap.slice(formula.as_const()),
                Err(SourceHeapError::MissingAllocation)
            );
            assert_eq!(
                equal_heap.slice(atom.as_const()),
                Err(SourceHeapError::MissingAllocation)
            );
            assert_eq!(
                equal_heap.slice(hydrogen.as_const()),
                Err(SourceHeapError::MissingAllocation)
            );
        }

        let mut split_heap = SourceHeap::default();
        let sodium_atom = split_heap.allocate_model_storage(vec![11_u8]).unwrap();
        let sodium_h = split_heap.allocate_model_storage(vec![1_i8]).unwrap();
        let sodium_formula = split_heap
            .allocate_model_storage(vec![b'H' as i8, b'N' as i8, b'a' as i8, 0])
            .unwrap();
        let reconnected_source = split_heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                nAtom: sodium_atom,
                nNum_H: sodium_h,
                szHillFormula: sodium_formula,
                ..INChI::default()
            }])
            .unwrap();
        let mut split_input = InpInChI::default();
        split_input.pInpInChI[INCHI_REC as usize][TAUT_YES as usize] = reconnected_source;
        split_input.nNumComponents[INCHI_REC as usize][TAUT_YES as usize] = 1;
        let mut split_counts = [0, 1];
        assert_eq!(
            output_requested(
                &mut split_heap,
                &INPUT_PARMS {
                    bReadInChIOptions: READ_INCHI_SPLIT_OUTPUT as i32,
                    bTautFlags: (TG_FLAG_DISCONNECT_COORD | TG_FLAG_RECONNECT_COORD) as INCHI_MODE,
                    ..INPUT_PARMS::default()
                },
                &STRUCT_DATA {
                    num_components: [0, 1],
                    ..STRUCT_DATA::default()
                },
                &mut split_input,
                &mut split_counts,
                &[MODE_PIXH_UNDEFINED, MODE_PIXH_UNDEFINED],
            ),
            Ok(0)
        );
        assert_eq!(split_counts, [0, 0]);
        assert!(split_input.pInpInChI[INCHI_REC as usize][TAUT_YES as usize].is_null());
        assert_eq!(
            split_heap.slice(reconnected_source.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            split_heap.slice(sodium_atom.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            split_heap.slice(sodium_h.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            split_heap.slice(sodium_formula.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__ichiread__binchihasreconnectedmetal__line_1993() {
        let mut heap = SourceHeap::default();
        assert_eq!(bInChIHasReconnectedMetal(&heap, None), Ok(0));
        for inchi in [
            INChI {
                bDeleted: 1,
                nNumberOfAtoms: 1,
                ..INChI::default()
            },
            INChI::default(),
            INChI {
                nNumberOfAtoms: -1,
                ..INChI::default()
            },
        ] {
            assert_eq!(bInChIHasReconnectedMetal(&heap, Some(&inchi)), Ok(0));
        }

        let carbon = heap.allocate(vec![6_u8]).unwrap();
        let sodium = heap.allocate(vec![11_u8]).unwrap();
        let sodium_carbon = heap.allocate(vec![11_u8, 6]).unwrap();
        let zero_h = heap.allocate(vec![0_i8]).unwrap();
        let one_h = heap.allocate(vec![1_i8]).unwrap();
        let carbon_inchi = INChI {
            nNumberOfAtoms: 1,
            nAtom: carbon,
            ..INChI::default()
        };
        assert_eq!(bInChIHasReconnectedMetal(&heap, Some(&carbon_inchi)), Ok(0));

        let mut metal = INChI {
            nNumberOfAtoms: 1,
            nAtom: sodium,
            ..INChI::default()
        };
        assert_eq!(bInChIHasReconnectedMetal(&heap, Some(&metal)), Ok(0));
        metal.nNum_H = zero_h;
        assert_eq!(bInChIHasReconnectedMetal(&heap, Some(&metal)), Ok(0));
        metal.nNum_H = one_h;
        assert_eq!(bInChIHasReconnectedMetal(&heap, Some(&metal)), Ok(1));
        metal.nNumberOfAtoms = 2;
        metal.nAtom = sodium_carbon;
        metal.nNum_H = SourceMutPointer::null();
        assert_eq!(bInChIHasReconnectedMetal(&heap, Some(&metal)), Ok(1));
    }
    use crate::source_types::{INCHI_IOS_STRING, inp_ATOM};

    fn c_string(heap: &mut SourceHeap, text: &str) -> SourceConstPointer<i8> {
        heap.allocate_model_storage(
            text.bytes()
                .chain(std::iter::once(0))
                .map(|byte| byte as i8)
                .collect(),
        )
        .unwrap()
        .as_const()
    }

    fn log_text(heap: &SourceHeap, log: &INCHI_IOSTREAM) -> String {
        let bytes = heap.slice(log.s.pStr.as_const()).unwrap();
        String::from_utf8(
            bytes[..log.s.nUsedLength as usize]
                .iter()
                .map(|byte| *byte as u8)
                .collect(),
        )
        .unwrap()
    }

    #[test]
    fn source_port__ichiread__preparesaveoptbits__line_10771() {
        let mut heap = SourceHeap::default();
        let mut log = INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        let header = c_string(&mut heap, "Header");
        let empty_header = c_string(&mut heap, "");

        let mut ip = INPUT_PARMS {
            bINChIOutputOptions: INCHI_OUT_SAVEOPT as i32 | 0x40,
            nMode: REQ_MODE_BASIC as INCHI_MODE,
            bTautFlags: TG_FLAG_RECONNECT_COORD as INCHI_MODE,
            ..INPUT_PARMS::default()
        };
        let mut save = 0xa5;
        PrepareSaveOptBits(
            &mut heap,
            &mut ip,
            Some(&mut log),
            -7,
            header,
            0,
            0xff,
            &mut save,
        )
        .unwrap();
        assert_eq!(ip.bINChIOutputOptions, 0x40);
        assert_eq!(ip.nMode, REQ_MODE_BASIC as INCHI_MODE);
        assert_eq!(ip.bTautFlags, TG_FLAG_RECONNECT_COORD as INCHI_MODE);
        assert_eq!(save, 0xa5);
        assert_eq!(
            log_text(&heap, &log),
            "Warning: ignore SaveOpt request for SaveOpt-less input, Header\n"
        );

        log.s = INCHI_IOS_STRING::default();
        ip = INPUT_PARMS {
            nMode: (REQ_MODE_BASIC | REQ_MODE_STEREO) as INCHI_MODE,
            bTautFlags: TG_FLAG_RECONNECT_COORD as INCHI_MODE,
            ..INPUT_PARMS::default()
        };
        save = 0;
        PrepareSaveOptBits(
            &mut heap,
            &mut ip,
            Some(&mut log),
            i64::MAX,
            SourceConstPointer::null(),
            1,
            0xff,
            &mut save,
        )
        .unwrap();
        assert_eq!(ip.nMode, (REQ_MODE_BASIC | REQ_MODE_STEREO) as INCHI_MODE);
        assert_eq!(ip.bTautFlags, TG_FLAG_RECONNECT_COORD as INCHI_MODE);
        assert_eq!(save, 0xff);
        assert_eq!(log.s.nUsedLength, 0);

        ip = INPUT_PARMS {
            nMode: (REQ_MODE_BASIC | REQ_MODE_STEREO) as INCHI_MODE,
            bTautFlags: TG_FLAG_RECONNECT_COORD as INCHI_MODE,
            ..INPUT_PARMS::default()
        };
        save = 0x08;
        PrepareSaveOptBits(
            &mut heap,
            &mut ip,
            Some(&mut log),
            -9,
            empty_header,
            -1,
            0,
            &mut save,
        )
        .unwrap();
        assert_eq!(ip.nMode, REQ_MODE_STEREO as INCHI_MODE);
        assert_eq!(ip.bTautFlags, 0);
        assert_eq!(save, 0x08);
        assert_eq!(
            log_text(&heap, &log),
            concat!(
                "Warning: input created w/o RecMet - ignoring RecMet request, Structure -9\n",
                "Warning: input created w/o FixedH - ignoring FixedH request, Structure -9\n"
            )
        );

        ip = INPUT_PARMS::default();
        save = 0x08;
        PrepareSaveOptBits(
            &mut heap,
            &mut ip,
            None,
            0,
            SourceConstPointer::null(),
            1,
            0xff,
            &mut save,
        )
        .unwrap();
        assert_eq!(save, 0xf8);

        assert_eq!(SAVE_OPT_PT_06_00 as u8, 0);
        assert_eq!(SAVE_OPT_PT_39_00 as u8, 0);
        assert_eq!(SAVE_OPT_PT_13_00 as u8, 0);
        assert_eq!(SAVE_OPT_PT_18_00 as u8, 0);
    }

    #[test]
    fn source_port__ichiread__getinchistatereaderr__line_536() {
        let expected = [
            "MOBILE_H_FORMULA",
            "MOBILE_H_CONNECTIONS",
            "MOBILE_H",
            "MOBILE_H_CHARGE",
            "MOBILE_H_PROTONS",
            "MOBILE_H_SP2",
            "MOBILE_H_SP3",
            "MOBILE_H_SP3_/m",
            "MOBILE_H_SP3_/s",
            "MOBILE_H_ISO_LAYER_FORK",
            "MOBILE_H_ISO_ATOMS",
            "MOBILE_H_ISO_EXCH_H",
            "MOBILE_H_ISO_SP2",
            "MOBILE_H_ISO_SP3",
            "MOBILE_H_ISO_SP3_/m",
            "MOBILE_H_ISO_SP3_/s",
            "FIXED_H_LAYER_FORK",
            "FIXED_H_FORMULA",
            "FIXED_H",
            "FIXED_H_CHARGE",
            "FIXED_H_SP2",
            "FIXED_H_SP3",
            "FIXED_H_SP3_/m",
            "FIXED_H_SP3_/s",
            "FIXED_H_PERMUTATION",
            "FIXED_H_ISO_LAYER_FORK",
            "FIXED_H_ISO_ATOMS",
            "FIXED_H_ISO_LAYER",
            "FIXED_H_ISO_SP2",
            "FIXED_H_ISO_SP3",
            "FIXED_H_ISO_SP3_m",
            "FIXED_H_ISO_SP3_s",
            "FIXED_H_ISO_PERMUTATION",
            "RECONNECTED_LAYER_FORK",
            "RECONNECTED_FORMULA",
            "MATERIAL_BALANCE",
            "POLYMER_LAYER",
        ];
        let mut heap = SourceHeap::default();
        let buffer = heap.allocate_model_storage(vec![99_i8; 128]).unwrap();
        for (stat, expected) in expected.iter().enumerate() {
            for (offset, suffix) in [(0, ""), (100, ", Reconnected layer")] {
                heap.slice_mut(buffer).unwrap().fill(99);
                getInchiStateReadErr(&mut heap, stat as i32 + offset, buffer).unwrap();
                let output = heap.slice(buffer.as_const()).unwrap();
                let expected = format!("{expected}{suffix}");
                assert_eq!(
                    &output[..expected.len()],
                    &expected.bytes().map(|byte| byte as i8).collect::<Vec<_>>()
                );
                assert_eq!(output[expected.len()], 0);
                assert!(output[expected.len() + 1..].iter().all(|byte| *byte == 99));
            }
        }
        for (stat, expected) in [
            (-1, "Unknown Error"),
            (37, "Unknown Error"),
            (99, "Unknown Error"),
            (137, "Unknown Error, Reconnected layer"),
            (i32::MAX, "Unknown Error, Reconnected layer"),
        ] {
            heap.slice_mut(buffer).unwrap().fill(99);
            getInchiStateReadErr(&mut heap, stat, buffer).unwrap();
            let output = heap.slice(buffer.as_const()).unwrap();
            assert_eq!(
                &output[..=expected.len()],
                &expected
                    .bytes()
                    .chain(std::iter::once(0))
                    .map(|byte| byte as i8)
                    .collect::<Vec<_>>()
            );
        }

        let short = heap.allocate_model_storage(vec![99_i8; 3]).unwrap();
        assert_eq!(
            getInchiStateReadErr(&mut heap, 0, short),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(heap.slice(short.as_const()).unwrap(), &[99; 3]);
        assert_eq!(
            getInchiStateReadErr(&mut heap, 0, SourceMutPointer::null()),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichiread__getinchierrname__line_562() {
        assert_eq!(getInchiErrName(RI_ERR_ALLOC), b"Allocation failed\0");
        assert_eq!(getInchiErrName(RI_ERR_PROGR), b"Program error\0");
        assert_eq!(getInchiErrName(RI_ERR_SYNTAX), b"Syntax error\0");
        assert_eq!(getInchiErrName(RI_ERR_EOL), b"End of line\0");
        for n_err in [i32::MIN, -9, -5, 0, 1, i32::MAX] {
            assert_eq!(getInchiErrName(n_err), b"Unknown error\0");
        }
    }

    #[test]
    fn source_port__ichiread__treaterrorsinreadinchistring__line_10901() {
        let stream = || INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };

        let mut heap = SourceHeap::default();
        let mut out = stream();
        let mut log = stream();
        let header = c_string(&mut heap, "EOF header").as_mut();
        let cur_header = c_string(&mut heap, "EOF current").as_mut();
        let atom = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let mut input = InpInChI {
            atom,
            num_atoms: 1,
            ..InpInChI::default()
        };
        let mut header_slot = header;
        let mut cur_header_slot = cur_header;
        let (mut num_inp, mut num_errors, mut num_processed) = (9, 8, 7);
        TreatErrorsInReadInChIString(
            &mut heap,
            RI_ERR_EOF as i32,
            0,
            0,
            &INPUT_PARMS::default(),
            Some(&mut out),
            Some(&mut log),
            &mut num_inp,
            &mut num_errors,
            &mut num_processed,
            &mut header_slot,
            &mut cur_header_slot,
            &mut input,
        )
        .unwrap();
        assert_eq!((num_inp, num_errors, num_processed), (9, 8, 7));
        assert_eq!(out.s.nUsedLength, 0);
        assert_eq!(log.s.nUsedLength, 0);
        assert!(header_slot.is_null());
        assert!(cur_header_slot.is_null());
        assert_eq!(input, InpInChI::default());
        for freed in [
            heap.slice(header.as_const()).map(|_| ()),
            heap.slice(cur_header.as_const()).map(|_| ()),
            heap.slice(atom.as_const()).map(|_| ()),
        ] {
            assert_eq!(freed, Err(SourceHeapError::MissingAllocation));
        }

        let mut heap = SourceHeap::default();
        let mut out = stream();
        let mut log = stream();
        let header = c_string(&mut heap, "Record 7").as_mut();
        let cur_header = c_string(&mut heap, "Current 7").as_mut();
        let mut header_slot = header;
        let mut cur_header_slot = cur_header;
        let mut input = InpInChI::default();
        let ip = INPUT_PARMS {
            nInputType: tagInputType_INPUT_INCHI,
            bReadInChIOptions: READ_INCHI_TO_STRUCTURE as i32,
            bINChIOutputOptions: INCHI_OUT_STDINCHI as i32,
            bINChIOutputOptions2: INCHI_OUT_INCHI_GEN_ERROR as i32,
            ..INPUT_PARMS::default()
        };
        let (mut num_inp, mut num_errors, mut num_processed) = (6, 2, 3);
        TreatErrorsInReadInChIString(
            &mut heap,
            RI_ERR_SYNTAX,
            RI_ERR_SYNTAX,
            5,
            &ip,
            Some(&mut out),
            Some(&mut log),
            &mut num_inp,
            &mut num_errors,
            &mut num_processed,
            &mut header_slot,
            &mut cur_header_slot,
            &mut input,
        )
        .unwrap();
        assert_eq!((num_inp, num_errors, num_processed), (7, 3, 4));
        assert_eq!(
            log_text(&heap, &log),
            "\nRecord 7 Syntax error (-2) in MOBILE_H_SP2 (5)\n"
        );
        assert_eq!(log_text(&heap, &out), "Record 7\nInChI=1S//\n");
        assert!(header_slot.is_null());
        assert!(cur_header_slot.is_null());
        assert_eq!(input, InpInChI::default());

        let mut heap = SourceHeap::default();
        let mut out = stream();
        let mut log = stream();
        let mut header_slot = SourceMutPointer::null();
        let mut cur_header_slot = SourceMutPointer::null();
        let mut input = InpInChI::default();
        let ip = INPUT_PARMS {
            bINChIOutputOptions2: INCHI_OUT_INCHI_GEN_ERROR as i32,
            ..INPUT_PARMS::default()
        };
        let (mut num_inp, mut num_errors, mut num_processed) = (0, 0, 0);
        TreatErrorsInReadInChIString(
            &mut heap,
            RI_ERR_SYNTAX,
            RI_ERR_SYNTAX,
            36,
            &ip,
            Some(&mut out),
            Some(&mut log),
            &mut num_inp,
            &mut num_errors,
            &mut num_processed,
            &mut header_slot,
            &mut cur_header_slot,
            &mut input,
        )
        .unwrap();
        assert_eq!((num_inp, num_errors, num_processed), (1, 0, 1));
        assert_eq!(
            log_text(&heap, &log),
            "Structure: 1 Skipping polymer InChI for conversion of InChI to InChI\n"
        );
        assert_eq!(log_text(&heap, &out), "Structure: 1\nInChI=1//\n");
        assert_eq!(input, InpInChI::default());
    }
}
