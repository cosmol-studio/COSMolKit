use std::cmp::Ordering;

use crate::source::base::ichi_io::inchi_ios_print;
use crate::source::base::ichierr::{AddErrorMessage, ErrMsg};
use crate::source::base::ichimake::{CompINChINonTaut2, CompINChITaut2};
use crate::source::base::ichiprt1::OutputINChI2;
use crate::source::base::mol_fmt2::MolfileSaveCopy;
use crate::source::base::runichi2::{eprint_bytes, sdf_label_value, source_array_c_bytes};
use crate::source::base::strutil::{Free_INChI, Free_INChI_Aux};
use crate::source::base::util::inchi_free;
use crate::source_types::{
    _IS_ERROR, _IS_FATAL, _IS_WARNING, AMBIGUOUS_STEREO_ATOM, AMBIGUOUS_STEREO_ATOM_ISO, AMBIGUOUS_STEREO_BOND,
    AMBIGUOUS_STEREO_BOND_ISO, CANON_GLOBALS, COMP_ATOM_DATA, CT_OUT_OF_RAM, CT_USER_QUIT_ERR, FILE, INCHI_BAS,
    INCHI_CLOCK, INCHI_IOS_STRING, INCHI_IOSTREAM, INCHI_MODE, INCHI_NUM, INCHI_OUT_EMBED_REC, INCHI_OUT_PLAIN_TEXT,
    INCHI_OUT_PLAIN_TEXT_COMMENTS, INCHI_OUT_PRINT_OPTIONS, INCHI_OUT_TABBED_OUTPUT, INCHI_SORT, INChI, INChI_Stereo,
    INP_ATOM_DATA, INPUT_PARMS, NORM_CANON_FLAGS, ORIG_ATOM_DATA, ORIG_STRUCT, OUT_TN, PINChI_Aux2, PINChI2,
    STRUCT_DATA, SourceConstPointer, SourceFormatArgument, SourceHeap, SourceHeapError, SourceMutPointer, SourceVaList,
    TAUT_NON, TAUT_NUM, TAUT_YES, TG_FLAG_DISCONNECT_COORD_DONE,
};

fn sort_inchi_rows(heap: &mut SourceHeap, rows: &mut [INCHI_SORT], tautomeric: bool) -> Result<(), SourceHeapError> {
    let mut error = None;
    rows.sort_unstable_by(|left, right| {
        if error.is_some() {
            return Ordering::Equal;
        }
        let result = if tautomeric {
            CompINChITaut2(heap, left, right)
        } else {
            CompINChINonTaut2(heap, left, right)
        };
        match result {
            Ok(value) => value.cmp(&0),
            Err(source_error) => {
                error = Some(source_error);
                Ordering::Equal
            }
        }
    });
    error.map_or(Ok(()), Err)
}

#[allow(non_snake_case)]
pub(crate) fn bIsStructChiral(
    heap: &SourceHeap,
    inchi_components: [SourceMutPointer<PINChI2>; INCHI_NUM as usize],
    component_counts: [i32; INCHI_NUM as usize],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1271 bIsStructChiral
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: bIsStructChiral
    // INCHI✔️❌: int bIsStructChiral( PINChI2 *pINChI2[INCHI_NUM], int num_components[] )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, k;
    // INCHI✔️❌:     INChI *pINChI;
    // INCHI✔️❌:     INChI_Stereo *Stereo;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (j = 0; j < INCHI_NUM; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* disconnected / reconnected */
    // INCHI✔️❌:         if (!num_components[j])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = 0; i < num_components[j]; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* i-th component */
    // INCHI✔️❌:             for (k = 0; k < TAUT_NUM; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* mobile/immobile H */
    // INCHI✔️❌:                 if (( pINChI = pINChI2[j][i][k] ) &&
    // INCHI✔️❌:                       !pINChI->bDeleted                &&
    // INCHI✔️❌:                       pINChI->nNumberOfAtoms > 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (( Stereo = pINChI->Stereo ) &&
    // INCHI✔️❌:                          Stereo->t_parity &&
    // INCHI✔️❌:                          Stereo->nNumberOfStereoCenters > 0 &&
    // INCHI✔️❌:                          Stereo->nCompInv2Abs)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return 1; /* inversion changed stereo */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (( Stereo = pINChI->StereoIsotopic ) &&
    // INCHI✔️❌:                          Stereo->t_parity &&
    // INCHI✔️❌:                          Stereo->nNumberOfStereoCenters > 0 &&
    // INCHI✔️❌:                          Stereo->nCompInv2Abs)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return 1; /* inversion changed stereo */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: bIsStructChiral
    // END INCHI C FUNCTION: bIsStructChiral
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: bIsStructChiral
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; INCHI_NUM=2 and TAUT_NUM=2.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: bIsStructChiral

    for domain in 0..INCHI_NUM as usize {
        let count = component_counts[domain];
        if count == 0 {
            continue;
        }
        for component in 0..count.max(0) {
            let component = usize::try_from(component).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let row = heap
                .slice(inchi_components[domain].as_const())?
                .get(component)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            for tautomer in 0..TAUT_NUM as usize {
                let inchi_pointer = row[tautomer];
                if inchi_pointer.is_null() {
                    continue;
                }
                let inchi = heap
                    .slice(inchi_pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if inchi.bDeleted != 0 || inchi.nNumberOfAtoms <= 0 {
                    continue;
                }
                for stereo_pointer in [inchi.Stereo, inchi.StereoIsotopic] {
                    if stereo_pointer.is_null() {
                        continue;
                    }
                    let stereo = heap
                        .slice(stereo_pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if !stereo.t_parity.is_null() && stereo.nNumberOfStereoCenters > 0 && stereo.nCompInv2Abs != 0 {
                        return Ok(1);
                    }
                }
            }
        }
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn SortAndPrintINChI(
    heap: &mut SourceHeap,
    canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    out_file: &mut INCHI_IOSTREAM,
    mut string_buffer: Option<&mut INCHI_IOS_STRING>,
    log_file: &mut INCHI_IOSTREAM,
    input_parameters: &INPUT_PARMS,
    original_atom_data: SourceConstPointer<ORIG_ATOM_DATA>,
    _prepared_atom_data: SourceConstPointer<ORIG_ATOM_DATA>,
    _composite_normalized_data: Option<&mut [[COMP_ATOM_DATA; TAUT_NUM as usize + 1]; INCHI_NUM as usize]>,
    original_structure: SourceConstPointer<ORIG_STRUCT>,
    component_counts: &[i32; INCHI_NUM as usize],
    non_tautomeric_counts: &[i32; INCHI_NUM as usize],
    tautomeric_counts: &[i32; INCHI_NUM as usize],
    taut_flags: &mut [INCHI_MODE; INCHI_NUM as usize],
    taut_flags_done: &mut [INCHI_MODE; INCHI_NUM as usize],
    normalization_flags: &NORM_CANON_FLAGS,
    input_structure_number: i64,
    inchi_components: [SourceMutPointer<PINChI2>; INCHI_NUM as usize],
    aux_components: [SourceMutPointer<PINChI_Aux2>; INCHI_NUM as usize],
    sort_print_flags: SourceMutPointer<i32>,
    save_option_bits: u8,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:102 SortAndPrintINChI
    // INCHI✔️❌: int SortAndPrintINChI( CANON_GLOBALS            *pCG,
    // INCHI✔️❌:                        INCHI_IOSTREAM           *out_file,
    // INCHI✔️❌:                        INCHI_IOS_STRING         *strbuf,
    // INCHI✔️❌:                        INCHI_IOSTREAM           *log_file,
    // INCHI✔️❌:                        INPUT_PARMS              *ip,
    // INCHI✔️❌:                        ORIG_ATOM_DATA           *orig_inp_data,
    // INCHI✔️❌:                        ORIG_ATOM_DATA           *prep_inp_data,
    // INCHI✔️❌:                        COMP_ATOM_DATA           composite_norm_data[INCHI_NUM][TAUT_NUM + 1],
    // INCHI✔️❌:                        ORIG_STRUCT              *pOrigStruct,
    // INCHI✔️❌:                        int                      num_components[INCHI_NUM],
    // INCHI✔️❌:                        int                      num_non_taut[INCHI_NUM],
    // INCHI✔️❌:                        int                      num_taut[INCHI_NUM],
    // INCHI✔️❌:                        INCHI_MODE               bTautFlags[INCHI_NUM],
    // INCHI✔️❌:                        INCHI_MODE               bTautFlagsDone[INCHI_NUM],
    // INCHI✔️❌:                        NORM_CANON_FLAGS         *pncFlags,
    // INCHI✔️❌:                        long                     num_inp,
    // INCHI✔️❌:                        PINChI2                  *pINChI[INCHI_NUM],
    // INCHI✔️❌:                        PINChI_Aux2              *pINChI_Aux[INCHI_NUM],
    // INCHI✔️❌:                        int                      *pSortPrintINChIFlags,
    // INCHI✔️❌:                        unsigned char            save_opt_bits )
    // INCHI✔️❌: {
    // INCHI✔️❌:     INCHI_SORT *pINChISort[INCHI_NUM][TAUT_NUM];
    // INCHI✔️❌:     int j, i, k, k1, ret, ret2, iINChI, max_num_components; /* djb-rwth: ignoring LLVM warning: variable used */
    // INCHI✔️❌:     int INCHI_basic_or_INCHI_reconnected;
    // INCHI✔️❌:     /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     int bDisconnectedCoord = ( 0 != ( bTautFlagsDone[0] & TG_FLAG_DISCONNECT_COORD_DONE ) );
    // INCHI✔️❌:     int bINChIOutputOptions0, bCurOption, bINChIOutputOptionsCur, bEmbedReconnected;
    // INCHI✔️❌:     static const char szAnnHdr[] = "InChI ANNOTATED CONTENTS";
    // INCHI✔️❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔️❌:     int ikflag = 0;
    // INCHI✔️❌:     out_file->type = INCHI_IOS_TYPE_STRING;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         Note:
    // INCHI✔️❌:
    // INCHI✔️❌:         pINChI[INCHI_BAS] refers to either disconnected or original structure;
    // INCHI✔️❌:             num_components[INCHI_BAS] > 0 if there was input structure
    // INCHI✔️❌:         pINChI[INCHI_REC] refers to the reconnected structure,
    // INCHI✔️❌:             and only if the input structure has been disconnected, that is,
    // INCHI✔️❌:             num_components[INCHI_REC] > 0
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = 1;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < INCHI_NUM; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (k = 0; k < TAUT_NUM; k++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             bTautFlags[i] |= pncFlags->bTautFlags[i][k];
    // INCHI✔️❌:             bTautFlagsDone[i] |= pncFlags->bTautFlagsDone[i][k];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:
    // INCHI✔️❌:     max_num_components = 0;
    // INCHI✔️❌:     for (j = 0; j < INCHI_NUM; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (max_num_components < num_components[j])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             max_num_components = num_components[j];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (max_num_components <= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         max_num_components = 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (j = 0, i = 0; j < INCHI_NUM; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (num_components[j])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (k1 = 0; k1 < TAUT_NUM; k1++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pINChISort[j][k1] =
    // INCHI✔️❌:                     (INCHI_SORT *) inchi_calloc( max_num_components,
    // INCHI✔️❌:                                                sizeof( pINChISort[0][0][0] ) );
    // INCHI✔️❌:                 i += !pINChISort[j][k1]; /* number of failed allocatons */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (k1 = 0; k1 < TAUT_NUM; k1++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pINChISort[j][k1] = NULL; /* keep BC happy */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = CT_OUT_OF_RAM;
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (j = 0; j < INCHI_NUM; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (!num_components[j])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         iINChI = j;
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( OUTPUT_CONNECTED_METAL_ONLY == 1 ) /* test: output connected as the only one INChI */
    // INCHI✔️❌:         if (INCHI_BAS == j && num_components[INCHI_REC])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             j = INCHI_REC;
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:         /* j = INCHI_BAS; <- for debug only */
    // INCHI✔️❌:         /* for only normal or disconnected coord compounds */
    // INCHI✔️❌:         /* (j=0=INCHI_BAS => normal or disconnected, j=1=INCHI_REC => reconnected */
    // INCHI✔️❌:
    // INCHI✔️❌:         for (k1 = 0; k1 < TAUT_NUM; k1++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (i = 0; i < num_components[j]; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (k = 0; k < TAUT_NUM; k++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     pINChISort[j][k1][i].pINChI[k] = pINChI[j][i][k];
    // INCHI✔️❌:                     pINChISort[j][k1][i].pINChI_Aux[k] = pINChI_Aux[j][i][k];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 pINChISort[j][k1][i].ord_number = i;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Sort component INChIs */
    // INCHI✔️❌:         for (k1 = 0; k1 < TAUT_NUM; k1++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             switch (k1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 case TAUT_NON:
    // INCHI✔️❌:                     qsort( pINChISort[j][k1],
    // INCHI✔️❌:                            num_components[j],
    // INCHI✔️❌:                            sizeof( pINChISort[0][0][0] ),
    // INCHI✔️❌:                            CompINChINonTaut2 );
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case TAUT_YES:
    // INCHI✔️❌:                     qsort( pINChISort[j][k1],
    // INCHI✔️❌:                            num_components[j],
    // INCHI✔️❌:                            sizeof( pINChISort[0][0][0] ),
    // INCHI✔️❌:                            CompINChITaut2 );
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌: #ifndef COMPILE_ANSI_ONLY
    // INCHI✔️❌:         /* Find equivalent and wINChI display order;
    // INCHI✔️❌:         use requested in ip->bCompareComponents comparison */
    // INCHI✔️❌:
    // INCHI✔️❌:         ret = SaveEquComponentsInfoAndSortOrder( iINChI, pINChISort[j],
    // INCHI✔️❌:                                                   num_components,
    // INCHI✔️❌:                                                   orig_inp_data, prep_inp_data,
    // INCHI✔️❌: #if ( FIX_DALKE_BUGS == 1 )
    // INCHI✔️❌:                                                   composite_norm_data ? composite_norm_data[j] : NULL,
    // INCHI✔️❌: #else
    // INCHI✔️❌:                                                   composite_norm_data[j],
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                                                   ip->bCompareComponents );
    // INCHI✔️❌:         if (RETURNED_ERROR( ret ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = 0;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     } /* j */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!( ip->bINChIOutputOptions & INCHI_OUT_PRINT_OPTIONS ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Prepare InChI from the structures obtained by
    // INCHI✔️❌:            reversing InChI for returning to the caller */
    // INCHI✔️❌:         for (j = 0; j < INCHI_NUM; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!num_components[j])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* pINChI[iINCHI][iComponent][bTaut] */
    // INCHI✔️❌:             /* j  = disconnected/connected */
    // INCHI✔️❌:             /* k1 = sort order for Mobile or Fixed H */
    // INCHI✔️❌:
    // INCHI✔️❌:             k1 = TAUT_YES; /* in Mobile H order */
    // INCHI✔️❌:             /* store components in Mobile H order */
    // INCHI✔️❌:
    // INCHI✔️❌:             for (i = 0; i < num_components[j] && i < max_num_components; i++) /* djb-rwth: fixing undefined index value / buffer overflow */
    // INCHI✔️❌:             {
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (pINChISort[j][k1][i].pINChI[TAUT_NON] &&
    // INCHI✔️❌:                     !pINChISort[j][k1][i].pINChI[TAUT_YES])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* make sure Mobile-H is always present */
    // INCHI✔️❌:                     for (k = 0; k < TAUT_NUM; k++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         pINChI[j][i][k] = pINChISort[j][k1][i].pINChI[ALT_TAUT( k )];
    // INCHI✔️❌:                         pINChI_Aux[j][i][k] = pINChISort[j][k1][i].pINChI_Aux[ALT_TAUT( k )];
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     for (k = 0; k < TAUT_NUM; k++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         pINChI[j][i][k] = pINChISort[j][k1][i].pINChI[k];
    // INCHI✔️❌:                         pINChI_Aux[j][i][k] = pINChISort[j][k1][i].pINChI_Aux[k];
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Print InChI string(s) */
    // INCHI✔️❌:         bINChIOutputOptions0 = ip->bINChIOutputOptions & ~INCHI_OUT_PRINT_OPTIONS;
    // INCHI✔️❌:         bEmbedReconnected = ip->bINChIOutputOptions & INCHI_OUT_EMBED_REC;
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = 1; i <= 2; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             bCurOption = INCHI_OUT_PLAIN_TEXT;
    // INCHI✔️❌:             if (i == 2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 bCurOption = INCHI_OUT_PLAIN_TEXT_COMMENTS;
    // INCHI✔️❌:                 /* continue; */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (!( ip->bINChIOutputOptions & bCurOption ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             bINChIOutputOptionsCur = bINChIOutputOptions0 | bCurOption;
    // INCHI✔️❌:             if (i == 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* output INChI */
    // INCHI✔️❌:                 bINChIOutputOptionsCur |= bEmbedReconnected;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else if (i == 2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* output annotation */
    // INCHI✔️❌:                 inchi_ios_print( out_file, "\n==== %s ====\n", szAnnHdr );
    // INCHI✔️❌:                 ITRACE_( "\n==== %s ====\n", szAnnHdr );
    // INCHI✔️❌:                 bINChIOutputOptionsCur |= bEmbedReconnected;
    // INCHI✔️❌:                 bINChIOutputOptionsCur &= ~INCHI_OUT_TABBED_OUTPUT;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌: #if 0 /*^^^#ifdef TARGET_LIB_FOR_WINCHI*/
    // INCHI✔️❌:             inchi_ios_free_str( out_file );
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:             INCHI_basic_or_INCHI_reconnected = INCHI_BAS;
    // INCHI✔️❌:
    // INCHI✔️❌:             ret2 = OutputINChI2( pCG,
    // INCHI✔️❌:                                  strbuf,
    // INCHI✔️❌:                                  pINChISort,
    // INCHI✔️❌:                                  INCHI_basic_or_INCHI_reconnected,
    // INCHI✔️❌:                                  orig_inp_data,
    // INCHI✔️❌:                                  pOrigStruct,
    // INCHI✔️❌:                                  ip,
    // INCHI✔️❌:                                  bDisconnectedCoord,
    // INCHI✔️❌:                                  OUT_TN,
    // INCHI✔️❌:                                  bINChIOutputOptionsCur,
    // INCHI✔️❌:                                  num_components,
    // INCHI✔️❌:                                  num_non_taut,
    // INCHI✔️❌:                                  num_taut,
    // INCHI✔️❌:                                  out_file,
    // INCHI✔️❌:                                  log_file,
    // INCHI✔️❌:                                  num_inp,
    // INCHI✔️❌:                                  pSortPrintINChIFlags,
    // INCHI✔️❌:                                  save_opt_bits );
    // INCHI✔️❌:
    // INCHI✔️❌:             ret &= ret2;
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔️❌:             /* always calculate InChIKey */
    // INCHI✔️❌:             winchi_calc_inchikey( ret, &ikflag, ip, out_file, log_file );
    // INCHI✔️❌:             /*push_to_winchi_text_window( out_file );
    // INCHI✔️❌:             inchi_ios_flush( out_file );*/
    // INCHI✔️❌:
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:             if (!ret)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         } /* i */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:
    // INCHI✔️❌:     for (j = 0; j < INCHI_NUM; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (k1 = 0; k1 < TAUT_NUM; k1++) /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (pINChISort[j][k1])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_free( pINChISort[j][k1] );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = ret ? 0 : _IS_FATAL;
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: SortAndPrintINChI

    let disconnected_coordinates = i32::from(taut_flags_done[0] & u64::from(TG_FLAG_DISCONNECT_COORD_DONE) != 0);
    for representation in 0..INCHI_NUM as usize {
        for tautomer in 0..TAUT_NUM as usize {
            taut_flags[representation] |= normalization_flags.bTautFlags[representation][tautomer];
            taut_flags_done[representation] |= normalization_flags.bTautFlagsDone[representation][tautomer];
        }
    }

    let mut max_components = component_counts.iter().copied().max().unwrap_or(0);
    if max_components <= 0 {
        max_components = 1;
    }
    let allocation_len = usize::try_from(max_components).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let mut sort_rows = [SourceMutPointer::<INCHI_SORT>::null(); 4];
    let mut allocation_failures = 0_i32;
    for representation in 0..INCHI_NUM as usize {
        if component_counts[representation] != 0 {
            for tautomer in 0..TAUT_NUM as usize {
                let index = representation * TAUT_NUM as usize + tautomer;
                match heap.allocate(vec![INCHI_SORT::default(); allocation_len]) {
                    Ok(pointer) => sort_rows[index] = pointer,
                    Err(SourceHeapError::AllocationFailed) => {
                        allocation_failures = allocation_failures.wrapping_add(1);
                    }
                    Err(error) => {
                        for pointer in sort_rows {
                            if !pointer.is_null() {
                                let _ = heap.free(pointer);
                            }
                        }
                        return Err(error);
                    }
                }
            }
        }
    }

    let mut sort_table = SourceMutPointer::<SourceMutPointer<INCHI_SORT>>::null();
    let execution = if allocation_failures != 0 {
        Ok(CT_OUT_OF_RAM)
    } else {
        (|| -> Result<i32, SourceHeapError> {
            let mut ret = 1_i32;
            for representation in 0..INCHI_NUM as usize {
                let component_count = component_counts[representation];
                if component_count == 0 {
                    continue;
                }
                let count = usize::try_from(component_count).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let source_inchi = heap
                    .slice(inchi_components[representation].as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                let source_aux = heap
                    .slice(aux_components[representation].as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();

                for tautomer_order in 0..TAUT_NUM as usize {
                    let row_pointer = sort_rows[representation * TAUT_NUM as usize + tautomer_order];
                    let row = heap
                        .slice_mut(row_pointer)?
                        .get_mut(..count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    for component in 0..count {
                        for tautomer in 0..TAUT_NUM as usize {
                            row[component].pINChI[tautomer] = source_inchi[component][tautomer];
                            row[component].pINChI_Aux[tautomer] = source_aux[component][tautomer];
                        }
                        row[component].ord_number = component as i16;
                    }
                }

                for tautomer_order in 0..TAUT_NUM as usize {
                    let row_pointer = sort_rows[representation * TAUT_NUM as usize + tautomer_order];
                    let mut row = heap
                        .slice(row_pointer.as_const())?
                        .get(..count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .to_vec();
                    sort_inchi_rows(heap, &mut row, tautomer_order == TAUT_YES as usize)?;
                    heap.slice_mut(row_pointer)?
                        .get_mut(..count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone_from_slice(&row);
                }
            }

            if input_parameters.bINChIOutputOptions & INCHI_OUT_PRINT_OPTIONS as i32 == 0 {
                for representation in 0..INCHI_NUM as usize {
                    let component_count = component_counts[representation];
                    if component_count == 0 {
                        continue;
                    }
                    let count = usize::try_from(component_count.min(max_components))
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    let mobile_pointer = sort_rows[representation * TAUT_NUM as usize + TAUT_YES as usize];
                    let mobile_rows = heap
                        .slice(mobile_pointer.as_const())?
                        .get(..count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .to_vec();
                    let mut reordered_inchi = heap
                        .slice(inchi_components[representation].as_const())?
                        .get(..count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .to_vec();
                    let mut reordered_aux = heap
                        .slice(aux_components[representation].as_const())?
                        .get(..count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .to_vec();
                    for component in 0..count {
                        let swap_tautomer_slots = !mobile_rows[component].pINChI[TAUT_NON as usize].is_null()
                            && mobile_rows[component].pINChI[TAUT_YES as usize].is_null();
                        for tautomer in 0..TAUT_NUM as usize {
                            let source_tautomer = if swap_tautomer_slots {
                                TAUT_NUM as usize - 1 - tautomer
                            } else {
                                tautomer
                            };
                            reordered_inchi[component][tautomer] = mobile_rows[component].pINChI[source_tautomer];
                            reordered_aux[component][tautomer] = mobile_rows[component].pINChI_Aux[source_tautomer];
                        }
                    }
                    heap.slice_mut(inchi_components[representation])?
                        .get_mut(..count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone_from_slice(&reordered_inchi);
                    heap.slice_mut(aux_components[representation])?
                        .get_mut(..count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone_from_slice(&reordered_aux);
                }
            } else {
                let base_options = input_parameters.bINChIOutputOptions & !(INCHI_OUT_PRINT_OPTIONS as i32);
                let embed_reconnected = input_parameters.bINChIOutputOptions & INCHI_OUT_EMBED_REC as i32;
                sort_table = heap.allocate_model_storage(sort_rows.to_vec())?;
                for pass in 1..=2 {
                    let current_option = if pass == 1 {
                        INCHI_OUT_PLAIN_TEXT as i32
                    } else {
                        INCHI_OUT_PLAIN_TEXT_COMMENTS as i32
                    };
                    if input_parameters.bINChIOutputOptions & current_option == 0 {
                        continue;
                    }
                    let mut current_options = base_options | current_option | embed_reconnected;
                    if pass == 2 {
                        let format = heap
                            .allocate_model_storage(b"\n==== %s ====\n\0".iter().map(|byte| *byte as i8).collect())?;
                        let header = heap.allocate_model_storage(
                            b"InChI ANNOTATED CONTENTS\0".iter().map(|byte| *byte as i8).collect(),
                        )?;
                        let print_result = inchi_ios_print(
                            heap,
                            Some(out_file),
                            stdout,
                            format.as_const(),
                            &SourceVaList {
                                arguments: vec![SourceFormatArgument::Bytes(header.as_const())],
                                ..SourceVaList::default()
                            },
                        );
                        let format_free = heap.free(format);
                        let header_free = heap.free(header);
                        print_result?;
                        format_free?;
                        header_free?;
                        current_options &= !(INCHI_OUT_TABBED_OUTPUT as i32);
                    }
                    let output_result = OutputINChI2(
                        heap,
                        canonical_globals,
                        string_buffer.as_deref_mut(),
                        sort_table,
                        INCHI_BAS as i32,
                        original_atom_data,
                        original_structure,
                        input_parameters,
                        disconnected_coordinates,
                        OUT_TN as i32,
                        current_options,
                        component_counts,
                        non_tautomeric_counts,
                        tautomeric_counts,
                        out_file,
                        log_file,
                        input_structure_number as i32,
                        sort_print_flags,
                        save_option_bits,
                        stdout,
                    )?;
                    ret &= output_result;
                    if ret == 0 {
                        break;
                    }
                }
            }
            Ok(ret)
        })()
    };

    let mut cleanup_error = None;
    if !sort_table.is_null() {
        if let Err(error) = heap.free(sort_table) {
            cleanup_error = Some(error);
        }
    }
    for pointer in sort_rows {
        if !pointer.is_null() {
            if let Err(error) = heap.free(pointer) {
                cleanup_error.get_or_insert(error);
            }
        }
    }
    let ret = execution?;
    if let Some(error) = cleanup_error {
        return Err(error);
    }
    Ok(if ret != 0 { 0 } else { _IS_FATAL as i32 })
}

#[allow(non_snake_case)]
pub(crate) fn FreeAllINChIArrays(
    heap: &mut SourceHeap,
    p_inchi: &mut [SourceMutPointer<PINChI2>; INCHI_NUM as usize],
    p_inchi_aux: &mut [SourceMutPointer<PINChI_Aux2>; INCHI_NUM as usize],
    num_components: &mut [i32; INCHI_NUM as usize],
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1322 FreeAllINChIArrays
    // INCHI✔️❌: void FreeAllINChIArrays( PINChI2 *pINChI[INCHI_NUM],
    // INCHI✔️❌:                          PINChI_Aux2 *pINChI_Aux[INCHI_NUM],
    // INCHI✔️❌:                          int num_components[INCHI_NUM] )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int k;
    // INCHI✔️❌:     for (k = 0; k < INCHI_NUM; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int nk = num_components[k];
    // INCHI✔️❌:
    // INCHI✔️❌:         FreeINChIArrays( pINChI[k], pINChI_Aux[k], num_components[k] );
    // INCHI✔️❌:
    // INCHI✔️❌:         num_components[k] = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nk &&                /* added check for nk: 2013-12-15 IPl */
    // INCHI✔️❌:              pINChI[k])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( pINChI[k] );
    // INCHI✔️❌:             pINChI[k] = NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nk &&                /* added check for nk: 2013-12-15 IPl */
    // INCHI✔️❌:              pINChI_Aux[k])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( pINChI_Aux[k] );
    // INCHI✔️❌:             pINChI_Aux[k] = NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: FreeAllINChIArrays

    for representation in 0..INCHI_NUM as usize {
        let old_count = num_components[representation];
        FreeINChIArrays(heap, p_inchi[representation], p_inchi_aux[representation], old_count)?;
        num_components[representation] = 0;
        if old_count != 0 && !p_inchi[representation].is_null() {
            inchi_free(heap, p_inchi[representation])?;
            p_inchi[representation] = SourceMutPointer::null();
        }
        if old_count != 0 && !p_inchi_aux[representation].is_null() {
            inchi_free(heap, p_inchi_aux[representation])?;
            p_inchi_aux[representation] = SourceMutPointer::null();
        }
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn FreeINChIArrays(
    heap: &mut SourceHeap,
    p_inchi: SourceMutPointer<PINChI2>,
    p_inchi_aux: SourceMutPointer<PINChI_Aux2>,
    num_components: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1357 FreeINChIArrays
    // INCHI✔️❌: void FreeINChIArrays( PINChI2 *pINChI, PINChI_Aux2 *pINChI_Aux, int num_components )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, k;
    // INCHI✔️❌:     /* Release allocated memory */
    // INCHI✔️❌:     if (pINChI)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = 0; i < num_components; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (k = 0; k < TAUT_NUM; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 Free_INChI( &pINChI[i][k] );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (pINChI_Aux)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (i = 0; i < num_components; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (k = 0; k < TAUT_NUM; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 Free_INChI_Aux( &pINChI_Aux[i][k] );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: FreeINChIArrays

    let count = if num_components > 0 {
        usize::try_from(num_components).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
    } else {
        0
    };
    if !p_inchi.is_null() {
        let mut rows = heap
            .slice(p_inchi.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec();
        for row in &mut rows {
            for owner in row.iter_mut().take(TAUT_NUM as usize) {
                Free_INChI(heap, owner)?;
            }
        }
        heap.slice_mut(p_inchi)?
            .get_mut(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone_from_slice(&rows);
    }
    if !p_inchi_aux.is_null() {
        let mut rows = heap
            .slice(p_inchi_aux.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec();
        for row in &mut rows {
            for owner in row.iter_mut().take(TAUT_NUM as usize) {
                Free_INChI_Aux(heap, owner)?;
            }
        }
        heap.slice_mut(p_inchi_aux)?
            .get_mut(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone_from_slice(&rows);
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn GetProcessingWarningsOneInChI(
    heap: &SourceHeap,
    p_inchi: &INChI,
    input_normalized_data: &INP_ATOM_DATA,
    error_buffer: &mut [i8],
    no_warnings: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1574 GetProcessingWarningsOneInChI
    // INCHI✔️✔️: int GetProcessingWarningsOneInChI( INChI *pINChI,
    // INCHI✔️✔️:                                    INP_ATOM_DATA *inp_norm_data,
    // INCHI✔️✔️:                                    char *pStrErrStruct,
    // INCHI✔️✔️:                                    int bNoWarnings )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int j;
    // INCHI✔️✔️:     int nAmbiguousStereoAtoms, nAmbiguousStereoBonds;
    // INCHI✔️✔️:     nAmbiguousStereoAtoms = 0;
    // INCHI✔️✔️:     nAmbiguousStereoBonds = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (inp_norm_data->at)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         for (j = 0; j < pINChI->nNumberOfAtoms; j++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (inp_norm_data->at[j].bAmbiguousStereo & ( AMBIGUOUS_STEREO_ATOM | AMBIGUOUS_STEREO_ATOM_ISO ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 nAmbiguousStereoAtoms++;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             if (inp_norm_data->at[j].bAmbiguousStereo & ( AMBIGUOUS_STEREO_BOND | AMBIGUOUS_STEREO_BOND_ISO ))
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 nAmbiguousStereoBonds++;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (nAmbiguousStereoAtoms)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (!bNoWarnings)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 WarningMessage( pStrErrStruct, "Ambiguous stereo:" );
    // INCHI✔️✔️:                 WarningMessage( pStrErrStruct, "center(s)" );
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (nAmbiguousStereoBonds)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (!bNoWarnings)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 WarningMessage( pStrErrStruct, "Ambiguous stereo:" );
    // INCHI✔️✔️:                 WarningMessage( pStrErrStruct, "bond(s)" );
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return ( nAmbiguousStereoAtoms || nAmbiguousStereoBonds );
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: GetProcessingWarningsOneInChI

    let mut ambiguous_stereo_atoms = 0_i32;
    let mut ambiguous_stereo_bonds = 0_i32;
    if !input_normalized_data.at.is_null() {
        let count = if p_inchi.nNumberOfAtoms > 0 {
            usize::try_from(p_inchi.nNumberOfAtoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
        } else {
            0
        };
        let atoms = heap
            .slice(input_normalized_data.at.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for atom in atoms {
            let ambiguity = i32::from(atom.bAmbiguousStereo);
            if ambiguity & (AMBIGUOUS_STEREO_ATOM as i32 | AMBIGUOUS_STEREO_ATOM_ISO as i32) != 0 {
                ambiguous_stereo_atoms += 1;
            }
            if ambiguity & (AMBIGUOUS_STEREO_BOND as i32 | AMBIGUOUS_STEREO_BOND_ISO as i32) != 0 {
                ambiguous_stereo_bonds += 1;
            }
        }
        if ambiguous_stereo_atoms != 0 && no_warnings == 0 {
            let heading = b"Ambiguous stereo:\0".map(|byte| byte as i8);
            let centers = b"center(s)\0".map(|byte| byte as i8);
            AddErrorMessage(Some(error_buffer), Some(&heading))?;
            AddErrorMessage(Some(error_buffer), Some(&centers))?;
        }
        if ambiguous_stereo_bonds != 0 && no_warnings == 0 {
            let heading = b"Ambiguous stereo:\0".map(|byte| byte as i8);
            let bonds = b"bond(s)\0".map(|byte| byte as i8);
            AddErrorMessage(Some(error_buffer), Some(&heading))?;
            AddErrorMessage(Some(error_buffer), Some(&bonds))?;
        }
    }
    Ok(i32::from(ambiguous_stereo_atoms != 0 || ambiguous_stereo_bonds != 0))
}

#[allow(non_snake_case)]
pub(crate) fn GetProcessingWarningsOneComponentInChI(
    heap: &SourceHeap,
    current_inchi: &[SourceMutPointer<INChI>; TAUT_NUM as usize],
    input_normalized_data: &[INP_ATOM_DATA; TAUT_NUM as usize],
    structure_data: &mut STRUCT_DATA,
    no_warnings: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1550 GetProcessingWarningsOneComponentInChI
    // INCHI✔️✔️: int GetProcessingWarningsOneComponentInChI( INChI *cur_INChI[],
    // INCHI✔️✔️:                                             INP_ATOM_DATA **inp_norm_data,
    // INCHI✔️✔️:                                             STRUCT_DATA *sd,
    // INCHI✔️✔️:                                             int bNoWarnings )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i, ret = 0;
    // INCHI✔️✔️:     for (i = 0; i < TAUT_NUM; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (cur_INChI[i] && cur_INChI[i]->nNumberOfAtoms > 0)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             ret |= GetProcessingWarningsOneInChI( cur_INChI[i],
    // INCHI✔️✔️:                                                   inp_norm_data[i],
    // INCHI✔️✔️:                                                   sd->pStrErrStruct,
    // INCHI✔️✔️:                                                   bNoWarnings );
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return ret;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: GetProcessingWarningsOneComponentInChI

    let mut result = 0_i32;
    for representation in 0..TAUT_NUM as usize {
        if !current_inchi[representation].is_null() {
            let inchi = heap
                .slice(current_inchi[representation].as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if inchi.nNumberOfAtoms > 0 {
                result |= GetProcessingWarningsOneInChI(
                    heap,
                    inchi,
                    &input_normalized_data[representation],
                    &mut structure_data.pStrErrStruct,
                    no_warnings,
                )?;
            }
        }
    }
    Ok(result)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn TreatErrorsInCreateOneComponentINChI(
    heap: &mut SourceHeap,
    structure_data: &mut STRUCT_DATA,
    input_parameters: &INPUT_PARMS,
    _original_atom_data: &ORIG_ATOM_DATA,
    component_index: i32,
    input_structure_number: i64,
    _input_file: Option<&mut INCHI_IOSTREAM>,
    log_file: Option<&mut INCHI_IOSTREAM>,
    _output_file: Option<&mut INCHI_IOSTREAM>,
    _problem_file: Option<&mut INCHI_IOSTREAM>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1387 TreatErrorsInCreateOneComponentINChI
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: TreatErrorsInCreateOneComponentINChI
    // INCHI✔️❌: int TreatErrorsInCreateOneComponentINChI( STRUCT_DATA *sd,
    // INCHI✔️❌:                                           INPUT_PARMS    *ip,
    // INCHI✔️❌:                                           ORIG_ATOM_DATA *orig_inp_data,
    // INCHI✔️❌:                                           int i, long num_inp,
    // INCHI✔️❌:                                           INCHI_IOSTREAM *inp_file,
    // INCHI✔️❌:                                           INCHI_IOSTREAM *log_file,
    // INCHI✔️❌:                                           INCHI_IOSTREAM *out_file,
    // INCHI✔️❌:                                           INCHI_IOSTREAM *prb_file )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (sd->nErrorCode)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         AddErrorMessage( sd->pStrErrStruct, ErrMsg( sd->nErrorCode ) );
    // INCHI✔️❌:         inchi_ios_eprint( log_file,
    // INCHI✔️❌:                           "Error %d (%s) structure #%ld component %d.%s%s%s%s\n",
    // INCHI✔️❌:                           sd->nErrorCode, sd->pStrErrStruct,
    // INCHI✔️❌:                           num_inp, i + 1, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
    // INCHI✔️❌:         sd->nErrorType = ( sd->nErrorCode == CT_OUT_OF_RAM || sd->nErrorCode == CT_USER_QUIT_ERR )
    // INCHI✔️❌:             ? _IS_FATAL
    // INCHI✔️❌:             : _IS_ERROR;
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔️❌:         if (( ip->bINChIOutputOptions & INCHI_OUT_WINCHI_WINDOW ) &&
    // INCHI✔️❌:             ( ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sd->nErrorType = ProcessStructError( out_file, log_file,
    // INCHI✔️❌:                                                  sd->pStrErrStruct, sd->nErrorType,
    // INCHI✔️❌:                                                  num_inp, ip );
    // INCHI✔️❌:             /*  Save the problem structure */
    // INCHI✔️❌:             if (prb_file->f &&
    // INCHI✔️❌:                  0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd &&
    // INCHI✔️❌:                  !ip->bSaveAllGoodStructsAsProblem)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 MolfileSaveCopy( inp_file, sd->fPtrStart, sd->fPtrEnd, prb_file->f, num_inp );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  Save the problem structure */
    // INCHI✔️❌:             if (sd->nErrorCode &&
    // INCHI✔️❌:                  prb_file->f &&
    // INCHI✔️❌:                  0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd &&
    // INCHI✔️❌:                  !ip->bSaveAllGoodStructsAsProblem)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 MolfileSaveCopy( inp_file, sd->fPtrStart, sd->fPtrEnd, prb_file->f, num_inp );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: /* #ifndef TARGET_API_LIB */
    // INCHI✔️❌: #if ( !defined( TARGET_API_LIB ) && !defined(TARGET_EXE_STANDALONE) )
    // INCHI✔️❌:     /*  print the logfile record */
    // INCHI✔️❌:     if (log_file->f && log_file->f != stderr && ( sd->ulStructTime >= 1000 || sd->nErrorCode ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         fprintf( log_file->f, "%10lu msec structure #%ld.%s%s%s%s (%d component%s, %d atom%s, error=%d).\n",
    // INCHI✔️❌:                 sd->ulStructTime, num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ),
    // INCHI✔️❌:                 orig_inp_data->num_components, orig_inp_data->num_components == 1 ? "" : "s",
    // INCHI✔️❌:                 orig_inp_data->num_inp_atoms, orig_inp_data->num_inp_atoms == 1 ? "" : "s", sd->nErrorCode );
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     return sd->nErrorType;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: TreatErrorsInCreateOneComponentINChI
    // END INCHI C FUNCTION: TreatErrorsInCreateOneComponentINChI
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: TreatErrorsInCreateOneComponentINChI
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; TARGET_LIB_FOR_WINCHI is undefined.
    // INCHI✔️❌: Both problem-file branches and the !TARGET_API_LIB timing record are inactive.
    // INCHI✔️❌: SourceHeap-backed formatting adds temporary allocations versus direct C varargs formatting.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: TreatErrorsInCreateOneComponentINChI
    if structure_data.nErrorCode != 0 {
        let error_message = ErrMsg(structure_data.nErrorCode);
        let mut c_error_message = error_message.bytes().map(|byte| byte as i8).collect::<Vec<_>>();
        c_error_message.push(0);
        AddErrorMessage(Some(&mut structure_data.pStrErrStruct), Some(&c_error_message))?;

        let error_text = source_array_c_bytes(&structure_data.pStrErrStruct)?;
        let mut log = b"Error ".to_vec();
        log.extend_from_slice(structure_data.nErrorCode.to_string().as_bytes());
        log.extend_from_slice(b" (");
        log.extend_from_slice(&error_text);
        log.extend_from_slice(b") structure #");
        log.extend_from_slice(input_structure_number.to_string().as_bytes());
        log.extend_from_slice(b" component ");
        log.extend_from_slice(component_index.wrapping_add(1).to_string().as_bytes());
        log.push(b'.');
        log.extend_from_slice(&sdf_label_value(heap, input_parameters)?);
        log.push(b'\n');
        eprint_bytes(heap, log_file, &log)?;

        structure_data.nErrorType =
            if structure_data.nErrorCode == CT_OUT_OF_RAM || structure_data.nErrorCode == CT_USER_QUIT_ERR {
                _IS_FATAL as i32
            } else {
                _IS_ERROR as i32
            };
    }
    Ok(structure_data.nErrorType)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn TreatCreateINChIWarning(
    heap: &mut SourceHeap,
    structure_data: &mut STRUCT_DATA,
    input_parameters: &INPUT_PARMS,
    _original_atom_data: &ORIG_ATOM_DATA,
    input_structure_number: i64,
    input_file: Option<&mut INCHI_IOSTREAM>,
    log_file: Option<&mut INCHI_IOSTREAM>,
    _output_file: Option<&mut INCHI_IOSTREAM>,
    problem_file: Option<&mut INCHI_IOSTREAM>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1453 TreatCreateINChIWarning
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: TreatCreateINChIWarning
    // INCHI✔️❌: int TreatCreateINChIWarning( STRUCT_DATA    *sd,
    // INCHI✔️❌:                              INPUT_PARMS    *ip,
    // INCHI✔️❌:                              ORIG_ATOM_DATA *orig_inp_data,
    // INCHI✔️❌:                              long           num_inp,
    // INCHI✔️❌:                              INCHI_IOSTREAM *inp_file,
    // INCHI✔️❌:                              INCHI_IOSTREAM *log_file,
    // INCHI✔️❌:                              INCHI_IOSTREAM *out_file,
    // INCHI✔️❌:                              INCHI_IOSTREAM *prb_file )
    // INCHI✔️❌: {
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( bRELEASE_VERSION == 0 && (EXTR_FLAGS || EXTR_MASK) )
    // INCHI✔️❌:     if (EXTR_MASK ? ( ( sd->bExtract & EXTR_MASK ) == EXTR_FLAGS )
    // INCHI✔️❌:                     : ( sd->bExtract & EXTR_FLAGS )
    // INCHI✔️❌:        )
    // INCHI✔️❌:     {
    // INCHI✔️❌:         char szMsg[64];
    // INCHI✔️❌:         sprintf( szMsg, "ExtractStruct.code=0x%X", sd->bExtract );
    // INCHI✔️❌:         if (!ip->bNoWarnings)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             WarningMessage( sd->pStrErrStruct, szMsg );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!sd->nErrorCode && sd->pStrErrStruct[0])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint( log_file, "Warning (%s) structure #%ld.%s%s%s%s\n",
    // INCHI✔️❌:             sd->pStrErrStruct, num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
    // INCHI✔️❌:         sd->nErrorType = _IS_WARNING;
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔️❌:         if (( ip->bINChIOutputOptions & INCHI_OUT_WINCHI_WINDOW ) &&
    // INCHI✔️❌:             ( ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sd->nErrorType = ProcessStructError( out_file,
    // INCHI✔️❌:                                                  log_file,
    // INCHI✔️❌:                                                  sd->pStrErrStruct,
    // INCHI✔️❌:                                                  sd->nErrorType,
    // INCHI✔️❌:                                                  num_inp,
    // INCHI✔️❌:                                                  ip );
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         /*  save the structure as a problem structure if requested */
    // INCHI✔️❌:         if (ip->bSaveWarningStructsAsProblem && !ip->bSaveAllGoodStructsAsProblem &&
    // INCHI✔️❌:              prb_file->f && 0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd)
    // INCHI✔️❌:         {   /* djb-rwth: addressing coverity ID #499545 -- return values handled properly */
    // INCHI✔️❌:             MolfileSaveCopy( inp_file,
    // INCHI✔️❌:                              sd->fPtrStart,
    // INCHI✔️❌:                              sd->fPtrEnd,
    // INCHI✔️❌:                              prb_file->f,
    // INCHI✔️❌:                              num_inp );
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( bRELEASE_VERSION == 0 )
    // INCHI✔️❌:         /*  otherwise extract the structure as a problem structure if requested */
    // INCHI✔️❌:         else
    // INCHI✔️❌:             if (( EXTR_MASK ? ( ( sd->bExtract & EXTR_MASK ) == EXTR_FLAGS ) : ( sd->bExtract & EXTR_FLAGS ) )
    // INCHI✔️❌:                             && !ip->bSaveAllGoodStructsAsProblem &&
    // INCHI✔️❌:                                prb_file->f &&
    // INCHI✔️❌:                                0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 MolfileSaveCopy( inp_file->f,
    // INCHI✔️❌:                                  sd->fPtrStart,
    // INCHI✔️❌:                                  sd->fPtrEnd,
    // INCHI✔️❌:                                  prb_file->f,
    // INCHI✔️❌:                                  num_inp );
    // INCHI✔️❌:             }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( bRELEASE_VERSION != 1 && bOUTPUT_ONE_STRUCT_TIME == 1 )
    // INCHI✔️❌: #ifndef TARGET_API_LIB
    // INCHI✔️❌:     if (log_file && log_file != stderr)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         fprintf( log_file, "%10lu msec structure %1dD #%ld.%s%s%s%s (%d component%s, %d atom%s, error=%d).\n",
    // INCHI✔️❌:                 sd->ulStructTime, orig_inp_data->num_dimensions, num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ),
    // INCHI✔️❌:                 orig_inp_data->num_components, orig_inp_data->num_components == 1 ? "" : "s",
    // INCHI✔️❌:                 orig_inp_data->num_inp_atoms, orig_inp_data->num_inp_atoms == 1 ? "" : "s", sd->nErrorCode );
    // INCHI✔️❌:     }
    // INCHI✔️❌: #else
    // INCHI✔️❌:     if (log_file)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint( log_file, "%10lu msec structure %1dD #%ld.%s%s%s%s (%d component%s, %d atom%s, error=%d).\n",
    // INCHI✔️❌:                 sd->ulStructTime, orig_inp_data->num_dimensions, num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ),
    // INCHI✔️❌:                 orig_inp_data->num_components, orig_inp_data->num_components == 1 ? "" : "s",
    // INCHI✔️❌:                 orig_inp_data->num_inp_atoms, orig_inp_data->num_inp_atoms == 1 ? "" : "s", sd->nErrorCode );
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     return sd->nErrorType;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: TreatCreateINChIWarning
    // END INCHI C FUNCTION: TreatCreateINChIWarning
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: TreatCreateINChIWarning
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; bRELEASE_VERSION=1.
    // INCHI✔️❌: EXTR_FLAGS=0; EXTR_MASK=0; TARGET_LIB_FOR_WINCHI is undefined;
    // INCHI✔️❌: bOUTPUT_ONE_STRUCT_TIME is forced to 0 for release builds.
    // INCHI✔️❌: SourceHeap-backed formatting adds temporary allocations versus direct C varargs formatting.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: TreatCreateINChIWarning

    if structure_data.nErrorCode == 0 && structure_data.pStrErrStruct[0] != 0 {
        let warning_text = source_array_c_bytes(&structure_data.pStrErrStruct)?;
        let mut log = b"Warning (".to_vec();
        log.extend_from_slice(&warning_text);
        log.extend_from_slice(b") structure #");
        log.extend_from_slice(input_structure_number.to_string().as_bytes());
        log.push(b'.');
        log.extend_from_slice(&sdf_label_value(heap, input_parameters)?);
        log.push(b'\n');
        eprint_bytes(heap, log_file, &log)?;
        structure_data.nErrorType = _IS_WARNING as i32;

        let problem_output = problem_file
            .as_ref()
            .map(|stream| stream.f)
            .unwrap_or_else(SourceMutPointer::null);
        if input_parameters.bSaveWarningStructsAsProblem != 0
            && input_parameters.bSaveAllGoodStructsAsProblem == 0
            && !problem_output.is_null()
            && structure_data.fPtrStart >= 0
            && structure_data.fPtrStart < structure_data.fPtrEnd
        {
            let _ = MolfileSaveCopy(
                heap,
                input_file,
                structure_data.fPtrStart,
                structure_data.fPtrEnd,
                problem_output,
                input_structure_number,
            )?;
        }
    }

    Ok(structure_data.nErrorType)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn DisplayTheWholeStructure(
    _canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    _clock: SourceMutPointer<INCHI_CLOCK>,
    _structure_data: SourceMutPointer<STRUCT_DATA>,
    _input_parameters: SourceMutPointer<INPUT_PARMS>,
    _title: SourceMutPointer<i8>,
    _input_file: SourceMutPointer<INCHI_IOSTREAM>,
    _log_file: SourceMutPointer<INCHI_IOSTREAM>,
    _original_atom_data: SourceMutPointer<ORIG_ATOM_DATA>,
    _input_structure_number: i64,
    _inchi_index: i32,
    _show_structure: i32,
    _inchi_library_flag: i32,
) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1204 DisplayTheWholeStructure
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: DisplayTheWholeStructure
    // INCHI✔️✔️: int DisplayTheWholeStructure( struct tagCANON_GLOBALS *pCG,
    // INCHI✔️✔️:                               struct tagINCHI_CLOCK   *ic,
    // INCHI✔️✔️:                               STRUCT_DATA             *sd,
    // INCHI✔️✔️:                               INPUT_PARMS             *ip,
    // INCHI✔️✔️:                               char                    *szTitle,
    // INCHI✔️✔️:                               INCHI_IOSTREAM          *inp_file,
    // INCHI✔️✔️:                               INCHI_IOSTREAM          *log_file,
    // INCHI✔️✔️:                               ORIG_ATOM_DATA          *orig_inp_data,
    // INCHI✔️✔️:                               long                    num_inp,
    // INCHI✔️✔️:                               int                     iINChI,
    // INCHI✔️✔️:                               int                     bShowStruct,
    // INCHI✔️✔️:                               int                     bINCHI_LIB_Flag )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: DisplayTheWholeStructure
    // END INCHI C FUNCTION: DisplayTheWholeStructure
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: DisplayTheWholeStructure
    // INCHI✔️✔️: COMPILE_ANSI_ONLY selects the dummy implementation; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️✔️: The GUI implementation at runichi4.c:994 is excluded by the active preprocessor branch.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: DisplayTheWholeStructure
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{
        _IS_WARNING, AMBIGUOUS_STEREO_ERROR, INCHI_IOS_TYPE_FILE, INCHI_IOS_TYPE_STRING, INCHI_OUT_NO_AUX_INFO,
        INChI_Aux, SourceFile, inp_ATOM,
    };

    fn buffer(heap: &mut SourceHeap, size: usize) -> INCHI_IOS_STRING {
        INCHI_IOS_STRING {
            pStr: heap.allocate_model_storage(vec![0_i8; size]).unwrap(),
            nAllocatedLength: size as i32,
            nPtr: size as i32,
            ..INCHI_IOS_STRING::default()
        }
    }

    fn stream(heap: &mut SourceHeap) -> INCHI_IOSTREAM {
        INCHI_IOSTREAM {
            s: buffer(heap, 4096),
            type_: INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        }
    }

    fn text(heap: &SourceHeap, stream: &INCHI_IOSTREAM) -> String {
        let bytes = heap.slice(stream.s.pStr.as_const()).unwrap();
        let end = bytes.iter().position(|byte| *byte == 0).unwrap();
        bytes[..end].iter().map(|byte| *byte as u8 as char).collect()
    }

    fn inchi(heap: &mut SourceHeap, charge: i32) -> SourceMutPointer<INChI> {
        let value = INChI {
            nTotalCharge: charge,
            nNumberOfAtoms: 1,
            szHillFormula: heap.allocate_model_storage(vec![b'C' as i8, 0]).unwrap(),
            nAtom: heap.allocate_model_storage(vec![6_u8]).unwrap(),
            nNum_H: heap.allocate_model_storage(vec![0_i8]).unwrap(),
            ..INChI::default()
        };
        heap.allocate_model_storage(vec![value]).unwrap()
    }

    #[test]
    fn source_port__runichi4__sortandprintinchi__line_102() {
        let mut heap = SourceHeap::default();
        let left_fixed = inchi(&mut heap, 1);
        let left_mobile = inchi(&mut heap, 1);
        let right_fixed = inchi(&mut heap, -1);
        let right_mobile = inchi(&mut heap, -1);
        let component_rows = heap
            .allocate_model_storage(vec![[left_fixed, left_mobile], [right_fixed, right_mobile]])
            .unwrap();
        let aux_left = heap.allocate_model_storage(vec![INChI_Aux::default()]).unwrap();
        let aux_right = heap.allocate_model_storage(vec![INChI_Aux::default()]).unwrap();
        let aux_rows = heap
            .allocate_model_storage(vec![
                [aux_left, SourceMutPointer::null()],
                [aux_right, SourceMutPointer::null()],
            ])
            .unwrap();
        let mut output = stream(&mut heap);
        let mut log = stream(&mut heap);
        let mut composite = std::array::from_fn(|_| std::array::from_fn(|_| COMP_ATOM_DATA::default()));
        let mut taut_flags = [0, 0];
        let mut taut_flags_done = [TG_FLAG_DISCONNECT_COORD_DONE as u64, 0];
        let normalization = NORM_CANON_FLAGS {
            bTautFlags: [[1, 2], [4, 8]],
            bTautFlagsDone: [[16, 32], [64, 128]],
            ..NORM_CANON_FLAGS::default()
        };
        assert_eq!(
            SortAndPrintINChI(
                &mut heap,
                SourceMutPointer::null(),
                &mut output,
                None,
                &mut log,
                &INPUT_PARMS::default(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                Some(&mut composite),
                SourceConstPointer::null(),
                &[2, 0],
                &[2, 0],
                &[2, 0],
                &mut taut_flags,
                &mut taut_flags_done,
                &normalization,
                1,
                [component_rows, SourceMutPointer::null()],
                [aux_rows, SourceMutPointer::null()],
                SourceMutPointer::null(),
                0,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        let reordered = heap.slice(component_rows.as_const()).unwrap();
        assert_eq!(reordered[0][TAUT_YES as usize], right_mobile);
        assert_eq!(reordered[1][TAUT_YES as usize], left_mobile);
        assert_eq!(taut_flags, [3, 12]);
        assert_eq!(
            taut_flags_done,
            [TG_FLAG_DISCONNECT_COORD_DONE as u64 | 16 | 32, 64 | 128]
        );

        let only_fixed = inchi(&mut heap, 0);
        let fixed_rows = heap
            .allocate_model_storage(vec![[only_fixed, SourceMutPointer::null()]])
            .unwrap();
        let fixed_aux = heap
            .allocate_model_storage(vec![[aux_left, SourceMutPointer::null()]])
            .unwrap();
        assert_eq!(
            SortAndPrintINChI(
                &mut heap,
                SourceMutPointer::null(),
                &mut output,
                None,
                &mut log,
                &INPUT_PARMS::default(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                Some(&mut composite),
                SourceConstPointer::null(),
                &[1, 0],
                &[1, 0],
                &[0, 0],
                &mut [0, 0],
                &mut [0, 0],
                &NORM_CANON_FLAGS::default(),
                2,
                [fixed_rows, SourceMutPointer::null()],
                [fixed_aux, SourceMutPointer::null()],
                SourceMutPointer::null(),
                0,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        let swapped = heap.slice(fixed_rows.as_const()).unwrap()[0];
        assert_eq!(swapped[TAUT_YES as usize], only_fixed);
        assert!(swapped[TAUT_NON as usize].is_null());

        let mut failing_heap = SourceHeap::default();
        failing_heap.fail_after_allocations(0);
        let mut failing_output = INCHI_IOSTREAM::default();
        let mut failing_log = INCHI_IOSTREAM::default();
        let mut failing_composite = std::array::from_fn(|_| std::array::from_fn(|_| COMP_ATOM_DATA::default()));
        assert_eq!(
            SortAndPrintINChI(
                &mut failing_heap,
                SourceMutPointer::null(),
                &mut failing_output,
                None,
                &mut failing_log,
                &INPUT_PARMS::default(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                Some(&mut failing_composite),
                SourceConstPointer::null(),
                &[1, 0],
                &[1, 0],
                &[1, 0],
                &mut [0, 0],
                &mut [0, 0],
                &NORM_CANON_FLAGS::default(),
                3,
                [SourceMutPointer::null(), SourceMutPointer::null()],
                [SourceMutPointer::null(), SourceMutPointer::null()],
                SourceMutPointer::null(),
                0,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(failing_heap.source_allocation_calls(), 2);

        let mut printing_heap = SourceHeap::default();
        let mut printing_output = stream(&mut printing_heap);
        let mut printing_log = stream(&mut printing_heap);
        let mut printing_buffer = buffer(&mut printing_heap, 4096);
        let mut printing_composite = std::array::from_fn(|_| std::array::from_fn(|_| COMP_ATOM_DATA::default()));
        let flags = printing_heap.allocate_model_storage(vec![0_i32]).unwrap();
        let print_options =
            (INCHI_OUT_NO_AUX_INFO | INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS | INCHI_OUT_TABBED_OUTPUT)
                as i32;
        assert_eq!(
            SortAndPrintINChI(
                &mut printing_heap,
                SourceMutPointer::null(),
                &mut printing_output,
                Some(&mut printing_buffer),
                &mut printing_log,
                &INPUT_PARMS {
                    bINChIOutputOptions: print_options,
                    ..INPUT_PARMS::default()
                },
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                Some(&mut printing_composite),
                SourceConstPointer::null(),
                &[0, 0],
                &[0, 0],
                &[0, 0],
                &mut [0, 0],
                &mut [0, 0],
                &NORM_CANON_FLAGS::default(),
                i64::from(i32::MAX) + 1,
                [SourceMutPointer::null(), SourceMutPointer::null()],
                [SourceMutPointer::null(), SourceMutPointer::null()],
                flags,
                0,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(
            text(&printing_heap, &printing_output),
            concat!("InChI=1//\n", "\n==== InChI ANNOTATED CONTENTS ====\n", "\nInChI=\n1\n")
        );
        assert_eq!(text(&printing_heap, &printing_log), "");
    }

    #[test]
    fn source_port__runichi4__freeinchiarrays__line_1357() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            FreeINChIArrays(&mut heap, SourceMutPointer::null(), SourceMutPointer::null(), i32::MIN,),
            Ok(())
        );

        let released_inchi = heap.allocate(vec![INChI::default()]).unwrap();
        let retained_inchi = heap
            .allocate(vec![INChI {
                nRefCount: 1,
                ..INChI::default()
            }])
            .unwrap();
        let inchi_rows = heap
            .allocate(vec![
                [released_inchi, retained_inchi],
                [SourceMutPointer::null(), SourceMutPointer::null()],
            ])
            .unwrap();
        let released_aux = heap.allocate(vec![INChI_Aux::default()]).unwrap();
        let retained_aux = heap
            .allocate(vec![INChI_Aux {
                nRefCount: 1,
                ..INChI_Aux::default()
            }])
            .unwrap();
        let aux_rows = heap
            .allocate(vec![
                [released_aux, retained_aux],
                [SourceMutPointer::null(), SourceMutPointer::null()],
            ])
            .unwrap();
        FreeINChIArrays(&mut heap, inchi_rows, aux_rows, 2).unwrap();
        let resulting_inchi = heap.slice(inchi_rows.as_const()).unwrap();
        assert!(resulting_inchi[0][TAUT_NON as usize].is_null());
        assert_eq!(resulting_inchi[0][TAUT_YES as usize], retained_inchi);
        assert_eq!(heap.slice(retained_inchi.as_const()).unwrap()[0].nRefCount, 0);
        let resulting_aux = heap.slice(aux_rows.as_const()).unwrap();
        assert!(resulting_aux[0][TAUT_NON as usize].is_null());
        assert_eq!(resulting_aux[0][TAUT_YES as usize], retained_aux);
        assert_eq!(heap.slice(retained_aux.as_const()).unwrap()[0].nRefCount, 0);
        assert_eq!(
            heap.slice(released_inchi.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(released_aux.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__runichi4__freeallinchiarrays__line_1322() {
        let mut heap = SourceHeap::default();
        let inchi = heap.allocate(vec![INChI::default()]).unwrap();
        let inchi_rows = heap.allocate(vec![[inchi, SourceMutPointer::null()]]).unwrap();
        let aux = heap.allocate(vec![INChI_Aux::default()]).unwrap();
        let aux_rows = heap.allocate(vec![[aux, SourceMutPointer::null()]]).unwrap();
        let zero_inchi_rows = heap
            .allocate(vec![[SourceMutPointer::null(), SourceMutPointer::null()]])
            .unwrap();
        let zero_aux_rows = heap
            .allocate(vec![[SourceMutPointer::null(), SourceMutPointer::null()]])
            .unwrap();
        let mut top_inchi = [inchi_rows, zero_inchi_rows];
        let mut top_aux = [aux_rows, zero_aux_rows];
        let mut counts = [1, 0];
        FreeAllINChIArrays(&mut heap, &mut top_inchi, &mut top_aux, &mut counts).unwrap();
        assert_eq!(counts, [0, 0]);
        assert!(top_inchi[0].is_null());
        assert!(top_aux[0].is_null());
        assert_eq!(top_inchi[1], zero_inchi_rows);
        assert_eq!(top_aux[1], zero_aux_rows);
        assert!(heap.slice(zero_inchi_rows.as_const()).is_ok());
        assert!(heap.slice(zero_aux_rows.as_const()).is_ok());
        assert_eq!(heap.slice(inchi.as_const()), Err(SourceHeapError::MissingAllocation));
        assert_eq!(heap.slice(aux.as_const()), Err(SourceHeapError::MissingAllocation));

        let negative_inchi = heap
            .allocate(vec![[SourceMutPointer::null(), SourceMutPointer::null()]])
            .unwrap();
        let negative_aux = heap
            .allocate(vec![[SourceMutPointer::null(), SourceMutPointer::null()]])
            .unwrap();
        let mut negative_top_inchi = [negative_inchi, SourceMutPointer::null()];
        let mut negative_top_aux = [negative_aux, SourceMutPointer::null()];
        let mut negative_counts = [i32::MIN, 0];
        FreeAllINChIArrays(
            &mut heap,
            &mut negative_top_inchi,
            &mut negative_top_aux,
            &mut negative_counts,
        )
        .unwrap();
        assert!(negative_top_inchi[0].is_null());
        assert!(negative_top_aux[0].is_null());
        assert_eq!(
            heap.slice(negative_inchi.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(negative_aux.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__runichi4__getprocessingwarningsoneinchi__line_1574() {
        fn error_text(buffer: &[i8]) -> String {
            let end = buffer.iter().position(|byte| *byte == 0).unwrap();
            buffer[..end].iter().map(|byte| *byte as u8 as char).collect()
        }

        fn run(
            heap: &SourceHeap,
            atoms: SourceMutPointer<inp_ATOM>,
            atom_count: i32,
            no_warnings: i32,
        ) -> (i32, String) {
            let mut errors = vec![0_i8; crate::source_types::STR_ERR_LEN as usize];
            let result = GetProcessingWarningsOneInChI(
                heap,
                &INChI {
                    nNumberOfAtoms: atom_count,
                    ..INChI::default()
                },
                &INP_ATOM_DATA {
                    at: atoms,
                    ..INP_ATOM_DATA::default()
                },
                &mut errors,
                no_warnings,
            )
            .unwrap();
            (result, error_text(&errors))
        }

        let mut heap = SourceHeap::default();
        assert_eq!(run(&heap, SourceMutPointer::null(), 4, 0), (0, String::new()));

        let atoms = heap
            .allocate_model_storage(vec![
                inp_ATOM {
                    bAmbiguousStereo: AMBIGUOUS_STEREO_ATOM as i8,
                    ..inp_ATOM::default()
                },
                inp_ATOM {
                    bAmbiguousStereo: AMBIGUOUS_STEREO_ATOM_ISO as i8,
                    ..inp_ATOM::default()
                },
                inp_ATOM {
                    bAmbiguousStereo: AMBIGUOUS_STEREO_BOND as i8,
                    ..inp_ATOM::default()
                },
                inp_ATOM {
                    bAmbiguousStereo: AMBIGUOUS_STEREO_BOND_ISO as i8,
                    ..inp_ATOM::default()
                },
                inp_ATOM {
                    bAmbiguousStereo: (AMBIGUOUS_STEREO_ATOM | AMBIGUOUS_STEREO_BOND) as i8,
                    ..inp_ATOM::default()
                },
                inp_ATOM {
                    bAmbiguousStereo: AMBIGUOUS_STEREO_ERROR as i8,
                    ..inp_ATOM::default()
                },
            ])
            .unwrap();

        assert_eq!(run(&heap, atoms, 0, 0), (0, String::new()));
        assert_eq!(run(&heap, atoms, -1, 0), (0, String::new()));
        assert_eq!(run(&heap, atoms, 1, 0), (1, "Ambiguous stereo: center(s)".to_owned()));
        assert_eq!(
            run(&heap, atoms.offset(1).unwrap(), 1, 0),
            (1, "Ambiguous stereo: center(s)".to_owned())
        );
        assert_eq!(
            run(&heap, atoms.offset(2).unwrap(), 1, 0),
            (1, "Ambiguous stereo: bond(s)".to_owned())
        );
        assert_eq!(
            run(&heap, atoms.offset(3).unwrap(), 1, 0),
            (1, "Ambiguous stereo: bond(s)".to_owned())
        );
        assert_eq!(
            run(&heap, atoms, 5, 0),
            (1, "Ambiguous stereo: center(s); bond(s)".to_owned())
        );
        assert_eq!(run(&heap, atoms.offset(5).unwrap(), 1, 0), (0, String::new()));
        assert_eq!(run(&heap, atoms, 5, 1), (1, String::new()));
        assert_eq!(run(&heap, atoms, 5, i32::MIN), (1, String::new()));
    }

    #[test]
    fn source_port__runichi4__getprocessingwarningsonecomponentinchi__line_1550() {
        fn error_text(buffer: &[i8]) -> String {
            let end = buffer.iter().position(|byte| *byte == 0).unwrap();
            buffer[..end].iter().map(|byte| *byte as u8 as char).collect()
        }

        let mut heap = SourceHeap::default();
        let mut structure = STRUCT_DATA::default();
        assert_eq!(
            GetProcessingWarningsOneComponentInChI(
                &heap,
                &[SourceMutPointer::null(), SourceMutPointer::null()],
                &std::array::from_fn(|_| INP_ATOM_DATA::default()),
                &mut structure,
                0,
            ),
            Ok(0)
        );

        let empty_inchi = heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 0,
                ..INChI::default()
            }])
            .unwrap();
        assert_eq!(
            GetProcessingWarningsOneComponentInChI(
                &heap,
                &[empty_inchi, SourceMutPointer::null()],
                &std::array::from_fn(|_| INP_ATOM_DATA::default()),
                &mut structure,
                0,
            ),
            Ok(0)
        );

        let atom = heap
            .allocate_model_storage(vec![inp_ATOM {
                bAmbiguousStereo: AMBIGUOUS_STEREO_ATOM_ISO as i8,
                ..inp_ATOM::default()
            }])
            .unwrap();
        let bond = heap
            .allocate_model_storage(vec![inp_ATOM {
                bAmbiguousStereo: AMBIGUOUS_STEREO_BOND_ISO as i8,
                ..inp_ATOM::default()
            }])
            .unwrap();
        let normalized = [
            INP_ATOM_DATA {
                at: atom,
                ..INP_ATOM_DATA::default()
            },
            INP_ATOM_DATA {
                at: bond,
                ..INP_ATOM_DATA::default()
            },
        ];
        let first = heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                ..INChI::default()
            }])
            .unwrap();
        let second = heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                ..INChI::default()
            }])
            .unwrap();

        assert_eq!(
            GetProcessingWarningsOneComponentInChI(&heap, &[first, second], &normalized, &mut structure, 0,),
            Ok(1)
        );
        assert_eq!(
            error_text(&structure.pStrErrStruct),
            "Ambiguous stereo: center(s); bond(s)"
        );

        let mut no_warning_structure = STRUCT_DATA::default();
        assert_eq!(
            GetProcessingWarningsOneComponentInChI(
                &heap,
                &[SourceMutPointer::null(), second],
                &normalized,
                &mut no_warning_structure,
                1,
            ),
            Ok(1)
        );
        assert_eq!(error_text(&no_warning_structure.pStrErrStruct), "");
    }

    #[test]
    fn source_port__runichi4__displaythewholestructure__line_1204() {
        assert_eq!(
            DisplayTheWholeStructure(
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                i64::MIN,
                i32::MIN,
                i32::MAX,
                -1,
            ),
            0
        );

        let mut heap = SourceHeap::default();
        let canonical_globals = heap.allocate_model_storage(vec![CANON_GLOBALS::default()]).unwrap();
        let clock = heap.allocate_model_storage(vec![INCHI_CLOCK::default()]).unwrap();
        let structure = heap
            .allocate_model_storage(vec![STRUCT_DATA {
                bUserQuit: 17,
                ..STRUCT_DATA::default()
            }])
            .unwrap();
        let parameters = heap
            .allocate_model_storage(vec![INPUT_PARMS {
                bDisplay: 23,
                ..INPUT_PARMS::default()
            }])
            .unwrap();
        let title = heap.allocate_model_storage(vec![b'X' as i8, 0]).unwrap();
        let input_file = heap.allocate_model_storage(vec![INCHI_IOSTREAM::default()]).unwrap();
        let log_file = heap.allocate_model_storage(vec![INCHI_IOSTREAM::default()]).unwrap();
        let original = heap
            .allocate_model_storage(vec![ORIG_ATOM_DATA {
                num_inp_atoms: 31,
                ..ORIG_ATOM_DATA::default()
            }])
            .unwrap();

        assert_eq!(
            DisplayTheWholeStructure(
                canonical_globals,
                clock,
                structure,
                parameters,
                title,
                input_file,
                log_file,
                original,
                i64::MAX,
                i32::MAX,
                i32::MIN,
                1,
            ),
            0
        );
        assert_eq!(heap.slice(structure.as_const()).unwrap()[0].bUserQuit, 17);
        assert_eq!(heap.slice(parameters.as_const()).unwrap()[0].bDisplay, 23);
        assert_eq!(heap.slice(original.as_const()).unwrap()[0].num_inp_atoms, 31);
        assert_eq!(heap.slice(title.as_const()).unwrap(), &[b'X' as i8, 0]);
    }

    #[test]
    fn source_port__runichi4__treaterrorsincreateonecomponentinchi__line_1387() {
        let mut no_error_heap = SourceHeap::default();
        let mut no_error_log = stream(&mut no_error_heap);
        let mut no_error_structure = STRUCT_DATA {
            nErrorType: _IS_WARNING as i32,
            ..STRUCT_DATA::default()
        };
        assert_eq!(
            TreatErrorsInCreateOneComponentINChI(
                &mut no_error_heap,
                &mut no_error_structure,
                &INPUT_PARMS::default(),
                &ORIG_ATOM_DATA::default(),
                0,
                1,
                None,
                Some(&mut no_error_log),
                None,
                None,
            ),
            Ok(_IS_WARNING as i32)
        );
        assert_eq!(text(&no_error_heap, &no_error_log), "");

        let mut ordinary_heap = SourceHeap::default();
        let label = ordinary_heap
            .allocate_model_storage(vec![b'I' as i8, b'D' as i8, 0])
            .unwrap();
        let value = ordinary_heap
            .allocate_model_storage(vec![b'4' as i8, b'2' as i8, 0])
            .unwrap();
        let parameters = INPUT_PARMS {
            pSdfLabel: label,
            pSdfValue: value,
            ..INPUT_PARMS::default()
        };
        let mut ordinary_log = stream(&mut ordinary_heap);
        let mut ordinary_structure = STRUCT_DATA {
            nErrorCode: crate::source_types::CT_ATOMCOUNT_ERR,
            ..STRUCT_DATA::default()
        };
        ordinary_structure.pStrErrStruct[..6]
            .copy_from_slice(&[b'p' as i8, b'r' as i8, b'i' as i8, b'o' as i8, b'r' as i8, 0]);
        assert_eq!(
            TreatErrorsInCreateOneComponentINChI(
                &mut ordinary_heap,
                &mut ordinary_structure,
                &parameters,
                &ORIG_ATOM_DATA::default(),
                2,
                -7,
                None,
                Some(&mut ordinary_log),
                None,
                None,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(
            text(&ordinary_heap, &ordinary_log),
            "Error -30011 (prior; ATOMCOUNT_ERR) structure #-7 component 3. ID=42\n"
        );
        let ordinary_error = ordinary_structure
            .pStrErrStruct
            .iter()
            .take_while(|byte| **byte != 0)
            .map(|byte| *byte as u8 as char)
            .collect::<String>();
        assert_eq!(ordinary_error, "prior; ATOMCOUNT_ERR");

        for fatal_code in [CT_OUT_OF_RAM, CT_USER_QUIT_ERR] {
            let mut fatal_heap = SourceHeap::default();
            let mut fatal_log = stream(&mut fatal_heap);
            let mut fatal_structure = STRUCT_DATA {
                nErrorCode: fatal_code,
                ..STRUCT_DATA::default()
            };
            assert_eq!(
                TreatErrorsInCreateOneComponentINChI(
                    &mut fatal_heap,
                    &mut fatal_structure,
                    &INPUT_PARMS::default(),
                    &ORIG_ATOM_DATA::default(),
                    0,
                    i64::MAX,
                    None,
                    Some(&mut fatal_log),
                    None,
                    None,
                ),
                Ok(_IS_FATAL as i32)
            );
            let expected_message = ErrMsg(fatal_code);
            assert_eq!(
                text(&fatal_heap, &fatal_log),
                format!(
                    "Error {fatal_code} ({expected_message}) structure #{} component 1.\n",
                    i64::MAX
                )
            );
        }

        let mut duplicate_heap = SourceHeap::default();
        let mut duplicate_log = stream(&mut duplicate_heap);
        let mut duplicate_structure = STRUCT_DATA {
            nErrorCode: CT_OUT_OF_RAM,
            ..STRUCT_DATA::default()
        };
        duplicate_structure.pStrErrStruct[..11].copy_from_slice(&[
            b'O' as i8, b'u' as i8, b't' as i8, b' ' as i8, b'o' as i8, b'f' as i8, b' ' as i8, b'R' as i8, b'A' as i8,
            b'M' as i8, 0,
        ]);
        assert_eq!(
            TreatErrorsInCreateOneComponentINChI(
                &mut duplicate_heap,
                &mut duplicate_structure,
                &INPUT_PARMS::default(),
                &ORIG_ATOM_DATA::default(),
                0,
                0,
                None,
                Some(&mut duplicate_log),
                None,
                None,
            ),
            Ok(_IS_FATAL as i32)
        );
        assert_eq!(
            text(&duplicate_heap, &duplicate_log),
            "Error -30002 (Out of RAM) structure #0 component 1.\n"
        );
    }

    #[test]
    fn source_port__runichi4__treatcreateinchiwarning__line_1453() {
        fn warning_structure(error_code: i32) -> STRUCT_DATA {
            let mut structure = STRUCT_DATA {
                nErrorCode: error_code,
                nErrorType: 37,
                fPtrStart: 0,
                fPtrEnd: 10,
                ..STRUCT_DATA::default()
            };
            structure.pStrErrStruct[..5].copy_from_slice(&[b'c' as i8, b'a' as i8, b'r' as i8, b'e' as i8, 0]);
            structure
        }

        let mut no_text_heap = SourceHeap::default();
        let mut no_text_log = stream(&mut no_text_heap);
        let mut no_text_structure = STRUCT_DATA {
            nErrorType: 41,
            ..STRUCT_DATA::default()
        };
        assert_eq!(
            TreatCreateINChIWarning(
                &mut no_text_heap,
                &mut no_text_structure,
                &INPUT_PARMS::default(),
                &ORIG_ATOM_DATA::default(),
                i64::MIN,
                None,
                Some(&mut no_text_log),
                None,
                None,
            ),
            Ok(41)
        );
        assert_eq!(text(&no_text_heap, &no_text_log), "");

        let mut error_heap = SourceHeap::default();
        let mut error_log = stream(&mut error_heap);
        let mut error_structure = warning_structure(9);
        assert_eq!(
            TreatCreateINChIWarning(
                &mut error_heap,
                &mut error_structure,
                &INPUT_PARMS::default(),
                &ORIG_ATOM_DATA::default(),
                i64::MAX,
                None,
                Some(&mut error_log),
                None,
                None,
            ),
            Ok(37)
        );
        assert_eq!(text(&error_heap, &error_log), "");

        let mut labeled_heap = SourceHeap::default();
        let label = labeled_heap
            .allocate_model_storage(vec![b'I' as i8, b'D' as i8, 0])
            .unwrap();
        let value = labeled_heap
            .allocate_model_storage(vec![b'4' as i8, b'2' as i8, 0])
            .unwrap();
        let parameters = INPUT_PARMS {
            pSdfLabel: label,
            pSdfValue: value,
            bNoWarnings: 1,
            ..INPUT_PARMS::default()
        };
        let mut labeled_log = stream(&mut labeled_heap);
        let mut labeled_structure = warning_structure(0);
        assert_eq!(
            TreatCreateINChIWarning(
                &mut labeled_heap,
                &mut labeled_structure,
                &parameters,
                &ORIG_ATOM_DATA::default(),
                i64::MIN,
                None,
                Some(&mut labeled_log),
                None,
                None,
            ),
            Ok(_IS_WARNING as i32)
        );
        assert_eq!(
            text(&labeled_heap, &labeled_log),
            format!("Warning (care) structure #{}. ID=42\n", i64::MIN)
        );

        let mut save_heap = SourceHeap::default();
        let input_file = save_heap
            .allocate_model_storage(vec![SourceFile {
                bytes: b"name\nbody\nrest\n".to_vec(),
                ..SourceFile::default()
            }])
            .unwrap();
        let problem_output = save_heap.allocate_model_storage(vec![SourceFile::default()]).unwrap();
        let mut input_stream = INCHI_IOSTREAM {
            f: input_file,
            type_: INCHI_IOS_TYPE_FILE as i32,
            ..INCHI_IOSTREAM::default()
        };
        let mut problem_stream = INCHI_IOSTREAM {
            f: problem_output,
            type_: INCHI_IOS_TYPE_FILE as i32,
            ..INCHI_IOSTREAM::default()
        };
        let mut save_log = stream(&mut save_heap);
        let mut save_structure = warning_structure(0);
        let save_parameters = INPUT_PARMS {
            bSaveWarningStructsAsProblem: 1,
            ..INPUT_PARMS::default()
        };
        assert_eq!(
            TreatCreateINChIWarning(
                &mut save_heap,
                &mut save_structure,
                &save_parameters,
                &ORIG_ATOM_DATA::default(),
                -17,
                Some(&mut input_stream),
                Some(&mut save_log),
                None,
                Some(&mut problem_stream),
            ),
            Ok(_IS_WARNING as i32)
        );
        assert_eq!(text(&save_heap, &save_log), "Warning (care) structure #-17.\n");
        assert_eq!(
            save_heap
                .slice(problem_output.as_const())
                .unwrap()
                .first()
                .unwrap()
                .bytes,
            b"#-17/name\nbody\n".to_vec()
        );

        for (save_warning, save_all, has_problem_file, start, end) in [
            (0, 0, true, 0, 10),
            (1, 1, true, 0, 10),
            (1, 0, false, 0, 10),
            (1, 0, true, -1, 10),
            (1, 0, true, 10, 10),
        ] {
            let mut heap = SourceHeap::default();
            let input = heap
                .allocate_model_storage(vec![SourceFile {
                    bytes: b"name\nbody\n".to_vec(),
                    ..SourceFile::default()
                }])
                .unwrap();
            let output = heap.allocate_model_storage(vec![SourceFile::default()]).unwrap();
            let mut input_stream = INCHI_IOSTREAM {
                f: input,
                type_: INCHI_IOS_TYPE_FILE as i32,
                ..INCHI_IOSTREAM::default()
            };
            let mut problem_stream = INCHI_IOSTREAM {
                f: if has_problem_file {
                    output
                } else {
                    SourceMutPointer::null()
                },
                type_: INCHI_IOS_TYPE_FILE as i32,
                ..INCHI_IOSTREAM::default()
            };
            let mut log = stream(&mut heap);
            let mut structure = warning_structure(0);
            structure.fPtrStart = start;
            structure.fPtrEnd = end;
            let parameters = INPUT_PARMS {
                bSaveWarningStructsAsProblem: save_warning,
                bSaveAllGoodStructsAsProblem: save_all,
                ..INPUT_PARMS::default()
            };
            assert_eq!(
                TreatCreateINChIWarning(
                    &mut heap,
                    &mut structure,
                    &parameters,
                    &ORIG_ATOM_DATA::default(),
                    3,
                    Some(&mut input_stream),
                    Some(&mut log),
                    None,
                    Some(&mut problem_stream),
                ),
                Ok(_IS_WARNING as i32)
            );
            assert!(heap.slice(output.as_const()).unwrap().first().unwrap().bytes.is_empty());
        }
    }

    #[test]
    fn source_port__runichi4__bisstructchiral__line_1271() {
        fn component_array(
            heap: &mut SourceHeap,
            first: SourceMutPointer<INChI>,
            second: SourceMutPointer<INChI>,
        ) -> SourceMutPointer<PINChI2> {
            heap.allocate_model_storage(vec![[first, second]]).unwrap()
        }

        let mut heap = SourceHeap::default();
        assert_eq!(
            bIsStructChiral(
                &heap,
                [SourceMutPointer::null(), SourceMutPointer::null()],
                [0, i32::MIN],
            ),
            Ok(0)
        );

        let deleted = heap
            .allocate_model_storage(vec![INChI {
                bDeleted: 1,
                nNumberOfAtoms: 3,
                ..INChI::default()
            }])
            .unwrap();
        let zero_atoms = heap.allocate_model_storage(vec![INChI::default()]).unwrap();
        let skipped = component_array(&mut heap, deleted, zero_atoms);
        assert_eq!(
            bIsStructChiral(&heap, [skipped, SourceMutPointer::null()], [1, 0],),
            Ok(0)
        );

        for stereo in [
            INChI_Stereo {
                nNumberOfStereoCenters: 1,
                nCompInv2Abs: 1,
                ..INChI_Stereo::default()
            },
            INChI_Stereo {
                t_parity: heap.allocate(vec![1_i8]).unwrap(),
                nNumberOfStereoCenters: 0,
                nCompInv2Abs: 1,
                ..INChI_Stereo::default()
            },
            INChI_Stereo {
                t_parity: heap.allocate(vec![1_i8]).unwrap(),
                nNumberOfStereoCenters: 1,
                nCompInv2Abs: 0,
                ..INChI_Stereo::default()
            },
        ] {
            let stereo = heap.allocate_model_storage(vec![stereo]).unwrap();
            let inchi = heap
                .allocate_model_storage(vec![INChI {
                    nNumberOfAtoms: 1,
                    Stereo: stereo,
                    ..INChI::default()
                }])
                .unwrap();
            let components = component_array(&mut heap, inchi, SourceMutPointer::null());
            assert_eq!(
                bIsStructChiral(&heap, [components, SourceMutPointer::null()], [1, 0],),
                Ok(0)
            );
        }

        let parity = heap.allocate(vec![i8::MIN]).unwrap();
        let valid_stereo = heap
            .allocate_model_storage(vec![INChI_Stereo {
                t_parity: parity,
                nNumberOfStereoCenters: i32::MAX,
                nCompInv2Abs: i32::MIN,
                ..INChI_Stereo::default()
            }])
            .unwrap();
        let ordinary = heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                Stereo: valid_stereo,
                ..INChI::default()
            }])
            .unwrap();
        let ordinary_components = component_array(&mut heap, SourceMutPointer::null(), ordinary);
        assert_eq!(
            bIsStructChiral(&heap, [ordinary_components, SourceMutPointer::null()], [1, 0],),
            Ok(1)
        );

        let isotopic = heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                StereoIsotopic: valid_stereo,
                ..INChI::default()
            }])
            .unwrap();
        let empty_component = [SourceMutPointer::null(), SourceMutPointer::null()];
        let isotopic_components = heap
            .allocate_model_storage(vec![empty_component, [isotopic, SourceMutPointer::null()]])
            .unwrap();
        assert_eq!(
            bIsStructChiral(&heap, [SourceMutPointer::null(), isotopic_components], [0, 2],),
            Ok(1)
        );
    }
}
