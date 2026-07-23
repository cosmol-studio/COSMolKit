use std::cmp::Ordering;

use crate::source::base::ichi_io::inchi_ios_print;
use crate::source::base::ichimake::{CompINChINonTaut2, CompINChITaut2};
use crate::source::base::ichiprt1::OutputINChI2;
use crate::source::base::strutil::{Free_INChI, Free_INChI_Aux};
use crate::source::base::util::inchi_free;
use crate::source_types::{
    _IS_FATAL, CANON_GLOBALS, COMP_ATOM_DATA, CT_OUT_OF_RAM, FILE, INCHI_BAS, INCHI_IOS_STRING,
    INCHI_IOSTREAM, INCHI_MODE, INCHI_NUM, INCHI_OUT_EMBED_REC, INCHI_OUT_PLAIN_TEXT,
    INCHI_OUT_PLAIN_TEXT_COMMENTS, INCHI_OUT_PRINT_OPTIONS, INCHI_OUT_TABBED_OUTPUT, INCHI_SORT,
    INPUT_PARMS, NORM_CANON_FLAGS, ORIG_ATOM_DATA, ORIG_STRUCT, OUT_TN, PINChI_Aux2, PINChI2,
    SourceConstPointer, SourceFormatArgument, SourceHeap, SourceHeapError, SourceMutPointer,
    SourceVaList, TAUT_NON, TAUT_NUM, TAUT_YES, TG_FLAG_DISCONNECT_COORD_DONE,
};

fn sort_inchi_rows(
    heap: &mut SourceHeap,
    rows: &mut [INCHI_SORT],
    tautomeric: bool,
) -> Result<(), SourceHeapError> {
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
    _composite_normalized_data: Option<
        &mut [[COMP_ATOM_DATA; TAUT_NUM as usize + 1]; INCHI_NUM as usize],
    >,
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

    let disconnected_coordinates =
        i32::from(taut_flags_done[0] & u64::from(TG_FLAG_DISCONNECT_COORD_DONE) != 0);
    for representation in 0..INCHI_NUM as usize {
        for tautomer in 0..TAUT_NUM as usize {
            taut_flags[representation] |= normalization_flags.bTautFlags[representation][tautomer];
            taut_flags_done[representation] |=
                normalization_flags.bTautFlagsDone[representation][tautomer];
        }
    }

    let mut max_components = component_counts.iter().copied().max().unwrap_or(0);
    if max_components <= 0 {
        max_components = 1;
    }
    let allocation_len =
        usize::try_from(max_components).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
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
                let count = usize::try_from(component_count)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
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
                    let row_pointer =
                        sort_rows[representation * TAUT_NUM as usize + tautomer_order];
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
                    let row_pointer =
                        sort_rows[representation * TAUT_NUM as usize + tautomer_order];
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
                    let mobile_pointer =
                        sort_rows[representation * TAUT_NUM as usize + TAUT_YES as usize];
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
                        let swap_tautomer_slots = !mobile_rows[component].pINChI[TAUT_NON as usize]
                            .is_null()
                            && mobile_rows[component].pINChI[TAUT_YES as usize].is_null();
                        for tautomer in 0..TAUT_NUM as usize {
                            let source_tautomer = if swap_tautomer_slots {
                                TAUT_NUM as usize - 1 - tautomer
                            } else {
                                tautomer
                            };
                            reordered_inchi[component][tautomer] =
                                mobile_rows[component].pINChI[source_tautomer];
                            reordered_aux[component][tautomer] =
                                mobile_rows[component].pINChI_Aux[source_tautomer];
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
                let base_options =
                    input_parameters.bINChIOutputOptions & !(INCHI_OUT_PRINT_OPTIONS as i32);
                let embed_reconnected =
                    input_parameters.bINChIOutputOptions & INCHI_OUT_EMBED_REC as i32;
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
                        let format = heap.allocate_model_storage(
                            b"\n==== %s ====\n\0"
                                .iter()
                                .map(|byte| *byte as i8)
                                .collect(),
                        )?;
                        let header = heap.allocate_model_storage(
                            b"InChI ANNOTATED CONTENTS\0"
                                .iter()
                                .map(|byte| *byte as i8)
                                .collect(),
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
        FreeINChIArrays(
            heap,
            p_inchi[representation],
            p_inchi_aux[representation],
            old_count,
        )?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{INCHI_IOS_TYPE_STRING, INCHI_OUT_NO_AUX_INFO, INChI, INChI_Aux};

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
        bytes[..end]
            .iter()
            .map(|byte| *byte as u8 as char)
            .collect()
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
        let aux_left = heap
            .allocate_model_storage(vec![INChI_Aux::default()])
            .unwrap();
        let aux_right = heap
            .allocate_model_storage(vec![INChI_Aux::default()])
            .unwrap();
        let aux_rows = heap
            .allocate_model_storage(vec![
                [aux_left, SourceMutPointer::null()],
                [aux_right, SourceMutPointer::null()],
            ])
            .unwrap();
        let mut output = stream(&mut heap);
        let mut log = stream(&mut heap);
        let mut composite =
            std::array::from_fn(|_| std::array::from_fn(|_| COMP_ATOM_DATA::default()));
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
        let mut failing_composite =
            std::array::from_fn(|_| std::array::from_fn(|_| COMP_ATOM_DATA::default()));
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
        let mut printing_composite =
            std::array::from_fn(|_| std::array::from_fn(|_| COMP_ATOM_DATA::default()));
        let flags = printing_heap.allocate_model_storage(vec![0_i32]).unwrap();
        let print_options = (INCHI_OUT_NO_AUX_INFO
            | INCHI_OUT_PLAIN_TEXT
            | INCHI_OUT_PLAIN_TEXT_COMMENTS
            | INCHI_OUT_TABBED_OUTPUT) as i32;
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
            concat!(
                "InChI=1//\n",
                "\n==== InChI ANNOTATED CONTENTS ====\n",
                "\nInChI=\n1\n"
            )
        );
        assert_eq!(text(&printing_heap, &printing_log), "");
    }

    #[test]
    fn source_port__runichi4__freeinchiarrays__line_1357() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            FreeINChIArrays(
                &mut heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                i32::MIN,
            ),
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
        assert_eq!(
            heap.slice(retained_inchi.as_const()).unwrap()[0].nRefCount,
            0
        );
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
        let inchi_rows = heap
            .allocate(vec![[inchi, SourceMutPointer::null()]])
            .unwrap();
        let aux = heap.allocate(vec![INChI_Aux::default()]).unwrap();
        let aux_rows = heap
            .allocate(vec![[aux, SourceMutPointer::null()]])
            .unwrap();
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
        assert_eq!(
            heap.slice(inchi.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(aux.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

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
}
