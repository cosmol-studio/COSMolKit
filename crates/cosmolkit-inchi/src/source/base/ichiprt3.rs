use crate::source::base::ichi_io::inchi_strbuf_printf;
use crate::source::base::ichimak2::MakeHillFormulaString;
use crate::source::base::ichiprt2::{
    Eql_INChI_Aux_Equ, Eql_INChI_Isotopic, Eql_INChI_Stereo, EqlOrigInfo, MakeCRVString, MakeCtString, MakeCtStringNew,
    MakeDelim, MakeEqStr, MakeEquString, MakeHString, MakeIsoAtomString, MakeIsoTautString, MakeMult, MakeStereoString,
    MakeTautString, bHasEquString, bHasOrigInfo,
};
use crate::source::base::util::{inchi_calloc, inchi_free};
use crate::source_types::{
    INCHI_IOS_STRING, INCHI_SORT, OUT_N1, OUT_NN, OUT_NT, OUT_T1, OUT_TN, SourceConstPointer, SourceFormatArgument,
    SourceHeap, SourceHeapError, SourceMutPointer, SourceVaList, TAUT_NON, TAUT_YES,
};

fn get_ii(heap: &SourceHeap, output_type: i32, sort: &INCHI_SORT) -> Result<Option<usize>, SourceHeapError> {
    // BEGIN INCHI ACTIVE MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.h:69 GET_II dependencies
    // INCHI✔️✔️: #define HAS_T(S)  (S->pINChI[TAUT_YES] && S->pINChI[TAUT_YES]->nNumberOfAtoms)
    // INCHI✔️✔️: #define HAS_N(S)  (S->pINChI[TAUT_NON] && S->pINChI[TAUT_NON]->nNumberOfAtoms)
    // INCHI✔️✔️:
    // INCHI✔️✔️: /* S->pINChI[TAUT_YES] has tautomeric info: */
    // INCHI✔️✔️: #define HAS_TT(S) (S->pINChI[TAUT_YES] && S->pINChI[TAUT_YES]->nNumberOfAtoms && S->pINChI[TAUT_YES]->lenTautomer>0)
    // INCHI✔️✔️: /* S->pINChI[TAUT_YES] has non-taitomeric info: */
    // INCHI✔️✔️: #define HAS_TN(S) (S->pINChI[TAUT_YES] && S->pINChI[TAUT_YES]->nNumberOfAtoms && !S->pINChI[TAUT_YES]->lenTautomer)
    // INCHI✔️✔️: /* S->pINChI[TAUT_NON] has non-tautomeric info: */
    // INCHI✔️✔️: #define HAS_NN(S) (S->pINChI[TAUT_NON] && S->pINChI[TAUT_NON]->nNumberOfAtoms && !S->pINChI[TAUT_NON]->lenTautomer)
    // INCHI✔️✔️: #define GET_II(M,S) ((M==OUT_N1)?              (HAS_TN(S)? TAUT_YES : HAS_NN(S)? TAUT_NON : -1): \
    // INCHI✔️✔️:                      (M==OUT_T1 || M==OUT_TN)? (HAS_T(S) ? TAUT_YES : HAS_N(S) ? TAUT_NON : -1): \
    // INCHI✔️✔️:                      (M==OUT_NN)?              (HAS_NN(S)? TAUT_NON : HAS_TN(S)? TAUT_YES : -1): \
    // INCHI✔️✔️:                      (M==OUT_NT)?              ((HAS_TT(S) && HAS_NN(S))       ? TAUT_NON : -1) : -1)
    // END INCHI ACTIVE MACRO: GET_II dependencies

    let non = sort.pINChI[TAUT_NON as usize];
    let taut = sort.pINChI[TAUT_YES as usize];
    let non_value = if non.is_null() {
        None
    } else {
        Some(
            heap.slice(non.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
    };
    let taut_value = if taut.is_null() {
        None
    } else {
        Some(
            heap.slice(taut.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
    };
    let has_n = non_value.is_some_and(|inchi| inchi.nNumberOfAtoms != 0);
    let has_t = taut_value.is_some_and(|inchi| inchi.nNumberOfAtoms != 0);
    let has_nn = non_value.is_some_and(|inchi| inchi.nNumberOfAtoms != 0 && inchi.lenTautomer == 0);
    let has_tn = taut_value.is_some_and(|inchi| inchi.nNumberOfAtoms != 0 && inchi.lenTautomer == 0);
    let has_tt = taut_value.is_some_and(|inchi| inchi.nNumberOfAtoms != 0 && inchi.lenTautomer > 0);
    Ok(if output_type == OUT_N1 as i32 {
        has_tn
            .then_some(TAUT_YES as usize)
            .or_else(|| has_nn.then_some(TAUT_NON as usize))
    } else if output_type == OUT_T1 as i32 || output_type == OUT_TN as i32 {
        has_t
            .then_some(TAUT_YES as usize)
            .or_else(|| has_n.then_some(TAUT_NON as usize))
    } else if output_type == OUT_NN as i32 {
        has_nn
            .then_some(TAUT_NON as usize)
            .or_else(|| has_tn.then_some(TAUT_YES as usize))
    } else if output_type == OUT_NT as i32 {
        (has_tt && has_nn).then_some(TAUT_NON as usize)
    } else {
        None
    })
}

fn selected_inchi(
    heap: &SourceHeap,
    sort_pointer: SourceConstPointer<INCHI_SORT>,
    output_type: i32,
) -> Result<SourceMutPointer<crate::source_types::INChI>, SourceHeapError> {
    let sort = heap
        .slice(sort_pointer)?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    Ok(get_ii(heap, output_type, sort)?.map_or(SourceMutPointer::null(), |index| sort.pINChI[index]))
}

fn c_string<'a>(heap: &'a SourceHeap, pointer: SourceConstPointer<i8>) -> Result<&'a [i8], SourceHeapError> {
    let bytes = heap.slice(pointer)?;
    let length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    Ok(&bytes[..length])
}

#[allow(non_snake_case)]
pub(crate) fn str_HillFormula(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    number_of_components: i32,
    use_multipliers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:60 str_HillFormula
    // INCHI✔️❌: int str_HillFormula( INCHI_SORT       *pINChISort,
    // INCHI✔️❌:                      INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                      int               *bOverflow,
    // INCHI✔️❌:                      int               bOutType,
    // INCHI✔️❌:                      int               num_components,
    // INCHI✔️❌:                      int               bUseMulipliers )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int          i, ii, nUsedLength0;
    // INCHI✔️❌:     INCHI_SORT   *is, *is0;
    // INCHI✔️❌:     INChI        *pINChI, *pINChI_Prev;
    // INCHI✔️❌:     int          mult, eq2prev, bNext;
    // INCHI✔️❌:
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!( is0 = pINChISort ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return nUsedLength0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     i = 0;
    // INCHI✔️❌:     pINChI_Prev = ( 0 <= ( ii = GET_II( bOutType, is0 ) ) ) ? is0->pINChI[ii] : NULL;
    // INCHI✔️❌:     mult = 0;
    // INCHI✔️❌:     bNext = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* For each connected component */
    // INCHI✔️❌:     for (i++; i <= num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         pINChI = ( i < num_components &&
    // INCHI✔️❌:             ( is = is0 + i, 0 <= ( ii = GET_II( bOutType, is ) ) ) )
    // INCHI✔️❌:             ? is->pINChI[ii]
    // INCHI✔️❌:             : NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:         eq2prev = bUseMulipliers &&
    // INCHI✔️❌:             pINChI &&
    // INCHI✔️❌:             pINChI_Prev &&
    // INCHI✔️❌:             pINChI->szHillFormula &&
    // INCHI✔️❌:             pINChI_Prev->szHillFormula &&
    // INCHI✔️❌:             pINChI->szHillFormula[0] &&
    // INCHI✔️❌:             !strcmp( pINChI_Prev->szHillFormula, pINChI->szHillFormula );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (eq2prev)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             mult++; /* mult = (number of non-empty equal items)-1 */
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int ok;
    // INCHI✔️❌:             if (bNext++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 MakeDelim( ".", strbuf, bOverflow );
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             ok = pINChI_Prev &&
    // INCHI✔️❌:                 pINChI_Prev->szHillFormula &&
    // INCHI✔️❌:                 pINChI_Prev->szHillFormula[0];
    // INCHI✔️❌:
    // INCHI✔️❌:             if (ok)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 MakeMult( mult + 1, "", strbuf, 0, bOverflow );
    // INCHI✔️❌:                 MakeHillFormulaString( pINChI_Prev->szHillFormula,
    // INCHI✔️❌:                     strbuf,
    // INCHI✔️❌:                     bOverflow );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         pINChI_Prev = pINChI;
    // INCHI✔️❌:         mult = 0; /* we do not know whether the item is empty */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( strbuf->nUsedLength - nUsedLength0 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: str_HillFormula

    let initial_length = string_buffer.nUsedLength;
    if sorted_inchi.is_null() {
        return Ok(initial_length);
    }
    let mut previous = selected_inchi(heap, sorted_inchi.as_const(), output_type)?;
    let dot = heap.allocate_model_storage(vec![b'.' as i8, 0])?;
    let empty = heap.allocate_model_storage(vec![0_i8])?;
    let mut multiplier = 0_i32;
    let mut next = 0_i32;
    let mut index = 1_i32;
    while index <= number_of_components {
        let current = if index < number_of_components {
            selected_inchi(heap, sorted_inchi.as_const().offset(i64::from(index))?, output_type)?
        } else {
            SourceMutPointer::null()
        };
        let equal_to_previous = if use_multipliers != 0 && !current.is_null() && !previous.is_null() {
            let current_inchi = heap
                .slice(current.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let previous_inchi = heap
                .slice(previous.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            !current_inchi.szHillFormula.is_null()
                && !previous_inchi.szHillFormula.is_null()
                && !c_string(heap, current_inchi.szHillFormula.as_const())?.is_empty()
                && c_string(heap, current_inchi.szHillFormula.as_const())?
                    == c_string(heap, previous_inchi.szHillFormula.as_const())?
        } else {
            false
        };
        if equal_to_previous {
            multiplier = multiplier.wrapping_add(1);
            index = index.wrapping_add(1);
            continue;
        }
        if next != 0 {
            MakeDelim(heap, dot.as_const(), string_buffer, overflow)?;
        }
        next = next.wrapping_add(1);
        if !previous.is_null() {
            let formula = heap
                .slice(previous.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .szHillFormula;
            if !formula.is_null() && !c_string(heap, formula.as_const())?.is_empty() {
                MakeMult(
                    heap,
                    multiplier.wrapping_add(1),
                    empty.as_const(),
                    string_buffer,
                    0,
                    overflow,
                )?;
                MakeHillFormulaString(heap, formula.as_const(), string_buffer, overflow)?;
            }
        }
        previous = current;
        multiplier = 0;
        index = index.wrapping_add(1);
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_HillFormula2(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    sorted_inchi2: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    number_of_components: i32,
    use_multipliers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:137 str_HillFormula2
    // INCHI✔️❌: int str_HillFormula2( INCHI_SORT       *pINChISort,     /* non-taut */
    // INCHI✔️❌:                       INCHI_SORT       *pINChISort2,    /* taut     */
    // INCHI✔️❌:                       INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                       int              *bOverflow,
    // INCHI✔️❌:                       int              bOutType,
    // INCHI✔️❌:                       int              num_components,
    // INCHI✔️❌:                       int              bUseMulipliers )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int          i, ii, ii2, nUsedLength0;
    // INCHI✔️❌:     INCHI_SORT   *is, *is2, *is0, *is20;
    // INCHI✔️❌:     INChI        *pINChI, *pINChI_Prev, *pINChI_Taut, *pINChI_Taut_Prev;
    // INCHI✔️❌:     int          mult, eq2prev, bNext, bEqToTaut;
    // INCHI✔️❌:
    // INCHI✔️❌:     is = NULL;
    // INCHI✔️❌:     is2 = NULL;
    // INCHI✔️❌:     is0 = pINChISort;
    // INCHI✔️❌:     is20 = pINChISort2;
    // INCHI✔️❌:     i = 0;
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     pINChI_Prev = ( 0 <= ( ii = GET_II( bOutType, is0 ) ) ) ? is0->pINChI[ii] : NULL;
    // INCHI✔️❌:     pINChI_Taut_Prev = ( 0 <= ( ii2 = GET_II( OUT_T1, is20 ) ) ) ? is20->pINChI[ii2] : NULL;
    // INCHI✔️❌:     mult = 0;
    // INCHI✔️❌:     bNext = 0;
    // INCHI✔️❌:     bEqToTaut = 1;
    // INCHI✔️❌:     bEqToTaut = bEqToTaut                 &&
    // INCHI✔️❌:         pINChI_Prev                       &&
    // INCHI✔️❌:         pINChI_Taut_Prev &&
    // INCHI✔️❌:         !pINChI_Taut_Prev->bDeleted       &&
    // INCHI✔️❌:         pINChI_Prev->szHillFormula        &&
    // INCHI✔️❌:         pINChI_Taut_Prev->szHillFormula   &&
    // INCHI✔️❌:         !strcmp( pINChI_Prev->szHillFormula, pINChI_Taut_Prev->szHillFormula );
    // INCHI✔️❌:
    // INCHI✔️❌:     /* For each connected component    */
    // INCHI✔️❌:     for (i++; i <= num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         pINChI = ( i < num_components && ( is = is0 + i, 0 <= ( ii = GET_II( bOutType, is ) ) ) ) ? is->pINChI[ii] : NULL;
    // INCHI✔️❌:         pINChI_Taut = ( i < num_components && ( is2 = is20 + i, 0 <= ( ii2 = GET_II( OUT_T1, is2 ) ) ) ) ? is2->pINChI[ii2] : NULL;
    // INCHI✔️❌:         if (bEqToTaut && ( pINChI || pINChI_Taut ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             bEqToTaut = pINChI && pINChI_Taut && !pINChI_Taut->bDeleted &&
    // INCHI✔️❌:                 pINChI->szHillFormula && pINChI_Taut->szHillFormula     &&
    // INCHI✔️❌:                 !strcmp( pINChI->szHillFormula, pINChI_Taut->szHillFormula );
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         eq2prev = bUseMulipliers        &&
    // INCHI✔️❌:             pINChI && pINChI_Prev       &&
    // INCHI✔️❌:             pINChI->szHillFormula       &&
    // INCHI✔️❌:             pINChI_Prev->szHillFormula  &&
    // INCHI✔️❌:             pINChI->szHillFormula[0]    &&
    // INCHI✔️❌:             !strcmp( pINChI_Prev->szHillFormula, pINChI->szHillFormula );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (eq2prev)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             mult++; /* mult = (number of non-empty equal items)-1 */
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (bNext++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 MakeDelim( ".", strbuf, bOverflow );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (pINChI_Prev && pINChI_Prev->szHillFormula && pINChI_Prev->szHillFormula[0])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 MakeMult( mult + 1, "", strbuf, 0, bOverflow );
    // INCHI✔️❌:                 MakeHillFormulaString( pINChI_Prev->szHillFormula, strbuf, bOverflow );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         pINChI_Prev = pINChI;
    // INCHI✔️❌:         mult = 0; /* we do not know whether the item is empty */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (bEqToTaut)
    // INCHI✔️❌:         strbuf->nUsedLength = nUsedLength0;
    // INCHI✔️❌:
    // INCHI✔️❌:     {
    // INCHI✔️❌:         strbuf->pStr[strbuf->nUsedLength] = '\0';
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( strbuf->nUsedLength - nUsedLength0 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: str_HillFormula2

    fn equal_formula(
        heap: &SourceHeap,
        fixed: SourceMutPointer<crate::source_types::INChI>,
        taut: SourceMutPointer<crate::source_types::INChI>,
    ) -> Result<bool, SourceHeapError> {
        if fixed.is_null() || taut.is_null() {
            return Ok(false);
        }
        let fixed = heap
            .slice(fixed.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let taut = heap
            .slice(taut.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        Ok(taut.bDeleted == 0
            && !fixed.szHillFormula.is_null()
            && !taut.szHillFormula.is_null()
            && c_string(heap, fixed.szHillFormula.as_const())? == c_string(heap, taut.szHillFormula.as_const())?)
    }

    let initial_length = string_buffer.nUsedLength;
    let mut previous = selected_inchi(heap, sorted_inchi.as_const(), output_type)?;
    let mut taut_previous = selected_inchi(heap, sorted_inchi2.as_const(), OUT_T1 as i32)?;
    let mut equal_to_taut = equal_formula(heap, previous, taut_previous)?;
    let dot = heap.allocate_model_storage(vec![b'.' as i8, 0])?;
    let empty = heap.allocate_model_storage(vec![0_i8])?;
    let mut multiplier = 0_i32;
    let mut next = 0_i32;
    let mut index = 1_i32;
    while index <= number_of_components {
        let current = if index < number_of_components {
            selected_inchi(heap, sorted_inchi.as_const().offset(i64::from(index))?, output_type)?
        } else {
            SourceMutPointer::null()
        };
        let taut_current = if index < number_of_components {
            selected_inchi(heap, sorted_inchi2.as_const().offset(i64::from(index))?, OUT_T1 as i32)?
        } else {
            SourceMutPointer::null()
        };
        if equal_to_taut && (!current.is_null() || !taut_current.is_null()) {
            equal_to_taut = equal_formula(heap, current, taut_current)?;
        }
        let equal_to_previous = if use_multipliers != 0 && !current.is_null() && !previous.is_null() {
            let current_value = heap
                .slice(current.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let previous_value = heap
                .slice(previous.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            !current_value.szHillFormula.is_null()
                && !previous_value.szHillFormula.is_null()
                && !c_string(heap, current_value.szHillFormula.as_const())?.is_empty()
                && c_string(heap, current_value.szHillFormula.as_const())?
                    == c_string(heap, previous_value.szHillFormula.as_const())?
        } else {
            false
        };
        if equal_to_previous {
            multiplier = multiplier.wrapping_add(1);
            index = index.wrapping_add(1);
            continue;
        }
        if next != 0 {
            MakeDelim(heap, dot.as_const(), string_buffer, overflow)?;
        }
        next = next.wrapping_add(1);
        if !previous.is_null() {
            let formula = heap
                .slice(previous.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .szHillFormula;
            if !formula.is_null() && !c_string(heap, formula.as_const())?.is_empty() {
                MakeMult(
                    heap,
                    multiplier.wrapping_add(1),
                    empty.as_const(),
                    string_buffer,
                    0,
                    overflow,
                )?;
                MakeHillFormulaString(heap, formula.as_const(), string_buffer, overflow)?;
            }
        }
        previous = current;
        taut_previous = taut_current;
        let _ = taut_previous;
        multiplier = 0;
        index = index.wrapping_add(1);
    }
    if equal_to_taut {
        string_buffer.nUsedLength = initial_length;
    }
    let used = usize::try_from(string_buffer.nUsedLength).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    *heap
        .slice_mut(string_buffer.pStr)?
        .get_mut(used)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case)]
pub(crate) fn str_Connections(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    atom_mode: i32,
    number_of_components: i32,
    use_multipliers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:224 str_Connections
    // INCHI✔️❌: int str_Connections( CANON_GLOBALS    *pCG,
    // INCHI✔️❌:                      INCHI_SORT       *pINChISort,
    // INCHI✔️❌:                      INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                      int              *bOverflow,
    // INCHI✔️❌:                      int              bOutType,
    // INCHI✔️❌:                      int              ATOM_MODE,
    // INCHI✔️❌:                      int              num_components,
    // INCHI✔️❌:                      int              bUseMulipliers )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int          i, ii, nUsedLength0;
    // INCHI✔️❌:     INCHI_SORT   *is, *is0;
    // INCHI✔️❌:     INChI        *pINChI, *pINChI_Prev;
    // INCHI✔️❌:     int          mult, eq2prev, bNext, nNumEmpty;
    // INCHI✔️❌:
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:     if (!( is0 = pINChISort ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return nUsedLength0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     i = 0;
    // INCHI✔️❌:     pINChI_Prev = ( 0 <= ( ii = GET_II( bOutType, is0 ) ) ) ? is0->pINChI[ii] : NULL;
    // INCHI✔️❌:     is = NULL;
    // INCHI✔️❌:     mult = 0;
    // INCHI✔️❌:     bNext = 0;
    // INCHI✔️❌:     nNumEmpty = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* For each connected component...    */
    // INCHI✔️❌:     for (i++; i <= num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         pINChI = ( i < num_components && ( is = is0 + i, 0 <= ( ii = GET_II( bOutType, is ) ) ) ) ? is->pINChI[ii] : NULL;
    // INCHI✔️❌:         eq2prev = bUseMulipliers                                &&
    // INCHI✔️❌:             pINChI && pINChI_Prev && pINChI->lenConnTable > 1   &&
    // INCHI✔️❌:             pINChI_Prev->lenConnTable == pINChI->lenConnTable   &&
    // INCHI✔️❌:             !memcmp( pINChI_Prev->nConnTable, pINChI->nConnTable,
    // INCHI✔️❌:                 pINChI_Prev->lenConnTable * sizeof( pINChI->nConnTable[0] ) );
    // INCHI✔️❌:         if (eq2prev)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             mult++; /* mult = (number of non-empty equal items)-1 */
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (pINChI_Prev)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (bNext++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (pINChI_Prev && pINChI_Prev->lenConnTable > 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
    // INCHI✔️❌:                     MakeCtStringNew( pCG, pINChI_Prev->nConnTable,
    // INCHI✔️❌:                         pINChI_Prev->lenConnTable,
    // INCHI✔️❌:                         0, NULL,
    // INCHI✔️❌:                         pINChI_Prev->nNumberOfAtoms,
    // INCHI✔️❌:                         strbuf, ATOM_MODE, bOverflow );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nNumEmpty++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         pINChI_Prev = pINChI;
    // INCHI✔️❌:         mult = 0; /* we do not know whether the item is empty */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (nNumEmpty == num_components && strbuf->nUsedLength > nUsedLength0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         strbuf->nUsedLength = nUsedLength0;
    // INCHI✔️❌:         strbuf->pStr[strbuf->nUsedLength] = '\0';
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( nNumEmpty == num_components && tot_len > tot_len_inp ) {
    // INCHI✔️❌:     tot_len = tot_len_inp;
    // INCHI✔️❌:     strbuf->pStr[tot_len] = '\0';
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( strbuf->nUsedLength - nUsedLength0 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: str_Connections

    let initial_length = string_buffer.nUsedLength;
    if sorted_inchi.is_null() {
        return Ok(initial_length);
    }
    let mut previous = selected_inchi(heap, sorted_inchi.as_const(), output_type)?;
    let component_delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let star = heap.allocate_model_storage(vec![b'*' as i8, 0])?;
    let mut multiplier = 0_i32;
    let mut next = 0_i32;
    let mut number_empty = 0_i32;
    let mut index = 1_i32;
    while index <= number_of_components {
        let current = if index < number_of_components {
            selected_inchi(heap, sorted_inchi.as_const().offset(i64::from(index))?, output_type)?
        } else {
            SourceMutPointer::null()
        };
        let equal_to_previous = if use_multipliers != 0 && !current.is_null() && !previous.is_null() {
            let current_inchi = heap
                .slice(current.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let previous_inchi = heap
                .slice(previous.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if current_inchi.lenConnTable > 1 && previous_inchi.lenConnTable == current_inchi.lenConnTable {
                let length =
                    usize::try_from(current_inchi.lenConnTable).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                heap.slice(current_inchi.nConnTable.as_const())?
                    .get(..length)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    == heap
                        .slice(previous_inchi.nConnTable.as_const())?
                        .get(..length)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
            } else {
                false
            }
        } else {
            false
        };
        if equal_to_previous {
            multiplier = multiplier.wrapping_add(1);
            index = index.wrapping_add(1);
            continue;
        }
        if !previous.is_null() {
            if next != 0 {
                MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
            }
            next = next.wrapping_add(1);
            let previous_inchi = heap
                .slice(previous.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if previous_inchi.lenConnTable > 1 {
                let connection_table = previous_inchi.nConnTable;
                let connection_length = previous_inchi.lenConnTable;
                let atom_count = previous_inchi.nNumberOfAtoms;
                MakeMult(
                    heap,
                    multiplier.wrapping_add(1),
                    star.as_const(),
                    string_buffer,
                    0,
                    overflow,
                )?;
                MakeCtStringNew(
                    heap,
                    connection_table,
                    connection_length,
                    0,
                    SourceConstPointer::null(),
                    atom_count,
                    string_buffer,
                    atom_mode,
                    overflow,
                )?;
            } else {
                number_empty = number_empty.wrapping_add(1);
            }
        }
        previous = current;
        multiplier = 0;
        index = index.wrapping_add(1);
    }
    if number_empty == number_of_components && string_buffer.nUsedLength > initial_length {
        string_buffer.nUsedLength = initial_length;
        let used = usize::try_from(initial_length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        heap.slice_mut(string_buffer.pStr)?[used] = 0;
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case)]
pub(crate) fn str_H_atoms(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    atom_mode: i32,
    taut_mode: i32,
    number_of_components: i32,
    use_multipliers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:309 str_H_atoms
    // INCHI✔️❌: int str_H_atoms( INCHI_SORT       *pINChISort,
    // INCHI✔️❌:                  INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                  int               *bOverflow,
    // INCHI✔️❌:                  int               bOutType,
    // INCHI✔️❌:                  int               ATOM_MODE,
    // INCHI✔️❌:                  int               TAUT_MODE,
    // INCHI✔️❌:                  int               num_components,
    // INCHI✔️❌:                  int               bUseMulipliers )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int          i, j, ii, len_H, nUsedLength0;
    // INCHI✔️❌:     INCHI_SORT   *is, *is0;
    // INCHI✔️❌:     INChI        *pINChI, *pINChI_Prev;
    // INCHI✔️❌:     int          mult, eq2prev, bNext, bNotEmpty, nNumEmpty;
    // INCHI✔️❌:
    // INCHI✔️❌:     nNumEmpty = 0;
    // INCHI✔️❌:     is0 = pINChISort;
    // INCHI✔️❌:     is = NULL;
    // INCHI✔️❌:     i = 0;
    // INCHI✔️❌:     pINChI_Prev = ( 0 <= ( ii = GET_II( bOutType, is0 ) ) ) ? is0->pINChI[ii] : NULL;
    // INCHI✔️❌:     mult = 0;
    // INCHI✔️❌:     bNext = 0;
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* For each connected component...    */
    // INCHI✔️❌:     for (i++; i <= num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         pINChI = ( i < num_components && ( is = is0 + i, 0 <= ( ii = GET_II( bOutType, is ) ) ) ) ? is->pINChI[ii] : NULL;
    // INCHI✔️❌:         /*========== compare to previous ============*/
    // INCHI✔️❌:         eq2prev = bUseMulipliers &&
    // INCHI✔️❌:             pINChI && pINChI_Prev && ( pINChI->nNumberOfAtoms > 0 || pINChI->lenTautomer > 1 ) &&
    // INCHI✔️❌:             pINChI_Prev->nNumberOfAtoms == pINChI->nNumberOfAtoms &&
    // INCHI✔️❌:             ( !pINChI_Prev->nNumberOfAtoms || !memcmp( pINChI_Prev->nNum_H, pINChI->nNum_H,
    // INCHI✔️❌:                 pINChI_Prev->nNumberOfAtoms * sizeof( pINChI->nNum_H[0] ) ) ) &&
    // INCHI✔️❌:             !CompareTautNonIsoPartOfINChI( pINChI_Prev, pINChI );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (eq2prev && pINChI_Prev->lenTautomer <= 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* make sure it is not empty */
    // INCHI✔️❌:             eq2prev = 0;
    // INCHI✔️❌:             for (j = 0; j < pINChI_Prev->nNumberOfAtoms; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (pINChI_Prev->nNum_H[j])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     eq2prev = 1;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (eq2prev)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             mult++; /* mult = (number of non-empty equal items)-1 */
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (pINChI_Prev)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* delimiter */
    // INCHI✔️❌:                 if (bNext++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* verify non-empty */
    // INCHI✔️❌:                 bNotEmpty = 0;
    // INCHI✔️❌:                 if (pINChI_Prev)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     bNotEmpty = ( pINChI_Prev->lenTautomer > 1 );
    // INCHI✔️❌:                     if (!bNotEmpty)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         for (j = 0; j < pINChI_Prev->nNumberOfAtoms; j++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (pINChI_Prev->nNum_H[j])
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 bNotEmpty = 1;
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (bNotEmpty)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
    // INCHI✔️❌:                     /* H-atoms */
    // INCHI✔️❌:                     len_H = MakeHString( 0, pINChI_Prev->nNum_H, pINChI_Prev->nNumberOfAtoms,
    // INCHI✔️❌:                         strbuf, ATOM_MODE, bOverflow );
    // INCHI✔️❌:                     /*  tautomeric groups */
    // INCHI✔️❌:                     MakeTautString( pINChI_Prev->nTautomer, pINChI_Prev->lenTautomer, ( 0 != len_H ),
    // INCHI✔️❌:                         strbuf, TAUT_MODE, bOverflow );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nNumEmpty++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         pINChI_Prev = pINChI;
    // INCHI✔️❌:         mult = 0; /* we do not know whether the item is empty */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (nNumEmpty == num_components && strbuf->nUsedLength > nUsedLength0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         strbuf->nUsedLength = nUsedLength0;
    // INCHI✔️❌:         strbuf->pStr[nUsedLength0] = '\0';
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( nNumEmpty == num_components && tot_len > tot_len_inp ) {
    // INCHI✔️❌:     tot_len = tot_len_inp;
    // INCHI✔️❌:     strbuf->pStr[tot_len] = '\0';
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( strbuf->nUsedLength - nUsedLength0 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: str_H_atoms

    let initial_length = string_buffer.nUsedLength;
    let mut previous = selected_inchi(heap, sorted_inchi.as_const(), output_type)?;
    let component_delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let star = heap.allocate_model_storage(vec![b'*' as i8, 0])?;
    let mut multiplier = 0_i32;
    let mut next = 0_i32;
    let mut empty_count = 0_i32;
    let mut index = 1_i32;
    while index <= number_of_components {
        let current = if index < number_of_components {
            selected_inchi(heap, sorted_inchi.as_const().offset(i64::from(index))?, output_type)?
        } else {
            SourceMutPointer::null()
        };
        let mut equal = false;
        if use_multipliers != 0 && !current.is_null() && !previous.is_null() {
            let current_value = heap
                .slice(current.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let previous_value = heap
                .slice(previous.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            if (current_value.nNumberOfAtoms > 0 || current_value.lenTautomer > 1)
                && previous_value.nNumberOfAtoms == current_value.nNumberOfAtoms
            {
                let hydrogens_equal = if previous_value.nNumberOfAtoms == 0 {
                    true
                } else {
                    let count = usize::try_from(previous_value.nNumberOfAtoms)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    heap.slice(previous_value.nNum_H.as_const())?
                        .get(..count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        == heap
                            .slice(current_value.nNum_H.as_const())?
                            .get(..count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                };
                equal = hydrogens_equal
                    && crate::source::base::ichimake::CompareTautNonIsoPartOfINChI(
                        heap,
                        &previous_value,
                        &current_value,
                    )? == 0;
            }
            if equal && previous_value.lenTautomer <= 1 {
                let count = usize::try_from(previous_value.nNumberOfAtoms)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                equal = heap
                    .slice(previous_value.nNum_H.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .iter()
                    .any(|value| *value != 0);
            }
        }
        if equal {
            multiplier = multiplier.wrapping_add(1);
            index = index.wrapping_add(1);
            continue;
        }
        if !previous.is_null() {
            if next != 0 {
                MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
            }
            next = next.wrapping_add(1);
            let value = heap
                .slice(previous.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let count = usize::try_from(value.nNumberOfAtoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let nonempty = if value.lenTautomer > 1 {
                true
            } else if count == 0 {
                false
            } else {
                heap.slice(value.nNum_H.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .iter()
                    .any(|hydrogen| *hydrogen != 0)
            };
            if nonempty {
                MakeMult(
                    heap,
                    multiplier.wrapping_add(1),
                    star.as_const(),
                    string_buffer,
                    0,
                    overflow,
                )?;
                let h_length = MakeHString(
                    heap,
                    0,
                    value.nNum_H.as_const(),
                    value.nNumberOfAtoms,
                    string_buffer,
                    atom_mode,
                    overflow,
                )?;
                MakeTautString(
                    heap,
                    value.nTautomer,
                    value.lenTautomer,
                    i32::from(h_length != 0),
                    string_buffer,
                    taut_mode,
                    overflow,
                )?;
            } else {
                empty_count = empty_count.wrapping_add(1);
            }
        }
        previous = current;
        multiplier = 0;
        index = index.wrapping_add(1);
    }
    if empty_count == number_of_components && string_buffer.nUsedLength > initial_length {
        string_buffer.nUsedLength = initial_length;
        heap.slice_mut(string_buffer.pStr)?
            [usize::try_from(initial_length).map_err(|_| SourceHeapError::PointerOutOfBounds)?] = 0;
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

fn append_signed_decimal(
    heap: &mut SourceHeap,
    string_buffer: &mut INCHI_IOS_STRING,
    value: i32,
) -> Result<(), SourceHeapError> {
    let format = heap.allocate_model_storage(vec![b'%' as i8, b'+' as i8, b'd' as i8, 0])?;
    match inchi_strbuf_printf(
        heap,
        Some(string_buffer),
        format.as_const(),
        &SourceVaList {
            arguments: vec![SourceFormatArgument::Signed(i64::from(value))],
            ..SourceVaList::default()
        },
    ) {
        Ok(_) | Err(SourceHeapError::AllocationFailed) => Ok(()),
        Err(error) => Err(error),
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_FixedH_atoms(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    atom_mode: i32,
    number_of_components: i32,
    use_multipliers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:618 str_FixedH_atoms
    // INCHI✔️❌: int str_FixedH_atoms( INCHI_SORT       *pINChISort,
    // INCHI✔️❌:                       INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                       int              *bOverflow,
    // INCHI✔️❌:                       int              bOutType,
    // INCHI✔️❌:                       int              ATOM_MODE,
    // INCHI✔️❌:                       int              num_components,
    // INCHI✔️❌:                       int              bUseMulipliers )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int          i, j, ii, nNumEmpty, nUsedLength0;
    // INCHI✔️❌:     INCHI_SORT   *is, *is0;
    // INCHI✔️❌:     INChI        *pINChI, *pINChI_Prev;
    // INCHI✔️❌:     int          mult, eq2prev, bNext, bNotEmpty;
    // INCHI✔️❌:
    // INCHI✔️❌:     is = NULL;
    // INCHI✔️❌:     is0 = pINChISort;
    // INCHI✔️❌:     i = 0;
    // INCHI✔️❌:     pINChI_Prev = ( 0 <= ( ii = GET_II( bOutType, is0 ) ) ) ? is0->pINChI[ii] : NULL;
    // INCHI✔️❌:     mult = 0;
    // INCHI✔️❌:     bNext = 0;
    // INCHI✔️❌:     nNumEmpty = 0;
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* For each connected component...    */
    // INCHI✔️❌:     for (i++; i <= num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* only non-tautomeric representation of tautomeric */
    // INCHI✔️❌:         pINChI = ( i < num_components && ( is = is0 + i, 0 <= ( ii = GET_II( bOutType, is ) ) ) ) ? is->pINChI[ii] : NULL;
    // INCHI✔️❌:         /*================ compare fixed H to previous =====================*/
    // INCHI✔️❌:         eq2prev = bUseMulipliers &&
    // INCHI✔️❌:             pINChI && pINChI_Prev && pINChI->nNumberOfAtoms > 0 &&
    // INCHI✔️❌:             pINChI_Prev->nNumberOfAtoms == pINChI->nNumberOfAtoms &&
    // INCHI✔️❌:             !memcmp( pINChI_Prev->nNum_H_fixed, pINChI->nNum_H_fixed,
    // INCHI✔️❌:                 pINChI_Prev->nNumberOfAtoms * sizeof( pINChI->nNum_H_fixed[0] ) );
    // INCHI✔️❌:         if (eq2prev)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* make sure it is not empty */
    // INCHI✔️❌:             eq2prev = 0;
    // INCHI✔️❌:             for (j = 0; j < pINChI_Prev->nNumberOfAtoms; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (pINChI_Prev->nNum_H_fixed[j])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     eq2prev = 1;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (eq2prev)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             mult++; /* mult = (number of non-empty equal items)-1 */
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* print pINChI_Prev */
    // INCHI✔️❌:             /* delimiter */
    // INCHI✔️❌:             if (bNext++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (pINChI_Prev)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* verify it is not empty */
    // INCHI✔️❌:                 bNotEmpty = 0;
    // INCHI✔️❌:                 for (j = 0; j < pINChI_Prev->nNumberOfAtoms; j++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (pINChI_Prev->nNum_H_fixed[j])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         bNotEmpty = 1;
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (bNotEmpty)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
    // INCHI✔️❌:                     /* H-atoms-fixed */
    // INCHI✔️❌:                     MakeHString( 0, pINChI_Prev->nNum_H_fixed,
    // INCHI✔️❌:                         pINChI_Prev->nNumberOfAtoms,
    // INCHI✔️❌:                         strbuf, ATOM_MODE, bOverflow );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nNumEmpty++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         pINChI_Prev = pINChI;
    // INCHI✔️❌:         mult = 0; /* we do not know whether the item is empty */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (nNumEmpty == num_components && strbuf->nUsedLength > nUsedLength0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         strbuf->nUsedLength = nUsedLength0;
    // INCHI✔️❌:         strbuf->pStr[nUsedLength0] = '\0';
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( nNumEmpty == num_components && tot_len > tot_len_inp ) {
    // INCHI✔️❌:     tot_len = tot_len_inp;
    // INCHI✔️❌:     strbuf->pStr[tot_len] = '\0';
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( strbuf->nUsedLength - nUsedLength0 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: str_FixedH_atoms

    fn fixed_h<'a>(
        heap: &'a SourceHeap,
        pointer: SourceMutPointer<crate::source_types::INChI>,
    ) -> Result<&'a [i8], SourceHeapError> {
        let value = heap
            .slice(pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let count = usize::try_from(value.nNumberOfAtoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        heap.slice(value.nNum_H_fixed.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }
    let initial_length = string_buffer.nUsedLength;
    let mut previous = selected_inchi(heap, sorted_inchi.as_const(), output_type)?;
    let delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let star = heap.allocate_model_storage(vec![b'*' as i8, 0])?;
    let mut multiplier = 0_i32;
    let mut next = 0_i32;
    let mut empty_count = 0_i32;
    let mut index = 1_i32;
    while index <= number_of_components {
        let current = if index < number_of_components {
            selected_inchi(heap, sorted_inchi.as_const().offset(i64::from(index))?, output_type)?
        } else {
            SourceMutPointer::null()
        };
        let equal = if use_multipliers != 0 && !current.is_null() && !previous.is_null() {
            let current_value = heap
                .slice(current.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let previous_value = heap
                .slice(previous.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            current_value.nNumberOfAtoms > 0
                && previous_value.nNumberOfAtoms == current_value.nNumberOfAtoms
                && fixed_h(heap, previous)? == fixed_h(heap, current)?
                && fixed_h(heap, previous)?.iter().any(|value| *value != 0)
        } else {
            false
        };
        if equal {
            multiplier = multiplier.wrapping_add(1);
            index = index.wrapping_add(1);
            continue;
        }
        if next != 0 {
            MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
        }
        next = next.wrapping_add(1);
        if !previous.is_null() {
            let hydrogens = fixed_h(heap, previous)?;
            if hydrogens.iter().any(|value| *value != 0) {
                let value = heap
                    .slice(previous.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let pointer = value.nNum_H_fixed;
                let atom_count = value.nNumberOfAtoms;
                MakeMult(
                    heap,
                    multiplier.wrapping_add(1),
                    star.as_const(),
                    string_buffer,
                    0,
                    overflow,
                )?;
                MakeHString(
                    heap,
                    0,
                    pointer.as_const(),
                    atom_count,
                    string_buffer,
                    atom_mode,
                    overflow,
                )?;
            } else {
                empty_count = empty_count.wrapping_add(1);
            }
        }
        previous = current;
        multiplier = 0;
        index = index.wrapping_add(1);
    }
    if empty_count == number_of_components && string_buffer.nUsedLength > initial_length {
        string_buffer.nUsedLength = initial_length;
        let used = usize::try_from(initial_length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        heap.slice_mut(string_buffer.pStr)?[used] = 0;
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_Charge2(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    sorted_inchi2: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    number_of_components: i32,
    second_non_taut_pass: i32,
    omit_repetitions: i32,
    use_multipliers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:427 str_Charge2
    // INCHI✔️❌: int str_Charge2( INCHI_SORT       *pINChISort,
    // INCHI✔️❌:                  INCHI_SORT       *pINChISort2,
    // INCHI✔️❌:                  INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                  int              *bOverflow,
    // INCHI✔️❌:                  int              bOutType,
    // INCHI✔️❌:                  int              num_components,
    // INCHI✔️❌:                  int              bSecondNonTautPass,
    // INCHI✔️❌:                  int              bOmitRepetitions,
    // INCHI✔️❌:                  int              bUseMulipliers )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int          i, ii, ii2, nUsedLength0;
    // INCHI✔️❌:     INCHI_SORT   *is, *is2, *is0, *is20;
    // INCHI✔️❌:     INChI        *pINChI, *pINChI_Prev, *pINChI_Taut, *pINChI_Taut_Prev;
    // INCHI✔️❌:     int         nTotalCharge, nTotalCharge_Prev, nTotalCharge_Taut; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     int          mult, eq2prev, eq2taut, eq2tautPrev, bNext;
    // INCHI✔️❌:     const char  *pPrevEquStr, *pCurrEquStr;
    // INCHI✔️❌:     int         multPrevEquStr;
    // INCHI✔️❌:     pINChI_Taut = NULL;
    // INCHI✔️❌:     pINChI_Prev = NULL;
    // INCHI✔️❌:     pINChI_Taut_Prev = NULL;
    // INCHI✔️❌:     mult = 0;
    // INCHI✔️❌:     bNext = 0;
    // INCHI✔️❌:     is = NULL;
    // INCHI✔️❌:     is2 = NULL;
    // INCHI✔️❌:     is0 = pINChISort;
    // INCHI✔️❌:     is20 = bSecondNonTautPass ? pINChISort2 : NULL;
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     eq2tautPrev = 1; /* pINChI_Prev (previous pINChI) does not exist */
    // INCHI✔️❌:     pPrevEquStr = NULL; /*, *pCurrEquStr;*/
    // INCHI✔️❌:     multPrevEquStr = 0;
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* For each connected component...    */
    // INCHI✔️❌:     for (i = 0; i <= num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         /* 1st (taut) pass: bOutType=OUT_TN  ; 2nd (non-taut pass) bOutType=OUT_NT */
    // INCHI✔️❌:         pINChI = ( i < num_components && ( is = is0 + i, 0 <= ( ii = GET_II( bOutType, is ) ) ) ) ? is->pINChI[ii] : NULL;
    // INCHI✔️❌:         /*================ compare sp3 to previous =====================*/
    // INCHI✔️❌:         if (bSecondNonTautPass)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* component that was output on the 1st pass */
    // INCHI✔️❌:             pINChI_Taut = ( i < num_components && ( is2 = is20 + i, 0 <= ( ii2 = GET_II( OUT_T1, is2 ) ) ) ) ? is2->pINChI[ii2] : NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /*========= if bSecondNonTautPass then compare non-iso non-taut stereo to non-iso taut ========*/
    // INCHI✔️❌:         eq2taut = 0;
    // INCHI✔️❌:         if (!eq2taut && bSecondNonTautPass && bOmitRepetitions)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             eq2taut = pINChI && pINChI_Taut && !pINChI_Taut->bDeleted &&
    // INCHI✔️❌:                 ( nTotalCharge = pINChI->nTotalCharge ) && ( nTotalCharge_Taut = pINChI_Taut->nTotalCharge ) &&
    // INCHI✔️❌:                 nTotalCharge == nTotalCharge_Taut;
    // INCHI✔️❌:             eq2taut = eq2taut ? ( iiEQU | iitNONTAUT ) : 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (eq2taut)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* we may be here only in case of the second (non-taut) pass */
    // INCHI✔️❌:             /* current non-taut stereo has been found to be same as tautomeric */
    // INCHI✔️❌:             if (pINChI_Prev && pINChI_Prev->nNumberOfAtoms)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* previous component exists; output it */
    // INCHI✔️❌:                 if (bNext++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if ((nTotalCharge_Prev = pINChI_Prev->nTotalCharge)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
    // INCHI✔️❌:                     inchi_strbuf_printf( strbuf, "%+d", nTotalCharge_Prev );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (pINChI_Taut_Prev && pINChI_Taut_Prev->nNumberOfAtoms && !pINChI_Taut_Prev->bDeleted)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* previous non-taut component exists only in taut list */
    // INCHI✔️❌:                     if (bNext++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* we have found pINChI->nTotalCharge same as in pINChI_Taut */
    // INCHI✔️❌:             /* output this (current) equivalence as '*', that is, same as tautomeric */
    // INCHI✔️❌:             /* that was printed on the 1st pass. */
    // INCHI✔️❌:
    // INCHI✔️❌:             pCurrEquStr = EquString( eq2taut );
    // INCHI✔️❌:             if (multPrevEquStr && pPrevEquStr)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (pCurrEquStr && !strcmp( pCurrEquStr, pPrevEquStr ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     multPrevEquStr++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* new EqStr is different; output it */
    // INCHI✔️❌:                     if (bNext++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
    // INCHI✔️❌:                     pPrevEquStr = pCurrEquStr;
    // INCHI✔️❌:                     multPrevEquStr = 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pPrevEquStr = pCurrEquStr;
    // INCHI✔️❌:                 multPrevEquStr = 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             pINChI_Prev = NULL; /* pINChI_Prev sp2 does not exist since */
    // INCHI✔️❌:             pINChI_Taut_Prev = NULL; /* pINChI has just been printed */
    // INCHI✔️❌:             mult = 0;
    // INCHI✔️❌:             eq2tautPrev = 1;     /* pINChI_Prev sp2 does not exist */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (eq2tautPrev)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* at this point pINChI_Prev does not exist; however, pINChI */
    // INCHI✔️❌:                 /*might have been discovered and it is different from pINChI_Taut */
    // INCHI✔️❌:                 if (multPrevEquStr && pPrevEquStr)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* new EqStr is different; output it */
    // INCHI✔️❌:                     if (bNext++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
    // INCHI✔️❌:                     pPrevEquStr = NULL;
    // INCHI✔️❌:                     multPrevEquStr = 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 eq2tautPrev = 0;
    // INCHI✔️❌:                 pINChI_Prev = pINChI;
    // INCHI✔️❌:                 pINChI_Taut_Prev = pINChI_Taut;
    // INCHI✔️❌:                 mult = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* check whether pINChI and pINChI_Prev have non-zero identical stereo sp3 */
    // INCHI✔️❌:                 /*================ compare sp3 to previous =====================*/
    // INCHI✔️❌:                 eq2prev = bUseMulipliers &&
    // INCHI✔️❌:                     pINChI && pINChI_Prev &&
    // INCHI✔️❌:                     ( nTotalCharge = pINChI->nTotalCharge ) && ( nTotalCharge_Prev = pINChI_Prev->nTotalCharge ) &&
    // INCHI✔️❌:                     nTotalCharge == nTotalCharge_Prev;
    // INCHI✔️❌:                 if (eq2prev)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     mult++; /* mult = (number of non-empty equal items)-1 */
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (bNext++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (pINChI_Prev && pINChI_Prev->nNumberOfAtoms)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if ((nTotalCharge_Prev = pINChI_Prev->nTotalCharge)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* pINChI_Prev exists and has charge info */
    // INCHI✔️❌:                             MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
    // INCHI✔️❌:                             inchi_strbuf_printf( strbuf, "%+d", nTotalCharge_Prev );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         /* else charge is not present in pINChI_Prev */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* djb-rwth: removing redundant code */
    // INCHI✔️❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI✔️❌:                         if (!(bSecondNonTautPass && pINChI_Taut_Prev && pINChI_Taut_Prev->nNumberOfAtoms && !pINChI_Taut_Prev->bDeleted))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             int stop = 1;   /* <BRKPT> */
    // INCHI✔️❌:                         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 pINChI_Prev = pINChI;
    // INCHI✔️❌:                 pINChI_Taut_Prev = pINChI_Taut;
    // INCHI✔️❌:                 mult = 0; /* we do not know whether the item is empty */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( strbuf->nUsedLength - nUsedLength0 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: str_Charge2

    fn clone_inchi(
        heap: &SourceHeap,
        pointer: SourceMutPointer<crate::source_types::INChI>,
    ) -> Result<Option<crate::source_types::INChI>, SourceHeapError> {
        if pointer.is_null() {
            Ok(None)
        } else {
            Ok(Some(
                heap.slice(pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            ))
        }
    }

    fn equivalence_pointer(heap: &mut SourceHeap, value: &str) -> Result<SourceConstPointer<i8>, SourceHeapError> {
        heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain(std::iter::once(0)).collect())
            .map(SourceMutPointer::as_const)
    }

    let initial_length = string_buffer.nUsedLength;
    let component_delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let star = heap.allocate_model_storage(vec![b'*' as i8, 0])?;
    let mut previous = SourceMutPointer::null();
    let mut taut_previous = SourceMutPointer::null();
    let mut multiplier = 0_i32;
    let mut next = 0_i32;
    let mut previous_equivalence: Option<&'static str> = None;
    let mut previous_equivalence_multiplier = 0_i32;
    let mut equal_to_taut_previous = 1_i32;
    let mut index = 0_i32;

    while index <= number_of_components {
        let current = if index < number_of_components {
            selected_inchi(heap, sorted_inchi.as_const().offset(i64::from(index))?, output_type)?
        } else {
            SourceMutPointer::null()
        };
        let taut_current = if second_non_taut_pass != 0 && index < number_of_components {
            selected_inchi(heap, sorted_inchi2.as_const().offset(i64::from(index))?, OUT_T1 as i32)?
        } else {
            SourceMutPointer::null()
        };
        let current_value = clone_inchi(heap, current)?;
        let taut_current_value = clone_inchi(heap, taut_current)?;
        let equal_to_taut = second_non_taut_pass != 0
            && omit_repetitions != 0
            && current_value.as_ref().is_some_and(|value| value.nTotalCharge != 0)
            && taut_current_value
                .as_ref()
                .is_some_and(|value| value.nTotalCharge != 0 && value.bDeleted == 0)
            && current_value.as_ref().map(|value| value.nTotalCharge)
                == taut_current_value.as_ref().map(|value| value.nTotalCharge);

        if equal_to_taut {
            let previous_value = clone_inchi(heap, previous)?;
            if previous_value.as_ref().is_some_and(|value| value.nNumberOfAtoms != 0) {
                if next != 0 {
                    MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                if let Some(charge) = previous_value
                    .as_ref()
                    .map(|value| value.nTotalCharge)
                    .filter(|charge| *charge != 0)
                {
                    MakeMult(
                        heap,
                        multiplier.wrapping_add(1),
                        star.as_const(),
                        string_buffer,
                        0,
                        overflow,
                    )?;
                    append_signed_decimal(heap, string_buffer, charge)?;
                }
            } else {
                let taut_previous_value = clone_inchi(heap, taut_previous)?;
                if taut_previous_value
                    .as_ref()
                    .is_some_and(|value| value.nNumberOfAtoms != 0 && value.bDeleted == 0)
                {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                }
            }

            let current_equivalence = crate::source::base::ichiprt1::EquString(
                (crate::source_types::iiEQU | crate::source_types::iitNONTAUT) as i32,
            );
            if previous_equivalence_multiplier != 0 && previous_equivalence.is_some() {
                if previous_equivalence == Some(current_equivalence) {
                    previous_equivalence_multiplier = previous_equivalence_multiplier.wrapping_add(1);
                } else {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let pointer = equivalence_pointer(heap, previous_equivalence.unwrap())?;
                    MakeEqStr(heap, pointer, previous_equivalence_multiplier, string_buffer, overflow)?;
                    previous_equivalence = Some(current_equivalence);
                    previous_equivalence_multiplier = 1;
                }
            } else {
                previous_equivalence = Some(current_equivalence);
                previous_equivalence_multiplier = 1;
            }
            previous = SourceMutPointer::null();
            taut_previous = SourceMutPointer::null();
            multiplier = 0;
            equal_to_taut_previous = 1;
        } else if equal_to_taut_previous != 0 {
            if previous_equivalence_multiplier != 0 {
                if let Some(value) = previous_equivalence {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let pointer = equivalence_pointer(heap, value)?;
                    MakeEqStr(heap, pointer, previous_equivalence_multiplier, string_buffer, overflow)?;
                    previous_equivalence = None;
                    previous_equivalence_multiplier = 0;
                }
            }
            equal_to_taut_previous = 0;
            previous = current;
            taut_previous = taut_current;
            multiplier = 0;
        } else {
            let previous_value = clone_inchi(heap, previous)?;
            let equal_to_previous = use_multipliers != 0
                && current_value.as_ref().is_some_and(|value| value.nTotalCharge != 0)
                && previous_value.as_ref().is_some_and(|value| value.nTotalCharge != 0)
                && current_value.as_ref().map(|value| value.nTotalCharge)
                    == previous_value.as_ref().map(|value| value.nTotalCharge);
            if equal_to_previous {
                multiplier = multiplier.wrapping_add(1);
                index = index.wrapping_add(1);
                continue;
            }
            if next != 0 {
                MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
            }
            next = next.wrapping_add(1);
            if previous_value.as_ref().is_some_and(|value| value.nNumberOfAtoms != 0) {
                if let Some(charge) = previous_value
                    .as_ref()
                    .map(|value| value.nTotalCharge)
                    .filter(|charge| *charge != 0)
                {
                    MakeMult(
                        heap,
                        multiplier.wrapping_add(1),
                        star.as_const(),
                        string_buffer,
                        0,
                        overflow,
                    )?;
                    append_signed_decimal(heap, string_buffer, charge)?;
                }
            }
            previous = current;
            taut_previous = taut_current;
            multiplier = 0;
        }
        index = index.wrapping_add(1);
    }

    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_Sp2(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    sorted_inchi2: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    taut_mode: i32,
    number_of_components: i32,
    second_non_taut_pass: i32,
    omit_repetitions: i32,
    use_multipliers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:725 str_Sp2
    // INCHI✔️❌: int str_Sp2( INCHI_SORT       *pINChISort,
    // INCHI✔️❌:              INCHI_SORT       *pINChISort2,
    // INCHI✔️❌:              INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:              int              *bOverflow,
    // INCHI✔️❌:              int              bOutType,
    // INCHI✔️❌:              int              TAUT_MODE,
    // INCHI✔️❌:              int              num_components,
    // INCHI✔️❌:              int              bSecondNonTautPass,
    // INCHI✔️❌:              int              bOmitRepetitions,
    // INCHI✔️❌:              int              bUseMulipliers )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int          i, ii, ii2, nUsedLength0;
    // INCHI✔️❌:     INCHI_SORT   *is, *is2, *is0, *is20;
    // INCHI✔️❌:     INChI        *pINChI, *pINChI_Prev, *pINChI_Taut, *pINChI_Taut_Prev;
    // INCHI✔️❌:     INChI_Stereo *Stereo, *Stereo_Prev, *Stereo_Taut, *Stereo_Taut_Prev;
    // INCHI✔️❌:     int          mult, eq2prev, eq2taut, eq2tautPrev, bNext;
    // INCHI✔️❌:     const char  *pPrevEquStr, *pCurrEquStr;
    // INCHI✔️❌:     int         multPrevEquStr;
    // INCHI✔️❌:
    // INCHI✔️❌:     pINChI_Taut = NULL;
    // INCHI✔️❌:     pINChI_Prev = NULL;
    // INCHI✔️❌:     pINChI_Taut_Prev = NULL;
    // INCHI✔️❌:     mult = 0;
    // INCHI✔️❌:     bNext = 0;
    // INCHI✔️❌:     is = NULL;
    // INCHI✔️❌:     is2 = NULL;
    // INCHI✔️❌:     is0 = pINChISort;
    // INCHI✔️❌:     is20 = bSecondNonTautPass ? pINChISort2 : NULL;
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     eq2tautPrev = 1; /* pINChI_Prev (previous pINChI) does not exist */
    // INCHI✔️❌:     pPrevEquStr = NULL; /*, *pCurrEquStr;*/
    // INCHI✔️❌:     multPrevEquStr = 0;
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* For each connected component ... */
    // INCHI✔️❌:     for (i = 0; i <= num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         /* 1st (taut) pass: bOutType=OUT_TN  ; 2nd (non-taut pass) bOutType=OUT_NT */
    // INCHI✔️❌:         pINChI = ( i < num_components && ( is = is0 + i, 0 <= ( ii = GET_II( bOutType, is ) ) ) ) ? is->pINChI[ii] : NULL;
    // INCHI✔️❌:         /*================ compare sp2 to previous =====================*/
    // INCHI✔️❌:         if (bSecondNonTautPass)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* component that was output on the 1st pass */
    // INCHI✔️❌:             pINChI_Taut = ( i < num_components && ( is2 = is20 + i, 0 <= ( ii2 = GET_II( OUT_T1, is2 ) ) ) ) ? is2->pINChI[ii2] : NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /*========= if bSecondNonTautPass then compare non-iso non-taut stereo to non-iso taut ========*/
    // INCHI✔️❌:         eq2taut = 0;
    // INCHI✔️❌: #if ( FIX_EMPTY_LAYER_BUG == 1 )
    // INCHI✔️❌:         if (!eq2taut && bSecondNonTautPass && bOmitRepetitions && pINChI && pINChI_Taut)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             Stereo = pINChI->Stereo;
    // INCHI✔️❌:             Stereo_Taut = pINChI_Taut->Stereo;
    // INCHI✔️❌:             eq2taut = Stereo && Stereo_Taut &&
    // INCHI✔️❌:                 Eql_INChI_Stereo( Stereo, EQL_SP2, Stereo_Taut, EQL_SP2, 0 );
    // INCHI✔️❌:             eq2taut = eq2taut ? ( iiSTEREO | iitNONTAUT ) : 0;
    // INCHI✔️❌:
    // INCHI✔️❌:             if (!eq2taut &&
    // INCHI✔️❌:                 !Eql_INChI_Stereo( Stereo, EQL_SP2, NULL, EQL_EXISTS, 0 ) &&
    // INCHI✔️❌:                 Eql_INChI_Stereo( Stereo_Taut, EQL_SP2, NULL, EQL_EXISTS, 0 ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 eq2taut = iiEmpty; /* the current is empty while the preceding (taut) is not */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌: #else
    // INCHI✔️❌:         if (!eq2taut && bSecondNonTautPass && bOmitRepetitions)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             eq2taut = pINChI && pINChI_Taut &&
    // INCHI✔️❌:                 ( Stereo = pINChI->Stereo ) && ( Stereo_Taut = pINChI_Taut->Stereo ) &&
    // INCHI✔️❌:                 Eql_INChI_Stereo( Stereo, EQL_SP2, Stereo_Taut, EQL_SP2, 0 );
    // INCHI✔️❌:             eq2taut = eq2taut ? ( iiSTEREO | iitNONTAUT ) : 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         if (eq2taut)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* we may be here only in case of the second (non-taut) pass */
    // INCHI✔️❌:             /* current non-taut stereo has been found to be same as tautomeric */
    // INCHI✔️❌:             if (pINChI_Prev && pINChI_Prev->nNumberOfAtoms)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* previous component exists; output it */
    // INCHI✔️❌:                 if (bNext++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (( Stereo_Prev = pINChI_Prev->Stereo ) && Stereo_Prev->nNumberOfStereoBonds > 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
    // INCHI✔️❌:
    // INCHI✔️❌:                     MakeStereoString( Stereo_Prev->nBondAtom1, Stereo_Prev->nBondAtom2,
    // INCHI✔️❌:                         Stereo_Prev->b_parity,
    // INCHI✔️❌:                         0, Stereo_Prev->nNumberOfStereoBonds,
    // INCHI✔️❌:                         strbuf, TAUT_MODE, bOverflow );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (pINChI_Taut_Prev && pINChI_Taut_Prev->nNumberOfAtoms)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* previous non-taut component exists only in taut list */
    // INCHI✔️❌:                     if (bNext++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* we have found pINChI->Stereo sp2 same as in pINChI_Taut */
    // INCHI✔️❌:             /* output this (current) equivalence as '*', that is, same as tautomeric */
    // INCHI✔️❌:             /* that was printed on the 1st pass. */
    // INCHI✔️❌:             pCurrEquStr = EquString( eq2taut );
    // INCHI✔️❌:             if (multPrevEquStr && pPrevEquStr)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (pCurrEquStr && !strcmp( pCurrEquStr, pPrevEquStr ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     multPrevEquStr++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* new EqStr is different; output the previous one */
    // INCHI✔️❌:                     if (bNext++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
    // INCHI✔️❌:                     pPrevEquStr = pCurrEquStr;
    // INCHI✔️❌:                     multPrevEquStr = 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pPrevEquStr = pCurrEquStr;
    // INCHI✔️❌:                 multPrevEquStr = 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             pINChI_Prev = NULL; /* pINChI_Prev sp2 does not exist since */
    // INCHI✔️❌:             pINChI_Taut_Prev = NULL; /* pINChI has just been printed */
    // INCHI✔️❌:             mult = 0;
    // INCHI✔️❌:             eq2tautPrev = 1;     /* pINChI_Prev sp2 does not exist */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (eq2tautPrev)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* at this point pINChI_Prev does not exist; however, pINChI */
    // INCHI✔️❌:                 /*might have been discovered and it is different from pINChI_Taut */
    // INCHI✔️❌:                 if (multPrevEquStr && pPrevEquStr)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* new EqStr is different; output it */
    // INCHI✔️❌:                     if (bNext++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
    // INCHI✔️❌:                     pPrevEquStr = NULL;
    // INCHI✔️❌:                     multPrevEquStr = 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 eq2tautPrev = 0;
    // INCHI✔️❌:                 pINChI_Prev = pINChI;
    // INCHI✔️❌:                 pINChI_Taut_Prev = pINChI_Taut;
    // INCHI✔️❌:                 mult = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* check whether pINChI and pINChI_Prev have non-zero identical stereo sp2 */
    // INCHI✔️❌:                 eq2prev = bUseMulipliers &&
    // INCHI✔️❌:                     pINChI && pINChI_Prev &&
    // INCHI✔️❌:                     ( Stereo = pINChI->Stereo ) && ( Stereo_Prev = pINChI_Prev->Stereo ) &&
    // INCHI✔️❌:                     Eql_INChI_Stereo( Stereo, EQL_SP2, Stereo_Prev, EQL_SP2, 0 );
    // INCHI✔️❌:                 if (eq2prev)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     mult++; /* mult = (number of non-empty equal items)-1 */
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* pINChI sp2 info is either different or trivial. Output pINChI_Prev anyway */
    // INCHI✔️❌:                     if (bNext++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (pINChI_Prev && pINChI_Prev->nNumberOfAtoms)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (( Stereo_Prev = pINChI_Prev->Stereo ) && Stereo_Prev->nNumberOfStereoBonds > 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* pINChI_Prev exists and has sp2 info */
    // INCHI✔️❌:                             MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
    // INCHI✔️❌:
    // INCHI✔️❌:                             MakeStereoString( Stereo_Prev->nBondAtom1, Stereo_Prev->nBondAtom2,
    // INCHI✔️❌:                                 Stereo_Prev->b_parity,
    // INCHI✔️❌:                                 0, Stereo_Prev->nNumberOfStereoBonds,
    // INCHI✔️❌:                                 strbuf, TAUT_MODE, bOverflow );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         /* else sp2 info is not present in pINChI_Prev */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (bSecondNonTautPass && pINChI_Taut_Prev && pINChI_Taut_Prev->nNumberOfAtoms)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (( Stereo_Taut_Prev = pINChI_Taut_Prev->Stereo ) && Stereo_Taut_Prev->nNumberOfStereoBonds > 0)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* since pINChI_Prev does not exist, pINChI_Taut_Prev is non-tautomeric */
    // INCHI✔️❌:                                 /* and it has non-trivial sp2 info */
    // INCHI✔️❌:                                 /*
    // INCHI✔️❌:                                 tot_len += MakeDelim( sIdenticalValues, strbuf, bOverflow);
    // INCHI✔️❌:                                 */
    // INCHI✔️❌:                                 ;/* pINChI_Taut_Prev sp2 info was output in the main stereo section */
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 ; /* pINChI_Taut_Prev exists and has not sp2 info */
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             int stop = 1;   /* <BRKPT> */
    // INCHI✔️❌:                         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 pINChI_Prev = pINChI;
    // INCHI✔️❌:                 pINChI_Taut_Prev = pINChI_Taut;
    // INCHI✔️❌:                 mult = 0; /* we do not know whether the item is empty */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }    /* end of for each connected component ... */
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( strbuf->nUsedLength - nUsedLength0 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: str_Sp2

    fn clone_inchi(
        heap: &SourceHeap,
        pointer: SourceMutPointer<crate::source_types::INChI>,
    ) -> Result<Option<crate::source_types::INChI>, SourceHeapError> {
        if pointer.is_null() {
            Ok(None)
        } else {
            Ok(Some(
                heap.slice(pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            ))
        }
    }

    fn clone_stereo(
        heap: &SourceHeap,
        inchi: Option<&crate::source_types::INChI>,
    ) -> Result<Option<crate::source_types::INChI_Stereo>, SourceHeapError> {
        let Some(inchi) = inchi else {
            return Ok(None);
        };
        if inchi.Stereo.is_null() {
            Ok(None)
        } else {
            Ok(Some(
                heap.slice(inchi.Stereo.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            ))
        }
    }

    fn equivalence_pointer(heap: &mut SourceHeap, value: &str) -> Result<SourceConstPointer<i8>, SourceHeapError> {
        heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain(std::iter::once(0)).collect())
            .map(SourceMutPointer::as_const)
    }

    fn output_stereo(
        heap: &mut SourceHeap,
        inchi: Option<&crate::source_types::INChI>,
        multiplier: i32,
        star: SourceConstPointer<i8>,
        string_buffer: &mut INCHI_IOS_STRING,
        taut_mode: i32,
        overflow: &mut i32,
    ) -> Result<(), SourceHeapError> {
        let Some(inchi) = inchi.filter(|inchi| inchi.nNumberOfAtoms != 0) else {
            return Ok(());
        };
        let Some(stereo) = clone_stereo(heap, Some(inchi))? else {
            return Ok(());
        };
        if stereo.nNumberOfStereoBonds > 0 {
            MakeMult(heap, multiplier, star, string_buffer, 0, overflow)?;
            MakeStereoString(
                heap,
                stereo.nBondAtom1.as_const(),
                stereo.nBondAtom2.as_const(),
                stereo.b_parity.as_const(),
                0,
                stereo.nNumberOfStereoBonds,
                string_buffer,
                taut_mode,
                overflow,
            )?;
        }
        Ok(())
    }

    let initial_length = string_buffer.nUsedLength;
    let component_delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let star = heap.allocate_model_storage(vec![b'*' as i8, 0])?;
    let mut previous = SourceMutPointer::null();
    let mut taut_previous = SourceMutPointer::null();
    let mut multiplier = 0_i32;
    let mut next = 0_i32;
    let mut previous_equivalence: Option<&'static str> = None;
    let mut previous_equivalence_multiplier = 0_i32;
    let mut equal_to_taut_previous = 1_i32;
    let mut index = 0_i32;

    while index <= number_of_components {
        let current = if index < number_of_components {
            selected_inchi(heap, sorted_inchi.as_const().offset(i64::from(index))?, output_type)?
        } else {
            SourceMutPointer::null()
        };
        let taut_current = if second_non_taut_pass != 0 && index < number_of_components {
            selected_inchi(heap, sorted_inchi2.as_const().offset(i64::from(index))?, OUT_T1 as i32)?
        } else {
            SourceMutPointer::null()
        };
        let current_value = clone_inchi(heap, current)?;
        let taut_current_value = clone_inchi(heap, taut_current)?;
        let current_stereo = clone_stereo(heap, current_value.as_ref())?;
        let taut_current_stereo = clone_stereo(heap, taut_current_value.as_ref())?;
        let equal_to_taut = second_non_taut_pass != 0
            && omit_repetitions != 0
            && current_value.is_some()
            && taut_current_value.is_some()
            && current_stereo.is_some()
            && taut_current_stereo.is_some()
            && Eql_INChI_Stereo(
                heap,
                current_stereo.as_ref(),
                crate::source_types::EQL_SP2 as i32,
                taut_current_stereo.as_ref(),
                crate::source_types::EQL_SP2 as i32,
                0,
            )? != 0;

        if equal_to_taut {
            let previous_value = clone_inchi(heap, previous)?;
            if previous_value.as_ref().is_some_and(|inchi| inchi.nNumberOfAtoms != 0) {
                if next != 0 {
                    MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                output_stereo(
                    heap,
                    previous_value.as_ref(),
                    multiplier.wrapping_add(1),
                    star.as_const(),
                    string_buffer,
                    taut_mode,
                    overflow,
                )?;
            } else {
                let taut_previous_value = clone_inchi(heap, taut_previous)?;
                if taut_previous_value
                    .as_ref()
                    .is_some_and(|inchi| inchi.nNumberOfAtoms != 0)
                {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                }
            }
            let current_equivalence = crate::source::base::ichiprt1::EquString(
                (crate::source_types::iiSTEREO | crate::source_types::iitNONTAUT) as i32,
            );
            if previous_equivalence_multiplier != 0 && previous_equivalence.is_some() {
                if previous_equivalence == Some(current_equivalence) {
                    previous_equivalence_multiplier = previous_equivalence_multiplier.wrapping_add(1);
                } else {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let pointer = equivalence_pointer(heap, previous_equivalence.unwrap())?;
                    MakeEqStr(heap, pointer, previous_equivalence_multiplier, string_buffer, overflow)?;
                    previous_equivalence = Some(current_equivalence);
                    previous_equivalence_multiplier = 1;
                }
            } else {
                previous_equivalence = Some(current_equivalence);
                previous_equivalence_multiplier = 1;
            }
            previous = SourceMutPointer::null();
            taut_previous = SourceMutPointer::null();
            multiplier = 0;
            equal_to_taut_previous = 1;
        } else if equal_to_taut_previous != 0 {
            if previous_equivalence_multiplier != 0 {
                if let Some(value) = previous_equivalence {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let pointer = equivalence_pointer(heap, value)?;
                    MakeEqStr(heap, pointer, previous_equivalence_multiplier, string_buffer, overflow)?;
                    previous_equivalence = None;
                    previous_equivalence_multiplier = 0;
                }
            }
            equal_to_taut_previous = 0;
            previous = current;
            taut_previous = taut_current;
            multiplier = 0;
        } else {
            let previous_value = clone_inchi(heap, previous)?;
            let previous_stereo = clone_stereo(heap, previous_value.as_ref())?;
            let equal_to_previous = use_multipliers != 0
                && current_value.is_some()
                && previous_value.is_some()
                && current_stereo.is_some()
                && previous_stereo.is_some()
                && Eql_INChI_Stereo(
                    heap,
                    current_stereo.as_ref(),
                    crate::source_types::EQL_SP2 as i32,
                    previous_stereo.as_ref(),
                    crate::source_types::EQL_SP2 as i32,
                    0,
                )? != 0;
            if equal_to_previous {
                multiplier = multiplier.wrapping_add(1);
                index = index.wrapping_add(1);
                continue;
            }
            if next != 0 {
                MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
            }
            next = next.wrapping_add(1);
            output_stereo(
                heap,
                previous_value.as_ref(),
                multiplier.wrapping_add(1),
                star.as_const(),
                string_buffer,
                taut_mode,
                overflow,
            )?;
            previous = current;
            taut_previous = taut_current;
            multiplier = 0;
        }
        index = index.wrapping_add(1);
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_Sp3(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    sorted_inchi2: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    taut_mode: i32,
    number_of_components: i32,
    mut relative_or_racemic: i32,
    second_non_taut_pass: i32,
    omit_repetitions: i32,
    use_multipliers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:957 str_Sp3
    // INCHI✔️❌: int str_Sp3( INCHI_SORT       *pINChISort,
    // INCHI✔️❌:              INCHI_SORT       *pINChISort2,
    // INCHI✔️❌:              INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:              int              *bOverflow,
    // INCHI✔️❌:              int              bOutType,
    // INCHI✔️❌:              int              TAUT_MODE,
    // INCHI✔️❌:              int              num_components,
    // INCHI✔️❌:              int              bRelRac,
    // INCHI✔️❌:              int              bSecondNonTautPass,
    // INCHI✔️❌:              int              bOmitRepetitions,
    // INCHI✔️❌:              int              bUseMulipliers )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int          i, ii, ii2, nUsedLength0;
    // INCHI✔️❌:     INCHI_SORT   *is, *is2, *is0, *is20;
    // INCHI✔️❌:     INChI        *pINChI, *pINChI_Prev, *pINChI_Taut, *pINChI_Taut_Prev;
    // INCHI✔️❌:     INChI_Stereo *Stereo, *Stereo_Prev, *Stereo_Taut, *Stereo_Taut_Prev;
    // INCHI✔️❌:     int          mult, eq2prev, eq2taut, eq2tautPrev, bNext;
    // INCHI✔️❌:     const char  *pPrevEquStr, *pCurrEquStr;
    // INCHI✔️❌:     int         multPrevEquStr;
    // INCHI✔️❌:     pINChI_Taut = NULL;
    // INCHI✔️❌:     pINChI_Prev = NULL;
    // INCHI✔️❌:     pINChI_Taut_Prev = NULL;
    // INCHI✔️❌:     mult = 0;
    // INCHI✔️❌:     bNext = 0;
    // INCHI✔️❌:     is = NULL;
    // INCHI✔️❌:     is2 = NULL;
    // INCHI✔️❌:     is0 = pINChISort;
    // INCHI✔️❌:     is20 = bSecondNonTautPass ? pINChISort2 : NULL;
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     eq2tautPrev = 1; /* pINChI_Prev (previous pINChI) does not exist */
    // INCHI✔️❌:     pPrevEquStr = NULL; /*, *pCurrEquStr;*/
    // INCHI✔️❌:     multPrevEquStr = 0;
    // INCHI✔️❌: #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    // INCHI✔️❌: #else
    // INCHI✔️❌:     bRelRac = 0;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* For each connected component...    */
    // INCHI✔️❌:     for (i = 0; i <= num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         /* 1st (taut) pass: bOutType=OUT_TN  ; 2nd (non-taut pass) bOutType=OUT_NT */
    // INCHI✔️❌:         pINChI = ( i < num_components && ( is = is0 + i, 0 <= ( ii = GET_II( bOutType, is ) ) ) ) ? is->pINChI[ii] : NULL;
    // INCHI✔️❌:         /*================ compare sp3 to previous =====================*/
    // INCHI✔️❌:         if (bSecondNonTautPass)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* component that was output on the 1st pass */
    // INCHI✔️❌:             pINChI_Taut = ( i < num_components && ( is2 = is20 + i, 0 <= ( ii2 = GET_II( OUT_T1, is2 ) ) ) ) ? is2->pINChI[ii2] : NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /*========= if bSecondNonTautPass then compare non-iso non-taut stereo to non-iso taut ========*/
    // INCHI✔️❌:         eq2taut = 0;
    // INCHI✔️❌: #if ( FIX_EMPTY_LAYER_BUG == 1 )
    // INCHI✔️❌:         if (!eq2taut && bSecondNonTautPass && bOmitRepetitions && pINChI && pINChI_Taut)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             Stereo = pINChI->Stereo;
    // INCHI✔️❌:             Stereo_Taut = pINChI_Taut->Stereo;
    // INCHI✔️❌:             eq2taut = Stereo && Stereo_Taut &&
    // INCHI✔️❌:                 Eql_INChI_Stereo( Stereo, EQL_SP3, Stereo_Taut, EQL_SP3, bRelRac );
    // INCHI✔️❌:             eq2taut = eq2taut ? ( iiSTEREO | iitNONTAUT ) : 0;
    // INCHI✔️❌:             if (!eq2taut &&
    // INCHI✔️❌:                 !Eql_INChI_Stereo( Stereo, EQL_SP3, NULL, EQL_EXISTS, 0 ) &&
    // INCHI✔️❌:                 Eql_INChI_Stereo( Stereo_Taut, EQL_SP3, NULL, EQL_EXISTS, 0 ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 eq2taut = iiEmpty; /* the current is empty while the preceding (taut) is not */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌: #else
    // INCHI✔️❌:         if (!eq2taut && bSecondNonTautPass && bOmitRepetitions)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             eq2taut = pINChI && pINChI_Taut &&
    // INCHI✔️❌:                 ( Stereo = pINChI->Stereo ) && ( Stereo_Taut = pINChI_Taut->Stereo ) &&
    // INCHI✔️❌:                 Eql_INChI_Stereo( Stereo, EQL_SP3, Stereo_Taut, EQL_SP3, bRelRac );
    // INCHI✔️❌:             eq2taut = eq2taut ? ( iiSTEREO | iitNONTAUT ) : 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         if (eq2taut)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* we may be here only in case of the second (non-taut) pass */
    // INCHI✔️❌:             /* current non-taut stereo has been found to be same as tautomeric */
    // INCHI✔️❌:             if (pINChI_Prev && pINChI_Prev->nNumberOfAtoms)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* previous component exists; output it */
    // INCHI✔️❌:                 if (bNext++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (( Stereo_Prev = pINChI_Prev->Stereo ) && Stereo_Prev->nNumberOfStereoCenters > 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* non-empty item */
    // INCHI✔️❌:                     MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
    // INCHI✔️❌:                     MakeStereoString( Stereo_Prev->nNumber, NULL,
    // INCHI✔️❌:                                       Stereo_Prev->t_parity,0,
    // INCHI✔️❌:                                       Stereo_Prev->nNumberOfStereoCenters,
    // INCHI✔️❌:                                       strbuf, TAUT_MODE, bOverflow );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:                 if (pINChI_Taut_Prev && pINChI_Taut_Prev->nNumberOfAtoms)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* previous non-taut component exists only in taut list */
    // INCHI✔️❌:                     if (bNext++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             /* we have found pINChI->Stereo sp3 same as in pINChI_Taut */
    // INCHI✔️❌:             /* output this (current) equivalence as '*', that is, same as tautomeric */
    // INCHI✔️❌:             /* that was printed on the 1st pass. */
    // INCHI✔️❌:
    // INCHI✔️❌:             pCurrEquStr = EquString( eq2taut );
    // INCHI✔️❌:             if (multPrevEquStr && pPrevEquStr)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (pCurrEquStr && !strcmp( pCurrEquStr, pPrevEquStr ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     multPrevEquStr++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* new EqStr is different; output it */
    // INCHI✔️❌:                     if (bNext++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
    // INCHI✔️❌:                     pPrevEquStr = pCurrEquStr;
    // INCHI✔️❌:                     multPrevEquStr = 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pPrevEquStr = pCurrEquStr;
    // INCHI✔️❌:                 multPrevEquStr = 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             pINChI_Prev = NULL; /* pINChI_Prev sp2 does not exist since */
    // INCHI✔️❌:             pINChI_Taut_Prev = NULL; /* pINChI has just been printed */
    // INCHI✔️❌:             mult = 0;
    // INCHI✔️❌:             eq2tautPrev = 1;     /* pINChI_Prev sp2 does not exist */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (eq2tautPrev)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* at this point pINChI_Prev does not exist; however, pINChI */
    // INCHI✔️❌:                 /*might have been discovered and it is different from pINChI_Taut */
    // INCHI✔️❌:                 if (multPrevEquStr && pPrevEquStr)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* new EqStr is different; output it */
    // INCHI✔️❌:                     if (bNext++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
    // INCHI✔️❌:                     pPrevEquStr = NULL;
    // INCHI✔️❌:                     multPrevEquStr = 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 eq2tautPrev = 0;
    // INCHI✔️❌:                 pINChI_Prev = pINChI;
    // INCHI✔️❌:                 pINChI_Taut_Prev = pINChI_Taut;
    // INCHI✔️❌:                 mult = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* check whether pINChI and pINChI_Prev have non-zero identical stereo sp3 */
    // INCHI✔️❌:                 /*================ compare sp3 to previous =====================*/
    // INCHI✔️❌:                 eq2prev = bUseMulipliers &&
    // INCHI✔️❌:                     pINChI && pINChI_Prev &&
    // INCHI✔️❌:                     ( Stereo = pINChI->Stereo ) && ( Stereo_Prev = pINChI_Prev->Stereo ) &&
    // INCHI✔️❌:                     Eql_INChI_Stereo( Stereo, EQL_SP3, Stereo_Prev, EQL_SP3, bRelRac );
    // INCHI✔️❌:                 if (eq2prev)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     mult++; /* mult = (number of non-empty equal items)-1 */
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (bNext++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (pINChI_Prev && pINChI_Prev->nNumberOfAtoms)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (( Stereo_Prev = pINChI_Prev->Stereo ) && Stereo_Prev->nNumberOfStereoCenters > bRelRac)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* pINChI_Prev exists and has sp3 info */
    // INCHI✔️❌:                             MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
    // INCHI✔️❌:
    // INCHI✔️❌:                             MakeStereoString( Stereo_Prev->nNumber, NULL, Stereo_Prev->t_parity,
    // INCHI✔️❌:                                 0, Stereo_Prev->nNumberOfStereoCenters,
    // INCHI✔️❌:                                 strbuf, TAUT_MODE, bOverflow );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         /* else sp3 info is not present in pINChI_Prev */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (bSecondNonTautPass && pINChI_Taut_Prev && pINChI_Taut_Prev->nNumberOfAtoms)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (( Stereo_Taut_Prev = pINChI_Taut_Prev->Stereo ) && Stereo_Taut_Prev->nNumberOfStereoCenters > bRelRac)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* since pINChI_Prev does not exist, pINChI_Taut_Prev is non-tautomeric */
    // INCHI✔️❌:                                 /* and it has non-trivial sp3 info. This info has already been printed in the main section */
    // INCHI✔️❌:                                 /*
    // INCHI✔️❌:                                 tot_len += MakeDelim( sIdenticalValues, strbuf, bOverflow);
    // INCHI✔️❌:                                 */
    // INCHI✔️❌:                                 ; /* pINChI_Taut_Prev sp3 info was output in the main stereo section */
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 ; /* pINChI_Taut_Prev exists and has not sp3 info */
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             int stop = 1;   /* <BRKPT> */
    // INCHI✔️❌:                         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 pINChI_Prev = pINChI;
    // INCHI✔️❌:                 pINChI_Taut_Prev = pINChI_Taut;
    // INCHI✔️❌:                 mult = 0; /* we do not know whether the item is empty */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( strbuf->nUsedLength - nUsedLength0 );
    // INCHI✔️❌: }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: /****************************************************************************
    // INCHI✔️❌:   Output abs stero inversion substring of the whole structure InChI string
    // END INCHI C FUNCTION: str_Sp3

    fn clone_inchi(
        heap: &SourceHeap,
        pointer: SourceMutPointer<crate::source_types::INChI>,
    ) -> Result<Option<crate::source_types::INChI>, SourceHeapError> {
        if pointer.is_null() {
            Ok(None)
        } else {
            Ok(Some(
                heap.slice(pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            ))
        }
    }

    fn clone_stereo(
        heap: &SourceHeap,
        inchi: Option<&crate::source_types::INChI>,
    ) -> Result<Option<crate::source_types::INChI_Stereo>, SourceHeapError> {
        let Some(inchi) = inchi else {
            return Ok(None);
        };
        if inchi.Stereo.is_null() {
            Ok(None)
        } else {
            Ok(Some(
                heap.slice(inchi.Stereo.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            ))
        }
    }

    fn equivalence_pointer(heap: &mut SourceHeap, value: &str) -> Result<SourceConstPointer<i8>, SourceHeapError> {
        heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain(std::iter::once(0)).collect())
            .map(SourceMutPointer::as_const)
    }

    fn output_stereo(
        heap: &mut SourceHeap,
        inchi: Option<&crate::source_types::INChI>,
        multiplier: i32,
        star: SourceConstPointer<i8>,
        string_buffer: &mut INCHI_IOS_STRING,
        taut_mode: i32,
        overflow: &mut i32,
        relative_or_racemic: i32,
    ) -> Result<(), SourceHeapError> {
        let Some(inchi) = inchi.filter(|inchi| inchi.nNumberOfAtoms != 0) else {
            return Ok(());
        };
        let Some(stereo) = clone_stereo(heap, Some(inchi))? else {
            return Ok(());
        };
        if stereo.nNumberOfStereoCenters > relative_or_racemic {
            MakeMult(heap, multiplier, star, string_buffer, 0, overflow)?;
            MakeStereoString(
                heap,
                stereo.nNumber.as_const(),
                SourceConstPointer::null(),
                stereo.t_parity.as_const(),
                0,
                stereo.nNumberOfStereoCenters,
                string_buffer,
                taut_mode,
                overflow,
            )?;
        }
        Ok(())
    }

    relative_or_racemic = 0;
    let initial_length = string_buffer.nUsedLength;
    let component_delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let star = heap.allocate_model_storage(vec![b'*' as i8, 0])?;
    let mut previous = SourceMutPointer::null();
    let mut taut_previous = SourceMutPointer::null();
    let mut multiplier = 0_i32;
    let mut next = 0_i32;
    let mut previous_equivalence: Option<&'static str> = None;
    let mut previous_equivalence_multiplier = 0_i32;
    let mut equal_to_taut_previous = 1_i32;
    let mut index = 0_i32;

    while index <= number_of_components {
        let current = if index < number_of_components {
            selected_inchi(heap, sorted_inchi.as_const().offset(i64::from(index))?, output_type)?
        } else {
            SourceMutPointer::null()
        };
        let taut_current = if second_non_taut_pass != 0 && index < number_of_components {
            selected_inchi(heap, sorted_inchi2.as_const().offset(i64::from(index))?, OUT_T1 as i32)?
        } else {
            SourceMutPointer::null()
        };
        let current_value = clone_inchi(heap, current)?;
        let taut_current_value = clone_inchi(heap, taut_current)?;
        let current_stereo = clone_stereo(heap, current_value.as_ref())?;
        let taut_current_stereo = clone_stereo(heap, taut_current_value.as_ref())?;
        let equal_to_taut = second_non_taut_pass != 0
            && omit_repetitions != 0
            && current_value.is_some()
            && taut_current_value.is_some()
            && current_stereo.is_some()
            && taut_current_stereo.is_some()
            && Eql_INChI_Stereo(
                heap,
                current_stereo.as_ref(),
                crate::source_types::EQL_SP3 as i32,
                taut_current_stereo.as_ref(),
                crate::source_types::EQL_SP3 as i32,
                relative_or_racemic,
            )? != 0;

        if equal_to_taut {
            let previous_value = clone_inchi(heap, previous)?;
            if previous_value.as_ref().is_some_and(|inchi| inchi.nNumberOfAtoms != 0) {
                if next != 0 {
                    MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                output_stereo(
                    heap,
                    previous_value.as_ref(),
                    multiplier.wrapping_add(1),
                    star.as_const(),
                    string_buffer,
                    taut_mode,
                    overflow,
                    relative_or_racemic,
                )?;
            } else {
                let taut_previous_value = clone_inchi(heap, taut_previous)?;
                if taut_previous_value
                    .as_ref()
                    .is_some_and(|inchi| inchi.nNumberOfAtoms != 0)
                {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                }
            }
            let current_equivalence = crate::source::base::ichiprt1::EquString(
                (crate::source_types::iiSTEREO | crate::source_types::iitNONTAUT) as i32,
            );
            if previous_equivalence_multiplier != 0 && previous_equivalence.is_some() {
                if previous_equivalence == Some(current_equivalence) {
                    previous_equivalence_multiplier = previous_equivalence_multiplier.wrapping_add(1);
                } else {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let pointer = equivalence_pointer(heap, previous_equivalence.unwrap())?;
                    MakeEqStr(heap, pointer, previous_equivalence_multiplier, string_buffer, overflow)?;
                    previous_equivalence = Some(current_equivalence);
                    previous_equivalence_multiplier = 1;
                }
            } else {
                previous_equivalence = Some(current_equivalence);
                previous_equivalence_multiplier = 1;
            }
            previous = SourceMutPointer::null();
            taut_previous = SourceMutPointer::null();
            multiplier = 0;
            equal_to_taut_previous = 1;
        } else if equal_to_taut_previous != 0 {
            if previous_equivalence_multiplier != 0 {
                if let Some(value) = previous_equivalence {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let pointer = equivalence_pointer(heap, value)?;
                    MakeEqStr(heap, pointer, previous_equivalence_multiplier, string_buffer, overflow)?;
                    previous_equivalence = None;
                    previous_equivalence_multiplier = 0;
                }
            }
            equal_to_taut_previous = 0;
            previous = current;
            taut_previous = taut_current;
            multiplier = 0;
        } else {
            let previous_value = clone_inchi(heap, previous)?;
            let previous_stereo = clone_stereo(heap, previous_value.as_ref())?;
            let equal_to_previous = use_multipliers != 0
                && current_value.is_some()
                && previous_value.is_some()
                && current_stereo.is_some()
                && previous_stereo.is_some()
                && Eql_INChI_Stereo(
                    heap,
                    current_stereo.as_ref(),
                    crate::source_types::EQL_SP3 as i32,
                    previous_stereo.as_ref(),
                    crate::source_types::EQL_SP3 as i32,
                    relative_or_racemic,
                )? != 0;
            if equal_to_previous {
                multiplier = multiplier.wrapping_add(1);
                index = index.wrapping_add(1);
                continue;
            }
            if next != 0 {
                MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
            }
            next = next.wrapping_add(1);
            output_stereo(
                heap,
                previous_value.as_ref(),
                multiplier.wrapping_add(1),
                star.as_const(),
                string_buffer,
                taut_mode,
                overflow,
                relative_or_racemic,
            )?;
            previous = current;
            taut_previous = taut_current;
            multiplier = 0;
        }
        index = index.wrapping_add(1);
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case)]
pub(crate) fn str_StereoAbsInv(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    number_of_components: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:1191 str_StereoAbsInv
    // INCHI✔️❌: int str_StereoAbsInv( INCHI_SORT       *pINChISort,
    // INCHI✔️❌:                       INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                       int              *bOverflow,
    // INCHI✔️❌:                       int              bOutType,
    // INCHI✔️❌:                       int              num_components )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int          i, j, ii, nUsedLength0;
    // INCHI✔️❌:     INCHI_SORT   *is, *is0;
    // INCHI✔️❌:     INChI_Stereo *Stereo;
    // INCHI✔️❌:     INChI        *pINChI;
    // INCHI✔️❌:
    // INCHI✔️❌:     is = NULL;
    // INCHI✔️❌:     is0 = pINChISort;
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* For each connected component...    */
    // INCHI✔️❌:     for (i = 0; !*bOverflow && i < num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         is = is0 + i;
    // INCHI✔️❌:         pINChI = ( 0 <= ( ii = GET_II( bOutType, is ) ) ) ? is->pINChI[ii] : NULL;
    // INCHI✔️❌:         if (pINChI && ( Stereo = pINChI->Stereo ) && ( j = Stereo->nCompInv2Abs ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             MakeDelim( j < 0 ? "1" : "0", strbuf, bOverflow );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             MakeDelim( ".", strbuf, bOverflow );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( strbuf->nUsedLength - nUsedLength0 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: str_StereoAbsInv

    let initial_length = string_buffer.nUsedLength;
    let inverted = heap.allocate_model_storage(vec![b'1' as i8, 0])?;
    let absolute = heap.allocate_model_storage(vec![b'0' as i8, 0])?;
    let absent = heap.allocate_model_storage(vec![b'.' as i8, 0])?;
    let mut index = 0_i32;
    while *overflow == 0 && index < number_of_components {
        let inchi = selected_inchi(heap, sorted_inchi.as_const().offset(i64::from(index))?, output_type)?;
        let inversion = if inchi.is_null() {
            None
        } else {
            let inchi = heap
                .slice(inchi.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if inchi.Stereo.is_null() {
                None
            } else {
                Some(
                    heap.slice(inchi.Stereo.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCompInv2Abs,
                )
            }
        };
        let delimiter = match inversion {
            Some(value) if value < 0 => inverted.as_const(),
            Some(value) if value > 0 => absolute.as_const(),
            _ => absent.as_const(),
        };
        MakeDelim(heap, delimiter, string_buffer, overflow)?;
        index = index.wrapping_add(1);
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case)]
pub(crate) fn str_IsoStereoAbsInv(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    number_of_components: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:2058 str_IsoStereoAbsInv
    // INCHI✔️❌: int str_IsoStereoAbsInv( INCHI_SORT       *pINChISort,
    // INCHI✔️❌:                          INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                          int              *bOverflow,
    // INCHI✔️❌:                          int              bOutType,
    // INCHI✔️❌:                          int              num_components )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int          i, j, ii, nUsedLength0;
    // INCHI✔️❌:     INCHI_SORT   *is, *is0;
    // INCHI✔️❌:     INChI_Stereo *Stereo;
    // INCHI✔️❌:     INChI        *pINChI;
    // INCHI✔️❌:
    // INCHI✔️❌:     is = NULL;
    // INCHI✔️❌:     is0 = pINChISort;
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* For each connected component...    */
    // INCHI✔️❌:     for (i = 0; !*bOverflow && i < num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         is = is0 + i;
    // INCHI✔️❌:         pINChI = ( 0 <= ( ii = GET_II( bOutType, is ) ) )
    // INCHI✔️❌:             ? is->pINChI[ii]
    // INCHI✔️❌:             : NULL;
    // INCHI✔️❌:         if (pINChI &&
    // INCHI✔️❌:             ( Stereo = pINChI->StereoIsotopic ) &&
    // INCHI✔️❌:             ( j = Stereo->nCompInv2Abs ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             MakeDelim( j < 0 ? "1" : "0", strbuf, bOverflow );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             MakeDelim( ".", strbuf, bOverflow );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( strbuf->nUsedLength - nUsedLength0 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: str_IsoStereoAbsInv

    let initial_length = string_buffer.nUsedLength;
    let inverted = heap.allocate_model_storage(vec![b'1' as i8, 0])?;
    let absolute = heap.allocate_model_storage(vec![b'0' as i8, 0])?;
    let absent = heap.allocate_model_storage(vec![b'.' as i8, 0])?;
    let mut index = 0_i32;
    while *overflow == 0 && index < number_of_components {
        let inchi = selected_inchi(heap, sorted_inchi.as_const().offset(i64::from(index))?, output_type)?;
        let inversion = if inchi.is_null() {
            None
        } else {
            let inchi = heap
                .slice(inchi.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if inchi.StereoIsotopic.is_null() {
                None
            } else {
                Some(
                    heap.slice(inchi.StereoIsotopic.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCompInv2Abs,
                )
            }
        };
        let delimiter = match inversion {
            Some(value) if value < 0 => inverted.as_const(),
            Some(value) if value > 0 => absolute.as_const(),
            _ => absent.as_const(),
        };
        MakeDelim(heap, delimiter, string_buffer, overflow)?;
        index = index.wrapping_add(1);
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_IsoAtoms(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    sorted_inchi2: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    taut_mode: i32,
    number_of_components: i32,
    abc_numbers: i32,
    second_non_taut_pass: i32,
    omit_repetitions: i32,
    use_multipliers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:1229 str_IsoAtoms
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int str_IsoAtoms( INCHI_SORT *pINChISort,
                      INCHI_SORT *pINChISort2,
                      INCHI_IOS_STRING *strbuf,
                      int *bOverflow,
                      int bOutType,
                      int TAUT_MODE,
                      int num_components,
                      int bAbcNumbers,
                      int bSecondNonTautPass,
                      int bOmitRepetitions,
                      int bUseMulipliers )
    {
        int          i, ii, ii2, nUsedLength0;
        INCHI_SORT   *is, *is2, *is0, *is20;
        INChI        *pINChI, *pINChI_Prev, *pINChI_Taut, *pINChI_Taut_Prev;
        int          mult, eq2prev, eq2taut, eq2tautPrev, bNext;
        const char  *pPrevEquStr, *pCurrEquStr;
        int         multPrevEquStr;
        pINChI_Taut = NULL;
        pINChI_Prev = NULL;
        pINChI_Taut_Prev = NULL;
        mult = 0;
        bNext = 0;
        is = NULL;
        is2 = NULL;
        is0 = pINChISort;
        is20 = bSecondNonTautPass ? pINChISort2 : NULL;
        /* djb-rwth: removing redundant code */
        eq2tautPrev = 1; /* pINChI_Prev (previous pINChI) does not exist */
        pPrevEquStr = NULL; /*, *pCurrEquStr;*/
        multPrevEquStr = 0;
        nUsedLength0 = strbuf->nUsedLength;

        /* For each connected component...    */
        for (i = 0; i <= num_components; i++)
        {

            /* 1st (taut) pass: bOutType=OUT_TN  ; 2nd (non-taut pass) bOutType=OUT_NT */
            pINChI = ( i < num_components && ( is = is0 + i, 0 <= ( ii = GET_II( bOutType, is ) ) ) ) ? is->pINChI[ii] : NULL;
            /*================ compare isotopic info to previous component =====================*/
            if (bSecondNonTautPass)
            {
                /* component that was output on the 1st pass */
                pINChI_Taut = ( i < num_components && ( is2 = is20 + i, 0 <= ( ii2 = GET_II( OUT_T1, is2 ) ) ) ) ? is2->pINChI[ii2] : NULL;
            }
            /*========= if bSecondNonTautPass then compare iso non-taut to taut non-iso ========*/
            eq2taut = 0;
            if (!eq2taut && bSecondNonTautPass && bOmitRepetitions)
            {
                eq2taut = Eql_INChI_Isotopic( pINChI, pINChI_Taut );
                eq2taut = eq2taut ? ( iiNUMB | iitNONTAUT ) : 0;
            }
            if (eq2taut)
            {
                /* we may be here only in case of the second (non-taut) pass */
                /* current non-taut isotopic info has been found to be same as current tautomeric */
                if (pINChI_Prev && pINChI_Prev->nNumberOfAtoms)
                {
                    /* previous component exists; output it */
                    if (bNext++)
                    {
                        MakeDelim( sCompDelim, strbuf, bOverflow );
                    }
                    if (pINChI_Prev && ( pINChI_Prev->nNumberOfIsotopicAtoms > 0 ||
                        pINChI_Prev->nNumberOfIsotopicTGroups > 0 ))
                    {
                        /* non-empty item */
                        MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
                        /*  Isotopic atoms */
                        if (pINChI_Prev->nNumberOfIsotopicAtoms > 0/* && nStrLen-tot_len > 2*/ && !*bOverflow)
                        { /* dereferenced bOverflow 2004-06-07 */
                            MakeIsoAtomString( pINChI_Prev->IsotopicAtom,
                                pINChI_Prev->nNumberOfIsotopicAtoms,
                                strbuf,
                                TAUT_MODE, bOverflow );
                        }
                        /*  Isotopic tautomeric groups */
                        if (pINChI_Prev->nNumberOfIsotopicTGroups > 0 &&
                            /*nStrLen-tot_len > 3 && */
                            !*bOverflow)
                        {
                            MakeDelim( bAbcNumbers ? ITEM_DELIMETER : "(", strbuf, bOverflow );
                            MakeIsoTautString( pINChI_Prev->IsotopicTGroup, pINChI_Prev->nNumberOfIsotopicTGroups,
                                strbuf, TAUT_MODE, bOverflow );
                            if (!bAbcNumbers)
                            {
                                MakeDelim( ")", strbuf, bOverflow );
                            }
                        }
                    }
                }
                else
                {
                    if (pINChI_Taut_Prev && pINChI_Taut_Prev->nNumberOfAtoms)
                    {
                        /* previous non-taut component exists only in taut list */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                    }
                }
                /* we have found pINChI isotopic info to be same as in pINChI_Taut */
                /* output this (current) equivalence as '*', that is, same as tautomeric */
                /* that was printed on the 1st pass. */
                pCurrEquStr = EquString( eq2taut );
                if (multPrevEquStr && pPrevEquStr)
                {
                    if (pCurrEquStr && !strcmp( pCurrEquStr, pPrevEquStr ))
                    {
                        multPrevEquStr++;
                    }
                    else
                    {
                        /* new EqStr is different; output it */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                        pPrevEquStr = pCurrEquStr;
                        multPrevEquStr = 1;
                    }
                }
                else
                {
                    pPrevEquStr = pCurrEquStr;
                    multPrevEquStr = 1;
                }
                pINChI_Prev = NULL; /* pINChI_Prev isotopic info does not exist since */
                pINChI_Taut_Prev = NULL; /* pINChI has just been printed */
                mult = 0;
                eq2tautPrev = 1;     /* pINChI_Prev isotopic info does not exist */
            }
            else
                if (eq2tautPrev)
                {
                    /* at this point pINChI_Prev does not exist; however, pINChI */
                    /* might have been discovered and it is different from pINChI_Taut */
                    if (multPrevEquStr && pPrevEquStr)
                    {
                        /* new EqStr is different; output it */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                        pPrevEquStr = NULL;
                        multPrevEquStr = 0;
                    }
                    eq2tautPrev = 0;
                    pINChI_Prev = pINChI;
                    pINChI_Taut_Prev = pINChI_Taut;
                    mult = 0;
                }
                else
                {
                    /*================ compare iso composition to previous =====================*/
                    /* check whether pINChI and pINChI_Prev have non-zero identical isotopic info */
                    eq2prev = bUseMulipliers && Eql_INChI_Isotopic( pINChI, pINChI_Prev );
                    if (eq2prev)
                    {
                        mult++; /* mult = (number of non-empty equal items)-1 */
                        continue;
                    }
                    else
                    {
                        /* pINChI isotopic info is either different or empty. Output pINChI_Prev anyway */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        if (pINChI_Prev && pINChI_Prev->nNumberOfAtoms)
                        {
                            if (( pINChI_Prev->nNumberOfIsotopicAtoms > 0 ||
                                  pINChI_Prev->nNumberOfIsotopicTGroups > 0 ))
                            {
                                /* pINChI_Prev exists and has isotopic info */
                                MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
                                /*  Isotopic atoms */
                                if (pINChI_Prev->nNumberOfIsotopicAtoms > 0 &&
                                    /*nStrLen-tot_len > 2 && */
                                    !*bOverflow)
                                {
                                    MakeIsoAtomString( pINChI_Prev->IsotopicAtom,
                                        pINChI_Prev->nNumberOfIsotopicAtoms,
                                        strbuf, TAUT_MODE, bOverflow );
                                }
                                /*  Isotopic tautomeric groups */
                                if (pINChI_Prev->nNumberOfIsotopicTGroups > 0 &&
                                    /*nStrLen-tot_len > 3 && */
                                    !*bOverflow)
                                {
                                    MakeDelim( bAbcNumbers ? ITEM_DELIMETER : "(", strbuf, bOverflow );
                                    MakeIsoTautString( pINChI_Prev->IsotopicTGroup, pINChI_Prev->nNumberOfIsotopicTGroups,
                                        strbuf, TAUT_MODE, bOverflow );
                                    if (!bAbcNumbers)
                                    {
                                        MakeDelim( ")", strbuf, bOverflow );
                                    }
                                }
                            }
                            /* else isotopic info is not present in pINChI_Prev */
                        }
                        else
                        {
                            if (bSecondNonTautPass && pINChI_Taut_Prev && pINChI_Taut_Prev->nNumberOfAtoms)
                            {
                                if (( pINChI_Taut_Prev->nNumberOfIsotopicAtoms > 0 ||
                                      pINChI_Taut_Prev->nNumberOfIsotopicTGroups > 0 ))
                                {
                                    /* since pINChI_Prev does not exist, pINChI_Taut_Prev is non-tautomeric */
                                    /* and it has non-trivial isotopic info */
                                    /*
                                    tot_len += MakeDelim( sIdenticalValues, strbuf, bOverflow);
                                    */
                                    ;/* pINChI_Taut_Prev isotopic info was output in the main isotopic section */
                                }
                                else
                                {
                                    ; /* pINChI_Taut_Prev exists and has not isotopic info */
                                }
                            }
    #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
                            else
                            {
                                int stop = 1;   /* <BRKPT> */
                            }
    #endif
                        }
                    }
                    /* Fix17: moved here 2004-10-08 */
                    pINChI_Prev = pINChI;
                    pINChI_Taut_Prev = pINChI_Taut;
                    mult = 0; /* we do not know whether the item is empty */
                }
            /* Fix17: moved from here 2004-10-08
            pINChI_Prev = pINChI;
            pINChI_Taut_Prev = pINChI_Taut;
            mult = 0;
            */
        }

        return ( strbuf->nUsedLength - nUsedLength0 );
    }
        */
    // END INCHI C FUNCTION: str_IsoAtoms

    fn clone_inchi(
        heap: &SourceHeap,
        pointer: SourceMutPointer<crate::source_types::INChI>,
    ) -> Result<Option<crate::source_types::INChI>, SourceHeapError> {
        if pointer.is_null() {
            Ok(None)
        } else {
            Ok(Some(
                heap.slice(pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            ))
        }
    }

    fn equivalence_pointer(heap: &mut SourceHeap, value: &str) -> Result<SourceConstPointer<i8>, SourceHeapError> {
        heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain(std::iter::once(0)).collect())
            .map(SourceMutPointer::as_const)
    }

    #[allow(clippy::too_many_arguments)]
    fn output_isotopic(
        heap: &mut SourceHeap,
        inchi: Option<&crate::source_types::INChI>,
        multiplier: i32,
        star: SourceConstPointer<i8>,
        open_group: SourceConstPointer<i8>,
        close_group: SourceConstPointer<i8>,
        abc_numbers: i32,
        string_buffer: &mut INCHI_IOS_STRING,
        taut_mode: i32,
        overflow: &mut i32,
    ) -> Result<(), SourceHeapError> {
        let Some(inchi) = inchi.filter(|inchi| inchi.nNumberOfAtoms != 0) else {
            return Ok(());
        };
        if inchi.nNumberOfIsotopicAtoms <= 0 && inchi.nNumberOfIsotopicTGroups <= 0 {
            return Ok(());
        }
        MakeMult(heap, multiplier, star, string_buffer, 0, overflow)?;
        if inchi.nNumberOfIsotopicAtoms > 0 && *overflow == 0 {
            MakeIsoAtomString(
                heap,
                inchi.IsotopicAtom.as_const(),
                inchi.nNumberOfIsotopicAtoms,
                string_buffer,
                taut_mode,
                overflow,
            )?;
        }
        if inchi.nNumberOfIsotopicTGroups > 0 && *overflow == 0 {
            MakeDelim(heap, open_group, string_buffer, overflow)?;
            MakeIsoTautString(
                heap,
                inchi.IsotopicTGroup.as_const(),
                inchi.nNumberOfIsotopicTGroups,
                string_buffer,
                taut_mode,
                overflow,
            )?;
            if abc_numbers == 0 {
                MakeDelim(heap, close_group, string_buffer, overflow)?;
            }
        }
        Ok(())
    }

    let initial_length = string_buffer.nUsedLength;
    let component_delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let star = heap.allocate_model_storage(vec![b'*' as i8, 0])?;
    let open_group = heap.allocate_model_storage(if abc_numbers != 0 {
        vec![b',' as i8, 0]
    } else {
        vec![b'(' as i8, 0]
    })?;
    let close_group = heap.allocate_model_storage(vec![b')' as i8, 0])?;
    let mut previous = SourceMutPointer::null();
    let mut taut_previous = SourceMutPointer::null();
    let mut multiplier = 0_i32;
    let mut next = 0_i32;
    let mut previous_equivalence: Option<&'static str> = None;
    let mut previous_equivalence_multiplier = 0_i32;
    let mut equal_to_taut_previous = 1_i32;
    let mut index = 0_i32;

    while index <= number_of_components {
        let current = if index < number_of_components {
            selected_inchi(heap, sorted_inchi.as_const().offset(i64::from(index))?, output_type)?
        } else {
            SourceMutPointer::null()
        };
        let taut_current = if second_non_taut_pass != 0 && index < number_of_components {
            selected_inchi(heap, sorted_inchi2.as_const().offset(i64::from(index))?, OUT_T1 as i32)?
        } else {
            SourceMutPointer::null()
        };
        let current_value = clone_inchi(heap, current)?;
        let taut_current_value = clone_inchi(heap, taut_current)?;
        let equal_to_taut = second_non_taut_pass != 0
            && omit_repetitions != 0
            && Eql_INChI_Isotopic(heap, current_value.as_ref(), taut_current_value.as_ref())? != 0;

        if equal_to_taut {
            let previous_value = clone_inchi(heap, previous)?;
            if previous_value.as_ref().is_some_and(|inchi| inchi.nNumberOfAtoms != 0) {
                if next != 0 {
                    MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                output_isotopic(
                    heap,
                    previous_value.as_ref(),
                    multiplier.wrapping_add(1),
                    star.as_const(),
                    open_group.as_const(),
                    close_group.as_const(),
                    abc_numbers,
                    string_buffer,
                    taut_mode,
                    overflow,
                )?;
            } else {
                let taut_previous_value = clone_inchi(heap, taut_previous)?;
                if taut_previous_value
                    .as_ref()
                    .is_some_and(|inchi| inchi.nNumberOfAtoms != 0)
                {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                }
            }
            let current_equivalence = crate::source::base::ichiprt1::EquString(
                (crate::source_types::iiNUMB | crate::source_types::iitNONTAUT) as i32,
            );
            if previous_equivalence_multiplier != 0 && previous_equivalence.is_some() {
                if previous_equivalence == Some(current_equivalence) {
                    previous_equivalence_multiplier = previous_equivalence_multiplier.wrapping_add(1);
                } else {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let pointer = equivalence_pointer(heap, previous_equivalence.unwrap())?;
                    MakeEqStr(heap, pointer, previous_equivalence_multiplier, string_buffer, overflow)?;
                    previous_equivalence = Some(current_equivalence);
                    previous_equivalence_multiplier = 1;
                }
            } else {
                previous_equivalence = Some(current_equivalence);
                previous_equivalence_multiplier = 1;
            }
            previous = SourceMutPointer::null();
            taut_previous = SourceMutPointer::null();
            multiplier = 0;
            equal_to_taut_previous = 1;
        } else if equal_to_taut_previous != 0 {
            if previous_equivalence_multiplier != 0 {
                if let Some(value) = previous_equivalence {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let pointer = equivalence_pointer(heap, value)?;
                    MakeEqStr(heap, pointer, previous_equivalence_multiplier, string_buffer, overflow)?;
                    previous_equivalence = None;
                    previous_equivalence_multiplier = 0;
                }
            }
            equal_to_taut_previous = 0;
            previous = current;
            taut_previous = taut_current;
            multiplier = 0;
        } else {
            let previous_value = clone_inchi(heap, previous)?;
            let equal_to_previous =
                use_multipliers != 0 && Eql_INChI_Isotopic(heap, current_value.as_ref(), previous_value.as_ref())? != 0;
            if equal_to_previous {
                multiplier = multiplier.wrapping_add(1);
                index = index.wrapping_add(1);
                continue;
            }
            if next != 0 {
                MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
            }
            next = next.wrapping_add(1);
            output_isotopic(
                heap,
                previous_value.as_ref(),
                multiplier.wrapping_add(1),
                star.as_const(),
                open_group.as_const(),
                close_group.as_const(),
                abc_numbers,
                string_buffer,
                taut_mode,
                overflow,
            )?;
            previous = current;
            taut_previous = taut_current;
            multiplier = 0;
        }
        index = index.wrapping_add(1);
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_IsoSp2(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    sorted_inchi2: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    taut_mode: i32,
    number_of_components: i32,
    second_non_taut_pass: i32,
    omit_repetitions: i32,
    use_multipliers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:1479 str_IsoSp2
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int str_IsoSp2( INCHI_SORT       *pINChISort,
                    INCHI_SORT       *pINChISort2,
                    INCHI_IOS_STRING *strbuf,
                    int              *bOverflow,
                    int              bOutType,
                    int              TAUT_MODE,
                    int              num_components,
                    int              bSecondNonTautPass,
                    int              bOmitRepetitions,
                    int              bUseMulipliers )
    {
        int          i, ii, ii2, nUsedLength0;
        INCHI_SORT   *is, *is2, *is0, *is20;
        INChI        *pINChI, *pINChI_Prev, *pINChI_Taut, *pINChI_Taut_Prev;
        INChI_Stereo *Stereo, *Stereo_Prev, *Stereo_Taut, *Stereo_Taut_Prev;
        int          mult, eq2prev, eq2taut, eq2tautPrev = 1, bNext; /* djb-rwth: initialisation required to avoid garbage values */
        const char  *pPrevEquStr, *pCurrEquStr;
        int         multPrevEquStr;
        pINChI_Taut = NULL;
        pINChI_Prev = NULL;
        pINChI_Taut_Prev = NULL;
        mult = 0;
        bNext = 0;
        is = NULL;
        is2 = NULL;
        is0 = pINChISort;
        is20 = bSecondNonTautPass ? pINChISort2 : NULL;
        /* djb-rwth: removing redundant code */
        pPrevEquStr = NULL; /*, *pCurrEquStr;*/
        multPrevEquStr = 0;
        nUsedLength0 = strbuf->nUsedLength;

        /* For each connected component...    */
        for (i = 0; i <= num_components; i++)
        {
            /* 1st (taut) pass: bOutType=OUT_TN  ; 2nd (non-taut pass) bOutType=OUT_NT */
            pINChI = ( i < num_components && ( is = is0 + i, 0 <= ( ii = GET_II( bOutType, is ) ) ) ) ? is->pINChI[ii] : NULL;
            /*================ compare sp2 to previous =====================*/
            if (bSecondNonTautPass)
            {
                /* component that was output on the 1st pass */
                pINChI_Taut = ( i < num_components && ( is2 = is20 + i, 0 <= ( ii2 = GET_II( OUT_T1, is2 ) ) ) ) ? is2->pINChI[ii2] : NULL;
            }
            eq2taut = 0;
            /*========= if bSecondNonTautPass then compare iso non-taut stereo to other stereo ========*/
            if (bSecondNonTautPass && bOmitRepetitions)
            {
                /* compare non-tautomeric isotopic to:
                *   a) non-tautomeric non-isotopic
                *   b) tautomeric non-isotopic
                *   c) tautomeric isotopic
                */
                /* a) compare non-tautomeric isotopic to non-tautomeric non-isotopic */
                if (!eq2taut)
                {
                    eq2taut = pINChI &&
                        /* non-taut isotopic */                  /* non-taut non-isotopic */
                        ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Taut = pINChI->Stereo ) &&
                        Eql_INChI_Stereo( Stereo, EQL_SP2, Stereo_Taut, EQL_SP2, 0 );
                    /* stereo     isotopic non-taut =  non-taut (stereo) */
                    eq2taut = eq2taut ? ( iiSTEREO | iitISO | iitNONTAUT | iiEq2NONTAUT ) : 0;
                }
                /* b) compare non-tautomeric isotopic to tautomeric non-isotopic */
                if (!eq2taut)
                {
                    eq2taut = pINChI && pINChI_Taut &&
                        /* non-taut isotopic */                  /* taut non-isotopic */
                        ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Taut = pINChI_Taut->Stereo ) &&
                        Eql_INChI_Stereo( Stereo, EQL_SP2, Stereo_Taut, EQL_SP2, 0 );
                    /* stereo     isotopic non-taut =  taut (stereo) */
                    eq2taut = eq2taut ? ( iiSTEREO | iitISO | iitNONTAUT ) : 0;
                }
                /* c) compare non-tautomeric isotopic to tautomeric isotopic */
                if (!eq2taut && bSecondNonTautPass && bOmitRepetitions)
                {
                    eq2taut = pINChI && pINChI_Taut &&
                        ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Taut = pINChI_Taut->StereoIsotopic ) &&
                        Eql_INChI_Stereo( Stereo, EQL_SP2, Stereo_Taut, EQL_SP2, 0 );
                    /* stereo     isotopic non-taut =  isotopic taut (stereo) */
                    eq2taut = eq2taut ? ( iiSTEREO | iitISO | iitNONTAUT | iiEq2ISO ) : 0;
                }
    #if ( FIX_EMPTY_LAYER_BUG == 1 )
                if (!eq2taut && pINChI && !( ( Stereo = pINChI->StereoIsotopic ) &&
                    Eql_INChI_Stereo( Stereo, EQL_SP2, NULL, EQL_EXISTS, 0 ) ))
                {
                    /* component has no stereo; check whether it has stereo in the preceding layers */
                    if (pINChI_Taut && ( Stereo_Taut = pINChI_Taut->Stereo ) && /* F is not empty */
                        Eql_INChI_Stereo( Stereo_Taut, EQL_SP2, NULL, EQL_EXISTS, 0 ) ||
                        !( pINChI_Taut && ( Stereo_Taut = pINChI_Taut->Stereo ) &&  /* M is empty and ... */
                            Eql_INChI_Stereo( Stereo_Taut, EQL_SP2, NULL, EQL_EXISTS, 0 ) ) &&
                            ( pINChI_Taut && ( Stereo_Taut = pINChI_Taut->StereoIsotopic ) &&  /* ... MI is not empty */
                                Eql_INChI_Stereo( Stereo_Taut, EQL_SP2, NULL, EQL_EXISTS, 0 ) ))
                    {

                        eq2taut = iiEmpty; /* the component has stereo in the preceding layer  */
                    }
                }
    #endif
            }
            else
            {
                /*========= if not bSecondNonTautPass then compare iso taut stereo to non-iso taut ========*/
                if (!bSecondNonTautPass && bOmitRepetitions)
                {
                    /* compare tautomeric isotopic to tautomeric non-isotopic */
                    if (!eq2taut)
                    {
                        eq2taut = pINChI &&
                            /* taut isotopic */                  /* taut non-isotopic */
                            ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Taut = pINChI->Stereo ) &&
                            Eql_INChI_Stereo( Stereo, EQL_SP2, Stereo_Taut, EQL_SP2, 0 );
                        /* stereo     isotopic taut =  taut (stereo) */
                        eq2taut = eq2taut ? ( iiSTEREO | iitISO ) : 0;
    #if ( FIX_EMPTY_LAYER_BUG == 1 )
                        if (!eq2taut && pINChI && !( ( Stereo = pINChI->StereoIsotopic ) &&
                                                     Eql_INChI_Stereo( Stereo, EQL_SP2, NULL, EQL_EXISTS, 0 ) ))
                        {
                            /* component has no MI stereo; check whether it has stereo in the preceding layer M */
                            if (( Stereo_Taut = pINChI->Stereo ) &&
                                Eql_INChI_Stereo( Stereo_Taut, EQL_SP2, NULL, EQL_EXISTS, 0 ))
                            {
                                eq2taut = iiEmpty; /* the component has stereo in the preceding layer  */
                            }
                        }
    #endif
                    }
                }
            }
            if (eq2taut)
            {
                /* we may be here only in case of the current layer found equal in another layer the same component */
                if (pINChI_Prev && pINChI_Prev->nNumberOfAtoms)
                {
                    /* previous component exists; output it before output the current component */
                    if (bNext++)
                    {
                        MakeDelim( sCompDelim, strbuf, bOverflow );
                    }
                    if (( Stereo_Prev = pINChI_Prev->StereoIsotopic ) && Stereo_Prev->nNumberOfStereoBonds > 0)
                    {
                        MakeMult( mult + 1, "*", strbuf, 0, bOverflow );

                        MakeStereoString( Stereo_Prev->nBondAtom1, Stereo_Prev->nBondAtom2,
                            Stereo_Prev->b_parity,
                            0, Stereo_Prev->nNumberOfStereoBonds,
                            strbuf, TAUT_MODE, bOverflow );
                    }
                }
                else
                {
                    if (pINChI_Taut_Prev && pINChI_Taut_Prev->nNumberOfAtoms)
                    {
                        /* previous non-taut component exists only in taut list */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        /* do not output stereo of non-tautomeric in non-taut layer: it has been output in the main layer */
                    }
                }
                /* we have found another (previously printed) layer of the current component equal to this layer */
                /* output this (current) equivalence mark = EquString(eq2taut) */
                pCurrEquStr = EquString( eq2taut );
                if (multPrevEquStr && pPrevEquStr)
                {
                    if (pCurrEquStr && !strcmp( pCurrEquStr, pPrevEquStr ))
                    {
                        multPrevEquStr++;
                    }
                    else
                    {
                        /* new EqStr is different; output it */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                        pPrevEquStr = pCurrEquStr;
                        multPrevEquStr = 1;
                    }
                }
                else
                {
                    pPrevEquStr = pCurrEquStr;
                    multPrevEquStr = 1;
                }
                pINChI_Prev = NULL; /* pINChI_Prev sp2 does not exist since */
                pINChI_Taut_Prev = NULL; /* pINChI has just been printed */
                mult = 0;
                eq2tautPrev = 1;     /* pINChI_Prev and pINChI_Taut_Prev have already been */
            }
            else
            {
                if (eq2tautPrev)
                {
                    /* at this point pINChI_Prev does not exist; however, pINChI */
                    /*might have been discovered and it is different from pINChI_Taut */
                    if (multPrevEquStr && pPrevEquStr)
                    {
                        /* new EqStr is different; output it */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                        pPrevEquStr = NULL;
                        multPrevEquStr = 0;
                    }
                    eq2tautPrev = 0;
                    pINChI_Prev = pINChI;
                    pINChI_Taut_Prev = pINChI_Taut;
                    mult = 0;
                }
                else
                {
                    /* current layer is different from previously printed layers of the current component */
                    /* compare the current layer to this layer of the previous component: */
                    /* check whether pINChI and pINChI_Prev have non-zero identical stereo sp2 */
                    /*================ compare iso sp2 to previous =====================*/
                    eq2prev = bUseMulipliers &&
                        pINChI && pINChI_Prev &&
                        ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Prev = pINChI_Prev->StereoIsotopic ) &&
                        Eql_INChI_Stereo( Stereo, EQL_SP2, Stereo_Prev, EQL_SP2, 0 );
                    if (eq2prev)
                    {
                        mult++; /* mult = (number of non-empty equal items)-1 */
                        continue;
                    }
                    else
                    {
                        /* the current layer is different from this layer of the previous component */
                        /* therefore print the current layer */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        if (pINChI_Prev && pINChI_Prev->nNumberOfAtoms)
                        {
                            if (( Stereo_Prev = pINChI_Prev->StereoIsotopic ) && Stereo_Prev->nNumberOfStereoBonds > 0)
                            {
                                MakeMult( mult + 1, "*", strbuf, 0, bOverflow );

                                MakeStereoString( Stereo_Prev->nBondAtom1, Stereo_Prev->nBondAtom2,
                                    Stereo_Prev->b_parity,
                                    0, Stereo_Prev->nNumberOfStereoBonds,
                                    strbuf, TAUT_MODE, bOverflow );
                            }
                            /* else sp2 info is not present in pINChI_Prev */
                        }
                        else
                        {/* do not print pINChI_Prev because it either do not exist of have already been printed */
                            if (bSecondNonTautPass && pINChI_Taut_Prev && pINChI_Taut_Prev->nNumberOfAtoms)
                            {
                                if (( Stereo_Taut_Prev = pINChI_Taut_Prev->StereoIsotopic ) && Stereo_Taut_Prev->nNumberOfStereoBonds > 0)
                                {
                                    /* since pINChI_Prev does not exist, pINChI_Taut_Prev is non-tautomeric */
                                    /* and it has non-trivial sp2 info */
                                    /*
                                    tot_len += MakeDelim( sIdenticalValues, strbuf, bOverflow);
                                    */
                                    ;/* pINChI_Taut_Prev sp3 info was output in the main stereo section */
                                }
                                else
                                {
                                    ; /* pINChI_Taut_Prev exists and has not sp2 info */
                                }
                            }
    #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
                            else
                            {
                                int stop = 1;   /* <BRKPT> */
                            }
    #endif
                        }
                    }
                    pINChI_Prev = pINChI;
                    pINChI_Taut_Prev = pINChI_Taut;
                    mult = 0; /* we do not know whether the item is empty */
                }
            }
        }

        return ( strbuf->nUsedLength - nUsedLength0 );
    }
        */
    // END INCHI C FUNCTION: str_IsoSp2

    fn clone_inchi(
        heap: &SourceHeap,
        pointer: SourceMutPointer<crate::source_types::INChI>,
    ) -> Result<Option<crate::source_types::INChI>, SourceHeapError> {
        if pointer.is_null() {
            Ok(None)
        } else {
            Ok(Some(
                heap.slice(pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            ))
        }
    }

    fn clone_stereo(
        heap: &SourceHeap,
        inchi: Option<&crate::source_types::INChI>,
        isotopic: bool,
    ) -> Result<Option<crate::source_types::INChI_Stereo>, SourceHeapError> {
        let Some(inchi) = inchi else {
            return Ok(None);
        };
        let pointer = if isotopic { inchi.StereoIsotopic } else { inchi.Stereo };
        if pointer.is_null() {
            Ok(None)
        } else {
            Ok(Some(
                heap.slice(pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            ))
        }
    }

    fn equivalence_pointer(heap: &mut SourceHeap, value: &str) -> Result<SourceConstPointer<i8>, SourceHeapError> {
        heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain(std::iter::once(0)).collect())
            .map(SourceMutPointer::as_const)
    }

    fn stereo_equal(
        heap: &SourceHeap,
        left: Option<&crate::source_types::INChI_Stereo>,
        right: Option<&crate::source_types::INChI_Stereo>,
    ) -> Result<bool, SourceHeapError> {
        Ok(left.is_some()
            && right.is_some()
            && Eql_INChI_Stereo(
                heap,
                left,
                crate::source_types::EQL_SP2 as i32,
                right,
                crate::source_types::EQL_SP2 as i32,
                0,
            )? != 0)
    }

    fn output_stereo(
        heap: &mut SourceHeap,
        inchi: Option<&crate::source_types::INChI>,
        multiplier: i32,
        star: SourceConstPointer<i8>,
        string_buffer: &mut INCHI_IOS_STRING,
        taut_mode: i32,
        overflow: &mut i32,
    ) -> Result<(), SourceHeapError> {
        let Some(inchi) = inchi.filter(|inchi| inchi.nNumberOfAtoms != 0) else {
            return Ok(());
        };
        let Some(stereo) = clone_stereo(heap, Some(inchi), true)? else {
            return Ok(());
        };
        if stereo.nNumberOfStereoBonds > 0 {
            MakeMult(heap, multiplier, star, string_buffer, 0, overflow)?;
            MakeStereoString(
                heap,
                stereo.nBondAtom1.as_const(),
                stereo.nBondAtom2.as_const(),
                stereo.b_parity.as_const(),
                0,
                stereo.nNumberOfStereoBonds,
                string_buffer,
                taut_mode,
                overflow,
            )?;
        }
        Ok(())
    }

    let initial_length = string_buffer.nUsedLength;
    let component_delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let star = heap.allocate_model_storage(vec![b'*' as i8, 0])?;
    let mut previous = SourceMutPointer::null();
    let mut taut_previous = SourceMutPointer::null();
    let mut multiplier = 0_i32;
    let mut next = 0_i32;
    let mut previous_equivalence: Option<&'static str> = None;
    let mut previous_equivalence_multiplier = 0_i32;
    let mut equal_to_taut_previous = 1_i32;
    let mut index = 0_i32;

    while index <= number_of_components {
        let current = if index < number_of_components {
            selected_inchi(heap, sorted_inchi.as_const().offset(i64::from(index))?, output_type)?
        } else {
            SourceMutPointer::null()
        };
        let taut_current = if second_non_taut_pass != 0 && index < number_of_components {
            selected_inchi(heap, sorted_inchi2.as_const().offset(i64::from(index))?, OUT_T1 as i32)?
        } else {
            SourceMutPointer::null()
        };
        let current_value = clone_inchi(heap, current)?;
        let taut_current_value = clone_inchi(heap, taut_current)?;
        let current_iso = clone_stereo(heap, current_value.as_ref(), true)?;
        let current_main = clone_stereo(heap, current_value.as_ref(), false)?;
        let taut_main = clone_stereo(heap, taut_current_value.as_ref(), false)?;
        let taut_iso = clone_stereo(heap, taut_current_value.as_ref(), true)?;

        let mut equivalence = 0_i32;
        if second_non_taut_pass != 0 && omit_repetitions != 0 {
            if stereo_equal(heap, current_iso.as_ref(), current_main.as_ref())? {
                equivalence = (crate::source_types::iiSTEREO
                    | crate::source_types::iitISO
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2NONTAUT) as i32;
            } else if stereo_equal(heap, current_iso.as_ref(), taut_main.as_ref())? {
                equivalence = (crate::source_types::iiSTEREO
                    | crate::source_types::iitISO
                    | crate::source_types::iitNONTAUT) as i32;
            } else if stereo_equal(heap, current_iso.as_ref(), taut_iso.as_ref())? {
                equivalence = (crate::source_types::iiSTEREO
                    | crate::source_types::iitISO
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2ISO) as i32;
            }
        } else if second_non_taut_pass == 0
            && omit_repetitions != 0
            && stereo_equal(heap, current_iso.as_ref(), current_main.as_ref())?
        {
            equivalence = (crate::source_types::iiSTEREO | crate::source_types::iitISO) as i32;
        }

        if equivalence != 0 {
            let previous_value = clone_inchi(heap, previous)?;
            if previous_value.as_ref().is_some_and(|inchi| inchi.nNumberOfAtoms != 0) {
                if next != 0 {
                    MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                output_stereo(
                    heap,
                    previous_value.as_ref(),
                    multiplier.wrapping_add(1),
                    star.as_const(),
                    string_buffer,
                    taut_mode,
                    overflow,
                )?;
            } else {
                let taut_previous_value = clone_inchi(heap, taut_previous)?;
                if taut_previous_value
                    .as_ref()
                    .is_some_and(|inchi| inchi.nNumberOfAtoms != 0)
                {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                }
            }
            let current_equivalence = crate::source::base::ichiprt1::EquString(equivalence);
            if previous_equivalence_multiplier != 0 && previous_equivalence.is_some() {
                if previous_equivalence == Some(current_equivalence) {
                    previous_equivalence_multiplier = previous_equivalence_multiplier.wrapping_add(1);
                } else {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let pointer = equivalence_pointer(heap, previous_equivalence.unwrap())?;
                    MakeEqStr(heap, pointer, previous_equivalence_multiplier, string_buffer, overflow)?;
                    previous_equivalence = Some(current_equivalence);
                    previous_equivalence_multiplier = 1;
                }
            } else {
                previous_equivalence = Some(current_equivalence);
                previous_equivalence_multiplier = 1;
            }
            previous = SourceMutPointer::null();
            taut_previous = SourceMutPointer::null();
            multiplier = 0;
            equal_to_taut_previous = 1;
        } else if equal_to_taut_previous != 0 {
            if previous_equivalence_multiplier != 0 {
                if let Some(value) = previous_equivalence {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let pointer = equivalence_pointer(heap, value)?;
                    MakeEqStr(heap, pointer, previous_equivalence_multiplier, string_buffer, overflow)?;
                    previous_equivalence = None;
                    previous_equivalence_multiplier = 0;
                }
            }
            equal_to_taut_previous = 0;
            previous = current;
            taut_previous = taut_current;
            multiplier = 0;
        } else {
            let previous_value = clone_inchi(heap, previous)?;
            let previous_iso = clone_stereo(heap, previous_value.as_ref(), true)?;
            let equal_to_previous =
                use_multipliers != 0 && stereo_equal(heap, current_iso.as_ref(), previous_iso.as_ref())?;
            if equal_to_previous {
                multiplier = multiplier.wrapping_add(1);
                index = index.wrapping_add(1);
                continue;
            }
            if next != 0 {
                MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
            }
            next = next.wrapping_add(1);
            output_stereo(
                heap,
                previous_value.as_ref(),
                multiplier.wrapping_add(1),
                star.as_const(),
                string_buffer,
                taut_mode,
                overflow,
            )?;
            previous = current;
            taut_previous = taut_current;
            multiplier = 0;
        }
        index = index.wrapping_add(1);
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_AuxTgroupEqu(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    taut_mode: i32,
    number_of_components: i32,
    use_multipliers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:3953 str_AuxTgroupEqu
    // INCHI✔️❌: int str_AuxTgroupEqu( INCHI_SORT       *pINChISort,
    // INCHI✔️❌:                       INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                       int              *bOverflow,
    // INCHI✔️❌:                       int              bOutType,
    // INCHI✔️❌:                       int              TAUT_MODE,
    // INCHI✔️❌:                       int              num_components,
    // INCHI✔️❌:                       int              bUseMulipliers )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int          i, ii, nUsedLength0;
    // INCHI✔️❌:     INCHI_SORT   *is, *is0;
    // INCHI✔️❌:     INChI_Aux    *pINChI_Aux, *pINChI_Aux_Prev;
    // INCHI✔️❌:     int          mult, eq2prev, bNext;
    // INCHI✔️❌:
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     is0 = pINChISort;
    // INCHI✔️❌:     is = NULL;
    // INCHI✔️❌:     i = 0;
    // INCHI✔️❌:     pINChI_Aux_Prev = ( 0 <= ( ii = GET_II( bOutType, is0 ) ) ) ? is0->pINChI_Aux[ii] : NULL;
    // INCHI✔️❌:     mult = 0;
    // INCHI✔️❌:     bNext = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* For each connected component...    */
    // INCHI✔️❌:     for (i++; i <= num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         pINChI_Aux = ( i < num_components && ( is = is0 + i, 0 <= ( ii = GET_II( bOutType, is ) ) ) ) ? is->pINChI_Aux[ii] : NULL;
    // INCHI✔️❌:         eq2prev = bUseMulipliers &&
    // INCHI✔️❌:             Eql_INChI_Aux_Equ( pINChI_Aux, EQL_EQU_TG, pINChI_Aux_Prev, EQL_EQU_TG );
    // INCHI✔️❌:         if (eq2prev)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             mult++; /* mult = (number of non-empty equal items)-1 */
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (bNext++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (pINChI_Aux_Prev && pINChI_Aux_Prev->nNumberOfTGroups &&
    // INCHI✔️❌:                 bHasEquString( pINChI_Aux_Prev->nConstitEquTGroupNumbers, pINChI_Aux_Prev->nNumberOfTGroups ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
    // INCHI✔️❌:                 MakeEquString( pINChI_Aux_Prev->nConstitEquTGroupNumbers,
    // INCHI✔️❌:                     pINChI_Aux_Prev->nNumberOfTGroups,
    // INCHI✔️❌:                     0, strbuf,
    // INCHI✔️❌:                     TAUT_MODE, bOverflow );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         pINChI_Aux_Prev = pINChI_Aux;
    // INCHI✔️❌:         mult = 0; /* we do not know whether the item is empty */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( strbuf->nUsedLength - nUsedLength0 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: str_AuxTgroupEqu

    let initial_length = string_buffer.nUsedLength;
    if sorted_inchi.is_null() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let aux_at = |heap: &SourceHeap, index: i32| -> Result<Option<crate::source_types::INChI_Aux>, SourceHeapError> {
        let sort = heap
            .slice(sorted_inchi.as_const().offset(i64::from(index))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let Some(selected) = get_ii(heap, output_type, sort)? else {
            return Ok(None);
        };
        let pointer = sort.pINChI_Aux[selected];
        if pointer.is_null() {
            Ok(None)
        } else {
            Ok(heap.slice(pointer.as_const())?.first().cloned())
        }
    };
    let mut previous = aux_at(&*heap, 0)?;
    let delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let star = heap.allocate_model_storage(vec![b'*' as i8, 0])?;
    let mut multiplier = 0_i32;
    let mut next = 0_i32;
    for index in 1..=number_of_components {
        let current = if index < number_of_components {
            aux_at(&*heap, index)?
        } else {
            None
        };
        if use_multipliers != 0
            && Eql_INChI_Aux_Equ(
                &*heap,
                current.as_ref(),
                crate::source_types::EQL_EQU_TG as i32,
                previous.as_ref(),
                crate::source_types::EQL_EQU_TG as i32,
            )? != 0
        {
            multiplier = multiplier.wrapping_add(1);
            continue;
        }
        if next != 0 {
            MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
        }
        next = next.wrapping_add(1);
        if let Some(aux) = previous.as_ref() {
            if aux.nNumberOfTGroups != 0
                && bHasEquString(&*heap, aux.nConstitEquTGroupNumbers.as_const(), aux.nNumberOfTGroups)? != 0
            {
                MakeMult(
                    heap,
                    multiplier.wrapping_add(1),
                    star.as_const(),
                    string_buffer,
                    0,
                    overflow,
                )?;
                MakeEquString(
                    heap,
                    aux.nConstitEquTGroupNumbers.as_const(),
                    aux.nNumberOfTGroups,
                    0,
                    string_buffer,
                    taut_mode,
                    overflow,
                )?;
            }
        }
        previous = current;
        multiplier = 0;
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_AuxChargeRadVal(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    taut_mode: i32,
    number_of_components: i32,
    use_multipliers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:4013 str_AuxChargeRadVal
    // INCHI✔️❌: int str_AuxChargeRadVal( INCHI_SORT       *pINChISort,
    // INCHI✔️❌:     INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:     int               *bOverflow,
    // INCHI✔️❌:     int               bOutType,
    // INCHI✔️❌:     int               TAUT_MODE,
    // INCHI✔️❌:     int               num_components,
    // INCHI✔️❌:     int               bUseMulipliers )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int          i, ii, nUsedLength0;
    // INCHI✔️❌:     INCHI_SORT   *is, *is0;
    // INCHI✔️❌:     INChI_Aux    *pINChI_Aux, *pINChI_Aux_Prev;
    // INCHI✔️❌:     int          mult, eq2prev, bNext;
    // INCHI✔️❌:
    // INCHI✔️❌:     pINChI_Aux_Prev = NULL;
    // INCHI✔️❌:     mult = 0;
    // INCHI✔️❌:     bNext = 0;
    // INCHI✔️❌:     is = NULL;
    // INCHI✔️❌:     is0 = pINChISort;
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* For each connected component...    */
    // INCHI✔️❌:     for (i = 0; i <= num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         /* 1st (taut) pass: bOutType=OUT_TN  ; 2nd (non-taut pass) bOutType=OUT_NT */
    // INCHI✔️❌:         pINChI_Aux = ( i < num_components && ( is = is0 + i, 0 <= ( ii = GET_II( bOutType, is ) ) ) ) ? is->pINChI_Aux[ii] : NULL;
    // INCHI✔️❌:         /* check whether pINChI_Aux and pINChI_Aux_Prev have identical info */
    // INCHI✔️❌:         eq2prev = bUseMulipliers &&
    // INCHI✔️❌:             EqlOrigInfo( pINChI_Aux, pINChI_Aux_Prev );
    // INCHI✔️❌:         if (eq2prev)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* eq. info is same and non-trivial */
    // INCHI✔️❌:             mult++; /* mult = (number of non-empty equal items)-1 */
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* pINChI_Aux info is either different or trivial. Output pINChI_Aux_Prev anyway */
    // INCHI✔️❌:                 if (bNext++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (pINChI_Aux_Prev && pINChI_Aux_Prev->nNumberOfAtoms)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (bHasOrigInfo( pINChI_Aux_Prev->OrigInfo, pINChI_Aux_Prev->nNumberOfAtoms ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* pINChI_Aux_Prev exists and has orig. info info */
    // INCHI✔️❌:                         MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
    // INCHI✔️❌:                         MakeCRVString( pINChI_Aux_Prev->OrigInfo,
    // INCHI✔️❌:                             pINChI_Aux_Prev->nNumberOfAtoms,
    // INCHI✔️❌:                             0, strbuf, TAUT_MODE, bOverflow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ; /* pINChI_Aux_Prev exists and has only trivial info */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI❌❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI❌❌:                 else
    // INCHI❌❌:                 {
    // INCHI❌❌:                     int stop = 1;   /* <BRKPT> */
    // INCHI❌❌:                 }
    // INCHI❌❌: #endif
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         pINChI_Aux_Prev = pINChI_Aux;
    // INCHI✔️❌:         mult = 0; /* we do not know whether the item is empty */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( strbuf->nUsedLength - nUsedLength0 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: str_AuxChargeRadVal

    let initial_length = string_buffer.nUsedLength;
    let aux_at = |heap: &SourceHeap, index: i32| -> Result<Option<crate::source_types::INChI_Aux>, SourceHeapError> {
        let sort = heap
            .slice(sorted_inchi.as_const().offset(i64::from(index))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let Some(selected) = get_ii(heap, output_type, sort)? else {
            return Ok(None);
        };
        let pointer = sort.pINChI_Aux[selected];
        if pointer.is_null() {
            Ok(None)
        } else {
            Ok(heap.slice(pointer.as_const())?.first().cloned())
        }
    };
    let delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let star = heap.allocate_model_storage(vec![b'*' as i8, 0])?;
    let mut previous: Option<crate::source_types::INChI_Aux> = None;
    let mut multiplier = 0_i32;
    let mut next = 0_i32;
    let mut index = 0_i32;
    while index <= number_of_components {
        let current = if index < number_of_components {
            aux_at(&*heap, index)?
        } else {
            None
        };
        if use_multipliers != 0 && EqlOrigInfo(&*heap, current.as_ref(), previous.as_ref())? != 0 {
            multiplier = multiplier.wrapping_add(1);
            index = index.wrapping_add(1);
            continue;
        }
        if index != 0 {
            if next != 0 {
                MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
            }
            next = next.wrapping_add(1);
            if let Some(aux) = previous.as_ref() {
                if aux.nNumberOfAtoms != 0 && bHasOrigInfo(&*heap, aux.OrigInfo.as_const(), aux.nNumberOfAtoms)? != 0 {
                    MakeMult(
                        heap,
                        multiplier.wrapping_add(1),
                        star.as_const(),
                        string_buffer,
                        0,
                        overflow,
                    )?;
                    MakeCRVString(
                        heap,
                        aux.OrigInfo.as_const(),
                        aux.nNumberOfAtoms,
                        0,
                        string_buffer,
                        taut_mode,
                        overflow,
                    )?;
                }
            }
        }
        previous = current;
        multiplier = 0;
        index = index.wrapping_add(1);
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_AuxIsoTgroupEqu(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    taut_mode: i32,
    number_of_components: i32,
    omit_repetitions: i32,
    use_multipliers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:4236 str_AuxIsoTgroupEqu
    // INCHI✔️❌: complete configured source frame follows verbatim; the non-release debug block is inactive.
    /*
    int str_AuxIsoTgroupEqu( INCHI_SORT       *pINChISort,
        INCHI_IOS_STRING *strbuf,
        int              *bOverflow,
        int              bOutType,
        int              TAUT_MODE,
        int              num_components,
        int              bOmitRepetitions,
        int              bUseMulipliers )
    {
        int          i, ii, nUsedLength0;
        INCHI_SORT   *is, *is0;
        INChI_Aux    *pINChI_Aux, *pINChI_Aux_Prev;
        int          mult, eq2prev, eq2taut, eq2tautPrev, bNext;
        const char  *pPrevEquStr, *pCurrEquStr;
        int         multPrevEquStr;
        pINChI_Aux = NULL;
        pINChI_Aux_Prev = NULL;
        mult = 0;
        bNext = 0;
        is = NULL;
        is0 = pINChISort;
        /* djb-rwth: removing redundant code */
        eq2tautPrev = 1; /* pINChI_Aux_Prev (previous pINChI_Aux) does not exist */
        pPrevEquStr = NULL; /*, *pCurrEquStr;*/
        multPrevEquStr = 0;
        nUsedLength0 = strbuf->nUsedLength;

        /* For each connected component...    */
        for (i = 0; i <= num_components; i++)
        {
            /* 1st (taut) pass: bOutType=OUT_TN  ; 2nd (non-taut pass) bOutType=OUT_NT */
            pINChI_Aux = ( i < num_components && ( is = is0 + i, 0 <= ( ii = GET_II( bOutType, is ) ) ) ) ? is->pINChI_Aux[ii] : NULL;
            /*================ compare iso non-taut equivalence info to non-iso taut ========*/
            eq2taut = 0;
            if (bOmitRepetitions && pINChI_Aux && pINChI_Aux->bIsIsotopic)
            {
                /**************************************************
                * compare isotopic tautomeric equivalence to:
                *    a) non-isotopic tautomeric
                */
                /* compare isotopic t-group equivalence to non-isotopic */
                eq2taut = Eql_INChI_Aux_Equ( pINChI_Aux, EQL_EQU_TG | EQL_EQU_ISO, pINChI_Aux, EQL_EQU_TG );
                /* equ   taut-isotopic = tautomeric, same as for isotopic atom equivalence info*/
                eq2taut = eq2taut ? ( iiEQU | iitISO ) : 0;
            }
            if (eq2taut)
            {
                /* current isotopic t-group equivalence has been found to be same as non-isotopic */
                if (pINChI_Aux_Prev && pINChI_Aux_Prev->nNumberOfAtoms)
                {
                    /* previous component exists */
                    if (bNext++)
                    {
                        MakeDelim( sCompDelim, strbuf, bOverflow );
                    }
                    if (bHasEquString( pINChI_Aux_Prev->nConstitEquIsotopicTGroupNumbers, pINChI_Aux_Prev->nNumberOfTGroups ))
                    {
                        /* output previous component(s) equivalence since it was found to be non-trivial */
                        MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
                        MakeEquString( pINChI_Aux_Prev->nConstitEquIsotopicTGroupNumbers, pINChI_Aux_Prev->nNumberOfTGroups, 0,
                            strbuf, TAUT_MODE, bOverflow );
                    }
                    else
                    {
                        ; /* pINChI_Aux_Prev exists and does not have non-trivial t-group equivalence info */
                    }
                }
                /* we have found pINChI_Aux->pINChI_Aux->nConstitEquIsotopicTGroupNumbers same as in pINChI_Aux->nConstitEquTGroupNumbers */
                pCurrEquStr = EquString( eq2taut );
                if (multPrevEquStr && pPrevEquStr)
                {
                    if (pCurrEquStr && !strcmp( pCurrEquStr, pPrevEquStr ))
                    {
                        multPrevEquStr++;
                    }
                    else
                    {
                        /* new EqStr is different; output it */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                        pPrevEquStr = pCurrEquStr;
                        multPrevEquStr = 1;
                    }
                }
                else
                {
                    pPrevEquStr = pCurrEquStr;
                    multPrevEquStr = 1;
                }
                pINChI_Aux_Prev = NULL; /* pINChI_Aux_Prev has already been output */
                mult = 0;
                eq2tautPrev = 1;
            }
            else
            {
                if (eq2tautPrev)
                {
                    /* at this point pINChI_Aux_Prev does not exist; however, pINChI_Aux */
                    /* might have been discovered and it may be different from non-isotopic */
                    if (multPrevEquStr && pPrevEquStr)
                    {
                        /* new EqStr is different; output it */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                        pPrevEquStr = NULL;
                        multPrevEquStr = 0;
                    }
                    eq2tautPrev = 0;
                    pINChI_Aux_Prev = pINChI_Aux;
                    mult = 0;
                }
                else
                {
                    /* check whether pINChI_Aux and pINChI_Aux_Prev have identical non-trivial isotopic t-group equivalence info */
                    eq2prev = bUseMulipliers && Eql_INChI_Aux_Equ( pINChI_Aux, EQL_EQU_TG | EQL_EQU_ISO, pINChI_Aux_Prev, EQL_EQU_TG | EQL_EQU_ISO );
                    if (eq2prev)
                    {
                        /* eq. info is same and non-trivial */
                        mult++; /* mult = (number of non-empty equal items)-1 */
                        continue;
                    }
                    else
                    {
                        /* pINChI_Aux eq. info is either different or trivial. Output pINChI_Aux_Prev anyway */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        if (pINChI_Aux_Prev && pINChI_Aux_Prev->nNumberOfAtoms)
                        {
                            if (bHasEquString( pINChI_Aux_Prev->nConstitEquIsotopicTGroupNumbers, pINChI_Aux_Prev->nNumberOfTGroups ))
                            {
                                /* pINChI_Aux_Prev exists and has equivalence info */
                                MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
                                MakeEquString( pINChI_Aux_Prev->nConstitEquIsotopicTGroupNumbers, pINChI_Aux_Prev->nNumberOfTGroups, 0,
                                    strbuf, TAUT_MODE, bOverflow );
                            }
                            else
                            {
                                ; /* pINChI_Aux_Prev exists and has only trivial equivalence info */
                            }
                        }
    #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
                        else
                        {
                            int stop = 1;   /* <BRKPT> */
                        }
    #endif
                    }
                    pINChI_Aux_Prev = pINChI_Aux;
                    mult = 0; /* we do not know whether the item is empty */
                }
            }
        }

        return ( strbuf->nUsedLength - nUsedLength0 );
    }
    */
    // END INCHI C FUNCTION: str_AuxIsoTgroupEqu

    let initial_length = string_buffer.nUsedLength;
    let delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let star = heap.allocate_model_storage(vec![b'*' as i8, 0])?;
    let mut previous_aux: Option<crate::source_types::INChI_Aux> = None;
    let mut previous_equivalence: Option<&'static str> = None;
    let mut previous_equivalence_multiplier = 0_i32;
    let mut multiplier = 0_i32;
    let mut next = 0_i32;
    let mut equal_to_taut_previous = true;
    let equivalence_mode = (crate::source_types::EQL_EQU_TG | crate::source_types::EQL_EQU_ISO) as i32;

    for component in 0..=number_of_components {
        let current_aux = if component < number_of_components {
            let sort = heap
                .slice(sorted_inchi.as_const().offset(i64::from(component))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let selected = get_ii(heap, output_type, sort)?;
            let pointer = selected.map_or(SourceMutPointer::null(), |slot| sort.pINChI_Aux[slot]);
            if pointer.is_null() {
                None
            } else {
                heap.slice(pointer.as_const())?.first().cloned()
            }
        } else {
            None
        };
        let mut equivalence = 0_i32;
        if omit_repetitions != 0
            && current_aux.as_ref().is_some_and(|aux| aux.bIsIsotopic != 0)
            && Eql_INChI_Aux_Equ(
                heap,
                current_aux.as_ref(),
                equivalence_mode,
                current_aux.as_ref(),
                crate::source_types::EQL_EQU_TG as i32,
            )? != 0
        {
            equivalence = (crate::source_types::iiEQU | crate::source_types::iitISO) as i32;
        }

        if equivalence != 0 {
            if let Some(aux) = previous_aux.as_ref()
                && aux.nNumberOfAtoms != 0
            {
                if next != 0 {
                    MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                if bHasEquString(
                    heap,
                    aux.nConstitEquIsotopicTGroupNumbers.as_const(),
                    aux.nNumberOfTGroups,
                )? != 0
                {
                    MakeMult(
                        heap,
                        multiplier.wrapping_add(1),
                        star.as_const(),
                        string_buffer,
                        0,
                        overflow,
                    )?;
                    MakeEquString(
                        heap,
                        aux.nConstitEquIsotopicTGroupNumbers.as_const(),
                        aux.nNumberOfTGroups,
                        0,
                        string_buffer,
                        taut_mode,
                        overflow,
                    )?;
                }
            }
            let current_equivalence = crate::source::base::ichiprt1::EquString(equivalence);
            if previous_equivalence_multiplier != 0 {
                if previous_equivalence == Some(current_equivalence) {
                    previous_equivalence_multiplier = previous_equivalence_multiplier.wrapping_add(1);
                } else {
                    if next != 0 {
                        MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let value = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let value =
                        heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain([0]).collect())?;
                    MakeEqStr(
                        heap,
                        value.as_const(),
                        previous_equivalence_multiplier,
                        string_buffer,
                        overflow,
                    )?;
                    previous_equivalence = Some(current_equivalence);
                    previous_equivalence_multiplier = 1;
                }
            } else {
                previous_equivalence = Some(current_equivalence);
                previous_equivalence_multiplier = 1;
            }
            previous_aux = None;
            multiplier = 0;
            equal_to_taut_previous = true;
        } else if equal_to_taut_previous {
            if previous_equivalence_multiplier != 0 {
                if next != 0 {
                    MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                let value = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
                let value = heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain([0]).collect())?;
                MakeEqStr(
                    heap,
                    value.as_const(),
                    previous_equivalence_multiplier,
                    string_buffer,
                    overflow,
                )?;
                previous_equivalence = None;
                previous_equivalence_multiplier = 0;
            }
            equal_to_taut_previous = false;
            previous_aux = current_aux;
            multiplier = 0;
        } else {
            let equal_to_previous = use_multipliers != 0
                && Eql_INChI_Aux_Equ(
                    heap,
                    current_aux.as_ref(),
                    equivalence_mode,
                    previous_aux.as_ref(),
                    equivalence_mode,
                )? != 0;
            if equal_to_previous {
                multiplier = multiplier.wrapping_add(1);
                continue;
            }
            if next != 0 {
                MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
            }
            next = next.wrapping_add(1);
            if let Some(aux) = previous_aux.as_ref()
                && aux.nNumberOfAtoms != 0
                && bHasEquString(
                    heap,
                    aux.nConstitEquIsotopicTGroupNumbers.as_const(),
                    aux.nNumberOfTGroups,
                )? != 0
            {
                MakeMult(
                    heap,
                    multiplier.wrapping_add(1),
                    star.as_const(),
                    string_buffer,
                    0,
                    overflow,
                )?;
                MakeEquString(
                    heap,
                    aux.nConstitEquIsotopicTGroupNumbers.as_const(),
                    aux.nNumberOfTGroups,
                    0,
                    string_buffer,
                    taut_mode,
                    overflow,
                )?;
            }
            previous_aux = current_aux;
            multiplier = 0;
        }
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_AuxEqu(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    sorted_inchi2: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    taut_mode: i32,
    number_of_components: i32,
    second_non_taut_pass: i32,
    omit_repetitions: i32,
    use_multipliers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:2106 str_AuxEqu
    // INCHI✔️❌: int str_AuxEqu( INCHI_SORT       *pINChISort,
    // INCHI✔️❌:                 INCHI_SORT       *pINChISort2,
    // INCHI✔️❌:                 INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                 int              *bOverflow,
    // INCHI✔️❌:                 int              bOutType,
    // INCHI✔️❌:                 int              TAUT_MODE,
    // INCHI✔️❌:                 int              num_components,
    // INCHI✔️❌:                 int              bSecondNonTautPass,
    // INCHI✔️❌:                 int              bOmitRepetitions,
    // INCHI✔️❌:                 int              bUseMulipliers )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int          i, ii, ii2, nUsedLength0;
    // INCHI✔️❌:     INCHI_SORT   *is, *is2, *is0, *is20;
    // INCHI✔️❌:     INChI_Aux    *pINChI_Aux = NULL, *pINChI_Aux_Prev, *pINChI_Aux_Taut, *pINChI_Aux_Taut_Prev;
    // INCHI✔️❌:     int          mult, eq2prev, eq2taut, eq2tautPrev, bNext;
    // INCHI✔️❌:     const char  *pPrevEquStr, *pCurrEquStr;
    // INCHI✔️❌:     int         multPrevEquStr;
    // INCHI✔️❌:     pINChI_Aux_Prev = NULL;
    // INCHI✔️❌:     pINChI_Aux_Taut = NULL;
    // INCHI✔️❌:     pINChI_Aux_Taut_Prev = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     mult = 0;
    // INCHI✔️❌:     bNext = 0;
    // INCHI✔️❌:     is = NULL;
    // INCHI✔️❌:     is2 = NULL;
    // INCHI✔️❌:     is0 = pINChISort;
    // INCHI✔️❌:     is20 = bSecondNonTautPass ? pINChISort2 : NULL;
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     eq2tautPrev = 1; /* pINChI_Aux_Prev (previous pINChI_Aux) does not exist */
    // INCHI✔️❌:     pPrevEquStr = NULL; /*, *pCurrEquStr;*/
    // INCHI✔️❌:     multPrevEquStr = 0;
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* For each connected component...    */
    // INCHI✔️❌:     for (i = 0; i <= num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         /* 1st (taut) pass: bOutType=OUT_TN  ; 2nd (non-taut pass) bOutType=OUT_NT */
    // INCHI✔️❌:         pINChI_Aux = ( i < num_components && ( is = is0 + i, 0 <= ( ii = GET_II( bOutType, is ) ) ) ) ? is->pINChI_Aux[ii] : NULL;
    // INCHI✔️❌:         if (bSecondNonTautPass)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* component that was output on the 1st pass */
    // INCHI✔️❌:             pINChI_Aux_Taut = ( i < num_components && ( is2 = is20 + i, 0 <= ( ii2 = GET_II( OUT_T1, is2 ) ) ) ) ? is2->pINChI_Aux[ii2] : NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /*================ compare non-iso non-taut equivalence info to non-iso taut ========*/
    // INCHI✔️❌:         eq2taut = bSecondNonTautPass && bOmitRepetitions &&
    // INCHI✔️❌:             Eql_INChI_Aux_Equ( pINChI_Aux, EQL_EQU, pINChI_Aux_Taut, EQL_EQU );
    // INCHI✔️❌:         eq2taut = eq2taut ? ( iiEQU | iitNONTAUT ) : 0;
    // INCHI✔️❌:         if (eq2taut)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* we may be here only in case of the second (non-taut) pass */
    // INCHI✔️❌:             /* current non-taut equivalence has been found to be same as tautomeric */
    // INCHI✔️❌:             if (pINChI_Aux_Prev && pINChI_Aux_Prev->nNumberOfAtoms)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* previous component exists */
    // INCHI✔️❌:                 if (bNext++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (bHasEquString( pINChI_Aux_Prev->nConstitEquNumbers, pINChI_Aux_Prev->nNumberOfAtoms ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* output previous component(s) equivalence since it was found to be non-trivial */
    // INCHI✔️❌:                     MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
    // INCHI✔️❌:                     MakeEquString( pINChI_Aux_Prev->nConstitEquNumbers,
    // INCHI✔️❌:                         pINChI_Aux_Prev->nNumberOfAtoms,
    // INCHI✔️❌:                         0, strbuf, TAUT_MODE, bOverflow );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     ; /* pINChI_Aux_Prev exists and has only trivial equivalence info */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:                 if (pINChI_Aux_Taut_Prev && pINChI_Aux_Taut_Prev->nNumberOfAtoms)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* previous non-taut component exists only in taut list */
    // INCHI✔️❌:                     if (bNext++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             /* we have found pINChI_Aux->nConstitEquNumbers same as in pINChI_Aux_Taut */
    // INCHI✔️❌:             /* output this (current) equivalence as '*', that is, same as tautomeric */
    // INCHI✔️❌:             /* that was printed on the 1st pass. */
    // INCHI✔️❌:             pCurrEquStr = EquString( eq2taut );
    // INCHI✔️❌:             if (multPrevEquStr && pPrevEquStr)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (pCurrEquStr && !strcmp( pCurrEquStr, pPrevEquStr ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     multPrevEquStr++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* new EqStr is different; output it */
    // INCHI✔️❌:                     if (bNext++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
    // INCHI✔️❌:                     pPrevEquStr = pCurrEquStr;
    // INCHI✔️❌:                     multPrevEquStr = 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pPrevEquStr = pCurrEquStr;
    // INCHI✔️❌:                 multPrevEquStr = 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             pINChI_Aux_Prev = NULL; /* pINChI_Aux_Prev does not exist since */
    // INCHI✔️❌:             pINChI_Aux_Taut_Prev = NULL; /* pINChI_Aux has just been printed */
    // INCHI✔️❌:             mult = 0;
    // INCHI✔️❌:             eq2tautPrev = 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (eq2tautPrev)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* at this point pINChI_Aux_Prev does not exist; however, pINChI_Aux */
    // INCHI✔️❌:                 /*might have been discovered and it is different from pINChI_Aux_Taut */
    // INCHI✔️❌:                 if (multPrevEquStr && pPrevEquStr)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* new EqStr is different; output it */
    // INCHI✔️❌:                     if (bNext++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
    // INCHI✔️❌:                     pPrevEquStr = NULL;
    // INCHI✔️❌:                     multPrevEquStr = 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 eq2tautPrev = 0;
    // INCHI✔️❌:                 pINChI_Aux_Prev = pINChI_Aux;
    // INCHI✔️❌:                 pINChI_Aux_Taut_Prev = pINChI_Aux_Taut;
    // INCHI✔️❌:                 mult = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* check whether pINChI_Aux and pINChI_Aux_Prev have identical non-trivial equivalence info */
    // INCHI✔️❌:                 eq2prev = bUseMulipliers &&
    // INCHI✔️❌:                     Eql_INChI_Aux_Equ( pINChI_Aux, EQL_EQU, pINChI_Aux_Prev, EQL_EQU );
    // INCHI✔️❌:                 if (eq2prev)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* eq. info is same and non-trivial */
    // INCHI✔️❌:                     mult++; /* mult = (number of non-empty equal items)-1 */
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* pINChI_Aux eq. info is either different or trivial. Output pINChI_Aux_Prev anyway */
    // INCHI✔️❌:                     if (bNext++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (pINChI_Aux_Prev && pINChI_Aux_Prev->nNumberOfAtoms)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (bHasEquString( pINChI_Aux_Prev->nConstitEquNumbers, pINChI_Aux_Prev->nNumberOfAtoms ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* pINChI_Aux_Prev exists and has equivalence info */
    // INCHI✔️❌:                             MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
    // INCHI✔️❌:                             MakeEquString( pINChI_Aux_Prev->nConstitEquNumbers, pINChI_Aux_Prev->nNumberOfAtoms, 0,
    // INCHI✔️❌:                                 strbuf, TAUT_MODE, bOverflow );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             ; /* pINChI_Aux_Prev exists and has only trivial equivalence info */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (bSecondNonTautPass && pINChI_Aux_Taut_Prev && pINChI_Aux_Taut_Prev->nNumberOfAtoms)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (bHasEquString( pINChI_Aux_Taut_Prev->nConstitEquNumbers, pINChI_Aux_Taut_Prev->nNumberOfAtoms ))
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 /* since pINChI_Aux_Prev does not exist, pINChI_Aux_Taut_Prev is non-tautomeric */
    // INCHI✔️❌:                                 /* and it has non-trivial equivalence info. This info has already been printed in the main section  */
    // INCHI✔️❌:                                 /*
    // INCHI✔️❌:                                 MakeDelim( sIdenticalValues, strbuf, bOverflow);
    // INCHI✔️❌:                                 */
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 ; /* pINChI_Aux_Taut_Prev exists and has only trivial equivalence info */
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌: #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             int stop = 1;   /* <BRKPT> */
    // INCHI✔️❌:                         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 pINChI_Aux_Prev = pINChI_Aux;
    // INCHI✔️❌:                 pINChI_Aux_Taut_Prev = pINChI_Aux_Taut;
    // INCHI✔️❌:                 mult = 0; /* we do not know whether the item is empty */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( strbuf->nUsedLength - nUsedLength0 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: str_AuxEqu

    let initial_length = string_buffer.nUsedLength;
    let delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let mut previous_aux: Option<crate::source_types::INChI_Aux> = None;
    let mut previous_taut_aux: Option<crate::source_types::INChI_Aux> = None;
    let mut previous_equivalence: Option<&'static str> = None;
    let mut previous_multiplier = 0_i32;
    let mut multiplier = 0_i32;
    let mut next = 0_i32;
    let mut equal_to_taut_previous = 1_i32;
    let aux_for = |heap: &SourceHeap,
                   sort_pointer: SourceConstPointer<INCHI_SORT>,
                   output: i32|
     -> Result<Option<crate::source_types::INChI_Aux>, SourceHeapError> {
        let sort = heap
            .slice(sort_pointer)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let Some(index) = get_ii(heap, output, &sort)? else {
            return Ok(None);
        };
        let pointer = sort.pINChI_Aux[index];
        if pointer.is_null() {
            return Ok(None);
        }
        Ok(Some(
            heap.slice(pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone(),
        ))
    };

    for component in 0..=number_of_components {
        let current_aux = if component < number_of_components && !sorted_inchi.is_null() {
            let sort_pointer = sorted_inchi.as_const().offset(i64::from(component))?;
            aux_for(&*heap, sort_pointer, output_type)?
        } else {
            None
        };
        let current_taut_aux =
            if second_non_taut_pass != 0 && component < number_of_components && !sorted_inchi2.is_null() {
                let sort_pointer = sorted_inchi2.as_const().offset(i64::from(component))?;
                aux_for(&*heap, sort_pointer, OUT_T1 as i32)?
            } else {
                None
            };

        let equal_to_taut = second_non_taut_pass != 0
            && omit_repetitions != 0
            && Eql_INChI_Aux_Equ(
                &*heap,
                current_aux.as_ref(),
                crate::source_types::EQL_EQU as i32,
                current_taut_aux.as_ref(),
                crate::source_types::EQL_EQU as i32,
            )? != 0;
        if equal_to_taut {
            if previous_aux.as_ref().is_some_and(|aux| aux.nNumberOfAtoms != 0) {
                if next != 0 {
                    MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                let aux = previous_aux.as_ref().unwrap();
                if bHasEquString(&*heap, aux.nConstitEquNumbers.as_const(), aux.nNumberOfAtoms)? != 0 {
                    let star = heap.allocate_model_storage(vec![b'*' as i8, 0])?;
                    MakeMult(
                        heap,
                        multiplier.wrapping_add(1),
                        star.as_const(),
                        string_buffer,
                        0,
                        overflow,
                    )?;
                    MakeEquString(
                        heap,
                        aux.nConstitEquNumbers.as_const(),
                        aux.nNumberOfAtoms,
                        0,
                        string_buffer,
                        taut_mode,
                        overflow,
                    )?;
                }
            } else if previous_taut_aux.as_ref().is_some_and(|aux| aux.nNumberOfAtoms != 0) {
                if next != 0 {
                    MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
            }
            let current_equivalence = crate::source::base::ichiprt1::EquString(
                (crate::source_types::iiEQU | crate::source_types::iitNONTAUT) as i32,
            );
            if previous_multiplier != 0 {
                if previous_equivalence == Some(current_equivalence) {
                    previous_multiplier = previous_multiplier.wrapping_add(1);
                } else {
                    if next != 0 {
                        MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let value = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let value =
                        heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain([0]).collect())?;
                    MakeEqStr(heap, value.as_const(), previous_multiplier, string_buffer, overflow)?;
                    previous_equivalence = Some(current_equivalence);
                    previous_multiplier = 1;
                }
            } else {
                previous_equivalence = Some(current_equivalence);
                previous_multiplier = 1;
            }
            previous_aux = None;
            previous_taut_aux = None;
            multiplier = 0;
            equal_to_taut_previous = 1;
            continue;
        }

        if equal_to_taut_previous != 0 {
            if previous_multiplier != 0 {
                if next != 0 {
                    MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                let value = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
                let value = heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain([0]).collect())?;
                MakeEqStr(heap, value.as_const(), previous_multiplier, string_buffer, overflow)?;
                previous_equivalence = None;
                previous_multiplier = 0;
            }
            equal_to_taut_previous = 0;
            previous_aux = current_aux;
            previous_taut_aux = current_taut_aux;
            multiplier = 0;
            continue;
        }

        let equal_to_previous = use_multipliers != 0
            && Eql_INChI_Aux_Equ(
                &*heap,
                current_aux.as_ref(),
                crate::source_types::EQL_EQU as i32,
                previous_aux.as_ref(),
                crate::source_types::EQL_EQU as i32,
            )? != 0;
        if equal_to_previous {
            multiplier = multiplier.wrapping_add(1);
            continue;
        }
        if next != 0 {
            MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
        }
        next = next.wrapping_add(1);
        if let Some(aux) = previous_aux.as_ref() {
            if aux.nNumberOfAtoms != 0
                && bHasEquString(&*heap, aux.nConstitEquNumbers.as_const(), aux.nNumberOfAtoms)? != 0
            {
                let star = heap.allocate_model_storage(vec![b'*' as i8, 0])?;
                MakeMult(
                    heap,
                    multiplier.wrapping_add(1),
                    star.as_const(),
                    string_buffer,
                    0,
                    overflow,
                )?;
                MakeEquString(
                    heap,
                    aux.nConstitEquNumbers.as_const(),
                    aux.nNumberOfAtoms,
                    0,
                    string_buffer,
                    taut_mode,
                    overflow,
                )?;
            }
        }
        previous_aux = current_aux;
        previous_taut_aux = current_taut_aux;
        multiplier = 0;
    }
    if previous_multiplier != 0 {
        if next != 0 {
            MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
        }
        let value = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
        let value = heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain([0]).collect())?;
        MakeEqStr(heap, value.as_const(), previous_multiplier, string_buffer, overflow)?;
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_AuxIsoEqu(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    sorted_inchi2: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    taut_mode: i32,
    number_of_components: i32,
    second_non_taut_pass: i32,
    omit_repetitions: i32,
    use_multipliers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:2962 str_AuxIsoEqu
    // INCHI✔️❌: complete configured source frame follows verbatim; the non-release debug block is inactive.
    /*
    int str_AuxIsoEqu( INCHI_SORT       *pINChISort,
                       INCHI_SORT       *pINChISort2,
                       INCHI_IOS_STRING *strbuf,
                       int              *bOverflow,
                       int              bOutType,
                       int              TAUT_MODE,
                       int              num_components,
                       int              bSecondNonTautPass,
                       int              bOmitRepetitions,
                       int              bUseMulipliers )
    {
        int          i, ii, ii2, nUsedLength0;
        INCHI_SORT   *is, *is2, *is0, *is20;
        INChI_Aux    *pINChI_Aux, *pINChI_Aux_Prev, *pINChI_Aux_Taut, *pINChI_Aux_Taut_Prev;
        int          mult, eq2prev, eq2taut, eq2tautPrev, bNext;
        const char  *pPrevEquStr, *pCurrEquStr;
        int         multPrevEquStr;
        pINChI_Aux = NULL;
        pINChI_Aux_Prev = NULL;
        pINChI_Aux_Taut = NULL;
        pINChI_Aux_Taut_Prev = NULL;
        mult = 0;
        bNext = 0;
        is = NULL;
        is2 = NULL;
        is0 = pINChISort;
        is20 = bSecondNonTautPass ? pINChISort2 : NULL;
        /* djb-rwth: removing redundant code */
        eq2tautPrev = 1; /* pINChI_Aux_Prev (previous pINChI_Aux) does not exist */
        pPrevEquStr = NULL; /*, *pCurrEquStr;*/
        multPrevEquStr = 0;
        nUsedLength0 = strbuf->nUsedLength;

        /* For each connected component...    */
        for (i = 0; i <= num_components; i++)
        {

            /* 1st (taut) pass: bOutType=OUT_TN  ; 2nd (non-taut pass) bOutType=OUT_NT */
            pINChI_Aux = ( i < num_components && ( is = is0 + i, 0 <= ( ii = GET_II( bOutType, is ) ) ) ) ? is->pINChI_Aux[ii] : NULL;
            if (bSecondNonTautPass)
            {
                /* component that was output on the 1st pass */
                pINChI_Aux_Taut = ( i < num_components && ( is2 = is20 + i, 0 <= ( ii2 = GET_II( OUT_T1, is2 ) ) ) ) ? is2->pINChI_Aux[ii2] : NULL;
            }
            /*================ compare iso non-taut equivalence info to non-iso taut ========*/
            eq2taut = 0;
            if (bSecondNonTautPass && bOmitRepetitions && pINChI_Aux && pINChI_Aux->bIsIsotopic)
            {
                /**************************************************
                * compare isotopic non-tautomeric equivalence to:
                *    a) tautomeric
                *    b) non-tautomeric
                *    c) isotopic tautomeric
                */
                if (!eq2taut)
                {
                    /* compare isotopic non-tautomeric equivalence to tautomeric */
                    eq2taut = Eql_INChI_Aux_Equ( pINChI_Aux, EQL_EQU_ISO, pINChI_Aux_Taut, EQL_EQU );
                    /* equ   non-taut     isotopic = tautomeric*/
                    eq2taut = eq2taut ? ( iiEQU | iitNONTAUT | iitISO ) : 0;
                }
                if (!eq2taut)
                {
                    /* compare isotopic non-tautomeric equivalence to non-tautomeric */
                    eq2taut = Eql_INChI_Aux_Equ( pINChI_Aux, EQL_EQU_ISO, pINChI_Aux, EQL_EQU );
                    /* equ  non-taut    isotopic = non-tautomeric*/
                    eq2taut = eq2taut ? ( iiEQU | iitNONTAUT | iitISO | iiEq2NONTAUT ) : 0;
                }
                if (!eq2taut)
                {
                    /* compare isotopic non-tautomeric equivalence to isotopic tautomeric */
                    eq2taut = Eql_INChI_Aux_Equ( pINChI_Aux, EQL_EQU_ISO, pINChI_Aux_Taut, EQL_EQU_ISO );
                    /* equ   non-taut     isotopic = isotopic tautomeric*/
                    eq2taut = eq2taut ? ( iiEQU | iitNONTAUT | iitISO | iiEq2ISO ) : 0;
                }
            }
            else
            {
                if (!bSecondNonTautPass && bOmitRepetitions && pINChI_Aux && pINChI_Aux->bIsIsotopic)
                {
                    /**************************************************
                    * compare isotopic tautomeric equivalence to:
                    *    a) non-isotopic tautomeric
                    */
                    if (!eq2taut)
                    {
                        /* compare isotopic tautomeric equivalence to tautomeric */
                        eq2taut = Eql_INChI_Aux_Equ( pINChI_Aux, EQL_EQU_ISO, pINChI_Aux, EQL_EQU );
                        /* equ   taut-isotopic = tautomeric*/
                        eq2taut = eq2taut ? ( iiEQU | iitISO ) : 0;
                    }
                }
            }
            if (eq2taut)
            {
                /* we may be here only in case of the second (non-taut) pass */
                /* current non-taut equivalence has been found to be same as tautomeric */
                if (pINChI_Aux_Prev && pINChI_Aux_Prev->nNumberOfAtoms)
                {
                    /* previous component exists */
                    if (bNext++)
                    {
                        MakeDelim( sCompDelim, strbuf, bOverflow );
                    }
                    if (bHasEquString( pINChI_Aux_Prev->nConstitEquIsotopicNumbers, pINChI_Aux_Prev->nNumberOfAtoms ))
                    {
                        /* output previous component(s) equivalence since it was found to be non-trivial */
                        MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
                        MakeEquString( pINChI_Aux_Prev->nConstitEquIsotopicNumbers, pINChI_Aux_Prev->nNumberOfAtoms, 0,
                            strbuf, TAUT_MODE, bOverflow );
                    }
                    else
                    {
                        ; /* pINChI_Aux_Prev exists and has only trivial equivalence info */
                    }
                }
                else
                {
                    if (pINChI_Aux_Taut_Prev && pINChI_Aux_Taut_Prev->nNumberOfAtoms)
                    {
                        /* previous non-taut component exists only in taut list */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                    }
                }
                /* we have found pINChI_Aux->pINChI_Aux->nConstitEquIsotopicNumbers same as in pINChI_Aux_Taut */
                /* output this (current) equivalence as '*', that is, same as tautomeric */
                /* that was printed on the 1st pass. */
                pCurrEquStr = EquString( eq2taut );
                if (multPrevEquStr && pPrevEquStr)
                {
                    if (pCurrEquStr && !strcmp( pCurrEquStr, pPrevEquStr ))
                    {
                        multPrevEquStr++;
                    }
                    else
                    {
                        /* new EqStr is different; output it */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                        pPrevEquStr = pCurrEquStr;
                        multPrevEquStr = 1;
                    }
                }
                else
                {
                    pPrevEquStr = pCurrEquStr;
                    multPrevEquStr = 1;
                }
                pINChI_Aux_Prev = NULL; /* pINChI_Aux_Prev does not exist since */
                pINChI_Aux_Taut_Prev = NULL; /* pINChI_Aux has just been printed */
                mult = 0;
                eq2tautPrev = 1;
            }
            else
            {
                if (eq2tautPrev)
                {
                    /* at this point pINChI_Aux_Prev does not exist; however, pINChI_Aux */
                    /*might have been discovered and it is different from pINChI_Aux_Taut */
                    if (multPrevEquStr && pPrevEquStr)
                    {
                        /* new EqStr is different; output it */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                        pPrevEquStr = NULL;
                        multPrevEquStr = 0;
                    }
                    eq2tautPrev = 0;
                    pINChI_Aux_Prev = pINChI_Aux;
                    pINChI_Aux_Taut_Prev = pINChI_Aux_Taut;
                    mult = 0;
                }
                else
                {
                    /* check whether pINChI_Aux and pINChI_Aux_Prev have identical non-trivial equivalence info */
                    eq2prev = bUseMulipliers && Eql_INChI_Aux_Equ( pINChI_Aux, EQL_EQU_ISO, pINChI_Aux_Prev, EQL_EQU_ISO );
                    if (eq2prev)
                    {
                        /* eq. info is same and non-trivial */
                        mult++; /* mult = (number of non-empty equal items)-1 */
                        continue;
                    }
                    else
                    {
                        /* pINChI_Aux eq. info is either different or trivial. Output pINChI_Aux_Prev anyway */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        if (pINChI_Aux_Prev && pINChI_Aux_Prev->nNumberOfAtoms)
                        {
                            if (bHasEquString( pINChI_Aux_Prev->nConstitEquIsotopicNumbers, pINChI_Aux_Prev->nNumberOfAtoms ))
                            {
                                /* pINChI_Aux_Prev exists and has equivalence info */
                                MakeMult( mult + 1, "*", strbuf, 0, bOverflow );
                                MakeEquString( pINChI_Aux_Prev->nConstitEquIsotopicNumbers, pINChI_Aux_Prev->nNumberOfAtoms, 0,
                                    strbuf, TAUT_MODE, bOverflow );
                            }
                            else
                            {
                                ; /* pINChI_Aux_Prev exists and has only trivial equivalence info */
                            }
                        }
                        else
                        {
                            if (bSecondNonTautPass && pINChI_Aux_Taut_Prev && pINChI_Aux_Taut_Prev->nNumberOfAtoms)
                            {
                                if (bHasEquString( pINChI_Aux_Taut_Prev->nConstitEquIsotopicNumbers, pINChI_Aux_Taut_Prev->nNumberOfAtoms ))
                                {
                                    /* since pINChI_Aux_Prev does not exist, pINChI_Aux_Taut_Prev is non-tautomeric */
                                    /* and it has non-trivial equivalence info. This info has already been printed in the main section  */
                                    /*
                                    MakeDelim( sIdenticalValues, strbuf, bOverflow);
                                    */
                                }
                                else
                                {
                                    ; /* pINChI_Aux_Taut_Prev exists and has only trivial equivalence info */
                                }
                            }
    #if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
                            else
                            {
                                int stop = 1;   /* <BRKPT> */
                            }
    #endif
                        }
                    }
                    pINChI_Aux_Prev = pINChI_Aux;
                    pINChI_Aux_Taut_Prev = pINChI_Aux_Taut;
                    mult = 0; /* we do not know whether the item is empty */
                }
            }
        }

        return ( strbuf->nUsedLength - nUsedLength0 );
    }
    */
    // END INCHI C FUNCTION: str_AuxIsoEqu

    let initial_length = string_buffer.nUsedLength;
    let delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let star = heap.allocate_model_storage(vec![b'*' as i8, 0])?;
    let mut previous_aux: Option<crate::source_types::INChI_Aux> = None;
    let mut previous_taut_aux: Option<crate::source_types::INChI_Aux> = None;
    let mut previous_equivalence: Option<&'static str> = None;
    let mut previous_equivalence_multiplier = 0_i32;
    let mut multiplier = 0_i32;
    let mut next = 0_i32;
    let mut equal_to_taut_previous = true;
    let aux_for = |heap: &SourceHeap,
                   sort_pointer: SourceConstPointer<INCHI_SORT>,
                   output: i32|
     -> Result<Option<crate::source_types::INChI_Aux>, SourceHeapError> {
        let sort = heap
            .slice(sort_pointer)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let Some(slot) = get_ii(heap, output, sort)? else {
            return Ok(None);
        };
        let pointer = sort.pINChI_Aux[slot];
        if pointer.is_null() {
            return Ok(None);
        }
        Ok(Some(
            heap.slice(pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone(),
        ))
    };

    for component in 0..=number_of_components {
        let current_aux = if component < number_of_components && !sorted_inchi.is_null() {
            aux_for(heap, sorted_inchi.as_const().offset(i64::from(component))?, output_type)?
        } else {
            None
        };
        let current_taut_aux =
            if second_non_taut_pass != 0 && component < number_of_components && !sorted_inchi2.is_null() {
                aux_for(
                    heap,
                    sorted_inchi2.as_const().offset(i64::from(component))?,
                    OUT_T1 as i32,
                )?
            } else {
                None
            };
        let equal_equivalence = |left: Option<&crate::source_types::INChI_Aux>,
                                 left_mode: i32,
                                 right: Option<&crate::source_types::INChI_Aux>,
                                 right_mode: i32|
         -> Result<bool, SourceHeapError> {
            Ok(Eql_INChI_Aux_Equ(heap, left, left_mode, right, right_mode)? != 0)
        };
        let mut equivalence = 0_i32;
        if second_non_taut_pass != 0
            && omit_repetitions != 0
            && current_aux.as_ref().is_some_and(|aux| aux.bIsIsotopic != 0)
        {
            if equal_equivalence(
                current_aux.as_ref(),
                crate::source_types::EQL_EQU_ISO as i32,
                current_taut_aux.as_ref(),
                crate::source_types::EQL_EQU as i32,
            )? {
                equivalence =
                    (crate::source_types::iiEQU | crate::source_types::iitNONTAUT | crate::source_types::iitISO) as i32;
            } else if equal_equivalence(
                current_aux.as_ref(),
                crate::source_types::EQL_EQU_ISO as i32,
                current_aux.as_ref(),
                crate::source_types::EQL_EQU as i32,
            )? {
                equivalence = (crate::source_types::iiEQU
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iitISO
                    | crate::source_types::iiEq2NONTAUT) as i32;
            } else if equal_equivalence(
                current_aux.as_ref(),
                crate::source_types::EQL_EQU_ISO as i32,
                current_taut_aux.as_ref(),
                crate::source_types::EQL_EQU_ISO as i32,
            )? {
                equivalence = (crate::source_types::iiEQU
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iitISO
                    | crate::source_types::iiEq2ISO) as i32;
            }
        } else if second_non_taut_pass == 0
            && omit_repetitions != 0
            && current_aux.as_ref().is_some_and(|aux| aux.bIsIsotopic != 0)
            && equal_equivalence(
                current_aux.as_ref(),
                crate::source_types::EQL_EQU_ISO as i32,
                current_aux.as_ref(),
                crate::source_types::EQL_EQU as i32,
            )?
        {
            equivalence = (crate::source_types::iiEQU | crate::source_types::iitISO) as i32;
        }

        if equivalence != 0 {
            if let Some(aux) = previous_aux.as_ref()
                && aux.nNumberOfAtoms != 0
            {
                if next != 0 {
                    MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                if bHasEquString(heap, aux.nConstitEquIsotopicNumbers.as_const(), aux.nNumberOfAtoms)? != 0 {
                    MakeMult(
                        heap,
                        multiplier.wrapping_add(1),
                        star.as_const(),
                        string_buffer,
                        0,
                        overflow,
                    )?;
                    MakeEquString(
                        heap,
                        aux.nConstitEquIsotopicNumbers.as_const(),
                        aux.nNumberOfAtoms,
                        0,
                        string_buffer,
                        taut_mode,
                        overflow,
                    )?;
                }
            } else if previous_aux.is_none() && previous_taut_aux.as_ref().is_some_and(|aux| aux.nNumberOfAtoms != 0) {
                if next != 0 {
                    MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
            }

            let current_equivalence = crate::source::base::ichiprt1::EquString(equivalence);
            if previous_equivalence_multiplier != 0 {
                if previous_equivalence == Some(current_equivalence) {
                    previous_equivalence_multiplier = previous_equivalence_multiplier.wrapping_add(1);
                } else {
                    if next != 0 {
                        MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let value = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let value =
                        heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain([0]).collect())?;
                    MakeEqStr(
                        heap,
                        value.as_const(),
                        previous_equivalence_multiplier,
                        string_buffer,
                        overflow,
                    )?;
                    previous_equivalence = Some(current_equivalence);
                    previous_equivalence_multiplier = 1;
                }
            } else {
                previous_equivalence = Some(current_equivalence);
                previous_equivalence_multiplier = 1;
            }
            previous_aux = None;
            previous_taut_aux = None;
            multiplier = 0;
            equal_to_taut_previous = true;
            continue;
        }

        if equal_to_taut_previous {
            if previous_equivalence_multiplier != 0 {
                if next != 0 {
                    MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                let value = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
                let value = heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain([0]).collect())?;
                MakeEqStr(
                    heap,
                    value.as_const(),
                    previous_equivalence_multiplier,
                    string_buffer,
                    overflow,
                )?;
                previous_equivalence = None;
                previous_equivalence_multiplier = 0;
            }
            equal_to_taut_previous = false;
            previous_aux = current_aux;
            previous_taut_aux = current_taut_aux;
            multiplier = 0;
            continue;
        }

        let equal_to_previous = use_multipliers != 0
            && equal_equivalence(
                current_aux.as_ref(),
                crate::source_types::EQL_EQU_ISO as i32,
                previous_aux.as_ref(),
                crate::source_types::EQL_EQU_ISO as i32,
            )?;
        if equal_to_previous {
            multiplier = multiplier.wrapping_add(1);
            continue;
        }
        if next != 0 {
            MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
        }
        next = next.wrapping_add(1);
        if let Some(aux) = previous_aux.as_ref()
            && aux.nNumberOfAtoms != 0
            && bHasEquString(heap, aux.nConstitEquIsotopicNumbers.as_const(), aux.nNumberOfAtoms)? != 0
        {
            MakeMult(
                heap,
                multiplier.wrapping_add(1),
                star.as_const(),
                string_buffer,
                0,
                overflow,
            )?;
            MakeEquString(
                heap,
                aux.nConstitEquIsotopicNumbers.as_const(),
                aux.nNumberOfAtoms,
                0,
                string_buffer,
                taut_mode,
                overflow,
            )?;
        }
        previous_aux = current_aux;
        previous_taut_aux = current_taut_aux;
        multiplier = 0;
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_AuxInvIsoSp3(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    sorted_inchi2: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    taut_mode: i32,
    number_of_components: i32,
    second_non_taut_pass: i32,
    omit_repetitions: i32,
    use_multipliers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:3213 str_AuxInvIsoSp3
    // INCHI✔️❌: complete configured source frame follows verbatim; inactive FIX_EMPTY_LAYER_BUG and non-release debug blocks are excluded.
    /*
    int str_AuxInvIsoSp3( INCHI_SORT       *pINChISort,
                          INCHI_SORT       *pINChISort2,
                          INCHI_IOS_STRING *strbuf,
                          int              *bOverflow,
                          int              bOutType,
                          int              TAUT_MODE,
                          int              num_components,
                          int              bSecondNonTautPass,
                          int              bOmitRepetitions,
                          int              bUseMulipliers )
    {
        int          i, ii, ii2, nUsedLength0;
        INCHI_SORT   *is, *is2, *is0, *is20;
        INChI        *pINChI, *pINChI_Prev, *pINChI_Taut, *pINChI_Taut_Prev;
        INChI_Stereo *Stereo, *Stereo_Prev, *Stereo_Taut, *Stereo_Taut_Prev;
        int          mult, eq2prev, eq2taut, eq2tautPrev, bNext;
        const char  *pPrevEquStr, *pCurrEquStr;
        int         multPrevEquStr;
        /********************************
        inverted isotopic sp3
        *********************************/
        pINChI_Taut = NULL;
        pINChI_Prev = NULL;
        pINChI_Taut_Prev = NULL;
        mult = 0;
        bNext = 0;
        is = NULL;
        is2 = NULL;
        is0 = pINChISort;
        is20 = bSecondNonTautPass ? pINChISort2 : NULL;
        /* djb-rwth: removing redundant code */
        eq2tautPrev = 1; /* pINChI_Prev (previous pINChI) does not exist */
        pPrevEquStr = NULL; /*, *pCurrEquStr;*/
        multPrevEquStr = 0;
        nUsedLength0 = strbuf->nUsedLength;

        /* For each connected component...    */
        for (i = 0; i <= num_components; i++)
        {

            /* 1st (taut) pass: bOutType=OUT_TN  ; 2nd (non-taut pass) bOutType=OUT_NT */
            pINChI = ( i < num_components && ( is = is0 + i, 0 <= ( ii = GET_II( bOutType, is ) ) ) ) ? is->pINChI[ii] : NULL;
            /*================ compare sp2 to previous =====================*/
            if (bSecondNonTautPass)
            {
                /* component that was output on the 1st pass */
                pINChI_Taut = ( i < num_components && ( is2 = is20 + i, 0 <= ( ii2 = GET_II( OUT_T1, is2 ) ) ) ) ? is2->pINChI[ii2] : NULL;
            }
            eq2taut = 0;
            /*========= if bSecondNonTautPass then compare iso non-taut stereo to other stereo ========*/
            if (bSecondNonTautPass && bOmitRepetitions && pINChI && pINChI->nNumberOfIsotopicAtoms + pINChI->nNumberOfIsotopicTGroups > 0)
            {
                /* compare non-tautomeric isotopic inverted to:
                *   a) tautomeric inverted
                *   b) *non-tautomeric inverted
                *   c) *isotopic tautomeric inverted
                *   d) Inverted(tautomeric)
                *   e) *Inverted(tautomeric isotopic)
                *   f) Inverted(non-tautomeric)
                *   g) *Inverted(non-tautomeric isotopic)
                */
                /* a) compare non-tautomeric isotopic inverted to tautomeric inverted */
                if (!eq2taut)
                {
                    eq2taut = pINChI && pINChI_Taut &&
                        /* non-taut inverted */          /* taut invertedc */
                        ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Taut = pINChI_Taut->Stereo ) &&
                        Eql_INChI_Stereo( Stereo, EQL_SP3_INV, Stereo_Taut, EQL_SP3_INV, 0 );
                    /* stereo-inv    isotopic  non-taut =  taut (stereo-inv) */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iitNONTAUT ) : 0;
                }
                /* b) compare non-tautomeric isotopic inverted to non-tautomeric inverted */
                if (!eq2taut)
                {
                    eq2taut = pINChI &&                    /* it is non-taut non-iso stereo */
                        ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Taut = pINChI->Stereo ) &&
                        Eql_INChI_Stereo( Stereo, EQL_SP3_INV, Stereo_Taut, EQL_SP3_INV, 0 );
                    /* stereo-inv    isotopic non-taut =  non-taut stereo-inv */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iitNONTAUT | iiEq2NONTAUT ) : 0;
                }
                /* c) compare non-tautomeric isotopic inverted to isotopic tautomeric inverted */
                if (!eq2taut)
                {
                    eq2taut = pINChI && pINChI_Taut &&
                        /* non-taut iso. inverted */             /* taut iso. inverted */
                        ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Taut = pINChI_Taut->StereoIsotopic ) &&
                        Eql_INChI_Stereo( Stereo, EQL_SP3_INV, Stereo_Taut, EQL_SP3_INV, 0 );
                    /* stereo-inv    isotopic  non-taut =  taut iso. stereo-inv */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iitNONTAUT | iiEq2ISO ) : 0;
                }
                /* d) compare non-tautomeric inverted to Inverted(tautomeric stereo) */
                if (!eq2taut)
                {
                    eq2taut = pINChI && pINChI_Taut &&
                        ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Taut = pINChI_Taut->Stereo ) &&
                        Eql_INChI_Stereo( Stereo, EQL_SP3_INV, Stereo_Taut, EQL_SP3_INV, 0 );
                    /* stereo-inv   isotopic  non-taut =  Inv(non-iso taut stereo) */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iitNONTAUT | iiEq2INV ) : 0;
                }
                /* e) compare non-tautomeric inverted to Inverted(isotopic tautomeric stereo) */
                if (!eq2taut)
                {
                    eq2taut = pINChI && pINChI_Taut &&
                        ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Taut = pINChI_Taut->StereoIsotopic ) &&
                        Eql_INChI_Stereo( Stereo, EQL_SP3_INV, Stereo_Taut, EQL_SP3, 0 );
                    /* stereo-inv   isotopic  non-taut =  Inv(iso taut stereo) */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iitNONTAUT | iiEq2INV | iiEq2ISO ) : 0;
                }
                /* f) compare non-tautomeric isotopic inverted to Inverted(non-tautomeric stereo) */
                if (!eq2taut)
                {
                    eq2taut = pINChI &&                    /* it is non-taut non-iso stereo */
                        ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Taut = pINChI->Stereo ) &&
                        Eql_INChI_Stereo( Stereo, EQL_SP3_INV, Stereo_Taut, EQL_SP3, 0 );
                    /* stereo-inv   isotopic    non-taut =  Inv(non-taut stereo) */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iitNONTAUT | iiEq2INV | iiEq2NONTAUT ) : 0;
                }
                /* g) compare non-tautomeric isotopic inverted to Inverted(non-tautomeric isotopic stereo) */
                if (!eq2taut)
                {
                    eq2taut = pINChI &&                    /* it is non-taut non-iso stereo */
                        ( Stereo = pINChI->StereoIsotopic ) &&
                        Eql_INChI_Stereo( Stereo, EQL_SP3_INV, Stereo, EQL_SP3, 0 );
                    /* stereo-inv    isotopic  non-taut =   Inv( iso non-taut stereo) */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iitNONTAUT | iiEq2INV | iiEq2ISO | iiEq2NONTAUT ) : 0;
                }
            }
            else
            {
                /*========= if not bSecondNonTautPass then compare inv taut stereo to various stereo ========*/
                if (!bSecondNonTautPass && bOmitRepetitions && pINChI &&
                    ( pINChI->nNumberOfIsotopicAtoms > 0 ||
                      pINChI->nNumberOfIsotopicTGroups > 0 ||
                      (pINChI->nPossibleLocationsOfIsotopicH && pINChI->nPossibleLocationsOfIsotopicH[0] > 1) )) /* djb-rwth: addressing LLVM warning */
                {
                    /* compare tautomeric isotopic stereo-inverted to:
                    *    a) tautomeric stereo-inverted
                    *    b) Inverted(tautomeric stereo)
                    *    c) Inverted(tautomeric isotopic stereo)
                    */
                    /* a) compare tautomeric isotopic stereo-inverted to tautomeric stereo-inverted */
                    if (!eq2taut)
                    {
                        eq2taut = pINChI &&
                            ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Taut = pINChI->Stereo ) &&
                            Eql_INChI_Stereo( Stereo, EQL_SP3_INV, Stereo_Taut, EQL_SP3_INV, 0 );
                        /* stereo-inv  isotopic taut =  taut stereo-inv */
                        eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO ) : 0;
                    }
                    /* b) compare tautomeric isotopic stereo-inverted to Inverted(tautomeric stereo) */
                    if (!eq2taut)
                    {
                        eq2taut = pINChI &&
                            ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Taut = pINChI->Stereo ) &&
                            Eql_INChI_Stereo( Stereo, EQL_SP3_INV, Stereo_Taut, EQL_SP3, 0 );
                        /* stereo-inv   isotopic taut =  Inv(taut stereo) */
                        eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iiEq2INV ) : 0;
                    }
                    /* c) compare tautomeric isotopic stereo-inverted to Inverted(tautomeric isotopic stereo) */
                    if (!eq2taut)
                    {
                        eq2taut = pINChI &&
                            ( Stereo = pINChI->StereoIsotopic ) &&
                            Eql_INChI_Stereo( Stereo, EQL_SP3_INV, Stereo, EQL_SP3, 0 );
                        /* stereo-inv   isotopic taut =  Inv(taut iso stereo) */
                        eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iiEq2INV | iiEq2ISO ) : 0;
                    }
                }
            }
            if (eq2taut)
            {
                /* we may be here only in case of the current layer found equal in another layer the same component */
                if (pINChI_Prev && pINChI_Prev->nNumberOfAtoms)
                {
                    /* previous component exists; output it before output the current component */
                    if (bNext++)
                    {
                        MakeDelim( sCompDelim, strbuf, bOverflow );
                    }
                    if (( Stereo_Prev = pINChI_Prev->StereoIsotopic ) && Stereo_Prev->nNumberOfStereoCenters > 0)
                    {
                        MakeMult( mult + 1, "*", strbuf, 0, bOverflow );

                        MakeStereoString( Stereo_Prev->nNumber, NULL, Stereo_Prev->t_parityInv,
                            0, Stereo_Prev->nNumberOfStereoCenters,
                            strbuf, TAUT_MODE, bOverflow );
                    }
                }
                else
                {
                    if (pINChI_Taut_Prev && pINChI_Taut_Prev->nNumberOfAtoms)
                    {
                        /* previous non-taut component exists only in taut list */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        /* do not output stereo of non-tautomeric in non-taut layer: it has been output in the main layer */
                    }
                }
                /* we have found another (previously printed) layer of the current component equal to this layer */
                /* output this (current) equivalence mark = EquString(eq2taut) */
                pCurrEquStr = EquString( eq2taut );
                if (multPrevEquStr && pPrevEquStr)
                {
                    if (pCurrEquStr && !strcmp( pCurrEquStr, pPrevEquStr ))
                    {
                        multPrevEquStr++;
                    }
                    else
                    {
                        /* new EqStr is different; output it */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                        pPrevEquStr = pCurrEquStr;
                        multPrevEquStr = 1;
                    }
                }
                else
                {
                    pPrevEquStr = pCurrEquStr;
                    multPrevEquStr = 1;
                }
                pINChI_Prev = NULL; /* pINChI_Prev sp2 does not exist since */
                pINChI_Taut_Prev = NULL; /* pINChI has just been printed */
                mult = 0;
                eq2tautPrev = 1;     /* pINChI_Prev and pINChI_Taut_Prev have already been output */
            }
            else
            {
                if (eq2tautPrev)
                {
                    /* at this point pINChI_Prev does not exist; however, pINChI */
                    /*might have been discovered and it is different from pINChI_Taut */
                    if (multPrevEquStr && pPrevEquStr)
                    {
                        /* new EqStr is different; output it */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                        pPrevEquStr = NULL;
                        multPrevEquStr = 0;
                    }
                    eq2tautPrev = 0;
                    pINChI_Prev = pINChI;
                    pINChI_Taut_Prev = pINChI_Taut;
                    mult = 0;
                }
                else
                {
                    /* current layer is different from previously printed layers of the current component */
                    /* compare the current layer to this layer of the previous component: */
                    /* check whether pINChI and pINChI_Prev have non-zero identical stereo sp2 */
                    /*================ compare iso sp3 to previous =====================*/
                    eq2prev = bUseMulipliers &&
                        pINChI && pINChI->nNumberOfIsotopicAtoms + pINChI->nNumberOfIsotopicTGroups > 0 &&
                        pINChI_Prev && pINChI_Prev->nNumberOfIsotopicAtoms + pINChI_Prev->nNumberOfIsotopicTGroups > 0 &&
                        /* do both have stereo? */
                        ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Prev = pINChI_Prev->StereoIsotopic ) &&
                        /* is their inverted stereo same? */
                        Eql_INChI_Stereo( Stereo, EQL_SP3_INV, Stereo_Prev, EQL_SP3_INV, 0 );
                    if (eq2prev)
                    {
                        mult++; /* mult = (number of non-empty equal items)-1 */
                        continue;
                    }
                    else
                    {
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        if (pINChI_Prev && pINChI_Prev->nNumberOfAtoms && pINChI_Prev->nNumberOfIsotopicAtoms + pINChI_Prev->nNumberOfIsotopicTGroups > 0)
                        {
                            if (( Stereo_Prev = pINChI_Prev->StereoIsotopic ) &&
                                Stereo_Prev->nNumberOfStereoCenters > 0 && Stereo_Prev->nCompInv2Abs)
                            {
                                MakeMult( mult + 1, "*", strbuf, 0, bOverflow );

                                MakeStereoString( Stereo_Prev->nNumberInv, NULL, Stereo_Prev->t_parityInv,
                                    0, Stereo_Prev->nNumberOfStereoCenters,
                                    strbuf, TAUT_MODE, bOverflow );
                            }
                            /* else sp3 info is not present in pINChI_Prev */
                        }
                        else
                        {
                            /* do not print pINChI_Prev because it either do not exist of have already been printed */
                            if (bSecondNonTautPass && pINChI_Taut_Prev && pINChI_Taut_Prev->nNumberOfAtoms)
                            {
                                if (( Stereo_Taut_Prev = pINChI_Taut_Prev->StereoIsotopic ) &&
                                    Stereo_Taut_Prev->nNumberOfStereoCenters > 0 && Stereo_Taut_Prev->nCompInv2Abs)
                                {
                                    /* since pINChI_Prev does not exist, pINChI_Taut_Prev is non-tautomeric */
                                    /* and it has non-trivial inv sp3 info. It has already been printed in the main section */
                                    /*
                                    tot_len += MakeDelim( sIdenticalValues, strbuf, bOverflow);
                                    */
                                    ;/* pINChI_Taut_Prev sp3 info was output in the main stereo section */
                                }
                                else
                                {
                                    ; /* pINChI_Taut_Prev exists and has not sp3 info */
                                }
                            }
                        }
                    }
                    pINChI_Prev = pINChI;
                    mult = 0; /* we do not know whether the item is empty */
                }
            }
        }

        return ( strbuf->nUsedLength - nUsedLength0 );
    }
    */
    // END INCHI C FUNCTION: str_AuxInvIsoSp3

    fn clone_inchi(
        heap: &SourceHeap,
        pointer: SourceMutPointer<crate::source_types::INChI>,
    ) -> Result<Option<crate::source_types::INChI>, SourceHeapError> {
        if pointer.is_null() {
            Ok(None)
        } else {
            Ok(Some(
                heap.slice(pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            ))
        }
    }

    fn clone_stereo(
        heap: &SourceHeap,
        inchi: Option<&crate::source_types::INChI>,
        isotopic: bool,
    ) -> Result<Option<crate::source_types::INChI_Stereo>, SourceHeapError> {
        let Some(inchi) = inchi else { return Ok(None) };
        let pointer = if isotopic { inchi.StereoIsotopic } else { inchi.Stereo };
        if pointer.is_null() {
            Ok(None)
        } else {
            Ok(Some(
                heap.slice(pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            ))
        }
    }

    fn stereo_equal(
        heap: &SourceHeap,
        left: Option<&crate::source_types::INChI_Stereo>,
        left_mode: i32,
        right: Option<&crate::source_types::INChI_Stereo>,
        right_mode: i32,
    ) -> Result<bool, SourceHeapError> {
        Ok(left.is_some() && Eql_INChI_Stereo(heap, left, left_mode, right, right_mode, 0)? != 0)
    }

    fn equivalence_pointer(heap: &mut SourceHeap, value: &str) -> Result<SourceConstPointer<i8>, SourceHeapError> {
        heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain(std::iter::once(0)).collect())
            .map(SourceMutPointer::as_const)
    }

    fn has_isotopic_layer(
        heap: &SourceHeap,
        inchi: Option<&crate::source_types::INChI>,
        include_possible_h: bool,
    ) -> Result<bool, SourceHeapError> {
        let Some(inchi) = inchi else { return Ok(false) };
        if inchi
            .nNumberOfIsotopicAtoms
            .wrapping_add(inchi.nNumberOfIsotopicTGroups)
            > 0
        {
            return Ok(true);
        }
        if !include_possible_h || inchi.nPossibleLocationsOfIsotopicH.is_null() {
            return Ok(false);
        }
        Ok(heap
            .slice(inchi.nPossibleLocationsOfIsotopicH.as_const())?
            .first()
            .is_some_and(|value| *value > 1))
    }

    fn output_stereo(
        heap: &mut SourceHeap,
        inchi: Option<&crate::source_types::INChI>,
        multiplier: i32,
        star: SourceConstPointer<i8>,
        string_buffer: &mut INCHI_IOS_STRING,
        taut_mode: i32,
        overflow: &mut i32,
        equivalence_flush: bool,
    ) -> Result<(), SourceHeapError> {
        let Some(inchi) = inchi.filter(|inchi| {
            inchi.nNumberOfAtoms != 0
                && inchi
                    .nNumberOfIsotopicAtoms
                    .wrapping_add(inchi.nNumberOfIsotopicTGroups)
                    > 0
        }) else {
            return Ok(());
        };
        let Some(stereo) = clone_stereo(heap, Some(inchi), true)? else {
            return Ok(());
        };
        if stereo.nNumberOfStereoCenters > 0 && (equivalence_flush || stereo.nCompInv2Abs != 0) {
            MakeMult(heap, multiplier, star, string_buffer, 0, overflow)?;
            MakeStereoString(
                heap,
                if equivalence_flush {
                    stereo.nNumber.as_const()
                } else {
                    stereo.nNumberInv.as_const()
                },
                SourceConstPointer::null(),
                stereo.t_parityInv.as_const(),
                0,
                stereo.nNumberOfStereoCenters,
                string_buffer,
                taut_mode,
                overflow,
            )?;
        }
        Ok(())
    }

    let initial_length = string_buffer.nUsedLength;
    let component_delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let star = heap.allocate_model_storage(vec![b'*' as i8, 0])?;
    let mut previous = SourceMutPointer::null();
    let mut taut_previous = SourceMutPointer::null();
    let mut multiplier = 0_i32;
    let mut next = 0_i32;
    let mut previous_equivalence: Option<&'static str> = None;
    let mut previous_equivalence_multiplier = 0_i32;
    let mut equal_to_taut_previous = 1_i32;
    let mut index = 0_i32;

    while index <= number_of_components {
        let current = if index < number_of_components {
            selected_inchi(heap, sorted_inchi.as_const().offset(i64::from(index))?, output_type)?
        } else {
            SourceMutPointer::null()
        };
        let taut_current = if second_non_taut_pass != 0 && index < number_of_components {
            selected_inchi(heap, sorted_inchi2.as_const().offset(i64::from(index))?, OUT_T1 as i32)?
        } else {
            SourceMutPointer::null()
        };
        let current_value = clone_inchi(heap, current)?;
        let taut_value = clone_inchi(heap, taut_current)?;
        let current_iso = clone_stereo(heap, current_value.as_ref(), true)?;
        let current_normal = clone_stereo(heap, current_value.as_ref(), false)?;
        let taut_iso = clone_stereo(heap, taut_value.as_ref(), true)?;
        let taut_normal = clone_stereo(heap, taut_value.as_ref(), false)?;
        let mut equivalence = 0_i32;

        if second_non_taut_pass != 0
            && omit_repetitions != 0
            && has_isotopic_layer(heap, current_value.as_ref(), false)?
        {
            if stereo_equal(
                heap,
                current_iso.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
                taut_normal.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
            )? {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iitNONTAUT) as i32;
            } else if stereo_equal(
                heap,
                current_iso.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
                current_normal.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
            )? {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2NONTAUT) as i32;
            } else if stereo_equal(
                heap,
                current_iso.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
                taut_iso.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
            )? {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2ISO) as i32;
            } else if stereo_equal(
                heap,
                current_iso.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
                taut_normal.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
            )? {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2INV) as i32;
            } else if stereo_equal(
                heap,
                current_iso.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
                taut_iso.as_ref(),
                crate::source_types::EQL_SP3 as i32,
            )? {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2INV
                    | crate::source_types::iiEq2ISO) as i32;
            } else if stereo_equal(
                heap,
                current_iso.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
                current_normal.as_ref(),
                crate::source_types::EQL_SP3 as i32,
            )? {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2INV
                    | crate::source_types::iiEq2NONTAUT) as i32;
            } else if stereo_equal(
                heap,
                current_iso.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
                current_iso.as_ref(),
                crate::source_types::EQL_SP3 as i32,
            )? {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2INV
                    | crate::source_types::iiEq2ISO
                    | crate::source_types::iiEq2NONTAUT) as i32;
            }
        } else if second_non_taut_pass == 0
            && omit_repetitions != 0
            && has_isotopic_layer(heap, current_value.as_ref(), true)?
        {
            if stereo_equal(
                heap,
                current_iso.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
                current_normal.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
            )? {
                equivalence = (crate::source_types::iiSTEREO_INV | crate::source_types::iitISO) as i32;
            } else if stereo_equal(
                heap,
                current_iso.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
                current_normal.as_ref(),
                crate::source_types::EQL_SP3 as i32,
            )? {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iiEq2INV) as i32;
            } else if stereo_equal(
                heap,
                current_iso.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
                current_iso.as_ref(),
                crate::source_types::EQL_SP3 as i32,
            )? {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iiEq2INV
                    | crate::source_types::iiEq2ISO) as i32;
            }
        }

        if equivalence != 0 {
            let previous_value = clone_inchi(heap, previous)?;
            if previous_value.as_ref().is_some_and(|inchi| inchi.nNumberOfAtoms != 0) {
                if next != 0 {
                    MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                output_stereo(
                    heap,
                    previous_value.as_ref(),
                    multiplier.wrapping_add(1),
                    star.as_const(),
                    string_buffer,
                    taut_mode,
                    overflow,
                    true,
                )?;
            } else {
                let taut_previous_value = clone_inchi(heap, taut_previous)?;
                if taut_previous_value
                    .as_ref()
                    .is_some_and(|inchi| inchi.nNumberOfAtoms != 0)
                {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                }
            }
            let current_equivalence = crate::source::base::ichiprt1::EquString(equivalence);
            if previous_equivalence_multiplier != 0 && previous_equivalence.is_some() {
                if previous_equivalence == Some(current_equivalence) {
                    previous_equivalence_multiplier = previous_equivalence_multiplier.wrapping_add(1);
                } else {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let pointer = equivalence_pointer(heap, previous_equivalence.unwrap())?;
                    MakeEqStr(heap, pointer, previous_equivalence_multiplier, string_buffer, overflow)?;
                    previous_equivalence = Some(current_equivalence);
                    previous_equivalence_multiplier = 1;
                }
            } else {
                previous_equivalence = Some(current_equivalence);
                previous_equivalence_multiplier = 1;
            }
            previous = SourceMutPointer::null();
            taut_previous = SourceMutPointer::null();
            multiplier = 0;
            equal_to_taut_previous = 1;
        } else if equal_to_taut_previous != 0 {
            if previous_equivalence_multiplier != 0 {
                if next != 0 {
                    MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                let pointer =
                    equivalence_pointer(heap, previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?)?;
                MakeEqStr(heap, pointer, previous_equivalence_multiplier, string_buffer, overflow)?;
                previous_equivalence = None;
                previous_equivalence_multiplier = 0;
            }
            equal_to_taut_previous = 0;
            previous = current;
            taut_previous = taut_current;
            multiplier = 0;
        } else {
            let previous_value = clone_inchi(heap, previous)?;
            let previous_iso = clone_stereo(heap, previous_value.as_ref(), true)?;
            let equal_to_previous = use_multipliers != 0
                && has_isotopic_layer(heap, current_value.as_ref(), false)?
                && has_isotopic_layer(heap, previous_value.as_ref(), false)?
                && stereo_equal(
                    heap,
                    current_iso.as_ref(),
                    crate::source_types::EQL_SP3_INV as i32,
                    previous_iso.as_ref(),
                    crate::source_types::EQL_SP3_INV as i32,
                )?;
            if equal_to_previous {
                multiplier = multiplier.wrapping_add(1);
                index = index.wrapping_add(1);
                continue;
            }
            if next != 0 {
                MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
            }
            next = next.wrapping_add(1);
            output_stereo(
                heap,
                previous_value.as_ref(),
                multiplier.wrapping_add(1),
                star.as_const(),
                string_buffer,
                taut_mode,
                overflow,
                false,
            )?;
            previous = current;
            multiplier = 0;
        }
        index = index.wrapping_add(1);
    }

    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_AuxInvSp3(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    sorted_inchi2: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    taut_mode: i32,
    number_of_components: i32,
    second_non_taut_pass: i32,
    omit_repetitions: i32,
    use_multipliers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:2315 str_AuxInvSp3
    // INCHI✔️❌: complete configured active source frame follows verbatim; FIX_EMPTY_LAYER_BUG and non-release debug blocks are inactive.
    /*
    int str_AuxInvSp3( INCHI_SORT       *pINChISort,
                       INCHI_SORT       *pINChISort2,
                       INCHI_IOS_STRING *strbuf,
                       int              *bOverflow,
                       int              bOutType,
                       int              TAUT_MODE,
                       int              num_components,
                       int              bSecondNonTautPass,
                       int              bOmitRepetitions,
                       int              bUseMulipliers )
    {
        int          i, ii, ii2, nUsedLength0;
        INCHI_SORT   *is, *is2, *is0, *is20;
        INChI        *pINChI, *pINChI_Prev, *pINChI_Taut, *pINChI_Taut_Prev;
        INChI_Stereo *Stereo, *Stereo_Prev, *Stereo_Taut, *Stereo_Taut_Prev;
        int          mult, eq2prev, eq2taut, eq2tautPrev, bNext;
        const char  *pPrevEquStr, *pCurrEquStr;
        int         multPrevEquStr;
        /***************
        inverted sp3
        ****************/
        pINChI_Taut = NULL;
        pINChI_Prev = NULL;
        pINChI_Taut_Prev = NULL;
        mult = 0;
        bNext = 0;
        is = NULL;
        is2 = NULL;
        is0 = pINChISort;
        is20 = bSecondNonTautPass ? pINChISort2 : NULL;
        /* djb-rwth: removing redundant code */
        eq2tautPrev = 1; /* pINChI_Prev (previous pINChI) does not exist */
        pPrevEquStr = NULL; /*, *pCurrEquStr;*/
        multPrevEquStr = 0;
        nUsedLength0 = strbuf->nUsedLength;

        /* For each connected component...    */
        for (i = 0; i <= num_components; i++)
        {

            /* 1st (taut) pass: bOutType=OUT_TN  ; 2nd (non-taut pass) bOutType=OUT_NT */
            pINChI = ( i < num_components && ( is = is0 + i, 0 <= ( ii = GET_II( bOutType, is ) ) ) ) ? is->pINChI[ii] : NULL;
            /*================ compare sp2 to previous =====================*/
            if (bSecondNonTautPass)
            {
                /* component that was output on the 1st pass */
                pINChI_Taut = ( i < num_components && ( is2 = is20 + i, 0 <= ( ii2 = GET_II( OUT_T1, is2 ) ) ) ) ? is2->pINChI[ii2] : NULL;
            }
            eq2taut = 0;
            /*========= if bSecondNonTautPass then compare iso non-taut stereo to other stereo ========*/
            if (bSecondNonTautPass && bOmitRepetitions)
            {
                /* compare non-tautomeric inverted to:
                *   a) tautomeric inverted
                *   b) Inverted(tautomeric)
                *   c) Inverted(non-tautomeric)
                */
                /* a) compare non-tautomeric inverted to tautomeric inverted */
                if (!eq2taut)
                {
                    eq2taut = pINChI && pINChI_Taut &&
                        /* non-taut inverted */          /* taut invertedc */
                        ( Stereo = pINChI->Stereo ) && ( Stereo_Taut = pINChI_Taut->Stereo ) &&
                        Eql_INChI_Stereo( Stereo, EQL_SP3_INV, Stereo_Taut, EQL_SP3_INV, 0 );
                    /* stereo-inv      non-taut =  taut (stereo-inv) */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iitNONTAUT ) : 0;
                }
                /* b) compare non-tautomeric inverted to Inverted(tautomeric stereo) */
                if (!eq2taut)
                {
                    eq2taut = pINChI && pINChI_Taut &&
                        ( Stereo = pINChI->Stereo ) && ( Stereo_Taut = pINChI_Taut->Stereo ) &&
                        Eql_INChI_Stereo( Stereo, EQL_SP3_INV, Stereo_Taut, EQL_SP3, 0 );
                    /* stereo-inv    non-taut =  Inv(taut stereo) */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iitNONTAUT | iiEq2INV ) : 0;
                }
                /* c) compare non-tautomeric inverted to Inverted(non-tautomeric stereo) */
                if (!eq2taut)
                {
                    eq2taut = pINChI &&
                        ( Stereo = pINChI->Stereo ) &&
                        Eql_INChI_Stereo( Stereo, EQL_SP3_INV, Stereo, EQL_SP3, 0 );
                    /* stereo-inv    non-taut =  Inv(non-taut stereo) */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iitNONTAUT | iiEq2INV | iiEq2NONTAUT ) : 0;
                }
            }
            else
                /*========= if not bSecondNonTautPass then compare inv taut stereo to various taut stereo ========*/
                if (!bSecondNonTautPass && bOmitRepetitions)
                {
                    /* compare tautomeric inverted to Invetred(tautomeric) */
                    if (!eq2taut)
                    {
                        eq2taut = pINChI &&
                            ( Stereo = pINChI->Stereo ) &&
                            Eql_INChI_Stereo( Stereo, EQL_SP3_INV, Stereo, EQL_SP3, 0 );
                        /* stereo     isotopic taut =  taut (stereo) */
                        eq2taut = eq2taut ? ( iiSTEREO_INV | iiEq2INV ) : 0;
                    }
                }
            if (eq2taut)
            {
                /* we may be here only in case of the current layer found equal in another layer the same component */
                if (pINChI_Prev && pINChI_Prev->nNumberOfAtoms)
                {
                    /* previous component exists; output it before output the current component */
                    if (bNext++)
                    {
                        MakeDelim( sCompDelim, strbuf, bOverflow );
                    }
                    if (( Stereo_Prev = pINChI_Prev->Stereo ) && Stereo_Prev->nNumberOfStereoCenters > 0)
                    {
                        MakeMult( mult + 1, "*", strbuf, 0, bOverflow );

                        MakeStereoString( Stereo_Prev->nNumber, NULL, Stereo_Prev->t_parityInv,
                            0, Stereo_Prev->nNumberOfStereoCenters,
                            strbuf, TAUT_MODE, bOverflow );
                    }
                }
                else
                    if (pINChI_Taut_Prev && pINChI_Taut_Prev->nNumberOfAtoms)
                    {
                        /* previous non-taut component exists only in taut list */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        /* do not output stereo of non-tautomeric in non-taut layer: it has been output in the main layer */
                    }
                /* we have found another (previously printed) layer of the current component equal to this layer */
                /* output this (current) equivalence mark = EquString(eq2taut) */
                pCurrEquStr = EquString( eq2taut );
                if (multPrevEquStr && pPrevEquStr)
                {
                    if (pCurrEquStr && !strcmp( pCurrEquStr, pPrevEquStr ))
                    {
                        multPrevEquStr++;
                    }
                    else
                    {
                        /* new EqStr is different; output it */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                        pPrevEquStr = pCurrEquStr;
                        multPrevEquStr = 1;
                    }
                }
                else
                {
                    pPrevEquStr = pCurrEquStr;
                    multPrevEquStr = 1;
                }
                pINChI_Prev = NULL; /* pINChI_Prev sp2 does not exist since */
                pINChI_Taut_Prev = NULL; /* pINChI has just been printed */
                mult = 0;
                eq2tautPrev = 1;     /* pINChI_Prev and pINChI_Taut_Prev have already been output */
            }
            else
                if (eq2tautPrev)
                {
                    /* at this point pINChI_Prev does not exist; however, pINChI */
                    /*might have been discovered and it is different from pINChI_Taut */
                    if (multPrevEquStr && pPrevEquStr)
                    {
                        /* new EqStr is different; output it */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                        pPrevEquStr = NULL;
                        multPrevEquStr = 0;
                    }
                    eq2tautPrev = 0;
                    pINChI_Prev = pINChI;
                    pINChI_Taut_Prev = pINChI_Taut;
                    mult = 0;
                }
                else
                {
                    /* current layer is different from previously printed layers of the current component */
                    /* compare the current layer to this layer of the previous component: */
                    /* check whether pINChI and pINChI_Prev have non-zero identical stereo sp2 */
                    /*================ compare iso sp3 to previous =====================*/
                    eq2prev = bUseMulipliers &&
                        pINChI && pINChI_Prev &&
                        /* do both have stereo? */
                        ( Stereo = pINChI->Stereo ) && ( Stereo_Prev = pINChI_Prev->Stereo ) &&
                        /* is their inverted stereo same? */
                        Eql_INChI_Stereo( Stereo, EQL_SP3_INV, Stereo_Prev, EQL_SP3_INV, 0 );
                    if (eq2prev)
                    {
                        mult++; /* mult = (number of non-empty equal items)-1 */
                        continue;
                    }
                    else
                    {
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        if (pINChI_Prev && pINChI_Prev->nNumberOfAtoms)
                        {
                            if (( Stereo_Prev = pINChI_Prev->Stereo ) &&
                                Stereo_Prev->nNumberOfStereoCenters > 0 && Stereo_Prev->nCompInv2Abs)
                            {
                                MakeMult( mult + 1, "*", strbuf, 0, bOverflow );

                                MakeStereoString( Stereo_Prev->nNumberInv, NULL, Stereo_Prev->t_parityInv,
                                    0, Stereo_Prev->nNumberOfStereoCenters,
                                    strbuf, TAUT_MODE, bOverflow );
                            }
                            /* else sp3 info is not present in pINChI_Prev */
                        }
                        else
                            /* do not print pINChI_Prev because it either do not exist of have already been printed */
                            if (bSecondNonTautPass && pINChI_Taut_Prev && pINChI_Taut_Prev->nNumberOfAtoms)
                            {
                                if (( Stereo_Taut_Prev = pINChI_Taut_Prev->Stereo ) &&
                                    Stereo_Taut_Prev->nNumberOfStereoCenters > 0 && Stereo_Taut_Prev->nCompInv2Abs)
                                {
                                    /* since pINChI_Prev does not exist, pINChI_Taut_Prev is non-tautomeric */
                                    /* and it has non-trivial inv sp3 info. It has already been printed in the main section */
                                    /*
                                    MakeDelim( sIdenticalValues, strbuf, bOverflow);
                                    */
                                    ;/* pINChI_Taut_Prev sp3 info was output in the main stereo section */
                                }
                                else
                                {
                                    ; /* pINChI_Taut_Prev exists and has not sp3 info */
                                }
                            }
                    }
                    pINChI_Prev = pINChI;
                    mult = 0; /* we do not know whether the item is empty */
                }
        }

        return ( strbuf->nUsedLength - nUsedLength0 );
    }
    */
    // END INCHI C FUNCTION: str_AuxInvSp3

    fn clone_inchi(
        heap: &SourceHeap,
        pointer: SourceMutPointer<crate::source_types::INChI>,
    ) -> Result<Option<crate::source_types::INChI>, SourceHeapError> {
        if pointer.is_null() {
            Ok(None)
        } else {
            Ok(Some(
                heap.slice(pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            ))
        }
    }

    fn clone_stereo(
        heap: &SourceHeap,
        inchi: Option<&crate::source_types::INChI>,
    ) -> Result<Option<crate::source_types::INChI_Stereo>, SourceHeapError> {
        let Some(inchi) = inchi else { return Ok(None) };
        if inchi.Stereo.is_null() {
            Ok(None)
        } else {
            Ok(Some(
                heap.slice(inchi.Stereo.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            ))
        }
    }

    fn stereo_equal(
        heap: &SourceHeap,
        left: Option<&crate::source_types::INChI_Stereo>,
        left_mode: i32,
        right: Option<&crate::source_types::INChI_Stereo>,
        right_mode: i32,
    ) -> Result<bool, SourceHeapError> {
        Ok(left.is_some() && Eql_INChI_Stereo(heap, left, left_mode, right, right_mode, 0)? != 0)
    }

    fn equivalence_pointer(heap: &mut SourceHeap, value: &str) -> Result<SourceConstPointer<i8>, SourceHeapError> {
        heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain(std::iter::once(0)).collect())
            .map(SourceMutPointer::as_const)
    }

    fn output_stereo(
        heap: &mut SourceHeap,
        inchi: Option<&crate::source_types::INChI>,
        multiplier: i32,
        star: SourceConstPointer<i8>,
        string_buffer: &mut INCHI_IOS_STRING,
        taut_mode: i32,
        overflow: &mut i32,
        equivalence_flush: bool,
    ) -> Result<(), SourceHeapError> {
        let Some(inchi) = inchi.filter(|inchi| inchi.nNumberOfAtoms != 0) else {
            return Ok(());
        };
        let Some(stereo) = clone_stereo(heap, Some(inchi))? else {
            return Ok(());
        };
        if stereo.nNumberOfStereoCenters > 0 && (equivalence_flush || stereo.nCompInv2Abs != 0) {
            MakeMult(heap, multiplier, star, string_buffer, 0, overflow)?;
            let numbers = if equivalence_flush {
                stereo.nNumber.as_const()
            } else {
                stereo.nNumberInv.as_const()
            };
            MakeStereoString(
                heap,
                numbers,
                SourceConstPointer::null(),
                stereo.t_parityInv.as_const(),
                0,
                stereo.nNumberOfStereoCenters,
                string_buffer,
                taut_mode,
                overflow,
            )?;
        }
        Ok(())
    }

    let initial_length = string_buffer.nUsedLength;
    let component_delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let star = heap.allocate_model_storage(vec![b'*' as i8, 0])?;
    let mut previous = SourceMutPointer::null();
    let mut taut_previous = SourceMutPointer::null();
    let mut multiplier = 0_i32;
    let mut next = 0_i32;
    let mut previous_equivalence: Option<&'static str> = None;
    let mut previous_equivalence_multiplier = 0_i32;
    let mut equal_to_taut_previous = 1_i32;
    let mut index = 0_i32;

    while index <= number_of_components {
        let current = if index < number_of_components {
            selected_inchi(heap, sorted_inchi.as_const().offset(i64::from(index))?, output_type)?
        } else {
            SourceMutPointer::null()
        };
        let taut_current = if second_non_taut_pass != 0 && index < number_of_components {
            selected_inchi(heap, sorted_inchi2.as_const().offset(i64::from(index))?, OUT_T1 as i32)?
        } else {
            SourceMutPointer::null()
        };
        let current_value = clone_inchi(heap, current)?;
        let taut_current_value = clone_inchi(heap, taut_current)?;
        let current_stereo = clone_stereo(heap, current_value.as_ref())?;
        let taut_stereo = clone_stereo(heap, taut_current_value.as_ref())?;
        let mut equivalence = 0_i32;
        if second_non_taut_pass != 0 && omit_repetitions != 0 {
            if stereo_equal(
                heap,
                current_stereo.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
                taut_stereo.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
            )? {
                equivalence = (crate::source_types::iiSTEREO_INV | crate::source_types::iitNONTAUT) as i32;
            } else if stereo_equal(
                heap,
                current_stereo.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
                taut_stereo.as_ref(),
                crate::source_types::EQL_SP3 as i32,
            )? {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2INV) as i32;
            } else if stereo_equal(
                heap,
                current_stereo.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
                current_stereo.as_ref(),
                crate::source_types::EQL_SP3 as i32,
            )? {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2INV
                    | crate::source_types::iiEq2NONTAUT) as i32;
            }
        } else if second_non_taut_pass == 0
            && omit_repetitions != 0
            && stereo_equal(
                heap,
                current_stereo.as_ref(),
                crate::source_types::EQL_SP3_INV as i32,
                current_stereo.as_ref(),
                crate::source_types::EQL_SP3 as i32,
            )?
        {
            equivalence = (crate::source_types::iiSTEREO_INV | crate::source_types::iiEq2INV) as i32;
        }

        if equivalence != 0 {
            let previous_value = clone_inchi(heap, previous)?;
            if previous_value.as_ref().is_some_and(|inchi| inchi.nNumberOfAtoms != 0) {
                if next != 0 {
                    MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                output_stereo(
                    heap,
                    previous_value.as_ref(),
                    multiplier.wrapping_add(1),
                    star.as_const(),
                    string_buffer,
                    taut_mode,
                    overflow,
                    true,
                )?;
            } else {
                let taut_previous_value = clone_inchi(heap, taut_previous)?;
                if taut_previous_value
                    .as_ref()
                    .is_some_and(|inchi| inchi.nNumberOfAtoms != 0)
                {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                }
            }
            let current_equivalence = crate::source::base::ichiprt1::EquString(equivalence);
            if previous_equivalence_multiplier != 0 && previous_equivalence.is_some() {
                if previous_equivalence == Some(current_equivalence) {
                    previous_equivalence_multiplier = previous_equivalence_multiplier.wrapping_add(1);
                } else {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let pointer = equivalence_pointer(heap, previous_equivalence.unwrap())?;
                    MakeEqStr(heap, pointer, previous_equivalence_multiplier, string_buffer, overflow)?;
                    previous_equivalence = Some(current_equivalence);
                    previous_equivalence_multiplier = 1;
                }
            } else {
                previous_equivalence = Some(current_equivalence);
                previous_equivalence_multiplier = 1;
            }
            previous = SourceMutPointer::null();
            taut_previous = SourceMutPointer::null();
            multiplier = 0;
            equal_to_taut_previous = 1;
        } else if equal_to_taut_previous != 0 {
            if previous_equivalence_multiplier != 0 {
                if next != 0 {
                    MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                let pointer =
                    equivalence_pointer(heap, previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?)?;
                MakeEqStr(heap, pointer, previous_equivalence_multiplier, string_buffer, overflow)?;
                previous_equivalence = None;
                previous_equivalence_multiplier = 0;
            }
            equal_to_taut_previous = 0;
            previous = current;
            taut_previous = taut_current;
            multiplier = 0;
        } else {
            let previous_value = clone_inchi(heap, previous)?;
            let previous_stereo = clone_stereo(heap, previous_value.as_ref())?;
            let equal_to_previous = use_multipliers != 0
                && current_value.is_some()
                && previous_value.is_some()
                && stereo_equal(
                    heap,
                    current_stereo.as_ref(),
                    crate::source_types::EQL_SP3_INV as i32,
                    previous_stereo.as_ref(),
                    crate::source_types::EQL_SP3_INV as i32,
                )?;
            if equal_to_previous {
                multiplier = multiplier.wrapping_add(1);
                index = index.wrapping_add(1);
                continue;
            }
            if next != 0 {
                MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
            }
            next = next.wrapping_add(1);
            output_stereo(
                heap,
                previous_value.as_ref(),
                multiplier.wrapping_add(1),
                star.as_const(),
                string_buffer,
                taut_mode,
                overflow,
                false,
            )?;
            previous = current;
            multiplier = 0;
        }
        index = index.wrapping_add(1);
    }

    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_AuxInvSp3Numb(
    heap: &mut SourceHeap,
    _canonical_globals: SourceMutPointer<crate::source_types::CANON_GLOBALS>,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    _sorted_inchi2: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    taut_mode: i32,
    number_of_components: i32,
    second_non_taut_pass: i32,
    omit_repetitions: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:2584 str_AuxInvSp3Numb
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int str_AuxInvSp3Numb( CANON_GLOBALS    *pCG,
                           INCHI_SORT       *pINChISort,
                           INCHI_SORT       *pINChISort2,
                           INCHI_IOS_STRING *strbuf,
                           int              *bOverflow,
                           int              bOutType,
                           int              TAUT_MODE,
                           int              num_components,
                           int              bSecondNonTautPass,
                           int               bOmitRepetitions )
    {
        int          i, ii, ii2, nUsedLength0;
        INCHI_SORT   *is, *is0 /*, *is2*/;
        INChI        *pINChI, *pINChI_Taut;
        INChI_Aux    *pINChI_Aux, *pINChI_Aux_Taut; /* djb-rwth: removing redundant variables */
        INChI_Stereo *Stereo, *Stereo_Taut;
        int          eq2taut, bNext;
        const char  *pPrevEquStr, *pCurrEquStr;
        int         multPrevEquStr;
        /**************************************************
        * specificity of numbering: there is no previous *
        * component because no repetition is possible    *
        **************************************************/
        pINChI = NULL;
        pINChI_Taut = NULL;
        pINChI_Aux = NULL;
        pINChI_Aux_Taut = NULL;
        /* djb-rwth: removing redundant code */
        bNext = 0;
        is = NULL;
        is0 = pINChISort;
        /*is2         = bSecondNonTautPass? pINChISort2 : NULL;*/
        /* djb-rwth: removing redundant code */
        pPrevEquStr = NULL; /*, *pCurrEquStr;*/
        multPrevEquStr = 0;
        nUsedLength0 = strbuf->nUsedLength;

        /* For each connected component...    */
        for (i = 0; i < num_components; i++)
        {

            /* 1st (taut) pass: bOutType=OUT_TN  ; 2nd (non-taut pass) bOutType=OUT_NT */
            is = is0 + i;
            pINChI = ( 0 <= ( ii = GET_II( bOutType, is ) ) ) ? is->pINChI[ii] : NULL;
            pINChI_Aux = pINChI ? is->pINChI_Aux[ii] : NULL;
            /*================ to compare to previously printed =====================*/
            if (bSecondNonTautPass)
            {
                /* component that was printed on the 1st pass */
                pINChI_Taut = ( 0 <= ( ii2 = GET_II( OUT_T1, is ) ) ) ? is->pINChI[ii2] : NULL;
                pINChI_Aux_Taut = pINChI_Taut ? is->pINChI_Aux[ii2] : NULL;
            }

            eq2taut = 0;
            /*========= if bSecondNonTautPass then compare inv non-taut stereo to other stereo ========*/
            if (bSecondNonTautPass && bOmitRepetitions &&
                pINChI && ( Stereo = pINChI->Stereo ) && Stereo->nCompInv2Abs)
            {
                /* compare non-tautomeric inverted stereo numbering to:
                *   a) tautomeric numbering
                *   b) non-tautomeric numbering
                *   c) tautomeric inverted stereo numbering
                */
                /* a) compare non-tautomeric inverted stereo numbering to tautomeric numbering */
                if (!eq2taut)
                {
                    eq2taut = pINChI_Taut &&
                        Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_INV, pINChI_Aux_Taut, EQL_NUM );
                    /* stereo-inv     numbering  non-taut =  taut numbering */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iiNUMB | iitNONTAUT ) : 0;
                }
                /* b) compare non-tautomeric inverted stereo numbering to non-tautomeric numbering */
                if (!eq2taut)
                {
                    eq2taut =
                        Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_INV, pINChI_Aux, EQL_NUM );
                    /* stereo-inv     numb.     non-taut =  non-taut numbering */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iiNUMB | iitNONTAUT | iiEq2NONTAUT ) : 0;
                }
                /* c) compare non-tautomeric inverted stereo numbering to tautomeric inverted stereo numbering */
                if (!eq2taut)
                {
                    eq2taut = pINChI_Taut &&
                        ( Stereo_Taut = pINChI_Taut->Stereo ) && Stereo_Taut->nCompInv2Abs &&
                        Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_INV, pINChI_Aux_Taut, EQL_NUM_INV );
                    /* stereo-inv     numb.     non-taut =  taut inv stereo numbering */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iiNUMB | iitNONTAUT | iiEq2INV ) : 0;
                }
            }
            else
            {
                /*========= if not bSecondNonTautPass then compare inv taut stereo numb to taut numb ========*/
                if (!bSecondNonTautPass && bOmitRepetitions &&
                    pINChI && ( Stereo = pINChI->Stereo ) && Stereo->nCompInv2Abs)
                {
                    /* compare tautomeric inverted stereo numbering to tautomeric numbering */
                    if (!eq2taut)
                    {
                        eq2taut =
                            Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_INV, pINChI_Aux, EQL_NUM );
                        /* stereo-inv     numbering  (taut) =  taut numbering */
                        eq2taut = eq2taut ? ( iiSTEREO_INV | iiNUMB ) : 0;
                    }
                }
            }
            if (eq2taut)
            {
                /* we have found another (previously printed) layer of the current component equal to this layer */
                /* output this (current) equivalence mark = EquString(eq2taut) */
                pCurrEquStr = EquString( eq2taut );
                if (multPrevEquStr && pPrevEquStr)
                {
                    if (pCurrEquStr && !strcmp( pCurrEquStr, pPrevEquStr ))
                    {
                        multPrevEquStr++;
                    }
                    else
                    {
                        /* new EqStr is different; output it */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                        pPrevEquStr = pCurrEquStr;
                        multPrevEquStr = 1;
                    }
                }
                else
                {
                    pPrevEquStr = pCurrEquStr;
                    multPrevEquStr = 1;
                }
            }
            else
            {
                /* current layer is different from previously printed layers of the current component */
                if (multPrevEquStr && pPrevEquStr)
                {
                    /* new EqStr is different; output it */
                    if (bNext++)
                    {
                        MakeDelim( sCompDelim, strbuf, bOverflow );
                    }
                    MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                    pPrevEquStr = NULL;
                    multPrevEquStr = 0;
                }
                if (bNext++)
                {
                    MakeDelim( sCompDelim, strbuf, bOverflow );
                }
                if (pINChI && pINChI_Aux && pINChI_Aux->nNumberOfAtoms &&
                    ( Stereo = pINChI->Stereo ) && Stereo->nNumberOfStereoCenters &&
                    Stereo->nCompInv2Abs && pINChI_Aux->nOrigAtNosInCanonOrdInv)
                {
                    MakeCtString( pCG, pINChI_Aux->nOrigAtNosInCanonOrdInv,
                        pINChI_Aux->nNumberOfAtoms, 0, NULL, 0,
                        strbuf, TAUT_MODE, bOverflow );
                }
                /* else inv stereo info is not present in pINChI */
            }
        }

        if (multPrevEquStr && pPrevEquStr)
        {
            /* the new EqStr of the last item has not been printed; output it now */
            if (bNext++)
            {
                MakeDelim( sCompDelim, strbuf, bOverflow );
            }
            MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
            pPrevEquStr = NULL;
            /* djb-rwth: removing redundant code */
        }

        return ( strbuf->nUsedLength - nUsedLength0 );
    }
    */
    // END INCHI C FUNCTION: str_AuxInvSp3Numb

    let initial_length = string_buffer.nUsedLength;
    if sorted_inchi.is_null() {
        return Ok(initial_length);
    }
    let delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let mut previous_equivalence: Option<&'static str> = None;
    let mut previous_multiplier = 0_i32;
    let mut next = 0_i32;

    for index in 0..number_of_components {
        let sort = heap
            .slice(sorted_inchi.as_const().offset(i64::from(index))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let selected = get_ii(heap, output_type, &sort)?;
        let current_inchi_pointer = selected.map_or(SourceMutPointer::null(), |slot| sort.pINChI[slot]);
        let current_aux_pointer = selected.map_or(SourceMutPointer::null(), |slot| sort.pINChI_Aux[slot]);
        let taut = if second_non_taut_pass != 0 {
            get_ii(heap, OUT_T1 as i32, &sort)?
        } else {
            None
        };
        let taut_inchi_pointer = taut.map_or(SourceMutPointer::null(), |slot| sort.pINChI[slot]);
        let taut_aux_pointer = taut.map_or(SourceMutPointer::null(), |slot| sort.pINChI_Aux[slot]);
        let current_inchi = if current_inchi_pointer.is_null() {
            None
        } else {
            Some(
                heap.slice(current_inchi_pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            )
        };
        let taut_inchi = if taut_inchi_pointer.is_null() {
            None
        } else {
            Some(
                heap.slice(taut_inchi_pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            )
        };
        let current_aux = if current_aux_pointer.is_null() {
            None
        } else {
            Some(
                heap.slice(current_aux_pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            )
        };
        let taut_aux = if taut_aux_pointer.is_null() {
            None
        } else {
            Some(
                heap.slice(taut_aux_pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            )
        };
        let current_stereo = if let Some(inchi) = current_inchi.as_ref() {
            if inchi.Stereo.is_null() {
                None
            } else {
                Some(
                    heap.slice(inchi.Stereo.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone(),
                )
            }
        } else {
            None
        };
        let taut_stereo = if let Some(inchi) = taut_inchi.as_ref() {
            if inchi.Stereo.is_null() {
                None
            } else {
                Some(
                    heap.slice(inchi.Stereo.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone(),
                )
            }
        } else {
            None
        };
        let has_inversion = current_stereo.as_ref().is_some_and(|stereo| stereo.nCompInv2Abs != 0);
        let equal_numbers = |heap: &SourceHeap,
                             left: Option<&crate::source_types::INChI_Aux>,
                             left_mode: i32,
                             right: Option<&crate::source_types::INChI_Aux>,
                             right_mode: i32|
         -> Result<bool, SourceHeapError> {
            Ok(crate::source::base::ichiprt2::Eql_INChI_Aux_Num(heap, left, left_mode, right, right_mode)? != 0)
        };
        let mut equivalence = 0_i32;
        if second_non_taut_pass != 0 && omit_repetitions != 0 && current_inchi.is_some() && has_inversion {
            if taut_inchi.is_some()
                && equal_numbers(
                    heap,
                    current_aux.as_ref(),
                    crate::source_types::EQL_NUM_INV as i32,
                    taut_aux.as_ref(),
                    crate::source_types::EQL_NUM as i32,
                )?
            {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iiNUMB
                    | crate::source_types::iitNONTAUT) as i32;
            } else if equal_numbers(
                heap,
                current_aux.as_ref(),
                crate::source_types::EQL_NUM_INV as i32,
                current_aux.as_ref(),
                crate::source_types::EQL_NUM as i32,
            )? {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iiNUMB
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2NONTAUT) as i32;
            } else if taut_inchi.is_some()
                && taut_stereo.as_ref().is_some_and(|stereo| stereo.nCompInv2Abs != 0)
                && equal_numbers(
                    heap,
                    current_aux.as_ref(),
                    crate::source_types::EQL_NUM_INV as i32,
                    taut_aux.as_ref(),
                    crate::source_types::EQL_NUM_INV as i32,
                )?
            {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iiNUMB
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2INV) as i32;
            }
        } else if second_non_taut_pass == 0
            && omit_repetitions != 0
            && current_inchi.is_some()
            && has_inversion
            && equal_numbers(
                heap,
                current_aux.as_ref(),
                crate::source_types::EQL_NUM_INV as i32,
                current_aux.as_ref(),
                crate::source_types::EQL_NUM as i32,
            )?
        {
            equivalence = (crate::source_types::iiSTEREO_INV | crate::source_types::iiNUMB) as i32;
        }

        if equivalence != 0 {
            let current_equivalence = crate::source::base::ichiprt1::EquString(equivalence);
            if previous_multiplier != 0 {
                if previous_equivalence == Some(current_equivalence) {
                    previous_multiplier = previous_multiplier.wrapping_add(1);
                } else {
                    if next != 0 {
                        MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let value = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let value =
                        heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain([0]).collect())?;
                    MakeEqStr(heap, value.as_const(), previous_multiplier, string_buffer, overflow)?;
                    previous_equivalence = Some(current_equivalence);
                    previous_multiplier = 1;
                }
            } else {
                previous_equivalence = Some(current_equivalence);
                previous_multiplier = 1;
            }
        } else {
            if previous_multiplier != 0 {
                if next != 0 {
                    MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                let value = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
                let value = heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain([0]).collect())?;
                MakeEqStr(heap, value.as_const(), previous_multiplier, string_buffer, overflow)?;
                previous_equivalence = None;
                previous_multiplier = 0;
            }
            if next != 0 {
                MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
            }
            next = next.wrapping_add(1);
            if let (Some(inchi), Some(aux), Some(stereo)) =
                (current_inchi.as_ref(), current_aux.as_ref(), current_stereo.as_ref())
            {
                if inchi.nNumberOfAtoms != 0
                    && aux.nNumberOfAtoms != 0
                    && stereo.nNumberOfStereoCenters != 0
                    && stereo.nCompInv2Abs != 0
                    && !aux.nOrigAtNosInCanonOrdInv.is_null()
                {
                    MakeCtString(
                        heap,
                        aux.nOrigAtNosInCanonOrdInv,
                        aux.nNumberOfAtoms,
                        0,
                        SourceConstPointer::null(),
                        0,
                        string_buffer,
                        taut_mode,
                        overflow,
                    )?;
                }
            }
        }
    }

    if previous_multiplier != 0 {
        if next != 0 {
            MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
        }
        let value = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
        let value = heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain([0]).collect())?;
        MakeEqStr(heap, value.as_const(), previous_multiplier, string_buffer, overflow)?;
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_AuxInvIsoSp3Numb(
    heap: &mut SourceHeap,
    _canonical_globals: SourceMutPointer<crate::source_types::CANON_GLOBALS>,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    _sorted_inchi2: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    taut_mode: i32,
    number_of_components: i32,
    second_non_taut_pass: i32,
    omit_repetitions: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:3575 str_AuxInvIsoSp3Numb
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int str_AuxInvIsoSp3Numb( CANON_GLOBALS    *pCG,
                              INCHI_SORT       *pINChISort,
                              INCHI_SORT       *pINChISort2,
                              INCHI_IOS_STRING *strbuf,
                              int              *bOverflow,
                              int              bOutType,
                              int              TAUT_MODE,
                              int              num_components,
                              int              bSecondNonTautPass,
                              int              bOmitRepetitions )
    {
        int          i, ii, ii2, nUsedLength0;
        INCHI_SORT   *is, *is0 /*, *is2*/;
        INChI        *pINChI, *pINChI_Taut;
        INChI_Aux    *pINChI_Aux, *pINChI_Aux_Taut; /* djb-rwth: removing redundant variables */
        INChI_Stereo *Stereo, *Stereo_Taut;
        int          eq2taut, bNext;
        const char  *pPrevEquStr, *pCurrEquStr;
        int         multPrevEquStr;
        /**************************************************
        * specificity of numbering: there is no previous *
        * component because no repetition is possible    *
        **************************************************/
        pINChI = NULL;
        pINChI_Taut = NULL;
        pINChI_Aux = NULL;
        pINChI_Aux_Taut = NULL;
        /* djb-rwth: removing redundant code */
        bNext = 0;
        is = NULL;
        /* is2         = NULL;*/
        is0 = pINChISort;
        /* is20        = bSecondNonTautPass? pINChISort2 : NULL;*/
        /* djb-rwth: removing redundant code */
        pPrevEquStr = NULL; /*, *pCurrEquStr;*/
        multPrevEquStr = 0;
        nUsedLength0 = strbuf->nUsedLength;

        /* For each connected component...    */
        for (i = 0; i < num_components; i++)
        {

            /* 1st (taut) pass: bOutType=OUT_TN  ; 2nd (non-taut pass) bOutType=OUT_NT */
            is = is0 + i;
            pINChI = ( 0 <= ( ii = GET_II( bOutType, is ) ) ) ? is->pINChI[ii] : NULL;
            pINChI_Aux = pINChI ? is->pINChI_Aux[ii] : NULL;
            /*================ to compare to previously printed =====================*/
            if (bSecondNonTautPass)
            {
                /* component that was printed on the 1st pass */
                pINChI_Taut = ( 0 <= ( ii2 = GET_II( OUT_T1, is ) ) ) ? is->pINChI[ii2] : NULL;
                pINChI_Aux_Taut = pINChI_Taut ? is->pINChI_Aux[ii2] : NULL;
            }
            eq2taut = 0;
            /*========= if bSecondNonTautPass then compare iso non-taut stereo to other stereo ========*/
            if (bSecondNonTautPass && bOmitRepetitions && pINChI && pINChI_Aux && pINChI_Aux->bIsIsotopic &&
                ( Stereo = pINChI->StereoIsotopic ) && Stereo->nCompInv2Abs &&
                pINChI_Aux->nNumberOfAtoms > 0 && pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv)
            {
                /* compare isotopic non-tautomeric inverted stereo numbering to:
                *   a) tautomeric numbering
                *   b) non-tautomeric numbering
                *   c) *tautomeric isotopic numbering
                *   d) *non-tautomeric isotopic numbering
                *   e) tautomeric inverted stereo numbering
                *   f) *non-tautomeric inverted stereo numbering
                *   g) tautomeric isotopic inverted stereo numbering
                */
                /* a) compare isotopic non-tautomeric inverted stereo numbering to tautomeric numbering */
                if (!eq2taut)
                {
                    eq2taut = pINChI_Taut &&
                        Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_INV | EQL_NUM_ISO, pINChI_Aux_Taut, EQL_NUM );
                    /* stereo-inv   isotopic numbering  non-taut =  taut numbering */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iiNUMB | iitNONTAUT ) : 0;
                }
                /* b) compare isotopic non-tautomeric inverted stereo numbering to non-tautomeric numbering */
                if (!eq2taut)
                {
                    eq2taut =
                        Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_INV | EQL_NUM_ISO, pINChI_Aux, EQL_NUM );
                    /* stereo-inv    isotopic   numb.    non-taut =  non-taut numbering */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iiNUMB | iitNONTAUT | iiEq2NONTAUT ) : 0;
                }
                /* c) compare isotopic non-tautomeric inverted stereo numbering to tautomeric isotopic numbering */
                if (!eq2taut)
                {
                    eq2taut = pINChI_Taut &&
                        Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_INV | EQL_NUM_ISO, pINChI_Aux_Taut, EQL_NUM_ISO );
                    /* stereo-inv   isotopic   numb.     non-taut =  taut iso numbering */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iiNUMB | iitNONTAUT | iiEq2ISO ) : 0;
                }
                /* d) compare isotopic non-tautomeric inverted stereo numbering to non-tautomeric isotopic numbering */
                if (!eq2taut)
                {
                    eq2taut = Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_INV | EQL_NUM_ISO, pINChI_Aux, EQL_NUM_ISO );
                    /* stereo-inv   isotopic   numb.     non-taut =  non-taut isotopic numbering */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iiNUMB | iitNONTAUT | iiEq2NONTAUT | iiEq2ISO ) : 0;
                }
                /* e) compare isotopic non-tautomeric inverted stereo numbering to tautomeric inverted stereo numbering */
                if (!eq2taut)
                {
                    eq2taut = pINChI_Taut && pINChI_Aux_Taut &&
                        ( Stereo_Taut = pINChI_Taut->Stereo ) && Stereo_Taut->nCompInv2Abs &&
                        Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_INV | EQL_NUM_ISO, pINChI_Aux_Taut, EQL_NUM_INV );
                    /* stereo-inv   isotopic numbering  non-taut =  stereo-inv taut numbering */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iiNUMB | iitNONTAUT | iiEq2INV ) : 0;
                }
                /* f) compare isotopic non-tautomeric inverted stereo numbering to non-tautomeric inverted stereo numbering */
                if (!eq2taut)
                {
                    eq2taut =
                        ( Stereo_Taut = pINChI->StereoIsotopic ) && Stereo_Taut->nCompInv2Abs &&
                        Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_INV | EQL_NUM_ISO, pINChI_Aux, EQL_NUM_INV );
                    /* stereo-inv   isotopic numbering  non-taut =  stereo-inv non-taut numbering */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iiNUMB | iitNONTAUT | iiEq2INV | iiEq2NONTAUT ) : 0;
                }

                /* g) compare isotopic non-tautomeric inverted stereo numbering to tautomeric isotopic inverted stereo numbering */
                if (!eq2taut)
                {
                    eq2taut = pINChI_Taut &&
                        ( Stereo_Taut = pINChI_Taut->StereoIsotopic ) && Stereo_Taut->nCompInv2Abs &&
                        Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_INV | EQL_NUM_ISO, pINChI_Aux_Taut, EQL_NUM_INV | EQL_NUM_ISO );
                    /* stereo-inv   isotopic numbering  non-taut =  stereo-inv iso taut numbering */
                    eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iiNUMB | iitNONTAUT | iiEq2INV | iiEq2ISO ) : 0;
                }
            }
            else
            {
                /*========= if not bSecondNonTautPass then compare inv taut stereo numb to taut numb ========*/
                if (!bSecondNonTautPass && bOmitRepetitions && pINChI && pINChI_Aux && pINChI_Aux->bIsIsotopic &&
                    ( Stereo = pINChI->StereoIsotopic ) && Stereo->nCompInv2Abs &&
                    pINChI_Aux->nNumberOfAtoms > 0 && pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv)
                {
                    /* compare isotopic tautomeric inverted stereo numbering to:
                    *   a) tautomeric numbering
                    *   b) tautomeric isotopic numbering
                    *   c) tautomeric inverted stereo numbering
                    */
                    /* a) compare isotopic tautomeric inverted stereo numbering to tautomeric numbering */
                    if (!eq2taut)
                    {
                        eq2taut = Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_INV | EQL_NUM_ISO, pINChI_Aux, EQL_NUM );
                        /* stereo-inv   isotopic numbering  (taut) =  taut numbering */
                        eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iiNUMB ) : 0;
                    }
                    /* b) compare isotopic tautomeric inverted stereo numbering to tautomeric isotopic numbering */
                    if (!eq2taut)
                    {
                        eq2taut = Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_INV | EQL_NUM_ISO, pINChI_Aux, EQL_NUM_ISO );
                        /* stereo-inv   isotopic numbering(taut) =  isotopic taut numbering */
                        eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iiNUMB | iiEq2ISO ) : 0;
                    }
                    /* b) compare isotopic tautomeric inverted stereo numbering to tautomeric inverted stereo numbering */
                    if (!eq2taut)
                    {
                        eq2taut = pINChI->Stereo && Stereo->nCompInv2Abs &&
                            Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_INV | EQL_NUM_ISO, pINChI_Aux, EQL_NUM_INV ); /* djb-rwth: removing redundant code */
                        /* stereo-inv   isotopic numbering  (taut) =  taut stereo-inv numbering */
                        eq2taut = eq2taut ? ( iiSTEREO_INV | iitISO | iiNUMB | iiEq2INV ) : 0;
                    }
                }
            }
            if (eq2taut)
            {
                /* we have found another (previously printed) layer of the current component equal to this layer */
                /* output this (current) equivalence mark = EquString(eq2taut) */
                pCurrEquStr = EquString( eq2taut );
                if (multPrevEquStr && pPrevEquStr)
                {
                    if (pCurrEquStr && !strcmp( pCurrEquStr, pPrevEquStr ))
                    {
                        multPrevEquStr++;
                    }
                    else
                    {
                        /* new EqStr is different; output it */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                        pPrevEquStr = pCurrEquStr;
                        multPrevEquStr = 1;
                    }
                }
                else
                {
                    pPrevEquStr = pCurrEquStr;
                    multPrevEquStr = 1;
                }
            }
            else
            {
                /* current layer is different from previously printed layers of the current component */
                if (multPrevEquStr && pPrevEquStr)
                {
                    /* new EqStr is different; output it */
                    if (bNext++)
                    {
                        MakeDelim( sCompDelim, strbuf, bOverflow );
                    }
                    MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                    pPrevEquStr = NULL;
                    multPrevEquStr = 0;
                }
                if (bNext++)
                {
                    MakeDelim( sCompDelim, strbuf, bOverflow );
                }
                if (pINChI && pINChI_Aux && pINChI_Aux->bIsIsotopic && pINChI_Aux->nNumberOfAtoms &&
                    ( Stereo = pINChI->StereoIsotopic ) && Stereo->nNumberOfStereoCenters &&
                    Stereo->nCompInv2Abs && pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv)
                {
                    MakeCtString( pCG, pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv,
                        pINChI_Aux->nNumberOfAtoms, 0, NULL, 0,
                        strbuf, TAUT_MODE, bOverflow );
                }
                /* else isotopic inv stereo info is not present in pINChI */
            }
        }
        if (multPrevEquStr && pPrevEquStr)
        {
            /* the new EqStr of the last item has not been printed; output it now */
            if (bNext++)
            {
                MakeDelim( sCompDelim, strbuf, bOverflow );
            }
            MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
            pPrevEquStr = NULL;
            /* djb-rwth: removing redundant code */
        }

        return ( strbuf->nUsedLength - nUsedLength0 );
    }
    */
    // END INCHI C FUNCTION: str_AuxInvIsoSp3Numb

    let initial_length = string_buffer.nUsedLength;
    let delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let mut previous_equivalence: Option<&'static str> = None;
    let mut previous_multiplier = 0_i32;
    let mut next = 0_i32;
    let iso_inverted_mode = (crate::source_types::EQL_NUM_INV | crate::source_types::EQL_NUM_ISO) as i32;

    for index in 0..number_of_components {
        let sort = heap
            .slice(sorted_inchi.as_const().offset(i64::from(index))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let selected = get_ii(heap, output_type, &sort)?;
        let current_inchi_pointer = selected.map_or(SourceMutPointer::null(), |slot| sort.pINChI[slot]);
        let current_aux_pointer = selected.map_or(SourceMutPointer::null(), |slot| sort.pINChI_Aux[slot]);
        let taut = if second_non_taut_pass != 0 {
            get_ii(heap, OUT_T1 as i32, &sort)?
        } else {
            None
        };
        let taut_inchi_pointer = taut.map_or(SourceMutPointer::null(), |slot| sort.pINChI[slot]);
        let taut_aux_pointer = taut.map_or(SourceMutPointer::null(), |slot| sort.pINChI_Aux[slot]);
        let load_inchi = |pointer: SourceMutPointer<crate::source_types::INChI>|
         -> Result<Option<crate::source_types::INChI>, SourceHeapError> {
            if pointer.is_null() {
                Ok(None)
            } else {
                Ok(Some(
                    heap.slice(pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone(),
                ))
            }
        };
        let load_aux = |pointer: SourceMutPointer<crate::source_types::INChI_Aux>|
         -> Result<Option<crate::source_types::INChI_Aux>, SourceHeapError> {
            if pointer.is_null() {
                Ok(None)
            } else {
                Ok(Some(
                    heap.slice(pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone(),
                ))
            }
        };
        let current_inchi = load_inchi(current_inchi_pointer)?;
        let taut_inchi = load_inchi(taut_inchi_pointer)?;
        let current_aux = load_aux(current_aux_pointer)?;
        let taut_aux = load_aux(taut_aux_pointer)?;
        let load_stereo = |inchi: Option<&crate::source_types::INChI>,
                           isotopic: bool|
         -> Result<Option<crate::source_types::INChI_Stereo>, SourceHeapError> {
            let Some(inchi) = inchi else { return Ok(None) };
            let pointer = if isotopic { inchi.StereoIsotopic } else { inchi.Stereo };
            if pointer.is_null() {
                Ok(None)
            } else {
                Ok(Some(
                    heap.slice(pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone(),
                ))
            }
        };
        let current_iso_stereo = load_stereo(current_inchi.as_ref(), true)?;
        let taut_normal_stereo = load_stereo(taut_inchi.as_ref(), false)?;
        let taut_iso_stereo = load_stereo(taut_inchi.as_ref(), true)?;
        let eligible = current_inchi.is_some()
            && current_aux.as_ref().is_some_and(|aux| {
                aux.bIsIsotopic != 0 && aux.nNumberOfAtoms > 0 && !aux.nIsotopicOrigAtNosInCanonOrdInv.is_null()
            })
            && current_iso_stereo
                .as_ref()
                .is_some_and(|stereo| stereo.nCompInv2Abs != 0);
        let equal_numbers = |left: Option<&crate::source_types::INChI_Aux>,
                             left_mode: i32,
                             right: Option<&crate::source_types::INChI_Aux>,
                             right_mode: i32|
         -> Result<bool, SourceHeapError> {
            Ok(crate::source::base::ichiprt2::Eql_INChI_Aux_Num(heap, left, left_mode, right, right_mode)? != 0)
        };
        let mut equivalence = 0_i32;
        if second_non_taut_pass != 0 && omit_repetitions != 0 && eligible {
            if taut_inchi.is_some()
                && equal_numbers(
                    current_aux.as_ref(),
                    iso_inverted_mode,
                    taut_aux.as_ref(),
                    crate::source_types::EQL_NUM as i32,
                )?
            {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iiNUMB
                    | crate::source_types::iitNONTAUT) as i32;
            } else if equal_numbers(
                current_aux.as_ref(),
                iso_inverted_mode,
                current_aux.as_ref(),
                crate::source_types::EQL_NUM as i32,
            )? {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iiNUMB
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2NONTAUT) as i32;
            } else if taut_inchi.is_some()
                && equal_numbers(
                    current_aux.as_ref(),
                    iso_inverted_mode,
                    taut_aux.as_ref(),
                    crate::source_types::EQL_NUM_ISO as i32,
                )?
            {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iiNUMB
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2ISO) as i32;
            } else if equal_numbers(
                current_aux.as_ref(),
                iso_inverted_mode,
                current_aux.as_ref(),
                crate::source_types::EQL_NUM_ISO as i32,
            )? {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iiNUMB
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2NONTAUT
                    | crate::source_types::iiEq2ISO) as i32;
            } else if taut_inchi.is_some()
                && taut_normal_stereo
                    .as_ref()
                    .is_some_and(|stereo| stereo.nCompInv2Abs != 0)
                && equal_numbers(
                    current_aux.as_ref(),
                    iso_inverted_mode,
                    taut_aux.as_ref(),
                    crate::source_types::EQL_NUM_INV as i32,
                )?
            {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iiNUMB
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2INV) as i32;
            } else if current_iso_stereo
                .as_ref()
                .is_some_and(|stereo| stereo.nCompInv2Abs != 0)
                && equal_numbers(
                    current_aux.as_ref(),
                    iso_inverted_mode,
                    current_aux.as_ref(),
                    crate::source_types::EQL_NUM_INV as i32,
                )?
            {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iiNUMB
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2INV
                    | crate::source_types::iiEq2NONTAUT) as i32;
            } else if taut_inchi.is_some()
                && taut_iso_stereo.as_ref().is_some_and(|stereo| stereo.nCompInv2Abs != 0)
                && equal_numbers(
                    current_aux.as_ref(),
                    iso_inverted_mode,
                    taut_aux.as_ref(),
                    iso_inverted_mode,
                )?
            {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iiNUMB
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2INV
                    | crate::source_types::iiEq2ISO) as i32;
            }
        } else if second_non_taut_pass == 0 && omit_repetitions != 0 && eligible {
            if equal_numbers(
                current_aux.as_ref(),
                iso_inverted_mode,
                current_aux.as_ref(),
                crate::source_types::EQL_NUM as i32,
            )? {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iiNUMB) as i32;
            } else if equal_numbers(
                current_aux.as_ref(),
                iso_inverted_mode,
                current_aux.as_ref(),
                crate::source_types::EQL_NUM_ISO as i32,
            )? {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iiNUMB
                    | crate::source_types::iiEq2ISO) as i32;
            } else if current_inchi.as_ref().is_some_and(|inchi| !inchi.Stereo.is_null())
                && equal_numbers(
                    current_aux.as_ref(),
                    iso_inverted_mode,
                    current_aux.as_ref(),
                    crate::source_types::EQL_NUM_INV as i32,
                )?
            {
                equivalence = (crate::source_types::iiSTEREO_INV
                    | crate::source_types::iitISO
                    | crate::source_types::iiNUMB
                    | crate::source_types::iiEq2INV) as i32;
            }
        }

        if equivalence != 0 {
            let current_equivalence = crate::source::base::ichiprt1::EquString(equivalence);
            if previous_multiplier != 0 {
                if previous_equivalence == Some(current_equivalence) {
                    previous_multiplier = previous_multiplier.wrapping_add(1);
                } else {
                    if next != 0 {
                        MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let value = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let value =
                        heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain([0]).collect())?;
                    MakeEqStr(heap, value.as_const(), previous_multiplier, string_buffer, overflow)?;
                    previous_equivalence = Some(current_equivalence);
                    previous_multiplier = 1;
                }
            } else {
                previous_equivalence = Some(current_equivalence);
                previous_multiplier = 1;
            }
        } else {
            if previous_multiplier != 0 {
                if next != 0 {
                    MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                let value = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
                let value = heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain([0]).collect())?;
                MakeEqStr(heap, value.as_const(), previous_multiplier, string_buffer, overflow)?;
                previous_equivalence = None;
                previous_multiplier = 0;
            }
            if next != 0 {
                MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
            }
            next = next.wrapping_add(1);
            if let (Some(aux), Some(stereo)) = (current_aux.as_ref(), current_iso_stereo.as_ref())
                && aux.bIsIsotopic != 0
                && aux.nNumberOfAtoms != 0
                && stereo.nNumberOfStereoCenters != 0
                && stereo.nCompInv2Abs != 0
                && !aux.nIsotopicOrigAtNosInCanonOrdInv.is_null()
            {
                MakeCtString(
                    heap,
                    aux.nIsotopicOrigAtNosInCanonOrdInv,
                    aux.nNumberOfAtoms,
                    0,
                    SourceConstPointer::null(),
                    0,
                    string_buffer,
                    taut_mode,
                    overflow,
                )?;
            }
        }
    }

    if previous_multiplier != 0 {
        if next != 0 {
            MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
        }
        let value = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
        let value = heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain([0]).collect())?;
        MakeEqStr(heap, value.as_const(), previous_multiplier, string_buffer, overflow)?;
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_AuxIsoNumb(
    heap: &mut SourceHeap,
    _canonical_globals: SourceMutPointer<crate::source_types::CANON_GLOBALS>,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    _sorted_inchi2: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    taut_mode: i32,
    number_of_components: i32,
    second_non_taut_pass: i32,
    omit_repetitions: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:2769 str_AuxIsoNumb
    // INCHI✔️❌: complete configured source frame follows verbatim; AUX_ISO_NUMB_BUG_FIX is the active branch.
    /*
    int str_AuxIsoNumb( CANON_GLOBALS    *pCG,
                        INCHI_SORT       *pINChISort,
                        INCHI_SORT       *pINChISort2,
                        INCHI_IOS_STRING *strbuf,
                        int              *bOverflow,
                        int              bOutType,
                        int              TAUT_MODE,
                        int              num_components,
                        int              bSecondNonTautPass,
                        int              bOmitRepetitions )
    {
        int          i, ii, ii2, nUsedLength0;
        INCHI_SORT   *is, *is0 /*, *is2*/;
        /* djb-rwth: removing redundant variables */
        INChI_Aux    *pINChI_Aux, *pINChI_Aux_Taut; /* djb-rwth: removing redundant variables */
        int          eq2taut, bNext;
        const char  *pPrevEquStr, *pCurrEquStr;
        int         multPrevEquStr;
        /**************************************************
        * specificity of numbering: there is no previous *
        * component because no repetition is possible    *
        **************************************************/
        /* djb-rwth: removing redundant code */
        pINChI_Aux = NULL;
        pINChI_Aux_Taut = NULL;
        /* djb-rwth: removing redundant code */
        bNext = 0;
        is = NULL;
        is0 = pINChISort;
        /*is2         = bSecondNonTautPass? pINChISort2 : NULL;*/
        /* djb-rwth: removing redundant code */
        pPrevEquStr = NULL; /*, *pCurrEquStr;*/
        multPrevEquStr = 0;
        nUsedLength0 = strbuf->nUsedLength;

        /* For each connected component...    */
        for (i = 0; i < num_components; i++)
        {

            /* 1st (taut) pass: bOutType=OUT_TN  ; 2nd (non-taut pass) bOutType=OUT_NT */
            is = is0 + i;
            pINChI_Aux = ( i < num_components && 0 <= ( ii = GET_II( bOutType, is ) ) ) ? is->pINChI_Aux[ii] : NULL;
            /*================ to compare to previously printed =====================*/
            if (bSecondNonTautPass)
            {
                pINChI_Aux_Taut = ( 0 <= ( ii2 = GET_II( OUT_T1, is ) ) ) ? is->pINChI_Aux[ii2] : NULL;
            }
            eq2taut = 0;
            /*========= if bSecondNonTautPass then compare iso non-taut numb to other numb ========*/
            if (bSecondNonTautPass && bOmitRepetitions && pINChI_Aux && pINChI_Aux->bIsIsotopic)
            {
                /* compare non-tautomeric isotopic numbering to:
                *   a) tautomeric numbering
                *   b) non-tautomeric numbering
                *   c) tautomeric isotopic numbering
                */
                /* a) compare non-tautomeric isotopic numbering to tautomeric numbering */
                if (!eq2taut)
                {
                    eq2taut = Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_ISO, pINChI_Aux_Taut, EQL_NUM );
                    /* numbering  non-taut isotopic =  taut numbering */
                    eq2taut = eq2taut ? ( iiNUMB | iitNONTAUT | iitISO ) : 0;
                }
                /* b) compare non-tautomeric isotopic numbering to non-tautomeric numbering */
                if (!eq2taut)
                {
                    eq2taut = Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_ISO, pINChI_Aux, EQL_NUM );
                    /* numbering  non-taut isotopic =  non-taut numbering */
                    eq2taut = eq2taut ? ( iiNUMB | iitNONTAUT | iitISO | iiEq2NONTAUT ) : 0;
                }
                /* c) compare non-tautomeric isotopic numbering to tautomeric isotopic numbering */
                if (!eq2taut)
                {
                    eq2taut = Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_ISO, pINChI_Aux_Taut, EQL_NUM_ISO );
                    /* numbering  non-taut isotopic =  taut isotopic numbering */
                    eq2taut = eq2taut ? ( iiNUMB | iitNONTAUT | iitISO | iiEq2ISO ) : 0;
                }
            }
            else
            {
                /*    2011-10-28
                Fix bug in src:cano mapping of  atoms printed to AuxInfo
                Reported by Sandor Mark on 2011-10-25 in inchi-discuss
                See http://sourceforge.net/p/inchi/mailman/message/28292914/ also
                */
    #define AUX_ISO_NUMB_BUG_FIX
    #ifdef AUX_ISO_NUMB_BUG_FIX
                /* Bug-fixed version */
                /*========= if not bSecondNonTautPass then compare isotopic taut numbering to taut numbering ========*/
                if (!bSecondNonTautPass && bOmitRepetitions && pINChI_Aux && pINChI_Aux->bIsIsotopic)
                {
                    /* compare tautomeric isotopic numbering to tautomeric non-isotopic numbering */
                    if (!eq2taut)
                    {
                        eq2taut = Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_ISO, pINChI_Aux, EQL_NUM );
                        /*  numbering isotopic (taut) =  taut numbering */
                        eq2taut = eq2taut ? ( iiNUMB | iitISO ) : 0;
                    }
                }
    #else
                /* Original (buggy) version */

                /*========= if not bSecondNonTautPass then compare inv taut stereo numb to taut numb ========*/
                if (!bSecondNonTautPass && bOmitRepetitions && pINChI_Aux && pINChI_Aux->bIsIsotopic)
                {
                    /* compare tautomeric isotopic numbering to tautomeric non-isotopic numbering */
                    if (!eq2taut)
                    {
                        eq2taut = Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM_ISO, pINChI_Aux, EQL_NUM_ISO );
                        /* stereo-inv     numbering  (taut) =  taut numbering */
                        eq2taut = eq2taut ? ( iiSTEREO_INV | iiNUMB ) : 0;
                    }
                }
    #endif
            }

            if (eq2taut)
            {
                /* we have found another (previously printed) layer of the current component equal to this layer */
                /* output this (current) equivalence mark = EquString(eq2taut) */
                pCurrEquStr = EquString( eq2taut );
                if (multPrevEquStr && pPrevEquStr)
                {
                    if (pCurrEquStr && !strcmp( pCurrEquStr, pPrevEquStr ))
                    {
                        multPrevEquStr++;
                    }
                    else
                    {
                        /* new EqStr is different; output it */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                        pPrevEquStr = pCurrEquStr;
                        multPrevEquStr = 1;
                    }
                }
                else
                {
                    pPrevEquStr = pCurrEquStr;
                    multPrevEquStr = 1;
                }
            }
            else
            {
                /* current layer is different from previously printed layers of the current component */
                if (multPrevEquStr && pPrevEquStr)
                {
                    /* new EqStr is different; output it */
                    if (bNext++)
                    {
                        MakeDelim( sCompDelim, strbuf, bOverflow );
                    }
                    MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                    pPrevEquStr = NULL;
                    multPrevEquStr = 0;
                }
                if (bNext++)
                {
                    MakeDelim( sCompDelim, strbuf, bOverflow );
                }
                if (pINChI_Aux && pINChI_Aux->nNumberOfAtoms && pINChI_Aux->bIsIsotopic &&
                    pINChI_Aux->nIsotopicOrigAtNosInCanonOrd)
                {
                    MakeCtString( pCG, pINChI_Aux->nIsotopicOrigAtNosInCanonOrd,
                        pINChI_Aux->nNumberOfAtoms,
                        0, NULL, 0, strbuf, TAUT_MODE, bOverflow );
                }
                /* else isotopic numbering is not present in pINChI */
            }
        }

        if (multPrevEquStr && pPrevEquStr)
        {
            /* the new EqStr of the last item has not been printed; output it now */
            if (bNext++)
            {
                MakeDelim( sCompDelim, strbuf, bOverflow );
            }
            MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
            pPrevEquStr = NULL;
            /* djb-rwth: removing redundant code */
        }

        return ( strbuf->nUsedLength - nUsedLength0 );
    }
    */
    // END INCHI C FUNCTION: str_AuxIsoNumb

    let initial_length = string_buffer.nUsedLength;
    let delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let mut previous_equivalence: Option<&'static str> = None;
    let mut previous_multiplier = 0_i32;
    let mut next = 0_i32;

    for index in 0..number_of_components {
        let sort = heap
            .slice(sorted_inchi.as_const().offset(i64::from(index))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let selected = get_ii(heap, output_type, &sort)?;
        let current_aux_pointer = selected.map_or(SourceMutPointer::null(), |slot| sort.pINChI_Aux[slot]);
        let taut_aux_pointer = if second_non_taut_pass != 0 {
            get_ii(heap, OUT_T1 as i32, &sort)?.map_or(SourceMutPointer::null(), |slot| sort.pINChI_Aux[slot])
        } else {
            SourceMutPointer::null()
        };
        let current_aux = if current_aux_pointer.is_null() {
            None
        } else {
            Some(
                heap.slice(current_aux_pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            )
        };
        let taut_aux = if taut_aux_pointer.is_null() {
            None
        } else {
            Some(
                heap.slice(taut_aux_pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            )
        };
        let equal_numbers = |left: Option<&crate::source_types::INChI_Aux>,
                             left_mode: i32,
                             right: Option<&crate::source_types::INChI_Aux>,
                             right_mode: i32|
         -> Result<bool, SourceHeapError> {
            Ok(crate::source::base::ichiprt2::Eql_INChI_Aux_Num(heap, left, left_mode, right, right_mode)? != 0)
        };
        let mut equivalence = 0_i32;
        if second_non_taut_pass != 0
            && omit_repetitions != 0
            && current_aux.as_ref().is_some_and(|aux| aux.bIsIsotopic != 0)
        {
            if equal_numbers(
                current_aux.as_ref(),
                crate::source_types::EQL_NUM_ISO as i32,
                taut_aux.as_ref(),
                crate::source_types::EQL_NUM as i32,
            )? {
                equivalence = (crate::source_types::iiNUMB
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iitISO) as i32;
            } else if equal_numbers(
                current_aux.as_ref(),
                crate::source_types::EQL_NUM_ISO as i32,
                current_aux.as_ref(),
                crate::source_types::EQL_NUM as i32,
            )? {
                equivalence = (crate::source_types::iiNUMB
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iitISO
                    | crate::source_types::iiEq2NONTAUT) as i32;
            } else if equal_numbers(
                current_aux.as_ref(),
                crate::source_types::EQL_NUM_ISO as i32,
                taut_aux.as_ref(),
                crate::source_types::EQL_NUM_ISO as i32,
            )? {
                equivalence = (crate::source_types::iiNUMB
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iitISO
                    | crate::source_types::iiEq2ISO) as i32;
            }
        } else if second_non_taut_pass == 0
            && omit_repetitions != 0
            && current_aux.as_ref().is_some_and(|aux| aux.bIsIsotopic != 0)
            && equal_numbers(
                current_aux.as_ref(),
                crate::source_types::EQL_NUM_ISO as i32,
                current_aux.as_ref(),
                crate::source_types::EQL_NUM as i32,
            )?
        {
            equivalence = (crate::source_types::iiNUMB | crate::source_types::iitISO) as i32;
        }

        if equivalence != 0 {
            let current_equivalence = crate::source::base::ichiprt1::EquString(equivalence);
            if previous_multiplier != 0 {
                if previous_equivalence == Some(current_equivalence) {
                    previous_multiplier = previous_multiplier.wrapping_add(1);
                } else {
                    if next != 0 {
                        MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let value = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let value =
                        heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain([0]).collect())?;
                    MakeEqStr(heap, value.as_const(), previous_multiplier, string_buffer, overflow)?;
                    previous_equivalence = Some(current_equivalence);
                    previous_multiplier = 1;
                }
            } else {
                previous_equivalence = Some(current_equivalence);
                previous_multiplier = 1;
            }
        } else {
            if previous_multiplier != 0 {
                if next != 0 {
                    MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                let value = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
                let value = heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain([0]).collect())?;
                MakeEqStr(heap, value.as_const(), previous_multiplier, string_buffer, overflow)?;
                previous_equivalence = None;
                previous_multiplier = 0;
            }
            if next != 0 {
                MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
            }
            next = next.wrapping_add(1);
            if let Some(aux) = current_aux.as_ref()
                && aux.nNumberOfAtoms != 0
                && aux.bIsIsotopic != 0
                && !aux.nIsotopicOrigAtNosInCanonOrd.is_null()
            {
                MakeCtString(
                    heap,
                    aux.nIsotopicOrigAtNosInCanonOrd,
                    aux.nNumberOfAtoms,
                    0,
                    SourceConstPointer::null(),
                    0,
                    string_buffer,
                    taut_mode,
                    overflow,
                )?;
            }
        }
    }

    if previous_multiplier != 0 {
        if next != 0 {
            MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
        }
        let value = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
        let value = heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain([0]).collect())?;
        MakeEqStr(heap, value.as_const(), previous_multiplier, string_buffer, overflow)?;
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_AuxNumb(
    heap: &mut SourceHeap,
    _canonical_globals: SourceMutPointer<crate::source_types::CANON_GLOBALS>,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    _sorted_inchi2: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    taut_mode: i32,
    number_of_components: i32,
    second_non_taut_pass: i32,
    omit_repetitions: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:3818 str_AuxNumb
    // INCHI✔️❌: int str_AuxNumb( CANON_GLOBALS    *pCG,
    // INCHI✔️❌:                  INCHI_SORT       *pINChISort,
    // INCHI✔️❌:                  INCHI_SORT       *pINChISort2,
    // INCHI✔️❌:                  INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                  int              *bOverflow,
    // INCHI✔️❌:                  int              bOutType,
    // INCHI✔️❌:                  int              TAUT_MODE,
    // INCHI✔️❌:                  int              num_components,
    // INCHI✔️❌:                  int              bSecondNonTautPass,
    // INCHI✔️❌:                  int              bOmitRepetitions )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int          i, ii, ii2, nUsedLength0;
    // INCHI✔️❌:     INCHI_SORT   *is, *is0 /*, *is2*/;
    // INCHI✔️❌:     INChI        *pINChI, *pINChI_Taut = NULL;
    // INCHI✔️❌:     INChI_Aux    *pINChI_Aux, *pINChI_Aux_Taut = NULL;
    // INCHI✔️❌:     int          eq2taut, bNext;
    // INCHI✔️❌:     const char  *pPrevEquStr, *pCurrEquStr;
    // INCHI✔️❌:     int         multPrevEquStr;
    // INCHI✔️❌:     bNext = 0;
    // INCHI✔️❌:     /*is2         = bSecondNonTautPass? pINChISort2 : NULL;*/
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     pPrevEquStr = NULL; /*, *pCurrEquStr;*/
    // INCHI✔️❌:     multPrevEquStr = 0;
    // INCHI✔️❌:     is = NULL;
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!( is0 = pINChISort ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return nUsedLength0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* For each connected component...    */
    // INCHI✔️❌:     for (i = 0; i < num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* 1st (taut) pass: bOutType=OUT_TN  ; 2nd (non-taut pass) bOutType=OUT_NT */
    // INCHI✔️❌:         is = is0 + i;
    // INCHI✔️❌:         pINChI = ( 0 <= ( ii = GET_II( bOutType, is ) ) ) ? is->pINChI[ii] : NULL;
    // INCHI✔️❌:         pINChI_Aux = pINChI ? is->pINChI_Aux[ii] : NULL;
    // INCHI✔️❌:         /*================ to compare to previously printed =====================*/
    // INCHI✔️❌:         if (bSecondNonTautPass)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* component that was printed on the 1st pass */
    // INCHI✔️❌:             pINChI_Taut = ( 0 <= ( ii2 = GET_II( OUT_T1, is ) ) ) ? is->pINChI[ii2] : NULL;
    // INCHI✔️❌:             pINChI_Aux_Taut = pINChI_Taut ? is->pINChI_Aux[ii2] : NULL;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         eq2taut = 0;
    // INCHI✔️❌:         /*========= if bSecondNonTautPass then compare iso non-taut stereo to other stereo ========*/
    // INCHI✔️❌:         if (bSecondNonTautPass && bOmitRepetitions && pINChI && pINChI_Aux && pINChI_Aux->nNumberOfAtoms > 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* compare non-tautomeric numbering to:
    // INCHI✔️❌:             *   a) tautomeric numbering
    // INCHI✔️❌:             */
    // INCHI✔️❌:             /* a) compare non-tautomeric numbering to tautomeric numbering */
    // INCHI✔️❌:             if (!eq2taut)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 eq2taut = pINChI_Taut && !pINChI_Taut->bDeleted &&
    // INCHI✔️❌:                     Eql_INChI_Aux_Num( pINChI_Aux, EQL_NUM, pINChI_Aux_Taut, EQL_NUM );
    // INCHI✔️❌:                 /* numbering  non-taut =  taut numbering */
    // INCHI✔️❌:                 eq2taut = eq2taut ? ( iiNUMB | iitNONTAUT ) : 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (eq2taut)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* we have found another (previously printed) layer of the current component equal to this layer */
    // INCHI✔️❌:             /* output this (current) equivalence mark = EquString(eq2taut) */
    // INCHI✔️❌:             pCurrEquStr = EquString( eq2taut );
    // INCHI✔️❌:             if (multPrevEquStr && pPrevEquStr)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (pCurrEquStr && !strcmp( pCurrEquStr, pPrevEquStr ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     multPrevEquStr++;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* new EqStr is different; output it */
    // INCHI✔️❌:                     if (bNext++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
    // INCHI✔️❌:                     pPrevEquStr = pCurrEquStr;
    // INCHI✔️❌:                     multPrevEquStr = 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 pPrevEquStr = pCurrEquStr;
    // INCHI✔️❌:                 multPrevEquStr = 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* current layer is different from previously printed layers of the current component */
    // INCHI✔️❌:             if (multPrevEquStr && pPrevEquStr)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* new EqStr is different; output it */
    // INCHI✔️❌:                 if (bNext++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
    // INCHI✔️❌:                 pPrevEquStr = NULL;
    // INCHI✔️❌:                 multPrevEquStr = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (bNext++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (pINChI && pINChI_Aux && pINChI_Aux->nNumberOfAtoms)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 MakeCtString( pCG, pINChI_Aux->nOrigAtNosInCanonOrd,
    // INCHI✔️❌:                     pINChI_Aux->nNumberOfAtoms, 0, NULL, 0,
    // INCHI✔️❌:                     strbuf, TAUT_MODE, bOverflow );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (multPrevEquStr && pPrevEquStr)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* the new EqStr of the last item has not been printed; output it now */
    // INCHI✔️❌:         if (bNext++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             MakeDelim( sCompDelim, strbuf, bOverflow );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
    // INCHI✔️❌:         pPrevEquStr = NULL;
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( strbuf->nUsedLength - nUsedLength0 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: str_AuxNumb

    let initial_length = string_buffer.nUsedLength;
    if sorted_inchi.is_null() {
        return Ok(initial_length);
    }
    let delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let mut previous_equivalence: Option<&'static str> = None;
    let mut previous_multiplier = 0_i32;
    let mut next = 0_i32;

    for index in 0..number_of_components {
        let index_pointer = sorted_inchi.as_const().offset(i64::from(index))?;
        let sort = heap
            .slice(index_pointer)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let current_index = get_ii(heap, output_type, sort)?;
        let current_inchi = current_index.map_or(SourceMutPointer::null(), |value| sort.pINChI[value]);
        let current_aux = current_index.map_or(SourceMutPointer::null(), |value| sort.pINChI_Aux[value]);
        let taut_index = if second_non_taut_pass != 0 {
            get_ii(heap, OUT_T1 as i32, sort)?
        } else {
            None
        };
        let taut_inchi = taut_index.map_or(SourceMutPointer::null(), |value| sort.pINChI[value]);
        let taut_aux = taut_index.map_or(SourceMutPointer::null(), |value| sort.pINChI_Aux[value]);

        let equal_to_taut = if second_non_taut_pass != 0 && omit_repetitions != 0 {
            if current_inchi.is_null() || current_aux.is_null() || taut_inchi.is_null() {
                false
            } else {
                let taut_inchi_value = heap
                    .slice(taut_inchi.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if taut_inchi_value.bDeleted != 0 {
                    false
                } else {
                    let current_aux_value = heap
                        .slice(current_aux.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if current_aux_value.nNumberOfAtoms <= 0 {
                        false
                    } else {
                        crate::source::base::ichiprt2::Eql_INChI_Aux_Num(
                            &*heap,
                            Some(current_aux_value),
                            crate::source_types::EQL_NUM as i32,
                            if taut_aux.is_null() {
                                None
                            } else {
                                Some(
                                    heap.slice(taut_aux.as_const())?
                                        .first()
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                                )
                            },
                            crate::source_types::EQL_NUM as i32,
                        )? != 0
                    }
                }
            }
        } else {
            false
        };
        if equal_to_taut {
            let current_equivalence = crate::source::base::ichiprt1::EquString(
                (crate::source_types::iiNUMB | crate::source_types::iitNONTAUT) as i32,
            );
            if previous_multiplier != 0 {
                if previous_equivalence == Some(current_equivalence) {
                    previous_multiplier = previous_multiplier.wrapping_add(1);
                } else {
                    if next != 0 {
                        MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let previous = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let previous_pointer =
                        heap.allocate_model_storage(previous.bytes().map(|byte| byte as i8).chain([0]).collect())?;
                    MakeEqStr(
                        heap,
                        previous_pointer.as_const(),
                        previous_multiplier,
                        string_buffer,
                        overflow,
                    )?;
                    previous_equivalence = Some(current_equivalence);
                    previous_multiplier = 1;
                }
            } else {
                previous_equivalence = Some(current_equivalence);
                previous_multiplier = 1;
            }
        } else {
            if previous_multiplier != 0 {
                if next != 0 {
                    MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                let equivalence = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
                let equivalence_pointer =
                    heap.allocate_model_storage(equivalence.bytes().map(|byte| byte as i8).chain([0]).collect())?;
                MakeEqStr(
                    heap,
                    equivalence_pointer.as_const(),
                    previous_multiplier,
                    string_buffer,
                    overflow,
                )?;
                previous_equivalence = None;
                previous_multiplier = 0;
            }
            if next != 0 {
                MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
            }
            next = next.wrapping_add(1);
            if !current_inchi.is_null() && !current_aux.is_null() {
                let aux = heap
                    .slice(current_aux.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if aux.nNumberOfAtoms != 0 {
                    MakeCtString(
                        heap,
                        aux.nOrigAtNosInCanonOrd,
                        aux.nNumberOfAtoms,
                        0,
                        SourceConstPointer::null(),
                        0,
                        string_buffer,
                        taut_mode,
                        overflow,
                    )?;
                }
            }
        }
    }
    if previous_multiplier != 0 {
        if next != 0 {
            MakeDelim(heap, delimiter.as_const(), string_buffer, overflow)?;
        }
        next = next.wrapping_add(1);
        let equivalence = previous_equivalence.ok_or(SourceHeapError::PointerOutOfBounds)?;
        let equivalence_pointer =
            heap.allocate_model_storage(equivalence.bytes().map(|byte| byte as i8).chain([0]).collect())?;
        MakeEqStr(
            heap,
            equivalence_pointer.as_const(),
            previous_multiplier,
            string_buffer,
            overflow,
        )?;
    }
    let _ = next;
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_IsoSp3(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    sorted_inchi2: SourceMutPointer<INCHI_SORT>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    output_type: i32,
    taut_mode: i32,
    number_of_components: i32,
    mut relative_or_racemic: i32,
    second_non_taut_pass: i32,
    omit_repetitions: i32,
    use_multipliers: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:1768 str_IsoSp3
    // INCHI✔️❌: complete configured source frame follows verbatim; inactive FIX_EMPTY_LAYER_BUG blocks remain in the official file.
    /*
    int str_IsoSp3( INCHI_SORT       *pINChISort,
                    INCHI_SORT       *pINChISort2,
                    INCHI_IOS_STRING *strbuf,
                    int              *bOverflow,
                    int              bOutType,
                    int              TAUT_MODE,
                    int              num_components,
                    int              bRelRac,
                    int              bSecondNonTautPass,
                    int              bOmitRepetitions,
                    int              bUseMulipliers )
    {
        int          i, ii, ii2, nUsedLength0;
        INCHI_SORT   *is, *is2, *is0, *is20;
        INChI        *pINChI, *pINChI_Prev, *pINChI_Taut, *pINChI_Taut_Prev;
        INChI_Stereo *Stereo, *Stereo_Prev, *Stereo_Taut, *Stereo_Taut_Prev;
        int          mult, eq2prev, eq2taut, eq2tautPrev, bNext;
        const char  *pPrevEquStr, *pCurrEquStr;
        int         multPrevEquStr;
        pINChI_Taut = NULL;
        pINChI_Prev = NULL;
        pINChI_Taut_Prev = NULL;
        mult = 0;
        bNext = 0;
        is = NULL;
        is2 = NULL;
        is0 = pINChISort;
        is20 = bSecondNonTautPass ? pINChISort2 : NULL;
        /* djb-rwth: removing redundant code */
        eq2tautPrev = 1; /* pINChI_Prev (previous pINChI) does not exist */
        pPrevEquStr = NULL; /*, *pCurrEquStr;*/
        multPrevEquStr = 0;
    #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    #else
        bRelRac = 0;
    #endif
        nUsedLength0 = strbuf->nUsedLength;

        /* For each connected component...    */
        for (i = 0; i <= num_components; i++)
        {

            /* 1st (taut) pass: bOutType=OUT_TN  ; 2nd (non-taut pass) bOutType=OUT_NT */
            pINChI = ( i < num_components && ( is = is0 + i, 0 <= ( ii = GET_II( bOutType, is ) ) ) ) ? is->pINChI[ii] : NULL;
            /*================ compare sp2 to previous =====================*/
            if (bSecondNonTautPass)
            {
                /* component that was output on the 1st pass */
                pINChI_Taut = ( i < num_components && ( is2 = is20 + i, 0 <= ( ii2 = GET_II( OUT_T1, is2 ) ) ) ) ? is2->pINChI[ii2] : NULL;
            }
            eq2taut = 0;
            /*========= if bSecondNonTautPass then compare iso non-taut stereo to other stereo ========*/
            if (bSecondNonTautPass && bOmitRepetitions)
            {
                /* compare non-tautomeric isotopic to:
                *   a) non-tautomeric non-isotopic
                *   b) tautomeric non-isotopic
                *   c) tautomeric isotopic
                */
                /* a) compare non-tautomeric isotopic to non-tautomeric non-isotopic */
                if (!eq2taut)
                {
                    eq2taut = pINChI && /* non-taut isotopic */                  /* non-taut non-isotopic */
                        ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Taut = pINChI->Stereo ) &&
                        Eql_INChI_Stereo( Stereo, EQL_SP3, Stereo_Taut, EQL_SP3, bRelRac );
                    /* stereo     isotopic non-taut =  non-taut (stereo) */
                    eq2taut = eq2taut ? ( iiSTEREO | iitISO | iitNONTAUT | iiEq2NONTAUT ) : 0;
                }
                /* b) compare non-tautomeric isotopic to tautomeric non-isotopic */
                if (!eq2taut)
                {
                    eq2taut = pINChI && pINChI_Taut &&
                        /* non-taut isotopic */                  /* taut non-isotopic */
                        ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Taut = pINChI_Taut->Stereo ) &&
                        Eql_INChI_Stereo( Stereo, EQL_SP3, Stereo_Taut, EQL_SP3, bRelRac );
                    /* stereo     isotopic non-taut =  taut (stereo) */
                    eq2taut = eq2taut ? ( iiSTEREO | iitISO | iitNONTAUT ) : 0;
                }
                /* c) compare non-tautomeric isotopic to tautomeric isotopic */
                if (!eq2taut && bSecondNonTautPass && bOmitRepetitions)
                {
                    eq2taut = pINChI && pINChI_Taut &&
                        ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Taut = pINChI_Taut->StereoIsotopic ) &&
                        Eql_INChI_Stereo( Stereo, EQL_SP3, Stereo_Taut, EQL_SP3, bRelRac );
                    /* stereo     isotopic non-taut =  isotopic taut (stereo) */
                    eq2taut = eq2taut ? ( iiSTEREO | iitISO | iitNONTAUT | iiEq2ISO ) : 0;
                }
            }
            else
            {
                /*========= if not bSecondNonTautPass then compare iso taut stereo to non-iso taut ========*/
                if (!bSecondNonTautPass && bOmitRepetitions)
                {
                    /* compare tautomeric isotopic to tautomeric non-isotopic */
                    if (!eq2taut)
                    {
                        eq2taut = pINChI &&
                            /* taut isotopic */                  /* taut non-isotopic */
                            ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Taut = pINChI->Stereo ) &&
                            Eql_INChI_Stereo( Stereo, EQL_SP3, Stereo_Taut, EQL_SP3, bRelRac );
                        /* stereo     isotopic taut =  taut (stereo) */
                        eq2taut = eq2taut ? ( iiSTEREO | iitISO ) : 0;
                    }
                }
            }
            if (eq2taut)
            {
                /* we may be here only in case of the current layer found equal in another layer the same component */
                if (pINChI_Prev && pINChI_Prev->nNumberOfAtoms)
                {
                    /* previous component exists; output it before output the current component */
                    if (bNext++)
                    {
                        MakeDelim( sCompDelim, strbuf, bOverflow );
                    }
                    if (( Stereo_Prev = pINChI_Prev->StereoIsotopic ) && Stereo_Prev->nNumberOfStereoCenters > bRelRac)
                    {
                        MakeMult( mult + 1, "*", strbuf, 0, bOverflow );

                        MakeStereoString( Stereo_Prev->nNumber, NULL, Stereo_Prev->t_parity,
                            0, Stereo_Prev->nNumberOfStereoCenters,
                            strbuf, TAUT_MODE, bOverflow );
                    }
                }
                else
                {
                    if (pINChI_Taut_Prev && pINChI_Taut_Prev->nNumberOfAtoms)
                    {
                        /* previous non-taut component exists only in taut list */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        /* do not output stereo of non-tautomeric in non-taut layer: it has been output in the main layer */
                    }
                }
                /* we have found another (previously printed) layer of the current component equal to this layer */
                /* output this (current) equivalence mark = EquString(eq2taut) */
                pCurrEquStr = EquString( eq2taut );
                if (multPrevEquStr && pPrevEquStr)
                {
                    if (pCurrEquStr && !strcmp( pCurrEquStr, pPrevEquStr ))
                    {
                        multPrevEquStr++;
                    }
                    else
                    {
                        /* new EqStr is different; output it */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                        pPrevEquStr = pCurrEquStr;
                        multPrevEquStr = 1;
                    }
                }
                else
                {
                    pPrevEquStr = pCurrEquStr;
                    multPrevEquStr = 1;
                }
                pINChI_Prev = NULL; /* pINChI_Prev sp2 does not exist since */
                pINChI_Taut_Prev = NULL; /* pINChI has just been printed */
                mult = 0;
                eq2tautPrev = 1;     /* pINChI_Prev and pINChI_Taut_Prev have already been output */
            }
            else
            {
                if (eq2tautPrev)
                {
                    /* at this point pINChI_Prev does not exist; however, pINChI */
                    /*might have been discovered and it is different from pINChI_Taut */
                    if (multPrevEquStr && pPrevEquStr)
                    {
                        /* new EqStr is different; output it */
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        MakeEqStr( pPrevEquStr, multPrevEquStr, strbuf, bOverflow );
                        pPrevEquStr = NULL;
                        multPrevEquStr = 0;
                    }
                    eq2tautPrev = 0;
                    pINChI_Prev = pINChI;
                    pINChI_Taut_Prev = pINChI_Taut;
                    mult = 0;
                }
                else
                {
                    /* current layer is different from previously printed layers of the current component */
                    /* compare the current layer to this layer of the previous component: */
                    /* check whether pINChI and pINChI_Prev have non-zero identical stereo sp2 */
                    /*================ compare iso sp3 to previous =====================*/
                    eq2prev = bUseMulipliers && pINChI && pINChI_Prev &&
                        ( Stereo = pINChI->StereoIsotopic ) && ( Stereo_Prev = pINChI_Prev->StereoIsotopic ) &&
                        Eql_INChI_Stereo( Stereo, EQL_SP3, Stereo_Prev, EQL_SP3, bRelRac );
                    if (eq2prev)
                    {
                        mult++; /* mult = (number of non-empty equal items)-1 */
                        continue;
                    }
                    else
                    {
                        if (bNext++)
                        {
                            MakeDelim( sCompDelim, strbuf, bOverflow );
                        }
                        if (pINChI_Prev && pINChI_Prev->nNumberOfAtoms)
                        {
                            if (( Stereo_Prev = pINChI_Prev->StereoIsotopic ) && Stereo_Prev->nNumberOfStereoCenters > bRelRac)
                            {
                                MakeMult( mult + 1, "*", strbuf, 0, bOverflow );

                                MakeStereoString( Stereo_Prev->nNumber, NULL, Stereo_Prev->t_parity,
                                    0, Stereo_Prev->nNumberOfStereoCenters,
                                    strbuf, TAUT_MODE, bOverflow );
                            }
                            /* else sp3 info is not present in pINChI_Prev */
                        }
                        else
                            /* do not print pINChI_Prev because it either do not exist of have already been printed */
                            if (bSecondNonTautPass && pINChI_Taut_Prev && pINChI_Taut_Prev->nNumberOfAtoms)
                            {
                                if (( Stereo_Taut_Prev = pINChI_Taut_Prev->StereoIsotopic ) && Stereo_Taut_Prev->nNumberOfStereoCenters > bRelRac)
                                {
                                    /* since pINChI_Prev does not exist, pINChI_Taut_Prev is non-tautomeric */
                                    /* and it has non-trivial sp2 info */
                                    /*
                                    tot_len += MakeDelim( sIdenticalValues, strbuf, bOverflow);
                                    */
                                    ;/* pINChI_Taut_Prev sp3 info was output in the main stereo section */
                                }
                                else
                                {
                                    ; /* pINChI_Taut_Prev exists and has not sp3 info */
                                }
                            }
                    }
                    pINChI_Prev = pINChI;
                    pINChI_Taut_Prev = pINChI_Taut;
                    mult = 0; /* we do not know whether the item is empty */
                }
            }
        }

        return ( strbuf->nUsedLength - nUsedLength0 );
    }
        */
    // END INCHI C FUNCTION: str_IsoSp3

    fn clone_inchi(
        heap: &SourceHeap,
        pointer: SourceMutPointer<crate::source_types::INChI>,
    ) -> Result<Option<crate::source_types::INChI>, SourceHeapError> {
        if pointer.is_null() {
            Ok(None)
        } else {
            Ok(Some(
                heap.slice(pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            ))
        }
    }

    fn clone_stereo(
        heap: &SourceHeap,
        inchi: Option<&crate::source_types::INChI>,
        isotopic: bool,
    ) -> Result<Option<crate::source_types::INChI_Stereo>, SourceHeapError> {
        let Some(inchi) = inchi else {
            return Ok(None);
        };
        let pointer = if isotopic { inchi.StereoIsotopic } else { inchi.Stereo };
        if pointer.is_null() {
            Ok(None)
        } else {
            Ok(Some(
                heap.slice(pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone(),
            ))
        }
    }

    fn equivalence_pointer(heap: &mut SourceHeap, value: &str) -> Result<SourceConstPointer<i8>, SourceHeapError> {
        heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain(std::iter::once(0)).collect())
            .map(SourceMutPointer::as_const)
    }

    fn stereo_equal(
        heap: &SourceHeap,
        left: Option<&crate::source_types::INChI_Stereo>,
        right: Option<&crate::source_types::INChI_Stereo>,
        relative_or_racemic: i32,
    ) -> Result<bool, SourceHeapError> {
        Ok(left.is_some()
            && right.is_some()
            && Eql_INChI_Stereo(
                heap,
                left,
                crate::source_types::EQL_SP3 as i32,
                right,
                crate::source_types::EQL_SP3 as i32,
                relative_or_racemic,
            )? != 0)
    }

    fn output_stereo(
        heap: &mut SourceHeap,
        inchi: Option<&crate::source_types::INChI>,
        multiplier: i32,
        star: SourceConstPointer<i8>,
        string_buffer: &mut INCHI_IOS_STRING,
        taut_mode: i32,
        overflow: &mut i32,
        relative_or_racemic: i32,
    ) -> Result<(), SourceHeapError> {
        let Some(inchi) = inchi.filter(|inchi| inchi.nNumberOfAtoms != 0) else {
            return Ok(());
        };
        let Some(stereo) = clone_stereo(heap, Some(inchi), true)? else {
            return Ok(());
        };
        if stereo.nNumberOfStereoCenters > relative_or_racemic {
            MakeMult(heap, multiplier, star, string_buffer, 0, overflow)?;
            MakeStereoString(
                heap,
                stereo.nNumber.as_const(),
                SourceConstPointer::null(),
                stereo.t_parity.as_const(),
                0,
                stereo.nNumberOfStereoCenters,
                string_buffer,
                taut_mode,
                overflow,
            )?;
        }
        Ok(())
    }

    relative_or_racemic = 0;
    let initial_length = string_buffer.nUsedLength;
    let component_delimiter = heap.allocate_model_storage(vec![b';' as i8, 0])?;
    let star = heap.allocate_model_storage(vec![b'*' as i8, 0])?;
    let mut previous = SourceMutPointer::null();
    let mut taut_previous = SourceMutPointer::null();
    let mut multiplier = 0_i32;
    let mut next = 0_i32;
    let mut previous_equivalence: Option<&'static str> = None;
    let mut previous_equivalence_multiplier = 0_i32;
    let mut equal_to_taut_previous = 1_i32;
    let mut index = 0_i32;

    while index <= number_of_components {
        let current = if index < number_of_components {
            selected_inchi(heap, sorted_inchi.as_const().offset(i64::from(index))?, output_type)?
        } else {
            SourceMutPointer::null()
        };
        let taut_current = if second_non_taut_pass != 0 && index < number_of_components {
            selected_inchi(heap, sorted_inchi2.as_const().offset(i64::from(index))?, OUT_T1 as i32)?
        } else {
            SourceMutPointer::null()
        };
        let current_value = clone_inchi(heap, current)?;
        let taut_current_value = clone_inchi(heap, taut_current)?;
        let current_iso = clone_stereo(heap, current_value.as_ref(), true)?;
        let current_main = clone_stereo(heap, current_value.as_ref(), false)?;
        let taut_main = clone_stereo(heap, taut_current_value.as_ref(), false)?;
        let taut_iso = clone_stereo(heap, taut_current_value.as_ref(), true)?;
        let mut equivalence = 0_i32;
        if second_non_taut_pass != 0 && omit_repetitions != 0 {
            if stereo_equal(heap, current_iso.as_ref(), current_main.as_ref(), relative_or_racemic)? {
                equivalence = (crate::source_types::iiSTEREO
                    | crate::source_types::iitISO
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2NONTAUT) as i32;
            } else if stereo_equal(heap, current_iso.as_ref(), taut_main.as_ref(), relative_or_racemic)? {
                equivalence = (crate::source_types::iiSTEREO
                    | crate::source_types::iitISO
                    | crate::source_types::iitNONTAUT) as i32;
            } else if stereo_equal(heap, current_iso.as_ref(), taut_iso.as_ref(), relative_or_racemic)? {
                equivalence = (crate::source_types::iiSTEREO
                    | crate::source_types::iitISO
                    | crate::source_types::iitNONTAUT
                    | crate::source_types::iiEq2ISO) as i32;
            }
        } else if second_non_taut_pass == 0
            && omit_repetitions != 0
            && stereo_equal(heap, current_iso.as_ref(), current_main.as_ref(), relative_or_racemic)?
        {
            equivalence = (crate::source_types::iiSTEREO | crate::source_types::iitISO) as i32;
        }

        if equivalence != 0 {
            let previous_value = clone_inchi(heap, previous)?;
            if previous_value.as_ref().is_some_and(|inchi| inchi.nNumberOfAtoms != 0) {
                if next != 0 {
                    MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                }
                next = next.wrapping_add(1);
                output_stereo(
                    heap,
                    previous_value.as_ref(),
                    multiplier.wrapping_add(1),
                    star.as_const(),
                    string_buffer,
                    taut_mode,
                    overflow,
                    relative_or_racemic,
                )?;
            } else {
                let taut_previous_value = clone_inchi(heap, taut_previous)?;
                if taut_previous_value
                    .as_ref()
                    .is_some_and(|inchi| inchi.nNumberOfAtoms != 0)
                {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                }
            }
            let current_equivalence = crate::source::base::ichiprt1::EquString(equivalence);
            if previous_equivalence_multiplier != 0 && previous_equivalence.is_some() {
                if previous_equivalence == Some(current_equivalence) {
                    previous_equivalence_multiplier = previous_equivalence_multiplier.wrapping_add(1);
                } else {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let pointer = equivalence_pointer(heap, previous_equivalence.unwrap())?;
                    MakeEqStr(heap, pointer, previous_equivalence_multiplier, string_buffer, overflow)?;
                    previous_equivalence = Some(current_equivalence);
                    previous_equivalence_multiplier = 1;
                }
            } else {
                previous_equivalence = Some(current_equivalence);
                previous_equivalence_multiplier = 1;
            }
            previous = SourceMutPointer::null();
            taut_previous = SourceMutPointer::null();
            multiplier = 0;
            equal_to_taut_previous = 1;
        } else if equal_to_taut_previous != 0 {
            if previous_equivalence_multiplier != 0 {
                if let Some(value) = previous_equivalence {
                    if next != 0 {
                        MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
                    }
                    next = next.wrapping_add(1);
                    let pointer = equivalence_pointer(heap, value)?;
                    MakeEqStr(heap, pointer, previous_equivalence_multiplier, string_buffer, overflow)?;
                    previous_equivalence = None;
                    previous_equivalence_multiplier = 0;
                }
            }
            equal_to_taut_previous = 0;
            previous = current;
            taut_previous = taut_current;
            multiplier = 0;
        } else {
            let previous_value = clone_inchi(heap, previous)?;
            let previous_iso = clone_stereo(heap, previous_value.as_ref(), true)?;
            let equal_to_previous = use_multipliers != 0
                && stereo_equal(heap, current_iso.as_ref(), previous_iso.as_ref(), relative_or_racemic)?;
            if equal_to_previous {
                multiplier = multiplier.wrapping_add(1);
                index = index.wrapping_add(1);
                continue;
            }
            if next != 0 {
                MakeDelim(heap, component_delimiter.as_const(), string_buffer, overflow)?;
            }
            next = next.wrapping_add(1);
            output_stereo(
                heap,
                previous_value.as_ref(),
                multiplier.wrapping_add(1),
                star.as_const(),
                string_buffer,
                taut_mode,
                overflow,
                relative_or_racemic,
            )?;
            previous = current;
            taut_previous = taut_current;
            multiplier = 0;
        }
        index = index.wrapping_add(1);
    }
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn bin_AuxTautTrans(
    heap: &mut SourceHeap,
    sorted_inchi: SourceMutPointer<INCHI_SORT>,
    sorted_inchi2: SourceMutPointer<INCHI_SORT>,
    transposed_non_taut: &mut SourceMutPointer<u16>,
    transposed_taut: &mut SourceMutPointer<u16>,
    output_type: i32,
    number_of_components: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:4091 bin_AuxTautTrans
    // INCHI✔️❌: int bin_AuxTautTrans( INCHI_SORT *pINChISort,
    // INCHI✔️❌:     INCHI_SORT *pINChISort2,
    // INCHI✔️❌:     AT_NUMB    **pTrans_n,
    // INCHI✔️❌:     AT_NUMB    **pTrans_s,
    // INCHI✔️❌:     int        bOutType,
    // INCHI✔️❌:     int        num_components )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int          i, ii, ii2, ret = 0;
    // INCHI✔️❌:     INCHI_SORT   *is, *is2, *is0, *is20;
    // INCHI✔️❌:     INChI        *pINChI, *pINChI_Taut;
    // INCHI✔️❌:     AT_NUMB     *nTrans_n = NULL;
    // INCHI✔️❌:     AT_NUMB     *nTrans_s = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     is0 = pINChISort;
    // INCHI✔️❌:     is20 = pINChISort2;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* djb-rwth: rewritten to avoid memory leaks -- updated 25/09/2025 */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Pass 1: save new non-taut numbering */
    // INCHI✔️❌:     /* For each connected component...    */
    // INCHI✔️❌:     for (i = 0; i < num_components; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         is = is0 + i;
    // INCHI✔️❌:         is2 = is20 + i;
    // INCHI✔️❌:         pINChI = (0 <= (ii = GET_II(bOutType, is))) ? is->pINChI[ii] : NULL;
    // INCHI✔️❌:         pINChI_Taut = (0 <= (ii2 = GET_II(OUT_T1, is2))) ? is2->pINChI[ii2] : NULL;
    // INCHI✔️❌:         if (pINChI && pINChI->nNumberOfAtoms > 0 &&
    // INCHI✔️❌:             pINChI_Taut && pINChI_Taut->nNumberOfAtoms > 0 &&
    // INCHI✔️❌:             /* different components save equal new ord. numbers: */
    // INCHI✔️❌:             is->ord_number != is2->ord_number)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!nTrans_n)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nTrans_n = (AT_NUMB*)inchi_calloc((long long)num_components + 1, sizeof(nTrans_n[0]));
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (!nTrans_s)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nTrans_s = (AT_NUMB*)inchi_calloc((long long)num_components + 1, sizeof(nTrans_s[0]));
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (nTrans_n && nTrans_s)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* new ordering number for original non-tautomeric component number is->ord_number */
    // INCHI✔️❌:                 nTrans_n[is->ord_number] = i + 1; /*nTrans_t[is2->ord_number] =*/
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (nTrans_n && nTrans_s)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Pass 2: get new taut numbering, retrieve new non-taut and save the transposition */
    // INCHI✔️❌:         for (i = 0; i < num_components; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             is = is0 + i;
    // INCHI✔️❌:             is2 = is20 + i;
    // INCHI✔️❌:             pINChI = (0 <= (ii = GET_II(bOutType, is))) ? is->pINChI[ii] : NULL;
    // INCHI✔️❌:             pINChI_Taut = (0 <= (ii2 = GET_II(OUT_T1, is2))) ? is2->pINChI[ii2] : NULL;
    // INCHI✔️❌:             if (pINChI && pINChI->nNumberOfAtoms > 0 &&
    // INCHI✔️❌:                 pINChI_Taut && pINChI_Taut->nNumberOfAtoms > 0 &&
    // INCHI✔️❌:                 is->ord_number != is2->ord_number &&
    // INCHI✔️❌:                 nTrans_n[is2->ord_number])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* nTrans_n[is2->ord_number] is new ordering number of
    // INCHI✔️❌:                 the non-taut representation of the tautomeric component
    // INCHI✔️❌:                 that has new ord number i+1 and orig ordering number is2->ord_number.
    // INCHI✔️❌:                 Old numbers start from 0, new start from 1
    // INCHI✔️❌:                 */
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* n = nTrans_s[t]: taut component #t is in position #n of the non-taut representation */
    // INCHI✔️❌:                 nTrans_s[i + 1] = nTrans_n[is2->ord_number];
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         *pTrans_n = nTrans_n;
    // INCHI✔️❌:         *pTrans_s = nTrans_s;
    // INCHI✔️❌:         ret = 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (nTrans_n)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free(nTrans_n);
    // INCHI✔️❌:             ret = -1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (nTrans_s)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free(nTrans_s);
    // INCHI✔️❌:             ret = -1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: bin_AuxTautTrans

    fn has_atoms(
        heap: &SourceHeap,
        pointer: SourceMutPointer<crate::source_types::INChI>,
    ) -> Result<bool, SourceHeapError> {
        if pointer.is_null() {
            return Ok(false);
        }
        Ok(heap
            .slice(pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .nNumberOfAtoms
            > 0)
    }

    let mut non_taut_transposition = SourceMutPointer::null();
    let mut taut_transposition = SourceMutPointer::null();
    let allocation_count = if number_of_components >= 0 {
        u64::try_from(i64::from(number_of_components) + 1)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?
    } else {
        0
    };
    let mut index = 0_i32;
    while index < number_of_components {
        let sort = heap
            .slice(sorted_inchi.as_const().offset(i64::from(index))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let sort2 = heap
            .slice(sorted_inchi2.as_const().offset(i64::from(index))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let inchi =
            get_ii(heap, output_type, &sort)?.map_or(SourceMutPointer::null(), |selected| sort.pINChI[selected]);
        let taut =
            get_ii(heap, OUT_T1 as i32, &sort2)?.map_or(SourceMutPointer::null(), |selected| sort2.pINChI[selected]);
        if has_atoms(heap, inchi)? && has_atoms(heap, taut)? && sort.ord_number != sort2.ord_number {
            if non_taut_transposition.is_null() {
                non_taut_transposition = match inchi_calloc(
                    heap,
                    allocation_count,
                    u64::try_from(std::mem::size_of::<u16>()).map_err(|_| SourceHeapError::AllocationSizeOverflow)?,
                ) {
                    Ok(pointer) => pointer,
                    Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                    Err(error) => return Err(error),
                };
            }
            if taut_transposition.is_null() {
                taut_transposition = match inchi_calloc(
                    heap,
                    allocation_count,
                    u64::try_from(std::mem::size_of::<u16>()).map_err(|_| SourceHeapError::AllocationSizeOverflow)?,
                ) {
                    Ok(pointer) => pointer,
                    Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                    Err(error) => return Err(error),
                };
            }
            if !non_taut_transposition.is_null() && !taut_transposition.is_null() {
                let order = usize::try_from(sort.ord_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                *heap
                    .slice_mut(non_taut_transposition)?
                    .get_mut(order)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = index.wrapping_add(1) as u16;
            }
        }
        index = index.wrapping_add(1);
    }

    if !non_taut_transposition.is_null() && !taut_transposition.is_null() {
        index = 0;
        while index < number_of_components {
            let sort = heap
                .slice(sorted_inchi.as_const().offset(i64::from(index))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let sort2 = heap
                .slice(sorted_inchi2.as_const().offset(i64::from(index))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let inchi =
                get_ii(heap, output_type, &sort)?.map_or(SourceMutPointer::null(), |selected| sort.pINChI[selected]);
            let taut = get_ii(heap, OUT_T1 as i32, &sort2)?
                .map_or(SourceMutPointer::null(), |selected| sort2.pINChI[selected]);
            let source_order = usize::try_from(sort2.ord_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let transposed = *heap
                .slice(non_taut_transposition.as_const())?
                .get(source_order)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if has_atoms(heap, inchi)?
                && has_atoms(heap, taut)?
                && sort.ord_number != sort2.ord_number
                && transposed != 0
            {
                let destination =
                    usize::try_from(index.wrapping_add(1)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                *heap
                    .slice_mut(taut_transposition)?
                    .get_mut(destination)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = transposed;
            }
            index = index.wrapping_add(1);
        }
        *transposed_non_taut = non_taut_transposition;
        *transposed_taut = taut_transposition;
        return Ok(1);
    }

    let mut result = 0;
    if !non_taut_transposition.is_null() {
        inchi_free(heap, non_taut_transposition)?;
        result = -1;
    }
    if !taut_transposition.is_null() {
        inchi_free(heap, taut_transposition)?;
        result = -1;
    }
    Ok(result)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn str_AuxTautTrans(
    heap: &mut SourceHeap,
    transposed_non_taut: SourceMutPointer<u16>,
    transposed_taut: SourceMutPointer<u16>,
    string_buffer: &mut INCHI_IOS_STRING,
    overflow: &mut i32,
    taut_mode: i32,
    number_of_components: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:4187 str_AuxTautTrans
    // INCHI✔️❌: int str_AuxTautTrans( CANON_GLOBALS     *pCG,
    // INCHI✔️❌:     AT_NUMB           *nTrans_n,
    // INCHI✔️❌:     AT_NUMB           *nTrans_s,
    // INCHI✔️❌:     INCHI_IOS_STRING  *strbuf,
    // INCHI✔️❌:     int               *bOverflow,
    // INCHI✔️❌:     int               TAUT_MODE,
    // INCHI✔️❌:     int               num_components )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, k, len, j, nUsedLength0;
    // INCHI✔️❌:
    // INCHI✔️❌:     nUsedLength0 = strbuf->nUsedLength;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nTrans_n && nTrans_s)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* print the transposition, cycle after cycle */
    // INCHI✔️❌:         for (i = 1; i <= num_components; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (nTrans_s[i])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* get one cycle of the transposition */
    // INCHI✔️❌:                 for (j = i, len = 0; ( k = nTrans_s[j] ); j = k, len++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nTrans_n[len] = j; /* save the transposition */
    // INCHI✔️❌:                     nTrans_s[j] = 0; /* clear used element to avoid repetitions */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* print one cycle of the transposition */
    // INCHI✔️❌:                 MakeDelim( "(", strbuf, bOverflow );
    // INCHI✔️❌:                 MakeCtString( pCG, nTrans_n, len, 0, NULL, 0,
    // INCHI✔️❌:                     strbuf, TAUT_MODE, bOverflow );
    // INCHI✔️❌:                 MakeDelim( ")", strbuf, bOverflow );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (nTrans_n)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( nTrans_n );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (nTrans_s)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( nTrans_s );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ( strbuf->nUsedLength - nUsedLength0 );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: str_AuxTautTrans

    let initial_length = string_buffer.nUsedLength;
    if !transposed_non_taut.is_null() && !transposed_taut.is_null() {
        let opening = heap.allocate_model_storage(vec![b'(' as i8, 0])?;
        let closing = heap.allocate_model_storage(vec![b')' as i8, 0])?;
        let mut index = 1_i32;
        while index <= number_of_components {
            let index_usize = usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let first = *heap
                .slice(transposed_taut.as_const())?
                .get(index_usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if first != 0 {
                let mut current = index;
                let mut length = 0_i32;
                loop {
                    let current_usize = usize::try_from(current).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let next = *heap
                        .slice(transposed_taut.as_const())?
                        .get(current_usize)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if next == 0 {
                        break;
                    }
                    let length_usize = usize::try_from(length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    *heap
                        .slice_mut(transposed_non_taut)?
                        .get_mut(length_usize)
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = current as u16;
                    *heap
                        .slice_mut(transposed_taut)?
                        .get_mut(current_usize)
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
                    current = i32::from(next);
                    length = length.wrapping_add(1);
                }
                MakeDelim(heap, opening.as_const(), string_buffer, overflow)?;
                MakeCtString(
                    heap,
                    transposed_non_taut,
                    length,
                    0,
                    SourceConstPointer::null(),
                    0,
                    string_buffer,
                    taut_mode,
                    overflow,
                )?;
                MakeDelim(heap, closing.as_const(), string_buffer, overflow)?;
            }
            index = index.wrapping_add(1);
        }
        inchi_free(heap, opening)?;
        inchi_free(heap, closing)?;
    }
    inchi_free(heap, transposed_non_taut)?;
    inchi_free(heap, transposed_taut)?;
    Ok(string_buffer.nUsedLength.wrapping_sub(initial_length))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{INChI, INChI_Aux, INChI_Stereo, OUT_N1, OUT_T1};

    fn text(heap: &mut SourceHeap, value: &str) -> SourceMutPointer<i8> {
        heap.allocate_model_storage(value.bytes().map(|byte| byte as i8).chain(std::iter::once(0)).collect())
            .unwrap()
    }

    fn inchi(heap: &mut SourceHeap, formula: &str, atoms: i32, tautomer: i32) -> SourceMutPointer<INChI> {
        let formula = text(heap, formula);
        heap.allocate_model_storage(vec![INChI {
            nNumberOfAtoms: atoms,
            lenTautomer: tautomer,
            szHillFormula: formula,
            ..INChI::default()
        }])
        .unwrap()
    }

    fn hydrogen_inchi(heap: &mut SourceHeap, hydrogens: &[i8], tautomer: &[u16]) -> SourceMutPointer<INChI> {
        let hydrogen_pointer = if hydrogens.is_empty() {
            SourceMutPointer::null()
        } else {
            heap.allocate_model_storage(hydrogens.to_vec()).unwrap()
        };
        let tautomer_pointer = if tautomer.is_empty() {
            SourceMutPointer::null()
        } else {
            heap.allocate_model_storage(tautomer.to_vec()).unwrap()
        };
        heap.allocate_model_storage(vec![INChI {
            nNumberOfAtoms: i32::try_from(hydrogens.len()).unwrap(),
            nNum_H: hydrogen_pointer,
            lenTautomer: i32::try_from(tautomer.len()).unwrap(),
            nTautomer: tautomer_pointer,
            ..INChI::default()
        }])
        .unwrap()
    }

    fn output_buffer(heap: &mut SourceHeap, prefix: &str) -> INCHI_IOS_STRING {
        let mut bytes = vec![81_i8; 128];
        for (target, source) in bytes.iter_mut().zip(prefix.bytes()) {
            *target = source as i8;
        }
        bytes[prefix.len()] = 0;
        INCHI_IOS_STRING {
            pStr: heap.allocate_model_storage(bytes).unwrap(),
            nAllocatedLength: 128,
            nUsedLength: prefix.len() as i32,
            nPtr: 32,
        }
    }

    fn visible(heap: &SourceHeap, buffer: &INCHI_IOS_STRING) -> String {
        let bytes = heap.slice(buffer.pStr.as_const()).unwrap();
        let nul = bytes.iter().position(|byte| *byte == 0).unwrap();
        String::from_utf8(bytes[..nul].iter().map(|byte| *byte as u8).collect()).unwrap()
    }

    #[test]
    fn source_port__ichiprt3__str_auxequ__line_2106() {
        let mut heap = SourceHeap::default();
        let mut output = output_buffer(&mut heap, "pre:");
        let mut overflow = 0;
        assert_eq!(
            str_AuxEqu(
                &mut heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                0,
                0,
                0,
                0,
            ),
            Ok(0)
        );
        assert_eq!(visible(&heap, &output), "pre:");

        let numbers = heap.allocate_model_storage(vec![1_u16, 1, 3]).unwrap();
        let aux = heap
            .allocate_model_storage(vec![crate::source_types::INChI_Aux {
                nNumberOfAtoms: 3,
                nConstitEquNumbers: numbers,
                ..crate::source_types::INChI_Aux::default()
            }])
            .unwrap();
        let inchi = heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 3,
                ..INChI::default()
            }])
            .unwrap();
        let sort = heap
            .allocate_model_storage(vec![crate::source_types::INCHI_SORT {
                pINChI: [inchi, SourceMutPointer::null()],
                pINChI_Aux: [aux, SourceMutPointer::null()],
                ..crate::source_types::INCHI_SORT::default()
            }])
            .unwrap();
        output = output_buffer(&mut heap, "");
        overflow = 0;
        assert_eq!(
            str_AuxEqu(
                &mut heap,
                sort,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                1,
                0,
                0,
                0,
            ),
            Ok(8)
        );
        assert_eq!(visible(&heap, &output), "(1,2)(3)");
        assert_eq!(overflow, 0);
    }

    #[test]
    fn source_port__ichiprt3__str_auxisoequ__line_2962() {
        use crate::source_types::{INChI_Aux, OUT_NT, OUT_T1, TAUT_NON, TAUT_YES};

        fn aux(heap: &mut SourceHeap, normal: &[u16], isotopic: Option<&[u16]>) -> SourceMutPointer<INChI_Aux> {
            let normal_pointer = heap.allocate_model_storage(normal.to_vec()).unwrap();
            let isotopic_pointer = isotopic.map_or(SourceMutPointer::null(), |values| {
                heap.allocate_model_storage(values.to_vec()).unwrap()
            });
            heap.allocate_model_storage(vec![INChI_Aux {
                nNumberOfAtoms: i32::try_from(normal.len()).unwrap(),
                bIsIsotopic: i32::from(isotopic.is_some()),
                nConstitEquNumbers: normal_pointer,
                nConstitEquIsotopicNumbers: isotopic_pointer,
                ..INChI_Aux::default()
            }])
            .unwrap()
        }

        fn taut_sort(heap: &mut SourceHeap, normal: &[u16], isotopic: Option<&[u16]>) -> INCHI_SORT {
            let aux = aux(heap, normal, isotopic);
            let inchi = heap
                .allocate_model_storage(vec![INChI {
                    nNumberOfAtoms: i32::try_from(normal.len()).unwrap(),
                    ..INChI::default()
                }])
                .unwrap();
            let mut sort = INCHI_SORT::default();
            sort.pINChI[TAUT_YES as usize] = inchi;
            sort.pINChI_Aux[TAUT_YES as usize] = aux;
            sort
        }

        fn non_taut_sort(heap: &mut SourceHeap, normal: &[u16], isotopic: Option<&[u16]>) -> INCHI_SORT {
            let taut_aux = aux(heap, &[1, 2, 3], None);
            let taut_inchi = heap
                .allocate_model_storage(vec![INChI {
                    nNumberOfAtoms: 3,
                    lenTautomer: 1,
                    ..INChI::default()
                }])
                .unwrap();
            let non_aux = aux(heap, normal, isotopic);
            let non_inchi = heap
                .allocate_model_storage(vec![INChI {
                    nNumberOfAtoms: i32::try_from(normal.len()).unwrap(),
                    ..INChI::default()
                }])
                .unwrap();
            let mut sort = INCHI_SORT::default();
            sort.pINChI[TAUT_YES as usize] = taut_inchi;
            sort.pINChI_Aux[TAUT_YES as usize] = taut_aux;
            sort.pINChI[TAUT_NON as usize] = non_inchi;
            sort.pINChI_Aux[TAUT_NON as usize] = non_aux;
            sort
        }

        fn run(
            heap: &mut SourceHeap,
            current: Vec<INCHI_SORT>,
            taut: Vec<INCHI_SORT>,
            output_type: i32,
            second_pass: i32,
            omit_repetitions: i32,
            use_multipliers: i32,
        ) -> (i32, String, i32) {
            let count = i32::try_from(current.len()).unwrap();
            let current = if current.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(current).unwrap()
            };
            let taut = if taut.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(taut).unwrap()
            };
            let mut output = output_buffer(heap, "pre:");
            let mut overflow = 0;
            let length = str_AuxIsoEqu(
                heap,
                current,
                taut,
                &mut output,
                &mut overflow,
                output_type,
                0,
                count,
                second_pass,
                omit_repetitions,
                use_multipliers,
            )
            .unwrap();
            (length, visible(heap, &output), overflow)
        }

        let mut heap = SourceHeap::default();
        assert_eq!(
            run(&mut heap, vec![], vec![], OUT_T1 as i32, 0, 0, 0),
            (0, "pre:".into(), 0)
        );

        let concrete = taut_sort(&mut heap, &[1, 2, 3], Some(&[1, 1, 3]));
        assert_eq!(
            run(&mut heap, vec![concrete], vec![], OUT_T1 as i32, 0, 0, 0,),
            (8, "pre:(1,2)(3)".into(), 0)
        );

        let first_equal = taut_sort(&mut heap, &[1, 1, 3], Some(&[1, 1, 3]));
        assert_eq!(
            run(&mut heap, vec![first_equal], vec![], OUT_T1 as i32, 0, 1, 0,),
            (1, "pre:m".into(), 0)
        );

        let current = non_taut_sort(&mut heap, &[1, 2, 2], Some(&[1, 1, 3]));
        let taut = taut_sort(&mut heap, &[1, 1, 3], Some(&[1, 2, 2]));
        assert_eq!(
            run(&mut heap, vec![current], vec![taut], OUT_NT as i32, 1, 1, 0,),
            (1, "pre:m".into(), 0)
        );

        let current = non_taut_sort(&mut heap, &[1, 1, 3], Some(&[1, 1, 3]));
        let taut = taut_sort(&mut heap, &[1, 2, 2], Some(&[1, 3, 3]));
        assert_eq!(
            run(&mut heap, vec![current], vec![taut], OUT_NT as i32, 1, 1, 0,),
            (1, "pre:n".into(), 0)
        );

        let current = non_taut_sort(&mut heap, &[1, 2, 2], Some(&[1, 1, 3]));
        let taut = taut_sort(&mut heap, &[1, 3, 3], Some(&[1, 1, 3]));
        assert_eq!(
            run(&mut heap, vec![current], vec![taut], OUT_NT as i32, 1, 1, 0,),
            (1, "pre:M".into(), 0)
        );

        let repeated = taut_sort(&mut heap, &[1, 2, 3], Some(&[1, 1, 3]));
        assert_eq!(
            run(
                &mut heap,
                vec![repeated.clone(), repeated],
                vec![],
                OUT_T1 as i32,
                0,
                0,
                1,
            ),
            (10, "pre:2*(1,2)(3)".into(), 0)
        );
    }

    #[test]
    fn source_port__ichiprt3__str_auxtgroupequ__line_3953() {
        let mut heap = SourceHeap::default();
        let numbers = heap.allocate_model_storage(vec![1_u16, 1]).unwrap();
        let aux = heap
            .allocate_model_storage(vec![crate::source_types::INChI_Aux {
                nNumberOfAtoms: 2,
                nNumberOfTGroups: 2,
                nConstitEquTGroupNumbers: numbers,
                ..crate::source_types::INChI_Aux::default()
            }])
            .unwrap();
        let inchi = heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 2,
                ..INChI::default()
            }])
            .unwrap();
        let sort = heap
            .allocate_model_storage(vec![crate::source_types::INCHI_SORT {
                pINChI: [inchi, SourceMutPointer::null()],
                pINChI_Aux: [aux, SourceMutPointer::null()],
                ..crate::source_types::INCHI_SORT::default()
            }])
            .unwrap();
        let mut output = output_buffer(&mut heap, "");
        let mut overflow = 0;
        assert_eq!(
            str_AuxTgroupEqu(&mut heap, sort, &mut output, &mut overflow, OUT_T1 as i32, 0, 1, 0,),
            Ok(5)
        );
        assert_eq!(visible(&heap, &output), "(1,2)");
        assert_eq!(overflow, 0);
    }

    #[test]
    fn source_port__ichiprt3__str_auxchargeradval__line_4013() {
        fn component(
            heap: &mut SourceHeap,
            info: Vec<crate::source_types::ORIG_INFO>,
            tautomeric_slot: bool,
        ) -> INCHI_SORT {
            let count = i32::try_from(info.len()).unwrap();
            let original_info = if info.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(info).unwrap()
            };
            let aux = heap
                .allocate_model_storage(vec![INChI_Aux {
                    nNumberOfAtoms: count,
                    OrigInfo: original_info,
                    ..INChI_Aux::default()
                }])
                .unwrap();
            let inchi = heap
                .allocate_model_storage(vec![INChI {
                    nNumberOfAtoms: count,
                    lenTautomer: i32::from(tautomeric_slot),
                    ..INChI::default()
                }])
                .unwrap();
            let slot = usize::from(tautomeric_slot);
            let mut sort = INCHI_SORT::default();
            sort.pINChI[slot] = inchi;
            sort.pINChI_Aux[slot] = aux;
            sort
        }

        fn run(
            heap: &mut SourceHeap,
            components: Vec<INCHI_SORT>,
            output_type: i32,
            mode: i32,
            use_multipliers: i32,
            initial_overflow: i32,
        ) -> (Result<i32, SourceHeapError>, String, i32) {
            let count = i32::try_from(components.len()).unwrap();
            let components = if components.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(components).unwrap()
            };
            let mut output = output_buffer(heap, "pre:");
            let mut overflow = initial_overflow;
            let result = str_AuxChargeRadVal(
                heap,
                components,
                &mut output,
                &mut overflow,
                output_type,
                mode,
                count,
                use_multipliers,
            );
            (result, visible(heap, &output), overflow)
        }

        let info = |charge, radical, valence| crate::source_types::ORIG_INFO {
            cCharge: charge,
            cRadical: radical,
            cUnusualValence: valence,
        };
        let mut heap = SourceHeap::default();
        assert_eq!(
            run(&mut heap, vec![], OUT_T1 as i32, 0, 0, 0),
            (Ok(0), "pre:".into(), 0)
        );

        let charged = component(&mut heap, vec![info(1, 0, 0)], true);
        assert_eq!(
            run(&mut heap, vec![charged], OUT_T1 as i32, 0, 0, 0),
            (Ok(3), "pre:1+1".into(), 0)
        );

        let sparse = component(
            &mut heap,
            vec![crate::source_types::ORIG_INFO::default(), info(0, 1, 0)],
            true,
        );
        assert_eq!(
            run(
                &mut heap,
                vec![sparse],
                OUT_T1 as i32,
                crate::source_types::CT_MODE_ABC_NUMBERS as i32,
                0,
                0,
            ),
            (Ok(3), "pre:B.d".into(), 0)
        );

        let first = component(&mut heap, vec![info(1, 0, 0)], true);
        let repeated = first.clone();
        assert_eq!(
            run(&mut heap, vec![first, repeated], OUT_T1 as i32, 0, 1, 0,),
            (Ok(5), "pre:2*1+1".into(), 0)
        );

        let charged = component(&mut heap, vec![info(1, 0, 0)], true);
        let radical = component(&mut heap, vec![info(0, 1, 0)], true);
        assert_eq!(
            run(&mut heap, vec![charged, radical], OUT_T1 as i32, 0, 1, 0,),
            (Ok(6), "pre:1+1;1d".into(), 0)
        );

        let trivial = component(&mut heap, vec![crate::source_types::ORIG_INFO::default()], true);
        let charged = component(&mut heap, vec![info(1, 0, 0)], true);
        assert_eq!(
            run(&mut heap, vec![trivial, charged], OUT_T1 as i32, 0, 0, 0,),
            (Ok(4), "pre:;1+1".into(), 0)
        );

        let non_taut = component(&mut heap, vec![info(-2, 0, 3)], false);
        assert_eq!(
            run(&mut heap, vec![non_taut], crate::source_types::OUT_NN as i32, 0, 0, 0,),
            (Ok(5), "pre:1-2.3".into(), 0)
        );

        let selected = component(&mut heap, vec![info(1, 0, 0)], true);
        assert_eq!(run(&mut heap, vec![selected], -1, 0, 0, 0), (Ok(0), "pre:".into(), 0));
        let held = component(&mut heap, vec![info(1, 0, 0)], true);
        assert_eq!(
            run(&mut heap, vec![held], OUT_T1 as i32, 0, 0, 7),
            (Ok(0), "pre:".into(), 7)
        );

        let mut output = output_buffer(&mut heap, "held");
        let mut overflow = 0;
        assert_eq!(
            str_AuxChargeRadVal(
                &mut heap,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                -1,
                0,
            ),
            Ok(0)
        );
        assert!(matches!(
            str_AuxChargeRadVal(
                &mut heap,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                1,
                0,
            ),
            Err(SourceHeapError::NullPointer | SourceHeapError::PointerOutOfBounds)
        ));
    }

    #[test]
    fn source_port__ichiprt3__str_auxisotgroupequ__line_4236() {
        fn component(heap: &mut SourceHeap, normal: &[u16], isotopic: &[u16]) -> INCHI_SORT {
            let normal = heap.allocate_model_storage(normal.to_vec()).unwrap();
            let isotopic = heap.allocate_model_storage(isotopic.to_vec()).unwrap();
            let aux = heap
                .allocate_model_storage(vec![INChI_Aux {
                    nNumberOfAtoms: 2,
                    nNumberOfTGroups: 2,
                    bIsIsotopic: 1,
                    nConstitEquTGroupNumbers: normal,
                    nConstitEquIsotopicTGroupNumbers: isotopic,
                    ..INChI_Aux::default()
                }])
                .unwrap();
            let inchi = heap
                .allocate_model_storage(vec![INChI {
                    nNumberOfAtoms: 2,
                    ..INChI::default()
                }])
                .unwrap();
            INCHI_SORT {
                pINChI: [inchi, SourceMutPointer::null()],
                pINChI_Aux: [aux, SourceMutPointer::null()],
                ..INCHI_SORT::default()
            }
        }

        fn run(
            heap: &mut SourceHeap,
            components: Vec<INCHI_SORT>,
            omit_repetitions: i32,
            use_multipliers: i32,
        ) -> (i32, String, i32) {
            let count = i32::try_from(components.len()).unwrap();
            let components = if components.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(components).unwrap()
            };
            let mut output = output_buffer(heap, "pre:");
            let mut overflow = 0;
            let length = str_AuxIsoTgroupEqu(
                heap,
                components,
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                count,
                omit_repetitions,
                use_multipliers,
            )
            .unwrap();
            (length, visible(heap, &output), overflow)
        }

        let mut heap = SourceHeap::default();
        assert_eq!(run(&mut heap, vec![], 0, 0), (0, "pre:".into(), 0));

        let concrete = component(&mut heap, &[1, 2], &[1, 1]);
        assert_eq!(run(&mut heap, vec![concrete], 0, 0), (5, "pre:(1,2)".into(), 0));

        let equal = component(&mut heap, &[1, 1], &[1, 1]);
        assert_eq!(run(&mut heap, vec![equal], 1, 0), (1, "pre:m".into(), 0));

        let equal_a = component(&mut heap, &[1, 1], &[1, 1]);
        let equal_b = equal_a.clone();
        assert_eq!(run(&mut heap, vec![equal_a, equal_b], 1, 0), (2, "pre:2m".into(), 0));

        let concrete_a = component(&mut heap, &[1, 2], &[1, 1]);
        let concrete_b = concrete_a.clone();
        assert_eq!(
            run(&mut heap, vec![concrete_a, concrete_b], 0, 1),
            (7, "pre:2*(1,2)".into(), 0)
        );

        let pending = component(&mut heap, &[1, 1], &[1, 1]);
        let following = component(&mut heap, &[1, 2], &[1, 1]);
        assert_eq!(
            run(&mut heap, vec![pending, following], 1, 0),
            (7, "pre:m;(1,2)".into(), 0)
        );

        let trivial = component(&mut heap, &[1, 2], &[1, 2]);
        let following = component(&mut heap, &[1, 2], &[1, 1]);
        assert_eq!(
            run(&mut heap, vec![trivial, following], 0, 0),
            (6, "pre:;(1,2)".into(), 0)
        );
    }

    #[test]
    fn source_port__ichiprt3__str_hillformula__line_60() {
        let mut heap = SourceHeap::default();
        let mut output = output_buffer(&mut heap, "pre:");
        let mut overflow = 0;
        assert_eq!(
            str_HillFormula(
                &mut heap,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                3,
                1,
            ),
            Ok(4)
        );
        assert_eq!(visible(&heap, &output), "pre:");

        let c2h6_a = inchi(&mut heap, "C2H6", 2, 1);
        let c2h6_b = inchi(&mut heap, "C2H6", 2, 1);
        let h2o = inchi(&mut heap, "H2O", 1, 1);
        let sorts = heap
            .allocate_model_storage(vec![
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), c2h6_a],
                    ..INCHI_SORT::default()
                },
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), c2h6_b],
                    ..INCHI_SORT::default()
                },
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), h2o],
                    ..INCHI_SORT::default()
                },
            ])
            .unwrap();
        output = output_buffer(&mut heap, "");
        assert_eq!(
            str_HillFormula(&mut heap, sorts, &mut output, &mut overflow, OUT_T1 as i32, 3, 1),
            Ok(9)
        );
        assert_eq!(visible(&heap, &output), "2C2H6.H2O");

        output = output_buffer(&mut heap, "");
        assert_eq!(
            str_HillFormula(&mut heap, sorts, &mut output, &mut overflow, OUT_T1 as i32, 3, 0),
            Ok(13)
        );
        assert_eq!(visible(&heap, &output), "C2H6.C2H6.H2O");

        let empty = inchi(&mut heap, "", 1, 0);
        let methane_non = inchi(&mut heap, "CH4", 1, 0);
        let methane_taut_non = inchi(&mut heap, "wrong", 1, 1);
        let selection = heap
            .allocate_model_storage(vec![
                INCHI_SORT {
                    pINChI: [empty, SourceMutPointer::null()],
                    ..INCHI_SORT::default()
                },
                INCHI_SORT {
                    pINChI: [methane_non, methane_taut_non],
                    ..INCHI_SORT::default()
                },
            ])
            .unwrap();
        output = output_buffer(&mut heap, "");
        assert_eq!(
            str_HillFormula(&mut heap, selection, &mut output, &mut overflow, OUT_N1 as i32, 2, 1),
            Ok(4)
        );
        assert_eq!(visible(&heap, &output), ".CH4");

        output = output_buffer(&mut heap, "");
        overflow = 1;
        assert_eq!(
            str_HillFormula(&mut heap, sorts, &mut output, &mut overflow, OUT_T1 as i32, 3, 1),
            Ok(0)
        );
        assert_eq!(visible(&heap, &output), "");
    }

    #[test]
    fn source_port__ichiprt3__str_connections__line_224() {
        fn connected(heap: &mut SourceHeap, table: &[u16], atoms: i32) -> SourceMutPointer<INChI> {
            let table = heap.allocate_model_storage(table.to_vec()).unwrap();
            heap.allocate_model_storage(vec![INChI {
                nNumberOfAtoms: atoms,
                lenConnTable: heap.slice(table.as_const()).unwrap().len() as i32,
                nConnTable: table,
                ..INChI::default()
            }])
            .unwrap()
        }

        let mut heap = SourceHeap::default();
        let mut output = output_buffer(&mut heap, "pre");
        let mut overflow = 0;
        assert_eq!(
            str_Connections(
                &mut heap,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                2,
                1,
            ),
            Ok(3)
        );
        assert_eq!(visible(&heap, &output), "pre");

        let first = connected(&mut heap, &[2, 1, 3, 2], 3);
        let second = connected(&mut heap, &[2, 1, 3, 2], 3);
        let sorts = heap
            .allocate_model_storage(vec![
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), first],
                    ..INCHI_SORT::default()
                },
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), second],
                    ..INCHI_SORT::default()
                },
            ])
            .unwrap();
        output = output_buffer(&mut heap, "");
        assert_eq!(
            str_Connections(&mut heap, sorts, &mut output, &mut overflow, OUT_T1 as i32, 0, 2, 1,),
            Ok(7)
        );
        assert_eq!(visible(&heap, &output), "2*1-2-3");

        output = output_buffer(&mut heap, "");
        assert_eq!(
            str_Connections(&mut heap, sorts, &mut output, &mut overflow, OUT_T1 as i32, 0, 2, 0,),
            Ok(11)
        );
        assert_eq!(visible(&heap, &output), "1-2-3;1-2-3");

        let empty1 = connected(&mut heap, &[1], 1);
        let empty2 = connected(&mut heap, &[1], 1);
        let empty_sorts = heap
            .allocate_model_storage(vec![
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), empty1],
                    ..INCHI_SORT::default()
                },
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), empty2],
                    ..INCHI_SORT::default()
                },
            ])
            .unwrap();
        output = output_buffer(&mut heap, "x");
        assert_eq!(
            str_Connections(
                &mut heap,
                empty_sorts,
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                2,
                0,
            ),
            Ok(0)
        );
        assert_eq!(visible(&heap, &output), "x");
    }

    #[test]
    fn source_port__ichiprt3__str_h_atoms__line_309() {
        let mut heap = SourceHeap::default();
        let h = hydrogen_inchi(&mut heap, &[1, 1, 0], &[]);
        let taut = hydrogen_inchi(&mut heap, &[0], &[1, 3, 1, 0, 5]);
        let both = hydrogen_inchi(&mut heap, &[1], &[1, 3, 1, 0, 5]);
        let empty1 = hydrogen_inchi(&mut heap, &[0], &[]);
        let empty2 = hydrogen_inchi(&mut heap, &[0], &[]);

        for (inchi, expected) in [(h, "1-2H"), (taut, "(H,5)"), (both, "1H,(H,5)")] {
            let sorts = heap
                .allocate_model_storage(vec![INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), inchi],
                    ..INCHI_SORT::default()
                }])
                .unwrap();
            let mut output = output_buffer(&mut heap, "");
            let mut overflow = 0;
            assert_eq!(
                str_H_atoms(&mut heap, sorts, &mut output, &mut overflow, OUT_T1 as i32, 0, 0, 1, 1,),
                Ok(i32::try_from(expected.len()).unwrap())
            );
            assert_eq!(visible(&heap, &output), expected);
        }

        let duplicate_sorts = heap
            .allocate_model_storage(vec![
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), both],
                    ..INCHI_SORT::default()
                },
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), both],
                    ..INCHI_SORT::default()
                },
            ])
            .unwrap();
        let mut output = output_buffer(&mut heap, "");
        let mut overflow = 0;
        assert_eq!(
            str_H_atoms(
                &mut heap,
                duplicate_sorts,
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                0,
                2,
                1,
            ),
            Ok(10)
        );
        assert_eq!(visible(&heap, &output), "2*1H,(H,5)");

        let mixed = heap
            .allocate_model_storage(vec![
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), empty1],
                    ..INCHI_SORT::default()
                },
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), both],
                    ..INCHI_SORT::default()
                },
            ])
            .unwrap();
        output = output_buffer(&mut heap, "");
        assert_eq!(
            str_H_atoms(&mut heap, mixed, &mut output, &mut overflow, OUT_T1 as i32, 0, 0, 2, 1,),
            Ok(9)
        );
        assert_eq!(visible(&heap, &output), ";1H,(H,5)");

        let empties = heap
            .allocate_model_storage(vec![
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), empty1],
                    ..INCHI_SORT::default()
                },
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), empty2],
                    ..INCHI_SORT::default()
                },
            ])
            .unwrap();
        output = output_buffer(&mut heap, "pre");
        assert_eq!(
            str_H_atoms(
                &mut heap,
                empties,
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                0,
                2,
                1,
            ),
            Ok(0)
        );
        assert_eq!(visible(&heap, &output), "pre");

        assert_eq!(
            str_H_atoms(
                &mut heap,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                0,
                1,
                1,
            ),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichiprt3__str_charge2__line_427() {
        fn charged(
            heap: &mut SourceHeap,
            charge: i32,
            atoms: i32,
            tautomer_length: i32,
            deleted: i32,
        ) -> SourceMutPointer<INChI> {
            heap.allocate_model_storage(vec![INChI {
                nNumberOfAtoms: atoms,
                lenTautomer: tautomer_length,
                nTotalCharge: charge,
                bDeleted: deleted,
                ..INChI::default()
            }])
            .unwrap()
        }

        fn taut_sorts(heap: &mut SourceHeap, values: &[SourceMutPointer<INChI>]) -> SourceMutPointer<INCHI_SORT> {
            heap.allocate_model_storage(
                values
                    .iter()
                    .map(|value| INCHI_SORT {
                        pINChI: [SourceMutPointer::null(), *value],
                        ..INCHI_SORT::default()
                    })
                    .collect(),
            )
            .unwrap()
        }

        fn non_taut_sorts(
            heap: &mut SourceHeap,
            values: &[(SourceMutPointer<INChI>, SourceMutPointer<INChI>)],
        ) -> SourceMutPointer<INCHI_SORT> {
            heap.allocate_model_storage(
                values
                    .iter()
                    .map(|(non_taut, taut)| INCHI_SORT {
                        pINChI: [*non_taut, *taut],
                        ..INCHI_SORT::default()
                    })
                    .collect(),
            )
            .unwrap()
        }

        let mut heap = SourceHeap::default();
        let plus_a = charged(&mut heap, 1, 1, 0, 0);
        let plus_b = charged(&mut heap, 1, 1, 0, 0);
        let zero = charged(&mut heap, 0, 1, 0, 0);
        let minus = charged(&mut heap, -2, 1, 0, 0);
        let first_pass = taut_sorts(&mut heap, &[plus_a, plus_b, zero, minus]);
        let mut overflow = 0;

        let mut output = output_buffer(&mut heap, "");
        assert_eq!(
            str_Charge2(
                &mut heap,
                first_pass,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                4,
                0,
                0,
                1,
            ),
            Ok(8)
        );
        assert_eq!(visible(&heap, &output), "2*+1;;-2");

        output = output_buffer(&mut heap, "");
        assert_eq!(
            str_Charge2(
                &mut heap,
                first_pass,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                4,
                0,
                0,
                0,
            ),
            Ok(9)
        );
        assert_eq!(visible(&heap, &output), "+1;+1;;-2");

        let all_zero = taut_sorts(&mut heap, &[zero, zero, zero]);
        output = output_buffer(&mut heap, "pre:");
        assert_eq!(
            str_Charge2(
                &mut heap,
                all_zero,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                3,
                0,
                0,
                1,
            ),
            Ok(2)
        );
        assert_eq!(visible(&heap, &output), "pre:;;");

        let non_a = charged(&mut heap, 1, 1, 0, 0);
        let non_b = charged(&mut heap, 1, 1, 0, 0);
        let taut_a = charged(&mut heap, 1, 1, 1, 0);
        let taut_b = charged(&mut heap, 1, 1, 1, 0);
        let non_pass = non_taut_sorts(&mut heap, &[(non_a, taut_a), (non_b, taut_b)]);
        let taut_pass = taut_sorts(&mut heap, &[taut_a, taut_b]);
        output = output_buffer(&mut heap, "");
        assert_eq!(
            str_Charge2(
                &mut heap,
                non_pass,
                taut_pass,
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                2,
                1,
                1,
                1,
            ),
            Ok(2)
        );
        assert_eq!(visible(&heap, &output), "2m");

        let explicit = charged(&mut heap, 2, 1, 0, 0);
        let different_taut = charged(&mut heap, 3, 1, 1, 0);
        let mixed_non = non_taut_sorts(&mut heap, &[(non_a, taut_a), (explicit, different_taut)]);
        let mixed_taut = taut_sorts(&mut heap, &[taut_a, different_taut]);
        output = output_buffer(&mut heap, "");
        assert_eq!(
            str_Charge2(
                &mut heap,
                mixed_non,
                mixed_taut,
                &mut output,
                &mut overflow,
                OUT_NT as i32,
                2,
                1,
                1,
                1,
            ),
            Ok(4)
        );
        assert_eq!(visible(&heap, &output), "m;+2");

        let missing_non = SourceMutPointer::null();
        let only_taut = charged(&mut heap, -1, 1, 1, 0);
        let later_non = charged(&mut heap, 2, 1, 0, 0);
        let later_taut = charged(&mut heap, 3, 1, 1, 0);
        let missing_then_explicit = non_taut_sorts(&mut heap, &[(missing_non, only_taut), (later_non, later_taut)]);
        let missing_taut_pass = taut_sorts(&mut heap, &[only_taut, later_taut]);
        output = output_buffer(&mut heap, "");
        assert_eq!(
            str_Charge2(
                &mut heap,
                missing_then_explicit,
                missing_taut_pass,
                &mut output,
                &mut overflow,
                OUT_NT as i32,
                2,
                1,
                1,
                1,
            ),
            Ok(3)
        );
        assert_eq!(visible(&heap, &output), ";+2");

        let deleted_taut = charged(&mut heap, 1, 1, 1, 1);
        let deleted_non = non_taut_sorts(&mut heap, &[(non_a, deleted_taut)]);
        let deleted_taut_pass = taut_sorts(&mut heap, &[deleted_taut]);
        output = output_buffer(&mut heap, "");
        assert_eq!(
            str_Charge2(
                &mut heap,
                deleted_non,
                deleted_taut_pass,
                &mut output,
                &mut overflow,
                OUT_NT as i32,
                1,
                1,
                1,
                1,
            ),
            Ok(2)
        );
        assert_eq!(visible(&heap, &output), "+1");

        output = output_buffer(&mut heap, "unchanged");
        assert!(matches!(
            str_Charge2(
                &mut heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                1,
                0,
                0,
                1,
            ),
            Err(SourceHeapError::NullPointer) | Err(SourceHeapError::PointerOutOfBounds)
        ));
        assert_eq!(visible(&heap, &output), "unchanged");

        let mut failure_heap = SourceHeap::default();
        let charged_value = charged(&mut failure_heap, 1, 1, 0, 0);
        let failure_sorts = taut_sorts(&mut failure_heap, &[charged_value]);
        let mut failure_output = INCHI_IOS_STRING {
            pStr: failure_heap.allocate_model_storage(vec![0_i8]).unwrap(),
            nAllocatedLength: 1,
            nUsedLength: 0,
            nPtr: 0,
        };
        let mut failure_overflow = 0;
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            str_Charge2(
                &mut failure_heap,
                failure_sorts,
                SourceMutPointer::null(),
                &mut failure_output,
                &mut failure_overflow,
                OUT_T1 as i32,
                1,
                0,
                0,
                1,
            ),
            Ok(0)
        );
        assert_eq!(visible(&failure_heap, &failure_output), "");
        assert_eq!(failure_overflow, 0);
    }

    #[test]
    fn source_port__ichiprt3__str_isosp3__line_1768() {
        fn stereo(
            heap: &mut SourceHeap,
            center: u16,
            parity: i8,
        ) -> SourceMutPointer<crate::source_types::INChI_Stereo> {
            let centers = heap.allocate_model_storage(vec![center]).unwrap();
            let parities = heap.allocate_model_storage(vec![parity]).unwrap();
            heap.allocate_model_storage(vec![crate::source_types::INChI_Stereo {
                nNumberOfStereoCenters: 1,
                nNumber: centers,
                t_parity: parities,
                ..crate::source_types::INChI_Stereo::default()
            }])
            .unwrap()
        }

        fn value(
            heap: &mut SourceHeap,
            main: Option<(u16, i8)>,
            isotopic: Option<(u16, i8)>,
            tautomer_length: i32,
        ) -> SourceMutPointer<INChI> {
            let main = main.map_or(SourceMutPointer::null(), |(center, parity)| {
                stereo(heap, center, parity)
            });
            let isotopic = isotopic.map_or(SourceMutPointer::null(), |(center, parity)| {
                stereo(heap, center, parity)
            });
            heap.allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                lenTautomer: tautomer_length,
                Stereo: main,
                StereoIsotopic: isotopic,
                ..INChI::default()
            }])
            .unwrap()
        }

        fn taut_sorts(heap: &mut SourceHeap, values: &[SourceMutPointer<INChI>]) -> SourceMutPointer<INCHI_SORT> {
            heap.allocate_model_storage(
                values
                    .iter()
                    .map(|value| INCHI_SORT {
                        pINChI: [SourceMutPointer::null(), *value],
                        ..INCHI_SORT::default()
                    })
                    .collect(),
            )
            .unwrap()
        }

        fn non_taut_sort(
            heap: &mut SourceHeap,
            non: SourceMutPointer<INChI>,
            taut: SourceMutPointer<INChI>,
        ) -> SourceMutPointer<INCHI_SORT> {
            heap.allocate_model_storage(vec![INCHI_SORT {
                pINChI: [non, taut],
                ..INCHI_SORT::default()
            }])
            .unwrap()
        }

        #[allow(clippy::too_many_arguments)]
        fn assert_output(
            heap: &mut SourceHeap,
            sorts: SourceMutPointer<INCHI_SORT>,
            sorts2: SourceMutPointer<INCHI_SORT>,
            output_type: i32,
            components: i32,
            relative_or_racemic: i32,
            second_pass: i32,
            omit: i32,
            multipliers: i32,
            expected: &str,
        ) {
            let mut output = output_buffer(heap, "");
            let mut overflow = 0;
            assert_eq!(
                str_IsoSp3(
                    heap,
                    sorts,
                    sorts2,
                    &mut output,
                    &mut overflow,
                    output_type,
                    0,
                    components,
                    relative_or_racemic,
                    second_pass,
                    omit,
                    multipliers,
                ),
                Ok(i32::try_from(expected.len()).unwrap())
            );
            assert_eq!(visible(heap, &output), expected);
            assert_eq!(overflow, 0);
        }

        let mut heap = SourceHeap::default();
        let iso_a = value(&mut heap, None, Some((1, 1)), 1);
        let iso_a_copy = value(&mut heap, None, Some((1, 1)), 1);
        let empty = value(&mut heap, None, None, 1);
        let iso_b = value(&mut heap, None, Some((3, 2)), 1);
        let first = taut_sorts(&mut heap, &[iso_a, iso_a_copy, empty, iso_b]);
        assert_output(
            &mut heap,
            first,
            SourceMutPointer::null(),
            OUT_T1 as i32,
            4,
            i32::MAX,
            0,
            0,
            1,
            "2*1-;;3+",
        );

        let equal = value(&mut heap, Some((1, 1)), Some((1, 1)), 1);
        let equal_sorts = taut_sorts(&mut heap, &[equal]);
        assert_output(
            &mut heap,
            equal_sorts,
            SourceMutPointer::null(),
            OUT_T1 as i32,
            1,
            1,
            0,
            1,
            1,
            "m",
        );

        let internal_taut = value(&mut heap, None, None, 1);
        let external_other = value(&mut heap, Some((9, 1)), None, 1);
        let external_other_sorts = taut_sorts(&mut heap, &[external_other]);
        let non_equal_main = value(&mut heap, Some((1, 1)), Some((1, 1)), 0);
        let non_equal_main_sorts = non_taut_sort(&mut heap, non_equal_main, internal_taut);
        assert_output(
            &mut heap,
            non_equal_main_sorts,
            external_other_sorts,
            OUT_NT as i32,
            1,
            1,
            1,
            1,
            1,
            "n",
        );

        let non_distinct = value(&mut heap, Some((3, 2)), Some((1, 1)), 0);
        let non_distinct_sorts = non_taut_sort(&mut heap, non_distinct, internal_taut);
        let external_main = value(&mut heap, Some((1, 1)), None, 1);
        let external_main_sorts = taut_sorts(&mut heap, &[external_main]);
        assert_output(
            &mut heap,
            non_distinct_sorts,
            external_main_sorts,
            OUT_NT as i32,
            1,
            0,
            1,
            1,
            1,
            "m",
        );

        let external_iso = value(&mut heap, Some((9, 1)), Some((1, 1)), 1);
        let external_iso_sorts = taut_sorts(&mut heap, &[external_iso]);
        assert_output(
            &mut heap,
            non_distinct_sorts,
            external_iso_sorts,
            OUT_NT as i32,
            1,
            0,
            1,
            1,
            1,
            "M",
        );

        let mut output = output_buffer(&mut heap, "held");
        let mut overflow = 1;
        assert_eq!(
            str_IsoSp3(
                &mut heap,
                first,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                4,
                1,
                0,
                0,
                1,
            ),
            Ok(0)
        );
        assert_eq!(visible(&heap, &output), "held");
    }

    #[test]
    fn source_port__ichiprt3__str_isosp2__line_1479() {
        fn stereo(
            heap: &mut SourceHeap,
            atom1: u16,
            atom2: u16,
            parity: i8,
        ) -> SourceMutPointer<crate::source_types::INChI_Stereo> {
            let first = heap.allocate_model_storage(vec![atom1]).unwrap();
            let second = heap.allocate_model_storage(vec![atom2]).unwrap();
            let parity = heap.allocate_model_storage(vec![parity]).unwrap();
            heap.allocate_model_storage(vec![crate::source_types::INChI_Stereo {
                nNumberOfStereoBonds: 1,
                nBondAtom1: first,
                nBondAtom2: second,
                b_parity: parity,
                ..crate::source_types::INChI_Stereo::default()
            }])
            .unwrap()
        }

        fn value(
            heap: &mut SourceHeap,
            main: Option<(u16, u16, i8)>,
            isotopic: Option<(u16, u16, i8)>,
            tautomer_length: i32,
        ) -> SourceMutPointer<INChI> {
            let main = main.map_or(SourceMutPointer::null(), |(a, b, p)| stereo(heap, a, b, p));
            let isotopic = isotopic.map_or(SourceMutPointer::null(), |(a, b, p)| stereo(heap, a, b, p));
            heap.allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                lenTautomer: tautomer_length,
                Stereo: main,
                StereoIsotopic: isotopic,
                ..INChI::default()
            }])
            .unwrap()
        }

        fn taut_sorts(heap: &mut SourceHeap, values: &[SourceMutPointer<INChI>]) -> SourceMutPointer<INCHI_SORT> {
            heap.allocate_model_storage(
                values
                    .iter()
                    .map(|value| INCHI_SORT {
                        pINChI: [SourceMutPointer::null(), *value],
                        ..INCHI_SORT::default()
                    })
                    .collect(),
            )
            .unwrap()
        }

        fn non_taut_sorts(
            heap: &mut SourceHeap,
            values: &[(SourceMutPointer<INChI>, SourceMutPointer<INChI>)],
        ) -> SourceMutPointer<INCHI_SORT> {
            heap.allocate_model_storage(
                values
                    .iter()
                    .map(|(non, taut)| INCHI_SORT {
                        pINChI: [*non, *taut],
                        ..INCHI_SORT::default()
                    })
                    .collect(),
            )
            .unwrap()
        }

        #[allow(clippy::too_many_arguments)]
        fn assert_output(
            heap: &mut SourceHeap,
            sorts: SourceMutPointer<INCHI_SORT>,
            sorts2: SourceMutPointer<INCHI_SORT>,
            output_type: i32,
            components: i32,
            second_pass: i32,
            omit: i32,
            multipliers: i32,
            expected: &str,
        ) {
            let mut output = output_buffer(heap, "");
            let mut overflow = 0;
            assert_eq!(
                str_IsoSp2(
                    heap,
                    sorts,
                    sorts2,
                    &mut output,
                    &mut overflow,
                    output_type,
                    0,
                    components,
                    second_pass,
                    omit,
                    multipliers,
                ),
                Ok(i32::try_from(expected.len()).unwrap())
            );
            assert_eq!(visible(heap, &output), expected);
            assert_eq!(overflow, 0);
        }

        let mut heap = SourceHeap::default();
        let iso_a = value(&mut heap, None, Some((1, 2, 1)), 1);
        let iso_a_copy = value(&mut heap, None, Some((1, 2, 1)), 1);
        let empty = value(&mut heap, None, None, 1);
        let iso_b = value(&mut heap, None, Some((3, 4, 2)), 1);
        let first = taut_sorts(&mut heap, &[iso_a, iso_a_copy, empty, iso_b]);
        assert_output(
            &mut heap,
            first,
            SourceMutPointer::null(),
            OUT_T1 as i32,
            4,
            0,
            0,
            1,
            "2*1-2-;;3-4+",
        );

        let main_equal = value(&mut heap, Some((1, 2, 1)), Some((1, 2, 1)), 1);
        let main_equal_sorts = taut_sorts(&mut heap, &[main_equal]);
        assert_output(
            &mut heap,
            main_equal_sorts,
            SourceMutPointer::null(),
            OUT_T1 as i32,
            1,
            0,
            1,
            1,
            "m",
        );

        let internal_taut = value(&mut heap, None, None, 1);
        let external_b = value(&mut heap, Some((5, 6, 1)), None, 1);
        let non_a = value(&mut heap, Some((1, 2, 1)), Some((1, 2, 1)), 0);
        let non_a_sorts = non_taut_sorts(&mut heap, &[(non_a, internal_taut)]);
        let external_b_sorts = taut_sorts(&mut heap, &[external_b]);
        assert_output(&mut heap, non_a_sorts, external_b_sorts, OUT_NT as i32, 1, 1, 1, 1, "n");

        let non_b = value(&mut heap, Some((3, 4, 2)), Some((1, 2, 1)), 0);
        let external_main = value(&mut heap, Some((1, 2, 1)), None, 1);
        let non_b_sorts = non_taut_sorts(&mut heap, &[(non_b, internal_taut)]);
        let external_main_sorts = taut_sorts(&mut heap, &[external_main]);
        assert_output(
            &mut heap,
            non_b_sorts,
            external_main_sorts,
            OUT_NT as i32,
            1,
            1,
            1,
            1,
            "m",
        );

        let external_iso = value(&mut heap, Some((5, 6, 1)), Some((1, 2, 1)), 1);
        let external_iso_sorts = taut_sorts(&mut heap, &[external_iso]);
        assert_output(
            &mut heap,
            non_b_sorts,
            external_iso_sorts,
            OUT_NT as i32,
            1,
            1,
            1,
            1,
            "M",
        );

        let mut output = output_buffer(&mut heap, "held");
        let mut overflow = 1;
        assert_eq!(
            str_IsoSp2(
                &mut heap,
                first,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                4,
                0,
                0,
                1,
            ),
            Ok(0)
        );
        assert_eq!(visible(&heap, &output), "held");
    }

    #[test]
    fn source_port__ichiprt3__str_isoatoms__line_1229() {
        fn isotopic_inchi(
            heap: &mut SourceHeap,
            atom: Option<crate::source_types::INChI_IsotopicAtom>,
            group: Option<crate::source_types::INChI_IsotopicTGroup>,
            tautomer_length: i32,
        ) -> SourceMutPointer<INChI> {
            let (atom_count, atom_pointer) = atom.map_or((0, SourceMutPointer::null()), |value| {
                (1, heap.allocate_model_storage(vec![value]).unwrap())
            });
            let (group_count, group_pointer) = group.map_or((0, SourceMutPointer::null()), |value| {
                (1, heap.allocate_model_storage(vec![value]).unwrap())
            });
            heap.allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                lenTautomer: tautomer_length,
                nNumberOfIsotopicAtoms: atom_count,
                IsotopicAtom: atom_pointer,
                nNumberOfIsotopicTGroups: group_count,
                IsotopicTGroup: group_pointer,
                ..INChI::default()
            }])
            .unwrap()
        }

        fn taut_sorts(heap: &mut SourceHeap, values: &[SourceMutPointer<INChI>]) -> SourceMutPointer<INCHI_SORT> {
            heap.allocate_model_storage(
                values
                    .iter()
                    .map(|value| INCHI_SORT {
                        pINChI: [SourceMutPointer::null(), *value],
                        ..INCHI_SORT::default()
                    })
                    .collect(),
            )
            .unwrap()
        }

        fn non_taut_sorts(
            heap: &mut SourceHeap,
            values: &[(SourceMutPointer<INChI>, SourceMutPointer<INChI>)],
        ) -> SourceMutPointer<INCHI_SORT> {
            heap.allocate_model_storage(
                values
                    .iter()
                    .map(|(non_taut, taut)| INCHI_SORT {
                        pINChI: [*non_taut, *taut],
                        ..INCHI_SORT::default()
                    })
                    .collect(),
            )
            .unwrap()
        }

        #[allow(clippy::too_many_arguments)]
        fn assert_output(
            heap: &mut SourceHeap,
            sorts: SourceMutPointer<INCHI_SORT>,
            sorts2: SourceMutPointer<INCHI_SORT>,
            output_type: i32,
            components: i32,
            mode: i32,
            abc: i32,
            second_pass: i32,
            omit: i32,
            multipliers: i32,
            expected: &str,
        ) {
            let mut output = output_buffer(heap, "");
            let mut overflow = 0;
            assert_eq!(
                str_IsoAtoms(
                    heap,
                    sorts,
                    sorts2,
                    &mut output,
                    &mut overflow,
                    output_type,
                    mode,
                    components,
                    abc,
                    second_pass,
                    omit,
                    multipliers,
                ),
                Ok(i32::try_from(expected.len()).unwrap())
            );
            assert_eq!(visible(heap, &output), expected);
            assert_eq!(overflow, 0);
        }

        let atom_a = crate::source_types::INChI_IsotopicAtom {
            nAtomNumber: 1,
            nIsoDifference: 1,
            ..crate::source_types::INChI_IsotopicAtom::default()
        };
        let atom_b = crate::source_types::INChI_IsotopicAtom {
            nAtomNumber: 2,
            nIsoDifference: -2,
            ..crate::source_types::INChI_IsotopicAtom::default()
        };
        let group = crate::source_types::INChI_IsotopicTGroup {
            nTGroupNumber: 1,
            nNum_H: 1,
            ..crate::source_types::INChI_IsotopicTGroup::default()
        };

        let mut heap = SourceHeap::default();
        let a = isotopic_inchi(&mut heap, Some(atom_a.clone()), None, 1);
        let a_copy = isotopic_inchi(&mut heap, Some(atom_a), None, 1);
        let empty = isotopic_inchi(&mut heap, None, None, 1);
        let b_group = isotopic_inchi(&mut heap, Some(atom_b), Some(group), 1);
        let first_pass = taut_sorts(&mut heap, &[a, a_copy, empty, b_group]);
        assert_output(
            &mut heap,
            first_pass,
            SourceMutPointer::null(),
            OUT_T1 as i32,
            4,
            0,
            0,
            0,
            0,
            1,
            "2*1+0;;2-2(1H)",
        );
        assert_output(
            &mut heap,
            first_pass,
            SourceMutPointer::null(),
            OUT_T1 as i32,
            4,
            0,
            0,
            0,
            0,
            0,
            "1+0;1+0;;2-2(1H)",
        );
        let abc_sorts = taut_sorts(&mut heap, &[b_group]);
        assert_output(
            &mut heap,
            abc_sorts,
            SourceMutPointer::null(),
            OUT_T1 as i32,
            1,
            crate::source_types::CT_MODE_ABC_NUMBERS as i32,
            1,
            0,
            0,
            1,
            "B-2,A1h",
        );

        let non_a = isotopic_inchi(
            &mut heap,
            Some(crate::source_types::INChI_IsotopicAtom {
                nAtomNumber: 1,
                nIsoDifference: 1,
                ..crate::source_types::INChI_IsotopicAtom::default()
            }),
            None,
            0,
        );
        let non_a_value = heap.slice(non_a.as_const()).unwrap().to_vec();
        let non_a_copy = heap.allocate_model_storage(non_a_value).unwrap();
        let taut_a = a;
        let taut_a_copy = a_copy;
        let equal_non = non_taut_sorts(&mut heap, &[(non_a, taut_a), (non_a_copy, taut_a_copy)]);
        let equal_taut = taut_sorts(&mut heap, &[taut_a, taut_a_copy]);
        assert_output(&mut heap, equal_non, equal_taut, OUT_NT as i32, 2, 0, 0, 1, 1, 1, "2m");

        let mut output = output_buffer(&mut heap, "held");
        let mut overflow = 1;
        assert_eq!(
            str_IsoAtoms(
                &mut heap,
                first_pass,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                4,
                0,
                0,
                0,
                1,
            ),
            Ok(0)
        );
        assert_eq!(visible(&heap, &output), "held");
    }

    #[test]
    fn source_port__ichiprt3__str_sp2__line_725() {
        fn stereo_inchi(
            heap: &mut SourceHeap,
            atom1: &[u16],
            atom2: &[u16],
            parity: &[i8],
            tautomer_length: i32,
        ) -> SourceMutPointer<INChI> {
            let n_bond_atom1 = heap.allocate_model_storage(atom1.to_vec()).unwrap();
            let n_bond_atom2 = heap.allocate_model_storage(atom2.to_vec()).unwrap();
            let b_parity = heap.allocate_model_storage(parity.to_vec()).unwrap();
            let stereo = heap
                .allocate_model_storage(vec![crate::source_types::INChI_Stereo {
                    nNumberOfStereoBonds: i32::try_from(atom1.len()).unwrap(),
                    nBondAtom1: n_bond_atom1,
                    nBondAtom2: n_bond_atom2,
                    b_parity,
                    ..crate::source_types::INChI_Stereo::default()
                }])
                .unwrap();
            heap.allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                lenTautomer: tautomer_length,
                Stereo: stereo,
                ..INChI::default()
            }])
            .unwrap()
        }

        fn taut_sorts(heap: &mut SourceHeap, values: &[SourceMutPointer<INChI>]) -> SourceMutPointer<INCHI_SORT> {
            heap.allocate_model_storage(
                values
                    .iter()
                    .map(|value| INCHI_SORT {
                        pINChI: [SourceMutPointer::null(), *value],
                        ..INCHI_SORT::default()
                    })
                    .collect(),
            )
            .unwrap()
        }

        fn non_taut_sorts(
            heap: &mut SourceHeap,
            values: &[(SourceMutPointer<INChI>, SourceMutPointer<INChI>)],
        ) -> SourceMutPointer<INCHI_SORT> {
            heap.allocate_model_storage(
                values
                    .iter()
                    .map(|(non_taut, taut)| INCHI_SORT {
                        pINChI: [*non_taut, *taut],
                        ..INCHI_SORT::default()
                    })
                    .collect(),
            )
            .unwrap()
        }

        fn assert_output(
            heap: &mut SourceHeap,
            sorts: SourceMutPointer<INCHI_SORT>,
            sorts2: SourceMutPointer<INCHI_SORT>,
            output_type: i32,
            components: i32,
            second_pass: i32,
            omit: i32,
            multipliers: i32,
            expected: &str,
        ) {
            let mut output = output_buffer(heap, "");
            let mut overflow = 0;
            assert_eq!(
                str_Sp2(
                    heap,
                    sorts,
                    sorts2,
                    &mut output,
                    &mut overflow,
                    output_type,
                    0,
                    components,
                    second_pass,
                    omit,
                    multipliers,
                ),
                Ok(i32::try_from(expected.len()).unwrap())
            );
            assert_eq!(visible(heap, &output), expected);
            assert_eq!(overflow, 0);
        }

        let mut heap = SourceHeap::default();
        let stereo_a = stereo_inchi(&mut heap, &[1], &[2], &[1], 1);
        let stereo_a_copy = stereo_inchi(&mut heap, &[1], &[2], &[1], 1);
        let empty = stereo_inchi(&mut heap, &[], &[], &[], 1);
        let stereo_b = stereo_inchi(&mut heap, &[3], &[4], &[2], 1);
        let first_pass = taut_sorts(&mut heap, &[stereo_a, stereo_a_copy, empty, stereo_b]);
        assert_output(
            &mut heap,
            first_pass,
            SourceMutPointer::null(),
            OUT_T1 as i32,
            4,
            0,
            0,
            1,
            "2*1-2-;;3-4+",
        );
        assert_output(
            &mut heap,
            first_pass,
            SourceMutPointer::null(),
            OUT_T1 as i32,
            4,
            0,
            0,
            0,
            "1-2-;1-2-;;3-4+",
        );

        let non_a = stereo_inchi(&mut heap, &[1], &[2], &[1], 0);
        let non_a_copy = stereo_inchi(&mut heap, &[1], &[2], &[1], 0);
        let taut_a = stereo_inchi(&mut heap, &[1], &[2], &[1], 1);
        let taut_a_copy = stereo_inchi(&mut heap, &[1], &[2], &[1], 1);
        let equal_non = non_taut_sorts(&mut heap, &[(non_a, taut_a), (non_a_copy, taut_a_copy)]);
        let equal_taut = taut_sorts(&mut heap, &[taut_a, taut_a_copy]);
        assert_output(&mut heap, equal_non, equal_taut, OUT_NT as i32, 2, 1, 1, 1, "2m");

        let non_b = stereo_inchi(&mut heap, &[3], &[4], &[2], 0);
        let taut_other = stereo_inchi(&mut heap, &[5], &[6], &[1], 1);
        let mixed_non = non_taut_sorts(&mut heap, &[(non_a, taut_a), (non_b, taut_other)]);
        let mixed_taut = taut_sorts(&mut heap, &[taut_a, taut_other]);
        assert_output(&mut heap, mixed_non, mixed_taut, OUT_NT as i32, 2, 1, 1, 1, "m;3-4+");

        let missing_then_b = non_taut_sorts(&mut heap, &[(SourceMutPointer::null(), taut_a), (non_b, taut_other)]);
        let missing_taut = taut_sorts(&mut heap, &[taut_a, taut_other]);
        assert_output(
            &mut heap,
            missing_then_b,
            missing_taut,
            OUT_NT as i32,
            2,
            1,
            1,
            1,
            ";3-4+",
        );

        let mut output = output_buffer(&mut heap, "held");
        let mut overflow = 1;
        assert_eq!(
            str_Sp2(
                &mut heap,
                first_pass,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                4,
                0,
                0,
                1,
            ),
            Ok(0)
        );
        assert_eq!(visible(&heap, &output), "held");
    }

    #[test]
    fn source_port__ichiprt3__str_sp3__line_957() {
        fn stereo_inchi(
            heap: &mut SourceHeap,
            centers: &[u16],
            parity: &[i8],
            tautomer_length: i32,
        ) -> SourceMutPointer<INChI> {
            let n_number = heap.allocate_model_storage(centers.to_vec()).unwrap();
            let t_parity = heap.allocate_model_storage(parity.to_vec()).unwrap();
            let stereo = heap
                .allocate_model_storage(vec![crate::source_types::INChI_Stereo {
                    nNumberOfStereoCenters: i32::try_from(centers.len()).unwrap(),
                    nNumber: n_number,
                    t_parity,
                    ..crate::source_types::INChI_Stereo::default()
                }])
                .unwrap();
            heap.allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                lenTautomer: tautomer_length,
                Stereo: stereo,
                ..INChI::default()
            }])
            .unwrap()
        }

        fn taut_sorts(heap: &mut SourceHeap, values: &[SourceMutPointer<INChI>]) -> SourceMutPointer<INCHI_SORT> {
            heap.allocate_model_storage(
                values
                    .iter()
                    .map(|value| INCHI_SORT {
                        pINChI: [SourceMutPointer::null(), *value],
                        ..INCHI_SORT::default()
                    })
                    .collect(),
            )
            .unwrap()
        }

        fn non_taut_sorts(
            heap: &mut SourceHeap,
            values: &[(SourceMutPointer<INChI>, SourceMutPointer<INChI>)],
        ) -> SourceMutPointer<INCHI_SORT> {
            heap.allocate_model_storage(
                values
                    .iter()
                    .map(|(non_taut, taut)| INCHI_SORT {
                        pINChI: [*non_taut, *taut],
                        ..INCHI_SORT::default()
                    })
                    .collect(),
            )
            .unwrap()
        }

        fn assert_output(
            heap: &mut SourceHeap,
            sorts: SourceMutPointer<INCHI_SORT>,
            sorts2: SourceMutPointer<INCHI_SORT>,
            output_type: i32,
            components: i32,
            second_pass: i32,
            omit: i32,
            multipliers: i32,
            expected: &str,
        ) {
            let mut output = output_buffer(heap, "");
            let mut overflow = 0;
            assert_eq!(
                str_Sp3(
                    heap,
                    sorts,
                    sorts2,
                    &mut output,
                    &mut overflow,
                    output_type,
                    0,
                    components,
                    1,
                    second_pass,
                    omit,
                    multipliers,
                ),
                Ok(i32::try_from(expected.len()).unwrap())
            );
            assert_eq!(visible(heap, &output), expected);
            assert_eq!(overflow, 0);
        }

        let mut heap = SourceHeap::default();
        let stereo_a = stereo_inchi(&mut heap, &[1], &[1], 1);
        let stereo_a_copy = stereo_inchi(&mut heap, &[1], &[1], 1);
        let empty = stereo_inchi(&mut heap, &[], &[], 1);
        let stereo_b = stereo_inchi(&mut heap, &[3], &[2], 1);
        let first_pass = taut_sorts(&mut heap, &[stereo_a, stereo_a_copy, empty, stereo_b]);
        assert_output(
            &mut heap,
            first_pass,
            SourceMutPointer::null(),
            OUT_T1 as i32,
            4,
            0,
            0,
            1,
            "2*1-;;3+",
        );
        assert_output(
            &mut heap,
            first_pass,
            SourceMutPointer::null(),
            OUT_T1 as i32,
            4,
            0,
            0,
            0,
            "1-;1-;;3+",
        );

        let non_a = stereo_inchi(&mut heap, &[1], &[1], 0);
        let non_a_copy = stereo_inchi(&mut heap, &[1], &[1], 0);
        let taut_a = stereo_inchi(&mut heap, &[1], &[1], 1);
        let taut_a_copy = stereo_inchi(&mut heap, &[1], &[1], 1);
        let equal_non = non_taut_sorts(&mut heap, &[(non_a, taut_a), (non_a_copy, taut_a_copy)]);
        let equal_taut = taut_sorts(&mut heap, &[taut_a, taut_a_copy]);
        assert_output(&mut heap, equal_non, equal_taut, OUT_NT as i32, 2, 1, 1, 1, "2m");

        let non_b = stereo_inchi(&mut heap, &[3], &[2], 0);
        let taut_other = stereo_inchi(&mut heap, &[5], &[1], 1);
        let mixed_non = non_taut_sorts(&mut heap, &[(non_a, taut_a), (non_b, taut_other)]);
        let mixed_taut = taut_sorts(&mut heap, &[taut_a, taut_other]);
        assert_output(&mut heap, mixed_non, mixed_taut, OUT_NT as i32, 2, 1, 1, 1, "m;3+");

        let missing_then_b = non_taut_sorts(&mut heap, &[(SourceMutPointer::null(), taut_a), (non_b, taut_other)]);
        let missing_taut = taut_sorts(&mut heap, &[taut_a, taut_other]);
        assert_output(
            &mut heap,
            missing_then_b,
            missing_taut,
            OUT_NT as i32,
            2,
            1,
            1,
            1,
            ";3+",
        );

        let mut output = output_buffer(&mut heap, "held");
        let mut overflow = 1;
        assert_eq!(
            str_Sp3(
                &mut heap,
                first_pass,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                4,
                1,
                0,
                0,
                1,
            ),
            Ok(0)
        );
        assert_eq!(visible(&heap, &output), "held");
    }

    #[test]
    fn source_port__ichiprt3__bin_auxtauttrans__line_4091() {
        fn inchi(heap: &mut SourceHeap, atoms: i32) -> SourceMutPointer<INChI> {
            heap.allocate_model_storage(vec![INChI {
                nNumberOfAtoms: atoms,
                lenTautomer: 1,
                ..INChI::default()
            }])
            .unwrap()
        }

        fn sorts(
            heap: &mut SourceHeap,
            value: SourceMutPointer<INChI>,
            orders: &[i16],
        ) -> SourceMutPointer<INCHI_SORT> {
            heap.allocate_model_storage(
                orders
                    .iter()
                    .map(|order| INCHI_SORT {
                        pINChI: [SourceMutPointer::null(), value],
                        ord_number: *order,
                        ..INCHI_SORT::default()
                    })
                    .collect(),
            )
            .unwrap()
        }

        let mut heap = SourceHeap::default();
        let present = inchi(&mut heap, 1);
        let non = sorts(&mut heap, present, &[2, 0, 1]);
        let taut = sorts(&mut heap, present, &[0, 1, 2]);
        let old_non = heap.allocate_model_storage(vec![91_u16]).unwrap();
        let old_taut = heap.allocate_model_storage(vec![92_u16]).unwrap();
        let mut transposed_non = old_non;
        let mut transposed_taut = old_taut;
        assert_eq!(
            bin_AuxTautTrans(
                &mut heap,
                non,
                taut,
                &mut transposed_non,
                &mut transposed_taut,
                OUT_T1 as i32,
                3,
            ),
            Ok(1)
        );
        assert_eq!(&heap.slice(transposed_non.as_const()).unwrap()[..4], &[2, 3, 1, 0]);
        assert_eq!(&heap.slice(transposed_taut.as_const()).unwrap()[..4], &[0, 2, 3, 1]);

        let equal = sorts(&mut heap, present, &[0, 1, 2]);
        transposed_non = old_non;
        transposed_taut = old_taut;
        assert_eq!(
            bin_AuxTautTrans(
                &mut heap,
                equal,
                taut,
                &mut transposed_non,
                &mut transposed_taut,
                OUT_T1 as i32,
                3,
            ),
            Ok(0)
        );
        assert_eq!(transposed_non, old_non);
        assert_eq!(transposed_taut, old_taut);

        let mut first_fail_heap = SourceHeap::default();
        let present = inchi(&mut first_fail_heap, 1);
        let non = sorts(&mut first_fail_heap, present, &[1]);
        let taut = sorts(&mut first_fail_heap, present, &[0]);
        let old_non = first_fail_heap.allocate_model_storage(vec![71_u16]).unwrap();
        let old_taut = first_fail_heap.allocate_model_storage(vec![72_u16]).unwrap();
        let mut transposed_non = old_non;
        let mut transposed_taut = old_taut;
        first_fail_heap.fail_after_allocations(0);
        assert_eq!(
            bin_AuxTautTrans(
                &mut first_fail_heap,
                non,
                taut,
                &mut transposed_non,
                &mut transposed_taut,
                OUT_T1 as i32,
                1,
            ),
            Ok(-1)
        );
        assert_eq!(transposed_non, old_non);
        assert_eq!(transposed_taut, old_taut);

        let mut one_fail_heap = SourceHeap::default();
        let present = inchi(&mut one_fail_heap, 1);
        let non = sorts(&mut one_fail_heap, present, &[1]);
        let taut = sorts(&mut one_fail_heap, present, &[0]);
        let old_non = one_fail_heap.allocate_model_storage(vec![81_u16]).unwrap();
        let old_taut = one_fail_heap.allocate_model_storage(vec![82_u16]).unwrap();
        let mut transposed_non = old_non;
        let mut transposed_taut = old_taut;
        one_fail_heap.fail_after_allocations(1);
        assert_eq!(
            bin_AuxTautTrans(
                &mut one_fail_heap,
                non,
                taut,
                &mut transposed_non,
                &mut transposed_taut,
                OUT_T1 as i32,
                1,
            ),
            Ok(-1)
        );
        assert_eq!(transposed_non, old_non);
        assert_eq!(transposed_taut, old_taut);

        let mut untouched_non = old_non;
        let mut untouched_taut = old_taut;
        assert_eq!(
            bin_AuxTautTrans(
                &mut one_fail_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut untouched_non,
                &mut untouched_taut,
                OUT_T1 as i32,
                -1,
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichiprt3__str_isostereoabsinv__line_2058() {
        fn isotopic_stereo_inchi(heap: &mut SourceHeap, inversion: Option<i32>) -> SourceMutPointer<INChI> {
            let stereo = inversion.map_or(SourceMutPointer::null(), |value| {
                heap.allocate_model_storage(vec![crate::source_types::INChI_Stereo {
                    nCompInv2Abs: value,
                    ..crate::source_types::INChI_Stereo::default()
                }])
                .unwrap()
            });
            heap.allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                lenTautomer: 1,
                StereoIsotopic: stereo,
                ..INChI::default()
            }])
            .unwrap()
        }

        let mut heap = SourceHeap::default();
        let negative = isotopic_stereo_inchi(&mut heap, Some(i32::MIN));
        let positive = isotopic_stereo_inchi(&mut heap, Some(i32::MAX));
        let zero = isotopic_stereo_inchi(&mut heap, Some(0));
        let no_stereo = isotopic_stereo_inchi(&mut heap, None);
        let sorts = heap
            .allocate_model_storage(vec![
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), negative],
                    ..INCHI_SORT::default()
                },
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), positive],
                    ..INCHI_SORT::default()
                },
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), zero],
                    ..INCHI_SORT::default()
                },
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), no_stereo],
                    ..INCHI_SORT::default()
                },
                INCHI_SORT::default(),
            ])
            .unwrap();
        let mut output = output_buffer(&mut heap, "pre:");
        let mut overflow = 0;
        assert_eq!(
            str_IsoStereoAbsInv(&mut heap, sorts, &mut output, &mut overflow, OUT_T1 as i32, 5,),
            Ok(5)
        );
        assert_eq!(visible(&heap, &output), "pre:10...");

        output = output_buffer(&mut heap, "held");
        overflow = 1;
        assert_eq!(
            str_IsoStereoAbsInv(
                &mut heap,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                i32::MAX,
            ),
            Ok(0)
        );
        assert_eq!(visible(&heap, &output), "held");

        overflow = 0;
        assert_eq!(
            str_IsoStereoAbsInv(
                &mut heap,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                -1,
            ),
            Ok(0)
        );
        assert!(matches!(
            str_IsoStereoAbsInv(
                &mut heap,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                1,
            ),
            Err(SourceHeapError::NullPointer)
        ));
    }

    #[test]
    fn source_port__ichiprt3__str_stereoabsinv__line_1191() {
        fn stereo_inchi(heap: &mut SourceHeap, inversion: Option<i32>) -> SourceMutPointer<INChI> {
            let stereo = inversion.map_or(SourceMutPointer::null(), |value| {
                heap.allocate_model_storage(vec![crate::source_types::INChI_Stereo {
                    nCompInv2Abs: value,
                    ..crate::source_types::INChI_Stereo::default()
                }])
                .unwrap()
            });
            heap.allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                lenTautomer: 1,
                Stereo: stereo,
                ..INChI::default()
            }])
            .unwrap()
        }

        let mut heap = SourceHeap::default();
        let negative = stereo_inchi(&mut heap, Some(i32::MIN));
        let positive = stereo_inchi(&mut heap, Some(i32::MAX));
        let zero = stereo_inchi(&mut heap, Some(0));
        let no_stereo = stereo_inchi(&mut heap, None);
        let sorts = heap
            .allocate_model_storage(vec![
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), negative],
                    ..INCHI_SORT::default()
                },
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), positive],
                    ..INCHI_SORT::default()
                },
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), zero],
                    ..INCHI_SORT::default()
                },
                INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), no_stereo],
                    ..INCHI_SORT::default()
                },
                INCHI_SORT::default(),
            ])
            .unwrap();
        let mut output = output_buffer(&mut heap, "pre:");
        let mut overflow = 0;
        assert_eq!(
            str_StereoAbsInv(&mut heap, sorts, &mut output, &mut overflow, OUT_T1 as i32, 5,),
            Ok(5)
        );
        assert_eq!(visible(&heap, &output), "pre:10...");

        output = output_buffer(&mut heap, "held");
        overflow = 1;
        assert_eq!(
            str_StereoAbsInv(
                &mut heap,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                i32::MAX,
            ),
            Ok(0)
        );
        assert_eq!(visible(&heap, &output), "held");

        overflow = 0;
        assert_eq!(
            str_StereoAbsInv(
                &mut heap,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                -1,
            ),
            Ok(0)
        );
        assert!(matches!(
            str_StereoAbsInv(
                &mut heap,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                1,
            ),
            Err(SourceHeapError::NullPointer)
        ));
    }

    #[test]
    fn source_port__ichiprt3__str_auxtauttrans__line_4187() {
        fn run(
            mapping: Option<&[u16]>,
            mode: i32,
            components: i32,
            initial_overflow: i32,
        ) -> (i32, String, i32, bool, bool) {
            let mut heap = SourceHeap::default();
            let work = mapping.map_or(SourceMutPointer::null(), |values| {
                heap.allocate_model_storage(vec![91_u16; values.len()]).unwrap()
            });
            let transposition = mapping.map_or(SourceMutPointer::null(), |values| {
                heap.allocate_model_storage(values.to_vec()).unwrap()
            });
            let mut output = output_buffer(&mut heap, "pre:");
            let mut overflow = initial_overflow;
            let result = str_AuxTautTrans(
                &mut heap,
                work,
                transposition,
                &mut output,
                &mut overflow,
                mode,
                components,
            )
            .unwrap();
            let text = visible(&heap, &output);
            let work_freed = heap.slice(work.as_const()).is_err();
            let transposition_freed = heap.slice(transposition.as_const()).is_err();
            (result, text, overflow, work_freed, transposition_freed)
        }

        assert_eq!(
            run(Some(&[0, 2, 1, 4, 3]), 0, 4, 0),
            (10, "pre:(1,2)(3,4)".to_owned(), 0, true, true)
        );
        assert_eq!(
            run(
                Some(&[0, 2, 1, 4, 3]),
                crate::source_types::CT_MODE_ABC_NUMBERS as i32,
                4,
                0,
            ),
            (8, "pre:(AB)(CD)".to_owned(), 0, true, true)
        );
        assert_eq!(run(Some(&[0, 2, 0]), 0, 2, 0), (3, "pre:(1)".to_owned(), 0, true, true));
        assert_eq!(run(Some(&[0, 2, 1]), 0, 2, 1), (0, "pre:".to_owned(), 1, true, true));
        assert_eq!(run(None, 0, i32::MIN, 0), (0, "pre:".to_owned(), 0, true, true));

        let mut heap = SourceHeap::default();
        let only_work = heap.allocate_model_storage(vec![7_u16]).unwrap();
        let mut output = output_buffer(&mut heap, "");
        let mut overflow = 0;
        assert_eq!(
            str_AuxTautTrans(
                &mut heap,
                only_work,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                0,
                i32::MAX,
            ),
            Ok(0)
        );
        assert!(heap.slice(only_work.as_const()).is_err());
    }

    #[test]
    fn source_port__ichiprt3__str_hillformula2__line_137() {
        fn formula_inchi(heap: &mut SourceHeap, formula: &str, tautomer: i32, deleted: i32) -> SourceMutPointer<INChI> {
            let formula = text(heap, formula);
            heap.allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                lenTautomer: tautomer,
                bDeleted: deleted,
                szHillFormula: formula,
                ..INChI::default()
            }])
            .unwrap()
        }
        fn run(
            fixed: &[(&str, i32)],
            taut: &[(&str, i32)],
            multipliers: i32,
            initial_overflow: i32,
        ) -> (i32, String, i32) {
            let mut heap = SourceHeap::default();
            let fixed_values: Vec<_> = fixed
                .iter()
                .map(|(formula, deleted)| formula_inchi(&mut heap, formula, 0, *deleted))
                .collect();
            let taut_values: Vec<_> = taut
                .iter()
                .map(|(formula, deleted)| formula_inchi(&mut heap, formula, 1, *deleted))
                .collect();
            let fixed_sorts = heap
                .allocate_model_storage(
                    fixed_values
                        .iter()
                        .map(|value| INCHI_SORT {
                            pINChI: [*value, SourceMutPointer::null()],
                            ..INCHI_SORT::default()
                        })
                        .collect(),
                )
                .unwrap();
            let taut_sorts = heap
                .allocate_model_storage(
                    taut_values
                        .iter()
                        .map(|value| INCHI_SORT {
                            pINChI: [SourceMutPointer::null(), *value],
                            ..INCHI_SORT::default()
                        })
                        .collect(),
                )
                .unwrap();
            let mut output = output_buffer(&mut heap, "pre:");
            let mut overflow = initial_overflow;
            let result = str_HillFormula2(
                &mut heap,
                fixed_sorts,
                taut_sorts,
                &mut output,
                &mut overflow,
                OUT_NN as i32,
                i32::try_from(fixed.len()).unwrap(),
                multipliers,
            )
            .unwrap();
            (result, visible(&heap, &output), overflow)
        }

        assert_eq!(
            run(&[("C2H6", 0), ("H2O", 0)], &[("C2H6", 0), ("H2O", 0)], 1, 0),
            (0, "pre:".to_owned(), 0)
        );
        assert_eq!(run(&[("C2H6", 0)], &[("CH4", 0)], 1, 0), (4, "pre:C2H6".to_owned(), 0));
        assert_eq!(
            run(&[("C", 0), ("C", 0), ("O", 0)], &[("N", 0), ("N", 0), ("N", 0)], 1, 0),
            (4, "pre:2C.O".to_owned(), 0)
        );
        assert_eq!(
            run(&[("", 0), ("O", 0)], &[("N", 0), ("N", 0)], 0, 0),
            (2, "pre:.O".to_owned(), 0)
        );
        assert_eq!(run(&[("C", 0)], &[("C", 1)], 1, 0), (1, "pre:C".to_owned(), 0));
        assert_eq!(run(&[("C", 0)], &[("N", 0)], 1, 7), (0, "pre:".to_owned(), 7));
    }

    #[test]
    fn source_port__ichiprt3__str_fixedh_atoms__line_618() {
        fn run(values: &[&[i8]], multipliers: i32, initial_overflow: i32) -> (i32, String, i32) {
            let mut heap = SourceHeap::default();
            let mut inchis = Vec::new();
            for hydrogens in values {
                let pointer = heap.allocate_model_storage(hydrogens.to_vec()).unwrap();
                inchis.push(
                    heap.allocate_model_storage(vec![INChI {
                        nNumberOfAtoms: i32::try_from(hydrogens.len()).unwrap(),
                        lenTautomer: 0,
                        nNum_H_fixed: pointer,
                        ..INChI::default()
                    }])
                    .unwrap(),
                );
            }
            let sorts = heap
                .allocate_model_storage(
                    inchis
                        .iter()
                        .map(|value| INCHI_SORT {
                            pINChI: [*value, SourceMutPointer::null()],
                            ..INCHI_SORT::default()
                        })
                        .collect(),
                )
                .unwrap();
            let mut output = output_buffer(&mut heap, "pre:");
            let mut overflow = initial_overflow;
            let result = str_FixedH_atoms(
                &mut heap,
                sorts,
                &mut output,
                &mut overflow,
                OUT_NN as i32,
                0,
                i32::try_from(values.len()).unwrap(),
                multipliers,
            )
            .unwrap();
            (result, visible(&heap, &output), overflow)
        }

        assert_eq!(run(&[&[1, 0, 2]], 1, 0), (6, "pre:1H,3H2".to_owned(), 0));
        assert_eq!(run(&[&[1, 0, 2], &[1, 0, 2]], 1, 0), (8, "pre:2*1H,3H2".to_owned(), 0));
        assert_eq!(run(&[&[0, 0], &[0, 0]], 1, 0), (0, "pre:".to_owned(), 0));
        assert_eq!(run(&[&[0], &[1]], 0, 0), (3, "pre:;1H".to_owned(), 0));
        assert_eq!(run(&[&[1]], 1, 5), (0, "pre:".to_owned(), 5));
    }

    #[test]
    fn source_port__ichiprt3__str_auxnumb__line_3818() {
        use crate::source_types::{INChI, INChI_Aux, OUT_NT, OUT_T1, TAUT_NON, TAUT_YES};

        fn component(heap: &mut SourceHeap, taut_numbers: &[u16], non_taut_numbers: Option<&[u16]>) -> INCHI_SORT {
            let taut_number_pointer = heap.allocate_model_storage(taut_numbers.to_vec()).unwrap();
            let taut_aux = heap
                .allocate_model_storage(vec![INChI_Aux {
                    nNumberOfAtoms: taut_numbers.len() as i32,
                    nOrigAtNosInCanonOrd: taut_number_pointer,
                    ..INChI_Aux::default()
                }])
                .unwrap();
            let taut_inchi = heap
                .allocate_model_storage(vec![INChI {
                    nNumberOfAtoms: taut_numbers.len() as i32,
                    lenTautomer: i32::from(non_taut_numbers.is_some()),
                    ..INChI::default()
                }])
                .unwrap();
            let mut sort = INCHI_SORT::default();
            sort.pINChI[TAUT_YES as usize] = taut_inchi;
            sort.pINChI_Aux[TAUT_YES as usize] = taut_aux;
            if let Some(non_taut_numbers) = non_taut_numbers {
                let numbers = heap.allocate_model_storage(non_taut_numbers.to_vec()).unwrap();
                let aux = heap
                    .allocate_model_storage(vec![INChI_Aux {
                        nNumberOfAtoms: non_taut_numbers.len() as i32,
                        nOrigAtNosInCanonOrd: numbers,
                        ..INChI_Aux::default()
                    }])
                    .unwrap();
                let inchi = heap
                    .allocate_model_storage(vec![INChI {
                        nNumberOfAtoms: non_taut_numbers.len() as i32,
                        lenTautomer: 0,
                        ..INChI::default()
                    }])
                    .unwrap();
                sort.pINChI[TAUT_NON as usize] = inchi;
                sort.pINChI_Aux[TAUT_NON as usize] = aux;
            }
            sort
        }

        let mut heap = SourceHeap::default();
        let mut output = output_buffer(&mut heap, "pre:");
        let mut overflow = 0;
        assert_eq!(
            str_AuxNumb(
                &mut heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                1,
                0,
                0,
            ),
            Ok(4)
        );
        assert_eq!(visible(&heap, &output), "pre:");

        let ordinary = component(&mut heap, &[1, 2, 3], None);
        let ordinary_sorts = heap.allocate_model_storage(vec![ordinary]).unwrap();
        let mut output = output_buffer(&mut heap, "pre:");
        assert_eq!(
            str_AuxNumb(
                &mut heap,
                SourceMutPointer::null(),
                ordinary_sorts,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                1,
                0,
                0,
            ),
            Ok(5)
        );
        assert_eq!(visible(&heap, &output), "pre:1,2,3");

        let equal = component(&mut heap, &[1, 2, 3], Some(&[1, 2, 3]));
        let equal_sorts = heap.allocate_model_storage(vec![equal.clone()]).unwrap();
        let mut output = output_buffer(&mut heap, "pre:");
        assert_eq!(
            str_AuxNumb(
                &mut heap,
                SourceMutPointer::null(),
                equal_sorts,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                1,
                1,
                1,
            ),
            Ok(1)
        );
        assert_eq!(visible(&heap, &output), "pre:m");

        let repeated = heap.allocate_model_storage(vec![equal.clone(), equal]).unwrap();
        let mut output = output_buffer(&mut heap, "pre:");
        assert_eq!(
            str_AuxNumb(
                &mut heap,
                SourceMutPointer::null(),
                repeated,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_NT as i32,
                0,
                2,
                1,
                1,
            ),
            Ok(2)
        );
        assert_eq!(visible(&heap, &output), "pre:2m");

        let different = component(&mut heap, &[1, 2, 3], Some(&[3, 2, 1]));
        let same_again = component(&mut heap, &[1, 2, 3], Some(&[1, 2, 3]));
        let mixed = heap.allocate_model_storage(vec![same_again, different]).unwrap();
        let mut output = output_buffer(&mut heap, "pre:");
        assert_eq!(
            str_AuxNumb(
                &mut heap,
                SourceMutPointer::null(),
                mixed,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_NT as i32,
                0,
                2,
                1,
                1,
            ),
            Ok(7)
        );
        assert_eq!(visible(&heap, &output), "pre:m;3,2,1");
    }

    #[test]
    fn source_port__ichiprt3__str_auxisonumb__line_2769() {
        use crate::source_types::{INChI, INChI_Aux, OUT_NT, OUT_T1, TAUT_NON, TAUT_YES};

        fn aux(heap: &mut SourceHeap, normal: &[u16], isotopic: Option<&[u16]>) -> SourceMutPointer<INChI_Aux> {
            let number_of_atoms = i32::try_from(normal.len()).unwrap();
            let normal = heap.allocate_model_storage(normal.to_vec()).unwrap();
            let isotopic_pointer = isotopic.map_or(SourceMutPointer::null(), |values| {
                heap.allocate_model_storage(values.to_vec()).unwrap()
            });
            heap.allocate_model_storage(vec![INChI_Aux {
                nNumberOfAtoms: number_of_atoms,
                bIsIsotopic: i32::from(isotopic.is_some()),
                nOrigAtNosInCanonOrd: normal,
                nIsotopicOrigAtNosInCanonOrd: isotopic_pointer,
                ..INChI_Aux::default()
            }])
            .unwrap()
        }

        fn component(
            heap: &mut SourceHeap,
            taut_normal: &[u16],
            taut_isotopic: Option<&[u16]>,
            non_taut: Option<(&[u16], Option<&[u16]>)>,
        ) -> INCHI_SORT {
            let taut_aux = aux(heap, taut_normal, taut_isotopic);
            let taut_inchi = heap
                .allocate_model_storage(vec![INChI {
                    nNumberOfAtoms: i32::try_from(taut_normal.len()).unwrap(),
                    lenTautomer: i32::from(non_taut.is_some()),
                    ..INChI::default()
                }])
                .unwrap();
            let mut sort = INCHI_SORT::default();
            sort.pINChI[TAUT_YES as usize] = taut_inchi;
            sort.pINChI_Aux[TAUT_YES as usize] = taut_aux;
            if let Some((normal, isotopic)) = non_taut {
                let non_aux = aux(heap, normal, isotopic);
                let non_inchi = heap
                    .allocate_model_storage(vec![INChI {
                        nNumberOfAtoms: i32::try_from(normal.len()).unwrap(),
                        lenTautomer: 0,
                        ..INChI::default()
                    }])
                    .unwrap();
                sort.pINChI[TAUT_NON as usize] = non_inchi;
                sort.pINChI_Aux[TAUT_NON as usize] = non_aux;
            }
            sort
        }

        fn run(
            heap: &mut SourceHeap,
            components: Vec<INCHI_SORT>,
            output_type: i32,
            second_pass: i32,
            omit_repetitions: i32,
        ) -> (i32, String, i32) {
            let count = i32::try_from(components.len()).unwrap();
            let sorts = if components.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(components).unwrap()
            };
            let mut output = output_buffer(heap, "pre:");
            let mut overflow = 0;
            let length = str_AuxIsoNumb(
                heap,
                SourceMutPointer::null(),
                sorts,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                output_type,
                0,
                count,
                second_pass,
                omit_repetitions,
            )
            .unwrap();
            (length, visible(heap, &output), overflow)
        }

        let mut heap = SourceHeap::default();
        assert_eq!(run(&mut heap, vec![], OUT_T1 as i32, 0, 0), (0, "pre:".into(), 0));

        let concrete = component(&mut heap, &[1, 2, 3], Some(&[3, 1, 2]), None);
        assert_eq!(
            run(&mut heap, vec![concrete], OUT_T1 as i32, 0, 0),
            (5, "pre:3,1,2".into(), 0)
        );

        let missing = component(&mut heap, &[1, 2, 3], None, None);
        assert_eq!(
            run(&mut heap, vec![missing], OUT_T1 as i32, 0, 0),
            (0, "pre:".into(), 0)
        );

        let first_pass_equal = component(&mut heap, &[1, 2, 3], Some(&[1, 2, 3]), None);
        assert_eq!(
            run(&mut heap, vec![first_pass_equal], OUT_T1 as i32, 0, 1,),
            (1, "pre:m".into(), 0)
        );

        let equal_taut_normal = component(
            &mut heap,
            &[1, 2, 3],
            Some(&[3, 2, 1]),
            Some((&[3, 1, 2], Some(&[1, 2, 3]))),
        );
        assert_eq!(
            run(&mut heap, vec![equal_taut_normal], OUT_NT as i32, 1, 1,),
            (1, "pre:m".into(), 0)
        );

        let equal_non_taut_normal = component(
            &mut heap,
            &[3, 2, 1],
            Some(&[2, 1, 3]),
            Some((&[1, 2, 3], Some(&[1, 2, 3]))),
        );
        assert_eq!(
            run(&mut heap, vec![equal_non_taut_normal], OUT_NT as i32, 1, 1,),
            (1, "pre:n".into(), 0)
        );

        let equal_taut_isotopic = component(
            &mut heap,
            &[3, 2, 1],
            Some(&[1, 2, 3]),
            Some((&[2, 3, 1], Some(&[1, 2, 3]))),
        );
        assert_eq!(
            run(&mut heap, vec![equal_taut_isotopic], OUT_NT as i32, 1, 1,),
            (1, "pre:M".into(), 0)
        );

        let repeated = component(
            &mut heap,
            &[1, 2, 3],
            Some(&[3, 2, 1]),
            Some((&[3, 1, 2], Some(&[1, 2, 3]))),
        );
        let different_equivalence = component(
            &mut heap,
            &[3, 2, 1],
            Some(&[2, 1, 3]),
            Some((&[1, 2, 3], Some(&[1, 2, 3]))),
        );
        let concrete = component(
            &mut heap,
            &[1, 2, 3],
            Some(&[2, 3, 1]),
            Some((&[2, 1, 3], Some(&[3, 1, 2]))),
        );
        assert_eq!(
            run(
                &mut heap,
                vec![repeated.clone(), repeated, different_equivalence, concrete],
                OUT_NT as i32,
                1,
                1,
            ),
            (10, "pre:2m;n;3,1,2".into(), 0)
        );
    }

    #[test]
    fn source_port__ichiprt3__str_auxinvisosp3__line_3213() {
        #[derive(Clone, Copy)]
        struct StereoSpec {
            number: u16,
            parity: i8,
            inverted_number: u16,
            inverted_parity: i8,
        }

        fn stereo(heap: &mut SourceHeap, spec: StereoSpec) -> SourceMutPointer<INChI_Stereo> {
            let number = heap.allocate_model_storage(vec![spec.number]).unwrap();
            let parity = heap.allocate_model_storage(vec![spec.parity]).unwrap();
            let inverted_number = heap.allocate_model_storage(vec![spec.inverted_number]).unwrap();
            let inverted_parity = heap.allocate_model_storage(vec![spec.inverted_parity]).unwrap();
            heap.allocate_model_storage(vec![INChI_Stereo {
                nNumberOfStereoCenters: 1,
                nNumber: number,
                t_parity: parity,
                nNumberInv: inverted_number,
                t_parityInv: inverted_parity,
                nCompInv2Abs: 1,
                ..INChI_Stereo::default()
            }])
            .unwrap()
        }

        fn inchi(
            heap: &mut SourceHeap,
            normal: Option<StereoSpec>,
            isotopic: Option<StereoSpec>,
            isotope_count: i32,
            possible_isotopic_h: bool,
        ) -> SourceMutPointer<INChI> {
            let possible = if possible_isotopic_h {
                heap.allocate_model_storage(vec![2_u16]).unwrap()
            } else {
                SourceMutPointer::null()
            };
            let normal = normal.map_or(SourceMutPointer::null(), |value| stereo(heap, value));
            let isotopic = isotopic.map_or(SourceMutPointer::null(), |value| stereo(heap, value));
            heap.allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                nNumberOfIsotopicAtoms: isotope_count,
                nPossibleLocationsOfIsotopicH: possible,
                Stereo: normal,
                StereoIsotopic: isotopic,
                ..INChI::default()
            }])
            .unwrap()
        }

        fn sort(inchi: SourceMutPointer<INChI>) -> INCHI_SORT {
            INCHI_SORT {
                pINChI: [inchi, SourceMutPointer::null()],
                ..INCHI_SORT::default()
            }
        }

        fn output_buffer(heap: &mut SourceHeap) -> INCHI_IOS_STRING {
            INCHI_IOS_STRING {
                pStr: heap.allocate_model_storage(vec![0_i8; 256]).unwrap(),
                nAllocatedLength: 256,
                ..INCHI_IOS_STRING::default()
            }
        }

        fn visible(heap: &SourceHeap, output: &INCHI_IOS_STRING) -> String {
            let bytes = heap.slice(output.pStr.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            bytes[..end].iter().map(|byte| *byte as u8 as char).collect()
        }

        fn run(
            heap: &mut SourceHeap,
            current: Vec<INCHI_SORT>,
            taut: Vec<INCHI_SORT>,
            second_pass: i32,
            omit_repetitions: i32,
            use_multipliers: i32,
        ) -> (i32, String, i32) {
            let count = i32::try_from(current.len()).unwrap();
            let current = if current.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(current).unwrap()
            };
            let taut = if taut.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(taut).unwrap()
            };
            let mut output = output_buffer(heap);
            let mut overflow = 0;
            let length = str_AuxInvIsoSp3(
                heap,
                current,
                taut,
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                count,
                second_pass,
                omit_repetitions,
                use_multipliers,
            )
            .unwrap();
            (length, visible(heap, &output), overflow)
        }

        const ISO: StereoSpec = StereoSpec {
            number: 8,
            parity: 1,
            inverted_number: 2,
            inverted_parity: 2,
        };
        let mut heap = SourceHeap::default();
        assert_eq!(run(&mut heap, vec![], vec![], 0, 0, 0), (0, "".into(), 0));

        let concrete = inchi(&mut heap, None, Some(ISO), 1, false);
        assert_eq!(
            run(&mut heap, vec![sort(concrete)], vec![], 0, 0, 0),
            (2, "2+".into(), 0)
        );
        assert_eq!(
            run(&mut heap, vec![sort(concrete), sort(concrete)], vec![], 0, 0, 1,),
            (4, "2*2+".into(), 0)
        );

        let first_m = inchi(
            &mut heap,
            Some(StereoSpec {
                inverted_number: 2,
                inverted_parity: 2,
                ..ISO
            }),
            Some(ISO),
            1,
            false,
        );
        assert_eq!(run(&mut heap, vec![sort(first_m)], vec![], 0, 1, 0).1, "m");

        let first_im = inchi(
            &mut heap,
            Some(StereoSpec {
                number: 2,
                parity: 1,
                inverted_number: 3,
                inverted_parity: 2,
            }),
            Some(ISO),
            1,
            false,
        );
        assert_eq!(run(&mut heap, vec![sort(first_im)], vec![], 0, 1, 0).1, "im");

        let self_inverted = StereoSpec {
            number: 2,
            parity: 1,
            ..ISO
        };
        let first_i_m = inchi(&mut heap, None, Some(self_inverted), 0, true);
        assert_eq!(run(&mut heap, vec![sort(first_i_m)], vec![], 0, 1, 0).1, "iM");

        let second_a = inchi(&mut heap, None, Some(ISO), 1, false);
        let taut_a = inchi(
            &mut heap,
            Some(StereoSpec {
                inverted_number: 2,
                inverted_parity: 2,
                ..ISO
            }),
            None,
            0,
            false,
        );
        assert_eq!(run(&mut heap, vec![sort(second_a)], vec![sort(taut_a)], 1, 1, 0).1, "m");

        let second_b = inchi(
            &mut heap,
            Some(StereoSpec {
                inverted_number: 2,
                inverted_parity: 2,
                ..ISO
            }),
            Some(ISO),
            1,
            false,
        );
        let taut_different = inchi(
            &mut heap,
            Some(StereoSpec {
                inverted_number: 5,
                ..ISO
            }),
            None,
            0,
            false,
        );
        assert_eq!(
            run(&mut heap, vec![sort(second_b)], vec![sort(taut_different)], 1, 1, 0,).1,
            "n"
        );

        let second_c = inchi(&mut heap, None, Some(ISO), 1, false);
        let taut_c = inchi(
            &mut heap,
            Some(StereoSpec {
                inverted_number: 5,
                ..ISO
            }),
            Some(ISO),
            1,
            false,
        );
        assert_eq!(run(&mut heap, vec![sort(second_c)], vec![sort(taut_c)], 1, 1, 0).1, "M");

        let taut_e = inchi(
            &mut heap,
            Some(StereoSpec {
                inverted_number: 5,
                ..ISO
            }),
            Some(StereoSpec {
                number: 2,
                parity: 1,
                inverted_number: 6,
                inverted_parity: 2,
            }),
            1,
            false,
        );
        assert_eq!(
            run(&mut heap, vec![sort(second_c)], vec![sort(taut_e)], 1, 1, 0).1,
            "iM"
        );

        let second_f = inchi(
            &mut heap,
            Some(StereoSpec {
                number: 2,
                parity: 1,
                inverted_number: 7,
                inverted_parity: 2,
            }),
            Some(ISO),
            1,
            false,
        );
        assert_eq!(
            run(&mut heap, vec![sort(second_f)], vec![sort(taut_different)], 1, 1, 0,).1,
            "in"
        );

        let second_g = inchi(&mut heap, None, Some(self_inverted), 1, false);
        assert_eq!(
            run(&mut heap, vec![sort(second_g)], vec![sort(taut_different)], 1, 1, 0,).1,
            "iN"
        );

        let no_isotope_layer = inchi(&mut heap, None, Some(ISO), 0, false);
        assert_eq!(
            run(&mut heap, vec![sort(no_isotope_layer)], vec![], 0, 1, 0),
            (0, "".into(), 0)
        );
    }

    #[test]
    fn source_port__ichiprt3__str_auxinvisosp3numb__line_3575() {
        use crate::source_types::{OUT_NT, TAUT_NON, TAUT_YES};

        #[derive(Clone, Copy)]
        struct Numbering {
            normal: [u16; 2],
            isotopic: [u16; 2],
            inverted: [u16; 2],
            isotopic_inverted: [u16; 2],
            isotopic_flag: i32,
            normal_stereo: bool,
            isotopic_stereo: bool,
            isotopic_centers: i32,
        }

        const A: [u16; 2] = [1, 2];
        const B: [u16; 2] = [2, 1];
        const C: [u16; 2] = [1, 3];
        const D: [u16; 2] = [3, 1];
        const TARGET: [u16; 2] = [4, 2];
        const DIFFERENT: Numbering = Numbering {
            normal: A,
            isotopic: B,
            inverted: C,
            isotopic_inverted: TARGET,
            isotopic_flag: 1,
            normal_stereo: true,
            isotopic_stereo: true,
            isotopic_centers: 1,
        };

        fn stereo(heap: &mut SourceHeap, centers: i32) -> SourceMutPointer<INChI_Stereo> {
            heap.allocate_model_storage(vec![INChI_Stereo {
                nNumberOfStereoCenters: centers,
                nCompInv2Abs: 1,
                ..INChI_Stereo::default()
            }])
            .unwrap()
        }

        fn layer(
            heap: &mut SourceHeap,
            numbering: Numbering,
            tautomer_length: i32,
        ) -> (SourceMutPointer<INChI>, SourceMutPointer<INChI_Aux>) {
            let normal = heap.allocate_model_storage(numbering.normal.to_vec()).unwrap();
            let isotopic = heap.allocate_model_storage(numbering.isotopic.to_vec()).unwrap();
            let inverted = heap.allocate_model_storage(numbering.inverted.to_vec()).unwrap();
            let isotopic_inverted = heap
                .allocate_model_storage(numbering.isotopic_inverted.to_vec())
                .unwrap();
            let aux = heap
                .allocate_model_storage(vec![INChI_Aux {
                    nNumberOfAtoms: 2,
                    bIsIsotopic: numbering.isotopic_flag,
                    nOrigAtNosInCanonOrd: normal,
                    nIsotopicOrigAtNosInCanonOrd: isotopic,
                    nOrigAtNosInCanonOrdInv: inverted,
                    nIsotopicOrigAtNosInCanonOrdInv: isotopic_inverted,
                    ..INChI_Aux::default()
                }])
                .unwrap();
            let normal_stereo = numbering
                .normal_stereo
                .then(|| stereo(heap, 1))
                .unwrap_or_else(SourceMutPointer::null);
            let isotopic_stereo = numbering
                .isotopic_stereo
                .then(|| stereo(heap, numbering.isotopic_centers))
                .unwrap_or_else(SourceMutPointer::null);
            let inchi = heap
                .allocate_model_storage(vec![INChI {
                    nNumberOfAtoms: 2,
                    lenTautomer: tautomer_length,
                    Stereo: normal_stereo,
                    StereoIsotopic: isotopic_stereo,
                    ..INChI::default()
                }])
                .unwrap();
            (inchi, aux)
        }

        fn first_pass_component(heap: &mut SourceHeap, numbering: Numbering) -> INCHI_SORT {
            let (inchi, aux) = layer(heap, numbering, 0);
            let mut sort = INCHI_SORT::default();
            sort.pINChI[TAUT_YES as usize] = inchi;
            sort.pINChI_Aux[TAUT_YES as usize] = aux;
            sort
        }

        fn second_pass_component(heap: &mut SourceHeap, current: Numbering, taut: Numbering) -> INCHI_SORT {
            let (taut_inchi, taut_aux) = layer(heap, taut, 1);
            let (current_inchi, current_aux) = layer(heap, current, 0);
            let mut sort = INCHI_SORT::default();
            sort.pINChI[TAUT_YES as usize] = taut_inchi;
            sort.pINChI_Aux[TAUT_YES as usize] = taut_aux;
            sort.pINChI[TAUT_NON as usize] = current_inchi;
            sort.pINChI_Aux[TAUT_NON as usize] = current_aux;
            sort
        }

        fn run(
            heap: &mut SourceHeap,
            components: Vec<INCHI_SORT>,
            output_type: i32,
            second_pass: i32,
            omit_repetitions: i32,
        ) -> (i32, String, i32) {
            let count = i32::try_from(components.len()).unwrap();
            let components = if components.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(components).unwrap()
            };
            let mut output = output_buffer(heap, "");
            let mut overflow = 0;
            let length = str_AuxInvIsoSp3Numb(
                heap,
                SourceMutPointer::null(),
                components,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                output_type,
                0,
                count,
                second_pass,
                omit_repetitions,
            )
            .unwrap();
            (length, visible(heap, &output), overflow)
        }

        let mut heap = SourceHeap::default();
        assert_eq!(run(&mut heap, vec![], OUT_T1 as i32, 0, 0), (0, "".into(), 0));

        let concrete = first_pass_component(&mut heap, DIFFERENT);
        assert_eq!(
            run(&mut heap, vec![concrete], OUT_T1 as i32, 0, 0),
            (3, "4,2".into(), 0)
        );

        let missing = first_pass_component(
            &mut heap,
            Numbering {
                isotopic_flag: 0,
                ..DIFFERENT
            },
        );
        assert_eq!(run(&mut heap, vec![missing], OUT_T1 as i32, 0, 1), (0, "".into(), 0));

        let first_a = first_pass_component(
            &mut heap,
            Numbering {
                normal: TARGET,
                ..DIFFERENT
            },
        );
        assert_eq!(run(&mut heap, vec![first_a], OUT_T1 as i32, 0, 1).1, "m");

        let first_b = first_pass_component(
            &mut heap,
            Numbering {
                isotopic: TARGET,
                ..DIFFERENT
            },
        );
        assert_eq!(run(&mut heap, vec![first_b], OUT_T1 as i32, 0, 1).1, "M");

        let first_c = first_pass_component(
            &mut heap,
            Numbering {
                inverted: TARGET,
                ..DIFFERENT
            },
        );
        assert_eq!(run(&mut heap, vec![first_c], OUT_T1 as i32, 0, 1).1, "im");

        let current = DIFFERENT;
        let second_a = second_pass_component(
            &mut heap,
            current,
            Numbering {
                normal: TARGET,
                ..DIFFERENT
            },
        );
        assert_eq!(run(&mut heap, vec![second_a], OUT_NT as i32, 1, 1).1, "m");

        let second_b = second_pass_component(
            &mut heap,
            Numbering {
                normal: TARGET,
                ..current
            },
            DIFFERENT,
        );
        assert_eq!(run(&mut heap, vec![second_b], OUT_NT as i32, 1, 1).1, "n");

        let second_c = second_pass_component(
            &mut heap,
            current,
            Numbering {
                isotopic: TARGET,
                ..DIFFERENT
            },
        );
        assert_eq!(run(&mut heap, vec![second_c], OUT_NT as i32, 1, 1).1, "M");

        let second_d = second_pass_component(
            &mut heap,
            Numbering {
                isotopic: TARGET,
                ..current
            },
            DIFFERENT,
        );
        assert_eq!(run(&mut heap, vec![second_d], OUT_NT as i32, 1, 1).1, "N");

        let second_e = second_pass_component(
            &mut heap,
            current,
            Numbering {
                inverted: TARGET,
                ..DIFFERENT
            },
        );
        assert_eq!(run(&mut heap, vec![second_e], OUT_NT as i32, 1, 1).1, "im");

        let second_f = second_pass_component(
            &mut heap,
            Numbering {
                inverted: TARGET,
                ..current
            },
            DIFFERENT,
        );
        assert_eq!(run(&mut heap, vec![second_f], OUT_NT as i32, 1, 1).1, "in");

        let second_g = second_pass_component(
            &mut heap,
            current,
            Numbering {
                isotopic_inverted: TARGET,
                ..DIFFERENT
            },
        );
        assert_eq!(run(&mut heap, vec![second_g], OUT_NT as i32, 1, 1).1, "iM");

        let repeated_a = first_pass_component(
            &mut heap,
            Numbering {
                normal: TARGET,
                ..DIFFERENT
            },
        );
        let repeated_b = repeated_a.clone();
        assert_eq!(
            run(&mut heap, vec![repeated_a, repeated_b], OUT_T1 as i32, 0, 1,).1,
            "2m"
        );

        let transition_m = first_pass_component(
            &mut heap,
            Numbering {
                normal: TARGET,
                ..DIFFERENT
            },
        );
        let transition_im = first_pass_component(
            &mut heap,
            Numbering {
                inverted: TARGET,
                ..DIFFERENT
            },
        );
        let transition_concrete = first_pass_component(&mut heap, DIFFERENT);
        assert_eq!(
            run(
                &mut heap,
                vec![transition_m, transition_im, transition_concrete],
                OUT_T1 as i32,
                0,
                1,
            ),
            (8, "m;im;4,2".into(), 0)
        );
    }

    #[test]
    fn source_port__ichiprt3__str_auxinvsp3__line_2315() {
        fn stereo_inchi(heap: &mut SourceHeap, atoms: i32) -> SourceMutPointer<INChI> {
            let number = heap.allocate_model_storage(vec![1_u16]).unwrap();
            let parity = heap.allocate_model_storage(vec![1_i8]).unwrap();
            let number_inv = heap.allocate_model_storage(vec![1_u16]).unwrap();
            let parity_inv = heap.allocate_model_storage(vec![2_i8]).unwrap();
            let stereo = heap
                .allocate_model_storage(vec![INChI_Stereo {
                    nNumberOfStereoCenters: 1,
                    nNumber: number,
                    t_parity: parity,
                    nNumberInv: number_inv,
                    t_parityInv: parity_inv,
                    nCompInv2Abs: 1,
                    ..INChI_Stereo::default()
                }])
                .unwrap();
            heap.allocate_model_storage(vec![INChI {
                nNumberOfAtoms: atoms,
                Stereo: stereo,
                ..INChI::default()
            }])
            .unwrap()
        }

        fn output_buffer(heap: &mut SourceHeap) -> INCHI_IOS_STRING {
            INCHI_IOS_STRING {
                pStr: heap.allocate_model_storage(vec![0_i8; 128]).unwrap(),
                nAllocatedLength: 128,
                ..INCHI_IOS_STRING::default()
            }
        }

        fn visible(heap: &SourceHeap, output: &INCHI_IOS_STRING) -> String {
            let bytes = heap.slice(output.pStr.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            bytes[..end].iter().map(|byte| *byte as u8 as char).collect()
        }

        let mut heap = SourceHeap::default();
        let inchi = stereo_inchi(&mut heap, 1);
        let sorts = heap
            .allocate_model_storage(vec![INCHI_SORT {
                pINChI: [inchi, SourceMutPointer::null()],
                ..INCHI_SORT::default()
            }])
            .unwrap();
        let mut output = output_buffer(&mut heap);
        let mut overflow = 0;
        assert_eq!(
            str_AuxInvSp3(
                &mut heap,
                sorts,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                1,
                0,
                0,
                0,
            ),
            Ok(2)
        );
        assert_eq!(visible(&heap, &output), "1+");

        let mut output = output_buffer(&mut heap);
        assert_eq!(
            str_AuxInvSp3(
                &mut heap,
                sorts,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                1,
                0,
                1,
                0,
            ),
            Ok(2)
        );
        assert_eq!(visible(&heap, &output), "im");

        let repeated = heap
            .allocate_model_storage(vec![
                INCHI_SORT {
                    pINChI: [inchi, SourceMutPointer::null()],
                    ..INCHI_SORT::default()
                },
                INCHI_SORT {
                    pINChI: [inchi, SourceMutPointer::null()],
                    ..INCHI_SORT::default()
                },
            ])
            .unwrap();
        let mut output = output_buffer(&mut heap);
        assert_eq!(
            str_AuxInvSp3(
                &mut heap,
                repeated,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                2,
                0,
                1,
                1,
            ),
            Ok(3)
        );
        assert_eq!(visible(&heap, &output), "2im");

        let mut output = output_buffer(&mut heap);
        assert_eq!(
            str_AuxInvSp3(
                &mut heap,
                sorts,
                sorts,
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                1,
                1,
                1,
                0,
            ),
            Ok(1)
        );
        assert_eq!(visible(&heap, &output), "m");

        let mut empty = output_buffer(&mut heap);
        assert_eq!(
            str_AuxInvSp3(
                &mut heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut empty,
                &mut overflow,
                OUT_T1 as i32,
                0,
                0,
                0,
                0,
                0,
            ),
            Ok(0)
        );
        assert_eq!(visible(&heap, &empty), "");
    }

    #[test]
    fn source_port__ichiprt3__str_auxinvsp3numb__line_2584() {
        fn output_buffer(heap: &mut SourceHeap) -> INCHI_IOS_STRING {
            INCHI_IOS_STRING {
                pStr: heap.allocate_model_storage(vec![0_i8; 128]).unwrap(),
                nAllocatedLength: 128,
                ..INCHI_IOS_STRING::default()
            }
        }
        fn visible(heap: &SourceHeap, output: &INCHI_IOS_STRING) -> String {
            let bytes = heap.slice(output.pStr.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            bytes[..end].iter().map(|byte| *byte as u8 as char).collect()
        }
        fn component(heap: &mut SourceHeap, normal: &[u16], inverted: &[u16], with_inversion: bool) -> INCHI_SORT {
            let normal_numbers = heap.allocate_model_storage(normal.to_vec()).unwrap();
            let inverted_numbers = heap.allocate_model_storage(inverted.to_vec()).unwrap();
            let aux = heap
                .allocate_model_storage(vec![INChI_Aux {
                    nNumberOfAtoms: normal.len() as i32,
                    nOrigAtNosInCanonOrd: normal_numbers,
                    nOrigAtNosInCanonOrdInv: if with_inversion {
                        inverted_numbers
                    } else {
                        SourceMutPointer::null()
                    },
                    ..INChI_Aux::default()
                }])
                .unwrap();
            let stereo = heap
                .allocate_model_storage(vec![INChI_Stereo {
                    nNumberOfStereoCenters: if with_inversion { 1 } else { 0 },
                    nCompInv2Abs: i32::from(with_inversion),
                    ..INChI_Stereo::default()
                }])
                .unwrap();
            let inchi = heap
                .allocate_model_storage(vec![INChI {
                    nNumberOfAtoms: normal.len() as i32,
                    Stereo: stereo,
                    ..INChI::default()
                }])
                .unwrap();
            INCHI_SORT {
                pINChI: [inchi, SourceMutPointer::null()],
                pINChI_Aux: [aux, SourceMutPointer::null()],
                ..INCHI_SORT::default()
            }
        }

        let mut heap = SourceHeap::default();
        let ordinary_component = component(&mut heap, &[1, 2], &[2, 1], true);
        let ordinary = heap.allocate_model_storage(vec![ordinary_component]).unwrap();
        let mut output = output_buffer(&mut heap);
        let mut overflow = 0;
        assert_eq!(
            str_AuxInvSp3Numb(
                &mut heap,
                SourceMutPointer::null(),
                ordinary,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                1,
                0,
                0,
            ),
            Ok(3)
        );
        assert_eq!(visible(&heap, &output), "2,1");

        let equal_component = component(&mut heap, &[1, 2], &[1, 2], true);
        let equal = heap.allocate_model_storage(vec![equal_component]).unwrap();
        let mut output = output_buffer(&mut heap);
        assert_eq!(
            str_AuxInvSp3Numb(
                &mut heap,
                SourceMutPointer::null(),
                equal,
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                1,
                0,
                1,
            ),
            Ok(1)
        );
        assert_eq!(visible(&heap, &output), "m");

        let mut output = output_buffer(&mut heap);
        assert_eq!(
            str_AuxInvSp3Numb(
                &mut heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &mut output,
                &mut overflow,
                OUT_T1 as i32,
                0,
                0,
                0,
                0,
            ),
            Ok(0)
        );
        assert_eq!(visible(&heap, &output), "");
    }
}
