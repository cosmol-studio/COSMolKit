use crate::source::base::ichiprt2::{Eql_INChI_Stereo, inchi_strtol};
use crate::source::base::ichisort::{CreateNeighListFromLinearCT, FreeNeighList, insertions_sort};
use crate::source::base::util::{inchi_calloc, inchi_free};
use crate::source_types::{
    _IS_ERROR, _IS_FATAL, _IS_OKAY, _IS_WARNING, AT_NUMB, AT_RANK, EQL_SP2, EQL_SP3,
    INCHI_FLAG_RAC_STEREO, INCHI_FLAG_REL_STEREO, INCHI_MODE, INCHI_SORT, INChI, INChI_Aux,
    INChI_Stereo, INPUT_PARMS, MAX_ATOMS, SourceConstPointer, SourceHeap, SourceHeapError,
    SourceMutPointer, TAUT_NON, TAUT_YES, tagDiffINChILayers_DIFL_F, tagDiffINChILayers_DIFL_FI,
    tagDiffINChILayers_DIFL_M, tagDiffINChILayers_DIFL_MI, tagDiffINChISegments_DIFS_b_SBONDS,
    tagDiffINChISegments_DIFS_f_FORMULA, tagDiffINChISegments_DIFS_h_H_ATOMS,
    tagDiffINChISegments_DIFS_i_IATOMS, tagDiffINChISegments_DIFS_m_SP3INV,
    tagDiffINChISegments_DIFS_o_TRANSP, tagDiffINChISegments_DIFS_q_CHARGE,
    tagDiffINChISegments_DIFS_s_STYPE, tagDiffINChISegments_DIFS_t_SATOMS,
    tagMarkDiff_DIFV_BOTH_EMPTY, tagMarkDiff_DIFV_EQL2PRECED, tagMarkDiff_DIFV_FI_EQ_MI,
    tagMarkDiff_DIFV_IS_EMPTY, tagMarkDiff_DIFV_NEQ2PRECED, tagMarkDiff_DIFV_OUTPUT_OMIT_F,
};

#[allow(non_snake_case)]
pub(crate) fn GetElementAndCount(
    heap: &mut SourceHeap,
    formula: &mut SourceConstPointer<i8>,
    element: &mut [i8],
    count: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:173 GetElementAndCount
    // INCHI✔️❌: int GetElementAndCount(const char** f, char* szEl, int* count)
    // INCHI✔️❌: {
    // INCHI✔️❌:     const char* p = *f;
    // INCHI✔️❌:     char* q;
    // INCHI✔️❌:     int   i = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* djb-rwth: fixing oss-fuzz issue #37224 */
    // INCHI✔️❌:     if (p && *p)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (isupper(UCINT * p))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             szEl[i++] = *p++;
    // INCHI✔️❌:             if (*p && islower(UCINT * p))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 szEl[i++] = *p++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             szEl[i] = '\0';
    // INCHI✔️❌:             if (1 == i && szEl[0] == 'C')
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 szEl[0] = 'A'; /*  less than any element: */
    // INCHI✔️❌:                 /*  carbon-containing compounds should be first */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (*p && isdigit(UCINT * p))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *count = strtol(p, &q, 10);
    // INCHI✔️❌:                 p = q;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *count = 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             *f = p; /*  next element; */
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return -1; /*  not a chemical formula */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* v. 1.06 Changed "Zz" to "Zzz" as "Zz" is valid symbol now */
    // INCHI✔️❌:     strcpy(szEl, "Zzz");
    // INCHI✔️❌:     /*strcpy( szEl, "Zz" );*/
    // INCHI✔️❌:     /*  zero termination 'element' is larger than any other element */
    // INCHI✔️❌:     *count = 99999;         /* zero termination 'element count' is larger than any other count */
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: GetElementAndCount

    if !formula.is_null() {
        let bytes = heap.slice(*formula)?;
        let length = bytes
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        if length > 0 {
            let first = bytes[0] as u8;
            if first.is_ascii_uppercase() {
                let mut consumed = 1_usize;
                let mut symbol = [0_i8; 3];
                symbol[0] = bytes[0];
                if consumed < length && (bytes[consumed] as u8).is_ascii_lowercase() {
                    symbol[1] = bytes[consumed];
                    consumed += 1;
                }
                if consumed == 1 && symbol[0] as u8 == b'C' {
                    symbol[0] = b'A' as i8;
                }
                let has_digits = consumed < length && (bytes[consumed] as u8).is_ascii_digit();
                let symbol_length = consumed + 1;
                element
                    .get_mut(..symbol_length)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .copy_from_slice(&symbol[..symbol_length]);
                let number = (*formula).offset(
                    i64::try_from(consumed).map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                )?;
                if has_digits {
                    let mut end = SourceConstPointer::null();
                    *count = inchi_strtol(heap, number, Some(&mut end), 10)? as i32;
                    *formula = end;
                } else {
                    *count = 1;
                    *formula = number;
                }
                return Ok(1);
            }
            return Ok(-1);
        }
    }
    element
        .get_mut(..4)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .copy_from_slice(&[b'Z' as i8, b'z' as i8, b'z' as i8, 0]);
    *count = 99999;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CompareHillFormulas(
    heap: &mut SourceHeap,
    mut formula1: SourceConstPointer<i8>,
    mut formula2: SourceConstPointer<i8>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:241 CompareHillFormulas
    // INCHI✔️❌: int CompareHillFormulas(const char* f1, const char* f2)
    // INCHI✔️❌: {
    // INCHI✔️❌:     char szEl1[4], szEl2[4];
    // INCHI✔️❌:     int  count1, count2, ret1, ret2, ret;
    // INCHI✔️❌:
    // INCHI✔️❌:     do
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret1 = GetElementAndCount(&f1, szEl1, &count1);
    // INCHI✔️❌:         ret2 = GetElementAndCount(&f2, szEl2, &count2);
    // INCHI✔️❌:         if (0 <= ret1 && 0 <= ret2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((ret = strcmp(szEl1, szEl2))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret; /*  lexicographic order, string termination > any character */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = count2 - count1)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret; /*  inverse atom count order */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0; /*  program error <BRKPT> */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     } while (0 < ret1 && 0 < ret2);
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompareHillFormulas

    loop {
        let mut element1 = [0_i8; 4];
        let mut element2 = [0_i8; 4];
        let mut count1 = 0;
        let mut count2 = 0;
        let ret1 = GetElementAndCount(heap, &mut formula1, &mut element1, &mut count1)?;
        let ret2 = GetElementAndCount(heap, &mut formula2, &mut element2, &mut count2)?;
        if ret1 < 0 || ret2 < 0 {
            return Ok(0);
        }
        for index in 0..element1.len() {
            let ret = i32::from(element1[index] as u8) - i32::from(element2[index] as u8);
            if ret != 0 {
                return Ok(ret);
            }
            if element1[index] == 0 {
                break;
            }
        }
        let ret = count2.wrapping_sub(count1);
        if ret != 0 {
            return Ok(ret);
        }
        if ret1 <= 0 || ret2 <= 0 {
            return Ok(0);
        }
    }
}

#[allow(non_snake_case)]
pub(crate) fn CompareHillFormulasNoH(
    heap: &mut SourceHeap,
    mut formula1: SourceConstPointer<i8>,
    mut formula2: SourceConstPointer<i8>,
    num_h1: &mut i32,
    num_h2: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:272 CompareHillFormulasNoH
    // INCHI✔️❌: int CompareHillFormulasNoH(const char* f1,
    // INCHI✔️❌:     const char* f2,
    // INCHI✔️❌:     int* num_H1,
    // INCHI✔️❌:     int* num_H2)
    // INCHI✔️❌: {
    // INCHI✔️❌:     char szEl1[4], szEl2[4];
    // INCHI✔️❌:     int  count1, count2, ret1, ret2, ret;
    // INCHI✔️❌:
    // INCHI✔️❌:     do
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret1 = GetElementAndCount(&f1, szEl1, &count1);
    // INCHI✔️❌:         if (0 < ret1 && szEl1[0] == 'H' && !szEl1[1])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *num_H1 += count1;
    // INCHI✔️❌:             ret1 = GetElementAndCount(&f1, szEl1, &count1);
    // INCHI✔️❌:         }
    // INCHI✔️❌:         ret2 = GetElementAndCount(&f2, szEl2, &count2);
    // INCHI✔️❌:         if (0 < ret2 && szEl2[0] == 'H' && !szEl2[1])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *num_H2 += count2;
    // INCHI✔️❌:             ret2 = GetElementAndCount(&f2, szEl2, &count2);
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (0 <= ret1 && 0 <= ret2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((ret = strcmp(szEl1, szEl2))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret; /*  lexicographic order, string termination > any character */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = count2 - count1)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret; /*  inverse atom count order */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0; /*  program error <BRKPT> */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     } while (0 < ret1 && 0 < ret2);
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompareHillFormulasNoH

    loop {
        let mut element1 = [0_i8; 4];
        let mut element2 = [0_i8; 4];
        let mut count1 = 0;
        let mut count2 = 0;
        let mut ret1 = GetElementAndCount(heap, &mut formula1, &mut element1, &mut count1)?;
        if ret1 > 0 && element1[0] as u8 == b'H' && element1[1] == 0 {
            *num_h1 = num_h1
                .checked_add(count1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            ret1 = GetElementAndCount(heap, &mut formula1, &mut element1, &mut count1)?;
        }
        let mut ret2 = GetElementAndCount(heap, &mut formula2, &mut element2, &mut count2)?;
        if ret2 > 0 && element2[0] as u8 == b'H' && element2[1] == 0 {
            *num_h2 = num_h2
                .checked_add(count2)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            ret2 = GetElementAndCount(heap, &mut formula2, &mut element2, &mut count2)?;
        }
        if ret1 < 0 || ret2 < 0 {
            return Ok(0);
        }
        let mut ret = 0;
        for index in 0..4 {
            ret = i32::from(element1[index] as u8) - i32::from(element2[index] as u8);
            if ret != 0 || element1[index] == 0 {
                break;
            }
        }
        if ret != 0 {
            return Ok(ret);
        }
        let ret = count2 - count1;
        if ret != 0 {
            return Ok(ret);
        }
        if ret1 <= 0 || ret2 <= 0 {
            return Ok(0);
        }
    }
}

#[allow(non_snake_case)]
pub(crate) fn CompareTautNonIsoPartOfINChI(
    heap: &SourceHeap,
    inchi1: &INChI,
    inchi2: &INChI,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:316 CompareTautNonIsoPartOfINChI
    // INCHI✔️❌: int CompareTautNonIsoPartOfINChI(const INChI* i1, const INChI* i2)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int len1, len2, ret, i;
    // INCHI✔️❌:
    // INCHI✔️❌:     len1 = i1->lenTautomer > 0 && i1->nTautomer[0] ? i1->lenTautomer : 0;
    // INCHI✔️❌:     len2 = i2->lenTautomer > 0 && i2->nTautomer[0] ? i2->lenTautomer : 0;
    // INCHI✔️❌:     if ((ret = len2 - len1)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (i = 0; i < len1; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((ret = (int)i2->nTautomer[i] - (int)i1->nTautomer[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompareTautNonIsoPartOfINChI

    let len1 = if inchi1.lenTautomer > 0 && heap.slice(inchi1.nTautomer.as_const())?[0] != 0 {
        inchi1.lenTautomer
    } else {
        0
    };
    let len2 = if inchi2.lenTautomer > 0 && heap.slice(inchi2.nTautomer.as_const())?[0] != 0 {
        inchi2.lenTautomer
    } else {
        0
    };
    let ret = len2 - len1;
    if ret != 0 {
        return Ok(ret);
    }
    if len1 > 0 {
        let count = usize::try_from(len1).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let tautomer1 = heap
            .slice(inchi1.nTautomer.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let tautomer2 = heap
            .slice(inchi2.nTautomer.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for index in 0..count {
            let ret = i32::from(tautomer2[index]) - i32::from(tautomer1[index]);
            if ret != 0 {
                return Ok(ret);
            }
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CompINChITautVsNonTaut(
    heap: &mut SourceHeap,
    p1: &INCHI_SORT,
    p2: &INCHI_SORT,
    compare_isotopic: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:341 CompINChITautVsNonTaut
    // INCHI✔️❌: int CompINChITautVsNonTaut(const INCHI_SORT* p1,
    // INCHI✔️❌:     const INCHI_SORT* p2,
    // INCHI✔️❌:     int bCompareIsotopic)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret, num, i, num_H1, num_H2;
    // INCHI✔️❌:
    // INCHI✔️❌:     const INChI* i1 = NULL; /* Mobile-H layers in Mobile-H sorting order */
    // INCHI✔️❌:     const INChI* i2 = NULL; /* Fixed-H  layers in Fixed-H  sorting order */
    // INCHI✔️❌:
    // INCHI✔️❌:     int   n1;               /* TAUT_YES if tautomeric i1 exists, otherwise TAUT_NON */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* INChI_Stereo *Stereo1, *Stereo2; */
    // INCHI✔️❌:
    // INCHI✔️❌:     n1 = (p1->pINChI[TAUT_YES] && p1->pINChI[TAUT_YES]->nNumberOfAtoms) ? TAUT_YES : TAUT_NON;
    // INCHI✔️❌:
    // INCHI✔️❌:     i1 = p1->pINChI[n1];
    // INCHI✔️❌:     i2 = (n1 == TAUT_YES && p2->pINChI[TAUT_NON] &&
    // INCHI✔️❌:         p2->pINChI[TAUT_NON]->nNumberOfAtoms) ? p2->pINChI[TAUT_NON] : (const INChI*)NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* non-deleted-non-empty < deleted < empty */
    // INCHI✔️❌:     if (i1 && !i2)
    // INCHI✔️❌:         return 0;   /* non-empty is the smallest (first) */
    // INCHI✔️❌:     if (!i1 && i2)
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     if (!i1 && !i2)
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     if (i1->bDeleted)
    // INCHI✔️❌:         return 1;    /* deleted is the largest (last) among non-empty */
    // INCHI✔️❌:     if (i2->bDeleted)
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i1->nNumberOfAtoms > 0 && !i2->nNumberOfAtoms)
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* i2 = i2;         djb-rwth: an obviously useless statement */
    // INCHI✔️❌:
    // INCHI✔️❌:     num_H1 = num_H2 = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* do not compare terminal H */
    // INCHI✔️❌:     if ((ret = CompareHillFormulasNoH(i1->szHillFormula, i2->szHillFormula, &num_H1, &num_H2))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret;  /* lexicographic order except the shorter one is greater (last): CH2O < CH2; C3XX < C2XX */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         compare non-isotopic non-tautomeric part
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* compare number of atoms (excluding terminal H) */
    // INCHI✔️❌:     if ((ret = i2->nNumberOfAtoms - i1->nNumberOfAtoms)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         return ret; /*  more atoms first */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  compare elements  (excluding terminal H) */
    // INCHI✔️❌:     num = i1->nNumberOfAtoms;
    // INCHI✔️❌:     for (i = 0; i < num; i++)
    // INCHI✔️❌:     { /* should always be equal if Hill formulas are same */
    // INCHI✔️❌:         if ((ret = (int)i2->nAtom[i] - (int)i1->nAtom[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             return ret; /* greater periodic number first */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         compare connection tables
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if ((ret = i2->lenConnTable - i1->lenConnTable)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         return ret; /* longer connection table first */
    // INCHI✔️❌:     num = i2->lenConnTable;
    // INCHI✔️❌:     for (i = 0; i < num; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((ret = (int)i2->nConnTable[i] - (int)i1->nConnTable[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             return ret; /* greater connection table first */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*
    // INCHI✔️❌:       compare total number of H (inverse: H3 < H2 )
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if ((ret = num_H2 - num_H1)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     /*
    // INCHI✔️❌:       compare non-tautomeric num_H: N < NH3 < NH2 < NH
    // INCHI✔️❌:     */
    // INCHI✔️❌:     num = i1->nNumberOfAtoms;
    // INCHI✔️❌:     for (i = 0; i < num; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (i2->nNum_H[i] != i1->nNum_H[i])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return !i2->nNum_H[i] ? 1 :  /* no H first */
    // INCHI✔️❌:                 !i1->nNum_H[i] ? -1 :
    // INCHI✔️❌:                 (int)i2->nNum_H[i] - (int)i1->nNum_H[i];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         compare non-isotopic tautomeric part
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if ((ret = CompareTautNonIsoPartOfINChI(i1, i2))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( ret = i2->lenTautomer - i1->lenTautomer )
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     num = inchi_min( i2->lenTautomer, i1->lenTautomer );
    // INCHI✔️❌:     for ( i = 0; i < num; i ++ ) {
    // INCHI✔️❌:         if ( ret = (int)i2->nTautomer[i] - (int)i1->nTautomer[i] )
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         at this point both components are either tautomeric
    // INCHI✔️❌:         or non-tautomeric
    // INCHI✔️❌:      */
    // INCHI✔️❌:
    // INCHI✔️❌:      /*
    // INCHI✔️❌:          non-tautomeric "fixed H" specific
    // INCHI✔️❌:      */
    // INCHI✔️❌:     if ( /*TAUT_NON == bTaut && (i2 &&*/ i2->nNum_H_fixed) /* djb-rwth: fixing coverity ID #499493 */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* first, compare non-tautomeric chem. formulas -- they may be different */
    // INCHI✔️❌:         /* secondly, compare fixed-H distribution */
    // INCHI✔️❌:         if (i2->nNum_H_fixed)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num = i2->nNumberOfAtoms;
    // INCHI✔️❌:             for (i = 0; i < num; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i2->nNum_H_fixed[i] != 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         compare non-isotopic stereo
    // INCHI✔️❌:     */
    // INCHI✔️❌:     ret = CompareInchiStereo(i1->Stereo, i1->nFlags, i2->Stereo, i2->nFlags);
    // INCHI✔️❌:     if (ret)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         do not switch back to tautomeric i1, i2
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* -- how to switch back --
    // INCHI✔️❌:     if ( i1t ) {
    // INCHI✔️❌:         i1  = i1t;
    // INCHI✔️❌:         i1t = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if ( i2t ) {
    // INCHI✔️❌:         i2  = i2t;
    // INCHI✔️❌:         i2t = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:          compare isotopic non-tautomeric part
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if (bCompareIsotopic)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((ret = i2->nNumberOfIsotopicAtoms - i1->nNumberOfIsotopicAtoms)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         num = i1->nNumberOfIsotopicAtoms;
    // INCHI✔️❌:
    // INCHI✔️❌:         /*  compare isotopic atoms */
    // INCHI✔️❌:         for (i = 0; i < num; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nAtomNumber - (int)i1->IsotopicAtom[i].nAtomNumber)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nIsoDifference - (int)i1->IsotopicAtom[i].nIsoDifference)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* compare isotopic H */
    // INCHI✔️❌:         /* if tautomeric comparison mode then here are compared only non-tautomeric H */
    // INCHI✔️❌:         for (i = 0; i < num; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nNum_T - (int)i1->IsotopicAtom[i].nNum_T)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nNum_D - (int)i1->IsotopicAtom[i].nNum_D)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nNum_H - (int)i1->IsotopicAtom[i].nNum_H)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* compare isotopic tautomeric part */
    // INCHI✔️❌:         if ((ret = i2->nNumberOfIsotopicTGroups || i1->nNumberOfIsotopicTGroups)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         num = i1->nNumberOfIsotopicTGroups;
    // INCHI✔️❌:         for ( i = 0; i < num; i ++ ) {
    // INCHI✔️❌:             if ( ret = (int)i2->IsotopicTGroup[i].nTGroupNumber - (int)i1->IsotopicTGroup[i].nTGroupNumber )
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             if ( ret = (int)i2->IsotopicTGroup[i].nNum_T - (int)i1->IsotopicTGroup[i].nNum_T )
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             if ( ret = (int)i2->IsotopicTGroup[i].nNum_D - (int)i1->IsotopicTGroup[i].nNum_D )
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             if ( ret = (int)i2->IsotopicTGroup[i].nNum_H - (int)i1->IsotopicTGroup[i].nNum_H )
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         */
    // INCHI✔️❌:
    // INCHI✔️❌:         /* compare isotopic stereo */
    // INCHI✔️❌:         ret = CompareInchiStereo(i1->StereoIsotopic, i1->nFlags,
    // INCHI✔️❌:             i2->StereoIsotopic, i2->nFlags);
    // INCHI✔️❌:         if (ret)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         compare charges: non-charged first, then in order of
    // INCHI✔️❌:         ascending charges (negative first)
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i2->nTotalCharge && i1->nTotalCharge)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  both are charged; smaller charges first */
    // INCHI✔️❌:         ret = (int)i1->nTotalCharge - (int)i2->nTotalCharge;
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if ((ret = (i1->nTotalCharge ? 1 : 0) - (i2->nTotalCharge ? 1 : 0))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  only one is charged; uncharged first */
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* stable sort */
    // INCHI✔️❌:     /*ret = p1->ord_number - p2->ord_number;*/
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompINChITautVsNonTaut

    let taut_yes = TAUT_YES as usize;
    let taut_non = TAUT_NON as usize;
    let p1_taut = p1.pINChI[taut_yes];
    let n1 = if !p1_taut.is_null() && heap.slice(p1_taut.as_const())?[0].nNumberOfAtoms != 0 {
        taut_yes
    } else {
        taut_non
    };
    let i1_pointer = p1.pINChI[n1];
    let p2_non = p2.pINChI[taut_non];
    let i2_pointer = if n1 == taut_yes
        && !p2_non.is_null()
        && heap.slice(p2_non.as_const())?[0].nNumberOfAtoms != 0
    {
        p2_non
    } else {
        crate::source_types::SourceMutPointer::null()
    };
    if i1_pointer.is_null() || i2_pointer.is_null() {
        return Ok(0);
    }
    let i1 = heap.slice(i1_pointer.as_const())?[0].clone();
    let i2 = heap.slice(i2_pointer.as_const())?[0].clone();
    if i1.bDeleted != 0 {
        return Ok(1);
    }
    if i2.bDeleted != 0 {
        return Ok(-1);
    }
    if i1.nNumberOfAtoms > 0 && i2.nNumberOfAtoms == 0 {
        return Ok(0);
    }

    let mut num_h1 = 0;
    let mut num_h2 = 0;
    let ret = CompareHillFormulasNoH(
        heap,
        i1.szHillFormula.as_const(),
        i2.szHillFormula.as_const(),
        &mut num_h1,
        &mut num_h2,
    )?;
    if ret != 0 {
        return Ok(ret);
    }
    let ret = i2.nNumberOfAtoms.wrapping_sub(i1.nNumberOfAtoms);
    if ret != 0 {
        return Ok(ret);
    }
    let atom_count = usize::try_from(i1.nNumberOfAtoms.max(0))
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let atoms1 = heap
        .slice(i1.nAtom.as_const())?
        .get(..atom_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let atoms2 = heap
        .slice(i2.nAtom.as_const())?
        .get(..atom_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for index in 0..atom_count {
        let ret = i32::from(atoms2[index]) - i32::from(atoms1[index]);
        if ret != 0 {
            return Ok(ret);
        }
    }
    let ret = i2.lenConnTable.wrapping_sub(i1.lenConnTable);
    if ret != 0 {
        return Ok(ret);
    }
    let connection_count = usize::try_from(i2.lenConnTable.max(0))
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if connection_count > 0 {
        let connections1 = heap
            .slice(i1.nConnTable.as_const())?
            .get(..connection_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let connections2 = heap
            .slice(i2.nConnTable.as_const())?
            .get(..connection_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for index in 0..connection_count {
            let ret = i32::from(connections2[index]) - i32::from(connections1[index]);
            if ret != 0 {
                return Ok(ret);
            }
        }
    }
    let ret = num_h2.wrapping_sub(num_h1);
    if ret != 0 {
        return Ok(ret);
    }
    let hydrogens1 = heap
        .slice(i1.nNum_H.as_const())?
        .get(..atom_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let hydrogens2 = heap
        .slice(i2.nNum_H.as_const())?
        .get(..atom_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for index in 0..atom_count {
        if hydrogens2[index] != hydrogens1[index] {
            return Ok(if hydrogens2[index] == 0 {
                1
            } else if hydrogens1[index] == 0 {
                -1
            } else {
                i32::from(hydrogens2[index]) - i32::from(hydrogens1[index])
            });
        }
    }
    let ret = CompareTautNonIsoPartOfINChI(heap, &i1, &i2)?;
    if ret != 0 {
        return Ok(ret);
    }
    if !i2.nNum_H_fixed.is_null() {
        let fixed = heap
            .slice(i2.nNum_H_fixed.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if fixed.iter().any(|value| *value != 0) {
            return Ok(1);
        }
    }
    let stereo1 = if i1.Stereo.is_null() {
        None
    } else {
        Some(heap.slice(i1.Stereo.as_const())?[0].clone())
    };
    let stereo2 = if i2.Stereo.is_null() {
        None
    } else {
        Some(heap.slice(i2.Stereo.as_const())?[0].clone())
    };
    let ret = CompareInchiStereo(
        heap,
        stereo1.as_ref(),
        i1.nFlags,
        stereo2.as_ref(),
        i2.nFlags,
    )?;
    if ret != 0 {
        return Ok(ret);
    }

    if compare_isotopic != 0 {
        let ret = i2
            .nNumberOfIsotopicAtoms
            .wrapping_sub(i1.nNumberOfIsotopicAtoms);
        if ret != 0 {
            return Ok(ret);
        }
        let isotope_count = usize::try_from(i1.nNumberOfIsotopicAtoms.max(0))
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if isotope_count > 0 {
            let isotopes1 = heap
                .slice(i1.IsotopicAtom.as_const())?
                .get(..isotope_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let isotopes2 = heap
                .slice(i2.IsotopicAtom.as_const())?
                .get(..isotope_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            for index in 0..isotope_count {
                for ret in [
                    i32::from(isotopes2[index].nAtomNumber)
                        - i32::from(isotopes1[index].nAtomNumber),
                    i32::from(isotopes2[index].nIsoDifference)
                        - i32::from(isotopes1[index].nIsoDifference),
                ] {
                    if ret != 0 {
                        return Ok(ret);
                    }
                }
            }
            for index in 0..isotope_count {
                for ret in [
                    i32::from(isotopes2[index].nNum_T) - i32::from(isotopes1[index].nNum_T),
                    i32::from(isotopes2[index].nNum_D) - i32::from(isotopes1[index].nNum_D),
                    i32::from(isotopes2[index].nNum_H) - i32::from(isotopes1[index].nNum_H),
                ] {
                    if ret != 0 {
                        return Ok(ret);
                    }
                }
            }
        }
        if i2.nNumberOfIsotopicTGroups != 0 || i1.nNumberOfIsotopicTGroups != 0 {
            return Ok(1);
        }
        let stereo1 = if i1.StereoIsotopic.is_null() {
            None
        } else {
            Some(heap.slice(i1.StereoIsotopic.as_const())?[0].clone())
        };
        let stereo2 = if i2.StereoIsotopic.is_null() {
            None
        } else {
            Some(heap.slice(i2.StereoIsotopic.as_const())?[0].clone())
        };
        let ret = CompareInchiStereo(
            heap,
            stereo1.as_ref(),
            i1.nFlags,
            stereo2.as_ref(),
            i2.nFlags,
        )?;
        if ret != 0 {
            return Ok(ret);
        }
    }
    if i2.nTotalCharge != 0 && i1.nTotalCharge != 0 {
        return Ok(i1.nTotalCharge.wrapping_sub(i2.nTotalCharge));
    }
    Ok(i32::from(i1.nTotalCharge != 0) - i32::from(i2.nTotalCharge != 0))
}

#[allow(non_snake_case)]
pub(crate) fn GetSp3RelRacAbs(inchi: Option<&INChI>, stereo: Option<&INChI_Stereo>) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:593 GetSp3RelRacAbs
    // INCHI✔️❌: int GetSp3RelRacAbs(const INChI* pINChI, INChI_Stereo* Stereo)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int nRet = SP3_NONE;
    // INCHI✔️❌:     if (pINChI && !pINChI->bDeleted && Stereo && 0 < Stereo->nNumberOfStereoCenters)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (0 != Stereo->nCompInv2Abs)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (pINChI->nFlags & INCHI_FLAG_REL_STEREO)
    // INCHI✔️❌:             {
    // INCHI✔️❌: #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    // INCHI✔️❌:                 if (1 < Stereo->nNumberOfStereoCenters)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nRet = SP3_REL;
    // INCHI✔️❌:                 }
    // INCHI✔️❌: #else
    // INCHI✔️❌:                 nRet = SP3_REL;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:                 if (pINChI->nFlags & INCHI_FLAG_RAC_STEREO)
    // INCHI✔️❌:                 {
    // INCHI✔️❌: #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    // INCHI✔️❌:                     if (1 < Stereo->nNumberOfStereoCenters)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         nRet = SP3_REL;
    // INCHI✔️❌:                     }
    // INCHI✔️❌: #else
    // INCHI✔️❌:                     nRet = SP3_RAC;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nRet = SP3_ABS;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌: #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    // INCHI✔️❌:             if (!((pINChI->nFlags & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) && 1 == Stereo->nNumberOfStereoCenters))
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nRet = SP3_ONLY; /*  SP3_NONE if relative stereo and 1 stereocenter */
    // INCHI✔️❌:             }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return nRet;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: GetSp3RelRacAbs

    const SP3_NONE: i32 = 0;
    const SP3_ONLY: i32 = 1;
    const SP3_ABS: i32 = 2;
    const SP3_REL: i32 = 4;

    let (Some(inchi), Some(stereo)) = (inchi, stereo) else {
        return SP3_NONE;
    };
    if inchi.bDeleted != 0 || stereo.nNumberOfStereoCenters <= 0 {
        return SP3_NONE;
    }
    let relative_or_racemic = INCHI_MODE::from(INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO);
    if stereo.nCompInv2Abs != 0 {
        if inchi.nFlags & INCHI_MODE::from(INCHI_FLAG_REL_STEREO) != 0 {
            return if stereo.nNumberOfStereoCenters > 1 {
                SP3_REL
            } else {
                SP3_NONE
            };
        }
        if inchi.nFlags & INCHI_MODE::from(INCHI_FLAG_RAC_STEREO) != 0 {
            return if stereo.nNumberOfStereoCenters > 1 {
                SP3_REL
            } else {
                SP3_NONE
            };
        }
        return SP3_ABS;
    }
    if inchi.nFlags & relative_or_racemic != 0 && stereo.nNumberOfStereoCenters == 1 {
        SP3_NONE
    } else {
        SP3_ONLY
    }
}

#[allow(non_snake_case)]
pub(crate) fn CompINChILayers(
    heap: &mut SourceHeap,
    p1: &INCHI_SORT,
    p2: &INCHI_SORT,
    difference_segments: &mut [[i8; 11]; 4],
    fix_transposed_charge_bug: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:645 CompINChILayers
    // INCHI✔️❌: int CompINChILayers(const INCHI_SORT* p1,
    // INCHI✔️❌:     const INCHI_SORT* p2,
    // INCHI✔️❌:     char sDifSegs[][DIFS_LENGTH],
    // INCHI✔️❌:     int bFixTranspChargeBug)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0, num, i, num_H1, num_H2;
    // INCHI✔️❌:
    // INCHI✔️❌:     const INChI* i1 = NULL; /* Mobile-H layers in Mobile-H sorting order */
    // INCHI✔️❌:     const INChI* i2 = NULL; /* Fixed-H  layers in Fixed-H  sorting order */
    // INCHI✔️❌:
    // INCHI✔️❌:     int   n1;               /* TAUT_YES if tautomeric i1 exists, otherwise TAUT_NON */
    // INCHI✔️❌:
    // INCHI✔️❌:     INChI_Stereo* Stereo1, * Stereo2;
    // INCHI✔️❌:     INChI_Stereo* IsoStereo1, * IsoStereo2;
    // INCHI✔️❌:     int bRelRac[DIFL_LENGTH];
    // INCHI✔️❌:     char* psDifSegs;
    // INCHI✔️❌:
    // INCHI✔️❌:     n1 = (p1->pINChI[TAUT_YES] && p1->pINChI[TAUT_YES]->nNumberOfAtoms) ? TAUT_YES : TAUT_NON;
    // INCHI✔️❌:
    // INCHI✔️❌:     i1 = p1->pINChI[n1];
    // INCHI✔️❌:     i2 = (n1 == TAUT_YES && p2->pINChI[TAUT_NON] &&
    // INCHI✔️❌:         p2->pINChI[TAUT_NON]->nNumberOfAtoms) ? p2->pINChI[TAUT_NON] : (const INChI*)NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     num_H1 = num_H2 = 0;
    // INCHI✔️❌:     memset(bRelRac, DIFV_BOTH_EMPTY, sizeof(bRelRac)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /f    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     if (i1 && !i1->bDeleted && i1->szHillFormula && i1->szHillFormula[0])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_f_FORMULA] |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         if (i2 && !i2->bDeleted && i2->szHillFormula && i2->szHillFormula[0])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!CompareHillFormulasNoH(i1->szHillFormula, i2->szHillFormula, &num_H1, &num_H2) &&
    // INCHI✔️❌:                 num_H1 == num_H2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 sDifSegs[DIFL_F][DIFS_f_FORMULA] |= DIFV_EQL2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 sDifSegs[DIFL_F][DIFS_f_FORMULA] |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sDifSegs[DIFL_F][DIFS_f_FORMULA] |= i2 ? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_f_FORMULA] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         if (i2 && !i2->bDeleted && i2->szHillFormula && i2->szHillFormula[0])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sDifSegs[DIFL_F][DIFS_f_FORMULA] |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sDifSegs[DIFL_F][DIFS_f_FORMULA] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /c    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     if (i1 && !i1->bDeleted && i1->lenConnTable > 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_f_FORMULA] |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_f_FORMULA] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /h    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /* M: H atoms */
    // INCHI✔️❌:     if (i1 && !i1->bDeleted)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         num_H1 = (i1->lenTautomer > 0 && i1->nTautomer && i1->nTautomer[0]) ? 1 : 0; /* number of t-groups */
    // INCHI✔️❌:         if (!num_H1 && i1->nNum_H)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (i = 0; i < i1->nNumberOfAtoms; i++)
    // INCHI✔️❌:             { /* immobile H */
    // INCHI✔️❌:                 if (i1->nNum_H[i])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     num_H1 = 1;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_h_H_ATOMS] |= num_H1 ? DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_h_H_ATOMS] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* F: fixed mobile H */
    // INCHI✔️❌:     if (i2 && !i2->bDeleted && i2->nNum_H_fixed)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         num_H2 = 0;
    // INCHI✔️❌:         if (i1 && !i1->bDeleted)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (i = 0; i < i1->nNumberOfAtoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i2->nNum_H_fixed[i])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     num_H2 = 1;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         sDifSegs[DIFL_F][DIFS_h_H_ATOMS] |= num_H2 ? DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_F][DIFS_h_H_ATOMS] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* MI: exchangable isotopic H: see OutputINChI1(), num_iso_H[] */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /q    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_F][DIFS_q_CHARGE];
    // INCHI✔️❌:     if (i1 && !i1->bDeleted)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (i1->nTotalCharge)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sDifSegs[DIFL_M][DIFS_q_CHARGE] |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sDifSegs[DIFL_M][DIFS_q_CHARGE] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i2 && !i2->bDeleted)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i1->nTotalCharge)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i1->nTotalCharge == i2->nTotalCharge)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (i2->nTotalCharge)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         *psDifSegs |= DIFV_IS_EMPTY;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i2->nTotalCharge)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!i2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (bFixTranspChargeBug == 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* bug explanation:
    // INCHI✔️❌:
    // INCHI✔️❌:                     component #1 is tautomeric, component #2 is not
    // INCHI✔️❌:                     Mobile-H(#2) > Mobile-H(#1)
    // INCHI✔️❌:                     Fixed-H(#2) = Mobile-H(#2) < Fixed-H(#1)
    // INCHI✔️❌:
    // INCHI✔️❌:                     Layer       first_charge   second_charge
    // INCHI✔️❌:
    // INCHI✔️❌:                     Mobile-H    0    (comp#1)  -1 (comp#2)
    // INCHI✔️❌:                     Fixed-H     none (comp#2)  -1 (comp#1)
    // INCHI✔️❌:
    // INCHI✔️❌:                     v1.01 charge compared decided that charge layers are same and omitted Fixed-H /q layer
    // INCHI✔️❌:
    // INCHI✔️❌:                     Solution: when component permutation is detected AND fixed-H component does not exist,
    // INCHI✔️❌:                     compare Mobile-H charge [0 (comp#1) in the example] to the charge of Mobile-H [-1 (comp#2)]
    // INCHI✔️❌:                     of the component that has none Fixed-H charge
    // INCHI✔️❌:                     */
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* Fixed-H i2 is empty because Fixed-H struct is same as Mobile-H */
    // INCHI✔️❌:                     if (p1->ord_number != p2->ord_number && /* component order in Fixed-H is different from Mobile-H */
    // INCHI✔️❌:                         n1 == TAUT_YES && p2->pINChI[TAUT_YES] && !p2->pINChI[TAUT_YES]->bDeleted &&
    // INCHI✔️❌:                         p2->pINChI[TAUT_YES]->nNumberOfAtoms)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         int i2_nTotalCharge = p2->pINChI[TAUT_YES]->nTotalCharge;
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (i1->nTotalCharge)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (i1->nTotalCharge == i2_nTotalCharge)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 if (i2_nTotalCharge)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     *psDifSegs |= DIFV_IS_EMPTY;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (i2_nTotalCharge)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         *psDifSegs |= i1->nTotalCharge ? DIFV_EQL2PRECED : DIFV_BOTH_EMPTY;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else /* if (bFixTranspChargeBug==1) */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     *psDifSegs |= i1->nTotalCharge ? DIFV_EQL2PRECED : DIFV_BOTH_EMPTY;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else /* if ( !i2 ) { */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* i2 && i2->bDeleted */
    // INCHI✔️❌:                 *psDifSegs |= i1->nTotalCharge ? DIFV_IS_EMPTY : DIFV_BOTH_EMPTY;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_q_CHARGE] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         if (i2 && !i2->bDeleted)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i2->nTotalCharge)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 sDifSegs[DIFL_F][DIFS_q_CHARGE] |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 sDifSegs[DIFL_F][DIFS_q_CHARGE] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*************** stereo *****************/
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i1 && !i1->bDeleted)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Stereo1 = i1->Stereo;
    // INCHI✔️❌:         IsoStereo1 = i1->StereoIsotopic;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Stereo1 = NULL;
    // INCHI✔️❌:         IsoStereo1 = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (i2 && !i2->bDeleted)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Stereo2 = i2->Stereo;
    // INCHI✔️❌:         IsoStereo2 = i2->StereoIsotopic;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Stereo2 = NULL;
    // INCHI✔️❌:         IsoStereo2 = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /b    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /* M double bond stereo */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_M][DIFS_b_SBONDS];
    // INCHI✔️❌:     if (Stereo1 && Stereo1->nNumberOfStereoBonds)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* F double bond stereo */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_F][DIFS_b_SBONDS];
    // INCHI✔️❌:     if (Stereo2 && Stereo2->nNumberOfStereoBonds)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Stereo1 && Stereo1->nNumberOfStereoBonds)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (Eql_INChI_Stereo(Stereo1, EQL_SP2, Stereo2, EQL_SP2, 0))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Stereo1 && Stereo1->nNumberOfStereoBonds)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= i2 ? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* MI double bond stereo */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_MI][DIFS_b_SBONDS];
    // INCHI✔️❌:     if (IsoStereo1 && IsoStereo1->nNumberOfStereoBonds)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Eql_INChI_Stereo(IsoStereo1, EQL_SP2, Stereo1, EQL_SP2, 0))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Stereo1 && Stereo1->nNumberOfStereoBonds)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* FI double bond stereo */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_FI][DIFS_b_SBONDS];
    // INCHI✔️❌:     if (IsoStereo2 && IsoStereo2->nNumberOfStereoBonds)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Eql_INChI_Stereo(IsoStereo2, EQL_SP2, Stereo2, EQL_SP2, 0))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!(Stereo1 && Stereo1->nNumberOfStereoBonds) &&
    // INCHI✔️❌:                 !(Stereo2 && Stereo2->nNumberOfStereoBonds) &&
    // INCHI✔️❌:                 Eql_INChI_Stereo(IsoStereo2, EQL_SP2, IsoStereo1, EQL_SP2, 0))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_FI_EQ_MI;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* the solution table for FI stereo,
    // INCHI✔️❌:            in case of FI stereo is empty
    // INCHI✔️❌:            E = segment is empty, NE = not empty
    // INCHI✔️❌:            +==============================+
    // INCHI✔️❌:            | M   | MI  | F   |  result    |
    // INCHI✔️❌:            +=====+=====+=====+============+
    // INCHI✔️❌:            | E   | E   | E   | both empty |
    // INCHI✔️❌:            +-----+-----+-----+------------+
    // INCHI✔️❌:            | NE  | E   | E   | both empty |
    // INCHI✔️❌:            +-----+-----+-----+------------+
    // INCHI✔️❌:            | E   | NE  | E   | is empty   |
    // INCHI✔️❌:            +-----+-----+-----+------------+
    // INCHI✔️❌:            | NE  | NE  | E   | both empty |
    // INCHI✔️❌:            +-----+-----+-----+------------+
    // INCHI✔️❌:            | E   | E   | NE  | is empty   |
    // INCHI✔️❌:            +-----+-----+-----+------------+
    // INCHI✔️❌:            | NE  | E   | NE  | is empty   |
    // INCHI✔️❌:            +-----+-----+-----+------------+
    // INCHI✔️❌:            | E   | NE  | NE  | is empty   |
    // INCHI✔️❌:            +-----+-----+-----+------------+
    // INCHI✔️❌:            | NE  | NE  | ME  | is empty   |
    // INCHI✔️❌:            +==============================+
    // INCHI✔️❌:         */
    // INCHI✔️❌:         if (Stereo2 && Stereo2->nNumberOfStereoBonds)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:             if (IsoStereo1 && IsoStereo1->nNumberOfStereoBonds &&
    // INCHI✔️❌:                 !(Stereo1 && Stereo1->nNumberOfStereoBonds)
    // INCHI✔️❌:                 )
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= i2 ? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:             }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*==================================*/
    // INCHI✔️❌:     /*====     /t, /m, /s for M   ======*/
    // INCHI✔️❌:     /*==================================*/
    // INCHI✔️❌:     /* M sp3 stereo */
    // INCHI✔️❌:     bRelRac[DIFL_M] = GetSp3RelRacAbs(i1, Stereo1);       /* Mobile-H */
    // INCHI✔️❌:     bRelRac[DIFL_MI] = GetSp3RelRacAbs(i1, IsoStereo1);
    // INCHI✔️❌:     bRelRac[DIFL_F] = GetSp3RelRacAbs(i2, Stereo2);       /* Fixed-H */
    // INCHI✔️❌:     bRelRac[DIFL_FI] = GetSp3RelRacAbs(i2, IsoStereo2);
    // INCHI✔️❌:     if (SP3_NONE != bRelRac[DIFL_M])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_t_SATOMS] |= (bRelRac[DIFL_M] & SP3_ANY) ? DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_m_SP3INV] |= (bRelRac[DIFL_M] & SP3_ABS) ? DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_s_STYPE] |= (bRelRac[DIFL_M] & SP3_TYPE) ? DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_t_SATOMS] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_m_SP3INV] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_s_STYPE] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /t    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /* F sp3 stereo */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_F][DIFS_t_SATOMS];
    // INCHI✔️❌:     if (SP3_ANY & bRelRac[DIFL_F])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Eql_INChI_Stereo(Stereo2, EQL_SP3, Stereo1, EQL_SP3, 0))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (SP3_ANY & bRelRac[DIFL_M])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= i2 ? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* MI sp3 stereo */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_MI][DIFS_t_SATOMS];
    // INCHI✔️❌:     if (SP3_ANY & bRelRac[DIFL_MI])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Eql_INChI_Stereo(IsoStereo1, EQL_SP3, Stereo1, EQL_SP3, 0))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (SP3_ANY & bRelRac[DIFL_M])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* FI sp3 stereo */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_FI][DIFS_t_SATOMS];
    // INCHI✔️❌:     if (SP3_ANY & bRelRac[DIFL_FI])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Eql_INChI_Stereo(IsoStereo2, EQL_SP3, Stereo2, EQL_SP3, 0))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!(SP3_ANY & bRelRac[DIFL_M]) &&
    // INCHI✔️❌:                 !(SP3_ANY & bRelRac[DIFL_F]) &&
    // INCHI✔️❌:                 Eql_INChI_Stereo(IsoStereo2, EQL_SP3, IsoStereo1, EQL_SP3, 0))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_FI_EQ_MI;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else /* similar to /b */
    // INCHI✔️❌:         if ((SP3_ANY & bRelRac[DIFL_F]))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((SP3_ANY & bRelRac[DIFL_MI]) && !(SP3_ANY & bRelRac[DIFL_M]))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= i2 ? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /m    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /* F sp3 abs stereo inversion */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_F][DIFS_m_SP3INV];
    // INCHI✔️❌:     if (bRelRac[DIFL_F] & SP3_ABS)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* the order of || operands below is critically important: || is not a commutative operation */
    // INCHI✔️❌:         if (!(bRelRac[DIFL_M] & SP3_ABS) || Stereo2->nCompInv2Abs != Stereo1->nCompInv2Abs)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (bRelRac[DIFL_M] & SP3_ABS)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= i2 ? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* MI sp3 abs stereo inversion */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_MI][DIFS_m_SP3INV];
    // INCHI✔️❌:     if (SP3_ABS & bRelRac[DIFL_MI])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((SP3_ABS & bRelRac[DIFL_M]) && IsoStereo1->nCompInv2Abs == Stereo1->nCompInv2Abs)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (SP3_ABS & bRelRac[DIFL_M])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* FI sp3 abs stereo inversion */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_FI][DIFS_m_SP3INV];
    // INCHI✔️❌:     if (SP3_ABS & bRelRac[DIFL_FI])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((SP3_ABS & bRelRac[DIFL_F]) && IsoStereo2->nCompInv2Abs == Stereo2->nCompInv2Abs)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!(SP3_ABS & bRelRac[DIFL_M]) &&
    // INCHI✔️❌:                 !(SP3_ABS & bRelRac[DIFL_F]) &&
    // INCHI✔️❌:                 (SP3_ABS & bRelRac[DIFL_MI]) && /* make sure IsoStereo1 != NULL */
    // INCHI✔️❌:                 IsoStereo2->nCompInv2Abs == IsoStereo1->nCompInv2Abs)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_FI_EQ_MI;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* similar to /b */
    // INCHI✔️❌:         /* the order of || operands below is critically important: || is no a commutative operation */
    // INCHI✔️❌:         if ((SP3_ABS & bRelRac[DIFL_F]))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((SP3_ABS & bRelRac[DIFL_MI]) && !(SP3_ABS & bRelRac[DIFL_M]))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= i2 ? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /s    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /* F sp3 stereo type */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_F][DIFS_s_STYPE];
    // INCHI✔️❌:     if (bRelRac[DIFL_F] & SP3_TYPE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((bRelRac[DIFL_F] & SP3_TYPE) == (bRelRac[DIFL_M] & SP3_TYPE))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (bRelRac[DIFL_M] & SP3_TYPE)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= i2 ? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* MI sp3 stereo type */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_MI][DIFS_s_STYPE];
    // INCHI✔️❌:     if (SP3_TYPE & bRelRac[DIFL_MI])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((SP3_TYPE & bRelRac[DIFL_MI]) == (SP3_TYPE & bRelRac[DIFL_M]))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (SP3_TYPE & bRelRac[DIFL_M])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* FI sp3 stereo type */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_FI][DIFS_s_STYPE];
    // INCHI✔️❌:     if (SP3_TYPE & bRelRac[DIFL_FI])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((SP3_TYPE & bRelRac[DIFL_FI]) == (SP3_TYPE & bRelRac[DIFL_F]))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!(SP3_TYPE & bRelRac[DIFL_M]) &&
    // INCHI✔️❌:                 !(SP3_TYPE & bRelRac[DIFL_F]) &&
    // INCHI✔️❌:                 (SP3_TYPE & bRelRac[DIFL_MI]))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_FI_EQ_MI;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* similar to /b */
    // INCHI✔️❌:         /* the order of || operands below is critically important: || is not a commutative operation */
    // INCHI✔️❌:         if ((SP3_TYPE & bRelRac[DIFL_F]))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((SP3_TYPE & bRelRac[DIFL_MI]) && !(SP3_TYPE & bRelRac[DIFL_M]))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= i2 ? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /o    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     if (p1 && p2 && p1->ord_number != p2->ord_number)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_F][DIFS_o_TRANSP] |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /i    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:
    // INCHI✔️❌:     /* M isotopic atoms */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_MI][DIFS_i_IATOMS];
    // INCHI✔️❌:     if (i1 && !i1->bDeleted && (i1->nNumberOfIsotopicAtoms || i1->nNumberOfIsotopicTGroups))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* F isotopic atoms */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_FI][DIFS_i_IATOMS];
    // INCHI✔️❌:     if (i2 && !i2->bDeleted)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (i2->nNumberOfIsotopicAtoms || i2->nNumberOfIsotopicTGroups)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!i1 || i1->bDeleted ||
    // INCHI✔️❌:                 i2->nNumberOfIsotopicAtoms != i1->nNumberOfIsotopicAtoms ||
    // INCHI✔️❌:                 i2->nNumberOfIsotopicTGroups != i1->nNumberOfIsotopicTGroups)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 int diff;
    // INCHI✔️❌:                 num = i1->nNumberOfIsotopicAtoms;
    // INCHI✔️❌:                 diff = 0;
    // INCHI✔️❌:                 for (i = 0; i < num; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* compare isotopic atoms */
    // INCHI✔️❌:                     if ((diff = (int)i2->IsotopicAtom[i].nAtomNumber - (int)i1->IsotopicAtom[i].nAtomNumber)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if ((diff = (int)i2->IsotopicAtom[i].nIsoDifference - (int)i1->IsotopicAtom[i].nIsoDifference)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /* compare isotopic H */
    // INCHI✔️❌:                     if ((diff = (int)i2->IsotopicAtom[i].nNum_T - (int)i1->IsotopicAtom[i].nNum_T)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if ((diff = (int)i2->IsotopicAtom[i].nNum_D - (int)i1->IsotopicAtom[i].nNum_D)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if ((diff = (int)i2->IsotopicAtom[i].nNum_H - (int)i1->IsotopicAtom[i].nNum_H)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (!diff)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     num = i1->nNumberOfIsotopicTGroups;
    // INCHI✔️❌:                     for (i = 0; i < num; i++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if ((diff = (int)i2->IsotopicTGroup[i].nTGroupNumber - (int)i1->IsotopicTGroup[i].nTGroupNumber)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if ((diff = (int)i2->IsotopicTGroup[i].nNum_T - (int)i1->IsotopicTGroup[i].nNum_T)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if ((diff = (int)i2->IsotopicTGroup[i].nNum_D - (int)i1->IsotopicTGroup[i].nNum_D)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             return diff;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if ((diff = (int)i2->IsotopicTGroup[i].nNum_H - (int)i1->IsotopicTGroup[i].nNum_H)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 *psDifSegs |= diff ? DIFV_NEQ2PRECED : DIFV_FI_EQ_MI;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i1 && !i1->bDeleted && (i1->nNumberOfIsotopicAtoms || i1->nNumberOfIsotopicTGroups))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_IS_EMPTY;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (!i2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i1 && !i1->bDeleted && (i1->nNumberOfIsotopicAtoms || i1->nNumberOfIsotopicTGroups))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // INCHI✔️❌:
    // END INCHI C FUNCTION: CompINChILayers

    const M: usize = tagDiffINChILayers_DIFL_M as usize;
    const MI: usize = tagDiffINChILayers_DIFL_MI as usize;
    const F: usize = tagDiffINChILayers_DIFL_F as usize;
    const FI: usize = tagDiffINChILayers_DIFL_FI as usize;
    const FORMULA: usize = tagDiffINChISegments_DIFS_f_FORMULA as usize;
    const H_ATOMS: usize = tagDiffINChISegments_DIFS_h_H_ATOMS as usize;
    const CHARGE: usize = tagDiffINChISegments_DIFS_q_CHARGE as usize;
    const SBONDS: usize = tagDiffINChISegments_DIFS_b_SBONDS as usize;
    const SATOMS: usize = tagDiffINChISegments_DIFS_t_SATOMS as usize;
    const SP3INV: usize = tagDiffINChISegments_DIFS_m_SP3INV as usize;
    const STYPE: usize = tagDiffINChISegments_DIFS_s_STYPE as usize;
    const IATOMS: usize = tagDiffINChISegments_DIFS_i_IATOMS as usize;
    const TRANSP: usize = tagDiffINChISegments_DIFS_o_TRANSP as usize;
    const BOTH_EMPTY: i8 = tagMarkDiff_DIFV_BOTH_EMPTY as i8;
    const EQUAL: i8 = tagMarkDiff_DIFV_EQL2PRECED as i8;
    const NOT_EQUAL: i8 = tagMarkDiff_DIFV_NEQ2PRECED as i8;
    const IS_EMPTY: i8 = tagMarkDiff_DIFV_IS_EMPTY as i8;
    const FI_EQ_MI: i8 = tagMarkDiff_DIFV_FI_EQ_MI as i8;
    const SP3_NONE: i32 = 0;
    const SP3_ONLY: i32 = 1;
    const SP3_ABS: i32 = 2;
    const SP3_REL: i32 = 4;
    const SP3_RAC: i32 = 8;
    const SP3_TYPE: i32 = SP3_ABS | SP3_REL | SP3_RAC;
    const SP3_ANY: i32 = SP3_TYPE | SP3_ONLY;

    let taut_yes = TAUT_YES as usize;
    let taut_non = TAUT_NON as usize;
    let p1_taut = p1.pINChI[taut_yes];
    let n1 = if !p1_taut.is_null() && heap.slice(p1_taut.as_const())?[0].nNumberOfAtoms != 0 {
        taut_yes
    } else {
        taut_non
    };
    let i1_pointer = p1.pINChI[n1];
    let p2_non = p2.pINChI[taut_non];
    let i2_pointer = if n1 == taut_yes
        && !p2_non.is_null()
        && heap.slice(p2_non.as_const())?[0].nNumberOfAtoms != 0
    {
        p2_non
    } else {
        Default::default()
    };
    let i1 = if i1_pointer.is_null() {
        None
    } else {
        Some(heap.slice(i1_pointer.as_const())?[0].clone())
    };
    let i2 = if i2_pointer.is_null() {
        None
    } else {
        Some(heap.slice(i2_pointer.as_const())?[0].clone())
    };
    let i1_live = i1.as_ref().filter(|inchi| inchi.bDeleted == 0);
    let i2_live = i2.as_ref().filter(|inchi| inchi.bDeleted == 0);

    let formula1_present = if let Some(inchi) = i1_live {
        !inchi.szHillFormula.is_null() && heap.slice(inchi.szHillFormula.as_const())?[0] != 0
    } else {
        false
    };
    let formula2_present = if let Some(inchi) = i2_live {
        !inchi.szHillFormula.is_null() && heap.slice(inchi.szHillFormula.as_const())?[0] != 0
    } else {
        false
    };
    let mut num_h1 = 0;
    let mut num_h2 = 0;
    if formula1_present {
        difference_segments[M][FORMULA] |= NOT_EQUAL;
        if formula2_present {
            let first = i1_live.ok_or(SourceHeapError::NullPointer)?;
            let second = i2_live.ok_or(SourceHeapError::NullPointer)?;
            if CompareHillFormulasNoH(
                heap,
                first.szHillFormula.as_const(),
                second.szHillFormula.as_const(),
                &mut num_h1,
                &mut num_h2,
            )? == 0
                && num_h1 == num_h2
            {
                difference_segments[F][FORMULA] |= EQUAL;
            } else {
                difference_segments[F][FORMULA] |= NOT_EQUAL;
            }
        } else {
            difference_segments[F][FORMULA] |= if i2.is_some() { IS_EMPTY } else { EQUAL };
        }
    } else {
        difference_segments[M][FORMULA] |= BOTH_EMPTY;
        difference_segments[F][FORMULA] |= if formula2_present {
            NOT_EQUAL
        } else {
            BOTH_EMPTY
        };
    }

    difference_segments[M][FORMULA] |= if i1_live.is_some_and(|inchi| inchi.lenConnTable > 1) {
        NOT_EQUAL
    } else {
        BOTH_EMPTY
    };

    let mobile_h_present = if let Some(inchi) = i1_live {
        let tautomeric = inchi.lenTautomer > 0
            && !inchi.nTautomer.is_null()
            && heap.slice(inchi.nTautomer.as_const())?[0] != 0;
        if tautomeric {
            true
        } else if !inchi.nNum_H.is_null() {
            let count = usize::try_from(inchi.nNumberOfAtoms.max(0))
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            heap.slice(inchi.nNum_H.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .iter()
                .any(|value| *value != 0)
        } else {
            false
        }
    } else {
        false
    };
    difference_segments[M][H_ATOMS] |= if mobile_h_present {
        NOT_EQUAL
    } else {
        BOTH_EMPTY
    };
    let fixed_h_present = if let (Some(first), Some(second)) = (i1_live, i2_live) {
        if second.nNum_H_fixed.is_null() {
            false
        } else {
            let count = usize::try_from(first.nNumberOfAtoms.max(0))
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            heap.slice(second.nNum_H_fixed.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .iter()
                .any(|value| *value != 0)
        }
    } else {
        false
    };
    difference_segments[F][H_ATOMS] |= if fixed_h_present {
        NOT_EQUAL
    } else {
        BOTH_EMPTY
    };

    if let Some(first) = i1_live {
        difference_segments[M][CHARGE] |= if first.nTotalCharge != 0 {
            NOT_EQUAL
        } else {
            BOTH_EMPTY
        };
        if let Some(second) = i2_live {
            difference_segments[F][CHARGE] |= if first.nTotalCharge != 0 {
                if first.nTotalCharge == second.nTotalCharge {
                    EQUAL
                } else if second.nTotalCharge != 0 {
                    NOT_EQUAL
                } else {
                    IS_EMPTY
                }
            } else if second.nTotalCharge != 0 {
                NOT_EQUAL
            } else {
                BOTH_EMPTY
            };
        } else if i2.is_none() {
            let fixed_mark = if fix_transposed_charge_bug == 1
                && p1.ord_number != p2.ord_number
                && n1 == taut_yes
            {
                let p2_taut = p2.pINChI[taut_yes];
                if !p2_taut.is_null() {
                    let second_taut = &heap.slice(p2_taut.as_const())?[0];
                    if second_taut.bDeleted == 0 && second_taut.nNumberOfAtoms != 0 {
                        if first.nTotalCharge != 0 {
                            if first.nTotalCharge == second_taut.nTotalCharge {
                                EQUAL
                            } else if second_taut.nTotalCharge != 0 {
                                NOT_EQUAL
                            } else {
                                IS_EMPTY
                            }
                        } else if second_taut.nTotalCharge != 0 {
                            NOT_EQUAL
                        } else {
                            BOTH_EMPTY
                        }
                    } else if first.nTotalCharge != 0 {
                        EQUAL
                    } else {
                        BOTH_EMPTY
                    }
                } else if first.nTotalCharge != 0 {
                    EQUAL
                } else {
                    BOTH_EMPTY
                }
            } else if first.nTotalCharge != 0 {
                EQUAL
            } else {
                BOTH_EMPTY
            };
            difference_segments[F][CHARGE] |= fixed_mark;
        } else {
            difference_segments[F][CHARGE] |= if first.nTotalCharge != 0 {
                IS_EMPTY
            } else {
                BOTH_EMPTY
            };
        }
    } else {
        difference_segments[M][CHARGE] |= BOTH_EMPTY;
        if let Some(second) = i2_live {
            difference_segments[F][CHARGE] |= if second.nTotalCharge != 0 {
                NOT_EQUAL
            } else {
                BOTH_EMPTY
            };
        }
    }

    let stereo1 = if let Some(inchi) = i1_live.filter(|inchi| !inchi.Stereo.is_null()) {
        Some(heap.slice(inchi.Stereo.as_const())?[0].clone())
    } else {
        None
    };
    let iso_stereo1 = if let Some(inchi) = i1_live.filter(|inchi| !inchi.StereoIsotopic.is_null()) {
        Some(heap.slice(inchi.StereoIsotopic.as_const())?[0].clone())
    } else {
        None
    };
    let stereo2 = if let Some(inchi) = i2_live.filter(|inchi| !inchi.Stereo.is_null()) {
        Some(heap.slice(inchi.Stereo.as_const())?[0].clone())
    } else {
        None
    };
    let iso_stereo2 = if let Some(inchi) = i2_live.filter(|inchi| !inchi.StereoIsotopic.is_null()) {
        Some(heap.slice(inchi.StereoIsotopic.as_const())?[0].clone())
    } else {
        None
    };
    let has_bonds = |stereo: Option<&INChI_Stereo>| {
        stereo.is_some_and(|stereo| stereo.nNumberOfStereoBonds != 0)
    };
    let s1 = stereo1.as_ref();
    let is1 = iso_stereo1.as_ref();
    let s2 = stereo2.as_ref();
    let is2 = iso_stereo2.as_ref();

    difference_segments[M][SBONDS] |= if has_bonds(s1) { NOT_EQUAL } else { BOTH_EMPTY };
    difference_segments[F][SBONDS] |= if has_bonds(s2) {
        if has_bonds(s1) {
            if Eql_INChI_Stereo(heap, s1, EQL_SP2 as i32, s2, EQL_SP2 as i32, 0)? != 0 {
                EQUAL
            } else {
                NOT_EQUAL
            }
        } else {
            NOT_EQUAL
        }
    } else if has_bonds(s1) {
        if i2.is_some() { IS_EMPTY } else { EQUAL }
    } else {
        BOTH_EMPTY
    };
    difference_segments[MI][SBONDS] |= if has_bonds(is1) {
        if Eql_INChI_Stereo(heap, is1, EQL_SP2 as i32, s1, EQL_SP2 as i32, 0)? != 0 {
            EQUAL
        } else {
            NOT_EQUAL
        }
    } else if has_bonds(s1) {
        EQUAL
    } else {
        BOTH_EMPTY
    };
    difference_segments[FI][SBONDS] |= if has_bonds(is2) {
        if Eql_INChI_Stereo(heap, is2, EQL_SP2 as i32, s2, EQL_SP2 as i32, 0)? != 0 {
            EQUAL
        } else if !has_bonds(s1)
            && !has_bonds(s2)
            && Eql_INChI_Stereo(heap, is2, EQL_SP2 as i32, is1, EQL_SP2 as i32, 0)? != 0
        {
            FI_EQ_MI
        } else {
            NOT_EQUAL
        }
    } else if has_bonds(s2) {
        EQUAL
    } else if has_bonds(is1) && !has_bonds(s1) {
        if i2.is_some() { IS_EMPTY } else { EQUAL }
    } else {
        BOTH_EMPTY
    };

    let rel_rac = [
        GetSp3RelRacAbs(i1_live, s1),
        GetSp3RelRacAbs(i1_live, is1),
        GetSp3RelRacAbs(i2_live, s2),
        GetSp3RelRacAbs(i2_live, is2),
    ];
    if rel_rac[M] != SP3_NONE {
        difference_segments[M][SATOMS] |= if rel_rac[M] & SP3_ANY != 0 {
            NOT_EQUAL
        } else {
            BOTH_EMPTY
        };
        difference_segments[M][SP3INV] |= if rel_rac[M] & SP3_ABS != 0 {
            NOT_EQUAL
        } else {
            BOTH_EMPTY
        };
        difference_segments[M][STYPE] |= if rel_rac[M] & SP3_TYPE != 0 {
            NOT_EQUAL
        } else {
            BOTH_EMPTY
        };
    } else {
        difference_segments[M][SATOMS] |= BOTH_EMPTY;
        difference_segments[M][SP3INV] |= BOTH_EMPTY;
        difference_segments[M][STYPE] |= BOTH_EMPTY;
    }

    difference_segments[F][SATOMS] |= if rel_rac[F] & SP3_ANY != 0 {
        if Eql_INChI_Stereo(heap, s2, EQL_SP3 as i32, s1, EQL_SP3 as i32, 0)? != 0 {
            EQUAL
        } else {
            NOT_EQUAL
        }
    } else if rel_rac[M] & SP3_ANY != 0 {
        if i2.is_some() { IS_EMPTY } else { EQUAL }
    } else {
        BOTH_EMPTY
    };
    difference_segments[MI][SATOMS] |= if rel_rac[MI] & SP3_ANY != 0 {
        if Eql_INChI_Stereo(heap, is1, EQL_SP3 as i32, s1, EQL_SP3 as i32, 0)? != 0 {
            EQUAL
        } else {
            NOT_EQUAL
        }
    } else if rel_rac[M] & SP3_ANY != 0 {
        EQUAL
    } else {
        BOTH_EMPTY
    };
    difference_segments[FI][SATOMS] |= if rel_rac[FI] & SP3_ANY != 0 {
        if Eql_INChI_Stereo(heap, is2, EQL_SP3 as i32, s2, EQL_SP3 as i32, 0)? != 0 {
            EQUAL
        } else if rel_rac[M] & SP3_ANY == 0
            && rel_rac[F] & SP3_ANY == 0
            && Eql_INChI_Stereo(heap, is2, EQL_SP3 as i32, is1, EQL_SP3 as i32, 0)? != 0
        {
            FI_EQ_MI
        } else {
            NOT_EQUAL
        }
    } else if rel_rac[F] & SP3_ANY != 0 {
        EQUAL
    } else if rel_rac[MI] & SP3_ANY != 0 && rel_rac[M] & SP3_ANY == 0 {
        if i2.is_some() { IS_EMPTY } else { EQUAL }
    } else {
        BOTH_EMPTY
    };

    difference_segments[F][SP3INV] |= if rel_rac[F] & SP3_ABS != 0 {
        if rel_rac[M] & SP3_ABS == 0
            || s2.ok_or(SourceHeapError::NullPointer)?.nCompInv2Abs
                != s1.ok_or(SourceHeapError::NullPointer)?.nCompInv2Abs
        {
            NOT_EQUAL
        } else {
            EQUAL
        }
    } else if rel_rac[M] & SP3_ABS != 0 {
        if i2.is_some() { IS_EMPTY } else { EQUAL }
    } else {
        BOTH_EMPTY
    };
    difference_segments[MI][SP3INV] |= if rel_rac[MI] & SP3_ABS != 0 {
        if rel_rac[M] & SP3_ABS != 0
            && is1.ok_or(SourceHeapError::NullPointer)?.nCompInv2Abs
                == s1.ok_or(SourceHeapError::NullPointer)?.nCompInv2Abs
        {
            EQUAL
        } else {
            NOT_EQUAL
        }
    } else if rel_rac[M] & SP3_ABS != 0 {
        EQUAL
    } else {
        BOTH_EMPTY
    };
    difference_segments[FI][SP3INV] |= if rel_rac[FI] & SP3_ABS != 0 {
        if rel_rac[F] & SP3_ABS != 0
            && is2.ok_or(SourceHeapError::NullPointer)?.nCompInv2Abs
                == s2.ok_or(SourceHeapError::NullPointer)?.nCompInv2Abs
        {
            EQUAL
        } else if rel_rac[M] & SP3_ABS == 0
            && rel_rac[F] & SP3_ABS == 0
            && rel_rac[MI] & SP3_ABS != 0
            && is2.ok_or(SourceHeapError::NullPointer)?.nCompInv2Abs
                == is1.ok_or(SourceHeapError::NullPointer)?.nCompInv2Abs
        {
            FI_EQ_MI
        } else {
            NOT_EQUAL
        }
    } else if rel_rac[F] & SP3_ABS != 0 {
        EQUAL
    } else if rel_rac[MI] & SP3_ABS != 0 && rel_rac[M] & SP3_ABS == 0 {
        if i2.is_some() { IS_EMPTY } else { EQUAL }
    } else {
        BOTH_EMPTY
    };

    difference_segments[F][STYPE] |= if rel_rac[F] & SP3_TYPE != 0 {
        if rel_rac[F] & SP3_TYPE == rel_rac[M] & SP3_TYPE {
            EQUAL
        } else {
            NOT_EQUAL
        }
    } else if rel_rac[M] & SP3_TYPE != 0 {
        if i2.is_some() { IS_EMPTY } else { EQUAL }
    } else {
        BOTH_EMPTY
    };
    difference_segments[MI][STYPE] |= if rel_rac[MI] & SP3_TYPE != 0 {
        if rel_rac[MI] & SP3_TYPE == rel_rac[M] & SP3_TYPE {
            EQUAL
        } else {
            NOT_EQUAL
        }
    } else if rel_rac[M] & SP3_TYPE != 0 {
        EQUAL
    } else {
        BOTH_EMPTY
    };
    difference_segments[FI][STYPE] |= if rel_rac[FI] & SP3_TYPE != 0 {
        if rel_rac[FI] & SP3_TYPE == rel_rac[F] & SP3_TYPE {
            EQUAL
        } else if rel_rac[M] & SP3_TYPE == 0
            && rel_rac[F] & SP3_TYPE == 0
            && rel_rac[MI] & SP3_TYPE != 0
        {
            FI_EQ_MI
        } else {
            NOT_EQUAL
        }
    } else if rel_rac[F] & SP3_TYPE != 0 {
        EQUAL
    } else if rel_rac[MI] & SP3_TYPE != 0 && rel_rac[M] & SP3_TYPE == 0 {
        if i2.is_some() { IS_EMPTY } else { EQUAL }
    } else {
        BOTH_EMPTY
    };

    if p1.ord_number != p2.ord_number {
        difference_segments[F][TRANSP] |= NOT_EQUAL;
    }

    let isotopic1_present = i1_live.is_some_and(|inchi| {
        inchi.nNumberOfIsotopicAtoms != 0 || inchi.nNumberOfIsotopicTGroups != 0
    });
    difference_segments[MI][IATOMS] |= if isotopic1_present {
        NOT_EQUAL
    } else {
        BOTH_EMPTY
    };
    if let Some(second) = i2_live {
        let isotopic2_present =
            second.nNumberOfIsotopicAtoms != 0 || second.nNumberOfIsotopicTGroups != 0;
        if isotopic2_present {
            let Some(first) = i1_live else {
                difference_segments[FI][IATOMS] |= NOT_EQUAL;
                return Ok(0);
            };
            if second.nNumberOfIsotopicAtoms != first.nNumberOfIsotopicAtoms
                || second.nNumberOfIsotopicTGroups != first.nNumberOfIsotopicTGroups
            {
                difference_segments[FI][IATOMS] |= NOT_EQUAL;
            } else {
                let atom_count = usize::try_from(first.nNumberOfIsotopicAtoms.max(0))
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let mut different = false;
                if atom_count > 0 {
                    let atoms1 = heap
                        .slice(first.IsotopicAtom.as_const())?
                        .get(..atom_count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let atoms2 = heap
                        .slice(second.IsotopicAtom.as_const())?
                        .get(..atom_count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    different = atoms1 != atoms2;
                }
                if !different {
                    let group_count = usize::try_from(first.nNumberOfIsotopicTGroups.max(0))
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    if group_count > 0 {
                        let groups1 = heap
                            .slice(first.IsotopicTGroup.as_const())?
                            .get(..group_count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let groups2 = heap
                            .slice(second.IsotopicTGroup.as_const())?
                            .get(..group_count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        for index in 0..group_count {
                            if groups2[index].nTGroupNumber != groups1[index].nTGroupNumber
                                || groups2[index].nNum_T != groups1[index].nNum_T
                            {
                                different = true;
                                break;
                            }
                            let d_difference =
                                i32::from(groups2[index].nNum_D) - i32::from(groups1[index].nNum_D);
                            if d_difference != 0 {
                                return Ok(d_difference);
                            }
                            if groups2[index].nNum_H != groups1[index].nNum_H {
                                different = true;
                                break;
                            }
                        }
                    }
                }
                difference_segments[FI][IATOMS] |= if different { NOT_EQUAL } else { FI_EQ_MI };
            }
        } else if isotopic1_present {
            difference_segments[FI][IATOMS] |= IS_EMPTY;
        }
    } else if i2.is_none() {
        difference_segments[FI][IATOMS] |= if isotopic1_present { EQUAL } else { BOTH_EMPTY };
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CompareInchiStereo(
    heap: &SourceHeap,
    stereo1: Option<&INChI_Stereo>,
    flags1: INCHI_MODE,
    stereo2: Option<&INChI_Stereo>,
    flags2: INCHI_MODE,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:1607 CompareInchiStereo
    // INCHI✔️❌: int CompareInchiStereo(INChI_Stereo* Stereo1,
    // INCHI✔️❌:     INCHI_MODE nFlags1,
    // INCHI✔️❌:     INChI_Stereo* Stereo2,
    // INCHI✔️❌:     INCHI_MODE nFlags2)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, num, ret;
    // INCHI✔️❌:     if (Stereo2 && Stereo1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  compare stereogenic bonds */
    // INCHI✔️❌:         num = inchi_min(Stereo2->nNumberOfStereoBonds, Stereo1->nNumberOfStereoBonds);
    // INCHI✔️❌:         for (i = 0; i < num; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((ret = (int)Stereo2->nBondAtom1[i] - (int)Stereo1->nBondAtom1[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)Stereo2->nBondAtom2[i] - (int)Stereo1->nBondAtom2[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)Stereo2->b_parity[i] - (int)Stereo1->b_parity[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if ((ret = (int)Stereo2->nNumberOfStereoBonds - (int)Stereo1->nNumberOfStereoBonds)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /*  compare stereogenic atoms */
    // INCHI✔️❌: #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    // INCHI✔️❌:         if (((nFlags1 | nFlags2) & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) &&
    // INCHI✔️❌:             1 == Stereo2->nNumberOfStereoCenters &&
    // INCHI✔️❌:             1 == Stereo1->nNumberOfStereoCenters)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ; /*  do not compare single stereocenters in case of relative stereo */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num = inchi_min(Stereo2->nNumberOfStereoCenters, Stereo1->nNumberOfStereoCenters);
    // INCHI✔️❌:             for (i = 0; i < num; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if ((ret = (int)Stereo2->nNumber[i] - (int)Stereo1->nNumber[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return ret;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if ((ret = (int)Stereo2->t_parity[i] - (int)Stereo1->t_parity[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return ret;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)Stereo2->nNumberOfStereoCenters - (int)Stereo1->nNumberOfStereoCenters)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             /*  compare stereo-abs-is-inverted flags  for non-relative, non-racemic */
    // INCHI✔️❌:             if (!((nFlags1 | nFlags2) & (INCHI_FLAG_RAC_STEREO | INCHI_FLAG_REL_STEREO)))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if ((ret = (Stereo2->nCompInv2Abs < 0) - (Stereo1->nCompInv2Abs < 0))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return ret;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Stereo2 && (Stereo2->nNumberOfStereoBonds > 0 || Stereo2->nNumberOfStereoCenters > 0
    // INCHI✔️❌: #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    // INCHI✔️❌:             && /*  do not compare single stereocenters in case of relative stereo */
    // INCHI✔️❌:             !((nFlags2 & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) &&
    // INCHI✔️❌:                 1 == Stereo2->nNumberOfStereoCenters
    // INCHI✔️❌:                 )
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (Stereo1 && (Stereo1->nNumberOfStereoBonds > 0 || Stereo1->nNumberOfStereoCenters > 0
    // INCHI✔️❌: #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    // INCHI✔️❌:                 && /*  do not compare single stereocenters in case of relative stereo */
    // INCHI✔️❌:                 !((nFlags1 & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) &&
    // INCHI✔️❌:                     1 == Stereo1->nNumberOfStereoCenters
    // INCHI✔️❌:                     )
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                 ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompareInchiStereo

    if let (Some(stereo1), Some(stereo2)) = (stereo1, stereo2) {
        let bond_count = stereo1
            .nNumberOfStereoBonds
            .min(stereo2.nNumberOfStereoBonds);
        if bond_count > 0 {
            let count =
                usize::try_from(bond_count).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let atom11 = heap
                .slice(stereo1.nBondAtom1.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let atom12 = heap
                .slice(stereo2.nBondAtom1.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let atom21 = heap
                .slice(stereo1.nBondAtom2.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let atom22 = heap
                .slice(stereo2.nBondAtom2.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let parity1 = heap
                .slice(stereo1.b_parity.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let parity2 = heap
                .slice(stereo2.b_parity.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            for index in 0..count {
                let ret = i32::from(atom12[index]) - i32::from(atom11[index]);
                if ret != 0 {
                    return Ok(ret);
                }
                let ret = i32::from(atom22[index]) - i32::from(atom21[index]);
                if ret != 0 {
                    return Ok(ret);
                }
                let ret = i32::from(parity2[index]) - i32::from(parity1[index]);
                if ret != 0 {
                    return Ok(ret);
                }
            }
        }
        let ret = stereo2.nNumberOfStereoBonds - stereo1.nNumberOfStereoBonds;
        if ret != 0 {
            return Ok(ret);
        }
        let center_count = stereo1
            .nNumberOfStereoCenters
            .min(stereo2.nNumberOfStereoCenters);
        if center_count > 0 {
            let count = usize::try_from(center_count)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let number1 = heap
                .slice(stereo1.nNumber.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let number2 = heap
                .slice(stereo2.nNumber.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let parity1 = heap
                .slice(stereo1.t_parity.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let parity2 = heap
                .slice(stereo2.t_parity.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            for index in 0..count {
                let ret = i32::from(number2[index]) - i32::from(number1[index]);
                if ret != 0 {
                    return Ok(ret);
                }
                let ret = i32::from(parity2[index]) - i32::from(parity1[index]);
                if ret != 0 {
                    return Ok(ret);
                }
            }
        }
        let ret = stereo2.nNumberOfStereoCenters - stereo1.nNumberOfStereoCenters;
        if ret != 0 {
            return Ok(ret);
        }
        if (flags1 | flags2) & INCHI_MODE::from(INCHI_FLAG_RAC_STEREO | INCHI_FLAG_REL_STEREO) == 0
        {
            let ret = i32::from(stereo2.nCompInv2Abs < 0) - i32::from(stereo1.nCompInv2Abs < 0);
            if ret != 0 {
                return Ok(ret);
            }
        }
    } else {
        if stereo2.is_some_and(|stereo| {
            stereo.nNumberOfStereoBonds > 0 || stereo.nNumberOfStereoCenters > 0
        }) {
            return Ok(1);
        }
        if stereo1.is_some_and(|stereo| {
            stereo.nNumberOfStereoBonds > 0 || stereo.nNumberOfStereoCenters > 0
        }) {
            return Ok(-1);
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CompINChI2(
    heap: &mut SourceHeap,
    p1: &INCHI_SORT,
    p2: &INCHI_SORT,
    b_taut: u32,
    compare_isotopic: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:1712 CompINChI2
    // INCHI✔️❌: int CompINChI2(const INCHI_SORT* p1,
    // INCHI✔️❌:     const INCHI_SORT* p2,
    // INCHI✔️❌:     int bTaut,
    // INCHI✔️❌:     int bCompareIsotopic)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret, num, i, num_H1, num_H2;
    // INCHI✔️❌:
    // INCHI✔️❌:     const INChI* i1 = NULL; /* tautomeric if exists, otherwise non-tautomeric */
    // INCHI✔️❌:     const INChI* i2 = NULL; /* tautomeric if exists, otherwise non-tautomeric */
    // INCHI✔️❌:
    // INCHI✔️❌:     int   n1;               /* TAUT_YES if tautomeric i1 exists, otherwise TAUT_NON */
    // INCHI✔️❌:     int   n2;               /* TAUT_YES if tautomeric i2 exists, otherwise TAUT_NON */
    // INCHI✔️❌:
    // INCHI✔️❌:     const INChI* i1n = NULL; /* non-tautomeric if both tautomeric AND non-tautomeric exist */
    // INCHI✔️❌:     const INChI* i2n = NULL; /* non-tautomeric if both tautomeric AND non-tautomeric exist */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*const INChI *i1t = NULL;*/ /* temp for i1 if both tautomeric AND non-tautomeric exist */
    // INCHI✔️❌:     /*const INChI *i2t = NULL;*/ /* temp for i2 if both tautomeric AND non-tautomeric exist */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* INChI_Stereo *Stereo1, *Stereo2; */
    // INCHI✔️❌:
    // INCHI✔️❌:     n1 = (p1->pINChI[TAUT_YES] && p1->pINChI[TAUT_YES]->nNumberOfAtoms) ? TAUT_YES : TAUT_NON;
    // INCHI✔️❌:     n2 = (p2->pINChI[TAUT_YES] && p2->pINChI[TAUT_YES]->nNumberOfAtoms) ? TAUT_YES : TAUT_NON;
    // INCHI✔️❌:
    // INCHI✔️❌:     i1 = p1->pINChI[n1];
    // INCHI✔️❌:     i1n = (n1 == TAUT_YES && p1->pINChI[TAUT_NON] &&
    // INCHI✔️❌:         p1->pINChI[TAUT_NON]->nNumberOfAtoms) ? p1->pINChI[TAUT_NON] : (const INChI*)NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     i2 = p2->pINChI[n2];
    // INCHI✔️❌:     i2n = (n2 == TAUT_YES && p2->pINChI[TAUT_NON] &&
    // INCHI✔️❌:         p2->pINChI[TAUT_NON]->nNumberOfAtoms) ? p2->pINChI[TAUT_NON] : (const INChI*)NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* non-deleted-non-empty < deleted < empty */
    // INCHI✔️❌:     if (i1 && !i2)
    // INCHI✔️❌:         return -1;   /* non-empty is the smallest (first) */
    // INCHI✔️❌:     if (!i1 && i2)
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     if (!i1 && !i2)
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     if (i1->bDeleted && !i2->bDeleted)
    // INCHI✔️❌:         return 1;    /* deleted is the largest (last) among non-empty */
    // INCHI✔️❌:     if (!i1->bDeleted && i2->bDeleted)
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:
    // INCHI✔️❌:     num_H1 = num_H2 = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* do not compare terminal H */
    // INCHI✔️❌:     if ((ret = CompareHillFormulasNoH(i1->szHillFormula, i2->szHillFormula, &num_H1, &num_H2))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret;  /* lexicographic order except the shorter one is greater (last): CH2O < CH2; C3XX < C2XX */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*********************************************************
    // INCHI✔️❌:             compare non-isotopic non-tautomeric part
    // INCHI✔️❌:      *********************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:      /* compare number of atoms (excluding terminal H) */
    // INCHI✔️❌:     if ((ret = i2->nNumberOfAtoms - i1->nNumberOfAtoms)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         return ret; /*  more atoms first */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  compare elements  (excluding terminal H) */
    // INCHI✔️❌:     num = i1->nNumberOfAtoms;
    // INCHI✔️❌:     for (i = 0; i < num; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* should always be equal if Hill formulas are same */
    // INCHI✔️❌:         if ((ret = (int)i2->nAtom[i] - (int)i1->nAtom[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             return ret; /* greater periodic number first */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /**********************************************************
    // INCHI✔️❌:         compare connection tables
    // INCHI✔️❌:     ***********************************************************/
    // INCHI✔️❌:     if ((ret = i2->lenConnTable - i1->lenConnTable)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         return ret; /* longer connection table first */
    // INCHI✔️❌:     num = i2->lenConnTable;
    // INCHI✔️❌:     for (i = 0; i < num; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((ret = (int)i2->nConnTable[i] - (int)i1->nConnTable[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             return ret; /* greater connection table first */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*********************************************************
    // INCHI✔️❌:       compare compare total number of H (inverse: H3 < H2 )
    // INCHI✔️❌:     **********************************************************/
    // INCHI✔️❌:     if ((ret = num_H2 - num_H1)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*********************************************************
    // INCHI✔️❌:       compare non-tautomeric num_H: N < NH3 < NH2 < NH
    // INCHI✔️❌:     **********************************************************/
    // INCHI✔️❌:     num = i1->nNumberOfAtoms;
    // INCHI✔️❌:     for (i = 0; i < num; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (i2->nNum_H[i] != i1->nNum_H[i])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return !i2->nNum_H[i] ? 1 :  /* no H first */
    // INCHI✔️❌:                 !i1->nNum_H[i] ? -1 :
    // INCHI✔️❌:                 (int)i2->nNum_H[i] - (int)i1->nNum_H[i];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*********************************************************
    // INCHI✔️❌:          compare non-isotopic tautomeric part
    // INCHI✔️❌:      *********************************************************/
    // INCHI✔️❌:     if ((ret = CompareTautNonIsoPartOfINChI(i1, i2))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( ret = i2->lenTautomer - i1->lenTautomer )
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     num = inchi_min( i2->lenTautomer, i1->lenTautomer );
    // INCHI✔️❌:     for ( i = 0; i < num; i ++ ) {
    // INCHI✔️❌:         if ( ret = (int)i2->nTautomer[i] - (int)i1->nTautomer[i] )
    // INCHI✔️❌:            return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*********************************************************
    // INCHI✔️❌:      *                                                       *
    // INCHI✔️❌:      *  at this point both components are either tautomeric  *
    // INCHI✔️❌:      *  or non-tautomeric                                    *
    // INCHI✔️❌:      *                                                       *
    // INCHI✔️❌:      *********************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:      /*********************************************************
    // INCHI✔️❌:         non-tautomeric "fixed H" specific
    // INCHI✔️❌:       *********************************************************/
    // INCHI✔️❌:     if (TAUT_NON == bTaut && ((i1n && i1n->nNum_H_fixed) || (i2n && i2n->nNum_H_fixed))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* first, compare non-tautomeric chem. formulas -- they may be different */
    // INCHI✔️❌:         const char* f1 = (i1n /*&& i1n->nNum_H_fixed*/) ? i1n->szHillFormula : i1->szHillFormula;
    // INCHI✔️❌:         const char* f2 = (i2n /*&& i2n->nNum_H_fixed*/) ? i2n->szHillFormula : i2->szHillFormula;
    // INCHI✔️❌:         if (f1 && f2 && (ret = CompareHillFormulas(f1, f2)))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* secondly, compare fixed-H distribution */
    // INCHI✔️❌:         if (i1n && i1n->nNum_H_fixed && i2n && i2n->nNum_H_fixed)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num = inchi_min(i1n->nNumberOfAtoms, i2n->nNumberOfAtoms);
    // INCHI✔️❌:             for (i = 0; i < num; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i2n->nNum_H_fixed[i] != i1n->nNum_H_fixed[i])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return !i2n->nNum_H_fixed[i] ? 1 : /* no fixed H first */
    // INCHI✔️❌:                         !i1n->nNum_H_fixed[i] ? -1 :
    // INCHI✔️❌:                         (int)i2n->nNum_H_fixed[i] - (int)i1n->nNum_H_fixed[i];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)i2n->nNumberOfAtoms - (int)i1n->nNumberOfAtoms)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret; /* should not happen <BRKPT> */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i1n && i1n->nNum_H_fixed)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 num = i1n->nNumberOfAtoms;
    // INCHI✔️❌:                 for (i = 0; i < num; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* added 2004-05-04 */
    // INCHI✔️❌:                     if (i1n->nNum_H_fixed[i])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return -1; /* i1n->nNum_H_fixed[i] > 0? -1:1;*/
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* p1 is tautomeric, p2 is not tautomeric; this must have been detected earlier */
    // INCHI✔️❌:                 /*return -1;*/ /* has fixed H first *//* <BRKPT> */ /* removed 2004-05-04 */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 num = i2n->nNumberOfAtoms;
    // INCHI✔️❌:                 for (i = 0; i < num; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* added 2004-05-04 */
    // INCHI✔️❌:                     if (i2n->nNum_H_fixed[i])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return 1; /* i2n->nNum_H_fixed[i] > 0? 1:-1;*/
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* p2 is tautomeric, p1 is not tautomeric; this must have been detected earlier */
    // INCHI✔️❌:                 /*return 1; */ /* has fixed H first *//* <BRKPT> */ /* removed 2004-05-04 */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*************************************************************************
    // INCHI✔️❌:               if requested non-tautomeric comparison then
    // INCHI✔️❌:               prepare to compare non-taut non-isotopic stereo, etc.
    // INCHI✔️❌:      *************************************************************************/
    // INCHI✔️❌:     if (TAUT_NON == bTaut)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (i1n)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*i1t = i1;*/
    // INCHI✔️❌:             i1 = i1n;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i2n)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*i2t = i2;*/
    // INCHI✔️❌:             i2 = i2n;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*********************************************************
    // INCHI✔️❌:         compare non-isotopic stereo
    // INCHI✔️❌:      *********************************************************/
    // INCHI✔️❌:     ret = CompareInchiStereo(i1->Stereo, i1->nFlags, i2->Stereo, i2->nFlags);
    // INCHI✔️❌:     if (ret)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*******************************************************
    // INCHI✔️❌:         do not switch back to tautomeric i1, i2
    // INCHI✔️❌:      *******************************************************/
    // INCHI✔️❌:      /* -- how to switch back --
    // INCHI✔️❌:      if ( i1t ) {
    // INCHI✔️❌:          i1  = i1t;
    // INCHI✔️❌:          i1t = NULL;
    // INCHI✔️❌:      }
    // INCHI✔️❌:      if ( i2t ) {
    // INCHI✔️❌:          i2  = i2t;
    // INCHI✔️❌:          i2t = NULL;
    // INCHI✔️❌:      }
    // INCHI✔️❌:      */
    // INCHI✔️❌:
    // INCHI✔️❌:      /******************************************************
    // INCHI✔️❌:           compare isotopic non-tautomeric part
    // INCHI✔️❌:       ******************************************************/
    // INCHI✔️❌:     if (bCompareIsotopic)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((ret = i2->nNumberOfIsotopicAtoms - i1->nNumberOfIsotopicAtoms)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         num = i1->nNumberOfIsotopicAtoms;
    // INCHI✔️❌:         /*  compare isotopic atoms */
    // INCHI✔️❌:         for (i = 0; i < num; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nAtomNumber - (int)i1->IsotopicAtom[i].nAtomNumber)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nIsoDifference - (int)i1->IsotopicAtom[i].nIsoDifference)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* compare isotopic H */
    // INCHI✔️❌:         /* if tautomeric comparison mode then here are compared only non-tautomeric H */
    // INCHI✔️❌:         for (i = 0; i < num; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nNum_T - (int)i1->IsotopicAtom[i].nNum_T)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nNum_D - (int)i1->IsotopicAtom[i].nNum_D)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nNum_H - (int)i1->IsotopicAtom[i].nNum_H)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /*****************************************************
    // INCHI✔️❌:              compare isotopic tautomeric part
    // INCHI✔️❌:          *****************************************************/
    // INCHI✔️❌:         if ((ret = i2->nNumberOfIsotopicTGroups - i1->nNumberOfIsotopicTGroups)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         num = i1->nNumberOfIsotopicTGroups;
    // INCHI✔️❌:         for (i = 0; i < num; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicTGroup[i].nTGroupNumber - (int)i1->IsotopicTGroup[i].nTGroupNumber)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicTGroup[i].nNum_T - (int)i1->IsotopicTGroup[i].nNum_T)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicTGroup[i].nNum_D - (int)i1->IsotopicTGroup[i].nNum_D)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicTGroup[i].nNum_H - (int)i1->IsotopicTGroup[i].nNum_H)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /****************************************************
    // INCHI✔️❌:             compare isotopic stereo
    // INCHI✔️❌:          ****************************************************/
    // INCHI✔️❌:         ret = CompareInchiStereo(i1->StereoIsotopic, i1->nFlags, i2->StereoIsotopic, i2->nFlags);
    // INCHI✔️❌:         if (ret)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /**********************************************************
    // INCHI✔️❌:         compare charges: non-charged first, then in order of
    // INCHI✔️❌:         ascending charges (negative first)
    // INCHI✔️❌:     ***********************************************************/
    // INCHI✔️❌:     if (i2->nTotalCharge && i1->nTotalCharge)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  both are charged; smaller charges first */
    // INCHI✔️❌:         ret = (int)i1->nTotalCharge - (int)i2->nTotalCharge;
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if ((ret = (i1->nTotalCharge ? 1 : 0) - (i2->nTotalCharge ? 1 : 0))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  only one is charged; uncharged first */
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* stable sort */
    // INCHI✔️❌:     /*ret = p1->ord_number - p2->ord_number;*/
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompINChI2

    let taut_yes = TAUT_YES as usize;
    let taut_non = TAUT_NON as usize;
    let p1_taut = p1.pINChI[taut_yes];
    let p2_taut = p2.pINChI[taut_yes];
    let n1 = if !p1_taut.is_null() && heap.slice(p1_taut.as_const())?[0].nNumberOfAtoms != 0 {
        taut_yes
    } else {
        taut_non
    };
    let n2 = if !p2_taut.is_null() && heap.slice(p2_taut.as_const())?[0].nNumberOfAtoms != 0 {
        taut_yes
    } else {
        taut_non
    };
    let i1_pointer = p1.pINChI[n1];
    let i2_pointer = p2.pINChI[n2];
    let mut i1 = if i1_pointer.is_null() {
        None
    } else {
        Some(heap.slice(i1_pointer.as_const())?[0].clone())
    };
    let mut i2 = if i2_pointer.is_null() {
        None
    } else {
        Some(heap.slice(i2_pointer.as_const())?[0].clone())
    };
    let p1_non = p1.pINChI[taut_non];
    let p2_non = p2.pINChI[taut_non];
    let i1n = if n1 == taut_yes
        && !p1_non.is_null()
        && heap.slice(p1_non.as_const())?[0].nNumberOfAtoms != 0
    {
        Some(heap.slice(p1_non.as_const())?[0].clone())
    } else {
        None
    };
    let i2n = if n2 == taut_yes
        && !p2_non.is_null()
        && heap.slice(p2_non.as_const())?[0].nNumberOfAtoms != 0
    {
        Some(heap.slice(p2_non.as_const())?[0].clone())
    } else {
        None
    };

    match (&i1, &i2) {
        (Some(_), None) => return Ok(-1),
        (None, Some(_)) => return Ok(1),
        (None, None) => return Ok(0),
        (Some(first), Some(second)) if first.bDeleted != 0 && second.bDeleted == 0 => return Ok(1),
        (Some(first), Some(second)) if first.bDeleted == 0 && second.bDeleted != 0 => {
            return Ok(-1);
        }
        _ => {}
    }
    let first = i1.as_ref().ok_or(SourceHeapError::NullPointer)?;
    let second = i2.as_ref().ok_or(SourceHeapError::NullPointer)?;
    let mut num_h1 = 0;
    let mut num_h2 = 0;
    let ret = CompareHillFormulasNoH(
        heap,
        first.szHillFormula.as_const(),
        second.szHillFormula.as_const(),
        &mut num_h1,
        &mut num_h2,
    )?;
    if ret != 0 {
        return Ok(ret);
    }
    let ret = second.nNumberOfAtoms.wrapping_sub(first.nNumberOfAtoms);
    if ret != 0 {
        return Ok(ret);
    }
    let atom_count = usize::try_from(first.nNumberOfAtoms.max(0))
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if atom_count > 0 {
        let atoms1 = heap
            .slice(first.nAtom.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let atoms2 = heap
            .slice(second.nAtom.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for index in 0..atom_count {
            let ret = i32::from(atoms2[index]) - i32::from(atoms1[index]);
            if ret != 0 {
                return Ok(ret);
            }
        }
    }
    let ret = second.lenConnTable.wrapping_sub(first.lenConnTable);
    if ret != 0 {
        return Ok(ret);
    }
    let connection_count = usize::try_from(second.lenConnTable.max(0))
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if connection_count > 0 {
        let connections1 = heap
            .slice(first.nConnTable.as_const())?
            .get(..connection_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let connections2 = heap
            .slice(second.nConnTable.as_const())?
            .get(..connection_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for index in 0..connection_count {
            let ret = i32::from(connections2[index]) - i32::from(connections1[index]);
            if ret != 0 {
                return Ok(ret);
            }
        }
    }
    let ret = num_h2.wrapping_sub(num_h1);
    if ret != 0 {
        return Ok(ret);
    }
    if atom_count > 0 {
        let hydrogens1 = heap
            .slice(first.nNum_H.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let hydrogens2 = heap
            .slice(second.nNum_H.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for index in 0..atom_count {
            if hydrogens2[index] != hydrogens1[index] {
                return Ok(if hydrogens2[index] == 0 {
                    1
                } else if hydrogens1[index] == 0 {
                    -1
                } else {
                    i32::from(hydrogens2[index]) - i32::from(hydrogens1[index])
                });
            }
        }
    }
    let ret = CompareTautNonIsoPartOfINChI(heap, first, second)?;
    if ret != 0 {
        return Ok(ret);
    }

    if b_taut == TAUT_NON
        && (i1n
            .as_ref()
            .is_some_and(|value| !value.nNum_H_fixed.is_null())
            || i2n
                .as_ref()
                .is_some_and(|value| !value.nNum_H_fixed.is_null()))
    {
        let formula1 = i1n.as_ref().unwrap_or(first).szHillFormula;
        let formula2 = i2n.as_ref().unwrap_or(second).szHillFormula;
        if !formula1.is_null() && !formula2.is_null() {
            let ret = CompareHillFormulas(heap, formula1.as_const(), formula2.as_const())?;
            if ret != 0 {
                return Ok(ret);
            }
        }
        let fixed1 = i1n.as_ref().filter(|value| !value.nNum_H_fixed.is_null());
        let fixed2 = i2n.as_ref().filter(|value| !value.nNum_H_fixed.is_null());
        match (fixed1, fixed2) {
            (Some(left), Some(right)) => {
                let count = usize::try_from(left.nNumberOfAtoms.min(right.nNumberOfAtoms).max(0))
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let values1 = heap
                    .slice(left.nNum_H_fixed.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let values2 = heap
                    .slice(right.nNum_H_fixed.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                for index in 0..count {
                    if values2[index] != values1[index] {
                        return Ok(if values2[index] == 0 {
                            1
                        } else if values1[index] == 0 {
                            -1
                        } else {
                            i32::from(values2[index]) - i32::from(values1[index])
                        });
                    }
                }
                let ret = right.nNumberOfAtoms.wrapping_sub(left.nNumberOfAtoms);
                if ret != 0 {
                    return Ok(ret);
                }
            }
            (Some(left), None) => {
                let count = usize::try_from(left.nNumberOfAtoms.max(0))
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                if heap
                    .slice(left.nNum_H_fixed.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .iter()
                    .any(|value| *value != 0)
                {
                    return Ok(-1);
                }
            }
            (None, Some(right)) => {
                let count = usize::try_from(right.nNumberOfAtoms.max(0))
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                if heap
                    .slice(right.nNum_H_fixed.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .iter()
                    .any(|value| *value != 0)
                {
                    return Ok(1);
                }
            }
            (None, None) => {}
        }
    }

    if b_taut == TAUT_NON {
        if let Some(fixed) = i1n {
            i1 = Some(fixed);
        }
        if let Some(fixed) = i2n {
            i2 = Some(fixed);
        }
    }
    let first = i1.as_ref().ok_or(SourceHeapError::NullPointer)?;
    let second = i2.as_ref().ok_or(SourceHeapError::NullPointer)?;
    let stereo1 = if first.Stereo.is_null() {
        None
    } else {
        Some(heap.slice(first.Stereo.as_const())?[0].clone())
    };
    let stereo2 = if second.Stereo.is_null() {
        None
    } else {
        Some(heap.slice(second.Stereo.as_const())?[0].clone())
    };
    let ret = CompareInchiStereo(
        heap,
        stereo1.as_ref(),
        first.nFlags,
        stereo2.as_ref(),
        second.nFlags,
    )?;
    if ret != 0 {
        return Ok(ret);
    }

    if compare_isotopic != 0 {
        let ret = second
            .nNumberOfIsotopicAtoms
            .wrapping_sub(first.nNumberOfIsotopicAtoms);
        if ret != 0 {
            return Ok(ret);
        }
        let isotope_count = usize::try_from(first.nNumberOfIsotopicAtoms.max(0))
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if isotope_count > 0 {
            let isotopes1 = heap
                .slice(first.IsotopicAtom.as_const())?
                .get(..isotope_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let isotopes2 = heap
                .slice(second.IsotopicAtom.as_const())?
                .get(..isotope_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            for index in 0..isotope_count {
                for ret in [
                    i32::from(isotopes2[index].nAtomNumber)
                        - i32::from(isotopes1[index].nAtomNumber),
                    i32::from(isotopes2[index].nIsoDifference)
                        - i32::from(isotopes1[index].nIsoDifference),
                ] {
                    if ret != 0 {
                        return Ok(ret);
                    }
                }
            }
            for index in 0..isotope_count {
                for ret in [
                    i32::from(isotopes2[index].nNum_T) - i32::from(isotopes1[index].nNum_T),
                    i32::from(isotopes2[index].nNum_D) - i32::from(isotopes1[index].nNum_D),
                    i32::from(isotopes2[index].nNum_H) - i32::from(isotopes1[index].nNum_H),
                ] {
                    if ret != 0 {
                        return Ok(ret);
                    }
                }
            }
        }
        let ret = second
            .nNumberOfIsotopicTGroups
            .wrapping_sub(first.nNumberOfIsotopicTGroups);
        if ret != 0 {
            return Ok(ret);
        }
        let group_count = usize::try_from(first.nNumberOfIsotopicTGroups.max(0))
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if group_count > 0 {
            let groups1 = heap
                .slice(first.IsotopicTGroup.as_const())?
                .get(..group_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let groups2 = heap
                .slice(second.IsotopicTGroup.as_const())?
                .get(..group_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            for index in 0..group_count {
                for ret in [
                    i32::from(groups2[index].nTGroupNumber)
                        - i32::from(groups1[index].nTGroupNumber),
                    i32::from(groups2[index].nNum_T) - i32::from(groups1[index].nNum_T),
                    i32::from(groups2[index].nNum_D) - i32::from(groups1[index].nNum_D),
                    i32::from(groups2[index].nNum_H) - i32::from(groups1[index].nNum_H),
                ] {
                    if ret != 0 {
                        return Ok(ret);
                    }
                }
            }
        }
        let stereo1 = if first.StereoIsotopic.is_null() {
            None
        } else {
            Some(heap.slice(first.StereoIsotopic.as_const())?[0].clone())
        };
        let stereo2 = if second.StereoIsotopic.is_null() {
            None
        } else {
            Some(heap.slice(second.StereoIsotopic.as_const())?[0].clone())
        };
        let ret = CompareInchiStereo(
            heap,
            stereo1.as_ref(),
            first.nFlags,
            stereo2.as_ref(),
            second.nFlags,
        )?;
        if ret != 0 {
            return Ok(ret);
        }
    }

    if second.nTotalCharge != 0 && first.nTotalCharge != 0 {
        return Ok(first.nTotalCharge.wrapping_sub(second.nTotalCharge));
    }
    Ok(i32::from(first.nTotalCharge != 0) - i32::from(second.nTotalCharge != 0))
}

#[allow(non_snake_case)]
pub(crate) fn CompINChINonTaut2(
    heap: &mut SourceHeap,
    p1: &INCHI_SORT,
    p2: &INCHI_SORT,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2042 CompINChINonTaut2
    // INCHI✔️❌: int CompINChINonTaut2(const void* p1, const void* p2)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret;
    // INCHI✔️❌:     ret = CompINChI2((const INCHI_SORT*)p1, (const INCHI_SORT*)p2, TAUT_NON, 1);
    // INCHI✔️❌: #if ( CANON_FIXH_TRANS == 1 )
    // INCHI✔️❌:     if (!ret)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* to obtain canonical transposition 2004-05-10 */
    // INCHI✔️❌:         ret = CompINChI2((const INCHI_SORT*)p1, (const INCHI_SORT*)p2, TAUT_YES, 1);
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     if (!ret)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* stable sort */
    // INCHI✔️❌:         ret = ((const INCHI_SORT*)p1)->ord_number - ((const INCHI_SORT*)p2)->ord_number;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompINChINonTaut2

    let mut ret = CompINChI2(heap, p1, p2, TAUT_NON, 1)?;
    if ret == 0 {
        ret = CompINChI2(heap, p1, p2, TAUT_YES, 1)?;
    }
    if ret == 0 {
        ret = i32::from(p1.ord_number) - i32::from(p2.ord_number);
    }
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn CompINChITaut2(
    heap: &mut SourceHeap,
    p1: &INCHI_SORT,
    p2: &INCHI_SORT,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2064 CompINChITaut2
    // INCHI✔️❌: int CompINChITaut2(const void* p1, const void* p2)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret;
    // INCHI✔️❌:     ret = CompINChI2((const INCHI_SORT*)p1, (const INCHI_SORT*)p2, TAUT_YES, 1);
    // INCHI✔️❌: #if ( CANON_FIXH_TRANS == 1 )
    // INCHI✔️❌:     if (!ret)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* to obtain canonical transposition 2004-05-10 */
    // INCHI✔️❌:         ret = CompINChI2((const INCHI_SORT*)p1, (const INCHI_SORT*)p2, TAUT_NON, 1);
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     if (!ret)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* stable sort */
    // INCHI✔️❌:         ret = ((const INCHI_SORT*)p1)->ord_number - ((const INCHI_SORT*)p2)->ord_number;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompINChITaut2

    let mut ret = CompINChI2(heap, p1, p2, TAUT_YES, 1)?;
    if ret == 0 {
        ret = CompINChI2(heap, p1, p2, TAUT_NON, 1)?;
    }
    if ret == 0 {
        ret = i32::from(p1.ord_number) - i32::from(p2.ord_number);
    }
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn INChI_SegmentAction(c_dif_segs: i8) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:1486 INChI_SegmentAction
    // INCHI✔️✔️: int INChI_SegmentAction(char cDifSegs)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     if (!(cDifSegs & DIFV_OUTPUT_OMIT_F))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return INCHI_SEGM_OMIT;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if ((cDifSegs & DIFV_OUTPUT_EMPTY_T) && !(cDifSegs & DIFV_OUTPUT_EMPTY_F))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return INCHI_SEGM_EMPTY;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if ((cDifSegs & DIFV_OUTPUT_FILL_T))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return INCHI_SEGM_FILL;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return INCHI_SEGM_OMIT; /* the control flow shoul never reach this point */
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: INChI_SegmentAction

    let value = i32::from(c_dif_segs);
    if value & crate::source_types::tagMarkDiff_DIFV_OUTPUT_OMIT_F as i32 == 0 {
        return crate::source_types::tagINChISegmAction_INCHI_SEGM_OMIT as i32;
    }
    if value & crate::source_types::tagMarkDiff_DIFV_OUTPUT_EMPTY_T as i32 != 0
        && value & crate::source_types::tagMarkDiff_DIFV_OUTPUT_EMPTY_F as i32 == 0
    {
        return crate::source_types::tagINChISegmAction_INCHI_SEGM_EMPTY as i32;
    }
    if value & crate::source_types::tagMarkDiff_DIFV_OUTPUT_FILL_T as i32 != 0 {
        return crate::source_types::tagINChISegmAction_INCHI_SEGM_FILL as i32;
    }
    crate::source_types::tagINChISegmAction_INCHI_SEGM_OMIT as i32
}

#[allow(non_snake_case)]
pub(crate) fn MarkUnusedAndEmptyLayers(s_dif_segs: &mut [[i8; 11]; 4]) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:1525 MarkUnusedAndEmptyLayers
    // INCHI✔️❌: int MarkUnusedAndEmptyLayers(char sDifSegs[][DIFS_LENGTH])
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, nLayer, sBits, nFirstSegm;
    // INCHI✔️❌: #define nFirstFmlSegm   DIFS_f_FORMULA
    // INCHI✔️❌: #define nFirstIsoSegm   DIFS_i_IATOMS
    // INCHI✔️❌:     /* FI */
    // INCHI✔️❌:     nLayer = DIFL_FI;
    // INCHI✔️❌:     nFirstSegm = nFirstIsoSegm;
    // INCHI✔️❌:     sBits = 0;
    // INCHI✔️❌:     for (i = 0; i < DIFS_idf_LENGTH; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sBits |= sDifSegs[nLayer][i];
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (!(sBits & DIFV_OUTPUT_OMIT_F))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Omit the FI layer */
    // INCHI✔️❌:         memset(sDifSegs[nLayer], DIFV_BOTH_EMPTY, DIFS_idf_LENGTH); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (sDifSegs[nLayer][nFirstSegm] == DIFV_BOTH_EMPTY ||
    // INCHI✔️❌:             !(sDifSegs[nLayer][nFirstSegm] & DIFV_OUTPUT_OMIT_F))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sDifSegs[nLayer][nFirstSegm] = DIFV_IS_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* MI */
    // INCHI✔️❌:     nLayer = DIFL_MI;
    // INCHI✔️❌:     nFirstSegm = nFirstIsoSegm;
    // INCHI✔️❌:     sBits = 0;
    // INCHI✔️❌:     for (i = 0; i < DIFS_idf_LENGTH; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sBits |= sDifSegs[nLayer][i];
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (!(sBits & DIFV_OUTPUT_OMIT_F))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Omit the MI layer */
    // INCHI✔️❌:         memset(sDifSegs[nLayer], DIFV_BOTH_EMPTY, DIFS_idf_LENGTH); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (sDifSegs[nLayer][nFirstSegm] == DIFV_BOTH_EMPTY ||
    // INCHI✔️❌:             !(sDifSegs[nLayer][nFirstSegm] & DIFV_OUTPUT_OMIT_F))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sDifSegs[nLayer][nFirstSegm] = DIFV_IS_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* F */
    // INCHI✔️❌:     nLayer = DIFL_F;
    // INCHI✔️❌:     nFirstSegm = nFirstFmlSegm;
    // INCHI✔️❌:     sBits = 0;
    // INCHI✔️❌:     for (i = 0; i < DIFS_idf_LENGTH; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sBits |= sDifSegs[nLayer][i];
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (!(sBits & DIFV_OUTPUT_OMIT_F) &&
    // INCHI✔️❌:         sDifSegs[DIFL_FI][nFirstIsoSegm] == DIFV_BOTH_EMPTY)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Omit the F layer: no non-iotopic and no isotopic segments */
    // INCHI✔️❌:         memset(sDifSegs[nLayer], DIFV_BOTH_EMPTY, DIFS_idf_LENGTH); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {    /* do not omit fixed-H layer */
    // INCHI✔️❌:         if (sDifSegs[nLayer][nFirstSegm] == DIFV_BOTH_EMPTY ||
    // INCHI✔️❌:             !(sDifSegs[nLayer][nFirstSegm] & DIFV_OUTPUT_OMIT_F))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sDifSegs[nLayer][nFirstSegm] = DIFV_IS_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* M -- leave as it is */
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌:
    // INCHI✔️❌: #undef nFirstFmlSegm
    // INCHI✔️❌: #undef nFirstIsoSegm
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MarkUnusedAndEmptyLayers

    let first_iso = tagDiffINChISegments_DIFS_i_IATOMS as usize;
    let output_omit = tagMarkDiff_DIFV_OUTPUT_OMIT_F as i8;
    for layer in [
        tagDiffINChILayers_DIFL_FI as usize,
        tagDiffINChILayers_DIFL_MI as usize,
    ] {
        let bits = s_dif_segs[layer]
            .iter()
            .fold(0_i8, |bits, value| bits | value);
        if bits & output_omit == 0 {
            s_dif_segs[layer].fill(tagMarkDiff_DIFV_BOTH_EMPTY as i8);
        } else if s_dif_segs[layer][first_iso] == tagMarkDiff_DIFV_BOTH_EMPTY as i8
            || s_dif_segs[layer][first_iso] & output_omit == 0
        {
            s_dif_segs[layer][first_iso] = tagMarkDiff_DIFV_IS_EMPTY as i8;
        }
    }
    let fixed = tagDiffINChILayers_DIFL_F as usize;
    let formula = tagDiffINChISegments_DIFS_f_FORMULA as usize;
    let bits = s_dif_segs[fixed]
        .iter()
        .fold(0_i8, |bits, value| bits | value);
    if bits & output_omit == 0
        && s_dif_segs[tagDiffINChILayers_DIFL_FI as usize][first_iso]
            == tagMarkDiff_DIFV_BOTH_EMPTY as i8
    {
        s_dif_segs[fixed].fill(tagMarkDiff_DIFV_BOTH_EMPTY as i8);
    } else if s_dif_segs[fixed][formula] == tagMarkDiff_DIFV_BOTH_EMPTY as i8
        || s_dif_segs[fixed][formula] & output_omit == 0
    {
        s_dif_segs[fixed][formula] = tagMarkDiff_DIFV_IS_EMPTY as i8;
    }
    0
}

#[allow(non_snake_case)]
pub(crate) fn CompareReversedStereoINChI(
    heap: &SourceHeap,
    s1: Option<&INChI_Stereo>,
    s2: Option<&INChI_Stereo>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2555 CompareReversedStereoINChI
    // INCHI✔️❌: int CompareReversedStereoINChI(INChI_Stereo* s1/* InChI from reversed struct */,
    // INCHI✔️❌:     INChI_Stereo* s2 /* input InChI */
    // INCHI✔️❌: )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (s1 == NULL && s2 == NULL)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if ((s1 == NULL) ^ (s2 == NULL))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         INChI_Stereo* s = s1 ? s1 : s2;
    // INCHI✔️❌:         if (s->nNumberOfStereoCenters || s->nNumberOfStereoBonds)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 20; /* Diff: Missing Stereo */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (s1->nNumberOfStereoCenters != s2->nNumberOfStereoCenters)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 21;      /* Diff: Number of sp3 stereocenters */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (s1->nNumberOfStereoCenters > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (memcmp(s1->nNumber, s2->nNumber, s1->nNumberOfStereoCenters * sizeof(s1->nNumber[0])))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 22;  /* Diff: sp3 stereocenter locations */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (memcmp(s1->t_parity, s2->t_parity, s1->nNumberOfStereoCenters * sizeof(s1->t_parity[0])))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 23;  /* Diff: sp3 stereocenter parities */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (s1->nCompInv2Abs != s2->nCompInv2Abs && s1->nCompInv2Abs && s2->nCompInv2Abs)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 24;  /* Diff: sp3 inversion */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         if ( s1->nNumberInv && s2->nNumberInv ) {
    // INCHI✔️❌:             if ( memcmp( s1->nNumberInv, s2->nNumberInv, s1->nNumberOfStereoCenters*sizeof(s1->nNumber[0]) ) )
    // INCHI✔️❌:                 return 25;
    // INCHI✔️❌:             if ( memcmp( s1->t_parityInv, s2->t_parityInv, s1->nNumberOfStereoCenters*sizeof(s1->t_parity[0]) ) )
    // INCHI✔️❌:                 return 26;
    // INCHI✔️❌:             if ( s1->nCompInv2Abs != s2->nCompInv2Abs ||
    // INCHI✔️❌:                  s1->bTrivialInv  != s2->bTrivialInv ) {
    // INCHI✔️❌:                 return 27;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         } else
    // INCHI✔️❌:         if ( s1->nNumberInv || s2->nNumberInv ) {
    // INCHI✔️❌:             return 28;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (s1->nNumberOfStereoBonds != s2->nNumberOfStereoBonds)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 25;      /* Diff: Number of stereobonds */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (s1->nNumberOfStereoBonds > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (memcmp(s1->nBondAtom1, s2->nBondAtom1, s1->nNumberOfStereoBonds * sizeof(s1->nBondAtom1[0])))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 26; /* Diff: Stereobond 1st atom locations */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (memcmp(s1->nBondAtom2, s2->nBondAtom2, s1->nNumberOfStereoBonds * sizeof(s1->nBondAtom2[0])))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 27; /* Diff: Stereobond 2nd atom locations */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (memcmp(s1->b_parity, s2->b_parity, s1->nNumberOfStereoBonds * sizeof(s1->b_parity[0])))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 28; /* Diff: Stereobond parities */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompareReversedStereoINChI

    let (Some(s1), Some(s2)) = (s1, s2) else {
        let Some(stereo) = s1.or(s2) else {
            return Ok(0);
        };
        return Ok(
            if stereo.nNumberOfStereoCenters != 0 || stereo.nNumberOfStereoBonds != 0 {
                20
            } else {
                0
            },
        );
    };
    if s1.nNumberOfStereoCenters != s2.nNumberOfStereoCenters {
        return Ok(21);
    }
    if s1.nNumberOfStereoCenters > 0 {
        let count = usize::try_from(s1.nNumberOfStereoCenters)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if heap
            .slice(s1.nNumber.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            != heap
                .slice(s2.nNumber.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        {
            return Ok(22);
        }
        if heap
            .slice(s1.t_parity.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            != heap
                .slice(s2.t_parity.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        {
            return Ok(23);
        }
        if s1.nCompInv2Abs != s2.nCompInv2Abs && s1.nCompInv2Abs != 0 && s2.nCompInv2Abs != 0 {
            return Ok(24);
        }
    }
    if s1.nNumberOfStereoBonds != s2.nNumberOfStereoBonds {
        return Ok(25);
    }
    if s1.nNumberOfStereoBonds > 0 {
        let count = usize::try_from(s1.nNumberOfStereoBonds)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if heap
            .slice(s1.nBondAtom1.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            != heap
                .slice(s2.nBondAtom1.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        {
            return Ok(26);
        }
        if heap
            .slice(s1.nBondAtom2.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            != heap
                .slice(s2.nBondAtom2.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        {
            return Ok(27);
        }
        if heap
            .slice(s1.b_parity.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            != heap
                .slice(s2.b_parity.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        {
            return Ok(28);
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CompareReversedINChI(
    heap: &SourceHeap,
    i1: Option<&INChI>,
    i2: Option<&INChI>,
    a1: Option<&INChI_Aux>,
    a2: Option<&INChI_Aux>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2936 CompareReversedINChI
    // INCHI✔️❌: int CompareReversedINChI(INChI* i1 /* InChI from reversed struct */,
    // INCHI✔️❌:     INChI* i2 /* input InChI */,
    // INCHI✔️❌:     INChI_Aux* a1,
    // INCHI✔️❌:     INChI_Aux* a2)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i1 == NULL && i2 == NULL)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if ((i1 == NULL) ^ (i2 == NULL))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 1; /* Diff: Missing InChI */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i1->nErrorCode == i2->nErrorCode)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (i1->nErrorCode)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 2; /* Diff: Error codes */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (i1->bDeleted != i2->bDeleted)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 1; /* Diff: Missing InChI */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (i1->nNumberOfAtoms != i2->nNumberOfAtoms)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 3;  /* Diff: Num. atoms */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i1->nNumberOfAtoms > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (memcmp(i1->nAtom, i2->nAtom, i1->nNumberOfAtoms * sizeof(i1->nAtom[0])))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 4; /* Diff: Elements */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (strcmp(i1->szHillFormula, i2->szHillFormula))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 7; /* Diff: Hill Formulas */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (memcmp(i1->nNum_H, i2->nNum_H, i1->nNumberOfAtoms * sizeof(i1->nNum_H[0])))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i1->lenConnTable > 1 || i2->lenConnTable > 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 5; /* Diff: H Locations (mobile H present) */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 6; /* Diff: H Locations (no mobile H) */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* fixed H */
    // INCHI✔️❌:         if (i1->nNum_H_fixed || i2->nNum_H_fixed)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int bHasFixedH1 = 0, bHasFixedH2 = 0, i, j1, j2;
    // INCHI✔️❌:             if (i1->nNum_H_fixed)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (i = 0; i < i1->nNumberOfAtoms; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (i1->nNum_H_fixed[i])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         bHasFixedH1++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (i2->nNum_H_fixed)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (i = 0; i < i2->nNumberOfAtoms; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (i2->nNum_H_fixed[i])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         bHasFixedH2++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* count the differences */
    // INCHI✔️❌:             j1 = j2 = 0;
    // INCHI✔️❌:             if (bHasFixedH1 && !bHasFixedH2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (i = 0; i < i1->nNumberOfAtoms; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (i1->nNum_H_fixed[i] > 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         j1++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                         if (i1->nNum_H_fixed[i] < 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             j2++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 return 18; /* Diff: Extra Fixed-H */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (!bHasFixedH1 && bHasFixedH2)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     for (i = j1 = j2 = 0; i < i1->nNumberOfAtoms; i++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (0 > i2->nNum_H_fixed[i])
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             j1++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                             if (0 < i2->nNum_H_fixed[i])
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 j2++;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     return 19; /* Diff: Missed Fixed-H */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (bHasFixedH1 && bHasFixedH2 &&
    // INCHI✔️❌:                         memcmp(i1->nNum_H_fixed, i2->nNum_H_fixed, i1->nNumberOfAtoms * sizeof(i1->nNum_H_fixed[0])))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         for (i = j1 = j2 = 0; i < i1->nNumberOfAtoms; i++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (i1->nNum_H_fixed[i] > i2->nNum_H_fixed[i])
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 j1++;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 if (i1->nNum_H_fixed[i] < i2->nNum_H_fixed[i])
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     j2++;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             ret = (j1 && j2) ? 20 : j1 ? 18 : j2 ? 19 : 0;
    // INCHI✔️❌:             if (ret)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret; /* 20 => Diff: NotEql Fixed-H */
    // INCHI✔️❌:                 /* 19 => Diff: Missed Fixed-H (i1 has less) */
    // INCHI✔️❌:                 /* 18 => Diff: Extra Fixed-H  (i1 has more) */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i1->lenConnTable != i2->lenConnTable)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 8; /* Diff: Connections length */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (i1->lenConnTable > 0 && memcmp(i1->nConnTable, i2->nConnTable, i1->lenConnTable * sizeof(i1->nConnTable[0])))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 9; /* Diff: Connections */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* output special cases: different number of t-groups, different sizes of t-groups, different endpoints */
    // INCHI✔️❌:     if (i1->lenTautomer != i2->lenTautomer && (i1->lenTautomer > 1 || i2->lenTautomer > 1))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 10; /* Diff: Mobile groups length */ /* in isotopic or deprotonated cases i1->lenTautomer == 1 && i1->nTautomer[0] = 0 */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if ((i1->lenTautomer > 1 && i2->lenTautomer > 1) &&
    // INCHI✔️❌:         memcmp(i1->nTautomer, i2->nTautomer, i1->lenTautomer * sizeof(i1->nTautomer[0])))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 11; /* Diff: Mobile groups */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i1->nNumberOfIsotopicAtoms != i2->nNumberOfIsotopicAtoms)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 12; /* Diff: Isotopic atoms number */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (i1->nNumberOfIsotopicAtoms > 0 && memcmp(i1->IsotopicAtom, i2->IsotopicAtom, i1->nNumberOfIsotopicAtoms * sizeof(i1->IsotopicAtom[0])))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 13; /* Diff: Isotopic atoms */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (i1->nTotalCharge != i2->nTotalCharge)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 14; /* Diff: Charge */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         if ( i1->nNumberOfIsotopicTGroups != i2->nNumberOfIsotopicTGroups )
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 14;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if ( i1->nNumberOfIsotopicTGroups > 0 && memcmp( i1->IsotopicTGroup, i2->IsotopicTGroup, i1->nNumberOfIsotopicTGroups*sizeof(i1->IsotopicTGroup[0]) ) )
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 15;
    // INCHI✔️❌:             }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (a1 && a2)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (a1->nNumRemovedProtons != a2->nNumRemovedProtons)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 16; /* Diff: Number of removed protons */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (memcmp(a1->nNumRemovedIsotopicH, a2->nNumRemovedIsotopicH, sizeof(a1->nNumRemovedIsotopicH)))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 17; /* Diff: Removed isotopic H */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         if ( i1->nPossibleLocationsOfIsotopicH && i2->nPossibleLocationsOfIsotopicH ) {
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ( i1->nPossibleLocationsOfIsotopicH[0] != i2->nPossibleLocationsOfIsotopicH[0] ||
    // INCHI✔️❌:                  memcmp(i1->nPossibleLocationsOfIsotopicH, i2->nPossibleLocationsOfIsotopicH,
    // INCHI✔️❌:                         sizeof(i1->nPossibleLocationsOfIsotopicH[0])*i1->nPossibleLocationsOfIsotopicH[0]) )
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 18;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         } else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ( !i1->nPossibleLocationsOfIsotopicH != !i2->nPossibleLocationsOfIsotopicH ) {
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 19;}
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* ret = 20..31 => 40..51 */
    // INCHI✔️❌:     if ((ret = CompareReversedStereoINChI(i1->Stereo, i2->Stereo))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret + 20;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* ret = 40..51 => 60..71 */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!i2->StereoIsotopic && i2->Stereo && i1->StereoIsotopic &&
    // INCHI✔️❌:         0 < (i1->StereoIsotopic->nNumberOfStereoBonds + i1->StereoIsotopic->nNumberOfStereoCenters) &&
    // INCHI✔️❌:         0 == CompareReversedStereoINChI(i1->StereoIsotopic, i2->Stereo))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* InChI from reversed structure does not contain fully duplicated isotopic stereo */
    // INCHI✔️❌:         ;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((ret = CompareReversedStereoINChI(i1->StereoIsotopic, i2->StereoIsotopic))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return ret + 40;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompareReversedINChI

    let (Some(i1), Some(i2)) = (i1, i2) else {
        return Ok(if i1.is_none() && i2.is_none() { 0 } else { 1 });
    };
    if i1.nErrorCode == i2.nErrorCode {
        if i1.nErrorCode != 0 {
            return Ok(0);
        }
    } else {
        return Ok(2);
    }
    if i1.bDeleted != i2.bDeleted {
        return Ok(1);
    }
    if i1.nNumberOfAtoms != i2.nNumberOfAtoms {
        return Ok(3);
    }
    if i1.nNumberOfAtoms > 0 {
        let count = usize::try_from(i1.nNumberOfAtoms)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let atoms1 = heap
            .slice(i1.nAtom.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let atoms2 = heap
            .slice(i2.nAtom.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if atoms1 != atoms2 {
            return Ok(4);
        }
        let formula1 = heap.slice(i1.szHillFormula.as_const())?;
        let formula2 = heap.slice(i2.szHillFormula.as_const())?;
        let end1 = formula1
            .iter()
            .position(|value| *value == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        let end2 = formula2
            .iter()
            .position(|value| *value == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        if formula1[..end1] != formula2[..end2] {
            return Ok(7);
        }
        let hydrogens1 = heap
            .slice(i1.nNum_H.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let hydrogens2 = heap
            .slice(i2.nNum_H.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if hydrogens1 != hydrogens2 {
            return Ok(if i1.lenConnTable > 1 || i2.lenConnTable > 1 {
                5
            } else {
                6
            });
        }
        if !i1.nNum_H_fixed.is_null() || !i2.nNum_H_fixed.is_null() {
            let fixed1 = if i1.nNum_H_fixed.is_null() {
                None
            } else {
                Some(
                    heap.slice(i1.nNum_H_fixed.as_const())?
                        .get(..count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                )
            };
            let fixed2 = if i2.nNum_H_fixed.is_null() {
                None
            } else {
                Some(
                    heap.slice(i2.nNum_H_fixed.as_const())?
                        .get(..count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                )
            };
            let has1 = fixed1.map_or(0, |values| {
                values.iter().filter(|value| **value != 0).count()
            });
            let has2 = fixed2.map_or(0, |values| {
                values.iter().filter(|value| **value != 0).count()
            });
            if has1 != 0 && has2 == 0 {
                return Ok(18);
            }
            if has1 == 0 && has2 != 0 {
                return Ok(19);
            }
            if has1 != 0 && has2 != 0 {
                let fixed1 = fixed1.ok_or(SourceHeapError::NullPointer)?;
                let fixed2 = fixed2.ok_or(SourceHeapError::NullPointer)?;
                if fixed1 != fixed2 {
                    let mut greater = false;
                    let mut less = false;
                    for (&left, &right) in fixed1.iter().zip(fixed2) {
                        greater |= left > right;
                        less |= left < right;
                    }
                    let ret = if greater && less {
                        20
                    } else if greater {
                        18
                    } else if less {
                        19
                    } else {
                        0
                    };
                    if ret != 0 {
                        return Ok(ret);
                    }
                }
            }
        }
    }
    if i1.lenConnTable != i2.lenConnTable {
        return Ok(8);
    }
    if i1.lenConnTable > 0 {
        let count =
            usize::try_from(i1.lenConnTable).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let connections1 = heap
            .slice(i1.nConnTable.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let connections2 = heap
            .slice(i2.nConnTable.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if connections1 != connections2 {
            return Ok(9);
        }
    }
    if i1.lenTautomer != i2.lenTautomer && (i1.lenTautomer > 1 || i2.lenTautomer > 1) {
        return Ok(10);
    }
    if i1.lenTautomer > 1 && i2.lenTautomer > 1 {
        let count =
            usize::try_from(i1.lenTautomer).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let tautomer1 = heap
            .slice(i1.nTautomer.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let tautomer2 = heap
            .slice(i2.nTautomer.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if tautomer1 != tautomer2 {
            return Ok(11);
        }
    }
    if i1.nNumberOfIsotopicAtoms != i2.nNumberOfIsotopicAtoms {
        return Ok(12);
    }
    if i1.nNumberOfIsotopicAtoms > 0 {
        let count = usize::try_from(i1.nNumberOfIsotopicAtoms)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let isotopic_atoms1 = heap
            .slice(i1.IsotopicAtom.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let isotopic_atoms2 = heap
            .slice(i2.IsotopicAtom.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if isotopic_atoms1 != isotopic_atoms2 {
            return Ok(13);
        }
    }
    if i1.nTotalCharge != i2.nTotalCharge {
        return Ok(14);
    }
    if let (Some(a1), Some(a2)) = (a1, a2) {
        if a1.nNumRemovedProtons != a2.nNumRemovedProtons {
            return Ok(16);
        }
        if a1.nNumRemovedIsotopicH != a2.nNumRemovedIsotopicH {
            return Ok(17);
        }
    }
    let stereo1 = if i1.Stereo.is_null() {
        None
    } else {
        Some(&heap.slice(i1.Stereo.as_const())?[0])
    };
    let stereo2 = if i2.Stereo.is_null() {
        None
    } else {
        Some(&heap.slice(i2.Stereo.as_const())?[0])
    };
    let ret = CompareReversedStereoINChI(heap, stereo1, stereo2)?;
    if ret != 0 {
        return Ok(ret + 20);
    }
    let isotopic1 = if i1.StereoIsotopic.is_null() {
        None
    } else {
        Some(&heap.slice(i1.StereoIsotopic.as_const())?[0])
    };
    let isotopic2 = if i2.StereoIsotopic.is_null() {
        None
    } else {
        Some(&heap.slice(i2.StereoIsotopic.as_const())?[0])
    };
    if isotopic2.is_none()
        && stereo2.is_some()
        && isotopic1.is_some()
        && isotopic1
            .is_some_and(|stereo| stereo.nNumberOfStereoBonds + stereo.nNumberOfStereoCenters > 0)
        && CompareReversedStereoINChI(heap, isotopic1, stereo2)? == 0
    {
        return Ok(0);
    }
    let ret = CompareReversedStereoINChI(heap, isotopic1, isotopic2)?;
    Ok(if ret != 0 { ret + 40 } else { 0 })
}

#[allow(non_snake_case)]
pub(crate) fn GetInpStructErrorType(
    input_parameters: Option<&INPUT_PARMS>,
    error: i32,
    structure_error: Option<&[i8]>,
    num_input_atoms: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2480 GetInpStructErrorType
    // INCHI✔️❌: int GetInpStructErrorType(INPUT_PARMS* ip,
    // INCHI✔️❌:     int err,
    // INCHI✔️❌:     char* pStrErrStruct,
    // INCHI✔️❌:     int num_inp_atoms)
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (err && err == 9)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return _IS_ERROR; /*  sdfile bypassed to $$$$ */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (err && err < 30)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return _IS_FATAL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (num_inp_atoms <= 0 || err)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (98 == err && 0 == num_inp_atoms && ip->bAllowEmptyStructure)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return _IS_WARNING;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return _IS_ERROR;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (pStrErrStruct[0])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return _IS_WARNING;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return _IS_OKAY;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: GetInpStructErrorType

    if error != 0 && error == 9 {
        return Ok(_IS_ERROR as i32);
    }
    if error != 0 && error < 30 {
        return Ok(_IS_FATAL as i32);
    }
    if num_input_atoms <= 0 || error != 0 {
        if error == 98
            && num_input_atoms == 0
            && input_parameters
                .ok_or(SourceHeapError::NullPointer)?
                .bAllowEmptyStructure
                != 0
        {
            return Ok(_IS_WARNING as i32);
        }
        return Ok(_IS_ERROR as i32);
    }
    let first = structure_error
        .and_then(|error| error.first())
        .ok_or(SourceHeapError::NullPointer)?;
    if *first != 0 {
        return Ok(_IS_WARNING as i32);
    }
    Ok(_IS_OKAY as i32)
}

pub(crate) fn mystrrev(
    heap: &mut SourceHeap,
    string: SourceMutPointer<i8>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2090 mystrrev
    // INCHI✔️✔️: void mystrrev(char* p)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     char c, * q = p;
    // INCHI✔️✔️:     while (*q++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     q -= 2; /*  pointer to the last character */
    // INCHI✔️✔️:     while (p < q)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         c = *q;  /*  swap */
    // INCHI✔️✔️:         *q-- = *p;
    // INCHI✔️✔️:         *p++ = c;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: mystrrev

    let bytes = heap.slice_mut(string)?;
    let length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    bytes[..length].reverse();
    Ok(())
}

struct OrderStruct<'a> {
    dfs_number: &'a [AT_NUMB],
    number_of_descendants: &'a [AT_NUMB],
    current_atom: i32,
}

#[allow(non_snake_case)]
fn CompareDfsDescendants4CT(
    neighbor1: AT_RANK,
    neighbor2: AT_RANK,
    order: &OrderStruct<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2119 CompareDfsDescendants4CT
    // INCHI✔️✔️: static int CompareDfsDescendants4CT(const void* a1, const void* a2, void* p)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     OrderStruct* os = (OrderStruct*)p;
    // INCHI✔️✔️:     int neigh1 = (int)*(const AT_RANK*)a1;
    // INCHI✔️✔️:     int neigh2 = (int)*(const AT_RANK*)a2;
    // INCHI✔️✔️:     if (neigh1 > MAX_ATOMS)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (neigh2 > MAX_ATOMS)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (neigh2 > MAX_ATOMS)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         AT_RANK nCurDfsNumber = os->m_gDfs4CT_nDfsNumber[os->m_gDfs4CT_nCurrentAtom];
    // INCHI✔️✔️:         int nDesc1 = nCurDfsNumber > os->m_gDfs4CT_nDfsNumber[neigh1] ?
    // INCHI✔️✔️:             0 : (int)os->m_gDfs4CT_nNumDescendants[neigh1];
    // INCHI✔️✔️:         int nDesc2 = nCurDfsNumber > os->m_gDfs4CT_nDfsNumber[neigh2] ?
    // INCHI✔️✔️:             0 : (int)os->m_gDfs4CT_nNumDescendants[neigh2];
    // INCHI✔️✔️:         int ret;
    // INCHI✔️✔️:         if ((ret = nDesc1 - nDesc2)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return ret;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         return  (int)neigh1 - (int)neigh2; /*  canon. numbers difference */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: CompareDfsDescendants4CT

    if u32::from(neighbor1) > MAX_ATOMS {
        return Ok(if u32::from(neighbor2) > MAX_ATOMS {
            0
        } else {
            1
        });
    }
    if u32::from(neighbor2) > MAX_ATOMS {
        return Ok(-1);
    }
    let current =
        usize::try_from(order.current_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let current_dfs = *order
        .dfs_number
        .get(current)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let first = usize::from(neighbor1);
    let second = usize::from(neighbor2);
    let descendants1 = if current_dfs
        > *order
            .dfs_number
            .get(first)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
    {
        0
    } else {
        i32::from(
            *order
                .number_of_descendants
                .get(first)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
    };
    let descendants2 = if current_dfs
        > *order
            .dfs_number
            .get(second)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
    {
        0
    } else {
        i32::from(
            *order
                .number_of_descendants
                .get(second)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
    };
    let result = descendants1 - descendants2;
    if result != 0 {
        return Ok(result);
    }
    Ok(i32::from(neighbor1) - i32::from(neighbor2))
}

#[allow(non_snake_case)]
pub(crate) fn GetDfsOrder4CT(
    heap: &mut SourceHeap,
    linear_ct: SourceMutPointer<AT_NUMB>,
    length_ct: i32,
    number_hydrogens: SourceConstPointer<i8>,
    number_of_atoms: i32,
    ct_mode: i32,
) -> Result<SourceMutPointer<AT_NUMB>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2154 GetDfsOrder4CT
    // INCHI✔️✔️: AT_NUMB* GetDfsOrder4CT(CANON_GLOBALS* pCG,
    // INCHI✔️✔️:     AT_NUMB* LinearCT,
    // INCHI✔️✔️:     int nLenCT,
    // INCHI✔️✔️:     S_CHAR* nNum_H,
    // INCHI✔️✔️:     int num_atoms,
    // INCHI✔️✔️:     int nCtMode)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     AT_NUMB* nStackAtom = NULL;
    // INCHI✔️✔️:     int         nTopStackAtom = -1;
    // INCHI✔️✔️:     AT_NUMB* nNumDescendants = NULL; /*  number of descendants incl. closures and the atom itself */
    // INCHI✔️✔️:     AT_NUMB* nDfsNumber = NULL;
    // INCHI✔️✔️:     S_CHAR* cNeighNumb = NULL;
    // INCHI✔️✔️:     NEIGH_LIST* nl = NULL;
    // INCHI✔️✔️:     AT_NUMB     nDfs;
    // INCHI✔️✔️:     int         i, j, u, k, start, num_rings, nTotOutputStringLen;
    // INCHI✔️✔️:     AT_NUMB* nOutputString = NULL, cDelim;
    // INCHI✔️✔️:     int         bCtPredecessors = (nCtMode & CT_MODE_PREDECESSORS);
    // INCHI✔️✔️:     OrderStruct os;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*  allocate arrays */
    // INCHI✔️✔️:     nStackAtom = (AT_NUMB*)inchi_malloc(num_atoms * sizeof(nStackAtom[0]));
    // INCHI✔️✔️:     nNumDescendants = (AT_NUMB*)inchi_malloc(num_atoms * sizeof(nNumDescendants[0]));
    // INCHI✔️✔️:     nDfsNumber = (AT_NUMB*)inchi_malloc(num_atoms * sizeof(nDfsNumber[0]));
    // INCHI✔️✔️:     cNeighNumb = (S_CHAR*)inchi_malloc(num_atoms * sizeof(cNeighNumb[0]));
    // INCHI✔️✔️:     nl = CreateNeighListFromLinearCT(LinearCT, nLenCT, num_atoms);
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*  check allocation */
    // INCHI✔️✔️:     if (!nStackAtom || !nNumDescendants || !nDfsNumber || !cNeighNumb || !nl)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* ret = CT_OUT_OF_RAM; */ /*  program error */ /*   <BRKPT> */
    // INCHI✔️✔️:         goto exit_function;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (bCtPredecessors)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         start = 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /*  find DFS start vertex (atom) */
    // INCHI✔️✔️:         for (i = 1, start = 0; i < num_atoms; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (nl[i][0] < nl[start][0])
    // INCHI✔️✔️:             { /*  index = nRank-1 */
    // INCHI✔️✔️:                 start = i;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*
    // INCHI✔️✔️:       vertex information:
    // INCHI✔️✔️:         1. Number of (forward edges) + (back edges, first visit -- ring closures): nl[i][0]
    // INCHI✔️✔️:         2. Number of vertices traversed from this vertex, including the vertex:    nNumDescendants[i]
    // INCHI✔️✔️:         3. Each edge information:
    // INCHI✔️✔️:            a. forward edge (0) or back edge (1) indicator: nDfsNumber[i] > nDfsNumber[neigh]
    // INCHI✔️✔️:            b. neighbor at another end of the edge neigh = nl[i][k+1], k < i
    // INCHI✔️✔️:
    // INCHI✔️✔️:         Total per edge: 2 + 2*(number of edges)
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* DFS initiation */
    // INCHI✔️✔️:     u = start; /* start atom */
    // INCHI✔️✔️:     nDfs = 0;
    // INCHI✔️✔️:     nTopStackAtom = -1;
    // INCHI✔️✔️:     memset(nDfsNumber, 0, num_atoms * sizeof(nDfsNumber[0])); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:     memset(nNumDescendants, 0, num_atoms * sizeof(nNumDescendants[0])); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:     memset(cNeighNumb, 0, num_atoms * sizeof(cNeighNumb[0])); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:     /*  push the start atom on the stack */
    // INCHI✔️✔️:     nDfsNumber[u] = ++nDfs;
    // INCHI✔️✔️:     if (bCtPredecessors)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nNumDescendants[u] = 0; /* atom #1 has no predecessor */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nNumDescendants[u] = 1; /* count itself as a descendant */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
    // INCHI✔️✔️:     /* nNumStartChildren = 0; */
    // INCHI✔️✔️:     num_rings = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* DFS */
    // INCHI✔️✔️:     do
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* advance */
    // INCHI✔️✔️:         while (i = (int)nStackAtom[nTopStackAtom], j = (int)cNeighNumb[i] + 1, (int)nl[i][0] >= j)
    // INCHI✔️✔️:             /*while ( (int)nl[i=nStackAtom[nTopStackAtom]][0] >= (j = (int)cNeighNumb[i]+1) )*/
    // INCHI✔️✔️:
    // INCHI✔️✔️:             /* replaced due to missing sequence point; undefined behavior, pointed by Geoffrey Hutchison */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             cNeighNumb[i]++;
    // INCHI✔️✔️:             u = (int)nl[i][j]; /*  jth neighbor of the vertex i */
    // INCHI✔️✔️:             if (!nDfsNumber[u])
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 /* tree edge, 1st visit -- advance */
    // INCHI✔️✔️:                 /* put unexplored vertex u on the stack for further examination */
    // INCHI✔️✔️:                 nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
    // INCHI✔️✔️:                 nDfsNumber[u] = ++nDfs;
    // INCHI✔️✔️:                 if (bCtPredecessors)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     nNumDescendants[u] = i + 1; /* predecessor's rank */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     nNumDescendants[u]++; /* count atom u as its descendant */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (nTopStackAtom && u != (int)nStackAtom[nTopStackAtom - 1] &&
    // INCHI✔️✔️:                     /* back edge: u is not a predecessor of i */
    // INCHI✔️✔️:                     nDfsNumber[u] < nDfsNumber[i])
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /* Back edge, 1st visit: u is an ancestor of i (ring closure) */
    // INCHI✔️✔️:                     if (!bCtPredecessors)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         nNumDescendants[i]++; /* count closures as descendants */
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     num_rings++;          /* count ring closures */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     nl[i][j] = MAX_ATOMS + 1; /* back edge, 2nd visit: mark as deleted */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         cNeighNumb[i] = 0; /* all neighbors of the ith atom have been
    // INCHI✔️✔️:                               traversed; resore the neighbor counter */
    // INCHI✔️✔️:                               /* back up */
    // INCHI✔️✔️:         if (!bCtPredecessors && nTopStackAtom /* that is, i != start */)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             u = (int)nStackAtom[nTopStackAtom - 1]; /* predecessor of i */
    // INCHI✔️✔️:             nNumDescendants[u] += nNumDescendants[i]; /* add descendants */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     } while (--nTopStackAtom >= 0);
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Sort the neighbors in ascending order so that:
    // INCHI✔️✔️:        primary key   = number of descendants in the DFS tree; closure neighbor is 0
    // INCHI✔️✔️:        secondary key = canonical number (here vertex number = canonical number - 1)
    // INCHI✔️✔️:      */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     os.m_gDfs4CT_nDfsNumber = nDfsNumber;
    // INCHI✔️✔️:     os.m_gDfs4CT_nNumDescendants = nNumDescendants;
    // INCHI✔️✔️:     os.m_gDfs4CT_nCurrentAtom = -1;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* sorting; deleted will be the last neighbors */
    // INCHI✔️✔️:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (nl[i][0] > 1)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             os.m_gDfs4CT_nCurrentAtom = i;
    // INCHI✔️✔️:             insertions_sort(&os, &nl[i][1], nl[i][0], sizeof(nl[i][1]), CompareDfsDescendants4CT);
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         /* reduce number of neighbors to exclude deleted */
    // INCHI✔️✔️:         for (k = 0; k < nl[i][0] && nl[i][k + 1] <= MAX_ATOMS; k++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             ;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         nl[i][0] = k;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     nTotOutputStringLen = 3 * (num_atoms + num_rings + 1); /*  last 3 elements are a 'zero termination' */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (bCtPredecessors)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if ((nOutputString = (AT_RANK*)inchi_calloc(nTotOutputStringLen, sizeof(nOutputString[0])))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             cDelim = '-';
    // INCHI✔️✔️:             for (u = 0, k = -3; u < num_atoms; u++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 k += 3;
    // INCHI✔️✔️:                 if (k + 6 > nTotOutputStringLen)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     goto exit_error;  /* program error */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 nOutputString[k] = nNumDescendants[u] ? nNumDescendants[u] : MAX_ATOMS + 1;
    // INCHI✔️✔️:                 nOutputString[k + 1] = nNum_H ? 16 + nNum_H[u] : 0;
    // INCHI✔️✔️:                 nOutputString[k + 2] = k ? ',' : '\0';
    // INCHI✔️✔️:                 for (j = 1; j <= nl[u][0] && nDfsNumber[u] > nDfsNumber[i = nl[u][j]]; j++)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /* closures */
    // INCHI✔️✔️:                     k += 3;
    // INCHI✔️✔️:                     if (k + 6 > nTotOutputStringLen)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         goto exit_error;  /* program error */
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     nOutputString[k] = i + 1;  /* closure */
    // INCHI✔️✔️:                     nOutputString[k + 1] = 0;
    // INCHI✔️✔️:                     nOutputString[k + 2] = cDelim;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (nNumDescendants)
    // INCHI✔️✔️:         {  /* do not need anymore */
    // INCHI✔️✔️:             inchi_free(nNumDescendants);
    // INCHI✔️✔️:             nNumDescendants = NULL;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /*
    // INCHI✔️✔️:             the output string contains:
    // INCHI✔️✔️:               (num_atoms) atoms for the DFS (spanning) tree
    // INCHI✔️✔️:               (num_atoms-1) delimiters for the DFS (spanning) tree
    // INCHI✔️✔️:               1 character for each atom that has 1 terminal hydrogen atoms
    // INCHI✔️✔️:               2 characters  for each atom that has 2-9 terminal hydrogen atoms
    // INCHI✔️✔️:               3 characters  for each atom that has 10-99 terminal hydrogen atoms, etc.
    // INCHI✔️✔️:               (num_rings) atoms for the ring closures
    // INCHI✔️✔️:               (num_rings) delimiters for the ring closures
    // INCHI✔️✔️:         */
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if ((nOutputString = (AT_RANK*)inchi_calloc(nTotOutputStringLen, sizeof(nOutputString[0])))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             u = start; /*  start atom */
    // INCHI✔️✔️:             nTopStackAtom = -1;
    // INCHI✔️✔️:             memset(cNeighNumb, 0, num_atoms * sizeof(cNeighNumb[0])); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:             /*  push the start atom on the stack */
    // INCHI✔️✔️:             nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
    // INCHI✔️✔️:             /*  output the starting atom */
    // INCHI✔️✔️:             k = 0;
    // INCHI✔️✔️:             nOutputString[k] = u + 1;
    // INCHI✔️✔️:             nOutputString[k + 1] = nNum_H ? 16 + nNum_H[u] : 0;
    // INCHI✔️✔️:             nOutputString[k + 2] = '\0';
    // INCHI✔️✔️:
    // INCHI✔️✔️:             do
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 /* advance */
    // INCHI✔️✔️:                 while (i = (int)nStackAtom[nTopStackAtom], j = (int)cNeighNumb[i] + 1, i < num_atoms && (int)nl[i][0] >= j) /* djb-rwth: additional condition to avoid buffer overruns */
    // INCHI✔️✔️:                     /*while ( (int)nl[i=nStackAtom[nTopStackAtom]][0] >= (j = (int)cNeighNumb[i]+1) )*/
    // INCHI✔️✔️:                     /* replaced due to missing sequence point; undefined behavior, reported by Geoffrey Hutchison */
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     k += 3;
    // INCHI✔️✔️:                     if (k + 6 > nTotOutputStringLen)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         goto exit_error;  /* program error */
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     cNeighNumb[i]++;
    // INCHI✔️✔️:                     u = (int)nl[i][j]; /* neighbor */
    // INCHI✔️✔️:
    // INCHI✔️✔️:                     /* output neighbor's canonical number */
    // INCHI✔️✔️:                     nOutputString[k] = u + 1;
    // INCHI✔️✔️:
    // INCHI✔️✔️:                     if (nDfsNumber[u] > nDfsNumber[i])
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         /* tree edge, 1st visit -- advance */
    // INCHI✔️✔️:                         /* put 'unexplored' vertex u on the stack */
    // INCHI✔️✔️:                         nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
    // INCHI✔️✔️:
    // INCHI✔️✔️:                         /* output neighbor's number of H */
    // INCHI✔️✔️:                         nOutputString[k + 1] = nNum_H ? 16 + nNum_H[u] : 0;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     else
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         nOutputString[k + 1] = 0;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     /* output a delimiter preceding the neighbor */
    // INCHI✔️✔️:                     if (1 < nl[i][0])
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         if (j == 1)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             cDelim = '(';
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         else
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             if (j == nl[i][0])
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 cDelim = ')';
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                             else
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 cDelim = ',';
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     else
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         cDelim = '-';
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     nOutputString[k + 2] = cDelim;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 if ((i >= 0) && (i < num_atoms)) /* djb-rwth: fixing coverity ID #499483 */
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     cNeighNumb[i] = 0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 /* back up: nothing else to do */
    // INCHI✔️✔️:             } while (--nTopStackAtom >= 0);
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     goto exit_function;
    // INCHI✔️✔️:
    // INCHI✔️✔️: exit_error:
    // INCHI✔️✔️:     if (nOutputString)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inchi_free(nOutputString);
    // INCHI✔️✔️:         nOutputString = NULL;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️: exit_function:
    // INCHI✔️✔️:     if (nStackAtom)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inchi_free(nStackAtom);
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (nNumDescendants)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inchi_free(nNumDescendants);
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (nDfsNumber)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inchi_free(nDfsNumber);
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (cNeighNumb)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inchi_free(cNeighNumb);
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (nl)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         FreeNeighList(nl);
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return nOutputString;
    // INCHI✔️✔️: }
    // INCHI✔️✔️:
    // END INCHI C FUNCTION: GetDfsOrder4CT

    let atom_count =
        usize::try_from(number_of_atoms).map_err(|_| SourceHeapError::UnsupportedSourceBehavior)?;
    if atom_count == 0 {
        return Err(SourceHeapError::UnsupportedSourceBehavior);
    }
    let allocate_u16 = |heap: &mut SourceHeap| match inchi_calloc::<u16>(
        heap,
        atom_count as u64,
        std::mem::size_of::<u16>() as u64,
    ) {
        Ok(pointer) => Ok(pointer),
        Err(SourceHeapError::AllocationFailed) => Ok(SourceMutPointer::null()),
        Err(error) => Err(error),
    };
    let stack = allocate_u16(heap)?;
    let mut number_of_descendants = allocate_u16(heap)?;
    let dfs_number = allocate_u16(heap)?;
    let neighbor_number = match inchi_calloc::<i8>(heap, atom_count as u64, 1) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
        Err(error) => return Err(error),
    };
    let neighbor_lists =
        CreateNeighListFromLinearCT(heap, linear_ct.as_const(), length_ct, number_of_atoms)?;
    let predecessors = ct_mode & crate::source_types::CT_MODE_PREDECESSORS as i32;

    let computation = (|| -> Result<SourceMutPointer<AT_NUMB>, SourceHeapError> {
        if stack.is_null()
            || number_of_descendants.is_null()
            || dfs_number.is_null()
            || neighbor_number.is_null()
            || neighbor_lists.is_null()
        {
            return Ok(SourceMutPointer::null());
        }

        let mut start = 0_usize;
        if predecessors == 0 {
            for atom in 1..atom_count {
                let pointers = heap.slice(neighbor_lists.as_const())?;
                let valence = heap.slice(pointers[atom].as_const())?[0];
                let start_valence = heap.slice(pointers[start].as_const())?[0];
                if valence < start_valence {
                    start = atom;
                }
            }
        }

        heap.slice_mut(dfs_number)?.fill(0);
        heap.slice_mut(number_of_descendants)?.fill(0);
        heap.slice_mut(neighbor_number)?.fill(0);
        let mut dfs = 1_u16;
        heap.slice_mut(dfs_number)?[start] = dfs;
        heap.slice_mut(number_of_descendants)?[start] = if predecessors != 0 { 0 } else { 1 };
        heap.slice_mut(stack)?[0] = start as u16;
        let mut stack_top = 0_i32;
        let mut number_of_rings = 0_i32;

        loop {
            let mut atom = usize::from(heap.slice(stack.as_const())?[stack_top as usize]);
            loop {
                let neighbor_index = i32::from(heap.slice(neighbor_number.as_const())?[atom]) + 1;
                let list_pointer = heap.slice(neighbor_lists.as_const())?[atom];
                if i32::from(heap.slice(list_pointer.as_const())?[0]) < neighbor_index {
                    break;
                }
                heap.slice_mut(neighbor_number)?[atom] =
                    heap.slice(neighbor_number.as_const())?[atom].wrapping_add(1);
                let adjacent =
                    usize::from(heap.slice(list_pointer.as_const())?[neighbor_index as usize]);
                if heap.slice(dfs_number.as_const())?[adjacent] == 0 {
                    stack_top += 1;
                    heap.slice_mut(stack)?[stack_top as usize] = adjacent as u16;
                    dfs = dfs.wrapping_add(1);
                    heap.slice_mut(dfs_number)?[adjacent] = dfs;
                    if predecessors != 0 {
                        heap.slice_mut(number_of_descendants)?[adjacent] = (atom + 1) as u16;
                    } else {
                        heap.slice_mut(number_of_descendants)?[adjacent] =
                            heap.slice(number_of_descendants.as_const())?[adjacent].wrapping_add(1);
                    }
                } else if stack_top != 0
                    && adjacent
                        != usize::from(heap.slice(stack.as_const())?[(stack_top - 1) as usize])
                    && heap.slice(dfs_number.as_const())?[adjacent]
                        < heap.slice(dfs_number.as_const())?[atom]
                {
                    if predecessors == 0 {
                        heap.slice_mut(number_of_descendants)?[atom] =
                            heap.slice(number_of_descendants.as_const())?[atom].wrapping_add(1);
                    }
                    number_of_rings = number_of_rings.wrapping_add(1);
                } else {
                    heap.slice_mut(list_pointer)?[neighbor_index as usize] = (MAX_ATOMS + 1) as u16;
                }
                atom = usize::from(heap.slice(stack.as_const())?[stack_top as usize]);
            }
            heap.slice_mut(neighbor_number)?[atom] = 0;
            if predecessors == 0 && stack_top != 0 {
                let parent = usize::from(heap.slice(stack.as_const())?[(stack_top - 1) as usize]);
                heap.slice_mut(number_of_descendants)?[parent] = heap
                    .slice(number_of_descendants.as_const())?[parent]
                    .wrapping_add(heap.slice(number_of_descendants.as_const())?[atom]);
            }
            stack_top -= 1;
            if stack_top < 0 {
                break;
            }
        }

        for atom in 0..atom_count {
            let list_pointer = heap.slice(neighbor_lists.as_const())?[atom];
            let neighbor_count = usize::from(heap.slice(list_pointer.as_const())?[0]);
            if neighbor_count > 1 {
                heap.with_slice_mut_and_heap(list_pointer.offset(1)?, |neighbors, read_heap| {
                    let records = bytemuck::cast_slice_mut::<u16, u8>(
                        neighbors
                            .get_mut(..neighbor_count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    let order = OrderStruct {
                        dfs_number: read_heap.slice(dfs_number.as_const())?,
                        number_of_descendants: read_heap.slice(number_of_descendants.as_const())?,
                        current_atom: atom as i32,
                    };
                    insertions_sort(records, neighbor_count, 2, &mut |first, second| {
                        CompareDfsDescendants4CT(
                            u16::from_ne_bytes([first[0], first[1]]),
                            u16::from_ne_bytes([second[0], second[1]]),
                            &order,
                        )
                    })?;
                    Ok(())
                })?;
            }
            let list = heap.slice_mut(list_pointer)?;
            let mut retained = 0_usize;
            while retained < neighbor_count && u32::from(list[retained + 1]) <= MAX_ATOMS {
                retained += 1;
            }
            list[0] = retained as u16;
        }

        let total_length = number_of_atoms
            .checked_add(number_of_rings)
            .and_then(|value| value.checked_add(1))
            .and_then(|value| value.checked_mul(3))
            .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
        let output_length = usize::try_from(total_length)
            .map_err(|_| SourceHeapError::UnsupportedSourceBehavior)?;
        let output = match inchi_calloc::<u16>(heap, output_length as u64, 2) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => return Ok(SourceMutPointer::null()),
            Err(error) => return Err(error),
        };

        if predecessors != 0 {
            let mut position = 0_usize;
            for atom in 0..atom_count {
                if position + 6 > output_length {
                    inchi_free(heap, output)?;
                    return Ok(SourceMutPointer::null());
                }
                let predecessor = heap.slice(number_of_descendants.as_const())?[atom];
                heap.slice_mut(output)?[position] = if predecessor != 0 {
                    predecessor
                } else {
                    (MAX_ATOMS + 1) as u16
                };
                heap.slice_mut(output)?[position + 1] = if number_hydrogens.is_null() {
                    0
                } else {
                    (16_i32 + i32::from(heap.slice(number_hydrogens)?[atom])) as u16
                };
                heap.slice_mut(output)?[position + 2] = if position != 0 { b',' as u16 } else { 0 };
                let list_pointer = heap.slice(neighbor_lists.as_const())?[atom];
                let list_count = usize::from(heap.slice(list_pointer.as_const())?[0]);
                let mut list_index = 1_usize;
                while list_index <= list_count {
                    let adjacent = heap.slice(list_pointer.as_const())?[list_index];
                    if heap.slice(dfs_number.as_const())?[atom]
                        <= heap.slice(dfs_number.as_const())?[usize::from(adjacent)]
                    {
                        break;
                    }
                    position += 3;
                    if position + 6 > output_length {
                        inchi_free(heap, output)?;
                        return Ok(SourceMutPointer::null());
                    }
                    heap.slice_mut(output)?[position] = adjacent.wrapping_add(1);
                    heap.slice_mut(output)?[position + 1] = 0;
                    heap.slice_mut(output)?[position + 2] = b'-' as u16;
                    list_index += 1;
                }
                position += 3;
            }
        } else {
            inchi_free(heap, number_of_descendants)?;
            number_of_descendants = SourceMutPointer::null();
            heap.slice_mut(neighbor_number)?.fill(0);
            heap.slice_mut(stack)?[0] = start as u16;
            stack_top = 0;
            let mut position = 0_usize;
            heap.slice_mut(output)?[0] = (start + 1) as u16;
            heap.slice_mut(output)?[1] = if number_hydrogens.is_null() {
                0
            } else {
                (16_i32 + i32::from(heap.slice(number_hydrogens)?[start])) as u16
            };
            heap.slice_mut(output)?[2] = 0;
            loop {
                let atom = usize::from(heap.slice(stack.as_const())?[stack_top as usize]);
                let list_pointer = heap.slice(neighbor_lists.as_const())?[atom];
                let list_count = usize::from(heap.slice(list_pointer.as_const())?[0]);
                let neighbor_index =
                    usize::try_from(i32::from(heap.slice(neighbor_number.as_const())?[atom]))
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?
                        + 1;
                if atom < atom_count && neighbor_index <= list_count {
                    position += 3;
                    if position + 6 > output_length {
                        inchi_free(heap, output)?;
                        return Ok(SourceMutPointer::null());
                    }
                    heap.slice_mut(neighbor_number)?[atom] =
                        heap.slice(neighbor_number.as_const())?[atom].wrapping_add(1);
                    let adjacent =
                        usize::from(heap.slice(list_pointer.as_const())?[neighbor_index]);
                    heap.slice_mut(output)?[position] = (adjacent + 1) as u16;
                    if heap.slice(dfs_number.as_const())?[adjacent]
                        > heap.slice(dfs_number.as_const())?[atom]
                    {
                        stack_top += 1;
                        heap.slice_mut(stack)?[stack_top as usize] = adjacent as u16;
                        heap.slice_mut(output)?[position + 1] = if number_hydrogens.is_null() {
                            0
                        } else {
                            (16_i32 + i32::from(heap.slice(number_hydrogens)?[adjacent])) as u16
                        };
                    } else {
                        heap.slice_mut(output)?[position + 1] = 0;
                    }
                    heap.slice_mut(output)?[position + 2] = if list_count > 1 {
                        if neighbor_index == 1 {
                            b'(' as u16
                        } else if neighbor_index == list_count {
                            b')' as u16
                        } else {
                            b',' as u16
                        }
                    } else {
                        b'-' as u16
                    };
                } else {
                    if atom < atom_count {
                        heap.slice_mut(neighbor_number)?[atom] = 0;
                    }
                    stack_top -= 1;
                    if stack_top < 0 {
                        break;
                    }
                }
            }
        }
        Ok(output)
    })();

    inchi_free(heap, stack)?;
    inchi_free(heap, number_of_descendants)?;
    inchi_free(heap, dfs_number)?;
    inchi_free(heap, neighbor_number)?;
    FreeNeighList(heap, neighbor_lists)?;
    computation
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_port__ichimake__mystrrev__line_2090() {
        for (input, expected) in [
            ("", ""),
            ("a", "a"),
            ("ab", "ba"),
            ("abc", "cba"),
            ("abcdef", "fedcba"),
        ] {
            let mut heap = SourceHeap::default();
            let mut bytes: Vec<i8> = input.bytes().map(|byte| byte as i8).collect();
            bytes.push(0);
            bytes.push(79);
            let string = heap.allocate_model_storage(bytes).unwrap();
            assert_eq!(mystrrev(&mut heap, string), Ok(()));
            let bytes = heap.slice(string.as_const()).unwrap();
            assert_eq!(
                &bytes[..expected.len()],
                &expected.bytes().map(|byte| byte as i8).collect::<Vec<_>>()
            );
            assert_eq!(bytes[expected.len()], 0);
            assert_eq!(bytes[expected.len() + 1], 79);
        }

        let mut heap = SourceHeap::default();
        let unterminated = heap.allocate_model_storage(vec![b'a' as i8]).unwrap();
        assert_eq!(
            mystrrev(&mut heap, unterminated),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(
            mystrrev(&mut heap, SourceMutPointer::null()),
            Err(SourceHeapError::NullPointer)
        );
    }
    use crate::source_types::{
        INChI_IsotopicAtom, INChI_IsotopicTGroup, tagMarkDiff_DIFV_NEQ2PRECED,
    };

    fn inchi_fixture(heap: &mut SourceHeap) -> INChI {
        INChI {
            nNumberOfAtoms: 2,
            szHillFormula: heap.allocate(vec![b'C' as i8, b'2' as i8, 0]).unwrap(),
            nAtom: heap.allocate(vec![6_u8, 6]).unwrap(),
            lenConnTable: 2,
            nConnTable: heap.allocate(vec![1_u16, 2]).unwrap(),
            lenTautomer: 2,
            nTautomer: heap.allocate(vec![2_u16, 1]).unwrap(),
            nNum_H: heap.allocate(vec![3_i8, 3]).unwrap(),
            nNumberOfIsotopicAtoms: 1,
            IsotopicAtom: heap
                .allocate(vec![INChI_IsotopicAtom {
                    nAtomNumber: 1,
                    ..INChI_IsotopicAtom::default()
                }])
                .unwrap(),
            ..INChI::default()
        }
    }

    #[test]
    fn source_port__ichimake__markunusedandemptylayers__line_1525() {
        let main = [9_i8; 11];
        let mut omitted = [main, [1_i8; 11], [1_i8; 11], [1_i8; 11]];
        assert_eq!(MarkUnusedAndEmptyLayers(&mut omitted), 0);
        assert_eq!(omitted[0], main);
        assert_eq!(omitted[1], [0; 11]);
        assert_eq!(omitted[2], [0; 11]);
        assert_eq!(omitted[3], [0; 11]);

        let mut exposed = [[0_i8; 11]; 4];
        exposed[0] = main;
        exposed[tagDiffINChILayers_DIFL_FI as usize][1] = tagMarkDiff_DIFV_NEQ2PRECED as i8;
        exposed[tagDiffINChILayers_DIFL_MI as usize][2] = tagMarkDiff_DIFV_IS_EMPTY as i8;
        assert_eq!(MarkUnusedAndEmptyLayers(&mut exposed), 0);
        assert_eq!(exposed[0], main);
        assert_eq!(
            exposed[tagDiffINChILayers_DIFL_FI as usize]
                [tagDiffINChISegments_DIFS_i_IATOMS as usize],
            tagMarkDiff_DIFV_IS_EMPTY as i8
        );
        assert_eq!(
            exposed[tagDiffINChILayers_DIFL_MI as usize]
                [tagDiffINChISegments_DIFS_i_IATOMS as usize],
            tagMarkDiff_DIFV_IS_EMPTY as i8
        );
        assert_eq!(
            exposed[tagDiffINChILayers_DIFL_F as usize]
                [tagDiffINChISegments_DIFS_f_FORMULA as usize],
            tagMarkDiff_DIFV_IS_EMPTY as i8
        );

        let mut preserved = [[0_i8; 11]; 4];
        preserved[0] = main;
        preserved[tagDiffINChILayers_DIFL_FI as usize]
            [tagDiffINChISegments_DIFS_i_IATOMS as usize] = tagMarkDiff_DIFV_OUTPUT_OMIT_F as i8;
        preserved[tagDiffINChILayers_DIFL_MI as usize]
            [tagDiffINChISegments_DIFS_i_IATOMS as usize] = tagMarkDiff_DIFV_OUTPUT_OMIT_F as i8;
        preserved[tagDiffINChILayers_DIFL_F as usize]
            [tagDiffINChISegments_DIFS_f_FORMULA as usize] = tagMarkDiff_DIFV_OUTPUT_OMIT_F as i8;
        let before = preserved;
        assert_eq!(MarkUnusedAndEmptyLayers(&mut preserved), 0);
        assert_eq!(preserved, before);
    }

    #[test]
    fn source_port__ichimake__compareinchistereo__line_1607() {
        let mut heap = SourceHeap::default();
        assert_eq!(CompareInchiStereo(&heap, None, 0, None, 0), Ok(0));
        let empty = INChI_Stereo::default();
        assert_eq!(CompareInchiStereo(&heap, Some(&empty), 0, None, 0), Ok(0));
        let nonempty = stereo_fixture(&mut heap, &[2], &[1], 1, &[3], &[4], &[1]);
        assert_eq!(
            CompareInchiStereo(&heap, None, 0, Some(&nonempty), 0),
            Ok(1)
        );
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, None, 0),
            Ok(-1)
        );

        let equal = stereo_fixture(&mut heap, &[2], &[1], 1, &[3], &[4], &[1]);
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, Some(&equal), 0),
            Ok(0)
        );
        let mut changed = stereo_fixture(&mut heap, &[2], &[1], 1, &[5], &[4], &[1]);
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, Some(&changed), 0),
            Ok(2)
        );
        changed = stereo_fixture(&mut heap, &[2], &[1], 1, &[3], &[6], &[1]);
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, Some(&changed), 0),
            Ok(2)
        );
        changed = stereo_fixture(&mut heap, &[2], &[1], 1, &[3], &[4], &[3]);
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, Some(&changed), 0),
            Ok(2)
        );
        changed = equal.clone();
        changed.nNumberOfStereoBonds = 2;
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, Some(&changed), 0),
            Ok(1)
        );

        changed = stereo_fixture(&mut heap, &[5], &[1], 1, &[3], &[4], &[1]);
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, Some(&changed), 0),
            Ok(3)
        );
        changed = stereo_fixture(&mut heap, &[2], &[3], 1, &[3], &[4], &[1]);
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, Some(&changed), 0),
            Ok(2)
        );
        changed = equal.clone();
        changed.nNumberOfStereoCenters = 2;
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, Some(&changed), 0),
            Ok(1)
        );

        changed = equal.clone();
        changed.nCompInv2Abs = -1;
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, Some(&changed), 0),
            Ok(1)
        );
        assert_eq!(
            CompareInchiStereo(
                &heap,
                Some(&nonempty),
                INCHI_MODE::from(INCHI_FLAG_REL_STEREO),
                Some(&changed),
                0,
            ),
            Ok(0)
        );
        assert_eq!(
            CompareInchiStereo(
                &heap,
                Some(&nonempty),
                0,
                Some(&changed),
                INCHI_MODE::from(INCHI_FLAG_RAC_STEREO),
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichimake__getelementandcount__line_173() {
        let mut heap = SourceHeap::default();
        let formula = heap
            .allocate(vec![b'C' as i8, b'2' as i8, b'H' as i8, b'6' as i8, 0])
            .unwrap();
        let mut output = [99_i8; 6];
        let mut cursor = formula.as_const();
        let mut count = -7;
        assert_eq!(
            GetElementAndCount(&mut heap, &mut cursor, &mut output, &mut count),
            Ok(1)
        );
        assert_eq!(count, 2);
        assert_eq!(cursor, formula.as_const().offset(2).unwrap());
        assert_eq!(output, [65, 0, 99, 99, 99, 99]);
        assert_eq!(
            GetElementAndCount(&mut heap, &mut cursor, &mut output, &mut count),
            Ok(1)
        );
        assert_eq!(count, 6);
        assert_eq!(cursor, formula.as_const().offset(4).unwrap());
        assert_eq!(output, [72, 0, 99, 99, 99, 99]);
        assert_eq!(
            GetElementAndCount(&mut heap, &mut cursor, &mut output, &mut count),
            Ok(0)
        );
        assert_eq!(count, 99999);
        assert_eq!(cursor, formula.as_const().offset(4).unwrap());
        assert_eq!(output, [90, 122, 122, 0, 99, 99]);

        let chlorine = heap
            .allocate(vec![
                b'C' as i8, b'l' as i8, b'1' as i8, b'2' as i8, b'X' as i8, 0,
            ])
            .unwrap();
        cursor = chlorine.as_const();
        count = 0;
        assert_eq!(
            GetElementAndCount(&mut heap, &mut cursor, &mut output, &mut count),
            Ok(1)
        );
        assert_eq!(count, 12);
        assert_eq!(cursor, chlorine.as_const().offset(4).unwrap());
        assert_eq!(&output[..4], &[67, 108, 0, 0]);

        let nitrogen = heap.allocate(vec![b'N' as i8, 0]).unwrap();
        cursor = nitrogen.as_const();
        assert_eq!(
            GetElementAndCount(&mut heap, &mut cursor, &mut output, &mut count),
            Ok(1)
        );
        assert_eq!(count, 1);
        assert_eq!(cursor, nitrogen.as_const().offset(1).unwrap());

        let invalid = heap.allocate(vec![b'c' as i8, b'2' as i8, 0]).unwrap();
        cursor = invalid.as_const();
        count = 77;
        output.fill(88);
        assert_eq!(
            GetElementAndCount(&mut heap, &mut cursor, &mut output, &mut count),
            Ok(-1)
        );
        assert_eq!(cursor, invalid.as_const());
        assert_eq!(count, 77);
        assert_eq!(output, [88; 6]);

        cursor = SourceConstPointer::null();
        count = 0;
        assert_eq!(
            GetElementAndCount(&mut heap, &mut cursor, &mut output, &mut count),
            Ok(0)
        );
        assert!(cursor.is_null());
        assert_eq!(count, 99999);
        assert_eq!(&output[..5], &[90, 122, 122, 0, 88]);

        let overflow = heap
            .allocate(
                b"N9223372036854775808X\0"
                    .iter()
                    .map(|byte| *byte as i8)
                    .collect(),
            )
            .unwrap();
        cursor = overflow.as_const();
        assert_eq!(
            GetElementAndCount(&mut heap, &mut cursor, &mut output, &mut count),
            Ok(1)
        );
        assert_eq!(count, -1);
        assert_eq!(cursor, overflow.as_const().offset(20).unwrap());
    }

    #[test]
    fn source_port__ichimake__comparehillformulas__line_241() {
        let mut heap = SourceHeap::default();
        let mut formula = |bytes: &[u8]| {
            let mut bytes: Vec<i8> = bytes.iter().map(|byte| *byte as i8).collect();
            bytes.push(0);
            heap.allocate(bytes).unwrap().as_const()
        };
        let c2h6 = formula(b"C2H6");
        let c3h6 = formula(b"C3H6");
        let c2n = formula(b"C2N");
        let c2o = formula(b"C2O");
        let oxygen = formula(b"O");
        let empty = formula(b"");
        let invalid = formula(b"c2");
        let max_count = formula(b"N2147483647");
        let wrapped_count = formula(b"N2147483648");

        assert_eq!(CompareHillFormulas(&mut heap, c2h6, c2h6), Ok(0));
        assert_eq!(CompareHillFormulas(&mut heap, c2h6, c3h6), Ok(1));
        assert_eq!(CompareHillFormulas(&mut heap, c3h6, c2h6), Ok(-1));
        assert_eq!(CompareHillFormulas(&mut heap, c2n, c2o), Ok(-1));
        assert_eq!(CompareHillFormulas(&mut heap, c2o, c2n), Ok(1));
        assert_eq!(CompareHillFormulas(&mut heap, oxygen, empty), Ok(-11));
        assert_eq!(CompareHillFormulas(&mut heap, empty, oxygen), Ok(11));
        assert_eq!(CompareHillFormulas(&mut heap, invalid, oxygen), Ok(0));
        assert_eq!(
            CompareHillFormulas(&mut heap, SourceConstPointer::null(), empty),
            Ok(0)
        );
        assert_eq!(
            CompareHillFormulas(&mut heap, max_count, wrapped_count),
            Ok(1)
        );
    }

    #[test]
    fn source_port__ichimake__compinchi2__line_1712() {
        fn inchi(heap: &mut SourceHeap, formula: &[u8], atoms: &[u8], hydrogens: &[i8]) -> INChI {
            let mut formula_bytes: Vec<i8> = formula.iter().map(|byte| *byte as i8).collect();
            formula_bytes.push(0);
            INChI {
                nNumberOfAtoms: atoms.len() as i32,
                szHillFormula: heap.allocate(formula_bytes).unwrap(),
                nAtom: heap.allocate(atoms.to_vec()).unwrap(),
                nNum_H: heap.allocate(hydrogens.to_vec()).unwrap(),
                ..INChI::default()
            }
        }

        fn compare(
            heap: &mut SourceHeap,
            left: INChI,
            right: INChI,
            taut: u32,
            isotopic: i32,
        ) -> i32 {
            let mut first = INCHI_SORT::default();
            let mut second = INCHI_SORT::default();
            first.pINChI[TAUT_NON as usize] = heap.allocate(vec![left]).unwrap();
            second.pINChI[TAUT_NON as usize] = heap.allocate(vec![right]).unwrap();
            CompINChI2(heap, &first, &second, taut, isotopic).unwrap()
        }

        let mut heap = SourceHeap::default();
        let empty1 = INCHI_SORT::default();
        let empty2 = INCHI_SORT::default();
        assert_eq!(CompINChI2(&mut heap, &empty1, &empty2, TAUT_NON, 1), Ok(0));

        let live = inchi(&mut heap, b"C", &[6], &[0]);
        let mut only_left = INCHI_SORT::default();
        only_left.pINChI[TAUT_NON as usize] = heap.allocate(vec![live.clone()]).unwrap();
        assert_eq!(
            CompINChI2(&mut heap, &only_left, &empty2, TAUT_NON, 1),
            Ok(-1)
        );
        assert_eq!(
            CompINChI2(&mut heap, &empty1, &only_left, TAUT_NON, 1),
            Ok(1)
        );

        let mut deleted = live.clone();
        deleted.bDeleted = 1;
        assert_eq!(
            compare(&mut heap, deleted.clone(), live.clone(), TAUT_NON, 1),
            1
        );
        assert_eq!(compare(&mut heap, live.clone(), deleted, TAUT_NON, 1), -1);

        let nitrogen = inchi(&mut heap, b"N", &[7], &[0]);
        assert!(compare(&mut heap, live.clone(), nitrogen, TAUT_NON, 1) < 0);
        let c2 = inchi(&mut heap, b"C2", &[6, 6], &[0, 0]);
        assert_eq!(compare(&mut heap, live.clone(), c2, TAUT_NON, 1), 1);

        let mut atom2 = live.clone();
        atom2.nAtom = heap.allocate(vec![7]).unwrap();
        assert_eq!(compare(&mut heap, live.clone(), atom2, TAUT_NON, 1), 1);
        let mut conn1 = live.clone();
        conn1.lenConnTable = 1;
        conn1.nConnTable = heap.allocate(vec![1_u16]).unwrap();
        let mut conn2 = live.clone();
        conn2.lenConnTable = 1;
        conn2.nConnTable = heap.allocate(vec![2_u16]).unwrap();
        assert_eq!(
            compare(&mut heap, conn1.clone(), conn2.clone(), TAUT_NON, 1),
            1
        );
        conn2.lenConnTable = 2;
        conn2.nConnTable = heap.allocate(vec![1_u16, 2]).unwrap();
        assert_eq!(compare(&mut heap, conn1, conn2, TAUT_NON, 1), 1);

        let ch2 = inchi(&mut heap, b"CH2", &[6], &[0]);
        let ch3 = inchi(&mut heap, b"CH3", &[6], &[0]);
        assert_eq!(compare(&mut heap, ch2, ch3, TAUT_NON, 1), 1);
        let h1 = inchi(&mut heap, b"C", &[6], &[1]);
        let h2 = inchi(&mut heap, b"C", &[6], &[2]);
        assert_eq!(compare(&mut heap, h1, h2, TAUT_NON, 1), 1);

        let mut taut1 = live.clone();
        taut1.lenTautomer = 1;
        taut1.nTautomer = heap.allocate(vec![1_u16]).unwrap();
        let mut taut2 = live.clone();
        taut2.lenTautomer = 1;
        taut2.nTautomer = heap.allocate(vec![2_u16]).unwrap();
        assert_eq!(compare(&mut heap, taut1, taut2, TAUT_YES, 1), 1);

        let mut fixed1 = live.clone();
        fixed1.nNum_H_fixed = heap.allocate(vec![1_i8]).unwrap();
        let mut fixed2 = live.clone();
        fixed2.nNum_H_fixed = heap.allocate(vec![2_i8]).unwrap();
        let mut fixed_sort1 = INCHI_SORT::default();
        fixed_sort1.pINChI[TAUT_YES as usize] = heap.allocate(vec![live.clone()]).unwrap();
        fixed_sort1.pINChI[TAUT_NON as usize] = heap.allocate(vec![fixed1]).unwrap();
        let mut fixed_sort2 = INCHI_SORT::default();
        fixed_sort2.pINChI[TAUT_YES as usize] = heap.allocate(vec![live.clone()]).unwrap();
        fixed_sort2.pINChI[TAUT_NON as usize] = heap.allocate(vec![fixed2]).unwrap();
        assert_eq!(
            CompINChI2(&mut heap, &fixed_sort1, &fixed_sort2, TAUT_NON, 1),
            Ok(1)
        );

        let mut stereo_inchi = live.clone();
        let stereo = INChI_Stereo {
            nNumberOfStereoBonds: 1,
            nBondAtom1: heap.allocate(vec![1_u16]).unwrap(),
            nBondAtom2: heap.allocate(vec![2_u16]).unwrap(),
            b_parity: heap.allocate(vec![1_i8]).unwrap(),
            ..INChI_Stereo::default()
        };
        stereo_inchi.Stereo = heap.allocate(vec![stereo]).unwrap();
        assert_eq!(
            compare(&mut heap, live.clone(), stereo_inchi, TAUT_NON, 1),
            1
        );

        let mut isotope1 = live.clone();
        isotope1.nNumberOfIsotopicAtoms = 1;
        isotope1.IsotopicAtom = heap
            .allocate(vec![crate::source_types::INChI_IsotopicAtom {
                nAtomNumber: 1,
                nIsoDifference: 1,
                nNum_H: 1,
                nNum_D: 1,
                nNum_T: 1,
            }])
            .unwrap();
        let mut isotope2 = isotope1.clone();
        isotope2.IsotopicAtom = heap
            .allocate(vec![crate::source_types::INChI_IsotopicAtom {
                nAtomNumber: 1,
                nIsoDifference: 1,
                nNum_H: 1,
                nNum_D: 1,
                nNum_T: 2,
            }])
            .unwrap();
        assert_eq!(
            compare(&mut heap, isotope1.clone(), isotope2.clone(), TAUT_NON, 0),
            0
        );
        assert_eq!(compare(&mut heap, isotope1, isotope2, TAUT_NON, 1), 1);

        let mut group1 = live.clone();
        group1.nNumberOfIsotopicTGroups = 1;
        group1.IsotopicTGroup = heap
            .allocate(vec![crate::source_types::INChI_IsotopicTGroup {
                nTGroupNumber: 1,
                nNum_H: 1,
                nNum_D: 1,
                nNum_T: 1,
            }])
            .unwrap();
        let mut group2 = group1.clone();
        group2.IsotopicTGroup = heap
            .allocate(vec![crate::source_types::INChI_IsotopicTGroup {
                nTGroupNumber: 2,
                nNum_H: 1,
                nNum_D: 1,
                nNum_T: 1,
            }])
            .unwrap();
        assert_eq!(compare(&mut heap, group1, group2, TAUT_NON, 1), 1);

        let mut charged1 = live.clone();
        charged1.nTotalCharge = -1;
        let mut charged2 = live.clone();
        charged2.nTotalCharge = 1;
        assert_eq!(compare(&mut heap, charged1, charged2, TAUT_NON, 1), -2);
        let mut charged = live.clone();
        charged.nTotalCharge = 1;
        assert_eq!(compare(&mut heap, charged, live, TAUT_NON, 1), 1);
    }

    #[test]
    fn source_port__ichimake__compinchinontaut2__line_2042() {
        fn inchi(heap: &mut SourceHeap, charge: i32) -> INChI {
            INChI {
                nTotalCharge: charge,
                nNumberOfAtoms: 1,
                szHillFormula: heap.allocate(vec![b'C' as i8, 0]).unwrap(),
                nAtom: heap.allocate(vec![6_u8]).unwrap(),
                nNum_H: heap.allocate(vec![0_i8]).unwrap(),
                ..INChI::default()
            }
        }

        let mut heap = SourceHeap::default();
        let mut first = INCHI_SORT::default();
        let mut second = INCHI_SORT::default();
        let charged_left = inchi(&mut heap, -1);
        let charged_right = inchi(&mut heap, 1);
        first.pINChI[TAUT_NON as usize] = heap.allocate(vec![charged_left]).unwrap();
        second.pINChI[TAUT_NON as usize] = heap.allocate(vec![charged_right]).unwrap();
        assert_eq!(CompINChINonTaut2(&mut heap, &first, &second), Ok(-2));

        let mut stable1 = INCHI_SORT {
            ord_number: i16::MAX,
            ..INCHI_SORT::default()
        };
        let mut stable2 = INCHI_SORT {
            ord_number: i16::MIN,
            ..INCHI_SORT::default()
        };
        let equal1 = inchi(&mut heap, 0);
        let equal2 = inchi(&mut heap, 0);
        stable1.pINChI[TAUT_NON as usize] = heap.allocate(vec![equal1]).unwrap();
        stable2.pINChI[TAUT_NON as usize] = heap.allocate(vec![equal2]).unwrap();
        assert_eq!(CompINChINonTaut2(&mut heap, &stable1, &stable2), Ok(65535));

        let mut mobile_stereo = inchi(&mut heap, 0);
        let bond_atom1 = heap.allocate(vec![1_u16]).unwrap();
        let bond_atom2 = heap.allocate(vec![2_u16]).unwrap();
        let bond_parity = heap.allocate(vec![1_i8]).unwrap();
        mobile_stereo.Stereo = heap
            .allocate(vec![INChI_Stereo {
                nNumberOfStereoBonds: 1,
                nBondAtom1: bond_atom1,
                nBondAtom2: bond_atom2,
                b_parity: bond_parity,
                ..INChI_Stereo::default()
            }])
            .unwrap();
        let mut fallback1 = INCHI_SORT::default();
        let mut fallback2 = INCHI_SORT::default();
        let mobile_plain = inchi(&mut heap, 0);
        let fixed_plain1 = inchi(&mut heap, 0);
        let fixed_plain2 = inchi(&mut heap, 0);
        fallback1.pINChI[TAUT_YES as usize] = heap.allocate(vec![mobile_plain]).unwrap();
        fallback2.pINChI[TAUT_YES as usize] = heap.allocate(vec![mobile_stereo]).unwrap();
        fallback1.pINChI[TAUT_NON as usize] = heap.allocate(vec![fixed_plain1]).unwrap();
        fallback2.pINChI[TAUT_NON as usize] = heap.allocate(vec![fixed_plain2]).unwrap();
        assert_eq!(CompINChINonTaut2(&mut heap, &fallback1, &fallback2), Ok(1));
    }

    #[test]
    fn source_port__ichimake__compinchitaut2__line_2064() {
        fn inchi(heap: &mut SourceHeap, charge: i32) -> INChI {
            INChI {
                nTotalCharge: charge,
                nNumberOfAtoms: 1,
                szHillFormula: heap.allocate(vec![b'C' as i8, 0]).unwrap(),
                nAtom: heap.allocate(vec![6_u8]).unwrap(),
                nNum_H: heap.allocate(vec![0_i8]).unwrap(),
                ..INChI::default()
            }
        }

        let mut heap = SourceHeap::default();
        let mut first = INCHI_SORT::default();
        let mut second = INCHI_SORT::default();
        let charged_left = inchi(&mut heap, -1);
        let charged_right = inchi(&mut heap, 1);
        first.pINChI[TAUT_YES as usize] = heap.allocate(vec![charged_left]).unwrap();
        second.pINChI[TAUT_YES as usize] = heap.allocate(vec![charged_right]).unwrap();
        assert_eq!(CompINChITaut2(&mut heap, &first, &second), Ok(-2));

        let mut stable1 = INCHI_SORT {
            ord_number: i16::MAX,
            ..INCHI_SORT::default()
        };
        let mut stable2 = INCHI_SORT {
            ord_number: i16::MIN,
            ..INCHI_SORT::default()
        };
        let equal1 = inchi(&mut heap, 0);
        let equal2 = inchi(&mut heap, 0);
        stable1.pINChI[TAUT_YES as usize] = heap.allocate(vec![equal1]).unwrap();
        stable2.pINChI[TAUT_YES as usize] = heap.allocate(vec![equal2]).unwrap();
        assert_eq!(CompINChITaut2(&mut heap, &stable1, &stable2), Ok(65535));

        let mut fixed_stereo = inchi(&mut heap, 0);
        let bond_atom1 = heap.allocate(vec![1_u16]).unwrap();
        let bond_atom2 = heap.allocate(vec![2_u16]).unwrap();
        let bond_parity = heap.allocate(vec![1_i8]).unwrap();
        fixed_stereo.Stereo = heap
            .allocate(vec![INChI_Stereo {
                nNumberOfStereoBonds: 1,
                nBondAtom1: bond_atom1,
                nBondAtom2: bond_atom2,
                b_parity: bond_parity,
                ..INChI_Stereo::default()
            }])
            .unwrap();
        let mut fallback1 = INCHI_SORT::default();
        let mut fallback2 = INCHI_SORT::default();
        let mobile_plain1 = inchi(&mut heap, 0);
        let mobile_plain2 = inchi(&mut heap, 0);
        let fixed_plain = inchi(&mut heap, 0);
        fallback1.pINChI[TAUT_YES as usize] = heap.allocate(vec![mobile_plain1]).unwrap();
        fallback2.pINChI[TAUT_YES as usize] = heap.allocate(vec![mobile_plain2]).unwrap();
        fallback1.pINChI[TAUT_NON as usize] = heap.allocate(vec![fixed_plain]).unwrap();
        fallback2.pINChI[TAUT_NON as usize] = heap.allocate(vec![fixed_stereo]).unwrap();
        assert_eq!(CompINChITaut2(&mut heap, &fallback1, &fallback2), Ok(1));
    }

    #[test]
    fn source_port__ichimake__comparehillformulasnoh__line_272() {
        let mut heap = SourceHeap::default();
        let mut formula = |bytes: &[u8]| {
            let mut bytes: Vec<i8> = bytes.iter().map(|byte| *byte as i8).collect();
            bytes.push(0);
            heap.allocate(bytes).unwrap().as_const()
        };
        let c2h6 = formula(b"C2H6");
        let c2h4 = formula(b"C2H4");
        let c2n = formula(b"C2N");
        let c3n = formula(b"C3N");
        let c2o = formula(b"C2O");
        let h2o = formula(b"H2O");
        let oxygen = formula(b"O");
        let hydrogen = formula(b"H2");
        let empty = formula(b"");
        let invalid = formula(b"c2");
        let mut h1 = 10;
        let mut h2 = 20;
        assert_eq!(
            CompareHillFormulasNoH(&mut heap, c2h6, c2h4, &mut h1, &mut h2),
            Ok(0)
        );
        assert_eq!((h1, h2), (16, 24));

        h1 = 0;
        h2 = 0;
        assert_eq!(
            CompareHillFormulasNoH(&mut heap, c2n, c3n, &mut h1, &mut h2),
            Ok(1)
        );
        assert_eq!(
            CompareHillFormulasNoH(&mut heap, c2n, c2o, &mut h1, &mut h2),
            Ok(-1)
        );
        assert_eq!(
            CompareHillFormulasNoH(&mut heap, h2o, oxygen, &mut h1, &mut h2),
            Ok(0)
        );
        assert_eq!((h1, h2), (2, 0));
        assert_eq!(
            CompareHillFormulasNoH(&mut heap, hydrogen, empty, &mut h1, &mut h2),
            Ok(0)
        );
        assert_eq!(h1, 4);
        assert_eq!(
            CompareHillFormulasNoH(&mut heap, empty, empty, &mut h1, &mut h2),
            Ok(0)
        );
        assert_eq!(
            CompareHillFormulasNoH(&mut heap, invalid, oxygen, &mut h1, &mut h2),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichimake__comparetautnonisopartofinchi__line_316() {
        let mut heap = SourceHeap::default();
        let empty1 = INChI::default();
        let empty2 = INChI {
            lenTautomer: -1,
            ..INChI::default()
        };
        assert_eq!(CompareTautNonIsoPartOfINChI(&heap, &empty1, &empty2), Ok(0));

        let zero_first1 = heap.allocate(vec![0_u16, 9]).unwrap();
        let zero_first2 = heap.allocate(vec![0_u16, 7]).unwrap();
        let ignored1 = INChI {
            lenTautomer: 2,
            nTautomer: zero_first1,
            ..INChI::default()
        };
        let ignored2 = INChI {
            lenTautomer: 2,
            nTautomer: zero_first2,
            ..INChI::default()
        };
        assert_eq!(
            CompareTautNonIsoPartOfINChI(&heap, &ignored1, &ignored2),
            Ok(0)
        );

        let values1 = heap.allocate(vec![2_u16, 3]).unwrap();
        let values2 = heap.allocate(vec![2_u16, 3]).unwrap();
        let mut left = INChI {
            lenTautomer: 2,
            nTautomer: values1,
            ..INChI::default()
        };
        let mut right = INChI {
            lenTautomer: 2,
            nTautomer: values2,
            ..INChI::default()
        };
        assert_eq!(CompareTautNonIsoPartOfINChI(&heap, &left, &right), Ok(0));
        right.lenTautomer = 3;
        assert_eq!(CompareTautNonIsoPartOfINChI(&heap, &left, &right), Ok(1));
        right.lenTautomer = 2;
        right.nTautomer = heap.allocate(vec![2_u16, 5]).unwrap();
        assert_eq!(CompareTautNonIsoPartOfINChI(&heap, &left, &right), Ok(2));
        left.nTautomer = heap.allocate(vec![4_u16, 3]).unwrap();
        right.nTautomer = heap.allocate(vec![2_u16, 3]).unwrap();
        assert_eq!(CompareTautNonIsoPartOfINChI(&heap, &left, &right), Ok(-2));
    }

    fn comp_taut_fixture(heap: &mut SourceHeap) -> (INChI, INChI) {
        let make = |heap: &mut SourceHeap| INChI {
            nNumberOfAtoms: 2,
            szHillFormula: heap.allocate(vec![b'C' as i8, b'2' as i8, 0]).unwrap(),
            nAtom: heap.allocate(vec![6_u8, 6]).unwrap(),
            lenConnTable: 2,
            nConnTable: heap.allocate(vec![1_u16, 2]).unwrap(),
            lenTautomer: 2,
            nTautomer: heap.allocate(vec![2_u16, 1]).unwrap(),
            nNum_H: heap.allocate(vec![1_i8, 1]).unwrap(),
            ..INChI::default()
        };
        (make(heap), make(heap))
    }

    fn comp_taut_pair(
        heap: &mut SourceHeap,
        i1: INChI,
        i2: INChI,
        compare_isotopic: i32,
    ) -> Result<i32, SourceHeapError> {
        let i1 = heap.allocate(vec![i1]).unwrap();
        let i2 = heap.allocate(vec![i2]).unwrap();
        let mut p1 = INCHI_SORT::default();
        let mut p2 = INCHI_SORT::default();
        p1.pINChI[TAUT_YES as usize] = i1;
        p2.pINChI[TAUT_NON as usize] = i2;
        CompINChITautVsNonTaut(heap, &p1, &p2, compare_isotopic)
    }

    fn assert_comp_taut<F>(expected: i32, compare_isotopic: i32, update: F)
    where
        F: FnOnce(&mut SourceHeap, &mut INChI, &mut INChI),
    {
        let mut heap = SourceHeap::default();
        let (mut i1, mut i2) = comp_taut_fixture(&mut heap);
        update(&mut heap, &mut i1, &mut i2);
        assert_eq!(
            comp_taut_pair(&mut heap, i1, i2, compare_isotopic),
            Ok(expected)
        );
    }

    #[test]
    fn source_port__ichimake__compinchitautvsnontaut__line_341() {
        let mut heap = SourceHeap::default();
        let (i1, i2) = comp_taut_fixture(&mut heap);
        let i1_pointer = heap.allocate(vec![i1.clone()]).unwrap();
        let i2_pointer = heap.allocate(vec![i2.clone()]).unwrap();
        let mut p1 = INCHI_SORT::default();
        let mut p2 = INCHI_SORT::default();
        p1.pINChI[TAUT_NON as usize] = i1_pointer;
        p2.pINChI[TAUT_NON as usize] = i2_pointer;
        assert_eq!(CompINChITautVsNonTaut(&mut heap, &p1, &p2, 1), Ok(0));

        p1 = INCHI_SORT::default();
        p1.pINChI[TAUT_YES as usize] = i1_pointer;
        p2 = INCHI_SORT::default();
        assert_eq!(CompINChITautVsNonTaut(&mut heap, &p1, &p2, 1), Ok(0));
        assert_eq!(
            CompINChITautVsNonTaut(&mut heap, &INCHI_SORT::default(), &INCHI_SORT::default(), 1,),
            Ok(0)
        );

        assert_comp_taut(1, 0, |_, i1, _| i1.bDeleted = 1);
        assert_comp_taut(-1, 0, |_, _, i2| i2.bDeleted = 1);
        assert_comp_taut(1, 0, |heap, i1, i2| {
            i1.szHillFormula = heap
                .allocate(vec![b'C' as i8, b'2' as i8, b'N' as i8, 0])
                .unwrap();
            i2.szHillFormula = heap
                .allocate(vec![b'C' as i8, b'3' as i8, b'N' as i8, 0])
                .unwrap();
        });
        assert_comp_taut(1, 0, |heap, i1, i2| {
            i1.szHillFormula = heap
                .allocate(vec![b'C' as i8, b'2' as i8, b'H' as i8, b'2' as i8, 0])
                .unwrap();
            i2.szHillFormula = heap
                .allocate(vec![b'C' as i8, b'2' as i8, b'H' as i8, b'3' as i8, 0])
                .unwrap();
        });
        assert_comp_taut(1, 0, |_, _, i2| i2.nNumberOfAtoms = 3);
        assert_comp_taut(1, 0, |heap, _, i2| {
            i2.nAtom = heap.allocate(vec![6_u8, 7]).unwrap();
        });
        assert_comp_taut(1, 0, |_, _, i2| i2.lenConnTable = 3);
        assert_comp_taut(2, 0, |heap, _, i2| {
            i2.nConnTable = heap.allocate(vec![1_u16, 4]).unwrap();
        });
        assert_comp_taut(0, 0, |_, i1, i2| {
            i1.lenConnTable = 0;
            i1.nConnTable = Default::default();
            i2.lenConnTable = 0;
            i2.nConnTable = Default::default();
        });

        assert_comp_taut(1, 0, |heap, _, i2| {
            i2.nNum_H = heap.allocate(vec![0_i8, 1]).unwrap();
        });
        assert_comp_taut(-1, 0, |heap, i1, _| {
            i1.nNum_H = heap.allocate(vec![0_i8, 1]).unwrap();
        });
        assert_comp_taut(2, 0, |heap, _, i2| {
            i2.nNum_H = heap.allocate(vec![3_i8, 1]).unwrap();
        });
        assert_comp_taut(1, 0, |_, _, i2| i2.lenTautomer = 3);
        assert_comp_taut(2, 0, |heap, _, i2| {
            i2.nTautomer = heap.allocate(vec![2_u16, 3]).unwrap();
        });
        assert_comp_taut(0, 0, |heap, _, i2| {
            i2.nNum_H_fixed = heap.allocate(vec![0_i8, 0]).unwrap();
        });
        assert_comp_taut(1, 0, |heap, _, i2| {
            i2.nNum_H_fixed = heap.allocate(vec![0_i8, 1]).unwrap();
        });
        assert_comp_taut(1, 0, |heap, _, i2| {
            let stereo = stereo_fixture(heap, &[1], &[1], 0, &[], &[], &[]);
            i2.Stereo = heap.allocate(vec![stereo]).unwrap();
        });

        assert_comp_taut(0, 0, |_, _, i2| i2.nNumberOfIsotopicAtoms = 1);
        assert_comp_taut(1, 1, |_, _, i2| i2.nNumberOfIsotopicAtoms = 1);
        for (field, expected) in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)] {
            assert_comp_taut(expected, 1, |heap, i1, i2| {
                let left = INChI_IsotopicAtom {
                    nAtomNumber: 1,
                    ..INChI_IsotopicAtom::default()
                };
                let mut right = left.clone();
                match field {
                    0 => right.nAtomNumber += 1,
                    1 => right.nIsoDifference += 2,
                    2 => right.nNum_T += 3,
                    3 => right.nNum_D += 4,
                    _ => right.nNum_H += 5,
                }
                i1.nNumberOfIsotopicAtoms = 1;
                i2.nNumberOfIsotopicAtoms = 1;
                i1.IsotopicAtom = heap.allocate(vec![left]).unwrap();
                i2.IsotopicAtom = heap.allocate(vec![right]).unwrap();
            });
        }
        assert_comp_taut(1, 1, |_, i1, _| i1.nNumberOfIsotopicTGroups = -7);
        assert_comp_taut(1, 1, |heap, _, i2| {
            let stereo = stereo_fixture(heap, &[1], &[1], 0, &[], &[], &[]);
            i2.StereoIsotopic = heap.allocate(vec![stereo]).unwrap();
        });

        assert_comp_taut(-3, 1, |_, i1, i2| {
            i1.nTotalCharge = -2;
            i2.nTotalCharge = 1;
        });
        assert_comp_taut(1, 1, |_, i1, _| i1.nTotalCharge = 1);
        assert_comp_taut(-1, 1, |_, _, i2| i2.nTotalCharge = 1);
        assert_comp_taut(0, 1, |_, _, _| {});
    }

    #[test]
    fn source_port__ichimake__getsp3relracabs__line_593() {
        const SP3_NONE: i32 = 0;
        const SP3_ONLY: i32 = 1;
        const SP3_ABS: i32 = 2;
        const SP3_REL: i32 = 4;

        let plain = INChI::default();
        let mut stereo = INChI_Stereo::default();
        assert_eq!(GetSp3RelRacAbs(None, None), SP3_NONE);
        assert_eq!(GetSp3RelRacAbs(Some(&plain), None), SP3_NONE);
        assert_eq!(GetSp3RelRacAbs(None, Some(&stereo)), SP3_NONE);
        assert_eq!(GetSp3RelRacAbs(Some(&plain), Some(&stereo)), SP3_NONE);

        stereo.nNumberOfStereoCenters = 1;
        stereo.nCompInv2Abs = 1;
        let deleted = INChI {
            bDeleted: 1,
            ..INChI::default()
        };
        assert_eq!(GetSp3RelRacAbs(Some(&deleted), Some(&stereo)), SP3_NONE);
        assert_eq!(GetSp3RelRacAbs(Some(&plain), Some(&stereo)), SP3_ABS);

        let relative = INChI {
            nFlags: INCHI_MODE::from(INCHI_FLAG_REL_STEREO),
            ..INChI::default()
        };
        let racemic = INChI {
            nFlags: INCHI_MODE::from(INCHI_FLAG_RAC_STEREO),
            ..INChI::default()
        };
        assert_eq!(GetSp3RelRacAbs(Some(&relative), Some(&stereo)), SP3_NONE);
        assert_eq!(GetSp3RelRacAbs(Some(&racemic), Some(&stereo)), SP3_NONE);
        stereo.nNumberOfStereoCenters = 2;
        assert_eq!(GetSp3RelRacAbs(Some(&relative), Some(&stereo)), SP3_REL);
        assert_eq!(GetSp3RelRacAbs(Some(&racemic), Some(&stereo)), SP3_REL);

        stereo.nCompInv2Abs = 0;
        stereo.nNumberOfStereoCenters = 1;
        assert_eq!(GetSp3RelRacAbs(Some(&plain), Some(&stereo)), SP3_ONLY);
        assert_eq!(GetSp3RelRacAbs(Some(&relative), Some(&stereo)), SP3_NONE);
        assert_eq!(GetSp3RelRacAbs(Some(&racemic), Some(&stereo)), SP3_NONE);
        stereo.nNumberOfStereoCenters = 2;
        assert_eq!(GetSp3RelRacAbs(Some(&relative), Some(&stereo)), SP3_ONLY);
        assert_eq!(GetSp3RelRacAbs(Some(&racemic), Some(&stereo)), SP3_ONLY);
    }

    fn layers_inchi(heap: &mut SourceHeap) -> INChI {
        let stereo = stereo_fixture(heap, &[1, 2], &[1, 2], 1, &[1], &[2], &[1]);
        let isotopic_stereo = stereo_fixture(heap, &[1, 2], &[1, 2], 1, &[1], &[2], &[1]);
        INChI {
            nTotalCharge: 1,
            nNumberOfAtoms: 2,
            szHillFormula: heap.allocate(vec![b'C' as i8, b'2' as i8, 0]).unwrap(),
            nAtom: heap.allocate(vec![6_u8, 6]).unwrap(),
            lenConnTable: 2,
            nConnTable: heap.allocate(vec![1_u16, 2]).unwrap(),
            lenTautomer: 1,
            nTautomer: heap.allocate(vec![0_u16]).unwrap(),
            nNum_H: heap.allocate(vec![1_i8, 0]).unwrap(),
            nNum_H_fixed: heap.allocate(vec![1_i8, 0]).unwrap(),
            nNumberOfIsotopicAtoms: 1,
            IsotopicAtom: heap
                .allocate(vec![INChI_IsotopicAtom {
                    nAtomNumber: 1,
                    nIsoDifference: 1,
                    nNum_H: 1,
                    nNum_D: 1,
                    nNum_T: 1,
                }])
                .unwrap(),
            nNumberOfIsotopicTGroups: 1,
            IsotopicTGroup: heap
                .allocate(vec![INChI_IsotopicTGroup {
                    nTGroupNumber: 1,
                    nNum_H: 1,
                    nNum_D: 1,
                    nNum_T: 1,
                }])
                .unwrap(),
            Stereo: heap.allocate(vec![stereo]).unwrap(),
            StereoIsotopic: heap.allocate(vec![isotopic_stereo]).unwrap(),
            ..INChI::default()
        }
    }

    fn run_layers(
        heap: &mut SourceHeap,
        first: Option<INChI>,
        second: Option<INChI>,
        fix_charge: i32,
    ) -> (i32, [[i8; 11]; 4]) {
        let mut p1 = INCHI_SORT::default();
        let mut p2 = INCHI_SORT::default();
        if let Some(first) = first {
            p1.pINChI[TAUT_YES as usize] = heap.allocate(vec![first]).unwrap();
        }
        if let Some(second) = second {
            p2.pINChI[TAUT_NON as usize] = heap.allocate(vec![second]).unwrap();
        }
        let mut segments = [[0_i8; 11]; 4];
        let result = CompINChILayers(heap, &p1, &p2, &mut segments, fix_charge).unwrap();
        (result, segments)
    }

    #[test]
    fn source_port__ichimake__compinchilayers__line_645() {
        const EQL: i8 = tagMarkDiff_DIFV_EQL2PRECED as i8;
        const NEQ: i8 = tagMarkDiff_DIFV_NEQ2PRECED as i8;
        const EMPTY: i8 = tagMarkDiff_DIFV_IS_EMPTY as i8;
        const FI_MI: i8 = tagMarkDiff_DIFV_FI_EQ_MI as i8;

        let mut heap = SourceHeap::default();
        assert_eq!(run_layers(&mut heap, None, None, 0), (0, [[0; 11]; 4]));

        let first = layers_inchi(&mut heap);
        let second = layers_inchi(&mut heap);
        let (result, segments) = run_layers(&mut heap, Some(first), Some(second), 0);
        assert_eq!(result, 0);
        let mut expected = [[0_i8; 11]; 4];
        expected[tagDiffINChILayers_DIFL_M as usize]
            [tagDiffINChISegments_DIFS_f_FORMULA as usize] = NEQ;
        expected[tagDiffINChILayers_DIFL_F as usize]
            [tagDiffINChISegments_DIFS_f_FORMULA as usize] = EQL;
        expected[tagDiffINChILayers_DIFL_M as usize]
            [tagDiffINChISegments_DIFS_h_H_ATOMS as usize] = NEQ;
        expected[tagDiffINChILayers_DIFL_F as usize]
            [tagDiffINChISegments_DIFS_h_H_ATOMS as usize] = NEQ;
        expected[tagDiffINChILayers_DIFL_M as usize][tagDiffINChISegments_DIFS_q_CHARGE as usize] =
            NEQ;
        expected[tagDiffINChILayers_DIFL_F as usize][tagDiffINChISegments_DIFS_q_CHARGE as usize] =
            EQL;
        expected[tagDiffINChILayers_DIFL_M as usize][tagDiffINChISegments_DIFS_b_SBONDS as usize] =
            NEQ;
        expected[tagDiffINChILayers_DIFL_F as usize][tagDiffINChISegments_DIFS_b_SBONDS as usize] =
            EQL;
        expected[tagDiffINChILayers_DIFL_MI as usize]
            [tagDiffINChISegments_DIFS_b_SBONDS as usize] = EQL;
        expected[tagDiffINChILayers_DIFL_FI as usize]
            [tagDiffINChISegments_DIFS_b_SBONDS as usize] = EQL;
        for segment in [
            tagDiffINChISegments_DIFS_t_SATOMS,
            tagDiffINChISegments_DIFS_m_SP3INV,
            tagDiffINChISegments_DIFS_s_STYPE,
        ] {
            expected[tagDiffINChILayers_DIFL_M as usize][segment as usize] = NEQ;
            expected[tagDiffINChILayers_DIFL_F as usize][segment as usize] = EQL;
            expected[tagDiffINChILayers_DIFL_MI as usize][segment as usize] = EQL;
            expected[tagDiffINChILayers_DIFL_FI as usize][segment as usize] = EQL;
        }
        expected[tagDiffINChILayers_DIFL_MI as usize]
            [tagDiffINChISegments_DIFS_i_IATOMS as usize] = NEQ;
        expected[tagDiffINChILayers_DIFL_FI as usize]
            [tagDiffINChISegments_DIFS_i_IATOMS as usize] = FI_MI;
        assert_eq!(segments, expected);

        let first = layers_inchi(&mut heap);
        let mut second = layers_inchi(&mut heap);
        second.szHillFormula = heap.allocate(vec![b'C' as i8, b'3' as i8, 0]).unwrap();
        second.nTotalCharge = 0;
        second.nNum_H_fixed = heap.allocate(vec![0_i8, 0]).unwrap();
        second.Stereo = Default::default();
        second.StereoIsotopic = Default::default();
        second.nNumberOfIsotopicAtoms = 0;
        second.nNumberOfIsotopicTGroups = 0;
        let (_, changed) = run_layers(&mut heap, Some(first), Some(second), 0);
        assert_eq!(changed[tagDiffINChILayers_DIFL_F as usize][0], NEQ);
        assert_eq!(changed[tagDiffINChILayers_DIFL_F as usize][2], 0);
        assert_eq!(changed[tagDiffINChILayers_DIFL_F as usize][3], EMPTY);
        assert_eq!(changed[tagDiffINChILayers_DIFL_F as usize][5], EMPTY);
        assert_eq!(changed[tagDiffINChILayers_DIFL_FI as usize][9], EMPTY);

        let first = layers_inchi(&mut heap);
        let mut deleted = layers_inchi(&mut heap);
        deleted.bDeleted = 1;
        let (_, missing) = run_layers(&mut heap, Some(first), Some(deleted), 0);
        assert_eq!(missing[tagDiffINChILayers_DIFL_F as usize][0], EMPTY);
        assert_eq!(missing[tagDiffINChILayers_DIFL_F as usize][3], EMPTY);

        let first = layers_inchi(&mut heap);
        let mut second = layers_inchi(&mut heap);
        second.IsotopicTGroup = heap
            .allocate(vec![INChI_IsotopicTGroup {
                nTGroupNumber: 1,
                nNum_H: 1,
                nNum_D: 4,
                nNum_T: 1,
            }])
            .unwrap();
        let (difference, _) = run_layers(&mut heap, Some(first), Some(second), 0);
        assert_eq!(difference, 3);

        let mut first = layers_inchi(&mut heap);
        first.nTotalCharge = 0;
        let first_pointer = heap.allocate(vec![first]).unwrap();
        let mut p1 = INCHI_SORT::default();
        p1.pINChI[TAUT_YES as usize] = first_pointer;
        p1.ord_number = 1;
        let mut p2_taut = layers_inchi(&mut heap);
        p2_taut.nTotalCharge = -1;
        let mut p2 = INCHI_SORT::default();
        p2.pINChI[TAUT_YES as usize] = heap.allocate(vec![p2_taut]).unwrap();
        p2.ord_number = 2;
        let mut fixed = [[0_i8; 11]; 4];
        assert_eq!(CompINChILayers(&mut heap, &p1, &p2, &mut fixed, 1), Ok(0));
        assert_eq!(fixed[tagDiffINChILayers_DIFL_F as usize][3], NEQ);
        assert_eq!(fixed[tagDiffINChILayers_DIFL_F as usize][10], NEQ);
    }

    fn stereo_fixture(
        heap: &mut SourceHeap,
        centers: &[u16],
        center_parity: &[i8],
        inversion: i32,
        bond_atom_1: &[u16],
        bond_atom_2: &[u16],
        bond_parity: &[i8],
    ) -> INChI_Stereo {
        INChI_Stereo {
            nNumberOfStereoCenters: centers.len() as i32,
            nNumber: heap.allocate(centers.to_vec()).unwrap(),
            t_parity: heap.allocate(center_parity.to_vec()).unwrap(),
            nCompInv2Abs: inversion,
            nNumberOfStereoBonds: bond_atom_1.len() as i32,
            nBondAtom1: heap.allocate(bond_atom_1.to_vec()).unwrap(),
            nBondAtom2: heap.allocate(bond_atom_2.to_vec()).unwrap(),
            b_parity: heap.allocate(bond_parity.to_vec()).unwrap(),
            ..INChI_Stereo::default()
        }
    }

    #[test]
    fn source_port__ichimake__comparereversedstereoinchi__line_2555() {
        let mut heap = SourceHeap::default();
        assert_eq!(CompareReversedStereoINChI(&heap, None, None), Ok(0));
        let empty = INChI_Stereo::default();
        assert_eq!(CompareReversedStereoINChI(&heap, Some(&empty), None), Ok(0));
        let center_only = INChI_Stereo {
            nNumberOfStereoCenters: 1,
            ..INChI_Stereo::default()
        };
        assert_eq!(
            CompareReversedStereoINChI(&heap, None, Some(&center_only)),
            Ok(20)
        );
        let bond_only = INChI_Stereo {
            nNumberOfStereoBonds: 1,
            ..INChI_Stereo::default()
        };
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&bond_only), None),
            Ok(20)
        );

        let base = stereo_fixture(&mut heap, &[1, 3], &[1, 2], 1, &[2], &[4], &[2]);
        let equal = stereo_fixture(&mut heap, &[1, 3], &[1, 2], 1, &[2], &[4], &[2]);
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&equal)),
            Ok(0)
        );

        let mut changed = equal.clone();
        changed.nNumberOfStereoCenters = 1;
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&changed)),
            Ok(21)
        );
        changed = stereo_fixture(&mut heap, &[1, 4], &[1, 2], 1, &[2], &[4], &[2]);
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&changed)),
            Ok(22)
        );
        changed = stereo_fixture(&mut heap, &[1, 3], &[1, 1], 1, &[2], &[4], &[2]);
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&changed)),
            Ok(23)
        );
        changed = stereo_fixture(&mut heap, &[1, 3], &[1, 2], 2, &[2], &[4], &[2]);
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&changed)),
            Ok(24)
        );
        changed.nCompInv2Abs = 0;
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&changed)),
            Ok(0)
        );

        changed = equal.clone();
        changed.nNumberOfStereoBonds = 0;
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&changed)),
            Ok(25)
        );
        changed = stereo_fixture(&mut heap, &[1, 3], &[1, 2], 1, &[3], &[4], &[2]);
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&changed)),
            Ok(26)
        );
        changed = stereo_fixture(&mut heap, &[1, 3], &[1, 2], 1, &[2], &[5], &[2]);
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&changed)),
            Ok(27)
        );
        changed = stereo_fixture(&mut heap, &[1, 3], &[1, 2], 1, &[2], &[4], &[1]);
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&changed)),
            Ok(28)
        );
    }

    #[test]
    fn source_port__ichimake__comparereversedinchi__line_2936() {
        let mut heap = SourceHeap::default();
        assert_eq!(CompareReversedINChI(&heap, None, None, None, None), Ok(0));
        let base = inchi_fixture(&mut heap);
        let equal = inchi_fixture(&mut heap);
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), None, None, None),
            Ok(1)
        );
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&equal), None, None),
            Ok(0)
        );

        let mut left = base.clone();
        let mut right = equal.clone();
        left.nErrorCode = 7;
        right.nErrorCode = 7;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(0)
        );
        right.nErrorCode = 8;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(2)
        );

        right = equal.clone();
        right.bDeleted = 1;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(1)
        );
        right = equal.clone();
        right.nNumberOfAtoms = 1;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(3)
        );
        right = equal.clone();
        right.nAtom = heap.allocate(vec![6_u8, 8]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(4)
        );
        right = equal.clone();
        right.szHillFormula = heap.allocate(vec![b'C' as i8, b'H' as i8, 0]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(7)
        );
        right = equal.clone();
        right.nNum_H = heap.allocate(vec![3_i8, 2]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(5)
        );
        left = base.clone();
        left.lenConnTable = 1;
        right = equal.clone();
        right.lenConnTable = 1;
        right.nNum_H = heap.allocate(vec![2_i8, 3]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(6)
        );

        left = base.clone();
        right = equal.clone();
        left.nNum_H_fixed = heap.allocate(vec![1_i8, 0]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(18)
        );
        left = base.clone();
        right.nNum_H_fixed = heap.allocate(vec![0_i8, -1]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(19)
        );
        left.nNum_H_fixed = heap.allocate(vec![2_i8, 0]).unwrap();
        right.nNum_H_fixed = heap.allocate(vec![1_i8, 1]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(20)
        );
        left.nNum_H_fixed = heap.allocate(vec![2_i8, 1]).unwrap();
        right.nNum_H_fixed = heap.allocate(vec![1_i8, 1]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(18)
        );
        left.nNum_H_fixed = heap.allocate(vec![1_i8, 1]).unwrap();
        right.nNum_H_fixed = heap.allocate(vec![2_i8, 1]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(19)
        );

        right = equal.clone();
        right.lenConnTable = 1;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(8)
        );
        right = equal.clone();
        right.nConnTable = heap.allocate(vec![1_u16, 3]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(9)
        );
        right = equal.clone();
        right.lenTautomer = 3;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(10)
        );
        right = equal.clone();
        right.nTautomer = heap.allocate(vec![2_u16, 2]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(11)
        );
        right = equal.clone();
        right.nNumberOfIsotopicAtoms = 0;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(12)
        );
        right = equal.clone();
        right.IsotopicAtom = heap
            .allocate(vec![INChI_IsotopicAtom {
                nAtomNumber: 2,
                ..INChI_IsotopicAtom::default()
            }])
            .unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(13)
        );
        right = equal.clone();
        right.nTotalCharge = -1;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(14)
        );

        let aux1 = INChI_Aux::default();
        let mut aux2 = INChI_Aux {
            nNumRemovedProtons: 1,
            ..INChI_Aux::default()
        };
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&equal), Some(&aux1), Some(&aux2)),
            Ok(16)
        );
        aux2.nNumRemovedProtons = 0;
        aux2.nNumRemovedIsotopicH = [0, 1, 0];
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&equal), Some(&aux1), Some(&aux2)),
            Ok(17)
        );
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&equal), Some(&aux1), None),
            Ok(0)
        );

        let stereo = stereo_fixture(&mut heap, &[1], &[1], 1, &[], &[], &[]);
        let stereo_pointer = heap.allocate(vec![stereo]).unwrap();
        left = base.clone();
        left.Stereo = stereo_pointer;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&equal), None, None),
            Ok(40)
        );
        left = base.clone();
        left.StereoIsotopic = stereo_pointer;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&equal), None, None),
            Ok(60)
        );
        left.Stereo = stereo_pointer;
        right = equal.clone();
        right.Stereo = stereo_pointer;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(0)
        );

        left = base.clone();
        right = equal.clone();
        left.lenTautomer = 1;
        right.lenTautomer = 0;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichimake__getinpstructerrortype__line_2480() {
        for error in [i32::MIN, -1, 1, 8, 10, 29] {
            assert_eq!(
                GetInpStructErrorType(None, error, None, i32::MAX).unwrap(),
                _IS_FATAL as i32
            );
        }
        assert_eq!(
            GetInpStructErrorType(None, 9, None, 0).unwrap(),
            _IS_ERROR as i32
        );
        for error in [30, 31, 97, 98, i32::MAX] {
            assert_eq!(
                GetInpStructErrorType(None, error, None, 1).unwrap(),
                _IS_ERROR as i32
            );
        }
        assert_eq!(
            GetInpStructErrorType(None, 0, None, 0).unwrap(),
            _IS_ERROR as i32
        );
        assert_eq!(
            GetInpStructErrorType(None, 0, None, -1).unwrap(),
            _IS_ERROR as i32
        );

        let mut parameters = INPUT_PARMS::default();
        assert_eq!(
            GetInpStructErrorType(Some(&parameters), 98, None, 0).unwrap(),
            _IS_ERROR as i32
        );
        parameters.bAllowEmptyStructure = 1;
        assert_eq!(
            GetInpStructErrorType(Some(&parameters), 98, None, 0).unwrap(),
            _IS_WARNING as i32
        );
        assert_eq!(
            GetInpStructErrorType(None, 98, None, 0),
            Err(SourceHeapError::NullPointer)
        );

        assert_eq!(
            GetInpStructErrorType(None, 0, Some(&[0]), 1).unwrap(),
            _IS_OKAY as i32
        );
        assert_eq!(
            GetInpStructErrorType(None, 0, Some(&[b'x' as i8, 0]), 1).unwrap(),
            _IS_WARNING as i32
        );
        assert_eq!(
            GetInpStructErrorType(None, 0, None, 1),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichimake__comparedfsdescendants4ct__line_2119() {
        let dfs = [1_u16, 2, 3, 4];
        let descendants = [4_u16, 3, 1, 2];
        let deleted = (MAX_ATOMS + 1) as u16;
        let mut order = OrderStruct {
            dfs_number: &dfs,
            number_of_descendants: &descendants,
            current_atom: 0,
        };

        assert_eq!(
            CompareDfsDescendants4CT(deleted, deleted + 1, &order),
            Ok(0)
        );
        assert_eq!(CompareDfsDescendants4CT(deleted, 1, &order), Ok(1));
        assert_eq!(CompareDfsDescendants4CT(1, deleted, &order), Ok(-1));
        assert_eq!(CompareDfsDescendants4CT(1, 2, &order), Ok(2));
        assert_eq!(CompareDfsDescendants4CT(2, 3, &order), Ok(-1));

        order.current_atom = 3;
        assert_eq!(CompareDfsDescendants4CT(0, 1, &order), Ok(-1));
        assert_eq!(CompareDfsDescendants4CT(2, 1, &order), Ok(1));
        order.current_atom = -1;
        assert_eq!(
            CompareDfsDescendants4CT(0, 1, &order),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichimake__getdfsorder4ct__line_2154() {
        let mut heap = SourceHeap::default();
        let single = heap.allocate_model_storage(vec![1_u16]).unwrap();
        let output =
            GetDfsOrder4CT(&mut heap, single, 1, SourceConstPointer::null(), 1, 0).unwrap();
        assert_eq!(heap.slice(output.as_const()).unwrap(), &[1, 0, 0, 0, 0, 0]);
        inchi_free(&mut heap, output).unwrap();

        let chain = heap.allocate_model_storage(vec![2_u16, 1, 3, 2]).unwrap();
        let hydrogens = heap.allocate_model_storage(vec![1_i8, 2, 3]).unwrap();
        let output = GetDfsOrder4CT(&mut heap, chain, 4, hydrogens.as_const(), 3, 0).unwrap();
        assert_eq!(
            heap.slice(output.as_const()).unwrap(),
            &[1, 17, 0, 2, 18, b'-' as u16, 3, 19, b'-' as u16, 0, 0, 0]
        );
        inchi_free(&mut heap, output).unwrap();

        let output = GetDfsOrder4CT(
            &mut heap,
            chain,
            4,
            hydrogens.as_const(),
            3,
            crate::source_types::CT_MODE_PREDECESSORS as i32,
        )
        .unwrap();
        assert_eq!(
            heap.slice(output.as_const()).unwrap(),
            &[
                (MAX_ATOMS + 1) as u16,
                17,
                0,
                1,
                18,
                b',' as u16,
                2,
                19,
                b',' as u16,
                0,
                0,
                0,
            ]
        );
        inchi_free(&mut heap, output).unwrap();

        for failure_after in 0..8 {
            let mut failing_heap = SourceHeap::default();
            let chain = failing_heap
                .allocate_model_storage(vec![2_u16, 1, 3, 2])
                .unwrap();
            failing_heap.fail_after_allocations(failure_after);
            let output = GetDfsOrder4CT(
                &mut failing_heap,
                chain,
                4,
                SourceConstPointer::null(),
                3,
                0,
            )
            .unwrap();
            assert!(output.is_null(), "allocation failure {failure_after}");
        }
    }

    #[test]
    fn source_port__ichimake__inchi_segmentaction__line_1486() {
        for value in i8::MIN..=i8::MAX {
            let promoted = i32::from(value);
            let expected = if promoted & crate::source_types::tagMarkDiff_DIFV_OUTPUT_OMIT_F as i32
                == 0
            {
                crate::source_types::tagINChISegmAction_INCHI_SEGM_OMIT as i32
            } else if promoted & crate::source_types::tagMarkDiff_DIFV_OUTPUT_EMPTY_T as i32 != 0
                && promoted & crate::source_types::tagMarkDiff_DIFV_OUTPUT_EMPTY_F as i32 == 0
            {
                crate::source_types::tagINChISegmAction_INCHI_SEGM_EMPTY as i32
            } else if promoted & crate::source_types::tagMarkDiff_DIFV_OUTPUT_FILL_T as i32 != 0 {
                crate::source_types::tagINChISegmAction_INCHI_SEGM_FILL as i32
            } else {
                crate::source_types::tagINChISegmAction_INCHI_SEGM_OMIT as i32
            };
            assert_eq!(INChI_SegmentAction(value), expected, "value={value}");
        }
    }
}
