use crate::source::base::ichicano::{InchiTimeElapsed, InchiTimeGet};
use crate::source::base::ichimake::{CompINChITaut2, CompareHillFormulasNoH, CompareReversedINChI};
use crate::source::base::ichiread::{
    bInpInchiComponentDeleted, bInpInchiComponentExists, bRevInchiComponentDeleted,
    bRevInchiComponentExists, insertions_sort_AT_NUMB,
};
use crate::source::base::ichirvr4::{
    AddRemIsoProtonsInRestrStruct, AddRemProtonsInRestrStruct, OneInChI2Atom,
};
use crate::source::base::ichitaut::free_t_group_info;
use crate::source::base::mol2atom::FreeExtOrigAtData;
use crate::source::base::runichi4::FreeAllINChIArrays;
use crate::source::base::strutil::Free_INChI_Members;
use crate::source::base::util::{inchi_calloc, inchi_free, inchi_malloc};
use crate::source_types::{
    _IS_OKAY, _IS_WARNING, AB_PARITY_UNDF, AT_NUMB, CANON_GLOBALS, COMPONENT_REM_PROTONS,
    CT_USER_QUIT_ERR, I2A_FLAG_FIXEDH, I2A_FLAG_RECMET, ICR, ICR_MAX_DIFF_FIXED_H,
    ICR_MAX_ENDP_IN1_ONLY, ICR_MAX_ENDP_IN2_ONLY, ICR_MAX_SB_IN1_ONLY, ICR_MAX_SB_IN2_ONLY,
    ICR_MAX_SB_UNDF, ICR_MAX_SC_IN1_ONLY, ICR_MAX_SC_IN2_ONLY, ICR_MAX_SC_UNDF, INCHI_BAS,
    INCHI_CLOCK, INCHI_MODE, INCHI_NUM, INCHI_REC, INCHI_SORT, INChI, INChI_Aux, INChI_Stereo,
    INPUT_PARMS, InpInChI, MAX_NUM_STEREO_ATOM_NEIGH, MAX_NUM_STEREO_BONDS, NO_VALUE_INT,
    NUM_H_ISOTOPES, REQ_MODE_BASIC, RI_ERR_ALLOC, RI_ERR_PROGR, SRM, STRUCT_DATA,
    SourceConstPointer, SourceHeap, SourceHeapError, SourceMutPointer, StrFromINChI, TAUT_NON,
    TAUT_NUM, TAUT_YES, clock_t, inchiTime, inp_ATOM, tagInchiCompareDiffBits_INCHIDIFF_ATOMS,
    tagInchiCompareDiffBits_INCHIDIFF_CHARGE, tagInchiCompareDiffBits_INCHIDIFF_COMP_HLAYER,
    tagInchiCompareDiffBits_INCHIDIFF_COMP_NUMBER, tagInchiCompareDiffBits_INCHIDIFF_CON_LEN,
    tagInchiCompareDiffBits_INCHIDIFF_CON_TBL, tagInchiCompareDiffBits_INCHIDIFF_DIFF_TG_ENDP,
    tagInchiCompareDiffBits_INCHIDIFF_EXTRA_TG_ENDP, tagInchiCompareDiffBits_INCHIDIFF_ISO_AT,
    tagInchiCompareDiffBits_INCHIDIFF_LESS_FH, tagInchiCompareDiffBits_INCHIDIFF_LESS_H,
    tagInchiCompareDiffBits_INCHIDIFF_MISS_TG_ENDP, tagInchiCompareDiffBits_INCHIDIFF_MOB_ISO_H,
    tagInchiCompareDiffBits_INCHIDIFF_MOBH_PROTONS, tagInchiCompareDiffBits_INCHIDIFF_MORE_FH,
    tagInchiCompareDiffBits_INCHIDIFF_MORE_H, tagInchiCompareDiffBits_INCHIDIFF_MULTIPLE_TG,
    tagInchiCompareDiffBits_INCHIDIFF_NO_TAUT, tagInchiCompareDiffBits_INCHIDIFF_NUM_AT,
    tagInchiCompareDiffBits_INCHIDIFF_NUM_EL, tagInchiCompareDiffBits_INCHIDIFF_NUM_ISO_AT,
    tagInchiCompareDiffBits_INCHIDIFF_NUM_TG, tagInchiCompareDiffBits_INCHIDIFF_POSITION_H,
    tagInchiCompareDiffBits_INCHIDIFF_PROBLEM, tagInchiCompareDiffBits_INCHIDIFF_REM_ISO_H,
    tagInchiCompareDiffBits_INCHIDIFF_REM_PROT, tagInchiCompareDiffBits_INCHIDIFF_SB_EXTRA,
    tagInchiCompareDiffBits_INCHIDIFF_SB_EXTRA_UNDF, tagInchiCompareDiffBits_INCHIDIFF_SB_MISS,
    tagInchiCompareDiffBits_INCHIDIFF_SB_MISS_UNDF, tagInchiCompareDiffBits_INCHIDIFF_SB_PARITY,
    tagInchiCompareDiffBits_INCHIDIFF_SC_EXTRA, tagInchiCompareDiffBits_INCHIDIFF_SC_EXTRA_UNDF,
    tagInchiCompareDiffBits_INCHIDIFF_SC_INV, tagInchiCompareDiffBits_INCHIDIFF_SC_MISS,
    tagInchiCompareDiffBits_INCHIDIFF_SC_MISS_UNDF, tagInchiCompareDiffBits_INCHIDIFF_SC_PARITY,
    tagInchiCompareDiffBits_INCHIDIFF_SINGLE_TG, tagInchiCompareDiffBits_INCHIDIFF_STR2INCHI_ERR,
    tagInchiCompareDiffBits_INCHIDIFF_TG, tagInchiCompareDiffBits_INCHIDIFF_WRONG_TAUT,
    tagtagCompareInchiMsgGroupID_IDGRP_CHARGE, tagtagCompareInchiMsgGroupID_IDGRP_COMP,
    tagtagCompareInchiMsgGroupID_IDGRP_CONV_ERR, tagtagCompareInchiMsgGroupID_IDGRP_ERR,
    tagtagCompareInchiMsgGroupID_IDGRP_H, tagtagCompareInchiMsgGroupID_IDGRP_HLAYER,
    tagtagCompareInchiMsgGroupID_IDGRP_ISO_AT, tagtagCompareInchiMsgGroupID_IDGRP_ISO_H,
    tagtagCompareInchiMsgGroupID_IDGRP_MOB_GRP, tagtagCompareInchiMsgGroupID_IDGRP_PROTONS,
    tagtagCompareInchiMsgGroupID_IDGRP_SB, tagtagCompareInchiMsgGroupID_IDGRP_SC,
    tagtagCompareInchiMsgGroupID_IDGRP_ZERO,
};

// BEGIN INCHI ACTIVE MACRO CONFIGURATION: util.h pseudo-element numbers
// INCHI✔❌: #define EL_NUMBER_ZY ((U_CHAR)119)
// INCHI✔❌: #define EL_NUMBER_ZZ ((U_CHAR)120)
// END INCHI ACTIVE MACRO CONFIGURATION: util.h pseudo-element numbers
const EL_NUMBER_ZY: u8 = 119;
const EL_NUMBER_ZZ: u8 = 120;

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn RemoveFixHInChIIdentical2MobH(
    heap: &mut SourceHeap,
    pOneInput: &mut InpInChI,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:256 RemoveFixHInChIIdentical2MobH
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
void RemoveFixHInChIIdentical2MobH(InpInChI* pOneInput)
{
    int iInchiRec, cur_num_comp, k;

    /* eliminate Fixed-H InChI that are exactly same as the corresponding Mobile-H structures */
    for (iInchiRec = 0; iInchiRec < INCHI_NUM; iInchiRec++)
    {
        cur_num_comp = inchi_min(pOneInput->nNumComponents[iInchiRec][TAUT_YES],
            pOneInput->nNumComponents[iInchiRec][TAUT_NON]);
        for (k = 0; k < cur_num_comp; k++)
        {
            if (!CompareReversedINChI(pOneInput->pInpInChI[iInchiRec][TAUT_YES] + k,
                pOneInput->pInpInChI[iInchiRec][TAUT_NON] + k, NULL, NULL))
            {
                Free_INChI_Members(pOneInput->pInpInChI[iInchiRec][TAUT_NON] + k);
                memset(pOneInput->pInpInChI[iInchiRec][TAUT_NON] + k, 0, sizeof(pOneInput->pInpInChI[0][0][0])); /* djb-rwth: memset_s C11/Annex K variant? */
            }
        }
    }
}
    */
    // END INCHI C FUNCTION: RemoveFixHInChIIdentical2MobH
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: RemoveFixHInChIIdentical2MobH
    // INCHI✔️❌: READ_INCHI_STRING=1 includes this production helper.
    // INCHI✔️❌: inchi_min is reproduced as the active header macro.
    // INCHI✔️❌: COMPILE_ANSI_ONLY and TARGET_API_LIB do not alter this function body.
    // END INCHI ACTIVE MACRO CONFIGURATION: RemoveFixHInChIIdentical2MobH

    for representation in 0..INCHI_NUM as usize {
        let component_count = pOneInput.nNumComponents[representation][TAUT_YES as usize]
            .min(pOneInput.nNumComponents[representation][TAUT_NON as usize]);
        for component in 0..component_count.max(0) {
            let mobile_base = pOneInput.pInpInChI[representation][TAUT_YES as usize];
            let fixed_base = pOneInput.pInpInChI[representation][TAUT_NON as usize];
            if mobile_base.is_null() || fixed_base.is_null() {
                return Err(SourceHeapError::NullPointer);
            }
            let mobile_pointer = mobile_base.offset(i64::from(component))?;
            let fixed_pointer = fixed_base.offset(i64::from(component))?;
            let mobile = heap
                .slice(mobile_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let fixed = heap
                .slice(fixed_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            if CompareReversedINChI(heap, Some(&mobile), Some(&fixed), None, None)? == 0 {
                Free_INChI_Members(heap, fixed_pointer)?;
                *heap
                    .slice_mut(fixed_pointer)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = INChI::default();
            }
        }
    }
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MarkDisconectedIdenticalToReconnected(
    heap: &mut SourceHeap,
    pOneInput: &mut InpInChI,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:283 MarkDisconectedIdenticalToReconnected
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int MarkDisconectedIdenticalToReconnected(InpInChI* pOneInput)
{
    int k1, k2, num_marked = 0;
    int k1max = inchi_max(pOneInput->nNumComponents[INCHI_BAS][TAUT_YES],
        pOneInput->nNumComponents[INCHI_BAS][TAUT_NON]);

    for (k1 = 0; k1 < k1max; k1++)
    {
        int k2max = inchi_max(pOneInput->nNumComponents[INCHI_REC][TAUT_YES],
            pOneInput->nNumComponents[INCHI_REC][TAUT_NON]);

        for (k2 = 0; k2 < k2max; k2++)
        {
            int eqM = k1 < pOneInput->nNumComponents[INCHI_BAS][TAUT_YES]
                &&
                k2 < pOneInput->nNumComponents[INCHI_REC][TAUT_YES]
                &&
                !pOneInput->pInpInChI[INCHI_REC][TAUT_YES][k2].nLink /* already linked */
                &&
                !pOneInput->pInpInChI[INCHI_BAS][TAUT_YES][k1].bDeleted
                &&
                pOneInput->pInpInChI[INCHI_BAS][TAUT_YES][k1].nNumberOfAtoms
                &&
                pOneInput->pInpInChI[INCHI_BAS][TAUT_YES][k1].nNumberOfAtoms ==
                pOneInput->pInpInChI[INCHI_REC][TAUT_YES][k2].nNumberOfAtoms
                &&
                !pOneInput->pInpInChI[INCHI_REC][TAUT_YES][k2].bDeleted
                &&
                !CompareReversedINChI(pOneInput->pInpInChI[INCHI_REC][TAUT_YES] + k2,
                    pOneInput->pInpInChI[INCHI_BAS][TAUT_YES] + k1,
                    NULL, NULL);

            int isF1 = k1 < pOneInput->nNumComponents[INCHI_BAS][TAUT_NON]
                &&
                0 == pOneInput->pInpInChI[INCHI_BAS][TAUT_NON][k1].bDeleted
                &&
                0 < pOneInput->pInpInChI[INCHI_BAS][TAUT_NON][k1].nNumberOfAtoms
                ;

            int isF2 = k2 < pOneInput->nNumComponents[INCHI_REC][TAUT_NON]
                &&
                0 == pOneInput->pInpInChI[INCHI_REC][TAUT_NON][k2].bDeleted
                &&
                0 < pOneInput->pInpInChI[INCHI_REC][TAUT_NON][k2].nNumberOfAtoms
                ;

            int eqF = isF1
                &&
                isF2
                &&
                !pOneInput->pInpInChI[INCHI_REC][TAUT_NON][k2].nLink
                &&
                pOneInput->pInpInChI[INCHI_BAS][TAUT_NON][k1].nNumberOfAtoms ==
                pOneInput->pInpInChI[INCHI_REC][TAUT_NON][k2].nNumberOfAtoms
                &&
                !CompareReversedINChI(pOneInput->pInpInChI[INCHI_REC][TAUT_NON] + k2,
                    pOneInput->pInpInChI[INCHI_BAS][TAUT_NON] + k1,
                    NULL, NULL)
                ;

            if (eqM && ((!isF1 && !isF2) || eqF)) /* djb-rwth: addressing LLVM warning */
            {
                pOneInput->pInpInChI[INCHI_BAS][TAUT_YES][k1].nLink = -(k2 + 1);
                pOneInput->pInpInChI[INCHI_REC][TAUT_YES][k2].nLink = (k1 + 1);
                if (eqF)
                {
                    pOneInput->pInpInChI[INCHI_BAS][TAUT_NON][k1].nLink = -(k2 + 1);
                    pOneInput->pInpInChI[INCHI_REC][TAUT_NON][k2].nLink = (k1 + 1);
                }
                num_marked++;
                break;
                /* equal InChI has been deleted from the disconnected layer, get next k1 */
            }
        }
    }

    return num_marked;
}
    */
    // END INCHI C FUNCTION: MarkDisconectedIdenticalToReconnected
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: MarkDisconectedIdenticalToReconnected
    // INCHI✔️❌: READ_INCHI_STRING=1 includes this production helper.
    // INCHI✔️❌: inchi_max is reproduced as the active header macro.
    // INCHI✔️❌: COMPILE_ANSI_ONLY and TARGET_API_LIB do not alter this function body.
    // END INCHI ACTIVE MACRO CONFIGURATION: MarkDisconectedIdenticalToReconnected

    let pointer = |input: &InpInChI,
                   representation: usize,
                   mobile_h: usize,
                   component: i32|
     -> Result<SourceMutPointer<INChI>, SourceHeapError> {
        let base = input.pInpInChI[representation][mobile_h];
        if base.is_null() {
            return Err(SourceHeapError::NullPointer);
        }
        base.offset(i64::from(component))
    };
    let disconnected_count = pOneInput.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize]
        .max(pOneInput.nNumComponents[INCHI_BAS as usize][TAUT_NON as usize]);
    let mut num_marked = 0_i32;
    for disconnected in 0..disconnected_count.max(0) {
        let reconnected_count = pOneInput.nNumComponents[INCHI_REC as usize][TAUT_YES as usize]
            .max(pOneInput.nNumComponents[INCHI_REC as usize][TAUT_NON as usize]);
        for reconnected in 0..reconnected_count.max(0) {
            let mobile1 = if disconnected
                < pOneInput.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize]
            {
                Some(pointer(
                    pOneInput,
                    INCHI_BAS as usize,
                    TAUT_YES as usize,
                    disconnected,
                )?)
            } else {
                None
            };
            let mobile2 = if reconnected
                < pOneInput.nNumComponents[INCHI_REC as usize][TAUT_YES as usize]
            {
                Some(pointer(
                    pOneInput,
                    INCHI_REC as usize,
                    TAUT_YES as usize,
                    reconnected,
                )?)
            } else {
                None
            };
            let mut eq_mobile = false;
            if let (Some(mobile1_pointer), Some(mobile2_pointer)) = (mobile1, mobile2) {
                let mobile1_value = heap
                    .slice(mobile1_pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone();
                let mobile2_value = heap
                    .slice(mobile2_pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone();
                if mobile2_value.nLink == 0
                    && mobile1_value.bDeleted == 0
                    && mobile1_value.nNumberOfAtoms != 0
                    && mobile1_value.nNumberOfAtoms == mobile2_value.nNumberOfAtoms
                    && mobile2_value.bDeleted == 0
                    && CompareReversedINChI(
                        heap,
                        Some(&mobile2_value),
                        Some(&mobile1_value),
                        None,
                        None,
                    )? == 0
                {
                    eq_mobile = true;
                }
            }

            let fixed1 = if disconnected
                < pOneInput.nNumComponents[INCHI_BAS as usize][TAUT_NON as usize]
            {
                Some(pointer(
                    pOneInput,
                    INCHI_BAS as usize,
                    TAUT_NON as usize,
                    disconnected,
                )?)
            } else {
                None
            };
            let fixed1_value = if let Some(fixed1) = fixed1 {
                Some(
                    heap.slice(fixed1.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone(),
                )
            } else {
                None
            };
            let is_fixed1 = fixed1_value
                .as_ref()
                .is_some_and(|value| value.bDeleted == 0 && value.nNumberOfAtoms > 0);
            let fixed2 = if reconnected
                < pOneInput.nNumComponents[INCHI_REC as usize][TAUT_NON as usize]
            {
                Some(pointer(
                    pOneInput,
                    INCHI_REC as usize,
                    TAUT_NON as usize,
                    reconnected,
                )?)
            } else {
                None
            };
            let fixed2_value = if let Some(fixed2) = fixed2 {
                Some(
                    heap.slice(fixed2.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone(),
                )
            } else {
                None
            };
            let is_fixed2 = fixed2_value
                .as_ref()
                .is_some_and(|value| value.bDeleted == 0 && value.nNumberOfAtoms > 0);
            let eq_fixed = if is_fixed1 && is_fixed2 {
                let fixed1_value = fixed1_value
                    .as_ref()
                    .ok_or(SourceHeapError::NullPointer)?;
                let fixed2_value = fixed2_value
                    .as_ref()
                    .ok_or(SourceHeapError::NullPointer)?;
                fixed2_value.nLink == 0
                    && fixed1_value.nNumberOfAtoms == fixed2_value.nNumberOfAtoms
                    && CompareReversedINChI(
                        heap,
                        Some(fixed2_value),
                        Some(fixed1_value),
                        None,
                        None,
                    )? == 0
            } else {
                false
            };
            if eq_mobile && ((!is_fixed1 && !is_fixed2) || eq_fixed) {
                let disconnected_link = reconnected.wrapping_add(1).wrapping_neg();
                let reconnected_link = disconnected.wrapping_add(1);
                heap.slice_mut(mobile1.ok_or(SourceHeapError::NullPointer)?)?[0].nLink =
                    disconnected_link;
                heap.slice_mut(mobile2.ok_or(SourceHeapError::NullPointer)?)?[0].nLink =
                    reconnected_link;
                if eq_fixed {
                    heap.slice_mut(fixed1.ok_or(SourceHeapError::NullPointer)?)?[0].nLink =
                        disconnected_link;
                    heap.slice_mut(fixed2.ok_or(SourceHeapError::NullPointer)?)?[0].nLink =
                        reconnected_link;
                }
                num_marked = num_marked.wrapping_add(1);
                break;
            }
        }
    }
    Ok(num_marked)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn AddOneMsg(
    szMsg: &mut [i8],
    used_len: i32,
    tot_len: i32,
    szAddMsg: &[i8],
    szDelim: Option<&[i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:3138 AddOneMsg
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int AddOneMsg(char* szMsg,
    int used_len,
    int tot_len,
    const char* szAddMsg,
    const char* szDelim)
{
    const char ellip[] = "...";
    int len = (int)strlen(szAddMsg);
    int len_delim = (used_len && szDelim) ? strlen(szDelim) : 0;
    int len_to_copy;
    if (len + len_delim + used_len < tot_len)
    {
        if (len_delim)
        {
            strcpy(szMsg + used_len, szDelim);
            used_len += len_delim;
        }
        strcpy(szMsg + used_len, szAddMsg);
        used_len += len;
    }
    else
    {
        if ((len_to_copy = (tot_len - used_len - len_delim - (int)sizeof(ellip))) > 10)
        {
            if (len_delim)
            {
                strcpy(szMsg + used_len, szDelim);
                used_len += len_delim;
            }
            strncpy(szMsg + used_len, szAddMsg, len_to_copy);
            used_len += len_to_copy;
            strcpy(szMsg + used_len, ellip);
            used_len += sizeof(ellip) - 1;
        }
    }

    return used_len;
}
    */
    // END INCHI C FUNCTION: AddOneMsg
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: AddOneMsg
    // INCHI✔️❌: TARGET_API_LIB and COMPILE_ANSI_ONLY do not alter this function body.
    // INCHI✔️❌: sizeof(ellip) is 4, including its terminating NUL.
    // END INCHI ACTIVE MACRO CONFIGURATION: AddOneMsg

    let c_len = |value: &[i8]| -> Result<i32, SourceHeapError> {
        let length = value
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        i32::try_from(length).map_err(|_| SourceHeapError::SourceIntegerOverflow)
    };
    let copy_c_string = |target: &mut [i8], offset: i32, value: &[i8]| -> Result<(), SourceHeapError> {
        let start = usize::try_from(offset).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let length = usize::try_from(c_len(value)?)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let end = start
            .checked_add(length)
            .and_then(|end| end.checked_add(1))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let destination = target
            .get_mut(start..end)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        destination.copy_from_slice(&value[..=length]);
        Ok(())
    };

    let len = c_len(szAddMsg)?;
    let len_delim = if used_len != 0 {
        szDelim.map(c_len).transpose()?.unwrap_or(0)
    } else {
        0
    };
    let mut used_len = used_len;
    if len.wrapping_add(len_delim).wrapping_add(used_len) < tot_len {
        if len_delim != 0 {
            copy_c_string(
                szMsg,
                used_len,
                szDelim.ok_or(SourceHeapError::NullPointer)?,
            )?;
            used_len = used_len.wrapping_add(len_delim);
        }
        copy_c_string(szMsg, used_len, szAddMsg)?;
        used_len = used_len.wrapping_add(len);
    } else {
        let len_to_copy = tot_len
            .wrapping_sub(used_len)
            .wrapping_sub(len_delim)
            .wrapping_sub(4);
        if len_to_copy > 10 {
            if len_delim != 0 {
                copy_c_string(
                    szMsg,
                    used_len,
                    szDelim.ok_or(SourceHeapError::NullPointer)?,
                )?;
                used_len = used_len.wrapping_add(len_delim);
            }
            let start = usize::try_from(used_len)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let count = usize::try_from(len_to_copy)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let end = start
                .checked_add(count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let destination = szMsg
                .get_mut(start..end)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let source_length = usize::try_from(len)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            for (index, byte) in destination.iter_mut().enumerate() {
                *byte = if index < source_length {
                    szAddMsg[index]
                } else {
                    0
                };
            }
            used_len = used_len.wrapping_add(len_to_copy);
            copy_c_string(szMsg, used_len, &[b'.' as i8, b'.' as i8, b'.' as i8, 0])?;
            used_len = used_len.wrapping_add(3);
        }
    }
    Ok(used_len)
}

// BEGIN INCHI C DATA FRAME: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:3074 CompareInchiMsgsGroup/CompareInchiMsgs
// INCHI✔️❌: complete source frame follows verbatim.
/*
/* message group names */
const CMP_INCHI_MSG_GROUP CompareInchiMsgsGroup[] =
{
    {IDGRP_ERR,     " Error:"},
    {IDGRP_H,       " Hydrogens:"},
    {IDGRP_MOB_GRP, " Mobile-H groups:"},
    {IDGRP_ISO_AT,  " Isotopic:"},
    {IDGRP_CHARGE,  " Charge(s):"},
    {IDGRP_PROTONS, " Proton balance:"},
    {IDGRP_ISO_H,   " Exchangeable isotopic H:"},
    {IDGRP_SC,      " Stereo centers/allenes:"},
    {IDGRP_SB,      " Stereobonds/cumulenes:"},
    {IDGRP_HLAYER,  " Fixed-H layer:"},
    {IDGRP_COMP,    " Number of components:"},
    {IDGRP_CONV_ERR," Conversion encountered:"},
    {IDGRP_ZERO,    ""}
};



/* messages */
const CMP_INCHI_MSG  CompareInchiMsgs[] =
{
    {INCHIDIFF_PROBLEM      ,IDGRP_ERR,     " Wrong result"                   }, /*0x00000001,  severe: at least one InChI does not exist */
    {INCHIDIFF_POSITION_H   ,IDGRP_H,       " Locations or number"            }, /*0x00000002,  difference in non-taut {Mobile-H} or all H {Fixed-H} location/number */
    {INCHIDIFF_MORE_FH      ,IDGRP_H,       " Fixed-H"                        }, /*0x00000004,  extra fixed H */
    {INCHIDIFF_LESS_FH      ,IDGRP_H,       " Fixed-H"                        }, /*0x00000004,  missing fixed H */
    {INCHIDIFF_MORE_H       ,IDGRP_H,       " Number"                         }, /*0x00000008,  formulas differ in number of H */
    {INCHIDIFF_LESS_H       ,IDGRP_H,       " Number"                         }, /*0x00000008,  formulas differ in number of H */
    {INCHIDIFF_NO_TAUT      ,IDGRP_MOB_GRP, " Missing"                        }, /*0x00000010,  restored structure has no taut groups while the original InChI has some */
    {INCHIDIFF_WRONG_TAUT   ,IDGRP_MOB_GRP, " Falsely present"                }, /*0x00000020,  restored has tautomerism while the original does not have it */
    {INCHIDIFF_SINGLE_TG    ,IDGRP_MOB_GRP, " One instead of multiple"        }, /*0x00000040,  restored has 1 taut. group while the original InChI has multiple tg */
    {INCHIDIFF_MULTIPLE_TG  ,IDGRP_MOB_GRP, " Multiple instead of one"        }, /*0x00000080,  restored has multiple tg while the original InChI has only one tg */
    {INCHIDIFF_EXTRA_TG_ENDP,IDGRP_MOB_GRP, " Attachment points"              }, /*0x00000100,  extra tautomeric endpoint{s} in restored structure */
    {INCHIDIFF_MISS_TG_ENDP ,IDGRP_MOB_GRP, " Attachment points"              }, /*0x00000100,  one or more tg endpoint is not in the restored structure */
    {INCHIDIFF_DIFF_TG_ENDP ,IDGRP_MOB_GRP, " Attachment points"              }, /*0x00000100,  lists of tg endpoints are different */
    {INCHIDIFF_NUM_TG       ,IDGRP_MOB_GRP, " Number"                         }, /*0x00000200,  different number of tautomeric groups */
    {INCHIDIFF_TG           ,IDGRP_MOB_GRP, " Do not match"                   }, /*0x00000200,  different tautomeric groups */
    {INCHIDIFF_NUM_ISO_AT   ,IDGRP_ISO_AT,  " Atoms do not match"             }, /*0x00000400,  ?severe: restored struct. has different number of isotopic atoms */
    {INCHIDIFF_ISO_AT       ,IDGRP_ISO_AT,  " Atoms do not match"             }, /*0x00000400,  ?severe: restored struct. has different locations/isotopes of isotopic atoms */
    {INCHIDIFF_REM_ISO_H    ,IDGRP_ISO_H,   " Does not match for a component" }, /*0x00000800,  isotopic H removed */
    {INCHIDIFF_MOB_ISO_H    ,IDGRP_ISO_H,   " Do not match"                   }, /*0x00001000,  different number of mobile exchangeable isotopic H */
    {INCHIDIFF_CHARGE       ,IDGRP_CHARGE,  " Do not match"                   }, /*0x00002000,  restored structure has different charge */
    {INCHIDIFF_REM_PROT     ,IDGRP_PROTONS, " Does not match for a component" }, /*0x00004000,  proton{s} removed/added from the restored structure */
    {INCHIDIFF_MOBH_PROTONS ,IDGRP_PROTONS, " Does not match"                 }, /*0x00008000,  different proton balance */
    {INCHIDIFF_SC_INV       ,IDGRP_SC,      " Falsely inverted"               }, /*0x00010000,  restores structure has different inversion stereocenter mark */
    {INCHIDIFF_SC_PARITY    ,IDGRP_SC,      " Wrong parity"                   }, /*0x00020000,  restored structure has stereoatoms or allenes with different parity */
    {INCHIDIFF_SC_EXTRA_UNDF,IDGRP_SC,      " Extra undefined"                }, /*0x00040000,  restored structure has extra undefined stereocenter{s} */
    {INCHIDIFF_SC_EXTRA     ,IDGRP_SC,      " Extra known"                    }, /*0x00080000,  restored structure has extra stereocenter{s} */
    {INCHIDIFF_SC_MISS_UNDF ,IDGRP_SC,      " Missing undefined"              }, /*0x00100000,  restored structure has not some undefined stereocenter{s} */
    {INCHIDIFF_SC_MISS      ,IDGRP_SC,      " Missing known"                  }, /*0x00200000,  restored structure has not some stereocenters that are not undefined */
    {INCHIDIFF_SB_PARITY    ,IDGRP_SB,      " Wrong parity"                   }, /*0x00400000,  restored structure has stereobonds or cumulenes with different parity */
    {INCHIDIFF_SB_EXTRA_UNDF,IDGRP_SB,      " Extra undefined"                }, /*0x00800000,  restored structure has extra undefined stereobond{s} */
    {INCHIDIFF_SB_EXTRA     ,IDGRP_SB,      " Missing known"                  }, /*0x01000000,  restored structure has extra stereobond{s} */
    {INCHIDIFF_SB_MISS_UNDF ,IDGRP_SB,      " Missing undefined"              }, /*0x02000000,  restored structure has not some undefined stereocenters */
    {INCHIDIFF_SB_MISS      ,IDGRP_SB,      " Missing known"                  }, /*0x04000000,  restored structure has not some stereobonds that are not undefined */
    {INCHIDIFF_COMP_HLAYER  ,IDGRP_HLAYER,  " Missing or extra"               }, /*0x08000000,  Restored component has Mobile-H layer instead of both Mobile-H & Fixed-H or both instead of one */
    {INCHIDIFF_COMP_NUMBER  ,IDGRP_COMP,    " Does not match"                 }, /*0x10000000,  wrong number of components */
    {INCHIDIFF_STR2INCHI_ERR,IDGRP_CONV_ERR," Error"                          },  /*0x20000000   Restored structure to InChI conversion error */
    {INCHIDIFF_ZERO         ,IDGRP_ZERO,    ""                                }
};
*/
// END INCHI C DATA FRAME: CompareInchiMsgsGroup/CompareInchiMsgs
const COMPARE_INCHI_MESSAGE_GROUPS: &[(u32, &str)] = &[
    (tagtagCompareInchiMsgGroupID_IDGRP_ERR, " Error:"),
    (tagtagCompareInchiMsgGroupID_IDGRP_H, " Hydrogens:"),
    (
        tagtagCompareInchiMsgGroupID_IDGRP_MOB_GRP,
        " Mobile-H groups:",
    ),
    (tagtagCompareInchiMsgGroupID_IDGRP_ISO_AT, " Isotopic:"),
    (tagtagCompareInchiMsgGroupID_IDGRP_CHARGE, " Charge(s):"),
    (
        tagtagCompareInchiMsgGroupID_IDGRP_PROTONS,
        " Proton balance:",
    ),
    (
        tagtagCompareInchiMsgGroupID_IDGRP_ISO_H,
        " Exchangeable isotopic H:",
    ),
    (
        tagtagCompareInchiMsgGroupID_IDGRP_SC,
        " Stereo centers/allenes:",
    ),
    (
        tagtagCompareInchiMsgGroupID_IDGRP_SB,
        " Stereobonds/cumulenes:",
    ),
    (tagtagCompareInchiMsgGroupID_IDGRP_HLAYER, " Fixed-H layer:"),
    (
        tagtagCompareInchiMsgGroupID_IDGRP_COMP,
        " Number of components:",
    ),
    (
        tagtagCompareInchiMsgGroupID_IDGRP_CONV_ERR,
        " Conversion encountered:",
    ),
    (tagtagCompareInchiMsgGroupID_IDGRP_ZERO, ""),
];
const COMPARE_INCHI_MESSAGES: &[(u32, u32, &str)] = &[
    (
        tagInchiCompareDiffBits_INCHIDIFF_PROBLEM,
        tagtagCompareInchiMsgGroupID_IDGRP_ERR,
        " Wrong result",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_POSITION_H,
        tagtagCompareInchiMsgGroupID_IDGRP_H,
        " Locations or number",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_MORE_FH,
        tagtagCompareInchiMsgGroupID_IDGRP_H,
        " Fixed-H",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_LESS_FH,
        tagtagCompareInchiMsgGroupID_IDGRP_H,
        " Fixed-H",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_MORE_H,
        tagtagCompareInchiMsgGroupID_IDGRP_H,
        " Number",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_LESS_H,
        tagtagCompareInchiMsgGroupID_IDGRP_H,
        " Number",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_NO_TAUT,
        tagtagCompareInchiMsgGroupID_IDGRP_MOB_GRP,
        " Missing",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_WRONG_TAUT,
        tagtagCompareInchiMsgGroupID_IDGRP_MOB_GRP,
        " Falsely present",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_SINGLE_TG,
        tagtagCompareInchiMsgGroupID_IDGRP_MOB_GRP,
        " One instead of multiple",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_MULTIPLE_TG,
        tagtagCompareInchiMsgGroupID_IDGRP_MOB_GRP,
        " Multiple instead of one",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_EXTRA_TG_ENDP,
        tagtagCompareInchiMsgGroupID_IDGRP_MOB_GRP,
        " Attachment points",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_MISS_TG_ENDP,
        tagtagCompareInchiMsgGroupID_IDGRP_MOB_GRP,
        " Attachment points",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_DIFF_TG_ENDP,
        tagtagCompareInchiMsgGroupID_IDGRP_MOB_GRP,
        " Attachment points",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_NUM_TG,
        tagtagCompareInchiMsgGroupID_IDGRP_MOB_GRP,
        " Number",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_TG,
        tagtagCompareInchiMsgGroupID_IDGRP_MOB_GRP,
        " Do not match",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_NUM_ISO_AT,
        tagtagCompareInchiMsgGroupID_IDGRP_ISO_AT,
        " Atoms do not match",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_ISO_AT,
        tagtagCompareInchiMsgGroupID_IDGRP_ISO_AT,
        " Atoms do not match",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_REM_ISO_H,
        tagtagCompareInchiMsgGroupID_IDGRP_ISO_H,
        " Does not match for a component",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_MOB_ISO_H,
        tagtagCompareInchiMsgGroupID_IDGRP_ISO_H,
        " Do not match",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_CHARGE,
        tagtagCompareInchiMsgGroupID_IDGRP_CHARGE,
        " Do not match",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_REM_PROT,
        tagtagCompareInchiMsgGroupID_IDGRP_PROTONS,
        " Does not match for a component",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_MOBH_PROTONS,
        tagtagCompareInchiMsgGroupID_IDGRP_PROTONS,
        " Does not match",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_SC_INV,
        tagtagCompareInchiMsgGroupID_IDGRP_SC,
        " Falsely inverted",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_SC_PARITY,
        tagtagCompareInchiMsgGroupID_IDGRP_SC,
        " Wrong parity",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_SC_EXTRA_UNDF,
        tagtagCompareInchiMsgGroupID_IDGRP_SC,
        " Extra undefined",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_SC_EXTRA,
        tagtagCompareInchiMsgGroupID_IDGRP_SC,
        " Extra known",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_SC_MISS_UNDF,
        tagtagCompareInchiMsgGroupID_IDGRP_SC,
        " Missing undefined",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_SC_MISS,
        tagtagCompareInchiMsgGroupID_IDGRP_SC,
        " Missing known",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_SB_PARITY,
        tagtagCompareInchiMsgGroupID_IDGRP_SB,
        " Wrong parity",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_SB_EXTRA_UNDF,
        tagtagCompareInchiMsgGroupID_IDGRP_SB,
        " Extra undefined",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_SB_EXTRA,
        tagtagCompareInchiMsgGroupID_IDGRP_SB,
        " Missing known",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_SB_MISS_UNDF,
        tagtagCompareInchiMsgGroupID_IDGRP_SB,
        " Missing undefined",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_SB_MISS,
        tagtagCompareInchiMsgGroupID_IDGRP_SB,
        " Missing known",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_COMP_HLAYER,
        tagtagCompareInchiMsgGroupID_IDGRP_HLAYER,
        " Missing or extra",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_COMP_NUMBER,
        tagtagCompareInchiMsgGroupID_IDGRP_COMP,
        " Does not match",
    ),
    (
        tagInchiCompareDiffBits_INCHIDIFF_STR2INCHI_ERR,
        tagtagCompareInchiMsgGroupID_IDGRP_CONV_ERR,
        " Error",
    ),
    (0, tagtagCompareInchiMsgGroupID_IDGRP_ZERO, ""),
];

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn FillOutCompareMessage(
    szMsg: &mut [i8],
    nLenMsg: i32,
    bits: &[INCHI_MODE; TAUT_NUM as usize],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:3179 FillOutCompareMessage
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int FillOutCompareMessage(char* szMsg, int nLenMsg, INCHI_MODE bits[])
{
    int bMobileH, k, n, len = (int)strlen(szMsg);
    int iPrevGrpIdx, iCurGrpIdx, bFound;
    INCHI_MODE bit;
    static const char* hdr = " Problems/mismatches:";
    char szOneMsg[256];

    int mismatch = 0;

    if (bits[TAUT_YES] || bits[TAUT_NON])
    {

        mismatch = -1;

        if (!strstr(szMsg, hdr))
        {
            len = AddOneMsg(szMsg, len, nLenMsg, hdr, NULL);
        }

        for (bMobileH = TAUT_YES; 0 <= bMobileH; bMobileH--)
        {
            /*      bMobileH = TAUT_YES, TAUT_NON */

            if (bits[bMobileH])
            {
                strcpy(szOneMsg, bMobileH == TAUT_YES ? " Mobile-H(" : " Fixed-H(");
                len = AddOneMsg(szMsg, len, nLenMsg, szOneMsg, NULL);
            }

            bit = 1;
            iPrevGrpIdx = -1;

            do
            {
                if (bit & bits[bMobileH])
                {
                    /* search for the message */
                    bFound = 0;
                    for (k = 0; CompareInchiMsgs[k].nBit != INCHIDIFF_ZERO && !bFound; k++)
                    {
                        if (bit & (INCHI_MODE)CompareInchiMsgs[k].nBit)
                        {
                            /* message found */
                            for (n = 0; CompareInchiMsgsGroup[n].nGroupID != IDGRP_ZERO; n++)
                            {
                                if (CompareInchiMsgsGroup[n].nGroupID == CompareInchiMsgs[k].nGroupID)
                                {
                                    iCurGrpIdx = n;
                                    if (iCurGrpIdx != iPrevGrpIdx)
                                    {
                                        if (iPrevGrpIdx >= 0)
                                        {
                                            len = AddOneMsg(szMsg, len, nLenMsg, ";", NULL);
                                        }
                                        len = AddOneMsg(szMsg, len, nLenMsg, CompareInchiMsgsGroup[iCurGrpIdx].szGroupName, NULL);
                                    }
                                    len = AddOneMsg(szMsg, len, nLenMsg, CompareInchiMsgs[k].szMsg, iCurGrpIdx == iPrevGrpIdx ? "," : NULL);
                                    iPrevGrpIdx = iCurGrpIdx;
                                    bFound = 1;
                                    break;
                                }
                            }
                        }
                    }
                }
                bit <<= 1;
            } while (bit);

            if (bits[bMobileH])
            {
                len = AddOneMsg(szMsg, len, nLenMsg, ")", NULL);
            }

        }
    }

    return mismatch; /*len; */
}
    */
    // END INCHI C FUNCTION: FillOutCompareMessage
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FillOutCompareMessage
    // INCHI✔️❌: READ_INCHI_STRING=1 includes this production helper.
    // INCHI✔️❌: INCHI_MODE is unsigned 64-bit, so bit <<= 1 terminates after bit 63.
    // INCHI✔️❌: TARGET_API_LIB and COMPILE_ANSI_ONLY do not alter this function body.
    // END INCHI ACTIVE MACRO CONFIGURATION: FillOutCompareMessage

    let mut len = i32::try_from(
        szMsg
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    )
    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let c_text = |text: &str| -> Result<[i8; 256], SourceHeapError> {
        if text.len() >= 256 {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let mut value = [0_i8; 256];
        for (target, source) in value.iter_mut().zip(text.bytes()) {
            *target = source as i8;
        }
        Ok(value)
    };
    let header = " Problems/mismatches:";
    if bits[TAUT_YES as usize] != 0 || bits[TAUT_NON as usize] != 0 {
        let initial_length = usize::try_from(len)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let header_c = c_text(header)?;
        if !szMsg[..initial_length]
            .windows(header.len())
            .any(|window| window == &header_c[..header.len()])
        {
            len = AddOneMsg(szMsg, len, nLenMsg, &header_c, None)?;
        }
        for mobile_h in [TAUT_YES as usize, TAUT_NON as usize] {
            if bits[mobile_h] != 0 {
                let layer = c_text(if mobile_h == TAUT_YES as usize {
                    " Mobile-H("
                } else {
                    " Fixed-H("
                })?;
                len = AddOneMsg(szMsg, len, nLenMsg, &layer, None)?;
            }
            let mut bit = 1 as INCHI_MODE;
            let mut previous_group = None;
            loop {
                if bit & bits[mobile_h] != 0 {
                    if let Some((_, group_id, message)) = COMPARE_INCHI_MESSAGES
                        .iter()
                        .find(|(message_bit, _, _)| {
                            *message_bit != 0 && bit & INCHI_MODE::from(*message_bit) != 0
                        })
                    {
                        if let Some((group_index, (_, group_name))) = COMPARE_INCHI_MESSAGE_GROUPS
                            .iter()
                            .enumerate()
                            .find(|(_, (candidate, _))| *candidate == *group_id)
                        {
                            if previous_group != Some(group_index) {
                                if previous_group.is_some() {
                                    let semicolon = c_text(";")?;
                                    len = AddOneMsg(szMsg, len, nLenMsg, &semicolon, None)?;
                                }
                                let group = c_text(group_name)?;
                                len = AddOneMsg(szMsg, len, nLenMsg, &group, None)?;
                            }
                            let message = c_text(message)?;
                            let delimiter = if previous_group == Some(group_index) {
                                Some(c_text(",")?)
                            } else {
                                None
                            };
                            len = AddOneMsg(
                                szMsg,
                                len,
                                nLenMsg,
                                &message,
                                delimiter.as_ref().map(|value| value.as_slice()),
                            )?;
                            previous_group = Some(group_index);
                        }
                    }
                }
                bit = bit.wrapping_shl(1);
                if bit == 0 {
                    break;
                }
            }
            if bits[mobile_h] != 0 {
                let close = c_text(")")?;
                len = AddOneMsg(szMsg, len, nLenMsg, &close, None)?;
            }
        }
        Ok(-1)
    } else {
        Ok(0)
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn SetUpSrm(pSrm: &mut SRM) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:364 SetUpSrm
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
void SetUpSrm(SRM* pSrm)
{
    /* structure restore parms !!!!! */
    memset(pSrm, 0, sizeof(pSrm[0])); /* djb-rwth: memset_s C11/Annex K variant? */
    pSrm->bFixStereoBonds = FIX_STEREO_BOND_ORDER;
    pSrm->nMetal2EndpointMinBondOrder = 1;
    pSrm->nMetal2EndpointInitEdgeFlow = 0;

    if (METAL_FREE_CHARGE_VAL == 1)
    {
        pSrm->bMetalAddFlower = 1;
        /* the next 3 parameters: */
        /* 0, 0, 0 => all bonds 0, no init radical on metal */
        /* 0, 0, 1 => all bonds 0,    init radical on metal */
        /* 0, 1, 0 => wrong */
        /* 0, 1, 1 => all bonds 1, no init radical on metal */
        /* 1, 0, 1 => min bond order 1, all bonds to metal have order 1 */
        /* 1, 1, 0 => wrong */
        /* 1, 1, 1 => wrong */
        pSrm->nMetalMinBondOrder = 0;
        pSrm->nMetalInitEdgeFlow = 1;
        pSrm->nMetalInitBondOrder = 1;
        pSrm->bStereoRemovesMetalFlag = pSrm->bFixStereoBonds;
        pSrm->nMetalFlowerParam_D = 16;
        pSrm->nMetalMaxCharge_D = 16;
    }
    else
    {
        pSrm->bMetalAddFlower = 0;
        pSrm->nMetalMinBondOrder = 1;
        pSrm->nMetalInitEdgeFlow = 0;
        pSrm->nMetalInitBondOrder = 1;
        pSrm->bStereoRemovesMetalFlag = pSrm->bFixStereoBonds;
        pSrm->nMetalFlowerParam_D = 16;
        pSrm->nMetalMaxCharge_D = 0;
    }
    pSrm->nMetal2EndpointInitBondOrder = pSrm->nMetal2EndpointMinBondOrder
        + pSrm->nMetal2EndpointInitEdgeFlow;
}
    */
    // END INCHI C FUNCTION: SetUpSrm
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: SetUpSrm
    // INCHI✔️❌: ichirvrs.h defines FIX_STEREO_BOND_ORDER=0.
    // INCHI✔️❌: ichirvrs.h defines METAL_FREE_CHARGE_VAL=1, selecting the first branch.
    // INCHI✔️❌: TARGET_API_LIB and COMPILE_ANSI_ONLY do not alter this function body.
    // END INCHI ACTIVE MACRO CONFIGURATION: SetUpSrm

    *pSrm = SRM::default();
    pSrm.bFixStereoBonds = 0;
    pSrm.nMetal2EndpointMinBondOrder = 1;
    pSrm.nMetal2EndpointInitEdgeFlow = 0;
    pSrm.bMetalAddFlower = 1;
    pSrm.nMetalMinBondOrder = 0;
    pSrm.nMetalInitEdgeFlow = 1;
    pSrm.nMetalInitBondOrder = 1;
    pSrm.bStereoRemovesMetalFlag = pSrm.bFixStereoBonds;
    pSrm.nMetalFlowerParam_D = 16;
    pSrm.nMetalMaxCharge_D = 16;
    pSrm.nMetal2EndpointInitBondOrder = pSrm
        .nMetal2EndpointMinBondOrder
        .wrapping_add(pSrm.nMetal2EndpointInitEdgeFlow);
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MergeStructureComponents(
    heap: &mut SourceHeap,
    _ip: &INPUT_PARMS,
    _sd: &mut STRUCT_DATA,
    _num_inp: i64,
    _szCurHdr: SourceMutPointer<i8>,
    _pSrm: &SRM,
    _bReqNonTaut: i32,
    pStruct: &[[SourceMutPointer<StrFromINChI>; TAUT_NUM as usize]; INCHI_NUM as usize],
    pOneInput: &mut InpInChI,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:406 MergeStructureComponents
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int MergeStructureComponents(ICHICONST INPUT_PARMS* ip,
    STRUCT_DATA* sd,
    long num_inp,
    char* szCurHdr,
    ICHICONST SRM* pSrm,
    int bReqNonTaut,
    StrFromINChI* pStruct[INCHI_NUM][TAUT_NUM],
    InpInChI* pOneInput)
{
    int iInchiRec, iMobileH, iAlternH, num_components, tot_just_atoms, tot_atoms, cur_nA, cur_nH; /* djb-rwth: removing redundant variables */
    int k, i, j, ret, iCurAtomOffs, iNxtAtomOffs, iCurDelHOffs, iNxtDelHOffs, len, len2, iShiftH, icomp;
    int* nAtomOffs = NULL, * nDelHOffs = NULL;
    StrFromINChI* pStruct1;
    inp_ATOM* at = NULL, * a;

    ret = 0;
    pOneInput->num_atoms = 0;
    /* select highest detail level */
    if ((num_components = pOneInput->nNumComponents[INCHI_REC][TAUT_NON])) /* djb-rwth: addressing LLVM warning */
    {
        iInchiRec = INCHI_REC;
        iMobileH = TAUT_NON;
    }
    else
    {
        if ((num_components = pOneInput->nNumComponents[INCHI_REC][TAUT_YES])) /* djb-rwth: addressing LLVM warning */
        {
            iInchiRec = INCHI_REC;
            iMobileH = TAUT_YES;
        }
        else
        {
            if ((num_components = pOneInput->nNumComponents[INCHI_BAS][TAUT_NON])) /* djb-rwth: addressing LLVM warning */
            {
                iInchiRec = INCHI_BAS;
                iMobileH = TAUT_NON;
            }
            else
            {
                if ((num_components = pOneInput->nNumComponents[INCHI_BAS][TAUT_YES])) /* djb-rwth: addressing LLVM warning */
                {
                    iInchiRec = INCHI_BAS;
                    iMobileH = TAUT_YES;
                }
                else
                {
                    return 0; /* no components available */
                }
            }
        }
    }

    nAtomOffs = (int*)inchi_malloc(((long long)num_components + 1) * sizeof(nAtomOffs[0])); /* djb-rwth: cast operator added */
    nDelHOffs = (int*)inchi_malloc(((long long)num_components + 1) * sizeof(nDelHOffs[0])); /* djb-rwth: cast operator added */
    if (!nAtomOffs || !nDelHOffs)
    {
        ret = RI_ERR_ALLOC;
        goto exit_function;
    }
    /* count number of atoms and removed H */
    /* djb-rwth: removing redundant code */
    iAlternH = (iMobileH == TAUT_NON && pOneInput->nNumComponents[iInchiRec][TAUT_YES]) ? TAUT_YES : -1;
    nAtomOffs[0] = nDelHOffs[0] = 0;
    for (k = 0; k < num_components; k++)
    {
        pStruct1 = pStruct[iInchiRec][iMobileH][k].num_atoms ? pStruct[iInchiRec][iMobileH] + k :
            iAlternH >= 0 &&
            pStruct[iInchiRec][iAlternH][k].num_atoms ? pStruct[iInchiRec][iAlternH] + k : NULL;
        if (!pStruct1 || !pStruct1->at2 || !pStruct1->num_atoms || pStruct1->bDeleted)
        {
            cur_nA = cur_nH = 0;
        }
        else
        {
            cur_nA = pStruct1->num_atoms;
            cur_nH = pStruct1->num_deleted_H;
        }
        nAtomOffs[k + 1] = nAtomOffs[k] + cur_nA;
        nDelHOffs[k + 1] = nDelHOffs[k] + cur_nH;
    }
    tot_just_atoms = nAtomOffs[num_components];
    /* shift all H to the end */
    for (k = 0; k <= num_components; k++)
    {
        nDelHOffs[k] += tot_just_atoms;
    }
    tot_atoms = nDelHOffs[num_components];

    /* merge atoms together: 1. Allocate */
    at = (inp_ATOM*)inchi_malloc(((long long)tot_atoms + 1) * sizeof(at[0])); /* djb-rwth: cast operator added */
    if (NULL == at)
    {
        ret = RI_ERR_ALLOC;
        goto exit_function;
    }
    if (!tot_atoms)
    {
        ret = 0;
        goto exit_function; /* empty structure */
    }
    /* merge atoms together: 2. Copy */
    for (k = 0; k < num_components; k++)
    {
        pStruct1 = pStruct[iInchiRec][iMobileH][k].num_atoms ? pStruct[iInchiRec][iMobileH] + k :
            iAlternH >= 0 &&
            pStruct[iInchiRec][iAlternH][k].num_atoms ? pStruct[iInchiRec][iAlternH] + k : NULL;
        if ((len = nAtomOffs[k + 1] - nAtomOffs[k]) && pStruct1) /* djb-rwth: addressing LLVM warning; fixing coverity ID #499555 */ 
        {
            memcpy(at + nAtomOffs[k], pStruct1->at2, len * sizeof(at[0]));
            if ((len2 = nDelHOffs[k + 1] - nDelHOffs[k])) /* djb-rwth: addressing LLVM warning */
            {
                memcpy(at + nDelHOffs[k], pStruct1->at2 + len, len2 * sizeof(at[0]));
            }
        }
    }
    /* merge atoms together: 3. Update atom numbers */
    icomp = 0;
    for (k = 0; k < num_components; k++)
    {
        iCurAtomOffs = nAtomOffs[k];
        iNxtAtomOffs = nAtomOffs[k + 1];
        iCurDelHOffs = nDelHOffs[k];
        iNxtDelHOffs = nDelHOffs[k + 1];
        len = nAtomOffs[k + 1] - nAtomOffs[k]; /* number of atoms in a component excluding explicit H */
        iShiftH = iCurDelHOffs - len;
        if (!len)
        {
            continue;
        }
        icomp++; /* current component number */
        /* update atoms */
        for (i = iCurAtomOffs; i < iNxtAtomOffs; i++)
        {

            a = at + i;

            a->endpoint = 0;
            a->bAmbiguousStereo = 0;
            a->at_type = 0;
            a->bCutVertex = 0;
            a->bUsed0DParity = 0;
            a->cFlags = 0;
            a->nBlockSystem = 0;
            a->nNumAtInRingSystem = 0;
            a->nRingSystem = 0;
            /* djb-rwth: addressing coverity ID #499524 -- initialisation with at */

            for (j = 0; j < a->valence; j++)
            {
                if (a->neighbor[j] < len)
                {
                    a->neighbor[j] += iCurAtomOffs; /* atom */
                }
                else
                {
                    a->neighbor[j] += iShiftH;      /* explicit H */
                }
            }
            a->orig_at_number += iCurAtomOffs;
            a->component = icomp;
            if (a->p_parity)
            {
                for (j = 0; j < MAX_NUM_STEREO_ATOM_NEIGH; j++)
                {
                    if (a->p_orig_at_num[j] <= len)
                    {
                        /* originally, orig_at_num = atom_index+1, therefore <= instead of < */
                        a->p_orig_at_num[j] += iCurAtomOffs;
                    }
                    else
                    {
                        a->p_orig_at_num[j] += iShiftH;
                    }
                }
            }
            for (j = 0; j < MAX_NUM_STEREO_BONDS && a->sb_parity[j]; j++)
            {
                if (a->sn_orig_at_num[j] <= len)
                {
                    /* originally, orig_at_num = atom_index+1, therefore <= instead of < */
                    a->sn_orig_at_num[j] += iCurAtomOffs;
                }
                else
                {
                    a->sn_orig_at_num[j] += iShiftH;
                }
            }
        }
        /* update fixed-H */
        for (i = iCurDelHOffs; i < iNxtDelHOffs; i++)
        {
            a = at + i;
            a->neighbor[0] += iCurAtomOffs;
            a->orig_at_number += iShiftH;
        }
    }
    /* save the results */
    pOneInput->atom = at;
    pOneInput->num_atoms = tot_atoms;
    at = NULL;

exit_function:
    if (at)        inchi_free(at);  /* in case of failure */
    if (nAtomOffs) inchi_free(nAtomOffs);
    if (nDelHOffs) inchi_free(nDelHOffs);
    return ret;
}
    */
    // END INCHI C FUNCTION: MergeStructureComponents
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: MergeStructureComponents
    // INCHI✔️❌: READ_INCHI_STRING=1 includes this production helper.
    // INCHI✔️❌: inchi_malloc/inchi_free resolve to the GCC/Linux libc allocation macros.
    // INCHI✔️❌: TARGET_API_LIB and COMPILE_ANSI_ONLY do not alter this function body.
    // END INCHI ACTIVE MACRO CONFIGURATION: MergeStructureComponents

    pOneInput.num_atoms = 0;
    let (representation, mobile_h, num_components) = if pOneInput.nNumComponents
        [INCHI_REC as usize][TAUT_NON as usize]
        != 0
    {
        (
            INCHI_REC as usize,
            TAUT_NON as usize,
            pOneInput.nNumComponents[INCHI_REC as usize][TAUT_NON as usize],
        )
    } else if pOneInput.nNumComponents[INCHI_REC as usize][TAUT_YES as usize] != 0 {
        (
            INCHI_REC as usize,
            TAUT_YES as usize,
            pOneInput.nNumComponents[INCHI_REC as usize][TAUT_YES as usize],
        )
    } else if pOneInput.nNumComponents[INCHI_BAS as usize][TAUT_NON as usize] != 0 {
        (
            INCHI_BAS as usize,
            TAUT_NON as usize,
            pOneInput.nNumComponents[INCHI_BAS as usize][TAUT_NON as usize],
        )
    } else if pOneInput.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] != 0 {
        (
            INCHI_BAS as usize,
            TAUT_YES as usize,
            pOneInput.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize],
        )
    } else {
        return Ok(0);
    };
    if num_components < 0 {
        return Err(SourceHeapError::SourceIntegerOverflow);
    }
    let count = usize::try_from(num_components)
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let allocation_count = count
        .checked_add(1)
        .ok_or(SourceHeapError::AllocationElementCountOutOfRange)?;
    let atom_offsets = match heap.allocate(vec![0_i32; allocation_count]) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
        Err(error) => return Err(error),
    };
    let hydrogen_offsets = match heap.allocate(vec![0_i32; allocation_count]) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
        Err(error) => {
            if !atom_offsets.is_null() {
                inchi_free(heap, atom_offsets)?;
            }
            return Err(error);
        }
    };
    if atom_offsets.is_null() || hydrogen_offsets.is_null() {
        if !atom_offsets.is_null() {
            inchi_free(heap, atom_offsets)?;
        }
        if !hydrogen_offsets.is_null() {
            inchi_free(heap, hydrogen_offsets)?;
        }
        return Ok(RI_ERR_ALLOC);
    }

    let alternate_h = if mobile_h == TAUT_NON as usize
        && pOneInput.nNumComponents[representation][TAUT_YES as usize] != 0
    {
        Some(TAUT_YES as usize)
    } else {
        None
    };
    let select_structure = |heap: &SourceHeap,
                            component: usize|
     -> Result<Option<StrFromINChI>, SourceHeapError> {
        let selected_base = pStruct[representation][mobile_h];
        if selected_base.is_null() {
            return Err(SourceHeapError::NullPointer);
        }
        let selected = heap
            .slice(selected_base.offset(component as i64)?.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        if selected.num_atoms != 0 {
            return Ok(Some(selected));
        }
        if let Some(alternate_h) = alternate_h {
            let alternate_base = pStruct[representation][alternate_h];
            if alternate_base.is_null() {
                return Err(SourceHeapError::NullPointer);
            }
            let alternate = heap
                .slice(alternate_base.offset(component as i64)?.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            if alternate.num_atoms != 0 {
                return Ok(Some(alternate));
            }
        }
        Ok(None)
    };

    let mut atoms = SourceMutPointer::null();
    let mut merged_atom_count = 0_i32;
    let execution = (|| -> Result<i32, SourceHeapError> {
        heap.slice_mut(atom_offsets)?[0] = 0;
        heap.slice_mut(hydrogen_offsets)?[0] = 0;
        for component in 0..count {
            let structure = select_structure(heap, component)?;
            let (atoms, hydrogens) = structure.as_ref().map_or((0, 0), |structure| {
                if structure.at2.is_null()
                    || structure.num_atoms == 0
                    || structure.bDeleted != 0
                {
                    (0, 0)
                } else {
                    (structure.num_atoms, structure.num_deleted_H)
                }
            });
            let previous_atoms = heap.slice(atom_offsets.as_const())?[component];
            let previous_hydrogens = heap.slice(hydrogen_offsets.as_const())?[component];
            heap.slice_mut(atom_offsets)?[component + 1] = previous_atoms.wrapping_add(atoms);
            heap.slice_mut(hydrogen_offsets)?[component + 1] =
                previous_hydrogens.wrapping_add(hydrogens);
        }
        let total_regular_atoms = heap.slice(atom_offsets.as_const())?[count];
        for component in 0..=count {
            let value = heap.slice(hydrogen_offsets.as_const())?[component];
            heap.slice_mut(hydrogen_offsets)?[component] =
                value.wrapping_add(total_regular_atoms);
        }
        merged_atom_count = heap.slice(hydrogen_offsets.as_const())?[count];
        if merged_atom_count < 0 {
            return Err(SourceHeapError::SourceIntegerOverflow);
        }
        let atom_count = usize::try_from(merged_atom_count)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
            .checked_add(1)
            .ok_or(SourceHeapError::AllocationElementCountOutOfRange)?;
        atoms = match heap.allocate(vec![inp_ATOM::default(); atom_count]) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => return Ok(RI_ERR_ALLOC),
            Err(error) => return Err(error),
        };
        if merged_atom_count == 0 {
            return Ok(0);
        }

        for component in 0..count {
            let structure = select_structure(heap, component)?;
            let regular_start = heap.slice(atom_offsets.as_const())?[component];
            let regular_end = heap.slice(atom_offsets.as_const())?[component + 1];
            let deleted_start = heap.slice(hydrogen_offsets.as_const())?[component];
            let deleted_end = heap.slice(hydrogen_offsets.as_const())?[component + 1];
            let regular_len = regular_end.wrapping_sub(regular_start);
            if regular_len != 0 {
                let structure = structure.ok_or(SourceHeapError::NullPointer)?;
                if regular_len < 0 || deleted_end.wrapping_sub(deleted_start) < 0 {
                    return Err(SourceHeapError::SourceIntegerOverflow);
                }
                for index in 0..regular_len {
                    let source = heap
                        .slice(structure.at2.offset(i64::from(index))?.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    let target = atoms.offset(i64::from(regular_start.wrapping_add(index)))?;
                    *heap
                        .slice_mut(target)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = source;
                }
                let deleted_len = deleted_end.wrapping_sub(deleted_start);
                for index in 0..deleted_len {
                    let source_index = regular_len.wrapping_add(index);
                    let source = heap
                        .slice(structure.at2.offset(i64::from(source_index))?.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    let target = atoms.offset(i64::from(deleted_start.wrapping_add(index)))?;
                    *heap
                        .slice_mut(target)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)? = source;
                }
            }
        }

        let mut current_component = 0_i32;
        for component in 0..count {
            let regular_start = heap.slice(atom_offsets.as_const())?[component];
            let regular_end = heap.slice(atom_offsets.as_const())?[component + 1];
            let deleted_start = heap.slice(hydrogen_offsets.as_const())?[component];
            let deleted_end = heap.slice(hydrogen_offsets.as_const())?[component + 1];
            let regular_len = regular_end.wrapping_sub(regular_start);
            let shift_h = deleted_start.wrapping_sub(regular_len);
            if regular_len == 0 {
                continue;
            }
            current_component = current_component.wrapping_add(1);
            for index in regular_start..regular_end {
                let pointer = atoms.offset(i64::from(index))?;
                let atom = heap
                    .slice_mut(pointer)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                atom.endpoint = 0;
                atom.bAmbiguousStereo = 0;
                atom.at_type = 0;
                atom.bCutVertex = 0;
                atom.bUsed0DParity = 0;
                atom.cFlags = 0;
                atom.nBlockSystem = 0;
                atom.nNumAtInRingSystem = 0;
                atom.nRingSystem = 0;
                if atom.valence > 0 {
                    for neighbor in 0..usize::try_from(atom.valence)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                    {
                        let value = atom
                            .neighbor
                            .get_mut(neighbor)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        *value = if i32::from(*value) < regular_len {
                            i32::from(*value).wrapping_add(regular_start) as AT_NUMB
                        } else {
                            i32::from(*value).wrapping_add(shift_h) as AT_NUMB
                        };
                    }
                }
                atom.orig_at_number =
                    i32::from(atom.orig_at_number).wrapping_add(regular_start) as AT_NUMB;
                atom.component = current_component as AT_NUMB;
                if atom.p_parity != 0 {
                    for stereo in 0..MAX_NUM_STEREO_ATOM_NEIGH as usize {
                        atom.p_orig_at_num[stereo] = if i32::from(atom.p_orig_at_num[stereo])
                            <= regular_len
                        {
                            i32::from(atom.p_orig_at_num[stereo]).wrapping_add(regular_start)
                                as AT_NUMB
                        } else {
                            i32::from(atom.p_orig_at_num[stereo]).wrapping_add(shift_h) as AT_NUMB
                        };
                    }
                }
                for stereo in 0..MAX_NUM_STEREO_BONDS as usize {
                    if atom.sb_parity[stereo] == 0 {
                        break;
                    }
                    atom.sn_orig_at_num[stereo] = if i32::from(atom.sn_orig_at_num[stereo])
                        <= regular_len
                    {
                        i32::from(atom.sn_orig_at_num[stereo]).wrapping_add(regular_start)
                            as AT_NUMB
                    } else {
                        i32::from(atom.sn_orig_at_num[stereo]).wrapping_add(shift_h) as AT_NUMB
                    };
                }
            }
            for index in deleted_start..deleted_end {
                let pointer = atoms.offset(i64::from(index))?;
                let atom = heap
                    .slice_mut(pointer)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                atom.neighbor[0] =
                    i32::from(atom.neighbor[0]).wrapping_add(regular_start) as AT_NUMB;
                atom.orig_at_number =
                    i32::from(atom.orig_at_number).wrapping_add(shift_h) as AT_NUMB;
            }
        }
        Ok(0)
    })();

    let transfer = matches!(execution, Ok(0)) && !atoms.is_null() && merged_atom_count != 0;
    if !atoms.is_null() && !transfer {
        inchi_free(heap, atoms)?;
    }
    let free_atoms = inchi_free(heap, atom_offsets);
    let free_hydrogens = inchi_free(heap, hydrogen_offsets);
    if let Err(error) = free_atoms {
        return Err(error);
    }
    if let Err(error) = free_hydrogens {
        return Err(error);
    }
    let result = execution?;
    if transfer {
        pOneInput.atom = atoms;
        pOneInput.num_atoms = merged_atom_count;
    }
    Ok(result)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn InChI2Atom(
    heap: &mut SourceHeap,
    ic: SourceMutPointer<INCHI_CLOCK>,
    pCG: &mut CANON_GLOBALS,
    ip: &INPUT_PARMS,
    sd: &mut STRUCT_DATA,
    szCurHdr: SourceConstPointer<i8>,
    num_inp: i64,
    pStruct: &mut StrFromINChI,
    iComponent: i32,
    iAtNoOffset: i32,
    bI2A_Flag: i32,
    bHasSomeFixedH: i32,
    one_input: &mut InpInChI,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:101 InChI2Atom
    // INCHI✔❌: complete source frame follows verbatim.
    /*
    int InChI2Atom(INCHI_CLOCK* ic,
        CANON_GLOBALS* pCG,
        ICHICONST INPUT_PARMS* ip,
        STRUCT_DATA* sd,
        const char* szCurHdr,
        long num_inp,
        StrFromINChI* pStruct,
        int iComponent,
        int iAtNoOffset,
        int  bI2A_Flag,
        int bHasSomeFixedH,
        InpInChI* OneInput)
    {
        int ret = 0;

        int iINChI = (bI2A_Flag & I2A_FLAG_RECMET) ? INCHI_REC : INCHI_BAS;
        int bMobileH = (bI2A_Flag & I2A_FLAG_FIXEDH) ? TAUT_NON : TAUT_YES;

        INChI* pInChI[TAUT_NUM];



        memset(pInChI, 0, sizeof(pInChI)); /* djb-rwth: memset_s C11/Annex K variant? */


        /* disconnected or reconnected */

        if (iINChI == INCHI_REC)
        {
            if (!OneInput->nNumComponents[iINChI][TAUT_YES])
            {
                iINChI = INCHI_BAS;
            }
        }
        if (iComponent >= OneInput->nNumComponents[iINChI][TAUT_YES])
        {
            return 0; /* component does not exist */
        }

        /* mobile or fixed H */
        pStruct->bFixedHExists = 0;
        if (bMobileH == TAUT_NON)
        {
            if (!OneInput->nNumComponents[iINChI][bMobileH])
            {
                /* only one InChI exists (no mobile H) */
                bMobileH = TAUT_YES;
            }
        }
        if (iComponent >= OneInput->nNumComponents[iINChI][bMobileH])
        {
            return 0; /* component does not exist */
        }

        /* pointer to the InChI that is going to be reversed */
        pInChI[0] = &OneInput->pInpInChI[iINChI][bMobileH][iComponent];
        pStruct->bMobileH = bMobileH;
        pStruct->iINCHI = iINChI;

        /* deleted component only in case Mobile-H and compound contains only protons */
        if (pInChI[0]->bDeleted)
        {
            return 0; /* deleted component, presumably H(+) */
        }

        if (bMobileH == TAUT_NON && OneInput->nNumProtons[iINChI][TAUT_YES].pNumProtons)
        {
            pStruct->nNumRemovedProtonsMobHInChI =
                OneInput->nNumProtons[iINChI][TAUT_YES].pNumProtons[iComponent].nNumRemovedProtons;
        }

        if (bMobileH == TAUT_NON ||
            (bMobileH == TAUT_YES &&
                OneInput->pInpInChI[iINChI][TAUT_NON] &&
                OneInput->pInpInChI[iINChI][TAUT_NON][iComponent].nNumberOfAtoms > 0 &&
                !OneInput->pInpInChI[iINChI][TAUT_NON][iComponent].bDeleted) /* djb-rwth: addressing LLVM warning */
            )
        {
            pStruct->bFixedHExists = 1;
        }

        if (bMobileH == TAUT_NON &&
            iComponent < OneInput->nNumComponents[iINChI][TAUT_YES] &&
            OneInput->pInpInChI[iINChI][TAUT_YES] &&
            OneInput->pInpInChI[iINChI][TAUT_YES][iComponent].nNumberOfAtoms > 0 &&
            !OneInput->pInpInChI[iINChI][TAUT_YES][iComponent].bDeleted
            )
        {
            /* pointer to the Mobile-H InChI if we are reversing Fixed-H InChI */
            pInChI[1] = &OneInput->pInpInChI[iINChI][TAUT_YES][iComponent];
        }

        pStruct->num_inp_actual = OneInput->num_inp;

        /* Intercept and correct non-polymer Zz to Zy if applicable */
        if (OneInput->polymer)
        {
            int a, k, new_num;
            OAD_Polymer* p = OneInput->polymer;
            pStruct->n_pzz = 0;
            pStruct->n_zy = 0;
            for (a = 0; a < pInChI[0]->nNumberOfAtoms; a++)
            {
                int aglob = iAtNoOffset + a + 1;
                if (pInChI[0]->nAtom[a] == EL_NUMBER_ZZ)
                {
                    new_num = EL_NUMBER_ZY; /* Zy */
                    for (k = 0; k < p->n; k++)
                    {
                        if ((aglob == p->units[k]->cap1) || (aglob == p->units[k]->cap2))
                        {
                            new_num = EL_NUMBER_ZZ;
                            break;
                        }
                    }
                    pInChI[0]->nAtom[a] = new_num;
                    if (new_num == EL_NUMBER_ZY)
                    {
                        pStruct->n_zy++;
                    }
                    else if (new_num == EL_NUMBER_ZZ)
                    {
                        pStruct->n_pzz++;
                    }
                }
            }
        }

        ret = OneInChI2Atom(ic, pCG, ip, sd, szCurHdr, num_inp, pStruct,
            iComponent, 0 /* iAtNoOffset*/, bHasSomeFixedH, pInChI);

        /* djb-rwth: fixing oss-fuzz issue #66758, #30283 */
        if (pStruct->at && pStruct->at2)
        {
            int a;
            for (a = 0; a < pInChI[0]->nNumberOfAtoms; a++)
            {
                if (pInChI[0]->nAtom[a] == EL_NUMBER_ZY)
                {
                    pInChI[0]->nAtom[a] = EL_NUMBER_ZZ;
                    pStruct->at[a].el_number = EL_NUMBER_ZZ;
                    strcpy(pStruct->at[a].elname, "Zz");
                    pStruct->at2[a].el_number = EL_NUMBER_ZZ;
                    strcpy(pStruct->at2[a].elname, "Zz");
                    pStruct->n_zy--;
                    pStruct->n_pzz++;
                }
            }
        }

        return ret; /* same interpretation as in ProcessOneStructure ??? */
    }
    */
    // END INCHI C FUNCTION: InChI2Atom
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: InChI2Atom
    // INCHI✔❌: READ_INCHI_STRING=1 includes this production function.
    // INCHI✔❌: COMPILE_ANSI_ONLY and TARGET_API_LIB do not alter this function body.
    // END INCHI ACTIVE MACRO CONFIGURATION: InChI2Atom

    let mut i_inchi = if bI2A_Flag & I2A_FLAG_RECMET as i32 != 0 {
        INCHI_REC as usize
    } else {
        INCHI_BAS as usize
    };
    let mut mobile_h = if bI2A_Flag & I2A_FLAG_FIXEDH as i32 != 0 {
        TAUT_NON as usize
    } else {
        TAUT_YES as usize
    };
    let mut selected = [SourceMutPointer::<INChI>::null(); 2];

    if i_inchi == INCHI_REC as usize && one_input.nNumComponents[i_inchi][TAUT_YES as usize] == 0 {
        i_inchi = INCHI_BAS as usize;
    }
    if iComponent >= one_input.nNumComponents[i_inchi][TAUT_YES as usize] {
        return Ok(0);
    }

    pStruct.bFixedHExists = 0;
    if mobile_h == TAUT_NON as usize && one_input.nNumComponents[i_inchi][mobile_h] == 0 {
        mobile_h = TAUT_YES as usize;
    }
    if iComponent >= one_input.nNumComponents[i_inchi][mobile_h] {
        return Ok(0);
    }

    selected[0] = one_input.pInpInChI[i_inchi][mobile_h].offset(i64::from(iComponent))?;
    pStruct.bMobileH = mobile_h as i8;
    pStruct.iINCHI = i_inchi as i8;
    let selected_value = heap
        .slice(selected[0].as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    if selected_value.bDeleted != 0 {
        return Ok(0);
    }

    if mobile_h == TAUT_NON as usize {
        let proton_components = one_input.nNumProtons[i_inchi][TAUT_YES as usize].pNumProtons;
        if !proton_components.is_null() {
            let proton = heap
                .slice(proton_components.offset(i64::from(iComponent))?.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            pStruct.nNumRemovedProtonsMobHInChI = i32::from(proton.nNumRemovedProtons);
        }
    }

    if mobile_h == TAUT_NON as usize {
        pStruct.bFixedHExists = 1;
    } else {
        let fixed_components = one_input.pInpInChI[i_inchi][TAUT_NON as usize];
        if !fixed_components.is_null() {
            let fixed = heap
                .slice(fixed_components.offset(i64::from(iComponent))?.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if fixed.nNumberOfAtoms > 0 && fixed.bDeleted == 0 {
                pStruct.bFixedHExists = 1;
            }
        }
    }

    if mobile_h == TAUT_NON as usize
        && iComponent < one_input.nNumComponents[i_inchi][TAUT_YES as usize]
    {
        let mobile_components = one_input.pInpInChI[i_inchi][TAUT_YES as usize];
        if !mobile_components.is_null() {
            let component = mobile_components.offset(i64::from(iComponent))?;
            let mobile = heap
                .slice(component.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if mobile.nNumberOfAtoms > 0 && mobile.bDeleted == 0 {
                selected[1] = component;
            }
        }
    }

    pStruct.num_inp_actual = one_input.num_inp;

    if !one_input.polymer.is_null() {
        pStruct.n_pzz = 0;
        pStruct.n_zy = 0;
        let polymer = heap
            .slice(one_input.polymer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        for atom_index in 0..selected_value.nNumberOfAtoms {
            let atom_pointer = selected_value.nAtom.offset(i64::from(atom_index))?;
            let atom_number = *heap
                .slice(atom_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let global_atom = iAtNoOffset
                .checked_add(atom_index)
                .and_then(|value| value.checked_add(1))
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            if atom_number == EL_NUMBER_ZZ {
                let mut new_number = EL_NUMBER_ZY;
                for unit_index in 0..polymer.n {
                    let unit_pointer = *heap
                        .slice(polymer.units.offset(i64::from(unit_index))?.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let unit = heap
                        .slice(unit_pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if global_atom == unit.cap1 || global_atom == unit.cap2 {
                        new_number = EL_NUMBER_ZZ;
                        break;
                    }
                }
                *heap
                    .slice_mut(atom_pointer)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = new_number;
                if new_number == EL_NUMBER_ZY {
                    pStruct.n_zy = pStruct
                        .n_zy
                        .checked_add(1)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                } else if new_number == EL_NUMBER_ZZ {
                    pStruct.n_pzz = pStruct
                        .n_pzz
                        .checked_add(1)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                }
            }
        }
    }

    let result = OneInChI2Atom(
        heap,
        ic,
        pCG,
        ip,
        sd,
        szCurHdr,
        num_inp,
        pStruct,
        iComponent,
        0,
        bHasSomeFixedH,
        selected,
        clock_result,
    );
    if !pStruct.at.is_null() && !pStruct.at2.is_null() {
        for atom_index in 0..selected_value.nNumberOfAtoms {
            let atom_pointer = selected_value.nAtom.offset(i64::from(atom_index))?;
            let atom_number = *heap
                .slice(atom_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if atom_number == EL_NUMBER_ZY {
                *heap
                    .slice_mut(atom_pointer)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = EL_NUMBER_ZZ;
                for output in [pStruct.at, pStruct.at2] {
                    let atom = heap
                        .slice_mut(output.offset(i64::from(atom_index))?)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    atom.el_number = EL_NUMBER_ZZ;
                    atom.elname[0] = b'Z' as i8;
                    atom.elname[1] = b'z' as i8;
                    atom.elname[2] = 0;
                }
                pStruct.n_zy = pStruct
                    .n_zy
                    .checked_sub(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                pStruct.n_pzz = pStruct
                    .n_pzz
                    .checked_add(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            }
        }
    }

    result
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn AllInchiToStructure(
    heap: &mut SourceHeap,
    ic: SourceMutPointer<INCHI_CLOCK>,
    pCG: &mut CANON_GLOBALS,
    ip_inp: &INPUT_PARMS,
    sd_inp: &STRUCT_DATA,
    num_inp: i64,
    szCurHdr: SourceConstPointer<i8>,
    pSrm: SourceConstPointer<SRM>,
    bHasSomeFixedH: i32,
    pStruct: &mut [[SourceMutPointer<StrFromINChI>; 2]; 2],
    pOneInput: &mut InpInChI,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1042 AllInchiToStructure
    // INCHI✔❌: complete source frame follows verbatim.
    /*
    int AllInchiToStructure(INCHI_CLOCK* ic,
        CANON_GLOBALS* pCG,
        ICHICONST INPUT_PARMS* ip_inp,
        STRUCT_DATA* sd_inp,
        long num_inp,
        char* szCurHdr,
        ICHICONST SRM* pSrm,
        int bHasSomeFixedH,
        StrFromINChI* pStruct[INCHI_NUM][TAUT_NUM],
        InpInChI* pOneInput)
    {
        int iInchiRec, iMobileH, cur_num_comp, bCurI2A_Flag, k, ret, num_err;
        int iAtNoOffset;
        INPUT_PARMS* ip, ip_loc;
        STRUCT_DATA* sd, sd_loc;
        long          ulProcessingTime = 0;
        inchiTime     ulTStart;

        InchiTimeGet(&ulTStart);
        ip = &ip_loc;
        *ip = *ip_inp;
        sd = &sd_loc;
        memset(sd, 0, sizeof(*sd)); /* djb-rwth: memset_s C11/Annex K variant? */
        sd->ulStructTime = sd_inp->ulStructTime;
        ret = 0;
        num_err = 0;
        for (iInchiRec = 0; iInchiRec < INCHI_NUM; iInchiRec++)
        {
            /* Disconnected/Connected */
            for (iMobileH = 0; iMobileH < TAUT_NUM; iMobileH++)
            {
                /* Mobile/Fixed H */
                cur_num_comp = pOneInput->nNumComponents[iInchiRec][iMobileH];
                if (!cur_num_comp)
                {
                    continue;
                }
                /* allocate memory for all existing components */
                pStruct[iInchiRec][iMobileH] = (StrFromINChI*)inchi_calloc(cur_num_comp, sizeof(pStruct[0][0][0]));
                if (!pStruct[iInchiRec][iMobileH])
                {
                    ret = RI_ERR_ALLOC;
                    goto exit_error;
                }
                /* set conversion mode */
                bCurI2A_Flag = (iMobileH ? 0 : I2A_FLAG_FIXEDH) | (iInchiRec ? I2A_FLAG_RECMET : 0);
                if (iMobileH)
                {
                    ip->nMode &= ~REQ_MODE_BASIC;
                }
                else
                {
                    ip->nMode |= REQ_MODE_BASIC;
                }
                /* InChI --> structure conversion for all components except duplicated */
                iAtNoOffset = 0;
                for (k = 0; k < cur_num_comp; k++)
                {
                    /* components */
                    if ((!iMobileH && !pOneInput->pInpInChI[iInchiRec][iMobileH][k].nNumberOfAtoms) ||
                        pOneInput->pInpInChI[iInchiRec][iMobileH][k].bDeleted ||
                        pOneInput->pInpInChI[iInchiRec][iMobileH][k].nLink < 0) /* djb-rwth: addressing LLVM warning */
                    {
                        pStruct[iInchiRec][iMobileH][k].nLink = pOneInput->pInpInChI[iInchiRec][iMobileH][k].nLink;
                        pStruct[iInchiRec][iMobileH][k].bDeleted = pOneInput->pInpInChI[iInchiRec][iMobileH][k].bDeleted;
                        continue; /* do not create a structure out of an unavailable
                                     Fixed-H InChI or out of the one present in Reconnected layer */

    #ifdef NEVER  /* a wrong attempt to process deleted components here */
                        if (pStruct[iInchiRec][iMobileH][k].nLink = pOneInput->pInpInChI[iInchiRec][iMobileH][k].nLink)
                        {
                            continue; /* do not create a structure out of an unavailable
                                         Fixed-H InChI or out of the one present in Reconnected layer */
                        }
                        else
                            if (iMobileH && pOneInput->pInpInChI[iInchiRec][iMobileH][k].nNumberOfAtoms &&
                                pOneInput->pInpInChI[iInchiRec][iMobileH][k].bDeleted &&
                                pOneInput->pInpInChI[iInchiRec][iMobileH][0].bDeleted)
                            {
                                /* all components are protons */
                                ;
                            }
                            else
                            {
                                continue;
                            }
    #endif

                    }
                    if (bHasSomeFixedH && iMobileH && k < pOneInput->nNumComponents[iInchiRec][TAUT_NON] &&
                        pOneInput->pInpInChI[iInchiRec][TAUT_NON][k].nNumberOfAtoms)
                    {
                        continue; /* do not process Mobile-H if Fixed-H is requested and exists */
                    }
                    pStruct[iInchiRec][iMobileH][k].pSrm = pSrm;
                    pStruct[iInchiRec][iMobileH][k].iInchiRec = iInchiRec;
                    pStruct[iInchiRec][iMobileH][k].iMobileH = iMobileH;

                    /****************************************************/
                    /*                                                  */
                    /* Convert InChI of one component into a Structure  */
                    /*                                                  */
                    /****************************************************/

                    ret = InChI2Atom(ic, pCG, ip, sd, szCurHdr, num_inp, pStruct[iInchiRec][iMobileH] + k, k,
                        iAtNoOffset /* 0*/, bCurI2A_Flag, bHasSomeFixedH, pOneInput);
                    pStruct[iInchiRec][iMobileH][k].nLink = pOneInput->pInpInChI[iInchiRec][iMobileH][k].nLink;
                    if (ret < 0)
                    {
    #if ( bRELEASE_VERSION != 1 )
    #ifndef TARGET_API_LIB
                        /* !!! Conversion Error -- Ignore for now !!! */
                        fprintf(stdout, "%ld %s Conversion failed: %d, %c%c comp %d\n",
                            num_inp, szCurHdr ? szCurHdr : "Struct", ret, iInchiRec ? 'R' : 'D', iMobileH ? 'M' : 'F', k + 1);
    #endif
    #endif
                        if (ret == CT_USER_QUIT_ERR)
                        {
                            goto exit_error;
                        }
                        pStruct[iInchiRec][iMobileH][k].nError = ret;
                        ret = 0; /* force to ignore the errors for now !!!! */
                        num_err++;
                    }
                    iAtNoOffset += pOneInput->pInpInChI[iInchiRec][iMobileH][k].nNumberOfAtoms;
                } /* k-th component */
            }
        }

    exit_error:
        ulProcessingTime += InchiTimeElapsed(ic, &ulTStart);
        sd->ulStructTime += ulProcessingTime;

        return ret < 0 ? ret : num_err;
    }
    */
    // END INCHI C FUNCTION: AllInchiToStructure
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: AllInchiToStructure
    // INCHI✔❌: READ_INCHI_STRING=1 includes this production function.
    // INCHI✔❌: NEVER is undefined, excluding the wrong deleted-component attempt.
    // INCHI✔❌: TARGET_API_LIB is defined, excluding the diagnostic fprintf branch.
    // INCHI✔❌: COMPILE_ANSI_ONLY does not alter this function body.
    // END INCHI ACTIVE MACRO CONFIGURATION: AllInchiToStructure

    let mut start = inchiTime::default();
    InchiTimeGet(&mut start, clock_result);
    let mut ip = ip_inp.clone();
    let mut sd = STRUCT_DATA {
        ulStructTime: sd_inp.ulStructTime,
        ..STRUCT_DATA::default()
    };
    let mut ret = 0_i32;
    let mut num_err = 0_i32;

    'all_layers: for i_inchi_rec in 0..2_usize {
        for i_mobile_h in 0..2_usize {
            let component_count = pOneInput.nNumComponents[i_inchi_rec][i_mobile_h];
            if component_count == 0 {
                continue;
            }
            let structures = match inchi_calloc::<StrFromINChI>(
                heap,
                component_count as u64,
                std::mem::size_of::<StrFromINChI>() as u64,
            ) {
                Ok(pointer) => pointer,
                Err(_) => {
                    ret = RI_ERR_ALLOC;
                    break 'all_layers;
                }
            };
            pStruct[i_inchi_rec][i_mobile_h] = structures;

            let current_flags = (if i_mobile_h != 0 {
                0
            } else {
                I2A_FLAG_FIXEDH as i32
            }) | if i_inchi_rec != 0 {
                I2A_FLAG_RECMET as i32
            } else {
                0
            };
            if i_mobile_h != 0 {
                ip.nMode &= !u64::from(REQ_MODE_BASIC);
            } else {
                ip.nMode |= u64::from(REQ_MODE_BASIC);
            }

            let mut atom_offset = 0_i32;
            for component_index in 0..component_count {
                let input_pointer = pOneInput.pInpInChI[i_inchi_rec][i_mobile_h]
                    .offset(i64::from(component_index))?;
                let input_component = heap
                    .slice(input_pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone();
                let structure_pointer = structures.offset(i64::from(component_index))?;

                if (i_mobile_h == TAUT_NON as usize && input_component.nNumberOfAtoms == 0)
                    || input_component.bDeleted != 0
                    || input_component.nLink < 0
                {
                    let structure = heap
                        .slice_mut(structure_pointer)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    structure.nLink = input_component.nLink;
                    structure.bDeleted = input_component.bDeleted as i8;
                    continue;
                }

                if bHasSomeFixedH != 0
                    && i_mobile_h != 0
                    && component_index < pOneInput.nNumComponents[i_inchi_rec][TAUT_NON as usize]
                {
                    let fixed_pointer = pOneInput.pInpInChI[i_inchi_rec][TAUT_NON as usize]
                        .offset(i64::from(component_index))?;
                    let fixed_component = heap
                        .slice(fixed_pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if fixed_component.nNumberOfAtoms != 0 {
                        continue;
                    }
                }

                // The validated in-place borrow is Rust's representation of
                // the source `pStruct[...] + k` pointer. It remains live for
                // exactly the source mutation interval.
                let conversion_result =
                    heap.with_slice_mut_and_heap_mut(structure_pointer, |structures, heap| {
                        let structure = structures
                            .first_mut()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        // INCHI✔️✔️: pStruct[iInchiRec][iMobileH][k].pSrm = pSrm;
                        structure.pSrm = pSrm;
                        // INCHI✔️✔️: pStruct[iInchiRec][iMobileH][k].iInchiRec = iInchiRec;
                        structure.iInchiRec = i_inchi_rec as i8;
                        // INCHI✔️✔️: pStruct[iInchiRec][iMobileH][k].iMobileH = iMobileH;
                        structure.iMobileH = i_mobile_h as i8;
                        // INCHI✔️✔️: ret = InChI2Atom(ic, pCG, ip, sd, szCurHdr, num_inp,
                        // INCHI✔️✔️:     pStruct[iInchiRec][iMobileH] + k, k,
                        // INCHI✔️✔️:     iAtNoOffset, bCurI2A_Flag, bHasSomeFixedH, pOneInput);
                        let conversion = InChI2Atom(
                            heap,
                            ic,
                            pCG,
                            &ip,
                            &mut sd,
                            szCurHdr,
                            num_inp,
                            structure,
                            component_index,
                            atom_offset,
                            current_flags,
                            bHasSomeFixedH,
                            pOneInput,
                            clock_result,
                        )?;
                        // INCHI✔️✔️: pStruct[iInchiRec][iMobileH][k].nLink =
                        // INCHI✔️✔️:     pOneInput->pInpInChI[iInchiRec][iMobileH][k].nLink;
                        structure.nLink = input_component.nLink;
                        Ok(conversion)
                    });
                ret = conversion_result?;
                if ret < 0 {
                    if ret == CT_USER_QUIT_ERR {
                        break 'all_layers;
                    }
                    heap.slice_mut(structure_pointer)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nError = ret;
                    ret = 0;
                    num_err = num_err
                        .checked_add(1)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                }
                atom_offset = atom_offset
                    .checked_add(input_component.nNumberOfAtoms)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            }
        }
    }

    let mut clock = heap
        .slice(ic.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    let processing_time = InchiTimeElapsed(&mut clock, Some(&start), clock_result);
    *heap
        .slice_mut(ic)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = clock;
    sd.ulStructTime = sd.ulStructTime.wrapping_add(processing_time as u64);

    Ok(if ret < 0 { ret } else { num_err })
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn AddProtonAndIsoHBalanceToMobHStruct(
    heap: &mut SourceHeap,
    ic: SourceMutPointer<INCHI_CLOCK>,
    pCG: SourceMutPointer<CANON_GLOBALS>,
    ip: &INPUT_PARMS,
    sd: &mut STRUCT_DATA,
    num_inp: i64,
    bHasSomeFixedH: i32,
    szCurHdr: SourceMutPointer<i8>,
    pStruct: &[[SourceMutPointer<StrFromINChI>; TAUT_NUM as usize]; INCHI_NUM as usize],
    pOneInput: &InpInChI,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1180 AddProtonAndIsoHBalanceToMobHStruct
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int AddProtonAndIsoHBalanceToMobHStruct(INCHI_CLOCK* ic,
    CANON_GLOBALS* pCG,
    ICHICONST INPUT_PARMS* ip,
    STRUCT_DATA* sd,
    long num_inp,
    int bHasSomeFixedH,
    char* szCurHdr,
    StrFromINChI* pStruct[INCHI_NUM][TAUT_NUM],
    InpInChI* pOneInput)
{
    COMPONENT_REM_PROTONS nToBeRemovedByNormFromRevrs[INCHI_NUM];
    int                   nRemovedByNormFromRevrs[INCHI_NUM];
    int                   nRemovedByRevrs[INCHI_NUM];

    int   nDeltaFromDisconnected = 0, nRemovedProtonsByNormFromRevrs, nRemovedProtonsByRevrs; /* djb-rwth: removing redundant variables */
    NUM_H nIsoDeltaFromDisconnected[NUM_H_ISOTOPES];
    int iInchiRec, i, k, k1, ret = 0;
    int  nChargeInChI, nChargeRevrs;

    if (bHasSomeFixedH)
    {
        return 0; /* 2005-03-01 */
    }

    /* num protons removed by InChI Normalization from the original structure */
    for (i = 0; i < INCHI_NUM; i++)
    {
        nToBeRemovedByNormFromRevrs[i].nNumRemovedProtons = pOneInput->nNumProtons[i][TAUT_YES].nNumRemovedProtons;
        for (k = 0; k < NUM_H_ISOTOPES; k++)
        {
            nToBeRemovedByNormFromRevrs[i].nNumRemovedIsotopicH[k] = pOneInput->nNumProtons[i][TAUT_YES].nNumRemovedIsotopicH[k];
        }
    }
    /* accumulate here num. protons removed by the normalization from the reversed structure */
    nRemovedByNormFromRevrs[INCHI_BAS] =
        nRemovedByNormFromRevrs[INCHI_REC] = 0;
    nRemovedByRevrs[INCHI_REC] =
        nRemovedByRevrs[INCHI_BAS] = 0;
    /* protons added/removed by InChI Normalization to/from Restored Structure might have been added by StructureRestore */

    for (iInchiRec = 0; iInchiRec < INCHI_NUM; iInchiRec++)
    {
        for (k = 0; k < pOneInput->nNumComponents[iInchiRec][TAUT_YES]; k++)
        {
            if (!bInpInchiComponentExists(pOneInput, iInchiRec, TAUT_YES, k))
            {
                continue;
            }
            nRemovedProtonsByNormFromRevrs = 0; /* Num protons removed from the Restored Structure by InChI Normalization */
            nRemovedProtonsByRevrs = 0; /* Num protons removed by the Reconstruction from the Restored Structure */
            if (iInchiRec == INCHI_REC || (iInchiRec == INCHI_BAS && (pStruct[iInchiRec][TAUT_YES][k].nLink) >= 0)) /* djb-rwth: addressing LLVM warning; removing redundant code */
            {

                REV_INCHI* pRevInChI = &pStruct[iInchiRec][TAUT_YES][k].RevInChI;
                INChI_Aux** pINChI_Aux2 = pRevInChI->pINChI_Aux[iInchiRec][0]; /* component 0*/
                INChI** pINChI_Revr = pRevInChI->pINChI[iInchiRec][0];
                INChI* pINChI_Orig = pOneInput->pInpInChI[iInchiRec][TAUT_YES] + k;
                nChargeRevrs = pINChI_Revr ? pINChI_Revr[TAUT_YES]->nTotalCharge : NO_VALUE_INT;
                nChargeInChI = pINChI_Orig->nTotalCharge;
                if (pINChI_Aux2)
                {
                    nRemovedProtonsByNormFromRevrs = pINChI_Aux2[TAUT_YES]->nNumRemovedProtons;
                }
                nRemovedProtonsByRevrs = pStruct[iInchiRec][TAUT_YES][k].nNumRemovedProtonsByRevrs;
                pStruct[iInchiRec][TAUT_YES][k].nChargeRevrs = nChargeRevrs;
                pStruct[iInchiRec][TAUT_YES][k].nChargeInChI = nChargeInChI;
            }
            else
            {
                if (0 <= (k1 = -(1 + pStruct[iInchiRec][TAUT_YES][k].nLink)))
                {
                    REV_INCHI* pRevInChI = &pStruct[INCHI_REC][TAUT_YES][k1].RevInChI;
                    INChI_Aux** pINChI_Aux2 = pRevInChI->pINChI_Aux[INCHI_BAS][0]; /* component 0 */
                    INChI** pINChI_Revr = pRevInChI->pINChI[INCHI_BAS][0];
                    INChI* pINChI_Orig = pOneInput->pInpInChI[INCHI_REC][TAUT_YES] + k1;
                    nChargeRevrs = pINChI_Revr ? pINChI_Revr[TAUT_YES]->nTotalCharge : NO_VALUE_INT;
                    nChargeInChI = pINChI_Orig->nTotalCharge;
                    if (pINChI_Aux2)
                    {
                        nRemovedProtonsByNormFromRevrs = pINChI_Aux2[TAUT_YES]->nNumRemovedProtons;
                    }
                    /* this component cannot be disconnected because it is same as in reconnected layer */
                    nRemovedProtonsByRevrs = pStruct[INCHI_REC][TAUT_YES][k1].nNumRemovedProtonsByRevrs;
                    pStruct[iInchiRec][TAUT_YES][k1].nChargeRevrs = nChargeRevrs;
                    pStruct[iInchiRec][TAUT_YES][k1].nChargeInChI = nChargeInChI;
                }
            }
            /* how many protons (to be removed by InChI Normalization) to add =
            (proton balance in InChI} -
            {number of protons known to be removed by InChI Normalization from Reconstructed structure} */
            nToBeRemovedByNormFromRevrs[iInchiRec].nNumRemovedProtons -= nRemovedProtonsByNormFromRevrs;
            nRemovedByNormFromRevrs[iInchiRec] += nRemovedProtonsByNormFromRevrs;
            nRemovedByRevrs[iInchiRec] += nRemovedProtonsByRevrs;
            pStruct[iInchiRec][TAUT_YES][k].nRemovedProtonsByNormFromRevrs = nRemovedProtonsByNormFromRevrs;
        }
    }

    /* Since fixed-H layer is missing we need to add proton balance to the components */
    memset(nIsoDeltaFromDisconnected, 0, sizeof(nIsoDeltaFromDisconnected)); /* djb-rwth: memset_s C11/Annex K variant? */
    for (iInchiRec = INCHI_REC; INCHI_BAS <= iInchiRec; iInchiRec--)
    {
        /*
        if ( !pOneInput->nNumComponents[iInchiRec][TAUT_NON] &&
              pOneInput->nNumComponents[iInchiRec][TAUT_YES] ) {
        */
        int bHasRecMobH = (iInchiRec == INCHI_BAS && pOneInput->nNumComponents[INCHI_REC][TAUT_YES]);
        /* bHasRecMobH means all components that could not be disconnected are in reconnected part */
        if (iInchiRec == INCHI_BAS)
        {
            /* second pass: common structures have been changed */
            nToBeRemovedByNormFromRevrs[INCHI_BAS].nNumRemovedProtons += nDeltaFromDisconnected;
        }
        /* after proton removal InChI is recalculated */

        ret = AddRemProtonsInRestrStruct(ic, pCG, ip, sd, num_inp, bHasSomeFixedH, pStruct[iInchiRec][TAUT_YES],
            pOneInput->nNumComponents[iInchiRec][TAUT_YES],
            bHasRecMobH ? pStruct[INCHI_REC][TAUT_YES] : NULL,
            bHasRecMobH ? pOneInput->nNumComponents[INCHI_REC][TAUT_YES] : 0,
            &nToBeRemovedByNormFromRevrs[iInchiRec].nNumRemovedProtons,
            (iInchiRec == INCHI_REC) ? &nDeltaFromDisconnected : NULL);
        if (ret < 0)
        {
            goto exit_function;
        }
        /* djb-rwth: removing redundant code */
       /*
       }
       */
    }

    /* if fixed-H layer is missing then we need to add isotopic exchangeable proton balance to the components */
    for (iInchiRec = INCHI_REC; INCHI_BAS <= iInchiRec; iInchiRec--)
    {
        /*
        if ( !pOneInput->nNumComponents[iInchiRec][TAUT_NON] &&
              pOneInput->nNumComponents[iInchiRec][TAUT_YES] ) {
        */
        int bHasRecMobH = (iInchiRec == INCHI_BAS && pOneInput->nNumComponents[INCHI_REC][TAUT_YES]);
        /* bHasRecMobH means all components that could not be disconnected are in reconnected part */
        if (iInchiRec == INCHI_BAS)
        {
            /* second pass: common structures have been changed */
            for (k = 0; k < NUM_H_ISOTOPES; k++)
            {
                nToBeRemovedByNormFromRevrs[INCHI_BAS].nNumRemovedIsotopicH[k] += nIsoDeltaFromDisconnected[k];
            }
        }

        /* after proton removal InChI is recalculated */
        ret = AddRemIsoProtonsInRestrStruct(ic, pCG, ip, sd, num_inp, bHasSomeFixedH, pStruct[iInchiRec][TAUT_YES],
            pOneInput->nNumComponents[iInchiRec][TAUT_YES],
            bHasRecMobH ? pStruct[INCHI_REC][TAUT_YES] : NULL,
            bHasRecMobH ? pOneInput->nNumComponents[INCHI_REC][TAUT_YES] : 0,
            nToBeRemovedByNormFromRevrs[iInchiRec].nNumRemovedIsotopicH,
            (iInchiRec == INCHI_REC) ? nIsoDeltaFromDisconnected : NULL);

        if (ret < 0)
        {
            goto exit_function;
        }
        /* djb-rwth: removing redundant code */
        /*
        }
        */
    }

exit_function:

    return ret;
}
    */
    // END INCHI C FUNCTION: AddProtonAndIsoHBalanceToMobHStruct
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: AddProtonAndIsoHBalanceToMobHStruct
    // INCHI✔️❌: READ_INCHI_STRING=1 includes this function; TARGET_API_LIB removes only inactive diagnostic output outside this body.
    // INCHI✔️❌: NUM_H is signed 16-bit; INCHI_NUM=2, TAUT_NUM=2, NUM_H_ISOTOPES=3, NO_VALUE_INT=9999.
    // INCHI✔️❌: SourceHeap checked pointer ownership adds overhead versus direct C pointer arithmetic.
    // END INCHI ACTIVE MACRO CONFIGURATION: AddProtonAndIsoHBalanceToMobHStruct

    let _ = szCurHdr;
    if bHasSomeFixedH != 0 {
        return Ok(0);
    }

    let mut nToBeRemovedByNormFromRevrs: [COMPONENT_REM_PROTONS; INCHI_NUM as usize] =
        std::array::from_fn(|representation| COMPONENT_REM_PROTONS {
            nNumRemovedProtons: pOneInput.nNumProtons[representation][TAUT_YES as usize]
                .nNumRemovedProtons,
            nNumRemovedIsotopicH: pOneInput.nNumProtons[representation][TAUT_YES as usize]
                .nNumRemovedIsotopicH,
        });
    let mut nRemovedByNormFromRevrs = [0_i32; INCHI_NUM as usize];
    let mut nRemovedByRevrs = [0_i32; INCHI_NUM as usize];

    for iInchiRec in 0..INCHI_NUM as usize {
        let component_count = pOneInput.nNumComponents[iInchiRec][TAUT_YES as usize];
        let mut k = 0_i32;
        while k < component_count {
            if bInpInchiComponentExists(
                heap,
                pOneInput,
                iInchiRec as i32,
                TAUT_YES as i32,
                k,
            )? == 0
            {
                k = k.wrapping_add(1);
                continue;
            }

            let structure_pointer = pStruct[iInchiRec][TAUT_YES as usize].offset(i64::from(k))?;
            let structure = heap
                .slice(structure_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let mut nRemovedProtonsByNormFromRevrs = 0_i32;
            let mut nRemovedProtonsByRevrs = 0_i32;
            if iInchiRec == INCHI_REC as usize
                || (iInchiRec == INCHI_BAS as usize && structure.nLink >= 0)
            {
                let reversed = &structure.RevInChI;
                let reversed_rows = reversed.pINChI[iInchiRec];
                let nChargeRevrs = if reversed_rows.is_null() {
                    NO_VALUE_INT as i32
                } else {
                    let row = heap
                        .slice(reversed_rows.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    heap.slice(row[TAUT_YES as usize].as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nTotalCharge
                };
                let original_pointer = pOneInput.pInpInChI[iInchiRec][TAUT_YES as usize]
                    .offset(i64::from(k))?;
                let nChargeInChI = heap
                    .slice(original_pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nTotalCharge;
                let aux_rows = reversed.pINChI_Aux[iInchiRec];
                if !aux_rows.is_null() {
                    let row = heap
                        .slice(aux_rows.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    nRemovedProtonsByNormFromRevrs = i32::from(
                        heap.slice(row[TAUT_YES as usize].as_const())?
                            .first()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nNumRemovedProtons,
                    );
                }
                nRemovedProtonsByRevrs = structure.nNumRemovedProtonsByRevrs;
                let target = heap
                    .slice_mut(structure_pointer)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                target.nChargeRevrs = nChargeRevrs;
                target.nChargeInChI = nChargeInChI;
            } else {
                let k1 = structure.nLink.wrapping_add(1).wrapping_neg();
                if k1 >= 0 {
                    let rec_pointer = pStruct[INCHI_REC as usize][TAUT_YES as usize]
                        .offset(i64::from(k1))?;
                    let rec_structure = heap
                        .slice(rec_pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    let reversed = &rec_structure.RevInChI;
                    let reversed_rows = reversed.pINChI[INCHI_BAS as usize];
                    let nChargeRevrs = if reversed_rows.is_null() {
                        NO_VALUE_INT as i32
                    } else {
                        let row = heap
                            .slice(reversed_rows.as_const())?
                            .first()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        heap.slice(row[TAUT_YES as usize].as_const())?
                            .first()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nTotalCharge
                    };
                    let original_pointer = pOneInput.pInpInChI[INCHI_REC as usize]
                        [TAUT_YES as usize]
                        .offset(i64::from(k1))?;
                    let nChargeInChI = heap
                        .slice(original_pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nTotalCharge;
                    let aux_rows = reversed.pINChI_Aux[INCHI_BAS as usize];
                    if !aux_rows.is_null() {
                        let row = heap
                            .slice(aux_rows.as_const())?
                            .first()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        nRemovedProtonsByNormFromRevrs = i32::from(
                            heap.slice(row[TAUT_YES as usize].as_const())?
                                .first()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .nNumRemovedProtons,
                        );
                    }
                    nRemovedProtonsByRevrs = rec_structure.nNumRemovedProtonsByRevrs;
                    let target_pointer = pStruct[iInchiRec][TAUT_YES as usize]
                        .offset(i64::from(k1))?;
                    let target = heap
                        .slice_mut(target_pointer)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    target.nChargeRevrs = nChargeRevrs;
                    target.nChargeInChI = nChargeInChI;
                }
            }

            let proton_balance = &mut nToBeRemovedByNormFromRevrs[iInchiRec].nNumRemovedProtons;
            *proton_balance = i32::from(*proton_balance)
                .wrapping_sub(nRemovedProtonsByNormFromRevrs) as _;
            nRemovedByNormFromRevrs[iInchiRec] = nRemovedByNormFromRevrs[iInchiRec]
                .wrapping_add(nRemovedProtonsByNormFromRevrs);
            nRemovedByRevrs[iInchiRec] =
                nRemovedByRevrs[iInchiRec].wrapping_add(nRemovedProtonsByRevrs);
            heap.slice_mut(structure_pointer)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .nRemovedProtonsByNormFromRevrs = nRemovedProtonsByNormFromRevrs;
            k = k.wrapping_add(1);
        }
    }

    let mut nDeltaFromDisconnected = 0_i32;
    let mut nIsoDeltaFromDisconnected = [0; NUM_H_ISOTOPES as usize];
    let mut ret = 0_i32;
    for iInchiRec in (INCHI_BAS as usize..=INCHI_REC as usize).rev() {
        let bHasRecMobH = iInchiRec == INCHI_BAS as usize
            && pOneInput.nNumComponents[INCHI_REC as usize][TAUT_YES as usize] != 0;
        if iInchiRec == INCHI_BAS as usize {
            let balance = &mut nToBeRemovedByNormFromRevrs[INCHI_BAS as usize]
                .nNumRemovedProtons;
            *balance = i32::from(*balance).wrapping_add(nDeltaFromDisconnected) as _;
        }
        if nToBeRemovedByNormFromRevrs[iInchiRec].nNumRemovedProtons == 0 {
            ret = 0;
            continue;
        }
        let rec_snapshot = if bHasRecMobH {
            let count = usize::try_from(
                pOneInput.nNumComponents[INCHI_REC as usize][TAUT_YES as usize],
            )
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            Some(
                heap.slice(pStruct[INCHI_REC as usize][TAUT_YES as usize].as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec(),
            )
        } else {
            None
        };
        let structure_pointer = pStruct[iInchiRec][TAUT_YES as usize];
        let component_count = pOneInput.nNumComponents[iInchiRec][TAUT_YES as usize];
        let balance = &mut nToBeRemovedByNormFromRevrs[iInchiRec].nNumRemovedProtons;
        ret = if structure_pointer.is_null() {
            AddRemProtonsInRestrStruct(
                heap,
                ic,
                pCG,
                ip,
                sd,
                num_inp,
                bHasSomeFixedH,
                &mut [],
                component_count,
                rec_snapshot.as_deref(),
                if bHasRecMobH {
                    pOneInput.nNumComponents[INCHI_REC as usize][TAUT_YES as usize]
                } else {
                    0
                },
                balance,
                if iInchiRec == INCHI_REC as usize {
                    Some(&mut nDeltaFromDisconnected)
                } else {
                    None
                },
                clock_result,
            )?
        } else {
            heap.with_slice_mut_and_heap_mut(structure_pointer, |structures, heap| {
                AddRemProtonsInRestrStruct(
                    heap,
                    ic,
                    pCG,
                    ip,
                    sd,
                    num_inp,
                    bHasSomeFixedH,
                    structures,
                    component_count,
                    rec_snapshot.as_deref(),
                    if bHasRecMobH {
                        pOneInput.nNumComponents[INCHI_REC as usize][TAUT_YES as usize]
                    } else {
                        0
                    },
                    balance,
                    if iInchiRec == INCHI_REC as usize {
                        Some(&mut nDeltaFromDisconnected)
                    } else {
                        None
                    },
                    clock_result,
                )
            })?
        };
        if ret < 0 {
            return Ok(ret);
        }
    }

    for iInchiRec in (INCHI_BAS as usize..=INCHI_REC as usize).rev() {
        let bHasRecMobH = iInchiRec == INCHI_BAS as usize
            && pOneInput.nNumComponents[INCHI_REC as usize][TAUT_YES as usize] != 0;
        if iInchiRec == INCHI_BAS as usize {
            for isotope in 0..NUM_H_ISOTOPES as usize {
                let balance = &mut nToBeRemovedByNormFromRevrs[INCHI_BAS as usize]
                    .nNumRemovedIsotopicH[isotope];
                *balance = i32::from(*balance)
                    .wrapping_add(i32::from(nIsoDeltaFromDisconnected[isotope])) as _;
            }
        }
        let mut isotopic_balance_present = 0_i32;
        for balance in nToBeRemovedByNormFromRevrs[iInchiRec].nNumRemovedIsotopicH {
            isotopic_balance_present |= i32::from(balance);
        }
        if isotopic_balance_present == 0 {
            ret = 0;
            continue;
        }
        let rec_snapshot = if bHasRecMobH {
            let count = usize::try_from(
                pOneInput.nNumComponents[INCHI_REC as usize][TAUT_YES as usize],
            )
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            Some(
                heap.slice(pStruct[INCHI_REC as usize][TAUT_YES as usize].as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec(),
            )
        } else {
            None
        };
        let structure_pointer = pStruct[iInchiRec][TAUT_YES as usize];
        let component_count = pOneInput.nNumComponents[iInchiRec][TAUT_YES as usize];
        let balance = &mut nToBeRemovedByNormFromRevrs[iInchiRec].nNumRemovedIsotopicH;
        ret = if structure_pointer.is_null() {
            AddRemIsoProtonsInRestrStruct(
                heap,
                ic,
                pCG,
                ip,
                sd,
                num_inp,
                bHasSomeFixedH,
                &mut [],
                component_count,
                rec_snapshot.as_deref(),
                if bHasRecMobH {
                    pOneInput.nNumComponents[INCHI_REC as usize][TAUT_YES as usize]
                } else {
                    0
                },
                balance,
                if iInchiRec == INCHI_REC as usize {
                    Some(&mut nIsoDeltaFromDisconnected)
                } else {
                    None
                },
                clock_result,
            )?
        } else {
            heap.with_slice_mut_and_heap_mut(structure_pointer, |structures, heap| {
                AddRemIsoProtonsInRestrStruct(
                    heap,
                    ic,
                    pCG,
                    ip,
                    sd,
                    num_inp,
                    bHasSomeFixedH,
                    structures,
                    component_count,
                    rec_snapshot.as_deref(),
                    if bHasRecMobH {
                        pOneInput.nNumComponents[INCHI_REC as usize][TAUT_YES as usize]
                    } else {
                        0
                    },
                    balance,
                    if iInchiRec == INCHI_REC as usize {
                        Some(&mut nIsoDeltaFromDisconnected)
                    } else {
                        None
                    },
                    clock_result,
                )
            })?
        };
        if ret < 0 {
            return Ok(ret);
        }
    }

    Ok(ret)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn FreeStrFromINChI(
    heap: &mut SourceHeap,
    pStruct: &mut [[SourceMutPointer<StrFromINChI>; TAUT_NUM as usize]; INCHI_NUM as usize],
    nNumComponents: &[[i32; TAUT_NUM as usize]; INCHI_NUM as usize],
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1353 FreeStrFromINChI
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
void FreeStrFromINChI(StrFromINChI* pStruct[INCHI_NUM][TAUT_NUM],
    int nNumComponents[INCHI_NUM][TAUT_NUM])
{
    int iInchiRec, iMobileH, cur_num_comp, k, j;
    StrFromINChI* pStruct1;
    for (iInchiRec = 0; iInchiRec < INCHI_NUM; iInchiRec++)
    {
        for (iMobileH = 0; iMobileH < TAUT_NUM; iMobileH++)
        {
            cur_num_comp = nNumComponents[iInchiRec][iMobileH];
            if (!cur_num_comp || !(pStruct1 = pStruct[iInchiRec][iMobileH]))
            {
                continue;
            }
            for (k = 0; k < cur_num_comp; k++)
            {
                if (pStruct1[k].at)
                {
                    inchi_free(pStruct1[k].at);
                }
                if (pStruct1[k].at2)
                {
                    inchi_free(pStruct1[k].at2);
                }
                if (pStruct1[k].st)
                {
                    inchi_free(pStruct1[k].st);
                }
                if (pStruct1[k].pVA)
                {
                    inchi_free(pStruct1[k].pVA);
                }
                /*
                if ( pStruct1[k].ti.t_group ) {
                    inchi_free( pStruct1[k].ti.t_group );
                }
                */
                if (pStruct1[k].One_ti.t_group) {
                    inchi_free(pStruct1[k].One_ti.t_group); /* ricrogz: fixing memory leak */
                }
                if (pStruct1[k].pXYZ)
                {
                    inchi_free(pStruct1[k].pXYZ); /* djb-rwth: unresolved issue -- revision required? -- false positive as this function just does the clean-up job */
                }
                /*==== begin ====*/
                free_t_group_info(&pStruct1[k].ti);
                if (pStruct1[k].endpoint)
                {
                    inchi_free(pStruct1[k].endpoint); /* djb-rwth: unresolved issue -- revision required? -- false positive as this function just does the clean-up job */
                }
                if (pStruct1[k].fixed_H)
                {
                    inchi_free(pStruct1[k].fixed_H); /* djb-rwth: unresolved issue -- revision required? -- false positive as this function just does the clean-up job */
                }
                for (j = 0; j < TAUT_NUM; j++)
                {
                    if (pStruct1[k].nAtno2Canon[j])
                        inchi_free(pStruct1[k].nAtno2Canon[j]);
                    if (pStruct1[k].nCanon2Atno[j])
                        inchi_free(pStruct1[k].nCanon2Atno[j]);
                }
                /*===== end ======*/
                /*  free INChI memory */
                FreeAllINChIArrays(pStruct1[k].RevInChI.pINChI,
                    pStruct1[k].RevInChI.pINChI_Aux,
                    pStruct1[k].RevInChI.num_components);
#ifdef NEVER
                /* don't do that: these are just pointers to OneInput structure members */
                Free_INChI(&pStruct1[k].pINChI);
                Free_INChI_Aux(&pStruct1[k].pINChI_Aux);
                if (pStruct1[k].inp_norm_data)
                {
                    FreeInpAtomData(pStruct1[k].inp_norm_data);
                    inchi_free(pStruct1[k].inp_norm_data);
                }
#endif
            }
            inchi_free(pStruct[iInchiRec][iMobileH]);
            pStruct[iInchiRec][iMobileH] = NULL;
        }
    }
}
    */
    // END INCHI C FUNCTION: FreeStrFromINChI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FreeStrFromINChI
    // INCHI✔️❌: READ_INCHI_STRING=1 includes this function; NEVER is undefined, excluding borrowed-member cleanup.
    // INCHI✔️❌: COMPILE_ANSI_ONLY and TARGET_API_LIB do not alter the active body.
    // INCHI✔️❌: SourceHeap checked frees add overhead versus direct C free calls.
    // END INCHI ACTIVE MACRO CONFIGURATION: FreeStrFromINChI

    for iInchiRec in 0..INCHI_NUM as usize {
        for iMobileH in 0..TAUT_NUM as usize {
            let cur_num_comp = nNumComponents[iInchiRec][iMobileH];
            let structures = pStruct[iInchiRec][iMobileH];
            if cur_num_comp == 0 || structures.is_null() {
                continue;
            }
            let components = if cur_num_comp > 0 {
                let count = usize::try_from(cur_num_comp)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                heap.slice(structures.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec()
            } else {
                Vec::new()
            };
            for mut structure in components {
                for pointer in [structure.at, structure.at2] {
                    if !pointer.is_null() {
                        inchi_free(heap, pointer)?;
                    }
                }
                if !structure.st.is_null() {
                    inchi_free(heap, structure.st)?;
                }
                if !structure.pVA.is_null() {
                    inchi_free(heap, structure.pVA)?;
                }
                if !structure.One_ti.t_group.is_null() {
                    inchi_free(heap, structure.One_ti.t_group)?;
                }
                if !structure.pXYZ.is_null() {
                    inchi_free(heap, structure.pXYZ)?;
                }
                free_t_group_info(heap, Some(&mut structure.ti))?;
                if !structure.endpoint.is_null() {
                    inchi_free(heap, structure.endpoint)?;
                }
                if !structure.fixed_H.is_null() {
                    inchi_free(heap, structure.fixed_H)?;
                }
                for j in 0..TAUT_NUM as usize {
                    if !structure.nAtno2Canon[j].is_null() {
                        inchi_free(heap, structure.nAtno2Canon[j])?;
                    }
                    if !structure.nCanon2Atno[j].is_null() {
                        inchi_free(heap, structure.nCanon2Atno[j])?;
                    }
                }
                FreeAllINChIArrays(
                    heap,
                    &mut structure.RevInChI.pINChI,
                    &mut structure.RevInChI.pINChI_Aux,
                    &mut structure.RevInChI.num_components,
                )?;
            }
            inchi_free(heap, structures)?;
            pStruct[iInchiRec][iMobileH] = SourceMutPointer::null();
        }
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn FreeInpInChI(
    heap: &mut SourceHeap,
    p_one_input: &mut InpInChI,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1438 FreeInpInChI
    // INCHI✔❌: void FreeInpInChI(InpInChI* pOneInput)
    // INCHI✔❌: {
    // INCHI✔❌:     int iINChI, k, j;
    // INCHI✔❌:     for (iINChI = 0; iINChI < INCHI_NUM; iINChI++)
    // INCHI✔❌:     {
    // INCHI✔❌:         for (j = 0; j < TAUT_NUM; j++)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (pOneInput->pInpInChI[iINChI][j])
    // INCHI✔❌:             {
    // INCHI✔❌:                 for (k = 0; k < pOneInput->nNumComponents[iINChI][j]; k++)
    // INCHI✔❌:                 {
    // INCHI✔❌: #if (FIX_OSS_FUZZ_25734_28139 == 1)
    // INCHI✔❌:                     U_CHAR* k_nAtom = (&pOneInput->pInpInChI[iINChI][j][k])->nAtom;
    // INCHI✔❌:                     AT_NUMB* k_nConnTable = (&pOneInput->pInpInChI[iINChI][j][k])->nConnTable;
    // INCHI✔❌:                     AT_NUMB* k_nTautomer = (&pOneInput->pInpInChI[iINChI][j][k])->nTautomer;
    // INCHI✔❌:                     S_CHAR* k_nNum_H = (&pOneInput->pInpInChI[iINChI][j][k])->nNum_H;
    // INCHI✔❌:                     S_CHAR* k_nNum_H_fixed = (&pOneInput->pInpInChI[iINChI][j][k])->nNum_H_fixed;
    // INCHI✔❌:                     char* k_szHillFormula = (&pOneInput->pInpInChI[iINChI][j][k])->szHillFormula;
    // INCHI✔❌:                     AT_NUMB* k_nPossibleLocationsOfIsotopicH = (&pOneInput->pInpInChI[iINChI][j][k])->nPossibleLocationsOfIsotopicH;
    // INCHI✔❌:                     INChI_IsotopicAtom* k_IsotopicAtom = (&pOneInput->pInpInChI[iINChI][j][k])->IsotopicAtom;
    // INCHI✔❌:                     INChI_IsotopicTGroup* k_IsotopicTGroup = (&pOneInput->pInpInChI[iINChI][j][k])->IsotopicTGroup;
    // INCHI✔❌:                     INChI_Stereo* k_Stereo = (&pOneInput->pInpInChI[iINChI][j][k])->Stereo;
    // INCHI✔❌:                     INChI_Stereo* k_StereoIsotopic = (&pOneInput->pInpInChI[iINChI][j][k])->StereoIsotopic;
    // INCHI✔❌:
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌:                     Free_INChI_Members(&pOneInput->pInpInChI[iINChI][j][k]);
    // INCHI✔❌:
    // INCHI✔❌: #if (FIX_OSS_FUZZ_25734_28139 == 1)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         /* prevent erroneous repeated freeing in copied pInpInChIp[][][kk] */
    // INCHI✔❌:                         int kk;
    // INCHI✔❌:                         for (kk = k + 1; kk < pOneInput->nNumComponents[iINChI][j]; kk++)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             if (k_nAtom == (&pOneInput->pInpInChI[iINChI][j][kk])->nAtom)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->nAtom = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             if (k_nConnTable == (&pOneInput->pInpInChI[iINChI][j][kk])->nConnTable)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->nConnTable = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             if (k_nTautomer == (&pOneInput->pInpInChI[iINChI][j][kk])->nTautomer)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->nTautomer = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             if (k_nNum_H == (&pOneInput->pInpInChI[iINChI][j][kk])->nNum_H)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->nNum_H = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             if (k_nNum_H_fixed == (&pOneInput->pInpInChI[iINChI][j][kk])->nNum_H_fixed)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->nNum_H_fixed = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             if (k_szHillFormula == (&pOneInput->pInpInChI[iINChI][j][kk])->szHillFormula)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->szHillFormula = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             if (k_nPossibleLocationsOfIsotopicH == (&pOneInput->pInpInChI[iINChI][j][kk])->nPossibleLocationsOfIsotopicH)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->nPossibleLocationsOfIsotopicH = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:
    // INCHI✔❌:                             if (k_IsotopicAtom == (&pOneInput->pInpInChI[iINChI][j][kk])->IsotopicAtom)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->IsotopicAtom = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:
    // INCHI✔❌:                             if (k_IsotopicTGroup == (&pOneInput->pInpInChI[iINChI][j][kk])->IsotopicTGroup)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->IsotopicTGroup = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:
    // INCHI✔❌:                             if (k_Stereo == (&pOneInput->pInpInChI[iINChI][j][kk])->Stereo)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->Stereo = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:
    // INCHI✔❌:                             if (k_StereoIsotopic == (&pOneInput->pInpInChI[iINChI][j][kk])->StereoIsotopic)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 (&pOneInput->pInpInChI[iINChI][j][kk])->StereoIsotopic = NULL;
    // INCHI✔❌:                             }
    // INCHI✔❌:
    // INCHI✔❌:                         }
    // INCHI✔❌:                     }
    // INCHI✔❌: #endif
    // INCHI✔❌:
    // INCHI✔❌:                 }
    // INCHI✔❌:                 inchi_free(pOneInput->pInpInChI[iINChI][j]);
    // INCHI✔❌:                 pOneInput->pInpInChI[iINChI][j] = NULL;
    // INCHI✔❌:             }
    // INCHI✔❌:             if (pOneInput->nNumProtons[iINChI][j].pNumProtons)
    // INCHI✔❌:             {
    // INCHI✔❌:                 inchi_free(pOneInput->nNumProtons[iINChI][j].pNumProtons);
    // INCHI✔❌:                 pOneInput->nNumProtons[iINChI][j].pNumProtons = NULL;
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:     if (pOneInput->atom)
    // INCHI✔❌:     {
    // INCHI✔❌:         inchi_free(pOneInput->atom);
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     FreeExtOrigAtData(pOneInput->polymer, pOneInput->v3000);
    // INCHI✔❌:
    // INCHI✔❌:     memset(pOneInput, 0, sizeof(*pOneInput)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌: }
    // END INCHI C FUNCTION: FreeInpInChI

    for i_inchi in 0..2_usize {
        for j in 0..2_usize {
            let components = p_one_input.pInpInChI[i_inchi][j];
            let component_count = p_one_input.nNumComponents[i_inchi][j];
            if !components.is_null() {
                for k in 0..component_count {
                    let k = i64::from(k);
                    let component = components.offset(k)?;
                    let snapshot = heap
                        .slice(component.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    Free_INChI_Members(heap, component)?;
                    for kk in (k + 1)..i64::from(component_count) {
                        let later = components.offset(kk)?;
                        let value = heap
                            .slice_mut(later)?
                            .first_mut()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if snapshot.nAtom == value.nAtom {
                            value.nAtom = SourceMutPointer::null();
                        }
                        if snapshot.nConnTable == value.nConnTable {
                            value.nConnTable = SourceMutPointer::null();
                        }
                        if snapshot.nTautomer == value.nTautomer {
                            value.nTautomer = SourceMutPointer::null();
                        }
                        if snapshot.nNum_H == value.nNum_H {
                            value.nNum_H = SourceMutPointer::null();
                        }
                        if snapshot.nNum_H_fixed == value.nNum_H_fixed {
                            value.nNum_H_fixed = SourceMutPointer::null();
                        }
                        if snapshot.szHillFormula == value.szHillFormula {
                            value.szHillFormula = SourceMutPointer::null();
                        }
                        if snapshot.nPossibleLocationsOfIsotopicH
                            == value.nPossibleLocationsOfIsotopicH
                        {
                            value.nPossibleLocationsOfIsotopicH = SourceMutPointer::null();
                        }
                        if snapshot.IsotopicAtom == value.IsotopicAtom {
                            value.IsotopicAtom = SourceMutPointer::null();
                        }
                        if snapshot.IsotopicTGroup == value.IsotopicTGroup {
                            value.IsotopicTGroup = SourceMutPointer::null();
                        }
                        if snapshot.Stereo == value.Stereo {
                            value.Stereo = SourceMutPointer::null();
                        }
                        if snapshot.StereoIsotopic == value.StereoIsotopic {
                            value.StereoIsotopic = SourceMutPointer::null();
                        }
                    }
                }
                inchi_free(heap, components)?;
                p_one_input.pInpInChI[i_inchi][j] = SourceMutPointer::null();
            }
            let protons = p_one_input.nNumProtons[i_inchi][j].pNumProtons;
            if !protons.is_null() {
                inchi_free(heap, protons)?;
                p_one_input.nNumProtons[i_inchi][j].pNumProtons = SourceMutPointer::null();
            }
        }
    }
    if !p_one_input.atom.is_null() {
        inchi_free(heap, p_one_input.atom)?;
    }
    FreeExtOrigAtData(heap, p_one_input.polymer, p_one_input.v3000)?;
    *p_one_input = InpInChI::default();
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn CompareAllOrigInchiToRevInChI(
    heap: &mut SourceHeap,
    pStruct: &[[SourceMutPointer<StrFromINChI>; TAUT_NUM as usize]; INCHI_NUM as usize],
    pOneInput: &mut InpInChI,
    bReqNonTaut: i32,
    num_inp: i64,
    szCurHdr: SourceMutPointer<i8>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1548 CompareAllOrigInchiToRevInChI
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int CompareAllOrigInchiToRevInChI(StrFromINChI* pStruct[INCHI_NUM][TAUT_NUM],
    InpInChI* pOneInput,
    int bReqNonTaut,
    long num_inp,
    char* szCurHdr)
{
    int i, iInchiRec, iMobileH, iMobileHpStruct, num_components, iComponent, ret = 0; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    COMPONENT_REM_PROTONS nCurRemovedProtons, nNumRemovedProtons;
    INChI* pInChI[TAUT_NUM];
    INCHI_MODE  CompareInchiFlags[TAUT_NUM];
    memset(pOneInput->CompareInchiFlags[0], 0, sizeof(pOneInput->CompareInchiFlags[0])); /* djb-rwth: memset_s C11/Annex K variant? */
    memset(&nNumRemovedProtons, 0, sizeof(nNumRemovedProtons)); /* djb-rwth: memset_s C11/Annex K variant? */

    /* do we have reconnected InChI ?*/
    iInchiRec = INCHI_REC;
    iMobileH = TAUT_NON;
    if (!pOneInput->nNumComponents[iInchiRec][TAUT_YES] && !pOneInput->nNumComponents[iInchiRec][TAUT_NON])
    {
        iInchiRec = INCHI_BAS;
    }
    /* do we have Mobile or Fixed-H ? */
    if (!pOneInput->nNumComponents[iInchiRec][TAUT_NON] || !bReqNonTaut)
    {
        iMobileH = TAUT_YES;  /* index for pOneInput */
    }
    /* if a restored structure has Fixed-H InChI then its mobile-H restored InChI is in Fixed-H pStruct */
    num_components = pOneInput->nNumComponents[iInchiRec][iMobileH];
    for (iComponent = 0; iComponent < num_components; iComponent++)
    {
        int bMobileH = iMobileH;
        pInChI[0] = pInChI[1] = NULL;
        if (pOneInput->pInpInChI[iInchiRec][bMobileH][iComponent].nNumberOfAtoms &&
            !pOneInput->pInpInChI[iInchiRec][bMobileH][iComponent].bDeleted)
        {
            /* the requested InChI layer exists */
            pInChI[0] = &pOneInput->pInpInChI[iInchiRec][bMobileH][iComponent];
            if (bMobileH == TAUT_NON)
            {
                pInChI[1] = &pOneInput->pInpInChI[iInchiRec][TAUT_YES][iComponent];
            }
        }
        else
        {
            if (bMobileH == TAUT_NON &&
                pOneInput->pInpInChI[iInchiRec][TAUT_YES][iComponent].nNumberOfAtoms &&
                !pOneInput->pInpInChI[iInchiRec][TAUT_YES][iComponent].bDeleted)
            {
                /* the requested Fixed-H InChI layer does not exist; however, the Mobile-H does exist */
                bMobileH = TAUT_YES; /* only Mobile-H is available */
                pInChI[0] = &pOneInput->pInpInChI[iInchiRec][bMobileH][iComponent];
            }
        }
        memset(CompareInchiFlags, 0, sizeof(CompareInchiFlags)); /* djb-rwth: memset_s C11/Annex K variant? */
        memset(&nCurRemovedProtons, 0, sizeof(nCurRemovedProtons)); /* djb-rwth: memset_s C11/Annex K variant? */
        iMobileHpStruct =
#if ( bRELEASE_VERSION == 0 )
#ifndef TARGET_API_LIB
            /* legacy: reproduce old output */
            OldPrintCompareOneOrigInchiToRevInChI(pStruct[iInchiRec][bMobileH] + iComponent, pInChI, bMobileH,
                iComponent, num_inp, szCurHdr);
#endif
#endif
        /* one component comparison result bits */
        ret = CompareOneOrigInchiToRevInChI(pStruct[iInchiRec][bMobileH] + iComponent, pInChI, bMobileH, iComponent,
            num_inp, szCurHdr, &nCurRemovedProtons, CompareInchiFlags);
        if (ret >= 0)
        {
            /* no errors encountered -> accumulate removed protons from individual Mobile-H layers of components */
            nNumRemovedProtons.nNumRemovedProtons += nCurRemovedProtons.nNumRemovedProtons;
            for (i = 0; i < NUM_H_ISOTOPES; i++)
            {
                nNumRemovedProtons.nNumRemovedIsotopicH[i] += nCurRemovedProtons.nNumRemovedIsotopicH[i];
            }
            /* accumulate compare bits */
            for (i = 0; i < TAUT_NUM; i++)
            {
                pOneInput->CompareInchiFlags[0][i] |= CompareInchiFlags[i];
            }
        }
        else
        {
            goto exit_function;
        }
    }
    if (iMobileH == TAUT_YES)
    {
        if (pOneInput->nNumProtons[iInchiRec][iMobileH].pNumProtons)
        {
            ret = RI_ERR_PROGR; /* in Mobile-H case proton balances are split between compoments */
        }
        else
        {
            /*   num removed protons in orig. InChI      num removed protons in restored InChi */
            if (nNumRemovedProtons.nNumRemovedProtons != pOneInput->nNumProtons[iInchiRec][iMobileH].nNumRemovedProtons)
            {
                /* restored structure InChI has less or more removed protons */
                pOneInput->CompareInchiFlags[0][TAUT_YES] |= INCHIDIFF_MOBH_PROTONS;
#if ( bRELEASE_VERSION == 0 )
                /* debug output only */
                {
                    int num_H_AddedByRevrs = pOneInput->nNumProtons[iInchiRec][iMobileH].nNumRemovedProtons
                        - nNumRemovedProtons.nNumRemovedProtons;
                    fprintf(stdout, "COMPARE_INCHI: %ld: %s %cM: Proton balance (Diff: %d, RevrsRem=%d)\n",
                        num_inp, szCurHdr ? szCurHdr : "Struct", iInchiRec ? 'R' : 'D',
                        pOneInput->nNumProtons[iInchiRec][iMobileH].nNumRemovedProtons, num_H_AddedByRevrs);
                }
#endif
            }
            for (i = 0; i < NUM_H_ISOTOPES; i++)
            {
                if (nNumRemovedProtons.nNumRemovedIsotopicH[i] != pOneInput->nNumProtons[iInchiRec][TAUT_YES].nNumRemovedIsotopicH[i])
                {
                    pOneInput->CompareInchiFlags[0][TAUT_YES] |= INCHIDIFF_MOB_ISO_H;
#if ( bRELEASE_VERSION == 0 )
                    /* debug output only */
                    {
                        int num_H_AddedByRevrs = pOneInput->nNumProtons[iInchiRec][TAUT_YES].nNumRemovedIsotopicH[i]
                            - nNumRemovedProtons.nNumRemovedIsotopicH[i];
                        fprintf(stdout, "COMPARE_INCHI: %ld: %s %cM: Iso Xchg %dH balance (Diff: %d, RevrsRem=%d)\n",
                            num_inp, szCurHdr ? szCurHdr : "Struct", iInchiRec ? 'R' : 'D', i + 1,
                            pOneInput->nNumProtons[iInchiRec][TAUT_YES].nNumRemovedIsotopicH[i], num_H_AddedByRevrs);
                    }
#endif
                }
            }
        }
    }

exit_function:

    return ret;
}
    */
    // END INCHI C FUNCTION: CompareAllOrigInchiToRevInChI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CompareAllOrigInchiToRevInChI
    // INCHI✔️❌: bRELEASE_VERSION=1 and TARGET_API_LIB exclude OldPrint/debug output.
    // INCHI✔️❌: COMPILE_ANSI_ONLY does not alter the active production body.
    // INCHI✔️❌: Component clones avoid holding SourceHeap borrows across mutable callee calls.
    // END INCHI ACTIVE MACRO CONFIGURATION: CompareAllOrigInchiToRevInChI

    pOneInput.CompareInchiFlags[0] = [0; TAUT_NUM as usize];
    let mut total_removed = COMPONENT_REM_PROTONS::default();
    let mut inchi_representation = INCHI_REC as usize;
    let mut mobile_h = TAUT_NON as usize;
    if pOneInput.nNumComponents[inchi_representation][TAUT_YES as usize] == 0
        && pOneInput.nNumComponents[inchi_representation][TAUT_NON as usize] == 0
    {
        inchi_representation = INCHI_BAS as usize;
    }
    if pOneInput.nNumComponents[inchi_representation][TAUT_NON as usize] == 0
        || bReqNonTaut == 0
    {
        mobile_h = TAUT_YES as usize;
    }
    let num_components = pOneInput.nNumComponents[inchi_representation][mobile_h];
    let mut ret = 0_i32;
    if num_components > 0 {
        let count = usize::try_from(num_components)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        for component in 0..count {
            let mut selected_mobile = mobile_h;
            let mut input = [SourceMutPointer::null(); TAUT_NUM as usize];
            let requested_base = pOneInput.pInpInChI[inchi_representation][selected_mobile];
            if requested_base.is_null() {
                return Err(SourceHeapError::NullPointer);
            }
            let requested_pointer = requested_base.offset(component as i64)?;
            let requested = heap
                .slice(requested_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if requested.nNumberOfAtoms != 0 && requested.bDeleted == 0 {
                input[0] = requested_pointer;
                if selected_mobile == TAUT_NON as usize {
                    let mobile_base =
                        pOneInput.pInpInChI[inchi_representation][TAUT_YES as usize];
                    if mobile_base.is_null() {
                        return Err(SourceHeapError::NullPointer);
                    }
                    input[1] = mobile_base.offset(component as i64)?;
                }
            } else if selected_mobile == TAUT_NON as usize {
                let mobile_base = pOneInput.pInpInChI[inchi_representation][TAUT_YES as usize];
                if mobile_base.is_null() {
                    return Err(SourceHeapError::NullPointer);
                }
                let mobile_pointer = mobile_base.offset(component as i64)?;
                let mobile = heap
                    .slice(mobile_pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if mobile.nNumberOfAtoms != 0 && mobile.bDeleted == 0 {
                    selected_mobile = TAUT_YES as usize;
                    input[0] = mobile_pointer;
                }
            }
            let structure_base = pStruct[inchi_representation][selected_mobile];
            if structure_base.is_null() {
                return Err(SourceHeapError::NullPointer);
            }
            let structure = heap
                .slice(structure_base.offset(component as i64)?.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let mut component_flags = [0; TAUT_NUM as usize];
            let mut current_removed = COMPONENT_REM_PROTONS::default();
            ret = CompareOneOrigInchiToRevInChI(
                heap,
                Some(&structure),
                input,
                selected_mobile as i32,
                component as i32,
                num_inp,
                szCurHdr,
                &mut current_removed,
                &mut component_flags,
            )?;
            if ret >= 0 {
                total_removed.nNumRemovedProtons = total_removed
                    .nNumRemovedProtons
                    .wrapping_add(current_removed.nNumRemovedProtons);
                for isotope in 0..NUM_H_ISOTOPES as usize {
                    total_removed.nNumRemovedIsotopicH[isotope] =
                        total_removed.nNumRemovedIsotopicH[isotope]
                            .wrapping_add(current_removed.nNumRemovedIsotopicH[isotope]);
                }
                for tautomer in 0..TAUT_NUM as usize {
                    pOneInput.CompareInchiFlags[0][tautomer] |= component_flags[tautomer];
                }
            } else {
                return Ok(ret);
            }
        }
    }
    if mobile_h == TAUT_YES as usize {
        let proton_balance = &pOneInput.nNumProtons[inchi_representation][mobile_h];
        if !proton_balance.pNumProtons.is_null() {
            ret = RI_ERR_PROGR;
        } else {
            if total_removed.nNumRemovedProtons != proton_balance.nNumRemovedProtons {
                pOneInput.CompareInchiFlags[0][TAUT_YES as usize] |=
                    tagInchiCompareDiffBits_INCHIDIFF_MOBH_PROTONS as INCHI_MODE;
            }
            for isotope in 0..NUM_H_ISOTOPES as usize {
                if total_removed.nNumRemovedIsotopicH[isotope]
                    != proton_balance.nNumRemovedIsotopicH[isotope]
                {
                    pOneInput.CompareInchiFlags[0][TAUT_YES as usize] |=
                        tagInchiCompareDiffBits_INCHIDIFF_MOB_ISO_H as INCHI_MODE;
                }
            }
        }
    }
    Ok(ret)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn CompareAllDisconnectedOrigInchiToRevInChI(
    heap: &mut SourceHeap,
    pStruct: &[[SourceMutPointer<StrFromINChI>; TAUT_NUM as usize]; INCHI_NUM as usize],
    pOneInput: &mut InpInChI,
    bHasSomeFixedH: i32,
    _num_inp: i64,
    _szCurHdr: SourceMutPointer<i8>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1683 CompareAllDisconnectedOrigInchiToRevInChI
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int CompareAllDisconnectedOrigInchiToRevInChI(StrFromINChI* pStruct[INCHI_NUM][TAUT_NUM],
    InpInChI* pOneInput,
    int bHasSomeFixedH,
    long num_inp,
    char* szCurHdr)
{
    int i, k, m, n, iInChI, iMobileH, bMobileH, ifk;
    int num_components_D, num_components_R;
    int nNumCompHaveSeparateProtons_D, nNumCompHaveSeparateProtons_R;
    int num_fragments_D, num_fragments_R, num_fragments_DR, num_fragments, iComponent, ret;
    int ifInChI, ifMobileH, bfMobileH, nLink;
    COMPONENT_REM_PROTONS nNumRemovedProtons_D;     /* removed from the disconnected layer of the Input InChI */
    COMPONENT_REM_PROTONS nNumRemovedProtons_D_all; /* if only totals are avalable */
    COMPONENT_REM_PROTONS nNumRemovedProtons_R; /* removed from disconnected layer of the reconstructed struct */
    COMPONENT_REM_PROTONS nNumRemovedProtons_R_all;
    INCHI_MODE  CompareInchiFlags[TAUT_NUM];
    StrFromINChI* pStruct1;
    INChI_Aux* pINChI_Aux;
    INCHI_SORT* pINChISort1 = NULL; /* from reversed structure */
    INCHI_SORT* pINChISort2 = NULL; /* original input InChI */
    int        nNumNonTaut1 = 0, nNumNonTaut2 = 0;

    ret = 0;
    memset(pOneInput->CompareInchiFlags[1], 0, sizeof(pOneInput->CompareInchiFlags[1])); /* djb-rwth: memset_s C11/Annex K variant? */

    /* count components that are not subject to disconnection */
    if (!pOneInput->nNumComponents[INCHI_REC][TAUT_YES] &&
        !pOneInput->nNumComponents[INCHI_REC][TAUT_NON])
    {
        return 0; /* nothing to do */
    }

    memset(&nNumRemovedProtons_D, 0, sizeof(nNumRemovedProtons_D)); /* djb-rwth: memset_s C11/Annex K variant? */
    memset(&nNumRemovedProtons_R, 0, sizeof(nNumRemovedProtons_R)); /* djb-rwth: memset_s C11/Annex K variant? */
    memset(&nNumRemovedProtons_D_all, 0, sizeof(nNumRemovedProtons_D_all)); /* djb-rwth: memset_s C11/Annex K variant? */
    memset(&nNumRemovedProtons_R_all, 0, sizeof(nNumRemovedProtons_R_all)); /* djb-rwth: memset_s C11/Annex K variant? */
    memset(CompareInchiFlags, 0, sizeof(CompareInchiFlags)); /* djb-rwth: memset_s C11/Annex K variant? */

    num_components_D = inchi_max(pOneInput->nNumComponents[INCHI_BAS][TAUT_YES],
        pOneInput->nNumComponents[INCHI_BAS][TAUT_NON]);
    num_components_R = inchi_max(pOneInput->nNumComponents[INCHI_REC][TAUT_YES],
        pOneInput->nNumComponents[INCHI_REC][TAUT_NON]);

    /***********************************************************************************************/
    /* InpInChI: count fragments -- disconnected components that do not match reconnected          */
    /* Accumulate removed H and isotopic H from ALL Fixed-H disconnected components except deleted */
    /* This segment collects info from the original InChI                                          */
    /***********************************************************************************************/
    /*---- Original InChI ----*/

    num_fragments_D = 0;
    iInChI = INCHI_BAS;
    iMobileH = bHasSomeFixedH ? !pOneInput->nNumComponents[iInChI][TAUT_NON] : TAUT_YES;
    nNumCompHaveSeparateProtons_D = 0;

    /* in case of Mobile-H components here are the proton totals from the original InChI disconn. layer */
    nNumRemovedProtons_D.nNumRemovedProtons = pOneInput->nNumProtons[iInChI][TAUT_YES].nNumRemovedProtons;
    memcpy(nNumRemovedProtons_D.nNumRemovedIsotopicH,
        pOneInput->nNumProtons[iInChI][TAUT_YES].nNumRemovedIsotopicH,
        sizeof(nNumRemovedProtons_D.nNumRemovedIsotopicH)); /* total for the disconnected layer */
    for (k = 0; k < num_components_D; k++)
    {
        bMobileH = iMobileH;
        if (!bInpInchiComponentExists(pOneInput, iInChI, bMobileH, k))
        {
            if (bInpInchiComponentExists(pOneInput, iInChI, TAUT_YES, k))
            {
                bMobileH = TAUT_YES;
            }
            else
            {
                continue; /* component is missing ??? */
            }
        }
        if (0 > (nLink = pOneInput->pInpInChI[iInChI][bMobileH][k].nLink))
        {
            /* component in Disconnected layer is linked to the identical one in the Reconnected layer */
            if (pOneInput->nNumProtons[INCHI_REC][TAUT_YES].pNumProtons)
            {
                nNumCompHaveSeparateProtons_D++;
                nLink = -(1 + nLink);
                nNumRemovedProtons_D.nNumRemovedProtons += pOneInput->nNumProtons[INCHI_REC][TAUT_YES].pNumProtons[nLink].nNumRemovedProtons;
                for (m = 0; m < NUM_H_ISOTOPES; m++)
                {
                    nNumRemovedProtons_D.nNumRemovedIsotopicH[m] += pOneInput->nNumProtons[INCHI_REC][TAUT_YES].pNumProtons[nLink].nNumRemovedIsotopicH[m];
                }
            }
            continue; /* same as reconnected */
        }
        /* component in the reconnected layer that was disconnected */
        nNumNonTaut2 += (bMobileH == TAUT_NON);
        if (pOneInput->nNumProtons[iInChI][TAUT_YES].pNumProtons)
        {
            nNumCompHaveSeparateProtons_D++;
            nNumRemovedProtons_D.nNumRemovedProtons += pOneInput->nNumProtons[iInChI][TAUT_YES].pNumProtons[k].nNumRemovedProtons;
            for (m = 0; m < NUM_H_ISOTOPES; m++)
            {
                nNumRemovedProtons_D.nNumRemovedIsotopicH[m] += pOneInput->nNumProtons[iInChI][TAUT_YES].pNumProtons[k].nNumRemovedIsotopicH[m];
            }
        }
        num_fragments_D++; /* number of disconnected fragments from original reconnected structure */
    }
    /* in case of Mobile-H components here are the proton totals from the original InChI */
    /*
    nNumRemovedProtons_D_all.nNumRemovedProtons = pOneInput->nNumProtons[iInChI][TAUT_YES].nNumRemovedProtons;
    memcpy( nNumRemovedProtons_D_all.nNumRemovedIsotopicH,
            pOneInput->nNumProtons[iInChI][TAUT_YES].nNumRemovedIsotopicH,
            sizeof(nNumRemovedProtons_D_all.nNumRemovedIsotopicH) );

    */
    /****************************************************************************************************/
    /* count fragments in reconstructed reconnected structure                                           */
    /* accumulate removed H and isotopic H from ALL reconstructed reconnected components except deleted */
    /* This segment collects info from the reconstructed structure InChI                                */
    /****************************************************************************************************/
    /*---- InChI from the reconstructed reconnected structure ----*/
    num_fragments_R = 0;
    iInChI = INCHI_REC;
    iMobileH = bHasSomeFixedH ? !pOneInput->nNumComponents[iInChI][TAUT_NON] : TAUT_YES;
    nNumCompHaveSeparateProtons_R = 0;
    for (k = 0; k < num_components_R; k++)
    {
        bMobileH = iMobileH;
        if (!bInpInchiComponentExists(pOneInput, iInChI, bMobileH, k))
        {
            if (bInpInchiComponentExists(pOneInput, iInChI, TAUT_YES, k))
            {
                bMobileH = TAUT_YES;
            }
            else
            {
                continue; /* component is missing ??? (Deleted proton in Mobile-H layer) */
            }
        }
        if (0 < pOneInput->pInpInChI[iInChI][bMobileH][k].nLink)
        {
            /* this reconstructed reconnected component was NOT DISCONNECTED */
            /* same component is in the disconnected layer, it has no metal atoms or is an isolated metal atom */
            pStruct1 = pStruct[iInChI][bMobileH] + k;
            ifMobileH = TAUT_YES;  /* Mobile-H Aux_Info contains number removed protons */
            ifInChI = INCHI_BAS; /* this component cannot be reconnected */
            ifk = 0;         /* 0th component since it is InChI of a single component */
            /* The statement in the following line is *WRONG*, component number mixed with bMobileH:  */
            /* in RevInchi, when only Mobile-H is present then its only non-NULL InChI has index 0==TAUT_NON */
            if (bRevInchiComponentExists(pStruct1, ifInChI, ifMobileH, ifk))
            {
                /* count protons */
                pINChI_Aux = pStruct1->RevInChI.pINChI_Aux[ifInChI][ifk][ifMobileH];
                if (pINChI_Aux)
                {
                    nNumRemovedProtons_R.nNumRemovedProtons += pINChI_Aux->nNumRemovedProtons;
                    for (m = 0; m < NUM_H_ISOTOPES; m++)
                    {
                        nNumRemovedProtons_R.nNumRemovedIsotopicH[m] += pINChI_Aux->nNumRemovedIsotopicH[m];
                    }
                }
            }
            nNumCompHaveSeparateProtons_R += bRevInchiComponentExists(pStruct1, ifInChI, ALT_TAUT(ifMobileH), ifk);
            continue; /* same as disconnected, has no metal atoms */
        }
        /* this reconstructed reconnected component WAS DISCONNECTED; check its fragments */
        /* it does not have same component in the disconnected layer */
        pStruct1 = pStruct[iInChI][bMobileH] + k;
        num_fragments = pStruct1->RevInChI.num_components[INCHI_BAS];
        ifInChI = INCHI_BAS; /* disconnected layer */
        ifMobileH = bHasSomeFixedH ? TAUT_NON : TAUT_YES;
        for (ifk = 0; ifk < num_fragments; ifk++)
        {
            bfMobileH = ifMobileH;
            if (!bRevInchiComponentExists(pStruct1, ifInChI, bfMobileH, ifk))
            {
                if (bRevInchiComponentExists(pStruct1, ifInChI, TAUT_YES, ifk))
                {
                    bfMobileH = TAUT_YES;
                }
                else
                {
                    continue; /* fragment does not exist ??? */
                }
            }
            nNumNonTaut1 += (bfMobileH == TAUT_NON);
            nNumCompHaveSeparateProtons_R += (bfMobileH == TAUT_NON);
            /* count protons from fragments made by metal disconnection */
            pINChI_Aux = pStruct1->RevInChI.pINChI_Aux[ifInChI][ifk][TAUT_YES];
            if (pINChI_Aux)
            {
                nNumRemovedProtons_R.nNumRemovedProtons += pINChI_Aux->nNumRemovedProtons;
                for (m = 0; m < NUM_H_ISOTOPES; m++)
                {
                    nNumRemovedProtons_R.nNumRemovedIsotopicH[m] += pINChI_Aux->nNumRemovedIsotopicH[m];
                }
            }
            num_fragments_R++; /* number of disconnected fragments from reconstructed reconnected structure */
        }
    }

    /*---------------- special treatment of the last reconstructed component -----------------*/
    /*---------------- this may contain separate protons added by the reconstruction ---------*/
    k = num_components_R - 1;
    pStruct1 = pStruct[iInChI][iMobileH] + k;
    if (iMobileH == TAUT_YES && !bHasSomeFixedH &&
        bInpInchiComponentDeleted(pOneInput, iInChI, iMobileH, k) &&
        (num_fragments = pStruct1->RevInChI.num_components[INCHI_BAS]))
    {

        ifInChI = INCHI_BAS; /* disconnected layer */
        ifMobileH = TAUT_YES;
        for (ifk = 0; ifk < num_fragments; ifk++)
        {
            bfMobileH = ifMobileH;
            if (!bRevInchiComponentDeleted(pStruct1, ifInChI, bfMobileH, ifk))
            {
                continue; /* fragment does exist ??? Should not happen */
            }
            /*
            nNumNonTaut1           += (bfMobileH == TAUT_NON);
            nNumCompHaveSeparateProtons_R += (bfMobileH == TAUT_NON);
            */
            /* count protons from fragments made by metal disconnection */
            pINChI_Aux = pStruct1->RevInChI.pINChI_Aux[ifInChI][ifk][TAUT_YES];
            if (pINChI_Aux)
            {
                nNumRemovedProtons_R.nNumRemovedProtons += pINChI_Aux->nNumRemovedProtons;
                for (m = 0; m < NUM_H_ISOTOPES; m++)
                {
                    nNumRemovedProtons_R.nNumRemovedIsotopicH[m] += pINChI_Aux->nNumRemovedIsotopicH[m];
                }
            }
            /*num_fragments_R ++;*/ /* number of disconnected fragments from reconstructed reconnected structure */
        }
    }

    num_fragments_DR = inchi_max(num_fragments_D, num_fragments_R);
    /* in case of correct reconstruction, num_fragments_D, num_fragments_R */
    if (!num_fragments_DR)
    {
        return 0; /* no component was disconnected */
    }
    if (num_fragments_D != num_fragments_R)
    {
        for (i = 0; i < TAUT_NUM; i++)
        {
            if (pOneInput->nNumComponents[INCHI_BAS][i])
            {
                pOneInput->CompareInchiFlags[1][i] |= INCHIDIFF_PROBLEM;
            }
        }
        return 1; /* severe error */
    }

    pINChISort1 = (INCHI_SORT*)inchi_calloc(num_fragments_DR, sizeof(pINChISort1[0]));
    pINChISort2 = (INCHI_SORT*)inchi_calloc(num_fragments_DR, sizeof(pINChISort2[0]));
    if (!pINChISort1 || !pINChISort2)
    {
        ret = RI_ERR_ALLOC;
        goto exit_function;
    }

    /* accumulate original InChI of fragments -- disconnected components that do not match reconnected */
    iInChI = INCHI_BAS;
    iMobileH = bHasSomeFixedH ? !pOneInput->nNumComponents[iInChI][TAUT_NON] : TAUT_YES;
    for (k = n = 0; k < num_components_D; k++)
    {
        bMobileH = iMobileH;
        if (!bInpInchiComponentExists(pOneInput, iInChI, bMobileH, k))
        {
            if (bInpInchiComponentExists(pOneInput, iInChI, TAUT_YES, k))
            {
                bMobileH = TAUT_YES;
            }
            else
            {
                continue; /* component is missing ??? (Deleted proton in Mobile-H layer) */
            }
        }
        if (0 > pOneInput->pInpInChI[iInChI][bMobileH][k].nLink)
        {
            continue; /* same as reconnected */
        }
        /* the component exists in disconnected layer of the orig. InChI only: it is a fragment */
        pINChISort2[n].pINChI[bMobileH] = pOneInput->pInpInChI[iInChI][bMobileH] + k;
        if (bMobileH == TAUT_NON &&
            (bInpInchiComponentExists(pOneInput, iInChI, TAUT_YES, k) ||
                bInpInchiComponentDeleted(pOneInput, iInChI, TAUT_YES, k)))
        {
            pINChISort2[n].pINChI[TAUT_YES] = pOneInput->pInpInChI[iInChI][TAUT_YES] + k;
        }
        /* the last sort key is a number of removed protons */
        pINChISort2[n].ord_number = pOneInput->nNumProtons[iInChI][TAUT_YES].pNumProtons ?
            pOneInput->nNumProtons[iInChI][TAUT_YES].pNumProtons[k].nNumRemovedProtons : 0;
        pINChISort2[n].n1 = k;  /* orig. InChI disconnected layer component number */
        pINChISort2[n].n2 = -1; /* no fragment index */
        n++;
    }

    /* accumulate fragments from the reconstructed structure */
    iInChI = INCHI_REC;
    iMobileH = bHasSomeFixedH ? !pOneInput->nNumComponents[iInChI][TAUT_NON] : TAUT_YES;
    for (k = n = 0; k < num_components_R; k++)
    {
        bMobileH = iMobileH;
        if (!bInpInchiComponentExists(pOneInput, iInChI, bMobileH, k))
        {
            if (bInpInchiComponentExists(pOneInput, iInChI, TAUT_YES, k))
            {
                bMobileH = TAUT_YES;
            }
            else
            {
                continue; /* component is missing ??? (Deleted proton in Mobile-H layer) */
            }
        }
        /* the reconstructed structure */
        if (0 < pOneInput->pInpInChI[iInChI][bMobileH][k].nLink)
        {
            continue; /* same as disconnected, has no metal atoms */
        }
        /* this reconstructed structure was disconnected */
        pStruct1 = pStruct[iInChI][bMobileH] + k;
        num_fragments = pStruct1->RevInChI.num_components[INCHI_BAS];
        ifInChI = INCHI_BAS;
        ifMobileH = bHasSomeFixedH ? TAUT_NON : TAUT_YES;
        for (i = 0; i < num_fragments; i++)
        {
            bfMobileH = ifMobileH;
            if (!bRevInchiComponentExists(pStruct1, ifInChI, bfMobileH, i))
            {
                if (bRevInchiComponentExists(pStruct1, ifInChI, TAUT_YES, i))
                {
                    bfMobileH = TAUT_YES;
                }
                else
                {
                    continue; /* component is missing ??? */
                }
            }
            pINChISort1[n].pINChI[bfMobileH] = pStruct1->RevInChI.pINChI[ifInChI][i][bfMobileH];
            if (bfMobileH == TAUT_NON /*&& bRevInchiComponentExists( pStruct1, ifInChI, TAUT_YES, i )*/)
            {
                pINChISort1[n].pINChI[TAUT_YES] = pStruct1->RevInChI.pINChI[ifInChI][i][TAUT_YES];
                /* remove Fixed-H InChI if is is identical to Mobile-H */
                /* do it exactly same way the identical components were removed from InpInChI */
                if (!CompareReversedINChI(pINChISort1[n].pINChI[bfMobileH],
                    pINChISort1[n].pINChI[TAUT_YES], NULL, NULL))
                {
                    pINChISort1[n].pINChI[bfMobileH] = NULL; /* remove Fixed-H layer */
                }
                else
                {
                    pINChISort1[n].ord_number = pStruct1->RevInChI.pINChI_Aux[ifInChI][i][TAUT_YES]->nNumRemovedProtons;
                }
            }

            pINChISort1[n].n1 = k;  /* reconstructed reconnected structure component index */
            pINChISort1[n].n2 = i;  /* index of a fragment made out of this component */
            n++;
        }
    }

    /* sort fragment InChI before comparing them */
    qsort(pINChISort1, num_fragments_D, sizeof(pINChISort1[0]), CompINChITaut2);
    qsort(pINChISort2, num_fragments_R, sizeof(pINChISort2[0]), CompINChITaut2);

    /* compare fragments -- components present in disconnected layer only */
    for (iComponent = 0; iComponent < num_fragments_DR; iComponent++)
    {
        INChI* pInChI1[TAUT_NUM]; /* from reversed structure */
        INChI* pInChI2[TAUT_NUM]; /* original input InChI */
        for (i = 0; i < TAUT_NUM; i++)
        {
            pInChI1[i] = pINChISort1[iComponent].pINChI[i];
            pInChI2[i] = pINChISort2[iComponent].pINChI[i];
        }
        CompareTwoPairsOfInChI(pInChI1, pInChI2, !bHasSomeFixedH, CompareInchiFlags);
    }

    if ( /*nNumNonTaut1 && nNumNonTaut2 &&*/ bHasSomeFixedH)
    {
        if (nNumCompHaveSeparateProtons_D || nNumCompHaveSeparateProtons_R)
        {
            /* for each component, compare number removed protons */
            /* comparison does not make sense if Disconnected Fixed-H layer is not present */
            for (iComponent = 0; iComponent < num_fragments_DR; iComponent++)
            {
                NUM_H   nNumRemovedIsotopicH1[NUM_H_ISOTOPES];
                NUM_H   nNumRemovedIsotopicH2[NUM_H_ISOTOPES];

                memset(nNumRemovedIsotopicH1, 0, sizeof(nNumRemovedIsotopicH1)); /* djb-rwth: memset_s C11/Annex K variant? */
                memset(nNumRemovedIsotopicH2, 0, sizeof(nNumRemovedIsotopicH2)); /* djb-rwth: memset_s C11/Annex K variant? */
                /* compare removed protons */
                if (pINChISort1[iComponent].ord_number != pINChISort2[iComponent].ord_number)
                {
                    CompareInchiFlags[TAUT_YES] |= INCHIDIFF_MOBH_PROTONS; /* diff number of removed protons */
                }
                /* also compare removed isotopic atoms H */
                k = pINChISort2[iComponent].n1; /* input InChI, OneInput */
                if (pOneInput->nNumProtons[INCHI_BAS][TAUT_YES].pNumProtons)
                {
                    memcpy(nNumRemovedIsotopicH2,
                        pOneInput->nNumProtons[INCHI_BAS][TAUT_YES].pNumProtons[k].nNumRemovedIsotopicH,
                        sizeof(nNumRemovedIsotopicH2));
                }
                /* get fragments of reconstructed structure removed protons info */
                k = pINChISort1[iComponent].n1; /* restored component number */
                i = pINChISort1[iComponent].n2; /* subcomponent number */
                iInChI = INCHI_REC;
                iMobileH = bHasSomeFixedH ? !pOneInput->nNumComponents[iInChI][TAUT_NON] : TAUT_YES;
                bMobileH = iMobileH;
                if (!bInpInchiComponentExists(pOneInput, iInChI, bMobileH, k))
                {
                    if (bInpInchiComponentExists(pOneInput, iInChI, TAUT_YES, k))
                    {
                        bMobileH = TAUT_YES;
                    }
                    else
                    {
                        goto compare_iso_H;
                    }
                }
                if (pOneInput->pInpInChI[iInChI][bMobileH][k].nLink)
                {
                    continue;
                    /*
                    ret = RI_ERR_PROGR;
                    goto exit_function;
                    */
                }
                pStruct1 = pStruct[iInChI][bMobileH] + k;
                num_fragments = pStruct1->RevInChI.num_components[INCHI_BAS];
                ifInChI = INCHI_BAS;
                ifMobileH = bHasSomeFixedH ? TAUT_NON : TAUT_YES;
                if (i < num_fragments)
                {
                    bfMobileH = ifMobileH;
                    if (!bRevInchiComponentExists(pStruct1, ifInChI, bfMobileH, i))
                    {
                        if (!bRevInchiComponentExists(pStruct1, ifInChI, TAUT_YES, i))
                            /* djb-rwth: removing redundant code */
                        {
                            goto compare_iso_H;
                        }
                    }
                    memcpy(nNumRemovedIsotopicH1,
                        pStruct1->RevInChI.pINChI_Aux[ifInChI][i][TAUT_YES]->nNumRemovedIsotopicH,
                        sizeof(nNumRemovedIsotopicH1));
                }
            compare_iso_H:
                if (memcmp(nNumRemovedIsotopicH1, nNumRemovedIsotopicH2, sizeof(nNumRemovedIsotopicH1)))
                {
                    CompareInchiFlags[TAUT_YES] |= INCHIDIFF_REM_ISO_H;
                }
            }
        }
    }
    else /*if ( !nNumNonTaut1 && !nNumNonTaut2 || !bHasSomeFixedH )*/
    {
        /* compare totals for removed protons and isotopic H */
        if (pOneInput->nNumProtons[INCHI_BAS][TAUT_YES].nNumRemovedProtons !=
            nNumRemovedProtons_R.nNumRemovedProtons)
        {
            CompareInchiFlags[TAUT_YES] |= INCHIDIFF_MOBH_PROTONS;
        }
        if (memcmp(pOneInput->nNumProtons[INCHI_BAS][TAUT_YES].nNumRemovedIsotopicH,
            nNumRemovedProtons_R.nNumRemovedIsotopicH,
            sizeof(nNumRemovedProtons_R.nNumRemovedIsotopicH)))
        {
            CompareInchiFlags[TAUT_YES] |= INCHIDIFF_REM_ISO_H;
        }
    }

    if (!nNumNonTaut1 == !nNumNonTaut2)
    {
        ; /* difference if(nNumNonTaut1 != nNumNonTaut2) will be caught in InChI comparison */
    }
    else
    {
        if (nNumNonTaut1)
        {
            /* reconstructed has Fixed-H while the original has not: extra Fixed-H layer */
            CompareInchiFlags[TAUT_YES] |= INCHIDIFF_WRONG_TAUT;
        }
        else
        {
            /* the original InChI has Fixed-H while the reconstructed one has not: missing Fixed-H layer */
            CompareInchiFlags[TAUT_YES] |= INCHIDIFF_NO_TAUT;
        }
    }
    for (i = 0; i < TAUT_NUM; i++)
    {
        pOneInput->CompareInchiFlags[1][i] |= CompareInchiFlags[i];
    }

    /* compare totals */
    if (nNumRemovedProtons_R.nNumRemovedProtons != nNumRemovedProtons_D.nNumRemovedProtons)
    {
        CompareInchiFlags[TAUT_YES] |= INCHIDIFF_MOBH_PROTONS; /* diff number of removed protons */
    }
    if (memcmp(nNumRemovedProtons_R.nNumRemovedIsotopicH,
        nNumRemovedProtons_D.nNumRemovedIsotopicH,
        sizeof(nNumRemovedProtons_D.nNumRemovedIsotopicH)))
    {
        CompareInchiFlags[TAUT_YES] |= INCHIDIFF_REM_ISO_H;
    }

exit_function:

    if (pINChISort1)
    {
        inchi_free(pINChISort1);
    }
    if (pINChISort2)
    {
        inchi_free(pINChISort2);
    }

    return ret;
}
    */
    // END INCHI C FUNCTION: CompareAllDisconnectedOrigInchiToRevInChI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CompareAllDisconnectedOrigInchiToRevInChI
    // INCHI✔️❌: READ_INCHI_STRING=1 includes this production helper.
    // INCHI✔️❌: COMPILE_ANSI_ONLY and TARGET_API_LIB do not alter this function body.
    // INCHI✔️❌: inchi_max and ALT_TAUT are reproduced as active header macros.
    // INCHI✔️❌: qsort uses the already ported CompINChITaut2 comparator.
    // END INCHI ACTIVE MACRO CONFIGURATION: CompareAllDisconnectedOrigInchiToRevInChI

    let component_pointer = |input: &InpInChI,
                             representation: usize,
                             mobile_h: usize,
                             component: i32|
     -> Result<SourceMutPointer<INChI>, SourceHeapError> {
        if component < 0 {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let base = input.pInpInChI[representation][mobile_h];
        if base.is_null() {
            return Err(SourceHeapError::NullPointer);
        }
        base.offset(i64::from(component))
    };
    let structure_value = |heap: &SourceHeap,
                           base: SourceMutPointer<StrFromINChI>,
                           component: i32|
     -> Result<StrFromINChI, SourceHeapError> {
        if base.is_null() {
            return Err(SourceHeapError::NullPointer);
        }
        if component < 0 {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        Ok(heap
            .slice(base.offset(i64::from(component))?.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone())
    };
    let rev_inchi_pointer = |heap: &SourceHeap,
                             structure: &StrFromINChI,
                             representation: usize,
                             component: i32,
                             mobile_h: usize|
     -> Result<SourceMutPointer<INChI>, SourceHeapError> {
        if component < 0 {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let rows = structure.RevInChI.pINChI[representation];
        if rows.is_null() {
            return Err(SourceHeapError::NullPointer);
        }
        Ok(*heap
            .slice(rows.as_const())?
            .get(usize::try_from(component).map_err(|_| SourceHeapError::SourceIntegerOverflow)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .get(mobile_h)
            .ok_or(SourceHeapError::PointerOutOfBounds)?)
    };
    let rev_aux_pointer = |heap: &SourceHeap,
                           structure: &StrFromINChI,
                           representation: usize,
                           component: i32,
                           mobile_h: usize|
     -> Result<SourceMutPointer<INChI_Aux>, SourceHeapError> {
        if component < 0 {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let rows = structure.RevInChI.pINChI_Aux[representation];
        if rows.is_null() {
            return Err(SourceHeapError::NullPointer);
        }
        Ok(*heap
            .slice(rows.as_const())?
            .get(usize::try_from(component).map_err(|_| SourceHeapError::SourceIntegerOverflow)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .get(mobile_h)
            .ok_or(SourceHeapError::PointerOutOfBounds)?)
    };
    let add_removed = |total: &mut COMPONENT_REM_PROTONS, value: &COMPONENT_REM_PROTONS| {
        total.nNumRemovedProtons = total
            .nNumRemovedProtons
            .wrapping_add(value.nNumRemovedProtons);
        for isotope in 0..NUM_H_ISOTOPES as usize {
            total.nNumRemovedIsotopicH[isotope] = total.nNumRemovedIsotopicH[isotope]
                .wrapping_add(value.nNumRemovedIsotopicH[isotope]);
        }
    };
    let add_aux = |total: &mut COMPONENT_REM_PROTONS, value: &INChI_Aux| {
        total.nNumRemovedProtons = total
            .nNumRemovedProtons
            .wrapping_add(value.nNumRemovedProtons);
        for isotope in 0..NUM_H_ISOTOPES as usize {
            total.nNumRemovedIsotopicH[isotope] = total.nNumRemovedIsotopicH[isotope]
                .wrapping_add(value.nNumRemovedIsotopicH[isotope]);
        }
    };

    pOneInput.CompareInchiFlags[1] = [0; TAUT_NUM as usize];
    if pOneInput.nNumComponents[INCHI_REC as usize][TAUT_YES as usize] == 0
        && pOneInput.nNumComponents[INCHI_REC as usize][TAUT_NON as usize] == 0
    {
        return Ok(0);
    }

    let num_components_d = pOneInput.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize]
        .max(pOneInput.nNumComponents[INCHI_BAS as usize][TAUT_NON as usize]);
    let num_components_r = pOneInput.nNumComponents[INCHI_REC as usize][TAUT_YES as usize]
        .max(pOneInput.nNumComponents[INCHI_REC as usize][TAUT_NON as usize]);
    let mut removed_d = COMPONENT_REM_PROTONS {
        nNumRemovedProtons: pOneInput.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize]
            .nNumRemovedProtons,
        nNumRemovedIsotopicH: pOneInput.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize]
            .nNumRemovedIsotopicH,
    };
    let mut removed_r = COMPONENT_REM_PROTONS::default();
    let mut compare_flags = [0; TAUT_NUM as usize];
    let mut num_fragments_d = 0_i32;
    let input_mobile_d = if bHasSomeFixedH != 0 {
        i32::from(pOneInput.nNumComponents[INCHI_BAS as usize][TAUT_NON as usize] == 0)
    } else {
        TAUT_YES as i32
    };
    let mut separate_d = 0_i32;
    let mut non_taut_2 = 0_i32;
    for component in 0..num_components_d.max(0) {
        let mut mobile_h = input_mobile_d;
        if bInpInchiComponentExists(
            heap,
            pOneInput,
            INCHI_BAS as i32,
            mobile_h,
            component,
        )? == 0
        {
            if bInpInchiComponentExists(
                heap,
                pOneInput,
                INCHI_BAS as i32,
                TAUT_YES as i32,
                component,
            )? != 0
            {
                mobile_h = TAUT_YES as i32;
            } else {
                continue;
            }
        }
        let pointer = component_pointer(
            pOneInput,
            INCHI_BAS as usize,
            mobile_h as usize,
            component,
        )?;
        let value = heap
            .slice(pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        if value.nLink < 0 {
            let protons = pOneInput.nNumProtons[INCHI_REC as usize][TAUT_YES as usize].pNumProtons;
            if !protons.is_null() {
                separate_d = separate_d.wrapping_add(1);
                let link = value.nLink.wrapping_add(1).wrapping_neg();
                let proton = heap
                    .slice(protons.as_const())?
                    .get(usize::try_from(link).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone();
                add_removed(&mut removed_d, &proton);
            }
            continue;
        }
        non_taut_2 = non_taut_2.wrapping_add(i32::from(mobile_h == TAUT_NON as i32));
        let protons = pOneInput.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize].pNumProtons;
        if !protons.is_null() {
            separate_d = separate_d.wrapping_add(1);
            let proton = heap
                .slice(protons.as_const())?
                .get(usize::try_from(component).map_err(|_| SourceHeapError::SourceIntegerOverflow)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            add_removed(&mut removed_d, &proton);
        }
        num_fragments_d = num_fragments_d.wrapping_add(1);
    }

    let input_mobile_r = if bHasSomeFixedH != 0 {
        i32::from(pOneInput.nNumComponents[INCHI_REC as usize][TAUT_NON as usize] == 0)
    } else {
        TAUT_YES as i32
    };
    let mut num_fragments_r = 0_i32;
    let mut separate_r = 0_i32;
    let mut non_taut_1 = 0_i32;
    for component in 0..num_components_r.max(0) {
        let mut mobile_h = input_mobile_r;
        if bInpInchiComponentExists(
            heap,
            pOneInput,
            INCHI_REC as i32,
            mobile_h,
            component,
        )? == 0
        {
            if bInpInchiComponentExists(
                heap,
                pOneInput,
                INCHI_REC as i32,
                TAUT_YES as i32,
                component,
            )? != 0
            {
                mobile_h = TAUT_YES as i32;
            } else {
                continue;
            }
        }
        let input_pointer = component_pointer(
            pOneInput,
            INCHI_REC as usize,
            mobile_h as usize,
            component,
        )?;
        let input_value = heap
            .slice(input_pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let structure = structure_value(
            heap,
            pStruct[INCHI_REC as usize][mobile_h as usize],
            component,
        )?;
        if input_value.nLink > 0 {
            if bRevInchiComponentExists(
                heap,
                Some(&structure),
                INCHI_BAS as i32,
                TAUT_YES as i32,
                0,
            )? != 0
            {
                let aux = rev_aux_pointer(
                    heap,
                    &structure,
                    INCHI_BAS as usize,
                    0,
                    TAUT_YES as usize,
                )?;
                if !aux.is_null() {
                    let value = heap
                        .slice(aux.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    add_aux(&mut removed_r, &value);
                }
            }
            separate_r = separate_r.wrapping_add(bRevInchiComponentExists(
                heap,
                Some(&structure),
                INCHI_BAS as i32,
                TAUT_NON as i32,
                0,
            )?);
            continue;
        }
        let fragments = structure.RevInChI.num_components[INCHI_BAS as usize];
        let preferred = if bHasSomeFixedH != 0 {
            TAUT_NON as i32
        } else {
            TAUT_YES as i32
        };
        for fragment in 0..fragments.max(0) {
            let mut fragment_mobile = preferred;
            if bRevInchiComponentExists(
                heap,
                Some(&structure),
                INCHI_BAS as i32,
                fragment_mobile,
                fragment,
            )? == 0
            {
                if bRevInchiComponentExists(
                    heap,
                    Some(&structure),
                    INCHI_BAS as i32,
                    TAUT_YES as i32,
                    fragment,
                )? != 0
                {
                    fragment_mobile = TAUT_YES as i32;
                } else {
                    continue;
                }
            }
            non_taut_1 =
                non_taut_1.wrapping_add(i32::from(fragment_mobile == TAUT_NON as i32));
            separate_r =
                separate_r.wrapping_add(i32::from(fragment_mobile == TAUT_NON as i32));
            let aux = rev_aux_pointer(
                heap,
                &structure,
                INCHI_BAS as usize,
                fragment,
                TAUT_YES as usize,
            )?;
            if !aux.is_null() {
                let value = heap
                    .slice(aux.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone();
                add_aux(&mut removed_r, &value);
            }
            num_fragments_r = num_fragments_r.wrapping_add(1);
        }
    }

    let last_component = num_components_r.wrapping_sub(1);
    let last_structure = structure_value(
        heap,
        pStruct[INCHI_REC as usize][input_mobile_r as usize],
        last_component,
    )?;
    if input_mobile_r == TAUT_YES as i32
        && bHasSomeFixedH == 0
        && bInpInchiComponentDeleted(
            heap,
            pOneInput,
            INCHI_REC as i32,
            input_mobile_r,
            last_component,
        )? != 0
    {
        let fragments = last_structure.RevInChI.num_components[INCHI_BAS as usize];
        if fragments != 0 {
            for fragment in 0..fragments.max(0) {
                if bRevInchiComponentDeleted(
                    heap,
                    Some(&last_structure),
                    INCHI_BAS as i32,
                    TAUT_YES as i32,
                    fragment,
                )? == 0
                {
                    continue;
                }
                let aux = rev_aux_pointer(
                    heap,
                    &last_structure,
                    INCHI_BAS as usize,
                    fragment,
                    TAUT_YES as usize,
                )?;
                if !aux.is_null() {
                    let value = heap
                        .slice(aux.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    add_aux(&mut removed_r, &value);
                }
            }
        }
    }

    let num_fragments = num_fragments_d.max(num_fragments_r);
    if num_fragments == 0 {
        return Ok(0);
    }
    if num_fragments_d != num_fragments_r {
        for mobile_h in 0..TAUT_NUM as usize {
            if pOneInput.nNumComponents[INCHI_BAS as usize][mobile_h] != 0 {
                pOneInput.CompareInchiFlags[1][mobile_h] |=
                    tagInchiCompareDiffBits_INCHIDIFF_PROBLEM as INCHI_MODE;
            }
        }
        return Ok(1);
    }

    let count = u64::try_from(num_fragments)
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
    let source_size = std::mem::size_of::<INCHI_SORT>() as u64;
    let sort1 = match inchi_calloc::<INCHI_SORT>(heap, count, source_size) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
        Err(error) => return Err(error),
    };
    let sort2 = match inchi_calloc::<INCHI_SORT>(heap, count, source_size) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
        Err(error) => {
            if !sort1.is_null() {
                inchi_free(heap, sort1)?;
            }
            return Err(error);
        }
    };
    if sort1.is_null() || sort2.is_null() {
        if !sort1.is_null() {
            inchi_free(heap, sort1)?;
        }
        if !sort2.is_null() {
            inchi_free(heap, sort2)?;
        }
        return Ok(RI_ERR_ALLOC);
    }

    let execution = (|| -> Result<i32, SourceHeapError> {
        let mut row_index = 0_usize;
        for component in 0..num_components_d.max(0) {
            let mut mobile_h = input_mobile_d;
            if bInpInchiComponentExists(
                heap,
                pOneInput,
                INCHI_BAS as i32,
                mobile_h,
                component,
            )? == 0
            {
                if bInpInchiComponentExists(
                    heap,
                    pOneInput,
                    INCHI_BAS as i32,
                    TAUT_YES as i32,
                    component,
                )? != 0
                {
                    mobile_h = TAUT_YES as i32;
                } else {
                    continue;
                }
            }
            let pointer = component_pointer(
                pOneInput,
                INCHI_BAS as usize,
                mobile_h as usize,
                component,
            )?;
            let input_value = heap
                .slice(pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            if input_value.nLink < 0 {
                continue;
            }
            let mut row = INCHI_SORT::default();
            row.pINChI[mobile_h as usize] = pointer;
            if mobile_h == TAUT_NON as i32
                && (bInpInchiComponentExists(
                    heap,
                    pOneInput,
                    INCHI_BAS as i32,
                    TAUT_YES as i32,
                    component,
                )? != 0
                    || bInpInchiComponentDeleted(
                        heap,
                        pOneInput,
                        INCHI_BAS as i32,
                        TAUT_YES as i32,
                        component,
                    )? != 0)
            {
                row.pINChI[TAUT_YES as usize] = component_pointer(
                    pOneInput,
                    INCHI_BAS as usize,
                    TAUT_YES as usize,
                    component,
                )?;
            }
            let protons =
                pOneInput.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize].pNumProtons;
            if !protons.is_null() {
                row.ord_number = heap
                    .slice(protons.as_const())?
                    .get(
                        usize::try_from(component)
                            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNumRemovedProtons;
            }
            row.n1 = component as i16;
            row.n2 = -1;
            *heap
                .slice_mut(sort2)?
                .get_mut(row_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = row;
            row_index += 1;
        }

        row_index = 0;
        for component in 0..num_components_r.max(0) {
            let mut mobile_h = input_mobile_r;
            if bInpInchiComponentExists(
                heap,
                pOneInput,
                INCHI_REC as i32,
                mobile_h,
                component,
            )? == 0
            {
                if bInpInchiComponentExists(
                    heap,
                    pOneInput,
                    INCHI_REC as i32,
                    TAUT_YES as i32,
                    component,
                )? != 0
                {
                    mobile_h = TAUT_YES as i32;
                } else {
                    continue;
                }
            }
            let input_pointer = component_pointer(
                pOneInput,
                INCHI_REC as usize,
                mobile_h as usize,
                component,
            )?;
            let input_value = heap
                .slice(input_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            if input_value.nLink > 0 {
                continue;
            }
            let structure = structure_value(
                heap,
                pStruct[INCHI_REC as usize][mobile_h as usize],
                component,
            )?;
            let fragments = structure.RevInChI.num_components[INCHI_BAS as usize];
            let preferred = if bHasSomeFixedH != 0 {
                TAUT_NON as i32
            } else {
                TAUT_YES as i32
            };
            for fragment in 0..fragments.max(0) {
                let mut fragment_mobile = preferred;
                if bRevInchiComponentExists(
                    heap,
                    Some(&structure),
                    INCHI_BAS as i32,
                    fragment_mobile,
                    fragment,
                )? == 0
                {
                    if bRevInchiComponentExists(
                        heap,
                        Some(&structure),
                        INCHI_BAS as i32,
                        TAUT_YES as i32,
                        fragment,
                    )? != 0
                    {
                        fragment_mobile = TAUT_YES as i32;
                    } else {
                        continue;
                    }
                }
                let mut row = INCHI_SORT::default();
                row.pINChI[fragment_mobile as usize] = rev_inchi_pointer(
                    heap,
                    &structure,
                    INCHI_BAS as usize,
                    fragment,
                    fragment_mobile as usize,
                )?;
                if fragment_mobile == TAUT_NON as i32 {
                    row.pINChI[TAUT_YES as usize] = rev_inchi_pointer(
                        heap,
                        &structure,
                        INCHI_BAS as usize,
                        fragment,
                        TAUT_YES as usize,
                    )?;
                    let fixed = if row.pINChI[TAUT_NON as usize].is_null() {
                        None
                    } else {
                        heap.slice(row.pINChI[TAUT_NON as usize].as_const())?
                            .first()
                            .cloned()
                    };
                    let mobile = if row.pINChI[TAUT_YES as usize].is_null() {
                        None
                    } else {
                        heap.slice(row.pINChI[TAUT_YES as usize].as_const())?
                            .first()
                            .cloned()
                    };
                    if CompareReversedINChI(
                        heap,
                        fixed.as_ref(),
                        mobile.as_ref(),
                        None,
                        None,
                    )? == 0
                    {
                        row.pINChI[TAUT_NON as usize] = SourceMutPointer::null();
                    } else {
                        let aux = rev_aux_pointer(
                            heap,
                            &structure,
                            INCHI_BAS as usize,
                            fragment,
                            TAUT_YES as usize,
                        )?;
                        row.ord_number = heap
                            .slice(aux.as_const())?
                            .first()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nNumRemovedProtons;
                    }
                }
                row.n1 = component as i16;
                row.n2 = fragment as i16;
                *heap
                    .slice_mut(sort1)?
                    .get_mut(row_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = row;
                row_index += 1;
            }
        }

        let count = usize::try_from(num_fragments)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let mut rows1 = heap
            .slice(sort1.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec();
        let mut rows2 = heap
            .slice(sort2.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec();
        let mut sort_error = None;
        rows1.sort_unstable_by(|left, right| {
            if sort_error.is_some() {
                return std::cmp::Ordering::Equal;
            }
            match CompINChITaut2(heap, left, right) {
                Ok(value) => value.cmp(&0),
                Err(error) => {
                    sort_error = Some(error);
                    std::cmp::Ordering::Equal
                }
            }
        });
        if let Some(error) = sort_error.take() {
            return Err(error);
        }
        rows2.sort_unstable_by(|left, right| {
            if sort_error.is_some() {
                return std::cmp::Ordering::Equal;
            }
            match CompINChITaut2(heap, left, right) {
                Ok(value) => value.cmp(&0),
                Err(error) => {
                    sort_error = Some(error);
                    std::cmp::Ordering::Equal
                }
            }
        });
        if let Some(error) = sort_error {
            return Err(error);
        }

        for component in 0..count {
            let _ = CompareTwoPairsOfInChI(
                heap,
                rows1[component].pINChI,
                rows2[component].pINChI,
                i32::from(bHasSomeFixedH == 0),
                &mut compare_flags,
            )?;
        }

        if bHasSomeFixedH != 0 {
            if separate_d != 0 || separate_r != 0 {
                for component in 0..count {
                    let mut isotopic1 = [0_i16; NUM_H_ISOTOPES as usize];
                    let mut isotopic2 = [0_i16; NUM_H_ISOTOPES as usize];
                    if rows1[component].ord_number != rows2[component].ord_number {
                        compare_flags[TAUT_YES as usize] |=
                            tagInchiCompareDiffBits_INCHIDIFF_MOBH_PROTONS as INCHI_MODE;
                    }
                    let input_component = i32::from(rows2[component].n1);
                    let protons = pOneInput.nNumProtons[INCHI_BAS as usize]
                        [TAUT_YES as usize]
                        .pNumProtons;
                    if !protons.is_null() {
                        isotopic2 = heap
                            .slice(protons.as_const())?
                            .get(
                                usize::try_from(input_component)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nNumRemovedIsotopicH;
                    }
                    let restored_component = i32::from(rows1[component].n1);
                    let fragment = i32::from(rows1[component].n2);
                    let mut mobile_h = input_mobile_r;
                    if bInpInchiComponentExists(
                        heap,
                        pOneInput,
                        INCHI_REC as i32,
                        mobile_h,
                        restored_component,
                    )? == 0
                    {
                        if bInpInchiComponentExists(
                            heap,
                            pOneInput,
                            INCHI_REC as i32,
                            TAUT_YES as i32,
                            restored_component,
                        )? != 0
                        {
                            mobile_h = TAUT_YES as i32;
                        } else {
                            if isotopic1 != isotopic2 {
                                compare_flags[TAUT_YES as usize] |=
                                    tagInchiCompareDiffBits_INCHIDIFF_REM_ISO_H as INCHI_MODE;
                            }
                            continue;
                        }
                    }
                    let input_pointer = component_pointer(
                        pOneInput,
                        INCHI_REC as usize,
                        mobile_h as usize,
                        restored_component,
                    )?;
                    let input_value = heap
                        .slice(input_pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    if input_value.nLink != 0 {
                        continue;
                    }
                    let structure = structure_value(
                        heap,
                        pStruct[INCHI_REC as usize][mobile_h as usize],
                        restored_component,
                    )?;
                    if fragment < structure.RevInChI.num_components[INCHI_BAS as usize] {
                        let preferred = TAUT_NON as i32;
                        if bRevInchiComponentExists(
                            heap,
                            Some(&structure),
                            INCHI_BAS as i32,
                            preferred,
                            fragment,
                        )? != 0
                            || bRevInchiComponentExists(
                                heap,
                                Some(&structure),
                                INCHI_BAS as i32,
                                TAUT_YES as i32,
                                fragment,
                            )? != 0
                        {
                            let aux = rev_aux_pointer(
                                heap,
                                &structure,
                                INCHI_BAS as usize,
                                fragment,
                                TAUT_YES as usize,
                            )?;
                            isotopic1 = heap
                                .slice(aux.as_const())?
                                .first()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .nNumRemovedIsotopicH;
                        }
                    }
                    if isotopic1 != isotopic2 {
                        compare_flags[TAUT_YES as usize] |=
                            tagInchiCompareDiffBits_INCHIDIFF_REM_ISO_H as INCHI_MODE;
                    }
                }
            }
        } else {
            if pOneInput.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize]
                .nNumRemovedProtons
                != removed_r.nNumRemovedProtons
            {
                compare_flags[TAUT_YES as usize] |=
                    tagInchiCompareDiffBits_INCHIDIFF_MOBH_PROTONS as INCHI_MODE;
            }
            if pOneInput.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize]
                .nNumRemovedIsotopicH
                != removed_r.nNumRemovedIsotopicH
            {
                compare_flags[TAUT_YES as usize] |=
                    tagInchiCompareDiffBits_INCHIDIFF_REM_ISO_H as INCHI_MODE;
            }
        }

        if i32::from(non_taut_1 == 0) != i32::from(non_taut_2 == 0) {
            if non_taut_1 != 0 {
                compare_flags[TAUT_YES as usize] |=
                    tagInchiCompareDiffBits_INCHIDIFF_WRONG_TAUT as INCHI_MODE;
            } else {
                compare_flags[TAUT_YES as usize] |=
                    tagInchiCompareDiffBits_INCHIDIFF_NO_TAUT as INCHI_MODE;
            }
        }
        for mobile_h in 0..TAUT_NUM as usize {
            pOneInput.CompareInchiFlags[1][mobile_h] |= compare_flags[mobile_h];
        }

        if removed_r.nNumRemovedProtons != removed_d.nNumRemovedProtons {
            compare_flags[TAUT_YES as usize] |=
                tagInchiCompareDiffBits_INCHIDIFF_MOBH_PROTONS as INCHI_MODE;
        }
        if removed_r.nNumRemovedIsotopicH != removed_d.nNumRemovedIsotopicH {
            compare_flags[TAUT_YES as usize] |=
                tagInchiCompareDiffBits_INCHIDIFF_REM_ISO_H as INCHI_MODE;
        }
        Ok(0)
    })();

    let free1 = inchi_free(heap, sort1);
    let free2 = inchi_free(heap, sort2);
    match (execution, free1, free2) {
        (Err(error), _, _) => Err(error),
        (Ok(_), Err(error), _) | (Ok(_), Ok(()), Err(error)) => Err(error),
        (Ok(value), Ok(()), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn CompareTwoPairsOfInChI(
    heap: &mut SourceHeap,
    pInChI1: [SourceMutPointer<INChI>; TAUT_NUM as usize],
    pInChI2: [SourceMutPointer<INChI>; TAUT_NUM as usize],
    _bMobileH: i32,
    CompareInchiFlags: &mut [INCHI_MODE; TAUT_NUM as usize],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:2203 CompareTwoPairsOfInChI
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int CompareTwoPairsOfInChI(INChI* pInChI1[TAUT_NUM],
    INChI* pInChI2[TAUT_NUM],
    int bMobileH,
    INCHI_MODE CompareInchiFlags[])
{
    int iMobileH, err = 0;
    INCHI_MODE cmp;
    for (iMobileH = 0; iMobileH < TAUT_NUM; iMobileH++)
    {
        if (!pInChI1[iMobileH] != !pInChI2[iMobileH])
        {
            if (iMobileH == TAUT_NON &&
                pInChI1[TAUT_YES] && pInChI2[TAUT_YES]) /* djb-rwth: condition corrected */
            {
                CompareInchiFlags[iMobileH] |= INCHIDIFF_COMP_HLAYER;
            }
            else
            {
                CompareInchiFlags[iMobileH] |= INCHIDIFF_COMP_NUMBER;
            }
            continue;
        }
        if (pInChI1[iMobileH] && pInChI2[iMobileH])
        {
            cmp = CompareReversedINChI3(pInChI1[iMobileH], pInChI2[iMobileH], NULL, NULL, &err);
            if (cmp)
            {
                CompareInchiFlags[iMobileH] |= cmp;
            }
        }
    }

    return err;
}
    */
    // END INCHI C FUNCTION: CompareTwoPairsOfInChI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CompareTwoPairsOfInChI
    // INCHI✔️❌: READ_INCHI_STRING=1 includes this production helper.
    // INCHI✔️❌: COMPILE_ANSI_ONLY and TARGET_API_LIB do not alter this function body.
    // INCHI✔️❌: Shallow INChI clones avoid retaining SourceHeap borrows across the mutable callee.
    // END INCHI ACTIVE MACRO CONFIGURATION: CompareTwoPairsOfInChI

    let mut err = 0_i32;
    for iMobileH in 0..TAUT_NUM as usize {
        if pInChI1[iMobileH].is_null() != pInChI2[iMobileH].is_null() {
            if iMobileH == TAUT_NON as usize
                && !pInChI1[TAUT_YES as usize].is_null()
                && !pInChI2[TAUT_YES as usize].is_null()
            {
                CompareInchiFlags[iMobileH] |=
                    tagInchiCompareDiffBits_INCHIDIFF_COMP_HLAYER as INCHI_MODE;
            } else {
                CompareInchiFlags[iMobileH] |=
                    tagInchiCompareDiffBits_INCHIDIFF_COMP_NUMBER as INCHI_MODE;
            }
            continue;
        }
        if !pInChI1[iMobileH].is_null() && !pInChI2[iMobileH].is_null() {
            let inchi1 = heap
                .slice(pInChI1[iMobileH].as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let inchi2 = heap
                .slice(pInChI2[iMobileH].as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let cmp = CompareReversedINChI3(
                heap,
                Some(&inchi1),
                Some(&inchi2),
                None,
                None,
                &mut err,
            )?;
            if cmp != 0 {
                CompareInchiFlags[iMobileH] |= cmp;
            }
        }
    }
    Ok(err)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn CompareOneOrigInchiToRevInChI(
    heap: &mut SourceHeap,
    pStruct: Option<&StrFromINChI>,
    pInChI: [SourceMutPointer<INChI>; TAUT_NUM as usize],
    bMobileH: i32,
    _iComponent: i32,
    _num_inp: i64,
    _szCurHdr: SourceMutPointer<i8>,
    nCurRemovedProtons: &mut COMPONENT_REM_PROTONS,
    CompareInchiFlags: &mut [INCHI_MODE; TAUT_NUM as usize],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:2240 CompareOneOrigInchiToRevInChI
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int CompareOneOrigInchiToRevInChI(StrFromINChI* pStruct,
    INChI* pInChI[TAUT_NUM],
    int bMobileH,
    int iComponent,
    long num_inp,
    char* szCurHdr,
    COMPONENT_REM_PROTONS* nCurRemovedProtons,
    INCHI_MODE CompareInchiFlags[])
{
    int ret, err = 0;
    INCHI_MODE cmp;

    if (pStruct) /* djb-rwth: fixing coverity ID #499601 */
    {
        ret = pStruct->RevInChI.nRetVal;
    }
    else
    {
        ret = -1;
    }

    if ((ret == _IS_OKAY || ret == _IS_WARNING) && pStruct) /* djb-rwth: fixing a NULL pointer dereference */
    {
        /* ignore bMobileH for now */
        int i, i0, b /* created type */, b0 /* requested type*/, j, k;
        /* pINChI[iINCHI][iComponent][bTaut] */
        /* i0 = requested Rec/Disconnected: 1/0 */
        /* i  = what InChI creaded out of the restored structure */
        /* b0 = requested Mobile/Fixed-H: 1/0 */
        /* b  = what InChI creaded out of the restored structure */
        i = i0 = pStruct->iINCHI;
        b = b0 = pStruct->iMobileH;
        if (i == INCHI_REC && !pStruct->RevInChI.num_components[i])
        {
            i = INCHI_BAS;
        }
        if (b == TAUT_NON && (!pStruct->RevInChI.pINChI[i] ||
            !pStruct->RevInChI.pINChI[i][0][b] ||
            !pStruct->RevInChI.pINChI[i][0][b]->nNumberOfAtoms))
        {
            b = TAUT_YES;
        }
        if (pStruct->bDeleted && (!pInChI[0] || pInChI[0]->bDeleted))
        {
            return 0;
        }

        if (pStruct->RevInChI.pINChI[i]) /* djb-rwth: fixing a NULL pointer dereference */
        {
            if ((pStruct->RevInChI.num_components[i] > 1 &&
                !pStruct->RevInChI.pINChI[i][1][b]->bDeleted) ||
                pStruct->RevInChI.num_components[i] < 1) /* djb-rwth: addressing LLVM warning */
            {
                CompareInchiFlags[bMobileH] |= INCHIDIFF_COMP_NUMBER;
            }
        }
        if (b != b0 || b != bMobileH || b0 != bMobileH || i > i0)
        {
            /* do not print messages about TAUT_YES instead of TAUT_NON */
            CompareInchiFlags[bMobileH] |= INCHIDIFF_COMP_HLAYER;
        }

        if (pStruct->RevInChI.num_components[i] && pStruct->RevInChI.pINChI[i]) /* djb-rwth: fixing a NULL pointer dereference */
        {
            /* compare InChI from restored structure; '0' in [i][0][b] is the first component */
            if (b == TAUT_YES && pStruct->RevInChI.pINChI[i][0][b]->bDeleted && ( !pInChI[0] || pInChI[0]->bDeleted ))
            {
                /* the 1st component is made out of proton(s) and the input component is missing or also a proton */
                cmp = 0;
            }
            else
            {
                cmp = CompareReversedINChI3( pStruct->RevInChI.pINChI[i][0][b], pInChI[0], NULL, NULL, &err );
                if (cmp)
                {
                    CompareInchiFlags[bMobileH] |= cmp;
                }
            }
            if (b == b0 && b == TAUT_NON)
            {
                if ((pStruct->RevInChI.pINChI[i][0][TAUT_YES] &&
                    !pStruct->RevInChI.pINChI[i][0][TAUT_YES]->bDeleted) ||
                    (pInChI[1] && !pInChI[1]->bDeleted)) /* djb-rwth: addressing LLVM warnings */
                {
                    /* in addition to fixed-H also compare mobile-H InChI */
                    cmp = CompareReversedINChI3(pStruct->RevInChI.pINChI[i][0][TAUT_YES], pInChI[1], NULL, NULL, &err);
                    if (cmp)
                    {
                        CompareInchiFlags[TAUT_YES] |= cmp;
                    }
                }
                /* compare removed H */
                if (pStruct->nNumRemovedProtonsMobHInChI != pStruct->RevInChI.pINChI_Aux[i][0][TAUT_YES]->nNumRemovedProtons)
                {
                    CompareInchiFlags[TAUT_YES] |= INCHIDIFF_MOBH_PROTONS;
                }
            }
            memset(nCurRemovedProtons, 0, sizeof(*nCurRemovedProtons)); /* djb-rwth: memset_s C11/Annex K variant? */
            for (k = 0; k < pStruct->RevInChI.num_components[i]; k++)
            {
                if (!k || pStruct->RevInChI.pINChI[i][k][TAUT_YES]->bDeleted)
                {
                    /* get removed protons from the 1st component; add othere only if they are deleted protons */
                    nCurRemovedProtons->nNumRemovedProtons += pStruct->RevInChI.pINChI_Aux[i][k][TAUT_YES]->nNumRemovedProtons;
                    for (j = 0; j < NUM_H_ISOTOPES; j++)
                    {
                        nCurRemovedProtons->nNumRemovedIsotopicH[j] += pStruct->RevInChI.pINChI_Aux[i][k][TAUT_YES]->nNumRemovedIsotopicH[j];
                    }
                }
            }
        }
    }
    else
    {
        CompareInchiFlags[bMobileH] |= INCHIDIFF_STR2INCHI_ERR;
    }

    return err;
}
    */
    // END INCHI C FUNCTION: CompareOneOrigInchiToRevInChI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CompareOneOrigInchiToRevInChI
    // INCHI✔️❌: READ_INCHI_STRING=1 includes this production helper.
    // INCHI✔️❌: COMPILE_ANSI_ONLY and TARGET_API_LIB do not alter this function body.
    // INCHI✔️❌: Shallow structure clones avoid holding SourceHeap borrows across mutable callee calls.
    // END INCHI ACTIVE MACRO CONFIGURATION: CompareOneOrigInchiToRevInChI

    let mobile_index = usize::try_from(bMobileH)
        .ok()
        .filter(|index| *index < TAUT_NUM as usize)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let input = pInChI.map(|pointer| -> Result<Option<INChI>, SourceHeapError> {
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
    });
    let [input0, input1] = input;
    let input0 = input0?;
    let input1 = input1?;
    let ret = pStruct.map_or(-1, |structure| structure.RevInChI.nRetVal);
    let mut err = 0_i32;
    if (ret == _IS_OKAY as i32 || ret == _IS_WARNING as i32) && pStruct.is_some() {
        let structure = pStruct.ok_or(SourceHeapError::NullPointer)?;
        let i0 = i32::from(structure.iINCHI);
        let b0 = i32::from(structure.iMobileH);
        let mut i = i0;
        let mut b = b0;
        let i0_index = usize::try_from(i0)
            .ok()
            .filter(|index| *index < INCHI_NUM as usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if i == INCHI_REC as i32 && structure.RevInChI.num_components[i0_index] == 0 {
            i = INCHI_BAS as i32;
        }
        let i_index = usize::try_from(i)
            .ok()
            .filter(|index| *index < INCHI_NUM as usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let b0_index = usize::try_from(b0)
            .ok()
            .filter(|index| *index < TAUT_NUM as usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let rows_pointer = structure.RevInChI.pINChI[i_index];
        if b == TAUT_NON as i32 {
            let fixed_missing = if rows_pointer.is_null() {
                true
            } else {
                let rows = heap.slice(rows_pointer.as_const())?;
                let row = rows.first().ok_or(SourceHeapError::PointerOutOfBounds)?;
                let pointer = row[b0_index];
                if pointer.is_null() {
                    true
                } else {
                    heap.slice(pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nNumberOfAtoms
                        == 0
                }
            };
            if fixed_missing {
                b = TAUT_YES as i32;
            }
        }
        let b_index = usize::try_from(b)
            .ok()
            .filter(|index| *index < TAUT_NUM as usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if structure.bDeleted != 0 && input0.as_ref().is_none_or(|inchi| inchi.bDeleted != 0) {
            return Ok(0);
        }
        if !rows_pointer.is_null() {
            let rows = heap.slice(rows_pointer.as_const())?;
            if structure.RevInChI.num_components[i_index] > 1 {
                let row = rows.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?;
                let pointer = row[b_index];
                let second = if pointer.is_null() {
                    return Err(SourceHeapError::NullPointer);
                } else {
                    heap.slice(pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                };
                if second.bDeleted == 0 {
                    CompareInchiFlags[mobile_index] |=
                        tagInchiCompareDiffBits_INCHIDIFF_COMP_NUMBER as INCHI_MODE;
                }
            } else if structure.RevInChI.num_components[i_index] < 1 {
                CompareInchiFlags[mobile_index] |=
                    tagInchiCompareDiffBits_INCHIDIFF_COMP_NUMBER as INCHI_MODE;
            }
        }
        if b != b0 || b != bMobileH || b0 != bMobileH || i > i0 {
            CompareInchiFlags[mobile_index] |=
                tagInchiCompareDiffBits_INCHIDIFF_COMP_HLAYER as INCHI_MODE;
        }
        if structure.RevInChI.num_components[i_index] != 0 && !rows_pointer.is_null() {
            let first_row = *heap
                .slice(rows_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let generated_pointer = first_row[b_index];
            let generated = if generated_pointer.is_null() {
                return Err(SourceHeapError::NullPointer);
            } else {
                heap.slice(generated_pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone()
            };
            let skip_deleted_mobile = b == TAUT_YES as i32
                && generated.bDeleted != 0
                && input0.as_ref().is_none_or(|inchi| inchi.bDeleted != 0);
            if !skip_deleted_mobile {
                let cmp = CompareReversedINChI3(
                    heap,
                    Some(&generated),
                    input0.as_ref(),
                    None,
                    None,
                    &mut err,
                )?;
                if cmp != 0 {
                    CompareInchiFlags[mobile_index] |= cmp;
                }
            }
            if b == b0 && b == TAUT_NON as i32 {
                let generated_mobile_pointer = first_row[TAUT_YES as usize];
                let generated_mobile = if generated_mobile_pointer.is_null() {
                    return Err(SourceHeapError::NullPointer);
                } else {
                    heap.slice(generated_mobile_pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone()
                };
                if generated_mobile.bDeleted == 0
                    || input1.as_ref().is_some_and(|inchi| inchi.bDeleted == 0)
                {
                    let cmp = CompareReversedINChI3(
                        heap,
                        Some(&generated_mobile),
                        input1.as_ref(),
                        None,
                        None,
                        &mut err,
                    )?;
                    if cmp != 0 {
                        CompareInchiFlags[TAUT_YES as usize] |= cmp;
                    }
                }
                let aux_rows_pointer = structure.RevInChI.pINChI_Aux[i_index];
                if aux_rows_pointer.is_null() {
                    return Err(SourceHeapError::NullPointer);
                }
                let aux_row = *heap
                    .slice(aux_rows_pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let aux_pointer = aux_row[TAUT_YES as usize];
                let aux = if aux_pointer.is_null() {
                    return Err(SourceHeapError::NullPointer);
                } else {
                    heap.slice(aux_pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                };
                if structure.nNumRemovedProtonsMobHInChI != i32::from(aux.nNumRemovedProtons) {
                    CompareInchiFlags[TAUT_YES as usize] |=
                        tagInchiCompareDiffBits_INCHIDIFF_MOBH_PROTONS as INCHI_MODE;
                }
            }
            *nCurRemovedProtons = COMPONENT_REM_PROTONS::default();
            let component_count = structure.RevInChI.num_components[i_index];
            if component_count > 0 {
                let count = usize::try_from(component_count)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let inchi_rows = heap
                    .slice(rows_pointer.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                let aux_rows_pointer = structure.RevInChI.pINChI_Aux[i_index];
                if aux_rows_pointer.is_null() {
                    return Err(SourceHeapError::NullPointer);
                }
                let aux_rows = heap
                    .slice(aux_rows_pointer.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                for k in 0..count {
                    let mobile_pointer = inchi_rows[k][TAUT_YES as usize];
                    let mobile = if mobile_pointer.is_null() {
                        return Err(SourceHeapError::NullPointer);
                    } else {
                        heap.slice(mobile_pointer.as_const())?
                            .first()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                    };
                    if k == 0 || mobile.bDeleted != 0 {
                        let aux_pointer = aux_rows[k][TAUT_YES as usize];
                        let aux = if aux_pointer.is_null() {
                            return Err(SourceHeapError::NullPointer);
                        } else {
                            heap.slice(aux_pointer.as_const())?
                                .first()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                        };
                        nCurRemovedProtons.nNumRemovedProtons = nCurRemovedProtons
                            .nNumRemovedProtons
                            .wrapping_add(aux.nNumRemovedProtons);
                        for isotope in 0..NUM_H_ISOTOPES as usize {
                            nCurRemovedProtons.nNumRemovedIsotopicH[isotope] =
                                nCurRemovedProtons.nNumRemovedIsotopicH[isotope]
                                    .wrapping_add(aux.nNumRemovedIsotopicH[isotope]);
                        }
                    }
                }
            }
        }
    } else {
        CompareInchiFlags[mobile_index] |=
            tagInchiCompareDiffBits_INCHIDIFF_STR2INCHI_ERR as INCHI_MODE;
    }
    Ok(err)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn CompareReversedStereoINChI3(
    heap: &SourceHeap,
    s1: Option<&INChI_Stereo>,
    s2: Option<&INChI_Stereo>,
    picr: &mut ICR,
) -> Result<INCHI_MODE, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:2362 CompareReversedStereoINChI3
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
INCHI_MODE CompareReversedStereoINChI3(INChI_Stereo* s1,
    /* InChI from reversed struct   */
    INChI_Stereo* s2,
    /* input InChI                  */
    ICR* picr)
{
    int ret = 0;
    int j1, j2, num_dif, num_extra_undf, num_miss_undf, num_in1_only, num_in2_only; /* djb-rwth: removing redundant variables */
    int bAddSb = !(picr->num_sb_undef_in1_only + picr->num_sb_in1_only + picr->num_sb_in2_only);
    int bAddSc = !(picr->num_sc_undef_in1_only + picr->num_sc_in1_only + picr->num_sc_in2_only);

    int nNumSc1 = s1 ? s1->nNumberOfStereoCenters : 0;
    int nNumSc2 = s2 ? s2->nNumberOfStereoCenters : 0;
    int nNumSb1 = s1 ? s1->nNumberOfStereoBonds : 0;
    int nNumSb2 = s2 ? s2->nNumberOfStereoBonds : 0;

    if ((nNumSc1 || nNumSc2) &&
        (nNumSc1 != nNumSc2 ||
            memcmp(s1->nNumber, s2->nNumber, nNumSc1 * sizeof(s1->nNumber[0])) ||
            memcmp(s1->t_parity, s2->t_parity, nNumSc1 * sizeof(s1->t_parity[0])))) /* djb-rwth: condition corrected */
    {

        num_dif = num_extra_undf = num_miss_undf = num_in1_only = num_in2_only = 0; /* djb-rwth: removing redundant code */
        for (j1 = j2 = 0; j1 < nNumSc1 && j2 < nNumSc2; )
        {
            if (s1->nNumber[j1] == s2->nNumber[j2])
            {
                if (s1->t_parity[j1] != s2->t_parity[j2])
                    /* djb-rwth: removing redundant code */
                {
                    num_dif++;
                }
                j1++;
                j2++;
            }
            else
            {
                if (s1->nNumber[j1] < s2->nNumber[j2])
                {
                    num_in1_only++;
                    if (s1->t_parity[j1] == AB_PARITY_UNDF)
                    {
                        num_extra_undf++;
                    }
                    if (bAddSc)
                    {
                        if (picr->num_sc_in1_only < ICR_MAX_SC_IN1_ONLY)
                            picr->sc_in1_only[picr->num_sc_in1_only++] = j1;
                        if (s1->t_parity[j1] == AB_PARITY_UNDF)
                        {
                            if (picr->num_sc_undef_in1_only < ICR_MAX_SC_UNDF)
                                picr->sc_undef_in1_only[picr->num_sc_undef_in1_only++] = j1;
                        }
                    }
                    j1++;
                }
                else
                {
                    num_in2_only++;
                    if (s2->t_parity[j2] == AB_PARITY_UNDF)
                    {
                        num_miss_undf++;
                    }
                    if (bAddSc)
                    {
                        if (picr->num_sc_in2_only < ICR_MAX_SC_IN2_ONLY)
                            picr->sc_in2_only[picr->num_sc_in2_only++] = j2;
                        if (s2->t_parity[j2] == AB_PARITY_UNDF)
                        {
                            if (picr->num_sc_undef_in2_only < ICR_MAX_SC_UNDF)
                                picr->sc_undef_in2_only[picr->num_sc_undef_in2_only++] = j1;
                        }
                    }
                    j2++;
                }
            }
        }

        while (j1 < nNumSc1)
        {
            if (s1->t_parity[j1] == AB_PARITY_UNDF)
            {
                num_extra_undf++;
            }
            num_in1_only++;
            if (bAddSc)
            {
                if (picr->num_sc_in1_only < ICR_MAX_SC_IN1_ONLY)
                    picr->sc_in1_only[picr->num_sc_in1_only++] = j1;
                if (s1->t_parity[j1] == AB_PARITY_UNDF)
                {
                    if (picr->num_sc_undef_in1_only < ICR_MAX_SC_UNDF)
                        picr->sc_undef_in1_only[picr->num_sc_undef_in1_only++] = j1;
                }
            }
            j1++;
        }

        while (j2 < nNumSc2)
        {
            if (s2->t_parity[j2] == AB_PARITY_UNDF)
            {
                num_miss_undf++;
            }
            num_in2_only++;
            if (bAddSc)
            {
                if (picr->num_sc_in2_only < ICR_MAX_SC_IN2_ONLY)
                    picr->sc_in2_only[picr->num_sc_in2_only++] = j2;
            }
            j2++;
        }

        if (num_dif)
        {
            ret |= INCHIDIFF_SC_PARITY;
        }
        if (num_in1_only)
        {
            if (num_extra_undf)
            {
                ret |= INCHIDIFF_SC_EXTRA_UNDF;
            }
            if (num_in1_only != num_extra_undf)
            {
                ret |= INCHIDIFF_SC_EXTRA;
            }
        }
        if (num_in2_only)
        {
            if (num_miss_undf)
            {
                ret |= INCHIDIFF_SC_MISS_UNDF;
            }
            if (num_in2_only != num_miss_undf)
            {
                ret |= INCHIDIFF_SC_MISS;
            }
        }
    }

    if (s1 && s2 && (s2->nCompInv2Abs != 2) && s1->nCompInv2Abs != s2->nCompInv2Abs && s1->nCompInv2Abs && s2->nCompInv2Abs)
    {
        ret |= INCHIDIFF_SC_INV; /* 2007-07-13 DT: added (s2->nCompInv2Abs != 2) to fix bug reoprted by Yerin on 2007/02/28 */
        /* Bug description: falsely reported "Stereo centers/allenes: Falsely inverted" for /S2 or /S3 */
    }

    if ((nNumSb1 || nNumSb2) &&
        (nNumSb1 != nNumSb2 ||
            memcmp(s1->nBondAtom1, s2->nBondAtom1, nNumSb1 * sizeof(s1->nBondAtom1[0])) ||
            memcmp(s1->nBondAtom2, s2->nBondAtom2, nNumSb1 * sizeof(s1->nBondAtom2[0])) ||
            memcmp(s1->b_parity, s2->b_parity, nNumSb1 * sizeof(s1->b_parity[0]))))
    {

        num_dif = num_extra_undf = num_miss_undf = num_in1_only = num_in2_only = 0; /* djb-rwth: removing redundant code */
        for (j1 = j2 = 0; j1 < nNumSb1 && j2 < nNumSb2; )
        {
            if (s1->nBondAtom1[j1] == s2->nBondAtom1[j2] &&
                s1->nBondAtom2[j1] == s2->nBondAtom2[j2])
            {
                if (s1->b_parity[j1] != s2->b_parity[j2])
                    /* djb-rwth: removing redundant code */
                {
                    num_dif++;
                }
                j1++;
                j2++;
            }
            else
            {
                if (s1->nBondAtom1[j1] < s2->nBondAtom1[j2] ||
                    (s1->nBondAtom1[j1] == s2->nBondAtom1[j2] && s1->nBondAtom2[j1] < s2->nBondAtom2[j2])) /* djb-rwth: addressing LLVM warning */
                {
                    num_in1_only++;
                    if (s1->b_parity[j1] == AB_PARITY_UNDF)
                    {
                        num_extra_undf++;
                    }
                    if (bAddSb)
                    {
                        if (picr->num_sb_in1_only < ICR_MAX_SB_IN1_ONLY)
                            picr->sb_in1_only[picr->num_sb_in1_only++] = j1;
                        if (s1->b_parity[j1] == AB_PARITY_UNDF)
                        {
                            if (picr->num_sb_undef_in1_only < ICR_MAX_SB_UNDF)
                                picr->sb_undef_in1_only[picr->num_sb_undef_in1_only++] = j1;
                        }
                    }
                    j1++;
                }
                else
                {
                    num_in2_only++;
                    if (s2->b_parity[j2] == AB_PARITY_UNDF)
                    {
                        num_miss_undf++;
                    }
                    if (bAddSb)
                    {
                        if (picr->num_sb_in2_only < ICR_MAX_SB_IN2_ONLY)
                            picr->sb_in2_only[picr->num_sb_in2_only++] = j2;
                        if (s2->b_parity[j2] == AB_PARITY_UNDF)
                        {
                            if (picr->num_sb_undef_in2_only < ICR_MAX_SB_UNDF)
                                picr->sb_undef_in2_only[picr->num_sb_undef_in2_only++] = j1;
                        }
                    }
                    j2++;
                }
            }
        }
        while (j1 < nNumSb1)
        {
            num_in1_only++;
            if (s1->b_parity[j1] == AB_PARITY_UNDF)
            {
                num_extra_undf++;
            }
            if (bAddSb)
            {
                if (picr->num_sb_in1_only < ICR_MAX_SB_IN1_ONLY)
                    picr->sb_in1_only[picr->num_sb_in1_only++] = j1;
                if (s1->b_parity[j1] == AB_PARITY_UNDF)
                {
                    if (picr->num_sb_undef_in1_only < ICR_MAX_SB_UNDF)
                        picr->sb_undef_in1_only[picr->num_sb_undef_in1_only++] = j1;
                }
            }
            j1++;
        }
        while (j2 < nNumSb2)
        {
            num_in2_only++;
            if (s2->b_parity[j2] == AB_PARITY_UNDF)
            {
                num_miss_undf++;
            }
            if (bAddSb)
            {
                if (picr->num_sb_in2_only < ICR_MAX_SB_IN2_ONLY)
                    picr->sb_in2_only[picr->num_sb_in2_only++] = j2;
                if (s2->b_parity[j2] == AB_PARITY_UNDF)
                {
                    if (picr->num_sb_undef_in2_only < ICR_MAX_SB_UNDF)
                        picr->sb_undef_in2_only[picr->num_sb_undef_in2_only++] = j1;
                }
            }
            j2++;
        }
        if (num_dif)
        {
            ret |= INCHIDIFF_SB_PARITY;
        }
        if (num_in1_only)
        {
            if (num_extra_undf)
            {
                ret |= INCHIDIFF_SB_EXTRA_UNDF;
            }
            if (num_in1_only != num_extra_undf)
            {
                ret |= INCHIDIFF_SB_EXTRA;
            }
        }
        if (num_in2_only)
        {
            if (num_miss_undf)
            {
                ret |= INCHIDIFF_SB_MISS_UNDF;
            }
            if (num_in2_only != num_miss_undf)
            {
                ret |= INCHIDIFF_SB_MISS;
            }
        }
    }

    return ret;
}
    */
    // END INCHI C FUNCTION: CompareReversedStereoINChI3
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CompareReversedStereoINChI3
    // INCHI✔️❌: READ_INCHI_STRING=1 includes this production helper.
    // INCHI✔️❌: COMPILE_ANSI_ONLY and TARGET_API_LIB do not alter this function body.
    // INCHI✔️❌: Checked SourceHeap slices add overhead versus direct memcmp and pointer indexing.
    // END INCHI ACTIVE MACRO CONFIGURATION: CompareReversedStereoINChI3

    let add_sb = picr.num_sb_undef_in1_only
        .wrapping_add(picr.num_sb_in1_only)
        .wrapping_add(picr.num_sb_in2_only) == 0;
    let add_sc = picr.num_sc_undef_in1_only
        .wrapping_add(picr.num_sc_in1_only)
        .wrapping_add(picr.num_sc_in2_only) == 0;
    let num_sc1 = s1.map_or(0, |stereo| stereo.nNumberOfStereoCenters);
    let num_sc2 = s2.map_or(0, |stereo| stereo.nNumberOfStereoCenters);
    let num_sb1 = s1.map_or(0, |stereo| stereo.nNumberOfStereoBonds);
    let num_sb2 = s2.map_or(0, |stereo| stereo.nNumberOfStereoBonds);
    let mut ret: INCHI_MODE = 0;

    if num_sc1 != 0 || num_sc2 != 0 {
        let count1 = usize::try_from(num_sc1)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let count2 = usize::try_from(num_sc2)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let (numbers1, parities1): (&[AT_NUMB], &[i8]) = if count1 == 0 {
            (&[], &[])
        } else {
            let stereo = s1.ok_or(SourceHeapError::NullPointer)?;
            (
                heap.slice(stereo.nNumber.as_const())?
                    .get(..count1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                heap.slice(stereo.t_parity.as_const())?
                    .get(..count1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
        };
        let (numbers2, parities2): (&[AT_NUMB], &[i8]) = if count2 == 0 {
            (&[], &[])
        } else {
            let stereo = s2.ok_or(SourceHeapError::NullPointer)?;
            (
                heap.slice(stereo.nNumber.as_const())?
                    .get(..count2)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                heap.slice(stereo.t_parity.as_const())?
                    .get(..count2)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
        };
        if num_sc1 != num_sc2 || numbers1 != numbers2 || parities1 != parities2 {
            let (mut j1, mut j2) = (0_usize, 0_usize);
            let (mut num_dif, mut num_extra_undf, mut num_miss_undf) = (0_i32, 0_i32, 0_i32);
            let (mut num_in1_only, mut num_in2_only) = (0_i32, 0_i32);
            while j1 < count1 && j2 < count2 {
                if numbers1[j1] == numbers2[j2] {
                    if parities1[j1] != parities2[j2] {
                        num_dif = num_dif.wrapping_add(1);
                    }
                    j1 += 1;
                    j2 += 1;
                } else if numbers1[j1] < numbers2[j2] {
                    num_in1_only = num_in1_only.wrapping_add(1);
                    let undefined = parities1[j1] == AB_PARITY_UNDF as i8;
                    if undefined {
                        num_extra_undf = num_extra_undf.wrapping_add(1);
                    }
                    if add_sc {
                        if picr.num_sc_in1_only < 0 || picr.num_sc_undef_in1_only < 0 {
                            return Err(SourceHeapError::SourceIntegerOverflow);
                        }
                        if picr.num_sc_in1_only < ICR_MAX_SC_IN1_ONLY as i32 {
                            picr.sc_in1_only[picr.num_sc_in1_only as usize] = j1 as AT_NUMB;
                            picr.num_sc_in1_only = picr.num_sc_in1_only.wrapping_add(1);
                        }
                        if undefined && picr.num_sc_undef_in1_only < ICR_MAX_SC_UNDF as i32 {
                            picr.sc_undef_in1_only[picr.num_sc_undef_in1_only as usize] = j1 as AT_NUMB;
                            picr.num_sc_undef_in1_only = picr.num_sc_undef_in1_only.wrapping_add(1);
                        }
                    }
                    j1 += 1;
                } else {
                    num_in2_only = num_in2_only.wrapping_add(1);
                    let undefined = parities2[j2] == AB_PARITY_UNDF as i8;
                    if undefined {
                        num_miss_undf = num_miss_undf.wrapping_add(1);
                    }
                    if add_sc {
                        if picr.num_sc_in2_only < 0 || picr.num_sc_undef_in2_only < 0 {
                            return Err(SourceHeapError::SourceIntegerOverflow);
                        }
                        if picr.num_sc_in2_only < ICR_MAX_SC_IN2_ONLY as i32 {
                            picr.sc_in2_only[picr.num_sc_in2_only as usize] = j2 as AT_NUMB;
                            picr.num_sc_in2_only = picr.num_sc_in2_only.wrapping_add(1);
                        }
                        if undefined && picr.num_sc_undef_in2_only < ICR_MAX_SC_UNDF as i32 {
                            picr.sc_undef_in2_only[picr.num_sc_undef_in2_only as usize] = j1 as AT_NUMB;
                            picr.num_sc_undef_in2_only = picr.num_sc_undef_in2_only.wrapping_add(1);
                        }
                    }
                    j2 += 1;
                }
            }
            while j1 < count1 {
                let undefined = parities1[j1] == AB_PARITY_UNDF as i8;
                if undefined {
                    num_extra_undf = num_extra_undf.wrapping_add(1);
                }
                num_in1_only = num_in1_only.wrapping_add(1);
                if add_sc {
                    if picr.num_sc_in1_only < 0 || picr.num_sc_undef_in1_only < 0 {
                        return Err(SourceHeapError::SourceIntegerOverflow);
                    }
                    if picr.num_sc_in1_only < ICR_MAX_SC_IN1_ONLY as i32 {
                        picr.sc_in1_only[picr.num_sc_in1_only as usize] = j1 as AT_NUMB;
                        picr.num_sc_in1_only = picr.num_sc_in1_only.wrapping_add(1);
                    }
                    if undefined && picr.num_sc_undef_in1_only < ICR_MAX_SC_UNDF as i32 {
                        picr.sc_undef_in1_only[picr.num_sc_undef_in1_only as usize] = j1 as AT_NUMB;
                        picr.num_sc_undef_in1_only = picr.num_sc_undef_in1_only.wrapping_add(1);
                    }
                }
                j1 += 1;
            }
            while j2 < count2 {
                let undefined = parities2[j2] == AB_PARITY_UNDF as i8;
                if undefined {
                    num_miss_undf = num_miss_undf.wrapping_add(1);
                }
                num_in2_only = num_in2_only.wrapping_add(1);
                if add_sc {
                    if picr.num_sc_in2_only < 0 {
                        return Err(SourceHeapError::SourceIntegerOverflow);
                    }
                    if picr.num_sc_in2_only < ICR_MAX_SC_IN2_ONLY as i32 {
                        picr.sc_in2_only[picr.num_sc_in2_only as usize] = j2 as AT_NUMB;
                        picr.num_sc_in2_only = picr.num_sc_in2_only.wrapping_add(1);
                    }
                }
                j2 += 1;
            }
            if num_dif != 0 {
                ret |= tagInchiCompareDiffBits_INCHIDIFF_SC_PARITY as INCHI_MODE;
            }
            if num_in1_only != 0 {
                if num_extra_undf != 0 {
                    ret |= tagInchiCompareDiffBits_INCHIDIFF_SC_EXTRA_UNDF as INCHI_MODE;
                }
                if num_in1_only != num_extra_undf {
                    ret |= tagInchiCompareDiffBits_INCHIDIFF_SC_EXTRA as INCHI_MODE;
                }
            }
            if num_in2_only != 0 {
                if num_miss_undf != 0 {
                    ret |= tagInchiCompareDiffBits_INCHIDIFF_SC_MISS_UNDF as INCHI_MODE;
                }
                if num_in2_only != num_miss_undf {
                    ret |= tagInchiCompareDiffBits_INCHIDIFF_SC_MISS as INCHI_MODE;
                }
            }
        }
    }
    if let (Some(stereo1), Some(stereo2)) = (s1, s2) {
        if stereo2.nCompInv2Abs != 2
            && stereo1.nCompInv2Abs != stereo2.nCompInv2Abs
            && stereo1.nCompInv2Abs != 0
            && stereo2.nCompInv2Abs != 0
        {
            ret |= tagInchiCompareDiffBits_INCHIDIFF_SC_INV as INCHI_MODE;
        }
    }

    if num_sb1 != 0 || num_sb2 != 0 {
        let count1 = usize::try_from(num_sb1)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let count2 = usize::try_from(num_sb2)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let (atom11, atom12, parity1): (&[AT_NUMB], &[AT_NUMB], &[i8]) = if count1 == 0 {
            (&[], &[], &[])
        } else {
            let stereo = s1.ok_or(SourceHeapError::NullPointer)?;
            (
                heap.slice(stereo.nBondAtom1.as_const())?
                    .get(..count1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                heap.slice(stereo.nBondAtom2.as_const())?
                    .get(..count1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                heap.slice(stereo.b_parity.as_const())?
                    .get(..count1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
        };
        let (atom21, atom22, parity2): (&[AT_NUMB], &[AT_NUMB], &[i8]) = if count2 == 0 {
            (&[], &[], &[])
        } else {
            let stereo = s2.ok_or(SourceHeapError::NullPointer)?;
            (
                heap.slice(stereo.nBondAtom1.as_const())?
                    .get(..count2)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                heap.slice(stereo.nBondAtom2.as_const())?
                    .get(..count2)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                heap.slice(stereo.b_parity.as_const())?
                    .get(..count2)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
        };
        if num_sb1 != num_sb2 || atom11 != atom21 || atom12 != atom22 || parity1 != parity2 {
            let (mut j1, mut j2) = (0_usize, 0_usize);
            let (mut num_dif, mut num_extra_undf, mut num_miss_undf) = (0_i32, 0_i32, 0_i32);
            let (mut num_in1_only, mut num_in2_only) = (0_i32, 0_i32);
            while j1 < count1 && j2 < count2 {
                if atom11[j1] == atom21[j2] && atom12[j1] == atom22[j2] {
                    if parity1[j1] != parity2[j2] {
                        num_dif = num_dif.wrapping_add(1);
                    }
                    j1 += 1;
                    j2 += 1;
                } else if atom11[j1] < atom21[j2]
                    || (atom11[j1] == atom21[j2] && atom12[j1] < atom22[j2])
                {
                    num_in1_only = num_in1_only.wrapping_add(1);
                    let undefined = parity1[j1] == AB_PARITY_UNDF as i8;
                    if undefined {
                        num_extra_undf = num_extra_undf.wrapping_add(1);
                    }
                    if add_sb {
                        if picr.num_sb_in1_only < 0 || picr.num_sb_undef_in1_only < 0 {
                            return Err(SourceHeapError::SourceIntegerOverflow);
                        }
                        if picr.num_sb_in1_only < ICR_MAX_SB_IN1_ONLY as i32 {
                            picr.sb_in1_only[picr.num_sb_in1_only as usize] = j1 as AT_NUMB;
                            picr.num_sb_in1_only = picr.num_sb_in1_only.wrapping_add(1);
                        }
                        if undefined && picr.num_sb_undef_in1_only < ICR_MAX_SB_UNDF as i32 {
                            picr.sb_undef_in1_only[picr.num_sb_undef_in1_only as usize] = j1 as AT_NUMB;
                            picr.num_sb_undef_in1_only = picr.num_sb_undef_in1_only.wrapping_add(1);
                        }
                    }
                    j1 += 1;
                } else {
                    num_in2_only = num_in2_only.wrapping_add(1);
                    let undefined = parity2[j2] == AB_PARITY_UNDF as i8;
                    if undefined {
                        num_miss_undf = num_miss_undf.wrapping_add(1);
                    }
                    if add_sb {
                        if picr.num_sb_in2_only < 0 || picr.num_sb_undef_in2_only < 0 {
                            return Err(SourceHeapError::SourceIntegerOverflow);
                        }
                        if picr.num_sb_in2_only < ICR_MAX_SB_IN2_ONLY as i32 {
                            picr.sb_in2_only[picr.num_sb_in2_only as usize] = j2 as AT_NUMB;
                            picr.num_sb_in2_only = picr.num_sb_in2_only.wrapping_add(1);
                        }
                        if undefined && picr.num_sb_undef_in2_only < ICR_MAX_SB_UNDF as i32 {
                            picr.sb_undef_in2_only[picr.num_sb_undef_in2_only as usize] = j1 as AT_NUMB;
                            picr.num_sb_undef_in2_only = picr.num_sb_undef_in2_only.wrapping_add(1);
                        }
                    }
                    j2 += 1;
                }
            }
            while j1 < count1 {
                num_in1_only = num_in1_only.wrapping_add(1);
                let undefined = parity1[j1] == AB_PARITY_UNDF as i8;
                if undefined {
                    num_extra_undf = num_extra_undf.wrapping_add(1);
                }
                if add_sb {
                    if picr.num_sb_in1_only < 0 || picr.num_sb_undef_in1_only < 0 {
                        return Err(SourceHeapError::SourceIntegerOverflow);
                    }
                    if picr.num_sb_in1_only < ICR_MAX_SB_IN1_ONLY as i32 {
                        picr.sb_in1_only[picr.num_sb_in1_only as usize] = j1 as AT_NUMB;
                        picr.num_sb_in1_only = picr.num_sb_in1_only.wrapping_add(1);
                    }
                    if undefined && picr.num_sb_undef_in1_only < ICR_MAX_SB_UNDF as i32 {
                        picr.sb_undef_in1_only[picr.num_sb_undef_in1_only as usize] = j1 as AT_NUMB;
                        picr.num_sb_undef_in1_only = picr.num_sb_undef_in1_only.wrapping_add(1);
                    }
                }
                j1 += 1;
            }
            while j2 < count2 {
                num_in2_only = num_in2_only.wrapping_add(1);
                let undefined = parity2[j2] == AB_PARITY_UNDF as i8;
                if undefined {
                    num_miss_undf = num_miss_undf.wrapping_add(1);
                }
                if add_sb {
                    if picr.num_sb_in2_only < 0 || picr.num_sb_undef_in2_only < 0 {
                        return Err(SourceHeapError::SourceIntegerOverflow);
                    }
                    if picr.num_sb_in2_only < ICR_MAX_SB_IN2_ONLY as i32 {
                        picr.sb_in2_only[picr.num_sb_in2_only as usize] = j2 as AT_NUMB;
                        picr.num_sb_in2_only = picr.num_sb_in2_only.wrapping_add(1);
                    }
                    if undefined && picr.num_sb_undef_in2_only < ICR_MAX_SB_UNDF as i32 {
                        picr.sb_undef_in2_only[picr.num_sb_undef_in2_only as usize] = j1 as AT_NUMB;
                        picr.num_sb_undef_in2_only = picr.num_sb_undef_in2_only.wrapping_add(1);
                    }
                }
                j2 += 1;
            }
            if num_dif != 0 {
                ret |= tagInchiCompareDiffBits_INCHIDIFF_SB_PARITY as INCHI_MODE;
            }
            if num_in1_only != 0 {
                if num_extra_undf != 0 {
                    ret |= tagInchiCompareDiffBits_INCHIDIFF_SB_EXTRA_UNDF as INCHI_MODE;
                }
                if num_in1_only != num_extra_undf {
                    ret |= tagInchiCompareDiffBits_INCHIDIFF_SB_EXTRA as INCHI_MODE;
                }
            }
            if num_in2_only != 0 {
                if num_miss_undf != 0 {
                    ret |= tagInchiCompareDiffBits_INCHIDIFF_SB_MISS_UNDF as INCHI_MODE;
                }
                if num_in2_only != num_miss_undf {
                    ret |= tagInchiCompareDiffBits_INCHIDIFF_SB_MISS as INCHI_MODE;
                }
            }
        }
    }
    Ok(ret)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn CompareReversedINChI3(
    heap: &mut SourceHeap,
    i1: Option<&INChI>,
    i2: Option<&INChI>,
    a1: Option<&INChI_Aux>,
    a2: Option<&INChI_Aux>,
    err: &mut i32,
) -> Result<INCHI_MODE, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:2644 CompareReversedINChI3
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
INCHI_MODE CompareReversedINChI3(INChI* i1 /* InChI from reversed struct */,
    INChI* i2 /* input InChI */,
    INChI_Aux* a1,
    INChI_Aux* a2,
    int* err)
{
    INCHI_MODE ret = 0;
    INChI_Stereo* Stereo1 = NULL, * Stereo2 = NULL;
    int  n1, n2, m, j, j1, j2, ret2, num_H1, num_H2;
    ICR icr;
    ICR* picr = &icr;

    *err = 0;

    memset(picr, 0, sizeof(*picr)); /* djb-rwth: memset_s C11/Annex K variant? */

    if (i1 == NULL && i2 == NULL)
    {
        return 0;
    }
    if ((i1 == NULL) ^ (i2 == NULL))
    {
        ret |= INCHIDIFF_PROBLEM; /* one InChI exists while another doesn't */
        goto exit_function;
    }

    if (i1->nErrorCode == i2->nErrorCode)
    {
        if (i1->nErrorCode)
        {
            ret |= INCHIDIFF_PROBLEM; /* both InChI have same error codes */
            goto exit_function;
        }
    }
    else
    {
        ret |= INCHIDIFF_PROBLEM; /* at least one InChI has an error code */
        goto exit_function;
    }

    if (i1->nNumberOfAtoms != i2->nNumberOfAtoms)
    {
        ret |= INCHIDIFF_NUM_AT;
        goto exit_function;
    }
    if (i1->nNumberOfAtoms > 0)
    {
        if (memcmp(i1->nAtom, i2->nAtom, i1->nNumberOfAtoms * sizeof(i1->nAtom[0])))
        {
            ret |= INCHIDIFF_ATOMS;
            goto exit_function;
        }
        /* INCHIDIFF_NON_TAUT_H,  INCHIDIFF_MORE_FH, INCHIDIFF_LESS_FH */
        if (memcmp(i1->nNum_H, i2->nNum_H, i1->nNumberOfAtoms * sizeof(i1->nNum_H[0])))
        {
            ret |= INCHIDIFF_POSITION_H;
            for (j1 = 0; j1 < i1->nNumberOfAtoms; j1++)
            {
                if (i1->nNum_H[j1] != i2->nNum_H[j1] && picr->num_diff_pos_H < ICR_MAX_DIFF_FIXED_H)
                {
                    picr->diff_pos_H_at[picr->num_diff_pos_H] = j1;
                    picr->diff_pos_H_nH[picr->num_diff_pos_H] = i1->nNum_H[j1] - i2->nNum_H[j1];
                    picr->num_diff_pos_H++;
                }
            }
        }
        /* fixed H */
        if (i1->nNum_H_fixed || i2->nNum_H_fixed)
        {
            int bHasFixedH1 = 0, bHasFixedH2 = 0, i;
            if (i1->nNum_H_fixed)
            {
                for (i = 0; i < i1->nNumberOfAtoms; i++)
                {
                    if (i1->nNum_H_fixed[i])
                    {
                        bHasFixedH1++;
                    }
                }
            }
            if (i2->nNum_H_fixed)
            {
                for (i = 0; i < i2->nNumberOfAtoms; i++)
                {
                    if (i2->nNum_H_fixed[i])
                    {
                        bHasFixedH2++;
                    }
                }
            }
            if (bHasFixedH1 && !bHasFixedH2)
            {
                for (i = j = 0; i < i1->nNumberOfAtoms; i++)
                {
                    if (i1->nNum_H_fixed[i])
                    {
                        if (j < ICR_MAX_DIFF_FIXED_H)
                        {
                            picr->fixed_H_at1_more[j] = i;
                            picr->fixed_H_nH1_more[j] = i1->nNum_H_fixed[i];
                            j++;
                        }
                    }
                }
                picr->num_fixed_H1_more = j;
                ret |= INCHIDIFF_MORE_FH; /* Extra Fixed-H */
            }
            else
            {
                if (!bHasFixedH1 && bHasFixedH2)
                {
                    for (i = j = 0; i < i2->nNumberOfAtoms; i++)
                    {
                        if (i2->nNum_H_fixed[i])
                        {
                            if (j < ICR_MAX_DIFF_FIXED_H)
                            {
                                picr->fixed_H_at2_more[j] = i;
                                picr->fixed_H_nH2_more[j] = i2->nNum_H_fixed[i];
                                j++;
                            }
                        }
                    }
                    picr->num_fixed_H2_more = j;
                    ret |= INCHIDIFF_LESS_FH; /* Missed Fixed-H */
                }
                else
                {
                    if (bHasFixedH1 && bHasFixedH2 &&
                        memcmp(i1->nNum_H_fixed, i2->nNum_H_fixed, i1->nNumberOfAtoms * sizeof(i1->nNum_H_fixed[0])))
                    {
                        for (i = j1 = j2 = 0; i < i1->nNumberOfAtoms; i++)
                        {
                            if (i1->nNum_H_fixed[i] > i2->nNum_H_fixed[i])
                            {
                                if (j1 < ICR_MAX_DIFF_FIXED_H)
                                {
                                    picr->fixed_H_at1_more[j1] = i;
                                    picr->fixed_H_nH1_more[j1] = i1->nNum_H_fixed[i] - i2->nNum_H_fixed[i];
                                    j1++;
                                }
                            }
                            else
                                if (i1->nNum_H_fixed[i] < i2->nNum_H_fixed[i])
                                {
                                    if (j2 < ICR_MAX_DIFF_FIXED_H)
                                    {
                                        picr->fixed_H_at2_more[j2] = i;
                                        picr->fixed_H_nH2_more[j2] = i2->nNum_H_fixed[i] - i1->nNum_H_fixed[i];
                                        j2++;
                                    }
                                }
                        }
                        ret |= (j1 ? INCHIDIFF_MORE_FH : 0) | (j2 ? INCHIDIFF_LESS_FH : 0);
                        picr->num_fixed_H1_more = j1;
                        picr->num_fixed_H2_more = j2;
                    }
                }
            }
        }
    }

    /* compare formulas and H */
    num_H1 = 0;
    num_H2 = 0;
    ret2 = CompareHillFormulasNoH(i1->szHillFormula, i2->szHillFormula, &num_H1, &num_H2);
    picr->tot_num_H1 = num_H1;
    picr->tot_num_H2 = num_H2;
    if (ret2)
    {
        ret |= INCHIDIFF_NUM_EL;
        goto exit_function;
    }
    if (num_H1 > num_H2)
    {
        ret |= INCHIDIFF_MORE_H;
    }
    if (num_H1 < num_H2)
    {
        ret |= INCHIDIFF_LESS_H;
    }

    if (i1->lenConnTable != i2->lenConnTable)
    {
        ret |= INCHIDIFF_CON_LEN;
        goto exit_function;
    }
    else
    {
        if (i1->lenConnTable > 0 && memcmp(i1->nConnTable, i2->nConnTable, i1->lenConnTable * sizeof(i1->nConnTable[0])))
        {
            ret |= INCHIDIFF_CON_TBL;
            goto exit_function;
        }
    }
    /* output special cases: different number of t-groups, different sizes of t-groups, different endpoints */
    /* in isotopic or deprotonated cases i1->lenTautomer == 1 && i1->nTautomer[0] = 0 */
        /*
    if ( i1->lenTautomer != i2->lenTautomer && (i1->lenTautomer > 1 || i2->lenTautomer > 1) )
    {
        ret |=  INCHIDIFF_TAUT_LEN;
    }
    */

    /* compare number of t-groups */
    n1 = i1->lenTautomer ? i1->nTautomer[0] : 0;
    n2 = i2->lenTautomer ? i2->nTautomer[0] : 0;
    if (!n1 && n2)
    {
        ret |= INCHIDIFF_NO_TAUT;
    }
    else
    {
        if (n1 && !n2)
        {
            ret |= INCHIDIFF_WRONG_TAUT;
        }
        else
        {
            if (n1 == 1 && n2 > 1)
            {
                ret |= INCHIDIFF_SINGLE_TG;
            }
            else
            {
                if (n1 > 1 && n2 == 1)
                {
                    ret |= INCHIDIFF_MULTIPLE_TG;
                }
                else
                {
                    if (n1 != n2)
                    {
                        ret |= INCHIDIFF_NUM_TG;
                    }
                }
            }
        }
    }

    if (n1 || n2)
    {
        /* number of endpoints */
        int num1 = 0, num2 = 0, num_M1 = 0, num_M2 = 0;
        int len, num_eq, num_in1_only, num_in2_only;
        AT_NUMB* pe1 = (AT_NUMB*)inchi_malloc(((long long)i1->lenTautomer + 1) * sizeof(pe1[0])); /* djb-rwth: cast operator added */
        AT_NUMB* pe2 = (AT_NUMB*)inchi_malloc(((long long)i2->lenTautomer + 1) * sizeof(pe2[0])); /* djb-rwth: cast operator added */
        num_H1 = num_H2 = 0;
        /* collect endpoints, H, (-) */
        if (!pe1 || !pe2)
        {
            if (pe1) inchi_free(pe1);
            if (pe2) inchi_free(pe2);
            *err = RI_ERR_ALLOC; /* allocation error */
            goto exit_function;
        }
        for (m = 1; m < i1->lenTautomer; m += len)
        {
            len = i1->nTautomer[m++];
            num_H1 += i1->nTautomer[m];
            num_M1 += i1->nTautomer[m + 1];
            for (j = 2; j < len; j++)
            {
                pe1[num1++] = i1->nTautomer[m + j];
            }
        }
        for (m = 1; m < i2->lenTautomer; m += len)
        {
            len = i2->nTautomer[m++];
            num_H2 += i2->nTautomer[m];
            num_M2 += i2->nTautomer[m + 1];
            for (j = 2; j < len; j++)
            {
                pe2[num2++] = i2->nTautomer[m + j];
            }
        }
        picr->num_taut_H1 = num_H1;
        picr->num_taut_H2 = num_H2;
        picr->num_taut_M1 = num_M1;
        picr->num_taut_M2 = num_M2;
        /* sort endpoints */
        insertions_sort_AT_NUMB(pe1, num1);
        insertions_sort_AT_NUMB(pe2, num2);
        /* compare */
        /*
        if ( num1 < num2 ) {
            ret |= INCHIDIFF_LESS_TG_ENDP;
        } else
        if ( num1 > num2 ) {
            ret |= INCHIDIFF_MORE_TG_ENDP;
        }
        */
        /* compare all */
        num_eq = num_in1_only = num_in2_only = 0;
        for (j1 = j2 = 0; j1 < num1 && j2 < num2; )
        {
            if (pe1[j1] == pe2[j2])
            {
                j1++;
                j2++;
                num_eq++;
            }
            else
            {
                if (pe1[j1] < pe2[j1])
                {
                    if (picr->num_endp_in1_only < ICR_MAX_ENDP_IN1_ONLY)
                    {
                        picr->endp_in1_only[picr->num_endp_in1_only++] = pe1[j1];
                    }
                    j1++;
                    num_in1_only++;
                }
                else
                {
                    if (picr->num_endp_in2_only < ICR_MAX_ENDP_IN2_ONLY)
                    {
                        picr->endp_in2_only[picr->num_endp_in2_only++] = pe2[j2];
                    }
                    j2++;
                    num_in2_only++;
                }
            }
        }
        while (j1 < num1)
        {
            if (picr->num_endp_in1_only < ICR_MAX_ENDP_IN1_ONLY)
            {
                picr->endp_in1_only[picr->num_endp_in1_only++] = pe1[j1];
            }
            j1++;
            num_in1_only++;
        }
        while (j2 < num2)
        {
            if (picr->num_endp_in2_only < ICR_MAX_ENDP_IN2_ONLY)
            {
                picr->endp_in2_only[picr->num_endp_in2_only++] = pe2[j2];
            }
            j2++;
            num_in2_only++;
        }
        if (num_in1_only)
        {
            ret |= INCHIDIFF_EXTRA_TG_ENDP;
        }
        if (num_in2_only)
        {
            ret |= INCHIDIFF_MISS_TG_ENDP;
        }
        if (!num_in1_only && !num_in2_only && num_eq)
        {
            ; /* same t-groups endpoints */
        }
        else
        {
            ret |= INCHIDIFF_DIFF_TG_ENDP;
        }
        inchi_free(pe1);
        inchi_free(pe2);
    }

    if ((i1->lenTautomer > 1 && i2->lenTautomer > 1) &&
        (i1->lenTautomer != i2->lenTautomer ||
            memcmp(i1->nTautomer, i2->nTautomer, i1->lenTautomer * sizeof(i1->nTautomer[0]))))
        ret |= INCHIDIFF_TG;

    if (i1->nNumberOfIsotopicAtoms != i2->nNumberOfIsotopicAtoms)
    {
        ret |= INCHIDIFF_NUM_ISO_AT;
    }
    else
    {
        if (i1->nNumberOfIsotopicAtoms > 0 && memcmp(i1->IsotopicAtom, i2->IsotopicAtom, i1->nNumberOfIsotopicAtoms * sizeof(i1->IsotopicAtom[0])))
        {
            ret |= INCHIDIFF_ISO_AT;
        }
    }
    if (i1->nTotalCharge != i2->nTotalCharge)
        ret |= INCHIDIFF_CHARGE;
    if (a1 && a1->nNumRemovedProtons && (!a2 || a2->nNumRemovedProtons != a1->nNumRemovedProtons))
    {
        ret |= INCHIDIFF_REM_PROT;
    }
    if (a1 && (!a2 ||
        a2->nNumRemovedIsotopicH[0] != a1->nNumRemovedIsotopicH[0] ||
        a2->nNumRemovedIsotopicH[1] != a1->nNumRemovedIsotopicH[1] ||
        a2->nNumRemovedIsotopicH[2] != a1->nNumRemovedIsotopicH[2]))
    {
        ret |= INCHIDIFF_REM_ISO_H;
    }

    /*
    if ( i1->nPossibleLocationsOfIsotopicH && i2->nPossibleLocationsOfIsotopicH ) {
        if ( i1->nPossibleLocationsOfIsotopicH[0] != i2->nPossibleLocationsOfIsotopicH[0] ||
             memcmp(i1->nPossibleLocationsOfIsotopicH, i2->nPossibleLocationsOfIsotopicH,
                    sizeof(i1->nPossibleLocationsOfIsotopicH[0])*i1->nPossibleLocationsOfIsotopicH[0]) )
            return 18;
    } else
    if ( !i1->nPossibleLocationsOfIsotopicH != !i2->nPossibleLocationsOfIsotopicH ) {
        return 19;
    }
    */
    if (i1->StereoIsotopic &&
        i1->StereoIsotopic->nNumberOfStereoBonds + i1->StereoIsotopic->nNumberOfStereoCenters)
    {
        Stereo1 = i1->StereoIsotopic;
    }
    else
    {
        Stereo1 = i1->Stereo;
    }
    if (i2->StereoIsotopic &&
        i2->StereoIsotopic->nNumberOfStereoBonds + i2->StereoIsotopic->nNumberOfStereoCenters)
    {
        Stereo2 = i2->StereoIsotopic;
    }
    else
    {
        Stereo2 = i2->Stereo;
    }
    ret |= CompareReversedStereoINChI3(Stereo1, Stereo2, picr);

exit_function:
    picr->flags = ret;

    return ret;
}
    */
    // END INCHI C FUNCTION: CompareReversedINChI3
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CompareReversedINChI3
    // INCHI✔️❌: READ_INCHI_STRING=1 includes this production helper.
    // INCHI✔️❌: COMPILE_ANSI_ONLY and TARGET_API_LIB do not alter this function body.
    // INCHI✔️❌: inchi_malloc/free resolve to the active GCC/Linux mode.h macros.
    // INCHI✔️❌: Rust-owned endpoint staging and checked SourceHeap access add allocation and bounds-check overhead.
    // END INCHI ACTIVE MACRO CONFIGURATION: CompareReversedINChI3

    *err = 0;
    let mut picr = ICR::default();
    let mut ret = 0_u32;

    macro_rules! finish {
        () => {{
            picr.flags = INCHI_MODE::from(ret);
            return Ok(INCHI_MODE::from(ret));
        }};
    }

    let (i1, i2) = match (i1, i2) {
        (None, None) => return Ok(0),
        (Some(i1), Some(i2)) => (i1, i2),
        _ => {
            ret |= tagInchiCompareDiffBits_INCHIDIFF_PROBLEM;
            finish!();
        }
    };

    if i1.nErrorCode == i2.nErrorCode {
        if i1.nErrorCode != 0 {
            ret |= tagInchiCompareDiffBits_INCHIDIFF_PROBLEM;
            finish!();
        }
    } else {
        ret |= tagInchiCompareDiffBits_INCHIDIFF_PROBLEM;
        finish!();
    }

    if i1.nNumberOfAtoms != i2.nNumberOfAtoms {
        ret |= tagInchiCompareDiffBits_INCHIDIFF_NUM_AT;
        finish!();
    }
    if i1.nNumberOfAtoms > 0 {
        let atom_count = usize::try_from(i1.nNumberOfAtoms)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let atoms1 = heap
            .slice(i1.nAtom.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let atoms2 = heap
            .slice(i2.nAtom.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if atoms1 != atoms2 {
            ret |= tagInchiCompareDiffBits_INCHIDIFF_ATOMS;
            finish!();
        }

        let hydrogens1 = heap
            .slice(i1.nNum_H.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let hydrogens2 = heap
            .slice(i2.nNum_H.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if hydrogens1 != hydrogens2 {
            ret |= tagInchiCompareDiffBits_INCHIDIFF_POSITION_H;
            for index in 0..atom_count {
                if hydrogens1[index] != hydrogens2[index]
                    && picr.num_diff_pos_H < ICR_MAX_DIFF_FIXED_H as i32
                {
                    let output = picr.num_diff_pos_H as usize;
                    picr.diff_pos_H_at[output] = index as AT_NUMB;
                    picr.diff_pos_H_nH[output] = hydrogens1[index].wrapping_sub(hydrogens2[index]);
                    picr.num_diff_pos_H = picr.num_diff_pos_H.wrapping_add(1);
                }
            }
        }

        if !i1.nNum_H_fixed.is_null() || !i2.nNum_H_fixed.is_null() {
            let fixed1 = if i1.nNum_H_fixed.is_null() {
                None
            } else {
                Some(
                    heap.slice(i1.nNum_H_fixed.as_const())?
                        .get(..atom_count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                )
            };
            let fixed2 = if i2.nNum_H_fixed.is_null() {
                None
            } else {
                Some(
                    heap.slice(i2.nNum_H_fixed.as_const())?
                        .get(..atom_count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                )
            };
            let has_fixed1 = fixed1.map_or(0_i32, |values| {
                values.iter().fold(0_i32, |count, value| {
                    count.wrapping_add(i32::from(*value != 0))
                })
            });
            let has_fixed2 = fixed2.map_or(0_i32, |values| {
                values.iter().fold(0_i32, |count, value| {
                    count.wrapping_add(i32::from(*value != 0))
                })
            });

            if has_fixed1 != 0 && has_fixed2 == 0 {
                let values = fixed1.expect("nonzero count requires a source array");
                let mut count = 0_i32;
                for (index, &value) in values.iter().enumerate() {
                    if value != 0 && count < ICR_MAX_DIFF_FIXED_H as i32 {
                        picr.fixed_H_at1_more[count as usize] = index as AT_NUMB;
                        picr.fixed_H_nH1_more[count as usize] = value;
                        count = count.wrapping_add(1);
                    }
                }
                picr.num_fixed_H1_more = count;
                ret |= tagInchiCompareDiffBits_INCHIDIFF_MORE_FH;
            } else if has_fixed1 == 0 && has_fixed2 != 0 {
                let values = fixed2.expect("nonzero count requires a source array");
                let mut count = 0_i32;
                for (index, &value) in values.iter().enumerate() {
                    if value != 0 && count < ICR_MAX_DIFF_FIXED_H as i32 {
                        picr.fixed_H_at2_more[count as usize] = index as AT_NUMB;
                        picr.fixed_H_nH2_more[count as usize] = value;
                        count = count.wrapping_add(1);
                    }
                }
                picr.num_fixed_H2_more = count;
                ret |= tagInchiCompareDiffBits_INCHIDIFF_LESS_FH;
            } else if has_fixed1 != 0 && has_fixed2 != 0 {
                let values1 = fixed1.expect("nonzero count requires a source array");
                let values2 = fixed2.expect("nonzero count requires a source array");
                if values1 != values2 {
                    let mut count1 = 0_i32;
                    let mut count2 = 0_i32;
                    for index in 0..atom_count {
                        if values1[index] > values2[index] {
                            if count1 < ICR_MAX_DIFF_FIXED_H as i32 {
                                picr.fixed_H_at1_more[count1 as usize] = index as AT_NUMB;
                                picr.fixed_H_nH1_more[count1 as usize] =
                                    values1[index].wrapping_sub(values2[index]);
                                count1 = count1.wrapping_add(1);
                            }
                        } else if values1[index] < values2[index]
                            && count2 < ICR_MAX_DIFF_FIXED_H as i32
                        {
                            picr.fixed_H_at2_more[count2 as usize] = index as AT_NUMB;
                            picr.fixed_H_nH2_more[count2 as usize] =
                                values2[index].wrapping_sub(values1[index]);
                            count2 = count2.wrapping_add(1);
                        }
                    }
                    if count1 != 0 {
                        ret |= tagInchiCompareDiffBits_INCHIDIFF_MORE_FH;
                    }
                    if count2 != 0 {
                        ret |= tagInchiCompareDiffBits_INCHIDIFF_LESS_FH;
                    }
                    picr.num_fixed_H1_more = count1;
                    picr.num_fixed_H2_more = count2;
                }
            }
        }
    }

    let mut total_h1 = 0_i32;
    let mut total_h2 = 0_i32;
    let formula_difference = CompareHillFormulasNoH(
        heap,
        i1.szHillFormula.as_const(),
        i2.szHillFormula.as_const(),
        &mut total_h1,
        &mut total_h2,
    )?;
    picr.tot_num_H1 = total_h1;
    picr.tot_num_H2 = total_h2;
    if formula_difference != 0 {
        ret |= tagInchiCompareDiffBits_INCHIDIFF_NUM_EL;
        finish!();
    }
    if total_h1 > total_h2 {
        ret |= tagInchiCompareDiffBits_INCHIDIFF_MORE_H;
    }
    if total_h1 < total_h2 {
        ret |= tagInchiCompareDiffBits_INCHIDIFF_LESS_H;
    }

    if i1.lenConnTable != i2.lenConnTable {
        ret |= tagInchiCompareDiffBits_INCHIDIFF_CON_LEN;
        finish!();
    }
    if i1.lenConnTable > 0 {
        let connection_count =
            usize::try_from(i1.lenConnTable).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let connections1 = heap
            .slice(i1.nConnTable.as_const())?
            .get(..connection_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let connections2 = heap
            .slice(i2.nConnTable.as_const())?
            .get(..connection_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if connections1 != connections2 {
            ret |= tagInchiCompareDiffBits_INCHIDIFF_CON_TBL;
            finish!();
        }
    }

    let num_groups1 = if i1.lenTautomer != 0 {
        *heap
            .slice(i1.nTautomer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)? as i32
    } else {
        0
    };
    let num_groups2 = if i2.lenTautomer != 0 {
        *heap
            .slice(i2.nTautomer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)? as i32
    } else {
        0
    };
    if num_groups1 == 0 && num_groups2 != 0 {
        ret |= tagInchiCompareDiffBits_INCHIDIFF_NO_TAUT;
    } else if num_groups1 != 0 && num_groups2 == 0 {
        ret |= tagInchiCompareDiffBits_INCHIDIFF_WRONG_TAUT;
    } else if num_groups1 == 1 && num_groups2 > 1 {
        ret |= tagInchiCompareDiffBits_INCHIDIFF_SINGLE_TG;
    } else if num_groups1 > 1 && num_groups2 == 1 {
        ret |= tagInchiCompareDiffBits_INCHIDIFF_MULTIPLE_TG;
    } else if num_groups1 != num_groups2 {
        ret |= tagInchiCompareDiffBits_INCHIDIFF_NUM_TG;
    }

    if num_groups1 != 0 || num_groups2 != 0 {
        let byte_count1 = u64::try_from(
            (i64::from(i1.lenTautomer) + 1)
                .checked_mul(std::mem::size_of::<AT_NUMB>() as i64)
                .ok_or(SourceHeapError::AllocationSizeOverflow)?,
        )
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        let byte_count2 = u64::try_from(
            (i64::from(i2.lenTautomer) + 1)
                .checked_mul(std::mem::size_of::<AT_NUMB>() as i64)
                .ok_or(SourceHeapError::AllocationSizeOverflow)?,
        )
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;

        let allocation1 = inchi_malloc(heap, byte_count1);
        let allocation2 = inchi_malloc(heap, byte_count2);
        let (endpoint_allocation1, endpoint_allocation2) = match (allocation1, allocation2) {
            (Ok(pointer1), Ok(pointer2)) => (pointer1, pointer2),
            (result1, result2) => {
                if let Ok(pointer) = result1.as_ref() {
                    inchi_free(heap, *pointer)?;
                }
                if let Ok(pointer) = result2.as_ref() {
                    inchi_free(heap, *pointer)?;
                }
                if matches!(&result1, Err(SourceHeapError::AllocationFailed))
                    || matches!(&result2, Err(SourceHeapError::AllocationFailed))
                {
                    *err = RI_ERR_ALLOC;
                    finish!();
                }
                return Err(result1.err().or_else(|| result2.err()).ok_or(
                    SourceHeapError::AllocationFailed,
                )?);
            }
        };

        let endpoint_result = (|| -> Result<
            (Vec<AT_NUMB>, usize, Vec<AT_NUMB>, usize, i32, i32, i32, i32),
            SourceHeapError,
        > {
                let tautomer_count1 = usize::try_from(i1.lenTautomer)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let tautomer_count2 = usize::try_from(i2.lenTautomer)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let tautomer1 = if tautomer_count1 == 0 {
                    Vec::new()
                } else {
                    heap.slice(i1.nTautomer.as_const())?
                        .get(..tautomer_count1)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .to_vec()
                };
                let tautomer2 = if tautomer_count2 == 0 {
                    Vec::new()
                } else {
                    heap.slice(i2.nTautomer.as_const())?
                        .get(..tautomer_count2)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .to_vec()
                };

                // INCHI✔️❌: AT_NUMB* pe1 = (AT_NUMB*)inchi_malloc(((long long)i1->lenTautomer + 1) * sizeof(pe1[0]));
                // INCHI✔️❌: AT_NUMB* pe2 = (AT_NUMB*)inchi_malloc(((long long)i2->lenTautomer + 1) * sizeof(pe2[0]));
                // Keep the complete source allocation extents, not only the initialized
                // endpoint prefixes. The later official `pe2[j1]` access is measured
                // against this allocation. `inchi_malloc`'s SourceHeap model supplies
                // zero-filled bytes for the otherwise uninitialized tail.
                let endpoint_extent1 = tautomer_count1
                    .checked_add(1)
                    .ok_or(SourceHeapError::AllocationSizeOverflow)?;
                let endpoint_extent2 = tautomer_count2
                    .checked_add(1)
                    .ok_or(SourceHeapError::AllocationSizeOverflow)?;
                let mut endpoints1 = vec![0; endpoint_extent1];
                let mut endpoints2 = vec![0; endpoint_extent2];
                let mut endpoint_count1 = 0_usize;
                let mut endpoint_count2 = 0_usize;
                let mut taut_h1 = 0_i32;
                let mut taut_h2 = 0_i32;
                let mut taut_m1 = 0_i32;
                let mut taut_m2 = 0_i32;

                let mut m = 1_usize;
                while m < tautomer_count1 {
                    let length = usize::from(tautomer1[m]);
                    m += 1;
                    let record = tautomer1
                        .get(m..m.saturating_add(length))
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    taut_h1 = taut_h1.wrapping_add(i32::from(
                        *record.first().ok_or(SourceHeapError::PointerOutOfBounds)?,
                    ));
                    taut_m1 = taut_m1.wrapping_add(i32::from(
                        *record.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?,
                    ));
                    let source_endpoints = record
                        .get(2..length)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let next_count = endpoint_count1
                        .checked_add(source_endpoints.len())
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                    endpoints1
                        .get_mut(endpoint_count1..next_count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .copy_from_slice(source_endpoints);
                    endpoint_count1 = next_count;
                    m = m
                        .checked_add(length)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                }

                m = 1;
                while m < tautomer_count2 {
                    let length = usize::from(tautomer2[m]);
                    m += 1;
                    let record = tautomer2
                        .get(m..m.saturating_add(length))
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    taut_h2 = taut_h2.wrapping_add(i32::from(
                        *record.first().ok_or(SourceHeapError::PointerOutOfBounds)?,
                    ));
                    taut_m2 = taut_m2.wrapping_add(i32::from(
                        *record.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?,
                    ));
                    let source_endpoints = record
                        .get(2..length)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let next_count = endpoint_count2
                        .checked_add(source_endpoints.len())
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                    endpoints2
                        .get_mut(endpoint_count2..next_count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .copy_from_slice(source_endpoints);
                    endpoint_count2 = next_count;
                    m = m
                        .checked_add(length)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                }

                let count1 = i32::try_from(endpoint_count1)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let count2 = i32::try_from(endpoint_count2)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                insertions_sort_AT_NUMB(&mut endpoints1, count1)?;
                insertions_sort_AT_NUMB(&mut endpoints2, count2)?;
                Ok((
                    endpoints1,
                    endpoint_count1,
                    endpoints2,
                    endpoint_count2,
                    taut_h1,
                    taut_h2,
                    taut_m1,
                    taut_m2,
                ))
            })();

        let free1 = inchi_free(heap, endpoint_allocation1);
        let free2 = inchi_free(heap, endpoint_allocation2);
        free1?;
        free2?;
        let (
            endpoints1,
            endpoint_count1,
            endpoints2,
            endpoint_count2,
            taut_h1,
            taut_h2,
            taut_m1,
            taut_m2,
        ) = endpoint_result?;

        picr.num_taut_H1 = taut_h1;
        picr.num_taut_H2 = taut_h2;
        picr.num_taut_M1 = taut_m1;
        picr.num_taut_M2 = taut_m2;

        let mut index1 = 0_usize;
        let mut index2 = 0_usize;
        let mut num_equal = 0_i32;
        let mut num_in1_only = 0_i32;
        let mut num_in2_only = 0_i32;
        while index1 < endpoint_count1 && index2 < endpoint_count2 {
            if endpoints1[index1] == endpoints2[index2] {
                index1 += 1;
                index2 += 1;
                num_equal = num_equal.wrapping_add(1);
            // INCHI✔️❌: if (pe1[j1] < pe2[j1])
            // This intentionally uses `index1`, exactly as the official source
            // does, while retaining checked access to the complete allocation.
            } else if endpoints1[index1]
                < *endpoints2
                    .get(index1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
            {
                if picr.num_endp_in1_only < ICR_MAX_ENDP_IN1_ONLY as i32 {
                    picr.endp_in1_only[picr.num_endp_in1_only as usize] = endpoints1[index1];
                    picr.num_endp_in1_only = picr.num_endp_in1_only.wrapping_add(1);
                }
                index1 += 1;
                num_in1_only = num_in1_only.wrapping_add(1);
            } else {
                if picr.num_endp_in2_only < ICR_MAX_ENDP_IN2_ONLY as i32 {
                    picr.endp_in2_only[picr.num_endp_in2_only as usize] = endpoints2[index2];
                    picr.num_endp_in2_only = picr.num_endp_in2_only.wrapping_add(1);
                }
                index2 += 1;
                num_in2_only = num_in2_only.wrapping_add(1);
            }
        }
        while index1 < endpoint_count1 {
            if picr.num_endp_in1_only < ICR_MAX_ENDP_IN1_ONLY as i32 {
                picr.endp_in1_only[picr.num_endp_in1_only as usize] = endpoints1[index1];
                picr.num_endp_in1_only = picr.num_endp_in1_only.wrapping_add(1);
            }
            index1 += 1;
            num_in1_only = num_in1_only.wrapping_add(1);
        }
        while index2 < endpoint_count2 {
            if picr.num_endp_in2_only < ICR_MAX_ENDP_IN2_ONLY as i32 {
                picr.endp_in2_only[picr.num_endp_in2_only as usize] = endpoints2[index2];
                picr.num_endp_in2_only = picr.num_endp_in2_only.wrapping_add(1);
            }
            index2 += 1;
            num_in2_only = num_in2_only.wrapping_add(1);
        }
        if num_in1_only != 0 {
            ret |= tagInchiCompareDiffBits_INCHIDIFF_EXTRA_TG_ENDP;
        }
        if num_in2_only != 0 {
            ret |= tagInchiCompareDiffBits_INCHIDIFF_MISS_TG_ENDP;
        }
        if num_in1_only != 0 || num_in2_only != 0 || num_equal == 0 {
            ret |= tagInchiCompareDiffBits_INCHIDIFF_DIFF_TG_ENDP;
        }
    }

    if i1.lenTautomer > 1 && i2.lenTautomer > 1 {
        let raw_tautomer_difference = if i1.lenTautomer != i2.lenTautomer {
            true
        } else {
            let count = usize::try_from(i1.lenTautomer)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            heap.slice(i1.nTautomer.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                != heap
                    .slice(i2.nTautomer.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
        };
        if raw_tautomer_difference {
            ret |= tagInchiCompareDiffBits_INCHIDIFF_TG;
        }
    }

    if i1.nNumberOfIsotopicAtoms != i2.nNumberOfIsotopicAtoms {
        ret |= tagInchiCompareDiffBits_INCHIDIFF_NUM_ISO_AT;
    } else if i1.nNumberOfIsotopicAtoms > 0 {
        let count = usize::try_from(i1.nNumberOfIsotopicAtoms)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if heap
            .slice(i1.IsotopicAtom.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            != heap
                .slice(i2.IsotopicAtom.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        {
            ret |= tagInchiCompareDiffBits_INCHIDIFF_ISO_AT;
        }
    }
    if i1.nTotalCharge != i2.nTotalCharge {
        ret |= tagInchiCompareDiffBits_INCHIDIFF_CHARGE;
    }
    if let Some(aux1) = a1 {
        if aux1.nNumRemovedProtons != 0
            && a2.is_none_or(|aux2| aux2.nNumRemovedProtons != aux1.nNumRemovedProtons)
        {
            ret |= tagInchiCompareDiffBits_INCHIDIFF_REM_PROT;
        }
        if a2.is_none_or(|aux2| {
            aux2.nNumRemovedIsotopicH[0] != aux1.nNumRemovedIsotopicH[0]
                || aux2.nNumRemovedIsotopicH[1] != aux1.nNumRemovedIsotopicH[1]
                || aux2.nNumRemovedIsotopicH[2] != aux1.nNumRemovedIsotopicH[2]
        }) {
            ret |= tagInchiCompareDiffBits_INCHIDIFF_REM_ISO_H;
        }
    }

    let isotopic_stereo1 = if i1.StereoIsotopic.is_null() {
        None
    } else {
        heap.slice(i1.StereoIsotopic.as_const())?.first()
    };
    let regular_stereo1 = if i1.Stereo.is_null() {
        None
    } else {
        heap.slice(i1.Stereo.as_const())?.first()
    };
    let stereo1 = if isotopic_stereo1.is_some_and(|stereo| {
        stereo
            .nNumberOfStereoBonds
            .wrapping_add(stereo.nNumberOfStereoCenters)
            != 0
    }) {
        isotopic_stereo1
    } else {
        regular_stereo1
    };

    let isotopic_stereo2 = if i2.StereoIsotopic.is_null() {
        None
    } else {
        heap.slice(i2.StereoIsotopic.as_const())?.first()
    };
    let regular_stereo2 = if i2.Stereo.is_null() {
        None
    } else {
        heap.slice(i2.Stereo.as_const())?.first()
    };
    let stereo2 = if isotopic_stereo2.is_some_and(|stereo| {
        stereo
            .nNumberOfStereoBonds
            .wrapping_add(stereo.nNumberOfStereoCenters)
            != 0
    }) {
        isotopic_stereo2
    } else {
        regular_stereo2
    };

    ret |= CompareReversedStereoINChI3(heap, stereo1, stereo2, &mut picr)? as u32;
    picr.flags = INCHI_MODE::from(ret);
    Ok(INCHI_MODE::from(ret))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source::base::ichirvr1::MakeOneInChIOutOfStrFromINChI;
    use crate::source_types::{
        ALL_TC_GROUPS, COMPONENT_REM_PROTONS, INChI, INChI_Aux, INChI_IsotopicAtom,
        INChI_IsotopicTGroup, INChI_Stereo, OAD_Polymer, OAD_PolymerUnit, OAD_V3000,
        REQ_MODE_NON_ISO, REQ_MODE_TAUT, RI_ERR_ALLOC, SRM, T_GROUP, TG_FLAG_FIX_ISO_FIXEDH_BUG,
        TG_FLAG_FIX_TERM_H_CHRG_BUG, VAL_AT, XYZ_COORD, inp_ATOM, inp_ATOM_STEREO,
    };

    #[test]
    fn source_port__ichirvr7__comparetwopairsofinchi__line_2203() {
        let mut heap = SourceHeap::default();
        let nulls = [SourceMutPointer::null(); TAUT_NUM as usize];
        let mut flags = [7 as INCHI_MODE, 11 as INCHI_MODE];
        assert_eq!(
            CompareTwoPairsOfInChI(&mut heap, nulls, nulls, i32::MIN, &mut flags),
            Ok(0)
        );
        assert_eq!(flags, [7, 11]);

        let equal1_value = reversed_inchi3_fixture(&mut heap);
        let equal1 = heap.allocate_model_storage(vec![equal1_value]).unwrap();
        let equal2_value = reversed_inchi3_fixture(&mut heap);
        let equal2 = heap.allocate_model_storage(vec![equal2_value]).unwrap();
        flags = [0; TAUT_NUM as usize];
        assert_eq!(
            CompareTwoPairsOfInChI(
                &mut heap,
                [equal1, equal1],
                [equal2, equal2],
                i32::MAX,
                &mut flags,
            ),
            Ok(0)
        );
        assert_eq!(flags, [0, 0]);

        flags = [0; TAUT_NUM as usize];
        assert_eq!(
            CompareTwoPairsOfInChI(
                &mut heap,
                [SourceMutPointer::null(), equal1],
                [equal2, equal2],
                0,
                &mut flags,
            ),
            Ok(0)
        );
        assert_eq!(
            flags,
            [
                tagInchiCompareDiffBits_INCHIDIFF_COMP_HLAYER as INCHI_MODE,
                0,
            ]
        );

        flags = [13 as INCHI_MODE, 0];
        assert_eq!(
            CompareTwoPairsOfInChI(
                &mut heap,
                [SourceMutPointer::null(), SourceMutPointer::null()],
                [equal2, equal2],
                -99,
                &mut flags,
            ),
            Ok(0)
        );
        assert_eq!(
            flags,
            [
                13 | tagInchiCompareDiffBits_INCHIDIFF_COMP_NUMBER as INCHI_MODE,
                tagInchiCompareDiffBits_INCHIDIFF_COMP_NUMBER as INCHI_MODE,
            ]
        );

        let mut changed = reversed_inchi3_fixture(&mut heap);
        changed.nNumberOfAtoms = changed.nNumberOfAtoms.wrapping_add(1);
        let changed = heap.allocate_model_storage(vec![changed]).unwrap();
        flags = [0; TAUT_NUM as usize];
        assert_eq!(
            CompareTwoPairsOfInChI(
                &mut heap,
                [equal1, SourceMutPointer::null()],
                [changed, SourceMutPointer::null()],
                1,
                &mut flags,
            ),
            Ok(0)
        );
        assert_eq!(
            flags,
            [tagInchiCompareDiffBits_INCHIDIFF_NUM_AT as INCHI_MODE, 0]
        );

        let mut tautomer1 = reversed_inchi3_fixture(&mut heap);
        tautomer1.nTautomer = heap
            .allocate_model_storage(vec![1_u16, 4, 2, 1, 5, 1])
            .unwrap();
        tautomer1.lenTautomer = 6;
        let tautomer1 = heap.allocate_model_storage(vec![tautomer1]).unwrap();
        let mut tautomer2 = reversed_inchi3_fixture(&mut heap);
        tautomer2.nTautomer = heap
            .allocate_model_storage(vec![1_u16, 4, 1, 0, 2, 5])
            .unwrap();
        tautomer2.lenTautomer = 6;
        let tautomer2 = heap.allocate_model_storage(vec![tautomer2]).unwrap();
        let baseline_allocations = heap.live_allocation_count();
        heap.fail_after_allocations(0);
        flags = [0; TAUT_NUM as usize];
        assert_eq!(
            CompareTwoPairsOfInChI(
                &mut heap,
                [SourceMutPointer::null(), tautomer1],
                [SourceMutPointer::null(), tautomer2],
                0,
                &mut flags,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(flags, [0, 0]);
        assert_eq!(heap.source_allocation_calls(), 2);
        assert_eq!(heap.live_allocation_count(), baseline_allocations);
    }

    fn stereo3_fixture(
        heap: &mut SourceHeap,
        centers: &[(AT_NUMB, i8)],
        inversion: i32,
        bonds: &[(AT_NUMB, AT_NUMB, i8)],
    ) -> INChI_Stereo {
        let center_numbers = (!centers.is_empty())
            .then(|| {
                heap.allocate_model_storage(centers.iter().map(|entry| entry.0).collect())
                    .unwrap()
            })
            .unwrap_or_default();
        let center_parities = (!centers.is_empty())
            .then(|| {
                heap.allocate_model_storage(centers.iter().map(|entry| entry.1).collect())
                    .unwrap()
            })
            .unwrap_or_default();
        let bond_atom1 = (!bonds.is_empty())
            .then(|| {
                heap.allocate_model_storage(bonds.iter().map(|entry| entry.0).collect())
                    .unwrap()
            })
            .unwrap_or_default();
        let bond_atom2 = (!bonds.is_empty())
            .then(|| {
                heap.allocate_model_storage(bonds.iter().map(|entry| entry.1).collect())
                    .unwrap()
            })
            .unwrap_or_default();
        let bond_parities = (!bonds.is_empty())
            .then(|| {
                heap.allocate_model_storage(bonds.iter().map(|entry| entry.2).collect())
                    .unwrap()
            })
            .unwrap_or_default();
        INChI_Stereo {
            nNumberOfStereoCenters: centers.len() as i32,
            nNumber: center_numbers,
            t_parity: center_parities,
            nCompInv2Abs: inversion,
            nNumberOfStereoBonds: bonds.len() as i32,
            nBondAtom1: bond_atom1,
            nBondAtom2: bond_atom2,
            b_parity: bond_parities,
            ..INChI_Stereo::default()
        }
    }

    struct CompareDisconnectedFixture {
        heap: SourceHeap,
        structures: [[SourceMutPointer<StrFromINChI>; TAUT_NUM as usize]; INCHI_NUM as usize],
        input: InpInChI,
        generated: SourceMutPointer<INChI>,
        original: SourceMutPointer<INChI>,
        reconnected: SourceMutPointer<INChI>,
        aux: SourceMutPointer<INChI_Aux>,
        rows: SourceMutPointer<[SourceMutPointer<INChI>; TAUT_NUM as usize]>,
        structure: SourceMutPointer<StrFromINChI>,
    }

    fn compare_disconnected_fixture() -> CompareDisconnectedFixture {
        let mut heap = SourceHeap::default();
        let generated_value = reversed_inchi3_fixture(&mut heap);
        let original_value = generated_value.clone();
        let generated = heap.allocate_model_storage(vec![generated_value]).unwrap();
        let original = heap.allocate_model_storage(vec![original_value]).unwrap();
        let reconnected = heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                ..INChI::default()
            }])
            .unwrap();
        let aux = heap
            .allocate_model_storage(vec![INChI_Aux::default()])
            .unwrap();
        let rows = heap
            .allocate_model_storage(vec![[SourceMutPointer::null(), generated]])
            .unwrap();
        let aux_rows = heap
            .allocate_model_storage(vec![[SourceMutPointer::null(), aux]])
            .unwrap();
        let mut structure_value = StrFromINChI {
            num_atoms: 1,
            ..StrFromINChI::default()
        };
        structure_value.RevInChI.pINChI[INCHI_BAS as usize] = rows;
        structure_value.RevInChI.pINChI_Aux[INCHI_BAS as usize] = aux_rows;
        structure_value.RevInChI.num_components[INCHI_BAS as usize] = 1;
        let structure = heap.allocate_model_storage(vec![structure_value]).unwrap();
        let mut structures = [[SourceMutPointer::null(); 2]; 2];
        structures[INCHI_REC as usize][TAUT_YES as usize] = structure;
        let mut input = InpInChI::default();
        input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize] = original;
        input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 1;
        input.pInpInChI[INCHI_REC as usize][TAUT_YES as usize] = reconnected;
        input.nNumComponents[INCHI_REC as usize][TAUT_YES as usize] = 1;
        CompareDisconnectedFixture {
            heap,
            structures,
            input,
            generated,
            original,
            reconnected,
            aux,
            rows,
            structure,
        }
    }

    #[test]
    fn source_port__ichirvr7__removefixhinchiidentical2mobh__line_256() {
        let mut empty_heap = SourceHeap::default();
        let mut negative = InpInChI::default();
        negative.nNumComponents[INCHI_BAS as usize] = [-1, i32::MAX];
        negative.nNumComponents[INCHI_REC as usize] = [i32::MIN, i32::MAX];
        assert_eq!(
            RemoveFixHInChIIdentical2MobH(&mut empty_heap, &mut negative),
            Ok(())
        );

        let mut invalid = InpInChI::default();
        invalid.nNumComponents[INCHI_BAS as usize] = [1, 1];
        assert_eq!(
            RemoveFixHInChIIdentical2MobH(&mut empty_heap, &mut invalid),
            Err(SourceHeapError::NullPointer)
        );

        let mut heap = SourceHeap::default();
        let mobile0 = reversed_inchi3_fixture(&mut heap);
        let fixed0 = reversed_inchi3_fixture(&mut heap);
        let freed_atom = fixed0.nAtom;
        let mobile1 = reversed_inchi3_fixture(&mut heap);
        let mut fixed1 = reversed_inchi3_fixture(&mut heap);
        fixed1.nTotalCharge = 1;
        let retained_atom = fixed1.nAtom;
        let mobile2 = reversed_inchi3_fixture(&mut heap);
        let fixed2 = reversed_inchi3_fixture(&mut heap);
        let beyond_min_atom = fixed2.nAtom;
        let mobile_bas = heap
            .allocate_model_storage(vec![mobile0, mobile1, mobile2])
            .unwrap();
        let fixed_bas = heap
            .allocate_model_storage(vec![fixed0, fixed1.clone(), fixed2.clone()])
            .unwrap();

        let mobile_rec_value = reversed_inchi3_fixture(&mut heap);
        let fixed_rec_value = reversed_inchi3_fixture(&mut heap);
        let freed_rec_atom = fixed_rec_value.nAtom;
        let mobile_rec = heap.allocate_model_storage(vec![mobile_rec_value]).unwrap();
        let fixed_rec = heap.allocate_model_storage(vec![fixed_rec_value]).unwrap();
        let mut input = InpInChI::default();
        input.pInpInChI[INCHI_BAS as usize] = [fixed_bas, mobile_bas];
        input.nNumComponents[INCHI_BAS as usize] = [2, 3];
        input.pInpInChI[INCHI_REC as usize] = [fixed_rec, mobile_rec];
        input.nNumComponents[INCHI_REC as usize] = [1, 1];
        assert_eq!(RemoveFixHInChIIdentical2MobH(&mut heap, &mut input), Ok(()));
        let fixed_values = heap.slice(fixed_bas.as_const()).unwrap();
        assert_eq!(fixed_values[0], INChI::default());
        assert_eq!(fixed_values[1], fixed1);
        assert_eq!(fixed_values[2], fixed2);
        assert_eq!(
            heap.slice(freed_atom.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert!(heap.slice(retained_atom.as_const()).is_ok());
        assert!(heap.slice(beyond_min_atom.as_const()).is_ok());
        assert_eq!(
            heap.slice(fixed_rec.as_const()).unwrap()[0],
            INChI::default()
        );
        assert_eq!(
            heap.slice(freed_rec_atom.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__ichirvr7__markdisconectedidenticaltoreconnected__line_283() {
        let mut heap = SourceHeap::default();
        let mut empty = InpInChI::default();
        empty.nNumComponents[INCHI_BAS as usize] = [i32::MIN, -1];
        empty.nNumComponents[INCHI_REC as usize] = [i32::MAX, i32::MIN];
        assert_eq!(
            MarkDisconectedIdenticalToReconnected(&mut heap, &mut empty),
            Ok(0)
        );

        let mut invalid = InpInChI::default();
        invalid.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 1;
        invalid.nNumComponents[INCHI_REC as usize][TAUT_YES as usize] = 1;
        assert_eq!(
            MarkDisconectedIdenticalToReconnected(&mut heap, &mut invalid),
            Err(SourceHeapError::NullPointer)
        );

        let disconnected0 = reversed_inchi3_fixture(&mut heap);
        let mut disconnected1 = reversed_inchi3_fixture(&mut heap);
        disconnected1.nTotalCharge = 1;
        let mut reconnected0 = reversed_inchi3_fixture(&mut heap);
        reconnected0.nTotalCharge = 1;
        let reconnected1 = reversed_inchi3_fixture(&mut heap);
        let disconnected_mobile = heap
            .allocate_model_storage(vec![disconnected0, disconnected1])
            .unwrap();
        let reconnected_mobile = heap
            .allocate_model_storage(vec![reconnected0, reconnected1])
            .unwrap();
        let mut crossed = InpInChI::default();
        crossed.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize] = disconnected_mobile;
        crossed.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 2;
        crossed.pInpInChI[INCHI_REC as usize][TAUT_YES as usize] = reconnected_mobile;
        crossed.nNumComponents[INCHI_REC as usize][TAUT_YES as usize] = 2;
        assert_eq!(
            MarkDisconectedIdenticalToReconnected(&mut heap, &mut crossed),
            Ok(2)
        );
        assert_eq!(
            heap.slice(disconnected_mobile.as_const())
                .unwrap()
                .iter()
                .map(|value| value.nLink)
                .collect::<Vec<_>>(),
            [-2, -1]
        );
        assert_eq!(
            heap.slice(reconnected_mobile.as_const())
                .unwrap()
                .iter()
                .map(|value| value.nLink)
                .collect::<Vec<_>>(),
            [2, 1]
        );
        assert_eq!(
            MarkDisconectedIdenticalToReconnected(&mut heap, &mut crossed),
            Ok(0)
        );

        let mut fixed_heap = SourceHeap::default();
        let disconnected_mobile_value = reversed_inchi3_fixture(&mut fixed_heap);
        let reconnected_mobile_value = reversed_inchi3_fixture(&mut fixed_heap);
        let disconnected_fixed_value = reversed_inchi3_fixture(&mut fixed_heap);
        let reconnected_fixed_value = reversed_inchi3_fixture(&mut fixed_heap);
        let disconnected_mobile = fixed_heap
            .allocate_model_storage(vec![disconnected_mobile_value])
            .unwrap();
        let reconnected_mobile = fixed_heap
            .allocate_model_storage(vec![reconnected_mobile_value])
            .unwrap();
        let disconnected_fixed = fixed_heap
            .allocate_model_storage(vec![disconnected_fixed_value])
            .unwrap();
        let reconnected_fixed = fixed_heap
            .allocate_model_storage(vec![reconnected_fixed_value])
            .unwrap();
        let mut fixed = InpInChI::default();
        fixed.pInpInChI[INCHI_BAS as usize] = [disconnected_fixed, disconnected_mobile];
        fixed.pInpInChI[INCHI_REC as usize] = [reconnected_fixed, reconnected_mobile];
        fixed.nNumComponents = [[[1_i32; 2]; 2][0], [[1_i32; 2]; 2][1]];
        assert_eq!(
            MarkDisconectedIdenticalToReconnected(&mut fixed_heap, &mut fixed),
            Ok(1)
        );
        assert_eq!(
            fixed_heap.slice(disconnected_mobile.as_const()).unwrap()[0].nLink,
            -1
        );
        assert_eq!(
            fixed_heap.slice(reconnected_mobile.as_const()).unwrap()[0].nLink,
            1
        );
        assert_eq!(
            fixed_heap.slice(disconnected_fixed.as_const()).unwrap()[0].nLink,
            -1
        );
        assert_eq!(
            fixed_heap.slice(reconnected_fixed.as_const()).unwrap()[0].nLink,
            1
        );

        let mut rejected_heap = SourceHeap::default();
        let mobile1_value = reversed_inchi3_fixture(&mut rejected_heap);
        let mobile2_value = reversed_inchi3_fixture(&mut rejected_heap);
        let fixed1_value = reversed_inchi3_fixture(&mut rejected_heap);
        let mut fixed2_value = reversed_inchi3_fixture(&mut rejected_heap);
        fixed2_value.nTotalCharge = 1;
        let mobile1 = rejected_heap
            .allocate_model_storage(vec![mobile1_value])
            .unwrap();
        let mobile2 = rejected_heap
            .allocate_model_storage(vec![mobile2_value])
            .unwrap();
        let fixed1 = rejected_heap
            .allocate_model_storage(vec![fixed1_value])
            .unwrap();
        let fixed2 = rejected_heap
            .allocate_model_storage(vec![fixed2_value])
            .unwrap();
        let mut rejected = InpInChI::default();
        rejected.pInpInChI[INCHI_BAS as usize] = [fixed1, mobile1];
        rejected.pInpInChI[INCHI_REC as usize] = [fixed2, mobile2];
        rejected.nNumComponents = [[1, 1], [1, 1]];
        assert_eq!(
            MarkDisconectedIdenticalToReconnected(&mut rejected_heap, &mut rejected),
            Ok(0)
        );
        rejected_heap.slice_mut(fixed2).unwrap()[0].bDeleted = 1;
        assert_eq!(
            MarkDisconectedIdenticalToReconnected(&mut rejected_heap, &mut rejected),
            Ok(0)
        );
        rejected.nNumComponents[INCHI_REC as usize][TAUT_NON as usize] = 0;
        assert_eq!(
            MarkDisconectedIdenticalToReconnected(&mut rejected_heap, &mut rejected),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichirvr7__addonemsg__line_3138() {
        let bytes = |value: &[u8]| value.iter().map(|byte| *byte as i8).collect::<Vec<_>>();

        let mut initial = [99_i8; 12];
        assert_eq!(
            AddOneMsg(&mut initial, 0, 12, &bytes(b"abc\0"), Some(&[b',' as i8]),),
            Ok(3)
        );
        assert_eq!(&initial[..4], &bytes(b"abc\0"));
        assert!(initial[4..].iter().all(|byte| *byte == 99));

        let mut complete = [99_i8; 10];
        complete[..2].copy_from_slice(&bytes(b"a\0"));
        assert_eq!(
            AddOneMsg(&mut complete, 1, 6, &bytes(b"bc\0"), Some(&bytes(b", \0")),),
            Ok(5)
        );
        assert_eq!(&complete[..6], &bytes(b"a, bc\0"));
        assert!(complete[6..].iter().all(|byte| *byte == 99));

        let mut exact = [99_i8; 10];
        exact[..2].copy_from_slice(&bytes(b"a\0"));
        let exact_before = exact;
        assert_eq!(
            AddOneMsg(&mut exact, 1, 5, &bytes(b"bc\0"), Some(&bytes(b", \0")),),
            Ok(1)
        );
        assert_eq!(exact, exact_before);

        let long = bytes(b"abcdefghijklmnopqrstuvwxyz\0");
        let mut threshold = [99_i8; 24];
        threshold[..2].copy_from_slice(&bytes(b"a\0"));
        let threshold_before = threshold;
        assert_eq!(
            AddOneMsg(&mut threshold, 1, 17, &long, Some(&bytes(b", \0")),),
            Ok(1)
        );
        assert_eq!(threshold, threshold_before);

        let mut truncated = [99_i8; 24];
        truncated[..2].copy_from_slice(&bytes(b"a\0"));
        assert_eq!(
            AddOneMsg(&mut truncated, 1, 18, &long, Some(&bytes(b", \0")),),
            Ok(17)
        );
        assert_eq!(&truncated[..18], &bytes(b"a, abcdefghijk...\0"));
        assert!(truncated[18..].iter().all(|byte| *byte == 99));

        let mut negative_total = [99_i8; 8];
        let negative_before = negative_total;
        assert_eq!(
            AddOneMsg(&mut negative_total, 0, -1, &bytes(b"abc\0"), None),
            Ok(0)
        );
        assert_eq!(negative_total, negative_before);
        assert_eq!(
            AddOneMsg(&mut [99_i8; 8], -1, 8, &bytes(b"abc\0"), None),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            AddOneMsg(&mut [99_i8; 8], 0, 8, &bytes(b"abc"), None),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut partial = [99_i8; 4];
        partial[..2].copy_from_slice(&bytes(b"a\0"));
        assert_eq!(
            AddOneMsg(
                &mut partial,
                1,
                20,
                &bytes(b"long message\0"),
                Some(&bytes(b", \0")),
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(partial, [b'a' as i8, b',' as i8, b' ' as i8, 0]);
    }

    #[test]
    fn source_port__ichirvr7__filloutcomparemessage__line_3179() {
        let text = |buffer: &[i8]| {
            let length = buffer.iter().position(|byte| *byte == 0).unwrap();
            String::from_utf8(buffer[..length].iter().map(|byte| *byte as u8).collect()).unwrap()
        };

        let mut unchanged = [99_i8; 32];
        unchanged[0] = 0;
        let before = unchanged;
        assert_eq!(FillOutCompareMessage(&mut unchanged, 32, &[0, 0]), Ok(0));
        assert_eq!(unchanged, before);

        let cases = [
            (1_u64, " Error: Wrong result"),
            (2, " Hydrogens: Locations or number"),
            (4, " Hydrogens: Fixed-H"),
            (8, " Hydrogens: Number"),
            (16, " Mobile-H groups: Missing"),
            (32, " Mobile-H groups: Falsely present"),
            (64, " Mobile-H groups: One instead of multiple"),
            (128, " Mobile-H groups: Multiple instead of one"),
            (256, " Mobile-H groups: Attachment points"),
            (512, " Mobile-H groups: Number"),
            (1024, " Isotopic: Atoms do not match"),
            (
                2048,
                " Exchangeable isotopic H: Does not match for a component",
            ),
            (4096, " Exchangeable isotopic H: Do not match"),
            (8192, " Charge(s): Do not match"),
            (16384, " Proton balance: Does not match for a component"),
            (32768, " Proton balance: Does not match"),
            (65536, " Stereo centers/allenes: Falsely inverted"),
            (131072, " Stereo centers/allenes: Wrong parity"),
            (262144, " Stereo centers/allenes: Extra undefined"),
            (524288, " Stereo centers/allenes: Extra known"),
            (1048576, " Stereo centers/allenes: Missing undefined"),
            (2097152, " Stereo centers/allenes: Missing known"),
            (4194304, " Stereobonds/cumulenes: Wrong parity"),
            (8388608, " Stereobonds/cumulenes: Extra undefined"),
            (16777216, " Stereobonds/cumulenes: Missing known"),
            (33554432, " Stereobonds/cumulenes: Missing undefined"),
            (67108864, " Stereobonds/cumulenes: Missing known"),
            (134217728, " Fixed-H layer: Missing or extra"),
            (268435456, " Number of components: Does not match"),
            (536870912, " Conversion encountered: Error"),
        ];
        for (bit, expected) in cases {
            let mut buffer = [99_i8; 512];
            buffer[0] = 0;
            assert_eq!(FillOutCompareMessage(&mut buffer, 512, &[0, bit]), Ok(-1));
            assert_eq!(
                text(&buffer),
                format!(" Problems/mismatches: Mobile-H({expected})")
            );
            let nul = buffer.iter().position(|byte| *byte == 0).unwrap();
            assert!(buffer[nul + 1..].iter().all(|byte| *byte == 99));
        }

        let mut grouped = [99_i8; 512];
        grouped[0] = 0;
        let grouped_bits = 2 | 4 | 8 | 16 | 8192;
        assert_eq!(
            FillOutCompareMessage(&mut grouped, 512, &[grouped_bits, grouped_bits]),
            Ok(-1)
        );
        assert_eq!(
            text(&grouped),
            " Problems/mismatches: Mobile-H( Hydrogens: Locations or number, Fixed-H, Number; Mobile-H groups: Missing; Charge(s): Do not match) Fixed-H( Hydrogens: Locations or number, Fixed-H, Number; Mobile-H groups: Missing; Charge(s): Do not match)"
        );

        let mut existing = [0_i8; 256];
        let prefix = b"prefix Problems/mismatches: suffix";
        for (target, source) in existing.iter_mut().zip(prefix) {
            *target = *source as i8;
        }
        assert_eq!(FillOutCompareMessage(&mut existing, 256, &[0, 1]), Ok(-1));
        assert_eq!(
            text(&existing),
            "prefix Problems/mismatches: suffix Mobile-H( Error: Wrong result)"
        );

        let mut unknown = [99_i8; 128];
        unknown[0] = 0;
        assert_eq!(
            FillOutCompareMessage(&mut unknown, 128, &[0, 1_u64 << 63]),
            Ok(-1)
        );
        assert_eq!(text(&unknown), " Problems/mismatches: Mobile-H()");

        let mut truncated = [99_i8; 64];
        truncated[0] = 0;
        assert_eq!(FillOutCompareMessage(&mut truncated, 5, &[0, 1]), Ok(-1));
        assert_eq!(text(&truncated), ")");
        assert_eq!(
            FillOutCompareMessage(&mut [99_i8; 8], 8, &[0, 0]),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichirvr7__setupsrm__line_364() {
        let mut restore_mode = SRM {
            bMetalAddFlower: -1,
            nMetalMinBondOrder: -2,
            nMetalInitEdgeFlow: -3,
            nMetalInitBondOrder: -4,
            nMetal2EndpointMinBondOrder: -5,
            nMetal2EndpointInitBondOrder: -6,
            nMetal2EndpointInitEdgeFlow: -7,
            nMetalFlowerParam_D: -8,
            nMetalMaxCharge_D: -9,
            bStereoRemovesMetalFlag: -10,
            bFixStereoBonds: -11,
        };
        SetUpSrm(&mut restore_mode);
        assert_eq!(
            restore_mode,
            SRM {
                bMetalAddFlower: 1,
                nMetalMinBondOrder: 0,
                nMetalInitEdgeFlow: 1,
                nMetalInitBondOrder: 1,
                nMetal2EndpointMinBondOrder: 1,
                nMetal2EndpointInitBondOrder: 1,
                nMetal2EndpointInitEdgeFlow: 0,
                nMetalFlowerParam_D: 16,
                nMetalMaxCharge_D: 16,
                bStereoRemovesMetalFlag: 0,
                bFixStereoBonds: 0,
            }
        );
    }

    #[test]
    fn source_port__ichirvr7__mergestructurecomponents__line_406() {
        fn call_merge(
            heap: &mut SourceHeap,
            structures: &[[SourceMutPointer<StrFromINChI>; 2]; 2],
            input: &mut InpInChI,
        ) -> Result<i32, SourceHeapError> {
            MergeStructureComponents(
                heap,
                &INPUT_PARMS::default(),
                &mut STRUCT_DATA::default(),
                i64::MIN,
                SourceMutPointer::null(),
                &SRM::default(),
                i32::MAX,
                structures,
                input,
            )
        }

        fn single_branch(representation: usize, mobile_h: usize, marker: u8) {
            let mut heap = SourceHeap::default();
            let atom = heap
                .allocate_model_storage(vec![inp_ATOM {
                    el_number: marker,
                    orig_at_number: 1,
                    ..inp_ATOM::default()
                }])
                .unwrap();
            let structure = heap
                .allocate_model_storage(vec![StrFromINChI {
                    at2: atom,
                    num_atoms: 1,
                    ..StrFromINChI::default()
                }])
                .unwrap();
            let mut structures = [[SourceMutPointer::null(); 2]; 2];
            structures[representation][mobile_h] = structure;
            let mut input = InpInChI::default();
            input.nNumComponents[representation][mobile_h] = 1;
            assert_eq!(call_merge(&mut heap, &structures, &mut input), Ok(0));
            assert_eq!(input.num_atoms, 1);
            assert_eq!(
                heap.slice(input.atom.as_const()).unwrap()[0].el_number,
                marker
            );
        }

        let mut empty_heap = SourceHeap::default();
        let old_empty = empty_heap
            .allocate_model_storage(vec![inp_ATOM {
                el_number: 99,
                ..inp_ATOM::default()
            }])
            .unwrap();
        empty_heap.trace_source_allocations();
        let mut empty_input = InpInChI {
            atom: old_empty,
            num_atoms: 17,
            ..InpInChI::default()
        };
        assert_eq!(
            call_merge(
                &mut empty_heap,
                &[[SourceMutPointer::null(); 2]; 2],
                &mut empty_input,
            ),
            Ok(0)
        );
        assert_eq!((empty_input.atom, empty_input.num_atoms), (old_empty, 0));
        assert_eq!(empty_heap.source_allocation_calls(), 0);

        single_branch(INCHI_REC as usize, TAUT_NON as usize, 11);
        single_branch(INCHI_REC as usize, TAUT_YES as usize, 12);
        single_branch(INCHI_BAS as usize, TAUT_NON as usize, 13);
        single_branch(INCHI_BAS as usize, TAUT_YES as usize, 14);

        let mut priority_heap = SourceHeap::default();
        let mut priority_structures = [[SourceMutPointer::null(); 2]; 2];
        let mut priority_input = InpInChI::default();
        for (representation, mobile_h, marker) in [
            (INCHI_BAS as usize, TAUT_YES as usize, 21_u8),
            (INCHI_BAS as usize, TAUT_NON as usize, 22_u8),
            (INCHI_REC as usize, TAUT_YES as usize, 23_u8),
            (INCHI_REC as usize, TAUT_NON as usize, 24_u8),
        ] {
            let atom = priority_heap
                .allocate_model_storage(vec![inp_ATOM {
                    el_number: marker,
                    ..inp_ATOM::default()
                }])
                .unwrap();
            priority_structures[representation][mobile_h] = priority_heap
                .allocate_model_storage(vec![StrFromINChI {
                    at2: atom,
                    num_atoms: 1,
                    ..StrFromINChI::default()
                }])
                .unwrap();
            priority_input.nNumComponents[representation][mobile_h] = 1;
        }
        assert_eq!(
            call_merge(
                &mut priority_heap,
                &priority_structures,
                &mut priority_input,
            ),
            Ok(0)
        );
        assert_eq!(
            priority_heap.slice(priority_input.atom.as_const()).unwrap()[0].el_number,
            24
        );

        let mut fallback_heap = SourceHeap::default();
        let fallback_atom = fallback_heap
            .allocate_model_storage(vec![inp_ATOM {
                el_number: 31,
                ..inp_ATOM::default()
            }])
            .unwrap();
        let fixed = fallback_heap
            .allocate_model_storage(vec![StrFromINChI::default()])
            .unwrap();
        let mobile = fallback_heap
            .allocate_model_storage(vec![StrFromINChI {
                at2: fallback_atom,
                num_atoms: 1,
                ..StrFromINChI::default()
            }])
            .unwrap();
        let mut fallback_structures = [[SourceMutPointer::null(); 2]; 2];
        fallback_structures[INCHI_REC as usize][TAUT_NON as usize] = fixed;
        fallback_structures[INCHI_REC as usize][TAUT_YES as usize] = mobile;
        let mut fallback_input = InpInChI::default();
        fallback_input.nNumComponents[INCHI_REC as usize][TAUT_NON as usize] = 1;
        fallback_input.nNumComponents[INCHI_REC as usize][TAUT_YES as usize] = 1;
        assert_eq!(
            call_merge(
                &mut fallback_heap,
                &fallback_structures,
                &mut fallback_input,
            ),
            Ok(0)
        );
        assert_eq!(
            fallback_heap.slice(fallback_input.atom.as_const()).unwrap()[0].el_number,
            31
        );

        let mut merge_heap = SourceHeap::default();
        let mut first = inp_ATOM {
            el_number: 41,
            valence: 2,
            orig_at_number: 1,
            endpoint: 10,
            bAmbiguousStereo: 11,
            at_type: 12,
            bCutVertex: 13,
            bUsed0DParity: 14,
            cFlags: 15,
            nBlockSystem: 16,
            nNumAtInRingSystem: 17,
            nRingSystem: 18,
            p_parity: 1,
            p_orig_at_num: [1, 2, 3, AT_NUMB::MAX],
            sb_parity: [1, 1, 0],
            sn_orig_at_num: [2, 3, 7],
            ..inp_ATOM::default()
        };
        first.neighbor[0] = 1;
        first.neighbor[1] = 2;
        let second = inp_ATOM {
            el_number: 42,
            valence: -1,
            orig_at_number: AT_NUMB::MAX,
            p_orig_at_num: [9, 8, 7, 6],
            sn_orig_at_num: [5, 4, 3],
            ..inp_ATOM::default()
        };
        let mut first_h = inp_ATOM {
            el_number: 1,
            orig_at_number: 3,
            ..inp_ATOM::default()
        };
        first_h.neighbor[0] = 0;
        let first_atoms = merge_heap
            .allocate_model_storage(vec![first, second, first_h])
            .unwrap();

        let mut third = inp_ATOM {
            el_number: 43,
            valence: 2,
            orig_at_number: AT_NUMB::MAX,
            p_parity: 1,
            p_orig_at_num: [1, 2, 0, AT_NUMB::MAX],
            sb_parity: [1, 0, 1],
            sn_orig_at_num: [1, 2, 9],
            ..inp_ATOM::default()
        };
        third.neighbor[0] = 0;
        third.neighbor[1] = 1;
        let mut second_h = inp_ATOM {
            el_number: 1,
            orig_at_number: 2,
            ..inp_ATOM::default()
        };
        second_h.neighbor[0] = 0;
        let second_atoms = merge_heap
            .allocate_model_storage(vec![third, second_h])
            .unwrap();
        let ignored_atom = merge_heap
            .allocate_model_storage(vec![inp_ATOM {
                el_number: 88,
                ..inp_ATOM::default()
            }])
            .unwrap();
        let components = merge_heap
            .allocate_model_storage(vec![
                StrFromINChI {
                    at2: first_atoms,
                    num_atoms: 2,
                    num_deleted_H: 1,
                    ..StrFromINChI::default()
                },
                StrFromINChI {
                    at2: second_atoms,
                    num_atoms: 1,
                    num_deleted_H: 1,
                    ..StrFromINChI::default()
                },
                StrFromINChI {
                    at2: ignored_atom,
                    num_atoms: 1,
                    bDeleted: 1,
                    ..StrFromINChI::default()
                },
                StrFromINChI {
                    num_atoms: 1,
                    ..StrFromINChI::default()
                },
            ])
            .unwrap();
        let mut merge_structures = [[SourceMutPointer::null(); 2]; 2];
        merge_structures[INCHI_BAS as usize][TAUT_YES as usize] = components;
        let mut merge_input = InpInChI::default();
        merge_input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 4;
        assert_eq!(
            call_merge(&mut merge_heap, &merge_structures, &mut merge_input),
            Ok(0)
        );
        assert_eq!(merge_input.num_atoms, 5);
        let merged = merge_heap.slice(merge_input.atom.as_const()).unwrap();
        assert_eq!(
            merged[..5]
                .iter()
                .map(|atom| atom.el_number)
                .collect::<Vec<_>>(),
            [41, 42, 43, 1, 1]
        );
        assert_eq!((merged[0].neighbor[0], merged[0].neighbor[1]), (1, 3));
        assert_eq!(merged[0].p_orig_at_num, [1, 2, 4, 0]);
        assert_eq!(merged[0].sn_orig_at_num, [2, 4, 7]);
        assert_eq!((merged[0].orig_at_number, merged[0].component), (1, 1));
        assert_eq!(
            (
                merged[0].endpoint,
                merged[0].bAmbiguousStereo,
                merged[0].at_type,
                merged[0].bCutVertex,
                merged[0].bUsed0DParity,
                merged[0].cFlags,
                merged[0].nBlockSystem,
                merged[0].nNumAtInRingSystem,
                merged[0].nRingSystem,
            ),
            (0, 0, 0, 0, 0, 0, 0, 0, 0)
        );
        assert_eq!(merged[1].orig_at_number, AT_NUMB::MAX);
        assert_eq!(merged[1].p_orig_at_num, [9, 8, 7, 6]);
        assert_eq!(merged[1].sn_orig_at_num, [5, 4, 3]);
        assert_eq!((merged[2].neighbor[0], merged[2].neighbor[1]), (2, 4));
        assert_eq!((merged[2].orig_at_number, merged[2].component), (1, 2));
        assert_eq!(merged[2].p_orig_at_num, [3, 5, 2, 2]);
        assert_eq!(merged[2].sn_orig_at_num, [3, 2, 9]);
        assert_eq!((merged[3].neighbor[0], merged[3].orig_at_number), (0, 4));
        assert_eq!((merged[4].neighbor[0], merged[4].orig_at_number), (2, 5));

        let mut zero_heap = SourceHeap::default();
        let zero_atom = zero_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let zero_structure = zero_heap
            .allocate_model_storage(vec![StrFromINChI {
                at2: zero_atom,
                num_atoms: 1,
                bDeleted: 1,
                ..StrFromINChI::default()
            }])
            .unwrap();
        let old_zero = zero_heap
            .allocate_model_storage(vec![inp_ATOM {
                el_number: 77,
                ..inp_ATOM::default()
            }])
            .unwrap();
        let mut zero_structures = [[SourceMutPointer::null(); 2]; 2];
        zero_structures[INCHI_BAS as usize][TAUT_YES as usize] = zero_structure;
        let mut zero_input = InpInChI {
            atom: old_zero,
            num_atoms: 9,
            ..InpInChI::default()
        };
        zero_input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 1;
        let zero_baseline = zero_heap.live_allocation_count();
        zero_heap.trace_source_allocations();
        assert_eq!(
            call_merge(&mut zero_heap, &zero_structures, &mut zero_input),
            Ok(0)
        );
        assert_eq!((zero_input.atom, zero_input.num_atoms), (old_zero, 0));
        assert_eq!(zero_heap.source_allocation_calls(), 3);
        assert_eq!(zero_heap.live_allocation_count(), zero_baseline);

        for failure_after in 0..3_u64 {
            let mut heap = SourceHeap::default();
            let atom = heap
                .allocate_model_storage(vec![inp_ATOM::default()])
                .unwrap();
            let structure = heap
                .allocate_model_storage(vec![StrFromINChI {
                    at2: atom,
                    num_atoms: 1,
                    ..StrFromINChI::default()
                }])
                .unwrap();
            let old = heap
                .allocate_model_storage(vec![inp_ATOM {
                    el_number: 66,
                    ..inp_ATOM::default()
                }])
                .unwrap();
            let mut structures = [[SourceMutPointer::null(); 2]; 2];
            structures[INCHI_REC as usize][TAUT_YES as usize] = structure;
            let mut input = InpInChI {
                atom: old,
                num_atoms: 8,
                ..InpInChI::default()
            };
            input.nNumComponents[INCHI_REC as usize][TAUT_YES as usize] = 1;
            let baseline = heap.live_allocation_count();
            heap.fail_after_allocations(failure_after);
            assert_eq!(
                call_merge(&mut heap, &structures, &mut input),
                Ok(RI_ERR_ALLOC)
            );
            assert_eq!((input.atom, input.num_atoms), (old, 0));
            assert_eq!(
                heap.source_allocation_calls(),
                if failure_after < 2 { 2 } else { 3 }
            );
            assert_eq!(heap.live_allocation_count(), baseline);
        }

        let mut malformed_heap = SourceHeap::default();
        let malformed_atom = malformed_heap
            .allocate_model_storage(vec![inp_ATOM {
                valence: 21,
                ..inp_ATOM::default()
            }])
            .unwrap();
        let malformed_structure = malformed_heap
            .allocate_model_storage(vec![StrFromINChI {
                at2: malformed_atom,
                num_atoms: 1,
                ..StrFromINChI::default()
            }])
            .unwrap();
        let old_malformed = malformed_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let mut malformed_structures = [[SourceMutPointer::null(); 2]; 2];
        malformed_structures[INCHI_REC as usize][TAUT_YES as usize] = malformed_structure;
        let mut malformed_input = InpInChI {
            atom: old_malformed,
            num_atoms: 3,
            ..InpInChI::default()
        };
        malformed_input.nNumComponents[INCHI_REC as usize][TAUT_YES as usize] = 1;
        let malformed_baseline = malformed_heap.live_allocation_count();
        assert_eq!(
            call_merge(
                &mut malformed_heap,
                &malformed_structures,
                &mut malformed_input,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            (malformed_input.atom, malformed_input.num_atoms),
            (old_malformed, 0)
        );
        assert_eq!(malformed_heap.live_allocation_count(), malformed_baseline);

        let mut negative_heap = SourceHeap::default();
        let old_negative = negative_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let mut negative_input = InpInChI {
            atom: old_negative,
            num_atoms: 2,
            ..InpInChI::default()
        };
        negative_input.nNumComponents[INCHI_REC as usize][TAUT_NON as usize] = -1;
        assert_eq!(
            call_merge(
                &mut negative_heap,
                &[[SourceMutPointer::null(); 2]; 2],
                &mut negative_input,
            ),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(
            (negative_input.atom, negative_input.num_atoms),
            (old_negative, 0)
        );
    }

    fn compare_disconnected(
        fixture: &mut CompareDisconnectedFixture,
        fixed: i32,
    ) -> Result<i32, SourceHeapError> {
        CompareAllDisconnectedOrigInchiToRevInChI(
            &mut fixture.heap,
            &fixture.structures,
            &mut fixture.input,
            fixed,
            i64::MIN,
            SourceMutPointer::null(),
        )
    }

    #[test]
    fn source_port__ichirvr7__comparealldisconnectedoriginchitorevinchi__line_1683() {
        let mut empty_heap = SourceHeap::default();
        let empty_structures = [[SourceMutPointer::null(); 2]; 2];
        let mut empty_input = InpInChI::default();
        empty_input.CompareInchiFlags = [[INCHI_MODE::MAX; 2]; 2];
        assert_eq!(
            CompareAllDisconnectedOrigInchiToRevInChI(
                &mut empty_heap,
                &empty_structures,
                &mut empty_input,
                i32::MIN,
                i64::MAX,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(empty_input.CompareInchiFlags[0], [INCHI_MODE::MAX; 2]);
        assert_eq!(empty_input.CompareInchiFlags[1], [0, 0]);

        let mut no_fragments = compare_disconnected_fixture();
        no_fragments.heap.slice_mut(no_fragments.original).unwrap()[0].nLink = -1;
        no_fragments
            .heap
            .slice_mut(no_fragments.reconnected)
            .unwrap()[0]
            .nLink = 1;
        no_fragments.heap.slice_mut(no_fragments.structure).unwrap()[0]
            .RevInChI
            .num_components[INCHI_BAS as usize] = 0;
        assert_eq!(compare_disconnected(&mut no_fragments, 0), Ok(0));

        let mut count_mismatch = compare_disconnected_fixture();
        count_mismatch
            .heap
            .slice_mut(count_mismatch.reconnected)
            .unwrap()[0]
            .nLink = 1;
        count_mismatch
            .heap
            .slice_mut(count_mismatch.structure)
            .unwrap()[0]
            .RevInChI
            .num_components[INCHI_BAS as usize] = 0;
        assert_eq!(compare_disconnected(&mut count_mismatch, 0), Ok(1));
        assert_eq!(
            count_mismatch.input.CompareInchiFlags[1],
            [0, tagInchiCompareDiffBits_INCHIDIFF_PROBLEM as INCHI_MODE]
        );

        let mut equal = compare_disconnected_fixture();
        assert_eq!(compare_disconnected(&mut equal, 0), Ok(0));
        assert_eq!(equal.input.CompareInchiFlags[1], [0, 0]);

        let mut changed = compare_disconnected_fixture();
        changed.heap.slice_mut(changed.generated).unwrap()[0].nNumberOfAtoms = 3;
        assert_eq!(compare_disconnected(&mut changed, 0), Ok(0));
        assert_eq!(
            changed.input.CompareInchiFlags[1],
            [0, tagInchiCompareDiffBits_INCHIDIFF_NUM_AT as INCHI_MODE]
        );

        let mut proton_difference = compare_disconnected_fixture();
        proton_difference
            .heap
            .slice_mut(proton_difference.aux)
            .unwrap()[0]
            .nNumRemovedProtons = 3;
        proton_difference
            .heap
            .slice_mut(proton_difference.aux)
            .unwrap()[0]
            .nNumRemovedIsotopicH = [1, 0, 2];
        assert_eq!(compare_disconnected(&mut proton_difference, 0), Ok(0));
        assert_eq!(
            proton_difference.input.CompareInchiFlags[1],
            [
                0,
                (tagInchiCompareDiffBits_INCHIDIFF_MOBH_PROTONS
                    | tagInchiCompareDiffBits_INCHIDIFF_REM_ISO_H) as INCHI_MODE,
            ]
        );

        let mut deleted_last = compare_disconnected_fixture();
        deleted_last.input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize] =
            SourceMutPointer::null();
        deleted_last.input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 0;
        deleted_last
            .heap
            .slice_mut(deleted_last.reconnected)
            .unwrap()[0]
            .bDeleted = 1;
        deleted_last.heap.slice_mut(deleted_last.generated).unwrap()[0].bDeleted = 1;
        deleted_last.heap.slice_mut(deleted_last.aux).unwrap()[0].nNumRemovedProtons = 5;
        assert_eq!(compare_disconnected(&mut deleted_last, 0), Ok(0));
        assert_eq!(deleted_last.input.CompareInchiFlags[1], [0, 0]);

        let mut fixed = compare_disconnected_fixture();
        let mut fixed_value = fixed.heap.slice(fixed.generated.as_const()).unwrap()[0].clone();
        fixed_value.nTotalCharge = 1;
        let fixed_generated = fixed
            .heap
            .allocate_model_storage(vec![fixed_value.clone()])
            .unwrap();
        let fixed_original = fixed
            .heap
            .allocate_model_storage(vec![fixed_value])
            .unwrap();
        let fixed_reconnected = fixed
            .heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                ..INChI::default()
            }])
            .unwrap();
        fixed.heap.slice_mut(fixed.rows).unwrap()[0][TAUT_NON as usize] = fixed_generated;
        fixed.input.pInpInChI[INCHI_BAS as usize][TAUT_NON as usize] = fixed_original;
        fixed.input.nNumComponents[INCHI_BAS as usize][TAUT_NON as usize] = 1;
        fixed.input.pInpInChI[INCHI_REC as usize][TAUT_NON as usize] = fixed_reconnected;
        fixed.input.nNumComponents[INCHI_REC as usize][TAUT_NON as usize] = 1;
        fixed.structures[INCHI_REC as usize][TAUT_NON as usize] = fixed.structure;
        assert_eq!(compare_disconnected(&mut fixed, 1), Ok(0));
        assert_eq!(fixed.input.CompareInchiFlags[1], [0, 0]);

        let original_protons = fixed
            .heap
            .allocate_model_storage(vec![COMPONENT_REM_PROTONS {
                nNumRemovedProtons: 1,
                nNumRemovedIsotopicH: [1, 2, 3],
            }])
            .unwrap();
        fixed.input.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize].pNumProtons =
            original_protons;
        fixed.heap.slice_mut(fixed.aux).unwrap()[0].nNumRemovedProtons = 2;
        fixed.heap.slice_mut(fixed.aux).unwrap()[0].nNumRemovedIsotopicH = [1, 9, 3];
        assert_eq!(compare_disconnected(&mut fixed, 1), Ok(0));
        let proton_bits = (tagInchiCompareDiffBits_INCHIDIFF_MOBH_PROTONS
            | tagInchiCompareDiffBits_INCHIDIFF_REM_ISO_H) as INCHI_MODE;
        assert_eq!(
            fixed.input.CompareInchiFlags[1][TAUT_YES as usize] & proton_bits,
            proton_bits
        );

        let mut allocation = compare_disconnected_fixture();
        let baseline_allocations = allocation.heap.live_allocation_count();
        allocation.heap.fail_after_allocations(0);
        assert_eq!(compare_disconnected(&mut allocation, 0), Ok(RI_ERR_ALLOC));
        assert_eq!(allocation.heap.source_allocation_calls(), 2);
        assert_eq!(
            allocation.heap.live_allocation_count(),
            baseline_allocations
        );
    }

    struct CompareOneFixture {
        structure: StrFromINChI,
        input: [SourceMutPointer<INChI>; TAUT_NUM as usize],
        generated: Vec<[SourceMutPointer<INChI>; TAUT_NUM as usize]>,
        aux: Vec<[SourceMutPointer<INChI_Aux>; TAUT_NUM as usize]>,
    }

    fn compare_one_fixture(heap: &mut SourceHeap, components: usize) -> CompareOneFixture {
        let mut generated = Vec::new();
        let mut aux = Vec::new();
        for index in 0..components {
            let fixed = reversed_inchi3_fixture(heap);
            let mobile = reversed_inchi3_fixture(heap);
            generated.push([
                heap.allocate_model_storage(vec![fixed]).unwrap(),
                heap.allocate_model_storage(vec![mobile]).unwrap(),
            ]);
            let removed = INChI_Aux {
                nNumRemovedProtons: index as i16 + 1,
                nNumRemovedIsotopicH: [index as i16 + 2, index as i16 + 3, index as i16 + 4],
                ..INChI_Aux::default()
            };
            aux.push([
                heap.allocate_model_storage(vec![INChI_Aux::default()])
                    .unwrap(),
                heap.allocate_model_storage(vec![removed]).unwrap(),
            ]);
        }
        let input_fixed = reversed_inchi3_fixture(heap);
        let input_mobile = reversed_inchi3_fixture(heap);
        let input = [
            heap.allocate_model_storage(vec![input_fixed]).unwrap(),
            heap.allocate_model_storage(vec![input_mobile]).unwrap(),
        ];
        let rows = if components == 0 {
            heap.allocate_model_storage(Vec::<[SourceMutPointer<INChI>; 2]>::new())
                .unwrap()
        } else {
            heap.allocate_model_storage(generated.clone()).unwrap()
        };
        let aux_rows = if components == 0 {
            heap.allocate_model_storage(Vec::<[SourceMutPointer<INChI_Aux>; 2]>::new())
                .unwrap()
        } else {
            heap.allocate_model_storage(aux.clone()).unwrap()
        };
        let mut structure = StrFromINChI {
            iINCHI: INCHI_BAS as i8,
            iMobileH: TAUT_YES as i8,
            RevInChI: crate::source_types::REV_INCHI {
                nRetVal: _IS_OKAY as i32,
                ..crate::source_types::REV_INCHI::default()
            },
            ..StrFromINChI::default()
        };
        structure.RevInChI.pINChI[INCHI_BAS as usize] = rows;
        structure.RevInChI.pINChI_Aux[INCHI_BAS as usize] = aux_rows;
        structure.RevInChI.num_components[INCHI_BAS as usize] = components as i32;
        CompareOneFixture {
            structure,
            input,
            generated,
            aux,
        }
    }

    #[test]
    fn source_port__ichirvr7__compareoneoriginchitorevinchi__line_2240() {
        fn call(
            heap: &mut SourceHeap,
            structure: Option<&StrFromINChI>,
            input: [SourceMutPointer<INChI>; 2],
            mobile_h: i32,
            removed: &mut COMPONENT_REM_PROTONS,
            flags: &mut [INCHI_MODE; 2],
        ) -> Result<i32, SourceHeapError> {
            CompareOneOrigInchiToRevInChI(
                heap,
                structure,
                input,
                mobile_h,
                17,
                23,
                SourceMutPointer::null(),
                removed,
                flags,
            )
        }

        let mut heap = SourceHeap::default();
        let mut removed = COMPONENT_REM_PROTONS {
            nNumRemovedProtons: 99,
            nNumRemovedIsotopicH: [99; 3],
        };
        let mut flags = [0; 2];
        assert_eq!(
            call(
                &mut heap,
                None,
                [SourceMutPointer::null(); 2],
                TAUT_YES as i32,
                &mut removed,
                &mut flags,
            ),
            Ok(0)
        );
        assert_eq!(
            flags[TAUT_YES as usize],
            tagInchiCompareDiffBits_INCHIDIFF_STR2INCHI_ERR as INCHI_MODE
        );
        assert_eq!(removed.nNumRemovedProtons, 99);

        let mut fixture = compare_one_fixture(&mut heap, 1);
        let mut error_structure = fixture.structure.clone();
        error_structure.RevInChI.nRetVal = crate::source_types::_IS_ERROR as i32;
        flags = [0; 2];
        assert_eq!(
            call(
                &mut heap,
                Some(&error_structure),
                fixture.input,
                TAUT_NON as i32,
                &mut removed,
                &mut flags,
            ),
            Ok(0)
        );
        assert_eq!(
            flags[TAUT_NON as usize],
            tagInchiCompareDiffBits_INCHIDIFF_STR2INCHI_ERR as INCHI_MODE
        );

        flags = [0; 2];
        removed = COMPONENT_REM_PROTONS::default();
        assert_eq!(
            call(
                &mut heap,
                Some(&fixture.structure),
                fixture.input,
                TAUT_YES as i32,
                &mut removed,
                &mut flags,
            ),
            Ok(0)
        );
        assert_eq!(flags, [0, 0]);
        assert_eq!(
            removed,
            COMPONENT_REM_PROTONS {
                nNumRemovedProtons: 1,
                nNumRemovedIsotopicH: [2, 3, 4],
            }
        );

        heap.slice_mut(fixture.input[0]).unwrap()[0].nAtom =
            heap.allocate_model_storage(vec![6_u8, 7]).unwrap();
        flags = [0; 2];
        call(
            &mut heap,
            Some(&fixture.structure),
            fixture.input,
            TAUT_YES as i32,
            &mut removed,
            &mut flags,
        )
        .unwrap();
        assert_eq!(
            flags[TAUT_YES as usize],
            tagInchiCompareDiffBits_INCHIDIFF_PROBLEM as INCHI_MODE
        );

        let mut deleted = compare_one_fixture(&mut heap, 1);
        deleted.structure.bDeleted = 1;
        heap.slice_mut(deleted.input[0]).unwrap()[0].bDeleted = 1;
        flags = [0; 2];
        removed.nNumRemovedProtons = 77;
        assert_eq!(
            call(
                &mut heap,
                Some(&deleted.structure),
                deleted.input,
                TAUT_YES as i32,
                &mut removed,
                &mut flags,
            ),
            Ok(0)
        );
        assert_eq!(flags, [0, 0]);
        assert_eq!(removed.nNumRemovedProtons, 77);

        let mut skipped = compare_one_fixture(&mut heap, 1);
        heap.slice_mut(skipped.generated[0][TAUT_YES as usize])
            .unwrap()[0]
            .bDeleted = 1;
        heap.slice_mut(skipped.input[0]).unwrap()[0].bDeleted = 1;
        flags = [0; 2];
        call(
            &mut heap,
            Some(&skipped.structure),
            skipped.input,
            TAUT_YES as i32,
            &mut removed,
            &mut flags,
        )
        .unwrap();
        assert_eq!(flags, [0, 0]);

        let mut fallback = compare_one_fixture(&mut heap, 1);
        fallback.structure.iMobileH = TAUT_NON as i8;
        heap.slice_mut(fallback.generated[0][TAUT_NON as usize])
            .unwrap()[0]
            .nNumberOfAtoms = 0;
        flags = [0; 2];
        call(
            &mut heap,
            Some(&fallback.structure),
            fallback.input,
            TAUT_NON as i32,
            &mut removed,
            &mut flags,
        )
        .unwrap();
        assert_ne!(
            flags[TAUT_NON as usize] & tagInchiCompareDiffBits_INCHIDIFF_COMP_HLAYER as INCHI_MODE,
            0
        );

        let mut fixed = compare_one_fixture(&mut heap, 1);
        fixed.structure.iMobileH = TAUT_NON as i8;
        fixed.structure.nNumRemovedProtonsMobHInChI = 9;
        heap.slice_mut(fixed.input[1]).unwrap()[0].nAtom =
            heap.allocate_model_storage(vec![6_u8, 8]).unwrap();
        flags = [0; 2];
        call(
            &mut heap,
            Some(&fixed.structure),
            fixed.input,
            TAUT_NON as i32,
            &mut removed,
            &mut flags,
        )
        .unwrap();
        assert_eq!(flags[TAUT_NON as usize], 0);
        assert_eq!(
            flags[TAUT_YES as usize],
            tagInchiCompareDiffBits_INCHIDIFF_PROBLEM as INCHI_MODE
                | tagInchiCompareDiffBits_INCHIDIFF_MOBH_PROTONS as INCHI_MODE
        );

        let mut multiple = compare_one_fixture(&mut heap, 2);
        flags = [0; 2];
        call(
            &mut heap,
            Some(&multiple.structure),
            multiple.input,
            TAUT_YES as i32,
            &mut removed,
            &mut flags,
        )
        .unwrap();
        assert_ne!(
            flags[TAUT_YES as usize] & tagInchiCompareDiffBits_INCHIDIFF_COMP_NUMBER as INCHI_MODE,
            0
        );
        assert_eq!(removed.nNumRemovedProtons, 1);
        heap.slice_mut(multiple.generated[1][TAUT_YES as usize])
            .unwrap()[0]
            .bDeleted = 1;
        flags = [0; 2];
        call(
            &mut heap,
            Some(&multiple.structure),
            multiple.input,
            TAUT_YES as i32,
            &mut removed,
            &mut flags,
        )
        .unwrap();
        assert_eq!(
            flags[TAUT_YES as usize] & tagInchiCompareDiffBits_INCHIDIFF_COMP_NUMBER as INCHI_MODE,
            0
        );
        assert_eq!(removed.nNumRemovedProtons, 3);
        assert_eq!(removed.nNumRemovedIsotopicH, [5, 7, 9]);

        let zero = compare_one_fixture(&mut heap, 0);
        flags = [0; 2];
        call(
            &mut heap,
            Some(&zero.structure),
            zero.input,
            TAUT_YES as i32,
            &mut removed,
            &mut flags,
        )
        .unwrap();
        assert_eq!(
            flags[TAUT_YES as usize],
            tagInchiCompareDiffBits_INCHIDIFF_COMP_NUMBER as INCHI_MODE
        );

        let mut reconnected = compare_one_fixture(&mut heap, 1);
        reconnected.structure.iINCHI = INCHI_REC as i8;
        flags = [0; 2];
        call(
            &mut heap,
            Some(&reconnected.structure),
            reconnected.input,
            TAUT_YES as i32,
            &mut removed,
            &mut flags,
        )
        .unwrap();
        assert_eq!(flags, [0, 0]);

        assert_eq!(
            call(
                &mut heap,
                Some(&fixture.structure),
                fixture.input,
                -1,
                &mut removed,
                &mut flags,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        fixture.structure.RevInChI.pINChI[INCHI_BAS as usize] = SourceMutPointer::null();
        assert_eq!(
            call(
                &mut heap,
                Some(&fixture.structure),
                fixture.input,
                TAUT_YES as i32,
                &mut removed,
                &mut flags,
            ),
            Ok(0)
        );
    }

    struct CompareAllFixture {
        structures: [[SourceMutPointer<StrFromINChI>; TAUT_NUM as usize]; INCHI_NUM as usize],
        input: InpInChI,
        one: CompareOneFixture,
    }

    fn compare_all_fixture(heap: &mut SourceHeap) -> CompareAllFixture {
        let one = compare_one_fixture(heap, 1);
        let structure_pointer = heap
            .allocate_model_storage(vec![one.structure.clone()])
            .unwrap();
        let mut structures = [[SourceMutPointer::null(); 2]; 2];
        structures[INCHI_BAS as usize][TAUT_YES as usize] = structure_pointer;
        let mut input = InpInChI::default();
        input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize] = one.input[0];
        input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 1;
        input.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize].nNumRemovedProtons = 1;
        input.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize].nNumRemovedIsotopicH = [2, 3, 4];
        CompareAllFixture {
            structures,
            input,
            one,
        }
    }

    #[test]
    fn source_port__ichirvr7__comparealloriginchitorevinchi__line_1548() {
        let mut heap = SourceHeap::default();
        let mut fixture = compare_all_fixture(&mut heap);
        fixture.input.CompareInchiFlags = [[INCHI_MODE::MAX; 2]; 2];
        assert_eq!(
            CompareAllOrigInchiToRevInChI(
                &mut heap,
                &fixture.structures,
                &mut fixture.input,
                0,
                17,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(fixture.input.CompareInchiFlags[0], [0, 0]);
        assert_eq!(fixture.input.CompareInchiFlags[1], [INCHI_MODE::MAX; 2]);

        fixture.input.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize].nNumRemovedProtons = 9;
        fixture.input.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize].nNumRemovedIsotopicH =
            [2, 8, 4];
        assert_eq!(
            CompareAllOrigInchiToRevInChI(
                &mut heap,
                &fixture.structures,
                &mut fixture.input,
                0,
                0,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(
            fixture.input.CompareInchiFlags[0][TAUT_YES as usize],
            tagInchiCompareDiffBits_INCHIDIFF_MOBH_PROTONS as INCHI_MODE
                | tagInchiCompareDiffBits_INCHIDIFF_MOB_ISO_H as INCHI_MODE
        );

        fixture.input.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize].pNumProtons = heap
            .allocate_model_storage(vec![COMPONENT_REM_PROTONS::default()])
            .unwrap();
        assert_eq!(
            CompareAllOrigInchiToRevInChI(
                &mut heap,
                &fixture.structures,
                &mut fixture.input,
                0,
                0,
                SourceMutPointer::null(),
            ),
            Ok(RI_ERR_PROGR)
        );

        let mut fixed = compare_one_fixture(&mut heap, 1);
        fixed.structure.iMobileH = TAUT_NON as i8;
        fixed.structure.nNumRemovedProtonsMobHInChI = 1;
        let fixed_structure = heap
            .allocate_model_storage(vec![fixed.structure.clone()])
            .unwrap();
        let mut fixed_structures = [[SourceMutPointer::null(); 2]; 2];
        fixed_structures[INCHI_BAS as usize][TAUT_NON as usize] = fixed_structure;
        let mut fixed_input = InpInChI::default();
        fixed_input.pInpInChI[INCHI_BAS as usize][TAUT_NON as usize] = fixed.input[0];
        fixed_input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize] = fixed.input[1];
        fixed_input.nNumComponents[INCHI_BAS as usize][TAUT_NON as usize] = 1;
        fixed_input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 1;
        assert_eq!(
            CompareAllOrigInchiToRevInChI(
                &mut heap,
                &fixed_structures,
                &mut fixed_input,
                1,
                0,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(fixed_input.CompareInchiFlags[0], [0, 0]);

        let fallback = compare_one_fixture(&mut heap, 1);
        heap.slice_mut(fallback.input[0]).unwrap()[0].nNumberOfAtoms = 0;
        let fallback_structure = heap
            .allocate_model_storage(vec![fallback.structure.clone()])
            .unwrap();
        let mut fallback_structures = [[SourceMutPointer::null(); 2]; 2];
        fallback_structures[INCHI_BAS as usize][TAUT_YES as usize] = fallback_structure;
        let mut fallback_input = InpInChI::default();
        fallback_input.pInpInChI[INCHI_BAS as usize][TAUT_NON as usize] = fallback.input[0];
        fallback_input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize] = fallback.input[1];
        fallback_input.nNumComponents[INCHI_BAS as usize][TAUT_NON as usize] = 1;
        fallback_input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 1;
        fallback_input.nNumProtons[INCHI_BAS as usize][TAUT_NON as usize].nNumRemovedProtons = 1;
        assert_eq!(
            CompareAllOrigInchiToRevInChI(
                &mut heap,
                &fallback_structures,
                &mut fallback_input,
                1,
                0,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );

        let reconnected = compare_one_fixture(&mut heap, 1);
        let mut reconnected_structure_value = reconnected.structure.clone();
        reconnected_structure_value.iINCHI = INCHI_REC as i8;
        reconnected_structure_value.RevInChI.pINChI[INCHI_REC as usize] =
            reconnected_structure_value.RevInChI.pINChI[INCHI_BAS as usize];
        reconnected_structure_value.RevInChI.pINChI_Aux[INCHI_REC as usize] =
            reconnected_structure_value.RevInChI.pINChI_Aux[INCHI_BAS as usize];
        reconnected_structure_value.RevInChI.num_components[INCHI_REC as usize] = 1;
        let reconnected_structure = heap
            .allocate_model_storage(vec![reconnected_structure_value])
            .unwrap();
        let mut reconnected_structures = [[SourceMutPointer::null(); 2]; 2];
        reconnected_structures[INCHI_REC as usize][TAUT_YES as usize] = reconnected_structure;
        let mut reconnected_input = InpInChI::default();
        reconnected_input.pInpInChI[INCHI_REC as usize][TAUT_YES as usize] = reconnected.input[0];
        reconnected_input.nNumComponents[INCHI_REC as usize][TAUT_YES as usize] = 1;
        reconnected_input.nNumProtons[INCHI_REC as usize][TAUT_YES as usize].nNumRemovedProtons = 1;
        reconnected_input.nNumProtons[INCHI_REC as usize][TAUT_YES as usize].nNumRemovedIsotopicH =
            [2, 3, 4];
        assert_eq!(
            CompareAllOrigInchiToRevInChI(
                &mut heap,
                &reconnected_structures,
                &mut reconnected_input,
                0,
                0,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );

        let first = compare_one_fixture(&mut heap, 1);
        let second = compare_one_fixture(&mut heap, 1);
        let structure_rows = heap
            .allocate_model_storage(vec![first.structure.clone(), second.structure.clone()])
            .unwrap();
        let first_input = heap.slice(first.input[0].as_const()).unwrap()[0].clone();
        let second_input = heap.slice(second.input[0].as_const()).unwrap()[0].clone();
        let input_rows = heap
            .allocate_model_storage(vec![first_input, second_input])
            .unwrap();
        let mut multiple_structures = [[SourceMutPointer::null(); 2]; 2];
        multiple_structures[INCHI_BAS as usize][TAUT_YES as usize] = structure_rows;
        let mut multiple_input = InpInChI::default();
        multiple_input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize] = input_rows;
        multiple_input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 2;
        multiple_input.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize].nNumRemovedProtons = 2;
        multiple_input.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize].nNumRemovedIsotopicH =
            [4, 6, 8];
        assert_eq!(
            CompareAllOrigInchiToRevInChI(
                &mut heap,
                &multiple_structures,
                &mut multiple_input,
                0,
                0,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(multiple_input.CompareInchiFlags[0], [0, 0]);

        let mut allocation = compare_all_fixture(&mut heap);
        let tautomer = heap
            .allocate_model_storage(vec![1_u16, 3, 0, 0, 1])
            .unwrap();
        heap.slice_mut(allocation.one.generated[0][TAUT_YES as usize])
            .unwrap()[0]
            .nTautomer = tautomer;
        heap.slice_mut(allocation.one.generated[0][TAUT_YES as usize])
            .unwrap()[0]
            .lenTautomer = 5;
        let input_tautomer = heap
            .allocate_model_storage(vec![1_u16, 3, 0, 0, 1])
            .unwrap();
        heap.slice_mut(allocation.one.input[0]).unwrap()[0].nTautomer = input_tautomer;
        heap.slice_mut(allocation.one.input[0]).unwrap()[0].lenTautomer = 5;
        heap.fail_after_allocations(0);
        assert_eq!(
            CompareAllOrigInchiToRevInChI(
                &mut heap,
                &allocation.structures,
                &mut allocation.input,
                0,
                0,
                SourceMutPointer::null(),
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(heap.source_allocation_calls(), 2);

        let mut malformed = compare_all_fixture(&mut heap);
        malformed.structures[INCHI_BAS as usize][TAUT_YES as usize] = SourceMutPointer::null();
        assert_eq!(
            CompareAllOrigInchiToRevInChI(
                &mut heap,
                &malformed.structures,
                &mut malformed.input,
                0,
                0,
                SourceMutPointer::null(),
            ),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichirvr7__comparereversedstereoinchi3__line_2362() {
        let mut heap = SourceHeap::default();
        let mut comparison = ICR::default();
        assert_eq!(
            CompareReversedStereoINChI3(&heap, None, None, &mut comparison),
            Ok(0)
        );
        let empty = INChI_Stereo::default();
        assert_eq!(
            CompareReversedStereoINChI3(&heap, Some(&empty), None, &mut comparison),
            Ok(0)
        );

        let equal1 = stereo3_fixture(&mut heap, &[(1, 1), (3, 2)], 1, &[(2, 4, 2)]);
        let equal2 = stereo3_fixture(&mut heap, &[(1, 1), (3, 2)], 1, &[(2, 4, 2)]);
        assert_eq!(
            CompareReversedStereoINChI3(&heap, Some(&equal1), Some(&equal2), &mut ICR::default(),),
            Ok(0)
        );

        let parity1 = stereo3_fixture(&mut heap, &[(1, 1)], 0, &[]);
        let parity2 = stereo3_fixture(&mut heap, &[(1, 2)], 0, &[]);
        assert_eq!(
            CompareReversedStereoINChI3(&heap, Some(&parity1), Some(&parity2), &mut ICR::default(),),
            Ok(tagInchiCompareDiffBits_INCHIDIFF_SC_PARITY as INCHI_MODE)
        );

        let extra_centers =
            stereo3_fixture(&mut heap, &[(1, AB_PARITY_UNDF as i8), (2, 1)], 0, &[]);
        let mut extra_result = ICR::default();
        assert_eq!(
            CompareReversedStereoINChI3(&heap, Some(&extra_centers), None, &mut extra_result,),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_SC_EXTRA_UNDF
                | tagInchiCompareDiffBits_INCHIDIFF_SC_EXTRA) as INCHI_MODE)
        );
        assert_eq!(extra_result.num_sc_in1_only, 2);
        assert_eq!(&extra_result.sc_in1_only[..2], &[0, 1]);
        assert_eq!(extra_result.num_sc_undef_in1_only, 1);
        assert_eq!(extra_result.sc_undef_in1_only[0], 0);

        let missing_centers =
            stereo3_fixture(&mut heap, &[(1, AB_PARITY_UNDF as i8), (2, 1)], 0, &[]);
        let mut missing_result = ICR::default();
        assert_eq!(
            CompareReversedStereoINChI3(&heap, None, Some(&missing_centers), &mut missing_result,),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_SC_MISS_UNDF
                | tagInchiCompareDiffBits_INCHIDIFF_SC_MISS) as INCHI_MODE)
        );
        assert_eq!(missing_result.num_sc_in2_only, 2);
        assert_eq!(&missing_result.sc_in2_only[..2], &[0, 1]);
        assert_eq!(missing_result.num_sc_undef_in2_only, 0);

        let merge1 = stereo3_fixture(&mut heap, &[(2, 1)], 0, &[]);
        let merge2 = stereo3_fixture(&mut heap, &[(1, AB_PARITY_UNDF as i8)], 0, &[]);
        let mut merge_result = ICR::default();
        assert_eq!(
            CompareReversedStereoINChI3(&heap, Some(&merge1), Some(&merge2), &mut merge_result,),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_SC_EXTRA
                | tagInchiCompareDiffBits_INCHIDIFF_SC_MISS_UNDF) as INCHI_MODE)
        );
        assert_eq!(merge_result.sc_undef_in2_only[0], 0);

        let mut no_add = ICR {
            num_sc_in1_only: 1,
            sc_in1_only: [17; ICR_MAX_SC_IN1_ONLY as usize],
            ..ICR::default()
        };
        CompareReversedStereoINChI3(&heap, Some(&extra_centers), None, &mut no_add).unwrap();
        assert_eq!(no_add.num_sc_in1_only, 1);
        assert_eq!(no_add.sc_in1_only, [17; ICR_MAX_SC_IN1_ONLY as usize]);

        let many_centers: Vec<_> = (0..34).map(|index| (index as AT_NUMB, 1)).collect();
        let many = stereo3_fixture(&mut heap, &many_centers, 0, &[]);
        let mut capped = ICR::default();
        CompareReversedStereoINChI3(&heap, Some(&many), None, &mut capped).unwrap();
        assert_eq!(capped.num_sc_in1_only, ICR_MAX_SC_IN1_ONLY as i32);
        assert_eq!(capped.sc_in1_only[31], 31);

        let inverted1 = stereo3_fixture(&mut heap, &[], 1, &[]);
        let inverted2 = stereo3_fixture(&mut heap, &[], -1, &[]);
        assert_eq!(
            CompareReversedStereoINChI3(
                &heap,
                Some(&inverted1),
                Some(&inverted2),
                &mut ICR::default(),
            ),
            Ok(tagInchiCompareDiffBits_INCHIDIFF_SC_INV as INCHI_MODE)
        );
        let suppressed = stereo3_fixture(&mut heap, &[], 2, &[]);
        assert_eq!(
            CompareReversedStereoINChI3(
                &heap,
                Some(&inverted1),
                Some(&suppressed),
                &mut ICR::default(),
            ),
            Ok(0)
        );

        let bond1 = stereo3_fixture(
            &mut heap,
            &[],
            0,
            &[(2, 2, AB_PARITY_UNDF as i8), (4, 4, 1)],
        );
        let bond2 = stereo3_fixture(
            &mut heap,
            &[],
            0,
            &[(1, 1, AB_PARITY_UNDF as i8), (4, 4, 2)],
        );
        let mut bond_result = ICR::default();
        assert_eq!(
            CompareReversedStereoINChI3(&heap, Some(&bond1), Some(&bond2), &mut bond_result,),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_SB_EXTRA_UNDF
                | tagInchiCompareDiffBits_INCHIDIFF_SB_MISS_UNDF
                | tagInchiCompareDiffBits_INCHIDIFF_SB_PARITY) as INCHI_MODE)
        );
        assert_eq!(bond_result.num_sb_in1_only, 1);
        assert_eq!(bond_result.num_sb_undef_in1_only, 1);
        assert_eq!(bond_result.num_sb_in2_only, 1);
        assert_eq!(bond_result.num_sb_undef_in2_only, 1);
        assert_eq!(bond_result.sb_undef_in2_only[0], 0);

        let extra_bonds = stereo3_fixture(
            &mut heap,
            &[],
            0,
            &[(1, 1, AB_PARITY_UNDF as i8), (2, 2, 1)],
        );
        assert_eq!(
            CompareReversedStereoINChI3(&heap, Some(&extra_bonds), None, &mut ICR::default(),),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_SB_EXTRA_UNDF
                | tagInchiCompareDiffBits_INCHIDIFF_SB_EXTRA) as INCHI_MODE)
        );
        assert_eq!(
            CompareReversedStereoINChI3(&heap, None, Some(&extra_bonds), &mut ICR::default(),),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_SB_MISS_UNDF
                | tagInchiCompareDiffBits_INCHIDIFF_SB_MISS) as INCHI_MODE)
        );

        let second_atom1 = stereo3_fixture(&mut heap, &[], 0, &[(1, 2, 1)]);
        let second_atom2 = stereo3_fixture(&mut heap, &[], 0, &[(1, 3, 1)]);
        assert_eq!(
            CompareReversedStereoINChI3(
                &heap,
                Some(&second_atom1),
                Some(&second_atom2),
                &mut ICR::default(),
            ),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_SB_EXTRA
                | tagInchiCompareDiffBits_INCHIDIFF_SB_MISS) as INCHI_MODE)
        );

        let negative = INChI_Stereo {
            nNumberOfStereoCenters: -1,
            ..INChI_Stereo::default()
        };
        assert_eq!(
            CompareReversedStereoINChI3(&heap, Some(&negative), None, &mut ICR::default(),),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        let malformed = INChI_Stereo {
            nNumberOfStereoBonds: 1,
            ..INChI_Stereo::default()
        };
        assert_eq!(
            CompareReversedStereoINChI3(&heap, Some(&malformed), None, &mut ICR::default(),),
            Err(SourceHeapError::NullPointer)
        );
    }

    fn reversed_inchi3_fixture(heap: &mut SourceHeap) -> INChI {
        INChI {
            nNumberOfAtoms: 2,
            szHillFormula: heap
                .allocate_model_storage(vec![b'C' as i8, b'2' as i8, b'H' as i8, b'6' as i8, 0])
                .unwrap(),
            nAtom: heap.allocate_model_storage(vec![6_u8, 6]).unwrap(),
            lenConnTable: 2,
            nConnTable: heap.allocate_model_storage(vec![1_u16, 2]).unwrap(),
            nNum_H: heap.allocate_model_storage(vec![3_i8, 3]).unwrap(),
            ..INChI::default()
        }
    }

    #[test]
    fn source_port__ichirvr7__comparereversedinchi3__line_2644() {
        fn compare(
            heap: &mut SourceHeap,
            left: Option<&INChI>,
            right: Option<&INChI>,
            aux_left: Option<&INChI_Aux>,
            aux_right: Option<&INChI_Aux>,
        ) -> Result<(u32, i32), SourceHeapError> {
            let mut error = 77;
            CompareReversedINChI3(heap, left, right, aux_left, aux_right, &mut error)
                .map(|difference| (difference as u32, error))
        }

        let mut heap = SourceHeap::default();
        assert_eq!(compare(&mut heap, None, None, None, None), Ok((0, 0)));
        let left = reversed_inchi3_fixture(&mut heap);
        let right = reversed_inchi3_fixture(&mut heap);
        assert_eq!(
            compare(&mut heap, Some(&left), None, None, None),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_PROBLEM, 0))
        );
        assert_eq!(
            compare(&mut heap, Some(&left), Some(&right), None, None),
            Ok((0, 0))
        );

        let mut error_left = left.clone();
        let mut error_right = right.clone();
        error_left.nErrorCode = 4;
        error_right.nErrorCode = 4;
        assert_eq!(
            compare(&mut heap, Some(&error_left), Some(&error_right), None, None),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_PROBLEM, 0))
        );
        error_right.nErrorCode = 5;
        assert_eq!(
            compare(&mut heap, Some(&error_left), Some(&error_right), None, None),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_PROBLEM, 0))
        );

        let mut changed = reversed_inchi3_fixture(&mut heap);
        changed.nNumberOfAtoms = 1;
        assert_eq!(
            compare(&mut heap, Some(&left), Some(&changed), None, None),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_NUM_AT, 0))
        );
        changed = reversed_inchi3_fixture(&mut heap);
        changed.nAtom = heap.allocate_model_storage(vec![6_u8, 7]).unwrap();
        assert_eq!(
            compare(&mut heap, Some(&left), Some(&changed), None, None),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_ATOMS, 0))
        );

        let mut hydrogen_left = reversed_inchi3_fixture(&mut heap);
        let mut hydrogen_right = reversed_inchi3_fixture(&mut heap);
        hydrogen_left.nNum_H = heap.allocate_model_storage(vec![i8::MIN, 4]).unwrap();
        hydrogen_right.nNum_H = heap.allocate_model_storage(vec![1_i8, 2]).unwrap();
        hydrogen_left.nNum_H_fixed = heap.allocate_model_storage(vec![2_i8, -1]).unwrap();
        hydrogen_right.nNum_H_fixed = heap.allocate_model_storage(vec![1_i8, 3]).unwrap();
        assert_eq!(
            compare(
                &mut heap,
                Some(&hydrogen_left),
                Some(&hydrogen_right),
                None,
                None,
            ),
            Ok((
                tagInchiCompareDiffBits_INCHIDIFF_POSITION_H
                    | tagInchiCompareDiffBits_INCHIDIFF_MORE_FH
                    | tagInchiCompareDiffBits_INCHIDIFF_LESS_FH,
                0,
            ))
        );
        let mut fixed_only = reversed_inchi3_fixture(&mut heap);
        fixed_only.nNum_H_fixed = heap.allocate_model_storage(vec![0_i8, -3]).unwrap();
        assert_eq!(
            compare(&mut heap, Some(&fixed_only), Some(&right), None, None),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_MORE_FH, 0))
        );
        assert_eq!(
            compare(&mut heap, Some(&right), Some(&fixed_only), None, None),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_LESS_FH, 0))
        );

        let mut fewer_h = reversed_inchi3_fixture(&mut heap);
        fewer_h.szHillFormula = heap
            .allocate_model_storage(vec![b'C' as i8, b'2' as i8, b'H' as i8, b'4' as i8, 0])
            .unwrap();
        assert_eq!(
            compare(&mut heap, Some(&left), Some(&fewer_h), None, None),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_MORE_H, 0))
        );
        assert_eq!(
            compare(&mut heap, Some(&fewer_h), Some(&left), None, None),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_LESS_H, 0))
        );
        let mut different_element = reversed_inchi3_fixture(&mut heap);
        different_element.szHillFormula = heap
            .allocate_model_storage(vec![b'N' as i8, b'2' as i8, b'H' as i8, b'6' as i8, 0])
            .unwrap();
        assert_eq!(
            compare(&mut heap, Some(&left), Some(&different_element), None, None,),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_NUM_EL, 0))
        );

        let mut connection = reversed_inchi3_fixture(&mut heap);
        connection.lenConnTable = 1;
        assert_eq!(
            compare(&mut heap, Some(&left), Some(&connection), None, None),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_CON_LEN, 0))
        );
        connection = reversed_inchi3_fixture(&mut heap);
        connection.nConnTable = heap.allocate_model_storage(vec![1_u16, 3]).unwrap();
        assert_eq!(
            compare(&mut heap, Some(&left), Some(&connection), None, None),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_CON_TBL, 0))
        );

        let mut tautomer_left = reversed_inchi3_fixture(&mut heap);
        tautomer_left.nTautomer = heap
            .allocate_model_storage(vec![1_u16, 4, 2, 1, 5, 1])
            .unwrap();
        tautomer_left.lenTautomer = 6;
        let mut tautomer_right = reversed_inchi3_fixture(&mut heap);
        tautomer_right.nTautomer = heap
            .allocate_model_storage(vec![1_u16, 4, 1, 0, 2, 5])
            .unwrap();
        tautomer_right.lenTautomer = 6;
        assert_eq!(
            compare(
                &mut heap,
                Some(&tautomer_left),
                Some(&tautomer_right),
                None,
                None,
            ),
            Ok((
                tagInchiCompareDiffBits_INCHIDIFF_DIFF_TG_ENDP
                    | tagInchiCompareDiffBits_INCHIDIFF_TG,
                0,
            ))
        );

        let one_group = |heap: &mut SourceHeap| {
            let mut value = reversed_inchi3_fixture(heap);
            value.nTautomer = heap
                .allocate_model_storage(vec![1_u16, 3, 0, 0, 1])
                .unwrap();
            value.lenTautomer = 5;
            value
        };
        let two_groups = |heap: &mut SourceHeap| {
            let mut value = reversed_inchi3_fixture(heap);
            value.nTautomer = heap
                .allocate_model_storage(vec![2_u16, 3, 0, 0, 1, 3, 0, 0, 2])
                .unwrap();
            value.lenTautomer = 9;
            value
        };
        let three_groups = |heap: &mut SourceHeap| {
            let mut value = reversed_inchi3_fixture(heap);
            value.nTautomer = heap
                .allocate_model_storage(vec![3_u16, 3, 0, 0, 1, 3, 0, 0, 2, 3, 0, 0, 3])
                .unwrap();
            value.lenTautomer = 13;
            value
        };
        let one = one_group(&mut heap);
        let two = two_groups(&mut heap);
        let three = three_groups(&mut heap);
        assert_eq!(
            compare(&mut heap, Some(&one), Some(&two), None, None),
            Ok((
                tagInchiCompareDiffBits_INCHIDIFF_SINGLE_TG
                    | tagInchiCompareDiffBits_INCHIDIFF_DIFF_TG_ENDP
                    | tagInchiCompareDiffBits_INCHIDIFF_TG,
                0,
            ))
        );
        assert_eq!(
            compare(&mut heap, Some(&two), Some(&one), None, None),
            Ok((
                tagInchiCompareDiffBits_INCHIDIFF_MULTIPLE_TG
                    | tagInchiCompareDiffBits_INCHIDIFF_DIFF_TG_ENDP
                    | tagInchiCompareDiffBits_INCHIDIFF_TG,
                0,
            ))
        );
        assert_eq!(
            compare(&mut heap, Some(&two), Some(&three), None, None),
            Ok((
                tagInchiCompareDiffBits_INCHIDIFF_NUM_TG
                    | tagInchiCompareDiffBits_INCHIDIFF_DIFF_TG_ENDP,
                0,
            ))
        );
        assert_eq!(
            compare(&mut heap, Some(&right), Some(&one), None, None),
            Ok((
                tagInchiCompareDiffBits_INCHIDIFF_NO_TAUT
                    | tagInchiCompareDiffBits_INCHIDIFF_DIFF_TG_ENDP,
                0,
            ))
        );
        assert_eq!(
            compare(&mut heap, Some(&one), Some(&right), None, None),
            Ok((
                tagInchiCompareDiffBits_INCHIDIFF_WRONG_TAUT
                    | tagInchiCompareDiffBits_INCHIDIFF_DIFF_TG_ENDP,
                0,
            ))
        );

        let baseline_allocations = heap.live_allocation_count();
        heap.fail_after_allocations(0);
        assert_eq!(
            compare(
                &mut heap,
                Some(&tautomer_left),
                Some(&tautomer_right),
                None,
                None,
            ),
            Ok((0, RI_ERR_ALLOC))
        );
        assert_eq!(heap.source_allocation_calls(), 2);
        assert_eq!(heap.live_allocation_count(), baseline_allocations);
        heap.fail_after_allocations(1);
        assert_eq!(
            compare(
                &mut heap,
                Some(&tautomer_left),
                Some(&tautomer_right),
                None,
                None,
            ),
            Ok((0, RI_ERR_ALLOC))
        );
        assert_eq!(heap.source_allocation_calls(), 2);
        assert_eq!(heap.live_allocation_count(), baseline_allocations);

        let mut isotope_left = reversed_inchi3_fixture(&mut heap);
        isotope_left.nNumberOfIsotopicAtoms = 1;
        isotope_left.IsotopicAtom = heap
            .allocate_model_storage(vec![INChI_IsotopicAtom {
                nAtomNumber: 1,
                ..INChI_IsotopicAtom::default()
            }])
            .unwrap();
        let mut isotope_right = reversed_inchi3_fixture(&mut heap);
        assert_eq!(
            compare(
                &mut heap,
                Some(&isotope_left),
                Some(&isotope_right),
                None,
                None,
            ),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_NUM_ISO_AT, 0))
        );
        isotope_right.nNumberOfIsotopicAtoms = 1;
        isotope_right.IsotopicAtom = heap
            .allocate_model_storage(vec![INChI_IsotopicAtom {
                nAtomNumber: 2,
                ..INChI_IsotopicAtom::default()
            }])
            .unwrap();
        assert_eq!(
            compare(
                &mut heap,
                Some(&isotope_left),
                Some(&isotope_right),
                None,
                None,
            ),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_ISO_AT, 0))
        );

        let mut charged = reversed_inchi3_fixture(&mut heap);
        charged.nTotalCharge = -1;
        assert_eq!(
            compare(&mut heap, Some(&charged), Some(&right), None, None),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_CHARGE, 0))
        );
        let aux = INChI_Aux {
            nNumRemovedProtons: 2,
            nNumRemovedIsotopicH: [1, 2, 3],
            ..INChI_Aux::default()
        };
        assert_eq!(
            compare(&mut heap, Some(&left), Some(&right), Some(&aux), None),
            Ok((
                tagInchiCompareDiffBits_INCHIDIFF_REM_PROT
                    | tagInchiCompareDiffBits_INCHIDIFF_REM_ISO_H,
                0,
            ))
        );
        assert_eq!(
            compare(
                &mut heap,
                Some(&left),
                Some(&right),
                Some(&INChI_Aux::default()),
                None,
            ),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_REM_ISO_H, 0))
        );

        let regular_left = stereo3_fixture(&mut heap, &[(1, 1)], 0, &[]);
        let regular_right = stereo3_fixture(&mut heap, &[(1, 2)], 0, &[]);
        let mut stereo_left = reversed_inchi3_fixture(&mut heap);
        stereo_left.Stereo = heap.allocate_model_storage(vec![regular_left]).unwrap();
        stereo_left.StereoIsotopic = heap
            .allocate_model_storage(vec![INChI_Stereo::default()])
            .unwrap();
        let mut stereo_right = reversed_inchi3_fixture(&mut heap);
        stereo_right.Stereo = heap.allocate_model_storage(vec![regular_right]).unwrap();
        stereo_right.StereoIsotopic = heap
            .allocate_model_storage(vec![INChI_Stereo::default()])
            .unwrap();
        assert_eq!(
            compare(
                &mut heap,
                Some(&stereo_left),
                Some(&stereo_right),
                None,
                None,
            ),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_SC_PARITY, 0))
        );
        let isotopic_left = stereo3_fixture(&mut heap, &[(2, 1)], 0, &[]);
        let isotopic_right = stereo3_fixture(&mut heap, &[(2, 2)], 0, &[]);
        stereo_left.StereoIsotopic = heap.allocate_model_storage(vec![isotopic_left]).unwrap();
        stereo_right.StereoIsotopic = heap.allocate_model_storage(vec![isotopic_right]).unwrap();
        assert_eq!(
            compare(
                &mut heap,
                Some(&stereo_left),
                Some(&stereo_right),
                None,
                None,
            ),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_SC_PARITY, 0))
        );
    }

    #[test]
    fn compare_reversed_inchi_uses_the_source_allocated_extent_for_pe2_j1() {
        let mut heap = SourceHeap::default();
        let mut left = reversed_inchi3_fixture(&mut heap);
        left.nTautomer = heap
            .allocate_model_storage(vec![1_u16, 4, 0, 0, 1, 2])
            .unwrap();
        left.lenTautomer = 6;
        let mut right = reversed_inchi3_fixture(&mut heap);
        right.nTautomer = heap
            .allocate_model_storage(vec![1_u16, 3, 0, 0, 3])
            .unwrap();
        right.lenTautomer = 5;
        let mut error = 77;

        assert_eq!(
            CompareReversedINChI3(&mut heap, Some(&left), Some(&right), None, None, &mut error),
            Ok((tagInchiCompareDiffBits_INCHIDIFF_EXTRA_TG_ENDP
                | tagInchiCompareDiffBits_INCHIDIFF_MISS_TG_ENDP
                | tagInchiCompareDiffBits_INCHIDIFF_DIFF_TG_ENDP
                | tagInchiCompareDiffBits_INCHIDIFF_TG) as INCHI_MODE)
        );
        assert_eq!(error, 0);
    }

    #[test]
    fn source_port__ichirvr7__inchi2atom__line_101() {
        struct Fixture {
            heap: SourceHeap,
            clock: SourceMutPointer<INCHI_CLOCK>,
            parameters: INPUT_PARMS,
            structure: StrFromINChI,
            input: InpInChI,
        }

        fn fixture() -> Fixture {
            let mut heap = SourceHeap::default();
            let clock = heap
                .allocate_model_storage(vec![INCHI_CLOCK::default()])
                .unwrap();
            let restore_mode = heap.allocate_model_storage(vec![SRM::default()]).unwrap();
            let mut carbon = inp_ATOM {
                el_number: 6,
                num_H: 4,
                orig_at_number: 1,
                component: 1,
                ..inp_ATOM::default()
            };
            carbon.elname[0] = b'C' as i8;
            let atoms = heap.allocate_model_storage(vec![carbon]).unwrap();
            let scratch = heap
                .allocate_model_storage(vec![inp_ATOM::default()])
                .unwrap();
            let parameters = INPUT_PARMS {
                nMode: u64::from(REQ_MODE_TAUT | REQ_MODE_NON_ISO),
                bTautFlags: u64::from(TG_FLAG_FIX_ISO_FIXEDH_BUG | TG_FLAG_FIX_TERM_H_CHRG_BUG),
                ..INPUT_PARMS::default()
            };
            let mut generated = StrFromINChI {
                num_atoms: 1,
                bMobileH: TAUT_YES as i8,
                iMobileH: TAUT_YES as i8,
                ..StrFromINChI::default()
            };
            assert_eq!(
                MakeOneInChIOutOfStrFromINChI(
                    &mut heap,
                    &mut CANON_GLOBALS::default(),
                    clock,
                    &parameters,
                    &mut STRUCT_DATA::default(),
                    &mut generated,
                    atoms,
                    scratch,
                    &ALL_TC_GROUPS::default(),
                    0,
                ),
                Ok(0)
            );
            let component = generated.pOneINChI[0];
            assert!(!component.is_null());
            generated.pOneINChI[0] = SourceMutPointer::null();
            let mut input = InpInChI {
                num_inp: 91,
                ..InpInChI::default()
            };
            input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize] = component;
            input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 1;
            Fixture {
                heap,
                clock,
                parameters,
                structure: StrFromINChI {
                    pSrm: restore_mode.as_const(),
                    ..StrFromINChI::default()
                },
                input,
            }
        }

        fn invoke(
            fixture: &mut Fixture,
            flags: i32,
            atom_offset: i32,
        ) -> Result<i32, SourceHeapError> {
            InChI2Atom(
                &mut fixture.heap,
                fixture.clock,
                &mut CANON_GLOBALS::default(),
                &fixture.parameters,
                &mut STRUCT_DATA::default(),
                SourceConstPointer::null(),
                17,
                &mut fixture.structure,
                0,
                atom_offset,
                flags,
                0,
                &mut fixture.input,
                0,
            )
        }

        let mut missing = InpInChI::default();
        let mut missing_structure = StrFromINChI {
            bFixedHExists: 7,
            ..StrFromINChI::default()
        };
        assert_eq!(
            InChI2Atom(
                &mut SourceHeap::default(),
                SourceMutPointer::null(),
                &mut CANON_GLOBALS::default(),
                &INPUT_PARMS::default(),
                &mut STRUCT_DATA::default(),
                SourceConstPointer::null(),
                0,
                &mut missing_structure,
                0,
                0,
                0,
                0,
                &mut missing,
                0,
            ),
            Ok(0)
        );
        assert_eq!(missing_structure.bFixedHExists, 7);

        let mut deleted_heap = SourceHeap::default();
        let deleted = deleted_heap
            .allocate_model_storage(vec![INChI {
                bDeleted: 1,
                ..INChI::default()
            }])
            .unwrap();
        let mut deleted_input = InpInChI::default();
        deleted_input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize] = deleted;
        deleted_input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 1;
        let mut deleted_structure = StrFromINChI {
            bFixedHExists: 9,
            ..StrFromINChI::default()
        };
        assert_eq!(
            InChI2Atom(
                &mut deleted_heap,
                SourceMutPointer::null(),
                &mut CANON_GLOBALS::default(),
                &INPUT_PARMS::default(),
                &mut STRUCT_DATA::default(),
                SourceConstPointer::null(),
                0,
                &mut deleted_structure,
                0,
                0,
                (I2A_FLAG_RECMET | I2A_FLAG_FIXEDH) as i32,
                0,
                &mut deleted_input,
                0,
            ),
            Ok(0)
        );
        assert_eq!(deleted_structure.bFixedHExists, 0);
        assert_eq!(deleted_structure.iINCHI, INCHI_BAS as i8);
        assert_eq!(deleted_structure.bMobileH, TAUT_YES as i8);

        let mut fixed_missing = InpInChI::default();
        fixed_missing.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 2;
        fixed_missing.nNumComponents[INCHI_BAS as usize][TAUT_NON as usize] = 1;
        let mut fixed_missing_structure = StrFromINChI {
            bFixedHExists: 11,
            ..StrFromINChI::default()
        };
        assert_eq!(
            InChI2Atom(
                &mut SourceHeap::default(),
                SourceMutPointer::null(),
                &mut CANON_GLOBALS::default(),
                &INPUT_PARMS::default(),
                &mut STRUCT_DATA::default(),
                SourceConstPointer::null(),
                0,
                &mut fixed_missing_structure,
                1,
                0,
                I2A_FLAG_FIXEDH as i32,
                0,
                &mut fixed_missing,
                0,
            ),
            Ok(0)
        );
        assert_eq!(fixed_missing_structure.bFixedHExists, 0);

        let mut negative = fixture();
        assert_eq!(
            InChI2Atom(
                &mut negative.heap,
                negative.clock,
                &mut CANON_GLOBALS::default(),
                &negative.parameters,
                &mut STRUCT_DATA::default(),
                SourceConstPointer::null(),
                0,
                &mut negative.structure,
                -1,
                0,
                0,
                0,
                &mut negative.input,
                0,
            ),
            Err(SourceHeapError::PointerOffsetOverflow)
        );

        let mut fixed = fixture();
        let component_value = fixed
            .heap
            .slice(fixed.input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize].as_const())
            .unwrap()[0]
            .clone();
        let fixed_component = fixed
            .heap
            .allocate_model_storage(vec![component_value])
            .unwrap();
        fixed.input.pInpInChI[INCHI_BAS as usize][TAUT_NON as usize] = fixed_component;
        fixed.input.nNumComponents[INCHI_BAS as usize][TAUT_NON as usize] = 1;
        let protons = fixed
            .heap
            .allocate_model_storage(vec![COMPONENT_REM_PROTONS {
                nNumRemovedProtons: 13,
                ..COMPONENT_REM_PROTONS::default()
            }])
            .unwrap();
        fixed.input.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize].pNumProtons = protons;
        fixed.heap.fail_after_allocations(0);
        assert_eq!(
            invoke(&mut fixed, I2A_FLAG_FIXEDH as i32, 23),
            Ok(RI_ERR_ALLOC)
        );
        assert_eq!(fixed.structure.bMobileH, TAUT_NON as i8);
        assert_eq!(fixed.structure.bFixedHExists, 1);
        assert_eq!(fixed.structure.nNumRemovedProtonsMobHInChI, 13);
        assert_eq!(fixed.structure.num_inp_actual, 91);

        let mut mobile_with_fixed = fixture();
        let fixed_value = mobile_with_fixed
            .heap
            .slice(
                mobile_with_fixed.input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize].as_const(),
            )
            .unwrap()[0]
            .clone();
        let fixed_pointer = mobile_with_fixed
            .heap
            .allocate_model_storage(vec![fixed_value])
            .unwrap();
        mobile_with_fixed.input.pInpInChI[INCHI_BAS as usize][TAUT_NON as usize] = fixed_pointer;
        mobile_with_fixed.input.nNumComponents[INCHI_BAS as usize][TAUT_NON as usize] = 1;
        mobile_with_fixed.heap.fail_after_allocations(0);
        assert_eq!(invoke(&mut mobile_with_fixed, 0, 0), Ok(RI_ERR_ALLOC));
        assert_eq!(mobile_with_fixed.structure.bFixedHExists, 1);

        let mut failed_polymer = fixture();
        let source_atom = failed_polymer
            .heap
            .slice(failed_polymer.input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize].as_const())
            .unwrap()[0]
            .nAtom;
        failed_polymer.heap.slice_mut(source_atom).unwrap()[0] = EL_NUMBER_ZZ;
        failed_polymer.input.polymer = failed_polymer
            .heap
            .allocate_model_storage(vec![OAD_Polymer::default()])
            .unwrap();
        failed_polymer.heap.fail_after_allocations(0);
        assert_eq!(invoke(&mut failed_polymer, 0, 40), Ok(RI_ERR_ALLOC));
        assert_eq!(
            failed_polymer.heap.slice(source_atom.as_const()).unwrap()[0],
            EL_NUMBER_ZY
        );
        assert_eq!(failed_polymer.structure.n_zy, 1);
        assert_eq!(failed_polymer.structure.n_pzz, 0);

        let mut capped = fixture();
        let capped_source_atom = capped
            .heap
            .slice(capped.input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize].as_const())
            .unwrap()[0]
            .nAtom;
        capped.heap.slice_mut(capped_source_atom).unwrap()[0] = EL_NUMBER_ZZ;
        let unit = capped
            .heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                cap1: -1,
                cap2: 41,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let units = capped.heap.allocate_model_storage(vec![unit]).unwrap();
        capped.input.polymer = capped
            .heap
            .allocate_model_storage(vec![OAD_Polymer {
                units,
                n: 1,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        capped.heap.fail_after_allocations(0);
        assert_eq!(invoke(&mut capped, 0, 40), Ok(RI_ERR_ALLOC));
        assert_eq!(
            capped.heap.slice(capped_source_atom.as_const()).unwrap()[0],
            EL_NUMBER_ZZ
        );
        assert_eq!(capped.structure.n_zy, 0);
        assert_eq!(capped.structure.n_pzz, 1);

        let mut restored = fixture();
        let restored_source_atom = restored
            .heap
            .slice(restored.input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize].as_const())
            .unwrap()[0]
            .nAtom;
        restored.heap.slice_mut(restored_source_atom).unwrap()[0] = EL_NUMBER_ZZ;
        restored.input.polymer = restored
            .heap
            .allocate_model_storage(vec![OAD_Polymer::default()])
            .unwrap();
        assert_eq!(invoke(&mut restored, 0, 40), Ok(1));
        assert_eq!(
            restored
                .heap
                .slice(restored_source_atom.as_const())
                .unwrap()[0],
            EL_NUMBER_ZZ
        );
        assert_eq!(restored.structure.n_zy, 0);
        assert_eq!(restored.structure.n_pzz, 1);
        for output in [restored.structure.at, restored.structure.at2] {
            let atom = &restored.heap.slice(output.as_const()).unwrap()[0];
            assert_eq!(atom.el_number, EL_NUMBER_ZZ);
            assert_eq!(&atom.elname[..3], &[b'Z' as i8, b'z' as i8, 0]);
        }
    }

    #[test]
    fn source_port__ichirvr7__allinchitostructure__line_1042() {
        struct GeneratedFixture {
            heap: SourceHeap,
            clock: SourceMutPointer<INCHI_CLOCK>,
            parameters: INPUT_PARMS,
            restore_mode: SourceConstPointer<SRM>,
            component: SourceMutPointer<INChI>,
        }

        fn generated_fixture() -> GeneratedFixture {
            let mut heap = SourceHeap::default();
            let clock = heap
                .allocate_model_storage(vec![INCHI_CLOCK::default()])
                .unwrap();
            let restore_mode = heap.allocate_model_storage(vec![SRM::default()]).unwrap();
            let mut carbon = inp_ATOM {
                el_number: 6,
                num_H: 4,
                orig_at_number: 1,
                component: 1,
                ..inp_ATOM::default()
            };
            carbon.elname[0] = b'C' as i8;
            let atoms = heap.allocate_model_storage(vec![carbon]).unwrap();
            let scratch = heap
                .allocate_model_storage(vec![inp_ATOM::default()])
                .unwrap();
            let parameters = INPUT_PARMS {
                nMode: u64::from(REQ_MODE_TAUT | REQ_MODE_NON_ISO | REQ_MODE_BASIC),
                bTautFlags: u64::from(TG_FLAG_FIX_ISO_FIXEDH_BUG | TG_FLAG_FIX_TERM_H_CHRG_BUG),
                ..INPUT_PARMS::default()
            };
            let mut generated = StrFromINChI {
                num_atoms: 1,
                bMobileH: TAUT_YES as i8,
                iMobileH: TAUT_YES as i8,
                ..StrFromINChI::default()
            };
            assert_eq!(
                MakeOneInChIOutOfStrFromINChI(
                    &mut heap,
                    &mut CANON_GLOBALS::default(),
                    clock,
                    &parameters,
                    &mut STRUCT_DATA::default(),
                    &mut generated,
                    atoms,
                    scratch,
                    &ALL_TC_GROUPS::default(),
                    0,
                ),
                Ok(0)
            );
            let component = generated.pOneINChI[0];
            assert!(!component.is_null());
            generated.pOneINChI[0] = SourceMutPointer::null();
            GeneratedFixture {
                heap,
                clock,
                parameters,
                restore_mode: restore_mode.as_const(),
                component,
            }
        }

        let mut empty_heap = SourceHeap::default();
        let empty_clock = empty_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let sentinel = empty_heap
            .allocate_model_storage(vec![StrFromINChI::default()])
            .unwrap();
        let mut empty_outputs = [[SourceMutPointer::null(); 2]; 2];
        empty_outputs[1][1] = sentinel;
        let empty_parameters = INPUT_PARMS {
            nMode: 0x1234_5678,
            ..INPUT_PARMS::default()
        };
        let empty_data = STRUCT_DATA {
            ulStructTime: 73,
            nErrorCode: 19,
            ..STRUCT_DATA::default()
        };
        assert_eq!(
            AllInchiToStructure(
                &mut empty_heap,
                empty_clock,
                &mut CANON_GLOBALS::default(),
                &empty_parameters,
                &empty_data,
                5,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                0,
                &mut empty_outputs,
                &mut InpInChI::default(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(empty_outputs[1][1], sentinel);
        assert_eq!(empty_parameters.nMode, 0x1234_5678);
        assert_eq!(empty_data.ulStructTime, 73);
        assert_eq!(empty_data.nErrorCode, 19);

        let mut allocation_heap = SourceHeap::default();
        let allocation_clock = allocation_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let mut allocation_input = InpInChI::default();
        allocation_input.nNumComponents[INCHI_BAS as usize][TAUT_NON as usize] = 1;
        let mut allocation_outputs = [[SourceMutPointer::null(); 2]; 2];
        allocation_heap.fail_after_allocations(0);
        assert_eq!(
            AllInchiToStructure(
                &mut allocation_heap,
                allocation_clock,
                &mut CANON_GLOBALS::default(),
                &INPUT_PARMS::default(),
                &STRUCT_DATA::default(),
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                0,
                &mut allocation_outputs,
                &mut allocation_input,
                0,
            ),
            Ok(RI_ERR_ALLOC)
        );
        assert!(allocation_outputs[0][0].is_null());

        let mut negative_heap = SourceHeap::default();
        let negative_clock = negative_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let mut negative_input = InpInChI::default();
        negative_input.nNumComponents[INCHI_BAS as usize][TAUT_NON as usize] = -1;
        let mut negative_outputs = [[SourceMutPointer::null(); 2]; 2];
        assert_eq!(
            AllInchiToStructure(
                &mut negative_heap,
                negative_clock,
                &mut CANON_GLOBALS::default(),
                &INPUT_PARMS::default(),
                &STRUCT_DATA::default(),
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                0,
                &mut negative_outputs,
                &mut negative_input,
                0,
            ),
            Ok(RI_ERR_ALLOC)
        );

        let mut skipped_heap = SourceHeap::default();
        let skipped_clock = skipped_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let skipped_components = skipped_heap
            .allocate_model_storage(vec![
                INChI {
                    nNumberOfAtoms: 0,
                    nLink: 7,
                    ..INChI::default()
                },
                INChI {
                    nNumberOfAtoms: 1,
                    bDeleted: 258,
                    nLink: 8,
                    ..INChI::default()
                },
                INChI {
                    nNumberOfAtoms: 1,
                    nLink: -3,
                    ..INChI::default()
                },
            ])
            .unwrap();
        let mut skipped_input = InpInChI::default();
        skipped_input.pInpInChI[INCHI_BAS as usize][TAUT_NON as usize] = skipped_components;
        skipped_input.nNumComponents[INCHI_BAS as usize][TAUT_NON as usize] = 3;
        let mut skipped_outputs = [[SourceMutPointer::null(); 2]; 2];
        assert_eq!(
            AllInchiToStructure(
                &mut skipped_heap,
                skipped_clock,
                &mut CANON_GLOBALS::default(),
                &INPUT_PARMS::default(),
                &STRUCT_DATA::default(),
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                0,
                &mut skipped_outputs,
                &mut skipped_input,
                0,
            ),
            Ok(0)
        );
        let skipped = skipped_heap
            .slice(skipped_outputs[0][TAUT_NON as usize].as_const())
            .unwrap();
        assert_eq!((skipped[0].nLink, skipped[0].bDeleted), (7, 0));
        assert_eq!((skipped[1].nLink, skipped[1].bDeleted), (8, 2));
        assert_eq!((skipped[2].nLink, skipped[2].bDeleted), (-3, 0));

        let mut fixed_skip = generated_fixture();
        let fixed_value = fixed_skip
            .heap
            .slice(fixed_skip.component.as_const())
            .unwrap()[0]
            .clone();
        let fixed_component = fixed_skip
            .heap
            .allocate_model_storage(vec![INChI {
                bDeleted: 1,
                nLink: 12,
                ..fixed_value
            }])
            .unwrap();
        let mut fixed_skip_input = InpInChI::default();
        fixed_skip_input.pInpInChI[INCHI_BAS as usize][TAUT_NON as usize] = fixed_component;
        fixed_skip_input.nNumComponents[INCHI_BAS as usize][TAUT_NON as usize] = 1;
        fixed_skip_input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize] = fixed_skip.component;
        fixed_skip_input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 1;
        let mut fixed_skip_outputs = [[SourceMutPointer::null(); 2]; 2];
        assert_eq!(
            AllInchiToStructure(
                &mut fixed_skip.heap,
                fixed_skip.clock,
                &mut CANON_GLOBALS::default(),
                &fixed_skip.parameters,
                &STRUCT_DATA::default(),
                0,
                SourceConstPointer::null(),
                fixed_skip.restore_mode,
                1,
                &mut fixed_skip_outputs,
                &mut fixed_skip_input,
                0,
            ),
            Ok(0)
        );
        let fixed_target = &fixed_skip
            .heap
            .slice(fixed_skip_outputs[0][TAUT_NON as usize].as_const())
            .unwrap()[0];
        assert_eq!((fixed_target.nLink, fixed_target.bDeleted), (12, 1));
        assert_eq!(
            fixed_skip
                .heap
                .slice(fixed_skip_outputs[0][TAUT_YES as usize].as_const())
                .unwrap()[0],
            StrFromINChI::default()
        );

        let mut counted = generated_fixture();
        let mut counted_input = InpInChI::default();
        counted_input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize] = counted.component;
        counted_input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 1;
        let mut counted_outputs = [[SourceMutPointer::null(); 2]; 2];
        counted.heap.fail_after_allocations(1);
        assert_eq!(
            AllInchiToStructure(
                &mut counted.heap,
                counted.clock,
                &mut CANON_GLOBALS::default(),
                &counted.parameters,
                &STRUCT_DATA::default(),
                0,
                SourceConstPointer::null(),
                counted.restore_mode,
                0,
                &mut counted_outputs,
                &mut counted_input,
                0,
            ),
            Ok(1)
        );
        let counted_target = &counted
            .heap
            .slice(counted_outputs[0][TAUT_YES as usize].as_const())
            .unwrap()[0];
        assert_eq!(counted_target.nError, RI_ERR_ALLOC);
        assert_eq!(counted_target.iInchiRec, INCHI_BAS as i8);
        assert_eq!(counted_target.iMobileH, TAUT_YES as i8);

        let mut complete = generated_fixture();
        let original_mode = complete.parameters.nMode;
        let mut complete_input = InpInChI {
            num_inp: 101,
            ..InpInChI::default()
        };
        complete_input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize] = complete.component;
        complete_input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 1;
        let mut complete_outputs = [[SourceMutPointer::null(); 2]; 2];
        let complete_data = STRUCT_DATA {
            ulStructTime: 211,
            nErrorCode: 37,
            ..STRUCT_DATA::default()
        };
        assert_eq!(
            AllInchiToStructure(
                &mut complete.heap,
                complete.clock,
                &mut CANON_GLOBALS::default(),
                &complete.parameters,
                &complete_data,
                41,
                SourceConstPointer::null(),
                complete.restore_mode,
                0,
                &mut complete_outputs,
                &mut complete_input,
                0,
            ),
            Ok(0)
        );
        assert_eq!(complete.parameters.nMode, original_mode);
        assert_eq!(complete_data.ulStructTime, 211);
        assert_eq!(complete_data.nErrorCode, 37);
        let complete_target = &complete
            .heap
            .slice(complete_outputs[0][TAUT_YES as usize].as_const())
            .unwrap()[0];
        assert_eq!(complete_target.pSrm, complete.restore_mode);
        assert_eq!(complete_target.iInchiRec, INCHI_BAS as i8);
        assert_eq!(complete_target.iMobileH, TAUT_YES as i8);
        assert_eq!(complete_target.num_inp_actual, 101);
        assert!(!complete_target.at.is_null());
        assert!(!complete_target.at2.is_null());
        assert_eq!(complete_target.nError, 0);
    }

    #[test]
    fn source_port__ichirvr7__addprotonandisohbalancetomobhstruct__line_1180() {
        fn call(
            heap: &mut SourceHeap,
            fixed_h: i32,
            structures: &[[SourceMutPointer<StrFromINChI>; 2]; 2],
            input: &InpInChI,
        ) -> Result<i32, SourceHeapError> {
            AddProtonAndIsoHBalanceToMobHStruct(
                heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                &INPUT_PARMS::default(),
                &mut STRUCT_DATA::default(),
                23,
                fixed_h,
                SourceMutPointer::null(),
                structures,
                input,
                0,
            )
        }

        fn add_layer(
            heap: &mut SourceHeap,
            structure: &mut StrFromINChI,
            representation: usize,
            charge: Option<i32>,
            removed: i16,
        ) {
            if let Some(charge) = charge {
                let inchi = heap
                    .allocate_model_storage(vec![INChI {
                        nTotalCharge: charge,
                        ..INChI::default()
                    }])
                    .unwrap();
                structure.RevInChI.pINChI[representation] = heap
                    .allocate_model_storage(vec![[SourceMutPointer::null(), inchi]])
                    .unwrap();
            }
            let aux = heap
                .allocate_model_storage(vec![INChI_Aux {
                    nNumRemovedProtons: removed,
                    ..INChI_Aux::default()
                }])
                .unwrap();
            structure.RevInChI.pINChI_Aux[representation] = heap
                .allocate_model_storage(vec![[SourceMutPointer::null(), aux]])
                .unwrap();
            structure.RevInChI.num_components[representation] = 1;
        }

        let null_structures = [[SourceMutPointer::null(); 2]; 2];
        let mut fixed_input = InpInChI::default();
        fixed_input.nNumComponents = [[i32::MAX; 2]; 2];
        fixed_input.nNumProtons[0][1].nNumRemovedProtons = i16::MIN;
        let fixed_before = fixed_input.clone();
        assert_eq!(
            call(
                &mut SourceHeap::default(),
                1,
                &null_structures,
                &fixed_input,
            ),
            Ok(0)
        );
        assert_eq!(fixed_input, fixed_before);

        let mut no_components = InpInChI::default();
        no_components.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize].nNumRemovedProtons = 7;
        no_components.nNumProtons[INCHI_REC as usize][TAUT_YES as usize].nNumRemovedIsotopicH =
            [11, 13, 17];
        assert_eq!(
            call(
                &mut SourceHeap::default(),
                0,
                &null_structures,
                &no_components,
            ),
            Ok(0)
        );

        let mut direct_heap = SourceHeap::default();
        let original = direct_heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                nTotalCharge: -4,
                ..INChI::default()
            }])
            .unwrap();
        let mut direct = StrFromINChI {
            nNumRemovedProtonsByRevrs: 7,
            nChargeRevrs: 101,
            nChargeInChI: 103,
            ..StrFromINChI::default()
        };
        add_layer(&mut direct_heap, &mut direct, INCHI_BAS as usize, None, 3);
        let direct = direct_heap.allocate_model_storage(vec![direct]).unwrap();
        let mut direct_structures = null_structures;
        direct_structures[INCHI_BAS as usize][TAUT_YES as usize] = direct;
        let mut direct_input = InpInChI::default();
        direct_input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize] = original;
        direct_input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 1;
        direct_input.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize].nNumRemovedProtons = 3;
        assert_eq!(
            call(&mut direct_heap, 0, &direct_structures, &direct_input),
            Ok(0)
        );
        let direct_result = &direct_heap.slice(direct.as_const()).unwrap()[0];
        assert_eq!(
            (
                direct_result.nChargeRevrs,
                direct_result.nChargeInChI,
                direct_result.nRemovedProtonsByNormFromRevrs,
                direct_result.nNumRemovedProtonsByRevrs,
            ),
            (NO_VALUE_INT as i32, -4, 3, 7)
        );

        let mut deleted_heap = SourceHeap::default();
        let deleted = deleted_heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                bDeleted: 1,
                ..INChI::default()
            }])
            .unwrap();
        let mut deleted_input = InpInChI::default();
        deleted_input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize] = deleted;
        deleted_input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 1;
        assert_eq!(
            call(&mut deleted_heap, 0, &null_structures, &deleted_input,),
            Ok(0)
        );

        let mut mapped_heap = SourceHeap::default();
        let base_original = mapped_heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                nTotalCharge: 3,
                ..INChI::default()
            }])
            .unwrap();
        let rec_original = mapped_heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                nTotalCharge: 17,
                ..INChI::default()
            }])
            .unwrap();
        let mut rec = StrFromINChI {
            nNumRemovedProtonsByRevrs: 29,
            ..StrFromINChI::default()
        };
        add_layer(&mut mapped_heap, &mut rec, INCHI_BAS as usize, Some(23), 5);
        add_layer(&mut mapped_heap, &mut rec, INCHI_REC as usize, Some(19), 0);
        let rec = mapped_heap.allocate_model_storage(vec![rec]).unwrap();
        let base = mapped_heap
            .allocate_model_storage(vec![StrFromINChI {
                nLink: -1,
                ..StrFromINChI::default()
            }])
            .unwrap();
        let mut mapped_structures = null_structures;
        mapped_structures[INCHI_BAS as usize][TAUT_YES as usize] = base;
        mapped_structures[INCHI_REC as usize][TAUT_YES as usize] = rec;
        let mut mapped_input = InpInChI::default();
        mapped_input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize] = base_original;
        mapped_input.pInpInChI[INCHI_REC as usize][TAUT_YES as usize] = rec_original;
        mapped_input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 1;
        mapped_input.nNumComponents[INCHI_REC as usize][TAUT_YES as usize] = 1;
        mapped_input.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize].nNumRemovedProtons = 5;
        assert_eq!(
            call(&mut mapped_heap, 0, &mapped_structures, &mapped_input),
            Ok(0)
        );
        let base_result = &mapped_heap.slice(base.as_const()).unwrap()[0];
        assert_eq!(
            (
                base_result.nChargeRevrs,
                base_result.nChargeInChI,
                base_result.nRemovedProtonsByNormFromRevrs,
            ),
            (23, 17, 5)
        );
        let rec_result = &mapped_heap.slice(rec.as_const()).unwrap()[0];
        assert_eq!(
            (
                rec_result.nChargeRevrs,
                rec_result.nChargeInChI,
                rec_result.nRemovedProtonsByNormFromRevrs,
            ),
            (19, 17, 0)
        );

        let mut error_heap = SourceHeap::default();
        let error_original = error_heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                ..INChI::default()
            }])
            .unwrap();
        let atoms = error_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let mut error = StrFromINChI {
            at2: atoms,
            num_atoms: 1,
            bPostProcessed: 4,
            ..StrFromINChI::default()
        };
        add_layer(&mut error_heap, &mut error, INCHI_BAS as usize, None, 0);
        let error = error_heap.allocate_model_storage(vec![error]).unwrap();
        let mut error_structures = null_structures;
        error_structures[INCHI_BAS as usize][TAUT_YES as usize] = error;
        let mut error_input = InpInChI::default();
        error_input.pInpInChI[INCHI_BAS as usize][TAUT_YES as usize] = error_original;
        error_input.nNumComponents[INCHI_BAS as usize][TAUT_YES as usize] = 1;
        error_input.nNumProtons[INCHI_BAS as usize][TAUT_YES as usize].nNumRemovedIsotopicH =
            [0, 0, -1];
        assert_eq!(
            call(&mut error_heap, 0, &error_structures, &error_input),
            Ok(crate::source_types::RI_ERR_PROGR)
        );
        assert_eq!(
            error_heap.slice(error.as_const()).unwrap()[0].bPostProcessed,
            4 | crate::source_types::RI_ERR_PROGR
        );
    }

    #[test]
    fn source_port__ichirvr7__freestrfrominchi__line_1353() {
        let mut heap = SourceHeap::default();
        let borrowed = heap
            .allocate_model_storage(vec![INChI {
                nTotalCharge: 71,
                ..INChI::default()
            }])
            .unwrap();

        let mut structure = StrFromINChI {
            at: heap
                .allocate_model_storage(vec![inp_ATOM::default()])
                .unwrap(),
            at2: heap
                .allocate_model_storage(vec![inp_ATOM::default()])
                .unwrap(),
            st: heap
                .allocate_model_storage(vec![inp_ATOM_STEREO::default()])
                .unwrap(),
            pVA: heap
                .allocate_model_storage(vec![VAL_AT::default()])
                .unwrap(),
            pXYZ: heap
                .allocate_model_storage(vec![XYZ_COORD::default()])
                .unwrap(),
            endpoint: heap.allocate_model_storage(vec![1_u16]).unwrap(),
            fixed_H: heap.allocate_model_storage(vec![1_i8]).unwrap(),
            pOneINChI: [borrowed, SourceMutPointer::null()],
            ..StrFromINChI::default()
        };
        structure.One_ti.t_group = heap
            .allocate_model_storage(vec![T_GROUP::default()])
            .unwrap();
        structure.ti.t_group = heap
            .allocate_model_storage(vec![T_GROUP::default()])
            .unwrap();
        structure.ti.nEndpointAtomNumber = heap.allocate_model_storage(vec![2_u16]).unwrap();
        structure.ti.tGroupNumber = heap.allocate_model_storage(vec![3_u16]).unwrap();
        structure.ti.nIsotopicEndpointAtomNumber =
            heap.allocate_model_storage(vec![4_u16]).unwrap();
        for index in 0..2 {
            structure.nAtno2Canon[index] =
                heap.allocate_model_storage(vec![index as u16 + 5]).unwrap();
            structure.nCanon2Atno[index] =
                heap.allocate_model_storage(vec![index as u16 + 7]).unwrap();
        }

        let reversed_inchi = heap.allocate_model_storage(vec![INChI::default()]).unwrap();
        structure.RevInChI.pINChI[INCHI_BAS as usize] = heap
            .allocate_model_storage(vec![[SourceMutPointer::null(), reversed_inchi]])
            .unwrap();
        let reversed_aux = heap
            .allocate_model_storage(vec![INChI_Aux::default()])
            .unwrap();
        structure.RevInChI.pINChI_Aux[INCHI_BAS as usize] = heap
            .allocate_model_storage(vec![[SourceMutPointer::null(), reversed_aux]])
            .unwrap();
        structure.RevInChI.num_components[INCHI_BAS as usize] = 1;

        let owned = heap
            .allocate_model_storage(vec![structure, StrFromINChI::default()])
            .unwrap();
        let retained_zero_count = heap
            .allocate_model_storage(vec![StrFromINChI::default()])
            .unwrap();
        let freed_negative_count = heap
            .allocate_model_storage(vec![StrFromINChI::default()])
            .unwrap();
        let mut structures = [
            [owned, retained_zero_count],
            [freed_negative_count, SourceMutPointer::null()],
        ];
        let counts = [[2, 0], [-1, 1]];

        FreeStrFromINChI(&mut heap, &mut structures, &counts).unwrap();
        assert_eq!(structures[0][0], SourceMutPointer::null());
        assert_eq!(structures[0][1], retained_zero_count);
        assert_eq!(structures[1][0], SourceMutPointer::null());
        assert_eq!(structures[1][1], SourceMutPointer::null());
        assert_eq!(heap.slice(borrowed.as_const()).unwrap()[0].nTotalCharge, 71);
        assert_eq!(
            heap.slice(retained_zero_count.as_const()).unwrap(),
            &[StrFromINChI::default()]
        );
        assert_eq!(heap.live_allocation_count(), 2);

        let mut no_op_heap = SourceHeap::default();
        let retained = no_op_heap
            .allocate_model_storage(vec![StrFromINChI::default()])
            .unwrap();
        let mut no_op = [[SourceMutPointer::null(); 2]; 2];
        no_op[1][1] = retained;
        FreeStrFromINChI(&mut no_op_heap, &mut no_op, &[[0; 2]; 2]).unwrap();
        assert_eq!(no_op[1][1], retained);
        assert_eq!(no_op_heap.live_allocation_count(), 1);
    }

    #[test]
    fn source_port__ichirvr7__freeinpinchi__line_1438() {
        let mut heap = SourceHeap::default();
        let mut empty = InpInChI::default();
        FreeInpInChI(&mut heap, &mut empty).unwrap();
        assert_eq!(empty, InpInChI::default());

        let hill = heap.allocate(vec![b'C' as i8, 0]).unwrap();
        let atom_member = heap.allocate(vec![6_u8]).unwrap();
        let connection = heap.allocate(vec![1_u16]).unwrap();
        let tautomer = heap.allocate(vec![2_u16]).unwrap();
        let num_h = heap.allocate(vec![1_i8]).unwrap();
        let num_h_fixed = heap.allocate(vec![2_i8]).unwrap();
        let possible_h = heap.allocate(vec![3_u16]).unwrap();
        let isotopic_atom = heap.allocate(vec![INChI_IsotopicAtom::default()]).unwrap();
        let isotopic_t_group = heap
            .allocate(vec![INChI_IsotopicTGroup::default()])
            .unwrap();
        let stereo = heap.allocate(vec![INChI_Stereo::default()]).unwrap();
        let stereo_isotopic = heap.allocate(vec![INChI_Stereo::default()]).unwrap();
        let shared = INChI {
            szHillFormula: hill,
            nAtom: atom_member,
            nConnTable: connection,
            nTautomer: tautomer,
            nNum_H: num_h,
            nNum_H_fixed: num_h_fixed,
            IsotopicAtom: isotopic_atom,
            IsotopicTGroup: isotopic_t_group,
            Stereo: stereo,
            StereoIsotopic: stereo_isotopic,
            nPossibleLocationsOfIsotopicH: possible_h,
            ..INChI::default()
        };
        let components = heap.allocate(vec![shared.clone(), shared]).unwrap();
        let second_components = heap.allocate(vec![INChI::default()]).unwrap();
        let proton_00 = heap
            .allocate(vec![COMPONENT_REM_PROTONS::default()])
            .unwrap();
        let proton_01 = heap
            .allocate(vec![COMPONENT_REM_PROTONS::default()])
            .unwrap();
        let proton_10 = heap
            .allocate(vec![COMPONENT_REM_PROTONS::default()])
            .unwrap();
        let proton_11 = heap
            .allocate(vec![COMPONENT_REM_PROTONS::default()])
            .unwrap();
        let atoms = heap.allocate(vec![inp_ATOM::default()]).unwrap();
        let polymer = heap.allocate(vec![OAD_Polymer::default()]).unwrap();
        let v3000 = heap.allocate(vec![OAD_V3000::default()]).unwrap();

        let mut input = InpInChI::default();
        input.pInpInChI[0][0] = components;
        input.nNumComponents[0][0] = 2;
        input.pInpInChI[1][1] = second_components;
        input.nNumComponents[1][1] = 1;
        input.nNumProtons[0][0].pNumProtons = proton_00;
        input.nNumProtons[0][1].pNumProtons = proton_01;
        input.nNumProtons[1][0].pNumProtons = proton_10;
        input.nNumProtons[1][1].pNumProtons = proton_11;
        input.atom = atoms;
        input.num_atoms = 1;
        input.num_inp = 77;
        input.polymer = polymer;
        input.v3000 = v3000;

        FreeInpInChI(&mut heap, &mut input).unwrap();
        assert_eq!(input, InpInChI::default());
        for freed in [
            heap.slice(hill.as_const()).map(|_| ()),
            heap.slice(atom_member.as_const()).map(|_| ()),
            heap.slice(connection.as_const()).map(|_| ()),
            heap.slice(tautomer.as_const()).map(|_| ()),
            heap.slice(num_h.as_const()).map(|_| ()),
            heap.slice(num_h_fixed.as_const()).map(|_| ()),
            heap.slice(possible_h.as_const()).map(|_| ()),
            heap.slice(isotopic_atom.as_const()).map(|_| ()),
            heap.slice(isotopic_t_group.as_const()).map(|_| ()),
            heap.slice(stereo.as_const()).map(|_| ()),
            heap.slice(stereo_isotopic.as_const()).map(|_| ()),
            heap.slice(components.as_const()).map(|_| ()),
            heap.slice(second_components.as_const()).map(|_| ()),
            heap.slice(proton_00.as_const()).map(|_| ()),
            heap.slice(proton_01.as_const()).map(|_| ()),
            heap.slice(proton_10.as_const()).map(|_| ()),
            heap.slice(proton_11.as_const()).map(|_| ()),
            heap.slice(atoms.as_const()).map(|_| ()),
            heap.slice(polymer.as_const()).map(|_| ()),
            heap.slice(v3000.as_const()).map(|_| ()),
        ] {
            assert_eq!(freed, Err(SourceHeapError::MissingAllocation));
        }
    }
}
