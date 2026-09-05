use std::mem::size_of;

use crate::source::base::ichicans::{
    FillOutStereoParities, InvertStereo, SetCtToIsotopicStereo, SetCtToNonIsotopicStereo,
    SwitchAtomStereoAndIsotopicStereo,
};
use crate::source::base::ichimap1::{CompareLinCtStereo, CurTreeAlloc, CurTreeFree, CurTreeSetPos};
use crate::source::base::ichimap2::{SortedEquInfoToRanks, SortedRanksToEquInfo};
use crate::source::base::ichimap4::map_stereo_bonds4;
use crate::source::base::ichisort::{
    CompNeighborsAT_NUMBER, CompNeighborsATNumberContext, CompRank, CompRanksOrd, FreeNeighList,
    inchi_qsort, insertions_sort,
};
use crate::source::base::ichitaut::SortTautomerGroupsAndEndpoints;
use crate::source::base::util::{inchi_calloc, inchi_free};
use crate::source_types::{
    AB_PARITY_UNDF, AB_PARITY_UNKN, AT_FLAG_ISO_H_POINT, AT_ISO_TGROUP, AT_ISOTOPIC, AT_NUMB,
    AT_RANK, AT_TAUTOMER, ATOM_SIZES, CANON_GLOBALS, CANON_STAT, CMODE_ISO, CMODE_ISO_STEREO,
    CMODE_NOEQ_STEREO, CMODE_STEREO, CT_ATOMID, CT_ATOMID_DONTINCLUDE, CT_CALC_STEREO_ERR,
    CT_CANON_ERR, CT_ERR_MAX, CT_ERR_MIN, CT_LEN_MISMATCH, CT_OUT_OF_RAM, CT_OVERFLOW,
    CT_STEREOCOUNT_ERR, CUR_TREE, INCHI_CLOCK, MAX_NUM_STEREO_BONDS, MAXVAL,
    REQ_MODE_DIFF_UU_STEREO, SourceConstPointer, SourceHeap, SourceHeapError, SourceMutPointer,
    T_GROUP_HDR_LEN, T_GROUP_INFO, T_NUM_NO_ISOTOPIC, TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE,
    TG_FLAG_FOUND_ISOTOPIC_H_DONE, TGSO_SYMM_IORDER, TGSO_SYMM_IRANK, TGSO_SYMM_RANK, clock_t,
    inchiTime, sp_ATOM,
};

fn source_copy<T: Clone + 'static>(
    heap: &mut SourceHeap,
    destination: SourceMutPointer<T>,
    source: SourceConstPointer<T>,
    count: usize,
) -> Result<(), SourceHeapError> {
    heap.with_slice_mut_and_heap(destination, |destination, heap| {
        let destination = destination
            .get_mut(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let source = heap
            .slice(source)?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        destination.clone_from_slice(source);
        Ok(())
    })
}

fn source_returned_error(value: i32) -> bool {
    CT_ERR_MIN <= value && value <= CT_ERR_MAX
}

#[allow(non_snake_case)]
pub(crate) fn InchiClock(clock_result: clock_t) -> clock_t {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:111 InchiClock
    // INCHI✔❌: static clock_t InchiClock( void )
    // INCHI✔❌: {
    // INCHI✔❌:     clock_t c = clock( );
    // INCHI✔❌:     if (c != (clock_t) -1)
    // INCHI✔❌:     {
    // INCHI✔❌:         return c;
    // INCHI✔❌:     }
    // INCHI✔❌:     return 0;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: InchiClock

    // clock() is a libc environment boundary. Its clock_t result is injected
    // explicitly so production code does not replace process CPU time with a
    // wall-clock heuristic or a native FFI call.
    if clock_result != -1 { clock_result } else { 0 }
}

#[allow(non_snake_case)]
pub(crate) fn InchiTimeGet(tick_end: &mut inchiTime, clock_result: clock_t) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:146 InchiTimeGet
    // INCHI✔❌: void InchiTimeGet( inchiTime *TickEnd )
    // INCHI✔❌: {
    // INCHI✔❌:     TickEnd->clockTime = InchiClock( );
    // INCHI✔❌: }
    // END INCHI C FUNCTION: InchiTimeGet

    tick_end.clockTime = InchiClock(clock_result);
}

#[allow(non_snake_case)]
pub(crate) fn FillMaxMinClock(ic: &mut INCHI_CLOCK) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:128 FillMaxMinClock
    // INCHI✔❌: static void FillMaxMinClock( INCHI_CLOCK *ic )
    // INCHI✔❌: {
    // INCHI✔❌:     if (!ic->m_MaxPositiveClock)
    // INCHI✔❌:     {
    // INCHI✔❌:         clock_t valPos = 0, val1 = 1;
    // INCHI✔❌:         while (0 < ( ( val1 <<= 1 ), ( val1 |= 1 ) )) /* djb-rwth: ignoring GH issue #59.3 -- LLONG_MIN/LLONG_MAX not found in <limits.h> on Linux */
    // INCHI✔❌:         {
    // INCHI✔❌:             valPos = val1;
    // INCHI✔❌:         }
    // INCHI✔❌:         ic->m_MaxPositiveClock = valPos;
    // INCHI✔❌:         ic->m_MinNegativeClock = -valPos;
    // INCHI✔❌:         ic->m_HalfMaxPositiveClock = ic->m_MaxPositiveClock / 2;
    // INCHI✔❌:         ic->m_HalfMinNegativeClock = ic->m_MinNegativeClock / 2;
    // INCHI✔❌:     }
    // INCHI✔❌: }
    // END INCHI C FUNCTION: FillMaxMinClock

    if ic.m_MaxPositiveClock == 0 {
        let mut val_pos = 0_i64;
        let mut val1 = 1_i64;
        loop {
            val1 = val1.wrapping_shl(1);
            val1 |= 1;
            if val1 <= 0 {
                break;
            }
            val_pos = val1;
        }
        ic.m_MaxPositiveClock = val_pos;
        ic.m_MinNegativeClock = -val_pos;
        ic.m_HalfMaxPositiveClock = ic.m_MaxPositiveClock / 2;
        ic.m_HalfMinNegativeClock = ic.m_MinNegativeClock / 2;
    }
}

#[allow(non_snake_case)]
pub(crate) fn InchiTimeMsecDiff(
    ic: &mut INCHI_CLOCK,
    tick_end: Option<&inchiTime>,
    tick_start: Option<&inchiTime>,
) -> i64 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:151 InchiTimeMsecDiff
    // INCHI✔❌: long InchiTimeMsecDiff( INCHI_CLOCK *ic, inchiTime *TickEnd, inchiTime *TickStart )
    // INCHI✔❌: {
    // INCHI✔❌:     if (FullMaxClock > 0)
    // INCHI✔❌:     {
    // INCHI✔❌:         clock_t delta;
    // INCHI✔❌:         if (!TickEnd || !TickStart)
    // INCHI✔❌:             return 0;
    // INCHI✔❌:         /* clock_t is unsigned */
    // INCHI✔❌:         if (TickEnd->clockTime > TickStart->clockTime)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (TickEnd->clockTime > HalfMaxClock &&
    // INCHI✔❌:                  TickEnd->clockTime - TickStart->clockTime > HalfMaxClock)
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* overflow in TickStart->clockTime, actually TickStart->clockTime was later */
    // INCHI✔❌:                 delta = ( FullMaxClock - TickEnd->clockTime ) + TickStart->clockTime;
    // INCHI✔❌:                 return -INCHI_MSEC( delta );
    // INCHI✔❌:             }
    // INCHI✔❌:             delta = TickEnd->clockTime - TickStart->clockTime;
    // INCHI✔❌:             return INCHI_MSEC( delta );
    // INCHI✔❌:         }
    // INCHI✔❌:         else if (TickEnd->clockTime < TickStart->clockTime)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (TickStart->clockTime > HalfMaxClock &&
    // INCHI✔❌:                  TickStart->clockTime - TickEnd->clockTime > HalfMaxClock)
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* overflow in TickEnd->clockTime, actually TickEnd->clockTime was later */
    // INCHI✔❌:                 delta = ( FullMaxClock - TickStart->clockTime ) + TickEnd->clockTime;
    // INCHI✔❌:                 return INCHI_MSEC( delta );
    // INCHI✔❌:             }
    // INCHI✔❌:             delta = TickStart->clockTime - TickEnd->clockTime;
    // INCHI✔❌:             return -INCHI_MSEC( delta );
    // INCHI✔❌:         }
    // INCHI✔❌:         return 0; /* TickEnd->clockTime == TickStart->clockTime */
    // INCHI✔❌:     }
    // INCHI✔❌:     else
    // INCHI✔❌:     {
    // INCHI✔❌:         /* may happen under Win32 only where clock_t is SIGNED long */
    // INCHI✔❌:         clock_t delta;
    // INCHI✔❌:         FillMaxMinClock( ic );
    // INCHI✔❌:         if (!TickEnd || !TickStart)
    // INCHI✔❌:             return 0;
    // INCHI✔❌:         if ((TickEnd->clockTime >= 0 && TickStart->clockTime >= 0) ||
    // INCHI✔❌:              (TickEnd->clockTime <= 0 && TickStart->clockTime <= 0)) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔❌:         {
    // INCHI✔❌:             delta = TickEnd->clockTime - TickStart->clockTime;
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:             if (TickEnd->clockTime >= ic->m_HalfMaxPositiveClock &&
    // INCHI✔❌:                  TickStart->clockTime <= ic->m_HalfMinNegativeClock)
    // INCHI✔❌:             {
    // INCHI✔❌: /* end is earlier than start */
    // INCHI✔❌:                 delta = ( ic->m_MaxPositiveClock - TickEnd->clockTime ) + ( TickStart->clockTime - ic->m_MinNegativeClock );
    // INCHI✔❌:                 delta = -delta;
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:                 if (TickEnd->clockTime <= ic->m_HalfMinNegativeClock &&
    // INCHI✔❌:                      TickStart->clockTime >= ic->m_HalfMaxPositiveClock)
    // INCHI✔❌:                 {
    // INCHI✔❌: /* start was earlier than end */
    // INCHI✔❌:                     delta = ( ic->m_MaxPositiveClock - TickStart->clockTime ) + ( TickEnd->clockTime - ic->m_MinNegativeClock );
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                      /* there was no overflow, clock passed zero */
    // INCHI✔❌:                     delta = TickEnd->clockTime - TickStart->clockTime;
    // INCHI✔❌:                 }
    // INCHI✔❌:         return INCHI_MSEC( delta );
    // INCHI✔❌:     }
    // INCHI✔❌: }
    // END INCHI C FUNCTION: InchiTimeMsecDiff

    FillMaxMinClock(ic);
    let (Some(tick_end), Some(tick_start)) = (tick_end, tick_start) else {
        return 0;
    };
    let end = tick_end.clockTime;
    let start = tick_start.clockTime;
    let delta = if (end >= 0 && start >= 0) || (end <= 0 && start <= 0) {
        end - start
    } else if end >= ic.m_HalfMaxPositiveClock && start <= ic.m_HalfMinNegativeClock {
        -((ic.m_MaxPositiveClock - end) + (start - ic.m_MinNegativeClock))
    } else if end <= ic.m_HalfMinNegativeClock && start >= ic.m_HalfMaxPositiveClock {
        (ic.m_MaxPositiveClock - start) + (end - ic.m_MinNegativeClock)
    } else {
        end - start
    };
    ((1000.0 / 1_000_000.0) * delta as f64) as i64
}

#[allow(non_snake_case)]
pub(crate) fn InchiTimeElapsed(
    ic: &mut INCHI_CLOCK,
    tick_start: Option<&inchiTime>,
    clock_result: clock_t,
) -> i64 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:223 InchiTimeElapsed
    // INCHI✔❌: long InchiTimeElapsed( INCHI_CLOCK *ic, inchiTime *TickStart )
    // INCHI✔❌: {
    // INCHI✔❌:     inchiTime TickEnd;
    // INCHI✔❌:     if (!TickStart)
    // INCHI✔❌:         return 0;
    // INCHI✔❌:     InchiTimeGet( &TickEnd );
    // INCHI✔❌:     return InchiTimeMsecDiff( ic, &TickEnd, TickStart );
    // INCHI✔❌: }
    // END INCHI C FUNCTION: InchiTimeElapsed

    let Some(tick_start) = tick_start else {
        return 0;
    };
    let mut tick_end = inchiTime::default();
    InchiTimeGet(&mut tick_end, clock_result);
    InchiTimeMsecDiff(ic, Some(&tick_end), Some(tick_start))
}

#[allow(non_snake_case)]
pub(crate) fn InchiTimeAddMsec(
    ic: &mut INCHI_CLOCK,
    tick_end: Option<&mut inchiTime>,
    number_msec: u64,
) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:234 InchiTimeAddMsec
    // INCHI✔️✔️: void InchiTimeAddMsec( INCHI_CLOCK *ic, inchiTime *TickEnd, unsigned long nNumMsec )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     clock_t delta;
    // INCHI✔️✔️:     if (!TickEnd)
    // INCHI✔️✔️:         return;
    // INCHI✔️✔️:     if (FullMaxClock > 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* clock_t is unsigned */
    // INCHI✔️✔️:         delta = INCHI_CLOCK_T( nNumMsec );
    // INCHI✔️✔️:         TickEnd->clockTime += delta;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* may happen under Win32 only where clock_t is SIGNED long */
    // INCHI✔️✔️:         /* clock_t is unsigned */
    // INCHI✔️✔️:         FillMaxMinClock( ic );
    // INCHI✔️✔️:         delta = INCHI_CLOCK_T( nNumMsec );
    // INCHI✔️✔️:         TickEnd->clockTime += delta;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: InchiTimeAddMsec
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: InchiTimeAddMsec
    // INCHI✔️✔️: #define INCHI_CLOCK_T(X) (clock_t)( (double)(X) / 1000.0 * (double)CLOCKS_PER_SEC )
    // INCHI✔️✔️: CLOCKS_PER_SEC == 1000000; sizeof(unsigned long) == 8; clock_t == signed i64; FullMaxClock == -1
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: InchiTimeAddMsec

    let Some(tick_end) = tick_end else {
        return;
    };
    FillMaxMinClock(ic);
    let delta = (number_msec as f64 / 1000.0 * 1_000_000.0) as clock_t;
    tick_end.clockTime = tick_end.clockTime.wrapping_add(delta);
}

#[allow(non_snake_case)]
pub(crate) fn bInchiTimeIsOver(
    ic: &mut INCHI_CLOCK,
    tick_start: Option<&inchiTime>,
    clock_result: clock_t,
) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:257 bInchiTimeIsOver
    // INCHI✔❌: int bInchiTimeIsOver( INCHI_CLOCK *ic, inchiTime *TickStart )
    // INCHI✔❌: {
    // INCHI✔❌:     if (FullMaxClock > 0)
    // INCHI✔❌:     {
    // INCHI✔❌:         clock_t clockCurrTime;
    // INCHI✔❌:         if (!TickStart)
    // INCHI✔❌:             return 0;
    // INCHI✔❌:         clockCurrTime = InchiClock( );
    // INCHI✔❌:         /* clock_t is unsigned */
    // INCHI✔❌:         if (TickStart->clockTime > clockCurrTime)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (TickStart->clockTime > HalfMaxClock &&
    // INCHI✔❌:                  TickStart->clockTime - clockCurrTime > HalfMaxClock)
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* overflow in clockCurrTime, actually clockCurrTime was later */
    // INCHI✔❌:                 return 1;
    // INCHI✔❌:             }
    // INCHI✔❌:             return 0;
    // INCHI✔❌:         }
    // INCHI✔❌:         else if (TickStart->clockTime < clockCurrTime)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (clockCurrTime > HalfMaxClock &&
    // INCHI✔❌:                  clockCurrTime - TickStart->clockTime > HalfMaxClock)
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* overflow in TickStart->clockTime, actually TickStart->clockTime was later */
    // INCHI✔❌:                 return 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             return 1;
    // INCHI✔❌:         }
    // INCHI✔❌:         return 0; /* TickStart->clockTime == clockCurrTime */
    // INCHI✔❌:     }
    // INCHI✔❌:     else
    // INCHI✔❌:     {
    // INCHI✔❌:         /* may happen under Win32 only where clock_t is SIGNED long */
    // INCHI✔❌:         clock_t clockCurrTime;
    // INCHI✔❌:         FillMaxMinClock( ic );
    // INCHI✔❌:         if (!TickStart)
    // INCHI✔❌:             return 0;
    // INCHI✔❌:         clockCurrTime = InchiClock( );
    // INCHI✔❌:         if ((clockCurrTime >= 0 && TickStart->clockTime >= 0) ||
    // INCHI✔❌:              (clockCurrTime <= 0 && TickStart->clockTime <= 0)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:         {
    // INCHI✔❌:             return ( clockCurrTime > TickStart->clockTime );
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:             if (clockCurrTime >= ic->m_HalfMaxPositiveClock &&
    // INCHI✔❌:                  TickStart->clockTime <= ic->m_HalfMinNegativeClock)
    // INCHI✔❌:             {
    // INCHI✔❌: /* curr is earlier than start */
    // INCHI✔❌:                 return 0;
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:                 if (clockCurrTime <= ic->m_HalfMinNegativeClock &&
    // INCHI✔❌:                      TickStart->clockTime >= ic->m_HalfMaxPositiveClock)
    // INCHI✔❌:                 {
    // INCHI✔❌: /* start was earlier than curr */
    // INCHI✔❌:                     return 1;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                      /* there was no overflow, clock passed zero */
    // INCHI✔❌:                     return ( clockCurrTime > TickStart->clockTime );
    // INCHI✔❌:                 }
    // INCHI✔❌:     }
    // INCHI✔❌: }
    // END INCHI C FUNCTION: bInchiTimeIsOver

    FillMaxMinClock(ic);
    let Some(tick_start) = tick_start else {
        return 0;
    };
    let clock_curr_time = InchiClock(clock_result);
    let start = tick_start.clockTime;
    if (clock_curr_time >= 0 && start >= 0) || (clock_curr_time <= 0 && start <= 0) {
        (clock_curr_time > start) as i32
    } else if clock_curr_time >= ic.m_HalfMaxPositiveClock && start <= ic.m_HalfMinNegativeClock {
        0
    } else if clock_curr_time <= ic.m_HalfMinNegativeClock && start >= ic.m_HalfMaxPositiveClock {
        1
    } else {
        (clock_curr_time > start) as i32
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn FillIsotopicAtLinearCT(
    heap: &mut SourceHeap,
    num_atoms: i32,
    at: SourceConstPointer<sp_ATOM>,
    nAtomNumber: SourceConstPointer<AT_RANK>,
    LinearCTIsotopic: SourceMutPointer<AT_ISOTOPIC>,
    nMaxLenLinearCTIsotopic: i32,
    pnLenLinearCTIsotopic: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:780 FillIsotopicAtLinearCT
    // INCHI✔️❌: int FillIsotopicAtLinearCT( int num_atoms,
    // INCHI✔️❌:                             sp_ATOM* at,
    // INCHI✔️❌:                             const AT_RANK *nAtomNumber,
    // INCHI✔️❌:                             AT_ISOTOPIC *LinearCTIsotopic,
    // INCHI✔️❌:                             int nMaxLenLinearCTIsotopic,
    // INCHI✔️❌:                             int *pnLenLinearCTIsotopic )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /* at[i].init_rank = initial ranks before canonizing */
    // INCHI✔️❌:     /* nRank[i]  = new ordering number for atoms: nRank=1,2,.. */
    // INCHI✔️❌:     /* nAtomNumber[r] = orig. atom number= 0,1,...  for r = nRank-1  */
    // INCHI✔️❌:     /* nRank[nAtomNumber[r]] = r; r = 0,1,... */
    // INCHI✔️❌:     /* nAtomNumber[nRank[i]-1] = i; */
    // INCHI✔️❌:
    // INCHI✔️❌:     int  i, k, rank;
    // INCHI✔️❌:     int  nLinearCTIsotopicLen = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* the following parts of the "name" should be compared */
    // INCHI✔️❌:     /* after the connection table comparison is done */
    // INCHI✔️❌:     /* to avoid wrong difference sign. So, these parts */
    // INCHI✔️❌:     /* go to a separate buffers. */
    // INCHI✔️❌:     if (LinearCTIsotopic && nMaxLenLinearCTIsotopic > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         memset( LinearCTIsotopic, 0, nMaxLenLinearCTIsotopic * sizeof( LinearCTIsotopic[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* rank = nRank[nAtomNumber[rank-1]] -- proposed atoms canon. numbers */
    // INCHI✔️❌:     for (rank = 1; rank <= num_atoms; rank++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         i = (int) nAtomNumber[rank - 1];  /* current atom */
    // INCHI✔️❌:
    // INCHI✔️❌:         /****************************************************
    // INCHI✔️❌:              add isotopic atom info to LinearCTIsotopic
    // INCHI✔️❌:         *****************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:         /* if the atom itself is not isotopic then add it only if */
    // INCHI✔️❌:         /* the atom is not an endpoint AND has attached T or D or 1H.  */
    // INCHI✔️❌:         k = ( !at[i].endpoint && !( at[i].cFlags & AT_FLAG_ISO_H_POINT ) && ( at[i].num_iso_H[0] || at[i].num_iso_H[1] || at[i].num_iso_H[2] ) );
    // INCHI✔️❌:         if (at[i].iso_atw_diff || k)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (CHECK_OVERFLOW( nLinearCTIsotopicLen, nMaxLenLinearCTIsotopic ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return CT_OVERFLOW;  /*  <BRKPT> */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             LinearCTIsotopic[nLinearCTIsotopicLen].at_num = (AT_RANK) rank;
    // INCHI✔️❌:             LinearCTIsotopic[nLinearCTIsotopicLen].iso_atw_diff = at[i].iso_atw_diff;
    // INCHI✔️❌:             LinearCTIsotopic[nLinearCTIsotopicLen].num_1H = (NUM_H) ( k ? at[i].num_iso_H[0] : 0 );
    // INCHI✔️❌:             LinearCTIsotopic[nLinearCTIsotopicLen].num_D = (NUM_H) ( k ? at[i].num_iso_H[1] : 0 );
    // INCHI✔️❌:             LinearCTIsotopic[nLinearCTIsotopicLen].num_T = (NUM_H) ( k ? at[i].num_iso_H[2] : 0 );
    // INCHI✔️❌:             nLinearCTIsotopicLen++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     } /* end of cycle over all atoms. */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (LinearCTIsotopic)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (*pnLenLinearCTIsotopic)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (*pnLenLinearCTIsotopic != nLinearCTIsotopicLen)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return CT_LEN_MISMATCH;  /*  <BRKPT> */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *pnLenLinearCTIsotopic = nLinearCTIsotopicLen;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Return value: >0 => OK */
    // INCHI✔️❌:     return nLinearCTIsotopicLen;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: FillIsotopicAtLinearCT

    if LinearCTIsotopic.is_null() || nMaxLenLinearCTIsotopic <= 0 {
        return Ok(0);
    }
    let capacity = usize::try_from(nMaxLenLinearCTIsotopic)
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    heap.slice_mut(LinearCTIsotopic)?
        .get_mut(..capacity)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .fill(AT_ISOTOPIC::default());

    let mut output_len = 0_i32;
    if num_atoms > 0 {
        for rank in 1..=num_atoms {
            let rank_index = usize::try_from(rank - 1).unwrap();
            let atom_number = usize::from(
                *heap
                    .slice(nAtomNumber)?
                    .get(rank_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            let atom = heap
                .slice(at)?
                .get(atom_number)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let has_isotopic_h = atom.endpoint == 0
                && atom.cFlags & AT_FLAG_ISO_H_POINT as i8 == 0
                && atom.num_iso_H.iter().any(|&count| count != 0);
            if atom.iso_atw_diff != 0 || has_isotopic_h {
                if output_len >= nMaxLenLinearCTIsotopic {
                    return Ok(CT_OVERFLOW);
                }
                let record = heap
                    .slice_mut(LinearCTIsotopic)?
                    .get_mut(usize::try_from(output_len).unwrap())
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                record.at_num = rank as AT_RANK;
                record.iso_atw_diff = atom.iso_atw_diff as _;
                record.num_1H = if has_isotopic_h {
                    atom.num_iso_H[0] as _
                } else {
                    0
                };
                record.num_D = if has_isotopic_h {
                    atom.num_iso_H[1] as _
                } else {
                    0
                };
                record.num_T = if has_isotopic_h {
                    atom.num_iso_H[2] as _
                } else {
                    0
                };
                output_len = output_len.wrapping_add(1);
            }
        }
    }
    if *pnLenLinearCTIsotopic != 0 {
        if *pnLenLinearCTIsotopic != output_len {
            return Ok(CT_LEN_MISMATCH);
        }
    } else {
        *pnLenLinearCTIsotopic = output_len;
    }
    Ok(output_len)
}

#[allow(non_snake_case)]
pub(crate) fn GetCanonLengths(
    heap: &SourceHeap,
    num_at: i32,
    at: SourceConstPointer<sp_ATOM>,
    s: &mut ATOM_SIZES,
    t_group_info: Option<&T_GROUP_INFO>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:418 GetCanonLengths
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int GetCanonLengths( int num_at,
                         sp_ATOM* at,
                         ATOM_SIZES *s,
                         T_GROUP_INFO *t_group_info )
    {
        /* include taut. groups as additional "atoms" to the connection table 07-22-2002 */

        int  i, nNumCT, nNumBonds, nNumTBonds = 0, nNumDblBondsStereo = 0, nNumAsymCarbStereo = 0, nNumIsotopic = 0;
        T_GROUP *t_group = ( s->nLenLinearCTTautomer && t_group_info ) ? t_group_info->t_group : NULL;

        for (nNumBonds = 0, i = 0; i < num_at; i++)
        {
            nNumBonds += at[i].valence;
            if (at[i].iso_sort_key)
            {
                nNumIsotopic++;  /* not including tautomeric endpoints that are isotopic only due to mobile atoms */
            }

            if (at[i].parity > 0)
            {
                /* ignore hydrogen isotope parities in at[i].parity2 */
                int j = 0, nStereoBondsToAtom = 0;  /* number of stereo double bonds at this atom */
                int k;
                for (; j < MAX_NUM_STEREO_BONDS && ( k = at[i].stereo_bond_neighbor[j] ); j++)
                {
                    nStereoBondsToAtom += ( at[k - 1].parity > 0 );
                }
                nNumDblBondsStereo += nStereoBondsToAtom;
                nNumAsymCarbStereo += !j;
            }
        }

        nNumDblBondsStereo /= 2;
        nNumBonds /= 2;

        s->nLenBonds = inchi_max( s->nLenBonds, nNumBonds );
        nNumCT = nNumBonds; /* total number of neighbors in the CT */

    #if ( CT_ATOMID != CT_ATOMID_DONTINCLUDE )
        nNumCT += num_at;
    #endif

        s->nLenCTAtOnly = inchi_max( s->nLenCTAtOnly, nNumCT );

        if (t_group)
        {
            for (i = 0; i < t_group_info->num_t_groups; i++)
            {
                nNumTBonds += t_group[i].nNumEndpoints;
            }
            nNumCT += nNumTBonds;

    #if ( CT_ATOMID != CT_ATOMID_DONTINCLUDE )
            nNumCT += t_group_info->num_t_groups;
    #endif
        }

        nNumCT = inchi_max( 1, nNumCT ); /* keep GetBaseCanonRanking() happy */
        s->nLenCT = inchi_max( s->nLenCT, nNumCT );
        s->nLenIsotopic = inchi_max( s->nLenIsotopic, nNumIsotopic );
        s->nLenLinearCTStereoDble = inchi_max( s->nLenLinearCTStereoDble, nNumDblBondsStereo );
        s->nLenLinearCTStereoCarb = inchi_max( s->nLenLinearCTStereoCarb, nNumAsymCarbStereo );

        if (t_group_info)
        {
            s->nLenIsotopicEndpoints = inchi_max( s->nLenIsotopicEndpoints, t_group_info->nNumIsotopicEndpoints );
        }

        return 0;
    }
    */
    // END INCHI C FUNCTION: GetCanonLengths

    let atom_count = usize::try_from(num_at).unwrap_or(0);
    let atoms = if atom_count == 0 {
        &[][..]
    } else {
        heap.slice(at)?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
    };
    let mut n_num_bonds = 0_i32;
    let mut n_num_dbl_bonds_stereo = 0_i32;
    let mut n_num_asym_carb_stereo = 0_i32;
    let mut n_num_isotopic = 0_i32;
    for atom in atoms {
        n_num_bonds = n_num_bonds.wrapping_add(i32::from(atom.valence));
        if atom.iso_sort_key != 0 {
            n_num_isotopic = n_num_isotopic.wrapping_add(1);
        }
        if atom.parity > 0 {
            let mut j = 0_usize;
            let mut n_stereo_bonds_to_atom = 0_i32;
            while j < MAX_NUM_STEREO_BONDS as usize {
                let opposite = atom.stereo_bond_neighbor[j];
                if opposite == 0 {
                    break;
                }
                let opposite_atom = atoms
                    .get(usize::from(opposite - 1))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                n_stereo_bonds_to_atom =
                    n_stereo_bonds_to_atom.wrapping_add(i32::from(opposite_atom.parity > 0));
                j += 1;
            }
            n_num_dbl_bonds_stereo = n_num_dbl_bonds_stereo.wrapping_add(n_stereo_bonds_to_atom);
            n_num_asym_carb_stereo = n_num_asym_carb_stereo.wrapping_add(i32::from(j == 0));
        }
    }
    n_num_dbl_bonds_stereo /= 2;
    n_num_bonds /= 2;

    s.nLenBonds = s.nLenBonds.max(n_num_bonds);
    let mut n_num_ct = n_num_bonds;
    if CT_ATOMID != CT_ATOMID_DONTINCLUDE {
        n_num_ct = n_num_ct.wrapping_add(num_at);
    }
    s.nLenCTAtOnly = s.nLenCTAtOnly.max(n_num_ct);

    let t_group_pointer = if s.nLenLinearCTTautomer != 0 {
        t_group_info.map(|info| info.t_group)
    } else {
        None
    };
    if let (Some(info), Some(pointer)) = (t_group_info, t_group_pointer)
        && !pointer.is_null()
    {
        let group_count = usize::try_from(info.num_t_groups).unwrap_or(0);
        let groups = heap
            .slice(pointer.as_const())?
            .get(..group_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut n_num_t_bonds = 0_i32;
        for group in groups {
            n_num_t_bonds = n_num_t_bonds.wrapping_add(i32::from(group.nNumEndpoints));
        }
        n_num_ct = n_num_ct.wrapping_add(n_num_t_bonds);
        if CT_ATOMID != CT_ATOMID_DONTINCLUDE {
            n_num_ct = n_num_ct.wrapping_add(info.num_t_groups);
        }
    }

    n_num_ct = n_num_ct.max(1);
    s.nLenCT = s.nLenCT.max(n_num_ct);
    s.nLenIsotopic = s.nLenIsotopic.max(n_num_isotopic);
    s.nLenLinearCTStereoDble = s.nLenLinearCTStereoDble.max(n_num_dbl_bonds_stereo);
    s.nLenLinearCTStereoCarb = s.nLenLinearCTStereoCarb.max(n_num_asym_carb_stereo);
    if let Some(info) = t_group_info {
        s.nLenIsotopicEndpoints = s.nLenIsotopicEndpoints.max(info.nNumIsotopicEndpoints);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn DeAllocateCS(
    heap: &mut SourceHeap,
    pCS: &mut CANON_STAT,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:491 DeAllocateCS
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int DeAllocateCS( CANON_STAT *pCS )
    {
    #define LOCAL_FREE( X) do{if(X){inchi_free( X); X=NULL;}}while(0)

        /* connection table */
        LOCAL_FREE( pCS->LinearCT );
        LOCAL_FREE( pCS->nCanonOrd );
        LOCAL_FREE( pCS->nSymmRank );
        LOCAL_FREE( pCS->nNum_H );
        LOCAL_FREE( pCS->nNum_H_fixed );
        LOCAL_FREE( pCS->nExchgIsoH );
        /* isotopic */
        LOCAL_FREE( pCS->LinearCTIsotopic );
        LOCAL_FREE( pCS->nSymmRankIsotopic );
        LOCAL_FREE( pCS->nCanonOrdIsotopic );
        /* isotopic tautomeric */
        LOCAL_FREE( pCS->LinearCTIsotopicTautomer );
        LOCAL_FREE( pCS->nCanonOrdIsotopicTaut );
        LOCAL_FREE( pCS->nSymmRankIsotopicTaut );
        /* stereo */
        LOCAL_FREE( pCS->LinearCTStereoDble );
        LOCAL_FREE( pCS->LinearCTStereoCarb );
        LOCAL_FREE( pCS->LinearCTStereoDbleInv );
        LOCAL_FREE( pCS->LinearCTStereoCarbInv );
        LOCAL_FREE( pCS->nCanonOrdStereo );
        LOCAL_FREE( pCS->nCanonOrdStereoInv );
        LOCAL_FREE( pCS->nCanonOrdStereoTaut );
        /* isotopic stereo */
        LOCAL_FREE( pCS->LinearCTIsotopicStereoDble );
        LOCAL_FREE( pCS->LinearCTIsotopicStereoCarb );
        LOCAL_FREE( pCS->LinearCTIsotopicStereoDbleInv );
        LOCAL_FREE( pCS->LinearCTIsotopicStereoCarbInv );
        LOCAL_FREE( pCS->bRankUsedForStereo );
        LOCAL_FREE( pCS->bAtomUsedForStereo );

        LOCAL_FREE( pCS->nCanonOrdIsotopicStereo );
        LOCAL_FREE( pCS->nCanonOrdIsotopicStereoInv );
        LOCAL_FREE( pCS->nCanonOrdIsotopicStereoTaut );
        /* tautomeric part of the connection table */
        LOCAL_FREE( pCS->LinearCTTautomer );
        LOCAL_FREE( pCS->nCanonOrdTaut );
        LOCAL_FREE( pCS->nSymmRankTaut );

        LOCAL_FREE( pCS->LinearCT2 );

        /* for establishing constitutional equivalence */
        LOCAL_FREE( pCS->nPrevAtomNumber );

        FreeNeighList( pCS->NeighList );
        pCS->NeighList = NULL;

        /* set zero lengths */
        pCS->nMaxLenLinearCTStereoDble = 0;
        pCS->nLenLinearCTStereoDble = 0;
        pCS->nMaxLenLinearCTStereoCarb = 0;
        pCS->nLenLinearCTStereoCarb = 0;
        pCS->nMaxLenLinearCTIsotopicStereoDble = 0;
        pCS->nLenLinearCTIsotopicStereoDble = 0;
        pCS->nMaxLenLinearCTIsotopicStereoCarb = 0;
        pCS->nLenLinearCTIsotopicStereoCarb = 0;
        pCS->nMaxLenLinearCTTautomer = 0;
        pCS->nLenLinearCTTautomer = 0;
        pCS->nMaxLenLinearCTIsotopic = 0;
        pCS->nLenLinearCTIsotopic = 0;
        pCS->nMaxLenLinearCTIsotopicTautomer = 0;
        pCS->nLenLinearCTIsotopicTautomer = 0;

        /* set canon numbering lengths to zero */
        pCS->nLenCanonOrd = 0;
        pCS->nLenCanonOrdIsotopic = 0;
        pCS->nLenCanonOrdIsotopicTaut = 0;
        pCS->nLenCanonOrdStereo = 0;
        pCS->nLenCanonOrdStereoTaut = 0;
        pCS->nLenCanonOrdIsotopicStereo = 0;
        pCS->nLenCanonOrdIsotopicStereoTaut = 0;
        pCS->nLenCanonOrdTaut = 0;

        return 0;

    #undef LOCAL_FREE
    }
    */
    // END INCHI C FUNCTION: DeAllocateCS

    macro_rules! local_free {
        ($field:ident) => {
            if !pCS.$field.is_null() {
                inchi_free(heap, pCS.$field)?;
                pCS.$field = SourceMutPointer::null();
            }
        };
    }

    local_free!(LinearCT);
    local_free!(nCanonOrd);
    local_free!(nSymmRank);
    local_free!(nNum_H);
    local_free!(nNum_H_fixed);
    local_free!(nExchgIsoH);
    local_free!(LinearCTIsotopic);
    local_free!(nSymmRankIsotopic);
    local_free!(nCanonOrdIsotopic);
    local_free!(LinearCTIsotopicTautomer);
    local_free!(nCanonOrdIsotopicTaut);
    local_free!(nSymmRankIsotopicTaut);
    local_free!(LinearCTStereoDble);
    local_free!(LinearCTStereoCarb);
    local_free!(LinearCTStereoDbleInv);
    local_free!(LinearCTStereoCarbInv);
    local_free!(nCanonOrdStereo);
    local_free!(nCanonOrdStereoInv);
    local_free!(nCanonOrdStereoTaut);
    local_free!(LinearCTIsotopicStereoDble);
    local_free!(LinearCTIsotopicStereoCarb);
    local_free!(LinearCTIsotopicStereoDbleInv);
    local_free!(LinearCTIsotopicStereoCarbInv);
    local_free!(bRankUsedForStereo);
    local_free!(bAtomUsedForStereo);
    local_free!(nCanonOrdIsotopicStereo);
    local_free!(nCanonOrdIsotopicStereoInv);
    local_free!(nCanonOrdIsotopicStereoTaut);
    local_free!(LinearCTTautomer);
    local_free!(nCanonOrdTaut);
    local_free!(nSymmRankTaut);
    local_free!(LinearCT2);
    local_free!(nPrevAtomNumber);

    FreeNeighList(heap, pCS.NeighList)?;
    pCS.NeighList = SourceMutPointer::null();

    pCS.nMaxLenLinearCTStereoDble = 0;
    pCS.nLenLinearCTStereoDble = 0;
    pCS.nMaxLenLinearCTStereoCarb = 0;
    pCS.nLenLinearCTStereoCarb = 0;
    pCS.nMaxLenLinearCTIsotopicStereoDble = 0;
    pCS.nLenLinearCTIsotopicStereoDble = 0;
    pCS.nMaxLenLinearCTIsotopicStereoCarb = 0;
    pCS.nLenLinearCTIsotopicStereoCarb = 0;
    pCS.nMaxLenLinearCTTautomer = 0;
    pCS.nLenLinearCTTautomer = 0;
    pCS.nMaxLenLinearCTIsotopic = 0;
    pCS.nLenLinearCTIsotopic = 0;
    pCS.nMaxLenLinearCTIsotopicTautomer = 0;
    pCS.nLenLinearCTIsotopicTautomer = 0;
    pCS.nLenCanonOrd = 0;
    pCS.nLenCanonOrdIsotopic = 0;
    pCS.nLenCanonOrdIsotopicTaut = 0;
    pCS.nLenCanonOrdStereo = 0;
    pCS.nLenCanonOrdStereoTaut = 0;
    pCS.nLenCanonOrdIsotopicStereo = 0;
    pCS.nLenCanonOrdIsotopicStereoTaut = 0;
    pCS.nLenCanonOrdTaut = 0;
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn AllocateCS(
    heap: &mut SourceHeap,
    pCS: &mut CANON_STAT,
    num_at: i32,
    num_at_tg: i32,
    nLenCT: i32,
    nLenCTAtOnly: i32,
    nLenLinearCTStereoDble: i32,
    nLenLinearCTIsotopicStereoDble: i32,
    nLenLinearCTStereoCarb: i32,
    nLenLinearCTIsotopicStereoCarb: i32,
    nLenLinearCTTautomer: i32,
    nLenLinearCTIsotopicTautomer: i32,
    nLenIsotopic: i32,
    nMode: crate::source_types::INCHI_MODE,
    pBCN: SourceMutPointer<crate::source_types::BCN>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:575 AllocateCS
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int AllocateCS( CANON_STAT *pCS,
                    int num_at,
                    int num_at_tg,
                    int nLenCT,
                    int nLenCTAtOnly,
                    int nLenLinearCTStereoDble,
                    int nLenLinearCTIsotopicStereoDble,
                    int nLenLinearCTStereoCarb,
                    int nLenLinearCTIsotopicStereoCarb,
                    int nLenLinearCTTautomer,
                    int nLenLinearCTIsotopicTautomer,
                    int nLenIsotopic,
                    INCHI_MODE nMode,
                    BCN *pBCN )
    {
    #define pCS_CALLOC( PTR,TYPE,LEN) (pCS->PTR=(TYPE*)inchi_calloc( (size_t)(LEN),sizeof(*pCS->PTR)))

        int num_err = 0;
        int num_t_groups = num_at_tg - num_at;

        pCS->nMode = nMode;

        /* connection table */
        if (( nMode & CMODE_CT ) && nLenCT > 0)
        {
            num_err += !pCS_CALLOC( LinearCT, AT_NUMB, nLenCT );
            pCS->nMaxLenLinearCT =
                pCS->nLenLinearCT = nLenCT;
            pCS->nLenLinearCTAtOnly = nLenCTAtOnly;
            num_err += !pCS_CALLOC( nCanonOrd, AT_RANK, num_at_tg );
            num_err += !pCS_CALLOC( nSymmRank, AT_RANK, num_at_tg );
            if (pBCN)
            {
                num_err += !pCS_CALLOC( nNum_H, S_CHAR, num_at );
                num_err += !pCS_CALLOC( nNum_H_fixed, S_CHAR, num_at );
                num_err += !pCS_CALLOC( nExchgIsoH, S_CHAR, num_at );
            }
        }

        /* isotopic */
        if (( nMode & CMODE_ISO ) && nLenIsotopic > 0)
        {
            num_err += !pCS_CALLOC( LinearCTIsotopic, AT_ISOTOPIC, nLenIsotopic );
            pCS->nMaxLenLinearCTIsotopic =
                pCS->nLenLinearCTIsotopic = nLenIsotopic;
        }

        /* isotopic tautomeric */
        if (( nMode & CMODE_ISO ) && CANON_MODE_TAUT == ( nMode & CANON_MODE_TAUT ))
        {
            if (nLenLinearCTIsotopicTautomer > 0)
            {
                num_err += !pCS_CALLOC( LinearCTIsotopicTautomer, AT_ISO_TGROUP, nLenLinearCTIsotopicTautomer );
                pCS->nMaxLenLinearCTIsotopicTautomer =
                    pCS->nLenLinearCTIsotopicTautomer = nLenLinearCTIsotopicTautomer;
            }
            if (num_t_groups > 0)
            {
                num_err += !pCS_CALLOC( nCanonOrdIsotopicTaut, AT_RANK, num_t_groups );
                num_err += !pCS_CALLOC( nSymmRankIsotopicTaut, AT_RANK, num_t_groups );
            }
        }
        /* isotopic atoms & t-groups */
        if (( nMode & CMODE_ISO ) /*&& nLenIsotopic > 0*/ ||
            (( nMode & CMODE_ISO ) && CANON_MODE_TAUT == ( nMode & CANON_MODE_TAUT ) && nLenLinearCTIsotopicTautomer > 0
            )) /* djb-rwth: addressing LLVM warning */
        {
            num_err += !pCS_CALLOC( nSymmRankIsotopic, AT_RANK, num_at_tg );
            num_err += !pCS_CALLOC( nCanonOrdIsotopic, AT_RANK, num_at_tg );
        }
        /* stereo */
        if (( nMode & CMODE_STEREO ) && nLenLinearCTStereoDble > 0)
        {
            num_err += !pCS_CALLOC( LinearCTStereoDble, AT_STEREO_DBLE, nLenLinearCTStereoDble );
            num_err += !pCS_CALLOC( LinearCTStereoDbleInv, AT_STEREO_DBLE, nLenLinearCTStereoDble );
            pCS->nLenLinearCTStereoDbleInv =
                pCS->nMaxLenLinearCTStereoDble =
                pCS->nLenLinearCTStereoDble = nLenLinearCTStereoDble;
        }
        if (( nMode & CMODE_STEREO ) && nLenLinearCTStereoCarb > 0)
        {
            num_err += !pCS_CALLOC( LinearCTStereoCarb, AT_STEREO_CARB, nLenLinearCTStereoCarb );
            num_err += !pCS_CALLOC( LinearCTStereoCarbInv, AT_STEREO_CARB, nLenLinearCTStereoCarb );
            pCS->nLenLinearCTStereoCarbInv =
                pCS->nMaxLenLinearCTStereoCarb =
                pCS->nLenLinearCTStereoCarb = nLenLinearCTStereoCarb;
        }
        if (( nMode & CMODE_STEREO ) && ( nLenLinearCTStereoDble > 0 || nLenLinearCTStereoCarb > 0 ))
        {
            num_err += !pCS_CALLOC( nCanonOrdStereo, AT_RANK, num_at_tg );
            num_err += !pCS_CALLOC( nCanonOrdStereoInv, AT_RANK, num_at_tg );
            if (( nMode & CMODE_TAUT ) && nLenLinearCTTautomer > 0 && num_t_groups > 0)
            {
                num_err += !pCS_CALLOC( nCanonOrdStereoTaut, AT_RANK, num_t_groups );
            }
        }
        /* isotopic stereo */
        if (( nMode & CMODE_ISO_STEREO ) && nLenLinearCTIsotopicStereoDble > 0)
        {
            num_err += !pCS_CALLOC( LinearCTIsotopicStereoDble, AT_STEREO_DBLE, nLenLinearCTIsotopicStereoDble );
            num_err += !pCS_CALLOC( LinearCTIsotopicStereoDbleInv, AT_STEREO_DBLE, nLenLinearCTIsotopicStereoDble );
            pCS->nLenLinearCTIsotopicStereoDbleInv =
                pCS->nMaxLenLinearCTIsotopicStereoDble =
                pCS->nLenLinearCTIsotopicStereoDble = nLenLinearCTIsotopicStereoDble;
        }
        if (( nMode & CMODE_ISO_STEREO ) && nLenLinearCTIsotopicStereoCarb > 0)
        {
            num_err += !pCS_CALLOC( LinearCTIsotopicStereoCarb, AT_STEREO_CARB, nLenLinearCTIsotopicStereoCarb );
            num_err += !pCS_CALLOC( LinearCTIsotopicStereoCarbInv, AT_STEREO_CARB, nLenLinearCTIsotopicStereoCarb );
            pCS->nLenLinearCTIsotopicStereoCarbInv =
                pCS->nMaxLenLinearCTIsotopicStereoCarb =
                pCS->nLenLinearCTIsotopicStereoCarb = nLenLinearCTIsotopicStereoCarb;
        }
        if (( nMode & CMODE_ISO_STEREO ) && ( nLenLinearCTIsotopicStereoDble > 0 || nLenLinearCTIsotopicStereoCarb > 0 ))
        {
            num_err += !pCS_CALLOC( nCanonOrdIsotopicStereo, AT_RANK, num_at_tg );
            num_err += !pCS_CALLOC( nCanonOrdIsotopicStereoInv, AT_RANK, num_at_tg );
            if (( nMode & CMODE_TAUT ) && nLenLinearCTTautomer > 0 && num_t_groups > 0)
            {
                num_err += !pCS_CALLOC( nCanonOrdIsotopicStereoTaut, AT_RANK, num_t_groups );
            }
        }
        if ((( nMode & CMODE_STEREO ) && ( nLenLinearCTStereoDble > 0 || nLenLinearCTStereoCarb > 0 )) ||
            (( nMode & CMODE_ISO_STEREO ) && ( nLenLinearCTIsotopicStereoDble > 0 || nLenLinearCTIsotopicStereoCarb > 0 ))) /* djb-rwth: addressing LLVM warning */
        {
            num_err += !pCS_CALLOC( bRankUsedForStereo, S_CHAR, num_at );
            num_err += !pCS_CALLOC( bAtomUsedForStereo, S_CHAR, num_at );
        }
        /* tautomeric part of the connection table */
        if (( nMode & CMODE_CT ) && ( nMode & CMODE_TAUT ) && nLenLinearCTTautomer > 0)
        {
            num_err += !pCS_CALLOC( LinearCTTautomer, AT_TAUTOMER, nLenLinearCTTautomer );
            pCS->nMaxLenLinearCTTautomer =
                pCS->nLenLinearCTTautomer = nLenLinearCTTautomer;
            if (num_t_groups > 0)
            {
                num_err += !pCS_CALLOC( nCanonOrdTaut, AT_RANK, num_t_groups );
                num_err += !pCS_CALLOC( nSymmRankTaut, AT_RANK, num_t_groups );
            }
        }

        if (nMode & CMODE_CT)
        {
            num_err += !pCS_CALLOC( LinearCT2, AT_NUMB, nLenCT );
        }

            /* for establishing constitutional equivalence */
        num_err += !pCS_CALLOC( nPrevAtomNumber, AT_RANK, num_at_tg );

        /* set canon numbering lengths to zero */
        pCS->nLenCanonOrd = 0;
        pCS->nLenCanonOrdIsotopic = 0;
        pCS->nLenCanonOrdIsotopicTaut = 0;
        pCS->nLenCanonOrdStereo = 0;
        pCS->nLenCanonOrdStereoTaut = 0;
        pCS->nLenCanonOrdIsotopicStereo = 0;
        pCS->nLenCanonOrdIsotopicStereoTaut = 0;
        pCS->nLenCanonOrdTaut = 0;


        if (num_err)
        {
            DeAllocateCS( pCS );
            return CT_OUT_OF_RAM;  /*  <BRKPT> */
        }
        return 0;

    #undef pCS_CALLOC
    }
    */
    // END INCHI C FUNCTION: AllocateCS

    let mut num_err = 0_i32;
    let num_t_groups = num_at_tg.wrapping_sub(num_at);
    pCS.nMode = nMode;

    macro_rules! calloc_field {
        ($field:ident, $type:ty, $length:expr) => {{
            let source_count = (i64::from($length)) as u64;
            match inchi_calloc::<$type>(heap, source_count, size_of::<$type>() as u64) {
                Ok(pointer) => pCS.$field = pointer,
                Err(_) => {
                    pCS.$field = SourceMutPointer::null();
                    num_err = num_err.wrapping_add(1);
                }
            }
        }};
    }

    if nMode & u64::from(crate::source_types::CMODE_CT) != 0 && nLenCT > 0 {
        calloc_field!(LinearCT, AT_NUMB, nLenCT);
        pCS.nMaxLenLinearCT = nLenCT;
        pCS.nLenLinearCT = nLenCT;
        pCS.nLenLinearCTAtOnly = nLenCTAtOnly;
        calloc_field!(nCanonOrd, AT_RANK, num_at_tg);
        calloc_field!(nSymmRank, AT_RANK, num_at_tg);
        if !pBCN.is_null() {
            calloc_field!(nNum_H, i8, num_at);
            calloc_field!(nNum_H_fixed, i8, num_at);
            calloc_field!(nExchgIsoH, i8, num_at);
        }
    }
    if nMode & u64::from(CMODE_ISO) != 0 && nLenIsotopic > 0 {
        calloc_field!(LinearCTIsotopic, AT_ISOTOPIC, nLenIsotopic);
        pCS.nMaxLenLinearCTIsotopic = nLenIsotopic;
        pCS.nLenLinearCTIsotopic = nLenIsotopic;
    }
    if nMode & u64::from(CMODE_ISO) != 0
        && u64::from(crate::source_types::CANON_MODE_TAUT)
            == nMode & u64::from(crate::source_types::CANON_MODE_TAUT)
    {
        if nLenLinearCTIsotopicTautomer > 0 {
            calloc_field!(
                LinearCTIsotopicTautomer,
                AT_ISO_TGROUP,
                nLenLinearCTIsotopicTautomer
            );
            pCS.nMaxLenLinearCTIsotopicTautomer = nLenLinearCTIsotopicTautomer;
            pCS.nLenLinearCTIsotopicTautomer = nLenLinearCTIsotopicTautomer;
        }
        if num_t_groups > 0 {
            calloc_field!(nCanonOrdIsotopicTaut, AT_RANK, num_t_groups);
            calloc_field!(nSymmRankIsotopicTaut, AT_RANK, num_t_groups);
        }
    }
    if nMode & u64::from(CMODE_ISO) != 0 {
        calloc_field!(nSymmRankIsotopic, AT_RANK, num_at_tg);
        calloc_field!(nCanonOrdIsotopic, AT_RANK, num_at_tg);
    }
    if nMode & u64::from(CMODE_STEREO) != 0 && nLenLinearCTStereoDble > 0 {
        calloc_field!(
            LinearCTStereoDble,
            crate::source_types::AT_STEREO_DBLE,
            nLenLinearCTStereoDble
        );
        calloc_field!(
            LinearCTStereoDbleInv,
            crate::source_types::AT_STEREO_DBLE,
            nLenLinearCTStereoDble
        );
        pCS.nLenLinearCTStereoDbleInv = nLenLinearCTStereoDble;
        pCS.nMaxLenLinearCTStereoDble = nLenLinearCTStereoDble;
        pCS.nLenLinearCTStereoDble = nLenLinearCTStereoDble;
    }
    if nMode & u64::from(CMODE_STEREO) != 0 && nLenLinearCTStereoCarb > 0 {
        calloc_field!(
            LinearCTStereoCarb,
            crate::source_types::AT_STEREO_CARB,
            nLenLinearCTStereoCarb
        );
        calloc_field!(
            LinearCTStereoCarbInv,
            crate::source_types::AT_STEREO_CARB,
            nLenLinearCTStereoCarb
        );
        pCS.nLenLinearCTStereoCarbInv = nLenLinearCTStereoCarb;
        pCS.nMaxLenLinearCTStereoCarb = nLenLinearCTStereoCarb;
        pCS.nLenLinearCTStereoCarb = nLenLinearCTStereoCarb;
    }
    if nMode & u64::from(CMODE_STEREO) != 0
        && (nLenLinearCTStereoDble > 0 || nLenLinearCTStereoCarb > 0)
    {
        calloc_field!(nCanonOrdStereo, AT_RANK, num_at_tg);
        calloc_field!(nCanonOrdStereoInv, AT_RANK, num_at_tg);
        if nMode & u64::from(crate::source_types::CMODE_TAUT) != 0
            && nLenLinearCTTautomer > 0
            && num_t_groups > 0
        {
            calloc_field!(nCanonOrdStereoTaut, AT_RANK, num_t_groups);
        }
    }
    if nMode & u64::from(CMODE_ISO_STEREO) != 0 && nLenLinearCTIsotopicStereoDble > 0 {
        calloc_field!(
            LinearCTIsotopicStereoDble,
            crate::source_types::AT_STEREO_DBLE,
            nLenLinearCTIsotopicStereoDble
        );
        calloc_field!(
            LinearCTIsotopicStereoDbleInv,
            crate::source_types::AT_STEREO_DBLE,
            nLenLinearCTIsotopicStereoDble
        );
        pCS.nLenLinearCTIsotopicStereoDbleInv = nLenLinearCTIsotopicStereoDble;
        pCS.nMaxLenLinearCTIsotopicStereoDble = nLenLinearCTIsotopicStereoDble;
        pCS.nLenLinearCTIsotopicStereoDble = nLenLinearCTIsotopicStereoDble;
    }
    if nMode & u64::from(CMODE_ISO_STEREO) != 0 && nLenLinearCTIsotopicStereoCarb > 0 {
        calloc_field!(
            LinearCTIsotopicStereoCarb,
            crate::source_types::AT_STEREO_CARB,
            nLenLinearCTIsotopicStereoCarb
        );
        calloc_field!(
            LinearCTIsotopicStereoCarbInv,
            crate::source_types::AT_STEREO_CARB,
            nLenLinearCTIsotopicStereoCarb
        );
        pCS.nLenLinearCTIsotopicStereoCarbInv = nLenLinearCTIsotopicStereoCarb;
        pCS.nMaxLenLinearCTIsotopicStereoCarb = nLenLinearCTIsotopicStereoCarb;
        pCS.nLenLinearCTIsotopicStereoCarb = nLenLinearCTIsotopicStereoCarb;
    }
    if nMode & u64::from(CMODE_ISO_STEREO) != 0
        && (nLenLinearCTIsotopicStereoDble > 0 || nLenLinearCTIsotopicStereoCarb > 0)
    {
        calloc_field!(nCanonOrdIsotopicStereo, AT_RANK, num_at_tg);
        calloc_field!(nCanonOrdIsotopicStereoInv, AT_RANK, num_at_tg);
        if nMode & u64::from(crate::source_types::CMODE_TAUT) != 0
            && nLenLinearCTTautomer > 0
            && num_t_groups > 0
        {
            calloc_field!(nCanonOrdIsotopicStereoTaut, AT_RANK, num_t_groups);
        }
    }
    if (nMode & u64::from(CMODE_STEREO) != 0
        && (nLenLinearCTStereoDble > 0 || nLenLinearCTStereoCarb > 0))
        || (nMode & u64::from(CMODE_ISO_STEREO) != 0
            && (nLenLinearCTIsotopicStereoDble > 0 || nLenLinearCTIsotopicStereoCarb > 0))
    {
        calloc_field!(bRankUsedForStereo, i8, num_at);
        calloc_field!(bAtomUsedForStereo, i8, num_at);
    }
    if nMode & u64::from(crate::source_types::CMODE_CT) != 0
        && nMode & u64::from(crate::source_types::CMODE_TAUT) != 0
        && nLenLinearCTTautomer > 0
    {
        calloc_field!(LinearCTTautomer, AT_TAUTOMER, nLenLinearCTTautomer);
        pCS.nMaxLenLinearCTTautomer = nLenLinearCTTautomer;
        pCS.nLenLinearCTTautomer = nLenLinearCTTautomer;
        if num_t_groups > 0 {
            calloc_field!(nCanonOrdTaut, AT_RANK, num_t_groups);
            calloc_field!(nSymmRankTaut, AT_RANK, num_t_groups);
        }
    }
    if nMode & u64::from(crate::source_types::CMODE_CT) != 0 {
        calloc_field!(LinearCT2, AT_NUMB, nLenCT);
    }
    calloc_field!(nPrevAtomNumber, AT_RANK, num_at_tg);

    pCS.nLenCanonOrd = 0;
    pCS.nLenCanonOrdIsotopic = 0;
    pCS.nLenCanonOrdIsotopicTaut = 0;
    pCS.nLenCanonOrdStereo = 0;
    pCS.nLenCanonOrdStereoTaut = 0;
    pCS.nLenCanonOrdIsotopicStereo = 0;
    pCS.nLenCanonOrdIsotopicStereoTaut = 0;
    pCS.nLenCanonOrdTaut = 0;

    if num_err != 0 {
        DeAllocateCS(heap, pCS)?;
        return Ok(CT_OUT_OF_RAM);
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn FillTautLinearCT2(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    num_atoms: i32,
    num_at_tg: i32,
    bIsoTaut: i32,
    nRank: SourceConstPointer<AT_RANK>,
    nAtomNumber: SourceConstPointer<AT_RANK>,
    nSymmRank: SourceConstPointer<AT_RANK>,
    nRankIso: SourceConstPointer<AT_RANK>,
    nAtomNumberIso: SourceConstPointer<AT_RANK>,
    nSymmRankIso: SourceConstPointer<AT_RANK>,
    LinearCTTautomer: SourceMutPointer<AT_TAUTOMER>,
    nMaxLenLinearCTTautomer: i32,
    pnLenLinearCTTautomer: &mut i32,
    LinearCTIsotopicTautomer: SourceMutPointer<AT_ISO_TGROUP>,
    nMaxLenLinearCTIsotopicTautomer: i32,
    pnLenLinearCTIsotopicTautomer: &mut i32,
    t_group_info: Option<&mut T_GROUP_INFO>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:858 FillTautLinearCT2
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // END INCHI C FUNCTION: FillTautLinearCT2

    let mut len = 0_i32;
    let Some(info) = t_group_info else {
        return Ok(0);
    };
    if num_at_tg <= num_atoms || info.num_t_groups == 0 {
        return Ok(0);
    }
    let group_count =
        usize::try_from(info.num_t_groups).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom_count = usize::try_from(num_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let total_count =
        usize::try_from(num_at_tg).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let offset = num_atoms as AT_RANK;
    for (group_index, pseudo_index) in (atom_count..total_count).enumerate() {
        let group_number = heap
            .slice(nAtomNumber)?
            .get(pseudo_index)
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .wrapping_sub(offset);
        let symmetry = heap
            .slice(nSymmRank)?
            .get(pseudo_index)
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .wrapping_sub(offset);
        let work = heap.slice_mut(info.tGroupNumber)?;
        *work
            .get_mut(group_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = group_number;
        *work
            .get_mut(TGSO_SYMM_RANK as usize * group_count + group_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)? = symmetry;
        if bIsoTaut != 0 {
            let iso_group = heap
                .slice(nAtomNumberIso)?
                .get(pseudo_index)
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .wrapping_sub(offset);
            let iso_symmetry = heap
                .slice(nSymmRankIso)?
                .get(pseudo_index)
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .wrapping_sub(offset);
            let work = heap.slice_mut(info.tGroupNumber)?;
            *work
                .get_mut(TGSO_SYMM_IORDER as usize * group_count + group_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = iso_group;
            *work
                .get_mut(TGSO_SYMM_IRANK as usize * group_count + group_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = iso_symmetry;
        }
    }

    pCG.m_pn_RankForSort = nRank;
    for group_index in 0..group_count {
        let group = heap
            .slice(info.t_group.as_const())?
            .get(group_index)
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let start = usize::from(group.nFirstEndpointAtNoPos);
        let count = usize::from(group.nNumEndpoints);
        let mut endpoints = heap
            .slice(info.nEndpointAtomNumber.as_const())?
            .get(start..start + count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec();
        inchi_qsort(
            bytemuck::cast_slice_mut(&mut endpoints),
            count,
            size_of::<AT_NUMB>(),
            &mut |left, right| {
                CompRank(
                    heap,
                    AT_NUMB::from_ne_bytes([left[0], left[1]]),
                    AT_NUMB::from_ne_bytes([right[0], right[1]]),
                    pCG,
                )
            },
        )?;
        heap.slice_mut(info.nEndpointAtomNumber)?[start..start + count].copy_from_slice(&endpoints);
    }

    let mut max_len = 0_i32;
    if nMaxLenLinearCTTautomer != 0 {
        max_len = (T_GROUP_HDR_LEN as i32)
            .wrapping_mul(info.num_t_groups)
            .wrapping_add(info.nNumEndpoints)
            .wrapping_add(1);
        if max_len > nMaxLenLinearCTTautomer {
            return Ok(CT_OVERFLOW);
        }
    }
    for canonical_index in 0..group_count {
        let group_number = usize::from(
            *heap
                .slice(info.tGroupNumber.as_const())?
                .get(canonical_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let group = heap
            .slice(info.t_group.as_const())?
            .get(group_number)
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let required = len
            .wrapping_add(T_GROUP_HDR_LEN as i32)
            .wrapping_add(i32::from(group.nNumEndpoints));
        if required >= max_len {
            return Ok(CT_OVERFLOW);
        }
        let output = heap.slice_mut(LinearCTTautomer)?;
        output[usize::try_from(len).unwrap()] = group.nNumEndpoints;
        len = len.wrapping_add(1);
        for value in &group.num[..T_NUM_NO_ISOTOPIC as usize] {
            output[usize::try_from(len).unwrap()] = *value;
            len = len.wrapping_add(1);
        }
        let start = usize::from(group.nFirstEndpointAtNoPos);
        for endpoint_index in 0..usize::from(group.nNumEndpoints) {
            let atom_number = usize::from(
                *heap
                    .slice(info.nEndpointAtomNumber.as_const())?
                    .get(start + endpoint_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            let rank = *heap
                .slice(nRank)?
                .get(atom_number)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            heap.slice_mut(LinearCTTautomer)?[usize::try_from(len).unwrap()] = rank;
            len = len.wrapping_add(1);
        }
    }
    if nMaxLenLinearCTTautomer != 0 {
        heap.slice_mut(LinearCTTautomer)?[usize::try_from(len).unwrap()] = 0;
        len = len.wrapping_add(1);
        if len != max_len {
            len = len.wrapping_neg();
        } else if *pnLenLinearCTTautomer != 0 && *pnLenLinearCTTautomer != len {
            return Ok(CT_LEN_MISMATCH);
        } else {
            *pnLenLinearCTTautomer = len;
        }
    } else {
        *pnLenLinearCTTautomer = 0;
    }

    let mut len_iso = 0_i32;
    if nMaxLenLinearCTIsotopicTautomer != 0 && info.nNumIsotopicEndpoints == 0 {
        for canonical_index in 0..group_count {
            let group_number = usize::from(
                *heap
                    .slice(info.tGroupNumber.as_const())?
                    .get(TGSO_SYMM_IORDER as usize * group_count + canonical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            let group = heap
                .slice(info.t_group.as_const())?
                .get(group_number)
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if group.iWeight == 0 {
                continue;
            }
            if len_iso >= nMaxLenLinearCTIsotopicTautomer {
                return Ok(CT_OVERFLOW);
            }
            let record = heap
                .slice_mut(LinearCTIsotopicTautomer)?
                .get_mut(usize::try_from(len_iso).unwrap())
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            for index in T_NUM_NO_ISOTOPIC as usize..group.num.len() {
                record.num[index - T_NUM_NO_ISOTOPIC as usize] = group.num[index];
            }
            record.tgroup_num = (canonical_index as AT_NUMB).wrapping_add(1);
            len_iso = len_iso.wrapping_add(1);
        }
    }
    if nMaxLenLinearCTIsotopicTautomer != 0 {
        if *pnLenLinearCTIsotopicTautomer != 0 && *pnLenLinearCTIsotopicTautomer != len_iso {
            return Ok(CT_LEN_MISMATCH);
        }
        *pnLenLinearCTIsotopicTautomer = len_iso;
    } else {
        *pnLenLinearCTIsotopicTautomer = 0;
    }
    let _ = nRankIso;
    Ok(len)
}

fn update_full_linear_ct_value(
    heap: &mut SourceHeap,
    linear_ct: SourceMutPointer<AT_NUMB>,
    length: &mut i32,
    maximum_length: i32,
    value: AT_NUMB,
    compare: &mut i32,
) -> Result<Option<i32>, SourceHeapError> {
    if *length >= maximum_length {
        return Ok(Some(CT_OVERFLOW));
    }
    let index = usize::try_from(*length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if *compare != 0 {
        let previous = *heap
            .slice(linear_ct.as_const())?
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if value > previous {
            return Ok(Some(1));
        }
        *compare = i32::from(value == previous);
    }
    *heap
        .slice_mut(linear_ct)?
        .get_mut(index)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = value;
    *length = length.wrapping_add(1);
    Ok(None)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn UpdateFullLinearCT(
    heap: &mut SourceHeap,
    num_atoms: i32,
    num_at_tg: i32,
    at: SourceConstPointer<sp_ATOM>,
    nRank: SourceMutPointer<AT_RANK>,
    nAtomNumber: SourceMutPointer<AT_RANK>,
    pCS: &mut CANON_STAT,
    pCG: &mut CANON_GLOBALS,
    bFirstTime: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:1050 UpdateFullLinearCT
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap accesses add overhead.
    /*
    int UpdateFullLinearCT( int num_atoms,
                            int num_at_tg,
                            sp_ATOM* at,
                            AT_RANK *nRank,
                            AT_RANK *nAtomNumber,
                            CANON_STAT* pCS,
                            CANON_GLOBALS *pCG,
                            int bFirstTime )
    {
        /* at[i].init_rank = initial ranks before canonizing */
        /* nRank[i]  = new ordering number for atoms: nRank=1,2,.. */
        /* nAtomNumber[r] = orig. atom number= 0,1,...  for r = nRank-1  */
        /* nRank[nAtomNumber[r]] = r; r = 0,1,... */
        /* nAtomNumber[nRank[i]-1] = i; */
    
        AT_NUMB nNeighborNumber[MAXVAL];
        int  i, j, k, num_neigh, rank, bCompare; /*, nRetVal; */
    
        T_GROUP_INFO *t_group_info = NULL;
        T_GROUP      *t_group = NULL;
        AT_NUMB      *nEndpointAtomNumber = NULL;
    
        int  nCTLen = 0, nCTLenAtOnly = 0;
        
    
        AT_NUMB         r_neigh;
        AT_NUMB        *LinearCT = pCS->LinearCT;
    
        /* the following parts of the "name" should be compared */
        /* after the connection table comparison is done */
        /* to avoid wrong difference sign. So, these parts */
        /* go to a separate buffers. */
        /* -- currently not used at all at all -- */
    
    #if CT_ATOMID != CT_ATOMID_DONTINCLUDE
        AT_NUMB          r0_at_type;
    #endif
    
        num_neigh = 0; /* Moved from above 2024-09-01 DT; djb-rwth: num_neigh initialisation added */
    
        bCompare = bFirstTime ? 0 : 1;
    
        if (num_at_tg > num_atoms)
        {
            t_group_info = pCS->t_group_info;
            t_group = t_group_info->t_group;
        }
        else
        {
            t_group_info = NULL;
            t_group = NULL;
        }
    
        /**********************************************************************/
        /*                                                                    */
        /*    CYCLE 1: FILL OUT  CONNECTION TABLE(S) FOR ALL ATOMS            */
        /*      ** NOT INCLUDING ISOTOPIC ATOMS AND 1H, 2H(D), 3H(T) **       */
        /*                                                                    */
        /* rank = nRank[nAtomNumber[rank-1]] -- proposed atoms canon. numbers */
        /**********************************************************************/
        for (rank = 1; rank <= num_atoms; rank++)
        {
            i = (int) nAtomNumber[rank - 1];  /* current atom */
    
    #if ( CT_ATOMID == CT_ATOMID_IS_CURRANK )
            r0_at_type = (AT_NUMB) rank; /* current Rank */
    #else
    #if ( CT_ATOMID == CT_ATOMID_IS_INITRANK )
            r0_at_type = (AT_NUMB) at[i].init_rank; /* chemical + neighborhood ID */
    #else
    #if ( CT_ATOMID == CT_ATOMID_DONTINCLUDE )
    #else
    #error Undefined or wrong definition of CT_ATOMID
    #endif
    #endif
    #endif
    
            /* add atom to the CT */
    #if ( CT_ATOMID != CT_ATOMID_DONTINCLUDE )
            if (CHECK_OVERFLOW( nCTLen, pCS->nMaxLenLinearCT ))
            {
                return CT_OVERFLOW;  /*  <BRKPT> */
            }
            COMPARE_WITH_CT( LinearCT, nCTLen, r0_at_type, bCompare );
    #endif
    
            /*******************************************************
                 add neighbors and (if required) bonds to CT
            ********************************************************/
    
            /* sort neighbors */
            num_neigh = at[i].valence;
            for (k = 0; k < num_neigh; k++)
            {
                nNeighborNumber[k] = (AT_NUMB) k;
            }
            pCG->m_pNeighborsForSort = at[i].neighbor;
            pCG->m_pn_RankForSort = nRank;
            insertions_sort( pCG, nNeighborNumber, (size_t) num_neigh, sizeof( nNeighborNumber[0] ), CompNeighborsAT_NUMBER );
    
            for (k = 0; k < num_neigh; k++)
            {
                /* rank = (new current atom Rank) */
                if ((int) ( r_neigh = (AT_NUMB) nRank[(int) at[i].neighbor[(int) nNeighborNumber[k]]] )
                                                                  CT_NEIGH_SMALLER_THAN rank)
                {
                    if (CHECK_OVERFLOW( nCTLen, pCS->nMaxLenLinearCT ))
                    {
                        return CT_OVERFLOW;  /*  <BRKPT> */
                    }
                    COMPARE_WITH_CT( LinearCT, nCTLen, r_neigh, bCompare );
                }
            }
    
            /* add CT row delimiter */
        } /* end of cycle over all atoms. */
    
        nCTLenAtOnly = nCTLen;
    
        /**************************************************************
    
                    Tautomeric groups 07-22-2002
    
        ***************************************************************/
    
        for (rank = num_atoms + 1; rank <= num_at_tg; rank++)
        {
            j = (int) nAtomNumber[rank - 1];  /* current "atom" */
            i = j - num_atoms;             /* current t-group */
    
    #if ( CT_ATOMID == CT_ATOMID_IS_CURRANK )
            r0_at_type = (AT_NUMB) rank; /* current Rank */
    #else
    #if ( CT_ATOMID == CT_ATOMID_IS_INITRANK )
            r0_at_type = (AT_NUMB) rank; /* current Rank or  (AT_NUMB)at[i].init_rank; ==> chemical + neighborhood ID */
    #else
    #if ( CT_ATOMID == CT_ATOMID_DONTINCLUDE )
    #else
    #error Undefined or wrong definition of CT_ATOMID
    #endif
    #endif
    #endif
    
            /* add atom to the CT */
    #if ( CT_ATOMID != CT_ATOMID_DONTINCLUDE )
            if (CHECK_OVERFLOW( nCTLen, pCS->nMaxLenLinearCT ))
                return CT_OVERFLOW;  /*  <BRKPT> */
            COMPARE_WITH_CT( LinearCT, nCTLen, r0_at_type, bCompare );
    #endif
    
            /*******************************************************
                  add neighbors and (if required) bonds to CT
            ********************************************************/
    
            /* sort endpoints */
            if (t_group_info) 
                nEndpointAtomNumber = t_group_info->nEndpointAtomNumber + (int) t_group[i].nFirstEndpointAtNoPos; /* djb-rwth: fixing a NULL pointer dereference */
            pCG->m_pn_RankForSort = nRank;
            if (t_group + i) /* djb-rwth: ignoring GCC warning */
                num_neigh = (int)t_group[i].nNumEndpoints;
            insertions_sort( pCG, nEndpointAtomNumber, (size_t) num_neigh, sizeof( nEndpointAtomNumber[0] ), CompRank );
    
            for (k = 0; k < num_neigh; k++)
            {
                /* rank = (new current atom Rank) */
                if (nEndpointAtomNumber + k) /* djb-rwth: fixing a NULL pointer dereference; ignoring GCC warning */
                {
                    if ((int)(r_neigh = (AT_NUMB)nRank[(int)nEndpointAtomNumber[k]]) CT_NEIGH_SMALLER_THAN rank) 
                    {
                        if (CHECK_OVERFLOW(nCTLen, pCS->nMaxLenLinearCT))
                        {
                            return CT_OVERFLOW;  /*  <BRKPT> */
                        }
                        COMPARE_WITH_CT(LinearCT, nCTLen, r_neigh, bCompare);
                    }
                }
            }
        } /* end of cycle over all tautomeric groups. */
    
        /* compare bonds types */
        /* compare elements */
    
        if (LinearCT)
        {
    
            if (pCS->nLenLinearCT)
            {
                if (pCS->nLenLinearCT != nCTLen)
                {
                    return CT_LEN_MISMATCH;  /*  <BRKPT> */
                }
            }
            else
            {
                pCS->nLenLinearCT = nCTLen;
            }
    
            if (pCS->nLenLinearCT)
            {
                if (pCS->nLenLinearCTAtOnly != nCTLenAtOnly)
                {
                    return CT_LEN_MISMATCH;  /*  <BRKPT> */
                }
            }
            else
            {
                pCS->nLenLinearCTAtOnly = nCTLenAtOnly;
            }
        }
    
        /* Return: 0=> identical CT; -1=> new CT is smaller than the previous one */
        return ( bCompare - 1 );
    }
    
    
    */
    // END INCHI C FUNCTION: UpdateFullLinearCT
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: UpdateFullLinearCT
    // INCHI✔️❌: #define CT_ATOMID CT_ATOMID_IS_CURRANK
    // INCHI✔️❌: #define CT_NEIGH_SMALLER_THAN <
    // INCHI✔️❌: #define CT_NEIGH_INCREASE
    // INCHI✔️❌: #define CT_GREATER_THAN >
    // INCHI✔️❌: #define CHECK_OVERFLOW(Len, Maxlen) ((Len) >= (Maxlen))
    // END INCHI ACTIVE MACRO CONFIGURATION: UpdateFullLinearCT

    let atom_count =
        usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let total_count =
        usize::try_from(num_at_tg).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if total_count < atom_count {
        return Err(SourceHeapError::SourceIntegerOverflow);
    }
    let mut ct_length = 0_i32;
    let mut compare = i32::from(bFirstTime == 0);
    let linear_ct = pCS.LinearCT;
    let mut neighbor_numbers = [0_u16; MAXVAL as usize];
    let t_group_info = if total_count > atom_count {
        Some(
            heap.slice(pCS.t_group_info.as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
    } else {
        None
    };

    for rank in 1..=atom_count {
        let atom_index = usize::from(
            *heap
                .slice(nAtomNumber.as_const())?
                .get(rank - 1)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let atom = heap
            .slice(at)?
            .get(atom_index)
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if let Some(result) = update_full_linear_ct_value(
            heap,
            linear_ct,
            &mut ct_length,
            pCS.nMaxLenLinearCT,
            rank as AT_NUMB,
            &mut compare,
        )? {
            return Ok(result);
        }

        let neighbor_count = usize::try_from(i32::from(atom.valence).max(0))
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if neighbor_count > neighbor_numbers.len() {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        for (index, value) in neighbor_numbers[..neighbor_count].iter_mut().enumerate() {
            *value = index as AT_NUMB;
        }
        pCG.m_pn_RankForSort = nRank.as_const();
        if neighbor_count > 1 {
            let ranks = heap.slice(nRank.as_const())?;
            let bytes =
                bytemuck::cast_slice_mut::<AT_NUMB, u8>(&mut neighbor_numbers[..neighbor_count]);
            insertions_sort(
                bytes,
                neighbor_count,
                size_of::<AT_NUMB>(),
                &mut |first, second| {
                    let first = AT_NUMB::from_ne_bytes(
                        first
                            .try_into()
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    );
                    let second = AT_NUMB::from_ne_bytes(
                        second
                            .try_into()
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    );
                    CompNeighborsAT_NUMBER(
                        first,
                        second,
                        CompNeighborsATNumberContext::Slices {
                            neighbors: &atom.neighbor,
                            ranks,
                        },
                    )
                },
            )?;
        }
        for &neighbor_order in &neighbor_numbers[..neighbor_count] {
            let neighbor = atom.neighbor[usize::from(neighbor_order)];
            let neighbor_rank = *heap
                .slice(nRank.as_const())?
                .get(usize::from(neighbor))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if usize::from(neighbor_rank) < rank
                && let Some(result) = update_full_linear_ct_value(
                    heap,
                    linear_ct,
                    &mut ct_length,
                    pCS.nMaxLenLinearCT,
                    neighbor_rank,
                    &mut compare,
                )?
            {
                return Ok(result);
            }
        }
    }
    let atom_only_length = ct_length;

    for rank in atom_count + 1..=total_count {
        let pseudo_atom = usize::from(
            *heap
                .slice(nAtomNumber.as_const())?
                .get(rank - 1)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let group_index = pseudo_atom
            .checked_sub(atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if let Some(result) = update_full_linear_ct_value(
            heap,
            linear_ct,
            &mut ct_length,
            pCS.nMaxLenLinearCT,
            rank as AT_NUMB,
            &mut compare,
        )? {
            return Ok(result);
        }

        let info = t_group_info.as_ref().ok_or(SourceHeapError::NullPointer)?;
        let group = heap
            .slice(info.t_group.as_const())?
            .get(group_index)
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let endpoint_count = usize::from(group.nNumEndpoints);
        let endpoint_pointer = info
            .nEndpointAtomNumber
            .offset(i64::from(group.nFirstEndpointAtNoPos))?;
        pCG.m_pn_RankForSort = nRank.as_const();
        if endpoint_count > 1 {
            heap.with_slice_mut_and_heap(endpoint_pointer, |endpoints, heap| {
                let endpoints = endpoints
                    .get_mut(..endpoint_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let bytes = bytemuck::cast_slice_mut::<AT_NUMB, u8>(endpoints);
                insertions_sort(
                    bytes,
                    endpoint_count,
                    size_of::<AT_NUMB>(),
                    &mut |first, second| {
                        let first = AT_NUMB::from_ne_bytes(
                            first
                                .try_into()
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        );
                        let second = AT_NUMB::from_ne_bytes(
                            second
                                .try_into()
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        );
                        CompRank(heap, first, second, pCG)
                    },
                )
                .map(|_| ())
            })?;
        }
        for index in 0..endpoint_count {
            let endpoint = heap.slice(endpoint_pointer.as_const())?[index];
            let endpoint_rank = *heap
                .slice(nRank.as_const())?
                .get(usize::from(endpoint))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if usize::from(endpoint_rank) < rank
                && let Some(result) = update_full_linear_ct_value(
                    heap,
                    linear_ct,
                    &mut ct_length,
                    pCS.nMaxLenLinearCT,
                    endpoint_rank,
                    &mut compare,
                )?
            {
                return Ok(result);
            }
        }
    }

    if !linear_ct.is_null() {
        if pCS.nLenLinearCT != 0 {
            if pCS.nLenLinearCT != ct_length {
                return Ok(CT_LEN_MISMATCH);
            }
        } else {
            pCS.nLenLinearCT = ct_length;
        }
        if pCS.nLenLinearCT != 0 {
            if pCS.nLenLinearCTAtOnly != atom_only_length {
                return Ok(CT_LEN_MISMATCH);
            }
        } else {
            pCS.nLenLinearCTAtOnly = atom_only_length;
        }
    }
    Ok(compare - 1)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn FixCanonEquivalenceInfo(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    num_at_tg: i32,
    nSymmRank: SourceMutPointer<AT_RANK>,
    nCurrRank: SourceMutPointer<AT_RANK>,
    nTempRank: SourceMutPointer<AT_RANK>,
    nAtomNumber: SourceMutPointer<AT_NUMB>,
    bChanged: Option<&mut i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:1278 FixCanonEquivalenceInfo
    // INCHI✔️✔️: int FixCanonEquivalenceInfo( CANON_GLOBALS *pCG,
    // INCHI✔️✔️:                              int num_at_tg,
    // INCHI✔️✔️:                              AT_RANK *nSymmRank,
    // INCHI✔️✔️:                              AT_RANK *nCurrRank,
    // INCHI✔️✔️:                              AT_RANK *nTempRank,
    // INCHI✔️✔️:                              AT_NUMB *nAtomNumber,
    // INCHI✔️✔️:                              int *bChanged )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int nNumDiffRanks, bChangeSymmRank, bChangeCurrRank = 0;
    // INCHI✔️✔️:     /* sort equivalence information */
    // INCHI✔️✔️:     /*
    // INCHI✔️✔️:     int i;
    // INCHI✔️✔️:     for ( i = 0; i < num_at_tg; i ++ )
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nAtomNumber[i] = i;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     pCG->m_pn_RankForSort = nSymmRank; /* minimal class representatives: min ranks for equiv. atoms */
    // INCHI✔️✔️:     inchi_qsort( pCG, nAtomNumber, num_at_tg, sizeof( nAtomNumber[0] ), CompRanksOrd );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* convert equivalence information nSymmRank[] into ranks array nTempRank[] */
    // INCHI✔️✔️:     /* eq. info contains min. possible ranks for eq. atoms; nCurrRank contains max. possible ranks */
    // INCHI✔️✔️:     nNumDiffRanks = SortedEquInfoToRanks( nSymmRank/*inp*/, nTempRank/*out*/, nAtomNumber, num_at_tg, &bChangeSymmRank );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* check whether nCurrRank is same as new initial ranks calculated from nSymmRank[] */
    // INCHI✔️✔️:     bChangeCurrRank = memcmp( nCurrRank, nTempRank, num_at_tg * sizeof( nTempRank[0] ) );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*-----------------------------------------------------------------------
    // INCHI✔️✔️:     if ( bChangeSymmRank || bChangeCurrRank ) {
    // INCHI✔️✔️:          This is the case when the initial equitable partitioning does not produce
    // INCHI✔️✔️:          constitutionally equivalent classes of atoms.
    // INCHI✔️✔️:          Rebuild nSymmRank[] according to the new nCurrRank[] := nTempRank[]
    // INCHI✔️✔️:          For such structures the found canonical numbers of the constitutionally equivalent atoms
    // INCHI✔️✔️:          are not contiguous (see nCanonRank and nSymmRank examples below). Here arrays
    // INCHI✔️✔️:          nCurrRank, nAtomNumber, and nSymmRank are changed so that later the
    // INCHI✔️✔️:          contiguous canonical numbers for equivalent atoms can be obtained
    // INCHI✔️✔️:          (see GetCanonRanking under
    // INCHI✔️✔️:          "III. Get final canonical numbering (no stereo, no isotopic)".
    // INCHI✔️✔️:
    // INCHI✔️✔️:          Example: for CAS=37520-11-9 (ID=21247: Ethane, 1,2-dicyclopropyl-),
    // INCHI✔️✔️:
    // INCHI✔️✔️:                      the numbers are the "final canon. numbers, nCanonRank"
    // INCHI✔️✔️:           1
    // INCHI✔️✔️:
    // INCHI✔️✔️:           HC   7    5         3
    // INCHI✔️✔️:            | \
    // INCHI✔️✔️:            |  >CH--CH2        CH
    // INCHI✔️✔️:            | /       \      / |
    // INCHI✔️✔️:           HC        H2C--CH<  |
    // INCHI✔️✔️:                             \ |
    // INCHI✔️✔️:           2          6    8   CH
    // INCHI✔️✔️:
    // INCHI✔️✔️:                               4
    // INCHI✔️✔️:
    // INCHI✔️✔️:          the arrays (arranged according to ordering in nAtomNumberTemp) are:
    // INCHI✔️✔️:                                  before SortedEquInfoToRanks  after SortedRanksToEquInfo
    // INCHI✔️✔️:          orig. atom nos.,nAtomNumberTemp:  {4 5 6 7 0 1 2 3}   {4 5 6 7 0 1 2 3}
    // INCHI✔️✔️:          order numbers for sorted  ranks:  {0 1 2 3 4 5 6 7}   {0 1 2 3 4 5 6 7}
    // INCHI✔️✔️:          canonical numbering, nCanonRank:  {1 2 5 6 3 4 7 8}   {1 2 5 6 3 4 7 8}
    // INCHI✔️✔️:          constit. equivalence, nSymmRank:  {1 1 1 1 3 3 7 7}   {1 1 1 1 5 5 7 7} used later
    // INCHI✔️✔️:          initial equivalence,  nCurrRank:  {6 6 6 6 6 6 8 8}   {4 4 4 4 6 6 8 8} used later
    // INCHI✔️✔️:          initial numbering,  nAtomNumber:  {2 3 4 7 0 1 6 7}   {0 1 2 3 4 5 6 7} used later
    // INCHI✔️✔️:          final, no stereo, no isotopic, after  III. GetCanonRanking:
    // INCHI✔️✔️:          final canon. numbers, nCanonRank:                     {1 2 3 4 5 6 7 8} final
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     ----------------------------------------------------------------------------------*/
    // INCHI✔️✔️:     if (bChangeCurrRank)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         memcpy(nCurrRank, nTempRank, num_at_tg * sizeof(nCurrRank[0]));
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (bChangeSymmRank)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         SortedRanksToEquInfo( nSymmRank/*out*/, nTempRank/*inp*/, nAtomNumber, num_at_tg );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (bChanged)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         *bChanged = ( 0 != bChangeSymmRank ) | 2 * ( 0 != bChangeCurrRank );
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return nNumDiffRanks;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: FixCanonEquivalenceInfo

    let count = usize::try_from(num_at_tg).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    pCG.m_pn_RankForSort = nSymmRank.as_const();
    if count > 1 {
        heap.with_slice_mut_and_heap(nAtomNumber, |atom_numbers, heap| {
            let atom_numbers = atom_numbers
                .get_mut(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let byte_len = count
                .checked_mul(size_of::<AT_NUMB>())
                .ok_or(SourceHeapError::AllocationSizeOverflow)?;
            // AT_NUMB and the comparator's AT_RANK records are the same source uint16_t type.
            let bytes = unsafe {
                std::slice::from_raw_parts_mut(atom_numbers.as_mut_ptr().cast::<u8>(), byte_len)
            };
            inchi_qsort(bytes, count, size_of::<AT_NUMB>(), &mut |first, second| {
                let first = AT_RANK::from_ne_bytes(
                    first
                        .try_into()
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                );
                let second = AT_RANK::from_ne_bytes(
                    second
                        .try_into()
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                );
                CompRanksOrd(heap, first, second, pCG)
            })
        })?;
    }

    let mut bChangeSymmRank = 0_i32;
    let nNumDiffRanks = SortedEquInfoToRanks(
        heap,
        nSymmRank,
        nTempRank,
        nAtomNumber,
        num_at_tg,
        Some(&mut bChangeSymmRank),
    )?;

    let bChangeCurrRank = if count == 0 {
        false
    } else {
        let current = heap
            .slice(nCurrRank.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let temporary = heap
            .slice(nTempRank.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        current != temporary
    };
    if bChangeCurrRank {
        heap.with_slice_mut_and_heap(nCurrRank, |current, heap| {
            let current = current
                .get_mut(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let temporary = heap
                .slice(nTempRank.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            current.copy_from_slice(temporary);
            Ok(())
        })?;
    }
    if bChangeSymmRank != 0 {
        SortedRanksToEquInfo(heap, nSymmRank, nTempRank, nAtomNumber, num_at_tg)?;
    }
    if let Some(changed) = bChanged {
        *changed = i32::from(bChangeSymmRank != 0) | 2 * i32::from(bChangeCurrRank);
    }
    Ok(nNumDiffRanks)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn Canon_INChI3(
    heap: &mut SourceHeap,
    ic: &mut INCHI_CLOCK,
    clock_result: clock_t,
    user_action: Option<fn() -> i32>,
    console_quit: Option<fn() -> i32>,
    num_atoms: i32,
    num_at_tg: i32,
    at: SourceMutPointer<sp_ATOM>,
    pCS: &mut CANON_STAT,
    pCG: &mut CANON_GLOBALS,
    mut nMode: u32,
    bTautFtcn: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:1370 Canon_INChI3
    // INCHI✔️❌: complete raw source frame follows verbatim; Rust executes the
    // configured SEPARATE_CANON_CALLS=0 and FIX_STEREOCOUNT_ERR=1 branches.
    // The C source has undefined null dereferences when only nSymmRank qmalloc
    // fails and when NOEQ first-CT descriptor mapping writes nSymmStereo. Those
    // states do not define C behavior; Rust reports OOM for the former.
    /*
    int Canon_INChI3( INCHI_CLOCK *ic,
                     int num_atoms,
                     int num_at_tg,
                     sp_ATOM* at,
                     CANON_STAT* pCS,
                     CANON_GLOBALS *pCG,
                     INCHI_MODE nMode,
                     int bTautFtcn )
    {
    
    /****************************************************************
    
    0.    Initiation, Prepare initial ranks for GetCanonRanking()
    
    I.    Find constitutionally equivalent atoms and possibly canonical numbering
    I.1      Set tautomer=On, stereo=isotopic=Off
    I.2      GetCanonRanking(): Find constitutionally equivalent atoms and possibly canonical numbering
    1.3      Fix canonical equivalence info if needed (if the fix is needed then the numbering is not canonical)
    
    II.   Get final non-isotopic canonical numbering. Simultaneously obtain non-minimal isotopic and stereo CTs
             GetCanonRanking() with pCS->bKeepSymmRank = 1
             FillOutStereoParities() (create initial stereo descriptors)
             save non-isotopic canonicalization final results
             hide isotopic and tautomeric results (for historical reasons only)
    
    
    III.  Find constitutionally equivalent isotopic atoms (for isotopic stereo canonicalization)
    III.1    Allocate more memory
    III.2    fill allocated memory with the initial data
    III.3    duplicate, save old and add isotopic info to the new pCS->t_group_info
    III.4    Prepare initial isotopic ranks for GetCanonRanking()
    III.5    GetCanonRanking() to Find constitutionally equivalent ISOTOPIC atoms and tautomer groups
    III.6    Fix canonical isotopic equivalence information and derive ranks out of it
    
    IV.      Prepare a second Rank/AtomNumber Stack for mapping.
    
    V.    Optimize isotopic part (optimized)
             map_isotopic_atoms2()
             save isotopic canonical numbering
    
    VI.   Optimize stereo descriptors (optimized)
             map_stereo_bonds4()
    
    
    VII. Optimize isotopic stereo descriptors (optimized)
             SwitchAtomStereoAndIsotopicStereo()
             SetCtToIsotopicStereo()
             FillOutStereoParities()
             SetUseAtomForStereo()
             map_stereo_bonds4()
    
             SwitchAtomStereoAndIsotopicStereo()
             SetCtToNonIsotopicStereo()
    
    
    
    
    *****************************************************************/
    
        int     nRet = 0, i, n;
    
    
        /********************************************************
                  input non-stereo canonical info
         ********************************************************/
        BCN            *pBCN = pCS->pBCN;
        FTCN           *ftcn = pBCN->ftcn + bTautFtcn;
    
        /********************************************************
                  set mode flags
         ********************************************************/
    
        /* tautomeric structure */
        int bTaut = ( num_at_tg > num_atoms ) && pCS->t_group_info && pCS->t_group_info->num_t_groups && pCS->t_group_info->t_group;
    
        /* special case: induced by exchangable isotopic H inequivalence of atoms in formally non-tautomeric structure */
        int bIsoXchgH = pCS->t_group_info && pCS->t_group_info->nNumIsotopicEndpoints > 1 &&
            pCS->t_group_info->nIsotopicEndpointAtomNumber && pCS->t_group_info->nIsotopicEndpointAtomNumber[0] &&
            ( pCS->t_group_info->bTautFlagsDone & ( TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE ) )
            /* && (ftcn->nCanonFlags & CANON_FLAG_ISO_TAUT_DIFF)*/;
        int bHasIsotopicCanonData = ( ftcn->PartitionCtIso.AtNumber && ftcn->PartitionCtIso.Rank && ftcn->nSymmRankCtIso );
    
        /* bHasIsotopicCanonData==0 means
         *       (1) No isotopic atoms in the component OR
         *       (2) the component has only exchangable isotopic H that do not change canonical numbering and equivalence.
         */
        T_GROUP_INFO *t_group_info1 = bTaut ? pCS->t_group_info : NULL;
    
        /*int bIsoXchgH = t_group_info1 && t_group_info1->nNumIsotopicEndpoints && t_group_info1->nIsotopicEndpointAtomNumber;*/
        /* isotopic canonicalization */
        int bCanonIsotopic = bHasIsotopicCanonData && ( nMode & CMODE_ISO ) && ( pCS->LinearCTIsotopic || pCS->LinearCTIsotopicTautomer || bIsoXchgH );
    
        /* stereo canonicalization */
        int bCanonStereo = ( nMode & CMODE_STEREO ) && ( pCS->LinearCTStereoDble || pCS->LinearCTStereoCarb );
    
        /* stereo isotopic canonicalization */
        int bCanonIsoStereo = bHasIsotopicCanonData && ( nMode & CMODE_ISO_STEREO ) && ( pCS->LinearCTIsotopicStereoDble || pCS->LinearCTIsotopicStereoCarb ) && bCanonIsotopic;
        int bIsoTaut = ( bTaut && bCanonIsotopic );
    
        int            bIgnoreIsotopicInputGroups;
        int            bIgnoreIsotopicInputAtoms;
    
        AT_RANK       **pRankStack1 = pBCN->pRankStack;
        int             nRankStackLen = pBCN->nMaxLenRankStack;
        int             num_max = pBCN->num_max;     /* allocation lengths in *pRankStack1[] */
        NEIGH_LIST     *NeighList = ftcn->NeighList;
    
        int             nNumCurrRanks = 0;
        AT_RANK        *nTempRank = NULL;
    
        AT_RANK        *nSymmRank = NULL;
    
        AT_RANK        *nAtomNumber = NULL;
        AT_RANK        *nRank = NULL;
    
        AT_RANK       **pRankStack2 = NULL;
        AT_RANK        *nCanonRankStereo = NULL;
        AT_RANK        *nCanonRankStereoInv = NULL;
        AT_RANK        *nSymmStereo = NULL;
    
        AT_RANK        *nCanonRankIsotopicStereo = NULL;
        AT_RANK        *nCanonRankIsotopicStereoInv = NULL;
    
        CUR_TREE *cur_tree = NULL;
        CUR_TREE CurrentTree;
    
    
        /*AT_ISO_TGROUP  *LinearCTIsotopicTautomer  = NULL; */
    
    
        CANON_STAT   CS2;
        CANON_STAT* pCS2 = &CS2;
    
        inchiTime   ulStartTime, ulEndTime;
        /*=========== Mode Bits (low 8 bits, bit 0 is Least Significant Bit) ===========
    
          Mode      Bits       Description
           '0' c    0          Only one connection table canonicalization
           '1' C    1          Recalculate CT using fixed nSymmRank
           '2' i    1|2        Isotopic canonicalization (internal)
           '3' I    1|2|4      Isotopic canonicalization (output)
           '4' s    1|8        Stereo canonicalization
           '5' S    1|2|4|16   Stereo isotopic canonicalization
           '6' A    1|2|4|8|16 Output All
    
          --- high 8 bits ----
          --- obsolete, only historical interest. ------
          1-2 : 0 => at[i].init_rank from Morgan+NeighList
                1 => at[i].init_rank from Atom Invariants
                2 => at[i].init_rank from nSymmRank[]
                     (at[i].init_rank is included in LinearCT
                           depending on CT_ATOMID definition)
          3   : 1 => Get Stereo canonical info
          4   : 1 => Get Isotopic canonical info
          5   : 1 => Get Charge/Radical canonical info
        ==================================================================*/
        /*int             nOutputMode = 0;*/ /* obsolete */
    
    
        int bSwitchedAtomToIsotopic = 0;
    
    
        /* vABParityUnknown holds actual value of an internal constant signifying       */
        /* unknown parity: either the same as for undefined parity (default==standard)  */
        /*  or a specific one (non-std; requested by SLUUD switch).                     */
        int vABParityUnknown = AB_PARITY_UNDF;
        if (0 != ( nMode & REQ_MODE_DIFF_UU_STEREO ))
        {
            /* Make labels for unknown and undefined stereo different */
            vABParityUnknown = AB_PARITY_UNKN;
        }
    
        InchiTimeGet( &ulStartTime );
    
        *pCS2 = *pCS;  /* save input information and pointers to allocated memory */
    
        /* Set "ignore isotopic differences in tautomer groups" true */
        if (bTaut)
        {
            /* save request for isotopic tautomeric groups */
            bIgnoreIsotopicInputGroups = pCS->t_group_info->bIgnoreIsotopic;
            pCS->t_group_info->bIgnoreIsotopic = 1;
        }
        else
        {
            bIgnoreIsotopicInputGroups = 1;
        }
        /* Save request for isotopic name */
        bIgnoreIsotopicInputAtoms = pCS->bIgnoreIsotopic;
        /* set "ignore isotopic differences in atoms" true */
        pCS->bIgnoreIsotopic = 1;
    
    
        /* Save non-isotopic and isotopic canonicalization results */
        pCS->nCanonFlags = ftcn->nCanonFlags;
        /* 1. non-isotopic */
    
        /* linear CT, H */
        memcpy(pCS->LinearCT, ftcn->LinearCt, ftcn->nLenLinearCt * sizeof(pCS->LinearCT[0]));
        if (pCS->nNum_H && ftcn->nNumH)
        {
            for (i = 0; i < num_atoms; i++)
            {
                pCS->nNum_H[i] = /*(S_CHAR)*/( CHAR_MASK & ftcn->nNumH[i] );
            }
        }
        if (pCS->nNum_H_fixed && ftcn->nNumHFixH)
        {
            for (i = 0; i < num_atoms; i++)
            {
                pCS->nNum_H_fixed[i] = /*(S_CHAR)*/( CHAR_MASK & ftcn->nNumHFixH[i] );
            }
        }
        pCS->nLenLinearCT = ftcn->nLenLinearCt;
        pCS->nLenLinearCTAtOnly = ftcn->nLenLinearCtAtOnly;
    
        /* save non-isotopic atoms equivalence and numbering */
        if (pCS->nSymmRank)
        {
            memcpy(pCS->nSymmRank, ftcn->nSymmRankCt, num_at_tg * sizeof(pCS->nSymmRank[0]));
        }
        if (pCS->nCanonOrd)
        {
            memcpy(pCS->nCanonOrd, ftcn->PartitionCt.AtNumber, num_at_tg * sizeof(pCS->nCanonOrd[0]));
            pCS->nLenCanonOrd = num_atoms;
        }
        if (ftcn->iso_exchg_atnos && pCS->nExchgIsoH)
        {
            for (i = 0; i < num_atoms; i++)
            {
                pCS->nExchgIsoH[i] = !ftcn->iso_exchg_atnos[i]; /* (pCS->nExchgIsoH[i]==1) => tautomeric or hetero atoms that may exchange isotopic H */
            }
        }
        /* 2. isotopic */
    
        if (bCanonIsotopic)
        {
            /* linear CT, num_H are same as non-isotopic */
            /* save atoms equivalence and numbering */
            if (pCS->nSymmRankIsotopic)
            {
                memcpy(pCS->nSymmRankIsotopic, ftcn->nSymmRankCtIso, num_at_tg * sizeof(pCS->nSymmRankIsotopic[0]));
            }
            if (pCS->nCanonOrdIsotopic)
            {
                memcpy(pCS->nCanonOrdIsotopic, ftcn->PartitionCtIso.AtNumber, num_at_tg * sizeof(pCS->nCanonOrdIsotopic[0]));
                pCS->nLenCanonOrdIsotopic = num_at_tg;
            }
    
            nRet = FillIsotopicAtLinearCT( num_atoms, at,
                                           ftcn->PartitionCtIso.AtNumber,
                                           pCS->LinearCTIsotopic,
                                           pCS->nMaxLenLinearCTIsotopic,
                                           &pCS->nLenLinearCTIsotopic );
    
            if (RETURNED_ERROR( nRet ))
            {
                goto exit_function;
            }
            if (nRet < 0)
            {
                nRet = CT_TAUCOUNT_ERR;
                goto exit_function;
            }
        }
        else
        {
            pCS->nMaxLenLinearCTIsotopic = 0;
            pCS->nMaxLenLinearCTIsotopicTautomer = 0;
        }
    
        /* fill out tautomeric groups, isotopic and non-isotopic tautomeric CT and t_group_info1->tGroupNumber */
        if (bTaut)
        {
            bIsoTaut = bIsoTaut && ftcn->PartitionCtIso.Rank &&
                ftcn->PartitionCtIso.AtNumber && ftcn->nSymmRankCtIso;
    
            nRet = FillTautLinearCT2( pCG, num_atoms, num_at_tg, bIsoTaut,
                    ftcn->PartitionCt.Rank, ftcn->PartitionCt.AtNumber, ftcn->nSymmRankCt,
                    ftcn->PartitionCtIso.Rank, ftcn->PartitionCtIso.AtNumber, ftcn->nSymmRankCtIso,
                    pCS->LinearCTTautomer, pCS->nMaxLenLinearCTTautomer, &pCS->nLenLinearCTTautomer,
                    pCS->LinearCTIsotopicTautomer, pCS->nMaxLenLinearCTIsotopicTautomer, &pCS->nLenLinearCTIsotopicTautomer,
                    t_group_info1 );
    
            if (RETURNED_ERROR( nRet ))
            {
                goto exit_function;
            }
            if (nRet <= 0)
            {
                nRet = CT_TAUCOUNT_ERR;
                goto exit_function;
            }
            else
            {
                /* tautomeric groups: save non-isotopic symmetry & t_group order */
                int num_t_groups = t_group_info1->num_t_groups;
                AT_NUMB *tGroupNumber = t_group_info1->tGroupNumber;
                AT_NUMB *tSymmRank = tGroupNumber + TGSO_SYMM_RANK*num_t_groups;
                if (pCS->nSymmRankTaut)
                {
                    memcpy(pCS->nSymmRankTaut, tSymmRank, num_t_groups * sizeof(pCS->nSymmRank[0])); /* fixed 5-23-02 */
                }
                if (pCS->nCanonOrdTaut)
                {
                    memcpy(pCS->nCanonOrdTaut, tGroupNumber, num_t_groups * sizeof(pCS->nCanonOrdTaut[0]));
                    pCS->nLenCanonOrdTaut = num_t_groups;
                }
                if (bCanonIsotopic /*&& pCS->nLenLinearCTIsotopicTautomer*/)
                {
                    /* tautomeric groups: save isotopic symmetry & t_group order */
                    /*AT_NUMB ntRankOffset       = (AT_RANK)num_atoms;*/
                    AT_NUMB *tiSymmRank = tGroupNumber + TGSO_SYMM_IRANK*(long long)num_t_groups; /* djb-rwth: cast operator added */
                    AT_NUMB *tiGroupNumber = tGroupNumber + TGSO_SYMM_IORDER*(long long)num_t_groups; /* djb-rwth: cast operator added */
                    if (pCS->nSymmRankIsotopicTaut)
                    {
                        memcpy(pCS->nSymmRankIsotopicTaut, tiSymmRank, num_t_groups * sizeof(pCS->nSymmRankIsotopicTaut[0]));
                    }
                    memcpy(pCS->nCanonOrdIsotopicTaut, tiGroupNumber, num_t_groups * sizeof(pCS->nCanonOrdIsotopicTaut[0]));
                    pCS->nLenCanonOrdIsotopicTaut = num_t_groups;
                }
            }
        }
        /* save connection table if requested */
        if (pCS->LinearCT2)
        {
            memcpy(pCS->LinearCT2, pCS->LinearCT, sizeof(pCS->LinearCT2[0])* pCS->nLenLinearCT);
            pCS->nLenLinearCT2 = pCS->nLenLinearCT;
            pCS->nLenLinearCTAtOnly2 = pCS->nLenLinearCTAtOnly;
        }
    
        if (num_atoms <= 1)
        {
            bCanonStereo = 0;  /* a sinle atom + possibly terminal hydrogen atoms */
            if (num_atoms < 1 || !at[0].parity2)
            {
                bCanonIsoStereo = 0;  /*  structure; for example Cl- or CH4 */
            }
        }
    
        if (!bCanonStereo && !( bCanonIsotopic && bCanonIsoStereo ))
        {
            goto exit_function; /* skip stereo canonicalization */
        }
    
    
    
        /**********************************************************
                         Mode
        ***********************************************************/
        nMode = nMode & CANON_MODE_MASK;
    
        /* memory allocation */
    
        nAtomNumber = (AT_RANK *) qmalloc( num_max * sizeof( *nAtomNumber ) );
        nRank = (AT_RANK *) qmalloc( num_max * sizeof( *nRank ) );
        nTempRank = (AT_RANK *) qmalloc( num_max * sizeof( *nTempRank ) );
        nSymmRank = (AT_RANK *) qmalloc( num_max * sizeof( *nSymmRank ) );
    
        /***********************************************
                    0.1 Initialization
        ************************************************/
    
    
        if (!NeighList || !nAtomNumber || !nTempRank ||
             !nRank || !pCS->LinearCT)
        {
            nRet = CT_OUT_OF_RAM;  /* program error */  /*  <BRKPT> */
            goto exit_function;
        }
    
        pCS->NeighList = NeighList;
    
        *pCS2 = *pCS;  /* save input information and pointers to allocated memory */
    
        if (!( nMode & CMODE_NOEQ_STEREO ) && ( bCanonStereo || bCanonIsoStereo ))
        {
            /* will be used to discover vertex equivalences in stereo canonicalization */
            memset( &CurrentTree, 0, sizeof( CurrentTree ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            cur_tree = &CurrentTree;
        }
    
    
        pCS->bCmpStereo = 0;
        pCS->bCmpIsotopicStereo = 0;
    
    
        if (bCanonStereo || bCanonIsoStereo)
        {
            int ii, nn;
    
            /* stereo or isotopic canonicalization: we need a second set of ranks for mapping */
            /* (isotopic atoms or stereo can only increase nNumCurrRanks) */
            pRankStack2 = (AT_RANK **) inchi_calloc( nRankStackLen, sizeof( AT_RANK * ) );
            if (pRankStack2)
            {
                /* prepare for ranks reuse */
                for (nn = 2; nn < nRankStackLen && pRankStack1[nn]; nn++)
                {
                    pRankStack1[nn][0] = 0; /* means ranks have to be calculated */
                }
                /* reuse memory to reduce number of allocations: */
                /* move last half of pointers from pRankStack1 to pRankStack2 */
                /* The first 2 elements will be assigned separately */
                if (( nn = ( nn - 2 ) / 2 ) > 0)
                {
                    for (ii = 2 + nn; ii < nRankStackLen && pRankStack1[ii]; ii++)
                    {
                        pRankStack2[ii - nn] = pRankStack1[ii];
                        pRankStack1[ii] = NULL;
                    }
                }
            }
            else
            {
                nRet = CT_OUT_OF_RAM;  /*  <BRKPT> */
                goto exit_function; /* program error */
            }
        }
    
        if (bCanonStereo)
        {
    
           /* *pCS2 = *pCS; */ /* save input information and pointers to allocated memory */
    
            /* initial ranking for non-isotopic mapping */
            memcpy(nAtomNumber, ftcn->PartitionCt.AtNumber, num_at_tg * sizeof(nAtomNumber[0]));
            memcpy(nRank, ftcn->PartitionCt.Rank, num_at_tg * sizeof(nRank[0]));
            memcpy(nSymmRank, ftcn->nSymmRankCt, num_at_tg * sizeof(nSymmRank[0]));
    
            /* nSymmRank changes if canonical numbers of constitutionally equivalent atoms are not contiguous */
            nNumCurrRanks = FixCanonEquivalenceInfo( pCG, num_at_tg,
                                                     nSymmRank /* in&out*/,
                                                     nRank, nTempRank /* out */,
                                                     nAtomNumber /* in&out */,
                                                     NULL );
            /* atom numbers in canonical order */
            memcpy(pCS->nPrevAtomNumber, ftcn->PartitionCt.AtNumber, num_at_tg * sizeof(nAtomNumber[0]));
    
            /* fill stereo part of the connection table with initial (not optimized) parities */
            /* input
            pCS->LinearCTStereoDble
            pCS->LinearCTStereoCarb
            pCS->nMaxLenLinearCTStereoCarb
            pCS->nMaxLenLinearCTStereoDble
            */
    
            nRet = FillOutStereoParities( at, num_atoms, ftcn->PartitionCt.Rank, ftcn->PartitionCt.AtNumber,
                                          nRank, nAtomNumber, pCS, pCG, 0 /* bIsotopic */ );
    
            /* output
            pCS->LinearCTStereoDble
            pCS->LinearCTStereoCarb
            pCS2->nLenLinearCTStereoCarb
            pCS2->nLenLinearCTStereoDble
            */
    
            if (RETURNED_ERROR( nRet ))
            {
                goto exit_function;
            }
            if (nRet < 0)
            {
                nRet = CT_STEREOCOUNT_ERR;
                goto exit_function;
            }
    
            /***************************************************************
             *
             *  VI. Optimize non-isotopic stereo descriptors (optimized)
             *
             ***************************************************************/
    
            /* allocate memory for stereo canonicalization */
    
            if (!nCanonRankStereo)
            {
                nCanonRankStereo = (AT_RANK *) qmalloc( num_max * sizeof( *nCanonRankStereo ) );
            }
            if (!nSymmStereo && !( nMode & CMODE_NOEQ_STEREO ))
            {
                nSymmStereo = (AT_RANK *) qmalloc( ( (long long)num_max + 1 ) * sizeof( *nSymmStereo ) ); /* djb-rwth: cast operator added */
            }
    
            if (!( nMode & CMODE_NOEQ_STEREO ) && 0 > CurTreeAlloc( cur_tree, num_at_tg ))
            {
                nRet = CT_OUT_OF_RAM;  /*  <BRKPT> */
                goto exit_function;
            }
            /* check allocations and assign first 2 elements of pRankStack2 */
            if (pRankStack1 && pRankStack2 &&
                 nCanonRankStereo &&
                 /* nCurrRankStereo  && nAtomNumberCurrStereo &&*/
                 ( nSymmStereo || ( nMode & CMODE_NOEQ_STEREO ) ))
            {
                pRankStack1[0] = pRankStack2[0] = nRank;
                pRankStack1[1] = pRankStack2[1] = nAtomNumber;
            }
            else
            {
                nRet = CT_OUT_OF_RAM;  /*  <BRKPT> */
                goto exit_function;
            }
    
            /****************************************************************
             *
             *  VI-A. Optimize non-isotopic non-inverted stereo descriptors
             *
             ****************************************************************/
    
            /* set the 1st ranks in the rest of the stack to zero: prepare for ranks reuse */
            for (n = 2; n < nRankStackLen && pRankStack1[n]; n++)
            {
                pRankStack1[n][0] = 0; /* means ranks have to be recalculated */
            }
            /* set the 1st ranks to zero: prepare for ranks reuse */
            for (n = 2; n < nRankStackLen && pRankStack2[n]; n++)
            {
                pRankStack2[n][0] = 0; /* means ranks have to be recalculated */
            }
    
            /* for debugging or statistics */
            pCS->lNumBreakTies =
                pCS->lNumNeighListIter =
                pCS->lNumTotCT =
                pCS->lNumDecreasedCT =
                pCS->lNumRejectedCT =
                pCS->lNumEqualCT = 0;
            pCS->bKeepSymmRank = 0;
            pCS->bFirstCT = 1; /* To fill out nCanonRankStereo[] in map_stero_atoms2() */
    
            /******************************************************************************
                 nCanonRank contains input canonical numbering
                 nCanonRankStereo will be filled with a transposition of canonical numbering
                   which (1) keeps connection table unchanged and
                         (2) provides minimal stereo descriptors in
                             pCS->LinearCTStereoDble (length=pCS->nLenLinearCTStereoDble)
                             pCS->LinearCTStereoCarb (length=pCS->nLenLinearCTStereoCarb)
             */
    
            nRet = map_stereo_bonds4( ic, pCG, at, num_atoms, num_at_tg, num_max, 0,
                                       ftcn->PartitionCt.Rank, ftcn->PartitionCt.AtNumber,
                                       nCanonRankStereo, nSymmRank,
                                       pRankStack1, pRankStack2,
                                       nTempRank, nNumCurrRanks,nSymmStereo,
                                       NeighList, pCS, cur_tree, 0 /* nNumMappedBonds */,
                                       vABParityUnknown );
    
            if (RETURNED_ERROR( nRet ))
            {
                if (nRet == CT_TIMEOUT_ERR)
                {
                    goto exit_function;
                }
                else
                {
                    goto exit_function; /* program error */
                }
            }
            else
            {
                int bFailed = 0;
                if (!nRet)
                {
                    /* djb-rwth: removing redundant code */
                    pCS2->nLenLinearCTStereoCarb =
                        pCS->nLenLinearCTStereoCarb = -abs( pCS->nLenLinearCTStereoCarb );
                    pCS2->nLenLinearCTStereoDble =
                        pCS->nLenLinearCTStereoDble = -abs( pCS->nLenLinearCTStereoDble );
                    nRet = CT_STEREOCOUNT_ERR;  /*  <BRKPT> */
                    goto exit_function; /* program error */
                }
                else
                {
                    /* save non-isotopic lengths */
                    pCS2->nLenLinearCTStereoDble = pCS->nLenLinearCTStereoDble;
                    pCS2->nLenLinearCTStereoCarb = pCS->nLenLinearCTStereoCarb;
                    /* djb-rwth: removing redundant code */
                }
    
                /* save stereo canonical numbering */
                if (pCS->nCanonOrdStereo)
                {
                    for (i = 0; i < num_at_tg; i++) /* djb-rwth: removing redundant code */
                    {
                        if (nCanonRankStereo[i] && (int) nCanonRankStereo[i] <= num_at_tg)
                        {
                            pCS->nCanonOrdStereo[(int) nCanonRankStereo[i] - 1] = (AT_NUMB) i;
                        }
                        else
                        {
                            bFailed++;
                        }
                    }
                    pCS->nLenCanonOrdStereo = ( bFailed ) ? -num_atoms : num_atoms;
                }
    
                /* save stereo tautomer groups numbering */
                if (bTaut && pCS->nCanonOrdStereoTaut)
                {
                    if (0 < ( nRet = SortTautomerGroupsAndEndpoints( pCG, t_group_info1, num_atoms, num_at_tg, nCanonRankStereo ) ))
                    {
                        /*non-isotopic contains symmetry ranks */
                        int num_t_groups = t_group_info1->num_t_groups;
                        AT_NUMB *tGroupNumber = t_group_info1->tGroupNumber;
                        /*AT_NUMB *tiSymmRank        = tGroupNumber + TGSO_SYMM_IRANK*num_t_groups; */
                        memcpy(pCS->nCanonOrdStereoTaut, tGroupNumber, num_t_groups * sizeof(pCS->nCanonOrdStereoTaut[0]));
                        pCS->nLenCanonOrdStereoTaut = ( bFailed ) ? -num_t_groups
                            : num_t_groups;
                    }
                    else
                    {
                        if (RETURNED_ERROR(nRet))
                        {
                            goto exit_function;
                        }
                    }
                    /* djb-rwth: removing redundant code */
                    /*SortTautomerGroupsAndEndpoints( t_group_info1, nCanonRank ); */ /* ??? return to non-isotopic canonical numbering */
                }
            }
    
    
            /****************************************************
             *
             *  VI-B. Optimize INVERTED stereo descriptors
             *
             ****************************************************/
            if (!nCanonRankStereoInv)
            {
                nCanonRankStereoInv = (AT_RANK *) qmalloc( num_max * sizeof( *nCanonRankStereoInv ) );
            }
            if (!nCanonRankStereoInv)
            {
                nRet = CT_OUT_OF_RAM;  /*  <BRKPT> */
                goto exit_function;
            }
    
            /* copy previous non-isotopic stereo canonicalization results to Inv initial data */
            /* assign pointers */
            pCS->LinearCTStereoDble = pCS2->LinearCTStereoDbleInv;
            pCS->LinearCTStereoCarb = pCS2->LinearCTStereoCarbInv;
    
            /* copy the lengths */
            pCS2->nLenLinearCTStereoDbleInv =
                pCS->nLenLinearCTStereoDbleInv =
                pCS->nLenLinearCTStereoDble = pCS2->nLenLinearCTStereoDble;
    
            pCS2->nLenLinearCTStereoCarbInv =
                pCS->nLenLinearCTStereoCarbInv =
                pCS->nLenLinearCTStereoCarb = pCS2->nLenLinearCTStereoCarb;
    
            if (pCS->nLenLinearCTStereoDble > 0 || pCS->nLenLinearCTStereoCarb > 0)
            {
                /* copy previous results, the canonical stereo CT */
                memcpy(pCS->LinearCTStereoDble, pCS2->LinearCTStereoDble, pCS->nLenLinearCTStereoDble * sizeof(pCS->LinearCTStereoDble[0]));
                memcpy(pCS->LinearCTStereoCarb, pCS2->LinearCTStereoCarb, pCS->nLenLinearCTStereoCarb * sizeof(pCS->LinearCTStereoCarb[0]));
            }
            memcpy(nCanonRankStereoInv, nCanonRankStereo, num_max * sizeof(nCanonRankStereoInv[0]));
            if (pCS->nCanonOrdStereoInv && pCS->nCanonOrdStereo)
            {
                /* in case there is nothing to invert */
                memcpy(pCS->nCanonOrdStereoInv, pCS->nCanonOrdStereo, num_at_tg * sizeof(pCS->nCanonOrdStereoInv[0]));
            }
    
            /******************************
             *
             * Invert stereo
             *
             ******************************/
    
            /*********************************************************************************
             * Create initial approximation for the minimization of the stereo descriptors:
             *  invert stereogenic atom parities, one parity in each allene, all parities in
             *  pCS->LinearCTStereoCarb and allene parities in pCS->nLenLinearCTStereoDble
             */
            nRet = InvertStereo( at, num_at_tg, nCanonRankStereo, nTempRank, pCS, 1 /* bInvertLinearCTStereo */ );
            if (RETURNED_ERROR( nRet ))
            {
                goto exit_function;
            }
            else if (nRet > 0)
            {
                /* InvertStereo() has done some changes */
                /* djb-rwth: removing redundant code */
                /* FillOutStereoParities() has already been called to fill out these 2 LinearCTs */
    
                /* set the 1st ranks in the rest of the stack to zero: prepare for ranks reuse */
                for (n = 2; n < nRankStackLen && pRankStack1[n]; n++)
                {
                    pRankStack1[n][0] = 0; /* means ranks have to be recalculated */
                }
                /* set the 1st ranks to zero: prepare for ranks reuse */
                for (n = 2; n < nRankStackLen && pRankStack2[n]; n++)
                {
                    pRankStack2[n][0] = 0; /* means ranks have to be recalculated */
                }
    
                /* for debugging or statistics */
                pCS->lNumBreakTies =
                    pCS->lNumNeighListIter =
                    pCS->lNumTotCT =
                    pCS->lNumDecreasedCT =
                    pCS->lNumRejectedCT =
                    pCS->lNumEqualCT = 0;
                pCS->bKeepSymmRank = 0;
                pCS->bFirstCT = 1; /* To fill out nCanonRankStereo[] in map_stero_atoms2() */
    
    #ifdef FIX_STEREOCOUNT_ERR
                CurTreeSetPos( cur_tree, 0 );
    #endif
    
                /******************************************************************************
                     ftcn->PartitionCt.Rank contains input canonical numbering
                     nCanonRankStereoInv will be filled with a transposition of canonical numbering
                       which (1) keeps connection table unchanged and
                             (2) provides minimal stereo descriptors in
                                 pCS->LinearCTStereoDble (length=pCS->nLenLinearCTStereoDble)
                                 pCS->LinearCTStereoCarb (length=pCS->nLenLinearCTStereoCarb)
                 ******************************************************************************/
    
                nRet = map_stereo_bonds4( ic,pCG, at, num_atoms, num_at_tg, num_max, 0,
                                          ftcn->PartitionCt.Rank,
                                          ftcn->PartitionCt.AtNumber,
                                          nCanonRankStereoInv,
                                          nSymmRank, pRankStack1, pRankStack2,
                                          nTempRank, nNumCurrRanks, nSymmStereo,
                                          NeighList, pCS, cur_tree, 0,
                                          vABParityUnknown );
    
                if (RETURNED_ERROR( nRet ))
                {
                    if (nRet == CT_TIMEOUT_ERR)
                    {
                        goto exit_function;
                    }
                    else
                    {
                        goto exit_function; /* program error */
                    }
                }
                else
                {
                    int bFailed = 0;
                    if (!nRet)
                    {
                        /* djb-rwth: removing redundant code */
                        pCS2->nLenLinearCTStereoCarb =
                            pCS->nLenLinearCTStereoCarb = -abs( pCS->nLenLinearCTStereoCarb );
                        pCS2->nLenLinearCTStereoDble =
                            pCS->nLenLinearCTStereoDble = -abs( pCS->nLenLinearCTStereoDble );
                        nRet = CT_STEREOCOUNT_ERR;  /*  <BRKPT> */
                        goto exit_function; /* program error */
                    }
    
                    /* save non-isotopic pointers & lengths for INVERTED stereo */
                    pCS->nLenLinearCTStereoDbleInv =
                        pCS2->nLenLinearCTStereoDbleInv = pCS->nLenLinearCTStereoDble;
                    pCS->nLenLinearCTStereoCarbInv =
                        pCS2->nLenLinearCTStereoCarbInv = pCS->nLenLinearCTStereoCarb;
    
                        /* restore pointers and lengths to non-inverted stereo    */
                        /*  -- this is needed for InvertStereo() back, see below  */
                    pCS->LinearCTStereoDble = pCS2->LinearCTStereoDble;
                    pCS->LinearCTStereoCarb = pCS2->LinearCTStereoCarb;
                    pCS->nLenLinearCTStereoDble = pCS2->nLenLinearCTStereoDble;
                    pCS->nLenLinearCTStereoCarb = pCS2->nLenLinearCTStereoCarb;
                    /* consistency check */
                    if (pCS->nLenLinearCTStereoDbleInv != pCS->nLenLinearCTStereoDble ||
                         pCS->nLenLinearCTStereoCarbInv != pCS->nLenLinearCTStereoCarb)
                    {
                        nRet = CT_CALC_STEREO_ERR;
                        goto exit_function; /* program error */
                    }
    
                    /******************************
                     *
                     * Invert stereo back
                     *
                     ******************************
                     *  (make sure that pointers
                     *  pCS->LinearCTStereoCarb,
                     *  pCS->LinearCTStereoDble
                     *  and corresponding lengths
                     *  have been restored)
                     ******************************/
                /*********************************************************************************
                 *  invert only stereogenic atom parities and one parity in each allene, DO NOT
                 *  change parities in pCS->LinearCTStereoCarb and pCS->nLenLinearCTStereoDble
                 */
                    nRet = InvertStereo( at,
                                         num_at_tg,
                                         nCanonRankStereo,
                                         nTempRank,
                                         pCS,
                                         0 );
    
                    if (RETURNED_ERROR( nRet ))
                    {
                        goto exit_function;
                    }
                    nRet = 0;
    
    
                    /* save stereo canonical numbering */
                    if (pCS->nCanonOrdStereoInv)
                    {
                        for (i = 0; i < num_at_tg; i++) /* djb-rwth: removing redundant code */
                        {
                            if (nCanonRankStereoInv[i] && (int) nCanonRankStereoInv[i] <= num_at_tg)
                            {
                                pCS->nCanonOrdStereoInv[(int) nCanonRankStereoInv[i] - 1] = (AT_NUMB) i;
                            }
                            else
                            {
                                bFailed++;
                            }
                        }
                        pCS->nLenCanonOrdStereo = ( bFailed ) ? -num_atoms : num_atoms;
                    }
    
                    /* compare inverted and non-inverted stereo */
                    pCS->bCmpStereo = 2 + CompareLinCtStereo(
                                             pCS->LinearCTStereoDbleInv, pCS->nLenLinearCTStereoDbleInv,
                                             pCS->LinearCTStereoCarbInv, pCS->nLenLinearCTStereoCarbInv,
                                             pCS->LinearCTStereoDble, pCS->nLenLinearCTStereoDble,
                                             pCS->LinearCTStereoCarb, pCS->nLenLinearCTStereoCarb
                    );
                }
            }
            else if (0 == nRet)
            {
                /* nothing has been done, restore pointers and lengths for stereo */
                pCS->LinearCTStereoDble = pCS2->LinearCTStereoDble;
                pCS->LinearCTStereoCarb = pCS2->LinearCTStereoCarb;
                pCS->nLenLinearCTStereoDble = pCS2->nLenLinearCTStereoDble;
                pCS->nLenLinearCTStereoCarb = pCS2->nLenLinearCTStereoCarb;
            }
        }
    
        /* restore "ignore isotopic differences in tautomer groups" */
        if (bTaut)
        {
            /* save request for isotopic tautomeric groups */
            pCS->t_group_info->bIgnoreIsotopic = bIgnoreIsotopicInputGroups;
        }
    
        /* restore request for isotopic name */
        pCS->bIgnoreIsotopic = bIgnoreIsotopicInputAtoms;
    
        if (bCanonIsoStereo && bCanonIsotopic)
        {
    
            /****************************************************************
             *
             *   VII. Optimize isotopic stereo descriptors (optimized)
             *
             ****************************************************************/
            /*
            pCS->LinearCTIsotopic     = NULL;
            */
    
            /* Initial ranking for isotopic mapping */
            memcpy(nAtomNumber, ftcn->PartitionCtIso.AtNumber, num_at_tg * sizeof(nAtomNumber[0]));
            memcpy(nRank, ftcn->PartitionCtIso.Rank, num_at_tg * sizeof(nRank[0]));
            memcpy(nSymmRank, ftcn->nSymmRankCtIso, num_at_tg * sizeof(nSymmRank[0]));
            /* nSymmRank will change if canonical numbers of of constitutionally equivalent atoms are not contiguous */
            nNumCurrRanks = FixCanonEquivalenceInfo( pCG, num_at_tg, nSymmRank /* in&out*/,
                                               nRank, nTempRank /* out */, nAtomNumber /* in&out */, NULL );
            memcpy(pCS->nPrevAtomNumber, ftcn->PartitionCtIso.AtNumber, num_at_tg * sizeof(nAtomNumber[0]));
            /* Allocate memory for optimized stereo canonicalization */
            /* for stereo canonical numbering to be found. */
            if (!nCanonRankIsotopicStereo)
            {
                nCanonRankIsotopicStereo = (AT_RANK *) qmalloc( num_max * sizeof( *nCanonRankIsotopicStereo ) );
            }
            if (!nSymmStereo && !( nMode & CMODE_NOEQ_STEREO ))
            {
                nSymmStereo = (AT_RANK *) qmalloc( ( (long long)num_max + 1 ) * sizeof( *nSymmStereo ) ); /* djb-rwth: cast operator added */
            }
            if (!( nMode & CMODE_NOEQ_STEREO ) && CurTreeAlloc( cur_tree, num_at_tg ))
            {
                nRet = CT_OUT_OF_RAM;  /*  <BRKPT> */
                goto exit_function;
            }
    
            /* Check allocations and assign first 2 elements of pRankStack2 */
            if (pRankStack1 && pRankStack2 &&
                 nCanonRankIsotopicStereo &&
                 ( nSymmStereo || ( nMode & CMODE_NOEQ_STEREO ) ))
            {
    
                pRankStack1[0] = pRankStack2[0] = nRank; /* pRankStack1[0,1] shall be unchanged */
                pRankStack1[1] = pRankStack2[1] = nAtomNumber;
            }
            else
            {
                nRet = CT_OUT_OF_RAM;  /*  <BRKPT> */
                goto exit_function;
            }
    
            /******************************************************************
               Important: fill out a list of stereo atoms and bonds including
               those which are stereo due to isotopic atoms only and create
               LinearCT stereo descriptors for the canonical numbering
             ******************************************************************/
    
    
            /* at[] has certain members for non-isotopic and isotopic stereo; switch them */
            SwitchAtomStereoAndIsotopicStereo( at, num_atoms, &bSwitchedAtomToIsotopic );
    
            /* Prepare stereo connection tables' pointers */
            SetCtToIsotopicStereo( pCS, pCS2 );
    
            nRet = FillOutStereoParities( at, num_atoms,
                                          ftcn->PartitionCtIso.Rank,
                                          ftcn->PartitionCtIso.AtNumber,
                                          nRank, nAtomNumber,
                                          pCS, pCG, 1 /* bIsotopic */ );
    
    
            if (RETURNED_ERROR( nRet ))
            {
                goto exit_function;  /* program error */
            }
            else
            {
                if (!nRet)
                {
                    /* no isotopic stereo */
                    pCS2->nLenLinearCTIsotopicStereoDble = pCS->nLenLinearCTIsotopicStereoDble = 0;
                    pCS2->nLenLinearCTIsotopicStereoCarb = pCS->nLenLinearCTIsotopicStereoCarb = 0;
                    pCS->nLenCanonOrdIsotopicStereo = 0;
                    pCS->nLenCanonOrdIsotopicStereoTaut = 0;
                    pCS2->nLenLinearCTIsotopicStereoDbleInv = pCS->nLenLinearCTIsotopicStereoDbleInv = 0;
                    pCS2->nLenLinearCTIsotopicStereoCarbInv = pCS->nLenLinearCTIsotopicStereoCarbInv = 0;
                    goto bypass_isotopic_stereo;
                }
            }
    
            /* djb-rwth: removing redundant code */
    
    
            /*************************************************************
             *
             *  VII-A. Optimize non-inverted isotopic stereo descriptors
             *
             *************************************************************/
    
            /* set the 1st ranks in the rest of the stack to zero: prepare for ranks reuse */
            for (n = 2; n < nRankStackLen && pRankStack1[n]; n++)
            {
                pRankStack1[n][0] = 0; /* means ranks have to be recalculated */
            }
            /* set the 1st ranks to zero: prepare for ranks reuse */
            for (n = 2; n < nRankStackLen && pRankStack2[n]; n++)
            {
                pRankStack2[n][0] = 0; /* means ranks have to be recalculated */
            }
    
            /* for debugging or statistics */
            pCS->lNumBreakTies =
                pCS->lNumNeighListIter =
                pCS->lNumTotCT =
                pCS->lNumDecreasedCT =
                pCS->lNumRejectedCT =
                pCS->lNumEqualCT = 0;
            pCS->bKeepSymmRank = 0;
            pCS->bFirstCT = 1; /* To fill out nCanonRankStereo[] in map_stero_atoms2() */
    
            /**************************************************************************************
              nCanonRankIsotopic contains input canonical numbering
              nCanonRankIsotopicStereo will be filled with a transposition of canonical numbering
                that  (1) keeps connection table unchanged and
                      (2) provides minimal stereo descriptors in
                          pCS->LinearCTStereoDble (length=pCS->nLenLinearCTStereoDble)
                          pCS->LinearCTStereoCarb (length=pCS->nLenLinearCTStereoCarb)
            ***************************************************************************************/
    
            nRet = map_stereo_bonds4( ic, pCG,at, num_atoms, num_at_tg, num_max, 0,
                                      ftcn->PartitionCtIso.Rank,
                                      ftcn->PartitionCtIso.AtNumber,
                                      nCanonRankIsotopicStereo,
                                      nSymmRank, pRankStack1, pRankStack2,
                                      nTempRank, nNumCurrRanks, nSymmStereo,
                                      NeighList, pCS, cur_tree,
                                      0, vABParityUnknown );
    
            if (RETURNED_ERROR( nRet ))
            {
                goto exit_function;
            }
            else
            {
                int bFailed = 0;
    
                if (!nRet)
                {
                    /* djb-rwth: removing redundant code */
                    pCS2->nLenLinearCTIsotopicStereoDble =
                        pCS->nLenLinearCTIsotopicStereoDble = -abs( pCS->nLenLinearCTStereoDble );
                    pCS2->nLenLinearCTIsotopicStereoCarb =
                        pCS->nLenLinearCTIsotopicStereoCarb = -abs( pCS->nLenLinearCTStereoCarb );
                    nRet = CT_STEREOCOUNT_ERR;  /*  <BRKPT> */
                    goto exit_function; /* program error */
                }
                else
                {
                    /* save isotopic lengths */
                    pCS->nLenLinearCTIsotopicStereoDble =
                        pCS2->nLenLinearCTIsotopicStereoDble = pCS->nLenLinearCTStereoDble;
                    pCS->nLenLinearCTIsotopicStereoCarb =
                        pCS2->nLenLinearCTIsotopicStereoCarb = pCS->nLenLinearCTStereoCarb;
    
                        /* save stereo canonical numbering */
                    if (pCS->nCanonOrdIsotopicStereo)
                    {
                        for (i = 0; i < num_at_tg; i++) /* djb-rwth: removing redundant code */
                        {
                            if (nCanonRankIsotopicStereo[i] && (int) nCanonRankIsotopicStereo[i] <= num_at_tg)
                            {
                                pCS->nCanonOrdIsotopicStereo[(int) nCanonRankIsotopicStereo[i] - 1] = (AT_NUMB) i;
                            }
                            else
                            {
                                bFailed++;
                            }
                        }
                        pCS->nLenCanonOrdIsotopicStereo = bFailed ? -num_atoms : num_atoms;
                    }
    
                    /* save stereo tautomer groups numbering */
                    if (pCS->nCanonOrdIsotopicStereoTaut)
                    {
                        if (0 < ( nRet = SortTautomerGroupsAndEndpoints( pCG, t_group_info1, num_atoms, num_at_tg, nCanonRankIsotopicStereo ) ))
                        {
                            /*non-isotopic contains symmetry ranks */
                            int num_t_groups = t_group_info1->num_t_groups;
                            AT_NUMB *tGroupNumber = t_group_info1->tGroupNumber;
                            /*AT_NUMB *tiSymmRank        = tGroupNumber + TGSO_SYMM_IRANK*num_t_groups; */
                            memcpy(pCS->nCanonOrdIsotopicStereoTaut, tGroupNumber, num_t_groups * sizeof(pCS->nCanonOrdIsotopicStereoTaut[0]));
                            pCS->nLenCanonOrdIsotopicStereoTaut = bFailed ? -num_t_groups : num_t_groups;
    
                            /*SortTautomerGroupsAndEndpoints( t_group_info1, nCanonRank ); */ /* ??? return to non-isotopic canonical numbering */
                        }
                        else
                        {
                            if (RETURNED_ERROR(nRet))
                            {
                                goto exit_function;
                            }
                        }
                        /* djb-rwth: removing redundant code */
                    }
                }
            }
    
            /**********************************************************
             *
             *  VII-B. Optimize INVERTED isotopic stereo descriptors
             *
             **********************************************************/
            if (!nCanonRankIsotopicStereoInv)
                nCanonRankIsotopicStereoInv = (AT_RANK *) qmalloc( num_max * sizeof( *nCanonRankIsotopicStereoInv ) );
            if (!nCanonRankIsotopicStereoInv)
            {
                nRet = CT_OUT_OF_RAM;  /*  <BRKPT> */
                goto exit_function;
            }
    
            /* copy previous isotopic stereo canonicalization results to Inv initial data */
            /* assign pointers */
            pCS->LinearCTStereoDble = pCS2->LinearCTIsotopicStereoDbleInv; /*  enable stereo */
            pCS->LinearCTStereoCarb = pCS2->LinearCTIsotopicStereoCarbInv;
    
            /* copy the lengths */
            pCS2->nLenLinearCTIsotopicStereoDbleInv =
                pCS->nLenLinearCTStereoDbleInv =
                pCS->nLenLinearCTStereoDble = pCS2->nLenLinearCTIsotopicStereoDble;
    
            pCS2->nLenLinearCTIsotopicStereoCarbInv =
                pCS->nLenLinearCTStereoCarbInv =
                pCS->nLenLinearCTStereoCarb = pCS2->nLenLinearCTIsotopicStereoCarb;
    
            if (pCS->nLenLinearCTStereoDble > 0 || pCS->nLenLinearCTStereoCarb > 0)
            {
                /* copy previous results, the canonical stereo CT */
                memcpy(pCS->LinearCTStereoDble, pCS2->LinearCTIsotopicStereoDble, pCS->nLenLinearCTStereoDble * sizeof(pCS->LinearCTStereoDble[0]));
                memcpy(pCS->LinearCTStereoCarb, pCS2->LinearCTIsotopicStereoCarb, pCS->nLenLinearCTStereoCarb * sizeof(pCS->LinearCTStereoCarb[0]));
            }
            memcpy(nCanonRankIsotopicStereoInv, nCanonRankIsotopicStereo, num_max * sizeof(nCanonRankIsotopicStereoInv[0]));
            if (pCS->nCanonOrdIsotopicStereoInv && pCS->nCanonOrdIsotopicStereo)
            {
                /* in case there is nothing to invert */
                memcpy(pCS->nCanonOrdIsotopicStereoInv, pCS->nCanonOrdIsotopicStereo, num_at_tg * sizeof(pCS->nCanonOrdIsotopicStereoInv[0]));
            }
    
            /******************************
             *
             * Invert isotopic stereo
             *
             ******************************/
    
            /*********************************************************************************
             * Create initial approximation for the minimization of the stereo descriptors:
             *  invert stereogenic atom parities, one parity in each allene, all parities in
             *  pCS->LinearCTStereoCarb and allene parities in pCS->nLenLinearCTStereoDble
             */
    
            nRet = InvertStereo( at,
                                 num_at_tg,
                                 nCanonRankIsotopicStereo,
                                 nTempRank,
                                 pCS,
                                 1 );
    
            if (RETURNED_ERROR( nRet ))
            {
                goto exit_function;
            }
            else if (nRet > 0)
            {
                /* InvertStereo() has done some changes */
                /* djb-rwth: removing redundant code */
                /* FillOutStereoParities() has already been called to fill out these 2 LinearCTs */
    
                /* set the 1st ranks in the rest of the stack to zero: prepare for ranks reuse */
                for (n = 2; n < nRankStackLen && pRankStack1[n]; n++)
                {
                    pRankStack1[n][0] = 0; /* means ranks have to be recalculated */
                }
                /* set the 1st ranks to zero: prepare for ranks reuse */
                for (n = 2; n < nRankStackLen && pRankStack2[n]; n++)
                {
                    pRankStack2[n][0] = 0; /* means ranks have to be recalculated */
                }
                /* for debugging or statistics */
                pCS->lNumBreakTies =
                    pCS->lNumNeighListIter =
                    pCS->lNumTotCT =
                    pCS->lNumDecreasedCT =
                    pCS->lNumRejectedCT =
                    pCS->lNumEqualCT = 0;
                pCS->bKeepSymmRank = 0;
                pCS->bFirstCT = 1; /* To fill out nCanonRankStereo[] in map_stero_atoms2() */
    
                /**************************************************************************************
                  nCanonRankIsotopic contains input canonical numbering
                  nCanonRankIsotopicStereo will be filled with a transposition of canonical numbering
                    that  (1) keeps connection table unchanged and
                          (2) provides minimal stereo descriptors in
                              pCS->LinearCTStereoDble (length=pCS->nLenLinearCTStereoDble)
                              pCS->LinearCTStereoCarb (length=pCS->nLenLinearCTStereoCarb)
                */
                nRet = map_stereo_bonds4( ic, pCG,
                                          at,
                                          num_atoms,
                                          num_at_tg,
                                          num_max,
                                          0,
                                          ftcn->PartitionCtIso.Rank,
                                          ftcn->PartitionCtIso.AtNumber,
                                          nCanonRankIsotopicStereoInv,
                                          nSymmRank,
                                          pRankStack1,
                                          pRankStack2,
                                          nTempRank,
                                          nNumCurrRanks,
                                          nSymmStereo,
                                          NeighList,
                                          pCS,
                                          cur_tree,
                                          0,
                                          vABParityUnknown );
    
                if (RETURNED_ERROR( nRet ))
                {
                    if (nRet == CT_TIMEOUT_ERR)
                        goto exit_function;
                    else
                        goto exit_function; /* program error */
                }
                else
                {
                    int bFailed = 0;
    
                    if (!nRet)
                    {
                        /* djb-rwth: removing redundant code */
                        pCS2->nLenLinearCTIsotopicStereoDble =
                            pCS->nLenLinearCTIsotopicStereoDble = -abs( pCS->nLenLinearCTStereoDble );
                        pCS2->nLenLinearCTIsotopicStereoCarb =
                            pCS->nLenLinearCTIsotopicStereoCarb = -abs( pCS->nLenLinearCTStereoCarb );
                        nRet = CT_STEREOCOUNT_ERR;  /*  <BRKPT> */
                        goto exit_function; /* program error */
                    }
    
                    /* save isotopic pointers & lengths for INVERTED stereo */
    
                    /* save isotopic lengths */
                    pCS->nLenLinearCTIsotopicStereoDbleInv =
                        pCS2->nLenLinearCTIsotopicStereoDbleInv = pCS->nLenLinearCTStereoDble;
                    pCS->nLenLinearCTIsotopicStereoCarbInv =
                        pCS2->nLenLinearCTIsotopicStereoCarbInv = pCS->nLenLinearCTStereoCarb;
    
                        /* restore pointers and lengths to non-inverted isotopic stereo */
                        /*  -- this is needed for InvertStereo() back, see below        */
                    pCS->LinearCTStereoDble = pCS2->LinearCTIsotopicStereoDble;
                    pCS->LinearCTStereoCarb = pCS2->LinearCTIsotopicStereoCarb;
                    pCS->nLenLinearCTStereoDble = pCS2->nLenLinearCTIsotopicStereoDble;
                    pCS->nLenLinearCTStereoCarb = pCS2->nLenLinearCTIsotopicStereoCarb;
    
                    /* consistency check */
                    if (pCS->nLenLinearCTIsotopicStereoDbleInv != pCS->nLenLinearCTIsotopicStereoDble ||
                         pCS->nLenLinearCTIsotopicStereoCarbInv != pCS->nLenLinearCTIsotopicStereoCarb)
                    {
                        nRet = CT_CALC_STEREO_ERR;
                        goto exit_function; /* program error */
                    }
    
                    /******************************
                     *
                     * Invert stereo back
                     *
                     ******************************
                     *  (make sure that pointers
                     *  pCS->LinearCTStereoCarb,
                     *  pCS->LinearCTStereoDble
                     *  and corresponding lengths
                     *  have been restored)
                     ******************************/
    
                    nRet = InvertStereo( at, num_at_tg, nCanonRankIsotopicStereo,
                                         nTempRank, pCS, 0 );
    
                    if (RETURNED_ERROR( nRet ))
                    {
                        goto exit_function;
                    }
                    nRet = 0;
    
                    /* save stereo canonical numbering */
                    if (pCS->nCanonOrdIsotopicStereoInv)
                    {
                        for (i = 0; i < num_at_tg; i++) /* djb-rwth: removing redundant code */
                        {
                            if (nCanonRankIsotopicStereoInv[i] && (int) nCanonRankIsotopicStereoInv[i] <= num_at_tg)
                            {
                                pCS->nCanonOrdIsotopicStereoInv[(int) nCanonRankIsotopicStereoInv[i] - 1] = (AT_NUMB) i;
                            }
                            else
                            {
                                bFailed++;
                            }
                        }
                        pCS->nLenCanonOrdIsotopicStereo = bFailed ? -num_atoms : num_atoms;
                    }
    
                    /* compare inverted and non-inverted isotopic stereo */
                    pCS->bCmpIsotopicStereo = 2 + CompareLinCtStereo(
                                             pCS->LinearCTIsotopicStereoDbleInv, pCS->nLenLinearCTIsotopicStereoDbleInv,
                                             pCS->LinearCTIsotopicStereoCarbInv, pCS->nLenLinearCTIsotopicStereoCarbInv,
                                             pCS->LinearCTIsotopicStereoDble, pCS->nLenLinearCTIsotopicStereoDble,
                                             pCS->LinearCTIsotopicStereoCarb, pCS->nLenLinearCTIsotopicStereoCarb
                    );
                }
            }
            else if (0 == nRet)
            {
                /* nothing has been done, restore pointers and lengths for stereo */
                pCS->LinearCTStereoDble = pCS2->LinearCTIsotopicStereoDble;
                pCS->LinearCTStereoCarb = pCS2->LinearCTIsotopicStereoCarb;
                pCS->nLenLinearCTStereoDble = pCS2->nLenLinearCTIsotopicStereoDble;
                pCS->nLenLinearCTStereoCarb = pCS2->nLenLinearCTIsotopicStereoCarb;
            }
    
    
        bypass_isotopic_stereo:;  /* ???       */
    
            pCS->LinearCTIsotopic = pCS2->LinearCTIsotopic;
        }
    
    
    
    exit_function:
    
        if (bSwitchedAtomToIsotopic)
        {
            SwitchAtomStereoAndIsotopicStereo( at, num_atoms, &bSwitchedAtomToIsotopic );
            SetCtToNonIsotopicStereo( pCS, pCS2 ); /* ??? */
        }
    
        /* restore non-isotopic connection table */
        if (pCS->LinearCT2)
        {
            inchi_swap( (char*) &pCS->LinearCT, (char*) &pCS->LinearCT2, sizeof( pCS->LinearCT ) );
            inchi_swap( (char*) &pCS->nLenLinearCT, (char*) &pCS->nLenLinearCT2, sizeof( pCS->nLenLinearCT ) );
            inchi_swap( (char*) &pCS->nLenLinearCTAtOnly, (char*) &pCS->nLenLinearCTAtOnly2, sizeof( pCS->nLenLinearCTAtOnly ) );
        }
    
        /* free memory */
        i = 2;
        if (pRankStack1)
        {
            pRankStack1[0] =
                pRankStack1[1] = NULL; /* deallocated separately */
            for (; i < nRankStackLen && pRankStack1[i]; i++)
            {
                ;
            }
        }
        if (pRankStack1 && pRankStack2)
        {
            for (n = 2; n < nRankStackLen && pRankStack2[n]; n++)
            {
                if (i < nRankStackLen - 1)
                {
                    pRankStack1[i++] = pRankStack2[n];
                }
                else
                {
                    inchi_free( pRankStack2[n] );
                }
            }
        }
    
        inchi_free(pRankStack2); /* djb-rwth: fixing coverity ID #499631 */
    
        pCS->NeighList = NULL; /* keep the pointer in pBCN->ftcn[bTaut].NeighList for further deallocation */
        qfree( nAtomNumber );
        qfree( nTempRank );
        qfree( nRank );
        qfree( nSymmRank );
    
        qfree( nSymmStereo );
        CurTreeFree( cur_tree );
        /* memory leak fix */
        /*
        qfree ( nCurrRankIsotopicStereo );
        qfree ( nAtomNumberCurrIsotopicStereo);
        */
        qfree( nCanonRankIsotopicStereo );
        qfree( nCanonRankIsotopicStereoInv );
    
        qfree( nCanonRankStereo );
        qfree( nCanonRankStereoInv );
    
        InchiTimeGet( &ulEndTime );
    
        pCS->lTotalTime = InchiTimeMsecDiff( ic, &ulEndTime, &ulStartTime );
    
        return ( nRet >= -1 ) ? num_atoms : nRet;
            /* cannot easily get number of ranks for now */
    }
    */

    let mut start_time = inchiTime::default();
    InchiTimeGet(&mut start_time, clock_result);
    let bcn = heap
        .slice(pCS.pBCN.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let ftcn_index = usize::try_from(bTautFtcn).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let ftcn = bcn
        .ftcn
        .get(ftcn_index)
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut tgi = if pCS.t_group_info.is_null() {
        None
    } else {
        Some(
            heap.slice(pCS.t_group_info.as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
    };
    let b_taut = num_at_tg > num_atoms
        && tgi
            .as_ref()
            .is_some_and(|value| value.num_t_groups != 0 && !value.t_group.is_null());
    let b_iso_xchg_h = if let Some(value) = tgi.as_ref() {
        value.nNumIsotopicEndpoints > 1
            && !value.nIsotopicEndpointAtomNumber.is_null()
            && heap
                .slice(value.nIsotopicEndpointAtomNumber.as_const())?
                .first()
                .is_some_and(|endpoint| *endpoint != 0)
            && value.bTautFlagsDone
                & u64::from(TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE)
                != 0
    } else {
        false
    };
    let has_isotopic_canon_data = !ftcn.PartitionCtIso.AtNumber.is_null()
        && !ftcn.PartitionCtIso.Rank.is_null()
        && !ftcn.nSymmRankCtIso.is_null();
    let mut canon_isotopic = has_isotopic_canon_data
        && nMode & CMODE_ISO != 0
        && (!pCS.LinearCTIsotopic.is_null()
            || !pCS.LinearCTIsotopicTautomer.is_null()
            || b_iso_xchg_h);
    let mut canon_stereo = nMode & CMODE_STEREO != 0
        && (!pCS.LinearCTStereoDble.is_null() || !pCS.LinearCTStereoCarb.is_null());
    let mut canon_iso_stereo = has_isotopic_canon_data
        && nMode & CMODE_ISO_STEREO != 0
        && (!pCS.LinearCTIsotopicStereoDble.is_null()
            || !pCS.LinearCTIsotopicStereoCarb.is_null())
        && canon_isotopic;
    let mut iso_taut = b_taut && canon_isotopic;
    let mut switched_atom_to_isotopic = 0_i32;
    let parity_unknown = if nMode & REQ_MODE_DIFF_UU_STEREO != 0 {
        AB_PARITY_UNKN as i32
    } else {
        AB_PARITY_UNDF as i32
    };
    let mut saved = pCS.clone();
    let mut current_tree = CUR_TREE::default();
    let mut use_tree = false;

    let mut atom_number = SourceMutPointer::<AT_RANK>::null();
    let mut rank = SourceMutPointer::<AT_RANK>::null();
    let mut temp_rank = SourceMutPointer::<AT_RANK>::null();
    let mut symm_rank = SourceMutPointer::<AT_RANK>::null();
    let mut rank_stack2 = SourceMutPointer::<SourceMutPointer<AT_RANK>>::null();
    let mut canon_rank_stereo = SourceMutPointer::<AT_RANK>::null();
    let mut canon_rank_stereo_inv = SourceMutPointer::<AT_RANK>::null();
    let mut symm_stereo = SourceMutPointer::<AT_RANK>::null();
    let mut canon_rank_iso_stereo = SourceMutPointer::<AT_RANK>::null();
    let mut canon_rank_iso_stereo_inv = SourceMutPointer::<AT_RANK>::null();

    let status = (|| -> Result<i32, SourceHeapError> {
        let status = 'exit_function: {
            let mut ignore_iso_groups = 1_i32;
            if b_taut {
                let value = tgi.as_mut().ok_or(SourceHeapError::NullPointer)?;
                ignore_iso_groups = value.bIgnoreIsotopic;
                value.bIgnoreIsotopic = 1;
                heap.slice_mut(pCS.t_group_info)?[0].bIgnoreIsotopic = 1;
            }
            let ignore_iso_atoms = pCS.bIgnoreIsotopic;
            pCS.bIgnoreIsotopic = 1;
            pCS.nCanonFlags = ftcn.nCanonFlags;
            let linear_count = usize::try_from(ftcn.nLenLinearCt)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            source_copy(heap, pCS.LinearCT, ftcn.LinearCt.as_const(), linear_count)?;
            let atom_count = usize::try_from(num_atoms)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let atom_tg_count = usize::try_from(num_at_tg)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            if !pCS.nNum_H.is_null() && !ftcn.nNumH.is_null() {
                let input = heap
                    .slice(ftcn.nNumH.as_const())?
                    .get(..atom_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                let output = heap
                    .slice_mut(pCS.nNum_H)?
                    .get_mut(..atom_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                for (output, input) in output.iter_mut().zip(input) {
                    *output = (input as u16 & 0xff) as u8 as i8;
                }
            }
            if !pCS.nNum_H_fixed.is_null() && !ftcn.nNumHFixH.is_null() {
                let input = heap
                    .slice(ftcn.nNumHFixH.as_const())?
                    .get(..atom_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                let output = heap
                    .slice_mut(pCS.nNum_H_fixed)?
                    .get_mut(..atom_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                for (output, input) in output.iter_mut().zip(input) {
                    *output = (input as u16 & 0xff) as u8 as i8;
                }
            }
            pCS.nLenLinearCT = ftcn.nLenLinearCt;
            pCS.nLenLinearCTAtOnly = ftcn.nLenLinearCtAtOnly;
            if !pCS.nSymmRank.is_null() {
                source_copy(heap, pCS.nSymmRank, ftcn.nSymmRankCt.as_const(), atom_tg_count)?;
            }
            if !pCS.nCanonOrd.is_null() {
                source_copy(
                    heap,
                    pCS.nCanonOrd,
                    ftcn.PartitionCt.AtNumber.as_const(),
                    atom_tg_count,
                )?;
                pCS.nLenCanonOrd = num_atoms;
            }
            if !ftcn.iso_exchg_atnos.is_null() && !pCS.nExchgIsoH.is_null() {
                let input = heap
                    .slice(ftcn.iso_exchg_atnos.as_const())?
                    .get(..atom_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                let output = heap
                    .slice_mut(pCS.nExchgIsoH)?
                    .get_mut(..atom_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                for (output, input) in output.iter_mut().zip(input) {
                    *output = i8::from(input == 0);
                }
            }
            let mut nret = 0_i32;
            if canon_isotopic {
                if !pCS.nSymmRankIsotopic.is_null() {
                    source_copy(
                        heap,
                        pCS.nSymmRankIsotopic,
                        ftcn.nSymmRankCtIso.as_const(),
                        atom_tg_count,
                    )?;
                }
                if !pCS.nCanonOrdIsotopic.is_null() {
                    source_copy(
                        heap,
                        pCS.nCanonOrdIsotopic,
                        ftcn.PartitionCtIso.AtNumber.as_const(),
                        atom_tg_count,
                    )?;
                    pCS.nLenCanonOrdIsotopic = num_at_tg;
                }
                nret = FillIsotopicAtLinearCT(
                    heap,
                    num_atoms,
                    at.as_const(),
                    ftcn.PartitionCtIso.AtNumber.as_const(),
                    pCS.LinearCTIsotopic,
                    pCS.nMaxLenLinearCTIsotopic,
                    &mut pCS.nLenLinearCTIsotopic,
                )?;
                if source_returned_error(nret) {
                    break 'exit_function nret;
                }
                if nret < 0 {
                    break 'exit_function CT_STEREOCOUNT_ERR.wrapping_add(5);
                }
            } else {
                pCS.nMaxLenLinearCTIsotopic = 0;
                pCS.nMaxLenLinearCTIsotopicTautomer = 0;
            }
            if b_taut {
                iso_taut = iso_taut
                    && !ftcn.PartitionCtIso.Rank.is_null()
                    && !ftcn.PartitionCtIso.AtNumber.is_null()
                    && !ftcn.nSymmRankCtIso.is_null();
                nret = FillTautLinearCT2(
                    heap,
                    pCG,
                    num_atoms,
                    num_at_tg,
                    i32::from(iso_taut),
                    ftcn.PartitionCt.Rank.as_const(),
                    ftcn.PartitionCt.AtNumber.as_const(),
                    ftcn.nSymmRankCt.as_const(),
                    ftcn.PartitionCtIso.Rank.as_const(),
                    ftcn.PartitionCtIso.AtNumber.as_const(),
                    ftcn.nSymmRankCtIso.as_const(),
                    pCS.LinearCTTautomer,
                    pCS.nMaxLenLinearCTTautomer,
                    &mut pCS.nLenLinearCTTautomer,
                    pCS.LinearCTIsotopicTautomer,
                    pCS.nMaxLenLinearCTIsotopicTautomer,
                    &mut pCS.nLenLinearCTIsotopicTautomer,
                    tgi.as_mut(),
                )?;
                if let Some(value) = tgi.as_ref() {
                    heap.slice_mut(pCS.t_group_info)?[0] = value.clone();
                }
                if source_returned_error(nret) {
                    break 'exit_function nret;
                }
                if nret <= 0 {
                    break 'exit_function CT_STEREOCOUNT_ERR.wrapping_add(5);
                }
                let value = tgi.as_ref().ok_or(SourceHeapError::NullPointer)?;
                let group_count = usize::try_from(value.num_t_groups)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if !pCS.nSymmRankTaut.is_null() {
                    source_copy(
                        heap,
                        pCS.nSymmRankTaut,
                        value.tGroupNumber.offset(value.num_t_groups as i64)?.as_const(),
                        group_count,
                    )?;
                }
                if !pCS.nCanonOrdTaut.is_null() {
                    source_copy(
                        heap,
                        pCS.nCanonOrdTaut,
                        value.tGroupNumber.as_const(),
                        group_count,
                    )?;
                    pCS.nLenCanonOrdTaut = value.num_t_groups;
                }
                if canon_isotopic {
                    if !pCS.nSymmRankIsotopicTaut.is_null() {
                        source_copy(
                            heap,
                            pCS.nSymmRankIsotopicTaut,
                            value
                                .tGroupNumber
                                .offset(i64::from(value.num_t_groups) * 3)?
                                .as_const(),
                            group_count,
                        )?;
                    }
                    source_copy(
                        heap,
                        pCS.nCanonOrdIsotopicTaut,
                        value
                            .tGroupNumber
                            .offset(i64::from(value.num_t_groups) * 2)?
                            .as_const(),
                        group_count,
                    )?;
                    pCS.nLenCanonOrdIsotopicTaut = value.num_t_groups;
                }
            }
            if !pCS.LinearCT2.is_null() {
                source_copy(
                    heap,
                    pCS.LinearCT2,
                    pCS.LinearCT.as_const(),
                    usize::try_from(pCS.nLenLinearCT)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )?;
                pCS.nLenLinearCT2 = pCS.nLenLinearCT;
                pCS.nLenLinearCTAtOnly2 = pCS.nLenLinearCTAtOnly;
            }
            if num_atoms <= 1 {
                canon_stereo = false;
                if num_atoms < 1
                    || heap
                        .slice(at.as_const())?
                        .first()
                        .is_none_or(|atom| atom.parity2 == 0)
                {
                    canon_iso_stereo = false;
                }
            }
            if !canon_stereo && !(canon_isotopic && canon_iso_stereo) {
                break 'exit_function nret;
            }

            nMode &= 0x00ff;
            let num_max = usize::try_from(bcn.num_max)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let allocate_rank = |heap: &mut SourceHeap| heap.allocate(vec![0_u16; num_max]);
            let atom_number_result = allocate_rank(heap);
            let rank_result = allocate_rank(heap);
            let temp_rank_result = allocate_rank(heap);
            let symm_rank_result = allocate_rank(heap);
            let source_malloc_result = |result| match result {
                Ok(pointer) => Ok(pointer),
                Err(SourceHeapError::AllocationFailed) => Ok(SourceMutPointer::null()),
                Err(error) => Err(error),
            };
            atom_number = source_malloc_result(atom_number_result)?;
            rank = source_malloc_result(rank_result)?;
            temp_rank = source_malloc_result(temp_rank_result)?;
            symm_rank = source_malloc_result(symm_rank_result)?;
            if ftcn.NeighList.is_null()
                || atom_number.is_null()
                || temp_rank.is_null()
                || rank.is_null()
                // The C condition omits nSymmRank, then dereferences it. Treat the
                // allocation-failure-only undefined path as out-of-memory so the
                // three successful sibling allocations still reach source cleanup.
                || symm_rank.is_null()
                || pCS.LinearCT.is_null()
            {
                break 'exit_function CT_OUT_OF_RAM;
            }
            pCS.NeighList = ftcn.NeighList;
            saved = pCS.clone();
            if nMode & CMODE_NOEQ_STEREO == 0 && (canon_stereo || canon_iso_stereo) {
                use_tree = true;
            }
            pCS.bCmpStereo = 0;
            pCS.bCmpIsotopicStereo = 0;
            if canon_stereo || canon_iso_stereo {
                let stack_len = usize::try_from(bcn.nMaxLenRankStack)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                rank_stack2 = match heap.allocate(vec![SourceMutPointer::null(); stack_len]) {
                    Ok(value) => value,
                    Err(SourceHeapError::AllocationFailed) => break 'exit_function CT_OUT_OF_RAM,
                    Err(error) => return Err(error),
                };
                let mut nn = 2_usize;
                while nn < stack_len {
                    let pointer = heap.slice(bcn.pRankStack.as_const())?[nn];
                    if pointer.is_null() {
                        break;
                    }
                    heap.slice_mut(pointer)?[0] = 0;
                    nn += 1;
                }
                let moved = nn.saturating_sub(2) / 2;
                if moved > 0 {
                    for index in 2 + moved..nn {
                        let pointer = heap.slice(bcn.pRankStack.as_const())?[index];
                        heap.slice_mut(rank_stack2)?[index - moved] = pointer;
                        heap.slice_mut(bcn.pRankStack)?[index] = SourceMutPointer::null();
                    }
                }
            }

            let mut num_current_ranks = 0_i32;
            if canon_stereo {
                source_copy(
                    heap,
                    atom_number,
                    ftcn.PartitionCt.AtNumber.as_const(),
                    atom_tg_count,
                )?;
                source_copy(
                    heap,
                    rank,
                    ftcn.PartitionCt.Rank.as_const(),
                    atom_tg_count,
                )?;
                source_copy(
                    heap,
                    symm_rank,
                    ftcn.nSymmRankCt.as_const(),
                    atom_tg_count,
                )?;
                num_current_ranks = FixCanonEquivalenceInfo(
                    heap,
                    pCG,
                    num_at_tg,
                    symm_rank,
                    rank,
                    temp_rank,
                    atom_number,
                    None,
                )?;
                source_copy(
                    heap,
                    pCS.nPrevAtomNumber,
                    ftcn.PartitionCt.AtNumber.as_const(),
                    atom_tg_count,
                )?;
                nret = FillOutStereoParities(
                    heap,
                    at,
                    num_atoms,
                    ftcn.PartitionCt.Rank,
                    ftcn.PartitionCt.AtNumber,
                    rank.as_const(),
                    atom_number.as_const(),
                    pCS,
                    pCG,
                    0,
                )?;
                if source_returned_error(nret) {
                    break 'exit_function nret;
                }
                if nret < 0 {
                    break 'exit_function CT_STEREOCOUNT_ERR;
                }
                canon_rank_stereo = match allocate_rank(heap) {
                    Ok(value) => value,
                    Err(SourceHeapError::AllocationFailed) => break 'exit_function CT_OUT_OF_RAM,
                    Err(error) => return Err(error),
                };
                if nMode & CMODE_NOEQ_STEREO == 0 {
                    symm_stereo = match heap.allocate(vec![0_u16; num_max + 1]) {
                        Ok(value) => value,
                        Err(SourceHeapError::AllocationFailed) => {
                            break 'exit_function CT_OUT_OF_RAM;
                        }
                        Err(error) => return Err(error),
                    };
                    if CurTreeAlloc(heap, Some(&mut current_tree), num_at_tg)? < 0 {
                        break 'exit_function CT_OUT_OF_RAM;
                    }
                }
                if rank_stack2.is_null()
                    || canon_rank_stereo.is_null()
                    || (symm_stereo.is_null() && nMode & CMODE_NOEQ_STEREO == 0)
                {
                    break 'exit_function CT_OUT_OF_RAM;
                }
                heap.slice_mut(bcn.pRankStack)?[0] = rank;
                heap.slice_mut(rank_stack2)?[0] = rank;
                heap.slice_mut(bcn.pRankStack)?[1] = atom_number;
                heap.slice_mut(rank_stack2)?[1] = atom_number;
                for stack in [bcn.pRankStack, rank_stack2] {
                    let mut index = 2_usize;
                    let stack_len = usize::try_from(bcn.nMaxLenRankStack)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    while index < stack_len {
                        let pointer = heap.slice(stack.as_const())?[index];
                        if pointer.is_null() {
                            break;
                        }
                        heap.slice_mut(pointer)?[0] = 0;
                        index += 1;
                    }
                }
                pCS.lNumBreakTies = 0;
                pCS.lNumNeighListIter = 0;
                pCS.lNumTotCT = 0;
                pCS.lNumDecreasedCT = 0;
                pCS.lNumRejectedCT = 0;
                pCS.lNumEqualCT = 0;
                pCS.bKeepSymmRank = 0;
                pCS.bFirstCT = 1;
                nret = map_stereo_bonds4(
                    heap,
                    ic,
                    clock_result,
                    user_action,
                    console_quit,
                    pCG,
                    at,
                    num_atoms,
                    num_at_tg,
                    bcn.num_max,
                    0,
                    ftcn.PartitionCt.Rank,
                    ftcn.PartitionCt.AtNumber,
                    canon_rank_stereo,
                    symm_rank,
                    bcn.pRankStack,
                    rank_stack2,
                    temp_rank,
                    num_current_ranks,
                    symm_stereo,
                    ftcn.NeighList,
                    pCS,
                    if use_tree { Some(&mut current_tree) } else { None },
                    0,
                    parity_unknown,
                )?;
                if source_returned_error(nret) {
                    break 'exit_function nret;
                }
                let mut failed = 0_i32;
                if nret == 0 {
                    pCS.nLenLinearCTStereoCarb = -pCS.nLenLinearCTStereoCarb.abs();
                    saved.nLenLinearCTStereoCarb = pCS.nLenLinearCTStereoCarb;
                    pCS.nLenLinearCTStereoDble = -pCS.nLenLinearCTStereoDble.abs();
                    saved.nLenLinearCTStereoDble = pCS.nLenLinearCTStereoDble;
                    break 'exit_function CT_STEREOCOUNT_ERR;
                }
                saved.nLenLinearCTStereoDble = pCS.nLenLinearCTStereoDble;
                saved.nLenLinearCTStereoCarb = pCS.nLenLinearCTStereoCarb;
                if !pCS.nCanonOrdStereo.is_null() {
                    let ranks = heap
                        .slice(canon_rank_stereo.as_const())?
                        .get(..atom_tg_count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .to_vec();
                    let order = heap.slice_mut(pCS.nCanonOrdStereo)?;
                    for (index, value) in ranks.into_iter().enumerate() {
                        if value != 0 && i32::from(value) <= num_at_tg {
                            order[usize::from(value) - 1] = index as AT_NUMB;
                        } else {
                            failed += 1;
                        }
                    }
                    pCS.nLenCanonOrdStereo = if failed != 0 { -num_atoms } else { num_atoms };
                }
                if b_taut && !pCS.nCanonOrdStereoTaut.is_null() {
                    nret = SortTautomerGroupsAndEndpoints(
                        heap,
                        pCG,
                        tgi.as_ref().ok_or(SourceHeapError::NullPointer)?,
                        num_atoms,
                        num_at_tg,
                        canon_rank_stereo,
                    )?;
                    if nret > 0 {
                        let value = tgi.as_ref().ok_or(SourceHeapError::NullPointer)?;
                        let count = usize::try_from(value.num_t_groups)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        source_copy(
                            heap,
                            pCS.nCanonOrdStereoTaut,
                            value.tGroupNumber.as_const(),
                            count,
                        )?;
                        pCS.nLenCanonOrdStereoTaut = if failed != 0 {
                            -value.num_t_groups
                        } else {
                            value.num_t_groups
                        };
                    } else if source_returned_error(nret) {
                        break 'exit_function nret;
                    }
                }

                canon_rank_stereo_inv = match allocate_rank(heap) {
                    Ok(value) => value,
                    Err(SourceHeapError::AllocationFailed) => break 'exit_function CT_OUT_OF_RAM,
                    Err(error) => return Err(error),
                };
                pCS.LinearCTStereoDble = saved.LinearCTStereoDbleInv;
                pCS.LinearCTStereoCarb = saved.LinearCTStereoCarbInv;
                pCS.nLenLinearCTStereoDble = saved.nLenLinearCTStereoDble;
                pCS.nLenLinearCTStereoDbleInv = saved.nLenLinearCTStereoDble;
                saved.nLenLinearCTStereoDbleInv = saved.nLenLinearCTStereoDble;
                pCS.nLenLinearCTStereoCarb = saved.nLenLinearCTStereoCarb;
                pCS.nLenLinearCTStereoCarbInv = saved.nLenLinearCTStereoCarb;
                saved.nLenLinearCTStereoCarbInv = saved.nLenLinearCTStereoCarb;
                if pCS.nLenLinearCTStereoDble > 0 {
                    source_copy(
                        heap,
                        pCS.LinearCTStereoDble,
                        saved.LinearCTStereoDble.as_const(),
                        usize::try_from(pCS.nLenLinearCTStereoDble)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )?;
                }
                if pCS.nLenLinearCTStereoCarb > 0 {
                    source_copy(
                        heap,
                        pCS.LinearCTStereoCarb,
                        saved.LinearCTStereoCarb.as_const(),
                        usize::try_from(pCS.nLenLinearCTStereoCarb)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )?;
                }
                source_copy(
                    heap,
                    canon_rank_stereo_inv,
                    canon_rank_stereo.as_const(),
                    num_max,
                )?;
                if !pCS.nCanonOrdStereoInv.is_null() && !pCS.nCanonOrdStereo.is_null() {
                    source_copy(
                        heap,
                        pCS.nCanonOrdStereoInv,
                        pCS.nCanonOrdStereo.as_const(),
                        atom_tg_count,
                    )?;
                }
                nret = InvertStereo(
                    heap,
                    at,
                    num_at_tg,
                    canon_rank_stereo.as_const(),
                    temp_rank,
                    pCS,
                    1,
                )?;
                if source_returned_error(nret) {
                    break 'exit_function nret;
                }
                if nret > 0 {
                    for stack in [bcn.pRankStack, rank_stack2] {
                        let mut index = 2_usize;
                        while index < usize::try_from(bcn.nMaxLenRankStack)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?
                        {
                            let pointer = heap.slice(stack.as_const())?[index];
                            if pointer.is_null() {
                                break;
                            }
                            heap.slice_mut(pointer)?[0] = 0;
                            index += 1;
                        }
                    }
                    pCS.lNumBreakTies = 0;
                    pCS.lNumNeighListIter = 0;
                    pCS.lNumTotCT = 0;
                    pCS.lNumDecreasedCT = 0;
                    pCS.lNumRejectedCT = 0;
                    pCS.lNumEqualCT = 0;
                    pCS.bKeepSymmRank = 0;
                    pCS.bFirstCT = 1;
                    CurTreeSetPos(if use_tree { Some(&mut current_tree) } else { None }, 0);
                    nret = map_stereo_bonds4(
                        heap,
                        ic,
                        clock_result,
                        user_action,
                        console_quit,
                        pCG,
                        at,
                        num_atoms,
                        num_at_tg,
                        bcn.num_max,
                        0,
                        ftcn.PartitionCt.Rank,
                        ftcn.PartitionCt.AtNumber,
                        canon_rank_stereo_inv,
                        symm_rank,
                        bcn.pRankStack,
                        rank_stack2,
                        temp_rank,
                        num_current_ranks,
                        symm_stereo,
                        ftcn.NeighList,
                        pCS,
                        if use_tree { Some(&mut current_tree) } else { None },
                        0,
                        parity_unknown,
                    )?;
                    if source_returned_error(nret) {
                        break 'exit_function nret;
                    }
                    failed = 0;
                    if nret == 0 {
                        pCS.nLenLinearCTStereoCarb = -pCS.nLenLinearCTStereoCarb.abs();
                        saved.nLenLinearCTStereoCarb = pCS.nLenLinearCTStereoCarb;
                        pCS.nLenLinearCTStereoDble = -pCS.nLenLinearCTStereoDble.abs();
                        saved.nLenLinearCTStereoDble = pCS.nLenLinearCTStereoDble;
                        break 'exit_function CT_STEREOCOUNT_ERR;
                    }
                    pCS.nLenLinearCTStereoDbleInv = pCS.nLenLinearCTStereoDble;
                    saved.nLenLinearCTStereoDbleInv = pCS.nLenLinearCTStereoDble;
                    pCS.nLenLinearCTStereoCarbInv = pCS.nLenLinearCTStereoCarb;
                    saved.nLenLinearCTStereoCarbInv = pCS.nLenLinearCTStereoCarb;
                    pCS.LinearCTStereoDble = saved.LinearCTStereoDble;
                    pCS.LinearCTStereoCarb = saved.LinearCTStereoCarb;
                    pCS.nLenLinearCTStereoDble = saved.nLenLinearCTStereoDble;
                    pCS.nLenLinearCTStereoCarb = saved.nLenLinearCTStereoCarb;
                    if pCS.nLenLinearCTStereoDbleInv != pCS.nLenLinearCTStereoDble
                        || pCS.nLenLinearCTStereoCarbInv != pCS.nLenLinearCTStereoCarb
                    {
                        break 'exit_function CT_CALC_STEREO_ERR;
                    }
                    nret = InvertStereo(
                        heap,
                        at,
                        num_at_tg,
                        canon_rank_stereo.as_const(),
                        temp_rank,
                        pCS,
                        0,
                    )?;
                    if source_returned_error(nret) {
                        break 'exit_function nret;
                    }
                    nret = 0;
                    if !pCS.nCanonOrdStereoInv.is_null() {
                        let ranks = heap
                            .slice(canon_rank_stereo_inv.as_const())?
                            .get(..atom_tg_count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .to_vec();
                        let order = heap.slice_mut(pCS.nCanonOrdStereoInv)?;
                        for (index, value) in ranks.into_iter().enumerate() {
                            if value != 0 && i32::from(value) <= num_at_tg {
                                order[usize::from(value) - 1] = index as AT_NUMB;
                            } else {
                                failed += 1;
                            }
                        }
                        pCS.nLenCanonOrdStereo = if failed != 0 { -num_atoms } else { num_atoms };
                    }
                    pCS.bCmpStereo = 2
                        + CompareLinCtStereo(
                            heap,
                            pCS.LinearCTStereoDbleInv.as_const(),
                            pCS.nLenLinearCTStereoDbleInv,
                            pCS.LinearCTStereoCarbInv.as_const(),
                            pCS.nLenLinearCTStereoCarbInv,
                            pCS.LinearCTStereoDble.as_const(),
                            pCS.nLenLinearCTStereoDble,
                            pCS.LinearCTStereoCarb.as_const(),
                            pCS.nLenLinearCTStereoCarb,
                        )?;
                } else {
                    pCS.LinearCTStereoDble = saved.LinearCTStereoDble;
                    pCS.LinearCTStereoCarb = saved.LinearCTStereoCarb;
                    pCS.nLenLinearCTStereoDble = saved.nLenLinearCTStereoDble;
                    pCS.nLenLinearCTStereoCarb = saved.nLenLinearCTStereoCarb;
                }
            }
            if b_taut {
                let value = tgi.as_mut().ok_or(SourceHeapError::NullPointer)?;
                value.bIgnoreIsotopic = ignore_iso_groups;
                heap.slice_mut(pCS.t_group_info)?[0].bIgnoreIsotopic = ignore_iso_groups;
            }
            pCS.bIgnoreIsotopic = ignore_iso_atoms;

            if canon_iso_stereo && canon_isotopic {
                source_copy(
                    heap,
                    atom_number,
                    ftcn.PartitionCtIso.AtNumber.as_const(),
                    atom_tg_count,
                )?;
                source_copy(
                    heap,
                    rank,
                    ftcn.PartitionCtIso.Rank.as_const(),
                    atom_tg_count,
                )?;
                source_copy(
                    heap,
                    symm_rank,
                    ftcn.nSymmRankCtIso.as_const(),
                    atom_tg_count,
                )?;
                num_current_ranks = FixCanonEquivalenceInfo(
                    heap,
                    pCG,
                    num_at_tg,
                    symm_rank,
                    rank,
                    temp_rank,
                    atom_number,
                    None,
                )?;
                source_copy(
                    heap,
                    pCS.nPrevAtomNumber,
                    ftcn.PartitionCtIso.AtNumber.as_const(),
                    atom_tg_count,
                )?;
                canon_rank_iso_stereo = match allocate_rank(heap) {
                    Ok(value) => value,
                    Err(SourceHeapError::AllocationFailed) => break 'exit_function CT_OUT_OF_RAM,
                    Err(error) => return Err(error),
                };
                if symm_stereo.is_null() && nMode & CMODE_NOEQ_STEREO == 0 {
                    symm_stereo = match heap.allocate(vec![0_u16; num_max + 1]) {
                        Ok(value) => value,
                        Err(SourceHeapError::AllocationFailed) => {
                            break 'exit_function CT_OUT_OF_RAM;
                        }
                        Err(error) => return Err(error),
                    };
                }
                if nMode & CMODE_NOEQ_STEREO == 0
                    && CurTreeAlloc(heap, Some(&mut current_tree), num_at_tg)? != 0
                {
                    break 'exit_function CT_OUT_OF_RAM;
                }
                if rank_stack2.is_null()
                    || canon_rank_iso_stereo.is_null()
                    || (symm_stereo.is_null() && nMode & CMODE_NOEQ_STEREO == 0)
                {
                    break 'exit_function CT_OUT_OF_RAM;
                }
                heap.slice_mut(bcn.pRankStack)?[0] = rank;
                heap.slice_mut(rank_stack2)?[0] = rank;
                heap.slice_mut(bcn.pRankStack)?[1] = atom_number;
                heap.slice_mut(rank_stack2)?[1] = atom_number;
                SwitchAtomStereoAndIsotopicStereo(
                    heap,
                    at,
                    num_atoms,
                    &mut switched_atom_to_isotopic,
                )?;
                SetCtToIsotopicStereo(pCS, &saved);
                nret = FillOutStereoParities(
                    heap,
                    at,
                    num_atoms,
                    ftcn.PartitionCtIso.Rank,
                    ftcn.PartitionCtIso.AtNumber,
                    rank.as_const(),
                    atom_number.as_const(),
                    pCS,
                    pCG,
                    1,
                )?;
                if source_returned_error(nret) {
                    break 'exit_function nret;
                }
                if nret == 0 {
                    pCS.nLenLinearCTIsotopicStereoDble = 0;
                    saved.nLenLinearCTIsotopicStereoDble = 0;
                    pCS.nLenLinearCTIsotopicStereoCarb = 0;
                    saved.nLenLinearCTIsotopicStereoCarb = 0;
                    pCS.nLenCanonOrdIsotopicStereo = 0;
                    pCS.nLenCanonOrdIsotopicStereoTaut = 0;
                    pCS.nLenLinearCTIsotopicStereoDbleInv = 0;
                    saved.nLenLinearCTIsotopicStereoDbleInv = 0;
                    pCS.nLenLinearCTIsotopicStereoCarbInv = 0;
                    saved.nLenLinearCTIsotopicStereoCarbInv = 0;
                } else {
                    for stack in [bcn.pRankStack, rank_stack2] {
                        let mut index = 2_usize;
                        while index < usize::try_from(bcn.nMaxLenRankStack)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?
                        {
                            let pointer = heap.slice(stack.as_const())?[index];
                            if pointer.is_null() {
                                break;
                            }
                            heap.slice_mut(pointer)?[0] = 0;
                            index += 1;
                        }
                    }
                    pCS.lNumBreakTies = 0;
                    pCS.lNumNeighListIter = 0;
                    pCS.lNumTotCT = 0;
                    pCS.lNumDecreasedCT = 0;
                    pCS.lNumRejectedCT = 0;
                    pCS.lNumEqualCT = 0;
                    pCS.bKeepSymmRank = 0;
                    pCS.bFirstCT = 1;
                    nret = map_stereo_bonds4(
                        heap,
                        ic,
                        clock_result,
                        user_action,
                        console_quit,
                        pCG,
                        at,
                        num_atoms,
                        num_at_tg,
                        bcn.num_max,
                        0,
                        ftcn.PartitionCtIso.Rank,
                        ftcn.PartitionCtIso.AtNumber,
                        canon_rank_iso_stereo,
                        symm_rank,
                        bcn.pRankStack,
                        rank_stack2,
                        temp_rank,
                        num_current_ranks,
                        symm_stereo,
                        ftcn.NeighList,
                        pCS,
                        if use_tree { Some(&mut current_tree) } else { None },
                        0,
                        parity_unknown,
                    )?;
                    if source_returned_error(nret) {
                        break 'exit_function nret;
                    }
                    let mut failed = 0_i32;
                    if nret == 0 {
                        pCS.nLenLinearCTIsotopicStereoDble =
                            -pCS.nLenLinearCTStereoDble.abs();
                        saved.nLenLinearCTIsotopicStereoDble =
                            pCS.nLenLinearCTIsotopicStereoDble;
                        pCS.nLenLinearCTIsotopicStereoCarb =
                            -pCS.nLenLinearCTStereoCarb.abs();
                        saved.nLenLinearCTIsotopicStereoCarb =
                            pCS.nLenLinearCTIsotopicStereoCarb;
                        break 'exit_function CT_STEREOCOUNT_ERR;
                    }
                    pCS.nLenLinearCTIsotopicStereoDble = pCS.nLenLinearCTStereoDble;
                    saved.nLenLinearCTIsotopicStereoDble = pCS.nLenLinearCTStereoDble;
                    pCS.nLenLinearCTIsotopicStereoCarb = pCS.nLenLinearCTStereoCarb;
                    saved.nLenLinearCTIsotopicStereoCarb = pCS.nLenLinearCTStereoCarb;
                    if !pCS.nCanonOrdIsotopicStereo.is_null() {
                        let ranks = heap
                            .slice(canon_rank_iso_stereo.as_const())?
                            .get(..atom_tg_count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .to_vec();
                        let order = heap.slice_mut(pCS.nCanonOrdIsotopicStereo)?;
                        for (index, value) in ranks.into_iter().enumerate() {
                            if value != 0 && i32::from(value) <= num_at_tg {
                                order[usize::from(value) - 1] = index as AT_NUMB;
                            } else {
                                failed += 1;
                            }
                        }
                        pCS.nLenCanonOrdIsotopicStereo =
                            if failed != 0 { -num_atoms } else { num_atoms };
                    }
                    if !pCS.nCanonOrdIsotopicStereoTaut.is_null() {
                        nret = SortTautomerGroupsAndEndpoints(
                            heap,
                            pCG,
                            tgi.as_ref().ok_or(SourceHeapError::NullPointer)?,
                            num_atoms,
                            num_at_tg,
                            canon_rank_iso_stereo,
                        )?;
                        if nret > 0 {
                            let value = tgi.as_ref().ok_or(SourceHeapError::NullPointer)?;
                            let count = usize::try_from(value.num_t_groups)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            source_copy(
                                heap,
                                pCS.nCanonOrdIsotopicStereoTaut,
                                value.tGroupNumber.as_const(),
                                count,
                            )?;
                            pCS.nLenCanonOrdIsotopicStereoTaut = if failed != 0 {
                                -value.num_t_groups
                            } else {
                                value.num_t_groups
                            };
                        } else if source_returned_error(nret) {
                            break 'exit_function nret;
                        }
                    }

                    canon_rank_iso_stereo_inv = match allocate_rank(heap) {
                        Ok(value) => value,
                        Err(SourceHeapError::AllocationFailed) => {
                            break 'exit_function CT_OUT_OF_RAM;
                        }
                        Err(error) => return Err(error),
                    };
                    pCS.LinearCTStereoDble = saved.LinearCTIsotopicStereoDbleInv;
                    pCS.LinearCTStereoCarb = saved.LinearCTIsotopicStereoCarbInv;
                    pCS.nLenLinearCTStereoDble = saved.nLenLinearCTIsotopicStereoDble;
                    pCS.nLenLinearCTStereoDbleInv = saved.nLenLinearCTIsotopicStereoDble;
                    pCS.nLenLinearCTIsotopicStereoDbleInv =
                        saved.nLenLinearCTIsotopicStereoDble;
                    saved.nLenLinearCTIsotopicStereoDbleInv =
                        saved.nLenLinearCTIsotopicStereoDble;
                    pCS.nLenLinearCTStereoCarb = saved.nLenLinearCTIsotopicStereoCarb;
                    pCS.nLenLinearCTStereoCarbInv = saved.nLenLinearCTIsotopicStereoCarb;
                    pCS.nLenLinearCTIsotopicStereoCarbInv =
                        saved.nLenLinearCTIsotopicStereoCarb;
                    saved.nLenLinearCTIsotopicStereoCarbInv =
                        saved.nLenLinearCTIsotopicStereoCarb;
                    if pCS.nLenLinearCTStereoDble > 0 {
                        source_copy(
                            heap,
                            pCS.LinearCTStereoDble,
                            saved.LinearCTIsotopicStereoDble.as_const(),
                            usize::try_from(pCS.nLenLinearCTStereoDble)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )?;
                    }
                    if pCS.nLenLinearCTStereoCarb > 0 {
                        source_copy(
                            heap,
                            pCS.LinearCTStereoCarb,
                            saved.LinearCTIsotopicStereoCarb.as_const(),
                            usize::try_from(pCS.nLenLinearCTStereoCarb)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )?;
                    }
                    source_copy(
                        heap,
                        canon_rank_iso_stereo_inv,
                        canon_rank_iso_stereo.as_const(),
                        num_max,
                    )?;
                    if !pCS.nCanonOrdIsotopicStereoInv.is_null()
                        && !pCS.nCanonOrdIsotopicStereo.is_null()
                    {
                        source_copy(
                            heap,
                            pCS.nCanonOrdIsotopicStereoInv,
                            pCS.nCanonOrdIsotopicStereo.as_const(),
                            atom_tg_count,
                        )?;
                    }
                    nret = InvertStereo(
                        heap,
                        at,
                        num_at_tg,
                        canon_rank_iso_stereo.as_const(),
                        temp_rank,
                        pCS,
                        1,
                    )?;
                    if source_returned_error(nret) {
                        break 'exit_function nret;
                    }
                    if nret > 0 {
                        for stack in [bcn.pRankStack, rank_stack2] {
                            let mut index = 2_usize;
                            while index < usize::try_from(bcn.nMaxLenRankStack)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?
                            {
                                let pointer = heap.slice(stack.as_const())?[index];
                                if pointer.is_null() {
                                    break;
                                }
                                heap.slice_mut(pointer)?[0] = 0;
                                index += 1;
                            }
                        }
                        pCS.lNumBreakTies = 0;
                        pCS.lNumNeighListIter = 0;
                        pCS.lNumTotCT = 0;
                        pCS.lNumDecreasedCT = 0;
                        pCS.lNumRejectedCT = 0;
                        pCS.lNumEqualCT = 0;
                        pCS.bKeepSymmRank = 0;
                        pCS.bFirstCT = 1;
                        nret = map_stereo_bonds4(
                            heap,
                            ic,
                            clock_result,
                            user_action,
                            console_quit,
                            pCG,
                            at,
                            num_atoms,
                            num_at_tg,
                            bcn.num_max,
                            0,
                            ftcn.PartitionCtIso.Rank,
                            ftcn.PartitionCtIso.AtNumber,
                            canon_rank_iso_stereo_inv,
                            symm_rank,
                            bcn.pRankStack,
                            rank_stack2,
                            temp_rank,
                            num_current_ranks,
                            symm_stereo,
                            ftcn.NeighList,
                            pCS,
                            if use_tree { Some(&mut current_tree) } else { None },
                            0,
                            parity_unknown,
                        )?;
                        if source_returned_error(nret) {
                            break 'exit_function nret;
                        }
                        failed = 0;
                        if nret == 0 {
                            pCS.nLenLinearCTIsotopicStereoDble =
                                -pCS.nLenLinearCTStereoDble.abs();
                            saved.nLenLinearCTIsotopicStereoDble =
                                pCS.nLenLinearCTIsotopicStereoDble;
                            pCS.nLenLinearCTIsotopicStereoCarb =
                                -pCS.nLenLinearCTStereoCarb.abs();
                            saved.nLenLinearCTIsotopicStereoCarb =
                                pCS.nLenLinearCTIsotopicStereoCarb;
                            break 'exit_function CT_STEREOCOUNT_ERR;
                        }
                        pCS.nLenLinearCTIsotopicStereoDbleInv = pCS.nLenLinearCTStereoDble;
                        saved.nLenLinearCTIsotopicStereoDbleInv = pCS.nLenLinearCTStereoDble;
                        pCS.nLenLinearCTIsotopicStereoCarbInv = pCS.nLenLinearCTStereoCarb;
                        saved.nLenLinearCTIsotopicStereoCarbInv = pCS.nLenLinearCTStereoCarb;
                        pCS.LinearCTStereoDble = saved.LinearCTIsotopicStereoDble;
                        pCS.LinearCTStereoCarb = saved.LinearCTIsotopicStereoCarb;
                        pCS.nLenLinearCTStereoDble = saved.nLenLinearCTIsotopicStereoDble;
                        pCS.nLenLinearCTStereoCarb = saved.nLenLinearCTIsotopicStereoCarb;
                        if pCS.nLenLinearCTIsotopicStereoDbleInv
                            != pCS.nLenLinearCTIsotopicStereoDble
                            || pCS.nLenLinearCTIsotopicStereoCarbInv
                                != pCS.nLenLinearCTIsotopicStereoCarb
                        {
                            break 'exit_function CT_CALC_STEREO_ERR;
                        }
                        nret = InvertStereo(
                            heap,
                            at,
                            num_at_tg,
                            canon_rank_iso_stereo.as_const(),
                            temp_rank,
                            pCS,
                            0,
                        )?;
                        if source_returned_error(nret) {
                            break 'exit_function nret;
                        }
                        nret = 0;
                        if !pCS.nCanonOrdIsotopicStereoInv.is_null() {
                            let ranks = heap
                                .slice(canon_rank_iso_stereo_inv.as_const())?
                                .get(..atom_tg_count)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .to_vec();
                            let order = heap.slice_mut(pCS.nCanonOrdIsotopicStereoInv)?;
                            for (index, value) in ranks.into_iter().enumerate() {
                                if value != 0 && i32::from(value) <= num_at_tg {
                                    order[usize::from(value) - 1] = index as AT_NUMB;
                                } else {
                                    failed += 1;
                                }
                            }
                            pCS.nLenCanonOrdIsotopicStereo =
                                if failed != 0 { -num_atoms } else { num_atoms };
                        }
                        pCS.bCmpIsotopicStereo = 2
                            + CompareLinCtStereo(
                                heap,
                                pCS.LinearCTIsotopicStereoDbleInv.as_const(),
                                pCS.nLenLinearCTIsotopicStereoDbleInv,
                                pCS.LinearCTIsotopicStereoCarbInv.as_const(),
                                pCS.nLenLinearCTIsotopicStereoCarbInv,
                                pCS.LinearCTIsotopicStereoDble.as_const(),
                                pCS.nLenLinearCTIsotopicStereoDble,
                                pCS.LinearCTIsotopicStereoCarb.as_const(),
                                pCS.nLenLinearCTIsotopicStereoCarb,
                            )?;
                    } else {
                        pCS.LinearCTStereoDble = saved.LinearCTIsotopicStereoDble;
                        pCS.LinearCTStereoCarb = saved.LinearCTIsotopicStereoCarb;
                        pCS.nLenLinearCTStereoDble = saved.nLenLinearCTIsotopicStereoDble;
                        pCS.nLenLinearCTStereoCarb = saved.nLenLinearCTIsotopicStereoCarb;
                    }
                }
                pCS.LinearCTIsotopic = saved.LinearCTIsotopic;
            }

            if b_taut {
                let value = tgi.as_mut().ok_or(SourceHeapError::NullPointer)?;
                value.bIgnoreIsotopic = ignore_iso_groups;
                heap.slice_mut(pCS.t_group_info)?[0].bIgnoreIsotopic = ignore_iso_groups;
            }
            pCS.bIgnoreIsotopic = ignore_iso_atoms;
            break 'exit_function nret;
        };
        Ok(status)
    })()?;

    if switched_atom_to_isotopic != 0 {
        SwitchAtomStereoAndIsotopicStereo(heap, at, num_atoms, &mut switched_atom_to_isotopic)?;
        SetCtToNonIsotopicStereo(pCS, &saved);
    }
    if !pCS.LinearCT2.is_null() {
        std::mem::swap(&mut pCS.LinearCT, &mut pCS.LinearCT2);
        std::mem::swap(&mut pCS.nLenLinearCT, &mut pCS.nLenLinearCT2);
        std::mem::swap(&mut pCS.nLenLinearCTAtOnly, &mut pCS.nLenLinearCTAtOnly2);
    }
    let stack_len = usize::try_from(bcn.nMaxLenRankStack)
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mut first_free = 2_usize;
    if !bcn.pRankStack.is_null() {
        if stack_len > 0 {
            heap.slice_mut(bcn.pRankStack)?[0] = SourceMutPointer::null();
        }
        if stack_len > 1 {
            heap.slice_mut(bcn.pRankStack)?[1] = SourceMutPointer::null();
        }
        while first_free < stack_len
            && !heap.slice(bcn.pRankStack.as_const())?[first_free].is_null()
        {
            first_free += 1;
        }
    }
    if !bcn.pRankStack.is_null() && !rank_stack2.is_null() {
        let mut index = 2_usize;
        while index < stack_len && !heap.slice(rank_stack2.as_const())?[index].is_null() {
            let pointer = heap.slice(rank_stack2.as_const())?[index];
            if first_free < stack_len.saturating_sub(1) {
                heap.slice_mut(bcn.pRankStack)?[first_free] = pointer;
                first_free += 1;
            } else {
                inchi_free(heap, pointer)?;
            }
            index += 1;
        }
    }
    inchi_free(heap, rank_stack2)?;
    pCS.NeighList = SourceMutPointer::null();
    for pointer in [atom_number, temp_rank, rank, symm_rank, symm_stereo] {
        inchi_free(heap, pointer)?;
    }
    CurTreeFree(heap, if use_tree { Some(&mut current_tree) } else { None })?;
    for pointer in [
        canon_rank_iso_stereo,
        canon_rank_iso_stereo_inv,
        canon_rank_stereo,
        canon_rank_stereo_inv,
    ] {
        inchi_free(heap, pointer)?;
    }
    let mut end_time = inchiTime::default();
    InchiTimeGet(&mut end_time, clock_result);
    pCS.lTotalTime = InchiTimeMsecDiff(ic, Some(&end_time), Some(&start_time));
    // END INCHI C FUNCTION: Canon_INChI3
    Ok(if status >= -1 { num_atoms } else { status })
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn Canon_INChI(
    heap: &mut SourceHeap,
    ic: &mut INCHI_CLOCK,
    clock_result: clock_t,
    user_action: Option<fn() -> i32>,
    console_quit: Option<fn() -> i32>,
    num_atoms: i32,
    num_at_tg: i32,
    at: SourceMutPointer<sp_ATOM>,
    pCS: &mut CANON_STAT,
    pCG: &mut CANON_GLOBALS,
    nMode: u32,
    bTautFtcn: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:2727 Canon_INChI
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
    int Canon_INChI( INCHI_CLOCK *ic,
                    int num_atoms,
                    int num_at_tg,
                    sp_ATOM* at,
                    CANON_STAT* pCS,
                    CANON_GLOBALS *pCG,
                    INCHI_MODE nMode,
                    int bTautFtcn )
    {
        if (pCS->pBCN && !pCS->NeighList)
        {
            return Canon_INChI3( ic, num_atoms, num_at_tg, at, pCS, pCG, nMode, bTautFtcn );
        }

        return CT_CANON_ERR;
    }
    */
    // END INCHI C FUNCTION: Canon_INChI

    if !pCS.pBCN.is_null() && pCS.NeighList.is_null() {
        return Canon_INChI3(
            heap,
            ic,
            clock_result,
            user_action,
            console_quit,
            num_atoms,
            num_at_tg,
            at,
            pCS,
            pCG,
            nMode,
            bTautFtcn,
        );
    }
    Ok(CT_CANON_ERR)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source::base::ichisort::CreateNeighList;
    use crate::source_types::{BCN, CMODE_REDNDNT_STEREO, CT_USER_QUIT_ERR, FTCN, Partition};

    #[test]
    fn source_port__ichicano__updatefulllinearct__line_1050() {
        fn atom(neighbors: &[u16]) -> sp_ATOM {
            let mut value = sp_ATOM {
                valence: neighbors.len() as i8,
                ..sp_ATOM::default()
            };
            value.neighbor[..neighbors.len()].copy_from_slice(neighbors);
            value
        }
        fn fixture(
            ct: Vec<u16>,
            max: i32,
            length: i32,
            atom_length: i32,
        ) -> (
            SourceHeap,
            SourceMutPointer<sp_ATOM>,
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_RANK>,
            CANON_STAT,
        ) {
            let mut heap = SourceHeap::default();
            let atoms = heap
                .allocate_model_storage(vec![atom(&[1, 2]), atom(&[0]), atom(&[0])])
                .unwrap();
            let ranks = heap.allocate_model_storage(vec![3_u16, 1, 2]).unwrap();
            let order = heap.allocate_model_storage(vec![1_u16, 2, 0]).unwrap();
            let linear = heap.allocate_model_storage(ct).unwrap();
            let canon = CANON_STAT {
                LinearCT: linear,
                nMaxLenLinearCT: max,
                nLenLinearCT: length,
                nLenLinearCTAtOnly: atom_length,
                ..CANON_STAT::default()
            };
            (heap, atoms, ranks, order, canon)
        }
        fn call(
            heap: &mut SourceHeap,
            atoms: SourceMutPointer<sp_ATOM>,
            ranks: SourceMutPointer<AT_RANK>,
            order: SourceMutPointer<AT_RANK>,
            canon: &mut CANON_STAT,
            first: i32,
        ) -> Result<i32, SourceHeapError> {
            UpdateFullLinearCT(
                heap,
                3,
                3,
                atoms.as_const(),
                ranks,
                order,
                canon,
                &mut CANON_GLOBALS::default(),
                first,
            )
        }

        let (mut heap, atoms, ranks, order, mut canon) = fixture(vec![99; 8], 8, 0, 5);
        assert_eq!(call(&mut heap, atoms, ranks, order, &mut canon, 1), Ok(-1));
        assert_eq!(canon.nLenLinearCT, 5);
        assert_eq!(
            heap.slice(canon.LinearCT.as_const()).unwrap(),
            &[1, 2, 3, 1, 2, 99, 99, 99]
        );
        assert_eq!(call(&mut heap, atoms, ranks, order, &mut canon, 0), Ok(0));

        let (mut heap, atoms, ranks, order, mut greater) = fixture(vec![0; 5], 5, 5, 5);
        assert_eq!(call(&mut heap, atoms, ranks, order, &mut greater, 0), Ok(1));
        assert_eq!(heap.slice(greater.LinearCT.as_const()).unwrap()[0], 0);

        let (mut heap, atoms, ranks, order, mut smaller) =
            fixture(vec![2, 99, 99, 99, 99], 5, 5, 5);
        assert_eq!(
            call(&mut heap, atoms, ranks, order, &mut smaller, 0),
            Ok(-1)
        );
        assert_eq!(
            heap.slice(smaller.LinearCT.as_const()).unwrap(),
            &[1, 2, 3, 1, 2]
        );

        let (mut heap, atoms, ranks, order, mut overflow) = fixture(vec![90; 5], 2, 0, 5);
        assert_eq!(
            call(&mut heap, atoms, ranks, order, &mut overflow, 1),
            Ok(CT_OVERFLOW)
        );
        assert_eq!(
            &heap.slice(overflow.LinearCT.as_const()).unwrap()[..3],
            &[1, 2, 90]
        );
        for (length, atom_length) in [(4, 5), (5, 4)] {
            let (mut heap, atoms, ranks, order, mut wrong) =
                fixture(vec![0; 5], 5, length, atom_length);
            assert_eq!(
                call(&mut heap, atoms, ranks, order, &mut wrong, 1),
                Ok(CT_LEN_MISMATCH)
            );
        }

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![atom(&[]), atom(&[])])
            .unwrap();
        let ranks = heap.allocate_model_storage(vec![1_u16, 2, 3]).unwrap();
        let order = heap.allocate_model_storage(vec![0_u16, 1, 2]).unwrap();
        let endpoints = heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
        let groups = heap
            .allocate_model_storage(vec![crate::source_types::T_GROUP {
                nNumEndpoints: 2,
                ..crate::source_types::T_GROUP::default()
            }])
            .unwrap();
        let group_info = heap
            .allocate_model_storage(vec![T_GROUP_INFO {
                t_group: groups,
                nEndpointAtomNumber: endpoints,
                num_t_groups: 1,
                nNumEndpoints: 2,
                ..T_GROUP_INFO::default()
            }])
            .unwrap();
        let linear = heap.allocate_model_storage(vec![0_u16; 5]).unwrap();
        let mut canon = CANON_STAT {
            LinearCT: linear,
            nMaxLenLinearCT: 5,
            nLenLinearCTAtOnly: 2,
            t_group_info: group_info,
            ..CANON_STAT::default()
        };
        assert_eq!(
            UpdateFullLinearCT(
                &mut heap,
                2,
                3,
                atoms.as_const(),
                ranks,
                order,
                &mut canon,
                &mut CANON_GLOBALS::default(),
                1,
            ),
            Ok(-1)
        );
        assert_eq!(heap.slice(linear.as_const()).unwrap(), &[1, 2, 3, 1, 2]);
        assert_eq!(heap.slice(endpoints.as_const()).unwrap(), &[0, 1]);
    }

    #[test]
    fn source_port__ichicano__getcanonlengths__line_418() {
        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![
                sp_ATOM {
                    valence: 1,
                    iso_sort_key: 9,
                    parity: 1,
                    stereo_bond_neighbor: [2, 0, 0],
                    ..sp_ATOM::default()
                },
                sp_ATOM {
                    valence: 1,
                    parity: 1,
                    stereo_bond_neighbor: [1, 0, 0],
                    ..sp_ATOM::default()
                },
                sp_ATOM {
                    iso_sort_key: -3,
                    parity: 1,
                    ..sp_ATOM::default()
                },
                sp_ATOM::default(),
            ])
            .unwrap();
        let groups = heap
            .allocate_model_storage(vec![
                crate::source_types::T_GROUP {
                    nNumEndpoints: 7,
                    ..crate::source_types::T_GROUP::default()
                },
                crate::source_types::T_GROUP {
                    nNumEndpoints: 9,
                    ..crate::source_types::T_GROUP::default()
                },
            ])
            .unwrap();
        let group_info = T_GROUP_INFO {
            t_group: groups,
            num_t_groups: 2,
            nNumIsotopicEndpoints: 6,
            ..T_GROUP_INFO::default()
        };

        let mut sizes = ATOM_SIZES {
            nLenLinearCTTautomer: 1,
            ..ATOM_SIZES::default()
        };
        assert_eq!(
            GetCanonLengths(&heap, 4, atoms.as_const(), &mut sizes, Some(&group_info)),
            Ok(0)
        );
        assert_eq!(sizes.nLenBonds, 1);
        assert_eq!(sizes.nLenCTAtOnly, 5);
        assert_eq!(sizes.nLenCT, 23);
        assert_eq!(sizes.nLenIsotopic, 2);
        assert_eq!(sizes.nLenLinearCTStereoDble, 1);
        assert_eq!(sizes.nLenLinearCTStereoCarb, 1);
        assert_eq!(sizes.nLenIsotopicEndpoints, 6);

        let mut no_taut_sizes = ATOM_SIZES::default();
        assert_eq!(
            GetCanonLengths(
                &heap,
                4,
                atoms.as_const(),
                &mut no_taut_sizes,
                Some(&group_info),
            ),
            Ok(0)
        );
        assert_eq!(no_taut_sizes.nLenCT, 5);
        assert_eq!(no_taut_sizes.nLenIsotopicEndpoints, 6);

        let mut preserved = ATOM_SIZES {
            nLenCT: 20,
            nLenBonds: 20,
            nLenIsotopic: 20,
            nLenCTAtOnly: 20,
            nLenLinearCTStereoDble: 20,
            nLenLinearCTStereoCarb: 20,
            nLenIsotopicEndpoints: 20,
            ..ATOM_SIZES::default()
        };
        let expected = preserved.clone();
        assert_eq!(
            GetCanonLengths(&heap, 0, SourceConstPointer::null(), &mut preserved, None,),
            Ok(0)
        );
        assert_eq!(preserved, expected);

        let mut empty = ATOM_SIZES::default();
        assert_eq!(
            GetCanonLengths(&heap, 0, SourceConstPointer::null(), &mut empty, None,),
            Ok(0)
        );
        assert_eq!(empty.nLenCT, 1);
    }

    #[test]
    fn source_port__ichicano__deallocatecs__line_491() {
        let mut heap = SourceHeap::default();
        macro_rules! allocation {
            () => {
                heap.allocate_model_storage(vec![Default::default()])
                    .unwrap()
            };
        }

        let neighbor_atoms = allocation!();
        let neighbor_list = heap.allocate_model_storage(vec![neighbor_atoms]).unwrap();
        let retained_time = allocation!();
        let retained_group_info = allocation!();
        let retained_bcn = allocation!();
        let mut state = CANON_STAT {
            LinearCT: allocation!(),
            nCanonOrd: allocation!(),
            nSymmRank: allocation!(),
            nNum_H: allocation!(),
            nNum_H_fixed: allocation!(),
            nExchgIsoH: allocation!(),
            LinearCTIsotopic: allocation!(),
            nSymmRankIsotopic: allocation!(),
            nCanonOrdIsotopic: allocation!(),
            LinearCTIsotopicTautomer: allocation!(),
            nCanonOrdIsotopicTaut: allocation!(),
            nSymmRankIsotopicTaut: allocation!(),
            LinearCTStereoDble: allocation!(),
            LinearCTStereoCarb: allocation!(),
            LinearCTStereoDbleInv: allocation!(),
            LinearCTStereoCarbInv: allocation!(),
            nCanonOrdStereo: allocation!(),
            nCanonOrdStereoInv: allocation!(),
            nCanonOrdStereoTaut: allocation!(),
            LinearCTIsotopicStereoDble: allocation!(),
            LinearCTIsotopicStereoCarb: allocation!(),
            LinearCTIsotopicStereoDbleInv: allocation!(),
            LinearCTIsotopicStereoCarbInv: allocation!(),
            bRankUsedForStereo: allocation!(),
            bAtomUsedForStereo: allocation!(),
            nCanonOrdIsotopicStereo: allocation!(),
            nCanonOrdIsotopicStereoInv: allocation!(),
            nCanonOrdIsotopicStereoTaut: allocation!(),
            LinearCTTautomer: allocation!(),
            nCanonOrdTaut: allocation!(),
            nSymmRankTaut: allocation!(),
            LinearCT2: allocation!(),
            nPrevAtomNumber: allocation!(),
            NeighList: neighbor_list,
            ulTimeOutTime: retained_time,
            t_group_info: retained_group_info,
            pBCN: retained_bcn,
            nMaxLenLinearCTStereoDble: 1,
            nLenLinearCTStereoDble: 2,
            nMaxLenLinearCTStereoCarb: 3,
            nLenLinearCTStereoCarb: 4,
            nMaxLenLinearCTIsotopicStereoDble: 5,
            nLenLinearCTIsotopicStereoDble: 6,
            nMaxLenLinearCTIsotopicStereoCarb: 7,
            nLenLinearCTIsotopicStereoCarb: 8,
            nMaxLenLinearCTTautomer: 9,
            nLenLinearCTTautomer: 10,
            nMaxLenLinearCTIsotopic: 11,
            nLenLinearCTIsotopic: 12,
            nMaxLenLinearCTIsotopicTautomer: 13,
            nLenLinearCTIsotopicTautomer: 14,
            nLenCanonOrd: 15,
            nLenCanonOrdIsotopic: 16,
            nLenCanonOrdIsotopicTaut: 17,
            nLenCanonOrdStereo: 18,
            nLenCanonOrdStereoTaut: 19,
            nLenCanonOrdIsotopicStereo: 20,
            nLenCanonOrdIsotopicStereoTaut: 21,
            nLenCanonOrdTaut: 22,
            nLenLinearCT: 77,
            nLenLinearCTStereoDbleInv: 78,
            bFirstCT: 79,
            ..CANON_STAT::default()
        };
        assert_eq!(heap.live_allocation_count(), 38);

        assert_eq!(DeAllocateCS(&mut heap, &mut state), Ok(0));
        let expected = CANON_STAT {
            ulTimeOutTime: retained_time,
            t_group_info: retained_group_info,
            pBCN: retained_bcn,
            nLenLinearCT: 77,
            nLenLinearCTStereoDbleInv: 78,
            bFirstCT: 79,
            ..CANON_STAT::default()
        };
        assert_eq!(state, expected);
        assert_eq!(heap.live_allocation_count(), 3);

        assert_eq!(DeAllocateCS(&mut heap, &mut state), Ok(0));
        assert_eq!(state, expected);
        assert_eq!(heap.live_allocation_count(), 3);
    }

    #[test]
    fn source_port__ichicano__allocatecs__line_575() {
        const MODE: u64 = crate::source_types::CMODE_CT as u64
            | CMODE_ISO as u64
            | CMODE_STEREO as u64
            | CMODE_ISO_STEREO as u64
            | crate::source_types::CMODE_TAUT as u64;

        fn allocate_all(
            heap: &mut SourceHeap,
            state: &mut CANON_STAT,
            bcn: SourceMutPointer<BCN>,
        ) -> Result<i32, SourceHeapError> {
            AllocateCS(heap, state, 3, 5, 7, 6, 2, 3, 4, 5, 6, 7, 8, MODE, bcn)
        }

        let mut heap = SourceHeap::default();
        let bcn = heap.allocate_model_storage(vec![BCN::default()]).unwrap();
        heap.trace_source_allocations();
        let mut state = CANON_STAT::default();
        assert_eq!(allocate_all(&mut heap, &mut state, bcn), Ok(0));
        assert_eq!(heap.source_allocation_calls(), 33);
        assert_eq!(heap.live_allocation_count(), 34);
        macro_rules! assert_allocated {
            ($field:ident) => {
                assert!(!state.$field.is_null(), stringify!($field));
            };
        }
        assert_allocated!(LinearCT);
        assert_allocated!(nCanonOrd);
        assert_allocated!(nSymmRank);
        assert_allocated!(nNum_H);
        assert_allocated!(nNum_H_fixed);
        assert_allocated!(nExchgIsoH);
        assert_allocated!(LinearCTIsotopic);
        assert_allocated!(LinearCTIsotopicTautomer);
        assert_allocated!(nCanonOrdIsotopicTaut);
        assert_allocated!(nSymmRankIsotopicTaut);
        assert_allocated!(nSymmRankIsotopic);
        assert_allocated!(nCanonOrdIsotopic);
        assert_allocated!(LinearCTStereoDble);
        assert_allocated!(LinearCTStereoDbleInv);
        assert_allocated!(LinearCTStereoCarb);
        assert_allocated!(LinearCTStereoCarbInv);
        assert_allocated!(nCanonOrdStereo);
        assert_allocated!(nCanonOrdStereoInv);
        assert_allocated!(nCanonOrdStereoTaut);
        assert_allocated!(LinearCTIsotopicStereoDble);
        assert_allocated!(LinearCTIsotopicStereoDbleInv);
        assert_allocated!(LinearCTIsotopicStereoCarb);
        assert_allocated!(LinearCTIsotopicStereoCarbInv);
        assert_allocated!(nCanonOrdIsotopicStereo);
        assert_allocated!(nCanonOrdIsotopicStereoInv);
        assert_allocated!(nCanonOrdIsotopicStereoTaut);
        assert_allocated!(bRankUsedForStereo);
        assert_allocated!(bAtomUsedForStereo);
        assert_allocated!(LinearCTTautomer);
        assert_allocated!(nCanonOrdTaut);
        assert_allocated!(nSymmRankTaut);
        assert_allocated!(LinearCT2);
        assert_allocated!(nPrevAtomNumber);
        assert_eq!(heap.slice(state.LinearCT.as_const()).unwrap(), &[0; 7]);
        assert_eq!(heap.slice(state.nCanonOrd.as_const()).unwrap(), &[0; 5]);
        assert_eq!(heap.slice(state.nNum_H.as_const()).unwrap(), &[0; 3]);
        let isotopic = heap.slice(state.LinearCTIsotopic.as_const()).unwrap();
        assert_eq!(isotopic.len(), 8);
        assert!(
            isotopic
                .iter()
                .all(|value| value == &AT_ISOTOPIC::default())
        );
        let isotopic_tautomer = heap
            .slice(state.LinearCTIsotopicTautomer.as_const())
            .unwrap();
        assert_eq!(isotopic_tautomer.len(), 7);
        assert!(
            isotopic_tautomer
                .iter()
                .all(|value| value == &AT_ISO_TGROUP::default())
        );
        let stereo_double = heap.slice(state.LinearCTStereoDble.as_const()).unwrap();
        assert_eq!(stereo_double.len(), 2);
        assert!(
            stereo_double
                .iter()
                .all(|value| value == &crate::source_types::AT_STEREO_DBLE::default())
        );
        let stereo_atom = heap.slice(state.LinearCTStereoCarb.as_const()).unwrap();
        assert_eq!(stereo_atom.len(), 4);
        assert!(
            stereo_atom
                .iter()
                .all(|value| value == &crate::source_types::AT_STEREO_CARB::default())
        );
        let isotopic_stereo_double = heap
            .slice(state.LinearCTIsotopicStereoDble.as_const())
            .unwrap();
        assert_eq!(isotopic_stereo_double.len(), 3);
        assert!(
            isotopic_stereo_double
                .iter()
                .all(|value| value == &crate::source_types::AT_STEREO_DBLE::default())
        );
        let isotopic_stereo_atom = heap
            .slice(state.LinearCTIsotopicStereoCarb.as_const())
            .unwrap();
        assert_eq!(isotopic_stereo_atom.len(), 5);
        assert!(
            isotopic_stereo_atom
                .iter()
                .all(|value| value == &crate::source_types::AT_STEREO_CARB::default())
        );
        assert_eq!(
            heap.slice(state.LinearCTTautomer.as_const()).unwrap(),
            &[0; 6]
        );
        assert_eq!(
            heap.slice(state.nPrevAtomNumber.as_const()).unwrap(),
            &[0; 5]
        );
        assert_eq!(state.nMode, MODE);
        assert_eq!(state.nMaxLenLinearCT, 7);
        assert_eq!(state.nLenLinearCT, 7);
        assert_eq!(state.nLenLinearCTAtOnly, 6);
        assert_eq!(state.nMaxLenLinearCTIsotopic, 8);
        assert_eq!(state.nLenLinearCTIsotopic, 8);
        assert_eq!(state.nMaxLenLinearCTIsotopicTautomer, 7);
        assert_eq!(state.nLenLinearCTIsotopicTautomer, 7);
        assert_eq!(state.nLenLinearCTStereoDbleInv, 2);
        assert_eq!(state.nMaxLenLinearCTStereoDble, 2);
        assert_eq!(state.nLenLinearCTStereoDble, 2);
        assert_eq!(state.nLenLinearCTStereoCarbInv, 4);
        assert_eq!(state.nMaxLenLinearCTStereoCarb, 4);
        assert_eq!(state.nLenLinearCTStereoCarb, 4);
        assert_eq!(state.nLenLinearCTIsotopicStereoDbleInv, 3);
        assert_eq!(state.nMaxLenLinearCTIsotopicStereoDble, 3);
        assert_eq!(state.nLenLinearCTIsotopicStereoDble, 3);
        assert_eq!(state.nLenLinearCTIsotopicStereoCarbInv, 5);
        assert_eq!(state.nMaxLenLinearCTIsotopicStereoCarb, 5);
        assert_eq!(state.nLenLinearCTIsotopicStereoCarb, 5);
        assert_eq!(state.nMaxLenLinearCTTautomer, 6);
        assert_eq!(state.nLenLinearCTTautomer, 6);
        assert_eq!(DeAllocateCS(&mut heap, &mut state), Ok(0));
        assert_eq!(heap.live_allocation_count(), 1);

        let failed_state = CANON_STAT {
            nMode: MODE,
            nMaxLenLinearCT: 7,
            nLenLinearCT: 7,
            nLenLinearCTAtOnly: 6,
            nLenLinearCTStereoDbleInv: 2,
            nLenLinearCTStereoCarbInv: 4,
            nLenLinearCTIsotopicStereoDbleInv: 3,
            nLenLinearCTIsotopicStereoCarbInv: 5,
            ..CANON_STAT::default()
        };
        for failure_ordinal in 0..33_u64 {
            let mut failure_heap = SourceHeap::default();
            let failure_bcn = failure_heap
                .allocate_model_storage(vec![BCN::default()])
                .unwrap();
            failure_heap.fail_after_allocations(failure_ordinal);
            let mut failure_state = CANON_STAT::default();
            assert_eq!(
                allocate_all(&mut failure_heap, &mut failure_state, failure_bcn),
                Ok(CT_OUT_OF_RAM),
                "allocation ordinal {failure_ordinal}"
            );
            assert_eq!(
                failure_heap.source_allocation_calls(),
                33,
                "allocation ordinal {failure_ordinal} did not continue"
            );
            assert_eq!(
                failure_heap.live_allocation_count(),
                1,
                "allocation ordinal {failure_ordinal} leaked"
            );
            assert_eq!(
                failure_state, failed_state,
                "allocation ordinal {failure_ordinal} cleanup state"
            );
        }

        let mut minimal_heap = SourceHeap::default();
        minimal_heap.trace_source_allocations();
        let mut minimal_state = CANON_STAT::default();
        assert_eq!(
            AllocateCS(
                &mut minimal_heap,
                &mut minimal_state,
                3,
                5,
                7,
                6,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                0,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(minimal_heap.source_allocation_calls(), 1);
        assert_eq!(minimal_heap.live_allocation_count(), 1);
        assert!(!minimal_state.nPrevAtomNumber.is_null());
        let minimal_prev = minimal_state.nPrevAtomNumber;
        minimal_state.nPrevAtomNumber = SourceMutPointer::null();
        assert_eq!(minimal_state, CANON_STAT::default());
        minimal_state.nPrevAtomNumber = minimal_prev;
        assert_eq!(DeAllocateCS(&mut minimal_heap, &mut minimal_state), Ok(0));

        let mut negative_heap = SourceHeap::default();
        negative_heap.trace_source_allocations();
        let mut negative_state = CANON_STAT::default();
        assert_eq!(
            AllocateCS(
                &mut negative_heap,
                &mut negative_state,
                0,
                -1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                SourceMutPointer::null(),
            ),
            Ok(CT_OUT_OF_RAM)
        );
        assert_eq!(negative_heap.source_allocation_calls(), 0);
        assert_eq!(negative_heap.live_allocation_count(), 0);
        assert_eq!(negative_state, CANON_STAT::default());
    }

    #[test]
    fn source_port__ichicano__canon_inchi__line_2727() {
        let mut empty_heap = SourceHeap::default();
        assert_eq!(
            Canon_INChI(
                &mut empty_heap,
                &mut INCHI_CLOCK::default(),
                0,
                None,
                None,
                0,
                0,
                SourceMutPointer::null(),
                &mut CANON_STAT::default(),
                &mut CANON_GLOBALS::default(),
                0,
                0,
            ),
            Ok(CT_CANON_ERR)
        );

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![sp_ATOM::default()])
            .unwrap();
        let linear_input = heap.allocate_model_storage(vec![1_u16]).unwrap();
        let atom_order = heap.allocate_model_storage(vec![0_u16]).unwrap();
        let rank = heap.allocate_model_storage(vec![1_u16]).unwrap();
        let symmetry = heap.allocate_model_storage(vec![1_u16]).unwrap();
        let rank_stack = heap
            .allocate_model_storage(vec![
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ])
            .unwrap();
        let neighbors = CreateNeighList(&mut heap, 1, 1, atoms.as_const(), 0, None).unwrap();
        let bcn = heap
            .allocate_model_storage(vec![BCN {
                pRankStack: rank_stack,
                nMaxLenRankStack: 3,
                num_max: 1,
                num_at_tg: 1,
                num_atoms: 1,
                ftcn: [
                    FTCN {
                        LinearCt: linear_input,
                        nLenLinearCt: 1,
                        nLenLinearCtAtOnly: 1,
                        PartitionCt: Partition {
                            AtNumber: atom_order,
                            Rank: rank,
                        },
                        nSymmRankCt: symmetry,
                        NeighList: neighbors,
                        ..FTCN::default()
                    },
                    FTCN::default(),
                ],
                ..BCN::default()
            }])
            .unwrap();
        let output = heap.allocate_model_storage(vec![9_u16]).unwrap();
        let mut state = CANON_STAT {
            pBCN: bcn,
            LinearCT: output,
            ..CANON_STAT::default()
        };
        assert_eq!(
            Canon_INChI(
                &mut heap,
                &mut INCHI_CLOCK::default(),
                0,
                None,
                None,
                1,
                1,
                atoms,
                &mut state,
                &mut CANON_GLOBALS::default(),
                0,
                0,
            ),
            Ok(1)
        );
        assert_eq!(heap.slice(output.as_const()).unwrap(), &[1]);
        assert!(state.NeighList.is_null());

        state.NeighList = neighbors;
        heap.slice_mut(output).unwrap()[0] = 7;
        assert_eq!(
            Canon_INChI(
                &mut heap,
                &mut INCHI_CLOCK::default(),
                0,
                None,
                None,
                1,
                1,
                atoms,
                &mut state,
                &mut CANON_GLOBALS::default(),
                0,
                0,
            ),
            Ok(CT_CANON_ERR)
        );
        assert_eq!(heap.slice(output.as_const()).unwrap(), &[7]);
    }

    #[test]
    fn source_port__ichicano__canon_inchi3__line_1370() {
        fn request_quit() -> i32 {
            1
        }

        fn fixture(
            heap: &mut SourceHeap,
            stereo: bool,
        ) -> (
            SourceMutPointer<sp_ATOM>,
            CANON_STAT,
            SourceMutPointer<SourceMutPointer<AT_RANK>>,
        ) {
            let atoms = heap
                .allocate_model_storage(vec![sp_ATOM::default(); 2])
                .unwrap();
            let linear = heap.allocate_model_storage(vec![2_u16, 1]).unwrap();
            let atom_number = heap.allocate_model_storage(vec![0_u16, 1]).unwrap();
            let ranks = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
            let symmetry = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
            let hydrogens = heap.allocate_model_storage(vec![257_i16, -1]).unwrap();
            let fixed_hydrogens = heap.allocate_model_storage(vec![258_i16, 511]).unwrap();
            let scratch_rank = heap.allocate_model_storage(vec![7_u16; 2]).unwrap();
            let scratch_order = heap.allocate_model_storage(vec![8_u16; 2]).unwrap();
            let rank_stack = heap
                .allocate_model_storage(vec![
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    scratch_rank,
                    scratch_order,
                    SourceMutPointer::null(),
                ])
                .unwrap();
            let neighbors = CreateNeighList(&mut *heap, 2, 2, atoms.as_const(), 0, None).unwrap();
            let ftcn = FTCN {
                LinearCt: linear,
                nLenLinearCt: 2,
                nLenLinearCtAtOnly: 1,
                PartitionCt: Partition {
                    AtNumber: atom_number,
                    Rank: ranks,
                },
                nSymmRankCt: symmetry,
                nNumH: hydrogens,
                nNumHFixH: fixed_hydrogens,
                NeighList: neighbors,
                ..FTCN::default()
            };
            let bcn = heap
                .allocate_model_storage(vec![BCN {
                    pRankStack: rank_stack,
                    nMaxLenRankStack: 5,
                    num_max: 2,
                    num_at_tg: 2,
                    num_atoms: 2,
                    ftcn: [ftcn, FTCN::default()],
                    ..BCN::default()
                }])
                .unwrap();
            let output = heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
            let output2 = heap.allocate_model_storage(vec![9_u16; 2]).unwrap();
            let output_h = heap.allocate_model_storage(vec![0_i8; 2]).unwrap();
            let output_fixed_h = heap.allocate_model_storage(vec![0_i8; 2]).unwrap();
            let output_symmetry = heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
            let output_order = heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
            let stereo_double = if stereo {
                heap.allocate_model_storage(vec![Default::default(); 2])
                    .unwrap()
            } else {
                SourceMutPointer::null()
            };
            let stereo_double_inv = if stereo {
                heap.allocate_model_storage(vec![Default::default(); 2])
                    .unwrap()
            } else {
                SourceMutPointer::null()
            };
            let stereo_carb = if stereo {
                heap.allocate_model_storage(vec![Default::default(); 2])
                    .unwrap()
            } else {
                SourceMutPointer::null()
            };
            let stereo_carb_inv = if stereo {
                heap.allocate_model_storage(vec![Default::default(); 2])
                    .unwrap()
            } else {
                SourceMutPointer::null()
            };
            let ranks_used = heap.allocate_model_storage(vec![0_i8; 2]).unwrap();
            let atoms_used = heap.allocate_model_storage(vec![0_i8; 2]).unwrap();
            let previous_atoms = heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
            (
                atoms,
                CANON_STAT {
                    pBCN: bcn,
                    LinearCT: output,
                    LinearCT2: output2,
                    nNum_H: output_h,
                    nNum_H_fixed: output_fixed_h,
                    nSymmRank: output_symmetry,
                    nCanonOrd: output_order,
                    LinearCTStereoDble: stereo_double,
                    LinearCTStereoDbleInv: stereo_double_inv,
                    LinearCTStereoCarb: stereo_carb,
                    LinearCTStereoCarbInv: stereo_carb_inv,
                    nMaxLenLinearCTStereoDble: 2,
                    nMaxLenLinearCTStereoCarb: 2,
                    bRankUsedForStereo: ranks_used,
                    bAtomUsedForStereo: atoms_used,
                    nPrevAtomNumber: previous_atoms,
                    ..CANON_STAT::default()
                },
                rank_stack,
            )
        }

        fn enable_isotopic_stereo(
            heap: &mut SourceHeap,
            atoms: SourceMutPointer<sp_ATOM>,
            state: &mut CANON_STAT,
            parity2: i8,
        ) {
            let atom = &mut heap.slice_mut(atoms).unwrap()[0];
            atom.parity = 7;
            atom.parity2 = parity2;
            atom.final_parity = 8;
            atom.final_parity2 = 2;
            atom.stereo_atom_parity = 9;
            atom.stereo_atom_parity2 = parity2;
            atom.bHasStereoOrEquToStereo = 10;
            atom.bHasStereoOrEquToStereo2 = i8::from(parity2 != 0);
            let isotope_order = heap.allocate_model_storage(vec![0_u16, 1]).unwrap();
            let isotope_rank = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
            let isotope_symmetry = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
            heap.slice_mut(state.pBCN).unwrap()[0].ftcn[0].PartitionCtIso = Partition {
                AtNumber: isotope_order,
                Rank: isotope_rank,
            };
            heap.slice_mut(state.pBCN).unwrap()[0].ftcn[0].nSymmRankCtIso = isotope_symmetry;
            state.LinearCTIsotopic = heap
                .allocate_model_storage(vec![AT_ISOTOPIC::default(); 2])
                .unwrap();
            state.nMaxLenLinearCTIsotopic = 2;
            state.nSymmRankIsotopic = heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
            state.nCanonOrdIsotopic = heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
            state.LinearCTIsotopicStereoDble = heap
                .allocate_model_storage(vec![Default::default(); 2])
                .unwrap();
            state.LinearCTIsotopicStereoDbleInv = heap
                .allocate_model_storage(vec![Default::default(); 2])
                .unwrap();
            state.LinearCTIsotopicStereoCarb = heap
                .allocate_model_storage(vec![Default::default(); 2])
                .unwrap();
            state.LinearCTIsotopicStereoCarbInv = heap
                .allocate_model_storage(vec![Default::default(); 2])
                .unwrap();
            state.nMaxLenLinearCTIsotopicStereoDble = 2;
            state.nMaxLenLinearCTIsotopicStereoCarb = 2;
        }

        fn double_bond_fixture(
            heap: &mut SourceHeap,
            z_product: i8,
            with_center: bool,
        ) -> (SourceMutPointer<sp_ATOM>, CANON_STAT) {
            let mut atoms = vec![sp_ATOM::default(); 6];
            atoms[0].valence = 1;
            atoms[0].neighbor[0] = 4;
            if with_center {
                atoms[0].parity = crate::source_types::BEST_PARITY as i8;
                atoms[0].stereo_atom_parity = crate::source_types::AB_PARITY_CALC as i8;
            }
            atoms[1].valence = 1;
            atoms[1].neighbor[0] = 4;
            atoms[2].valence = 1;
            atoms[2].neighbor[0] = 5;
            atoms[3].valence = 1;
            atoms[3].neighbor[0] = 5;
            atoms[4].valence = 3;
            atoms[4].neighbor[..3].copy_from_slice(&[5, 0, 1]);
            atoms[4].parity = crate::source_types::BEST_PARITY as i8;
            atoms[4].stereo_bond_neighbor[0] = 6;
            atoms[4].stereo_bond_z_prod[0] = z_product;
            atoms[4].stereo_bond_parity[0] = crate::source_types::AB_PARITY_CALC as i8;
            atoms[5].valence = 3;
            atoms[5].neighbor[..3].copy_from_slice(&[4, 2, 3]);
            atoms[5].parity = crate::source_types::BEST_PARITY as i8;
            atoms[5].stereo_bond_neighbor[0] = 5;
            atoms[5].stereo_bond_z_prod[0] = z_product;
            atoms[5].stereo_bond_parity[0] = crate::source_types::AB_PARITY_CALC as i8;
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let order = heap.allocate_model_storage((0_u16..6).collect()).unwrap();
            let rank = heap.allocate_model_storage((1_u16..=6).collect()).unwrap();
            let symmetry = heap.allocate_model_storage((1_u16..=6).collect()).unwrap();
            let linear_input = heap.allocate_model_storage((1_u16..=6).collect()).unwrap();
            let hydrogens = heap.allocate_model_storage(vec![0_i16; 6]).unwrap();
            let scratch_rank = heap.allocate_model_storage(vec![7_u16; 6]).unwrap();
            let scratch_order = heap.allocate_model_storage(vec![8_u16; 6]).unwrap();
            let rank_stack = heap
                .allocate_model_storage(vec![
                    SourceMutPointer::null(),
                    SourceMutPointer::null(),
                    scratch_rank,
                    scratch_order,
                    SourceMutPointer::null(),
                ])
                .unwrap();
            let neighbors = CreateNeighList(heap, 6, 6, atoms.as_const(), 0, None).unwrap();
            let bcn = heap
                .allocate_model_storage(vec![BCN {
                    pRankStack: rank_stack,
                    nMaxLenRankStack: 5,
                    num_max: 6,
                    num_at_tg: 6,
                    num_atoms: 6,
                    ftcn: [
                        FTCN {
                            LinearCt: linear_input,
                            nLenLinearCt: 6,
                            nLenLinearCtAtOnly: 6,
                            PartitionCt: Partition {
                                AtNumber: order,
                                Rank: rank,
                            },
                            nSymmRankCt: symmetry,
                            nNumH: hydrogens,
                            NeighList: neighbors,
                            ..FTCN::default()
                        },
                        FTCN::default(),
                    ],
                    ..BCN::default()
                }])
                .unwrap();
            let linear = heap.allocate_model_storage(vec![0_u16; 6]).unwrap();
            let normal_double = heap
                .allocate_model_storage(vec![Default::default(); 6])
                .unwrap();
            let inverse_double = heap
                .allocate_model_storage(vec![Default::default(); 6])
                .unwrap();
            let normal_carb = heap
                .allocate_model_storage(vec![Default::default(); 6])
                .unwrap();
            let inverse_carb = heap
                .allocate_model_storage(vec![Default::default(); 6])
                .unwrap();
            let state = CANON_STAT {
                pBCN: bcn,
                LinearCT: linear,
                nNum_H: heap.allocate_model_storage(vec![0_i8; 6]).unwrap(),
                nSymmRank: heap.allocate_model_storage(vec![0_u16; 6]).unwrap(),
                nCanonOrd: heap.allocate_model_storage(vec![0_u16; 6]).unwrap(),
                LinearCTStereoDble: normal_double,
                LinearCTStereoDbleInv: inverse_double,
                LinearCTStereoCarb: normal_carb,
                LinearCTStereoCarbInv: inverse_carb,
                nMaxLenLinearCTStereoDble: 6,
                nMaxLenLinearCTStereoCarb: 6,
                bRankUsedForStereo: heap.allocate_model_storage(vec![0_i8; 6]).unwrap(),
                bAtomUsedForStereo: heap.allocate_model_storage(vec![0_i8; 6]).unwrap(),
                nPrevAtomNumber: heap.allocate_model_storage(vec![0_u16; 6]).unwrap(),
                nCanonOrdStereo: heap.allocate_model_storage(vec![0_u16; 6]).unwrap(),
                nCanonOrdStereoInv: heap.allocate_model_storage(vec![0_u16; 6]).unwrap(),
                nMode: u64::from(CMODE_STEREO | CMODE_REDNDNT_STEREO),
                ..CANON_STAT::default()
            };
            (atoms, state)
        }

        fn enable_double_bond_isotopic_stereo(
            heap: &mut SourceHeap,
            atoms: SourceMutPointer<sp_ATOM>,
            state: &mut CANON_STAT,
        ) {
            for atom in &mut heap.slice_mut(atoms).unwrap()[4..=5] {
                atom.parity2 = atom.parity;
                atom.stereo_bond_neighbor2 = atom.stereo_bond_neighbor;
                atom.stereo_bond_ord2 = atom.stereo_bond_ord;
                atom.stereo_bond_z_prod2 = atom.stereo_bond_z_prod;
                atom.stereo_bond_parity2 = atom.stereo_bond_parity;
            }
            let center = &mut heap.slice_mut(atoms).unwrap()[0];
            center.parity2 = center.parity;
            center.stereo_atom_parity2 = center.stereo_atom_parity;
            let isotope_order = heap.allocate_model_storage((0_u16..6).collect()).unwrap();
            let isotope_rank = heap.allocate_model_storage((1_u16..=6).collect()).unwrap();
            let isotope_symmetry = heap.allocate_model_storage((1_u16..=6).collect()).unwrap();
            heap.slice_mut(state.pBCN).unwrap()[0].ftcn[0].PartitionCtIso = Partition {
                AtNumber: isotope_order,
                Rank: isotope_rank,
            };
            heap.slice_mut(state.pBCN).unwrap()[0].ftcn[0].nSymmRankCtIso = isotope_symmetry;
            state.LinearCTIsotopic = heap
                .allocate_model_storage(vec![AT_ISOTOPIC::default(); 6])
                .unwrap();
            state.nMaxLenLinearCTIsotopic = 6;
            state.nSymmRankIsotopic = heap.allocate_model_storage(vec![0_u16; 6]).unwrap();
            state.nCanonOrdIsotopic = heap.allocate_model_storage(vec![0_u16; 6]).unwrap();
            state.LinearCTIsotopicStereoDble = heap
                .allocate_model_storage(vec![Default::default(); 6])
                .unwrap();
            state.LinearCTIsotopicStereoDbleInv = heap
                .allocate_model_storage(vec![Default::default(); 6])
                .unwrap();
            state.LinearCTIsotopicStereoCarb = heap
                .allocate_model_storage(vec![Default::default(); 6])
                .unwrap();
            state.LinearCTIsotopicStereoCarbInv = heap
                .allocate_model_storage(vec![Default::default(); 6])
                .unwrap();
            state.nMaxLenLinearCTIsotopicStereoDble = 6;
            state.nMaxLenLinearCTIsotopicStereoCarb = 6;
            state.nCanonOrdIsotopicStereo = heap.allocate_model_storage(vec![0_u16; 6]).unwrap();
            state.nCanonOrdIsotopicStereoInv = heap.allocate_model_storage(vec![0_u16; 6]).unwrap();
            state.nMode =
                u64::from(CMODE_STEREO | CMODE_ISO | CMODE_ISO_STEREO | CMODE_REDNDNT_STEREO);
        }

        let mut heap = SourceHeap::default();
        let (atoms, mut state, rank_stack) = fixture(&mut heap, false);
        let initial_stack = heap.slice(rank_stack.as_const()).unwrap().to_vec();
        let original_linear = state.LinearCT;
        let original_linear2 = state.LinearCT2;
        assert_eq!(
            Canon_INChI3(
                &mut heap,
                &mut INCHI_CLOCK::default(),
                100,
                None,
                None,
                2,
                2,
                atoms,
                &mut state,
                &mut CANON_GLOBALS::default(),
                0,
                0,
            ),
            Ok(2)
        );
        assert_eq!(state.LinearCT, original_linear2);
        assert_eq!(state.LinearCT2, original_linear);
        assert_eq!(heap.slice(state.LinearCT.as_const()).unwrap(), &[2, 1]);
        assert_eq!(heap.slice(state.nNum_H.as_const()).unwrap(), &[1, -1]);
        assert_eq!(heap.slice(state.nNum_H_fixed.as_const()).unwrap(), &[2, -1]);
        assert_eq!(heap.slice(state.nSymmRank.as_const()).unwrap(), &[1, 2]);
        assert_eq!(heap.slice(state.nCanonOrd.as_const()).unwrap(), &[0, 1]);
        assert_eq!(state.nLenLinearCT, 2);
        assert_eq!(state.nLenLinearCTAtOnly, 1);
        assert_eq!(state.nLenCanonOrd, 2);
        assert_eq!(state.lTotalTime, 0);
        assert_eq!(
            heap.slice(rank_stack.as_const()).unwrap(),
            &[
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                initial_stack[2],
                initial_stack[3],
                SourceMutPointer::null(),
            ]
        );

        let mut isotope_heap = SourceHeap::default();
        let (isotope_atoms, mut isotope_state, isotope_stack) = fixture(&mut isotope_heap, false);
        {
            let atoms = isotope_heap.slice_mut(isotope_atoms).unwrap();
            atoms[0].num_iso_H = [1, 2, 3];
            atoms[1].iso_atw_diff = -2;
        }
        let isotope_atom_number = isotope_heap.allocate_model_storage(vec![1_u16, 0]).unwrap();
        let isotope_rank = isotope_heap.allocate_model_storage(vec![2_u16, 1]).unwrap();
        let isotope_symmetry = isotope_heap.allocate_model_storage(vec![2_u16, 1]).unwrap();
        let exchange_h = isotope_heap.allocate_model_storage(vec![0_i8, 7]).unwrap();
        isotope_heap.slice_mut(isotope_state.pBCN).unwrap()[0].ftcn[0].PartitionCtIso = Partition {
            AtNumber: isotope_atom_number,
            Rank: isotope_rank,
        };
        isotope_heap.slice_mut(isotope_state.pBCN).unwrap()[0].ftcn[0].nSymmRankCtIso =
            isotope_symmetry;
        isotope_heap.slice_mut(isotope_state.pBCN).unwrap()[0].ftcn[0].iso_exchg_atnos = exchange_h;
        let isotope_output = isotope_heap
            .allocate_model_storage(vec![AT_ISOTOPIC::default(); 2])
            .unwrap();
        let isotope_symmetry_output = isotope_heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
        let isotope_order_output = isotope_heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
        let exchange_h_output = isotope_heap.allocate_model_storage(vec![9_i8; 2]).unwrap();
        isotope_state.LinearCTIsotopic = isotope_output;
        isotope_state.nMaxLenLinearCTIsotopic = 2;
        isotope_state.nSymmRankIsotopic = isotope_symmetry_output;
        isotope_state.nCanonOrdIsotopic = isotope_order_output;
        isotope_state.nExchgIsoH = exchange_h_output;
        isotope_heap.trace_source_allocations();
        assert_eq!(
            Canon_INChI3(
                &mut isotope_heap,
                &mut INCHI_CLOCK::default(),
                150,
                None,
                None,
                2,
                2,
                isotope_atoms,
                &mut isotope_state,
                &mut CANON_GLOBALS::default(),
                CMODE_ISO,
                0,
            ),
            Ok(2)
        );
        assert_eq!(isotope_heap.source_allocation_calls(), 0);
        assert_eq!(
            isotope_heap.slice(isotope_output.as_const()).unwrap(),
            &[
                AT_ISOTOPIC {
                    at_num: 1,
                    num_1H: 0,
                    num_D: 0,
                    num_T: 0,
                    iso_atw_diff: -2,
                },
                AT_ISOTOPIC {
                    at_num: 2,
                    num_1H: 1,
                    num_D: 2,
                    num_T: 3,
                    iso_atw_diff: 0,
                },
            ]
        );
        assert_eq!(isotope_state.nLenLinearCTIsotopic, 2);
        assert_eq!(isotope_state.nLenCanonOrdIsotopic, 2);
        assert_eq!(
            isotope_heap
                .slice(isotope_symmetry_output.as_const())
                .unwrap(),
            &[2, 1]
        );
        assert_eq!(
            isotope_heap.slice(isotope_order_output.as_const()).unwrap(),
            &[1, 0]
        );
        assert_eq!(
            isotope_heap.slice(exchange_h_output.as_const()).unwrap(),
            &[1, 0]
        );
        assert_eq!(
            isotope_heap.slice(isotope_stack.as_const()).unwrap()[..2],
            [SourceMutPointer::null(), SourceMutPointer::null()]
        );

        let mut exchange_only_heap = SourceHeap::default();
        let (exchange_only_atoms, mut exchange_only_state, _) =
            fixture(&mut exchange_only_heap, false);
        let exchange_only_order = exchange_only_heap
            .allocate_model_storage(vec![0_u16, 1])
            .unwrap();
        let exchange_only_rank = exchange_only_heap
            .allocate_model_storage(vec![1_u16, 2])
            .unwrap();
        let exchange_only_symmetry = exchange_only_heap
            .allocate_model_storage(vec![1_u16, 2])
            .unwrap();
        exchange_only_heap
            .slice_mut(exchange_only_state.pBCN)
            .unwrap()[0]
            .ftcn[0]
            .PartitionCtIso = Partition {
            AtNumber: exchange_only_order,
            Rank: exchange_only_rank,
        };
        exchange_only_heap
            .slice_mut(exchange_only_state.pBCN)
            .unwrap()[0]
            .ftcn[0]
            .nSymmRankCtIso = exchange_only_symmetry;
        let isotopic_endpoints = exchange_only_heap
            .allocate_model_storage(vec![1_u16, 0])
            .unwrap();
        exchange_only_state.t_group_info = exchange_only_heap
            .allocate_model_storage(vec![T_GROUP_INFO {
                nIsotopicEndpointAtomNumber: isotopic_endpoints,
                nNumIsotopicEndpoints: 2,
                bTautFlagsDone: u64::from(TG_FLAG_FOUND_ISOTOPIC_H_DONE),
                ..T_GROUP_INFO::default()
            }])
            .unwrap();
        exchange_only_state.nSymmRankIsotopic = exchange_only_heap
            .allocate_model_storage(vec![0_u16; 2])
            .unwrap();
        exchange_only_state.nCanonOrdIsotopic = exchange_only_heap
            .allocate_model_storage(vec![0_u16; 2])
            .unwrap();
        exchange_only_state.nMaxLenLinearCTIsotopic = 17;
        assert_eq!(
            Canon_INChI3(
                &mut exchange_only_heap,
                &mut INCHI_CLOCK::default(),
                160,
                None,
                None,
                2,
                2,
                exchange_only_atoms,
                &mut exchange_only_state,
                &mut CANON_GLOBALS::default(),
                CMODE_ISO,
                0,
            ),
            Ok(2)
        );
        assert_eq!(exchange_only_state.nMaxLenLinearCTIsotopic, 17);
        assert_eq!(exchange_only_state.nLenCanonOrdIsotopic, 2);
        assert_eq!(
            exchange_only_heap
                .slice(exchange_only_state.nCanonOrdIsotopic.as_const())
                .unwrap(),
            &[0, 1]
        );

        let mut isotope_overflow_heap = SourceHeap::default();
        let (isotope_overflow_atoms, mut isotope_overflow_state, _) =
            fixture(&mut isotope_overflow_heap, false);
        for (index, atom) in isotope_overflow_heap
            .slice_mut(isotope_overflow_atoms)
            .unwrap()
            .iter_mut()
            .enumerate()
        {
            atom.iso_atw_diff = index as i8 + 1;
        }
        let isotope_overflow_order = isotope_overflow_heap
            .allocate_model_storage(vec![0_u16, 1])
            .unwrap();
        let isotope_overflow_rank = isotope_overflow_heap
            .allocate_model_storage(vec![1_u16, 2])
            .unwrap();
        let isotope_overflow_symmetry = isotope_overflow_heap
            .allocate_model_storage(vec![1_u16, 2])
            .unwrap();
        isotope_overflow_heap
            .slice_mut(isotope_overflow_state.pBCN)
            .unwrap()[0]
            .ftcn[0]
            .PartitionCtIso = Partition {
            AtNumber: isotope_overflow_order,
            Rank: isotope_overflow_rank,
        };
        isotope_overflow_heap
            .slice_mut(isotope_overflow_state.pBCN)
            .unwrap()[0]
            .ftcn[0]
            .nSymmRankCtIso = isotope_overflow_symmetry;
        let isotope_overflow_output = isotope_overflow_heap
            .allocate_model_storage(vec![AT_ISOTOPIC::default()])
            .unwrap();
        isotope_overflow_state.LinearCTIsotopic = isotope_overflow_output;
        isotope_overflow_state.nMaxLenLinearCTIsotopic = 1;
        assert_eq!(
            Canon_INChI3(
                &mut isotope_overflow_heap,
                &mut INCHI_CLOCK::default(),
                165,
                None,
                None,
                2,
                2,
                isotope_overflow_atoms,
                &mut isotope_overflow_state,
                &mut CANON_GLOBALS::default(),
                CMODE_ISO,
                0,
            ),
            Ok(CT_OVERFLOW)
        );
        assert_eq!(
            isotope_overflow_heap
                .slice(isotope_overflow_output.as_const())
                .unwrap()[0]
                .iso_atw_diff,
            1
        );
        assert!(isotope_overflow_state.NeighList.is_null());

        let mut taut_heap = SourceHeap::default();
        let taut_atoms = taut_heap
            .allocate_model_storage(vec![sp_ATOM::default(); 2])
            .unwrap();
        let taut_rank = taut_heap
            .allocate_model_storage(vec![2_u16, 1, 3, 4])
            .unwrap();
        let taut_order = taut_heap
            .allocate_model_storage(vec![1_u16, 0, 3, 2])
            .unwrap();
        let taut_symmetry = taut_heap
            .allocate_model_storage(vec![2_u16, 2, 3, 4])
            .unwrap();
        let taut_iso_rank = taut_heap
            .allocate_model_storage(vec![2_u16, 1, 3, 4])
            .unwrap();
        let taut_iso_order = taut_heap
            .allocate_model_storage(vec![1_u16, 0, 2, 3])
            .unwrap();
        let taut_iso_symmetry = taut_heap
            .allocate_model_storage(vec![2_u16, 2, 3, 4])
            .unwrap();
        let taut_groups = taut_heap
            .allocate_model_storage(vec![
                crate::source_types::T_GROUP {
                    num: [5, 1, 10, 11, 12],
                    iWeight: 1,
                    nGroupNumber: 1,
                    nNumEndpoints: 2,
                    nFirstEndpointAtNoPos: 0,
                    ..crate::source_types::T_GROUP::default()
                },
                crate::source_types::T_GROUP {
                    num: [7, 2, 20, 21, 22],
                    iWeight: 1,
                    nGroupNumber: 2,
                    nNumEndpoints: 1,
                    nFirstEndpointAtNoPos: 2,
                    ..crate::source_types::T_GROUP::default()
                },
            ])
            .unwrap();
        let taut_endpoints = taut_heap.allocate_model_storage(vec![0_u16, 1, 1]).unwrap();
        let taut_group_work = taut_heap.allocate_model_storage(vec![99_u16; 8]).unwrap();
        let taut_info = T_GROUP_INFO {
            t_group: taut_groups,
            nEndpointAtomNumber: taut_endpoints,
            tGroupNumber: taut_group_work,
            nNumEndpoints: 3,
            num_t_groups: 2,
            max_num_t_groups: 2,
            bIgnoreIsotopic: 7,
            ..T_GROUP_INFO::default()
        };
        let taut_neighbors = CreateNeighList(
            &mut taut_heap,
            2,
            4,
            taut_atoms.as_const(),
            0,
            Some(&taut_info),
        )
        .unwrap();
        let taut_linear_input = taut_heap.allocate_model_storage(vec![2_u16, 1]).unwrap();
        let taut_h = taut_heap.allocate_model_storage(vec![0_i16; 2]).unwrap();
        let taut_rank_scratch = taut_heap.allocate_model_storage(vec![7_u16; 4]).unwrap();
        let taut_order_scratch = taut_heap.allocate_model_storage(vec![8_u16; 4]).unwrap();
        let taut_rank_stack = taut_heap
            .allocate_model_storage(vec![
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                taut_rank_scratch,
                taut_order_scratch,
                SourceMutPointer::null(),
            ])
            .unwrap();
        let taut_bcn = taut_heap
            .allocate_model_storage(vec![BCN {
                pRankStack: taut_rank_stack,
                nMaxLenRankStack: 5,
                num_max: 4,
                num_at_tg: 4,
                num_atoms: 2,
                ftcn: [
                    FTCN {
                        NeighList: taut_neighbors,
                        LinearCt: taut_linear_input,
                        nLenLinearCt: 2,
                        nLenLinearCtAtOnly: 2,
                        PartitionCt: Partition {
                            AtNumber: taut_order,
                            Rank: taut_rank,
                        },
                        nSymmRankCt: taut_symmetry,
                        nNumH: taut_h,
                        PartitionCtIso: Partition {
                            AtNumber: taut_iso_order,
                            Rank: taut_iso_rank,
                        },
                        nSymmRankCtIso: taut_iso_symmetry,
                        ..FTCN::default()
                    },
                    FTCN::default(),
                ],
                ..BCN::default()
            }])
            .unwrap();
        let taut_info_pointer = taut_heap.allocate_model_storage(vec![taut_info]).unwrap();
        let taut_linear = taut_heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
        let taut_linear_output = taut_heap.allocate_model_storage(vec![88_u16; 10]).unwrap();
        let taut_isotope_atoms = taut_heap
            .allocate_model_storage(vec![AT_ISOTOPIC::default(); 2])
            .unwrap();
        let taut_isotope_groups = taut_heap
            .allocate_model_storage(vec![AT_ISO_TGROUP::default(); 2])
            .unwrap();
        let taut_symmetry_output = taut_heap.allocate_model_storage(vec![0_u16; 4]).unwrap();
        let taut_order_output = taut_heap.allocate_model_storage(vec![0_u16; 4]).unwrap();
        let taut_group_symmetry_output = taut_heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
        let taut_group_order_output = taut_heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
        let taut_iso_symmetry_output = taut_heap.allocate_model_storage(vec![0_u16; 4]).unwrap();
        let taut_iso_order_output = taut_heap.allocate_model_storage(vec![0_u16; 4]).unwrap();
        let taut_iso_group_symmetry_output =
            taut_heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
        let taut_iso_group_order_output = taut_heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
        let mut taut_state = CANON_STAT {
            pBCN: taut_bcn,
            t_group_info: taut_info_pointer,
            LinearCT: taut_linear,
            LinearCTTautomer: taut_linear_output,
            nMaxLenLinearCTTautomer: 10,
            LinearCTIsotopic: taut_isotope_atoms,
            nMaxLenLinearCTIsotopic: 2,
            LinearCTIsotopicTautomer: taut_isotope_groups,
            nMaxLenLinearCTIsotopicTautomer: 2,
            nSymmRank: taut_symmetry_output,
            nCanonOrd: taut_order_output,
            nSymmRankTaut: taut_group_symmetry_output,
            nCanonOrdTaut: taut_group_order_output,
            nSymmRankIsotopic: taut_iso_symmetry_output,
            nCanonOrdIsotopic: taut_iso_order_output,
            nSymmRankIsotopicTaut: taut_iso_group_symmetry_output,
            nCanonOrdIsotopicTaut: taut_iso_group_order_output,
            bIgnoreIsotopic: 5,
            ..CANON_STAT::default()
        };
        taut_heap.trace_source_allocations();
        assert_eq!(
            Canon_INChI3(
                &mut taut_heap,
                &mut INCHI_CLOCK::default(),
                175,
                None,
                None,
                2,
                4,
                taut_atoms,
                &mut taut_state,
                &mut CANON_GLOBALS::default(),
                CMODE_ISO,
                0,
            ),
            Ok(2)
        );
        assert_eq!(taut_heap.source_allocation_calls(), 0);
        assert_eq!(
            taut_heap.slice(taut_linear_output.as_const()).unwrap(),
            &[1, 7, 2, 1, 2, 5, 1, 1, 2, 0]
        );
        assert_eq!(
            taut_heap.slice(taut_isotope_groups.as_const()).unwrap(),
            &[
                AT_ISO_TGROUP {
                    tgroup_num: 1,
                    num: [10, 11, 12],
                },
                AT_ISO_TGROUP {
                    tgroup_num: 2,
                    num: [20, 21, 22],
                },
            ]
        );
        assert_eq!(
            taut_heap
                .slice(taut_group_symmetry_output.as_const())
                .unwrap(),
            &[1, 2]
        );
        assert_eq!(
            taut_heap.slice(taut_group_order_output.as_const()).unwrap(),
            &[1, 0]
        );
        assert_eq!(
            taut_heap
                .slice(taut_iso_group_symmetry_output.as_const())
                .unwrap(),
            &[1, 2]
        );
        assert_eq!(
            taut_heap
                .slice(taut_iso_group_order_output.as_const())
                .unwrap(),
            &[0, 1]
        );
        assert_eq!(taut_state.nLenLinearCTTautomer, 10);
        assert_eq!(taut_state.nLenLinearCTIsotopicTautomer, 2);
        assert_eq!(taut_state.nLenCanonOrdTaut, 2);
        assert_eq!(taut_state.nLenCanonOrdIsotopicTaut, 2);
        assert_eq!(taut_state.bIgnoreIsotopic, 1);
        assert_eq!(
            taut_heap.slice(taut_info_pointer.as_const()).unwrap()[0].bIgnoreIsotopic,
            1
        );

        let iso_stereo_mode = CMODE_ISO | CMODE_ISO_STEREO | CMODE_NOEQ_STEREO;
        let mut gated_heap = SourceHeap::default();
        let (gated_atoms, mut gated_state, _) = fixture(&mut gated_heap, false);
        enable_isotopic_stereo(&mut gated_heap, gated_atoms, &mut gated_state, 0);
        let gated_atoms_before = gated_heap.slice(gated_atoms.as_const()).unwrap().to_vec();
        gated_heap.trace_source_allocations();
        assert_eq!(
            Canon_INChI3(
                &mut gated_heap,
                &mut INCHI_CLOCK::default(),
                180,
                None,
                None,
                1,
                1,
                gated_atoms,
                &mut gated_state,
                &mut CANON_GLOBALS::default(),
                iso_stereo_mode,
                0,
            ),
            Ok(1)
        );
        assert_eq!(gated_heap.source_allocation_calls(), 0);
        assert_eq!(
            gated_heap.slice(gated_atoms.as_const()).unwrap(),
            gated_atoms_before.as_slice()
        );

        let mut entered_heap = SourceHeap::default();
        let (entered_atoms, mut entered_state, _) = fixture(&mut entered_heap, false);
        enable_isotopic_stereo(&mut entered_heap, entered_atoms, &mut entered_state, 1);
        let entered_atoms_before = entered_heap
            .slice(entered_atoms.as_const())
            .unwrap()
            .to_vec();
        let normal_double_before = entered_state.LinearCTStereoDble;
        let normal_carb_before = entered_state.LinearCTStereoCarb;
        entered_heap.trace_source_allocations();
        assert_eq!(
            Canon_INChI3(
                &mut entered_heap,
                &mut INCHI_CLOCK::default(),
                185,
                Some(request_quit),
                None,
                1,
                1,
                entered_atoms,
                &mut entered_state,
                &mut CANON_GLOBALS::default(),
                iso_stereo_mode,
                0,
            ),
            Ok(CT_USER_QUIT_ERR)
        );
        assert!(entered_heap.source_allocation_calls() >= 6);
        assert_eq!(
            entered_heap.slice(entered_atoms.as_const()).unwrap(),
            entered_atoms_before.as_slice()
        );
        assert_eq!(entered_state.LinearCTStereoDble, normal_double_before);
        assert_eq!(entered_state.LinearCTStereoCarb, normal_carb_before);
        assert!(entered_state.NeighList.is_null());

        let mut descriptor_probe_heap = SourceHeap::default();
        let (descriptor_probe_atoms, mut descriptor_probe_state) =
            double_bond_fixture(&mut descriptor_probe_heap, 0, false);
        let descriptor_probe_ftcn = descriptor_probe_heap
            .slice(descriptor_probe_state.pBCN.as_const())
            .unwrap()[0]
            .ftcn[0]
            .clone();
        assert_eq!(
            FillOutStereoParities(
                &mut descriptor_probe_heap,
                descriptor_probe_atoms,
                6,
                descriptor_probe_ftcn.PartitionCt.Rank,
                descriptor_probe_ftcn.PartitionCt.AtNumber,
                descriptor_probe_ftcn.PartitionCt.Rank.as_const(),
                descriptor_probe_ftcn.PartitionCt.AtNumber.as_const(),
                &mut descriptor_probe_state,
                &mut CANON_GLOBALS::default(),
                0,
            ),
            Ok(1)
        );
        assert_eq!(descriptor_probe_state.nLenLinearCTStereoDble, 1);

        let mut descriptor_heap = SourceHeap::default();
        let (descriptor_atoms, mut descriptor_state) =
            double_bond_fixture(&mut descriptor_heap, 50, true);
        let mut descriptor_atoms_expected = descriptor_heap
            .slice(descriptor_atoms.as_const())
            .unwrap()
            .to_vec();
        for atom in &mut descriptor_atoms_expected[4..=5] {
            atom.stereo_bond_parity[0] =
                (crate::source_types::KNOWN_PARITIES_EQL | crate::source_types::WORSE_PARITY) as i8;
            atom.bHasStereoOrEquToStereo = 1;
        }
        descriptor_atoms_expected[0].stereo_atom_parity =
            (crate::source_types::KNOWN_PARITIES_EQL | crate::source_types::BEST_PARITY) as i8;
        descriptor_atoms_expected[0].bHasStereoOrEquToStereo = 1;
        let normal_double = descriptor_state.LinearCTStereoDble;
        let inverse_double = descriptor_state.LinearCTStereoDbleInv;
        let normal_carb = descriptor_state.LinearCTStereoCarb;
        let inverse_carb = descriptor_state.LinearCTStereoCarbInv;
        assert_eq!(
            Canon_INChI3(
                &mut descriptor_heap,
                &mut INCHI_CLOCK::default(),
                190,
                None,
                None,
                6,
                6,
                descriptor_atoms,
                &mut descriptor_state,
                &mut CANON_GLOBALS::default(),
                CMODE_STEREO | CMODE_REDNDNT_STEREO,
                0,
            ),
            Ok(6)
        );
        assert_eq!(descriptor_state.nLenLinearCTStereoDble, 1);
        assert_eq!(descriptor_state.nLenLinearCTStereoDbleInv, 1);
        assert_eq!(descriptor_state.nLenLinearCTStereoCarb, 1);
        assert_eq!(descriptor_state.nLenLinearCTStereoCarbInv, 1);
        let forward = descriptor_heap.slice(normal_double.as_const()).unwrap()[0].clone();
        let inverse = descriptor_heap.slice(inverse_double.as_const()).unwrap()[0].clone();
        assert_eq!((forward.at_num1, forward.at_num2), (6, 5));
        assert_eq!((inverse.at_num1, inverse.at_num2), (6, 5));
        assert_eq!(forward.parity, crate::source_types::WORSE_PARITY as u8);
        assert_eq!(inverse.parity, crate::source_types::WORSE_PARITY as u8);
        let forward_center = descriptor_heap.slice(normal_carb.as_const()).unwrap()[0].clone();
        let inverse_center = descriptor_heap.slice(inverse_carb.as_const()).unwrap()[0].clone();
        assert_eq!(forward_center.at_num, 1);
        assert_eq!(inverse_center.at_num, 1);
        assert_eq!(
            forward_center.parity,
            crate::source_types::BEST_PARITY as u8
        );
        assert_eq!(
            inverse_center.parity,
            crate::source_types::WORSE_PARITY as u8
        );
        assert_eq!(descriptor_state.nLenCanonOrdStereo, 6);
        assert_eq!(
            descriptor_heap
                .slice(descriptor_state.nCanonOrdStereo.as_const())
                .unwrap(),
            &[0, 1, 2, 3, 4, 5]
        );
        assert_eq!(
            descriptor_heap
                .slice(descriptor_state.nCanonOrdStereoInv.as_const())
                .unwrap(),
            &[0, 1, 2, 3, 4, 5]
        );
        assert_eq!(
            descriptor_heap.slice(descriptor_atoms.as_const()).unwrap(),
            descriptor_atoms_expected.as_slice()
        );
        assert!(descriptor_state.NeighList.is_null());

        let mut isotope_descriptor_heap = SourceHeap::default();
        let (isotope_descriptor_atoms, mut isotope_descriptor_state) =
            double_bond_fixture(&mut isotope_descriptor_heap, 50, true);
        enable_double_bond_isotopic_stereo(
            &mut isotope_descriptor_heap,
            isotope_descriptor_atoms,
            &mut isotope_descriptor_state,
        );
        let isotope_double = isotope_descriptor_state.LinearCTIsotopicStereoDble;
        let isotope_double_inv = isotope_descriptor_state.LinearCTIsotopicStereoDbleInv;
        let isotope_carb = isotope_descriptor_state.LinearCTIsotopicStereoCarb;
        let isotope_carb_inv = isotope_descriptor_state.LinearCTIsotopicStereoCarbInv;
        assert_eq!(
            Canon_INChI3(
                &mut isotope_descriptor_heap,
                &mut INCHI_CLOCK::default(),
                195,
                None,
                None,
                6,
                6,
                isotope_descriptor_atoms,
                &mut isotope_descriptor_state,
                &mut CANON_GLOBALS::default(),
                CMODE_STEREO | CMODE_ISO | CMODE_ISO_STEREO | CMODE_REDNDNT_STEREO,
                0,
            ),
            Ok(6)
        );
        assert_eq!(isotope_descriptor_state.nLenLinearCTIsotopic, 0);
        assert_eq!(isotope_descriptor_state.nLenLinearCTIsotopicStereoDble, 1);
        assert_eq!(
            isotope_descriptor_state.nLenLinearCTIsotopicStereoDbleInv,
            1
        );
        assert_eq!(isotope_descriptor_state.nLenLinearCTIsotopicStereoCarb, 1);
        assert_eq!(
            isotope_descriptor_state.nLenLinearCTIsotopicStereoCarbInv,
            1
        );
        let expected_isotope_double = crate::source_types::AT_STEREO_DBLE {
            at_num1: 6,
            at_num2: 5,
            parity: crate::source_types::WORSE_PARITY as u8,
        };
        assert_eq!(
            isotope_descriptor_heap
                .slice(isotope_double.as_const())
                .unwrap()[0],
            expected_isotope_double
        );
        assert_eq!(
            isotope_descriptor_heap
                .slice(isotope_double_inv.as_const())
                .unwrap()[0],
            crate::source_types::AT_STEREO_DBLE {
                parity: crate::source_types::WORSE_PARITY as u8,
                ..expected_isotope_double
            }
        );
        assert_eq!(
            isotope_descriptor_heap
                .slice(isotope_carb.as_const())
                .unwrap()[0],
            crate::source_types::AT_STEREO_CARB {
                at_num: 1,
                parity: crate::source_types::BEST_PARITY as u8,
            }
        );
        assert_eq!(
            isotope_descriptor_heap
                .slice(isotope_carb_inv.as_const())
                .unwrap()[0],
            crate::source_types::AT_STEREO_CARB {
                at_num: 1,
                parity: crate::source_types::WORSE_PARITY as u8,
            }
        );
        assert_eq!(isotope_descriptor_state.nLenCanonOrdIsotopicStereo, 6);
        assert_eq!(
            isotope_descriptor_heap
                .slice(isotope_descriptor_state.nCanonOrdIsotopicStereo.as_const())
                .unwrap(),
            &[0, 1, 2, 3, 4, 5]
        );
        assert_eq!(
            isotope_descriptor_heap
                .slice(
                    isotope_descriptor_state
                        .nCanonOrdIsotopicStereoInv
                        .as_const()
                )
                .unwrap(),
            &[0, 1, 2, 3, 4, 5]
        );
        let isotope_descriptor_atoms = isotope_descriptor_heap
            .slice(isotope_descriptor_atoms.as_const())
            .unwrap();
        assert_eq!(
            isotope_descriptor_atoms[0].stereo_atom_parity,
            (crate::source_types::KNOWN_PARITIES_EQL | crate::source_types::BEST_PARITY) as i8
        );
        assert_eq!(
            isotope_descriptor_atoms[0].stereo_atom_parity2,
            (crate::source_types::KNOWN_PARITIES_EQL | crate::source_types::BEST_PARITY) as i8
        );
        assert_eq!(isotope_descriptor_atoms[0].bHasStereoOrEquToStereo, 1);
        assert_eq!(isotope_descriptor_atoms[0].bHasStereoOrEquToStereo2, 1);
        for atom in &isotope_descriptor_atoms[4..=5] {
            let expected_parity =
                (crate::source_types::KNOWN_PARITIES_EQL | crate::source_types::WORSE_PARITY) as i8;
            assert_eq!(atom.stereo_bond_parity[0], expected_parity);
            assert_eq!(atom.stereo_bond_parity2[0], expected_parity);
            assert_eq!(atom.bHasStereoOrEquToStereo, 1);
            assert_eq!(atom.bHasStereoOrEquToStereo2, 1);
        }
        assert!(isotope_descriptor_state.NeighList.is_null());

        let mut failure_heap = SourceHeap::default();
        let (failure_atoms, mut failure_state, failure_stack) = fixture(&mut failure_heap, true);
        let initial_failure_stack = failure_heap
            .slice(failure_stack.as_const())
            .unwrap()
            .to_vec();
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            Canon_INChI3(
                &mut failure_heap,
                &mut INCHI_CLOCK::default(),
                200,
                None,
                None,
                2,
                2,
                failure_atoms,
                &mut failure_state,
                &mut CANON_GLOBALS::default(),
                CMODE_STEREO,
                0,
            ),
            Ok(CT_OUT_OF_RAM)
        );
        assert!(failure_state.NeighList.is_null());
        assert_eq!(
            failure_heap.slice(failure_stack.as_const()).unwrap(),
            &[
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                initial_failure_stack[2],
                initial_failure_stack[3],
                SourceMutPointer::null(),
            ]
        );

        let mut noeq_heap = SourceHeap::default();
        let (noeq_atoms, mut noeq_state, noeq_stack) = fixture(&mut noeq_heap, true);
        noeq_state.nMode = u64::from(CMODE_STEREO | CMODE_NOEQ_STEREO);
        noeq_heap.trace_source_allocations();
        let noeq_result = Canon_INChI3(
            &mut noeq_heap,
            &mut INCHI_CLOCK::default(),
            300,
            Some(request_quit),
            None,
            2,
            2,
            noeq_atoms,
            &mut noeq_state,
            &mut CANON_GLOBALS::default(),
            CMODE_STEREO | CMODE_NOEQ_STEREO,
            0,
        );
        let noeq_allocations = noeq_heap.source_allocation_calls();
        assert_eq!(noeq_result, Ok(CT_USER_QUIT_ERR));
        assert!(noeq_state.NeighList.is_null());
        assert!(noeq_heap.slice(noeq_stack.as_const()).unwrap()[0].is_null());
        assert!(noeq_heap.slice(noeq_stack.as_const()).unwrap()[1].is_null());

        let mut redundant_heap = SourceHeap::default();
        let (redundant_atoms, mut redundant_state, redundant_stack) =
            fixture(&mut redundant_heap, true);
        redundant_state.nMode = u64::from(CMODE_STEREO | CMODE_REDNDNT_STEREO);
        redundant_heap.trace_source_allocations();
        let redundant_result = Canon_INChI3(
            &mut redundant_heap,
            &mut INCHI_CLOCK::default(),
            400,
            Some(request_quit),
            None,
            2,
            2,
            redundant_atoms,
            &mut redundant_state,
            &mut CANON_GLOBALS::default(),
            CMODE_STEREO | CMODE_REDNDNT_STEREO,
            0,
        );
        let redundant_allocations = redundant_heap.source_allocation_calls();
        assert_eq!(redundant_result, Ok(CT_USER_QUIT_ERR));
        assert!(redundant_state.NeighList.is_null());
        assert!(redundant_heap.slice(redundant_stack.as_const()).unwrap()[0].is_null());
        assert!(redundant_heap.slice(redundant_stack.as_const()).unwrap()[1].is_null());
        assert_eq!(noeq_allocations, 7);
        assert_eq!(redundant_allocations, 9);

        for failed_allocation in 0..noeq_allocations {
            let mut allocation_heap = SourceHeap::default();
            let (allocation_atoms, mut allocation_state, allocation_stack) =
                fixture(&mut allocation_heap, true);
            let baseline_allocations = allocation_heap.live_allocation_count();
            let original_stack = allocation_heap
                .slice(allocation_stack.as_const())
                .unwrap()
                .to_vec();
            allocation_heap.fail_after_allocations(failed_allocation);
            let expected_status = if failed_allocation == 5 {
                CT_STEREOCOUNT_ERR
            } else {
                CT_OUT_OF_RAM
            };
            assert_eq!(
                Canon_INChI3(
                    &mut allocation_heap,
                    &mut INCHI_CLOCK::default(),
                    500,
                    Some(request_quit),
                    None,
                    2,
                    2,
                    allocation_atoms,
                    &mut allocation_state,
                    &mut CANON_GLOBALS::default(),
                    CMODE_STEREO | CMODE_NOEQ_STEREO,
                    0,
                ),
                Ok(expected_status),
                "source allocation ordinal {failed_allocation}"
            );
            assert_eq!(
                allocation_heap.live_allocation_count(),
                baseline_allocations,
                "source allocation ordinal {failed_allocation} leaked"
            );
            assert_eq!(
                allocation_heap.slice(allocation_stack.as_const()).unwrap(),
                original_stack.as_slice(),
                "source allocation ordinal {failed_allocation} changed rank-stack ownership"
            );
            assert!(allocation_state.NeighList.is_null());
        }
    }

    #[test]
    fn source_port__ichicano__fillisotopicatlinearct__line_780() {
        let mut empty_heap = SourceHeap::default();
        let mut untouched_len = 41;
        assert_eq!(
            FillIsotopicAtLinearCT(
                &mut empty_heap,
                7,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                9,
                &mut untouched_len,
            ),
            Ok(0)
        );
        assert_eq!(untouched_len, 41);

        let mut atoms = vec![sp_ATOM::default(); 5];
        atoms[1].iso_atw_diff = -2;
        atoms[1].num_iso_H = [1, 2, 3];
        atoms[2].endpoint = 1;
        atoms[2].num_iso_H = [4, 5, 6];
        atoms[3].iso_atw_diff = 1;
        atoms[3].cFlags = AT_FLAG_ISO_H_POINT as i8;
        atoms[3].num_iso_H = [7, 8, 9];
        atoms[4].num_iso_H = [10, 11, 12];

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(atoms).unwrap();
        let order = heap
            .allocate_model_storage(vec![4_u16, 0, 3, 2, 1])
            .unwrap();
        let sentinel = AT_ISOTOPIC {
            at_num: 99,
            num_1H: 99,
            num_D: 99,
            num_T: 99,
            iso_atw_diff: 99,
        };
        let output = heap
            .allocate_model_storage(vec![sentinel.clone(); 5])
            .unwrap();
        let mut output_len = 0;
        assert_eq!(
            FillIsotopicAtLinearCT(
                &mut heap,
                5,
                atoms.as_const(),
                order.as_const(),
                output,
                5,
                &mut output_len,
            ),
            Ok(3)
        );
        assert_eq!(output_len, 3);
        assert_eq!(
            heap.slice(output.as_const()).unwrap(),
            &[
                AT_ISOTOPIC {
                    at_num: 1,
                    num_1H: 10,
                    num_D: 11,
                    num_T: 12,
                    iso_atw_diff: 0,
                },
                AT_ISOTOPIC {
                    at_num: 3,
                    num_1H: 0,
                    num_D: 0,
                    num_T: 0,
                    iso_atw_diff: 1,
                },
                AT_ISOTOPIC {
                    at_num: 5,
                    num_1H: 1,
                    num_D: 2,
                    num_T: 3,
                    iso_atw_diff: -2,
                },
                AT_ISOTOPIC::default(),
                AT_ISOTOPIC::default(),
            ]
        );

        heap.slice_mut(output).unwrap().fill(sentinel.clone());
        let mut zero_capacity_len = 17;
        assert_eq!(
            FillIsotopicAtLinearCT(
                &mut heap,
                5,
                atoms.as_const(),
                order.as_const(),
                output,
                0,
                &mut zero_capacity_len,
            ),
            Ok(0)
        );
        assert_eq!(zero_capacity_len, 17);
        assert!(
            heap.slice(output.as_const())
                .unwrap()
                .iter()
                .all(|item| item == &sentinel)
        );

        let short_output = heap
            .allocate_model_storage(vec![sentinel.clone(); 2])
            .unwrap();
        let mut overflow_len = 77;
        assert_eq!(
            FillIsotopicAtLinearCT(
                &mut heap,
                5,
                atoms.as_const(),
                order.as_const(),
                short_output,
                2,
                &mut overflow_len,
            ),
            Ok(CT_OVERFLOW)
        );
        assert_eq!(overflow_len, 77);
        assert_eq!(
            heap.slice(short_output.as_const()).unwrap(),
            &[
                AT_ISOTOPIC {
                    at_num: 1,
                    num_1H: 10,
                    num_D: 11,
                    num_T: 12,
                    iso_atw_diff: 0,
                },
                AT_ISOTOPIC {
                    at_num: 3,
                    num_1H: 0,
                    num_D: 0,
                    num_T: 0,
                    iso_atw_diff: 1,
                },
            ]
        );

        let mismatch_output = heap
            .allocate_model_storage(vec![sentinel.clone(); 5])
            .unwrap();
        let mut mismatch_len = 2;
        assert_eq!(
            FillIsotopicAtLinearCT(
                &mut heap,
                5,
                atoms.as_const(),
                order.as_const(),
                mismatch_output,
                5,
                &mut mismatch_len,
            ),
            Ok(CT_LEN_MISMATCH)
        );
        assert_eq!(mismatch_len, 2);
        assert_eq!(heap.slice(mismatch_output.as_const()).unwrap()[2].at_num, 5);

        let negative_output = heap.allocate_model_storage(vec![sentinel.clone()]).unwrap();
        let mut negative_len = 0;
        assert_eq!(
            FillIsotopicAtLinearCT(
                &mut heap,
                -1,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                negative_output,
                1,
                &mut negative_len,
            ),
            Ok(0)
        );
        assert_eq!(negative_len, 0);
        assert_eq!(
            heap.slice(negative_output.as_const()).unwrap(),
            &[AT_ISOTOPIC::default()]
        );
    }

    #[test]
    fn source_port__ichicano__filltautlinearct2__line_858() {
        fn fixture(
            isotopic_endpoints: i32,
        ) -> (
            SourceHeap,
            CANON_GLOBALS,
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_TAUTOMER>,
            SourceMutPointer<AT_ISO_TGROUP>,
            T_GROUP_INFO,
        ) {
            let mut heap = SourceHeap::default();
            let rank = heap.allocate_model_storage(vec![2, 1, 3, 4]).unwrap();
            let atom_number = heap.allocate_model_storage(vec![1, 0, 3, 2]).unwrap();
            let symmetry = heap.allocate_model_storage(vec![2, 2, 3, 4]).unwrap();
            let rank_iso = heap.allocate_model_storage(vec![2, 1, 3, 4]).unwrap();
            let atom_number_iso = heap.allocate_model_storage(vec![1, 0, 2, 3]).unwrap();
            let symmetry_iso = heap.allocate_model_storage(vec![2, 2, 3, 4]).unwrap();
            let groups = heap
                .allocate_model_storage(vec![
                    crate::source_types::T_GROUP {
                        num: [5, 1, 10, 11, 12],
                        iWeight: 1,
                        nGroupNumber: 1,
                        nNumEndpoints: 2,
                        nFirstEndpointAtNoPos: 0,
                        ..crate::source_types::T_GROUP::default()
                    },
                    crate::source_types::T_GROUP {
                        num: [7, 2, 20, 21, 22],
                        iWeight: 1,
                        nGroupNumber: 2,
                        nNumEndpoints: 1,
                        nFirstEndpointAtNoPos: 2,
                        ..crate::source_types::T_GROUP::default()
                    },
                ])
                .unwrap();
            let endpoints = heap.allocate_model_storage(vec![0, 1, 1]).unwrap();
            let group_work = heap.allocate_model_storage(vec![99; 8]).unwrap();
            let linear = heap.allocate_model_storage(vec![88; 10]).unwrap();
            let linear_iso = heap
                .allocate_model_storage(vec![AT_ISO_TGROUP::default(); 2])
                .unwrap();
            let info = T_GROUP_INFO {
                t_group: groups,
                nEndpointAtomNumber: endpoints,
                tGroupNumber: group_work,
                nNumEndpoints: 3,
                num_t_groups: 2,
                max_num_t_groups: 2,
                nNumIsotopicEndpoints: isotopic_endpoints,
                ..T_GROUP_INFO::default()
            };
            (
                heap,
                CANON_GLOBALS::default(),
                rank,
                atom_number,
                symmetry,
                rank_iso,
                atom_number_iso,
                symmetry_iso,
                linear,
                linear_iso,
                info,
            )
        }

        let mut empty_heap = SourceHeap::default();
        let mut normal_len = 17;
        let mut isotope_len = 19;
        assert_eq!(
            FillTautLinearCT2(
                &mut empty_heap,
                &mut CANON_GLOBALS::default(),
                2,
                2,
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                0,
                &mut normal_len,
                SourceMutPointer::null(),
                0,
                &mut isotope_len,
                None,
            ),
            Ok(0)
        );
        assert_eq!((normal_len, isotope_len), (17, 19));

        let (
            mut heap,
            mut globals,
            rank,
            atom_number,
            symmetry,
            rank_iso,
            atom_number_iso,
            symmetry_iso,
            linear,
            linear_iso,
            mut info,
        ) = fixture(0);
        let mut normal_len = 0;
        let mut isotope_len = 0;
        assert_eq!(
            FillTautLinearCT2(
                &mut heap,
                &mut globals,
                2,
                4,
                1,
                rank.as_const(),
                atom_number.as_const(),
                symmetry.as_const(),
                rank_iso.as_const(),
                atom_number_iso.as_const(),
                symmetry_iso.as_const(),
                linear,
                10,
                &mut normal_len,
                linear_iso,
                2,
                &mut isotope_len,
                Some(&mut info),
            ),
            Ok(10)
        );
        assert_eq!((normal_len, isotope_len), (10, 2));
        assert_eq!(
            heap.slice(info.tGroupNumber.as_const()).unwrap(),
            &[1, 0, 1, 2, 0, 1, 1, 2]
        );
        assert_eq!(
            heap.slice(info.nEndpointAtomNumber.as_const()).unwrap(),
            &[1, 0, 1]
        );
        assert_eq!(
            heap.slice(linear.as_const()).unwrap(),
            &[1, 7, 2, 1, 2, 5, 1, 1, 2, 0]
        );
        assert_eq!(
            heap.slice(linear_iso.as_const()).unwrap(),
            &[
                AT_ISO_TGROUP {
                    tgroup_num: 1,
                    num: [10, 11, 12],
                },
                AT_ISO_TGROUP {
                    tgroup_num: 2,
                    num: [20, 21, 22],
                },
            ]
        );
        assert_eq!(globals.m_pn_RankForSort, rank.as_const());

        let (
            mut overflow_heap,
            mut overflow_globals,
            overflow_rank,
            overflow_atom_number,
            overflow_symmetry,
            overflow_rank_iso,
            overflow_atom_number_iso,
            overflow_symmetry_iso,
            overflow_linear,
            overflow_linear_iso,
            mut overflow_info,
        ) = fixture(0);
        let mut overflow_normal_len = 31;
        let mut overflow_iso_len = 32;
        assert_eq!(
            FillTautLinearCT2(
                &mut overflow_heap,
                &mut overflow_globals,
                2,
                4,
                1,
                overflow_rank.as_const(),
                overflow_atom_number.as_const(),
                overflow_symmetry.as_const(),
                overflow_rank_iso.as_const(),
                overflow_atom_number_iso.as_const(),
                overflow_symmetry_iso.as_const(),
                overflow_linear,
                9,
                &mut overflow_normal_len,
                overflow_linear_iso,
                2,
                &mut overflow_iso_len,
                Some(&mut overflow_info),
            ),
            Ok(CT_OVERFLOW)
        );
        assert_eq!((overflow_normal_len, overflow_iso_len), (31, 32));
        assert_eq!(
            overflow_heap
                .slice(overflow_info.nEndpointAtomNumber.as_const())
                .unwrap(),
            &[1, 0, 1]
        );

        let (
            mut iso_overflow_heap,
            mut iso_overflow_globals,
            iso_overflow_rank,
            iso_overflow_atom_number,
            iso_overflow_symmetry,
            iso_overflow_rank_iso,
            iso_overflow_atom_number_iso,
            iso_overflow_symmetry_iso,
            iso_overflow_linear,
            iso_overflow_linear_iso,
            mut iso_overflow_info,
        ) = fixture(0);
        let mut iso_overflow_normal_len = 0;
        let mut iso_overflow_len = 0;
        assert_eq!(
            FillTautLinearCT2(
                &mut iso_overflow_heap,
                &mut iso_overflow_globals,
                2,
                4,
                1,
                iso_overflow_rank.as_const(),
                iso_overflow_atom_number.as_const(),
                iso_overflow_symmetry.as_const(),
                iso_overflow_rank_iso.as_const(),
                iso_overflow_atom_number_iso.as_const(),
                iso_overflow_symmetry_iso.as_const(),
                iso_overflow_linear,
                10,
                &mut iso_overflow_normal_len,
                iso_overflow_linear_iso,
                1,
                &mut iso_overflow_len,
                Some(&mut iso_overflow_info),
            ),
            Ok(CT_OVERFLOW)
        );
        assert_eq!(iso_overflow_normal_len, 10);
        assert_eq!(iso_overflow_len, 0);
        assert_eq!(
            iso_overflow_heap
                .slice(iso_overflow_linear_iso.as_const())
                .unwrap()[0],
            AT_ISO_TGROUP {
                tgroup_num: 1,
                num: [10, 11, 12],
            }
        );
    }

    #[test]
    fn source_port__ichicano__fixcanonequivalenceinfo__line_1278() {
        fn allocate_case(
            heap: &mut SourceHeap,
            symmetry: &[AT_RANK],
            current: &[AT_RANK],
            order: &[AT_NUMB],
        ) -> (
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_RANK>,
            SourceMutPointer<AT_NUMB>,
        ) {
            (
                heap.allocate_model_storage(symmetry.to_vec()).unwrap(),
                heap.allocate_model_storage(current.to_vec()).unwrap(),
                heap.allocate_model_storage(vec![99; symmetry.len()])
                    .unwrap(),
                heap.allocate_model_storage(order.to_vec()).unwrap(),
            )
        }

        let mut stable_heap = SourceHeap::default();
        let (stable_symmetry, stable_current, stable_temporary, stable_order) = allocate_case(
            &mut stable_heap,
            &[1, 1, 3, 3],
            &[2, 2, 4, 4],
            &[3, 1, 2, 0],
        );
        let mut stable_globals = CANON_GLOBALS::default();
        let mut stable_changed = 9;
        assert_eq!(
            FixCanonEquivalenceInfo(
                &mut stable_heap,
                &mut stable_globals,
                4,
                stable_symmetry,
                stable_current,
                stable_temporary,
                stable_order,
                Some(&mut stable_changed),
            ),
            Ok(2)
        );
        assert_eq!(stable_changed, 0);
        assert_eq!(
            stable_heap.slice(stable_order.as_const()).unwrap(),
            &[0, 1, 2, 3]
        );
        assert_eq!(
            stable_heap.slice(stable_symmetry.as_const()).unwrap(),
            &[1, 1, 3, 3]
        );
        assert_eq!(
            stable_heap.slice(stable_current.as_const()).unwrap(),
            &[2, 2, 4, 4]
        );
        assert_eq!(
            stable_heap.slice(stable_temporary.as_const()).unwrap(),
            &[2, 2, 4, 4]
        );
        assert_eq!(stable_globals.m_pn_RankForSort, stable_symmetry.as_const());

        let mut current_only_heap = SourceHeap::default();
        let (symmetry, current, temporary, order) = allocate_case(
            &mut current_only_heap,
            &[1, 1, 3, 3],
            &[0, 0, 0, 0],
            &[0, 1, 2, 3],
        );
        let mut globals = CANON_GLOBALS::default();
        let mut changed = 0;
        assert_eq!(
            FixCanonEquivalenceInfo(
                &mut current_only_heap,
                &mut globals,
                4,
                symmetry,
                current,
                temporary,
                order,
                Some(&mut changed),
            ),
            Ok(2)
        );
        assert_eq!(changed, 2);
        assert_eq!(
            current_only_heap.slice(current.as_const()).unwrap(),
            &[2, 2, 4, 4]
        );
        assert_eq!(
            current_only_heap.slice(symmetry.as_const()).unwrap(),
            &[1, 1, 3, 3]
        );

        let source_symmetry = [1_u16, 1, 1, 1, 3, 3, 7, 7];
        let rebuilt_ranks = [4_u16, 4, 4, 4, 6, 6, 8, 8];
        let scrambled_order = [6_u16, 2, 5, 0, 7, 3, 1, 4];
        let mut both_heap = SourceHeap::default();
        let (symmetry, current, temporary, order) = allocate_case(
            &mut both_heap,
            &source_symmetry,
            &[6, 6, 6, 6, 6, 6, 8, 8],
            &scrambled_order,
        );
        let mut globals = CANON_GLOBALS::default();
        let mut changed = 0;
        assert_eq!(
            FixCanonEquivalenceInfo(
                &mut both_heap,
                &mut globals,
                8,
                symmetry,
                current,
                temporary,
                order,
                Some(&mut changed),
            ),
            Ok(3)
        );
        assert_eq!(changed, 3);
        assert_eq!(
            both_heap.slice(order.as_const()).unwrap(),
            &[0, 1, 2, 3, 4, 5, 6, 7]
        );
        assert_eq!(
            both_heap.slice(temporary.as_const()).unwrap(),
            &rebuilt_ranks
        );
        assert_eq!(both_heap.slice(current.as_const()).unwrap(), &rebuilt_ranks);
        assert_eq!(
            both_heap.slice(symmetry.as_const()).unwrap(),
            &[1, 1, 1, 1, 5, 5, 7, 7]
        );

        let mut symmetry_only_heap = SourceHeap::default();
        let (symmetry, current, temporary, order) = allocate_case(
            &mut symmetry_only_heap,
            &source_symmetry,
            &rebuilt_ranks,
            &scrambled_order,
        );
        let mut globals = CANON_GLOBALS::default();
        let mut changed = 0;
        assert_eq!(
            FixCanonEquivalenceInfo(
                &mut symmetry_only_heap,
                &mut globals,
                8,
                symmetry,
                current,
                temporary,
                order,
                Some(&mut changed),
            ),
            Ok(3)
        );
        assert_eq!(changed, 1);
        assert_eq!(
            symmetry_only_heap.slice(symmetry.as_const()).unwrap(),
            &[1, 1, 1, 1, 5, 5, 7, 7]
        );
        assert_eq!(
            symmetry_only_heap.slice(current.as_const()).unwrap(),
            &rebuilt_ranks
        );

        let mut empty_heap = SourceHeap::default();
        let mut empty_globals = CANON_GLOBALS::default();
        let mut empty_changed = 7;
        assert_eq!(
            FixCanonEquivalenceInfo(
                &mut empty_heap,
                &mut empty_globals,
                0,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                Some(&mut empty_changed),
            ),
            Ok(1)
        );
        assert_eq!(empty_changed, 0);

        let mut error_heap = SourceHeap::default();
        let (short_symmetry, current, temporary, order) =
            allocate_case(&mut error_heap, &[1], &[2, 2], &[0, 1]);
        let mut error_globals = CANON_GLOBALS::default();
        let mut error_changed = 6;
        assert_eq!(
            FixCanonEquivalenceInfo(
                &mut error_heap,
                &mut error_globals,
                2,
                short_symmetry,
                current,
                temporary,
                order,
                Some(&mut error_changed),
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(error_changed, 6);
        assert_eq!(error_globals.m_pn_RankForSort, short_symmetry.as_const());
    }

    #[test]
    fn source_port__ichicano__inchiclock__line_111() {
        assert_eq!(InchiClock(-1), 0);
        assert_eq!(InchiClock(0), 0);
        assert_eq!(InchiClock(1), 1);
        assert_eq!(InchiClock(i64::MIN), i64::MIN);
        assert_eq!(InchiClock(i64::MAX), i64::MAX);
    }

    #[test]
    fn source_port__ichicano__inchitimeget__line_146() {
        let mut tick_end = inchiTime { clockTime: 17 };

        InchiTimeGet(&mut tick_end, -1);
        assert_eq!(tick_end.clockTime, 0);

        for clock_result in [0, 1, i64::MIN, i64::MAX] {
            tick_end.clockTime = 17;
            InchiTimeGet(&mut tick_end, clock_result);
            assert_eq!(tick_end.clockTime, clock_result);
        }
    }

    #[test]
    fn source_port__ichicano__fillmaxminclock__line_128() {
        let mut clock = INCHI_CLOCK {
            m_MaxPositiveClock: 0,
            m_MinNegativeClock: 11,
            m_HalfMaxPositiveClock: 12,
            m_HalfMinNegativeClock: 13,
        };
        FillMaxMinClock(&mut clock);
        assert_eq!(clock.m_MaxPositiveClock, i64::MAX);
        assert_eq!(clock.m_MinNegativeClock, -i64::MAX);
        assert_eq!(clock.m_HalfMaxPositiveClock, i64::MAX / 2);
        assert_eq!(clock.m_HalfMinNegativeClock, -i64::MAX / 2);

        let unchanged = INCHI_CLOCK {
            m_MaxPositiveClock: 7,
            m_MinNegativeClock: 8,
            m_HalfMaxPositiveClock: 9,
            m_HalfMinNegativeClock: 10,
        };
        let mut clock = unchanged.clone();
        FillMaxMinClock(&mut clock);
        assert_eq!(clock, unchanged);
    }

    #[test]
    fn source_port__ichicano__inchitimemsecdiff__line_151() {
        let mut clock = INCHI_CLOCK::default();
        assert_eq!(InchiTimeMsecDiff(&mut clock, None, None), 0);
        assert_eq!(clock.m_MaxPositiveClock, i64::MAX);

        let tick = |clock_time| inchiTime {
            clockTime: clock_time,
        };
        let cases = [
            (2_500_000, 1_000_000, 1_500),
            (1_000_000, 2_500_000, -1_500),
            (-1_000_000, -2_500_000, 1_500),
            (-2_500_000, -1_000_000, -1_500),
            (1_500_000, -500_000, 2_000),
            (-500_000, 1_500_000, -2_000),
            (0, 0, 0),
        ];
        for (end, start, expected) in cases {
            assert_eq!(
                InchiTimeMsecDiff(&mut clock, Some(&tick(end)), Some(&tick(start))),
                expected
            );
        }

        assert_eq!(
            InchiTimeMsecDiff(
                &mut clock,
                Some(&tick(i64::MAX - 2_000)),
                Some(&tick(-i64::MAX + 3_000)),
            ),
            -5
        );
        assert_eq!(
            InchiTimeMsecDiff(
                &mut clock,
                Some(&tick(-i64::MAX + 3_000)),
                Some(&tick(i64::MAX - 2_000)),
            ),
            5
        );

        let mut initialized = INCHI_CLOCK::default();
        assert_eq!(InchiTimeMsecDiff(&mut initialized, Some(&tick(1)), None), 0);
        assert_eq!(initialized.m_MaxPositiveClock, i64::MAX);
    }

    #[test]
    fn source_port__ichicano__inchitimeelapsed__line_223() {
        let mut clock = INCHI_CLOCK::default();
        assert_eq!(InchiTimeElapsed(&mut clock, None, 9_000_000), 0);
        assert_eq!(clock, INCHI_CLOCK::default());

        let start = inchiTime {
            clockTime: 1_000_000,
        };
        assert_eq!(InchiTimeElapsed(&mut clock, Some(&start), 2_500_000), 1_500);
        assert_eq!(clock.m_MaxPositiveClock, i64::MAX);

        assert_eq!(InchiTimeElapsed(&mut clock, Some(&start), -1), -1_000);

        let wrapped_start = inchiTime {
            clockTime: i64::MAX - 2_000,
        };
        assert_eq!(
            InchiTimeElapsed(&mut clock, Some(&wrapped_start), -i64::MAX + 3_000,),
            5
        );
    }

    #[test]
    fn source_port__ichicano__inchitimeaddmsec__line_234() {
        let mut untouched_clock = INCHI_CLOCK::default();
        InchiTimeAddMsec(&mut untouched_clock, None, u64::MAX);
        assert_eq!(untouched_clock, INCHI_CLOCK::default());

        for (milliseconds, initial, expected) in [
            (0_u64, 17_i64, 17_i64),
            (1, 17, 1_017),
            (999, -999_000, 0),
            (1_001, -1_000, 999_999),
            (u64::from(u32::MAX), 0, 4_294_967_295_000),
        ] {
            let mut clock = INCHI_CLOCK::default();
            let mut tick = inchiTime { clockTime: initial };
            InchiTimeAddMsec(&mut clock, Some(&mut tick), milliseconds);
            assert_eq!(tick.clockTime, expected, "milliseconds={milliseconds}");
            assert_eq!(clock.m_MaxPositiveClock, i64::MAX);
            assert_eq!(clock.m_MinNegativeClock, -i64::MAX);
            assert_eq!(clock.m_HalfMaxPositiveClock, i64::MAX / 2);
            assert_eq!(clock.m_HalfMinNegativeClock, -i64::MAX / 2);
        }

        let initialized = INCHI_CLOCK {
            m_MaxPositiveClock: 101,
            m_MinNegativeClock: -101,
            m_HalfMaxPositiveClock: 50,
            m_HalfMinNegativeClock: -50,
        };
        let mut clock = initialized.clone();
        let mut tick = inchiTime { clockTime: 5 };
        InchiTimeAddMsec(&mut clock, Some(&mut tick), 2);
        assert_eq!(clock, initialized);
        assert_eq!(tick.clockTime, 2_005);
    }

    #[test]
    fn source_port__ichicano__binchitimeisover__line_257() {
        let tick = |clock_time| inchiTime {
            clockTime: clock_time,
        };

        let mut clock = INCHI_CLOCK::default();
        assert_eq!(bInchiTimeIsOver(&mut clock, None, 9), 0);
        assert_eq!(clock.m_MaxPositiveClock, i64::MAX);
        assert_eq!(clock.m_MinNegativeClock, -i64::MAX);
        assert_eq!(clock.m_HalfMaxPositiveClock, i64::MAX / 2);
        assert_eq!(clock.m_HalfMinNegativeClock, -i64::MAX / 2);

        let cases = [
            (0, 0, 0),
            (2, 1, 1),
            (1, 2, 0),
            (-1, -2, 1),
            (-2, -1, 0),
            (1, -1, 1),
            (-1, 1, 0),
            (i64::MAX - 1, -i64::MAX + 1, 0),
            (-i64::MAX + 1, i64::MAX - 1, 1),
            (i64::MAX / 2 - 1, -i64::MAX / 2 + 1, 1),
            (-i64::MAX / 2 + 1, i64::MAX / 2 - 1, 0),
            (-1, -1, 1),
        ];
        for (clock_result, start, expected) in cases {
            let mut clock = INCHI_CLOCK::default();
            assert_eq!(
                bInchiTimeIsOver(&mut clock, Some(&tick(start)), clock_result),
                expected,
                "clock_result={clock_result}, start={start}"
            );
        }

        let mut clock = INCHI_CLOCK::default();
        assert_eq!(bInchiTimeIsOver(&mut clock, Some(&tick(-1)), -1), 1);
    }
}
