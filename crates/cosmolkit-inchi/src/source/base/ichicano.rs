use crate::source_types::{INCHI_CLOCK, clock_t, inchiTime};

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

#[cfg(test)]
mod tests {
    use super::*;

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
