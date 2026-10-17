use super::clock_t;

/// Private sampling transport. Values preserve the existing deterministic
/// source-unit inputs; production entrypoints use LiveSourceClock.
pub(crate) trait SourceClock: Copy {
    fn sample(self) -> clock_t;
}

impl SourceClock for i64 {
    #[inline]
    fn sample(self) -> clock_t { self }
}

impl SourceClock for i32 {
    #[inline]
    fn sample(self) -> clock_t { i64::from(self) }
}

#[derive(Clone, Copy)]
pub(crate) struct LiveSourceClock;

impl SourceClock for LiveSourceClock {
    fn sample(self) -> clock_t {
        // VERBATIM SOURCE FRAGMENT: ichicano.c113..119
        //     clock_t c = clock( );
        //     if (c != (clock_t) -1)
        //     {
        //         return c;
        //     }
        //     return 0;
        // INCHI✔️❗: Native samples whole-process user+system CPU time.
        // The existing LP64 port uses 1_000_000 ticks/second. Failure is -1
        // here so InchiClock retains the source's failure-to-zero mapping.
        // Human 2026-10-17 explicitly selects elapsed web-time on Wasm:
        // waiting counts there, unlike the native CPU timeout contract.
        // INCHI-PERF-PLATFORM-CLOCK: O(1), no allocation, cloning or scans.
        // One backend query per reached source sample; Wasm initializes one
        // shared epoch so separate generation stages use the same origin.
        // Backend/conversion constants remain OPEN in performance_pending.md.
        #[cfg(not(target_family = "wasm"))]
        let micros = match cpu_time::ProcessTime::try_now() {
            Ok(time) => time.as_duration().as_micros(),
            Err(_) => return -1,
        };
        #[cfg(target_family = "wasm")]
        let micros = {
            static EPOCH: std::sync::OnceLock<web_time::Instant> = std::sync::OnceLock::new();
            let now = web_time::Instant::now();
            let epoch = EPOCH.get_or_init(|| now);
            now.duration_since(*epoch).as_micros()
        };
        clock_t::try_from(micros).unwrap_or(-1)
    }
}
