use std::sync::Mutex;

const MULTIPLIER: u32 = 48_271;
const MODULUS: u32 = 2_147_483_647;

#[derive(Debug, Clone)]
pub struct RdkitRandomEngine {
    state: u32,
}

static RDKIT_RANDOM_GENERATOR: Mutex<RdkitRandomEngine> =
    Mutex::new(RdkitRandomEngine { state: 42 });

impl RdkitRandomEngine {
    pub fn from_seed(seed: u32) -> Self {
        // BEGIN BOOST RANDOM FUNCTION linear_congruential_engine arithmetic constructor
        // Boost❗✔️: BOOST_RANDOM_DETAIL_ARITHMETIC_CONSTRUCTOR(linear_congruential_engine,
        // Boost❗✔️:                                                IntType, x0)
        // Boost❗✔️: { seed(x0); }
        // END BOOST RANDOM FUNCTION linear_congruential_engine arithmetic constructor
        // Behavior: the fixed minstd_rand specialization stores its normalized
        // seed before the first draw. The seed implementation below preserves
        // the source's zero-to-one normalization.
        // Complexity: construction and normalization are constant-time and use
        // no allocation, matching the scalar source constructor.
        let mut engine = Self { state: 1 };
        engine.seed(seed);
        engine
    }

    fn seed(&mut self, seed: u32) {
        // BEGIN BOOST RANDOM FUNCTION linear_congruential_engine::seed
        // Boost❗✔️: BOOST_RANDOM_DETAIL_ARITHMETIC_SEED(linear_congruential_engine, IntType, x0_)
        // Boost❗✔️: {
        // Boost❗✔️:     // Work around a msvc 12/14 optimizer bug, which causes
        // Boost❗✔️:     // the line _x = 1 to run unconditionally sometimes.
        // Boost❗✔️:     // Creating a local copy seems to work around it.
        // Boost❗✔️:     IntType x0 = x0_;
        // Boost❗✔️:     // wrap _x if it doesn't fit in the destination
        // Boost❗✔️:     if(modulus == 0) {
        // Boost❗✔️:         _x = x0;
        // Boost❗✔️:     } else {
        // Boost❗✔️:         _x = x0 % modulus;
        // Boost❗✔️:     }
        // Boost❗✔️:     // handle negative seeds
        // Boost❗✔️:     if(_x < 0) {
        // Boost❗✔️:         _x += modulus;
        // Boost❗✔️:     }
        // Boost❗✔️:     // adjust to the correct range
        // Boost❗✔️:     if(increment == 0 && _x == 0) {
        // Boost❗✔️:         _x = 1;
        // Boost❗✔️:     }
        // Boost❗✔️:     BOOST_ASSERT(_x >= (min)());
        // Boost❗✔️:     BOOST_ASSERT(_x <= (max)());
        // Boost❗✔️: }
        // END BOOST RANDOM FUNCTION linear_congruential_engine::seed
        // Behavior: IntType is u32 here, so signed negative-seed correction is
        // unreachable; source remainder and zero-to-one normalization remain.
        // Complexity: one u32 remainder and a conditional, with no allocation.
        let state = seed % MODULUS;
        self.state = if state == 0 { 1 } else { state };
    }

    pub fn next_u32(&mut self) -> u32 {
        // BEGIN BOOST RANDOM FUNCTION linear_congruential_engine::operator()
        // Boost❗✔️: typedef linear_congruential_engine<uint32_t, 48271, 0, 2147483647> minstd_rand;
        // Boost❗✔️: IntType operator()()
        // Boost❗✔️: {
        // Boost❗✔️:     _x = const_mod<IntType, m>::mult_add(a, _x, c);
        // Boost❗✔️:     return _x;
        // Boost❗✔️: }
        // END BOOST RANDOM FUNCTION linear_congruential_engine::operator()
        // BEGIN BOOST RANDOM HELPER const_mod::mult_add
        // Boost❗✔️: static IntType mult_add(IntType a, IntType x, IntType c)
        // Boost❗✔️: {
        // Boost❗✔️:   if(((unsigned_m() - 1) & unsigned_m()) == 0)
        // Boost❗✔️:     return (unsigned_type(a) * unsigned_type(x) + unsigned_type(c)) & (unsigned_m() - 1);
        // Boost❗✔️:   else if(a == 0)
        // Boost❗✔️:     return c;
        // Boost❗✔️:   else if(m <= (traits::const_max-c)/a) {  // i.e. a*m+c <= max
        // Boost❗✔️:     IntType suppress_warnings = (m == 0);
        // Boost❗✔️:     BOOST_ASSERT(suppress_warnings == 0);
        // Boost❗✔️:     return (a*x+c) % (m + suppress_warnings);
        // Boost❗✔️:   } else
        // Boost❗✔️:     return add(mult(a, x), c);
        // Boost❗✔️: }
        // END BOOST RANDOM HELPER const_mod::mult_add
        // BEGIN BOOST RANDOM HELPER const_mod::add
        // Boost❗✔️: static IntType add(IntType x, IntType c)
        // Boost❗✔️: {
        // Boost❗✔️:   if(((unsigned_m() - 1) & unsigned_m()) == 0)
        // Boost❗✔️:     return (unsigned_type(x) + unsigned_type(c)) & (unsigned_m() - 1);
        // Boost❗✔️:   else if(c == 0)
        // Boost❗✔️:     return x;
        // Boost❗✔️:   else if(x < m - c)
        // Boost❗✔️:     return x + c;
        // Boost❗✔️:   else
        // Boost❗✔️:     return x - (m - c);
        // Boost❗✔️: }
        // END BOOST RANDOM HELPER const_mod::add
        // BEGIN BOOST RANDOM HELPER const_mod::mult
        // Boost❗✔️: static IntType mult(IntType a, IntType x)
        // Boost❗✔️: {
        // Boost❗✔️:   if(((unsigned_m() - 1) & unsigned_m()) == 0)
        // Boost❗✔️:     return unsigned_type(a) * unsigned_type(x) & (unsigned_m() - 1);
        // Boost❗✔️:   else if(a == 0)
        // Boost❗✔️:     return 0;
        // Boost❗✔️:   else if(a == 1)
        // Boost❗✔️:     return x;
        // Boost❗✔️:   else if(m <= traits::const_max/a)      // i.e. a*m <= max
        // Boost❗✔️:     return mult_small(a, x);
        // Boost❗✔️:   else if(traits::is_signed && (m%a < m/a))
        // Boost❗✔️:     return mult_schrage(a, x);
        // Boost❗✔️:   else
        // Boost❗✔️:     return mult_general(a, x);
        // Boost❗✔️: }
        // END BOOST RANDOM HELPER const_mod::mult
        // BEGIN BOOST RANDOM HELPER const_mod::mult_general
        // Boost❗✔️: static IntType mult_general(IntType a, IntType b)
        // Boost❗✔️: {
        // Boost❗✔️:   IntType suppress_warnings = (m == 0);
        // Boost❗✔️:   BOOST_ASSERT(suppress_warnings == 0);
        // Boost❗✔️:   IntType modulus = m + suppress_warnings;
        // Boost❗✔️:   BOOST_ASSERT(modulus == m);
        // Boost❗✔️:   if(::boost::uintmax_t(modulus) <=
        // Boost❗✔️:       (::std::numeric_limits< ::boost::uintmax_t>::max)() / modulus)
        // Boost❗✔️:   {
        // Boost❗✔️:     return static_cast<IntType>(boost::uintmax_t(a) * b % modulus);
        // Boost❗✔️:   } else {
        // Boost❗✔️:     return static_cast<IntType>(detail::mulmod(a, b, modulus));
        // Boost❗✔️:   }
        // Boost❗✔️: }
        // END BOOST RANDOM HELPER const_mod::mult_general
        // Behavior: for this u32 specialization, mult_add selects add(mult, 0),
        // mult selects mult_general, and its u64 product branch is proven safe
        // by the fixed source bounds recorded in the lane receipt.
        // Complexity: one widened multiply, fixed-width remainder, assignment
        // and return; O(1), allocation-free, with no per-draw locking.
        let next = (u64::from(MULTIPLIER) * u64::from(self.state)) % u64::from(MODULUS);
        self.state = next as u32;
        self.state
    }
    /// Draw from Boost's integral-engine uniform [0,1) distribution.
    pub fn next_unit_f64(&mut self) -> f64 {
        // BEGIN BOOST CPP FUNCTION boost::random::detail::generate_uniform_real integral path (uniform_real_distribution.hpp)
        // RDKit❗✔️: typedef typename Engine::result_type base_result;
        // RDKit❗✔️: result_type numerator = static_cast<T>(subtract<base_result>()(eng(), (eng.min)()));
        // RDKit❗✔️: result_type divisor = static_cast<T>(subtract<base_result>()((eng.max)(), (eng.min)())) + 1;
        // RDKit❗✔️: T result = numerator / divisor * (max_value - min_value) + min_value;
        // RDKit❗✔️: if(result < max_value) return result;
        // END BOOST CPP FUNCTION boost::random::detail::generate_uniform_real integral path
        // Fixed source minstd_rand bounds are 1..=2147483646: subtracting
        // the minimum gives 0..=2147483645 and divisor2147483646. Every
        // result is below1, so the source retry is unreachable. One draw and
        // fixed scalar subtraction/division, O(1), no allocation or locking.
        (f64::from(self.next_u32()) - 1.0) / (f64::from(MODULUS) - 1.0)
    }
}

/// A borrowed draw handle for one outer RDKit-compatible random operation.
pub struct RdkitRandomGenerator<'a> {
    engine: &'a mut RdkitRandomEngine,
}

impl RdkitRandomGenerator<'_> {
    pub fn next_unit_f64(&mut self) -> f64 {
        self.engine.next_unit_f64()
    }

    /// Advance the shared minstd_rand stream by one value.
    pub fn next_u32(&mut self) -> u32 {
        self.engine.next_u32()
    }
}

/// Borrow the process random stream for one outer operation.
///
/// The closure must not re-enter this function while the borrow is active.
pub fn with_rdkit_random_generator<T>(
    seed: i32,
    f: impl FnOnce(&mut RdkitRandomGenerator<'_>) -> T,
) -> T {
    // BEGIN RDKIT RANDOM FUNCTION getRandomGenerator process lifetime and seed gate
    // RDKit❗❌: static rng_type generator(42u);
    // RDKit❗❌: rng_type &getRandomGenerator(int seed) {
    // RDKit❗❌:   if (seed > 0) {
    // RDKit❗❌:     generator.seed(seed);
    // RDKit❗❌:   }
    // RDKit❗❌:   return generator;
    // RDKit❗❌: }
    // END RDKIT RANDOM FUNCTION getRandomGenerator process lifetime and seed gate
    // Behavior: one process-owned minstd_rand begins at seed42. Each call
    // optionally reseeds only for a positive i32 before lending the same state
    // to the entire outer operation; returned errors and panics do not roll it
    // back. Poison recovery keeps the valid advanced engine state.
    // Complexity: one Mutex acquisition per outer call (extra O(1) sync cost
    // versus the unsynchronized C++ global); draws take no lock or allocation.
    let mut engine = RDKIT_RANDOM_GENERATOR
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    if seed > 0 {
        engine.seed(seed as u32);
    }
    let mut generator = RdkitRandomGenerator {
        engine: &mut engine,
    };
    f(&mut generator)
}

#[cfg(test)]
mod tests {
    use super::{MODULUS, RdkitRandomEngine, with_rdkit_random_generator};

    fn first_eight(seed: u32) -> [u32; 8] {
        let mut engine = RdkitRandomEngine::from_seed(seed);
        std::array::from_fn(|_| engine.next_u32())
    }

    #[test]
    fn rdkit_rng_engine_seed_42_matches_pinned_boost_sequence() {
        assert_eq!(
            first_eight(42),
            [
                2_027_382,
                1_226_992_407,
                551_494_037,
                961_371_815,
                1_404_753_842,
                2_076_553_157,
                1_350_734_175,
                1_538_354_858,
            ]
        );
    }

    #[test]
    fn rdkit_rng_engine_seed_1_matches_pinned_boost_sequence() {
        assert_eq!(
            first_eight(1),
            [
                48_271,
                182_605_794,
                1_291_394_886,
                1_914_720_637,
                2_078_669_041,
                407_355_683,
                1_105_902_161,
                854_716_505,
            ]
        );
    }

    #[test]
    fn rdkit_rng_engine_seed_modulus_minus_one_matches_pinned_boost_sequence() {
        assert_eq!(
            first_eight(MODULUS - 1),
            [
                2_147_435_376,
                1_964_877_853,
                856_088_761,
                232_763_010,
                68_814_606,
                1_740_127_964,
                1_041_581_486,
                1_292_767_142,
            ]
        );
    }

    #[test]
    fn rdkit_rng_engine_zero_and_modulus_seeds_normalize_to_one() {
        let expected = first_eight(1);
        assert_eq!(first_eight(0), expected);
        assert_eq!(first_eight(MODULUS), expected);
    }

    #[test]
    fn rdkit_rng_engine_repeated_seeds_repeat_and_draws_advance() {
        let expected = first_eight(42);
        assert_eq!(first_eight(42), expected);

        let mut engine = RdkitRandomEngine::from_seed(42);
        for output in expected {
            let next = engine.next_u32();
            assert_eq!(next, output);
            assert!((1..MODULUS).contains(&next));
        }
        assert_eq!(engine.next_u32(), 90_320_905);
    }

    #[test]
    fn rdkit_rng_engine_reaches_both_source_range_endpoints() {
        assert_eq!(RdkitRandomEngine::from_seed(1_899_818_559).next_u32(), 1);
        assert_eq!(
            RdkitRandomEngine::from_seed(247_665_088).next_u32(),
            MODULUS - 1
        );
    }

    #[test]
    fn rdkit_rng_stream_seed_acquisition_error_and_unwind_continuity() {
        use std::panic::{AssertUnwindSafe, catch_unwind};

        // This is the only test that touches the process stream, so the
        // initial-seed assertion cannot depend on parallel test scheduling.
        assert_eq!(
            with_rdkit_random_generator(-1, |rng| rng.next_u32()),
            2_027_382
        );
        assert_eq!(
            with_rdkit_random_generator(0, |rng| rng.next_u32()),
            1_226_992_407
        );
        assert_eq!(
            with_rdkit_random_generator(-9, |rng| rng.next_u32()),
            551_494_037
        );
        assert_eq!(
            with_rdkit_random_generator(i32::MIN, |rng| rng.next_u32()),
            961_371_815
        );

        // A no-draw closure leaves the process state untouched.
        with_rdkit_random_generator(0, |_| ());
        assert_eq!(
            with_rdkit_random_generator(0, |rng| rng.next_u32()),
            1_404_753_842
        );

        // Positive reseeds consume no draw; i32::MAX is the modulus and
        // therefore normalizes to the same engine state as seed one.
        with_rdkit_random_generator(1, |_| ());
        assert_eq!(with_rdkit_random_generator(0, |rng| rng.next_u32()), 48_271);
        assert_eq!(
            with_rdkit_random_generator(i32::MAX, |rng| rng.next_u32()),
            48_271
        );

        // Returning an error preserves draws already consumed by the closure.
        let returned_error = with_rdkit_random_generator(42, |rng| {
            assert_eq!(rng.next_u32(), 2_027_382);
            Err::<(), _>("retained error")
        });
        assert_eq!(returned_error, Err("retained error"));
        assert_eq!(
            with_rdkit_random_generator(0, |rng| rng.next_u32()),
            1_226_992_407
        );

        // Unwinding poisons the lock but does not reset the valid advanced
        // engine; later acquisitions recover and continue the same sequence.
        let panic_result = catch_unwind(AssertUnwindSafe(|| {
            with_rdkit_random_generator(42, |rng| {
                assert_eq!(rng.next_u32(), 2_027_382);
                panic!("intentional RNG continuity regression");
            });
        }));
        assert!(panic_result.is_err());
        assert_eq!(
            with_rdkit_random_generator(-1, |rng| rng.next_u32()),
            1_226_992_407
        );
        assert_eq!(
            with_rdkit_random_generator(-1, |rng| rng.next_u32()),
            551_494_037
        );
    }
}
