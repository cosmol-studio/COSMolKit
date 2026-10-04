const STATE_WORDS: usize = 4;
const SHIFT_SIZE: usize = 2;
const REDUNDANT_BITS: usize = 31;
const TWIST_MASK: u32 = 0x9908_b0df;
const TEMPERING_B: u32 = 0x9d2c_5680;
const TEMPERING_C: u32 = 0xefc6_0000;
const INITIALIZATION_MULTIPLIER: u32 = 1_812_433_253;
const DEFAULT_SEED: u32 = 5_489;

/// The four-word Boost MT specialization used by RDKit's fingerprint writer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct BoostFingerprintRng {
    state: [u32; STATE_WORDS],
    index: usize,
}

impl BoostFingerprintRng {
    pub(crate) fn new(seed: u32) -> Self {
        // BEGIN BOOST RANDOM FUNCTION mersenne_twister arithmetic constructor [Boost.Random 1.85.0@09594050246adcf944be2b2c8dd2b9901083de79; include/boost/random/mersenne_twister.hpp:112-114,650-651; include/boost/random/detail/seed.hpp:20,54-55,63,103-105,111]
        // Boost❗🔝:     BOOST_RANDOM_DETAIL_ARITHMETIC_CONSTRUCTOR(mersenne_twister_engine,
        // Boost❗🔝:                                                UIntType, value)
        // Boost❗🔝:     { seed(value); }
        // Boost❗🔝:     BOOST_RANDOM_DETAIL_ARITHMETIC_CONSTRUCTOR(mersenne_twister, UIntType, val)
        // Boost❗🔝:     { seed(val); }
        // Boost❗🔝: #if !defined(BOOST_NO_SFINAE) && !defined(__SUNPRO_CC) && !defined(BOOST_BORLANDC)
        // Boost❗🔝: #define BOOST_RANDOM_DETAIL_ARITHMETIC_CONSTRUCTOR(Self, T, x)  \
        // Boost❗🔝:     explicit Self(const T& x)
        // Boost❗🔝: #else
        // Boost❗🔝: #define BOOST_RANDOM_DETAIL_ARITHMETIC_CONSTRUCTOR(Self, T, x)  \
        // Boost❗🔝:     explicit Self(const T& x) { boost_random_constructor_impl(x, ::boost::true_type()); }\
        // Boost❗🔝:     void boost_random_constructor_impl(const T& x, ::boost::true_type)
        // Boost❗🔝: #endif
        // END BOOST RANDOM FUNCTION mersenne_twister arithmetic constructor
        // The pinned deprecated wrapper binds d to all-one UIntType and f to
        // 1812433253; its final template argument v is not forwarded. The
        // fixed array keeps the same four-word state without heap allocation.
        let mut engine = Self {
            state: [0; STATE_WORDS],
            index: 0,
        };
        engine.seed(seed);
        engine
    }

    pub(crate) fn seed(&mut self, value: u32) {
        // BEGIN BOOST RANDOM FUNCTION mersenne_twister::seed [Boost.Random@09594050246adcf944be2b2c8dd2b9901083de79; include/boost/random/mersenne_twister.hpp:660-661]
        // Boost❗🔝:     BOOST_RANDOM_DETAIL_ARITHMETIC_SEED(mersenne_twister, UIntType, val)
        // Boost❗🔝:     { base_type::seed(val); }
        // END BOOST RANDOM FUNCTION mersenne_twister::seed
        // BEGIN BOOST RANDOM FUNCTION mersenne_twister_engine::seed [Boost.Random@09594050246adcf944be2b2c8dd2b9901083de79; include/boost/random/detail/seed.hpp:20,57-58,63,107-109,111; include/boost/random/mersenne_twister.hpp:141-156]
        // Boost❗🔝: #if !defined(BOOST_NO_SFINAE) && !defined(__SUNPRO_CC) && !defined(BOOST_BORLANDC)
        // Boost❗🔝: #define BOOST_RANDOM_DETAIL_ARITHMETIC_SEED(Self, T, x) \
        // Boost❗🔝:     void seed(const T& x)
        // Boost❗🔝: #else
        // Boost❗🔝: #define BOOST_RANDOM_DETAIL_ARITHMETIC_SEED(Self, T, x) \
        // Boost❗🔝:     void seed(const T& x) { boost_random_seed_impl(x, ::boost::true_type()); }\
        // Boost❗🔝:     void boost_random_seed_impl(const T& x, ::boost::true_type)
        // Boost❗🔝: #endif
        // Boost❗🔝:     BOOST_RANDOM_DETAIL_ARITHMETIC_SEED(mersenne_twister_engine, UIntType, value)
        // Boost❗🔝:     {
        // Boost❗🔝:         // New seeding algorithm from
        // Boost❗🔝:         // http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html
        // Boost❗🔝:         // In the previous versions, MSBs of the seed affected only MSBs of the
        // Boost❗🔝:         // state x[].
        // Boost❗🔝:         const UIntType mask = (max)();
        // Boost❗🔝:         x[0] = value & mask;
        // Boost❗🔝:         for (i = 1; i < n; i++) {
        // Boost❗🔝:             // See Knuth "The Art of Computer Programming"
        // Boost❗🔝:             // Vol. 2, 3rd ed., page 106
        // Boost❗🔝:             x[i] = (f * (x[i-1] ^ (x[i-1] >> (w-2))) + i) & mask;
        // Boost❗🔝:         }
        // Boost❗🔝:
        // Boost❗🔝:         normalize_state();
        // Boost❗🔝:     }
        // END BOOST RANDOM FUNCTION mersenne_twister_engine::seed
        let mask = u32::MAX;
        self.state[0] = value & mask;
        let mut index = 1;
        while index < STATE_WORDS {
            let previous = self.state[index - 1];
            let folded = previous ^ (previous >> 30);
            self.state[index] = INITIALIZATION_MULTIPLIER
                .wrapping_mul(folded)
                .wrapping_add(index as u32)
                & mask;
            index += 1;
        }
        self.index = index;
        self.normalize_state();
    }

    fn normalize_state(&mut self) {
        // BEGIN BOOST RANDOM FUNCTION mersenne_twister_engine::normalize_state [Boost.Random 1.85.0@09594050246adcf944be2b2c8dd2b9901083de79; include/boost/random/mersenne_twister.hpp:346-363]
        // Boost❗🔝:     void normalize_state()
        // Boost❗🔝:     {
        // Boost❗🔝:         const UIntType upper_mask = (~static_cast<UIntType>(0)) << r;
        // Boost❗🔝:         const UIntType lower_mask = ~upper_mask;
        // Boost❗🔝:         UIntType y0 = x[m-1] ^ x[n-1];
        // Boost❗🔝:         if(y0 & (static_cast<UIntType>(1) << (w-1))) {
        // Boost❗🔝:             y0 = ((y0 ^ a) << 1) | 1;
        // Boost❗🔝:         } else {
        // Boost❗🔝:             y0 = y0 << 1;
        // Boost❗🔝:         }
        // Boost❗🔝:         x[0] = (x[0] & upper_mask) | (y0 & lower_mask);
        // Boost❗🔝:
        // Boost❗🔝:         // fix up the state if it's all zeroes.
        // Boost❗🔝:         for(std::size_t j = 0; j < n; ++j) {
        // Boost❗🔝:             if(x[j] != 0) return;
        // Boost❗🔝:         }
        // Boost❗🔝:         x[0] = static_cast<UIntType>(1) << (w-1);
        // Boost❗🔝:     }
        // END BOOST RANDOM FUNCTION mersenne_twister_engine::normalize_state
        let upper_mask = (!0_u32) << REDUNDANT_BITS;
        let lower_mask = !upper_mask;
        let mut y0 = self.state[SHIFT_SIZE - 1] ^ self.state[STATE_WORDS - 1];
        if y0 & (1_u32 << 31) != 0 {
            y0 = ((y0 ^ TWIST_MASK).wrapping_shl(1)) | 1;
        } else {
            y0 = y0.wrapping_shl(1);
        }
        self.state[0] = (self.state[0] & upper_mask) | (y0 & lower_mask);
        for word in &self.state {
            if *word != 0 {
                return;
            }
        }
        self.state[0] = 1_u32 << 31;
    }

    fn twist(&mut self) {
        // BEGIN BOOST RANDOM FUNCTION mersenne_twister_engine::twist [Boost.Random 1.85.0@09594050246adcf944be2b2c8dd2b9901083de79; include/boost/random/mersenne_twister.hpp:530-574]
        // Boost❗🔝: template<class UIntType,
        // Boost❗🔝:          std::size_t w, std::size_t n, std::size_t m, std::size_t r,
        // Boost❗🔝:          UIntType a, std::size_t u, UIntType d, std::size_t s,
        // Boost❗🔝:          UIntType b, std::size_t t,
        // Boost❗🔝:          UIntType c, std::size_t l, UIntType f>
        // Boost❗🔝: void
        // Boost❗🔝: mersenne_twister_engine<UIntType,w,n,m,r,a,u,d,s,b,t,c,l,f>::twist()
        // Boost❗🔝: {
        // Boost❗🔝:     const UIntType upper_mask = (~static_cast<UIntType>(0)) << r;
        // Boost❗🔝:     const UIntType lower_mask = ~upper_mask;
        // Boost❗🔝:
        // Boost❗🔝:     const std::size_t unroll_factor = 6;
        // Boost❗🔝:     const std::size_t unroll_extra1 = (n-m) % unroll_factor;
        // Boost❗🔝:     const std::size_t unroll_extra2 = (m-1) % unroll_factor;
        // Boost❗🔝:
        // Boost❗🔝:     // split loop to avoid costly modulo operations
        // Boost❗🔝:     {  // extra scope for MSVC brokenness w.r.t. for scope
        // Boost❗🔝:         for(std::size_t j = 0; j < n-m-unroll_extra1; j++) {
        // Boost❗🔝:             UIntType y = (x[j] & upper_mask) | (x[j+1] & lower_mask);
        // Boost❗🔝:             x[j] = x[j+m] ^ (y >> 1) ^ ((x[j+1]&1) * a);
        // Boost❗🔝:         }
        // Boost❗🔝:     }
        // Boost❗🔝:     {
        // Boost❗🔝:         for(std::size_t j = n-m-unroll_extra1; j < n-m; j++) {
        // Boost❗🔝:             UIntType y = (x[j] & upper_mask) | (x[j+1] & lower_mask);
        // Boost❗🔝:             x[j] = x[j+m] ^ (y >> 1) ^ ((x[j+1]&1) * a);
        // Boost❗🔝:         }
        // Boost❗🔝:     }
        // Boost❗🔝:     {
        // Boost❗🔝:         for(std::size_t j = n-m; j < n-1-unroll_extra2; j++) {
        // Boost❗🔝:             UIntType y = (x[j] & upper_mask) | (x[j+1] & lower_mask);
        // Boost❗🔝:             x[j] = x[j-(n-m)] ^ (y >> 1) ^ ((x[j+1]&1) * a);
        // Boost❗🔝:         }
        // Boost❗🔝:     }
        // Boost❗🔝:     {
        // Boost❗🔝:         for(std::size_t j = n-1-unroll_extra2; j < n-1; j++) {
        // Boost❗🔝:             UIntType y = (x[j] & upper_mask) | (x[j+1] & lower_mask);
        // Boost❗🔝:             x[j] = x[j-(n-m)] ^ (y >> 1) ^ ((x[j+1]&1) * a);
        // Boost❗🔝:         }
        // Boost❗🔝:     }
        // Boost❗🔝:     // last iteration
        // Boost❗🔝:     UIntType y = (x[n-1] & upper_mask) | (x[0] & lower_mask);
        // Boost❗🔝:     x[n-1] = x[m-1] ^ (y >> 1) ^ ((x[0]&1) * a);
        // Boost❗🔝:     i = 0;
        // Boost❗🔝: }
        // END BOOST RANDOM FUNCTION mersenne_twister_engine::twist
        let upper_mask = (!0_u32) << REDUNDANT_BITS;
        let lower_mask = !upper_mask;
        let unroll_factor = 6;
        let unroll_extra1 = (STATE_WORDS - SHIFT_SIZE) % unroll_factor;
        let unroll_extra2 = (SHIFT_SIZE - 1) % unroll_factor;

        for index in 0..STATE_WORDS - SHIFT_SIZE - unroll_extra1 {
            let y = (self.state[index] & upper_mask) | (self.state[index + 1] & lower_mask);
            self.state[index] = self.state[index + SHIFT_SIZE]
                ^ (y >> 1)
                ^ ((self.state[index + 1] & 1) * TWIST_MASK);
        }
        for index in STATE_WORDS - SHIFT_SIZE - unroll_extra1..STATE_WORDS - SHIFT_SIZE {
            let y = (self.state[index] & upper_mask) | (self.state[index + 1] & lower_mask);
            self.state[index] = self.state[index + SHIFT_SIZE]
                ^ (y >> 1)
                ^ ((self.state[index + 1] & 1) * TWIST_MASK);
        }
        for index in STATE_WORDS - SHIFT_SIZE..STATE_WORDS - 1 - unroll_extra2 {
            let y = (self.state[index] & upper_mask) | (self.state[index + 1] & lower_mask);
            self.state[index] = self.state[index - (STATE_WORDS - SHIFT_SIZE)]
                ^ (y >> 1)
                ^ ((self.state[index + 1] & 1) * TWIST_MASK);
        }
        for index in STATE_WORDS - 1 - unroll_extra2..STATE_WORDS - 1 {
            let y = (self.state[index] & upper_mask) | (self.state[index + 1] & lower_mask);
            self.state[index] = self.state[index - (STATE_WORDS - SHIFT_SIZE)]
                ^ (y >> 1)
                ^ ((self.state[index + 1] & 1) * TWIST_MASK);
        }
        let y = (self.state[STATE_WORDS - 1] & upper_mask) | (self.state[0] & lower_mask);
        self.state[STATE_WORDS - 1] =
            self.state[SHIFT_SIZE - 1] ^ (y >> 1) ^ ((self.state[0] & 1) * TWIST_MASK);
        self.index = 0;
    }

    pub(crate) fn next_u32(&mut self) -> u32 {
        // BEGIN BOOST RANDOM FUNCTION mersenne_twister_engine::operator() [Boost.Random 1.85.0@09594050246adcf944be2b2c8dd2b9901083de79; include/boost/random/mersenne_twister.hpp:577-596]
        // Boost❗🔝: template<class UIntType,
        // Boost❗🔝:          std::size_t w, std::size_t n, std::size_t m, std::size_t r,
        // Boost❗🔝:          UIntType a, std::size_t u, UIntType d, std::size_t s,
        // Boost❗🔝:          UIntType b, std::size_t t,
        // Boost❗🔝:          UIntType c, std::size_t l, UIntType f>
        // Boost❗🔝: inline typename
        // Boost❗🔝: mersenne_twister_engine<UIntType,w,n,m,r,a,u,d,s,b,t,c,l,f>::result_type
        // Boost❗🔝: mersenne_twister_engine<UIntType,w,n,m,r,a,u,d,s,b,t,c,l,f>::operator()()
        // Boost❗🔝: {
        // Boost❗🔝:     if(i == n)
        // Boost❗🔝:         twist();
        // Boost❗🔝:     // Step 4
        // Boost❗🔝:     UIntType z = x[i];
        // Boost❗🔝:     ++i;
        // Boost❗🔝:     z ^= ((z >> u) & d);
        // Boost❗🔝:     z ^= ((z << s) & b);
        // Boost❗🔝:     z ^= ((z << t) & c);
        // Boost❗🔝:     z ^= (z >> l);
        // Boost❗🔝:     return z;
        // Boost❗🔝: }
        // END BOOST RANDOM FUNCTION mersenne_twister_engine::operator()
        if self.index == STATE_WORDS {
            self.twist();
        }
        let mut value = self.state[self.index];
        self.index += 1;
        value ^= (value >> 11) & u32::MAX;
        value ^= value.wrapping_shl(7) & TEMPERING_B;
        value ^= value.wrapping_shl(15) & TEMPERING_C;
        value ^= value >> 18;
        value
    }
}

impl Default for BoostFingerprintRng {
    fn default() -> Self {
        // BEGIN BOOST RANDOM FUNCTION mersenne_twister_engine default constructor and seed [Boost.Random 1.85.0@09594050246adcf944be2b2c8dd2b9901083de79; include/boost/random/mersenne_twister.hpp:89,104-107,132-133,647,654]
        // Boost❗🔝:     BOOST_STATIC_CONSTANT(UIntType, default_seed = 5489u);
        // Boost❗🔝:     /**
        // Boost❗🔝:      * Constructs a @c mersenne_twister_engine and calls @c seed().
        // Boost❗🔝:      */
        // Boost❗🔝:     mersenne_twister_engine() { seed(); }
        // Boost❗🔝:     /** Calls @c seed(default_seed). */
        // Boost❗🔝:     void seed() { seed(default_seed); }
        // Boost❗🔝:     mersenne_twister() {}
        // Boost❗🔝:     void seed() { base_type::seed(); }
        // END BOOST RANDOM FUNCTION mersenne_twister_engine default constructor and seed
        Self::new(DEFAULT_SEED)
    }
}

/// The integral engine surface consumed by Boost's uniform-int distribution.
/// The production implementation is the pinned four-word fingerprint MT;
/// tests may supply finite literal sequences to exercise generic branches.
pub(super) trait BoostUniformIntEngine {
    fn range_min(&self) -> u32;
    fn range_max(&self) -> u32;
    fn next_word(&mut self) -> u32;
}

impl BoostUniformIntEngine for BoostFingerprintRng {
    fn range_min(&self) -> u32 {
        0
    }

    fn range_max(&self) -> u32 {
        u32::MAX
    }

    fn next_word(&mut self) -> u32 {
        self.next_u32()
    }
}

/// Stateless inclusive integer distribution used for additional Morgan bits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct BoostUniformIntDistribution {
    min_value: i32,
    max_value: i32,
}

impl BoostUniformIntDistribution {
    fn new(min_value: i32, max_value: i32) -> Self {
        // BEGIN BOOST RANDOM FUNCTION uniform_int constructor [Boost.Random 1.85.0@09594050246adcf944be2b2c8dd2b9901083de79; include/boost/random/uniform_int.hpp:55-63]
        // Boost❗✔️:     /**
        // Boost❗✔️:      * Constructs a uniform_int object. @c min and @c max are
        // Boost❗✔️:      * the parameters of the distribution.
        // Boost❗✔️:      *
        // Boost❗✔️:      * Requires: min <= max
        // Boost❗✔️:      */
        // Boost❗✔️:     explicit uniform_int(IntType min_arg = 0, IntType max_arg = 9)
        // Boost❗✔️:       : base_type(min_arg, max_arg)
        // Boost❗✔️:     {}
        // END BOOST RANDOM FUNCTION uniform_int constructor
        // BEGIN BOOST RANDOM FUNCTION uniform_int_distribution constructor [Boost.Random 1.85.0@09594050246adcf944be2b2c8dd2b9901083de79; include/boost/random/uniform_int_distribution.hpp:326-338]
        // Boost❗✔️:     /**
        // Boost❗✔️:      * Constructs a uniform_int_distribution. @c min and @c max are
        // Boost❗✔️:      * the parameters of the distribution.
        // Boost❗✔️:      *
        // Boost❗✔️:      * Requires: min <= max
        // Boost❗✔️:      */
        // Boost❗✔️:     explicit uniform_int_distribution(
        // Boost❗✔️:         IntType min_arg = 0,
        // Boost❗✔️:         IntType max_arg = (std::numeric_limits<IntType>::max)())
        // Boost❗✔️:       : _min(min_arg), _max(max_arg)
        // Boost❗✔️:     {
        // Boost❗✔️:         BOOST_ASSERT(min_arg <= max_arg);
        // Boost❗✔️:     }
        // END BOOST RANDOM FUNCTION uniform_int_distribution constructor
        debug_assert!(min_value <= max_value);
        Self {
            min_value,
            max_value,
        }
    }

    pub(super) fn rdkit_additional_bit() -> Self {
        // RDKit explicitly constructs boost::uniform_int<>(0, INT_MAX),
        // overriding the deprecated wrapper's default upper endpoint of 9.
        Self::new(0, i32::MAX)
    }

    pub(super) fn sample<Engine: BoostUniformIntEngine>(&self, engine: &mut Engine) -> i32 {
        // BEGIN BOOST RANDOM FUNCTION uniform_int_distribution::operator() [Boost.Random 1.85.0@09594050246adcf944be2b2c8dd2b9901083de79; include/boost/random/uniform_int.hpp:77-81; include/boost/random/uniform_int_distribution.hpp:368-371]
        // Boost❗✔️:     template<class Engine>
        // Boost❗✔️:     IntType operator()(Engine& eng) const
        // Boost❗✔️:     {
        // Boost❗✔️:         return static_cast<const base_type&>(*this)(eng);
        // Boost❗✔️:     }
        // Boost❗✔️:     /** Returns an integer uniformly distributed in the range [min, max]. */
        // Boost❗✔️:     template<class Engine>
        // Boost❗✔️:     result_type operator()(Engine& eng) const
        // Boost❗✔️:     { return detail::generate_uniform_int(eng, _min, _max); }
        // END BOOST RANDOM FUNCTION uniform_int_distribution::operator()
        // BEGIN BOOST RANDOM FUNCTION variate_generator::operator()/engine [Boost.Random 1.85.0@09594050246adcf944be2b2c8dd2b9901083de79; include/boost/random/variate_generator.hpp:69-73,83-87; include/boost/random/detail/ptr_helper.hpp:35-43]
        // Boost❗✔️:     variate_generator(Engine e, Distribution d)
        // Boost❗✔️:       : _eng(e), _dist(d) { }
        // Boost❗✔️:
        // Boost❗✔️:     /** Returns: distribution()(engine()) */
        // Boost❗✔️:     result_type operator()() { return _dist(engine()); }
        // Boost❗✔️:     engine_value_type& engine() { return helper_type::ref(_eng); }
        // Boost❗✔️:     /**
        // Boost❗✔️:      * Returns: A reference to the associated uniform random number generator.
        // Boost❗✔️:      */
        // Boost❗✔️:     const engine_value_type& engine() const { return helper_type::ref(_eng); }
        // Boost❗✔️: template<class T>
        // Boost❗✔️: struct ptr_helper<T&>
        // Boost❗✔️: {
        // Boost❗✔️:   typedef T value_type;
        // Boost❗✔️:   typedef T& reference_type;
        // Boost❗✔️:   typedef T& rvalue_type;
        // Boost❗✔️:   static reference_type ref(T& r) { return r; }
        // Boost❗✔️:   static const T& ref(const T& r) { return r; }
        // Boost❗✔️: };
        // END BOOST RANDOM FUNCTION variate_generator::operator()/engine
        // The pinned RDKit caller at FingerprintGenerator.cpp:383 instantiates
        // variate_generator<rng_type &, distrib_type>, selecting ptr_helper<T&>
        // so each draw advances the same borrowed engine.
        // The source engine is borrowed and the distribution is stateless; this
        // call advances the same engine object held for the outer operation.
        let range = signed_i32_range(self.min_value, self.max_value);
        let offset = generate_uniform_u32(engine, range);
        add_u32_offset_to_i32(self.min_value, offset)
    }
}

fn signed_i32_range(min_value: i32, max_value: i32) -> u32 {
    // BEGIN BOOST RANDOM FUNCTION signed_unsigned_tools::subtract<signed> [Boost.Random 1.85.0@09594050246adcf944be2b2c8dd2b9901083de79; include/boost/random/detail/signed_unsigned_tools.hpp:37-51]
    // Boost❗✔️: template<class T>
    // Boost❗✔️: struct subtract<T, /* signed */ true>
    // Boost❗✔️: {
    // Boost❗✔️:   typedef typename boost::random::traits::make_unsigned_or_unbounded<T>::type result_type;
    // Boost❗✔️:   result_type operator()(T x, T y)
    // Boost❗✔️:   {
    // Boost❗✔️:     if (y >= 0)   // because x >= y, it follows that x >= 0, too
    // Boost❗✔️:       return result_type(x) - result_type(y);
    // Boost❗✔️:     if (x >= 0)   // y < 0
    // Boost❗✔️:       // avoid the nasty two's complement case for y == min()
    // Boost❗✔️:       return result_type(x) + result_type(-(y+1)) + 1;
    // Boost❗✔️:     // both x and y are negative: no signed overflow
    // Boost❗✔️:     return result_type(x - y);
    // Boost❗✔️:   }
    // Boost❗✔️: };
    // END BOOST RANDOM FUNCTION signed_unsigned_tools::subtract<signed>
    // A widened subtraction has the same unsigned result for every valid i32
    // interval, including i32::MIN..=i32::MAX, without signed overflow.
    debug_assert!(min_value <= max_value);
    (i64::from(max_value) - i64::from(min_value)) as u32
}

fn add_u32_offset_to_i32(min_value: i32, offset: u32) -> i32 {
    // BEGIN BOOST RANDOM FUNCTION signed_unsigned_tools::add<unsigned,signed> [Boost.Random 1.85.0@09594050246adcf944be2b2c8dd2b9901083de79; include/boost/random/detail/signed_unsigned_tools.hpp:67-82]
    // Boost❗✔️: template<class T1, class T2>
    // Boost❗✔️: struct add<T1, T2, /* signed */ true>
    // Boost❗✔️: {
    // Boost❗✔️:   typedef T2 result_type;
    // Boost❗✔️:   result_type operator()(T1 x, T2 y)
    // Boost❗✔️:   {
    // Boost❗✔️:     if (y >= 0)
    // Boost❗✔️:       return T2(x) + y;
    // Boost❗✔️:     // y < 0
    // Boost❗✔️:     if (x > T1(-(y+1)))  // result >= 0 after subtraction
    // Boost❗✔️:       // avoid the nasty two's complement edge case for y == min()
    // Boost❗✔️:       return T2(x - T1(-(y+1)) - 1);
    // Boost❗✔️:     // abs(x) < abs(y), thus T2 able to represent x
    // Boost❗✔️:     return T2(x) + y;
    // Boost❗✔️:   }
    // Boost❗✔️: };
    // END BOOST RANDOM FUNCTION signed_unsigned_tools::add<unsigned,signed>
    // Source-valid offsets keep the mathematical sum inside i32; i64 makes
    // the source's signed-minimum avoidance explicit without changing it.
    (i64::from(min_value) + i64::from(offset)) as i32
}

fn generate_uniform_u32<Engine: BoostUniformIntEngine>(engine: &mut Engine, range: u32) -> u32 {
    // Boost.Random 1.85.0, detail::generate_uniform_int for an integral
    // engine and u32 range_type. The source's signed/unsigned helper overloads
    // are represented by signed_i32_range/add_u32_offset_to_i32 above.
    // BEGIN BOOST RANDOM FUNCTION detail::generate_uniform_int integral [Boost.Random 1.85.0@09594050246adcf944be2b2c8dd2b9901083de79; include/boost/random/uniform_int_distribution.hpp:49-228]
    // Boost❗✔️: template<class Engine, class T>
    // Boost❗✔️: T generate_uniform_int(
    // Boost❗✔️:     Engine& eng, T min_value, T max_value,
    // Boost❗✔️:     boost::true_type /** is_integral<Engine::result_type> */)
    // Boost❗✔️: {
    // Boost❗✔️:     typedef T result_type;
    // Boost❗✔️:     typedef typename boost::random::traits::make_unsigned_or_unbounded<T>::type range_type;
    // Boost❗✔️:     typedef typename Engine::result_type base_result;
    // Boost❗✔️:     // ranges are always unsigned or unbounded
    // Boost❗✔️:     typedef typename boost::random::traits::make_unsigned_or_unbounded<base_result>::type base_unsigned;
    // Boost❗✔️:     const range_type range = random::detail::subtract<result_type>()(max_value, min_value);
    // Boost❗✔️:     const base_result bmin = (eng.min)();
    // Boost❗✔️:     const base_unsigned brange =
    // Boost❗✔️:       random::detail::subtract<base_result>()((eng.max)(), (eng.min)());
    // Boost❗✔️:
    // Boost❗✔️:     if(range == 0) {
    // Boost❗✔️:       return min_value;
    // Boost❗✔️:     } else if(brange == range) {
    // Boost❗✔️:       // this will probably never happen in real life
    // Boost❗✔️:       // basically nothing to do; just take care we don't overflow / underflow
    // Boost❗✔️:       base_unsigned v = random::detail::subtract<base_result>()(eng(), bmin);
    // Boost❗✔️:       return random::detail::add<base_unsigned, result_type>()(v, min_value);
    // Boost❗✔️:     } else if(brange < range) {
    // Boost❗✔️:       // use rejection method to handle things like 0..3 --> 0..4
    // Boost❗✔️:       for(;;) {
    // Boost❗✔️:         // concatenate several invocations of the base RNG
    // Boost❗✔️:         // take extra care to avoid overflows
    // Boost❗✔️:
    // Boost❗✔️:         //  limit == floor((range+1)/(brange+1))
    // Boost❗✔️:         //  Therefore limit*(brange+1) <= range+1
    // Boost❗✔️:         range_type limit;
    // Boost❗✔️:         if(range == (std::numeric_limits<range_type>::max)()) {
    // Boost❗✔️:           limit = range/(range_type(brange)+1);
    // Boost❗✔️:           if(range % (range_type(brange)+1) == range_type(brange))
    // Boost❗✔️:             ++limit;
    // Boost❗✔️:         } else {
    // Boost❗✔️:           limit = (range+1)/(range_type(brange)+1);
    // Boost❗✔️:         }
    // Boost❗✔️:
    // Boost❗✔️:         // We consider "result" as expressed to base (brange+1):
    // Boost❗✔️:         // For every power of (brange+1), we determine a random factor
    // Boost❗✔️:         range_type result = range_type(0);
    // Boost❗✔️:         range_type mult = range_type(1);
    // Boost❗✔️:
    // Boost❗✔️:         // loop invariants:
    // Boost❗✔️:         //  result < mult
    // Boost❗✔️:         //  mult <= range
    // Boost❗✔️:         while(mult <= limit) {
    // Boost❗✔️:           // Postcondition: result <= range, thus no overflow
    // Boost❗✔️:           //
    // Boost❗✔️:           // limit*(brange+1)<=range+1                   def. of limit       (1)
    // Boost❗✔️:           // eng()-bmin<=brange                          eng() post.         (2)
    // Boost❗✔️:           // and mult<=limit.                            loop condition      (3)
    // Boost❗✔️:           // Therefore mult*(eng()-bmin+1)<=range+1      by (1),(2),(3)      (4)
    // Boost❗✔️:           // Therefore mult*(eng()-bmin)+mult<=range+1   rearranging (4)     (5)
    // Boost❗✔️:           // result<mult                                 loop invariant      (6)
    // Boost❗✔️:           // Therefore result+mult*(eng()-bmin)<range+1  by (5), (6)         (7)
    // Boost❗✔️:           //
    // Boost❗✔️:           // Postcondition: result < mult*(brange+1)
    // Boost❗✔️:           //
    // Boost❗✔️:           // result<mult                                 loop invariant      (1)
    // Boost❗✔️:           // eng()-bmin<=brange                          eng() post.         (2)
    // Boost❗✔️:           // Therefore result+mult*(eng()-bmin) <
    // Boost❗✔️:           //           mult+mult*(eng()-bmin)            by (1)              (3)
    // Boost❗✔️:           // Therefore result+(eng()-bmin)*mult <
    // Boost❗✔️:           //           mult+mult*brange                  by (2), (3)         (4)
    // Boost❗✔️:           // Therefore result+(eng()-bmin)*mult <
    // Boost❗✔️:           //           mult*(brange+1)                   by (4)
    // Boost❗✔️:           result += static_cast<range_type>(static_cast<range_type>(random::detail::subtract<base_result>()(eng(), bmin)) * mult);
    // Boost❗✔️:
    // Boost❗✔️:           // equivalent to (mult * (brange+1)) == range+1, but avoids overflow.
    // Boost❗✔️:           if(mult * range_type(brange) == range - mult + 1) {
    // Boost❗✔️:               // The destination range is an integer power of
    // Boost❗✔️:               // the generator's range.
    // Boost❗✔️:               return(result);
    // Boost❗✔️:           }
    // Boost❗✔️:
    // Boost❗✔️:           // Postcondition: mult <= range
    // Boost❗✔️:           //
    // Boost❗✔️:           // limit*(brange+1)<=range+1                   def. of limit       (1)
    // Boost❗✔️:           // mult<=limit                                 loop condition      (2)
    // Boost❗✔️:           // Therefore mult*(brange+1)<=range+1          by (1), (2)         (3)
    // Boost❗✔️:           // mult*(brange+1)!=range+1                    preceding if        (4)
    // Boost❗✔️:           // Therefore mult*(brange+1)<range+1           by (3), (4)         (5)
    // Boost❗✔️:           //
    // Boost❗✔️:           // Postcondition: result < mult
    // Boost❗✔️:           //
    // Boost❗✔️:           // See the second postcondition on the change to result.
    // Boost❗✔️:           mult *= range_type(brange)+range_type(1);
    // Boost❗✔️:         }
    // Boost❗✔️:         // loop postcondition: range/mult < brange+1
    // Boost❗✔️:         //
    // Boost❗✔️:         // mult > limit                                  loop condition      (1)
    // Boost❗✔️:         // Suppose range/mult >= brange+1                Assumption          (2)
    // Boost❗✔️:         // range >= mult*(brange+1)                      by (2)              (3)
    // Boost❗✔️:         // range+1 > mult*(brange+1)                     by (3)              (4)
    // Boost❗✔️:         // range+1 > (limit+1)*(brange+1)                by (1), (4)         (5)
    // Boost❗✔️:         // (range+1)/(brange+1) > limit+1                by (5)              (6)
    // Boost❗✔️:         // limit < floor((range+1)/(brange+1))           by (6)              (7)
    // Boost❗✔️:         // limit==floor((range+1)/(brange+1))            def. of limit       (8)
    // Boost❗✔️:         // not (2)                                       reductio            (9)
    // Boost❗✔️:         //
    // Boost❗✔️:         // loop postcondition: (range/mult)*mult+(mult-1) >= range
    // Boost❗✔️:         //
    // Boost❗✔️:         // (range/mult)*mult + range%mult == range       identity            (1)
    // Boost❗✔️:         // range%mult < mult                             def. of %           (2)
    // Boost❗✔️:         // (range/mult)*mult+mult > range                by (1), (2)         (3)
    // Boost❗✔️:         // (range/mult)*mult+(mult-1) >= range           by (3)              (4)
    // Boost❗✔️:         //
    // Boost❗✔️:         // Note that the maximum value of result at this point is (mult-1),
    // Boost❗✔️:         // so after this final step, we generate numbers that can be
    // Boost❗✔️:         // at least as large as range.  We have to really careful to avoid
    // Boost❗✔️:         // overflow in this final addition and in the rejection.  Anything
    // Boost❗✔️:         // that overflows is larger than range and can thus be rejected.
    // Boost❗✔️:
    // Boost❗✔️:         // range/mult < brange+1  -> no endless loop
    // Boost❗✔️:         range_type result_increment =
    // Boost❗✔️:             generate_uniform_int(
    // Boost❗✔️:                 eng,
    // Boost❗✔️:                 static_cast<range_type>(0),
    // Boost❗✔️:                 static_cast<range_type>(range/mult),
    // Boost❗✔️:                 boost::true_type());
    // Boost❗✔️:         if(std::numeric_limits<range_type>::is_bounded && ((std::numeric_limits<range_type>::max)() / mult < result_increment)) {
    // Boost❗✔️:           // The multiplcation would overflow.  Reject immediately.
    // Boost❗✔️:           continue;
    // Boost❗✔️:         }
    // Boost❗✔️:         result_increment *= mult;
    // Boost❗✔️:         // unsigned integers are guaranteed to wrap on overflow.
    // Boost❗✔️:         result += result_increment;
    // Boost❗✔️:         if(result < result_increment) {
    // Boost❗✔️:           // The addition overflowed.  Reject.
    // Boost❗✔️:           continue;
    // Boost❗✔️:         }
    // Boost❗✔️:         if(result > range) {
    // Boost❗✔️:           // Too big.  Reject.
    // Boost❗✔️:           continue;
    // Boost❗✔️:         }
    // Boost❗✔️:         return random::detail::add<range_type, result_type>()(result, min_value);
    // Boost❗✔️:       }
    // Boost❗✔️:     } else {                   // brange > range
    // Boost❗✔️: #ifdef BOOST_NO_CXX11_EXPLICIT_CONVERSION_OPERATORS
    // Boost❗✔️:       typedef typename conditional<
    // Boost❗✔️:          std::numeric_limits<range_type>::is_specialized && std::numeric_limits<base_unsigned>::is_specialized
    // Boost❗✔️:          && (std::numeric_limits<range_type>::digits >= std::numeric_limits<base_unsigned>::digits),
    // Boost❗✔️:          range_type, base_unsigned>::type mixed_range_type;
    // Boost❗✔️: #else
    // Boost❗✔️:       typedef base_unsigned mixed_range_type;
    // Boost❗✔️: #endif
    // Boost❗✔️:
    // Boost❗✔️:       mixed_range_type bucket_size;
    // Boost❗✔️:       // it's safe to add 1 to range, as long as we cast it first,
    // Boost❗✔️:       // because we know that it is less than brange.  However,
    // Boost❗✔️:       // we do need to be careful not to cause overflow by adding 1
    // Boost❗✔️:       // to brange.  We use mixed_range_type throughout for mixed
    // Boost❗✔️:       // arithmetic between base_unsigned and range_type - in the case
    // Boost❗✔️:       // that range_type has more bits than base_unsigned it is always
    // Boost❗✔️:       // safe to use range_type for this albeit it may be more effient
    // Boost❗✔️:       // to use base_unsigned.  The latter is a narrowing conversion though
    // Boost❗✔️:       // which may be disallowed if range_type is a multiprecision type
    // Boost❗✔️:       // and there are no explicit converison operators.
    // Boost❗✔️:
    // Boost❗✔️:       if(brange == (std::numeric_limits<base_unsigned>::max)()) {
    // Boost❗✔️:         bucket_size = static_cast<mixed_range_type>(brange) / (static_cast<mixed_range_type>(range)+1);
    // Boost❗✔️:         if(static_cast<mixed_range_type>(brange) % (static_cast<mixed_range_type>(range)+1) == static_cast<mixed_range_type>(range)) {
    // Boost❗✔️:           ++bucket_size;
    // Boost❗✔️:         }
    // Boost❗✔️:       } else {
    // Boost❗✔️:         bucket_size = static_cast<mixed_range_type>(brange + 1) / (static_cast<mixed_range_type>(range)+1);
    // Boost❗✔️:       }
    // Boost❗✔️:       for(;;) {
    // Boost❗✔️:         mixed_range_type result =
    // Boost❗✔️:           random::detail::subtract<base_result>()(eng(), bmin);
    // Boost❗✔️:         result /= bucket_size;
    // Boost❗✔️:         // result and range are non-negative, and result is possibly larger
    // Boost❗✔️:         // than range, so the cast is safe
    // Boost❗✔️:         if(result <= static_cast<mixed_range_type>(range))
    // Boost❗✔️:           return random::detail::add<mixed_range_type, result_type>()(result, min_value);
    // Boost❗✔️:       }
    // Boost❗✔️:     }
    // Boost❗✔️: }
    // END BOOST RANDOM FUNCTION detail::generate_uniform_int integral
    let base_min = engine.range_min();
    let base_range = engine.range_max().wrapping_sub(base_min);

    if range == 0 {
        return 0;
    }

    if base_range == range {
        return engine.next_word().wrapping_sub(base_min);
    }

    if base_range < range {
        loop {
            let limit = if range == u32::MAX {
                let divisor = base_range + 1;
                let mut limit = range / divisor;
                if range % divisor == base_range {
                    limit += 1;
                }
                limit
            } else {
                (range + 1) / (base_range + 1)
            };

            let mut result = 0_u32;
            let mut multiplier = 1_u32;
            while multiplier <= limit {
                let digit = engine.next_word().wrapping_sub(base_min);
                result = result.wrapping_add(digit.wrapping_mul(multiplier));
                if multiplier.wrapping_mul(base_range)
                    == range.wrapping_sub(multiplier).wrapping_add(1)
                {
                    return result;
                }
                multiplier = multiplier.wrapping_mul(base_range + 1);
            }

            let mut result_increment = generate_uniform_u32(engine, range / multiplier);
            if u32::MAX / multiplier < result_increment {
                continue;
            }
            result_increment = result_increment.wrapping_mul(multiplier);
            result = result.wrapping_add(result_increment);
            if result < result_increment {
                continue;
            }
            if result > range {
                continue;
            }
            return result;
        }
    }

    // For the fixed u32 engine, range_type and mixed_range_type are u32.
    let bucket_size = if base_range == u32::MAX {
        let divisor = range + 1;
        let mut bucket_size = base_range / divisor;
        if base_range % divisor == range {
            bucket_size += 1;
        }
        bucket_size
    } else {
        (base_range + 1) / (range + 1)
    };
    loop {
        let result = engine.next_word().wrapping_sub(base_min) / bucket_size;
        if result <= range {
            return result;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct LiteralEngine<'a> {
        min: u32,
        max: u32,
        words: &'a [u32],
        index: usize,
    }

    impl<'a> LiteralEngine<'a> {
        fn new(min: u32, max: u32, words: &'a [u32]) -> Self {
            assert!(min <= max);
            Self {
                min,
                max,
                words,
                index: 0,
            }
        }

        fn consumed(&self) -> usize {
            self.index
        }
    }

    impl BoostUniformIntEngine for LiteralEngine<'_> {
        fn range_min(&self) -> u32 {
            self.min
        }

        fn range_max(&self) -> u32 {
            self.max
        }

        fn next_word(&mut self) -> u32 {
            let word = *self
                .words
                .get(self.index)
                .expect("literal engine exhausted");
            assert!((self.min..=self.max).contains(&word));
            self.index += 1;
            word
        }
    }

    fn first_eight(seed: u32) -> [u32; 8] {
        let mut engine = BoostFingerprintRng::new(seed);
        std::array::from_fn(|_| engine.next_u32())
    }

    #[test]
    fn fingerprint_morgan_r01_parameters_and_seed_recurrence() {
        assert_eq!(STATE_WORDS, 4);
        assert_eq!(SHIFT_SIZE, 2);
        assert_eq!(REDUNDANT_BITS, 31);
        assert_eq!(TWIST_MASK, 0x9908_b0df);
        assert_eq!(TEMPERING_B, 0x9d2c_5680);
        assert_eq!(TEMPERING_C, 0xefc6_0000);
        assert_eq!(INITIALIZATION_MULTIPLIER, 1_812_433_253);
        assert_eq!(DEFAULT_SEED, 5_489);

        // Literal post-seed states independently follow the pinned u32
        // recurrence and normalize_state equations for both MSB branches.
        let expected_states = [
            (0, [0x6295_9680, 0x0000_0001, 0x6c07_8967, 0x714a_cb41]),
            (1, [0x18ba_0471, 0x6c07_8966, 0xdd52_54a5, 0xb952_3b81]),
            (42, [0x228c_8b84, 0xb93c_8a93, 0x7101_4437, 0xe87a_cf51]),
            (
                DEFAULT_SEED,
                [0x7693_c9bf, 0x4d98_ee96, 0xaf25_f095, 0xafd9_ba96],
            ),
            (
                0x7fff_ffff,
                [0x0733_88a5, 0xa7f0_ed37, 0x0971_f2eb, 0x7d61_99ba],
            ),
            (
                0x8000_0000,
                [0xbb0e_9c34, 0x580f_12cb, 0x8a86_83b4, 0x0588_5cd1],
            ),
            (
                u32::MAX,
                [0xafb4_e08d, 0x4fe1_da6d, 0xeaf2_f89e, 0xc133_1af4],
            ),
        ];

        for (seed, expected_state) in expected_states {
            let engine = BoostFingerprintRng::new(seed);
            assert_eq!(engine.state, expected_state, "seed {seed:#010x}");
            assert_eq!(engine.index, STATE_WORDS, "seed {seed:#010x}");
        }
    }

    #[test]
    fn fingerprint_morgan_r01_zero_high_bit_and_all_ones_words() {
        let cases = [
            (
                0,
                [
                    0xb25b_b249,
                    0xfaf6_4da7,
                    0x4ec0_2f2a,
                    0xad14_82f0,
                    0xe5f0_ada2,
                    0x3895_4814,
                    0xb5ea_7577,
                    0x4338_2910,
                ],
            ),
            (
                0x8000_0000,
                [
                    0xd166_6be0,
                    0x8dcf_331b,
                    0x6913_b0d9,
                    0x92ca_5746,
                    0x4fe3_fc18,
                    0x0653_7727,
                    0xb623_8692,
                    0x28bf_2625,
                ],
            ),
            (
                u32::MAX,
                [
                    0x1a48_3231,
                    0xb2f2_91cb,
                    0xbf95_818d,
                    0x0b21_7867,
                    0x6b8c_42b9,
                    0xe95e_5a85,
                    0xdcf0_0eb8,
                    0xbd51_a262,
                ],
            ),
        ];

        for (seed, expected) in cases {
            assert_eq!(first_eight(seed), expected, "seed {seed:#010x}");
        }
    }

    #[test]
    fn fingerprint_morgan_r01_default_and_reseed_paths() {
        assert_eq!(
            std::array::from_fn::<_, 8, _>({
                let mut engine = BoostFingerprintRng::default();
                move |_| engine.next_u32()
            }),
            [
                0x73b2_f002,
                0x462d_dd76,
                0x9514_a97b,
                0x393f_8d45,
                0xc38a_92ab,
                0xcc4d_2ecf,
                0x79fd_3d33,
                0x10a7_add8,
            ]
        );

        let mut reseeded = BoostFingerprintRng::new(0);
        let _ = reseeded.next_u32();
        reseeded.seed(42);
        assert_eq!(
            std::array::from_fn::<_, 4, _>(|_| reseeded.next_u32()),
            [0x704b_3dc5, 0x5215_b14b, 0x46c5_072b, 0xb44a_98be]
        );
    }

    #[test]
    fn fingerprint_morgan_r01_is_not_mt19937_or_core_minstd_rand() {
        let writer_words: [u32; 4] = {
            let mut engine = BoostFingerprintRng::new(42);
            std::array::from_fn(|_| engine.next_u32())
        };

        // Pinned Boost mt19937 uses 624 words; core::MinStdRand has its own
        // source-tested seed-42 sequence. Neither is the RDKit writer engine.
        assert_eq!(
            writer_words,
            [0x704b_3dc5, 0x5215_b14b, 0x46c5_072b, 0xb44a_98be]
        );
        assert_ne!(
            writer_words,
            [0x5fe1_dc66, 0xcbea_3db3, 0xf362_035c, 0x2ef5_950e]
        );
        assert_ne!(
            writer_words,
            [2_027_382, 1_226_992_407, 551_494_037, 961_371_815]
        );
    }

    #[test]
    fn fingerprint_morgan_r01_normalization_repairs_zero_state() {
        let mut engine = BoostFingerprintRng::new(42);
        let index = engine.index;
        engine.state = [0; STATE_WORDS];

        engine.normalize_state();

        assert_eq!(engine.state, [0x8000_0000, 0, 0, 0]);
        assert_eq!(engine.index, index);
    }

    #[test]
    fn fingerprint_morgan_r02_range_endpoints_and_zero_width() {
        let writer_words = [0, 1, 0xffff_fffe, u32::MAX];
        let mut writer_engine = LiteralEngine::new(0, u32::MAX, &writer_words);
        let writer_distribution = BoostUniformIntDistribution::rdkit_additional_bit();
        let writer_results =
            std::array::from_fn::<_, 4, _>(|_| writer_distribution.sample(&mut writer_engine));
        assert_eq!(writer_results, [0, 0, i32::MAX, i32::MAX]);
        assert_eq!(writer_engine.consumed(), 4);

        // The full signed i32 interval takes Boost's equal base/destination
        // range path; the literal endpoints exercise its safe signed add.
        let mut signed_engine = LiteralEngine::new(0, u32::MAX, &writer_words);
        let signed_distribution = BoostUniformIntDistribution::new(i32::MIN, i32::MAX);
        let signed_results =
            std::array::from_fn::<_, 4, _>(|_| signed_distribution.sample(&mut signed_engine));
        assert_eq!(
            signed_results,
            [i32::MIN, i32::MIN + 1, i32::MAX - 1, i32::MAX]
        );
        assert_eq!(signed_engine.consumed(), 4);

        // Boost returns the minimum before asking the engine for a word when
        // the requested range is zero.
        let mut unused_engine = LiteralEngine::new(0, 0, &[]);
        let singleton_distribution = BoostUniformIntDistribution::new(-7, -7);
        assert_eq!(singleton_distribution.sample(&mut unused_engine), -7);
        assert_eq!(unused_engine.consumed(), 0);
    }

    #[test]
    fn fingerprint_morgan_r02_rejection_and_unsigned_overflow() {
        // For range 0..=4 with a 0..=1 engine, the first two digits plus a
        // recursive increment produce 7 and are rejected; the next attempt
        // yields the fixed result 2.
        let rejection_words = [1, 1, 1, 0, 1, 0];
        let mut rejection_engine = LiteralEngine::new(0, 1, &rejection_words);
        let rejection_distribution = BoostUniformIntDistribution::new(0, 4);
        assert_eq!(rejection_distribution.sample(&mut rejection_engine), 2);
        assert_eq!(rejection_engine.consumed(), 6);

        // Twenty base-3 digits of 2 form 3^20-1. The recursive increment 1
        // adds 3^20, wrapping u32 to 2_678_601_505; Boost detects result <
        // increment, rejects that draw sequence, then accepts the zero row.
        let mut overflow_words = [0_u32; 42];
        overflow_words[..20].fill(2);
        overflow_words[20] = 1;
        assert_eq!(3_486_784_400_u32.wrapping_add(3_486_784_401), 2_678_601_505);
        assert!(2_678_601_505_u32 < 3_486_784_401_u32);

        let mut overflow_engine = LiteralEngine::new(0, 2, &overflow_words);
        let full_signed_distribution = BoostUniformIntDistribution::new(i32::MIN, i32::MAX);
        assert_eq!(
            full_signed_distribution.sample(&mut overflow_engine),
            i32::MIN
        );
        assert_eq!(overflow_engine.consumed(), 42);
    }

    #[test]
    fn fingerprint_morgan_r02_bits_per_feature_seeded_stream_matrix() {
        // One engine and distribution live across each source outer call.
        // Each environment reseeds that same engine, then consumes exactly
        // bitsPerFeature-1 draws. The repeated seed after two other seeds must
        // reproduce the same ordered suffix.
        let environment_seeds = [42, 0x8000_0000, u32::MAX, 42];
        let cases: [(u32, [&[i32]; 4]); 3] = [
            (1, [&[], &[], &[], &[]]),
            (
                2,
                [
                    &[0x3825_9ee2],
                    &[0x68b3_35f0],
                    &[0x0d24_1918],
                    &[0x3825_9ee2],
                ],
            ),
            (
                3,
                [
                    &[0x3825_9ee2, 0x290a_d8a5],
                    &[0x68b3_35f0, 0x46e7_998d],
                    &[0x0d24_1918, 0x5979_48e5],
                    &[0x3825_9ee2, 0x290a_d8a5],
                ],
            ),
        ];

        for (bits_per_feature, expected_rows) in cases {
            let mut actual_rows: [Vec<i32>; 4] = std::array::from_fn(|_| Vec::new());
            if bits_per_feature > 1 {
                let mut engine = BoostFingerprintRng::new(42);
                let distribution = BoostUniformIntDistribution::rdkit_additional_bit();
                for (row, seed) in environment_seeds.into_iter().enumerate() {
                    engine.seed(seed);
                    for _bit_number in 1..bits_per_feature {
                        actual_rows[row].push(distribution.sample(&mut engine));
                    }
                }
            }

            for row in 0..environment_seeds.len() {
                assert_eq!(
                    actual_rows[row].as_slice(),
                    expected_rows[row],
                    "bitsPerFeature={bits_per_feature}, environment seed={:#010x}",
                    environment_seeds[row]
                );
            }
            assert_eq!(actual_rows[0], actual_rows[3]);
        }
    }
}
