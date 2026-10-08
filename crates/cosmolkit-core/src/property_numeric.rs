//! Source arithmetic reads for canonical property values.
use cosmolkit_model::{PropertyText, PropertyValue, PropertyValueKind};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PropertyIntReadError {
    #[error("unsigned value {value} causes positive_overflow converting to signed int")]
    UnsignedOverflow { value: u32 },
    #[error("bad_any_cast reading {kind:?} as signed int")]
    InvalidKind { kind: PropertyValueKind },
    #[error(
        "bad_any_cast reading string {value:?}: decimal magnitude {magnitude} with negative={negative} exceeds signed int range"
    )]
    SignedTextOverflow {
        value: PropertyText,
        magnitude: u32,
        negative: bool,
    },
    #[error("bad_any_cast reading string {value:?} as signed int: {source}")]
    Lexical {
        value: PropertyText,
        #[source]
        source: UIntLexicalReadError,
    },
}

#[doc(hidden)]
pub fn property_value_to_int(value: &PropertyValue) -> Result<i32, PropertyIntReadError> {
    // RDKit✔️✔️: void getProp(const std::string_view key, T &res) const {
    // RDKit✔️✔️:     d_props.getVal(key, res);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️:   void getVal(const std::string_view what, T &res) const {
    // RDKit✔️✔️:     res = getVal<T>(what);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️:   T getVal(const std::string_view what) const {
    // RDKit✔️✔️:     for (auto &data : _data) {
    // RDKit✔️✔️:       if (data.key == what) {
    // RDKit✔️✔️:         return from_rdvalue<T>(data.val);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     throw KeyErrorException(what);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: template <class T>
    // RDKit✔️✔️: typename boost::enable_if<boost::is_arithmetic<T>, T>::type from_rdvalue(
    // RDKit✔️✔️:     RDValue_cast_t arg) {
    // RDKit✔️✔️:   T res;
    // RDKit✔️✔️:   if (arg.getTag() == RDTypeTag::StringTag) {
    // RDKit✔️✔️:     Utils::LocaleSwitcher ls;
    // RDKit✔️✔️:     try {
    // RDKit✔️✔️:       res = rdvalue_cast<T>(arg);
    // RDKit✔️✔️:     } catch (const std::bad_any_cast &exc) {
    // RDKit✔️✔️:       try {
    // RDKit✔️✔️: 	std::string val = rdvalue_cast<std::string>(arg);
    // RDKit✔️✔️: 	// trim only the right characters, this mimics how SD values
    // RDKit✔️✔️: 	//  work on read, they will be trimmed by the MolFile parser
    // RDKit✔️✔️: 	boost::trim_right(val);
    // RDKit✔️✔️:         res = boost::lexical_cast<T>(val);
    // RDKit✔️✔️:       } catch (...) {
    // RDKit✔️✔️:         throw exc;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res = rdvalue_cast<T>(arg);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: template <>
    // RDKit✔️✔️: inline int rdvalue_cast<int>(RDValue_cast_t v) {
    // RDKit✔️✔️:   if (rdvalue_is<int>(v)) {
    // RDKit✔️✔️:     return v.value.i;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (rdvalue_is<unsigned int>(v)) {
    // RDKit✔️✔️:     return boost::numeric_cast<int>(v.value.u);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   throw std::bad_any_cast();
    // RDKit✔️✔️: }
    // Boost✔️✔️: template <typename Type>
    // Boost✔️✔️:             bool shr_signed(Type& output) {
    // Boost✔️✔️:                 if (start == finish) return false;
    // Boost✔️✔️:                 CharT const minus = lcast_char_constants<CharT>::minus;
    // Boost✔️✔️:                 CharT const plus = lcast_char_constants<CharT>::plus;
    // Boost✔️✔️:                 typedef BOOST_DEDUCED_TYPENAME make_unsigned<Type>::type utype;
    // Boost✔️✔️:                 utype out_tmp = 0;
    // Boost✔️✔️:                 bool const has_minus = Traits::eq(minus, *start);
    // Boost✔️✔️:
    // Boost✔️✔️:                 /* We won`t use `start' any more, so no need in decrementing it after */
    // Boost✔️✔️:                 if (has_minus || Traits::eq(plus, *start)) {
    // Boost✔️✔️:                     ++start;
    // Boost✔️✔️:                 }
    // Boost✔️✔️:
    // Boost✔️✔️:                 bool succeed = lcast_ret_unsigned<Traits, utype, CharT>(out_tmp, start, finish).convert();
    // Boost✔️✔️:                 if (has_minus) {
    // Boost✔️✔️:                     utype const comp_val = (static_cast<utype>(1) << std::numeric_limits<Type>::digits);
    // Boost✔️✔️:                     succeed = succeed && out_tmp<=comp_val;
    // Boost✔️✔️:                     output = static_cast<Type>(0u - out_tmp);
    // Boost✔️✔️:                 } else {
    // Boost✔️✔️:                     utype const comp_val = static_cast<utype>((std::numeric_limits<Type>::max)());
    // Boost✔️✔️:                     succeed = succeed && out_tmp<=comp_val;
    // Boost✔️✔️:                     output = static_cast<Type>(out_tmp);
    // Boost✔️✔️:                 }
    // Boost✔️✔️:                 return succeed;
    // Boost✔️✔️:             }
    // This narrow value converter implements the reached from_rdvalue helper.
    // Canonical MODEL readers supply byte-key presence/required missing errors;
    // dictionary forwarding does not duplicate conversion or live runtime access.
    // C-locale string reads trim only the right source-classified space bytes,
    // then use the sign/magnitude helpers below; source rethrows bad_any_cast
    // on all lexical failures. Structural errors retain that category plus
    // original counted bytes and failure detail. Typed numeric-cast overflow
    // remains distinct. No arbitrary string decoding or failure-to-default.
    // Complexity: O(1) typed scalars; O(bytes) trim/decimal scan with constant
    // success-path storage. Borrowed trimming avoids the source string clone.
    match value {
        PropertyValue::Int(value) => Ok(*value),
        PropertyValue::UInt(value) => i32::try_from(*value)
            .map_err(|_| PropertyIntReadError::UnsignedOverflow { value: *value }),
        PropertyValue::String(value) => {
            // Boost✔️✔️:             template <typename Type>
            // Boost✔️✔️:             bool shr_signed(Type& output) {
            // Boost✔️✔️:                 if (start == finish) return false;
            // Boost✔️✔️:                 CharT const minus = lcast_char_constants<CharT>::minus;
            // Boost✔️✔️:                 CharT const plus = lcast_char_constants<CharT>::plus;
            // Boost✔️✔️:                 typedef BOOST_DEDUCED_TYPENAME make_unsigned<Type>::type utype;
            // Boost✔️✔️:                 utype out_tmp = 0;
            // Boost✔️✔️:                 bool const has_minus = Traits::eq(minus, *start);
            // Boost✔️✔️:
            // Boost✔️✔️:                 /* We won`t use `start' any more, so no need in decrementing it after */
            // Boost✔️✔️:                 if (has_minus || Traits::eq(plus, *start)) {
            // Boost✔️✔️:                     ++start;
            // Boost✔️✔️:                 }
            // Boost✔️✔️:
            // Boost✔️✔️:                 bool succeed = lcast_ret_unsigned<Traits, utype, CharT>(out_tmp, start, finish).convert();
            // Boost✔️✔️:                 if (has_minus) {
            // Boost✔️✔️:                     utype const comp_val = (static_cast<utype>(1) << std::numeric_limits<Type>::digits);
            // Boost✔️✔️:                     succeed = succeed && out_tmp<=comp_val;
            // Boost✔️✔️:                     output = static_cast<Type>(0u - out_tmp);
            // Boost✔️✔️:                 } else {
            // Boost✔️✔️:                     utype const comp_val = static_cast<utype>((std::numeric_limits<Type>::max)());
            // Boost✔️✔️:                     succeed = succeed && out_tmp<=comp_val;
            // Boost✔️✔️:                     output = static_cast<Type>(out_tmp);
            // Boost✔️✔️:                 }
            // Boost✔️✔️:                 return succeed;
            // Boost✔️✔️:             }
            // The source reads a sign and an unsigned decimal magnitude before
            // applying the signed bound. Borrow raw counted bytes throughout;
            // C-locale right trimming is unchanged, not a text decoding step.
            let bytes = trim_source_c_locale_right(value.as_bytes());
            let negative = bytes.first() == Some(&b'-');
            let start = usize::from(negative || bytes.first() == Some(&b'+'));
            let magnitude = parse_decimal_magnitude(&bytes[start..], start).map_err(|source| {
                PropertyIntReadError::Lexical {
                    value: value.clone(),
                    source,
                }
            })?;
            let maximum = if negative {
                1u32 << 31
            } else {
                i32::MAX as u32
            };
            if magnitude > maximum {
                return Err(PropertyIntReadError::SignedTextOverflow {
                    value: value.clone(),
                    magnitude,
                    negative,
                });
            }
            Ok(if negative {
                magnitude.wrapping_neg() as i32
            } else {
                magnitude as i32
            })
        }
        other => Err(PropertyIntReadError::InvalidKind { kind: other.kind() }),
    }
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum UIntLexicalReadError {
    #[error("empty unsigned decimal magnitude")]
    Empty,
    #[error("invalid decimal character {byte} at byte {position}")]
    Character { position: usize, byte: u8 },
    #[error("unsigned decimal magnitude exceeds u32")]
    Overflow,
}

fn trim_source_c_locale_right(bytes: &[u8]) -> &[u8] {
    // Boost✔️🔝: template<typename SequenceT>
    // Boost✔️🔝:         inline void trim_right(SequenceT& Input, const std::locale& Loc=std::locale())
    // Boost✔️🔝:         {
    // Boost✔️🔝:             ::boost::algorithm::trim_right_if(
    // Boost✔️🔝:                 Input,
    // Boost✔️🔝:                 is_space(Loc) );
    // Boost✔️🔝:         }
    // Boost✔️🔝: template<typename SequenceT, typename PredicateT>
    // Boost✔️🔝:         inline void trim_right_if(SequenceT& Input, PredicateT IsSpace)
    // Boost✔️🔝:         {
    // Boost✔️🔝:             Input.erase(
    // Boost✔️🔝:                 ::boost::algorithm::detail::trim_end(
    // Boost✔️🔝:                     ::boost::begin(Input),
    // Boost✔️🔝:                     ::boost::end(Input),
    // Boost✔️🔝:                     IsSpace ),
    // Boost✔️🔝:                 ::boost::end(Input)
    // Boost✔️🔝:                 );
    // Boost✔️🔝:         }
    // Boost✔️🔝: template< typename ForwardIteratorT, typename PredicateT >
    // Boost✔️🔝:             inline ForwardIteratorT trim_end(
    // Boost✔️🔝:                 ForwardIteratorT InBegin,
    // Boost✔️🔝:                 ForwardIteratorT InEnd,
    // Boost✔️🔝:                 PredicateT IsSpace )
    // Boost✔️🔝:             {
    // Boost✔️🔝:                 typedef BOOST_STRING_TYPENAME
    // Boost✔️🔝:                     std::iterator_traits<ForwardIteratorT>::iterator_category category;
    // Boost✔️🔝:
    // Boost✔️🔝:                 return ::boost::algorithm::detail::trim_end_iter_select( InBegin, InEnd, IsSpace, category() );
    // Boost✔️🔝:             }
    // Boost✔️🔝:             template< typename ForwardIteratorT, typename PredicateT >
    // Boost✔️🔝:             inline ForwardIteratorT trim_end_iter_select(
    // Boost✔️🔝:                 ForwardIteratorT InBegin,
    // Boost✔️🔝:                 ForwardIteratorT InEnd,
    // Boost✔️🔝:                 PredicateT IsSpace,
    // Boost✔️🔝:                 std::bidirectional_iterator_tag )
    // Boost✔️🔝:             {
    // Boost✔️🔝:                 for( ForwardIteratorT It=InEnd; It!=InBegin;  )
    // Boost✔️🔝:                 {
    // Boost✔️🔝:                     if ( !IsSpace(*(--It)) )
    // Boost✔️🔝:                         return ++It;
    // Boost✔️🔝:                 }
    // Boost✔️🔝:
    // Boost✔️🔝:                 return InBegin;
    // Boost✔️🔝:             }
    // Boost✔️🔝:             struct is_classifiedF :
    // Boost✔️🔝:                 public predicate_facade<is_classifiedF>
    // Boost✔️🔝:             {
    // Boost✔️🔝:                 // Boost.ResultOf support
    // Boost✔️🔝:                 typedef bool result_type;
    // Boost✔️🔝:
    // Boost✔️🔝:                 // Constructor from a locale
    // Boost✔️🔝:                 is_classifiedF(std::ctype_base::mask Type, std::locale const & Loc = std::locale()) :
    // Boost✔️🔝:                     m_Type(Type), m_Locale(Loc) {}
    // Boost✔️🔝:                 // Operation
    // Boost✔️🔝:                 template<typename CharT>
    // Boost✔️🔝:                 bool operator()( CharT Ch ) const
    // Boost✔️🔝:                 {
    // Boost✔️🔝:                     return std::use_facet< std::ctype<CharT> >(m_Locale).is( m_Type, Ch );
    // Boost✔️🔝:                 }
    // Boost✔️🔝:
    // Boost✔️🔝:                 #if defined(BOOST_BORLANDC) && (BOOST_BORLANDC >= 0x560) && (BOOST_BORLANDC <= 0x582) && !defined(_USE_OLD_RW_STL)
    // Boost✔️🔝:                     template<>
    // Boost✔️🔝:                     bool operator()( char const Ch ) const
    // Boost✔️🔝:                     {
    // Boost✔️🔝:                         return std::use_facet< std::ctype<char> >(m_Locale).is( m_Type, Ch );
    // Boost✔️🔝:                     }
    // Boost✔️🔝:                 #endif
    // Boost✔️🔝:
    // Boost✔️🔝:             private:
    // Boost✔️🔝:                 std::ctype_base::mask m_Type;
    // Boost✔️🔝:                 std::locale m_Locale;
    // Boost✔️🔝:             };
    // Boost✔️🔝:         inline detail::is_classifiedF
    // Boost✔️🔝:         is_space(const std::locale& Loc=std::locale())
    // Boost✔️🔝:         {
    // Boost✔️🔝:             return detail::is_classifiedF(std::ctype_base::space, Loc);
    // Boost✔️🔝:         }
    // The primary C-locale ctype table classifies exactly 0x09..=0x0d and
    // 0x20 as space. In particular vertical tab 0x0b is space; Rust's
    // trim_ascii_end omits it and therefore was not a source-equivalent trim.
    // Boost std::string uses the bidirectional helper. Keep its reverse-only
    // trim boundary with a borrowed slice instead of cloning/erasing a string.
    // Source witness: Step0136-boost181-trim-source/{full-read,c-locale-source}.json.
    // Behavior: no leading/interior trimming, NUL or high-byte classification.
    // Complexity: O(trailing bytes), O(1) state, no success-path allocation.
    let end = bytes
        .iter()
        .rposition(|byte| !matches!(*byte, b'\t'..=b'\r' | b' '))
        .map_or(0, |position| position + 1);
    &bytes[..end]
}

fn parse_decimal_magnitude(bytes: &[u8], offset: usize) -> Result<u32, UIntLexicalReadError> {
    parse_decimal_magnitude_with_limit(bytes, offset, u64::from(u32::MAX)).map(|value| value as u32)
}

fn parse_decimal_magnitude_with_limit(
    bytes: &[u8],
    offset: usize,
    maximum: u64,
) -> Result<u64, UIntLexicalReadError> {
    // BEGIN BOOST CPP INPUT CONVERTER lcast_ret_unsigned
    // Boost✔️✔️:         template <class Traits, class T, class CharT>
    // Boost✔️✔️:         class lcast_ret_unsigned: boost::noncopyable {
    // Boost✔️✔️:             bool m_multiplier_overflowed;
    // Boost✔️✔️:             T m_multiplier;
    // Boost✔️✔️:             T& m_value;
    // Boost✔️✔️:             const CharT* const m_begin;
    // Boost✔️✔️:             const CharT* m_end;
    // Boost✔️✔️:
    // Boost✔️✔️:         public:
    // Boost✔️✔️:             lcast_ret_unsigned(T& value, const CharT* const begin, const CharT* end) BOOST_NOEXCEPT
    // Boost✔️✔️:                 : m_multiplier_overflowed(false), m_multiplier(1), m_value(value), m_begin(begin), m_end(end)
    // Boost✔️✔️:             {
    // Boost✔️✔️: #ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
    // Boost✔️✔️:                 BOOST_STATIC_ASSERT(!std::numeric_limits<T>::is_signed);
    // Boost✔️✔️:
    // Boost✔️✔️:                 // GCC when used with flag -std=c++0x may not have std::numeric_limits
    // Boost✔️✔️:                 // specializations for __int128 and unsigned __int128 types.
    // Boost✔️✔️:                 // Try compilation with -std=gnu++0x or -std=gnu++11.
    // Boost✔️✔️:                 //
    // Boost✔️✔️:                 // http://gcc.gnu.org/bugzilla/show_bug.cgi?id=40856
    // Boost✔️✔️:                 BOOST_STATIC_ASSERT_MSG(std::numeric_limits<T>::is_specialized,
    // Boost✔️✔️:                     "std::numeric_limits are not specialized for integral type passed to boost::lexical_cast"
    // Boost✔️✔️:                 );
    // Boost✔️✔️: #endif
    // Boost✔️✔️:             }
    // Boost✔️✔️:
    // Boost✔️✔️:             inline bool convert() {
    // Boost✔️✔️:                 CharT const czero = lcast_char_constants<CharT>::zero;
    // Boost✔️✔️:                 --m_end;
    // Boost✔️✔️:                 m_value = static_cast<T>(0);
    // Boost✔️✔️:
    // Boost✔️✔️:                 if (m_begin > m_end || *m_end < czero || *m_end >= czero + 10)
    // Boost✔️✔️:                     return false;
    // Boost✔️✔️:                 m_value = static_cast<T>(*m_end - czero);
    // Boost✔️✔️:                 --m_end;
    // Boost✔️✔️:
    // Boost✔️✔️: #ifdef BOOST_LEXICAL_CAST_ASSUME_C_LOCALE
    // Boost✔️✔️:                 return main_convert_loop();
    // Boost✔️✔️: #else
    // Boost✔️✔️:                 std::locale loc;
    // Boost✔️✔️:                 if (loc == std::locale::classic()) {
    // Boost✔️✔️:                     return main_convert_loop();
    // Boost✔️✔️:                 }
    // Boost✔️✔️:
    // Boost❌❌:                 typedef std::numpunct<CharT> numpunct;
    // Boost❌❌:                 numpunct const& np = BOOST_USE_FACET(numpunct, loc);
    // Boost❌❌:                 std::string const& grouping = np.grouping();
    // Boost❌❌:                 std::string::size_type const grouping_size = grouping.size();
    // Boost❌❌:
    // Boost❌❌:                 /* According to Programming languages - C++
    // Boost❌❌:                  * we MUST check for correct grouping
    // Boost❌❌:                  */
    // Boost❌❌:                 if (!grouping_size || grouping[0] <= 0) {
    // Boost❌❌:                     return main_convert_loop();
    // Boost❌❌:                 }
    // Boost❌❌:
    // Boost❌❌:                 unsigned char current_grouping = 0;
    // Boost❌❌:                 CharT const thousands_sep = np.thousands_sep();
    // Boost❌❌:                 char remained = static_cast<char>(grouping[current_grouping] - 1);
    // Boost❌❌:
    // Boost❌❌:                 for (;m_end >= m_begin; --m_end)
    // Boost❌❌:                 {
    // Boost❌❌:                     if (remained) {
    // Boost❌❌:                         if (!main_convert_iteration()) {
    // Boost❌❌:                             return false;
    // Boost❌❌:                         }
    // Boost❌❌:                         --remained;
    // Boost❌❌:                     } else {
    // Boost❌❌:                         if ( !Traits::eq(*m_end, thousands_sep) ) //|| begin == end ) return false;
    // Boost❌❌:                         {
    // Boost❌❌:                             /*
    // Boost❌❌:                              * According to Programming languages - C++
    // Boost❌❌:                              * Digit grouping is checked. That is, the positions of discarded
    // Boost❌❌:                              * separators is examined for consistency with
    // Boost❌❌:                              * use_facet<numpunct<charT> >(loc ).grouping()
    // Boost❌❌:                              *
    // Boost❌❌:                              * BUT what if there is no separators at all and grouping()
    // Boost❌❌:                              * is not empty? Well, we have no extraced separators, so we
    // Boost❌❌:                              * won`t check them for consistency. This will allow us to
    // Boost❌❌:                              * work with "C" locale from other locales
    // Boost❌❌:                              */
    // Boost❌❌:                             return main_convert_loop();
    // Boost❌❌:                         } else {
    // Boost❌❌:                             if (m_begin == m_end) return false;
    // Boost❌❌:                             if (current_grouping < grouping_size - 1) ++current_grouping;
    // Boost❌❌:                             remained = grouping[current_grouping];
    // Boost❌❌:                         }
    // Boost❌❌:                     }
    // Boost❌❌:                 } /*for*/
    // Boost❌❌:
    // Boost❌❌:                 return true;
    // Boost✔️✔️: #endif
    // Boost✔️✔️:             }
    // Boost✔️✔️:
    // Boost✔️✔️:         private:
    // Boost✔️✔️:             // Iteration that does not care about grouping/separators and assumes that all
    // Boost✔️✔️:             // input characters are digits
    // Boost✔️✔️:             inline bool main_convert_iteration() BOOST_NOEXCEPT {
    // Boost✔️✔️:                 CharT const czero = lcast_char_constants<CharT>::zero;
    // Boost✔️✔️:                 T const maxv = (std::numeric_limits<T>::max)();
    // Boost✔️✔️:
    // Boost✔️✔️:                 m_multiplier_overflowed = m_multiplier_overflowed || (maxv/10 < m_multiplier);
    // Boost✔️✔️:                 m_multiplier = static_cast<T>(m_multiplier * 10);
    // Boost✔️✔️:
    // Boost✔️✔️:                 T const dig_value = static_cast<T>(*m_end - czero);
    // Boost✔️✔️:                 T const new_sub_value = static_cast<T>(m_multiplier * dig_value);
    // Boost✔️✔️:
    // Boost✔️✔️:                 // We must correctly handle situations like `000000000000000000000000000001`.
    // Boost✔️✔️:                 // So we take care of overflow only if `dig_value` is not '0'.
    // Boost✔️✔️:                 if (*m_end < czero || *m_end >= czero + 10  // checking for correct digit
    // Boost✔️✔️:                     || (dig_value && (                      // checking for overflow of ...
    // Boost✔️✔️:                         m_multiplier_overflowed                             // ... multiplier
    // Boost✔️✔️:                         || static_cast<T>(maxv / dig_value) < m_multiplier  // ... subvalue
    // Boost✔️✔️:                         || static_cast<T>(maxv - new_sub_value) < m_value   // ... whole expression
    // Boost✔️✔️:                     ))
    // Boost✔️✔️:                 ) return false;
    // Boost✔️✔️:
    // Boost✔️✔️:                 m_value = static_cast<T>(m_value + new_sub_value);
    // Boost✔️✔️:
    // Boost✔️✔️:                 return true;
    // Boost✔️✔️:             }
    // Boost✔️✔️:
    // Boost✔️✔️:             bool main_convert_loop() BOOST_NOEXCEPT {
    // Boost✔️✔️:                 for ( ; m_end >= m_begin; --m_end) {
    // Boost✔️✔️:                     if (!main_convert_iteration()) {
    // Boost✔️✔️:                         return false;
    // Boost✔️✔️:                     }
    // Boost✔️✔️:                 }
    // Boost✔️✔️:
    // Boost✔️✔️:                 return true;
    // Boost✔️✔️:             }
    // Boost✔️✔️:         };
    // END BOOST CPP INPUT CONVERTER lcast_ret_unsigned
    // Behavior: the modeled C-locale branch accepts a nonempty decimal run
    // iff its magnitude fits its actual unsigned source width. Arbitrarily many leading zeros are accepted.
    // Source reverse weighting and this forward checked accumulation have the
    // same success/value domain; failure returns the enclosing bad_any_cast,
    // with additional structural diagnostic detail retained locally.
    // Nonclassic locale grouping is independent and not supplied by this
    // source C-locale chemical path. No Unicode digit/whitespace acceptance.
    // Complexity: one digit pass, O(1) state, no allocation or text decoding.
    if bytes.is_empty() {
        return Err(UIntLexicalReadError::Empty);
    }
    let mut magnitude = 0u64;
    for (position, &byte) in bytes.iter().enumerate() {
        if !byte.is_ascii_digit() {
            return Err(UIntLexicalReadError::Character {
                position: position + offset,
                byte,
            });
        }
        magnitude = magnitude
            .checked_mul(10)
            .and_then(|value| value.checked_add(u64::from(byte - b'0')))
            .filter(|value| *value <= maximum)
            .ok_or(UIntLexicalReadError::Overflow)?;
    }
    Ok(magnitude)
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PropertyUIntReadError {
    #[error("signed value {value} causes negative_overflow converting to unsigned int")]
    Negative { value: i32 },
    #[error("bad_any_cast reading {kind:?} as unsigned int")]
    InvalidKind { kind: PropertyValueKind },
    #[error("bad_any_cast reading string {value:?} as unsigned int: {source}")]
    Lexical {
        value: PropertyText,
        #[source]
        source: UIntLexicalReadError,
    },
}

#[doc(hidden)]
pub fn property_value_to_uint(value: &PropertyValue) -> Result<u32, PropertyUIntReadError> {
    // RDKit✔️✔️: void getProp(const std::string_view key, T &res) const {
    // RDKit✔️✔️:     d_props.getVal(key, res);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️:   void getVal(const std::string_view what, T &res) const {
    // RDKit✔️✔️:     res = getVal<T>(what);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️:   T getVal(const std::string_view what) const {
    // RDKit✔️✔️:     for (auto &data : _data) {
    // RDKit✔️✔️:       if (data.key == what) {
    // RDKit✔️✔️:         return from_rdvalue<T>(data.val);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     throw KeyErrorException(what);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: template <class T>
    // RDKit✔️✔️: typename boost::enable_if<boost::is_arithmetic<T>, T>::type from_rdvalue(
    // RDKit✔️✔️:     RDValue_cast_t arg) {
    // RDKit✔️✔️:   T res;
    // RDKit✔️✔️:   if (arg.getTag() == RDTypeTag::StringTag) {
    // RDKit✔️✔️:     Utils::LocaleSwitcher ls;
    // RDKit✔️✔️:     try {
    // RDKit✔️✔️:       res = rdvalue_cast<T>(arg);
    // RDKit✔️✔️:     } catch (const std::bad_any_cast &exc) {
    // RDKit✔️✔️:       try {
    // RDKit✔️✔️: 	std::string val = rdvalue_cast<std::string>(arg);
    // RDKit✔️✔️: 	// trim only the right characters, this mimics how SD values
    // RDKit✔️✔️: 	//  work on read, they will be trimmed by the MolFile parser
    // RDKit✔️✔️: 	boost::trim_right(val);
    // RDKit✔️✔️:         res = boost::lexical_cast<T>(val);
    // RDKit✔️✔️:       } catch (...) {
    // RDKit✔️✔️:         throw exc;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res = rdvalue_cast<T>(arg);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: template <>
    // RDKit✔️✔️: inline unsigned int rdvalue_cast<unsigned int>(RDValue_cast_t v) {
    // RDKit✔️✔️:   if (rdvalue_is<unsigned int>(v)) {
    // RDKit✔️✔️:     return v.value.u;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (rdvalue_is<int>(v)) {
    // RDKit✔️✔️:     return boost::numeric_cast<unsigned int>(v.value.i);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   throw std::bad_any_cast();
    // RDKit✔️✔️: }
    // Boost✔️✔️: template <typename Type>
    // Boost✔️✔️:             bool shr_unsigned(Type& output) {
    // Boost✔️✔️:                 if (start == finish) return false;
    // Boost✔️✔️:                 CharT const minus = lcast_char_constants<CharT>::minus;
    // Boost✔️✔️:                 CharT const plus = lcast_char_constants<CharT>::plus;
    // Boost✔️✔️:                 bool const has_minus = Traits::eq(minus, *start);
    // Boost✔️✔️:
    // Boost✔️✔️:                 /* We won`t use `start' any more, so no need in decrementing it after */
    // Boost✔️✔️:                 if (has_minus || Traits::eq(plus, *start)) {
    // Boost✔️✔️:                     ++start;
    // Boost✔️✔️:                 }
    // Boost✔️✔️:
    // Boost✔️✔️:                 bool const succeed = lcast_ret_unsigned<Traits, Type, CharT>(output, start, finish).convert();
    // Boost✔️✔️:
    // Boost✔️✔️:                 if (has_minus) {
    // Boost✔️✔️:                     output = static_cast<Type>(0u - output);
    // Boost✔️✔️:                 }
    // Boost✔️✔️:
    // Boost✔️✔️:                 return succeed;
    // Boost✔️✔️:             }
    // This narrow value converter implements the reached from_rdvalue helper.
    // Canonical MODEL readers supply byte-key presence/required missing errors;
    // dictionary forwarding does not duplicate conversion or live runtime access.
    // C-locale string reads trim only the right source-classified space bytes,
    // then use the sign/magnitude helpers below; source rethrows bad_any_cast
    // on all lexical failures. Structural errors retain that category plus
    // original counted bytes and failure detail. Typed numeric-cast overflow
    // remains distinct. No arbitrary string decoding or failure-to-default.
    // Complexity: O(1) typed scalars; O(bytes) trim/decimal scan with constant
    // success-path storage. Borrowed trimming avoids the source string clone.
    match value {
        PropertyValue::UInt(value) => Ok(*value),
        PropertyValue::Int(value) => {
            u32::try_from(*value).map_err(|_| PropertyUIntReadError::Negative { value: *value })
        }
        PropertyValue::String(value) => {
            let bytes = trim_source_c_locale_right(value.as_bytes());
            let negative = bytes.first() == Some(&b'-');
            let start = usize::from(negative || bytes.first() == Some(&b'+'));
            let parsed = parse_decimal_magnitude(&bytes[start..], start).map(|magnitude| {
                if negative {
                    magnitude.wrapping_neg()
                } else {
                    magnitude
                }
            });
            parsed.map_err(|source| PropertyUIntReadError::Lexical {
                value: value.clone(),
                source,
            })
        }
        other => Err(PropertyUIntReadError::InvalidKind { kind: other.kind() }),
    }
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PropertyDoubleReadError {
    #[error("bad_any_cast reading {kind:?} as double")]
    InvalidKind { kind: PropertyValueKind },
    #[error("bad_any_cast reading double property {value:?}: {source}")]
    Lexical {
        value: PropertyText,
        #[source]
        source: DoubleLexicalReadError,
    },
}

#[derive(Debug, thiserror::Error)]
pub enum UnsignedStreamArrayError {
    #[error("unsigned stream array first value is unavailable")]
    MissingFirstValue,
    #[error("unsigned stream array allocation failed: {0}")]
    Allocation(#[from] std::collections::TryReserveError),
}

use std::cmp::Ordering;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SourceNumericRoundingMode {
    NearestEven,
    Downward,
    Upward,
    TowardZero,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct SourceNumericStreamState {
    pub eof: bool,
    pub fail: bool,
    pub bad: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum DoubleLexicalReadErrorKind {
    Syntax,
    TrailingByte,
    Overflow,
    FloatingEnvironmentNotModeled,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DoubleLexicalReadError {
    pub kind: DoubleLexicalReadErrorKind,
    pub consumed: usize,
    pub state: SourceNumericStreamState,
    pub value_bits: Option<u64>,
}

impl std::fmt::Display for DoubleLexicalReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "bad lexical cast: source type value could not be interpreted as target ({:?}, byte {})",
            self.kind, self.consumed
        )
    }
}
impl std::error::Error for DoubleLexicalReadError {}

#[derive(Clone, Debug, PartialEq)]
pub struct SourceDoubleExtraction {
    pub value: f64,
    pub consumed: usize,
    pub state: SourceNumericStreamState,
    pub failure: Option<DoubleLexicalReadErrorKind>,
}

// Values are exact unsigned integers in radix 2^32. The largest comparison
// has fewer than 4,700 bits: 10^1092 times a binary64 coefficient/exponent.
// Every exponent here is bounded after full lexical/exponent cancellation.
#[derive(Clone, Debug)]
struct Binary64Integer(Vec<u32>);

impl Binary64Integer {
    fn from_u64(value: u64) -> Self {
        if value >> 32 != 0 {
            Self(vec![value as u32, (value >> 32) as u32])
        } else {
            Self(vec![value as u32])
        }
    }
    fn multiply_add(&mut self, factor: u64, add: u64) {
        let mut carry = u128::from(add);
        for limb in &mut self.0 {
            let value = u128::from(*limb) * u128::from(factor) + carry;
            *limb = value as u32;
            carry = value >> 32;
        }
        while carry != 0 {
            self.0.push(carry as u32);
            carry >>= 32;
        }
        while self.0.len() > 1 && self.0.last() == Some(&0) {
            self.0.pop();
        }
    }
    fn shift_left(&mut self, shift: usize) {
        if self.0.len() == 1 && self.0[0] == 0 {
            return;
        }
        let words = shift / 32;
        let remainder = shift % 32;
        if remainder != 0 {
            let mut carry = 0_u64;
            for limb in &mut self.0 {
                let value = (u64::from(*limb) << remainder) | carry;
                *limb = value as u32;
                carry = value >> 32;
            }
            if carry != 0 {
                self.0.push(carry as u32);
            }
        }
        if words != 0 {
            let old = self.0.len();
            self.0.resize(old + words, 0);
            self.0.copy_within(0..old, words);
            self.0[..words].fill(0);
        }
    }
    fn compare(&self, other: &Self) -> Ordering {
        self.0
            .len()
            .cmp(&other.0.len())
            .then_with(|| self.0.iter().rev().cmp(other.0.iter().rev()))
    }
}

struct ExactDecimal {
    numerator: Binary64Integer,
    denominator: Binary64Integer,
}
impl ExactDecimal {
    fn new(digits: &[u8], normalized_exponent: i128) -> Self {
        let mut numerator = Binary64Integer::from_u64(0);
        for &digit in digits {
            numerator.multiply_add(10, u64::from(digit - b'0'));
        }
        let mut denominator = Binary64Integer::from_u64(1);
        let power = normalized_exponent - digits.len() as i128 + 1;
        if power >= 0 {
            for _ in 0..power {
                numerator.multiply_add(10, 0);
            }
        } else {
            for _ in 0..-power {
                denominator.multiply_add(10, 0);
            }
        }
        Self {
            numerator,
            denominator,
        }
    }
    fn compare_binary(&self, coefficient: u64, exponent: i32) -> Ordering {
        let mut left = self.numerator.clone();
        let mut right = self.denominator.clone();
        right.multiply_add(coefficient, 0);
        if exponent >= 0 {
            right.shift_left(exponent as usize);
        } else {
            left.shift_left((-exponent) as usize);
        }
        left.compare(&right)
    }
}

const MAX_FINITE_BITS: u64 = 0x7fef_ffff_ffff_ffff;
const INFINITY_BITS: u64 = 0x7ff0_0000_0000_0000;
const SIGN_BIT: u64 = 1 << 63;

fn binary64_parts(bits: u64) -> (u64, i32) {
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & ((1_u64 << 52) - 1);
    if exponent == 0 {
        (fraction, -1074)
    } else {
        ((1_u64 << 52) | fraction, exponent - 1023 - 52)
    }
}

#[rustfmt::skip]
fn round_away(
    negative: bool,
    last_digit_odd: bool,
    half_bit: bool,
    more_bits: bool,
    mode: SourceNumericRoundingMode,
) -> bool {
    // /* Handle floating-point rounding mode within libc.
    //    Copyright (C) 2012-2026 Free Software Foundation, Inc.
    //    This file is part of the GNU C Library.
    // 
    //    The GNU C Library is free software; you can redistribute it and/or
    //    modify it under the terms of the GNU Lesser General Public
    //    License as published by the Free Software Foundation; either
    //    version 2.1 of the License, or (at your option) any later version.
    // 
    //    The GNU C Library is distributed in the hope that it will be useful,
    //    but WITHOUT ANY WARRANTY; without even the implied warranty of
    //    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    //    Lesser General Public License for more details.
    // 
    //    You should have received a copy of the GNU Lesser General Public
    //    License along with the GNU C Library; if not, see
    //    <https://www.gnu.org/licenses/>.  */
    // 
    // #ifndef _ROUNDING_MODE_H
    // #define _ROUNDING_MODE_H	1
    // 
    // #include <fenv.h>
    // #include <stdbool.h>
    // #include <stdlib.h>
    // 
    // /* Get the architecture-specific definition of how to determine the
    //    rounding mode in libc.  This header must also define the FE_*
    //    macros for any standard rounding modes the architecture does not
    //    have in <fenv.h>, to arbitrary distinct values.  */
    // #include <get-rounding-mode.h>
    // 
    // /* Return true if a number should be rounded away from zero in
    //    rounding mode MODE, false otherwise.  NEGATIVE is true if the
    //    number is negative, false otherwise.  LAST_DIGIT_ODD is true if the
    //    last digit of the truncated value (last bit for binary) is odd,
    //    false otherwise.  HALF_BIT is true if the number is at least half
    //    way from the truncated value to the next value with the
    //    least-significant digit in the same place, false otherwise.
    //    MORE_BITS is true if the number is not exactly equal to the
    //    truncated value or the half-way value, false otherwise.  */
    // 
    // static bool
    // round_away (bool negative, bool last_digit_odd, bool half_bit, bool more_bits,
    // 	    int mode)
    // {
    //   switch (mode)
    //     {
    //     case FE_DOWNWARD:
    //       return negative && (half_bit || more_bits);
    // 
    //     case FE_TONEAREST:
    //       return half_bit && (last_digit_odd || more_bits);
    // 
    //     case FE_TOWARDZERO:
    //       return false;
    // 
    //     case FE_UPWARD:
    //       return !negative && (half_bit || more_bits);
    // 
    //     default:
    //       abort ();
    //     }
    // }
    // 
    // #endif /* rounding-mode.h */
    // glibc❗✔️: all four source branches are selected by explicit state.
    match mode {
        SourceNumericRoundingMode::Downward => negative && (half_bit || more_bits),
        SourceNumericRoundingMode::NearestEven => half_bit && (last_digit_odd || more_bits),
        SourceNumericRoundingMode::TowardZero => false,
        SourceNumericRoundingMode::Upward => !negative && (half_bit || more_bits),
    }
}

#[rustfmt::skip]
fn decimal_binary64_bits(
    digits: &[u8],
    exponent: i128,
    negative: bool,
    mode: SourceNumericRoundingMode,
) -> u64 {
    // /* Convert string representing a number to float value, using given locale.
    //    Copyright (C) 1997-2026 Free Software Foundation, Inc.
    //    This file is part of the GNU C Library.
    // 
    //    The GNU C Library is free software; you can redistribute it and/or
    //    modify it under the terms of the GNU Lesser General Public
    //    License as published by the Free Software Foundation; either
    //    version 2.1 of the License, or (at your option) any later version.
    // 
    //    The GNU C Library is distributed in the hope that it will be useful,
    //    but WITHOUT ANY WARRANTY; without even the implied warranty of
    //    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    //    Lesser General Public License for more details.
    // 
    //    You should have received a copy of the GNU Lesser General Public
    //    License along with the GNU C Library; if not, see
    //    <https://www.gnu.org/licenses/>.  */
    // 
    // #include <bits/floatn.h>
    // 
    // #ifdef FLOAT
    // # define BUILD_DOUBLE 0
    // #else
    // # define BUILD_DOUBLE 1
    // #endif
    // 
    // #if BUILD_DOUBLE
    // # if __HAVE_FLOAT64 && !__HAVE_DISTINCT_FLOAT64
    // #  define strtof64_l __hide_strtof64_l
    // #  define wcstof64_l __hide_wcstof64_l
    // # endif
    // # if __HAVE_FLOAT32X && !__HAVE_DISTINCT_FLOAT32X
    // #  define strtof32x_l __hide_strtof32x_l
    // #  define wcstof32x_l __hide_wcstof32x_l
    // # endif
    // #endif
    // 
    // #include <locale.h>
    // 
    // extern double ____strtod_l_internal (const char *, char **, int, locale_t);
    // 
    // /* Configuration part.  These macros are defined by `strtold.c',
    //    `strtof.c', `wcstod.c', `wcstold.c', and `wcstof.c' to produce the
    //    `long double' and `float' versions of the reader.  */
    // #ifndef FLOAT
    // # include <math_ldbl_opt.h>
    // # define FLOAT		double
    // # define FLT		DBL
    // # ifdef USE_WIDE_CHAR
    // #  define STRTOF	wcstod_l
    // #  define __STRTOF	__wcstod_l
    // #  define STRTOF_NAN	__wcstod_nan
    // # else
    // #  define STRTOF	strtod_l
    // #  define __STRTOF	__strtod_l
    // #  define STRTOF_NAN	__strtod_nan
    // # endif
    // # define MPN2FLOAT	__mpn_construct_double
    // # define FLOAT_HUGE_VAL	HUGE_VAL
    // #endif
    // /* End of configuration part.  */
    // 
    // 
    // #include <ctype.h>
    // #include <errno.h>
    // #include <float.h>
    // #include "../locale/localeinfo.h"
    // #include <math.h>
    // #include <math-barriers.h>
    // #include <math-narrow-eval.h>
    // #include <stdlib.h>
    // #include <string.h>
    // #include <stdint.h>
    // #include <rounding-mode.h>
    // #include <tininess.h>
    // #include <stdbit.h>
    // 
    // /* The gmp headers need some configuration frobs.  */
    // #define HAVE_ALLOCA 1
    // 
    // /* Include gmp-mparam.h first, such that definitions of _SHORT_LIMB
    //    and _LONG_LONG_LIMB in it can take effect into gmp.h.  */
    // #include <gmp-mparam.h>
    // #include <gmp.h>
    // #include "gmp-impl.h"
    // #include "fpioconst.h"
    // 
    // #include <assert.h>
    // 
    // 
    // /* We use this code for the extended locale handling where the
    //    function gets as an additional argument the locale which has to be
    //    used.  To access the values we have to redefine the _NL_CURRENT and
    //    _NL_CURRENT_WORD macros.  */
    // #undef _NL_CURRENT
    // #define _NL_CURRENT(category, item) \
    //   (current->values[_NL_ITEM_INDEX (item)].string)
    // #undef _NL_CURRENT_WORD
    // #define _NL_CURRENT_WORD(category, item) \
    //   ((uint32_t) current->values[_NL_ITEM_INDEX (item)].word)
    // 
    // #if defined _LIBC || defined HAVE_WCHAR_H
    // # include <wchar.h>
    // #endif
    // 
    // #ifdef USE_WIDE_CHAR
    // # include <wctype.h>
    // # define STRING_TYPE wchar_t
    // # define CHAR_TYPE wint_t
    // # define L_(Ch) L##Ch
    // # define ISSPACE(Ch) __iswspace_l ((Ch), loc)
    // # define ISDIGIT(Ch) __iswdigit_l ((Ch), loc)
    // # define ISXDIGIT(Ch) __iswxdigit_l ((Ch), loc)
    // # define TOLOWER(Ch) __towlower_l ((Ch), loc)
    // # define TOLOWER_C(Ch) __towlower_l ((Ch), _nl_C_locobj_ptr)
    // # define STRNCASECMP(S1, S2, N) \
    //   __wcsncasecmp_l ((S1), (S2), (N), _nl_C_locobj_ptr)
    // #else
    // # define STRING_TYPE char
    // # define CHAR_TYPE char
    // # define L_(Ch) Ch
    // # define ISSPACE(Ch) __isspace_l ((Ch), loc)
    // # define ISDIGIT(Ch) __isdigit_l ((Ch), loc)
    // # define ISXDIGIT(Ch) __isxdigit_l ((Ch), loc)
    // # define TOLOWER(Ch) __tolower_l ((Ch), loc)
    // # define TOLOWER_C(Ch) __tolower_l ((Ch), _nl_C_locobj_ptr)
    // # define STRNCASECMP(S1, S2, N) \
    //   __strncasecmp_l ((S1), (S2), (N), _nl_C_locobj_ptr)
    // #endif
    // 
    // 
    // /* Constants we need from float.h; select the set for the FLOAT precision.  */
    // #define MANT_DIG	PASTE(FLT,_MANT_DIG)
    // #define	DIG		PASTE(FLT,_DIG)
    // #define	MAX_EXP		PASTE(FLT,_MAX_EXP)
    // #define	MIN_EXP		PASTE(FLT,_MIN_EXP)
    // #define MAX_10_EXP	PASTE(FLT,_MAX_10_EXP)
    // #define MIN_10_EXP	PASTE(FLT,_MIN_10_EXP)
    // #define MAX_VALUE	PASTE(FLT,_MAX)
    // #define MIN_VALUE	PASTE(FLT,_MIN)
    // 
    // /* Extra macros required to get FLT expanded before the pasting.  */
    // #define PASTE(a,b)	PASTE1(a,b)
    // #define PASTE1(a,b)	a##b
    // 
    // /* Function to construct a floating point number from an MP integer
    //    containing the fraction bits, a base 2 exponent, and a sign flag.  */
    // extern FLOAT MPN2FLOAT (mp_srcptr mpn, int exponent, int negative);
    // 
    // 
    // /* Definitions according to limb size used.  */
    // #if	BITS_PER_MP_LIMB == 32
    // # define MAX_DIG_PER_LIMB	9
    // # define MAX_FAC_PER_LIMB	1000000000UL
    // #elif	BITS_PER_MP_LIMB == 64
    // # define MAX_DIG_PER_LIMB	19
    // # define MAX_FAC_PER_LIMB	10000000000000000000ULL
    // #else
    // # error "mp_limb_t size " BITS_PER_MP_LIMB "not accounted for"
    // #endif
    // 
    // extern const mp_limb_t _tens_in_limb[MAX_DIG_PER_LIMB + 1];
    // 
    // 
    // #ifndef	howmany
    // #define	howmany(x,y)		(((x)+((y)-1))/(y))
    // #endif
    // #define SWAP(x, y)		({ typeof(x) _tmp = x; x = y; y = _tmp; })
    // 
    // #define	RETURN_LIMB_SIZE		howmany (MANT_DIG, BITS_PER_MP_LIMB)
    // 
    // #define RETURN(val,end)							      \
    //     do { if (endptr != NULL) *endptr = (STRING_TYPE *) (end);		      \
    // 	 return val; } while (0)
    // 
    // /* Maximum size necessary for mpn integers to hold floating point
    //    numbers.  The largest number we need to hold is 10^n where 2^-n is
    //    1/4 ulp of the smallest representable value (that is, n = MANT_DIG
    //    - MIN_EXP + 2).  Approximate using 10^3 < 2^10.  */
    // #define	MPNSIZE		(howmany (1 + ((MANT_DIG - MIN_EXP + 2) * 10) / 3, \
    // 				  BITS_PER_MP_LIMB) + 2)
    // /* Declare an mpn integer variable that big.  */
    // #define	MPN_VAR(name)	mp_limb_t name[MPNSIZE]; mp_size_t name##size
    // /* Copy an mpn integer value.  */
    // #define MPN_ASSIGN(dst, src) \
    // 	memcpy (dst, src, (dst##size = src##size) * sizeof (mp_limb_t))
    // 
    // 
    // /* Set errno and return an overflowing value with sign specified by
    //    NEGATIVE.  */
    // static FLOAT
    // overflow_value (int negative)
    // {
    //   __set_errno (ERANGE);
    //   FLOAT result = math_narrow_eval ((negative ? -MAX_VALUE : MAX_VALUE)
    // 				   * MAX_VALUE);
    //   return result;
    // }
    // 
    // 
    // /* Set errno and return an underflowing value with sign specified by
    //    NEGATIVE.  */
    // static FLOAT
    // underflow_value (int negative)
    // {
    //   __set_errno (ERANGE);
    //   FLOAT result = math_narrow_eval ((negative ? -MIN_VALUE : MIN_VALUE)
    // 				   * MIN_VALUE);
    //   return result;
    // }
    // 
    // 
    // /* Return a floating point number of the needed type according to the given
    //    multi-precision number after possible rounding.  */
    // static FLOAT
    // round_and_return (mp_limb_t *retval, intmax_t exponent, int negative,
    // 		  mp_limb_t round_limb, mp_size_t round_bit, int more_bits)
    // {
    //   int mode = get_rounding_mode ();
    // 
    //   if (exponent < MIN_EXP - 1)
    //     {
    //       if (exponent < MIN_EXP - 1 - MANT_DIG)
    // 	return underflow_value (negative);
    // 
    //       mp_size_t shift = MIN_EXP - 1 - exponent;
    //       bool is_tiny = true;
    //       bool old_half_bit = (round_limb & (((mp_limb_t) 1) << round_bit)) != 0;
    // 
    //       more_bits |= (round_limb & ((((mp_limb_t) 1) << round_bit) - 1)) != 0;
    //       if (shift == MANT_DIG)
    // 	/* This is a special case to handle the very seldom case where
    // 	   the mantissa will be empty after the shift.  */
    // 	{
    // 	  int i;
    // 
    // 	  round_limb = retval[RETURN_LIMB_SIZE - 1];
    // 	  round_bit = (MANT_DIG - 1) % BITS_PER_MP_LIMB;
    // 	  for (i = 0; i < RETURN_LIMB_SIZE - 1; ++i)
    // 	    more_bits |= retval[i] != 0;
    // 	  MPN_ZERO (retval, RETURN_LIMB_SIZE);
    // 	}
    //       else if (shift >= BITS_PER_MP_LIMB)
    // 	{
    // 	  int i;
    // 
    // 	  round_limb = retval[(shift - 1) / BITS_PER_MP_LIMB];
    // 	  round_bit = (shift - 1) % BITS_PER_MP_LIMB;
    // 	  for (i = 0; i < (shift - 1) / BITS_PER_MP_LIMB; ++i)
    // 	    more_bits |= retval[i] != 0;
    // 	  more_bits |= ((round_limb & ((((mp_limb_t) 1) << round_bit) - 1))
    // 			!= 0);
    // 
    // 	  /* __mpn_rshift requires 0 < shift < BITS_PER_MP_LIMB.  */
    // 	  if ((shift % BITS_PER_MP_LIMB) != 0)
    // 	    (void) __mpn_rshift (retval, &retval[shift / BITS_PER_MP_LIMB],
    // 			         RETURN_LIMB_SIZE - (shift / BITS_PER_MP_LIMB),
    // 			         shift % BITS_PER_MP_LIMB);
    // 	  else
    // 	    for (i = 0; i < RETURN_LIMB_SIZE - (shift / BITS_PER_MP_LIMB); i++)
    // 	      retval[i] = retval[i + (shift / BITS_PER_MP_LIMB)];
    // 	  MPN_ZERO (&retval[RETURN_LIMB_SIZE - (shift / BITS_PER_MP_LIMB)],
    // 		    shift / BITS_PER_MP_LIMB);
    // 	}
    //       else if (shift > 0)
    // 	{
    // 	  if (TININESS_AFTER_ROUNDING && shift == 1)
    // 	    {
    // 	      /* Whether the result counts as tiny depends on whether,
    // 		 after rounding to the normal precision, it still has
    // 		 a subnormal exponent.  */
    // 	      mp_limb_t retval_normal[RETURN_LIMB_SIZE];
    // 	      if (round_away (negative,
    // 			      (retval[0] & 1) != 0,
    // 			      (round_limb
    // 			       & (((mp_limb_t) 1) << round_bit)) != 0,
    // 			      (more_bits
    // 			       || ((round_limb
    // 				    & ((((mp_limb_t) 1) << round_bit) - 1))
    // 				   != 0)),
    // 			      mode))
    // 		{
    // 		  mp_limb_t cy = __mpn_add_1 (retval_normal, retval,
    // 					      RETURN_LIMB_SIZE, 1);
    // 
    // 		  if (((MANT_DIG % BITS_PER_MP_LIMB) == 0 && cy)
    // 		      || ((MANT_DIG % BITS_PER_MP_LIMB) != 0
    // 			  && ((retval_normal[RETURN_LIMB_SIZE - 1]
    // 			       & (((mp_limb_t) 1)
    // 				  << (MANT_DIG % BITS_PER_MP_LIMB)))
    // 			      != 0)))
    // 		    is_tiny = false;
    // 		}
    // 	    }
    // 	  round_limb = retval[0];
    // 	  round_bit = shift - 1;
    // 	  (void) __mpn_rshift (retval, retval, RETURN_LIMB_SIZE, shift);
    // 	}
    //       more_bits |= old_half_bit;
    //       /* This is a hook for the m68k long double format, where the
    // 	 exponent bias is the same for normalized and denormalized
    // 	 numbers.  */
    // #ifndef DENORM_EXP
    // # define DENORM_EXP (MIN_EXP - 2)
    // #endif
    //       exponent = DENORM_EXP;
    //       if (is_tiny
    // 	  && ((round_limb & (((mp_limb_t) 1) << round_bit)) != 0
    // 	      || more_bits
    // 	      || (round_limb & ((((mp_limb_t) 1) << round_bit) - 1)) != 0))
    // 	{
    // 	  __set_errno (ERANGE);
    // 	  FLOAT force_underflow = MIN_VALUE * MIN_VALUE;
    // 	  math_force_eval (force_underflow);
    // 	}
    //     }
    // 
    //   if (exponent >= MAX_EXP)
    //     goto overflow;
    // 
    //   bool half_bit = (round_limb & (((mp_limb_t) 1) << round_bit)) != 0;
    //   bool more_bits_nonzero
    //     = (more_bits
    //        || (round_limb & ((((mp_limb_t) 1) << round_bit) - 1)) != 0);
    //   if (round_away (negative,
    // 		  (retval[0] & 1) != 0,
    // 		  half_bit,
    // 		  more_bits_nonzero,
    // 		  mode))
    //     {
    //       mp_limb_t cy = __mpn_add_1 (retval, retval, RETURN_LIMB_SIZE, 1);
    // 
    //       if (((MANT_DIG % BITS_PER_MP_LIMB) == 0 && cy)
    // 	  || ((MANT_DIG % BITS_PER_MP_LIMB) != 0
    // 	      && (retval[RETURN_LIMB_SIZE - 1]
    // 		  & (((mp_limb_t) 1) << (MANT_DIG % BITS_PER_MP_LIMB))) != 0))
    // 	{
    // 	  ++exponent;
    // 	  (void) __mpn_rshift (retval, retval, RETURN_LIMB_SIZE, 1);
    // 	  retval[RETURN_LIMB_SIZE - 1]
    // 	    |= ((mp_limb_t) 1) << ((MANT_DIG - 1) % BITS_PER_MP_LIMB);
    // 	}
    //       else if (exponent == DENORM_EXP
    // 	       && (retval[RETURN_LIMB_SIZE - 1]
    // 		   & (((mp_limb_t) 1) << ((MANT_DIG - 1) % BITS_PER_MP_LIMB)))
    // 	       != 0)
    // 	  /* The number was denormalized but now normalized.  */
    // 	exponent = MIN_EXP - 1;
    //     }
    // 
    //   if (exponent >= MAX_EXP)
    //   overflow:
    //     return overflow_value (negative);
    // 
    //   if (half_bit || more_bits_nonzero)
    //     {
    //       FLOAT force_inexact = (FLOAT) 1 + MIN_VALUE;
    //       math_force_eval (force_inexact);
    //     }
    //   return MPN2FLOAT (retval, exponent, negative);
    // }
    // 
    // 
    // /* Read a multi-precision integer starting at STR with exactly DIGCNT digits
    //    into N.  Return the size of the number limbs in NSIZE at the first
    //    character od the string that is not part of the integer as the function
    //    value.  If the EXPONENT is small enough to be taken as an additional
    //    factor for the resulting number (see code) multiply by it.  */
    // static const STRING_TYPE *
    // str_to_mpn (const STRING_TYPE *str, int digcnt, mp_limb_t *n, mp_size_t *nsize,
    // 	    intmax_t *exponent
    // #ifndef USE_WIDE_CHAR
    // 	    , const char *decimal, size_t decimal_len, const char *thousands
    // #endif
    // 
    // 	    )
    // {
    //   /* Number of digits for actual limb.  */
    //   int cnt = 0;
    //   mp_limb_t low = 0;
    //   mp_limb_t start;
    // 
    //   *nsize = 0;
    //   assert (digcnt > 0);
    //   do
    //     {
    //       if (cnt == MAX_DIG_PER_LIMB)
    // 	{
    // 	  if (*nsize == 0)
    // 	    {
    // 	      n[0] = low;
    // 	      *nsize = 1;
    // 	    }
    // 	  else
    // 	    {
    // 	      mp_limb_t cy;
    // 	      cy = __mpn_mul_1 (n, n, *nsize, MAX_FAC_PER_LIMB);
    // 	      cy += __mpn_add_1 (n, n, *nsize, low);
    // 	      if (cy != 0)
    // 		{
    // 		  assert (*nsize < MPNSIZE);
    // 		  n[*nsize] = cy;
    // 		  ++(*nsize);
    // 		}
    // 	    }
    // 	  cnt = 0;
    // 	  low = 0;
    // 	}
    // 
    //       /* There might be thousands separators or radix characters in
    // 	 the string.  But these all can be ignored because we know the
    // 	 format of the number is correct and we have an exact number
    // 	 of characters to read.  */
    // #ifdef USE_WIDE_CHAR
    //       if (*str < L'0' || *str > L'9')
    // 	++str;
    // #else
    //       if (*str < '0' || *str > '9')
    // 	{
    // 	  int inner = 0;
    // 	  if (thousands != NULL && *str == *thousands
    // 	      && ({ for (inner = 1; thousands[inner] != '\0'; ++inner)
    // 		      if (thousands[inner] != str[inner])
    // 			break;
    // 		    thousands[inner] == '\0'; }))
    // 	    str += inner;
    // 	  else
    // 	    str += decimal_len;
    // 	}
    // #endif
    //       low = low * 10 + *str++ - L_('0');
    //       ++cnt;
    //     }
    //   while (--digcnt > 0);
    // 
    //   if (*exponent > 0 && *exponent <= MAX_DIG_PER_LIMB - cnt)
    //     {
    //       low *= _tens_in_limb[*exponent];
    //       start = _tens_in_limb[cnt + *exponent];
    //       *exponent = 0;
    //     }
    //   else
    //     start = _tens_in_limb[cnt];
    // 
    //   if (*nsize == 0)
    //     {
    //       n[0] = low;
    //       *nsize = 1;
    //     }
    //   else
    //     {
    //       mp_limb_t cy;
    //       cy = __mpn_mul_1 (n, n, *nsize, start);
    //       cy += __mpn_add_1 (n, n, *nsize, low);
    //       if (cy != 0)
    // 	{
    // 	  assert (*nsize < MPNSIZE);
    // 	  n[(*nsize)++] = cy;
    // 	}
    //     }
    // 
    //   return str;
    // }
    // 
    // 
    // /* Shift {PTR, SIZE} COUNT bits to the left, and fill the vacated bits
    //    with the COUNT most significant bits of LIMB.
    // 
    //    Implemented as a macro, so that __builtin_constant_p works even at -O0.
    // 
    //    Tege doesn't like this macro so I have to write it here myself. :)
    //    --drepper */
    // #define __mpn_lshift_1(ptr, size, count, limb) \
    //   do									\
    //     {									\
    //       mp_limb_t *__ptr = (ptr);						\
    //       if (__builtin_constant_p (count) && count == BITS_PER_MP_LIMB)	\
    // 	{								\
    // 	  mp_size_t i;							\
    // 	  for (i = (size) - 1; i > 0; --i)				\
    // 	    __ptr[i] = __ptr[i - 1];					\
    // 	  __ptr[0] = (limb);						\
    // 	}								\
    //       else								\
    // 	{								\
    // 	  /* We assume count > 0 && count < BITS_PER_MP_LIMB here.  */	\
    // 	  unsigned int __count = (count);				\
    // 	  (void) __mpn_lshift (__ptr, __ptr, size, __count);		\
    // 	  __ptr[0] |= (limb) >> (BITS_PER_MP_LIMB - __count);		\
    // 	}								\
    //     }									\
    //   while (0)
    // 
    // 
    // #define INTERNAL(x) INTERNAL1(x)
    // #define INTERNAL1(x) __##x##_internal
    // #ifndef ____STRTOF_INTERNAL
    // # define ____STRTOF_INTERNAL INTERNAL (__STRTOF)
    // #endif
    // 
    // /* This file defines a function to check for correct grouping.  */
    // #include "grouping.h"
    // 
    // 
    // /* Return a floating point number with the value of the given string NPTR.
    //    Set *ENDPTR to the character after the last used one.  If the number is
    //    smaller than the smallest representable number, set `errno' to ERANGE and
    //    return 0.0.  If the number is too big to be represented, set `errno' to
    //    ERANGE and return HUGE_VAL with the appropriate sign.  */
    // FLOAT
    // ____STRTOF_INTERNAL (const STRING_TYPE *nptr, STRING_TYPE **endptr, int group,
    // 		     locale_t loc)
    // {
    //   int negative;			/* The sign of the number.  */
    //   MPN_VAR (num);		/* MP representation of the number.  */
    //   intmax_t exponent;		/* Exponent of the number.  */
    // 
    //   /* Numbers starting `0X' or `0x' have to be processed with base 16.  */
    //   int base = 10;
    // 
    //   /* When we have to compute fractional digits we form a fraction with a
    //      second multi-precision number (and we sometimes need a second for
    //      temporary results).  */
    //   MPN_VAR (den);
    // 
    //   /* Representation for the return value.  */
    //   mp_limb_t retval[RETURN_LIMB_SIZE];
    //   /* Number of bits currently in result value.  */
    //   int bits;
    // 
    //   /* Running pointer after the last character processed in the string.  */
    //   const STRING_TYPE *cp, *tp;
    //   /* Start of significant part of the number.  */
    //   const STRING_TYPE *startp, *start_of_digits;
    //   /* Points at the character following the integer and fractional digits.  */
    //   const STRING_TYPE *expp;
    //   /* Total number of digit and number of digits in integer part.  */
    //   size_t dig_no, int_no, lead_zero;
    //   /* Contains the last character read.  */
    //   CHAR_TYPE c;
    // 
    // /* We should get wint_t from <stddef.h>, but not all GCC versions define it
    //    there.  So define it ourselves if it remains undefined.  */
    // #ifndef _WINT_T
    //   typedef unsigned int wint_t;
    // #endif
    //   /* The radix character of the current locale.  */
    // #ifdef USE_WIDE_CHAR
    //   wchar_t decimal;
    // #else
    //   const char *decimal;
    //   size_t decimal_len;
    // #endif
    //   /* The thousands character of the current locale.  */
    // #ifdef USE_WIDE_CHAR
    //   wchar_t thousands = L'\0';
    // #else
    //   const char *thousands = NULL;
    // #endif
    //   /* The numeric grouping specification of the current locale,
    //      in the format described in <locale.h>.  */
    //   const char *grouping;
    //   /* Used in several places.  */
    //   int cnt;
    // 
    //   struct __locale_data *current = loc->__locales[LC_NUMERIC];
    // 
    //   if (__glibc_unlikely (group))
    //     {
    //       grouping = _NL_CURRENT (LC_NUMERIC, GROUPING);
    //       if (*grouping <= 0 || *grouping == CHAR_MAX)
    // 	grouping = NULL;
    //       else
    // 	{
    // 	  /* Figure out the thousands separator character.  */
    // #ifdef USE_WIDE_CHAR
    // 	  thousands = _NL_CURRENT_WORD (LC_NUMERIC,
    // 					_NL_NUMERIC_THOUSANDS_SEP_WC);
    // 	  if (thousands == L'\0')
    // 	    grouping = NULL;
    // #else
    // 	  thousands = _NL_CURRENT (LC_NUMERIC, THOUSANDS_SEP);
    // 	  if (*thousands == '\0')
    // 	    {
    // 	      thousands = NULL;
    // 	      grouping = NULL;
    // 	    }
    // #endif
    // 	}
    //     }
    //   else
    //     grouping = NULL;
    // 
    //   /* Find the locale's decimal point character.  */
    // #ifdef USE_WIDE_CHAR
    //   decimal = _NL_CURRENT_WORD (LC_NUMERIC, _NL_NUMERIC_DECIMAL_POINT_WC);
    //   assert (decimal != L'\0');
    // # define decimal_len 1
    // #else
    //   decimal = _NL_CURRENT (LC_NUMERIC, DECIMAL_POINT);
    //   decimal_len = strlen (decimal);
    //   assert (decimal_len > 0);
    // #endif
    // 
    //   /* Prepare number representation.  */
    //   exponent = 0;
    //   negative = 0;
    //   bits = 0;
    // 
    //   /* Parse string to get maximal legal prefix.  We need the number of
    //      characters of the integer part, the fractional part and the exponent.  */
    //   cp = nptr - 1;
    //   /* Ignore leading white space.  */
    //   do
    //     c = *++cp;
    //   while (ISSPACE (c));
    // 
    //   /* Get sign of the result.  */
    //   if (c == L_('-'))
    //     {
    //       negative = 1;
    //       c = *++cp;
    //     }
    //   else if (c == L_('+'))
    //     c = *++cp;
    // 
    //   /* Return 0.0 if no legal string is found.
    //      No character is used even if a sign was found.  */
    // #ifdef USE_WIDE_CHAR
    //   if (c == (wint_t) decimal
    //       && (wint_t) cp[1] >= L'0' && (wint_t) cp[1] <= L'9')
    //     {
    //       /* We accept it.  This funny construct is here only to indent
    // 	 the code correctly.  */
    //     }
    // #else
    //   for (cnt = 0; decimal[cnt] != '\0'; ++cnt)
    //     if (cp[cnt] != decimal[cnt])
    //       break;
    //   if (decimal[cnt] == '\0' && cp[cnt] >= '0' && cp[cnt] <= '9')
    //     {
    //       /* We accept it.  This funny construct is here only to indent
    // 	 the code correctly.  */
    //     }
    // #endif
    //   else if (c < L_('0') || c > L_('9'))
    //     {
    //       /* Check for `INF' or `INFINITY'.  */
    //       CHAR_TYPE lowc = TOLOWER_C (c);
    // 
    //       if (lowc == L_('i') && STRNCASECMP (cp, L_("inf"), 3) == 0)
    // 	{
    // 	  /* Return +/- infinity.  */
    // 	  if (endptr != NULL)
    // 	    *endptr = (STRING_TYPE *)
    // 		      (cp + (STRNCASECMP (cp + 3, L_("inity"), 5) == 0
    // 			     ? 8 : 3));
    // 
    // 	  return negative ? -FLOAT_HUGE_VAL : FLOAT_HUGE_VAL;
    // 	}
    // 
    //       if (lowc == L_('n') && STRNCASECMP (cp, L_("nan"), 3) == 0)
    // 	{
    // 	  /* Return NaN.  */
    // 	  FLOAT retval = NAN;
    // 
    // 	  cp += 3;
    // 
    // 	  /* Match `(n-char-sequence-digit)'.  */
    // 	  if (*cp == L_('('))
    // 	    {
    // 	      const STRING_TYPE *startp = cp;
    // 	      STRING_TYPE *endp;
    // 	      retval = STRTOF_NAN (cp + 1, &endp, L_(')'));
    // 	      if (*endp == L_(')'))
    // 		/* Consume the closing parenthesis.  */
    // 		cp = endp + 1;
    // 	      else
    // 		/* Only match the NAN part.  */
    // 		cp = startp;
    // 	    }
    // 
    // 	  if (endptr != NULL)
    // 	    *endptr = (STRING_TYPE *) cp;
    // 
    // 	  return negative ? -retval : retval;
    // 	}
    // 
    //       /* It is really a text we do not recognize.  */
    //       RETURN (0.0, nptr);
    //     }
    // 
    //   /* First look whether we are faced with a hexadecimal number.  */
    //   if (c == L_('0') && TOLOWER (cp[1]) == L_('x'))
    //     {
    //       /* Okay, it is a hexa-decimal number.  Remember this and skip
    // 	 the characters.  BTW: hexadecimal numbers must not be
    // 	 grouped.  */
    //       base = 16;
    //       cp += 2;
    //       c = *cp;
    //       grouping = NULL;
    //     }
    // 
    //   /* Record the start of the digits, in case we will check their grouping.  */
    //   start_of_digits = startp = cp;
    // 
    //   /* Ignore leading zeroes.  This helps us to avoid useless computations.  */
    // #ifdef USE_WIDE_CHAR
    //   while (c == L'0' || ((wint_t) thousands != L'\0' && c == (wint_t) thousands))
    //     c = *++cp;
    // #else
    //   if (__glibc_likely (thousands == NULL))
    //     while (c == '0')
    //       c = *++cp;
    //   else
    //     {
    //       /* We also have the multibyte thousands string.  */
    //       while (1)
    // 	{
    // 	  if (c != '0')
    // 	    {
    // 	      for (cnt = 0; thousands[cnt] != '\0'; ++cnt)
    // 		if (thousands[cnt] != cp[cnt])
    // 		  break;
    // 	      if (thousands[cnt] != '\0')
    // 		break;
    // 	      cp += cnt - 1;
    // 	    }
    // 	  c = *++cp;
    // 	}
    //     }
    // #endif
    // 
    //   /* If no other digit but a '0' is found the result is 0.0.
    //      Return current read pointer.  */
    //   CHAR_TYPE lowc = TOLOWER (c);
    //   if (!((c >= L_('0') && c <= L_('9'))
    // 	|| (base == 16 && lowc >= L_('a') && lowc <= L_('f'))
    // 	|| (
    // #ifdef USE_WIDE_CHAR
    // 	    c == (wint_t) decimal
    // #else
    // 	    ({ for (cnt = 0; decimal[cnt] != '\0'; ++cnt)
    // 		 if (decimal[cnt] != cp[cnt])
    // 		   break;
    // 	       decimal[cnt] == '\0'; })
    // #endif
    // 	    /* '0x.' alone is not a valid hexadecimal number.
    // 	       '.' alone is not valid either, but that has been checked
    // 	       already earlier.  */
    // 	    && (base != 16
    // 		|| cp != start_of_digits
    // 		|| (cp[decimal_len] >= L_('0') && cp[decimal_len] <= L_('9'))
    // 		|| ({ CHAR_TYPE lo = TOLOWER (cp[decimal_len]);
    // 		      lo >= L_('a') && lo <= L_('f'); })))
    // 	|| (base == 16 && (cp != start_of_digits
    // 			   && lowc == L_('p')))
    // 	|| (base != 16 && lowc == L_('e'))))
    //     {
    // #ifdef USE_WIDE_CHAR
    //       tp = __correctly_grouped_prefixwc (start_of_digits, cp, thousands,
    // 					 grouping);
    // #else
    //       tp = __correctly_grouped_prefixmb (start_of_digits, cp, thousands,
    // 					 grouping);
    // #endif
    //       /* If TP is at the start of the digits, there was no correctly
    // 	 grouped prefix of the string; so no number found.  */
    //       RETURN (negative ? -0.0 : 0.0,
    // 	      tp == start_of_digits ? (base == 16 ? cp - 1 : nptr) : tp);
    //     }
    // 
    //   /* Remember first significant digit and read following characters until the
    //      decimal point, exponent character or any non-FP number character.  */
    //   startp = cp;
    //   dig_no = 0;
    //   while (1)
    //     {
    //       if ((c >= L_('0') && c <= L_('9'))
    // 	  || (base == 16
    // 	      && ({ CHAR_TYPE lo = TOLOWER (c);
    // 		    lo >= L_('a') && lo <= L_('f'); })))
    // 	++dig_no;
    //       else
    // 	{
    // #ifdef USE_WIDE_CHAR
    // 	  if (__builtin_expect ((wint_t) thousands == L'\0', 1)
    // 	      || c != (wint_t) thousands)
    // 	    /* Not a digit or separator: end of the integer part.  */
    // 	    break;
    // #else
    // 	  if (__glibc_likely (thousands == NULL))
    // 	    break;
    // 	  else
    // 	    {
    // 	      for (cnt = 0; thousands[cnt] != '\0'; ++cnt)
    // 		if (thousands[cnt] != cp[cnt])
    // 		  break;
    // 	      if (thousands[cnt] != '\0')
    // 		break;
    // 	      cp += cnt - 1;
    // 	    }
    // #endif
    // 	}
    //       c = *++cp;
    //     }
    // 
    //   if (__builtin_expect (grouping != NULL, 0) && cp > start_of_digits)
    //     {
    //       /* Check the grouping of the digits.  */
    // #ifdef USE_WIDE_CHAR
    //       tp = __correctly_grouped_prefixwc (start_of_digits, cp, thousands,
    // 					 grouping);
    // #else
    //       tp = __correctly_grouped_prefixmb (start_of_digits, cp, thousands,
    // 					 grouping);
    // #endif
    //       if (cp != tp)
    // 	{
    // 	  /* Less than the entire string was correctly grouped.  */
    // 
    // 	  if (tp == start_of_digits)
    // 	    /* No valid group of numbers at all: no valid number.  */
    // 	    RETURN (0.0, nptr);
    // 
    // 	  if (tp < startp)
    // 	    /* The number is validly grouped, but consists
    // 	       only of zeroes.  The whole value is zero.  */
    // 	    RETURN (negative ? -0.0 : 0.0, tp);
    // 
    // 	  /* Recompute DIG_NO so we won't read more digits than
    // 	     are properly grouped.  */
    // 	  cp = tp;
    // 	  dig_no = 0;
    // 	  for (tp = startp; tp < cp; ++tp)
    // 	    if (*tp >= L_('0') && *tp <= L_('9'))
    // 	      ++dig_no;
    // 
    // 	  int_no = dig_no;
    // 	  lead_zero = 0;
    // 
    // 	  goto number_parsed;
    // 	}
    //     }
    // 
    //   /* We have the number of digits in the integer part.  Whether these
    //      are all or any is really a fractional digit will be decided
    //      later.  */
    //   int_no = dig_no;
    //   lead_zero = int_no == 0 ? (size_t) -1 : 0;
    // 
    //   /* Read the fractional digits.  A special case are the 'american
    //      style' numbers like `16.' i.e. with decimal point but without
    //      trailing digits.  */
    //   if (
    // #ifdef USE_WIDE_CHAR
    //       c == (wint_t) decimal
    // #else
    //       ({ for (cnt = 0; decimal[cnt] != '\0'; ++cnt)
    // 	   if (decimal[cnt] != cp[cnt])
    // 	     break;
    // 	 decimal[cnt] == '\0'; })
    // #endif
    //       )
    //     {
    //       cp += decimal_len;
    //       c = *cp;
    //       while ((c >= L_('0') && c <= L_('9'))
    // 	     || (base == 16 && ({ CHAR_TYPE lo = TOLOWER (c);
    // 				  lo >= L_('a') && lo <= L_('f'); })))
    // 	{
    // 	  if (c != L_('0') && lead_zero == (size_t) -1)
    // 	    lead_zero = dig_no - int_no;
    // 	  ++dig_no;
    // 	  c = *++cp;
    // 	}
    //     }
    //   assert (dig_no <= (uintmax_t) INTMAX_MAX);
    // 
    //   /* Remember start of exponent (if any).  */
    //   expp = cp;
    // 
    //   /* Read exponent.  */
    //   lowc = TOLOWER (c);
    //   if ((base == 16 && lowc == L_('p'))
    //       || (base != 16 && lowc == L_('e')))
    //     {
    //       int exp_negative = 0;
    // 
    //       c = *++cp;
    //       if (c == L_('-'))
    // 	{
    // 	  exp_negative = 1;
    // 	  c = *++cp;
    // 	}
    //       else if (c == L_('+'))
    // 	c = *++cp;
    // 
    //       if (c >= L_('0') && c <= L_('9'))
    // 	{
    // 	  intmax_t exp_limit;
    // 
    // 	  /* Get the exponent limit. */
    // 	  if (base == 16)
    // 	    {
    // 	      if (exp_negative)
    // 		{
    // 		  assert (int_no <= (uintmax_t) (INTMAX_MAX
    // 						 + MIN_EXP - MANT_DIG) / 4);
    // 		  exp_limit = -MIN_EXP + MANT_DIG + 4 * (intmax_t) int_no;
    // 		}
    // 	      else
    // 		{
    // 		  if (int_no)
    // 		    {
    // 		      assert (lead_zero == 0
    // 			      && int_no <= (uintmax_t) INTMAX_MAX / 4);
    // 		      exp_limit = MAX_EXP - 4 * (intmax_t) int_no + 3;
    // 		    }
    // 		  else if (lead_zero == (size_t) -1)
    // 		    {
    // 		      /* The number is zero and this limit is
    // 			 arbitrary.  */
    // 		      exp_limit = MAX_EXP + 3;
    // 		    }
    // 		  else
    // 		    {
    // 		      assert (lead_zero
    // 			      <= (uintmax_t) (INTMAX_MAX - MAX_EXP - 3) / 4);
    // 		      exp_limit = (MAX_EXP
    // 				   + 4 * (intmax_t) lead_zero
    // 				   + 3);
    // 		    }
    // 		}
    // 	    }
    // 	  else
    // 	    {
    // 	      if (exp_negative)
    // 		{
    // 		  assert (int_no
    // 			  <= (uintmax_t) (INTMAX_MAX + MIN_10_EXP - MANT_DIG));
    // 		  exp_limit = -MIN_10_EXP + MANT_DIG + (intmax_t) int_no;
    // 		}
    // 	      else
    // 		{
    // 		  if (int_no)
    // 		    {
    // 		      assert (lead_zero == 0
    // 			      && int_no <= (uintmax_t) INTMAX_MAX);
    // 		      exp_limit = MAX_10_EXP - (intmax_t) int_no + 1;
    // 		    }
    // 		  else if (lead_zero == (size_t) -1)
    // 		    {
    // 		      /* The number is zero and this limit is
    // 			 arbitrary.  */
    // 		      exp_limit = MAX_10_EXP + 1;
    // 		    }
    // 		  else
    // 		    {
    // 		      assert (lead_zero
    // 			      <= (uintmax_t) (INTMAX_MAX - MAX_10_EXP - 1));
    // 		      exp_limit = MAX_10_EXP + (intmax_t) lead_zero + 1;
    // 		    }
    // 		}
    // 	    }
    // 
    // 	  if (exp_limit < 0)
    // 	    exp_limit = 0;
    // 
    // 	  do
    // 	    {
    // 	      if (__builtin_expect ((exponent > exp_limit / 10
    // 				     || (exponent == exp_limit / 10
    // 					 && c - L_('0') > exp_limit % 10)), 0))
    // 		/* The exponent is too large/small to represent a valid
    // 		   number.  */
    // 		{
    // 		  FLOAT result;
    // 
    // 		  /* We have to take care for special situation: a joker
    // 		     might have written "0.0e100000" which is in fact
    // 		     zero.  */
    // 		  if (lead_zero == (size_t) -1)
    // 		    result = negative ? -0.0 : 0.0;
    // 		  else
    // 		    {
    // 		      /* Overflow or underflow.  */
    // 		      result = (exp_negative
    // 				? underflow_value (negative)
    // 				: overflow_value (negative));
    // 		    }
    // 
    // 		  /* Accept all following digits as part of the exponent.  */
    // 		  do
    // 		    ++cp;
    // 		  while (*cp >= L_('0') && *cp <= L_('9'));
    // 
    // 		  RETURN (result, cp);
    // 		  /* NOTREACHED */
    // 		}
    // 
    // 	      exponent *= 10;
    // 	      exponent += c - L_('0');
    // 
    // 	      c = *++cp;
    // 	    }
    // 	  while (c >= L_('0') && c <= L_('9'));
    // 
    // 	  if (exp_negative)
    // 	    exponent = -exponent;
    // 	}
    //       else
    // 	cp = expp;
    //     }
    // 
    //   /* We don't want to have to work with trailing zeroes after the radix.  */
    //   if (dig_no > int_no)
    //     {
    //       while (expp[-1] == L_('0'))
    // 	{
    // 	  --expp;
    // 	  --dig_no;
    // 	}
    //       assert (dig_no >= int_no);
    //     }
    // 
    //   if (dig_no == int_no && dig_no > 0 && exponent < 0)
    //     do
    //       {
    // 	while (! (base == 16 ? ISXDIGIT (expp[-1]) : ISDIGIT (expp[-1])))
    // 	  --expp;
    // 
    // 	if (expp[-1] != L_('0'))
    // 	  break;
    // 
    // 	--expp;
    // 	--dig_no;
    // 	--int_no;
    // 	exponent += base == 16 ? 4 : 1;
    //       }
    //     while (dig_no > 0 && exponent < 0);
    // 
    //  number_parsed:
    // 
    //   /* The whole string is parsed.  Store the address of the next character.  */
    //   if (endptr)
    //     *endptr = (STRING_TYPE *) cp;
    // 
    //   if (dig_no == 0)
    //     return negative ? -0.0 : 0.0;
    // 
    //   if (lead_zero)
    //     {
    //       /* Find the decimal point */
    // #ifdef USE_WIDE_CHAR
    //       while (*startp != decimal)
    // 	++startp;
    // #else
    //       while (1)
    // 	{
    // 	  if (*startp == decimal[0])
    // 	    {
    // 	      for (cnt = 1; decimal[cnt] != '\0'; ++cnt)
    // 		if (decimal[cnt] != startp[cnt])
    // 		  break;
    // 	      if (decimal[cnt] == '\0')
    // 		break;
    // 	    }
    // 	  ++startp;
    // 	}
    // #endif
    //       startp += lead_zero + decimal_len;
    //       assert (lead_zero <= (base == 16
    // 			    ? (uintmax_t) INTMAX_MAX / 4
    // 			    : (uintmax_t) INTMAX_MAX));
    //       assert (lead_zero <= (base == 16
    // 			    ? ((uintmax_t) exponent
    // 			       - (uintmax_t) INTMAX_MIN) / 4
    // 			    : ((uintmax_t) exponent - (uintmax_t) INTMAX_MIN)));
    //       exponent -= base == 16 ? 4 * (intmax_t) lead_zero : (intmax_t) lead_zero;
    //       dig_no -= lead_zero;
    //     }
    // 
    //   /* If the BASE is 16 we can use a simpler algorithm.  */
    //   if (base == 16)
    //     {
    //       static const int nbits[16] = { 0, 1, 2, 2, 3, 3, 3, 3,
    // 				     4, 4, 4, 4, 4, 4, 4, 4 };
    //       int idx = (MANT_DIG - 1) / BITS_PER_MP_LIMB;
    //       int pos = (MANT_DIG - 1) % BITS_PER_MP_LIMB;
    //       mp_limb_t val;
    // 
    //       while (!ISXDIGIT (*startp))
    // 	++startp;
    //       while (*startp == L_('0'))
    // 	++startp;
    //       if (ISDIGIT (*startp))
    // 	val = *startp++ - L_('0');
    //       else
    // 	val = 10 + TOLOWER (*startp++) - L_('a');
    //       bits = nbits[val];
    //       /* We cannot have a leading zero.  */
    //       assert (bits != 0);
    // 
    //       if (pos + 1 >= 4 || pos + 1 >= bits)
    // 	{
    // 	  /* We don't have to care for wrapping.  This is the normal
    // 	     case so we add the first clause in the `if' expression as
    // 	     an optimization.  It is a compile-time constant and so does
    // 	     not cost anything.  */
    // 	  retval[idx] = val << (pos - bits + 1);
    // 	  pos -= bits;
    // 	}
    //       else
    // 	{
    // 	  retval[idx--] = val >> (bits - pos - 1);
    // 	  retval[idx] = val << (BITS_PER_MP_LIMB - (bits - pos - 1));
    // 	  pos = BITS_PER_MP_LIMB - 1 - (bits - pos - 1);
    // 	}
    // 
    //       /* Adjust the exponent for the bits we are shifting in.  */
    //       assert (int_no <= (uintmax_t) (exponent < 0
    // 				     ? (INTMAX_MAX - bits + 1) / 4
    // 				     : (INTMAX_MAX - exponent - bits + 1) / 4));
    //       exponent += bits - 1 + ((intmax_t) int_no - 1) * 4;
    // 
    //       while (--dig_no > 0 && idx >= 0)
    // 	{
    // 	  if (!ISXDIGIT (*startp))
    // 	    startp += decimal_len;
    // 	  if (ISDIGIT (*startp))
    // 	    val = *startp++ - L_('0');
    // 	  else
    // 	    val = 10 + TOLOWER (*startp++) - L_('a');
    // 
    // 	  if (pos + 1 >= 4)
    // 	    {
    // 	      retval[idx] |= val << (pos - 4 + 1);
    // 	      pos -= 4;
    // 	    }
    // 	  else
    // 	    {
    // 	      retval[idx--] |= val >> (4 - pos - 1);
    // 	      val <<= BITS_PER_MP_LIMB - (4 - pos - 1);
    // 	      if (idx < 0)
    // 		{
    // 		  int rest_nonzero = 0;
    // 		  while (--dig_no > 0)
    // 		    {
    // 		      if (*startp != L_('0'))
    // 			{
    // 			  rest_nonzero = 1;
    // 			  break;
    // 			}
    // 		      startp++;
    // 		    }
    // 		  return round_and_return (retval, exponent, negative, val,
    // 					   BITS_PER_MP_LIMB - 1, rest_nonzero);
    // 		}
    // 
    // 	      retval[idx] = val;
    // 	      pos = BITS_PER_MP_LIMB - 1 - (4 - pos - 1);
    // 	    }
    // 	}
    // 
    //       /* We ran out of digits.  */
    //       MPN_ZERO (retval, idx);
    // 
    //       return round_and_return (retval, exponent, negative, 0, 0, 0);
    //     }
    // 
    //   /* Now we have the number of digits in total and the integer digits as well
    //      as the exponent and its sign.  We can decide whether the read digits are
    //      really integer digits or belong to the fractional part; i.e. we normalize
    //      123e-2 to 1.23.  */
    //   {
    //     intmax_t incr = (exponent < 0
    // 		     ? MAX (-(intmax_t) int_no, exponent)
    // 		     : MIN ((intmax_t) dig_no - (intmax_t) int_no, exponent));
    //     int_no += incr;
    //     exponent -= incr;
    //   }
    // 
    //   if (__glibc_unlikely (exponent > MAX_10_EXP + 1 - (intmax_t) int_no))
    //     return overflow_value (negative);
    // 
    //   /* 10^(MIN_10_EXP-1) is not normal.  Thus, 10^(MIN_10_EXP-1) /
    //      2^MANT_DIG is below half the least subnormal, so anything with a
    //      base-10 exponent less than the base-10 exponent (which is
    //      MIN_10_EXP - 1 - ceil(MANT_DIG*log10(2))) of that value
    //      underflows.  DIG is floor((MANT_DIG-1)log10(2)), so an exponent
    //      below MIN_10_EXP - (DIG + 3) underflows.  But EXPONENT is
    //      actually an exponent multiplied only by a fractional part, not an
    //      integer part, so an exponent below MIN_10_EXP - (DIG + 2)
    //      underflows.  */
    //   if (__glibc_unlikely (exponent < MIN_10_EXP - (DIG + 2)))
    //     return underflow_value (negative);
    // 
    //   if (int_no > 0)
    //     {
    //       /* Read the integer part as a multi-precision number to NUM.  */
    //       startp = str_to_mpn (startp, int_no, num, &numsize, &exponent
    // #ifndef USE_WIDE_CHAR
    // 			   , decimal, decimal_len, thousands
    // #endif
    // 			   );
    // 
    //       if (exponent > 0)
    // 	{
    // 	  /* We now multiply the gained number by the given power of ten.  */
    // 	  mp_limb_t *psrc = num;
    // 	  mp_limb_t *pdest = den;
    // 	  int expbit = 1;
    // 	  const struct mp_power *ttab = &_fpioconst_pow10[0];
    // 
    // 	  do
    // 	    {
    // 	      if ((exponent & expbit) != 0)
    // 		{
    // 		  size_t size = ttab->arraysize - _FPIO_CONST_OFFSET;
    // 		  mp_limb_t cy;
    // 		  exponent ^= expbit;
    // 
    // 		  /* FIXME: not the whole multiplication has to be
    // 		     done.  If we have the needed number of bits we
    // 		     only need the information whether more non-zero
    // 		     bits follow.  */
    // 		  if (numsize >= ttab->arraysize - _FPIO_CONST_OFFSET)
    // 		    cy = __mpn_mul (pdest, psrc, numsize,
    // 				    &__tens[ttab->arrayoff
    // 					   + _FPIO_CONST_OFFSET],
    // 				    size);
    // 		  else
    // 		    cy = __mpn_mul (pdest, &__tens[ttab->arrayoff
    // 						  + _FPIO_CONST_OFFSET],
    // 				    size, psrc, numsize);
    // 		  numsize += size;
    // 		  if (cy == 0)
    // 		    --numsize;
    // 		  (void) SWAP (psrc, pdest);
    // 		}
    // 	      expbit <<= 1;
    // 	      ++ttab;
    // 	    }
    // 	  while (exponent != 0);
    // 
    // 	  if (psrc == den)
    // 	    memcpy (num, den, numsize * sizeof (mp_limb_t));
    // 	}
    // 
    //       /* Determine how many bits of the result we already have.  */
    //       bits = stdc_leading_zeros (num[numsize - 1]);
    //       bits = numsize * BITS_PER_MP_LIMB - bits;
    // 
    //       /* Now we know the exponent of the number in base two.
    // 	 Check it against the maximum possible exponent.  */
    //       if (__glibc_unlikely (bits > MAX_EXP))
    // 	return overflow_value (negative);
    // 
    //       /* We have already the first BITS bits of the result.  Together with
    // 	 the information whether more non-zero bits follow this is enough
    // 	 to determine the result.  */
    //       if (bits > MANT_DIG)
    // 	{
    // 	  int i;
    // 	  const mp_size_t least_idx = (bits - MANT_DIG) / BITS_PER_MP_LIMB;
    // 	  const mp_size_t least_bit = (bits - MANT_DIG) % BITS_PER_MP_LIMB;
    // 	  const mp_size_t round_idx = least_bit == 0 ? least_idx - 1
    // 						     : least_idx;
    // 	  const mp_size_t round_bit = least_bit == 0 ? BITS_PER_MP_LIMB - 1
    // 						     : least_bit - 1;
    // 
    // 	  if (least_bit == 0)
    // 	    memcpy (retval, &num[least_idx],
    // 		    RETURN_LIMB_SIZE * sizeof (mp_limb_t));
    // 	  else
    // 	    {
    // 	      for (i = least_idx; i < numsize - 1; ++i)
    // 		retval[i - least_idx] = (num[i] >> least_bit)
    // 					| (num[i + 1]
    // 					   << (BITS_PER_MP_LIMB - least_bit));
    // 	      if (i - least_idx < RETURN_LIMB_SIZE)
    // 		retval[RETURN_LIMB_SIZE - 1] = num[i] >> least_bit;
    // 	    }
    // 
    // 	  /* Check whether any limb beside the ones in RETVAL are non-zero.  */
    // 	  for (i = 0; num[i] == 0; ++i)
    // 	    ;
    // 
    // 	  return round_and_return (retval, bits - 1, negative,
    // 				   num[round_idx], round_bit,
    // 				   int_no < dig_no || i < round_idx);
    // 	  /* NOTREACHED */
    // 	}
    //       else if (dig_no == int_no)
    // 	{
    // 	  const mp_size_t target_bit = (MANT_DIG - 1) % BITS_PER_MP_LIMB;
    // 	  const mp_size_t is_bit = (bits - 1) % BITS_PER_MP_LIMB;
    // 
    // 	  if (target_bit == is_bit)
    // 	    {
    // 	      memcpy (&retval[RETURN_LIMB_SIZE - numsize], num,
    // 		      numsize * sizeof (mp_limb_t));
    // 	      /* FIXME: the following loop can be avoided if we assume a
    // 		 maximal MANT_DIG value.  */
    // 	      MPN_ZERO (retval, RETURN_LIMB_SIZE - numsize);
    // 	    }
    // 	  else if (target_bit > is_bit)
    // 	    {
    // 	      (void) __mpn_lshift (&retval[RETURN_LIMB_SIZE - numsize],
    // 				   num, numsize, target_bit - is_bit);
    // 	      /* FIXME: the following loop can be avoided if we assume a
    // 		 maximal MANT_DIG value.  */
    // 	      MPN_ZERO (retval, RETURN_LIMB_SIZE - numsize);
    // 	    }
    // 	  else
    // 	    {
    // 	      mp_limb_t cy;
    // 	      assert (numsize < RETURN_LIMB_SIZE);
    // 
    // 	      cy = __mpn_rshift (&retval[RETURN_LIMB_SIZE - numsize],
    // 				 num, numsize, is_bit - target_bit);
    // 	      retval[RETURN_LIMB_SIZE - numsize - 1] = cy;
    // 	      /* FIXME: the following loop can be avoided if we assume a
    // 		 maximal MANT_DIG value.  */
    // 	      MPN_ZERO (retval, RETURN_LIMB_SIZE - numsize - 1);
    // 	    }
    // 
    // 	  return round_and_return (retval, bits - 1, negative, 0, 0, 0);
    // 	  /* NOTREACHED */
    // 	}
    // 
    //       /* Store the bits we already have.  */
    //       memcpy (retval, num, numsize * sizeof (mp_limb_t));
    // #if RETURN_LIMB_SIZE > 1
    //       if (numsize < RETURN_LIMB_SIZE)
    // # if RETURN_LIMB_SIZE == 2
    // 	retval[numsize] = 0;
    // # else
    // 	MPN_ZERO (retval + numsize, RETURN_LIMB_SIZE - numsize);
    // # endif
    // #endif
    //     }
    // 
    //   /* We have to compute at least some of the fractional digits.  */
    //   {
    //     /* We construct a fraction and the result of the division gives us
    //        the needed digits.  The denominator is 1.0 multiplied by the
    //        exponent of the lowest digit; i.e. 0.123 gives 123 / 1000 and
    //        123e-6 gives 123 / 1000000.  */
    // 
    //     int expbit;
    //     int neg_exp;
    //     int more_bits;
    //     int need_frac_digits;
    //     mp_limb_t cy;
    //     mp_limb_t *psrc = den;
    //     mp_limb_t *pdest = num;
    //     const struct mp_power *ttab = &_fpioconst_pow10[0];
    // 
    //     assert (dig_no > int_no
    // 	    && exponent <= 0
    // 	    && exponent >= MIN_10_EXP - (DIG + 2));
    // 
    //     /* We need to compute MANT_DIG - BITS fractional bits that lie
    //        within the mantissa of the result, the following bit for
    //        rounding, and to know whether any subsequent bit is 0.
    //        Computing a bit with value 2^-n means looking at n digits after
    //        the decimal point.  */
    //     if (bits > 0)
    //       {
    // 	/* The bits required are those immediately after the point.  */
    // 	assert (int_no > 0 && exponent == 0);
    // 	need_frac_digits = 1 + MANT_DIG - bits;
    //       }
    //     else
    //       {
    // 	/* The number is in the form .123eEXPONENT.  */
    // 	assert (int_no == 0 && *startp != L_('0'));
    // 	/* The number is at least 10^(EXPONENT-1), and 10^3 <
    // 	   2^10.  */
    // 	int neg_exp_2 = ((1 - exponent) * 10) / 3 + 1;
    // 	/* The number is at least 2^-NEG_EXP_2.  We need up to
    // 	   MANT_DIG bits following that bit.  */
    // 	need_frac_digits = neg_exp_2 + MANT_DIG;
    // 	/* However, we never need bits beyond 1/4 ulp of the smallest
    // 	   representable value.  (That 1/4 ulp bit is only needed to
    // 	   determine tinyness on machines where tinyness is determined
    // 	   after rounding.)  */
    // 	if (need_frac_digits > MANT_DIG - MIN_EXP + 2)
    // 	  need_frac_digits = MANT_DIG - MIN_EXP + 2;
    // 	/* At this point, NEED_FRAC_DIGITS is the total number of
    // 	   digits needed after the point, but some of those may be
    // 	   leading 0s.  */
    // 	need_frac_digits += exponent;
    // 	/* Any cases underflowing enough that none of the fractional
    // 	   digits are needed should have been caught earlier (such
    // 	   cases are on the order of 10^-n or smaller where 2^-n is
    // 	   the least subnormal).  */
    // 	assert (need_frac_digits > 0);
    //       }
    // 
    //     if (need_frac_digits > (intmax_t) dig_no - (intmax_t) int_no)
    //       need_frac_digits = (intmax_t) dig_no - (intmax_t) int_no;
    // 
    //     if ((intmax_t) dig_no > (intmax_t) int_no + need_frac_digits)
    //       {
    // 	dig_no = int_no + need_frac_digits;
    // 	more_bits = 1;
    //       }
    //     else
    //       more_bits = 0;
    // 
    //     neg_exp = (intmax_t) dig_no - (intmax_t) int_no - exponent;
    // 
    //     /* Construct the denominator.  */
    //     densize = 0;
    //     expbit = 1;
    //     do
    //       {
    // 	if ((neg_exp & expbit) != 0)
    // 	  {
    // 	    mp_limb_t cy;
    // 	    neg_exp ^= expbit;
    // 
    // 	    if (densize == 0)
    // 	      {
    // 		densize = ttab->arraysize - _FPIO_CONST_OFFSET;
    // 		memcpy (psrc, &__tens[ttab->arrayoff + _FPIO_CONST_OFFSET],
    // 			densize * sizeof (mp_limb_t));
    // 	      }
    // 	    else
    // 	      {
    // 		cy = __mpn_mul (pdest, &__tens[ttab->arrayoff
    // 					      + _FPIO_CONST_OFFSET],
    // 				ttab->arraysize - _FPIO_CONST_OFFSET,
    // 				psrc, densize);
    // 		densize += ttab->arraysize - _FPIO_CONST_OFFSET;
    // 		if (cy == 0)
    // 		  --densize;
    // 		(void) SWAP (psrc, pdest);
    // 	      }
    // 	  }
    // 	expbit <<= 1;
    // 	++ttab;
    //       }
    //     while (neg_exp != 0);
    // 
    //     if (psrc == num)
    //       memcpy (den, num, densize * sizeof (mp_limb_t));
    // 
    //     /* Read the fractional digits from the string.  */
    //     (void) str_to_mpn (startp, dig_no - int_no, num, &numsize, &exponent
    // #ifndef USE_WIDE_CHAR
    // 		       , decimal, decimal_len, thousands
    // #endif
    // 		       );
    // 
    //     /* We now have to shift both numbers so that the highest bit in the
    //        denominator is set.  In the same process we copy the numerator to
    //        a high place in the array so that the division constructs the wanted
    //        digits.  This is done by a "quasi fix point" number representation.
    // 
    //        num:   ddddddddddd . 0000000000000000000000
    // 	      |--- m ---|
    //        den:                            ddddddddddd      n >= m
    // 				       |--- n ---|
    //      */
    // 
    //     cnt = stdc_leading_zeros (den[densize - 1]);
    // 
    // 
    //     if (cnt > 0)
    //       {
    // 	/* Don't call `mpn_shift' with a count of zero since the specification
    // 	   does not allow this.  */
    // 	(void) __mpn_lshift (den, den, densize, cnt);
    // 	cy = __mpn_lshift (num, num, numsize, cnt);
    // 	if (cy != 0)
    // 	  num[numsize++] = cy;
    //       }
    // 
    //     /* Now we are ready for the division.  But it is not necessary to
    //        do a full multi-precision division because we only need a small
    //        number of bits for the result.  So we do not use __mpn_divmod
    //        here but instead do the division here by hand and stop whenever
    //        the needed number of bits is reached.  The code itself comes
    //        from the GNU MP Library by Torbj\"orn Granlund.  */
    // 
    //     exponent = bits;
    // 
    //     switch (densize)
    //       {
    //       case 1:
    // 	{
    // 	  mp_limb_t d, n, quot;
    // 	  int used = 0;
    // 
    // 	  n = num[0];
    // 	  d = den[0];
    // 	  assert (numsize == 1 && n < d);
    // 
    // 	  do
    // 	    {
    // 	      udiv_qrnnd (quot, n, n, 0, d);
    // 
    // #define got_limb							      \
    // 	      if (bits == 0)						      \
    // 		{							      \
    // 		  int cnt = stdc_leading_zeros (quot);			      \
    // 		  exponent -= cnt;					      \
    // 		  if (BITS_PER_MP_LIMB - cnt > MANT_DIG)		      \
    // 		    {							      \
    // 		      used = MANT_DIG + cnt;				      \
    // 		      retval[0] = quot >> (BITS_PER_MP_LIMB - used);	      \
    // 		      bits = MANT_DIG + 1;				      \
    // 		    }							      \
    // 		  else							      \
    // 		    {							      \
    // 		      /* Note that we only clear the second element.  */      \
    // 		      /* The conditional is determined at compile time.  */   \
    // 		      if (RETURN_LIMB_SIZE > 1)				      \
    // 			retval[1] = 0;					      \
    // 		      retval[0] = quot;					      \
    // 		      bits = -cnt;					      \
    // 		    }							      \
    // 		}							      \
    // 	      else if (bits + BITS_PER_MP_LIMB <= MANT_DIG)		      \
    // 		__mpn_lshift_1 (retval, RETURN_LIMB_SIZE, BITS_PER_MP_LIMB,   \
    // 				quot);					      \
    // 	      else							      \
    // 		{							      \
    // 		  used = MANT_DIG - bits;				      \
    // 		  if (used > 0)						      \
    // 		    __mpn_lshift_1 (retval, RETURN_LIMB_SIZE, used, quot);    \
    // 		}							      \
    // 	      bits += BITS_PER_MP_LIMB
    // 
    // 	      got_limb;
    // 	    }
    // 	  while (bits <= MANT_DIG);
    // 
    // 	  return round_and_return (retval, exponent - 1, negative,
    // 				   quot, BITS_PER_MP_LIMB - 1 - used,
    // 				   more_bits || n != 0);
    // 	}
    //       case 2:
    // 	{
    // 	  mp_limb_t d0, d1, n0, n1;
    // 	  mp_limb_t quot = 0;
    // 	  int used = 0;
    // 
    // 	  d0 = den[0];
    // 	  d1 = den[1];
    // 
    // 	  if (numsize < densize)
    // 	    {
    // 	      if (num[0] >= d1)
    // 		{
    // 		  /* The numerator of the number occupies fewer bits than
    // 		     the denominator but the one limb is bigger than the
    // 		     high limb of the numerator.  */
    // 		  n1 = 0;
    // 		  n0 = num[0];
    // 		}
    // 	      else
    // 		{
    // 		  if (bits <= 0)
    // 		    exponent -= BITS_PER_MP_LIMB;
    // 		  else
    // 		    {
    // 		      if (bits + BITS_PER_MP_LIMB <= MANT_DIG)
    // 			__mpn_lshift_1 (retval, RETURN_LIMB_SIZE,
    // 					BITS_PER_MP_LIMB, 0);
    // 		      else
    // 			{
    // 			  used = MANT_DIG - bits;
    // 			  if (used > 0)
    // 			    __mpn_lshift_1 (retval, RETURN_LIMB_SIZE, used, 0);
    // 			}
    // 		      bits += BITS_PER_MP_LIMB;
    // 		    }
    // 		  n1 = num[0];
    // 		  n0 = 0;
    // 		}
    // 	    }
    // 	  else
    // 	    {
    // 	      n1 = num[1];
    // 	      n0 = num[0];
    // 	    }
    // 
    // 	  while (bits <= MANT_DIG)
    // 	    {
    // 	      mp_limb_t r;
    // 
    // 	      if (n1 == d1)
    // 		{
    // 		  /* QUOT should be either 111..111 or 111..110.  We need
    // 		     special treatment of this rare case as normal division
    // 		     would give overflow.  */
    // 		  quot = ~(mp_limb_t) 0;
    // 
    // 		  r = n0 + d1;
    // 		  if (r < d1)	/* Carry in the addition?  */
    // 		    {
    // 		      add_ssaaaa (n1, n0, r - d0, 0, 0, d0);
    // 		      goto have_quot;
    // 		    }
    // 		  n1 = d0 - (d0 != 0);
    // 		  n0 = -d0;
    // 		}
    // 	      else
    // 		{
    // 		  udiv_qrnnd (quot, r, n1, n0, d1);
    // 		  umul_ppmm (n1, n0, d0, quot);
    // 		}
    // 
    // 	    q_test:
    // 	      if (n1 > r || (n1 == r && n0 > 0))
    // 		{
    // 		  /* The estimated QUOT was too large.  */
    // 		  --quot;
    // 
    // 		  sub_ddmmss (n1, n0, n1, n0, 0, d0);
    // 		  r += d1;
    // 		  if (r >= d1)	/* If not carry, test QUOT again.  */
    // 		    goto q_test;
    // 		}
    // 	      sub_ddmmss (n1, n0, r, 0, n1, n0);
    // 
    // 	    have_quot:
    // 	      got_limb;
    // 	    }
    // 
    // 	  return round_and_return (retval, exponent - 1, negative,
    // 				   quot, BITS_PER_MP_LIMB - 1 - used,
    // 				   more_bits || n1 != 0 || n0 != 0);
    // 	}
    //       default:
    // 	{
    // 	  int i;
    // 	  mp_limb_t cy, dX, d1, n0, n1;
    // 	  mp_limb_t quot = 0;
    // 	  int used = 0;
    // 
    // 	  dX = den[densize - 1];
    // 	  d1 = den[densize - 2];
    // 
    // 	  /* The division does not work if the upper limb of the two-limb
    // 	     numerator is greater than or equal to the denominator.  */
    // 	  if (__mpn_cmp (num, &den[densize - numsize], numsize) >= 0)
    // 	    num[numsize++] = 0;
    // 
    // 	  if (numsize < densize)
    // 	    {
    // 	      mp_size_t empty = densize - numsize;
    // 	      int i;
    // 
    // 	      if (bits <= 0)
    // 		exponent -= empty * BITS_PER_MP_LIMB;
    // 	      else
    // 		{
    // 		  if (bits + empty * BITS_PER_MP_LIMB <= MANT_DIG)
    // 		    {
    // 		      /* We make a difference here because the compiler
    // 			 cannot optimize the `else' case that good and
    // 			 this reflects all currently used FLOAT types
    // 			 and GMP implementations.  */
    // #if RETURN_LIMB_SIZE <= 2
    // 		      assert (empty == 1);
    // 		      __mpn_lshift_1 (retval, RETURN_LIMB_SIZE,
    // 				      BITS_PER_MP_LIMB, 0);
    // #else
    // 		      for (i = RETURN_LIMB_SIZE - 1; i >= empty; --i)
    // 			retval[i] = retval[i - empty];
    // 		      while (i >= 0)
    // 			retval[i--] = 0;
    // #endif
    // 		    }
    // 		  else
    // 		    {
    // 		      used = MANT_DIG - bits;
    // 		      if (used >= BITS_PER_MP_LIMB)
    // 			{
    // 			  int i;
    // 			  (void) __mpn_lshift (&retval[used
    // 						       / BITS_PER_MP_LIMB],
    // 					       retval,
    // 					       (RETURN_LIMB_SIZE
    // 						- used / BITS_PER_MP_LIMB),
    // 					       used % BITS_PER_MP_LIMB);
    // 			  for (i = used / BITS_PER_MP_LIMB - 1; i >= 0; --i)
    // 			    retval[i] = 0;
    // 			}
    // 		      else if (used > 0)
    // 			__mpn_lshift_1 (retval, RETURN_LIMB_SIZE, used, 0);
    // 		    }
    // 		  bits += empty * BITS_PER_MP_LIMB;
    // 		}
    // 	      for (i = numsize; i > 0; --i)
    // 		num[i + empty] = num[i - 1];
    // 	      MPN_ZERO (num, empty + 1);
    // 	    }
    // 	  else
    // 	    {
    // 	      int i;
    // 	      assert (numsize == densize);
    // 	      for (i = numsize; i > 0; --i)
    // 		num[i] = num[i - 1];
    // 	      num[0] = 0;
    // 	    }
    // 
    // 	  den[densize] = 0;
    // 	  n0 = num[densize];
    // 
    // 	  while (bits <= MANT_DIG)
    // 	    {
    // 	      if (n0 == dX)
    // 		/* This might over-estimate QUOT, but it's probably not
    // 		   worth the extra code here to find out.  */
    // 		quot = ~(mp_limb_t) 0;
    // 	      else
    // 		{
    // 		  mp_limb_t r;
    // 
    // 		  udiv_qrnnd (quot, r, n0, num[densize - 1], dX);
    // 		  umul_ppmm (n1, n0, d1, quot);
    // 
    // 		  while (n1 > r || (n1 == r && n0 > num[densize - 2]))
    // 		    {
    // 		      --quot;
    // 		      r += dX;
    // 		      if (r < dX) /* I.e. "carry in previous addition?" */
    // 			break;
    // 		      n1 -= n0 < d1;
    // 		      n0 -= d1;
    // 		    }
    // 		}
    // 
    // 	      /* Possible optimization: We already have (q * n0) and (1 * n1)
    // 		 after the calculation of QUOT.  Taking advantage of this, we
    // 		 could make this loop make two iterations less.  */
    // 
    // 	      cy = __mpn_submul_1 (num, den, densize + 1, quot);
    // 
    // 	      if (num[densize] != cy)
    // 		{
    // 		  cy = __mpn_add_n (num, num, den, densize);
    // 		  assert (cy != 0);
    // 		  --quot;
    // 		}
    // 	      n0 = num[densize] = num[densize - 1];
    // 	      for (i = densize - 1; i > 0; --i)
    // 		num[i] = num[i - 1];
    // 	      num[0] = 0;
    // 
    // 	      got_limb;
    // 	    }
    // 
    // 	  for (i = densize; i >= 0 && num[i] == 0; --i)
    // 	    ;
    // 	  return round_and_return (retval, exponent - 1, negative,
    // 				   quot, BITS_PER_MP_LIMB - 1 - used,
    // 				   more_bits || i >= 0);
    // 	}
    //       }
    //   }
    // 
    //   /* NOTREACHED */
    // }
    // #if defined _LIBC && !defined USE_WIDE_CHAR
    // libc_hidden_def (____STRTOF_INTERNAL)
    // #endif
    // 
    // 
    // /* External user entry point.  */
    // 
    // FLOAT
    // #ifdef weak_function
    // weak_function
    // #endif
    // __STRTOF (const STRING_TYPE *nptr, STRING_TYPE **endptr, locale_t loc)
    // {
    //   return ____STRTOF_INTERNAL (nptr, endptr, 0, loc);
    // }
    // #if defined _LIBC
    // libc_hidden_def (__STRTOF)
    // libc_hidden_ver (__STRTOF, STRTOF)
    // #endif
    // static_weak_alias (__STRTOF, STRTOF)
    // 
    // #ifdef LONG_DOUBLE_COMPAT
    // # if LONG_DOUBLE_COMPAT(libc, GLIBC_2_1)
    // #  ifdef USE_WIDE_CHAR
    // compat_symbol (libc, __wcstod_l, __wcstold_l, GLIBC_2_1);
    // #  else
    // compat_symbol (libc, __strtod_l, __strtold_l, GLIBC_2_1);
    // #  endif
    // # endif
    // # if LONG_DOUBLE_COMPAT(libc, GLIBC_2_3)
    // #  ifdef USE_WIDE_CHAR
    // compat_symbol (libc, wcstod_l, wcstold_l, GLIBC_2_3);
    // #  else
    // compat_symbol (libc, strtod_l, strtold_l, GLIBC_2_3);
    // #  endif
    // # endif
    // #endif
    // 
    // #if BUILD_DOUBLE
    // # if __HAVE_FLOAT64 && !__HAVE_DISTINCT_FLOAT64
    // #  undef strtof64_l
    // #  undef wcstof64_l
    // #  ifdef USE_WIDE_CHAR
    // weak_alias (__wcstod_l, wcstof64_l)
    // #  else
    // weak_alias (__strtod_l, strtof64_l)
    // #  endif
    // # endif
    // # if __HAVE_FLOAT32X && !__HAVE_DISTINCT_FLOAT32X
    // #  undef strtof32x_l
    // #  undef wcstof32x_l
    // #  ifdef USE_WIDE_CHAR
    // weak_alias (__wcstod_l, wcstof32x_l)
    // #  else
    // weak_alias (__strtod_l, strtof32x_l)
    // #  endif
    // # endif
    // #endif
    // glibc❗❌: exact rational comparisons implement the complete modeled
    // decimal result, including directed rounding, without floating arithmetic.
    // The fixed binary64 search is an explicit constant-factor cost difference
    // from GNU's limb extraction; no source algorithm-equivalence claim.
    // Any binary64 endpoint has <=767 significant decimal digits; adjacent
    // midpoints have <=768. The complete significant suffix is scanned. A
    // nonzero suffix beyond 768 is represented by one sticky digit at 769:
    // both numbers are strictly on the same side of every endpoint/midpoint.
    // Historical loaded GNU sources/global locale custody remains unresolved.
    let magnitude = if exponent >= 309 {
        if round_away(negative, true, true, true, mode) {
            INFINITY_BITS
        } else {
            MAX_FINITE_BITS
        }
    } else if exponent < -324 {
        if round_away(negative, false, false, true, mode) {
            1
        } else {
            0
        }
    } else {
        let exact = ExactDecimal::new(digits, exponent);
        let mut low = 0_u64;
        let mut high = MAX_FINITE_BITS;
        while low < high {
            let middle = low + (high - low + 1) / 2;
            let (coefficient, shift) = binary64_parts(middle);
            if exact.compare_binary(coefficient, shift) == Ordering::Less {
                high = middle - 1;
            } else {
                low = middle;
            }
        }
        let (coefficient, shift) = binary64_parts(low);
        if exact.compare_binary(coefficient, shift) == Ordering::Equal {
            low
        } else {
            let (next_coefficient, next_shift) = if low == MAX_FINITE_BITS {
                (1_u64 << 53, 971)
            } else {
                binary64_parts(low + 1)
            };
            let midpoint_coefficient = coefficient + (next_coefficient << (next_shift - shift));
            let midpoint = exact.compare_binary(midpoint_coefficient, shift - 1);
            let half_bit = midpoint != Ordering::Less;
            let more_bits = midpoint != Ordering::Equal;
            if round_away(negative, low & 1 != 0, half_bit, more_bits, mode) {
                if low == MAX_FINITE_BITS {
                    INFINITY_BITS
                } else {
                    low + 1
                }
            } else {
                low
            }
        }
    };
    magnitude | if negative { SIGN_BIT } else { 0 }
}

fn ascii_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
}

/// C-locale formatted double extraction. File grammar and warning assembly
/// stay with their owning IO parser. This extracts decimal, not C atof/hex.
#[rustfmt::skip]
pub fn source_extract_double(
    bytes: &[u8],
    skip_whitespace: bool,
    mode: SourceNumericRoundingMode,
) -> SourceDoubleExtraction {
    // // Locale support -*- C++ -*-
    // 
    // // Copyright (C) 1997-2024 Free Software Foundation, Inc.
    // //
    // // This file is part of the GNU ISO C++ Library.  This library is free
    // // software; you can redistribute it and/or modify it under the
    // // terms of the GNU General Public License as published by the
    // // Free Software Foundation; either version 3, or (at your option)
    // // any later version.
    // 
    // // This library is distributed in the hope that it will be useful,
    // // but WITHOUT ANY WARRANTY; without even the implied warranty of
    // // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    // // GNU General Public License for more details.
    // 
    // // Under Section 7 of GPL version 3, you are granted additional
    // // permissions described in the GCC Runtime Library Exception, version
    // // 3.1, as published by the Free Software Foundation.
    // 
    // // You should have received a copy of the GNU General Public License and
    // // a copy of the GCC Runtime Library Exception along with this program;
    // // see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
    // // <http://www.gnu.org/licenses/>.
    // 
    // /** @file bits/locale_facets.tcc
    //  *  This is an internal header file, included by other library headers.
    //  *  Do not attempt to use it directly. @headername{locale}
    //  */
    // 
    // #ifndef _LOCALE_FACETS_TCC
    // #define _LOCALE_FACETS_TCC 1
    // 
    // #pragma GCC system_header
    // 
    // namespace std _GLIBCXX_VISIBILITY(default)
    // {
    // _GLIBCXX_BEGIN_NAMESPACE_VERSION
    // 
    //   // Routine to access a cache for the facet.  If the cache didn't
    //   // exist before, it gets constructed on the fly.
    //   template<typename _Facet>
    //     struct __use_cache
    //     {
    //       const _Facet*
    //       operator() (const locale& __loc) const;
    //     };
    // 
    //   // Specializations.
    //   template<typename _CharT>
    //     struct __use_cache<__numpunct_cache<_CharT> >
    //     {
    //       const __numpunct_cache<_CharT>*
    //       operator() (const locale& __loc) const
    //       {
    // 	const size_t __i = numpunct<_CharT>::id._M_id();
    // 	const locale::facet** __caches = __loc._M_impl->_M_caches;
    // 	if (!__caches[__i])
    // 	  {
    // 	    __numpunct_cache<_CharT>* __tmp = 0;
    // 	    __try
    // 	      {
    // 		__tmp = new __numpunct_cache<_CharT>;
    // 		__tmp->_M_cache(__loc);
    // 	      }
    // 	    __catch(...)
    // 	      {
    // 		delete __tmp;
    // 		__throw_exception_again;
    // 	      }
    // 	    __loc._M_impl->_M_install_cache(__tmp, __i);
    // 	  }
    // 	return static_cast<const __numpunct_cache<_CharT>*>(__caches[__i]);
    //       }
    //     };
    // 
    //   template<typename _CharT>
    //     void
    //     __numpunct_cache<_CharT>::_M_cache(const locale& __loc)
    //     {
    //       const numpunct<_CharT>& __np = use_facet<numpunct<_CharT> >(__loc);
    // 
    //       char* __grouping = 0;
    //       _CharT* __truename = 0;
    //       _CharT* __falsename = 0;
    //       __try
    // 	{
    // 	  const string& __g = __np.grouping();
    // 	  _M_grouping_size = __g.size();
    // 	  __grouping = new char[_M_grouping_size];
    // 	  __g.copy(__grouping, _M_grouping_size);
    // 	  _M_use_grouping = (_M_grouping_size
    // 			     && static_cast<signed char>(__grouping[0]) > 0
    // 			     && (__grouping[0]
    // 				 != __gnu_cxx::__numeric_traits<char>::__max));
    // 
    // 	  const basic_string<_CharT>& __tn = __np.truename();
    // 	  _M_truename_size = __tn.size();
    // 	  __truename = new _CharT[_M_truename_size];
    // 	  __tn.copy(__truename, _M_truename_size);
    // 
    // 	  const basic_string<_CharT>& __fn = __np.falsename();
    // 	  _M_falsename_size = __fn.size();
    // 	  __falsename = new _CharT[_M_falsename_size];
    // 	  __fn.copy(__falsename, _M_falsename_size);
    // 
    // 	  _M_decimal_point = __np.decimal_point();
    // 	  _M_thousands_sep = __np.thousands_sep();
    // 
    // 	  const ctype<_CharT>& __ct = use_facet<ctype<_CharT> >(__loc);
    // 	  __ct.widen(__num_base::_S_atoms_out,
    // 		     __num_base::_S_atoms_out
    // 		     + __num_base::_S_oend, _M_atoms_out);
    // 	  __ct.widen(__num_base::_S_atoms_in,
    // 		     __num_base::_S_atoms_in
    // 		     + __num_base::_S_iend, _M_atoms_in);
    // 
    // 	  _M_grouping = __grouping;
    // 	  _M_truename = __truename;
    // 	  _M_falsename = __falsename;
    // 	  _M_allocated = true;
    // 	}
    //       __catch(...)
    // 	{
    // 	  delete [] __grouping;
    // 	  delete [] __truename;
    // 	  delete [] __falsename;
    // 	  __throw_exception_again;
    // 	}
    //     }
    // 
    //   // Used by both numeric and monetary facets.
    //   // Check to make sure that the __grouping_tmp string constructed in
    //   // money_get or num_get matches the canonical grouping for a given
    //   // locale.
    //   // __grouping_tmp is parsed L to R
    //   // 1,222,444 == __grouping_tmp of "\1\3\3"
    //   // __grouping is parsed R to L
    //   // 1,222,444 == __grouping of "\3" == "\3\3\3"
    //   _GLIBCXX_PURE bool
    //   __verify_grouping(const char* __grouping, size_t __grouping_size,
    // 		    const string& __grouping_tmp) throw ();
    // 
    // _GLIBCXX_BEGIN_NAMESPACE_LDBL
    // 
    //   template<typename _CharT, typename _InIter>
    //     _GLIBCXX_DEFAULT_ABI_TAG
    //     _InIter
    //     num_get<_CharT, _InIter>::
    //     _M_extract_float(_InIter __beg, _InIter __end, ios_base& __io,
    // 		     ios_base::iostate& __err, string& __xtrc) const
    //     {
    //       typedef char_traits<_CharT>			__traits_type;
    //       typedef __numpunct_cache<_CharT>                  __cache_type;
    //       __use_cache<__cache_type> __uc;
    //       const locale& __loc = __io._M_getloc();
    //       const __cache_type* __lc = __uc(__loc);
    //       const _CharT* __lit = __lc->_M_atoms_in;
    //       char_type __c = char_type();
    // 
    //       // True if __beg becomes equal to __end.
    //       bool __testeof = __beg == __end;
    // 
    //       // First check for sign.
    //       if (!__testeof)
    // 	{
    // 	  __c = *__beg;
    // 	  const bool __plus = __c == __lit[__num_base::_S_iplus];
    // 	  if ((__plus || __c == __lit[__num_base::_S_iminus])
    // 	      && !(__lc->_M_use_grouping && __c == __lc->_M_thousands_sep)
    // 	      && !(__c == __lc->_M_decimal_point))
    // 	    {
    // 	      __xtrc += __plus ? '+' : '-';
    // 	      if (++__beg != __end)
    // 		__c = *__beg;
    // 	      else
    // 		__testeof = true;
    // 	    }
    // 	}
    // 
    //       // Next, look for leading zeros.
    //       bool __found_mantissa = false;
    //       int __sep_pos = 0;
    //       while (!__testeof)
    // 	{
    // 	  if ((__lc->_M_use_grouping && __c == __lc->_M_thousands_sep)
    // 	      || __c == __lc->_M_decimal_point)
    // 	    break;
    // 	  else if (__c == __lit[__num_base::_S_izero])
    // 	    {
    // 	      if (!__found_mantissa)
    // 		{
    // 		  __xtrc += '0';
    // 		  __found_mantissa = true;
    // 		}
    // 	      ++__sep_pos;
    // 
    // 	      if (++__beg != __end)
    // 		__c = *__beg;
    // 	      else
    // 		__testeof = true;
    // 	    }
    // 	  else
    // 	    break;
    // 	}
    // 
    //       // Only need acceptable digits for floating point numbers.
    //       bool __found_dec = false;
    //       bool __found_sci = false;
    //       string __found_grouping;
    //       if (__lc->_M_use_grouping)
    // 	__found_grouping.reserve(32);
    //       const char_type* __lit_zero = __lit + __num_base::_S_izero;
    // 
    //       if (!__lc->_M_allocated)
    // 	// "C" locale
    // 	while (!__testeof)
    // 	  {
    // 	    const int __digit = _M_find(__lit_zero, 10, __c);
    // 	    if (__digit != -1)
    // 	      {
    // 		__xtrc += '0' + __digit;
    // 		__found_mantissa = true;
    // 	      }
    // 	    else if (__c == __lc->_M_decimal_point
    // 		     && !__found_dec && !__found_sci)
    // 	      {
    // 		__xtrc += '.';
    // 		__found_dec = true;
    // 	      }
    // 	    else if ((__c == __lit[__num_base::_S_ie] 
    // 		      || __c == __lit[__num_base::_S_iE])
    // 		     && !__found_sci && __found_mantissa)
    // 	      {
    // 		// Scientific notation.
    // 		__xtrc += 'e';
    // 		__found_sci = true;
    // 		
    // 		// Remove optional plus or minus sign, if they exist.
    // 		if (++__beg != __end)
    // 		  {
    // 		    __c = *__beg;
    // 		    const bool __plus = __c == __lit[__num_base::_S_iplus];
    // 		    if (__plus || __c == __lit[__num_base::_S_iminus])
    // 		      __xtrc += __plus ? '+' : '-';
    // 		    else
    // 		      continue;
    // 		  }
    // 		else
    // 		  {
    // 		    __testeof = true;
    // 		    break;
    // 		  }
    // 	      }
    // 	    else
    // 	      break;
    // 
    // 	    if (++__beg != __end)
    // 	      __c = *__beg;
    // 	    else
    // 	      __testeof = true;
    // 	  }
    //       else
    // 	while (!__testeof)
    // 	  {
    // 	    // According to 22.2.2.1.2, p8-9, first look for thousands_sep
    // 	    // and decimal_point.
    // 	    if (__lc->_M_use_grouping && __c == __lc->_M_thousands_sep)
    // 	      {
    // 		if (!__found_dec && !__found_sci)
    // 		  {
    // 		    // NB: Thousands separator at the beginning of a string
    // 		    // is a no-no, as is two consecutive thousands separators.
    // 		    if (__sep_pos)
    // 		      {
    // 			__found_grouping += static_cast<char>(__sep_pos);
    // 			__sep_pos = 0;
    // 		      }
    // 		    else
    // 		      {
    // 			// NB: __convert_to_v will not assign __v and will
    // 			// set the failbit.
    // 			__xtrc.clear();
    // 			break;
    // 		      }
    // 		  }
    // 		else
    // 		  break;
    // 	      }
    // 	    else if (__c == __lc->_M_decimal_point)
    // 	      {
    // 		if (!__found_dec && !__found_sci)
    // 		  {
    // 		    // If no grouping chars are seen, no grouping check
    // 		    // is applied. Therefore __found_grouping is adjusted
    // 		    // only if decimal_point comes after some thousands_sep.
    // 		    if (__found_grouping.size())
    // 		      __found_grouping += static_cast<char>(__sep_pos);
    // 		    __xtrc += '.';
    // 		    __found_dec = true;
    // 		  }
    // 		else
    // 		  break;
    // 	      }
    // 	    else
    // 	      {
    // 		const char_type* __q =
    // 		  __traits_type::find(__lit_zero, 10, __c);
    // 		if (__q)
    // 		  {
    // 		    __xtrc += '0' + (__q - __lit_zero);
    // 		    __found_mantissa = true;
    // 		    ++__sep_pos;
    // 		  }
    // 		else if ((__c == __lit[__num_base::_S_ie] 
    // 			  || __c == __lit[__num_base::_S_iE])
    // 			 && !__found_sci && __found_mantissa)
    // 		  {
    // 		    // Scientific notation.
    // 		    if (__found_grouping.size() && !__found_dec)
    // 		      __found_grouping += static_cast<char>(__sep_pos);
    // 		    __xtrc += 'e';
    // 		    __found_sci = true;
    // 		    
    // 		    // Remove optional plus or minus sign, if they exist.
    // 		    if (++__beg != __end)
    // 		      {
    // 			__c = *__beg;
    // 			const bool __plus = __c == __lit[__num_base::_S_iplus];
    // 			if ((__plus || __c == __lit[__num_base::_S_iminus])
    // 			    && !(__lc->_M_use_grouping
    // 				 && __c == __lc->_M_thousands_sep)
    // 			    && !(__c == __lc->_M_decimal_point))
    // 		      __xtrc += __plus ? '+' : '-';
    // 			else
    // 			  continue;
    // 		      }
    // 		    else
    // 		      {
    // 			__testeof = true;
    // 			break;
    // 		      }
    // 		  }
    // 		else
    // 		  break;
    // 	      }
    // 	    
    // 	    if (++__beg != __end)
    // 	      __c = *__beg;
    // 	    else
    // 	      __testeof = true;
    // 	  }
    // 
    //       // Digit grouping is checked. If grouping and found_grouping don't
    //       // match, then get very very upset, and set failbit.
    //       if (__found_grouping.size())
    //         {
    //           // Add the ending grouping if a decimal or 'e'/'E' wasn't found.
    // 	  if (!__found_dec && !__found_sci)
    // 	    __found_grouping += static_cast<char>(__sep_pos);
    // 
    //           if (!std::__verify_grouping(__lc->_M_grouping, 
    // 				      __lc->_M_grouping_size,
    // 				      __found_grouping))
    // 	    __err = ios_base::failbit;
    //         }
    // 
    //       return __beg;
    //     }
    // 
    //   template<typename _CharT, typename _InIter>
    //     template<typename _ValueT>
    //       _GLIBCXX_DEFAULT_ABI_TAG
    //       _InIter
    //       num_get<_CharT, _InIter>::
    //       _M_extract_int(_InIter __beg, _InIter __end, ios_base& __io,
    // 		     ios_base::iostate& __err, _ValueT& __v) const
    //       {
    //         typedef char_traits<_CharT>			    __traits_type;
    // 	using __gnu_cxx::__add_unsigned;
    // 	typedef typename __add_unsigned<_ValueT>::__type    __unsigned_type;
    // 	typedef __numpunct_cache<_CharT>                    __cache_type;
    // 	__use_cache<__cache_type> __uc;
    // 	const locale& __loc = __io._M_getloc();
    // 	const __cache_type* __lc = __uc(__loc);
    // 	const _CharT* __lit = __lc->_M_atoms_in;
    // 	char_type __c = char_type();
    // 
    // 	// NB: Iff __basefield == 0, __base can change based on contents.
    // 	const ios_base::fmtflags __basefield = __io.flags()
    // 	                                       & ios_base::basefield;
    // 	const bool __oct = __basefield == ios_base::oct;
    // 	int __base = __oct ? 8 : (__basefield == ios_base::hex ? 16 : 10);
    // 
    // 	// True if __beg becomes equal to __end.
    // 	bool __testeof = __beg == __end;
    // 
    // 	// First check for sign.
    // 	bool __negative = false;
    // 	if (!__testeof)
    // 	  {
    // 	    __c = *__beg;
    // 	    __negative = __c == __lit[__num_base::_S_iminus];
    // 	    if ((__negative || __c == __lit[__num_base::_S_iplus])
    // 		&& !(__lc->_M_use_grouping && __c == __lc->_M_thousands_sep)
    // 		&& !(__c == __lc->_M_decimal_point))
    // 	      {
    // 		if (++__beg != __end)
    // 		  __c = *__beg;
    // 		else
    // 		  __testeof = true;
    // 	      }
    // 	  }
    // 
    // 	// Next, look for leading zeros and check required digits
    // 	// for base formats.
    // 	bool __found_zero = false;
    // 	int __sep_pos = 0;
    // 	while (!__testeof)
    // 	  {
    // 	    if ((__lc->_M_use_grouping && __c == __lc->_M_thousands_sep)
    // 		|| __c == __lc->_M_decimal_point)
    // 	      break;
    // 	    else if (__c == __lit[__num_base::_S_izero] 
    // 		     && (!__found_zero || __base == 10))
    // 	      {
    // 		__found_zero = true;
    // 		++__sep_pos;
    // 		if (__basefield == 0)
    // 		  __base = 8;
    // 		if (__base == 8)
    // 		  __sep_pos = 0;
    // 	      }
    // 	    else if (__found_zero
    // 		     && (__c == __lit[__num_base::_S_ix]
    // 			 || __c == __lit[__num_base::_S_iX]))
    // 	      {
    // 		if (__basefield == 0)
    // 		  __base = 16;
    // 		if (__base == 16)
    // 		  {
    // 		    __found_zero = false;
    // 		    __sep_pos = 0;
    // 		  }
    // 		else
    // 		  break;
    // 	      }
    // 	    else
    // 	      break;
    // 
    // 	    if (++__beg != __end)
    // 	      {
    // 		__c = *__beg;
    // 		if (!__found_zero)
    // 		  break;
    // 	      }
    // 	    else
    // 	      __testeof = true;
    // 	  }
    // 	
    // 	// At this point, base is determined. If not hex, only allow
    // 	// base digits as valid input.
    // 	const size_t __len = (__base == 16 ? __num_base::_S_iend
    // 			      - __num_base::_S_izero : __base);
    // 
    // 	// Extract.
    // 	typedef __gnu_cxx::__numeric_traits<_ValueT> __num_traits;
    // 	string __found_grouping;
    // 	if (__lc->_M_use_grouping)
    // 	  __found_grouping.reserve(32);
    // 	bool __testfail = false;
    // 	bool __testoverflow = false;
    // 	const __unsigned_type __max =
    // 	  (__negative && __num_traits::__is_signed)
    // 	  ? -static_cast<__unsigned_type>(__num_traits::__min)
    // 	  : __num_traits::__max;
    // 	const __unsigned_type __smax = __max / __base;
    // 	__unsigned_type __result = 0;
    // 	int __digit = 0;
    // 	const char_type* __lit_zero = __lit + __num_base::_S_izero;
    // 
    // 	if (!__lc->_M_allocated)
    // 	  // "C" locale
    // 	  while (!__testeof)
    // 	    {
    // 	      __digit = _M_find(__lit_zero, __len, __c);
    // 	      if (__digit == -1)
    // 		break;
    // 	      
    // 	      if (__result > __smax)
    // 		__testoverflow = true;
    // 	      else
    // 		{
    // 		  __result *= __base;
    // 		  __testoverflow |= __result > __max - __digit;
    // 		  __result += __digit;
    // 		  ++__sep_pos;
    // 		}
    // 	      
    // 	      if (++__beg != __end)
    // 		__c = *__beg;
    // 	      else
    // 		__testeof = true;
    // 	    }
    // 	else
    // 	  while (!__testeof)
    // 	    {
    // 	      // According to 22.2.2.1.2, p8-9, first look for thousands_sep
    // 	      // and decimal_point.
    // 	      if (__lc->_M_use_grouping && __c == __lc->_M_thousands_sep)
    // 		{
    // 		  // NB: Thousands separator at the beginning of a string
    // 		  // is a no-no, as is two consecutive thousands separators.
    // 		  if (__sep_pos)
    // 		    {
    // 		      __found_grouping += static_cast<char>(__sep_pos);
    // 		      __sep_pos = 0;
    // 		    }
    // 		  else
    // 		    {
    // 		      __testfail = true;
    // 		      break;
    // 		    }
    // 		}
    // 	      else if (__c == __lc->_M_decimal_point)
    // 		break;
    // 	      else
    // 		{
    // 		  const char_type* __q =
    // 		    __traits_type::find(__lit_zero, __len, __c);
    // 		  if (!__q)
    // 		    break;
    // 		  
    // 		  __digit = __q - __lit_zero;
    // 		  if (__digit > 15)
    // 		    __digit -= 6;
    // 		  if (__result > __smax)
    // 		    __testoverflow = true;
    // 		  else
    // 		    {
    // 		      __result *= __base;
    // 		      __testoverflow |= __result > __max - __digit;
    // 		      __result += __digit;
    // 		      ++__sep_pos;
    // 		    }
    // 		}
    // 	      
    // 	      if (++__beg != __end)
    // 		__c = *__beg;
    // 	      else
    // 		__testeof = true;
    // 	    }
    // 	
    // 	// Digit grouping is checked. If grouping and found_grouping don't
    // 	// match, then get very very upset, and set failbit.
    // 	if (__found_grouping.size())
    // 	  {
    // 	    // Add the ending grouping.
    // 	    __found_grouping += static_cast<char>(__sep_pos);
    // 
    // 	    if (!std::__verify_grouping(__lc->_M_grouping,
    // 					__lc->_M_grouping_size,
    // 					__found_grouping))
    // 	      __err = ios_base::failbit;
    // 	  }
    // 
    // 	// _GLIBCXX_RESOLVE_LIB_DEFECTS
    // 	// 23. Num_get overflow result.
    // 	if ((!__sep_pos && !__found_zero && !__found_grouping.size())
    // 	    || __testfail)
    // 	  {
    // 	    __v = 0;
    // 	    __err = ios_base::failbit;
    // 	  }
    // 	else if (__testoverflow)
    // 	  {
    // 	    if (__negative && __num_traits::__is_signed)
    // 	      __v = __num_traits::__min;
    // 	    else
    // 	      __v = __num_traits::__max;
    // 	    __err = ios_base::failbit;
    // 	  }
    // 	else
    // 	  __v = __negative ? -__result : __result;
    // 
    // 	if (__testeof)
    // 	  __err |= ios_base::eofbit;
    // 	return __beg;
    //       }
    // 
    //   // _GLIBCXX_RESOLVE_LIB_DEFECTS
    //   // 17.  Bad bool parsing
    //   template<typename _CharT, typename _InIter>
    //     _InIter
    //     num_get<_CharT, _InIter>::
    //     do_get(iter_type __beg, iter_type __end, ios_base& __io,
    //            ios_base::iostate& __err, bool& __v) const
    //     {
    //       if (!(__io.flags() & ios_base::boolalpha))
    //         {
    // 	  // Parse bool values as long.
    //           // NB: We can't just call do_get(long) here, as it might
    //           // refer to a derived class.
    // 	  long __l = -1;
    //           __beg = _M_extract_int(__beg, __end, __io, __err, __l);
    // 	  if (__l == 0 || __l == 1)
    // 	    __v = bool(__l);
    // 	  else
    // 	    {
    // 	      // _GLIBCXX_RESOLVE_LIB_DEFECTS
    // 	      // 23. Num_get overflow result.
    // 	      __v = true;
    // 	      __err = ios_base::failbit;
    // 	      if (__beg == __end)
    // 		__err |= ios_base::eofbit;
    // 	    }
    //         }
    //       else
    //         {
    // 	  // Parse bool values as alphanumeric.
    // 	  typedef __numpunct_cache<_CharT>  __cache_type;
    // 	  __use_cache<__cache_type> __uc;
    // 	  const locale& __loc = __io._M_getloc();
    // 	  const __cache_type* __lc = __uc(__loc);
    // 
    // 	  bool __testf = true;
    // 	  bool __testt = true;
    // 	  bool __donef = __lc->_M_falsename_size == 0;
    // 	  bool __donet = __lc->_M_truename_size == 0;
    // 	  bool __testeof = false;
    // 	  size_t __n = 0;
    // 	  while (!__donef || !__donet)
    // 	    {
    // 	      if (__beg == __end)
    // 		{
    // 		  __testeof = true;
    // 		  break;
    // 		}
    // 
    // 	      const char_type __c = *__beg;
    // 
    // 	      if (!__donef)
    // 		__testf = __c == __lc->_M_falsename[__n];
    // 
    // 	      if (!__testf && __donet)
    // 		break;
    // 
    // 	      if (!__donet)
    // 		__testt = __c == __lc->_M_truename[__n];
    // 
    // 	      if (!__testt && __donef)
    // 		break;
    // 
    // 	      if (!__testt && !__testf)
    // 		break;
    // 
    // 	      ++__n;
    // 	      ++__beg;
    // 
    // 	      __donef = !__testf || __n >= __lc->_M_falsename_size;
    // 	      __donet = !__testt || __n >= __lc->_M_truename_size;
    // 	    }
    // 	  if (__testf && __n == __lc->_M_falsename_size && __n)
    // 	    {
    // 	      __v = false;
    // 	      if (__testt && __n == __lc->_M_truename_size)
    // 		__err = ios_base::failbit;
    // 	      else
    // 		__err = __testeof ? ios_base::eofbit : ios_base::goodbit;
    // 	    }
    // 	  else if (__testt && __n == __lc->_M_truename_size && __n)
    // 	    {
    // 	      __v = true;
    // 	      __err = __testeof ? ios_base::eofbit : ios_base::goodbit;
    // 	    }
    // 	  else
    // 	    {
    // 	      // _GLIBCXX_RESOLVE_LIB_DEFECTS
    // 	      // 23. Num_get overflow result.
    // 	      __v = false;
    // 	      __err = ios_base::failbit;
    // 	      if (__testeof)
    // 		__err |= ios_base::eofbit;
    // 	    }
    // 	}
    //       return __beg;
    //     }
    // 
    //   template<typename _CharT, typename _InIter>
    //     _InIter
    //     num_get<_CharT, _InIter>::
    //     do_get(iter_type __beg, iter_type __end, ios_base& __io,
    // 	   ios_base::iostate& __err, float& __v) const
    //     {
    //       string __xtrc;
    //       __xtrc.reserve(32);
    //       __beg = _M_extract_float(__beg, __end, __io, __err, __xtrc);
    //       std::__convert_to_v(__xtrc.c_str(), __v, __err, _S_get_c_locale());
    //       if (__beg == __end)
    // 	__err |= ios_base::eofbit;
    //       return __beg;
    //     }
    // 
    //   template<typename _CharT, typename _InIter>
    //     _InIter
    //     num_get<_CharT, _InIter>::
    //     do_get(iter_type __beg, iter_type __end, ios_base& __io,
    //            ios_base::iostate& __err, double& __v) const
    //     {
    //       string __xtrc;
    //       __xtrc.reserve(32);
    //       __beg = _M_extract_float(__beg, __end, __io, __err, __xtrc);
    //       std::__convert_to_v(__xtrc.c_str(), __v, __err, _S_get_c_locale());
    //       if (__beg == __end)
    // 	__err |= ios_base::eofbit;
    //       return __beg;
    //     }
    // 
    // #if defined _GLIBCXX_LONG_DOUBLE_COMPAT && defined __LONG_DOUBLE_128__
    //   template<typename _CharT, typename _InIter>
    //     _InIter
    //     num_get<_CharT, _InIter>::
    //     __do_get(iter_type __beg, iter_type __end, ios_base& __io,
    // 	     ios_base::iostate& __err, double& __v) const
    //     {
    //       string __xtrc;
    //       __xtrc.reserve(32);
    //       __beg = _M_extract_float(__beg, __end, __io, __err, __xtrc);
    //       std::__convert_to_v(__xtrc.c_str(), __v, __err, _S_get_c_locale());
    //       if (__beg == __end)
    // 	__err |= ios_base::eofbit;
    //       return __beg;
    //     }
    // #endif
    // 
    //   template<typename _CharT, typename _InIter>
    //     _InIter
    //     num_get<_CharT, _InIter>::
    //     do_get(iter_type __beg, iter_type __end, ios_base& __io,
    //            ios_base::iostate& __err, long double& __v) const
    //     {
    //       string __xtrc;
    //       __xtrc.reserve(32);
    //       __beg = _M_extract_float(__beg, __end, __io, __err, __xtrc);
    //       std::__convert_to_v(__xtrc.c_str(), __v, __err, _S_get_c_locale());
    //       if (__beg == __end)
    // 	__err |= ios_base::eofbit;
    //       return __beg;
    //     }
    // 
    //   template<typename _CharT, typename _InIter>
    //     _InIter
    //     num_get<_CharT, _InIter>::
    //     do_get(iter_type __beg, iter_type __end, ios_base& __io,
    //            ios_base::iostate& __err, void*& __v) const
    //     {
    //       // Prepare for hex formatted input.
    //       typedef ios_base::fmtflags        fmtflags;
    //       const fmtflags __fmt = __io.flags();
    //       __io.flags((__fmt & ~ios_base::basefield) | ios_base::hex);
    // 
    //       typedef __gnu_cxx::__conditional_type<(sizeof(void*)
    // 					     <= sizeof(unsigned long)),
    // 	unsigned long, unsigned long long>::__type _UIntPtrType;       
    // 
    //       _UIntPtrType __ul;
    //       __beg = _M_extract_int(__beg, __end, __io, __err, __ul);
    // 
    //       // Reset from hex formatted input.
    //       __io.flags(__fmt);
    // 
    //       __v = reinterpret_cast<void*>(__ul);
    //       return __beg;
    //     }
    // 
    // #if defined _GLIBCXX_LONG_DOUBLE_ALT128_COMPAT \
    //       && defined __LONG_DOUBLE_IEEE128__
    //   template<typename _CharT, typename _InIter>
    //     _InIter
    //     num_get<_CharT, _InIter>::
    //     __do_get(iter_type __beg, iter_type __end, ios_base& __io,
    // 	     ios_base::iostate& __err, __ibm128& __v) const
    //     {
    //       string __xtrc;
    //       __xtrc.reserve(32);
    //       __beg = _M_extract_float(__beg, __end, __io, __err, __xtrc);
    //       std::__convert_to_v(__xtrc.c_str(), __v, __err, _S_get_c_locale());
    //       if (__beg == __end)
    // 	__err |= ios_base::eofbit;
    //       return __beg;
    //     }
    // #endif
    // 
    //   // For use by integer and floating-point types after they have been
    //   // converted into a char_type string.
    //   template<typename _CharT, typename _OutIter>
    //     void
    //     num_put<_CharT, _OutIter>::
    //     _M_pad(_CharT __fill, streamsize __w, ios_base& __io,
    // 	   _CharT* __new, const _CharT* __cs, int& __len) const
    //     {
    //       // [22.2.2.2.2] Stage 3.
    //       // If necessary, pad.
    //       __pad<_CharT, char_traits<_CharT> >::_S_pad(__io, __fill, __new,
    // 						  __cs, __w, __len);
    //       __len = static_cast<int>(__w);
    //     }
    // 
    // _GLIBCXX_END_NAMESPACE_LDBL
    // 
    //   template<typename _CharT, typename _ValueT>
    //     int
    //     __int_to_char(_CharT* __bufend, _ValueT __v, const _CharT* __lit,
    // 		  ios_base::fmtflags __flags, bool __dec)
    //     {
    //       _CharT* __buf = __bufend;
    //       if (__builtin_expect(__dec, true))
    // 	{
    // 	  // Decimal.
    // 	  do
    // 	    {
    // 	      *--__buf = __lit[(__v % 10) + __num_base::_S_odigits];
    // 	      __v /= 10;
    // 	    }
    // 	  while (__v != 0);
    // 	}
    //       else if ((__flags & ios_base::basefield) == ios_base::oct)
    // 	{
    // 	  // Octal.
    // 	  do
    // 	    {
    // 	      *--__buf = __lit[(__v & 0x7) + __num_base::_S_odigits];
    // 	      __v >>= 3;
    // 	    }
    // 	  while (__v != 0);
    // 	}
    //       else
    // 	{
    // 	  // Hex.
    // 	  const bool __uppercase = __flags & ios_base::uppercase;
    // 	  const int __case_offset = __uppercase ? __num_base::_S_oudigits
    // 	                                        : __num_base::_S_odigits;
    // 	  do
    // 	    {
    // 	      *--__buf = __lit[(__v & 0xf) + __case_offset];
    // 	      __v >>= 4;
    // 	    }
    // 	  while (__v != 0);
    // 	}
    //       return __bufend - __buf;
    //     }
    // 
    // _GLIBCXX_BEGIN_NAMESPACE_LDBL
    // 
    //   template<typename _CharT, typename _OutIter>
    //     void
    //     num_put<_CharT, _OutIter>::
    //     _M_group_int(const char* __grouping, size_t __grouping_size, _CharT __sep,
    // 		 ios_base&, _CharT* __new, _CharT* __cs, int& __len) const
    //     {
    //       _CharT* __p = std::__add_grouping(__new, __sep, __grouping,
    // 					__grouping_size, __cs, __cs + __len);
    //       __len = __p - __new;
    //     }
    //   
    //   template<typename _CharT, typename _OutIter>
    //     template<typename _ValueT>
    //       _OutIter
    //       num_put<_CharT, _OutIter>::
    //       _M_insert_int(_OutIter __s, ios_base& __io, _CharT __fill,
    // 		    _ValueT __v) const
    //       {
    // 	using __gnu_cxx::__add_unsigned;
    // 	typedef typename __add_unsigned<_ValueT>::__type __unsigned_type;
    // 	typedef __numpunct_cache<_CharT>	             __cache_type;
    // 	__use_cache<__cache_type> __uc;
    // 	const locale& __loc = __io._M_getloc();
    // 	const __cache_type* __lc = __uc(__loc);
    // 	const _CharT* __lit = __lc->_M_atoms_out;
    // 	const ios_base::fmtflags __flags = __io.flags();
    // 
    // 	// Long enough to hold hex, dec, and octal representations.
    // 	const int __ilen = 5 * sizeof(_ValueT);
    // 	_CharT* __cs = static_cast<_CharT*>(__builtin_alloca(sizeof(_CharT)
    // 							     * __ilen));
    // 
    // 	// [22.2.2.2.2] Stage 1, numeric conversion to character.
    // 	// Result is returned right-justified in the buffer.
    // 	const ios_base::fmtflags __basefield = __flags & ios_base::basefield;
    // 	const bool __dec = (__basefield != ios_base::oct
    // 			    && __basefield != ios_base::hex);
    // 	const __unsigned_type __u = ((__v > 0 || !__dec)
    // 				     ? __unsigned_type(__v)
    // 				     : -__unsigned_type(__v));
    //  	int __len = __int_to_char(__cs + __ilen, __u, __lit, __flags, __dec);
    // 	__cs += __ilen - __len;
    // 
    // 	// Add grouping, if necessary.
    // 	if (__lc->_M_use_grouping)
    // 	  {
    // 	    // Grouping can add (almost) as many separators as the number
    // 	    // of digits + space is reserved for numeric base or sign.
    // 	    _CharT* __cs2 = static_cast<_CharT*>(__builtin_alloca(sizeof(_CharT)
    // 								  * (__len + 1)
    // 								  * 2));
    // 	    _M_group_int(__lc->_M_grouping, __lc->_M_grouping_size,
    // 			 __lc->_M_thousands_sep, __io, __cs2 + 2, __cs, __len);
    // 	    __cs = __cs2 + 2;
    // 	  }
    // 
    // 	// Complete Stage 1, prepend numeric base or sign.
    // 	if (__builtin_expect(__dec, true))
    // 	  {
    // 	    // Decimal.
    // 	    if (__v >= 0)
    // 	      {
    // 		if (bool(__flags & ios_base::showpos)
    // 		    && __gnu_cxx::__numeric_traits<_ValueT>::__is_signed)
    // 		  *--__cs = __lit[__num_base::_S_oplus], ++__len;
    // 	      }
    // 	    else
    // 	      *--__cs = __lit[__num_base::_S_ominus], ++__len;
    // 	  }
    // 	else if (bool(__flags & ios_base::showbase) && __v)
    // 	  {
    // 	    if (__basefield == ios_base::oct)
    // 	      *--__cs = __lit[__num_base::_S_odigits], ++__len;
    // 	    else
    // 	      {
    // 		// 'x' or 'X'
    // 		const bool __uppercase = __flags & ios_base::uppercase;
    // 		*--__cs = __lit[__num_base::_S_ox + __uppercase];
    // 		// '0'
    // 		*--__cs = __lit[__num_base::_S_odigits];
    // 		__len += 2;
    // 	      }
    // 	  }
    // 
    // 	// Pad.
    // 	const streamsize __w = __io.width();
    // 	if (__w > static_cast<streamsize>(__len))
    // 	  {
    // 	    _CharT* __cs3 = static_cast<_CharT*>(__builtin_alloca(sizeof(_CharT)
    // 								  * __w));
    // 	    _M_pad(__fill, __w, __io, __cs3, __cs, __len);
    // 	    __cs = __cs3;
    // 	  }
    // 	__io.width(0);
    // 
    // 	// [22.2.2.2.2] Stage 4.
    // 	// Write resulting, fully-formatted string to output iterator.
    // 	return std::__write(__s, __cs, __len);
    //       }
    // 
    //   template<typename _CharT, typename _OutIter>
    //     void
    //     num_put<_CharT, _OutIter>::
    //     _M_group_float(const char* __grouping, size_t __grouping_size,
    // 		   _CharT __sep, const _CharT* __p, _CharT* __new,
    // 		   _CharT* __cs, int& __len) const
    //     {
    //       // _GLIBCXX_RESOLVE_LIB_DEFECTS
    //       // 282. What types does numpunct grouping refer to?
    //       // Add grouping, if necessary.
    //       const int __declen = __p ? __p - __cs : __len;
    //       _CharT* __p2 = std::__add_grouping(__new, __sep, __grouping,
    // 					 __grouping_size,
    // 					 __cs, __cs + __declen);
    // 
    //       // Tack on decimal part.
    //       int __newlen = __p2 - __new;
    //       if (__p)
    // 	{
    // 	  char_traits<_CharT>::copy(__p2, __p, __len - __declen);
    // 	  __newlen += __len - __declen;
    // 	}
    //       __len = __newlen;
    //     }
    // 
    //   // The following code uses vsnprintf (or vsprintf(), when
    //   // _GLIBCXX_USE_C99_STDIO is not defined) to convert floating point
    //   // values for insertion into a stream.  An optimization would be to
    //   // replace them with code that works directly on a wide buffer and
    //   // then use __pad to do the padding.  It would be good to replace
    //   // them anyway to gain back the efficiency that C++ provides by
    //   // knowing up front the type of the values to insert.  Also, sprintf
    //   // is dangerous since may lead to accidental buffer overruns.  This
    //   // implementation follows the C++ standard fairly directly as
    //   // outlined in 22.2.2.2 [lib.locale.num.put]
    //   template<typename _CharT, typename _OutIter>
    //     template<typename _ValueT>
    //       _OutIter
    //       num_put<_CharT, _OutIter>::
    //       _M_insert_float(_OutIter __s, ios_base& __io, _CharT __fill, char __mod,
    // 		       _ValueT __v) const
    //       {
    // 	typedef __numpunct_cache<_CharT>                __cache_type;
    // 	__use_cache<__cache_type> __uc;
    // 	const locale& __loc = __io._M_getloc();
    // 	const __cache_type* __lc = __uc(__loc);
    // 
    // 	// Use default precision if out of range.
    // 	const streamsize __prec = __io.precision() < 0 ? 6 : __io.precision();
    // 
    // 	const int __max_digits =
    // 	  __gnu_cxx::__numeric_traits<_ValueT>::__digits10;
    // 
    // 	// [22.2.2.2.2] Stage 1, numeric conversion to character.
    // 	int __len;
    // 	// Long enough for the max format spec.
    // 	char __fbuf[16];
    // 	__num_base::_S_format_float(__io, __fbuf, __mod);
    // 
    // #if _GLIBCXX_USE_C99_STDIO && !_GLIBCXX_HAVE_BROKEN_VSNPRINTF
    // 	// Precision is always used except for hexfloat format.
    // 	const bool __use_prec =
    // 	  (__io.flags() & ios_base::floatfield) != ios_base::floatfield;
    // 
    // 	// First try a buffer perhaps big enough (most probably sufficient
    // 	// for non-ios_base::fixed outputs)
    // 	int __cs_size = __max_digits * 3;
    // 	char* __cs = static_cast<char*>(__builtin_alloca(__cs_size));
    // 	if (__use_prec)
    // 	  __len = std::__convert_from_v(_S_get_c_locale(), __cs, __cs_size,
    // 					__fbuf, __prec, __v);
    // 	else
    // 	  __len = std::__convert_from_v(_S_get_c_locale(), __cs, __cs_size,
    // 					__fbuf, __v);
    // 
    // 	// If the buffer was not large enough, try again with the correct size.
    // 	if (__len >= __cs_size)
    // 	  {
    // 	    __cs_size = __len + 1;
    // 	    __cs = static_cast<char*>(__builtin_alloca(__cs_size));
    // 	    if (__use_prec)
    // 	      __len = std::__convert_from_v(_S_get_c_locale(), __cs, __cs_size,
    // 					    __fbuf, __prec, __v);
    // 	    else
    // 	      __len = std::__convert_from_v(_S_get_c_locale(), __cs, __cs_size,
    // 					    __fbuf, __v);
    // 	  }
    // #else
    // 	// Consider the possibility of long ios_base::fixed outputs
    // 	const bool __fixed = __io.flags() & ios_base::fixed;
    // 	const int __max_exp =
    // 	  __gnu_cxx::__numeric_traits<_ValueT>::__max_exponent10;
    // 
    // 	// The size of the output string is computed as follows.
    // 	// ios_base::fixed outputs may need up to __max_exp + 1 chars
    // 	// for the integer part + __prec chars for the fractional part
    // 	// + 3 chars for sign, decimal point, '\0'. On the other hand,
    // 	// for non-fixed outputs __max_digits * 2 + __prec chars are
    // 	// largely sufficient.
    // 	const int __cs_size = __fixed ? __max_exp + __prec + 4
    // 	                              : __max_digits * 2 + __prec;
    // 	char* __cs = static_cast<char*>(__builtin_alloca(__cs_size));
    // 	__len = std::__convert_from_v(_S_get_c_locale(), __cs, 0, __fbuf, 
    // 				      __prec, __v);
    // #endif
    // 
    // 	// [22.2.2.2.2] Stage 2, convert to char_type, using correct
    // 	// numpunct.decimal_point() values for '.' and adding grouping.
    // 	const ctype<_CharT>& __ctype = use_facet<ctype<_CharT> >(__loc);
    // 	
    // 	_CharT* __ws = static_cast<_CharT*>(__builtin_alloca(sizeof(_CharT)
    // 							     * __len));
    // 	__ctype.widen(__cs, __cs + __len, __ws);
    // 	
    // 	// Replace decimal point.
    // 	_CharT* __wp = 0;
    // 	const char* __p = char_traits<char>::find(__cs, __len, '.');
    // 	if (__p)
    // 	  {
    // 	    __wp = __ws + (__p - __cs);
    // 	    *__wp = __lc->_M_decimal_point;
    // 	  }
    // 	
    // 	// Add grouping, if necessary.
    // 	// N.B. Make sure to not group things like 2e20, i.e., no decimal
    // 	// point, scientific notation.
    // 	if (__lc->_M_use_grouping
    // 	    && (__wp || __len < 3 || (__cs[1] <= '9' && __cs[2] <= '9'
    // 				      && __cs[1] >= '0' && __cs[2] >= '0')))
    // 	  {
    // 	    // Grouping can add (almost) as many separators as the
    // 	    // number of digits, but no more.
    // 	    _CharT* __ws2 = static_cast<_CharT*>(__builtin_alloca(sizeof(_CharT)
    // 								  * __len * 2));
    // 	    
    // 	    streamsize __off = 0;
    // 	    if (__cs[0] == '-' || __cs[0] == '+')
    // 	      {
    // 		__off = 1;
    // 		__ws2[0] = __ws[0];
    // 		__len -= 1;
    // 	      }
    // 	    
    // 	    _M_group_float(__lc->_M_grouping, __lc->_M_grouping_size,
    // 			   __lc->_M_thousands_sep, __wp, __ws2 + __off,
    // 			   __ws + __off, __len);
    // 	    __len += __off;
    // 	    
    // 	    __ws = __ws2;
    // 	  }
    // 
    // 	// Pad.
    // 	const streamsize __w = __io.width();
    // 	if (__w > static_cast<streamsize>(__len))
    // 	  {
    // 	    _CharT* __ws3 = static_cast<_CharT*>(__builtin_alloca(sizeof(_CharT)
    // 								  * __w));
    // 	    _M_pad(__fill, __w, __io, __ws3, __ws, __len);
    // 	    __ws = __ws3;
    // 	  }
    // 	__io.width(0);
    // 	
    // 	// [22.2.2.2.2] Stage 4.
    // 	// Write resulting, fully-formatted string to output iterator.
    // 	return std::__write(__s, __ws, __len);
    //       }
    //   
    //   template<typename _CharT, typename _OutIter>
    //     _OutIter
    //     num_put<_CharT, _OutIter>::
    //     do_put(iter_type __s, ios_base& __io, char_type __fill, bool __v) const
    //     {
    //       const ios_base::fmtflags __flags = __io.flags();
    //       if ((__flags & ios_base::boolalpha) == 0)
    //         {
    //           const long __l = __v;
    //           __s = _M_insert_int(__s, __io, __fill, __l);
    //         }
    //       else
    //         {
    // 	  typedef __numpunct_cache<_CharT>              __cache_type;
    // 	  __use_cache<__cache_type> __uc;
    // 	  const locale& __loc = __io._M_getloc();
    // 	  const __cache_type* __lc = __uc(__loc);
    // 
    // 	  const _CharT* __name = __v ? __lc->_M_truename
    // 	                             : __lc->_M_falsename;
    // 	  int __len = __v ? __lc->_M_truename_size
    // 	                  : __lc->_M_falsename_size;
    // 
    // 	  const streamsize __w = __io.width();
    // 	  if (__w > static_cast<streamsize>(__len))
    // 	    {
    // 	      const streamsize __plen = __w - __len;
    // 	      _CharT* __ps
    // 		= static_cast<_CharT*>(__builtin_alloca(sizeof(_CharT)
    // 							* __plen));
    // 
    // 	      char_traits<_CharT>::assign(__ps, __plen, __fill);
    // 	      __io.width(0);
    // 
    // 	      if ((__flags & ios_base::adjustfield) == ios_base::left)
    // 		{
    // 		  __s = std::__write(__s, __name, __len);
    // 		  __s = std::__write(__s, __ps, __plen);
    // 		}
    // 	      else
    // 		{
    // 		  __s = std::__write(__s, __ps, __plen);
    // 		  __s = std::__write(__s, __name, __len);
    // 		}
    // 	      return __s;
    // 	    }
    // 	  __io.width(0);
    // 	  __s = std::__write(__s, __name, __len);
    // 	}
    //       return __s;
    //     }
    // 
    //   template<typename _CharT, typename _OutIter>
    //     _OutIter
    //     num_put<_CharT, _OutIter>::
    //     do_put(iter_type __s, ios_base& __io, char_type __fill, double __v) const
    //     { return _M_insert_float(__s, __io, __fill, char(), __v); }
    // 
    // #if defined _GLIBCXX_LONG_DOUBLE_COMPAT && defined __LONG_DOUBLE_128__
    //   template<typename _CharT, typename _OutIter>
    //     _OutIter
    //     num_put<_CharT, _OutIter>::
    //     __do_put(iter_type __s, ios_base& __io, char_type __fill, double __v) const
    //     { return _M_insert_float(__s, __io, __fill, char(), __v); }
    // #endif
    // 
    //   template<typename _CharT, typename _OutIter>
    //     _OutIter
    //     num_put<_CharT, _OutIter>::
    //     do_put(iter_type __s, ios_base& __io, char_type __fill,
    // 	   long double __v) const
    //     { return _M_insert_float(__s, __io, __fill, 'L', __v); }
    // 
    //   template<typename _CharT, typename _OutIter>
    //     _OutIter
    //     num_put<_CharT, _OutIter>::
    //     do_put(iter_type __s, ios_base& __io, char_type __fill,
    //            const void* __v) const
    //     {
    //       const ios_base::fmtflags __flags = __io.flags();
    //       const ios_base::fmtflags __fmt = ~(ios_base::basefield
    // 					 | ios_base::uppercase);
    //       __io.flags((__flags & __fmt) | (ios_base::hex | ios_base::showbase));
    // 
    //       typedef __gnu_cxx::__conditional_type<(sizeof(const void*)
    // 					     <= sizeof(unsigned long)),
    // 	unsigned long, unsigned long long>::__type _UIntPtrType;       
    // 
    //       __s = _M_insert_int(__s, __io, __fill,
    // 			  reinterpret_cast<_UIntPtrType>(__v));
    //       __io.flags(__flags);
    //       return __s;
    //     }
    // 
    // #if defined _GLIBCXX_LONG_DOUBLE_ALT128_COMPAT \
    //       && defined __LONG_DOUBLE_IEEE128__
    //   template<typename _CharT, typename _OutIter>
    //     _OutIter
    //     num_put<_CharT, _OutIter>::
    //     __do_put(iter_type __s, ios_base& __io, char_type __fill,
    // 	     __ibm128 __v) const
    //     { return _M_insert_float(__s, __io, __fill, 'L', __v); }
    // #endif
    // _GLIBCXX_END_NAMESPACE_LDBL
    // 
    //   // Construct correctly padded string, as per 22.2.2.2.2
    //   // Assumes
    //   // __newlen > __oldlen
    //   // __news is allocated for __newlen size
    // 
    //   // NB: Of the two parameters, _CharT can be deduced from the
    //   // function arguments. The other (_Traits) has to be explicitly specified.
    //   template<typename _CharT, typename _Traits>
    //     void
    //     __pad<_CharT, _Traits>::_S_pad(ios_base& __io, _CharT __fill,
    // 				   _CharT* __news, const _CharT* __olds,
    // 				   streamsize __newlen, streamsize __oldlen)
    //     {
    //       const size_t __plen = static_cast<size_t>(__newlen - __oldlen);
    //       const ios_base::fmtflags __adjust = __io.flags() & ios_base::adjustfield;
    // 
    //       // Padding last.
    //       if (__adjust == ios_base::left)
    // 	{
    // 	  _Traits::copy(__news, __olds, __oldlen);
    // 	  _Traits::assign(__news + __oldlen, __plen, __fill);
    // 	  return;
    // 	}
    // 
    //       size_t __mod = 0;
    //       if (__adjust == ios_base::internal)
    // 	{
    // 	  // Pad after the sign, if there is one.
    // 	  // Pad after 0[xX], if there is one.
    // 	  // Who came up with these rules, anyway? Jeeze.
    //           const locale& __loc = __io._M_getloc();
    // 	  const ctype<_CharT>& __ctype = use_facet<ctype<_CharT> >(__loc);
    // 
    // 	  if (__ctype.widen('-') == __olds[0]
    // 	      || __ctype.widen('+') == __olds[0])
    // 	    {
    // 	      __news[0] = __olds[0];
    // 	      __mod = 1;
    // 	      ++__news;
    // 	    }
    // 	  else if (__ctype.widen('0') == __olds[0]
    // 		   && __oldlen > 1
    // 		   && (__ctype.widen('x') == __olds[1]
    // 		       || __ctype.widen('X') == __olds[1]))
    // 	    {
    // 	      __news[0] = __olds[0];
    // 	      __news[1] = __olds[1];
    // 	      __mod = 2;
    // 	      __news += 2;
    // 	    }
    // 	  // else Padding first.
    // 	}
    //       _Traits::assign(__news, __plen, __fill);
    //       _Traits::copy(__news + __plen, __olds + __mod, __oldlen - __mod);
    //     }
    // 
    //   template<typename _CharT>
    //     _CharT*
    //     __add_grouping(_CharT* __s, _CharT __sep,
    // 		   const char* __gbeg, size_t __gsize,
    // 		   const _CharT* __first, const _CharT* __last)
    //     {
    //       size_t __idx = 0;
    //       size_t __ctr = 0;
    // 
    //       while (__last - __first > __gbeg[__idx]
    // 	     && static_cast<signed char>(__gbeg[__idx]) > 0
    // 	     && __gbeg[__idx] != __gnu_cxx::__numeric_traits<char>::__max)
    // 	{
    // 	  __last -= __gbeg[__idx];
    // 	  __idx < __gsize - 1 ? ++__idx : ++__ctr;
    // 	}
    // 
    //       while (__first != __last)
    // 	*__s++ = *__first++;
    // 
    //       while (__ctr--)
    // 	{
    // 	  *__s++ = __sep;	  
    // 	  for (char __i = __gbeg[__idx]; __i > 0; --__i)
    // 	    *__s++ = *__first++;
    // 	}
    // 
    //       while (__idx--)
    // 	{
    // 	  *__s++ = __sep;	  
    // 	  for (char __i = __gbeg[__idx]; __i > 0; --__i)
    // 	    *__s++ = *__first++;
    // 	}
    // 
    //       return __s;
    //     }
    // 
    //   // Inhibit implicit instantiations for required instantiations,
    //   // which are defined via explicit instantiations elsewhere.
    // #if _GLIBCXX_EXTERN_TEMPLATE
    //   extern template class _GLIBCXX_NAMESPACE_CXX11 numpunct<char>;
    //   extern template class _GLIBCXX_NAMESPACE_CXX11 numpunct_byname<char>;
    //   extern template class _GLIBCXX_NAMESPACE_LDBL num_get<char>;
    //   extern template class _GLIBCXX_NAMESPACE_LDBL num_put<char>;
    //   extern template class ctype_byname<char>;
    // 
    //   extern template
    //     const ctype<char>*
    //     __try_use_facet<ctype<char> >(const locale&) _GLIBCXX_NOTHROW;
    // 
    //   extern template
    //     const numpunct<char>*
    //     __try_use_facet<numpunct<char> >(const locale&) _GLIBCXX_NOTHROW;
    // 
    //   extern template
    //     const num_put<char>*
    //     __try_use_facet<num_put<char> >(const locale&) _GLIBCXX_NOTHROW;
    // 
    //   extern template
    //     const num_get<char>*
    //     __try_use_facet<num_get<char> >(const locale&) _GLIBCXX_NOTHROW;
    // 
    //   extern template
    //     const ctype<char>&
    //     use_facet<ctype<char> >(const locale&);
    // 
    //   extern template
    //     const numpunct<char>&
    //     use_facet<numpunct<char> >(const locale&);
    // 
    //   extern template
    //     const num_put<char>&
    //     use_facet<num_put<char> >(const locale&);
    // 
    //   extern template
    //     const num_get<char>&
    //     use_facet<num_get<char> >(const locale&);
    // 
    //   extern template
    //     bool
    //     has_facet<ctype<char> >(const locale&);
    // 
    //   extern template
    //     bool
    //     has_facet<numpunct<char> >(const locale&);
    // 
    //   extern template
    //     bool
    //     has_facet<num_put<char> >(const locale&);
    // 
    //   extern template
    //     bool
    //     has_facet<num_get<char> >(const locale&);
    // 
    // #ifdef _GLIBCXX_USE_WCHAR_T
    //   extern template class _GLIBCXX_NAMESPACE_CXX11 numpunct<wchar_t>;
    //   extern template class _GLIBCXX_NAMESPACE_CXX11 numpunct_byname<wchar_t>;
    //   extern template class _GLIBCXX_NAMESPACE_LDBL num_get<wchar_t>;
    //   extern template class _GLIBCXX_NAMESPACE_LDBL num_put<wchar_t>;
    //   extern template class ctype_byname<wchar_t>;
    // 
    //   extern template
    //     const ctype<wchar_t>*
    //     __try_use_facet<ctype<wchar_t> >(const locale&) _GLIBCXX_NOTHROW;
    // 
    //   extern template
    //     const numpunct<wchar_t>*
    //     __try_use_facet<numpunct<wchar_t> >(const locale&) _GLIBCXX_NOTHROW;
    // 
    //   extern template
    //     const num_put<wchar_t>*
    //     __try_use_facet<num_put<wchar_t> >(const locale&) _GLIBCXX_NOTHROW;
    // 
    //   extern template
    //     const num_get<wchar_t>*
    //     __try_use_facet<num_get<wchar_t> >(const locale&) _GLIBCXX_NOTHROW;
    // 
    //   extern template
    //     const ctype<wchar_t>&
    //     use_facet<ctype<wchar_t> >(const locale&);
    // 
    //   extern template
    //     const numpunct<wchar_t>&
    //     use_facet<numpunct<wchar_t> >(const locale&);
    // 
    //   extern template
    //     const num_put<wchar_t>&
    //     use_facet<num_put<wchar_t> >(const locale&);
    // 
    //   extern template
    //     const num_get<wchar_t>&
    //     use_facet<num_get<wchar_t> >(const locale&);
    // 
    //   extern template
    //     bool
    //     has_facet<ctype<wchar_t> >(const locale&);
    // 
    //   extern template
    //     bool
    //     has_facet<numpunct<wchar_t> >(const locale&);
    // 
    //   extern template
    //     bool
    //     has_facet<num_put<wchar_t> >(const locale&);
    // 
    //   extern template
    //     bool
    //     has_facet<num_get<wchar_t> >(const locale&);
    // #endif
    // #endif
    // 
    // _GLIBCXX_END_NAMESPACE_VERSION
    // } // namespace
    // 
    // #endif
    // // Wrapper for underlying C-language localization -*- C++ -*-
    // 
    // // Copyright (C) 2001-2024 Free Software Foundation, Inc.
    // //
    // // This file is part of the GNU ISO C++ Library.  This library is free
    // // software; you can redistribute it and/or modify it under the
    // // terms of the GNU General Public License as published by the
    // // Free Software Foundation; either version 3, or (at your option)
    // // any later version.
    // 
    // // This library is distributed in the hope that it will be useful,
    // // but WITHOUT ANY WARRANTY; without even the implied warranty of
    // // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    // // GNU General Public License for more details.
    // 
    // // Under Section 7 of GPL version 3, you are granted additional
    // // permissions described in the GCC Runtime Library Exception, version
    // // 3.1, as published by the Free Software Foundation.
    // 
    // // You should have received a copy of the GNU General Public License and
    // // a copy of the GCC Runtime Library Exception along with this program;
    // // see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
    // // <http://www.gnu.org/licenses/>.
    // 
    // //
    // // ISO C++ 14882: 22.8  Standard locale categories.
    // //
    // 
    // // Written by Benjamin Kosnik <bkoz@redhat.com>
    // 
    // #include <locale>
    // #include <stdexcept>
    // #include <limits>
    // #include <algorithm>
    // #include <langinfo.h>
    // #include <bits/c++locale_internal.h>
    // 
    // #include <backward/auto_ptr.h>
    // 
    // namespace std _GLIBCXX_VISIBILITY(default)
    // {
    // _GLIBCXX_BEGIN_NAMESPACE_VERSION
    // 
    //   template<>
    //     void
    //     __convert_to_v(const char* __s, float& __v, ios_base::iostate& __err,
    // 		   const __c_locale& __cloc) throw()
    //     {
    //       char* __sanity;
    //       __v = __strtof_l(__s, &__sanity, __cloc);
    // 
    //       // _GLIBCXX_RESOLVE_LIB_DEFECTS
    //       // 23. Num_get overflow result.
    //       if (__sanity == __s || *__sanity != '\0')
    // 	{
    // 	  __v = 0.0f;
    // 	  __err = ios_base::failbit;
    // 	}
    //       else if (__v == numeric_limits<float>::infinity())
    // 	{
    // 	  __v = numeric_limits<float>::max();
    // 	  __err = ios_base::failbit;
    // 	}
    //       else if (__v == -numeric_limits<float>::infinity())
    // 	{
    // 	  __v = -numeric_limits<float>::max();
    // 	  __err = ios_base::failbit;
    // 	}
    //     }
    // 
    //   template<>
    //     void
    //     __convert_to_v(const char* __s, double& __v, ios_base::iostate& __err,
    // 		   const __c_locale& __cloc) throw()
    //     {
    //       char* __sanity;
    //       __v = __strtod_l(__s, &__sanity, __cloc);
    // 
    //       // _GLIBCXX_RESOLVE_LIB_DEFECTS
    //       // 23. Num_get overflow result.
    //       if (__sanity == __s || *__sanity != '\0')
    // 	{
    // 	  __v = 0.0;
    // 	  __err = ios_base::failbit;
    // 	}
    //       else if (__v == numeric_limits<double>::infinity())
    // 	{
    // 	  __v = numeric_limits<double>::max();
    // 	  __err = ios_base::failbit;
    // 	}
    //       else if (__v == -numeric_limits<double>::infinity())
    // 	{
    // 	  __v = -numeric_limits<double>::max();
    // 	  __err = ios_base::failbit;
    // 	}
    //     }
    // 
    //   template<>
    //     void
    //     __convert_to_v(const char* __s, long double& __v, ios_base::iostate& __err,
    // 		   const __c_locale& __cloc) throw()
    //     {
    //       char* __sanity;
    // #if __GLIBC__ > 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ > 2)
    //       // Prefer strtold_l, as __strtold_l isn't prototyped in more recent
    //       // glibc versions.
    //       __v = strtold_l(__s, &__sanity, __cloc);
    // #else
    //       __v = __strtold_l(__s, &__sanity, __cloc);
    // #endif
    // 
    //       // _GLIBCXX_RESOLVE_LIB_DEFECTS
    //       // 23. Num_get overflow result.
    //       if (__sanity == __s || *__sanity != '\0')
    // 	{
    // 	  __v = 0.0l;
    // 	  __err = ios_base::failbit;
    // 	}
    //       else if (__v == numeric_limits<long double>::infinity())
    // 	{
    // 	  __v = numeric_limits<long double>::max();
    // 	  __err = ios_base::failbit;
    // 	}
    //       else if (__v == -numeric_limits<long double>::infinity())
    // 	{
    // 	  __v = -numeric_limits<long double>::max();
    // 	  __err = ios_base::failbit;
    // 	}
    //     }
    // 
    //   void
    //   locale::facet::_S_create_c_locale(__c_locale& __cloc, const char* __s,
    // 				    __c_locale __old)
    //   {
    //     __cloc = __newlocale(1 << LC_ALL, __s, __old);
    //     if (!__cloc)
    //       {
    // 	// This named locale is not supported by the underlying OS.
    // 	__throw_runtime_error(__N("locale::facet::_S_create_c_locale "
    // 				  "name not valid"));
    //       }
    //   }
    // 
    //   void
    //   locale::facet::_S_destroy_c_locale(__c_locale& __cloc)
    //   {
    //     if (__cloc && _S_get_c_locale() != __cloc)
    //       __freelocale(__cloc);
    //   }
    // 
    //   __c_locale
    //   locale::facet::_S_clone_c_locale(__c_locale& __cloc) throw()
    //   { return __duplocale(__cloc); }
    // 
    //   __c_locale
    //   locale::facet::_S_lc_ctype_c_locale(__c_locale __cloc, const char* __s)
    //   {
    //     __c_locale __dup = __duplocale(__cloc);
    //     if (__dup == __c_locale(0))
    //       __throw_runtime_error(__N("locale::facet::_S_lc_ctype_c_locale "
    // 				"duplocale error"));
    // #if __GLIBC__ > 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ > 2)
    //     __c_locale __changed = __newlocale(LC_CTYPE_MASK, __s, __dup);
    // #else
    //     __c_locale __changed = __newlocale(1 << LC_CTYPE, __s, __dup);
    // #endif
    //     if (__changed == __c_locale(0))
    //       {
    // 	__freelocale(__dup);
    // 	__throw_runtime_error(__N("locale::facet::_S_lc_ctype_c_locale "
    // 				  "newlocale error"));
    //       }
    //     return __changed;
    //   }
    // 
    //   struct _CatalogIdComp
    //   {
    //     bool
    //     operator()(messages_base::catalog __cat, const Catalog_info* __info) const
    //     { return __cat < __info->_M_id; }
    // 
    //     bool
    //     operator()(const Catalog_info* __info, messages_base::catalog __cat) const
    //     { return __info->_M_id < __cat; }
    //   };
    // 
    //   Catalogs::~Catalogs()
    //   {
    //     for (vector<Catalog_info*>::iterator __it = _M_infos.begin();
    // 	 __it != _M_infos.end(); ++__it)
    //       delete *__it;
    //   }
    // 
    //   messages_base::catalog
    //   Catalogs::_M_add(const char* __domain, locale __l)
    //   {
    //     __gnu_cxx::__scoped_lock lock(_M_mutex);
    // 
    //     // The counter is not likely to roll unless catalogs keep on being
    //     // opened/closed which is consider as an application mistake for the
    //     // moment.
    //     if (_M_catalog_counter == numeric_limits<messages_base::catalog>::max())
    //       return -1;
    // 
    //     auto_ptr<Catalog_info> info(new Catalog_info(_M_catalog_counter++,
    // 						 __domain, __l));
    // 
    //     // Check if we managed to allocate memory for domain.
    //     if (!info->_M_domain)
    //       return -1;
    // 
    //     _M_infos.push_back(info.get());
    //     return info.release()->_M_id;
    //   }
    // 
    //   void
    //   Catalogs::_M_erase(messages_base::catalog __c)
    //   {
    //     __gnu_cxx::__scoped_lock lock(_M_mutex);
    // 
    //     vector<Catalog_info*>::iterator __res =
    //       lower_bound(_M_infos.begin(), _M_infos.end(), __c, _CatalogIdComp());
    //     if (__res == _M_infos.end() || (*__res)->_M_id != __c)
    //       return;
    // 
    //     delete *__res;
    //     _M_infos.erase(__res);
    // 
    //     // Just in case closed catalog was the last open.
    //     if (__c == _M_catalog_counter - 1)
    //       --_M_catalog_counter;
    //   }
    // 
    //   const Catalog_info*
    //   Catalogs::_M_get(messages_base::catalog __c) const
    //   {
    //     __gnu_cxx::__scoped_lock lock(_M_mutex);
    // 
    //     vector<Catalog_info*>::const_iterator __res =
    //       lower_bound(_M_infos.begin(), _M_infos.end(), __c, _CatalogIdComp());
    // 
    //     if (__res != _M_infos.end() && (*__res)->_M_id == __c)
    //       return *__res;
    // 
    //     return 0;
    //   }
    // 
    //   Catalogs&
    //   get_catalogs()
    //   {
    //     static Catalogs __catalogs;
    //     return __catalogs;
    //   }
    // 
    // _GLIBCXX_END_NAMESPACE_VERSION
    // } // namespace
    // 
    // namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
    // {
    // _GLIBCXX_BEGIN_NAMESPACE_VERSION
    // 
    //   const char* const category_names[6 + _GLIBCXX_NUM_CATEGORIES] =
    //     {
    //       "LC_CTYPE",
    //       "LC_NUMERIC",
    //       "LC_TIME",
    //       "LC_COLLATE",
    //       "LC_MONETARY",
    //       "LC_MESSAGES",
    //       "LC_PAPER",
    //       "LC_NAME",
    //       "LC_ADDRESS",
    //       "LC_TELEPHONE",
    //       "LC_MEASUREMENT",
    //       "LC_IDENTIFICATION"
    //     };
    // 
    // _GLIBCXX_END_NAMESPACE_VERSION
    // } // namespace
    // 
    // namespace std _GLIBCXX_VISIBILITY(default)
    // {
    // _GLIBCXX_BEGIN_NAMESPACE_VERSION
    // 
    //   const char* const* const locale::_S_categories = __gnu_cxx::category_names;
    // 
    // _GLIBCXX_END_NAMESPACE_VERSION
    // } // namespace
    // 
    // // XXX GLIBCXX_ABI Deprecated
    // #ifdef _GLIBCXX_LONG_DOUBLE_COMPAT
    // #pragma GCC diagnostic ignored "-Wattribute-alias"
    // #define _GLIBCXX_LDBL_COMPAT(dbl, ldbl) \
    //   extern "C" void ldbl (void) __attribute__ ((alias (#dbl)))
    // _GLIBCXX_LDBL_COMPAT(_ZSt14__convert_to_vIdEvPKcRT_RSt12_Ios_IostateRKP15__locale_struct, _ZSt14__convert_to_vIeEvPKcRT_RSt12_Ios_IostateRKP15__locale_struct);
    // #endif // _GLIBCXX_LONG_DOUBLE_COMPAT
    // // istream classes -*- C++ -*-
    // 
    // // Copyright (C) 1997-2024 Free Software Foundation, Inc.
    // //
    // // This file is part of the GNU ISO C++ Library.  This library is free
    // // software; you can redistribute it and/or modify it under the
    // // terms of the GNU General Public License as published by the
    // // Free Software Foundation; either version 3, or (at your option)
    // // any later version.
    // 
    // // This library is distributed in the hope that it will be useful,
    // // but WITHOUT ANY WARRANTY; without even the implied warranty of
    // // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    // // GNU General Public License for more details.
    // 
    // // Under Section 7 of GPL version 3, you are granted additional
    // // permissions described in the GCC Runtime Library Exception, version
    // // 3.1, as published by the Free Software Foundation.
    // 
    // // You should have received a copy of the GNU General Public License and
    // // a copy of the GCC Runtime Library Exception along with this program;
    // // see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
    // // <http://www.gnu.org/licenses/>.
    // 
    // /** @file bits/istream.tcc
    //  *  This is an internal header file, included by other library headers.
    //  *  Do not attempt to use it directly. @headername{istream}
    //  */
    // 
    // //
    // // ISO C++ 14882: 27.6.1  Input streams
    // //
    // 
    // #ifndef _ISTREAM_TCC
    // #define _ISTREAM_TCC 1
    // 
    // #pragma GCC system_header
    // 
    // #include <bits/cxxabi_forced.h>
    // 
    // namespace std _GLIBCXX_VISIBILITY(default)
    // {
    // _GLIBCXX_BEGIN_NAMESPACE_VERSION
    // 
    //   template<typename _CharT, typename _Traits>
    //     basic_istream<_CharT, _Traits>::sentry::
    //     sentry(basic_istream<_CharT, _Traits>& __in, bool __noskip) : _M_ok(false)
    //     {
    //       ios_base::iostate __err = ios_base::goodbit;
    //       if (__in.good())
    // 	{
    // 	  __try
    // 	    {
    // 	      if (__in.tie())
    // 		__in.tie()->flush();
    // 	      if (!__noskip && bool(__in.flags() & ios_base::skipws))
    // 		{
    // 		  const __int_type __eof = traits_type::eof();
    // 		  __streambuf_type* __sb = __in.rdbuf();
    // 		  __int_type __c = __sb->sgetc();
    // 
    // 		  const __ctype_type& __ct = __check_facet(__in._M_ctype);
    // 		  while (!traits_type::eq_int_type(__c, __eof)
    // 			 && __ct.is(ctype_base::space,
    // 				    traits_type::to_char_type(__c)))
    // 		    __c = __sb->snextc();
    // 
    // 		  // _GLIBCXX_RESOLVE_LIB_DEFECTS
    // 		  // 195. Should basic_istream::sentry's constructor ever
    // 		  // set eofbit?
    // 		  if (traits_type::eq_int_type(__c, __eof))
    // 		    __err |= ios_base::eofbit;
    // 		}
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      __in._M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { __in._M_setstate(ios_base::badbit); }
    // 	}
    // 
    //       if (__in.good() && __err == ios_base::goodbit)
    // 	_M_ok = true;
    //       else
    // 	{
    // 	  __err |= ios_base::failbit;
    // 	  __in.setstate(__err);
    // 	}
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     template<typename _ValueT>
    //       basic_istream<_CharT, _Traits>&
    //       basic_istream<_CharT, _Traits>::
    //       _M_extract(_ValueT& __v)
    //       {
    // 	sentry __cerb(*this, false);
    // 	if (__cerb)
    // 	  {
    // 	    ios_base::iostate __err = ios_base::goodbit;
    // 	    __try
    // 	      {
    // #ifndef _GLIBCXX_LONG_DOUBLE_ALT128_COMPAT
    // 		const __num_get_type& __ng = __check_facet(this->_M_num_get);
    // #else
    // 		const __num_get_type& __ng
    // 		  = use_facet<__num_get_type>(this->_M_ios_locale);
    // #endif
    // 		__ng.get(*this, 0, *this, __err, __v);
    // 	      }
    // 	    __catch(__cxxabiv1::__forced_unwind&)
    // 	      {
    // 		this->_M_setstate(ios_base::badbit);
    // 		__throw_exception_again;
    // 	      }
    // 	    __catch(...)
    // 	      { this->_M_setstate(ios_base::badbit); }
    // 	    if (__err)
    // 	      this->setstate(__err);
    // 	  }
    // 	return *this;
    //       }
    // 
    //   template<typename _CharT, typename _Traits>
    //     basic_istream<_CharT, _Traits>&
    //     basic_istream<_CharT, _Traits>::
    //     operator>>(short& __n)
    //     {
    //       // _GLIBCXX_RESOLVE_LIB_DEFECTS
    //       // 118. basic_istream uses nonexistent num_get member functions.
    //       sentry __cerb(*this, false);
    //       if (__cerb)
    // 	{
    // 	  ios_base::iostate __err = ios_base::goodbit;
    // 	  __try
    // 	    {
    // 	      long __l;
    // #ifndef _GLIBCXX_LONG_DOUBLE_ALT128_COMPAT
    // 	      const __num_get_type& __ng = __check_facet(this->_M_num_get);
    // #else
    // 	      const __num_get_type& __ng
    // 		= use_facet<__num_get_type>(this->_M_ios_locale);
    // #endif
    // 	      __ng.get(*this, 0, *this, __err, __l);
    // 
    // 	      // _GLIBCXX_RESOLVE_LIB_DEFECTS
    // 	      // 696. istream::operator>>(int&) broken.
    // 	      if (__l < __gnu_cxx::__numeric_traits<short>::__min)
    // 		{
    // 		  __err |= ios_base::failbit;
    // 		  __n = __gnu_cxx::__numeric_traits<short>::__min;
    // 		}
    // 	      else if (__l > __gnu_cxx::__numeric_traits<short>::__max)
    // 		{
    // 		  __err |= ios_base::failbit;
    // 		  __n = __gnu_cxx::__numeric_traits<short>::__max;
    // 		}
    // 	      else
    // 		__n = short(__l);
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { this->_M_setstate(ios_base::badbit); }
    // 	  if (__err)
    // 	    this->setstate(__err);
    // 	}
    //       return *this;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     basic_istream<_CharT, _Traits>&
    //     basic_istream<_CharT, _Traits>::
    //     operator>>(int& __n)
    //     {
    //       // _GLIBCXX_RESOLVE_LIB_DEFECTS
    //       // 118. basic_istream uses nonexistent num_get member functions.
    //       sentry __cerb(*this, false);
    //       if (__cerb)
    // 	{
    // 	  ios_base::iostate __err = ios_base::goodbit;
    // 	  __try
    // 	    {
    // 	      long __l;
    // #ifndef _GLIBCXX_LONG_DOUBLE_ALT128_COMPAT
    // 	      const __num_get_type& __ng = __check_facet(this->_M_num_get);
    // #else
    // 	      const __num_get_type& __ng
    // 		= use_facet<__num_get_type>(this->_M_ios_locale);
    // #endif
    // 	      __ng.get(*this, 0, *this, __err, __l);
    // 
    // 	      // _GLIBCXX_RESOLVE_LIB_DEFECTS
    // 	      // 696. istream::operator>>(int&) broken.
    // 	      if (__l < __gnu_cxx::__numeric_traits<int>::__min)
    // 		{
    // 		  __err |= ios_base::failbit;
    // 		  __n = __gnu_cxx::__numeric_traits<int>::__min;
    // 		}
    // 	      else if (__l > __gnu_cxx::__numeric_traits<int>::__max)
    // 		{
    // 		  __err |= ios_base::failbit;	      
    // 		  __n = __gnu_cxx::__numeric_traits<int>::__max;
    // 		}
    // 	      else
    // 		__n = int(__l);
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { this->_M_setstate(ios_base::badbit); }
    // 	  if (__err)
    // 	    this->setstate(__err);
    // 	}
    //       return *this;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     basic_istream<_CharT, _Traits>&
    //     basic_istream<_CharT, _Traits>::
    //     operator>>(__streambuf_type* __sbout)
    //     {
    //       ios_base::iostate __err = ios_base::goodbit;
    //       sentry __cerb(*this, false);
    //       if (__cerb && __sbout)
    // 	{
    // 	  __try
    // 	    {
    // 	      bool __ineof;
    // 	      if (!__copy_streambufs_eof(this->rdbuf(), __sbout, __ineof))
    // 		__err |= ios_base::failbit;
    // 	      if (__ineof)
    // 		__err |= ios_base::eofbit;
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::failbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { this->_M_setstate(ios_base::failbit); }
    // 	}
    //       else if (!__sbout)
    // 	__err |= ios_base::failbit;
    //       if (__err)
    // 	this->setstate(__err);
    //       return *this;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     typename basic_istream<_CharT, _Traits>::int_type
    //     basic_istream<_CharT, _Traits>::
    //     get(void)
    //     {
    //       const int_type __eof = traits_type::eof();
    //       int_type __c = __eof;
    //       _M_gcount = 0;
    //       ios_base::iostate __err = ios_base::goodbit;
    //       sentry __cerb(*this, true);
    //       if (__cerb)
    // 	{
    // 	  __try
    // 	    {
    // 	      __c = this->rdbuf()->sbumpc();
    // 	      // 27.6.1.1 paragraph 3
    // 	      if (!traits_type::eq_int_type(__c, __eof))
    // 		_M_gcount = 1;
    // 	      else
    // 		__err |= ios_base::eofbit;
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { this->_M_setstate(ios_base::badbit); }
    // 	}
    //       if (!_M_gcount)
    // 	__err |= ios_base::failbit;
    //       if (__err)
    // 	this->setstate(__err);
    //       return __c;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     basic_istream<_CharT, _Traits>&
    //     basic_istream<_CharT, _Traits>::
    //     get(char_type& __c)
    //     {
    //       _M_gcount = 0;
    //       ios_base::iostate __err = ios_base::goodbit;
    //       sentry __cerb(*this, true);
    //       if (__cerb)
    // 	{
    // 	  __try
    // 	    {
    // 	      const int_type __cb = this->rdbuf()->sbumpc();
    // 	      // 27.6.1.1 paragraph 3
    // 	      if (!traits_type::eq_int_type(__cb, traits_type::eof()))
    // 		{
    // 		  _M_gcount = 1;
    // 		  __c = traits_type::to_char_type(__cb);
    // 		}
    // 	      else
    // 		__err |= ios_base::eofbit;
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { this->_M_setstate(ios_base::badbit); }
    // 	}
    //       if (!_M_gcount)
    // 	__err |= ios_base::failbit;
    //       if (__err)
    // 	this->setstate(__err);
    //       return *this;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     basic_istream<_CharT, _Traits>&
    //     basic_istream<_CharT, _Traits>::
    //     get(char_type* __s, streamsize __n, char_type __delim)
    //     {
    //       _M_gcount = 0;
    //       ios_base::iostate __err = ios_base::goodbit;
    //       sentry __cerb(*this, true);
    //       if (__cerb)
    // 	{
    // 	  __try
    // 	    {
    // 	      const int_type __idelim = traits_type::to_int_type(__delim);
    // 	      const int_type __eof = traits_type::eof();
    // 	      __streambuf_type* __sb = this->rdbuf();
    // 	      int_type __c = __sb->sgetc();
    // 
    // 	      while (_M_gcount + 1 < __n
    // 		     && !traits_type::eq_int_type(__c, __eof)
    // 		     && !traits_type::eq_int_type(__c, __idelim))
    // 		{
    // 		  *__s++ = traits_type::to_char_type(__c);
    // 		  ++_M_gcount;
    // 		  __c = __sb->snextc();
    // 		}
    // 	      if (traits_type::eq_int_type(__c, __eof))
    // 		__err |= ios_base::eofbit;
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { this->_M_setstate(ios_base::badbit); }
    // 	}
    //       // _GLIBCXX_RESOLVE_LIB_DEFECTS
    //       // 243. get and getline when sentry reports failure.
    //       if (__n > 0)
    // 	*__s = char_type();
    //       if (!_M_gcount)
    // 	__err |= ios_base::failbit;
    //       if (__err)
    // 	this->setstate(__err);
    //       return *this;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     basic_istream<_CharT, _Traits>&
    //     basic_istream<_CharT, _Traits>::
    //     get(__streambuf_type& __sb, char_type __delim)
    //     {
    //       _M_gcount = 0;
    //       ios_base::iostate __err = ios_base::goodbit;
    //       sentry __cerb(*this, true);
    //       if (__cerb)
    // 	{
    // 	  __try
    // 	    {
    // 	      const int_type __idelim = traits_type::to_int_type(__delim);
    // 	      const int_type __eof = traits_type::eof();
    // 	      __streambuf_type* __this_sb = this->rdbuf();
    // 	      int_type __c = __this_sb->sgetc();
    // 	      char_type __c2 = traits_type::to_char_type(__c);
    // 	      unsigned long long __gcount = 0;
    // 
    // 	      while (!traits_type::eq_int_type(__c, __eof)
    // 		     && !traits_type::eq_int_type(__c, __idelim)
    // 		     && !traits_type::eq_int_type(__sb.sputc(__c2), __eof))
    // 		{
    // 		  ++__gcount;
    // 		  __c = __this_sb->snextc();
    // 		  __c2 = traits_type::to_char_type(__c);
    // 		}
    // 	      if (traits_type::eq_int_type(__c, __eof))
    // 		__err |= ios_base::eofbit;
    // 	      // _GLIBCXX_RESOLVE_LIB_DEFECTS
    // 	      // 3464. istream::gcount() can overflow
    // 	      if (__gcount <= __gnu_cxx::__numeric_traits<streamsize>::__max)
    // 		_M_gcount = __gcount;
    // 	      else
    // 		_M_gcount = __gnu_cxx::__numeric_traits<streamsize>::__max;
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { this->_M_setstate(ios_base::badbit); }
    // 	}
    //       if (!_M_gcount)
    // 	__err |= ios_base::failbit;
    //       if (__err)
    // 	this->setstate(__err);
    //       return *this;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     basic_istream<_CharT, _Traits>&
    //     basic_istream<_CharT, _Traits>::
    //     getline(char_type* __s, streamsize __n, char_type __delim)
    //     {
    //       _M_gcount = 0;
    //       ios_base::iostate __err = ios_base::goodbit;
    //       sentry __cerb(*this, true);
    //       if (__cerb)
    //         {
    //           __try
    //             {
    //               const int_type __idelim = traits_type::to_int_type(__delim);
    //               const int_type __eof = traits_type::eof();
    //               __streambuf_type* __sb = this->rdbuf();
    //               int_type __c = __sb->sgetc();
    // 
    //               while (_M_gcount + 1 < __n
    //                      && !traits_type::eq_int_type(__c, __eof)
    //                      && !traits_type::eq_int_type(__c, __idelim))
    //                 {
    //                   *__s++ = traits_type::to_char_type(__c);
    //                   __c = __sb->snextc();
    //                   ++_M_gcount;
    //                 }
    //               if (traits_type::eq_int_type(__c, __eof))
    //                 __err |= ios_base::eofbit;
    //               else
    //                 {
    //                   if (traits_type::eq_int_type(__c, __idelim))
    //                     {
    //                       __sb->sbumpc();
    //                       ++_M_gcount;
    //                     }
    //                   else
    //                     __err |= ios_base::failbit;
    //                 }
    //             }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    //           __catch(...)
    //             { this->_M_setstate(ios_base::badbit); }
    //         }
    //       // _GLIBCXX_RESOLVE_LIB_DEFECTS
    //       // 243. get and getline when sentry reports failure.
    //       if (__n > 0)
    // 	*__s = char_type();
    //       if (!_M_gcount)
    //         __err |= ios_base::failbit;
    //       if (__err)
    //         this->setstate(__err);
    //       return *this;
    //     }
    // 
    //   // We provide three overloads, since the first two are much simpler
    //   // than the general case. Also, the latter two can thus adopt the
    //   // same "batchy" strategy used by getline above.
    //   template<typename _CharT, typename _Traits>
    //     basic_istream<_CharT, _Traits>&
    //     basic_istream<_CharT, _Traits>::
    //     ignore(void)
    //     {
    //       _M_gcount = 0;
    //       sentry __cerb(*this, true);
    //       if (__cerb)
    // 	{
    // 	  ios_base::iostate __err = ios_base::goodbit;
    // 	  __try
    // 	    {
    // 	      const int_type __eof = traits_type::eof();
    // 	      __streambuf_type* __sb = this->rdbuf();
    // 
    // 	      if (traits_type::eq_int_type(__sb->sbumpc(), __eof))
    // 		__err |= ios_base::eofbit;
    // 	      else
    // 		_M_gcount = 1;
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { this->_M_setstate(ios_base::badbit); }
    // 	  if (__err)
    // 	    this->setstate(__err);
    // 	}
    //       return *this;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     basic_istream<_CharT, _Traits>&
    //     basic_istream<_CharT, _Traits>::
    //     ignore(streamsize __n)
    //     {
    //       _M_gcount = 0;
    //       sentry __cerb(*this, true);
    //       if (__cerb && __n > 0)
    //         {
    //           ios_base::iostate __err = ios_base::goodbit;
    //           __try
    //             {
    //               const int_type __eof = traits_type::eof();
    //               __streambuf_type* __sb = this->rdbuf();
    //               int_type __c = __sb->sgetc();
    // 
    // 	      // N.B. On LFS-enabled platforms streamsize is still 32 bits
    // 	      // wide: if we want to implement the standard mandated behavior
    // 	      // for n == max() (see 27.6.1.3/24) we are at risk of signed
    // 	      // integer overflow: thus these contortions. Also note that,
    // 	      // by definition, when more than 2G chars are actually ignored,
    // 	      // _M_gcount (the return value of gcount, that is) cannot be
    // 	      // really correct, being unavoidably too small.
    // 	      bool __large_ignore = false;
    // 	      while (true)
    // 		{
    // 		  while (_M_gcount < __n
    // 			 && !traits_type::eq_int_type(__c, __eof))
    // 		    {
    // 		      ++_M_gcount;
    // 		      __c = __sb->snextc();
    // 		    }
    // 		  if (__n == __gnu_cxx::__numeric_traits<streamsize>::__max
    // 		      && !traits_type::eq_int_type(__c, __eof))
    // 		    {
    // 		      _M_gcount =
    // 			__gnu_cxx::__numeric_traits<streamsize>::__min;
    // 		      __large_ignore = true;
    // 		    }
    // 		  else
    // 		    break;
    // 		}
    // 
    // 	      if (__n == __gnu_cxx::__numeric_traits<streamsize>::__max)
    // 		{
    // 		  if (__large_ignore)
    // 		    _M_gcount = __gnu_cxx::__numeric_traits<streamsize>::__max;
    // 
    // 		  if (traits_type::eq_int_type(__c, __eof))
    // 		    __err |= ios_base::eofbit;
    // 		}
    // 	      else if (_M_gcount < __n)
    // 		{
    // 		  if (traits_type::eq_int_type(__c, __eof))
    // 		    __err |= ios_base::eofbit;
    // 		}
    //             }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    //           __catch(...)
    //             { this->_M_setstate(ios_base::badbit); }
    //           if (__err)
    //             this->setstate(__err);
    //         }
    //       return *this;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     basic_istream<_CharT, _Traits>&
    //     basic_istream<_CharT, _Traits>::
    //     ignore(streamsize __n, int_type __delim)
    //     {
    //       _M_gcount = 0;
    //       sentry __cerb(*this, true);
    //       if (__cerb && __n > 0)
    //         {
    //           ios_base::iostate __err = ios_base::goodbit;
    //           __try
    //             {
    //               const int_type __eof = traits_type::eof();
    //               __streambuf_type* __sb = this->rdbuf();
    //               int_type __c = __sb->sgetc();
    // 
    // 	      // See comment above.
    // 	      bool __large_ignore = false;
    // 	      while (true)
    // 		{
    // 		  while (_M_gcount < __n
    // 			 && !traits_type::eq_int_type(__c, __eof)
    // 			 && !traits_type::eq_int_type(__c, __delim))
    // 		    {
    // 		      ++_M_gcount;
    // 		      __c = __sb->snextc();
    // 		    }
    // 		  if (__n == __gnu_cxx::__numeric_traits<streamsize>::__max
    // 		      && !traits_type::eq_int_type(__c, __eof)
    // 		      && !traits_type::eq_int_type(__c, __delim))
    // 		    {
    // 		      _M_gcount =
    // 			__gnu_cxx::__numeric_traits<streamsize>::__min;
    // 		      __large_ignore = true;
    // 		    }
    // 		  else
    // 		    break;
    // 		}
    // 
    // 	      if (__n == __gnu_cxx::__numeric_traits<streamsize>::__max)
    // 		{
    // 		  if (__large_ignore)
    // 		    _M_gcount = __gnu_cxx::__numeric_traits<streamsize>::__max;
    // 
    // 		  if (traits_type::eq_int_type(__c, __eof))
    // 		    __err |= ios_base::eofbit;
    // 		  else
    // 		    {
    // 		      if (_M_gcount != __n)
    // 			++_M_gcount;
    // 		      __sb->sbumpc();
    // 		    }
    // 		}
    // 	      else if (_M_gcount < __n) // implies __c == __delim or EOF
    // 		{
    // 		  if (traits_type::eq_int_type(__c, __eof))
    // 		    __err |= ios_base::eofbit;
    // 		  else
    // 		    {
    // 		      ++_M_gcount;
    // 		      __sb->sbumpc();
    // 		    }
    // 		}
    //             }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    //           __catch(...)
    //             { this->_M_setstate(ios_base::badbit); }
    //           if (__err)
    //             this->setstate(__err);
    //         }
    //       return *this;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     typename basic_istream<_CharT, _Traits>::int_type
    //     basic_istream<_CharT, _Traits>::
    //     peek(void)
    //     {
    //       int_type __c = traits_type::eof();
    //       _M_gcount = 0;
    //       sentry __cerb(*this, true);
    //       if (__cerb)
    // 	{
    // 	  ios_base::iostate __err = ios_base::goodbit;
    // 	  __try
    // 	    {
    // 	      __c = this->rdbuf()->sgetc();
    // 	      if (traits_type::eq_int_type(__c, traits_type::eof()))
    // 		__err |= ios_base::eofbit;
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { this->_M_setstate(ios_base::badbit); }
    // 	  if (__err)
    // 	    this->setstate(__err);
    // 	}
    //       return __c;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     basic_istream<_CharT, _Traits>&
    //     basic_istream<_CharT, _Traits>::
    //     read(char_type* __s, streamsize __n)
    //     {
    //       _M_gcount = 0;
    //       sentry __cerb(*this, true);
    //       if (__cerb)
    // 	{
    // 	  ios_base::iostate __err = ios_base::goodbit;
    // 	  __try
    // 	    {
    // 	      _M_gcount = this->rdbuf()->sgetn(__s, __n);
    // 	      if (_M_gcount != __n)
    // 		__err |= (ios_base::eofbit | ios_base::failbit);
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { this->_M_setstate(ios_base::badbit); }
    // 	  if (__err)
    // 	    this->setstate(__err);
    // 	}
    //       return *this;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     streamsize
    //     basic_istream<_CharT, _Traits>::
    //     readsome(char_type* __s, streamsize __n)
    //     {
    //       _M_gcount = 0;
    //       sentry __cerb(*this, true);
    //       if (__cerb)
    // 	{
    // 	  ios_base::iostate __err = ios_base::goodbit;
    // 	  __try
    // 	    {
    // 	      // Cannot compare int_type with streamsize generically.
    // 	      const streamsize __num = this->rdbuf()->in_avail();
    // 	      if (__num > 0)
    // 		_M_gcount = this->rdbuf()->sgetn(__s, std::min(__num, __n));
    // 	      else if (__num == -1)
    // 		__err |= ios_base::eofbit;
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { this->_M_setstate(ios_base::badbit); }
    // 	  if (__err)
    // 	    this->setstate(__err);
    // 	}
    //       return _M_gcount;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     basic_istream<_CharT, _Traits>&
    //     basic_istream<_CharT, _Traits>::
    //     putback(char_type __c)
    //     {
    //       // _GLIBCXX_RESOLVE_LIB_DEFECTS
    //       // 60. What is a formatted input function?
    //       _M_gcount = 0;
    //       // Clear eofbit per N3168.
    //       this->clear(this->rdstate() & ~ios_base::eofbit);
    //       sentry __cerb(*this, true);
    //       if (__cerb)
    // 	{
    // 	  ios_base::iostate __err = ios_base::goodbit;
    // 	  __try
    // 	    {
    // 	      const int_type __eof = traits_type::eof();
    // 	      __streambuf_type* __sb = this->rdbuf();
    // 	      if (!__sb
    // 		  || traits_type::eq_int_type(__sb->sputbackc(__c), __eof))
    // 		__err |= ios_base::badbit;
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { this->_M_setstate(ios_base::badbit); }
    // 	  if (__err)
    // 	    this->setstate(__err);
    // 	}
    //       return *this;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     basic_istream<_CharT, _Traits>&
    //     basic_istream<_CharT, _Traits>::
    //     unget(void)
    //     {
    //       // _GLIBCXX_RESOLVE_LIB_DEFECTS
    //       // 60. What is a formatted input function?
    //       _M_gcount = 0;
    //       // Clear eofbit per N3168.
    //       this->clear(this->rdstate() & ~ios_base::eofbit);
    //       sentry __cerb(*this, true);
    //       if (__cerb)
    // 	{
    // 	  ios_base::iostate __err = ios_base::goodbit;
    // 	  __try
    // 	    {
    // 	      const int_type __eof = traits_type::eof();
    // 	      __streambuf_type* __sb = this->rdbuf();
    // 	      if (!__sb
    // 		  || traits_type::eq_int_type(__sb->sungetc(), __eof))
    // 		__err |= ios_base::badbit;
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { this->_M_setstate(ios_base::badbit); }
    // 	  if (__err)
    // 	    this->setstate(__err);
    // 	}
    //       return *this;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     int
    //     basic_istream<_CharT, _Traits>::
    //     sync(void)
    //     {
    //       // _GLIBCXX_RESOLVE_LIB_DEFECTS
    //       // DR60.  Do not change _M_gcount.
    //       int __ret = -1;
    //       sentry __cerb(*this, true);
    //       if (__cerb)
    // 	{
    // 	  ios_base::iostate __err = ios_base::goodbit;
    // 	  __try
    // 	    {
    // 	      __streambuf_type* __sb = this->rdbuf();
    // 	      if (__sb)
    // 		{
    // 		  if (__sb->pubsync() == -1)
    // 		    __err |= ios_base::badbit;
    // 		  else
    // 		    __ret = 0;
    // 		}
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { this->_M_setstate(ios_base::badbit); }
    // 	  if (__err)
    // 	    this->setstate(__err);
    // 	}
    //       return __ret;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     typename basic_istream<_CharT, _Traits>::pos_type
    //     basic_istream<_CharT, _Traits>::
    //     tellg(void)
    //     {
    //       // _GLIBCXX_RESOLVE_LIB_DEFECTS
    //       // DR60.  Do not change _M_gcount.
    //       pos_type __ret = pos_type(-1);
    //       sentry __cerb(*this, true);
    //       if (__cerb)
    // 	{
    // 	  __try
    // 	    {
    // 	      if (!this->fail())
    // 		__ret = this->rdbuf()->pubseekoff(0, ios_base::cur,
    // 						  ios_base::in);
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { this->_M_setstate(ios_base::badbit); }
    // 	}
    //       return __ret;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     basic_istream<_CharT, _Traits>&
    //     basic_istream<_CharT, _Traits>::
    //     seekg(pos_type __pos)
    //     {
    //       // _GLIBCXX_RESOLVE_LIB_DEFECTS
    //       // DR60.  Do not change _M_gcount.
    //       // Clear eofbit per N3168.
    //       this->clear(this->rdstate() & ~ios_base::eofbit);
    //       sentry __cerb(*this, true);
    //       if (__cerb)
    // 	{
    // 	  ios_base::iostate __err = ios_base::goodbit;
    // 	  __try
    // 	    {
    // 	      if (!this->fail())
    // 		{
    // 		  // 136.  seekp, seekg setting wrong streams?
    // 		  const pos_type __p = this->rdbuf()->pubseekpos(__pos,
    // 								 ios_base::in);
    // 		  
    // 		  // 129.  Need error indication from seekp() and seekg()
    // 		  if (__p == pos_type(off_type(-1)))
    // 		    __err |= ios_base::failbit;
    // 		}
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { this->_M_setstate(ios_base::badbit); }
    // 	  if (__err)
    // 	    this->setstate(__err);
    // 	}
    //       return *this;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     basic_istream<_CharT, _Traits>&
    //     basic_istream<_CharT, _Traits>::
    //     seekg(off_type __off, ios_base::seekdir __dir)
    //     {
    //       // _GLIBCXX_RESOLVE_LIB_DEFECTS
    //       // DR60.  Do not change _M_gcount.
    //       // Clear eofbit per N3168.
    //       this->clear(this->rdstate() & ~ios_base::eofbit);
    //       sentry __cerb(*this, true);
    //       if (__cerb)
    // 	{
    // 	  ios_base::iostate __err = ios_base::goodbit;
    // 	  __try
    // 	    {
    // 	      if (!this->fail())
    // 		{
    // 		  // 136.  seekp, seekg setting wrong streams?
    // 		  const pos_type __p = this->rdbuf()->pubseekoff(__off, __dir,
    // 								 ios_base::in);
    // 	      
    // 		  // 129.  Need error indication from seekp() and seekg()
    // 		  if (__p == pos_type(off_type(-1)))
    // 		    __err |= ios_base::failbit;
    // 		}
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      this->_M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { this->_M_setstate(ios_base::badbit); }
    // 	  if (__err)
    // 	    this->setstate(__err);
    // 	}
    //       return *this;
    //     }
    // 
    //   // 27.6.1.2.3 Character extraction templates
    //   template<typename _CharT, typename _Traits>
    //     basic_istream<_CharT, _Traits>&
    //     operator>>(basic_istream<_CharT, _Traits>& __in, _CharT& __c)
    //     {
    //       typedef basic_istream<_CharT, _Traits>		__istream_type;
    //       typedef typename __istream_type::int_type         __int_type;
    // 
    //       typename __istream_type::sentry __cerb(__in, false);
    //       if (__cerb)
    // 	{
    // 	  ios_base::iostate __err = ios_base::goodbit;
    // 	  __try
    // 	    {
    // 	      const __int_type __cb = __in.rdbuf()->sbumpc();
    // 	      if (!_Traits::eq_int_type(__cb, _Traits::eof()))
    // 		__c = _Traits::to_char_type(__cb);
    // 	      else
    // 		__err |= (ios_base::eofbit | ios_base::failbit);
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      __in._M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { __in._M_setstate(ios_base::badbit); }
    // 	  if (__err)
    // 	    __in.setstate(__err);
    // 	}
    //       return __in;
    //     }
    // 
    //   template<typename _CharT, typename _Traits>
    //     void
    //     __istream_extract(basic_istream<_CharT, _Traits>& __in, _CharT* __s,
    // 		      streamsize __num)
    //     {
    //       typedef basic_istream<_CharT, _Traits>		__istream_type;
    //       typedef basic_streambuf<_CharT, _Traits>          __streambuf_type;
    //       typedef typename _Traits::int_type		int_type;
    //       typedef _CharT					char_type;
    //       typedef ctype<_CharT>				__ctype_type;
    // 
    //       streamsize __extracted = 0;
    //       ios_base::iostate __err = ios_base::goodbit;
    //       typename __istream_type::sentry __cerb(__in, false);
    //       if (__cerb)
    // 	{
    // 	  __try
    // 	    {
    // 	      // Figure out how many characters to extract.
    // 	      streamsize __width = __in.width();
    // 	      if (0 < __width && __width < __num)
    // 		__num = __width;
    // 
    // 	      const __ctype_type& __ct = use_facet<__ctype_type>(__in.getloc());
    // 
    // 	      const int_type __eof = _Traits::eof();
    // 	      __streambuf_type* __sb = __in.rdbuf();
    // 	      int_type __c = __sb->sgetc();
    // 
    // 	      while (__extracted < __num - 1
    // 		     && !_Traits::eq_int_type(__c, __eof)
    // 		     && !__ct.is(ctype_base::space,
    // 				 _Traits::to_char_type(__c)))
    // 		{
    // 		  *__s++ = _Traits::to_char_type(__c);
    // 		  ++__extracted;
    // 		  __c = __sb->snextc();
    // 		}
    // 
    // 	      if (__extracted < __num - 1
    // 		  && _Traits::eq_int_type(__c, __eof))
    // 		__err |= ios_base::eofbit;
    // 
    // 	      // _GLIBCXX_RESOLVE_LIB_DEFECTS
    // 	      // 68.  Extractors for char* should store null at end
    // 	      *__s = char_type();
    // 	      __in.width(0);
    // 	    }
    // 	  __catch(__cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      __in._M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch(...)
    // 	    { __in._M_setstate(ios_base::badbit); }
    // 	}
    //       if (!__extracted)
    // 	__err |= ios_base::failbit;
    //       if (__err)
    // 	__in.setstate(__err);
    //     }
    // 
    //   // 27.6.1.4 Standard basic_istream manipulators
    //   template<typename _CharT, typename _Traits>
    //     basic_istream<_CharT, _Traits>&
    //     ws(basic_istream<_CharT, _Traits>& __in)
    //     {
    //       typedef basic_istream<_CharT, _Traits>		__istream_type;
    //       typedef basic_streambuf<_CharT, _Traits>          __streambuf_type;
    //       typedef typename __istream_type::int_type		__int_type;
    //       typedef ctype<_CharT>				__ctype_type;
    // 
    //       // _GLIBCXX_RESOLVE_LIB_DEFECTS
    //       // 451. behavior of std::ws
    //       typename __istream_type::sentry __cerb(__in, true);
    //       if (__cerb)
    // 	{
    // 	  ios_base::iostate __err = ios_base::goodbit;
    // 	  __try
    // 	    {
    // 	      const __ctype_type& __ct = use_facet<__ctype_type>(__in.getloc());
    // 	      const __int_type __eof = _Traits::eof();
    // 	      __streambuf_type* __sb = __in.rdbuf();
    // 	      __int_type __c = __sb->sgetc();
    // 
    // 	      while (true)
    // 		{
    // 		  if (_Traits::eq_int_type(__c, __eof))
    // 		    {
    // 		      __err = ios_base::eofbit;
    // 		      break;
    // 		    }
    // 		  if (!__ct.is(ctype_base::space, _Traits::to_char_type(__c)))
    // 		    break;
    // 		  __c = __sb->snextc();
    // 		}
    // 	    }
    // 	  __catch (const __cxxabiv1::__forced_unwind&)
    // 	    {
    // 	      __in._M_setstate(ios_base::badbit);
    // 	      __throw_exception_again;
    // 	    }
    // 	  __catch (...)
    // 	    {
    // 	      __in._M_setstate(ios_base::badbit);
    // 	    }
    // 	  if (__err)
    // 	    __in.setstate(__err);
    // 	}
    //       return __in;
    //     }
    // 
    //   // Inhibit implicit instantiations for required instantiations,
    //   // which are defined via explicit instantiations elsewhere.
    // #if _GLIBCXX_EXTERN_TEMPLATE
    //   extern template class basic_istream<char>;
    //   extern template istream& ws(istream&);
    //   extern template istream& operator>>(istream&, char&);
    //   extern template istream& operator>>(istream&, unsigned char&);
    //   extern template istream& operator>>(istream&, signed char&);
    // 
    //   extern template istream& istream::_M_extract(unsigned short&);
    //   extern template istream& istream::_M_extract(unsigned int&);  
    //   extern template istream& istream::_M_extract(long&);
    //   extern template istream& istream::_M_extract(unsigned long&);
    //   extern template istream& istream::_M_extract(bool&);
    // #ifdef _GLIBCXX_USE_LONG_LONG
    //   extern template istream& istream::_M_extract(long long&);
    //   extern template istream& istream::_M_extract(unsigned long long&);
    // #endif
    //   extern template istream& istream::_M_extract(float&);
    //   extern template istream& istream::_M_extract(double&);
    //   extern template istream& istream::_M_extract(long double&);
    //   extern template istream& istream::_M_extract(void*&);
    // 
    //   extern template class basic_iostream<char>;
    // 
    // #ifdef _GLIBCXX_USE_WCHAR_T
    //   extern template class basic_istream<wchar_t>;
    //   extern template wistream& ws(wistream&);
    //   extern template wistream& operator>>(wistream&, wchar_t&);
    //   extern template void __istream_extract(wistream&, wchar_t*, streamsize);
    // 
    //   extern template wistream& wistream::_M_extract(unsigned short&);
    //   extern template wistream& wistream::_M_extract(unsigned int&);  
    //   extern template wistream& wistream::_M_extract(long&);
    //   extern template wistream& wistream::_M_extract(unsigned long&);
    //   extern template wistream& wistream::_M_extract(bool&);
    // #ifdef _GLIBCXX_USE_LONG_LONG
    //   extern template wistream& wistream::_M_extract(long long&);
    //   extern template wistream& wistream::_M_extract(unsigned long long&);
    // #endif
    //   extern template wistream& wistream::_M_extract(float&);
    //   extern template wistream& wistream::_M_extract(double&);
    //   extern template wistream& wistream::_M_extract(long double&);
    //   extern template wistream& wistream::_M_extract(void*&);
    // 
    //   extern template class basic_iostream<wchar_t>;
    // #endif
    // #endif
    // 
    // _GLIBCXX_END_NAMESPACE_VERSION
    // } // namespace std
    // 
    // #endif
    // GNU❗✔️: counted C-locale decimal stream extraction; explicit EOF/fail
    // and produced value preserve the source distinction from lexical cast.
    let mut cursor = 0;
    if skip_whitespace {
        while bytes.get(cursor).is_some_and(|b| ascii_space(*b)) {
            cursor += 1;
        }
    }
    let negative = bytes.get(cursor) == Some(&b'-');
    if matches!(bytes.get(cursor), Some(b'+' | b'-')) {
        cursor += 1;
    }
    let mut seen_point = false;
    let mut integer_digits = 0_usize;
    let mut digit_count = 0_usize;
    let mut first_significant = None;
    let mut significant = Vec::with_capacity(769);
    let mut sticky = false;
    while let Some(&byte) = bytes.get(cursor) {
        match byte {
            digit @ b'0'..=b'9' => {
                if !seen_point {
                    integer_digits += 1;
                }
                if first_significant.is_none() && digit != b'0' {
                    first_significant = Some(digit_count);
                }
                if first_significant.is_some() {
                    if significant.len() < 768 {
                        significant.push(digit);
                    } else {
                        sticky |= digit != b'0';
                    }
                }
                digit_count += 1;
            }
            b'.' if !seen_point => seen_point = true,
            _ => break,
        }
        cursor += 1;
    }
    let mut exponent = 0_i128;
    let mut syntax_failure = digit_count == 0;
    if digit_count != 0 && matches!(bytes.get(cursor), Some(b'e' | b'E')) {
        cursor += 1;
        let exponent_negative = bytes.get(cursor) == Some(&b'-');
        if matches!(bytes.get(cursor), Some(b'+' | b'-')) {
            cursor += 1;
        }
        let exponent_start = cursor;
        // Input displacement cannot exceed input length. Saturation beyond
        // that length +4096 cannot cancel into binary64's boundary interval.
        let limit = bytes.len() as i128 + 4096;
        while let Some(&digit) = bytes.get(cursor) {
            if !digit.is_ascii_digit() {
                break;
            }
            exponent = (10 * exponent + i128::from(digit - b'0')).min(limit);
            cursor += 1;
        }
        syntax_failure |= cursor == exponent_start;
        if exponent_negative {
            exponent = -exponent;
        }
    }
    let mut state = SourceNumericStreamState {
        eof: cursor == bytes.len(),
        fail: syntax_failure,
        bad: false,
    };
    let mut failure = syntax_failure.then_some(DoubleLexicalReadErrorKind::Syntax);
    let bits = if syntax_failure {
        0
    } else if let Some(first) = first_significant {
        exponent += integer_digits as i128 - first as i128 - 1;
        if sticky {
            significant.push(b'1');
        } else {
            while significant.len() > 1 && significant.last() == Some(&b'0') {
                significant.pop();
            }
        }
        let bits = decimal_binary64_bits(&significant, exponent, negative, mode);
        if bits & !SIGN_BIT == INFINITY_BITS {
            state.fail = true;
            failure = Some(DoubleLexicalReadErrorKind::Overflow);
            MAX_FINITE_BITS | (bits & SIGN_BIT)
        } else {
            bits
        }
    } else {
        if negative { SIGN_BIT } else { 0 }
    };
    SourceDoubleExtraction {
        value: f64::from_bits(bits),
        consumed: cursor,
        state,
        failure,
    }
}

/// Boost's counted-byte lexical-cast double, with no whitespace skipping.
#[rustfmt::skip]
pub fn source_lexical_double_with_rounding(
    bytes: &[u8],
    mode: SourceNumericRoundingMode,
) -> Result<f64, DoubleLexicalReadError> {
    // BEGIN BOOST185 VERBATIM FUNCTION 0
    /*
    template <typename Target, typename Source>
    inline Target lexical_cast(const Source &arg)
    {
        Target result = Target();

        if (!boost::conversion::detail::try_lexical_convert(arg, result)) {
            boost::conversion::detail::throw_bad_cast<Source, Target>();
        }

        return result;
    }
    */
    // END BOOST185 VERBATIM FUNCTION 0
    // BEGIN BOOST185 VERBATIM FUNCTION 1
    /*
        template <class CharT, class T>
        inline bool parse_inf_nan_impl(const CharT* begin, const CharT* end, T& value
            , const CharT* lc_NAN, const CharT* lc_nan
            , const CharT* lc_INFINITY, const CharT* lc_infinity
            , const CharT opening_brace, const CharT closing_brace) noexcept
        {
            if (begin == end) return false;
            const CharT minus = lcast_char_constants<CharT>::minus;
            const CharT plus = lcast_char_constants<CharT>::plus;
            const int inifinity_size = 8; // == sizeof("infinity") - 1

            /* Parsing +/- */
            bool const has_minus = (*begin == minus);
            if (has_minus || *begin == plus) {
                ++ begin;
            }

            if (end - begin < 3) return false;
            if (lc_iequal(begin, lc_nan, lc_NAN, 3)) {
                begin += 3;
                if (end != begin) {
                    /* It is 'nan(...)' or some bad input*/

                    if (end - begin < 2) return false; // bad input
                    -- end;
                    if (*begin != opening_brace || *end != closing_brace) return false; // bad input
                }

                if( !has_minus ) value = std::numeric_limits<T>::quiet_NaN();
                else value = boost::core::copysign(std::numeric_limits<T>::quiet_NaN(), static_cast<T>(-1));
                return true;
            } else if (
                ( /* 'INF' or 'inf' */
                  end - begin == 3      // 3 == sizeof('inf') - 1
                  && lc_iequal(begin, lc_infinity, lc_INFINITY, 3)
                )
                ||
                ( /* 'INFINITY' or 'infinity' */
                  end - begin == inifinity_size
                  && lc_iequal(begin, lc_infinity, lc_INFINITY, inifinity_size)
                )
             )
            {
                if( !has_minus ) value = std::numeric_limits<T>::infinity();
                else value = -std::numeric_limits<T>::infinity();
                return true;
            }

            return false;
        }
    */
    // END BOOST185 VERBATIM FUNCTION 1
    // BEGIN BOOST185 VERBATIM FUNCTION 2
    /*
            template<typename InputStreamable>
            bool shr_using_base_class(InputStreamable& output)
            {
                static_assert(
                    !boost::is_pointer<InputStreamable>::value,
                    "boost::lexical_cast can not convert to pointers"
                );

    #if defined(BOOST_NO_STRINGSTREAM) || defined(BOOST_NO_STD_LOCALE)
                static_assert(boost::is_same<char, CharT>::value,
                    "boost::lexical_cast can not convert, because your STL library does not "
                    "support such conversions. Try updating it."
                );
    #endif

    #if defined(BOOST_NO_STRINGSTREAM)
                std::istrstream stream(start, static_cast<std::istrstream::streamsize>(finish - start));
    #else
                typedef detail::lcast::buffer_t<CharT, Traits> buffer_t;
                buffer_t buf;
                // Usually `istream` and `basic_istream` do not modify
                // content of buffer; `buffer_t` assures that this is true
                buf.setbuf(const_cast<CharT*>(start), static_cast<typename buffer_t::streamsize>(finish - start));
    #if defined(BOOST_NO_STD_LOCALE)
                std::istream stream(&buf);
    #else
                std::basic_istream<CharT, Traits> stream(&buf);
    #endif // BOOST_NO_STD_LOCALE
    #endif // BOOST_NO_STRINGSTREAM

    #ifndef BOOST_NO_EXCEPTIONS
                stream.exceptions(std::ios::badbit);
                try {
    #endif
                stream.unsetf(std::ios::skipws);
                lcast_set_precision(stream, static_cast<InputStreamable*>(0));

                return (stream >> output)
                    && (stream.get() == Traits::eof());

    #ifndef BOOST_NO_EXCEPTIONS
                } catch (const ::std::ios_base::failure& /*f*/) {
                    return false;
                }
    #endif
            }
        */
    // END BOOST185 VERBATIM FUNCTION 2
    // BEGIN BOOST185 VERBATIM FUNCTION 3
    /*
        template <class T>
        bool float_types_converter_internal(T& output) {
            if (parse_inf_nan(start, finish, output)) return true;
            bool const return_value = shr_using_base_class(output);

            /* Some compilers and libraries successfully
             * parse 'inf', 'INFINITY', '1.0E', '1.0E-'...
             * We are trying to provide a unified behaviour,
             * so we just forbid such conversions (as some
             * of the most popular compilers/libraries do)
             * */
            CharT const minus = lcast_char_constants<CharT>::minus;
            CharT const plus = lcast_char_constants<CharT>::plus;
            CharT const capital_e = lcast_char_constants<CharT>::capital_e;
            CharT const lowercase_e = lcast_char_constants<CharT>::lowercase_e;
            if ( return_value &&
                 (
                    Traits::eq(*(finish-1), lowercase_e)                   // 1.0e
                    || Traits::eq(*(finish-1), capital_e)                  // 1.0E
                    || Traits::eq(*(finish-1), minus)                      // 1.0e- or 1.0E-
                    || Traits::eq(*(finish-1), plus)                       // 1.0e+ or 1.0E+
                 )
            ) return false;

            return return_value;
        }
    */
    // END BOOST185 VERBATIM FUNCTION 3
    // BEGIN BOOST185 VERBATIM FUNCTION 4
    /*
        template <class CharT, class T>
        bool parse_inf_nan(const CharT* begin, const CharT* end, T& value) noexcept {
            return parse_inf_nan_impl(begin, end, value
                               , "NAN", "nan"
                               , "INFINITY", "infinity"
                               , '(', ')');
        }
    */
    // END BOOST185 VERBATIM FUNCTION 4
    // BEGIN BOOST185 VERBATIM FUNCTION 5
    /*
        template <class CharT>
        bool lc_iequal(const CharT* val, const CharT* lcase, const CharT* ucase, unsigned int len) noexcept {
            for( unsigned int i=0; i < len; ++i ) {
                if ( val[i] != lcase[i] && val[i] != ucase[i] ) return false;
            }

            return true;
        }
    */
    // END BOOST185 VERBATIM FUNCTION 5
    // Boost❗✔️: special literals precede formatted extraction and retain
    // sign without interpreting NaN payload bytes. The extraction/get/EOF
    // sequence retains error value/state rather than replacing it by default.
    if let Some(value) = source_special_double_literal(bytes) { return Ok(value); }
    source_lexical_decimal_with_rounding(bytes, mode)
}

fn source_lexical_decimal_with_rounding(
    bytes: &[u8],
    mode: SourceNumericRoundingMode,
) -> Result<f64, DoubleLexicalReadError> {
    let mut extraction = source_extract_double(bytes, false, mode);
    if extraction.state.fail {
        return Err(DoubleLexicalReadError {
            kind: extraction
                .failure
                .expect("failed extraction carries its source cause"),
            consumed: extraction.consumed,
            state: extraction.state,
            value_bits: Some(extraction.value.to_bits()),
        });
    }
    if extraction.consumed != bytes.len() {
        // stream.get() consumes exactly one extra byte after extraction.
        extraction.consumed += 1;
        return Err(DoubleLexicalReadError {
            kind: DoubleLexicalReadErrorKind::TrailingByte,
            consumed: extraction.consumed,
            state: extraction.state,
            value_bits: Some(extraction.value.to_bits()),
        });
    }
    Ok(extraction.value)
}

#[rustfmt::skip]
fn source_special_double_literal(bytes: &[u8]) -> Option<f64> {
    // BEGIN BOOST185 VERBATIM FUNCTION 1
    /*
        template <class CharT, class T>
        inline bool parse_inf_nan_impl(const CharT* begin, const CharT* end, T& value
            , const CharT* lc_NAN, const CharT* lc_nan
            , const CharT* lc_INFINITY, const CharT* lc_infinity
            , const CharT opening_brace, const CharT closing_brace) noexcept
        {
            if (begin == end) return false;
            const CharT minus = lcast_char_constants<CharT>::minus;
            const CharT plus = lcast_char_constants<CharT>::plus;
            const int inifinity_size = 8; // == sizeof("infinity") - 1

            /* Parsing +/- */
            bool const has_minus = (*begin == minus);
            if (has_minus || *begin == plus) {
                ++ begin;
            }

            if (end - begin < 3) return false;
            if (lc_iequal(begin, lc_nan, lc_NAN, 3)) {
                begin += 3;
                if (end != begin) {
                    /* It is 'nan(...)' or some bad input*/

                    if (end - begin < 2) return false; // bad input
                    -- end;
                    if (*begin != opening_brace || *end != closing_brace) return false; // bad input
                }

                if( !has_minus ) value = std::numeric_limits<T>::quiet_NaN();
                else value = boost::core::copysign(std::numeric_limits<T>::quiet_NaN(), static_cast<T>(-1));
                return true;
            } else if (
                ( /* 'INF' or 'inf' */
                  end - begin == 3      // 3 == sizeof('inf') - 1
                  && lc_iequal(begin, lc_infinity, lc_INFINITY, 3)
                )
                ||
                ( /* 'INFINITY' or 'infinity' */
                  end - begin == inifinity_size
                  && lc_iequal(begin, lc_infinity, lc_INFINITY, inifinity_size)
                )
             )
            {
                if( !has_minus ) value = std::numeric_limits<T>::infinity();
                else value = -std::numeric_limits<T>::infinity();
                return true;
            }

            return false;
        }
    */
    // END BOOST185 VERBATIM FUNCTION 1
    let negative = bytes.first() == Some(&b'-');
    let body = if matches!(bytes.first(), Some(b'+' | b'-')) {
        &bytes[1..]
    } else {
        bytes
    };
    if body.eq_ignore_ascii_case(b"inf") || body.eq_ignore_ascii_case(b"infinity") {
        return Some(f64::from_bits(
            INFINITY_BITS | if negative { SIGN_BIT } else { 0 },
        ));
    }
    if body.len() >= 3
        && body[..3].eq_ignore_ascii_case(b"nan")
        && (body.len() == 3 || (body.len() >= 5 && body[3] == b'(' && body.last() == Some(&b')')))
    {
        return Some(f64::from_bits(
            0x7ff8_0000_0000_0000 | if negative { SIGN_BIT } else { 0 },
        ));
    }

    None
}

/// Observe the source floating-point rounding state. No assumed native default.
#[rustfmt::skip]
fn source_current_rounding_mode() -> Result<SourceNumericRoundingMode, DoubleLexicalReadError> {
    // /* Determine floating-point rounding mode within libc.  Generic version.
    //    Copyright (C) 2012-2026 Free Software Foundation, Inc.
    //    This file is part of the GNU C Library.
    // 
    //    The GNU C Library is free software; you can redistribute it and/or
    //    modify it under the terms of the GNU Lesser General Public
    //    License as published by the Free Software Foundation; either
    //    version 2.1 of the License, or (at your option) any later version.
    // 
    //    The GNU C Library is distributed in the hope that it will be useful,
    //    but WITHOUT ANY WARRANTY; without even the implied warranty of
    //    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    //    Lesser General Public License for more details.
    // 
    //    You should have received a copy of the GNU Lesser General Public
    //    License along with the GNU C Library; if not, see
    //    <https://www.gnu.org/licenses/>.  */
    // 
    // #ifndef _GET_ROUNDING_MODE_H
    // #define _GET_ROUNDING_MODE_H	1
    // 
    // #include <fpu_control.h>
    // #include <stdlib.h>
    // 
    // /* Define values for FE_* modes not defined for this architecture.  */
    // #ifdef FE_DOWNWARD
    // # define ORIG_FE_DOWNWARD FE_DOWNWARD
    // #else
    // # define ORIG_FE_DOWNWARD 0
    // #endif
    // #ifdef FE_TONEAREST
    // # define ORIG_FE_TONEAREST FE_TONEAREST
    // #else
    // # define ORIG_FE_TONEAREST 0
    // #endif
    // #ifdef FE_TOWARDZERO
    // # define ORIG_FE_TOWARDZERO FE_TOWARDZERO
    // #else
    // # define ORIG_FE_TOWARDZERO 0
    // #endif
    // #ifdef FE_UPWARD
    // # define ORIG_FE_UPWARD FE_UPWARD
    // #else
    // # define ORIG_FE_UPWARD 0
    // #endif
    // #define FE_CONSTRUCT_DISTINCT_VALUE(X, Y, Z) \
    //   ((((X) & 1) | ((Y) & 2) | ((Z) & 4)) ^ 7)
    // #ifndef FE_DOWNWARD
    // # define FE_DOWNWARD FE_CONSTRUCT_DISTINCT_VALUE (ORIG_FE_TONEAREST,	\
    // 						  ORIG_FE_TOWARDZERO,	\
    // 						  ORIG_FE_UPWARD)
    // #endif
    // #ifndef FE_TONEAREST
    // # define FE_TONEAREST FE_CONSTRUCT_DISTINCT_VALUE (FE_DOWNWARD,		\
    // 						   ORIG_FE_TOWARDZERO,	\
    // 						   ORIG_FE_UPWARD)
    // #endif
    // #ifndef FE_TOWARDZERO
    // # define FE_TOWARDZERO FE_CONSTRUCT_DISTINCT_VALUE (FE_DOWNWARD,	\
    // 						    FE_TONEAREST,	\
    // 						    ORIG_FE_UPWARD)
    // #endif
    // #ifndef FE_UPWARD
    // # define FE_UPWARD FE_CONSTRUCT_DISTINCT_VALUE (FE_DOWNWARD,	\
    // 						FE_TONEAREST,	\
    // 						FE_TOWARDZERO)
    // #endif
    // 
    // /* Return the floating-point rounding mode.  */
    // 
    // static inline int
    // get_rounding_mode (void)
    // {
    // #if (defined _FPU_RC_DOWN			\
    //      || defined _FPU_RC_NEAREST			\
    //      || defined _FPU_RC_ZERO			\
    //      || defined _FPU_RC_UP)
    //   fpu_control_t fc;
    //   const fpu_control_t mask = (0
    // # ifdef _FPU_RC_DOWN
    // 			      | _FPU_RC_DOWN
    // # endif
    // # ifdef _FPU_RC_NEAREST
    // 			      | _FPU_RC_NEAREST
    // # endif
    // # ifdef _FPU_RC_ZERO
    // 			      | _FPU_RC_ZERO
    // # endif
    // # ifdef _FPU_RC_UP
    // 			      | _FPU_RC_UP
    // # endif
    // 			      );
    // 
    //   _FPU_GETCW (fc);
    //   switch (fc & mask)
    //     {
    // # ifdef _FPU_RC_DOWN
    //     case _FPU_RC_DOWN:
    //       return FE_DOWNWARD;
    // # endif
    // 
    // # ifdef _FPU_RC_NEAREST
    //     case _FPU_RC_NEAREST:
    //       return FE_TONEAREST;
    // # endif
    // 
    // # ifdef _FPU_RC_ZERO
    //     case _FPU_RC_ZERO:
    //       return FE_TOWARDZERO;
    // # endif
    // 
    // # ifdef _FPU_RC_UP
    //     case _FPU_RC_UP:
    //       return FE_UPWARD;
    // # endif
    // 
    //     default:
    //       abort ();
    //     }
    // #else
    //   return FE_TONEAREST;
    // #endif
    // }
    // 
    // #endif /* get-rounding-mode.h */
    // /* FPU control word bits.  x86 version.
    //    Copyright (C) 1993-2026 Free Software Foundation, Inc.
    //    This file is part of the GNU C Library.
    // 
    //    The GNU C Library is free software; you can redistribute it and/or
    //    modify it under the terms of the GNU Lesser General Public
    //    License as published by the Free Software Foundation; either
    //    version 2.1 of the License, or (at your option) any later version.
    // 
    //    The GNU C Library is distributed in the hope that it will be useful,
    //    but WITHOUT ANY WARRANTY; without even the implied warranty of
    //    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    //    Lesser General Public License for more details.
    // 
    //    You should have received a copy of the GNU Lesser General Public
    //    License along with the GNU C Library; if not, see
    //    <https://www.gnu.org/licenses/>.  */
    // 
    // #ifndef _FPU_CONTROL_H
    // #define _FPU_CONTROL_H	1
    // 
    // /* Note that this file sets on x86-64 only the x87 FPU, it does not
    //    touch the SSE unit.  */
    // 
    // /* Here is the dirty part. Set up your 387 through the control word
    //  * (cw) register.
    //  *
    //  *     15-13    12  11-10  9-8     7-6     5    4    3    2    1    0
    //  * | reserved | IC | RC  | PC | reserved | PM | UM | OM | ZM | DM | IM
    //  *
    //  * IM: Invalid operation mask
    //  * DM: Denormalized operand mask
    //  * ZM: Zero-divide mask
    //  * OM: Overflow mask
    //  * UM: Underflow mask
    //  * PM: Precision (inexact result) mask
    //  *
    //  * Mask bit is 1 means no interrupt.
    //  *
    //  * PC: Precision control
    //  * 11 - round to extended precision
    //  * 10 - round to double precision
    //  * 00 - round to single precision
    //  *
    //  * RC: Rounding control
    //  * 00 - rounding to nearest
    //  * 01 - rounding down (toward - infinity)
    //  * 10 - rounding up (toward + infinity)
    //  * 11 - rounding toward zero
    //  *
    //  * IC: Infinity control
    //  * That is for 8087 and 80287 only.
    //  *
    //  * The hardware default is 0x037f which we use.
    //  */
    // 
    // #include <features.h>
    // 
    // /* masking of interrupts */
    // #define _FPU_MASK_IM  0x01
    // #define _FPU_MASK_DM  0x02
    // #define _FPU_MASK_ZM  0x04
    // #define _FPU_MASK_OM  0x08
    // #define _FPU_MASK_UM  0x10
    // #define _FPU_MASK_PM  0x20
    // 
    // /* precision control */
    // #define _FPU_EXTENDED 0x300	/* libm requires double extended precision.  */
    // #define _FPU_DOUBLE   0x200
    // #define _FPU_SINGLE   0x0
    // 
    // /* rounding control */
    // #define _FPU_RC_NEAREST 0x0    /* RECOMMENDED */
    // #define _FPU_RC_DOWN    0x400
    // #define _FPU_RC_UP      0x800
    // #define _FPU_RC_ZERO    0xC00
    // 
    // #define _FPU_RESERVED 0xF0C0  /* Reserved bits in cw */
    // 
    // 
    // /* The fdlibm code requires strict IEEE double precision arithmetic,
    //    and no interrupts for exceptions, rounding to nearest.  */
    // 
    // #define _FPU_DEFAULT  0x037f
    // 
    // /* IEEE:  same as above.  */
    // #define _FPU_IEEE     0x037f
    // 
    // /* Type of the control word.  */
    // typedef unsigned int fpu_control_t __attribute__ ((__mode__ (__HI__)));
    // 
    // /* Macros for accessing the hardware control word.  "*&" is used to
    //    work around a bug in older versions of GCC.  __volatile__ is used
    //    to support combination of writing the control register and reading
    //    it back.  Without __volatile__, the old value may be used for reading
    //    back under compiler optimization.
    // 
    //    Note that the use of these macros is not sufficient anymore with
    //    recent hardware nor on x86-64.  Some floating point operations are
    //    executed in the SSE/SSE2 engines which have their own control and
    //    status register.  */
    // #define _FPU_GETCW(cw) __asm__ __volatile__ ("fnstcw %0" : "=m" (*&cw))
    // #define _FPU_SETCW(cw) __asm__ __volatile__ ("fldcw %0" : : "m" (*&cw))
    // 
    // /* Default control word set at startup.  */
    // extern fpu_control_t __fpu_control;
    // 
    // #endif	/* fpu_control.h */
    // /* Return current rounding direction.
    //    Copyright (C) 1997-2026 Free Software Foundation, Inc.
    //    This file is part of the GNU C Library.
    // 
    //    The GNU C Library is free software; you can redistribute it and/or
    //    modify it under the terms of the GNU Lesser General Public
    //    License as published by the Free Software Foundation; either
    //    version 2.1 of the License, or (at your option) any later version.
    // 
    //    The GNU C Library is distributed in the hope that it will be useful,
    //    but WITHOUT ANY WARRANTY; without even the implied warranty of
    //    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    //    Lesser General Public License for more details.
    // 
    //    You should have received a copy of the GNU Lesser General Public
    //    License along with the GNU C Library; if not, see
    //    <https://www.gnu.org/licenses/>.  */
    // 
    // #include <fenv.h>
    // 
    // int
    // __fegetround (void)
    // {
    //   int cw;
    //   /* We only check the x87 FPU unit.  The SSE unit should be the same
    //      - and if it's not the same there's no way to signal it.  */
    // 
    //   __asm__ ("fnstcw %0" : "=m" (cw));
    // 
    //   return cw & 0xc00;
    // }
    // libm_hidden_def (__fegetround)
    // static_weak_alias (__fegetround, fegetround)
    // libm_hidden_weak (fegetround)
    // glibc❗✔️: native x86 source queries x87, including when MXCSR differs.
    // Reading the control register leaves state and exception flags unchanged.
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        let mut control = 0_u16;
        // SAFETY: fnstcw writes one aligned, initialized two-byte stack object;
        // it does not change the x87 environment, touch unowned memory, or call
        // foreign code. Match the source get_rounding_mode/_FPU_GETCW primitive.
        unsafe {
            std::arch::asm!("fnstcw [{destination}]", destination = in(reg) &mut control, options(nostack, preserves_flags));
        }
        return Ok(match (control >> 10) & 3 {
            0 => SourceNumericRoundingMode::NearestEven,
            1 => SourceNumericRoundingMode::Downward,
            2 => SourceNumericRoundingMode::Upward,
            _ => SourceNumericRoundingMode::TowardZero,
        });
    }
    // /* Determine floating-point rounding mode within libc.  AArch64 version.
    // 
    //    Copyright (C) 2012-2026 Free Software Foundation, Inc.
    // 
    //    This file is part of the GNU C Library.
    // 
    //    The GNU C Library is free software; you can redistribute it and/or
    //    modify it under the terms of the GNU Lesser General Public
    //    License as published by the Free Software Foundation; either
    //    version 2.1 of the License, or (at your option) any later version.
    // 
    //    The GNU C Library is distributed in the hope that it will be useful,
    //    but WITHOUT ANY WARRANTY; without even the implied warranty of
    //    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    //    Lesser General Public License for more details.
    // 
    //    You should have received a copy of the GNU Lesser General Public
    //    License along with the GNU C Library; if not, see
    //    <https://www.gnu.org/licenses/>.  */
    // 
    // #ifndef _AARCH64_GET_ROUNDING_MODE_H
    // #define _AARCH64_GET_ROUNDING_MODE_H	1
    // 
    // #include <fenv.h>
    // #include <fpu_control.h>
    // 
    // /* Return the floating-point rounding mode.  */
    // 
    // static inline int
    // get_rounding_mode (void)
    // {
    //   fpu_control_t fpcr;
    // 
    //   _FPU_GETCW (fpcr);
    //   return fpcr & _FPU_FPCR_RM_MASK;
    // }
    // 
    // #endif /* get-rounding-mode.h */
    // /* Copyright (C) 1996-2026 Free Software Foundation, Inc.
    // 
    //    This file is part of the GNU C Library.
    // 
    //    The GNU C Library is free software; you can redistribute it and/or
    //    modify it under the terms of the GNU Lesser General Public License as
    //    published by the Free Software Foundation; either version 2.1 of the
    //    License, or (at your option) any later version.
    // 
    //    The GNU C Library is distributed in the hope that it will be useful,
    //    but WITHOUT ANY WARRANTY; without even the implied warranty of
    //    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    //    Lesser General Public License for more details.
    // 
    //    You should have received a copy of the GNU Lesser General Public
    //    License along with the GNU C Library; if not, see
    //    <https://www.gnu.org/licenses/>.  */
    // 
    // #ifndef _AARCH64_FPU_CONTROL_H
    // #define _AARCH64_FPU_CONTROL_H
    // 
    // #include <features.h>
    // #include <sys/types.h>
    // 
    // /* Macros for accessing the FPCR and FPSR.  */
    // 
    // #if __GNUC_PREREQ (6,0)
    // # define _FPU_GETCW(fpcr) (fpcr = __builtin_aarch64_get_fpcr ())
    // # define _FPU_SETCW(fpcr) __builtin_aarch64_set_fpcr (fpcr)
    // # define _FPU_GETFPSR(fpsr) (fpsr = __builtin_aarch64_get_fpsr ())
    // # define _FPU_SETFPSR(fpsr) __builtin_aarch64_set_fpsr (fpsr)
    // #else
    // # define _FPU_GETCW(fpcr)					\
    //   ({ 								\
    //    __uint64_t __fpcr;						\
    //    __asm__ __volatile__ ("mrs	%0, fpcr" : "=r" (__fpcr));	\
    //    fpcr = __fpcr;						\
    //   })
    // 
    // # define _FPU_SETCW(fpcr)					\
    //   ({								\
    //    __uint64_t __fpcr = fpcr;					\
    //    __asm__ __volatile__ ("msr	fpcr, %0" : : "r" (__fpcr));    \
    //   })
    // 
    // # define _FPU_GETFPSR(fpsr)					\
    //   ({								\
    //    __uint64_t __fpsr;						\
    //    __asm__ __volatile__ ("mrs	%0, fpsr" : "=r" (__fpsr));	\
    //    fpsr = __fpsr;						\
    //   })
    // 
    // # define _FPU_SETFPSR(fpsr)					\
    //   ({								\
    //    __uint64_t __fpsr = fpsr;					\
    //    __asm__ __volatile__ ("msr	fpsr, %0" : : "r" (__fpsr));    \
    //   })
    // #endif
    // 
    // /* Reserved bits should be preserved when modifying register
    //    contents. These two masks indicate which bits in each of FPCR and
    //    FPSR should not be changed.  */
    // 
    // #define _FPU_RESERVED		0xfe0fe0f8
    // #define _FPU_FPSR_RESERVED	0x0fffffe0
    // 
    // #define _FPU_DEFAULT		0x00000000
    // #define _FPU_FPSR_DEFAULT	0x00000000
    // 
    // /* Layout of FPCR and FPSR:
    // 
    //    |       |       |       |       |       |       |       |
    //    0 0 0 0 1 1 1 0 0 0 0 0 1 0 0 0 1 1 1 0 0 0 0 0 1 1 1 0 0 0 0 0
    //    s s s s s                                       s     s s s s s
    //              c c c c c c c               c c c c c
    //    N Z C V Q A D F R R S S S L L L I U U I U O D I I U U I U O D I
    //            C H N Z M M T T B E E E D N N X F F Z O D N N X F F Z O
    //              P     O O R R Z N N N E K K E E E E E C K K C C C C C
    //                    D D I I P
    //                    E E D D
    //                        E E
    //  */
    // 
    // #define _FPU_FPCR_RM_MASK  0xc00000
    // 
    // #define _FPU_FPCR_MASK_IXE 0x1000
    // #define _FPU_FPCR_MASK_UFE 0x0800
    // #define _FPU_FPCR_MASK_OFE 0x0400
    // #define _FPU_FPCR_MASK_DZE 0x0200
    // #define _FPU_FPCR_MASK_IOE 0x0100
    // 
    // #define _FPU_FPCR_IEEE                       \
    //   (_FPU_DEFAULT  | _FPU_FPCR_MASK_IXE	     \
    //    | _FPU_FPCR_MASK_UFE | _FPU_FPCR_MASK_OFE \
    //    | _FPU_FPCR_MASK_DZE | _FPU_FPCR_MASK_IOE)
    // 
    // #define _FPU_FPSR_IEEE 0
    // 
    // typedef unsigned int fpu_control_t;
    // typedef unsigned int fpu_fpsr_t;
    // 
    // /* Default control word set at startup.  */
    // extern fpu_control_t __fpu_control;
    // 
    // #endif
    // glibc❗✔️: AArch64 source observes the two FPCR rounding bits.
    #[cfg(target_arch = "aarch64")]
    {
        let control: u64;
        // SAFETY: mrs only observes the thread-local FPCR; it does not alter
        // rounding, flags, memory, or floating-point register contents.
        unsafe {
            std::arch::asm!("mrs {control}, fpcr", control = out(reg) control, options(nomem, nostack, preserves_flags));
        }
        return Ok(match (control >> 22) & 3 {
            0 => SourceNumericRoundingMode::NearestEven,
            1 => SourceNumericRoundingMode::Upward,
            2 => SourceNumericRoundingMode::Downward,
            _ => SourceNumericRoundingMode::TowardZero,
        });
    }
    // WebAssembly numerics specifies round-to-nearest-ties-to-even and has no
    // mutable hardware rounding mode. This is a target numerical state,
    // not a claim about an unidentified historical native GNU environment.
    // https://webassembly.github.io/spec/core/exec/numerics.html#floating-point-operations
    #[cfg(any(target_arch = "wasm32", target_arch = "wasm64"))]
    {
        return Ok(SourceNumericRoundingMode::NearestEven);
    }
    // glibc❌❌: other native hardware environment access is not modeled here.
    // Do not substitute a conventional nearest mode for an unknown target.
    #[cfg(not(any(
        target_arch = "x86",
        target_arch = "x86_64",
        target_arch = "aarch64",
        target_arch = "wasm32",
        target_arch = "wasm64"
    )))]
    {
        Err(DoubleLexicalReadError {
            kind: DoubleLexicalReadErrorKind::FloatingEnvironmentNotModeled,
            consumed: 0,
            state: SourceNumericStreamState::default(),
            value_bits: None,
        })
    }
}

/// Canonical counted-byte C-locale conversion using current source mode.
pub fn source_lexical_double(bytes: &[u8]) -> Result<f64, DoubleLexicalReadError> {
    if let Some(value) = source_special_double_literal(bytes) {
        return Ok(value);
    }
    source_lexical_decimal_with_rounding(bytes, source_current_rounding_mode()?)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn bits(s: &[u8], mode: SourceNumericRoundingMode) -> u64 {
        source_lexical_double_with_rounding(s, mode)
            .unwrap()
            .to_bits()
    }
    #[test]
    fn complete_literals_sign_payload_and_counted_bytes() {
        for (s, b) in [
            (b"-INFINITY".as_slice(), 0xfff0_0000_0000_0000),
            (b"+inf".as_slice(), INFINITY_BITS),
            (b"NaN(a\0b)".as_slice(), 0x7ff8_0000_0000_0000),
            (b"-nan(whatever)".as_slice(), 0xfff8_0000_0000_0000),
        ] {
            assert_eq!(bits(s, SourceNumericRoundingMode::NearestEven), b);
        }
        for s in [
            b"nanx".as_slice(),
            b"infx",
            b" nan",
            b"1\0",
            b"1 ",
            b"0x1p0",
            b"1e",
            b"1e+",
            b".",
            b"+",
        ] {
            assert!(
                source_lexical_double_with_rounding(s, SourceNumericRoundingMode::NearestEven)
                    .is_err(),
                "{s:?}"
            );
        }
    }
    #[test]
    fn complete_decimal_general_length_cancellation() {
        let mut s = vec![b'1'];
        s.extend(std::iter::repeat_n(b'0', 700000));
        s.extend_from_slice(b"e-700000");
        assert_eq!(
            bits(&s, SourceNumericRoundingMode::NearestEven),
            1.0_f64.to_bits()
        );
        let mut s = b"0.".to_vec();
        s.extend(std::iter::repeat_n(b'0', 700000));
        s.extend_from_slice(b"1e700001");
        assert_eq!(
            bits(&s, SourceNumericRoundingMode::NearestEven),
            1.0_f64.to_bits()
        );
        let mut s = b"-0e+".to_vec();
        s.extend(std::iter::repeat_n(b'9', 10000));
        assert_eq!(bits(&s, SourceNumericRoundingMode::Upward), SIGN_BIT);
    }
    #[test]
    fn complete_all_rounding_modes_and_signs() {
        let half = b"1.00000000000000011102230246251565404236316680908203125";
        for (mode, p, n) in [
            (
                SourceNumericRoundingMode::NearestEven,
                0x3ff0_0000_0000_0000,
                0xbff0_0000_0000_0000,
            ),
            (
                SourceNumericRoundingMode::TowardZero,
                0x3ff0_0000_0000_0000,
                0xbff0_0000_0000_0000,
            ),
            (
                SourceNumericRoundingMode::Upward,
                0x3ff0_0000_0000_0001,
                0xbff0_0000_0000_0000,
            ),
            (
                SourceNumericRoundingMode::Downward,
                0x3ff0_0000_0000_0000,
                0xbff0_0000_0000_0001,
            ),
        ] {
            assert_eq!(bits(half, mode), p);
            let mut minus = vec![b'-'];
            minus.extend_from_slice(half);
            assert_eq!(bits(&minus, mode), n);
        }
        for mode in [
            SourceNumericRoundingMode::NearestEven,
            SourceNumericRoundingMode::TowardZero,
            SourceNumericRoundingMode::Upward,
            SourceNumericRoundingMode::Downward,
        ] {
            assert_eq!(bits(b"2", mode), 2.0_f64.to_bits());
            assert_eq!(
                bits(b"1e-400", mode),
                u64::from(mode == SourceNumericRoundingMode::Upward)
            );
            assert_eq!(
                bits(b"-1e-400", mode),
                SIGN_BIT | u64::from(mode == SourceNumericRoundingMode::Downward)
            );
            let positive = source_lexical_double_with_rounding(b"1e400", mode);
            assert_eq!(
                positive.is_err(),
                matches!(
                    mode,
                    SourceNumericRoundingMode::NearestEven | SourceNumericRoundingMode::Upward
                )
            );
            let negative = source_lexical_double_with_rounding(b"-1e400", mode);
            assert_eq!(
                negative.is_err(),
                matches!(
                    mode,
                    SourceNumericRoundingMode::NearestEven | SourceNumericRoundingMode::Downward
                )
            );
        }
    }
    #[test]
    fn complete_stream_failure_value_cursor_and_eof_state() {
        let mode = SourceNumericRoundingMode::NearestEven;
        let r = source_extract_double(b"  1.5 tail", true, mode);
        assert_eq!(r.value, 1.5);
        assert_eq!(r.consumed, 5);
        assert_eq!(r.state, SourceNumericStreamState::default());
        let r = source_extract_double(b"1e+ tail", false, mode);
        assert_eq!(r.value.to_bits(), 0);
        assert_eq!(r.consumed, 3);
        assert!(r.state.fail);
        assert!(!r.state.eof);
        let r = source_extract_double(b"1e309", false, mode);
        assert_eq!(r.value.to_bits(), MAX_FINITE_BITS);
        assert_eq!(
            r.state,
            SourceNumericStreamState {
                eof: true,
                fail: true,
                bad: false
            }
        );
        let r = source_extract_double(b"1e-400", false, mode);
        assert_eq!(r.value.to_bits(), 0);
        assert_eq!(
            r.state,
            SourceNumericStreamState {
                eof: true,
                fail: false,
                bad: false
            }
        );
        let e = source_lexical_double_with_rounding(b"12abc", mode).unwrap_err();
        assert_eq!(e.consumed, 3);
        assert_eq!(e.value_bits, Some(12.0_f64.to_bits()));
        assert_eq!(e.kind, DoubleLexicalReadErrorKind::TrailingByte);
    }
    #[test]
    fn complete_subnormal_and_maximum_midpoints() {
        let mode = SourceNumericRoundingMode::NearestEven;
        assert_eq!(bits(b"0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000024703282292062327208828439643411068618252990130716238221279284125033775363510437593264991818081799618989828234772285886546332835517796989819938739800539093906315035659515570226392290858392449105184435931802849936536152500319370457678249219365623669863658480757001585769269903706311928279558551332927834338409351978015531246597263579574622766465272827220056374006485499977096599470454020828166226237857393450736339007967761930577506740176324673600968951340535537458516661134223766678604162159680461914467291840300530057530849048765391711386591646239524912623653881879636239373280423891018672348497668235089863388587925628302755995657524455507255189313690836254779186948667994968324049705821028513185451396213837722826145437693412532098591327667236328125",mode),0);
        assert_eq!(bits(b"0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000247032822920623272088284396434110686182529901307162382212792841250337753635104375932649918180817996189898282347722858865463328355177969898199387398005390939063150356595155702263922908583924491051844359318028499365361525003193704576782492193656236698636584807570015857692699037063119282795585513329278343384093519780155312465972635795746227664652728272200563740064854999770965994704540208281662262378573934507363390079677619305775067401763246736009689513405355374585166611342237666786041621596804619144672918403005300575308490487653917113865916462395249126236538818796362393732804238910186723484976682350898633885879256283027559956575244555072551893136908362547791869486679949683240497058210285131854513962138377228261454376934125320985913276672363281251",mode),1);
        assert_eq!(bits(b"0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000247032822920623272088284396434110686182529901307162382212792841250337753635104375932649918180817996189898282347722858865463328355177969898199387398005390939063150356595155702263922908583924491051844359318028499365361525003193704576782492193656236698636584807570015857692699037063119282795585513329278343384093519780155312465972635795746227664652728272200563740064854999770965994704540208281662262378573934507363390079677619305775067401763246736009689513405355374585166611342237666786041621596804619144672918403005300575308490487653917113865916462395249126236538818796362393732804238910186723484976682350898633885879256283027559956575244555072551893136908362547791869486679949683240497058210285131854513962138377228261454376934125320985913276672363281249",mode),0);
        assert!(source_lexical_double_with_rounding(b"179769313486231580793728971405303415079934132710037826936173778980444968292764750946649017977587207096330286416692887910946555547851940402630657488671505820681908902000708383676273854845817711531764475730270069855571366959622842914819860834936475292719074168444365510704342711559699508093042880177904174497792",mode).is_err());
        assert_eq!(bits(b"179769313486231580793728971405303415079934132710037826936173778980444968292764750946649017977587207096330286416692887910946555547851940402630657488671505820681908902000708383676273854845817711531764475730270069855571366959622842914819860834936475292719074168444365510704342711559699508093042880177904174497791",mode),MAX_FINITE_BITS);
        assert_eq!(bits(b"0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000024703282292062327208828439643411068618252990130716238221279284125033775363510437593264991818081799618989828234772285886546332835517796989819938739800539093906315035659515570226392290858392449105184435931802849936536152500319370457678249219365623669863658480757001585769269903706311928279558551332927834338409351978015531246597263579574622766465272827220056374006485499977096599470454020828166226237857393450736339007967761930577506740176324673600968951340535537458516661134223766678604162159680461914467291840300530057530849048765391711386591646239524912623653881879636239373280423891018672348497668235089863388587925628302755995657524455507255189313690836254779186948667994968324049705821028513185451396213837722826145437693412532098591327667236328125",SourceNumericRoundingMode::Upward),1);
        assert_eq!(bits(b"0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000024703282292062327208828439643411068618252990130716238221279284125033775363510437593264991818081799618989828234772285886546332835517796989819938739800539093906315035659515570226392290858392449105184435931802849936536152500319370457678249219365623669863658480757001585769269903706311928279558551332927834338409351978015531246597263579574622766465272827220056374006485499977096599470454020828166226237857393450736339007967761930577506740176324673600968951340535537458516661134223766678604162159680461914467291840300530057530849048765391711386591646239524912623653881879636239373280423891018672348497668235089863388587925628302755995657524455507255189313690836254779186948667994968324049705821028513185451396213837722826145437693412532098591327667236328125",SourceNumericRoundingMode::TowardZero),0);
    }
    #[test]
    fn complete_768_digit_sticky_boundary() {
        let midpoint = b"1.00000000000000011102230246251565404236316680908203125";
        let mut above = midpoint.to_vec();
        above.extend(std::iter::repeat_n(b'0', 10000));
        above.push(b'1');
        assert_eq!(
            bits(&above, SourceNumericRoundingMode::NearestEven),
            0x3ff0_0000_0000_0001
        );
        assert_eq!(
            bits(&above, SourceNumericRoundingMode::TowardZero),
            0x3ff0_0000_0000_0000
        );
        let mut exact = midpoint.to_vec();
        exact.extend(std::iter::repeat_n(b'0', 10000));
        assert_eq!(
            bits(&exact, SourceNumericRoundingMode::NearestEven),
            0x3ff0_0000_0000_0000
        );
        assert_eq!(
            bits(&exact, SourceNumericRoundingMode::Upward),
            0x3ff0_0000_0000_0001
        );
    }
}

#[doc(hidden)]
pub fn property_value_to_double(value: &PropertyValue) -> Result<f64, PropertyDoubleReadError> {
    // RDKit❗✔️: template <class T>
    // RDKit❗✔️: typename boost::enable_if<boost::is_arithmetic<T>, T>::type from_rdvalue(
    // RDKit❗✔️:     RDValue_cast_t arg) {
    // RDKit❗✔️:   T res;
    // RDKit❗✔️:   if (arg.getTag() == RDTypeTag::StringTag) {
    // RDKit❗✔️:     Utils::LocaleSwitcher ls;
    // RDKit❗✔️:     try {
    // RDKit❗✔️:       res = rdvalue_cast<T>(arg);
    // RDKit❗✔️:     } catch (const std::bad_any_cast &exc) {
    // RDKit❗✔️:       try {
    // RDKit❗✔️: 	std::string val = rdvalue_cast<std::string>(arg);
    // RDKit❗✔️: 	// trim only the right characters, this mimics how SD values
    // RDKit❗✔️: 	//  work on read, they will be trimmed by the MolFile parser
    // RDKit❗✔️: 	boost::trim_right(val);
    // RDKit❗✔️:         res = boost::lexical_cast<T>(val);
    // RDKit❗✔️:       } catch (...) {
    // RDKit❗✔️:         throw exc;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     res = rdvalue_cast<T>(arg);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // RDKit❗✔️: template <>
    // RDKit❗✔️: inline double rdvalue_cast<double>(RDValue_cast_t v) {
    // RDKit❗✔️:   if (rdvalue_is<double>(v)) {
    // RDKit❗✔️:     return v.value.d;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (rdvalue_is<float>(v)) {
    // RDKit❗✔️:     return v.value.f;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   throw std::bad_any_cast();
    // RDKit❗✔️: }

    match value {
        PropertyValue::Double(value) => Ok(*value),
        PropertyValue::String(value) => {
            let end = value
                .as_bytes()
                .iter()
                .rposition(|byte| !matches!(byte, b' ' | b'\t' | b'\r' | b'\n' | 0x0b | 0x0c))
                .map_or(0, |index| index + 1);
            source_lexical_double(&value.as_bytes()[..end]).map_err(|source| {
                PropertyDoubleReadError::Lexical {
                    value: value.clone(),
                    source,
                }
            })
        }
        value => Err(PropertyDoubleReadError::InvalidKind { kind: value.kind() }),
    }
}

#[doc(hidden)]
pub fn source_field_double(input: &[u8]) -> Result<f64, DoubleLexicalReadError> {
    // RDKit❗✔️: RDKIT_FILEPARSERS_EXPORT inline std::string_view strip(
    // RDKit❗✔️:     std::string_view orig, std::string stripChars = " \t\r\n") {
    // RDKit❗✔️:   std::string_view res = orig;
    // RDKit❗✔️:   auto start = res.find_first_not_of(stripChars);
    // RDKit❗✔️:   if (start != std::string_view::npos) {
    // RDKit❗✔️:     auto end = res.find_last_not_of(stripChars) + 1;
    // RDKit❗✔️:     res = res.substr(start, end - start);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     res = "";
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // RDKit❗✔️:
    // RDKit❗✔️: template <typename T>
    // RDKit❗✔️: T stripSpacesAndCast(std::string_view input, bool acceptSpaces = false) {
    // RDKit❗✔️:   auto trimmed = strip(input, " ");
    // RDKit❗✔️:   if (acceptSpaces && trimmed.empty()) {
    // RDKit❗✔️:     return 0;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     return boost::lexical_cast<T>(trimmed);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️: template <typename T>
    // RDKit❗✔️: T stripSpacesAndCast(const std::string &input, bool acceptSpaces = false) {
    // RDKit❗✔️:   return stripSpacesAndCast<T>(std::string_view(input.c_str()), acceptSpaces);
    // RDKit❗✔️: }

    // std::string overload creates string_view(input.c_str()), so only this
    // source field boundary ends at NUL. Then strip only ordinary spaces.
    let end = input
        .iter()
        .position(|byte| *byte == 0)
        .unwrap_or(input.len());
    let input = &input[..end];
    let start = input
        .iter()
        .position(|byte| *byte != b' ')
        .unwrap_or(input.len());
    let end = input
        .iter()
        .rposition(|byte| *byte != b' ')
        .map_or(start, |index| index + 1);
    source_lexical_double(&input[start..end])
}

#[doc(hidden)]
pub fn source_unsigned_stream_read(
    bytes: &[u8],
    position: &mut usize,
    failed: &mut bool,
) -> Option<u32> {
    // libstdc++❗✔️:       num_get<_CharT, _InIter>::
    // libstdc++❗✔️:       _M_extract_int(_InIter __beg, _InIter __end, ios_base& __io,
    // libstdc++❗✔️: 		     ios_base::iostate& __err, _ValueT& __v) const
    // libstdc++❗✔️:       {
    // libstdc++❗✔️:         typedef char_traits<_CharT>			    __traits_type;
    // libstdc++❗✔️: 	using __gnu_cxx::__add_unsigned;
    // libstdc++❗✔️: 	typedef typename __add_unsigned<_ValueT>::__type    __unsigned_type;
    // libstdc++❗✔️: 	typedef __numpunct_cache<_CharT>                    __cache_type;
    // libstdc++❗✔️: 	__use_cache<__cache_type> __uc;
    // libstdc++❗✔️: 	const locale& __loc = __io._M_getloc();
    // libstdc++❗✔️: 	const __cache_type* __lc = __uc(__loc);
    // libstdc++❗✔️: 	const _CharT* __lit = __lc->_M_atoms_in;
    // libstdc++❗✔️: 	char_type __c = char_type();
    // libstdc++❗✔️:
    // libstdc++❗✔️: 	// NB: Iff __basefield == 0, __base can change based on contents.
    // libstdc++❗✔️: 	const ios_base::fmtflags __basefield = __io.flags()
    // libstdc++❗✔️: 	                                       & ios_base::basefield;
    // libstdc++❗✔️: 	const bool __oct = __basefield == ios_base::oct;
    // libstdc++❗✔️: 	int __base = __oct ? 8 : (__basefield == ios_base::hex ? 16 : 10);
    // libstdc++❗✔️:
    // libstdc++❗✔️: 	// True if __beg becomes equal to __end.
    // libstdc++❗✔️: 	bool __testeof = __beg == __end;
    // libstdc++❗✔️:
    // libstdc++❗✔️: 	// First check for sign.
    // libstdc++❗✔️: 	bool __negative = false;
    // libstdc++❗✔️: 	if (!__testeof)
    // libstdc++❗✔️: 	  {
    // libstdc++❗✔️: 	    __c = *__beg;
    // libstdc++❗✔️: 	    __negative = __c == __lit[__num_base::_S_iminus];
    // libstdc++❗✔️: 	    if ((__negative || __c == __lit[__num_base::_S_iplus])
    // libstdc++❗✔️: 		&& !(__lc->_M_use_grouping && __c == __lc->_M_thousands_sep)
    // libstdc++❗✔️: 		&& !(__c == __lc->_M_decimal_point))
    // libstdc++❗✔️: 	      {
    // libstdc++❗✔️: 		if (++__beg != __end)
    // libstdc++❗✔️: 		  __c = *__beg;
    // libstdc++❗✔️: 		else
    // libstdc++❗✔️: 		  __testeof = true;
    // libstdc++❗✔️: 	      }
    // libstdc++❗✔️: 	  }
    // libstdc++❗✔️:
    // libstdc++❗✔️: 	// Next, look for leading zeros and check required digits
    // libstdc++❗✔️: 	// for base formats.
    // libstdc++❗✔️: 	bool __found_zero = false;
    // libstdc++❗✔️: 	int __sep_pos = 0;
    // libstdc++❗✔️: 	while (!__testeof)
    // libstdc++❗✔️: 	  {
    // libstdc++❗✔️: 	    if ((__lc->_M_use_grouping && __c == __lc->_M_thousands_sep)
    // libstdc++❗✔️: 		|| __c == __lc->_M_decimal_point)
    // libstdc++❗✔️: 	      break;
    // libstdc++❗✔️: 	    else if (__c == __lit[__num_base::_S_izero]
    // libstdc++❗✔️: 		     && (!__found_zero || __base == 10))
    // libstdc++❗✔️: 	      {
    // libstdc++❗✔️: 		__found_zero = true;
    // libstdc++❗✔️: 		++__sep_pos;
    // libstdc++❗✔️: 		if (__basefield == 0)
    // libstdc++❗✔️: 		  __base = 8;
    // libstdc++❗✔️: 		if (__base == 8)
    // libstdc++❗✔️: 		  __sep_pos = 0;
    // libstdc++❗✔️: 	      }
    // libstdc++❗✔️: 	    else if (__found_zero
    // libstdc++❗✔️: 		     && (__c == __lit[__num_base::_S_ix]
    // libstdc++❗✔️: 			 || __c == __lit[__num_base::_S_iX]))
    // libstdc++❗✔️: 	      {
    // libstdc++❗✔️: 		if (__basefield == 0)
    // libstdc++❗✔️: 		  __base = 16;
    // libstdc++❗✔️: 		if (__base == 16)
    // libstdc++❗✔️: 		  {
    // libstdc++❗✔️: 		    __found_zero = false;
    // libstdc++❗✔️: 		    __sep_pos = 0;
    // libstdc++❗✔️: 		  }
    // libstdc++❗✔️: 		else
    // libstdc++❗✔️: 		  break;
    // libstdc++❗✔️: 	      }
    // libstdc++❗✔️: 	    else
    // libstdc++❗✔️: 	      break;
    // libstdc++❗✔️:
    // libstdc++❗✔️: 	    if (++__beg != __end)
    // libstdc++❗✔️: 	      {
    // libstdc++❗✔️: 		__c = *__beg;
    // libstdc++❗✔️: 		if (!__found_zero)
    // libstdc++❗✔️: 		  break;
    // libstdc++❗✔️: 	      }
    // libstdc++❗✔️: 	    else
    // libstdc++❗✔️: 	      __testeof = true;
    // libstdc++❗✔️: 	  }
    // libstdc++❗✔️:
    // libstdc++❗✔️: 	// At this point, base is determined. If not hex, only allow
    // libstdc++❗✔️: 	// base digits as valid input.
    // libstdc++❗✔️: 	const size_t __len = (__base == 16 ? __num_base::_S_iend
    // libstdc++❗✔️: 			      - __num_base::_S_izero : __base);
    // libstdc++❗✔️:
    // libstdc++❗✔️: 	// Extract.
    // libstdc++❗✔️: 	typedef __gnu_cxx::__numeric_traits<_ValueT> __num_traits;
    // libstdc++❗✔️: 	string __found_grouping;
    // libstdc++❗✔️: 	if (__lc->_M_use_grouping)
    // libstdc++❗✔️: 	  __found_grouping.reserve(32);
    // libstdc++❗✔️: 	bool __testfail = false;
    // libstdc++❗✔️: 	bool __testoverflow = false;
    // libstdc++❗✔️: 	const __unsigned_type __max =
    // libstdc++❗✔️: 	  (__negative && __num_traits::__is_signed)
    // libstdc++❗✔️: 	  ? -static_cast<__unsigned_type>(__num_traits::__min)
    // libstdc++❗✔️: 	  : __num_traits::__max;
    // libstdc++❗✔️: 	const __unsigned_type __smax = __max / __base;
    // libstdc++❗✔️: 	__unsigned_type __result = 0;
    // libstdc++❗✔️: 	int __digit = 0;
    // libstdc++❗✔️: 	const char_type* __lit_zero = __lit + __num_base::_S_izero;
    // libstdc++❗✔️:
    // libstdc++❗✔️: 	if (!__lc->_M_allocated)
    // libstdc++❗✔️: 	  // "C" locale
    // libstdc++❗✔️: 	  while (!__testeof)
    // libstdc++❗✔️: 	    {
    // libstdc++❗✔️: 	      __digit = _M_find(__lit_zero, __len, __c);
    // libstdc++❗✔️: 	      if (__digit == -1)
    // libstdc++❗✔️: 		break;
    // libstdc++❗✔️:
    // libstdc++❗✔️: 	      if (__result > __smax)
    // libstdc++❗✔️: 		__testoverflow = true;
    // libstdc++❗✔️: 	      else
    // libstdc++❗✔️: 		{
    // libstdc++❗✔️: 		  __result *= __base;
    // libstdc++❗✔️: 		  __testoverflow |= __result > __max - __digit;
    // libstdc++❗✔️: 		  __result += __digit;
    // libstdc++❗✔️: 		  ++__sep_pos;
    // libstdc++❗✔️: 		}
    // libstdc++❗✔️:
    // libstdc++❗✔️: 	      if (++__beg != __end)
    // libstdc++❗✔️: 		__c = *__beg;
    // libstdc++❗✔️: 	      else
    // libstdc++❗✔️: 		__testeof = true;
    // libstdc++❗✔️: 	    }
    // libstdc++❗✔️: 	else
    // libstdc++❗✔️: 	  while (!__testeof)
    // libstdc++❗✔️: 	    {
    // libstdc++❗✔️: 	      // According to 22.2.2.1.2, p8-9, first look for thousands_sep
    // libstdc++❗✔️: 	      // and decimal_point.
    // libstdc++❗✔️: 	      if (__lc->_M_use_grouping && __c == __lc->_M_thousands_sep)
    // libstdc++❗✔️: 		{
    // libstdc++❗✔️: 		  // NB: Thousands separator at the beginning of a string
    // libstdc++❗✔️: 		  // is a no-no, as is two consecutive thousands separators.
    // libstdc++❗✔️: 		  if (__sep_pos)
    // libstdc++❗✔️: 		    {
    // libstdc++❗✔️: 		      __found_grouping += static_cast<char>(__sep_pos);
    // libstdc++❗✔️: 		      __sep_pos = 0;
    // libstdc++❗✔️: 		    }
    // libstdc++❗✔️: 		  else
    // libstdc++❗✔️: 		    {
    // libstdc++❗✔️: 		      __testfail = true;
    // libstdc++❗✔️: 		      break;
    // libstdc++❗✔️: 		    }
    // libstdc++❗✔️: 		}
    // libstdc++❗✔️: 	      else if (__c == __lc->_M_decimal_point)
    // libstdc++❗✔️: 		break;
    // libstdc++❗✔️: 	      else
    // libstdc++❗✔️: 		{
    // libstdc++❗✔️: 		  const char_type* __q =
    // libstdc++❗✔️: 		    __traits_type::find(__lit_zero, __len, __c);
    // libstdc++❗✔️: 		  if (!__q)
    // libstdc++❗✔️: 		    break;
    // libstdc++❗✔️:
    // libstdc++❗✔️: 		  __digit = __q - __lit_zero;
    // libstdc++❗✔️: 		  if (__digit > 15)
    // libstdc++❗✔️: 		    __digit -= 6;
    // libstdc++❗✔️: 		  if (__result > __smax)
    // libstdc++❗✔️: 		    __testoverflow = true;
    // libstdc++❗✔️: 		  else
    // libstdc++❗✔️: 		    {
    // libstdc++❗✔️: 		      __result *= __base;
    // libstdc++❗✔️: 		      __testoverflow |= __result > __max - __digit;
    // libstdc++❗✔️: 		      __result += __digit;
    // libstdc++❗✔️: 		      ++__sep_pos;
    // libstdc++❗✔️: 		    }
    // libstdc++❗✔️: 		}
    // libstdc++❗✔️:
    // libstdc++❗✔️: 	      if (++__beg != __end)
    // libstdc++❗✔️: 		__c = *__beg;
    // libstdc++❗✔️: 	      else
    // libstdc++❗✔️: 		__testeof = true;
    // libstdc++❗✔️: 	    }
    // libstdc++❗✔️:
    // libstdc++❗✔️: 	// Digit grouping is checked. If grouping and found_grouping don't
    // libstdc++❗✔️: 	// match, then get very very upset, and set failbit.
    // libstdc++❗✔️: 	if (__found_grouping.size())
    // libstdc++❗✔️: 	  {
    // libstdc++❗✔️: 	    // Add the ending grouping.
    // libstdc++❗✔️: 	    __found_grouping += static_cast<char>(__sep_pos);
    // libstdc++❗✔️:
    // libstdc++❗✔️: 	    if (!std::__verify_grouping(__lc->_M_grouping,
    // libstdc++❗✔️: 					__lc->_M_grouping_size,
    // libstdc++❗✔️: 					__found_grouping))
    // libstdc++❗✔️: 	      __err = ios_base::failbit;
    // libstdc++❗✔️: 	  }
    // libstdc++❗✔️:
    // libstdc++❗✔️: 	// _GLIBCXX_RESOLVE_LIB_DEFECTS
    // libstdc++❗✔️: 	// 23. Num_get overflow result.
    // libstdc++❗✔️: 	if ((!__sep_pos && !__found_zero && !__found_grouping.size())
    // libstdc++❗✔️: 	    || __testfail)
    // libstdc++❗✔️: 	  {
    // libstdc++❗✔️: 	    __v = 0;
    // libstdc++❗✔️: 	    __err = ios_base::failbit;
    // libstdc++❗✔️: 	  }
    // libstdc++❗✔️: 	else if (__testoverflow)
    // libstdc++❗✔️: 	  {
    // libstdc++❗✔️: 	    if (__negative && __num_traits::__is_signed)
    // libstdc++❗✔️: 	      __v = __num_traits::__min;
    // libstdc++❗✔️: 	    else
    // libstdc++❗✔️: 	      __v = __num_traits::__max;
    // libstdc++❗✔️: 	    __err = ios_base::failbit;
    // libstdc++❗✔️: 	  }
    // libstdc++❗✔️: 	else
    // libstdc++❗✔️: 	  __v = __negative ? -__result : __result;
    // libstdc++❗✔️:
    // libstdc++❗✔️: 	if (__testeof)
    // libstdc++❗✔️: 	  __err |= ios_base::eofbit;
    // libstdc++❗✔️: 	return __beg;
    // libstdc++❗✔️:       }
    // libstdc++❗✔️:

    // EOF/sentry failure retains the caller's previous value. A malformed
    // numeric token assigns zero; magnitude overflow assigns UINT_MAX.
    // Signed input wraps only after successful unsigned magnitude extraction.
    // O(consumed bytes), no text decode, temporary string or allocation.
    if *failed {
        return None;
    }
    while bytes
        .get(*position)
        .is_some_and(|byte| matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c))
    {
        *position += 1;
    }
    if *position == bytes.len() {
        *failed = true;
        return None;
    }
    let negative = bytes.get(*position) == Some(&b'-');
    if matches!(bytes.get(*position), Some(b'+' | b'-')) {
        *position += 1;
    }
    let start = *position;
    let mut value = 0u32;
    let mut overflow = false;
    while let Some(byte) = bytes.get(*position).filter(|byte| byte.is_ascii_digit()) {
        match value
            .checked_mul(10)
            .and_then(|value| value.checked_add(u32::from(*byte - b'0')))
        {
            Some(next) if !overflow => value = next,
            _ => overflow = true,
        }
        *position += 1;
    }
    if *position == start {
        *failed = true;
        return Some(0);
    }
    if overflow {
        *failed = true;
        return Some(u32::MAX);
    }
    Some(if negative {
        value.wrapping_neg()
    } else {
        value
    })
}

#[doc(hidden)]
pub fn source_unsigned_stream_array(
    bytes: &[u8],
) -> Result<(Vec<u32>, bool, bool), UnsignedStreamArrayError> {
    // RDKit❗✔️: template <class T>
    // RDKit❗✔️: std::vector<T> ParseV3000Array(std::stringstream &stream, int maxV,
    // RDKit❗✔️:                                bool strictParsing) {
    // RDKit❗✔️:   auto paren = stream.get();  // discard parentheses
    // RDKit❗✔️:   if (paren != '(') {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:         << "WARNING: first character of V3000 array is not '('" << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int count = 0;
    // RDKit❗✔️:   stream >> count;
    // RDKit❗✔️:   std::vector<T> values;
    // RDKit❗✔️:   if (maxV >= 0 && count > static_cast<unsigned int>(maxV)) {
    // RDKit❗✔️:     SGroupWarnOrThrow(strictParsing, "invalid count value");
    // RDKit❗✔️:     return values;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   values.reserve(count);
    // RDKit❗✔️:   T value;
    // RDKit❗✔️:   for (unsigned i = 0; i < count; ++i) {
    // RDKit❗✔️:     stream >> value;
    // RDKit❗✔️:     values.push_back(value);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   paren = stream.get();  // discard parentheses
    // RDKit❗✔️:   if (paren != ')') {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:         << "WARNING: final character of V3000 array is not ')'" << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return values;
    // RDKit❗✔️: }

    // Reuse the sole numeric stream primitive. Count initializes to zero,
    // conversion failure assigns 0/MAX, subsequent sentry failures retain the
    // last assigned element. No token filtering or invented length limit.
    // Source's uninitialized first element after sentry failure has no defined
    // value; report that structural failure instead of fabricating one.
    let first = bytes.first() == Some(&b'(');
    let mut position = usize::from(!bytes.is_empty());
    let mut failed = bytes.is_empty();
    let count = source_unsigned_stream_read(bytes, &mut position, &mut failed).unwrap_or(0);
    let mut values = Vec::new();
    values.try_reserve_exact(count as usize)?;
    let mut previous = None;
    for _ in 0..count {
        if let Some(value) = source_unsigned_stream_read(bytes, &mut position, &mut failed) {
            previous = Some(value);
        }
        values.push(previous.ok_or(UnsignedStreamArrayError::MissingFirstValue)?);
    }
    let last = !failed && bytes.get(position) == Some(&b')');
    Ok((values, first, last))
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PropertyULongReadError {
    #[error("signed value {value} causes negative_overflow converting to uint64_t")]
    Negative { value: i32 },
    #[error("bad_any_cast reading {kind:?} as uint64_t")]
    InvalidKind { kind: PropertyValueKind },
    #[error(
        "bad_any_cast reading string {value:?} as uint64_t: unsigned decimal magnitude exceeds u64"
    )]
    LexicalOverflow { value: PropertyText },
    #[error("bad_any_cast reading string {value:?} as uint64_t: {source}")]
    Lexical {
        value: PropertyText,
        #[source]
        source: UIntLexicalReadError,
    },
}

/// Reached source size_t getter on the pinned unsigned-long64 ABI.
/// Generic Any-backed uint64_t values are supplied by explicit source metadata;
/// every existing represented scalar tag follows the actual source specialization.
#[doc(hidden)]
pub fn property_value_to_ulong(value: &PropertyValue) -> Result<u64, PropertyULongReadError> {
    // RDKit❗✔️: template <class T>
    // RDKit❗✔️: typename boost::enable_if<boost::is_arithmetic<T>, T>::type from_rdvalue(
    // RDKit❗✔️:     RDValue_cast_t arg) {
    // RDKit❗✔️:   T res;
    // RDKit❗✔️:   if (arg.getTag() == RDTypeTag::StringTag) {
    // RDKit❗✔️:     Utils::LocaleSwitcher ls;
    // RDKit❗✔️:     try {
    // RDKit❗✔️:       res = rdvalue_cast<T>(arg);
    // RDKit❗✔️:     } catch (const std::bad_any_cast &exc) {
    // RDKit❗✔️:       try {
    // RDKit❗✔️: 	std::string val = rdvalue_cast<std::string>(arg);
    // RDKit❗✔️: 	// trim only the right characters, this mimics how SD values
    // RDKit❗✔️: 	//  work on read, they will be trimmed by the MolFile parser
    // RDKit❗✔️: 	boost::trim_right(val);
    // RDKit❗✔️:         res = boost::lexical_cast<T>(val);
    // RDKit❗✔️:       } catch (...) {
    // RDKit❗✔️:         throw exc;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     res = rdvalue_cast<T>(arg);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // RDKit❗✔️:
    // RDKit❗✔️: template <>
    // RDKit❗✔️: inline std::uint64_t rdvalue_cast<std::uint64_t>(RDValue_cast_t v) {
    // RDKit❗✔️:   if (rdvalue_is<unsigned int>(v)) {
    // RDKit❗✔️:     return static_cast<std::uint64_t>(v.value.u);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (rdvalue_is<int>(v)) {
    // RDKit❗✔️:     return boost::numeric_cast<std::uint64_t>(v.value.i);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (rdvalue_is<std::any>(v)) {
    // RDKit❗✔️:     return std::any_cast<std::uint64_t>(*v.ptrCast<std::any>());
    // RDKit❗✔️:   }
    // RDKit❗✔️:   throw std::bad_any_cast();
    // RDKit❗✔️: }
    // Boost✔️✔️:             template <typename Type>
    // Boost✔️✔️:             bool shr_unsigned(Type& output) {
    // Boost✔️✔️:                 if (start == finish) return false;
    // Boost✔️✔️:                 CharT const minus = lcast_char_constants<CharT>::minus;
    // Boost✔️✔️:                 CharT const plus = lcast_char_constants<CharT>::plus;
    // Boost✔️✔️:                 bool const has_minus = Traits::eq(minus, *start);
    // Boost✔️✔️:
    // Boost✔️✔️:                 /* We won`t use `start' any more, so no need in decrementing it after */
    // Boost✔️✔️:                 if (has_minus || Traits::eq(plus, *start)) {
    // Boost✔️✔️:                     ++start;
    // Boost✔️✔️:                 }
    // Boost✔️✔️:
    // Boost✔️✔️:                 bool const succeed = lcast_ret_unsigned<Traits, Type, CharT>(output, start, finish).convert();
    // Boost✔️✔️:
    // Boost✔️✔️:                 if (has_minus) {
    // Boost✔️✔️:                     output = static_cast<Type>(0u - output);
    // Boost✔️✔️:                 }
    // Boost✔️✔️:
    // Boost✔️✔️:                 return succeed;
    // Boost✔️✔️:             }
    // Native UInt widens; signed negative typed values throw numeric-cast
    // errors. Counted strings use right-only C-locale trimming and optional
    // sign; textual minus is unsigned wrapping after full-width magnitude.
    // The single bounded unsigned parser below serves u32 and u64, preserving
    // all prior u32 limits/errors; no alternate chemistry parser or heuristic.
    // O(bytes), constant scalar state, no successful-string allocation/UTF8
    // decoding, versus native string copy and reverse weighted conversion.
    match value {
        PropertyValue::UInt(value) => Ok(u64::from(*value)),
        PropertyValue::Int(value) => {
            u64::try_from(*value).map_err(|_| PropertyULongReadError::Negative { value: *value })
        }
        PropertyValue::String(value) => {
            let bytes = trim_source_c_locale_right(value.as_bytes());
            let negative = bytes.first() == Some(&b'-');
            let start = usize::from(negative || bytes.first() == Some(&b'+'));
            let magnitude = parse_decimal_magnitude_with_limit(&bytes[start..], start, u64::MAX)
                .map_err(|source| {
                    if source == UIntLexicalReadError::Overflow {
                        PropertyULongReadError::LexicalOverflow {
                            value: value.clone(),
                        }
                    } else {
                        PropertyULongReadError::Lexical {
                            value: value.clone(),
                            source,
                        }
                    }
                })?;
            Ok(if negative {
                magnitude.wrapping_neg()
            } else {
                magnitude
            })
        }
        other => Err(PropertyULongReadError::InvalidKind { kind: other.kind() }),
    }
}

#[cfg(test)]
mod source590_ulong_tests {
    use super::*;
    #[test]
    fn source590_ulong_preserves_full_width_and_source_tag_asymmetry() {
        assert_eq!(
            property_value_to_ulong(&PropertyValue::UInt(u32::MAX)),
            Ok(u64::from(u32::MAX))
        );
        assert!(matches!(
            property_value_to_ulong(&PropertyValue::Int(-1)),
            Err(PropertyULongReadError::Negative { value: -1 })
        ));
        for (text, expected) in [
            ("4294967296", 4294967296),
            ("18446744073709551615", u64::MAX),
            ("-1", u64::MAX),
            ("-18446744073709551615", 1),
            ("+0001\x0b ", 1),
        ] {
            assert_eq!(
                property_value_to_ulong(&PropertyValue::String(text.into())),
                Ok(expected),
                "{text:?}"
            );
        }
        assert!(matches!(
            property_value_to_ulong(&PropertyValue::String("18446744073709551616".into())),
            Err(PropertyULongReadError::LexicalOverflow { .. })
        ));
        assert!(property_value_to_ulong(&PropertyValue::String(" 1".into())).is_err());
        assert!(property_value_to_ulong(&PropertyValue::Bool(true)).is_err());
        assert!(
            property_value_to_ulong(&PropertyValue::String(PropertyText::from_bytes(&[
                b'1', 0, b'2'
            ])))
            .is_err()
        );
    }
    #[test]
    fn source590_shared_parser_keeps_existing_unsigned32_acceptance_boundary() {
        assert_eq!(
            property_value_to_uint(&PropertyValue::String("4294967295".into())),
            Ok(u32::MAX)
        );
        assert!(property_value_to_uint(&PropertyValue::String("4294967296".into())).is_err());
        assert_eq!(
            property_value_to_int(&PropertyValue::String("-2147483648".into())),
            Ok(i32::MIN)
        );
        assert!(property_value_to_int(&PropertyValue::String("2147483648".into())).is_err());
        assert_eq!(
            property_value_to_uint(&PropertyValue::String("-4294967295".into())),
            Ok(1)
        );
        assert!(property_value_to_uint(&PropertyValue::String("-4294967296".into())).is_err());
    }
}
