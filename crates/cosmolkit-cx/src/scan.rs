use crate::CxParseError;
use cosmolkit_model::PropertyText;

pub(crate) fn read_number(text: &[u8], cursor: &mut usize) -> Result<usize, CxParseError> {
    // RDKit✔️🔝: template <typename Iterator>
    // RDKit✔️🔝: bool read_int(Iterator &first, Iterator last, unsigned int &res) {
    // RDKit✔️🔝:   std::string num = "";
    // RDKit✔️🔝:   while (first <= last && *first >= '0' && *first <= '9') {
    // RDKit✔️🔝:     num += *first;
    // RDKit✔️🔝:     ++first;
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   if (num.empty()) {
    // RDKit✔️🔝:     return false;
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   res = boost::lexical_cast<unsigned int>(num);
    // RDKit✔️🔝:   return true;
    // RDKit✔️🔝: }
    // Actual reached Boost 1.81 unsigned conversion helper:
    // Boost✔️✔️:     namespace detail // lcast_ret_unsigned
    // Boost✔️✔️:     {
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
    // Boost❌❌: #else
    // Boost❌❌:                 std::locale loc;
    // Boost❌❌:                 if (loc == std::locale::classic()) {
    // Boost❌❌:                     return main_convert_loop();
    // Boost❌❌:                 }
    // Boost❌❌:
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
    // Boost❌❌: #endif
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
    // Boost✔️✔️:     }
    // Boost✔️✔️: } // namespace boost
    // Boost✔️✔️:
    // Behavior: the source helper checks each decimal digit and rejects a
    // value above unsigned-int max; zero leading digits remain valid even when
    // its multiplier overflows. The private Rust byte-digit helper uses checked
    // multiplication/addition and accepts that same 0..=u32::MAX set for digit-only
    // input, including arbitrarily many leading zeros. No sign/separator reaches
    // conversion. Return/project an error without an output value on failure.
    // Complexity: both scan and conversion are linear. Borrowing the digit run
    // avoids source num's O(n) allocation while preserving full cursor advance.
    // The Rust scanner parses the same single digit run directly from the
    // borrowed input slice, avoiding the source's incrementally grown string
    // while preserving one forward pass and unsigned-int range semantics.
    let start = *cursor;
    while text.get(*cursor).is_some_and(u8::is_ascii_digit) {
        *cursor += 1;
    }
    // Match num.empty() before conversion. Source byte iterators can point
    // to any byte, with no text-validity requirement.
    // No digits leave the iterator/output untouched; conversion failures occur
    // only after the complete digit run has advanced the iterator.
    if *cursor == start {
        return Err(CxParseError::new(start, "invalid CX integer"));
    }
    parse_decimal_u32(&text[start..*cursor])
        .map(|value| value as usize)
        .map_err(|_| CxParseError::new(start, "invalid CX integer"))
}

pub(crate) fn read_pair(
    text: &[u8],
    cursor: &mut usize,
    separator: u8,
) -> Result<(usize, usize), CxParseError> {
    // RDKit✔️✔️: template <typename Iterator>
    // RDKit✔️✔️: bool read_int_pair(Iterator &first, Iterator last, unsigned int &n1,
    // RDKit✔️✔️:                    unsigned int &n2, char sep = '.') {
    // RDKit✔️✔️:   if (!read_int(first, last, n1)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (first >= last || *first != sep) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   return read_int(first, last, n2);
    // RDKit✔️✔️: }
    let first = read_number(text, cursor)?;
    expect_byte(text, cursor, separator)?;
    let second = read_number(text, cursor)?;
    Ok((first, second))
}

pub(crate) fn parse_delimited_number_list(
    text: &[u8],
    cursor: &mut usize,
    separator: u8,
) -> Result<Vec<usize>, CxParseError> {
    // RDKit✔️🔝: template <typename Iterator>
    // RDKit✔️🔝: bool read_int_list(Iterator &first, Iterator last,
    // RDKit✔️🔝:                    std::vector<unsigned int> &res, char sep = ',') {
    // RDKit✔️🔝:   while (1) {
    // RDKit✔️🔝:     std::string num = "";
    // RDKit✔️🔝:     while (first <= last && *first >= '0' && *first <= '9') {
    // RDKit✔️🔝:       num += *first;
    // RDKit✔️🔝:       ++first;
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:     if (!num.empty()) {
    // RDKit✔️🔝:       res.push_back(boost::lexical_cast<unsigned int>(num));
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:     if (first >= last || *first != sep) {
    // RDKit✔️🔝:       break;
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:     ++first;
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   return true;
    // RDKit✔️🔝: }
    // Each number reuses `read_number`, so the Rust path avoids one temporary
    // owned string per list element while retaining linear scan and order.
    let mut values = Vec::new();
    loop {
        if text.get(*cursor).is_some_and(u8::is_ascii_digit) {
            values.push(read_number(text, cursor)?);
        }
        if text.get(*cursor) != Some(&separator) {
            break;
        }
        *cursor += 1;
    }
    Ok(values)
}

pub(crate) fn read_text_to(
    text: &[u8],
    cursor: &mut usize,
    delimiters: &[u8],
) -> Result<PropertyText, CxParseError> {
    // RDKit✔️✔️: template <typename Iterator>
    // RDKit✔️✔️: std::string read_text_to(Iterator &first, Iterator last, std::string delims) {
    // RDKit✔️✔️:   std::string res = "";
    // RDKit✔️✔️:   Iterator start = first;
    // RDKit✔️✔️:   // EFF: there are certainly faster ways to do this
    // RDKit✔️✔️:   while (first <= last && delims.find_first_of(*first) == std::string::npos) {
    // RDKit✔️✔️:     if (*first == '&' && std::distance(first, last) > 2 &&
    // RDKit✔️✔️:         *(first + 1) == '#') {
    // RDKit✔️✔️:       // escaped char
    // RDKit✔️✔️:       if (start != first) {
    // RDKit✔️✔️:         res += std::string(start, first);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       Iterator next = first + 2;
    // RDKit✔️✔️:       while (next != last && *next >= '0' && *next <= '9') {
    // RDKit✔️✔️:         ++next;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (next == last || *next != ';') {
    // RDKit✔️✔️:         throw RDKit::SmilesParseException(
    // RDKit✔️✔️:             "failure parsing CXSMILES extensions: quoted block not terminated "
    // RDKit✔️✔️:             "with ';'");
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (next > first + 2) {
    // RDKit✔️✔️:         std::string blk = std::string(first + 2, next);
    // RDKit✔️✔️:         res += (char)(boost::lexical_cast<int>(blk));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       first = next + 1;
    // RDKit✔️✔️:       start = first;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (start != first) {
    // RDKit✔️✔️:     res += std::string(start, first);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Reached Boost 1.81 signed-int conversion, digit-only CX input:
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
    // Reached unsigned-magnitude conversion helper; C locale branch only:
    // Boost✔️✔️:     namespace detail // lcast_ret_unsigned
    // Boost✔️✔️:     {
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
    // Boost❌❌: #else
    // Boost❌❌:                 std::locale loc;
    // Boost❌❌:                 if (loc == std::locale::classic()) {
    // Boost❌❌:                     return main_convert_loop();
    // Boost❌❌:                 }
    // Boost❌❌:
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
    // Boost❌❌: #endif
    // Boost✔️✔️:             }
    // Boost✔️✔️:
    // Boost✔️✔️:         private:
    // Boost❌❌:             // Iteration that does not care about grouping/separators and assumes that all
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
    // Boost✔️✔️:     }
    // Boost✔️✔️: } // namespace boost
    // Boost✔️✔️:
    // Boost✔️✔️: #endif // BOOST_LEXICAL_CAST_DETAIL_LCAST_UNSIGNED_CONVERTERS_HPP
    // Boost✔️✔️:
    // Behavior: counted source bytes remain counted bytes; escaped digit-only
    // signed-int conversion retains its low eight bits, including NUL and high
    // bytes. The distance > 2 guard keeps a final literal "&#" unchanged.
    // Complexity: one forward scan and amortized buffer growth, matching the
    // source; borrowed digit conversion removes its temporary blk allocation.
    let mut result = PropertyText::new();
    let mut segment_start = *cursor;
    while let Some(&byte) = text.get(*cursor) {
        if delimiters.contains(&byte) {
            break;
        }
        if byte == b'&' && text.len() - *cursor > 2 && text.get(*cursor + 1) == Some(&b'#') {
            result.extend_bytes(&text[segment_start..*cursor]);
            let entity_start = *cursor;
            let mut next = *cursor + 2;
            while text.get(next).is_some_and(u8::is_ascii_digit) {
                next += 1;
            }
            if text.get(next) != Some(&b';') {
                return Err(CxParseError::new(
                    entity_start,
                    "failure parsing CXSMILES extensions: quoted block not terminated with ';'",
                ));
            }
            if next > entity_start + 2 {
                let value = parse_decimal_u32(&text[entity_start + 2..next])
                    .and_then(|value| i32::try_from(value).map_err(|_| ()))
                    .map_err(|_| CxParseError::new(entity_start, "invalid CX character code"))?;
                result.push_byte(value as u8);
            }
            *cursor = next + 1;
            segment_start = *cursor;
        } else {
            *cursor += 1;
        }
    }
    result.extend_bytes(&text[segment_start..*cursor]);
    Ok(result)
}

pub(crate) fn expect_byte(
    text: &[u8],
    cursor: &mut usize,
    expected: u8,
) -> Result<(), CxParseError> {
    if text.get(*cursor) == Some(&expected) {
        *cursor += 1;
        Ok(())
    } else {
        Err(CxParseError::new(
            *cursor,
            format!("expected '{}', found CX syntax mismatch", expected as char),
        ))
    }
}

fn parse_decimal_u32(digits: &[u8]) -> Result<u32, ()> {
    // Boost✔️✔️:     namespace detail // lcast_ret_unsigned
    // Boost✔️✔️:     {
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
    // Boost❌❌: #else
    // Boost❌❌:                 std::locale loc;
    // Boost❌❌:                 if (loc == std::locale::classic()) {
    // Boost❌❌:                     return main_convert_loop();
    // Boost❌❌:                 }
    // Boost❌❌:
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
    // Boost❌❌: #endif
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
    // Boost✔️✔️:     }
    // Boost✔️✔️: } // namespace boost
    // Boost✔️✔️:
    // Boost✔️✔️: #endif // BOOST_LEXICAL_CAST_DETAIL_LCAST_UNSIGNED_CONVERTERS_HPP
    // Boost✔️✔️:
    // Behavior: this caller supplies only a nonempty decimal digit span in C
    // locale. Checked forward accumulation accepts exactly the same u32 values
    // as the source reverse multiplier algorithm, including any leading zeros.
    // Other locale grouping is not modeled and no locale branch is inferred.
    // Complexity: one pass, constant storage; no string allocation or decoding.
    if digits.is_empty() {
        return Err(());
    }
    digits.iter().try_fold(0_u32, |value, &byte| {
        if !byte.is_ascii_digit() {
            return Err(());
        }
        value
            .checked_mul(10)
            .and_then(|value| value.checked_add(u32::from(byte - b'0')))
            .ok_or(())
    })
}
