// MMFF table lexical owners and structured parse errors.
//
// Source notice for the Boost-derived helpers added in this module:
// Copyright Kevlin Henney, 2000-2005.
// Copyright Alexander Nasonov, 2006-2010.
// Copyright Antony Polukhin, 2011-2024.
// Copyright John R. Bandela 2001.
// Distributed under the Boost Software License, Version 1.0.

use std::fmt;

/// Borrowed records from the pinned RDKit `getLine`/MMFF constructor loop.
pub(super) struct SourceLines<'a> {
    remaining: &'a str,
}

pub(super) fn source_lines(source: &str) -> SourceLines<'_> {
    SourceLines { remaining: source }
}

impl<'a> Iterator for SourceLines<'a> {
    type Item = &'a str;

    fn next(&mut self) -> Option<Self::Item> {
        // RDKit helper source, third_party/rdkit/Code/RDGeneral/StreamOps.h:337-345:
        // RDKit✔️🔝: inline std::string getLine(std::istream *inStream) {
        // RDKit✔️🔝:   std::string res;
        // RDKit✔️🔝:   std::getline(*inStream, res);
        // RDKit✔️🔝:   if (!res.empty() && (res.back() == '\r')) {
        // RDKit✔️🔝:     res.resize(res.length() - 1);
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:   return res;
        // RDKit✔️🔝: }
        // RDKit caller overload, third_party/rdkit/Code/RDGeneral/StreamOps.h:347-349:
        // RDKit✔️🔝: inline std::string getLine(std::istream &inStream) {
        // RDKit✔️🔝:   return getLine(&inStream);
        // RDKit✔️🔝: }
        // Both MMFF constructors (Params.cpp:45-94, 427-478) first read an
        // inLine, then process it only while `!(inStream.eof())`, and read
        // another line at the loop end. A non-LF-terminated final line sets
        // EOF during that read and is therefore never processed.
        // RDKit✔️🔝:   std::string inLine = RDKit::getLine(inStream);
        // RDKit✔️🔝:   while (!(inStream.eof())) {
        // RDKit✔️🔝:     inLine = RDKit::getLine(inStream);
        // RDKit✔️🔝:   }
        // Behavior review: find LF-delimited records, remove exactly one
        // trailing CR, and discard the unterminated suffix as the constructor
        // loop does. Complexity review: scan forward once and borrow slices;
        // unlike the source helper, do not construct/copy an owned string per
        // record, which avoids that per-line copying cost without changing bytes.
        let Some(line_end) = self.remaining.find('\n') else {
            self.remaining = "";
            return None;
        };

        let (record, after_lf) = self.remaining.split_at(line_end);
        self.remaining = &after_lf[1..];
        Some(record.strip_suffix('\r').unwrap_or(record))
    }
}

impl std::iter::FusedIterator for SourceLines<'_> {}

/// Tokenize the frozen MMFF table format: dropped tab delimiters, no kept
/// delimiters, and Boost's default empty-token dropping policy.
pub(super) fn tokenize_mmff_line<'a>(line: &'a str) -> impl Iterator<Item = &'a str> + 'a {
    // Boost.Tokenizer 1.85.0 source, include/boost/token_functions.hpp:424-557.
    // This ports the MMFF call-site specialization char_separator<char>("\t"):
    // the general keep-empty, kept-delimiter, and locale modes are not APIs here.
    // Boost✔️🔝:   enum empty_token_policy { drop_empty_tokens, keep_empty_tokens };
    // Boost✔️🔝:
    // Boost✔️🔝:   // The out of the box GCC 2.95 on cygwin does not have a char_traits class.
    // Boost✔️🔝:   template <typename Char,
    // Boost✔️🔝:     typename Tr = BOOST_DEDUCED_TYPENAME std::basic_string<Char>::traits_type >
    // Boost✔️🔝:   class char_separator
    // Boost✔️🔝:   {
    // Boost✔️🔝:     typedef tokenizer_detail::traits_extension<Tr> Traits;
    // Boost✔️🔝:     typedef std::basic_string<Char,Tr> string_type;
    // Boost✔️🔝:   public:
    // Boost✔️🔝:     explicit
    // Boost✔️🔝:     char_separator(const Char* dropped_delims,
    // Boost✔️🔝:                    const Char* kept_delims = 0,
    // Boost✔️🔝:                    empty_token_policy empty_tokens = drop_empty_tokens)
    // Boost✔️🔝:       : m_dropped_delims(dropped_delims),
    // Boost✔️🔝:         m_use_ispunct(false),
    // Boost✔️🔝:         m_use_isspace(false),
    // Boost✔️🔝:         m_empty_tokens(empty_tokens),
    // Boost✔️🔝:         m_output_done(false)
    // Boost✔️🔝:     {
    // Boost✔️🔝:       // Borland workaround
    // Boost✔️🔝:       if (kept_delims)
    // Boost✔️🔝:         m_kept_delims = kept_delims;
    // Boost✔️🔝:     }
    // Boost✔️🔝:
    // Boost✔️🔝:                 // use ispunct() for kept delimiters and isspace for dropped.
    // Boost✔️🔝:     explicit
    // Boost✔️🔝:     char_separator()
    // Boost✔️🔝:       : m_use_ispunct(true),
    // Boost✔️🔝:         m_use_isspace(true),
    // Boost✔️🔝:         m_empty_tokens(drop_empty_tokens),
    // Boost✔️🔝:         m_output_done(false) { }
    // Boost✔️🔝:
    // Boost✔️🔝:     void reset() { }
    // Boost✔️🔝:
    // Boost✔️🔝:     template <typename InputIterator, typename Token>
    // Boost✔️🔝:     bool operator()(InputIterator& next, InputIterator end, Token& tok)
    // Boost✔️🔝:     {
    // Boost✔️🔝:       typedef tokenizer_detail::assign_or_plus_equal<
    // Boost✔️🔝:         BOOST_DEDUCED_TYPENAME tokenizer_detail::get_iterator_category<
    // Boost✔️🔝:           InputIterator
    // Boost✔️🔝:         >::iterator_category
    // Boost✔️🔝:       > assigner;
    // Boost✔️🔝:
    // Boost✔️🔝:       assigner::clear(tok);
    // Boost✔️🔝:
    // Boost✔️🔝:       // skip past all dropped_delims
    // Boost✔️🔝:       if (m_empty_tokens == drop_empty_tokens)
    // Boost✔️🔝:         for (; next != end  && is_dropped(*next); ++next)
    // Boost✔️🔝:           { }
    // Boost✔️🔝:
    // Boost✔️🔝:       InputIterator start(next);
    // Boost✔️🔝:
    // Boost✔️🔝:       if (m_empty_tokens == drop_empty_tokens) {
    // Boost✔️🔝:
    // Boost✔️🔝:         if (next == end)
    // Boost✔️🔝:           return false;
    // Boost✔️🔝:
    // Boost✔️🔝:
    // Boost✔️🔝:         // if we are on a kept_delims move past it and stop
    // Boost✔️🔝:         if (is_kept(*next)) {
    // Boost✔️🔝:           assigner::plus_equal(tok,*next);
    // Boost✔️🔝:           ++next;
    // Boost✔️🔝:         } else
    // Boost✔️🔝:           // append all the non delim characters
    // Boost✔️🔝:           for (; next != end && !is_dropped(*next) && !is_kept(*next); ++next)
    // Boost✔️🔝:             assigner::plus_equal(tok,*next);
    // Boost✔️🔝:       }
    // Boost✔️🔝:       else { // m_empty_tokens == keep_empty_tokens
    // Boost✔️🔝:
    // Boost✔️🔝:         // Handle empty token at the end
    // Boost✔️🔝:         if (next == end)
    // Boost✔️🔝:         {
    // Boost✔️🔝:           if (m_output_done == false)
    // Boost✔️🔝:           {
    // Boost✔️🔝:             m_output_done = true;
    // Boost✔️🔝:             assigner::assign(start,next,tok);
    // Boost✔️🔝:             return true;
    // Boost✔️🔝:           }
    // Boost✔️🔝:           else
    // Boost✔️🔝:             return false;
    // Boost✔️🔝:         }
    // Boost✔️🔝:
    // Boost✔️🔝:         if (is_kept(*next)) {
    // Boost✔️🔝:           if (m_output_done == false)
    // Boost✔️🔝:             m_output_done = true;
    // Boost✔️🔝:           else {
    // Boost✔️🔝:             assigner::plus_equal(tok,*next);
    // Boost✔️🔝:             ++next;
    // Boost✔️🔝:             m_output_done = false;
    // Boost✔️🔝:           }
    // Boost✔️🔝:         }
    // Boost✔️🔝:         else if (m_output_done == false && is_dropped(*next)) {
    // Boost✔️🔝:           m_output_done = true;
    // Boost✔️🔝:         }
    // Boost✔️🔝:         else {
    // Boost✔️🔝:           if (is_dropped(*next))
    // Boost✔️🔝:             start=++next;
    // Boost✔️🔝:           for (; next != end && !is_dropped(*next) && !is_kept(*next); ++next)
    // Boost✔️🔝:             assigner::plus_equal(tok,*next);
    // Boost✔️🔝:           m_output_done = true;
    // Boost✔️🔝:         }
    // Boost✔️🔝:       }
    // Boost✔️🔝:       assigner::assign(start,next,tok);
    // Boost✔️🔝:       return true;
    // Boost✔️🔝:     }
    // Boost✔️🔝:
    // Boost✔️🔝:   private:
    // Boost✔️🔝:     string_type m_kept_delims;
    // Boost✔️🔝:     string_type m_dropped_delims;
    // Boost✔️🔝:     bool m_use_ispunct;
    // Boost✔️🔝:     bool m_use_isspace;
    // Boost✔️🔝:     empty_token_policy m_empty_tokens;
    // Boost✔️🔝:     bool m_output_done;
    // Boost✔️🔝:
    // Boost✔️🔝:     bool is_kept(Char E) const
    // Boost✔️🔝:     {
    // Boost✔️🔝:       if (m_kept_delims.length())
    // Boost✔️🔝:         return m_kept_delims.find(E) != string_type::npos;
    // Boost✔️🔝:       else if (m_use_ispunct) {
    // Boost✔️🔝:         return Traits::ispunct(E) != 0;
    // Boost✔️🔝:       } else
    // Boost✔️🔝:         return false;
    // Boost✔️🔝:     }
    // Boost✔️🔝:     bool is_dropped(Char E) const
    // Boost✔️🔝:     {
    // Boost✔️🔝:       if (m_dropped_delims.length())
    // Boost✔️🔝:         return m_dropped_delims.find(E) != string_type::npos;
    // Boost✔️🔝:       else if (m_use_isspace) {
    // Boost✔️🔝:         return Traits::isspace(E) != 0;
    // Boost✔️🔝:       } else
    // Boost✔️🔝:         return false;
    // Boost✔️🔝:     }
    // Boost✔️🔝:   };
    // Behavior review: the configured source drops only tab and drops empty
    // tokens; split plus this predicate preserves spaces, CR, NUL, and all other
    // cell bytes. Complexity review: one forward pass with borrowed slices and
    // constant state, avoiding source token string copies and any token Vec.
    line.split('\t').filter(|token| !token.is_empty())
}

/// Parse the frozen classic-locale unsigned MMFF table integer spelling.
pub(super) fn parse_mmff_u32(cell: &str) -> Option<u32> {
    // Boost lexical-cast source pin 02e5821ab32c45fad719829e9644e5d681c9ba0b:
    // converter_lexical_streams.hpp:496-514 (`shr_unsigned`):
    // Boost✔️✔️:         bool shr_unsigned(Type& output) {
    // Boost✔️✔️:             if (start == finish) return false;
    // Boost✔️✔️:             CharT const minus = lcast_char_constants<CharT>::minus;
    // Boost✔️✔️:             CharT const plus = lcast_char_constants<CharT>::plus;
    // Boost✔️✔️:             bool const has_minus = Traits::eq(minus, *start);
    // Boost✔️✔️:
    // Boost✔️✔️:             /* We won`t use `start' any more, so no need in decrementing it after */
    // Boost✔️✔️:             if (has_minus || Traits::eq(plus, *start)) {
    // Boost✔️✔️:                 ++start;
    // Boost✔️✔️:             }
    // Boost✔️✔️:
    // Boost✔️✔️:             bool const succeed = lcast_ret_unsigned<Traits, Type, CharT>(output, start, finish).convert();
    // Boost✔️✔️:
    // Boost✔️✔️:             if (has_minus) {
    // Boost✔️✔️:                 output = static_cast<Type>(0u - output);
    // Boost✔️✔️:             }
    // Boost✔️✔️:
    // Boost✔️✔️:             return succeed;
    // Boost✔️✔️:         }
    // lcast_unsigned_converters.hpp:160-294 (`lcast_ret_unsigned`):
    // Boost✔️✔️:         template <class Traits, class T, class CharT>
    // Boost✔️✔️:         class lcast_ret_unsigned: boost::noncopyable {
    // Boost✔️✔️:             bool m_multiplier_overflowed;
    // Boost✔️✔️:             T m_multiplier;
    // Boost✔️✔️:             T& m_value;
    // Boost✔️✔️:             const CharT* const m_begin;
    // Boost✔️✔️:             const CharT* m_end;
    // Boost✔️✔️:
    // Boost✔️✔️:         public:
    // Boost✔️✔️:             lcast_ret_unsigned(T& value, const CharT* const begin, const CharT* end) noexcept
    // Boost✔️✔️:                 : m_multiplier_overflowed(false), m_multiplier(1), m_value(value), m_begin(begin), m_end(end)
    // Boost✔️✔️:             {
    // Boost✔️✔️: #ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
    // Boost✔️✔️:                 static_assert(!std::numeric_limits<T>::is_signed, "");
    // Boost✔️✔️:
    // Boost✔️✔️:                 // GCC when used with flag -std=c++0x may not have std::numeric_limits
    // Boost✔️✔️:                 // specializations for __int128 and unsigned __int128 types.
    // Boost✔️✔️:                 // Try compilation with -std=gnu++0x or -std=gnu++11.
    // Boost✔️✔️:                 //
    // Boost✔️✔️:                 // http://gcc.gnu.org/bugzilla/show_bug.cgi?id=40856
    // Boost✔️✔️:                 static_assert(std::numeric_limits<T>::is_specialized,
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
    // Boost❌❌: #endif
    // Boost✔️✔️:             }
    // Boost✔️✔️:
    // Boost✔️✔️:         private:
    // Boost✔️✔️:             // Iteration that does not care about grouping/separators and assumes that all
    // Boost✔️✔️:             // input characters are digits
    // Boost✔️✔️: #if defined(__clang__) && (__clang_major__ > 3 || __clang_minor__ > 6)
    // Boost✔️✔️:             __attribute__((no_sanitize("unsigned-integer-overflow")))
    // Boost✔️✔️: #endif
    // Boost✔️✔️:             inline bool main_convert_iteration() noexcept {
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
    // Boost✔️✔️:             bool main_convert_loop() noexcept {
    // Boost✔️✔️:                 for ( ; m_end >= m_begin; --m_end) {
    // Boost✔️✔️:                     if (!main_convert_iteration()) {
    // Boost✔️✔️:                         return false;
    // Boost✔️✔️:                     }
    // Boost✔️✔️:                 }
    // Boost✔️✔️:
    // Boost✔️✔️:                 return true;
    // Boost✔️✔️:             }
    // Boost✔️✔️:         };
    // Behavior review: the configured source consumes one optional sign and
    // requires every remaining byte to be an ASCII decimal digit with a
    // checked u32 magnitude; minus wraps only after successful conversion.
    // The source's general non-classic locale grouping branch is anchored above
    // but remains outside this frozen C/classic-locale specialization.
    // Complexity review: accumulate once over the borrowed bytes with scalar
    // checked multiply/add state: O(n), O(1), and no cell clone. Leading zero
    // multiplication remains zero, so arbitrary zero prefixes stay valid.
    let bytes = cell.as_bytes();
    let (negative, digits) = match bytes.first().copied() {
        Some(b'-') => (true, &bytes[1..]),
        Some(b'+') => (false, &bytes[1..]),
        _ => (false, bytes),
    };

    if digits.is_empty() {
        return None;
    }

    let mut magnitude = 0_u32;
    for &byte in digits {
        if !byte.is_ascii_digit() {
            return None;
        }
        magnitude = magnitude
            .checked_mul(10)?
            .checked_add(u32::from(byte - b'0'))?;
    }

    Some(if negative {
        0_u32.wrapping_sub(magnitude)
    } else {
        magnitude
    })
}

/// Parse the pinned classic-locale binary64 MMFF cell spelling.
pub(super) fn parse_mmff_f64(cell: &str) -> Option<f64> {
    // Boost lexical_cast source pin 02e5821ab32c45fad719829e9644e5d681c9ba0b:
    // converter_lexical_streams.hpp::to_target_stream::shr_using_base_class:
    // Boost❗✔️:         template<typename InputStreamable>
    // Boost❗✔️:         bool shr_using_base_class(InputStreamable& output)
    // Boost❗✔️:         {
    // Boost❗✔️:             static_assert(
    // Boost❗✔️:                 !boost::is_pointer<InputStreamable>::value,
    // Boost❗✔️:                 "boost::lexical_cast can not convert to pointers"
    // Boost❗✔️:             );
    // Boost❗✔️:
    // Boost❗✔️: #if defined(BOOST_NO_STRINGSTREAM) || defined(BOOST_NO_STD_LOCALE)
    // Boost❗✔️:             static_assert(boost::is_same<char, CharT>::value,
    // Boost❗✔️:                 "boost::lexical_cast can not convert, because your STL library does not "
    // Boost❗✔️:                 "support such conversions. Try updating it."
    // Boost❗✔️:             );
    // Boost❗✔️: #endif
    // Boost❗✔️:
    // Boost❗✔️: #if defined(BOOST_NO_STRINGSTREAM)
    // Boost❗✔️:             std::istrstream stream(start, static_cast<std::istrstream::streamsize>(finish - start));
    // Boost❗✔️: #else
    // Boost❗✔️:             typedef detail::lcast::buffer_t<CharT, Traits> buffer_t;
    // Boost❗✔️:             buffer_t buf;
    // Boost❗✔️:             // Usually `istream` and `basic_istream` do not modify
    // Boost❗✔️:             // content of buffer; `buffer_t` assures that this is true
    // Boost❗✔️:             buf.setbuf(const_cast<CharT*>(start), static_cast<typename buffer_t::streamsize>(finish - start));
    // Boost❗✔️: #if defined(BOOST_NO_STD_LOCALE)
    // Boost❗✔️:             std::istream stream(&buf);
    // Boost❗✔️: #else
    // Boost❗✔️:             std::basic_istream<CharT, Traits> stream(&buf);
    // Boost❗✔️: #endif // BOOST_NO_STD_LOCALE
    // Boost❗✔️: #endif // BOOST_NO_STRINGSTREAM
    // Boost❗✔️:
    // Boost❗✔️: #ifndef BOOST_NO_EXCEPTIONS
    // Boost❗✔️:             stream.exceptions(std::ios::badbit);
    // Boost❗✔️:             try {
    // Boost❗✔️: #endif
    // Boost❗✔️:             stream.unsetf(std::ios::skipws);
    // Boost❗✔️:             lcast_set_precision(stream, static_cast<InputStreamable*>(0));
    // Boost❗✔️:
    // Boost❗✔️:             return (stream >> output)
    // Boost❗✔️:                 && (stream.get() == Traits::eof());
    // Boost❗✔️:
    // Boost❗✔️: #ifndef BOOST_NO_EXCEPTIONS
    // Boost❗✔️:             } catch (const ::std::ios_base::failure& /*f*/) {
    // Boost❗✔️:                 return false;
    // Boost❗✔️:             }
    // Boost❗✔️: #endif
    // Boost❗✔️:         }
    //
    // converter_lexical_streams.hpp::to_target_stream::float_types_converter_internal:
    // Boost❗✔️:         template <class T>
    // Boost❗✔️:         bool float_types_converter_internal(T& output) {
    // Boost❗✔️:             if (parse_inf_nan(start, finish, output)) return true;
    // Boost❗✔️:             bool const return_value = shr_using_base_class(output);
    // Boost❗✔️:
    // Boost❗✔️:             /* Some compilers and libraries successfully
    // Boost❗✔️:              * parse 'inf', 'INFINITY', '1.0E', '1.0E-'...
    // Boost❗✔️:              * We are trying to provide a unified behaviour,
    // Boost❗✔️:              * so we just forbid such conversions (as some
    // Boost❗✔️:              * of the most popular compilers/libraries do)
    // Boost❗✔️:              * */
    // Boost❗✔️:             CharT const minus = lcast_char_constants<CharT>::minus;
    // Boost❗✔️:             CharT const plus = lcast_char_constants<CharT>::plus;
    // Boost❗✔️:             CharT const capital_e = lcast_char_constants<CharT>::capital_e;
    // Boost❗✔️:             CharT const lowercase_e = lcast_char_constants<CharT>::lowercase_e;
    // Boost❗✔️:             if ( return_value &&
    // Boost❗✔️:                  (
    // Boost❗✔️:                     Traits::eq(*(finish-1), lowercase_e)                   // 1.0e
    // Boost❗✔️:                     || Traits::eq(*(finish-1), capital_e)                  // 1.0E
    // Boost❗✔️:                     || Traits::eq(*(finish-1), minus)                      // 1.0e- or 1.0E-
    // Boost❗✔️:                     || Traits::eq(*(finish-1), plus)                       // 1.0e+ or 1.0E+
    // Boost❗✔️:                  )
    // Boost❗✔️:             ) return false;
    // Boost❗✔️:
    // Boost❗✔️:             return return_value;
    // Boost❗✔️:         }
    //
    // converter_lexical_streams.hpp::to_target_stream::stream_out(double&):
    // Boost❗✔️:         bool stream_out(double& output) { return float_types_converter_internal(output); }
    //
    // inf_nan.hpp::lc_iequal:
    // Boost❗✔️:         template <class CharT>
    // Boost❗✔️:         bool lc_iequal(const CharT* val, const CharT* lcase, const CharT* ucase, unsigned int len) noexcept {
    // Boost❗✔️:             for( unsigned int i=0; i < len; ++i ) {
    // Boost❗✔️:                 if ( val[i] != lcase[i] && val[i] != ucase[i] ) return false;
    // Boost❗✔️:             }
    // Boost❗✔️:
    // Boost❗✔️:             return true;
    // Boost❗✔️:         }
    //
    // inf_nan.hpp::parse_inf_nan_impl:
    // Boost❗✔️:         template <class CharT, class T>
    // Boost❗✔️:         inline bool parse_inf_nan_impl(const CharT* begin, const CharT* end, T& value
    // Boost❗✔️:             , const CharT* lc_NAN, const CharT* lc_nan
    // Boost❗✔️:             , const CharT* lc_INFINITY, const CharT* lc_infinity
    // Boost❗✔️:             , const CharT opening_brace, const CharT closing_brace) noexcept
    // Boost❗✔️:         {
    // Boost❗✔️:             if (begin == end) return false;
    // Boost❗✔️:             const CharT minus = lcast_char_constants<CharT>::minus;
    // Boost❗✔️:             const CharT plus = lcast_char_constants<CharT>::plus;
    // Boost❗✔️:             const int inifinity_size = 8; // == sizeof("infinity") - 1
    // Boost❗✔️:
    // Boost❗✔️:             /* Parsing +/- */
    // Boost❗✔️:             bool const has_minus = (*begin == minus);
    // Boost❗✔️:             if (has_minus || *begin == plus) {
    // Boost❗✔️:                 ++ begin;
    // Boost❗✔️:             }
    // Boost❗✔️:
    // Boost❗✔️:             if (end - begin < 3) return false;
    // Boost❗✔️:             if (lc_iequal(begin, lc_nan, lc_NAN, 3)) {
    // Boost❗✔️:                 begin += 3;
    // Boost❗✔️:                 if (end != begin) {
    // Boost❗✔️:                     /* It is 'nan(...)' or some bad input*/
    // Boost❗✔️:
    // Boost❗✔️:                     if (end - begin < 2) return false; // bad input
    // Boost❗✔️:                     -- end;
    // Boost❗✔️:                     if (*begin != opening_brace || *end != closing_brace) return false; // bad input
    // Boost❗✔️:                 }
    // Boost❗✔️:
    // Boost❗✔️:                 if( !has_minus ) value = std::numeric_limits<T>::quiet_NaN();
    // Boost❗✔️:                 else value = boost::core::copysign(std::numeric_limits<T>::quiet_NaN(), static_cast<T>(-1));
    // Boost❗✔️:                 return true;
    // Boost❗✔️:             } else if (
    // Boost❗✔️:                 ( /* 'INF' or 'inf' */
    // Boost❗✔️:                   end - begin == 3      // 3 == sizeof('inf') - 1
    // Boost❗✔️:                   && lc_iequal(begin, lc_infinity, lc_INFINITY, 3)
    // Boost❗✔️:                 )
    // Boost❗✔️:                 ||
    // Boost❗✔️:                 ( /* 'INFINITY' or 'infinity' */
    // Boost❗✔️:                   end - begin == inifinity_size
    // Boost❗✔️:                   && lc_iequal(begin, lc_infinity, lc_INFINITY, inifinity_size)
    // Boost❗✔️:                 )
    // Boost❗✔️:              )
    // Boost❗✔️:             {
    // Boost❗✔️:                 if( !has_minus ) value = std::numeric_limits<T>::infinity();
    // Boost❗✔️:                 else value = -std::numeric_limits<T>::infinity();
    // Boost❗✔️:                 return true;
    // Boost❗✔️:             }
    // Boost❗✔️:
    // Boost❗✔️:             return false;
    // Boost❗✔️:         }
    //
    // inf_nan.hpp::parse_inf_nan(char):
    // Boost❗✔️:         template <class CharT, class T>
    // Boost❗✔️:         bool parse_inf_nan(const CharT* begin, const CharT* end, T& value) noexcept {
    // Boost❗✔️:             return parse_inf_nan_impl(begin, end, value
    // Boost❗✔️:                                , "NAN", "nan"
    // Boost❗✔️:                                , "INFINITY", "infinity"
    // Boost❗✔️:                                , '(', ')');
    // Boost❗✔️:         }
    //
    // Behavior review: parse source NaN/Inf tokens before decimal conversion.
    // Preserve its optional sign and require only the NaN payload's first and
    // final parentheses; payload bytes are ignored. Decimal parsing consumes
    // the complete borrowed cell, rejects overflow to infinity, and preserves
    // successful underflow as signed zero.
    // Complexity review: source uses one stream conversion; this code scans
    // borrowed bytes without allocation or cloning, then delegates decimal
    // conversion to the Rust binary64 parser in O(n) time.
    let bytes = cell.as_bytes();
    let (negative, unsigned) = match bytes.first().copied() {
        Some(b'-') => (true, &bytes[1..]),
        Some(b'+') => (false, &bytes[1..]),
        _ => (false, bytes),
    };

    let equals_ascii_case = |value: &[u8], expected: &[u8]| {
        value.len() == expected.len()
            && value
                .iter()
                .zip(expected)
                .all(|(actual, expected)| actual.eq_ignore_ascii_case(expected))
    };

    if equals_ascii_case(unsigned, b"nan")
        || (unsigned.len() >= 5
            && equals_ascii_case(&unsigned[..3], b"nan")
            && unsigned[3] == b'('
            && unsigned.last() == Some(&b')'))
    {
        return Some(f64::from_bits(if negative {
            0xfff8_0000_0000_0000
        } else {
            0x7ff8_0000_0000_0000
        }));
    }

    if equals_ascii_case(unsigned, b"inf") || equals_ascii_case(unsigned, b"infinity") {
        return Some(if negative {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        });
    }

    let value = cell.parse::<f64>().ok()?;
    value.is_finite().then_some(value)
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum MmffParamTable {
    Def,
    Prop,
    Pbci,
    Chg,
    Bond,
    Bndk,
    Stbn,
    Dfsb,
    Angle,
    Oop,
    Tor,
    Vdw,
    HerschbachLaurie,
    CovRadPauEle,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum MmffParamParseCause {
    InvalidUnsigned { cell: String },
    InvalidFloat { cell: String },
    MissingToken,
    EmptyProcessedLine,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct MmffParamParseError {
    pub(super) table: MmffParamTable,
    pub(super) line: usize,
    pub(super) column: usize,
    pub(super) cause: MmffParamParseCause,
}

impl fmt::Display for MmffParamParseError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let table = match self.table {
            MmffParamTable::Def => "Def",
            MmffParamTable::Prop => "Prop",
            MmffParamTable::Pbci => "PBCI",
            MmffParamTable::Chg => "Chg",
            MmffParamTable::Bond => "Bond",
            MmffParamTable::Bndk => "Bndk",
            MmffParamTable::Stbn => "Stbn",
            MmffParamTable::Dfsb => "Dfsb",
            MmffParamTable::Angle => "Angle",
            MmffParamTable::Oop => "Oop",
            MmffParamTable::Tor => "Tor",
            MmffParamTable::Vdw => "VdW",
            MmffParamTable::HerschbachLaurie => "HerschbachLaurie",
            MmffParamTable::CovRadPauEle => "CovRadPauEle",
        };
        write!(
            formatter,
            "MMFF {table} parse error at physical line {} token {}: ",
            self.line, self.column
        )?;
        match &self.cause {
            MmffParamParseCause::InvalidUnsigned { cell } => {
                write!(formatter, "invalid unsigned integer {cell:?}")
            }
            MmffParamParseCause::InvalidFloat { cell } => {
                write!(formatter, "invalid floating point number {cell:?}")
            }
            MmffParamParseCause::MissingToken => formatter.write_str("missing token"),
            MmffParamParseCause::EmptyProcessedLine => formatter.write_str("empty processed line"),
        }
    }
}

impl std::error::Error for MmffParamParseError {}

#[cfg(test)]
mod tests {
    use super::{parse_mmff_f64, parse_mmff_u32, source_lines, tokenize_mmff_line};

    #[test]
    fn mmff_dp_m02_source_lines_match_literal_records_and_borrow_input() {
        let cases: [(&str, &[&str], &[usize]); 10] = [
            ("", &[], &[]),
            ("x", &[], &[]),
            ("x\n", &["x"], &[0]),
            ("x\r\n", &["x"], &[0]),
            ("x\r\r\n", &["x\r"], &[0]),
            ("x\ny", &["x"], &[0]),
            ("\n", &[""], &[0]),
            ("x\n\n", &["x", ""], &[0, 2]),
            ("*\nrow", &["*"], &[0]),
            ("*\nrow\n", &["*", "row"], &[0, 2]),
        ];

        let mut actual_calls = 0;
        for (input, expected, expected_offsets) in cases {
            actual_calls += 1;
            let actual: Vec<_> = source_lines(input).collect();

            assert_eq!(actual.as_slice(), expected, "input {input:?}");
            assert_eq!(actual.len(), expected_offsets.len(), "input {input:?}");
            for (record, &offset) in actual.iter().zip(expected_offsets) {
                assert_eq!(
                    record.as_ptr(),
                    input.as_ptr().wrapping_add(offset),
                    "record {record:?} did not borrow the expected input range for {input:?}"
                );
                assert!(offset + record.len() <= input.len());
            }
        }

        assert_eq!(actual_calls, 10);
    }

    #[test]
    fn mmff_dp_m03_tab_tokens_match_literal_fields_and_borrow_input() {
        let cases: [(&str, &[&str], &[usize]); 11] = [
            ("", &[], &[]),
            ("\t", &[], &[]),
            ("\t\t", &[], &[]),
            ("a", &["a"], &[0]),
            ("\ta\t", &["a"], &[1]),
            ("a\t\tb", &["a", "b"], &[0, 3]),
            ("a b\tc", &["a b", "c"], &[0, 4]),
            ("*\ta", &["*", "a"], &[0, 2]),
            (" a\t1", &[" a", "1"], &[0, 3]),
            ("a\r\t1", &["a\r", "1"], &[0, 3]),
            ("a\0\t1", &["a\0", "1"], &[0, 3]),
        ];

        let mut actual_calls = 0;
        for (input, expected, expected_offsets) in cases {
            actual_calls += 1;
            let actual: Vec<_> = tokenize_mmff_line(input).collect();

            assert_eq!(actual.as_slice(), expected, "input {input:?}");
            assert_eq!(actual.len(), expected_offsets.len(), "input {input:?}");
            for (token, &offset) in actual.iter().zip(expected_offsets) {
                assert_eq!(
                    token.as_ptr(),
                    input.as_ptr().wrapping_add(offset),
                    "token {token:?} did not borrow the expected input range for {input:?}"
                );
                assert!(offset + token.len() <= input.len());
            }
        }

        assert_eq!(actual_calls, 11);
    }

    #[test]
    fn mmff_dp_m04_unsigned_sign_wrapping_and_checked_digits() {
        const SIGNS: [&str; 3] = ["", "+", "-"];
        let cases: [(&str, [Option<u32>; 3]); 24] = [
            ("0", [Some(0), Some(0), Some(0)]),
            ("00", [Some(0), Some(0), Some(0)]),
            ("1", [Some(1), Some(1), Some(4294967295)]),
            ("255", [Some(255), Some(255), Some(4294967041)]),
            ("256", [Some(256), Some(256), Some(4294967040)]),
            ("4294967295", [Some(4294967295), Some(4294967295), Some(1)]),
            ("4294967296", [None, None, None]),
            ("18446744073709551615", [None, None, None]),
            ("", [None, None, None]),
            (" ", [None, None, None]),
            (" 1", [None, None, None]),
            ("1 ", [None, None, None]),
            ("1\t", [None, None, None]),
            ("1\n", [None, None, None]),
            ("1\r", [None, None, None]),
            ("1x", [None, None, None]),
            ("1e2", [None, None, None]),
            ("1,000", [None, None, None]),
            ("0x10", [None, None, None]),
            ("1\0", [None, None, None]),
            ("1.0", [None, None, None]),
            ("é", [None, None, None]),
            ("++1", [None, None, None]),
            ("--1", [None, None, None]),
        ];

        let mut actual_calls = 0;
        for (magnitude, expected) in cases {
            for (sign_index, sign) in SIGNS.iter().enumerate() {
                let input = format!("{sign}{magnitude}");
                assert!(actual_calls < 76, "unexpected extra parser call");
                assert_eq!(
                    parse_mmff_u32(&input),
                    expected[sign_index],
                    "input {input:?}"
                );
                actual_calls += 1;
            }
        }
        assert_eq!(actual_calls, 72);

        let zero_prefix = "0".repeat(400);
        let long_cases = [
            (format!("{zero_prefix}0"), Some(0)),
            (format!("{zero_prefix}1"), Some(1)),
            (format!("{zero_prefix}4294967295"), Some(4294967295)),
            ("9".repeat(400), None),
        ];
        for (input, expected) in long_cases {
            assert!(actual_calls < 76, "unexpected extra parser call");
            assert_eq!(
                parse_mmff_u32(&input),
                expected,
                "long input length {}",
                input.len()
            );
            actual_calls += 1;
        }

        assert_eq!(actual_calls, 76);
    }

    #[test]
    fn mmff_pc_numeric_fixed_cells_match_source_bits_and_rejections() {
        let accepted: [(&str, u64); 12] = [
            ("0", 0x0000_0000_0000_0000),
            ("-0", 0x8000_0000_0000_0000),
            ("1.25", 0x3ff4_0000_0000_0000),
            ("-.5", 0xbfe0_0000_0000_0000),
            ("+1e2", 0x4059_0000_0000_0000),
            ("1e-9999", 0x0000_0000_0000_0000),
            ("-1e-9999", 0x8000_0000_0000_0000),
            ("inf", 0x7ff0_0000_0000_0000),
            ("-Infinity", 0xfff0_0000_0000_0000),
            ("nan()", 0x7ff8_0000_0000_0000),
            ("-nan(x)", 0xfff8_0000_0000_0000),
            ("NaN(a)b)", 0x7ff8_0000_0000_0000),
        ];
        let rejected = [
            "", " 1", "1 ", "1,2", "1x", "0x1", "1e", "++1", "nan(x", "infix", "1e9999", "-1e9999",
        ];

        // parse_mmff_f64 borrows each success cell; no successful parse clones it.
        let mut actual_calls = 0;
        for (cell, expected_bits) in accepted {
            assert_eq!(
                parse_mmff_f64(cell).map(f64::to_bits),
                Some(expected_bits),
                "accepted cell {cell:?}"
            );
            actual_calls += 1;
        }
        for cell in rejected {
            assert_eq!(parse_mmff_f64(cell), None, "rejected cell {cell:?}");
            actual_calls += 1;
        }
        assert_eq!(actual_calls, 24);
    }
}
