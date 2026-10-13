//! Decimal dispatch derives from musl `__floatscan` dispatcher and `scanexp` exponent scanner (1.2.5).
//!
//! Ported from `third_party/musl/src/internal/floatscan.c:36-60` and
//! `:426-507`, together with the selected binary64 `prec=1`/`pok=1` path.

use super::cursor::ScanCursor;
use super::decimal_float::decfloat;
use super::hex_float::hexfloat;

/// `(unsigned)c - '0' < 10U` for a byte-or-EOF value, matching the source's
/// unsigned digit test.
fn is_digit(c: i32) -> bool {
    (c.wrapping_sub(i32::from(b'0')) as u32) < 10
}

/// `isspace(c)` from musl `include/ctype.h`: `c == ' ' || (unsigned)c-'\t' < 5`.
fn is_space(c: i32) -> bool {
    c == i32::from(b' ') || (c.wrapping_sub(i32::from(b'\t')) as u32) < 5
}

/// musl `floatscan.c:36-60::scanexp` (1.2.5).
///
/// Returns `i64::MIN` for an absent exponent. The `x`/`y` accumulation guards
/// bound the value well inside `i32`/`i64` (see the integration report), and the
/// final digit loop only consumes remaining digits.
pub(super) fn scanexp(cursor: &mut ScanCursor, pok: bool) -> i64 {
    // musl✔️✔️: static long long scanexp(FILE *f, int pok)
    // musl✔️✔️: {
    // musl✔️✔️: 	int c;
    // musl✔️✔️: 	int x;
    // musl✔️✔️: 	long long y;
    // musl✔️✔️: 	int neg = 0;
    // musl✔️✔️:
    // musl✔️✔️: 	c = shgetc(f);
    // musl✔️✔️: 	if (c=='+' || c=='-') {
    // musl✔️✔️: 		neg = (c=='-');
    // musl✔️✔️: 		c = shgetc(f);
    // musl✔️✔️: 		if (c-'0'>=10U && pok) shunget(f);
    // musl✔️✔️: 	}
    // musl✔️✔️: 	if (c-'0'>=10U) {
    // musl✔️✔️: 		shunget(f);
    // musl✔️✔️: 		return LLONG_MIN;
    // musl✔️✔️: 	}
    // musl✔️✔️: 	for (x=0; c-'0'<10U && x<INT_MAX/10; c = shgetc(f))
    // musl✔️✔️: 		x = 10*x + c-'0';
    // musl✔️✔️: 	for (y=x; c-'0'<10U && y<LLONG_MAX/100; c = shgetc(f))
    // musl✔️✔️: 		y = 10*y + c-'0';
    // musl✔️✔️: 	for (; c-'0'<10U; c = shgetc(f));
    // musl✔️✔️: 	shunget(f);
    // musl✔️✔️: 	return neg ? -y : y;
    // musl✔️✔️: }
    const INT_MAX: i32 = i32::MAX;
    const LLONG_MAX: i64 = i64::MAX;

    let mut c = cursor.getc();
    let mut neg = false;
    if c == i32::from(b'+') || c == i32::from(b'-') {
        neg = c == i32::from(b'-');
        c = cursor.getc();
        if !is_digit(c) && pok {
            cursor.ungetc();
        }
    }
    if !is_digit(c) {
        cursor.ungetc();
        return i64::MIN;
    }
    let mut x: i32 = 0;
    while is_digit(c) && x < INT_MAX / 10 {
        x = 10 * x + (c - i32::from(b'0'));
        c = cursor.getc();
    }
    let mut y: i64 = i64::from(x);
    while is_digit(c) && y < LLONG_MAX / 100 {
        y = 10 * y + i64::from(c - i32::from(b'0'));
        c = cursor.getc();
    }
    while is_digit(c) {
        c = cursor.getc();
    }
    cursor.ungetc();
    if neg { -y } else { y }
}

/// musl `floatscan.c:426-507::__floatscan` (1.2.5), prec dispatch and
/// sign/infinity/zero handling with source glibc NaN/hex branches.
///
/// `prec=1` selects binary64 (`bits = DBL_MANT_DIG = 53`,
/// `emin = DBL_MIN_EXP - bits = -1074`). `prec=2` selects the retained binary64
/// long-double configuration, identical to `prec=1`. `prec=0` (f32) is part of
/// the dispatch but outside this port's scope; it is not selected by the V3000
/// call sites.
pub(crate) fn float_scan(cursor: &mut ScanCursor, prec: i32, pok: bool) -> f64 {
    // musl❗✔️: long double __floatscan(FILE *f, int prec, int pok)
    // musl❗✔️: {
    // musl❗✔️: 	int sign = 1;
    // musl❗✔️: 	size_t i;
    // musl❗✔️: 	int bits;
    // musl❗✔️: 	int emin;
    // musl❗✔️: 	int c;
    // musl❗✔️:
    // musl❗✔️: 	switch (prec) {
    // musl❗✔️: 	case 0:
    // musl❗✔️: 		bits = FLT_MANT_DIG;
    // musl❗✔️: 		emin = FLT_MIN_EXP-bits;
    // musl❗✔️: 		break;
    // musl❗✔️: 	case 1:
    // musl❗✔️: 		bits = DBL_MANT_DIG;
    // musl❗✔️: 		emin = DBL_MIN_EXP-bits;
    // musl❗✔️: 		break;
    // musl❗✔️: 	case 2:
    // musl❗✔️: 		bits = LDBL_MANT_DIG;
    // musl❗✔️: 		emin = LDBL_MIN_EXP-bits;
    // musl❗✔️: 		break;
    // musl❗✔️: 	default:
    // musl❗✔️: 		return 0;
    // musl❗✔️: 	}
    // musl❗✔️:
    // musl❗✔️: 	while (isspace((c=shgetc(f))));
    // musl❗✔️:
    // musl❗✔️: 	if (c=='+' || c=='-') {
    // musl❗✔️: 		sign -= 2*(c=='-');
    // musl❗✔️: 		c = shgetc(f);
    // musl❗✔️: 	}
    // musl❗✔️:
    // musl❗✔️: 	for (i=0; i<8 && (c|32)=="infinity"[i]; i++)
    // musl❗✔️: 		if (i<7) c = shgetc(f);
    // musl❗✔️: 	if (i==3 || i==8 || (i>3 && pok)) {
    // musl❗✔️: 		if (i!=8) {
    // musl❗✔️: 			shunget(f);
    // musl❗✔️: 			if (pok) for (; i>3; i--) shunget(f);
    // musl❗✔️: 		}
    // musl❗✔️: 		return sign * INFINITY;
    // musl❗✔️: 	}
    // musl❗✔️: 	if (!i) for (i=0; i<3 && (c|32)=="nan"[i]; i++)
    // musl❗✔️: 		if (i<2) c = shgetc(f);
    // musl❗✔️: 	if (i==3) {
    // musl❗✔️: 		if (shgetc(f) != '(') {
    // musl❗✔️: 			shunget(f);
    // musl❗✔️: 			return NAN;
    // musl❗✔️: 		}
    // musl❗✔️: 		for (i=1; ; i++) {
    // musl❗✔️: 			c = shgetc(f);
    // musl❗✔️: 			if (c-'0'<10U || c-'A'<26U || c-'a'<26U || c=='_')
    // musl❗✔️: 				continue;
    // musl❗✔️: 			if (c==')') return NAN;
    // musl❗✔️: 			shunget(f);
    // musl❗✔️: 			if (!pok) {
    // musl❗✔️: 				errno = EINVAL;
    // musl❗✔️: 				shlim(f, 0);
    // musl❗✔️: 				return 0;
    // musl❗✔️: 			}
    // musl❗✔️: 			while (i--) shunget(f);
    // musl❗✔️: 			return NAN;
    // musl❗✔️: 		}
    // musl❗✔️: 		return NAN;
    // musl❗✔️: 	}
    // musl❗✔️:
    // musl❗✔️: 	if (i) {
    // musl❗✔️: 		shunget(f);
    // musl❗✔️: 		errno = EINVAL;
    // musl❗✔️: 		shlim(f, 0);
    // musl❗✔️: 		return 0;
    // musl❗✔️: 	}
    // musl❗✔️:
    // musl❗✔️: 	if (c=='0') {
    // musl❗✔️: 		c = shgetc(f);
    // musl❗✔️: 		if ((c|32) == 'x')
    // musl❗✔️: 			return hexfloat(f, bits, emin, sign, pok);
    // musl❗✔️: 		shunget(f);
    // musl❗✔️: 		c = '0';
    // musl❗✔️: 	}
    // musl❗✔️:
    // musl❗✔️: 	return decfloat(f, c, bits, emin, sign, pok);
    // musl❗✔️: }
    const INFINITY: &[u8; 8] = b"infinity";
    const NAN_WORD: &[u8; 3] = b"nan";

    let (bits, emin) = match prec {
        0 => (24i32, -125 - 24),
        1 => (53i32, -1021 - 53),
        2 => (53i32, -1021 - 53),
        _ => return 0.0,
    };

    let mut sign = 1i32;
    let mut c = cursor.getc();
    while is_space(c) {
        c = cursor.getc();
    }

    if c == i32::from(b'+') || c == i32::from(b'-') {
        if c == i32::from(b'-') {
            sign = -1;
        }
        c = cursor.getc();
    }

    let mut i = 0usize;
    while i < 8 && (c | 32) == i32::from(INFINITY[i]) {
        if i < 7 {
            c = cursor.getc();
        }
        i += 1;
    }
    if i == 3 || i == 8 || (i > 3 && pok) {
        if i != 8 {
            cursor.ungetc();
            if pok {
                while i > 3 {
                    cursor.ungetc();
                    i -= 1;
                }
            }
        }
        return (sign as f64) * f64::INFINITY;
    }

    if i == 0 {
        while i < 3 && (c | 32) == i32::from(NAN_WORD[i]) {
            if i < 2 {
                c = cursor.getc();
            }
            i += 1;
        }
    }
    if i == 3 {
        return nan_value(cursor, sign);
    }

    if i != 0 {
        cursor.ungetc();
        cursor.reset_count();
        return 0.0;
    }

    if c == i32::from(b'0') {
        c = cursor.getc();
        if (c | 32) == i32::from(b'x') {
            return hexfloat(cursor, bits, emin, sign, pok);
        }
        cursor.ungetc();
        c = i32::from(b'0');
    }

    decfloat(cursor, c, bits, emin, sign, pok)
}

/// glibc NaN branch after the dispatch consumed `nan`.
fn nan_value(cursor: &mut ScanCursor, sign: i32) -> f64 {
    // glibc✔️✔️:       if (lowc == L_('n') && STRNCASECMP (cp, L_("nan"), 3) == 0)
    // glibc✔️✔️: 	{
    // glibc✔️✔️: 	  /* Return NaN.  */
    // glibc✔️✔️: 	  FLOAT retval = NAN;
    // glibc✔️✔️:
    // glibc✔️✔️: 	  cp += 3;
    // glibc✔️✔️:
    // glibc✔️✔️: 	  /* Match `(n-char-sequence-digit)'.  */
    // glibc✔️✔️: 	  if (*cp == L_('('))
    // glibc✔️✔️: 	    {
    // glibc✔️✔️: 	      const STRING_TYPE *startp = cp;
    // glibc✔️✔️: 	      STRING_TYPE *endp;
    // glibc✔️✔️: 	      retval = STRTOF_NAN (cp + 1, &endp, L_(')'));
    // glibc✔️✔️: 	      if (*endp == L_(')'))
    // glibc✔️✔️: 		/* Consume the closing parenthesis.  */
    // glibc✔️✔️: 		cp = endp + 1;
    // glibc✔️✔️: 	      else
    // glibc✔️✔️: 		/* Only match the NAN part.  */
    // glibc✔️✔️: 		cp = startp;
    // glibc✔️✔️: 	    }
    // glibc✔️✔️:
    // glibc✔️✔️: 	  if (endptr != NULL)
    // glibc✔️✔️: 	    *endptr = (STRING_TYPE *) cp;
    // glibc✔️✔️:
    // glibc✔️✔️: 	  return negative ? -retval : retval;
    // glibc✔️✔️: 	}
    // glibc✔️✔️:

    let bytes = cursor.remaining_bytes();
    let mut value = f64::NAN;
    if bytes.first() == Some(&b'(') {
        let (payload, end) = nan_payload(&bytes[1..]);
        value = payload;
        if bytes.get(end + 1) == Some(&b')') {
            cursor.advance_bytes(end + 2);
        }
    }
    if sign < 0 { -value } else { value }
}

fn nan_payload(bytes: &[u8]) -> (f64, usize) {
    // glibc✔️✔️: FLOAT
    // glibc✔️✔️: STRTOD_NAN (const STRING_TYPE *str, STRING_TYPE **endptr, STRING_TYPE endc)
    // glibc✔️✔️: {
    // glibc✔️✔️:   const STRING_TYPE *cp = str;
    // glibc✔️✔️:
    // glibc✔️✔️:   while ((*cp >= L_('0') && *cp <= L_('9'))
    // glibc✔️✔️: 	 || (*cp >= L_('A') && *cp <= L_('Z'))
    // glibc✔️✔️: 	 || (*cp >= L_('a') && *cp <= L_('z'))
    // glibc✔️✔️: 	 || *cp == L_('_'))
    // glibc✔️✔️:     ++cp;
    // glibc✔️✔️:
    // glibc✔️✔️:   FLOAT retval = NAN;
    // glibc✔️✔️:   if (*cp != endc)
    // glibc✔️✔️:     goto out;
    // glibc✔️✔️:
    // glibc✔️✔️:   /* This is a system-dependent way to specify the bitmask used for
    // glibc✔️✔️:      the NaN.  We expect it to be a number which is put in the
    // glibc✔️✔️:      mantissa of the number.  */
    // glibc✔️✔️:   STRING_TYPE *endp;
    // glibc✔️✔️:   unsigned long long int mant;
    // glibc✔️✔️:
    // glibc✔️✔️:   int save_errno = errno;
    // glibc✔️✔️:   mant = STRTOULL (str, &endp, 0);
    // glibc✔️✔️:   __set_errno (save_errno);
    // glibc✔️✔️:   if (endp == cp)
    // glibc✔️✔️:     SET_NAN_PAYLOAD (retval, mant);
    // glibc✔️✔️:
    // glibc✔️✔️:  out:
    // glibc✔️✔️:   if (endptr != NULL)
    // glibc✔️✔️:     *endptr = (STRING_TYPE *) cp;
    // glibc✔️✔️:   return retval;
    // glibc✔️✔️: }

    let mut cp = 0usize;
    while bytes
        .get(cp)
        .is_some_and(|c| c.is_ascii_alphanumeric() || *c == b'_')
    {
        cp += 1;
    }
    let mut value = f64::NAN;
    if bytes.get(cp) == Some(&b')') {
        let (mant, end) = nan_unsigned_prefix(bytes);
        if end == cp {
            value = set_nan_payload(mant);
        }
    }
    (value, cp)
}

fn nan_unsigned_prefix(bytes: &[u8]) -> (u64, usize) {
    // glibc✔️✔️: INT
    // glibc✔️✔️: INTERNAL (__strtol_l) (const STRING_TYPE *nptr, STRING_TYPE **endptr,
    // glibc✔️✔️: 		       int base, int group, bool bin_cst, locale_t loc)
    // glibc✔️✔️: {
    // glibc✔️✔️:   int negative;
    // glibc✔️✔️:   unsigned LONG int cutoff;
    // glibc✔️✔️:   unsigned int cutlim;
    // glibc✔️✔️:   unsigned LONG int i;
    // glibc✔️✔️:   const STRING_TYPE *s;
    // glibc✔️✔️:   UCHAR_TYPE c;
    // glibc✔️✔️:   const STRING_TYPE *save, *end;
    // glibc✔️✔️:   int overflow;
    // glibc✔️✔️: #ifndef USE_WIDE_CHAR
    // glibc✔️✔️:   size_t cnt;
    // glibc✔️✔️: #endif
    // glibc✔️✔️:
    // glibc✔️✔️: #ifdef USE_NUMBER_GROUPING
    // glibc✔️✔️:   struct __locale_data *current = loc->__locales[LC_NUMERIC];
    // glibc✔️✔️:   /* The thousands character of the current locale.  */
    // glibc✔️✔️: # ifdef USE_WIDE_CHAR
    // glibc✔️✔️:   wchar_t thousands = L'\0';
    // glibc✔️✔️: # else
    // glibc✔️✔️:   const char *thousands = NULL;
    // glibc✔️✔️:   size_t thousands_len = 0;
    // glibc✔️✔️: # endif
    // glibc✔️✔️:   /* The numeric grouping specification of the current locale,
    // glibc✔️✔️:      in the format described in <locale.h>.  */
    // glibc✔️✔️:   const char *grouping;
    // glibc✔️✔️:
    // glibc✔️✔️:   if (__glibc_unlikely (group))
    // glibc✔️✔️:     {
    // glibc✔️✔️:       grouping = _NL_CURRENT (LC_NUMERIC, GROUPING);
    // glibc✔️✔️:       if (*grouping <= 0 || *grouping == CHAR_MAX)
    // glibc✔️✔️: 	grouping = NULL;
    // glibc✔️✔️:       else
    // glibc✔️✔️: 	{
    // glibc✔️✔️: 	  /* Figure out the thousands separator character.  */
    // glibc✔️✔️: # ifdef USE_WIDE_CHAR
    // glibc✔️✔️: #  ifdef _LIBC
    // glibc✔️✔️: 	  thousands = _NL_CURRENT_WORD (LC_NUMERIC,
    // glibc✔️✔️: 					_NL_NUMERIC_THOUSANDS_SEP_WC);
    // glibc✔️✔️: #  endif
    // glibc✔️✔️: 	  if (thousands == L'\0')
    // glibc✔️✔️: 	    grouping = NULL;
    // glibc✔️✔️: # else
    // glibc✔️✔️: #  ifdef _LIBC
    // glibc✔️✔️: 	  thousands = _NL_CURRENT (LC_NUMERIC, THOUSANDS_SEP);
    // glibc✔️✔️: #  endif
    // glibc✔️✔️: 	  if (*thousands == '\0')
    // glibc✔️✔️: 	    {
    // glibc✔️✔️: 	      thousands = NULL;
    // glibc✔️✔️: 	      grouping = NULL;
    // glibc✔️✔️: 	    }
    // glibc✔️✔️: # endif
    // glibc✔️✔️: 	}
    // glibc✔️✔️:     }
    // glibc✔️✔️:   else
    // glibc✔️✔️:     grouping = NULL;
    // glibc✔️✔️: #endif
    // glibc✔️✔️:
    // glibc✔️✔️:   if (base < 0 || base == 1 || base > 36)
    // glibc✔️✔️:     {
    // glibc✔️✔️:       __set_errno (EINVAL);
    // glibc✔️✔️:       return 0;
    // glibc✔️✔️:     }
    // glibc✔️✔️:
    // glibc✔️✔️:   save = s = nptr;
    // glibc✔️✔️:
    // glibc✔️✔️:   /* Skip white space.  */
    // glibc✔️✔️:   while (ISSPACE (*s))
    // glibc✔️✔️:     ++s;
    // glibc✔️✔️:   if (__glibc_unlikely (*s == L_('\0')))
    // glibc✔️✔️:     goto noconv;
    // glibc✔️✔️:
    // glibc✔️✔️:   /* Check for a sign.  */
    // glibc✔️✔️:   negative = 0;
    // glibc✔️✔️:   if (*s == L_('-'))
    // glibc✔️✔️:     {
    // glibc✔️✔️:       negative = 1;
    // glibc✔️✔️:       ++s;
    // glibc✔️✔️:     }
    // glibc✔️✔️:   else if (*s == L_('+'))
    // glibc✔️✔️:     ++s;
    // glibc✔️✔️:
    // glibc✔️✔️:   /* Recognize number prefix and if BASE is zero, figure it out ourselves.  */
    // glibc✔️✔️:   if (*s == L_('0'))
    // glibc✔️✔️:     {
    // glibc✔️✔️:       if ((base == 0 || base == 16) && TOUPPER (s[1]) == L_('X'))
    // glibc✔️✔️: 	{
    // glibc✔️✔️: 	  s += 2;
    // glibc✔️✔️: 	  base = 16;
    // glibc✔️✔️: 	}
    // glibc✔️✔️:       else if (bin_cst && (base == 0 || base == 2) && TOUPPER (s[1]) == L_('B'))
    // glibc✔️✔️: 	{
    // glibc✔️✔️: 	  s += 2;
    // glibc✔️✔️: 	  base = 2;
    // glibc✔️✔️: 	}
    // glibc✔️✔️:       else if (base == 0)
    // glibc✔️✔️: 	base = 8;
    // glibc✔️✔️:     }
    // glibc✔️✔️:   else if (base == 0)
    // glibc✔️✔️:     base = 10;
    // glibc✔️✔️:
    // glibc✔️✔️:   /* Save the pointer so we can check later if anything happened.  */
    // glibc✔️✔️:   save = s;
    // glibc✔️✔️:
    // glibc✔️✔️: #ifdef USE_NUMBER_GROUPING
    // glibc✔️✔️:   if (base != 10)
    // glibc✔️✔️:     grouping = NULL;
    // glibc✔️✔️:
    // glibc✔️✔️:   if (__glibc_unlikely (grouping != NULL))
    // glibc✔️✔️:     {
    // glibc✔️✔️: # ifndef USE_WIDE_CHAR
    // glibc✔️✔️:       thousands_len = strlen (thousands);
    // glibc✔️✔️: # endif
    // glibc✔️✔️:
    // glibc✔️✔️:       /* Find the end of the digit string and check its grouping.  */
    // glibc✔️✔️:       end = s;
    // glibc✔️✔️:       if (
    // glibc✔️✔️: # ifdef USE_WIDE_CHAR
    // glibc✔️✔️: 	  *s != thousands
    // glibc✔️✔️: # else
    // glibc✔️✔️: 	  ({ for (cnt = 0; cnt < thousands_len; ++cnt)
    // glibc✔️✔️: 	       if (thousands[cnt] != end[cnt])
    // glibc✔️✔️: 		 break;
    // glibc✔️✔️: 	     cnt < thousands_len; })
    // glibc✔️✔️: # endif
    // glibc✔️✔️: 	  )
    // glibc✔️✔️: 	{
    // glibc✔️✔️: 	  for (c = *end; c != L_('\0'); c = *++end)
    // glibc✔️✔️: 	    if (((STRING_TYPE) c < L_('0') || (STRING_TYPE) c > L_('9'))
    // glibc✔️✔️: # ifdef USE_WIDE_CHAR
    // glibc✔️✔️: 		&& (wchar_t) c != thousands
    // glibc✔️✔️: # else
    // glibc✔️✔️: 		&& ({ for (cnt = 0; cnt < thousands_len; ++cnt)
    // glibc✔️✔️: 			if (thousands[cnt] != end[cnt])
    // glibc✔️✔️: 			  break;
    // glibc✔️✔️: 		      cnt < thousands_len; })
    // glibc✔️✔️: # endif
    // glibc✔️✔️: 		&& (!ISALPHA (c)
    // glibc✔️✔️: 		    || (int) (TOUPPER (c) - L_('A') + 10) >= base))
    // glibc✔️✔️: 	      break;
    // glibc✔️✔️:
    // glibc✔️✔️: # ifdef USE_WIDE_CHAR
    // glibc✔️✔️: 	  end = __correctly_grouped_prefixwc (s, end, thousands, grouping);
    // glibc✔️✔️: # else
    // glibc✔️✔️: 	  end = __correctly_grouped_prefixmb (s, end, thousands, grouping);
    // glibc✔️✔️: # endif
    // glibc✔️✔️: 	}
    // glibc✔️✔️:     }
    // glibc✔️✔️:   else
    // glibc✔️✔️: #endif
    // glibc✔️✔️:     end = NULL;
    // glibc✔️✔️:
    // glibc✔️✔️:   /* Avoid runtime division; lookup cutoff and limit.  */
    // glibc✔️✔️:   cutoff = cutoff_tab[base - 2];
    // glibc✔️✔️:   cutlim = cutlim_tab[base - 2];
    // glibc✔️✔️:
    // glibc✔️✔️:   overflow = 0;
    // glibc✔️✔️:   i = 0;
    // glibc✔️✔️:   c = *s;
    // glibc✔️✔️:   if (sizeof (long int) != sizeof (LONG int))
    // glibc✔️✔️:     {
    // glibc✔️✔️:       unsigned long int j = 0;
    // glibc✔️✔️:       unsigned long int jmax = jmax_tab[base - 2];
    // glibc✔️✔️:
    // glibc✔️✔️:       for (;c != L_('\0'); c = *++s)
    // glibc✔️✔️: 	{
    // glibc✔️✔️: 	  if (s == end)
    // glibc✔️✔️: 	    break;
    // glibc✔️✔️: 	  if (c >= L_('0') && c <= L_('9'))
    // glibc✔️✔️: 	    c -= L_('0');
    // glibc✔️✔️: #ifdef USE_NUMBER_GROUPING
    // glibc✔️✔️: # ifdef USE_WIDE_CHAR
    // glibc✔️✔️: 	  else if (grouping && (wchar_t) c == thousands)
    // glibc✔️✔️: 	    continue;
    // glibc✔️✔️: # else
    // glibc✔️✔️: 	  else if (thousands_len)
    // glibc✔️✔️: 	    {
    // glibc✔️✔️: 	      for (cnt = 0; cnt < thousands_len; ++cnt)
    // glibc✔️✔️: 		if (thousands[cnt] != s[cnt])
    // glibc✔️✔️: 		  break;
    // glibc✔️✔️: 	      if (cnt == thousands_len)
    // glibc✔️✔️: 		{
    // glibc✔️✔️: 		  s += thousands_len - 1;
    // glibc✔️✔️: 		  continue;
    // glibc✔️✔️: 		}
    // glibc✔️✔️: 	      if (ISALPHA (c))
    // glibc✔️✔️: 		c = TOUPPER (c) - L_('A') + 10;
    // glibc✔️✔️: 	      else
    // glibc✔️✔️: 		break;
    // glibc✔️✔️: 	    }
    // glibc✔️✔️: # endif
    // glibc✔️✔️: #endif
    // glibc✔️✔️: 	  else if (ISALPHA (c))
    // glibc✔️✔️: 	    c = TOUPPER (c) - L_('A') + 10;
    // glibc✔️✔️: 	  else
    // glibc✔️✔️: 	    break;
    // glibc✔️✔️: 	  if ((int) c >= base)
    // glibc✔️✔️: 	    break;
    // glibc✔️✔️: 	  /* Note that we never can have an overflow.  */
    // glibc✔️✔️: 	  else if (j >= jmax)
    // glibc✔️✔️: 	    {
    // glibc✔️✔️: 	      /* We have an overflow.  Now use the long representation.  */
    // glibc✔️✔️: 	      i = (unsigned LONG int) j;
    // glibc✔️✔️: 	      goto use_long;
    // glibc✔️✔️: 	    }
    // glibc✔️✔️: 	  else
    // glibc✔️✔️: 	    j = j * (unsigned long int) base + c;
    // glibc✔️✔️: 	}
    // glibc✔️✔️:
    // glibc✔️✔️:       i = (unsigned LONG int) j;
    // glibc✔️✔️:     }
    // glibc✔️✔️:   else
    // glibc✔️✔️:     for (;c != L_('\0'); c = *++s)
    // glibc✔️✔️:       {
    // glibc✔️✔️: 	if (s == end)
    // glibc✔️✔️: 	  break;
    // glibc✔️✔️: 	if (c >= L_('0') && c <= L_('9'))
    // glibc✔️✔️: 	  c -= L_('0');
    // glibc✔️✔️: #ifdef USE_NUMBER_GROUPING
    // glibc✔️✔️: # ifdef USE_WIDE_CHAR
    // glibc✔️✔️: 	else if (grouping && (wchar_t) c == thousands)
    // glibc✔️✔️: 	  continue;
    // glibc✔️✔️: # else
    // glibc✔️✔️: 	else if (thousands_len)
    // glibc✔️✔️: 	  {
    // glibc✔️✔️: 	    for (cnt = 0; cnt < thousands_len; ++cnt)
    // glibc✔️✔️: 	      if (thousands[cnt] != s[cnt])
    // glibc✔️✔️: 		break;
    // glibc✔️✔️: 	    if (cnt == thousands_len)
    // glibc✔️✔️: 	      {
    // glibc✔️✔️: 		s += thousands_len - 1;
    // glibc✔️✔️: 		continue;
    // glibc✔️✔️: 	      }
    // glibc✔️✔️: 	    if (ISALPHA (c))
    // glibc✔️✔️: 	      c = TOUPPER (c) - L_('A') + 10;
    // glibc✔️✔️: 	    else
    // glibc✔️✔️: 	      break;
    // glibc✔️✔️: 	  }
    // glibc✔️✔️: # endif
    // glibc✔️✔️: #endif
    // glibc✔️✔️: 	else if (ISALPHA (c))
    // glibc✔️✔️: 	  c = TOUPPER (c) - L_('A') + 10;
    // glibc✔️✔️: 	else
    // glibc✔️✔️: 	  break;
    // glibc✔️✔️: 	if ((int) c >= base)
    // glibc✔️✔️: 	  break;
    // glibc✔️✔️: 	/* Check for overflow.  */
    // glibc✔️✔️: 	if (i > cutoff || (i == cutoff && c > cutlim))
    // glibc✔️✔️: 	  overflow = 1;
    // glibc✔️✔️: 	else
    // glibc✔️✔️: 	  {
    // glibc✔️✔️: 	  use_long:
    // glibc✔️✔️: 	    i *= (unsigned LONG int) base;
    // glibc✔️✔️: 	    i += c;
    // glibc✔️✔️: 	  }
    // glibc✔️✔️:       }
    // glibc✔️✔️:
    // glibc✔️✔️:   /* Check if anything actually happened.  */
    // glibc✔️✔️:   if (s == save)
    // glibc✔️✔️:     goto noconv;
    // glibc✔️✔️:
    // glibc✔️✔️:   /* Store in ENDPTR the address of one character
    // glibc✔️✔️:      past the last character we converted.  */
    // glibc✔️✔️:   if (endptr != NULL)
    // glibc✔️✔️:     *endptr = (STRING_TYPE *) s;
    // glibc✔️✔️:
    // glibc✔️✔️: #if !UNSIGNED
    // glibc✔️✔️:   /* Check for a value that is within the range of
    // glibc✔️✔️:      `unsigned LONG int', but outside the range of `LONG int'.  */
    // glibc✔️✔️:   if (overflow == 0
    // glibc✔️✔️:       && i > (negative
    // glibc✔️✔️: 	      ? -((unsigned LONG int) (STRTOL_LONG_MIN + 1)) + 1
    // glibc✔️✔️: 	      : (unsigned LONG int) STRTOL_LONG_MAX))
    // glibc✔️✔️:     overflow = 1;
    // glibc✔️✔️: #endif
    // glibc✔️✔️:
    // glibc✔️✔️:   if (__glibc_unlikely (overflow))
    // glibc✔️✔️:     {
    // glibc✔️✔️:       __set_errno (ERANGE);
    // glibc✔️✔️: #if UNSIGNED
    // glibc✔️✔️:       return STRTOL_ULONG_MAX;
    // glibc✔️✔️: #else
    // glibc✔️✔️:       return negative ? STRTOL_LONG_MIN : STRTOL_LONG_MAX;
    // glibc✔️✔️: #endif
    // glibc✔️✔️:     }
    // glibc✔️✔️:
    // glibc✔️✔️:   /* Return the result of the appropriate sign.  */
    // glibc✔️✔️:   return negative ? -i : i;
    // glibc✔️✔️:
    // glibc✔️✔️: noconv:
    // glibc✔️✔️:   /* We must handle a special case here: the base is 0 or 16 and the
    // glibc✔️✔️:      first two characters are '0' and 'x', but the rest are no
    // glibc✔️✔️:      hexadecimal digits.  Likewise when the base is 0 or 2 and the
    // glibc✔️✔️:      first two characters are '0' and 'b', but the rest are no binary
    // glibc✔️✔️:      digits.  This is no error case.  We return 0 and ENDPTR points to
    // glibc✔️✔️:      the 'x' or 'b'.  */
    // glibc✔️✔️:   if (endptr != NULL)
    // glibc✔️✔️:     {
    // glibc✔️✔️:       if (save - nptr >= 2
    // glibc✔️✔️: 	  && (TOUPPER (save[-1]) == L_('X')
    // glibc✔️✔️: 	      || (bin_cst && TOUPPER (save[-1]) == L_('B')))
    // glibc✔️✔️: 	  && save[-2] == L_('0'))
    // glibc✔️✔️: 	*endptr = (STRING_TYPE *) &save[-1];
    // glibc✔️✔️:       else
    // glibc✔️✔️: 	/*  There was no number to convert.  */
    // glibc✔️✔️: 	*endptr = (STRING_TYPE *) nptr;
    // glibc✔️✔️:     }
    // glibc✔️✔️:
    // glibc✔️✔️:   return 0L;
    // glibc✔️✔️: }

    // Source reachable specialization: C ASCII n-char sequence means no
    // whitespace/sign, base=0, group=0, bin_cst=false, unsigned 64-bit.
    // Once overflow is observed it stays set while all legal digits consume.
    let mut s = 0usize;
    let base = if bytes.first() == Some(&b'0') {
        if bytes.get(1).is_some_and(|c| c.eq_ignore_ascii_case(&b'x')) {
            s = 2;
            16u64
        } else {
            8u64
        }
    } else {
        10u64
    };
    let save = s;
    let cutoff = u64::MAX / base;
    let cutlim = u64::MAX % base;
    let mut overflow = false;
    let mut value = 0u64;
    while let Some(&c) = bytes.get(s) {
        let d = match c {
            b'0'..=b'9' => u64::from(c - b'0'),
            b'A'..=b'Z' => u64::from(c - b'A' + 10),
            b'a'..=b'z' => u64::from(c - b'a' + 10),
            _ => break,
        };
        if d >= base {
            break;
        }
        if value > cutoff || (value == cutoff && d > cutlim) {
            overflow = true;
        } else {
            value = value * base + d;
        }
        s += 1;
    }
    if s == save {
        return (0, if save == 2 { 1 } else { 0 });
    }
    (if overflow { u64::MAX } else { value }, s)
}

fn set_nan_payload(mant: u64) -> f64 {
    // glibc✔️✔️: #define SET_NAN_PAYLOAD(flt, mant)			\
    // glibc✔️✔️:   do							\
    // glibc✔️✔️:     {							\
    // glibc✔️✔️:       union ieee754_double u;				\
    // glibc✔️✔️:       u.d = (flt);					\
    // glibc✔️✔️:       u.ieee_nan.mantissa0 = (mant) >> 32;		\
    // glibc✔️✔️:       u.ieee_nan.mantissa1 = (mant);			\
    // glibc✔️✔️:       if ((u.ieee.mantissa0 | u.ieee.mantissa1) != 0)	\
    // glibc✔️✔️: 	(flt) = u.d;					\
    // glibc✔️✔️:     }							\
    // glibc✔️✔️:   while (0)

    // ieee_nan has a separate quiet bit; assigning its 19+32 mantissa
    // bitfields truncates the unsigned payload to 51 bits and retains NAN's
    // quiet bit. Cost: constant integer state, no allocation or cloning.
    f64::from_bits(0x7ff8_0000_0000_0000 | (mant & 0x0007_ffff_ffff_ffff))
}

#[cfg(test)]
mod tests {
    use super::{is_digit, is_space, scanexp};
    use crate::numeric::cursor::ScanCursor;

    #[test]
    fn atof_port_scanexp_digit_classification() {
        assert!(is_digit(i32::from(b'0')));
        assert!(is_digit(i32::from(b'9')));
        assert!(!is_digit(i32::from(b'e')));
        assert!(!is_digit(0));
        assert!(!is_digit(-1));
        assert!(is_space(i32::from(b' ')));
        assert!(is_space(i32::from(b'\t')));
        assert!(is_space(i32::from(b'\n')));
        assert!(!is_space(0));
        assert!(!is_space(i32::from(b'0')));
    }

    fn scan(bytes: &[u8]) -> (i64, usize) {
        let mut cursor = ScanCursor::new(bytes);
        let value = scanexp(&mut cursor, true);
        (value, cursor.consumed())
    }

    #[test]
    fn atof_port_scanexp_signs_and_missing_digits() {
        assert_eq!(scan(b"+12"), (12, 3));
        assert_eq!(scan(b"-12"), (-12, 3));
        assert_eq!(scan(b"0"), (0, 1));
        // Missing digits: the exponent is absent and the cursor rewinds.
        assert_eq!(scan(b"e"), (i64::MIN, 0));
        assert_eq!(scan(b"+"), (i64::MIN, 0));
    }

    #[test]
    fn atof_port_scanexp_long_runs_bound_without_saturation() {
        // Fifty 9s far exceed i64; the source guard stops accumulating and only
        // consumes the remaining digits, so the value stays below LLONG_MAX
        // (about 9.22e17) instead of saturating.
        let token = [b'9'; 50];
        let (value, offset) = scan(&token);
        assert_eq!(offset, 50);
        assert!(value > 0);
        assert!(value < i64::MAX);
        assert_eq!(value, 99999999999999999);
    }

    #[test]
    fn atof_port_scanexp_thousand_digit_exponent_consumes_all() {
        let token = [b'9'; 1000];
        let (value, offset) = scan(&token);
        assert_eq!(offset, 1000);
        assert_eq!(value, 99999999999999999);
    }
}
