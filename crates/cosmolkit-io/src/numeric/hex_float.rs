//! glibc 2.43 hexadecimal binary64 conversion, C locale, FE_TONEAREST.
//! Source: stdlib/strtod_l.c and sysdeps/ieee754/dbl-64/mpn2dbl.c.
use super::cursor::ScanCursor;

fn digit(byte: u8) -> Option<u64> {
    match byte {
        b'0'..=b'9' => Some(u64::from(byte - b'0')),
        b'a'..=b'f' => Some(u64::from(byte - b'a' + 10)),
        b'A'..=b'F' => Some(u64::from(byte - b'A' + 10)),
        _ => None,
    }
}

/// Complete base-16 branch after the dispatch consumed the sign and `0x`.
/// The full corresponding upstream function is retained inline, including
/// its unreachable decimal/group/wide-character branches. This function
/// specializes base=16, group=0, narrow C locale, 64-bit limb, binary64.
/// It scans/normalizes borrowed bytes then packs integer mantissa/round/sticky
/// bits. Cost: linear bounded passes, constant state, no allocation or clone.
pub(super) fn hexfloat(cursor: &mut ScanCursor, bits: i32, emin: i32, sign: i32, pok: bool) -> f64 {
    // glibc❗✔️: FLOAT
    // glibc❗✔️: ____STRTOF_INTERNAL (const STRING_TYPE *nptr, STRING_TYPE **endptr, int group,
    // glibc❗✔️: 		     locale_t loc)
    // glibc❗✔️: {
    // glibc❗✔️:   int negative;			/* The sign of the number.  */
    // glibc❗✔️:   MPN_VAR (num);		/* MP representation of the number.  */
    // glibc❗✔️:   intmax_t exponent;		/* Exponent of the number.  */
    // glibc❗✔️:
    // glibc❗✔️:   /* Numbers starting `0X' or `0x' have to be processed with base 16.  */
    // glibc❗✔️:   int base = 10;
    // glibc❗✔️:
    // glibc❗✔️:   /* When we have to compute fractional digits we form a fraction with a
    // glibc❗✔️:      second multi-precision number (and we sometimes need a second for
    // glibc❗✔️:      temporary results).  */
    // glibc❗✔️:   MPN_VAR (den);
    // glibc❗✔️:
    // glibc❗✔️:   /* Representation for the return value.  */
    // glibc❗✔️:   mp_limb_t retval[RETURN_LIMB_SIZE];
    // glibc❗✔️:   /* Number of bits currently in result value.  */
    // glibc❗✔️:   int bits;
    // glibc❗✔️:
    // glibc❗✔️:   /* Running pointer after the last character processed in the string.  */
    // glibc❗✔️:   const STRING_TYPE *cp, *tp;
    // glibc❗✔️:   /* Start of significant part of the number.  */
    // glibc❗✔️:   const STRING_TYPE *startp, *start_of_digits;
    // glibc❗✔️:   /* Points at the character following the integer and fractional digits.  */
    // glibc❗✔️:   const STRING_TYPE *expp;
    // glibc❗✔️:   /* Total number of digit and number of digits in integer part.  */
    // glibc❗✔️:   size_t dig_no, int_no, lead_zero;
    // glibc❗✔️:   /* Contains the last character read.  */
    // glibc❗✔️:   CHAR_TYPE c;
    // glibc❗✔️:
    // glibc❗✔️: /* We should get wint_t from <stddef.h>, but not all GCC versions define it
    // glibc❗✔️:    there.  So define it ourselves if it remains undefined.  */
    // glibc❗✔️: #ifndef _WINT_T
    // glibc❗✔️:   typedef unsigned int wint_t;
    // glibc❗✔️: #endif
    // glibc❗✔️:   /* The radix character of the current locale.  */
    // glibc❗✔️: #ifdef USE_WIDE_CHAR
    // glibc❗✔️:   wchar_t decimal;
    // glibc❗✔️: #else
    // glibc❗✔️:   const char *decimal;
    // glibc❗✔️:   size_t decimal_len;
    // glibc❗✔️: #endif
    // glibc❗✔️:   /* The thousands character of the current locale.  */
    // glibc❗✔️: #ifdef USE_WIDE_CHAR
    // glibc❗✔️:   wchar_t thousands = L'\0';
    // glibc❗✔️: #else
    // glibc❗✔️:   const char *thousands = NULL;
    // glibc❗✔️: #endif
    // glibc❗✔️:   /* The numeric grouping specification of the current locale,
    // glibc❗✔️:      in the format described in <locale.h>.  */
    // glibc❗✔️:   const char *grouping;
    // glibc❗✔️:   /* Used in several places.  */
    // glibc❗✔️:   int cnt;
    // glibc❗✔️:
    // glibc❗✔️:   struct __locale_data *current = loc->__locales[LC_NUMERIC];
    // glibc❗✔️:
    // glibc❗✔️:   if (__glibc_unlikely (group))
    // glibc❗✔️:     {
    // glibc❗✔️:       grouping = _NL_CURRENT (LC_NUMERIC, GROUPING);
    // glibc❗✔️:       if (*grouping <= 0 || *grouping == CHAR_MAX)
    // glibc❗✔️: 	grouping = NULL;
    // glibc❗✔️:       else
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  /* Figure out the thousands separator character.  */
    // glibc❗✔️: #ifdef USE_WIDE_CHAR
    // glibc❗✔️: 	  thousands = _NL_CURRENT_WORD (LC_NUMERIC,
    // glibc❗✔️: 					_NL_NUMERIC_THOUSANDS_SEP_WC);
    // glibc❗✔️: 	  if (thousands == L'\0')
    // glibc❗✔️: 	    grouping = NULL;
    // glibc❗✔️: #else
    // glibc❗✔️: 	  thousands = _NL_CURRENT (LC_NUMERIC, THOUSANDS_SEP);
    // glibc❗✔️: 	  if (*thousands == '\0')
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      thousands = NULL;
    // glibc❗✔️: 	      grouping = NULL;
    // glibc❗✔️: 	    }
    // glibc❗✔️: #endif
    // glibc❗✔️: 	}
    // glibc❗✔️:     }
    // glibc❗✔️:   else
    // glibc❗✔️:     grouping = NULL;
    // glibc❗✔️:
    // glibc❗✔️:   /* Find the locale's decimal point character.  */
    // glibc❗✔️: #ifdef USE_WIDE_CHAR
    // glibc❗✔️:   decimal = _NL_CURRENT_WORD (LC_NUMERIC, _NL_NUMERIC_DECIMAL_POINT_WC);
    // glibc❗✔️:   assert (decimal != L'\0');
    // glibc❗✔️: # define decimal_len 1
    // glibc❗✔️: #else
    // glibc❗✔️:   decimal = _NL_CURRENT (LC_NUMERIC, DECIMAL_POINT);
    // glibc❗✔️:   decimal_len = strlen (decimal);
    // glibc❗✔️:   assert (decimal_len > 0);
    // glibc❗✔️: #endif
    // glibc❗✔️:
    // glibc❗✔️:   /* Prepare number representation.  */
    // glibc❗✔️:   exponent = 0;
    // glibc❗✔️:   negative = 0;
    // glibc❗✔️:   bits = 0;
    // glibc❗✔️:
    // glibc❗✔️:   /* Parse string to get maximal legal prefix.  We need the number of
    // glibc❗✔️:      characters of the integer part, the fractional part and the exponent.  */
    // glibc❗✔️:   cp = nptr - 1;
    // glibc❗✔️:   /* Ignore leading white space.  */
    // glibc❗✔️:   do
    // glibc❗✔️:     c = *++cp;
    // glibc❗✔️:   while (ISSPACE (c));
    // glibc❗✔️:
    // glibc❗✔️:   /* Get sign of the result.  */
    // glibc❗✔️:   if (c == L_('-'))
    // glibc❗✔️:     {
    // glibc❗✔️:       negative = 1;
    // glibc❗✔️:       c = *++cp;
    // glibc❗✔️:     }
    // glibc❗✔️:   else if (c == L_('+'))
    // glibc❗✔️:     c = *++cp;
    // glibc❗✔️:
    // glibc❗✔️:   /* Return 0.0 if no legal string is found.
    // glibc❗✔️:      No character is used even if a sign was found.  */
    // glibc❗✔️: #ifdef USE_WIDE_CHAR
    // glibc❗✔️:   if (c == (wint_t) decimal
    // glibc❗✔️:       && (wint_t) cp[1] >= L'0' && (wint_t) cp[1] <= L'9')
    // glibc❗✔️:     {
    // glibc❗✔️:       /* We accept it.  This funny construct is here only to indent
    // glibc❗✔️: 	 the code correctly.  */
    // glibc❗✔️:     }
    // glibc❗✔️: #else
    // glibc❗✔️:   for (cnt = 0; decimal[cnt] != '\0'; ++cnt)
    // glibc❗✔️:     if (cp[cnt] != decimal[cnt])
    // glibc❗✔️:       break;
    // glibc❗✔️:   if (decimal[cnt] == '\0' && cp[cnt] >= '0' && cp[cnt] <= '9')
    // glibc❗✔️:     {
    // glibc❗✔️:       /* We accept it.  This funny construct is here only to indent
    // glibc❗✔️: 	 the code correctly.  */
    // glibc❗✔️:     }
    // glibc❗✔️: #endif
    // glibc❗✔️:   else if (c < L_('0') || c > L_('9'))
    // glibc❗✔️:     {
    // glibc❗✔️:       /* Check for `INF' or `INFINITY'.  */
    // glibc❗✔️:       CHAR_TYPE lowc = TOLOWER_C (c);
    // glibc❗✔️:
    // glibc❗✔️:       if (lowc == L_('i') && STRNCASECMP (cp, L_("inf"), 3) == 0)
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  /* Return +/- infinity.  */
    // glibc❗✔️: 	  if (endptr != NULL)
    // glibc❗✔️: 	    *endptr = (STRING_TYPE *)
    // glibc❗✔️: 		      (cp + (STRNCASECMP (cp + 3, L_("inity"), 5) == 0
    // glibc❗✔️: 			     ? 8 : 3));
    // glibc❗✔️:
    // glibc❗✔️: 	  return negative ? -FLOAT_HUGE_VAL : FLOAT_HUGE_VAL;
    // glibc❗✔️: 	}
    // glibc❗✔️:
    // glibc❗✔️:       if (lowc == L_('n') && STRNCASECMP (cp, L_("nan"), 3) == 0)
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  /* Return NaN.  */
    // glibc❗✔️: 	  FLOAT retval = NAN;
    // glibc❗✔️:
    // glibc❗✔️: 	  cp += 3;
    // glibc❗✔️:
    // glibc❗✔️: 	  /* Match `(n-char-sequence-digit)'.  */
    // glibc❗✔️: 	  if (*cp == L_('('))
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      const STRING_TYPE *startp = cp;
    // glibc❗✔️: 	      STRING_TYPE *endp;
    // glibc❗✔️: 	      retval = STRTOF_NAN (cp + 1, &endp, L_(')'));
    // glibc❗✔️: 	      if (*endp == L_(')'))
    // glibc❗✔️: 		/* Consume the closing parenthesis.  */
    // glibc❗✔️: 		cp = endp + 1;
    // glibc❗✔️: 	      else
    // glibc❗✔️: 		/* Only match the NAN part.  */
    // glibc❗✔️: 		cp = startp;
    // glibc❗✔️: 	    }
    // glibc❗✔️:
    // glibc❗✔️: 	  if (endptr != NULL)
    // glibc❗✔️: 	    *endptr = (STRING_TYPE *) cp;
    // glibc❗✔️:
    // glibc❗✔️: 	  return negative ? -retval : retval;
    // glibc❗✔️: 	}
    // glibc❗✔️:
    // glibc❗✔️:       /* It is really a text we do not recognize.  */
    // glibc❗✔️:       RETURN (0.0, nptr);
    // glibc❗✔️:     }
    // glibc❗✔️:
    // glibc❗✔️:   /* First look whether we are faced with a hexadecimal number.  */
    // glibc❗✔️:   if (c == L_('0') && TOLOWER (cp[1]) == L_('x'))
    // glibc❗✔️:     {
    // glibc❗✔️:       /* Okay, it is a hexa-decimal number.  Remember this and skip
    // glibc❗✔️: 	 the characters.  BTW: hexadecimal numbers must not be
    // glibc❗✔️: 	 grouped.  */
    // glibc❗✔️:       base = 16;
    // glibc❗✔️:       cp += 2;
    // glibc❗✔️:       c = *cp;
    // glibc❗✔️:       grouping = NULL;
    // glibc❗✔️:     }
    // glibc❗✔️:
    // glibc❗✔️:   /* Record the start of the digits, in case we will check their grouping.  */
    // glibc❗✔️:   start_of_digits = startp = cp;
    // glibc❗✔️:
    // glibc❗✔️:   /* Ignore leading zeroes.  This helps us to avoid useless computations.  */
    // glibc❗✔️: #ifdef USE_WIDE_CHAR
    // glibc❗✔️:   while (c == L'0' || ((wint_t) thousands != L'\0' && c == (wint_t) thousands))
    // glibc❗✔️:     c = *++cp;
    // glibc❗✔️: #else
    // glibc❗✔️:   if (__glibc_likely (thousands == NULL))
    // glibc❗✔️:     while (c == '0')
    // glibc❗✔️:       c = *++cp;
    // glibc❗✔️:   else
    // glibc❗✔️:     {
    // glibc❗✔️:       /* We also have the multibyte thousands string.  */
    // glibc❗✔️:       while (1)
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  if (c != '0')
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      for (cnt = 0; thousands[cnt] != '\0'; ++cnt)
    // glibc❗✔️: 		if (thousands[cnt] != cp[cnt])
    // glibc❗✔️: 		  break;
    // glibc❗✔️: 	      if (thousands[cnt] != '\0')
    // glibc❗✔️: 		break;
    // glibc❗✔️: 	      cp += cnt - 1;
    // glibc❗✔️: 	    }
    // glibc❗✔️: 	  c = *++cp;
    // glibc❗✔️: 	}
    // glibc❗✔️:     }
    // glibc❗✔️: #endif
    // glibc❗✔️:
    // glibc❗✔️:   /* If no other digit but a '0' is found the result is 0.0.
    // glibc❗✔️:      Return current read pointer.  */
    // glibc❗✔️:   CHAR_TYPE lowc = TOLOWER (c);
    // glibc❗✔️:   if (!((c >= L_('0') && c <= L_('9'))
    // glibc❗✔️: 	|| (base == 16 && lowc >= L_('a') && lowc <= L_('f'))
    // glibc❗✔️: 	|| (
    // glibc❗✔️: #ifdef USE_WIDE_CHAR
    // glibc❗✔️: 	    c == (wint_t) decimal
    // glibc❗✔️: #else
    // glibc❗✔️: 	    ({ for (cnt = 0; decimal[cnt] != '\0'; ++cnt)
    // glibc❗✔️: 		 if (decimal[cnt] != cp[cnt])
    // glibc❗✔️: 		   break;
    // glibc❗✔️: 	       decimal[cnt] == '\0'; })
    // glibc❗✔️: #endif
    // glibc❗✔️: 	    /* '0x.' alone is not a valid hexadecimal number.
    // glibc❗✔️: 	       '.' alone is not valid either, but that has been checked
    // glibc❗✔️: 	       already earlier.  */
    // glibc❗✔️: 	    && (base != 16
    // glibc❗✔️: 		|| cp != start_of_digits
    // glibc❗✔️: 		|| (cp[decimal_len] >= L_('0') && cp[decimal_len] <= L_('9'))
    // glibc❗✔️: 		|| ({ CHAR_TYPE lo = TOLOWER (cp[decimal_len]);
    // glibc❗✔️: 		      lo >= L_('a') && lo <= L_('f'); })))
    // glibc❗✔️: 	|| (base == 16 && (cp != start_of_digits
    // glibc❗✔️: 			   && lowc == L_('p')))
    // glibc❗✔️: 	|| (base != 16 && lowc == L_('e'))))
    // glibc❗✔️:     {
    // glibc❗✔️: #ifdef USE_WIDE_CHAR
    // glibc❗✔️:       tp = __correctly_grouped_prefixwc (start_of_digits, cp, thousands,
    // glibc❗✔️: 					 grouping);
    // glibc❗✔️: #else
    // glibc❗✔️:       tp = __correctly_grouped_prefixmb (start_of_digits, cp, thousands,
    // glibc❗✔️: 					 grouping);
    // glibc❗✔️: #endif
    // glibc❗✔️:       /* If TP is at the start of the digits, there was no correctly
    // glibc❗✔️: 	 grouped prefix of the string; so no number found.  */
    // glibc❗✔️:       RETURN (negative ? -0.0 : 0.0,
    // glibc❗✔️: 	      tp == start_of_digits ? (base == 16 ? cp - 1 : nptr) : tp);
    // glibc❗✔️:     }
    // glibc❗✔️:
    // glibc❗✔️:   /* Remember first significant digit and read following characters until the
    // glibc❗✔️:      decimal point, exponent character or any non-FP number character.  */
    // glibc❗✔️:   startp = cp;
    // glibc❗✔️:   dig_no = 0;
    // glibc❗✔️:   while (1)
    // glibc❗✔️:     {
    // glibc❗✔️:       if ((c >= L_('0') && c <= L_('9'))
    // glibc❗✔️: 	  || (base == 16
    // glibc❗✔️: 	      && ({ CHAR_TYPE lo = TOLOWER (c);
    // glibc❗✔️: 		    lo >= L_('a') && lo <= L_('f'); })))
    // glibc❗✔️: 	++dig_no;
    // glibc❗✔️:       else
    // glibc❗✔️: 	{
    // glibc❗✔️: #ifdef USE_WIDE_CHAR
    // glibc❗✔️: 	  if (__builtin_expect ((wint_t) thousands == L'\0', 1)
    // glibc❗✔️: 	      || c != (wint_t) thousands)
    // glibc❗✔️: 	    /* Not a digit or separator: end of the integer part.  */
    // glibc❗✔️: 	    break;
    // glibc❗✔️: #else
    // glibc❗✔️: 	  if (__glibc_likely (thousands == NULL))
    // glibc❗✔️: 	    break;
    // glibc❗✔️: 	  else
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      for (cnt = 0; thousands[cnt] != '\0'; ++cnt)
    // glibc❗✔️: 		if (thousands[cnt] != cp[cnt])
    // glibc❗✔️: 		  break;
    // glibc❗✔️: 	      if (thousands[cnt] != '\0')
    // glibc❗✔️: 		break;
    // glibc❗✔️: 	      cp += cnt - 1;
    // glibc❗✔️: 	    }
    // glibc❗✔️: #endif
    // glibc❗✔️: 	}
    // glibc❗✔️:       c = *++cp;
    // glibc❗✔️:     }
    // glibc❗✔️:
    // glibc❗✔️:   if (__builtin_expect (grouping != NULL, 0) && cp > start_of_digits)
    // glibc❗✔️:     {
    // glibc❗✔️:       /* Check the grouping of the digits.  */
    // glibc❗✔️: #ifdef USE_WIDE_CHAR
    // glibc❗✔️:       tp = __correctly_grouped_prefixwc (start_of_digits, cp, thousands,
    // glibc❗✔️: 					 grouping);
    // glibc❗✔️: #else
    // glibc❗✔️:       tp = __correctly_grouped_prefixmb (start_of_digits, cp, thousands,
    // glibc❗✔️: 					 grouping);
    // glibc❗✔️: #endif
    // glibc❗✔️:       if (cp != tp)
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  /* Less than the entire string was correctly grouped.  */
    // glibc❗✔️:
    // glibc❗✔️: 	  if (tp == start_of_digits)
    // glibc❗✔️: 	    /* No valid group of numbers at all: no valid number.  */
    // glibc❗✔️: 	    RETURN (0.0, nptr);
    // glibc❗✔️:
    // glibc❗✔️: 	  if (tp < startp)
    // glibc❗✔️: 	    /* The number is validly grouped, but consists
    // glibc❗✔️: 	       only of zeroes.  The whole value is zero.  */
    // glibc❗✔️: 	    RETURN (negative ? -0.0 : 0.0, tp);
    // glibc❗✔️:
    // glibc❗✔️: 	  /* Recompute DIG_NO so we won't read more digits than
    // glibc❗✔️: 	     are properly grouped.  */
    // glibc❗✔️: 	  cp = tp;
    // glibc❗✔️: 	  dig_no = 0;
    // glibc❗✔️: 	  for (tp = startp; tp < cp; ++tp)
    // glibc❗✔️: 	    if (*tp >= L_('0') && *tp <= L_('9'))
    // glibc❗✔️: 	      ++dig_no;
    // glibc❗✔️:
    // glibc❗✔️: 	  int_no = dig_no;
    // glibc❗✔️: 	  lead_zero = 0;
    // glibc❗✔️:
    // glibc❗✔️: 	  goto number_parsed;
    // glibc❗✔️: 	}
    // glibc❗✔️:     }
    // glibc❗✔️:
    // glibc❗✔️:   /* We have the number of digits in the integer part.  Whether these
    // glibc❗✔️:      are all or any is really a fractional digit will be decided
    // glibc❗✔️:      later.  */
    // glibc❗✔️:   int_no = dig_no;
    // glibc❗✔️:   lead_zero = int_no == 0 ? (size_t) -1 : 0;
    // glibc❗✔️:
    // glibc❗✔️:   /* Read the fractional digits.  A special case are the 'american
    // glibc❗✔️:      style' numbers like `16.' i.e. with decimal point but without
    // glibc❗✔️:      trailing digits.  */
    // glibc❗✔️:   if (
    // glibc❗✔️: #ifdef USE_WIDE_CHAR
    // glibc❗✔️:       c == (wint_t) decimal
    // glibc❗✔️: #else
    // glibc❗✔️:       ({ for (cnt = 0; decimal[cnt] != '\0'; ++cnt)
    // glibc❗✔️: 	   if (decimal[cnt] != cp[cnt])
    // glibc❗✔️: 	     break;
    // glibc❗✔️: 	 decimal[cnt] == '\0'; })
    // glibc❗✔️: #endif
    // glibc❗✔️:       )
    // glibc❗✔️:     {
    // glibc❗✔️:       cp += decimal_len;
    // glibc❗✔️:       c = *cp;
    // glibc❗✔️:       while ((c >= L_('0') && c <= L_('9'))
    // glibc❗✔️: 	     || (base == 16 && ({ CHAR_TYPE lo = TOLOWER (c);
    // glibc❗✔️: 				  lo >= L_('a') && lo <= L_('f'); })))
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  if (c != L_('0') && lead_zero == (size_t) -1)
    // glibc❗✔️: 	    lead_zero = dig_no - int_no;
    // glibc❗✔️: 	  ++dig_no;
    // glibc❗✔️: 	  c = *++cp;
    // glibc❗✔️: 	}
    // glibc❗✔️:     }
    // glibc❗✔️:   assert (dig_no <= (uintmax_t) INTMAX_MAX);
    // glibc❗✔️:
    // glibc❗✔️:   /* Remember start of exponent (if any).  */
    // glibc❗✔️:   expp = cp;
    // glibc❗✔️:
    // glibc❗✔️:   /* Read exponent.  */
    // glibc❗✔️:   lowc = TOLOWER (c);
    // glibc❗✔️:   if ((base == 16 && lowc == L_('p'))
    // glibc❗✔️:       || (base != 16 && lowc == L_('e')))
    // glibc❗✔️:     {
    // glibc❗✔️:       int exp_negative = 0;
    // glibc❗✔️:
    // glibc❗✔️:       c = *++cp;
    // glibc❗✔️:       if (c == L_('-'))
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  exp_negative = 1;
    // glibc❗✔️: 	  c = *++cp;
    // glibc❗✔️: 	}
    // glibc❗✔️:       else if (c == L_('+'))
    // glibc❗✔️: 	c = *++cp;
    // glibc❗✔️:
    // glibc❗✔️:       if (c >= L_('0') && c <= L_('9'))
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  intmax_t exp_limit;
    // glibc❗✔️:
    // glibc❗✔️: 	  /* Get the exponent limit. */
    // glibc❗✔️: 	  if (base == 16)
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      if (exp_negative)
    // glibc❗✔️: 		{
    // glibc❗✔️: 		  assert (int_no <= (uintmax_t) (INTMAX_MAX
    // glibc❗✔️: 						 + MIN_EXP - MANT_DIG) / 4);
    // glibc❗✔️: 		  exp_limit = -MIN_EXP + MANT_DIG + 4 * (intmax_t) int_no;
    // glibc❗✔️: 		}
    // glibc❗✔️: 	      else
    // glibc❗✔️: 		{
    // glibc❗✔️: 		  if (int_no)
    // glibc❗✔️: 		    {
    // glibc❗✔️: 		      assert (lead_zero == 0
    // glibc❗✔️: 			      && int_no <= (uintmax_t) INTMAX_MAX / 4);
    // glibc❗✔️: 		      exp_limit = MAX_EXP - 4 * (intmax_t) int_no + 3;
    // glibc❗✔️: 		    }
    // glibc❗✔️: 		  else if (lead_zero == (size_t) -1)
    // glibc❗✔️: 		    {
    // glibc❗✔️: 		      /* The number is zero and this limit is
    // glibc❗✔️: 			 arbitrary.  */
    // glibc❗✔️: 		      exp_limit = MAX_EXP + 3;
    // glibc❗✔️: 		    }
    // glibc❗✔️: 		  else
    // glibc❗✔️: 		    {
    // glibc❗✔️: 		      assert (lead_zero
    // glibc❗✔️: 			      <= (uintmax_t) (INTMAX_MAX - MAX_EXP - 3) / 4);
    // glibc❗✔️: 		      exp_limit = (MAX_EXP
    // glibc❗✔️: 				   + 4 * (intmax_t) lead_zero
    // glibc❗✔️: 				   + 3);
    // glibc❗✔️: 		    }
    // glibc❗✔️: 		}
    // glibc❗✔️: 	    }
    // glibc❗✔️: 	  else
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      if (exp_negative)
    // glibc❗✔️: 		{
    // glibc❗✔️: 		  assert (int_no
    // glibc❗✔️: 			  <= (uintmax_t) (INTMAX_MAX + MIN_10_EXP - MANT_DIG));
    // glibc❗✔️: 		  exp_limit = -MIN_10_EXP + MANT_DIG + (intmax_t) int_no;
    // glibc❗✔️: 		}
    // glibc❗✔️: 	      else
    // glibc❗✔️: 		{
    // glibc❗✔️: 		  if (int_no)
    // glibc❗✔️: 		    {
    // glibc❗✔️: 		      assert (lead_zero == 0
    // glibc❗✔️: 			      && int_no <= (uintmax_t) INTMAX_MAX);
    // glibc❗✔️: 		      exp_limit = MAX_10_EXP - (intmax_t) int_no + 1;
    // glibc❗✔️: 		    }
    // glibc❗✔️: 		  else if (lead_zero == (size_t) -1)
    // glibc❗✔️: 		    {
    // glibc❗✔️: 		      /* The number is zero and this limit is
    // glibc❗✔️: 			 arbitrary.  */
    // glibc❗✔️: 		      exp_limit = MAX_10_EXP + 1;
    // glibc❗✔️: 		    }
    // glibc❗✔️: 		  else
    // glibc❗✔️: 		    {
    // glibc❗✔️: 		      assert (lead_zero
    // glibc❗✔️: 			      <= (uintmax_t) (INTMAX_MAX - MAX_10_EXP - 1));
    // glibc❗✔️: 		      exp_limit = MAX_10_EXP + (intmax_t) lead_zero + 1;
    // glibc❗✔️: 		    }
    // glibc❗✔️: 		}
    // glibc❗✔️: 	    }
    // glibc❗✔️:
    // glibc❗✔️: 	  if (exp_limit < 0)
    // glibc❗✔️: 	    exp_limit = 0;
    // glibc❗✔️:
    // glibc❗✔️: 	  do
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      if (__builtin_expect ((exponent > exp_limit / 10
    // glibc❗✔️: 				     || (exponent == exp_limit / 10
    // glibc❗✔️: 					 && c - L_('0') > exp_limit % 10)), 0))
    // glibc❗✔️: 		/* The exponent is too large/small to represent a valid
    // glibc❗✔️: 		   number.  */
    // glibc❗✔️: 		{
    // glibc❗✔️: 		  FLOAT result;
    // glibc❗✔️:
    // glibc❗✔️: 		  /* We have to take care for special situation: a joker
    // glibc❗✔️: 		     might have written "0.0e100000" which is in fact
    // glibc❗✔️: 		     zero.  */
    // glibc❗✔️: 		  if (lead_zero == (size_t) -1)
    // glibc❗✔️: 		    result = negative ? -0.0 : 0.0;
    // glibc❗✔️: 		  else
    // glibc❗✔️: 		    {
    // glibc❗✔️: 		      /* Overflow or underflow.  */
    // glibc❗✔️: 		      result = (exp_negative
    // glibc❗✔️: 				? underflow_value (negative)
    // glibc❗✔️: 				: overflow_value (negative));
    // glibc❗✔️: 		    }
    // glibc❗✔️:
    // glibc❗✔️: 		  /* Accept all following digits as part of the exponent.  */
    // glibc❗✔️: 		  do
    // glibc❗✔️: 		    ++cp;
    // glibc❗✔️: 		  while (*cp >= L_('0') && *cp <= L_('9'));
    // glibc❗✔️:
    // glibc❗✔️: 		  RETURN (result, cp);
    // glibc❗✔️: 		  /* NOTREACHED */
    // glibc❗✔️: 		}
    // glibc❗✔️:
    // glibc❗✔️: 	      exponent *= 10;
    // glibc❗✔️: 	      exponent += c - L_('0');
    // glibc❗✔️:
    // glibc❗✔️: 	      c = *++cp;
    // glibc❗✔️: 	    }
    // glibc❗✔️: 	  while (c >= L_('0') && c <= L_('9'));
    // glibc❗✔️:
    // glibc❗✔️: 	  if (exp_negative)
    // glibc❗✔️: 	    exponent = -exponent;
    // glibc❗✔️: 	}
    // glibc❗✔️:       else
    // glibc❗✔️: 	cp = expp;
    // glibc❗✔️:     }
    // glibc❗✔️:
    // glibc❗✔️:   /* We don't want to have to work with trailing zeroes after the radix.  */
    // glibc❗✔️:   if (dig_no > int_no)
    // glibc❗✔️:     {
    // glibc❗✔️:       while (expp[-1] == L_('0'))
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  --expp;
    // glibc❗✔️: 	  --dig_no;
    // glibc❗✔️: 	}
    // glibc❗✔️:       assert (dig_no >= int_no);
    // glibc❗✔️:     }
    // glibc❗✔️:
    // glibc❗✔️:   if (dig_no == int_no && dig_no > 0 && exponent < 0)
    // glibc❗✔️:     do
    // glibc❗✔️:       {
    // glibc❗✔️: 	while (! (base == 16 ? ISXDIGIT (expp[-1]) : ISDIGIT (expp[-1])))
    // glibc❗✔️: 	  --expp;
    // glibc❗✔️:
    // glibc❗✔️: 	if (expp[-1] != L_('0'))
    // glibc❗✔️: 	  break;
    // glibc❗✔️:
    // glibc❗✔️: 	--expp;
    // glibc❗✔️: 	--dig_no;
    // glibc❗✔️: 	--int_no;
    // glibc❗✔️: 	exponent += base == 16 ? 4 : 1;
    // glibc❗✔️:       }
    // glibc❗✔️:     while (dig_no > 0 && exponent < 0);
    // glibc❗✔️:
    // glibc❗✔️:  number_parsed:
    // glibc❗✔️:
    // glibc❗✔️:   /* The whole string is parsed.  Store the address of the next character.  */
    // glibc❗✔️:   if (endptr)
    // glibc❗✔️:     *endptr = (STRING_TYPE *) cp;
    // glibc❗✔️:
    // glibc❗✔️:   if (dig_no == 0)
    // glibc❗✔️:     return negative ? -0.0 : 0.0;
    // glibc❗✔️:
    // glibc❗✔️:   if (lead_zero)
    // glibc❗✔️:     {
    // glibc❗✔️:       /* Find the decimal point */
    // glibc❗✔️: #ifdef USE_WIDE_CHAR
    // glibc❗✔️:       while (*startp != decimal)
    // glibc❗✔️: 	++startp;
    // glibc❗✔️: #else
    // glibc❗✔️:       while (1)
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  if (*startp == decimal[0])
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      for (cnt = 1; decimal[cnt] != '\0'; ++cnt)
    // glibc❗✔️: 		if (decimal[cnt] != startp[cnt])
    // glibc❗✔️: 		  break;
    // glibc❗✔️: 	      if (decimal[cnt] == '\0')
    // glibc❗✔️: 		break;
    // glibc❗✔️: 	    }
    // glibc❗✔️: 	  ++startp;
    // glibc❗✔️: 	}
    // glibc❗✔️: #endif
    // glibc❗✔️:       startp += lead_zero + decimal_len;
    // glibc❗✔️:       assert (lead_zero <= (base == 16
    // glibc❗✔️: 			    ? (uintmax_t) INTMAX_MAX / 4
    // glibc❗✔️: 			    : (uintmax_t) INTMAX_MAX));
    // glibc❗✔️:       assert (lead_zero <= (base == 16
    // glibc❗✔️: 			    ? ((uintmax_t) exponent
    // glibc❗✔️: 			       - (uintmax_t) INTMAX_MIN) / 4
    // glibc❗✔️: 			    : ((uintmax_t) exponent - (uintmax_t) INTMAX_MIN)));
    // glibc❗✔️:       exponent -= base == 16 ? 4 * (intmax_t) lead_zero : (intmax_t) lead_zero;
    // glibc❗✔️:       dig_no -= lead_zero;
    // glibc❗✔️:     }
    // glibc❗✔️:
    // glibc❗✔️:   /* If the BASE is 16 we can use a simpler algorithm.  */
    // glibc❗✔️:   if (base == 16)
    // glibc❗✔️:     {
    // glibc❗✔️:       static const int nbits[16] = { 0, 1, 2, 2, 3, 3, 3, 3,
    // glibc❗✔️: 				     4, 4, 4, 4, 4, 4, 4, 4 };
    // glibc❗✔️:       int idx = (MANT_DIG - 1) / BITS_PER_MP_LIMB;
    // glibc❗✔️:       int pos = (MANT_DIG - 1) % BITS_PER_MP_LIMB;
    // glibc❗✔️:       mp_limb_t val;
    // glibc❗✔️:
    // glibc❗✔️:       while (!ISXDIGIT (*startp))
    // glibc❗✔️: 	++startp;
    // glibc❗✔️:       while (*startp == L_('0'))
    // glibc❗✔️: 	++startp;
    // glibc❗✔️:       if (ISDIGIT (*startp))
    // glibc❗✔️: 	val = *startp++ - L_('0');
    // glibc❗✔️:       else
    // glibc❗✔️: 	val = 10 + TOLOWER (*startp++) - L_('a');
    // glibc❗✔️:       bits = nbits[val];
    // glibc❗✔️:       /* We cannot have a leading zero.  */
    // glibc❗✔️:       assert (bits != 0);
    // glibc❗✔️:
    // glibc❗✔️:       if (pos + 1 >= 4 || pos + 1 >= bits)
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  /* We don't have to care for wrapping.  This is the normal
    // glibc❗✔️: 	     case so we add the first clause in the `if' expression as
    // glibc❗✔️: 	     an optimization.  It is a compile-time constant and so does
    // glibc❗✔️: 	     not cost anything.  */
    // glibc❗✔️: 	  retval[idx] = val << (pos - bits + 1);
    // glibc❗✔️: 	  pos -= bits;
    // glibc❗✔️: 	}
    // glibc❗✔️:       else
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  retval[idx--] = val >> (bits - pos - 1);
    // glibc❗✔️: 	  retval[idx] = val << (BITS_PER_MP_LIMB - (bits - pos - 1));
    // glibc❗✔️: 	  pos = BITS_PER_MP_LIMB - 1 - (bits - pos - 1);
    // glibc❗✔️: 	}
    // glibc❗✔️:
    // glibc❗✔️:       /* Adjust the exponent for the bits we are shifting in.  */
    // glibc❗✔️:       assert (int_no <= (uintmax_t) (exponent < 0
    // glibc❗✔️: 				     ? (INTMAX_MAX - bits + 1) / 4
    // glibc❗✔️: 				     : (INTMAX_MAX - exponent - bits + 1) / 4));
    // glibc❗✔️:       exponent += bits - 1 + ((intmax_t) int_no - 1) * 4;
    // glibc❗✔️:
    // glibc❗✔️:       while (--dig_no > 0 && idx >= 0)
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  if (!ISXDIGIT (*startp))
    // glibc❗✔️: 	    startp += decimal_len;
    // glibc❗✔️: 	  if (ISDIGIT (*startp))
    // glibc❗✔️: 	    val = *startp++ - L_('0');
    // glibc❗✔️: 	  else
    // glibc❗✔️: 	    val = 10 + TOLOWER (*startp++) - L_('a');
    // glibc❗✔️:
    // glibc❗✔️: 	  if (pos + 1 >= 4)
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      retval[idx] |= val << (pos - 4 + 1);
    // glibc❗✔️: 	      pos -= 4;
    // glibc❗✔️: 	    }
    // glibc❗✔️: 	  else
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      retval[idx--] |= val >> (4 - pos - 1);
    // glibc❗✔️: 	      val <<= BITS_PER_MP_LIMB - (4 - pos - 1);
    // glibc❗✔️: 	      if (idx < 0)
    // glibc❗✔️: 		{
    // glibc❗✔️: 		  int rest_nonzero = 0;
    // glibc❗✔️: 		  while (--dig_no > 0)
    // glibc❗✔️: 		    {
    // glibc❗✔️: 		      if (*startp != L_('0'))
    // glibc❗✔️: 			{
    // glibc❗✔️: 			  rest_nonzero = 1;
    // glibc❗✔️: 			  break;
    // glibc❗✔️: 			}
    // glibc❗✔️: 		      startp++;
    // glibc❗✔️: 		    }
    // glibc❗✔️: 		  return round_and_return (retval, exponent, negative, val,
    // glibc❗✔️: 					   BITS_PER_MP_LIMB - 1, rest_nonzero);
    // glibc❗✔️: 		}
    // glibc❗✔️:
    // glibc❗✔️: 	      retval[idx] = val;
    // glibc❗✔️: 	      pos = BITS_PER_MP_LIMB - 1 - (4 - pos - 1);
    // glibc❗✔️: 	    }
    // glibc❗✔️: 	}
    // glibc❗✔️:
    // glibc❗✔️:       /* We ran out of digits.  */
    // glibc❗✔️:       MPN_ZERO (retval, idx);
    // glibc❗✔️:
    // glibc❗✔️:       return round_and_return (retval, exponent, negative, 0, 0, 0);
    // glibc❗✔️:     }
    // glibc❗✔️:
    // glibc❗✔️:   /* Now we have the number of digits in total and the integer digits as well
    // glibc❗✔️:      as the exponent and its sign.  We can decide whether the read digits are
    // glibc❗✔️:      really integer digits or belong to the fractional part; i.e. we normalize
    // glibc❗✔️:      123e-2 to 1.23.  */
    // glibc❗✔️:   {
    // glibc❗✔️:     intmax_t incr = (exponent < 0
    // glibc❗✔️: 		     ? MAX (-(intmax_t) int_no, exponent)
    // glibc❗✔️: 		     : MIN ((intmax_t) dig_no - (intmax_t) int_no, exponent));
    // glibc❗✔️:     int_no += incr;
    // glibc❗✔️:     exponent -= incr;
    // glibc❗✔️:   }
    // glibc❗✔️:
    // glibc❗✔️:   if (__glibc_unlikely (exponent > MAX_10_EXP + 1 - (intmax_t) int_no))
    // glibc❗✔️:     return overflow_value (negative);
    // glibc❗✔️:
    // glibc❗✔️:   /* 10^(MIN_10_EXP-1) is not normal.  Thus, 10^(MIN_10_EXP-1) /
    // glibc❗✔️:      2^MANT_DIG is below half the least subnormal, so anything with a
    // glibc❗✔️:      base-10 exponent less than the base-10 exponent (which is
    // glibc❗✔️:      MIN_10_EXP - 1 - ceil(MANT_DIG*log10(2))) of that value
    // glibc❗✔️:      underflows.  DIG is floor((MANT_DIG-1)log10(2)), so an exponent
    // glibc❗✔️:      below MIN_10_EXP - (DIG + 3) underflows.  But EXPONENT is
    // glibc❗✔️:      actually an exponent multiplied only by a fractional part, not an
    // glibc❗✔️:      integer part, so an exponent below MIN_10_EXP - (DIG + 2)
    // glibc❗✔️:      underflows.  */
    // glibc❗✔️:   if (__glibc_unlikely (exponent < MIN_10_EXP - (DIG + 2)))
    // glibc❗✔️:     return underflow_value (negative);
    // glibc❗✔️:
    // glibc❗✔️:   if (int_no > 0)
    // glibc❗✔️:     {
    // glibc❗✔️:       /* Read the integer part as a multi-precision number to NUM.  */
    // glibc❗✔️:       startp = str_to_mpn (startp, int_no, num, &numsize, &exponent
    // glibc❗✔️: #ifndef USE_WIDE_CHAR
    // glibc❗✔️: 			   , decimal, decimal_len, thousands
    // glibc❗✔️: #endif
    // glibc❗✔️: 			   );
    // glibc❗✔️:
    // glibc❗✔️:       if (exponent > 0)
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  /* We now multiply the gained number by the given power of ten.  */
    // glibc❗✔️: 	  mp_limb_t *psrc = num;
    // glibc❗✔️: 	  mp_limb_t *pdest = den;
    // glibc❗✔️: 	  int expbit = 1;
    // glibc❗✔️: 	  const struct mp_power *ttab = &_fpioconst_pow10[0];
    // glibc❗✔️:
    // glibc❗✔️: 	  do
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      if ((exponent & expbit) != 0)
    // glibc❗✔️: 		{
    // glibc❗✔️: 		  size_t size = ttab->arraysize - _FPIO_CONST_OFFSET;
    // glibc❗✔️: 		  mp_limb_t cy;
    // glibc❗✔️: 		  exponent ^= expbit;
    // glibc❗✔️:
    // glibc❗✔️: 		  /* FIXME: not the whole multiplication has to be
    // glibc❗✔️: 		     done.  If we have the needed number of bits we
    // glibc❗✔️: 		     only need the information whether more non-zero
    // glibc❗✔️: 		     bits follow.  */
    // glibc❗✔️: 		  if (numsize >= ttab->arraysize - _FPIO_CONST_OFFSET)
    // glibc❗✔️: 		    cy = __mpn_mul (pdest, psrc, numsize,
    // glibc❗✔️: 				    &__tens[ttab->arrayoff
    // glibc❗✔️: 					   + _FPIO_CONST_OFFSET],
    // glibc❗✔️: 				    size);
    // glibc❗✔️: 		  else
    // glibc❗✔️: 		    cy = __mpn_mul (pdest, &__tens[ttab->arrayoff
    // glibc❗✔️: 						  + _FPIO_CONST_OFFSET],
    // glibc❗✔️: 				    size, psrc, numsize);
    // glibc❗✔️: 		  numsize += size;
    // glibc❗✔️: 		  if (cy == 0)
    // glibc❗✔️: 		    --numsize;
    // glibc❗✔️: 		  (void) SWAP (psrc, pdest);
    // glibc❗✔️: 		}
    // glibc❗✔️: 	      expbit <<= 1;
    // glibc❗✔️: 	      ++ttab;
    // glibc❗✔️: 	    }
    // glibc❗✔️: 	  while (exponent != 0);
    // glibc❗✔️:
    // glibc❗✔️: 	  if (psrc == den)
    // glibc❗✔️: 	    memcpy (num, den, numsize * sizeof (mp_limb_t));
    // glibc❗✔️: 	}
    // glibc❗✔️:
    // glibc❗✔️:       /* Determine how many bits of the result we already have.  */
    // glibc❗✔️:       bits = stdc_leading_zeros (num[numsize - 1]);
    // glibc❗✔️:       bits = numsize * BITS_PER_MP_LIMB - bits;
    // glibc❗✔️:
    // glibc❗✔️:       /* Now we know the exponent of the number in base two.
    // glibc❗✔️: 	 Check it against the maximum possible exponent.  */
    // glibc❗✔️:       if (__glibc_unlikely (bits > MAX_EXP))
    // glibc❗✔️: 	return overflow_value (negative);
    // glibc❗✔️:
    // glibc❗✔️:       /* We have already the first BITS bits of the result.  Together with
    // glibc❗✔️: 	 the information whether more non-zero bits follow this is enough
    // glibc❗✔️: 	 to determine the result.  */
    // glibc❗✔️:       if (bits > MANT_DIG)
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  int i;
    // glibc❗✔️: 	  const mp_size_t least_idx = (bits - MANT_DIG) / BITS_PER_MP_LIMB;
    // glibc❗✔️: 	  const mp_size_t least_bit = (bits - MANT_DIG) % BITS_PER_MP_LIMB;
    // glibc❗✔️: 	  const mp_size_t round_idx = least_bit == 0 ? least_idx - 1
    // glibc❗✔️: 						     : least_idx;
    // glibc❗✔️: 	  const mp_size_t round_bit = least_bit == 0 ? BITS_PER_MP_LIMB - 1
    // glibc❗✔️: 						     : least_bit - 1;
    // glibc❗✔️:
    // glibc❗✔️: 	  if (least_bit == 0)
    // glibc❗✔️: 	    memcpy (retval, &num[least_idx],
    // glibc❗✔️: 		    RETURN_LIMB_SIZE * sizeof (mp_limb_t));
    // glibc❗✔️: 	  else
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      for (i = least_idx; i < numsize - 1; ++i)
    // glibc❗✔️: 		retval[i - least_idx] = (num[i] >> least_bit)
    // glibc❗✔️: 					| (num[i + 1]
    // glibc❗✔️: 					   << (BITS_PER_MP_LIMB - least_bit));
    // glibc❗✔️: 	      if (i - least_idx < RETURN_LIMB_SIZE)
    // glibc❗✔️: 		retval[RETURN_LIMB_SIZE - 1] = num[i] >> least_bit;
    // glibc❗✔️: 	    }
    // glibc❗✔️:
    // glibc❗✔️: 	  /* Check whether any limb beside the ones in RETVAL are non-zero.  */
    // glibc❗✔️: 	  for (i = 0; num[i] == 0; ++i)
    // glibc❗✔️: 	    ;
    // glibc❗✔️:
    // glibc❗✔️: 	  return round_and_return (retval, bits - 1, negative,
    // glibc❗✔️: 				   num[round_idx], round_bit,
    // glibc❗✔️: 				   int_no < dig_no || i < round_idx);
    // glibc❗✔️: 	  /* NOTREACHED */
    // glibc❗✔️: 	}
    // glibc❗✔️:       else if (dig_no == int_no)
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  const mp_size_t target_bit = (MANT_DIG - 1) % BITS_PER_MP_LIMB;
    // glibc❗✔️: 	  const mp_size_t is_bit = (bits - 1) % BITS_PER_MP_LIMB;
    // glibc❗✔️:
    // glibc❗✔️: 	  if (target_bit == is_bit)
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      memcpy (&retval[RETURN_LIMB_SIZE - numsize], num,
    // glibc❗✔️: 		      numsize * sizeof (mp_limb_t));
    // glibc❗✔️: 	      /* FIXME: the following loop can be avoided if we assume a
    // glibc❗✔️: 		 maximal MANT_DIG value.  */
    // glibc❗✔️: 	      MPN_ZERO (retval, RETURN_LIMB_SIZE - numsize);
    // glibc❗✔️: 	    }
    // glibc❗✔️: 	  else if (target_bit > is_bit)
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      (void) __mpn_lshift (&retval[RETURN_LIMB_SIZE - numsize],
    // glibc❗✔️: 				   num, numsize, target_bit - is_bit);
    // glibc❗✔️: 	      /* FIXME: the following loop can be avoided if we assume a
    // glibc❗✔️: 		 maximal MANT_DIG value.  */
    // glibc❗✔️: 	      MPN_ZERO (retval, RETURN_LIMB_SIZE - numsize);
    // glibc❗✔️: 	    }
    // glibc❗✔️: 	  else
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      mp_limb_t cy;
    // glibc❗✔️: 	      assert (numsize < RETURN_LIMB_SIZE);
    // glibc❗✔️:
    // glibc❗✔️: 	      cy = __mpn_rshift (&retval[RETURN_LIMB_SIZE - numsize],
    // glibc❗✔️: 				 num, numsize, is_bit - target_bit);
    // glibc❗✔️: 	      retval[RETURN_LIMB_SIZE - numsize - 1] = cy;
    // glibc❗✔️: 	      /* FIXME: the following loop can be avoided if we assume a
    // glibc❗✔️: 		 maximal MANT_DIG value.  */
    // glibc❗✔️: 	      MPN_ZERO (retval, RETURN_LIMB_SIZE - numsize - 1);
    // glibc❗✔️: 	    }
    // glibc❗✔️:
    // glibc❗✔️: 	  return round_and_return (retval, bits - 1, negative, 0, 0, 0);
    // glibc❗✔️: 	  /* NOTREACHED */
    // glibc❗✔️: 	}
    // glibc❗✔️:
    // glibc❗✔️:       /* Store the bits we already have.  */
    // glibc❗✔️:       memcpy (retval, num, numsize * sizeof (mp_limb_t));
    // glibc❗✔️: #if RETURN_LIMB_SIZE > 1
    // glibc❗✔️:       if (numsize < RETURN_LIMB_SIZE)
    // glibc❗✔️: # if RETURN_LIMB_SIZE == 2
    // glibc❗✔️: 	retval[numsize] = 0;
    // glibc❗✔️: # else
    // glibc❗✔️: 	MPN_ZERO (retval + numsize, RETURN_LIMB_SIZE - numsize);
    // glibc❗✔️: # endif
    // glibc❗✔️: #endif
    // glibc❗✔️:     }
    // glibc❗✔️:
    // glibc❗✔️:   /* We have to compute at least some of the fractional digits.  */
    // glibc❗✔️:   {
    // glibc❗✔️:     /* We construct a fraction and the result of the division gives us
    // glibc❗✔️:        the needed digits.  The denominator is 1.0 multiplied by the
    // glibc❗✔️:        exponent of the lowest digit; i.e. 0.123 gives 123 / 1000 and
    // glibc❗✔️:        123e-6 gives 123 / 1000000.  */
    // glibc❗✔️:
    // glibc❗✔️:     int expbit;
    // glibc❗✔️:     int neg_exp;
    // glibc❗✔️:     int more_bits;
    // glibc❗✔️:     int need_frac_digits;
    // glibc❗✔️:     mp_limb_t cy;
    // glibc❗✔️:     mp_limb_t *psrc = den;
    // glibc❗✔️:     mp_limb_t *pdest = num;
    // glibc❗✔️:     const struct mp_power *ttab = &_fpioconst_pow10[0];
    // glibc❗✔️:
    // glibc❗✔️:     assert (dig_no > int_no
    // glibc❗✔️: 	    && exponent <= 0
    // glibc❗✔️: 	    && exponent >= MIN_10_EXP - (DIG + 2));
    // glibc❗✔️:
    // glibc❗✔️:     /* We need to compute MANT_DIG - BITS fractional bits that lie
    // glibc❗✔️:        within the mantissa of the result, the following bit for
    // glibc❗✔️:        rounding, and to know whether any subsequent bit is 0.
    // glibc❗✔️:        Computing a bit with value 2^-n means looking at n digits after
    // glibc❗✔️:        the decimal point.  */
    // glibc❗✔️:     if (bits > 0)
    // glibc❗✔️:       {
    // glibc❗✔️: 	/* The bits required are those immediately after the point.  */
    // glibc❗✔️: 	assert (int_no > 0 && exponent == 0);
    // glibc❗✔️: 	need_frac_digits = 1 + MANT_DIG - bits;
    // glibc❗✔️:       }
    // glibc❗✔️:     else
    // glibc❗✔️:       {
    // glibc❗✔️: 	/* The number is in the form .123eEXPONENT.  */
    // glibc❗✔️: 	assert (int_no == 0 && *startp != L_('0'));
    // glibc❗✔️: 	/* The number is at least 10^(EXPONENT-1), and 10^3 <
    // glibc❗✔️: 	   2^10.  */
    // glibc❗✔️: 	int neg_exp_2 = ((1 - exponent) * 10) / 3 + 1;
    // glibc❗✔️: 	/* The number is at least 2^-NEG_EXP_2.  We need up to
    // glibc❗✔️: 	   MANT_DIG bits following that bit.  */
    // glibc❗✔️: 	need_frac_digits = neg_exp_2 + MANT_DIG;
    // glibc❗✔️: 	/* However, we never need bits beyond 1/4 ulp of the smallest
    // glibc❗✔️: 	   representable value.  (That 1/4 ulp bit is only needed to
    // glibc❗✔️: 	   determine tinyness on machines where tinyness is determined
    // glibc❗✔️: 	   after rounding.)  */
    // glibc❗✔️: 	if (need_frac_digits > MANT_DIG - MIN_EXP + 2)
    // glibc❗✔️: 	  need_frac_digits = MANT_DIG - MIN_EXP + 2;
    // glibc❗✔️: 	/* At this point, NEED_FRAC_DIGITS is the total number of
    // glibc❗✔️: 	   digits needed after the point, but some of those may be
    // glibc❗✔️: 	   leading 0s.  */
    // glibc❗✔️: 	need_frac_digits += exponent;
    // glibc❗✔️: 	/* Any cases underflowing enough that none of the fractional
    // glibc❗✔️: 	   digits are needed should have been caught earlier (such
    // glibc❗✔️: 	   cases are on the order of 10^-n or smaller where 2^-n is
    // glibc❗✔️: 	   the least subnormal).  */
    // glibc❗✔️: 	assert (need_frac_digits > 0);
    // glibc❗✔️:       }
    // glibc❗✔️:
    // glibc❗✔️:     if (need_frac_digits > (intmax_t) dig_no - (intmax_t) int_no)
    // glibc❗✔️:       need_frac_digits = (intmax_t) dig_no - (intmax_t) int_no;
    // glibc❗✔️:
    // glibc❗✔️:     if ((intmax_t) dig_no > (intmax_t) int_no + need_frac_digits)
    // glibc❗✔️:       {
    // glibc❗✔️: 	dig_no = int_no + need_frac_digits;
    // glibc❗✔️: 	more_bits = 1;
    // glibc❗✔️:       }
    // glibc❗✔️:     else
    // glibc❗✔️:       more_bits = 0;
    // glibc❗✔️:
    // glibc❗✔️:     neg_exp = (intmax_t) dig_no - (intmax_t) int_no - exponent;
    // glibc❗✔️:
    // glibc❗✔️:     /* Construct the denominator.  */
    // glibc❗✔️:     densize = 0;
    // glibc❗✔️:     expbit = 1;
    // glibc❗✔️:     do
    // glibc❗✔️:       {
    // glibc❗✔️: 	if ((neg_exp & expbit) != 0)
    // glibc❗✔️: 	  {
    // glibc❗✔️: 	    mp_limb_t cy;
    // glibc❗✔️: 	    neg_exp ^= expbit;
    // glibc❗✔️:
    // glibc❗✔️: 	    if (densize == 0)
    // glibc❗✔️: 	      {
    // glibc❗✔️: 		densize = ttab->arraysize - _FPIO_CONST_OFFSET;
    // glibc❗✔️: 		memcpy (psrc, &__tens[ttab->arrayoff + _FPIO_CONST_OFFSET],
    // glibc❗✔️: 			densize * sizeof (mp_limb_t));
    // glibc❗✔️: 	      }
    // glibc❗✔️: 	    else
    // glibc❗✔️: 	      {
    // glibc❗✔️: 		cy = __mpn_mul (pdest, &__tens[ttab->arrayoff
    // glibc❗✔️: 					      + _FPIO_CONST_OFFSET],
    // glibc❗✔️: 				ttab->arraysize - _FPIO_CONST_OFFSET,
    // glibc❗✔️: 				psrc, densize);
    // glibc❗✔️: 		densize += ttab->arraysize - _FPIO_CONST_OFFSET;
    // glibc❗✔️: 		if (cy == 0)
    // glibc❗✔️: 		  --densize;
    // glibc❗✔️: 		(void) SWAP (psrc, pdest);
    // glibc❗✔️: 	      }
    // glibc❗✔️: 	  }
    // glibc❗✔️: 	expbit <<= 1;
    // glibc❗✔️: 	++ttab;
    // glibc❗✔️:       }
    // glibc❗✔️:     while (neg_exp != 0);
    // glibc❗✔️:
    // glibc❗✔️:     if (psrc == num)
    // glibc❗✔️:       memcpy (den, num, densize * sizeof (mp_limb_t));
    // glibc❗✔️:
    // glibc❗✔️:     /* Read the fractional digits from the string.  */
    // glibc❗✔️:     (void) str_to_mpn (startp, dig_no - int_no, num, &numsize, &exponent
    // glibc❗✔️: #ifndef USE_WIDE_CHAR
    // glibc❗✔️: 		       , decimal, decimal_len, thousands
    // glibc❗✔️: #endif
    // glibc❗✔️: 		       );
    // glibc❗✔️:
    // glibc❗✔️:     /* We now have to shift both numbers so that the highest bit in the
    // glibc❗✔️:        denominator is set.  In the same process we copy the numerator to
    // glibc❗✔️:        a high place in the array so that the division constructs the wanted
    // glibc❗✔️:        digits.  This is done by a "quasi fix point" number representation.
    // glibc❗✔️:
    // glibc❗✔️:        num:   ddddddddddd . 0000000000000000000000
    // glibc❗✔️: 	      |--- m ---|
    // glibc❗✔️:        den:                            ddddddddddd      n >= m
    // glibc❗✔️: 				       |--- n ---|
    // glibc❗✔️:      */
    // glibc❗✔️:
    // glibc❗✔️:     cnt = stdc_leading_zeros (den[densize - 1]);
    // glibc❗✔️:
    // glibc❗✔️:
    // glibc❗✔️:     if (cnt > 0)
    // glibc❗✔️:       {
    // glibc❗✔️: 	/* Don't call `mpn_shift' with a count of zero since the specification
    // glibc❗✔️: 	   does not allow this.  */
    // glibc❗✔️: 	(void) __mpn_lshift (den, den, densize, cnt);
    // glibc❗✔️: 	cy = __mpn_lshift (num, num, numsize, cnt);
    // glibc❗✔️: 	if (cy != 0)
    // glibc❗✔️: 	  num[numsize++] = cy;
    // glibc❗✔️:       }
    // glibc❗✔️:
    // glibc❗✔️:     /* Now we are ready for the division.  But it is not necessary to
    // glibc❗✔️:        do a full multi-precision division because we only need a small
    // glibc❗✔️:        number of bits for the result.  So we do not use __mpn_divmod
    // glibc❗✔️:        here but instead do the division here by hand and stop whenever
    // glibc❗✔️:        the needed number of bits is reached.  The code itself comes
    // glibc❗✔️:        from the GNU MP Library by Torbj\"orn Granlund.  */
    // glibc❗✔️:
    // glibc❗✔️:     exponent = bits;
    // glibc❗✔️:
    // glibc❗✔️:     switch (densize)
    // glibc❗✔️:       {
    // glibc❗✔️:       case 1:
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  mp_limb_t d, n, quot;
    // glibc❗✔️: 	  int used = 0;
    // glibc❗✔️:
    // glibc❗✔️: 	  n = num[0];
    // glibc❗✔️: 	  d = den[0];
    // glibc❗✔️: 	  assert (numsize == 1 && n < d);
    // glibc❗✔️:
    // glibc❗✔️: 	  do
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      udiv_qrnnd (quot, n, n, 0, d);
    // glibc❗✔️:
    // glibc❗✔️: #define got_limb							      \
    // glibc❗✔️: 	      if (bits == 0)						      \
    // glibc❗✔️: 		{							      \
    // glibc❗✔️: 		  int cnt = stdc_leading_zeros (quot);			      \
    // glibc❗✔️: 		  exponent -= cnt;					      \
    // glibc❗✔️: 		  if (BITS_PER_MP_LIMB - cnt > MANT_DIG)		      \
    // glibc❗✔️: 		    {							      \
    // glibc❗✔️: 		      used = MANT_DIG + cnt;				      \
    // glibc❗✔️: 		      retval[0] = quot >> (BITS_PER_MP_LIMB - used);	      \
    // glibc❗✔️: 		      bits = MANT_DIG + 1;				      \
    // glibc❗✔️: 		    }							      \
    // glibc❗✔️: 		  else							      \
    // glibc❗✔️: 		    {							      \
    // glibc❗✔️: 		      /* Note that we only clear the second element.  */      \
    // glibc❗✔️: 		      /* The conditional is determined at compile time.  */   \
    // glibc❗✔️: 		      if (RETURN_LIMB_SIZE > 1)				      \
    // glibc❗✔️: 			retval[1] = 0;					      \
    // glibc❗✔️: 		      retval[0] = quot;					      \
    // glibc❗✔️: 		      bits = -cnt;					      \
    // glibc❗✔️: 		    }							      \
    // glibc❗✔️: 		}							      \
    // glibc❗✔️: 	      else if (bits + BITS_PER_MP_LIMB <= MANT_DIG)		      \
    // glibc❗✔️: 		__mpn_lshift_1 (retval, RETURN_LIMB_SIZE, BITS_PER_MP_LIMB,   \
    // glibc❗✔️: 				quot);					      \
    // glibc❗✔️: 	      else							      \
    // glibc❗✔️: 		{							      \
    // glibc❗✔️: 		  used = MANT_DIG - bits;				      \
    // glibc❗✔️: 		  if (used > 0)						      \
    // glibc❗✔️: 		    __mpn_lshift_1 (retval, RETURN_LIMB_SIZE, used, quot);    \
    // glibc❗✔️: 		}							      \
    // glibc❗✔️: 	      bits += BITS_PER_MP_LIMB
    // glibc❗✔️:
    // glibc❗✔️: 	      got_limb;
    // glibc❗✔️: 	    }
    // glibc❗✔️: 	  while (bits <= MANT_DIG);
    // glibc❗✔️:
    // glibc❗✔️: 	  return round_and_return (retval, exponent - 1, negative,
    // glibc❗✔️: 				   quot, BITS_PER_MP_LIMB - 1 - used,
    // glibc❗✔️: 				   more_bits || n != 0);
    // glibc❗✔️: 	}
    // glibc❗✔️:       case 2:
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  mp_limb_t d0, d1, n0, n1;
    // glibc❗✔️: 	  mp_limb_t quot = 0;
    // glibc❗✔️: 	  int used = 0;
    // glibc❗✔️:
    // glibc❗✔️: 	  d0 = den[0];
    // glibc❗✔️: 	  d1 = den[1];
    // glibc❗✔️:
    // glibc❗✔️: 	  if (numsize < densize)
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      if (num[0] >= d1)
    // glibc❗✔️: 		{
    // glibc❗✔️: 		  /* The numerator of the number occupies fewer bits than
    // glibc❗✔️: 		     the denominator but the one limb is bigger than the
    // glibc❗✔️: 		     high limb of the numerator.  */
    // glibc❗✔️: 		  n1 = 0;
    // glibc❗✔️: 		  n0 = num[0];
    // glibc❗✔️: 		}
    // glibc❗✔️: 	      else
    // glibc❗✔️: 		{
    // glibc❗✔️: 		  if (bits <= 0)
    // glibc❗✔️: 		    exponent -= BITS_PER_MP_LIMB;
    // glibc❗✔️: 		  else
    // glibc❗✔️: 		    {
    // glibc❗✔️: 		      if (bits + BITS_PER_MP_LIMB <= MANT_DIG)
    // glibc❗✔️: 			__mpn_lshift_1 (retval, RETURN_LIMB_SIZE,
    // glibc❗✔️: 					BITS_PER_MP_LIMB, 0);
    // glibc❗✔️: 		      else
    // glibc❗✔️: 			{
    // glibc❗✔️: 			  used = MANT_DIG - bits;
    // glibc❗✔️: 			  if (used > 0)
    // glibc❗✔️: 			    __mpn_lshift_1 (retval, RETURN_LIMB_SIZE, used, 0);
    // glibc❗✔️: 			}
    // glibc❗✔️: 		      bits += BITS_PER_MP_LIMB;
    // glibc❗✔️: 		    }
    // glibc❗✔️: 		  n1 = num[0];
    // glibc❗✔️: 		  n0 = 0;
    // glibc❗✔️: 		}
    // glibc❗✔️: 	    }
    // glibc❗✔️: 	  else
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      n1 = num[1];
    // glibc❗✔️: 	      n0 = num[0];
    // glibc❗✔️: 	    }
    // glibc❗✔️:
    // glibc❗✔️: 	  while (bits <= MANT_DIG)
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      mp_limb_t r;
    // glibc❗✔️:
    // glibc❗✔️: 	      if (n1 == d1)
    // glibc❗✔️: 		{
    // glibc❗✔️: 		  /* QUOT should be either 111..111 or 111..110.  We need
    // glibc❗✔️: 		     special treatment of this rare case as normal division
    // glibc❗✔️: 		     would give overflow.  */
    // glibc❗✔️: 		  quot = ~(mp_limb_t) 0;
    // glibc❗✔️:
    // glibc❗✔️: 		  r = n0 + d1;
    // glibc❗✔️: 		  if (r < d1)	/* Carry in the addition?  */
    // glibc❗✔️: 		    {
    // glibc❗✔️: 		      add_ssaaaa (n1, n0, r - d0, 0, 0, d0);
    // glibc❗✔️: 		      goto have_quot;
    // glibc❗✔️: 		    }
    // glibc❗✔️: 		  n1 = d0 - (d0 != 0);
    // glibc❗✔️: 		  n0 = -d0;
    // glibc❗✔️: 		}
    // glibc❗✔️: 	      else
    // glibc❗✔️: 		{
    // glibc❗✔️: 		  udiv_qrnnd (quot, r, n1, n0, d1);
    // glibc❗✔️: 		  umul_ppmm (n1, n0, d0, quot);
    // glibc❗✔️: 		}
    // glibc❗✔️:
    // glibc❗✔️: 	    q_test:
    // glibc❗✔️: 	      if (n1 > r || (n1 == r && n0 > 0))
    // glibc❗✔️: 		{
    // glibc❗✔️: 		  /* The estimated QUOT was too large.  */
    // glibc❗✔️: 		  --quot;
    // glibc❗✔️:
    // glibc❗✔️: 		  sub_ddmmss (n1, n0, n1, n0, 0, d0);
    // glibc❗✔️: 		  r += d1;
    // glibc❗✔️: 		  if (r >= d1)	/* If not carry, test QUOT again.  */
    // glibc❗✔️: 		    goto q_test;
    // glibc❗✔️: 		}
    // glibc❗✔️: 	      sub_ddmmss (n1, n0, r, 0, n1, n0);
    // glibc❗✔️:
    // glibc❗✔️: 	    have_quot:
    // glibc❗✔️: 	      got_limb;
    // glibc❗✔️: 	    }
    // glibc❗✔️:
    // glibc❗✔️: 	  return round_and_return (retval, exponent - 1, negative,
    // glibc❗✔️: 				   quot, BITS_PER_MP_LIMB - 1 - used,
    // glibc❗✔️: 				   more_bits || n1 != 0 || n0 != 0);
    // glibc❗✔️: 	}
    // glibc❗✔️:       default:
    // glibc❗✔️: 	{
    // glibc❗✔️: 	  int i;
    // glibc❗✔️: 	  mp_limb_t cy, dX, d1, n0, n1;
    // glibc❗✔️: 	  mp_limb_t quot = 0;
    // glibc❗✔️: 	  int used = 0;
    // glibc❗✔️:
    // glibc❗✔️: 	  dX = den[densize - 1];
    // glibc❗✔️: 	  d1 = den[densize - 2];
    // glibc❗✔️:
    // glibc❗✔️: 	  /* The division does not work if the upper limb of the two-limb
    // glibc❗✔️: 	     numerator is greater than or equal to the denominator.  */
    // glibc❗✔️: 	  if (__mpn_cmp (num, &den[densize - numsize], numsize) >= 0)
    // glibc❗✔️: 	    num[numsize++] = 0;
    // glibc❗✔️:
    // glibc❗✔️: 	  if (numsize < densize)
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      mp_size_t empty = densize - numsize;
    // glibc❗✔️: 	      int i;
    // glibc❗✔️:
    // glibc❗✔️: 	      if (bits <= 0)
    // glibc❗✔️: 		exponent -= empty * BITS_PER_MP_LIMB;
    // glibc❗✔️: 	      else
    // glibc❗✔️: 		{
    // glibc❗✔️: 		  if (bits + empty * BITS_PER_MP_LIMB <= MANT_DIG)
    // glibc❗✔️: 		    {
    // glibc❗✔️: 		      /* We make a difference here because the compiler
    // glibc❗✔️: 			 cannot optimize the `else' case that good and
    // glibc❗✔️: 			 this reflects all currently used FLOAT types
    // glibc❗✔️: 			 and GMP implementations.  */
    // glibc❗✔️: #if RETURN_LIMB_SIZE <= 2
    // glibc❗✔️: 		      assert (empty == 1);
    // glibc❗✔️: 		      __mpn_lshift_1 (retval, RETURN_LIMB_SIZE,
    // glibc❗✔️: 				      BITS_PER_MP_LIMB, 0);
    // glibc❗✔️: #else
    // glibc❗✔️: 		      for (i = RETURN_LIMB_SIZE - 1; i >= empty; --i)
    // glibc❗✔️: 			retval[i] = retval[i - empty];
    // glibc❗✔️: 		      while (i >= 0)
    // glibc❗✔️: 			retval[i--] = 0;
    // glibc❗✔️: #endif
    // glibc❗✔️: 		    }
    // glibc❗✔️: 		  else
    // glibc❗✔️: 		    {
    // glibc❗✔️: 		      used = MANT_DIG - bits;
    // glibc❗✔️: 		      if (used >= BITS_PER_MP_LIMB)
    // glibc❗✔️: 			{
    // glibc❗✔️: 			  int i;
    // glibc❗✔️: 			  (void) __mpn_lshift (&retval[used
    // glibc❗✔️: 						       / BITS_PER_MP_LIMB],
    // glibc❗✔️: 					       retval,
    // glibc❗✔️: 					       (RETURN_LIMB_SIZE
    // glibc❗✔️: 						- used / BITS_PER_MP_LIMB),
    // glibc❗✔️: 					       used % BITS_PER_MP_LIMB);
    // glibc❗✔️: 			  for (i = used / BITS_PER_MP_LIMB - 1; i >= 0; --i)
    // glibc❗✔️: 			    retval[i] = 0;
    // glibc❗✔️: 			}
    // glibc❗✔️: 		      else if (used > 0)
    // glibc❗✔️: 			__mpn_lshift_1 (retval, RETURN_LIMB_SIZE, used, 0);
    // glibc❗✔️: 		    }
    // glibc❗✔️: 		  bits += empty * BITS_PER_MP_LIMB;
    // glibc❗✔️: 		}
    // glibc❗✔️: 	      for (i = numsize; i > 0; --i)
    // glibc❗✔️: 		num[i + empty] = num[i - 1];
    // glibc❗✔️: 	      MPN_ZERO (num, empty + 1);
    // glibc❗✔️: 	    }
    // glibc❗✔️: 	  else
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      int i;
    // glibc❗✔️: 	      assert (numsize == densize);
    // glibc❗✔️: 	      for (i = numsize; i > 0; --i)
    // glibc❗✔️: 		num[i] = num[i - 1];
    // glibc❗✔️: 	      num[0] = 0;
    // glibc❗✔️: 	    }
    // glibc❗✔️:
    // glibc❗✔️: 	  den[densize] = 0;
    // glibc❗✔️: 	  n0 = num[densize];
    // glibc❗✔️:
    // glibc❗✔️: 	  while (bits <= MANT_DIG)
    // glibc❗✔️: 	    {
    // glibc❗✔️: 	      if (n0 == dX)
    // glibc❗✔️: 		/* This might over-estimate QUOT, but it's probably not
    // glibc❗✔️: 		   worth the extra code here to find out.  */
    // glibc❗✔️: 		quot = ~(mp_limb_t) 0;
    // glibc❗✔️: 	      else
    // glibc❗✔️: 		{
    // glibc❗✔️: 		  mp_limb_t r;
    // glibc❗✔️:
    // glibc❗✔️: 		  udiv_qrnnd (quot, r, n0, num[densize - 1], dX);
    // glibc❗✔️: 		  umul_ppmm (n1, n0, d1, quot);
    // glibc❗✔️:
    // glibc❗✔️: 		  while (n1 > r || (n1 == r && n0 > num[densize - 2]))
    // glibc❗✔️: 		    {
    // glibc❗✔️: 		      --quot;
    // glibc❗✔️: 		      r += dX;
    // glibc❗✔️: 		      if (r < dX) /* I.e. "carry in previous addition?" */
    // glibc❗✔️: 			break;
    // glibc❗✔️: 		      n1 -= n0 < d1;
    // glibc❗✔️: 		      n0 -= d1;
    // glibc❗✔️: 		    }
    // glibc❗✔️: 		}
    // glibc❗✔️:
    // glibc❗✔️: 	      /* Possible optimization: We already have (q * n0) and (1 * n1)
    // glibc❗✔️: 		 after the calculation of QUOT.  Taking advantage of this, we
    // glibc❗✔️: 		 could make this loop make two iterations less.  */
    // glibc❗✔️:
    // glibc❗✔️: 	      cy = __mpn_submul_1 (num, den, densize + 1, quot);
    // glibc❗✔️:
    // glibc❗✔️: 	      if (num[densize] != cy)
    // glibc❗✔️: 		{
    // glibc❗✔️: 		  cy = __mpn_add_n (num, num, den, densize);
    // glibc❗✔️: 		  assert (cy != 0);
    // glibc❗✔️: 		  --quot;
    // glibc❗✔️: 		}
    // glibc❗✔️: 	      n0 = num[densize] = num[densize - 1];
    // glibc❗✔️: 	      for (i = densize - 1; i > 0; --i)
    // glibc❗✔️: 		num[i] = num[i - 1];
    // glibc❗✔️: 	      num[0] = 0;
    // glibc❗✔️:
    // glibc❗✔️: 	      got_limb;
    // glibc❗✔️: 	    }
    // glibc❗✔️:
    // glibc❗✔️: 	  for (i = densize; i >= 0 && num[i] == 0; --i)
    // glibc❗✔️: 	    ;
    // glibc❗✔️: 	  return round_and_return (retval, exponent - 1, negative,
    // glibc❗✔️: 				   quot, BITS_PER_MP_LIMB - 1 - used,
    // glibc❗✔️: 				   more_bits || i >= 0);
    // glibc❗✔️: 	}
    // glibc❗✔️:       }
    // glibc❗✔️:   }
    // glibc❗✔️:
    // glibc❗✔️:   /* NOTREACHED */
    // glibc❗✔️: }

    debug_assert_eq!((bits, emin, pok), (53, -1074, true));
    let bytes = cursor.remaining_bytes();
    let byte = |i: usize| bytes.get(i).copied().unwrap_or(0);
    let mut cp = 0usize;
    while byte(cp) == b'0' {
        cp += 1;
    }
    let start_of_digits = 0usize;
    let mut startp = cp;
    let c = byte(cp);
    if !(digit(c).is_some()
        || (c == b'.' && (cp != start_of_digits || digit(byte(cp + 1)).is_some()))
        || (cp != start_of_digits && c.eq_ignore_ascii_case(&b'p')))
    {
        if cp == start_of_digits {
            // No legal hex digit: source returns at the x of the 0x prefix.
            cursor.ungetc();
        } else {
            cursor.advance_bytes(cp);
        }
        return if sign < 0 { -0.0 } else { 0.0 };
    }
    let mut dig_no = 0usize;
    while digit(byte(cp)).is_some() {
        dig_no += 1;
        cp += 1;
    }
    let mut int_no = dig_no;
    let mut lead_zero = if int_no == 0 { usize::MAX } else { 0 };
    if byte(cp) == b'.' {
        cp += 1;
        while digit(byte(cp)).is_some() {
            if byte(cp) != b'0' && lead_zero == usize::MAX {
                lead_zero = dig_no - int_no;
            }
            dig_no += 1;
            cp += 1;
        }
    }
    let mut expp = cp;
    let mut exponent = 0i64;
    let mut early = None;
    if byte(cp).eq_ignore_ascii_case(&b'p') {
        cp += 1;
        let exp_negative = byte(cp) == b'-';
        if matches!(byte(cp), b'-' | b'+') {
            cp += 1;
        }
        if byte(cp).is_ascii_digit() {
            let exp_limit = if exp_negative {
                1074 + 4 * int_no as i64
            } else if int_no != 0 {
                1024 - 4 * int_no as i64 + 3
            } else if lead_zero == usize::MAX {
                1024 + 3
            } else {
                1024 + 4 * lead_zero as i64 + 3
            }
            .max(0);
            while byte(cp).is_ascii_digit() {
                let d = i64::from(byte(cp) - b'0');
                if exponent > exp_limit / 10 || (exponent == exp_limit / 10 && d > exp_limit % 10) {
                    early = Some(if lead_zero == usize::MAX {
                        0.0
                    } else if exp_negative {
                        f64::MIN_POSITIVE * f64::MIN_POSITIVE
                    } else {
                        f64::MAX * f64::MAX
                    });
                    while byte(cp).is_ascii_digit() {
                        cp += 1;
                    }
                    break;
                }
                exponent = exponent * 10 + d;
                cp += 1;
            }
            if exp_negative {
                exponent = -exponent;
            }
        } else {
            cp = expp;
        }
    }
    if let Some(value) = early {
        cursor.advance_bytes(cp);
        return if sign < 0 { -value } else { value };
    }
    if dig_no > int_no {
        while byte(expp - 1) == b'0' {
            expp -= 1;
            dig_no -= 1;
        }
    }
    if dig_no == int_no && dig_no > 0 && exponent < 0 {
        while dig_no > 0 && exponent < 0 {
            while digit(byte(expp - 1)).is_none() {
                expp -= 1;
            }
            if byte(expp - 1) != b'0' {
                break;
            }
            expp -= 1;
            dig_no -= 1;
            int_no -= 1;
            exponent += 4;
        }
    }
    cursor.advance_bytes(cp);
    if dig_no == 0 {
        return if sign < 0 { -0.0 } else { 0.0 };
    }
    if lead_zero != 0 {
        while byte(startp) != b'.' {
            startp += 1;
        }
        startp += lead_zero + 1;
        exponent -= 4 * lead_zero as i64;
        dig_no -= lead_zero;
    }
    while digit(byte(startp)).is_none() {
        startp += 1;
    }
    while byte(startp) == b'0' {
        startp += 1;
    }
    let mut val = digit(byte(startp)).expect("source first significant hex digit");
    startp += 1;
    let first_bits = 64 - val.leading_zeros() as i32;
    let mut retval = val << (53 - first_bits);
    let mut pos = 52 - first_bits;
    exponent += i64::from(first_bits - 1) + (int_no as i64 - 1) * 4;
    dig_no -= 1;
    while dig_no > 0 {
        if digit(byte(startp)).is_none() {
            startp += 1;
        }
        val = digit(byte(startp)).expect("source normalized hex digit");
        startp += 1;
        if pos + 1 >= 4 {
            retval |= val << (pos - 3);
            pos -= 4;
        } else {
            retval |= val >> (3 - pos);
            val = val.wrapping_shl((64 - (3 - pos)) as u32);
            let mut rest_nonzero = false;
            dig_no -= 1;
            while dig_no > 0 {
                if byte(startp) != b'0' {
                    rest_nonzero = true;
                    break;
                }
                startp += 1;
                dig_no -= 1;
            }
            return round_and_return(retval, exponent, sign < 0, val, 63, rest_nonzero);
        }
        dig_no -= 1;
    }
    round_and_return(retval, exponent, sign < 0, 0, 0, false)
}

fn round_and_return(
    mut retval: u64,
    mut exponent: i64,
    negative: bool,
    mut round_limb: u64,
    mut round_bit: u32,
    mut more_bits: bool,
) -> f64 {
    // glibc✔️✔️: static FLOAT
    // glibc✔️✔️: round_and_return (mp_limb_t *retval, intmax_t exponent, int negative,
    // glibc✔️✔️: 		  mp_limb_t round_limb, mp_size_t round_bit, int more_bits)
    // glibc✔️✔️: {
    // glibc✔️✔️:   int mode = get_rounding_mode ();
    // glibc✔️✔️:
    // glibc✔️✔️:   if (exponent < MIN_EXP - 1)
    // glibc✔️✔️:     {
    // glibc✔️✔️:       if (exponent < MIN_EXP - 1 - MANT_DIG)
    // glibc✔️✔️: 	return underflow_value (negative);
    // glibc✔️✔️:
    // glibc✔️✔️:       mp_size_t shift = MIN_EXP - 1 - exponent;
    // glibc✔️✔️:       bool is_tiny = true;
    // glibc✔️✔️:       bool old_half_bit = (round_limb & (((mp_limb_t) 1) << round_bit)) != 0;
    // glibc✔️✔️:
    // glibc✔️✔️:       more_bits |= (round_limb & ((((mp_limb_t) 1) << round_bit) - 1)) != 0;
    // glibc✔️✔️:       if (shift == MANT_DIG)
    // glibc✔️✔️: 	/* This is a special case to handle the very seldom case where
    // glibc✔️✔️: 	   the mantissa will be empty after the shift.  */
    // glibc✔️✔️: 	{
    // glibc✔️✔️: 	  int i;
    // glibc✔️✔️:
    // glibc✔️✔️: 	  round_limb = retval[RETURN_LIMB_SIZE - 1];
    // glibc✔️✔️: 	  round_bit = (MANT_DIG - 1) % BITS_PER_MP_LIMB;
    // glibc✔️✔️: 	  for (i = 0; i < RETURN_LIMB_SIZE - 1; ++i)
    // glibc✔️✔️: 	    more_bits |= retval[i] != 0;
    // glibc✔️✔️: 	  MPN_ZERO (retval, RETURN_LIMB_SIZE);
    // glibc✔️✔️: 	}
    // glibc✔️✔️:       else if (shift >= BITS_PER_MP_LIMB)
    // glibc✔️✔️: 	{
    // glibc✔️✔️: 	  int i;
    // glibc✔️✔️:
    // glibc✔️✔️: 	  round_limb = retval[(shift - 1) / BITS_PER_MP_LIMB];
    // glibc✔️✔️: 	  round_bit = (shift - 1) % BITS_PER_MP_LIMB;
    // glibc✔️✔️: 	  for (i = 0; i < (shift - 1) / BITS_PER_MP_LIMB; ++i)
    // glibc✔️✔️: 	    more_bits |= retval[i] != 0;
    // glibc✔️✔️: 	  more_bits |= ((round_limb & ((((mp_limb_t) 1) << round_bit) - 1))
    // glibc✔️✔️: 			!= 0);
    // glibc✔️✔️:
    // glibc✔️✔️: 	  /* __mpn_rshift requires 0 < shift < BITS_PER_MP_LIMB.  */
    // glibc✔️✔️: 	  if ((shift % BITS_PER_MP_LIMB) != 0)
    // glibc✔️✔️: 	    (void) __mpn_rshift (retval, &retval[shift / BITS_PER_MP_LIMB],
    // glibc✔️✔️: 			         RETURN_LIMB_SIZE - (shift / BITS_PER_MP_LIMB),
    // glibc✔️✔️: 			         shift % BITS_PER_MP_LIMB);
    // glibc✔️✔️: 	  else
    // glibc✔️✔️: 	    for (i = 0; i < RETURN_LIMB_SIZE - (shift / BITS_PER_MP_LIMB); i++)
    // glibc✔️✔️: 	      retval[i] = retval[i + (shift / BITS_PER_MP_LIMB)];
    // glibc✔️✔️: 	  MPN_ZERO (&retval[RETURN_LIMB_SIZE - (shift / BITS_PER_MP_LIMB)],
    // glibc✔️✔️: 		    shift / BITS_PER_MP_LIMB);
    // glibc✔️✔️: 	}
    // glibc✔️✔️:       else if (shift > 0)
    // glibc✔️✔️: 	{
    // glibc✔️✔️: 	  if (TININESS_AFTER_ROUNDING && shift == 1)
    // glibc✔️✔️: 	    {
    // glibc✔️✔️: 	      /* Whether the result counts as tiny depends on whether,
    // glibc✔️✔️: 		 after rounding to the normal precision, it still has
    // glibc✔️✔️: 		 a subnormal exponent.  */
    // glibc✔️✔️: 	      mp_limb_t retval_normal[RETURN_LIMB_SIZE];
    // glibc✔️✔️: 	      if (round_away (negative,
    // glibc✔️✔️: 			      (retval[0] & 1) != 0,
    // glibc✔️✔️: 			      (round_limb
    // glibc✔️✔️: 			       & (((mp_limb_t) 1) << round_bit)) != 0,
    // glibc✔️✔️: 			      (more_bits
    // glibc✔️✔️: 			       || ((round_limb
    // glibc✔️✔️: 				    & ((((mp_limb_t) 1) << round_bit) - 1))
    // glibc✔️✔️: 				   != 0)),
    // glibc✔️✔️: 			      mode))
    // glibc✔️✔️: 		{
    // glibc✔️✔️: 		  mp_limb_t cy = __mpn_add_1 (retval_normal, retval,
    // glibc✔️✔️: 					      RETURN_LIMB_SIZE, 1);
    // glibc✔️✔️:
    // glibc✔️✔️: 		  if (((MANT_DIG % BITS_PER_MP_LIMB) == 0 && cy)
    // glibc✔️✔️: 		      || ((MANT_DIG % BITS_PER_MP_LIMB) != 0
    // glibc✔️✔️: 			  && ((retval_normal[RETURN_LIMB_SIZE - 1]
    // glibc✔️✔️: 			       & (((mp_limb_t) 1)
    // glibc✔️✔️: 				  << (MANT_DIG % BITS_PER_MP_LIMB)))
    // glibc✔️✔️: 			      != 0)))
    // glibc✔️✔️: 		    is_tiny = false;
    // glibc✔️✔️: 		}
    // glibc✔️✔️: 	    }
    // glibc✔️✔️: 	  round_limb = retval[0];
    // glibc✔️✔️: 	  round_bit = shift - 1;
    // glibc✔️✔️: 	  (void) __mpn_rshift (retval, retval, RETURN_LIMB_SIZE, shift);
    // glibc✔️✔️: 	}
    // glibc✔️✔️:       more_bits |= old_half_bit;
    // glibc✔️✔️:       /* This is a hook for the m68k long double format, where the
    // glibc✔️✔️: 	 exponent bias is the same for normalized and denormalized
    // glibc✔️✔️: 	 numbers.  */
    // glibc✔️✔️: #ifndef DENORM_EXP
    // glibc✔️✔️: # define DENORM_EXP (MIN_EXP - 2)
    // glibc✔️✔️: #endif
    // glibc✔️✔️:       exponent = DENORM_EXP;
    // glibc✔️✔️:       if (is_tiny
    // glibc✔️✔️: 	  && ((round_limb & (((mp_limb_t) 1) << round_bit)) != 0
    // glibc✔️✔️: 	      || more_bits
    // glibc✔️✔️: 	      || (round_limb & ((((mp_limb_t) 1) << round_bit) - 1)) != 0))
    // glibc✔️✔️: 	{
    // glibc✔️✔️: 	  __set_errno (ERANGE);
    // glibc✔️✔️: 	  FLOAT force_underflow = MIN_VALUE * MIN_VALUE;
    // glibc✔️✔️: 	  math_force_eval (force_underflow);
    // glibc✔️✔️: 	}
    // glibc✔️✔️:     }
    // glibc✔️✔️:
    // glibc✔️✔️:   if (exponent >= MAX_EXP)
    // glibc✔️✔️:     goto overflow;
    // glibc✔️✔️:
    // glibc✔️✔️:   bool half_bit = (round_limb & (((mp_limb_t) 1) << round_bit)) != 0;
    // glibc✔️✔️:   bool more_bits_nonzero
    // glibc✔️✔️:     = (more_bits
    // glibc✔️✔️:        || (round_limb & ((((mp_limb_t) 1) << round_bit) - 1)) != 0);
    // glibc✔️✔️:   if (round_away (negative,
    // glibc✔️✔️: 		  (retval[0] & 1) != 0,
    // glibc✔️✔️: 		  half_bit,
    // glibc✔️✔️: 		  more_bits_nonzero,
    // glibc✔️✔️: 		  mode))
    // glibc✔️✔️:     {
    // glibc✔️✔️:       mp_limb_t cy = __mpn_add_1 (retval, retval, RETURN_LIMB_SIZE, 1);
    // glibc✔️✔️:
    // glibc✔️✔️:       if (((MANT_DIG % BITS_PER_MP_LIMB) == 0 && cy)
    // glibc✔️✔️: 	  || ((MANT_DIG % BITS_PER_MP_LIMB) != 0
    // glibc✔️✔️: 	      && (retval[RETURN_LIMB_SIZE - 1]
    // glibc✔️✔️: 		  & (((mp_limb_t) 1) << (MANT_DIG % BITS_PER_MP_LIMB))) != 0))
    // glibc✔️✔️: 	{
    // glibc✔️✔️: 	  ++exponent;
    // glibc✔️✔️: 	  (void) __mpn_rshift (retval, retval, RETURN_LIMB_SIZE, 1);
    // glibc✔️✔️: 	  retval[RETURN_LIMB_SIZE - 1]
    // glibc✔️✔️: 	    |= ((mp_limb_t) 1) << ((MANT_DIG - 1) % BITS_PER_MP_LIMB);
    // glibc✔️✔️: 	}
    // glibc✔️✔️:       else if (exponent == DENORM_EXP
    // glibc✔️✔️: 	       && (retval[RETURN_LIMB_SIZE - 1]
    // glibc✔️✔️: 		   & (((mp_limb_t) 1) << ((MANT_DIG - 1) % BITS_PER_MP_LIMB)))
    // glibc✔️✔️: 	       != 0)
    // glibc✔️✔️: 	  /* The number was denormalized but now normalized.  */
    // glibc✔️✔️: 	exponent = MIN_EXP - 1;
    // glibc✔️✔️:     }
    // glibc✔️✔️:
    // glibc✔️✔️:   if (exponent >= MAX_EXP)
    // glibc✔️✔️:   overflow:
    // glibc✔️✔️:     return overflow_value (negative);
    // glibc✔️✔️:
    // glibc✔️✔️:   if (half_bit || more_bits_nonzero)
    // glibc✔️✔️:     {
    // glibc✔️✔️:       FLOAT force_inexact = (FLOAT) 1 + MIN_VALUE;
    // glibc✔️✔️:       math_force_eval (force_inexact);
    // glibc✔️✔️:     }
    // glibc✔️✔️:   return MPN2FLOAT (retval, exponent, negative);
    // glibc✔️✔️: }

    // Binary64 has RETURN_LIMB_SIZE=1 and fixed FE_TONEAREST. errno and
    // floating-environment flags are not part of the detached reader API.
    // Source shifts first, preserving old round bits as sticky, then rounds
    // exactly once. Integer arithmetic prevents intermediate double rounding.
    if exponent < -1022 {
        if exponent < -1075 {
            return signed_bits(0, negative);
        }
        let shift = (-1022 - exponent) as u32;
        let old_half_bit = round_limb & (1u64 << round_bit) != 0;
        more_bits |= round_limb & ((1u64 << round_bit) - 1) != 0;
        if shift == 53 {
            round_limb = retval;
            round_bit = 52;
            retval = 0;
        } else {
            round_limb = retval;
            round_bit = shift - 1;
            retval >>= shift;
        }
        more_bits |= old_half_bit;
        exponent = -1023;
    }
    if exponent >= 1024 {
        return signed_bits(0x7ff0_0000_0000_0000, negative);
    }
    let half_bit = round_limb & (1u64 << round_bit) != 0;
    let more_bits_nonzero = more_bits || round_limb & ((1u64 << round_bit) - 1) != 0;
    if round_away(retval & 1 != 0, half_bit, more_bits_nonzero) {
        retval += 1;
        if retval & (1u64 << 53) != 0 {
            exponent += 1;
            retval >>= 1;
            retval |= 1u64 << 52;
        } else if exponent == -1023 && retval & (1u64 << 52) != 0 {
            exponent = -1022;
        }
    }
    if exponent >= 1024 {
        return signed_bits(0x7ff0_0000_0000_0000, negative);
    }
    construct_double(retval, exponent, negative)
}

fn round_away(last_digit_odd: bool, half_bit: bool, more_bits: bool) -> bool {
    // glibc✔️✔️: static bool
    // glibc✔️✔️: round_away (bool negative, bool last_digit_odd, bool half_bit, bool more_bits,
    // glibc✔️✔️: 	    int mode)
    // glibc✔️✔️: {
    // glibc✔️✔️:   switch (mode)
    // glibc✔️✔️:     {
    // glibc✔️✔️:     case FE_DOWNWARD:
    // glibc✔️✔️:       return negative && (half_bit || more_bits);
    // glibc✔️✔️:
    // glibc✔️✔️:     case FE_TONEAREST:
    // glibc✔️✔️:       return half_bit && (last_digit_odd || more_bits);
    // glibc✔️✔️:
    // glibc✔️✔️:     case FE_TOWARDZERO:
    // glibc✔️✔️:       return false;
    // glibc✔️✔️:
    // glibc✔️✔️:     case FE_UPWARD:
    // glibc✔️✔️:       return !negative && (half_bit || more_bits);
    // glibc✔️✔️:
    // glibc✔️✔️:     default:
    // glibc✔️✔️:       abort ();
    // glibc✔️✔️:     }
    // glibc✔️✔️: }

    // The existing pinned coordinate profile fixes FE_TONEAREST.
    half_bit && (last_digit_odd || more_bits)
}

fn construct_double(frac: u64, exponent: i64, negative: bool) -> f64 {
    // glibc✔️✔️: double
    // glibc✔️✔️: __mpn_construct_double (mp_srcptr frac_ptr, int expt, int negative)
    // glibc✔️✔️: {
    // glibc✔️✔️:   union ieee754_double u;
    // glibc✔️✔️:
    // glibc✔️✔️:   u.ieee.negative = negative;
    // glibc✔️✔️:   u.ieee.exponent = expt + IEEE754_DOUBLE_BIAS;
    // glibc✔️✔️: #if BITS_PER_MP_LIMB == 32
    // glibc✔️✔️:   u.ieee.mantissa1 = frac_ptr[0];
    // glibc✔️✔️:   u.ieee.mantissa0 = frac_ptr[1] & (((mp_limb_t) 1
    // glibc✔️✔️: 				     << (DBL_MANT_DIG - 32)) - 1);
    // glibc✔️✔️: #elif BITS_PER_MP_LIMB == 64
    // glibc✔️✔️:   u.ieee.mantissa1 = frac_ptr[0] & (((mp_limb_t) 1 << 32) - 1);
    // glibc✔️✔️:   u.ieee.mantissa0 = (frac_ptr[0] >> 32) & (((mp_limb_t) 1
    // glibc✔️✔️: 					     << (DBL_MANT_DIG - 32)) - 1);
    // glibc✔️✔️: #else
    // glibc✔️✔️:   # error "mp_limb size " BITS_PER_MP_LIMB "not accounted for"
    // glibc✔️✔️: #endif
    // glibc✔️✔️:
    // glibc✔️✔️:   return u.d;
    // glibc✔️✔️: }

    // Source bitfield assignments, RETURN_LIMB_SIZE=1, no float arithmetic.
    signed_bits(
        ((exponent + 1023) as u64) << 52 | (frac & 0x000f_ffff_ffff_ffff),
        negative,
    )
}

fn signed_bits(bits: u64, negative: bool) -> f64 {
    f64::from_bits(bits | if negative { 1u64 << 63 } else { 0 })
}
