//! Crippen logP/MR parameter-table owner (RDKit `Crippen.cpp`).

use std::sync::Arc;

use crate::{DescriptorComputedState, DescriptorError, DescriptorInput, DescriptorResult};
use cosmolkit_search::{
    QueryGraph, SearchTarget, SmartsParseParams, SubstructMatchParams,
    build_prepared_query_match_context, parse_smarts,
    try_get_substruct_matches_with_params_and_context,
};

#[cfg(test)]
thread_local! {
    /// Per-thread count of SMARTS compiles performed by the Crippen
    /// pattern builders (mirrors `patterns::PATTERN_COMPILES`). Only
    /// observable in tests.
    pub(crate) static CRIPPEN_SMARTS_COMPILES: std::cell::Cell<usize> =
        const { std::cell::Cell::new(0) };
}

/// The pinned default parameter table (RDKit Crippen.cpp:192-319
/// `defaultParamData`), committed VERBATIM as a small production asset:
/// the exact header line, all 110 data rows in source order with
/// tab-separated keep-empty fields, the intentional O12-before-O7 and
/// S2-before-S1-before-S3 order flips, the nine blank MR cells and the
/// five notes rows. No value is derived from element arithmetic.
pub(crate) const DEFAULT_PARAM_DATA: &str = r#"#ID	SMARTS	logP	MR	Notes/Questions
C1	[CH4]	0.1441	2.503	
C1	[CH3]C	0.1441	2.503	
C1	[CH2](C)C	0.1441	2.503	
C2	[CH](C)(C)C	0	2.433	
C2	[C](C)(C)(C)C	0	2.433	
C3	[CH3][N,O,P,S,F,Cl,Br,I]	-0.2035	2.753	
C3	[CH2X4]([N,O,P,S,F,Cl,Br,I])[A;!#1]	-0.2035	2.753	
C4	[CH1X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])[A;!#1]	-0.2051	2.731	
C4	[CH0X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])([A;!#1])[A;!#1]	-0.2051	2.731	
C5	[C]=[!C;A;!#1]	-0.2783	5.007	
C6	[CH2]=C	0.1551	3.513	
C6	[CH1](=C)[A;!#1]	0.1551	3.513	
C6	[CH0](=C)([A;!#1])[A;!#1]	0.1551	3.513	
C6	[C](=C)=C	0.1551	3.513	
C7	[CX2]#[A;!#1]	0.0017	3.888	
C8	[CH3]c	0.08452	2.464	
C9	[CH3]a	-0.1444	2.412	
C10	[CH2X4]a	-0.0516	2.488	
C11	[CHX4]a	0.1193	2.582	
C12	[CH0X4]a	-0.0967	2.576	
C13	[cH0]-[A;!C;!N;!O;!S;!F;!Cl;!Br;!I;!#1]	-0.5443	4.041	
C14	[c][#9]	0	3.257	
C15	[c][#17]	0.245	3.564	
C16	[c][#35]	0.198	3.18	
C17	[c][#53]	0	3.104	
C18	[cH]	0.1581	3.35	
C19	[c](:a)(:a):a	0.2955	4.346	
C20	[c](:a)(:a)-a	0.2713	3.904	
C21	[c](:a)(:a)-C	0.136	3.509	
C22	[c](:a)(:a)-N	0.4619	4.067	
C23	[c](:a)(:a)-O	0.5437	3.853	
C24	[c](:a)(:a)-S	0.1893	2.673	
C25	[c](:a)(:a)=[C,N,O]	-0.8186	3.135	
C26	[C](=C)(a)[A;!#1]	0.264	4.305	
C26	[C](=C)(c)a	0.264	4.305	
C26	[CH1](=C)a	0.264	4.305	
C26	[C]=c	0.264	4.305	
C27	[CX4][A;!C;!N;!O;!P;!S;!F;!Cl;!Br;!I;!#1]	0.2148	2.693	
CS	[#6]	0.08129	3.243	
H1	[#1][#6,#1]	0.123	1.057	
H2	[#1]O[CX4,c]	-0.2677	1.395	
H2	[#1]O[!#6;!#7;!#8;!#16]	-0.2677	1.395	
H2	[#1][!#6;!#7;!#8]	-0.2677	1.395	
H3	[#1][#7]	0.2142	0.9627	
H3	[#1]O[#7]	0.2142	0.9627	
H4	[#1]OC=[#6,#7,O,S]	0.298	1.805	
H4	[#1]O[O,S]	0.298	1.805	
HS	[#1]	0.1125	1.112	
N1	[NH2+0][A;!#1]	-1.019	2.262	
N2	[NH+0]([A;!#1])[A;!#1]	-0.7096	2.173	
N3	[NH2+0]a	-1.027	2.827	
N4	[NH1+0]([!#1;A,a])a	-0.5188	3	
N5	[NH+0]=[!#1;A,a]	0.08387	1.757	
N6	[N+0](=[!#1;A,a])[!#1;A,a]	0.1836	2.428	
N7	[N+0]([A;!#1])([A;!#1])[A;!#1]	-0.3187	1.839	
N8	[N+0](a)([!#1;A,a])[A;!#1]	-0.4458	2.819	
N8	[N+0](a)(a)a	-0.4458	2.819	
N9	[N+0]#[A;!#1]	0.01508	1.725	
N10	[NH3,NH2,NH;+,+2,+3]	-1.95		
N11	[n+0]	-0.3239	2.202	
N12	[n;+,+2,+3]	-1.119		
N13	[NH0;+,+2,+3]([A;!#1])([A;!#1])([A;!#1])[A;!#1]	-0.3396	0.2604	
N13	[NH0;+,+2,+3](=[A;!#1])([A;!#1])[!#1;A,a]	-0.3396	0.2604	
N13	[NH0;+,+2,+3](=[#6])=[#7]	-0.3396	0.2604	
N14	[N;+,+2,+3]#[A;!#1]	0.2887	3.359	
N14	[N;-,-2,-3]	0.2887	3.359	
N14	[N;+,+2,+3](=[N;-,-2,-3])=N	0.2887	3.359	
NS	[#7]	-0.4806	2.134	
O1	[o]	0.1552	1.08	
O2	[OH,OH2]	-0.2893	0.8238	
O3	[O]([A;!#1])[A;!#1]	-0.0684	1.085	
O4	[O](a)[!#1;A,a]	-0.4195	1.182	
O5	[O]=[#7,#8]	0.0335	3.367	
O5	[OX1;-,-2,-3][#7]	0.0335	3.367	
O6	[OX1;-,-2,-2][#16]	-0.3339	0.7774	
O6	[O;-0]=[#16;-0]	-0.3339	0.7774	
O12	[O-]C(=O)	-1.326		"order flip here intentional"
O7	[OX1;-,-2,-3][!#1;!N;!S]	-1.189	0	
O8	[O]=c	0.1788	3.135	
O9	[O]=[CH]C	-0.1526	0	
O9	[O]=C(C)([A;!#1])	-0.1526	0	
O9	[O]=[CH][N,O]	-0.1526	0	
O9	[O]=[CH2]	-0.1526	0	
O9	[O]=[CX2]=O	-0.1526	0	
O10	[O]=[CH]c	0.1129	0.2215	
O10	[O]=C([C,c])[a;!#1]	0.1129	0.2215	
O10	[O]=C(c)[A;!#1]	0.1129	0.2215	
O11	[O]=C([!#1;!#6])[!#1;!#6]	0.4833	0.389	
OS	[#8]	-0.1188	0.6865	
F	[#9-0]	0.4202	1.108	
Cl	[#17-0]	0.6895	5.853	
Br	[#35-0]	0.8456	8.927	
I	[#53-0]	0.8857	14.02	
Hal	[#9,#17,#35,#53;-]	-2.996		
Hal	[#53;+,+2,+3]	-2.996		
Hal	[+;#3,#11,#19,#37,#55]	-2.996		"Footnote h indicates these should be here?"
P	[#15]	0.8612	6.92	
S2	[S;-,-2,-3,-4,+1,+2,+3,+5,+6]	-0.0024	7.365	"Order flip here is intentional"
S2	[S-0]=[N,O,P,S]	-0.0024	7.365	"Expanded definition of (pseudo-)ionic S"
S1	[S;A]	0.6482	7.591	"Order flip here is intentional"
S3	[s;a]	0.6237	6.691	
Me1	[#3,#11,#19,#37,#55]	-0.3808	5.754	
Me1	[#4,#12,#20,#38,#56]	-0.3808	5.754	
Me1	[#5,#13,#31,#49,#81]	-0.3808	5.754	
Me1	[#14,#32,#50,#82]	-0.3808	5.754	
Me1	[#33,#51,#83]	-0.3808	5.754	
Me1	[#34,#52,#84]	-0.3808	5.754	
Me2	[#21,#22,#23,#24,#25,#26,#27,#28,#29,#30]	-0.0025		
Me2	[#39,#40,#41,#42,#43,#44,#45,#46,#47,#48]	-0.0025		
Me2	[#72,#73,#74,#75,#76,#77,#78,#79,#80]	-0.0025		
"#;

/// One parsed default-table row (the source `CrippenParams` minus the
/// SMARTS pattern, which is not compiled at this boundary).
#[derive(Clone, Debug, PartialEq)]
pub struct CrippenParamRow {
    /// 0-based position in the table (source `idx`).
    pub idx: u32,
    /// Source label (`C1`, `H2`, `O12`, ...).
    pub label: &'static str,
    /// Source SMARTS text, uncompiled.
    pub smarts: &'static str,
    /// logP contribution (blank cell parses to 0.0).
    pub logp: f64,
    /// MR contribution (blank or unparseable cell parses to 0.0).
    pub mr: f64,
}

/// Parse one Crippen LogP numeric cell (source Crippen.cpp:170-174).
///
/// Behavior review: an EMPTY cell stores 0.0; a valid cell parses to its
/// value; a MALFORMED cell FAILS with a typed error preserving the
/// offending cell — the source `boost::lexical_cast<double>` throw has
/// NO catch on the LogP side, so no fallback exists here.
///
/// Complexity review: O(cell length), zero allocation on the success
/// paths; the error path allocates the preserved cell payload.
pub(crate) fn parse_crippen_logp(cell: &str) -> Result<f64, DescriptorError> {
    // RDKit source (Crippen.cpp:170-174):
    //   if (*token != "") {
    //     paramObj.logp = boost::lexical_cast<double>(*token);
    //   } else {
    //     paramObj.logp = 0.0;
    //   }
    // RDKit✔️✔️:   if (*token != "") {
    // RDKit✔️✔️:     paramObj.logp = boost::lexical_cast<double>(*token);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     paramObj.logp = 0.0;
    // RDKit✔️✔️:   }
    // Typed mapping: nonempty cells delegate to the ONE shared Boost
    // numeric owner; its failure maps to the typed
    // `CrippenParamNumeric` error preserving the ORIGINAL cell — the
    // uncaught bad_lexical_cast throw, NO fallback. The error is
    // constructed LAZILY (ok_or_else), so the success path never
    // allocates the owned cell payload.
    if cell.is_empty() {
        Ok(0.0)
    } else {
        crippen_numeric_cell(cell).ok_or_else(|| DescriptorError::CrippenParamNumeric {
            field: "logp",
            cell: cell.to_string(),
        })
    }
}

/// ONE shared Crippen numeric-cell conversion for NONEMPTY cells,
/// reproducing the pinned Boost `lexical_cast` double semantics
/// (boostorg/lexical_cast commit 02e5821, BSL-1.0):
/// `parse_inf_nan_impl` (inf_nan.hpp:45-100) runs FIRST over borrowed
/// bytes; unmatched input falls to the full-string decimal primitive,
/// which must yield a FINITE value (overflow to Inf rejects; finite
/// underflow accepts as signed zero).
///
/// Behavior review: consumes at most one leading `+`/`-`; ASCII
/// case-insensitive EXACT `inf` (len 3) or `infinity` (len 8) after the
/// sign; `nan` with either nothing after or `(...)` checked ONLY at the
/// opening position and the final byte — payload bytes between are
/// ignored, never interpreted, restricted or normalized (embedded NUL
/// and multibyte bytes retained). NaN is the canonical quiet NaN:
/// positive bits 0x7ff8000000000000, negative 0xfff8000000000000
/// (source `copysign(quiet_NaN, -1)`); Inf/zero signs preserved.
/// Special-token failures (e.g. `nan(x`, `infix`) FALL THROUGH to the
/// decimal primitive exactly as Boost falls through to its stream
/// parser, which then rejects them. Never trims, never accepts numeric
/// prefixes, no second decimal parser, no fitting.
///
/// Complexity review: O(cell bytes) over borrowed slices; per-byte
/// ASCII-case comparisons with no lowercase allocation; no heap use on
/// ANY path (the LogP wrapper's typed error owns the original cell;
/// success and MR paths build nothing).
fn crippen_numeric_cell(cell: &str) -> Option<f64> {
    // Boost lexical_cast, inf_nan.hpp (boostorg/lexical_cast commit
    // 02e5821ab32c45fad719829e9644e5d681c9ba0b, BSL-1.0) — original
    // file notice:
    //   Copyright Kevlin Henney, 2000-2005.
    //   Copyright Alexander Nasonov, 2006-2010.
    //   Copyright Antony Polukhin, 2011-2024.
    //
    //   Distributed under the Boost Software License, Version 1.0. (See
    //   accompanying file LICENSE_1_0.txt or copy at
    //   http://www.boost.org/LICENSE_1_0.txt)
    //
    // Boost source (inf_nan.hpp:41-100, verbatim):
    //   template <class CharT>
    //   bool lc_iequal(const CharT* val, const CharT* lcase, const CharT* ucase, unsigned int len) noexcept {
    //       for( unsigned int i=0; i < len; ++i ) {
    //           if ( val[i] != lcase[i] && val[i] != ucase[i] ) return false;
    //       }
    //
    //       return true;
    //   }
    //
    //   /* Returns true and sets the correct value if found NaN or Inf. */
    //   template <class CharT, class T>
    //   inline bool parse_inf_nan_impl(const CharT* begin, const CharT* end, T& value
    //       , const CharT* lc_NAN, const CharT* lc_nan
    //       , const CharT* lc_INFINITY, const CharT* lc_infinity
    //       , const CharT opening_brace, const CharT closing_brace) noexcept
    //   {
    //       if (begin == end) return false;
    //       const CharT minus = lcast_char_constants<CharT>::minus;
    //       const CharT plus = lcast_char_constants<CharT>::plus;
    //       const int inifinity_size = 8; // == sizeof("infinity") - 1
    //
    //       /* Parsing +/- */
    //       bool const has_minus = (*begin == minus);
    //       if (has_minus || *begin == plus) {
    //           ++ begin;
    //       }
    //
    //       if (end - begin < 3) return false;
    //       if (lc_iequal(begin, lc_nan, lc_NAN, 3)) {
    //           begin += 3;
    //           if (end != begin) {
    //               /* It is 'nan(...)' or some bad input*/
    //
    //               if (end - begin < 2) return false; // bad input
    //               -- end;
    //               if (*begin != opening_brace || *end != closing_brace) return false; // bad input
    //           }
    //
    //           if( !has_minus ) value = std::numeric_limits<T>::quiet_NaN();
    //           else value = boost::core::copysign(std::numeric_limits<T>::quiet_NaN(), static_cast<T>(-1));
    //           return true;
    //       } else if (
    //           ( /* 'INF' or 'inf' */
    //             end - begin == 3      // 3 == sizeof('inf') - 1
    //             && lc_iequal(begin, lc_infinity, lc_INFINITY, 3)
    //           )
    //           ||
    //           ( /* 'INFINITY' or 'infinity' */
    //             end - begin == inifinity_size
    //             && lc_iequal(begin, lc_infinity, lc_INFINITY, inifinity_size)
    //           )
    //        )
    //       {
    //           if( !has_minus ) value = std::numeric_limits<T>::infinity();
    //           else value = -std::numeric_limits<T>::infinity();
    //           return true;
    //       }
    //
    //       return false;
    //   }
    // Boost✔️✔️: the borrowed-byte scan below reproduces the exact
    // block above — one leading +/-, exact-length case-insensitive
    // nan/inf/infinity recognition, the nan(...) brace checks at the
    // opening and final byte with payload bytes ignored, the
    // quiet_NaN/copysign value selection, and the false-return
    // fall-through to the decimal primitive (the Boost stream parser).
    // Behavior review: special tokens map to the pinned binary64 NaN
    // bit constants and +/- infinity; malformed input NEVER matches and
    // falls through unchanged. Complexity review: O(cell bytes) over
    // borrowed slices, per-byte ASCII-case compares, no heap use on any
    // path of this owner.
    // (CRIPPEN-ANCHOR2: the earlier second condensed marked copy of the
    // helper bodies was removed; this single marker is the only one.)
    let bytes = cell.as_bytes();
    let mut rest = bytes;
    let negative = match rest.first() {
        Some(b'+') => {
            rest = &rest[1..];
            false
        }
        Some(b'-') => {
            rest = &rest[1..];
            true
        }
        _ => false,
    };
    if rest.len() >= 3 {
        if rest[..3].eq_ignore_ascii_case(b"nan") {
            let tail = &rest[3..];
            let matched = tail.is_empty()
                || (tail.len() >= 2 && tail[0] == b'(' && tail[tail.len() - 1] == b')');
            if matched {
                return Some(f64::from_bits(if negative {
                    0xfff8_0000_0000_0000
                } else {
                    0x7ff8_0000_0000_0000
                }));
            }
        } else if rest.len() == 3 && rest.eq_ignore_ascii_case(b"inf")
            || rest.len() == 8 && rest.eq_ignore_ascii_case(b"infinity")
        {
            return Some(if negative {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            });
        }
    }
    // Unmatched input falls to the full-string decimal primitive; only
    // FINITE results accept (overflow to Inf rejects, matching the
    // pinned Boost conversion failure; finite underflow keeps its sign).
    let value = cell.parse::<f64>().ok()?;
    if value.is_infinite() {
        return None;
    }
    Some(value)
}

/// Parse one Crippen MR numeric cell (source Crippen.cpp:175-183).
///
/// Behavior review: an EMPTY cell stores 0.0; a valid cell parses to its
/// value; a MALFORMED cell ALSO stores 0.0 — the source's own
/// `catch (boost::bad_lexical_cast&)` fallback, which exists ONLY on
/// the MR side.
///
/// Complexity review: O(cell length), zero allocation on every path.
pub(crate) fn parse_crippen_mr(cell: &str) -> f64 {
    // RDKit source (Crippen.cpp:175-183):
    //   if (*token != "") {
    //     try {
    //       paramObj.mr = boost::lexical_cast<double>(*token);
    //     } catch (boost::bad_lexical_cast &) {
    //       paramObj.mr = 0.0;
    //     }
    //   } else {
    //     paramObj.mr = 0.0;
    //   }
    // RDKit✔️✔️:   if (*token != "") {
    // RDKit✔️✔️:     try {
    // RDKit✔️✔️:       paramObj.mr = boost::lexical_cast<double>(*token);
    // RDKit✔️✔️:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️✔️:       paramObj.mr = 0.0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     paramObj.mr = 0.0;
    // RDKit✔️✔️:   }
    // Typed mapping: nonempty cells delegate to the ONE shared Boost
    // numeric owner; its failure becomes +0.0 — the source's catch IS
    // the only fallback (POSITIVE zero regardless of input sign).
    if cell.is_empty() {
        0.0
    } else {
        crippen_numeric_cell(cell).unwrap_or(0.0)
    }
}

/// Compile every SMARTS cell of a Crippen parameter table in row order,
/// retaining one parsed graph per row (source constructor, Crippen.cpp:184-186).
///
/// Behavior review: EVERY data row's SMARTS is compiled exactly once AT
/// COLLECTION CONSTRUCTION, in row order, and retained; a row whose
/// SMARTS fails to compile keeps `None` — the source `SmartsToMol`
/// nullptr retained WITHOUT any check or error at acquisition (failure
/// would only surface at later match use). No numeric parsing here:
/// only the label/smarts cells are consumed.
///
/// Complexity review: O(rows) SMARTS compiles, one per row, each the
/// audited `parse_smarts` cost; storage is one `Option<Arc<QueryGraph>>`
/// per row; no map, no dedup — the source constructs one pattern per
/// row, not a shared flyweight.
fn compile_crippen_patterns(param_data: &str) -> Vec<Option<Arc<QueryGraph>>> {
    // RDKit source (Crippen.cpp:184-186):
    //   paramObj.dp_pattern =
    //       boost::shared_ptr<const ROMol>(SmartsToMol(paramObj.smarts));
    //   d_params.push_back(paramObj);
    // RDKit✔️✔️:   paramObj.dp_pattern =
    // RDKit✔️✔️:       boost::shared_ptr<const ROMol>(SmartsToMol(paramObj.smarts));
    // RDKit✔️✔️:   d_params.push_back(paramObj);
    // Typed mapping: SmartsToMol's nullptr-on-failure (unchecked by the
    // constructor) maps to `None` retained per row; success wraps the
    // parsed QueryGraph in a shared Arc. The rows themselves come from
    // the same '#'-only line scan the source loop performs.
    let mut patterns = Vec::new();
    // RDKit source (Crippen.cpp:156-159, the line selection):
    //   std::string inLine = RDKit::getLine(inStream);
    //   unsigned int idx = 0;
    //   while (!inStream.eof()) {
    //     if (inLine[0] != '#') {
    // RDKit✔️✔️:   std::string inLine = RDKit::getLine(inStream);
    // RDKit✔️✔️:   unsigned int idx = 0;
    // RDKit✔️✔️:   while (!inStream.eof()) {
    // RDKit✔️✔️:     if (inLine[0] != '#') {
    // Typed mapping: the FIRST line enters the loop like every other
    // line and is skipped ONLY when its first byte is '#' — a
    // headerless custom table therefore KEEPS its first data row
    // (Crippen.cpp:145-188). No header heuristic, no line
    // normalization; the committed DEFAULT table's header line begins
    // with '#', so its behavior is unchanged.
    for line in param_data.lines() {
        if line.starts_with('#') {
            continue;
        }
        let smarts = line.split('\t').nth(1).unwrap_or("");
        #[cfg(test)]
        CRIPPEN_SMARTS_COMPILES.with(|count| count.set(count.get() + 1));
        patterns.push(
            parse_smarts(smarts, &SmartsParseParams::default())
                .ok()
                .map(Arc::new),
        );
    }
    patterns
}

/// Test-only actual child-call input observation; no production storage.
#[cfg(test)]
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct HydrogenWorkingCopyObservation {
    pub atom_count: usize,
    pub bond_count: usize,
    pub atom_membership_count: usize,
    pub bond_membership_count: usize,
    pub original_atom_ids: Vec<cosmolkit_model::AtomId>,
    pub original_bond_ids: Vec<cosmolkit_model::BondId>,
    pub original_isotopes: Vec<Option<u16>>,
    pub atom_rings: Vec<Vec<cosmolkit_model::AtomId>>,
    pub bond_rings: Vec<Vec<cosmolkit_model::BondId>>,
    pub child_state: DescriptorComputedState,
}

#[cfg(test)]
thread_local! {
    /// Cold-kernel invocations (test observability only, mirrors
    /// `patterns::PATTERN_COMPILES`).
    pub(crate) static CRIPPEN_KERNEL_CALLS: std::cell::Cell<usize> =
        const { std::cell::Cell::new(0) };
    /// Prepared matching contexts built inside the cold kernel.
    pub(crate) static CRIPPEN_CONTEXT_CALLS: std::cell::Cell<usize> =
        const { std::cell::Cell::new(0) };
    /// Total packed words across cold invocations (actual-site:
    /// incremented by atom_needed.len() after the final mask).
    pub(crate) static CRIPPEN_PACKED_WORDS: std::cell::Cell<usize> =
        const { std::cell::Cell::new(0) };
    /// The last initial packed word of the most recent cold invocation
    /// (actual-site: overwritten from atom_needed.last() after the mask).
    pub(crate) static CRIPPEN_LAST_WORD: std::cell::Cell<Option<u64>> =
        const { std::cell::Cell::new(None) };
    /// Parameter-row visits immediately before each match call.
    pub(crate) static CRIPPEN_ROW_VISITS: std::cell::Cell<usize> =
        const { std::cell::Cell::new(0) };
    /// Label writes inside the actual Some(labels) assignment.
    pub(crate) static CRIPPEN_LABEL_WRITES: std::cell::Cell<usize> =
        const { std::cell::Cell::new(0) };
    /// The ACTUAL extended input observed immediately before the
    /// canonical child contribution call inside the hydrogenated
    /// working-copy owner (cfg(test) observability only).
    pub(crate) static HYDROGEN_WORKING_COPY_OBSERVATION:
        std::cell::RefCell<Option<HydrogenWorkingCopyObservation>> =
        const { std::cell::RefCell::new(None) };
}

/// The ONE private cold first-match kernel (source
/// `getCrippenAtomContribs` cold region, Crippen.cpp:56-84).
///
/// Behavior review: zero-initializes ONLY the numeric output rows;
/// pre-sized optional sinks are NOT cleared — source-untyped atoms
/// keep their caller sentinel values. The packed `Vec<u64>` bitset
/// reproduces `boost::dynamic_bitset<> atomNeeded` with
/// `n.div_ceil(64)` words and a masked final word (bits beyond the
/// atom count stay zero so the all-words-zero early break is exact);
/// first-match marking, the source loop order and the
/// `atomNeeded.none()` early break are unchanged; the match
/// parameters keep the source defaults uniquify=false,
/// recursionPossible=true (maxMatches default 1000). A stored default
/// pattern that is unexpectedly None returns the typed
/// `MissingCrippenDefaultPattern` error (all 110 are source-proven
/// Some; the Option store retains no parser cause to preserve) — it is
/// never silently skipped or recompiled.
///
/// Complexity review: O(table-rows-until-all-typed x per-row match
/// cost) with ONE borrowed prepared matching context per cold
/// invocation (the audited equivalent of the source reusing
/// ROMol-internal query state); per-match work is the O(1)
/// first-atom projection plus O(1) bit twiddling; allocations are the
/// bitset words plus the match vectors the matcher itself owns. Label
/// strings are written ONLY into caller-supplied sinks, reusing each
/// String's existing capacity via clear()+push_str with NO temporary
/// String; absent sinks cause NO label buffer/string allocation or
/// growth. HONEST limits (corrected per CRIPPEN-KERNEL-PROOF): a
/// requested sink whose capacity is INSUFFICIENT for a label may grow
/// (reallocation is then the Rust allocator's, not a new owned String
/// per write); a sufficiently sized sink's capacity and pointer are
/// retained through the write. This does NOT claim C++ SSO behavior
/// or whole-matcher allocation parity.
pub(crate) fn crippen_cold_kernel(
    input: &DescriptorInput<'_>,
    logp_rows: &mut [f64],
    mr_rows: &mut [f64],
    mut atom_types: Option<&mut Vec<u32>>,
    mut atom_labels: Option<&mut Vec<String>>,
) -> DescriptorResult<()> {
    #[cfg(test)]
    CRIPPEN_KERNEL_CALLS.with(|count| count.set(count.get() + 1));
    let n_atoms = input.topology().atoms.len();
    debug_assert_eq!(logp_rows.len(), n_atoms);
    debug_assert_eq!(mr_rows.len(), n_atoms);
    logp_rows.fill(0.0);
    mr_rows.fill(0.0);
    // RDKit source (Crippen.cpp:56-84, verbatim; the single marker
    // below covers this whole block — no condensed duplicate copy):
    //   boost::dynamic_bitset<> atomNeeded(mol.getNumAtoms());
    //   atomNeeded.set();
    //   const CrippenParamCollection *params = CrippenParamCollection::getParams();
    //   for (const auto &param : *params) {
    //     std::vector<MatchVectType> matches;
    //     SubstructMatch(mol, *(param.dp_pattern.get()), matches, false, true);
    //     for (std::vector<MatchVectType>::const_iterator matchIt = matches.begin();
    //          matchIt != matches.end(); ++matchIt) {
    //       int idx = (*matchIt)[0].second;
    //       if (atomNeeded[idx]) {
    //         atomNeeded[idx] = 0;
    //         logpContribs[idx] = param.logp;
    //         mrContribs[idx] = param.mr;
    //         if (atomTypes) {
    //           (*atomTypes)[idx] = param.idx;
    //         }
    //         if (atomTypeLabels) {
    //           (*atomTypeLabels)[idx] = param.label;
    //         }
    //       }
    //     }
    //     // no need to keep matching stuff if we already found all the atoms:
    //     if (atomNeeded.none()) {
    //       break;
    //     }
    //   }
    // RDKit✔️✔️: the packed `Vec<u64>` words with a masked final word
    // reproduce dynamic_bitset(nAtoms)/set(); the loop order, the
    // first-match marking, the optional-sink writes and the early
    // break map line-for-line; `(*matchIt)[0].second` is the FIRST
    // QUERY ATOM's target index (`atom_mapping[0]`).
    let words = n_atoms.div_ceil(64);
    let mut atom_needed = vec![u64::MAX; words];
    if n_atoms % 64 != 0 {
        atom_needed[words - 1] = (1u64 << (n_atoms % 64)) - 1;
    }
    #[cfg(test)]
    {
        CRIPPEN_PACKED_WORDS.with(|c| c.set(c.get() + atom_needed.len()));
        CRIPPEN_LAST_WORD.with(|c| c.set(atom_needed.last().copied()));
    }
    let rows: &[CrippenParamRow] = default_crippen_params();
    let patterns: &[Option<Arc<QueryGraph>>] = crippen_patterns();
    let target = SearchTarget::new(
        input.topology(),
        input.coordinates(),
        &input.topology().stereo_groups,
        Some(input.ring_info()),
        Some(input.valence()),
    );
    #[cfg(test)]
    CRIPPEN_CONTEXT_CALLS.with(|count| count.set(count.get() + 1));
    let context =
        build_prepared_query_match_context(input.topology(), input.ring_info(), input.valence())
            .map_err(|source| DescriptorError::Search {
                function: "crippen_cold_kernel",
                source: crate::DescriptorSearchCause::Context(source),
            })?;
    let match_params = SubstructMatchParams {
        uniquify: false,
        recursion_possible: true,
        max_matches: 1000,
        ..SubstructMatchParams::default()
    };
    for (row, pattern) in rows.iter().zip(patterns) {
        let Some(query) = pattern else {
            return Err(DescriptorError::MissingCrippenDefaultPattern { row: row.idx });
        };
        #[cfg(test)]
        CRIPPEN_ROW_VISITS.with(|count| count.set(count.get() + 1));
        let matches = try_get_substruct_matches_with_params_and_context(
            &target,
            query,
            &match_params,
            &context,
        )
        .map_err(|source| DescriptorError::Search {
            function: "crippen_cold_kernel",
            source: crate::DescriptorSearchCause::Match(source),
        })?;
        for m in &matches {
            let idx = m.atom_mapping[0];
            if atom_needed[idx / 64] & (1u64 << (idx % 64)) != 0 {
                atom_needed[idx / 64] &= !(1u64 << (idx % 64));
                logp_rows[idx] = row.logp;
                mr_rows[idx] = row.mr;
                if let Some(types) = atom_types.as_deref_mut() {
                    types[idx] = row.idx;
                }
                if let Some(labels) = atom_labels.as_deref_mut() {
                    let target = &mut labels[idx];
                    target.clear();
                    target.push_str(row.label);
                    #[cfg(test)]
                    CRIPPEN_LABEL_WRITES.with(|count| count.set(count.get() + 1));
                }
            }
        }
        if atom_needed.iter().all(|&word| word == 0) {
            break;
        }
    }
    Ok(())
}

/// The ONE private contribution guard (source
/// `getCrippenAtomContribs` entry, Crippen.cpp:35-86).
///
/// Behavior review: the source PRECONDITIONs on the OPTIONAL sink
/// sizes are checked FIRST — before ANY cache access; the warm arm
/// reads the LogP rows and then the MR rows, failing with the typed
/// `MissingCrippenMrContributions` when the MR rows are absent BEFORE
/// either row length is examined; a hit (both lengths equal the atom
/// count) returns COPIES of the cached rows with the optional sinks
/// UNTOUCHED and performs no chemistry, preparation, context or kernel
/// work; every other combination falls through to the ONE cold kernel
/// (which zero-initializes only the numeric rows and never clears the
/// pre-sized sinks), and BOTH row vectors are published to the
/// Crippen cache slot only after complete successful matching — never
/// scalars, never unrelated slots, and errors preserve prior state.
///
/// Complexity review: the optional-sink checks and the warm arm are
/// O(1) presence checks plus the cached-row CLONES (the source's
/// `getProp` vector copies) — no context, no valence, no kernel on a
/// hit; the cold arm is the kernel's audited cost plus the two
/// publication clones.
pub(crate) fn crippen_contribution_guard(
    input: &DescriptorInput<'_>,
    force: bool,
    state: &mut DescriptorComputedState,
    mut atom_types: Option<&mut Vec<u32>>,
    atom_labels: Option<&mut Vec<String>>,
) -> DescriptorResult<(Vec<f64>, Vec<f64>)> {
    let n_atoms = input.topology().atoms.len();
    // RDKit source (Crippen.cpp:36-42, verbatim):
    //   PRECONDITION(logpContribs.size() == mol.getNumAtoms() &&
    //                    mrContribs.size() == mol.getNumAtoms(),
    //                "bad result vector size");
    //   PRECONDITION((!atomTypes || atomTypes->size() == mol.getNumAtoms()),
    //                "bad atomTypes vector");
    //   PRECONDITION((!atomTypeLabels || atomTypeLabels->size() == mol.getNumAtoms()),
    //                "bad atomTypeLabels vector");
    // RDKit✔️✔️: the two OPTIONAL-sink PRECONDITIONs map to the typed
    // InvalidCrippenOptionalRows error checked BEFORE any cache access;
    // the numeric-row PRECONDITION holds by construction (the guard
    // allocates them itself below).
    if let Some(types) = atom_types.as_deref() {
        if types.len() != n_atoms {
            return Err(DescriptorError::InvalidCrippenOptionalRows {
                field: "atom_types",
                actual: types.len(),
                expected: n_atoms,
            });
        }
    }
    if let Some(labels) = atom_labels.as_deref() {
        if labels.len() != n_atoms {
            return Err(DescriptorError::InvalidCrippenOptionalRows {
                field: "atom_labels",
                actual: labels.len(),
                expected: n_atoms,
            });
        }
    }
    // RDKit source (Crippen.cpp:43-53, verbatim):
    //   if (!force && mol.hasProp(common_properties::_crippenLogPContribs)) {
    //     std::vector<double> tmpVect1, tmpVect2;
    //     mol.getProp(common_properties::_crippenLogPContribs, tmpVect1);
    //     mol.getProp(common_properties::_crippenMRContribs, tmpVect2);
    //     if (tmpVect1.size() == mol.getNumAtoms() &&
    //         tmpVect2.size() == mol.getNumAtoms()) {
    //       logpContribs = tmpVect1;
    //       mrContribs = tmpVect2;
    //       return;
    //     }
    //   }
    // RDKit✔️✔️: the MR-rows read failure is the typed
    // MissingCrippenMrContributions error raised BEFORE either length
    // check; a length-mismatched pair falls through to the cold kernel
    // exactly as the source falls through to recompute.
    if !force {
        if let Some(lp_rows) = state.crippen_slot().logp_rows.clone() {
            let mr_rows = state.crippen_slot().mr_rows.clone().ok_or(
                DescriptorError::MissingCrippenMrContributions {
                    function: "crippen_contributions",
                },
            )?;
            if lp_rows.len() == n_atoms && mr_rows.len() == n_atoms {
                return Ok((lp_rows, mr_rows));
            }
        }
    }
    // RDKit source (Crippen.cpp:56-85, verbatim): the cold region and
    // the row publication are the ONE cold kernel plus the two-field
    // publication below (see crippen_cold_kernel for the full block).
    // RDKit✔️✔️: publication writes BOTH row vectors after complete
    // successful matching only.
    let mut logp_rows = vec![0.0; n_atoms];
    let mut mr_rows = vec![0.0; n_atoms];
    crippen_cold_kernel(input, &mut logp_rows, &mut mr_rows, atom_types, atom_labels)?;
    let slot = state.crippen_slot_mut();
    slot.logp_rows = Some(logp_rows.clone());
    slot.mr_rows = Some(mr_rows.clone());
    Ok((logp_rows, mr_rows))
}

/// The canonical per-atom Crippen contribution result (the source's
/// `logpContribs`/`mrContribs` output vectors only — the optional
/// type/label projections travel through caller-owned sinks).
#[derive(Clone, Debug, PartialEq)]
pub struct CrippenContributions {
    /// Per-atom logP contributions.
    pub logp: Vec<f64>,
    /// Per-atom molar-refractivity contributions.
    pub molar_refractivity: Vec<f64>,
}

/// Canonical Crippen atom contributions (source
/// `getCrippenAtomContribs`, Crippen.cpp:35-86).
///
/// Behavior review: delegates the WHOLE entry — the optional-sink
/// PRECONDITIONs before any cache access, the non-includeHs-keyed warm
/// guard with the MR-missing precedence, the exact hit/cold
/// discrimination and the ONE cold kernel — to the shared private
/// [`crippen_contribution_guard`], and publishes ONLY the two
/// contribution fields into the Crippen cache slot after successful
/// matching (never scalars, never unrelated slots); the optional
/// type/label sinks travel through untouched as in the source.
///
/// Complexity review: the canonical entry adds zero work over the
/// guard (a result-struct move; the rows are the guard's owned
/// outputs).
///
/// # Exact canonical access
///
/// ```
/// use cosmolkit_descriptors::{
///     CrippenContributions, DescriptorComputedState, DescriptorInput,
///     DescriptorResult, crippen_contributions,
/// };
///
/// let contributions: fn(
///     &DescriptorInput<'_>,
///     bool,
///     &mut DescriptorComputedState,
///     Option<&mut Vec<u32>>,
///     Option<&mut Vec<String>>,
/// ) -> DescriptorResult<CrippenContributions> = crippen_contributions;
/// ```
///
/// The implementation module and its private helpers stay private:
///
/// ```compile_fail,E0603
/// // Private implementation-module access: `crippen` is not a public
/// // module; `crippen_cold_kernel` and `crippen_contribution_guard`
/// // are crate-private. This proof is module privacy only — not full
/// // operation/runtime isolation.
/// let f = cosmolkit_descriptors::crippen::crippen_cold_kernel;
/// let g = cosmolkit_descriptors::crippen::crippen_contribution_guard;
/// ```
pub fn crippen_contributions(
    input: &DescriptorInput<'_>,
    force: bool,
    state: &mut DescriptorComputedState,
    atom_types: Option<&mut Vec<u32>>,
    atom_labels: Option<&mut Vec<String>>,
) -> DescriptorResult<CrippenContributions> {
    // RDKit source (Crippen.cpp:35-42, verbatim): the entry
    // PRECONDITIONs; RDKit✔️✔️: enforced by the shared guard BEFORE any
    // cache access (see crippen_contribution_guard for the exact
    // mapping to InvalidCrippenOptionalRows).
    // RDKit source (Crippen.cpp:43-53, verbatim): the warm guard;
    // RDKit✔️✔️: the guard's non-includeHs-keyed read with the
    // MR-missing typed failure and exact-length hit.
    // RDKit source (Crippen.cpp:84-85, verbatim):
    //   mol.setProp(common_properties::_crippenLogPContribs, logpContribs, true);
    //   mol.setProp(common_properties::_crippenMRContribs, mrContribs, true);
    // RDKit✔️✔️: the guard publishes BOTH row vectors after complete
    // successful matching only; this canonical entry adds nothing else.
    let (logp, molar_refractivity) =
        crippen_contribution_guard(input, force, state, atom_types, atom_labels)?;
    Ok(CrippenContributions {
        logp,
        molar_refractivity,
    })
}

/// The ONE private scalar cache guard (source
/// `calcCrippenDescriptors` entry, Crippen.cpp:87-93).
///
/// Behavior review: `!force` with the LogP scalar present reads the MR
/// scalar INDEPENDENTLY (its absence is the typed
/// `MissingCrippenMr` failure, preserving the source's second getProp
/// read failure) and returns the pair — BEFORE any includeHs handling
/// or prepared validation. Every other combination returns `None`,
/// leaving the cold branch (hydrogen option + contribution owner +
/// ordered sums + two-scalar publication after success) to the
/// canonical totals owner. Unrelated cache fields are never touched.
///
/// Complexity review: O(1) presence checks and scalar copies on the
/// warm arm — NO cached-vector clones, NO valence/context/copy work;
/// the guard itself allocates nothing.
pub(crate) fn crippen_scalar_guard(
    force: bool,
    state: &DescriptorComputedState,
) -> Result<Option<(f64, f64)>, DescriptorError> {
    // RDKit source (Crippen.cpp:87-93, verbatim):
    //   void calcCrippenDescriptors(const ROMol &mol, double &logp, double &mr,
    //                               bool includeHs, bool force) {
    //     if (!force && mol.hasProp(common_properties::_crippenLogP)) {
    //       mol.getProp(common_properties::_crippenLogP, logp);
    //       mol.getProp(common_properties::_crippenMR, mr);
    //       return;
    //     }
    //   }
    // RDKit✔️✔️: the warm arm maps the two INDEPENDENT scalar reads —
    // the MR read failure is the typed MissingCrippenMr error; the
    // return-before-anything-else maps to returning the pair with no
    // includeHs or prepared work; every other combination falls to the
    // cold branch (None).
    if !force {
        if let Some(logp) = state.crippen_slot().logp {
            let mr = state
                .crippen_slot()
                .mr
                .ok_or(DescriptorError::MissingCrippenMr {
                    function: "crippen_totals",
                })?;
            return Ok(Some((logp, mr)));
        }
    }
    Ok(None)
}

/// The ONE private hydrogenated working-copy owner (source
/// `calcCrippenDescriptors` includeHs branch, Crippen.cpp:94-100 with
/// AddHs.cpp:533-538 ring retention).
///
/// Behavior review: ONE core `add_hydrogens_with_params` on CLONED
/// detached blocks (explicit_only=false, add_coords=false; remaining
/// source defaults); the original input blocks are never mutated and
/// are asserted bit-identical by the caller-side tests. The appended
/// atom/bond IDs are verified to preserve every original index (the
/// working copy is a pure extension). HydrogenError propagates
/// STRUCTURALLY through `DescriptorError::Hydrogens` with the
/// borrowed `Error::source` — never flattened to text/Unsupported. A
/// FRESH cleared child `DescriptorComputedState` carries the extended
/// computation (AddHs clears computed properties — AddHs.cpp:533-538),
/// and the child is discarded after the ordered sums. The ring
/// context is the ORIGINAL paired atom/bond ring rows transported at
/// the extended counts via the EXISTING
/// `ring_info_from_selected_rows` — NO SSSR recomputation; that
/// constructor's OtherOrUnknown find type is a PRIVATE Crippen
/// membership projection (prepared matcher reads memberships/rows,
/// not a find-type acquisition decision), recorded as a limitation.
/// One final valence assignment prepares the extended topology. The
/// ordered sums fold left-to-right exactly as the source iterates the
/// contribution vectors; the PARENT receives only the two scalar
/// fields after complete success (parent contributions unchanged).
///
/// Complexity review: one AddHs transform (cloned blocks — the
/// source's owned working copy), one final valence assignment on the
/// extended topology, one ring-row transport (row copies only, no
/// finding), the ONE cold kernel through the canonical contributions
/// owner, and two O(nAtoms) ordered summations; no second SSSR pass.
pub(crate) fn crippen_hydrogenated_totals(
    input: &DescriptorInput<'_>,
    force: bool,
    parent_state: &mut DescriptorComputedState,
) -> DescriptorResult<(f64, f64)> {
    // RDKit source (Crippen.cpp:94-100, verbatim):
    //   // this isn't as bad as it looks, we aren't actually going
    //   // to harm the molecule in any way!
    //   auto *workMol = const_cast<ROMol *>(&mol);
    //   if (includeHs) {
    //     workMol = MolOps::addHs(mol, false, false);
    //   }
    // RDKit✔️✔️: the owned result blocks ARE the working copy; the
    // ORIGINAL input blocks are only read; MolOps::addHs(mol, false,
    // false) maps to add_hydrogens_with_params with
    // explicit_only=false, add_coords=false and the remaining source
    // defaults.
    let added = cosmolkit_core::add_hydrogens_with_params(
        input.topology().clone(),
        input.coordinates().clone(),
        input.properties().clone(),
        &cosmolkit_core::AddHsParams {
            explicit_only: false,
            add_coords: false,
            ..cosmolkit_core::AddHsParams::default()
        },
    )
    .map_err(|source| DescriptorError::Hydrogens {
        function: "crippen_totals",
        source,
    })?;
    // RDKit source (AddHs.cpp:533-538, verbatim):
    //   void addHs(RWMol &mol, const AddHsParameters &params,
    //              const UINT_VECT *onlyOnAtoms) {
    //     // when we hit each atom, clear its computed properties
    //     // NOTE: it is essential that we not clear the ring info in the
    //     // molecule's computed properties.  We don't want to have to
    //     // regenerate that.  This caused Issue210 and Issue212:
    //     mol.clearComputedProps(false);
    // RDKit✔️✔️: the fresh cleared child state carries the "cleared
    // computed properties" semantics; the ORIGINAL paired ring rows
    // are RETAINED (transported, never regenerated) exactly as the
    // source keeps the ring info.
    let original_atoms = input.topology().atoms.len();
    let original_bonds = input.topology().bonds.len();
    debug_assert!(added.topology.atoms.len() >= original_atoms);
    debug_assert!(added.topology.bonds.len() >= original_bonds);
    // The appended IDs preserve every original index: the first
    // original_atoms atoms and original_bonds bonds keep their
    // identity AND position (pure extension).
    for (index, atom) in added.topology.atoms.iter().enumerate().take(original_atoms) {
        debug_assert_eq!(atom.id(), input.topology().atoms[index].id());
    }
    for (index, bond) in added.topology.bonds.iter().enumerate().take(original_bonds) {
        debug_assert_eq!(bond.id(), input.topology().bonds[index].id());
    }
    // Retained paired ring-row transport at the EXTENDED counts —
    // ring_info_from_selected_rows is the EXISTING core constructor;
    // its OtherOrUnknown find type is the recorded private-membership
    // limitation, not an acquisition decision.
    let transported_rings = cosmolkit_core::ring_info_from_selected_rows(
        added.topology.atoms.len(),
        added.topology.bonds.len(),
        input.ring_info().atom_rings(),
        input.ring_info().bond_rings(),
    )
    .map_err(|source| DescriptorError::Ring {
        function: "crippen_totals",
        source,
    })?;
    // RDKit source (Crippen.cpp:101-117, verbatim):
    //   std::vector<double> logpContribs(workMol->getNumAtoms());
    //   std::vector<double> mrContribs(workMol->getNumAtoms());
    //   getCrippenAtomContribs(*workMol, logpContribs, mrContribs, force);
    //   logp = 0.0;
    //   for (std::vector<double>::const_iterator iter = logpContribs.begin();
    //        iter != logpContribs.end(); ++iter) {
    //     logp += *iter;
    //   }
    //   mr = 0.0;
    //   for (std::vector<double>::const_iterator iter = mrContribs.begin();
    //        iter != mrContribs.end(); ++iter) {
    //     mr += *iter;
    //   }
    //   if (includeHs) {
    //     delete workMol;
    //   }
    // RDKit✔️✔️: ONE final valence assignment prepares the extended
    // topology; the child state is FRESH (AddHs cleared computed
    // props); the contributions come from the canonical owner; the
    // ordered sums fold left-to-right; the child is discarded; the
    // parent publishes BOTH scalars only after success.
    let assignment = crate::valence(&added.topology, "crippen_totals")?;
    let extended_input = DescriptorInput::new(
        &added.topology,
        &added.coordinates,
        &added.properties,
        &assignment,
        &transported_rings,
    );
    let mut child_state = DescriptorComputedState::default();
    #[cfg(test)]
    {
        let observation = HydrogenWorkingCopyObservation {
            atom_count: extended_input.topology().atoms.len(),
            bond_count: extended_input.topology().bonds.len(),
            atom_membership_count: extended_input.ring_info().atom_row_count(),
            bond_membership_count: extended_input.ring_info().bond_row_count(),
            original_atom_ids: extended_input
                .topology()
                .atoms
                .iter()
                .take(original_atoms)
                .map(|a| a.id())
                .collect(),
            original_bond_ids: extended_input
                .topology()
                .bonds
                .iter()
                .take(original_bonds)
                .map(|b| b.id())
                .collect(),
            original_isotopes: extended_input
                .topology()
                .atoms
                .iter()
                .take(original_atoms)
                .map(|atom| atom.isotope())
                .collect(),
            atom_rings: extended_input.ring_info().atom_rings().to_vec(),
            bond_rings: extended_input.ring_info().bond_rings().to_vec(),
            child_state: child_state.clone(),
        };
        HYDROGEN_WORKING_COPY_OBSERVATION.with(|slot| {
            *slot.borrow_mut() = Some(observation);
        });
    }
    let contribs = crippen_contributions(&extended_input, force, &mut child_state, None, None)?;
    let mut logp = 0.0_f64;
    for value in &contribs.logp {
        logp += *value;
    }
    let mut mr = 0.0_f64;
    for value in &contribs.molar_refractivity {
        mr += *value;
    }
    let slot = parent_state.crippen_slot_mut();
    slot.logp = Some(logp);
    slot.mr = Some(mr);
    Ok((logp, mr))
}

/// The retained, once-compiled pattern collection of the DEFAULT Crippen
/// table (source `getParams("")` -> the defaultParamData flyweight entry,
/// Crippen.cpp:140-153). Index-aligned with [`default_crippen_params`].
///
/// Behavior review: the default collection compiles ALL 110 ordered
/// patterns exactly ONCE per process (flyweight construct-once); every
/// later acquisition is a retained borrow — no per-molecule and no
/// per-parameter reparse.
///
/// Complexity review: one-time O(rows) build inside a `OnceLock`; each
/// call is one atomic-load slice borrow with no Arc clone and no map
/// lookup — the same cost class as the flyweight handle.
pub(crate) fn crippen_patterns() -> &'static [Option<Arc<QueryGraph>>] {
    // RDKit source (Crippen.cpp:140-144):
    //   const CrippenParamCollection *CrippenParamCollection::getParams(
    //       const std::string &paramData) {
    //     const CrippenParamCollection *res = &(param_flyweight(paramData).get());
    //     return res;
    //   }
    // RDKit✔️✔️:   const CrippenParamCollection *CrippenParamCollection::getParams(
    // RDKit✔️✔️:       const std::string &paramData) {
    // RDKit✔️✔️:     const CrippenParamCollection *res = &(param_flyweight(paramData).get());
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    // Typed mapping: the "" flyweight entry maps to this process-global
    // OnceLock (construct once, borrow forever).
    // RDKit source (Crippen.cpp:149-153):
    //   if (paramData == "") {
    //     params = defaultParamData;
    //   } else {
    //     params = paramData;
    //   }
    // RDKit✔️✔️:   if (paramData == "") {
    // RDKit✔️✔️:     params = defaultParamData;
    // RDKit✔️✔️:   }
    static DEFAULT_PATTERNS: std::sync::OnceLock<Vec<Option<Arc<QueryGraph>>>> =
        std::sync::OnceLock::new();
    DEFAULT_PATTERNS.get_or_init(|| compile_crippen_patterns(DEFAULT_PARAM_DATA))
}

/// Compile the pattern collection of a CUSTOM Crippen parameter table
/// (source constructor on a non-empty `paramData`, Crippen.cpp:145-188).
/// The CALLER owns the returned collection — this boundary keeps no
/// second global registry keyed by custom strings; reparsing happens
/// only if the caller drops the result and calls again.
///
/// Behavior/complexity reviews: identical to
/// [`compile_crippen_patterns`]'s per-row semantics; one build per call,
/// O(rows) compiles.
pub(crate) fn crippen_patterns_from(param_data: &str) -> Vec<Option<Arc<QueryGraph>>> {
    // RDKit source (Crippen.cpp:151-152):
    //   } else {
    //     params = paramData;
    //   }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     params = paramData;
    // RDKit✔️✔️:   }
    compile_crippen_patterns(param_data)
}

/// The canonical Crippen scalar-totals result (source `logp`/`mr`
/// out-parameters of `calcCrippenDescriptors`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrippenTotals {
    /// The Wildman-Crippen logP estimate.
    pub logp: f64,
    /// The Wildman-Crippen molar-refractivity estimate.
    pub molar_refractivity: f64,
}

/// Canonical Crippen totals (source `calcCrippenDescriptors`,
/// Crippen.cpp:87-123).
///
/// Behavior review: the scalar guard runs FIRST (warm pair or typed
/// MissingCrippenMr — before any includeHs handling or preparation);
/// the no-H arm computes the canonical contributions on the ORIGINAL
/// input/state with left-to-right ordered sums; the with-H arm is
/// the ONE hydrogenated working-copy owner (fresh child state, ring
/// transport, structural HydrogenError). BOTH scalar fields publish
/// to the parent Crippen slot ONLY after complete success; the
/// parent's contribution fields and unrelated slots are retained.
///
/// Complexity review: the guard is O(1); the no-H arm is the
/// contribution owner's audited cost plus two O(nAtoms) folds; the
/// with-H arm is the working-copy owner's audited cost.
///
/// # Exact canonical access
///
/// ```
/// use cosmolkit_descriptors::{
///     CrippenTotals, DescriptorComputedState, DescriptorInput,
///     DescriptorResult, crippen_totals,
/// };
///
/// let totals: fn(
///     &DescriptorInput<'_>,
///     bool,
///     bool,
///     &mut DescriptorComputedState,
/// ) -> DescriptorResult<CrippenTotals> = crippen_totals;
/// ```
///
/// ```compile_fail,E0603
/// // Private implementation-module access: `crippen` is not a public
/// // module; `crippen_scalar_guard` and `crippen_hydrogenated_totals`
/// // are crate-private. This proof is module privacy only — not full
/// // operation/runtime isolation.
/// let f = cosmolkit_descriptors::crippen::crippen_scalar_guard;
/// let g = cosmolkit_descriptors::crippen::crippen_hydrogenated_totals;
/// ```
pub fn crippen_totals(
    input: &DescriptorInput<'_>,
    include_hydrogens: bool,
    force: bool,
    state: &mut DescriptorComputedState,
) -> DescriptorResult<CrippenTotals> {
    // RDKit source (Crippen.cpp:87-93, verbatim): the scalar warm
    // guard (anchored in crippen_scalar_guard) runs BEFORE the
    // includeHs branch; RDKit✔️✔️: delegated to the ONE guard.
    if let Some((logp, mr)) = crippen_scalar_guard(force, state)? {
        return Ok(CrippenTotals {
            logp,
            molar_refractivity: mr,
        });
    }
    if include_hydrogens {
        let (logp, mr) = crippen_hydrogenated_totals(input, force, state)?;
        Ok(CrippenTotals {
            logp,
            molar_refractivity: mr,
        })
    } else {
        // RDKit source (Crippen.cpp:101-117, verbatim — anchored in
        // the contributions owner and the working-copy owner):
        // RDKit✔️✔️: no-H arm — canonical contributions on the
        // ORIGINAL input/state, left-to-right ordered sums, both
        // scalars published only after success.
        let contribs = crippen_contributions(input, force, state, None, None)?;
        let mut logp = 0.0_f64;
        for value in &contribs.logp {
            logp += *value;
        }
        let mut mr = 0.0_f64;
        for value in &contribs.molar_refractivity {
            mr += *value;
        }
        let slot = state.crippen_slot_mut();
        slot.logp = Some(logp);
        slot.mr = Some(mr);
        Ok(CrippenTotals {
            logp,
            molar_refractivity: mr,
        })
    }
}

/// The logP half of the Crippen totals pair (the source `calcClogP`
/// projection, Crippen.cpp:123-128).
///
/// Behavior review: a pure projection of the ONE shared totals engine —
/// the source declares both locals, calls `calcCrippenDescriptors(mol,
/// clogp, mr)` with the DEFAULT arguments (includeHs=true, force=false
/// per the header defaults), and returns only clogp; the mr local is
/// discarded, exactly as the unreturned result field here. The explicit
/// mutable state carries the source's non-includeHs-keyed cache.
///
/// Complexity review: the shared `crippen_totals` path plus one field
/// projection; no additional allocation beyond that owner's work.
///
/// # Exact canonical access
///
/// ```
/// use cosmolkit_descriptors::{
///     DescriptorComputedState, DescriptorInput, DescriptorResult, crippen_clogp,
/// };
///
/// let projection: fn(
///     &DescriptorInput<'_>,
///     &mut DescriptorComputedState,
/// ) -> DescriptorResult<f64> = crippen_clogp;
/// ```
///
/// ```compile_fail,E0603
/// // Private implementation-module path: the crate-root re-export is
/// // the only public route to the projection.
/// let f = cosmolkit_descriptors::crippen::crippen_clogp;
/// ```
pub fn crippen_clogp(
    input: &DescriptorInput<'_>,
    state: &mut DescriptorComputedState,
) -> DescriptorResult<f64> {
    // RDKit source (Crippen.cpp:123-128):
    //   double calcClogP(const ROMol &mol) {
    //     double clogp, mr;
    //     calcCrippenDescriptors(mol, clogp, mr);
    //     return clogp;
    //   }
    // RDKit✔️✔️:   double calcClogP(const ROMol &mol) {
    // RDKit✔️✔️:     double clogp, mr;
    // RDKit✔️✔️:     calcCrippenDescriptors(mol, clogp, mr);
    // RDKit✔️✔️:     return clogp;
    // RDKit✔️✔️:   }
    // Typed mapping: calcCrippenDescriptors omits BOTH arguments, so
    // the pinned HEADER defaults apply — includeHs=TRUE, force=false
    // (Crippen.h:72-75); the projection carries the caller's explicit
    // mutable state; the discarded mr local maps to the unreturned
    // half of the ONE canonical totals owner.
    crippen_totals(input, true, false, state).map(|totals| totals.logp)
}

/// The MR half of the Crippen totals pair (the source `calcMR`
/// projection, Crippen.cpp:130-133).
///
/// Behavior review: identical projection shape returning the mr half.
///
/// Complexity review: identical to [`crippen_clogp`].
///
/// # Exact canonical access
///
/// ```
/// use cosmolkit_descriptors::{
///     DescriptorComputedState, DescriptorInput, DescriptorResult, crippen_mr,
/// };
///
/// let projection: fn(
///     &DescriptorInput<'_>,
///     &mut DescriptorComputedState,
/// ) -> DescriptorResult<f64> = crippen_mr;
/// ```
///
/// ```compile_fail,E0603
/// // Private implementation-module path: the crate-root re-export is
/// // the only public route to the projection.
/// let f = cosmolkit_descriptors::crippen::crippen_mr;
/// ```
pub fn crippen_mr(
    input: &DescriptorInput<'_>,
    state: &mut DescriptorComputedState,
) -> DescriptorResult<f64> {
    // RDKit source (Crippen.cpp:130-133):
    //   double calcMR(const ROMol &mol) {
    //     double clogp, mr;
    //     calcCrippenDescriptors(mol, clogp, mr);
    //     return mr;
    //   }
    // RDKit✔️✔️:   double calcMR(const ROMol &mol) {
    // RDKit✔️✔️:     double clogp, mr;
    // RDKit✔️✔️:     calcCrippenDescriptors(mol, clogp, mr);
    // RDKit✔️✔️:     return mr;
    // RDKit✔️✔️:   }
    crippen_totals(input, true, false, state).map(|totals| totals.molar_refractivity)
}

/// The ordered default Crippen parameter table (the source
/// `CrippenParamCollection::getParams("")` selecting `defaultParamData`,
/// Crippen.cpp:151-190).
///
/// Behavior review: reproduces the default-table row projection — the
/// header line is consumed first (`getLine` before the loop); every
/// subsequent line is tab-tokenized with EMPTY TOKENS KEPT; lines
/// starting with `#` are skipped; `idx` is a 0-based sequence over data
/// rows; a blank logP cell stores 0.0 while a malformed one aborts
/// (the source `boost::lexical_cast` throw — unreachable for the
/// committed well-formed asset); a blank OR malformed MR cell stores
/// 0.0 (the source's explicit `bad_lexical_cast` catch); notes are
/// carried in the asset but discarded, exactly as the parser ignores
/// them; no row is reordered and no value is invented. SMARTS pattern
/// compilation is NOT part of this owner yet (the `SmartsToMol` line
/// stays unmodeled until the Crippen pattern units).
///
/// Complexity review: the table parses ONCE into a process-global
/// vector (the source flyweight's parse-once semantics); the per-row
/// work is O(line length) with two `&'static str` slices borrowed from
/// the committed literal and no per-row allocation beyond the vector
/// slot itself; every later call is one O(1) atomic-load slice access,
/// matching the flyweight's amortized zero-cost lookup.
pub fn default_crippen_params() -> &'static [CrippenParamRow] {
    // RDKit source (Crippen.cpp:158-190):
    //   std::string inLine = RDKit::getLine(inStream);
    //   unsigned int idx = 0;
    //   while (!inStream.eof()) {
    //     if (inLine[0] != '#') {
    //       CrippenParams paramObj;
    //       paramObj.idx = idx++;
    //       tokenizer tokens(inLine, tabSep);
    //       tokenizer::iterator token = tokens.begin();
    //
    //       paramObj.label = *token;
    //       ++token;
    //       paramObj.smarts = *token;
    //       ++token;
    //       if (*token != "") {
    //         paramObj.logp = boost::lexical_cast<double>(*token);
    //       } else {
    //         paramObj.logp = 0.0;
    //       }
    //       ++token;
    //       if (*token != "") {
    //         try {
    //           paramObj.mr = boost::lexical_cast<double>(*token);
    //         } catch (boost::bad_lexical_cast &) {
    //           paramObj.mr = 0.0;
    //         }
    //       } else {
    //         paramObj.mr = 0.0;
    //       }
    //       paramObj.dp_pattern =
    //           boost::shared_ptr<const ROMol>(SmartsToMol(paramObj.smarts));
    //       d_params.push_back(paramObj);
    //     }
    //     inLine = RDKit::getLine(inStream);
    //   }
    // RDKit✔️✔️:   std::string inLine = RDKit::getLine(inStream);
    // RDKit✔️✔️:   unsigned int idx = 0;
    // RDKit✔️✔️:   while (!inStream.eof()) {
    // RDKit✔️✔️:     if (inLine[0] != '#') {
    // RDKit✔️✔️:       CrippenParams paramObj;
    // RDKit✔️✔️:       paramObj.idx = idx++;
    // RDKit✔️✔️:       tokenizer tokens(inLine, tabSep);
    // RDKit✔️✔️:       tokenizer::iterator token = tokens.begin();
    // RDKit✔️✔️:       paramObj.label = *token;
    // RDKit✔️✔️:       ++token;
    // RDKit✔️✔️:       paramObj.smarts = *token;
    // RDKit✔️✔️:       ++token;
    // RDKit✔️✔️:       if (*token != "") {
    // RDKit✔️✔️:         paramObj.logp = boost::lexical_cast<double>(*token);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         paramObj.logp = 0.0;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       ++token;
    // RDKit✔️✔️:       if (*token != "") {
    // RDKit✔️✔️:         try {
    // RDKit✔️✔️:           paramObj.mr = boost::lexical_cast<double>(*token);
    // RDKit✔️✔️:         } catch (boost::bad_lexical_cast &) {
    // RDKit✔️✔️:           paramObj.mr = 0.0;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         paramObj.mr = 0.0;
    // RDKit✔️✔️:       }
    // (The dp_pattern acquisition lines 184-185 are modeled by
    // `compile_crippen_patterns` below — the C03 pattern owner; this
    // row accessor keeps the raw `smarts` text.)
    // RDKit✔️✔️:       d_params.push_back(paramObj);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     inLine = RDKit::getLine(inStream);
    // RDKit✔️✔️:   }
    // Typed mapping: the header line above is the first `getLine`; the
    // istringstream loop maps to `lines()` over the committed literal;
    // `split('\t')` keeps empty tokens exactly like tabSep with
    // keep_empty_tokens; the malformed-logP abort maps to `expect`
    // (unreachable for this committed asset); the MR try/catch maps to
    // `unwrap_or(0.0)`.
    static ROWS: std::sync::OnceLock<Vec<CrippenParamRow>> = std::sync::OnceLock::new();
    ROWS.get_or_init(|| {
        let mut rows = Vec::new();
        let mut idx = 0u32;
        for (line_no, line) in DEFAULT_PARAM_DATA.lines().enumerate() {
            if line_no == 0 || line.starts_with('#') {
                continue;
            }
            let mut fields = line.split('\t');
            let label = fields.next().unwrap_or("");
            let smarts = fields.next().unwrap_or("");
            let logp = parse_crippen_logp(fields.next().unwrap_or(""))
                .expect("committed Crippen asset has valid logP cells");
            let mr = parse_crippen_mr(fields.next().unwrap_or(""));
            rows.push(CrippenParamRow {
                idx,
                label,
                smarts,
                logp,
                mr,
            });
            idx += 1;
        }
        rows
    })
}
