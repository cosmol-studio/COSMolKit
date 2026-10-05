//! Private Gemmi CID selection foundations (BIO-SEL packet).
//!
//! Source pins: third_party/gemmi/include/gemmi/select.hpp (Selection and
//! its nested predicates), third_party/gemmi/src/select.cpp (CID grammar
//! helpers), third_party/gemmi/include/gemmi/util.hpp (`is_in_list`),
//! third_party/gemmi/include/gemmi/iterator.hpp (FilterProxy laziness
//! contract). This module is private: nothing here is a public selection
//! API, and per the packet's decision queue no numeric parsing/formatting
//! owner, custom-flag storage or name-conversion owner is implied.

/// Gemmi `is_in_list` (util.hpp:227-237): comma-separated membership with
/// the source's exact length shortcut and tail-segment semantics.
use crate::hierarchy::EntityKind;
use crate::source_ids::AltLocLabel;

/// Private selection-copying child (BIO-COPY packet): detached
/// `empty_copy`/`add_matching_children` ports reusing this module's
/// canonical `Selection` and borrowed row predicates.
#[path = "selection_copy.rs"]
mod selection_copy;
pub use selection_copy::{
    BioSelectionBlocks, SelectionCopyCause as BioSelectionCopyCause,
    SelectionCopyError as BioSelectionCopyError, copy_selection_blocks,
};

pub(crate) fn is_in_list(name: &[u8], list: &[u8]) -> bool {
    // Gemmi✔️✔️: inline bool is_in_list(const std::string& name, const std::string& list,
    // Gemmi✔️✔️:                        char sep=',') {
    // Gemmi✔️✔️:   if (name.length() >= list.length())
    // Gemmi✔️✔️:     return name == list;
    // Gemmi✔️✔️:   for (size_t start=0, end=0; end != std::string::npos; start=end+1) {
    // Gemmi✔️✔️:     end = list.find(sep, start);
    // Gemmi✔️✔️:     if (list.compare(start, end - start, name) == 0)
    // Gemmi✔️✔️:       return true;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   return false;
    // Gemmi✔️✔️: }
    // Behavior: a name at least as long as the whole list matches only by
    // full equality (it cannot be a shorter member); otherwise each
    // comma-delimited segment — including the final tail and any empty
    // segments produced by consecutive/leading/trailing commas — is
    // compared literally (byte equality, no case folding, no trimming).
    // Complexity: O(list bytes) single pass with no allocation, matching
    // the source's in-place scanning.
    if name.len() >= list.len() {
        return name == list;
    }
    let mut start = 0;
    loop {
        let end = list[start..]
            .iter()
            .position(|byte| *byte == b',')
            .map_or_else(|| list.len(), |offset| start + offset);
        if &list[start..end] == name {
            return true;
        }
        if end == list.len() {
            return false;
        }
        start = end + 1;
    }
}

/// Gemmi `Selection::List` (select.hpp:22-39): all/inverted/comma list.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SelectionList {
    /// `*`: match everything.
    pub all: bool,
    /// `!`: invert the comma-list membership result.
    pub inverted: bool,
    /// Comma-separated member list (no escaping).
    pub list: String,
}

impl SelectionList {
    /// Gemmi `Selection::List::has` (select.hpp:33-38).
    pub(crate) fn has(&self, name: &str) -> bool {
        // Gemmi✔️✔️:     bool has(const std::string& name) const {
        // Gemmi✔️✔️:       if (all)
        // Gemmi✔️✔️:         return true;
        // Gemmi✔️✔️:       bool found = is_in_list(name, list);
        // Gemmi✔️✔️:       return inverted ? !found : found;
        // Gemmi✔️✔️:     }
        // Behavior: `all` short-circuits to true for every name (even an
        // empty one, even when a list is also stored); otherwise membership
        // follows is_in_list literally, optionally inverted.
        // Complexity: constant for `all`, else the is_in_list single scan.
        if self.all {
            return true;
        }
        let found = is_in_list(name.as_bytes(), self.list.as_bytes());
        if self.inverted { !found } else { found }
    }

    /// The source `List` member-initializer state (`all = true`).
    pub(crate) const fn default_all() -> Self {
        Self {
            all: true,
            inverted: false,
            list: String::new(),
        }
    }

    /// Byte-level form of `has` for source single-char fields (altloc) whose
    /// raw byte is not guaranteed to be UTF-8; identical semantics.
    pub(crate) fn has_bytes(&self, name: &[u8]) -> bool {
        if self.all {
            return true;
        }
        let found = is_in_list(name, self.list.as_bytes());
        if self.inverted { !found } else { found }
    }
}

/// Gemmi `Selection::FlagList` (select.hpp:41-51): a pattern over single
/// byte flags. BIO storage has no per-row custom flag field (decision
/// queue), so the predicate itself takes the explicit scalar flag byte.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct SelectionFlagList {
    /// Empty pattern matches everything; otherwise a byte pattern with an
    /// optional leading `!` inverter.
    pub(crate) pattern: String,
}

impl SelectionFlagList {
    /// Gemmi `Selection::FlagList::has` (select.hpp:43-50).
    pub(crate) fn has(&self, flag: u8) -> bool {
        // Gemmi✔️✔️:     bool has(char flag) const {
        // Gemmi✔️✔️:       if (pattern.empty())
        // Gemmi✔️✔️:         return true;
        // Gemmi✔️✔️:       bool invert = (pattern[0] == '!');
        // Gemmi✔️✔️:       bool found = (pattern.find(flag, invert ? 1 : 0) != std::string::npos);
        // Gemmi✔️✔️:       return invert ? !found : found;
        // Gemmi✔️✔️:     }
        // Behavior: an empty pattern matches every flag byte; a pattern
        // starting with '!' inverts and never matches the '!' itself (the
        // search starts after it); otherwise plain byte containment of the
        // flag anywhere in the pattern (a '!' in a non-leading position is
        // an ordinary member).
        // Complexity: O(pattern bytes) single scan, no allocation.
        if self.pattern.is_empty() {
            return true;
        }
        let bytes = self.pattern.as_bytes();
        let invert = bytes[0] == b'!';
        let found = bytes[invert as usize..].contains(&flag);
        if invert { !found } else { found }
    }
}

/// Gemmi `Selection::SequenceId` (select.hpp:52-68): sequence-number
/// selector with INT_MIN/INT_MAX emptiness sentinels and an optional
/// insertion code where `b'*'` is the wildcard.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SelectionSequenceId {
    pub seqnum: i32,
    /// Blank (`b' '`), a letter, or the `b'*'` wildcard.
    pub icode: u8,
}

impl SelectionSequenceId {
    /// Gemmi `Selection::SequenceId::empty` (select.hpp:57-59).
    pub(crate) fn empty(&self) -> bool {
        // Gemmi✔️✔️:     bool empty() const {
        // Gemmi✔️✔️:       return seqnum == INT_MIN || seqnum == INT_MAX;
        // Gemmi✔️✔️:     }
        // Behavior: both range endpoints are "unset" markers.
        // Complexity: two integer comparisons.
        self.seqnum == i32::MIN || self.seqnum == i32::MAX
    }

    /// Gemmi `Selection::SequenceId::compare` (select.hpp:63-68) against a
    /// source `SeqId` carried as its raw number (INT_MIN when unset) and
    /// icode byte (blank when unset).
    pub(crate) fn compare(&self, seqid_num: i32, seqid_icode: u8) -> i32 {
        // Gemmi✔️✔️:     int compare(const SeqId& seqid) const {
        // Gemmi✔️✔️:       if (seqnum != *seqid.num)
        // Gemmi✔️✔️:         return seqnum < *seqid.num ? -1 : 1;
        // Gemmi✔️✔️:       if (icode != '*' && icode != seqid.icode)
        // Gemmi✔️✔️:         return icode < seqid.icode ? -1 : 1;
        // Gemmi✔️✔️:       return 0;
        // Gemmi✔️✔️:     }
        // Behavior: primary ordering by sequence number with strict
        // inequalities; a wildcard selector icode equals every source
        // icode; blank and letter icodes compare by byte order (blank
        // sorts before letters); equal numbers with equal-or-wildcard
        // icode compare 0.
        // Complexity: at most two integer/byte comparisons.
        if self.seqnum != seqid_num {
            return if self.seqnum < seqid_num { -1 } else { 1 };
        }
        if self.icode != b'*' && self.icode != seqid_icode {
            return if self.icode < seqid_icode { -1 } else { 1 };
        }
        0
    }
}

/// Gemmi `Selection::AtomInequality` (select.hpp:70-82) minus its numeric
/// payload parser (decision queue: single approved numeric owner).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SelectionAtomInequality {
    /// `b'q'` occupancy, `b'b'` B-factor; any other byte selects neither.
    pub property: u8,
    /// Negative: `<`, zero: `=`, positive: `>`.
    pub relation: i32,
    /// Comparison threshold (f64 in the source; Gemmi's Atom stores float,
    /// so the compared carrier is the f32-narrowed atom value).
    pub value: f64,
}

impl SelectionAtomInequality {
    /// Gemmi `Selection::AtomInequality::matches` (select.hpp:74-82) with
    /// the atom's occupancy and B-factor supplied as the source's f32
    /// carriers (BIO stores f64; narrow once, like Gemmi's float members).
    pub(crate) fn matches(&self, atom_occ: f32, atom_b_iso: f32) -> bool {
        // Gemmi✔️✔️:     bool matches(const Atom& a) const {
        // Gemmi✔️✔️:       double atom_value = 0.;
        // Gemmi✔️✔️:       if (property == 'q')
        // Gemmi✔️✔️:         atom_value = a.occ;
        // Gemmi✔️✔️:       else if (property == 'b')
        // Gemmi✔️✔️:         atom_value = a.b_iso;
        // Gemmi✔️✔️:       if (relation < 0)
        // Gemmi✔️✔️:         return atom_value < value;
        // Gemmi✔️✔️:       if (relation > 0)
        // Gemmi✔️✔️:         return atom_value > value;
        // Gemmi✔️✔️:       return atom_value == value;
        // Gemmi✔️✔️:     }
        // Behavior: an unknown property leaves atom_value at the initial
        // 0.0 — so `<` holds against positive thresholds, `>` against
        // negative ones, and `=` ONLY against exactly 0.0 (it does not
        // hold against positive thresholds); comparisons are IEEE double
        // semantics after the
        // f64 promotion of the f32 carriers — NaN never satisfies any
        // strict or equality relation; -0.0 == 0.0 holds; the payload
        // itself is not parsed here (numeric-owner decision pending).
        // Complexity: one property dispatch and one comparison; no
        // allocation.
        let mut atom_value = 0.0f64;
        if self.property == b'q' {
            atom_value = f64::from(atom_occ);
        } else if self.property == b'b' {
            atom_value = f64::from(atom_b_iso);
        }
        if self.relation < 0 {
            atom_value < self.value
        } else if self.relation > 0 {
            atom_value > self.value
        } else {
            atom_value == self.value
        }
    }
}

// (BIO-CID C02) `wrong_syntax`/`SelectionSyntaxError` (select.cpp:14-26)
// were relocated verbatim to `cosmolkit-io/src/bio_cid.rs`, the single IO
// owner of the CID lexical closure. BIO keeps no copy.

// (BIO-CID C02) `OmittedCidFields` and `determine_omitted_cid_fields`
// (select.cpp:28-41) were relocated verbatim to
// `cosmolkit-io/src/bio_cid.rs`.

/// Parser-side CID field list carrying the exact source byte payload:
/// the source stores `cid.substr(pos, end - pos)` — raw bytes with no
/// UTF-8 requirement (select.cpp:48). The matcher-side `SelectionList`
/// (UTF-8 `String`) is a separate carrier; converting between them is a
/// caller decision that must not happen silently inside this helper.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
/// Detached CID list vocabulary: a `[!]/[*]/list` field's exact bytes
/// after `make_cid_list`. Pub in `cosmolkit-bio` for the IO CID decoder
/// (BIO-CID C01 frozen boundary); the root crate does not re-export it.
pub struct SelectionCidList {
    /// `*`: match everything.
    pub all: bool,
    /// `!`: invert the comma-list membership result.
    pub inverted: bool,
    /// Exact member bytes as produced by the source `substr`, including
    /// multibyte UTF-8 names and every empty comma member.
    pub list: Vec<u8>,
}

// (BIO-CID C02) `make_cid_list` (select.cpp:43-57) was relocated
// verbatim to `cosmolkit-io/src/bio_cid.rs`.

// (BIO-CID C02) `SelectionSeqidRangeError` and `parse_cid_seqid`
// (select.cpp:95-115) were relocated verbatim to
// `cosmolkit-io/src/bio_cid.rs`.

/// Gemmi `El::END` (elem.hpp:14-26): one past the last ordinal (D=119).
/// Pub in `cosmolkit-bio`: the minimal immutable element-lookup boundary
/// frozen by BIO-CID C01 for the IO CID decoder.
pub const GEMMI_EL_END: usize = 120;

/// Gemmi `element_uppercase_name` table (elem.hpp:296-313), El-ordinal
/// order. Pub in `cosmolkit-bio` as canonical element vocabulary: the
/// relocated IO CID owner imports it for its moved pinned-name
/// regressions and any future name resolution; IO must not copy it.
/// The root crate does not re-export it.
pub const GEMMI_ELEMENT_NAMES: [&str; GEMMI_EL_END] = [
    "X", "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE", "NA", "MG", "AL", "SI", "P", "S",
    "CL", "AR", "K", "CA", "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE",
    "AS", "SE", "BR", "KR", "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD",
    "IN", "SN", "SB", "TE", "I", "XE", "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD",
    "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG",
    "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH", "PA", "U", "NP", "PU", "AM", "CM",
    "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG", "CN",
    "NH", "FL", "MC", "LV", "TS", "OG", "D",
];

/// Gemmi `is_metal_value` table (elem.hpp:98-130), first `El::END` (120)
/// entries — the END slot at index 120 is never read by the source loops.
///
/// Canonical classifier boundary (BIO-CID C01/FIX): the IO CID decoder's
/// relocated `parse_cid_elements` MUST import this exact owner for
/// `metals`/`nonmetals` group expansion (`cosmolkit_bio::GEMMI_IS_METAL`);
/// IO must not copy the table or reimplement a classifier. The root crate
/// must not re-export it. Values and source provenance are unchanged from
/// the original private table (verified byte-for-byte against elem.hpp).
pub const GEMMI_IS_METAL: [bool; GEMMI_EL_END] = GEMMI_IS_METAL_SRC;

const GEMMI_IS_METAL_SRC: [bool; GEMMI_EL_END] = [
    false, false, false, true, true, false, false, false, false, false, false, true, true, true,
    false, false, false, false, false, true, true, true, true, true, true, true, true, true, true,
    true, true, true, true, false, false, false, false, true, true, true, true, true, true, true,
    true, true, true, true, true, true, true, true, false, false, false, true, true, true, true,
    true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    true, true, true, true, true, true, true, true, true, true, false, false, true, true, true,
    true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    true, true, true, true, true, true, true, true, true, true, true, false, false, false,
];

/// Gemmi `impl::find_single_letter_element` (elem.hpp:320-336).
const fn find_single_letter_element(c: u8) -> u8 {
    match c {
        b'H' => 1,
        b'B' => 5,
        b'C' => 6,
        b'N' => 7,
        b'O' => 8,
        b'F' => 9,
        b'P' => 15,
        b'S' => 16,
        b'K' => 19,
        b'V' => 23,
        b'Y' => 39,
        b'I' => 53,
        b'W' => 74,
        b'U' => 92,
        b'D' => 119,
        _ => 0,
    }
}

/// Gemmi `element_name(El)` (elem.hpp:269-294) — the canonical DISPLAY
/// name table ("He", "Li", "Cl", ...), distinct from the parser's
/// uppercase `element_uppercase_name` vocabulary in
/// [`GEMMI_ELEMENT_NAMES`]. Ordinals follow the Gemmi El enum
/// (X=0..Og=118, D=119); `GEMMI_EL_END` (120) returns the source's
/// empty END string and any higher ordinal returns `None` (a safe Rust
/// boundary, not a claim about source-undefined behavior). Narrow owner
/// exported for the IO CID serializer; NOT re-exported through the root
/// crate.
pub fn gemmi_element_name(ordinal: usize) -> Option<&'static str> {
    // Gemmi✔️✔️: typedef const char elname_t[3];
    // Gemmi✔️✔️: inline const char* element_name(El el) {
    // Gemmi✔️✔️:   static constexpr elname_t names[] = {
    // Gemmi✔️✔️:     "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O", "F", "Ne",
    // Gemmi✔️✔️:     "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",
    // Gemmi✔️✔️:     "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co",
    // Gemmi✔️✔️:     "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    // Gemmi✔️✔️:     "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
    // Gemmi✔️✔️:     "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    // Gemmi✔️✔️:     "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu",
    // Gemmi✔️✔️:     "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    // Gemmi✔️✔️:     "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    // Gemmi✔️✔️:     "Tl", "Pb", "Bi", "Po", "At", "Rn",
    // Gemmi✔️✔️:     "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am",
    // Gemmi✔️✔️:     "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    // Gemmi✔️✔️:     "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    // Gemmi✔️✔️:     "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
    // Gemmi✔️✔️:     "D", ""
    // Gemmi✔️✔️:   };
    // Gemmi✔️✔️:   static_assert(static_cast<int>(El::Og) == 118, "Hmm");
    // Gemmi✔️✔️:   static_assert(names[118][0] == 'O', "Hmm");
    // Gemmi✔️✔️:   static_assert(sizeof(names) / sizeof(names[0]) ==
    // Gemmi✔️✔️:                 static_cast<int>(El::END) + 1, "Hmm");
    // Gemmi✔️✔️:   return names[static_cast<int>(el)];
    // Gemmi✔️✔️: }
    // (elem.hpp:267-294, verbatim.)
    //
    // Behavior review: exact display spellings in exact Gemmi ordinal
    // order — 121 entries covering X..Og (0..=118), D (119) and the END
    // placeholder "" (120), which select.cpp:314's element_name call
    // serializes for CID element masks ("He"/"Li"/"Cl", never the
    // uppercase parser vocabulary). Bounds beyond the table return None
    // instead of indexing undefined source memory. Complexity review:
    // one constant table lookup, no scan.
    const NAMES: [&str; 121] = [
        "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P",
        "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
        "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re",
        "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db",
        "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", "D", "",
    ];
    NAMES.get(ordinal).copied()
}

/// Gemmi `find_element` (elem.hpp:343-360) over a one/two-byte symbol;
/// returns the Gemmi El ordinal (0 = X).
pub fn gemmi_find_element(first: u8, second: u8) -> u8 {
    // Gemmi✔️✔️: inline El find_element(const char* symbol) {
    // Gemmi✔️✔️:   if (symbol == nullptr || symbol[0] == '\0')
    // Gemmi✔️✔️:     return El::X;
    // Gemmi✔️✔️:   char first = symbol[0] & ~0x20;  // lower -> upper, space -> NUL
    // Gemmi✔️✔️:   char second = symbol[1] & ~0x20;
    // Gemmi✔️✔️:   if (first == '\0')
    // Gemmi✔️✔️:     return impl::find_single_letter_element(second);
    // Gemmi✔️✔️:   if (second < 14)
    // Gemmi✔️✔️:     return impl::find_single_letter_element(first);
    // Gemmi✔️✔️:   elname_t* names = &element_uppercase_name(El::X);
    // Gemmi✔️✔️:   for (int i = 0; i != 120; ++i) {
    // Gemmi✔️✔️:     if (names[i][0] == first && names[i][1] == second)
    // Gemmi✔️✔️:       return static_cast<El>(i);
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   return El::X;
    // Gemmi✔️✔️: }
    // Behavior: `& ~0x20` folds lowercase to uppercase and maps space to
    // NUL; a NUL first byte delegates to the single-letter table on the
    // second byte; a second byte below 14 (NUL, controls, punctuation up
    // to '\r') takes the single-letter table on the first byte; otherwise
    // the two-letter uppercase name table decides, defaulting to X.
    // Complexity: at most 120 comparisons of two bytes each, no allocation.
    let first = first & !0x20;
    let second = second & !0x20;
    if first == 0 {
        return find_single_letter_element(second);
    }
    if second < 14 {
        return find_single_letter_element(first);
    }
    for (index, name) in GEMMI_ELEMENT_NAMES.iter().enumerate() {
        let name = name.as_bytes();
        if name[0] == first && name.get(1).copied().unwrap_or(0) == second {
            return index as u8;
        }
    }
    0
}

/// A Gemmi-ordinal element mask (select.hpp `std::vector<char>`).
pub type GemmiElementMask = [bool; GEMMI_EL_END];

// (BIO-CID C02) `parse_cid_elements` (select.cpp:59-93) and
// `has_inequality` (select.cpp:143-148) were relocated verbatim to
// `cosmolkit-io/src/bio_cid.rs`; the element ordinals, name lookup and
// metals classifier remain here and are imported there through the
// C01/FIX boundary.

/// Gemmi `EntityType` ordinal (metadata.hpp:206-214) used to index
/// `et_flags`: CK `EntityKind` declares the same five kinds in the same
/// order (Unknown, Polymer, NonPolymer, Branched, Water); the source's
/// sixth array slot is unreachable from both type systems.
pub(crate) const fn entity_kind_index(kind: EntityKind) -> usize {
    match kind {
        EntityKind::Unknown => 0,
        EntityKind::Polymer => 1,
        EntityKind::NonPolymer => 2,
        EntityKind::Branched => 3,
        EntityKind::Water => 4,
    }
}

/// Source-shaped borrowed Model input for the private predicates below.
/// Only fields the Gemmi predicates read are carried; converting detached
/// `BioModelRow` storage into these views (including any policy for a
/// missing source model number) is a decision-queue row-adapter item, not
/// part of this packet.
/// CK missing-representation boundary for model-row selection (BIO-ROWS
/// R03): the stored row has no source model number, so no source-backed
/// predicate value exists. Not a Gemmi behavior claim — Gemmi's
/// `Model::num` is always reader-assigned.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BioRowModelError {
    /// `BioModelRow::source_model_number()` is `None`.
    MissingModelNumber,
}

impl std::fmt::Display for BioRowModelError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BioRowModelError::MissingModelNumber => {
                write!(f, "model row has no source model number")
            }
        }
    }
}

impl std::error::Error for BioRowModelError {}

/// Test-only gate-call counters (BIO-ROWS-CURSOR): compiled solely under
/// `cfg(test)` to prove each parent gate runs once per visited parent row
/// in the cursor regressions. Never present in production builds.
#[cfg(test)]
pub(crate) mod gate_counters {
    use std::cell::Cell;

    // Each synchronous cursor test measures its own calls. Shared counters
    // include gates run by unrelated tests on Rust's parallel test threads.
    std::thread_local! {
        pub(crate) static MODEL: Cell<usize> = const { Cell::new(0) };
        pub(crate) static CHAIN: Cell<usize> = const { Cell::new(0) };
        pub(crate) static RESIDUE: Cell<usize> = const { Cell::new(0) };
    }

    pub(crate) fn reset() {
        MODEL.with(|count| count.set(0));
        CHAIN.with(|count| count.set(0));
        RESIDUE.with(|count| count.set(0));
    }
}

/// Model-row matching on an explicit present source model number
/// (BIO-ROWS R03).
///
/// Extracts the row's source number and applies the pinned model
/// predicate. `None` rows return the typed [`BioRowModelError`] instead of
/// guessing Gemmi's default `num = 0`, the row index, or any other value.
pub(crate) fn bio_model_row_matches(
    selection: &Selection,
    row: &crate::hierarchy::BioModelRow,
) -> Result<bool, BioRowModelError> {
    // Gemmi✔️✔️: int num = 0;                       // Model default (model.hpp)
    // Gemmi✔️✔️: int num = read_int(line+6, 8);      // PDB MODEL record
    // Gemmi✔️✔️: model = &st.find_or_add_model(num);
    // Gemmi✔️✔️: int num = (int) st.models.size() + 1; // implicit single model
    // Gemmi✔️✔️: {"pdbx_PDB_model_num", ...},        // mmCIF atom_site tag
    // Behavior: every pinned Gemmi reader ASSIGNS `Model::num` before any
    // selection can run (PDB MODEL serial via read_int(line+6,8), implicit
    // models via size+1, mmCIF via pdbx_PDB_model_num), so `matches(const
    // Model&)` always sees a present number. CK's `BioModelRow` stores
    // `Option<i32>`: `Some` carries the exact reader-assigned number
    // (PDB serial or size+1; CK's mmCIF grouping records null tokens as
    // 0, the reader's own recorded assignment); `None` has no
    // source-backed value at all, so matching fails closed with the typed
    // error rather than substituting Gemmi's field default, the row ID,
    // or 0 — a CK missing-representation boundary, not a Gemmi behavior
    // claim. The predicate itself reuses the ported `matches_model`
    // (mdl == 0 wildcard || exact int equality) with an empty chain
    // slice, which the model predicate never reads.
    // Complexity: O(1); one wrapped struct over the existing predicate.
    #[cfg(test)]
    gate_counters::MODEL.with(|count| count.set(count.get() + 1));
    let num = row
        .source_model_number()
        .ok_or(BioRowModelError::MissingModelNumber)?;
    Ok(selection.matches_model(&SourceModel { num, chains: &[] }))
}

/// CK missing-representation boundary for chain-row selection (BIO-ROWS
/// R04): the stored row carries neither an auth nor a label source chain
/// identity, so no reader-populated canonical name exists. Not a Gemmi
/// behavior claim — Gemmi's `Chain::name` is always reader-assigned.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BioRowChainError {
    /// `ChainSourceIds` has neither `auth_chain_id` nor `label_asym_id`.
    MissingCanonicalChainName,
}

impl std::fmt::Display for BioRowChainError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BioRowChainError::MissingCanonicalChainName => {
                write!(f, "chain row has no source chain identity")
            }
        }
    }
}

impl std::error::Error for BioRowChainError {}

/// Chain-row matching through the canonical reader-populated source name
/// (BIO-ROWS R04).
///
/// Precedence is derived from the pinned reader paths — never a fallback
/// guess — and a missing canonical identity returns the typed error
/// rather than a row ID, empty string, or other substitute.
pub(crate) fn bio_chain_row_matches(
    selection: &Selection,
    row: &crate::hierarchy::BioChainRow,
    input_format: crate::hierarchy::BioCoordinateFormat,
) -> Result<bool, BioRowChainError> {
    // Gemmi✔️✔️: std::string chain_name = read_string(line+20, 2);  // PDB reader
    // Gemmi✔️✔️: RowAccess asym_id(atom_table, kAuthAsymId, kLabelAsymId); // mmCIF
    // Gemmi✔️✔️: if (!chain || cif::as_string(asym_id.get(gap)) != chain->name) {
    // Gemmi✔️✔️:   model->chains.emplace_back(cif::as_string(asym_id.get(gap)));
    // Gemmi✔️✔️:   bool matches(const Chain& chain) const {
    // Gemmi✔️✔️:     return chain_ids.has(chain.name);
    // Gemmi✔️✔️:   }
    // Behavior: the mmCIF reader resolves the chain column through a
    // two-slot RowAccess with AUTH first and LABEL as the fallback
    // (kAuthAsymId before kLabelAsymId), and the resulting string becomes
    // Chain::name verbatim; the PDB reader takes the read_string view of
    // the raw two-column chain field (a blank PDB field is the EMPTY
    // canonical name, not a missing one). CK mirrors this precedence in
    // `ChainSourceIds`: `auth_chain_id` (PDB column / auth_asym_id) is
    // preferred; `label_asym_id` is used only when auth is absent. Pdb
    // provenance applies the exact read_string view to the stored raw
    // column; non-PDB provenance preserves the decoded auth value
    // verbatim, exactly as as_string does. With both identities absent
    // there is no reader-populated canonical name: the typed error is
    // returned. The derived name feeds the ported exact-membership
    // predicate with an empty residue slice (never read). The canonical
    // name feeds the ported exact-membership predicate with an empty
    // residue slice (never read). BIO-CID C23: the narrow borrowed
    // `auth_chain_id_ref` accessor resolves the name as a borrow — no
    // per-chain String allocation on any path.
    // Complexity: O(len(name)) membership probe over borrowed data;
    // auth, PDB-view, and label paths all borrow directly.
    #[cfg(test)]
    gate_counters::CHAIN.with(|count| count.set(count.get() + 1));
    let source = row.source();
    let name: &str = if let Some(auth) = source.auth_chain_id_ref() {
        if input_format == crate::hierarchy::BioCoordinateFormat::Pdb {
            crate::hierarchy::pdb_read_string_view(auth.as_str())
        } else {
            auth.as_str()
        }
    } else {
        source
            .label_asym_id()
            .ok_or(BioRowChainError::MissingCanonicalChainName)?
    };
    Ok(selection.matches_chain(&SourceChain {
        name,
        residues: &[],
    }))
}

/// Chain-row matching through the canonical source name (BIO-ROWS R04).

/// Residue source seqid adapter (BIO-ROWS R05): maps the row's stored
/// `Option<PdbSeqId>` onto the pinned `SeqId` absence sentinels —
/// `INT_MIN` number and blank (`b' '`) icode — preserving negative,
/// zero and positive numbers and letter insertion codes byte-exactly.
#[must_use]
pub(crate) fn bio_residue_seqid(row: &crate::hierarchy::BioResidueRow) -> (i32, u8) {
    // Gemmi✔️✔️: using OptionalNum = OptionalInt<INT_MIN>;
    // Gemmi✔️✔️: OptionalNum num;   // sequence number
    // Gemmi✔️✔️: char icode = ' ';  // insertion code
    // Gemmi✔️✔️: SeqId(int num_, char icode_) { num = num_; icode = icode_; }
    // Behavior: a missing seq_id maps to the source's own absence
    // sentinels (num = INT_MIN, icode = ' '), exactly as a default
    // SeqId reads; a present seq_id keeps its stored number (including
    // negative and zero values) and its insertion-code byte, with
    // `ins_code: None` meaning the blank icode. No defaulting to zero,
    // row index, or any other guessed value, and no label_seq_id
    // substitution (label numbering is a separate source identity).
    // Complexity: O(1), no allocation.
    match row.source().seq_id() {
        Some(seqid) => (seqid.seq_num(), seqid.ins_code().unwrap_or(b' ')),
        None => (i32::MIN, b' '),
    }
}

/// Direct borrowed residue-row predicate (BIO-ROWS R06).
///
/// Assembles the five source conjuncts from the row via the R02 logical
/// name view and the R05 seqid adapter, plus the source entity kind. The
/// custom flag remains an explicit caller-supplied byte: BIO rows store no
/// custom residue flag, and substituting het_flag/calc_flag is forbidden
/// (decision queue) — the caller passes the flag their context defines,
/// defaulting to the Gemmi unset byte.
pub(crate) fn bio_residue_row_matches(
    selection: &Selection,
    row: &crate::hierarchy::BioResidueRow,
    input_format: crate::hierarchy::BioCoordinateFormat,
    flag: u8,
) -> bool {
    // Gemmi✔️✔️:   bool matches(const Residue& res) const {
    // Gemmi✔️✔️:     return (entity_types.all || et_flags[(int)res.entity_type]) &&
    // Gemmi✔️✔️:            residue_names.has(res.name) &&
    // Gemmi✔️✔️:            from_seqid.compare(res.seqid) <= 0 &&
    // Gemmi✔️✔️:            to_seqid.compare(res.seqid) >= 0 &&
    // Gemmi✔️✔️:            residue_flags.has(res.flag);
    // Gemmi✔️✔️:   }
    // Behavior: the row's name goes through the R02 logical view (PDB
    // read_string trim; CIF-family verbatim), the seqid through the R05
    // adapter (INT_MIN/blank sentinels for a missing seq_id), and the
    // entity kind through the row's stored source entity_kind (equal to
    // Gemmi EntityType ordinals). The flag byte is the caller's explicit
    // input — never derived from het_flag or any other row field — and
    // feeds the exact five-conjunct ported predicate with an empty atom
    // slice (never read by the residue conjuncts).
    // Complexity: as the ported predicate; the views are borrowed, O(1).
    #[cfg(test)]
    gate_counters::RESIDUE.with(|count| count.set(count.get() + 1));
    let (seqid_num, icode) = bio_residue_seqid(row);
    let row_name = row.name();
    // Borrowed logical view held for the synchronous predicate call —
    // no per-row String allocation (BIO-C24-CLOSE removed the copy).
    let name = crate::hierarchy::residue_name_logical_view(&row_name, input_format);
    selection.matches_residue(&SourceResidue {
        entity_type: row.entity_kind(),
        name,
        seqid_num,
        icode,
        flag,
        atoms: &[],
    })
}

/// Atom element identity adapter (BIO-ROWS R07): maps the row's typed
/// `(Element, isotope)` state onto the Gemmi element ordinal used by the
/// selection mask, preserving deuterium identity.
#[must_use]
pub(crate) fn bio_atom_element_ordinal(row: &crate::hierarchy::BioAtomRow) -> u8 {
    // Gemmi✔️✔️: inline bool is_hydrogen(El el) { return el == El::H || el == El::D; }
    // Gemmi✔️✔️: if (alpha_up(name[0]) == 'D')    // PDB deuterium inference
    // Gemmi✔️✔️:   return El::D;
    // Behavior: Gemmi's element ordinals are the atomic numbers for
    // X..Og (0..=118) with D at 119. The CK PDB reader projects source
    // El::D to the approved H + isotope_mass_number Some(2)
    // representation, so H + Some(2) maps back to ordinal 119 exactly
    // (the reader's recorded projection, not a guess): deuterium keeps
    // its distinct selection identity and is not equated to plain H.
    // Every other (element, isotope) pair keeps the element ordinal —
    // isotopes of other elements are NOT distinct Gemmi Els — and D's
    // 119 is never conflated with an atomic number (no element 119
    // exists; the mask table ends at 119 exactly for D).
    // Complexity: O(1), no allocation.
    let element = row.element();
    if element.atomic_number() == 1 && row.isotope_mass_number() == Some(2) {
        119
    } else {
        element.atomic_number()
    }
}

/// Direct borrowed atom-row predicate (BIO-ROWS R08).
///
/// Feeds the ported atom conjuncts from the row via the R01 logical name
/// view, the R07 element-ordinal adapter (preserving deuterium 119), the
/// exact stored altloc byte, source `float` occupancy/B-factor values,
/// and the caller's explicit flag byte.
pub(crate) fn bio_atom_row_matches(
    selection: &Selection,
    row: &crate::hierarchy::BioAtomRow,
    input_format: crate::hierarchy::BioCoordinateFormat,
    flag: u8,
) -> bool {
    // Gemmi✔️✔️:   bool matches(const Atom& a) const {
    // Gemmi✔️✔️:     return atom_names.has(a.name) &&
    // Gemmi✔️✔️:            (elements.empty() || elements[a.element.ordinal()]) &&
    // Gemmi✔️✔️:            (altlocs.all || altlocs.has(std::string(a.altloc ? 1 : 0, a.altloc))) &&
    // Gemmi✔️✔️:            atom_flags.has(a.flag) &&
    // Gemmi✔️✔️:            std::all_of(atom_inequalities.begin(), atom_inequalities.end(),
    // Gemmi✔️✔️:                        [&](const AtomInequality& i) { return i.matches(a); });
    // Gemmi✔️✔️:   }
    // Behavior: the name goes through the R01 logical view (PDB
    // read_string trim; CIF-family verbatim), the element through the R07
    // ordinal adapter (H+Some(2) keeps the distinct D ordinal 119), and
    // the altloc is the exact stored label byte with `None` as the
    // source's 0 sentinel (empty one-byte name). Occupancy and B-factor
    // are the source `float` (f32) semantics: Gemmi stores `float occ`
    // and `float b_iso`, the CK row widens reader-parsed f32 values into
    // f64, and this adapter narrows back — bit-exact for reader-sourced
    // values, documented as the source-width view (inequalities compare
    // doubles after the f32 read in the source, matching this shape).
    // The flag byte is the caller's explicit input, never derived from
    // calc_flag or any other row field (decision queue).
    // Complexity: as the ported predicate; views borrowed, O(1).
    let row_name = row.name();
    // Borrowed logical view held for the synchronous predicate call —
    // no per-row String allocation (BIO-C24-CLOSE removed the copy).
    let name = crate::hierarchy::atom_name_logical_view(&row_name, input_format);
    selection.matches_atom(&SourceAtom {
        name,
        element: bio_atom_element_ordinal(row),
        altloc: row.altloc().map_or(0, AltLocLabel::value),
        flag,
        occ: row.occupancy() as f32,
        b_iso: row.b_iso() as f32,
    })
}

/// Typed original row IDs of a first-match hit (BIO-ROWS R09): the
/// chain/residue/atom table indices of the matching rows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct BioRowFirst {
    pub(crate) chain: crate::hierarchy::BioChainId,
    pub(crate) residue: crate::hierarchy::BioResidueId,
    pub(crate) atom: crate::hierarchy::BioAtomId,
}

/// Traversal-level error for row-based first-match (BIO-ROWS R09): the
/// row tables carry a missing source representation the source hierarchy
/// always has.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BioRowTraverseError {
    Model(BioRowModelError),
    Chain(BioRowChainError),
}

impl std::fmt::Display for BioRowTraverseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BioRowTraverseError::Model(e) => write!(f, "{e}"),
            BioRowTraverseError::Chain(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for BioRowTraverseError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        // Structured leaf retention (BIO-C24-CLOSE): the Model/Chain leaf
        // is the chain terminus (its own source is None).
        match self {
            Self::Model(cause) => Some(cause as &(dyn std::error::Error + 'static)),
            Self::Chain(cause) => Some(cause as &(dyn std::error::Error + 'static)),
        }
    }
}

/// Flat-hierarchy first-match traversal over actual BIO row spans
/// (BIO-ROWS R09), mirroring Gemmi `Selection::first_in_model` without
/// constructing `SourceModel` trees or copying rows.
///
/// Residue and atom custom flags are explicit caller-supplied bytes (the
/// BIO rows store no custom flags; het/calc substitution is forbidden).
pub(crate) fn bio_first_in_model_rows(
    selection: &Selection,
    model: &crate::hierarchy::BioModelRow,
    chains: &[crate::hierarchy::BioChainRow],
    residues: &[crate::hierarchy::BioResidueRow],
    atoms: &[crate::hierarchy::BioAtomRow],
    input_format: crate::hierarchy::BioCoordinateFormat,
    residue_flag: u8,
    atom_flag: u8,
) -> Result<Option<BioRowFirst>, BioRowTraverseError> {
    // Gemmi✔️✔️:   CRA first_in_model(Model& model) const {
    // Gemmi✔️✔️:     if (matches(model))
    // Gemmi✔️✔️:       for (Chain& chain : model.chains) {
    // Gemmi✔️✔️:         if (matches(chain))
    // Gemmi✔️✔️:           for (Residue& res : chain.residues) {
    // Gemmi✔️✔️:             if (matches(res))
    // Gemmi✔️✔️:               for (Atom& atom : res.atoms) {
    // Gemmi✔️✔️:                 if (matches(atom))
    // Gemmi✔️✔️:                   return {&chain, &res, &atom};
    // Gemmi✔️✔️:               }
    // Gemmi✔️✔️:           }
    // Gemmi✔️✔️:       }
    // Gemmi✔️✔️:     return {nullptr, nullptr, nullptr};
    // Gemmi✔️:   }
    // Behavior: identical gate order and short-circuiting over the flat
    // row tables: the model gate first (rejected model -> None), then the
    // model's chain_span slice, each chain's residue_span slice, and each
    // residue's atom_span slice, each row gated by the R03/R04/R06/R08
    // row predicates (which apply the logical-name views, seqid sentinels
    // and element ordinals). A result is returned ONLY when an atom
    // matches; a matching residue with no matching atom continues to the
    // next residue. The returned IDs are the ORIGINAL row table indices
    // (enumerate positions, noncontiguous source identities preserved);
    // no Source* trees, no row copies, no intermediate vectors.
    // Complexity: worst case one scan over the model's rows with early
    // exit; the per-row derived views are O(1).
    if !bio_model_row_matches(selection, model).map_err(BioRowTraverseError::Model)? {
        return Ok(None);
    }
    let chain_span = model.chain_span();
    for (chain_offset, chain) in chains
        [chain_span.start() as usize..(chain_span.start() + chain_span.len()) as usize]
        .iter()
        .enumerate()
    {
        if !bio_chain_row_matches(selection, chain, input_format)
            .map_err(BioRowTraverseError::Chain)?
        {
            continue;
        }
        let chain_id = crate::hierarchy::BioChainId::new(chain_span.start() + chain_offset as u32);
        let residue_span = chain.residue_span();
        for (residue_offset, residue) in residues
            [residue_span.start() as usize..(residue_span.start() + residue_span.len()) as usize]
            .iter()
            .enumerate()
        {
            if !bio_residue_row_matches(selection, residue, input_format, residue_flag) {
                continue;
            }
            let residue_id =
                crate::hierarchy::BioResidueId::new(residue_span.start() + residue_offset as u32);
            let atom_span = residue.atom_span();
            for (atom_offset, atom) in atoms
                [atom_span.start() as usize..(atom_span.start() + atom_span.len()) as usize]
                .iter()
                .enumerate()
            {
                if bio_atom_row_matches(selection, atom, input_format, atom_flag) {
                    return Ok(Some(BioRowFirst {
                        chain: chain_id,
                        residue: residue_id,
                        atom: crate::hierarchy::BioAtomId::new(
                            atom_span.start() + atom_offset as u32,
                        ),
                    }));
                }
            }
        }
    }
    Ok(None)
}

/// Typed original row IDs of a whole-structure first-match (BIO-ROWS
/// R10): the model table index plus the R09 triple.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct BioRowFirstHit {
    pub(crate) model: crate::hierarchy::BioModelId,
    pub(crate) cra: BioRowFirst,
}

/// Ordered whole-structure first-match over existing row tables
/// (BIO-ROWS R10), mirroring Gemmi `Selection::first` in stored source
/// order with explicit flag-provider inputs.
pub(crate) fn bio_first_rows(
    selection: &Selection,
    models: &[crate::hierarchy::BioModelRow],
    chains: &[crate::hierarchy::BioChainRow],
    residues: &[crate::hierarchy::BioResidueRow],
    atoms: &[crate::hierarchy::BioAtomRow],
    input_format: crate::hierarchy::BioCoordinateFormat,
    residue_flag: u8,
    atom_flag: u8,
) -> Result<Option<BioRowFirstHit>, BioRowTraverseError> {
    // Gemmi✔️✔️:   std::pair<Model*, CRA> first(Structure& st) const {
    // Gemmi✔️✔️:     for (Model& model : st.models) {
    // Gemmi✔️✔️:       CRA cra = first_in_model(model);
    // Gemmi✔️✔️:       if (cra.chain)
    // Gemmi✔️✔️:          return {&model, cra};
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:     return {nullptr, {nullptr, nullptr, nullptr}};
    // Gemmi✔️:   }
    // Behavior: models are visited in stored order — no sorting by model
    // number or any other identity — and a model counts as a hit only
    // when its R09 traversal yields an atom. Each model delegates to
    // `bio_first_in_model_rows` (which applies the model/chain gates and
    // the row predicates), so missing-representation errors propagate
    // structurally from whichever row fails. The returned model ID is the
    // ORIGINAL row table index. The source's null-model return maps to
    // `None`.
    // Complexity: worst case one R09 scan per model with early exit; no
    // allocation.
    for (model_offset, model) in models.iter().enumerate() {
        if let Some(hit) = bio_first_in_model_rows(
            selection,
            model,
            chains,
            residues,
            atoms,
            input_format,
            residue_flag,
            atom_flag,
        )? {
            return Ok(Some(BioRowFirstHit {
                model: crate::hierarchy::BioModelId::new(model_offset as u32),
                cra: hit,
            }));
        }
    }
    Ok(None)
}

/// Lazy selected-atom-ID iterator over the actual BIO row hierarchy
/// (BIO-ROWS R11), mirroring Gemmi `FilterProxy`/`Selection::atoms`
/// laziness: yields only accepted atoms in exact source order, with no
/// intermediate `Vec` and no hierarchy snapshot — just borrowed slices
/// and cursor positions.
///
/// Model and chain gates propagate their typed missing-representation
/// errors structurally as `Some(Err(..))`, after which the iterator is
/// exhausted (the error is terminal, matching fail-closed gating).
pub(crate) struct BioSelectedAtomIds<'a> {
    selection: &'a Selection,
    models: &'a [crate::hierarchy::BioModelRow],
    chains: &'a [crate::hierarchy::BioChainRow],
    residues: &'a [crate::hierarchy::BioResidueRow],
    atoms: &'a [crate::hierarchy::BioAtomRow],
    input_format: crate::hierarchy::BioCoordinateFormat,
    residue_flag: u8,
    atom_flag: u8,
    stage: CursorStage,
    model_index: usize,
    chain_index: usize,
    residue_index: usize,
    atom_index: usize,
    error: Option<BioRowTraverseError>,
}

/// Lazy parent-entry stage (BIO-ROWS-CURSOR): each parent gate runs
/// exactly once per visited parent row, mirroring the source's nested
/// loops — a parent is re-gated only after its children are exhausted
/// and the cursor advances to the NEXT row of that level.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CursorStage {
    /// Next action: gate `models[model_index]` (or finish).
    AtModel,
    /// Model accepted; gate `chains[chain_index]` (or advance model).
    AtChain,
    /// Chain accepted; gate `residues[residue_index]` (or advance chain).
    AtResidue,
    /// Residue accepted; walk its atom span one atom at a time.
    InResidueAtoms,
}

impl<'a> BioSelectedAtomIds<'a> {
    /// Constructs the iterator over the full model/chain/residue/atom
    /// tables (BIO-ROWS R11). The caller supplies the tables and the
    /// explicit residue/atom flag bytes.
    pub(crate) fn over_structure(
        selection: &'a Selection,
        models: &'a [crate::hierarchy::BioModelRow],
        chains: &'a [crate::hierarchy::BioChainRow],
        residues: &'a [crate::hierarchy::BioResidueRow],
        atoms: &'a [crate::hierarchy::BioAtomRow],
        input_format: crate::hierarchy::BioCoordinateFormat,
        residue_flag: u8,
        atom_flag: u8,
    ) -> Self {
        // Gemmi✔️✔️: template<typename Filter, typename Value>
        // Gemmi✔️✔️: struct FilterProxy {
        // Gemmi✔️✔️:   const Filter& filter;
        // Gemmi✔️✔️:   std::vector<Value>& vec;
        // Gemmi✔️✔️:   using iterator = FilterIter<Filter, std::vector<Value>, Value>;
        // Gemmi✔️✔️:   iterator begin() { return {{&filter, &vec, 0}}; }
        // Gemmi✔️✔️:   iterator end() { return {{&filter, &vec, vec.size()}}; }
        // Gemmi✔️✔️: };
        // Gemmi✔️✔️: FilterProxy<Selection, Atom> atoms(Residue& residue) const {
        // Gemmi✔️✔️:   return {*this, residue.atoms};
        // Gemmi✔️✔️: }
        // Behavior: pure cursor state over borrowed slices; construction
        // performs no pre-scan and no allocation. The source proxy's
        // constructor pre-skip happens lazily inside the first increment
        // chain, not eagerly here.
        // Complexity: O(1) construction.
        Self {
            selection,
            models,
            chains,
            residues,
            atoms,
            input_format,
            residue_flag,
            atom_flag,
            stage: CursorStage::AtModel,
            model_index: 0,
            chain_index: 0,
            residue_index: 0,
            atom_index: 0,
            error: None,
        }
    }
}

impl<'a> Iterator for BioSelectedAtomIds<'a> {
    type Item = Result<crate::hierarchy::BioAtomId, BioRowTraverseError>;

    fn next(&mut self) -> Option<Self::Item> {
        // Gemmi✔️✔️:   FilterIterPolicy(const Filter* filter, Vector* vec, std::size_t pos)
        // Gemmi✔️✔️:       : filter_(filter), vec_(vec), pos_(pos) {
        // Gemmi✔️✔️:     while (pos_ != vec_->size() && !matches(pos_))
        // Gemmi✔️✔️:       ++pos_;
        // Gemmi✔️✔️:   }
        // Gemmi✔️✔️:   bool matches(std::size_t p) const { return filter_->matches((*vec_)[p]); }
        // Gemmi✔️✔️:   void increment() { while (++pos_ < vec_->size() && !matches(pos_)) {} }
        // Gemmi✔️✔️:   CRA first_in_model(Model& model) const {
        // Gemmi✔️✔️:     if (matches(model))
        // Gemmi✔️✔️:       for (Chain& chain : model.chains) {
        // Gemmi✔️✔️:         if (matches(chain))
        // Gemmi✔️✔️:           for (Residue& res : chain.residues) {
        // Gemmi✔️✔️:             if (matches(res))
        // Gemmi✔️✔️:               for (Atom& atom : res.atoms) {
        // Gemmi✔️✔️:                 if (matches(atom))
        // Gemmi✔️✔️:                   return {&chain, &res, &atom};
        // Gemmi✔️✔️:               }
        // Gemmi✔️✔️:           }
        // Gemmi✔️✔️:       }
        // Gemmi✔️✔️:     return {nullptr, nullptr, nullptr};
        // Gemmi✔️✔️:   }
        // Behavior (BIO-ROWS-CURSOR correction): explicit parent-entry
        // stages reproduce the source's nested-loop shape — the model
        // gate runs once per visited model row, the chain gate once per
        // visited chain row OF AN ACCEPTED MODEL, the residue gate once
        // per visited residue row OF AN ACCEPTED CHAIN, and the atom
        // predicate once per atom of an accepted residue. A parent row is
        // never re-gated for a second atom of the same child: after
        // descending, the parent's stage is re-entered only when the
        // child span is exhausted and the parent index advances. Yields
        // accepted atoms in exact source order; a structural
        // missing-representation error is returned once and the iterator
        // is exhausted afterwards (sticky). No eager full scan, no
        // hierarchy copy, no predicate-result cache.
        // Complexity: amortized O(1) per next() call — every row gate is
        // executed at most once across the entire iteration (one gate per
        // visited parent row), restoring the source's per-parent gating.
        if self.error.is_some() {
            return None;
        }
        loop {
            match self.stage {
                CursorStage::AtModel => {
                    if self.model_index >= self.models.len() {
                        return None;
                    }
                    match bio_model_row_matches(self.selection, &self.models[self.model_index]) {
                        Ok(true) => {
                            self.chain_index =
                                self.models[self.model_index].chain_span().start() as usize;
                            self.stage = CursorStage::AtChain;
                        }
                        Ok(false) => {
                            self.model_index += 1;
                        }
                        Err(e) => {
                            self.error = Some(BioRowTraverseError::Model(e));
                            return self.error.map(Err);
                        }
                    }
                }
                CursorStage::AtChain => {
                    let span = self.models[self.model_index].chain_span();
                    let end = (span.start() + span.len()) as usize;
                    if self.chain_index >= end {
                        self.model_index += 1;
                        self.stage = CursorStage::AtModel;
                        continue;
                    }
                    match bio_chain_row_matches(
                        self.selection,
                        &self.chains[self.chain_index],
                        self.input_format,
                    ) {
                        Ok(true) => {
                            self.residue_index =
                                self.chains[self.chain_index].residue_span().start() as usize;
                            self.stage = CursorStage::AtResidue;
                        }
                        Ok(false) => {
                            self.chain_index += 1;
                        }
                        Err(e) => {
                            self.error = Some(BioRowTraverseError::Chain(e));
                            return self.error.map(Err);
                        }
                    }
                }
                CursorStage::AtResidue => {
                    let rspan = self.chains[self.chain_index].residue_span();
                    let rend = (rspan.start() + rspan.len()) as usize;
                    if self.residue_index >= rend {
                        self.chain_index += 1;
                        self.stage = CursorStage::AtChain;
                        continue;
                    }
                    if bio_residue_row_matches(
                        self.selection,
                        &self.residues[self.residue_index],
                        self.input_format,
                        self.residue_flag,
                    ) {
                        self.atom_index =
                            self.residues[self.residue_index].atom_span().start() as usize;
                        self.stage = CursorStage::InResidueAtoms;
                    } else {
                        self.residue_index += 1;
                    }
                }
                CursorStage::InResidueAtoms => {
                    let aspan = self.residues[self.residue_index].atom_span();
                    let aend = (aspan.start() + aspan.len()) as usize;
                    if self.atom_index >= aend {
                        self.residue_index += 1;
                        self.stage = CursorStage::AtResidue;
                        continue;
                    }
                    let atom_id = crate::hierarchy::BioAtomId::new(self.atom_index as u32);
                    let accepted = bio_atom_row_matches(
                        self.selection,
                        &self.atoms[self.atom_index],
                        self.input_format,
                        self.atom_flag,
                    );
                    self.atom_index += 1;
                    if accepted {
                        return Some(Ok(atom_id));
                    }
                }
            }
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) struct SourceModel<'a> {
    /// Gemmi `Model::num` (model.hpp): the source `int` identity.
    pub(crate) num: i32,
    pub(crate) chains: &'a [SourceChain<'a>],
}

/// Source-shaped borrowed Chain input (`Chain::name` plus children).
#[derive(Clone, Copy)]
pub(crate) struct SourceChain<'a> {
    /// Gemmi `Chain::name`: the exact source chain-name string.
    pub(crate) name: &'a str,
    pub(crate) residues: &'a [SourceResidue<'a>],
}

/// Source-shaped borrowed Residue input.
#[derive(Clone, Copy)]
pub(crate) struct SourceResidue<'a> {
    pub(crate) entity_type: EntityKind,
    pub(crate) name: &'a str,
    /// Gemmi `Residue::seqid` raw number (`INT_MIN` when unset).
    pub(crate) seqid_num: i32,
    /// Gemmi `Residue::seqid` icode byte (`b' '` when unset).
    pub(crate) icode: u8,
    /// Gemmi `Residue::flag`: explicit scalar custom flag byte; BIO rows
    /// store no custom flag and het_flag/calc_flag substitution is
    /// forbidden (decision queue).
    pub(crate) flag: u8,
    pub(crate) atoms: &'a [SourceAtom<'a>],
}

/// Source-shaped borrowed Atom input.
#[derive(Clone, Copy)]
pub(crate) struct SourceAtom<'a> {
    /// Gemmi `Atom::name`: exact logical atom name, unnormalized.
    pub(crate) name: &'a str,
    /// Gemmi `Atom::element.ordinal()`: El ordinal (X=0..Og=118, D=119).
    pub(crate) element: u8,
    /// Gemmi `Atom::altloc`: `0` when unset.
    pub(crate) altloc: u8,
    /// Gemmi `Atom::flag`: explicit scalar custom flag byte.
    pub(crate) flag: u8,
    pub(crate) occ: f32,
    pub(crate) b_iso: f32,
}

/// Source-shaped CRA (chain/residue/atom) reference triple, mirroring
/// Gemmi's `CRA` of optional ancestor pointers.
pub(crate) struct SourceCra<'a> {
    pub(crate) chain: Option<&'a SourceChain<'a>>,
    pub(crate) residue: Option<&'a SourceResidue<'a>>,
    pub(crate) atom: Option<&'a SourceAtom<'a>>,
}

/// Gemmi `Selection` (select.hpp:20-182) carrier. The CID parser is not
/// part of this packet: fields are constructed directly, with the source's
/// default member state mirrored by [`Default`].
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Selection {
    pub(crate) mdl: i32,
    pub(crate) chain_ids: SelectionList,
    pub(crate) from_seqid: SelectionSequenceId,
    pub(crate) to_seqid: SelectionSequenceId,
    pub(crate) residue_names: SelectionList,
    pub(crate) entity_types: SelectionList,
    /// `std::array<char, 6>` in the source; indexed by `entity_kind_index`.
    pub(crate) et_flags: [bool; 6],
    pub(crate) atom_names: SelectionList,
    /// Source empty vector (`*`) ↔ `None`.
    pub(crate) elements: Option<GemmiElementMask>,
    pub(crate) altlocs: SelectionList,
    pub(crate) residue_flags: SelectionFlagList,
    pub(crate) atom_flags: SelectionFlagList,
    pub(crate) atom_inequalities: Vec<SelectionAtomInequality>,
}

impl Default for Selection {
    /// Gemmi member initializers: `mdl = 0`, `List` defaults `all = true`
    /// (select.hpp:23), `from_seqid = {INT_MIN, '*'}`, `to_seqid =
    /// {INT_MAX, '*'}`, `et_flags{}` zero-fills and the element vector is
    /// empty; flags/inequalities start empty.
    fn default() -> Self {
        Self {
            mdl: 0,
            chain_ids: SelectionList::default_all(),
            from_seqid: SelectionSequenceId {
                seqnum: i32::MIN,
                icode: b'*',
            },
            to_seqid: SelectionSequenceId {
                seqnum: i32::MAX,
                icode: b'*',
            },
            residue_names: SelectionList::default_all(),
            entity_types: SelectionList::default_all(),
            et_flags: [false; 6],
            atom_names: SelectionList::default_all(),
            elements: None,
            altlocs: SelectionList::default_all(),
            residue_flags: SelectionFlagList::default(),
            atom_flags: SelectionFlagList::default(),
            atom_inequalities: Vec::new(),
        }
    }
}

impl Selection {
    /// Gemmi `Selection::matches(Model)` (select.hpp:107-109).
    pub(crate) fn matches_model(&self, model: &SourceModel<'_>) -> bool {
        // Gemmi✔️✔️:   bool matches(const Model& model) const {
        // Gemmi✔️✔️:     return mdl == 0 || mdl == model.num;
        // Gemmi✔️✔️:   }
        // Behavior: 0 selects all models; otherwise exact `int` equality
        // on the source model number. No default is invented for a missing
        // source number — the input already carries the Gemmi `Model::num`
        // int shape, and detached-row adaptation is a decision-queue item.
        // Complexity: one integer comparison, no allocation.
        self.mdl == 0 || self.mdl == model.num
    }

    /// Gemmi `Selection::matches(Chain)` (select.hpp:110-112).
    pub(crate) fn matches_chain(&self, chain: &SourceChain<'_>) -> bool {
        // Gemmi✔️✔️:   bool matches(const Chain& chain) const {
        // Gemmi✔️✔️:     return chain_ids.has(chain.name);
        // Gemmi✔️✔️:   }
        // Behavior: exact source chain-name string membership with the
        // List semantics (all/inverted/comma members); byte-for-byte
        // comparison, no case folding, no trimming, and no row-ID,
        // label/auth or normalized-name substitution — the canonical
        // name-conversion owner is a decision-queue item, so the input
        // carries the Gemmi `Chain::name` string exactly.
        // Complexity: the List::has single scan over the stored list.
        self.chain_ids.has(chain.name)
    }

    /// Gemmi `Selection::matches(Residue)` (select.hpp:113-120).
    pub(crate) fn matches_residue(&self, residue: &SourceResidue<'_>) -> bool {
        // Gemmi✔️✔️:   bool matches(const Residue& res) const {
        // Gemmi✔️✔️:     return (entity_types.all || et_flags[(int)res.entity_type]) &&
        // Gemmi✔️✔️:            residue_names.has(res.name) &&
        // Gemmi✔️✔️:            from_seqid.compare(res.seqid) <= 0 &&
        // Gemmi✔️✔️:            to_seqid.compare(res.seqid) >= 0 &&
        // Gemmi✔️✔️:            residue_flags.has(res.flag);
        // Gemmi✔️✔️:   }
        // Behavior: five conjuncts short-circuit in source order. The
        // entity-type test is the `all` gate OR the et_flags ordinal: CK
        // `EntityKind` declaration order equals Gemmi `EntityType` ordinals
        // 0..4, so `entity_kind_index` is the (int) cast; the sixth array
        // slot is unreachable from both type systems. Sequence bounds are
        // inclusive via the source compare (wildcard `*` selector icodes
        // and INT_MIN/INT_MAX sentinel endpoints behave exactly as in
        // compare). The flag byte is an explicit scalar input — BIO rows
        // carry no custom residue flag and het_flag/calc_flag substitution
        // is forbidden (decision queue).
        // Complexity: one ordinal index, two list scans, two seqid
        // compares, one flag scan; all constant-space, no allocation.
        (self.entity_types.all || self.et_flags[entity_kind_index(residue.entity_type)])
            && self.residue_names.has(residue.name)
            && self.from_seqid.compare(residue.seqid_num, residue.icode) <= 0
            && self.to_seqid.compare(residue.seqid_num, residue.icode) >= 0
            && self.residue_flags.has(residue.flag)
    }

    /// Gemmi `Selection::matches(Atom)` (select.hpp:121-130).
    pub(crate) fn matches_atom(&self, atom: &SourceAtom<'_>) -> bool {
        // Gemmi✔️✔️:   bool matches(const Atom& a) const {
        // Gemmi✔️✔️:     return atom_names.has(a.name) &&
        // Gemmi✔️✔️:            (elements.empty() || elements[a.element.ordinal()]) &&
        // Gemmi✔️✔️:            (altlocs.all || altlocs.has(std::string(a.altloc ? 1 : 0, a.altloc))) &&
        // Gemmi✔️✔️:            atom_flags.has(a.flag) &&
        // Gemmi✔️✔️:            std::all_of(atom_inequalities.begin(), atom_inequalities.end(),
        // Gemmi✔️✔️:                        [&](const AtomInequality& i) { return i.matches(a); });
        // Gemmi✔️✔️:   }
        // Behavior: five conjuncts short-circuit in source order. The
        // element mask is indexed by Gemmi El ordinal (X=0..Og=118; the D
        // slot 119 stays addressable exactly as parsed but cannot arise
        // from a CK `Element` atomic number); `None` is the source's empty
        // mask = no element restriction. The altloc comparison reproduces
        // `std::string(a.altloc ? 1 : 0, a.altloc)`: the byte 0 sentinel
        // becomes the empty name, any other byte a one-byte name compared
        // byte-for-byte (no UTF-8 assumption, no normalization). The flag
        // byte is an explicit scalar; every inequality must hold (all_of
        // over an empty vector is true).
        // Complexity: two list scans, one mask index, one inequality pass;
        // the one-byte altloc name is stack-local — no allocation.
        let altloc_name: [u8; 1] = [atom.altloc];
        let altloc_bytes: &[u8] = if atom.altloc == 0 { &[] } else { &altloc_name };
        self.atom_names.has(atom.name)
            && self
                .elements
                .as_ref()
                .is_none_or(|mask| mask[usize::from(atom.element)])
            && (self.altlocs.all || self.altlocs.has_bytes(altloc_bytes))
            && self.atom_flags.has(atom.flag)
            && self
                .atom_inequalities
                .iter()
                .all(|inequality| inequality.matches(atom.occ, atom.b_iso))
    }

    /// Gemmi `Selection::matches(CRA)` (select.hpp:131-135).
    pub(crate) fn matches_cra(&self, cra: &SourceCra<'_>) -> bool {
        // Gemmi✔️✔️:   bool matches(const CRA& cra) const {
        // Gemmi✔️✔️:     return (cra.chain == nullptr || matches(*cra.chain)) &&
        // Gemmi✔️✔️:            (cra.residue == nullptr || matches(*cra.residue)) &&
        // Gemmi✔️✔️:            (cra.atom == nullptr || matches(*cra.atom));
        // Gemmi✔️✔️:   }
        // Behavior: each absent ancestor (`nullptr` in the source, `None`
        // here) short-circuits its conjunct to true; each present ancestor
        // delegates to the corresponding level predicate, evaluated in
        // source order chain -> residue -> atom.
        // Complexity: at most three predicate evaluations, no allocation.
        (cra.chain.is_none_or(|chain| self.matches_chain(chain)))
            && (cra
                .residue
                .is_none_or(|residue| self.matches_residue(residue)))
            && (cra.atom.is_none_or(|atom| self.matches_atom(atom)))
    }

    /// Gemmi `Selection::models` (select.hpp:137-139) as a lazy ordered
    /// filter over borrowed source-shaped model rows.
    pub(crate) fn models<'a>(
        &'a self,
        models: &'a [SourceModel<'a>],
    ) -> SelectionFilterIter<'a, SourceModel<'a>, impl Fn(&SourceModel<'a>) -> bool + 'a> {
        // Gemmi✔️✔️:   FilterProxy<Selection, Model> models(Structure& st) const {
        // Gemmi✔️✔️:     return {*this, st.models};
        // Gemmi✔️✔️:   }
        // Behavior/complexity: see SelectionFilterIter; the proxy only
        // pairs the selection with the borrowed model slice.
        SelectionFilterIter::new(models, move |model| self.matches_model(model))
    }

    /// Gemmi `Selection::chains` (select.hpp:140-142).
    pub(crate) fn chains<'a>(
        &'a self,
        chains: &'a [SourceChain<'a>],
    ) -> SelectionFilterIter<'a, SourceChain<'a>, impl Fn(&SourceChain<'a>) -> bool + 'a> {
        // Gemmi✔️✔️:   FilterProxy<Selection, Chain> chains(Model& model) const {
        // Gemmi✔️✔️:     return {*this, model.chains};
        // Gemmi✔️✔️:   }
        SelectionFilterIter::new(chains, move |chain| self.matches_chain(chain))
    }

    /// Gemmi `Selection::residues` (select.hpp:143-145).
    pub(crate) fn residues<'a>(
        &'a self,
        residues: &'a [SourceResidue<'a>],
    ) -> SelectionFilterIter<'a, SourceResidue<'a>, impl Fn(&SourceResidue<'a>) -> bool + 'a> {
        // Gemmi✔️✔️:   FilterProxy<Selection, Residue> residues(Chain& chain) const {
        // Gemmi✔️✔️:     return {*this, chain.residues};
        // Gemmi✔️✔️:   }
        SelectionFilterIter::new(residues, move |residue| self.matches_residue(residue))
    }

    /// Gemmi `Selection::atoms` (select.hpp:146-148).
    pub(crate) fn atoms<'a>(
        &'a self,
        atoms: &'a [SourceAtom<'a>],
    ) -> SelectionFilterIter<'a, SourceAtom<'a>, impl Fn(&SourceAtom<'a>) -> bool + 'a> {
        // Gemmi✔️✔️:   FilterProxy<Selection, Atom> atoms(Residue& residue) const {
        // Gemmi✔️✔️:     return {*this, residue.atoms};
        // Gemmi✔️✔️:   }
        SelectionFilterIter::new(atoms, move |atom| self.matches_atom(atom))
    }

    /// Gemmi `Selection::first_in_model` (select.hpp:150-163).
    pub(crate) fn first_in_model<'m>(&self, model: &'m SourceModel<'m>) -> SourceCra<'m> {
        // Gemmi✔️✔️:   CRA first_in_model(Model& model) const {
        // Gemmi✔️✔️:     if (matches(model))
        // Gemmi✔️✔️:       for (Chain& chain : model.chains) {
        // Gemmi✔️✔️:         if (matches(chain))
        // Gemmi✔️✔️:           for (Residue& res : chain.residues) {
        // Gemmi✔️✔️:             if (matches(res))
        // Gemmi✔️✔️:               for (Atom& atom : res.atoms) {
        // Gemmi✔️✔️:                 if (matches(atom))
        // Gemmi✔️✔️:                   return {&chain, &res, &atom};
        // Gemmi✔️✔️:               }
        // Gemmi✔️✔️:           }
        // Gemmi✔️✔️:       }
        // Gemmi✔️✔️:     return {nullptr, nullptr, nullptr};
        // Gemmi✔️✔️:   }
        // Behavior: the model gate runs first (a rejected model short-
        // circuits to the all-null CRA); chains, residues and atoms are
        // scanned in exact source order, each gated by its predicate
        // before descending. A result is returned ONLY when an atom
        // matches: a matching residue whose atom loop completes without a
        // match does not return — the residue loop continues to the next
        // residue (the source has no residue-only result).
        // Complexity: worst case one full hierarchy scan with one
        // predicate call per row; early exit on the first matching atom;
        // no allocation.
        if self.matches_model(model) {
            for chain in model.chains {
                if self.matches_chain(chain) {
                    for residue in chain.residues {
                        if self.matches_residue(residue) {
                            for atom in residue.atoms {
                                if self.matches_atom(atom) {
                                    return SourceCra {
                                        chain: Some(chain),
                                        residue: Some(residue),
                                        atom: Some(atom),
                                    };
                                }
                            }
                        }
                    }
                }
            }
        }
        SourceCra {
            chain: None,
            residue: None,
            atom: None,
        }
    }

    /// Gemmi `Selection::first` (select.hpp:165-172).
    pub(crate) fn first<'m>(
        &self,
        models: &'m [SourceModel<'m>],
    ) -> Option<(&'m SourceModel<'m>, SourceCra<'m>)> {
        // Gemmi✔️✔️:   std::pair<Model*, CRA> first(Structure& st) const {
        // Gemmi✔️✔️:     for (Model& model : st.models) {
        // Gemmi✔️✔️:       CRA cra = first_in_model(model);
        // Gemmi✔️✔️:       if (cra.chain)
        // Gemmi✔️:          return {&model, cra};
        // Gemmi✔️✔️:     }
        // Gemmi✔️✔️:     return {nullptr, {nullptr, nullptr, nullptr}};
        // Gemmi✔️✔️:   }
        // Behavior: models are visited in stored order; a model counts as
        // a hit only when first_in_model yields a chain (i.e. an atom
        // matched). The source's null-model return maps to `None`. The
        // chain-presence check (`cra.chain`) is exactly first_in_model's
        // non-null signal — residue/atom are set whenever chain is.
        // Complexity: worst case one first_in_model scan per model; early
        // exit on the first model with a match; no allocation.
        for model in models {
            let cra = self.first_in_model(model);
            if cra.chain.is_some() {
                return Some((model, cra));
            }
        }
        None
    }
}

/// Gemmi `FilterIterPolicy`/`FilterProxy` (iterator.hpp:180-211): a lazy,
/// allocation-free forward iterator yielding only the rows the predicate
/// accepts, in exact source order.
pub(crate) struct SelectionFilterIter<'a, T: 'a, F: Fn(&T) -> bool + 'a> {
    rows: &'a [T],
    pos: usize,
    predicate: F,
}

impl<'a, T, F: Fn(&T) -> bool + 'a> SelectionFilterIter<'a, T, F> {
    pub(crate) fn new(rows: &'a [T], predicate: F) -> Self {
        // Gemmi✔️✔️:   FilterIterPolicy(const Filter* filter, Vector* vec, std::size_t pos)
        // Gemmi✔️✔️:       : filter_(filter), vec_(vec), pos_(pos) {
        // Gemmi✔️✔️:     while (pos_ != vec_->size() && !matches(pos_))
        // Gemmi✔️✔️:       ++pos_;
        // Gemmi✔️✔️:   }
        // Behavior: the source constructor pre-skips non-matching rows so
        // begin() lands on the first match (or the end position); this port
        // performs the identical skip lazily on the first `next()` call,
        // which yields the same sequence with no observable difference for
        // a forward scan. `end()` is the row count (vec.size()).
        // Complexity: O(1) construction, one scan with no allocation; each
        // row is tested exactly once across the whole iteration.
        Self {
            rows,
            pos: 0,
            predicate,
        }
    }
}

impl<'a, T, F: Fn(&T) -> bool + 'a> Iterator for SelectionFilterIter<'a, T, F> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        // Gemmi✔️✔️:   void increment() { while (++pos_ < vec_->size() && !matches(pos_)) {} }
        // Gemmi✔️✔️:   Value& dereference() { return (*vec_)[pos_]; }
        // Behavior: advance past non-matching rows and yield a borrow of
        // the next matching row; exhaustion returns None (pos == len,
        // equal to the source end iterator).
        // Complexity: amortized one predicate call per row, no allocation.
        while self.pos < self.rows.len() {
            let index = self.pos;
            self.pos += 1;
            if (self.predicate)(&self.rows[index]) {
                return Some(&self.rows[index]);
            }
        }
        None
    }
}

// ---- BIO-CID C01: canonical detached selection data ----

/// Structured exact-bytes error when a CID list carrier does not hold
/// valid UTF-8. The CID grammar only produces valid-UTF-8 list slices
/// (BIO-SEL C04 conversion rule), so this error is fail-closed evidence,
/// never a lossy or empty substitute.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CidListUtf8Error {
    /// Which selection member the invalid bytes came from.
    pub field: &'static str,
    /// The exact invalid bytes, preserved.
    pub bytes: Vec<u8>,
}

/// Owned construction vocabulary for `BioSelectionData`: exactly the
/// ten CID-reachable member groups of Gemmi's `Selection`. The CID
/// grammar never sets `residue_flags`/`atom_flags`, so flag members are
/// deliberately absent and stay source-empty after construction.
pub struct BioSelectionParts {
    // Gemmi❗✔️:   int mdl = 0;            // 0 = all
    // Gemmi❗✔️:   List chain_ids;
    // Gemmi❗✔️:   SequenceId from_seqid = {INT_MIN, '*'};
    // Gemmi❗✔️:   SequenceId to_seqid = {INT_MAX, '*'};
    // Gemmi❗✔️:   List residue_names;
    // Gemmi❗✔️:   List entity_types;
    // Gemmi❗✔️:   // array corresponding to enum EntityType
    // Gemmi❗✔️:   std::array<char, 6> et_flags;
    // Gemmi❗✔️:   List atom_names;
    // Gemmi❗✔️:   std::vector<char> elements;
    // Gemmi❗✔️:   List altlocs;
    // Gemmi❗✔️:   FlagList residue_flags;
    // Gemmi❗✔️:   FlagList atom_flags;
    // Gemmi❗✔️:   std::vector<AtomInequality> atom_inequalities;
    //
    // Behavior review: this Rust struct assembles exactly the ten
    // CID-reachable member groups with the same carrier semantics
    // (List = SelectionCidList carrier, SequenceId, et_flags, element
    // mask, AtomInequality); the two FlagList members have NO part here
    // because the CID grammar never writes flags (packet fixed
    // contract; select.hpp:103-104 members stay source-empty). This is a
    // construction vocabulary, not a copy of one C++ function, so no
    // single-line behavioral equivalence is claimed beyond the member
    // set. Complexity review: plain owned moves, no scan or clone.
    pub mdl: i32,
    pub chain_ids: SelectionCidList,
    pub from_seqid: SelectionSequenceId,
    pub to_seqid: SelectionSequenceId,
    pub residue_names: SelectionCidList,
    pub entity_types: SelectionCidList,
    /// Indexed by `entity_kind_index`, like `std::array<char, 6>`.
    pub et_flags: [bool; 6],
    pub atom_names: SelectionCidList,
    /// Source empty vector (`*`) ↔ `None`.
    pub elements: Option<GemmiElementMask>,
    pub altlocs: SelectionCidList,
    pub atom_inequalities: Vec<SelectionAtomInequality>,
}

impl Default for BioSelectionParts {
    /// Same member initializers as `Selection::default`: `mdl = 0`,
    /// `List` defaults `all = true` (select.hpp:23),
    /// `from_seqid = {INT_MIN, '*'}`, `to_seqid = {INT_MAX, '*'}`
    /// (select.hpp:94-95), `et_flags{}` zero-fills, `elements` empty.
    fn default() -> Self {
        Self {
            mdl: 0,
            chain_ids: SelectionCidList {
                all: true,
                inverted: false,
                list: Vec::new(),
            },
            from_seqid: SelectionSequenceId {
                seqnum: i32::MIN,
                icode: b'*',
            },
            to_seqid: SelectionSequenceId {
                seqnum: i32::MAX,
                icode: b'*',
            },
            residue_names: SelectionCidList {
                all: true,
                inverted: false,
                list: Vec::new(),
            },
            entity_types: SelectionCidList {
                all: true,
                inverted: false,
                list: Vec::new(),
            },
            et_flags: [false; 6],
            atom_names: SelectionCidList {
                all: true,
                inverted: false,
                list: Vec::new(),
            },
            elements: None,
            altlocs: SelectionCidList {
                all: true,
                inverted: false,
                list: Vec::new(),
            },
            atom_inequalities: Vec::new(),
        }
    }
}

/// Borrowed read-only view over a `BioSelectionData`'s members. Gemmi
/// code reads `Selection` members directly; this Rust view supplies the
/// same read access without exposing a mutable engine, so there is no
/// single C++ counterpart to anchor — the member set mirrors
/// `BioSelectionParts` one-to-one. No clone: every field borrows.
pub struct BioSelectionPartsView<'a> {
    pub mdl: i32,
    pub chain_ids: &'a SelectionList,
    pub from_seqid: &'a SelectionSequenceId,
    pub to_seqid: &'a SelectionSequenceId,
    pub residue_names: &'a SelectionList,
    pub entity_types: &'a SelectionList,
    pub et_flags: &'a [bool; 6],
    pub atom_names: &'a SelectionList,
    pub elements: Option<&'a GemmiElementMask>,
    pub altlocs: &'a SelectionList,
    pub atom_inequalities: &'a [SelectionAtomInequality],
}

/// Typed failure of the detached selected-atom query (BIO-CID C24):
/// a structural parent gate encountered a row without the canonical
/// source representation the pinned reader paths require (the terminal
/// lazy-cursor cause is retained privately; its message is surfaced
/// through Display). This is fail-closed evidence, never an empty-ID
/// substitute.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BioSelectionMatchError {
    cause: BioRowTraverseError,
}

impl BioSelectionMatchError {
    /// Wraps the terminal lazy-cursor failure (crate-internal cause).
    pub(crate) fn from_traverse_error(cause: BioRowTraverseError) -> Self {
        Self { cause }
    }
}

impl std::fmt::Display for BioSelectionMatchError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "bio selection match failed: {}", self.cause)
    }
}

impl std::error::Error for BioSelectionMatchError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        // Structured cause retention (BIO-C24-CLOSE): expose the stored
        // traverse error through the standard chain without any new
        // public field/type.
        Some(&self.cause as &(dyn std::error::Error + 'static))
    }
}

/// The one canonical detached compiled-selection state (BIO-CID C01).
/// Every construction path (CID decode in IO, root adapter) reuses this
/// single representation; predicates and traversal live on the inner
/// `Selection` engine and are never duplicated here.
pub struct BioSelectionData {
    selection: Selection,
}

impl Default for BioSelectionData {
    /// `Selection() = default` (select.hpp:107): all-wildcard state.
    fn default() -> Self {
        Self {
            selection: Selection::default(),
        }
    }
}

impl BioSelectionData {
    /// Assemble the canonical state from owned parts. CID list bytes are
    /// converted to the `SelectionList` string form exactly: only a
    /// valid-UTF-8 carrier is accepted, and any invalid bytes return the
    /// exact bytes with their member label (no lossy fallback, no empty
    /// substitute). Flag members start source-empty because the parts
    /// vocabulary cannot carry them.
    pub fn from_parts(parts: BioSelectionParts) -> Result<Self, CidListUtf8Error> {
        let selection = Selection {
            mdl: parts.mdl,
            chain_ids: cid_list_into_selection_list(parts.chain_ids, "chain_ids")?,
            from_seqid: parts.from_seqid,
            to_seqid: parts.to_seqid,
            residue_names: cid_list_into_selection_list(parts.residue_names, "residue_names")?,
            entity_types: cid_list_into_selection_list(parts.entity_types, "entity_types")?,
            et_flags: parts.et_flags,
            atom_names: cid_list_into_selection_list(parts.atom_names, "atom_names")?,
            elements: parts.elements,
            altlocs: cid_list_into_selection_list(parts.altlocs, "altlocs")?,
            residue_flags: SelectionFlagList::default(),
            atom_flags: SelectionFlagList::default(),
            atom_inequalities: parts.atom_inequalities,
        };
        Ok(Self { selection })
    }

    /// Borrowed read-only view over every CID-reachable member. Pure
    /// borrows of the inner state — no clone, no rebuild. This IS the
    /// authorized external read boundary (with `from_parts` on the
    /// construction side); it never exposes `Source*` inputs.
    ///
    /// ```
    /// // Positive control (BIO-C01-PROOF): default construction and the
    /// // borrowed view are externally usable — all-wildcard members.
    /// let data = cosmolkit_bio::BioSelectionData::default();
    /// let parts = data.parts();
    /// assert_eq!(parts.mdl, 0);
    /// assert_eq!(parts.from_seqid.seqnum, i32::MIN);
    /// ```
    pub fn parts(&self) -> BioSelectionPartsView<'_> {
        BioSelectionPartsView {
            mdl: self.selection.mdl,
            chain_ids: &self.selection.chain_ids,
            from_seqid: &self.selection.from_seqid,
            to_seqid: &self.selection.to_seqid,
            residue_names: &self.selection.residue_names,
            entity_types: &self.selection.entity_types,
            et_flags: &self.selection.et_flags,
            atom_names: &self.selection.atom_names,
            elements: self.selection.elements.as_ref(),
            altlocs: &self.selection.altlocs,
            atom_inequalities: &self.selection.atom_inequalities,
        }
    }

    /// Detached selected-atom-ID query (BIO-CID C24): collects the
    /// ORIGINAL atom IDs of every atom accepted by this selection over
    /// the borrowed BIO row tables, in exact source order, through the
    /// existing lazy cursor. Typed failure: a structural parent gate
    /// with a missing canonical source representation terminates the
    /// query with [`BioSelectionMatchError`].
    #[must_use]
    pub fn selected_bio_atom_ids(
        &self,
        models: &[crate::hierarchy::BioModelRow],
        chains: &[crate::hierarchy::BioChainRow],
        residues: &[crate::hierarchy::BioResidueRow],
        atoms: &[crate::hierarchy::BioAtomRow],
        input_format: crate::hierarchy::BioCoordinateFormat,
    ) -> Result<Vec<crate::hierarchy::BioAtomId>, BioSelectionMatchError> {
        // Gemmi✔️✔️:   bool matches(const Model& model) const {
        // Gemmi✔️✔️:     return mdl == 0 || mdl == model.num;
        // Gemmi✔️✔️:   }
        // Gemmi✔️✔️:   bool matches(const Chain& chain) const {
        // Gemmi✔️✔️:     return chain_ids.has(chain.name);
        // Gemmi✔️✔️:   }
        // Gemmi✔️✔️:   bool matches(const Residue& res) const {
        // Gemmi✔️✔️:     return (entity_types.all || et_flags[(int)res.entity_type]) &&
        // Gemmi✔️✔️:            residue_names.has(res.name) &&
        // Gemmi✔️✔️:            from_seqid.compare(res.seqid) <= 0 &&
        // Gemmi✔️✔️:            to_seqid.compare(res.seqid) >= 0 &&
        // Gemmi✔️✔️:            residue_flags.has(res.flag);
        // Gemmi✔️✔️:   }
        // Gemmi✔️✔️:   bool matches(const Atom& a) const {
        // Gemmi✔️✔️:     return atom_names.has(a.name) &&
        // Gemmi✔️✔️:   FilterProxy<Selection, Model> models(Structure& st) const {
        // Gemmi✔️✔️:   FilterProxy<Selection, Chain> chains(Model& model) const {
        // Gemmi✔️✔️:   FilterProxy<Selection, Residue> residues(Chain& chain) const {
        // Gemmi✔️✔️:   FilterProxy<Selection, Atom> atoms(Residue& residue) const {
        // (select.hpp:110-148; the nested model/chain/residue/atom gating
        // executes inside the reused lazy cursor, BIO-ROWS-CURSOR shape.)
        //
        // Behavior review: the nested filtering and source-order yield
        // come from the EXISTING lazy cursor (BioSelectedAtomIds) over the
        // borrowed tables — no structural copy, no hierarchy snapshot, no
        // second predicate engine. CID-decoded selections carry
        // source-EMPTY residue/atom flag patterns (BIO-CID C01 contract:
        // the CID grammar never writes flags), and an empty FlagList
        // pattern accepts EVERY flag value, so the explicit zero flag
        // bytes handed to the cursor are inert — no invented flags. A
        // terminal structural missing-representation error from a parent
        // gate is mapped into the typed BioSelectionMatchError after the
        // cursor's sticky-error exhaustion. Complexity review: linear
        // visited-row gating plus the predicate scans (each list
        // membership is a scan over the stored names), plus one result
        // Vec allocation; no per-row String allocation — the row adapters
        // borrow their logical name views (BIO-C24-CLOSE).
        BioSelectedAtomIds::over_structure(
            &self.selection,
            models,
            chains,
            residues,
            atoms,
            input_format,
            0,
            0,
        )
        .collect::<Result<Vec<_>, _>>()
        .map_err(BioSelectionMatchError::from_traverse_error)
    }

    /// `Selection::matches(Model)` (select.hpp:110-112).
    ///
    /// Crate-private: the `SourceModel` input is source-shaped borrowed
    /// BIO vocabulary and must never become an external (root or IO)
    /// parameter; external matching goes through the root thin adapters.
    ///
    /// Proof note (BIO-C01-PROOF): the lookup below is argument-free, so
    /// the ONLY possible failure is the method's own visibility.
    ///
    /// ```compile_fail,E0624
    /// // E0624: `matches_model` is private to `cosmolkit-bio`.
    /// let _ = cosmolkit_bio::BioSelectionData::matches_model;
    /// ```
    pub(crate) fn matches_model(&self, model: &SourceModel<'_>) -> bool {
        self.selection.matches_model(model)
    }

    /// `Selection::matches(Chain)` (select.hpp:113-115).
    ///
    /// Crate-private: the `SourceChain` input is source-shaped borrowed
    /// BIO vocabulary and must never become an external parameter.
    ///
    /// Proof note (BIO-C01-PROOF): argument-free lookup; only visibility
    /// can fail.
    ///
    /// ```compile_fail,E0624
    /// // E0624: `matches_chain` is private to `cosmolkit-bio`.
    /// let _ = cosmolkit_bio::BioSelectionData::matches_chain;
    /// ```
    pub(crate) fn matches_chain(&self, chain: &SourceChain<'_>) -> bool {
        self.selection.matches_chain(chain)
    }

    /// `Selection::matches(Residue)` (select.hpp:116-123).
    ///
    /// Crate-private: the `SourceResidue` input is source-shaped borrowed
    /// BIO vocabulary and must never become an external parameter.
    ///
    /// Proof note (BIO-C01-PROOF): argument-free lookup; only visibility
    /// can fail.
    ///
    /// ```compile_fail,E0624
    /// // E0624: `matches_residue` is private to `cosmolkit-bio`.
    /// let _ = cosmolkit_bio::BioSelectionData::matches_residue;
    /// ```
    pub(crate) fn matches_residue(&self, residue: &SourceResidue<'_>) -> bool {
        self.selection.matches_residue(residue)
    }

    /// `Selection::matches(Atom)` (select.hpp:124-132).
    ///
    /// Crate-private: the `SourceAtom` input is source-shaped borrowed
    /// BIO vocabulary and must never become an external parameter.
    ///
    /// Proof note (BIO-C01-PROOF): argument-free lookup; only visibility
    /// can fail.
    ///
    /// ```compile_fail,E0624
    /// // E0624: `matches_atom` is private to `cosmolkit-bio`.
    /// let _ = cosmolkit_bio::BioSelectionData::matches_atom;
    /// ```
    pub(crate) fn matches_atom(&self, atom: &SourceAtom<'_>) -> bool {
        self.selection.matches_atom(atom)
    }
}

/// Exact carrier conversion: CID list bytes must be valid UTF-8
/// (`String::from_utf8`), mirroring the C++ `std::string` payload that
/// already holds those bytes. Invalid bytes fail closed with the exact
/// bytes and member label.
fn cid_list_into_selection_list(
    list: SelectionCidList,
    field: &'static str,
) -> Result<SelectionList, CidListUtf8Error> {
    let text = String::from_utf8(list.list).map_err(|error| CidListUtf8Error {
        field,
        bytes: error.into_bytes(),
    })?;
    Ok(SelectionList {
        all: list.all,
        inverted: list.inverted,
        list: text,
    })
}

#[cfg(test)]
mod bio_cid_c01_tests {
    use super::{
        BioSelectionData, BioSelectionParts, CidListUtf8Error, SelectionAtomInequality,
        SelectionCidList, SelectionList, SelectionSequenceId, SourceAtom, SourceChain, SourceModel,
        SourceResidue,
    };

    fn cid_list(list: &str) -> SelectionCidList {
        SelectionCidList {
            all: false,
            inverted: false,
            list: list.as_bytes().to_vec(),
        }
    }

    // Classifier boundary (BIO-C01-FIX Step 6): the exported canonical
    // table must equal the pinned elem.hpp:98-130 table at EVERY one of
    // the 120 ordinals. The expected literal below is transcribed
    // independently from the C++ source layout (per-row element
    // comments kept), not from GEMMI_IS_METAL or any mask built from it.
    #[test]
    fn bio_cid_c01_metal_table_matches_pinned_source_all_120_ordinals() {
        #[rustfmt::skip]
        const EXPECTED_FROM_ELEM_HPP: [bool; 120] = [
            // X     H     He
            false, false, false,
            // Li  Be     B      C      N      O      F     Ne
            true,  true,  false, false, false, false, false, false,
            // Na  Mg    Al     Si     P      S      Cl     Ar
            true,  true,  true,  false, false, false, false, false,
            // K   Ca    Sc    Ti    V     Cr    Mn    Fe    Co    Ni    Cu    Zn
            true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
            // Ga  Ge    As     Se     Br     Kr
            true,  true,  false, false, false, false,
            // Rb  Sr    Y     Zr    Nb    Mo    Tc    Ru    Rh    Pd    Ag    Cd
            true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
            // In  Sn    Sb    Te      I     Xe
            true,  true,  true,  false, false, false,
            // Cs  Ba    La    Ce    Pr    Nd    Pm    Sm    Eu    Gd    Tb    Dy
            true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
            // Ho  Er    Tm    Yb    Lu    Hf    Ta    W     Re    Os    Ir    Pt
            true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
            // Au  Hg    Tl    Pb    Bi    Po    At     Rn
            true,  true,  true,  true,  true,  true,  false, false,
            // Fr  Ra    Ac    Th    Pa    U     Np    Pu    Am    Cm    Bk    Cf
            true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
            // Es  Fm    Md    No    Lr    Rf    Db    Sg    Bh    Hs    Mt    Ds
            true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
            // Rg  Cn    Nh    Fl    Mc    Lv    Ts     Og
            true,  true,  true,  true,  true,  true,  false, false,
            // D
            false,
        ];
        for ordinal in 0..120 {
            assert_eq!(
                super::GEMMI_IS_METAL[ordinal],
                EXPECTED_FROM_ELEM_HPP[ordinal],
                "is_metal mismatch at El ordinal {ordinal} (X=0, H=1, Og=118, D=119)"
            );
        }
        // The 121st source slot (El::END) is deliberately not part of the
        // exported first-120 boundary; the source loop never reads it.
        assert_eq!(super::GEMMI_EL_END, 120);
    }

    // Default wildcard: `Selection() = default` leaves every CID member
    // at its select.hpp:92-105 initializer.
    #[test]
    fn bio_cid_c01_default_is_all_wildcard() {
        let data = BioSelectionData::default();
        let view = data.parts();
        assert_eq!(view.mdl, 0);
        for list in [
            view.chain_ids,
            view.residue_names,
            view.entity_types,
            view.atom_names,
            view.altlocs,
        ] {
            assert!(list.all);
            assert!(!list.inverted);
            assert!(list.list.is_empty());
        }
        assert_eq!(view.from_seqid.seqnum, i32::MIN);
        assert_eq!(view.from_seqid.icode, b'*');
        assert_eq!(view.to_seqid.seqnum, i32::MAX);
        assert_eq!(view.to_seqid.icode, b'*');
        assert!(view.et_flags.iter().all(|flag| !flag));
        assert_eq!(view.elements, None);
        assert!(view.atom_inequalities.is_empty());
        // Wildcard model matching: mdl = 0 selects every model number.
        let models = [SourceModel {
            num: 1,
            chains: &[],
        }];
        assert!(data.matches_model(&models[0]));
    }

    // Field and order preservation: every part survives assembly with
    // its exact value, and the borrowed view reads the same storage.
    #[test]
    fn bio_cid_c01_parts_field_and_order_preservation() {
        let mut et_flags = [false; 6];
        et_flags[2] = true;
        let mut elements = [false; super::GEMMI_EL_END];
        elements[6] = true; // C ordinal 6
        let parts = BioSelectionParts {
            mdl: 3,
            chain_ids: cid_list("A,B"),
            from_seqid: SelectionSequenceId {
                seqnum: 10,
                icode: b' ',
            },
            to_seqid: SelectionSequenceId {
                seqnum: 20,
                icode: b'A',
            },
            residue_names: cid_list("ALA,GLY"),
            entity_types: cid_list("polymer"),
            et_flags,
            atom_names: cid_list("CA,N"),
            elements: Some(elements),
            altlocs: cid_list("A"),
            atom_inequalities: vec![SelectionAtomInequality {
                property: b'q',
                relation: -1,
                value: 0.5,
            }],
        };
        let data = BioSelectionData::from_parts(parts).expect("valid UTF-8 parts");
        let view = data.parts();
        assert_eq!(view.mdl, 3);
        assert_eq!(view.chain_ids.list, "A,B");
        assert!(!view.chain_ids.all);
        assert_eq!(view.from_seqid.seqnum, 10);
        assert_eq!(view.from_seqid.icode, b' ');
        assert_eq!(view.to_seqid.seqnum, 20);
        assert_eq!(view.to_seqid.icode, b'A');
        assert_eq!(view.residue_names.list, "ALA,GLY");
        assert_eq!(view.entity_types.list, "polymer");
        assert!(view.et_flags[2]);
        assert!(view.et_flags.iter().filter(|flag| **flag).count() == 1);
        assert_eq!(view.atom_names.list, "CA,N");
        let elements = view.elements.expect("elements present");
        assert!(elements[6]);
        assert_eq!(view.altlocs.list, "A");
        let inequality = &view.atom_inequalities[0];
        assert_eq!(inequality.property, b'q');
        assert_eq!(inequality.relation, -1);
        assert_eq!(inequality.value, 0.5);
    }

    // No clone in the borrowed view: the view borrows the exact inner
    // member storage (pointer identity), and matching works through the
    // assembled state.
    #[test]
    fn bio_cid_c01_borrowed_view_is_not_a_clone() {
        let parts = BioSelectionParts {
            chain_ids: cid_list("A"),
            ..BioSelectionParts::default()
        };
        let data = BioSelectionData::from_parts(parts).expect("valid UTF-8 parts");
        let view = data.parts();
        assert!(std::ptr::eq(
            view.chain_ids as *const SelectionList,
            &data.selection.chain_ids as *const SelectionList
        ));
        assert!(std::ptr::eq(
            view.atom_inequalities.as_ptr(),
            data.selection.atom_inequalities.as_ptr()
        ));
        // Matching runs through the canonical engine with the assembled
        // chain list (byte-exact membership, no case folding).
        let chain = SourceChain {
            name: "A",
            residues: &[],
        };
        assert!(data.matches_chain(&chain));
        let other = SourceChain {
            name: "a",
            residues: &[],
        };
        assert!(!data.matches_chain(&other));
    }

    // Valid caller UTF-8 list boundaries: multibyte members round-trip
    // exactly through assembly; invalid bytes fail closed with the
    // exact bytes and member label, never a lossy substitute.
    #[test]
    fn bio_cid_c01_utf8_list_boundaries() {
        let parts = BioSelectionParts {
            residue_names: cid_list("ALA,\u{3a9}"),
            ..BioSelectionParts::default()
        };
        let data = BioSelectionData::from_parts(parts).expect("valid UTF-8 multibyte");
        assert_eq!(data.parts().residue_names.list, "ALA,\u{3a9}");

        let invalid = BioSelectionParts {
            atom_names: SelectionCidList {
                all: false,
                inverted: false,
                list: vec![b'C', 0xff],
            },
            ..BioSelectionParts::default()
        };
        match BioSelectionData::from_parts(invalid) {
            Err(CidListUtf8Error { field, bytes }) => {
                assert_eq!(field, "atom_names");
                assert_eq!(bytes, vec![b'C', 0xff]);
            }
            Ok(_) => panic!("invalid UTF-8 must fail closed"),
        }
    }

    // CID-created custom flag patterns remain source-empty: the parts
    // vocabulary cannot carry flags, so both flag lists keep their
    // default (empty pattern => has() true) exactly like a
    // default-constructed source Selection.
    #[test]
    fn bio_cid_c01_custom_flags_stay_source_empty() {
        let data = BioSelectionData::from_parts(BioSelectionParts::default())
            .expect("default parts are valid");
        assert!(data.selection.residue_flags.pattern.is_empty());
        assert!(data.selection.atom_flags.pattern.is_empty());
        // A flagged residue still matches the wildcard state, matching
        // the source default FlagList::has (empty pattern = true).
        let residue = SourceResidue {
            name: "ALA",
            seqid_num: 1,
            icode: b' ',
            entity_type: crate::hierarchy::EntityKind::Unknown,
            flag: 0,
            atoms: &[],
        };
        assert!(data.matches_residue(&residue));
        let atom = SourceAtom {
            name: "CA",
            element: 6,
            altloc: 0,
            occ: 1.0,
            b_iso: 20.0,
            flag: 0,
        };
        assert!(data.matches_atom(&atom));
    }
}

#[cfg(test)]
mod tests {
    use super::{SelectionList, is_in_list};

    fn list(all: bool, inverted: bool, list: &str) -> SelectionList {
        SelectionList {
            all,
            inverted,
            list: list.to_string(),
        }
    }

    #[test]
    fn bio_sel_b01_all_matches_every_name_including_empty() {
        // select.hpp:34-35: `if (all) return true;` — before any list use.
        for name in ["", "A", "AB", "water", ","] {
            assert!(
                list(true, false, "").has(name),
                "all + empty list, {name:?}"
            );
            // A stored list is ignored entirely when all is set.
            assert!(list(true, false, "X,Y").has(name), "all + list, {name:?}");
            assert!(list(true, true, "X").has(name), "all + inverted, {name:?}");
        }
    }

    #[test]
    fn bio_sel_b01_empty_list_membership() {
        // is_in_list length shortcut: only a name at least as long as the
        // whole list can match by full equality; "" == "" is true.
        assert!(list(false, false, "").has(""));
        assert!(!list(false, false, "").has("A"));
        // Inversion flips both outcomes.
        assert!(!list(false, true, "").has(""));
        assert!(list(false, true, "").has("A"));
    }

    #[test]
    fn bio_sel_b01_single_and_multi_member_lists() {
        let single = list(false, false, "A");
        assert!(single.has("A"));
        assert!(!single.has("B"));
        assert!(!single.has(""));
        assert!(!single.has("AA"));

        let multi = list(false, false, "A,B");
        assert!(multi.has("A"));
        assert!(multi.has("B"));
        assert!(!multi.has("AB"));
        assert!(!multi.has("A,B,C"));
        // The whole string equals the list, so the length shortcut accepts it.
        assert!(multi.has("A,B"));

        let inverted = list(false, true, "A,B");
        assert!(!inverted.has("A"));
        assert!(!inverted.has("B"));
        assert!(inverted.has("C"));
    }

    #[test]
    fn bio_sel_b01_trailing_and_interior_empty_members() {
        // "A," has members "A" and "" (trailing empty segment).
        let trailing = list(false, false, "A,");
        assert!(trailing.has("A"));
        assert!(trailing.has(""));
        assert!(!trailing.has("B"));
        // ",," has three empty members.
        let empties = list(false, false, ",,");
        assert!(empties.has(""));
        assert!(!empties.has("A"));
        // Interior empty member: "A,,B" matches "".
        let interior = list(false, false, "A,,B");
        assert!(interior.has(""));
        assert!(interior.has("B"));
    }

    #[test]
    fn bio_sel_b01_literal_case_sensitive_comparison() {
        // std::string::compare is byte equality: no case folding.
        let lower = list(false, false, "abc");
        assert!(!lower.has("ABC"));
        assert!(!lower.has("Abc"));
        assert!(lower.has("abc"));
        let mixed = list(false, false, "aB");
        assert!(!mixed.has("Ab"));
        assert!(mixed.has("aB"));
        // Whitespace is significant: " A" != "A".
        let spaced = list(false, false, "A,B");
        assert!(!spaced.has(" A"));
        assert!(!spaced.has("A "));
    }

    #[test]
    fn bio_sel_b01_is_in_list_direct_length_shortcut() {
        // util.hpp:229-230 exercised directly: longer/equal names compare
        // as whole strings only.
        assert!(is_in_list(b"same", b"same"));
        assert!(!is_in_list(b"same1", b"same"));
        assert!(!is_in_list(b"longer name", b"short"));
        // Shorter names scan segments.
        assert!(is_in_list(b"mid", b"start,mid,end"));
        assert!(!is_in_list(b"start,mid", b"start,mid,end"));
    }

    #[test]
    fn bio_sel_b02_empty_pattern_matches_every_flag() {
        // select.hpp:44-45: empty pattern short-circuits to true.
        let empty = super::SelectionFlagList::default();
        for flag in [0u8, b' ', b'A', b'!', b'x'] {
            assert!(empty.has(flag), "flag {flag:?}");
        }
    }

    #[test]
    fn bio_sel_b02_positive_pattern_containment() {
        // find(flag) anywhere in the pattern; zero byte never matches a
        // text pattern; first-byte membership works without special casing.
        let single = super::SelectionFlagList {
            pattern: "A".to_string(),
        };
        assert!(single.has(b'A'));
        assert!(!single.has(b'B'));
        assert!(!single.has(b'!'));
        let multi = super::SelectionFlagList {
            pattern: "ABC".to_string(),
        };
        assert!(multi.has(b'A'));
        assert!(multi.has(b'B'));
        assert!(multi.has(b'C'));
        assert!(!multi.has(b'a'));
        assert!(!multi.has(0));
        // A '!' in a non-leading position is an ordinary member.
        let interior = super::SelectionFlagList {
            pattern: "X!".to_string(),
        };
        assert!(interior.has(b'!'));
        assert!(interior.has(b'X'));
        assert!(!interior.has(b'Y'));
    }

    #[test]
    fn bio_sel_b02_inverted_pattern_skips_leading_bang() {
        // select.hpp:47-48: search starts at index 1 after '!', and the
        // result is inverted; '!' itself is therefore never excluded.
        let inverted = super::SelectionFlagList {
            pattern: "!A".to_string(),
        };
        assert!(!inverted.has(b'A'));
        assert!(inverted.has(b'B'));
        assert!(inverted.has(b'!'));
        assert!(inverted.has(0));
        let inverted_multi = super::SelectionFlagList {
            pattern: "!AB".to_string(),
        };
        assert!(!inverted_multi.has(b'A'));
        assert!(!inverted_multi.has(b'B'));
        assert!(inverted_multi.has(b'C'));
        // A bare "!" pattern searches an empty member set and inverts it:
        // every flag matches.
        let bare = super::SelectionFlagList {
            pattern: "!".to_string(),
        };
        assert!(bare.has(b'A'));
        assert!(bare.has(b'!'));
    }

    #[test]
    fn bio_sel_b03_empty_sentinels_and_representable_numbers() {
        use super::SelectionSequenceId as Sid;
        // select.hpp:57-59: both INT_MIN and INT_MAX are "empty".
        assert!(
            Sid {
                seqnum: i32::MIN,
                icode: b'*'
            }
            .empty()
        );
        assert!(
            Sid {
                seqnum: i32::MAX,
                icode: b' '
            }
            .empty()
        );
        for seqnum in [-1, 0, 1, 42, i32::MIN + 1, i32::MAX - 1] {
            assert!(
                !Sid {
                    seqnum,
                    icode: b' '
                }
                .empty(),
                "{seqnum}"
            );
        }
    }

    #[test]
    fn bio_sel_b03_compare_orders_numbers_strictly() {
        use super::SelectionSequenceId as Sid;
        // Number comparison dominates; icode is irrelevant then.
        assert_eq!(
            Sid {
                seqnum: 5,
                icode: b'Z'
            }
            .compare(10, b' '),
            -1
        );
        assert_eq!(
            Sid {
                seqnum: 10,
                icode: b'A'
            }
            .compare(5, b' '),
            1
        );
        assert_eq!(
            Sid {
                seqnum: -3,
                icode: b' '
            }
            .compare(-3, b' '),
            0
        );
        assert_eq!(
            Sid {
                seqnum: 0,
                icode: b'*'
            }
            .compare(0, b'A'),
            0
        );
        // INT_MIN/INT_MAX are ordinary extremes in compare (no special
        // casing): INT_MIN sorts below everything.
        assert_eq!(
            Sid {
                seqnum: i32::MIN,
                icode: b' '
            }
            .compare(i32::MIN + 1, b' '),
            -1
        );
        assert_eq!(
            Sid {
                seqnum: i32::MAX,
                icode: b' '
            }
            .compare(i32::MAX - 1, b' '),
            1
        );
    }

    #[test]
    fn bio_sel_b03_compare_insertion_codes_and_wildcard() {
        use super::SelectionSequenceId as Sid;
        // Equal numbers: icode decides; wildcard b'*' equals everything.
        assert_eq!(
            Sid {
                seqnum: 7,
                icode: b'*'
            }
            .compare(7, b' '),
            0
        );
        assert_eq!(
            Sid {
                seqnum: 7,
                icode: b'*'
            }
            .compare(7, b'Z'),
            0
        );
        // Blank (0x20) sorts before letters.
        assert_eq!(
            Sid {
                seqnum: 7,
                icode: b' '
            }
            .compare(7, b'A'),
            -1
        );
        assert_eq!(
            Sid {
                seqnum: 7,
                icode: b'A'
            }
            .compare(7, b' '),
            1
        );
        assert_eq!(
            Sid {
                seqnum: 7,
                icode: b'A'
            }
            .compare(7, b'B'),
            -1
        );
        assert_eq!(
            Sid {
                seqnum: 7,
                icode: b'B'
            }
            .compare(7, b'A'),
            1
        );
        assert_eq!(
            Sid {
                seqnum: 7,
                icode: b'A'
            }
            .compare(7, b'A'),
            0
        );
    }
    #[test]
    fn bio_sel_b04_q_b_scalar_comparisons() {
        use super::SelectionAtomInequality as Ineq;
        // q reads occupancy, b reads B-factor; f32 carriers promote to
        // f64 exactly, so equality against a representable threshold holds.
        let occ = 0.5f32;
        let b_iso = 20.0f32;
        for (relation, expected_q, expected_b) in
            [(-1, false, false), (0, true, false), (1, false, true)]
        {
            assert_eq!(
                Ineq {
                    property: b'q',
                    relation,
                    value: 0.5
                }
                .matches(occ, b_iso),
                expected_q,
                "q relation {relation}"
            );
            assert_eq!(
                Ineq {
                    property: b'b',
                    relation,
                    value: 0.5
                }
                .matches(occ, b_iso),
                expected_b,
                "b relation {relation}"
            );
        }
        // Threshold crossing both directions on each property.
        assert!(
            Ineq {
                property: b'q',
                relation: -1,
                value: 0.75
            }
            .matches(occ, b_iso)
        );
        assert!(
            Ineq {
                property: b'q',
                relation: 1,
                value: 0.25
            }
            .matches(occ, b_iso)
        );
        assert!(
            Ineq {
                property: b'b',
                relation: 1,
                value: 19.5
            }
            .matches(occ, b_iso)
        );
        assert!(
            Ineq {
                property: b'b',
                relation: -1,
                value: 20.5
            }
            .matches(occ, b_iso)
        );
        // Negative and zero thresholds behave through the same branches.
        assert!(
            Ineq {
                property: b'b',
                relation: 1,
                value: -1.0
            }
            .matches(occ, b_iso)
        );
        assert!(
            Ineq {
                property: b'b',
                relation: 0,
                value: 0.0
            }
            .matches(0.0f32, 0.0f32)
        );
    }

    #[test]
    fn bio_sel_b04_unknown_property_uses_initial_zero() {
        use super::SelectionAtomInequality as Ineq;
        // select.hpp:75-79: neither 'q' nor 'b' leaves atom_value at 0.0;
        // case variants are distinct bytes and are unknown too.
        for property in [b'x', b'Q', b'B', 0, b' '] {
            // occupancy/B-factor values are ignored entirely.
            assert!(
                Ineq {
                    property,
                    relation: -1,
                    value: 1.0
                }
                .matches(99.0f32, 99.0f32),
                "property {property:?}"
            );
            assert!(
                Ineq {
                    property,
                    relation: 0,
                    value: 0.0
                }
                .matches(99.0f32, 99.0f32)
            );
            assert!(
                !Ineq {
                    property,
                    relation: 1,
                    value: 0.0
                }
                .matches(99.0f32, 99.0f32)
            );
            // -0.0 == 0.0 under IEEE equality.
            assert!(
                Ineq {
                    property,
                    relation: 0,
                    value: -0.0
                }
                .matches(99.0f32, 99.0f32)
            );
        }
    }

    #[test]
    fn bio_sel_b04_nan_and_signed_zero_semantics() {
        use super::SelectionAtomInequality as Ineq;
        let nan = f32::NAN;
        let value = 1.0f64;
        // NaN atom values satisfy no relation.
        assert!(
            !Ineq {
                property: b'q',
                relation: -1,
                value
            }
            .matches(nan, 0.0f32)
        );
        assert!(
            !Ineq {
                property: b'q',
                relation: 0,
                value
            }
            .matches(nan, 0.0f32)
        );
        assert!(
            !Ineq {
                property: b'q',
                relation: 1,
                value
            }
            .matches(nan, 0.0f32)
        );
        // A NaN threshold also never matches a finite carrier.
        assert!(
            !Ineq {
                property: b'b',
                relation: 0,
                value: f64::NAN
            }
            .matches(1.0f32, 20.0f32)
        );
        // Signed zero carriers compare equal to both zero spellings.
        assert!(
            Ineq {
                property: b'q',
                relation: 0,
                value: 0.0
            }
            .matches(-0.0f32, 0.0f32)
        );
        // -0.0 < 0.0 is false.
        assert!(
            !Ineq {
                property: b'q',
                relation: -1,
                value: 0.0
            }
            .matches(-0.0f32, 0.0f32)
        );
        // f32 carrier narrowing is exact: 0.1f32 promotes to a double
        // unequal to the decimal literal 0.1.
        assert_ne!(f64::from(0.1f32), 0.1f64);
        assert!(
            !Ineq {
                property: b'q',
                relation: 0,
                value: 0.1
            }
            .matches(0.1f32, 0.0f32)
        );
        assert!(
            Ineq {
                property: b'q',
                relation: 0,
                value: f64::from(0.1f32)
            }
            .matches(0.1f32, 0.0f32)
        );
    }
    // (BIO-CID C02) The lexical regressions bio_sel_b05 (wrong_syntax),
    // b06 (determine_omitted_cid_fields), b07 (make_cid_list), b08
    // (parse_cid_seqid), b09 (parse_cid_elements) and b10
    // (has_inequality) were relocated UNCHANGED into
    // `cosmolkit-io/src/bio_cid.rs` (same names/inputs/assertions) when
    // their owning functions moved. Predicate/carrier tests b01-b04 and
    // the Selection tests below stay here.
    #[test]
    fn bio_sel_b11_wildcard_and_exact_model_identity() {
        use super::{Selection, SourceModel};
        let no_chains: [super::SourceChain; 0] = [];
        let model = |num: i32| SourceModel {
            num,
            chains: &no_chains,
        };
        // Default Selection has mdl = 0 -> all models.
        let wildcard = Selection::default();
        for num in [-1, 0, 1, 2, 9999] {
            assert!(wildcard.matches_model(&model(num)), "num={num}");
        }
        // Exact identity: only the equal source number matches.
        let mut sel = Selection::default();
        sel.mdl = 1;
        assert!(sel.matches_model(&model(1)));
        for num in [-1, 0, 2, 9999] {
            assert!(!sel.matches_model(&model(num)), "num={num}");
        }
        // Negative model numbers are ordinary source identities.
        sel.mdl = -3;
        assert!(sel.matches_model(&model(-3)));
        assert!(!sel.matches_model(&model(3)));
        assert!(!sel.matches_model(&model(0)));
        // An explicit mdl = 0 stays wildcard even for num = 0.
        sel.mdl = 0;
        assert!(sel.matches_model(&model(0)));
        // i32 extremes are exact identities.
        sel.mdl = i32::MIN;
        assert!(sel.matches_model(&model(i32::MIN)));
        assert!(!sel.matches_model(&model(0)));
    }
    #[test]
    fn bio_sel_b12_chain_name_membership() {
        use super::{Selection, SourceChain};
        fn mk_chain<'a>(
            name: &'a str,
            residues: &'a [super::SourceResidue<'a>],
        ) -> SourceChain<'a> {
            SourceChain { name, residues }
        }
        let no_residues: [super::SourceResidue; 0] = [];
        // Direct helper calls keep each name's lifetime independent of
        // the residue-array borrow (no unifying closure).
        // Default wildcard list matches every name, including blank and
        // multi-character names.
        let wildcard = Selection::default();
        for name in ["", "A", "AB", "SEGA", "z"] {
            assert!(
                wildcard.matches_chain(&mk_chain(name, &no_residues)),
                "name={name:?}"
            );
        }
        // Positive comma list: literal byte membership, case-sensitive.
        let mut sel = Selection::default();
        sel.chain_ids = super::SelectionList {
            all: false,
            inverted: false,
            list: "A,B".to_string(),
        };
        assert!(sel.matches_chain(&mk_chain("A", &no_residues)));
        assert!(sel.matches_chain(&mk_chain("B", &no_residues)));
        assert!(!sel.matches_chain(&mk_chain("a", &no_residues)));
        assert!(!sel.matches_chain(&mk_chain("AB", &no_residues)));
        assert!(!sel.matches_chain(&mk_chain("", &no_residues)));
        // Inverted list: everything except members, blank included.
        sel.chain_ids.inverted = true;
        assert!(!sel.matches_chain(&mk_chain("A", &no_residues)));
        assert!(sel.matches_chain(&mk_chain("C", &no_residues)));
        assert!(sel.matches_chain(&mk_chain("", &no_residues)));
        // Duplicate source names behave identically per row (the predicate
        // is stateless); order and duplication are the traversal's concern.
        // The list is restored to non-inverted first.
        sel.chain_ids.inverted = false;
        let chains = [
            mk_chain("A", &no_residues),
            mk_chain("A", &no_residues),
            mk_chain("B", &no_residues),
        ];
        let matched: Vec<&str> = chains
            .iter()
            .filter(|chain| sel.matches_chain(chain))
            .map(|chain| chain.name)
            .collect();
        assert_eq!(matched, vec!["A", "A", "B"]);
    }
    #[test]
    fn bio_sel_b13_entity_type_gate_and_flags() {
        use super::Selection;
        use crate::hierarchy::EntityKind;
        let no_atoms: [super::SourceAtom; 0] = [];
        let res = |entity_type, name, num, icode, flag| super::SourceResidue {
            entity_type,
            name,
            seqid_num: num,
            icode,
            flag,
            atoms: &no_atoms,
        };
        // Default selection matches every residue of every entity kind.
        let wildcard = Selection::default();
        for kind in [
            EntityKind::Unknown,
            EntityKind::Polymer,
            EntityKind::NonPolymer,
            EntityKind::Branched,
            EntityKind::Water,
        ] {
            assert!(
                wildcard.matches_residue(&res(kind, "ALA", 1, b' ', 0)),
                "{kind:?}"
            );
        }
        // et_flags ordinal selection: only the flagged kind passes; the
        // `all` list stays true but et_flags is consulted only when the
        // entity_types list is non-all (source short-circuit).
        let mut sel = Selection::default();
        sel.entity_types = super::SelectionList {
            all: false,
            inverted: false,
            list: String::new(),
        };
        sel.et_flags = [false; 6];
        sel.et_flags[3] = true; // Branched
        assert!(sel.matches_residue(&res(EntityKind::Branched, "X", 1, b' ', 0)));
        assert!(!sel.matches_residue(&res(EntityKind::Polymer, "X", 1, b' ', 0)));
        assert!(!sel.matches_residue(&res(EntityKind::Water, "X", 1, b' ', 0)));
        // Explicit scalar flag: pattern selection without substitution.
        sel = Selection::default();
        sel.residue_flags.pattern = "s".to_string();
        assert!(sel.matches_residue(&res(EntityKind::Polymer, "ALA", 1, b' ', b's')));
        assert!(!sel.matches_residue(&res(EntityKind::Polymer, "ALA", 1, b' ', 0)));
        assert!(!sel.matches_residue(&res(EntityKind::Polymer, "ALA", 1, b' ', b'x')));
    }

    #[test]
    fn bio_sel_b13_name_and_seqid_bounds() {
        use super::Selection;
        use crate::hierarchy::EntityKind;
        let no_atoms: [super::SourceAtom; 0] = [];
        let res = |name, num, icode| super::SourceResidue {
            entity_type: EntityKind::Polymer,
            name,
            seqid_num: num,
            icode,
            flag: 0,
            atoms: &no_atoms,
        };
        // Name membership is literal and case-sensitive.
        let mut sel = Selection::default();
        sel.residue_names = super::SelectionList {
            all: false,
            inverted: false,
            list: "ALA,GLY".to_string(),
        };
        assert!(sel.matches_residue(&res("ALA", 5, b' ')));
        assert!(sel.matches_residue(&res("GLY", 5, b' ')));
        assert!(!sel.matches_residue(&res("ala", 5, b' ')));
        assert!(!sel.matches_residue(&res("ALAA", 5, b' ')));
        // Inclusive bounds: from 3/A to 7/B.
        sel = Selection::default();
        sel.from_seqid = super::SelectionSequenceId {
            seqnum: 3,
            icode: b'A',
        };
        sel.to_seqid = super::SelectionSequenceId {
            seqnum: 7,
            icode: b'B',
        };
        // In-range interiors and exact endpoints.
        assert!(sel.matches_residue(&res("X", 4, b' ')));
        assert!(sel.matches_residue(&res("X", 7, b'A')));
        assert!(sel.matches_residue(&res("X", 7, b'B')));
        // Below the lower endpoint: 3 itself only at icode >= 'A'.
        assert!(sel.matches_residue(&res("X", 3, b'A')));
        assert!(!sel.matches_residue(&res("X", 3, b' ')));
        assert!(!sel.matches_residue(&res("X", 2, b'Z')));
        // Above the upper endpoint: 7 only up to icode <= 'B'.
        assert!(!sel.matches_residue(&res("X", 7, b'C')));
        assert!(!sel.matches_residue(&res("X", 8, b' ')));
        // Wildcard selector icode on the endpoints spans all icodes.
        sel.from_seqid.icode = b'*';
        sel.to_seqid.icode = b'*';
        assert!(sel.matches_residue(&res("X", 3, b' ')));
        assert!(sel.matches_residue(&res("X", 7, b'Z')));
        assert!(!sel.matches_residue(&res("X", 2, b' ')));
        assert!(!sel.matches_residue(&res("X", 8, b' ')));
        // All five conjuncts must hold together: name AND range.
        sel.residue_names = super::SelectionList {
            all: false,
            inverted: false,
            list: "GLY".to_string(),
        };
        assert!(!sel.matches_residue(&res("ALA", 5, b' ')));
        assert!(sel.matches_residue(&res("GLY", 5, b' ')));
    }
    #[test]
    fn bio_sel_b14_name_element_altloc_flag() {
        use super::{Selection, SourceAtom};
        fn mk_atom<'a>(name: &'a str, element: u8, altloc: u8, flag: u8) -> SourceAtom<'a> {
            SourceAtom {
                name,
                element,
                altloc,
                flag,
                occ: 1.0,
                b_iso: 20.0,
            }
        }
        // Default selection matches any atom.
        assert!(Selection::default().matches_atom(&mk_atom("CA", 6, 0, 0)));
        // Literal name membership.
        let mut sel = Selection::default();
        sel.atom_names = super::SelectionList {
            all: false,
            inverted: false,
            list: "CA,N".to_string(),
        };
        assert!(sel.matches_atom(&mk_atom("CA", 6, 0, 0)));
        assert!(sel.matches_atom(&mk_atom("N", 7, 0, 0)));
        assert!(!sel.matches_atom(&mk_atom("ca", 6, 0, 0)));
        // Element mask on Gemmi ordinals: carbon=6 selected, nitrogen=7 not.
        sel = Selection::default();
        let mut mask = [false; super::GEMMI_EL_END];
        mask[6] = true;
        sel.elements = Some(mask);
        assert!(sel.matches_atom(&mk_atom("CA", 6, 0, 0)));
        assert!(!sel.matches_atom(&mk_atom("N", 7, 0, 0)));
        assert!(sel.matches_atom(&mk_atom("X", 0, 0, 0)) == false);
        // Altloc: 0 -> empty name; 'A' -> one-byte name; no normalization.
        sel = Selection::default();
        sel.altlocs = super::SelectionList {
            all: false,
            inverted: false,
            list: "A".to_string(),
        };
        assert!(sel.matches_atom(&mk_atom("CA", 6, b'A', 0)));
        assert!(!sel.matches_atom(&mk_atom("CA", 6, b'B', 0)));
        assert!(!sel.matches_atom(&mk_atom("CA", 6, 0, 0)));
        // An explicit empty altloc member selects unset altlocs (the source
        // compares against the empty string).
        sel.altlocs.list = ",A".to_string();
        assert!(sel.matches_atom(&mk_atom("CA", 6, 0, 0)));
        assert!(sel.matches_atom(&mk_atom("CA", 6, b'A', 0)));
        assert!(!sel.matches_atom(&mk_atom("CA", 6, b'B', 0)));
        // Inverted altloc list.
        sel.altlocs.inverted = true;
        assert!(!sel.matches_atom(&mk_atom("CA", 6, b'A', 0)));
        assert!(sel.matches_atom(&mk_atom("CA", 6, b'B', 0)));
        // Explicit scalar atom flag pattern.
        sel = Selection::default();
        sel.atom_flags.pattern = "d".to_string();
        assert!(sel.matches_atom(&mk_atom("CA", 6, 0, b'd')));
        assert!(!sel.matches_atom(&mk_atom("CA", 6, 0, 0)));
    }

    #[test]
    fn bio_sel_b14_inequalities_on_f32_carriers() {
        use super::{Selection, SelectionAtomInequality, SourceAtom};
        fn mk_qb_atom(occ: f32, b_iso: f32) -> SourceAtom<'static> {
            SourceAtom {
                name: "CA",
                element: 6,
                altloc: 0,
                flag: 0,
                occ,
                b_iso,
            }
        }
        let ineq = |property: u8, relation: i32, value: f64| SelectionAtomInequality {
            property,
            relation,
            value,
        };
        // Zero inequalities: vacuously true (source all_of on empty).
        let mut sel = Selection::default();
        assert!(sel.matches_atom(&mk_qb_atom(0.5, 30.0)));
        // Single q inequality.
        sel.atom_inequalities = vec![ineq(b'q', -1, 0.5)];
        assert!(sel.matches_atom(&mk_qb_atom(0.25, 30.0)));
        assert!(!sel.matches_atom(&mk_qb_atom(0.75, 30.0)));
        assert!(!sel.matches_atom(&mk_qb_atom(0.5, 30.0))); // strict < is false for equal
        sel.atom_inequalities = vec![ineq(b'q', 0, 0.5)];
        assert!(sel.matches_atom(&mk_qb_atom(0.5, 30.0)));
        sel.atom_inequalities = vec![ineq(b'q', 1, 0.5)];
        assert!(sel.matches_atom(&mk_qb_atom(0.75, 30.0)));
        assert!(!sel.matches_atom(&mk_qb_atom(0.25, 30.0)));
        // b inequality on the f32-narrowed carrier.
        sel.atom_inequalities = vec![ineq(b'b', 1, 20.0)];
        assert!(sel.matches_atom(&mk_qb_atom(0.5, 25.5)));
        assert!(!sel.matches_atom(&mk_qb_atom(0.5, 15.0)));
        // Multiple inequalities must all hold.
        sel.atom_inequalities = vec![ineq(b'q', 1, 0.3), ineq(b'b', -1, 40.0)];
        assert!(sel.matches_atom(&mk_qb_atom(0.5, 35.0)));
        assert!(!sel.matches_atom(&mk_qb_atom(0.2, 35.0)));
        assert!(!sel.matches_atom(&mk_qb_atom(0.5, 45.0)));
        // Unknown property leaves atom_value at 0.0 (source branch).
        sel.atom_inequalities = vec![ineq(b'z', 0, 0.0)];
        assert!(sel.matches_atom(&mk_qb_atom(0.5, 30.0)));
        sel.atom_inequalities = vec![ineq(b'z', 1, 0.0)];
        assert!(!sel.matches_atom(&mk_qb_atom(0.5, 30.0)));
        // NaN carrier satisfies no relation.
        sel.atom_inequalities = vec![ineq(b'q', 0, f64::NAN)];
        assert!(!sel.matches_atom(&mk_qb_atom(f32::NAN, 30.0)));
    }
    #[test]
    fn bio_sel_b15_ancestor_presence_and_mismatch() {
        use super::{Selection, SourceAtom, SourceChain, SourceCra, SourceResidue};
        use crate::hierarchy::EntityKind;
        let atoms = [SourceAtom {
            name: "CA",
            element: 6,
            altloc: 0,
            flag: 0,
            occ: 1.0,
            b_iso: 20.0,
        }];
        let residues = [SourceResidue {
            entity_type: EntityKind::Polymer,
            name: "GLY",
            seqid_num: 5,
            icode: b' ',
            flag: 0,
            atoms: &atoms,
        }];
        let chains = [SourceChain {
            name: "A",
            residues: &residues,
        }];
        let chain = &chains[0];
        let residue = &residues[0];
        let atom = &atoms[0];
        // A selection that rejects chain name A, residue name GLY and atom
        // name CA individually, built per case below.
        //
        // All eight presence combinations with a wildcard selection:
        // every combination matches.
        let wildcard = Selection::default();
        for (chain, residue, atom) in [
            (None, None, None),
            (Some(chain), None, None),
            (None, Some(residue), None),
            (None, None, Some(atom)),
            (Some(chain), Some(residue), None),
            (Some(chain), None, Some(atom)),
            (None, Some(residue), Some(atom)),
            (Some(chain), Some(residue), Some(atom)),
        ] {
            let cra = SourceCra {
                chain,
                residue,
                atom,
            };
            assert!(wildcard.matches_cra(&cra));
        }
        // Mismatch at each present level vetoes the whole conjunction.
        let mut sel = Selection::default();
        sel.chain_ids = super::SelectionList {
            all: false,
            inverted: false,
            list: "B".to_string(),
        };
        assert!(!sel.matches_cra(&SourceCra {
            chain: Some(chain),
            residue: None,
            atom: None,
        }));
        // Absent chain conjunct passes; the residue name still decides.
        sel = Selection::default();
        sel.residue_names = super::SelectionList {
            all: false,
            inverted: false,
            list: "ALA".to_string(),
        };
        assert!(!sel.matches_cra(&SourceCra {
            chain: None,
            residue: Some(residue),
            atom: None,
        }));
        sel = Selection::default();
        sel.atom_names = super::SelectionList {
            all: false,
            inverted: false,
            list: "N".to_string(),
        };
        assert!(!sel.matches_cra(&SourceCra {
            chain: None,
            residue: None,
            atom: Some(atom),
        }));
        // Match at the present level passes.
        sel.atom_names.list = "CA".to_string();
        assert!(sel.matches_cra(&SourceCra {
            chain: None,
            residue: None,
            atom: Some(atom),
        }));
    }
    #[test]
    fn bio_sel_b16_ordered_lazy_filtering() {
        use super::{Selection, SourceAtom, SourceChain, SourceModel, SourceResidue};
        use crate::hierarchy::EntityKind;
        let atoms = [
            SourceAtom {
                name: "CA",
                element: 6,
                altloc: 0,
                flag: 0,
                occ: 1.0,
                b_iso: 20.0,
            },
            SourceAtom {
                name: "N",
                element: 7,
                altloc: b'A',
                flag: 0,
                occ: 1.0,
                b_iso: 20.0,
            },
            SourceAtom {
                name: "CA",
                element: 6,
                altloc: b'B',
                flag: 0,
                occ: 0.5,
                b_iso: 35.0,
            },
        ];
        let residues = [
            SourceResidue {
                entity_type: EntityKind::Polymer,
                name: "GLY",
                seqid_num: 1,
                icode: b' ',
                flag: 0,
                atoms: &atoms,
            },
            SourceResidue {
                entity_type: EntityKind::Water,
                name: "HOH",
                seqid_num: 5,
                icode: b' ',
                flag: 0,
                atoms: &[],
            },
            SourceResidue {
                entity_type: EntityKind::Polymer,
                name: "GLY",
                seqid_num: 2,
                icode: b' ',
                flag: 0,
                atoms: &atoms,
            },
        ];
        let chains = [
            SourceChain {
                name: "A",
                residues: &residues,
            },
            SourceChain {
                name: "B",
                residues: &[],
            },
            SourceChain {
                name: "A",
                residues: &residues,
            },
        ];
        let models = [
            SourceModel {
                num: 1,
                chains: &chains,
            },
            SourceModel {
                num: 2,
                chains: &[],
            },
            SourceModel {
                num: 1,
                chains: &chains,
            },
        ];
        // Empty input: no iterations, no panic.
        let empty: [SourceModel; 0] = [];
        assert_eq!(Selection::default().models(&empty).count(), 0);
        // Exact retained source order with duplicate names, no synthesis
        // or sorting: model filter mdl=1 keeps 1,1 in positions 0,2.
        let mut sel = Selection::default();
        sel.mdl = 1;
        let kept: Vec<i32> = sel.models(&models).map(|model| model.num).collect();
        assert_eq!(kept, vec![1, 1]);
        // Chain-name filter over a mixed hierarchy retains order and
        // duplicates: A, A from positions 0, 2.
        sel = Selection::default();
        sel.chain_ids = super::SelectionList {
            all: false,
            inverted: false,
            list: "A".to_string(),
        };
        let kept: Vec<&str> = sel.chains(&chains).map(|chain| chain.name).collect();
        assert_eq!(kept, vec!["A", "A"]);
        // Residue filter (polymer GLY) keeps positions 0 and 2 in order.
        sel = Selection::default();
        sel.residue_names = super::SelectionList {
            all: false,
            inverted: false,
            list: "GLY".to_string(),
        };
        let kept: Vec<i32> = sel
            .residues(&residues)
            .map(|residue| residue.seqid_num)
            .collect();
        assert_eq!(kept, vec![1, 2]);
        // Atom filter (name CA) keeps positions 0 and 2 in order.
        sel = Selection::default();
        sel.atom_names = super::SelectionList {
            all: false,
            inverted: false,
            list: "CA".to_string(),
        };
        let kept: Vec<&str> = sel.atoms(&atoms).map(|atom| atom.name).collect();
        assert_eq!(kept, vec!["CA", "CA"]);
        // Laziness: an iterator over an empty selection still yields items
        // without pre-materializing; count() drives one predicate per row.
        assert_eq!(Selection::default().atoms(&atoms).count(), 3);
        // A rejected-parent-shaped composition: iterating a chain's atoms
        // through an empty-predicate filter is independent of whether any
        // chain matched; order is the row order inside the slice.
        sel = Selection::default();
        sel.chain_ids = super::SelectionList {
            all: false,
            inverted: false,
            list: "Z".to_string(),
        };
        assert_eq!(sel.chains(&chains).count(), 0);
        assert_eq!(sel.atoms(&atoms).count(), 3);
    }
    #[test]
    fn bio_sel_b17_first_match_traversal() {
        use super::{Selection, SourceAtom, SourceChain, SourceModel, SourceResidue};
        use crate::hierarchy::EntityKind;
        let atoms = [
            SourceAtom {
                name: "N",
                element: 7,
                altloc: 0,
                flag: 0,
                occ: 1.0,
                b_iso: 20.0,
            },
            SourceAtom {
                name: "CA",
                element: 6,
                altloc: 0,
                flag: 0,
                occ: 1.0,
                b_iso: 20.0,
            },
        ];
        let residues = [
            SourceResidue {
                entity_type: EntityKind::Polymer,
                name: "GLY",
                seqid_num: 1,
                icode: b' ',
                flag: 0,
                atoms: &atoms,
            },
            SourceResidue {
                entity_type: EntityKind::Polymer,
                name: "ALA",
                seqid_num: 2,
                icode: b' ',
                flag: 0,
                atoms: &atoms,
            },
        ];
        let chains = [
            SourceChain {
                name: "A",
                residues: &residues,
            },
            SourceChain {
                name: "B",
                residues: &residues,
            },
        ];
        let model = SourceModel {
            num: 1,
            chains: &chains,
        };
        // Empty model and empty slices.
        let empty_model = SourceModel {
            num: 1,
            chains: &[],
        };
        let none = Selection::default().first_in_model(&empty_model);
        assert!(none.chain.is_none() && none.residue.is_none() && none.atom.is_none());
        // First match overall: chain A, residue GLY, atom N (source order).
        let cra = Selection::default().first_in_model(&model);
        assert_eq!(cra.chain.unwrap().name, "A");
        assert_eq!(cra.residue.unwrap().name, "GLY");
        assert_eq!(cra.atom.unwrap().name, "N");
        // Middle match: skip chain A and residue GLY.
        let mut sel = Selection::default();
        sel.chain_ids = super::SelectionList {
            all: false,
            inverted: false,
            list: "B".to_string(),
        };
        sel.residue_names = super::SelectionList {
            all: false,
            inverted: false,
            list: "ALA".to_string(),
        };
        let cra = sel.first_in_model(&model);
        assert_eq!(cra.chain.unwrap().name, "B");
        assert_eq!(cra.residue.unwrap().seqid_num, 2);
        assert_eq!(cra.atom.unwrap().name, "N");
        // Last match: second atom.
        sel = Selection::default();
        sel.atom_names = super::SelectionList {
            all: false,
            inverted: false,
            list: "CA".to_string(),
        };
        let cra = sel.first_in_model(&model);
        assert_eq!(cra.chain.unwrap().name, "A");
        assert_eq!(cra.residue.unwrap().name, "GLY");
        assert_eq!(cra.atom.unwrap().name, "CA");
        // No match: rejected model gate.
        sel = Selection::default();
        sel.mdl = 2;
        let cra = sel.first_in_model(&model);
        assert!(cra.chain.is_none());
        // No match: rejected chain gate.
        sel = Selection::default();
        sel.chain_ids = super::SelectionList {
            all: false,
            inverted: false,
            list: "Z".to_string(),
        };
        let cra = sel.first_in_model(&model);
        assert!(cra.chain.is_none());
        // A matching residue with no matching atom does not return a
        // result: the scan continues (here into the next residue).
        let mixed_residues = [
            SourceResidue {
                entity_type: EntityKind::Polymer,
                name: "GLY",
                seqid_num: 1,
                icode: b' ',
                flag: 0,
                atoms: &atoms[..1],
            },
            SourceResidue {
                entity_type: EntityKind::Polymer,
                name: "ALA",
                seqid_num: 2,
                icode: b' ',
                flag: 0,
                atoms: &atoms,
            },
        ];
        let mixed_chains = [SourceChain {
            name: "A",
            residues: &mixed_residues,
        }];
        let mixed_model = SourceModel {
            num: 1,
            chains: &mixed_chains,
        };
        sel = Selection::default();
        sel.atom_names = super::SelectionList {
            all: false,
            inverted: false,
            list: "CA".to_string(),
        };
        let cra = sel.first_in_model(&mixed_model);
        assert_eq!(cra.residue.unwrap().name, "ALA");
        assert_eq!(cra.atom.unwrap().name, "CA");
        // Duplicate chain names: the first duplicate in source order wins.
        let dup_chains = [
            SourceChain {
                name: "A",
                residues: &residues,
            },
            SourceChain {
                name: "A",
                residues: &residues,
            },
        ];
        let dup_model = SourceModel {
            num: 1,
            chains: &dup_chains,
        };
        let cra = Selection::default().first_in_model(&dup_model);
        let first = cra.chain.unwrap() as *const _;
        let second = &dup_chains[1] as *const _;
        assert_ne!(first, second);
        assert_eq!(cra.chain.unwrap().name, "A");
    }
    #[test]
    fn bio_sel_b18_model_ordered_first_match() {
        use super::{Selection, SourceAtom, SourceChain, SourceModel, SourceResidue};
        use crate::hierarchy::EntityKind;
        let atoms = [SourceAtom {
            name: "CA",
            element: 6,
            altloc: 0,
            flag: 0,
            occ: 1.0,
            b_iso: 20.0,
        }];
        let residues = [SourceResidue {
            entity_type: EntityKind::Polymer,
            name: "GLY",
            seqid_num: 1,
            icode: b' ',
            flag: 0,
            atoms: &atoms,
        }];
        let chains = [SourceChain {
            name: "A",
            residues: &residues,
        }];
        // Empty models.
        let empty: [SourceModel; 0] = [];
        assert!(Selection::default().first(&empty).is_none());
        // Noncontiguous model numbers: stored order decides, not number
        // order (5 before 3 in the slice).
        let model_five = SourceModel {
            num: 5,
            chains: &chains,
        };
        let model_three = SourceModel {
            num: 3,
            chains: &chains,
        };
        // Build a slice of the two models for first().
        let model_slice = [model_five, model_three];
        let (model, cra) = Selection::default().first(&model_slice).unwrap();
        assert_eq!(model.num, 5);
        assert_eq!(cra.chain.unwrap().name, "A");
        // Earlier rejected, later matching: model 1 is rejected by mdl,
        // model 7 matches.
        let model_one = SourceModel {
            num: 1,
            chains: &chains,
        };
        let model_seven = SourceModel {
            num: 7,
            chains: &chains,
        };
        let models = [model_one, model_seven];
        let mut sel = Selection::default();
        sel.mdl = 7;
        let (model, cra) = sel.first(&models).unwrap();
        assert_eq!(model.num, 7);
        assert_eq!(cra.atom.unwrap().name, "CA");
        // Earlier empty, later populated: an empty model contributes no
        // hit even when it is not filtered out.
        let empty_chains: [SourceChain; 0] = [];
        let hollow = SourceModel {
            num: 2,
            chains: &empty_chains,
        };
        let models = [hollow, model_seven];
        let (model, _) = Selection::default().first(&models).unwrap();
        assert_eq!(model.num, 7);
        // No model matches: None.
        let models = [model_one, model_three];
        let mut sel = Selection::default();
        sel.mdl = 9;
        assert!(sel.first(&models).is_none());
    }
    // (BIO-CID C02) bio_sel_b05_source_message_bytes_ascii_and_offsets and
    // bio_sel_b05_source_message_bytes_multibyte_window_cuts (BIO-SEL-FIX
    // Step 12 byte-exact message regressions) were relocated UNCHANGED to
    // `cosmolkit-io/src/bio_cid.rs` with the `wrong_syntax` owner.
    #[test]
    fn bio_rows_r03_model_numbers_and_missing_representation() {
        use super::{BioRowModelError, Selection, bio_model_row_matches};
        use crate::hierarchy::{BioChainId, BioModelRow, BioRowSpan};

        fn row(span: BioRowSpan<BioChainId>, num: Option<i32>) -> BioModelRow {
            BioModelRow::new(span, num)
        }
        let span = BioRowSpan::<BioChainId>::new(0, 0).unwrap();
        let zero = BioRowSpan::<BioChainId>::new(0, 0).unwrap();

        // Some(negative/zero/noncontiguous positive): the row predicate
        // compares the exact reader-assigned int — no renumbering.
        let mut sel = Selection::default();
        sel.mdl = -3;
        assert_eq!(bio_model_row_matches(&sel, &row(span, Some(-3))), Ok(true));
        assert_eq!(bio_model_row_matches(&sel, &row(zero, Some(-2))), Ok(false));
        sel.mdl = 0;
        assert_eq!(bio_model_row_matches(&sel, &row(zero, Some(0))), Ok(true));
        sel.mdl = 7;
        assert_eq!(bio_model_row_matches(&sel, &row(zero, Some(7))), Ok(true));
        // Noncontiguous source numbering (2 then 9) matches by value, not
        // by order or index.
        let m2 = BioRowSpan::<BioChainId>::new(0, 1).unwrap();
        let m9 = BioRowSpan::<BioChainId>::new(1, 2).unwrap();
        sel.mdl = 9;
        assert_eq!(bio_model_row_matches(&sel, &row(m2, Some(2))), Ok(false));
        assert_eq!(bio_model_row_matches(&sel, &row(m9, Some(9))), Ok(true));

        // Wildcard availability order: mdl == 0 selects regardless of the
        // row's number (documented availability order — wildcard first,
        // exact equality only when a concrete number is selected).
        sel.mdl = 0;
        assert_eq!(bio_model_row_matches(&sel, &row(zero, Some(-3))), Ok(true));
        assert_eq!(bio_model_row_matches(&sel, &row(m9, Some(9))), Ok(true));

        // None: typed missing-representation error, never a guessed value.
        assert_eq!(
            bio_model_row_matches(&sel, &row(zero, None)),
            Err(BioRowModelError::MissingModelNumber)
        );
        sel.mdl = 1;
        assert_eq!(
            bio_model_row_matches(&sel, &row(zero, None)),
            Err(BioRowModelError::MissingModelNumber)
        );
    }
    #[test]
    fn bio_rows_r04_chain_identity_precedence_and_missing() {
        use super::{BioRowChainError, Selection, SelectionList, bio_chain_row_matches};
        use crate::hierarchy::{BioChainRow, BioCoordinateFormat, ChainKind};
        use crate::source_ids::{ChainSourceIds, PdbChainId, ResidueName};

        fn chain(auth: Option<&[u8]>, label: Option<&str>) -> BioChainRow {
            let span = crate::hierarchy::BioRowSpan::new(0, 0).unwrap();
            let source = ChainSourceIds::new(
                auth.and_then(PdbChainId::from_ascii),
                label.map(str::to_string),
            );
            let _ = ResidueName::from_ascii(b"ALA");
            BioChainRow::new(
                crate::hierarchy::BioModelId::new(0),
                None,
                span,
                ChainKind::Protein,
                source,
            )
        }
        fn sel(list: &str) -> Selection {
            let mut s = Selection::default();
            s.chain_ids = SelectionList {
                all: false,
                inverted: false,
                list: list.to_string(),
            };
            s
        }

        // PDB blank chain column: read_string view of a raw blank field is
        // the EMPTY canonical name (not missing) — it matches "" only.
        let pdb = BioCoordinateFormat::Pdb;
        let blank = chain(Some(b"  "), Some("L"));
        assert_eq!(bio_chain_row_matches(&sel("A,B"), &blank, pdb), Ok(false));
        let empty_sel = sel("");
        assert_eq!(bio_chain_row_matches(&empty_sel, &blank, pdb), Ok(true));
        // PDB multichar (2-column) chain id, padded: "AB".
        let padded = chain(Some(b" AB"), Some("LAB"));
        assert_eq!(bio_chain_row_matches(&sel("AB"), &padded, pdb), Ok(true));
        assert_eq!(bio_chain_row_matches(&sel(" A"), &padded, pdb), Ok(false));
        // Same stored auth bytes under mmCIF provenance: verbatim (no
        // column trim) — auth " AB" matches " AB", not "AB".
        let cif = BioCoordinateFormat::Mmcif;
        let cif_auth = chain(Some(b" AB"), Some("LAB"));
        assert_eq!(bio_chain_row_matches(&sel(" AB"), &cif_auth, cif), Ok(true));
        assert_eq!(bio_chain_row_matches(&sel("AB"), &cif_auth, cif), Ok(false));
        // mmCIF auth-vs-label differing values: auth is preferred even
        // when the label would also match.
        let differing = chain(Some(b"A"), Some("LAB"));
        assert_eq!(bio_chain_row_matches(&sel("A"), &differing, cif), Ok(true));
        assert_eq!(
            bio_chain_row_matches(&sel("LAB"), &differing, cif),
            Ok(false)
        );
        // Label fallback only when auth is absent.
        let label_only = chain(None, Some("LAB"));
        assert_eq!(
            bio_chain_row_matches(&sel("LAB"), &label_only, cif),
            Ok(true)
        );
        assert_eq!(
            bio_chain_row_matches(&sel("A"), &label_only, cif),
            Ok(false)
        );
        // Both markers missing: typed error, never a guessed substitute.
        let missing = chain(None, None);
        assert_eq!(
            bio_chain_row_matches(&sel("A"), &missing, cif),
            Err(BioRowChainError::MissingCanonicalChainName)
        );
        // Wildcard list still demands a canonical name first.
        let mut all_sel = Selection::default();
        all_sel.chain_ids = SelectionList {
            all: true,
            inverted: false,
            list: String::new(),
        };
        assert_eq!(
            bio_chain_row_matches(&all_sel, &missing, cif),
            Err(BioRowChainError::MissingCanonicalChainName)
        );
    }
    #[test]
    fn bio_cid_c23_borrowed_access_and_unchanged_r04_fixtures() {
        use super::{BioRowChainError, Selection, SelectionList, bio_chain_row_matches};
        use crate::hierarchy::{BioChainRow, BioCoordinateFormat, ChainKind};
        use crate::source_ids::{ChainSourceIds, PdbChainId, ResidueName};

        // Borrowed accessor: pointer identity — two calls share the same
        // stored PdbChainId storage; absent auth yields None (label-only
        // rows take the label branch without any auth borrow).
        let ids = ChainSourceIds::new(PdbChainId::from_ascii(b"AB"), Some("LAB".into()));
        let first = ids.auth_chain_id_ref().expect("auth present");
        let second = ids.auth_chain_id_ref().expect("auth present");
        assert!(std::ptr::eq(first, second));
        assert_eq!(first.as_str(), "AB");
        let label_only = ChainSourceIds::new(None, Some("L".into()));
        assert!(label_only.auth_chain_id_ref().is_none());

        fn chain(auth: Option<&[u8]>, label: Option<&str>) -> BioChainRow {
            let span = crate::hierarchy::BioRowSpan::new(0, 0).unwrap();
            let source = ChainSourceIds::new(
                auth.and_then(PdbChainId::from_ascii),
                label.map(str::to_string),
            );
            let _ = ResidueName::from_ascii(b"ALA");
            BioChainRow::new(
                crate::hierarchy::BioModelId::new(0),
                None,
                span,
                ChainKind::Protein,
                source,
            )
        }
        fn sel(list: &str) -> Selection {
            let mut s = Selection::default();
            s.chain_ids = SelectionList {
                all: false,
                inverted: false,
                list: list.to_string(),
            };
            s
        }

        // Unchanged R04 fixtures through the allocation-free path: auth
        // preferred with the PDB read_string view, verbatim auth under
        // mmCIF, label fallback, and the typed missing-identity error.
        let pdb = BioCoordinateFormat::Pdb;
        let cif = BioCoordinateFormat::Mmcif;
        assert_eq!(
            bio_chain_row_matches(&sel("AB"), &chain(Some(b" AB"), Some("LAB")), pdb),
            Ok(true)
        );
        assert_eq!(
            bio_chain_row_matches(&sel(" A"), &chain(Some(b" AB"), Some("LAB")), pdb),
            Ok(false)
        );
        assert_eq!(
            bio_chain_row_matches(&sel(" AB"), &chain(Some(b" AB"), Some("LAB")), cif),
            Ok(true)
        );
        assert_eq!(
            bio_chain_row_matches(&sel("LAB"), &chain(None, Some("LAB")), cif),
            Ok(true)
        );
        assert_eq!(
            bio_chain_row_matches(&sel("A"), &chain(None, None), cif),
            Err(BioRowChainError::MissingCanonicalChainName)
        );
    }

    #[test]
    fn bio_cid_element_name_all_display_spellings_and_bounds() {
        use super::gemmi_element_name;
        // Independently literal display table (elem.hpp:269-294), NOT a
        // production-table oracle: mixed-case display spellings.
        let expected: [&str; 121] = [
            "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P",
            "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
            "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
            "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce",
            "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
            "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
            "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
            "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl",
            "Mc", "Lv", "Ts", "Og", "D", "",
        ];
        for (ordinal, name) in expected.iter().enumerate() {
            assert_eq!(gemmi_element_name(ordinal), Some(*name), "{ordinal}");
        }
        // Boundary: END (120) is the empty entry; 121 is beyond the
        // modeled vocabulary (safe Rust boundary, not a source-undefined
        // behavior claim).
        assert_eq!(gemmi_element_name(120), Some(""));
        assert_eq!(gemmi_element_name(121), None);
        // Display vocabulary differs from the parser's UPPERCASE table
        // exactly on the two-letter spellings (He/Li/Cl class).
        assert_eq!(super::GEMMI_ELEMENT_NAMES[2], "HE");
        assert_eq!(gemmi_element_name(2), Some("He"));
        assert_eq!(super::GEMMI_ELEMENT_NAMES[17], "CL");
        assert_eq!(gemmi_element_name(17), Some("Cl"));
    }

    #[test]
    fn bio_cid_c24_selected_atom_ids_zero_one_many_and_order() {
        use super::{BioSelectionData, BioSelectionParts, SelectionCidList};
        use crate::hierarchy::{
            BioAtomRow, BioCalcFlag, BioChainRow, BioCoordinateFormat, BioModelRow, BioResidueRow,
            BioRowSpan, ChainKind, EntityKind,
        };
        use crate::residue::ResidueInfoKind;
        use crate::source_ids::{
            AtomName, AtomSourceIds, ChainSourceIds, PdbChainId, PdbSeqId, ResidueName,
            ResidueSourceIds,
        };
        use cosmolkit_types::Element;

        fn model(num: Option<i32>, start: u32, len: u32) -> BioModelRow {
            BioModelRow::new(BioRowSpan::new(start, len).unwrap(), num)
        }
        fn chain_row(auth: Option<&[u8]>, start: u32, len: u32) -> BioChainRow {
            BioChainRow::new(
                crate::hierarchy::BioModelId::new(0),
                None,
                BioRowSpan::new(start, len).unwrap(),
                ChainKind::Protein,
                ChainSourceIds::new(auth.map(|a| PdbChainId::from_ascii(a).unwrap()), None),
            )
        }
        fn residue_row(name: &[u8], start: u32, len: u32) -> BioResidueRow {
            BioResidueRow::new(
                crate::hierarchy::BioChainId::new(0),
                BioRowSpan::new(start, len).unwrap(),
                ResidueName::from_ascii(name).unwrap(),
                ResidueInfoKind::Unknown,
                EntityKind::Polymer,
                None,
                None,
                ResidueSourceIds::new(Some(PdbSeqId::new(1, None)), None, None, None, None)
                    .unwrap(),
                crate::hierarchy::BioSiftsUnpResidue::default(),
            )
        }
        fn atom_row(name: &[u8]) -> BioAtomRow {
            BioAtomRow::new(
                crate::hierarchy::BioResidueId::new(0),
                AtomName::from_ascii(name).unwrap(),
                Element::from_atomic_number(6).unwrap(),
                None,
                None,
                0,
                BioCalcFlag::default(),
                1.0,
                20.0,
                [0.0; 6],
                0,
                0.0,
                AtomSourceIds::default(),
            )
        }
        fn list(text: &str) -> SelectionCidList {
            SelectionCidList {
                all: false,
                inverted: false,
                list: text.as_bytes().to_vec(),
            }
        }
        fn list_wildcard() -> SelectionCidList {
            // Source List default: all = true.
            SelectionCidList {
                all: true,
                inverted: false,
                list: Vec::new(),
            }
        }
        fn data(mdl: i32, atoms: &str, residues: &str) -> BioSelectionData {
            BioSelectionData::from_parts(BioSelectionParts {
                mdl,
                chain_ids: SelectionCidList {
                    all: true,
                    inverted: false,
                    list: Vec::new(),
                },
                from_seqid: crate::SelectionSequenceId {
                    seqnum: i32::MIN,
                    icode: b'*',
                },
                to_seqid: crate::SelectionSequenceId {
                    seqnum: i32::MAX,
                    icode: b'*',
                },
                residue_names: if residues.is_empty() {
                    list_wildcard()
                } else {
                    list(residues)
                },
                entity_types: list_wildcard(),
                et_flags: [false; 6],
                atom_names: if atoms.is_empty() {
                    list_wildcard()
                } else {
                    list(atoms)
                },
                elements: None,
                altlocs: list_wildcard(),
                atom_inequalities: Vec::new(),
            })
            .unwrap()
        }
        fn ids(data: &BioSelectionData) -> Vec<u32> {
            let models = [model(Some(1), 0, 1), model(Some(2), 1, 1)];
            let chains = [chain_row(Some(b"A"), 0, 2), chain_row(Some(b"A"), 2, 1)];
            let residues = [
                residue_row(b"ALA", 0, 2),
                residue_row(b"GLY", 2, 2),
                residue_row(b"HOH", 4, 1),
            ];
            let atoms = [
                atom_row(b" CA"),
                atom_row(b" N"),
                atom_row(b" CA"),
                atom_row(b" O"),
                atom_row(b" O"),
            ];
            data.selected_bio_atom_ids(
                &models,
                &chains,
                &residues,
                &atoms,
                BioCoordinateFormat::Mmcif,
            )
            .unwrap()
            .iter()
            .map(|id| id.value())
            .collect()
        }

        // All (zero constraints): every original ID in source order.
        assert_eq!(ids(&data(0, "", "")), vec![0u32, 1, 2, 3, 4]);
        // Many, sparse: " O" members only.
        assert_eq!(ids(&data(0, " O", "")), vec![3u32, 4]);
        // Rejected parent (residue gate): water residue's atom excluded
        // in BOTH models.
        assert_eq!(ids(&data(0, "", "ALA,GLY")), vec![0u32, 1, 2, 3]);
        // Rejected parent (model gate): model 2's chain row (its own span)
        // carries only the HOH residue, so model-1 selection is the four
        // atoms of ALA+GLY; atom 4 never reached.
        assert_eq!(ids(&data(1, "", "")), vec![0u32, 1, 2, 3]);
        // Exactly one hit.
        assert_eq!(ids(&data(0, " N", "")), vec![1u32]);
        // Zero hits: no matching atom name.
        assert_eq!(ids(&data(0, "ZZ", "")), Vec::<u32>::new());
    }

    #[test]
    fn bio_cid_c24_missing_source_representation_and_cursor_parity() {
        use super::{
            BioRowChainError, BioSelectionData, BioSelectionMatchError, BioSelectionParts,
        };
        use crate::hierarchy::{
            BioAtomRow, BioCalcFlag, BioChainRow, BioCoordinateFormat, BioModelRow, BioResidueRow,
            BioRowSpan, ChainKind, EntityKind,
        };
        use crate::residue::ResidueInfoKind;
        use crate::source_ids::{
            AtomName, AtomSourceIds, ChainSourceIds, PdbSeqId, ResidueName, ResidueSourceIds,
        };
        use cosmolkit_types::Element;

        let models = [BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1))];
        let chains = [BioChainRow::new(
            crate::hierarchy::BioModelId::new(0),
            None,
            BioRowSpan::new(0, 1).unwrap(),
            ChainKind::Protein,
            ChainSourceIds::new(None, None),
        )];
        let residues = [BioResidueRow::new(
            crate::hierarchy::BioChainId::new(0),
            BioRowSpan::new(0, 1).unwrap(),
            ResidueName::from_ascii(b"ALA").unwrap(),
            ResidueInfoKind::Unknown,
            EntityKind::Polymer,
            None,
            None,
            ResidueSourceIds::new(Some(PdbSeqId::new(1, None)), None, None, None, None).unwrap(),
            crate::hierarchy::BioSiftsUnpResidue::default(),
        )];
        let atoms = [BioAtomRow::new(
            crate::hierarchy::BioResidueId::new(0),
            AtomName::from_ascii(b" CA").unwrap(),
            Element::from_atomic_number(6).unwrap(),
            None,
            None,
            0,
            BioCalcFlag::default(),
            1.0,
            20.0,
            [0.0; 6],
            0,
            0.0,
            AtomSourceIds::default(),
        )];
        let data = BioSelectionData::default();
        let error = data
            .selected_bio_atom_ids(
                &models,
                &chains,
                &residues,
                &atoms,
                BioCoordinateFormat::Mmcif,
            )
            .unwrap_err();
        // Same-crate inspection: the private cause is the terminal chain
        // identity failure (typed, fail-closed, no empty-ID substitute).
        assert_eq!(
            error,
            BioSelectionMatchError::from_traverse_error(super::BioRowTraverseError::Chain(
                BioRowChainError::MissingCanonicalChainName
            ))
        );
        // Lazy-cursor regression: the underlying cursor still yields the
        // same first item/error directly (unchanged behavior, parity with
        // the detached query's mapped failure).
        let default_selection = super::Selection::default();
        let mut cursor = super::BioSelectedAtomIds::over_structure(
            &default_selection,
            &models,
            &chains,
            &residues,
            &atoms,
            BioCoordinateFormat::Mmcif,
            0,
            0,
        );
        assert_eq!(
            cursor.next(),
            Some(Err(super::BioRowTraverseError::Chain(
                BioRowChainError::MissingCanonicalChainName
            )))
        );
        assert_eq!(cursor.next(), None);
    }

    #[test]
    fn bio_cid_c24_finite_selection_product_48_and_error_chains() {
        use super::{
            BioRowChainError, BioRowModelError, BioRowTraverseError, BioSelectionData,
            BioSelectionMatchError, BioSelectionParts, SelectionCidList,
        };
        use crate::hierarchy::{
            BioAtomRow, BioCalcFlag, BioChainRow, BioCoordinateFormat, BioModelRow, BioResidueRow,
            BioRowSpan, ChainKind, EntityKind,
        };
        use crate::residue::ResidueInfoKind;
        use crate::source_ids::{
            AtomName, AtomSourceIds, ChainSourceIds, PdbChainId, PdbSeqId, ResidueName,
            ResidueSourceIds,
        };
        use cosmolkit_types::Element;
        use std::error::Error as _;

        fn model(num: Option<i32>, start: u32, len: u32) -> BioModelRow {
            BioModelRow::new(BioRowSpan::new(start, len).unwrap(), num)
        }
        fn chain_row(model_index: u32, auth: Option<&[u8]>, start: u32, len: u32) -> BioChainRow {
            BioChainRow::new(
                crate::hierarchy::BioModelId::new(model_index),
                None,
                BioRowSpan::new(start, len).unwrap(),
                ChainKind::Protein,
                ChainSourceIds::new(auth.map(|a| PdbChainId::from_ascii(a).unwrap()), None),
            )
        }
        fn residue_row(chain_index: u32, name: &[u8], start: u32, len: u32) -> BioResidueRow {
            BioResidueRow::new(
                crate::hierarchy::BioChainId::new(chain_index),
                BioRowSpan::new(start, len).unwrap(),
                ResidueName::from_ascii(name).unwrap(),
                ResidueInfoKind::Unknown,
                EntityKind::Polymer,
                None,
                None,
                ResidueSourceIds::new(Some(PdbSeqId::new(1, None)), None, None, None, None)
                    .unwrap(),
                crate::hierarchy::BioSiftsUnpResidue::default(),
            )
        }
        fn atom_row(residue_index: u32, name: &[u8]) -> BioAtomRow {
            BioAtomRow::new(
                crate::hierarchy::BioResidueId::new(residue_index),
                AtomName::from_ascii(name).unwrap(),
                Element::from_atomic_number(6).unwrap(),
                None,
                None,
                0,
                BioCalcFlag::default(),
                1.0,
                20.0,
                [0.0; 6],
                0,
                0.0,
                AtomSourceIds::default(),
            )
        }
        let wildcard = || SelectionCidList {
            all: true,
            inverted: false,
            list: Vec::new(),
        };
        let plain = |text: &str| SelectionCidList {
            all: false,
            inverted: false,
            list: text.as_bytes().to_vec(),
        };
        fn data(
            mdl: i32,
            chain: SelectionCidList,
            residues: SelectionCidList,
            atoms: SelectionCidList,
        ) -> BioSelectionData {
            let all_list = || SelectionCidList {
                all: true,
                inverted: false,
                list: Vec::new(),
            };
            BioSelectionData::from_parts(BioSelectionParts {
                mdl,
                chain_ids: chain,
                from_seqid: crate::SelectionSequenceId {
                    seqnum: i32::MIN,
                    icode: b'*',
                },
                to_seqid: crate::SelectionSequenceId {
                    seqnum: i32::MAX,
                    icode: b'*',
                },
                residue_names: residues,
                entity_types: all_list(),
                et_flags: [false; 6],
                atom_names: atoms,
                elements: None,
                altlocs: all_list(),
                atom_inequalities: Vec::new(),
            })
            .unwrap()
        }
        // Consistent reverse parent IDs (BIO-C24-PARENT): chain rows point
        // at their owning models [0,1], residue rows at their chains
        // [0,0,1], atom rows at their residues [0,0,1,1,2].
        let fixture = || {
            (
                vec![model(Some(1), 0, 1), model(Some(2), 1, 1)],
                vec![
                    chain_row(0, Some(b"A"), 0, 2),
                    chain_row(1, Some(b"A"), 2, 1),
                ],
                vec![
                    residue_row(0, b"ALA", 0, 2),
                    residue_row(0, b"GLY", 2, 2),
                    residue_row(1, b"HOH", 4, 1),
                ],
                vec![
                    atom_row(0, b" CA"),
                    atom_row(0, b" N"),
                    atom_row(1, b" CA"),
                    atom_row(1, b" O"),
                    atom_row(2, b" O"),
                ],
            )
        };
        let cif = BioCoordinateFormat::Mmcif;

        // 24 model/atom/residue rows with chain = all, literal original
        // ID arrays derived independently from the fixture spans:
        // model 1 = ALA(atoms 0,1)+GLY(atoms 2,3); model 2 = HOH(atom 4).
        let rows: [(i32, &str, &str, &[u32]); 24] = [
            (0, "", "", &[0, 1, 2, 3, 4]),
            (0, "", "ALA,GLY", &[0, 1, 2, 3]),
            (0, " O", "", &[3, 4]),
            (0, " O", "ALA,GLY", &[3]),
            (0, "ZZ", "", &[]),
            (0, "ZZ", "ALA,GLY", &[]),
            (1, "", "", &[0, 1, 2, 3]),
            (1, "", "ALA,GLY", &[0, 1, 2, 3]),
            (1, " O", "", &[3]),
            (1, " O", "ALA,GLY", &[3]),
            (1, "ZZ", "", &[]),
            (1, "ZZ", "ALA,GLY", &[]),
            (2, "", "", &[4]),
            (2, "", "ALA,GLY", &[]),
            (2, " O", "", &[4]),
            (2, " O", "ALA,GLY", &[]),
            (2, "ZZ", "", &[]),
            (2, "ZZ", "ALA,GLY", &[]),
            (99, "", "", &[]),
            (99, "", "ALA,GLY", &[]),
            (99, " O", "", &[]),
            (99, " O", "ALA,GLY", &[]),
            (99, "ZZ", "", &[]),
            (99, "ZZ", "ALA,GLY", &[]),
        ];
        let filter = |text: &str| {
            if text.is_empty() {
                wildcard()
            } else {
                plain(text)
            }
        };
        let mut calls = 0_usize;
        for (mdl, atoms, residues, expected) in rows {
            for chain in [wildcard(), plain("Z")] {
                let (models, chains, residue_rows, atom_rows) = fixture();
                {
                    let borrowed = (&models, &chains, &residue_rows, &atom_rows);
                    assert_eq!(
                        borrowed
                            .1
                            .iter()
                            .map(|row| row.model_id().value())
                            .collect::<Vec<_>>(),
                        [0, 1]
                    );
                    assert_eq!(
                        borrowed
                            .2
                            .iter()
                            .map(|row| row.chain_id().value())
                            .collect::<Vec<_>>(),
                        [0, 0, 1]
                    );
                    assert_eq!(
                        borrowed
                            .3
                            .iter()
                            .map(|row| row.residue_id().value())
                            .collect::<Vec<_>>(),
                        [0, 0, 1, 1, 2]
                    );
                }
                let result = data(mdl, chain.clone(), filter(residues), filter(atoms))
                    .selected_bio_atom_ids(&models, &chains, &residue_rows, &atom_rows, cif);
                calls += 1;
                let expected_empty = chain.list == b"Z";
                if expected_empty {
                    // Rejected Z chains: literal empty results.
                    assert_eq!(result.unwrap(), Vec::<crate::hierarchy::BioAtomId>::new());
                } else {
                    let got: Vec<u32> = result.unwrap().iter().map(|id| id.value()).collect();
                    assert_eq!(got, expected, "{mdl}/{atoms}/{residues}");
                }
            }
        }
        assert_eq!(calls, 48);

        // Missing-model leaf (fixture spans/parents all consistent: one
        // model/chain/residue/atom, spans 0..1): typed cause, full
        // Error::source chain, leaf source None, Err rather than partial
        // Ok; ALL supplied tables unchanged afterwards (snapshotted and
        // compared, not a single length check).
        let models_bad = vec![model(None, 0, 1)];
        let chains_ok = vec![chain_row(0, Some(b"A"), 0, 1)];
        let residues_one = vec![residue_row(0, b"ALA", 0, 1)];
        let atoms_one = vec![atom_row(0, b" CA")];
        let snapshot = (
            format!("{models_bad:?}"),
            format!("{chains_ok:?}"),
            format!("{residues_one:?}"),
            format!("{atoms_one:?}"),
        );
        let err = BioSelectionData::default()
            .selected_bio_atom_ids(&models_bad, &chains_ok, &residues_one, &atoms_one, cif)
            .unwrap_err();
        let traverse = err
            .source()
            .and_then(|s| s.downcast_ref::<BioRowTraverseError>())
            .expect("traverse cause retained");
        assert_eq!(
            traverse,
            &BioRowTraverseError::Model(BioRowModelError::MissingModelNumber)
        );
        let leaf = traverse
            .source()
            .and_then(|s| s.downcast_ref::<BioRowModelError>())
            .expect("model leaf retained");
        assert_eq!(leaf, &BioRowModelError::MissingModelNumber);
        assert!(leaf.source().is_none());
        assert_eq!(
            snapshot,
            (
                format!("{models_bad:?}"),
                format!("{chains_ok:?}"),
                format!("{residues_one:?}"),
                format!("{atoms_one:?}"),
            )
        );

        // Immediate missing-chain leaf: model Some(1) so the MODEL gate
        // passes and the failing branch is the CHAIN gate (the C24-CLOSE
        // version wrongly reused the None model and never reached Chain).
        let models_one = vec![model(Some(1), 0, 1)];
        let chains_missing = vec![chain_row(0, None, 0, 1)];
        let snapshot = format!("{models_one:?}{chains_missing:?}{residues_one:?}{atoms_one:?}");
        let err = BioSelectionData::default()
            .selected_bio_atom_ids(&models_one, &chains_missing, &residues_one, &atoms_one, cif)
            .unwrap_err();
        let traverse = err
            .source()
            .and_then(|s| s.downcast_ref::<BioRowTraverseError>())
            .expect("traverse cause retained");
        assert_eq!(
            traverse,
            &BioRowTraverseError::Chain(BioRowChainError::MissingCanonicalChainName)
        );
        let leaf = traverse
            .source()
            .and_then(|s| s.downcast_ref::<BioRowChainError>())
            .expect("chain leaf retained");
        assert_eq!(leaf, &BioRowChainError::MissingCanonicalChainName);
        assert!(leaf.source().is_none());
        assert_eq!(
            snapshot,
            format!("{models_one:?}{chains_missing:?}{residues_one:?}{atoms_one:?}")
        );

        // Late missing chain after accepted atoms: model 1 (chain A,
        // ALA+GLY atoms 0..4) yields first, then model 2's identity-less
        // chain fails. The ACTUAL lazy cursor must emit original IDs
        // 0,1,2,3, then the terminal Chain error, then None (sticky); the
        // collected query returns Err (all-or-error), never a partial Ok.
        let (models, chains_late, residue_rows2, atom_rows2) = (
            vec![model(Some(1), 0, 1), model(Some(2), 1, 1)],
            vec![chain_row(0, Some(b"A"), 0, 2), chain_row(1, None, 2, 1)],
            fixture().2,
            fixture().3,
        );
        assert_eq!(
            chains_late
                .iter()
                .map(|row| row.model_id().value())
                .collect::<Vec<_>>(),
            [0, 1]
        );
        assert_eq!(
            residue_rows2
                .iter()
                .map(|row| row.chain_id().value())
                .collect::<Vec<_>>(),
            [0, 0, 1]
        );
        assert_eq!(
            atom_rows2
                .iter()
                .map(|row| row.residue_id().value())
                .collect::<Vec<_>>(),
            [0, 0, 1, 1, 2]
        );
        let snapshot = format!("{models:?}{chains_late:?}{residue_rows2:?}{atom_rows2:?}");
        let default_selection = super::Selection::default();
        let mut cursor = super::BioSelectedAtomIds::over_structure(
            &default_selection,
            &models,
            &chains_late,
            &residue_rows2,
            &atom_rows2,
            cif,
            0,
            0,
        );
        assert_eq!(cursor.next().map(|r| r.map(|id| id.value())), Some(Ok(0)));
        assert_eq!(cursor.next().map(|r| r.map(|id| id.value())), Some(Ok(1)));
        assert_eq!(cursor.next().map(|r| r.map(|id| id.value())), Some(Ok(2)));
        assert_eq!(cursor.next().map(|r| r.map(|id| id.value())), Some(Ok(3)));
        assert_eq!(
            cursor.next(),
            Some(Err(BioRowTraverseError::Chain(
                BioRowChainError::MissingCanonicalChainName
            )))
        );
        assert_eq!(cursor.next(), None);
        let result = BioSelectionData::default().selected_bio_atom_ids(
            &models,
            &chains_late,
            &residue_rows2,
            &atom_rows2,
            cif,
        );
        assert!(result.is_err());
        let late_error = result.unwrap_err();
        let traverse = late_error
            .source()
            .and_then(|s| s.downcast_ref::<BioRowTraverseError>())
            .expect("traverse cause retained");
        assert_eq!(
            traverse,
            &BioRowTraverseError::Chain(BioRowChainError::MissingCanonicalChainName)
        );
        assert_eq!(
            snapshot,
            format!("{models:?}{chains_late:?}{residue_rows2:?}{atom_rows2:?}")
        );
    }

    #[test]
    fn bio_rows_r05_seqid_adapter_sentinels_and_bounds() {
        use super::{Selection, SelectionSequenceId, bio_residue_seqid};
        use crate::hierarchy::{BioResidueRow, EntityKind};
        use crate::residue::ResidueInfoKind;
        use crate::source_ids::{PdbSeqId, ResidueName, ResidueSourceIds};

        fn residue(seq: Option<PdbSeqId>) -> BioResidueRow {
            let span = crate::hierarchy::BioRowSpan::new(0, 0).unwrap();
            let source = ResidueSourceIds::new(seq, None, None, None, None).unwrap();
            BioResidueRow::new(
                crate::hierarchy::BioChainId::new(0),
                span,
                ResidueName::from_ascii(b"ALA").unwrap(),
                ResidueInfoKind::Unknown,
                EntityKind::Polymer,
                None,
                None,
                source,
                crate::hierarchy::BioSiftsUnpResidue::default(),
            )
        }

        // Missing seq_id: the pinned absence sentinels INT_MIN + blank.
        assert_eq!(bio_residue_seqid(&residue(None)), (i32::MIN, b' '));
        // Blank / A / B codes preserved byte-exactly with their numbers.
        assert_eq!(
            bio_residue_seqid(&residue(Some(PdbSeqId::new(10, None)))),
            (10, b' ')
        );
        assert_eq!(
            bio_residue_seqid(&residue(Some(PdbSeqId::new(10, Some(b'A'))))),
            (10, b'A')
        );
        assert_eq!(
            bio_residue_seqid(&residue(Some(PdbSeqId::new(-7, Some(b'B'))))),
            (-7, b'B')
        );
        // Negative, zero, positive numbers — including extremes.
        assert_eq!(
            bio_residue_seqid(&residue(Some(PdbSeqId::new(i32::MIN, None)))),
            (i32::MIN, b' ')
        );
        assert_eq!(
            bio_residue_seqid(&residue(Some(PdbSeqId::new(0, None)))),
            (0, b' ')
        );
        assert_eq!(
            bio_residue_seqid(&residue(Some(PdbSeqId::new(i32::MAX, Some(b'Z'))))),
            (i32::MAX, b'Z')
        );

        // Wildcard bounds through SelectionSequenceId::compare: the `*`
        // wildcard applies to the ICODE only; numbers compare directly,
        // and range membership is from.compare <= 0 <= to.compare.
        let from = SelectionSequenceId {
            seqnum: 5,
            icode: b'*',
        };
        let to = SelectionSequenceId {
            seqnum: 20,
            icode: b'*',
        };
        // Equal number: the icode wildcard makes the comparison 0.
        assert_eq!(from.compare(5, b'A'), 0);
        assert_eq!(to.compare(20, b' '), 0);
        // Inside the range: from.compare < 0 < to.compare.
        assert_eq!(from.compare(10, b' '), -1);
        assert_eq!(to.compare(10, b' '), 1);
        // Outside the range on either side: the selector-perspective
        // comparison is self.seqnum < source ? -1 : 1, so a source below
        // the from-bound compares +1 and a source above the to-bound -1.
        assert_eq!(from.compare(4, b' '), 1);
        assert_eq!(to.compare(21, b' '), -1);
        // INT_MIN sentinel compares below-from as +1 from the from-bound.
        assert_eq!(from.compare(i32::MIN, b' '), 1);
    }
    #[test]
    fn bio_rows_r06_residue_predicate_finite_product() {
        use super::{
            Selection, SelectionFlagList, SelectionList, SelectionSequenceId,
            bio_residue_row_matches,
        };
        use crate::hierarchy::{BioCoordinateFormat, EntityKind};
        use crate::residue::ResidueInfoKind;
        use crate::source_ids::{PdbSeqId, ResidueName, ResidueSourceIds};

        fn residue(
            entity: EntityKind,
            name: &[u8],
            seq: Option<PdbSeqId>,
        ) -> crate::hierarchy::BioResidueRow {
            let span = crate::hierarchy::BioRowSpan::new(0, 0).unwrap();
            let source = ResidueSourceIds::new(seq, None, None, None, None).unwrap();
            crate::hierarchy::BioResidueRow::new(
                crate::hierarchy::BioChainId::new(0),
                span,
                ResidueName::from_ascii(name).unwrap(),
                ResidueInfoKind::Unknown,
                entity,
                None,
                None,
                source,
                crate::hierarchy::BioSiftsUnpResidue::default(),
            )
        }

        fn sel() -> Selection {
            Selection::default()
        }
        let pdb = BioCoordinateFormat::Pdb;
        let unset_flag = b' ';

        // Entity kind: only the selected type passes the first conjunct.
        let mut only_polymer = sel();
        only_polymer.entity_types.all = false;
        only_polymer.et_flags[super::entity_kind_index(EntityKind::Polymer)] = true;
        assert!(bio_residue_row_matches(
            &only_polymer,
            &residue(EntityKind::Polymer, b"ALA", Some(PdbSeqId::new(1, None))),
            pdb,
            unset_flag
        ));
        assert!(!bio_residue_row_matches(
            &only_polymer,
            &residue(EntityKind::Water, b"HOH", Some(PdbSeqId::new(1, None))),
            pdb,
            unset_flag
        ));

        // Name: PDB-padded stored name matches the trimmed spelling.
        let mut ala_only = sel();
        ala_only.residue_names = SelectionList {
            all: false,
            inverted: false,
            list: "ALA".to_string(),
        };
        assert!(bio_residue_row_matches(
            &ala_only,
            &residue(EntityKind::Polymer, b"ALA ", Some(PdbSeqId::new(1, None))),
            pdb,
            unset_flag
        ));
        // ...while the same stored bytes under mmCIF keep "ALA " and fail.
        assert!(!bio_residue_row_matches(
            &ala_only,
            &residue(EntityKind::Polymer, b"ALA ", Some(PdbSeqId::new(1, None))),
            BioCoordinateFormat::Mmcif,
            unset_flag
        ));

        // Range: inclusive number bounds with icode wildcards.
        let mut band = sel();
        band.from_seqid = SelectionSequenceId {
            seqnum: 5,
            icode: b'*',
        };
        band.to_seqid = SelectionSequenceId {
            seqnum: 20,
            icode: b'*',
        };
        assert!(bio_residue_row_matches(
            &band,
            &residue(
                EntityKind::Polymer,
                b"ALA",
                Some(PdbSeqId::new(5, Some(b'A')))
            ),
            pdb,
            unset_flag
        ));
        assert!(bio_residue_row_matches(
            &band,
            &residue(EntityKind::Polymer, b"ALA", Some(PdbSeqId::new(20, None))),
            pdb,
            unset_flag
        ));
        assert!(!bio_residue_row_matches(
            &band,
            &residue(EntityKind::Polymer, b"ALA", Some(PdbSeqId::new(21, None))),
            pdb,
            unset_flag
        ));
        // Missing seq_id sentinel INT_MIN falls below the band.
        assert!(!bio_residue_row_matches(
            &band,
            &residue(EntityKind::Polymer, b"ALA", None),
            pdb,
            unset_flag
        ));

        // Icode letter outside the wildcard range when bounds pin a code.
        let mut exact = sel();
        exact.from_seqid = SelectionSequenceId {
            seqnum: 7,
            icode: b'A',
        };
        exact.to_seqid = SelectionSequenceId {
            seqnum: 7,
            icode: b'A',
        };
        assert!(bio_residue_row_matches(
            &exact,
            &residue(
                EntityKind::Polymer,
                b"ALA",
                Some(PdbSeqId::new(7, Some(b'A')))
            ),
            pdb,
            unset_flag
        ));
        assert!(!bio_residue_row_matches(
            &exact,
            &residue(
                EntityKind::Polymer,
                b"ALA",
                Some(PdbSeqId::new(7, Some(b'B')))
            ),
            pdb,
            unset_flag
        ));

        // Flag: the caller byte feeds the flag conjunct — inverted list.
        let mut no_flag = sel();
        no_flag.residue_flags = SelectionFlagList {
            pattern: "!X".to_string(),
        };
        assert!(bio_residue_row_matches(
            &no_flag,
            &residue(EntityKind::Polymer, b"ALA", Some(PdbSeqId::new(1, None))),
            pdb,
            unset_flag
        ));
        assert!(!bio_residue_row_matches(
            &no_flag,
            &residue(EntityKind::Polymer, b"ALA", Some(PdbSeqId::new(1, None))),
            pdb,
            b'X'
        ));

        // Source-ordered rejection: entity type fails before any name /
        // range / flag conjunct is consulted (a Water row never matches a
        // Polymer-only selection regardless of the other fields).
        assert!(!bio_residue_row_matches(
            &only_polymer,
            &residue(EntityKind::Water, b"ALA", Some(PdbSeqId::new(1, None))),
            pdb,
            b'X'
        ));
    }
    #[test]
    fn bio_rows_r07_element_ordinal_identity() {
        use super::bio_atom_element_ordinal;
        use crate::hierarchy::{BioAtomRow, BioCalcFlag};
        use crate::source_ids::{AtomName, AtomSourceIds};
        use cosmolkit_types::Element;

        fn atom(element: Element, isotope: Option<u16>) -> BioAtomRow {
            BioAtomRow::new(
                crate::hierarchy::BioResidueId::new(0),
                AtomName::from_ascii(b" CA").unwrap(),
                element,
                isotope,
                None,
                0,
                BioCalcFlag::default(),
                1.0,
                0.0,
                [0.0; 6],
                0,
                0.0,
                AtomSourceIds::default(),
            )
        }

        // X and H map to their ordinals 0 and 1.
        assert_eq!(
            bio_atom_element_ordinal(&atom(Element::from_atomic_number(0).unwrap(), None)),
            0
        );
        assert_eq!(
            bio_atom_element_ordinal(&atom(Element::from_atomic_number(1).unwrap(), None)),
            1
        );
        // Deuterium: the reader's approved H + Some(2) projection maps
        // back to the distinct Gemmi D ordinal 119, never plain H.
        assert_eq!(
            bio_atom_element_ordinal(&atom(Element::from_atomic_number(1).unwrap(), Some(2))),
            119
        );
        // Tritium stays H (no distinct Gemmi El for mass 3).
        assert_eq!(
            bio_atom_element_ordinal(&atom(Element::from_atomic_number(1).unwrap(), Some(3))),
            1
        );
        // Other isotopes keep their element ordinal: 13C -> C(6).
        assert_eq!(
            bio_atom_element_ordinal(&atom(Element::from_atomic_number(6).unwrap(), Some(13))),
            6
        );
        // Representative metals: Fe 26, Zn 30.
        assert_eq!(
            bio_atom_element_ordinal(&atom(Element::from_atomic_number(26).unwrap(), None)),
            26
        );
        assert_eq!(
            bio_atom_element_ordinal(&atom(Element::from_atomic_number(30).unwrap(), None)),
            30
        );
        // Og 118 — the largest real atomic number — is NOT D.
        assert_eq!(
            bio_atom_element_ordinal(&atom(Element::from_atomic_number(118).unwrap(), None)),
            118
        );
        // Distinguishable through the B09 mask: D selects the D ordinal
        // and not Og. (BIO-CID C02: explicit ordinals; the lexical parser
        // lives in the IO owner — same mask as parse "[D]" at pos 1.)
        let mut d_mask = [false; super::GEMMI_EL_END];
        d_mask[119] = true;
        assert!(d_mask[119]);
        assert!(!d_mask[118]);
        assert!(!d_mask[1]);
        // And H selects only ordinal 1 — deuterium is excluded.
        let mut h_mask = [false; super::GEMMI_EL_END];
        h_mask[1] = true;
        assert!(h_mask[1]);
        assert!(!h_mask[119]);
    }
    #[test]
    fn bio_rows_r08_atom_predicate_product() {
        use super::{
            Selection, SelectionAtomInequality, SelectionFlagList, SelectionList,
            bio_atom_row_matches,
        };
        use crate::hierarchy::{BioAtomRow, BioCalcFlag, BioCoordinateFormat};
        use crate::source_ids::{AltLocLabel, AtomName, AtomSourceIds};
        use cosmolkit_types::Element;

        fn atom(
            name: &[u8],
            altloc: Option<u8>,
            occ: f64,
            b_iso: f64,
            element: Element,
        ) -> BioAtomRow {
            BioAtomRow::new(
                crate::hierarchy::BioResidueId::new(0),
                AtomName::from_ascii(name).unwrap(),
                element,
                None,
                altloc.map(AltLocLabel::new),
                0,
                BioCalcFlag::default(),
                occ,
                b_iso,
                [0.0; 6],
                0,
                0.0,
                AtomSourceIds::default(),
            )
        }
        let h = Element::from_atomic_number(1).unwrap();
        let o = Element::from_atomic_number(8).unwrap();
        let pdb = BioCoordinateFormat::Pdb;
        let cif = BioCoordinateFormat::Mmcif;
        let unset_flag = b' ';

        // PDB name view: padded stored name matches the trimmed spelling.
        let mut ca_only = Selection::default();
        ca_only.atom_names = SelectionList {
            all: false,
            inverted: false,
            list: "CA".to_string(),
        };
        assert!(bio_atom_row_matches(
            &ca_only,
            &atom(b" CA ", None, 1.0, 20.0, h),
            pdb,
            unset_flag
        ));
        // CIF name view: the same stored bytes stay " CA " and fail.
        assert!(!bio_atom_row_matches(
            &ca_only,
            &atom(b" CA ", None, 1.0, 20.0, h),
            cif,
            unset_flag
        ));
        // CIF logical quoted name " CA " matches itself.
        let mut ca_spaced = Selection::default();
        ca_spaced.atom_names = SelectionList {
            all: false,
            inverted: false,
            list: " CA ".to_string(),
        };
        assert!(bio_atom_row_matches(
            &ca_spaced,
            &atom(b" CA ", None, 1.0, 20.0, h),
            cif,
            unset_flag
        ));

        // Altloc: blank label is the empty name; a present label is its
        // exact byte — no cross-matching.
        let mut alt_a = Selection::default();
        alt_a.altlocs = SelectionList {
            all: false,
            inverted: false,
            list: "A".to_string(),
        };
        assert!(bio_atom_row_matches(
            &alt_a,
            &atom(b" CA", Some(b'A'), 1.0, 20.0, h),
            cif,
            unset_flag
        ));
        assert!(!bio_atom_row_matches(
            &alt_a,
            &atom(b" CA", Some(b'B'), 1.0, 20.0, h),
            cif,
            unset_flag
        ));
        assert!(!bio_atom_row_matches(
            &alt_a,
            &atom(b" CA", None, 1.0, 20.0, h),
            cif,
            unset_flag
        ));
        let mut alt_blank = Selection::default();
        alt_blank.altlocs = SelectionList {
            all: false,
            inverted: false,
            list: String::new(),
        };
        assert!(bio_atom_row_matches(
            &alt_blank,
            &atom(b" CA", None, 1.0, 20.0, h),
            cif,
            unset_flag
        ));
        assert!(!bio_atom_row_matches(
            &alt_blank,
            &atom(b" CA", Some(b'A'), 1.0, 20.0, h),
            cif,
            unset_flag
        ));

        // Flag: caller byte through the inverted pattern.
        let mut no_x = Selection::default();
        no_x.atom_flags = SelectionFlagList {
            pattern: "!X".to_string(),
        };
        assert!(bio_atom_row_matches(
            &no_x,
            &atom(b" CA", None, 1.0, 20.0, h),
            cif,
            unset_flag
        ));
        assert!(!bio_atom_row_matches(
            &no_x,
            &atom(b" CA", None, 1.0, 20.0, h),
            cif,
            b'X'
        ));

        // Occupancy/B-factor thresholds (strict source inequalities on the
        // source float values).
        let mut occ_gt = Selection::default();
        occ_gt.atom_inequalities = vec![SelectionAtomInequality {
            property: b'q',
            relation: 1,
            value: 0.5,
        }];
        assert!(bio_atom_row_matches(
            &occ_gt,
            &atom(b" CA", None, 1.0, 20.0, h),
            cif,
            unset_flag
        ));
        assert!(!bio_atom_row_matches(
            &occ_gt,
            &atom(b" CA", None, 0.5, 20.0, h),
            cif,
            unset_flag
        ));
        assert!(!bio_atom_row_matches(
            &occ_gt,
            &atom(b" CA", None, 0.25, 20.0, h),
            cif,
            unset_flag
        ));
        let mut b_lt = Selection::default();
        b_lt.atom_inequalities = vec![SelectionAtomInequality {
            property: b'b',
            relation: -1,
            value: 40.0,
        }];
        assert!(bio_atom_row_matches(
            &b_lt,
            &atom(b" CA", None, 1.0, 20.0, h),
            cif,
            unset_flag
        ));
        assert!(!bio_atom_row_matches(
            &b_lt,
            &atom(b" CA", None, 1.0, 40.0, h),
            cif,
            unset_flag
        ));

        // Element conjunct: [O] restricts; deuterium stays distinct.
        // (BIO-CID C02: explicit ordinals — O=8, H=1, D=119; same masks
        // the IO owner's parse_cid_elements("[O]"/"[H]"/"[D]", 1)
        // produces.)
        let mut oxy = Selection::default();
        let mut o_mask = [false; super::GEMMI_EL_END];
        o_mask[8] = true;
        oxy.elements = Some(o_mask);
        assert!(bio_atom_row_matches(
            &oxy,
            &atom(b" O", None, 1.0, 20.0, o),
            cif,
            unset_flag
        ));
        assert!(!bio_atom_row_matches(
            &oxy,
            &atom(b" CA", None, 1.0, 20.0, h),
            cif,
            unset_flag
        ));
        let d_atom = BioAtomRow::new(
            crate::hierarchy::BioResidueId::new(0),
            AtomName::from_ascii(b" D").unwrap(),
            h,
            Some(2),
            None,
            0,
            BioCalcFlag::default(),
            1.0,
            20.0,
            [0.0; 6],
            0,
            0.0,
            AtomSourceIds::default(),
        );
        let mut hyd = Selection::default();
        let mut h_mask = [false; super::GEMMI_EL_END];
        h_mask[1] = true;
        hyd.elements = Some(h_mask);
        assert!(!bio_atom_row_matches(&hyd, &d_atom, cif, unset_flag));
        let mut deu = Selection::default();
        let mut d_mask = [false; super::GEMMI_EL_END];
        d_mask[119] = true;
        deu.elements = Some(d_mask);
        assert!(bio_atom_row_matches(&deu, &d_atom, cif, unset_flag));
    }
    #[test]
    fn bio_rows_r09_first_in_model_traversal() {
        use super::{
            BioRowModelError, BioRowTraverseError, Selection, SelectionList,
            bio_first_in_model_rows,
        };
        use crate::hierarchy::{
            BioAtomRow, BioCalcFlag, BioChainRow, BioCoordinateFormat, BioEntityId, BioModelRow,
            BioResidueRow, BioRowSpan, ChainKind, EntityKind,
        };
        use crate::residue::ResidueInfoKind;
        use crate::source_ids::{
            AtomName, AtomSourceIds, ChainSourceIds, PdbChainId, PdbSeqId, ResidueName,
            ResidueSourceIds,
        };
        use cosmolkit_types::Element;

        fn chain_row(model_id: u32, auth: &[u8], start: u32, len: u32) -> BioChainRow {
            BioChainRow::new(
                crate::hierarchy::BioModelId::new(model_id),
                None,
                BioRowSpan::new(start, len).unwrap(),
                ChainKind::Protein,
                ChainSourceIds::new(Some(PdbChainId::from_ascii(auth).unwrap()), None),
            )
        }
        fn residue_row(chain: u32, name: &[u8], seq: i32, start: u32, len: u32) -> BioResidueRow {
            BioResidueRow::new(
                crate::hierarchy::BioChainId::new(chain),
                BioRowSpan::new(start, len).unwrap(),
                ResidueName::from_ascii(name).unwrap(),
                ResidueInfoKind::Unknown,
                EntityKind::Polymer,
                None,
                None,
                ResidueSourceIds::new(Some(PdbSeqId::new(seq, None)), None, None, None, None)
                    .unwrap(),
                crate::hierarchy::BioSiftsUnpResidue::default(),
            )
        }
        fn atom_row(residue: u32, name: &[u8], element: u8) -> BioAtomRow {
            BioAtomRow::new(
                crate::hierarchy::BioResidueId::new(residue),
                AtomName::from_ascii(name).unwrap(),
                Element::from_atomic_number(element).unwrap(),
                None,
                None,
                0,
                BioCalcFlag::default(),
                1.0,
                20.0,
                [0.0; 6],
                0,
                0.0,
                AtomSourceIds::default(),
            )
        }

        // Two chains (A: 2 residues, B: 1), three atoms in A's residues.
        let chains = vec![chain_row(0, b"A", 0, 2), chain_row(0, b"B", 2, 1)];
        let residues = vec![
            residue_row(0, b"ALA", 1, 0, 2),
            residue_row(0, b"GLY", 2, 2, 1),
            residue_row(1, b"HOH", 3, 3, 1),
        ];
        let atoms = vec![
            atom_row(0, b" N", 7),
            atom_row(0, b" CA", 6),
            atom_row(1, b" CA", 6),
            atom_row(2, b" O", 8),
        ];
        let model = BioModelRow::new(BioRowSpan::new(0, 2).unwrap(), Some(1));
        let cif = BioCoordinateFormat::Mmcif;
        let unset = b' ';

        // Everything wildcard: first atom in source order.
        let all = Selection::default();
        assert_eq!(
            bio_first_in_model_rows(&all, &model, &chains, &residues, &atoms, cif, unset, unset),
            Ok(Some(super::BioRowFirst {
                chain: crate::hierarchy::BioChainId::new(0),
                residue: crate::hierarchy::BioResidueId::new(0),
                atom: crate::hierarchy::BioAtomId::new(0),
            }))
        );

        // Middle hit: CA of the first residue's second atom (skips N).
        // PDB provenance: the read_string view trims the stored " CA".
        let pdb = BioCoordinateFormat::Pdb;
        let mut ca_only = Selection::default();
        ca_only.atom_names = SelectionList {
            all: false,
            inverted: false,
            list: "CA".to_string(),
        };
        assert_eq!(
            bio_first_in_model_rows(
                &ca_only, &model, &chains, &residues, &atoms, pdb, unset, unset
            ),
            Ok(Some(super::BioRowFirst {
                chain: crate::hierarchy::BioChainId::new(0),
                residue: crate::hierarchy::BioResidueId::new(0),
                atom: crate::hierarchy::BioAtomId::new(1),
            }))
        );
        // ...while under mmCIF the same stored " CA" stays verbatim and
        // the member "CA" does not match: no hit.
        assert_eq!(
            bio_first_in_model_rows(
                &ca_only, &model, &chains, &residues, &atoms, cif, unset, unset
            ),
            Ok(None)
        );

        // Later hit in chain B (noncontiguous identities preserved).
        let mut water_res = Selection::default();
        water_res.residue_names = SelectionList {
            all: false,
            inverted: false,
            list: "HOH".to_string(),
        };
        assert_eq!(
            bio_first_in_model_rows(
                &water_res, &model, &chains, &residues, &atoms, cif, unset, unset
            ),
            Ok(Some(super::BioRowFirst {
                chain: crate::hierarchy::BioChainId::new(1),
                residue: crate::hierarchy::BioResidueId::new(2),
                atom: crate::hierarchy::BioAtomId::new(3),
            }))
        );

        // No match at all.
        let mut zz = Selection::default();
        zz.atom_names = SelectionList {
            all: false,
            inverted: false,
            list: "ZZ".to_string(),
        };
        assert_eq!(
            bio_first_in_model_rows(&zz, &model, &chains, &residues, &atoms, cif, unset, unset),
            Ok(None)
        );

        // Rejected model gate: a different model number.
        let mut mdl2 = Selection::default();
        mdl2.mdl = 2;
        assert_eq!(
            bio_first_in_model_rows(&mdl2, &model, &chains, &residues, &atoms, cif, unset, unset),
            Ok(None)
        );

        // Rejected chain gate (chain list) skips to chain B's atom.
        let mut only_b = Selection::default();
        only_b.chain_ids = SelectionList {
            all: false,
            inverted: false,
            list: "B".to_string(),
        };
        assert_eq!(
            bio_first_in_model_rows(
                &only_b, &model, &chains, &residues, &atoms, cif, unset, unset
            ),
            Ok(Some(super::BioRowFirst {
                chain: crate::hierarchy::BioChainId::new(1),
                residue: crate::hierarchy::BioResidueId::new(2),
                atom: crate::hierarchy::BioAtomId::new(3),
            }))
        );

        // Empty parents: model with an empty chain span -> None.
        let empty_model = BioModelRow::new(BioRowSpan::new(0, 0).unwrap(), Some(1));
        assert_eq!(
            bio_first_in_model_rows(
                &all,
                &empty_model,
                &chains,
                &residues,
                &atoms,
                cif,
                unset,
                unset
            ),
            Ok(None)
        );

        // Missing source model number: typed traversal error.
        let missing_model = BioModelRow::new(BioRowSpan::new(0, 2).unwrap(), None);
        assert_eq!(
            bio_first_in_model_rows(
                &all,
                &missing_model,
                &chains,
                &residues,
                &atoms,
                cif,
                unset,
                unset
            ),
            Err(BioRowTraverseError::Model(
                BioRowModelError::MissingModelNumber
            ))
        );
        // Missing canonical chain identity: typed traversal error.
        let mut bad_chains = vec![chain_row(0, b"A", 0, 1)];
        bad_chains[0] = BioChainRow::new(
            crate::hierarchy::BioModelId::new(0),
            None,
            BioRowSpan::new(0, 1).unwrap(),
            ChainKind::Protein,
            ChainSourceIds::new(None, None),
        );
        let one_chain_model = BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1));
        assert_eq!(
            bio_first_in_model_rows(
                &all,
                &one_chain_model,
                &bad_chains,
                &residues,
                &atoms,
                cif,
                unset,
                unset
            ),
            Err(BioRowTraverseError::Chain(
                super::BioRowChainError::MissingCanonicalChainName
            ))
        );
        // ...but a chain WITHOUT a matching chain predicate is skipped
        // without consulting its identity at all? No: the source gates
        // matches(chain) first, and identity IS the match — a missing
        // identity still errors even when the list would not match.
    }
    #[test]
    fn bio_rows_r10_whole_structure_first() {
        use super::{Selection, SelectionList, bio_first_rows};
        use crate::hierarchy::{
            BioAtomRow, BioCalcFlag, BioChainRow, BioCoordinateFormat, BioModelRow, BioResidueRow,
            BioRowSpan, ChainKind, EntityKind,
        };
        use crate::residue::ResidueInfoKind;
        use crate::source_ids::{
            AtomName, AtomSourceIds, ChainSourceIds, PdbChainId, PdbSeqId, ResidueName,
            ResidueSourceIds,
        };
        use cosmolkit_types::Element;

        fn model(num: Option<i32>, start: u32, len: u32) -> BioModelRow {
            BioModelRow::new(BioRowSpan::new(start, len).unwrap(), num)
        }
        fn chain_row(auth: &[u8], start: u32, len: u32) -> BioChainRow {
            BioChainRow::new(
                crate::hierarchy::BioModelId::new(0),
                None,
                BioRowSpan::new(start, len).unwrap(),
                ChainKind::Protein,
                ChainSourceIds::new(Some(PdbChainId::from_ascii(auth).unwrap()), None),
            )
        }
        fn residue_row(name: &[u8], seq: i32, start: u32, len: u32) -> BioResidueRow {
            BioResidueRow::new(
                crate::hierarchy::BioChainId::new(0),
                BioRowSpan::new(start, len).unwrap(),
                ResidueName::from_ascii(name).unwrap(),
                ResidueInfoKind::Unknown,
                EntityKind::Polymer,
                None,
                None,
                ResidueSourceIds::new(Some(PdbSeqId::new(seq, None)), None, None, None, None)
                    .unwrap(),
                crate::hierarchy::BioSiftsUnpResidue::default(),
            )
        }
        fn atom_row(name: &[u8]) -> BioAtomRow {
            BioAtomRow::new(
                crate::hierarchy::BioResidueId::new(0),
                AtomName::from_ascii(name).unwrap(),
                Element::from_atomic_number(6).unwrap(),
                None,
                None,
                0,
                BioCalcFlag::default(),
                1.0,
                20.0,
                [0.0; 6],
                0,
                0.0,
                AtomSourceIds::default(),
            )
        }

        // Two models; model 2 (noncontiguous number 9) holds the only ZN.
        let models = vec![model(Some(1), 0, 1), model(Some(9), 1, 1)];
        let chains = vec![chain_row(b"A", 0, 1), chain_row(b"A", 1, 1)];
        let residues = vec![residue_row(b"ALA", 1, 0, 1), residue_row(b"GLY", 1, 1, 1)];
        let atoms = vec![atom_row(b" CA"), atom_row(b" ZN")];
        let cif = BioCoordinateFormat::Mmcif;
        let unset = b' ';

        // Wildcard: first atom of the FIRST model in stored order.
        let all = Selection::default();
        assert_eq!(
            bio_first_rows(&all, &models, &chains, &residues, &atoms, cif, unset, unset),
            Ok(Some(super::BioRowFirstHit {
                model: crate::hierarchy::BioModelId::new(0),
                cra: super::BioRowFirst {
                    chain: crate::hierarchy::BioChainId::new(0),
                    residue: crate::hierarchy::BioResidueId::new(0),
                    atom: crate::hierarchy::BioAtomId::new(0),
                },
            }))
        );

        // Later hit in the second model by atom name (mmCIF verbatim).
        let mut zn = Selection::default();
        zn.atom_names = SelectionList {
            all: false,
            inverted: false,
            list: " ZN".to_string(),
        };
        assert_eq!(
            bio_first_rows(&zn, &models, &chains, &residues, &atoms, cif, unset, unset),
            Ok(Some(super::BioRowFirstHit {
                model: crate::hierarchy::BioModelId::new(1),
                cra: super::BioRowFirst {
                    chain: crate::hierarchy::BioChainId::new(1),
                    residue: crate::hierarchy::BioResidueId::new(1),
                    atom: crate::hierarchy::BioAtomId::new(1),
                },
            }))
        );

        // No hit anywhere.
        let mut none = Selection::default();
        none.atom_names = SelectionList {
            all: false,
            inverted: false,
            list: "XX".to_string(),
        };
        assert_eq!(
            bio_first_rows(
                &none, &models, &chains, &residues, &atoms, cif, unset, unset
            ),
            Ok(None)
        );

        // Stored order, NOT sorted by model number: swap the stored order
        // so number 9 comes first; the wildcard now hits it first.
        let swapped = vec![model(Some(9), 1, 1), model(Some(1), 0, 1)];
        assert_eq!(
            bio_first_rows(
                &all, &swapped, &chains, &residues, &atoms, cif, unset, unset
            ),
            Ok(Some(super::BioRowFirstHit {
                model: crate::hierarchy::BioModelId::new(0),
                cra: super::BioRowFirst {
                    chain: crate::hierarchy::BioChainId::new(1),
                    residue: crate::hierarchy::BioResidueId::new(1),
                    atom: crate::hierarchy::BioAtomId::new(1),
                },
            }))
        );

        // Duplicate chain names across models: chain B appears in both;
        // a chain-B-only selection hits the first stored occurrence.
        let mut chain_b = Selection::default();
        chain_b.chain_ids = SelectionList {
            all: false,
            inverted: false,
            list: "A".to_string(),
        };
        assert_eq!(
            bio_first_rows(
                &chain_b, &models, &chains, &residues, &atoms, cif, unset, unset
            ),
            Ok(Some(super::BioRowFirstHit {
                model: crate::hierarchy::BioModelId::new(0),
                cra: super::BioRowFirst {
                    chain: crate::hierarchy::BioChainId::new(0),
                    residue: crate::hierarchy::BioResidueId::new(0),
                    atom: crate::hierarchy::BioAtomId::new(0),
                },
            }))
        );

        // Missing model number in ANY visited model errors structurally;
        // an earlier model without a match is fine, the error surfaces
        // when its gate is consulted.
        let mut broken = vec![model(Some(1), 0, 1), model(None, 1, 1)];
        broken[0] = model(Some(1), 0, 1);
        let mut later_only = Selection::default();
        later_only.mdl = 9;
        assert_eq!(
            bio_first_rows(
                &later_only,
                &broken,
                &chains,
                &residues,
                &atoms,
                cif,
                unset,
                unset
            ),
            Err(super::BioRowTraverseError::Model(
                super::BioRowModelError::MissingModelNumber
            ))
        );
    }
    #[test]
    fn bio_rows_r11_lazy_atom_id_iterator() {
        use super::{
            BioRowModelError, BioSelectedAtomIds, Selection, SelectionList, bio_model_row_matches,
        };
        use crate::hierarchy::{
            BioAtomRow, BioCalcFlag, BioChainRow, BioCoordinateFormat, BioModelRow, BioResidueRow,
            BioRowSpan, ChainKind, EntityKind,
        };
        use crate::residue::ResidueInfoKind;
        use crate::source_ids::{
            AtomName, AtomSourceIds, ChainSourceIds, PdbChainId, PdbSeqId, ResidueName,
            ResidueSourceIds,
        };
        use cosmolkit_types::Element;

        fn model(num: Option<i32>, start: u32, len: u32) -> BioModelRow {
            BioModelRow::new(BioRowSpan::new(start, len).unwrap(), num)
        }
        fn chain_row(auth: &[u8], start: u32, len: u32) -> BioChainRow {
            BioChainRow::new(
                crate::hierarchy::BioModelId::new(0),
                None,
                BioRowSpan::new(start, len).unwrap(),
                ChainKind::Protein,
                ChainSourceIds::new(Some(PdbChainId::from_ascii(auth).unwrap()), None),
            )
        }
        fn residue_row(name: &[u8], start: u32, len: u32) -> BioResidueRow {
            BioResidueRow::new(
                crate::hierarchy::BioChainId::new(0),
                BioRowSpan::new(start, len).unwrap(),
                ResidueName::from_ascii(name).unwrap(),
                ResidueInfoKind::Unknown,
                EntityKind::Polymer,
                None,
                None,
                ResidueSourceIds::new(Some(PdbSeqId::new(1, None)), None, None, None, None)
                    .unwrap(),
                crate::hierarchy::BioSiftsUnpResidue::default(),
            )
        }
        fn atom_row(name: &[u8]) -> BioAtomRow {
            BioAtomRow::new(
                crate::hierarchy::BioResidueId::new(0),
                AtomName::from_ascii(name).unwrap(),
                Element::from_atomic_number(6).unwrap(),
                None,
                None,
                0,
                BioCalcFlag::default(),
                1.0,
                20.0,
                [0.0; 6],
                0,
                0.0,
                AtomSourceIds::default(),
            )
        }

        // Model 1 (chain A: 2 residues) and model 2 (chain A: 1 residue).
        let models = vec![model(Some(1), 0, 1), model(Some(2), 1, 1)];
        let chains = vec![chain_row(b"A", 0, 2), chain_row(b"A", 2, 1)];
        let residues = vec![
            residue_row(b"ALA", 0, 2),
            residue_row(b"GLY", 2, 2),
            residue_row(b"HOH", 4, 1),
        ];
        let atoms = vec![
            atom_row(b" CA"),
            atom_row(b" N"),
            atom_row(b" CA"),
            atom_row(b" O"),
            atom_row(b" O"),
        ];
        let cif = BioCoordinateFormat::Mmcif;
        let unset = b' ';

        // All hits: every atom ID in exact source order.
        let all = Selection::default();
        let ids: Vec<_> = BioSelectedAtomIds::over_structure(
            &all, &models, &chains, &residues, &atoms, cif, unset, unset,
        )
        .map(|r| r.unwrap().value())
        .collect();
        assert_eq!(ids, vec![0u32, 1, 2, 3, 4]);

        // Sparse hits: " O" members only (mmCIF verbatim names).
        let mut oxy = Selection::default();
        oxy.atom_names = SelectionList {
            all: false,
            inverted: false,
            list: " O".to_string(),
        };
        let ids: Vec<_> = BioSelectedAtomIds::over_structure(
            &oxy, &models, &chains, &residues, &atoms, cif, unset, unset,
        )
        .map(|r| r.unwrap().value())
        .collect();
        assert_eq!(ids, vec![3u32, 4]);

        // Dense-but-partial: parent gating — water residue excluded by
        // name list; model 2's HOH atom 4 never consulted for its flag.
        let mut poly = Selection::default();
        poly.residue_names = SelectionList {
            all: false,
            inverted: false,
            list: "ALA,GLY".to_string(),
        };
        let ids: Vec<_> = BioSelectedAtomIds::over_structure(
            &poly, &models, &chains, &residues, &atoms, cif, unset, unset,
        )
        .map(|r| r.unwrap().value())
        .collect();
        assert_eq!(ids, vec![0u32, 1, 2, 3]);

        // Early drop: take(2) stops without exhausting the tables.
        let ids: Vec<_> = BioSelectedAtomIds::over_structure(
            &all, &models, &chains, &residues, &atoms, cif, unset, unset,
        )
        .take(2)
        .map(|r| r.unwrap().value())
        .collect();
        assert_eq!(ids, vec![0u32, 1]);

        // No hits at all.
        let mut none = Selection::default();
        none.atom_names = SelectionList {
            all: false,
            inverted: false,
            list: "ZZ".to_string(),
        };
        let mut iter = BioSelectedAtomIds::over_structure(
            &none, &models, &chains, &residues, &atoms, cif, unset, unset,
        );
        assert!(iter.next().is_none());

        // Terminal error from a missing model number in a VISITED model:
        // laziness means model 1's four atoms are yielded first, and the
        // error surfaces exactly when the broken model's gate is
        // consulted; afterwards the iterator is exhausted.
        let broken = vec![model(Some(1), 0, 1), model(None, 1, 1)];
        let mut err_iter = BioSelectedAtomIds::over_structure(
            &all, &broken, &chains, &residues, &atoms, cif, unset, unset,
        );
        for expected in 0u32..4 {
            assert_eq!(
                err_iter.next(),
                Some(Ok(crate::hierarchy::BioAtomId::new(expected)))
            );
        }
        assert_eq!(
            err_iter.next(),
            Some(Err(super::BioRowTraverseError::Model(
                BioRowModelError::MissingModelNumber
            )))
        );
        assert!(err_iter.next().is_none());

        // Empty hierarchy.
        let empty_iter =
            BioSelectedAtomIds::over_structure(&all, &[], &[], &[], &[], cif, unset, unset);
        assert_eq!(empty_iter.count(), 0);
    }
    #[test]
    fn bio_rows_cursor_gate_measurements_are_thread_local() {
        use super::{Selection, gate_counters};

        let row = crate::BioModelRow::new(crate::BioRowSpan::new(0, 0).unwrap(), Some(1));
        gate_counters::reset();
        assert!(super::bio_model_row_matches(&Selection::default(), &row).unwrap());
        std::thread::spawn(|| {
            let row = crate::BioModelRow::new(crate::BioRowSpan::new(0, 0).unwrap(), Some(1));
            gate_counters::reset();
            for _ in 0..2 {
                assert!(super::bio_model_row_matches(&Selection::default(), &row).unwrap());
            }
            assert_eq!(gate_counters::MODEL.with(std::cell::Cell::get), 2);
        })
        .join()
        .unwrap();
        assert_eq!(gate_counters::MODEL.with(std::cell::Cell::get), 1);
    }

    #[test]
    fn bio_rows_cursor_parent_gates_once_per_visited_row() {
        use super::{BioSelectedAtomIds, Selection, gate_counters};
        use crate::hierarchy::{
            BioAtomRow, BioCalcFlag, BioChainRow, BioCoordinateFormat, BioModelRow, BioResidueRow,
            BioRowSpan, ChainKind, EntityKind,
        };
        use crate::residue::ResidueInfoKind;
        use crate::source_ids::{
            AtomName, AtomSourceIds, ChainSourceIds, PdbChainId, PdbSeqId, ResidueName,
            ResidueSourceIds,
        };
        use cosmolkit_types::Element;

        fn model(num: i32, start: u32, len: u32) -> BioModelRow {
            BioModelRow::new(BioRowSpan::new(start, len).unwrap(), Some(num))
        }
        fn chain_row(auth: &[u8], start: u32, len: u32) -> BioChainRow {
            BioChainRow::new(
                crate::hierarchy::BioModelId::new(0),
                None,
                BioRowSpan::new(start, len).unwrap(),
                ChainKind::Protein,
                ChainSourceIds::new(Some(PdbChainId::from_ascii(auth).unwrap()), None),
            )
        }
        fn residue_row(name: &[u8], start: u32, len: u32) -> BioResidueRow {
            BioResidueRow::new(
                crate::hierarchy::BioChainId::new(0),
                BioRowSpan::new(start, len).unwrap(),
                ResidueName::from_ascii(name).unwrap(),
                ResidueInfoKind::Unknown,
                EntityKind::Polymer,
                None,
                None,
                ResidueSourceIds::new(Some(PdbSeqId::new(1, None)), None, None, None, None)
                    .unwrap(),
                crate::hierarchy::BioSiftsUnpResidue::default(),
            )
        }
        fn atom_row(name: &[u8]) -> BioAtomRow {
            BioAtomRow::new(
                crate::hierarchy::BioResidueId::new(0),
                AtomName::from_ascii(name).unwrap(),
                Element::from_atomic_number(6).unwrap(),
                None,
                None,
                0,
                BioCalcFlag::default(),
                1.0,
                20.0,
                [0.0; 6],
                0,
                0.0,
                AtomSourceIds::default(),
            )
        }

        // 2 models x (2 chains x 2 residues x 2 atoms) — 16 atoms; the
        // second model's second chain is rejected, mixing accepted and
        // rejected parents at every level.
        let models = vec![model(1, 0, 2), model(2, 2, 2)];
        let chains = vec![
            chain_row(b"A", 0, 2),
            chain_row(b"B", 2, 2),
            chain_row(b"A", 4, 2),
            chain_row(b"B", 6, 2),
        ];
        let residues: Vec<_> = (0..8).map(|i| residue_row(b"ALA", i * 2, 2)).collect();
        let atoms: Vec<_> = (0..16).map(|_| atom_row(b" CA")).collect();
        let cif = BioCoordinateFormat::Mmcif;
        let unset = b' ';

        // Wildcard selection: all 16 atoms, chain B rejected only in
        // model 2 via a targeted list is impossible here (both models
        // have B) — instead reject chain "B" everywhere and verify counts
        // against the visited rows.
        let mut only_a = Selection::default();
        only_a.chain_ids = super::SelectionList {
            all: false,
            inverted: false,
            list: "A".to_string(),
        };
        gate_counters::reset();
        let ids: Vec<_> = BioSelectedAtomIds::over_structure(
            &only_a, &models, &chains, &residues, &atoms, cif, unset, unset,
        )
        .map(|r| r.unwrap().value())
        .collect();
        assert_eq!(ids.len(), 8); // chains A of both models, 2 residues x 2 atoms each
        // Model gate: exactly once per model row (2).
        assert_eq!(gate_counters::MODEL.with(std::cell::Cell::get), 2);
        // Chain gate: once per visited chain row of accepted models (4).
        assert_eq!(gate_counters::CHAIN.with(std::cell::Cell::get), 4);
        // Residue gate: once per visited residue row of accepted chains (4).
        assert_eq!(gate_counters::RESIDUE.with(std::cell::Cell::get), 4);

        // Full wildcard: every parent accepted; with early drop at 3 atoms
        // the counts reflect the partial visitation only.
        let all = Selection::default();
        gate_counters::reset();
        let ids: Vec<_> = BioSelectedAtomIds::over_structure(
            &all, &models, &chains, &residues, &atoms, cif, unset, unset,
        )
        .take(3)
        .map(|r| r.unwrap().value())
        .collect();
        assert_eq!(ids, vec![0u32, 1, 2]);
        assert_eq!(gate_counters::MODEL.with(std::cell::Cell::get), 1);
        assert_eq!(gate_counters::CHAIN.with(std::cell::Cell::get), 1);
        // Residue 0 exhausted, residue 1 entered: two residue gates.
        assert_eq!(gate_counters::RESIDUE.with(std::cell::Cell::get), 2);

        // Sticky error: after a missing-model-number error, further next()
        // calls neither yield nor re-run gates on later rows.
        let broken = vec![
            model(1, 0, 1),
            BioModelRow::new(BioRowSpan::new(1, 1).unwrap(), None),
        ];
        let mut mdl2 = Selection::default();
        mdl2.mdl = 2;
        gate_counters::reset();
        let mut err_iter = BioSelectedAtomIds::over_structure(
            &mdl2, &broken, &chains, &residues, &atoms, cif, unset, unset,
        );
        // Model 1 is gated and rejected; the missing-number model errors
        // on the first call.
        assert_eq!(
            err_iter.next(),
            Some(Err(super::BioRowTraverseError::Model(
                super::BioRowModelError::MissingModelNumber
            )))
        );
        assert_eq!(gate_counters::MODEL.with(std::cell::Cell::get), 2);
        // Sticky: no further yield and no further gate on any row.
        assert!(err_iter.next().is_none());
        assert!(err_iter.next().is_none());
        assert_eq!(gate_counters::MODEL.with(std::cell::Cell::get), 2);

        // Empty siblings: an accepted chain with an empty residue span
        // advances without any residue gate.
        let empty_chains = vec![chain_row(b"A", 0, 0)];
        let one_model = vec![model(1, 0, 1)];
        gate_counters::reset();
        let ids: Vec<u32> = BioSelectedAtomIds::over_structure(
            &all,
            &one_model,
            &empty_chains,
            &residues,
            &atoms,
            cif,
            unset,
            unset,
        )
        .map(|r| r.unwrap().value())
        .collect();
        assert!(ids.is_empty());
        assert_eq!(gate_counters::CHAIN.with(std::cell::Cell::get), 1);
        assert_eq!(gate_counters::RESIDUE.with(std::cell::Cell::get), 0);
    }
}
