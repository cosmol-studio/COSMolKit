//! BIO-CID C02: the single IO owner of the CID lexical closure, relocated
//! mechanically and compile-closed from `cosmolkit-bio/src/selection.rs`
//! (BIO-CID packet fixed contract). Bodies, verbatim Gemmi anchors and
//! two-axis markers are unchanged; only the module home moved. BIO keeps
//! the detached carriers (`SelectionCidList`, `SelectionSequenceId`,
//! `GemmiElementMask`) and the canonical element vocabulary/classification
//! (`GEMMI_EL_END`, `GEMMI_ELEMENT_NAMES`, `GEMMI_IS_METAL`,
//! `gemmi_find_element`), imported here through the C01/FIX boundary —
//! no copied table, no second classifier, no second predicate engine.

// Relocated items are exercised by this module's regressions today; the
// production caller is the CID decoder's parse_cid stages (BIO-CID C09+).
// Narrow module-level allowance: nothing else lives in this module.
#![allow(dead_code)]

use cosmolkit_bio::{
    BioSelectionData, BioSelectionParts, GEMMI_EL_END, GEMMI_IS_METAL, GemmiElementMask,
    SelectionAtomInequality, SelectionCidList, SelectionList, SelectionSequenceId,
    gemmi_find_element,
};

/// Canonical %.9g number owner (BIO-CID C20) consumed by the CID
/// serializer helpers below.
use crate::cif::format_cif_f64;

/// Canonical element DISPLAY-name owner (Gemmi `element_name(El)`,
/// elem.hpp:269-294: "He"/"Li"/"Cl", END's empty entry) for the CID
/// element-mask serialization in bio_selection_to_cid below — distinct
/// from the parser's UPPERCASE `element_uppercase_name` vocabulary.
use cosmolkit_bio::gemmi_element_name;

/// Relocated-test vocabulary only (BIO-CID C22-NAME): the parser's
/// UPPERCASE name table stays available to the relocated regressions and
/// local tests via `super::`; production serialization never uses it.
#[cfg(test)]
use cosmolkit_bio::GEMMI_ELEMENT_NAMES;

/// Gemmi `wrong_syntax` (select.cpp:14-26) as a typed error retaining
/// the full CID, the byte offset and the optional info note.
/// Exact `wrong_syntax` evidence carrier (public parse-error vocabulary,
/// BIO-CID C26): fields stay private; read access goes through the
/// accessors below so no construction surface leaks.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SelectionSyntaxError {
    /// The complete selection string, retained verbatim.
    pub(crate) cid: String,
    /// Byte offset of the failure site (0 suppresses the nearby text).
    pub(crate) pos: usize,
    /// Optional source info suffix; owned because `make_cid_list` composes
    /// its ` ('X' in a list)` note at runtime (select.cpp:53-55).
    pub(crate) info: Option<String>,
}

impl SelectionSyntaxError {
    /// The complete selection string, verbatim.
    #[must_use]
    pub fn cid(&self) -> &str {
        &self.cid
    }

    /// Byte offset of the failure site (0 suppresses the nearby text).
    #[must_use]
    pub const fn pos(&self) -> usize {
        self.pos
    }

    /// Optional source info suffix (e.g. ` ('X' in a list)`).
    pub fn info(&self) -> Option<&str> {
        self.info.as_deref()
    }
}

impl std::fmt::Display for SelectionSyntaxError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Gemmi❗✔️: inline GEMMI_COLD void wrong_syntax(const std::string& cid, size_t pos,
        // Gemmi✔️✔️:                                    const char* info=nullptr) {
        // Gemmi✔️✔️:   std::string msg = "Invalid selection syntax";
        // Gemmi✔️✔️:   if (info)
        // Gemmi✔️✔️:     msg += info;
        // Gemmi✔️✔️:   if (pos != 0)
        // Gemmi✔️✔️:     cat_to(msg, " near \"", cid.substr(pos, 8), '"');
        // Gemmi✔️✔️:   cat_to(msg, ": ", cid);
        // Gemmi✔️✔️:   fail(msg);
        // Gemmi✔️✔️: }
        // Behavior (❗ not byte-equivalent for split windows): the message
        // is "Invalid selection syntax" + optional info + (pos != 0)
        // ` near "<up to 8 bytes from pos>"` + ": " + the full cid. The
        // nearby window is taken on BYTE boundaries exactly like
        // std::string::substr; because a Rust String cannot hold invalid
        // UTF-8, a window splitting a multibyte character renders through
        // from_utf8_lossy (U+FFFD) — a presentation difference from the
        // source bytes. Exact byte fidelity, including split windows, is
        // provided by `source_message_bytes` (BIO-SEL-FIX Step 10); this
        // Display must not be marked byte-equivalent. pos beyond the end
        // yields an empty window without panicking (the source would throw
        // out_of_range there; no caller passes such a pos).
        // Complexity: O(cid) formatting work, one bounded window copy.
        write!(f, "Invalid selection syntax")?;
        if let Some(info) = self.info.as_deref() {
            write!(f, "{info}")?;
        }
        if self.pos != 0 {
            let tail = self.cid.as_bytes().get(self.pos..).map_or(&[][..], |t| t);
            let window = &tail[..tail.len().min(8)];
            write!(f, " near \"{}\"", String::from_utf8_lossy(window))?;
        }
        write!(f, ": {}", self.cid)
    }
}

impl std::error::Error for SelectionSyntaxError {}

impl SelectionSyntaxError {
    /// The exact `wrong_syntax` message BYTES (select.cpp:14-26):
    /// `"Invalid selection syntax"` + info + (pos != 0) ` near "` +
    /// `cid.substr(pos, 8)` + `"` + `": "` + cid — with the nearby window
    /// taken on raw byte boundaries exactly like `std::string::substr`
    /// (a window splitting a multibyte character contributes its exact
    /// cut bytes). `Display` is a UTF-8 presentation that renders such
    /// cut windows through `from_utf8_lossy`; it is NOT byte-equivalent
    /// and must not be marked as such. No public error-format policy is
    /// introduced: this accessor is private to the CID owner module.
    pub(crate) fn source_message_bytes(&self) -> Vec<u8> {
        let mut message =
            Vec::with_capacity(self.cid.len() + 24 + self.info.as_ref().map_or(0, String::len));
        message.extend_from_slice(b"Invalid selection syntax");
        if let Some(info) = &self.info {
            message.extend_from_slice(info.as_bytes());
        }
        if self.pos != 0 {
            message.extend_from_slice(b" near \"");
            let cid_bytes = self.cid.as_bytes();
            if self.pos <= cid_bytes.len() {
                let stop = (self.pos + 8).min(cid_bytes.len());
                message.extend_from_slice(&cid_bytes[self.pos..stop]);
            }
            message.push(b'"');
        }
        message.extend_from_slice(b": ");
        message.extend_from_slice(self.cid.as_bytes());
        message
    }
}

/// Source-shaped constructor mirroring the `wrong_syntax` call sites.
pub(crate) fn wrong_syntax(cid: &str, pos: usize, info: Option<&str>) -> SelectionSyntaxError {
    SelectionSyntaxError {
        cid: cid.to_string(),
        pos,
        info: info.map(str::to_string),
    }
}

/// Gemmi `determine_omitted_cid_fields` result (select.cpp:28-41):
/// how many leading `/mdl/chn/res` fields the CID omits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum OmittedCidFields {
    /// `/...`: the model field is present.
    Model = 0,
    /// A bare chain name.
    Chain = 1,
    /// A leading residue field.
    Residue = 2,
    /// A leading atom field.
    Atom = 3,
}

/// Gemmi `determine_omitted_cid_fields` (select.cpp:28-41).
pub(crate) fn determine_omitted_cid_fields(cid: &str) -> OmittedCidFields {
    // Gemmi✔️✔️: inline int determine_omitted_cid_fields(const std::string& cid) {
    // Gemmi✔️✔️:   if (cid[0] == '/')
    // Gemmi✔️✔️:     return 0; // model
    // Gemmi✔️✔️:   if (std::isdigit(cid[0]) || cid[0] == '.' || cid[0] == '(' || cid[0] == '-')
    // Gemmi✔️✔️:     return 2; // residue
    // Gemmi✔️✔️:   size_t sep = cid.find_first_of("/([:;");
    // Gemmi✔️✔️:   if (sep == std::string::npos || cid[sep] == '/' || cid[sep] == ';')
    // Gemmi✔️✔️:     return 1; // chain
    // Gemmi✔️✔️:   if (cid[sep] == '(')
    // Gemmi✔️✔️:     return 2; // residue
    // Gemmi✔️✔️:   return 3;  // atom
    // Gemmi✔️✔️: }
    // Behavior: '/' leads to model; a digit, '.', '(' or '-' first byte
    // means the residue field is first; otherwise the first of
    // `/([:;` decides — absent, '/' or ';' means a bare chain name, '('
    // means a residue name list, and '['/':'/';' position implies an atom
    // field. An empty cid reads as byte 0 (std::string::operator[] at
    // size() returns '\0'): no branch matches the leading tests and no
    // separator exists, so it classifies as a (empty) chain name, exactly
    // like the source.
    // Complexity: O(len) worst case for the separator scan; no allocation.
    let first = cid.as_bytes().first().copied().unwrap_or(0);
    if first == b'/' {
        return OmittedCidFields::Model;
    }
    if first.is_ascii_digit() || first == b'.' || first == b'(' || first == b'-' {
        return OmittedCidFields::Residue;
    }
    let sep = cid
        .as_bytes()
        .iter()
        .position(|b| matches!(b, b'/' | b'(' | b'[' | b':' | b';'));
    match sep {
        None => OmittedCidFields::Chain,
        Some(index) => match cid.as_bytes()[index] {
            b'/' | b';' => OmittedCidFields::Chain,
            b'(' => OmittedCidFields::Residue,
            _ => OmittedCidFields::Atom,
        },
    }
}

/// Gemmi `make_cid_list` (select.cpp:43-57): parse a `[!]/[*]/list` field
/// between byte positions, rejecting disallowed punctuation with the
/// source's composed info note.
pub(crate) fn make_cid_list(
    cid: &str,
    pos: usize,
    end: usize,
    disallowed_chars: &str,
) -> Result<SelectionCidList, SelectionSyntaxError> {
    // Gemmi✔️✔️: inline Selection::List make_cid_list(const std::string& cid, size_t pos, size_t end,
    // Gemmi✔️✔️:                                      const char* disallowed_chars="-[]()!/*.:;") {
    // Gemmi✔️✔️:   Selection::List list;
    // Gemmi✔️✔️:   list.all = (cid[pos] == '*');
    // Gemmi✔️✔️:   list.inverted = (cid[pos] == '!');
    // Gemmi✔️✔️:   if (list.all || list.inverted)
    // Gemmi✔️✔️:     ++pos;
    // Gemmi✔️✔️:   list.list = cid.substr(pos, end - pos);
    // Gemmi✔️✔️:   // if a list have punctuation other than ',' something must be wrong
    // Gemmi✔️✔️:   size_t idx = list.list.find_first_of(disallowed_chars);
    // Gemmi✔️✔️:   if (idx != std::string::npos)
    // Gemmi✔️✔️:     wrong_syntax(cid, pos + idx, cat(" ('", list.list[idx], "' in a list)").c_str());
    // Gemmi✔️✔️:   return list;
    // Gemmi✔️✔️: }
    // Behavior: the byte at `pos` is inspected regardless of `end` (the
    // source indexes cid[pos] directly; pos == len reads the '\0' that
    // std::string::operator[](size()) returns); '*' or '!' consumes
    // exactly one leading byte and sets all/inverted; the member payload
    // is the BYTE range reproducing std::string::substr exactly — clamped
    // to the string end when end exceeds it, and returning the remainder
    // [start, len] when end < start because the source's `end - pos`
    // size_t subtraction wraps to a huge count that substr clamps (an
    // impossible argument from the pinned parse_cid callers, retained for
    // helper-level substr fidelity). Multibyte UTF-8 member names are
    // preserved byte-for-byte; a byte-boundary cut inside a multibyte
    // character (also impossible from the callers — all scanning bytes
    // are ASCII) yields the exact cut bytes, never a lossy conversion or
    // an empty substitute. The first member byte in the disallowed set
    // fails with the advanced position and the composed ` ('X' in a list)`
    // note (the offender byte is always ASCII because every disallowed
    // set used by the source is ASCII). start > len reproduces the
    // source's std::out_of_range throw as a panic: it is unreachable from
    // the pinned callers and no typed syntax error exists for it in the
    // source (documented in the receipt audit).
    // Complexity: O(field bytes) single scan plus one byte-vector copy of
    // the members, matching the source's substr and find_first_of.
    let bytes = cid.as_bytes();
    let all = bytes.get(pos) == Some(&b'*');
    let inverted = bytes.get(pos) == Some(&b'!');
    let start = if all || inverted { pos + 1 } else { pos };
    if start > bytes.len() {
        panic!(
            "make_cid_list: position {start} beyond cid length {} \
             (std::string::substr out_of_range; unreachable from parse_cid callers)",
            bytes.len()
        );
    }
    // substr(rpos, rcount): rcount clamps to the remaining length, and a
    // wrapped (end < start) count likewise selects [start, len).
    let stop = if end >= start {
        end.min(bytes.len())
    } else {
        bytes.len()
    };
    let list_bytes: &[u8] = if start >= stop {
        &[]
    } else {
        &bytes[start..stop]
    };
    if let Some(idx) = list_bytes
        .iter()
        .position(|byte| disallowed_chars.as_bytes().contains(byte))
    {
        let offender = list_bytes[idx];
        debug_assert!(offender.is_ascii());
        return Err(wrong_syntax(
            cid,
            start + idx,
            Some(&format!(" ('{}' in a list)", offender as char)),
        ));
    }
    Ok(SelectionCidList {
        all,
        inverted,
        list: list_bytes.to_vec(),
    })
}

/// Typed bounds failure for a CID sequence number that the source's
/// `strtol` → `int` path cannot represent without C implementation-defined
/// truncation; the packet mandates an error instead of wrapping.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SelectionSeqidRangeError {
    /// The offending numeric text exactly as it appeared.
    pub(crate) text: String,
}

impl std::fmt::Display for SelectionSeqidRangeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "CID sequence number out of representable range: {}",
            self.text
        )
    }
}

impl std::error::Error for SelectionSeqidRangeError {}

/// Gemmi `parse_cid_seqid` (select.cpp:95-115): parse `s1.i1`-shaped
/// sequence fields, advancing the shared byte position.
pub(crate) fn parse_cid_seqid(
    cid: &str,
    pos: &mut usize,
    default_seqnum: i32,
) -> Result<SelectionSequenceId, SelectionSeqidRangeError> {
    // Gemmi✔️✔️: inline Selection::SequenceId parse_cid_seqid(const std::string& cid, size_t& pos,
    // Gemmi✔️✔️:                                            int default_seqnum) {
    // Gemmi✔️✔️:   size_t initial_pos = pos;
    // Gemmi✔️✔️:   int seqnum = default_seqnum;
    // Gemmi✔️✔️:   char icode = ' ';
    // Gemmi✔️✔️:   if (cid[pos] == '*') {
    // Gemmi✔️✔️:     ++pos;
    // Gemmi✔️✔️:     icode = '*';
    // Gemmi✔️✔️:   } else if (std::isdigit(cid[pos]) || cid[pos] == '-') {
    // Gemmi✔️✔️:     char* endptr;
    // Gemmi✔️✔️:     seqnum = std::strtol(&cid[pos], &endptr, 10);
    // Gemmi✔️✔️:     pos = endptr - &cid[0];
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   if (cid[pos] == '.')
    // Gemmi✔️✔️:     ++pos;
    // Gemmi✔️✔️:   if (initial_pos != pos && (std::isalpha(cid[pos]) || cid[pos] == '*'))
    // Gemmi✔️✔️:     icode = cid[pos++];
    // Gemmi✔️✔️:   return {seqnum, icode};
    // Gemmi✔️✔️: }
    // Behavior: '*' takes the wildcard icode and leaves seqnum at the
    // caller's default (INT_MIN for from, INT_MAX for to). A digit or '-'
    // runs a strtol-shaped conversion: optional sign then at least one
    // digit; with no digit after the sign strtol performs no conversion,
    // returns 0 and leaves the end pointer at the start — so the seqnum
    // becomes 0 and the position does not advance. The optional '.' is
    // consumed unconditionally; an icode letter (or '*') is consumed only
    // when something before it was consumed. Packet-mandated divergence:
    // a magnitude outside i32 raises a typed range error instead of the
    // source's implementation-defined long→int truncation (strtol itself
    // clamps to LONG_MIN/LONG_MAX with errno set, which the source ignores).
    // cid[pos] at the string end reads as '\0' exactly like
    // std::string::operator[] at size().
    // Complexity: O(digits) single scan with no allocation except the
    // range-error text.
    let bytes = cid.as_bytes();
    let initial_pos = *pos;
    let mut seqnum = default_seqnum;
    let mut icode = b' ';
    let at = |index: usize| -> u8 { bytes.get(index).copied().unwrap_or(0) };
    let current = at(*pos);
    if current == b'*' {
        *pos += 1;
        icode = b'*';
    } else if current.is_ascii_digit() || current == b'-' {
        // C06: the single checked integer-prefix owner (C05) supplies the
        // strtol conversion — no ad-hoc accumulation here. The caller
        // guard means the scan starts on a digit or '-'; whitespace is
        // unreachable in this branch but the owner keeps full strtol
        // semantics for the model stage's reuse.
        let scan = scan_cid_int_prefix(&bytes[*pos..]);
        if !scan.converted {
            // strtol no-conversion (e.g. a lone '-'): value 0, position
            // unchanged (endptr == start).
            seqnum = 0;
        } else {
            if scan.value < i64::from(i32::MIN) || scan.value > i64::from(i32::MAX) {
                return Err(SelectionSeqidRangeError {
                    text: cid[*pos..*pos + scan.consumed].to_string(),
                });
            }
            seqnum = scan.value as i32;
            *pos += scan.consumed;
        }
    }
    if at(*pos) == b'.' {
        *pos += 1;
    }
    let tail = at(*pos);
    if initial_pos != *pos && (tail.is_ascii_alphabetic() || tail == b'*') {
        icode = tail;
        *pos += 1;
    }
    Ok(SelectionSequenceId { seqnum, icode })
}

/// Gemmi `parse_cid_elements` (select.cpp:59-93): parse an `[el,el...]`
/// mask after '[', up to and including the closing ']'. Relocated to IO
/// (BIO-CID C02); the element ordinals, name lookup and metals classifier
/// are imported from the canonical BIO vocabulary — nothing is copied.
pub(crate) fn parse_cid_elements(
    cid: &str,
    pos: usize,
) -> Result<Option<GemmiElementMask>, SelectionSyntaxError> {
    // Gemmi✔️✔️: inline void parse_cid_elements(const std::string& cid, size_t pos,
    // Gemmi✔️✔️:                                std::vector<char>& elements) {
    // Gemmi✔️✔️:   elements.clear();  // just in case
    // Gemmi✔️✔️:   if (cid[pos] == '*')
    // Gemmi✔️✔️:     return;
    // Gemmi✔️✔️:   bool inverted = false;
    // Gemmi✔️✔️:   if (cid[pos] == '!') {
    // Gemmi✔️✔️:     inverted = true;
    // Gemmi✔️✔️:     ++pos;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   elements.resize((size_t)El::END, char(inverted));
    // Gemmi✔️✔️:   for (;;) {
    // Gemmi✔️✔️:     size_t sep = cid.find_first_of(",]", pos);
    // Gemmi✔️✔️:     if (sep == pos || sep > pos + 2) {
    // Gemmi✔️✔️:       if (sep == pos + 6 && cid.compare(pos, 6, "metals", 6) == 0) {
    // Gemmi✔️✔️:         for (size_t i = 0; i < elements.size(); ++i)
    // Gemmi✔️✔️:           if (is_metal(static_cast<El>(i)))
    // Gemmi✔️✔️:             elements[i] = char(!inverted);
    // Gemmi✔️✔️:       } else if (sep == pos + 9 && cid.compare(pos, 9, "nonmetals", 9) == 0) {
    // Gemmi✔️✔️:         for (size_t i = 0; i < elements.size(); ++i)
    // Gemmi✔️✔️:           if (!is_metal(static_cast<El>(i)))
    // Gemmi✔️✔️:             elements[i] = char(!inverted);
    // Gemmi✔️✔️:       } else {
    // Gemmi✔️✔️:         wrong_syntax(cid, 0, " in [...]");
    // Gemmi✔️✔️:       }
    // Gemmi✔️✔️:     } else {
    // Gemmi✔️✔️:       char elem_str[2] = {cid[pos], sep > pos+1 ? cid[pos+1] : '\0'};
    // Gemmi✔️:       Element el = find_element(elem_str);
    // Gemmi✔️✔️:       if (el == El::X && (alpha_up(elem_str[0]) != 'X' || elem_str[1] != '\0'))
    // Gemmi✔️✔️:         wrong_syntax(cid, 0, " (invalid element in [...])");
    // Gemmi✔️✔️:       elements[el.ordinal()] = char(!inverted);
    // Gemmi✔️:     }
    // Gemmi✔️✔️:     if (cid[sep] == ']')
    // Gemmi✔️✔️:       break;
    // Gemmi✔️✔️:     pos = sep + 1;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: }
    // Behavior: '*' returns no mask (the source leaves the vector empty);
    // '!' pre-fills every ordinal with inverted and selected ordinals flip
    // to !inverted. Members are one or two bytes (a longer token must be
    // exactly `metals`/`nonmetals`, case-sensitive); a two-byte member goes
    // through find_element with its `& ~0x20` folding (so lowercase input
    // resolves like the source); an X result is accepted ONLY for a
    // literal single-letter X/x spelling. Mask entries use Gemmi El
    // ordinals (X=0, H..Og=1..118, D=119) with the pinned is_metal table.
    // Divergence note: without any ',' or ']' the source's pos = npos+1
    // wraps to 0 and loops forever; callers guarantee a ']' (the earlier
    // "no matching ']'" check), and this port fails closed with the
    // neighboring " in [...]" syntax error instead of looping.
    // Complexity: O(members x 120) worst case, one fixed-size mask.
    let bytes = cid.as_bytes();
    let at = |index: usize| -> u8 { bytes.get(index).copied().unwrap_or(0) };
    if at(pos) == b'*' {
        return Ok(None);
    }
    let mut inverted = false;
    let mut pos = pos;
    if at(pos) == b'!' {
        inverted = true;
        pos += 1;
    }
    let mut elements = [inverted; GEMMI_EL_END];
    loop {
        let sep = bytes[pos..]
            .iter()
            .position(|byte| *byte == b',' || *byte == b']')
            .map_or(None, |offset| Some(pos + offset));
        let Some(sep) = sep else {
            return Err(wrong_syntax(cid, 0, Some(" in [...]")));
        };
        if sep == pos || sep > pos + 2 {
            let token = &bytes[pos..sep];
            if sep == pos + 6 && token == b"metals" {
                for index in 0..GEMMI_EL_END {
                    if GEMMI_IS_METAL[index] {
                        elements[index] = !inverted;
                    }
                }
            } else if sep == pos + 9 && token == b"nonmetals" {
                for index in 0..GEMMI_EL_END {
                    if !GEMMI_IS_METAL[index] {
                        elements[index] = !inverted;
                    }
                }
            } else {
                return Err(wrong_syntax(cid, 0, Some(" in [...]")));
            }
        } else {
            let elem_first = at(pos);
            let elem_second = if sep > pos + 1 { at(pos + 1) } else { 0 };
            let ordinal = gemmi_find_element(elem_first, elem_second);
            if ordinal == 0 && ((elem_first & !0x20) != b'X' || elem_second != 0) {
                return Err(wrong_syntax(cid, 0, Some(" (invalid element in [...])")));
            }
            elements[ordinal as usize] = !inverted;
        }
        if at(sep) == b']' {
            break;
        }
        pos = sep + 1;
    }
    Ok(Some(elements))
}

/// Gemmi `has_inequality` (select.cpp:143-148): does the byte range
/// contain a relation token? The numeric payload is deliberately not
/// parsed here (decision queue: single approved numeric owner).
pub(crate) fn has_inequality(cid: &str, start: usize, end: usize) -> bool {
    // Gemmi✔️✔️: inline bool has_inequality(const std::string& cid, size_t start, size_t end) {
    // Gemmi✔️✔️:   for (size_t i = start; i < end; ++i)
    // Gemmi✔️✔️:     if (cid[i] == '<' || cid[i] == '=' || cid[i] == '>')
    // Gemmi✔️✔️:       return true;
    // Gemmi✔️✔️:   return false;
    // Gemmi✔️✔️: }
    // Behavior: plain byte scan over [start, end); '=' anywhere counts,
    // including inside what would be a malformed token — recognition only,
    // no consumption and no numeric parsing.
    // Complexity: O(end - start) bytes, no allocation.
    cid.as_bytes()
        .get(start..end)
        .is_some_and(|window| window.iter().any(|byte| matches!(byte, b'<' | b'=' | b'>')))
}

/// Result of the one shared checked decimal integer-prefix scan
/// (BIO-CID C05). The CID model and sequence stages both call
/// `std::strtol(..., 10)` on a suffix of the selection string.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct CidIntScan {
    /// strtol-shaped value on the pinned LP64 platform (`long` = i64):
    /// `0` when no conversion happened; `LONG_MIN`/`LONG_MAX` saturation
    /// when the digit magnitude exceeds the long range (strtol also sets
    /// errno, which both source call sites ignore). The packet-mandated
    /// typed i32 bounds check belongs to the callers — never a C
    /// long→int truncation.
    pub value: i64,
    /// Full consumed extent in bytes: leading is_space bytes, one
    /// optional sign, then EVERY decimal digit (strtol consumes the
    /// complete digit run even when saturating).
    pub consumed: usize,
    /// Whether strtol performed a conversion. `false` means endptr was
    /// left at the ORIGINAL scan start (consumed = 0), exactly like C
    /// strtol with no digits — whitespace/sign alone never converts.
    pub converted: bool,
}

/// The one checked decimal integer-prefix owner shared by the CID model
/// and sequence conversions (BIO-CID C05). Byte accumulation saturates
/// instead of overflowing, so arbitrarily long digit runs never wrap the
/// internal arithmetic while the consumed extent still spans every digit.
pub(crate) fn scan_cid_int_prefix(field: &[u8]) -> CidIntScan {
    // Gemmi✔️✔️:     seqnum = std::strtol(&cid[pos], &endptr, 10);
    // Gemmi✔️✔️:     sel.mdl = std::strtol(&cid[1], &endptr, 10);
    //
    // Behavior (strtol semantics, base 10, pinned LP64 `long` = i64):
    // optional leading is_space bytes (the pinned Gemmi atox.hpp table
    // 9-13 + 32, identical to the C isspace set for ASCII), one optional
    // '+'/'-' sign AFTER the whitespace, then a decimal digit run. With no
    // digit after whitespace/sign the conversion fails: value 0 and the
    // end pointer stays at the ORIGINAL start (not after the whitespace).
    // A digit run longer than the long range saturates to LONG_MAX/LONG_MIN
    // while STILL consuming every digit (errno is set by strtol and ignored
    // by both call sites). Exact 2^63 with '-' is LONG_MIN and therefore
    // NOT saturation; exact 2^63-1 is LONG_MAX. The typed i32 range error
    // is applied by the callers (C06/C09) from this scan — the owner never
    // performs the source's implementation-defined long→int truncation.
    // Complexity: single pass over the consumed bytes; no allocation.
    let mut index = 0usize;
    while index < field.len() && crate::bio_pdb::gemmi_is_space(field[index]) {
        index += 1;
    }
    let mut negative = false;
    if index < field.len() && (field[index] == b'+' || field[index] == b'-') {
        negative = field[index] == b'-';
        index += 1;
    }
    let digits_start = index;
    let long_min_magnitude = 1u128 << 63;
    let long_max_magnitude = (1u128 << 63) - 1;
    let mut magnitude: u128 = 0;
    let mut saturated = false;
    while index < field.len() && field[index].is_ascii_digit() {
        if !saturated {
            let digit = u128::from(field[index] - b'0');
            magnitude = match magnitude
                .checked_mul(10)
                .and_then(|scaled| scaled.checked_add(digit))
            {
                Some(value) => value,
                // Beyond u128: clamp; the extent scan continues over digits.
                None => {
                    saturated = true;
                    magnitude
                }
            };
            let bound = if negative {
                long_min_magnitude
            } else {
                long_max_magnitude
            };
            if magnitude > bound {
                saturated = true;
            }
        }
        index += 1;
    }
    if index == digits_start {
        // strtol no-conversion: value 0, end pointer back at the start.
        return CidIntScan {
            value: 0,
            consumed: 0,
            converted: false,
        };
    }
    let value = if saturated {
        if negative { i64::MIN } else { i64::MAX }
    } else if negative {
        // magnitude <= 2^63 by the bound check above; negation fits i64.
        (magnitude as i128).wrapping_neg() as i64
    } else {
        magnitude as i64
    };
    CidIntScan {
        value,
        consumed: index,
        converted: true,
    }
}

/// Parses a `q…`/`b…` occupancy/B-factor inequality payload (BIO-CID C08)
/// through the shared bio_numeric C-string owner, filling the canonical
/// BIO `SelectionAtomInequality` carrier. On success `pos` is left at the
/// consumed boundary (equal to `end`).
pub(crate) fn parse_atom_inequality(
    cid: &str,
    pos: &mut usize,
    end: usize,
) -> Result<SelectionAtomInequality, SelectionSyntaxError> {
    // Gemmi✔️✔️: inline Selection::AtomInequality parse_atom_inequality(const std::string& cid,
    // Gemmi✔️✔️:                                                        size_t pos, size_t end) {
    // Gemmi✔️✔️:   Selection::AtomInequality r;
    // Gemmi✔️✔️:   if (cid[pos] != 'q' && cid[pos] != 'b')
    // Gemmi✔️✔️:     wrong_syntax(cid, pos);
    // Gemmi✔️✔️:   r.property = cid[pos];
    // Gemmi✔️✔️:   ++pos;
    // Gemmi✔️✔️:   while (cid[pos] == ' ')
    // Gemmi✔️✔️:     ++pos;
    // Gemmi✔️✔️:   if (cid[pos] == '<')
    // Gemmi✔️✔️:     r.relation = -1;
    // Gemmi✔️✔️:   else if (cid[pos] == '>')
    // Gemmi✔️✔️:     r.relation = 1;
    // Gemmi✔️✔️:   else if (cid[pos] == '=')
    // Gemmi✔️✔️:     r.relation = 0;
    // Gemmi✔️✔️:   else
    // Gemmi✔️✔️:     wrong_syntax(cid, pos);
    // Gemmi✔️✔️:   ++pos;
    // Gemmi✔️✔️:   auto result = fast_from_chars(cid.c_str() + pos, r.value);
    // Gemmi✔️✔️:   if (result.ec != std::errc())
    // Gemmi✔️✔️:     wrong_syntax(cid, pos, " (expected number)");
    // Gemmi✔️✔️:   pos = size_t(result.ptr - cid.c_str());
    // Gemmi✔️✔️:   while (cid[pos] == ' ')
    // Gemmi✔️✔️:     ++pos;
    // Gemmi✔️✔️:   if (pos != end)
    // Gemmi✔️✔️:     wrong_syntax(cid, pos);
    // Gemmi✔️✔️:   return r;
    // Gemmi✔️✔️: }
    //
    // Behavior review: the two space loops read `cid[pos]` WITHOUT the
    // end bound (std::string::operator[] at size() is '\0'), so trailing
    // spaces beyond `end` are consumed and then rejected by the final
    // `pos != end` check; the mirror `at()` yields 0 past the last byte.
    // The numeric tail is the C-string `fast_from_chars` overload (atof.hpp
    // :24-30): leading is_space skip, one optional '+', then fast_float
    // from_chars to strlen — all inside the shared bio_numeric owner, and
    // its `consumed` is relative to the slice passed in, so `pos` advances
    // to the absolute `result.ptr` position. `ec != errc()` covers BOTH
    // invalid_argument and result_out_of_range: a syntactically valid but
    // out-of-range literal (1e999) is rejected with the same note even
    // though fast_float assigns a saturated value; only literal
    // inf/infinity/nan spellings (with optional sign) parse successfully.
    // Complexity: O(payload) single pass, no allocation.
    let bytes = cid.as_bytes();
    let at = |index: usize| -> u8 { bytes.get(index).copied().unwrap_or(0) };
    let mut result = SelectionAtomInequality {
        property: 0,
        relation: 0,
        value: 0.0,
    };
    let start = *pos;
    let mut scan = start;
    if at(scan) != b'q' && at(scan) != b'b' {
        return Err(wrong_syntax(cid, scan, None));
    }
    result.property = at(scan);
    scan += 1;
    while at(scan) == b' ' {
        scan += 1;
    }
    result.relation = match at(scan) {
        b'<' => -1,
        b'>' => 1,
        b'=' => 0,
        _ => return Err(wrong_syntax(cid, scan, None)),
    };
    scan += 1;
    let outcome = crate::bio_numeric::fast_from_chars_cstring_typed(&bytes[scan..], result.value);
    if outcome.category.is_some() {
        return Err(wrong_syntax(cid, scan, Some(" (expected number)")));
    }
    result.value = outcome.value;
    scan += outcome.consumed;
    while at(scan) == b' ' {
        scan += 1;
    }
    if scan != end {
        return Err(wrong_syntax(cid, scan, None));
    }
    *pos = scan;
    Ok(result)
}

/// Gemmi parse_cid preamble bypass (BIO-CID C09): an empty selection or
/// a lone `*` selects everything — no stage runs.
pub(crate) fn cid_selection_bypassed(cid: &str) -> bool {
    // Gemmi✔️✔️:   if (cid.empty() || (cid.size() == 1 && cid[0] == '*'))
    // Gemmi✔️✔️:     return;
    cid.is_empty() || (cid.len() == 1 && cid.as_bytes()[0] == b'*')
}

/// CID model stage (BIO-CID C09). Runs only when the leading field is the
/// model (`determine_omitted_cid_fields` == 0, i.e. the string starts with
/// `/`). Returns the model number and the stage's separator position for
/// the chain stage (`None` mirrors std::string::npos). No public decoder
/// is exposed until all stages are composed.
pub(crate) fn parse_cid_model_stage(
    cid: &str,
    omit: i32,
    semi: Option<usize>,
) -> Result<(i32, Option<usize>), BioSelectionParseError> {
    // Gemmi✔️✔️:   int omit = determine_omitted_cid_fields(cid);
    // Gemmi✔️✔️:   size_t sep = 0;
    // Gemmi✔️✔️:   size_t semi = cid.find(';');
    // Gemmi✔️✔️:   // model
    // Gemmi✔️✔️:   if (omit == 0) {
    // Gemmi✔️✔️:     sep = std::min(cid.find('/', 1), semi);
    // Gemmi✔️✔️:     if (sep != 1 && cid[1] != '*') {
    // Gemmi✔️✔️:       char* endptr;
    // Gemmi✔️✔️:       sel.mdl = std::strtol(&cid[1], &endptr, 10);
    // Gemmi✔️✔️:       size_t end_pos = endptr - &cid[0];
    // Gemmi✔️✔️:       if (end_pos != sep && end_pos != cid.size())
    // Gemmi✔️✔️:         wrong_syntax(cid, 0, " (at model number)");
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:   }
    //
    // Behavior review: omit != 0 leaves mdl at the Selection default 0
    // ("0 = all", select.hpp:92) and sep at the source's initial 0 for the
    // chain stage. omit == 0 means cid[0] == '/': sep is min(first '/'
    // from index 1, semi) with npos as None. sep == 1 (empty model field,
    // "//...") or cid[1] == '*' ("/*...") bypasses the conversion — mdl
    // stays 0 = all models. Otherwise strtol runs from cid[1] through the
    // shared C05 owner: a no-conversion tail yields mdl 0 and end_pos 1
    // (e.g. "/" alone: end_pos == cid.size()); a digit run outside the
    // i32 range is the packet-mandated typed range error instead of the
    // source's implementation-defined long->int truncation; end_pos must
    // equal sep or cid.size(), else the source's fixed " (at model
    // number)" note at position 0. Complexity: one scan over the model
    // bytes, no allocation.
    if omit != 0 {
        return Ok((0, Some(0)));
    }
    let bytes = cid.as_bytes();
    let slash = bytes[1..]
        .iter()
        .position(|byte| *byte == b'/')
        .map(|index| index + 1);
    let sep = match (slash, semi) {
        (Some(a), Some(b)) => Some(a.min(b)),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    };
    if sep != Some(1) && bytes.get(1) != Some(&b'*') {
        let scan = scan_cid_int_prefix(&bytes[1..]);
        let mdl = if !scan.converted {
            0
        } else {
            if scan.value < i64::from(i32::MIN) || scan.value > i64::from(i32::MAX) {
                return Err(BioSelectionParseError::SeqidRange(
                    SelectionSeqidRangeError {
                        text: cid[1..1 + scan.consumed].to_string(),
                    },
                ));
            }
            scan.value as i32
        };
        let end_pos = if scan.converted { 1 + scan.consumed } else { 1 };
        if Some(end_pos) != sep && end_pos != cid.len() {
            return Err(BioSelectionParseError::Syntax(wrong_syntax(
                cid,
                0,
                Some(" (at model number)"),
            )));
        }
        return Ok((mdl, sep));
    }
    Ok((0, sep))
}

/// CID chain stage (BIO-CID C10). Consumes the model stage's separator
/// and yields the chain list plus the residue stage's separator.
pub(crate) fn parse_cid_chain_stage(
    cid: &str,
    omit: i32,
    semi: Option<usize>,
    model_sep: Option<usize>,
) -> Result<(SelectionCidList, Option<usize>), BioSelectionParseError> {
    // Gemmi✔️✔️:   // chain
    // Gemmi✔️✔️:   if (omit <= 1 && sep < semi) {
    // Gemmi✔️✔️:     size_t pos = (sep == 0 ? 0 : sep + 1);
    // Gemmi✔️✔️:     sep = std::min(cid.find('/', pos), semi);
    // Gemmi✔️✔️:     // These characters are not really disallowed, but are unexpected.
    // Gemmi✔️✔️:     // "-" is expected, it's in chain IDs in bioassembly files from RCSB.
    // Gemmi✔️✔️:     const char* disallowed_chars = "[]()!/*.:;";
    // Gemmi✔️✔️:     sel.chain_ids = make_cid_list(cid, pos, sep, disallowed_chars);
    // Gemmi✔️✔️:   }
    //
    // Behavior review: `sep` is the model stage's output separator; on a
    // false guard the pinned parse_cid leaves sep UNCHANGED (identity-on-
    // skip), so the residue stage still sees Some(0) for residue-first /
    // atom-first CIDs whose chain stage was skipped and Some(2) for "/1;"
    // — the returned separator is the incoming model_sep, never None
    // unless the model stage itself ended at npos (None mirrors npos).
    // pos restarts after the separator (or 0 when the model stage never
    // ran, sep == 0). The new separator is min(first '/' from pos, semi)
    // for the residue stage. The chain disallowed set omits '-': literal
    // '-' chain names from RCSB bioassembly files stay legal (source
    // comment above). Empty member lists are retained verbatim: an empty
    // chain name is a valid blank chain, never an error. Complexity: one
    // scan plus the C04 list builder; no allocation beyond the list.
    if omit <= 1 && model_sep.unwrap_or(usize::MAX) < semi.unwrap_or(usize::MAX) {
        let pos = match model_sep {
            Some(0) | None => 0,
            Some(sep) => sep + 1,
        };
        let bytes = cid.as_bytes();
        let slash = bytes[pos..]
            .iter()
            .position(|byte| *byte == b'/')
            .map(|index| index + pos);
        let sep = match (slash, semi) {
            (Some(a), Some(b)) => Some(a.min(b)),
            (Some(a), None) => Some(a),
            (None, Some(b)) => Some(b),
            (None, None) => None,
        };
        let chain_ids = make_cid_list(cid, pos, sep.unwrap_or(cid.len()), "[]()!/*.:;")
            .map_err(BioSelectionParseError::Syntax)?;
        return Ok((chain_ids, sep));
    }
    Ok((
        SelectionCidList {
            all: true,
            inverted: false,
            list: Vec::new(),
        },
        model_sep,
    ))
}

/// CID residue stage (BIO-CID C11). Consumes the chain stage's
/// separator and yields both seqid bounds, the residue name list and the
/// atom stage's separator.
pub(crate) fn parse_cid_residue_stage(
    cid: &str,
    omit: i32,
    semi: Option<usize>,
    chain_sep: Option<usize>,
) -> Result<
    (
        SelectionSequenceId,
        SelectionSequenceId,
        SelectionCidList,
        Option<usize>,
    ),
    BioSelectionParseError,
> {
    // Gemmi✔️✔️:   // residue; MMDB CID syntax: s1.i1-s2.i2 or *(res).ic
    // Gemmi✔️✔️:   // In gemmi both 14.a and 14a are accepted.
    // Gemmi✔️✔️:   // *(ALA). and *(ALA) and (ALA). can be used instead of (ALA) for
    // Gemmi✔️✔️:   // compatibility with MMDB.
    // Gemmi✔️✔️:   if (omit <= 2 && sep < semi) {
    // Gemmi✔️✔️:     size_t pos = (sep == 0 ? 0 : sep + 1);
    // Gemmi✔️✔️:     if (cid[pos] != '(')
    // Gemmi✔️✔️:       sel.from_seqid = parse_cid_seqid(cid, pos, INT_MIN);
    // Gemmi✔️✔️:     if (cid[pos] == '(') {
    // Gemmi✔️✔️:       ++pos;
    // Gemmi✔️✔️:       size_t right_br = cid.find(')', pos);
    // Gemmi✔️✔️:       sel.residue_names = make_cid_list(cid, pos, right_br);
    // Gemmi✔️✔️:       pos = right_br + 1;
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:     // allow "(RES)." and "(RES).*" and "(RES)*"
    // Gemmi✔️✔️:     if (cid[pos] == '.')
    // Gemmi✔️✔️:       ++pos;
    // Gemmi✔️✔️:     if (cid[pos] == '*')
    // Gemmi✔️✔️:       ++pos;
    // Gemmi✔️✔️:     if (cid[pos] == '-') {
    // Gemmi✔️✔️:       ++pos;
    // Gemmi✔️✔️:       sel.to_seqid = parse_cid_seqid(cid, pos, INT_MAX);
    // Gemmi✔️✔️:     } else if (sel.from_seqid.seqnum != INT_MIN) {
    // Gemmi✔️✔️:       sel.to_seqid = sel.from_seqid;
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:     sep = pos;
    // Gemmi✔️✔️:     if (cid[sep] != '/' && cid[sep] != ';' && cid[sep] != '\0')
    // Gemmi✔️✔️:       wrong_syntax(cid, 0, " (at residue)");
    // Gemmi✔️✔️:   }
    //
    // Behavior review: the name-list test RE-READS cid[pos] AFTER the
    // from_seqid parse advanced pos ("14(?)" enters the name branch from
    // the advanced position, exactly like the source). A missing ')'
    // gives right_br = npos; make_cid_list then reads to the string end,
    // and pos = right_br + 1 WRAPS to 0 — preserved via wrapping_add — so
    // the dot/star/dash checks restart from the string start before the
    // terminal check fires. A point selection (seqnum != INT_MIN, no
    // '-'-bound) copies from_seqid into to_seqid. On a false guard the
    // pinned parse_cid leaves sep unchanged, so the incoming chain_sep is
    // returned as-is (C10-SEP identity-on-skip). Markers ✔️✔️: one scan
    // plus the C04/C06 owners; no allocation beyond the name list.
    let mut from_seqid = SelectionSequenceId {
        seqnum: i32::MIN,
        icode: b'*',
    };
    let mut to_seqid = SelectionSequenceId {
        seqnum: i32::MAX,
        icode: b'*',
    };
    if omit <= 2 && chain_sep.unwrap_or(usize::MAX) < semi.unwrap_or(usize::MAX) {
        let bytes = cid.as_bytes();
        let at = |index: usize| bytes.get(index).copied().unwrap_or(0);
        let mut pos = match chain_sep {
            Some(0) | None => 0,
            Some(sep) => sep + 1,
        };
        if at(pos) != b'(' {
            from_seqid = parse_cid_seqid(cid, &mut pos, i32::MIN)
                .map_err(BioSelectionParseError::SeqidRange)?;
        }
        let mut residue_names = SelectionCidList {
            all: true,
            inverted: false,
            list: Vec::new(),
        };
        if at(pos) == b'(' {
            pos += 1;
            let right_br = bytes[pos..]
                .iter()
                .position(|byte| *byte == b')')
                .map(|index| index + pos);
            residue_names = make_cid_list(cid, pos, right_br.unwrap_or(cid.len()), "-[]()!/*.:;")
                .map_err(BioSelectionParseError::Syntax)?;
            pos = right_br.unwrap_or(usize::MAX).wrapping_add(1);
        }
        if at(pos) == b'.' {
            pos += 1;
        }
        if at(pos) == b'*' {
            pos += 1;
        }
        if at(pos) == b'-' {
            pos += 1;
            to_seqid = parse_cid_seqid(cid, &mut pos, i32::MAX)
                .map_err(BioSelectionParseError::SeqidRange)?;
        } else if from_seqid.seqnum != i32::MIN {
            to_seqid = from_seqid;
        }
        let sep = pos;
        if at(sep) != b'/' && at(sep) != b';' && at(sep) != 0 {
            return Err(BioSelectionParseError::Syntax(wrong_syntax(
                cid,
                0,
                Some(" (at residue)"),
            )));
        }
        return Ok((from_seqid, to_seqid, residue_names, Some(sep)));
    }
    Ok((
        from_seqid,
        to_seqid,
        SelectionCidList {
            all: true,
            inverted: false,
            list: Vec::new(),
        },
        chain_sep,
    ))
}

/// CID atom stage (BIO-CID C12). Consumes the residue stage's separator
/// and yields atom names, element mask, altloc list, and whether the
/// pinned source's early return fired (no element/altloc parsing after).
pub(crate) fn parse_cid_atom_stage(
    cid: &str,
    semi: Option<usize>,
    residue_sep: Option<usize>,
) -> Result<
    (
        SelectionCidList,
        Option<GemmiElementMask>,
        SelectionCidList,
        bool,
    ),
    BioSelectionParseError,
> {
    // Gemmi✔️✔️:   // atom;  at[el]:aloc
    // Gemmi✔️✔️:   if (sep < std::min(cid.size(), semi)) {
    // Gemmi✔️✔️:     size_t pos = (sep == 0 ? 0 : sep + 1);
    // Gemmi✔️✔️:     size_t end = cid.find_first_of("[:;", pos);
    // Gemmi✔️✔️:     if (end != pos) {
    // Gemmi✔️✔️:       sel.atom_names = make_cid_list(cid, pos, end);
    // Gemmi✔️✔️:       // Chain name can be empty, but not atom name,
    // Gemmi✔️✔️:       // so we interpret empty atom name as *.
    // Gemmi✔️✔️:       if (!sel.atom_names.inverted && sel.atom_names.list.empty())
    // Gemmi✔️✔️:         sel.atom_names.all = true;
    // Gemmi✔️✔️:       if (end == std::string::npos)
    // Gemmi✔️✔️:         return;
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:     if (cid[end] == '[') {
    // Gemmi✔️✔️:       pos = end + 1;
    // Gemmi✔️✔️:       end = cid.find(']', pos);
    // Gemmi✔️✔️:       if (end == std::string::npos)
    // Gemmi✔️✔️:         wrong_syntax(cid, 0, " (no matching ']')");
    // Gemmi✔️✔️:       parse_cid_elements(cid, pos, sel.elements);
    // Gemmi✔️✔️:       ++end;
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:     if (cid[end] == ':') {
    // Gemmi✔️✔️:       pos = end + 1;
    // Gemmi✔️✔️:       sel.altlocs = make_cid_list(cid, pos, semi);
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:   }
    //
    // Behavior review: the guard compares against min(size, semi) — a
    // residue separator AT the string end or the semicolon skips the
    // whole stage (identity-on-skip defaults returned). Empty atom name
    // lists become wildcards ONLY when not inverted (source comment:
    // chain name can be empty, atom name cannot). end == npos after a
    // name list is the SOURCE EARLY RETURN: elements/altlocs stay at
    // defaults even if junk follows in a longer buffer. The element
    // bracket must find ']' (else " (no matching ']')" at pos 0); the
    // C07 owner reads to ']' itself. Altlocs run to the semicolon (npos
    // -> end of string). Markers ✔️✔️: one scan plus C04/C07 owners.
    let default_names = SelectionCidList {
        all: true,
        inverted: false,
        list: Vec::new(),
    };
    if residue_sep.unwrap_or(usize::MAX) < cid.len().min(semi.unwrap_or(usize::MAX)) {
        let bytes = cid.as_bytes();
        let at = |index: usize| bytes.get(index).copied().unwrap_or(0);
        let mut pos = match residue_sep {
            Some(0) | None => 0,
            Some(sep) => sep + 1,
        };
        let mut end = bytes[pos..]
            .iter()
            .position(|byte| b"[:;".contains(byte))
            .map(|index| index + pos);
        let mut atom_names = SelectionCidList {
            all: true,
            inverted: false,
            list: Vec::new(),
        };
        if end != Some(pos) {
            atom_names = make_cid_list(cid, pos, end.unwrap_or(cid.len()), "-[]()!/*.:;")
                .map_err(BioSelectionParseError::Syntax)?;
            // Empty atom name -> '*' only when not inverted.
            if !atom_names.inverted && atom_names.list.is_empty() {
                atom_names.all = true;
            }
            if end.is_none() {
                return Ok((
                    atom_names,
                    None,
                    SelectionCidList {
                        all: true,
                        inverted: false,
                        list: Vec::new(),
                    },
                    true,
                ));
            }
        }
        let mut elements = None;
        let mut altlocs = SelectionCidList {
            all: true,
            inverted: false,
            list: Vec::new(),
        };
        if at(end.unwrap_or(usize::MAX)) == b'[' {
            pos = end.unwrap() + 1;
            end = bytes[pos..]
                .iter()
                .position(|byte| *byte == b']')
                .map(|index| index + pos);
            if end.is_none() {
                return Err(BioSelectionParseError::Syntax(wrong_syntax(
                    cid,
                    0,
                    Some(" (no matching ']')"),
                )));
            }
            elements = parse_cid_elements(cid, pos).map_err(BioSelectionParseError::Syntax)?;
            end = Some(end.unwrap() + 1);
        }
        if at(end.unwrap_or(usize::MAX)) == b':' {
            pos = end.unwrap() + 1;
            altlocs = make_cid_list(cid, pos, semi.unwrap_or(cid.len()), "-[]()!/*.:;")
                .map_err(BioSelectionParseError::Syntax)?;
        }
        return Ok((atom_names, elements, altlocs, false));
    }
    Ok((
        default_names,
        None,
        SelectionCidList {
            all: true,
            inverted: false,
            list: Vec::new(),
        },
        false,
    ))
}

/// CID extension stage (BIO-CID C13): the semicolon loop after the atom
/// stage. Inequalities accumulate; entity-type lists are re-assigned each
/// iteration (last wins), with flags re-filled then re-set.
pub(crate) struct CidExtensionStage {
    pub inequalities: Vec<SelectionAtomInequality>,
    pub entity_types: Option<SelectionCidList>,
    pub polymer: bool,
    pub water: bool,
}

pub(crate) fn parse_cid_extensions(
    cid: &str,
    semi: Option<usize>,
) -> Result<CidExtensionStage, BioSelectionParseError> {
    // Gemmi✔️✔️:   // extensions after semicolon(s)
    // Gemmi✔️✔️:   while (semi < cid.size()) {
    // Gemmi✔️✔️:     size_t pos = semi + 1;
    // Gemmi✔️✔️:     while (cid[pos] == ' ')
    // Gemmi✔️✔️:       ++pos;
    // Gemmi✔️✔️:     semi = std::min(cid.find(';', pos), cid.size());
    // Gemmi✔️✔️:     size_t end = semi;
    // Gemmi✔️✔️:     while (end > pos && cid[end-1] == ' ')
    // Gemmi✔️✔️:       --end;
    // Gemmi✔️✔️:     if (has_inequality(cid, pos, end)) {
    // Gemmi✔️✔️:       sel.atom_inequalities.push_back(parse_atom_inequality(cid, pos, end));
    // Gemmi✔️✔️:     } else {
    // Gemmi✔️✔️:       sel.entity_types = make_cid_list(cid, pos, end);
    // Gemmi✔️✔️:       bool inv = sel.entity_types.inverted;
    // Gemmi✔️✔️:       std::fill(sel.et_flags.begin(), sel.et_flags.end(), char(inv));
    // Gemmi✔️✔️:       for (const std::string& item : split_str(sel.entity_types.list, ',')) {
    // Gemmi✔️✔️:         EntityType et = EntityType::Unknown;
    // Gemmi✔️✔️:         if (item == "polymer")
    // Gemmi✔️✔️:           et = EntityType::Polymer;
    // Gemmi✔️✔️:         else if (item == "solvent")
    // Gemmi✔️✔️:           et = EntityType::Water;
    // Gemmi✔️✔️:         else
    // Gemmi✔️✔️:           wrong_syntax(cid, 0, (" at " + item).c_str());
    // Gemmi✔️✔️:         sel.et_flags[(int)et] = char(!inv);
    // Gemmi✔️✔️:       }
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:   }
    //
    // Behavior review: EVERY loop iteration re-assigns entity_types (last
    // extension wins) and re-fills the flags to inv before setting matched
    // items to !inv, while atom_inequalities accumulate across iterations.
    // Leading spaces after ';' are skipped unbounded (operator[] reads 0
    // past the end); trailing spaces are trimmed back from the next ';'.
    // split_str always pushes the final remainder (util.hpp:146-161), so
    // an empty extension (";") yields the single item "" and errors with
    // the composed " at " note. No extension item explicitly matches
    // EntityType::Unknown, but the inversion fill still writes its flag
    // bit true for inverted lists.
    // Markers ✔️✔️: one scan per extension; items are ASCII-comma slices
    // of the validated list text, so no lossy conversion.
    let bytes = cid.as_bytes();
    let at = |index: usize| bytes.get(index).copied().unwrap_or(0);
    let mut inequalities = Vec::new();
    let mut entity_types = None;
    let mut polymer = false;
    let mut water = false;
    let mut semi = semi.unwrap_or(usize::MAX);
    while semi < cid.len() {
        let mut pos = semi + 1;
        while at(pos) == b' ' {
            pos += 1;
        }
        semi = bytes[pos..]
            .iter()
            .position(|byte| *byte == b';')
            .map(|index| index + pos)
            .unwrap_or(cid.len())
            .min(cid.len());
        let mut end = semi;
        while end > pos && at(end - 1) == b' ' {
            end -= 1;
        }
        if has_inequality(cid, pos, end) {
            let mut cursor = pos;
            inequalities.push(
                parse_atom_inequality(cid, &mut cursor, end)
                    .map_err(BioSelectionParseError::Syntax)?,
            );
        } else {
            let list = make_cid_list(cid, pos, end, "-[]()!/*.:;")
                .map_err(BioSelectionParseError::Syntax)?;
            let inv = list.inverted;
            polymer = inv;
            water = inv;
            // Caller invariant (BIO-C14-CLOSE Step2): `list.list` holds
            // only bytes copied by make_cid_list from the `cid: &str`
            // window, so it is valid UTF-8 by construction. Splitting at
            // ASCII ',' can never land inside a multibyte code point, so
            // each borrowed member keeps that validity. A failure here is
            // a production bug, not CID input state; there is no empty
            // substitution or lossy fallback path.
            let list_text = std::str::from_utf8(&list.list)
                .expect("CID list members are slices of the cid: &str input");
            for item in split_cid_list_items(list_text) {
                match item {
                    "polymer" => polymer = !inv,
                    "solvent" => water = !inv,
                    other => {
                        return Err(BioSelectionParseError::Syntax(wrong_syntax(
                            cid,
                            0,
                            Some(&format!(" at {other}")),
                        )));
                    }
                }
            }
            entity_types = Some(list);
        }
    }
    Ok(CidExtensionStage {
        inequalities,
        entity_types,
        polymer,
        water,
    })
}

/// Split comma list members exactly like gemmi split_str (always keeps
/// the final remainder, including an empty one). The verbatim source
/// anchors and the behavior/complexity review live inside the body
/// (BIO-C14-ANCHOR Step2; source-reproduction protocol §1).
fn split_cid_list_items(list: &str) -> impl Iterator<Item = &str> {
    // Gemmi✔️✔️: template<typename S>
    // Gemmi✔️✔️: void split_str_into(const std::string& str, S sep,
    // Gemmi✔️✔️:                     std::vector<std::string>& result) {
    // Gemmi✔️✔️:   std::size_t start = 0, end;
    // Gemmi✔️✔️:   while ((end = str.find(sep, start)) != std::string::npos) {
    // Gemmi✔️✔️:     result.emplace_back(str, start, end - start);
    // Gemmi✔️✔️:     start = end + impl::length(sep);
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   result.emplace_back(str, start);
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: template<typename S>
    // Gemmi✔️✔️: std::vector<std::string> split_str(const std::string& str, S sep) {
    // Gemmi✔️✔️:   std::vector<std::string> result;
    // Gemmi✔️✔️:   split_str_into(str, sep, result);
    // Gemmi✔️✔️:   return result;
    // Gemmi✔️✔️: }
    //
    // Behavior review: `str::split(',')` walks the same delimiter endpoints
    // as the pinned find loop — every delimiter position yields the member
    // before it (possibly empty) and the final remainder is always produced,
    // so `""` -> `[""]`, `"a,"` -> `["a", ""]`. Complexity review: one
    // forward scan; unlike the source vector<std::string>, each member here
    // is a BORROWED slice — zero per-member allocation, and the caller keeps
    // owning the one validated list text.
    list.split(',')
}

/// Complete CID decode (BIO-CID C14): the Selection ctor entry composed
/// from the C09-C13 stage owners. No parser duplication — every stage is
/// a call into its owner.
pub(crate) fn read_bio_selection_cid(
    cid: &str,
) -> Result<BioSelectionData, BioSelectionParseError> {
    // Gemmi✔️✔️:   Selection::Selection(const std::string& cid) {
    // Gemmi✔️✔️:     parse_cid(cid, *this);
    // Gemmi✔️✔️:   }
    //
    // parse_cid skeleton (select.cpp:131-257): the bypass preamble, one
    // `semi` scan, the C03 omission class, and the threaded `sep` that
    // every stage reassigns only at its pinned source assignment (C10-SEP
    // identity-on-skip). The atom stage's source early return (names to
    // end of string) skips the extension loop entirely. Every list window
    // is a substring of the `cid: &str`, so `from_parts` cannot observe
    // invalid UTF-8; its error type is unreachable here by construction.
    if cid_selection_bypassed(cid) {
        return Ok(BioSelectionData::default());
    }
    let semi = cid.as_bytes().iter().position(|byte| *byte == b';');
    let omit = match determine_omitted_cid_fields(cid) {
        OmittedCidFields::Model => 0,
        OmittedCidFields::Chain => 1,
        OmittedCidFields::Residue => 2,
        OmittedCidFields::Atom => 3,
    };
    let (mdl, model_sep) = parse_cid_model_stage(cid, omit, semi)?;
    let (chain_ids, chain_sep) = parse_cid_chain_stage(cid, omit, semi, model_sep)?;
    let (from_seqid, to_seqid, residue_names, residue_sep) =
        parse_cid_residue_stage(cid, omit, semi, chain_sep)?;
    let (atom_names, elements, altlocs, early_return) =
        parse_cid_atom_stage(cid, semi, residue_sep)?;
    let mut entity_types = SelectionCidList {
        all: true,
        inverted: false,
        list: Vec::new(),
    };
    let mut et_flags = [false; 6];
    let mut atom_inequalities = Vec::new();
    if !early_return {
        let ext = parse_cid_extensions(cid, semi)?;
        atom_inequalities = ext.inequalities;
        if let Some(list) = ext.entity_types {
            // std::fill(et_flags, inv) then matched items set !inv. No item
            // explicitly matches Unknown(0), but the fill still writes its
            // bit; matched items flip at Polymer(1)/Water(4).
            et_flags = [list.inverted; 6];
            et_flags[1] = ext.polymer;
            et_flags[4] = ext.water;
            entity_types = list;
        }
    }
    let parts = BioSelectionParts {
        mdl,
        chain_ids,
        from_seqid,
        to_seqid,
        residue_names,
        entity_types,
        et_flags,
        atom_names,
        elements,
        altlocs,
        atom_inequalities,
    };
    match BioSelectionData::from_parts(parts) {
        Ok(data) => Ok(data),
        // Unreachable by construction: every list window is a substring
        // of the `cid: &str` input, so no carrier can hold invalid UTF-8.
        // This is an invariant statement, not a fallback path.
        Err(_) => unreachable!("cid substrings are valid UTF-8"),
    }
}

/// CID List serialization (BIO-CID C15) against the canonical borrowed
/// BIO `SelectionList` exposed by `parts()`; anchors/reviews inside body.
pub(crate) fn cid_list_str(list: &SelectionList) -> String {
    // Gemmi✔️✔️:     std::string str() const {
    // Gemmi✔️✔️:       if (all)
    // Gemmi✔️✔️:         return "*";
    // Gemmi✔️✔️:       return inverted ? "!" + list : list;
    // Gemmi✔️✔️:     }
    // (select.hpp List::str, inline.)
    //
    // Behavior review: wildcard prints "*" regardless of members; an
    // inverted non-wildcard list prints "!" + the exact member text
    // (empty members and multibyte names preserved exactly through the
    // canonical String members); a plain non-wildcard list returns the
    // member text itself, even when empty. Complexity review: one output
    // allocation for the inverted case ("!" + members copied once — the
    // source "!" + list concatenation also allocates once; no evidence
    // is claimed about SSO internals); the plain case returns a copy of
    // the members, as the source returns its list by value. No Vec
    // staging and no byte/char transcoding: borrowed String members only.
    if list.all {
        return "*".to_string();
    }
    if list.inverted {
        let mut text = String::with_capacity(1 + list.list.len());
        text.push('!');
        text.push_str(&list.list);
        text
    } else {
        list.list.clone()
    }
}

/// CID SequenceId serialization (BIO-CID C15) returning exact source
/// bytes; anchors/reviews inside body.
pub(crate) fn seqid_str(seqid: &SelectionSequenceId) -> Vec<u8> {
    // Gemmi✔️✔️: std::string Selection::SequenceId::str() const {
    // Gemmi✔️✔️:   std::string s;
    // Gemmi✔️✔️:   if (!empty()) {
    // Gemmi✔️✔️:     s = std::to_string(seqnum);
    // Gemmi✔️✔️:     if (icode != '*') {
    // Gemmi✔️✔️:       s += '.';
    // Gemmi✔️✔️:       if (icode != ' ')
    // Gemmi✔️✔️:         s += icode;
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   return s;
    // Gemmi✔️✔️: }
    // with `empty()` (select.hpp:56-58):
    // Gemmi✔️✔️:     bool empty() const {
    // Gemmi✔️✔️:       return seqnum == INT_MIN || seqnum == INT_MAX;
    // Gemmi✔️✔️:     }
    //
    // Behavior review: BOTH sentinel numbers serialize to empty (source
    // emptiness, not a placeholder); a real number renders in full;
    // icode '*' adds no suffix, ' ' adds a dot only, ANY other u8 appends
    // its EXACT raw byte after the dot — the C++ `s += icode` appends a
    // raw char byte, and a Vec<u8> preserves it with no lossy/Unicode
    // substitution (the unapproved error-Display lossy exception was
    // removed here with BIO-C15-BYTES). Public String conversion is
    // deferred to the final decoded-CID boundary (C22) under explicit
    // ASCII-icode provenance. Complexity review: the i32 is rendered
    // once into the owned output (to_string -> into_bytes; no second
    // intermediate digit buffer — same single-render cost class as
    // std::to_string(int)); suffix appends are constant-time byte pushes
    // on that same owned buffer.
    if seqid.seqnum == i32::MIN || seqid.seqnum == i32::MAX {
        return Vec::new();
    }
    let mut bytes = seqid.seqnum.to_string().into_bytes();
    if seqid.icode != b'*' {
        bytes.push(b'.');
        if seqid.icode != b' ' {
            bytes.push(seqid.icode);
        }
    }
    bytes
}

/// CID AtomInequality serialization (BIO-CID C21) returning exact source
/// bytes through the canonical %.9g owner (C20); anchors/reviews inside
/// body.
pub(crate) fn inequality_str(inequality: &SelectionAtomInequality) -> Vec<u8> {
    // Gemmi✔️✔️: std::string Selection::AtomInequality::str() const {
    // Gemmi✔️✔️:   std::string r = ";";
    // Gemmi✔️✔️:   r += property;
    // Gemmi✔️✔️:   r += relation == 0 ? '=' : relation < 0 ? '<' : '>';
    // Gemmi✔️✔️:   r += to_str(value);
    // Gemmi✔️✔️:   return r;
    // Gemmi✔️✔️: }
    // (select.cpp:280-285, verbatim; `to_str(double)` is sprintf.hpp:36-40,
    // ported as the canonical %.9g owner format_cif_f64 — BIO-CID C20.)
    //
    // Behavior review: leading semicolon, then the property as its EXACT
    // raw byte (the C++ `r += property` appends a raw char byte; the
    // parser domain is b'q'/b'b', and any from_parts byte is preserved
    // with NO lossy/Unicode substitution, matching the C15-BYTES rule;
    // public String conversion happens only at the C22 decoded-CID
    // boundary), then '=' for relation 0, '<' for negative, '>' for
    // positive, then the value serialized by the canonical %.9g owner
    // (nine significant digits: rounding-induced non-roundtrips such as a
    // value whose decimal input had more digits serialize to the rounded
    // owner output, never to the original input text). Complexity review
    // (BIO-C21-FIX): ONE canonical-owner value String — the owner
    // format_cif_f64 is invoked exactly once and stored, the output Vec
    // reserves exactly 3 + value length once (semicolon, property byte,
    // relation char prefix), and the stored result is appended once; the
    // earlier version invoked the owner twice and underreserved a
    // two-byte prefix. These are Rust-only allocation facts — no source
    // std::string SSO/ABI allocation-parity claim.
    let value_text = format_cif_f64(inequality.value);
    let mut bytes = Vec::with_capacity(3 + value_text.len());
    bytes.push(b';');
    bytes.push(inequality.property);
    bytes.push(if inequality.relation == 0 {
        b'='
    } else if inequality.relation < 0 {
        b'<'
    } else {
        b'>'
    });
    bytes.extend_from_slice(value_text.as_bytes());
    bytes
}

/// Complete CID serialization (BIO-CID C22) over the borrowed predicate
/// view, composing the C15 list/seqid and C21 inequality serializers and
/// the canonical %.9g owner; anchors/reviews inside body.
pub(crate) fn bio_selection_to_cid(data: &BioSelectionData) -> String {
    // Gemmi✔️✔️: std::string Selection::str() const {
    // Gemmi✔️✔️:   std::string cid = "/";
    // Gemmi✔️✔️:   if (mdl != 0)
    // Gemmi✔️✔️:     cid += std::to_string(mdl);
    // Gemmi✔️✔️:   cid += '/';
    // Gemmi✔️✔️:   cid += chain_ids.str();
    // Gemmi✔️✔️:   cid += '/';
    // Gemmi✔️✔️:   cid += from_seqid.str();
    // Gemmi✔️✔️:   if (!residue_names.all) {
    // Gemmi✔️✔️:     cid += '(';
    // Gemmi✔️✔️:     cid += residue_names.str();
    // Gemmi✔️✔️:     cid += ')';
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   if ((!from_seqid.empty() || !to_seqid.empty()) &&
    // Gemmi✔️✔️:       (from_seqid.seqnum != to_seqid.seqnum || from_seqid.icode != to_seqid.icode)) {
    // Gemmi✔️✔️:     cid += '-';
    // Gemmi✔️✔️:     cid += to_seqid.str();
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   cid += '/';
    // Gemmi✔️✔️:   if (!atom_names.all)
    // Gemmi✔️✔️:     cid += atom_names.str();
    // Gemmi✔️✔️:   if (!elements.empty()) {
    // Gemmi✔️✔️:     cid += '[';
    // Gemmi✔️✔️:     bool inv = (std::count(elements.begin(), elements.end(), 1) > 64);
    // Gemmi✔️✔️:     if (inv)
    // Gemmi✔️✔️:       cid += '!';
    // Gemmi✔️✔️:     for (size_t i = 0; i < elements.size(); ++i)
    // Gemmi✔️✔️:       if (elements[i] != char(inv)) {
    // Gemmi✔️✔️:         cid += element_name(static_cast<El>(i));
    // Gemmi✔️✔️:         cid += ',';
    // Gemmi✔️✔️:       }
    // Gemmi✔️✔️:     cid.back() = ']';
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   if (!altlocs.all) {
    // Gemmi✔️✔️:     cid += ':';
    // Gemmi✔️✔️:     cid += altlocs.str();
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   if (!entity_types.all) {
    // Gemmi✔️✔️:     cid += ';';
    // Gemmi✔️✔️:     cid += entity_types.str();
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   for (const AtomInequality& ai : atom_inequalities)
    // Gemmi✔️✔️:     cid += ai.str();
    // Gemmi✔️✔️:   return cid;
    // Gemmi✔️✔️: }
    // (select.cpp:288-327, verbatim.)
    //
    // Behavior review: exact source field order (model, chain, residue
    // names in parentheses, to-seqid range only when the endpoints differ
    // and at least one is non-empty, atom names bare, element mask with the
    // STRICT count > 64 inversion threshold, altloc after ':', entity
    // after ';', inequalities in order). The cid.back() = ']' replacement
    // is reproduced for DEGENERATE masks too: when every mask bit equals
    // char(inv) nothing is appended and the '[' itself becomes ']',
    // exactly like the source. Byte assembly uses the C15 borrowed-list
    // serializer, the C15 raw-byte seqid serializer, and the C21
    // inequality serializer; the single String conversion happens HERE at
    // the final decoded-CID boundary under ASCII provenance — CID-decoded
    // selections carry only ASCII icodes/properties, so the conversion
    // cannot substitute bytes (fail-closed expect, C14-CLOSE precedent).
    // Complexity review: one growing owned byte vector with amortized
    // appends (source std::string growth class); list serializers borrow
    // the canonical view members (C15), the element pass is one linear
    // mask scan with a table lookup per selected element (source loop
    // shape), and each inequality costs its canonical-owner value String
    // (C20 Rust-only allocation evidence; no source SSO/ABI parity claim).
    let parts = data.parts();
    let seqid_empty =
        |seqid: &SelectionSequenceId| seqid.seqnum == i32::MIN || seqid.seqnum == i32::MAX;
    let mut cid: Vec<u8> = Vec::with_capacity(32);
    cid.push(b'/');
    if parts.mdl != 0 {
        cid.extend_from_slice(parts.mdl.to_string().as_bytes());
    }
    cid.push(b'/');
    cid.extend_from_slice(cid_list_str(parts.chain_ids).as_bytes());
    cid.push(b'/');
    cid.extend_from_slice(&seqid_str(parts.from_seqid));
    if !parts.residue_names.all {
        cid.push(b'(');
        cid.extend_from_slice(cid_list_str(parts.residue_names).as_bytes());
        cid.push(b')');
    }
    if (!seqid_empty(parts.from_seqid) || !seqid_empty(parts.to_seqid))
        && (parts.from_seqid.seqnum != parts.to_seqid.seqnum
            || parts.from_seqid.icode != parts.to_seqid.icode)
    {
        cid.push(b'-');
        cid.extend_from_slice(&seqid_str(parts.to_seqid));
    }
    cid.push(b'/');
    if !parts.atom_names.all {
        cid.extend_from_slice(cid_list_str(parts.atom_names).as_bytes());
    }
    if let Some(mask) = parts.elements {
        cid.push(b'[');
        let inverted = mask.iter().filter(|&&bit| bit).count() > 64;
        if inverted {
            cid.push(b'!');
        }
        for (index, &bit) in mask.iter().enumerate() {
            if bit != inverted {
                // Bounded ordinal invariant: mask indices are within the
                // [bool; 120] mask, so the display owner always resolves;
                // fail-closed evidence, never an invented default name.
                let name =
                    gemmi_element_name(index).expect("element mask ordinal is within GEMMI_EL_END");
                cid.extend_from_slice(name.as_bytes());
                cid.push(b',');
            }
        }
        // cid.back() = ']' — including degenerate masks, where nothing was
        // appended and the '[' itself is replaced by ']'.
        *cid.last_mut().expect("cid has at least the pushed '['") = b']';
    }
    if !parts.altlocs.all {
        cid.push(b':');
        cid.extend_from_slice(cid_list_str(parts.altlocs).as_bytes());
    }
    if !parts.entity_types.all {
        cid.push(b';');
        cid.extend_from_slice(cid_list_str(parts.entity_types).as_bytes());
    }
    for inequality in parts.atom_inequalities {
        cid.extend_from_slice(&inequality_str(inequality));
    }
    String::from_utf8(cid).expect("decoded-CID provenance: sequence/property bytes are ASCII")
}

/// Structured CID decode failure (BIO-CID C02 frozen boundary, corrected
/// by BIO-C01-FIX Step 12; public parse-error vocabulary since BIO-CID
/// C26 with private-field carriers). `Syntax` carries the exact-byte
/// carrier; `SeqidRange` is the approved typed integer-range variant.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BioSelectionParseError {
    /// Exact `wrong_syntax` evidence: full cid, byte offset, info note.
    Syntax(SelectionSyntaxError),
    /// Approved typed i32 bounds failure (offending text preserved).
    SeqidRange(SelectionSeqidRangeError),
}

impl BioSelectionParseError {
    /// The real recorded source byte offset, if one exists.
    /// `Syntax` → `Some(pos)` (the actual `wrong_syntax` position);
    /// `SeqidRange` → `None` — that carrier records only the offending
    /// text and no offset; a placeholder 0 must never be fabricated.
    pub(crate) fn offset(&self) -> Option<usize> {
        match self {
            Self::Syntax(error) => Some(error.pos),
            Self::SeqidRange(_) => None,
        }
    }

    /// Exact message bytes of the underlying carrier. For `Syntax` this
    /// is the byte-exact `wrong_syntax` message (including split-window
    /// cut bytes); for `SeqidRange` the carrier's own composed message.
    /// `Display` remains a lossy presentation only.
    pub(crate) fn source_message_bytes(&self) -> Vec<u8> {
        match self {
            Self::Syntax(error) => error.source_message_bytes(),
            Self::SeqidRange(error) => error.to_string().into_bytes(),
        }
    }
}

impl std::fmt::Display for BioSelectionParseError {
    /// Lossy presentation only (BIO-SEL-FIX B26 precedent); exact bytes
    /// come from `source_message_bytes`.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Syntax(error) => write!(f, "{error}"),
            Self::SeqidRange(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for BioSelectionParseError {}

/// Public CID decode entry (BIO-CID C27): the registered root
/// `BioSelection::from_cid` delegates here; the body IS the canonical
/// reader (select.cpp:261-263 `Selection::Selection(cid){parse_cid}`),
/// so no second parser exists.
pub fn read_bio_selection(cid: &str) -> Result<BioSelectionData, BioSelectionParseError> {
    read_bio_selection_cid(cid)
}

/// Public CID serialize entry (BIO-CID C28): the registered root
/// `BioSelection::to_cid` delegates here; the body IS the canonical
/// C22 serializer (select.cpp:288-327), so no second formatter exists.
pub fn write_bio_selection(data: &BioSelectionData) -> String {
    bio_selection_to_cid(data)
}

#[cfg(test)]
mod tests {
    use super::{
        BioSelectionData, BioSelectionParseError, BioSelectionParts, CidIntScan, OmittedCidFields,
        SelectionAtomInequality, SelectionCidList, SelectionList, SelectionSeqidRangeError,
        SelectionSequenceId, bio_selection_to_cid, cid_list_str, cid_selection_bypassed,
        determine_omitted_cid_fields, inequality_str, make_cid_list, parse_atom_inequality,
        parse_cid_atom_stage, parse_cid_chain_stage, parse_cid_extensions, parse_cid_model_stage,
        parse_cid_residue_stage, parse_cid_seqid, read_bio_selection_cid, scan_cid_int_prefix,
        seqid_str, wrong_syntax,
    };

    // ---- BIO-CID C13: extension loop. Expectations from select.cpp:229-257
    // + split_str (util.hpp:146-161, remainder always kept). ----

    #[test]
    fn bio_cid_c13_orders_repetitions_and_inequalities() {
        // Single entity extension sets its flag and clears the other.
        let ext = parse_cid_extensions(";polymer", Some(0)).unwrap();
        assert_eq!((ext.polymer, ext.water), (true, false));
        assert_eq!(ext.entity_types.as_ref().unwrap().list, b"polymer");
        assert!(ext.inequalities.is_empty());
        let ext = parse_cid_extensions(";solvent", Some(0)).unwrap();
        assert_eq!((ext.polymer, ext.water), (false, true));
        // Repetition: last entity extension wins and re-fills the flags.
        let ext = parse_cid_extensions(";polymer;solvent", Some(0)).unwrap();
        assert_eq!((ext.polymer, ext.water), (false, true));
        // Combined list members set both flags.
        let ext = parse_cid_extensions(";polymer,solvent", Some(0)).unwrap();
        assert_eq!((ext.polymer, ext.water), (true, true));
        // Inequalities accumulate across extensions in order.
        let ext = parse_cid_extensions(";q>1;b<2;q=3", Some(0)).unwrap();
        assert_eq!(ext.inequalities.len(), 3);
        assert_eq!(ext.inequalities[0].relation, 1);
        assert_eq!(ext.inequalities[1].relation, -1);
        assert_eq!(ext.inequalities[2].relation, 0);
        // A later entity extension still replaces an earlier entity list
        // while earlier inequalities stay.
        let ext = parse_cid_extensions(";q>1;polymer", Some(0)).unwrap();
        assert_eq!(ext.inequalities.len(), 1);
        assert_eq!(ext.entity_types.as_ref().unwrap().list, b"polymer");
        assert!(ext.polymer && !ext.water);
    }

    #[test]
    fn bio_cid_c13_inverted_entities_and_whitespace() {
        // Inverted list: flags are filled with inv then matched items set
        // !inv — "!polymer" means water (and unknown) stay inverted-on.
        let ext = parse_cid_extensions(";!polymer", Some(0)).unwrap();
        assert_eq!((ext.polymer, ext.water), (false, true));
        assert!(ext.entity_types.as_ref().unwrap().inverted);
        // Whitespace: leading spaces after ';' skipped; trailing spaces
        // before the next ';' trimmed from the list window. A space AFTER
        // a numeric value is consumed by the unbounded skip and then
        // fails `pos != end` (source-exact), so the inequality window
        // itself must be space-free at its end.
        let ext = parse_cid_extensions(";  polymer  ; q>1", Some(0)).unwrap();
        assert_eq!(ext.entity_types.as_ref().unwrap().list, b"polymer");
        assert!(ext.polymer && !ext.water);
        assert_eq!(ext.inequalities.len(), 1);
        let error = match parse_cid_extensions("; q>1 ", Some(0)) {
            Err(error) => error,
            Ok(stage) => panic!("expected error, got {:?}", (stage.polymer, stage.water)),
        };
        match error {
            BioSelectionParseError::Syntax(error) => {
                assert_eq!(error.pos, 6);
                assert_eq!(error.info.as_deref(), None);
            }
            other => panic!("expected Syntax, got {other:?}"),
        }
        // Spaces INSIDE a list member do not match the fixed spellings.}        // Spaces INSIDE a list member do not match the fixed spellings.
        let error = match parse_cid_extensions("; polymer , solvent", Some(0)) {
            Err(error) => error,
            Ok(stage) => panic!("expected error, got {:?}", (stage.polymer, stage.water)),
        };
        match error {
            BioSelectionParseError::Syntax(error) => {
                assert_eq!(error.pos, 0);
                assert_eq!(error.info.as_deref(), Some(" at polymer "));
            }
            other => panic!("expected Syntax, got {other:?}"),
        }
    }

    #[test]
    fn bio_cid_c13_empty_and_trailing_extensions() {
        // Empty extension: split_str keeps the empty remainder -> the
        // composed " at " note with an empty item.
        let error = match parse_cid_extensions(";", Some(0)) {
            Err(error) => error,
            Ok(stage) => panic!("expected error, got {:?}", (stage.polymer, stage.water)),
        };
        match error {
            BioSelectionParseError::Syntax(error) => {
                assert_eq!(error.pos, 0);
                assert_eq!(error.info.as_deref(), Some(" at "));
            }
            other => panic!("expected Syntax, got {other:?}"),
        }
        // Unknown entity name is rejected with the item echoed.
        let error = match parse_cid_extensions(";ligand", Some(0)) {
            Err(error) => error,
            Ok(stage) => panic!("expected error, got {:?}", (stage.polymer, stage.water)),
        };
        match error {
            BioSelectionParseError::Syntax(error) => {
                assert_eq!(error.info.as_deref(), Some(" at ligand"));
            }
            other => panic!("expected Syntax, got {other:?}"),
        }
        // No extension region (semi npos / at end): empty stage.
        let ext = parse_cid_extensions("A", None).unwrap();
        assert!(ext.inequalities.is_empty() && ext.entity_types.is_none());
        let ext = parse_cid_extensions("A;polymer", Some(1)).unwrap();
        assert!(ext.polymer && ext.inequalities.is_empty());
        // A lone trailing semicolon after an entity still re-parses the
        // empty window and errors.
        let error = match parse_cid_extensions(";polymer;", Some(0)) {
            Err(error) => error,
            Ok(stage) => panic!("expected error, got {:?}", (stage.polymer, stage.water)),
        };
        match error {
            BioSelectionParseError::Syntax(error) => {
                assert_eq!(error.info.as_deref(), Some(" at "));
            }
            other => panic!("expected Syntax, got {other:?}"),
        }
    }

    #[test]
    fn bio_sel_b05_offset_positions_and_window_lengths() {
        use super::wrong_syntax;
        // pos == 0: no nearby text at all (source `if (pos != 0)`).
        let start = wrong_syntax("/a/b/c", 0, None);
        assert_eq!(start.to_string(), "Invalid selection syntax: /a/b/c");
        assert_eq!(start.pos, 0);
        assert_eq!(start.cid, "/a/b/c");
        // Middle offset: exactly eight bytes of nearby text.
        let middle = wrong_syntax("/aaa/bbb/ccc", 5, None);
        assert_eq!(
            middle.to_string(),
            format!(
                "Invalid selection syntax near \"{}\": /aaa/bbb/ccc",
                "bbb/ccc"
            )
        );
        // End offset: shorter-than-eight remainder is taken whole.
        let end = wrong_syntax("/atom", 3, None);
        assert_eq!(
            end.to_string(),
            format!("Invalid selection syntax near \"{}\": /atom", "om")
        );
        // Longer remainder: still exactly eight bytes.
        let long = wrong_syntax("0123456789ABCDEFGH", 2, None);
        assert_eq!(
            long.to_string(),
            format!(
                "Invalid selection syntax near \"{}\": 0123456789ABCDEFGH",
                "23456789"
            )
        );
        // Optional info is appended right after the base text.
        let info = wrong_syntax("/a[b]", 3, Some(" in [...]"));
        assert_eq!(
            info.to_string(),
            format!("Invalid selection syntax in [...] near \"{}\": /a[b]", "b]")
        );
        assert_eq!(info.info.as_deref(), Some(" in [...]"));
        // pos beyond the end: empty window, no panic.
        let beyond = wrong_syntax("cid", 10, None);
        assert_eq!(
            beyond.to_string(),
            format!("Invalid selection syntax near \"\": cid")
        );
    }

    #[test]
    fn bio_sel_b05_non_ascii_byte_boundaries_do_not_panic() {
        use super::wrong_syntax;
        // A two-byte character split by the offset renders lossily; the
        // retained cid and byte offset stay exact.
        let cid = "/\u{e9}abc"; // bytes: '/' 0xC3 0xA9 'a' 'b' 'c'
        let split = wrong_syntax(cid, 2, None);
        assert_eq!(split.pos, 2);
        assert_eq!(split.cid, cid);
        let text = split.to_string();
        assert!(text.starts_with("Invalid selection syntax near \""));
        // The window begins with the lone trailing byte of the character.
        assert!(text.contains("\u{fffd}abc"));
        assert!(text.ends_with(&format!(": {cid}")));
        // A whole multibyte character inside the window renders as itself.
        let whole = wrong_syntax(cid, 1, None);
        assert!(whole.to_string().contains("\u{e9}abc\""));
        // An eight-byte window over multibyte content still cuts at eight
        // bytes (four two-byte characters).
        let wide = "\u{e9}\u{e8}\u{e7}\u{e6}\u{e5}\u{e4}";
        let cut = wrong_syntax(wide, 0, None);
        assert_eq!(cut.pos, 0);
        assert_eq!(cut.to_string(), format!("Invalid selection syntax: {wide}"));
    }
    #[test]
    fn bio_sel_b06_leading_character_branches() {
        use super::OmittedCidFields::*;
        use super::determine_omitted_cid_fields as omit;
        // '/' -> model (select.cpp:29-30).
        assert_eq!(omit("/1/A"), Model);
        assert_eq!(omit("/"), Model);
        // digit, '.', '(' or '-' -> residue (select.cpp:31-32).
        assert_eq!(omit("14/A"), Residue);
        assert_eq!(omit("0"), Residue);
        assert_eq!(omit("9.a"), Residue);
        assert_eq!(omit("."), Residue);
        assert_eq!(omit("(ALA)"), Residue);
        assert_eq!(omit("-5"), Residue);
        // Non-ASCII first byte: none of the leading tests match (byte
        // values outside the compared ranges), so the separator scan runs;
        // the '(' separator classifies as residue.
        assert_eq!(omit("\u{e9}A("), Residue);
        assert_eq!(omit("\u{e9}A["), Atom);
    }

    #[test]
    fn bio_sel_b06_separator_branches() {
        use super::OmittedCidFields::*;
        use super::determine_omitted_cid_fields as omit;
        // No separator at all -> chain.
        assert_eq!(omit("A"), Chain);
        assert_eq!(omit("abc"), Chain);
        // First separator '/' or ';' -> chain.
        assert_eq!(omit("A/14"), Chain);
        assert_eq!(omit("AB;"), Chain);
        // '(' separator -> residue.
        assert_eq!(omit("A(ALA)"), Residue);
        // '[' or ':' separator -> atom.
        assert_eq!(omit("CA[C]"), Atom);
        assert_eq!(omit("CA:A"), Atom);
        // The FIRST separator decides: '(' later than ':' is irrelevant.
        assert_eq!(omit("CA:A(B)"), Atom);
        assert_eq!(omit("A(B):C"), Residue);
    }

    #[test]
    fn bio_sel_b06_empty_input_reads_as_chain() {
        // std::string cid[0] on empty yields '\0' (no leading branch), and
        // find_first_of finds nothing: classification is a bare chain.
        use super::OmittedCidFields::*;
        use super::determine_omitted_cid_fields as omit;
        assert_eq!(omit(""), Chain);
    }
    #[test]
    fn bio_sel_b07_wildcard_and_inversion_consume_one_byte() {
        use super::make_cid_list;
        // '*' sets all and consumes itself: empty member bytes.
        let wildcard = make_cid_list("/x/*/y", 3, 4, "-[]()!/*.:;").unwrap();
        assert!(wildcard.all);
        assert!(!wildcard.inverted);
        assert_eq!(wildcard.list, b"");
        // '!' sets inverted and consumes itself.
        let inverted = make_cid_list("/x/!A/y", 3, 5, "-[]()!/*.:;").unwrap();
        assert!(!inverted.all);
        assert!(inverted.inverted);
        assert_eq!(inverted.list, b"A");
        // A bare list keeps both flags false.
        let plain = make_cid_list("/x/AB/y", 3, 5, "-[]()!/*.:;").unwrap();
        assert!(!plain.all && !plain.inverted);
        assert_eq!(plain.list, b"AB");
        // '!' followed by '*': inversion consumes the '!' and the literal
        // '*' member is disallowed punctuation — the source errors here
        // with the composed note (verified against select.cpp:48-55).
        let bang_star = make_cid_list("/x/!*", 3, 5, "-[]()!/*.:;").unwrap_err();
        assert_eq!(bang_star.pos, 4);
        assert_eq!(bang_star.info.as_deref(), Some(" ('*' in a list)"));
    }

    #[test]
    fn bio_sel_b07_comma_members_and_empty_fields() {
        use super::make_cid_list;
        // Commas are never disallowed: empty members survive byte-exactly.
        let commas = make_cid_list("A,,B", 0, 4, "-[]()!/*.:;").unwrap();
        assert_eq!(commas.list, b"A,,B");
        // An entirely empty field (pos == end) with a non-special byte at
        // pos yields an empty plain list.
        let empty = make_cid_list("AB", 0, 0, "-[]()!/*.:;").unwrap();
        assert_eq!(empty.list, b"");
        assert!(!empty.all && !empty.inverted);
    }

    #[test]
    fn bio_sel_b07_disallowed_separators_report_position_and_info() {
        use super::make_cid_list;
        // Each disallowed byte fails with the advanced position and the
        // composed note; a leading '!'/'*' shifts the reported offset.
        let cases: [(&str, usize, usize, &str, usize, char); 5] = [
            ("A-B", 0, 3, "-[]()!/*.:;", 1, '-'),
            ("A[B", 0, 3, "-[]()!/*.:;", 1, '['),
            ("A.B", 0, 3, "-[]()!/*.:;", 1, '.'),
            ("!A/B", 0, 3, "-[]()!/*.:;", 2, '/'),
            ("*A:B", 0, 3, "-[]()!/*.:;", 2, ':'),
        ];
        for (cid, pos, end, disallowed, expected_pos, offender) in cases {
            let error = make_cid_list(cid, pos, end, disallowed).unwrap_err();
            assert_eq!(error.pos, expected_pos, "{cid}");
            assert_eq!(error.cid, cid);
            assert_eq!(
                error.info.as_deref(),
                Some(format!(" ('{offender}' in a list)").as_str())
            );
        }
    }

    #[test]
    fn bio_sel_b07_trailing_input_beyond_end_is_ignored() {
        use super::make_cid_list;
        // Only bytes before `end` belong to the field.
        let field = make_cid_list("AB/x", 0, 2, "-[]()!/*.:;").unwrap();
        assert_eq!(field.list, b"AB");
        // end past the string clamps like substr to the whole remainder.
        let clamped = make_cid_list("AB", 0, 99, "-[]()!/*.:;").unwrap();
        assert_eq!(clamped.list, b"AB");
        // pos at the very end reads the '\0' sentinel: plain empty list.
        let at_end = make_cid_list("AB", 2, 99, "-[]()!/*.:;").unwrap();
        assert_eq!(at_end.list, b"");
        assert!(!at_end.all && !at_end.inverted);
    }

    // BIO-SEL-FIX Step 6: byte-payload regressions. Caller-reachable cases
    // are marked; byte-boundary cuts and wrapped counts are impossible
    // from the pinned parse_cid callers (receipt audit) and are exercised
    // at helper level for exact substr fidelity.
    #[test]
    fn bio_sel_b07_multibyte_names_preserved_byte_exactly() {
        use super::make_cid_list;
        // CALLER-REACHABLE: multibyte UTF-8 member names between ASCII
        // separators keep their exact bytes ("Ä" = 0xC3 0x84). The field
        // ends at the BYTE offset of the following '/' ("/x/" is bytes
        // 0-2, "Ä" occupies 3-4, ',' 5, 'B' 6, '/' 7).
        let field = make_cid_list("/x/Ä,B/y", 3, 7, "-[]()!/*.:;").unwrap();
        assert_eq!(field.list, b"\xC3\x84,B");
        assert!(!field.all && !field.inverted);
        // Inversion in front of multibyte members consumes one byte.
        let inverted = make_cid_list("/x/!Ä/y", 3, 6, "-[]()!/*.:;").unwrap();
        assert!(inverted.inverted);
        assert_eq!(inverted.list, b"\xC3\x84");
        // A multibyte member after a wildcard consumption at the boundary.
        let wildcard = make_cid_list("/x/*Ä/y", 3, 6, "-[]()!/*.:;").unwrap();
        assert!(wildcard.all);
        assert_eq!(wildcard.list, b"\xC3\x84");
    }

    #[test]
    fn bio_sel_b07_byte_boundary_cuts_yield_exact_cut_bytes() {
        use super::make_cid_list;
        // IMPOSSIBLE FROM CALLERS (audit): end splitting the two-byte "Ä"
        // yields exactly the first byte — no lossy conversion, no empty
        // substitution, no panic.
        let cut = make_cid_list("ÄB", 0, 1, "-[]()!/*.:;").unwrap();
        assert_eq!(cut.list, b"\xC3");
        assert!(!cut.all && !cut.inverted);
        // Cut after the second byte of "Ä" inside a longer field.
        let cut2 = make_cid_list("ÄÄ", 0, 3, "-[]()!/*.:;").unwrap();
        assert_eq!(cut2.list, b"\xC3\x84\xC3");
        // IMPOSSIBLE FROM CALLERS: finite end < start reproduces the
        // source's wrapped `end - pos` count, i.e. substr(start, huge)
        // clamps to the remainder [start, len].
        let wrapped = make_cid_list("ABCD", 3, 1, "-[]()!/*.:;").unwrap();
        assert_eq!(wrapped.list, b"D");
        // The same wrap after '!' consumption (end 2 < start 3).
        let bang = make_cid_list("AB!CD", 2, 2, "-[]()!/*.:;").unwrap();
        assert!(bang.inverted);
        assert_eq!(bang.list, b"CD");
        // end == start after consumption is an ordinary empty field:
        // substr(3, 3-3=0) = "" (no wrap, verified against the source
        // subtraction which uses the incremented pos).
        let bang_empty = make_cid_list("AB!CD", 2, 3, "-[]()!/*.:;").unwrap();
        assert!(bang_empty.inverted);
        assert_eq!(bang_empty.list, b"");
        // start > len would mirror substr's out_of_range (panics; see the
        // implementation note) — deliberately not constructed here.
    }
    #[test]
    fn bio_sel_b08_numbers_wildcard_and_defaults() {
        use super::parse_cid_seqid;
        // Plain numbers advance the position; defaults are INT_MIN/INT_MAX
        // shaped by the caller.
        let mut pos = 0usize;
        let positive = parse_cid_seqid("14/A", &mut pos, i32::MIN).unwrap();
        assert_eq!((positive.seqnum, positive.icode, pos), (14, b' ', 2));
        let mut pos = 0usize;
        let negative = parse_cid_seqid("-7x", &mut pos, i32::MAX).unwrap();
        assert_eq!((negative.seqnum, negative.icode, pos), (-7, b'x', 3));
        let mut pos = 0usize;
        let zero = parse_cid_seqid("0", &mut pos, 5).unwrap();
        assert_eq!((zero.seqnum, zero.icode, pos), (0, b' ', 1));
        // Wildcard: seqnum keeps the default, icode becomes '*', one byte.
        let mut pos = 0usize;
        let wildcard = parse_cid_seqid("*-", &mut pos, i32::MIN).unwrap();
        assert_eq!((wildcard.seqnum, wildcard.icode, pos), (i32::MIN, b'*', 1));
        // No recognized leading byte: nothing consumed, default retained.
        let mut pos = 0usize;
        let untouched = parse_cid_seqid("A", &mut pos, i32::MAX).unwrap();
        assert_eq!(
            (untouched.seqnum, untouched.icode, pos),
            (i32::MAX, b' ', 0)
        );
    }

    #[test]
    fn bio_sel_b08_dot_and_insertion_code_grammar() {
        use super::parse_cid_seqid;
        // Optional dot before the icode; both "14.a" and "14a" are accepted.
        let mut pos = 0usize;
        let dotted = parse_cid_seqid("14.a/", &mut pos, i32::MIN).unwrap();
        assert_eq!((dotted.seqnum, dotted.icode, pos), (14, b'a', 4));
        let mut pos = 0usize;
        let bare = parse_cid_seqid("14a/", &mut pos, i32::MIN).unwrap();
        assert_eq!((bare.seqnum, bare.icode, pos), (14, b'a', 3));
        // '*' icode after a number.
        let mut pos = 0usize;
        let star = parse_cid_seqid("14*", &mut pos, i32::MIN).unwrap();
        assert_eq!((star.seqnum, star.icode, pos), (14, b'*', 3));
        // A trailing dot alone is consumed; no icode follows.
        let mut pos = 0usize;
        let dot = parse_cid_seqid("14.", &mut pos, i32::MIN).unwrap();
        assert_eq!((dot.seqnum, dot.icode, pos), (14, b' ', 3));
        // Wildcard then dot then letter: the letter OVERWRITES the wildcard
        // icode (the source assigns icode unconditionally in that branch).
        let mut pos = 0usize;
        let over = parse_cid_seqid("*.a", &mut pos, i32::MIN).unwrap();
        assert_eq!((over.seqnum, over.icode, pos), (i32::MIN, b'a', 3));
        // Nothing consumed before a letter: icode NOT taken.
        let mut pos = 0usize;
        let skip = parse_cid_seqid("a", &mut pos, i32::MIN).unwrap();
        assert_eq!((skip.seqnum, skip.icode, pos), (i32::MIN, b' ', 0));
        // Digit then non-alpha: icode not taken.
        let mut pos = 0usize;
        let num_only = parse_cid_seqid("3:", &mut pos, i32::MIN).unwrap();
        assert_eq!((num_only.seqnum, num_only.icode, pos), (3, b' ', 1));
    }

    #[test]
    fn bio_sel_b08_malformed_and_overflow() {
        use super::parse_cid_seqid;
        // A lone '-' performs strtol's no-conversion: seqnum 0, position
        // unchanged (so the caller sees the '-' still ahead).
        let mut pos = 0usize;
        let dash = parse_cid_seqid("-x", &mut pos, i32::MIN).unwrap();
        assert_eq!((dash.seqnum, dash.icode, pos), (0, b' ', 0));
        // Typed bounds error instead of C long->int truncation.
        let mut pos = 0usize;
        let error = parse_cid_seqid("99999999999x", &mut pos, i32::MIN).unwrap_err();
        assert_eq!(error.text, "99999999999");
        assert!(error.to_string().contains("99999999999"));
        let mut pos = 0usize;
        let negative = parse_cid_seqid("-3000000000", &mut pos, i32::MIN).unwrap_err();
        assert_eq!(negative.text, "-3000000000");
        // i32 boundary values themselves are representable.
        let mut pos = 0usize;
        let max = parse_cid_seqid("2147483647", &mut pos, 0).unwrap();
        assert_eq!(max.seqnum, i32::MAX);
        let mut pos = 0usize;
        let min = parse_cid_seqid("-2147483648", &mut pos, 0).unwrap();
        assert_eq!(min.seqnum, i32::MIN);
    }
    #[test]
    fn bio_sel_b09_all_recognized_symbols_and_case_rules() {
        use super::GEMMI_ELEMENT_NAMES;
        use super::parse_cid_elements;
        // Every entry of the pinned name table resolves to exactly its own
        // ordinal, including X and D.
        for (ordinal, name) in GEMMI_ELEMENT_NAMES.iter().enumerate() {
            let cid = format!("[{name}]");
            let mask = parse_cid_elements(&cid, 1).unwrap().unwrap();
            assert_eq!(mask[ordinal], true, "{name}");
            let selected: Vec<usize> = mask
                .iter()
                .enumerate()
                .filter(|(_, value)| **value)
                .map(|(index, _)| index)
                .collect();
            assert_eq!(selected, vec![ordinal], "{name}");
        }
        // Lowercase folds through `& ~0x20` exactly like the source —
        // including the "cn"/"CN" -> Copernicium(112) collision.
        for (lower, ordinal) in [
            ("c", 6usize),
            ("n", 7usize),
            ("zn", 30usize),
            ("d", 119usize),
            ("CN", 112usize),
            ("cn", 112usize),
        ] {
            let cid = format!("[{lower}]");
            let mask = parse_cid_elements(&cid, 1).unwrap().unwrap();
            assert!(mask[ordinal], "{lower}");
        }
        // Multiple members accumulate.
        let mask = parse_cid_elements("[C,N,O]", 1).unwrap().unwrap();
        assert!(mask[6] && mask[7] && mask[8]);
        assert!(!mask[1]);
    }

    #[test]
    fn bio_sel_b09_metals_nonmetals_and_inversion() {
        use super::GEMMI_IS_METAL;
        use super::parse_cid_elements;
        let metals = parse_cid_elements("[metals]", 1).unwrap().unwrap();
        for (index, value) in metals.iter().enumerate() {
            assert_eq!(*value, GEMMI_IS_METAL[index], "metals at {index}");
        }
        let nonmetals = parse_cid_elements("[nonmetals]", 1).unwrap().unwrap();
        for (index, value) in nonmetals.iter().enumerate() {
            assert_eq!(*value, !GEMMI_IS_METAL[index], "nonmetals at {index}");
        }
        // Spot pins read directly from the pinned table rows: Ge(32) and
        // Sb(51) are metals ("arbitrary division"), and so is Po(84);
        // X(0), H(1), He(2), At(85), Rn(86), Ts(117), Og(118) and D(119)
        // are not.
        for metal in [32usize, 51, 26, 92, 84] {
            assert!(GEMMI_IS_METAL[metal]);
        }
        for nonmetal in [0usize, 1, 2, 85, 86, 117, 118, 119] {
            assert!(!GEMMI_IS_METAL[nonmetal]);
        }
        // Inversion pre-fills inverted and flips the selected group.
        let inverted = parse_cid_elements("[!metals]", 1).unwrap().unwrap();
        for (index, value) in inverted.iter().enumerate() {
            assert_eq!(*value, !GEMMI_IS_METAL[index], "!metals at {index}");
        }
        let inverted_single = parse_cid_elements("[!C]", 1).unwrap().unwrap();
        assert!(!inverted_single[6]);
        assert!(inverted_single[1] && inverted_single[8]);
    }

    #[test]
    fn bio_sel_b09_unknown_and_malformed_tokens() {
        use super::parse_cid_elements;
        // Unknown single letters and two-letter spellings error. NOTE:
        // "CN" is NOT unknown — find_element uppercases before the name
        // scan, so "CN"/"cn" resolves to Copernicium (ordinal 112); that
        // collision is asserted in the case-rules test.
        for token in ["Q", "A", "Zz", "Xx", "METALS", "metal"] {
            let cid = format!("[{token}]");
            let error = parse_cid_elements(&cid, 1).unwrap_err();
            assert_eq!(error.pos, 0, "{token}");
            assert_eq!(error.cid, cid);
            assert!(
                error.info.as_deref() == Some(" in [...]")
                    || error.info.as_deref() == Some(" (invalid element in [...])"),
                "{token}: {:?}",
                error.info
            );
        }
        // Exact info distinction: letters give the invalid-element note,
        // over-long non-keywords give the in-list note.
        assert_eq!(
            parse_cid_elements("[Q]", 1).unwrap_err().info.as_deref(),
            Some(" (invalid element in [...])")
        );
        assert_eq!(
            parse_cid_elements("[metal]", 1)
                .unwrap_err()
                .info
                .as_deref(),
            Some(" in [...]")
        );
        // Empty member list.
        assert_eq!(
            parse_cid_elements("[]", 1).unwrap_err().info.as_deref(),
            Some(" in [...]")
        );
        // A missing closing bracket cannot loop (source would; callers
        // guarantee it) — the port fails closed with the in-list note.
        assert_eq!(
            parse_cid_elements("[C", 1).unwrap_err().info.as_deref(),
            Some(" in [...]")
        );
        // Wildcard returns no mask at all.
        assert!(parse_cid_elements("[*]", 1).unwrap().is_none());
    }
    #[test]
    fn bio_sel_b10_relation_tokens_by_position() {
        use super::has_inequality;
        let cid = "ab<c=de>f";
        // Every relation byte position is recognized inside the window.
        for (start, end, expected) in [
            (0usize, 9usize, true),
            (0, 3, true),
            (3, 4, false),
            (4, 6, true),
            (7, 9, true),
            (0, 2, false),
            (4, 5, true),
            (5, 7, false),
        ] {
            assert_eq!(has_inequality(cid, start, end), expected, "{start}..{end}");
        }
        // A window that clips a relation byte excludes it (end is exclusive).
        assert!(!has_inequality(cid, 0, 2));
        assert!(has_inequality(cid, 0, 3));
    }

    #[test]
    fn bio_sel_b10_missing_relation_whitespace_and_extensions() {
        use super::has_inequality;
        // No relation byte at all -> false, even with other punctuation.
        assert!(!has_inequality("q 0.5", 0, 5));
        assert!(!has_inequality("polymer", 0, 7));
        assert!(!has_inequality("!A,B", 0, 4));
        // Whitespace around a relation still counts (the byte is present).
        assert!(has_inequality("b < 1", 0, 5));
        assert!(has_inequality("q=0.5", 0, 5));
        assert!(has_inequality("b>2 ;", 0, 4));
        // A window beyond the string end is simply empty (no panic).
        assert!(!has_inequality("abc", 1, 99));
    }

    // BIO-SEL-FIX Step 12 regressions relocated with the owner (BIO-CID
    // C02); expected byte vectors are composed independently from
    // select.cpp:14-26 — never from the Display implementation.
    #[test]
    fn bio_sel_b05_source_message_bytes_ascii_and_offsets() {
        use super::wrong_syntax;
        // Plain ASCII middle window, exactly eight bytes.
        let error = wrong_syntax("ABCDEFGHIJKLMN", 3, None);
        assert_eq!(error.pos, 3);
        let expected: Vec<u8> =
            b"Invalid selection syntax near \"DEFGHIJK\": ABCDEFGHIJKLMN".to_vec();
        assert_eq!(error.source_message_bytes(), expected);
        // Remainder shorter than eight bytes: the window clamps.
        let error = wrong_syntax("ABCDEF", 4, None);
        let expected: Vec<u8> = b"Invalid selection syntax near \"EF\": ABCDEF".to_vec();
        assert_eq!(error.source_message_bytes(), expected);
        // pos = 0 suppresses the near section entirely.
        let error = wrong_syntax("ABC", 0, None);
        let expected: Vec<u8> = b"Invalid selection syntax: ABC".to_vec();
        assert_eq!(error.source_message_bytes(), expected);
        // Empty info vs. present info: the note inserts verbatim.
        let error = wrong_syntax("ABC", 0, Some(" (at residue)"));
        let expected: Vec<u8> = b"Invalid selection syntax (at residue): ABC".to_vec();
        assert_eq!(error.source_message_bytes(), expected);
        // pos == len yields an empty window (legal substr), no panic.
        let error = wrong_syntax("ABC", 3, None);
        let expected: Vec<u8> = b"Invalid selection syntax near \"\": ABC".to_vec();
        assert_eq!(error.source_message_bytes(), expected);
        // The full original cid is retained even with multibyte content.
        let error = wrong_syntax("ÄÖ", 0, None);
        let mut expected: Vec<u8> = b"Invalid selection syntax: ".to_vec();
        expected.extend_from_slice("ÄÖ".as_bytes());
        assert_eq!(error.source_message_bytes(), expected);
        assert_eq!(error.cid, "ÄÖ");
    }

    #[test]
    fn bio_sel_b05_source_message_bytes_multibyte_window_cuts() {
        use super::wrong_syntax;
        // An eight-byte window starting inside multibyte content cuts
        // exactly eight raw bytes, splitting characters as the source's
        // byte substr does: "ÄÖÜ" = C3 84 C3 96 C3 9C (6 bytes); pos = 1
        // yields 84 C3 96 C3 9C + ... use a longer cid so the window is
        // full: "xÄÖÜÄÖÜÄÖÜ" with pos = 1 selects bytes 1..9 =
        // C3 84 C3 96 C3 9C C3 84 (four bytes of one char pair boundary
        // cut mid-character at both ends is impossible; here the leading
        // 'x' makes the start a cut).
        let cid = "xÄÖÜÄÖÜÄÖÜ";
        let error = wrong_syntax(cid, 1, None);
        let mut expected: Vec<u8> = b"Invalid selection syntax near \"".to_vec();
        expected.extend_from_slice(&cid.as_bytes()[1..9]);
        expected.extend_from_slice(b"\": ");
        expected.extend_from_slice(cid.as_bytes());
        assert_eq!(error.source_message_bytes(), expected);
        // pos at the start of a multibyte character (pos = 2, the first
        // byte of "Ö"): the eight-byte window [2, 10) covers exactly four
        // two-byte characters and stays boundary-aligned — contrast with
        // the pos = 1 window above, which starts mid-character. (pos = 0
        // would suppress the near section entirely.)
        let cid2 = "ÄÖÜÄÖÜÄÖÜ";
        let error2 = wrong_syntax(cid2, 2, Some(" in [...]"));
        let mut expected2: Vec<u8> = b"Invalid selection syntax in [...] near \"".to_vec();
        expected2.extend_from_slice(&cid2.as_bytes()[2..10]);
        expected2.extend_from_slice(b"\": ");
        expected2.extend_from_slice(cid2.as_bytes());
        assert_eq!(error2.source_message_bytes(), expected2);
        // A boundary-aligned eight-byte window renders identically in
        // Display (no lossy replacement needed).
        assert!(!error2.to_string().contains('\u{FFFD}'));
        assert_eq!(
            error2.to_string().as_bytes(),
            error2.source_message_bytes().as_slice()
        );
        // The byte window and the lossy Display differ exactly on split
        // characters (U+FFFD): pos = 2 starts inside the first "Ä" (a
        // continuation byte), so both window ends cut mid-character —
        // Display is presentation only, never byte-equivalent there.
        let error_cut = wrong_syntax(cid, 2, None);
        let mut expected_cut: Vec<u8> = b"Invalid selection syntax near \"".to_vec();
        expected_cut.extend_from_slice(&cid.as_bytes()[2..10]);
        expected_cut.extend_from_slice(b"\": ");
        expected_cut.extend_from_slice(cid.as_bytes());
        assert_eq!(error_cut.source_message_bytes(), expected_cut);
        let rendered = error_cut.to_string();
        assert!(rendered.contains('\u{FFFD}'));
        assert_ne!(
            rendered.as_bytes(),
            error_cut.source_message_bytes().as_slice()
        );
    }

    // ---- BIO-CID C03: determine_omitted_cid_fields closure. No change to
    // the relocated body was needed (receipt no-change evidence): every
    // initial-character class and the FIRST-separator precedence below are
    // derived independently from select.cpp:24-41. ----

    #[test]
    fn bio_cid_c03_initial_character_classes_and_delimiter_precedence() {
        use super::OmittedCidFields::*;
        use super::determine_omitted_cid_fields as omit;
        // Model: exactly a leading '/' (select.cpp:25-26).
        for cid in ["/", "/1/A", "/x"] {
            assert_eq!(omit(cid), Model, "{cid}");
        }
        // Residue by leading class: EVERY ASCII digit, '.', '(' and '-'
        // (select.cpp:28-29).
        for digit in b'0'..=b'9' {
            let cid = format!("{}A", digit as char);
            assert_eq!(omit(&cid), Residue, "digit {:?}", digit as char);
        }
        for lead in ['.', '(', '-'] {
            let cid = format!("{lead}A");
            assert_eq!(omit(&cid), Residue, "{cid}");
        }
        // Chain: no separator at all, or the FIRST separator is '/' or ';'
        // (select.cpp:30-32).
        for cid in ["A", "abc", "ABC", "+x", " x", "A/1", "AB;C"] {
            assert_eq!(omit(cid), Chain, "{cid}");
        }
        // Atom: the FIRST separator is '[' or ':' (select.cpp:35-36).
        for cid in ["CA[C]", "CA:C", "C[A]:x", "C:A(x)"] {
            assert_eq!(omit(cid), Atom, "{cid}");
        }
        // Delimiter precedence: the FIRST separator decides — '(' later
        // than '[' or ':' is irrelevant, and vice versa.
        assert_eq!(omit("A(B)[C]"), Residue); // '(' first
        assert_eq!(omit("A[B](C)"), Atom); // '[' first
        assert_eq!(omit("A:B(C)"), Atom); // ':' first
        assert_eq!(omit("A(B):C"), Residue); // '(' first
        // '/' dominating later '(' — still chain.
        assert_eq!(omit("A/1(B)"), Chain);
        // Empty cid: operator[] at size() yields '\0'; no branch matches,
        // no separator exists → chain (source-exact).
        assert_eq!(omit(""), Chain);
    }

    #[test]
    fn bio_cid_c03_multibyte_leading_names() {
        use super::OmittedCidFields::*;
        use super::determine_omitted_cid_fields as omit;
        // A valid multibyte (non-ASCII) leading name byte matches none of
        // the leading tests (byte values outside the compared ASCII
        // ranges), so classification comes from the first ASCII separator
        // exactly like the source's byte-level find_first_of.
        assert_eq!(omit("\u{3a9}"), Chain); // "\u{3a9}" alone: no separator
        assert_eq!(omit("\u{3a9}\u{3b1}"), Chain); // two multibyte chars
        assert_eq!(omit("\u{3a9}/1"), Chain);
        assert_eq!(omit("\u{3a9};x"), Chain);
        assert_eq!(omit("\u{3a9}(ALA)"), Residue);
        assert_eq!(omit("\u{3a9}[C]"), Atom);
        assert_eq!(omit("\u{3a9}:C"), Atom);
        // Mixed: multibyte member AFTER an ASCII first byte still uses the
        // first separator; the multibyte bytes never match separators.
        assert_eq!(omit("A\u{3a9}[C]"), Atom);
        assert_eq!(omit("A\u{3a9}(B)"), Residue);
    }

    // ---- BIO-CID C04: make_cid_list closure. No change to the relocated
    // body (receipt no-change evidence): substr clamp/wrap, one-byte
    // wildcard/inversion consumption and the composed punctuation note are
    // all source-exact. Expectations below derive from select.cpp:43-53. ----

    #[test]
    fn bio_cid_c04_every_default_disallowed_punctuation_errors() {
        use super::make_cid_list;
        // The DEFAULT disallowed set "-[]()!/*.:;" (select.cpp:45): each
        // single member made only of that byte fails at its own position
        // with the composed note. Positions: '!' and '*' are CONSUMED first
        // (source ++pos), so a lone '!'/'*' member reports the SUBSTR
        // content (empty), i.e. no offender — the source reads cid[pos]
        // BEFORE the increment for the flags, but the punctuation scan runs
        // on the payload; verify the real behavior instead: use two-byte
        // members where the offender survives consumption.
        for (offender, cid, pos, end, expected_pos) in [
            ('-', "A-B", 0usize, 3usize, 1usize),
            ('[', "A[B", 0, 3, 1),
            (']', "A]B", 0, 3, 1),
            ('(', "A(B", 0, 3, 1),
            (')', "A)B", 0, 3, 1),
            ('!', "A!B", 0, 3, 1),
            ('/', "A/B", 0, 3, 1),
            ('*', "A*B", 0, 3, 1),
            ('.', "A.B", 0, 3, 1),
            (':', "A:B", 0, 3, 1),
            (';', "A;B", 0, 3, 1),
        ] {
            let error = make_cid_list(cid, pos, end, "-[]()!/*.:;").unwrap_err();
            assert_eq!(error.pos, expected_pos, "{cid}");
            assert_eq!(error.cid, cid);
            assert_eq!(
                error.info.as_deref(),
                Some(format!(" ('{offender}' in a list)").as_str()),
                "{cid}"
            );
        }
        // A member of ONLY '!' or '*': the flag byte is consumed and the
        // payload is empty — no punctuation error (select.cpp:48-51 order:
        // flags read at pos, THEN ++pos, THEN substr of the remainder).
        let bang = make_cid_list("!", 0, 1, "-[]()!/*.:;").unwrap();
        assert!(bang.inverted && bang.list.is_empty());
        let star = make_cid_list("*", 0, 1, "-[]()!/*.:;").unwrap();
        assert!(star.all && star.list.is_empty());
        // Inverted member with a NON-offender payload: no error.
        let ok = make_cid_list("!AB", 0, 3, "-[]()!/*.:;").unwrap();
        assert!(ok.inverted && ok.list == b"AB");
    }

    #[test]
    fn bio_cid_c04_chain_stage_dash_exception_and_first_offender_wins() {
        use super::make_cid_list;
        // The chain-stage caller passes the default set WITHOUT '-' so a
        // literal '-' chain name is legal (select.cpp:196-197 stage passes
        // "[]()!/*.:;"); '-' then never errors under that set.
        let chain_set = "[]()!/*.:;";
        let dash = make_cid_list("-", 0, 1, chain_set).unwrap();
        assert_eq!(dash.list, b"-");
        assert!(!dash.all && !dash.inverted);
        // Every OTHER member of that reduced set still errors.
        for (offender, cid) in [('(', "A(B)"), (')', "A)"), ('!', "AB!C"), ('*', "AB*C")] {
            let error = make_cid_list(cid, 0, cid.len(), chain_set).unwrap_err();
            assert_eq!(error.pos, cid.find(offender).unwrap(), "{cid}");
            assert_eq!(
                error.info.as_deref(),
                Some(format!(" ('{offender}' in a list)").as_str()),
                "{cid}"
            );
        }
        // FIRST offender wins under either set (find_first_of).
        let first = make_cid_list("A.B[C", 0, 5, "-[]()!/*.:;").unwrap_err();
        assert_eq!(first.pos, 1);
        assert_eq!(first.info.as_deref(), Some(" ('.' in a list)"));
        // Commas are NEVER in any source set: comma members survive.
        let commas = make_cid_list("A,B,C", 0, 5, chain_set).unwrap();
        assert_eq!(commas.list, b"A,B,C");
    }

    #[test]
    fn bio_cid_c04_canonical_conversion_boundary_exact_utf8_only() {
        use super::make_cid_list;
        use cosmolkit_bio::{BioSelectionData, BioSelectionParts};
        // Caller-reachable fields are valid UTF-8; the exact-bytes carrier
        // converts losslessly through BioSelectionData::from_parts (C01):
        // multibyte members preserved, empty members preserved.
        let field = make_cid_list("/x/\u{3a9},A/y", 3, 7, "-[]()!/*.:;").unwrap();
        let parts = BioSelectionParts {
            chain_ids: field.clone(),
            ..BioSelectionParts::default()
        };
        let data = BioSelectionData::from_parts(parts).expect("valid UTF-8");
        assert_eq!(data.parts().chain_ids.list, "\u{3a9},A");
        // Empty member field converts to an empty String member list.
        let empty = make_cid_list("AB", 0, 0, "-[]()!/*.:;").unwrap();
        let parts = BioSelectionParts {
            residue_names: empty,
            ..BioSelectionParts::default()
        };
        let data = BioSelectionData::from_parts(parts).expect("empty is valid");
        assert!(data.parts().residue_names.list.is_empty());
        // Impossible-from-callers byte cut (helper-level substr fidelity)
        // stays exact bytes in the carrier and fails CLOSED at conversion —
        // never lossy, never an empty substitute.
        let cut = make_cid_list("\u{3a9}B", 0, 1, "-[]()!/*.:;").unwrap();
        assert_eq!(cut.list, b"\xCE");
        let parts = BioSelectionParts {
            atom_names: cut,
            ..BioSelectionParts::default()
        };
        // BioSelectionData has no Debug (private engine state), so the
        // fail-closed conversion is matched instead of unwrap_err.
        let error = match BioSelectionData::from_parts(parts) {
            Ok(_) => panic!("cut multibyte bytes must fail closed"),
            Err(error) => error,
        };
        assert_eq!(error.field, "atom_names");
        assert_eq!(error.bytes, vec![0xCE]);
    }

    // ---- BIO-CID C05: the shared checked decimal integer-prefix owner.
    // Expectations derive from the two pinned strtol call sites
    // (select.cpp:101/158) and C strtol semantics on LP64, never from the
    // implementation. ----

    #[test]
    fn bio_cid_c05_signs_whitespace_no_conversion_and_offsets() {
        use super::scan_cid_int_prefix;
        // Plain digits: value + full extent.
        assert_eq!(
            scan_cid_int_prefix(b"14/"),
            CidIntScan {
                value: 14,
                consumed: 2,
                converted: true
            }
        );
        // Signs.
        assert_eq!(
            scan_cid_int_prefix(b"-7x"),
            CidIntScan {
                value: -7,
                consumed: 2,
                converted: true
            }
        );
        assert_eq!(
            scan_cid_int_prefix(b"+3:"),
            CidIntScan {
                value: 3,
                consumed: 2,
                converted: true
            }
        );
        // Leading is_space bytes (atox.hpp 9-13 + 32) are consumed AND
        // counted when a conversion happens.
        for (field, consumed) in [
            (&b" 5"[..], 2usize),
            (&b"\t5"[..], 2),
            (&b"\n\r\x0b\x0c 5"[..], 6),
        ] {
            assert_eq!(
                scan_cid_int_prefix(field),
                CidIntScan {
                    value: 5,
                    consumed,
                    converted: true
                },
                "{field:?}"
            );
        }
        // Whitespace BETWEEN sign and digits is NOT skipped (strtol stops).
        assert_eq!(
            scan_cid_int_prefix(b"- 5"),
            CidIntScan {
                value: 0,
                consumed: 0,
                converted: false
            }
        );
        // No digits at all: no conversion — value 0, end pointer at the
        // ORIGINAL start (whitespace/sign never count alone).
        for field in [
            &b"x"[..],
            &b""[..],
            &b" "[..],
            &b"-"[..],
            &b"+"[..],
            &b" -"[..],
            &b" /A"[..],
        ] {
            assert_eq!(
                scan_cid_int_prefix(field),
                CidIntScan {
                    value: 0,
                    consumed: 0,
                    converted: false
                },
                "{field:?}"
            );
        }
        // '0' alone converts with extent 1 (unlike a lone sign).
        assert_eq!(
            scan_cid_int_prefix(b"0"),
            CidIntScan {
                value: 0,
                consumed: 1,
                converted: true
            }
        );
    }

    #[test]
    fn bio_cid_c05_i32_limits_long_saturation_and_huge_runs() {
        use super::scan_cid_int_prefix;
        // i32 boundary values are representable longs.
        assert_eq!(
            scan_cid_int_prefix(b"2147483647/"),
            CidIntScan {
                value: 2147483647,
                consumed: 10,
                converted: true
            }
        );
        assert_eq!(
            scan_cid_int_prefix(b"-2147483648/"),
            CidIntScan {
                value: -2147483648,
                consumed: 11,
                converted: true
            }
        );
        // Just outside i32: still a valid LONG — the owner reports the real
        // long value; the typed i32 error belongs to the callers (C06/C09).
        assert_eq!(scan_cid_int_prefix(b"2147483648;").value, 2147483648);
        assert_eq!(scan_cid_int_prefix(b"-2147483649;").value, -2147483649);
        // Exact long bounds: LONG_MAX = 2^63-1 and LONG_MIN = -2^63 are NOT
        // saturation; one digit beyond saturates to LONG_MAX/LONG_MIN while
        // STILL consuming every digit (extent = full run).
        assert_eq!(
            scan_cid_int_prefix(b"9223372036854775807x"),
            CidIntScan {
                value: i64::MAX,
                consumed: 19,
                converted: true
            }
        );
        assert_eq!(
            scan_cid_int_prefix(b"-9223372036854775808x"),
            CidIntScan {
                value: i64::MIN,
                consumed: 20,
                converted: true
            }
        );
        assert_eq!(
            scan_cid_int_prefix(b"9223372036854775808x"),
            CidIntScan {
                value: i64::MAX,
                consumed: 19,
                converted: true
            }
        );
        assert_eq!(
            scan_cid_int_prefix(b"-9223372036854775809x"),
            CidIntScan {
                value: i64::MIN,
                consumed: 20,
                converted: true
            }
        );
        // 40-digit and 400-digit runs with a suffix: saturation + FULL
        // consumed extent, no internal overflow anywhere.
        let forty = format!("{};q>2", "9".repeat(40));
        let scan = scan_cid_int_prefix(forty.as_bytes());
        assert_eq!(scan.value, i64::MAX);
        assert_eq!(scan.consumed, 40);
        let four_hundred = format!("{}x", "7".repeat(400));
        let scan = scan_cid_int_prefix(four_hundred.as_bytes());
        assert_eq!(scan.value, i64::MAX);
        assert_eq!(scan.consumed, 400);
        // Negative huge run saturates to LONG_MIN with full extent.
        let neg = format!("-{}y", "1".repeat(100));
        let scan = scan_cid_int_prefix(neg.as_bytes());
        assert_eq!(scan.value, i64::MIN);
        assert_eq!(scan.consumed, 101);
        // A zero-padded in-range value keeps its real magnitude.
        assert_eq!(scan_cid_int_prefix(b"0000000012;").value, 12);
        assert_eq!(scan_cid_int_prefix(b"0000000012;").consumed, 10);
    }

    // ---- BIO-CID C06: parse_cid_seqid routed through the C05 owner.
    // Expectations from select.cpp:95-115 with the packet-mandated typed
    // i32 range error; existing b08 coverage retained unchanged. ----

    #[test]
    fn bio_cid_c06_owner_routing_and_typed_range_errors() {
        use super::parse_cid_seqid;
        // Ordinary routing stays identical: value, icode and advanced pos.
        let mut pos = 0usize;
        let plain = parse_cid_seqid("1234A", &mut pos, i32::MIN).unwrap();
        assert_eq!((plain.seqnum, plain.icode, pos), (1234, b'A', 5));
        // Long-digit overflow beyond i32 (but within long): typed range
        // error whose text spans the FULL consumed digit run.
        let mut pos = 0usize;
        let error = parse_cid_seqid("99999999999x", &mut pos, i32::MIN).unwrap_err();
        assert_eq!(error.text, "99999999999");
        let mut pos = 0usize;
        let error = parse_cid_seqid("-3000000000", &mut pos, i32::MAX).unwrap_err();
        assert_eq!(error.text, "-3000000000");
        // A saturated-long 40-digit run still reports the whole run as the
        // offending text (strtol consumed every digit) — no truncation.
        let huge = format!("{};q>2", "9".repeat(40));
        let mut pos = 0usize;
        let error = parse_cid_seqid(&huge, &mut pos, 0).unwrap_err();
        assert_eq!(error.text, "9".repeat(40));
        // i32 bounds themselves are representable (no error).
        let mut pos = 0usize;
        let max = parse_cid_seqid("2147483647/", &mut pos, 0).unwrap();
        assert_eq!((max.seqnum, pos), (i32::MAX, 10));
        let mut pos = 0usize;
        let min = parse_cid_seqid("-2147483648/", &mut pos, 0).unwrap();
        assert_eq!((min.seqnum, pos), (i32::MIN, 11));
        // One beyond each bound: typed error with exact text.
        let mut pos = 0usize;
        assert_eq!(
            parse_cid_seqid("2147483648", &mut pos, 0).unwrap_err().text,
            "2147483648"
        );
        let mut pos = 0usize;
        assert_eq!(
            parse_cid_seqid("-2147483649", &mut pos, 0)
                .unwrap_err()
                .text,
            "-2147483649"
        );
        // Minus-without-digits: no conversion, seqnum 0, pos unchanged —
        // the caller still sees the '-'.
        let mut pos = 0usize;
        let dash = parse_cid_seqid("-x", &mut pos, i32::MIN).unwrap();
        assert_eq!((dash.seqnum, dash.icode, pos), (0, b' ', 0));
        // Default sentinels preserved when nothing converts.
        let mut pos = 0usize;
        let untouched = parse_cid_seqid("/A", &mut pos, i32::MAX).unwrap();
        assert_eq!(untouched.seqnum, i32::MAX);
        assert_eq!(pos, 0);
        // Wildcard '*' keeps the default seqnum and takes icode '*'; a
        // following letter OVERWRITES it (source assigns unconditionally).
        let mut pos = 0usize;
        let wildcard = parse_cid_seqid("*a", &mut pos, 7).unwrap();
        assert_eq!((wildcard.seqnum, wildcard.icode, pos), (7, b'a', 2));
        let mut pos = 0usize;
        let star = parse_cid_seqid("*", &mut pos, 7).unwrap();
        assert_eq!((star.seqnum, star.icode, pos), (7, b'*', 1));
    }

    // ---- BIO-CID C07: parse_cid_elements closure via the canonical BIO
    // boundary (no copied table; imports GEMMI_EL_END/GEMMI_IS_METAL/
    // gemmi_find_element). No change to the relocated body. Expectations
    // from select.cpp:59-93 + pinned elem.hpp. ----

    #[test]
    fn bio_cid_c07_inversions_groups_and_fixed_spellings() {
        use super::parse_cid_elements;
        // Both inversions x both group expansions, checked against the
        // canonical table at every ordinal (metals/nonmetals under '!').
        for (cid, expect_metal) in [("metals", true), ("!metals", false)] {
            let mask = parse_cid_elements(&format!("[{cid}]"), 1).unwrap().unwrap();
            for ordinal in 0..super::GEMMI_EL_END {
                assert_eq!(
                    mask[ordinal],
                    expect_metal == super::GEMMI_IS_METAL[ordinal],
                    "{cid} at {ordinal}"
                );
            }
        }
        // Nonmetals = complement at every ordinal, both inversions.
        for (cid, expect_metal) in [("nonmetals", false), ("!nonmetals", true)] {
            let mask = parse_cid_elements(&format!("[{cid}]"), 1).unwrap().unwrap();
            for ordinal in 0..super::GEMMI_EL_END {
                assert_eq!(
                    mask[ordinal],
                    expect_metal == super::GEMMI_IS_METAL[ordinal],
                    "{cid} at {ordinal}"
                );
            }
        }
        // Empty and degenerate groups: empty member, inversion alone,
        // inversion with an empty member — all fail with the in-list note.
        for cid in ["[]", "[!]", "[!,]", "[,]"] {
            assert_eq!(
                parse_cid_elements(cid, 1).unwrap_err().info.as_deref(),
                Some(" in [...]"),
                "{cid}"
            );
        }
        // Fixed spellings are CASE-SENSITIVE exact tokens: only the exact
        // lower-case "metals"/"nonmetals" expand; every other spelling is
        // an invalid element / in-list error, never a group.
        for cid in [
            "[METALS]",
            "[Metal]",
            "[metals ]",
            "[ metals]",
            "[nonMetals]",
        ] {
            assert!(
                parse_cid_elements(cid, 1).is_err(),
                "{cid} must not be a group"
            );
        }
        // D ordinal 119 resolves through the canonical table; wildcard
        // yields no mask; mixed members accumulate across groups.
        let mask = parse_cid_elements("[D,metals,C]", 1).unwrap().unwrap();
        assert!(mask[119] && mask[6]);
        assert_eq!(mask[1], super::GEMMI_IS_METAL[1]); // H is not a metal
        assert!(mask[3]); // Li is a metal via the group
        assert!(parse_cid_elements("[*]", 1).unwrap().is_none());
    }

    #[test]
    fn bio_cid_c06_no_conversion_matrix() {
        use super::parse_cid_seqid;
        // 36-row matrix, independently derived from select.cpp:92-106:
        // the digit-or-'-' guard ENTERS the strtol branch for a lone minus,
        // strtol then performs NO conversion (value 0, endptr == start), so
        // pos is unchanged, the '.' guard reads '-' (not '.'), the
        // insertion-code guard sees initial_pos == pos and never fires.
        // Six suffix forms x two positions x three defaults.
        let suffixes = ["-", "-A", "- ", "-/", "-.", "-.A"];
        let starts = [("", 0usize), ("/x/", 3usize)];
        let defaults = [i32::MIN, i32::MAX, 7];
        for (prefix, start) in starts {
            for suffix in suffixes {
                for default in defaults {
                    let cid = format!("{prefix}{suffix}");
                    let mut pos = start;
                    let parsed = parse_cid_seqid(&cid, &mut pos, default).unwrap();
                    assert_eq!(
                        (parsed.seqnum, parsed.icode, pos),
                        (0, b' ', start),
                        "{cid:?} default {default}"
                    );
                }
            }
        }
    }

    // ---- BIO-CID C08: parse_atom_inequality through the shared
    // bio_numeric C-string owner. Expectations derived from
    // select.cpp:107-141 + atof.hpp:24-30, never the implementation. ----

    #[test]
    fn bio_cid_c08_properties_relations_and_finite_values() {
        // q/b x all three relations x finite tails; pos advances exactly
        // to end on success.
        for (cid, property, relation, value) in [
            ("q<3.5", b'q', -1, 3.5),
            ("b<2", b'b', -1, 2.0),
            ("q=2", b'q', 0, 2.0),
            ("b=0.25", b'b', 0, 0.25),
            ("q>-1", b'q', 1, -1.0),
            ("b>1e3", b'b', 1, 1000.0),
        ] {
            let mut pos = 0usize;
            let parsed = parse_atom_inequality(cid, &mut pos, cid.len()).unwrap();
            assert_eq!(parsed.property, property, "{cid}");
            assert_eq!(parsed.relation, relation, "{cid}");
            assert_eq!(parsed.value, value, "{cid}");
            assert_eq!(pos, cid.len(), "{cid}");
        }
        // Relation trimming: spaces around the relation and the number.
        let mut pos = 0usize;
        let spaced = parse_atom_inequality("q  >  3", &mut pos, 7).unwrap();
        assert_eq!(
            (spaced.property, spaced.relation, spaced.value, pos),
            (b'q', 1, 3.0, 7)
        );
        // Boundary consumption stops at end (semicolon left for caller).
        let mut pos = 0usize;
        let bounded = parse_atom_inequality("q>1;q<2", &mut pos, 3).unwrap();
        assert_eq!((bounded.relation, bounded.value, pos), (1, 1.0, 3));
        // Signed zero keeps its sign.
        let mut pos = 0usize;
        let zero = parse_atom_inequality("b<-0", &mut pos, 4).unwrap();
        assert_eq!(zero.value, 0.0);
        assert!(zero.value.is_sign_negative());
    }

    #[test]
    fn bio_cid_c08_special_literals_and_range_rejections() {
        // Literal specials parse successfully (fast_float accepts
        // case-insensitive inf/infinity/nan with optional sign).
        let mut pos = 0usize;
        let inf = parse_atom_inequality("q>inf", &mut pos, 5).unwrap();
        assert_eq!(inf.value, f64::INFINITY);
        let mut pos = 0usize;
        let neg_inf = parse_atom_inequality("b=-Infinity", &mut pos, 11).unwrap();
        assert_eq!(neg_inf.value, f64::NEG_INFINITY);
        let mut pos = 0usize;
        let nan = parse_atom_inequality("q=nan", &mut pos, 5).unwrap();
        assert!(nan.value.is_nan());
        // Zero mantissa with huge exponent (0e999) stays representable:
        // fast_float's range test needs a NONZERO mantissa rounding to
        // zero (N11 row), so it parses to 0.0 successfully.
        let mut pos = 0usize;
        let zero_exp = parse_atom_inequality("q=0e999", &mut pos, 7).unwrap();
        assert_eq!((zero_exp.value, pos), (0.0, 7));
        // Nonzero-mantissa over/underflow: ec != errc() -> same
        // " (expected number)" note AT THE NUMBER POSITION even though
        // fast_float assigns a saturated value.
        for cid in ["q>1e999", "b<1e-400"] {
            let mut pos = 0usize;
            let error = parse_atom_inequality(cid, &mut pos, cid.len()).unwrap_err();
            assert_eq!(error.pos, 2, "{cid}");
            assert_eq!(error.info.as_deref(), Some(" (expected number)"), "{cid}");
        }
        // Invalid numbers: same note, number position.
        for cid in ["q>x", "b<", "q= ", "q>"] {
            let mut pos = 0usize;
            let error = parse_atom_inequality(cid, &mut pos, cid.len()).unwrap_err();
            assert_eq!(error.pos, 2, "{cid}");
            assert_eq!(error.info.as_deref(), Some(" (expected number)"), "{cid}");
        }
    }

    #[test]
    fn bio_cid_c08_grammar_errors_and_trailing_tokens() {
        // Wrong property at pos.
        let mut pos = 0usize;
        assert_eq!(
            parse_atom_inequality("x>1", &mut pos, 3).unwrap_err().pos,
            0
        );
        // Empty payload.
        let mut pos = 0usize;
        assert_eq!(parse_atom_inequality("", &mut pos, 0).unwrap_err().pos, 0);
        // Wrong relation byte at its own position (spaces skipped first).
        for (cid, bad) in [("q?1", 1usize), ("q 1", 2), ("q?", 1)] {
            let mut pos = 0usize;
            assert_eq!(
                parse_atom_inequality(cid, &mut pos, cid.len())
                    .unwrap_err()
                    .pos,
                bad,
                "{cid}"
            );
        }
        // Trailing junk inside end: number consumed, no space, pos != end.
        let mut pos = 0usize;
        let junk = parse_atom_inequality("q=1x", &mut pos, 4).unwrap_err();
        assert_eq!((junk.pos, junk.info.as_deref()), (3, None));
        // Trailing byte after spaces.
        let mut pos = 0usize;
        let spaced_junk = parse_atom_inequality("q=1 1", &mut pos, 5).unwrap_err();
        assert_eq!(spaced_junk.pos, 4);
        // Trailing spaces WITHIN end are consumed and accepted.
        let mut pos = 0usize;
        let ok = parse_atom_inequality("q=1  ", &mut pos, 5).unwrap();
        assert_eq!((ok.value, pos), (1.0, 5));
        let mut pos = 0usize;
        let one_space = parse_atom_inequality("q=1 ", &mut pos, 4).unwrap();
        assert_eq!(pos, 4);
        // End beyond the consumed payload is rejected at the rest position.
        let mut pos = 0usize;
        assert_eq!(
            parse_atom_inequality("q>1", &mut pos, 6).unwrap_err().pos,
            3
        );
    }

    // ---- BIO-CID C09: model stage. Expectations derived from
    // select.cpp:146-163 (omit dispatch, sep/semi bounds, sep==1 and
    // cid[1]=='*' bypasses, strtol end-pointer validation at pos 0) and
    // the C03 omit classes. ----

    #[test]
    fn bio_cid_c09_preamble_bypass_and_omission_dispatch() {
        // Empty or lone '*' bypasses every stage.
        assert!(cid_selection_bypassed(""));
        assert!(cid_selection_bypassed("*"));
        assert!(!cid_selection_bypassed("/"));
        assert!(!cid_selection_bypassed("**"));
        assert!(!cid_selection_bypassed("/1"));
        // omit != 0 (chain/residue/atom-first CID): stage skipped, mdl
        // stays the Selection default 0 (all) and sep stays 0.
        for (cid, omit) in [("A/1-5", 1), ("1-5", 2), ("(LYS)", 2), ("q>1", 3)] {
            assert_eq!(
                parse_cid_model_stage(cid, omit, None).unwrap(),
                (0, Some(0)),
                "{cid}"
            );
        }
    }

    #[test]
    fn bio_cid_c09_model_numbers_and_bounds() {
        // Leading '/' starts the model field: number, then '/' or ';'.
        assert_eq!(
            parse_cid_model_stage("/1/A", 0, None).unwrap(),
            (1, Some(2))
        );
        assert_eq!(
            parse_cid_model_stage("/12/A", 0, None).unwrap(),
            (12, Some(3))
        );
        assert_eq!(
            parse_cid_model_stage("/0/A", 0, None).unwrap(),
            (0, Some(2))
        );
        // Signed model runs to the string end (no separator survives).
        assert_eq!(parse_cid_model_stage("/-2", 0, None).unwrap(), (-2, None));
        assert_eq!(
            parse_cid_model_stage("/1;", 0, Some(2)).unwrap(),
            (1, Some(2))
        );
        // No-conversion tail: mdl 0, end_pos 1; "/" alone satisfies
        // end_pos == cid.size().
        assert_eq!(parse_cid_model_stage("/", 0, None).unwrap(), (0, None));
        // Empty model field ("//...") and wildcard ("/*...") bypass the
        // conversion with mdl 0 = all models; sep is still reported.
        assert_eq!(parse_cid_model_stage("//A", 0, None).unwrap(), (0, Some(1)));
        assert_eq!(
            parse_cid_model_stage("/*/A", 0, None).unwrap(),
            (0, Some(2))
        );
    }

    #[test]
    fn bio_cid_c09_model_errors_and_range() {
        // Invalid trailing bytes: end_pos equals neither sep nor size ->
        // the fixed note at position 0.
        for cid in ["/1x/A", "/1 2/A", "/12x"] {
            match parse_cid_model_stage(cid, 0, None).unwrap_err() {
                BioSelectionParseError::Syntax(error) => {
                    assert_eq!(
                        (error.pos, error.info.as_deref()),
                        (0, Some(" (at model number)")),
                        "{cid}"
                    );
                }
                _ => panic!("{cid}: expected Syntax"),
            }
        }
        // "/A": no conversion, end_pos 1 matches neither npos nor size 2.
        match parse_cid_model_stage("/A", 0, None).unwrap_err() {
            BioSelectionParseError::Syntax(error) => {
                assert_eq!(
                    (error.pos, error.info.as_deref()),
                    (0, Some(" (at model number)"))
                );
            }
            other => panic!("expected Syntax, got {other:?}"),
        }
        // Beyond i32: typed range error with the full digit run as text
        // (packet divergence from long->int truncation).
        match parse_cid_model_stage("/5000000000/A", 0, None).unwrap_err() {
            BioSelectionParseError::SeqidRange(error) => assert_eq!(error.text, "5000000000"),
            _ => panic!("expected SeqidRange"),
        }
    }

    // ---- BIO-CID C10: chain stage via C04 make_cid_list. Expectations
    // from select.cpp:165-171: guard omit<=1 && sep<semi, pos restart
    // after the model separator, min(slash,semi) bound, disallowed set
    // WITHOUT '-' (RCSB bioassembly chain names). ----

    #[test]
    fn bio_cid_c10_guard_positions_and_names() {
        // omit 1 (chain-first): pos 0, list to the '/' bound.
        let (list, sep) = parse_cid_chain_stage("A,B/C", 1, None, Some(0)).unwrap();
        assert_eq!(
            (list.all, list.inverted, list.list.as_slice()),
            (false, false, b"A,B".as_slice())
        );
        assert_eq!(sep, Some(3));
        // No further '/': list runs to the semicolon or the end.
        let (list, sep) = parse_cid_chain_stage("A", 1, None, Some(0)).unwrap();
        assert_eq!(list.list, b"A");
        assert_eq!(sep, None);
        let (list, sep) = parse_cid_chain_stage("A;q>1", 1, Some(1), Some(0)).unwrap();
        assert_eq!(list.list, b"A");
        assert_eq!(sep, Some(1));
        // After the model stage (omit 0): pos = model sep + 1.
        let (list, sep) = parse_cid_chain_stage("/1/A,B/3", 0, None, Some(2)).unwrap();
        assert_eq!(list.list, b"A,B");
        assert_eq!(sep, Some(6));
        // '-' stays legal in chain names.
        let (list, _) = parse_cid_chain_stage("-A/1", 1, None, Some(0)).unwrap();
        assert_eq!(list.list, b"-A");
        // Multibyte UTF-8 chain names survive byte-exact.
        let (list, _) = parse_cid_chain_stage("\u{3a9},A/1", 1, None, Some(0)).unwrap();
        assert_eq!(list.list, "\u{3a9},A".as_bytes());
        // Empty member: "///1" -> model sep 1 -> pos 2 -> '/' bound 2.
        let (list, sep) = parse_cid_chain_stage("///1", 0, None, Some(1)).unwrap();
        assert_eq!((list.all, list.list.as_slice()), (false, b"".as_slice()));
        assert_eq!(sep, Some(2));
    }

    #[test]
    fn bio_cid_c10_wildcard_inversion_and_skips() {
        // Wildcard and inverted members consume their flag byte.
        let (list, _) = parse_cid_chain_stage("*/1", 1, None, Some(0)).unwrap();
        assert!(list.all && list.list.is_empty());
        let (list, _) = parse_cid_chain_stage("!A,B/1", 1, None, Some(0)).unwrap();
        assert!(list.inverted && list.list == b"A,B");
        // Guard skip: model sep None (npos) never satisfies sep < semi;
        // the pinned parse_cid leaves sep unchanged (identity-on-skip).
        let (list, sep) = parse_cid_chain_stage("/1", 0, None, None).unwrap();
        assert!(list.all);
        assert_eq!(sep, None);
        // Guard skip: residue/atom-first CIDs (omit 2/3) keep Some(0).
        for cid in ["1-5", "(LYS)", "q>1"] {
            let (list, sep) = parse_cid_chain_stage(cid, 2, None, Some(0)).unwrap();
            assert!(list.all && sep == Some(0), "{cid}");
        }
        // omit <= 1 but sep == semi: stage skipped; "/1;" keeps Some(2).
        let (list, sep) = parse_cid_chain_stage("/1;", 0, Some(2), Some(2)).unwrap();
        assert!(list.all);
        assert_eq!(sep, Some(2));
        // Punctuation (other than '-' and ',') inside the member errors
        // through the C04 boundary with its composed note.
        let error = parse_cid_chain_stage("A(B/1", 1, None, Some(0)).unwrap_err();
        match error {
            BioSelectionParseError::Syntax(error) => {
                assert_eq!(error.pos, 1);
                assert_eq!(error.info.as_deref(), Some(" ('(' in a list)"));
            }
            other => panic!("expected Syntax, got {other:?}"),
        }
    }

    #[test]
    fn bio_cid_c08_complete_property_relation_matrix() {
        // 360 REAL parse_atom_inequality calls: 6 property/relation
        // combos x 60 numeric-tail rows (10 tokens x 5 suffixes at
        // whole-string end + 10 token-boundary rows with ";q>2" tails and
        // field end at the token end). Expectations are literal, derived
        // from select.cpp:107-141 + the pinned fast_float semantics.
        let relations: [(u8, u8, i32); 3] = [(b'q', b'<', -1), (b'q', b'=', 0), (b'q', b'>', 1)];
        let relations: Vec<(u8, u8, i32)> = relations
            .into_iter()
            .chain([(b'b', b'<', -1), (b'b', b'=', 0), (b'b', b'>', 1)])
            .collect();
        // (token, expected value class). First six parse; last four
        // always fail numeric conversion with the expected-number note.
        let tokens: [(&str, Option<f64>); 10] = [
            ("1.5", Some(1.5)),
            ("-0", Some(-0.0)),
            ("inf", Some(f64::INFINITY)),
            ("-Infinity", Some(f64::NEG_INFINITY)),
            ("nan", Some(f64::NAN)),
            ("0e999", Some(0.0)),
            ("1e999", None),
            ("1e-400", None),
            ("x", None),
            ("", None),
        ];
        let suffixes = ["", "   ", "\t", "x", ";q>2"];
        let mut rows = 0usize;
        for (property, relation_char, relation) in relations {
            for (token, expected) in tokens {
                let token_len = token.len();
                // (a) whole-string field end rows.
                for suffix in suffixes {
                    let cid = format!(
                        "{}{}{}{}",
                        property as char, relation_char as char, token, suffix
                    );
                    let mut pos = 0usize;
                    let result = parse_atom_inequality(&cid, &mut pos, cid.len());
                    rows += 1;
                    match expected {
                        Some(value) if suffix.is_empty() || suffix == "   " => {
                            let parsed = result.unwrap_or_else(|e| panic!("{cid:?}: {e:?}"));
                            assert_eq!(parsed.property, property, "{cid:?}");
                            assert_eq!(parsed.relation, relation, "{cid:?}");
                            if value.is_nan() {
                                assert!(parsed.value.is_nan(), "{cid:?}");
                            } else {
                                assert_eq!(parsed.value, value, "{cid:?}");
                                assert_eq!(
                                    parsed.value.is_sign_negative(),
                                    value.is_sign_negative(),
                                    "{cid:?}"
                                );
                            }
                            assert_eq!(pos, cid.len(), "{cid:?}");
                        }
                        Some(_) => {
                            // tab/x/';q>2' after a valid token: the number
                            // consumed exactly the token, the trailing byte
                            // is not a space -> final-equality error with
                            // NO info at the consumed-token end; pos
                            // unchanged.
                            match result.unwrap_err() {
                                super::SelectionSyntaxError { pos: at, info, .. } => {
                                    assert_eq!(at, 2 + token_len, "{cid:?}");
                                    assert_eq!(info.as_deref(), None, "{cid:?}");
                                }
                                other => panic!("{cid:?}: {other:?}"),
                            }
                            assert_eq!(pos, 0, "{cid:?}");
                        }
                        None => {
                            match result.unwrap_err() {
                                super::SelectionSyntaxError { pos: at, info, .. } => {
                                    assert_eq!(at, 2, "{cid:?}");
                                    assert_eq!(
                                        info.as_deref(),
                                        Some(" (expected number)"),
                                        "{cid:?}"
                                    );
                                }
                                other => panic!("{cid:?}: {other:?}"),
                            }
                            assert_eq!(pos, 0, "{cid:?}");
                        }
                    }
                }
                // (b) token-boundary row: ";q>2" tail, field end at the
                // token end.
                let cid = format!("{}{}{};q>2", property as char, relation_char as char, token);
                let field_end = 2 + token_len;
                let mut pos = 0usize;
                let result = parse_atom_inequality(&cid, &mut pos, field_end);
                rows += 1;
                match expected {
                    Some(value) => {
                        let parsed = result.unwrap_or_else(|e| panic!("{cid:?}: {e:?}"));
                        assert_eq!(parsed.property, property, "{cid:?}");
                        assert_eq!(parsed.relation, relation, "{cid:?}");
                        if value.is_nan() {
                            assert!(parsed.value.is_nan(), "{cid:?}");
                        } else {
                            assert_eq!(parsed.value, value, "{cid:?}");
                            assert_eq!(
                                parsed.value.is_sign_negative(),
                                value.is_sign_negative(),
                                "{cid:?}"
                            );
                        }
                        assert_eq!(pos, field_end, "{cid:?}");
                    }
                    None => {
                        match result.unwrap_err() {
                            super::SelectionSyntaxError { pos: at, info, .. } => {
                                assert_eq!(at, 2, "{cid:?}");
                                assert_eq!(info.as_deref(), Some(" (expected number)"), "{cid:?}");
                            }
                            other => panic!("{cid:?}: {other:?}"),
                        }
                        assert_eq!(pos, 0, "{cid:?}");
                    }
                }
            }
        }
        assert_eq!(rows, 360, "every property/relation x tail row must execute");
    }

    #[test]
    fn bio_cid_c10_skip_preserves_separator() {
        // All 36 helper states of omit {0,1,2,3} x sep {Some0,Some2,None}
        // x semi {Some0,Some2,None}; the 30 satisfying the literal pinned
        // false guard `omit <= 1 && sep < semi` (None = npos = MAX) must
        // return the DEFAULT chain list plus the INCOMING sep unchanged —
        // select.cpp:165 `sep` is only assigned inside the guarded block.
        let seps = [Some(0usize), Some(2usize), None];
        let semis = [Some(0usize), Some(2usize), None];
        let mut guard_false_rows = 0usize;
        for omit in 0..=3 {
            for sep in seps {
                for semi in semis {
                    let guard_true =
                        omit <= 1 && sep.unwrap_or(usize::MAX) < semi.unwrap_or(usize::MAX);
                    if !guard_true {
                        guard_false_rows += 1;
                        let (list, out) = parse_cid_chain_stage("unused", omit, semi, sep).unwrap();
                        assert!(
                            list.all && list.list.is_empty(),
                            "omit={omit} sep={sep:?} semi={semi:?}"
                        );
                        assert_eq!(out, sep, "omit={omit} sep={sep:?} semi={semi:?}");
                    }
                }
            }
        }
        assert_eq!(guard_false_rows, 30, "guard-false state count");
    }

    #[test]
    fn bio_cid_c10_c09_handoff_separator_threading() {
        // Real C09 -> C10 handoffs: omission and semicolon positions are
        // derived from the actual input via the pinned C03 classifier and
        // a source ';' scan; the expected outgoing sep comes from the
        // pinned source, never from the C10 output.
        let handoffs: [(&str, Option<usize>); 5] = [
            ("1-5", Some(0)),
            ("(LYS)", Some(0)),
            ("CA[C]", Some(0)),
            ("/1;", Some(2)),
            ("/1", None),
        ];
        for (cid, expected) in handoffs {
            let omit = match determine_omitted_cid_fields(cid) {
                OmittedCidFields::Model => 0,
                OmittedCidFields::Chain => 1,
                OmittedCidFields::Residue => 2,
                OmittedCidFields::Atom => 3,
            };
            let semi = cid.as_bytes().iter().position(|byte| *byte == b';');
            let (model, model_sep) = parse_cid_model_stage(cid, omit, semi).unwrap();
            let (list, chain_sep) = parse_cid_chain_stage(cid, omit, semi, model_sep).unwrap();
            let _ = model;
            assert!(list.all, "default chain list for {cid:?}");
            assert_eq!(chain_sep, expected, "{cid:?}");
        }
    }

    // ---- BIO-CID C11: residue stage. Expectations from select.cpp:178-203
    // (re-check '(' after the seqid advance, dot/star compatibility
    // single consumption, point-copy default, terminal '/;\0' check) and
    // select.hpp:94-96 defaults (MIN/'*' and MAX/'*'). ----

    fn default_to() -> SelectionSequenceId {
        SelectionSequenceId {
            seqnum: i32::MAX,
            icode: b'*',
        }
    }

    #[test]
    fn bio_cid_c11_ranges_icode_point_copy_and_names() {
        // Plain ranges across negative/zero/positive bounds.
        let (from, to, names, sep) = parse_cid_residue_stage("1-5", 2, None, Some(0)).unwrap();
        assert_eq!((from.seqnum, from.icode), (1, b' '));
        assert_eq!((to.seqnum, to.icode), (5, b' '));
        assert!(names.all);
        assert_eq!(sep, Some(3));
        let (from, to, _, sep) = parse_cid_residue_stage("-3--1", 2, None, Some(0)).unwrap();
        assert_eq!((from.seqnum, to.seqnum), (-3, -1));
        assert_eq!(sep, Some(5));
        // Insertion codes: "14.a-20.b" keeps both letters.
        let (from, to, _, sep) = parse_cid_residue_stage("14.a-20.b", 2, None, Some(0)).unwrap();
        assert_eq!((from.seqnum, from.icode), (14, b'a'));
        assert_eq!((to.seqnum, to.icode), (20, b'b'));
        assert_eq!(sep, Some(9));
        // Point selection copies from_seqid when no '-' bound follows.
        let (from, to, _, sep) = parse_cid_residue_stage("14", 2, None, Some(0)).unwrap();
        assert_eq!(from.seqnum, 14);
        assert_eq!(to.seqnum, 14);
        assert_eq!(sep, Some(2));
        // Name lists: plain, inverted, multi and the wildcard default pair.
        let (from, to, names, sep) = parse_cid_residue_stage("(ALA)", 2, None, Some(0)).unwrap();
        assert_eq!((from.seqnum, to.seqnum), (i32::MIN, i32::MAX));
        assert_eq!(
            (names.all, names.list.as_slice()),
            (false, b"ALA".as_slice())
        );
        assert_eq!(sep, Some(5));
        let (_, _, names, _) = parse_cid_residue_stage("(ALA,GLY)/1", 2, None, Some(0)).unwrap();
        assert_eq!(names.list, b"ALA,GLY");
        let (_, _, names, _) = parse_cid_residue_stage("(!TRP)/1", 2, None, Some(0)).unwrap();
        assert!(names.inverted && names.list == b"TRP");
        // Re-check quirk: a seqid head followed directly by '(' enters the
        // name branch from the ADVANCED position.
        let (from, _, names, sep) = parse_cid_residue_stage("14(A)", 2, None, Some(0)).unwrap();
        assert_eq!(from.seqnum, 14);
        assert_eq!(names.list, b"A");
        assert_eq!(sep, Some(5));
        // Chain-stage handoff position (pos = chain sep + 1).
        let (from, _, _, sep) = parse_cid_residue_stage("/1/14-20", 0, None, Some(2)).unwrap();
        assert_eq!(from.seqnum, 14);
        assert_eq!(sep, Some(8));
    }

    #[test]
    fn bio_cid_c11_mmdb_compatibility_forms_and_skips() {
        // *(ALA). / *(ALA).* / (ALA). / (ALA).* / (ALA)* all parse the
        // same name list; each dot/star byte is consumed at most once.
        for (cid, expected_sep) in [
            ("*(ALA).", 7),
            ("*(ALA).*", 8),
            ("*(ALA)*", 7),
            ("(ALA).", 6),
            ("(ALA).*", 7),
            ("(ALA)*", 6),
        ] {
            let (_, _, names, sep) = parse_cid_residue_stage(cid, 2, None, Some(0)).unwrap();
            assert_eq!(
                (names.all, names.list.as_slice()),
                (false, b"ALA".as_slice()),
                "{cid}"
            );
            assert_eq!(sep, Some(expected_sep), "{cid}");
        }
        // '14.a' and '14a' are equivalent seqid forms.
        let (from, ..) = parse_cid_residue_stage("14a", 2, None, Some(0)).unwrap();
        assert_eq!((from.seqnum, from.icode), (14, b'a'));
        // Guard skip (atom-first omit 3): identity-on-skip keeps chain_sep.
        for cid in ["CA[C]", "q>1"] {
            let (from, to, names, sep) = parse_cid_residue_stage(cid, 3, None, Some(0)).unwrap();
            assert_eq!((from.seqnum, to.seqnum), (i32::MIN, i32::MAX), "{cid}");
            assert!(names.all, "{cid}");
            assert_eq!(sep, Some(0), "{cid}");
        }
        // chain_sep None (npos) skips the stage and returns None.
        let (_, _, _, sep) = parse_cid_residue_stage("/1", 0, None, None).unwrap();
        assert_eq!(sep, None);
    }

    #[test]
    fn bio_cid_c11_errors_wrap_quirk_and_range() {
        // Trailing junk after the bound: " (at residue)" at position 0.
        // A SINGLE trailing character is a legal icode (14x, 20x), so the
        // junk rows carry a second trailing character.
        for cid in ["14xy", "14-20xy", "(ALA)x"] {
            match parse_cid_residue_stage(cid, 2, None, Some(0)).unwrap_err() {
                BioSelectionParseError::Syntax(error) => {
                    assert_eq!(
                        (error.pos, error.info.as_deref()),
                        (0, Some(" (at residue)")),
                        "{cid}"
                    );
                }
                _ => panic!("{cid}: expected Syntax"),
            }
        }
        // Missing ')': right_br npos -> make_cid_list reads to the string
        // end, pos = right_br + 1 WRAPS to 0, and the terminal check on
        // cid[0] = '(' fires the same " (at residue)" note.
        let error = parse_cid_residue_stage("(ALA", 2, None, Some(0)).unwrap_err();
        match error {
            BioSelectionParseError::Syntax(error) => {
                assert_eq!(
                    (error.pos, error.info.as_deref()),
                    (0, Some(" (at residue)"))
                );
            }
            other => panic!("expected Syntax, got {other:?}"),
        }
        // Range errors propagate from the C06 owner at both bounds.
        for cid in ["5000000000-1", "1-5000000000"] {
            match parse_cid_residue_stage(cid, 2, None, Some(0)).unwrap_err() {
                BioSelectionParseError::SeqidRange(error) => {
                    assert_eq!(error.text, "5000000000", "{cid}");
                }
                other => panic!("{cid}: expected SeqidRange, got {other:?}"),
            }
        }
        // Name-list punctuation through the C04 default set (with '-').
        let error = parse_cid_residue_stage("(A[B)", 2, None, Some(0)).unwrap_err();
        match error {
            BioSelectionParseError::Syntax(error) => {
                assert_eq!(error.pos, 2);
                assert_eq!(error.info.as_deref(), Some(" ('[' in a list)"));
            }
            other => panic!("expected Syntax, got {other:?}"),
        }
        let _ = default_to();
    }

    // ---- BIO-CID C12: atom stage. Expectations from select.cpp:204-227:
    // guard sep < min(size, semi); names to the first of "[:;"; the
    // empty-to-wildcard rule only when NOT inverted; end == npos is the
    // source EARLY RETURN; '[..]' via the C07 owner; ':' altlocs to the
    // semicolon. ----

    #[test]
    fn bio_cid_c12_atom_names_elements_and_altlocs() {
        // Plain name list to the string end: source EARLY RETURN fires.
        let (names, elements, altlocs, early) = parse_cid_atom_stage("CA", None, Some(0)).unwrap();
        assert_eq!(
            (names.all, names.list.as_slice()),
            (false, b"CA".as_slice())
        );
        assert!(elements.is_none() && altlocs.all && early);
        // Bracket terminator before the semicolon.
        let (names, elements, _, early) =
            parse_cid_atom_stage("CA[C];q>1", Some(5), Some(0)).unwrap();
        assert_eq!(names.list, b"CA");
        assert!(elements.is_some() && !early);
        // Element brackets and altlocs combine after names.
        let (_, elements, altlocs, early) =
            parse_cid_atom_stage("CA[C]:A,B", None, Some(0)).unwrap();
        assert!(elements.is_some());
        assert_eq!(altlocs.list, b"A,B");
        assert!(!early);
        // Altlocs with no names (end == pos) and no elements.
        let (names, elements, altlocs, early) = parse_cid_atom_stage(":A", None, Some(0)).unwrap();
        assert!(names.all && elements.is_none());
        assert_eq!(altlocs.list, b"A");
        assert!(!early);
        // Inverted altloc list keeps its inversion.
        let (_, _, altlocs, _) = parse_cid_atom_stage(":!A", None, Some(0)).unwrap();
        assert!(altlocs.inverted && altlocs.list == b"A");
        // Empty altloc member list stays empty (the wildcard rewrite
        // applies to ATOM NAMES only).
        let (_, _, altlocs, _) = parse_cid_atom_stage(":", None, Some(0)).unwrap();
        assert_eq!(
            (altlocs.all, altlocs.inverted, altlocs.list.as_slice()),
            (false, false, b"".as_slice())
        );
        // After the residue stage separator (pos = sep + 1).
        let (names, _, _, _) = parse_cid_atom_stage("/1/14/CA", None, Some(5)).unwrap();
        assert_eq!(names.list, b"CA");
    }

    #[test]
    fn bio_cid_c12_empty_name_rules_and_early_return() {
        // Empty INVERTED atom name list stays inverted-empty: the
        // empty-to-wildcard rule fires only when NOT inverted (source
        // comment: chain name can be empty, atom name cannot). A
        // NON-inverted empty list is unreachable through make_cid_list —
        // emptiness only arises after '!' (inverted) or '*' (already
        // all) flag consumption — so the rule is exercised via the
        // inverted negative here and the guard-skip below.
        let (names, elements, _, early) = parse_cid_atom_stage("!;q>1", Some(1), Some(0)).unwrap();
        assert_eq!(
            (names.all, names.inverted, names.list.as_slice()),
            (false, true, b"".as_slice())
        );
        assert!(elements.is_none() && !early);
        // sep == semi (empty atom field): the guard skips the stage.
        let (names, elements, altlocs, early) =
            parse_cid_atom_stage(";q>1", Some(0), Some(0)).unwrap();
        assert_eq!(
            (names.all, names.inverted, names.list.as_slice()),
            (true, false, b"".as_slice())
        );
        assert!(elements.is_none() && altlocs.all && !early);
        // Wildcard name to the string end: early return, no elements.
        let (names, elements, _, early) = parse_cid_atom_stage("*", None, Some(0)).unwrap();
        assert!(names.all && names.list.is_empty());
        assert!(elements.is_none() && early);
    }

    #[test]
    fn bio_cid_c12_guards_and_missing_bracket() {
        // Guard skip: residue sep at the string end (sep == size).
        let (names, elements, altlocs, _) = parse_cid_atom_stage("/1", None, Some(2)).unwrap();
        assert!(names.all && elements.is_none() && altlocs.all);
        // Guard skip: sep == semi.
        let (names, _, _, _) = parse_cid_atom_stage("A/1;", Some(3), Some(3)).unwrap();
        assert!(names.all);
        // Missing ']': fixed note at position 0.
        let error = parse_cid_atom_stage("CA[C", None, Some(0)).unwrap_err();
        match error {
            BioSelectionParseError::Syntax(error) => {
                assert_eq!(
                    (error.pos, error.info.as_deref()),
                    (0, Some(" (no matching ']')"))
                );
            }
            other => panic!("expected Syntax, got {other:?}"),
        }
        // Altloc list punctuation still runs through the C04 default set.
        let error = parse_cid_atom_stage("CA:C[X", None, Some(0)).unwrap_err();
        match error {
            BioSelectionParseError::Syntax(error) => {
                assert_eq!(error.pos, 4);
                assert_eq!(error.info.as_deref(), Some(" ('[' in a list)"));
            }
            other => panic!("expected Syntax, got {other:?}"),
        }
    }

    // ---- BIO-CID C14: composed decode via read_bio_selection_cid.
    // Expectations derived independently from select.cpp:131-257 stage
    // order and the select.hpp member defaults. ----

    #[test]
    fn bio_cid_c14_omission_classes_and_stage_presence() {
        // Model-first: every stage runs in order.
        let data = read_bio_selection_cid("/2/A,B/14-20/CA[C]:A").unwrap();
        let p = data.parts();
        assert_eq!(p.mdl, 2);
        assert_eq!(p.chain_ids.list, "A,B");
        assert_eq!((p.from_seqid.seqnum, p.to_seqid.seqnum), (14, 20));
        assert_eq!(p.atom_names.list, "CA");
        assert!(p.elements.is_some());
        assert_eq!(p.altlocs.list, "A");
        // Chain-first: model skipped, sep 0 threads into the chain stage.
        let data = read_bio_selection_cid("A/14/CA").unwrap();
        let p = data.parts();
        assert_eq!(p.mdl, 0);
        assert_eq!(p.chain_ids.list, "A");
        assert_eq!(p.from_seqid.seqnum, 14);
        assert_eq!(p.atom_names.list, "CA");
        // Residue-first: model and chain skipped (sep stays 0).
        let data = read_bio_selection_cid("14-20/CA").unwrap();
        let p = data.parts();
        assert_eq!(p.mdl, 0);
        assert!(p.chain_ids.all);
        assert_eq!((p.from_seqid.seqnum, p.to_seqid.seqnum), (14, 20));
        assert_eq!(p.atom_names.list, "CA");
        // Atom-first: only the atom stage runs.
        let data = read_bio_selection_cid("CA[C]").unwrap();
        let p = data.parts();
        assert!(p.chain_ids.all && p.residue_names.all);
        assert_eq!(p.atom_names.list, "CA");
        assert!(p.elements.is_some());
        // Bypass: empty and lone '*' decode to the all-wildcard default.
        for cid in ["", "*"] {
            let data = read_bio_selection_cid(cid).unwrap();
            let p = data.parts();
            assert_eq!(p.mdl, 0, "{cid}");
            assert!(p.chain_ids.all && p.atom_names.all, "{cid}");
        }
    }

    #[test]
    fn bio_cid_c14_extension_kinds_and_early_return() {
        // Entity extension: et_flags filled and set at Polymer(1)/Water(4).
        let data = read_bio_selection_cid("CA;polymer").unwrap();
        let p = data.parts();
        assert_eq!(p.entity_types.list, "polymer");
        assert_eq!(p.et_flags, &[false, true, false, false, false, false]);
        // Inverted entity: fill true, matched item false.
        let data = read_bio_selection_cid("CA;!solvent").unwrap();
        let p = data.parts();
        assert_eq!(p.et_flags, &[true, true, true, true, false, true]);
        // Inequality extension: atom_inequalities carried through.
        let data = read_bio_selection_cid("CA;q>1.5").unwrap();
        let p = data.parts();
        assert_eq!(p.atom_inequalities.len(), 1);
        assert_eq!(p.atom_inequalities[0].value, 1.5);
        // Atom names to the string end: source early return SKIPS the
        // extension loop even though a semicolon region exists in the
        // buffer sense — here "CA" has none, so assert no extensions.
        let data = read_bio_selection_cid("A/1/CA").unwrap();
        let p = data.parts();
        assert!(p.entity_types.all);
        assert!(p.atom_inequalities.is_empty());
        // Names to end with a trailing extension string in the same input
        // would need a ';' — but the name list to end consumes everything,
        // so an early return with a semicolon cannot occur in one cid.
    }

    #[test]
    fn bio_cid_c14_malformed_positions_compose() {
        // Each stage's fixed note surfaces unchanged through the entry.
        let cases = [
            ("/1x/A", " (at model number)"),
            ("A(B/1", " (at residue)"),
            ("14-20xy", " (at residue)"),
            ("CA[C", " (no matching ']')"),
            (";ligand", " at ligand"),
        ];
        for (cid, info) in cases {
            match read_bio_selection_cid(cid) {
                Err(BioSelectionParseError::Syntax(error)) => {
                    assert_eq!(error.info.as_deref(), Some(info), "{cid}");
                    assert!(error.cid == cid, "{cid}");
                }
                _ => panic!("{cid}: expected Syntax"),
            }
        }
        // Typed range errors from both numeric owners.
        match read_bio_selection_cid("5000000000-1") {
            Err(BioSelectionParseError::SeqidRange(error)) => {
                assert_eq!(error.text, "5000000000");
            }
            _ => panic!("expected SeqidRange"),
        }
    }

    #[test]
    fn bio_cid_c13_utf8_members_no_substitution() {
        // Multibyte unknown entity text appears EXACTLY in the source
        // message bytes — never an empty substitute from UTF-8 handling.
        let error = match parse_cid_extensions(";Ünknown", Some(0)) {
            Err(error) => error,
            Ok(stage) => panic!(
                "expected error, got {} inequalities",
                stage.inequalities.len()
            ),
        };
        match error {
            BioSelectionParseError::Syntax(error) => {
                assert_eq!(error.info.as_deref(), Some(" at Ünknown"));
                // Full wrong_syntax byte composition (pos 0: no window):
                // header + info + ": " + cid, byte-exact incl. multibyte.
                assert_eq!(
                    error.source_message_bytes(),
                    b"Invalid selection syntax at \xc3\x9cnknown: ;\xc3\x9cnknown".as_slice()
                );
            }
            other => panic!("expected Syntax, got {other:?}"),
        }
        // Empty comma members remain errors (split keeps the empty
        // remainder): ";polymer," -> items ["polymer", ""] -> " at ".
        let error = match parse_cid_extensions(";polymer,", Some(0)) {
            Err(error) => error,
            Ok(stage) => panic!("expected error, got {}", stage.polymer),
        };
        match error {
            BioSelectionParseError::Syntax(error) => {
                assert_eq!(error.pos, 0);
                assert_eq!(error.info.as_deref(), Some(" at "));
            }
            other => panic!("expected Syntax, got {other:?}"),
        }
    }

    #[test]
    fn bio_cid_c14_complete_decode_matrix() {
        // Six literal extension actions from the pinned source
        // (select.cpp:229-257), INDEPENDENT of any parse output:
        // entity lists replace state and re-fill et_flags to inv then
        // flip the matched item bit; inequalities append in order.
        struct Action {
            tag: &'static str,
            // (list text, inverted)
            entity: Option<(&'static str, bool)>,
            flags: [bool; 6],
            ineq: Option<(u8, i32, f64)>,
        }
        let actions = [
            Action {
                tag: "polymer",
                entity: Some(("polymer", false)),
                flags: [false, true, false, false, false, false],
                ineq: None,
            },
            Action {
                tag: "solvent",
                entity: Some(("solvent", false)),
                flags: [false, false, false, false, true, false],
                ineq: None,
            },
            Action {
                tag: "!polymer",
                entity: Some(("polymer", true)),
                flags: [true, false, true, true, true, true],
                ineq: None,
            },
            Action {
                tag: "!solvent",
                entity: Some(("solvent", true)),
                flags: [true, true, true, true, false, true],
                ineq: None,
            },
            Action {
                tag: "q>1.5",
                entity: None,
                flags: [false; 6],
                ineq: Some((b'q', 1, 1.5)),
            },
            Action {
                tag: "b<2",
                entity: None,
                flags: [false; 6],
                ineq: Some((b'b', -1, 2.0)),
            },
        ];
        // Frozen per-shape literal bases (derived from select.cpp stage
        // order + select.hpp member defaults, not from decoder output).
        struct Shape {
            cid: &'static str,
            mdl: i32,
            chain: Option<&'static str>,
            from: (i32, u8),
            to: (i32, u8),
            atoms: Option<&'static str>,
            elements: bool,
            altlocs: Option<&'static str>,
        }
        let shapes = [
            Shape {
                cid: "/2",
                mdl: 2,
                chain: None,
                from: (i32::MIN, b'*'),
                to: (i32::MAX, b'*'),
                atoms: None,
                elements: false,
                altlocs: None,
            },
            Shape {
                cid: "/2/A",
                mdl: 2,
                chain: Some("A"),
                from: (i32::MIN, b'*'),
                to: (i32::MAX, b'*'),
                atoms: None,
                elements: false,
                altlocs: None,
            },
            Shape {
                cid: "/2/A/14-20",
                mdl: 2,
                chain: Some("A"),
                from: (14, b' '),
                to: (20, b' '),
                atoms: None,
                elements: false,
                altlocs: None,
            },
            Shape {
                cid: "/2/A/14-20/CA[C]:B",
                mdl: 2,
                chain: Some("A"),
                from: (14, b' '),
                to: (20, b' '),
                atoms: Some("CA"),
                elements: true,
                altlocs: Some("B"),
            },
            Shape {
                cid: "A",
                mdl: 0,
                chain: Some("A"),
                from: (i32::MIN, b'*'),
                to: (i32::MAX, b'*'),
                atoms: None,
                elements: false,
                altlocs: None,
            },
            Shape {
                cid: "A/14-20",
                mdl: 0,
                chain: Some("A"),
                from: (14, b' '),
                to: (20, b' '),
                atoms: None,
                elements: false,
                altlocs: None,
            },
            Shape {
                cid: "A/14-20/CA[C]:B",
                mdl: 0,
                chain: Some("A"),
                from: (14, b' '),
                to: (20, b' '),
                atoms: Some("CA"),
                elements: true,
                altlocs: Some("B"),
            },
            Shape {
                cid: "14-20",
                mdl: 0,
                chain: None,
                from: (14, b' '),
                to: (20, b' '),
                atoms: None,
                elements: false,
                altlocs: None,
            },
            Shape {
                cid: "14-20/CA[C]:B",
                mdl: 0,
                chain: None,
                from: (14, b' '),
                to: (20, b' '),
                atoms: Some("CA"),
                elements: true,
                altlocs: Some("B"),
            },
            Shape {
                cid: "CA[C]:B",
                mdl: 0,
                chain: None,
                from: (i32::MIN, b'*'),
                to: (i32::MAX, b'*'),
                atoms: Some("CA"),
                elements: true,
                altlocs: Some("B"),
            },
        ];
        let mut rows = 0usize;
        for shape in &shapes {
            // 43 configs: no extension, each single, every ordered pair.
            let mut configs: Vec<Vec<usize>> = vec![Vec::new()];
            for first in 0..6 {
                configs.push(vec![first]);
                for second in 0..6 {
                    configs.push(vec![first, second]);
                }
            }
            assert_eq!(configs.len(), 43);
            for config in &configs {
                let mut cid = String::from(shape.cid);
                for &index in config {
                    cid.push(';');
                    cid.push_str(actions[index].tag);
                }
                // Compose the EXPECTED extension state purely from the
                // literal action table: last entity action wins for
                // entity_types and flags; inequalities accumulate in order.
                let mut entity = ("", false);
                let mut has_entity = false;
                let mut flags = [false; 6];
                let mut ineqs: Vec<(u8, i32, f64)> = Vec::new();
                for &index in config {
                    let action = &actions[index];
                    if let Some((list, inverted)) = action.entity {
                        has_entity = true;
                        entity = (list, inverted);
                        flags = action.flags;
                    }
                    if let Some(literal) = action.ineq {
                        ineqs.push(literal);
                    }
                }
                let data = read_bio_selection_cid(&cid).unwrap();
                let parts = data.parts();
                assert_eq!(parts.mdl, shape.mdl, "{cid}");
                match shape.chain {
                    Some(list) => {
                        assert!(!parts.chain_ids.all, "{cid}");
                        assert_eq!(parts.chain_ids.list, list, "{cid}");
                    }
                    None => {
                        assert!(parts.chain_ids.all, "{cid}");
                        assert_eq!(parts.chain_ids.list, "", "{cid}");
                    }
                }
                assert_eq!(parts.from_seqid.seqnum, shape.from.0, "{cid}");
                assert_eq!(parts.from_seqid.icode, shape.from.1, "{cid}");
                assert_eq!(parts.to_seqid.seqnum, shape.to.0, "{cid}");
                assert_eq!(parts.to_seqid.icode, shape.to.1, "{cid}");
                assert!(parts.residue_names.all, "{cid}");
                assert_eq!(parts.residue_names.list, "", "{cid}");
                match shape.atoms {
                    Some(list) => {
                        assert!(!parts.atom_names.all, "{cid}");
                        assert_eq!(parts.atom_names.list, list, "{cid}");
                    }
                    None => assert!(parts.atom_names.all, "{cid}"),
                }
                match shape.elements {
                    true => {
                        let mask = parts.elements.expect("carbon mask");
                        // Carbon ordinal 6 on the Gemmi element table;
                        // every other element bit stays clear.
                        assert!(mask[6], "{cid}");
                        assert_eq!(mask.iter().filter(|bit| **bit).count(), 1, "{cid}");
                    }
                    false => assert!(parts.elements.is_none(), "{cid}"),
                }
                match shape.altlocs {
                    Some(list) => {
                        assert!(!parts.altlocs.all, "{cid}");
                        assert_eq!(parts.altlocs.list, list, "{cid}");
                    }
                    None => assert!(parts.altlocs.all, "{cid}"),
                }
                if has_entity {
                    assert!(!parts.entity_types.all, "{cid}");
                    assert_eq!(parts.entity_types.list, entity.0, "{cid}");
                    assert_eq!(parts.entity_types.inverted, entity.1, "{cid}");
                } else {
                    assert!(parts.entity_types.all, "{cid}");
                    assert_eq!(parts.entity_types.list, "", "{cid}");
                }
                assert_eq!(parts.et_flags, &flags, "{cid}");
                assert_eq!(parts.atom_inequalities.len(), ineqs.len(), "{cid}");
                for (index, (property, relation, value)) in ineqs.iter().enumerate() {
                    assert_eq!(parts.atom_inequalities[index].property, *property, "{cid}");
                    assert_eq!(parts.atom_inequalities[index].relation, *relation, "{cid}");
                    assert_eq!(parts.atom_inequalities[index].value, *value, "{cid}");
                }
                rows += 1;
            }
        }
        assert_eq!(rows, 430);
    }

    // ---- BIO-CID C15: List::str and SequenceId::str serialization.
    // Expectations derived independently from select.hpp:27-31 List::str
    // and select.cpp:267-277 SequenceId::str (empty() = select.hpp:56-58). ----

    fn c15_list(all: bool, inverted: bool, list: &str) -> SelectionList {
        SelectionList {
            all,
            inverted,
            list: list.to_string(),
        }
    }

    #[test]
    fn bio_cid_c15_list_str_all_states_and_members() {
        // Wildcard wins regardless of members/inversion.
        assert_eq!(cid_list_str(&c15_list(true, false, "")), "*");
        assert_eq!(cid_list_str(&c15_list(true, true, "A")), "*");
        // Plain non-wildcard: exact member text, empty included.
        assert_eq!(cid_list_str(&c15_list(false, false, "A,B")), "A,B");
        assert_eq!(cid_list_str(&c15_list(false, false, "")), "");
        // Inverted non-wildcard: "!" + exact member text (empty too).
        assert_eq!(cid_list_str(&c15_list(false, true, "A,B")), "!A,B");
        assert_eq!(cid_list_str(&c15_list(false, true, "")), "!");
        // Empty members and multibyte names survive byte-exact.
        assert_eq!(cid_list_str(&c15_list(false, false, "a,,b")), "a,,b");
        assert_eq!(
            cid_list_str(&c15_list(false, true, "\u{3a9},A")),
            "!\u{3a9},A"
        );
    }

    #[test]
    fn bio_cid_c15_list_str_complete_state_matrix() {
        // 16 rows: all x inverted x four member texts. Literal outputs
        // fixed independently from select.hpp:27-31; the fixture member
        // text must remain UNCHANGED after the call (borrowed, not
        // consumed).
        let texts = ["", "A,B", "a,,b", "\u{3a9},A"];
        let mut rows = 0;
        for text in texts {
            for all in [false, true] {
                for inverted in [false, true] {
                    let list = c15_list(all, inverted, text);
                    let expected = if all {
                        "*".to_string()
                    } else if inverted {
                        format!("!{text}")
                    } else {
                        text.to_string()
                    };
                    assert_eq!(cid_list_str(&list), expected, "{all}/{inverted}/{text}");
                    assert_eq!(list.list, text, "members unchanged");
                    rows += 1;
                }
            }
        }
        assert_eq!(rows, 16);
    }

    #[test]
    fn bio_cid_c15_seqid_str_numbers_and_icodes() {
        // Both sentinels: source emptiness (empty bytes, no dot).
        assert_eq!(
            seqid_str(&SelectionSequenceId {
                seqnum: i32::MIN,
                icode: b'A'
            })
            .as_slice(),
            b""
        );
        assert_eq!(
            seqid_str(&SelectionSequenceId {
                seqnum: i32::MAX,
                icode: b' '
            })
            .as_slice(),
            b""
        );
        // negative/zero/positive x blank/star/letter insertion codes,
        // asserted as EXACT bytes (Vec<u8> output, raw suffix).
        let cases: [(i32, u8, &[u8]); 9] = [
            (-7, b' ', b"-7."),
            (-7, b'*', b"-7"),
            (-7, b'A', b"-7.A"),
            (0, b' ', b"0."),
            (0, b'*', b"0"),
            (0, b'z', b"0.z"),
            (5, b' ', b"5."),
            (5, b'*', b"5"),
            (5, b'K', b"5.K"),
        ];
        for (seqnum, icode, expected) in cases {
            assert_eq!(seqid_str(&SelectionSequenceId { seqnum, icode }), expected);
        }
    }

    #[test]
    fn bio_cid_c15_seqid_str_all_icode_bytes() {
        // 1536 rows: six literal numbers x every u8 icode. Expected
        // bytes built ONLY from independently fixed literals and the
        // source suffix rule — no production formatter, no lossy path.
        let numbers: [(i32, &[u8]); 6] = [
            (i32::MIN, b""), // sentinel: empty regardless of icode
            (i32::MAX, b""), // sentinel: empty regardless of icode
            (i32::MIN + 1, b"-2147483647"),
            (-7, b"-7"),
            (0, b"0"),
            (i32::MAX - 1, b"2147483646"),
        ];
        let mut rows = 0;
        for (seqnum, prefix) in numbers {
            for icode in 0u8..=255 {
                let out = seqid_str(&SelectionSequenceId { seqnum, icode });
                let mut expected = prefix.to_vec();
                if seqnum != i32::MIN && seqnum != i32::MAX {
                    if icode != b'*' {
                        expected.push(b'.');
                        if icode != b' ' {
                            expected.push(icode);
                        }
                    }
                }
                assert_eq!(out, expected, "seqnum {seqnum} icode {icode:#04x}");
                rows += 1;
            }
        }
        // Explicit raw high-byte assertions (no Unicode substitution).
        assert_eq!(
            seqid_str(&SelectionSequenceId {
                seqnum: 0,
                icode: 0x80
            })
            .as_slice(),
            b"0.\x80"
        );
        assert_eq!(
            seqid_str(&SelectionSequenceId {
                seqnum: -7,
                icode: 0xff
            })
            .as_slice(),
            b"-7.\xff"
        );
        assert_eq!(rows, 1536);
    }

    // ---- BIO-CID C22: Selection::str over the borrowed predicate view,
    // composing C15/C21 serializers. Expectations derived independently
    // from select.cpp:288-327; element-mask expectations compose the
    // pinned element-name table with test-side selection rules mirroring
    // the source loop (strict count > 64 inversion, cid.back() = ']'
    // including degenerate masks). ----

    fn c22_cid_list(all: bool, inverted: bool, text: &str) -> SelectionCidList {
        SelectionCidList {
            all,
            inverted,
            list: text.as_bytes().to_vec(),
        }
    }

    fn c22_cid_list_all() -> SelectionCidList {
        // Source List default: `bool all = true;` (select.hpp:23).
        SelectionCidList {
            all: true,
            inverted: false,
            list: Vec::new(),
        }
    }

    fn c22_data(
        mask: Option<[bool; 120]>,
        inequalities: Vec<SelectionAtomInequality>,
    ) -> BioSelectionData {
        BioSelectionData::from_parts(BioSelectionParts {
            mdl: 0,
            chain_ids: c22_cid_list_all(),
            from_seqid: SelectionSequenceId {
                seqnum: i32::MIN,
                icode: b'*',
            },
            to_seqid: SelectionSequenceId {
                seqnum: i32::MAX,
                icode: b'*',
            },
            residue_names: c22_cid_list_all(),
            entity_types: c22_cid_list_all(),
            et_flags: [false; 6],
            atom_names: c22_cid_list_all(),
            elements: mask,
            altlocs: c22_cid_list_all(),
            atom_inequalities: inequalities,
        })
        .unwrap()
    }

    fn to_cid_text(cid: &str) -> String {
        bio_selection_to_cid(&read_bio_selection_cid(cid).unwrap())
    }

    #[test]
    fn bio_cid_c22_decode_serialize_rows_and_nonroundtrips() {
        // Field order and defaults per select.cpp:288-327 with the
        // source List default `all = true` (select.hpp:23): chain_ids.str
        // emits "*" for the default wildcard chain list, from_seqid is
        // EMPTY (sentinel seqnum), and atom_names.all SUPPRESSES the
        // atom-name output entirely (select.cpp:310-311 `if
        // (!atom_names.all)`), so a fully-default selection prints
        // "//*//". Single seqids copy from->to (select.cpp:200-202) so no
        // range is emitted; "14A" re-serializes as "14.A" (icode dot) and
        // a >9-significant-digit value serializes through the rounded
        // %.9g owner; parse_cid_seqid defaults icode to blank (select.cpp:96)
        // so bare seqids re-serialize as "14."/"20." — all explicit
        // NON-roundtrips; repeated entity
        // extensions keep the last one; inequality order is preserved.
        let rows: [(&str, &str); 10] = [
            ("/", "//*//"),
            ("/2", "/2/*//"),
            ("/2/A", "/2/A//"),
            ("A", "//A//"),
            ("/2/A/14-20", "/2/A/14.-20./"),
            ("/2/A/14-20/CA[C]:B", "/2/A/14.-20./CA[C]:B"),
            ("/2/A/14A", "/2/A/14.A/"),
            (
                "/2/A/14-20/CA[C]:B;q>0.1000000005",
                "/2/A/14.-20./CA[C]:B;q>0.100000001",
            ),
            (
                "/2/A/14-20/CA[C]:B;polymer;!polymer",
                "/2/A/14.-20./CA[C]:B;!polymer",
            ),
            (
                "/2/A/14-20;q>1.5;b<2.5;q>1.5",
                "/2/A/14.-20./;q>1.5;b<2.5;q>1.5",
            ),
        ];
        for (input, expected) in rows {
            assert_eq!(to_cid_text(input), expected, "{input}");
        }
    }

    #[test]
    fn bio_cid_c22_element_threshold_and_degenerate_masks() {
        // Strict > 64 threshold: 64 selected elements list directly; 65
        // selected elements invert and list the 55-name complement.
        // Degenerate masks reproduce cid.back() = ']': all-false (no '!',
        // nothing appended) leaves a bare ']'; all-true pushes '!' which
        // the back-replacement turns into ']' leaving "[]".
        // Expected strings use DISPLAY — an INDEPENDENT literal display
        // table (elem.hpp:269-294 "He"/"Li"/"Cl" spellings), never a
        // production table (neither the uppercase parser vocabulary nor
        // gemmi_element_name).
        const DISPLAY: [&str; 121] = [
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
        let mut mask_64 = [false; 120];
        for bit in mask_64.iter_mut().take(64) {
            *bit = true;
        }
        let mut mask_65 = mask_64;
        mask_65[64] = true;
        let mut mask_carbon = [false; 120];
        mask_carbon[6] = true; // Gemmi El::C ordinal (DISPLAY[6]).

        let expected_64 = format!(
            "[{}]",
            (0..64).map(|i| DISPLAY[i]).collect::<Vec<_>>().join(",")
        );
        let expected_65 = format!(
            "[!{}]",
            (65..120).map(|i| DISPLAY[i]).collect::<Vec<_>>().join(",")
        );
        assert_eq!(
            bio_selection_to_cid(&c22_data(Some(mask_64), vec![])),
            format!("//*//{expected_64}")
        );
        assert_eq!(
            bio_selection_to_cid(&c22_data(Some(mask_65), vec![])),
            format!("//*//{expected_65}")
        );
        assert_eq!(
            bio_selection_to_cid(&c22_data(Some([false; 120]), vec![])),
            "//*//]"
        );
        assert_eq!(
            bio_selection_to_cid(&c22_data(Some([true; 120]), vec![])),
            "//*//[]"
        );
        assert_eq!(
            bio_selection_to_cid(&c22_data(Some(mask_carbon), vec![])),
            "//*//[C]"
        );
        // Literal anchors for the direct-listing shape, including the
        // DISPLAY spellings the uppercase table got wrong (He/Li/Cl).
        let mut mask_head = [false; 120];
        for index in 0..4 {
            mask_head[index] = true;
        }
        assert_eq!(
            bio_selection_to_cid(&c22_data(Some(mask_head), vec![])),
            "//*//[X,H,He,Li]"
        );
        let mut mask_chlorine = [false; 120];
        mask_chlorine[17] = true;
        assert_eq!(
            bio_selection_to_cid(&c22_data(Some(mask_chlorine), vec![])),
            "//*//[Cl]"
        );
    }

    #[test]
    fn bio_cid_c22_all_singleton_masks_and_multibyte_lists() {
        // Every one of the 120 modeled ordinals as a singleton mask,
        // expected via an INDEPENDENT literal display table: each mask
        // serializes as "[<display name>]" exactly (He/Li/Cl included).
        // Decoder-backed multibyte list preservation at the end: decoded
        // CIDs with multibyte UTF-8 chain/residue members re-serialize
        // the exact bytes (decoded-CID UTF-8 provenance, not from_parts
        // raw-byte generality).
        const DISPLAY: [&str; 121] = [
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
        for ordinal in 0..120_usize {
            let mut mask = [false; 120];
            mask[ordinal] = true;
            let expected = format!("//*//[{}]", DISPLAY[ordinal]);
            assert_eq!(
                bio_selection_to_cid(&c22_data(Some(mask), vec![])),
                expected,
                "ordinal {ordinal}"
            );
        }
        assert_eq!(to_cid_text("//\u{3a9},A"), "//\u{3a9},A//");
        assert_eq!(to_cid_text("//A/(\u{3a9})"), "//A/(\u{3a9})/");
    }

    #[test]
    fn bio_cid_c22_repeated_inequality_order_via_parts() {
        // atom_inequalities serialize in stored order, repetitions kept
        // (select.cpp:324-325), via the canonical %.9g owner.
        let inequalities = vec![
            SelectionAtomInequality {
                property: b'q',
                relation: 1,
                value: 1.5,
            },
            SelectionAtomInequality {
                property: b'b',
                relation: -1,
                value: 2.5,
            },
            SelectionAtomInequality {
                property: b'q',
                relation: 1,
                value: 1.5,
            },
        ];
        assert_eq!(
            bio_selection_to_cid(&c22_data(None, inequalities)),
            "//*//;q>1.5;b<2.5;q>1.5"
        );
    }

    // ---- BIO-CID C21: AtomInequality::str through the canonical %.9g
    // owner. Expectations derived independently from select.cpp:280-285
    // (';' + raw property byte + relation char + to_str %.9g value). ----

    // ---- BIO-CID C25: real-CID parse -> BioSelectionData -> C24 detached
    // query over REAL reader-populated tables (PDB/mmCIF fixtures + an
    // mmJSON inline fixture), closing the C14 Step82 deferred obligation.
    // Expectations derive independently from the Gemmi CID closure and the
    // fixture texts (see the C25 matrix in selection_receipt.md). ----

    fn c25_query(structure: &cosmolkit_bio::BioStructureData, cid: &str) -> Vec<u32> {
        let selection = read_bio_selection_cid(cid).unwrap();
        selection
            .selected_bio_atom_ids(
                structure.models(),
                structure.chains(),
                structure.residues(),
                structure.atoms(),
                structure.input_format(),
            )
            .unwrap()
            .iter()
            .map(|id| id.value())
            .collect()
    }

    fn c25_structure(
        text: &str,
        format: cosmolkit_bio::BioCoordinateFormat,
    ) -> cosmolkit_bio::BioStructureData {
        crate::bio_read::read_bio_structure(
            text,
            &crate::bio_read::BioReadParams {
                format,
                source_name: "c25".to_owned(),
                ..crate::bio_read::BioReadParams::default()
            },
        )
        .unwrap()
    }

    #[test]
    fn bio_cid_c25_pdb_fixture_end_to_end_and_nonroundtrips() {
        let text = include_str!("../../../testdata/bio/fixtures/gemmi_full_feature_sample.pdb");
        let structure = c25_structure(text, cosmolkit_bio::BioCoordinateFormat::Pdb);
        // Source order: 0=SG CYS 3, 1=SG CYS 4, 2=CA ALA 7, 3=O HOH 2,
        // 4=ZN1 ZN 1 altloc A (single model 1, chain A, occ 1.00, B 20).
        let rows: [(&str, &[u32]); 25] = [
            ("A/3-4", &[0, 1]),
            ("A/7", &[2]),
            ("A/3A", &[]), // no insertion code A row
            ("A/(ALA)", &[2]),
            ("A/(HOH)", &[3]),
            ("A/(CYS,ALA)", &[0, 1, 2]),
            ("A/(ALA).", &[2]), // MMDB compat dot
            ("A/(ALA)*", &[2]), // MMDB compat star
            ("A/*", &[0, 1, 2, 3, 4]),
            ("A", &[0, 1, 2, 3, 4]),
            ("/1/A", &[0, 1, 2, 3, 4]),
            ("/2/A", &[]), // absent model number
            ("B", &[]),
            ("A//CA[C]", &[2]),
            ("A//SG[S]", &[0, 1]),
            ("A//O[O]", &[3]),
            ("A//O[!C]", &[3]), // mask complement excludes carbon atoms
            ("A//ZZ[D]", &[]),  // no deuterium rows
            ("A//O[!H]", &[3]), // inverted mask: O is not hydrogen
            ("A//:A", &[]),     // no altloc rows: ZN1A is a 4-char NAME, alt blank
            ("A//:B", &[]),
            (";solvent", &[]), // PDB fixture leaves HOH entity kind untyped
            (";!polymer", &[0, 1, 2, 3, 4]), // inverted mask leaves Unknown set
            ("A/;q>0.5", &[0, 1, 2, 3, 4]),
            ("A/;b<10", &[]),
        ];
        for (cid, expected) in rows {
            assert_eq!(c25_query(&structure, cid), expected, "{cid}");
        }
        // Additional stage rows: two ordered inequalities + B factor.
        assert_eq!(c25_query(&structure, "A/3-4;q<2;b>10"), vec![0, 1]);
        assert_eq!(c25_query(&structure, "A/;b<25"), vec![0, 1, 2, 3, 4]);
        // Source non-roundtrips through the real decode->serialize path.
        assert_eq!(to_cid_text("A/3-4"), "//A/3.-4./");
        assert_eq!(to_cid_text("A/7"), "//A/7./");
        assert_eq!(to_cid_text("A/(ALA)."), "//A/(ALA)/");
        assert_eq!(to_cid_text("A"), "//A//");
        assert_eq!(to_cid_text(";solvent"), "//*//;solvent");
        assert_eq!(to_cid_text("A//:B"), "//A//:B");
    }

    #[test]
    fn bio_cid_c25_mmcif_and_mmjson_fixture_rows() {
        let cif_text = include_str!("../../../testdata/bio/fixtures/gemmi_full_feature_sample.cif");
        let cif = c25_structure(cif_text, cosmolkit_bio::BioCoordinateFormat::Mmcif);
        // mmCIF fixture: chain A (polymer entity 1), CYS auth_seq 101/102,
        // atoms SG,SG in source order; reader stores the AUTH seqids
        // (auth-precedence), so label-seqid CIDs match nothing.
        let cif_rows: [(&str, &[u32]); 10] = [
            ("X/101-102", &[0, 1]),
            ("X/1-2", &[]),
            ("X/(CYS)", &[0, 1]),
            ("X//SG[S]", &[0, 1]),
            ("X//O[O]", &[]),
            (";polymer", &[0, 1]),
            (";solvent", &[]),
            ("/1/X", &[0, 1]),
            ("/2/A", &[]),
            ("B", &[]),
        ];
        for (cid, expected) in cif_rows {
            assert_eq!(c25_query(&cif, cid), expected, "cif {cid}");
        }
        // mmJSON projection: the canonical two-atom GLY fixture (label
        // chain A, AUTH chain X — auth precedence names the chain X),
        // rows in array order: 0=C2, 1=N1, auth_seq 7, model 1.
        let json = r#"{"data_demo":{"entry":{"id":["DEMO"]},"atom_site":{"id":[2,1],"group_PDB":["ATOM","ATOM"],"type_symbol":["C","N"],"label_atom_id":["C2","N1"],"label_alt_id":[null,null],"label_comp_id":["GLY","GLY"],"label_asym_id":["A","A"],"label_entity_id":["1","1"],"label_seq_id":[7,7],"pdbx_PDB_ins_code":[null,null],"Cartn_x":[3.25,-1.5],"Cartn_y":[4,8],"Cartn_z":[5,9],"auth_seq_id":[7,7],"auth_comp_id":["GLY","GLY"],"auth_asym_id":["X","X"],"pdbx_PDB_model_num":[1,1]}}}"#;
        let mmjson = c25_structure(json, cosmolkit_bio::BioCoordinateFormat::Mmjson);
        let json_rows: [(&str, &[u32]); 7] = [
            ("X/7", &[0, 1]), // auth chain name + stored auth seqid
            ("A/7", &[]),     // label chain name never stored
            ("X/(GLY)", &[0, 1]),
            ("X//C2[C]", &[0]),
            ("X//N1[N]", &[1]),
            (";polymer", &[]), // mmjson fixture carries no _entity.type
            ("/2/X", &[]),
        ];
        for (cid, expected) in json_rows {
            assert_eq!(c25_query(&mmjson, cid), expected, "mmjson {cid}");
        }
    }

    #[test]
    fn bio_cid_c21_relations_properties_and_values() {
        // Three relations x q/b x finite/signed-zero/carry values, as
        // exact bytes. Rows are fixed literals — no production helper
        // generates the expectations.
        let rows: [(&SelectionAtomInequality, &[u8]); 18] = [
            (
                &SelectionAtomInequality {
                    property: b'q',
                    relation: -1,
                    value: 0.5,
                },
                b";q<0.5",
            ),
            (
                &SelectionAtomInequality {
                    property: b'q',
                    relation: 0,
                    value: 0.5,
                },
                b";q=0.5",
            ),
            (
                &SelectionAtomInequality {
                    property: b'q',
                    relation: 1,
                    value: 0.5,
                },
                b";q>0.5",
            ),
            (
                &SelectionAtomInequality {
                    property: b'b',
                    relation: -3,
                    value: 1.5,
                },
                b";b<1.5",
            ),
            (
                &SelectionAtomInequality {
                    property: b'b',
                    relation: 0,
                    value: 2.5,
                },
                b";b=2.5",
            ),
            (
                &SelectionAtomInequality {
                    property: b'b',
                    relation: 7,
                    value: 9.5,
                },
                b";b>9.5",
            ),
            (
                &SelectionAtomInequality {
                    property: b'q',
                    relation: 1,
                    value: -0.0,
                },
                b";q>-0",
            ),
            (
                &SelectionAtomInequality {
                    property: b'q',
                    relation: -1,
                    value: 0.0,
                },
                b";q<0",
            ),
            (
                &SelectionAtomInequality {
                    property: b'b',
                    relation: 0,
                    value: 1023.0,
                },
                b";b=1023",
            ),
            (
                &SelectionAtomInequality {
                    property: b'q',
                    relation: 1,
                    value: 999_999_999.5,
                },
                b";q>1e+09",
            ),
            (
                &SelectionAtomInequality {
                    property: b'b',
                    relation: -1,
                    value: 1.23456789e18,
                },
                b";b<1.23456789e+18",
            ),
            (
                &SelectionAtomInequality {
                    property: b'q',
                    relation: 0,
                    value: 1e-5,
                },
                b";q=1e-05",
            ),
            (
                &SelectionAtomInequality {
                    property: b'b',
                    relation: 1,
                    value: 5e-324,
                },
                b";b>4.94065646e-324",
            ),
            (
                &SelectionAtomInequality {
                    property: b'q',
                    relation: -1,
                    value: f64::MAX,
                },
                b";q<1.79769313e+308",
            ),
            (
                &SelectionAtomInequality {
                    property: b'b',
                    relation: 0,
                    value: f64::INFINITY,
                },
                b";b=Inf",
            ),
            (
                &SelectionAtomInequality {
                    property: b'q',
                    relation: 1,
                    value: f64::NEG_INFINITY,
                },
                b";q>-Inf",
            ),
            (
                &SelectionAtomInequality {
                    property: b'b',
                    relation: -1,
                    value: f64::NAN,
                },
                b";b<NaN",
            ),
            (
                &SelectionAtomInequality {
                    property: b'q',
                    relation: 0,
                    value: -f64::NAN,
                },
                b";q=-NaN",
            ),
        ];
        for (inequality, expected) in rows {
            assert_eq!(inequality_str(inequality), expected);
        }
    }

    #[test]
    fn bio_cid_c21_complete_product_108() {
        // Full Cartesian product: q/b x relations -1/0/+1 x all 18
        // literal value/output pairs = 108 REAL helper calls. Numeric
        // outputs are independently fixed literals (never a production
        // formatter); expected bytes concatenate the literal
        // semicolon/property/relation bytes with them.
        let values: [(f64, &[u8]); 18] = [
            (0.5, b"0.5"),
            (1.5, b"1.5"),
            (2.5, b"2.5"),
            (9.5, b"9.5"),
            (-0.0, b"-0"),
            (0.0, b"0"),
            (1023.0, b"1023"),
            (999_999_999.5, b"1e+09"),
            (1.234_567_89e18, b"1.23456789e+18"),
            (1e-5, b"1e-05"),
            (5e-324, b"4.94065646e-324"),
            (f64::MAX, b"1.79769313e+308"),
            (f64::INFINITY, b"Inf"),
            (f64::NEG_INFINITY, b"-Inf"),
            (f64::NAN, b"NaN"),
            (-f64::NAN, b"-NaN"),
            (0.100_000_0005, b"0.100000001"),
            (1.000_000_005e18, b"1.00000001e+18"),
        ];
        let mut calls = 0_usize;
        for property in [b'q', b'b'] {
            for relation in [-1_i32, 0, 1] {
                let relation_byte = if relation == 0 {
                    b'='
                } else if relation < 0 {
                    b'<'
                } else {
                    b'>'
                };
                for (value, numeric) in values {
                    let out = inequality_str(&SelectionAtomInequality {
                        property,
                        relation,
                        value,
                    });
                    let mut expected = vec![b';', property, relation_byte];
                    expected.extend_from_slice(numeric);
                    assert_eq!(out, expected, "{} {relation} {value:?}", property as char);
                    calls += 1;
                }
            }
        }
        assert_eq!(calls, 108);
    }

    #[test]
    fn bio_cid_c21_all_property_bytes_768() {
        // Every u8 property byte x three relation signs with the literal
        // 0.5 output = 768 REAL helper calls; the property byte must be
        // preserved EXACTLY (raw byte, no substitution) for all 256
        // values including >= 0x80.
        let mut calls = 0_usize;
        for property in 0_u8..=255 {
            for relation in [-1_i32, 0, 1] {
                let relation_byte = if relation == 0 {
                    b'='
                } else if relation < 0 {
                    b'<'
                } else {
                    b'>'
                };
                let out = inequality_str(&SelectionAtomInequality {
                    property,
                    relation,
                    value: 0.5,
                });
                assert_eq!(
                    out.as_slice(),
                    &[b';', property, relation_byte, b'0', b'.', b'5'],
                    "{property:#04x} {relation}"
                );
                calls += 1;
            }
        }
        assert_eq!(calls, 768);
    }

    #[test]
    fn bio_cid_c21_rounding_nonroundtrip_and_raw_property_bytes() {
        // Rounding-induced non-roundtrip: a parse input with more than
        // nine significant digits serializes to the ROUNDED owner output,
        // never to the original input text (0.1000000005 -> "0.100000001",
        // 1.000000005e18 -> "1.00000001e+18").
        assert_eq!(
            inequality_str(&SelectionAtomInequality {
                property: b'q',
                relation: 1,
                value: 0.1000000005_f64,
            }),
            b";q>0.100000001"
        );
        assert_eq!(
            inequality_str(&SelectionAtomInequality {
                property: b'b',
                relation: 0,
                value: 1.000_000_005e18,
            }),
            b";b=1.00000001e+18"
        );
        // Raw property bytes: from_parts vocabulary can carry any u8;
        // bytes >= 0x80 are preserved EXACTLY (no lossy/Unicode
        // substitution — same C15-BYTES rule as the icode).
        assert_eq!(
            inequality_str(&SelectionAtomInequality {
                property: 0x80,
                relation: 1,
                value: 0.5
            })
            .as_slice(),
            b";\x80>0.5"
        );
        assert_eq!(
            inequality_str(&SelectionAtomInequality {
                property: 0xff,
                relation: -1,
                value: 0.0
            })
            .as_slice(),
            b";\xff<0"
        );
    }

    // ---- BIO-CID C02: fixed error-window regressions for the frozen
    // BioSelectionParseError boundary (expectations derived from
    // select.cpp:14-26 wrong_syntax, not from Display). ----

    #[test]
    fn bio_cid_c02_error_offset_is_real_or_absent() {
        // Syntax keeps the REAL recorded position (zero and nonzero).
        let zero = wrong_syntax("/a", 0, None);
        let mid = wrong_syntax("/aaa/bbb", 5, Some(" note"));
        let syntax_zero = BioSelectionParseError::Syntax(zero);
        let syntax_mid = BioSelectionParseError::Syntax(mid);
        assert_eq!(syntax_zero.offset(), Some(0));
        assert_eq!(syntax_mid.offset(), Some(5));
        // SeqidRange records NO offset: None, never a fabricated 0.
        let range = BioSelectionParseError::SeqidRange(SelectionSeqidRangeError {
            text: "99999999999".to_string(),
        });
        assert_eq!(range.offset(), None);
    }

    #[test]
    fn bio_cid_c02_source_message_bytes_exact_composition() {
        // Zero position: no window at all; info included before nothing.
        let zero = BioSelectionParseError::Syntax(wrong_syntax("/a/b", 0, Some(" in [...]")));
        assert_eq!(
            zero.source_message_bytes(),
            b"Invalid selection syntax in [...]: /a/b".to_vec()
        );
        // Nonzero position with an eight-byte window and info.
        let windowed =
            BioSelectionParseError::Syntax(wrong_syntax("/aaa/bbb/ccc", 5, Some(" note")));
        assert_eq!(
            windowed.source_message_bytes(),
            b"Invalid selection syntax note near \"bbb/ccc\": /aaa/bbb/ccc".to_vec()
        );
        // UTF-8 SPLIT window: the exact cut bytes appear, no U+FFFD, no
        // dropped bytes (cid '/' + 0xC3 0xA9 + 'abc'; pos 2 cuts inside
        // the two-byte character and the window takes exactly 8 bytes:
        // 0xA9 'a' 'b' 'c' is only 4, so the whole remainder is used).
        let cid = "/\u{e9}abc";
        let split = BioSelectionParseError::Syntax(wrong_syntax(cid, 2, None));
        assert_eq!(
            split.source_message_bytes(),
            b"Invalid selection syntax near \"\xA9abc\": /\xC3\xA9abc".to_vec()
        );
        // An eight-byte window that cuts a multibyte character contributes
        // only the bytes before the cut (6x two-byte characters; pos 1
        // starts inside the first one, so the window is A9 C3 A8 C3 A7
        // C3 A6 C3 — exactly eight bytes, cut mid-character).
        let wide = "\u{e9}\u{e8}\u{e7}\u{e6}\u{e5}\u{e4}";
        let cut = BioSelectionParseError::Syntax(wrong_syntax(wide, 1, None));
        let bytes = cut.source_message_bytes();
        let expected: Vec<u8> = [
            b"Invalid selection syntax near \"" as &[u8],
            &[0xA9, 0xC3, 0xA8, 0xC3, 0xA7, 0xC3, 0xA6, 0xC3],
            b"\": ",
            wide.as_bytes(),
        ]
        .concat();
        assert_eq!(bytes, expected);
        // SeqidRange variant: the carrier's own composed message bytes.
        let range = BioSelectionParseError::SeqidRange(SelectionSeqidRangeError {
            text: "-3000000000".to_string(),
        });
        assert_eq!(
            range.source_message_bytes(),
            b"CID sequence number out of representable range: -3000000000".to_vec()
        );
    }

    #[test]
    fn bio_cid_c02_syntax_error_round_trips_through_the_enum() {
        // A make_cid_list punctuation error keeps its exact bytes/pos/info
        // through the enum exactly as the standalone carrier reports.
        let direct = make_cid_list("/x/!*", 3, 5, "-[]()!/*.:;").unwrap_err();
        let wrapped = BioSelectionParseError::Syntax(direct.clone());
        assert_eq!(wrapped.offset(), Some(direct.pos));
        assert_eq!(
            wrapped.source_message_bytes(),
            direct.source_message_bytes()
        );
        assert!(wrapped.to_string().starts_with("Invalid selection syntax"));
    }
}
