//! Private legacy mmCIF `_struct`/`_struct_keywords` category emitter
//! (BIO-TITLE-KEYWORDS1-32).
//!
//! Complete selected-category helper only: the coordinate composer is NOT
//! dispatched here and caller group dispatch remains ROOT-owned. Reuses the
//! ONE existing quote owner and the ONE existing category setter; no
//! copied quoting/setter algorithm.

use std::collections::BTreeMap;

use super::super::cif::{CifBlock, quote_cif_value};

pub(super) fn write_title_keywords_category(
    info: &BTreeMap<String, String>,
    entry_id_raw: &str,
    block: &mut CifBlock,
) {
    // BEGIN GEMMI CPP FUNCTION gemmi::update_mmcif_block title_keywords branch (to_mmcif.cpp:943-959)
    // Gemmi❌❌:   if (groups.title_keywords) {
    // Gemmi✔️✔️:     auto title = st.info.find("_struct.title");
    // Gemmi✔️✔️:     if (title != st.info.end()) {
    // Gemmi❗❌:       cif::ItemSpan span(block.items, "_struct.");
    // Gemmi❗❌:       span.set_pair("_struct.entry_id", id);
    // Gemmi❗❌:       span.set_pair(title->first, cif::quote(title->second));
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:     auto pdbx_keywords = st.info.find("_struct_keywords.pdbx_keywords");
    // Gemmi✔️✔️:     auto keywords = st.info.find("_struct_keywords.text");
    // Gemmi❗❌:     cif::ItemSpan span(block.items, "_struct_keywords.");
    // Gemmi✔️✔️:     if (pdbx_keywords != st.info.end() || keywords != st.info.end())
    // Gemmi❗❌:       span.set_pair("_struct_keywords.entry_id", id);
    // Gemmi✔️✔️:     if (pdbx_keywords != st.info.end())
    // Gemmi❗❌:       span.set_pair(pdbx_keywords->first, cif::quote(pdbx_keywords->second));
    // Gemmi✔️✔️:     if (keywords != st.info.end())
    // Gemmi❗❌:       span.set_pair(keywords->first, cif::quote(keywords->second));
    // Gemmi❌❌:   }
    // END GEMMI CPP FUNCTION
    //
    // Guard scope: ONLY the outer `if (groups.title_keywords)` caller
    // dispatch bit is outside this complete selected-category boundary and
    // is NOT implemented here (coordinate-only composer undispatched;
    // ROOT-owned open item) — that omission alone carries ❌❌. Every inner
    // branch of the source (title presence guard, keywords OR-guard and
    // both keyword branches) is modeled below.
    //
    // Helper bodies consulted (cross-file rule, NOT duplicated): the
    // ItemSpan prefix constructor + set_pair (cifdoc.hpp:613-645) are the
    // existing CifBlock::set_pair_in_category (first Pair with
    // case-insensitive equal tag rewritten in place with caller tag case
    // and fresh value position 0:0 keeping Pair.line on replacement; else
    // first Loop containing the tag replaced by a Pair — other loop
    // columns disappear, a SOURCE-defined loss; else new Pair emplaced at
    // the category span end). `cif::quote` (cifdoc.hpp:1180-1197) is the
    // existing pub quote_cif_value owner (ordinary/null branch, empty ->
    // `''`, quote priority ' then " then ;) — invoked ONCE per quoted
    // field, never copied. `st.info` is the case-SENSITIVE
    // std::map<std::string,std::string> (model.hpp:916); the current
    // source-state storage is the existing BTreeMap<String,String>
    // (structure_metadata.rs), used as-is via THREE exact-case
    // get_key_value lookups — not moved, redefined, or duplicated.
    //
    // Behavior review: PRESENCE alone controls publication (unlike the
    // database-status nonempty guard) — a present EMPTY title/keyword is
    // written as the quoted empty raw `''`. If T is present: `_struct.`
    // entry pair with the caller-supplied entry raw forwarded verbatim
    // (source `id` is quote-processed upstream; no re-quoting here), then
    // the ACTUAL FOUND T KEY (exact map case) with quote_cif_value(value).
    // Then for keywords: if P OR K is present, ONE `_struct_keywords.`
    // entry pair with the entry raw; then P's actual found key with quoted
    // value if present; then K's actual found key with quoted value if
    // present — exactly the source's THREE ordered set_pair calls after
    // the OR-guard. No trimming, case-folding, fallback, or deletion of
    // absent fields; absent keys leave their categories completely
    // untouched (mask0 leaves the WHOLE block unchanged). Infallible
    // existing-setter boundary: no invented Result/error path.
    //
    // Cost review (separate axes): the source constructs ONE `_struct.`
    // span shared by up to two set_pair calls and ONE `_struct_keywords.`
    // span shared by up to three; every delegated Rust
    // set_pair_in_category call RECOMPUTES its category span — up to FIVE
    // span recomputations versus TWO persistent constructions, a KNOWN
    // local extra scan (❌ on the span/set_pair lines above), not
    // unresolved ❗. The source constructs the `_struct_keywords.` span
    // even when NEITHER P nor K exists; this Rust owner performs NO scan
    // in that both-absent case — a local behavior-preserving improvement
    // for that path only (the unused span produced no observable effect),
    // NOT a blanket algorithm/performance promotion; the used-path
    // recomputation cost above stays ❌. Each quoted/raw value is a fresh
    // Rust heap String where std::string may keep short values in SSO
    // inline storage (ABI-dependent, unmeasured, separate caveat). The
    // three lookups are ✔️ on BOTH axes: BTreeMap get_key_value is, like
    // std::map::find, an ordered-tree O(log N) key-comparison lookup
    // (key byte-comparison cost included) returning the same
    // present/absent result — same lookup class, no invented O(1) claim.
    // No attempt is made to optimize or rewrite the accepted setter or
    // quote owner in this packet, and passing fixed regression frames is
    // not a claim that all arbitrary interleaved CIF layouts are
    // equivalent.
    let title = info.get_key_value("_struct.title");
    if let Some((title_key, title_value)) = title {
        block.set_pair_in_category(
            Some("_struct."),
            "_struct.entry_id",
            entry_id_raw.to_owned(),
        );
        block.set_pair_in_category(
            Some("_struct."),
            title_key,
            quote_cif_value(title_value.clone()),
        );
    }
    let pdbx_keywords = info.get_key_value("_struct_keywords.pdbx_keywords");
    let keywords = info.get_key_value("_struct_keywords.text");
    if pdbx_keywords.is_some() || keywords.is_some() {
        block.set_pair_in_category(
            Some("_struct_keywords."),
            "_struct_keywords.entry_id",
            entry_id_raw.to_owned(),
        );
    }
    if let Some((pdbx_key, pdbx_value)) = pdbx_keywords {
        block.set_pair_in_category(
            Some("_struct_keywords."),
            pdbx_key,
            quote_cif_value(pdbx_value.clone()),
        );
    }
    if let Some((keywords_key, keywords_value)) = keywords {
        block.set_pair_in_category(
            Some("_struct_keywords."),
            keywords_key,
            quote_cif_value(keywords_value.clone()),
        );
    }
}

#[cfg(test)]
mod bio_title_keywords_writer_tests {
    use super::write_title_keywords_category;
    use crate::cif::{CifBlock, CifCheckLevel, CifItem, read_cif_document};
    use std::collections::BTreeMap;

    const KEY_T: &str = "_struct.title";
    const KEY_P: &str = "_struct_keywords.pdbx_keywords";
    const KEY_K: &str = "_struct_keywords.text";
    const TAG_SE: &str = "_struct.entry_id";
    const TAG_KE: &str = "_struct_keywords.entry_id";

    const F0: &str = "data_f0\n_entry.id before\n_after.id tail\n";
    const F1: &str = "data_f1\n_entry.id before\n_struct.entry_id OLD\n_struct.title oldtitle\n_struct_keywords.entry_id OLD\n_struct_keywords.pdbx_keywords oldp\n_struct_keywords.text oldtext\n_after.id tail\n";

    fn parse_block(text: &str) -> CifBlock {
        read_cif_document(text, "title_keywords_fixture", CifCheckLevel::Syntax)
            .expect("frozen block state parses under the Syntax owner")
            .blocks()[0]
            .clone()
    }

    /// Frozen quoting profiles: decoded (t/p/k) and expected quoted raws.
    /// Expected raws are frozen LITERALS here, never derived by calling
    /// the quote helper during expectation construction.
    struct Quoting {
        name: &'static str,
        t: &'static str,
        p: &'static str,
        k: &'static str,
        qt: &'static str,
        qp: &'static str,
        qk: &'static str,
    }

    const QUOTINGS: [Quoting; 3] = [
        Quoting {
            name: "Q0",
            t: "",
            p: "",
            k: "",
            qt: "''",
            qp: "''",
            qk: "''",
        },
        Quoting {
            name: "Q1",
            t: "TITLE",
            p: "KEY",
            k: "TEXT",
            qt: "TITLE",
            qp: "KEY",
            qk: "TEXT",
        },
        Quoting {
            name: "Q2",
            t: "TITLE A",
            p: "KEY A",
            k: "TEXT A",
            qt: "'TITLE A'",
            qp: "'KEY A'",
            qk: "'TEXT A'",
        },
    ];

    fn mask_bits(mask: u8) -> (bool, bool, bool) {
        (mask & 1 != 0, mask & 2 != 0, mask & 4 != 0)
    }

    /// Build the frozen info map for (mask, quoting) and return its
    /// literal tuple list in BTreeMap order for whole-map verification.
    fn info_for(
        mask: u8,
        quoting: &Quoting,
    ) -> (BTreeMap<String, String>, Vec<(&'static str, &'static str)>) {
        let mut info = BTreeMap::new();
        for (shadow_key, shadow_value) in [
            ("_STRUCT.TITLE", "SHADOW"),
            ("_STRUCT_KEYWORDS.PDBX_KEYWORDS", "SHADOW"),
            ("_STRUCT_KEYWORDS.TEXT", "SHADOW"),
        ] {
            info.insert(shadow_key.to_string(), shadow_value.to_string());
        }
        let (has_t, has_p, has_k) = mask_bits(mask);
        if has_t {
            info.insert(KEY_T.to_string(), quoting.t.to_string());
        }
        if has_p {
            info.insert(KEY_P.to_string(), quoting.p.to_string());
        }
        if has_k {
            info.insert(KEY_K.to_string(), quoting.k.to_string());
        }
        info.insert("_unrelated.keep".to_string(), "metadata".to_string());
        // Sorted BTreeMap order (uppercase < lowercase): the three SHADOW
        // controls, then _struct.title, _struct_keywords.pdbx_keywords,
        // _struct_keywords.text, then _unrelated.keep.
        let mut literal = vec![
            ("_STRUCT.TITLE", "SHADOW"),
            ("_STRUCT_KEYWORDS.PDBX_KEYWORDS", "SHADOW"),
            ("_STRUCT_KEYWORDS.TEXT", "SHADOW"),
        ];
        if has_t {
            literal.push((KEY_T, quoting.t));
        }
        if has_p {
            literal.push((KEY_P, quoting.p));
        }
        if has_k {
            literal.push((KEY_K, quoting.k));
        }
        literal.push(("_unrelated.keep", "metadata"));
        (info, literal)
    }

    /// Complete literal map prerequisite inside the innermost loop.
    fn verify_info(info: &BTreeMap<String, String>, literal: &[(&str, &str)], label: &str) {
        let observed: Vec<(&str, &str)> = info
            .iter()
            .map(|(key, value)| (key.as_str(), value.as_str()))
            .collect();
        assert_eq!(observed, literal, "{label}: frozen info map literal");
    }

    /// Complete parsed-frame prerequisite for F0/F1.
    fn verify_parsed(block: &CifBlock, name: &str) {
        match name {
            "F0" => {
                assert_eq!(block.items().len(), 2, "F0 literal dimension");
                assert_eq!(block.name(), "f0");
                let CifItem::Pair(entry) = &block.items()[0] else {
                    panic!("F0 entry position");
                };
                assert_eq!(
                    (
                        entry.tag(),
                        entry.value().map(|v| v.raw()),
                        entry.line(),
                        entry.value().map(|v| (v.line(), v.column()))
                    ),
                    ("_entry.id", Some("before"), 2, Some((2, 11)))
                );
                let CifItem::Pair(tail) = &block.items()[1] else {
                    panic!("F0 tail position");
                };
                assert_eq!(
                    (
                        tail.tag(),
                        tail.value().map(|v| v.raw()),
                        tail.line(),
                        tail.value().map(|v| (v.line(), v.column()))
                    ),
                    ("_after.id", Some("tail"), 3, Some((3, 11)))
                );
            }
            "F1" => {
                assert_eq!(block.items().len(), 7, "F1 literal dimension");
                assert_eq!(block.name(), "f1");
                let expected: [(&str, &str, usize, (usize, usize)); 7] = [
                    ("_entry.id", "before", 2, (2, 11)),
                    (TAG_SE, "OLD", 3, (3, 18)),
                    (KEY_T, "oldtitle", 4, (4, 15)),
                    (TAG_KE, "OLD", 5, (5, 27)),
                    (KEY_P, "oldp", 6, (6, 32)),
                    (KEY_K, "oldtext", 7, (7, 23)),
                    ("_after.id", "tail", 8, (8, 11)),
                ];
                for (index, (tag, raw, line, position)) in expected.iter().enumerate() {
                    let CifItem::Pair(pair) = &block.items()[index] else {
                        panic!("F1 item {index} position");
                    };
                    assert_eq!(
                        (
                            pair.tag(),
                            pair.value().map(|v| v.raw()),
                            pair.line(),
                            pair.value().map(|v| (v.line(), v.column()))
                        ),
                        (*tag, Some(*raw), *line, Some(*position)),
                        "F1 item {index}"
                    );
                }
            }
            other => panic!("unknown frozen block state {other}"),
        }
    }

    #[test]
    fn bio_title_keywords_writer_product96() {
        let entries: [&str; 2] = ["'ENTRY A'", "?"];
        let blocks: [(&str, &str); 2] = [("F0", F0), ("F1", F1)];

        let mut calls = 0usize;
        let mut active_calls = 0usize;
        let mut noop_calls = 0usize;
        let mut observed_new_values = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        for quoting in &QUOTINGS {
            for entry_raw in entries {
                for (block_name, block_text) in blocks {
                    for mask in 0u8..8 {
                        let (has_t, has_p, has_k) = mask_bits(mask);
                        let label = format!("{}/mask{mask}/{entry_raw}/{block_name}", quoting.name);
                        let (info, literal) = info_for(mask, quoting);
                        verify_info(&info, &literal, &label);
                        let info_before = info.clone();
                        let entry_bytes = entry_raw.as_bytes().to_vec();

                        let mut block = parse_block(block_text);
                        verify_parsed(&block, block_name);
                        let before = block.clone();

                        write_title_keywords_category(&info, entry_raw, &mut block);
                        calls += 1;

                        // Same-call input preservation BEFORE output checks.
                        if info != info_before {
                            discrepancies.push(format!("{label}: info map mutated"));
                        }
                        if entry_raw.as_bytes() != entry_bytes {
                            discrepancies.push(format!("{label}: entry raw bytes mutated"));
                        }
                        if block.name() != before.name() {
                            discrepancies.push(format!("{label}: block name changed"));
                        }

                        // Duplicate/loop counts BEFORE missing-Pair branches.
                        let category_tags = [TAG_SE, KEY_T, TAG_KE, KEY_P, KEY_K];
                        let mut duplicates = [0usize; 5];
                        for (field, tag) in category_tags.iter().enumerate() {
                            duplicates[field] = block
                                .items()
                                .iter()
                                .filter(|item| {
                                    matches!(item, CifItem::Pair(p) if p.tag().eq_ignore_ascii_case(tag))
                                })
                                .count();
                            if duplicates[field] > 1 {
                                discrepancies.push(format!(
                                    "{label}: {tag} duplicates {}",
                                    duplicates[field]
                                ));
                            }
                        }
                        let category_loops = block
                            .items()
                            .iter()
                            .filter(|item| {
                                matches!(item, CifItem::Loop(l) if l.tags().first().is_some_and(
                                    |tag| {
                                        tag.starts_with("_struct.")
                                            || tag.starts_with("_struct_keywords.")
                                    }
                                ))
                            })
                            .count();
                        if category_loops != 0 {
                            discrepancies.push(format!(
                                "{label}: unexpected category loop {category_loops}"
                            ));
                        }

                        if mask == 0 {
                            // Neither category key exists (SHADOWs are
                            // case-sensitive misses): WHOLE block no-op.
                            if block != before {
                                discrepancies.push(format!("{label}: mask0 mutated whole block"));
                            }
                            noop_calls += 1;
                            continue;
                        }

                        active_calls += 1;
                        // Source-order written sequence for this mask.
                        let mut written: Vec<(&str, &str)> = Vec::new();
                        if has_t {
                            written.push((TAG_SE, entry_raw));
                            written.push((KEY_T, quoting.qt));
                        }
                        if has_p || has_k {
                            written.push((TAG_KE, entry_raw));
                        }
                        if has_p {
                            written.push((KEY_P, quoting.qp));
                        }
                        if has_k {
                            written.push((KEY_K, quoting.qk));
                        }

                        match block_name {
                            "F0" => {
                                // Appends after tail, in source order.
                                let expected_len = 2 + written.len();
                                if block.items().len() != expected_len {
                                    discrepancies.push(format!(
                                        "{label}: item count {} != {expected_len}",
                                        block.items().len()
                                    ));
                                }
                                if before.items().get(0) != block.items().get(0) {
                                    discrepancies.push(format!("{label}: unrelated entry changed"));
                                }
                                if before.items().get(1) != block.items().get(1) {
                                    discrepancies.push(format!("{label}: unrelated tail changed"));
                                }
                                for (offset, (tag, raw)) in written.iter().enumerate() {
                                    match block.items().get(2 + offset) {
                                        Some(CifItem::Pair(pair)) => {
                                            if pair.tag() != *tag {
                                                discrepancies.push(format!(
                                                    "{label}: append {offset} tag {}",
                                                    pair.tag()
                                                ));
                                            }
                                            if pair.value().map(|v| v.raw()) != Some(*raw) {
                                                discrepancies.push(format!(
                                                    "{label}: append {offset} raw {:?}",
                                                    pair.value().map(|v| v.raw())
                                                ));
                                            }
                                            if pair.line() != 0 {
                                                discrepancies.push(format!(
                                                    "{label}: append {offset} line {}",
                                                    pair.line()
                                                ));
                                            }
                                            if let Some(value) = pair.value() {
                                                if value.line() != 0 || value.column() != 0 {
                                                    discrepancies.push(format!(
                                                        "{label}: append {offset} value {}:{}",
                                                        value.line(),
                                                        value.column()
                                                    ));
                                                }
                                            }
                                        }
                                        _ => discrepancies
                                            .push(format!("{label}: no append at {}", 2 + offset)),
                                    }
                                }
                            }
                            _ => {
                                // F1: seven items, frozen order; only
                                // selected fields change; unselected stay
                                // WHOLE-equal parsed originals.
                                if block.items().len() != 7 {
                                    discrepancies.push(format!(
                                        "{label}: item count {} != 7",
                                        block.items().len()
                                    ));
                                }
                                if before.items().get(0) != block.items().get(0) {
                                    discrepancies.push(format!("{label}: unrelated entry changed"));
                                }
                                if before.items().get(6) != block.items().get(6) {
                                    discrepancies.push(format!("{label}: unrelated tail changed"));
                                }
                                let fields: [(usize, &str, bool, &str); 5] = [
                                    (1, TAG_SE, has_t, entry_raw),
                                    (2, KEY_T, has_t, quoting.qt),
                                    (3, TAG_KE, has_p || has_k, entry_raw),
                                    (4, KEY_P, has_p, quoting.qp),
                                    (5, KEY_K, has_k, quoting.qk),
                                ];
                                for (index, tag, selected, raw) in fields {
                                    if selected {
                                        match block.items().get(index) {
                                            Some(CifItem::Pair(pair)) => {
                                                if pair.tag() != tag {
                                                    discrepancies.push(format!(
                                                        "{label}: {tag} tag {}",
                                                        pair.tag()
                                                    ));
                                                }
                                                if pair.value().map(|v| v.raw()) != Some(raw) {
                                                    discrepancies.push(format!(
                                                        "{label}: {tag} raw {:?} != {raw}",
                                                        pair.value().map(|v| v.raw())
                                                    ));
                                                }
                                                // Replacement keeps parsed line.
                                                if pair.line() != (index as usize + 2) {
                                                    discrepancies.push(format!(
                                                        "{label}: {tag} line {}",
                                                        pair.line()
                                                    ));
                                                }
                                                if let Some(value) = pair.value() {
                                                    if value.line() != 0 || value.column() != 0 {
                                                        discrepancies.push(format!(
                                                            "{label}: {tag} value {}:{}",
                                                            value.line(),
                                                            value.column()
                                                        ));
                                                    }
                                                }
                                            }
                                            _ => discrepancies
                                                .push(format!("{label}: no {tag} pair at {index}")),
                                        }
                                    } else if block.items().get(index) != before.items().get(index)
                                    {
                                        discrepancies
                                            .push(format!("{label}: unselected {tag} changed"));
                                    }
                                }
                            }
                        }
                        // OBSERVED NEW-value census: count actual expected-
                        // tag pairs carrying this call's written raws.
                        for (tag, raw) in &written {
                            let count = block
                                .items()
                                .iter()
                                .filter(|item| {
                                    matches!(item, CifItem::Pair(p) if p
                                        .tag()
                                        .eq_ignore_ascii_case(tag)
                                        && p.value().map(|v| v.raw()) == Some(*raw))
                                })
                                .count();
                            if count != 1 {
                                discrepancies
                                    .push(format!("{label}: {tag} observed count {count}"));
                            }
                            observed_new_values += count;
                        }
                    }
                }
            }
        }

        assert_eq!(
            calls, 96,
            "exact 96 real write_title_keywords_category calls"
        );
        assert_eq!(active_calls, 84, "84 active publications");
        assert_eq!(noop_calls, 12, "12 whole-block no-ops (mask0)");
        assert_eq!(
            observed_new_values, 264,
            "observed NEW value census (mask profile [0,2,2,4,2,4,3,5] x12)"
        );
        assert!(
            discrepancies.is_empty(),
            "product discrepancies: {discrepancies:?}"
        );
    }

    #[test]
    fn bio_title_keywords_writer_sequence8() {
        // ONE parsed F1 per repetition; ordered (mask7,Q0),(mask0,Q1),
        // (mask1,Q2),(mask6,Q1) with entry `'ENTRY A'`. The PRE oracle is
        // LITERAL-ONLY: the first PRE is the complete parsed F1 frame;
        // later PREs use the previous PUBLISHED literal five-category
        // frame, advanced ONLY on an active publication (the mask0 no-op
        // advances nothing). Q0 writes quoted empty `''` for present
        // empty values — presence, not emptiness, controls publication.
        const ENTRY_RAW: &str = "'ENTRY A'";
        struct SeqStep {
            mask: u8,
            quoting: usize,
            // Five written (tag, raw) pairs in source order; empty slice
            // marks the no-op step.
            written: &'static [(&'static str, &'static str)],
        }
        const SEQ: [SeqStep; 4] = [
            SeqStep {
                mask: 7,
                quoting: 0,
                written: &[
                    (TAG_SE, "'ENTRY A'"),
                    (KEY_T, "''"),
                    (TAG_KE, "'ENTRY A'"),
                    (KEY_P, "''"),
                    (KEY_K, "''"),
                ],
            },
            SeqStep {
                mask: 0,
                quoting: 1,
                written: &[],
            },
            SeqStep {
                mask: 1,
                quoting: 2,
                written: &[(TAG_SE, "'ENTRY A'"), (KEY_T, "'TITLE A'")],
            },
            SeqStep {
                mask: 6,
                quoting: 1,
                written: &[(TAG_KE, "'ENTRY A'"), (KEY_P, "KEY"), (KEY_K, "TEXT")],
            },
        ];

        let mut calls = 0usize;
        let mut active_calls = 0usize;
        let mut noop_calls = 0usize;
        let mut emitted_values = 0usize;
        let mut observed_post_values = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        for repetition in 0..2 {
            let mut block = parse_block(F1);
            let original = block.clone();
            // Literal-only previous PUBLISHED five-category frame:
            // (struct-entry, title, keywords-entry, pdbx, text). Starts
            // at the parsed F1 raws for the first PRE.
            let mut previous_literal: [&str; 5] = ["OLD", "oldtitle", "OLD", "oldp", "oldtext"];
            for (step_index, step) in SEQ.iter().enumerate() {
                let quoting = &QUOTINGS[step.quoting];
                let label = format!(
                    "seq/mask{}/{}/{}/{}",
                    step.mask, quoting.name, step_index, repetition
                );
                // Complete literal PRE frame BEFORE every call.
                if step_index == 0 {
                    verify_parsed(&block, "F1");
                } else {
                    assert_eq!(block.items().len(), 7, "{label}: PRE items");
                    assert_eq!(block.name(), "f1", "{label}: PRE name");
                    let CifItem::Pair(unrelated) = &block.items()[0] else {
                        panic!("{label}: PRE unrelated entry");
                    };
                    assert_eq!(
                        (
                            unrelated.tag(),
                            unrelated.value().map(|v| v.raw()),
                            unrelated.line(),
                            unrelated.value().map(|v| (v.line(), v.column()))
                        ),
                        ("_entry.id", Some("before"), 2, Some((2, 11)))
                    );
                    let fields = [
                        (1usize, TAG_SE),
                        (2usize, KEY_T),
                        (3usize, TAG_KE),
                        (4usize, KEY_P),
                        (5usize, KEY_K),
                    ];
                    for (offset, (index, tag)) in fields.iter().enumerate() {
                        let CifItem::Pair(pair) = &block.items()[*index] else {
                            panic!("{label}: PRE field {tag}");
                        };
                        assert_eq!(pair.tag(), *tag, "{label}: PRE {tag} tag");
                        assert_eq!(
                            pair.value().map(|v| v.raw()),
                            Some(previous_literal[offset]),
                            "{label}: PRE {tag} raw"
                        );
                        assert_eq!(pair.line(), *index + 2, "{label}: PRE {tag} Pair line");
                        assert_eq!(
                            pair.value().map(|v| (v.line(), v.column())),
                            Some((0, 0)),
                            "{label}: PRE {tag} value position"
                        );
                    }
                    let CifItem::Pair(tail) = &block.items()[6] else {
                        panic!("{label}: PRE tail");
                    };
                    assert_eq!(
                        (
                            tail.tag(),
                            tail.value().map(|v| v.raw()),
                            tail.line(),
                            tail.value().map(|v| (v.line(), v.column()))
                        ),
                        ("_after.id", Some("tail"), 8, Some((8, 11)))
                    );
                }
                // Fresh per-call info map/entry/block snapshots.
                let (info, literal) = info_for(step.mask, quoting);
                verify_info(&info, &literal, &label);
                let info_before = info.clone();
                let entry_bytes = ENTRY_RAW.as_bytes().to_vec();
                let before = block.clone();

                write_title_keywords_category(&info, ENTRY_RAW, &mut block);
                calls += 1;

                // Same-call input preservation BEFORE output checks.
                if info != info_before {
                    discrepancies.push(format!("{label}: info map mutated"));
                }
                if ENTRY_RAW.as_bytes() != entry_bytes {
                    discrepancies.push(format!("{label}: entry raw bytes mutated"));
                }
                if block.name() != before.name() {
                    discrepancies.push(format!("{label}: block name changed"));
                }
                if before.items().get(0) != block.items().get(0) {
                    discrepancies.push(format!("{label}: unrelated entry changed"));
                }
                if before.items().get(6) != block.items().get(6) {
                    discrepancies.push(format!("{label}: unrelated tail changed"));
                }
                if block.items().get(0) != original.items().get(0) {
                    discrepancies.push(format!("{label}: original unrelated entry changed"));
                }

                if step.written.is_empty() {
                    // mask0: whole block retained byte-identically.
                    if block != before {
                        discrepancies.push(format!("{label}: no-op mutated whole block"));
                    }
                    noop_calls += 1;
                } else {
                    active_calls += 1;
                    if block.items().len() != 7 {
                        discrepancies
                            .push(format!("{label}: item count {} != 7", block.items().len()));
                    }
                    let selected: Vec<usize> = step
                        .written
                        .iter()
                        .map(|(tag, _)| match *tag {
                            TAG_SE => 1,
                            KEY_T => 2,
                            TAG_KE => 3,
                            KEY_P => 4,
                            _ => 5,
                        })
                        .collect();
                    for index in 1..=5usize {
                        let tag = [TAG_SE, KEY_T, TAG_KE, KEY_P, KEY_K][index - 1];
                        if selected.contains(&index) {
                            let expected = step
                                .written
                                .iter()
                                .find(|(written_tag, _)| *written_tag == tag)
                                .expect("selected field has a written pair")
                                .1;
                            match block.items().get(index) {
                                Some(CifItem::Pair(pair)) => {
                                    if pair.tag() != tag {
                                        discrepancies
                                            .push(format!("{label}: {tag} tag {}", pair.tag()));
                                    }
                                    if pair.value().map(|v| v.raw()) != Some(expected) {
                                        discrepancies.push(format!(
                                            "{label}: {tag} raw {:?} != {expected}",
                                            pair.value().map(|v| v.raw())
                                        ));
                                    }
                                    if pair.line() != index + 2 {
                                        discrepancies
                                            .push(format!("{label}: {tag} line {}", pair.line()));
                                    }
                                    if let Some(value) = pair.value() {
                                        if value.line() != 0 || value.column() != 0 {
                                            discrepancies.push(format!(
                                                "{label}: {tag} value {}:{}",
                                                value.line(),
                                                value.column()
                                            ));
                                        }
                                    }
                                    if pair.value().map(|v| v.raw()) == Some(expected) {
                                        emitted_values += 1;
                                    }
                                }
                                _ => {
                                    discrepancies.push(format!("{label}: no {tag} pair at {index}"))
                                }
                            }
                        } else if block.items().get(index) != before.items().get(index) {
                            discrepancies.push(format!("{label}: unselected {tag} changed"));
                        }
                    }
                    // Advance the literal oracle ONLY on active
                    // publication, from THIS call's frozen literals.
                    let mut next = previous_literal;
                    for (tag, raw) in step.written {
                        let index = match *tag {
                            TAG_SE => 0,
                            KEY_T => 1,
                            TAG_KE => 2,
                            KEY_P => 3,
                            _ => 4,
                        };
                        next[index] = raw;
                    }
                    previous_literal = next;
                }
                // Observed POST census for EVERY call (retained fields
                // count; distinct from the emitted census above).
                let values_now = block
                    .items()
                    .iter()
                    .filter(|item| {
                        matches!(item, CifItem::Pair(p) if [TAG_SE, KEY_T, TAG_KE, KEY_P, KEY_K]
                            .iter()
                            .any(|tag| p.tag().eq_ignore_ascii_case(tag)))
                    })
                    .count();
                observed_post_values += values_now;
            }
        }

        assert_eq!(calls, 8, "exact 8 sequence calls");
        assert_eq!(active_calls, 6, "six active publications");
        assert_eq!(noop_calls, 2, "two whole-block no-ops");
        assert_eq!(emitted_values, 20, "20 newly written values (5+0+2+3)x2");
        assert_eq!(observed_post_values, 40, "40 observed POST values");
        assert!(
            discrepancies.is_empty(),
            "sequence discrepancies: {discrepancies:?}"
        );
    }

    #[test]
    fn bio_title_keywords_writer_mixed4() {
        const MIXED: &str = "data_bad\n_entry.id before\nloop_\n_struct.entry_id\n_foreign.keep\nOLD stale\n_after.id tail\n";
        // (mask0,Q1) and (mask1,Q1) on the mixed foreign-tag loop, entry
        // `'ENTRY A'`. mask0 (no exact keys) retains the WHOLE parsed
        // block. mask1 (T present) drives the FIRST set_pair through the
        // mixed loop: the source's ItemSpan::set_pair loop-replacement
        // semantics REPLACE the whole loop (including the `_foreign.keep`
        // column) with the struct-entry Pair; the title Pair is then
        // inserted. The foreign column DISAPPEARS — an explicitly
        // SOURCE-defined information loss, not a structured error and not
        // a lossless claim. The boundary is infallible: no fabricated
        // InvalidLoop/source-chain assertions.
        const ENTRY_RAW: &str = "'ENTRY A'";
        let profiles: [(u8, usize); 2] = [(0, 1), (1, 1)]; // (mask, Q1)

        let mut calls = 0usize;
        let mut active_calls = 0usize;
        let mut noop_calls = 0usize;
        let mut observed_new_values = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        for (mask, quoting_index) in profiles {
            for repetition in 0..2 {
                let quoting = &QUOTINGS[quoting_index];
                let label = format!("mixed/mask{mask}/{}/{}", quoting.name, repetition);
                let (info, literal) = info_for(mask, quoting);
                verify_info(&info, &literal, &label);
                let info_before = info.clone();
                let entry_bytes = ENTRY_RAW.as_bytes().to_vec();

                let mut block = parse_block(MIXED);
                // Complete parsed prerequisites on the exact mixed fixture.
                assert_eq!(block.items().len(), 3, "{label}: items");
                assert_eq!(block.name(), "bad", "{label}: name");
                let CifItem::Pair(entry) = &block.items()[0] else {
                    panic!("{label}: entry position");
                };
                assert_eq!(
                    (
                        entry.tag(),
                        entry.value().map(|v| v.raw()),
                        entry.line(),
                        entry.value().map(|v| (v.line(), v.column()))
                    ),
                    ("_entry.id", Some("before"), 2, Some((2, 11)))
                );
                let CifItem::Loop(stale) = &block.items()[1] else {
                    panic!("{label}: loop position");
                };
                assert_eq!(stale.line(), 3, "{label}: loop line");
                assert_eq!(
                    stale.tags(),
                    [TAG_SE, "_foreign.keep"],
                    "{label}: loop tags"
                );
                assert_eq!(stale.width(), 2, "{label}: loop width");
                assert_eq!(stale.len(), 1, "{label}: loop rows");
                assert_eq!(
                    stale.values().iter().map(|v| v.raw()).collect::<Vec<_>>(),
                    ["OLD", "stale"],
                    "{label}: loop raws"
                );
                assert_eq!(
                    (stale.values()[0].line(), stale.values()[0].column()),
                    (6, 1)
                );
                assert_eq!(
                    (stale.values()[1].line(), stale.values()[1].column()),
                    (6, 5)
                );
                let CifItem::Pair(tail) = &block.items()[2] else {
                    panic!("{label}: tail position");
                };
                assert_eq!(
                    (
                        tail.tag(),
                        tail.value().map(|v| v.raw()),
                        tail.line(),
                        tail.value().map(|v| (v.line(), v.column()))
                    ),
                    ("_after.id", Some("tail"), 7, Some((7, 11)))
                );
                let before = block.clone();

                write_title_keywords_category(&info, ENTRY_RAW, &mut block);
                calls += 1;

                // Same-call input preservation BEFORE output checks.
                if info != info_before {
                    discrepancies.push(format!("{label}: info map mutated"));
                }
                if ENTRY_RAW.as_bytes() != entry_bytes {
                    discrepancies.push(format!("{label}: entry raw bytes mutated"));
                }
                if block.name() != before.name() {
                    discrepancies.push(format!("{label}: block name changed"));
                }

                // Foreign-column disappearance bookkeeping BEFORE branches.
                let foreign_pairs = block
                    .items()
                    .iter()
                    .filter(|item| {
                        matches!(item, CifItem::Pair(p) if p
                            .tag()
                            .eq_ignore_ascii_case("_foreign.keep"))
                    })
                    .count();
                let foreign_loops = block
                    .items()
                    .iter()
                    .filter(|item| {
                        matches!(item, CifItem::Loop(l) if l
                            .tags()
                            .iter()
                            .any(|tag| tag.eq_ignore_ascii_case("_foreign.keep")))
                    })
                    .count();

                if mask == 0 {
                    if block != before {
                        discrepancies.push(format!("{label}: mask0 mutated whole block"));
                    }
                    if foreign_pairs + foreign_loops != 1 {
                        discrepancies.push(format!("{label}: mask0 foreign column not retained"));
                    }
                    noop_calls += 1;
                    continue;
                }

                active_calls += 1;
                if block.items().len() != 4 {
                    discrepancies.push(format!("{label}: item count {} != 4", block.items().len()));
                }
                // SOURCE-defined loss: the foreign column is GONE after
                // the mixed loop is replaced by the struct-entry Pair.
                if foreign_pairs + foreign_loops != 0 {
                    discrepancies.push(format!("{label}: foreign column unexpectedly retained"));
                }
                if before.items().get(0) != block.items().get(0) {
                    discrepancies.push(format!("{label}: unrelated entry changed"));
                }
                if before.items().get(2) != block.items().get(3) {
                    discrepancies.push(format!("{label}: old tail 2->3 changed"));
                }
                for (index, tag, expected_raw) in
                    [(1usize, TAG_SE, ENTRY_RAW), (2usize, KEY_T, quoting.qt)]
                {
                    match block.items().get(index) {
                        Some(CifItem::Pair(pair)) => {
                            if pair.tag() != tag {
                                discrepancies.push(format!("{label}: {tag} tag {}", pair.tag()));
                            }
                            if pair.value().map(|v| v.raw()) != Some(expected_raw) {
                                discrepancies.push(format!(
                                    "{label}: {tag} raw {:?} != {expected_raw}",
                                    pair.value().map(|v| v.raw())
                                ));
                            }
                            if pair.line() != 0 {
                                discrepancies
                                    .push(format!("{label}: {tag} line {} != 0", pair.line()));
                            }
                            if let Some(value) = pair.value() {
                                if value.line() != 0 || value.column() != 0 {
                                    discrepancies.push(format!(
                                        "{label}: {tag} value {}:{}",
                                        value.line(),
                                        value.column()
                                    ));
                                }
                            }
                        }
                        _ => discrepancies.push(format!("{label}: no {tag} pair at {index}")),
                    }
                }
                // Observed NEW-value census from the actual block.
                for tag in [TAG_SE, KEY_T] {
                    let count = block
                        .items()
                        .iter()
                        .filter(|item| {
                            matches!(item, CifItem::Pair(p) if p.tag().eq_ignore_ascii_case(tag))
                        })
                        .count();
                    if count != 1 {
                        discrepancies.push(format!("{label}: {tag} count {count}"));
                    }
                    observed_new_values += count;
                }
            }
        }

        assert_eq!(calls, 4, "exact 4 mixed calls, all infallible successes");
        assert_eq!(active_calls, 2, "two active publications");
        assert_eq!(noop_calls, 2, "two whole-block no-ops");
        assert_eq!(observed_new_values, 4, "observed NEW value census");
        assert!(
            discrepancies.is_empty(),
            "mixed discrepancies: {discrepancies:?}"
        );
    }
}
