//! Private selected mmCIF `_pdbx_database_status` category emitter
//! (BIO-DATABASE-STATUS1-32).
//!
//! Selected-category helper only: the coordinate composer is NOT
//! dispatched here and caller group dispatch remains ROOT-owned. Reuses
//! the one existing category setter; no second implementation.

use std::collections::BTreeMap;

use super::super::cif::CifBlock;

pub(super) fn write_database_status_category(
    info: &BTreeMap<String, String>,
    entry_id_raw: &str,
    block: &mut CifBlock,
) {
    // BEGIN GEMMI CPP FUNCTION gemmi::update_mmcif_block database_status branch (to_mmcif.cpp:474-481)
    // Gemmi❌❌:   if (groups.database_status) {
    // Gemmi✔️✔️:     auto initial_date = st.info.find("_pdbx_database_status.recvd_initial_deposition_date");
    // Gemmi✔️✔️:     if (initial_date != st.info.end() && !initial_date->second.empty()) {
    // Gemmi❗❌:       cif::ItemSpan span(block.items, "_pdbx_database_status.");
    // Gemmi❗❌:       span.set_pair("_pdbx_database_status.entry_id", id);
    // Gemmi❗❌:       span.set_pair(initial_date->first, initial_date->second);
    // Gemmi✔️✔️:     }
    // Gemmi❌❌:   }
    // END GEMMI CPP FUNCTION
    //
    // Guard scope: ONLY the outer `if (groups.database_status)` caller
    // dispatch bit is outside this selected-category boundary and is NOT
    // implemented here (composer undispatched; ROOT-owned open item) —
    // that omission alone carries ❌❌. The inner exact-map lookup and the
    // absent-or-empty predicate ARE modeled below.
    //
    // Helper bodies consulted (cross-file rule, NOT duplicated): the
    // ItemSpan prefix constructor + set_pair (cifdoc.hpp:613-645) are the
    // existing CifBlock::set_pair_in_category (first Pair with
    // case-insensitive equal tag rewritten in place with caller tag case
    // and fresh value position 0:0 keeping Pair.line on replacement; else
    // first Loop containing the tag replaced by a Pair — other loop
    // columns disappear, a SOURCE-defined loss; else new Pair emplaced at
    // the category span end). `st.info` is the case-SENSITIVE
    // std::map<std::string,std::string> (model.hpp:916); the current
    // source-state storage is the existing BTreeMap<String,String>
    // (structure_metadata.rs) — used as-is via ONE get_key_value
    // exact-key lookup, not moved, redefined, or duplicated.
    //
    // Behavior review: absent key or empty mapped string => NO category
    // scan, set, clear, or fallback — the WHOLE block (including any
    // conflicting parsed `_pdbx_database_status` pairs or a mixed
    // foreign-tag loop) is left byte-identically unchanged. Present
    // NONEMPTY => exactly TWO ordered existing-setter calls with
    // Some("_pdbx_database_status."): first the fixed
    // `_pdbx_database_status.entry_id` tag with the caller-supplied
    // entry raw forwarded verbatim (the source `id` is quote-processed
    // upstream at the composer boundary; this private boundary applies
    // no quoting/decoding/validation of its own), then the ACTUAL FOUND
    // KEY (exact map case, not the search literal) with the found raw
    // date bytes verbatim — no quote owner call, no date parser, no
    // trimming, case-folding, or second map model. Nonempty whitespace
    // and raw quote characters are written as-is: this private raw-value
    // forwarding is NOT a promise that such raw input serializes to a
    // grammar-valid CIF document.
    //
    // Cost review (separate axes): the source constructs ONE persistent
    // ItemSpan shared by both set_pair calls; each delegated Rust
    // set_pair_in_category call RECOMPUTES the category span — TWO span
    // recomputations versus one persistent span, a KNOWN local extra
    // scan (❌ on the span/set_pair lines above), not unresolved ❗. Each
    // raw value is a fresh Rust heap String where std::string may keep
    // short values in SSO inline storage (ABI-dependent, unmeasured,
    // separate caveat). The lookup line is ✔️ on BOTH axes: BTreeMap
    // get_key_value is, like std::map::find, an ordered-tree O(log N)
    // key-comparison lookup (key byte-comparison cost included in that
    // count) returning the same present/absent result and the same
    // found key/value — same lookup class, no invented O(1) claim. No
    // attempt is made to optimize or rewrite the accepted setter in this
    // packet.
    if let Some((found_key, raw_date)) =
        info.get_key_value("_pdbx_database_status.recvd_initial_deposition_date")
        && !raw_date.is_empty()
    {
        block.set_pair_in_category(
            Some("_pdbx_database_status."),
            "_pdbx_database_status.entry_id",
            entry_id_raw.to_owned(),
        );
        block.set_pair_in_category(Some("_pdbx_database_status."), found_key, raw_date.clone());
    }
}

#[cfg(test)]
mod bio_database_status_writer_tests {
    use super::write_database_status_category;
    use crate::cif::{CifBlock, CifCheckLevel, CifItem, read_cif_document};
    use std::collections::BTreeMap;

    const KEY: &str = "_pdbx_database_status.recvd_initial_deposition_date";
    const ENTRY_TAG: &str = "_pdbx_database_status.entry_id";

    const F0: &str = "data_f0\n_entry.id before\n_after.id tail\n";
    const F1: &str = "data_f1\n_entry.id before\n_pdbx_database_status.entry_id OLD\n_pdbx_database_status.recvd_initial_deposition_date old\n_after.id tail\n";
    const F2: &str = "data_f2\n_entry.id before\nloop_\n_pdbx_database_status.entry_id\n_pdbx_database_status.recvd_initial_deposition_date\nOLD old\n_after.id tail\n";

    fn parse_block(text: &str) -> CifBlock {
        read_cif_document(text, "database_status_fixture", CifCheckLevel::Syntax)
            .expect("frozen block state parses under the Syntax owner")
            .blocks()[0]
            .clone()
    }

    /// Frozen info profiles: literal key/value tuples in BTreeMap order.
    /// `active` marks a present NONEMPTY exact-key entry; `date` is the
    /// literal raw the call must publish.
    struct Profile {
        name: &'static str,
        literal: &'static [(&'static str, &'static str)],
        active: bool,
        date: &'static str,
    }

    const PROFILES: [Profile; 6] = [
        Profile {
            name: "D0",
            literal: &[("_unrelated.keep", "metadata")],
            active: false,
            date: "",
        },
        Profile {
            name: "D1",
            literal: &[
                ("_pdbx_database_status.recvd_initial_deposition_date", ""),
                ("_unrelated.keep", "metadata"),
            ],
            active: false,
            date: "",
        },
        Profile {
            name: "D2",
            literal: &[
                (
                    "_pdbx_database_status.recvd_initial_deposition_date",
                    "2026-10-04",
                ),
                ("_unrelated.keep", "metadata"),
            ],
            active: true,
            date: "2026-10-04",
        },
        Profile {
            name: "D3",
            literal: &[
                ("_pdbx_database_status.recvd_initial_deposition_date", " "),
                ("_unrelated.keep", "metadata"),
            ],
            active: true,
            date: " ",
        },
        Profile {
            name: "D4",
            literal: &[
                (
                    "_pdbx_database_status.recvd_initial_deposition_date",
                    "'2026 10'",
                ),
                ("_unrelated.keep", "metadata"),
            ],
            active: true,
            date: "'2026 10'",
        },
        Profile {
            name: "D5",
            literal: &[
                (
                    "_PDBX_DATABASE_STATUS.RECVD_INITIAL_DEPOSITION_DATE",
                    "2026-10-04",
                ),
                ("_unrelated.keep", "metadata"),
            ],
            active: false,
            date: "",
        },
    ];

    /// Complete literal map prerequisite inside the innermost loop.
    fn verify_info(info: &BTreeMap<String, String>, literal: &[(&str, &str)], label: &str) {
        let observed: Vec<(&str, &str)> = info
            .iter()
            .map(|(key, value)| (key.as_str(), value.as_str()))
            .collect();
        assert_eq!(observed, literal, "{label}: frozen info map literal");
    }

    /// Complete parsed-frame prerequisite. Returns (F1 entry Pair line,
    /// F1 date Pair line) for replacement-line expectations.
    fn verify_parsed(block: &CifBlock, name: &str) -> (usize, usize) {
        match name {
            "F0" => {
                assert_eq!(block.items().len(), 2, "F0 literal dimension");
                assert_eq!(block.name(), "f0");
                let CifItem::Pair(entry) = &block.items()[0] else {
                    panic!("F0 entry position");
                };
                assert_eq!(
                    (entry.tag(), entry.value().map(|v| v.raw())),
                    ("_entry.id", Some("before"))
                );
                assert_eq!(entry.line(), 2, "F0 entry Pair line");
                assert_eq!(entry.value().map(|v| (v.line(), v.column())), Some((2, 11)));
                let CifItem::Pair(tail) = &block.items()[1] else {
                    panic!("F0 tail position");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|v| v.raw())),
                    ("_after.id", Some("tail"))
                );
                assert_eq!(tail.line(), 3, "F0 tail Pair line");
                assert_eq!(tail.value().map(|v| (v.line(), v.column())), Some((3, 11)));
                (0, 0)
            }
            "F1" => {
                assert_eq!(block.items().len(), 4, "F1 literal dimension");
                assert_eq!(block.name(), "f1");
                let CifItem::Pair(entry) = &block.items()[0] else {
                    panic!("F1 entry position");
                };
                assert_eq!(
                    (entry.tag(), entry.value().map(|v| v.raw())),
                    ("_entry.id", Some("before"))
                );
                assert_eq!(entry.line(), 2, "F1 entry Pair line");
                assert_eq!(entry.value().map(|v| (v.line(), v.column())), Some((2, 11)));
                let CifItem::Pair(category) = &block.items()[1] else {
                    panic!("F1 category position");
                };
                assert_eq!(
                    (category.tag(), category.value().map(|v| v.raw())),
                    (ENTRY_TAG, Some("OLD"))
                );
                let entry_line = category.line();
                assert_eq!(entry_line, 3, "F1 category Pair line");
                assert_eq!(
                    category.value().map(|v| (v.line(), v.column())),
                    Some((3, 32))
                );
                let CifItem::Pair(date) = &block.items()[2] else {
                    panic!("F1 date position");
                };
                assert_eq!(
                    (date.tag(), date.value().map(|v| v.raw())),
                    (KEY, Some("old"))
                );
                let date_line = date.line();
                assert_eq!(date_line, 4, "F1 date Pair line");
                assert_eq!(date.value().map(|v| (v.line(), v.column())), Some((4, 53)));
                let CifItem::Pair(tail) = &block.items()[3] else {
                    panic!("F1 tail position");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|v| v.raw())),
                    ("_after.id", Some("tail"))
                );
                assert_eq!(tail.line(), 5, "F1 tail Pair line");
                assert_eq!(tail.value().map(|v| (v.line(), v.column())), Some((5, 11)));
                (entry_line, date_line)
            }
            "F2" => {
                assert_eq!(block.items().len(), 3, "F2 literal dimension");
                assert_eq!(block.name(), "f2");
                let CifItem::Pair(entry) = &block.items()[0] else {
                    panic!("F2 entry position");
                };
                assert_eq!(
                    (entry.tag(), entry.value().map(|v| v.raw())),
                    ("_entry.id", Some("before"))
                );
                assert_eq!(entry.line(), 2, "F2 entry Pair line");
                assert_eq!(entry.value().map(|v| (v.line(), v.column())), Some((2, 11)));
                let CifItem::Loop(stale) = &block.items()[1] else {
                    panic!("F2 stale loop position");
                };
                assert_eq!(stale.line(), 3, "F2 loop line");
                assert_eq!(stale.tags(), [ENTRY_TAG, KEY], "F2 loop tags");
                assert_eq!(stale.width(), 2, "F2 loop width");
                assert_eq!(stale.len(), 1, "F2 loop rows");
                assert_eq!(
                    stale.values().iter().map(|v| v.raw()).collect::<Vec<_>>(),
                    ["OLD", "old"],
                    "F2 loop raws"
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
                    panic!("F2 tail position");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|v| v.raw())),
                    ("_after.id", Some("tail"))
                );
                assert_eq!(tail.line(), 7, "F2 tail Pair line");
                assert_eq!(tail.value().map(|v| (v.line(), v.column())), Some((7, 11)));
                (0, 0)
            }
            other => panic!("unknown frozen block state {other}"),
        }
    }

    #[test]
    fn bio_database_status_writer_product72() {
        let entries: [&str; 2] = ["'ENTRY A'", "?"];
        let blocks: [(&str, &str); 3] = [("F0", F0), ("F1", F1), ("F2", F2)];

        let mut calls = 0usize;
        let mut active_calls = 0usize;
        let mut noop_calls = 0usize;
        let mut observed_category_pairs = 0usize;
        let mut observed_new_values = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        for profile in &PROFILES {
            for entry_raw in entries {
                for (block_name, block_text) in blocks {
                    for repetition in 0..2 {
                        let label =
                            format!("{}/{entry_raw}/{block_name}/{repetition}", profile.name);
                        // Complete literal info map prerequisite.
                        let mut info = BTreeMap::new();
                        info.insert("_unrelated.keep".to_string(), "metadata".to_string());
                        if let Some((key, value)) = profile
                            .literal
                            .iter()
                            .find(|(key, _)| *key != "_unrelated.keep")
                        {
                            info.insert(key.to_string(), value.to_string());
                        }
                        verify_info(&info, profile.literal, &label);
                        let entry_bytes = entry_raw.as_bytes().to_vec();

                        let mut block = parse_block(block_text);
                        let f1_lines = verify_parsed(&block, block_name);
                        let info_before = info.clone();
                        let before = block.clone();

                        write_database_status_category(&info, entry_raw, &mut block);
                        calls += 1;

                        // Same-call input preservation BEFORE output checks.
                        if info != info_before {
                            discrepancies.push(format!("{label}: info map mutated"));
                        }
                        if entry_raw.as_bytes() != entry_bytes {
                            discrepancies.push(format!("{label}: entry raw bytes mutated"));
                        }

                        // Name/duplicate/loop counts BEFORE missing-Pair
                        // branches, for every call.
                        if block.name() != before.name() {
                            discrepancies.push(format!("{label}: block name changed"));
                        }
                        let duplicate_entry = block
                            .items()
                            .iter()
                            .filter(|item| {
                                matches!(item, CifItem::Pair(p) if p
                                    .tag()
                                    .eq_ignore_ascii_case(ENTRY_TAG))
                            })
                            .count();
                        let duplicate_date = block
                            .items()
                            .iter()
                            .filter(|item| {
                                matches!(item, CifItem::Pair(p) if p.tag().eq_ignore_ascii_case(KEY))
                            })
                            .count();
                        let category_loops = block
                            .items()
                            .iter()
                            .filter(|item| {
                                matches!(item, CifItem::Loop(l) if l
                                    .tags()
                                    .first()
                                    .is_some_and(|tag| tag
                                        .starts_with("_pdbx_database_status.")))
                            })
                            .count();

                        if !profile.active {
                            // D0/D1/D5: whole block preserved byte-identically.
                            if block != before {
                                discrepancies.push(format!("{label}: no-op mutated whole block"));
                            }
                            noop_calls += 1;
                            continue;
                        }

                        active_calls += 1;
                        // Active POST layout per frozen fixture.
                        let (items_len, entry_index, date_index, unrelated, tail) = match block_name
                        {
                            "F0" => (4usize, 2usize, 3usize, (0usize, 0usize), (1usize, 1usize)),
                            "F1" => (4, 1, 2, (0, 0), (3, 3)),
                            _ => (4, 1, 2, (0, 0), (2, 3)),
                        };
                        if block.items().len() != items_len {
                            discrepancies.push(format!(
                                "{label}: item count {} != {items_len}",
                                block.items().len()
                            ));
                        }
                        if duplicate_entry != 1 || duplicate_date != 1 || category_loops != 0 {
                            discrepancies.push(format!(
                                "{label}: duplicates {duplicate_entry}/{duplicate_date} loops {category_loops}"
                            ));
                        }
                        // Unrelated mapping first.
                        if before.items().get(unrelated.0) != block.items().get(unrelated.1) {
                            discrepancies.push(format!("{label}: unrelated entry changed"));
                        }
                        if before.items().get(tail.0) != block.items().get(tail.1) {
                            discrepancies.push(format!("{label}: unrelated tail changed"));
                        }
                        for (index, tag, expected_raw, expected_line) in [
                            (entry_index, ENTRY_TAG, entry_raw, {
                                if block_name == "F1" { f1_lines.0 } else { 0 }
                            }),
                            (date_index, KEY, profile.date, {
                                if block_name == "F1" { f1_lines.1 } else { 0 }
                            }),
                        ] {
                            match block.items().get(index) {
                                Some(CifItem::Pair(pair)) => {
                                    if pair.tag() != tag {
                                        discrepancies
                                            .push(format!("{label}: {tag} tag {}", pair.tag()));
                                    }
                                    if pair.value().map(|v| v.raw()) != Some(expected_raw) {
                                        discrepancies.push(format!(
                                            "{label}: {tag} raw {:?} != {expected_raw}",
                                            pair.value().map(|v| v.raw())
                                        ));
                                    }
                                    if pair.line() != expected_line {
                                        discrepancies.push(format!(
                                            "{label}: {tag} line {} != {expected_line}",
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
                                _ => {
                                    discrepancies.push(format!("{label}: no {tag} pair at {index}"))
                                }
                            }
                        }
                        // OBSERVED censuses from the actual block.
                        let pairs_now = block
                            .items()
                            .iter()
                            .filter(|item| {
                                matches!(item, CifItem::Pair(p) if p
                                    .tag()
                                    .starts_with("_pdbx_database_status."))
                            })
                            .count();
                        observed_category_pairs += pairs_now;
                        let values_now = block
                            .items()
                            .iter()
                            .filter_map(|item| match item {
                                CifItem::Pair(p)
                                    if p.tag().starts_with("_pdbx_database_status.") =>
                                {
                                    p.value().map(|v| v.raw())
                                }
                                _ => None,
                            })
                            .filter(|raw| *raw == entry_raw || *raw == profile.date)
                            .count();
                        observed_new_values += values_now;
                    }
                }
            }
        }

        assert_eq!(
            calls, 72,
            "exact 72 real write_database_status_category calls"
        );
        assert_eq!(active_calls, 36, "36 active publications");
        assert_eq!(noop_calls, 36, "36 source no-op whole-block proofs");
        assert_eq!(
            observed_category_pairs, 72,
            "observed active category pair census"
        );
        assert_eq!(observed_new_values, 72, "observed active raw value census");
        assert!(
            discrepancies.is_empty(),
            "product discrepancies: {discrepancies:?}"
        );
    }

    #[test]
    fn bio_database_status_writer_sequence8() {
        // ONE parsed F1 per repetition; ordered D2, D0, D4, D1 with entry
        // `'ENTRY A'`. D2/D4 are active publications (each writes the entry
        // pair and the date pair); D0/D1 are source no-ops that retain the
        // whole previous output. The PRE oracle is LITERAL-ONLY: the first
        // call's PRE is the complete parsed F1 frame; later PREs use the
        // previous PUBLISHED literal frame, advanced ONLY on an active
        // publication — never from actual output capture.
        struct Step {
            profile: &'static Profile,
        }
        let sequence: [Step; 4] = [
            Step {
                profile: &PROFILES[2],
            }, // D2
            Step {
                profile: &PROFILES[0],
            }, // D0
            Step {
                profile: &PROFILES[4],
            }, // D4
            Step {
                profile: &PROFILES[1],
            }, // D1
        ];
        const ENTRY_RAW: &str = "'ENTRY A'";

        let mut calls = 0usize;
        let mut active_calls = 0usize;
        let mut noop_calls = 0usize;
        let mut emitted_values = 0usize;
        let mut observed_post_values = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        for repetition in 0..2 {
            let mut block = parse_block(F1);
            let original = block.clone();
            // Literal-only previous PUBLISHED frame: (entry raw, date raw).
            // Starts at the parsed F1 category raws for the first PRE.
            let mut previous_literal: (&str, &str) = ("OLD", "old");
            for (step_index, step) in sequence.iter().enumerate() {
                let profile = step.profile;
                let label = format!(
                    "seq/{}/{}{}/{}",
                    profile.name, "D-seq-", step_index, repetition
                );
                // Complete literal PRE frame BEFORE every call.
                if step_index == 0 {
                    verify_parsed(&block, "F1");
                } else {
                    assert_eq!(block.items().len(), 4, "{label}: PRE items");
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
                    let CifItem::Pair(entry) = &block.items()[1] else {
                        panic!("{label}: PRE entry pair");
                    };
                    assert_eq!(entry.tag(), ENTRY_TAG, "{label}: PRE entry tag");
                    assert_eq!(
                        entry.value().map(|v| v.raw()),
                        Some(previous_literal.0),
                        "{label}: PRE entry raw"
                    );
                    assert_eq!(entry.line(), 3, "{label}: PRE entry Pair line");
                    assert_eq!(
                        entry.value().map(|v| (v.line(), v.column())),
                        Some((0, 0)),
                        "{label}: PRE entry value position"
                    );
                    let CifItem::Pair(date) = &block.items()[2] else {
                        panic!("{label}: PRE date pair");
                    };
                    assert_eq!(date.tag(), KEY, "{label}: PRE date tag");
                    assert_eq!(
                        date.value().map(|v| v.raw()),
                        Some(previous_literal.1),
                        "{label}: PRE date raw"
                    );
                    assert_eq!(date.line(), 4, "{label}: PRE date Pair line");
                    assert_eq!(
                        date.value().map(|v| (v.line(), v.column())),
                        Some((0, 0)),
                        "{label}: PRE date value position"
                    );
                    let CifItem::Pair(tail) = &block.items()[3] else {
                        panic!("{label}: PRE tail");
                    };
                    assert_eq!(
                        (
                            tail.tag(),
                            tail.value().map(|v| v.raw()),
                            tail.line(),
                            tail.value().map(|v| (v.line(), v.column()))
                        ),
                        ("_after.id", Some("tail"), 5, Some((5, 11)))
                    );
                }
                // Fresh per-call info map/entry/block snapshots.
                let mut info = BTreeMap::new();
                info.insert("_unrelated.keep".to_string(), "metadata".to_string());
                if let Some((key, value)) = profile
                    .literal
                    .iter()
                    .find(|(key, _)| *key != "_unrelated.keep")
                {
                    info.insert(key.to_string(), value.to_string());
                }
                verify_info(&info, profile.literal, &label);
                let info_before = info.clone();
                let entry_bytes = ENTRY_RAW.as_bytes().to_vec();
                let before = block.clone();

                write_database_status_category(&info, ENTRY_RAW, &mut block);
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

                if !profile.active {
                    // D0/D1: whole block retained byte-identically.
                    if block != before {
                        discrepancies.push(format!("{label}: no-op mutated whole block"));
                    }
                    noop_calls += 1;
                } else {
                    active_calls += 1;
                    // Active F1 layout: replacement keeps Pair.line 3/4,
                    // fresh value positions 0:0, unrelated/tail unchanged.
                    if block.items().len() != 4 {
                        discrepancies
                            .push(format!("{label}: item count {} != 4", block.items().len()));
                    }
                    for (index, tag, expected_raw, expected_line) in [
                        (1usize, ENTRY_TAG, ENTRY_RAW, 3usize),
                        (2usize, KEY, profile.date, 4usize),
                    ] {
                        match block.items().get(index) {
                            Some(CifItem::Pair(pair)) => {
                                if pair.tag() != tag {
                                    discrepancies
                                        .push(format!("{label}: {tag} tag {}", pair.tag()));
                                }
                                if pair.value().map(|v| v.raw()) != Some(expected_raw) {
                                    discrepancies.push(format!(
                                        "{label}: {tag} raw {:?} != {expected_raw}",
                                        pair.value().map(|v| v.raw())
                                    ));
                                }
                                if pair.line() != expected_line {
                                    discrepancies.push(format!(
                                        "{label}: {tag} line {} != {expected_line}",
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
                                // Emitted census: values THIS call writes.
                                if pair.value().map(|v| v.raw()) == Some(expected_raw) {
                                    emitted_values += 1;
                                }
                            }
                            _ => discrepancies.push(format!("{label}: no {tag} pair at {index}")),
                        }
                    }
                    if before.items().get(0) != block.items().get(0) {
                        discrepancies.push(format!("{label}: unrelated entry changed"));
                    }
                    if before.items().get(3) != block.items().get(3) {
                        discrepancies.push(format!("{label}: unrelated tail changed"));
                    }
                    let duplicates = block
                        .items()
                        .iter()
                        .filter(|item| {
                            matches!(item, CifItem::Pair(p) if p
                                .tag()
                                .starts_with("_pdbx_database_status."))
                        })
                        .count();
                    if duplicates != 2 {
                        discrepancies.push(format!("{label}: category pairs {duplicates}"));
                    }
                    // Advance the literal oracle ONLY on active publication.
                    previous_literal = (ENTRY_RAW, profile.date);
                }
                // Observed POST census for EVERY call (retained no-op rows
                // count; distinct from the emitted census above).
                let values_now = block
                    .items()
                    .iter()
                    .filter_map(|item| match item {
                        CifItem::Pair(p) if p.tag().starts_with("_pdbx_database_status.") => {
                            p.value().map(|v| v.raw())
                        }
                        _ => None,
                    })
                    .count();
                observed_post_values += values_now;
                // Original unrelated frame never disturbed across the run.
                if block.items().get(0) != original.items().get(0) {
                    discrepancies.push(format!("{label}: original unrelated entry changed"));
                }
            }
        }

        assert_eq!(calls, 8, "exact 8 sequence calls");
        assert_eq!(active_calls, 4, "four active publications");
        assert_eq!(noop_calls, 4, "four source no-ops");
        assert_eq!(emitted_values, 8, "8 newly written values");
        assert_eq!(observed_post_values, 16, "16 observed POST values");
        assert!(
            discrepancies.is_empty(),
            "sequence discrepancies: {discrepancies:?}"
        );
    }

    #[test]
    fn bio_database_status_writer_mixed4() {
        const MIXED: &str = "data_bad\n_entry.id before\nloop_\n_pdbx_database_status.entry_id\n_foreign.keep\nOLD stale\n_after.id tail\n";
        // D0 (absent key) and D2 (present nonempty) on the mixed
        // foreign-tag loop. D0 retains the WHOLE parsed block. D2's first
        // set_pair replaces the mixed loop by the entry Pair — the
        // `_foreign.keep` column DISAPPEARS: this is the source's
        // ItemSpan::set_pair loop-replacement semantics, an explicitly
        // SOURCE-defined information loss, not a structured error and not
        // a lossless claim. The setter boundary is infallible: all four
        // calls succeed; no InvalidLoop/source-chain assertions exist.
        let profiles: [&Profile; 2] = [&PROFILES[0], &PROFILES[2]]; // D0, D2
        const ENTRY_RAW: &str = "'ENTRY A'";

        let mut calls = 0usize;
        let mut active_calls = 0usize;
        let mut noop_calls = 0usize;
        let mut observed_values = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        for profile in profiles {
            for repetition in 0..2 {
                let label = format!("mixed/{}/{}", profile.name, repetition);
                let mut info = BTreeMap::new();
                info.insert("_unrelated.keep".to_string(), "metadata".to_string());
                if let Some((key, value)) = profile
                    .literal
                    .iter()
                    .find(|(key, _)| *key != "_unrelated.keep")
                {
                    info.insert(key.to_string(), value.to_string());
                }
                verify_info(&info, profile.literal, &label);
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
                    [ENTRY_TAG, "_foreign.keep"],
                    "{label}: loop tags"
                );
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

                write_database_status_category(&info, ENTRY_RAW, &mut block);
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
                        matches!(item, CifItem::Pair(p) if p.tag().eq_ignore_ascii_case("_foreign.keep"))
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

                if !profile.active {
                    if block != before {
                        discrepancies.push(format!("{label}: D0 mutated whole block"));
                    }
                    if foreign_pairs + foreign_loops != 1 {
                        discrepancies.push(format!("{label}: D0 foreign column not retained"));
                    }
                    noop_calls += 1;
                    continue;
                }

                active_calls += 1;
                if block.items().len() != 4 {
                    discrepancies.push(format!("{label}: item count {} != 4", block.items().len()));
                }
                // SOURCE-defined loss: the foreign column is GONE after the
                // mixed loop is replaced by the entry Pair.
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
                    [(1usize, ENTRY_TAG, ENTRY_RAW), (2usize, KEY, profile.date)]
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
                // Observed census from the actual block.
                let values_now = block
                    .items()
                    .iter()
                    .filter_map(|item| match item {
                        CifItem::Pair(p) if p.tag().starts_with("_pdbx_database_status.") => {
                            p.value().map(|v| v.raw())
                        }
                        _ => None,
                    })
                    .count();
                observed_values += values_now;
            }
        }

        assert_eq!(calls, 4, "exact 4 mixed calls, all infallible successes");
        assert_eq!(active_calls, 2, "two active publications");
        assert_eq!(noop_calls, 2, "two whole-block no-ops");
        assert_eq!(observed_values, 4, "observed active value census");
        assert!(
            discrepancies.is_empty(),
            "mixed discrepancies: {discrepancies:?}"
        );
    }
}
