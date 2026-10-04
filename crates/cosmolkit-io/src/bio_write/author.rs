//! Private selected mmCIF `_audit_author` category emitter
//! (BIO-AUTHOR-WRITE1-32).
//!
//! Selected-category helper only: the coordinate composer is NOT
//! dispatched here and caller group dispatch remains ROOT-owned. Reuses
//! the one existing loop initialization/row-append chain and the accepted
//! quote owner.

use super::super::cif::{CifBlock, CifReadError};

pub(super) fn write_author_category(
    authors: &[String],
    block: &mut CifBlock,
) -> Result<(), CifReadError> {
    // BEGIN GEMMI CPP FUNCTION gemmi::write_author_category selected branch (to_mmcif.cpp:483-488)
    // Gemmi❌❌:   if (groups.author && !st.meta.authors.empty()) {
    // Gemmi❗❌:     cif::Loop& loop = block.init_mmcif_loop("_audit_author.", {"pdbx_ordinal", "name"});
    // Gemmi✔️✔️:     int n = 0;
    // Gemmi✔️✔️:     for (const std::string& author : st.meta.authors)
    // Gemmi❗❌:       loop.add_row({std::to_string(++n), cif::quote(author)});
    // Gemmi❌❌:   }
    // END GEMMI CPP FUNCTION
    //
    // Guard scope: ONLY the caller `groups.author` bit is outside this
    // selected-category boundary and is NOT implemented here (composer
    // undispatched) — that omission alone carries ❌❌. The source's
    // `!st.meta.authors.empty()` predicate IS IN SCOPE and is modeled
    // below: the early return runs for an empty slice, so empty input
    // does not touch the block at all — no loop initialization, no
    // category validation, and a conflicting existing `_audit_author`
    // category (including a mixed foreign-tag loop) is preserved
    // byte-identically. The borrowed `authors: &[String]` maps directly
    // to the source's borrowed iteration over `st.meta.authors`
    // (BioMetadata's `pub authors: Vec<String>` is the same
    // vector-of-strings storage; no copy or ownership change).
    //
    // Helper bodies consulted (cross-file rule, NOT duplicated): the
    // initializer_list loop tags are the two existing
    // init_mmcif_loop("_audit_author.", ["pdbx_ordinal","name"]) strings;
    // `cif::quote` is the existing quote_cif_value owner (ordinary/
    // null/quote-priority branches, empty input becomes `''`);
    // add_row's InvalidLoop width error propagates unchanged via `?`.
    //
    // Behavior review: for a NONEMPTY supplied slice the category is
    // reinitialized through the ONE existing init owner (any stale
    // parsed category pairs/loop are replaced), then rows iterate the
    // SUPPLIED order with NO deduplication/sorting/normalization or
    // silent fallback; each row is ONE add_row whose first value is the
    // 1-based ordinal `(index + 1)` rendered by the existing integer
    // Display (the source's `std::to_string(++n)` — counter restarts at
    // 1 on EVERY invocation; no persisted ordinal state) and whose second
    // value is the author quoted EXACTLY ONCE via quote_cif_value.
    // Structural InvalidLoop errors from the delegated init propagate
    // unchanged. Valid inputs are the source-defined domain: author
    // count <= i32::MAX (the source `++int` beyond that has no defined
    // behavior; no overflow test/policy/guard is invented here). No new
    // numeric parser.
    //
    // Cost review (separate axes): ONE delegated category scan
    // (init_mmcif_loop) plus ONE O(n + total author text) iteration,
    // the same shape as the source. Per row the Rust path allocates a
    // temporary Vec plus two owned Strings (ordinal and quoted name)
    // where the source passes a stack initializer_list of std::string
    // values whose allocation count is SSO/ABI-dependent — a KNOWN
    // extra allocation, recorded as ❌ on the two marked lines, not a
    // default unresolved ❗. No whole-cost, public, lossless, default,
    // history or generic float-state claim.
    if authors.is_empty() {
        return Ok(());
    }
    let loop_ = block.init_mmcif_loop("_audit_author.", &["pdbx_ordinal", "name"])?;
    for (index, author) in authors.iter().enumerate() {
        loop_.add_row(vec![
            (index + 1).to_string(),
            crate::cif::quote_cif_value(author.clone()),
        ])?;
    }
    Ok(())
}

#[cfg(test)]
mod bio_author_writer_tests {
    use super::write_author_category;
    use crate::cif::{CifBlock, CifCheckLevel, CifItem, CifReadErrorKind, read_cif_document};

    const AUTHOR_TAGS: [&str; 2] = ["_audit_author.pdbx_ordinal", "_audit_author.name"];

    fn parse_block(text: &str) -> CifBlock {
        read_cif_document(text, "author_writer_fixture", CifCheckLevel::Syntax)
            .expect("frozen block state parses under the Syntax owner")
            .blocks()[0]
            .clone()
    }

    const F0: &str = "data_f0\n_entry.id before\n_after.id tail\n";
    const F1: &str = "data_f1\n_entry.id before\n_audit_author.pdbx_ordinal OLD\n_audit_author.name old\n_after.id tail\n";
    const F2: &str = "data_f2\n_entry.id before\nloop_\n_audit_author.pdbx_ordinal\n_audit_author.name\nOLD old\n_after.id tail\n";

    /// Setup-time parsed-frame prerequisite per frozen state with the
    /// full lexical-position contract (source lines/values 1-based).
    /// Returns the parsed Loop.line for F2 (retained by reinit).
    fn verify_block_state(block: &CifBlock, name: &str) -> usize {
        let items = block.items();
        let CifItem::Pair(entry) = &items[0] else {
            panic!("{name}: entry pair position");
        };
        assert_eq!(entry.tag(), "_entry.id", "{name}: entry tag");
        assert_eq!(
            entry.value().map(|v| v.raw()),
            Some("before"),
            "{name}: entry raw"
        );
        assert_eq!(entry.line(), 2, "{name}: entry line");
        assert_eq!(
            entry.value().map(|v| (v.line(), v.column())),
            Some((2, 11)),
            "{name}: entry value position"
        );
        match name {
            "F0" => {
                assert_eq!(items.len(), 2, "F0 dimension");
                assert_eq!(block.name(), "f0");
                let CifItem::Pair(tail) = &items[1] else {
                    panic!("F0 tail");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|v| v.raw())),
                    ("_after.id", Some("tail"))
                );
                assert_eq!(tail.line(), 3);
                assert_eq!(tail.value().map(|v| (v.line(), v.column())), Some((3, 11)));
                0
            }
            "F1" => {
                assert_eq!(items.len(), 4, "F1 dimension");
                assert_eq!(block.name(), "f1");
                let CifItem::Pair(old_ord) = &items[1] else {
                    panic!("F1 old ordinal");
                };
                assert_eq!(
                    (old_ord.tag(), old_ord.value().map(|v| v.raw())),
                    ("_audit_author.pdbx_ordinal", Some("OLD"))
                );
                assert_eq!(old_ord.line(), 3);
                assert_eq!(
                    old_ord.value().map(|v| (v.line(), v.column())),
                    Some((3, 28))
                );
                let CifItem::Pair(old_name) = &items[2] else {
                    panic!("F1 old name");
                };
                assert_eq!(
                    (old_name.tag(), old_name.value().map(|v| v.raw())),
                    ("_audit_author.name", Some("old"))
                );
                assert_eq!(old_name.line(), 4);
                assert_eq!(
                    old_name.value().map(|v| (v.line(), v.column())),
                    Some((4, 20))
                );
                let CifItem::Pair(tail) = &items[3] else {
                    panic!("F1 tail");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|v| v.raw())),
                    ("_after.id", Some("tail"))
                );
                assert_eq!(tail.line(), 5);
                assert_eq!(tail.value().map(|v| (v.line(), v.column())), Some((5, 11)));
                0
            }
            "F2" => {
                assert_eq!(items.len(), 3, "F2 dimension");
                assert_eq!(block.name(), "f2");
                let CifItem::Loop(stale) = &items[1] else {
                    panic!("F2 stale loop");
                };
                assert_eq!(stale.tags(), AUTHOR_TAGS, "F2 tags");
                assert_eq!(stale.line(), 3, "F2 loop line");
                let raws: Vec<&str> = stale.values().iter().map(|v| v.raw()).collect();
                assert_eq!(raws, ["OLD", "old"], "F2 raws");
                assert_eq!(stale.values()[0].line(), 6, "F2 first value line");
                assert_eq!(stale.values()[0].column(), 1);
                assert_eq!(stale.values()[1].line(), 6);
                assert_eq!(stale.values()[1].column(), 5);
                let CifItem::Pair(tail) = &items[2] else {
                    panic!("F2 tail");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|v| v.raw())),
                    ("_after.id", Some("tail"))
                );
                assert_eq!(tail.line(), 7);
                assert_eq!(tail.value().map(|v| (v.line(), v.column())), Some((7, 11)));
                stale.line()
            }
            other => panic!("unknown frozen state {other}"),
        }
    }

    /// Non-panicking output collection for one successful nonempty call.
    /// Returns (rows, values, complete loops) observed in the output.
    fn collect_output(
        block: &CifBlock,
        before: &CifBlock,
        expected_raws: &[&str],
        expected_loop_line: usize,
        state: &str,
        label: &str,
        discrepancies: &mut Vec<String>,
    ) -> (usize, usize, usize) {
        let items = block.items();
        if items.len() != 3 {
            discrepancies.push(format!("{label}: item count {} != 3", items.len()));
        }
        // name / duplicate-category checks run BEFORE any missing-loop branch.
        if block.name() != before.name() {
            discrepancies.push(format!("{label}: block name changed"));
        }
        // Independent checks must also run when no expected loop survived.
        for tag in AUTHOR_TAGS {
            let pair_count = items
                .iter()
                .filter(
                    |item| matches!(item, CifItem::Pair(p) if p.tag().eq_ignore_ascii_case(tag)),
                )
                .count();
            if pair_count != 0 {
                discrepancies.push(format!("{label}: stale pair {tag} x{pair_count}"));
            }
        }
        let (tail_before, tail_after) = match state {
            "F0" => (1, 1),
            "F1" => (3, 2),
            _ => (2, 2),
        };
        if before.items().get(0) != items.get(0) {
            discrepancies.push(format!("{label}: unrelated entry changed"));
        }
        if before.items().get(tail_before) != items.get(tail_after) {
            discrepancies.push(format!("{label}: unrelated tail changed"));
        }
        let loop_positions: Vec<usize> = items
            .iter()
            .enumerate()
            .filter(|(_, item)| {
                matches!(item, CifItem::Loop(row) if row
                    .tags()
                    .first()
                    .is_some_and(|tag| tag.starts_with("_audit_author.")))
            })
            .map(|(index, _)| index)
            .collect();
        if loop_positions.len() != 1 {
            discrepancies.push(format!(
                "{label}: audit_author loop count {}",
                loop_positions.len()
            ));
            return (0, 0, 0);
        }
        let loop_index = loop_positions[0];
        let expected_index = if state == "F0" { 2 } else { 1 };
        if loop_index != expected_index {
            discrepancies.push(format!(
                "{label}: loop index {loop_index} != {expected_index}"
            ));
        }
        let CifItem::Loop(row) = &items[loop_index] else {
            unreachable!()
        };
        if row.tags() != AUTHOR_TAGS {
            discrepancies.push(format!("{label}: tags {:?}", row.tags()));
        }
        if row.line() != expected_loop_line {
            discrepancies.push(format!(
                "{label}: loop line {} != {expected_loop_line}",
                row.line()
            ));
        }
        if row.values().len() != expected_raws.len() {
            discrepancies.push(format!(
                "{label}: value count {} != {}",
                row.values().len(),
                expected_raws.len()
            ));
        }
        for (index, expected) in expected_raws.iter().enumerate() {
            match row.values().get(index) {
                Some(value) => {
                    if value.raw() != *expected {
                        discrepancies.push(format!(
                            "{label}: value {index} {:?} != {expected}",
                            value.raw()
                        ));
                    }
                    if value.line() != 0 || value.column() != 0 {
                        discrepancies.push(format!(
                            "{label}: value {index} position {}:{}",
                            value.line(),
                            value.column()
                        ));
                    }
                }
                None => discrepancies.push(format!("{label}: missing value {index}")),
            }
        }
        let complete = row.tags() == AUTHOR_TAGS && row.values().len() % 2 == 0;
        if !complete {
            discrepancies.push(format!("{label}: incomplete two-column loop"));
        }
        (
            row.values().len() / 2,
            row.values().len(),
            usize::from(complete),
        )
    }

    #[test]
    fn bio_author_writer_product36() {
        let profiles: [(&str, Vec<&str>, Vec<&str>); 6] = [
            ("A0", Vec::new(), Vec::new()),
            ("A1", vec![""], vec!["1", "''"]),
            ("A2", vec!["A B"], vec!["1", "'A B'"]),
            ("A3", vec!["?"], vec!["1", "'?'"]),
            ("A4", vec!["A B", "A B"], vec!["1", "'A B'", "2", "'A B'"]),
            (
                "A5",
                vec!["1", "?", "A B"],
                vec!["1", "1", "2", "'?'", "3", "'A B'"],
            ),
        ];
        let blocks: [(&str, &str); 3] = [("F0", F0), ("F1", F1), ("F2", F2)];

        let mut calls = 0usize;
        let mut successes = 0usize;
        let mut empty_preserved = 0usize;
        let mut publications = 0usize;
        let mut rows_total = 0usize;
        let mut values_total = 0usize;
        let mut loops_total = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        for (profile_name, names, expected_raws) in profiles {
            for (state, text) in blocks {
                for repetition in 0..2 {
                    let label = format!("{profile_name}/{state}/{repetition}");
                    let mut block = parse_block(text);
                    let expected_loop_line = verify_block_state(&block, state);
                    // Literal decoded vector bytes/length/order prerequisites.
                    let authors: Vec<String> = names.iter().map(|name| name.to_string()).collect();
                    let decoded_bytes: Vec<&[u8]> = authors.iter().map(String::as_bytes).collect();
                    for (index, name) in names.iter().enumerate() {
                        assert_eq!(decoded_bytes[index], name.as_bytes(), "{label}: bytes");
                    }
                    assert_eq!(authors.len(), names.len(), "{label}: length");
                    let before_authors = authors.clone();
                    let before = block.clone();

                    let result = write_author_category(&authors, &mut block);
                    calls += 1;

                    // Input preservation AFTER the SAME Result, BEFORE output.
                    if authors != before_authors {
                        discrepancies.push(format!("{label}: authors mutated"));
                    }

                    let Some(()) = result.ok() else {
                        discrepancies.push(format!("{label}: unexpected error"));
                        continue;
                    };
                    successes += 1;
                    if names.is_empty() {
                        // Empty guard: WHOLE block equals its captured original.
                        if block != before {
                            discrepancies.push(format!("{label}: empty guard mutated block"));
                        }
                        empty_preserved += 1;
                        continue;
                    }
                    publications += 1;
                    let (rows, values, loops) = collect_output(
                        &block,
                        &before,
                        &expected_raws,
                        expected_loop_line,
                        state,
                        &label,
                        &mut discrepancies,
                    );
                    rows_total += rows;
                    values_total += values;
                    loops_total += loops;
                }
            }
        }

        assert_eq!(calls, 36, "exact 36 real calls");
        assert_eq!(successes, 36, "all 36 success");
        assert_eq!(empty_preserved, 6, "six empty whole-block preservations");
        assert_eq!(publications, 30, "30 active publications");
        assert_eq!(rows_total, 48, "observed row census");
        assert_eq!(values_total, 96, "observed value census");
        assert_eq!(loops_total, 30, "30 actual complete loops");
        assert!(
            discrepancies.is_empty(),
            "product discrepancies: {discrepancies:?}"
        );
    }

    #[test]
    fn bio_author_writer_sequence8() {
        // S0=A4 (2 rows), S1=A0 (whole block retained, SAME 2 rows — no
        // reset), S2=A5 (3 rows, ordinal resets to 1), S3=A1 (1 row,
        // ordinal resets to 1). Emitted censuses count only rows this
        // invocation publishes; observed censuses count POST rows in the
        // block after every call (retained S1 rows included there only).
        let sequences: [(&str, Vec<&str>, Vec<&str>); 4] = [
            ("S0", vec!["A B", "A B"], vec!["1", "'A B'", "2", "'A B'"]),
            ("S1", Vec::new(), Vec::new()),
            (
                "S2",
                vec!["1", "?", "A B"],
                vec!["1", "1", "2", "'?'", "3", "'A B'"],
            ),
            ("S3", vec![""], vec!["1", "''"]),
        ];

        let mut calls = 0usize;
        let mut successes = 0usize;
        let mut publications = 0usize;
        let mut emitted_rows = 0usize;
        let mut emitted_values = 0usize;
        let mut observed_rows = 0usize;
        let mut observed_values = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        for repetition in 0..2 {
            let mut block = parse_block(F1);
            // PRE oracle is ONLY the frozen source-literal raw vector of
            // the last PUBLISHED profile — never actual produced output.
            // After S0 (A4) the literal is A4's ["1","'A B'","2","'A B'"];
            // S1 (A0) preserves it UNCHANGED, so the frame before S2 is
            // still A4's literal (not the empty profile's []); after S2
            // (A5) it becomes A5's literal before S3.
            let mut previous_expected_literal: &[&str] = &[];
            for step in 0..4 {
                let (seq_name, names, expected_raws) = &sequences[step];
                let label = format!("seq/{seq_name}/{step}/{repetition}");
                // Complete literal previous PRE frame BEFORE each call.
                if step == 0 {
                    verify_block_state(&block, "F1");
                } else {
                    assert_eq!(block.name(), "f1", "{label}: PRE name");
                    assert_eq!(block.items().len(), 3, "{label}: PRE items");
                    let CifItem::Pair(entry) = &block.items()[0] else {
                        panic!("{label}: PRE entry");
                    };
                    assert_eq!(
                        (entry.tag(), entry.value().map(|v| v.raw())),
                        ("_entry.id", Some("before"))
                    );
                    assert_eq!(entry.line(), 2);
                    assert_eq!(entry.value().map(|v| (v.line(), v.column())), Some((2, 11)));
                    let CifItem::Loop(previous) = &block.items()[1] else {
                        panic!("{label}: PRE loop");
                    };
                    assert_eq!(previous.tags(), AUTHOR_TAGS, "{label}: PRE tags");
                    if previous.line() != 0 {
                        discrepancies.push(format!("{label}: PRE loop line"));
                    }
                    let actual: Vec<&str> = previous.values().iter().map(|v| v.raw()).collect();
                    if actual != previous_expected_literal {
                        discrepancies.push(format!("{label}: PRE raws {actual:?}"));
                    }
                    for (index, value) in previous.values().iter().enumerate() {
                        if value.line() != 0 || value.column() != 0 {
                            discrepancies.push(format!(
                                "{label}: PRE value {index} {}:{}",
                                value.line(),
                                value.column()
                            ));
                        }
                    }
                    let CifItem::Pair(tail) = &block.items()[2] else {
                        panic!("{label}: PRE tail");
                    };
                    assert_eq!(
                        (tail.tag(), tail.value().map(|v| v.raw())),
                        ("_after.id", Some("tail"))
                    );
                    assert_eq!(tail.line(), 5);
                    assert_eq!(tail.value().map(|v| (v.line(), v.column())), Some((5, 11)));
                }
                // Fresh per-call baselines.
                let authors: Vec<String> = names.iter().map(|name| name.to_string()).collect();
                assert_eq!(authors.len(), names.len(), "{label}: decoded length");
                for (actual, expected) in authors.iter().zip(names) {
                    assert_eq!(
                        actual.as_bytes(),
                        expected.as_bytes(),
                        "{label}: decoded bytes/order"
                    );
                }
                let before_authors = authors.clone();
                let before = block.clone();

                let result = write_author_category(&authors, &mut block);
                calls += 1;

                if authors != before_authors {
                    discrepancies.push(format!("{label}: authors mutated"));
                }
                let Some(()) = result.ok() else {
                    discrepancies.push(format!("{label}: unexpected error"));
                    continue;
                };
                successes += 1;

                if names.is_empty() {
                    // S1: whole block equals pre-call (rows retained).
                    if block != before {
                        discrepancies.push(format!("{label}: empty guard mutated block"));
                    }
                } else {
                    publications += 1;
                    let (rows, values, _) = collect_output(
                        &block,
                        &before,
                        expected_raws,
                        0,
                        if step == 0 { "F1" } else { "F1steady" },
                        &label,
                        &mut discrepancies,
                    );
                    emitted_rows += rows;
                    emitted_values += values;
                }
                // Observed POST census for EVERY call (retained rows count).
                let mut loop_rows = 0usize;
                let mut loop_values = 0usize;
                for item in block.items() {
                    if let CifItem::Loop(row) = item
                        && row
                            .tags()
                            .first()
                            .is_some_and(|tag| tag.starts_with("_audit_author."))
                    {
                        loop_values = row.values().len();
                        loop_rows = loop_values / 2;
                    }
                }
                // The PRE oracle advances ONLY on a published (nonempty)
                // invocation and takes the FROZEN literal of that profile;
                // the empty guard advances nothing (retention).
                if !names.is_empty() {
                    previous_expected_literal = expected_raws;
                }
                observed_rows += loop_rows;
                observed_values += loop_values;
            }
        }

        assert_eq!(calls, 8, "exact 8 sequence calls");
        assert_eq!(successes, 8, "eight successes");
        assert_eq!(publications, 6, "six active publications");
        assert_eq!(emitted_rows, 12, "12 emitted rows");
        assert_eq!(emitted_values, 24, "24 emitted values");
        assert_eq!(observed_rows, 16, "16 observed POST rows");
        assert_eq!(observed_values, 32, "32 observed POST values");
        assert!(
            discrepancies.is_empty(),
            "sequence discrepancies: {discrepancies:?}"
        );
    }

    #[test]
    fn bio_author_writer_mixed4() {
        const MIXED: &str = "data_bad\n_entry.id before\nloop_\n_audit_author.pdbx_ordinal\n_foreign.keep\nOLD stale\n_after.id tail\n";
        let mut calls = 0usize;
        let mut successes = 0usize;
        let mut errors_observed = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        for profile in ["A0", "A1"] {
            for repetition in 0..2 {
                let label = format!("mixed/{profile}/{repetition}");
                let mut block = parse_block(MIXED);
                // Full parsed prerequisites on the exact mixed fixture.
                assert_eq!(block.name(), "bad", "{label}: name");
                assert_eq!(block.items().len(), 3, "{label}: items");
                let CifItem::Pair(entry) = &block.items()[0] else {
                    panic!("{label}: entry");
                };
                assert_eq!(
                    (entry.tag(), entry.value().map(|v| v.raw())),
                    ("_entry.id", Some("before"))
                );
                assert_eq!(entry.line(), 2);
                assert_eq!(entry.value().map(|v| (v.line(), v.column())), Some((2, 11)));
                let CifItem::Loop(stale) = &block.items()[1] else {
                    panic!("{label}: loop position");
                };
                assert_eq!(stale.line(), 3, "{label}: loop line");
                assert_eq!(
                    stale.tags(),
                    ["_audit_author.pdbx_ordinal", "_foreign.keep"],
                    "{label}: loop tags"
                );
                assert_eq!(
                    stale.values().iter().map(|v| v.raw()).collect::<Vec<_>>(),
                    ["OLD", "stale"],
                    "{label}: loop raws"
                );
                assert_eq!(stale.values()[0].line(), 6);
                assert_eq!(stale.values()[0].column(), 1);
                assert_eq!(stale.values()[1].line(), 6);
                assert_eq!(stale.values()[1].column(), 5);
                let CifItem::Pair(tail) = &block.items()[2] else {
                    panic!("{label}: tail");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|v| v.raw())),
                    ("_after.id", Some("tail"))
                );
                assert_eq!(tail.line(), 7);
                assert_eq!(tail.value().map(|v| (v.line(), v.column())), Some((7, 11)));

                let authors: Vec<String> = if profile == "A0" {
                    Vec::new()
                } else {
                    vec![String::new()]
                };
                let expected_names: &[&str] = if profile == "A0" { &[] } else { &[""] };
                assert_eq!(
                    authors.len(),
                    expected_names.len(),
                    "{label}: decoded length"
                );
                for (actual, expected) in authors.iter().zip(expected_names) {
                    assert_eq!(
                        actual.as_bytes(),
                        expected.as_bytes(),
                        "{label}: decoded bytes/order"
                    );
                }
                let before_authors = authors.clone();
                let before = block.clone();

                let result = write_author_category(&authors, &mut block);
                calls += 1;

                // Whole block/input preservation after EVERY call BEFORE
                // any error matching.
                if authors != before_authors {
                    discrepancies.push(format!("{label}: authors mutated"));
                }
                if block != before {
                    discrepancies.push(format!("{label}: block mutated"));
                }

                match result.as_ref() {
                    Ok(()) => {
                        if profile != "A0" {
                            discrepancies.push(format!("{label}: unexpected success"));
                        } else {
                            successes += 1;
                        }
                    }
                    Err(error) => {
                        if profile != "A1" {
                            discrepancies.push(format!("{label}: unexpected error"));
                            continue;
                        }
                        errors_observed += 1;
                        if error.kind() != CifReadErrorKind::InvalidLoop {
                            discrepancies.push(format!("{label}: kind {:?}", error.kind()));
                        }
                        if error.source() != "cif" {
                            discrepancies.push(format!("{label}: source"));
                        }
                        if error.line() != 3 {
                            discrepancies.push(format!("{label}: line {}", error.line()));
                        }
                        if error.column() != 1 {
                            discrepancies.push(format!("{label}: column"));
                        }
                        if error.message() != "Tag _foreign.keep in loop with _audit_author." {
                            discrepancies.push(format!("{label}: message {}", error.message()));
                        }
                        // std::error::Error::source()==None per current
                        // CifReadError — existing sourceNone, not a
                        // source-chain/downcast identity claim.
                        if std::error::Error::source(error).is_some() {
                            discrepancies.push(format!("{label}: unexpected source chain"));
                        }
                    }
                }
            }
        }

        assert_eq!(calls, 4, "exact 4 mixed calls");
        assert_eq!(successes, 2, "two successes (A0 empty guard)");
        assert_eq!(errors_observed, 2, "two actual errors observed");
        assert!(
            discrepancies.is_empty(),
            "mixed discrepancies: {discrepancies:?}"
        );
    }
}
