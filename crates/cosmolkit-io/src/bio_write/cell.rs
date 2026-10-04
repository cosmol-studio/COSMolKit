//! Private mmCIF `_cell` parameter emitter (BIO-CELL-WRITE1-24).
//!
//! Not enabled in the coordinate-only composer; the public writer still
//! does NOT emit cell parameters. Owns only the six source-backed pairs.

use super::super::cif::{CifBlock, format_cif_f64};
use cosmolkit_bio::BioCrystalCell;

pub(super) fn write_cell_parameters(cell: BioCrystalCell, block: &mut CifBlock) {
    // BEGIN GEMMI CPP FUNCTION gemmi::write_cell_parameters (to_mmcif.cpp:296-303)
    // Gemmi✔️✔️: void write_cell_parameters(const UnitCell& cell, cif::ItemSpan& span) {
    // Gemmi❗❌:   span.set_pair("_cell.length_a",    to_str(cell.a));
    // Gemmi❗❌:   span.set_pair("_cell.length_b",    to_str(cell.b));
    // Gemmi❗❌:   span.set_pair("_cell.length_c",    to_str(cell.c));
    // Gemmi❗❌:   span.set_pair("_cell.angle_alpha", to_str(cell.alpha));
    // Gemmi❗❌:   span.set_pair("_cell.angle_beta",  to_str(cell.beta));
    // Gemmi❗❌:   span.set_pair("_cell.angle_gamma", to_str(cell.gamma));
    // Gemmi✔️✔️: }
    // END GEMMI CPP FUNCTION
    //
    // Helper bodies consulted (cross-file rule, NOT duplicated here):
    // ItemSpan prefix constructor + set_pair (cifdoc.hpp:613-645) are the
    // existing CifBlock::set_pair_in_category; to_str(double) "%.9g"
    // (sprintf.hpp:36-40) is crate::cif::format_cif_f64.
    //
    // Behavior review: exactly six pairs in source order — length_a,
    // length_b, length_c, angle_alpha, angle_beta, angle_gamma — each
    // written through the one existing category setter with the numeric
    // text of the delegated formatter. First case-insensitive pair match
    // is replaced in place (caller tag case, existing pair line kept);
    // a loop containing the tag is replaced by a fresh pair; otherwise
    // the pair is inserted at the category span end. No validation,
    // error, default or missing-value behavior is invented.
    //
    // Cost review (separate axes): control order and the six setter
    // dispatches match the source one-for-one. KNOWN local extra cost
    // (❌, not unresolved): the source computes its category span ONCE
    // and keeps that persistent ItemSpan across all six set_pair calls,
    // while the reused Rust helper recomputes the category span on each
    // of the six calls — six span computations versus ONE in the source,
    // a clear local scanning overhead. Separately, each numeric value is
    // a fresh Rust heap String where source std::string uses SSO inline
    // storage for short outputs (ABI-dependent, unmeasured). Pair placement
    // matches the inspected source; delegated numeric behavior remains qualified.
    // format_cif_f64 profile (❗❌ in cif.rs) governs all float text and
    // these fixed literals upgrade nothing about arbitrary output.
    block.set_pair_in_category(Some("_cell."), "_cell.length_a", format_cif_f64(cell.a));
    block.set_pair_in_category(Some("_cell."), "_cell.length_b", format_cif_f64(cell.b));
    block.set_pair_in_category(Some("_cell."), "_cell.length_c", format_cif_f64(cell.c));
    block.set_pair_in_category(
        Some("_cell."),
        "_cell.angle_alpha",
        format_cif_f64(cell.alpha),
    );
    block.set_pair_in_category(
        Some("_cell."),
        "_cell.angle_beta",
        format_cif_f64(cell.beta),
    );
    block.set_pair_in_category(
        Some("_cell."),
        "_cell.angle_gamma",
        format_cif_f64(cell.gamma),
    );
}

#[cfg(test)]
mod bio_cell_writer_tests {
    use super::write_cell_parameters;
    use crate::cif::{CifBlock, CifCheckLevel, CifItem, read_cif_document};
    use cosmolkit_bio::BioCrystalCell;

    const TAGS: [&str; 6] = [
        "_cell.length_a",
        "_cell.length_b",
        "_cell.length_c",
        "_cell.angle_alpha",
        "_cell.angle_beta",
        "_cell.angle_gamma",
    ];

    fn cell(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> BioCrystalCell {
        BioCrystalCell {
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
        }
    }

    fn cell_bits(cell: &BioCrystalCell) -> [u64; 6] {
        [
            cell.a.to_bits(),
            cell.b.to_bits(),
            cell.c.to_bits(),
            cell.alpha.to_bits(),
            cell.beta.to_bits(),
            cell.gamma.to_bits(),
        ]
    }

    /// ROOT-decoded independent literal bit tables (declaration order
    /// a,b,c,alpha,beta,gamma) — per-call prerequisites asserted inside
    /// the invocation loops against these constants, never against
    /// bits derived from the constructed input value.
    const C0_BITS: [u64; 6] = [
        0x3ff0000000000000,
        0x3ff0000000000000,
        0x3ff0000000000000,
        0x4056800000000000,
        0x4056800000000000,
        0x4056800000000000,
    ];
    const C1_BITS: [u64; 6] = [
        0x3ff4000000000000,
        0x4004000000000000,
        0x400e000000000000,
        0x4046800000000000,
        0x404e000000000000,
        0x405e000000000000,
    ];
    const C2_BITS: [u64; 6] = [
        0x8000000000000000,
        0x0000000000000000,
        0xbff0000000000000,
        0x8000000000000000,
        0x4056800000000000,
        0x4066800000000000,
    ];
    const C3_BITS: [u64; 6] = [
        0x3f1a36e2eb1c432d,
        0x3ee4f8b588e368f1,
        0x419d6f3454000000,
        0x41d26580b4c00000,
        0x3fd5555555555555,
        0xc00c000000000000,
    ];
    const C4_BITS: [u64; 6] = [
        0x7ff0000000000000,
        0xfff0000000000000,
        0x7ff8000000000000,
        0x0000000000000000,
        0x8000000000000000,
        0x4062c00000000000,
    ];
    const C5_BITS: [u64; 6] = [
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
    ];

    fn literal_bits_for(name: &str) -> [u64; 6] {
        match name {
            "C0" => C0_BITS,
            "C1" => C1_BITS,
            "C2" => C2_BITS,
            "C3" => C3_BITS,
            "C4" => C4_BITS,
            "C5" => C5_BITS,
            other => panic!("unknown frozen cell literal {other}"),
        }
    }

    fn parse_block(text: &str) -> CifBlock {
        read_cif_document(text, "cell_writer_fixture", CifCheckLevel::Syntax)
            .expect("frozen block state parses under the Syntax owner")
            .blocks()[0]
            .clone()
    }

    const F0: &str = "data_f0\n_entry.id before\n_after.id tail\n";
    const F1: &str = "data_f1\n_entry.id before\n_cell.entry_id E1\n_after.id tail\n";
    const F2: &str =
        "data_f2\n_entry.id before\n_CELL.LENGTH_A OLD\n_cell.Z_PDB 4\n_after.id tail\n";
    const F3: &str = "data_f3\n_entry.id before\nloop_\n_cell.length_a\n_foreign.keep\nOLD stale\n_after.id tail\n";

    /// Per-state frozen layout after one complete call.
    struct Layout {
        total_items: usize,
        /// Output index of each of the six numeric pairs, in field order.
        field_indices: [usize; 6],
        /// (pre-call index, post-call index) of the entry pair.
        entry: (usize, usize),
        /// (pre-call index, post-call index, tag, raw) of a preserved
        /// unrelated mid-category item.
        mid: Option<(usize, usize, &'static str, &'static str)>,
        /// (pre-call index, post-call index) of the tail pair.
        tail: (usize, usize),
    }

    fn layout_for(name: &str) -> Layout {
        match name {
            "F0" => Layout {
                total_items: 8,
                field_indices: [2, 3, 4, 5, 6, 7],
                entry: (0, 0),
                mid: None,
                tail: (1, 1),
            },
            "F1" => Layout {
                total_items: 9,
                field_indices: [2, 3, 4, 5, 6, 7],
                entry: (0, 0),
                mid: Some((1, 1, "_cell.entry_id", "E1")),
                tail: (2, 8),
            },
            "F2" => Layout {
                total_items: 9,
                // a replaces index 1; Z_PDB stays 2; b..gamma 3..7.
                field_indices: [1, 3, 4, 5, 6, 7],
                entry: (0, 0),
                mid: Some((2, 2, "_cell.Z_PDB", "4")),
                tail: (3, 8),
            },
            "F3" => Layout {
                total_items: 8,
                // Loop replaced by the a Pair at 1; b..gamma 2..6.
                field_indices: [1, 2, 3, 4, 5, 6],
                entry: (0, 0),
                mid: None,
                tail: (2, 7),
            },
            other => panic!("unknown frozen block state {other}"),
        }
    }

    /// Setup-time (pre-call) validation of the parsed block literal:
    /// actual tags/raws/positions/dimensions and the block name.
    /// Returns the existing _CELL.LENGTH_A pair line for F2. A failure
    /// here is an invalid fixture — a real setup failure.
    fn verify_block_state(block: &CifBlock, name: &str) -> Option<usize> {
        assert_eq!(block.name(), &format!("{}", name.to_ascii_lowercase()));
        let CifItem::Pair(entry) = &block.items()[0] else {
            panic!("{name} entry pair position");
        };
        assert_eq!(
            (entry.tag(), entry.value().map(|v| v.raw())),
            ("_entry.id", Some("before"))
        );
        match name {
            "F0" => {
                assert_eq!(block.items().len(), 2, "F0 literal dimension");
                let CifItem::Pair(tail) = &block.items()[1] else {
                    panic!("F0 tail pair position");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|v| v.raw())),
                    ("_after.id", Some("tail"))
                );
                None
            }
            "F1" => {
                assert_eq!(block.items().len(), 3, "F1 literal dimension");
                let CifItem::Pair(entry_id) = &block.items()[1] else {
                    panic!("F1 entry_id pair position");
                };
                assert_eq!(
                    (entry_id.tag(), entry_id.value().map(|v| v.raw())),
                    ("_cell.entry_id", Some("E1"))
                );
                let CifItem::Pair(tail) = &block.items()[2] else {
                    panic!("F1 tail pair position");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|v| v.raw())),
                    ("_after.id", Some("tail"))
                );
                None
            }
            "F2" => {
                assert_eq!(block.items().len(), 4, "F2 literal dimension");
                let CifItem::Pair(old) = &block.items()[1] else {
                    panic!("F2 old uppercase pair position");
                };
                assert_eq!(
                    (old.tag(), old.value().map(|v| v.raw())),
                    ("_CELL.LENGTH_A", Some("OLD"))
                );
                let line = old.line();
                let CifItem::Pair(z) = &block.items()[2] else {
                    panic!("F2 Z_PDB pair position");
                };
                assert_eq!(
                    (z.tag(), z.value().map(|v| v.raw())),
                    ("_cell.Z_PDB", Some("4"))
                );
                let CifItem::Pair(tail) = &block.items()[3] else {
                    panic!("F2 tail pair position");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|v| v.raw())),
                    ("_after.id", Some("tail"))
                );
                Some(line)
            }
            "F3" => {
                assert_eq!(block.items().len(), 3, "F3 literal dimension");
                let CifItem::Loop(stale) = &block.items()[1] else {
                    panic!("F3 stale loop position");
                };
                assert_eq!(
                    stale.tags(),
                    ["_cell.length_a", "_foreign.keep"],
                    "F3 stale tags"
                );
                assert_eq!(
                    stale.values().iter().map(|v| v.raw()).collect::<Vec<_>>(),
                    ["OLD", "stale"],
                    "F3 stale raws"
                );
                let CifItem::Pair(tail) = &block.items()[2] else {
                    panic!("F3 tail pair position");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|v| v.raw())),
                    ("_after.id", Some("tail"))
                );
                None
            }
            other => panic!("unknown frozen block state {other}"),
        }
    }

    /// Non-panicking collected output/preservation discrepancies for one
    /// completed call. Returns the OBSERVED numeric pair count.
    #[allow(clippy::too_many_arguments)]
    fn collect_output(
        block: &CifBlock,
        before: &CifBlock,
        layout: &Layout,
        expected_raws: [&str; 6],
        preserved_pair_line: Option<usize>,
        label: &str,
        discrepancies: &mut Vec<String>,
    ) -> usize {
        let items = block.items();
        if items.len() != layout.total_items {
            discrepancies.push(format!(
                "{label}: item count {} != {}",
                items.len(),
                layout.total_items
            ));
        }
        for (field, tag) in TAGS.iter().enumerate() {
            let expected_index = layout.field_indices[field];
            let actual = items.get(expected_index);
            // Duplicate detection runs for EVERY field regardless of
            // whether the selected index holds the expected Pair; a
            // missing field never skips this independent check.
            let duplicate_count = items
                .iter()
                .filter(|item| {
                    matches!(item, CifItem::Pair(pair) if pair.tag().eq_ignore_ascii_case(tag))
                })
                .count();
            if duplicate_count != 1 {
                discrepancies.push(format!(
                    "{label}: field {field} duplicate count {duplicate_count}"
                ));
            }
            let pair = match actual {
                Some(CifItem::Pair(pair)) => Some(pair),
                _ => None,
            };
            let Some(pair) = pair else {
                discrepancies.push(format!(
                    "{label}: field {field} no pair at {expected_index}"
                ));
                // Only Pair-subject checks are skipped when the Pair is
                // absent; the unrelated/name checks after the loop still
                // run for this call.
                continue;
            };
            if pair.tag() != *tag {
                discrepancies.push(format!(
                    "{label}: field {field} tag {} != {tag}",
                    pair.tag()
                ));
            }
            let value = pair.value();
            if value.map(|v| v.raw()) != Some(expected_raws[field]) {
                discrepancies.push(format!(
                    "{label}: field {field} raw {:?} != {}",
                    value.map(|v| v.raw()),
                    expected_raws[field]
                ));
            }
            if let Some(value) = value {
                if value.line() != 0 || value.column() != 0 {
                    discrepancies.push(format!(
                        "{label}: field {field} value position {}:{}",
                        value.line(),
                        value.column()
                    ));
                }
            }
            // Field 0 may have replaced an existing Pair: F2 keeps the
            // existing pair line; fresh inserts/loop replacement carry 0.
            if field == 0 {
                match preserved_pair_line {
                    Some(expected_line) => {
                        if pair.line() != expected_line {
                            discrepancies.push(format!(
                                "{label}: replaced pair line {} != {expected_line}",
                                pair.line()
                            ));
                        }
                    }
                    None => {
                        if pair.line() != 0 {
                            discrepancies
                                .push(format!("{label}: fresh pair line {} != 0", pair.line()));
                        }
                    }
                }
            } else if pair.line() != 0 {
                discrepancies.push(format!(
                    "{label}: field {field} fresh pair line {} != 0",
                    pair.line()
                ));
            }
        }
        // Block name unchanged by the call (independent of field state).
        if block.name() != before.name() {
            discrepancies.push(format!(
                "{label}: block name {} != {}",
                block.name(),
                before.name()
            ));
        }
        // Whole unrelated items equal the PRE-call clone at frozen positions.
        let check_unrelated =
            |before_index: usize, after_index: usize, discrepancies: &mut Vec<String>| {
                if before.items().get(before_index) != items.get(after_index) {
                    discrepancies.push(format!(
                        "{label}: unrelated item {before_index}->{after_index} changed"
                    ));
                }
            };
        check_unrelated(layout.entry.0, layout.entry.1, discrepancies);
        check_unrelated(layout.tail.0, layout.tail.1, discrepancies);
        if let Some((mid_before, mid_after, mid_tag, mid_raw)) = layout.mid {
            check_unrelated(mid_before, mid_after, discrepancies);
            let mid_ok = items.get(mid_after).is_some_and(|item| {
                matches!(item, CifItem::Pair(pair) if pair.tag() == mid_tag
                    && pair.value().is_some_and(|v| v.raw() == mid_raw))
            });
            if !mid_ok {
                discrepancies.push(format!("{label}: mid category item moved"));
            }
        }
        // OBSERVED numeric pair census from actual matching output items.
        items
            .iter()
            .filter(|item| {
                matches!(item, CifItem::Pair(pair) if TAGS
                    .iter()
                    .any(|tag| pair.tag().eq_ignore_ascii_case(tag)))
            })
            .count()
    }

    #[test]
    fn bio_cell_writer_product48() {
        let cells: [(&str, BioCrystalCell, [&str; 6]); 6] = [
            (
                "C0",
                cell(1.0, 1.0, 1.0, 90.0, 90.0, 90.0),
                ["1", "1", "1", "90", "90", "90"],
            ),
            (
                "C1",
                cell(1.25, 2.5, 3.75, 45.0, 60.0, 120.0),
                ["1.25", "2.5", "3.75", "45", "60", "120"],
            ),
            (
                "C2",
                cell(-0.0, 0.0, -1.0, -0.0, 90.0, 180.0),
                ["-0", "0", "-1", "-0", "90", "180"],
            ),
            (
                "C3",
                cell(0.0001, 0.00001, 123456789.0, 1234567891.0, 1.0 / 3.0, -3.5),
                [
                    "0.0001",
                    "1e-05",
                    "123456789",
                    "1.23456789e+09",
                    "0.333333333",
                    "-3.5",
                ],
            ),
            (
                "C4",
                cell(f64::INFINITY, f64::NEG_INFINITY, f64::NAN, 0.0, -0.0, 150.0),
                ["Inf", "-Inf", "NaN", "0", "-0", "150"],
            ),
            (
                "C5",
                cell(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                ["0", "0", "0", "0", "0", "0"],
            ),
        ];
        let block_states: [(&str, &str); 4] = [("F0", F0), ("F1", F1), ("F2", F2), ("F3", F3)];

        let mut calls = 0usize;
        let mut numeric_pairs = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        for (cell_name, cell_value, expected_raws) in cells {
            for (block_name, block_text) in block_states {
                let layout = layout_for(block_name);
                for repetition in 0..2 {
                    let label = format!("{cell_name}/{block_name}/{repetition}");
                    let mut block = parse_block(block_text);
                    // Per-call prerequisite: actual parsed literal, name
                    // and (F2) the existing pair line to preserve.
                    let preserved_pair_line = verify_block_state(&block, block_name);
                    // Per-call independent six-bit prerequisite asserted
                    // INSIDE the invocation loop against the ROOT-decoded
                    // literal table (never outer-loop derived bits).
                    let literal_bits = literal_bits_for(cell_name);
                    let input_bits = cell_bits(&cell_value);
                    assert_eq!(
                        input_bits, literal_bits,
                        "{label}: input bits match literal table"
                    );
                    let before = block.clone();

                    write_cell_parameters(cell_value, &mut block);
                    calls += 1;

                    // Input preservation: the six literal bits unchanged.
                    if cell_bits(&cell_value) != input_bits {
                        discrepancies.push(format!("{label}: input bits mutated"));
                    }
                    // Baseline name retained (output-side independent).
                    if block.name() != before.name() {
                        discrepancies.push(format!("{label}: block name changed"));
                    }
                    // Non-panicking output/preservation collection.
                    numeric_pairs += collect_output(
                        &block,
                        &before,
                        &layout,
                        expected_raws,
                        preserved_pair_line,
                        &label,
                        &mut discrepancies,
                    );
                }
            }
        }

        assert_eq!(calls, 48, "exact 48 real write_cell_parameters calls");
        assert_eq!(numeric_pairs, 288, "observed numeric pair census");
        assert!(
            discrepancies.is_empty(),
            "product discrepancies: {discrepancies:?}"
        );
    }

    #[test]
    fn bio_cell_writer_sequence8() {
        let sequence: [(&str, BioCrystalCell, [&str; 6]); 4] = [
            (
                "C0",
                cell(1.0, 1.0, 1.0, 90.0, 90.0, 90.0),
                ["1", "1", "1", "90", "90", "90"],
            ),
            (
                "C1",
                cell(1.25, 2.5, 3.75, 45.0, 60.0, 120.0),
                ["1.25", "2.5", "3.75", "45", "60", "120"],
            ),
            (
                "C2",
                cell(-0.0, 0.0, -1.0, -0.0, 90.0, 180.0),
                ["-0", "0", "-1", "-0", "90", "180"],
            ),
            (
                "C5",
                cell(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                ["0", "0", "0", "0", "0", "0"],
            ),
        ];
        let layout = layout_for("F1");
        // Steady-state layout for calls after the first in a repetition:
        // the six pairs already exist at 2..7 and the tail sits at 8.
        let steady = Layout {
            total_items: 9,
            field_indices: [2, 3, 4, 5, 6, 7],
            entry: (0, 0),
            mid: Some((1, 1, "_cell.entry_id", "E1")),
            tail: (8, 8),
        };

        let mut calls = 0usize;
        let mut numeric_fields = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        for repetition in 0..2 {
            let mut block = parse_block(F1);
            // Previous-call frozen output literals for the sequence PRE
            // frame: after C0 the frame carries C0's raws, after C1
            // C1's, after C2 C2's (C5 ends each repetition).
            let previous_raws: [&[&str; 6]; 4] = [
                &["1", "1", "1", "90", "90", "90"],
                &["1.25", "2.5", "3.75", "45", "60", "120"],
                &["-0", "0", "-1", "-0", "90", "180"],
                &["0", "0", "0", "0", "0", "0"],
            ];
            for (cell_name, cell_value, expected_raws) in sequence {
                let label = format!("seq/{cell_name}/{repetition}");
                // Fresh per-call baseline and prerequisites. Before the
                // first call of a repetition the frame is the parsed F1
                // literal; afterwards the COMPLETE frozen previous
                // output is verified — numeric Pair raw values AND
                // retained Pair metadata — plus the whole unrelated
                // frame (entry/entry_id/tail positions and values) and
                // the block name.
                if calls % 4 == 0 {
                    assert_eq!(verify_block_state(&block, "F1"), None);
                } else {
                    // Complete frozen PRE frame: previous call's raw
                    // literals at 2..7 with retained metadata.
                    assert_eq!(block.items().len(), 9, "{label}: frame");
                    let CifItem::Pair(entry_id) = &block.items()[1] else {
                        panic!("{label}: frame entry_id");
                    };
                    assert_eq!(
                        (entry_id.tag(), entry_id.value().map(|v| v.raw())),
                        ("_cell.entry_id", Some("E1"))
                    );
                    for (offset, tag) in TAGS.iter().enumerate() {
                        let index = 2 + offset;
                        let CifItem::Pair(pair) = &block.items()[index] else {
                            panic!("{label}: frame pair {tag}");
                        };
                        assert_eq!(pair.tag(), *tag, "{label}: frame tag {tag}");
                        assert_eq!(pair.line(), 0, "{label}: frame pair line {tag}");
                        let previous = previous_raws[(calls % 4) - 1];
                        assert_eq!(
                            pair.value().map(|v| v.raw()),
                            Some(previous[offset]),
                            "{label}: frame raw {tag}"
                        );
                        // Replacement pairs retain their producing
                        // metadata: values keep line 0 / column 0.
                        let value = pair.value().expect("{label}: frame value");
                        assert_eq!(
                            (value.line(), value.column()),
                            (0, 0),
                            "{label}: frame value position {tag}"
                        );
                    }
                    let CifItem::Pair(entry) = &block.items()[0] else {
                        panic!("{label}: frame entry");
                    };
                    assert_eq!(
                        (entry.tag(), entry.value().map(|v| v.raw())),
                        ("_entry.id", Some("before"))
                    );
                    let CifItem::Pair(tail) = &block.items()[8] else {
                        panic!("{label}: frame tail");
                    };
                    assert_eq!(
                        (tail.tag(), tail.value().map(|v| v.raw())),
                        ("_after.id", Some("tail"))
                    );
                    assert_eq!(block.name(), "f1", "{label}: frame name");
                }
                // Per-call independent six-bit prerequisite against the
                // ROOT-decoded literal table.
                let input_bits = cell_bits(&cell_value);
                assert_eq!(
                    input_bits,
                    literal_bits_for(cell_name),
                    "{label}: input bits match literal table"
                );
                let call_layout = if calls % 4 == 0 { &layout } else { &steady };
                let before = block.clone();

                write_cell_parameters(cell_value, &mut block);
                calls += 1;

                if cell_bits(&cell_value) != input_bits {
                    discrepancies.push(format!("{label}: input bits mutated"));
                }
                let observed = collect_output(
                    &block,
                    &before,
                    call_layout,
                    expected_raws,
                    None,
                    &label,
                    &mut discrepancies,
                );
                if observed != 6 {
                    discrepancies.push(format!("{label}: observed {observed} numeric fields"));
                }
                numeric_fields += observed;
                // No duplicate growth: the F1 layout stays fixed.
                if block.items().len() != layout.total_items {
                    discrepancies.push(format!(
                        "{label}: item growth {} != {}",
                        block.items().len(),
                        layout.total_items
                    ));
                }
            }
        }

        assert_eq!(calls, 8, "exact 8 sequence calls");
        assert_eq!(numeric_fields, 48, "observed numeric field census");
        assert!(
            discrepancies.is_empty(),
            "sequence discrepancies: {discrepancies:?}"
        );
    }
}
