//! Private active mmCIF `_cell` category emitter (BIO-CELL-CATEGORY1-28).
//!
//! Active-helper only: the coordinate composer is NOT dispatched here and
//! caller group dispatch remains ROOT-owned. Reuses the accepted six-field
//! owner and the one existing category setter; no second implementation.

use super::super::cif::CifBlock;
use cosmolkit_bio::BioCrystalInfo;

pub(super) fn write_cell_category(
    crystal: &BioCrystalInfo,
    entry_id_raw: &str,
    block: &mut CifBlock,
) {
    // BEGIN GEMMI CPP FUNCTION gemmi::write_cell_category active branch (to_mmcif.cpp:490-497)
    // Gemmi❌❌:   if (groups.cell) {
    // Gemmi❗❌:     cif::ItemSpan cell_span(block.items, "_cell.");
    // Gemmi❗❌:     cell_span.set_pair("_cell.entry_id", id);
    // Gemmi❗✔️:     write_cell_parameters(st.cell, cell_span);
    // Gemmi✔️🔝:     auto z_pdb = st.info.find("_cell.Z_PDB");
    // Gemmi✔️✔️:     if (z_pdb != st.info.end())
    // Gemmi❗❌:       cell_span.set_pair(z_pdb->first, z_pdb->second);
    // Gemmi❌❌:   }
    // END GEMMI CPP FUNCTION
    //
    // Guard scope: the outer `if (groups.cell)` line is ❌❌ because the
    // caller group dispatch is OUTSIDE this active-helper boundary — the
    // coordinate-only composer does not call this function and public
    // group/default design remains ROOT-owned. This helper is the
    // complete ACTIVE branch, not a public dispatch stub.
    //
    // Helper bodies consulted (cross-file rule, NOT duplicated): the
    // ItemSpan prefix constructor + set_pair (cifdoc.hpp:613-645) are the
    // existing CifBlock::set_pair_in_category; write_cell_parameters
    // (:296-303) is the accepted sibling owner (its own in-body markers
    // govern its six numeric lines); st.info.find maps to the typed
    // BioCrystalInfo::z_pdb() Option field — the exact `_cell.Z_PDB`
    // source key's existing typed storage, not a second property map.
    //
    // Behavior review: ONE setter for `_cell.entry_id` with
    // entry_id_raw verbatim (already CIF raw text, quotes preserved,
    // no filtering/decoding/trimming/quoting/validation/error); ONE
    // accepted write_cell_parameters(crystal.cell(), block) for the six
    // numeric pairs in source order; IFF z_pdb() is Some(raw), ONE
    // setter for `_cell.Z_PDB` with raw verbatim — Some("") writes the
    // empty raw string, None performs NO third call and does NOT delete
    // an existing parsed Z item.
    //
    // Cost review (separate axes): the source constructs ONE persistent
    // ItemSpan shared by 7-8 set_pair invocations while the reused Rust
    // setter recomputes the category span on EACH invocation — 7 or 8
    // span recomputations versus one construction, a KNOWN local extra
    // scan (❌ on the span/set_pair lines above), not unresolved. Each
    // raw value is a fresh Rust heap String clone where source
    // std::string uses SSO inline storage for short outputs
    // (ABI-dependent, unmeasured, separate caveat). The delegated
    // six-field numeric profile stays exactly as marked inside the
    // accepted owner (❗❌). LOOKUP-SCOPED 🔝 on the find line: the source
    // `st.info` is a std::map<std::string,std::string> (model.hpp:916),
    // so `st.info.find("_cell.Z_PDB")` is an O(log N) key-comparison
    // lookup (string-comparison cost separate); the Rust
    // `crystal.z_pdb()` Option<&str> typed-field read is a direct O(1)
    // field access that does not change semantics (same present/absent
    // result, same raw value) and is strictly better for THIS lookup
    // only — the marker is scoped to the lookup itself, NOT to the
    // composed span/numeric/string costs, which carry their own markers
    // above. No all-float or blanket cost parity claim.
    block.set_pair_in_category(Some("_cell."), "_cell.entry_id", entry_id_raw.to_owned());
    super::cell::write_cell_parameters(crystal.cell(), block);
    if let Some(raw) = crystal.z_pdb() {
        block.set_pair_in_category(Some("_cell."), "_cell.Z_PDB", raw.to_owned());
    }
}

#[cfg(test)]
mod bio_cell_category_tests {
    use super::write_cell_category;
    use crate::cif::{CifBlock, CifCheckLevel, CifItem, read_cif_document};
    use cosmolkit_bio::{BioCrystalCell, BioCrystalInfo, BioTransform};

    const NUMERIC_TAGS: [&str; 6] = [
        "_cell.length_a",
        "_cell.length_b",
        "_cell.length_c",
        "_cell.angle_alpha",
        "_cell.angle_beta",
        "_cell.angle_gamma",
    ];

    const C0_RAW: [&str; 6] = ["1", "1", "1", "90", "90", "90"];
    const C2_RAW: [&str; 6] = ["-0", "0", "-1", "-0", "90", "180"];

    fn c0() -> BioCrystalCell {
        BioCrystalCell {
            a: 1.0,
            b: 1.0,
            c: 1.0,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        }
    }
    fn c2() -> BioCrystalCell {
        BioCrystalCell {
            a: -0.0,
            b: 0.0,
            c: -1.0,
            alpha: -0.0,
            beta: 90.0,
            gamma: 180.0,
        }
    }

    fn crystal_for(cell: BioCrystalCell, z_pdb: Option<&str>) -> BioCrystalInfo {
        BioCrystalInfo::new(
            cell,
            Some("P 1".to_string()),
            z_pdb.map(str::to_string),
            BioTransform::identity(),
            BioTransform::identity(),
            false,
            0,
            Vec::new(),
        )
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

    fn crystal_float_bits(crystal: &BioCrystalInfo) -> Vec<u64> {
        let mut bits = cell_bits(&crystal.cell()).to_vec();
        for transform in [crystal.orthogonal(), crystal.fractional()]
            .into_iter()
            .chain(crystal.symmetry_images())
        {
            bits.extend(
                transform
                    .matrix()
                    .iter()
                    .flatten()
                    .map(|value| value.to_bits()),
            );
            bits.extend(transform.translation().iter().map(|value| value.to_bits()));
        }
        bits.push(crystal.volume().to_bits());
        bits.extend(
            crystal
                .reciprocal_lengths()
                .iter()
                .map(|value| value.to_bits()),
        );
        bits.extend(
            crystal
                .reciprocal_cosines()
                .iter()
                .map(|value| value.to_bits()),
        );
        bits
    }

    /// Per-call constructor-getter prerequisite verification (no invented
    /// property calculation or crystal validation).
    fn verify_crystal(
        crystal: &BioCrystalInfo,
        cell: BioCrystalCell,
        z_pdb: Option<&str>,
        label: &str,
    ) {
        let literal_bits = if cell.a.to_bits() == 0x3ff0000000000000 {
            [
                0x3ff0000000000000,
                0x3ff0000000000000,
                0x3ff0000000000000,
                0x4056800000000000,
                0x4056800000000000,
                0x4056800000000000,
            ]
        } else {
            [
                0x8000000000000000,
                0,
                0xbff0000000000000,
                0x8000000000000000,
                0x4056800000000000,
                0x4066800000000000,
            ]
        };
        assert_eq!(
            cell_bits(&crystal.cell()),
            literal_bits,
            "{label}: frozen cell bits"
        );
        assert_eq!(
            cell_bits(&crystal.cell()),
            cell_bits(&cell),
            "{label}: cell bits"
        );
        assert_eq!(
            crystal.space_group_hm(),
            Some("P 1"),
            "{label}: space group"
        );
        assert_eq!(crystal.z_pdb(), z_pdb, "{label}: z_pdb");
        assert_eq!(
            *crystal.orthogonal().matrix(),
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        );
        assert_eq!(*crystal.orthogonal().translation(), [0.0; 3]);
        assert_eq!(
            *crystal.fractional().matrix(),
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        );
        assert_eq!(*crystal.fractional().translation(), [0.0; 3]);
        assert_eq!(crystal.volume().to_bits(), 1.0f64.to_bits());
        assert_eq!(*crystal.reciprocal_lengths(), [1.0; 3]);
        assert_eq!(*crystal.reciprocal_cosines(), [0.0; 3]);
        assert!(!crystal.explicit_matrices());
        assert_eq!(crystal.cs_count(), 0);
        assert!(crystal.symmetry_images().is_empty());
    }

    fn parse_block(text: &str) -> CifBlock {
        read_cif_document(text, "cell_category_fixture", CifCheckLevel::Syntax)
            .expect("frozen block state parses under the Syntax owner")
            .blocks()[0]
            .clone()
    }

    const F0: &str = "data_f0\n_entry.id before\n_after.id tail\n";
    const F1: &str =
        "data_f1\n_entry.id before\n_cell.entry_id OLD\n_cell.Z_PDB 99\n_after.id tail\n";
    const F2: &str = "data_f2\n_entry.id before\nloop_\n_cell.entry_id\n_foreign.keep\nOLD stale\n_after.id tail\n";

    /// Per-state frozen POST layout for one completed call.
    struct Layout {
        total_items: usize,
        entry_index: usize,
        /// (pre, post) — None when the block has no entry to preserve.
        numeric_start: usize,
        z_index: Option<usize>,
        entry0: (usize, usize),
        tail: (usize, usize),
    }

    fn layout_for(name: &str, z_present_in_output: bool) -> Layout {
        match name {
            "F0" => Layout {
                total_items: if z_present_in_output { 10 } else { 9 },
                entry_index: 2,
                numeric_start: 3,
                z_index: z_present_in_output.then_some(9),
                entry0: (0, 0),
                tail: (1, 1),
            },
            "F1" => Layout {
                total_items: 10,
                entry_index: 1,
                numeric_start: 3,
                z_index: Some(2),
                entry0: (0, 0),
                tail: (3, 9),
            },
            "F2" => Layout {
                total_items: if z_present_in_output { 10 } else { 9 },
                entry_index: 1,
                numeric_start: 2,
                z_index: z_present_in_output.then_some(8),
                entry0: (0, 0),
                tail: (2, if z_present_in_output { 9 } else { 8 }),
            },
            other => panic!("unknown frozen block state {other}"),
        }
    }

    /// Setup-time parsed-frame prerequisite. Returns (F1 entry line, F1 Z
    /// line) when the state is F1 — the preserved parsed Pair lines.
    fn verify_block_state(block: &CifBlock, name: &str) -> (usize, usize) {
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
                let CifItem::Pair(tail) = &block.items()[1] else {
                    panic!("F0 tail position");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|v| v.raw())),
                    ("_after.id", Some("tail"))
                );
                (0, 0)
            }
            "F1" => {
                assert_eq!(block.items().len(), 4, "F1 literal dimension");
                assert_eq!(block.name(), "f1");
                let CifItem::Pair(entry) = &block.items()[1] else {
                    panic!("F1 parsed entry position");
                };
                assert_eq!(
                    (entry.tag(), entry.value().map(|v| v.raw())),
                    ("_cell.entry_id", Some("OLD"))
                );
                let entry_line = entry.line();
                assert_eq!(entry_line, 3, "F1 literal entry Pair line");
                let CifItem::Pair(z) = &block.items()[2] else {
                    panic!("F1 parsed Z position");
                };
                assert_eq!(
                    (z.tag(), z.value().map(|v| v.raw())),
                    ("_cell.Z_PDB", Some("99"))
                );
                let z_line = z.line();
                assert_eq!(z_line, 4, "F1 literal Z Pair line");
                let z_value = z.value().expect("F1 parsed Z value").clone();
                let CifItem::Pair(tail) = &block.items()[3] else {
                    panic!("F1 tail position");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|v| v.raw())),
                    ("_after.id", Some("tail"))
                );
                let _ = z_value;
                (entry_line, z_line)
            }
            "F2" => {
                assert_eq!(block.items().len(), 3, "F2 literal dimension");
                assert_eq!(block.name(), "f2");
                let CifItem::Loop(stale) = &block.items()[1] else {
                    panic!("F2 stale loop position");
                };
                assert_eq!(
                    stale.tags(),
                    ["_cell.entry_id", "_foreign.keep"],
                    "F2 stale tags"
                );
                assert_eq!(
                    stale.values().iter().map(|v| v.raw()).collect::<Vec<_>>(),
                    ["OLD", "stale"],
                    "F2 stale raws"
                );
                let CifItem::Pair(tail) = &block.items()[2] else {
                    panic!("F2 tail position");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|v| v.raw())),
                    ("_after.id", Some("tail"))
                );
                (0, 0)
            }
            other => panic!("unknown frozen block state {other}"),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn collect_output(
        block: &CifBlock,
        before: &CifBlock,
        layout: &Layout,
        entry_raw: &str,
        numeric_raws: [&str; 6],
        z_expected: Option<&str>,
        f1_parsed_lines: (usize, usize),
        state: &str,
        label: &str,
        discrepancies: &mut Vec<String>,
    ) -> (usize, usize, usize) {
        let items = block.items();
        if items.len() != layout.total_items {
            discrepancies.push(format!(
                "{label}: item count {} != {}",
                items.len(),
                layout.total_items
            ));
        }
        // ---- active entry pair ----
        let duplicate_entry = items
            .iter()
            .filter(|item| matches!(item, CifItem::Pair(p) if p.tag().eq_ignore_ascii_case("_cell.entry_id")))
            .count();
        if duplicate_entry != 1 {
            discrepancies.push(format!("{label}: entry duplicates {duplicate_entry}"));
        }
        match items.get(layout.entry_index) {
            Some(CifItem::Pair(pair)) => {
                if pair.tag() != "_cell.entry_id" {
                    discrepancies.push(format!("{label}: entry tag {}", pair.tag()));
                }
                if pair.value().map(|v| v.raw()) != Some(entry_raw) {
                    discrepancies.push(format!(
                        "{label}: entry raw {:?} != {entry_raw}",
                        pair.value().map(|v| v.raw())
                    ));
                }
                // Pair.line: F1 replacement keeps parsed line; F0/F2 fresh 0.
                let expected_entry_line = if state == "F1" { f1_parsed_lines.0 } else { 0 };
                if pair.line() != expected_entry_line {
                    discrepancies.push(format!(
                        "{label}: entry line {} != {expected_entry_line}",
                        pair.line()
                    ));
                }
                if let Some(value) = pair.value() {
                    if value.line() != 0 || value.column() != 0 {
                        discrepancies.push(format!(
                            "{label}: entry value position {}:{}",
                            value.line(),
                            value.column()
                        ));
                    }
                }
            }
            _ => discrepancies.push(format!("{label}: no entry pair at {}", layout.entry_index)),
        }
        // ---- six numeric pairs ----
        for (field, tag) in NUMERIC_TAGS.iter().enumerate() {
            let expected_index = layout.numeric_start + field;
            let duplicate_count = items
                .iter()
                .filter(
                    |item| matches!(item, CifItem::Pair(p) if p.tag().eq_ignore_ascii_case(tag)),
                )
                .count();
            if duplicate_count != 1 {
                discrepancies.push(format!(
                    "{label}: field {field} duplicates {duplicate_count}"
                ));
            }
            match items.get(expected_index) {
                Some(CifItem::Pair(pair)) => {
                    if pair.tag() != *tag {
                        discrepancies.push(format!("{label}: field {field} tag {}", pair.tag()));
                    }
                    if pair.value().map(|v| v.raw()) != Some(numeric_raws[field]) {
                        discrepancies.push(format!(
                            "{label}: field {field} raw {:?} != {}",
                            pair.value().map(|v| v.raw()),
                            numeric_raws[field]
                        ));
                    }
                    if pair.line() != 0 {
                        discrepancies.push(format!("{label}: field {field} line {}", pair.line()));
                    }
                    if let Some(value) = pair.value() {
                        if value.line() != 0 || value.column() != 0 {
                            discrepancies.push(format!(
                                "{label}: field {field} value {}:{}",
                                value.line(),
                                value.column()
                            ));
                        }
                    }
                }
                _ => discrepancies.push(format!(
                    "{label}: field {field} no pair at {expected_index}"
                )),
            }
        }
        // ---- Z pair ----
        let duplicate_z = items
            .iter()
            .filter(|item| matches!(item, CifItem::Pair(p) if p.tag().eq_ignore_ascii_case("_cell.Z_PDB")))
            .count();
        if duplicate_z > 1 {
            discrepancies.push(format!("{label}: Z duplicates {duplicate_z}"));
        }
        match (layout.z_index, z_expected) {
            (Some(index), Some(expected_raw)) => match items.get(index) {
                Some(CifItem::Pair(pair)) => {
                    if pair.tag() != "_cell.Z_PDB" {
                        discrepancies.push(format!("{label}: Z tag {}", pair.tag()));
                    }
                    if pair.value().map(|v| v.raw()) != Some(expected_raw) {
                        discrepancies.push(format!(
                            "{label}: Z raw {:?} != {expected_raw}",
                            pair.value().map(|v| v.raw())
                        ));
                    }
                    // F1 Z replacement keeps parsed line4; F0/F2 fresh 0.
                    let expected_z_line = if state == "F1" { f1_parsed_lines.1 } else { 0 };
                    if pair.line() != expected_z_line {
                        discrepancies.push(format!(
                            "{label}: Z line {} != {expected_z_line}",
                            pair.line()
                        ));
                    }
                    if let Some(value) = pair.value() {
                        if value.line() != 0 || value.column() != 0 {
                            discrepancies.push(format!(
                                "{label}: Z value {}:{}",
                                value.line(),
                                value.column()
                            ));
                        }
                    }
                }
                _ => discrepancies.push(format!("{label}: no Z pair at {index}")),
            },
            (Some(index), None) => {
                // F1 None: parsed "99" preserved WHOLE (line + value metadata).
                if items.get(index) != before.items().get(index) {
                    discrepancies.push(format!("{label}: F1 whole preserved Z changed"));
                }
                match items.get(index) {
                    Some(CifItem::Pair(pair)) => {
                        if pair.tag() != "_cell.Z_PDB"
                            || pair.value().map(|v| v.raw()) != Some("99")
                        {
                            discrepancies.push(format!("{label}: F1 preserved Z changed"));
                        }
                        if pair.line() != f1_parsed_lines.1 {
                            discrepancies.push(format!("{label}: F1 preserved Z line"));
                        }
                    }
                    _ => discrepancies.push(format!("{label}: F1 preserved Z missing at {index}")),
                }
            }
            (None, Some(_)) => discrepancies.push(format!("{label}: unexpected Z present")),
            (None, None) => {
                // Absence verified: no _cell.Z_PDB pair anywhere.
                if duplicate_z != 0 {
                    discrepancies.push(format!("{label}: unexpected Z pair on None arm"));
                }
            }
        }
        // ---- unrelated items + name ----
        if before.items().get(layout.entry0.0) != items.get(layout.entry0.1) {
            discrepancies.push(format!("{label}: unrelated entry0 changed"));
        }
        if before.items().get(layout.tail.0) != items.get(layout.tail.1) {
            discrepancies.push(format!("{label}: unrelated tail changed"));
        }
        if block.name() != before.name() {
            discrepancies.push(format!("{label}: block name changed"));
        }
        // ---- observed censuses ----
        let numeric = items
            .iter()
            .filter(|item| {
                matches!(item, CifItem::Pair(p) if NUMERIC_TAGS
                .iter()
                .any(|tag| p.tag().eq_ignore_ascii_case(tag)))
            })
            .count();
        let entry = items
            .iter()
            .filter(|item| matches!(item, CifItem::Pair(p) if p.tag().eq_ignore_ascii_case("_cell.entry_id")))
            .count();
        let z = duplicate_z;
        (numeric, entry, z)
    }

    #[test]
    fn bio_cell_category_product96() {
        let cells: [(&str, BioCrystalCell, [&str; 6]); 2] =
            [("C0", c0(), C0_RAW), ("C2", c2(), C2_RAW)];
        let z_states: [(&str, Option<&str>); 4] = [
            ("ZNone", None),
            ("ZEmpty", Some("")),
            ("Z4", Some("4")),
            ("Z004", Some("004")),
        ];
        let entries: [&str; 2] = ["ID", "'E 1'"];
        let blocks: [(&str, &str); 3] = [("F0", F0), ("F1", F1), ("F2", F2)];

        let mut calls = 0usize;
        let mut numeric_pairs = 0usize;
        let mut entry_pairs = 0usize;
        let mut z_pairs = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        for (cell_name, cell_value, numeric_raws) in cells {
            for (z_name, z_value) in z_states {
                for entry_raw in entries {
                    for (block_name, block_text) in blocks {
                        for repetition in 0..2 {
                            let label = format!(
                                "{cell_name}/{z_name}/{entry_raw}/{block_name}/{repetition}"
                            );
                            let mut block = parse_block(block_text);
                            let f1_lines = verify_block_state(&block, block_name);
                            let crystal = crystal_for(cell_value, z_value);
                            verify_crystal(&crystal, cell_value, z_value, &label);
                            let before = block.clone();
                            let crystal_before = crystal.clone();
                            let crystal_bits_before = crystal_float_bits(&crystal);
                            let entry_bytes = entry_raw.as_bytes().to_vec();

                            write_cell_category(&crystal, entry_raw, &mut block);
                            calls += 1;

                            // Input preservation BEFORE output collection.
                            if crystal_float_bits(&crystal) != crystal_bits_before {
                                discrepancies.push(format!(
                                    "{label}: accessible crystal float bits mutated"
                                ));
                            }
                            if cell_bits(&crystal.cell()) != cell_bits(&cell_value) {
                                discrepancies.push(format!("{label}: cell bits mutated"));
                            }
                            if crystal != crystal_before {
                                discrepancies.push(format!("{label}: crystal mutated"));
                            }
                            if crystal.z_pdb() != z_value {
                                discrepancies.push(format!("{label}: z_pdb mutated"));
                            }
                            if entry_raw.as_bytes() != entry_bytes {
                                discrepancies.push(format!("{label}: entry bytes mutated"));
                            }

                            let z_in_output = z_value.is_some() || block_name == "F1";
                            let layout = layout_for(block_name, z_in_output);
                            let (n, e, z) = collect_output(
                                &block,
                                &before,
                                &layout,
                                entry_raw,
                                numeric_raws,
                                z_value,
                                f1_lines,
                                block_name,
                                &label,
                                &mut discrepancies,
                            );
                            numeric_pairs += n;
                            entry_pairs += e;
                            z_pairs += z;
                        }
                    }
                }
            }
        }

        assert_eq!(calls, 96, "exact 96 real write_cell_category calls");
        assert_eq!(numeric_pairs, 576, "observed numeric pair census");
        assert_eq!(entry_pairs, 96, "observed entry pair census");
        assert_eq!(z_pairs, 80, "observed Z pair census");
        assert!(
            discrepancies.is_empty(),
            "product discrepancies: {discrepancies:?}"
        );
    }

    #[test]
    fn bio_cell_category_sequence8() {
        // Frozen sequence: each repetition runs four real calls on ONE
        // parsed F1 block. Z expected ["4","4","","004"] — the second
        // call's None RETAINS the previous Z item WHOLE.
        let sequence: [(&str, BioCrystalCell, &str, Option<&str>, [&str; 6], &str); 4] = [
            ("C0", c0(), "ID", Some("4"), C0_RAW, "4"),
            ("C2", c2(), "'E 1'", None, C2_RAW, "4"),
            ("C0", c0(), "ID", Some(""), C0_RAW, ""),
            ("C2", c2(), "'E 1'", Some("004"), C2_RAW, "004"),
        ];
        let layout = Layout {
            total_items: 10,
            entry_index: 1,
            numeric_start: 3,
            z_index: Some(2),
            entry0: (0, 0),
            tail: (3, 9),
        };

        let mut calls = 0usize;
        let mut numeric_fields = 0usize;
        let mut entry_pairs = 0usize;
        let mut z_pairs = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        for repetition in 0..2 {
            let mut block = parse_block(F1);
            let original_frame = block.clone();
            for step in 0..4 {
                let (cell_name, cell_value, entry_raw, z_value, numeric_raws, z_expected) =
                    sequence[step];
                let label = format!("seq/{cell_name}/{step}/{repetition}");
                // Fresh per-call crystal + baselines; complete PRE frame
                // verified before every call.
                let crystal = crystal_for(cell_value, z_value);
                verify_crystal(&crystal, cell_value, z_value, &label);
                let (entry_line, z_line) = if step == 0 {
                    verify_block_state(&block, "F1")
                } else {
                    // Complete frozen PRE frame from the previous output:
                    // tags/raws/positions/Pair+Value metadata/name/unrelated.
                    assert_eq!(block.items().len(), 10, "{label}: PRE frame");
                    assert_eq!(block.name(), "f1", "{label}: PRE name");
                    let CifItem::Pair(entry) = &block.items()[1] else {
                        panic!("{label}: PRE entry");
                    };
                    assert_eq!(entry.tag(), "_cell.entry_id", "{label}: PRE entry tag");
                    let prev_entry_raw = sequence[step - 1].2;
                    assert_eq!(
                        entry.value().map(|v| v.raw()),
                        Some(prev_entry_raw),
                        "{label}: PRE entry raw"
                    );
                    let entry_line = entry.line();
                    assert_eq!(entry_line, 3, "{label}: PRE entry Pair line");
                    assert_eq!(
                        entry.value().map(|v| (v.line(), v.column())),
                        Some((0, 0)),
                        "{label}: PRE entry value position"
                    );
                    let CifItem::Pair(z) = &block.items()[2] else {
                        panic!("{label}: PRE Z");
                    };
                    assert_eq!(z.tag(), "_cell.Z_PDB", "{label}: PRE Z tag");
                    assert_eq!(
                        z.value().map(|v| v.raw()),
                        Some(sequence[step - 1].5),
                        "{label}: PRE Z raw"
                    );
                    let z_line = z.line();
                    assert_eq!(z_line, 4, "{label}: PRE Z Pair line");
                    assert_eq!(
                        z.value().map(|v| (v.line(), v.column())),
                        Some((0, 0)),
                        "{label}: PRE Z value position"
                    );
                    for field in 0..6 {
                        let CifItem::Pair(pair) = &block.items()[3 + field] else {
                            panic!("{label}: PRE numeric field {field}");
                        };
                        assert_eq!(
                            pair.tag(),
                            NUMERIC_TAGS[field],
                            "{label}: PRE numeric tag {field}"
                        );
                        assert_eq!(
                            pair.value().map(|v| v.raw()),
                            Some(sequence[step - 1].4[field]),
                            "{label}: PRE numeric raw {field}"
                        );
                        assert_eq!(pair.line(), 0, "{label}: PRE numeric Pair line {field}");
                        assert_eq!(
                            pair.value().map(|v| (v.line(), v.column())),
                            Some((0, 0)),
                            "{label}: PRE numeric value position {field}"
                        );
                    }
                    assert_eq!(
                        &block.items()[0],
                        &original_frame.items()[0],
                        "{label}: PRE whole unrelated entry"
                    );
                    assert_eq!(
                        &block.items()[9],
                        &original_frame.items()[3],
                        "{label}: PRE whole unrelated tail"
                    );
                    let CifItem::Pair(unrelated_entry) = &block.items()[0] else {
                        panic!("{label}: PRE unrelated entry");
                    };
                    assert_eq!(
                        (
                            unrelated_entry.tag(),
                            unrelated_entry.value().map(|v| v.raw())
                        ),
                        ("_entry.id", Some("before"))
                    );
                    let CifItem::Pair(tail) = &block.items()[9] else {
                        panic!("{label}: PRE tail");
                    };
                    assert_eq!(
                        (tail.tag(), tail.value().map(|v| v.raw())),
                        ("_after.id", Some("tail"))
                    );
                    (entry_line, z_line)
                };
                // First call of a repetition maps the parsed tail (3) to
                // its output slot (9); steady calls already have it at 9.
                let call_layout = if step == 0 {
                    &layout
                } else {
                    &Layout {
                        total_items: 10,
                        entry_index: 1,
                        numeric_start: 3,
                        z_index: Some(2),
                        entry0: (0, 0),
                        tail: (9, 9),
                    }
                };
                let before = block.clone();
                let crystal_before = crystal.clone();
                let crystal_bits_before = crystal_float_bits(&crystal);
                let entry_bytes = entry_raw.as_bytes().to_vec();

                write_cell_category(&crystal, entry_raw, &mut block);
                calls += 1;

                // Same-call input preservation BEFORE output collection.
                if crystal_float_bits(&crystal) != crystal_bits_before {
                    discrepancies.push(format!("{label}: accessible crystal float bits mutated"));
                }
                if z_value.is_none() && block.items().get(2) != before.items().get(2) {
                    discrepancies.push(format!("{label}: whole previous Z changed on None"));
                }
                if cell_bits(&crystal.cell()) != cell_bits(&cell_value) {
                    discrepancies.push(format!("{label}: cell bits mutated"));
                }
                if crystal != crystal_before {
                    discrepancies.push(format!("{label}: crystal mutated"));
                }
                if entry_raw.as_bytes() != entry_bytes {
                    discrepancies.push(format!("{label}: entry bytes mutated"));
                }

                // Z expectation for THIS call: Some(raw) writes raw; None
                // retains the previous output's Z whole (raw above).
                let z_for_output = z_value.or(Some(z_expected));
                let (n, e, z) = collect_output(
                    &block,
                    &before,
                    call_layout,
                    entry_raw,
                    numeric_raws,
                    z_for_output,
                    (entry_line, z_line),
                    "F1",
                    &label,
                    &mut discrepancies,
                );
                if n != 6 {
                    discrepancies.push(format!("{label}: observed {n} numeric fields"));
                }
                numeric_fields += n;
                entry_pairs += e;
                z_pairs += z;
                // No growth: F1 layout fixed at 10 items.
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
        assert_eq!(entry_pairs, 8, "observed entry census");
        assert_eq!(z_pairs, 8, "observed Z census");
        assert!(
            discrepancies.is_empty(),
            "sequence discrepancies: {discrepancies:?}"
        );
    }
}
