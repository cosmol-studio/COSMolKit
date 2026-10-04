//! Private active mmCIF `_entity` category emitter (BIO-ENTITY-WRITE1-32).
//!
//! Active-helper only: the coordinate composer is NOT dispatched here and
//! caller group dispatch remains ROOT-owned. Reuses the one existing loop
//! initialization/row-append chain and the accepted qchain quoting owner.

use super::super::cif::CifBlock;
use cosmolkit_bio::{BioEntityRow, EntityKind};

use super::super::cif::CifReadError;

pub(super) fn write_entity_category(
    entities: &[BioEntityRow],
    block: &mut CifBlock,
) -> Result<(), CifReadError> {
    // BEGIN GEMMI CPP FUNCTION gemmi::write_entity_category active branch (to_mmcif.cpp:508-513)
    // Gemmi❌❌:   if (groups.entity) {
    // Gemmi❗❗:     cif::Loop& entity_loop = block.init_mmcif_loop("_entity.", {"id", "type"});
    // Gemmi✔️✔️:     for (const Entity& ent : st.entities)
    // Gemmi❗❌:       entity_loop.add_row({qchain(ent.name),
    // Gemmi❗❌:                            entity_type_to_string(ent.entity_type)});
    // Gemmi❌❌:   }
    // END GEMMI CPP FUNCTION
    //
    // Guard scope: the outer `if (groups.entity)` line and closing `}` are
    // ❌❌ because the caller group dispatch is OUTSIDE this private
    // active-helper boundary — the coordinate-only composer does not call
    // this function and public group/default design remains ROOT-owned.
    // This helper is the complete ACTIVE branch, not a dispatch stub.
    //
    // Helper bodies consulted (cross-file rule, NOT duplicated): qchain
    // (to_mmcif.cpp:57-59, `cif::quote(s)`) is value.rs::qchain which
    // delegates to quote_cif_value (empty becomes '' not ./?);
    // entity_type_to_string is reproduced below in entity_type_raw;
    // init_mmcif_loop/add_row are the existing CIF owners whose exact
    // InvalidLoop errors propagate unchanged.
    //
    // Behavior review: ONE init_mmcif_loop("_entity.", ["id","type"])
    // EVEN when entities is empty — the category reinitializes, so an
    // empty input still produces the SAME two-tag zero-row loop and any
    // stale parsed category is replaced. Rows iterate the SUPPLIED order
    // with no filtering/sorting/deduplication/re-resolution; each row is
    // ONE add_row whose first value quotes the source entity id EXACTLY
    // ONCE via qchain and whose second value is the bare type literal
    // (never quoted; Unknown is the bare `?`). Source Entity.name maps
    // to the existing source_entity_id storage, not a local row index.
    // No invented validation, defaults or error mapping.
    //
    // Cost review (separate axes): ONE delegated init scan plus ONE walk
    // of the supplied slice — same shape as the source. Per row the Rust
    // path allocates a temporary Vec plus two owned Strings (the quoted
    // id and the type literal) where the source passes an
    // initializer_list of std::string values. The temporary Vec is an
    // extra heap allocation; the short type String also lacks source SSO.
    // These known local costs are ❌, not unresolved ❗ or whole-cost parity.
    // Quoted-id allocation depends on length/source SSO. The type
    // match in entity_type_raw is a fixed O(1) lookup. No blanket
    // allocation-free, whole-serialization, lossless, public/group or
    // default claim.
    let entity_loop = block.init_mmcif_loop("_entity.", &["id", "type"])?;
    for entity in entities {
        entity_loop.add_row(vec![
            super::value::qchain(entity.source().source_entity_id()),
            entity_type_raw(entity.kind()).to_owned(),
        ])?;
    }
    Ok(())
}

/// Exact source `entity_type_to_string` mapping (enumstr.hpp:14-22).
fn entity_type_raw(kind: EntityKind) -> &'static str {
    // BEGIN GEMMI CPP FUNCTION gemmi::entity_type_to_string (enumstr.hpp:14-22)
    // Gemmi✔️✔️: inline const char* entity_type_to_string(EntityType entity_type) {
    // Gemmi✔️✔️:   switch (entity_type) {
    // Gemmi✔️✔️:     case EntityType::Polymer: return "polymer";
    // Gemmi✔️✔️:     case EntityType::Branched: return "branched";
    // Gemmi✔️✔️:     case EntityType::NonPolymer: return "non-polymer";
    // Gemmi✔️✔️:     case EntityType::Water: return "water";
    // Gemmi✔️✔️:     default /*EntityType::Unknown*/: return "?";
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: }
    // END GEMMI CPP FUNCTION
    //
    // Behavior review: the four named kinds map to their exact source
    // literals; Unknown maps to the bare null type `?` exactly as the
    // source default arm. No other variants exist in the Rust enum.
    // Cost review: a fixed O(1) match, identical to the source switch.
    match kind {
        EntityKind::Polymer => "polymer",
        EntityKind::Branched => "branched",
        EntityKind::NonPolymer => "non-polymer",
        EntityKind::Water => "water",
        EntityKind::Unknown => "?",
    }
}

#[cfg(test)]
mod bio_entity_writer_tests {
    use super::write_entity_category;
    use crate::cif::{CifBlock, CifCheckLevel, CifItem, CifReadErrorKind, read_cif_document};
    use cosmolkit_bio::{BioEntityRow, EntityKind, EntitySourceIds, PolymerKind};

    const ENTITY_TAGS: [&str; 2] = ["_entity.id", "_entity.type"];

    fn entity_row(kind: EntityKind, name: &str) -> BioEntityRow {
        BioEntityRow::new(
            kind,
            PolymerKind::Unknown,
            true,
            vec!["ALA".to_string(), "GLY,VAL".to_string()],
            Vec::new(),
            vec!["P12345".to_string()],
            vec!["A".to_string(), "A".to_string()],
            EntitySourceIds::new(name.to_string()),
        )
    }

    fn kind_raw(kind: EntityKind) -> &'static str {
        match kind {
            EntityKind::Polymer => "polymer",
            EntityKind::Branched => "branched",
            EntityKind::NonPolymer => "non-polymer",
            EntityKind::Water => "water",
            EntityKind::Unknown => "?",
        }
    }

    /// Per-call literal prerequisite verification of ALL getters against
    /// the frozen construction. EntityRow contains no floats.
    fn verify_entity(row: &BioEntityRow, kind: EntityKind, name: &str, label: &str) {
        assert_eq!(row.kind(), kind, "{label}: kind");
        assert_eq!(
            row.polymer_kind(),
            PolymerKind::Unknown,
            "{label}: polymer kind"
        );
        assert!(row.reflects_microhetero(), "{label}: microhetero");
        assert_eq!(row.full_sequence(), ["ALA", "GLY,VAL"], "{label}: sequence");
        assert!(row.dbrefs().is_empty(), "{label}: dbrefs");
        assert_eq!(row.sifts_unp_accessions(), ["P12345"], "{label}: sifts");
        assert_eq!(row.subchains(), ["A", "A"], "{label}: subchains");
        assert_eq!(row.source().source_entity_id(), name, "{label}: source id");
    }

    fn parse_block(text: &str) -> CifBlock {
        read_cif_document(text, "entity_writer_fixture", CifCheckLevel::Syntax)
            .expect("frozen block state parses under the Syntax owner")
            .blocks()[0]
            .clone()
    }

    const F0: &str = "data_f0\n_entry.id before\n_after.id tail\n";
    const F1: &str =
        "data_f1\n_entry.id before\n_entity.id OLD\n_entity.type old\n_after.id tail\n";
    const F2: &str =
        "data_f2\n_entry.id before\nloop_\n_entity.id\n_entity.type\nOLD old\n_after.id tail\n";

    fn pair_frame(block: &CifBlock, index: usize) -> Option<(&str, usize, &str, usize, usize)> {
        let CifItem::Pair(pair) = block.items().get(index)? else {
            return None;
        };
        let value = pair.value()?;
        Some((
            pair.tag(),
            pair.line(),
            value.raw(),
            value.line(),
            value.column(),
        ))
    }

    /// Setup-time parsed-frame prerequisite per state. Returns the parsed
    /// Loop.line for F2 (retained by the reinitialized loop).
    fn verify_block_state(block: &CifBlock, name: &str) -> usize {
        match name {
            "F0" => {
                assert_eq!(block.items().len(), 2, "F0 literal dimension");
                assert_eq!(block.name(), "f0");
                assert_eq!(
                    pair_frame(block, 0),
                    Some(("_entry.id", 2, "before", 2, 11))
                );
                assert_eq!(pair_frame(block, 1), Some(("_after.id", 3, "tail", 3, 11)));
                let CifItem::Pair(entry) = &block.items()[0] else {
                    panic!("F0 entry");
                };
                assert_eq!(
                    (entry.tag(), entry.value().map(|v| v.raw())),
                    ("_entry.id", Some("before"))
                );
                let CifItem::Pair(tail) = &block.items()[1] else {
                    panic!("F0 tail");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|v| v.raw())),
                    ("_after.id", Some("tail"))
                );
                0
            }
            "F1" => {
                assert_eq!(block.items().len(), 4, "F1 literal dimension");
                assert_eq!(block.name(), "f1");
                assert_eq!(
                    pair_frame(block, 0),
                    Some(("_entry.id", 2, "before", 2, 11))
                );
                assert_eq!(pair_frame(block, 1), Some(("_entity.id", 3, "OLD", 3, 12)));
                assert_eq!(
                    pair_frame(block, 2),
                    Some(("_entity.type", 4, "old", 4, 14))
                );
                assert_eq!(pair_frame(block, 3), Some(("_after.id", 5, "tail", 5, 11)));
                let CifItem::Pair(old_id) = &block.items()[1] else {
                    panic!("F1 old id");
                };
                assert_eq!(
                    (old_id.tag(), old_id.value().map(|v| v.raw())),
                    ("_entity.id", Some("OLD"))
                );
                let CifItem::Pair(old_type) = &block.items()[2] else {
                    panic!("F1 old type");
                };
                assert_eq!(
                    (old_type.tag(), old_type.value().map(|v| v.raw())),
                    ("_entity.type", Some("old"))
                );
                let CifItem::Pair(tail) = &block.items()[3] else {
                    panic!("F1 tail");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|v| v.raw())),
                    ("_after.id", Some("tail"))
                );
                0
            }
            "F2" => {
                assert_eq!(block.items().len(), 3, "F2 literal dimension");
                assert_eq!(block.name(), "f2");
                assert_eq!(
                    pair_frame(block, 0),
                    Some(("_entry.id", 2, "before", 2, 11))
                );
                assert_eq!(pair_frame(block, 2), Some(("_after.id", 7, "tail", 7, 11)));
                let CifItem::Loop(stale) = &block.items()[1] else {
                    panic!("F2 stale loop");
                };
                assert_eq!(stale.tags(), ["_entity.id", "_entity.type"], "F2 tags");
                assert_eq!(
                    stale.values().iter().map(|v| v.raw()).collect::<Vec<_>>(),
                    ["OLD", "old"],
                    "F2 raws"
                );
                assert_eq!(stale.line(), 3, "F2 literal loop line");
                assert_eq!(
                    stale
                        .values()
                        .iter()
                        .map(|v| (v.line(), v.column()))
                        .collect::<Vec<_>>(),
                    [(6, 1), (6, 5)],
                    "F2 literal value positions"
                );
                let CifItem::Pair(tail) = &block.items()[2] else {
                    panic!("F2 tail");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|v| v.raw())),
                    ("_after.id", Some("tail"))
                );
                3
            }
            other => panic!("unknown frozen state {other}"),
        }
    }

    /// Non-panicking output collection for one successful call. Returns
    /// (complete rows, values, actual category loops) observed counts.
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
        let expected_total = 3;
        if items.len() != expected_total {
            discrepancies.push(format!(
                "{label}: item count {} != {expected_total}",
                items.len()
            ));
        }
        // Find the single _entity loop (duplicate category check).
        let loop_positions: Vec<usize> = items
            .iter()
            .enumerate()
            .filter(|(_, item)| {
                matches!(item, CifItem::Loop(row) if row
                    .tags()
                    .first()
                    .is_some_and(|tag| tag.starts_with("_entity.")))
            })
            .map(|(index, _)| index)
            .collect();
        if loop_positions.len() != 1 {
            discrepancies.push(format!(
                "{label}: entity loop count {}",
                loop_positions.len()
            ));
        }
        // Independent checks must execute even when the loop is missing.
        for tag in ["_entity.id", "_entity.type"] {
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
        let (entry_before, entry_after, tail_before, tail_after) = match state {
            "F0" => (0, 0, 1, 1),
            "F1" => (0, 0, 3, 2),
            _ => (0, 0, 2, 2),
        };
        if before.items().get(entry_before) != items.get(entry_after) {
            discrepancies.push(format!("{label}: unrelated entry changed"));
        }
        if before.items().get(tail_before) != items.get(tail_after) {
            discrepancies.push(format!("{label}: unrelated tail changed"));
        }
        if block.name() != before.name() {
            discrepancies.push(format!("{label}: block name changed"));
        }
        let Some(&loop_index) = loop_positions.first() else {
            return (0, 0, 0);
        };
        let expected_index = match state {
            "F0" => 2,
            _ => 1,
        };
        if loop_index != expected_index {
            discrepancies.push(format!(
                "{label}: loop index {loop_index} != {expected_index}"
            ));
        }
        let Some(CifItem::Loop(row)) = items.get(loop_index) else {
            discrepancies.push(format!("{label}: loop item unavailable"));
            return (0, 0, loop_positions.len());
        };
        if row.tags() != ENTITY_TAGS {
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
        let complete = row.tags().len() == 2 && row.values().len() % 2 == 0;
        if !complete {
            discrepancies.push(format!("{label}: incomplete two-column row frame"));
        }
        let rows = if complete { row.values().len() / 2 } else { 0 };
        (rows, row.values().len(), loop_positions.len())
    }

    #[test]
    fn bio_entity_writer_product120() {
        let kinds = [
            EntityKind::Polymer,
            EntityKind::Branched,
            EntityKind::NonPolymer,
            EntityKind::Water,
            EntityKind::Unknown,
        ];
        // Independent source-derived literals, never the tested quote helper.
        let names = [("1", "1"), ("", "''"), ("A B", "'A B'"), ("?", "'?'")];
        let blocks: [(&str, &str); 3] = [("F0", F0), ("F1", F1), ("F2", F2)];

        let mut calls = 0usize;
        let mut rows_total = 0usize;
        let mut values_total = 0usize;
        let mut loops_total = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        for kind in kinds {
            for (name, raw_name) in names {
                for (state, text) in blocks {
                    for repetition in 0..2 {
                        let label = format!("{}/{name:?}/{state}/{repetition}", kind_raw(kind));
                        let mut block = parse_block(text);
                        let expected_loop_line = verify_block_state(&block, state);
                        let row = entity_row(kind, name);
                        verify_entity(&row, kind, name, &label);
                        let entities = vec![row];
                        let before_entities = entities.clone();
                        let before = block.clone();

                        let result = write_entity_category(&entities, &mut block);
                        calls += 1;

                        // Same-call input preservation BEFORE output checks.
                        if entities != before_entities {
                            discrepancies.push(format!("{label}: entities mutated"));
                        }

                        let Some(()) = result.ok() else {
                            discrepancies.push(format!("{label}: unexpected error"));
                            continue;
                        };
                        let expected_raws = [raw_name, kind_raw(kind)];
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
        }

        assert_eq!(calls, 120, "exact 120 real write_entity_category calls");
        assert_eq!(rows_total, 120, "observed row census");
        assert_eq!(values_total, 240, "observed value census");
        assert_eq!(loops_total, 120, "observed loop census");
        assert!(
            discrepancies.is_empty(),
            "product discrepancies: {discrepancies:?}"
        );
    }

    #[test]
    fn bio_entity_writer_empty_sequence14() {
        let mut calls = 0usize;
        let mut rows_total = 0usize;
        let mut values_total = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        // Part 1: empty input still produces the SAME two-tag zero-row loop.
        for (state, text) in [("F0", F0), ("F1", F1), ("F2", F2)] {
            for repetition in 0..2 {
                let label = format!("empty/{state}/{repetition}");
                let mut block = parse_block(text);
                let expected_loop_line = verify_block_state(&block, state);
                let before = block.clone();
                let entities: Vec<BioEntityRow> = Vec::new();
                let before_entities = entities.clone();

                let result = write_entity_category(&entities, &mut block);
                calls += 1;

                if entities != before_entities {
                    discrepancies.push(format!("{label}: empty entities mutated"));
                }

                let Some(()) = result.ok() else {
                    discrepancies.push(format!("{label}: unexpected error"));
                    continue;
                };
                let (rows, values, _) = collect_output(
                    &block,
                    &before,
                    &[],
                    expected_loop_line,
                    state,
                    &label,
                    &mut discrepancies,
                );
                if rows != 0 || values != 0 {
                    discrepancies.push(format!("{label}: empty input produced rows"));
                }
            }
        }

        // Part 2: four ordered sequence calls on ONE parsed F1 per repetition.
        let sequences: [(&str, Vec<(&str, EntityKind)>, Vec<&str>); 4] = [
            ("S0", Vec::new(), Vec::new()),
            (
                "S1",
                vec![("A B", EntityKind::Polymer), ("", EntityKind::Unknown)],
                vec!["'A B'", "polymer", "''", "?"],
            ),
            ("S2", vec![("1", EntityKind::Water)], vec!["1", "water"]),
            (
                "S3",
                vec![
                    ("1", EntityKind::Water),
                    ("1", EntityKind::NonPolymer),
                    ("?", EntityKind::Branched),
                ],
                vec!["1", "water", "1", "non-polymer", "'?'", "branched"],
            ),
        ];
        for repetition in 0..2 {
            let mut block = parse_block(F1);
            for step in 0..4 {
                let (seq_name, spec, expected_raws) = &sequences[step];
                let label = format!("seq/{seq_name}/{step}/{repetition}");
                // Complete frozen PRE frame before each call.
                if step == 0 {
                    verify_block_state(&block, "F1");
                } else {
                    if block.items().len() != 3 || block.name() != "f1" {
                        discrepancies.push(format!("{label}: PRE frame/name"));
                    }
                    if pair_frame(&block, 0) != Some(("_entry.id", 2, "before", 2, 11)) {
                        discrepancies.push(format!("{label}: PRE full entry"));
                    }
                    if pair_frame(&block, 2) != Some(("_after.id", 5, "tail", 5, 11)) {
                        discrepancies.push(format!("{label}: PRE full tail"));
                    }
                    if let Some(CifItem::Loop(row)) = block.items().get(1) {
                        if row.tags() != ENTITY_TAGS || row.line() != 0 {
                            discrepancies.push(format!("{label}: PRE tags/line"));
                        }
                        let (_, _, previous_raws) = &sequences[step - 1];
                        if row.values().iter().map(|v| v.raw()).collect::<Vec<_>>()
                            != *previous_raws
                        {
                            discrepancies.push(format!("{label}: PRE raws"));
                        }
                        if row
                            .values()
                            .iter()
                            .any(|v| v.line() != 0 || v.column() != 0)
                        {
                            discrepancies.push(format!("{label}: PRE value positions"));
                        }
                    } else {
                        discrepancies.push(format!("{label}: PRE missing loop"));
                    }
                }
                // Fresh per-call entities + baselines.
                let entities: Vec<BioEntityRow> = spec
                    .iter()
                    .map(|(name, kind)| entity_row(*kind, name))
                    .collect();
                for (index, (name, kind)) in spec.iter().enumerate() {
                    verify_entity(&entities[index], *kind, name, &label);
                }
                let before_entities = entities.clone();
                let before = block.clone();

                let result = write_entity_category(&entities, &mut block);
                calls += 1;

                if entities != before_entities {
                    discrepancies.push(format!("{label}: entities mutated"));
                }
                let Some(()) = result.ok() else {
                    discrepancies.push(format!("{label}: unexpected error"));
                    continue;
                };
                let (rows, values, _) = collect_output(
                    &block,
                    &before,
                    expected_raws,
                    0,
                    if step == 0 { "F1" } else { "F1steady" },
                    &label,
                    &mut discrepancies,
                );
                if rows != spec.len() {
                    discrepancies.push(format!("{label}: rows {rows} != {}", spec.len()));
                }
                rows_total += rows;
                values_total += values;
            }
        }

        assert_eq!(calls, 14, "exact 14 combined empty+sequence calls");
        assert_eq!(rows_total, 12, "observed row census (sequence only)");
        assert_eq!(values_total, 24, "observed value census");
        assert!(
            discrepancies.is_empty(),
            "empty/sequence discrepancies: {discrepancies:?}"
        );
    }

    #[test]
    fn bio_entity_writer_mixed_error4() {
        const MIXED: &str = "data_bad\n_entry.id before\nloop_\n_entity.id\n_foreign.keep\nOLD stale\n_after.id tail\n";
        let mut calls = 0usize;
        let mut errors_observed = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        for rows_state in [Vec::new(), vec![("1", EntityKind::Polymer)]] {
            for repetition in 0..2 {
                let label = format!(
                    "mixed/{}/{}",
                    if rows_state.is_empty() {
                        "empty"
                    } else {
                        "one"
                    },
                    repetition
                );
                let mut block = parse_block(MIXED);
                // Pre-prove the parsed mixed loop's exact identity.
                assert_eq!(block.name(), "bad", "{label}: name");
                assert_eq!(block.items().len(), 3, "{label}: dimension");
                assert_eq!(
                    pair_frame(&block, 0),
                    Some(("_entry.id", 2, "before", 2, 11))
                );
                assert_eq!(pair_frame(&block, 2), Some(("_after.id", 7, "tail", 7, 11)));
                let CifItem::Loop(stale) = &block.items()[1] else {
                    panic!("{label}: loop index");
                };
                assert_eq!(stale.line(), 3, "{label}: loop line");
                assert_eq!(
                    stale.tags(),
                    ["_entity.id", "_foreign.keep"],
                    "{label}: exact tags"
                );
                assert_eq!(
                    stale
                        .values()
                        .iter()
                        .map(|v| (v.line(), v.column()))
                        .collect::<Vec<_>>(),
                    [(6, 1), (6, 5)],
                    "{label}: literal value positions"
                );
                assert_eq!(
                    stale.values().iter().map(|v| v.raw()).collect::<Vec<_>>(),
                    ["OLD", "stale"],
                    "{label}: loop raws"
                );
                let entities: Vec<BioEntityRow> = rows_state
                    .iter()
                    .map(|(name, kind)| entity_row(*kind, name))
                    .collect();
                for (index, (name, kind)) in rows_state.iter().enumerate() {
                    verify_entity(&entities[index], *kind, name, &label);
                }
                let before_entities = entities.clone();
                let before = block.clone();

                let result = write_entity_category(&entities, &mut block);
                calls += 1;

                // Same-call preservation BEFORE error field matching.
                if entities != before_entities {
                    discrepancies.push(format!("{label}: entities mutated"));
                }
                if block != before {
                    discrepancies.push(format!("{label}: block mutated on error"));
                }

                // Borrow the error; assert ALL exact typed fields.
                let Some(error) = result.as_ref().err() else {
                    discrepancies.push(format!("{label}: expected error, got Ok"));
                    continue;
                };
                errors_observed += 1;
                if error.kind() != CifReadErrorKind::InvalidLoop {
                    discrepancies.push(format!("{label}: kind {:?}", error.kind()));
                }
                if error.source() != "cif" {
                    discrepancies.push(format!("{label}: source {}", error.source()));
                }
                if error.line() != 3 {
                    discrepancies.push(format!("{label}: line {}", error.line()));
                }
                if error.column() != 1 {
                    discrepancies.push(format!("{label}: column {}", error.column()));
                }
                if error.message() != "Tag _foreign.keep in loop with _entity." {
                    discrepancies.push(format!("{label}: message {}", error.message()));
                }
                // std::error::Error::source()==None per current CifReadError.
                if std::error::Error::source(error).is_some() {
                    discrepancies.push(format!("{label}: unexpected source chain"));
                }
            }
        }

        assert_eq!(calls, 4, "exact 4 mixed-error calls");
        assert_eq!(errors_observed, 4, "all four errors OBSERVED");
        assert!(
            discrepancies.is_empty(),
            "mixed-error discrepancies: {discrepancies:?}"
        );
    }
}
