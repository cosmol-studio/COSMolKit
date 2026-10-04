//! Private mmCIF `_database_PDB_matrix` (ORIGX) category emitter
//! (BIO-ORIGX-WRITE1-32).
//!
//! Complete selected-category helper only: the coordinate composer is NOT
//! dispatched here and caller group dispatch remains ROOT-owned. Reuses the
//! ONE existing category setter and the ONE existing double formatter; no
//! copied CIF mutation or conversion algorithm.

use cosmolkit_bio::BioTransform;

use super::super::cif::{CifBlock, format_cif_f64};

/// Source `nontrivial_origx` predicate. The two `is_identity` helper
/// bodies are copied here (cross-file rule); BioTransform has no public
/// `is_identity` and none is invented.
fn origx_is_nontrivial(has_origx: bool, origx: &BioTransform) -> bool {
    // BEGIN GEMMI CPP HELPERS Mat33::is_identity (math.hpp:232-236) + Transform::is_identity (math.hpp:408-410)
    // Gemmi✔️✔️:   bool is_identity() const {
    // Gemmi✔️✔️:     return a[0][0] == 1 && a[0][1] == 0 && a[0][2] == 0 &&
    // Gemmi✔️✔️:            a[1][0] == 0 && a[1][1] == 1 && a[1][2] == 0 &&
    // Gemmi✔️✔️:            a[2][0] == 0 && a[2][1] == 0 && a[2][2] == 1;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   bool is_identity() const {
    // Gemmi✔️✔️:     return mat.is_identity() && vec.x == 0. && vec.y == 0. && vec.z == 0.;
    // Gemmi✔️✔️:   }
    // END GEMMI CPP HELPERS
    //
    // Behavior review: the ordered nine matrix comparisons run first
    // (diagonal == 1, off-diagonal == 0), then the three translation
    // comparisons == 0, short-circuiting exactly as the source's
    // `mat.is_identity() && vec.x == 0. && ...` does. IEEE `==` semantics
    // are preserved verbatim — NOT bits, epsilon, or approximation:
    // off-diagonal and vector `-0.0` compare equal to zero (a signed-zero
    // transform IS the identity here), while a NaN component compares
    // equal to nothing, making is_identity false and the transform
    // nontrivial (ACTIVE), matching the source exactly. No public method
    // is added to BioTransform; the predicate reads the existing borrowed
    // `matrix()`/`translation()` views.
    //
    // Cost review: twelve fixed scalar comparisons with short-circuiting,
    // O(1), zero allocation — identical to the source's inlined
    // comparisons. ✔️ on both axes for the predicate itself.
    let matrix = origx.matrix();
    let translation = origx.translation();
    let matrix_is_identity = matrix[0][0] == 1.0
        && matrix[0][1] == 0.0
        && matrix[0][2] == 0.0
        && matrix[1][0] == 0.0
        && matrix[1][1] == 1.0
        && matrix[1][2] == 0.0
        && matrix[2][0] == 0.0
        && matrix[2][1] == 0.0
        && matrix[2][2] == 1.0;
    has_origx
        && !(matrix_is_identity
            && translation[0] == 0.0
            && translation[1] == 0.0
            && translation[2] == 0.0)
}

pub(super) fn write_origx_category(
    has_origx: bool,
    origx: &BioTransform,
    entry_id_raw: &str,
    block: &mut CifBlock,
) {
    // BEGIN GEMMI CPP FUNCTION gemmi::update_mmcif_block origx branch (to_mmcif.cpp:977-993)
    // Gemmi✔️✔️:   bool nontrivial_origx = st.has_origx && !st.origx.is_identity();
    // Gemmi❌❌:   if (groups.origx && nontrivial_origx) { // _database_PDB_matrix (ORIGX)
    // Gemmi❗❌:     cif::ItemSpan span(block.items, "_database_PDB_matrix.");
    // Gemmi❗❌:     span.set_pair("_database_PDB_matrix.entry_id", id);
    // Gemmi❗❌:     std::string tag_mat = "_database_PDB_matrix.origx[0][0]";
    // Gemmi❗❌:     std::string tag_vec = "_database_PDB_matrix.origx_vector[0]";
    // Gemmi✔️✔️:     for (int i = 0; i < 3; ++i) {
    // Gemmi❗❌:       tag_mat[27] += 1;  // origx[0] -> origx[1] -> origx[2]
    // Gemmi❗❌:       tag_vec[34] += 1;
    // Gemmi✔️✔️:       for (int j = 0; j < 3; ++j) {
    // Gemmi❗❌:         tag_mat[30] = '1' + j;
    // Gemmi❗❌:         span.set_pair(tag_mat, to_str(st.origx.mat[i][j]));
    // Gemmi✔️✔️:       }
    // Gemmi❗❌:       span.set_pair(tag_vec, to_str(st.origx.vec.at(i)));
    // Gemmi✔️✔️:     }
    // Gemmi❌❌:   }
    //
    // END GEMMI CPP FUNCTION
    //
    // Guard scope: line 977's `nontrivial_origx` IS implemented — by the
    // private `origx_is_nontrivial` predicate above, with both source
    // `is_identity` helper bodies copied inside it and separate reviews.
    // The predicate runs BEFORE any span/set: a false flag or an exact
    // IEEE identity leaves the WHOLE block untouched (a stale mixed loop
    // is preserved). ONLY the caller `groups.origx` bit on line 978 is
    // outside this selected-category boundary (composer undispatched,
    // ROOT-owned open) — that omission alone carries ❌❌ there. No
    // finite-input filter, sanitize, normalization or fallback exists.
    //
    // Helper bodies consulted (cross-file rule, NOT duplicated): the
    // ItemSpan prefix constructor + set_pair (cifdoc.hpp:613-645) are the
    // existing CifBlock::set_pair_in_category; `to_str` (sprintf.hpp:36-40)
    // is the existing pub format_cif_f64 owner — its own in-body markers
    // and qualifications govern it; small-integer/-0 proofs here do NOT
    // prove all-f64 formatter equivalence.
    //
    // Behavior review: the active branch writes EXACT13 fields in source
    // order — the raw entry_id first, then for each row i the three
    // matrix entries followed by that row's vector component. The source
    // tag increments (tag_mat[27] row digit, tag_mat[30] = '1' + j column
    // digit, tag_vec[34] vector digit) produce exactly the thirteen
    // frozen literal names; this owner constructs those SAME names by
    // fixed index (i + 1 / j + 1 in the literal format) — an exact
    // equivalent construction, not a second stored representation. The
    // borrowed transform and entry raw are preserved; the entry raw is
    // forwarded unchanged. Infallible delegated-setter boundary: no
    // invented Result/error/source-chain.
    //
    // Cost review (separate axes): the source constructs ONE persistent
    // ItemSpan shared by all thirteen set_pair calls; every delegated
    // Rust set_pair_in_category call RECOMPUTES the category span —
    // THIRTEEN span recomputations versus one persistent span, a KNOWN
    // local extra scan (❌ on the span/set_pair lines), not unresolved ❗
    // and not whole-cost parity. Tag construction: the source reuses two
    // mutated std::string buffers while this owner allocates a fresh
    // String per tag (13 per active call) — stated truthfully as ❌ on
    // the tag lines; short std::string contents may additionally live in
    // SSO inline storage (ABI-dependent, unmeasured, separate caveat).
    // Each formatted value is likewise a fresh String from the delegated
    // formatter. The bounded i/j loops are fixed 3x3+3 iterations on both
    // sides (✔️✔️). The identity predicate is O(1) (see its own review).
    // No state cloning occurs in production.
    if origx_is_nontrivial(has_origx, origx) {
        block.set_pair_in_category(
            Some("_database_PDB_matrix."),
            "_database_PDB_matrix.entry_id",
            entry_id_raw.to_owned(),
        );
        let matrix = origx.matrix();
        let translation = origx.translation();
        for i in 0..3usize {
            for j in 0..3usize {
                let tag = format!("_database_PDB_matrix.origx[{}][{}]", i + 1, j + 1);
                block.set_pair_in_category(
                    Some("_database_PDB_matrix."),
                    &tag,
                    format_cif_f64(matrix[i][j]),
                );
            }
            let vector_tag = format!("_database_PDB_matrix.origx_vector[{}]", i + 1);
            block.set_pair_in_category(
                Some("_database_PDB_matrix."),
                &vector_tag,
                format_cif_f64(translation[i]),
            );
        }
    }
}

#[cfg(test)]
mod bio_origx_writer_tests {
    use super::write_origx_category;
    use crate::cif::{CifBlock, CifCheckLevel, CifItem, read_cif_document};
    use cosmolkit_bio::BioTransform;

    const F0: &str = "data_f0\n_entry.id before\n_after.id tail\n";
    const F1: &str = "data_f1\n_entry.id before\n_database_PDB_matrix.entry_id OLD\n_database_PDB_matrix.origx[1][1] OLD\n_database_PDB_matrix.origx[1][2] OLD\n_database_PDB_matrix.origx[1][3] OLD\n_database_PDB_matrix.origx_vector[1] OLD\n_database_PDB_matrix.origx[2][1] OLD\n_database_PDB_matrix.origx[2][2] OLD\n_database_PDB_matrix.origx[2][3] OLD\n_database_PDB_matrix.origx_vector[2] OLD\n_database_PDB_matrix.origx[3][1] OLD\n_database_PDB_matrix.origx[3][2] OLD\n_database_PDB_matrix.origx[3][3] OLD\n_database_PDB_matrix.origx_vector[3] OLD\n_after.id tail\n";
    const F2: &str = "data_f2\n_entry.id before\nloop_\n_database_PDB_matrix.entry_id\n_foreign.keep\nOLD stale\n_after.id tail\n";

    const ENTRY_TAG: &str = "_database_PDB_matrix.entry_id";
    /// Twelve numeric tags in production order; ENTRY_TAG is separate.
    const TAGS: [&str; 12] = [
        "_database_PDB_matrix.origx[1][1]",
        "_database_PDB_matrix.origx[1][2]",
        "_database_PDB_matrix.origx[1][3]",
        "_database_PDB_matrix.origx_vector[1]",
        "_database_PDB_matrix.origx[2][1]",
        "_database_PDB_matrix.origx[2][2]",
        "_database_PDB_matrix.origx[2][3]",
        "_database_PDB_matrix.origx_vector[2]",
        "_database_PDB_matrix.origx[3][1]",
        "_database_PDB_matrix.origx[3][2]",
        "_database_PDB_matrix.origx[3][3]",
        "_database_PDB_matrix.origx_vector[3]",
    ];

    fn parse_block(text: &str) -> CifBlock {
        read_cif_document(text, "origx_fixture", CifCheckLevel::Syntax)
            .expect("frozen block state parses under the Syntax owner")
            .blocks()[0]
            .clone()
    }

    /// Frozen transforms with twelve-component literal bit prerequisites.
    struct Profile {
        name: &'static str,
        transform: BioTransform,
        bits: [u64; 12],
        /// Active raw literals for the 12 numeric fields (empty => identity
        /// transform, never active). Frozen LITERALS — never computed via
        /// the tested formatter or production function.
        raws: [&'static str; 12],
        active: bool,
    }

    const B1: u64 = 0x3ff0000000000000;
    const B0: u64 = 0x0000000000000000;
    const BN0: u64 = 0x8000000000000000;

    fn profiles() -> [Profile; 4] {
        [
            Profile {
                name: "I",
                transform: BioTransform::new(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                    [0.0, 0.0, 0.0],
                ),
                bits: [B1, B0, B0, B0, B1, B0, B0, B0, B1, B0, B0, B0],
                raws: [""; 12],
                active: false,
            },
            Profile {
                name: "Z",
                transform: BioTransform::new(
                    [[1.0, -0.0, -0.0], [-0.0, 1.0, -0.0], [-0.0, -0.0, 1.0]],
                    [-0.0, -0.0, -0.0],
                ),
                bits: [B1, BN0, BN0, BN0, B1, BN0, BN0, BN0, B1, BN0, BN0, BN0],
                raws: [""; 12],
                active: false,
            },
            Profile {
                name: "M",
                transform: BioTransform::new(
                    [[2.0, 5.0, 6.0], [7.0, 3.0, 8.0], [9.0, 10.0, 4.0]],
                    [11.0, 12.0, 13.0],
                ),
                bits: [
                    2.0f64.to_bits(),
                    5.0f64.to_bits(),
                    6.0f64.to_bits(),
                    7.0f64.to_bits(),
                    3.0f64.to_bits(),
                    8.0f64.to_bits(),
                    9.0f64.to_bits(),
                    10.0f64.to_bits(),
                    4.0f64.to_bits(),
                    11.0f64.to_bits(),
                    12.0f64.to_bits(),
                    13.0f64.to_bits(),
                ],
                raws: [
                    "2", "5", "6", "11", "7", "3", "8", "12", "9", "10", "4", "13",
                ],
                active: true,
            },
            Profile {
                name: "V",
                transform: BioTransform::new(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                    [-0.0, 0.0, 1.25],
                ),
                bits: [
                    B1,
                    B0,
                    B0,
                    B0,
                    B1,
                    B0,
                    B0,
                    B0,
                    B1,
                    BN0,
                    B0,
                    1.25f64.to_bits(),
                ],
                raws: [
                    "1", "0", "0", "-0", "0", "1", "0", "0", "0", "0", "1", "1.25",
                ],
                active: true,
            },
        ]
    }

    /// Twelve-component bit snapshot in frozen order (row-major matrix,
    /// then translation).
    fn transform_bits(transform: &BioTransform) -> [u64; 12] {
        let matrix = transform.matrix();
        let translation = transform.translation();
        [
            matrix[0][0].to_bits(),
            matrix[0][1].to_bits(),
            matrix[0][2].to_bits(),
            matrix[1][0].to_bits(),
            matrix[1][1].to_bits(),
            matrix[1][2].to_bits(),
            matrix[2][0].to_bits(),
            matrix[2][1].to_bits(),
            matrix[2][2].to_bits(),
            translation[0].to_bits(),
            translation[1].to_bits(),
            translation[2].to_bits(),
        ]
    }

    /// Per-call twelve-component literal bit prerequisite.
    fn verify_bits(profile: &Profile, label: &str) {
        assert_eq!(
            transform_bits(&profile.transform),
            profile.bits,
            "{label}: frozen twelve-component bits"
        );
    }

    /// Complete parsed-frame prerequisite for F0/F1/F2.
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
                assert_eq!(block.items().len(), 15, "F1 literal dimension");
                assert_eq!(block.name(), "f1");
                let CifItem::Pair(entry) = &block.items()[0] else {
                    panic!("F1 entry position");
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
                // Frozen parsed positions from literal tag length + 2:
                // entry_id (len29) 3:31; origx[i][j] (len32) col34;
                // origx_vector[i] (len36) col38; lines 3..15.
                let expected_columns: [usize; 13] =
                    [31, 34, 34, 34, 38, 34, 34, 34, 38, 34, 34, 34, 38];
                for (index, column) in expected_columns.iter().enumerate() {
                    let tag = if index == 0 {
                        ENTRY_TAG
                    } else {
                        TAGS[index - 1]
                    };
                    let CifItem::Pair(pair) = &block.items()[1 + index] else {
                        panic!("F1 category item {index} position");
                    };
                    assert_eq!(
                        (
                            pair.tag(),
                            pair.value().map(|v| v.raw()),
                            pair.line(),
                            pair.value().map(|v| (v.line(), v.column()))
                        ),
                        (tag, Some("OLD"), 3 + index, Some((3 + index, *column))),
                        "F1 category item {index}"
                    );
                }
                let CifItem::Pair(tail) = &block.items()[14] else {
                    panic!("F1 tail position");
                };
                assert_eq!(
                    (
                        tail.tag(),
                        tail.value().map(|v| v.raw()),
                        tail.line(),
                        tail.value().map(|v| (v.line(), v.column()))
                    ),
                    ("_after.id", Some("tail"), 16, Some((16, 11)))
                );
            }
            "F2" => {
                assert_eq!(block.items().len(), 3, "F2 literal dimension");
                assert_eq!(block.name(), "f2");
                let CifItem::Pair(entry) = &block.items()[0] else {
                    panic!("F2 entry position");
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
                    panic!("F2 stale loop position");
                };
                assert_eq!(stale.line(), 3, "F2 loop line");
                assert_eq!(stale.tags(), [ENTRY_TAG, "_foreign.keep"], "F2 loop tags");
                assert_eq!(stale.width(), 2, "F2 loop width");
                assert_eq!(stale.len(), 1, "F2 loop rows");
                assert_eq!(
                    stale.values().iter().map(|v| v.raw()).collect::<Vec<_>>(),
                    ["OLD", "stale"],
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
                    (
                        tail.tag(),
                        tail.value().map(|v| v.raw()),
                        tail.line(),
                        tail.value().map(|v| (v.line(), v.column()))
                    ),
                    ("_after.id", Some("tail"), 7, Some((7, 11)))
                );
            }
            other => panic!("unknown frozen block state {other}"),
        }
    }

    #[test]
    fn bio_origx_writer_product48() {
        let entries: [&str; 2] = ["'ENTRY A'", "?"];
        let blocks: [(&str, &str); 3] = [("F0", F0), ("F1", F1), ("F2", F2)];

        let mut calls = 0usize;
        let mut active_calls = 0usize;
        let mut noop_calls = 0usize;
        let mut observed_new_values = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        for profile in &profiles() {
            for has_origx in [false, true] {
                for entry_raw in entries {
                    for (block_name, block_text) in blocks {
                        let label =
                            format!("{}/{has_origx}/{entry_raw}/{block_name}", profile.name);
                        // Twelve-component bit prerequisite inside the
                        // innermost loop (includes all Z negative-zero bits).
                        verify_bits(profile, &label);
                        let entry_bytes = entry_raw.as_bytes().to_vec();
                        let transform_before_bits = transform_bits(&profile.transform);

                        let mut block = parse_block(block_text);
                        verify_parsed(&block, block_name);
                        let before = block.clone();

                        write_origx_category(has_origx, &profile.transform, entry_raw, &mut block);
                        calls += 1;

                        // Same-call input preservation BEFORE output checks.
                        if transform_bits(&profile.transform) != transform_before_bits {
                            discrepancies.push(format!("{label}: transform bits mutated"));
                        }
                        if entry_raw.as_bytes() != entry_bytes {
                            discrepancies.push(format!("{label}: entry raw bytes mutated"));
                        }
                        if block.name() != before.name() {
                            discrepancies.push(format!("{label}: block name changed"));
                        }

                        // Duplicate/loop counts BEFORE safe missing-field
                        // handling.
                        let mut duplicates = 0usize;
                        for tag in std::iter::once(&ENTRY_TAG).chain(TAGS.iter()) {
                            let count = block
                                .items()
                                .iter()
                                .filter(|item| {
                                    matches!(item, CifItem::Pair(p) if p
                                        .tag()
                                        .eq_ignore_ascii_case(tag))
                                })
                                .count();
                            if count > 1 {
                                discrepancies.push(format!("{label}: {tag} duplicates {count}"));
                            }
                            duplicates += count;
                        }
                        let category_loops = block
                            .items()
                            .iter()
                            .filter(|item| {
                                matches!(item, CifItem::Loop(l) if l
                                    .tags()
                                    .first()
                                    .is_some_and(|tag| tag
                                        .starts_with("_database_PDB_matrix.")))
                            })
                            .count();
                        if category_loops != 0 {
                            discrepancies
                                .push(format!("{label}: remaining category loop {category_loops}"));
                        }

                        let is_active = has_origx && profile.active;
                        if !is_active {
                            // False flag or exact IEEE identity: WHOLE block
                            // untouched (stale mixed loop preserved).
                            if block != before {
                                discrepancies.push(format!("{label}: no-op mutated whole block"));
                            }
                            noop_calls += 1;
                            continue;
                        }

                        active_calls += 1;
                        // Exact13 written sequence: (entry raw) + 12 numeric.
                        let written: Vec<(&str, &str)> = std::iter::once((ENTRY_TAG, entry_raw))
                            .chain(TAGS.iter().copied().zip(profile.raws.iter().copied()))
                            .collect();

                        match block_name {
                            "F0" => {
                                // Appends after tail at 2..14.
                                if block.items().len() != 15 {
                                    discrepancies.push(format!(
                                        "{label}: item count {} != 15",
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
                            "F1" => {
                                // 15 items; category 1..13 rewritten with
                                // retained Pair.line3..15/value0:0.
                                if block.items().len() != 15 {
                                    discrepancies.push(format!(
                                        "{label}: item count {} != 15",
                                        block.items().len()
                                    ));
                                }
                                if before.items().get(0) != block.items().get(0) {
                                    discrepancies.push(format!("{label}: unrelated entry changed"));
                                }
                                if before.items().get(14) != block.items().get(14) {
                                    discrepancies.push(format!("{label}: unrelated tail changed"));
                                }
                                for (offset, (tag, raw)) in written.iter().enumerate() {
                                    match block.items().get(1 + offset) {
                                        Some(CifItem::Pair(pair)) => {
                                            if pair.tag() != *tag {
                                                discrepancies.push(format!(
                                                    "{label}: {tag} tag {}",
                                                    pair.tag()
                                                ));
                                            }
                                            if pair.value().map(|v| v.raw()) != Some(*raw) {
                                                discrepancies.push(format!(
                                                    "{label}: {tag} raw {:?} != {raw}",
                                                    pair.value().map(|v| v.raw())
                                                ));
                                            }
                                            if pair.line() != 3 + offset {
                                                discrepancies.push(format!(
                                                    "{label}: {tag} line {} != {}",
                                                    pair.line(),
                                                    3 + offset
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
                                        _ => discrepancies.push(format!(
                                            "{label}: no {tag} pair at {}",
                                            1 + offset
                                        )),
                                    }
                                }
                            }
                            _ => {
                                // F2: loop replaced by entry Pair at 1; 12
                                // numeric at 2..13; tail 2->14 WHOLE.
                                if block.items().len() != 15 {
                                    discrepancies.push(format!(
                                        "{label}: item count {} != 15",
                                        block.items().len()
                                    ));
                                }
                                if before.items().get(0) != block.items().get(0) {
                                    discrepancies.push(format!("{label}: unrelated entry changed"));
                                }
                                if before.items().get(2) != block.items().get(14) {
                                    discrepancies.push(format!("{label}: old tail 2->14 changed"));
                                }
                                for (offset, (tag, raw)) in written.iter().enumerate() {
                                    match block.items().get(1 + offset) {
                                        Some(CifItem::Pair(pair)) => {
                                            if pair.tag() != *tag {
                                                discrepancies.push(format!(
                                                    "{label}: {tag} tag {}",
                                                    pair.tag()
                                                ));
                                            }
                                            if pair.value().map(|v| v.raw()) != Some(*raw) {
                                                discrepancies.push(format!(
                                                    "{label}: {tag} raw {:?} != {raw}",
                                                    pair.value().map(|v| v.raw())
                                                ));
                                            }
                                            if pair.line() != 0 {
                                                discrepancies.push(format!(
                                                    "{label}: {tag} line {} != 0",
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
                                        _ => discrepancies.push(format!(
                                            "{label}: no {tag} pair at {}",
                                            1 + offset
                                        )),
                                    }
                                }
                            }
                        }
                        // OBSERVED NEW-value census from the actual block:
                        // count complete expected-tag Pairs carrying the
                        // frozen raws written THIS call.
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

        assert_eq!(calls, 48, "exact 48 real write_origx_category calls");
        assert_eq!(active_calls, 12, "12 active publications (has=true, M/V)");
        assert_eq!(noop_calls, 36, "36 whole-block no-ops");
        assert_eq!(
            observed_new_values, 156,
            "observed 156 newly written values (12x13)"
        );
        assert!(
            discrepancies.is_empty(),
            "product discrepancies: {discrepancies:?}"
        );
    }
}
