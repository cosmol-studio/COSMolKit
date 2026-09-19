//! V3000 read behavior regressions.
//!
//! Named behavior slices map to pinned RDKit 2026.03.1
//! `MolFileParser.cpp` / `MolSGroupParsing.cpp` and are exercised through the
//! public detached reader; private helpers stay in their owning module.

use cosmolkit_io::{MolBlockRecord, SdfReadError, read_mol_block_detached};
use cosmolkit_model::{AtomQueryPredicate, QueryNode};

#[test]
fn review_utf8_attchord_byte_truncation_is_a_parse_error() {
    // RDKit substr(1, size - 2) cuts a byte. The Rust label must remain UTF-8;
    // invalid byte truncation is a structured error, not a panic or lossy label.
    let block = v3000_block(&["M  V30 1 C 0 0 0 0 ATTCHORD=(2 1 é"], 1);
    assert!(matches!(
        read_mol_block_detached(&block),
        Err(SdfReadError::Parse(_))
    ));
    let block = v3000_block(&["M  V30 1 C 0 0 0 0 ATTCHORD=(2 1 é)"], 1);
    let MolBlockRecord::Concrete { topology, .. } = read_mol_block_detached(&block).unwrap() else {
        panic!("concrete record expected");
    };
    assert_eq!(
        topology.atoms[0]
            .template_attachment_order()
            .unwrap()
            .entries()[0]
            .label(),
        "é"
    );
}

#[test]
fn review_utf8_unknown_blocks_and_collections_are_skipped() {
    // Pinned ParseV3000CTAB skips unknown blocks; parseEnhancedStereo skips
    // nonmatching collection lines. Neither operation slices Unicode strings.
    for extra in [
        "M  V30 BEGIN😀ABCDEFG\nM  V30 END UNKNOWN\n",
        "M  V30 BEGIN COLLECTION\nM  V30 MDLV30/STE😀 ATOMS=(1 1)\nM  V30 END COLLECTION\n",
    ] {
        let block = v3000_block(&["M  V30 1 C 0 0 0 0"], 1)
            .replace("M  V30 END CTAB", &format!("{extra}M  V30 END CTAB"));
        let MolBlockRecord::Concrete { topology, .. } = read_mol_block_detached(&block).unwrap()
        else {
            panic!("concrete record expected");
        };
        assert_eq!(topology.atoms.len(), 1);
        assert!(topology.stereo_groups.is_empty());
    }
}

fn v3000_block(atom_lines: &[&str], atom_count: usize) -> String {
    let mut block = String::from(concat!(
        "\n",
        "  COSMolKit\n",
        "\n",
        "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
        "M  V30 BEGIN CTAB\n",
    ));
    block.push_str(&format!("M  V30 COUNTS {atom_count} 0 0 0 0\n"));
    block.push_str("M  V30 BEGIN ATOM\n");
    for line in atom_lines {
        block.push_str(line);
        block.push('\n');
    }
    block.push_str("M  V30 END ATOM\nM  V30 END CTAB\nM  END\n");
    block
}

fn review_collection_record(blocks: &[&str]) -> Result<MolBlockRecord, SdfReadError> {
    let mut extra = String::new();
    for lines in blocks {
        extra.push_str("M  V30 BEGIN COLLECTION\n");
        for line in lines.split('\n').filter(|line| !line.is_empty()) {
            extra.push_str(&format!("M  V30 {line}\n"));
        }
        extra.push_str("M  V30 END COLLECTION\n");
    }
    let block = v3000_block(&["M  V30 1 C 0 0 0 0", "M  V30 2 O 1 0 0 0"], 2)
        .replace("M  V30 END CTAB", &format!("{extra}M  V30 END CTAB"));
    read_mol_block_detached(&block)
}

#[test]
fn review_collection_formatted_unsigned_extraction() {
    // RDKit 2026.03.1 parseEnhancedStereo: stringstream >> unsigned int.
    // Expected rows verified against the pinned wheel with sanitize=False.
    for (text, expected) in [
        ("1 1x", vec![0]),
        ("2 1+2", vec![0, 1]),
        ("2 1", vec![0, 0]),
        ("2 1 ", vec![0, 0]),
        ("1 -4294967295", vec![0]),
        ("1 +1", vec![0]),
        ("1 01", vec![0]),
        ("2 1\t2", vec![0, 1]),
        ("2 1\u{b}2", vec![0, 1]),
        ("1 1 2", vec![0]),
    ] {
        let line = format!("MDLV30/STEREL1 ATOMS=({text})");
        let MolBlockRecord::Concrete { topology, .. } = review_collection_record(&[&line]).unwrap()
        else {
            panic!("concrete record expected");
        };
        let actual: Vec<_> = topology.stereo_groups[0]
            .atoms()
            .iter()
            .map(|id| id.index())
            .collect();
        assert_eq!(actual, expected, "{text:?}");
    }
    for text in [
        "2 1x",
        "1 0x1",
        "1 4294967296",
        "1 -1",
        "1 ",
        "1 +",
        "1 3",
        "2 1\u{a0}2",
    ] {
        let line = format!("MDLV30/STEREL1 ATOMS=({text})");
        assert!(
            matches!(
                review_collection_record(&[&line]),
                Err(SdfReadError::Parse(_))
            ),
            "{text:?}"
        );
    }
}

#[test]
fn review_collection_only_nonempty_stereo_replaces_previous_groups() {
    // parseEnhancedStereo installs only when !groups.empty(). A recognized
    // group with zero members is still one group and must replace the old one.
    for (second, expected_id, expected_atoms) in [
        ("", 1, vec![0]),
        ("MDLV30/HILITE ATOMS=(1 2)", 1, vec![0]),
        ("MDLV30/STEREL2 ATOMS=(1 2)", 2, vec![1]),
        ("MDLV30/STEREL2 ATOMS=(0 )", 2, vec![]),
    ] {
        let MolBlockRecord::Concrete { topology, .. } =
            review_collection_record(&["MDLV30/STEREL1 ATOMS=(1 1)", second]).unwrap()
        else {
            panic!("concrete record expected");
        };
        assert_eq!(topology.stereo_groups.len(), 1);
        assert_eq!(topology.stereo_groups[0].id(), Some(expected_id));
        let actual: Vec<_> = topology.stereo_groups[0]
            .atoms()
            .iter()
            .map(|id| id.index())
            .collect();
        assert_eq!(actual, expected_atoms, "{second}");
    }
}

#[test]
fn v3k_lines_joins_single_continuation_without_separator() {
    // Physical line 7 is `M  V30 1 C 0 0 -`; its continuation marker drops the
    // trailing dash and concatenates the next payload with no separator.
    let block = v3000_block(&["M  V30 1 C 0 0 -", "M  V30 0 0"], 1);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("continued atom row")
    else {
        panic!("ordinary atom must produce concrete topology");
    };
    assert_eq!(topology.atoms.len(), 1);
    assert_eq!(topology.atoms[0].atomic_number(), 6);
}

#[test]
fn v3k_lines_joins_multiple_continuations_in_order() {
    // Payload fragments `1 C 0 ` + `0 ` + `0 0` concatenate to a valid row.
    let block = v3000_block(&["M  V30 1 C 0 -", "M  V30 0 -", "M  V30 0 0"], 1);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("multi-continued atom row")
    else {
        panic!("ordinary atom must produce concrete topology");
    };
    assert_eq!(topology.atoms.len(), 1);
    assert_eq!(topology.atoms[0].atomic_number(), 6);
}

#[test]
fn v3k_lines_rejects_unprefixed_continuation_with_line_number() {
    // The third physical ATOM line lacks the `M  V30 ` prefix.
    let block = v3000_block(&["M  V30 1 C 0 -", "not a v30 line", "M  V30 0 0"], 1);
    match read_mol_block_detached(&block) {
        Err(SdfReadError::Parse(message)) => {
            assert!(
                message.contains("Line 9 does not start with 'M  V30 '"),
                "unexpected diagnostic: {message}"
            );
        }
        other => panic!("expected a structured prefix failure, got {other:?}"),
    }
}

#[test]
fn v3k_lines_reports_eof_on_truncated_continuation() {
    // The final physical line ends with a continuation marker with no next
    // line; the failure must carry the missing physical line number.
    let mut block = v3000_block(&["M  V30 1 C 0 -"], 1);
    // Drop the generated `M  V30 END ATOM\n...` tail so the record ends on the
    // continuation marker.
    if let Some(index) = block.find("M  V30 END ATOM") {
        block.truncate(index);
    }
    match read_mol_block_detached(&block) {
        Err(SdfReadError::Parse(message)) => {
            assert!(
                message.contains("Line 9 does not start with 'M  V30 '"),
                "unexpected diagnostic: {message}"
            );
        }
        other => panic!("expected a structured truncation failure, got {other:?}"),
    }
}

/// Build a V3000 record whose ATOM/BOND blocks are emitted only when the
/// corresponding slice is non-empty, so zero-count block omission is exercised
/// with genuinely absent blocks rather than empty BEGIN/END pairs.
fn v3000_with_outer_and_blocks(
    outer: &str,
    inner_counts: &str,
    atom_lines: &[&str],
    bond_lines: &[&str],
) -> String {
    let mut block = String::new();
    block.push_str("\n  COSMolKit\n\n");
    block.push_str(&format!("{outer} V3000\n"));
    block.push_str("M  V30 BEGIN CTAB\n");
    block.push_str(&format!("M  V30 COUNTS {inner_counts}\n"));
    if !atom_lines.is_empty() {
        block.push_str("M  V30 BEGIN ATOM\n");
        for line in atom_lines {
            block.push_str(line);
            block.push('\n');
        }
        block.push_str("M  V30 END ATOM\n");
    }
    if !bond_lines.is_empty() {
        block.push_str("M  V30 BEGIN BOND\n");
        for line in bond_lines {
            block.push_str(line);
            block.push('\n');
        }
        block.push_str("M  V30 END BOND\n");
    }
    block.push_str("M  V30 END CTAB\nM  END\n");
    block
}

fn v3000_with_outer_and_inner(outer: &str, inner_counts: &str, atom_lines: &[&str]) -> String {
    v3000_with_outer_and_blocks(outer, inner_counts, atom_lines, &[])
}

#[test]
fn v3k_counts_rejects_nonzero_outer_counts_in_strict_mode() {
    let block = v3000_with_outer_and_inner("  1  1  0  0  0  0  0  0  0  0999", "0 0 0 0 0", &[]);
    match read_mol_block_detached(&block) {
        Err(SdfReadError::Parse(message)) => assert!(
            message.contains("should have 0s"),
            "unexpected diagnostic: {message}"
        ),
        other => panic!("expected structured outer-count failure, got {other:?}"),
    }
}

#[test]
fn v3k_counts_allows_nonzero_outer_counts_in_non_strict_mode() {
    let block = v3000_with_outer_and_inner("  1  1  0  0  0  0  0  0  0  0999", "0 0 0 0 0", &[]);
    let params = cosmolkit_io::MolBlockReadParams {
        strict_parsing: false,
        ..cosmolkit_io::MolBlockReadParams::default()
    };
    let result = cosmolkit_io::read_mol_block_detached_with_params(&block, params);
    assert!(
        result.is_ok(),
        "non-strict outer counts must parse: {result:?}"
    );
}

#[test]
fn v3k_counts_requires_two_inner_count_fields() {
    let missing = v3000_with_outer_and_inner("  0  0  0  0  0  0  0  0  0  0999", "0", &[]);
    assert!(
        matches!(read_mol_block_detached(&missing), Err(SdfReadError::Counts)),
        "a single COUNTS field must be a structured Counts error"
    );
    let two_fields = v3000_with_outer_and_inner("  0  0  0  0  0  0  0  0  0  0999", "0 0", &[]);
    assert!(read_mol_block_detached(&two_fields).is_ok());
}

#[test]
fn v3k_counts_omits_zero_count_blocks_without_panic() {
    let block = v3000_with_outer_and_inner("  0  0  0  0  0  0  0  0  0  0999", "0 0 0 0 0", &[]);
    let MolBlockRecord::Concrete {
        topology,
        coordinates,
        properties,
    } = read_mol_block_detached(&block).expect("zero-count record")
    else {
        panic!("zero-count record must be concrete");
    };
    assert!(topology.atoms.is_empty());
    assert!(topology.bonds.is_empty());
    // RDKit still constructs the declared `Conformer(0)`; every coordinate row
    // (across both dimensions) must be absent, matching the atom table.
    assert!(
        coordinates
            .conformers_2d
            .iter()
            .all(|conformer| conformer.coordinates().is_empty())
    );
    assert!(coordinates.conformers_3d.is_empty());
    assert_eq!(properties.prop("_MolFileChiralFlag"), Some("0"));
}

#[test]
fn v3k_counts_rejects_malformed_counts_without_panic() {
    let block = v3000_with_outer_and_inner("  0  0  0  0  0  0  0  0  0  0999", "abc 0 0 0 0", &[]);
    assert!(
        matches!(
            read_mol_block_detached(&block),
            Err(SdfReadError::Field { .. })
        ),
        "malformed count text must be a structured field error"
    );
}

#[test]
fn v3k_counts_leading_plus_is_zero_atoms() {
    // Pinned `FileParserUtils::toUnsigned` passes the text to
    // `std::from_chars`, which does not recognize a leading `+`; RDKit therefore
    // reads `+1` as 0 atoms, so no ATOM block is required.
    let block = v3000_with_outer_and_inner("  0  0  0  0  0  0  0  0  0  0999", "+1 0 0 0 0", &[]);
    let MolBlockRecord::Concrete {
        topology,
        properties,
        ..
    } = read_mol_block_detached(&block).expect("leading plus must parse as zero atoms")
    else {
        panic!("zero-atom record must be concrete");
    };
    assert!(topology.atoms.is_empty());
    assert_eq!(properties.prop("_MolFileChiralFlag"), Some("0"));
}

#[test]
fn v3k_counts_unsigned_overflow_is_zero_atoms() {
    // `std::from_chars` reports out-of-range without modifying its initialized
    // target, so `2^32` reads as 0 atoms and must not allocate.
    let block = v3000_with_outer_and_inner(
        "  0  0  0  0  0  0  0  0  0  0999",
        "4294967296 0 0 0 0",
        &[],
    );
    let MolBlockRecord::Concrete {
        topology,
        properties,
        ..
    } = read_mol_block_detached(&block).expect("overflow must parse as zero atoms")
    else {
        panic!("zero-atom record must be concrete");
    };
    assert!(topology.atoms.is_empty());
    assert_eq!(properties.prop("_MolFileChiralFlag"), Some("0"));
}

#[test]
fn v3k_counts_chiral_flag_unsigned_boundaries() {
    // The maximum `unsigned int` is accepted and preserved. It is applied to
    // the chiral flag only: using it as an atom count would be an allocation
    // probe, which this suite intentionally avoids.
    let max = v3000_with_outer_and_inner(
        "  0  0  0  0  0  0  0  0  0  0999",
        "0 0 0 0 4294967295",
        &[],
    );
    let MolBlockRecord::Concrete { properties, .. } =
        read_mol_block_detached(&max).expect("max chiral flag")
    else {
        panic!("zero-atom record must be concrete");
    };
    assert_eq!(properties.prop("_MolFileChiralFlag"), Some("4294967295"));

    let overflow = v3000_with_outer_and_inner(
        "  0  0  0  0  0  0  0  0  0  0999",
        "0 0 0 0 4294967296",
        &[],
    );
    let MolBlockRecord::Concrete { properties, .. } =
        read_mol_block_detached(&overflow).expect("overflow chiral flag")
    else {
        panic!("zero-atom record must be concrete");
    };
    assert_eq!(properties.prop("_MolFileChiralFlag"), Some("0"));
}

#[test]
fn v3k_counts_optional_fields_absent_and_present() {
    let outer = "  0  0  0  0  0  0  0  0  0  0999";
    let atom = "M  V30 1 C 0.0 0.0 0.0 0";

    let absent = v3000_with_outer_and_inner(outer, "1 0", &[atom]);
    let MolBlockRecord::Concrete {
        topology,
        coordinates,
        properties,
    } = read_mol_block_detached(&absent).expect("absent optional COUNTS fields")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(topology.atoms.len(), 1);
    assert_eq!(coordinates.conformers_2d.len(), 1);
    assert_eq!(coordinates.conformers_2d[0].coordinates().len(), 1);
    assert_eq!(properties.prop("_MolFileChiralFlag"), Some("0"));

    let present = v3000_with_outer_and_inner(outer, "1 0 0 0 1", &[atom]);
    let MolBlockRecord::Concrete { properties, .. } =
        read_mol_block_detached(&present).expect("present optional COUNTS fields")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(properties.prop("_MolFileChiralFlag"), Some("1"));
}

#[test]
fn v3k_counts_zero_bond_record_omits_bond_block() {
    // Declared bond count 0 with no BEGIN BOND block must remain valid and must
    // not fabricate bonds.
    let block = v3000_with_outer_and_blocks(
        "  0  0  0  0  0  0  0  0  0  0999",
        "1 0 0 0 0",
        &["M  V30 1 C 0.0 0.0 0.0 0"],
        &[],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("zero-bond record")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(topology.atoms.len(), 1);
    assert!(topology.bonds.is_empty());
}

#[test]
fn v3k_counts_lowercase_keyword_is_accepted() {
    // The counts line is uppercased before `COUNTS ` is matched, so a lowercase
    // keyword in the source block is valid.
    let block = concat!(
        "\n  COSMolKit\n\n",
        "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
        "M  V30 BEGIN CTAB\n",
        "M  V30 counts 1 0 0 0 0\n",
        "M  V30 BEGIN ATOM\n",
        "M  V30 1 C 0.0 0.0 0.0 0\n",
        "M  V30 END ATOM\n",
        "M  V30 END CTAB\n",
        "M  END\n",
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(block).expect("lowercase counts keyword")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(topology.atoms.len(), 1);
}

#[test]
fn v3k_counts_rejects_invalid_optional_fields() {
    let outer = "  0  0  0  0  0  0  0  0  0  0999";
    for invalid in ["0 0 x 0 0", "0 0 0 x 0", "0 0 0 0 x"] {
        let block = v3000_with_outer_and_inner(outer, invalid, &[]);
        assert!(
            matches!(
                read_mol_block_detached(&block),
                Err(SdfReadError::Field { .. })
            ),
            "invalid optional field {invalid:?} must be a structured field error"
        );
    }
}

fn v3000_block_with_line_ending(line_ending: &str) -> String {
    let lines = [
        "",
        "  COSMolKit",
        "",
        "  0  0  0  0  0  0  0  0  0  0999 V3000",
        "M  V30 BEGIN CTAB",
        "M  V30 COUNTS 1 0 0 0 0",
        "M  V30 BEGIN ATOM",
        "M  V30 1 C 1-",
        "M  V30 0 0 0 0",
        "M  V30 END ATOM",
        "M  V30 END CTAB",
        "M  END",
    ];
    let mut block = lines.join(line_ending);
    block.push_str(line_ending);
    block
}

#[test]
fn v3k_lines_joins_without_inserting_whitespace() {
    // `1-` + `0 0 0 0` concatenates to `10 0 0 0`, so the x coordinate is 10.
    // If a separator were inserted the row would read `1 0 0 0 0` (x = 1).
    let block = v3000_block_with_line_ending("\n");
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("continued split token")
    else {
        panic!("ordinary atom must produce concrete topology");
    };
    assert_eq!(coordinates.conformers_2d.len(), 1);
    assert_eq!(coordinates.conformers_2d[0].coordinates()[0][0], 10.0);
}

#[test]
fn v3k_lines_handles_crlf_line_endings() {
    let block = v3000_block_with_line_ending("\r\n");
    let MolBlockRecord::Concrete {
        topology,
        coordinates,
        ..
    } = read_mol_block_detached(&block).expect("CRLF record")
    else {
        panic!("ordinary atom must produce concrete topology");
    };
    assert_eq!(topology.atoms.len(), 1);
    assert_eq!(coordinates.conformers_2d.len(), 1);
    assert_eq!(coordinates.conformers_2d[0].coordinates()[0][0], 10.0);
}

const ZERO_OUTER: &str = "  0  0  0  0  0  0  0  0  0  0999";

#[test]
fn v3k_atom_rows_nonsequential_bookmarks_resolve_by_bookmark() {
    // Bookmarks 10/20/30 are not row positions; the bond must resolve them
    // through the bookmark map, not by using the parsed value as a row index.
    let block = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "3 1 0 0 0",
        &[
            "M  V30 10 C 0.0 0.0 0.0 0",
            "M  V30 20 O 1.0 0.0 0.0 0",
            "M  V30 30 N 2.0 0.0 0.0 0",
        ],
        &["M  V30 1 1 10 20"],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("nonsequential bookmarks")
    else {
        panic!("atom/bond record must be concrete");
    };
    assert_eq!(topology.atoms.len(), 3);
    assert_eq!(topology.bonds.len(), 1);
    assert_eq!(topology.bonds[0].begin().index(), 0);
    assert_eq!(topology.bonds[0].end().index(), 1);
}

#[test]
fn v3k_atom_rows_duplicate_bookmark_keeps_first_atom() {
    // `setAtomBookmark` appends and `getAtomWithBookmark` returns the first
    // inserted atom, so the duplicate bookmark 1 resolves to row 0.
    let block = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "3 1 0 0 0",
        &[
            "M  V30 1 C 0.0 0.0 0.0 0",
            "M  V30 1 O 1.0 0.0 0.0 0",
            "M  V30 2 N 2.0 0.0 0.0 0",
        ],
        &["M  V30 1 1 1 2"],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("duplicate bookmark keeps first")
    else {
        panic!("atom/bond record must be concrete");
    };
    assert_eq!(topology.atoms.len(), 3);
    assert_eq!(topology.bonds[0].begin().index(), 0);
    assert_eq!(topology.bonds[0].end().index(), 2);
    assert_eq!(topology.atoms[0].atomic_number(), 6);
}

#[test]
fn v3k_atom_rows_zero_and_malformed_bookmarks_are_zero() {
    // `std::from_chars` leaves the initialized target at 0 for `abc`; the first
    // bookmark-0 association (row 0) is the one a bond resolves to.
    let block = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "3 1 0 0 0",
        &[
            "M  V30 0 C 0.0 0.0 0.0 0",
            "M  V30 abc O 1.0 0.0 0.0 0",
            "M  V30 2 N 2.0 0.0 0.0 0",
        ],
        &["M  V30 1 1 0 2"],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("malformed bookmark reads as zero")
    else {
        panic!("atom/bond record must be concrete");
    };
    assert_eq!(topology.bonds[0].begin().index(), 0);
    assert_eq!(topology.bonds[0].end().index(), 2);
}

#[test]
fn v3k_atom_rows_prefix_and_signed_bookmarks_follow_from_chars() {
    // `1x` parses as 1 (maximal decimal-digit prefix); `+1` is not a recognized
    // sign and parses as 0, so bookmark 1 resolves to row 0 and bookmark 2 to
    // row 2, leaving the `+1` row unreferenced at index 1.
    let block = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "3 1 0 0 0",
        &[
            "M  V30 1x C 0.0 0.0 0.0 0",
            "M  V30 +1 O 1.0 0.0 0.0 0",
            "M  V30 2 N 2.0 0.0 0.0 0",
        ],
        &["M  V30 1 1 1 2"],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("prefix/sign bookmark")
    else {
        panic!("atom/bond record must be concrete");
    };
    assert_eq!(topology.bonds[0].begin().index(), 0);
    assert_eq!(topology.bonds[0].end().index(), 2);
}

#[test]
fn v3k_atom_rows_short_row_is_bad_atom_line() {
    let block =
        v3000_with_outer_and_blocks(ZERO_OUTER, "1 0 0 0 0", &["M  V30 1 C 0.0 0.0 0.0"], &[]);
    match read_mol_block_detached(&block) {
        Err(SdfReadError::Parse(message)) => assert!(
            message.contains("Bad atom line"),
            "unexpected diagnostic: {message}"
        ),
        other => panic!("expected a short-row failure, got {other:?}"),
    }
}

#[test]
fn v3k_atom_rows_early_end_atom_is_bad_atom_line() {
    // Two atoms declared but only one row before `END ATOM`; the marker line is
    // consumed as the second atom row and fails the required-field check.
    let block =
        v3000_with_outer_and_blocks(ZERO_OUTER, "2 0 0 0 0", &["M  V30 1 C 0.0 0.0 0.0 0"], &[]);
    match read_mol_block_detached(&block) {
        Err(SdfReadError::Parse(message)) => assert!(
            message.contains("Bad atom line"),
            "unexpected diagnostic: {message}"
        ),
        other => panic!("expected an early-END-ATOM failure, got {other:?}"),
    }
}

#[test]
fn v3k_atom_rows_late_end_atom_is_reported() {
    // One atom declared but two rows present; the extra row is read by the
    // `END ATOM` check and reported.
    let block = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "1 0 0 0 0",
        &["M  V30 1 C 0.0 0.0 0.0 0", "M  V30 2 O 1.0 0.0 0.0 0"],
        &[],
    );
    match read_mol_block_detached(&block) {
        Err(SdfReadError::Parse(message)) => assert!(
            message.contains("END ATOM"),
            "unexpected diagnostic: {message}"
        ),
        other => panic!("expected a late-END-ATOM failure, got {other:?}"),
    }
}

#[test]
fn v3k_atom_rows_missing_bookmark_reference_is_structured_error() {
    // RDKit's `getAtomWithBookmark` is a debug-only precondition; release
    // behavior is undefined, so the detached reader uses a structured error.
    let block = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "2 1 0 0 0",
        &["M  V30 1 C 0.0 0.0 0.0 0", "M  V30 2 O 1.0 0.0 0.0 0"],
        &["M  V30 1 1 1 99"],
    );
    match read_mol_block_detached(&block) {
        Err(SdfReadError::Field { kind, .. }) => {
            assert_eq!(kind, "bond end atom bookmark");
        }
        other => panic!("expected a missing-reference failure, got {other:?}"),
    }
}

#[test]
fn v3k_atom_numbers_prefix_suffix_and_exponent() {
    // A non-zero z selects 3D storage, so all three axes are observable.
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "2 0 0 0 0",
        &["M  V30 1 C 1.5xyz 2.5 0 0", "M  V30 2 O 1e2 2E-1 3.0e+1 0"],
    );
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("suffix/exponent coordinates")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(coordinates.conformers_3d.len(), 1);
    let rows = coordinates.conformers_3d[0].coordinates();
    assert_eq!(rows[0], [1.5, 2.5, 0.0]);
    assert_eq!(rows[1], [100.0, 0.2, 30.0]);
}

#[test]
fn v3k_atom_numbers_leading_signs_and_invalid_text() {
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "3 0 0 0 0",
        &[
            "M  V30 1 C +1.5 -2.5 0 0",
            "M  V30 2 O .5 0 0 0",
            "M  V30 3 N abc def 0 0",
        ],
    );
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("signs and invalid text")
    else {
        panic!("atom record must be concrete");
    };
    let rows = coordinates.conformers_2d[0].coordinates();
    assert_eq!(rows[0], [1.5, -2.5]);
    assert_eq!(rows[1], [0.5, 0.0]);
    assert_eq!(rows[2], [0.0, 0.0]);
}

#[test]
fn v3k_atom_numbers_nonfinite_coordinates_are_rejected() {
    for value in ["inf", "nan", "1e999"] {
        let block = v3000_with_outer_and_inner(
            ZERO_OUTER,
            "1 0 0 0 0",
            &[&format!("M  V30 1 C {value} 0 0 0")],
        );
        assert!(
            matches!(
                read_mol_block_detached(&block),
                Err(SdfReadError::Coordinates(_))
            ),
            "non-finite coordinate {value:?} must fail the finite-coordinate boundary"
        );
    }
}

#[test]
fn v3k_atom_numbers_atom_map_positive_only_with_prefix_text() {
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "5 0 0 0 0",
        &[
            "M  V30 1 C 0 0 0 5x",
            "M  V30 2 O 0 0 0 5.9",
            "M  V30 3 N 0 0 0 -3",
            "M  V30 4 F 0 0 0 0",
            "M  V30 5 Cl 0 0 0 +7",
        ],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("atom map conversions")
    else {
        panic!("atom record must be concrete");
    };
    let maps = topology
        .atoms
        .iter()
        .map(|atom| atom.atom_map())
        .collect::<Vec<_>>();
    assert_eq!(maps, vec![Some(5), Some(5), None, None, Some(7)]);
}

#[test]
fn v3k_atom_numbers_atom_map_c_locale_whitespace_and_no_conversion() {
    // Quoting keeps whitespace inside the V3000 token, after which RDKit's
    // C-locale atoi accepts SP/HT/VT/FF/CR. LF is covered at the private
    // helper boundary because it terminates a physical molfile line.
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "9 0 0 0 0",
        &[
            "M  V30 1 C 0 0 0 \" 7\"",
            "M  V30 2 C 0 0 0 \"\t7\"",
            "M  V30 3 C 0 0 0 \"\u{000b}7\"",
            "M  V30 4 C 0 0 0 \"\u{000c}7\"",
            "M  V30 5 C 0 0 0 \"\r7\"",
            "M  V30 6 C 0 0 0 \"\"",
            "M  V30 7 C 0 0 0 \"+\"",
            "M  V30 8 C 0 0 0 \"-\"",
            "M  V30 9 C 0 0 0 abc",
        ],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("C-locale map conversion")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(
        topology
            .atoms
            .iter()
            .map(|atom| atom.atom_map())
            .collect::<Vec<_>>(),
        vec![
            Some(7),
            Some(7),
            Some(7),
            Some(7),
            Some(7),
            None,
            None,
            None,
            None,
        ]
    );
}

#[test]
fn v3k_atom_numbers_atom_map_fixed_glibc_int_width_boundaries() {
    // Pinned RDKit 2026.03.1 on x86_64 glibc 2.43 was observed on every row.
    // Values outside C `int` range are fixed-environment regressions, not a
    // claim that C defines portable atoi behavior for those inputs.
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "9 0 0 0 0",
        &[
            "M  V30 1 C 0 0 0 2147483647",
            "M  V30 2 C 0 0 0 2147483648",
            "M  V30 3 C 0 0 0 4294967296",
            "M  V30 4 C 0 0 0 4294967297",
            "M  V30 5 C 0 0 0 9223372036854775807",
            "M  V30 6 C 0 0 0 99999999999999999999",
            "M  V30 7 C 0 0 0 -2147483649",
            "M  V30 8 C 0 0 0 -4294967295",
            "M  V30 9 C 0 0 0 -99999999999999999999",
        ],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("fixed glibc map boundaries")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(
        topology
            .atoms
            .iter()
            .map(|atom| atom.atom_map())
            .collect::<Vec<_>>(),
        vec![
            Some(i32::MAX as u32),
            None,
            None,
            Some(1),
            None,
            None,
            Some(i32::MAX as u32),
            Some(1),
            None,
        ]
    );
}

#[test]
fn v3k_atom_numbers_coordinate_row_alignment() {
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "3 0 0 0 0",
        &[
            "M  V30 1 C 0.0 1.0 0 0",
            "M  V30 2 O 2.0 3.0 0 0",
            "M  V30 3 N 4.0 5.0 0 0",
        ],
    );
    let MolBlockRecord::Concrete {
        topology,
        coordinates,
        ..
    } = read_mol_block_detached(&block).expect("coordinate row alignment")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(topology.atoms.len(), 3);
    assert_eq!(coordinates.conformers_2d[0].coordinates().len(), 3);
    assert_eq!(coordinates.conformers_2d[0].coordinates()[2], [4.0, 5.0]);
}

#[test]
fn v3k_symbols_elements_mixed_case_second_char_is_normalized() {
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "3 0 0 0 0",
        &[
            "M  V30 1 CL 0 0 0 0",
            "M  V30 2 BR 0 0 0 0",
            "M  V30 3 NA 0 0 0 0",
        ],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("mixed-case elements")
    else {
        panic!("atom record must be concrete");
    };
    let numbers = topology
        .atoms
        .iter()
        .map(|atom| atom.element().atomic_number())
        .collect::<Vec<_>>();
    assert_eq!(numbers, vec![17, 35, 11]);
}

#[test]
fn v3k_symbols_elements_deuterium_and_tritium_isotopes() {
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "2 0 0 0 0",
        &["M  V30 1 D 0 0 0 0", "M  V30 2 T 0 0 0 0"],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("D/T shorthand")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(topology.atoms[0].element().atomic_number(), 1);
    assert_eq!(topology.atoms[0].isotope(), Some(2));
    assert_eq!(topology.atoms[1].element().atomic_number(), 1);
    assert_eq!(topology.atoms[1].isotope(), Some(3));
}

#[test]
fn v3k_symbols_elements_dummy_labels_for_pol_and_mod() {
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "2 0 0 0 0",
        &["M  V30 1 Pol 0 0 0 0", "M  V30 2 Mod 0 0 0 0"],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("Pol/Mod dummy labels")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(topology.atoms[0].element().atomic_number(), 0);
    assert_eq!(topology.atoms[0].prop("dummyLabel"), Some("Pol"));
    assert_eq!(topology.atoms[1].element().atomic_number(), 0);
    assert_eq!(topology.atoms[1].prop("dummyLabel"), Some("Mod"));
}

#[test]
fn v3k_symbols_elements_unknown_symbol_strictness() {
    let unknown = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "2 0 0 0 0",
        &["M  V30 1 Zz 0 0 0 0", "M  V30 2 cl 0 0 0 0"],
    );
    match read_mol_block_detached(&unknown) {
        Err(SdfReadError::Parse(message)) => {
            assert!(
                message.contains("Element 'Zz' not found"),
                "unexpected diagnostic: {message}"
            );
        }
        other => panic!("expected a strict unknown-symbol failure, got {other:?}"),
    }

    let params = cosmolkit_io::MolBlockReadParams {
        strict_parsing: false,
        ..cosmolkit_io::MolBlockReadParams::default()
    };
    let MolBlockRecord::Concrete { topology, .. } =
        cosmolkit_io::read_mol_block_detached_with_params(&unknown, params)
            .expect("non-strict unknown symbols")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(topology.atoms[0].element().atomic_number(), 0);
    assert_eq!(topology.atoms[0].prop("dummyLabel"), Some("Zz"));
    assert_eq!(topology.atoms[1].element().atomic_number(), 0);
    assert_eq!(topology.atoms[1].prop("dummyLabel"), Some("cl"));
}

#[test]
fn v3k_symbols_elements_empty_symbol_row_is_malformed() {
    // The doubled-quote tokenizer rule preserves `""` as the symbol token; it is
    // not a valid element, so strict parsing reports a structured parse error
    // rather than fabricating an atom.
    let block = v3000_with_outer_and_inner(ZERO_OUTER, "1 0 0 0 0", &["M  V30 1 \"\" 0 0 0 0"]);
    assert!(
        matches!(read_mol_block_detached(&block), Err(SdfReadError::Parse(_))),
        "an empty quoted symbol must fail as a structured parse error"
    );
}

fn v3000_single_atom(token: &str) -> String {
    let atom = format!("M  V30 1 {token} 0 0 0 0");
    v3000_with_outer_and_inner(ZERO_OUTER, "1 0 0 0 0", &[&atom])
}

#[test]
fn v3k_symbols_lists_positive_with_spaces_and_empty_members() {
    // Quoting keeps `[C,  N,,]` as a single token; the inner spaces and empty
    // members are stripped/skipped.
    let block = v3000_single_atom("\"[C,  N,,]\"");
    let MolBlockRecord::Query(record) = read_mol_block_detached(&block).expect("atom list query")
    else {
        panic!("atom list must produce a query record");
    };
    let query_atom = record.query.atom(0).expect("query atom");
    let debug = format!("{:?}", query_atom.predicate());
    assert!(debug.contains("AtomicNumberIn([6, 7])"), "{debug}");
    assert_eq!(query_atom.atom().element().atomic_number(), 0);
    assert!(!query_atom.atom().no_implicit());
}

#[test]
fn v3k_symbols_lists_not_case_variants() {
    for token in ["NOT[C,N]", "not[C,N]"] {
        let block = v3000_single_atom(token);
        let MolBlockRecord::Query(record) = read_mol_block_detached(&block).expect("NOT atom list")
        else {
            panic!("NOT atom list must produce a query record");
        };
        let debug = format!(
            "{:?}",
            record.query.atom(0).expect("query atom").predicate()
        );
        assert!(
            debug.contains("AtomicNumberNotIn([6, 7])"),
            "{token}: {debug}"
        );
    }
}

#[test]
fn v3k_symbols_lists_not_on_non_list_is_rejected() {
    let block = v3000_single_atom("NOTC");
    match read_mol_block_detached(&block) {
        Err(SdfReadError::Parse(message)) => assert!(
            message.contains("NOT tokens only supported for atom lists"),
            "unexpected diagnostic: {message}"
        ),
        other => panic!("expected a NOT-on-non-list failure, got {other:?}"),
    }
}

#[test]
fn v3k_symbols_lists_all_empty_is_structurally_rejected() {
    for token in ["[]", "[,,]"] {
        let block = v3000_single_atom(token);
        assert!(
            matches!(read_mol_block_detached(&block), Err(SdfReadError::Parse(_))),
            "an all-empty atom list {token:?} must be a structured error"
        );
    }
}

#[test]
fn v3k_symbols_lists_missing_bracket_is_rejected() {
    let block = v3000_single_atom("[C,N");
    match read_mol_block_detached(&block) {
        Err(SdfReadError::Parse(message)) => {
            assert!(
                message.contains("Bad atom token"),
                "unexpected diagnostic: {message}"
            );
        }
        other => panic!("expected a malformed atom-list failure, got {other:?}"),
    }
}

#[test]
fn v3k_symbols_queries_complex_names_and_wildcard_are_queries_without_implicit_hydrogens() {
    // RDKit 2026.03.1 ParseV3000AtomSymbol routes all eight complex query
    // names through convertComplexNameToQuery and `*` through AtomNull; only
    // these symbol-query branches explicitly set noImplicit.
    let symbols = ["A", "AH", "Q", "QH", "X", "XH", "M", "MH", "*"];
    let atom_lines = symbols
        .iter()
        .enumerate()
        .map(|(index, symbol)| format!("M  V30 {} {symbol} 0 0 0 0", index + 1))
        .collect::<Vec<_>>();
    let atom_refs = atom_lines.iter().map(String::as_str).collect::<Vec<_>>();
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        &format!("{} 0 0 0 0", symbols.len()),
        &atom_refs,
    );
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&block).expect("complex V3000 query symbols")
    else {
        panic!("complex symbols must produce a query record");
    };
    for (index, symbol) in symbols.iter().enumerate() {
        let atom = record.query.atom(index).expect("query atom");
        assert!(atom.atom().no_implicit(), "{symbol}");
        assert_eq!(atom.atom().element().atomic_number(), 0, "{symbol}");
    }
    let wildcard = format!("{:?}", record.query.atom(8).unwrap().predicate());
    assert!(wildcard.contains("Any"), "{wildcard}");
}

#[test]
fn v3k_symbols_queries_r_groups_remain_concrete_with_source_metadata() {
    // ParseV3000AtomSymbol constructs ordinary atomic-number-zero atoms for
    // R/R#/R-labels, stores the full dummy label, and stores positive parsed
    // suffixes as isotopes. The source's string comparison also admits R100.
    let symbols = ["R", "R#", "R0", "R7", "R99", "R100"];
    let atom_lines = symbols
        .iter()
        .enumerate()
        .map(|(index, symbol)| format!("M  V30 {} {symbol} 0 0 0 0", index + 1))
        .collect::<Vec<_>>();
    let atom_refs = atom_lines.iter().map(String::as_str).collect::<Vec<_>>();
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        &format!("{} 0 0 0 0", symbols.len()),
        &atom_refs,
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("R-group symbols")
    else {
        panic!("R-group symbols must remain concrete atoms");
    };
    let expected_isotopes = [None, None, None, Some(7), Some(99), Some(100)];
    for ((atom, symbol), expected_isotope) in
        topology.atoms.iter().zip(symbols).zip(expected_isotopes)
    {
        assert_eq!(atom.element().atomic_number(), 0, "{symbol}");
        assert_eq!(atom.prop("dummyLabel"), Some(symbol), "{symbol}");
        assert_eq!(atom.isotope(), expected_isotope, "{symbol}");
        assert!(!atom.no_implicit(), "{symbol}");
    }
}

#[test]
fn v3k_symbols_queries_generic_labels_preserve_the_complete_source_key_set() {
    // This is the complete key set of RDKit 2026.03.1
    // GenericGroups::genericMatchers. Parsing only preserves the query atom
    // and atomLabel; it does not run or duplicate generic matching.
    let labels = [
        "Group",
        "G",
        "GroupH",
        "GH",
        "Group*",
        "G*",
        "GroupH*",
        "GH*",
        "Alkyl",
        "ALK",
        "AlkylH",
        "ALH",
        "Alkenyl",
        "AEL",
        "AlkenylH",
        "AEH",
        "Alkynyl",
        "AYL",
        "AlkynylH",
        "AYH",
        "Carbocyclic",
        "CBC",
        "CarbocyclicH",
        "CBH",
        "Carbocycloalkyl",
        "CAL",
        "CarbocycloalkylH",
        "CAH",
        "Carbocycloalkenyl",
        "CEL",
        "CarbocycloalkenylH",
        "CEH",
        "Carboaryl",
        "ARY",
        "CarboarylH",
        "ARH",
        "Cyclic",
        "CYC",
        "CyclicH",
        "CYH",
        "Acyclic",
        "ACY",
        "AcyclicH",
        "ACH",
        "Carboacyclic",
        "ABC",
        "CarboacyclicH",
        "ABH",
        "Heteroacyclic",
        "AHC",
        "HeteroacyclicH",
        "AHH",
        "Alkoxy",
        "AOX",
        "AlkoxyH",
        "AOH",
        "Heterocyclic",
        "CHC",
        "HeterocyclicH",
        "CHH",
        "Heteroaryl",
        "HAR",
        "HeteroarylH",
        "HAH",
        "NoCarbonRing",
        "CXX",
        "NoCarbonRingH",
        "CXH",
    ];
    let atom_lines = labels
        .iter()
        .enumerate()
        .map(|(index, label)| format!("M  V30 {} {label} 0 0 0 0", index + 1))
        .collect::<Vec<_>>();
    let atom_refs = atom_lines.iter().map(String::as_str).collect::<Vec<_>>();
    let block =
        v3000_with_outer_and_inner(ZERO_OUTER, &format!("{} 0 0 0 0", labels.len()), &atom_refs);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&block).expect("generic group symbol table")
    else {
        panic!("generic labels must produce a query record");
    };
    for (index, label) in labels.iter().enumerate() {
        let atom = record.query.atom(index).expect("generic query atom");
        assert_eq!(atom.atom().prop("atomLabel"), Some(*label), "{label}");
        assert_eq!(atom.atom().element().atomic_number(), 0, "{label}");
        assert!(!atom.atom().no_implicit(), "{label}");
        let debug = format!("{:?}", atom.predicate());
        assert!(debug.contains("AtomicNumber(0)"), "{label}: {debug}");
    }
}

#[test]
fn v3k_symbols_queries_negative_lookalikes_are_not_special_symbols() {
    // Source matching is case-sensitive and exact. R-label recognition uses
    // the pinned lexicographic condition; these tokens fall outside it.
    for symbol in ["a", "qh", "r7", "R9x", "R999", "group", "Alkylx", "ALKY"] {
        let block = v3000_single_atom(symbol);
        assert!(
            matches!(read_mol_block_detached(&block), Err(SdfReadError::Parse(_))),
            "lookalike {symbol:?} must use the normal strict element path"
        );
    }
}

#[test]
fn v3k_charge_concrete_atoms_preserve_signed_model_boundaries() {
    // RDKit 2026.03.1 ParseV3000AtomProps parses CHG with toInt and calls
    // setFormalCharge for atoms without a query. These values also exercise
    // both inclusive boundaries of the detached i8 charge representation.
    let block = v3000_block(
        &[
            "M  V30 1 C 0 0 0 0 CHG=7",
            "M  V30 2 N 0 0 0 0 CHG=-2",
            "M  V30 3 O 0 0 0 0 CHG=127",
            "M  V30 4 F 0 0 0 0 CHG=-128",
        ],
        4,
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("signed concrete charges")
    else {
        panic!("ordinary element rows must remain concrete");
    };
    let charges = topology
        .atoms
        .iter()
        .map(|atom| atom.formal_charge())
        .collect::<Vec<_>>();
    assert_eq!(charges, [7, -2, 127, -128]);
}

#[test]
fn v3k_charge_leading_plus_is_screened_but_not_converted_for_concrete_and_query_atoms() {
    // RDKit 2026.03.1 FileParserUtils::toInt permits '+' during character
    // screening, but signed std::from_chars does not accept a leading plus.
    // Its initially-zero destination is retained because the error is ignored.
    let concrete = v3000_block(&["M  V30 1 C 0 0 0 0 CHG=+7"], 1);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&concrete).expect("screened concrete leading plus")
    else {
        panic!("ordinary element row must remain concrete");
    };
    assert_eq!(topology.atoms[0].formal_charge(), 0);

    let query = v3000_block(&["M  V30 1 * 0 0 0 0 CHG=+7"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&query).expect("screened query leading plus")
    else {
        panic!("wildcard row must remain a query");
    };
    let atom = record.query.atom(0).expect("query atom");
    assert_eq!(atom.atom().formal_charge(), 0);
    let QueryNode::And(children) = atom.predicate() else {
        panic!("expandQuery must conjoin the formal-charge predicate");
    };
    assert_eq!(children.len(), 2);
    assert!(matches!(
        &children[0],
        QueryNode::Predicate(AtomQueryPredicate::Any)
    ));
    assert!(matches!(
        &children[1],
        QueryNode::Predicate(AtomQueryPredicate::FormalCharge(0))
    ));
}

#[test]
fn v3k_charge_query_atoms_expand_the_query_without_setting_concrete_charge() {
    // ParseV3000AtomProps uses expandQuery(makeAtomFormalChargeQuery(charge))
    // when the atom already has a query, leaving its concrete charge field at
    // the source default.
    let block = v3000_block(&["M  V30 1 * 0 0 0 0 CHG=-3"], 1);
    let MolBlockRecord::Query(record) = read_mol_block_detached(&block).expect("query charge")
    else {
        panic!("wildcard with CHG must remain a query record");
    };
    let atom = record.query.atom(0).expect("query atom");
    assert_eq!(atom.atom().formal_charge(), 0);
    let QueryNode::And(children) = atom.predicate() else {
        panic!("expandQuery must conjoin the formal-charge predicate");
    };
    assert_eq!(children.len(), 2);
    assert!(matches!(
        &children[0],
        QueryNode::Predicate(AtomQueryPredicate::Any)
    ));
    assert!(matches!(
        &children[1],
        QueryNode::Predicate(AtomQueryPredicate::FormalCharge(-3))
    ));
}

#[test]
fn v3k_charge_to_int_empty_invalid_and_source_overflow_behavior() {
    // FileParserUtils::toInt initializes its destination to zero and ignores
    // from_chars' result. Empty text and an out-of-int-range decimal therefore
    // produce zero, while the preceding character screen rejects other bytes.
    for (property, expected) in [
        ("CHG=", 0),
        ("CHG=2147483648", 0),
        ("CHG=--1", 0),
        ("\"CHG=   -7\"", -7),
        ("\"CHG=   \"", 0),
        ("\"CHG=7-\"", 7),
        ("\"CHG=7 8\"", 7),
        ("\"CHG=-2+9\"", -2),
    ] {
        let atom = format!("M  V30 1 C 0 0 0 0 {property}");
        let block = v3000_block(&[&atom], 1);
        let MolBlockRecord::Concrete { topology, .. } =
            read_mol_block_detached(&block).unwrap_or_else(|error| panic!("{property}: {error:?}"))
        else {
            panic!("{property}: concrete record expected");
        };
        assert_eq!(topology.atoms[0].formal_charge(), expected, "{property}");
    }

    for symbol in ["C", "*"] {
        for property in ["CHG=1x", "\"CHG=\t7\""] {
            let atom = format!("M  V30 1 {symbol} 0 0 0 0 {property}");
            let block = v3000_block(&[&atom], 1);
            assert!(matches!(
                read_mol_block_detached(&block),
                Err(SdfReadError::Field {
                    kind: "V3000 atom charge",
                    ..
                })
            ));
        }
    }
}

#[test]
fn v3k_charge_rejects_unrepresentable_concrete_and_query_values_without_narrowing() {
    // RDKit concrete atoms narrow through int8_t while its equality-query
    // target remains int. COSMolKit currently uses a checked i8 model for both
    // and deliberately rejects these rows instead of copying source narrowing
    // or silently narrowing a full-width query target.
    for (symbol, charge) in [("C", "128"), ("C", "-129"), ("*", "128"), ("*", "-129")] {
        let atom = format!("M  V30 1 {symbol} 0 0 0 0 CHG={charge}");
        let block = v3000_block(&[&atom], 1);
        assert!(
            matches!(
                read_mol_block_detached(&block),
                Err(SdfReadError::Unsupported(_))
            ),
            "{symbol} CHG={charge} must not narrow"
        );
    }
}

#[test]
fn v3k_radical_all_four_source_values_preserve_typed_atom_state() {
    // RDKit 2026.03.1 ParseV3000AtomProps treats RAD=0 as a no-op and maps
    // CTAB radical codes 1/2/3 to 2/1/2 radical electrons respectively.
    let block = v3000_block(
        &[
            "M  V30 1 C 0 0 0 0 RAD=0",
            "M  V30 2 N 0 0 0 0 RAD=1",
            "M  V30 3 O 0 0 0 0 RAD=2",
            "M  V30 4 F 0 0 0 0 RAD=3",
        ],
        4,
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("all V3000 RAD codes")
    else {
        panic!("ordinary element rows must remain concrete");
    };
    let radicals = topology
        .atoms
        .iter()
        .map(|atom| atom.radical_electrons())
        .collect::<Vec<_>>();
    assert_eq!(radicals, [0, 2, 1, 2]);
}

#[test]
fn v3k_radical_repeats_zero_noop_and_query_atoms_do_not_gain_query_chemistry() {
    // The source processes properties in order. RAD=0 executes `break`
    // without calling the setter, while later nonzero values overwrite the
    // atom field. Its explicit query FIXME does not add a query predicate.
    let concrete = v3000_block(
        &[
            "M  V30 1 C 0 0 0 0 RAD=1 RAD=0",
            "M  V30 2 N 0 0 0 0 RAD=2 RAD=3",
            "M  V30 3 O 0 0 0 0 RAD=+7",
            "M  V30 4 F 0 0 0 0 RAD=",
        ],
        4,
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&concrete).expect("repeated and no-conversion RAD values")
    else {
        panic!("ordinary element rows must remain concrete");
    };
    let radicals = topology
        .atoms
        .iter()
        .map(|atom| atom.radical_electrons())
        .collect::<Vec<_>>();
    assert_eq!(radicals, [2, 2, 0, 0]);

    let query = v3000_block(&["M  V30 1 * 0 0 0 0 RAD=2"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&query).expect("query atom RAD state")
    else {
        panic!("wildcard row must remain a query");
    };
    let atom = record.query.atom(0).expect("query atom");
    assert_eq!(atom.atom().radical_electrons(), 1);
    assert!(matches!(
        atom.predicate(),
        QueryNode::Predicate(AtomQueryPredicate::Any)
    ));
}

#[test]
fn v3k_radical_invalid_positive_negative_and_screened_text_are_structured_errors() {
    // ParseV3000AtomProps rejects parsed enum values outside 0..=3, while
    // FileParserUtils::toInt rejects bytes outside its screened character set.
    for value in ["-1", "4"] {
        let atom = format!("M  V30 1 C 0 0 0 0 RAD={value}");
        let block = v3000_block(&[&atom], 1);
        assert!(
            matches!(read_mol_block_detached(&block), Err(SdfReadError::Parse(_))),
            "RAD={value} must be an invalid-enum parse error"
        );
    }

    for symbol in ["C", "*"] {
        let atom = format!("M  V30 1 {symbol} 0 0 0 0 RAD=1x");
        let block = v3000_block(&[&atom], 1);
        assert!(matches!(
            read_mol_block_detached(&block),
            Err(SdfReadError::Field {
                kind: "V3000 radical",
                ..
            })
        ));
    }
}

#[test]
fn v3k_mass_integer_fractional_and_integer_conversion_edges_are_source_shaped() {
    // RDKit 2026.03.1 ParseV3000AtomProps tries FileParserUtils::toInt first
    // and only falls back to toDouble/floor after the integer character screen
    // throws. Thus leading plus and signed-int overflow retain toInt's zero,
    // while a decimal point selects the floating fallback.
    let block = v3000_block(
        &[
            "M  V30 1 C 0 0 0 0 MASS=13",
            "M  V30 2 N 0 0 0 0 MASS=13.9",
            "M  V30 3 O 0 0 0 0 MASS=+13",
            "M  V30 4 F 0 0 0 0 MASS=2147483648",
            "M  V30 5 Cl 0 0 0 0 \"MASS=   14.9\"",
            "M  V30 6 Br 0 0 0 0 MASS=1-2",
        ],
        6,
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("source-shaped V3000 MASS conversion")
    else {
        panic!("ordinary element rows must remain concrete");
    };
    let isotopes = topology
        .atoms
        .iter()
        .map(|atom| atom.isotope())
        .collect::<Vec<_>>();
    assert_eq!(
        isotopes,
        [Some(13), Some(13), Some(0), Some(0), Some(14), Some(1)]
    );
}

#[test]
fn v3k_mass_query_atoms_expand_typed_isotope_predicates() {
    // The source calls expandQuery(makeAtomIsotopeQuery(v)) for an existing
    // query atom. Assert the typed conjunction and leave concrete isotope
    // state unchanged instead of relying on a Debug rendering.
    for (property, expected) in [
        ("MASS=13", 13),
        ("MASS=13.9", 13),
        ("MASS=+13", 0),
        ("MASS=2147483648", 0),
        ("MASS=1-2", 1),
    ] {
        let atom_line = format!("M  V30 1 * 0 0 0 0 {property}");
        let block = v3000_block(&[&atom_line], 1);
        let MolBlockRecord::Query(record) =
            read_mol_block_detached(&block).unwrap_or_else(|error| panic!("{property}: {error:?}"))
        else {
            panic!("{property}: wildcard row must remain a query");
        };
        let atom = record.query.atom(0).expect("query atom");
        assert_eq!(atom.atom().isotope(), None, "{property}");
        let QueryNode::And(children) = atom.predicate() else {
            panic!("{property}: isotope expansion must conjoin the predicate");
        };
        assert_eq!(children.len(), 2, "{property}");
        assert!(matches!(
            &children[0],
            QueryNode::Predicate(AtomQueryPredicate::Any)
        ));
        assert!(matches!(
            &children[1],
            QueryNode::Predicate(AtomQueryPredicate::Isotope(value)) if *value == expected
        ));
    }
}

#[test]
fn v3k_mass_negative_and_malformed_values_are_parse_errors() {
    // Negative integers and negative floored fractions reach the source's
    // explicit v < 0 failure. Text rejected by both toInt and toDouble is
    // converted to the same source MASS parse failure.
    for symbol in ["C", "*"] {
        for value in ["-1", "-0.1", "nope", "1x"] {
            let atom_line = format!("M  V30 1 {symbol} 0 0 0 0 MASS={value}");
            let block = v3000_block(&[&atom_line], 1);
            assert!(
                matches!(read_mol_block_detached(&block), Err(SdfReadError::Parse(_))),
                "{symbol} MASS={value} must be a source parse failure"
            );
        }
    }
}

#[test]
fn v3k_mass_unrepresentable_model_values_are_checked_without_narrowing() {
    // Upstream concrete isotope storage is uint16_t but query equality keeps
    // an int target. COSMolKit's concrete and typed-query isotope models are
    // both u16, so representationally valid source integers above 65535 are
    // an explicit current model boundary, not silently narrowed parity.
    for symbol in ["C", "*"] {
        for value in ["65536", "2147483647"] {
            let atom_line = format!("M  V30 1 {symbol} 0 0 0 0 MASS={value}");
            let block = v3000_block(&[&atom_line], 1);
            assert!(
                matches!(
                    read_mol_block_detached(&block),
                    Err(SdfReadError::Unsupported(_))
                ),
                "{symbol} MASS={value} must not narrow into the u16 model"
            );
        }
    }
}

#[test]
fn v3k_atom_cfg_all_source_values_store_only_typed_parity_metadata() {
    // RDKit 2026.03.1 ParseV3000AtomProps treats CFG=0 as a no-op and stores
    // molParity for CFG=1/2/3; stereo perception belongs to later processing.
    let block = v3000_block(
        &[
            "M  V30 1 C 0 0 0 0 CFG=0",
            "M  V30 2 N 0 0 0 0 CFG=1",
            "M  V30 3 O 0 0 0 0 CFG=2",
            "M  V30 4 F 0 0 0 0 CFG=3",
        ],
        4,
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("all V3000 atom CFG values")
    else {
        panic!("ordinary element rows must remain concrete");
    };
    let parities = topology
        .atoms
        .iter()
        .map(|atom| atom.mol_parity())
        .collect::<Vec<_>>();
    assert_eq!(parities, [None, Some(1), Some(2), Some(3)]);

    let query = v3000_block(&["M  V30 1 * 0 0 0 0 CFG=2"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&query).expect("query atom CFG metadata")
    else {
        panic!("wildcard row must remain a query");
    };
    let atom = record.query.atom(0).expect("query atom");
    assert_eq!(atom.atom().mol_parity(), Some(2));
    assert_eq!(atom.atom().prop("molParity"), Some("2"));
    assert!(matches!(
        atom.predicate(),
        QueryNode::Predicate(AtomQueryPredicate::Any)
    ));
}

#[test]
fn v3k_atom_cfg_repeats_and_source_integer_conversion_preserve_order() {
    let block = v3000_block(
        &[
            "M  V30 1 C 0 0 0 0 CFG=1 CFG=0",
            "M  V30 2 N 0 0 0 0 CFG=1 CFG=2",
            "M  V30 3 O 0 0 0 0 CFG=+2",
            "M  V30 4 F 0 0 0 0 CFG=",
            "M  V30 5 Cl 0 0 0 0 CFG=1-2",
        ],
        5,
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("ordered CFG properties")
    else {
        panic!("ordinary element rows must remain concrete");
    };
    let parities = topology
        .atoms
        .iter()
        .map(|atom| atom.mol_parity())
        .collect::<Vec<_>>();
    assert_eq!(parities, [Some(1), Some(2), None, None, Some(1)]);
}

#[test]
fn v3k_atom_cfg_invalid_values_are_structured_errors() {
    for symbol in ["C", "*"] {
        for value in ["-1", "4"] {
            let atom_line = format!("M  V30 1 {symbol} 0 0 0 0 CFG={value}");
            let block = v3000_block(&[&atom_line], 1);
            assert!(matches!(
                read_mol_block_detached(&block),
                Err(SdfReadError::Parse(_))
            ));
        }
        let atom_line = format!("M  V30 1 {symbol} 0 0 0 0 CFG=1x");
        let block = v3000_block(&[&atom_line], 1);
        assert!(matches!(
            read_mol_block_detached(&block),
            Err(SdfReadError::Field {
                kind: "V3000 atom CFG",
                ..
            })
        ));
    }
}

#[test]
fn v3k_hcount_literal_zero_is_the_only_source_noop() {
    // RDKit 2026.03.1 ParseV3000AtomProps guards HCOUNT with the exact
    // textual comparison `val != "0"`; this spelling never converts a
    // concrete atom to a query and never expands an existing query atom.
    let concrete = v3000_block(&["M  V30 1 C 0 0 0 0 HCOUNT=0"], 1);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&concrete).expect("literal-zero concrete HCOUNT")
    else {
        panic!("literal HCOUNT=0 must leave a concrete atom concrete");
    };
    assert_eq!(topology.atoms[0].atomic_number(), 6);

    let wildcard = v3000_block(&["M  V30 1 * 0 0 0 0 HCOUNT=0"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&wildcard).expect("literal-zero wildcard HCOUNT")
    else {
        panic!("wildcard atom must remain a query");
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::predicate(AtomQueryPredicate::Any)
    );
}

#[test]
fn v3k_hcount_alternate_zero_and_nonpositive_values_expand_equality_zero() {
    // Every non-literal-zero spelling enters FileParserUtils::toInt. The
    // source normalizes -1 explicitly and sends every other nonpositive
    // result, including -2 and numeric zero text, to the same equality query.
    for value in ["00", "-1", "-2"] {
        for (symbol, base) in [
            ("C", AtomQueryPredicate::AtomicNumber(6)),
            ("*", AtomQueryPredicate::Any),
        ] {
            let atom_line = format!("M  V30 1 {symbol} 0 0 0 0 HCOUNT={value}");
            let block = v3000_block(&[&atom_line], 1);
            let MolBlockRecord::Query(record) = read_mol_block_detached(&block)
                .unwrap_or_else(|error| panic!("{symbol} HCOUNT={value}: {error:?}"))
            else {
                panic!("{symbol} HCOUNT={value} must produce query topology");
            };
            assert_eq!(
                record.query.atoms()[0].predicate(),
                &QueryNode::and(vec![
                    QueryNode::predicate(base.clone()),
                    QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCount(0)),
                ]),
                "{symbol} HCOUNT={value}"
            );
        }
    }
}

#[test]
fn v3k_hcount_positive_values_expand_less_equal_with_checked_model_width() {
    // The source constructs ATOM_LESSEQUAL_QUERY with the parsed C++ int.
    // COSMolKit preserves representable u8 targets and fails closed at 256
    // instead of silently narrowing the independent source query target.
    for value in [1_u8, 255] {
        for (symbol, base) in [
            ("C", AtomQueryPredicate::AtomicNumber(6)),
            ("*", AtomQueryPredicate::Any),
        ] {
            let atom_line = format!("M  V30 1 {symbol} 0 0 0 0 HCOUNT={value}");
            let block = v3000_block(&[&atom_line], 1);
            let MolBlockRecord::Query(record) = read_mol_block_detached(&block)
                .unwrap_or_else(|error| panic!("{symbol} HCOUNT={value}: {error:?}"))
            else {
                panic!("{symbol} HCOUNT={value} must produce query topology");
            };
            assert_eq!(
                record.query.atoms()[0].predicate(),
                &QueryNode::and(vec![
                    QueryNode::predicate(base.clone()),
                    QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCountLessEqual(value),),
                ]),
                "{symbol} HCOUNT={value}"
            );
        }
    }

    for symbol in ["C", "*"] {
        let atom_line = format!("M  V30 1 {symbol} 0 0 0 0 HCOUNT=256");
        let block = v3000_block(&[&atom_line], 1);
        assert!(matches!(
            read_mol_block_detached(&block),
            Err(SdfReadError::Unsupported(message))
                if message.contains("outside the detached implicit-hydrogen-count query model")
        ));
    }
}

#[test]
fn v3k_hcount_concrete_conversion_skips_zero_isotope_state() {
    // RDKit 2026.03.1 QueryAtom(const Atom&) guards isotope expansion with
    // `if (other.getIsotope())`. MASS=0 stores zero on the concrete Atom, but
    // conversion for HCOUNT must therefore retain only atomic number before
    // appending the H-count predicate. MASS=+13 follows source toInt
    // no-conversion semantics and stores the same zero value.
    let expected = QueryNode::and(vec![
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCountLessEqual(1)),
    ]);
    for properties in ["MASS=0 HCOUNT=1", "MASS=+13 HCOUNT=1"] {
        let atom_line = format!("M  V30 1 C 0 0 0 0 {properties}");
        let block = v3000_block(&[&atom_line], 1);
        let MolBlockRecord::Query(record) = read_mol_block_detached(&block)
            .unwrap_or_else(|error| panic!("{properties}: {error:?}"))
        else {
            panic!("{properties} must produce query topology");
        };
        assert_eq!(
            record.query.atoms()[0].predicate(),
            &expected,
            "{properties}"
        );
    }
}

#[test]
fn v3k_hcount_concrete_conversion_retains_nonzero_constructor_predicates() {
    // QueryAtom(const Atom&) expands nonzero isotope, formal charge, and
    // radical predicates in that order before ParseV3000AtomProps appends the
    // HCOUNT predicate. Assert the complete typed tree so the zero-isotope
    // correction cannot discard meaningful constructor state.
    let cases = [
        (
            "MASS=13 HCOUNT=1",
            QueryNode::and(vec![
                QueryNode::and(vec![
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
                ]),
                QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCountLessEqual(1)),
            ]),
        ),
        (
            "CHG=1 MASS=13 RAD=2 HCOUNT=1",
            QueryNode::and(vec![
                QueryNode::and(vec![
                    QueryNode::and(vec![
                        QueryNode::and(vec![
                            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                            QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
                        ]),
                        QueryNode::predicate(AtomQueryPredicate::FormalCharge(1)),
                    ]),
                    QueryNode::predicate(AtomQueryPredicate::NumRadicalElectrons(1)),
                ]),
                QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCountLessEqual(1)),
            ]),
        ),
    ];
    for (properties, expected) in cases {
        let atom_line = format!("M  V30 1 C 0 0 0 0 {properties}");
        let block = v3000_block(&[&atom_line], 1);
        let MolBlockRecord::Query(record) = read_mol_block_detached(&block)
            .unwrap_or_else(|error| panic!("{properties}: {error:?}"))
        else {
            panic!("{properties} must produce query topology");
        };
        assert_eq!(
            record.query.atoms()[0].predicate(),
            &expected,
            "{properties}"
        );
    }
}

#[test]
fn v3k_hcount_then_mass_zero_retains_source_query_expansion() {
    // Property order is observable: after HCOUNT has replaced the concrete
    // atom with a QueryAtom, the MASS branch directly expands Isotope(0).
    // This is distinct from constructor conversion of pre-existing zero state.
    let block = v3000_block(&["M  V30 1 C 0 0 0 0 HCOUNT=1 MASS=0"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&block).expect("HCOUNT before MASS=0")
    else {
        panic!("HCOUNT=1 MASS=0 must produce query topology");
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::and(vec![
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCountLessEqual(1)),
            ]),
            QueryNode::predicate(AtomQueryPredicate::Isotope(0)),
        ])
    );
}

fn v3000_single_atom_x(x: &str) -> String {
    let atom = format!("M  V30 1 C {x} 0 0 0");
    v3000_with_outer_and_inner(ZERO_OUTER, "1 0 0 0 0", &[&atom])
}

#[test]
fn v3k_atom_numbers_hexadecimal_forms_exact() {
    // A non-zero z selects 3D storage so all three axes are observable.
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "2 0 0 0 0",
        &[
            "M  V30 1 C 0x1p+1 0x1.8p+1 0x.8p0 0",
            "M  V30 2 O 0X1P+1 0x10 -0x1p+1 0",
        ],
    );
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("hexadecimal coordinates")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(coordinates.conformers_3d.len(), 1);
    let rows = coordinates.conformers_3d[0].coordinates();
    assert_eq!(rows[0], [2.0, 3.0, 0.5]);
    assert_eq!(rows[1], [2.0, 16.0, -2.0]);
}

#[test]
fn v3k_atom_numbers_hexadecimal_partial_and_prefix() {
    let cases = [
        ("0x10", 16.0),
        ("0x1", 1.0),
        ("0x1p", 1.0),
        ("0x1p-1", 0.5),
        ("0x.8p0", 0.5),
        ("0x", 0.0),
        ("0xg", 0.0),
    ];
    for (token, expected) in cases {
        let block = v3000_single_atom_x(token);
        let MolBlockRecord::Concrete { coordinates, .. } =
            read_mol_block_detached(&block).unwrap_or_else(|error| panic!("{token}: {error:?}"))
        else {
            panic!("{token}: atom record must be concrete");
        };
        assert_eq!(
            coordinates.conformers_2d[0].coordinates()[0][0],
            expected,
            "{token}"
        );
    }
}

#[test]
fn v3k_atom_numbers_signed_zero_and_subnormal_boundaries() {
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "4 0 0 0 0",
        &[
            "M  V30 1 C -0 0 0 0",
            "M  V30 2 O 0x1p-1074 0 0 0",
            "M  V30 3 N 0x1.8p-1075 0 0 0",
            "M  V30 4 F 0x1p-1075 0 0 0",
        ],
    );
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("subnormal boundaries")
    else {
        panic!("atom record must be concrete");
    };
    let rows = coordinates.conformers_2d[0].coordinates();
    assert_eq!(rows[0][0], 0.0);
    assert!(rows[0][0].is_sign_negative());
    // 2^-1074 is the smallest positive subnormal.
    assert_eq!(rows[1][0], f64::from_bits(1));
    // 1.5 * 2^-1075 rounds up to the smallest subnormal.
    assert_eq!(rows[2][0], f64::from_bits(1));
    // Exactly 2^-1075 is a tie and rounds to even (zero).
    assert_eq!(rows[3][0], 0.0);
}

#[test]
fn v3k_atom_numbers_overflow_and_nonfinite_rejected() {
    for token in ["0x1p+2000", "0x1p+1024", "1e400"] {
        let block = v3000_single_atom_x(token);
        assert!(
            matches!(
                read_mol_block_detached(&block),
                Err(SdfReadError::Coordinates(_))
            ),
            "overflowing coordinate {token:?} must fail the finite-coordinate boundary"
        );
    }
    // Gradual underflow to zero is a finite coordinate and is accepted.
    let block = v3000_single_atom_x("0x1p-2000");
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("underflow to zero")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(coordinates.conformers_2d[0].coordinates()[0][0], 0.0);
}

#[test]
fn v3k_symbols_elements_multibyte_symbol_does_not_panic() {
    // U+1F600 is four bytes long, so byte index 3 is not a character boundary.
    let block = v3000_single_atom("\u{1F600}");
    assert!(
        matches!(read_mol_block_detached(&block), Err(SdfReadError::Parse(_))),
        "strict unknown multibyte symbol must be a structured parse error"
    );
    let params = cosmolkit_io::MolBlockReadParams {
        strict_parsing: false,
        ..cosmolkit_io::MolBlockReadParams::default()
    };
    let MolBlockRecord::Concrete { topology, .. } =
        cosmolkit_io::read_mol_block_detached_with_params(&block, params)
            .expect("non-strict multibyte dummy")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(topology.atoms[0].element().atomic_number(), 0);
    assert_eq!(topology.atoms[0].prop("dummyLabel"), Some("\u{1F600}"));
}

#[test]
fn v3k_symbols_elements_multibyte_prefixes_do_not_panic() {
    // Each token has text before the four-byte code point, so the old
    // `symbol[..3]` prefix check split a code point; all must be structured.
    for token in ["N\u{1F600}", "NO\u{1F600}", "NOT\u{1F600}"] {
        let block = v3000_single_atom(token);
        assert!(
            read_mol_block_detached(&block).is_err(),
            "{token:?} must produce a structured error, not a panic"
        );
    }
}

#[test]
fn v3k_symbols_lists_multibyte_entry_does_not_panic() {
    let block = v3000_single_atom("\"[\u{1F600},C]\"");
    assert!(
        matches!(read_mol_block_detached(&block), Err(SdfReadError::Parse(_))),
        "an unknown multibyte list entry must be a structured parse error"
    );
}

#[test]
fn v3k_atom_numbers_atof_port_counterexamples() {
    // 0x10p followed by fifty 9s overflows to +inf and must fail the existing
    // finite-coordinate boundary (previously a silent 0.0).
    let overflow = format!("0x10p{}", "9".repeat(50));
    let block = v3000_single_atom_x(&overflow);
    assert!(
        matches!(
            read_mol_block_detached(&block),
            Err(SdfReadError::Coordinates(_))
        ),
        "huge positive hex exponent must be rejected as non-finite"
    );

    // 0x0.1p- followed by fifty 9s underflows to +0.0 and is accepted as a
    // finite coordinate (previously a spurious non-finite error).
    let underflow = format!("0x0.1p-{}", "9".repeat(50));
    let block = v3000_single_atom_x(&underflow);
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("huge negative hex exponent underflows to zero")
    else {
        panic!("atom record must be concrete");
    };
    let value = coordinates.conformers_2d[0].coordinates()[0][0];
    assert_eq!(value, 0.0);
    assert!(!value.is_sign_negative());
}

#[test]
fn v3k_atom_numbers_atof_port_finite_bits_and_signed_zero() {
    let cases = [
        ("0x1.0000000000001p0", f64::from_bits(0x3ff0_0000_0000_0001)),
        ("0x1p-1074", f64::from_bits(0x0000_0000_0000_0001)),
        ("0x1.8p+1", 3.0),
    ];
    for (token, expected) in cases {
        let block = v3000_single_atom_x(token);
        let MolBlockRecord::Concrete { coordinates, .. } =
            read_mol_block_detached(&block).unwrap_or_else(|error| panic!("{token}: {error:?}"))
        else {
            panic!("{token}: atom record must be concrete");
        };
        assert_eq!(
            coordinates.conformers_2d[0].coordinates()[0][0].to_bits(),
            expected.to_bits(),
            "{token}"
        );
    }

    // Negative zero stays negative through the detached reader.
    let block = v3000_single_atom_x("-0");
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("negative zero coordinate")
    else {
        panic!("atom record must be concrete");
    };
    let value = coordinates.conformers_2d[0].coordinates()[0][0];
    assert_eq!(value.to_bits(), 0x8000_0000_0000_0000);
}

#[test]
fn v3k_atom_numbers_atof_port_suffix_and_malformed_prefixes() {
    let cases = [
        ("12xyz", 12.0),
        ("0x", 0.0),
        ("abc", 0.0),
        ("1e", 1.0),
        ("  3.5", 3.5),
        ("1,5", 1.0),
    ];
    for (token, expected) in cases {
        let block = v3000_single_atom_x(token);
        let MolBlockRecord::Concrete { coordinates, .. } =
            read_mol_block_detached(&block).unwrap_or_else(|error| panic!("{token}: {error:?}"))
        else {
            panic!("{token}: atom record must be concrete");
        };
        assert_eq!(
            coordinates.conformers_2d[0].coordinates()[0][0],
            expected,
            "{token}"
        );
    }
}

#[test]
fn v3k_atom_numbers_atof_port_nan_rejected_and_alignment_unchanged() {
    for token in ["nan", "-nan", "0x1p+1024", "1e400"] {
        let block = v3000_single_atom_x(token);
        assert!(
            matches!(
                read_mol_block_detached(&block),
                Err(SdfReadError::Coordinates(_))
            ),
            "non-finite coordinate {token:?} must fail the finite-coordinate boundary"
        );
    }

    // Row alignment and atom ordering are unchanged by the conversion swap.
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "3 0 0 0 0",
        &[
            "M  V30 1 C 0x1p+1 1.0 0 0",
            "M  V30 2 O 2.0 0X1P+1 0 0",
            "M  V30 3 N 4.0 5.0 0 0",
        ],
    );
    let MolBlockRecord::Concrete {
        topology,
        coordinates,
        ..
    } = read_mol_block_detached(&block).expect("conversion row alignment")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(topology.atoms.len(), 3);
    let rows = coordinates.conformers_2d[0].coordinates();
    assert_eq!(rows.len(), 3);
    assert_eq!(rows[0], [2.0, 1.0]);
    assert_eq!(rows[1], [2.0, 2.0]);
    assert_eq!(rows[2], [4.0, 5.0]);
}

#[test]
fn v3k_atom_numbers_atof_port_musl_long_hex_source_value() {
    // Retained excluded-scope discrepancy, NOT RDKit acceptance evidence:
    // this musl-derived implementation rounds a long hexadecimal significand
    // one ULP below glibc/RDKit; see IO-atof-port.md. Narrowing the coordinate
    // contract must not erase this regression or claim all-input equivalence.
    let block = v3000_single_atom_x("0x1814D89994224121");
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("long hex significand")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(
        coordinates.conformers_2d[0].coordinates()[0][0].to_bits(),
        0x43b8_14d8_9994_2242
    );
}

#[test]
fn v3k_atom_numbers_coordinate_contract_decimal_and_missing_fields() {
    // Pinned ParseV3000AtomBlock uses raw atof (unlike V2000 toDouble), but
    // requires each coordinate/map token to exist before conversion.
    for (token, expected) in [
        ("+1.25e2", 125.0_f64),
        ("-1.25E-2", -0.0125),
        ("-0.000e+12", -0.0),
        ("1.25suffix", 1.25),
        ("abc", 0.0),
    ] {
        let MolBlockRecord::Concrete { coordinates, .. } =
            read_mol_block_detached(&v3000_single_atom_x(token)).expect("V3000 prefix")
        else {
            panic!("expected concrete molecule");
        };
        assert_eq!(
            coordinates.conformers_2d[0].coordinates()[0][0].to_bits(),
            expected.to_bits(),
            "{token}"
        );
    }
    for row in [
        "M  V30 1 C",
        "M  V30 1 C 1",
        "M  V30 1 C 1 2",
        "M  V30 1 C 1 2 3",
    ] {
        let block = v3000_with_outer_and_inner(ZERO_OUTER, "1 0 0 0 0", &[row]);
        let error = read_mol_block_detached(&block).expect_err("missing required token");
        assert!(matches!(error, SdfReadError::Parse(_)), "{row}: {error:?}");
        assert!(error.to_string().contains("Bad atom line"), "{error:?}");
    }
}
