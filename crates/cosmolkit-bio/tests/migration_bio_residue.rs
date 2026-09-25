use cosmolkit_bio::{
    AltLocLabel, AtomName, AtomSourceIds, ChainSourceIds, EntitySourceIds, PdbAtomSerial,
    PdbChainId, PdbSeqId, ResidueInfoKind, ResidueName, ResidueSequenceError, ResidueSourceIds,
    UNKNOWN_TABULATED_RESIDUE_INDEX, expand_one_letter, expand_one_letter_sequence,
    find_residue_info, find_residue_info_index, residue_code, residue_info, residue_info_checked,
};

#[derive(Debug)]
struct SourceResidue {
    name: String,
    kind: ResidueInfoKind,
    linking_type: u8,
    one_letter_code: char,
    hydrogen_count: u8,
    weight: f32,
}

fn residue_fixture_rows() -> Vec<SourceResidue> {
    let text = include_str!("../../../testdata/bio/fixtures/gemmi_residues/residues.tsv");
    let mut lines = text.lines();
    assert_eq!(
        lines.next(),
        Some("index\tname\tkind\tlinking_type\tone_letter_code\thydrogen_count\tweight_f32")
    );
    let rows: Vec<_> = lines
        .enumerate()
        .map(|(index, line)| {
            // Spaces are data: several source one-letter codes are a literal space.
            let fields: Vec<_> = line.split('\t').collect();
            assert_eq!(fields.len(), 7, "fixture row {index}");
            assert_eq!(fields[0].parse::<usize>().unwrap(), index);
            let kind = match fields[2] {
                "UNKNOWN" => ResidueInfoKind::Unknown,
                "AA" => ResidueInfoKind::Aa,
                "AAD" => ResidueInfoKind::Aad,
                "PAA" => ResidueInfoKind::Paa,
                "MAA" => ResidueInfoKind::Maa,
                "RNA" => ResidueInfoKind::Rna,
                "DNA" => ResidueInfoKind::Dna,
                "BUF" => ResidueInfoKind::Buf,
                "HOH" => ResidueInfoKind::Hoh,
                "PYR" => ResidueInfoKind::Pyr,
                "KET" => ResidueInfoKind::Ket,
                "ELS" => ResidueInfoKind::Els,
                other => panic!("unknown fixture kind {other}"),
            };
            let mut code = fields[4].chars();
            let one_letter_code = code.next().expect("one-letter code");
            assert!(code.next().is_none());
            SourceResidue {
                name: fields[1].to_owned(),
                kind,
                linking_type: fields[3].parse().expect("linking type"),
                one_letter_code,
                hydrogen_count: fields[5].parse().expect("hydrogen count"),
                weight: fields[6].parse().expect("f32 weight"),
            }
        })
        .collect();
    assert_eq!(
        rows.len(),
        368,
        "complete fixed table including unknown row"
    );
    rows
}

#[test]
fn bio_residue_all_368_rows_match_fixed_fixture() {
    let source = residue_fixture_rows();
    assert_eq!(source.len(), 368);
    for (index, expected) in source.iter().enumerate() {
        let actual = residue_info(index);
        assert_eq!(actual.code.as_u16() as usize, index, "code at row {index}");
        assert_eq!(actual.name, expected.name, "name at row {index}");
        assert_eq!(actual.kind, expected.kind, "kind at row {index}");
        assert_eq!(
            actual.linking_type, expected.linking_type,
            "linking type at row {index}"
        );
        assert_eq!(
            actual.one_letter_code, expected.one_letter_code,
            "one-letter code at row {index}"
        );
        assert_eq!(
            actual.hydrogen_count, expected.hydrogen_count,
            "hydrogen count at row {index}"
        );
        assert_eq!(
            actual.weight.to_bits(),
            expected.weight.to_bits(),
            "weight bits at row {index}"
        );
        assert_eq!(residue_info_checked(index), Some(actual));
        assert_eq!(
            actual.found(),
            expected.kind != ResidueInfoKind::Unknown,
            "found at row {index}"
        );
        assert_eq!(actual.is_water(), expected.kind == ResidueInfoKind::Hoh);
        assert_eq!(actual.is_dna(), expected.kind == ResidueInfoKind::Dna);
        assert_eq!(actual.is_rna(), expected.kind == ResidueInfoKind::Rna);
        assert_eq!(
            actual.is_nucleic_acid(),
            matches!(expected.kind, ResidueInfoKind::Dna | ResidueInfoKind::Rna)
        );
        assert_eq!(
            actual.is_amino_acid(),
            matches!(
                expected.kind,
                ResidueInfoKind::Aa
                    | ResidueInfoKind::Aad
                    | ResidueInfoKind::Paa
                    | ResidueInfoKind::Maa
            )
        );
        assert_eq!(
            actual.is_buffer_or_water(),
            matches!(expected.kind, ResidueInfoKind::Buf | ResidueInfoKind::Hoh)
        );
        assert_eq!(
            actual.is_standard(),
            (expected.one_letter_code as u32 & 0x20) == 0
        );
        assert_eq!(
            actual.fasta_code(),
            if actual.is_standard() {
                expected.one_letter_code
            } else {
                'X'
            }
        );
        assert_eq!(actual.is_peptide_linking(), expected.linking_type & 1 != 0);
        assert_eq!(actual.is_na_linking(), expected.linking_type & 2 != 0);
    }
    assert_eq!(residue_info_checked(368), None);
    assert_eq!(residue_info_checked(usize::MAX), None);
}

#[test]
fn bio_residue_lookup_preserves_source_lengths_aliases_and_case_rules() {
    let source = residue_fixture_rows();
    for (index, row) in source
        .iter()
        .enumerate()
        .take(UNKNOWN_TABULATED_RESIDUE_INDEX)
    {
        assert_eq!(
            find_residue_info_index(&row.name),
            index,
            "lookup {}",
            row.name
        );
        assert_eq!(find_residue_info(&row.name), residue_info(index));
        assert_eq!(residue_code(&row.name).as_u16() as usize, index);
    }
    for (alias, index) in [("TRY", 23), ("WAT", 154), ("H2O", 154)] {
        assert_eq!(find_residue_info_index(alias), index);
    }
    for spelling in ["ala", "aLa", "trY", "h2o", "0td"] {
        assert_ne!(
            find_residue_info_index(spelling),
            UNKNOWN_TABULATED_RESIDUE_INDEX
        );
    }
    assert_eq!(find_residue_info_index("a"), 327);
    assert_eq!(find_residue_info_index("k"), 334);
    for (spelling, index) in [
        ("DA", 335),
        ("+A", 335),
        ("DC", 336),
        ("+U", 340),
        ("DN", 341),
        ("AG", 342),
        ("ZN", 366),
    ] {
        assert_eq!(find_residue_info_index(spelling), index);
    }
    for unknown in ["", "AAAA", "Da", "da", "+a", "ag", "J", "1", "?", "🧬"] {
        assert_eq!(
            find_residue_info_index(unknown),
            UNKNOWN_TABULATED_RESIDUE_INDEX,
            "unexpected known lookup for {unknown:?}"
        );
    }
}

#[test]
fn bio_residue_one_letter_expansion_covers_aa_dna_rna_and_invalid_inputs() {
    let aa = [
        ('A', "ALA"),
        ('B', "ASX"),
        ('C', "CYS"),
        ('D', "ASP"),
        ('E', "GLU"),
        ('F', "PHE"),
        ('G', "GLY"),
        ('H', "HIS"),
        ('I', "ILE"),
        ('K', "LYS"),
        ('L', "LEU"),
        ('M', "MET"),
        ('N', "ASN"),
        ('O', "PYL"),
        ('P', "PRO"),
        ('Q', "GLN"),
        ('R', "ARG"),
        ('S', "SER"),
        ('T', "THR"),
        ('U', "SEC"),
        ('V', "VAL"),
        ('W', "TRP"),
        ('X', "UNK"),
        ('Y', "TYR"),
        ('Z', "GLX"),
    ];
    for (letter, name) in aa {
        assert_eq!(expand_one_letter(letter, ResidueInfoKind::Aa), Some(name));
        assert_eq!(
            expand_one_letter(letter.to_ascii_lowercase(), ResidueInfoKind::Aa),
            Some(name)
        );
    }
    assert_eq!(expand_one_letter('J', ResidueInfoKind::Aa), None);

    let dna = [
        ('A', "DA"),
        ('C', "DC"),
        ('G', "DG"),
        ('I', "DI"),
        ('N', "DN"),
        ('T', "DT"),
        ('U', "DU"),
    ];
    let rna = [
        ('A', "A"),
        ('C', "C"),
        ('G', "G"),
        ('I', "I"),
        ('N', "N"),
        ('U', "U"),
    ];
    for (letter, name) in dna {
        assert_eq!(expand_one_letter(letter, ResidueInfoKind::Dna), Some(name));
        assert_eq!(
            expand_one_letter(letter.to_ascii_lowercase(), ResidueInfoKind::Dna),
            Some(name)
        );
    }
    for (letter, name) in rna {
        assert_eq!(expand_one_letter(letter, ResidueInfoKind::Rna), Some(name));
        assert_eq!(
            expand_one_letter(letter.to_ascii_lowercase(), ResidueInfoKind::Rna),
            Some(name)
        );
    }
    assert_eq!(expand_one_letter('T', ResidueInfoKind::Rna), None);
    assert_eq!(expand_one_letter('B', ResidueInfoKind::Dna), None);
    assert_eq!(expand_one_letter('A', ResidueInfoKind::Unknown), None);
    assert_eq!(expand_one_letter('é', ResidueInfoKind::Aa), None);
}

#[test]
fn bio_residue_sequence_expansion_preserves_source_order_whitespace_and_errors() {
    assert_eq!(
        expand_one_letter_sequence(" A\tC\nG\u{000b}T\u{000c}(MSE)() ", ResidueInfoKind::Aa)
            .unwrap(),
        ["ALA", "CYS", "GLY", "THR", "MSE", ""]
    );
    assert_eq!(
        expand_one_letter_sequence("ac", ResidueInfoKind::Dna).unwrap(),
        ["DA", "DC"]
    );
    assert_eq!(
        expand_one_letter_sequence("au", ResidueInfoKind::Rna).unwrap(),
        ["A", "U"]
    );
    assert_eq!(
        expand_one_letter_sequence("A(MSE", ResidueInfoKind::Aa),
        Err(ResidueSequenceError::UnmatchedParenthesis)
    );
    assert_eq!(
        expand_one_letter_sequence("AJ", ResidueInfoKind::Aa),
        Err(ResidueSequenceError::UnexpectedLetter {
            kind: "peptide",
            letter: 'J',
            source_code: 74,
        })
    );
    assert_eq!(
        expand_one_letter_sequence("é", ResidueInfoKind::Aa),
        Err(ResidueSequenceError::UnexpectedLetter {
            kind: "peptide",
            letter: 'Ã',
            source_code: -61,
        })
    );
}

#[test]
fn bio_residue_source_identifiers_preserve_bytes_and_independent_namespaces() {
    let atom_name = AtomName::from_ascii(*b" CA ").unwrap();
    assert_eq!(atom_name.as_bytes(), b" CA ");
    assert_eq!(atom_name.as_str(), " CA ");
    assert!(AtomName::from_ascii([0xff, b'A', b' ', b' ']).is_none());

    let empty = PdbChainId::from_ascii(b"").unwrap();
    let auth = PdbChainId::from_ascii(b"A").unwrap();
    let label = PdbChainId::from_ascii(b"L123").unwrap();
    assert_eq!(empty.as_bytes(), b"");
    assert_eq!(auth.as_str(), "A");
    assert_eq!(label.as_str(), "L123");
    assert!(PdbChainId::from_ascii(b"ABCDE").is_none());
    assert!(PdbChainId::from_ascii(&[0xff]).is_none());

    let chain = ChainSourceIds::new(Some(auth), Some("label_asym_long".to_owned()));
    assert_eq!(chain.auth_chain_id(), Some(auth));
    assert_eq!(chain.label_asym_id(), Some("label_asym_long"));

    let seq = PdbSeqId::new(-12, Some(b'B'));
    let residue = ResidueSourceIds::new(
        Some(seq),
        Some(42),
        Some(*b"SEG "),
        Some("subchain-long-source-id".to_owned()),
        Some("entity-long-source-id".to_owned()),
    )
    .unwrap();
    assert_eq!(residue.seq_id(), Some(seq));
    assert_eq!(seq.seq_num(), -12);
    assert_eq!(seq.ins_code(), Some(b'B'));
    assert_eq!(residue.label_seq_id(), Some(42));
    assert_eq!(residue.segment_id(), Some(b"SEG "));
    assert_eq!(residue.subchain_id(), Some("subchain-long-source-id"));
    assert_eq!(residue.label_entity_id(), Some("entity-long-source-id"));
    assert!(ResidueSourceIds::new(None, None, Some([0xff; 4]), None, None).is_none());

    let serial = PdbAtomSerial::new(-7);
    assert_eq!(serial.value(), -7);
    assert_eq!(AtomSourceIds::new(Some(serial)).serial(), Some(serial));
    assert_eq!(AltLocLabel::new(0).value(), 0);
    assert_eq!(AltLocLabel::new(b'A').value(), b'A');

    let residue_name = ResidueName::from_ascii(b"MSE ").unwrap();
    assert_eq!(residue_name.as_bytes(), b"MSE ");
    assert_eq!(residue_name.as_str(), "MSE ");
    assert!(ResidueName::from_ascii(b"ABCDE").is_none());
    assert!(ResidueName::from_ascii(&[0xff]).is_none());

    let entity = EntitySourceIds::new("entity-123".to_owned());
    assert_eq!(entity.source_entity_id(), "entity-123");
}

#[test]
fn bio_residue_unknown_row_and_out_of_range_contract_are_distinct() {
    let unknown = residue_info(UNKNOWN_TABULATED_RESIDUE_INDEX);
    assert_eq!(unknown.name, "");
    assert_eq!(unknown.kind, ResidueInfoKind::Unknown);
    assert_eq!(unknown.weight.to_bits(), 0.0f32.to_bits());
    assert_eq!(find_residue_info("not-a-residue"), unknown);
    assert_eq!(
        residue_info_checked(UNKNOWN_TABULATED_RESIDUE_INDEX),
        Some(unknown)
    );
    assert_eq!(
        residue_info_checked(UNKNOWN_TABULATED_RESIDUE_INDEX + 1),
        None
    );
}
