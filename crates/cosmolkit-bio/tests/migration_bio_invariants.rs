use cosmolkit_bio::{
    AtomName, AtomSourceIds, BioAssembly, BioAssemblyGenerator, BioAssemblySpecialKind, BioAtomId,
    BioAtomRow, BioCalcFlag, BioChainId, BioChainRow, BioCoordinateBlock, BioCoordinateFormat,
    BioEntityId, BioEntityRow, BioModelId, BioModelRow, BioResidueId, BioResidueRow, BioRowSpan,
    BioSiftsUnpResidue, BioStructure, BioStructureError, BioStructureParts, ChainKind,
    ChainSourceIds, EntityKind, EntitySourceIds, PdbChainId, PolymerKind, ResidueInfoKind,
    ResidueName, ResidueSourceIds,
};
use cosmolkit_types::Element;

fn span<I>(start: u32, len: u32) -> BioRowSpan<I> {
    BioRowSpan::new(start, len).unwrap()
}

fn atom_name() -> AtomName {
    AtomName::from_ascii(*b" CA ").unwrap()
}

fn residue_name() -> ResidueName {
    ResidueName::from_ascii(b"ALA").unwrap()
}

fn entity(subchains: &[&str]) -> BioEntityRow {
    BioEntityRow::new(
        EntityKind::Polymer,
        PolymerKind::PeptideL,
        false,
        Vec::new(),
        Vec::new(),
        Vec::new(),
        subchains.iter().map(|value| (*value).to_owned()).collect(),
        EntitySourceIds::new("1".to_owned()),
    )
}

fn atom(residue_id: u32) -> BioAtomRow {
    BioAtomRow::new(
        BioResidueId::new(residue_id),
        atom_name(),
        Element::C,
        None,
        None,
        0,
        BioCalcFlag::NotSet,
        1.0,
        20.0,
        [0.0; 6],
        -1,
        0.0,
        AtomSourceIds::new(None),
    )
}

fn residue(
    chain_id: u32,
    atom_span: BioRowSpan<BioAtomId>,
    entity_id: Option<BioEntityId>,
    subchain: Option<&str>,
) -> BioResidueRow {
    BioResidueRow::new(
        BioChainId::new(chain_id),
        atom_span,
        residue_name(),
        ResidueInfoKind::Aa,
        EntityKind::Polymer,
        entity_id,
        None,
        ResidueSourceIds::new(None, None, None, subchain.map(str::to_owned), None).unwrap(),
        BioSiftsUnpResidue::default(),
    )
}

fn chain(
    model_id: u32,
    residue_span: BioRowSpan<BioResidueId>,
    entity_id: Option<BioEntityId>,
    auth: Option<&[u8]>,
    label: Option<&str>,
) -> BioChainRow {
    BioChainRow::new(
        BioModelId::new(model_id),
        entity_id,
        residue_span,
        ChainKind::Protein,
        ChainSourceIds::new(
            auth.map(|value| PdbChainId::from_ascii(value).unwrap()),
            label.map(str::to_owned),
        ),
    )
}

fn one_atom_parts() -> BioStructureParts {
    BioStructureParts {
        input_format: BioCoordinateFormat::Mmcif,
        models: vec![BioModelRow::new(span(0, 1), Some(1))],
        chains: vec![chain(
            0,
            span(0, 1),
            Some(BioEntityId::new(0)),
            Some(b"A"),
            Some("LONG_LABEL"),
        )],
        residues: vec![residue(
            0,
            span(0, 1),
            Some(BioEntityId::new(0)),
            Some("LONG_LABEL"),
        )],
        atoms: vec![atom(0)],
        entities: vec![entity(&["LONG_LABEL"])],
        connections: Vec::new(),
        cispeps: Vec::new(),
        mod_residues: Vec::new(),
        helices: Vec::new(),
        sheets: Vec::new(),
        metadata: Default::default(),
        source_state: Default::default(),
        coordinates: BioCoordinateBlock::new(vec![[1.0, 2.0, 3.0]]),
        crystal: None,
        ncs_operators: Vec::new(),
        assemblies: Vec::new(),
    }
}

fn assembly(chains: &[&str], subchains: &[&str]) -> BioAssembly {
    BioAssembly::new(
        "1".to_owned(),
        true,
        false,
        BioAssemblySpecialKind::NotApplicable,
        1,
        String::new(),
        String::new(),
        0.0,
        0.0,
        0.0,
        vec![BioAssemblyGenerator::new(
            chains.iter().map(|value| (*value).to_owned()).collect(),
            subchains.iter().map(|value| (*value).to_owned()).collect(),
            Vec::new(),
        )],
    )
}

fn assert_rejected_without_mutation(parts: &BioStructureParts, expected: BioStructureError) {
    let before = parts.clone();
    assert_eq!(BioStructure::validate_parts(parts), Err(expected.clone()));
    assert_eq!(BioStructure::validate_parts(parts), Err(expected));
    assert_eq!(*parts, before);
}

#[test]
fn bio_invariants_accept_consecutive_empty_parents_and_validate_all_entrypoints() {
    let parts = BioStructureParts {
        input_format: BioCoordinateFormat::Unknown,
        models: vec![
            BioModelRow::new(span(0, 0), Some(1)),
            BioModelRow::new(span(0, 0), Some(2)),
        ],
        chains: Vec::new(),
        residues: Vec::new(),
        atoms: Vec::new(),
        entities: Vec::new(),
        connections: Vec::new(),
        cispeps: Vec::new(),
        mod_residues: Vec::new(),
        helices: Vec::new(),
        sheets: Vec::new(),
        metadata: Default::default(),
        source_state: Default::default(),
        coordinates: BioCoordinateBlock::default(),
        crystal: None,
        ncs_operators: Vec::new(),
        assemblies: Vec::new(),
    };
    assert_eq!(BioStructure::validate_parts(&parts), Ok(()));
    let structure = BioStructure::from_parts(parts).unwrap();
    assert_eq!(structure.validate(), Ok(()));

    let mut parts = one_atom_parts();
    parts.chains = vec![
        chain(0, span(0, 0), None, None, None),
        chain(0, span(0, 0), None, None, None),
        chain(
            0,
            span(0, 1),
            Some(BioEntityId::new(0)),
            Some(b"A"),
            Some("LONG_LABEL"),
        ),
    ];
    parts.models[0] = BioModelRow::new(span(0, 3), Some(1));
    parts.residues[0] = residue(2, span(0, 1), Some(BioEntityId::new(0)), Some("LONG_LABEL"));
    assert_eq!(BioStructure::validate_parts(&parts), Ok(()));
}

#[test]
fn bio_invariants_reject_every_span_shape_and_overflow() {
    assert_eq!(
        BioRowSpan::<BioAtomId>::new(u32::MAX, 1),
        Err(BioStructureError::RowSpanOverflow {
            start: u32::MAX,
            len: 1,
        })
    );

    let mut gap = one_atom_parts();
    gap.models[0] = BioModelRow::new(span(1, 0), Some(1));
    assert_rejected_without_mutation(
        &gap,
        BioStructureError::NonContiguousSpan {
            table: "models->chains",
            expected_start: 0,
            actual_start: 1,
        },
    );

    let mut out_of_bounds = one_atom_parts();
    out_of_bounds.models[0] = BioModelRow::new(span(0, 2), Some(1));
    assert_rejected_without_mutation(
        &out_of_bounds,
        BioStructureError::RowSpanOutOfBounds {
            start: 0,
            len: 2,
            table_len: 1,
        },
    );

    let mut incomplete = one_atom_parts();
    incomplete.models.clear();
    assert_rejected_without_mutation(
        &incomplete,
        BioStructureError::IncompleteCoverage {
            table: "models->chains",
            covered: 0,
            table_len: 1,
        },
    );

    let mut residue_gap = one_atom_parts();
    residue_gap.chains[0] = chain(
        0,
        span(1, 0),
        Some(BioEntityId::new(0)),
        Some(b"A"),
        Some("LONG_LABEL"),
    );
    assert!(matches!(
        BioStructure::validate_parts(&residue_gap),
        Err(BioStructureError::NonContiguousSpan {
            table: "chains->residues",
            ..
        })
    ));

    let mut atom_gap = one_atom_parts();
    atom_gap.residues[0] = residue(0, span(1, 0), Some(BioEntityId::new(0)), Some("LONG_LABEL"));
    assert!(matches!(
        BioStructure::validate_parts(&atom_gap),
        Err(BioStructureError::NonContiguousSpan {
            table: "residues->atoms",
            ..
        })
    ));
}

#[test]
fn bio_invariants_reject_each_parent_mismatch() {
    let mut wrong_chain = one_atom_parts();
    wrong_chain.chains[0] = chain(
        1,
        span(0, 1),
        Some(BioEntityId::new(0)),
        Some(b"A"),
        Some("LONG_LABEL"),
    );
    assert_rejected_without_mutation(
        &wrong_chain,
        BioStructureError::ParentMismatch {
            table: "chains",
            index: 0,
        },
    );

    let mut wrong_residue = one_atom_parts();
    wrong_residue.residues[0] =
        residue(1, span(0, 1), Some(BioEntityId::new(0)), Some("LONG_LABEL"));
    assert_rejected_without_mutation(
        &wrong_residue,
        BioStructureError::ParentMismatch {
            table: "residues",
            index: 0,
        },
    );

    let mut wrong_atom = one_atom_parts();
    wrong_atom.atoms[0] = atom(1);
    assert_rejected_without_mutation(
        &wrong_atom,
        BioStructureError::ParentMismatch {
            table: "atoms",
            index: 0,
        },
    );
}

#[test]
fn bio_invariants_preserve_coordinate_bits_and_require_exact_row_count() {
    let nan = f64::from_bits(0x7ff8_0000_0000_1234);
    let mut parts = one_atom_parts();
    parts.coordinates = BioCoordinateBlock::new(vec![[-0.0, nan, f64::INFINITY]]);
    assert_eq!(BioStructure::validate_parts(&parts), Ok(()));
    let structure = BioStructure::from_parts(parts).unwrap();
    let position = structure.coordinates().positions()[0];
    assert_eq!(position[0].to_bits(), (-0.0_f64).to_bits());
    assert_eq!(position[1].to_bits(), nan.to_bits());
    assert_eq!(position[2].to_bits(), f64::INFINITY.to_bits());
    assert_eq!(structure.validate(), Ok(()));

    let mut missing = one_atom_parts();
    missing.coordinates = BioCoordinateBlock::default();
    assert_rejected_without_mutation(
        &missing,
        BioStructureError::CoordinateCountMismatch {
            atom_count: 1,
            coordinate_count: 0,
        },
    );
}

#[test]
fn bio_invariants_accept_deuterium_without_normalizing_the_row() {
    let mut parts = one_atom_parts();
    parts.atoms[0] = BioAtomRow::new(
        BioResidueId::new(0),
        atom_name(),
        Element::H,
        Some(2),
        None,
        0,
        BioCalcFlag::NotSet,
        1.0,
        20.0,
        [0.0; 6],
        -1,
        0.0,
        AtomSourceIds::new(None),
    );
    let before = parts.clone();

    assert_eq!(BioStructure::validate_parts(&parts), Ok(()));
    assert_eq!(parts, before);

    let structure = BioStructure::from_parts(parts).unwrap();
    assert_eq!(structure.atoms()[0].element(), Element::H);
    assert_eq!(structure.atoms()[0].isotope_mass_number(), Some(2));
    let rebuilt = structure.into_parts();
    assert_eq!(rebuilt.atoms[0].isotope_mass_number(), Some(2));
}

#[test]
fn bio_invariants_cover_optional_entity_and_exact_subchain_membership() {
    let mut no_reference = one_atom_parts();
    no_reference.chains[0] = chain(0, span(0, 1), None, Some(b"A"), Some("NOT_MEMBER"));
    no_reference.residues[0] = residue(0, span(0, 1), None, Some("NOT_MEMBER"));
    assert_eq!(BioStructure::validate_parts(&no_reference), Ok(()));

    let mut empty_subchain = one_atom_parts();
    empty_subchain.chains[0] = chain(
        0,
        span(0, 1),
        Some(BioEntityId::new(0)),
        Some(b"A"),
        Some(""),
    );
    empty_subchain.residues[0] = residue(0, span(0, 1), Some(BioEntityId::new(0)), Some(""));
    assert_eq!(BioStructure::validate_parts(&empty_subchain), Ok(()));

    let mut missing_entity = one_atom_parts();
    missing_entity.chains[0] = chain(
        0,
        span(0, 1),
        Some(BioEntityId::new(3)),
        Some(b"A"),
        Some("LONG_LABEL"),
    );
    assert_rejected_without_mutation(
        &missing_entity,
        BioStructureError::RowReferenceOutOfBounds {
            table: "entities",
            index: 3,
            table_len: 1,
        },
    );

    let mut wrong_membership = one_atom_parts();
    wrong_membership.residues[0] = residue(
        0,
        span(0, 1),
        Some(BioEntityId::new(0)),
        Some("CASE_SENSITIVE"),
    );
    assert_rejected_without_mutation(
        &wrong_membership,
        BioStructureError::EntitySubchainMismatch {
            entity_id: BioEntityId::new(0),
            subchain: "CASE_SENSITIVE".to_owned(),
        },
    );
}

#[test]
fn bio_invariants_validate_assembly_auth_and_label_references_exactly() {
    let mut valid = one_atom_parts();
    valid
        .assemblies
        .push(assembly(&["A", "A"], &["LONG_LABEL", "LONG_LABEL"]));
    assert_eq!(BioStructure::validate_parts(&valid), Ok(()));
    assert_eq!(valid.assemblies[0].generators[0].chains, ["A", "A"]);

    let mut missing_chain = one_atom_parts();
    missing_chain.assemblies.push(assembly(&["a"], &[]));
    assert_rejected_without_mutation(
        &missing_chain,
        BioStructureError::AssemblyReferenceMissing {
            assembly: cosmolkit_bio::BioAssemblyId::new(0),
            kind: "chain",
            value: "a".to_owned(),
        },
    );

    let mut missing_subchain = one_atom_parts();
    missing_subchain
        .assemblies
        .push(assembly(&[], &["long_label"]));
    assert_rejected_without_mutation(
        &missing_subchain,
        BioStructureError::AssemblyReferenceMissing {
            assembly: cosmolkit_bio::BioAssemblyId::new(0),
            kind: "subchain",
            value: "long_label".to_owned(),
        },
    );
}

#[test]
fn bio_invariants_source_identifier_boundaries_are_lossless_and_constructor_owned() {
    assert_eq!(PdbChainId::from_ascii(b"").unwrap().as_bytes(), b"");
    assert_eq!(PdbChainId::from_ascii(b"ABCD").unwrap().as_bytes(), b"ABCD");
    assert!(PdbChainId::from_ascii(b"ABCDE").is_none());
    assert!(PdbChainId::from_ascii(&[0xff]).is_none());

    assert_eq!(ResidueName::from_ascii(b"").unwrap().as_bytes(), b"");
    assert_eq!(
        ResidueName::from_ascii(b"ABCD").unwrap().as_bytes(),
        b"ABCD"
    );
    assert!(ResidueName::from_ascii(b"ABCDE").is_none());
    assert!(ResidueName::from_ascii(&[0xff]).is_none());
    assert!(AtomName::from_ascii(*b" CA ").is_some());
    assert!(AtomName::from_ascii([0xff, b'A', b' ', b' ']).is_none());
    assert!(ResidueSourceIds::new(None, None, Some(*b"SEG "), None, None).is_some());
    assert!(ResidueSourceIds::new(None, None, Some([0xff; 4]), None, None).is_none());

    let long = "非ASCII-mmCIF-label-that-is-longer-than-four";
    let chain_ids = ChainSourceIds::new(None, Some(long.to_owned()));
    assert_eq!(chain_ids.label_asym_id(), Some(long));
    let residue_ids = ResidueSourceIds::new(None, None, None, Some(long.to_owned()), None).unwrap();
    assert_eq!(residue_ids.subchain_id(), Some(long));
}
