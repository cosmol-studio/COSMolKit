use cosmolkit::{
    AltLocLabel, AtomName, AtomSourceIds, BINDING_CONTRACT, BindingDefault, BindingExposure,
    BindingItem, BindingKind, BindingOwner, BindingParity, BindingSupport, BindingTypeRole,
    BioAtomId, BioAtomRow, BioCalcFlag, BioChainId, BioChainRow, BioCoordinateBlock,
    BioCoordinateFormat, BioModelId, BioModelRow, BioResidueId, BioResidueRow, BioRowSpan,
    BioSiftsUnpResidue, BioStructure, BioStructureParts, ChainKind, ChainSourceIds, Element,
    EntityKind, PdbAtomSerial, PdbChainId, PdbSeqId, Protein, ProteinAtomRef, ProteinChainRef,
    ProteinProjectionError, ProteinResidueRef, ResidueCode, ResidueInfoKind, ResidueName,
    ResidueSourceIds, StateModel,
};

fn span<I>(start: u32, len: u32) -> BioRowSpan<I> {
    BioRowSpan::new(start, len).unwrap()
}

fn residue(chain_id: u32, atom_start: u32, name: &str, kind: ResidueInfoKind) -> BioResidueRow {
    BioResidueRow::new(
        BioChainId::new(chain_id),
        span(atom_start, 1),
        ResidueName::from_ascii(name.as_bytes()).unwrap(),
        kind,
        EntityKind::Polymer,
        None,
        Some(b'A'),
        ResidueSourceIds::new(
            Some(PdbSeqId::new(atom_start as i32 + 1, None)),
            Some(atom_start as i32 + 1),
            Some(*b"SEG "),
            Some("label-A".to_owned()),
            None,
        )
        .unwrap(),
        BioSiftsUnpResidue::new(Some(b'A'), 1, u16::try_from(atom_start + 1).unwrap()),
    )
}

fn atom(residue_id: u32, serial: i32, altloc: Option<u8>) -> BioAtomRow {
    BioAtomRow::new(
        BioResidueId::new(residue_id),
        AtomName::from_ascii(*b" CA ").unwrap(),
        Element::C,
        None,
        altloc.map(AltLocLabel::new),
        0,
        BioCalcFlag::NotSet,
        1.0,
        10.0,
        [0.0; 6],
        0,
        1.0,
        AtomSourceIds::new(Some(PdbAtomSerial::new(serial))),
    )
}

fn mixed_structure() -> BioStructure {
    BioStructure::from_parts(BioStructureParts {
        input_format: BioCoordinateFormat::Mmcif,
        models: vec![BioModelRow::new(span(0, 1), Some(7))],
        chains: vec![BioChainRow::new(
            BioModelId::new(0),
            None,
            span(0, 2),
            ChainKind::Mixed,
            ChainSourceIds::new(
                Some(PdbChainId::from_ascii(b"AUTH").unwrap()),
                Some("label-A".to_owned()),
            ),
        )],
        residues: vec![
            residue(0, 0, "ALA", ResidueInfoKind::Aa),
            residue(0, 1, "HOH", ResidueInfoKind::Hoh),
        ],
        atoms: vec![atom(0, 101, Some(b'B')), atom(1, 102, None)],
        entities: Vec::new(),
        connections: Vec::new(),
        cispeps: Vec::new(),
        mod_residues: Vec::new(),
        helices: Vec::new(),
        sheets: Vec::new(),
        metadata: Default::default(),
        source_state: Default::default(),
        coordinates: BioCoordinateBlock::new(vec![
            [-0.0, f64::from_bits(0x7ff8_0000_0000_0042), 0.0005],
            [9.0, 8.0, 7.0],
        ]),
        crystal: None,
        ncs_operators: Vec::new(),
        assemblies: Vec::new(),
    })
    .unwrap()
}

fn compact(value: &str) -> String {
    value
        .chars()
        .filter(|character| !character.is_whitespace())
        .collect()
}

#[test]
fn canonical_protein_signatures_and_types_are_available_from_the_facade() {
    fn assert_type<T>() {}

    assert_type::<Protein>();
    assert_type::<ProteinProjectionError>();
    assert_type::<ProteinChainRef<'static>>();
    assert_type::<ProteinResidueRef<'static>>();
    assert_type::<ProteinAtomRef<'static>>();

    let _: fn(&BioStructure) -> Result<Protein, ProteinProjectionError> = BioStructure::protein;
    let _: fn(&Protein) -> usize = Protein::num_atoms;
    let _: fn(ProteinChainRef<'static>) -> BioChainId = ProteinChainRef::id;
    let _: fn(ProteinResidueRef<'static>) -> BioResidueId = ProteinResidueRef::id;
    let _: fn(ProteinAtomRef<'static>) -> BioAtomId = ProteinAtomRef::id;
}

#[test]
fn public_projection_preserves_order_views_coordinate_bits_and_source_state() {
    let source = mixed_structure();
    let before_positions = source
        .coordinates()
        .positions()
        .iter()
        .map(|position| position.map(f64::to_bits))
        .collect::<Vec<_>>();
    let first = source.protein().unwrap();
    let second = source.protein().unwrap();

    assert_eq!(
        (
            first.num_models(),
            first.num_chains(),
            first.num_residues(),
            first.num_atoms(),
        ),
        (1, 1, 1, 1)
    );
    assert_eq!(
        (
            second.num_models(),
            second.num_chains(),
            second.num_residues(),
            second.num_atoms(),
        ),
        (1, 1, 1, 1)
    );
    assert_eq!(
        first.atoms()[0].position().map(f64::to_bits),
        second.atoms()[0].position().map(f64::to_bits)
    );
    assert_eq!(source.models().len(), 1);
    assert_eq!(source.chains().len(), 1);
    assert_eq!(source.residues().len(), 2);
    assert_eq!(source.atoms().len(), 2);
    assert_eq!(
        source
            .coordinates()
            .positions()
            .iter()
            .map(|position| position.map(f64::to_bits))
            .collect::<Vec<_>>(),
        before_positions
    );

    let chain = first.chain(0).unwrap();
    assert!(first.chain(1).is_none());
    assert_eq!(chain.id(), BioChainId::new(0));
    assert_eq!(chain.kind(), ChainKind::Protein);
    assert_eq!(chain.source().auth_chain_id().unwrap().as_str(), "AUTH");
    assert_eq!(chain.source().label_asym_id(), Some("label-A"));
    assert_eq!(chain.residues().len(), 1);
    assert_eq!(chain.atoms().len(), 1);

    let residue = first.residues()[0];
    assert_eq!(residue.id(), BioResidueId::new(0));
    assert_eq!(residue.name().as_str(), "ALA");
    assert_eq!(residue.code(), ResidueCode::ALA);
    assert_eq!(residue.one_letter_code(), 'A');
    assert_eq!(residue.fasta_code(), 'A');
    assert!(residue.is_standard());
    assert_eq!(residue.chain().id(), chain.id());
    assert_eq!(residue.atoms().len(), 1);

    let atom = first.atoms()[0];
    assert_eq!(atom.id(), BioAtomId::new(0));
    assert_eq!(atom.name().as_str(), " CA ");
    assert_eq!(atom.element(), Element::C);
    assert_eq!(atom.altloc().unwrap().value(), b'B');
    assert_eq!(atom.residue().id(), residue.id());
    assert_eq!(atom.row().source().serial().unwrap().value(), 101);
    let position = atom.position();
    assert_eq!(position[0].to_bits(), (-0.0_f64).to_bits());
    assert_eq!(position[1].to_bits(), 0x7ff8_0000_0000_0042);
    assert_eq!(position[2].to_bits(), 0.0005_f64.to_bits());

    assert!(std::ptr::eq(chain.row(), first.chains()[0].row()));
    assert!(std::ptr::eq(residue.row(), chain.residues()[0].row()));
    assert!(std::ptr::eq(atom.row(), residue.atoms()[0].row()));
}

#[test]
fn protein_binding_contract_matches_the_detached_public_surface() {
    let semantic_ids = [
        "types.Protein",
        "types.ProteinProjectionError",
        "types.ProteinChainRef",
        "types.ProteinResidueRef",
        "types.ProteinAtomRef",
        "BioStructure.protein",
        "Protein.num_models",
        "Protein.num_chains",
        "Protein.num_residues",
        "Protein.num_atoms",
        "Protein.chains",
        "Protein.chain",
        "Protein.residues",
        "Protein.atoms",
        "ProteinChainRef.id",
        "ProteinChainRef.row",
        "ProteinChainRef.kind",
        "ProteinChainRef.source",
        "ProteinChainRef.residues",
        "ProteinChainRef.atoms",
        "ProteinResidueRef.id",
        "ProteinResidueRef.row",
        "ProteinResidueRef.name",
        "ProteinResidueRef.kind",
        "ProteinResidueRef.info",
        "ProteinResidueRef.code",
        "ProteinResidueRef.one_letter_code",
        "ProteinResidueRef.fasta_code",
        "ProteinResidueRef.is_standard",
        "ProteinResidueRef.chain",
        "ProteinResidueRef.atoms",
        "ProteinAtomRef.id",
        "ProteinAtomRef.row",
        "ProteinAtomRef.name",
        "ProteinAtomRef.element",
        "ProteinAtomRef.altloc",
        "ProteinAtomRef.residue",
        "ProteinAtomRef.position",
    ];
    let native_accessors = [
        "ProteinChainRef.id",
        "ProteinChainRef.row",
        "ProteinResidueRef.id",
        "ProteinResidueRef.row",
        "ProteinAtomRef.id",
        "ProteinAtomRef.row",
    ];

    let rows = BINDING_CONTRACT
        .iter()
        .filter(|row| semantic_ids.contains(&row.semantic_id))
        .collect::<Vec<_>>();
    assert_eq!(rows.len(), semantic_ids.len());
    for semantic_id in semantic_ids {
        assert_eq!(
            rows.iter()
                .filter(|row| row.semantic_id == semantic_id)
                .count(),
            1,
            "duplicate or missing contract row {semantic_id}"
        );
    }

    for row in rows {
        assert_eq!(row.owner, BindingOwner::Type);
        assert_eq!(row.feature, "bio");
        assert_eq!(row.exposure, BindingExposure::Public);
        if native_accessors.contains(&row.semantic_id) {
            assert_eq!(row.support, BindingSupport::Supported);
            assert_eq!(row.parity, BindingParity::NotApplicable);
        } else {
            assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
            assert_eq!(row.parity, BindingParity::RequiredNow);
        }

        if row.semantic_id.starts_with("types.") {
            assert_eq!(row.item, BindingItem::Type);
            assert!(row.callable.is_none());
            let expected_role = if row.semantic_id == "types.ProteinProjectionError" {
                BindingTypeRole::Error
            } else {
                BindingTypeRole::Value
            };
            assert_eq!(row.type_role, Some(expected_role));
        } else {
            assert_eq!(row.item, BindingItem::Callable);
            assert!(row.type_role.is_none());
            let callable = row.callable.unwrap();
            assert_eq!(callable.kind, BindingKind::Instance);
            assert_eq!(callable.operation_semantic_id, None);
            let expected_state = if row.semantic_id.starts_with("ProteinChainRef.")
                || row.semantic_id.starts_with("ProteinResidueRef.")
                || row.semantic_id.starts_with("ProteinAtomRef.")
            {
                StateModel::ValueReturning
            } else {
                StateModel::ReadOnly
            };
            assert_eq!(callable.state_model, expected_state);
            if row.semantic_id == "Protein.chain" {
                assert_eq!(callable.parameters.len(), 1);
                assert_eq!(callable.parameters[0].name, "index");
                assert_eq!(compact(callable.parameters[0].type_name), "usize");
                assert_eq!(callable.parameters[0].default, BindingDefault::Required);
            } else {
                assert!(callable.parameters.is_empty());
            }
        }
    }
}
