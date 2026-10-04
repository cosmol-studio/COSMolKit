use cosmolkit::{
    AltLocLabel, AtomName, AtomSourceIds, BINDING_CONTRACT, BindingDefault, BindingItem,
    BindingKind, BindingOwner, BindingTypeRole, BioAtomId, BioAtomRow, BioCalcFlag, BioChainId,
    BioChainRow, BioCoordinateBlock, BioCoordinateFormat, BioModelId, BioModelRow, BioResidueId,
    BioResidueRow, BioRowSpan, BioSiftsUnpResidue, BioStructure, BioStructureParts, ChainKind,
    ChainSourceIds, Element, EntityKind, FunctionStatus, PdbAtomSerial, PdbChainId, PdbSeqId,
    Protein, ProteinAtomRef, ProteinChainRef, ProteinProjectionError, ProteinResidueRef,
    ProteinSelectionSummary, ResidueCode, ResidueInfoKind, ResidueName, ResidueSourceIds,
    StateModel,
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
        AtomName::from_ascii(b" CA ").unwrap(),
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
fn protein_input_format_preserves_all_source_formats() {
    let _: fn(&Protein) -> BioCoordinateFormat = Protein::input_format;
    for format in [
        BioCoordinateFormat::Unknown,
        BioCoordinateFormat::Detect,
        BioCoordinateFormat::Pdb,
        BioCoordinateFormat::Mmcif,
        BioCoordinateFormat::Mmjson,
        BioCoordinateFormat::ChemComp,
    ] {
        let mut parts = mixed_structure().into_parts();
        parts.input_format = format;
        let source = BioStructure::from_parts(parts).unwrap();
        let protein = source.protein().unwrap();
        let atoms = protein.as_bio_structure().atoms().as_ptr();
        assert_eq!(protein.input_format(), format);
        assert_eq!(
            protein.input_format(),
            protein.as_bio_structure().input_format()
        );
        assert_eq!(source.input_format(), format);
        assert_eq!(protein.clone().input_format(), format);
        assert_eq!(protein.as_bio_structure().atoms().as_ptr(), atoms);
    }
}

#[test]
fn protein_input_format_contract_is_exact() {
    let rows = BINDING_CONTRACT
        .iter()
        .filter(|row| row.semantic_id == "Protein.input_format")
        .collect::<Vec<_>>();
    assert_eq!(rows.len(), 1);
    let row = rows[0];
    assert_eq!(row.item, BindingItem::Callable);
    assert_eq!(row.owner, BindingOwner::Type);
    assert_eq!(row.feature, "cap-bio");
    assert_eq!(row.status, FunctionStatus::Experimental);
    assert_eq!(row.python_name, "input_format");
    assert_eq!(row.javascript_name, "inputFormat");
    let callable = row.callable.unwrap();
    assert_eq!(callable.kind, BindingKind::Instance);
    assert_eq!(callable.state_model, StateModel::ReadOnly);
    assert!(callable.parameters.is_empty());
    assert_eq!(callable.operation_semantic_id, None);
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
    let _: for<'a> fn(&'a Protein) -> &'a BioStructure = Protein::as_bio_structure;
    let _: fn(Protein) -> BioStructure = Protein::into_bio_structure;
    let _: fn(ProteinChainRef<'static>) -> BioChainId = ProteinChainRef::id;
    let _: fn(ProteinResidueRef<'static>) -> BioResidueId = ProteinResidueRef::id;
    let _: fn(ProteinAtomRef<'static>) -> BioAtomId = ProteinAtomRef::id;
}

#[test]
fn public_protein_structure_view_and_move_use_the_same_filtered_storage() {
    let source = mixed_structure();
    let protein = source.protein().unwrap();
    let first = protein.as_bio_structure();
    let second = protein.as_bio_structure();
    assert!(std::ptr::eq(first, second));
    assert_eq!(first.models().len(), 1);
    assert_eq!(first.chains().len(), 1);
    assert_eq!(first.residues().len(), 1);
    assert_eq!(first.atoms().len(), 1);
    assert_eq!(source.atoms().len(), 2);
    assert_eq!(first.atoms()[0].source().serial().unwrap().value(), 101);
    assert_eq!(first.residues()[0].name().as_str(), "ALA");
    assert_eq!(
        first.coordinates().positions()[0][0].to_bits(),
        (-0.0_f64).to_bits()
    );
    assert_eq!(
        first.coordinates().positions()[0][1].to_bits(),
        0x7ff8_0000_0000_0042
    );

    let row_ptr = first.atoms().as_ptr();
    let coordinate_ptr = first.coordinates().positions().as_ptr();
    let moved = protein.into_bio_structure();
    assert_eq!(moved.atoms().as_ptr(), row_ptr);
    assert_eq!(moved.coordinates().positions().as_ptr(), coordinate_ptr);
    assert_eq!(moved.atoms().len(), 1);
    let translated = moved.with_translated_coordinates([1.0, 0.0, 0.0]).unwrap();
    assert_eq!(translated.coordinates().positions()[0][0], 1.0);
    assert_eq!(
        moved.coordinates().positions()[0][0].to_bits(),
        (-0.0_f64).to_bits()
    );
}

#[test]
fn protein_selection_summary_counts_only_filtered_rows() {
    let _: fn(&Protein) -> ProteinSelectionSummary = Protein::selection_summary;
    let empty = BioStructure::from_parts(BioStructureParts {
        input_format: BioCoordinateFormat::Pdb,
        models: Vec::new(),
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
        coordinates: BioCoordinateBlock::new(Vec::new()),
        crystal: None,
        ncs_operators: Vec::new(),
        assemblies: Vec::new(),
    })
    .unwrap();
    let empty_protein = empty.protein().unwrap();
    assert_eq!(
        empty_protein.selection_summary(),
        ProteinSelectionSummary {
            chains: 0,
            residues: 0,
            atoms: 0
        }
    );

    let mixed = mixed_structure().protein().unwrap();
    assert_eq!(
        mixed.selection_summary(),
        ProteinSelectionSummary {
            chains: 1,
            residues: 1,
            atoms: 1
        }
    );
    assert_eq!(mixed.selection_summary().chains, mixed.num_chains());
    assert_eq!(mixed.selection_summary().residues, mixed.num_residues());
    assert_eq!(mixed.selection_summary().atoms, mixed.num_atoms());

    let two = BioStructure::from_parts(BioStructureParts {
        input_format: BioCoordinateFormat::Mmcif,
        models: vec![
            BioModelRow::new(span(0, 1), Some(1)),
            BioModelRow::new(span(1, 1), Some(2)),
        ],
        chains: vec![
            BioChainRow::new(
                BioModelId::new(0),
                None,
                span(0, 1),
                ChainKind::Protein,
                ChainSourceIds::new(None, None),
            ),
            BioChainRow::new(
                BioModelId::new(1),
                None,
                span(1, 1),
                ChainKind::Protein,
                ChainSourceIds::new(None, None),
            ),
        ],
        residues: vec![
            residue(0, 0, "ALA", ResidueInfoKind::Aa),
            residue(1, 1, "GLY", ResidueInfoKind::Aa),
        ],
        atoms: vec![atom(0, 1, None), atom(1, 2, None)],
        entities: Vec::new(),
        connections: Vec::new(),
        cispeps: Vec::new(),
        mod_residues: Vec::new(),
        helices: Vec::new(),
        sheets: Vec::new(),
        metadata: Default::default(),
        source_state: Default::default(),
        coordinates: BioCoordinateBlock::new(vec![[0.0; 3], [1.0; 3]]),
        crystal: None,
        ncs_operators: Vec::new(),
        assemblies: Vec::new(),
    })
    .unwrap()
    .protein()
    .unwrap();
    assert_eq!(two.num_models(), 2);
    assert_eq!(
        two.selection_summary(),
        ProteinSelectionSummary {
            chains: two.num_chains(),
            residues: two.num_residues(),
            atoms: two.num_atoms()
        }
    );
    assert_eq!(
        two.selection_summary(),
        ProteinSelectionSummary {
            chains: 2,
            residues: 2,
            atoms: 2
        }
    );

    for id in ["types.ProteinSelectionSummary", "Protein.selection_summary"] {
        let rows = BINDING_CONTRACT
            .iter()
            .filter(|row| row.semantic_id == id)
            .collect::<Vec<_>>();
        assert_eq!(rows.len(), 1, "{id}");
        assert_eq!(rows[0].status, FunctionStatus::Experimental);
        assert_eq!(rows[0].feature, "cap-bio");
    }
}

#[test]
fn protein_structure_conversion_registry_is_experimental_and_exact() {
    for (id, rust_name, python, javascript, state) in [
        (
            "Protein.as_bio_structure",
            "as_bio_structure",
            "as_bio_structure",
            "asBioStructure",
            StateModel::ReadOnly,
        ),
        (
            "Protein.into_bio_structure",
            "into_bio_structure",
            "into_bio_structure",
            "intoBioStructure",
            StateModel::ValueReturning,
        ),
    ] {
        let rows = BINDING_CONTRACT
            .iter()
            .filter(|row| row.semantic_id == id)
            .collect::<Vec<_>>();
        assert_eq!(rows.len(), 1, "{id}");
        let row = rows[0];
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Type);
        assert_eq!(row.feature, "cap-bio");
        assert_eq!(row.status, FunctionStatus::Experimental);
        assert_eq!(row.python_name, python);
        assert_eq!(row.javascript_name, javascript);
        assert!(row.rust_path.contains(rust_name));
        let callable = row.callable.unwrap();
        assert_eq!(callable.kind, BindingKind::Instance);
        assert_eq!(callable.state_model, state);
        assert!(callable.parameters.is_empty());
        assert_eq!(callable.operation_semantic_id, None);
    }
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
        first.atoms().next().unwrap().position().map(f64::to_bits),
        second.atoms().next().unwrap().position().map(f64::to_bits)
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

    let residue = first.residues().next().unwrap();
    assert_eq!(residue.id(), BioResidueId::new(0));
    assert_eq!(residue.name().as_str(), "ALA");
    assert_eq!(residue.code(), ResidueCode::ALA);
    assert_eq!(residue.one_letter_code(), 'A');
    assert_eq!(residue.fasta_code(), 'A');
    assert!(residue.is_standard());
    assert_eq!(residue.chain().id(), chain.id());
    assert_eq!(residue.atoms().len(), 1);

    let atom = first.atoms().next().unwrap();
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

    assert!(std::ptr::eq(
        chain.row(),
        first.chains().next().unwrap().row()
    ));
    assert!(std::ptr::eq(
        residue.row(),
        chain.residues().next().unwrap().row()
    ));
    assert!(std::ptr::eq(
        atom.row(),
        residue.atoms().next().unwrap().row()
    ));
}

#[test]
fn protein_binding_contract_matches_the_detached_public_surface() {
    let semantic_ids = [
        "types.Protein",
        "types.ProteinProjectionError",
        "types.ProteinChainRef",
        "types.ProteinResidueRef",
        "types.ProteinAtomRef",
        "types.ProteinChainIter",
        "types.ProteinResidueIter",
        "types.ProteinAtomIter",
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
        assert_eq!(row.feature, "cap-bio");
        assert_eq!(row.status, FunctionStatus::Experimental);
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
