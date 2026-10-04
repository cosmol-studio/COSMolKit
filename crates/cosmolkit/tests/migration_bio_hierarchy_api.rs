use cosmolkit::{
    AtomName, AtomSourceIds, BINDING_CONTRACT, BindingItem, BioAltLocGroupId, BioAssembly,
    BioAssemblyGenerator, BioAssemblyId, BioAssemblyOperator, BioAssemblySpecialKind, BioAtomId,
    BioAtomRow, BioCalcFlag, BioChainId, BioChainRow, BioCoordinateBlock, BioCoordinateFormat,
    BioCrystalCell, BioCrystalInfo, BioEntityDbRef, BioEntityId, BioEntityRow, BioModelId,
    BioModelRow, BioNcsOperator, BioResidueId, BioResidueRow, BioRowSpan, BioSiftsUnpResidue,
    BioStructure, BioStructureError, BioStructureParts, BioTransform, ChainKind, ChainSourceIds,
    Element, EntityKind, EntitySourceIds, FunctionStatus, PolymerKind, ResidueInfoKind,
    ResidueKind, ResidueName, ResidueSourceIds,
};

fn empty_parts() -> BioStructureParts {
    BioStructureParts {
        input_format: BioCoordinateFormat::Unknown,
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
        coordinates: BioCoordinateBlock::default(),
        crystal: None,
        ncs_operators: Vec::new(),
        assemblies: Vec::new(),
    }
}

#[test]
fn bio_registry_covers_existing_accessors_and_owned_conversion() {
    for name in [
        "from_parts",
        "validate_parts",
        "validate",
        "into_parts",
        "input_format",
        "models",
        "chains",
        "residues",
        "atoms",
        "entities",
        "connections",
        "cispeps",
        "mod_residues",
        "helices",
        "sheets",
        "metadata",
        "source_state",
        "coordinates",
        "crystal",
        "ncs_operators",
        "assemblies",
        "find_entity",
        "find_entity_of_subchain",
        "find_atom",
        "atom_by_altloc",
    ] {
        let id = format!("BioStructure.{name}");
        let rows: Vec<_> = BINDING_CONTRACT
            .iter()
            .filter(|row| row.semantic_id == id)
            .collect();
        assert_eq!(rows.len(), 1, "{id}");
        assert_eq!(rows[0].status, FunctionStatus::Experimental, "{id}");
        let expected = match name {
            "from_parts" | "validate_parts" => None,
            "into_parts" => Some(cosmolkit::BindingReceiver::Owned),
            _ => Some(cosmolkit::BindingReceiver::Shared),
        };
        assert_eq!(rows[0].callable.unwrap().receiver, expected, "{id}");
    }
    for row in BINDING_CONTRACT
        .iter()
        .filter(|row| row.feature == "cap-bio" && row.item == BindingItem::Callable)
    {
        assert_eq!(
            row.status,
            FunctionStatus::Experimental,
            "{}",
            row.semantic_id
        );
    }
    let structure = BioStructure::from_parts(empty_parts()).unwrap();
    let parts = structure.into_parts();
    BioStructure::validate_parts(&parts).unwrap();
}

#[test]
fn bio_public_source_state_metadata_is_borrowed_and_bit_exact() {
    let defaults = BioStructure::from_parts(empty_parts()).unwrap();
    assert_eq!(defaults.name(), "");
    assert!(!defaults.has_origx());
    assert_eq!(defaults.resolution().to_bits(), 0.0_f64.to_bits());
    assert_eq!(defaults.ter_status(), 0);
    assert!(std::ptr::eq(
        defaults.origx(),
        &defaults.source_state().origx
    ));

    let mut parts = empty_parts();
    parts.source_state.name = "filtered source".into();
    parts.source_state.resolution = f64::from_bits(0x8000_0000_0000_0000);
    parts.source_state.ter_status = b'e';
    let stored_unset = BioStructure::from_parts(parts.clone()).unwrap();
    assert!(!stored_unset.has_origx());
    assert!(std::ptr::eq(
        stored_unset.origx(),
        &stored_unset.source_state().origx
    ));
    parts.source_state.has_origx = true;
    let explicit = BioStructure::from_parts(parts).unwrap();
    assert_eq!(explicit.name(), "filtered source");
    assert!(explicit.has_origx());
    assert_eq!(explicit.resolution().to_bits(), (-0.0_f64).to_bits());
    assert_eq!(explicit.ter_status(), b'e');
    assert!(std::ptr::eq(
        explicit.origx(),
        &explicit.source_state().origx
    ));
    let protein = explicit.protein().unwrap();
    let view = protein.as_bio_structure();
    assert_eq!(view.name(), explicit.name());
    assert_eq!(view.has_origx(), explicit.has_origx());
    assert_eq!(view.resolution().to_bits(), explicit.resolution().to_bits());
    assert_eq!(view.ter_status(), explicit.ter_status());
    assert_eq!(view.origx(), explicit.origx());

    let _: for<'a> fn(&'a BioStructure) -> &'a str = BioStructure::name;
    let _: fn(&BioStructure) -> bool = BioStructure::has_origx;
    let _: for<'a> fn(&'a BioStructure) -> &'a BioTransform = BioStructure::origx;
    let _: fn(&BioStructure) -> f64 = BioStructure::resolution;
    let _: fn(&BioStructure) -> u8 = BioStructure::ter_status;
    for id in ["name", "has_origx", "origx", "resolution", "ter_status"] {
        let rows = BINDING_CONTRACT
            .iter()
            .filter(|row| row.semantic_id == format!("BioStructure.{id}"))
            .collect::<Vec<_>>();
        assert_eq!(rows.len(), 1, "{id}");
        assert_eq!(rows[0].status, FunctionStatus::Experimental);
        assert_eq!(rows[0].feature, "cap-bio");
    }
}

#[test]
fn bio_public_counts_use_authoritative_blocks_and_filtered_protein_view() {
    let empty = BioStructure::from_parts(empty_parts()).unwrap();
    assert_eq!(
        (
            empty.num_models(),
            empty.num_chains(),
            empty.num_residues(),
            empty.num_atoms(),
            empty.num_entities()
        ),
        (0, 0, 0, 0, 0)
    );

    let mut parts = empty_parts();
    parts.models = vec![
        BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1)),
        BioModelRow::new(BioRowSpan::new(1, 1).unwrap(), Some(2)),
    ];
    parts.chains = vec![
        BioChainRow::new(
            BioModelId::new(0),
            None,
            BioRowSpan::new(0, 1).unwrap(),
            ChainKind::Protein,
            ChainSourceIds::default(),
        ),
        BioChainRow::new(
            BioModelId::new(1),
            None,
            BioRowSpan::new(1, 1).unwrap(),
            ChainKind::WaterOnly,
            ChainSourceIds::default(),
        ),
    ];
    parts.residues = vec![
        BioResidueRow::new(
            BioChainId::new(0),
            BioRowSpan::new(0, 1).unwrap(),
            ResidueName::from_ascii(b"ALA").unwrap(),
            ResidueInfoKind::Aa,
            EntityKind::Polymer,
            None,
            None,
            ResidueSourceIds::default(),
            BioSiftsUnpResidue::default(),
        ),
        BioResidueRow::new(
            BioChainId::new(1),
            BioRowSpan::new(1, 1).unwrap(),
            ResidueName::from_ascii(b"HOH").unwrap(),
            ResidueInfoKind::Hoh,
            EntityKind::Water,
            None,
            None,
            ResidueSourceIds::default(),
            BioSiftsUnpResidue::default(),
        ),
    ];
    parts.atoms = (0..2)
        .map(|index| {
            BioAtomRow::new(
                BioResidueId::new(index),
                AtomName::from_ascii(b" O  ").unwrap(),
                Element::O,
                None,
                None,
                0,
                BioCalcFlag::NotSet,
                1.0,
                0.0,
                [0.0; 6],
                0,
                1.0,
                AtomSourceIds::default(),
            )
        })
        .collect();
    parts.coordinates = BioCoordinateBlock::new(vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
    parts.entities = vec![
        BioEntityRow::new(
            EntityKind::Polymer,
            PolymerKind::PeptideL,
            false,
            vec![],
            vec![],
            vec![],
            vec![],
            EntitySourceIds::new("1".into()),
        ),
        BioEntityRow::new(
            EntityKind::Water,
            PolymerKind::Unknown,
            false,
            vec![],
            vec![],
            vec![],
            vec![],
            EntitySourceIds::new("2".into()),
        ),
    ];
    let source = BioStructure::from_parts(parts).unwrap();
    assert_eq!(
        (
            source.num_models(),
            source.num_chains(),
            source.num_residues(),
            source.num_atoms(),
            source.num_entities()
        ),
        (2, 2, 2, 2, 2)
    );
    let protein = source.protein().unwrap();
    let view = protein.as_bio_structure();
    assert_eq!(view.num_models(), view.models().len());
    assert_eq!(view.num_chains(), protein.num_chains());
    assert_eq!(view.num_residues(), protein.num_residues());
    assert_eq!(view.num_atoms(), protein.num_atoms());
    assert_eq!(view.num_entities(), view.entities().len());
    assert_eq!(view.num_atoms(), 1);
    assert_eq!(source.num_atoms(), 2);
    assert_eq!(
        source.atom_position(BioAtomId::new(0)),
        Some([1.0, 2.0, 3.0])
    );
    assert_eq!(
        source.atom_position(BioAtomId::new(1)),
        Some([4.0, 5.0, 6.0])
    );
    assert_eq!(
        source.residue_atoms(BioResidueId::new(1)),
        Some(&source.atoms()[1..2])
    );
    assert_eq!(view.atom_position(BioAtomId::new(0)), Some([1.0, 2.0, 3.0]));
    for (name, method) in [
        (
            "num_models",
            BioStructure::num_models as fn(&BioStructure) -> usize,
        ),
        ("num_chains", BioStructure::num_chains),
        ("num_residues", BioStructure::num_residues),
        ("num_atoms", BioStructure::num_atoms),
        ("num_entities", BioStructure::num_entities),
    ] {
        let row = BINDING_CONTRACT
            .iter()
            .find(|row| row.semantic_id == format!("BioStructure.{name}"))
            .unwrap();
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.feature, "cap-bio");
        assert_eq!(row.status, FunctionStatus::Experimental);
        assert_eq!(
            method(&source),
            match name {
                "num_models" => 2,
                "num_chains" => 2,
                "num_residues" => 2,
                "num_atoms" => 2,
                _ => 2,
            }
        );
    }
}

#[test]
fn bio_public_navigation_preserves_bits_borrows_and_empty_residues() {
    let mut parts = empty_parts();
    parts.models = vec![BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(7))];
    parts.chains = vec![BioChainRow::new(
        BioModelId::new(0),
        None,
        BioRowSpan::new(0, 2).unwrap(),
        ChainKind::Protein,
        ChainSourceIds::default(),
    )];
    parts.residues = vec![
        BioResidueRow::new(
            BioChainId::new(0),
            BioRowSpan::new(0, 1).unwrap(),
            ResidueName::from_ascii(b"ALA").unwrap(),
            ResidueInfoKind::Aa,
            EntityKind::Polymer,
            None,
            None,
            ResidueSourceIds::default(),
            BioSiftsUnpResidue::default(),
        ),
        BioResidueRow::new(
            BioChainId::new(0),
            BioRowSpan::new(1, 0).unwrap(),
            ResidueName::from_ascii(b"GLY").unwrap(),
            ResidueInfoKind::Aa,
            EntityKind::Polymer,
            None,
            None,
            ResidueSourceIds::default(),
            BioSiftsUnpResidue::default(),
        ),
    ];
    parts.atoms = vec![BioAtomRow::new(
        BioResidueId::new(0),
        AtomName::from_ascii(b" CA ").unwrap(),
        Element::C,
        None,
        None,
        0,
        BioCalcFlag::NotSet,
        1.0,
        0.0,
        [0.0; 6],
        0,
        1.0,
        AtomSourceIds::default(),
    )];
    parts.coordinates = BioCoordinateBlock::new(vec![[-0.0, f64::from_bits(1), 3.0]]);
    let structure = BioStructure::from_parts(parts).unwrap();
    assert_eq!(
        structure
            .atom_position(BioAtomId::new(0))
            .unwrap()
            .map(f64::to_bits),
        [(-0.0_f64).to_bits(), 1, 3.0_f64.to_bits()]
    );
    assert_eq!(structure.atom_position(BioAtomId::new(1)), None);
    assert_eq!(structure.atom_position(BioAtomId::new(u32::MAX)), None);
    let first = structure.residue_atoms(BioResidueId::new(0)).unwrap();
    assert_eq!(first.len(), 1);
    assert!(std::ptr::eq(first.as_ptr(), structure.atoms().as_ptr()));
    assert_eq!(structure.residue_atoms(BioResidueId::new(1)), Some(&[][..]));
    assert_eq!(structure.residue_atoms(BioResidueId::new(2)), None);
    let _: fn(&BioStructure, BioAtomId) -> Option<[f64; 3]> = BioStructure::atom_position;
    let _: for<'a> fn(&'a BioStructure, BioResidueId) -> Option<&'a [BioAtomRow]> =
        BioStructure::residue_atoms;
    for name in ["atom_position", "residue_atoms"] {
        let row = BINDING_CONTRACT
            .iter()
            .find(|row| row.semantic_id == format!("BioStructure.{name}"))
            .unwrap();
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.feature, "cap-bio");
        assert_eq!(row.status, FunctionStatus::Experimental);
    }
}

#[test]
fn canonical_hierarchy_types_are_available_from_the_public_facade() {
    fn assert_type<T>() {}

    assert_type::<BioStructure>();
    assert_type::<BioStructureParts>();
    assert_type::<BioStructureError>();
    assert_type::<BioCoordinateFormat>();
    assert_type::<BioCalcFlag>();
    assert_type::<EntityKind>();
    assert_type::<PolymerKind>();
    assert_type::<ResidueKind>();
    assert_type::<ChainKind>();
    assert_type::<BioRowSpan<BioAtomId>>();
    assert_type::<BioAtomId>();
    assert_type::<BioResidueId>();
    assert_type::<BioChainId>();
    assert_type::<BioEntityId>();
    assert_type::<BioModelId>();
    assert_type::<BioAssemblyId>();
    assert_type::<BioAltLocGroupId>();
    assert_type::<BioAtomRow>();
    assert_type::<BioResidueRow>();
    assert_type::<BioChainRow>();
    assert_type::<BioEntityRow>();
    assert_type::<BioModelRow>();
    assert_type::<BioCoordinateBlock>();
    assert_type::<BioTransform>();
    assert_type::<BioCrystalCell>();
    assert_type::<BioCrystalInfo>();
    assert_type::<BioNcsOperator>();
    assert_type::<BioAssemblyOperator>();
    assert_type::<BioAssemblyGenerator>();
    assert_type::<BioAssemblySpecialKind>();
    assert_type::<BioAssembly>();
    assert_type::<cosmolkit::AltLocRequest>();
    assert_type::<BioEntityDbRef>();
    assert_type::<BioSiftsUnpResidue>();
}

#[test]
fn public_hierarchy_construction_is_validated_and_failure_is_structured() {
    let empty = BioStructure::from_parts(empty_parts()).unwrap();
    assert!(empty.models().is_empty());
    assert!(empty.chains().is_empty());
    assert!(empty.residues().is_empty());
    assert!(empty.atoms().is_empty());
    assert!(empty.entities().is_empty());
    assert!(empty.coordinates().positions().is_empty());
    assert!(empty.crystal().is_none());
    assert!(empty.ncs_operators().is_empty());
    assert!(empty.assemblies().is_empty());
    assert_eq!(empty.input_format(), BioCoordinateFormat::Unknown);
    assert_eq!(empty.into_parts(), empty_parts());

    let mut invalid = empty_parts();
    invalid.coordinates = BioCoordinateBlock::new(vec![[-0.0, f64::from_bits(1), f64::NAN]]);
    let unchanged = invalid.clone();
    assert_eq!(
        BioStructure::from_parts(invalid.clone()),
        Err(BioStructureError::CoordinateCountMismatch {
            atom_count: 0,
            coordinate_count: 1,
        })
    );
    assert_eq!(invalid.input_format, unchanged.input_format);
    assert_eq!(invalid.models, unchanged.models);
    assert_eq!(invalid.chains, unchanged.chains);
    assert_eq!(invalid.residues, unchanged.residues);
    assert_eq!(invalid.atoms, unchanged.atoms);
    assert_eq!(invalid.entities, unchanged.entities);
    assert_eq!(invalid.crystal, unchanged.crystal);
    assert_eq!(invalid.ncs_operators, unchanged.ncs_operators);
    assert_eq!(invalid.assemblies, unchanged.assemblies);
    assert_eq!(
        invalid.coordinates.positions().len(),
        unchanged.coordinates.positions().len()
    );
    for (actual, expected) in invalid
        .coordinates
        .positions()
        .iter()
        .zip(unchanged.coordinates.positions())
    {
        assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
    }
    assert_eq!(
        BioRowSpan::<BioAtomId>::new(u32::MAX, 1),
        Err(BioStructureError::RowSpanOverflow {
            start: u32::MAX,
            len: 1,
        })
    );
}

#[test]
fn hierarchy_binding_contract_matches_the_public_type_projection() {
    let source_backed = [
        "BioStructure",
        "BioStructureParts",
        "BioStructureError",
        "BioCoordinateFormat",
        "BioCalcFlag",
        "EntityKind",
        "PolymerKind",
        "ResidueKind",
        "ChainKind",
        "BioAtomRow",
        "BioResidueRow",
        "BioChainRow",
        "BioEntityRow",
        "BioModelRow",
        "BioCoordinateBlock",
        "BioTransform",
        "BioCrystalCell",
        "BioCrystalInfo",
        "BioNcsOperator",
        "BioAssemblyOperator",
        "BioAssemblyGenerator",
        "BioAssemblySpecialKind",
        "BioAssembly",
        "AltLocRequest",
        "BioEntityDbRef",
        "BioSiftsUnpResidue",
    ];
    let native_ids = [
        "BioRowSpan",
        "BioAtomId",
        "BioResidueId",
        "BioChainId",
        "BioEntityId",
        "BioModelId",
        "BioAssemblyId",
        "BioAltLocGroupId",
    ];

    for name in source_backed {
        let row = BINDING_CONTRACT
            .iter()
            .find(|row| row.semantic_id == format!("types.{name}"))
            .unwrap();
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.feature, "cap-bio");
        assert_eq!(row.status, FunctionStatus::Experimental);
        assert!(row.callable.is_none());
    }
    for name in native_ids {
        let row = BINDING_CONTRACT
            .iter()
            .find(|row| row.semantic_id == format!("types.{name}"))
            .unwrap();
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.feature, "cap-bio");
        assert_eq!(row.status, FunctionStatus::Experimental);
        assert!(row.callable.is_none());
    }

    assert!(!BINDING_CONTRACT.iter().any(|row| {
        row.callable
            .is_some_and(|callable| callable.operation_semantic_id == Some("BIO-hierarchy"))
    }));
}
