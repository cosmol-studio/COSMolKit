use cosmolkit::{
    BINDING_CONTRACT, BindingExposure, BindingItem, BindingParity, BindingSupport,
    BioAltLocGroupId, BioAssembly, BioAssemblyGenerator, BioAssemblyId, BioAssemblyOperator,
    BioAssemblySpecialKind, BioAtomId, BioAtomRow, BioCalcFlag, BioChainId, BioChainRow,
    BioCoordinateBlock, BioCoordinateFormat, BioCrystalCell, BioCrystalInfo, BioEntityDbRef,
    BioEntityId, BioEntityRow, BioModelId, BioModelRow, BioNcsOperator, BioResidueId,
    BioResidueRow, BioRowSpan, BioSiftsUnpResidue, BioStructure, BioStructureError,
    BioStructureParts, BioTransform, ChainKind, EntityKind, PolymerKind, ResidueKind,
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
        assert_eq!(row.feature, "bio");
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
        assert!(row.callable.is_none());
    }
    for name in native_ids {
        let row = BINDING_CONTRACT
            .iter()
            .find(|row| row.semantic_id == format!("types.{name}"))
            .unwrap();
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.feature, "bio");
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::Supported);
        assert_eq!(row.parity, BindingParity::NotApplicable);
        assert!(row.callable.is_none());
    }

    assert!(!BINDING_CONTRACT.iter().any(|row| {
        row.callable
            .is_some_and(|callable| callable.operation_semantic_id == Some("BIO-hierarchy"))
    }));
}
