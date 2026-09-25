use cosmolkit_bio::{
    AltLocLabel, AltLocRequest, AtomName, AtomSourceIds, BioAssembly, BioAssemblyGenerator,
    BioAssemblyOperator, BioAssemblySpecialKind, BioAtomId, BioAtomRow, BioCalcFlag, BioChainId,
    BioChainRow, BioCoordinateBlock, BioCoordinateFormat, BioCrystalCell, BioCrystalInfo,
    BioEntityDbRef, BioEntityId, BioEntityRow, BioModelId, BioModelRow, BioNcsOperator,
    BioResidueId, BioResidueRow, BioRowSpan, BioSiftsUnpResidue, BioStructure, BioStructureError,
    BioStructureParts, BioTransform, ChainKind, ChainSourceIds, EntityKind, EntitySourceIds,
    PdbAtomSerial, PdbChainId, PdbSeqId, PolymerKind, ResidueInfoKind, ResidueName,
    ResidueSourceIds, altloc_matches, is_same_conformer,
};
use cosmolkit_types::Element;

fn atom_name(value: &[u8; 4]) -> AtomName {
    AtomName::from_ascii(*value).unwrap()
}

fn residue_name(value: &str) -> ResidueName {
    ResidueName::from_ascii(value.as_bytes()).unwrap()
}

fn span<I>(start: u32, len: u32) -> BioRowSpan<I> {
    BioRowSpan::new(start, len).unwrap()
}

fn entity(source_id: &str, subchains: &[&str]) -> BioEntityRow {
    BioEntityRow::new(
        EntityKind::Polymer,
        PolymerKind::PeptideL,
        true,
        vec!["ALA,GLY".to_owned(), "SER".to_owned()],
        vec![BioEntityDbRef {
            db_name: "UNP".to_owned(),
            accession_code: "P00001".to_owned(),
            id_code: "ENTRY".to_owned(),
            isoform: String::new(),
            seq_begin: PdbSeqId::new(1, None),
            seq_end: PdbSeqId::new(2, None),
            db_begin: PdbSeqId::new(10, None),
            db_end: PdbSeqId::new(11, None),
            label_seq_begin: Some(1),
            label_seq_end: Some(2),
        }],
        vec!["P00001".to_owned()],
        subchains.iter().map(|value| (*value).to_owned()).collect(),
        EntitySourceIds::new(source_id.to_owned()),
    )
}

fn one_atom_parts(position: [f64; 3]) -> BioStructureParts {
    let auth = PdbChainId::from_ascii(b"A").unwrap();
    BioStructureParts {
        input_format: BioCoordinateFormat::Mmcif,
        models: vec![BioModelRow::new(span(0, 1), Some(7))],
        chains: vec![BioChainRow::new(
            BioModelId::new(0),
            Some(BioEntityId::new(0)),
            span(0, 1),
            ChainKind::Protein,
            ChainSourceIds::new(Some(auth), Some("LABEL_LONG".to_owned())),
        )],
        residues: vec![BioResidueRow::new(
            BioChainId::new(0),
            span(0, 1),
            residue_name("ALA"),
            ResidueInfoKind::Aa,
            EntityKind::Polymer,
            Some(BioEntityId::new(0)),
            Some(b'A'),
            ResidueSourceIds::new(
                Some(PdbSeqId::new(4, Some(b'B'))),
                Some(1),
                Some(*b"SEG "),
                Some("LABEL_LONG".to_owned()),
                Some("1".to_owned()),
            )
            .unwrap(),
            BioSiftsUnpResidue::new(Some(b'A'), 0, 10),
        )],
        atoms: vec![BioAtomRow::new(
            BioResidueId::new(0),
            atom_name(b" CA "),
            Element::C,
            None,
            0,
            BioCalcFlag::NotSet,
            1.0,
            20.0,
            [0.0; 6],
            -1,
            0.0,
            AtomSourceIds::new(Some(PdbAtomSerial::new(99))),
        )],
        entities: vec![entity("1", &["LABEL_LONG"])],
        coordinates: BioCoordinateBlock::new(vec![position]),
        crystal: None,
        ncs_operators: Vec::new(),
        assemblies: Vec::new(),
    }
}

#[test]
fn bio_hierarchy_validates_empty_single_and_multi_model_contiguous_spans() {
    let empty = BioStructure::from_parts(BioStructureParts {
        input_format: BioCoordinateFormat::Unknown,
        models: Vec::new(),
        chains: Vec::new(),
        residues: Vec::new(),
        atoms: Vec::new(),
        entities: Vec::new(),
        coordinates: BioCoordinateBlock::default(),
        crystal: None,
        ncs_operators: Vec::new(),
        assemblies: Vec::new(),
    })
    .unwrap();
    assert!(empty.models().is_empty());

    let single = BioStructure::from_parts(one_atom_parts([1.0, 2.0, 3.0])).unwrap();
    assert_eq!(single.models().len(), 1);
    assert_eq!(
        single.chains()[0].source().label_asym_id(),
        Some("LABEL_LONG")
    );
    assert_eq!(
        single.residues()[0].source().subchain_id(),
        Some("LABEL_LONG")
    );

    let mut invalid = one_atom_parts([0.0; 3]);
    invalid.models[0] = BioModelRow::new(span(1, 0), Some(7));
    assert!(matches!(
        BioStructure::from_parts(invalid),
        Err(BioStructureError::NonContiguousSpan { .. })
    ));

    let mut wrong_parent = one_atom_parts([0.0; 3]);
    wrong_parent.atoms[0] = BioAtomRow::new(
        BioResidueId::new(1),
        atom_name(b" CA "),
        Element::C,
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
    assert!(matches!(
        BioStructure::from_parts(wrong_parent),
        Err(BioStructureError::ParentMismatch { table: "atoms", .. })
    ));
    assert!(matches!(
        BioRowSpan::<BioAtomId>::new(u32::MAX, 1),
        Err(BioStructureError::RowSpanOverflow { .. })
    ));
}

#[test]
fn bio_hierarchy_preserves_coordinate_bits_and_rejects_misalignment() {
    let nan = f64::from_bits(0x7ff8_0000_0000_0042);
    let position = [-0.0, nan, f64::INFINITY];
    let structure = BioStructure::from_parts(one_atom_parts(position)).unwrap();
    let stored = structure.coordinates().positions()[0];
    assert_eq!(stored[0].to_bits(), (-0.0_f64).to_bits());
    assert_eq!(stored[1].to_bits(), nan.to_bits());
    assert_eq!(stored[2].to_bits(), f64::INFINITY.to_bits());
    let parts = structure.into_parts();
    assert_eq!(parts.coordinates.positions()[0][1].to_bits(), nan.to_bits());
    assert!(BioStructure::from_parts(parts).is_ok());

    let mut invalid = one_atom_parts([0.0; 3]);
    invalid.coordinates = BioCoordinateBlock::default();
    assert!(matches!(
        BioStructure::from_parts(invalid),
        Err(BioStructureError::CoordinateCountMismatch {
            atom_count: 1,
            coordinate_count: 0
        })
    ));
}

#[test]
fn bio_hierarchy_altloc_matrix_and_source_order_lookup_match_gemmi() {
    let a = Some(AltLocLabel::new(b'A'));
    let b = Some(AltLocLabel::new(b'B'));
    assert!(is_same_conformer(None, b));
    assert!(is_same_conformer(a, None));
    assert!(is_same_conformer(a, a));
    assert!(!is_same_conformer(a, b));
    assert!(altloc_matches(b, AltLocRequest::Any));
    assert!(altloc_matches(None, AltLocRequest::Exact(a)));
    assert!(altloc_matches(a, AltLocRequest::Exact(a)));
    assert!(!altloc_matches(b, AltLocRequest::Exact(a)));

    let mut parts = one_atom_parts([0.0; 3]);
    parts.residues[0] = BioResidueRow::new(
        BioChainId::new(0),
        span(0, 2),
        residue_name("ALA"),
        ResidueInfoKind::Aa,
        EntityKind::Polymer,
        Some(BioEntityId::new(0)),
        Some(b'A'),
        ResidueSourceIds::new(None, None, None, Some("LABEL_LONG".to_owned()), None).unwrap(),
        BioSiftsUnpResidue::default(),
    );
    parts.atoms = vec![
        BioAtomRow::new(
            BioResidueId::new(0),
            atom_name(b" CA "),
            Element::C,
            a,
            0,
            BioCalcFlag::NotSet,
            1.0,
            20.0,
            [0.0; 6],
            -1,
            0.0,
            AtomSourceIds::new(None),
        ),
        BioAtomRow::new(
            BioResidueId::new(0),
            atom_name(b" CA "),
            Element::C,
            b,
            0,
            BioCalcFlag::NotSet,
            1.0,
            20.0,
            [0.0; 6],
            -1,
            0.0,
            AtomSourceIds::new(None),
        ),
    ];
    parts.coordinates = BioCoordinateBlock::new(vec![[0.0; 3]; 2]);
    let structure = BioStructure::from_parts(parts).unwrap();
    assert_eq!(
        structure
            .find_atom(
                BioResidueId::new(0),
                atom_name(b" CA "),
                AltLocRequest::Any,
                Some(Element::C)
            )
            .unwrap()
            .0,
        BioAtomId::new(0)
    );
    assert_eq!(
        structure
            .atom_by_altloc(BioResidueId::new(0), atom_name(b" CA "), b)
            .unwrap()
            .0,
        BioAtomId::new(1)
    );
    assert!(matches!(
        structure.atom_by_altloc(BioResidueId::new(0), atom_name(b" N  "), None),
        Err(BioStructureError::AtomNotFound)
    ));
}

#[test]
fn bio_hierarchy_entity_lookup_preserves_exact_strings_order_and_metadata() {
    let mut parts = one_atom_parts([0.0; 3]);
    parts.entities.push(entity("1", &["SECOND"]));
    let structure = BioStructure::from_parts(parts).unwrap();
    assert_eq!(BioEntityRow::first_mon("ALA,GLY"), "ALA");
    assert_eq!(BioEntityRow::first_mon(",GLY"), "");
    assert_eq!(BioEntityRow::first_mon("ALA"), "ALA");
    assert_eq!(structure.find_entity("1").unwrap().0, BioEntityId::new(0));
    assert_eq!(
        structure.find_entity_of_subchain("LABEL_LONG").unwrap().0,
        BioEntityId::new(0)
    );
    assert!(structure.find_entity_of_subchain("").is_none());
    assert_eq!(structure.entities()[0].full_sequence(), ["ALA,GLY", "SER"]);
    assert_eq!(structure.entities()[0].dbrefs()[0].db_name, "UNP");
    assert_eq!(structure.entities()[0].sifts_unp_accessions(), ["P00001"]);
}

#[test]
fn bio_hierarchy_transform_and_crystal_branches_match_gemmi_order() {
    let first = BioTransform::new(
        [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]],
        [1.0, 2.0, 3.0],
    );
    let second = BioTransform::new(
        [[1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [5.0, 0.0, -1.0],
    );
    let point = [1.0, 2.0, 3.0];
    assert_eq!(
        first.combine(&second).apply(point),
        first.apply(second.apply(point))
    );

    let defaults = BioCrystalCell::default();
    assert_eq!([defaults.a, defaults.b, defaults.c], [1.0; 3]);
    let mut orthogonal = BioCrystalInfo::new(
        BioCrystalCell {
            a: 10.0,
            b: 20.0,
            c: 30.0,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        },
        Some("P 1".to_owned()),
        None,
        BioTransform::identity(),
        BioTransform::identity(),
        false,
        0,
        vec![BioTransform::identity()],
    );
    orthogonal.calculate_properties().unwrap();
    assert_eq!(orthogonal.volume(), 6000.0);
    assert!(orthogonal.is_crystal());
    assert_eq!(
        orthogonal.orthogonal().matrix()[0][1].to_bits(),
        0.0_f64.to_bits()
    );

    let preserved = BioTransform::new([[7.0; 3]; 3], [8.0; 3]);
    let mut explicit = BioCrystalInfo::new(
        BioCrystalCell {
            a: 2.0,
            b: 3.0,
            c: 4.0,
            alpha: 80.0,
            beta: 90.0,
            gamma: 100.0,
        },
        None,
        None,
        preserved,
        preserved,
        true,
        2,
        Vec::new(),
    );
    explicit.calculate_properties().unwrap();
    assert_eq!(*explicit.orthogonal(), preserved);
    assert_eq!(*explicit.fractional(), preserved);

    for angle in [0.0, -0.0] {
        let mut impossible = BioCrystalInfo::new(
            BioCrystalCell {
                a: 1.0,
                b: 1.0,
                c: 1.0,
                alpha: angle,
                beta: 90.0,
                gamma: 90.0,
            },
            None,
            None,
            BioTransform::identity(),
            BioTransform::identity(),
            false,
            0,
            Vec::new(),
        );
        assert_eq!(
            impossible.calculate_properties(),
            Err(BioStructureError::ImpossibleCrystalAngle)
        );
    }
    let mut source_rounding = BioCrystalInfo::new(
        BioCrystalCell {
            a: 1.0,
            b: 1.0,
            c: 1.0,
            alpha: 180.0,
            beta: 90.0,
            gamma: 90.0,
        },
        None,
        None,
        BioTransform::identity(),
        BioTransform::identity(),
        false,
        0,
        Vec::new(),
    );
    assert!(source_rounding.calculate_properties().is_ok());
}

#[test]
fn bio_hierarchy_assembly_ncs_and_reference_validation_preserve_order() {
    let transform = BioTransform::identity();
    let operator = BioAssemblyOperator::new(Some("1".to_owned()), None, transform);
    let assembly = BioAssembly::new(
        "bio1".to_owned(),
        true,
        false,
        BioAssemblySpecialKind::CompletePoint,
        2,
        "dimer".to_owned(),
        "author".to_owned(),
        f64::NAN,
        -0.0,
        1.5,
        vec![BioAssemblyGenerator::new(
            vec!["A".to_owned(), "A".to_owned()],
            vec!["LABEL_LONG".to_owned()],
            vec![operator.clone(), operator],
        )],
    );
    let mut parts = one_atom_parts([0.0; 3]);
    parts
        .ncs_operators
        .push(BioNcsOperator::new("ncs1".to_owned(), true, transform));
    parts.assemblies.push(assembly);
    let structure = BioStructure::from_parts(parts).unwrap();
    assert_eq!(structure.assemblies()[0].generators[0].chains, ["A", "A"]);
    assert!(structure.assemblies()[0].buried_surface_area.is_nan());
    assert_eq!(
        structure.assemblies()[0].surface_area.to_bits(),
        (-0.0_f64).to_bits()
    );
    assert_eq!(structure.ncs_operators()[0].id, "ncs1");

    let mut bad = one_atom_parts([0.0; 3]);
    bad.assemblies.push(BioAssembly::new(
        "bad".to_owned(),
        false,
        false,
        BioAssemblySpecialKind::NotApplicable,
        0,
        String::new(),
        String::new(),
        f64::NAN,
        f64::NAN,
        f64::NAN,
        vec![BioAssemblyGenerator::new(
            vec!["missing".to_owned()],
            Vec::new(),
            Vec::new(),
        )],
    ));
    assert!(matches!(
        BioStructure::from_parts(bad),
        Err(BioStructureError::AssemblyReferenceMissing { kind: "chain", .. })
    ));
}

#[test]
fn bio_hierarchy_all_source_enum_states_remain_distinct() {
    assert_eq!(
        [
            BioCoordinateFormat::Unknown,
            BioCoordinateFormat::Detect,
            BioCoordinateFormat::Pdb,
            BioCoordinateFormat::Mmcif,
            BioCoordinateFormat::Mmjson,
            BioCoordinateFormat::ChemComp
        ]
        .len(),
        6
    );
    assert_eq!(
        [
            BioCalcFlag::NotSet,
            BioCalcFlag::NoHydrogen,
            BioCalcFlag::Determined,
            BioCalcFlag::Calculated,
            BioCalcFlag::Dummy
        ]
        .len(),
        5
    );
    assert_eq!(
        [
            EntityKind::Unknown,
            EntityKind::Polymer,
            EntityKind::NonPolymer,
            EntityKind::Branched,
            EntityKind::Water
        ]
        .len(),
        5
    );
    assert_eq!(
        [
            PolymerKind::Unknown,
            PolymerKind::PeptideL,
            PolymerKind::PeptideD,
            PolymerKind::Dna,
            PolymerKind::Rna,
            PolymerKind::DnaRnaHybrid,
            PolymerKind::SaccharideD,
            PolymerKind::SaccharideL,
            PolymerKind::Pna,
            PolymerKind::CyclicPseudoPeptide,
            PolymerKind::Other
        ]
        .len(),
        11
    );
}
