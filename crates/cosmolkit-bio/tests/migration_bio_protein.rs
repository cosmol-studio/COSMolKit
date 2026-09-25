use cosmolkit_bio::{
    AltLocLabel, AtomName, AtomSourceIds, BioAssembly, BioAssemblySpecialKind, BioAtomId,
    BioAtomRow, BioCalcFlag, BioChainId, BioChainRow, BioCoordinateBlock, BioCoordinateFormat,
    BioCrystalCell, BioCrystalInfo, BioEntityId, BioEntityRow, BioModelId, BioModelRow,
    BioNcsOperator, BioResidueId, BioResidueRow, BioRowSpan, BioSiftsUnpResidue, BioStructure,
    BioStructureError, BioStructureParts, BioTransform, ChainKind, ChainSourceIds, EntityKind,
    EntitySourceIds, PdbAtomSerial, PdbChainId, PdbSeqId, PolymerKind, ProteinProjectionError,
    ResidueCode, ResidueInfoKind, ResidueKind, ResidueName, ResidueSourceIds,
};
use cosmolkit_types::Element;

fn span<I>(start: u32, len: u32) -> BioRowSpan<I> {
    BioRowSpan::new(start, len).unwrap()
}

fn atom_name(value: &[u8; 4]) -> AtomName {
    AtomName::from_ascii(*value).unwrap()
}

fn residue_name(value: &str) -> ResidueName {
    ResidueName::from_ascii(value.as_bytes()).unwrap()
}

fn chain_source(auth: &[u8], label: &str) -> ChainSourceIds {
    ChainSourceIds::new(
        Some(PdbChainId::from_ascii(auth).unwrap()),
        Some(label.to_owned()),
    )
}

fn entity(source_id: &str, sequence: &str, subchain: &str) -> BioEntityRow {
    BioEntityRow::new(
        EntityKind::Polymer,
        PolymerKind::PeptideL,
        true,
        vec![sequence.to_owned()],
        Vec::new(),
        vec![format!("UNP-{source_id}")],
        vec![subchain.to_owned()],
        EntitySourceIds::new(source_id.to_owned()),
    )
}

#[allow(clippy::too_many_arguments)]
fn residue(
    chain: u32,
    atom_start: u32,
    atom_len: u32,
    name: &str,
    info_kind: ResidueInfoKind,
    entity_id: Option<u32>,
    label_seq_id: Option<i32>,
    subchain: Option<&str>,
) -> BioResidueRow {
    BioResidueRow::new(
        BioChainId::new(chain),
        span(atom_start, atom_len),
        residue_name(name),
        info_kind,
        EntityKind::Polymer,
        entity_id.map(BioEntityId::new),
        Some(b'A'),
        ResidueSourceIds::new(
            Some(PdbSeqId::new(40 + chain as i32, Some(b'B'))),
            label_seq_id,
            Some(*b"SEG "),
            subchain.map(str::to_owned),
            Some(format!("entity-{chain}")),
        )
        .unwrap(),
        BioSiftsUnpResidue::new(Some(b'X'), 1, 321),
    )
}

fn atom(residue_id: u32, altloc: Option<u8>, serial: i32) -> BioAtomRow {
    BioAtomRow::new(
        BioResidueId::new(residue_id),
        atom_name(b" CA "),
        Element::C,
        None,
        altloc.map(AltLocLabel::new),
        -3,
        BioCalcFlag::Calculated,
        0.625,
        17.25,
        [1.0, -2.0, 3.0, -4.0, 5.0, -6.0],
        9,
        0.375,
        AtomSourceIds::new(Some(PdbAtomSerial::new(serial))),
    )
}

fn structure(parts: BioStructureParts) -> BioStructure {
    BioStructure::from_parts(parts).unwrap()
}

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

fn one_residue_parts(kind: ResidueInfoKind, name: &str) -> BioStructureParts {
    BioStructureParts {
        input_format: BioCoordinateFormat::Pdb,
        models: vec![BioModelRow::new(span(0, 1), Some(7))],
        chains: vec![BioChainRow::new(
            BioModelId::new(0),
            None,
            span(0, 1),
            ChainKind::Mixed,
            chain_source(b"A", "label-A"),
        )],
        residues: vec![residue(
            0,
            0,
            0,
            name,
            kind,
            None,
            Some(11),
            Some("label-A"),
        )],
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
fn protein_projection_classifies_exactly_four_gemmi_amino_acid_kinds() {
    let included = [
        (ResidueInfoKind::Aa, "ALA"),
        (ResidueInfoKind::Aad, "DAL"),
        (ResidueInfoKind::Paa, "PAA"),
        (ResidueInfoKind::Maa, "MAA"),
    ];
    for (kind, name) in included {
        let protein = structure(one_residue_parts(kind, name)).protein().unwrap();
        assert_eq!(protein.num_models(), 1, "included kind {kind:?}");
        assert_eq!(protein.num_chains(), 1, "included kind {kind:?}");
        assert_eq!(protein.num_residues(), 1, "included kind {kind:?}");
    }

    let excluded = [
        ResidueInfoKind::Unknown,
        ResidueInfoKind::Rna,
        ResidueInfoKind::Dna,
        ResidueInfoKind::Buf,
        ResidueInfoKind::Hoh,
        ResidueInfoKind::Pyr,
        ResidueInfoKind::Ket,
        ResidueInfoKind::Els,
    ];
    for kind in excluded {
        let protein = structure(one_residue_parts(kind, "ZZZ")).protein().unwrap();
        assert_eq!(protein.num_models(), 0, "excluded kind {kind:?}");
        assert_eq!(protein.num_chains(), 0, "excluded kind {kind:?}");
        assert_eq!(protein.num_residues(), 0, "excluded kind {kind:?}");
        assert_eq!(protein.num_atoms(), 0, "excluded kind {kind:?}");
    }
}

#[test]
fn protein_projection_rebuilds_dense_hierarchy_and_retains_empty_residues() {
    let parts = BioStructureParts {
        input_format: BioCoordinateFormat::Mmcif,
        models: vec![
            BioModelRow::new(span(0, 2), Some(10)),
            BioModelRow::new(span(2, 1), Some(20)),
            BioModelRow::new(span(3, 1), Some(30)),
        ],
        chains: vec![
            BioChainRow::new(
                BioModelId::new(0),
                None,
                span(0, 1),
                ChainKind::Dna,
                chain_source(b"D", "dna"),
            ),
            BioChainRow::new(
                BioModelId::new(0),
                None,
                span(1, 3),
                ChainKind::Mixed,
                chain_source(b"A", "protein-a"),
            ),
            BioChainRow::new(
                BioModelId::new(1),
                None,
                span(4, 1),
                ChainKind::WaterOnly,
                chain_source(b"W", "water"),
            ),
            BioChainRow::new(
                BioModelId::new(2),
                None,
                span(5, 1),
                ChainKind::Unknown,
                chain_source(b"A", "protein-a"),
            ),
        ],
        residues: vec![
            residue(0, 0, 0, "DA", ResidueInfoKind::Dna, None, None, Some("dna")),
            residue(
                1,
                0,
                1,
                "ALA",
                ResidueInfoKind::Aa,
                None,
                Some(1),
                Some("protein-a"),
            ),
            residue(
                1,
                1,
                1,
                "HOH",
                ResidueInfoKind::Hoh,
                None,
                None,
                Some("protein-a"),
            ),
            residue(
                1,
                2,
                0,
                "DAL",
                ResidueInfoKind::Aad,
                None,
                Some(2),
                Some("protein-a"),
            ),
            residue(
                2,
                2,
                0,
                "HOH",
                ResidueInfoKind::Hoh,
                None,
                None,
                Some("water"),
            ),
            residue(
                3,
                2,
                0,
                "PAA",
                ResidueInfoKind::Paa,
                None,
                Some(3),
                Some("protein-a"),
            ),
        ],
        atoms: vec![atom(1, None, 101), atom(2, None, 102)],
        entities: Vec::new(),
        connections: Vec::new(),
        cispeps: Vec::new(),
        mod_residues: Vec::new(),
        helices: Vec::new(),
        sheets: Vec::new(),
        metadata: Default::default(),
        source_state: Default::default(),
        coordinates: BioCoordinateBlock::new(vec![[1.0, 2.0, 3.0], [9.0, 9.0, 9.0]]),
        crystal: None,
        ncs_operators: Vec::new(),
        assemblies: Vec::new(),
    };
    let protein = structure(parts).protein().unwrap();

    assert_eq!((protein.num_models(), protein.num_chains()), (2, 2));
    assert_eq!((protein.num_residues(), protein.num_atoms()), (3, 1));
    let chains = protein.chains();
    assert_eq!(
        chains
            .iter()
            .map(|chain| chain.id().value())
            .collect::<Vec<_>>(),
        vec![0, 1]
    );
    assert_eq!(chains[0].row().model_id(), BioModelId::new(0));
    assert_eq!(chains[1].row().model_id(), BioModelId::new(1));
    assert_eq!(chains[0].source().auth_chain_id().unwrap().as_str(), "A");
    assert_eq!(chains[1].source().auth_chain_id().unwrap().as_str(), "A");
    assert_eq!(chains[0].kind(), ChainKind::Protein);
    assert_eq!(chains[0].residues().len(), 2);
    assert_eq!(chains[0].atoms().len(), 1);
    assert_eq!(chains[1].residues().len(), 1);
    assert!(chains[1].atoms().is_empty());

    let residues = protein.residues();
    assert_eq!(
        residues
            .iter()
            .map(|row| row.id().value())
            .collect::<Vec<_>>(),
        vec![0, 1, 2]
    );
    assert_eq!(residues[0].row().chain_id(), BioChainId::new(0));
    assert_eq!(residues[1].row().chain_id(), BioChainId::new(0));
    assert_eq!(residues[2].row().chain_id(), BioChainId::new(1));
    assert_eq!(residues[0].row().atom_span(), span(0, 1));
    assert_eq!(residues[1].row().atom_span(), span(1, 0));
    assert_eq!(residues[2].row().atom_span(), span(1, 0));
    assert_eq!(protein.atoms()[0].row().residue_id(), BioResidueId::new(0));
}

#[test]
fn protein_projection_preserves_atom_state_altloc_sources_and_coordinate_bits() {
    let nan = f64::from_bits(0x7ff8_0000_0000_0042);
    let mut parts = one_residue_parts(ResidueInfoKind::Aa, "ALA");
    parts.residues[0] = residue(
        0,
        0,
        3,
        "ALA",
        ResidueInfoKind::Aa,
        None,
        Some(7),
        Some("label-A"),
    );
    parts.atoms = vec![
        atom(0, Some(b'B'), 111),
        atom(0, None, 112),
        atom(0, Some(b'A'), 113),
    ];
    parts.coordinates = BioCoordinateBlock::new(vec![
        [-0.0, nan, f64::INFINITY],
        [f64::NEG_INFINITY, 0.0, 0.125],
        [4.0, 5.0, 6.0],
    ]);
    let protein = structure(parts).protein().unwrap();
    let atoms = protein.atoms();
    assert_eq!(atoms.len(), 3);
    assert_eq!(atoms[0].id(), BioAtomId::new(0));
    assert_eq!(atoms[0].name().as_str(), " CA ");
    assert_eq!(atoms[0].element(), Element::C);
    assert_eq!(atoms[0].altloc().unwrap().value(), b'B');
    assert_eq!(atoms[0].row().formal_charge(), -3);
    assert_eq!(atoms[0].row().calc_flag(), BioCalcFlag::Calculated);
    assert_eq!(atoms[0].row().occupancy().to_bits(), 0.625_f64.to_bits());
    assert_eq!(atoms[0].row().b_iso().to_bits(), 17.25_f64.to_bits());
    assert_eq!(atoms[0].row().anisou(), &[1.0, -2.0, 3.0, -4.0, 5.0, -6.0]);
    assert_eq!(atoms[0].row().tls_group_id(), 9);
    assert_eq!(atoms[0].row().fraction().to_bits(), 0.375_f64.to_bits());
    assert_eq!(atoms[0].row().source().serial().unwrap().value(), 111);
    let first = atoms[0].position();
    assert_eq!(first[0].to_bits(), (-0.0_f64).to_bits());
    assert_eq!(first[1].to_bits(), nan.to_bits());
    assert_eq!(first[2].to_bits(), f64::INFINITY.to_bits());
    assert_eq!(
        atoms[1].position()[0].to_bits(),
        f64::NEG_INFINITY.to_bits()
    );
    assert_eq!(atoms[2].position()[2].to_bits(), 6.0_f64.to_bits());
}

#[test]
fn protein_projection_preserves_isotope_state_through_dense_atom_and_residue_remapping() {
    let mut parts = one_residue_parts(ResidueInfoKind::Aa, "ALA");
    parts.chains[0] = BioChainRow::new(
        BioModelId::new(0),
        None,
        span(0, 2),
        ChainKind::Mixed,
        chain_source(b"A", "label-A"),
    );
    parts.residues = vec![
        residue(
            0,
            0,
            1,
            "HOH",
            ResidueInfoKind::Hoh,
            None,
            None,
            Some("label-A"),
        ),
        residue(
            0,
            1,
            2,
            "ALA",
            ResidueInfoKind::Aa,
            None,
            Some(7),
            Some("label-A"),
        ),
    ];
    let isotope_hydrogen = |name: &[u8; 4],
                            isotope_mass_number: Option<u16>,
                            altloc: Option<u8>,
                            formal_charge: i8,
                            occupancy: f64,
                            b_iso: f64,
                            anisou: [f64; 6],
                            tls_group_id: i16,
                            fraction: f64,
                            serial: i32| {
        BioAtomRow::new(
            BioResidueId::new(1),
            atom_name(name),
            Element::H,
            isotope_mass_number,
            altloc.map(AltLocLabel::new),
            formal_charge,
            BioCalcFlag::Calculated,
            occupancy,
            b_iso,
            anisou,
            tls_group_id,
            fraction,
            AtomSourceIds::new(Some(PdbAtomSerial::new(serial))),
        )
    };
    parts.atoms = vec![
        atom(0, None, 110),
        isotope_hydrogen(
            b" H1 ",
            None,
            Some(b'B'),
            0,
            0.625,
            17.25,
            [1.0, -2.0, 3.0, -4.0, 5.0, -6.0],
            9,
            0.375,
            111,
        ),
        isotope_hydrogen(
            b" D1 ",
            Some(2),
            Some(b'A'),
            1,
            0.875,
            23.5,
            [-1.0, 2.0, -3.0, 4.0, -5.0, 6.0],
            3,
            0.25,
            112,
        ),
    ];
    parts.coordinates =
        BioCoordinateBlock::new(vec![[9.0, 8.0, 7.0], [-0.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);

    let source = structure(parts);
    let source_before = source.clone();
    let protein = source.protein().unwrap();

    assert_eq!(source, source_before);
    assert_eq!(protein.num_residues(), 1);
    assert_eq!(protein.num_atoms(), 2);
    assert_eq!(protein.residues()[0].id(), BioResidueId::new(0));

    let atoms = protein.atoms();
    assert_eq!(atoms[0].id(), BioAtomId::new(0));
    assert_eq!(atoms[1].id(), BioAtomId::new(1));
    assert_eq!(atoms[0].row().residue_id(), BioResidueId::new(0));
    assert_eq!(atoms[1].row().residue_id(), BioResidueId::new(0));
    assert_eq!(atoms[0].row().element(), Element::H);
    assert_eq!(atoms[1].row().element(), Element::H);
    assert_eq!(atoms[0].row().isotope_mass_number(), None);
    assert_eq!(atoms[1].row().isotope_mass_number(), Some(2));
    assert_eq!(atoms[0].name().as_str(), " H1 ");
    assert_eq!(atoms[1].name().as_str(), " D1 ");
    assert_eq!(atoms[0].altloc().unwrap().value(), b'B');
    assert_eq!(atoms[1].altloc().unwrap().value(), b'A');
    assert_eq!(atoms[0].row().formal_charge(), 0);
    assert_eq!(atoms[1].row().formal_charge(), 1);
    assert_eq!(atoms[0].row().occupancy().to_bits(), 0.625_f64.to_bits());
    assert_eq!(atoms[1].row().occupancy().to_bits(), 0.875_f64.to_bits());
    assert_eq!(atoms[0].row().b_iso().to_bits(), 17.25_f64.to_bits());
    assert_eq!(atoms[1].row().b_iso().to_bits(), 23.5_f64.to_bits());
    assert_eq!(atoms[0].row().anisou(), &[1.0, -2.0, 3.0, -4.0, 5.0, -6.0]);
    assert_eq!(atoms[1].row().anisou(), &[-1.0, 2.0, -3.0, 4.0, -5.0, 6.0]);
    assert_eq!(atoms[0].row().tls_group_id(), 9);
    assert_eq!(atoms[1].row().tls_group_id(), 3);
    assert_eq!(atoms[0].row().fraction().to_bits(), 0.375_f64.to_bits());
    assert_eq!(atoms[1].row().fraction().to_bits(), 0.25_f64.to_bits());
    assert_eq!(atoms[0].row().source().serial().unwrap().value(), 111);
    assert_eq!(atoms[1].row().source().serial().unwrap().value(), 112);
    assert_eq!(atoms[0].position()[0].to_bits(), (-0.0_f64).to_bits());
    assert_eq!(atoms[0].position()[1].to_bits(), 2.0_f64.to_bits());
    assert_eq!(atoms[0].position()[2].to_bits(), 3.0_f64.to_bits());
    assert_eq!(atoms[1].position(), [4.0, 5.0, 6.0]);
}

#[test]
fn protein_projection_remaps_entities_and_preserves_auth_label_distinctions() {
    let mut parts = one_residue_parts(ResidueInfoKind::Aa, "ALA");
    parts.entities = vec![
        entity("unused", "XXX", "unused-label"),
        entity("kept", "ALA,GLY", "label-A"),
    ];
    parts.chains[0] = BioChainRow::new(
        BioModelId::new(0),
        Some(BioEntityId::new(1)),
        span(0, 1),
        ChainKind::Mixed,
        chain_source(b"AUTH", "label-A"),
    );
    parts.residues[0] = residue(
        0,
        0,
        0,
        "ALA",
        ResidueInfoKind::Aa,
        Some(1),
        Some(5),
        Some("label-A"),
    );
    let source = structure(parts.clone());
    let protein = source.protein().unwrap();
    let chain = protein.chain(0).unwrap();
    let residue = chain.residues()[0];
    assert_eq!(chain.row().entity_id(), Some(BioEntityId::new(0)));
    assert_eq!(residue.row().entity_id(), Some(BioEntityId::new(0)));
    assert_eq!(chain.source().auth_chain_id().unwrap().as_str(), "AUTH");
    assert_eq!(chain.source().label_asym_id(), Some("label-A"));
    assert_eq!(residue.row().source().subchain_id(), Some("label-A"));

    let mut changed_unused = parts.clone();
    changed_unused.entities[0] = entity("different-unused", "YYY", "other");
    assert_eq!(protein, structure(changed_unused).protein().unwrap());
    let mut changed_kept = parts;
    changed_kept.entities[1] = entity("kept", "SER", "label-A");
    assert_ne!(protein, structure(changed_kept).protein().unwrap());
}

#[test]
fn protein_projection_preserves_global_frame_state_and_drops_assemblies() {
    let mut base = one_residue_parts(ResidueInfoKind::Aa, "ALA");
    let crystal = BioCrystalInfo::new(
        BioCrystalCell {
            a: 12.0,
            b: 13.0,
            c: 14.0,
            alpha: 90.0,
            beta: 90.0,
            gamma: 120.0,
        },
        Some("P 1".to_owned()),
        Some("4".to_owned()),
        BioTransform::identity(),
        BioTransform::identity(),
        true,
        2,
        vec![BioTransform::identity()],
    );
    base.crystal = Some(crystal.clone());
    base.ncs_operators = vec![BioNcsOperator::new(
        "ncs-1".to_owned(),
        true,
        BioTransform::new(
            [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]],
            [1.0, 2.0, 3.0],
        ),
    )];
    base.assemblies = vec![BioAssembly::new(
        "assembly-1".to_owned(),
        true,
        false,
        BioAssemblySpecialKind::NotApplicable,
        1,
        "details".to_owned(),
        "software".to_owned(),
        1.0,
        2.0,
        3.0,
        Vec::new(),
    )];
    let projected = structure(base.clone()).protein().unwrap();

    let mut changed_format = base.clone();
    changed_format.input_format = BioCoordinateFormat::ChemComp;
    assert_ne!(projected, structure(changed_format).protein().unwrap());
    let mut changed_crystal = base.clone();
    changed_crystal.crystal = None;
    assert_ne!(projected, structure(changed_crystal).protein().unwrap());
    let mut changed_ncs = base.clone();
    changed_ncs.ncs_operators.clear();
    assert_ne!(projected, structure(changed_ncs).protein().unwrap());
    let mut changed_assembly = base;
    changed_assembly.assemblies.clear();
    assert_eq!(projected, structure(changed_assembly).protein().unwrap());
}

#[test]
fn protein_read_only_views_cover_counts_indices_parents_rows_and_residue_info() {
    let mut parts = one_residue_parts(ResidueInfoKind::Aa, "ALA");
    parts.residues[0] = residue(
        0,
        0,
        1,
        "ALA",
        ResidueInfoKind::Aa,
        None,
        Some(1),
        Some("label-A"),
    );
    parts.atoms = vec![atom(0, Some(b'A'), 17)];
    parts.coordinates = BioCoordinateBlock::new(vec![[7.0, 8.0, 9.0]]);
    let protein = structure(parts).protein().unwrap();
    assert_eq!((protein.num_models(), protein.num_chains()), (1, 1));
    assert_eq!((protein.num_residues(), protein.num_atoms()), (1, 1));
    assert!(protein.chain(1).is_none());

    let chain = protein.chain(0).unwrap();
    let chain_again = protein.chains()[0];
    assert!(std::ptr::eq(chain.row(), chain_again.row()));
    let residue = protein.residues()[0];
    assert!(std::ptr::eq(residue.row(), chain.residues()[0].row()));
    assert!(std::ptr::eq(residue.chain().row(), chain.row()));
    assert_eq!(residue.name().as_str(), "ALA");
    assert_eq!(residue.kind(), ResidueKind::AminoAcid);
    assert_eq!(residue.info().kind, ResidueInfoKind::Aa);
    assert_eq!(residue.code(), ResidueCode::ALA);
    assert_eq!(residue.one_letter_code(), 'A');
    assert_eq!(residue.fasta_code(), 'A');
    assert!(residue.is_standard());

    let atom = protein.atoms()[0];
    assert!(std::ptr::eq(atom.row(), residue.atoms()[0].row()));
    assert!(std::ptr::eq(atom.residue().row(), residue.row()));
    assert_eq!(atom.position(), [7.0, 8.0, 9.0]);
}

#[test]
fn protein_projection_empty_source_immutability_determinism_and_error_shape() {
    let empty_source = structure(empty_parts());
    let empty = empty_source.protein().unwrap();
    assert_eq!(
        (
            empty.num_models(),
            empty.num_chains(),
            empty.num_residues(),
            empty.num_atoms()
        ),
        (0, 0, 0, 0)
    );

    let excluded_source = structure(one_residue_parts(ResidueInfoKind::Hoh, "HOH"));
    let excluded = excluded_source.protein().unwrap();
    assert_eq!(
        (
            excluded.num_models(),
            excluded.num_chains(),
            excluded.num_residues(),
            excluded.num_atoms()
        ),
        (0, 0, 0, 0)
    );
    let mut empty_pdb_parts = empty_parts();
    empty_pdb_parts.input_format = BioCoordinateFormat::Pdb;
    assert_eq!(excluded, structure(empty_pdb_parts).protein().unwrap());
    assert_ne!(excluded, empty);

    let source = structure(one_residue_parts(ResidueInfoKind::Aa, "ALA"));
    let before = source.clone();
    let first = source.protein().unwrap();
    let second = source.protein().unwrap();
    assert_eq!(source, before);
    assert_eq!(first, second);

    let cause = BioStructureError::CoordinateCountMismatch {
        atom_count: 1,
        coordinate_count: 0,
    };
    let error = ProteinProjectionError::from(cause.clone());
    assert_eq!(error, ProteinProjectionError::Structure(cause));
    assert!(error.to_string().contains("protein projection failed"));
    assert!(std::error::Error::source(&error).is_some());
}
