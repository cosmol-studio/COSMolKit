use cosmolkit::{
    Conformer2D, Conformer3D, CoordinateBlock, Molecule, MoleculeProperties, OperationError,
};

#[test]
fn selected_bonds_cut_value_and_in_place_share_chemistry_and_isolate_peers() {
    use cosmolkit::{BondCut, BondId};
    let source = Molecule::from_smiles("CCO").unwrap();
    let cuts = [BondCut::at_bond(BondId::new(1), Some(101), Some(202)).unwrap()];
    let value = source.with_bonds_cut(&cuts).unwrap();
    let mut target = source.clone();
    target.cut_bonds_(&cuts).unwrap();
    assert_eq!(value, target);
    assert_eq!(value.to_smiles().unwrap().as_bytes(), b"[101*]CC.[202*]O");
    assert_eq!(source.to_smiles().unwrap().as_bytes(), b"CCO");
    assert_eq!(target.fragments().unwrap().len(), 2);
}

#[test]
fn invalid_later_selection_changes_neither_cut_receiver_nor_peer() {
    use cosmolkit::{BondCut, BondCutError, BondId};
    let mut target = Molecule::from_smiles("CCO").unwrap();
    let untouched = target.clone();
    let cuts = [
        BondCut::at_bond(BondId::new(0), None, None).unwrap(),
        BondCut::at_bond(BondId::new(2), None, None).unwrap(),
    ];
    let expected = OperationError::BondCut(BondCutError::BondOutOfRange {
        bond: BondId::new(2),
        count: 2,
    });
    assert_eq!(target.with_bonds_cut(&cuts).unwrap_err(), expected);
    assert_eq!(target.cut_bonds_(&cuts).unwrap_err(), expected);
    assert_eq!(target, untouched);
    assert_eq!(target.to_smiles().unwrap().as_bytes(), b"CCO");
    assert_eq!(untouched.to_smiles().unwrap().as_bytes(), b"CCO");
}

#[test]
fn sanitized_fragment_retains_the_ring_cache_used_by_writer_stereo() {
    // RDKit 2026.03.6 MolOps::getTheFrags -> sanitizeMol retains SymmSSSR.
    // Replacing it with FastFindRings changes macrocycle cis/trans writing.
    let input = r"OCCOCCOCCOCCO[Si]1(OCCOCCOCCOCCO)n2c3c4ccccc4c2/N=C2\N=C(/N=c4/c5ccccc5/c(n41)=N/C1=N/C(=N\3)c3ccccc31)c1ccccc12";
    let expected = r"OCCOCCOCCOCCO[Si]1(OCCOCCOCCOCCO)n2c3c4ccccc4c2/N=C2N=C(/N=c4/c5ccccc5/c(n41)=N/C1=NC(=N\3)/c3ccccc31)c1ccccc1\2";
    for suffix in ["", ".[Na+]"] {
        let source = Molecule::from_smiles(&format!("{input}{suffix}")).unwrap();
        let before = source.clone();
        let fragment = source.fragments().unwrap().remove(0);
        assert_eq!(
            fragment.to_smiles().unwrap().as_bytes(),
            expected.as_bytes()
        );
        assert_eq!(
            source
                .largest_fragment()
                .unwrap()
                .to_smiles()
                .unwrap()
                .as_bytes(),
            expected.as_bytes()
        );
        assert_eq!(source, before);
    }
}

#[test]
fn fragments_preserve_component_order_and_source_value() {
    let source = Molecule::from_smiles("CC.O.[Na+]").unwrap();
    let peer = source.clone();
    let fragments = source.fragments().unwrap();
    let values: Vec<_> = fragments.iter().map(|m| m.to_smiles().unwrap()).collect();
    assert_eq!(values, ["CC".into(), "O".into(), "[Na+]".into()]);
    assert_eq!(source, peer);
    assert_eq!(
        source.largest_fragment().unwrap().to_smiles().unwrap(),
        "CC".into()
    );
}

#[test]
fn largest_fragment_retains_last_tie_and_empty_error() {
    assert_eq!(
        Molecule::from_smiles("CC.OO")
            .unwrap()
            .largest_fragment()
            .unwrap()
            .to_smiles()
            .unwrap(),
        "OO".into()
    );
    assert!(Molecule::new().fragments().unwrap().is_empty());
    assert!(matches!(
        Molecule::new().largest_fragment(),
        Err(OperationError::EmptyFragments)
    ));
}

#[test]
fn fragments_preserve_stereo_and_project_all_coordinate_rows() {
    let graph = Molecule::from_smiles("N[C@@H](C)O.Cl")
        .unwrap()
        .topology()
        .clone();
    let n = graph.atoms.len();
    let source = Molecule::from_parts(
        graph,
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                8,
                (0..n).map(|i| [i as f64, -0.0]).collect(),
            )],
            conformers_3d: vec![Conformer3D::new(
                19,
                (0..n).map(|i| [i as f64, 2.0, 3.0]).collect(),
                true,
            )],
            ..Default::default()
        },
        MoleculeProperties::default(),
    )
    .unwrap();
    let peer = source.clone();
    let fragments = source.fragments().unwrap();
    assert_eq!(fragments.len(), 2);
    assert_eq!(fragments[0].to_smiles().unwrap(), "C[C@H](N)O".into());
    for (part, plain) in source
        .fragments_with_atom_mapping()
        .unwrap()
        .iter()
        .zip(&fragments)
    {
        assert_eq!(part.molecule(), plain);
        assert_eq!(
            part.molecule().to_smiles().unwrap(),
            plain.to_smiles().unwrap()
        );
    }
    assert_eq!(fragments[1].coordinates_2d().unwrap(), &[[4.0, -0.0]]);
    assert_eq!(fragments[1].conformers_3d()[0].id(), 19);
    assert_eq!(
        fragments[1].conformers_3d()[0].coordinates(),
        &[[4.0, 2.0, 3.0]]
    );
    assert_eq!(source, peer);
}

#[test]
fn fragments_sanitize_errors_preserve_source_and_cause() {
    use std::error::Error;
    let source = Molecule::from_smiles_with_params(
        "[CH5].O",
        &cosmolkit::SmilesParseParams {
            sanitize: false,
            ..Default::default()
        },
    )
    .unwrap();
    let peer = source.clone();
    let error = source.fragments().unwrap_err();
    assert!(matches!(error, OperationError::Fragments(_)));
    assert!(error.source().is_some());
    assert_eq!(source, peer);
    assert!(source.largest_fragment().is_err());
    assert_eq!(source, peer);
}

#[test]
fn identical_fragments_scatter_atom_data_into_source_order() {
    use cosmolkit::FragmentExtractionParams;
    let source = Molecule::from_parts(
        Molecule::from_smiles("C1.CC.C1")
            .unwrap()
            .topology()
            .clone(),
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                8,
                (0..4).map(|i| [i as f64 + 0.25, -1.0]).collect(),
            )],
            conformers_3d: vec![Conformer3D::new(
                19,
                (0..4).map(|i| [1.5, 2.5, i as f64]).collect(),
                true,
            )],
            ..Default::default()
        },
        MoleculeProperties::default(),
    )
    .unwrap();
    let parts = source.fragments_with_atom_mapping().unwrap();
    assert_eq!(parts.len(), 2);
    let mut restored = [0; 4];
    for ((part, plain), values) in parts
        .iter()
        .zip(source.fragments().unwrap())
        .zip([[10, 40], [20, 30]])
    {
        assert_eq!(part.molecule(), &plain);
        assert_eq!(
            part.molecule().num_atoms(),
            part.source_atom_indices().len()
        );
        for (atom, value) in part.source_atom_indices().iter().zip(values) {
            restored[atom.index()] = value;
        }
    }
    assert_eq!(restored, [10, 20, 30, 40]);
    for sanitize in [true, false] {
        let mapped = source
            .fragments_with_atom_mapping_with_params(&FragmentExtractionParams { sanitize })
            .unwrap();
        for (part, indices) in mapped.iter().zip([[0, 3], [1, 2]]) {
            assert_eq!(
                part.source_atom_indices(),
                indices.map(cosmolkit::AtomId::new)
            );
            let mol = part.molecule();
            assert_eq!(
                mol.coordinates_2d().unwrap(),
                indices.map(|i| source.coordinates_2d().unwrap()[i])
            );
            assert_eq!(mol.conformers_3d()[0].id(), 19);
            assert_eq!(
                mol.conformers_3d()[0].coordinates(),
                indices.map(|i| source.conformers_3d()[0].coordinates()[i])
            );
        }
    }
    assert!(
        Molecule::new()
            .fragments_with_atom_mapping()
            .unwrap()
            .is_empty()
    );
}

#[test]
fn mapped_extraction_can_leave_components_unsanitized() {
    use cosmolkit::{AtomId, FragmentExtractionParams, SmilesParseParams};
    let source = Molecule::from_smiles_with_params(
        "C.C(F)(F)(F)(F)F",
        &SmilesParseParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        },
    )
    .unwrap();
    let before = source.clone();
    let error = source.fragments().unwrap_err();
    assert!(matches!(&error, OperationError::Fragments(e) if e.component_index() == Some(1)));
    assert_eq!(source.fragments_with_atom_mapping().unwrap_err(), error);

    let parts = source
        .fragments_with_atom_mapping_with_params(&FragmentExtractionParams { sanitize: false })
        .unwrap();
    assert_eq!(parts.len(), 2);
    assert_eq!(parts[0].source_atom_indices(), &[AtomId::new(0)]);
    assert_eq!(
        parts[1].source_atom_indices(),
        [1, 2, 3, 4, 5, 6].map(AtomId::new)
    );
    assert_eq!(parts[1].molecule().num_bonds(), 5);
    #[cfg(feature = "cap-sanitize")]
    assert!(matches!(
        parts[1].molecule().sanitize(),
        Err(OperationError::Sanitize(
            cosmolkit::SanitizeError::Properties { .. }
        ))
    ));
    assert_eq!(source, before);
}

#[test]
fn identical_dummy_labels_map_to_distinct_atoms_in_the_cut_graph() {
    use cosmolkit::{AtomId, BondCut, BondId};
    let cut = BondCut::at_bond(BondId::new(0), Some(9), Some(9)).unwrap();
    let product = Molecule::from_smiles("CC")
        .unwrap()
        .with_bonds_cut(&[cut])
        .unwrap();
    let parts = product.fragments_with_atom_mapping().unwrap();
    assert_eq!(parts.len(), 2);
    for (part, indices) in parts.iter().zip([[0, 3], [1, 2]]) {
        assert_eq!(part.source_atom_indices(), indices.map(AtomId::new));
        assert_eq!(part.molecule().atoms()[1].isotope(), Some(9));
    }
}
