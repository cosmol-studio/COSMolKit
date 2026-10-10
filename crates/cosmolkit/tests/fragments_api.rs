use cosmolkit::{
    Conformer2D, Conformer3D, CoordinateBlock, Molecule, MoleculeProperties, OperationError,
};

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
