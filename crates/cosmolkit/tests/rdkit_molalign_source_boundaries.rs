#![cfg(feature = "cap-alignment")]
use cosmolkit::*;
fn one() -> Molecule {
    let m = Molecule::from_smiles("CC").unwrap();
    Molecule::from_parts(
        m.topology().clone(),
        CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(17, vec![[0., 0., 0.], [1., 0., 0.]], true)],
            ..Default::default()
        },
        m.properties().clone(),
    )
    .unwrap()
}
#[test]
fn empty_atom_indices_select_all_atoms_like_source_wrapper() {
    for count in [1, 2] {
        let source = one();
        let mut coordinates = CoordinateBlock {
            conformers_3d: source.conformers_3d().to_vec(),
            ..Default::default()
        };
        if count == 2 {
            coordinates.conformers_3d.push(Conformer3D::new(
                27,
                vec![[2., 1., 0.], [3., 1., 0.]],
                true,
            ));
        }
        let source = Molecule::from_parts(
            source.topology().clone(),
            coordinates,
            source.properties().clone(),
        )
        .unwrap();
        let params = ConformerAlignmentParameters {
            atom_indices: Some(Vec::new()),
            ..Default::default()
        };
        let (aligned, report) = source.with_aligned_conformers_with_params(&params).unwrap();
        assert_eq!(report.rmsds, vec![0.; count - 1]);
        for conformer in aligned.conformers_3d() {
            assert_eq!(conformer.coordinates(), &[[0., 0., 0.], [1., 0., 0.]]);
        }
        let mut inplace = source.clone();
        assert_eq!(
            inplace.align_conformers_with_params_(&params).unwrap(),
            report
        );
        assert_eq!(inplace.conformers_3d(), aligned.conformers_3d());
        let bad = ConformerAlignmentParameters {
            atom_indices: Some(vec![2]),
            ..Default::default()
        };
        assert!(matches!(
            source.with_aligned_conformers_with_params(&bad),
            Err(OperationError::Alignment(
                AlignmentError::ProbeAtomOutOfRange {
                    index: 2,
                    atom_count: 2
                }
            ))
        ));
    }
}
#[test]
fn one_conformer_does_not_validate_unused_weights() {
    let m = one();
    let p = ConformerAlignmentParameters {
        weights: Some(vec![1.]),
        ..Default::default()
    };
    let (after, report) = m.with_aligned_conformers_with_params(&p).unwrap();
    assert!(report.rmsds.is_empty());
    assert_eq!(after.conformers_3d(), m.conformers_3d());
    let p = ConformerAlignmentParameters {
        weights: Some(vec![1., 0.]),
        ..Default::default()
    };
    assert!(
        m.with_aligned_conformers_with_params(&p)
            .unwrap()
            .1
            .rmsds
            .is_empty()
    );
}
#[test]
fn zero_max_matches_keeps_all_three_source_paths() {
    let m = one();
    let best = BestAlignmentParameters {
        max_matches: 0,
        ..Default::default()
    };
    assert_eq!(m.best_rmsd_to_with_params(&m, &best).unwrap(), 0.);
    let coord = CoordinateRmsdParameters {
        max_matches: 0,
        ..Default::default()
    };
    assert_eq!(m.coordinate_rmsd_to_with_params(&m, &coord).unwrap(), 0.);
    let all = AllConformerRmsdParameters {
        max_matches: 0,
        ..Default::default()
    };
    assert!(
        m.all_conformer_best_rmsds_with_params(&all)
            .unwrap()
            .is_empty()
    );
}
#[test]
fn coordinate_negative_weight_is_not_alignment_positive_weight_policy() {
    let m = one();
    let p = CoordinateRmsdParameters {
        atom_maps: vec![vec![
            AlignmentAtomMap::new(0, 0),
            AlignmentAtomMap::new(1, 1),
        ]],
        weights: Some(vec![1., -1.]),
        ..Default::default()
    };
    assert_eq!(m.coordinate_rmsd_to_with_params(&m, &p).unwrap(), 0.);
}
#[test]
fn source_matching_error_precedes_missing_conformer_for_best_and_coordinate() {
    let c = Molecule::from_smiles("C").unwrap();
    let o = Molecule::from_smiles("O").unwrap();
    assert_eq!(
        c.best_alignment_to(&o),
        Err(AlignmentError::NoSubstructureMatch)
    );
    assert_eq!(
        c.coordinate_rmsd_to(&o),
        Err(AlignmentError::NoSubstructureMatch)
    );
    assert_eq!(
        c.alignment_transform_to(&o),
        Err(AlignmentError::NoConformers)
    );
}
