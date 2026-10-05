#![cfg(feature = "cap-alignment")]
use cosmolkit::*;
#[test]
fn six_original_constructors_retain_all_fields_and_defaults() {
    let map = AlignmentAtomMap::new(2, 7);
    assert_eq!((map.probe_atom, map.reference_atom), (2, 7));
    assert_eq!(
        AlignmentParameters::new(-1, -1, None, None, false, 50),
        AlignmentParameters::default()
    );
    assert_eq!(
        BestAlignmentParameters::new(-1, -1, vec![], None, false, 50, 1000000, true, true, 1),
        BestAlignmentParameters::default()
    );
    assert_eq!(
        CoordinateRmsdParameters::new(-1, -1, vec![], None, 1000000, true),
        CoordinateRmsdParameters::default()
    );
    assert_eq!(
        AllConformerRmsdParameters::new(vec![], None, 1000000, true, true, 1),
        AllConformerRmsdParameters::default()
    );
    assert_eq!(
        ConformerAlignmentParameters::new(None, None, None, false, 50),
        ConformerAlignmentParameters::default()
    );
    let p = BestAlignmentParameters::new(
        17,
        7,
        vec![vec![map]],
        Some(vec![2.]),
        true,
        0,
        0,
        false,
        false,
        -1,
    );
    assert_eq!(
        (
            p.probe_conformer_id,
            p.reference_conformer_id,
            p.max_iterations,
            p.max_matches,
            p.num_threads
        ),
        (17, 7, 0, 0, -1)
    );
    assert!(p.reflect);
    assert!(!p.ignore_hydrogens && !p.symmetrize_conjugated_terminal_groups);
    assert_eq!(p.weights, Some(vec![2.]));
    assert_eq!(p.atom_maps, vec![vec![map]]);
    let p = ConformerAlignmentParameters::new(
        Some(vec![2]),
        Some(vec![17, 7]),
        Some(vec![3.]),
        true,
        0,
    );
    assert_eq!(p.atom_indices, Some(vec![2]));
    assert_eq!(p.conformer_ids, Some(vec![17, 7]));
    assert_eq!(p.weights, Some(vec![3.]));
    assert!(p.reflect);
    assert_eq!(p.max_iterations, 0);
}
#[test]
fn public_value_getters_and_canonical_registry_are_consistent() {
    let t = AlignmentTransform {
        matrix: [
            [1., 0., 0., 0.],
            [0., 1., 0., 0.],
            [0., 0., 1., 0.],
            [0., 0., 0., 1.],
        ],
    };
    let r = AlignmentResult {
        rmsd: 2.,
        transform: t,
        atom_map: vec![AlignmentAtomMap::new(2, 7)],
    };
    assert_eq!(r.rmsd(), 2.);
    assert_eq!(r.transform().matrix(), &t.matrix);
    assert_eq!(r.atom_map(), r.atom_map.as_slice());
    let p = ConformerRmsd {
        probe_conformer_id: 17,
        reference_conformer_id: 7,
        rmsd: 3.,
    };
    assert_eq!(
        (p.probe_conformer_id(), p.reference_conformer_id(), p.rmsd()),
        (17, 7, 3.)
    );
    let report = ConformerAlignmentReport {
        rmsds: vec![3., 4.],
    };
    assert_eq!(report.rmsds(), &[3., 4.]);
    for typ in [
        "AlignmentAtomMap",
        "AlignmentParameters",
        "BestAlignmentParameters",
        "CoordinateRmsdParameters",
        "AllConformerRmsdParameters",
        "ConformerAlignmentParameters",
        "AlignmentResult",
        "AlignmentTransform",
        "ConformerRmsd",
        "ConformerAlignmentReport",
        "AlignmentError",
    ] {
        let row = BINDING_CONTRACT
            .iter()
            .find(|r| r.semantic_id == format!("types.{typ}"))
            .unwrap();
        assert_eq!(row.python_name, typ);
        assert_eq!(row.feature, "cap-alignment");
    }
    for typ in [
        "AlignmentAtomMap",
        "AlignmentParameters",
        "BestAlignmentParameters",
        "CoordinateRmsdParameters",
        "AllConformerRmsdParameters",
        "ConformerAlignmentParameters",
    ] {
        let row = BINDING_CONTRACT
            .iter()
            .find(|r| r.semantic_id == format!("{typ}.new"))
            .unwrap();
        assert!(row.callable.is_some());
    }
}
