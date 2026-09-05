use cosmolkit_core::{
    AlignmentAtomMap, AlignmentError, AlignmentParameters, AllConformerRmsdParameters, BestAlignmentParameters,
    Conformer3D, ConformerAlignmentParameters, ConformerAlignmentReport, CoordinateRmsdParameters, Molecule,
    OperationError, mol_transforms,
};

fn molecule_with_conformer(smiles: &str, id: usize, coordinates: Vec<[f64; 3]>) -> Molecule {
    let molecule = Molecule::from_smiles(smiles).expect("molecule graph");
    let mut builder = molecule.to_builder();
    builder
        .add_conformer(Conformer3D::new(id, coordinates, true))
        .expect("3D conformer");
    builder.build().expect("molecule with conformer")
}

fn identity_map(count: usize) -> Vec<AlignmentAtomMap> {
    (0..count)
        .map(|index| AlignmentAtomMap {
            probe_atom: index,
            reference_atom: index,
        })
        .collect()
}

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1.0e-10,
        "actual={actual}, expected={expected}"
    );
}

fn assert_coordinates_close(actual: &[[f64; 3]], expected: &[[f64; 3]]) {
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(expected) {
        for axis in 0..3 {
            assert_close(actual[axis], expected[axis]);
        }
    }
}

#[test]
fn operation_output_signatures_preserve_existing_operations_and_type_alignment_reports() {
    let _: fn(&Molecule, Vec<[f64; 3]>, bool) -> Result<Molecule, OperationError> = Molecule::with_only_3d_conformer;
    let _: fn(&Molecule, ConformerAlignmentParameters) -> Result<(Molecule, ConformerAlignmentReport), OperationError> =
        Molecule::with_aligned_conformers_with_params;
    let _: fn(&mut Molecule, ConformerAlignmentParameters) -> Result<ConformerAlignmentReport, OperationError> =
        Molecule::align_conformers_with_params_;
    let _: fn(
        &Molecule,
        &Molecule,
        &AlignmentParameters,
    ) -> Result<(Molecule, cosmolkit_core::AlignmentResult), OperationError> = Molecule::with_alignment_to;
    let _: fn(
        &mut Molecule,
        &Molecule,
        &AlignmentParameters,
    ) -> Result<cosmolkit_core::AlignmentResult, OperationError> = Molecule::align_to_;
}

#[test]
fn explicit_alignment_mutation_is_registered_value_style_and_in_place() {
    let probe = molecule_with_conformer("CCC", 7, vec![[3.0, -2.0, 1.0], [4.0, -2.0, 1.0], [3.0, 0.0, 1.0]]);
    let reference = molecule_with_conformer("CCC", 17, vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]]);
    let params = AlignmentParameters {
        probe_conformer_id: 7,
        reference_conformer_id: 17,
        atom_map: Some(identity_map(3)),
        ..AlignmentParameters::default()
    };
    let source_coordinates = probe.conformers_3d()[0].coordinates().to_vec();
    let (aligned, result) = probe
        .with_alignment_to(&reference, &params)
        .expect("registered value-style alignment");
    assert_close(result.rmsd, 0.0);
    assert_eq!(probe.conformers_3d()[0].coordinates(), source_coordinates);
    assert_coordinates_close(
        aligned.conformers_3d()[0].coordinates(),
        reference.conformers_3d()[0].coordinates(),
    );

    let mut in_place = probe;
    let result = in_place
        .align_to_(&reference, &params)
        .expect("registered in-place alignment");
    assert_close(result.rmsd, 0.0);
    assert_coordinates_close(
        in_place.conformers_3d()[0].coordinates(),
        reference.conformers_3d()[0].coordinates(),
    );
}

#[test]
fn coordinate_value_transforms_leave_the_source_unchanged() {
    let molecule = Molecule::from_smiles("CC").expect("ethane");
    let with_coordinates = molecule
        .with_only_3d_conformer(vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], true)
        .expect("3D conformer");
    let original_coordinates = with_coordinates.conformers_3d()[0].coordinates().to_vec();

    let transformed = mol_transforms::transform_conformer(
        with_coordinates.clone(),
        &[
            [1.0, 0.0, 0.0, 2.0],
            [0.0, 1.0, 0.0, 3.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        0,
    )
    .expect("coordinate transform");

    assert_eq!(
        with_coordinates.conformers_3d()[0].coordinates(),
        original_coordinates.as_slice()
    );
    assert_eq!(
        transformed.conformers_3d()[0].coordinates(),
        &[[2.0, 3.0, 0.0], [3.5, 3.0, 0.0]]
    );
}

#[test]
fn conformer_identity_is_preserved_by_value_transforms() {
    let molecule = Molecule::from_smiles("C").expect("methane graph");
    let molecule = molecule
        .with_only_3d_conformer(vec![[0.0, 0.0, 0.0]], true)
        .expect("3D conformer");
    let conformer = molecule.conformers_3d()[0].clone();
    let named = Conformer3D::new(17, conformer.coordinates().to_vec(), true);
    assert_eq!(named.id(), 17);
    assert_eq!(named.coordinates(), conformer.coordinates());
}

#[test]
fn mol_transforms_resolve_sparse_stored_conformer_ids() {
    let molecule = Molecule::from_smiles("C").expect("methane graph");
    let mut builder = molecule.to_builder();
    builder
        .add_conformer(Conformer3D::new(17, vec![[1.0, 2.0, 3.0]], true))
        .expect("named conformer");
    let molecule = builder.build().expect("molecule");
    let transformed = mol_transforms::transform_conformer(
        molecule,
        &[
            [1.0, 0.0, 0.0, 4.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        17,
    )
    .expect("transform sparse conformer ID");
    assert_eq!(transformed.conformers_3d()[0].id(), 17);
    assert_eq!(transformed.conformers_3d()[0].coordinates(), &[[5.0, 2.0, 3.0]]);
}

#[test]
fn best_alignment_accepts_source_thread_count_semantics() {
    let molecule = Molecule::from_smiles("C").expect("methane graph");
    let molecule = molecule
        .with_only_3d_conformer(vec![[0.0, 0.0, 0.0]], true)
        .expect("3D conformer");
    let params = BestAlignmentParameters {
        num_threads: 2,
        atom_maps: vec![vec![cosmolkit_core::AlignmentAtomMap {
            probe_atom: 0,
            reference_atom: 0,
        }]],
        ..BestAlignmentParameters::default()
    };
    assert_eq!(
        molecule
            .best_alignment_to(&molecule, &params)
            .expect("parallel best alignment")
            .rmsd,
        0.0
    );
}

#[test]
fn every_negative_conformer_selector_uses_the_first_stored_conformer() {
    let molecule = Molecule::from_smiles("C").expect("methane graph");
    let mut builder = molecule.to_builder();
    builder
        .add_conformer(Conformer3D::new(7, vec![[1.0, 2.0, 3.0]], true))
        .expect("first conformer");
    builder
        .add_conformer(Conformer3D::new(17, vec![[9.0, 9.0, 9.0]], true))
        .expect("second conformer");
    let molecule = builder.build().expect("molecule");
    let params = AlignmentParameters {
        probe_conformer_id: -8,
        reference_conformer_id: -2,
        atom_map: Some(vec![cosmolkit_core::AlignmentAtomMap {
            probe_atom: 0,
            reference_atom: 0,
        }]),
        ..AlignmentParameters::default()
    };
    assert_eq!(
        molecule
            .alignment_transform_to(&molecule, &params)
            .expect("negative IDs select first")
            .rmsd,
        0.0
    );
}

#[test]
fn explicit_map_transform_matches_rdkit_orientation_weights_and_zero_iterations() {
    let probe_coordinates = vec![[3.0, -2.0, 1.0], [4.0, -2.0, 1.0], [3.0, 0.0, 1.0]];
    let reference_coordinates = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
    let probe = molecule_with_conformer("CCC", 7, probe_coordinates.clone());
    let reference = molecule_with_conformer("CCC", 17, reference_coordinates.clone());
    let params = AlignmentParameters {
        probe_conformer_id: 7,
        reference_conformer_id: 17,
        atom_map: Some(identity_map(3)),
        weights: Some(vec![1.0, 2.0, 3.0]),
        max_iterations: 0,
        ..AlignmentParameters::default()
    };

    let result = probe
        .alignment_transform_to(&reference, &params)
        .expect("explicit-map transform");
    assert_close(result.rmsd, 0.0);
    assert_eq!(result.atom_map, identity_map(3));
    let expected = [
        [1.0, 0.0, 0.0, -3.0],
        [0.0, 1.0, 0.0, 2.0],
        [0.0, 0.0, 1.0, -1.0],
        [0.0, 0.0, 0.0, 1.0],
    ];
    for (actual_row, expected_row) in result.transform.matrix.iter().zip(expected) {
        for (actual, expected) in actual_row.iter().zip(expected_row) {
            assert_close(*actual, expected);
        }
    }
    assert_eq!(probe.conformers_3d()[0].coordinates(), probe_coordinates);
    assert_eq!(reference.conformers_3d()[0].coordinates(), reference_coordinates);
}

#[test]
fn coordinate_rmsd_is_weighted_in_the_existing_frame_and_does_not_align() {
    let probe_coordinates = vec![[3.0, -2.0, 1.0], [4.0, -2.0, 1.0], [3.0, 0.0, 1.0]];
    let reference_coordinates = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
    let probe = molecule_with_conformer("CCC", 7, probe_coordinates.clone());
    let reference = molecule_with_conformer("CCC", 17, reference_coordinates.clone());
    let params = CoordinateRmsdParameters {
        probe_conformer_id: 7,
        reference_conformer_id: 17,
        atom_maps: vec![identity_map(3)],
        weights: Some(vec![1.0, 2.0, 3.0]),
        ..CoordinateRmsdParameters::default()
    };

    let rmsd = probe
        .coordinate_rmsd_to(&reference, &params)
        .expect("coordinate-frame RMSD");
    assert_close(rmsd, 28.0_f64.sqrt());
    assert_eq!(probe.conformers_3d()[0].coordinates(), probe_coordinates);
    assert_eq!(reference.conformers_3d()[0].coordinates(), reference_coordinates);
}

#[test]
fn explicit_map_reflection_matches_rdkit_source_regression() {
    let reference = molecule_with_conformer(
        "CCCC",
        0,
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    );
    let probe = molecule_with_conformer(
        "CCCC",
        0,
        vec![[2.0, 2.0, 3.0], [3.0, 2.0, 3.0], [2.0, 2.0, 4.0], [2.0, 3.0, 3.0]],
    );
    let without_reflection = probe
        .alignment_transform_to(
            &reference,
            &AlignmentParameters {
                atom_map: Some(identity_map(4)),
                ..AlignmentParameters::default()
            },
        )
        .expect("proper-rotation alignment");
    assert_close(without_reflection.rmsd, 0.5);

    let with_reflection = probe
        .alignment_transform_to(
            &reference,
            &AlignmentParameters {
                atom_map: Some(identity_map(4)),
                reflect: true,
                ..AlignmentParameters::default()
            },
        )
        .expect("reflected alignment");
    assert_close(with_reflection.rmsd, 0.0);
}

#[test]
fn explicit_map_validation_reports_each_source_boundary_without_mutation() {
    let molecule = molecule_with_conformer("CC", 17, vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
    let original = molecule.conformers_3d()[0].coordinates().to_vec();

    let cases = [
        (Vec::new(), None, AlignmentError::EmptyAtomMap),
        (
            vec![AlignmentAtomMap {
                probe_atom: 2,
                reference_atom: 0,
            }],
            None,
            AlignmentError::ProbeAtomOutOfRange {
                index: 2,
                atom_count: 2,
            },
        ),
        (
            vec![AlignmentAtomMap {
                probe_atom: 0,
                reference_atom: 2,
            }],
            None,
            AlignmentError::ReferenceAtomOutOfRange {
                index: 2,
                atom_count: 2,
            },
        ),
        (
            identity_map(2),
            Some(vec![1.0]),
            AlignmentError::WeightCountMismatch {
                map_len: 2,
                weight_len: 1,
            },
        ),
        (
            identity_map(2),
            Some(vec![1.0, 0.0]),
            AlignmentError::NonPositiveWeight { index: 1 },
        ),
        (
            identity_map(2),
            Some(vec![1.0, f64::NAN]),
            AlignmentError::NonPositiveWeight { index: 1 },
        ),
    ];
    for (atom_map, weights, expected) in cases {
        let alignment_params = AlignmentParameters {
            probe_conformer_id: 17,
            reference_conformer_id: 17,
            atom_map: Some(atom_map),
            weights: weights.clone(),
            ..AlignmentParameters::default()
        };
        assert_eq!(
            molecule.alignment_transform_to(&molecule, &alignment_params),
            Err(expected.clone())
        );
        let rmsd_params = CoordinateRmsdParameters {
            probe_conformer_id: 17,
            reference_conformer_id: 17,
            atom_maps: vec![alignment_params.atom_map.expect("explicit map")],
            weights,
            ..CoordinateRmsdParameters::default()
        };
        assert_eq!(molecule.coordinate_rmsd_to(&molecule, &rmsd_params), Err(expected));
    }
    assert_eq!(molecule.conformers_3d()[0].coordinates(), original);
}

#[test]
fn automatic_best_alignment_uses_all_non_unique_matches() {
    let probe_coordinates = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
    let reference_coordinates = vec![[0.0, 2.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
    let probe = molecule_with_conformer("CCC", 7, probe_coordinates.clone());
    let reference = molecule_with_conformer("CCC", 17, reference_coordinates.clone());
    let result = probe
        .best_alignment_to(
            &reference,
            &BestAlignmentParameters {
                probe_conformer_id: 7,
                reference_conformer_id: 17,
                ignore_hydrogens: false,
                ..BestAlignmentParameters::default()
            },
        )
        .expect("automatic best alignment");

    assert_close(result.rmsd, 0.0);
    assert_eq!(result.atom_map.len(), 3);
    assert_eq!(probe.conformers_3d()[0].coordinates(), probe_coordinates);
    assert_eq!(reference.conformers_3d()[0].coordinates(), reference_coordinates);
}

#[test]
fn best_alignment_default_ignores_explicit_hydrogens_only_for_automatic_maps() {
    let molecule = molecule_with_conformer("C", 0, vec![[0.0, 0.0, 0.0]])
        .with_hydrogens()
        .expect("explicit hydrogens");
    let automatic = molecule
        .best_alignment_to(&molecule, &BestAlignmentParameters::default())
        .expect("automatic map");
    assert_eq!(automatic.atom_map.len(), 1);

    let explicit = molecule
        .best_alignment_to(
            &molecule,
            &BestAlignmentParameters {
                atom_maps: vec![identity_map(5)],
                ..BestAlignmentParameters::default()
            },
        )
        .expect("explicit map");
    assert_eq!(explicit.atom_map.len(), 5);
}

#[test]
fn best_alignment_preserves_source_strict_less_tie_selection() {
    let coordinates = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let molecule = molecule_with_conformer("CCC", 0, coordinates);
    let bad = vec![
        AlignmentAtomMap {
            probe_atom: 0,
            reference_atom: 0,
        },
        AlignmentAtomMap {
            probe_atom: 1,
            reference_atom: 2,
        },
        AlignmentAtomMap {
            probe_atom: 2,
            reference_atom: 1,
        },
    ];
    let forward = identity_map(3);
    let reverse = vec![
        AlignmentAtomMap {
            probe_atom: 0,
            reference_atom: 2,
        },
        AlignmentAtomMap {
            probe_atom: 1,
            reference_atom: 1,
        },
        AlignmentAtomMap {
            probe_atom: 2,
            reference_atom: 0,
        },
    ];
    let maps = vec![bad, forward.clone(), reverse.clone()];

    let serial = molecule
        .best_alignment_to(
            &molecule,
            &BestAlignmentParameters {
                atom_maps: maps.clone(),
                num_threads: 1,
                ..BestAlignmentParameters::default()
            },
        )
        .expect("serial best alignment");
    assert_eq!(serial.atom_map, forward);

    let parallel = molecule
        .best_alignment_to(
            &molecule,
            &BestAlignmentParameters {
                atom_maps: maps,
                num_threads: 2,
                ..BestAlignmentParameters::default()
            },
        )
        .expect("parallel best alignment");
    assert_eq!(parallel.atom_map, reverse);
    assert_close(serial.rmsd, parallel.rmsd);
}

#[test]
fn automatic_map_reports_no_match_and_zero_max_matches_like_source() {
    let carbon = molecule_with_conformer("C", 0, vec![[0.0, 0.0, 0.0]]);
    let oxygen = molecule_with_conformer("O", 0, vec![[0.0, 0.0, 0.0]]);
    assert_eq!(
        carbon.best_alignment_to(&oxygen, &BestAlignmentParameters::default()),
        Err(AlignmentError::NoSubstructureMatch)
    );
    assert_eq!(
        carbon.best_alignment_to(
            &carbon,
            &BestAlignmentParameters {
                max_matches: 0,
                ..BestAlignmentParameters::default()
            }
        ),
        Err(AlignmentError::NoSubstructureMatch)
    );
}

#[test]
fn coordinate_rmsd_chooses_best_candidate_without_alignment() {
    let probe = molecule_with_conformer("CCC", 0, vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    let reference = molecule_with_conformer("CCC", 0, vec![[2.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]);
    let rmsd = probe
        .coordinate_rmsd_to(
            &reference,
            &CoordinateRmsdParameters {
                atom_maps: vec![
                    identity_map(3),
                    vec![
                        AlignmentAtomMap {
                            probe_atom: 0,
                            reference_atom: 2,
                        },
                        AlignmentAtomMap {
                            probe_atom: 1,
                            reference_atom: 1,
                        },
                        AlignmentAtomMap {
                            probe_atom: 2,
                            reference_atom: 0,
                        },
                    ],
                ],
                ..CoordinateRmsdParameters::default()
            },
        )
        .expect("best coordinate-frame RMSD");
    assert_close(rmsd, 0.0);
}

#[test]
fn automatic_best_alignment_symmetrizes_conjugated_terminal_atoms_on_a_clone() {
    let probe_coordinates = vec![[-2.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, -0.3, 0.0]];
    let reference_coordinates = vec![[-2.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, -0.3, 0.0], [0.0, 1.0, 0.0]];
    let probe = molecule_with_conformer("CC(=O)[O-]", 0, probe_coordinates.clone());
    let reference = molecule_with_conformer("CC(=O)[O-]", 0, reference_coordinates.clone());

    let symmetrized = probe
        .best_alignment_to(&reference, &BestAlignmentParameters::default())
        .expect("symmetrized terminal-group alignment");
    assert_close(symmetrized.rmsd, 0.0);
    assert_eq!(
        symmetrized.atom_map,
        vec![
            AlignmentAtomMap {
                probe_atom: 0,
                reference_atom: 0
            },
            AlignmentAtomMap {
                probe_atom: 1,
                reference_atom: 1
            },
            AlignmentAtomMap {
                probe_atom: 2,
                reference_atom: 3
            },
            AlignmentAtomMap {
                probe_atom: 3,
                reference_atom: 2
            },
        ]
    );

    let unsymmetrized = probe
        .best_alignment_to(
            &reference,
            &BestAlignmentParameters {
                symmetrize_conjugated_terminal_groups: false,
                ..BestAlignmentParameters::default()
            },
        )
        .expect("ordinary topology alignment");
    assert!(unsymmetrized.rmsd > 1.0e-6);
    assert_eq!(probe.conformers_3d()[0].coordinates(), probe_coordinates);
    assert_eq!(reference.conformers_3d()[0].coordinates(), reference_coordinates);
}

#[test]
fn all_conformer_best_rmsds_preserve_source_triangular_order_and_ids() {
    let molecule = Molecule::from_smiles("CCC").expect("molecule graph");
    let mut builder = molecule.to_builder();
    builder
        .add_conformer(Conformer3D::new(
            7,
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
            true,
        ))
        .expect("first conformer");
    builder
        .add_conformer(Conformer3D::new(
            17,
            vec![[3.0, 4.0, 0.0], [4.0, 4.0, 0.0], [3.0, 6.0, 0.0]],
            true,
        ))
        .expect("second conformer");
    builder
        .add_conformer(Conformer3D::new(
            29,
            vec![[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [0.0, 2.0, 0.0]],
            true,
        ))
        .expect("third conformer");
    let molecule = builder.build().expect("multi-conformer molecule");
    let original: Vec<_> = molecule
        .conformers_3d()
        .iter()
        .map(|conformer| conformer.coordinates().to_vec())
        .collect();
    let params = AllConformerRmsdParameters {
        atom_maps: vec![identity_map(3)],
        num_threads: 1,
        ..AllConformerRmsdParameters::default()
    };
    let serial = molecule
        .all_conformer_best_rmsds(&params)
        .expect("serial all-conformer RMSDs");
    assert_eq!(
        serial
            .iter()
            .map(|entry| (entry.probe_conformer_id, entry.reference_conformer_id))
            .collect::<Vec<_>>(),
        vec![(17, 7), (29, 7), (29, 17)]
    );
    assert_close(serial[0].rmsd, 0.0);

    let parallel = molecule
        .all_conformer_best_rmsds(&AllConformerRmsdParameters {
            num_threads: 2,
            ..params
        })
        .expect("parallel all-conformer RMSDs");
    assert_eq!(parallel, serial);
    assert_eq!(
        molecule
            .conformers_3d()
            .iter()
            .map(|conformer| conformer.coordinates().to_vec())
            .collect::<Vec<_>>(),
        original
    );
}

#[test]
fn conformer_alignment_mutation_is_explicit_and_value_style_preserves_source() {
    let molecule = Molecule::from_smiles("CCC").expect("molecule graph");
    let mut builder = molecule.to_builder();
    builder
        .add_conformer(Conformer3D::new(
            7,
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
            true,
        ))
        .expect("reference conformer");
    builder
        .add_conformer(Conformer3D::new(
            17,
            vec![[3.0, 4.0, 5.0], [4.0, 4.0, 5.0], [3.0, 6.0, 5.0]],
            true,
        ))
        .expect("probe conformer");
    let molecule = builder.build().expect("multi-conformer molecule");
    let original = molecule.conformers_3d()[1].coordinates().to_vec();

    let (aligned, report) = molecule
        .with_aligned_conformers_with_params(ConformerAlignmentParameters {
            atom_indices: Some(vec![0, 1, 2]),
            conformer_ids: Some(vec![7, 17]),
            ..ConformerAlignmentParameters::default()
        })
        .expect("value-style conformer alignment");
    assert_eq!(report.rmsds.len(), 1);
    assert_close(report.rmsds[0], 0.0);
    assert_eq!(molecule.conformers_3d()[1].coordinates(), original);
    for (actual, expected) in aligned.conformers_3d()[1]
        .coordinates()
        .iter()
        .zip(aligned.conformers_3d()[0].coordinates())
    {
        for axis in 0..3 {
            assert_close(actual[axis], expected[axis]);
        }
    }

    let mut inplace = molecule.clone();
    let report = inplace
        .align_conformers_with_params_(ConformerAlignmentParameters::default())
        .expect("explicit in-place conformer alignment");
    assert_eq!(report.rmsds.len(), 1);
    assert_close(report.rmsds[0], 0.0);
    assert_ne!(inplace.conformers_3d()[1].coordinates(), original);
}

#[test]
fn conformer_alignment_report_follows_selected_source_order() {
    let molecule = Molecule::from_smiles("CCC").expect("molecule graph");
    let mut builder = molecule.to_builder();
    for (id, coordinates) in [
        (7, vec![[0.0, 0.0, 0.0], [1.4, 0.0, 0.0], [0.2, 1.7, 0.0]]),
        (17, vec![[2.0, 3.0, 1.0], [3.8, 3.1, 1.0], [2.1, 4.2, 1.4]]),
        (29, vec![[-1.0, 2.0, 0.5], [0.1, 2.2, 0.6], [-0.8, 4.1, 0.2]]),
        (41, vec![[9.0, 8.0, 7.0], [8.0, 7.0, 6.0], [7.0, 6.0, 5.0]]),
    ] {
        builder
            .add_conformer(Conformer3D::new(id, coordinates, true))
            .expect("named conformer");
    }
    let molecule = builder.build().expect("multi-conformer molecule");
    let untouched = molecule.conformers_3d()[3].coordinates().to_vec();
    let atom_map = identity_map(3);
    let expected = [7, 17].map(|probe_conformer_id| {
        molecule
            .alignment_transform_to(
                &molecule,
                &AlignmentParameters {
                    probe_conformer_id,
                    reference_conformer_id: 29,
                    atom_map: Some(atom_map.clone()),
                    ..AlignmentParameters::default()
                },
            )
            .expect("read-only selected-conformer alignment")
            .rmsd
    });

    let (aligned, report) = molecule
        .with_aligned_conformers_with_params(ConformerAlignmentParameters {
            atom_indices: Some(vec![0, 1, 2]),
            conformer_ids: Some(vec![29, 7, 17]),
            ..ConformerAlignmentParameters::default()
        })
        .expect("selected conformer alignment");

    assert_eq!(report.rmsds.len(), expected.len());
    for (actual, expected) in report.rmsds.iter().zip(expected) {
        assert_close(*actual, expected);
    }
    assert_eq!(aligned.conformers_3d()[3].coordinates(), untouched);
    assert_eq!(molecule.conformers_3d()[3].coordinates(), untouched);
}

#[test]
fn empty_conformer_selection_uses_all_conformers_like_the_source_wrapper() {
    let molecule = Molecule::from_smiles("CCC").expect("molecule graph");
    let mut builder = molecule.to_builder();
    builder
        .add_conformer(Conformer3D::new(
            7,
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            true,
        ))
        .expect("reference conformer");
    builder
        .add_conformer(Conformer3D::new(
            17,
            vec![[2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [4.0, 0.0, 0.0]],
            true,
        ))
        .expect("probe conformer");
    let molecule = builder.build().expect("multi-conformer molecule");
    let source_probe = molecule.conformers_3d()[1].coordinates().to_vec();
    let (aligned, report) = molecule
        .with_aligned_conformers_with_params(ConformerAlignmentParameters {
            conformer_ids: Some(Vec::new()),
            ..ConformerAlignmentParameters::default()
        })
        .expect("empty conformer selection");

    assert_eq!(report.rmsds, vec![0.0]);
    assert_eq!(
        aligned.conformers_3d()[1].coordinates(),
        aligned.conformers_3d()[0].coordinates()
    );
    assert_eq!(molecule.conformers_3d()[1].coordinates(), source_probe);
}

#[test]
fn conformer_alignment_errors_return_no_report_and_keep_complete_state() {
    let molecule = molecule_with_conformer("CCC", 7, vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]]);
    let before = molecule.conformers_3d()[0].coordinates().to_vec();
    let params = ConformerAlignmentParameters {
        conformer_ids: Some(vec![999]),
        ..ConformerAlignmentParameters::default()
    };

    let error = molecule
        .with_aligned_conformers_with_params(params.clone())
        .expect_err("missing conformer must return a structured operation error");
    assert!(matches!(
        error,
        OperationError::Alignment {
            source: AlignmentError::ConformerNotFound { id: 999 },
            ..
        }
    ));
    assert_eq!(molecule.conformers_3d()[0].coordinates(), before);

    let mut inplace = molecule;
    let error = inplace
        .align_conformers_with_params_(params)
        .expect_err("in-place missing conformer must not return a report");
    assert!(matches!(
        error,
        OperationError::Alignment {
            source: AlignmentError::ConformerNotFound { id: 999 },
            ..
        }
    ));
    assert_eq!(inplace.num_atoms(), 3);
    assert_eq!(inplace.conformers_3d()[0].coordinates(), before);
}
