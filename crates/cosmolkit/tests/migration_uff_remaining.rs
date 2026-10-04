use cosmolkit::{
    AromaticityError, Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondQueryPredicate,
    BondSpec, Element, Hybridization, HydrogenError, KekulizeError, Molecule, MoleculeProperties,
    OperationError, QueryNode, SanitizeError, SanitizeStage, SdfError, SdfGraph, SdfRecord,
    SmilesError, SmilesParseParams, UffConformerOptimizationParams, UffOptimizationErrorKind,
    ValenceError, ValenceModel, ValenceParams,
};
use cosmolkit_model::{Conformer2D, Conformer3D, CoordinateBlock};

#[test]
fn uff_remaining_parse_rejections_() {
    let cases = [
        (
            "line:112",
            "C[C@@]([C@@]12C)(C[C@H](O[C@@](O[C@H](CO)[C@@H](O)[C@@H]3O)([H])[C@@H]3O[C@@](O[C@@H](C)[C@H](O)[C@H]4O)([H])[C@@H]4O)[C@@]5([H])C6(C)C)[C@@](C[C@@H](O)[C@]1([H])[C@]([C@]7(O[C@@H](C(C)(O)C)CC7)C)([H])CC2)([H])[C@]5(CC[C@@H]6O)C)",
            "parse",
            "SMILES parsing failed: invalid SMILES syntax at byte 228: unmatched branch close",
        ),
        (
            "line:128",
            "Nc1ncnc2n(cnc12)[C@@H]3O[C@H](CN=[N]=N)[C@@H](O)[C@H]3O",
            "hydrogen",
            "SMILES hydrogen removal failed: hydrogen-removal sanitize failed: sanitization failed at Properties: property-cache valence assignment failed: Explicit valence for atom # 15 N, 4, is greater than permitted",
        ),
        (
            "line:130",
            "c1c(ccc2NC(CN=c(c21)(C)C)=O)O",
            "hydrogen",
            "SMILES hydrogen removal failed: hydrogen-removal sanitize failed: sanitization failed at Properties: property-cache valence assignment failed: Explicit valence for atom # 9 C, 5, is greater than permitted",
        ),
        (
            "line:135",
            "C=CC(C)=Cc1cc2[nH]c3c(ccc1C)c(C)c(CCC(=O)O)c-3=c-c1c(C)c3[nH]/c(c(C)c3-c=1CCC(=O)O)=c/c1[nH]/c(=c/2C)-c=1C",
            "hydrogen",
            "SMILES hydrogen removal failed: hydrogen-removal sanitize failed: sanitization failed at Kekulize: could not kekulize molecule; remaining atoms: [AtomId(5), AtomId(6), AtomId(7), AtomId(9), AtomId(10), AtomId(11), AtomId(12), AtomId(13), AtomId(15), AtomId(17), AtomId(26), AtomId(28), AtomId(31), AtomId(33), AtomId(40)]",
        ),
        (
            "line:136",
            "[C+]([F])([F])([F])([F])[F]",
            "hydrogen",
            "SMILES hydrogen removal failed: hydrogen-removal sanitize failed: sanitization failed at Properties: property-cache valence assignment failed: Explicit valence for atom # 0 C, 5, is greater than permitted",
        ),
        (
            "line:138",
            "[F][He+]([F])[F]",
            "hydrogen",
            "SMILES hydrogen removal failed: hydrogen-removal sanitize failed: sanitization failed at Properties: property-cache valence assignment failed: Explicit valence for atom # 1 He, 3, is greater than permitted",
        ),
        (
            "line:139",
            "Br[Br-]Br.C=1C=CC(=CC1)[N+](C)C",
            "hydrogen",
            "SMILES hydrogen removal failed: hydrogen-removal sanitize failed: sanitization failed at Properties: property-cache valence assignment failed: Explicit valence for atom # 1 Br, 2, is greater than permitted",
        ),
        (
            "line:140",
            "FCl(F)F",
            "hydrogen",
            "SMILES hydrogen removal failed: hydrogen-removal sanitize failed: sanitization failed at Properties: property-cache valence assignment failed: Explicit valence for atom # 1 Cl, 3, is greater than permitted",
        ),
        (
            "line:141",
            "C1CCn2ccc(NCCOCCNc3ccn(CC1)c1ccccc31)c1ccccc21",
            "hydrogen",
            "SMILES hydrogen removal failed: hydrogen-removal sanitize failed: sanitization failed at Kekulize: could not kekulize molecule; remaining atoms: [AtomId(4), AtomId(5), AtomId(6), AtomId(14), AtomId(15), AtomId(16), AtomId(20), AtomId(21), AtomId(22), AtomId(23), AtomId(24), AtomId(25), AtomId(26), AtomId(27), AtomId(28), AtomId(29), AtomId(30), AtomId(31)]",
        ),
        (
            "line:144",
            "[n+:26]1(-[c:27]2[c:30]([CH3:36])[cH:35][cH:41][cH:37][c:31]2[CH3:38])[al+:28][n+:32](-[c:39]2[c:42]([CH3:46])[cH:45][cH:49][cH:47][c:43]2[CH3:48])[c:40]([CH3:44])[cH-:33][c:29]1[CH3:34]",
            "parse",
            "SMILES parsing failed: invalid atom at byte 70: al+:28",
        ),
        (
            "line:149",
            "CC1=CSC2=C1CCN(C(=O)[C@@h]1CCN(C)C(=O)C@@HNC3=CC(=CC(N4C=CC=C4)=C3)C(=O)N1C)C2",
            "parse",
            "SMILES parsing failed: invalid atom at byte 20: C@@h",
        ),
        (
            "line:150",
            "C=CC1=CC=C(NC(=O)[C@@h]2CCN(C)C(=O)C@@HNC3=CC(=CC(OCC(F)(F)F)=C3)C(=O)N2C)C=C1",
            "parse",
            "SMILES parsing failed: invalid atom at byte 17: C@@h",
        ),
    ];

    assert_eq!(cases.len(), 12);
    for (id, smiles, expected_kind, expected_display) in cases {
        let error = match Molecule::from_smiles(smiles) {
            Ok(molecule) => panic!(
                "CK public constructor accepted source-rejected {id}: atoms={}, bonds={}",
                molecule.num_atoms(),
                molecule.num_bonds()
            ),
            Err(error) => error,
        };
        let actual_kind = match &error {
            SmilesError::Parse(_) => "parse",
            SmilesError::Hydrogen(_) => "hydrogen",
            SmilesError::Sanitize(_) => "sanitize",
            SmilesError::Stereo(_) => "stereo",
            SmilesError::Construction(_) => "construction",
        };
        assert_eq!(actual_kind, expected_kind, "source row {id}");
        assert_eq!(error.to_string(), expected_display, "source row {id}");
    }
}

#[test]
fn uff_remaining_source_error_matrix_() {
    #[derive(Debug, PartialEq, Eq)]
    enum Outcome {
        Accepted {
            atoms: usize,
            bonds: usize,
        },
        Parse,
        Properties {
            hydrogen_wrapper: bool,
            stage: SanitizeStage,
            atom: usize,
            atomic_number: u8,
            calculated: Option<i32>,
        },
        AromaticityValence {
            hydrogen_wrapper: bool,
            stage: SanitizeStage,
            atom: usize,
            atomic_number: u8,
            calculated: Option<i32>,
        },
        Kekulize {
            hydrogen_wrapper: bool,
            stage: SanitizeStage,
            problem_atoms: Vec<usize>,
        },
        Other(String),
    }

    fn classify_error(error: &SmilesError, remove_hydrogens: bool) -> Outcome {
        let sanitize_error = match (remove_hydrogens, error) {
            (_, SmilesError::Parse(_)) => return Outcome::Parse,
            (true, SmilesError::Hydrogen(HydrogenError::Sanitize(error))) => error,
            (false, SmilesError::Sanitize(error)) => error,
            _ => return Outcome::Other(error.to_string()),
        };

        match sanitize_error {
            SanitizeError::Properties { stage, .. } => {
                let mut source = std::error::Error::source(sanitize_error);
                let mut valence = None;
                while let Some(error) = source {
                    if let Some(ValenceError::InvalidValence {
                        atom,
                        atomic_number,
                        calculated,
                        ..
                    }) = error.downcast_ref::<ValenceError>()
                    {
                        valence = Some((atom.index(), *atomic_number, *calculated));
                        break;
                    }
                    source = std::error::Error::source(error);
                }
                match valence {
                    Some((atom, atomic_number, calculated)) => Outcome::Properties {
                        hydrogen_wrapper: remove_hydrogens,
                        stage: *stage,
                        atom,
                        atomic_number,
                        calculated,
                    },
                    None => Outcome::Other(error_chain_text(sanitize_error)),
                }
            }
            SanitizeError::Kekulize { stage, source } => match source {
                KekulizeError::NotKekulizable { problem_atoms } => Outcome::Kekulize {
                    hydrogen_wrapper: remove_hydrogens,
                    stage: *stage,
                    problem_atoms: problem_atoms.iter().map(|atom| atom.index()).collect(),
                },
                _ => Outcome::Other(error_chain_text(sanitize_error)),
            },
            SanitizeError::Aromaticity {
                stage,
                source:
                    AromaticityError::Valence(ValenceError::InvalidValence {
                        atom,
                        atomic_number,
                        calculated,
                        ..
                    }),
            } => Outcome::AromaticityValence {
                hydrogen_wrapper: remove_hydrogens,
                stage: *stage,
                atom: atom.index(),
                atomic_number: *atomic_number,
                calculated: *calculated,
            },
            _ => Outcome::Other(error_chain_text(sanitize_error)),
        }
    }

    fn error_chain_text(error: &(dyn std::error::Error + 'static)) -> String {
        let mut text = error.to_string();
        let mut source = std::error::Error::source(error);
        while let Some(error) = source {
            text.push_str(" -> ");
            text.push_str(&error.to_string());
            source = std::error::Error::source(error);
        }
        text
    }

    let cases = [
        (
            "line:112",
            "C[C@@]([C@@]12C)(C[C@H](O[C@@](O[C@H](CO)[C@@H](O)[C@@H]3O)([H])[C@@H]3O[C@@](O[C@@H](C)[C@H](O)[C@H]4O)([H])[C@@H]4O)[C@@]5([H])C6(C)C)[C@@](C[C@@H](O)[C@]1([H])[C@]([C@]7(O[C@@H](C(C)(O)C)CC7)C)([H])CC2)([H])[C@]5(CC[C@@H]6O)C)",
            None,
        ),
        (
            "line:128",
            "Nc1ncnc2n(cnc12)[C@@H]3O[C@H](CN=[N]=N)[C@@H](O)[C@H]3O",
            Some((21, 23)),
        ),
        ("line:130", "c1c(ccc2NC(CN=c(c21)(C)C)=O)O", Some((15, 16))),
        (
            "line:135",
            "C=CC(C)=Cc1cc2[nH]c3c(ccc1C)c(C)c(CCC(=O)O)c-3=c-c1c(C)c3[nH]/c(c(C)c3-c=1CCC(=O)O)=c/c1[nH]/c(=c/2C)-c=1C",
            Some((48, 53)),
        ),
        ("line:136", "[C+]([F])([F])([F])([F])[F]", Some((6, 5))),
        ("line:138", "[F][He+]([F])[F]", Some((4, 3))),
        (
            "line:139",
            "Br[Br-]Br.C=1C=CC(=CC1)[N+](C)C",
            Some((12, 11)),
        ),
        ("line:140", "FCl(F)F", Some((4, 3))),
        (
            "line:141",
            "C1CCn2ccc(NCCOCCNc3ccn(CC1)c1ccccc31)c1ccccc21",
            Some((32, 36)),
        ),
        (
            "line:144",
            "[n+:26]1(-[c:27]2[c:30]([CH3:36])[cH:35][cH:41][cH:37][c:31]2[CH3:38])[al+:28][n+:32](-[c:39]2[c:42]([CH3:46])[cH:45][cH:49][cH:47][c:43]2[CH3:48])[c:40]([CH3:44])[cH-:33][c:29]1[CH3:34]",
            None,
        ),
        (
            "line:149",
            "CC1=CSC2=C1CCN(C(=O)[C@@h]1CCN(C)C(=O)C@@HNC3=CC(=CC(N4C=CC=C4)=C3)C(=O)N1C)C2",
            None,
        ),
        (
            "line:150",
            "C=CC1=CC=C(NC(=O)[C@@h]2CCN(C)C(=O)C@@HNC3=CC(=CC(OCC(F)(F)F)=C3)C(=O)N2C)C=C1",
            None,
        ),
    ];

    let mut observations = Vec::with_capacity(48);
    let mut mismatches = Vec::new();
    for (id, smiles, counts) in cases {
        for sanitize in [false, true] {
            for remove_hydrogens in [false, true] {
                let params = SmilesParseParams {
                    sanitize,
                    remove_hydrogens,
                    ..SmilesParseParams::default()
                };
                let actual = match Molecule::from_smiles_with_params(smiles, &params) {
                    Ok(molecule) => Outcome::Accepted {
                        atoms: molecule.num_atoms(),
                        bonds: molecule.num_bonds(),
                    },
                    Err(error) => classify_error(&error, remove_hydrogens),
                };
                let expected = if matches!(id, "line:112" | "line:144" | "line:149" | "line:150") {
                    Outcome::Parse
                } else if !sanitize {
                    let (atoms, bonds) = counts.expect("accepted source row has frozen counts");
                    Outcome::Accepted { atoms, bonds }
                } else {
                    let hydrogen_wrapper = remove_hydrogens;
                    match id {
                        "line:128" => Outcome::Properties {
                            hydrogen_wrapper,
                            stage: SanitizeStage::Properties,
                            atom: 15,
                            atomic_number: 7,
                            calculated: Some(4),
                        },
                        "line:130" => Outcome::Properties {
                            hydrogen_wrapper,
                            stage: SanitizeStage::Properties,
                            atom: 9,
                            atomic_number: 6,
                            calculated: Some(5),
                        },
                        "line:136" => Outcome::Properties {
                            hydrogen_wrapper,
                            stage: SanitizeStage::Properties,
                            atom: 0,
                            atomic_number: 6,
                            calculated: Some(5),
                        },
                        "line:138" => Outcome::Properties {
                            hydrogen_wrapper,
                            stage: SanitizeStage::Properties,
                            atom: 1,
                            atomic_number: 2,
                            calculated: Some(3),
                        },
                        "line:139" => Outcome::Properties {
                            hydrogen_wrapper,
                            stage: SanitizeStage::Properties,
                            atom: 1,
                            atomic_number: 35,
                            calculated: Some(2),
                        },
                        "line:140" => Outcome::Properties {
                            hydrogen_wrapper,
                            stage: SanitizeStage::Properties,
                            atom: 1,
                            atomic_number: 17,
                            calculated: Some(3),
                        },
                        "line:135" => Outcome::Kekulize {
                            hydrogen_wrapper,
                            stage: SanitizeStage::Kekulize,
                            problem_atoms: vec![
                                5, 6, 7, 9, 10, 11, 12, 13, 15, 17, 26, 28, 31, 33, 40,
                            ],
                        },
                        "line:141" => Outcome::Kekulize {
                            hydrogen_wrapper,
                            stage: SanitizeStage::Kekulize,
                            problem_atoms: vec![
                                4, 5, 6, 14, 15, 16, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                            ],
                        },
                        _ => unreachable!("all frozen source rows are classified"),
                    }
                };
                if actual != expected {
                    mismatches.push(format!(
                        "{id}, sanitize={sanitize}, remove_hydrogens={remove_hydrogens}: expected {expected:?}, observed {actual:?}"
                    ));
                }
                observations.push((id, sanitize, remove_hydrogens, actual));
            }
        }
    }

    assert_eq!(
        observations.len(),
        48,
        "all constructor cells were executed"
    );
    assert!(
        mismatches.is_empty(),
        "source error matrix mismatches after all 48 calls:\n{}",
        mismatches.join("\n")
    );
}

#[test]
fn uff_remaining_all_product_() {
    const INITIAL: [[f64; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [1.9, 0.2, 0.0],
        [5.0, 1.0, 0.0],
        [6.7, 1.1, 0.3],
    ];
    const IDS: [usize; 3] = [7, 3, 11];
    const SCALES: [f64; 3] = [1.0, 1.1, 1.2];

    fn coordinate_bits(molecule: &Molecule) -> Vec<(usize, Vec<[u64; 3]>)> {
        molecule
            .conformers_3d()
            .iter()
            .map(|conformer| {
                (
                    conformer.id(),
                    conformer
                        .coordinates()
                        .iter()
                        .map(|point| point.map(f64::to_bits))
                        .collect(),
                )
            })
            .collect()
    }

    fn coordinates_2d_bits(molecule: &Molecule) -> Vec<[u64; 2]> {
        molecule
            .coordinates_2d()
            .expect("the fixed source fixture retains one 2D row")
            .iter()
            .map(|point| point.map(f64::to_bits))
            .collect()
    }

    fn fixture(count: usize, warm_rings: bool) -> Molecule {
        let chemistry = Molecule::from_smiles("CC.CC").expect("fixed UFF graph parses");
        let conformers_3d = (0..count)
            .map(|row| {
                Conformer3D::new(
                    IDS[row],
                    INITIAL
                        .map(|point| point.map(|coordinate| coordinate * SCALES[row]))
                        .to_vec(),
                    true,
                )
            })
            .collect();
        let molecule = Molecule::from_parts(
            chemistry.topology().clone(),
            CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(
                    3,
                    vec![[9.0, -0.0], [8.0, 1.0], [7.0, 2.0], [6.0, 3.0]],
                )],
                conformers_3d,
                ..CoordinateBlock::default()
            },
            chemistry.properties().clone(),
        )
        .expect("fixed source coordinates fit CC.CC")
        .with_assigned_valence()
        .expect("fixed UFF graph has source valence state");
        if warm_rings {
            molecule
                .with_assigned_rings()
                .expect("the fixed acyclic graph has valid ring state")
        } else {
            molecule
        }
    }

    let mut optimization_calls = 0usize;
    let mut zero_conformer_errors = 0usize;
    let mut successful_calls = 0usize;
    for count in [0, 1, 3] {
        for iterations in [0, 1, 200, 1000] {
            for threshold in [0.0, 10.0, 100.0] {
                for ignore_interfragment_interactions in [false, true] {
                    for warm_rings in [false, true] {
                        for arc_peer in [false, true] {
                            let source = fixture(count, warm_rings);
                            let peer = arc_peer.then(|| source.clone());
                            if let Some(peer) = &peer {
                                assert!(std::ptr::eq(source.topology(), peer.topology()));
                                assert!(std::ptr::eq(source.properties(), peer.properties()));
                                assert_eq!(coordinate_bits(&source), coordinate_bits(peer));
                                assert_eq!(coordinates_2d_bits(&source), coordinates_2d_bits(peer));
                            }

                            let source_coordinates = coordinate_bits(&source);
                            let source_2d = coordinates_2d_bits(&source);
                            let source_topology = source.topology();
                            let source_properties = source.properties();
                            assert!(
                                source
                                    .uff_has_all_molecule_params()
                                    .expect("the prepared fixture retains its valence cache")
                            );

                            let params = UffConformerOptimizationParams {
                                max_iterations: iterations,
                                vdw_threshold: threshold,
                                ignore_interfragment_interactions,
                            };
                            let result = source.with_uff_optimized_conformers_with_params(&params);
                            optimization_calls += 1;

                            assert_eq!(coordinate_bits(&source), source_coordinates);
                            assert_eq!(coordinates_2d_bits(&source), source_2d);
                            assert!(std::ptr::eq(source_topology, source.topology()));
                            assert!(std::ptr::eq(source_properties, source.properties()));
                            if let Some(peer) = &peer {
                                assert_eq!(coordinate_bits(peer), source_coordinates);
                                assert_eq!(coordinates_2d_bits(peer), source_2d);
                                assert!(std::ptr::eq(source.topology(), peer.topology()));
                                assert!(std::ptr::eq(source.properties(), peer.properties()));
                                assert!(
                                    peer.uff_has_all_molecule_params()
                                        .expect("the Arc peer retains its prepared valence cache")
                                );
                            }

                            if count == 0 {
                                zero_conformer_errors += 1;
                                let error = match result {
                                    Err(error) => error,
                                    Ok(_) => panic!(
                                        "source confId=-1 construction must fail without stored 3D rows"
                                    ),
                                };
                                let OperationError::UffOptimization(error) = error else {
                                    panic!(
                                        "zero-3D failure must retain its typed UFF error: {error:?}"
                                    );
                                };
                                assert_eq!(
                                    error.kind(),
                                    UffOptimizationErrorKind::ConformerOptimization
                                );
                                assert!(std::error::Error::source(&error).is_some());
                                continue;
                            }

                            successful_calls += 1;
                            let result = result
                                .expect("stored 3D rows optimize through the public operation");
                            let output = &result.molecule;
                            let expected_ids = &IDS[..count];
                            assert_eq!(
                                result
                                    .conformers
                                    .iter()
                                    .map(|row| row.conformer_id)
                                    .collect::<Vec<_>>(),
                                expected_ids,
                                "stored source row order for count={count}, iterations={iterations}, threshold={threshold}, ignore={ignore_interfragment_interactions}, warm_rings={warm_rings}, arc_peer={arc_peer}"
                            );
                            assert!(result.conformers.iter().all(|row| matches!(
                                row.status,
                                0 | 1
                            )
                                && row.energy.is_finite()));
                            assert_eq!(
                                output
                                    .conformers_3d()
                                    .iter()
                                    .map(|row| row.id())
                                    .collect::<Vec<_>>(),
                                expected_ids
                            );
                            assert!(std::ptr::eq(source_topology, output.topology()));
                            assert!(std::ptr::eq(source_properties, output.properties()));
                            assert_eq!(coordinates_2d_bits(output), source_2d);
                            assert!(
                                output
                                    .uff_has_all_molecule_params()
                                    .expect("the result retains its prepared valence cache")
                            );
                            if iterations == 0 {
                                assert_eq!(coordinate_bits(output), source_coordinates);
                            }
                            for conformer in output.conformers_3d() {
                                assert_eq!(conformer.coordinates().len(), INITIAL.len());
                                assert!(
                                    conformer
                                        .coordinates()
                                        .iter()
                                        .flatten()
                                        .all(|coordinate| coordinate.is_finite())
                                );
                            }
                            if let Some(peer) = &peer {
                                assert!(std::ptr::eq(peer.topology(), output.topology()));
                                assert!(std::ptr::eq(peer.properties(), output.properties()));
                                assert_eq!(coordinate_bits(peer), source_coordinates);
                                assert_eq!(coordinates_2d_bits(peer), source_2d);
                            }
                        }
                    }
                }
            }
        }
    }

    assert_eq!(
        optimization_calls, 288,
        "execute the complete frozen product"
    );
    assert_eq!(
        zero_conformer_errors, 96,
        "all zero-3D cells are real calls"
    );
    assert_eq!(successful_calls, 192, "all nonempty cells are real calls");
}

#[test]
fn uff_remaining_all_reference_() {
    const INITIAL: [[f64; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [1.9, 0.2, 0.0],
        [5.0, 1.0, 0.0],
        [6.7, 1.1, 0.3],
    ];
    const IDS: [usize; 3] = [7, 3, 11];
    const SCALES: [f64; 3] = [1.0, 1.1, 1.2];
    const ONE_ITERATION_COORDINATE_BITS: [[[u64; 3]; 4]; 3] = [
        [
            [4598641110794173855, 4584026890252664668, 0],
            [4609993285745583870, 4595327571088798126, 0],
            [
                4617482136253981157,
                4607221623104106851,
                4583184447793597021,
            ],
            [
                4619062929510853352,
                4607593574458665015,
                4598605487821677197,
            ],
        ],
        [
            [4601035628630809199, 4586362091530870878, 0],
            [4610250340215625428, 4595544038010938386, 0],
            [
                4618178990165372210,
                4607703489880483840,
                4586949088337815106,
            ],
            [
                4619683378490468170,
                4608057463604035831,
                4598767838013282391,
            ],
        ],
        [
            [4603054482820045727, 4588378527603826956, 0],
            [4610507394685666986, 4595760504933078645, 0],
            [
                4618875844076763260,
                4608185356656860828,
                4589570881368892312,
            ],
            [
                4620303827470082988,
                4608521352749406646,
                4598930188204887588,
            ],
        ],
    ];

    fn coordinate_bits(molecule: &Molecule) -> Vec<(usize, Vec<[u64; 3]>)> {
        molecule
            .conformers_3d()
            .iter()
            .map(|conformer| {
                (
                    conformer.id(),
                    conformer
                        .coordinates()
                        .iter()
                        .map(|point| point.map(f64::to_bits))
                        .collect(),
                )
            })
            .collect()
    }

    fn coordinates_2d_bits(molecule: &Molecule) -> Vec<[u64; 2]> {
        molecule
            .coordinates_2d()
            .expect("the reference fixture retains one 2D row")
            .iter()
            .map(|point| point.map(f64::to_bits))
            .collect()
    }

    fn fixture(count: usize) -> Molecule {
        let chemistry = Molecule::from_smiles("CC.CC").expect("fixed UFF graph parses");
        let conformers_3d = (0..count)
            .map(|row| {
                Conformer3D::new(
                    IDS[row],
                    INITIAL
                        .map(|point| point.map(|coordinate| coordinate * SCALES[row]))
                        .to_vec(),
                    true,
                )
            })
            .collect();
        Molecule::from_parts(
            chemistry.topology().clone(),
            CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(
                    3,
                    vec![[9.0, -0.0], [8.0, 1.0], [7.0, 2.0], [6.0, 3.0]],
                )],
                conformers_3d,
                ..CoordinateBlock::default()
            },
            chemistry.properties().clone(),
        )
        .expect("fixed source coordinates fit CC.CC")
        .with_assigned_valence()
        .expect("fixed UFF graph has source valence state")
    }

    let cases: [(usize, i32, &[u64]); 4] = [
        (1, 0, &[4634709622640018448]),
        (1, 1, &[4622575255082429943]),
        (
            3,
            0,
            &[
                4634709622640018448,
                4640306196202131730,
                4644374702382500982,
            ],
        ),
        (
            3,
            1,
            &[
                4622575255082429943,
                4628472156411023569,
                4632421107459294556,
            ],
        ),
    ];
    for (count, iterations, expected_energy_bits) in cases {
        let source = fixture(count);
        let input_coordinate_bits = coordinate_bits(&source);
        let input_2d_bits = coordinates_2d_bits(&source);
        let topology = source.topology();
        let properties = source.properties();
        assert!(
            source
                .uff_has_all_molecule_params()
                .expect("the fixed fixture retains its valence cache")
        );

        let result = source
            .with_uff_optimized_conformers_with_params(&UffConformerOptimizationParams {
                max_iterations: iterations,
                vdw_threshold: 10.0,
                ignore_interfragment_interactions: true,
            })
            .expect("the fixed pinned source cases construct and optimize");
        let output = &result.molecule;
        assert_eq!(
            result
                .conformers
                .iter()
                .map(|row| row.conformer_id)
                .collect::<Vec<_>>(),
            IDS[..count]
        );
        assert!(result.conformers.iter().all(|row| row.status == 1));
        assert_eq!(
            result
                .conformers
                .iter()
                .map(|row| row.energy.to_bits())
                .collect::<Vec<_>>(),
            expected_energy_bits
        );
        assert!(std::ptr::eq(topology, output.topology()));
        assert!(std::ptr::eq(properties, output.properties()));
        assert_eq!(coordinates_2d_bits(output), input_2d_bits);
        assert!(
            output
                .uff_has_all_molecule_params()
                .expect("the result retains its valence cache")
        );

        let actual_coordinates = coordinate_bits(output);
        let expected_coordinates = if iterations == 0 {
            input_coordinate_bits.clone()
        } else {
            IDS[..count]
                .iter()
                .copied()
                .zip(ONE_ITERATION_COORDINATE_BITS[..count].iter())
                .map(|(id, row)| (id, row.to_vec()))
                .collect()
        };
        assert_eq!(actual_coordinates, expected_coordinates);
        assert_eq!(coordinate_bits(&source), input_coordinate_bits);
        assert_eq!(coordinates_2d_bits(&source), input_2d_bits);
    }
}

#[test]
fn uff_remaining_all_errors_() {
    use std::error::Error;

    #[derive(Debug, PartialEq, Eq)]
    struct CoordinateSnapshot {
        two_d: Option<Vec<[u64; 2]>>,
        three_d: Vec<(usize, Vec<[u64; 3]>)>,
    }

    fn coordinate_snapshot(molecule: &Molecule) -> CoordinateSnapshot {
        CoordinateSnapshot {
            two_d: molecule.coordinates_2d().map(|coordinates| {
                coordinates
                    .iter()
                    .map(|point| point.map(f64::to_bits))
                    .collect()
            }),
            three_d: molecule
                .conformers_3d()
                .iter()
                .map(|conformer| {
                    (
                        conformer.id(),
                        conformer
                            .coordinates()
                            .iter()
                            .map(|point| point.map(f64::to_bits))
                            .collect(),
                    )
                })
                .collect(),
        }
    }

    fn assert_error_and_storage(
        source: &Molecule,
        peer: &Molecule,
        expected_terminal: &str,
    ) -> OperationError {
        let source_coordinates = coordinate_snapshot(source);
        let peer_coordinates = coordinate_snapshot(peer);
        let source_topology = source.topology();
        let source_properties = source.properties();
        let source_two_d_storage = source.coordinates_2d().map(|rows| rows.as_ptr());
        let source_three_d_storage = source
            .conformers_3d()
            .iter()
            .map(|row| row.coordinates().as_ptr())
            .collect::<Vec<_>>();
        let peer_two_d_storage = peer.coordinates_2d().map(|rows| rows.as_ptr());
        let peer_three_d_storage = peer
            .conformers_3d()
            .iter()
            .map(|row| row.coordinates().as_ptr())
            .collect::<Vec<_>>();

        assert!(std::ptr::eq(source.topology(), peer.topology()));
        assert!(std::ptr::eq(source.properties(), peer.properties()));
        assert_eq!(source_two_d_storage, peer_two_d_storage);
        assert_eq!(source_three_d_storage, peer_three_d_storage);
        assert_eq!(source_coordinates, peer_coordinates);

        let error = source
            .with_uff_optimized_conformers_with_params(&UffConformerOptimizationParams {
                max_iterations: 0,
                vdw_threshold: 10.0,
                ignore_interfragment_interactions: true,
            })
            .expect_err("the fixed source failure must remain an error");
        let OperationError::UffOptimization(optimization) = &error else {
            panic!("source failure lost the typed UFF operation error: {error:?}");
        };
        assert_eq!(
            optimization.kind(),
            UffOptimizationErrorKind::ConformerOptimization
        );
        let direct_source = Error::source(optimization)
            .expect("all-conformer operation retains its concrete source error");
        assert!(
            direct_source
                .downcast_ref::<cosmolkit_forcefields::UffConformerError>()
                .is_some()
        );
        let mut terminal_source = direct_source;
        while let Some(next) = Error::source(terminal_source) {
            terminal_source = next;
        }
        let terminal_debug = format!("{terminal_source:?}");
        assert!(
            terminal_debug.contains(expected_terminal),
            "expected terminal source {expected_terminal:?}, got {terminal_debug}; operation={error:?}"
        );

        assert_eq!(coordinate_snapshot(source), source_coordinates);
        assert_eq!(coordinate_snapshot(peer), peer_coordinates);
        assert!(std::ptr::eq(source_topology, source.topology()));
        assert!(std::ptr::eq(source_properties, source.properties()));
        assert_eq!(
            source.coordinates_2d().map(|rows| rows.as_ptr()),
            source_two_d_storage
        );
        assert_eq!(
            source
                .conformers_3d()
                .iter()
                .map(|row| row.coordinates().as_ptr())
                .collect::<Vec<_>>(),
            source_three_d_storage
        );
        assert_eq!(
            peer.coordinates_2d().map(|rows| rows.as_ptr()),
            peer_two_d_storage
        );
        assert_eq!(
            peer.conformers_3d()
                .iter()
                .map(|row| row.coordinates().as_ptr())
                .collect::<Vec<_>>(),
            peer_three_d_storage
        );
        assert!(std::ptr::eq(source.topology(), peer.topology()));
        assert!(std::ptr::eq(source.properties(), peer.properties()));
        assert_eq!(coordinate_snapshot(peer), peer_coordinates);
        error
    }

    fn star_topology(
        center_atomic_number: u8,
        center_hybridization: Hybridization,
        neighbor_atomic_number: u8,
        neighbor_count: usize,
    ) -> cosmolkit::TopologyBlock {
        let mut atoms = Vec::with_capacity(neighbor_count + 1);
        for row in 0..=neighbor_count {
            let atomic_number = if row == 0 {
                center_atomic_number
            } else {
                neighbor_atomic_number
            };
            let hybridization = if row == 0 {
                center_hybridization
            } else {
                Hybridization::Sp3
            };
            let element = Element::from_atomic_number(atomic_number)
                .expect("the fixed source element is in the element table");
            atoms.push(Atom::from_spec(
                AtomId::new(row),
                AtomSpec::new(element)
                    .with_hybridization(hybridization)
                    .with_no_implicit(true),
            ));
        }
        let bonds = (0..neighbor_count)
            .map(|row| {
                Bond::from_spec(
                    BondId::new(row),
                    BondSpec::new(AtomId::new(0), AtomId::new(row + 1), BondOrder::Single),
                )
            })
            .collect();
        cosmolkit::TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("the fixed source star has valid row identity and endpoints")
    }

    fn prepared_molecule(
        topology: cosmolkit::TopologyBlock,
        two_d: Vec<[f64; 2]>,
        three_d: Vec<[f64; 3]>,
        conformer_id: usize,
        strict_valence: bool,
    ) -> Molecule {
        Molecule::from_parts(
            topology,
            CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(12, two_d)],
                conformers_3d: vec![Conformer3D::new(conformer_id, three_d, true)],
                ..CoordinateBlock::default()
            },
            MoleculeProperties::default(),
        )
        .expect("fixed source geometry fits the detached molecule blocks")
        .with_assigned_valence_with_params(&ValenceParams {
            model: ValenceModel::RdkitLike,
            strict: strict_valence,
        })
        .expect("the selected source cache policy creates the fixed fixture state")
    }

    fn star_coordinates(neighbor_count: usize) -> (Vec<[f64; 2]>, Vec<[f64; 3]>) {
        let two_d = (0..=neighbor_count)
            .map(|row| [row as f64, -(row as f64)])
            .collect();
        let three_d = (0..=neighbor_count)
            .map(|row| [row as f64 * 1.4, (row % 2) as f64 * 0.7, row as f64 * 0.2])
            .collect();
        (two_d, three_d)
    }

    fn peer_of(molecule: &Molecule) -> Molecule {
        molecule.clone()
    }

    let connected = Molecule::from_smiles("CC.CC")
        .expect("fixed no-coordinate source topology parses")
        .with_assigned_valence()
        .expect("fixed no-coordinate source topology has valid cached valence");
    let no_three_d_peer = peer_of(&connected);
    assert_error_and_storage(
        &connected,
        &no_three_d_peer,
        "SelectedThreeDimensionalConformerNotFound { conformer_id: 0 }",
    );

    let two_d_only = Molecule::from_parts(
        connected.topology().clone(),
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                7,
                vec![[17.0, -3.5], [16.0, -2.5], [15.0, -1.5], [14.0, -0.5]],
            )],
            ..CoordinateBlock::default()
        },
        connected.properties().clone(),
    )
    .expect("fixed 2D-only source blocks are valid")
    .with_assigned_valence()
    .expect("fixed 2D-only source topology retains valid cached valence");
    let two_d_peer = peer_of(&two_d_only);
    assert_error_and_storage(
        &two_d_only,
        &two_d_peer,
        "SelectedThreeDimensionalConformerNotFound { conformer_id: 0 }",
    );
    assert_eq!(two_d_only.coordinates_2d().unwrap()[0], [17.0, -3.5]);

    let (invalid_two_d, invalid_three_d) = star_coordinates(5);
    let invalid_valence = prepared_molecule(
        star_topology(6, Hybridization::Sp3, 9, 5),
        invalid_two_d,
        invalid_three_d,
        4,
        false,
    );
    let invalid_valence_peer = peer_of(&invalid_valence);
    let invalid_valence_error = assert_error_and_storage(
        &invalid_valence,
        &invalid_valence_peer,
        "InvalidValence { atom: AtomId(0), atomic_number: 6, formal_charge: 0, phase: Explicit, calculated: Some(5), reason: \"greater than permitted\"",
    );
    let OperationError::UffOptimization(optimization) = &invalid_valence_error else {
        unreachable!("the shared error helper already checks the operation variant")
    };
    let mut source_cursor = Error::source(optimization).expect("typed conformer source remains");
    let mut fragment_copy_error = None;
    let mut invalid_valence_source = None;
    while let Some(next) = Error::source(source_cursor) {
        if let Some(fragment_error) = next.downcast_ref::<cosmolkit_core::MoleculeFragmentsError>()
        {
            fragment_copy_error = Some(fragment_error);
        }
        if let Some(valence_error) = next.downcast_ref::<ValenceError>() {
            invalid_valence_source = Some(valence_error);
        }
        source_cursor = next;
    }
    assert_eq!(
        fragment_copy_error
            .expect("UFF preserves the source sanitized-fragment copy failure")
            .component_index(),
        Some(0)
    );
    assert!(matches!(
        invalid_valence_source.expect("fragment final sanitize preserves its core valence cause"),
        ValenceError::InvalidValence {
            atom,
            atomic_number: 6,
            formal_charge: 0,
            phase: cosmolkit_core::ValencePhase::Explicit,
            calculated: Some(5),
            reason: "greater than permitted",
            ..
        } if *atom == AtomId::new(0)
    ));

    let (tbp_two_d, mut tbp_three_d) = star_coordinates(5);
    tbp_three_d[1] = tbp_three_d[0];
    let degenerate_tbp = prepared_molecule(
        star_topology(15, Hybridization::Sp3d, 9, 5),
        tbp_two_d,
        tbp_three_d,
        4,
        true,
    );
    let degenerate_peer = peer_of(&degenerate_tbp);
    assert_error_and_storage(
        &degenerate_tbp,
        &degenerate_peer,
        "SourceDirectionVectorBelowTolerance { center_atom_index: 0, neighbor_atom_index: 1 }",
    );
}

#[test]
fn uff_remaining_query_boundary_() {
    // Fixed MolBlock emitted for the source C~C case. V3000 bond type 8 is
    // RDKit's BondNull query and must stay query-bearing at the SDF boundary.
    const SOURCE_CC_QUERY_MOLBLOCK: &str = "\n     RDKit          3D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 1.150762 -2.890078 3.560818 0\nM  V30 2 C 0.480135 -2.063470 3.748954 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 8 1 2\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";

    let record = SdfRecord::from_sdf(SOURCE_CC_QUERY_MOLBLOCK)
        .expect("source C~C MolBlock remains a readable SDF query record");
    assert!(matches!(record.graph(), SdfGraph::Query(_)));

    let query = record.query_graph().expect("query graph is preserved");
    assert_eq!(query.num_atoms(), 2);
    assert_eq!(query.num_bonds(), 1);
    assert_eq!(query.bonds()[0].bond().order(), BondOrder::Unspecified);
    assert_eq!(
        query.bonds()[0].predicate(),
        &QueryNode::predicate(BondQueryPredicate::Any)
    );

    let concrete_error = match Molecule::from_sdf(SOURCE_CC_QUERY_MOLBLOCK) {
        Ok(_) => panic!("query-bearing C~C MolBlock must not be coerced to Molecule"),
        Err(error) => error,
    };
    assert!(matches!(concrete_error, SdfError::QueryRecord));
    assert_eq!(
        concrete_error.to_string(),
        "query-bearing SDF record cannot be represented as a concrete molecule"
    );
}
