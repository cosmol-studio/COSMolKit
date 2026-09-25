use cosmolkit::{
    CoordinateDimension, Molecule, SdfCoordinateMode, SdfError, SdfGraph, SdfReadParams, SdfRecord,
};

fn v2000_atom(x: f64, y: f64, z: f64, symbol: &str) -> String {
    format!("{x:>10.4}{y:>10.4}{z:>10.4} {symbol:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
}

fn v2000_one_atom(symbol: &str, z: f64) -> String {
    format!(
        "one\n  COSMolKit         2D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n{}\nM  END\n",
        v2000_atom(1.25, -2.5, z, symbol)
    )
}

fn v3000_one_atom(symbol: &str, z: &str, attachment: &str) -> String {
    format!(
        "one\n  COSMolKit         2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 {symbol} 1.25 -2.5 {z} 0{attachment}\nM  V30 END ATOM\nM  V30 END CTAB\nM  END\n"
    )
}

#[test]
fn sdf_public_defaults_concrete_both_versions_and_first_record() {
    let defaults = SdfReadParams::default();
    assert!(defaults.sanitize);
    assert!(defaults.remove_hydrogens);
    assert!(defaults.strict_parsing);
    assert!(!defaults.expand_attachment_points);
    assert!(defaults.process_property_lists);
    assert_eq!(defaults.coordinate_mode, SdfCoordinateMode::Preserve);

    for input in [v2000_one_atom("C", 0.0), v3000_one_atom("C", "0", "")] {
        let stream = format!(
            "{input}>  <ID>\nfirst\n\n>  <ID>\nsecond\n\n$$$$\n{}$$$$\n",
            v2000_one_atom("O", 0.0)
        );
        let record = SdfRecord::from_sdf(&stream).expect("first complete record");
        assert!(matches!(record.graph(), SdfGraph::Molecule(_)));
        assert_eq!(record.molecule().unwrap().topology().atoms.len(), 1);
        assert_eq!(
            record.data_fields(),
            &[
                ("ID".into(), "first".into()),
                ("ID".into(), "second".into())
            ]
        );
        assert_eq!(record.properties().sdf_data_fields(), record.data_fields());
        assert_eq!(
            record.source_coordinate_dim(),
            Some(CoordinateDimension::TwoD)
        );
        assert_eq!(
            record.molecule().unwrap().coordinates_2d(),
            Some(&[[1.25, -2.5]][..])
        );
        let direct = Molecule::from_sdf(&stream).expect("concrete-only reader");
        assert_eq!(direct.topology().atoms.len(), 1);
    }
}

#[test]
fn sdf_public_coordinate_modes_preserve_one_selected_conformer_and_signed_zero() {
    for input in [v2000_one_atom("C", -3.5), v3000_one_atom("C", "-3.5", "")] {
        for (mode, dimension) in [
            (SdfCoordinateMode::Preserve, CoordinateDimension::ThreeD),
            (SdfCoordinateMode::Require2D, CoordinateDimension::TwoD),
            (SdfCoordinateMode::Require3D, CoordinateDimension::ThreeD),
        ] {
            let record = SdfRecord::from_sdf_with_params(
                &input,
                &SdfReadParams {
                    coordinate_mode: mode,
                    ..SdfReadParams::default()
                },
            )
            .unwrap();
            assert_eq!(record.source_coordinate_dim(), Some(dimension));
            let molecule = record.molecule().unwrap();
            match dimension {
                CoordinateDimension::TwoD => {
                    assert_eq!(molecule.coordinates_2d(), Some(&[[1.25, -2.5]][..]));
                    assert!(molecule.conformers_3d().is_empty());
                }
                CoordinateDimension::ThreeD => {
                    assert!(molecule.coordinates_2d().is_none());
                    assert_eq!(
                        molecule.conformers_3d()[0].coordinates(),
                        &[[1.25, -2.5, -3.5]]
                    );
                }
            }
        }
    }

    let negative_zero = v3000_one_atom("C", "-0.0", "");
    let record = SdfRecord::from_sdf(&negative_zero).unwrap();
    assert_eq!(
        record.source_coordinate_dim(),
        Some(CoordinateDimension::TwoD)
    );
    let xyz = &record.molecule().unwrap().conformers_3d()[0];
    assert!(!xyz.is_3d());
    assert_eq!(xyz.coordinates()[0][2].to_bits(), (-0.0f64).to_bits());
}

#[test]
fn sdf_public_query_preservation_wrong_kind_and_concrete_rejection() {
    for input in [v2000_one_atom("*", 0.0), v3000_one_atom("*", "0", "")] {
        let record = SdfRecord::from_sdf_with_params(
            &input,
            &SdfReadParams {
                sanitize: false,
                remove_hydrogens: false,
                ..SdfReadParams::default()
            },
        )
        .unwrap();
        assert!(matches!(record.graph(), SdfGraph::Query(_)));
        assert_eq!(record.query_graph().unwrap().num_atoms(), 1);
        assert!(matches!(
            record.molecule(),
            Err(SdfError::WrongGraphKind {
                expected: "molecule",
                actual: "query_graph"
            })
        ));
        assert!(matches!(
            Molecule::from_sdf_with_params(
                &input,
                &SdfReadParams {
                    sanitize: false,
                    remove_hydrogens: false,
                    ..SdfReadParams::default()
                }
            ),
            Err(SdfError::QueryRecord)
        ));
    }
    let concrete = SdfRecord::from_sdf(&v2000_one_atom("C", 0.0)).unwrap();
    assert!(matches!(
        concrete.query_graph(),
        Err(SdfError::WrongGraphKind {
            expected: "query_graph",
            actual: "molecule"
        })
    ));
}

#[test]
fn sdf_public_attachment_promotion_and_noop_are_classified_after_post() {
    let plain = v3000_one_atom("C", "0", "");
    let promoted = v3000_one_atom("C", "0", " ATTCHPT=1");
    let params = SdfReadParams {
        expand_attachment_points: true,
        sanitize: false,
        remove_hydrogens: false,
        ..SdfReadParams::default()
    };
    assert!(matches!(
        SdfRecord::from_sdf_with_params(&plain, &params)
            .unwrap()
            .graph(),
        SdfGraph::Molecule(_)
    ));
    let result = SdfRecord::from_sdf_with_params(&promoted, &params).unwrap();
    assert!(matches!(result.graph(), SdfGraph::Query(_)));
    assert_eq!(result.query_graph().unwrap().num_atoms(), 2);
    assert!(matches!(
        Molecule::from_sdf_with_params(&promoted, &params),
        Err(SdfError::QueryRecord)
    ));
}

#[test]
fn sdf_public_sanitize_and_hydrogen_parameter_branches() {
    let input = format!(
        "hydrogen\n  COSMolKit         2D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n{}\n{}\n  1  2  1  0  0  0  0\nM  END\n",
        v2000_atom(0.0, 0.0, 0.0, "C"),
        v2000_atom(1.0, 0.0, 0.0, "H")
    );
    for (sanitize, remove_hydrogens, expected_atoms) in [
        (false, false, 2),
        (false, true, 2),
        (true, false, 2),
        (true, true, 1),
    ] {
        let params = SdfReadParams {
            sanitize,
            remove_hydrogens,
            ..SdfReadParams::default()
        };
        let molecule = Molecule::from_sdf_with_params(&input, &params).unwrap();
        assert_eq!(
            molecule.topology().atoms.len(),
            expected_atoms,
            "sanitize={sanitize}, remove_hydrogens={remove_hydrogens}"
        );
        assert_eq!(molecule.coordinates_2d().unwrap().len(), expected_atoms);
    }
}

#[test]
fn sdf_public_property_lists_and_typed_sgroup_survive_runtime_install() {
    let input = format!(
        "props\n  COSMolKit         2D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n{}\n{}\n  1  2  1  0  0  0  0\nM  STY  1   1 DAT\nM  SAL   1  1   1\nM  END\n>  <atom.prop.Label>\nleft right\n\n$$$$\n",
        v2000_atom(0.0, 0.0, 0.0, "C"),
        v2000_atom(1.0, 0.0, 0.0, "O")
    );
    let applied = SdfRecord::from_sdf(&input).unwrap();
    assert_eq!(
        applied.molecule().unwrap().topology().atoms[0].prop("Label"),
        Some("left")
    );
    assert_eq!(
        applied.molecule().unwrap().topology().atoms[1].prop("Label"),
        Some("right")
    );
    assert_eq!(applied.properties().sdf_property_lists().len(), 1);
    assert_eq!(applied.substance_groups().len(), 1);
    assert_eq!(
        applied
            .molecule()
            .unwrap()
            .topology()
            .substance_groups
            .len(),
        1
    );

    let raw = SdfRecord::from_sdf_with_params(
        &input,
        &SdfReadParams {
            process_property_lists: false,
            ..SdfReadParams::default()
        },
    )
    .unwrap();
    assert_eq!(raw.properties().prop("atom.prop.Label"), Some("left right"));
    assert!(raw.properties().sdf_property_lists().is_empty());
    assert_eq!(
        raw.molecule().unwrap().topology().atoms[0].prop("Label"),
        None
    );
    assert_eq!(raw.substance_groups().len(), 1);
}

#[test]
fn sdf_public_errors_are_structured_and_do_not_modify_existing_values() {
    assert!(matches!(SdfRecord::from_sdf(""), Err(SdfError::Read(_))));
    assert!(matches!(
        Molecule::from_sdf("not a molfile"),
        Err(SdfError::Read(_))
    ));
    let bad_list = format!(
        "bad-list\n  COSMolKit         2D\n\n  2  0  0  0  0  0  0  0  0  0999 V2000\n{}\n{}\nM  END\n>  <atom.prop.Label>\nonly-one\n\n$$$$\n",
        v2000_atom(0.0, 0.0, 0.0, "C"),
        v2000_atom(1.0, 0.0, 0.0, "O")
    );
    assert!(matches!(
        SdfRecord::from_sdf(&bad_list),
        Err(SdfError::Read(_))
    ));
    let raw = SdfRecord::from_sdf_with_params(
        &bad_list,
        &SdfReadParams {
            process_property_lists: false,
            ..SdfReadParams::default()
        },
    )
    .unwrap();
    assert_eq!(raw.properties().prop("atom.prop.Label"), Some("only-one"));
    let source = Molecule::from_sdf(&v2000_one_atom("C", 0.0)).unwrap();
    let cloned = source.clone();
    let _ = Molecule::from_sdf("not a molfile");
    assert_eq!(source, cloned);
}
