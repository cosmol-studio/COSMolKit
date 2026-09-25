use std::collections::BTreeMap;

use cosmolkit_model::AtomId;
use cosmolkit_smiles::{SmilesParseError, SmilesParseParams, parse_smiles};
use cosmolkit_types::ChiralTag;

#[test]
fn defaults_match_the_pinned_v2_constructor() {
    let params = SmilesParseParams::default();
    assert!(params.sanitize);
    assert!(params.allow_cxsmiles);
    assert!(params.strict_cxsmiles);
    assert!(params.parse_name);
    assert!(params.remove_hydrogens);
    assert!(!params.skip_cleanup);
    assert!(!params.debug_parse);
    assert!(params.replacements.is_empty());
}

#[test]
fn zero_isotope_is_canonical_in_detached_parser_rows() {
    for (input, expected) in [("[0C]", None), ("[0H]", None), ("[13C]", Some(13))] {
        let record = parse_smiles(input, &Default::default()).unwrap();
        assert_eq!(record.topology.atoms[0].isotope(), expected, "{input}");
    }
}

#[test]
fn stereo_finalization_prefers_two_d_coordinates_and_keeps_stored_state() {
    use cosmolkit_model::{Conformer2D, Conformer3D};
    use cosmolkit_smiles::finalize_smiles_stereo;
    use cosmolkit_types::{BondDirection as D, BondStereo};
    let params = SmilesParseParams::default();
    let mut record = parse_smiles("CC=CC |c:1| sample", &params).unwrap();
    // Conflicting conformers prove first-2D precedence; the core geometry
    // consumes lifted XY without replacing the original coordinate tables.
    record.coordinates.conformers_2d.push(Conformer2D::new(
        3,
        vec![[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
    ));
    record.coordinates.conformers_3d.push(Conformer3D::new(
        8,
        vec![
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, -1.0, 0.0],
        ],
        true,
    ));
    let source = record.clone();
    let output = finalize_smiles_stereo(record, &params).unwrap();
    assert_eq!(output.coordinates, source.coordinates);
    assert_eq!(output.properties.name(), Some("sample"));
    assert_eq!(source.properties.prop("_needsDetectBondStereo"), Some("1"));
    assert_eq!(output.properties.prop("_needsDetectBondStereo"), None);
    assert_eq!(output.topology.bonds[1].stereo(), BondStereo::Z);
    assert_eq!(
        output
            .topology
            .bonds
            .iter()
            .map(|bond| bond.direction())
            .collect::<Vec<_>>(),
        [D::EndUpRight, D::None, D::EndDownRight]
    );
}

#[test]
fn stereo_finalization_propagates_typed_coordinate_failure() {
    use cosmolkit_core::DoubleBondStereoError;
    use cosmolkit_model::Conformer2D;
    use cosmolkit_smiles::{SmilesStereoError, finalize_smiles_stereo};
    let params = SmilesParseParams::default();
    let mut record = parse_smiles("CC=CC |c:1|", &params).unwrap();
    record
        .coordinates
        .conformers_2d
        .push(Conformer2D::new(0, vec![[0.0, 0.0]]));
    let error = finalize_smiles_stereo(record, &params).unwrap_err();
    assert!(matches!(
        error,
        SmilesStereoError::Directions(DoubleBondStereoError::InvalidConformer(_))
    ));
}

#[test]
fn empty_and_whitespace_graphs_are_valid_detached_records() {
    for input in ["", " \t\r\n"] {
        let record = parse_smiles(input, &Default::default()).expect("empty graph");
        assert!(record.topology.atoms.is_empty(), "{input:?}");
        assert!(record.topology.bonds.is_empty(), "{input:?}");
        assert!(record.coordinates.conformers_2d.is_empty(), "{input:?}");
        assert!(record.coordinates.conformers_3d.is_empty(), "{input:?}");
        assert_eq!(record.properties.name(), None, "{input:?}");
    }
}

#[test]
fn ordered_replacements_repeat_until_a_full_pass_is_stable() {
    let params = SmilesParseParams {
        replacements: BTreeMap::from([
            ("{A}".to_owned(), "C1CC1".to_owned()),
            ("{Q}".to_owned(), "{X}CC{X}".to_owned()),
            ("{X}".to_owned(), "N".to_owned()),
        ]),
        ..Default::default()
    };
    let record = parse_smiles("C{A}C{Q}C", &params).expect("replacement parse");
    assert_eq!(record.topology.atoms.len(), 10);
    // Expanded source: CC1CC1CNCCNC. Its nine traversal bonds plus the
    // C1...C1 ring closure produce ten bond rows.
    assert_eq!(record.topology.bonds.len(), 10);

    let unchanged = parse_smiles("CO", &params).expect("unchanged input");
    assert_eq!(unchanged.topology.atoms.len(), 2);

    let invalid = SmilesParseParams {
        replacements: BTreeMap::from([("{BAD}".to_owned(), "?".to_owned())]),
        ..Default::default()
    };
    assert!(matches!(
        parse_smiles("C{BAD}", &invalid),
        Err(SmilesParseError::Unsupported {
            token: '?',
            offset: 1
        })
    ));
}

#[test]
fn nonconvergent_and_empty_replacement_keys_fail_structurally() {
    for replacements in [
        BTreeMap::from([("{A}".to_owned(), "{A}C".to_owned())]),
        BTreeMap::from([
            ("{A}".to_owned(), "{B}".to_owned()),
            ("{B}".to_owned(), "{A}".to_owned()),
        ]),
    ] {
        let error = parse_smiles(
            "C{A}",
            &SmilesParseParams {
                replacements,
                ..Default::default()
            },
        )
        .unwrap_err();
        assert!(matches!(error, SmilesParseError::ReplacementCycle { .. }));
    }

    let error = parse_smiles(
        "C",
        &SmilesParseParams {
            replacements: BTreeMap::from([(String::new(), "N".to_owned())]),
            ..Default::default()
        },
    )
    .unwrap_err();
    assert_eq!(error, SmilesParseError::EmptyReplacementKey);
}

#[test]
fn branches_rings_and_stereo_finish_in_source_order() {
    let record = parse_smiles("F[C@]1(Cl)CCCC1", &Default::default()).expect("ring stereo");
    assert_eq!(record.topology.atoms.len(), 7);
    assert_eq!(record.topology.bonds.len(), 7);
    assert_eq!(
        record.topology.atoms[1].chiral_tag(),
        ChiralTag::TetrahedralCcw
    );
    let closure = &record.topology.bonds[6];
    assert_eq!(
        (closure.begin(), closure.end()),
        (AtomId::new(6), AtomId::new(1))
    );

    let nested = parse_smiles("CC(C)(N)O", &Default::default()).expect("branches");
    assert_eq!(nested.topology.atoms.len(), 5);
    assert_eq!(nested.topology.adjacency.neighbors_of(1).len(), 4);
}

#[test]
fn cx_and_name_policies_are_atomic_and_explicit() {
    let record =
        parse_smiles("CC |$left;right$| sample name", &Default::default()).expect("CX and name");
    assert_eq!(record.topology.atoms[0].prop("atomLabel"), Some("left"));
    assert_eq!(record.properties.name(), Some("sample name"));

    let no_cx = parse_smiles(
        "CC |plain name|",
        &SmilesParseParams {
            allow_cxsmiles: false,
            ..Default::default()
        },
    )
    .expect("name without CX");
    assert_eq!(no_cx.properties.name(), Some("|plain name|"));

    let strict = parse_smiles("CC |rb:0:0|", &Default::default()).unwrap_err();
    assert!(matches!(strict, SmilesParseError::UnsupportedCx(_)));
    let recovered = parse_smiles(
        "CC |rb:0:0| trailing",
        &SmilesParseParams {
            strict_cxsmiles: false,
            ..Default::default()
        },
    )
    .expect("non-strict CX recovery");
    assert_eq!(recovered.properties.prop("_CXSMILES_Data"), Some(""));
    assert_eq!(recovered.properties.name(), None);
}

#[test]
fn cleanup_debug_and_post_parse_chemistry_flags_keep_their_layer_boundaries() {
    let cleaned = parse_smiles("CC", &Default::default()).expect("cleaned");
    assert_eq!(cleaned.topology.bonds[0].prop("_cxsmilesBondIdx"), None);

    let retained = parse_smiles(
        "CC",
        &SmilesParseParams {
            skip_cleanup: true,
            ..Default::default()
        },
    )
    .expect("uncleaned");
    assert_eq!(
        retained.topology.bonds[0].prop("_cxsmilesBondIdx"),
        Some("0")
    );

    let baseline = parse_smiles(
        "[H]C",
        &SmilesParseParams {
            sanitize: true,
            remove_hydrogens: true,
            ..Default::default()
        },
    )
    .expect("detached parse");
    let deferred = parse_smiles(
        "[H]C",
        &SmilesParseParams {
            sanitize: false,
            remove_hydrogens: false,
            debug_parse: true,
            ..Default::default()
        },
    )
    .expect("detached parse without post chemistry");
    assert_eq!(baseline, deferred);
    assert_eq!(baseline.topology.atoms.len(), 2);
}

#[test]
fn errors_preserve_offsets_ring_context_and_final_model_validation() {
    assert!(matches!(
        parse_smiles("C?N", &Default::default()),
        Err(SmilesParseError::Unsupported {
            token: '?',
            offset: 1
        })
    ));
    assert_eq!(
        parse_smiles("C7CC", &Default::default()).unwrap_err(),
        SmilesParseError::UnclosedRing { index: 7 }
    );
    let coordinate = parse_smiles("C |(1e309,0)|", &Default::default()).unwrap_err();
    assert!(
        matches!(coordinate, SmilesParseError::Model(message) if message.contains("non-finite"))
    );
}
