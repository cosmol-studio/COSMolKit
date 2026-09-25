use std::io::Cursor;

use cosmolkit_io::{
    MolBlockRecord, MolPostParams, SdfCoordinateMode, SdfDataReadParams, SdfGraphDataset,
    SdfGraphReader, SdfReadError, read_sdf_graph_record_detached_with_params,
    read_sdf_record_detached, read_sdf_record_detached_with_params, read_sdf_records_detached,
};
use cosmolkit_model::{CoordinateDimension, SdfPropertyListTarget};
use cosmolkit_types::ChiralTag;

fn v2000_atom(x: f64, y: f64, z: f64, symbol: &str) -> String {
    format!("{x:>10.4}{y:>10.4}{z:>10.4} {symbol:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
}

fn v2000_one_atom(title: &str, info: &str, x: f64, y: f64, z: f64) -> String {
    format!(
        "{title}\n{info}\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n{}\nM  END\n",
        v2000_atom(x, y, z, "C")
    )
}

fn v3000_one_atom(title: &str, info: &str, x: &str, y: &str, z: &str) -> String {
    format!(
        concat!(
            "{title}\n{info}\n\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 1 0 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 C {x} {y} {z} 0\n",
            "M  V30 END ATOM\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        ),
        title = title,
        info = info,
        x = x,
        y = y,
        z = z,
    )
}

fn params(coordinate_mode: SdfCoordinateMode) -> SdfDataReadParams {
    SdfDataReadParams {
        coordinate_mode,
        ..SdfDataReadParams::default()
    }
}

#[test]
fn sdf_parse_state_handoff_v2000_header_wedge_and_query_dimension() {
    // RDKit 2026.03.1 MolFromMolDataStream/calculate3dFlag: an exact 3D
    // info-header label with zero Z remains 3D unless the parsed bond made
    // chiralityPossible true; nonzero Z takes precedence over that label.
    for (symbol, query) in [("C", false), ("*", true)] {
        let block = v2000_one_atom("header", "  COSMolKit         3D", 0.0, 0.0, 0.0)
            .replace(" C   0", &format!(" {symbol:<3} 0"));
        for (mode, expected) in [
            (SdfCoordinateMode::Preserve, CoordinateDimension::ThreeD),
            (SdfCoordinateMode::Require2D, CoordinateDimension::TwoD),
            (SdfCoordinateMode::Require3D, CoordinateDimension::ThreeD),
        ] {
            let record = read_sdf_graph_record_detached_with_params(&block, params(mode))
                .expect("marked 3D all-zero-Z record");
            let dim = match record.mol_block {
                MolBlockRecord::Concrete { coordinates, .. } if !query => {
                    coordinates.source_coordinate_dim
                }
                MolBlockRecord::Query(record) if query => record.source_coordinate_dim,
                other => panic!("unexpected graph kind for {symbol}: {other:?}"),
            };
            assert_eq!(dim, Some(expected), "symbol={symbol}, mode={mode:?}");
        }
    }

    let atoms = format!(
        "{}\n{}\n",
        v2000_atom(0.0, 0.0, 0.0, "C"),
        v2000_atom(1.0, 0.0, 0.0, "F")
    );
    for (stereo, expected) in [
        (0, CoordinateDimension::ThreeD),
        (1, CoordinateDimension::TwoD),
        (6, CoordinateDimension::TwoD),
    ] {
        let block = format!(
            "wedge\n  COSMolKit         3D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n{atoms}  1  2  1  {stereo}  0  0  0\nM  END\n"
        );
        let record =
            read_sdf_graph_record_detached_with_params(&block, params(SdfCoordinateMode::Preserve))
                .expect("marked 3D V2000 bond record");
        let MolBlockRecord::Concrete { coordinates, .. } = record.mol_block else {
            panic!("expected concrete record")
        };
        assert_eq!(
            coordinates.source_coordinate_dim,
            Some(expected),
            "stereo={stereo}"
        );
    }
    let nonzero = v2000_one_atom("z", "  COSMolKit         2D", 0.0, 0.0, 0.01);
    let record =
        read_sdf_graph_record_detached_with_params(&nonzero, params(SdfCoordinateMode::Preserve))
            .unwrap();
    let MolBlockRecord::Concrete { coordinates, .. } = record.mol_block else {
        panic!("expected concrete record")
    };
    assert_eq!(
        coordinates.source_coordinate_dim,
        Some(CoordinateDimension::ThreeD)
    );
}

#[test]
fn sdf_parse_state_handoff_v3000_cfg_and_mode_matrix() {
    // ParseV3000BondBlock sets chiralityPossible only for CFG=1/3. CFG=2
    // gives UNKNOWN direction but does not change the 3D-header decision.
    for (cfg, detected) in [
        (0, CoordinateDimension::ThreeD),
        (1, CoordinateDimension::TwoD),
        (2, CoordinateDimension::ThreeD),
        (3, CoordinateDimension::TwoD),
    ] {
        let block = format!(
            concat!(
                "cfg\n  COSMolKit         3D\n\n",
                "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
                "M  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\n",
                "M  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\n",
                "M  V30 2 F 1 0 0 0\nM  V30 END ATOM\n",
                "M  V30 BEGIN BOND\nM  V30 1 1 1 2 CFG={}\n",
                "M  V30 END BOND\nM  V30 END CTAB\nM  END\n"
            ),
            cfg
        );
        for (mode, expected) in [
            (SdfCoordinateMode::Preserve, detected),
            (SdfCoordinateMode::Require2D, CoordinateDimension::TwoD),
            (SdfCoordinateMode::Require3D, CoordinateDimension::ThreeD),
        ] {
            for (symbol, query) in [("C", false), ("*", true)] {
                let input = block.replace("M  V30 1 C ", &format!("M  V30 1 {symbol} "));
                let record = read_sdf_graph_record_detached_with_params(&input, params(mode))
                    .expect("V3000 CFG SDF record");
                let dim = match record.mol_block {
                    MolBlockRecord::Concrete { coordinates, .. } if !query => {
                        coordinates.source_coordinate_dim
                    }
                    MolBlockRecord::Query(record) if query => record.source_coordinate_dim,
                    other => panic!("unexpected graph kind for {symbol}: {other:?}"),
                };
                assert_eq!(
                    dim,
                    Some(expected),
                    "CFG={cfg}, mode={mode:?}, symbol={symbol}"
                );
            }
        }
    }
}

#[test]
fn sdf_parse_state_post_handoff_v2000_preserves_parser_stereo_bit() {
    // Fixed RDKit 2026.03.1: 3D header + zero Z + V2000 wedge is detected
    // as 2D and gives CHI_TETRAHEDRAL_CCW, including a query-bearing atom.
    let atoms = [
        v2000_atom(0.0, 0.0, 0.0, "C"),
        v2000_atom(1.0, 0.0, 0.0, "F"),
        v2000_atom(0.0, 1.0, 0.0, "Cl"),
        v2000_atom(-1.0, 0.0, 0.0, "Br"),
        v2000_atom(0.0, -1.0, 0.0, "I"),
    ]
    .join("\n");
    let block = format!(
        "v2\n  COSMolKit         3D\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n{atoms}\n  1  2  1  1  0  0  0\n  1  3  1  0  0  0  0\n  1  4  1  0  0  0  0\n  1  5  1  0  0  0  0\nM  END\n"
    );
    for (symbol, query) in [("C", false), ("*", true)] {
        let input = block.replacen(" C   0", &format!(" {symbol:<3} 0"), 1);
        for (mode, expected_tag) in [
            (SdfCoordinateMode::Preserve, ChiralTag::TetrahedralCcw),
            (SdfCoordinateMode::Require2D, ChiralTag::TetrahedralCcw),
            (SdfCoordinateMode::Require3D, ChiralTag::Unspecified),
        ] {
            let parsed = read_sdf_graph_record_detached_with_params(&input, params(mode)).unwrap();
            let finished = parsed
                .finish_mol_post(MolPostParams {
                    sanitize: false,
                    remove_hs: false,
                    expand_attachment_points: false,
                })
                .unwrap();
            let tag = match finished.mol_block {
                MolBlockRecord::Concrete { topology, .. } if !query => {
                    topology.atoms[0].chiral_tag()
                }
                MolBlockRecord::Query(record) if query => record.query.atoms()[0].chiral_tag(),
                other => panic!("unexpected graph kind: {other:?}"),
            };
            assert_eq!(tag, expected_tag, "symbol={symbol}, mode={mode:?}");
        }
    }
}

#[test]
fn sdf_parse_state_post_handoff_v3000_cfg_modes_and_errors() {
    for (cfg, source_tag) in [
        (0, ChiralTag::Unspecified),
        (1, ChiralTag::TetrahedralCcw),
        (2, ChiralTag::Unspecified),
        (3, ChiralTag::TetrahedralCw),
    ] {
        let block = format!(
            "cfg\n  COSMolKit         3D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 5 4 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 F 1 0 0 0\nM  V30 3 Cl 0 1 0 0\nM  V30 4 Br -1 0 0 0\nM  V30 5 I 0 -1 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2 CFG={cfg}\nM  V30 2 1 1 3\nM  V30 3 1 1 4\nM  V30 4 1 1 5\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n"
        );
        for (symbol, query) in [("C", false), ("*", true)] {
            let input = block.replacen("M  V30 1 C ", &format!("M  V30 1 {symbol} "), 1);
            for mode in [
                SdfCoordinateMode::Preserve,
                SdfCoordinateMode::Require2D,
                SdfCoordinateMode::Require3D,
            ] {
                let parsed =
                    read_sdf_graph_record_detached_with_params(&input, params(mode)).unwrap();
                let finished = parsed
                    .finish_mol_post(MolPostParams {
                        sanitize: false,
                        remove_hs: false,
                        expand_attachment_points: false,
                    })
                    .unwrap();
                let tag = match finished.mol_block {
                    MolBlockRecord::Concrete { topology, .. } if !query => {
                        topology.atoms[0].chiral_tag()
                    }
                    MolBlockRecord::Query(record) if query => record.query.atoms()[0].chiral_tag(),
                    other => panic!("unexpected graph kind: {other:?}"),
                };
                let expected = if mode == SdfCoordinateMode::Require3D {
                    ChiralTag::Unspecified
                } else {
                    source_tag
                };
                assert_eq!(tag, expected, "CFG={cfg}, symbol={symbol}, mode={mode:?}");
            }
        }
    }
    let invalid = v3000_one_atom("bad", "", "0", "0", "0").replace("COUNTS 1 0", "COUNTS X 0");
    assert!(matches!(
        read_sdf_graph_record_detached_with_params(&invalid, params(SdfCoordinateMode::Preserve)),
        Err(SdfReadError::Field { .. }) | Err(SdfReadError::Parse(_))
    ));
}

#[test]
fn sdf_concrete_boundary_rejects_queries_structurally_and_preserves_record_context() {
    // CK's concrete/query distinction is a type contract, not an RDKit parse
    // failure: the full reader accepts the same input without dropping predicates.
    let concrete = format!("{}$$$$\n", v3000_one_atom("ordinary", "", "1", "2", "3"));
    let query = format!(
        "{}>  <ID>\nquery-id\n\n$$$$\n",
        v3000_one_atom("query", "", "1", "2", "3").replace("M  V30 1 C ", "M  V30 1 * ")
    );
    let record =
        read_sdf_graph_record_detached_with_params(&query, SdfDataReadParams::default()).unwrap();
    assert!(matches!(&record.mol_block, MolBlockRecord::Query(_)));
    assert_eq!(
        record.data_fields,
        [("ID".to_owned(), "query-id".to_owned())]
    );
    assert!(matches!(
        record.into_concrete(),
        Err(SdfReadError::QueryRecord)
    ));
    assert!(matches!(
        read_sdf_record_detached(&query),
        Err(SdfReadError::QueryRecord)
    ));
    let stream = format!("{concrete}{query}");
    match read_sdf_records_detached(&stream).unwrap_err() {
        SdfReadError::Record {
            index,
            byte_offset,
            source,
            ..
        } => {
            assert_eq!(index, 1);
            assert_eq!(byte_offset, concrete.len() as u64);
            assert!(matches!(*source, SdfReadError::QueryRecord));
        }
        error => panic!("expected indexed query type error, got {error:?}"),
    }
    let mut reader = SdfGraphReader::new(Cursor::new(stream));
    assert!(matches!(
        reader.next_record().unwrap().unwrap().mol_block,
        MolBlockRecord::Concrete { .. }
    ));
    assert!(matches!(
        reader.next_record().unwrap().unwrap().mol_block,
        MolBlockRecord::Query(_)
    ));
    assert!(reader.next_record().unwrap().is_none());
}

#[test]
fn sdf_read_v2000_default_preserves_detected_2d_and_ordered_fields() {
    // Pinned ForwardSDMolSupplier::_next parses exactly one MolBlock and then
    // readMolProps installs fields in encounter order. The default coordinate
    // policy retains calculate3dFlag's effective dimension.
    let input = format!(
        "{}>  <ID>\nfirst\n\n>  <NOTE>\nalpha\nbeta\n\n$$$$\n",
        v2000_one_atom("v2", "  COSMolKit         2D", 1.25, -2.5, 0.0)
    );
    let record = read_sdf_record_detached(&input).expect("V2000 SDF record");

    assert_eq!(record.topology.atoms.len(), 1);
    assert_eq!(
        record.coordinates.source_coordinate_dim,
        Some(CoordinateDimension::TwoD)
    );
    assert_eq!(record.coordinates.conformers_2d.len(), 1);
    assert!(record.coordinates.conformers_3d.is_empty());
    assert_eq!(
        record.coordinates.conformers_2d[0].coordinates(),
        &[[1.25, -2.5]]
    );
    assert_eq!(
        record.data_fields,
        [
            ("ID".to_owned(), "first".to_owned()),
            ("NOTE".to_owned(), "alpha\nbeta".to_owned()),
        ]
    );
    assert_eq!(record.properties.sdf_data_fields(), record.data_fields);
}

#[test]
fn sdf_read_v2000_coordinate_modes_select_and_convert_one_conformer() {
    let source_2d = v2000_one_atom("v2-2d", "  COSMolKit         2D", 1.0, 2.0, 0.0);
    let promoted =
        read_sdf_record_detached_with_params(&source_2d, params(SdfCoordinateMode::Require3D))
            .expect("promote V2000 XY to XYZ");
    assert!(promoted.coordinates.conformers_2d.is_empty());
    assert_eq!(
        promoted.coordinates.source_coordinate_dim,
        Some(CoordinateDimension::ThreeD)
    );
    assert_eq!(promoted.coordinates.conformers_3d.len(), 1);
    assert!(promoted.coordinates.conformers_3d[0].is_3d());
    assert_eq!(
        promoted.coordinates.conformers_3d[0].coordinates(),
        &[[1.0, 2.0, 0.0]]
    );

    let source_3d = v2000_one_atom("v2-3d", "  COSMolKit         3D", 1.0, 2.0, -3.5);
    let preserved =
        read_sdf_record_detached_with_params(&source_3d, params(SdfCoordinateMode::Preserve))
            .expect("preserve V2000 XYZ");
    assert_eq!(
        preserved.coordinates.source_coordinate_dim,
        Some(CoordinateDimension::ThreeD)
    );
    assert_eq!(
        preserved.coordinates.conformers_3d[0].coordinates(),
        &[[1.0, 2.0, -3.5]]
    );

    let projected =
        read_sdf_record_detached_with_params(&source_3d, params(SdfCoordinateMode::Require2D))
            .expect("project V2000 XYZ to XY");
    assert_eq!(
        projected.coordinates.source_coordinate_dim,
        Some(CoordinateDimension::TwoD)
    );
    assert_eq!(
        projected.coordinates.conformers_2d[0].coordinates(),
        &[[1.0, 2.0]]
    );
    assert!(projected.coordinates.conformers_3d.is_empty());
}

#[test]
fn sdf_read_v3000_modes_preserve_xyz_bits_and_apply_to_query_records() {
    // A negative-zero Z cannot be represented by the XY carrier without bit
    // loss. Preserve keeps the XYZ carrier and its independent effective-2D
    // flag; explicit modes select exactly one destination dimension.
    let source = v3000_one_atom("v3", "  COSMolKit         2D", "1e0", "-2.5", "-0.0");
    let preserved =
        read_sdf_record_detached_with_params(&source, params(SdfCoordinateMode::Preserve))
            .expect("preserve V3000 negative-zero Z");
    assert_eq!(
        preserved.coordinates.source_coordinate_dim,
        Some(CoordinateDimension::TwoD)
    );
    assert!(preserved.coordinates.conformers_2d.is_empty());
    let xyz = &preserved.coordinates.conformers_3d[0];
    assert!(!xyz.is_3d());
    assert_eq!(xyz.coordinates()[0][2].to_bits(), (-0.0f64).to_bits());

    let forced_2d =
        read_sdf_record_detached_with_params(&source, params(SdfCoordinateMode::Require2D))
            .expect("force V3000 XY");
    assert_eq!(
        forced_2d.coordinates.conformers_2d[0].coordinates(),
        &[[1.0, -2.5]]
    );
    assert!(forced_2d.coordinates.conformers_3d.is_empty());

    let query = source.replace("M  V30 1 C ", "M  V30 1 * ");
    let query =
        read_sdf_graph_record_detached_with_params(&query, params(SdfCoordinateMode::Require3D))
            .expect("force query record to 3D");
    let MolBlockRecord::Query(query) = query.mol_block else {
        panic!("wildcard V3000 atom must remain a query graph");
    };
    assert_eq!(
        query.source_coordinate_dim,
        Some(CoordinateDimension::ThreeD)
    );
    assert_eq!(query.query.conformers_3d().len(), 1);
    assert!(query.query.conformers_3d()[0].is_3d());
    assert_eq!(
        query.query.conformers_3d()[0].coordinates()[0][2].to_bits(),
        (-0.0f64).to_bits()
    );
}

#[test]
fn sdf_read_property_lists_preserve_raw_fields_and_obey_processing_switch() {
    let mol = format!(
        concat!(
            "props\n  COSMolKit         2D\n\n",
            "  2  1  0  0  0  0  0  0  0  0999 V2000\n",
            "{}\n{}\n",
            "  1  2  1  0  0  0  0\n",
            "M  END\n",
            ">  <atom.prop.Label>\nleft right\n\n",
            ">  <bond.iprop.OrderTag>\n7\n\n",
            "$$$$\n",
        ),
        v2000_atom(0.0, 0.0, 0.0, "C"),
        v2000_atom(1.0, 0.0, 0.0, "O"),
    );
    let applied = read_sdf_record_detached(&mol).expect("apply property lists");
    assert_eq!(applied.topology.atoms[0].prop("Label"), Some("left"));
    assert_eq!(applied.topology.atoms[1].prop("Label"), Some("right"));
    assert_eq!(applied.topology.bonds[0].prop("OrderTag"), Some("7"));
    assert_eq!(applied.properties.sdf_property_lists().len(), 2);
    assert_eq!(
        applied.properties.sdf_property_lists()[0].target(),
        SdfPropertyListTarget::Atom
    );
    assert_eq!(
        applied.properties.prop("atom.prop.Label"),
        Some("left right")
    );

    let raw_only = read_sdf_record_detached_with_params(
        &mol,
        SdfDataReadParams {
            process_property_lists: false,
            ..SdfDataReadParams::default()
        },
    )
    .expect("preserve raw lists without expansion");
    assert_eq!(
        raw_only.properties.prop("atom.prop.Label"),
        Some("left right")
    );
    assert_eq!(raw_only.topology.atoms[0].prop("Label"), None);
    assert_eq!(raw_only.topology.bonds[0].prop("OrderTag"), None);
    assert!(raw_only.properties.sdf_property_lists().is_empty());
}

#[test]
fn sdf_read_detached_boundary_does_not_run_sanitization_or_hydrogen_removal() {
    // MolFromMolDataStream parsing precedes finishMolProcessing. This detached
    // IO boundary intentionally returns the parsed six-atom/five-bond graph;
    // chemistry finalization belongs to the live cosmolkit wrapper.
    let mut atoms = String::new();
    atoms.push_str("M  V30 1 C 0 0 0 0\n");
    for index in 2..=6 {
        atoms.push_str(&format!("M  V30 {index} H {index} 0 0 0\n"));
    }
    let mut bonds = String::new();
    for index in 1..=5 {
        bonds.push_str(&format!("M  V30 {index} 1 1 {}\n", index + 1));
    }
    let input = format!(
        concat!(
            "unsanitized\n  COSMolKit         2D\n\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\nM  V30 COUNTS 6 5 0 0 0\n",
            "M  V30 BEGIN ATOM\n{}M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n{}M  V30 END BOND\n",
            "M  V30 END CTAB\nM  END\n$$$$\n",
        ),
        atoms, bonds
    );
    let record = read_sdf_record_detached(&input).expect("detached hypervalent graph");
    assert_eq!(record.topology.atoms.len(), 6);
    assert_eq!(record.topology.bonds.len(), 5);
    assert_eq!(
        record
            .topology
            .atoms
            .iter()
            .filter(|atom| atom.element().symbol() == "H")
            .count(),
        5
    );
}

#[test]
fn sdf_read_stream_handles_empty_single_multiple_crlf_and_line_start_delimiters() {
    assert_eq!(read_sdf_records_detached(""), Err(SdfReadError::Empty));

    let one = v2000_one_atom("one", "  COSMolKit         2D", 0.0, 0.0, 0.0);
    let only = read_sdf_records_detached(&one).expect("final record without delimiter");
    assert_eq!(only.len(), 1);
    assert_eq!(only[0].properties.name(), Some("one"));

    let first = format!("{one}>  <TEXT>\r\ninside $$$$ value\r\n\r\n$$$$ suffix\r\n");
    let second = format!(
        "{}$$$$\r\n",
        v2000_one_atom("two", "  COSMolKit         2D", 1.0, 0.0, 0.0).replace('\n', "\r\n")
    );
    let records = read_sdf_records_detached(&format!("{first}{second}"))
        .expect("mixed newline framed stream");
    assert_eq!(records.len(), 2);
    assert_eq!(
        records[0].data_fields,
        [("TEXT".to_owned(), "inside $$$$ value".to_owned())]
    );
    assert_eq!(records[1].properties.name(), Some("two"));
}

#[test]
fn sdf_graph_reader_consumes_failed_record_and_reports_stable_offsets() {
    let bad = "bad\n  COSMolKit\n\nnot-a-counts-line\n$$$$\n";
    let good = format!(
        "{}$$$$\n",
        v3000_one_atom("good", "  COSMolKit         2D", "1", "2", "0")
    );
    let stream = format!("{bad}{good}");
    let mut reader = SdfGraphReader::new(Cursor::new(stream.as_bytes()));

    let error = reader.next_record().expect_err("first record must fail");
    assert!(matches!(
        error,
        SdfReadError::Record {
            index: 0,
            byte_offset: 0,
            line_offset: 0,
            ..
        }
    ));
    assert_eq!(reader.records_consumed(), 1);
    assert_eq!(reader.bytes_consumed(), bad.len() as u64);
    assert_eq!(reader.lines_consumed(), bad.lines().count());

    let recovered = reader
        .next_record()
        .expect("reader remains aligned")
        .expect("second record exists");
    let MolBlockRecord::Concrete { properties, .. } = recovered.mol_block else {
        panic!("ordinary record must be concrete");
    };
    assert_eq!(properties.name(), Some("good"));
    assert_eq!(reader.records_consumed(), 2);
    assert!(reader.next_record().expect("EOF").is_none());
    assert!(reader.is_end());
}

#[test]
fn sdf_graph_dataset_preserves_metadata_direct_access_and_error_index() {
    let first = format!(
        "{}$$$$\n",
        v2000_one_atom("first", "  COSMolKit         2D", 0.0, 0.0, 0.0)
    );
    let second = format!(
        "{}>  <ID>\n2\n\n$$$$\n",
        v2000_one_atom("second", "  COSMolKit         2D", 1.0, 0.0, 0.0)
    );
    let file = tempfile::NamedTempFile::new().expect("temporary SDF");
    std::fs::write(file.path(), format!("{first}{second}")).expect("write SDF fixture");

    let dataset = SdfGraphDataset::open(file.path()).expect("index SDF");
    assert_eq!(dataset.len(), 2);
    assert!(!dataset.is_empty());
    assert_eq!(dataset.metadata(0).unwrap().index, 0);
    assert_eq!(dataset.metadata(0).unwrap().title.as_deref(), Some("first"));
    assert_eq!(dataset.metadata(1).unwrap().byte_offset, first.len() as u64);
    assert_eq!(
        dataset.metadata(1).unwrap().line_offset,
        first.lines().count()
    );
    assert_eq!(dataset.record_text(1).unwrap(), second);
    assert_eq!(
        dataset.record(1).unwrap().data_fields,
        [("ID".to_owned(), "2".to_owned())]
    );
    assert_eq!(
        dataset.record(2),
        Err(SdfReadError::RecordIndexOutOfRange {
            index: 2,
            record_count: 2,
        })
    );
}
