use cosmolkit_io::{
    MolBlockReadParams, MolBlockRecord, QueryMolBlockRecord, SdfReadError, read_mol_block_detached,
    read_mol_block_detached_with_params,
};
use cosmolkit_model::{AtomQueryPredicate, BondQueryPredicate, QueryNode};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo};

#[allow(clippy::too_many_arguments)]
fn atom_line(
    x: f64,
    y: f64,
    z: f64,
    symbol: &str,
    mass_diff: i32,
    charge: i32,
    parity: i32,
    h_count: i32,
    stereo_care: i32,
    total_valence: i32,
    reaction_role: i32,
    reaction_component: i32,
    atom_map: i32,
    inversion: i32,
    exact_change: i32,
) -> String {
    format!(
        "{x:>10.4}{y:>10.4}{z:>10.4} {symbol:<3}{mass_diff:>2}{charge:>3}{parity:>3}{h_count:>3}{stereo_care:>3}{total_valence:>3}{:>3}{reaction_role:>3}{reaction_component:>3}{atom_map:>3}{inversion:>3}{exact_change:>3}",
        0
    )
}

fn plain_atom(symbol: &str) -> String {
    atom_line(0.0, 0.0, 0.0, symbol, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
}

fn bond_line(
    begin: usize,
    end: usize,
    kind: u32,
    stereo: u32,
    topology: i32,
    reaction: i32,
) -> String {
    format!(
        "{begin:>3}{end:>3}{kind:>3}{stereo:>3}{:>3}{topology:>3}{reaction:>3}",
        0
    )
}

fn molblock(title: &str, atoms: &[String], bonds: &[String]) -> String {
    format!(
        "{title}\n  pinned          2D\ncomment\n{:>3}{:>3}  0  0  1  0  0  0  0  0999 V2000\n{}{}M  END\n",
        atoms.len(),
        bonds.len(),
        atoms
            .iter()
            .map(|line| format!("{line}\n"))
            .collect::<String>(),
        bonds
            .iter()
            .map(|line| format!("{line}\n"))
            .collect::<String>(),
    )
}

fn concrete(
    block: &str,
) -> (
    cosmolkit_model::TopologyBlock,
    cosmolkit_model::CoordinateBlock,
) {
    match read_mol_block_detached(block).expect("concrete V2000 record") {
        MolBlockRecord::Concrete {
            topology,
            coordinates,
            ..
        } => (topology, coordinates),
        MolBlockRecord::Query(_) => panic!("expected concrete V2000 record"),
    }
}

#[test]
fn coordinate_contract_v2000_screening_then_prefix_conversion() {
    // RDKit 2026.03.1 MolFileParser.cpp::FileParserUtils::toDouble screens
    // characters first, then calls atof; it does not require full consumption.
    // Observed with the glibc 2.43-2ubuntu2.4/C-locale reference profile.
    for (field, expected) in [
        (" 1.25     ", 1.25_f64),
        ("-0.0000   ", -0.0),
        ("          ", 0.0),
        ("1 2       ", 1.0),
        ("1.2.3     ", 1.2),
        ("--1       ", 0.0),
        ("1,5       ", 1.0),
    ] {
        assert_eq!(field.len(), 10);
        let mut atom = plain_atom("C");
        atom.replace_range(..10, field);
        let (_, coordinates) = concrete(&molblock("coordinate contract", &[atom], &[]));
        assert_eq!(
            coordinates.conformers_2d[0].coordinates()[0][0].to_bits(),
            expected.to_bits(),
            "{field:?}"
        );
    }
    for field in ["1e2", "\t1.0", "abc"] {
        let mut atom = plain_atom("C");
        atom.replace_range(..10, &format!("{field:<10}"));
        let error = read_mol_block_detached(&molblock("coordinate contract", &[atom], &[]))
            .expect_err("V2000 rejects source-forbidden characters before atof");
        assert!(
            matches!(error, SdfReadError::Parse(_)),
            "{field:?}: {error:?}"
        );
    }
}

fn query(block: &str) -> QueryMolBlockRecord {
    match read_mol_block_detached(block).expect("query V2000 record") {
        MolBlockRecord::Query(record) => record,
        MolBlockRecord::Concrete { .. } => panic!("expected query V2000 record"),
    }
}

fn atom_query_contains(node: &QueryNode<AtomQueryPredicate>, wanted: &AtomQueryPredicate) -> bool {
    match node {
        QueryNode::Predicate(found) => found == wanted,
        QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => children
            .iter()
            .any(|child| atom_query_contains(child, wanted)),
        QueryNode::Not(child) => atom_query_contains(child, wanted),
    }
}

#[test]
fn counts_version_empty_tables_and_truncated_sections_preserve_error_categories() {
    assert_eq!(read_mol_block_detached(""), Err(SdfReadError::Empty));
    assert_eq!(
        read_mol_block_detached("x\ninfo\ncomment\n  1\n"),
        Err(SdfReadError::Counts)
    );

    let bad_count = "x\ninfo\ncomment\n xx  0  0  0  0  0  0  0  0  0999 V2000\n";
    assert!(
        read_mol_block_detached(bad_count)
            .unwrap_err()
            .to_string()
            .contains("Cannot convert ' xx' to unsigned int on line 4")
    );

    let empty = molblock("empty", &[], &[]);
    let (topology, coordinates) = concrete(&empty);
    assert!(topology.atoms.is_empty());
    assert!(topology.bonds.is_empty());
    assert_eq!(coordinates.conformers_2d.len(), 1);
    assert!(coordinates.conformers_2d[0].coordinates().is_empty());

    let atom_eof = "x\ninfo\ncomment\n  1  0  0  0  0  0  0  0  0  0999 V2000";
    assert_eq!(
        read_mol_block_detached(atom_eof).unwrap_err(),
        SdfReadError::Parse("EOF hit while reading atoms".to_owned())
    );
    let bond_eof = format!(
        "x\ninfo\ncomment\n  1  1  0  0  0  0  0  0  0  0999 V2000\n{}",
        plain_atom("C")
    );
    assert_eq!(
        read_mol_block_detached(&bond_eof).unwrap_err(),
        SdfReadError::Parse("EOF hit while reading bonds".to_owned())
    );
}

#[test]
fn strict_version_symbol_and_fixed_width_rules_match_the_pinned_parser() {
    let malformed = "x\ninfo\ncomment\n  0  0  0  0  0  0  0  0  0  0999 X2000\nM  END\n";
    assert!(
        read_mol_block_detached(malformed)
            .unwrap_err()
            .to_string()
            .contains("CTAB version string invalid")
    );
    assert!(matches!(
        read_mol_block_detached_with_params(
            malformed,
            MolBlockReadParams {
                strict_parsing: false
            }
        ),
        Ok(MolBlockRecord::Concrete { .. })
    ));
    assert!(
        read_mol_block_detached(&malformed.replace("X2000", "V4000"))
            .unwrap_err()
            .to_string()
            .contains("Unsupported CTAB version: 'V4000'")
    );

    let short_atom = format!("{:>10}{:>10}{:>10} C", "0", "0", "0");
    assert_eq!(short_atom.len(), 32);
    let short = molblock("short", &[short_atom], &[]);
    assert!(
        read_mol_block_detached(&short)
            .unwrap_err()
            .to_string()
            .contains("Atom line too short")
    );
    let relaxed = read_mol_block_detached_with_params(
        &short,
        MolBlockReadParams {
            strict_parsing: false,
        },
    )
    .expect("non-strict 32-byte atom line");
    assert!(matches!(relaxed, MolBlockRecord::Concrete { .. }));

    let unknown = molblock("unknown", &[plain_atom("Zz")], &[]);
    assert!(
        read_mol_block_detached(&unknown)
            .unwrap_err()
            .to_string()
            .contains("Element 'Zz' not found")
    );
    let relaxed = read_mol_block_detached_with_params(
        &unknown,
        MolBlockReadParams {
            strict_parsing: false,
        },
    )
    .expect("non-strict unknown symbol");
    let MolBlockRecord::Concrete { topology, .. } = relaxed else {
        panic!("unknown non-strict symbol is a concrete dummy")
    };
    assert_eq!(topology.atoms[0].prop("dummyLabel"), Some("Zz"));

    let mut bad_atom = atom_line(0.0, 0.0, 0.0, "Zz", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    bad_atom.replace_range(34..36, "xx");
    let bad_mass_and_symbol = molblock("source error order", &[bad_atom], &[]);
    assert!(
        read_mol_block_detached(&bad_mass_and_symbol)
            .unwrap_err()
            .to_string()
            .contains("Cannot convert 'xx' to int on line 5")
    );
}

#[test]
fn ordinary_shorthand_dummy_and_rgroup_symbols_keep_source_state() {
    for (symbol, atomic_number, isotope, label) in [
        ("CL", 17, None, None),
        ("D", 1, Some(2), None),
        ("T", 1, Some(3), None),
        ("Pol", 0, None, Some("Pol")),
        ("Mod", 0, None, Some("Mod")),
        ("L", 0, None, None),
        ("LP", 0, None, None),
        ("R#", 0, None, Some("R#")),
        ("R0", 0, Some(0), Some("R0")),
        ("R12", 0, Some(12), Some("R12")),
        ("R99", 0, Some(99), Some("R99")),
    ] {
        let (topology, _) = concrete(&molblock(symbol, &[plain_atom(symbol)], &[]));
        let atom = &topology.atoms[0];
        assert_eq!(atom.element().atomic_number(), atomic_number, "{symbol}");
        assert_eq!(atom.isotope(), isotope, "{symbol}");
        assert_eq!(atom.prop("dummyLabel"), label, "{symbol}");
        assert!(!atom.no_implicit(), "{symbol}");
    }
}

#[test]
fn complex_wildcard_and_generic_symbols_remain_typed_queries() {
    for symbol in ["*", "R"] {
        let record = query(&molblock(symbol, &[plain_atom(symbol)], &[]));
        let atom = &record.query.atoms()[0];
        assert_eq!(
            atom.predicate(),
            &QueryNode::predicate(AtomQueryPredicate::Any),
            "{symbol}"
        );
        assert!(atom.atom().no_implicit(), "{symbol}");
        assert_eq!(
            atom.atom().prop("dummyLabel"),
            (symbol == "R").then_some("R"),
            "{symbol}"
        );
    }

    for symbol in ["Q", "QH", "A", "AH", "X", "XH", "M", "MH"] {
        let record = query(&molblock(symbol, &[plain_atom(symbol)], &[]));
        assert!(record.query.atoms()[0].atom().no_implicit(), "{symbol}");
    }

    for symbol in [
        "G", "GH", "G*", "GH*", "ALK", "ALH", "AEL", "AEH", "AYL", "AYH", "CBC", "CBH", "CAL",
        "CAH", "CEL", "CEH", "ARY", "ARH", "CYC", "CYH", "ACY", "ACH", "ABC", "ABH", "AHC", "AHH",
        "AOX", "AOH", "CHC", "CHH", "HAR", "HAH", "CXX", "CXH",
    ] {
        let record = query(&molblock(symbol, &[plain_atom(symbol)], &[]));
        let atom = &record.query.atoms()[0];
        assert_eq!(
            atom.predicate(),
            &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(0)),
            "{symbol}"
        );
        assert_eq!(atom.atom().prop("atomLabel"), Some(symbol), "{symbol}");
        assert!(!atom.atom().no_implicit(), "{symbol}");
    }
}

#[test]
fn mass_charge_hcount_and_optional_atom_fields_preserve_order_and_ranges() {
    let rich = atom_line(1.25, -2.5, 3.75, "C", 1, 3, 2, 2, 1, 4, 2, 3, 12, 1, 5);
    let record = query(&molblock("rich", &[rich], &[]));
    let atom = &record.query.atoms()[0];
    assert_eq!(atom.atom().formal_charge(), 1);
    assert_eq!(atom.atom().isotope(), Some(13));
    assert_eq!(atom.atom().mol_parity(), Some(2));
    assert_eq!(atom.atom().atom_map(), Some(12));
    assert_eq!(atom.atom().mol_inversion_flag(), Some(1));
    assert_eq!(atom.atom().prop("molStereoCare"), Some("1"));
    assert_eq!(atom.atom().prop("molTotValence"), Some("4"));
    assert_eq!(atom.atom().prop("molRxnRole"), Some("2"));
    assert_eq!(atom.atom().prop("molRxnComponent"), Some("3"));
    assert_eq!(atom.atom().prop("molRxnExactChange"), Some("5"));
    assert!(atom.atom().no_implicit());
    assert!(atom_query_contains(
        atom.predicate(),
        &AtomQueryPredicate::AtomicNumber(6)
    ));
    assert!(atom_query_contains(
        atom.predicate(),
        &AtomQueryPredicate::FormalCharge(1)
    ));
    assert!(atom_query_contains(
        atom.predicate(),
        &AtomQueryPredicate::ImplicitHydrogenCountLessEqual(1)
    ));
    assert_eq!(
        record.query.conformers_3d()[0].coordinates(),
        &[[1.25, -2.5, 3.75]]
    );

    for (symbol, mass, expected) in [("C", 1, 13), ("D", 0, 2), ("T", 0, 3)] {
        let atom = atom_line(0.0, 0.0, 0.0, symbol, mass, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let (topology, _) = concrete(&molblock("isotope", &[atom], &[]));
        assert_eq!(topology.atoms[0].isotope(), Some(expected));
    }
    let negative = atom_line(0.0, 0.0, 0.0, "H", -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    assert!(matches!(
        read_mol_block_detached(&molblock("negative", &[negative], &[])),
        Err(SdfReadError::Unsupported(message)) if message.contains("negative isotope")
    ));
}

#[test]
fn concrete_bond_types_ids_order_and_aromaticity_are_exact() {
    for (kind, expected) in [
        (0, BondOrder::Unspecified),
        (1, BondOrder::Single),
        (2, BondOrder::Double),
        (3, BondOrder::Triple),
        (4, BondOrder::Aromatic),
        (9, BondOrder::Dative),
    ] {
        let block = molblock(
            "bond",
            &[plain_atom("C"), plain_atom("N")],
            &[bond_line(2, 1, kind, 0, 0, 0)],
        );
        let (topology, _) = concrete(&block);
        let bond = &topology.bonds[0];
        assert_eq!(bond.id().index(), 0, "kind {kind}");
        assert_eq!((bond.begin().index(), bond.end().index()), (1, 0));
        assert_eq!(bond.order(), expected, "kind {kind}");
        assert_eq!(bond.is_aromatic(), kind == 4, "kind {kind}");
        assert_eq!(
            bond.prop("_MolFileBondType"),
            Some(kind.to_string().as_str())
        );
    }
}

#[test]
fn query_bond_types_and_topology_keep_complete_source_predicates() {
    for (kind, expected) in [
        (
            5,
            QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Double,
            ])),
        ),
        (
            6,
            QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Aromatic,
            ])),
        ),
        (
            7,
            QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Double,
                BondOrder::Aromatic,
            ])),
        ),
        (8, QueryNode::predicate(BondQueryPredicate::Any)),
        (42, QueryNode::predicate(BondQueryPredicate::Any)),
    ] {
        let record = query(&molblock(
            "query bond",
            &[plain_atom("C"), plain_atom("C")],
            &[bond_line(1, 2, kind, 0, 0, 0)],
        ));
        assert_eq!(record.query.bonds()[0].predicate(), &expected, "{kind}");
    }

    let ring_single = query(&molblock(
        "ring single",
        &[plain_atom("C"), plain_atom("C")],
        &[bond_line(1, 2, 1, 0, 1, 0)],
    ));
    assert_eq!(
        ring_single.query.bonds()[0].predicate(),
        &QueryNode::and(vec![
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
            QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
        ])
    );
    let nonring_any = query(&molblock(
        "nonring any",
        &[plain_atom("C"), plain_atom("C")],
        &[bond_line(1, 2, 8, 0, 2, 0)],
    ));
    assert_eq!(
        nonring_any.query.bonds()[0].predicate(),
        &QueryNode::predicate(BondQueryPredicate::IsInRing(false))
    );
}

#[test]
fn bond_stereo_reaction_optional_lexing_and_stereo_care_follow_source() {
    for (stereo, direction, bond_stereo) in [
        (1, BondDirection::BeginWedge, BondStereo::None),
        (3, BondDirection::EitherDouble, BondStereo::Any),
        (4, BondDirection::Unknown, BondStereo::None),
        (6, BondDirection::BeginDash, BondStereo::None),
    ] {
        let (topology, _) = concrete(&molblock(
            "stereo",
            &[plain_atom("C"), plain_atom("C")],
            &[bond_line(1, 2, 1, stereo, 0, 7)],
        ));
        assert_eq!(topology.bonds[0].direction(), direction, "{stereo}");
        assert_eq!(topology.bonds[0].stereo(), bond_stereo, "{stereo}");
        assert_eq!(topology.bonds[0].prop("molReactStatus"), Some("7"));
    }

    let lexical = "  1  2  1abc  0xyzbad".to_owned();
    let (topology, _) = concrete(&molblock(
        "optional lexical",
        &[plain_atom("C"), plain_atom("C")],
        &[lexical],
    ));
    assert_eq!(topology.bonds[0].direction(), BondDirection::None);
    assert_eq!(topology.bonds[0].prop("molReactStatus"), None);

    let care = atom_line(0.0, 0.0, 0.0, "C", 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0);
    let (topology, _) = concrete(&molblock(
        "care",
        &[care.clone(), care],
        &[bond_line(1, 2, 1, 0, 0, 0)],
    ));
    assert_eq!(topology.bonds[0].prop("molStereoCare"), Some("1"));
}

#[test]
fn bond_endpoint_topology_and_primary_field_errors_are_structured() {
    for (begin, end) in [(0, 1), (1, 0), (1, 3)] {
        let error = read_mol_block_detached(&molblock(
            "endpoint",
            &[plain_atom("C"), plain_atom("C")],
            &[bond_line(begin, end, 1, 0, 0, 0)],
        ))
        .unwrap_err();
        assert!(matches!(
            error,
            SdfReadError::Field {
                kind: "bond endpoint",
                ..
            }
        ));
    }
    let error = read_mol_block_detached(&molblock(
        "topology",
        &[plain_atom("C"), plain_atom("C")],
        &[bond_line(1, 2, 1, 0, 3, 0)],
    ))
    .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("Unrecognized bond topology specifier: 3")
    );

    for line in [" xx  2  1", "  1 xx  1", "  1  2 xx"] {
        let error = read_mol_block_detached(&molblock(
            "primary",
            &[plain_atom("C"), plain_atom("C")],
            &[line.to_owned()],
        ))
        .unwrap_err();
        assert!(error.to_string().contains("Cannot convert"), "{error}");
    }
}

#[test]
fn coordinates_preserve_row_order_dimension_and_crlf_fixed_width_input() {
    let atoms = [
        atom_line(1.0, 2.0, 0.0, "C", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        atom_line(-3.0, 4.0, 0.0, "O", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    ];
    let block = molblock("coords", &atoms, &[bond_line(1, 2, 1, 0, 0, 0)]);
    let (topology, coordinates) = concrete(&block.replace('\n', "\r\n"));
    assert_eq!(topology.atoms[0].element().symbol(), "C");
    assert_eq!(topology.atoms[1].element().symbol(), "O");
    assert_eq!(coordinates.conformers_2d.len(), 1);
    assert_eq!(
        coordinates.conformers_2d[0].coordinates(),
        &[[1.0, 2.0], [-3.0, 4.0]]
    );
    assert!(coordinates.conformers_3d.is_empty());
}
